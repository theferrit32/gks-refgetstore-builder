#!/usr/bin/env python
"""The source model: what a build is configured to read, and from where.

This is the declarative half of the builder -- the manifest schema
(``AssemblyConfig``, ``SeqsetConfig``, ``RecordExclusion``), the concrete file
list those entries resolve to (``ResolvedSource``), the provider checksum
manifests consulted while resolving, and the cache-path and digest helpers that
name and fingerprint the resulting files.

It deliberately knows nothing about ingestion. ``build_store`` (the engine) and
``build_lock`` (the record of what an engine run consumed and produced) both
depend on this module and not on each other's internals, which is what keeps
the module graph acyclic.
"""

from __future__ import annotations

import fnmatch
import hashlib
import logging
import re
import subprocess
import time
import tomllib
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator
from urllib.parse import urljoin, urlsplit

logger = logging.getLogger("sources")

NA_VALUES = {"na", "", "<NA>"}
SEQSET_FORMATS = ("fasta", "gbff", "lrg_zip")

# Suffixes appended to a source artifact to name its derived FASTA. Shared with
# build_lock.records_from_log, which recovers the source cache path by stripping
# one of these; a resolver that invents its own suffix breaks that mapping
# silently, so add it here rather than hardcoding it at either site. The
# functions that actually produce these names are registered in
# build_store.DERIVED_FASTA_RESOLVERS -- keep the two in step.
DERIVED_SUFFIXES = (".filtered.fa.gz", ".fasta")


@dataclass(frozen=True)
class RecordExclusion:
    """Declarative rule for dropping records from named file classes.

    Exclusion is by record *name*, never by inspecting sequence content. Ensembl
    releases 76-109 mark their N-padded alt/patch scaffolds with a ``CHR_``
    prefix, which is both sound and complete for that set, so a name rule is
    exact where an N-fraction heuristic would be a guess.
    """

    file_classes: tuple[str, ...]
    record_prefixes: tuple[str, ...]

    def applies_to(self, file_class: str | None) -> bool:
        return file_class is not None and file_class in self.file_classes

    def matches(self, record_name: str) -> bool:
        return record_name.startswith(self.record_prefixes)


@dataclass
class AssemblyConfig:
    namespace: str
    report_url: str
    fasta_url: str | None = None
    load_fasta: bool = True
    fasta_path: str | None = None
    checksum_manifest_url: str | None = None

    def __post_init__(self) -> None:
        if self.load_fasta and not self.fasta_url:
            raise ValueError(
                f"assembly {self.namespace!r}: fasta_url is required when "
                "load_fasta is true"
            )


@dataclass
class SeqsetConfig:
    """Flat FASTA (optionally sharded) where the accession is the header name.

    Used for NCBI RefSeq transcript/protein shards that have no sidecar
    assembly_report.txt. Each shard is ingested as its own collection; header
    names become ``namespace:<name>`` aliases on the resulting digests.
    """

    name: str
    namespace: str
    url_template: str | None = None
    shard_range: list[int] | None = None
    urls: list[str] | None = None
    url_pattern: str | None = None
    checksum_manifest_url: str | None = None
    checksum_manifest_urls: list[str] | None = None
    md5_manifest_urls: list[str] | None = None
    file_class: str | None = None
    file_classes: list[str] | None = None
    release: int | None = None
    rolling_namespace: str | None = None
    format: str = "fasta"
    exclude: dict | None = None
    resolved_sources: list["ResolvedSource"] | None = field(default=None, init=False)
    exclusion: "RecordExclusion | None" = field(default=None, init=False)

    def __post_init__(self) -> None:
        modes = sum(value is not None for value in
                    (self.url_template, self.urls, self.url_pattern))
        if modes != 1:
            raise ValueError(
                f"seqset {self.name!r}: set exactly one of url_template, urls, "
                "or url_pattern"
            )
        if self.shard_range is not None and self.url_template is None:
            raise ValueError(
                f"seqset {self.name!r}: shard_range is only valid with url_template"
            )
        if (self.url_pattern is None) != (self.checksum_manifest_url is None):
            raise ValueError(
                f"seqset {self.name!r}: url_pattern and checksum_manifest_url "
                "must be set together"
            )
        if self.url_pattern is not None:
            pattern = urlsplit(self.url_pattern)
            manifest = urlsplit(self.checksum_manifest_url or "")
            if pattern.scheme != "https" or manifest.scheme != "https":
                raise ValueError(f"seqset {self.name!r}: pattern sources require HTTPS")
            if not any(c in Path(pattern.path).name for c in "*?["):
                raise ValueError(f"seqset {self.name!r}: url_pattern basename needs a glob")
            if (pattern.scheme, pattern.netloc, str(Path(pattern.path).parent)) != (
                manifest.scheme, manifest.netloc, str(Path(manifest.path).parent)
            ):
                raise ValueError(
                    f"seqset {self.name!r}: pattern and checksum manifest must share a directory"
                )
        if self.checksum_manifest_urls is not None:
            if self.urls is None:
                raise ValueError(
                    f"seqset {self.name!r}: checksum_manifest_urls requires urls"
                )
            if len(self.checksum_manifest_urls) != len(self.urls):
                raise ValueError(
                    f"seqset {self.name!r}: checksum_manifest_urls must match urls"
                )
            if self.checksum_manifest_url is not None:
                raise ValueError(
                    f"seqset {self.name!r}: use only one checksum manifest mode"
                )
        if self.md5_manifest_urls is not None:
            if self.urls is None or len(self.md5_manifest_urls) != len(self.urls):
                raise ValueError(
                    f"seqset {self.name!r}: md5_manifest_urls must match urls"
                )
            if self.checksum_manifest_urls is not None:
                raise ValueError(
                    f"seqset {self.name!r}: use only one explicit checksum mode"
                )
        if self.file_class is not None and self.file_classes is not None:
            raise ValueError(
                f"seqset {self.name!r}: use only one file class mode"
            )
        source_count = len(self.urls) if self.urls is not None else 1
        if self.file_classes is not None and len(self.file_classes) != source_count:
            raise ValueError(
                f"seqset {self.name!r}: file_classes must match source URLs"
            )
        if (self.release is None) != (self.rolling_namespace is None):
            raise ValueError(
                f"seqset {self.name!r}: release and rolling_namespace must be set together"
            )
        if self.release is not None:
            expected = f"{self.rolling_namespace}-{self.release}"
            if self.namespace != expected:
                raise ValueError(
                    f"seqset {self.name!r}: release namespace must be {expected!r}"
                )
        if self.format not in SEQSET_FORMATS:
            raise ValueError(
                f"seqset {self.name!r}: format must be one of "
                f"{', '.join(repr(f) for f in SEQSET_FORMATS)}, "
                f"got {self.format!r}"
            )
        self.exclusion = self._parse_exclude()

    def _parse_exclude(self) -> "RecordExclusion | None":
        """Validate the ``exclude`` table and freeze it into a RecordExclusion."""
        if self.exclude is None:
            return None
        if not isinstance(self.exclude, dict):
            raise ValueError(f"seqset {self.name!r}: exclude must be a table")
        unknown = set(self.exclude) - {"file_classes", "record_prefixes"}
        if unknown:
            raise ValueError(
                f"seqset {self.name!r}: unknown exclude key(s) "
                f"{', '.join(sorted(repr(k) for k in unknown))}"
            )
        classes = tuple(self.exclude.get("file_classes") or ())
        prefixes = tuple(self.exclude.get("record_prefixes") or ())
        if not classes:
            raise ValueError(f"seqset {self.name!r}: exclude.file_classes is required")
        if not prefixes:
            raise ValueError(f"seqset {self.name!r}: exclude.record_prefixes is required")
        known = set(self.file_classes or ([self.file_class] if self.file_class else []))
        missing = [c for c in classes if c not in known]
        if missing:
            raise ValueError(
                f"seqset {self.name!r}: exclude.file_classes "
                f"{', '.join(repr(c) for c in missing)} not among this seqset's "
                f"file classes {sorted(known)}"
            )
        return RecordExclusion(file_classes=classes, record_prefixes=prefixes)

    def ingest_spec(self, file_class: str | None) -> dict | None:
        """Declared transformation applied to a source of ``file_class``.

        Recorded in the build lock so a config change that alters *what gets
        ingested* from *unchanged upstream bytes* is detectable. The lock's
        ``sha256`` pins the bytes read; this pins how they are interpreted.

        ``None`` means "ingested as published" -- no derived-format conversion
        and no record exclusion. Most sources are in that case, so returning
        ``None`` rather than an empty dict keeps the lock free of noise.

        The exclusion is included only when it applies to *this* file class:
        a seqset declares exclusions per class, and editing the rule for
        ``dna.toplevel`` must not mark its ``cdna`` sibling as changed.
        """
        spec: dict = {}
        if self.format != "fasta":
            spec["format"] = self.format
        if self.exclusion is not None and self.exclusion.applies_to(file_class):
            spec["exclude"] = {
                "record_prefixes": list(self.exclusion.record_prefixes)
            }
        return spec or None

    def file_class_for_index(self, index: int) -> str | None:
        """File class of the ``index``-th resolved source (0-based).

        ``iter_shard_urls`` drops per-index metadata, so callers that need the
        file class alongside the URL enumerate and ask here.
        """
        if self.file_classes is not None:
            return self.file_classes[index] if index < len(self.file_classes) else None
        return self.file_class

    def iter_shard_urls(self) -> Iterator[tuple[str, str]]:
        """Yield (shard_label, url) pairs. Shard label is "" for unsharded."""
        if self.resolved_sources is not None:
            for i, source in enumerate(self.resolved_sources, 1):
                yield str(i), source.url
            return
        if self.url_pattern is not None:
            raise RuntimeError(f"seqset {self.name!r}: pattern source has not been resolved")
        if self.urls is not None:
            for i, url in enumerate(self.urls, 1):
                yield str(i), url
            return
        if self.shard_range is None:
            yield "", self.url_template
            return
        lo, hi = self.shard_range
        for i in range(lo, hi + 1):
            yield str(i), self.url_template.replace("{shard}", str(i))


@dataclass(frozen=True)
class ResolvedSource:
    kind: str
    owner: str
    url: str
    upstream_md5: str | None = None
    provider_checksum: str | None = None
    provider_checksum_algorithm: str | None = None
    provider_checksum_blocks: int | None = None
    checksum_url: str | None = None
    file_class: str | None = None


_MD5_RECORD = re.compile(r"^([0-9A-Fa-f]{32})[ \t]+([^\s]+)$")


def parse_checksum_manifest(text: str) -> dict[str, str]:
    """Parse strict NCBI ``MD5 filename`` records keyed by safe basenames."""
    records: dict[str, str] = {}
    for line_number, raw in enumerate(text.splitlines(), 1):
        match = _MD5_RECORD.fullmatch(raw)
        if not match:
            raise ValueError(f"malformed checksum manifest line {line_number}: {raw!r}")
        md5, filename = match.groups()
        if (filename in (".", "..") or Path(filename).name != filename
                or "\\" in filename or any(char in filename for char in "*?[]")):
            raise ValueError(f"unsafe checksum manifest filename on line {line_number}: {filename!r}")
        if filename in records:
            raise ValueError(f"duplicate checksum manifest filename: {filename!r}")
        records[filename] = md5.lower()
    return records


def parse_ncbi_md5checksums(text: str) -> dict[str, str]:
    """Parse NCBI directory ``md5checksums.txt`` records by basename."""
    records: dict[str, str] = {}
    pattern = re.compile(r"^([0-9A-Fa-f]{32})[ \t]+[*]?(.+)$")
    for line_number, raw in enumerate(text.splitlines(), 1):
        match = pattern.fullmatch(raw.strip())
        if not match:
            raise ValueError(f"malformed NCBI checksum line {line_number}: {raw!r}")
        md5, raw_name = match.groups()
        name = raw_name[2:] if raw_name.startswith("./") else raw_name
        if not name or name in (".", "..") or "\\" in name:
            raise ValueError(f"unsafe NCBI checksum filename on line {line_number}")
        # Assembly manifests recurse into nested auxiliary directories, which
        # can reuse basenames. Configured FASTA/report inputs live at the
        # manifest directory root, so nested records are intentionally ignored.
        if "/" in name:
            continue
        basename = Path(name).name
        if basename in records and records[basename] != md5.lower():
            raise ValueError(f"conflicting NCBI checksum basename: {basename!r}")
        records[basename] = md5.lower()
    return records


_ENSEMBL_CHECKSUM_RECORD = re.compile(
    r"^([0-9]{1,5})[ \t]+([0-9]+)[ \t]+([^\s]+)$"
)


def parse_ensembl_checksum_manifest(text: str) -> dict[str, tuple[str, int]]:
    """Parse Ensembl's BSD-sum ``CHECKSUMS`` records."""
    records: dict[str, tuple[str, int]] = {}
    for line_number, raw in enumerate(text.splitlines(), 1):
        match = _ENSEMBL_CHECKSUM_RECORD.fullmatch(raw.strip())
        if not match:
            raise ValueError(
                f"malformed Ensembl checksum line {line_number}: {raw!r}"
            )
        checksum, blocks, filename = match.groups()
        if (filename in (".", "..") or Path(filename).name != filename
                or "\\" in filename or any(char in filename for char in "*?[]")):
            raise ValueError(
                f"unsafe Ensembl checksum filename on line {line_number}: {filename!r}"
            )
        if filename in records:
            raise ValueError(f"duplicate Ensembl checksum filename: {filename!r}")
        records[filename] = (checksum.zfill(5), int(blocks))
    return records


def bsd_sum_file(path: Path) -> tuple[str, int]:
    """Return Ensembl's 16-bit BSD checksum and 1024-byte block count."""
    try:
        result = subprocess.run(
            ["sum", str(path)], capture_output=True, text=True, check=True
        )
        checksum, blocks, *_ = result.stdout.split()
        return checksum.zfill(5), int(blocks)
    except (OSError, subprocess.SubprocessError, ValueError):
        # Portable fallback for minimal environments without POSIX ``sum``.
        pass
    checksum = 0
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            size += len(chunk)
            for value in chunk:
                checksum = ((checksum >> 1) | ((checksum & 1) << 15))
                checksum = (checksum + value) & 0xFFFF
    return f"{checksum:05d}", (size + 1023) // 1024


def natural_sort_key(value: str) -> tuple:
    return tuple((1, int(part)) if part.isdigit() else (0, part.casefold())
                 for part in re.split(r"(\d+)", value))


def _fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "gks-refgetstore-builder"})
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt == 5:
                raise
            delay = attempt * 2
            logger.warning(
                "checksum fetch failed (%s); retrying %s in %ds",
                exc, url, delay,
            )
            time.sleep(delay)
    raise AssertionError("unreachable")


def resolve_sources(
    assemblies: list[AssemblyConfig],
    seqsets: list[SeqsetConfig],
    fetch_text: Callable[[str], str] | None = None,
) -> list[ResolvedSource]:
    """Resolve every config entry once, sharing fetched manifests by URL."""
    manifest_cache: dict[str, dict[str, str]] = {}
    ncbi_manifest_cache: dict[str, dict[str, str]] = {}
    ensembl_manifest_cache: dict[str, dict[str, tuple[str, int]]] = {}
    if fetch_text is None:
        manifest_urls = {
            url
            for seqset in seqsets
            for url in (
                ([seqset.checksum_manifest_url]
                 if seqset.checksum_manifest_url else [])
                + (seqset.checksum_manifest_urls or [])
                + (seqset.md5_manifest_urls or [])
            )
        }
        manifest_urls.update(
            assembly.checksum_manifest_url
            for assembly in assemblies if assembly.checksum_manifest_url
        )
        with ThreadPoolExecutor(max_workers=6) as executor:
            fetched = dict(zip(
                sorted(manifest_urls), executor.map(_fetch_text, sorted(manifest_urls))
            ))
        fetch_text = fetched.__getitem__
    resolved: list[ResolvedSource] = []
    for seqset in seqsets:
        seq_sources: list[ResolvedSource] = []
        if seqset.url_pattern is not None:
            manifest_url = seqset.checksum_manifest_url
            assert manifest_url is not None
            if manifest_url not in manifest_cache:
                manifest_cache[manifest_url] = parse_checksum_manifest(fetch_text(manifest_url))
            basename_pattern = Path(urlsplit(seqset.url_pattern).path).name
            matches = [(name, md5) for name, md5 in manifest_cache[manifest_url].items()
                       if fnmatch.fnmatchcase(name, basename_pattern)]
            if not matches:
                raise ValueError(f"seqset {seqset.name!r}: url_pattern matched no manifest files")
            base = seqset.url_pattern.rsplit("/", 1)[0] + "/"
            for name, md5 in sorted(matches, key=lambda item: natural_sort_key(item[0])):
                seq_sources.append(ResolvedSource(
                    "seqset", seqset.name, urljoin(base, name), md5,
                    file_class=seqset.file_class,
                ))
        else:
            urls = list(seqset.iter_shard_urls())
            for index, (_, url) in enumerate(urls):
                md5_url = (
                    seqset.md5_manifest_urls[index]
                    if seqset.md5_manifest_urls is not None else None
                )
                checksum_url = (
                    seqset.checksum_manifest_urls[index]
                    if seqset.checksum_manifest_urls is not None else None
                )
                checksum = None
                blocks = None
                upstream_md5 = None
                if md5_url is not None:
                    if md5_url not in ncbi_manifest_cache:
                        ncbi_manifest_cache[md5_url] = parse_ncbi_md5checksums(
                            fetch_text(md5_url)
                        )
                    basename = Path(urlsplit(url).path).name
                    try:
                        upstream_md5 = ncbi_manifest_cache[md5_url][basename]
                    except KeyError as exc:
                        raise ValueError(
                            f"seqset {seqset.name!r}: {basename!r} absent from "
                            f"{md5_url}"
                        ) from exc
                if checksum_url is not None:
                    if checksum_url not in ensembl_manifest_cache:
                        ensembl_manifest_cache[checksum_url] = (
                            parse_ensembl_checksum_manifest(fetch_text(checksum_url))
                        )
                    basename = Path(urlsplit(url).path).name
                    try:
                        checksum, blocks = ensembl_manifest_cache[checksum_url][basename]
                    except KeyError as exc:
                        raise ValueError(
                            f"seqset {seqset.name!r}: {basename!r} absent from "
                            f"{checksum_url}"
                        ) from exc
                file_class = (
                    seqset.file_classes[index]
                    if seqset.file_classes is not None else seqset.file_class
                )
                seq_sources.append(ResolvedSource(
                    "seqset", seqset.name, url,
                    upstream_md5=upstream_md5,
                    provider_checksum=checksum,
                    provider_checksum_algorithm="bsd-sum" if checksum else None,
                    provider_checksum_blocks=blocks,
                    checksum_url=checksum_url or md5_url,
                    file_class=file_class,
                ))
        seqset.resolved_sources = seq_sources
        resolved.extend(seq_sources)
    for assembly in assemblies:
        checksums: dict[str, str] = {}
        if assembly.checksum_manifest_url is not None:
            checksum_url = assembly.checksum_manifest_url
            if checksum_url not in ncbi_manifest_cache:
                ncbi_manifest_cache[checksum_url] = parse_ncbi_md5checksums(
                    fetch_text(checksum_url)
                )
            checksums = ncbi_manifest_cache[checksum_url]
        if assembly.load_fasta and not assembly.fasta_path:
            assert assembly.fasta_url is not None
            basename = Path(urlsplit(assembly.fasta_url).path).name
            if assembly.checksum_manifest_url is not None and basename not in checksums:
                raise ValueError(
                    f"assembly {assembly.namespace!r}: {basename!r} absent from "
                    f"{assembly.checksum_manifest_url}"
                )
            resolved.append(ResolvedSource(
                "assembly_fasta", assembly.namespace, assembly.fasta_url,
                checksums.get(basename), checksum_url=assembly.checksum_manifest_url,
                file_class="genomic",
            ))
        report_basename = Path(urlsplit(assembly.report_url).path).name
        report_md5 = checksums.get(report_basename)
        resolved.append(ResolvedSource(
            "assembly_report", assembly.namespace, assembly.report_url,
            report_md5,
            checksum_url=(assembly.checksum_manifest_url if report_md5 else None),
            file_class="assembly_report",
        ))
    urls = [source.url for source in resolved]
    if len(urls) != len(set(urls)):
        duplicate = next(url for url in urls if urls.count(url) > 1)
        raise ValueError(f"duplicate resolved source URL: {duplicate}")
    return resolved


# ``apply_locked_sources`` lives in build_lock: it reads the lock through that
# module's accessors, and keeping it there leaves exactly one place that knows
# which lock schema carries concrete sources.


def load_config(
    path: Path,
) -> tuple[list[AssemblyConfig], list[SeqsetConfig]]:
    with path.open("rb") as fh:
        data = tomllib.load(fh)
    assemblies = [AssemblyConfig(**entry) for entry in data.get("assembly", [])]
    seqsets = [SeqsetConfig(**entry) for entry in data.get("seqset", [])]
    names = [entry.name for entry in seqsets]
    if len(names) != len(set(names)):
        duplicate = next(name for name in names if names.count(name) > 1)
        raise ValueError(f"duplicate seqset name: {duplicate!r}")
    release_groups: dict[tuple[str, int], list[SeqsetConfig]] = {}
    for entry in seqsets:
        if entry.release is not None and entry.rolling_namespace is not None:
            release_groups.setdefault(
                (entry.rolling_namespace, entry.release), []
            ).append(entry)
    expected_ensembl_classes = {"dna.toplevel", "cdna", "ncrna", "pep"}
    for (provider, release), entries in release_groups.items():
        classes = [
            file_class
            for entry in entries
                for file_class in (
                    entry.file_classes
                    or ([entry.file_class] if entry.file_class is not None else [])
                )
        ]
        if len(classes) != len(set(classes)):
            raise ValueError(
                f"{provider} release {release}: duplicate file class"
            )
        if provider == "ensembl" and set(classes) != expected_ensembl_classes:
            raise ValueError(
                f"ensembl release {release}: expected file classes "
                f"{sorted(expected_ensembl_classes)}, got {sorted(classes)}"
            )
    return assemblies, seqsets


def mirror_cache_path(cache_root: Path, url: str) -> Path:
    """Cache path that mirrors the source host + full path under ``cache_root``.

    ``https://ftp.ncbi.nlm.nih.gov/genomes/all/.../X_rna.fna.gz`` caches to
    ``<cache_root>/ftp.ncbi.nlm.nih.gov/genomes/all/.../X_rna.fna.gz``.

    Mirroring host + full path is collision-proof by construction: different
    annotation-release snapshots and per-patch assembly dirs routinely ship
    FASTAs with identical basenames (e.g. ``GCF_000001405.39_GRCh38.p13_rna.fna.gz``
    recurs across every dated AR109 release) but live at distinct paths, so
    caching by basename alone would silently reuse the wrong file. The layout is
    also human-navigable and lets the fetcher (``fetch_sources.py``) and the
    builder share one cache.
    """
    parts = urlsplit(url)
    return cache_root / parts.netloc / parts.path.lstrip("/")


def md5_file(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 - required provider checksum algorithm
    with path.open("rb") as stream:
        while chunk := stream.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    """SHA-256 of a file's bytes -- the lock's identity anchor for a source.

    Lives beside the source model rather than in ``build_lock`` because both
    sides need it: the lock records it, and ``build_store.ensure_download``
    compares cached bytes against what a lock pinned before deciding whether to
    re-fetch. ``build_lock`` re-exports it so lock callers need only that module.
    """
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def provider_checksum_matches(path: Path, source: ResolvedSource) -> bool:
    """Validate any provider checksum attached to a resolved source."""
    if source.provider_checksum_algorithm == "md5":
        return md5_file(path) == source.provider_checksum
    if source.upstream_md5 is not None:
        return md5_file(path) == source.upstream_md5
    if source.provider_checksum_algorithm == "bsd-sum":
        return bsd_sum_file(path) == (
            source.provider_checksum, source.provider_checksum_blocks
        )
    return source.provider_checksum is None
