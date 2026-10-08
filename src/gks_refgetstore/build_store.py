#!/usr/bin/env python
"""Build a custom gtars RefgetStore with vrs-python alias parity.

Ingests NCBI genomic FASTA files into a local RefgetStore and augments it with
per-sequence and per-collection aliases parsed from the matching NCBI assembly
reports. See README.md for details.
"""

from __future__ import annotations

import gzip
import os
import shutil
import zipfile
import logging
import re
import tempfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator

from gtars.refget import RefgetStore

from . import build_lock, fs_checks
from .sources import (NA_VALUES, AssemblyConfig, RecordExclusion,
                      ResolvedSource, SeqsetConfig, load_config, md5_file,
                      mirror_cache_path, natural_sort_key,
                      provider_checksum_matches, resolve_sources, sha256_file)

logger = logging.getLogger("build_store")

SQ_PREFIX = "SQ."

# Files imported concurrently per batched ingest call. Capped well below the
# core count because peak RSS grows with concurrency on the largest inputs
# (~0.45 GiB per additional in-flight Ensembl ``dna.toplevel``), and the
# throughput gain flattens out long before the memory cost does.
INGEST_JOBS_DEFAULT = min(8, os.cpu_count() or 4)

# Compression level for derived filtered FASTAs. Level 1 costs ~28s per Ensembl
# dna.toplevel and yields 3.15 GB -> 1.02 GB; level 6 costs ~375s to reach
# 0.88 GB, and the gzip module's own default (9) is slower still. These files are
# regenerable cache artifacts, so the fast setting is the right trade.
FILTERED_FASTA_COMPRESSLEVEL = 1


@dataclass
class AssemblyReport:
    """Parsed NCBI assembly_report.txt."""

    refseq_assembly_accession: str | None
    genbank_assembly_accession: str | None
    rows: list[dict[str, str]] = field(default_factory=list)


@dataclass
class AssemblyStats:
    namespace: str
    collection_digest: str | None = None
    rows_total: int = 0
    rows_resolved: int = 0
    rows_skipped_non_refseq: int = 0
    sequence_aliases_added: int = 0
    sequence_aliases_skipped: int = 0
    collection_aliases_added: int = 0
    warnings: int = 0


@dataclass
class SeqsetStats:
    name: str
    namespace: str
    shards_processed: int = 0
    sequences_ingested: int = 0
    sequence_aliases_added: int = 0
    sequence_aliases_skipped: int = 0
    records_excluded: int = 0
    warnings: int = 0


@dataclass
class DeferredAssemblyReport:
    """Assembly state retained until all configured sequences are ingested."""

    entry: AssemblyConfig
    report_path: Path
    stats: AssemblyStats
    collection_digest: str | None = None
    report: AssemblyReport | None = None


def free_gib(path: Path) -> float:
    """Free space in GiB on the filesystem holding ``path``."""
    usage = shutil.disk_usage(path)
    return usage.free / 2**30


# Rough store-size multiplier over the compressed bytes of the ingested
# sources. Sequences are 2-bit encoded and deduplicated across collections,
# but the per-sequence and per-collection indexes add overhead, and the whole
# thing lands in the same order of magnitude as the gzipped inputs. This is a
# ballpark for a preflight warning, not an accounting.
STORE_SIZE_FACTOR = 1.2


def check_free_space(cache_dir: Path, store_dir: Path,
                     sources: list[ResolvedSource], min_free_gb: float) -> None:
    """Warn about the projected store size and stop if the disk is already tight.

    Checked once up front rather than monitored during the build: a full build
    runs for hours, and finding out at the end that the volume filled is far
    worse than being told at the start that it will.
    """
    cached = sum(
        path.stat().st_size
        for path in (mirror_cache_path(cache_dir, s.url) for s in sources)
        if path.exists()
    ) / 2**30
    projected = cached * STORE_SIZE_FACTOR
    free = free_gib(store_dir if store_dir.exists() else store_dir.parent)
    logger.info(
        "disk preflight: %.1f GiB of sources -> ~%.0f GiB store (rough); "
        "%.1f GiB free", cached, projected, free,
    )
    if free < min_free_gb:
        raise SystemExit(
            f"only {free:.1f} GiB free, below --min-free-gb {min_free_gb}; "
            "free space or lower the threshold"
        )
    if free < projected:
        logger.warning(
            "projected store (~%.0f GiB) exceeds free space (%.1f GiB); the "
            "build may run out part-way", projected, free,
        )


def gbff_to_fasta(gbff_path: Path, out_path: Path) -> int:
    """Convert a (gzipped) GenBank flat file to FASTA, keyed by accession.version.

    Streams the GBFF: each record's ``VERSION`` line supplies the FASTA header
    (the accession.version, which becomes the ``namespace:<accn.ver>`` alias),
    and the ``ORIGIN``..``//`` block supplies the sequence. Only the accession
    and sequence are needed for refget ingestion, so all other GenBank fields
    are skipped. Returns the number of records written.

    Used for the curated-history source
    (``…_knownrefseq_rna.gbff.gz``), which is distributed only as GBFF and holds
    replaced/suppressed ``NM_``/``NR_`` versions available in no bulk FASTA.
    """
    opener = gzip.open if gbff_path.suffix == ".gz" else open
    records = 0
    version: str | None = None
    seq_parts: list[str] = []
    in_origin = False
    tmp = out_path.with_suffix(out_path.suffix + ".part")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with opener(gbff_path, "rt") as fh, tmp.open("w") as out:  # type: ignore[operator]
        for line in fh:
            if in_origin:
                if line.startswith("//"):
                    if version and seq_parts:
                        out.write(f">{version}\n{''.join(seq_parts)}\n")
                        records += 1
                    version, seq_parts, in_origin = None, [], False
                else:
                    # ORIGIN sequence lines: leading base-count, then 10-char
                    # blocks of letters separated by spaces.
                    seq_parts.append("".join(c for c in line if c.isalpha()))
            elif line.startswith("VERSION"):
                version = line.split()[1] if len(line.split()) > 1 else None
            elif line.startswith("ORIGIN"):
                in_origin = True
    tmp.replace(out_path)
    return records


def resolve_gbff_fasta(gbff_path: Path, force: bool) -> Path:
    """Return a cached FASTA rendering of ``gbff_path``, converting if needed."""
    out_path = gbff_path.parent / (gbff_path.name + ".fasta")
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        logger.info("using cached gbff->fasta %s", out_path)
        return out_path
    logger.info("converting gbff -> fasta %s", gbff_path)
    n = gbff_to_fasta(gbff_path, out_path)
    logger.info("  wrote %d records to %s", n, out_path)
    return out_path


LRG_MEMBER_RE = re.compile(r"^(?:.*/)?(LRG_\d+)\.fasta$")
_LRG_GENOMIC_HEADER_RE = re.compile(rb"^>LRG_\d+g\b")


def lrg_zip_to_fasta(zip_path: Path, out_path: Path) -> tuple[int, int]:
    """Concatenate every ``LRG_N.fasta`` member of an EBI LRG bundle into one FASTA.

    EBI publishes the public LRG set both as 1325 standalone per-locus FASTAs and
    as one aggregate zip of the same bytes. Ingesting the aggregate pins a single
    immutable artifact instead of 1325 URLs, so this renders it back into the flat
    FASTA gtars expects.

    Members are emitted byte-faithfully — headers are *not* rewritten, so the
    resulting collection digest is a pure function of the upstream bytes plus the
    documented ordering. Members are ordered by ``natural_sort_key`` on the locus
    id so the digest does not depend on the archive's physical member order, and a
    newline is inserted where a member does not end in one so records cannot fuse
    across a member boundary.

    Returns ``(members, records)``.
    """
    tmp = out_path.with_suffix(out_path.suffix + ".part")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    members: list[tuple[str, zipfile.ZipInfo]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(zip_path) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            match = LRG_MEMBER_RE.match(info.filename)
            if match is None:
                # A silent skip would let an upstream repack change the ingested
                # contents without any signal.
                raise ValueError(
                    f"{zip_path.name}: unexpected member {info.filename!r}"
                )
            locus = match.group(1)
            if locus in seen:
                raise ValueError(f"{zip_path.name}: duplicate member {locus}.fasta")
            seen.add(locus)
            members.append((locus, info))
        if not members:
            raise ValueError(f"{zip_path.name}: no LRG_N.fasta members")
        members.sort(key=lambda item: natural_sort_key(item[0]))

        records = 0
        with tmp.open("wb") as out:
            for locus, info in members:
                member_records = 0
                first_chunk = True
                last_byte = b""
                with archive.open(info) as handle:
                    while chunk := handle.read(1 << 20):
                        if first_chunk:
                            if not _LRG_GENOMIC_HEADER_RE.match(chunk):
                                raise ValueError(
                                    f"{zip_path.name}: {locus}.fasta does not start "
                                    "with a genomic LRG header"
                                )
                            first_chunk = False
                        if last_byte == b"\n" and chunk.startswith(b">"):
                            member_records += 1  # header split across the boundary
                        member_records += chunk.count(b"\n>")
                        out.write(chunk)
                        last_byte = chunk[-1:]
                if first_chunk:
                    raise ValueError(f"{zip_path.name}: {locus}.fasta is empty")
                member_records += 1  # the leading header, which has no preceding \n
                if member_records < 2:
                    raise ValueError(
                        f"{zip_path.name}: {locus}.fasta has {member_records} record(s), "
                        "expected a genomic record plus at least one transcript"
                    )
                if last_byte != b"\n":
                    out.write(b"\n")
                records += member_records
    tmp.replace(out_path)
    return len(members), records


def resolve_lrg_fasta(zip_path: Path, force: bool) -> Path:
    """Return a cached FASTA rendering of an LRG bundle, unpacking if needed."""
    out_path = zip_path.parent / (zip_path.name + ".fasta")
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        logger.info("using cached lrg->fasta %s", out_path)
        return out_path
    logger.info("converting lrg zip -> fasta %s", zip_path)
    members, records = lrg_zip_to_fasta(zip_path, out_path)
    logger.info(
        "  wrote %d records from %d members to %s", records, members, out_path
    )
    return out_path


# Per-format conversion to a FASTA the store can ingest. The suffix each
# resolver appends must be registered in sources.DERIVED_SUFFIXES: that is how
# build_lock.records_from_log recovers the originating source from a derived
# path, and an unregistered suffix breaks the mapping silently.
DERIVED_FASTA_RESOLVERS: dict[str, Callable[[Path, bool], Path]] = {
    "gbff": resolve_gbff_fasta,
    "lrg_zip": resolve_lrg_fasta,
}


def filtered_fasta_path(src: Path) -> Path:
    """Cache path of ``src``'s filtered rendering, a sibling of the artifact."""
    return src.parent / (src.name + ".filtered.fa.gz")


def excluded_log_path(src: Path) -> Path:
    """Cache path of the audit log listing records dropped from ``src``."""
    return src.parent / (src.name + ".excluded.tsv")


def filter_fasta_records(
    src: Path, out_path: Path, exclusion: RecordExclusion, audit_log: Path | None = None
) -> tuple[int, int]:
    """Copy ``src`` to ``out_path``, dropping records whose name matches a prefix.

    Streams gzip->gzip, copying kept records' lines verbatim so their digests are
    unaffected. Returns ``(kept, dropped)``.

    Every dropped record's name and line count is written to ``audit_log`` so an
    exclusion is inspectable after the fact rather than inferred from a shrunken
    store. Raises if the rule matches nothing: a declared exclusion that silently
    no-ops would quietly retain everything it was meant to remove.
    """
    opener = gzip.open if src.suffix == ".gz" else open
    tmp = out_path.with_suffix(out_path.suffix + ".part")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    kept = dropped = 0
    dropped_rows: list[tuple[str, int]] = []
    emit = True
    name = ""
    lines = 0
    with opener(src, "rb") as fh, gzip.open(  # type: ignore[operator]
        tmp, "wb", compresslevel=FILTERED_FASTA_COMPRESSLEVEL
    ) as out:
        for line in fh:
            if line.startswith(b">"):
                if not emit and name:
                    dropped_rows.append((name, lines))
                name = line[1:].split()[0].decode("utf-8", "replace") if len(line) > 1 else ""
                emit = not exclusion.matches(name)
                kept += emit
                dropped += not emit
                lines = 0
            lines += 1
            if emit:
                out.write(line)
    if not emit and name:
        dropped_rows.append((name, lines))
    if not dropped:
        raise ValueError(
            f"exclusion matched no records in {src.name} "
            f"(prefixes {', '.join(repr(p) for p in exclusion.record_prefixes)})"
        )
    if audit_log is not None:
        audit_tmp = audit_log.with_suffix(audit_log.suffix + ".part")
        with audit_tmp.open("w", encoding="utf-8") as fh:
            fh.write("record_name\tlines\n")
            for record_name, n in dropped_rows:
                fh.write(f"{record_name}\t{n}\n")
        audit_tmp.replace(audit_log)
    tmp.replace(out_path)
    return kept, dropped


def resolve_filtered_fasta(src: Path, exclusion: RecordExclusion, force: bool) -> Path:
    """Return a cached filtered rendering of ``src``, producing it if needed.

    Same cache contract as ``resolve_gbff_fasta``/``resolve_lrg_fasta``: a sibling
    of the downloaded artifact, reused when present and non-empty unless ``force``.
    """
    out_path = filtered_fasta_path(src)
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        logger.info("using cached filtered fasta %s", out_path)
        return out_path
    logger.info("filtering %s", src)
    kept, dropped = filter_fasta_records(
        src, out_path, exclusion, audit_log=excluded_log_path(src)
    )
    logger.info("  kept %d record(s), dropped %d -> %s", kept, dropped, out_path)
    return out_path


def _excluded_count(src: Path) -> int:
    """Rows in ``src``'s exclusion audit log, or 0 if it is absent."""
    audit = excluded_log_path(src)
    if not audit.exists():
        return 0
    with audit.open("r", encoding="utf-8") as fh:
        return max(0, sum(1 for _ in fh) - 1)  # minus the header


def prepare_filtered_sources(
    seqsets: list[SeqsetConfig], download_dir: Path, jobs: int, force: bool = False
) -> dict[str, Path]:
    """Materialize every declared filtered FASTA up front, concurrently.

    Downloads are already a preflight (``ensure_download`` runs over all sources
    before any seqset is processed); filtration is the same shape. Doing it here
    rather than inline lets the whole set run in parallel: zlib releases the GIL,
    and decompression dominates at roughly four minutes per Ensembl release, so
    threading turns hours into minutes.

    Returns a ``{filtered input path: filtered path}`` map. ``process_seqset``
    resolves the same paths independently, so this is a warm-up, not a handoff.

    The filter input is the **derived** FASTA, not the downloaded artifact, for
    formats that have one. ``process_seqset`` converts before it filters, and a
    seqset combining ``format="lrg_zip"`` with ``exclude`` would otherwise hand
    ``filter_fasta_records`` a ZIP here: it would parse no records, drop
    nothing, and raise "exclusion matched no records" -- blaming the rule for a
    path bug. Production never combines the two today, which is exactly why it
    would have gone unnoticed.
    """
    work: list[tuple[Path, RecordExclusion]] = []
    for entry in seqsets:
        if entry.exclusion is None:
            continue
        resolver = DERIVED_FASTA_RESOLVERS.get(entry.format)
        for index, (_label, url) in enumerate(entry.iter_shard_urls()):
            if not entry.exclusion.applies_to(entry.file_class_for_index(index)):
                continue
            target = mirror_cache_path(download_dir, url)
            if not (target.exists() and target.stat().st_size > 0):
                continue
            work.append((
                resolver(target, force) if resolver is not None else target,
                entry.exclusion,
            ))
    if not work:
        return {}
    logger.info("filtering %d source file(s) with %d job(s)", len(work), jobs)
    out: dict[str, Path] = {}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = {
            pool.submit(resolve_filtered_fasta, src, exclusion, force): src
            for src, exclusion in work
        }
        for future in as_completed(futures):
            src = futures[future]
            out[str(src)] = future.result()
    return out


def lrg_alias_names(name: str) -> tuple[str, ...]:
    """LRG_1g is EBI's record id; HGVS/ClinVar spell the genomic record LRG_1."""
    match = re.fullmatch(r"(LRG_\d+)g", name)
    return (name, match.group(1)) if match else (name,)


SEQSET_ALIAS_EXPANDERS: dict[str, Callable[[str], tuple[str, ...]]] = {
    "lrg_zip": lrg_alias_names,
}


def ensure_download(url: str, target: Path, force: bool,
                    expected_md5: str | None = None,
                    source: ResolvedSource | None = None,
                    locked_sha256: str | None = None) -> Path:
    """Make ``target`` hold ``url``'s bytes, preferring what a lock pinned.

    ``locked_sha256`` is the lock's SHA-256 for this file, when a lock covers
    it. A provider-checksum mismatch on an otherwise valid cached file means the
    provider republished in place, and for a lock-based tool the lock wins: the
    cached bytes are the only surviving copy of what the lock describes, so they
    are kept and the drift is reported. Re-downloading would silently destroy
    the baseline the lock exists to preserve -- upstream has already changed the
    MD5 of all 30 remaining ``mRNA_Prot`` shards, so this is live, not
    hypothetical. ``--force-download`` remains the override.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        if source is not None:
            checksum_ok = provider_checksum_matches(target, source)
        else:
            checksum_ok = expected_md5 is None or md5_file(target) == expected_md5
        if checksum_ok:
            logger.info("using cached %s", target)
            return target
        if locked_sha256 is not None and sha256_file(target) == locked_sha256:
            logger.warning(
                "provider republished %s in place; keeping the lock-pinned "
                "cached bytes (sha256=%s…). Re-lock with --force-lock to accept "
                "the new upstream, or --force-download to overwrite.",
                target, locked_sha256[:12],
            )
            return target
        logger.warning("cached provider checksum mismatch; re-downloading %s", target)
    logger.info("downloading %s -> %s", url, target)
    tmp = target.with_suffix(target.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)  # noqa: S310
    if expected_md5 is not None and md5_file(tmp) != expected_md5:
        tmp.unlink(missing_ok=True)
        raise ValueError(f"NCBI MD5 mismatch for {url}")
    if source is not None and not provider_checksum_matches(tmp, source):
        tmp.unlink(missing_ok=True)
        raise ValueError(f"provider checksum mismatch for {url}")
    tmp.replace(target)
    return target


def require_prepared_source(path: Path) -> Path:
    """Return a nonempty source prepared by preflight; never access the network."""
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"prepared source is missing or empty: {path}")
    return path


def parse_assembly_report(path: Path) -> AssemblyReport:
    """Parse an NCBI assembly_report.txt into a structured object.

    Header lines begin with ``#`` and include ``GenBank assembly accession`` /
    ``RefSeq assembly accession`` entries. The final ``#``-prefixed line is the
    column header for the data rows (a tab-separated table).
    """
    refseq_accn: str | None = None
    genbank_accn: str | None = None
    columns: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("#"):
                stripped = line.lstrip("#").strip()
                if stripped.startswith("GenBank assembly accession:"):
                    genbank_accn = stripped.split(":", 1)[1].strip() or None
                elif stripped.startswith("RefSeq assembly accession:"):
                    refseq_accn = stripped.split(":", 1)[1].strip() or None
                elif "\t" in stripped and "Sequence-Name" in stripped:
                    columns = stripped.split("\t")
                continue
            if columns is None:
                continue
            parts = line.split("\t")
            if len(parts) != len(columns):
                continue
            rows.append(dict(zip(columns, parts)))
    return AssemblyReport(
        refseq_assembly_accession=refseq_accn,
        genbank_assembly_accession=genbank_accn,
        rows=rows,
    )


def resolve_present_alias(store, namespace: str, alias: str, *, collection=False):
    """Metadata for an alias the caller already saw listed in ``namespace``.

    gtars returns ``None`` rather than raising when an alias points at a digest
    the store no longer holds. That happens whenever a collection is removed
    without reconciling the namespaces referencing it -- a partially applied
    ``sync``, for instance -- and dereferencing the ``None`` raises an
    ``AttributeError`` naming neither the alias nor the cause.
    """
    getter = (store.get_collection_metadata_by_alias if collection
              else store.get_sequence_metadata_by_alias)
    metadata = getter(namespace, alias)
    if metadata is None:
        raise ValueError(
            f"dangling alias {namespace}:{alias} is listed in the namespace but "
            "resolves to no stored entry. The namespace still references a "
            "removed collection; reconcile it before ingesting into it"
        )
    return metadata


def build_name_to_digest_map(
    store: RefgetStore, collection_digest: str
) -> dict[str, str]:
    store.load_collection(collection_digest)
    level2 = store.get_collection_level2(collection_digest)
    names = level2["names"]
    sequences = level2["sequences"]
    out: dict[str, str] = {}
    for name, seq_id in zip(names, sequences):
        digest = seq_id[len(SQ_PREFIX) :] if seq_id.startswith(SQ_PREFIX) else seq_id
        out[name] = digest
    return out


def build_indexed_name_to_digest_map(
    store: RefgetStore,
    required_names: set[str],
    pending_refseq: ReportAliasAccumulator | None = None,
) -> dict[str, str]:
    """Resolve requested RefSeq accessions through pending and stored indexes."""
    started = time.monotonic()
    out: dict[str, str] = {}
    for name in required_names:
        pending = (
            pending_refseq.pending_sequence_digest("refseq", name)
            if pending_refseq is not None else None
        )
        if pending is not None:
            out[name] = pending
            continue
        try:
            metadata = store.get_sequence_metadata_by_alias("refseq", name)
        except KeyError:
            metadata = None
        if metadata is not None:
            out[name] = metadata.sha512t24u
    logger.info(
        "report-only indexed lookup: resolved %d/%d RefSeq accessions in %.2fs",
        len(out),
        len(required_names),
        time.monotonic() - started,
    )
    return out


class ReportAliasAccumulator:
    """First-seen-wins batches for aliases derived from assembly reports."""

    def __init__(self, store: RefgetStore) -> None:
        self.store = store
        self.sequence: dict[str, dict[str, tuple[str, AssemblyStats]]] = {}
        self.collection: dict[str, dict[str, tuple[str, AssemblyStats]]] = {}
        self._existing_sequence: dict[str, set[str]] = {}
        self._existing_collection: dict[str, set[str]] = {}

    def _sequence_aliases(self, namespace: str) -> set[str]:
        if namespace not in self._existing_sequence:
            self._existing_sequence[namespace] = set(
                self.store.list_sequence_aliases(namespace) or []
            )
        return self._existing_sequence[namespace]

    def _collection_aliases(self, namespace: str) -> set[str]:
        if namespace not in self._existing_collection:
            self._existing_collection[namespace] = set(
                self.store.list_collection_aliases(namespace) or []
            )
        return self._existing_collection[namespace]

    def add_sequence(
        self, namespace: str, alias: str, digest: str, stats: AssemblyStats
    ) -> bool:
        pending = self.sequence.setdefault(namespace, {})
        if alias in pending:
            if pending[alias][0] != digest:
                raise ValueError(
                    f"immutable alias collision {namespace}:{alias}: "
                    f"{pending[alias][0]} != {digest}"
                )
            stats.sequence_aliases_skipped += 1
            return False
        if alias in self._sequence_aliases(namespace):
            metadata = resolve_present_alias(self.store, namespace, alias)
            if metadata.sha512t24u != digest:
                raise ValueError(
                    f"immutable alias collision {namespace}:{alias}: "
                    f"{metadata.sha512t24u} != {digest}"
                )
            stats.sequence_aliases_skipped += 1
            return False
        pending[alias] = (digest, stats)
        return True

    def add_collection(
        self, namespace: str, alias: str, digest: str, stats: AssemblyStats
    ) -> bool:
        pending = self.collection.setdefault(namespace, {})
        if alias in pending:
            if pending[alias][0] != digest:
                raise ValueError(
                    f"immutable collection alias collision {namespace}:{alias}"
                )
            return False
        if alias in self._collection_aliases(namespace):
            metadata = resolve_present_alias(
                self.store, namespace, alias, collection=True
            )
            if metadata.digest != digest:
                raise ValueError(
                    f"immutable collection alias collision {namespace}:{alias}"
                )
            return False
        pending[alias] = (digest, stats)
        return True

    def pending_sequence_digest(self, namespace: str, alias: str) -> str | None:
        entry = self.sequence.get(namespace, {}).get(alias)
        return entry[0] if entry else None

    @staticmethod
    def _write_tsv(entries: dict[str, tuple[str, AssemblyStats]]) -> Path:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
            for alias, (digest, _) in entries.items():
                tmp.write(f"{alias}\t{digest}\n")
            return Path(tmp.name)

    def flush(self) -> None:
        started = time.monotonic()
        sequence_count = 0
        collection_count = 0
        for namespace, entries in self.sequence.items():
            if not entries:
                continue
            path = self._write_tsv(entries)
            try:
                loaded = self.store.load_sequence_aliases(namespace, str(path))
            finally:
                path.unlink(missing_ok=True)
            expected = len(entries)
            if loaded != expected:
                raise RuntimeError(
                    f"sequence alias batch mismatch for {namespace}: "
                    f"planned {expected}, loaded {loaded}"
                )
            for _, stats in entries.values():
                stats.sequence_aliases_added += 1
            sequence_count += loaded

        for namespace, entries in self.collection.items():
            if not entries:
                continue
            path = self._write_tsv(entries)
            try:
                loaded = self.store.load_collection_aliases(namespace, str(path))
            finally:
                path.unlink(missing_ok=True)
            expected = len(entries)
            if loaded != expected:
                raise RuntimeError(
                    f"collection alias batch mismatch for {namespace}: "
                    f"planned {expected}, loaded {loaded}"
                )
            for _, stats in entries.values():
                stats.collection_aliases_added += 1
            collection_count += loaded
        logger.info(
            "report alias batch flush: %d sequence, %d collection aliases in %.2fs",
            sequence_count,
            collection_count,
            time.monotonic() - started,
        )


def add_aliases_for_row(
    aliases: ReportAliasAccumulator,
    row: dict[str, str],
    namespace: str,
    name_to_digest: dict[str, str],
    stats: AssemblyStats,
    collection_scoped: bool = True,
) -> None:
    """Add supported aliases from one NCBI assembly-report row.

    Args:
        aliases: Accumulator receiving planned sequence aliases.
        row: Parsed NCBI assembly-report row.
        namespace: Assembly namespace for sequence-name and UCSC aliases.
        name_to_digest: RefSeq accession-to-digest mapping for the source FASTA.
        stats: Mutable assembly counters updated for aliases, skips, and warnings.

    GenBank-only rows and rows absent from the store are skipped. Resolved rows
    receive assembly-scoped aliases plus global ``refseq`` and ``insdc`` aliases.
    """
    refseq_ac = row.get("RefSeq-Accn", "").strip()
    ucsc_name = row.get("UCSC-style-name", "").strip()
    genbank_ac = row.get("GenBank-Accn", "").strip()
    seq_name = row.get("Sequence-Name", "").strip()

    # Scope: we only alias sequences that are actually in the RefSeq dataset.
    # Rows with RefSeq-Accn=na are GenBank-only contigs that the NCBI
    # _genomic.fna.gz (RefSeq view) deliberately omits, so there is no
    # sequence data to alias to. Skip cleanly without a warning.
    if refseq_ac.lower() in NA_VALUES:
        stats.rows_skipped_non_refseq += 1
        return

    digest = name_to_digest.get(refseq_ac)
    if digest is None:
        stats.warnings += 1
        source_description = (
            "associated FASTA collection" if collection_scoped else "completed store"
        )
        logger.warning(
            "no digest for %s in namespace %s (%s does not contain it)",
            refseq_ac,
            namespace,
            source_description,
        )
        return

    stats.rows_resolved += 1
    # The NCBI Sequence-Name (e.g. "1", "MT", "HSCHR1_CTG1_UNLOCALIZED") is the
    # spelling seqrepo uses for its assembly namespaces, distinct from the UCSC
    # "chr1"/"chrM" form. Add it so seqrepo-style lookups (GRCh38:1) resolve.
    if seq_name and seq_name.lower() not in NA_VALUES:
        aliases.add_sequence(namespace, seq_name, digest, stats)
    if ucsc_name and ucsc_name.lower() not in NA_VALUES:
        aliases.add_sequence(namespace, ucsc_name, digest, stats)
    if genbank_ac and genbank_ac.lower() not in NA_VALUES:
        aliases.add_sequence(namespace, genbank_ac, digest, stats)
        aliases.add_sequence("insdc", genbank_ac, digest, stats)
    aliases.add_sequence("refseq", refseq_ac, digest, stats)


def add_collection_aliases(
    aliases: ReportAliasAccumulator,
    report: AssemblyReport,
    collection_digest: str,
    stats: AssemblyStats,
) -> None:
    if report.refseq_assembly_accession:
        aliases.add_collection(
            "refseq", report.refseq_assembly_accession, collection_digest, stats
        )
    if report.genbank_assembly_accession:
        aliases.add_collection(
            "insdc", report.genbank_assembly_accession, collection_digest, stats
        )


def ingest_assembly(
    store: RefgetStore,
    entry: AssemblyConfig,
    download_dir: Path,
    provenance: dict[str, dict] | None = None,
) -> DeferredAssemblyReport:
    """Ingest an assembly FASTA and retain its report for deferred application.

    Args:
        store: On-disk RefgetStore to update.
        entry: Manifest assembly configuration and source URLs.
        download_dir: Root of the mirrored source cache.
        provenance: Optional cache-path-to-collection metadata for the build lock.

    A full assembly FASTA produces a collection and provenance record. Reports
    are deliberately left untouched until every selected source is ingested.
    """
    stats = AssemblyStats(namespace=entry.namespace)
    logger.info("=== %s ===", entry.namespace)
    fasta_path: Path | None = None
    if entry.load_fasta:
        # Already resolved against the manifest's directory by load_config.
        if entry.resolved_fasta_path and entry.resolved_fasta_path.exists():
            fasta_path = entry.resolved_fasta_path
        if fasta_path is None:
            assert entry.fasta_url is not None
            fasta_path = require_prepared_source(
                mirror_cache_path(download_dir, entry.fasta_url)
            )
    report_path = require_prepared_source(
        mirror_cache_path(download_dir, entry.report_url)
    )

    coll_digest: str | None = None
    if fasta_path is not None:
        logger.info("ingesting %s", fasta_path)
        coll_meta, was_new = store.add_sequence_collection_from_fasta(str(fasta_path))
        coll_digest = coll_meta.digest
        stats.collection_digest = coll_digest
        if provenance is not None:
            provenance[str(mirror_cache_path(download_dir, entry.fasta_url))] = {
                "collection_digest": coll_digest,
                "n_sequences": coll_meta.n_sequences,
            }
        logger.info(
            "collection %s (%s, %d sequences)",
            coll_digest,
            "new" if was_new else "existing",
            coll_meta.n_sequences,
        )
    return DeferredAssemblyReport(entry, report_path, stats, coll_digest)


def parse_deferred_assembly_report(context: DeferredAssemblyReport) -> None:
    """Parse a deferred report after all sequence ingestion has completed."""
    logger.info("parsing %s", context.report_path)
    report = parse_assembly_report(context.report_path)
    context.report = report
    stats = context.stats
    stats.rows_total = len(report.rows)
    logger.info(
        "report: %d rows, refseq=%s genbank=%s",
        stats.rows_total,
        report.refseq_assembly_accession,
        report.genbank_assembly_accession,
    )


def required_report_only_accessions(
    contexts: list[DeferredAssemblyReport],
) -> set[str]:
    """Return RefSeq accessions needed by parsed report-only assemblies."""
    required: set[str] = set()
    for context in contexts:
        if context.collection_digest is not None:
            continue
        assert context.report is not None
        for row in context.report.rows:
            refseq_ac = row.get("RefSeq-Accn", "").strip()
            if refseq_ac.lower() not in NA_VALUES:
                required.add(refseq_ac)
    return required


def apply_assembly_report(
    store: RefgetStore,
    context: DeferredAssemblyReport,
    report_only_name_to_digest: dict[str, str],
    aliases: ReportAliasAccumulator,
) -> AssemblyStats:
    """Apply one parsed report, preserving collection-scoped resolution."""
    report = context.report
    assert report is not None
    stats = context.stats
    if context.collection_digest is not None:
        name_to_digest = build_name_to_digest_map(store, context.collection_digest)
    else:
        name_to_digest = report_only_name_to_digest

    for row in report.rows:
        add_aliases_for_row(
            aliases,
            row,
            context.entry.namespace,
            name_to_digest,
            stats,
            collection_scoped=context.collection_digest is not None,
        )

    if context.collection_digest is not None:
        add_collection_aliases(aliases, report, context.collection_digest, stats)

    return stats


def _ingest_serial(
    store: RefgetStore, paths: list[Path], stats: SeqsetStats
) -> list[tuple | None]:
    """Import one file at a time, warning on each failure. Used as a fallback."""
    results: list[tuple | None] = []
    for path in paths:
        try:
            results.append(store.add_sequence_collection_from_fasta(str(path)))
        except Exception as exc:  # noqa: BLE001
            logger.warning("ingest failed for %s: %s", path, exc)
            stats.warnings += 1
            results.append(None)
    return results


def ingest_fastas(
    store: RefgetStore, paths: list[Path], stats: SeqsetStats, jobs: int
) -> list[tuple | None]:
    """Import ``paths`` in one batched call, aligned to the input order.

    Returns one ``(metadata, was_new)`` per input path, or ``None`` where that
    file failed.

    Batching is what makes a full build tractable. gtars persists the global
    sequence index per import call, so importing N files one at a time is
    O(N^2) in store size — measured at 9.9x slower for identical work against a
    32 MB index versus an empty one, and the production index exceeds 90 MB.
    One call per seqset amortizes that persist (3.4x on eight shards) and
    ``jobs`` then parallelizes the decompress/digest work, which is what
    dominates the large Ensembl inputs (3.0x on ``dna.toplevel`` at jobs=3).

    Ordering is preserved and the resulting store is byte-identical to serial
    ingestion: sequences are content-addressed, so collection digests do not
    depend on import order. gtars documents ``jobs`` as affecting neither
    ordering nor the resulting store.

    A batch failure aborts the whole call without naming a file, so any failure
    falls back to serial ingestion to attribute the error to a specific input
    and preserve per-file warn-and-continue semantics. Re-importing an
    already-added collection is a no-op (``force=False`` skips duplicates), so
    the retry is safe after a partially applied batch.
    """
    if not paths:
        return []
    try:
        # An ImportReport, whose ``collections`` holds the per-file results in
        # expanded-input order; its counters are not needed here.
        results = store.add_sequence_collections_from_fastas(
            [str(p) for p in paths], jobs=jobs
        ).collections
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "batched ingest of %d file(s) failed (%s); retrying serially to "
            "isolate the offending input", len(paths), exc,
        )
        return _ingest_serial(store, paths, stats)
    # gtars de-duplicates its expanded input list. Distinct URLs cannot collide
    # in the mirrored cache, so a length mismatch means an assumption broke;
    # fall back rather than mis-align results to shards.
    if len(results) != len(paths):
        logger.warning(
            "batched ingest returned %d result(s) for %d input(s); retrying "
            "serially", len(results), len(paths),
        )
        return _ingest_serial(store, paths, stats)
    return list(results)


def process_seqset(
    store: RefgetStore,
    entry: SeqsetConfig,
    download_dir: Path,
    provenance: dict[str, dict] | None = None,
    refresh_derived: bool = False,
    alias_sink: dict[str, str] | None = None,
    jobs: int = INGEST_JOBS_DEFAULT,
) -> SeqsetStats:
    """Ingest a flat/sharded seqset and alias every header name.

    Runs in three phases:
      1. Resolve every shard to a readable FASTA (preflight-prepared source,
         plus any derived-format conversion) before touching the store.
      2. Import them all in ONE batched call — see ``ingest_fastas`` for why
         per-file imports are O(N^2) in store size.
      3. Walk each collection's level-2 contents; collect ``(name, digest)``.

    Aliases are written in bulk via ``store.load_sequence_aliases`` after all
    shards have been ingested. ``add_sequence_alias`` rewrites the full
    namespace TSV on every call (gtars AliasManager), so inserting N entries
    one at a time is O(N^2) — 21k RNA aliases per shard turns into hours.
    ``load_sequence_aliases`` merges a whole TSV into the in-memory alias
    map and triggers exactly one persist, dropping the per-alias cost from
    ~7.5 ms to ~5 µs.

    The seqset is the batching unit deliberately. Batching more widely (across
    releases, or the whole manifest) measured only ~4% faster, while
    ``process_release_groups`` depends on releases being ingested and published
    oldest-first so the rolling namespace ends up holding exactly the newest
    release.
    """
    stats = SeqsetStats(name=entry.name, namespace=entry.namespace)
    logger.info("=== seqset %s ===", entry.name)

    # alias -> digest. Repeated identical mappings are harmless; a conflicting
    # mapping is rejected so immutable namespaces never depend on source order.
    pending: dict[str, str] = {}

    # Phase 1: resolve every shard before importing any of them.
    resolved: list[tuple[Path, Path]] = []  # (source cache path, FASTA to ingest)
    for index, (shard_label, url) in enumerate(entry.iter_shard_urls()):
        target = mirror_cache_path(download_dir, url)
        logger.info(
            "shard %s%s",
            shard_label or "-",
            f" ({Path(url).name})" if shard_label else "",
        )
        try:
            fasta_path = require_prepared_source(target)
        except Exception as exc:  # noqa: BLE001
            logger.warning("prepared source unavailable for %s: %s", url, exc)
            stats.warnings += 1
            continue

        resolver = DERIVED_FASTA_RESOLVERS.get(entry.format)
        if resolver is not None:
            try:
                fasta_path = resolver(fasta_path, refresh_derived)
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "%s conversion failed for %s: %s", entry.format, fasta_path, exc
                )
                stats.warnings += 1
                continue

        file_class = entry.file_class_for_index(index)
        if entry.exclusion is not None and entry.exclusion.applies_to(file_class):
            unfiltered = fasta_path
            try:
                fasta_path = resolve_filtered_fasta(
                    unfiltered, entry.exclusion, refresh_derived
                )
            except Exception as exc:  # noqa: BLE001
                # Naming the file class and URL matters: a seqset can declare
                # exclusions for several classes, and "it failed somewhere in
                # this release" is not actionable.
                logger.warning(
                    "record exclusion failed for %s (file_class=%s, url=%s): %s",
                    fasta_path, file_class, url, exc,
                )
                stats.warnings += 1
                continue
            stats.records_excluded += _excluded_count(unfiltered)

        resolved.append((target, fasta_path))

    # Phase 2: one batched import for the whole seqset.
    if resolved:
        logger.info(
            "  importing %d file(s) (jobs=%d)", len(resolved), jobs
        )
    results = ingest_fastas(store, [f for _, f in resolved], stats, jobs)

    # Phase 3: attribute each per-file result back to its shard.
    for (target, fasta_path), result in zip(resolved, results):
        if result is None:
            continue
        coll_meta, was_new = result

        logger.info(
            "  collection %s (%s, %d sequences)",
            coll_meta.digest,
            "new" if was_new else "existing",
            coll_meta.n_sequences,
        )
        stats.shards_processed += 1
        stats.sequences_ingested += coll_meta.n_sequences
        if provenance is not None:
            # key by the source cache path (the downloaded artifact for derived
            # formats, not the converted .fasta) so it matches what
            # sources.resolve_sources enumerates at lock time
            provenance[str(target)] = {
                "collection_digest": coll_meta.digest,
                "n_sequences": coll_meta.n_sequences,
            }

        name_to_digest = build_name_to_digest_map(store, coll_meta.digest)
        expand = SEQSET_ALIAS_EXPANDERS.get(entry.format, lambda name: (name,))
        for name, digest in name_to_digest.items():
            for alias in expand(name):
                previous = pending.get(alias)
                if previous is not None and previous != digest:
                    raise ValueError(
                        f"seqset {entry.name!r}: alias {alias!r} maps to both "
                        f"{previous} and {digest}"
                    )
                pending[alias] = digest

    if alias_sink is not None:
        for alias, digest in pending.items():
            previous = alias_sink.get(alias)
            if previous is not None and previous != digest:
                raise ValueError(
                    f"release {entry.release}: alias {alias!r} maps to both "
                    f"{previous} and {digest}"
                )
            alias_sink[alias] = digest
        return stats

    if pending:
        # Drop anything already aliased in this namespace (may be left over
        # from prior runs or overlap with other seqsets pointing at the same
        # namespace). list_sequence_aliases is an O(n) in-memory scan.
        existing = set(store.list_sequence_aliases(entry.namespace) or [])
        new_aliases = {a: d for a, d in pending.items() if a not in existing}
        stats.sequence_aliases_skipped += len(pending) - len(new_aliases)

        if new_aliases:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=False
            ) as tmp:
                for alias, digest in new_aliases.items():
                    tmp.write(f"{alias}\t{digest}\n")
                tmp_path = Path(tmp.name)
            try:
                count = store.load_sequence_aliases(entry.namespace, str(tmp_path))
                stats.sequence_aliases_added += count
                logger.info(
                    "  %s: bulk-loaded %d new aliases into ns=%s",
                    entry.name,
                    count,
                    entry.namespace,
                )
            finally:
                tmp_path.unlink(missing_ok=True)

    return stats


def _write_alias_tsv(path: Path, aliases: dict[str, str]) -> None:
    """Atomically write a complete, deterministically ordered alias namespace."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=path.parent, delete=False
    ) as tmp:
        for alias, digest in sorted(aliases.items()):
            tmp.write(f"{alias}\t{digest}\n")
        tmp_path = Path(tmp.name)
    tmp_path.replace(path)


def load_immutable_aliases(
    store: RefgetStore, namespace: str, aliases: dict[str, str]
) -> tuple[int, int]:
    """Load immutable aliases, rejecting an existing conflicting mapping."""
    existing = set(store.list_sequence_aliases(namespace) or [])
    additions: dict[str, str] = {}
    for alias, digest in aliases.items():
        if alias not in existing:
            additions[alias] = digest
            continue
        metadata = resolve_present_alias(store, namespace, alias)
        if metadata.sha512t24u != digest:
            raise ValueError(
                f"immutable alias collision {namespace}:{alias}: "
                f"{metadata.sha512t24u} != {digest}"
            )
    if not additions:
        return 0, len(aliases)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
        for alias, digest in sorted(additions.items()):
            tmp.write(f"{alias}\t{digest}\n")
        tmp_path = Path(tmp.name)
    try:
        loaded = store.load_sequence_aliases(namespace, str(tmp_path))
    finally:
        tmp_path.unlink(missing_ok=True)
    if loaded != len(additions):
        raise RuntimeError(
            f"immutable alias batch mismatch for {namespace}: "
            f"planned {len(additions)}, loaded {loaded}"
        )
    return loaded, len(aliases) - loaded


def process_release_groups(
    store: RefgetStore,
    entries: list[SeqsetConfig],
    download_dir: Path,
    store_dir: Path,
    provenance: dict[str, dict] | None = None,
    refresh_derived: bool = False,
    jobs: int = INGEST_JOBS_DEFAULT,
) -> list[SeqsetStats]:
    """Ingest provider releases oldest-first and publish complete alias snapshots.

    Every source in a release is accumulated before aliases become visible.
    ``provider-N`` namespaces are immutable. The rolling provider namespace is
    atomically replaced with exactly the newest processed release, so retired
    names and changed mappings cannot leak forward.
    """
    grouped: dict[str, dict[int, list[SeqsetConfig]]] = {}
    for entry in entries:
        if entry.release is None or entry.rolling_namespace is None:
            raise ValueError(f"seqset {entry.name!r} is not release-scoped")
        grouped.setdefault(entry.rolling_namespace, {}).setdefault(
            entry.release, []
        ).append(entry)

    stats: list[SeqsetStats] = []
    for rolling_namespace in sorted(grouped):
        list_namespaces = getattr(store, "list_sequence_alias_namespaces", lambda: [])
        registered_namespaces = set(list_namespaces())
        for release in sorted(grouped[rolling_namespace]):
            release_entries = grouped[rolling_namespace][release]
            immutable_namespace = f"{rolling_namespace}-{release}"
            if any(entry.namespace != immutable_namespace for entry in release_entries):
                raise ValueError(
                    f"release {release}: all entries must use {immutable_namespace!r}"
                )
            pending: dict[str, str] = {}
            logger.info(
                "=== %s release %d (%d source groups) ===",
                rolling_namespace, release, len(release_entries),
            )
            for entry in release_entries:
                entry_stats = process_seqset(
                    store, entry, download_dir, provenance, refresh_derived,
                    alias_sink=pending, jobs=jobs,
                )
                expected_sources = sum(1 for _ in entry.iter_shard_urls())
                if (entry_stats.warnings
                        or entry_stats.shards_processed != expected_sources):
                    raise RuntimeError(
                        f"{immutable_namespace}: release source group {entry.name!r} "
                        f"processed {entry_stats.shards_processed}/{expected_sources} "
                        f"files with {entry_stats.warnings} warning(s); aliases not published"
                    )
                stats.append(entry_stats)
            added, skipped = load_immutable_aliases(
                store, immutable_namespace, pending
            )
            logger.info(
                "  %s: loaded=%d existing=%d", immutable_namespace, added, skipped
            )
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", rolling_namespace):
                raise ValueError(f"unsafe rolling namespace: {rolling_namespace!r}")
            rolling_path = (
                store_dir / "aliases" / "sequences" / f"{rolling_namespace}.tsv"
            )
            _write_alias_tsv(rolling_path, pending)
            if rolling_namespace not in registered_namespaces:
                loaded = store.load_sequence_aliases(
                    rolling_namespace, str(rolling_path)
                )
                if loaded != len(pending):
                    raise RuntimeError(
                        f"rolling alias registration mismatch for "
                        f"{rolling_namespace}: planned {len(pending)}, loaded {loaded}"
                    )
                registered_namespaces.add(rolling_namespace)
            logger.info(
                "  %s: replaced complete namespace with %d aliases",
                rolling_namespace, len(pending),
            )
    return stats


def print_summary(
    all_stats: list[AssemblyStats],
    seqset_stats: list[SeqsetStats],
    store: RefgetStore,
) -> None:
    print()
    print("Summary")
    print("-------")
    for s in all_stats:
        print(
            f"  {s.namespace:14s}  rows={s.rows_total:4d} "
            f"resolved={s.rows_resolved:4d} "
            f"non_refseq={s.rows_skipped_non_refseq:3d} "
            f"seq_aliases_added={s.sequence_aliases_added:4d} "
            f"skipped={s.sequence_aliases_skipped:4d} "
            f"coll_aliases_added={s.collection_aliases_added:2d} "
            f"warnings={s.warnings:3d}"
        )
    for s in seqset_stats:
        print(
            f"  {s.name:24s}  ns={s.namespace:8s} "
            f"shards={s.shards_processed:3d} "
            f"seqs={s.sequences_ingested:7d} "
            f"seq_aliases_added={s.sequence_aliases_added:7d} "
            f"skipped={s.sequence_aliases_skipped:6d} "
            f"excluded={s.records_excluded:4d} "
            f"warnings={s.warnings:3d}"
        )
    print()
    print("Store stats:", store.stats())
    print(
        "Sequence alias namespaces:",
        sorted(store.list_sequence_alias_namespaces()),
    )
    print(
        "Collection alias namespaces:",
        sorted(store.list_collection_alias_namespaces()),
    )


def report_lock_check(check, mode: str) -> None:
    """Log the pre-flight lock-check outcome, emphasis per mode."""
    if check.changed:
        logger.warning(
            "lock check: %d file(s) CHANGED vs lock (content drift):", len(check.changed)
        )
        for rel, lsha, asha, mut in check.changed[:20]:
            logger.warning(
                "  %s%s  lock=%s… actual=%s…",
                rel, " (mutable)" if mut else "", lsha[:12], asha[:12],
            )
    if check.missing:
        logger.warning(
            "lock check: %d file(s) not in cache (will be fetched, unchecked): %s",
            len(check.missing), ", ".join(check.missing[:5]) + (" …" if len(check.missing) > 5 else ""),
        )
    if check.set_differs:
        emit = logger.warning if mode == "strict" else logger.info
        emit(
            "lock check: file set differs (in_build_not_lock=%d, in_lock_not_build=%d)",
            len(check.new), len(check.only_in_lock),
        )
    logger.info(
        "lock check: matched=%d changed=%d new=%d missing=%d",
        len(check.matched), len(check.changed), len(check.new), len(check.missing),
    )


def iter_selected_assemblies(
    entries: list[AssemblyConfig], selected: str | None
) -> Iterator[AssemblyConfig]:
    if selected is None:
        yield from entries
        return
    found = False
    for entry in entries:
        if entry.namespace == selected:
            found = True
            yield entry
    if not found:
        raise SystemExit(f"no assembly with namespace {selected!r} in config")


def iter_selected_seqsets(
    entries: list[SeqsetConfig], selected: str | None
) -> Iterator[SeqsetConfig]:
    if selected is None:
        yield from entries
        return
    found = False
    for entry in entries:
        if entry.name == selected:
            found = True
            yield entry
    if not found:
        raise SystemExit(f"no seqset with name {selected!r} in config")


def run_build(args) -> int:
    """Execute the ``build`` subcommand and conditionally refresh its lock.

    Args:
        args: Parsed CLI options controlling source selection, cache/store paths,
            downloads, and build-lock policy.

    Full builds replace the lock with all manifest sources; selected builds merge
    only touched sources. Content drift or a strict file-set mismatch suppresses
    a lock write unless the caller explicitly forces a re-baseline.
    """
    fs_checks.preflight_store_dir(args)
    assemblies, seqsets = load_config(args.config)
    logger.info(
        "loaded %d assembly + %d seqset entries from %s",
        len(assemblies),
        len(seqsets),
        args.config,
    )

    args.cache_dir.mkdir(parents=True, exist_ok=True)

    existing_lock = build_lock.load_lock(args.lock) if args.lock.exists() else None
    if args.locked_sources:
        if existing_lock is None:
            raise SystemExit("--locked-sources requires an existing --lock")
        all_sources = build_lock.apply_locked_sources(seqsets, existing_lock)
    else:
        all_sources = resolve_sources(assemblies, seqsets)

    # If the user passes --assembly or --seqset, they implicitly only want
    # that kind — otherwise we run both.
    run_assemblies = not args.skip_assemblies and args.seqset is None
    run_seqsets = not args.skip_seqsets and args.assembly is None

    sel_assemblies = (
        list(iter_selected_assemblies(assemblies, args.assembly)) if run_assemblies else []
    )
    sel_seqsets = (
        list(iter_selected_seqsets(seqsets, args.seqset)) if run_seqsets else []
    )
    # The concrete files this run will touch (used for both the pre-flight check
    # and the partial add-to merge).
    assembly_owners = {entry.namespace for entry in sel_assemblies}
    seqset_owners = {entry.name for entry in sel_seqsets}
    run_sources = [source for source in all_sources
                   if (source.kind.startswith("assembly_") and source.owner in assembly_owners)
                   or (source.kind == "seqset" and source.owner in seqset_owners)]

    # Materialize and validate every input before opening or mutating the store.
    locked_by_url = ({record["url"]: record
                      for record in build_lock.lock_files(existing_lock)}
                     if existing_lock else {})
    for source in run_sources:
        target = mirror_cache_path(args.cache_dir, source.url)
        locked = locked_by_url.get(source.url)
        if args.locked_sources:
            if not locked or not target.exists() or target.stat().st_size == 0:
                raise SystemExit(f"locked source missing from cache: {source.url}")
            if sha256_file(target) != locked.get("sha256"):
                raise SystemExit(f"locked source SHA-256 mismatch: {source.url}")
        else:
            # The lock's sha256 lets ensure_download tell "the provider
            # republished" apart from "the cache is damaged", and keep the
            # pinned bytes in the first case instead of overwriting them.
            ensure_download(
                source.url, target, args.force_download, source.upstream_md5,
                source, locked_sha256=locked.get("sha256") if locked else None,
            )

    # Second preflight: every declared record exclusion, in parallel. Downloads
    # above and filtration here are both "make the inputs ready before touching
    # the store"; keeping filtration out of the per-seqset loop is what lets the
    # whole set run concurrently.
    prepare_filtered_sources(
        sel_seqsets, args.cache_dir,
        getattr(args, "filter_jobs", INGEST_JOBS_DEFAULT), args.force_download,
    )

    check: build_lock.LockCheck | None = None
    membership_check: build_lock.LockCheck | None = None
    if existing_lock is not None and args.lock_check_mode != "ignore":
        membership_check = build_lock.evaluate_sources_vs_lock(run_sources, existing_lock)
        build_files = [(str(mirror_cache_path(args.cache_dir, source.url)
                            .relative_to(args.cache_dir)),
                        mirror_cache_path(args.cache_dir, source.url))
                       for source in run_sources]
        check = build_lock.evaluate_cache_vs_lock(build_files, existing_lock)
        fatal = (membership_check.drift or check.changed or check.new or check.missing)
        if fatal and not args.force_lock:
            details = []
            if membership_check.set_differs:
                details.append("resolved URL membership changed")
            if membership_check.upstream_changed:
                details.append("upstream MD5 changed")
            if check.changed:
                details.append("cached SHA-256 changed")
            if check.new:
                details.append("lock has no cached SHA-256 baseline")
            if check.missing:
                details.append("resolved input is missing from cache")
            raise SystemExit(
                "source drift detected before ingestion: " + "; ".join(details)
            )

    check_free_space(args.cache_dir, args.store_dir, run_sources, args.min_free_gb)

    args.store_dir.mkdir(parents=True, exist_ok=True)
    store = RefgetStore.on_disk(str(args.store_dir))
    logger.info("opened store at %s (mode=%s)", args.store_dir, store.storage_mode)

    all_stats: list[AssemblyStats] = []
    seqset_stats: list[SeqsetStats] = []
    provenance: dict[str, dict] = {}
    deferred_reports: list[DeferredAssemblyReport] = []

    # Phase 1: make every selected sequence available before applying reports.
    for entry in sel_assemblies:
        deferred_reports.append(
            ingest_assembly(store, entry, args.cache_dir, provenance)
        )
    ordinary_seqsets = [entry for entry in sel_seqsets if entry.release is None]
    release_seqsets = [entry for entry in sel_seqsets if entry.release is not None]
    ordinary_aliases: dict[str, dict[str, str]] = {}
    for entry in ordinary_seqsets:
        entry_stats = process_seqset(
            store, entry, args.cache_dir, provenance,
            refresh_derived=args.force_download,
            alias_sink=ordinary_aliases.setdefault(entry.namespace, {}),
            jobs=args.ingest_jobs,
        )
        if entry_stats.warnings:
            raise RuntimeError(
                f"seqset {entry.name!r} completed with "
                f"{entry_stats.warnings} warning(s); refusing partial publication"
            )
        seqset_stats.append(entry_stats)
    for namespace, pending in ordinary_aliases.items():
        added, skipped = load_immutable_aliases(store, namespace, pending)
        logger.info(
            "immutable seqset namespace %s: loaded=%d existing=%d",
            namespace, added, skipped,
        )
    if release_seqsets:
        seqset_stats.extend(process_release_groups(
            store, release_seqsets, args.cache_dir, args.store_dir, provenance,
            refresh_derived=args.force_download,
            jobs=args.ingest_jobs,
        ))
        # The rolling namespace is replaced as a complete on-disk snapshot.
        # Reopen so subsequent summaries and lock metadata see that snapshot,
        # rather than the alias manager state from before its atomic replacement.
        store.write()
        store = RefgetStore.open_local(str(args.store_dir))
        store.pull_aliases()

    # Phase 2: apply reports only after every sequence source is ingested.
    # Report-only patch releases may reference exact accessions supplied by a
    # later assembly FASTA. For example, the GRCh37.p11/p12 reports reference
    # NW_003871055.3, whose sequence is supplied by the later GRCh37.p13 FASTA;
    # applying those reports in manifest order would omit their scoped aliases.
    for context in deferred_reports:
        parse_deferred_assembly_report(context)
    collection_reports = [
        context for context in deferred_reports
        if context.collection_digest is not None
    ]
    report_only_reports = [
        context for context in deferred_reports
        if context.collection_digest is None
    ]

    aliases = ReportAliasAccumulator(store)
    planning_started = time.monotonic()

    # Collection-backed reports stay strictly scoped to their associated FASTA.
    # Plan all of them first so later assemblies (notably GRCh37.p13) establish
    # pending RefSeq aliases that earlier report-only patches can reuse.
    for context in collection_reports:
        apply_assembly_report(store, context, {}, aliases)

    required_names = required_report_only_accessions(report_only_reports)
    report_only_name_to_digest = build_indexed_name_to_digest_map(
        store, required_names, aliases
    )
    for context in report_only_reports:
        apply_assembly_report(
            store, context, report_only_name_to_digest, aliases
        )
    logger.info(
        "assembly report planning completed in %.2fs", time.monotonic() - planning_started
    )
    aliases.flush()

    # Preserve manifest ordering in the human-readable summary.
    all_stats.extend(context.stats for context in deferred_reports)

    logger.info("persisting store to %s", args.store_dir)
    store.write()
    print_summary(all_stats, seqset_stats, store)

    # ---- lock write decision ----
    full_build = (
        args.assembly is None
        and args.seqset is None
        and not args.skip_assemblies
        and not args.skip_seqsets
    )
    # The write gate keys on CONTENT DRIFT (never silently overwrite a good hash);
    # under strict, a differing file SET also suppresses the write. Both are
    # overridable with --force-lock.
    if args.no_lock:
        logger.info("skipping build-lock write (--no-lock)")
    elif full_build:
        logger.info("writing build-lock (hashing cache) -> %s", args.lock)
        lock = build_lock.build_lock_dict(
            config_path=args.config, download_dir=args.cache_dir,
            resolved_sources=all_sources,
            collection_by_cachepath=provenance, store=store,
            seqsets=seqsets,
        )
        build_lock.validate_lock(lock)
        build_lock.write_lock(args.lock, lock)
        logger.info(
            "wrote build-lock: %d input file(s), %d collection(s), "
            "%d sequence(s)", len(build_lock.lock_files(lock)),
            lock["outputs"]["n_collections"], lock["outputs"]["n_sequences"],
        )
    else:
        logger.info("partial build; merging %d touched source(s) into %s",
                    len(run_sources), args.lock)
        merged = build_lock.merge_into_lock(
            existing_lock, config_path=args.config, download_dir=args.cache_dir,
            touched_sources=run_sources, collection_by_cachepath=provenance, store=store,
            seqsets=seqsets,
        )
        build_lock.validate_lock(merged)
        build_lock.write_lock(args.lock, merged)
        logger.info("build-lock now has %d input file(s)",
                    len(build_lock.lock_files(merged)))
    return 0
