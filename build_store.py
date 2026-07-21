#!/usr/bin/env python
"""Build a custom gtars RefgetStore with vrs-python alias parity.

Ingests NCBI genomic FASTA files into a local RefgetStore and augments it with
per-sequence and per-collection aliases parsed from the matching NCBI assembly
reports. See README.md for details.
"""

from __future__ import annotations

import gzip
import logging
import tempfile
import tomllib
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator
from urllib.parse import urlsplit

from gtars.refget import RefgetStore

logger = logging.getLogger("build_store")

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(__file__).parent / "sources.toml"
DEFAULT_STORE = Path(__file__).parent / "store"
DEFAULT_DOWNLOADS = Path(__file__).parent / "downloads"

SQ_PREFIX = "SQ."
NA_VALUES = {"na", "", "<NA>"}


@dataclass
class AssemblyConfig:
    namespace: str
    report_url: str
    fasta_url: str | None = None
    load_fasta: bool = True
    fasta_path: str | None = None

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
    format: str = "fasta"

    def __post_init__(self) -> None:
        if (self.url_template is None) == (self.urls is None):
            raise ValueError(
                f"seqset {self.name!r}: set exactly one of url_template or urls"
            )
        if self.urls is not None and self.shard_range is not None:
            raise ValueError(
                f"seqset {self.name!r}: shard_range is not valid with urls"
            )
        if self.format not in ("fasta", "gbff"):
            raise ValueError(
                f"seqset {self.name!r}: format must be 'fasta' or 'gbff', "
                f"got {self.format!r}"
            )

    def iter_shard_urls(self) -> Iterator[tuple[str, str]]:
        """Yield (shard_label, url) pairs. Shard label is "" for unsharded."""
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
    warnings: int = 0


def load_config(
    path: Path,
) -> tuple[list[AssemblyConfig], list[SeqsetConfig]]:
    with path.open("rb") as fh:
        data = tomllib.load(fh)
    assemblies = [AssemblyConfig(**entry) for entry in data.get("assembly", [])]
    seqsets = [SeqsetConfig(**entry) for entry in data.get("seqset", [])]
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


def iter_source_urls(
    assemblies: list[AssemblyConfig], seqsets: list[SeqsetConfig]
) -> Iterator[tuple[str, str, str]]:
    """Yield ``(kind, owner, url)`` for every remote file the config references.

    ``kind`` is ``seqset``/``assembly_fasta``/``assembly_report``; ``owner`` is
    the seqset name or assembly namespace. Used by ``fetch_sources.py`` to
    pre-populate the cache and by tooling to enumerate the full source set.
    """
    for s in seqsets:
        for _, url in s.iter_shard_urls():
            yield "seqset", s.name, url
    for a in assemblies:
        if a.load_fasta and not a.fasta_path:
            assert a.fasta_url is not None
            yield "assembly_fasta", a.namespace, a.fasta_url
        yield "assembly_report", a.namespace, a.report_url


def ensure_download(url: str, target: Path, force: bool) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        logger.info("using cached %s", target)
        return target
    logger.info("downloading %s -> %s", url, target)
    tmp = target.with_suffix(target.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)  # noqa: S310
    tmp.replace(target)
    return target


def resolve_fasta_source(
    entry: AssemblyConfig, download_dir: Path, force: bool
) -> Path | None:
    if not entry.load_fasta:
        return None
    if entry.fasta_path:
        local = (REPO_ROOT / entry.fasta_path).resolve()
        if local.exists():
            logger.info("using local fasta override %s", local)
            return local
        logger.warning("fasta_path %s not found; falling back to download", local)
    assert entry.fasta_url is not None
    return ensure_download(
        entry.fasta_url,
        mirror_cache_path(download_dir, entry.fasta_url),
        force,
    )


def resolve_report_source(
    entry: AssemblyConfig, download_dir: Path, force: bool
) -> Path:
    return ensure_download(
        entry.report_url,
        mirror_cache_path(download_dir, entry.report_url),
        force,
    )


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


def build_global_name_to_digest_map(store: RefgetStore) -> dict[str, str]:
    """Fallback map built by scanning every sequence in the store.

    Used when `load_fasta=false` so we can still resolve refseq accessions to
    a digest without having just ingested this assembly's collection. Names
    are globally unique-ish (seqrepo jungle store: 1224 collisions in 580k);
    we keep the first digest we see for each name.
    """
    out: dict[str, str] = {}
    for record in store.iter_sequences():
        md = record.metadata
        out.setdefault(md.name, md.sha512t24u)
    return out


def add_unique_alias(
    store: RefgetStore,
    namespace: str,
    alias: str,
    digest: str,
    stats: AssemblyStats,
) -> None:
    existing = store.get_sequence_by_alias(namespace, alias)
    if existing is not None:
        stats.sequence_aliases_skipped += 1
        return
    store.add_sequence_alias(namespace, alias, digest)
    stats.sequence_aliases_added += 1


def add_aliases_for_row(
    store: RefgetStore,
    row: dict[str, str],
    namespace: str,
    name_to_digest: dict[str, str],
    stats: AssemblyStats,
) -> None:
    """Add supported aliases from one NCBI assembly-report row.

    Args:
        store: RefgetStore receiving sequence aliases.
        row: Parsed NCBI assembly-report row.
        namespace: Assembly namespace for sequence-name and UCSC aliases.
        name_to_digest: RefSeq accession-to-digest mapping for the source FASTA.
        stats: Mutable assembly counters updated for aliases, skips, and warnings.

    GenBank-only rows and rows absent from the FASTA are skipped. Resolved rows
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
        logger.warning(
            "no digest for %s in namespace %s (FASTA does not contain it)",
            refseq_ac,
            namespace,
        )
        return

    stats.rows_resolved += 1
    # The NCBI Sequence-Name (e.g. "1", "MT", "HSCHR1_CTG1_UNLOCALIZED") is the
    # spelling seqrepo uses for its assembly namespaces, distinct from the UCSC
    # "chr1"/"chrM" form. Add it so seqrepo-style lookups (GRCh38:1) resolve.
    if seq_name and seq_name.lower() not in NA_VALUES:
        add_unique_alias(store, namespace, seq_name, digest, stats)
    if ucsc_name and ucsc_name.lower() not in NA_VALUES:
        add_unique_alias(store, namespace, ucsc_name, digest, stats)
    if genbank_ac and genbank_ac.lower() not in NA_VALUES:
        add_unique_alias(store, namespace, genbank_ac, digest, stats)
        add_unique_alias(store, "insdc", genbank_ac, digest, stats)
    add_unique_alias(store, "refseq", refseq_ac, digest, stats)


def add_collection_aliases(
    store: RefgetStore,
    report: AssemblyReport,
    collection_digest: str,
    stats: AssemblyStats,
) -> None:
    if report.refseq_assembly_accession:
        store.add_collection_alias(
            "refseq", report.refseq_assembly_accession, collection_digest
        )
        stats.collection_aliases_added += 1
    if report.genbank_assembly_accession:
        store.add_collection_alias(
            "insdc", report.genbank_assembly_accession, collection_digest
        )
        stats.collection_aliases_added += 1


def process_assembly(
    store: RefgetStore,
    entry: AssemblyConfig,
    download_dir: Path,
    force_download: bool,
    global_name_map_cache: dict[str, dict[str, str]],
    provenance: dict[str, dict] | None = None,
) -> AssemblyStats:
    """Ingest one assembly and materialize its assembly-report aliases.

    Args:
        store: On-disk RefgetStore to update.
        entry: Manifest assembly configuration and source URLs.
        download_dir: Root of the mirrored source cache.
        force_download: Re-fetch cached source files when true.
        global_name_map_cache: Reusable fallback map for report-only assemblies.
        provenance: Optional cache-path-to-collection metadata for the build lock.

    A full assembly FASTA produces a collection and provenance record. A
    report-only entry resolves aliases against the existing store instead.
    """
    stats = AssemblyStats(namespace=entry.namespace)
    logger.info("=== %s ===", entry.namespace)
    fasta_path = resolve_fasta_source(entry, download_dir, force_download)
    report_path = resolve_report_source(entry, download_dir, force_download)

    name_to_digest: dict[str, str]
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
        name_to_digest = build_name_to_digest_map(store, coll_digest)
    else:
        logger.info("load_fasta=false; using global name map")
        cached = global_name_map_cache.get("_global")
        if cached is None:
            cached = build_global_name_to_digest_map(store)
            global_name_map_cache["_global"] = cached
        name_to_digest = cached

    logger.info("parsing %s", report_path)
    report = parse_assembly_report(report_path)
    stats.rows_total = len(report.rows)
    logger.info(
        "report: %d rows, refseq=%s genbank=%s",
        stats.rows_total,
        report.refseq_assembly_accession,
        report.genbank_assembly_accession,
    )

    for row in report.rows:
        add_aliases_for_row(store, row, entry.namespace, name_to_digest, stats)

    if coll_digest is not None:
        add_collection_aliases(store, report, coll_digest, stats)

    return stats


def process_seqset(
    store: RefgetStore,
    entry: SeqsetConfig,
    download_dir: Path,
    force_download: bool,
    provenance: dict[str, dict] | None = None,
) -> SeqsetStats:
    """Ingest a flat/sharded seqset and alias every header name.

    For each shard FASTA:
      1. Download (cached).
      2. ``add_sequence_collection_from_fasta`` -> collection digest.
      3. Walk collection level-2 contents; collect ``(name, digest)`` pairs.

    Aliases are written in bulk via ``store.load_sequence_aliases`` after all
    shards have been ingested. ``add_sequence_alias`` rewrites the full
    namespace TSV on every call (gtars AliasManager), so inserting N entries
    one at a time is O(N^2) — 21k RNA aliases per shard turns into hours.
    ``load_sequence_aliases`` merges a whole TSV into the in-memory alias
    map and triggers exactly one persist, dropping the per-alias cost from
    ~7.5 ms to ~5 µs.
    """
    stats = SeqsetStats(name=entry.name, namespace=entry.namespace)
    logger.info("=== seqset %s ===", entry.name)

    # alias -> digest; later entries win on collision. setdefault() preserves
    # the first digest seen (matches the old add_unique_alias behavior).
    pending: dict[str, str] = {}

    for shard_label, url in entry.iter_shard_urls():
        target = mirror_cache_path(download_dir, url)
        logger.info(
            "shard %s%s",
            shard_label or "-",
            f" ({Path(url).name})" if shard_label else "",
        )
        try:
            fasta_path = ensure_download(url, target, force_download)
        except Exception as exc:  # noqa: BLE001
            logger.warning("download failed for %s: %s", url, exc)
            stats.warnings += 1
            continue

        if entry.format == "gbff":
            try:
                fasta_path = resolve_gbff_fasta(fasta_path, force_download)
            except Exception as exc:  # noqa: BLE001
                logger.warning("gbff conversion failed for %s: %s", fasta_path, exc)
                stats.warnings += 1
                continue

        try:
            coll_meta, was_new = store.add_sequence_collection_from_fasta(
                str(fasta_path)
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("ingest failed for %s: %s", fasta_path, exc)
            stats.warnings += 1
            continue

        logger.info(
            "  collection %s (%s, %d sequences)",
            coll_meta.digest,
            "new" if was_new else "existing",
            coll_meta.n_sequences,
        )
        stats.shards_processed += 1
        stats.sequences_ingested += coll_meta.n_sequences
        if provenance is not None:
            # key by the source cache path (the .gbff.gz for GBFF, not the
            # converted .fasta) so it matches iter_source_urls at lock time
            provenance[str(target)] = {
                "collection_digest": coll_meta.digest,
                "n_sequences": coll_meta.n_sequences,
            }

        name_to_digest = build_name_to_digest_map(store, coll_meta.digest)
        for name, digest in name_to_digest.items():
            pending.setdefault(name, digest)

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
    assemblies, seqsets = load_config(args.config)
    logger.info(
        "loaded %d assembly + %d seqset entries from %s",
        len(assemblies),
        len(seqsets),
        args.config,
    )

    args.store_dir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    store = RefgetStore.on_disk(str(args.store_dir))
    logger.info("opened store at %s (mode=%s)", args.store_dir, store.storage_mode)

    import build_lock  # lazy: avoids a circular import at module load

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
    run_sources = list(iter_source_urls(sel_assemblies, sel_seqsets))

    # ---- pre-flight lock check (before ingesting anything) ----
    existing_lock = build_lock.load_lock(args.lock) if args.lock.exists() else None
    check: build_lock.LockCheck | None = None
    if existing_lock is None:
        logger.info("no existing lock at %s; skipping pre-flight check", args.lock)
    elif args.lock_check_mode == "ignore":
        logger.info("lock-check-mode=ignore; skipping pre-flight check")
    else:
        build_files = [
            (
                str(mirror_cache_path(args.cache_dir, url).relative_to(args.cache_dir)),
                mirror_cache_path(args.cache_dir, url),
            )
            for _, _, url in run_sources
        ]
        check = build_lock.evaluate_cache_vs_lock(build_files, existing_lock)
        report_lock_check(check, args.lock_check_mode)

    all_stats: list[AssemblyStats] = []
    seqset_stats: list[SeqsetStats] = []
    global_name_map_cache: dict[str, dict[str, str]] = {}
    provenance: dict[str, dict] = {}

    for entry in sel_assemblies:
        all_stats.append(
            process_assembly(
                store, entry, args.cache_dir, args.force_download,
                global_name_map_cache, provenance,
            )
        )
    for entry in sel_seqsets:
        seqset_stats.append(
            process_seqset(
                store, entry, args.cache_dir, args.force_download, provenance,
            )
        )

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
    suppress = False
    if check is not None and not args.force_lock:
        if check.drift:
            suppress = True
            logger.warning(
                "content drift vs lock on %d file(s); NOT writing lock "
                "(use --force-lock to re-baseline)", len(check.changed),
            )
        elif args.lock_check_mode == "strict" and check.set_differs:
            suppress = True
            logger.warning(
                "strict: build file set differs from lock; NOT writing lock "
                "(use --lock-check-mode subset, or --force-lock to re-baseline)",
            )

    if args.no_lock:
        logger.info("skipping build-lock write (--no-lock)")
    elif suppress:
        logger.warning("build-lock NOT updated; preserving %s", args.lock)
    elif full_build:
        logger.info("writing build-lock (hashing cache) -> %s", args.lock)
        lock = build_lock.build_lock_dict(
            config_path=args.config, download_dir=args.cache_dir,
            assemblies=assemblies, seqsets=seqsets,
            collection_by_cachepath=provenance, store=store,
        )
        build_lock.write_lock(args.lock, lock)
        logger.info("wrote build-lock with %d sources", len(lock["sources"]))
    else:
        logger.info("partial build; merging %d touched source(s) into %s",
                    len(run_sources), args.lock)
        merged = build_lock.merge_into_lock(
            existing_lock, config_path=args.config, download_dir=args.cache_dir,
            touched_sources=run_sources, collection_by_cachepath=provenance, store=store,
        )
        build_lock.write_lock(args.lock, merged)
        logger.info("build-lock now has %d sources", len(merged["sources"]))
    return 0
