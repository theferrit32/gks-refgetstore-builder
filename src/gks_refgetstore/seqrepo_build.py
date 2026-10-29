"""Build a biocommons SeqRepo instance holding exactly what a RefgetStore holds.

The instance is the SeqRepo side of a like-for-like comparison: the same
sequences, reachable by the same identifiers, as a RefgetStore built from a
``sources.toml`` and its build lock.

Every write goes through the ``biocommons.seqrepo`` library. For each sequence
in the store this module calls ``SeqRepo(..., writeable=True).store(seq,
nsaliases)`` once, with every alias the store holds for that sequence, and
finally ``commit()``. That is the call ``seqrepo load`` makes per FASTA record;
``load`` itself attaches only one alias per record. Seqrepo computes the
``seq_id``, writes and bgzips the FASTA, builds the pysam indexes, inserts the
alias rows, and adds its digest aliases (``VMC``, ``MD5``, ``SHA1``,
``SEGUID``). This module only decides which aliases each sequence gets, and
checks that seqrepo filed each sequence under the store's digest.

Namespaces are passed under the store's own names. Seqrepo's namespace
translation stores ``refseq`` as ``NCBI`` and ``lrg`` as ``LRG``; ``insdc`` and
the assembly namespaces (``GRCh38``, ``GRCh38.p14``, ``NCBI36``, ...) are
stored literally. Lookups translate back, so ``refseq:NC_000001.11`` and
``lrg:LRG_1`` resolve in both backends.

The store must be exactly what the lock pins for the manifest. Selection uses
the build's own mechanism: :func:`build_lock.apply_locked_sources` plus the
manifest's owners. A manifest that is a subset of the lock (``--lock-check-mode
subset``, as for ``sources.genomic.toml`` against the full lock) selects only
its own sources' collections.

The instance is written into a temporary directory beside the target, verified,
and then renamed into place. ``--add-missing`` instead opens an existing
instance writeable, in place, and passes to ``store()`` only the store
sequences it lacks or holds without some of the store's aliases; seqrepo writes
them to a new bgzip file and leaves the existing files as they are. Clone the
instance first (``cp -c -R`` on APFS) to keep the original. Verification checks, in both directions, that the
sequence sets are equal, that every sequence read back from SeqRepo hashes to
its ``seq_id`` and has the store's length and MD5, and that the alias sets are
equal. A provenance record is written beside the instance as
``<instance>.provenance.json``; nothing outside the seqrepo library writes
inside the instance directory.

Requires the ``seqrepo`` extra (``uv sync --extra seqrepo``) and an htslib
``bgzip`` on ``PATH`` (or ``SEQREPO_BGZIP_PATH``), which seqrepo runs when it
commits each FASTA file.
"""

from __future__ import annotations

import json
import logging
import os
import tempfile
import time
from collections import Counter
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path

from . import build_lock, store_census
from .build_store import free_gib
from .sources import load_config, sha256_file

logger = logging.getLogger("seqrepo_build")

PROVENANCE_SCHEMA = "gks-refgetstore-seqrepo/1"
EXTRA_HINT = "uv sync --extra seqrepo"
# Alias namespaces seqrepo's store() adds to every sequence on its own.
DIGEST_NAMESPACES = frozenset({"VMC", "MD5", "SHA1", "SEGUID"})
PROGRESS_EVERY_RESIDUES = 500_000_000
# Free space wanted before writing: bgzip output measured at ~0.23 bytes per
# residue, plus room for the uncompressed FASTA seqrepo writes before each
# commit (up to 1e9 residues).
BYTES_PER_RESIDUE = 0.3
COMMIT_SCRATCH_GIB = 1.5


@dataclass
class SeqRepoLib:
    """The parts of biocommons.seqrepo this module calls, resolved once."""

    SeqRepo: type
    validate_aliases: object
    bgzip: str
    bgzip_version: str


def require_seqrepo() -> SeqRepoLib:
    """Import biocommons.seqrepo and locate bgzip, or stop with what to install."""
    try:
        from biocommons.seqrepo import SeqRepo
        from biocommons.seqrepo.fastadir import fabgz
        from biocommons.seqrepo.utils import validate_aliases
    except ImportError as exc:
        raise SystemExit(
            f"biocommons.seqrepo is not installed ({exc}). It is the optional "
            f"'seqrepo' extra: {EXTRA_HINT}"
        ) from exc
    try:
        exe = fabgz._find_bgzip()
        version = fabgz._get_bgzip_version(exe)
    except (RuntimeError, ValueError) as exc:
        raise SystemExit(
            f"seqrepo needs htslib's bgzip on PATH (or SEQREPO_BGZIP_PATH): {exc}"
        ) from exc
    return SeqRepoLib(SeqRepo, validate_aliases, exe, version)


# ------------------------------------------------------------------ preflight

def expected_collections(config: Path, lock: dict, mode: str) -> dict[str, int]:
    """Collections the lock attributes to the manifest's sources.

    ``strict`` additionally requires the manifest to own every source in the
    lock, i.e. the lock was written for this manifest.
    """
    assemblies, seqsets = load_config(config)
    sources = build_lock.apply_locked_sources(seqsets, lock)
    owned = build_lock.sources_owned_by(sources, assemblies, seqsets)
    if mode == "strict":
        others = sorted({source.owner for source in sources} -
                        {source.owner for source in owned})
        if others:
            raise SystemExit(
                f"the lock pins sources the manifest does not own ({len(others)} "
                f"owner(s), e.g. {others[:3]}); use --lock-check-mode subset for a "
                "manifest that selects part of a larger lock"
            )
    return build_lock.locked_collections_for(lock, owned)


def check_store_scope(store, expected: dict[str, int]) -> list[str]:
    """Problems with the store holding exactly ``expected`` and nothing else."""
    problems = []
    census = store_census.collection_census(store)
    for digest, n in sorted(expected.items()):
        if digest not in census:
            problems.append(f"collection {digest} is pinned by the lock but not in the store")
        elif census[digest] != n:
            problems.append(f"collection {digest} has {census[digest]} sequences; "
                            f"the lock records {n}")
    for digest in sorted(set(census) - set(expected)):
        problems.append(f"collection {digest} is in the store but not pinned for "
                        "this manifest")
    members: set[str] = set()
    for digest in expected:
        if digest in census:
            members.update(store_census.collection_members(store, digest))
    indexed = store_census.sequence_digests(store)
    if indexed != members:
        problems.append(
            f"the store indexes {len(indexed - members)} sequence(s) outside the "
            f"pinned collections and lacks {len(members - indexed)} of their members"
        )
    return problems


# ---------------------------------------------------------------------- write

@dataclass
class WriteStats:
    sequences: int = 0          # sequences new to the instance
    completed: int = 0          # sequences already present that gained aliases
    aliases_added: int = 0      # store aliases added (seqrepo's digest aliases excluded)
    residues: int = 0           # residues passed to store()
    aliases_by_namespace: Counter = field(default_factory=Counter)
    seconds: float = 0.0


def store_aliases(store, digest: str) -> list[dict[str, str]]:
    """The store's sequence aliases for ``digest``, in seqrepo's ``nsaliases`` form."""
    pairs = store.get_aliases_for_sequence(digest) or []
    return [{"namespace": ns, "alias": alias} for ns, alias in sorted(pairs)]


def digests_to_add(store, sr) -> list[str]:
    """Store sequences the instance lacks, or holds without all their aliases."""
    missing = {d for d in store_census.sequence_digests(store) if d not in sr.sequences}
    sr_triples, _, _ = _seqrepo_alias_triples(sr)
    lacking = {d for _, _, d in _store_alias_triples(store) - sr_triples}
    return sorted(missing | lacking)


def write_instance(lib: SeqRepoLib, store, instance_dir: Path, *,
                   add_missing: bool = False) -> WriteStats:
    """Store the sequences of ``store`` in a SeqRepo at ``instance_dir``.

    By default the instance is new and every sequence is stored. With
    ``add_missing``, the instance already exists and only the sequences it
    lacks, or holds without some of the store's aliases, are passed to
    ``store()``, which adds only what is not already there.
    """
    stats = WriteStats()
    started = time.monotonic()
    next_report = PROGRESS_EVERY_RESIDUES
    sr = lib.SeqRepo(str(instance_dir), writeable=True)
    if add_missing:
        digests = digests_to_add(store, sr)
        logger.info("%d store sequence(s) to add or complete", len(digests))
    else:
        digests = sorted(store_census.sequence_digests(store))
    for digest in digests:
        present = add_missing and digest in sr.sequences
        nsaliases = store_aliases(store, digest)
        lib.validate_aliases(nsaliases)
        seq = "".join(store_census.stream_sequence_chunks(store, digest))
        n_added, n_aliases = sr.store(seq, nsaliases)
        # store() files the sequence under seqrepo's own digest of it. For the
        # instance to mirror the store, that must be the store's digest.
        if n_added != (0 if present else 1) or digest not in sr.sequences:
            raise RuntimeError(
                f"seqrepo did not file sequence {digest} under that digest "
                f"(added {n_added}); its normalization (uppercase, '*' removed) "
                "changes this sequence"
            )
        if present:
            stats.completed += 1
        else:
            stats.sequences += 1
        stats.aliases_added += n_aliases
        stats.residues += len(seq)
        stats.aliases_by_namespace.update(a["namespace"] for a in nsaliases)
        del seq
        if stats.residues >= next_report:
            elapsed = time.monotonic() - started
            logger.info("stored %d/%d sequences, %.2f Gbp, %.0fs (%.1f Mbp/s)",
                        stats.sequences, len(digests), stats.residues / 1e9,
                        elapsed, stats.residues / max(elapsed, 1e-9) / 1e6)
            next_report += PROGRESS_EVERY_RESIDUES
    logger.info("committing")
    sr.commit()
    del sr
    stats.seconds = time.monotonic() - started
    logger.info("stored %d new sequence(s) and completed %d, %d residues, %d "
                "alias(es) added, in %.0fs", stats.sequences, stats.completed,
                stats.residues, stats.aliases_added, stats.seconds)
    return stats


# --------------------------------------------------------------------- verify

@dataclass
class VerifyResult:
    problems: list[str] = field(default_factory=list)
    counts: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not self.problems


def _store_alias_triples(store) -> set[tuple[str, str, str]]:
    """``(namespace, alias, digest)`` for every store sequence alias, read forward."""
    out = set()
    for ns in store.list_sequence_alias_namespaces():
        for alias in store.list_sequence_aliases(ns) or []:
            metadata_ = store.get_sequence_metadata_by_alias(ns, alias)
            if metadata_ is None:
                raise RuntimeError(f"store alias {ns}:{alias} resolves to nothing")
            out.add((ns, alias, store_census.strip_sq(metadata_.sha512t24u)))
    return out


def _seqrepo_alias_triples(sr) -> tuple[set[tuple[str, str, str]], Counter, int]:
    """Current seqrepo alias rows, in the namespaces clients query.

    ``find_aliases`` yields each row under its stored namespace and then once
    per translated namespace (``NCBI`` then ``refseq``). The last record of a
    row is the client-facing spelling, which is the store's namespace name.
    Digest namespaces are counted separately.
    """
    by_row: dict[int, dict] = {}
    db_namespace: dict[int, str] = {}
    for record in sr.aliases.find_aliases():
        row = record["seqalias_id"]
        db_namespace.setdefault(row, record["namespace"])
        by_row[row] = record
    triples = set()
    db_counts: Counter = Counter()
    for row, record in by_row.items():
        db_counts[db_namespace[row]] += 1
        if db_namespace[row] in DIGEST_NAMESPACES:
            continue
        triples.add((record["namespace"], record["alias"], record["seq_id"]))
    n_all = len({record["seqalias_id"]
                 for record in sr.aliases.find_aliases(current_only=False)})
    return triples, db_counts, n_all


def verify_instance(lib: SeqRepoLib, store, instance_dir: Path) -> VerifyResult:
    """Prove the instance holds exactly the store's sequences and aliases."""
    result = VerifyResult()
    sr = lib.SeqRepo(str(instance_dir))
    meta = store_census.sequence_metadata_by_digest(store)

    # Sequences: every seqrepo record read back, re-hashed, and matched to the
    # store's length and MD5.
    started = time.monotonic()
    seen: set[str] = set()
    residues = 0
    for record in sr.sequences:
        seq_id, seq = record["seq_id"], record["seq"]
        seen.add(seq_id)
        residues += len(seq)
        if store_census.sha512t24u([seq]) != seq_id:
            result.problems.append(f"sequence {seq_id} read back hashes differently")
        m = meta.get(seq_id)
        if m is None:
            continue
        if not (len(seq) == record["len"] == int(m.length)):
            result.problems.append(f"sequence {seq_id} length {len(seq)}, store "
                                   f"{m.length}")
        md5 = [a["alias"] for a in sr.aliases.find_aliases(seq_id=seq_id, namespace="MD5")]
        if md5 != [m.md5]:
            result.problems.append(f"sequence {seq_id} MD5 {md5}, store {m.md5}")
    only_seqrepo = seen - set(meta)
    only_store = set(meta) - seen
    if only_seqrepo or only_store:
        result.problems.append(
            f"sequence sets differ: {len(only_store)} only in the store "
            f"(e.g. {sorted(only_store)[:3]}), {len(only_seqrepo)} only in seqrepo "
            f"(e.g. {sorted(only_seqrepo)[:3]})"
        )
    result.counts["sequences"] = len(seen)
    result.counts["residues"] = residues
    result.counts["sequence_check_seconds"] = round(time.monotonic() - started, 1)

    # Aliases, both directions, as (client namespace, alias, digest).
    store_triples = _store_alias_triples(store)
    sr_triples, db_counts, n_all = _seqrepo_alias_triples(sr)
    for label, missing in (("store", store_triples - sr_triples),
                           ("seqrepo", sr_triples - store_triples)):
        if missing:
            result.problems.append(
                f"{len(missing)} alias(es) only in the {label}, e.g. "
                f"{sorted(missing)[:3]}"
            )
    n_current = sum(db_counts.values())
    if n_all != n_current:
        result.problems.append(f"{n_all - n_current} non-current alias row(s)")
    for ns in DIGEST_NAMESPACES:
        if db_counts[ns] != len(seen):
            result.problems.append(f"{db_counts[ns]} {ns} aliases for {len(seen)} sequences")
    result.counts["store_aliases"] = len(store_triples)
    result.counts["seqrepo_alias_rows"] = n_current
    result.counts["seqrepo_aliases_by_namespace"] = dict(sorted(db_counts.items()))
    return result


def assembly_name_check(store) -> dict:
    """Compare the store's assembly-name aliases with bioutils' assembly data.

    ``seqrepo add-assembly-names`` aliases each ``=``/``<>`` row of a bioutils
    assembly whose RefSeq accession is stored, by its sequence name and UCSC
    name. This reports, per assembly namespace the store has, whether those
    aliases are all present with the same sequence, and per bioutils assembly
    the store does not have, how many aliases ``add-assembly-names`` would add
    that the store lacks. Informational; nothing here fails the build.
    """
    import bioutils.assemblies

    store_ns = set(store.list_sequence_alias_namespaces())

    def refseq_digest(accession: str) -> str | None:
        m = store.get_sequence_metadata_by_alias("refseq", accession)
        return store_census.strip_sq(m.sha512t24u) if m is not None else None

    shared: dict[str, dict] = {}
    extra: dict[str, int] = {}
    for name, assembly in sorted(bioutils.assemblies.get_assemblies().items()):
        wanted = {}
        for row in assembly["sequences"]:
            if row["relationship"] not in ("=", "<>"):
                continue
            digest = refseq_digest(row["refseq_ac"])
            if digest is None:
                continue
            for alias in [row["name"], *row["aliases"]]:
                wanted[alias] = digest
        if name not in store_ns:
            if wanted:
                extra[name] = len(wanted)
            continue
        missing = differing = 0
        for alias, digest in wanted.items():
            m = store.get_sequence_metadata_by_alias(name, alias)
            if m is None:
                missing += 1
            elif store_census.strip_sq(m.sha512t24u) != digest:
                differing += 1
        shared[name] = {"bioutils_aliases": len(wanted), "missing": missing,
                        "differing": differing}
    return {"shared_namespaces": shared, "not_in_store": extra}


# ----------------------------------------------------------------- provenance

def _version(dist: str) -> str | None:
    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return None


def _instance_files(instance_dir: Path) -> list[dict]:
    return [
        {"path": str(path.relative_to(instance_dir)), "bytes": path.stat().st_size,
         "sha256": sha256_file(path)}
        for path in sorted(p for p in instance_dir.rglob("*") if p.is_file())
    ]


def provenance_record(args, lib: SeqRepoLib, store_dir: Path, instance_dir: Path,
                      expected: dict[str, int], stats: WriteStats | None,
                      verified: VerifyResult | None, names: dict | None) -> dict:
    rgstore = json.loads((store_dir / "rgstore.json").read_text())
    files = _instance_files(instance_dir)
    return {
        "schema": PROVENANCE_SCHEMA,
        "created_utc": datetime.now(UTC).isoformat(),
        "instance": instance_dir.name,
        "config": {"path": str(args.config), "sha256": sha256_file(args.config)},
        "lock": {"path": str(args.lock), "sha256": sha256_file(args.lock),
                 "check_mode": args.lock_check_mode,
                 "collections_selected": len(expected)},
        "store": {"path": str(store_dir.resolve()),
                  **{key: rgstore.get(key) for key in (
                      "sequences_digest", "collections_digest", "aliases_digest",
                      "modified")}},
        "build": {
            "git": build_lock._git_info(Path(args.config).resolve().parent),
            "gtars_version": build_lock._gtars_version(),
            "gtars_source": build_lock._gtars_source(),
            "biocommons.seqrepo": _version("biocommons.seqrepo"),
            "bioutils": _version("bioutils"),
            "pysam": _version("pysam"),
            "bgzip": {"path": lib.bgzip, "version": lib.bgzip_version},
        },
        "write": None if stats is None else {
            "mode": "add-missing" if args.add_missing else "new",
            "sequences": stats.sequences, "completed": stats.completed,
            "aliases_added": stats.aliases_added, "residues": stats.residues,
            "store_aliases_by_namespace": dict(sorted(stats.aliases_by_namespace.items())),
            "seconds": round(stats.seconds, 1),
        },
        "verification": None if verified is None else {
            "ok": verified.ok, "problems": verified.problems[:50],
            "counts": verified.counts,
        },
        "assembly_names_vs_bioutils": names,
        "files": files,
        "bytes": sum(f["bytes"] for f in files),
    }


# ------------------------------------------------------------------------ run

def run_seqrepo(args) -> int:
    """Execute the ``seqrepo`` subcommand."""
    from gtars.refget import RefgetStore

    lib = require_seqrepo()
    if not getattr(args, "verbose", False):
        # seqrepo logs every stored sequence and alias at INFO.
        logging.getLogger("biocommons.seqrepo").setLevel(logging.WARNING)

    store_dir: Path = args.store_dir
    target_dir: Path = args.target_dir
    instance = args.instance_name or store_dir.resolve().name
    final_dir = target_dir / instance
    if args.verify_only and args.add_missing:
        raise SystemExit("--verify-only and --add-missing are mutually exclusive")
    if args.verify_only or args.add_missing:
        if not (final_dir / "aliases.sqlite3").is_file():
            raise SystemExit(f"no SeqRepo instance at {final_dir}")
    elif final_dir.exists():
        raise SystemExit(f"{final_dir} already exists; choose a new --instance-name "
                         "or --target-dir")

    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    lock = build_lock.load_lock(args.lock)
    expected = expected_collections(args.config, lock, args.lock_check_mode)
    problems = check_store_scope(store, expected)
    if problems:
        for problem in problems:
            logger.error("%s", problem)
        raise SystemExit(f"{store_dir} is not the store {args.config} pins in "
                         f"{args.lock}: {len(problems)} problem(s)")
    logger.info("store matches the lock: %d collection(s) pinned for %s",
                len(expected), args.config)

    if args.verify_only:
        verified = verify_instance(lib, store, final_dir)
        _log_verification(verified)
        return 0 if verified.ok else 1

    residues = sum(int(m.length) for m in store.list_sequences())
    target_dir.mkdir(parents=True, exist_ok=True)
    need = residues * BYTES_PER_RESIDUE / 2**30 + COMMIT_SCRATCH_GIB
    free = free_gib(target_dir)
    logger.info("disk preflight: ~%.1f GiB needed, %.1f GiB free", need, free)
    if free < need:
        raise SystemExit(f"only {free:.1f} GiB free in {target_dir}; need ~{need:.1f}")

    if args.add_missing:
        # Written in place: seqrepo adds new bgzip files and alias rows to the
        # existing instance. Clone the instance first to keep the original.
        work = final_dir
        logger.info("adding to %s", final_dir)
    else:
        work = Path(tempfile.mkdtemp(prefix=f".{instance}.", dir=target_dir))
        # mkdtemp creates 0700; the published instance directory is world-readable.
        work.chmod(0o755)
        logger.info("writing %s via %s", final_dir, work)
    stats = write_instance(lib, store, work, add_missing=args.add_missing)

    verified = None
    if not args.no_verify:
        verified = verify_instance(lib, store, work)
        _log_verification(verified)
        if not verified.ok:
            raise SystemExit(f"verification failed; the unverified instance is at "
                             f"{work}")
    names = assembly_name_check(store)
    if work != final_dir:
        os.rename(work, final_dir)
    record = provenance_record(args, lib, store_dir, final_dir, expected, stats,
                               verified, names)
    record_path = target_dir / f"{instance}.provenance.json"
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    logger.info("instance %s: %d files, %.2f GB; provenance %s", final_dir,
                len(record["files"]), record["bytes"] / 1e9, record_path)
    return 0


def _log_verification(verified: VerifyResult) -> None:
    for problem in verified.problems:
        logger.error("verify: %s", problem)
    logger.info("verify: %s %s", "PASS" if verified.ok else "FAIL",
                json.dumps(verified.counts, sort_keys=True))
