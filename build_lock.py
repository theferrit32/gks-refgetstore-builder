#!/usr/bin/env python
"""Build-lock: a provenance record of exactly what a build consumed.

`build.lock.json` pins, for one build, every concrete source file the build read:
its URL, cached path, byte length, sha256 of the (compressed) bytes, and the
collection digest it produced in the store — plus build metadata (timestamp, git
commit, gtars version, the `sources.toml` hash, store counts).

Two roles:
  * **provenance / reproducibility** — the source of truth for "what build X used".
    `gks-refgetstore verify --lock` re-checks the cache against the pinned sha256s.
  * **file <-> refget-digest mapping** — because each ingested file becomes exactly
    one collection and the store persists collection membership, the lock's
    `url -> collection_digest` is all `provenance.py` needs to answer
    `file -> digests` / `digest -> files` on demand (no giant table stored).

Some sources are mutable upstream. Live builds treat any URL, provider checksum, or
cached-SHA drift as fatal before ingestion unless ``--force-lock`` explicitly
accepts a new baseline. ``--locked-sources`` instead performs no discovery and
requires the exact cached SHA-256 bytes recorded here.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

# build_store imports this module lazily (inside run_build) to avoid a circular
# import, so importing from it at module load is safe here.
import build_store
from build_store import (ResolvedSource, load_config, mirror_cache_path,
                         provider_checksum_matches, resolve_sources)

SCHEMA = "gks-refgetstore-build-lock/3"
V2_SCHEMA = "gks-refgetstore-build-lock/2"
V1_SCHEMA = "gks-refgetstore-build-lock/1"
CHUNK = 1 << 20

# Schemas whose records predate ``ingest_spec``. A missing spec in one of these
# is *unknown*, not *known to be absent* -- the distinction drives sync's
# migration behaviour, and it cannot be recovered from the record itself
# because "no transformation" is also serialized as null.
SCHEMAS_WITHOUT_INGEST_SPEC = (V1_SCHEMA, V2_SCHEMA)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def is_mutable(url: str) -> bool:
    """True for sources NCBI/Ensembl refresh in place (current, non-archived).

    Everything under the immutable archives (genomes/all, historical GBFF, older
    Ensembl releases) is stable; only the "current" endpoints drift."""
    return (
        "/mRNA_Prot/" in url            # current RefSeq transcript/protein shards
        or "/RefSeqGene/" in url         # current RefSeqGene
        or "/release-113/" in url        # pinned-current Ensembl release
    )


def _git_info(repo_root: Path) -> dict:
    def run(*args: str) -> str | None:
        try:
            out = subprocess.run(
                ["git", "-C", str(repo_root), *args],
                capture_output=True, text=True, check=True,
            )
            return out.stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None

    commit = run("rev-parse", "HEAD")
    status = run("status", "--porcelain")
    describe = run("describe", "--tags", "--always", "--dirty")
    return {
        "commit": commit,
        "describe": describe,
        "dirty": bool(status) if status is not None else None,
    }


def _gtars_version() -> str | None:
    try:
        from importlib.metadata import version
        return version("gtars")
    except Exception:  # noqa: BLE001
        return None


def records_from_log(log_path: Path) -> dict[str, dict]:
    """Reconstruct file -> collection from a build log (for backfilling a lock).

    Parses the gtars import line ``Added <digest> (<n> seqs) from <path> in ...``.
    For GBFF the ``<path>`` is the converted ``…gbff.gz.fasta`` sibling; callers
    key by the source cache path, so both the exact path and (if it ends
    ``.fasta``) the stripped path are recorded.
    """
    out: dict[str, dict] = {}
    with log_path.open() as fh:
        for line in fh:
            if not line.startswith("Added "):
                continue
            # Added <digest> (<n> seqs) from <path> in <t>s
            try:
                rest = line[len("Added "):]
                digest, rest = rest.split(" (", 1)
                n_seqs = int(rest.split(" seqs) from ", 1)[0])
                path = rest.split(" seqs) from ", 1)[1].rsplit(" in ", 1)[0].strip()
            except (ValueError, IndexError):
                continue
            rec = {"collection_digest": digest.strip(), "n_sequences": n_seqs}
            out[path] = rec
            # Derived FASTAs are named by appending a suffix to the source
            # artifact, so stripping it recovers the cache path the lock is
            # keyed on. The suffix list lives in build_store beside the
            # resolvers that produce these names; a resolver that invents its
            # own suffix without registering it there would break this mapping
            # silently, leaving collection_digest null for those sources.
            for suffix in build_store.DERIVED_SUFFIXES:
                if path.endswith(suffix):
                    out[path[: -len(suffix)]] = rec
                    break
    return out


def _rel(download_dir: Path, cache_path: Path) -> str:
    try:
        return str(cache_path.relative_to(download_dir))
    except ValueError:
        return str(cache_path)


def canonical_spec_sha256(spec: dict | None) -> str | None:
    """Stable digest of an ingest spec. ``None`` for "no transformation"."""
    if spec is None:
        return None
    canonical = json.dumps(spec, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def ingest_spec_for(
    source: ResolvedSource, seqset_by_name: dict[str, object]
) -> dict | None:
    """Ingest spec for one resolved source, or ``None`` if it has no transform.

    Assembly sources are always ``None``: they are ingested as published.
    """
    if source.kind != "seqset":
        return None
    entry = seqset_by_name.get(source.owner)
    if entry is None:
        return None
    return entry.ingest_spec(source.file_class)


def _seqset_by_name(seqsets) -> dict[str, object]:
    return {entry.name: entry for entry in (seqsets or [])}


def source_record(
    source: ResolvedSource, download_dir: Path,
    collection_by_cachepath: dict[str, dict], hash_files: bool = True,
    seqset_by_name: dict[str, object] | None = None,
) -> dict:
    """Build the lock record for one manifest-referenced remote source.

    ``url`` is mapped to its expected local cache location with
    :func:`build_store.mirror_cache_path`. The resulting record intentionally
    captures the local cache state at lock-generation time, rather than making
    a remote request:

    - ``kind`` identifies ``seqset``, ``assembly_fasta``, or
      ``assembly_report``; ``owner`` identifies the containing manifest entry.
    - ``cache_path`` is relative to ``download_dir`` when possible, while
      ``url`` remains the authoritative remote origin.
    - ``mutable`` is a policy classification derived from the URL. It denotes
      endpoints known to refresh in place; it is not an HTTP immutability
      guarantee.
    - ``present`` is true only when the cache file exists and has nonzero size.
      ``bytes`` and ``sha256`` describe that local file; ``present`` does not
      verify gzip integrity or that the remote object is still available.
    - ``collection_digest`` and ``n_sequences`` come from build-time
      provenance, keyed by the absolute cache path. They describe the gtars
      sequence collection produced from this source, not an individual sequence
      digest or the number of globally new deduplicated sequences in the store.

    Args:
        kind: Source role emitted by :func:`build_store.iter_source_urls`.
        owner: Assembly namespace or seqset name that references ``url``.
        url: Manifest URL for the cached source file.
        download_dir: Root of the mirrored download cache.
        collection_by_cachepath: Build provenance keyed by absolute cache path.
        hash_files: Compute a SHA-256 for present files; disable only for callers
            that intentionally need a metadata-only record.

    Returns:
        A JSON-serializable source entry for ``build.lock.json``. Missing cache
        files retain their identity fields but have null size and hash values.
        Collection provenance is populated independently when supplied by the
        build.
    """
    kind, owner, url = source.kind, source.owner, source.url
    cache_path = mirror_cache_path(download_dir, url)
    rec: dict = {
        "kind": kind, "owner": owner, "url": url,
        "cache_path": _rel(download_dir, cache_path), "mutable": is_mutable(url),
        "upstream_md5": source.upstream_md5,
        "provider_checksum": source.provider_checksum or source.upstream_md5,
        "provider_checksum_algorithm": (
            source.provider_checksum_algorithm
            or ("md5" if source.upstream_md5 else None)
        ),
        "provider_checksum_blocks": source.provider_checksum_blocks,
        "checksum_url": source.checksum_url,
        "file_class": source.file_class,
    }
    if cache_path.exists() and cache_path.stat().st_size > 0:
        rec["bytes"] = cache_path.stat().st_size
        rec["sha256"] = sha256_file(cache_path) if hash_files else None
        rec["present"] = True
    else:
        rec["bytes"] = None
        rec["sha256"] = None
        rec["present"] = False
    coll = collection_by_cachepath.get(str(cache_path))
    rec["collection_digest"] = coll["collection_digest"] if coll else None
    rec["n_sequences"] = coll["n_sequences"] if coll else None
    spec = ingest_spec_for(source, seqset_by_name or {})
    rec["ingest_spec"] = spec
    rec["ingest_spec_sha256"] = canonical_spec_sha256(spec)
    return rec


def _build_meta(config_path: Path, store) -> dict:
    try:
        stats = dict(store.stats())
    except Exception:  # noqa: BLE001
        stats = {}
    config_path = Path(config_path)
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git": _git_info(config_path.resolve().parent),
        "gtars_version": _gtars_version(),
        "sources_toml": {
            "path": str(config_path),
            "sha256": sha256_file(config_path) if config_path.exists() else None,
        },
        "store": {
            "n_sequences": stats.get("n_sequences"),
            "n_collections": stats.get("n_collections"),
        },
    }


def build_lock_dict(
    *,
    config_path: Path,
    download_dir: Path,
    resolved_sources: list[ResolvedSource],
    collection_by_cachepath: dict[str, dict],
    store,
    hash_files: bool = True,
    seqsets=None,
) -> dict:
    """Full-build lock: one entry per source in the whole manifest (replace)."""
    by_name = _seqset_by_name(seqsets)
    sources = [
        source_record(source, download_dir, collection_by_cachepath, hash_files, by_name)
        for source in resolved_sources
    ]
    return {"schema": SCHEMA, "build": _build_meta(config_path, store), "sources": sources}


def merge_into_lock(
    existing_lock: dict | None,
    *,
    config_path: Path,
    download_dir: Path,
    touched_sources: list[ResolvedSource],
    collection_by_cachepath: dict[str, dict],
    store,
    hash_files: bool = True,
    seqsets=None,
    dropped_scopes: set[tuple[str, str]] | None = None,
) -> dict:
    """Replace touched owner scopes while preserving every untouched scope.

    ``dropped_scopes`` are removed outright rather than rewritten -- used by
    ``sync`` when a source disappears from the manifest, so the lock stops
    claiming a collection the store no longer holds.
    """
    by_rel: dict[str, dict] = {}
    by_name = _seqset_by_name(seqsets)
    touched_scopes = {(source.kind, source.owner) for source in touched_sources}
    touched_scopes |= dropped_scopes or set()
    if existing_lock:
        for s in existing_lock.get("sources", []):
            if (s.get("kind"), s.get("owner")) in touched_scopes:
                continue
            by_rel[s["cache_path"]] = s
    for source in touched_sources:
        rec = source_record(
            source, download_dir, collection_by_cachepath, hash_files, by_name
        )
        by_rel[rec["cache_path"]] = rec
    sources = sorted(by_rel.values(), key=lambda r: (r.get("kind", ""), r["cache_path"]))
    return {"schema": SCHEMA, "build": _build_meta(config_path, store), "sources": sources}


@dataclass
class LockCheck:
    """Result of comparing the files a build would ingest against a lock."""
    matched: list[str] = field(default_factory=list)          # in lock, sha256 unchanged
    changed: list[tuple] = field(default_factory=list)        # (rel, lock_sha, actual_sha, mutable)
    new: list[str] = field(default_factory=list)              # present, not in lock (no baseline)
    missing: list[str] = field(default_factory=list)          # not in cache
    only_in_lock: list[str] = field(default_factory=list)     # in lock, not in this build
    upstream_changed: list[tuple[str, str | None, str | None]] = field(default_factory=list)
    legacy_schema: bool = False

    @property
    def drift(self) -> bool:
        return bool(self.changed or self.upstream_changed or self.set_differs or self.legacy_schema)

    @property
    def set_differs(self) -> bool:
        return bool(self.new or self.only_in_lock)


def evaluate_cache_vs_lock(build_files: list[tuple[str, Path]], lock: dict) -> LockCheck:
    """Per-file drift check (check B). ``build_files`` is (rel_cache_path, abs_path)
    for each file the build would ingest. Compares the intersection with the lock
    by sha256; classifies the rest. Does not consider set membership fatal -- that
    is the caller's policy (check C)."""
    lock_by_rel = {s["cache_path"]: s for s in lock.get("sources", [])}
    res = LockCheck()
    seen: set[str] = set()
    for rel, path in build_files:
        seen.add(rel)
        locked = lock_by_rel.get(rel)
        if not path.exists() or path.stat().st_size == 0:
            res.missing.append(rel)
            continue
        if locked and locked.get("sha256"):
            actual = sha256_file(path)
            if actual == locked["sha256"]:
                res.matched.append(rel)
            else:
                res.changed.append((rel, locked["sha256"], actual, bool(locked.get("mutable"))))
        else:
            res.new.append(rel)
    res.only_in_lock = sorted(set(lock_by_rel) - seen)
    return res


def evaluate_sources_vs_lock(sources: list[ResolvedSource], lock: dict) -> LockCheck:
    """Classify URL membership and provider checksums independent of cache state."""
    result = LockCheck(legacy_schema=lock.get("schema") != SCHEMA)
    scope = {(source.kind, source.owner) for source in sources}
    locked = {entry.get("url"): entry for entry in lock.get("sources", [])
              if (entry.get("kind"), entry.get("owner")) in scope}
    live = {source.url: source for source in sources}
    result.new = sorted(set(live) - set(locked))
    result.only_in_lock = sorted(set(locked) - set(live))
    for url in sorted(set(live) & set(locked)):
        old_checksum = (
            locked[url].get("upstream_md5")
            or locked[url].get("provider_checksum")
        )
        new_checksum = live[url].upstream_md5 or live[url].provider_checksum
        if old_checksum != new_checksum:
            result.upstream_changed.append((url, old_checksum, new_checksum))
    return result


@dataclass
class SyncPlan:
    """What ``sync`` would do to make a store match the current manifest.

    Classification is pure: it reads the manifest, the lock, and the cache, and
    touches no store. Every source lands in exactly one bucket.
    """

    unchanged: list[str] = field(default_factory=list)
    added: list[str] = field(default_factory=list)
    removed: list[dict] = field(default_factory=list)
    reingest: list[tuple[str, str]] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    legacy_schema: bool = False
    # Lock records for the reingest set, kept alongside the (rel, reason) pairs
    # because removal needs their collection_digest.
    reingest_records: list[dict] = field(default_factory=list)

    @property
    def touched(self) -> bool:
        return bool(self.added or self.removed or self.reingest)

    @property
    def collection_digests_to_remove(self) -> list[str]:
        """Distinct collections to drop, for ``removed`` plus ``reingest``.

        Distinct because content-identical sources share one collection: the 34
        Ensembl releases carrying padded scaffolds resolve to 8 digests. Nulls
        are dropped -- a source the lock never mapped has nothing to remove.
        """
        digests = {r.get("collection_digest") for r in self.removed}
        digests |= {r.get("collection_digest") for r in self.reingest_records}
        return sorted(d for d in digests if d)


def plan_sync(
    sources: list[ResolvedSource],
    lock: dict,
    download_dir: Path,
    seqsets=None,
    hash_files: bool = True,
) -> SyncPlan:
    """Classify every manifest source against the lock and the cache.

    | class       | condition                                        |
    |-------------|--------------------------------------------------|
    | ``unchanged``| in lock, sha256 matches, ingest spec matches     |
    | ``added``    | not in lock                                      |
    | ``removed``  | in lock, absent from the manifest                |
    | ``reingest`` | sha256 **or** ingest spec moved                  |
    | ``missing``  | referenced by the manifest, absent from cache    |

    A lock written before ``ingest_spec`` existed cannot say whether a
    transformation was applied, so an absent spec is treated as *unknown* and
    resolved conservatively: a source that now declares a transformation is
    re-ingested, one that declares none is left alone. That migrates an older
    lock without special-casing and without re-ingesting the whole manifest.
    """
    by_name = _seqset_by_name(seqsets)
    spec_pinned = lock.get("schema") not in SCHEMAS_WITHOUT_INGEST_SPEC
    lock_by_rel = {s["cache_path"]: s for s in lock.get("sources", [])}
    plan = SyncPlan(legacy_schema=not spec_pinned)
    seen: set[str] = set()

    for source in sources:
        cache_path = mirror_cache_path(download_dir, source.url)
        rel = _rel(download_dir, cache_path)
        seen.add(rel)
        locked = lock_by_rel.get(rel)
        if locked is None:
            plan.added.append(rel)
            continue
        if not cache_path.exists() or cache_path.stat().st_size == 0:
            plan.missing.append(rel)
            continue

        desired = canonical_spec_sha256(ingest_spec_for(source, by_name))
        if hash_files and locked.get("sha256"):
            if sha256_file(cache_path) != locked["sha256"]:
                plan.reingest.append((rel, "upstream bytes changed"))
                plan.reingest_records.append(locked)
                continue
        if spec_pinned:
            if locked.get("ingest_spec_sha256") != desired:
                plan.reingest.append((rel, "ingest spec changed"))
                plan.reingest_records.append(locked)
                continue
        elif desired is not None:
            plan.reingest.append((rel, "ingest spec not pinned by this lock"))
            plan.reingest_records.append(locked)
            continue
        plan.unchanged.append(rel)

    plan.removed = [
        rec for rel, rec in sorted(lock_by_rel.items()) if rel not in seen
    ]
    return plan


def write_lock(path: Path, lock: dict) -> None:
    path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")


def load_lock(path: Path) -> dict:
    lock = json.loads(Path(path).read_text(encoding="utf-8"))
    schema = lock.get("schema")
    if schema not in (SCHEMA, *SCHEMAS_WITHOUT_INGEST_SPEC):
        raise ValueError(f"unsupported build lock schema: {schema!r}")
    return lock


# --------------------------------------------------------------------------- verify

def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}PB"


def gzip_intact(path: Path) -> tuple[bool, str]:
    """``gzip -t`` — validates the decompressed CRC32 + length trailer."""
    proc = subprocess.run(["gzip", "-t", str(path)], capture_output=True, text=True)
    if proc.returncode == 0:
        return True, "gzip ok"
    return False, (proc.stderr.strip() or f"gzip -t exit {proc.returncode}")


def remote_size(url: str, timeout: int) -> int | None:
    """Content-Length from an HTTP HEAD, or None if unavailable."""
    req = urllib.request.Request(
        url, method="HEAD", headers={"User-Agent": "gks-refgetstore-builder"}
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
            cl = resp.headers.get("Content-Length")
            return int(cl) if cl is not None else None
    except (urllib.error.URLError, ValueError, TimeoutError, OSError):
        return None


def _verify_integrity(args) -> int:
    """Check (A): each manifest file present + (for .gz) gzip-intact; optionally
    compare local size to the server's Content-Length."""
    assemblies, seqsets = load_config(args.config)
    sources = resolve_sources(assemblies, seqsets)
    if args.limit:
        sources = sources[: args.limit]

    missing: list[str] = []
    corrupt: list[str] = []
    size_mismatch: list[str] = []
    checksum_mismatch: list[str] = []
    ok = 0
    for i, source in enumerate(sources, 1):
        url = source.url
        target = mirror_cache_path(args.cache_dir, url)
        rel = _rel(args.cache_dir, target)
        if not target.exists() or target.stat().st_size == 0:
            missing.append(rel)
            print(f"[{i}/{len(sources)}] MISSING  {rel}")
            continue
        size = target.stat().st_size
        problems: list[str] = []
        if target.suffix == ".gz":
            good, msg = gzip_intact(target)
            if not good:
                problems.append(f"gzip: {msg}")
                corrupt.append(rel)
        if not provider_checksum_matches(target, source):
            problems.append("provider checksum mismatch")
            checksum_mismatch.append(rel)
        if args.check_remote:
            rsize = remote_size(url, args.timeout)
            if rsize is None:
                problems.append("remote size unavailable")
            elif rsize != size:
                problems.append(f"size local={size} remote={rsize}")
                size_mismatch.append(rel)
        status = "ok" if not problems else "FAIL"
        detail = f"  ({'; '.join(problems)})" if problems else ""
        print(f"[{i}/{len(sources)}] {status:7} {_human(size):>9}  {rel}{detail}")
        if not problems:
            ok += 1

    print(f"\nchecked={len(sources)} ok={ok} missing={len(missing)} "
          f"corrupt={len(corrupt)} checksum_mismatch={len(checksum_mismatch)} "
          f"size_mismatch={len(size_mismatch)}", file=sys.stderr)
    for label, items in (("MISSING", missing), ("CORRUPT", corrupt),
                         ("PROVIDER CHECKSUM MISMATCH", checksum_mismatch),
                         ("SIZE MISMATCH", size_mismatch)):
        if items:
            print(f"{label}:", file=sys.stderr)
            for r in items:
                print(f"  {r}", file=sys.stderr)
    return 1 if (missing or corrupt or checksum_mismatch or size_mismatch) else 0


def _verify_against_lock(args) -> int:
    """Checks (B) per-file drift and (C, with --strict-set) set reproducibility,
    manifest-driven, using the shared evaluator."""
    lock = load_lock(args.lock)
    assemblies, seqsets = load_config(args.config)
    sources = resolve_sources(assemblies, seqsets)
    if args.limit:
        sources = sources[: args.limit]

    build_files: list[tuple[str, Path]] = []
    rel_to_path: dict[str, Path] = {}
    for source in sources:
        url = source.url
        target = mirror_cache_path(args.cache_dir, url)
        rel = _rel(args.cache_dir, target)
        build_files.append((rel, target))
        rel_to_path[rel] = target
    res = evaluate_cache_vs_lock(build_files, lock)

    corrupt: list[str] = []
    new_ok = 0
    for rel in res.new:
        p = rel_to_path[rel]
        if p.suffix == ".gz":
            good, msg = gzip_intact(p)
            if not good:
                corrupt.append(rel)
                print(f"CORRUPT  {rel}  ({msg})")
                continue
        new_ok += 1
    for rel in res.matched:
        print(f"MATCH    {rel}")
    for rel, lsha, asha, mut in res.changed:
        print(f"CHANGED  {rel}{' (mutable)' if mut else ''}  "
              f"lock={lsha[:16]}… actual={asha[:16]}…")

    print(f"\nlock={args.lock} manifest={args.config}", file=sys.stderr)
    print(f"per-file: matched={len(res.matched)} changed={len(res.changed)} "
          f"new={new_ok} missing={len(res.missing)} corrupt={len(corrupt)}", file=sys.stderr)
    print(f"set-diff: in_build_not_lock={len(res.new)} "
          f"in_lock_not_build={len(res.only_in_lock)}", file=sys.stderr)
    for label, items in (
        ("CHANGED (drift/corruption)", [r for r, *_ in res.changed]),
        ("MISSING", res.missing), ("CORRUPT", corrupt),
        ("in lock, not used by this build", res.only_in_lock),
    ):
        if items:
            print(f"{label}:", file=sys.stderr)
            for r in items:
                print(f"  {r}", file=sys.stderr)

    fail = bool(res.missing or res.changed or corrupt)
    if args.strict_set and res.set_differs:
        print("STRICT-SET: build and lock file sets differ", file=sys.stderr)
        fail = True
    return 1 if fail else 0


def run_verify(args) -> int:
    """Execute the ``verify`` subcommand (verify-only; never ingests)."""
    if args.lock:
        return _verify_against_lock(args)
    return _verify_integrity(args)


def run_lock(args) -> int:
    """Execute the ``lock`` subcommand: (re)generate a build.lock.json for a build
    that already ran, reconstructing file->collection from its log."""
    from gtars.refget import RefgetStore  # local import: keeps module light

    assemblies, seqsets = load_config(args.config)
    collection_by_cachepath = records_from_log(args.from_log)
    print(f"parsed {len(collection_by_cachepath)} collection record(s) from "
          f"{args.from_log}", file=sys.stderr)

    store = RefgetStore.open_local(str(args.store_dir))
    lock = build_lock_dict(
        config_path=args.config, download_dir=args.cache_dir,
        resolved_sources=resolve_sources(assemblies, seqsets),
        collection_by_cachepath=collection_by_cachepath, store=store,
        seqsets=seqsets,
    )
    fasta_sources = [s for s in lock["sources"] if s["kind"] != "assembly_report"]
    mapped = [s for s in fasta_sources if s["collection_digest"]]
    print(f"sources={len(lock['sources'])} fasta/gbff={len(fasta_sources)} "
          f"collection-mapped={len(mapped)}", file=sys.stderr)
    write_lock(args.out, lock)
    print(f"wrote {args.out}", file=sys.stderr)
    return 0
