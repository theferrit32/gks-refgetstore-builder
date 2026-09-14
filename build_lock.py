#!/usr/bin/env python
"""Build-lock: what a build consumed, and what it produced.

`build.lock.json` has exactly two jobs, and the schema is shaped around them:

* **`inputs`** — upstream reproducibility. Every concrete source file the build
  read: URL, cached path, byte length, SHA-256 of the compressed bytes, the
  provider's own checksum, and the ingest spec that says how those bytes were
  interpreted. Re-pulling these gives the same bytes, or the lock says so.
* **`outputs`** — local completeness. A census of the store the build produced:
  its sequence and collection counts, a canonical root over each digest set,
  and every collection with the source files that contributed to it.

``outputs`` is enumerated **from the store**, not synthesized from build
provenance. Provenance only supplies the ``from`` annotation, so a collection
whose contributor is unknown is recorded honestly as ``"from": []`` rather than
omitted. That distinction is what lets ``verify --store`` be a statement about
the store rather than a restatement of the build's intentions.

``outputs.collections[].from`` is a **list** because a collection can have many
contributing files. Content-addressing means two byte-identical sources produce
one collection, and filtering the Ensembl toplevels collapses 8 collections into
3; in the current lock 43 collections have two or more contributors and one has
twelve. The inverse -- one file to one collection -- does hold.

Some sources are mutable upstream. Live builds treat any URL, provider checksum,
or cached-SHA drift as fatal before ingestion unless ``--force-lock`` explicitly
accepts a new baseline. ``--locked-sources`` instead performs no discovery and
requires the exact cached SHA-256 bytes recorded here.
"""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import store_census
from sources import (DERIVED_SUFFIXES, ResolvedSource, SeqsetConfig,
                     load_config, mirror_cache_path, resolve_sources,
                     sha256_file)

logger = logging.getLogger("build_lock")

SCHEMA = "gks-refgetstore-build-lock/4"

# Every key a v4 input record carries. Spelled out so validate_lock can reject a
# record that lost a field in transit rather than silently reading None for it.
FILE_FIELDS = (
    "kind", "owner", "url", "cache_path", "mutable", "sha256", "bytes",
    "present", "provider_checksum", "provider_checksum_algorithm",
    "provider_checksum_blocks", "checksum_url", "file_class", "ingest_spec",
    "ingest_spec_sha256",
)

__all__ = [  # noqa: RUF022 - grouped by role, not alphabetized
    "SCHEMA", "sha256_file", "is_mutable",
    "lock_files", "lock_collections", "collection_by_file",
    "files_by_collection", "lock_sources_toml_sha256", "upstream_md5_of",
    "validate_lock", "load_lock", "write_lock",
    "build_lock_dict", "merge_into_lock", "store_outputs",
    "apply_locked_sources", "canonical_spec_sha256", "ingest_spec_for",
    "LockCheck", "evaluate_cache_vs_lock", "evaluate_sources_vs_lock",
    "SyncPlan", "plan_sync", "reingest_scope",
]


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
            # keyed on. The suffix list lives in sources.DERIVED_SUFFIXES,
            # alongside the model these paths describe; the resolvers that
            # produce the names are in build_store.DERIVED_FASTA_RESOLVERS. A
            # resolver that invents its own suffix without registering it would
            # break this mapping silently, leaving the collection's ``from``
            # empty.
            for suffix in DERIVED_SUFFIXES:
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


# ----------------------------------------------------------------- accessors

# Every reader goes through these. Raw ``lock["inputs"]["files"]`` subscripting
# is how the old schema's readers ended up silently producing wrong output when
# the shape moved -- two of them bypassed load_lock entirely.

def lock_files(lock: dict) -> list[dict]:
    """Every input file record."""
    return lock.get("inputs", {}).get("files", [])


def lock_collections(lock: dict) -> list[dict]:
    """Every output collection record."""
    return lock.get("outputs", {}).get("collections", [])


def lock_sources_toml_sha256(lock: dict) -> str | None:
    return lock.get("inputs", {}).get("sources_toml_sha256")


def files_by_collection(lock: dict) -> dict[str, list[str]]:
    """``collection digest -> contributing cache paths``, as recorded."""
    return {c["digest"]: list(c.get("from", [])) for c in lock_collections(lock)}


def collection_by_file(lock: dict) -> dict[str, str]:
    """``cache path -> collection digest``.

    The inverse of ``from``, and well-defined: one file yields exactly one
    collection even though one collection can have many files. A file listed
    under two collections is rejected by :func:`validate_lock`.
    """
    out: dict[str, str] = {}
    for collection in lock_collections(lock):
        for cache_path in collection.get("from", []):
            out[cache_path] = collection["digest"]
    return out


def upstream_md5_of(record: dict) -> str | None:
    """The NCBI MD5 for an input record, or ``None``.

    v3 stored this twice -- as ``upstream_md5`` and again as
    ``provider_checksum`` with algorithm ``md5``. v4 keeps one copy and
    reconstructs the other, which is lossless.
    """
    if record.get("provider_checksum_algorithm") == "md5":
        return record.get("provider_checksum")
    return None


# ---------------------------------------------------------------- validation

class LockError(ValueError):
    """A lock that is structurally unusable, as opposed to merely stale."""


def validate_lock(lock: dict) -> dict:
    """Check the invariants a v4 lock must satisfy. Returns it, or raises.

    Checks structure and internal consistency only. It says nothing about
    whether the cache or the store still match -- that is ``verify``'s job.
    """
    schema = lock.get("schema")
    if schema != SCHEMA:
        raise LockError(f"unsupported build lock schema: {schema!r}")
    for section in ("build", "inputs", "outputs"):
        if not isinstance(lock.get(section), dict):
            raise LockError(f"lock is missing its {section!r} section")

    files = lock_files(lock)
    seen_paths: set[str] = set()
    seen_urls: set[str] = set()
    for record in files:
        missing = [key for key in FILE_FIELDS if key not in record]
        if missing:
            raise LockError(
                f"input record {record.get('cache_path')!r} is missing "
                f"{', '.join(missing)}"
            )
        if record["cache_path"] in seen_paths:
            raise LockError(f"duplicate cache_path: {record['cache_path']!r}")
        if record["url"] in seen_urls:
            raise LockError(f"duplicate url: {record['url']!r}")
        seen_paths.add(record["cache_path"])
        seen_urls.add(record["url"])
        expected = canonical_spec_sha256(record["ingest_spec"])
        if record["ingest_spec_sha256"] != expected:
            raise LockError(
                f"ingest_spec_sha256 does not digest ingest_spec for "
                f"{record['cache_path']!r}"
            )

    outputs = lock["outputs"]
    collections = lock_collections(lock)
    digests = [c["digest"] for c in collections]
    if len(digests) != len(set(digests)):
        raise LockError("duplicate collection digest in outputs")
    if outputs.get("n_collections") != len(collections):
        raise LockError(
            f"outputs.n_collections={outputs.get('n_collections')!r} disagrees "
            f"with {len(collections)} enumerated collection(s)"
        )
    expected_root = store_census.digest_root(digests)
    if outputs.get("collections_root") != expected_root:
        raise LockError("outputs.collections_root does not match its collections")
    if not str(outputs.get("sequences_root", "")).startswith("sha256:"):
        raise LockError("outputs.sequences_root is not a sha256 root")

    attributed: set[str] = set()
    for collection in collections:
        for cache_path in collection.get("from", []):
            if cache_path not in seen_paths:
                raise LockError(
                    f"collection {collection['digest']} claims unknown "
                    f"contributor {cache_path!r}"
                )
            if cache_path in attributed:
                raise LockError(
                    f"{cache_path!r} is attributed to more than one collection"
                )
            attributed.add(cache_path)
    return lock


def write_lock(path: Path, lock: dict) -> None:
    Path(path).write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")


def load_lock(path: Path) -> dict:
    """Read and validate a lock. The only supported entry point.

    Pre-v4 locks are rejected outright rather than shimmed: the project is
    unpublished, and a reader that silently accepted a v2 ``sources`` list would
    report an empty store census as a passing verify.
    """
    return validate_lock(json.loads(Path(path).read_text(encoding="utf-8")))


# -------------------------------------------------------------- lock writing

def file_record(
    source: ResolvedSource, download_dir: Path, hash_files: bool = True,
    seqset_by_name: dict[str, object] | None = None,
) -> dict:
    """The ``inputs.files`` record for one manifest-referenced remote source.

    ``url`` is mapped to its expected local cache location with
    :func:`sources.mirror_cache_path`. The record intentionally captures the
    local cache state at lock-generation time rather than making a remote
    request:

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
    - ``ingest_spec`` pins how the bytes were interpreted. ``None`` means
      "ingested as published"; the lock's ``sha256`` pins the bytes read, this
      pins what was done with them.

    Which collection the file produced is recorded on the *collection*, in
    ``outputs.collections[].from``, not here.
    """
    cache_path = mirror_cache_path(download_dir, source.url)
    present = cache_path.exists() and cache_path.stat().st_size > 0
    spec = ingest_spec_for(source, seqset_by_name or {})
    return {
        "kind": source.kind,
        "owner": source.owner,
        "url": source.url,
        "cache_path": _rel(download_dir, cache_path),
        "mutable": is_mutable(source.url),
        "provider_checksum": source.provider_checksum or source.upstream_md5,
        "provider_checksum_algorithm": (
            source.provider_checksum_algorithm
            or ("md5" if source.upstream_md5 else None)
        ),
        "provider_checksum_blocks": source.provider_checksum_blocks,
        "checksum_url": source.checksum_url,
        "file_class": source.file_class,
        "bytes": cache_path.stat().st_size if present else None,
        "sha256": (sha256_file(cache_path) if present and hash_files else None),
        "present": present,
        "ingest_spec": spec,
        "ingest_spec_sha256": canonical_spec_sha256(spec),
    }


def _build_meta(config_path: Path) -> dict:
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git": _git_info(Path(config_path).resolve().parent),
        "gtars_version": _gtars_version(),
    }


def _inputs(config_path: Path, files: list[dict]) -> dict:
    config_path = Path(config_path)
    return {
        "sources_toml_sha256": (
            sha256_file(config_path) if config_path.exists() else None
        ),
        "files": files,
    }


def provenance_from_map(
    download_dir: Path, collection_by_cachepath: dict[str, dict]
) -> dict[str, set[str]]:
    """``collection digest -> {cache path}`` from one build's provenance."""
    out: dict[str, set[str]] = {}
    for abs_path, record in collection_by_cachepath.items():
        digest = record.get("collection_digest")
        if digest:
            out.setdefault(digest, set()).add(_rel(download_dir, Path(abs_path)))
    return out


def store_outputs(
    store, from_map: dict[str, set[str]] | None = None,
    known_paths: set[str] | None = None,
) -> dict:
    """Census the store: counts, digest-set roots, and every collection.

    Enumerated from the store. ``from_map`` only annotates; a collection absent
    from it gets ``"from": []``, which is the truthful record of "this store
    holds a collection no known source file explains".

    ``known_paths`` restricts the annotation to files the lock actually records
    as inputs. Build provenance can name a path that is not one: an assembly
    with a ``fasta_path`` override is ingested from a local file, and
    ``resolve_sources`` deliberately emits no remote source for it, so nothing
    pins it in ``inputs.files``. Listing it under ``from`` would make the lock
    reference an input it does not contain -- which ``validate_lock`` rejects,
    and rightly: ``from`` means "contributing *recorded* input files". The
    honest record for such a collection is an empty ``from``.

    Costs about four seconds on the production store -- roughly 2.5s to
    enumerate 1.8M sequences and 1.5s to sort and hash them.
    """
    from_map = from_map or {}
    if known_paths is not None:
        filtered: dict[str, set[str]] = {}
        dropped: set[str] = set()
        for digest, paths in from_map.items():
            kept = {p for p in paths if p in known_paths}
            dropped |= paths - kept
            filtered[digest] = kept
        if dropped:
            logger.info(
                "%d ingested file(s) are not recorded inputs and are omitted "
                "from collection attribution: %s",
                len(dropped), ", ".join(sorted(dropped)[:3])
                + (" …" if len(dropped) > 3 else ""),
            )
        from_map = filtered
    census = store_census.collection_census(store)
    sequences = store_census.sequence_digests(store)
    return {
        # len() rather than store.stats(), which returns these as strings. A
        # lock recording "1779052" compares unequal to the 1779052 it describes.
        "n_sequences": len(sequences),
        "n_collections": len(census),
        "sequences_root": store_census.digest_root(sequences),
        "collections_root": store_census.digest_root(census),
        "collections": [
            {
                "digest": digest,
                "n_sequences": census[digest],
                "from": sorted(from_map.get(digest, ())),
            }
            for digest in sorted(census)
        ],
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
    files = sorted(
        (file_record(source, download_dir, hash_files, by_name)
         for source in resolved_sources),
        key=lambda r: (r["kind"], r["cache_path"]),
    )
    return {
        "schema": SCHEMA,
        "build": _build_meta(config_path),
        "inputs": _inputs(config_path, files),
        "outputs": store_outputs(
            store, provenance_from_map(download_dir, collection_by_cachepath),
            known_paths={record["cache_path"] for record in files},
        ),
    }


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
    claiming an input the build no longer reads.

    ``outputs`` is re-censused from the store rather than merged, because the
    store is the authority on what it now contains. Only the ``from``
    attribution is carried forward, and only for files that survived the merge.
    """
    by_rel: dict[str, dict] = {}
    by_name = _seqset_by_name(seqsets)
    touched_scopes = {(source.kind, source.owner) for source in touched_sources}
    touched_scopes |= dropped_scopes or set()
    carried: dict[str, set[str]] = {}
    if existing_lock:
        retained = {
            record["cache_path"] for record in lock_files(existing_lock)
            if (record["kind"], record["owner"]) not in touched_scopes
        }
        for record in lock_files(existing_lock):
            if record["cache_path"] in retained:
                by_rel[record["cache_path"]] = record
        for digest, paths in files_by_collection(existing_lock).items():
            kept = {p for p in paths if p in retained}
            if kept:
                carried[digest] = kept
    for source in touched_sources:
        record = file_record(source, download_dir, hash_files, by_name)
        by_rel[record["cache_path"]] = record

    from_map = carried
    for digest, paths in provenance_from_map(
        download_dir, collection_by_cachepath
    ).items():
        from_map.setdefault(digest, set()).update(paths)

    files = sorted(by_rel.values(), key=lambda r: (r.get("kind", ""), r["cache_path"]))
    return {
        "schema": SCHEMA,
        "build": _build_meta(config_path),
        "inputs": _inputs(config_path, files),
        "outputs": store_outputs(
            store, from_map,
            known_paths={record["cache_path"] for record in files},
        ),
    }


# ------------------------------------------------------------ locked sources

def apply_locked_sources(
    seqsets: list[SeqsetConfig], lock: dict
) -> list[ResolvedSource]:
    """Use concrete lock entries without performing live discovery.

    Lives here rather than in ``sources`` so it reads the lock through the
    accessors, and so there is exactly one place that knows which schema carries
    concrete sources. Construction is by keyword: ``ResolvedSource`` has nine
    fields, six of them optional strings, and positional construction from
    ``dict.get`` calls silently reorders into a valid-looking source if either
    side is edited.
    """
    validate_lock(lock)
    sources = [
        ResolvedSource(
            kind=record["kind"],
            owner=record["owner"],
            url=record["url"],
            upstream_md5=upstream_md5_of(record),
            provider_checksum=record.get("provider_checksum"),
            provider_checksum_algorithm=record.get("provider_checksum_algorithm"),
            provider_checksum_blocks=record.get("provider_checksum_blocks"),
            checksum_url=record.get("checksum_url"),
            file_class=record.get("file_class"),
        )
        for record in lock_files(lock)
    ]
    by_owner: dict[str, list[ResolvedSource]] = {}
    for source in sources:
        if source.kind == "seqset":
            by_owner.setdefault(source.owner, []).append(source)
    for seqset in seqsets:
        locked = by_owner.get(seqset.name, [])
        if not locked:
            raise LockError(
                f"lock has no concrete sources for seqset {seqset.name!r}"
            )
        seqset.resolved_sources = locked
    return sources


# ---------------------------------------------------------------- comparison

@dataclass
class LockCheck:
    """Result of comparing the files a build would ingest against a lock."""
    matched: list[str] = field(default_factory=list)          # in lock, sha256 unchanged
    changed: list[tuple] = field(default_factory=list)        # (rel, lock_sha, actual_sha, mutable)
    new: list[str] = field(default_factory=list)              # present, not in lock (no baseline)
    missing: list[str] = field(default_factory=list)          # not in cache
    only_in_lock: list[str] = field(default_factory=list)     # in lock, not in this build
    upstream_changed: list[tuple[str, str | None, str | None]] = field(default_factory=list)

    @property
    def drift(self) -> bool:
        return bool(self.changed or self.upstream_changed or self.set_differs)

    @property
    def set_differs(self) -> bool:
        return bool(self.new or self.only_in_lock)


def evaluate_cache_vs_lock(build_files: list[tuple[str, Path]], lock: dict) -> LockCheck:
    """Per-file drift check. ``build_files`` is (rel_cache_path, abs_path) for
    each file the build would ingest. Compares the intersection with the lock by
    sha256; classifies the rest. Does not consider set membership fatal -- that
    is the caller's policy."""
    lock_by_rel = {record["cache_path"]: record for record in lock_files(lock)}
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
    result = LockCheck()
    scope = {(source.kind, source.owner) for source in sources}
    locked = {record["url"]: record for record in lock_files(lock)
              if (record["kind"], record["owner"]) in scope}
    live = {source.url: source for source in sources}
    result.new = sorted(set(live) - set(locked))
    result.only_in_lock = sorted(set(locked) - set(live))
    for url in sorted(set(live) & set(locked)):
        old_checksum = locked[url].get("provider_checksum")
        new_checksum = live[url].provider_checksum or live[url].upstream_md5
        if old_checksum != new_checksum:
            result.upstream_changed.append((url, old_checksum, new_checksum))
    return result


# --------------------------------------------------------------------- sync

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
    # Lock records for the reingest set, kept alongside the (rel, reason) pairs
    # because removal needs their collection.
    reingest_records: list[dict] = field(default_factory=list)
    # cache_path -> collection digest, from the lock this plan was built against.
    collection_by_file: dict[str, str] = field(default_factory=dict)

    @property
    def touched(self) -> bool:
        return bool(self.added or self.removed or self.reingest)

    @property
    def collection_digests_to_remove(self) -> list[str]:
        """Distinct collections to drop, for ``removed`` plus ``reingest``.

        Distinct because content-identical sources share one collection: the 34
        Ensembl releases carrying padded scaffolds resolve to 8 digests. Files
        the lock never attributed to a collection contribute nothing to remove.
        """
        touched = [record["cache_path"]
                   for record in (*self.removed, *self.reingest_records)]
        return sorted({
            digest for digest in
            (self.collection_by_file.get(rel) for rel in touched) if digest
        })


def plan_sync(
    sources: list[ResolvedSource],
    lock: dict,
    download_dir: Path,
    seqsets=None,
    hash_files: bool = True,
) -> SyncPlan:
    """Classify every manifest source against the lock and the cache.

    | class        | condition                                       |
    |--------------|-------------------------------------------------|
    | ``unchanged``| in lock, sha256 matches, ingest spec matches     |
    | ``added``    | not in lock                                      |
    | ``removed``  | in lock, absent from the manifest                |
    | ``reingest`` | sha256 **or** ingest spec moved                  |
    | ``missing``  | referenced by the manifest, absent from cache    |
    """
    by_name = _seqset_by_name(seqsets)
    lock_by_rel = {record["cache_path"]: record for record in lock_files(lock)}
    plan = SyncPlan(collection_by_file=collection_by_file(lock))
    seen: set[str] = set()

    for source in sources:
        cache_path = mirror_cache_path(download_dir, source.url)
        rel = _rel(download_dir, cache_path)
        seen.add(rel)
        # Cache presence is checked before lock membership: a source in
        # neither -- an upstream file published since the lock was written --
        # is *missing*, not *added*. sync does not download, so classifying it
        # as added would let it past the pre-flight and fail during ingest,
        # after the removals have already been applied.
        if not cache_path.exists() or cache_path.stat().st_size == 0:
            plan.missing.append(rel)
            continue
        locked = lock_by_rel.get(rel)
        if locked is None:
            plan.added.append(rel)
            continue

        desired = canonical_spec_sha256(ingest_spec_for(source, by_name))
        if hash_files and locked.get("sha256"):
            if sha256_file(cache_path) != locked["sha256"]:
                plan.reingest.append((rel, "upstream bytes changed"))
                plan.reingest_records.append(locked)
                continue
        if locked.get("ingest_spec_sha256") != desired:
            plan.reingest.append((rel, "ingest spec changed"))
            plan.reingest_records.append(locked)
            continue
        plan.unchanged.append(rel)

    plan.removed = [
        record for rel, record in sorted(lock_by_rel.items()) if rel not in seen
    ]
    return plan


def reingest_scope(lock: dict, cache_paths: list[str]) -> set[str]:
    """Cache paths that must be re-ingested when ``cache_paths`` are removed.

    Removal is per *collection*, so dropping a collection because one
    contributor changed also drops every other contributor's sequences. The
    re-ingest set is therefore the transitive closure: expand each touched file
    to its collection, then back out to that collection's whole ``from`` list.

    ``sync`` previously approximated this by expanding to every seqset feeding
    the same alias *namespace*. In the committed lock that approximation happens
    to be sufficient -- all 31 collections with contributors in more than one
    seqset are Ensembl releases, and the release group expansion covers them.
    But it is sufficient by coincidence, not by construction: nothing about
    content-addressing respects namespace boundaries, and a manifest edit as
    small as splitting one seqset's URL list in two would produce a shared
    collection the namespace expansion cannot reach. Following ``from`` is
    correct for the same reason in every case.
    """
    by_file = collection_by_file(lock)
    from_map = files_by_collection(lock)
    scope = set(cache_paths)
    for cache_path in cache_paths:
        digest = by_file.get(cache_path)
        if digest:
            scope.update(from_map.get(digest, ()))
    return scope


def owners_for_paths(lock: dict, cache_paths: set[str]) -> set[tuple[str, str]]:
    """``(kind, owner)`` scopes covering the given cache paths."""
    by_rel = {record["cache_path"]: record for record in lock_files(lock)}
    return {
        (by_rel[rel]["kind"], by_rel[rel]["owner"])
        for rel in cache_paths if rel in by_rel
    }


# --------------------------------------------------------------------- lock

def run_lock(args) -> int:
    """Execute the ``lock`` subcommand: (re)generate a build.lock.json for a
    build that already ran, reconstructing file->collection from its log."""
    import sys

    from gtars.refget import RefgetStore  # local import: keeps module light

    assemblies, seqsets = load_config(args.config)
    collection_by_cachepath = records_from_log(args.from_log)
    print(f"parsed {len(collection_by_cachepath)} collection record(s) from "
          f"{args.from_log}", file=sys.stderr)

    store = RefgetStore.open_local(str(args.store_dir))
    store.set_quiet(True)
    lock = build_lock_dict(
        config_path=args.config, download_dir=args.cache_dir,
        resolved_sources=resolve_sources(assemblies, seqsets),
        collection_by_cachepath=collection_by_cachepath, store=store,
        seqsets=seqsets,
    )
    outputs = lock["outputs"]
    attributed = sum(1 for c in outputs["collections"] if c["from"])
    print(f"inputs={len(lock_files(lock))} "
          f"collections={outputs['n_collections']} "
          f"attributed={attributed} sequences={outputs['n_sequences']}",
          file=sys.stderr)
    validate_lock(lock)
    write_lock(args.out, lock)
    print(f"wrote {args.out}", file=sys.stderr)
    return 0
