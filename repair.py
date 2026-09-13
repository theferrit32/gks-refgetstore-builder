#!/usr/bin/env python
"""Restore the state the build lock records, after ``verify`` finds damage.

Repair is deliberately narrow. It puts back **what the lock says should be
there** and nothing else:

* A cache file that is missing, truncated, mis-hashed, or gzip-corrupt is
  re-fetched, and the download is accepted **only if its SHA-256 matches the
  lock**. New upstream bytes are not a repair -- they are a ``sync`` or
  ``--force-lock`` decision, and silently accepting them would make repair a
  re-baselining tool that quietly destroys the drift the lock preserves.
* A collection the store has lost, or whose sequences have lost their payload
  files, is removed and re-ingested from the cache.

Two things it will not do, for reasons worth stating rather than discovering:

* **``store.seq_digest_mismatch`` is unrepairable.** Those sequences decode to
  something other than what their digest promises because gtars' encoder is
  lossy for them -- no ``U`` in the protein alphabet, and short proteins
  misdetected as nucleotide. Re-ingesting reproduces the defect exactly. The fix
  is upstream in gtars, and ``verify --deep``'s baseline is what will report it
  when it lands.
* **A root mismatch with no attributable cause is refused.** If the store's
  digest set differs from the lock's but no collection is missing and no payload
  is absent, the store has diverged in a way repair cannot name. That is what
  ``sync`` is for.

Ordering is not incidental: **cache repair always precedes store repair, and
aborts the whole run if it fails.** Store mutation removes collections before it
re-ingests, removal persists immediately, and gtars offers no rollback -- so
discovering a missing source file after the removal has already happened leaves
the store worse than it started.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import build_lock
import build_store
import fetch_sources
import verify
from build_store import ResolvedSource, sha256_file

# verify codes this module knows how to act on. A code absent from both tables
# is reported and left alone, which is the safe default for a tool that mutates.
CACHE_CODES = frozenset({
    "cache.missing", "cache.sha256_mismatch", "cache.gzip_corrupt",
    "cache.size_mismatch",
})
STORE_CODES = frozenset({
    "store.collection_missing", "store.collection_size_mismatch",
    "store.collection_unresolvable", "store.seq_file_missing",
})
UNREPAIRABLE = {
    "store.seq_digest_mismatch": (
        "the stored payload does not re-digest to its own digest. This is a "
        "gtars encoder limitation, not bit rot: re-ingesting reproduces it "
        "byte for byte. Track it in the --deep baseline; fix it upstream."
    ),
    "store.sequences_root_mismatch": (
        "the store's sequence set is not the lock's, and no missing collection "
        "or payload explains it. Repair restores; it does not reconcile. Use "
        "`gks-refgetstore sync` to make the store match the manifest."
    ),
    "store.collections_root_mismatch": (
        "the store's collection set is not the lock's, and no missing "
        "collection explains it. Use `gks-refgetstore sync`."
    ),
}


@dataclass
class RepairPlan:
    """What repair would do, before it does any of it."""

    cache_files: list[str] = field(default_factory=list)
    store_files: list[str] = field(default_factory=list)
    store_collections: list[str] = field(default_factory=list)
    assembly_owners: set[str] = field(default_factory=set)
    unrepairable: list[tuple[str, str]] = field(default_factory=list)
    ignored: list[tuple[str, str]] = field(default_factory=list)

    @property
    def touched(self) -> bool:
        return bool(self.cache_files or self.store_files)


def locked_source(record: dict) -> ResolvedSource:
    """A ``ResolvedSource`` for one lock record, by keyword.

    Only ``kind``/``owner``/``url`` are load-bearing here -- repair accepts
    bytes on the lock's SHA-256, not on a provider checksum, so the provider
    fields are carried for logging rather than validation.
    """
    return ResolvedSource(
        kind=record["kind"],
        owner=record["owner"],
        url=record["url"],
        upstream_md5=build_lock.upstream_md5_of(record),
        provider_checksum=record.get("provider_checksum"),
        provider_checksum_algorithm=record.get("provider_checksum_algorithm"),
        provider_checksum_blocks=record.get("provider_checksum_blocks"),
        checksum_url=record.get("checksum_url"),
        file_class=record.get("file_class"),
    )


def plan_repair(
    lock: dict, findings: list[verify.Finding], store_dir: Path | None = None,
) -> RepairPlan:
    """Dispatch verify findings onto the actions that address them.

    Store findings are resolved digest -> collection -> ``from`` -> source
    files, which is why the ``from`` list had to become a list: a collection
    removed because one contributor's payload vanished carries *every*
    contributor's sequences, so all of them must be re-ingested.
    """
    plan = RepairPlan()
    by_digest_collection: dict[str, dict] = {
        c["digest"]: c for c in build_lock.lock_collections(lock)
    }
    # digest -> the collection publishing it, for seq_file_missing findings.
    collection_of_sequence: dict[str, str] = {}
    damaged_collections: set[str] = set()

    for finding in findings:
        if finding.code in CACHE_CODES:
            plan.cache_files.append(finding.subject)
        elif finding.code in ("store.collection_missing",
                              "store.collection_size_mismatch",
                              "store.collection_unresolvable"):
            damaged_collections.add(finding.subject)
        elif finding.code == "store.seq_file_missing":
            collection_of_sequence.setdefault(finding.subject, "")
        elif finding.code in UNREPAIRABLE:
            plan.unrepairable.append((finding.code, finding.subject))
        elif finding.severity != verify.INFO:
            plan.ignored.append((finding.code, finding.subject))

    # A missing payload names a sequence, not a collection. Resolve it by asking
    # which locked collections publish that digest -- this is the one place
    # repair has to open the store during planning.
    if collection_of_sequence:
        if store_dir is None:
            # Without the store there is no way from a damaged sequence back to
            # the sources that produced it. Say so rather than silently
            # planning nothing for these findings.
            plan.ignored.extend(
                ("store.seq_file_missing", digest)
                for digest in sorted(collection_of_sequence)
            )
        else:
            damaged_collections |= _collections_publishing(
                lock, store_dir, set(collection_of_sequence)
            )

    plan.store_collections = sorted(damaged_collections)
    contributors: set[str] = set()
    for digest in plan.store_collections:
        collection = by_digest_collection.get(digest)
        if collection is None:
            continue
        contributors.update(collection["from"])
    # Expand through the lock so every co-contributor of a shared collection is
    # re-ingested, not just the one whose damage was noticed.
    plan.store_files = sorted(build_lock.reingest_scope(lock, sorted(contributors)))

    for kind, owner in build_lock.owners_for_paths(lock, set(plan.store_files)):
        if kind != "seqset":
            plan.assembly_owners.add(owner)
    plan.cache_files = sorted(set(plan.cache_files))
    return plan


def _collections_publishing(
    lock: dict, store_dir: Path, digests: set[str]
) -> set[str]:
    """Which locked collections contain any of ``digests``.

    Deliberately opens the store read-only rather than trusting the lock: the
    lock records collection *sizes*, not membership, so this is the only way to
    get from a damaged sequence back to the sources that produced it.
    """
    from gtars.refget import RefgetStore

    import store_census

    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    found: set[str] = set()
    for collection in build_lock.lock_collections(lock):
        try:
            members = set(store_census.collection_members(
                store, collection["digest"]))
        except Exception:  # noqa: BLE001
            continue
        if members & digests:
            found.add(collection["digest"])
    return found


def print_plan(plan: RepairPlan) -> None:
    print(f"cache files to re-fetch:      {len(plan.cache_files)}")
    for rel in plan.cache_files:
        print(f"  {rel}")
    print(f"collections to re-ingest:     {len(plan.store_collections)}")
    for digest in plan.store_collections:
        print(f"  {digest}")
    print(f"source files to re-ingest:    {len(plan.store_files)}")
    for rel in plan.store_files:
        print(f"  {rel}")
    if plan.assembly_owners:
        print(f"\nassembly namespaces affected: {len(plan.assembly_owners)}")
        for owner in sorted(plan.assembly_owners):
            print(f"  {owner}  -> re-run: gks-refgetstore build --assembly {owner}")
        print("  repair re-ingests seqsets only; assembly ingestion carries its")
        print("  own deferred-report sequencing, which `build` owns.")
    if plan.unrepairable:
        print(f"\nunrepairable findings: {len(plan.unrepairable)}")
        seen: set[str] = set()
        for code, subject in plan.unrepairable:
            if code not in seen:
                seen.add(code)
                print(f"  {code}: {UNREPAIRABLE[code]}")
            print(f"    {subject}")
    if plan.ignored:
        print(f"\nnot addressed by repair: {len(plan.ignored)}")
        for code, subject in plan.ignored[:20]:
            print(f"  {code}  {subject}")
        if len(plan.ignored) > 20:
            print(f"  ... and {len(plan.ignored) - 20} more")


def repair_cache(plan: RepairPlan, lock: dict, args) -> bool:
    """Re-fetch every damaged cache file, accepting only the locked bytes."""
    by_rel = {r["cache_path"]: r for r in build_lock.lock_files(lock)}
    queued: list[tuple[ResolvedSource, Path]] = []
    expected: dict[str, str] = {}
    for rel in plan.cache_files:
        record = by_rel.get(rel)
        if record is None or not record.get("sha256"):
            print(f"  cannot repair {rel}: the lock pins no SHA-256 for it",
                  file=sys.stderr)
            return False
        target = args.cache_dir / rel
        target.unlink(missing_ok=True)
        expected[record["url"]] = record["sha256"]
        queued.append((locked_source(record), target))

    def accept(path: Path, source: ResolvedSource) -> bool:
        return sha256_file(path) == expected[source.url]

    outcome = fetch_sources.fetch_many(
        queued, args.cache_dir, timeout=args.timeout, retries=args.retries,
        jobs=args.jobs, min_free_gb=args.min_free_gb, accept=accept,
    )
    print(f"  re-fetched {outcome.fetched}, failed {outcome.failed}, "
          f"{fetch_sources.human(outcome.total_bytes)}")
    for url, reason in outcome.failures:
        print(f"  FAILED {url}\n    {reason}", file=sys.stderr)
    return outcome.ok


def repair_store(plan: RepairPlan, lock: dict, args) -> bool:
    """Remove the damaged collections and re-ingest their source files.

    Hands a synthesized :class:`build_lock.SyncPlan` to
    ``store_sync.apply_plan``, which already sequences prepare-filter ->
    load_for_mutation -> remove -> re-ingest -> write -> reconcile in the only
    order gtars makes safe.
    """
    import store_sync
    from gtars.refget import RefgetStore

    _, seqsets = build_store.load_config(args.config)
    by_rel = {r["cache_path"]: r for r in build_lock.lock_files(lock)}
    records = [by_rel[rel] for rel in plan.store_files
               if rel in by_rel and by_rel[rel]["kind"] == "seqset"]
    if not records:
        print("  nothing to re-ingest")
        return True

    missing = [r["cache_path"] for r in records
               if not (args.cache_dir / r["cache_path"]).exists()]
    if missing:
        print("  refusing to mutate the store: source file(s) absent from the "
              "cache", file=sys.stderr)
        for rel in missing:
            print(f"    {rel}", file=sys.stderr)
        return False

    synthesized = build_lock.SyncPlan(
        reingest=[(r["cache_path"], "repair: restoring the locked state")
                  for r in records],
        reingest_records=records,
        collection_by_file=build_lock.collection_by_file(lock),
    )
    store = RefgetStore.open_local(str(args.store_dir))
    store.set_quiet(True)
    report, provenance = store_sync.apply_plan(
        store, args.store_dir, synthesized, seqsets, args.cache_dir, lock,
        ingest_jobs=args.ingest_jobs, filter_jobs=args.filter_jobs,
    )
    print(f"  n_sequences: {report.sequences_before} -> "
          f"{report.sequences_after_removal} (after removal) -> "
          f"{report.sequences_final}")
    print(f"  collections removed: {report.collections_removed}; "
          f"created: {len(report.collections_added)}")
    return True


def run_repair(args) -> int:
    """Execute the ``repair`` subcommand."""
    targets = {t for t in ("cache", "store") if getattr(args, t, False)}
    if not targets:
        targets = {"cache", "store"}

    lock = build_lock.load_lock(args.lock)

    print("verifying before repairing ...", flush=True)
    report = verify.VerifyReport()
    if "cache" in targets:
        report.merge(verify.verify_cache(lock, args.cache_dir, jobs=args.jobs))
    if "store" in targets:
        # L0 + L1 only. L2 costs ~15 minutes and nothing it finds is
        # repairable -- an encoder round-trip defect is reproduced exactly by
        # re-ingesting. Run `verify --store --deep` for that separately.
        report.merge(verify.verify_store(lock, args.store_dir))
    if not report.findings:
        print("\nnothing to repair: cache and store match the lock")
        return 0

    plan = plan_repair(lock, report.findings, args.store_dir)
    print()
    print_plan(plan)

    # A root mismatch that no missing collection or payload explains is a
    # divergence, not damage. Refuse rather than guess.
    unexplained = [
        (code, subject) for code, subject in plan.unrepairable
        if code.endswith("_root_mismatch")
    ]
    if unexplained and not plan.store_collections:
        print("\nrefusing to repair the store: its digest set differs from the "
              "lock's with no attributable cause.", file=sys.stderr)
        return 1

    if not plan.touched:
        print("\nnothing repair can act on")
        return 1 if report.failed else 0
    if not args.apply:
        print("\ndry run; re-run with --apply to repair")
        return 0

    # Cache first, and abort the whole run if it fails: store removal persists
    # immediately and has no rollback, so a source file that turns out to be
    # unavailable must be discovered before anything is removed.
    if plan.cache_files:
        print("\nrepairing cache ...", flush=True)
        if not repair_cache(plan, lock, args):
            print("\ncache repair failed; refusing to touch the store",
                  file=sys.stderr)
            return 1
    if plan.store_files:
        print("\nrepairing store ...", flush=True)
        if not repair_store(plan, lock, args):
            return 1
    if plan.assembly_owners:
        print("\nassembly namespaces still need `build --assembly`; see above.")
        return 1
    print("\nrepair complete; re-run `gks-refgetstore verify` to confirm")
    return 0
