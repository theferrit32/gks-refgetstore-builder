#!/usr/bin/env python
"""Apply manifest changes to an existing store as incremental adds and deletes.

A full rebuild re-encodes every sequence to change one release. ``sync``
classifies each source against the build lock, removes only the collections
whose ingest inputs moved, re-ingests those, and rewrites only the affected
alias namespaces -- leaving the rest of the store byte-identical.

Three gtars 0.9.2 behaviours shape this module, all established empirically
against the installed package:

1. ``remove_collection(remove_orphan_sequences=True)`` collects **nothing** on a
   lazily loaded store, silently. Collections open as metadata stubs, and orphan
   detection consults only loaded ones. :func:`load_for_mutation` is therefore
   mandatory, and :func:`remove_collections` refuses to run without it.

2. Removal **persists immediately** -- orphaned ``.seq`` files are unlinked and
   the indexes rewritten without waiting for ``write()``. There is no in-memory
   staging and no abort path, so classification must be complete and checked
   before the first removal.

3. Alias namespaces cannot be shrunk through the alias API.
   ``load_sequence_aliases`` merges rather than replaces, and
   ``remove_sequence_alias`` rewrites the whole namespace TSV per call (~2s on a
   360k-alias namespace). :func:`reconcile_namespace` writes the namespace file
   directly instead, which is durable across ``write()``.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import build_store
from build_store import (INGEST_JOBS_DEFAULT, SeqsetConfig, _write_alias_tsv,
                         prepare_filtered_sources, process_seqset)

logger = logging.getLogger("store_sync")


class MutationOrderError(RuntimeError):
    """Raised when a removal is attempted before collections were loaded."""


@dataclass
class SyncReport:
    collections_removed: int = 0
    sequences_before: int = 0
    sequences_after_removal: int = 0
    sequences_final: int = 0
    namespaces_reconciled: dict[str, tuple[int, int]] = field(default_factory=dict)
    seqsets_processed: list[str] = field(default_factory=list)
    collections_added: dict[str, int] = field(default_factory=dict)


def store_sequence_count(store) -> int:
    return int(store.stats()["n_sequences"])


def load_for_mutation(store) -> int:
    """Load every collection so orphan detection can see all references.

    Returns the number loaded. Costs ~14s and ~1.7 GiB of additional RSS on the
    production store; it loads collection *metadata*, not sequence payloads.
    """
    store.load_all_collections()
    loaded = int(store.stats()["n_collections_loaded"])
    logger.info("loaded %d collection(s) for mutation", loaded)
    return loaded


def remove_collections(store, digests: list[str]) -> int:
    """Remove collections and collect orphaned sequences.

    Removals must all precede any re-ingest: a sequence published by several
    collections only becomes an orphan once the last one referencing it is
    gone. Interleaving remove-and-ingest per source reclaims almost nothing.
    """
    stats = store.stats()
    if int(stats["n_collections_loaded"]) < int(stats["n_collections"]):
        raise MutationOrderError(
            "collections are not fully loaded: orphan collection would "
            "silently free nothing. Call load_for_mutation() first "
            f"(loaded {stats['n_collections_loaded']} of {stats['n_collections']})"
        )
    removed = 0
    for digest in digests:
        if store.remove_collection(digest, remove_orphan_sequences=True):
            removed += 1
            logger.info("  removed collection %s (n_sequences=%s)",
                        digest, store.stats()["n_sequences"])
        else:
            logger.warning("  collection %s not present; nothing removed", digest)
    return removed


def reconcile_namespace(
    store_dir: Path, namespace: str, desired: dict[str, str]
) -> tuple[int, int]:
    """Make ``namespace`` hold exactly ``desired``; returns (added, removed).

    Writes the namespace TSV directly rather than going through the alias API.
    ``remove_sequence_alias`` is O(namespace) per call, so reconciling ~10k
    aliases across the ``ensembl-N`` namespaces would take hours; one atomic
    write is O(namespace) once. A hand-written namespace file survives a later
    ``store.write()``.

    The caller must not ``write()`` the store afterwards: an open store's
    in-memory alias map does not know about the file.
    """
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", namespace):
        raise ValueError(f"unsafe namespace: {namespace!r}")
    path = store_dir / "aliases" / "sequences" / f"{namespace}.tsv"
    existing: dict[str, str] = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            alias, _, digest = line.partition("\t")
            existing[alias] = digest
    added = len(set(desired) - set(existing))
    removed = len(set(existing) - set(desired))
    _write_alias_tsv(path, desired)
    logger.info("  namespace %s: %d alias(es) (+%d, -%d)",
                namespace, len(desired), added, removed)
    return added, removed


def seqsets_for_namespaces(
    seqsets: list[SeqsetConfig], namespaces: set[str]
) -> list[SeqsetConfig]:
    """Every seqset feeding any of ``namespaces``.

    A release namespace is fed by several seqsets -- ``ensembl-100`` draws from
    ``dna.toplevel``, ``cdna``, ``ncrna`` and ``pep``. Reconciling it requires
    the complete desired alias map, so all contributors are reprocessed even
    when only one of them changed. Re-ingesting an unchanged collection is ~8x
    cheaper than a fresh one, which is what makes this affordable.
    """
    return [entry for entry in seqsets if entry.namespace in namespaces]


def apply_plan(
    store,
    store_dir: Path,
    plan,
    seqsets: list[SeqsetConfig],
    download_dir: Path,
    lock: dict,
    *,
    ingest_jobs: int = INGEST_JOBS_DEFAULT,
    filter_jobs: int = INGEST_JOBS_DEFAULT,
    force_download: bool = False,
) -> tuple[SyncReport, dict[str, dict]]:
    """Execute a :class:`build_lock.SyncPlan` against an open store.

    Returns the report and the build provenance (cache path -> collection),
    which the caller merges into the lock.

    The re-ingest scope is widened twice, for two independent reasons:

    1. **Collection coverage.** Removal is per collection, so dropping one
       because a single contributor changed also drops every *other*
       contributor's sequences. ``build_lock.reingest_scope`` follows the
       lock's ``from`` lists to recover them.
    2. **Alias coverage.** A rolling release namespace is fed by several
       seqsets, and reconciling it needs the complete desired alias map, so
       every contributor to a touched namespace is reprocessed.
    """
    import build_lock

    report = SyncReport()
    by_name = {entry.name: entry for entry in seqsets}
    touched_paths = (
        {record["cache_path"] for record in plan.removed}
        | set(plan.added)
        | {rel for rel, _ in plan.reingest}
    )
    # Expanding by collection is correct by construction; expanding by
    # namespace was correct by coincidence. In the committed lock every
    # collection with contributors in more than one seqset is an Ensembl
    # release group, which the namespace expansion below already covers -- but
    # content-addressing does not respect namespace boundaries, so that holds
    # only until the manifest grows a shared collection that crosses them.
    covered = build_lock.reingest_scope(lock, sorted(touched_paths))
    touched_owners = {
        owner for _kind, owner in build_lock.owners_for_paths(lock, covered)
    }
    touched_owners |= _owners_of(plan, seqsets, download_dir)
    affected = [by_name[n] for n in sorted(touched_owners) if n in by_name]
    namespaces = {entry.namespace for entry in affected}
    scope = seqsets_for_namespaces(seqsets, namespaces)
    logger.info(
        "touched files: %d; expanded to %d by collection membership; "
        "affected seqsets: %d; expanded to %d for namespace coverage",
        len(touched_paths), len(covered), len(affected), len(scope),
    )

    report.sequences_before = store_sequence_count(store)

    # Phase 1: materialize derived inputs before touching the store, so a
    # filter failure aborts before anything is removed.
    prepare_filtered_sources(scope, download_dir, filter_jobs, force_download)

    # Phase 2: every removal, before any ingest.
    load_for_mutation(store)
    report.collections_removed = remove_collections(
        store, plan.collection_digests_to_remove
    )
    report.sequences_after_removal = store_sequence_count(store)
    logger.info("after removal: n_sequences=%d (was %d)",
                report.sequences_after_removal, report.sequences_before)

    # Phase 3: re-ingest, collecting the complete desired alias map per
    # namespace. alias_sink makes process_seqset return the map instead of
    # writing it, so nothing is persisted add-only.
    provenance: dict[str, dict] = {}
    alias_maps: dict[str, dict[str, str]] = {}
    for entry in scope:
        sink = alias_maps.setdefault(entry.namespace, {})
        stats = process_seqset(
            store, entry, download_dir, provenance=provenance,
            alias_sink=sink, jobs=ingest_jobs,
        )
        report.seqsets_processed.append(entry.name)
        if stats.warnings:
            raise RuntimeError(
                f"seqset {entry.name!r} reported {stats.warnings} warning(s); "
                "refusing to continue a partially applied sync"
            )

    # Phase 4: persist sequences and indexes, THEN rewrite alias namespaces.
    store.write()
    report.sequences_final = store_sequence_count(store)
    # Content-addressing means the re-ingested collections need not map 1:1 to
    # the removed ones. Filtering the Ensembl toplevels collapses 8 collections
    # into 3, because the padding was most of what distinguished the release
    # groups. Report what was actually created rather than assuming symmetry.
    report.collections_added = {
        rec["collection_digest"]: rec["n_sequences"]
        for rec in provenance.values() if rec.get("collection_digest")
    }
    for namespace, desired in sorted(alias_maps.items()):
        report.namespaces_reconciled[namespace] = reconcile_namespace(
            store_dir, namespace, desired
        )
    return report, provenance


def print_plan(plan, limit: int = 40) -> None:
    """Human-readable classification, reason-first."""
    print(f"unchanged={len(plan.unchanged)}  added={len(plan.added)}  "
          f"removed={len(plan.removed)}  reingest={len(plan.reingest)}  "
          f"missing={len(plan.missing)}")
    for label, items in (("added", plan.added),
                         ("missing", plan.missing),
                         ("removed", [r["cache_path"] for r in plan.removed])):
        if items:
            print(f"\n{label}:")
            for rel in items[:limit]:
                print(f"  {rel}")
            if len(items) > limit:
                print(f"  ... and {len(items) - limit} more")
    if plan.reingest:
        print("\nreingest:")
        for rel, reason in plan.reingest[:limit]:
            print(f"  {rel}\n      reason: {reason}")
        if len(plan.reingest) > limit:
            print(f"  ... and {len(plan.reingest) - limit} more")
    digests = plan.collection_digests_to_remove
    print(f"\ncollections to remove: {len(digests)}")
    for digest in digests:
        print(f"  {digest}")


def run_sync(args) -> int:
    """Execute the ``sync`` subcommand."""
    from gtars.refget import RefgetStore

    import build_lock

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    assemblies, seqsets = build_store.load_config(args.config)
    lock = build_lock.load_lock(args.lock)
    sources = build_store.resolve_sources(assemblies, seqsets)
    plan = build_lock.plan_sync(
        sources, lock, args.cache_dir, seqsets=seqsets,
        hash_files=not args.no_hash,
    )
    print_plan(plan)

    if plan.missing:
        print("\nrefusing to sync: sources missing from the cache", flush=True)
        return 1
    if not plan.touched:
        print("\nstore already matches the manifest; nothing to do", flush=True)
        return 0
    if not args.apply:
        print("\ndry run; re-run with --apply to mutate the store", flush=True)
        return 0

    print(f"\napplying to {args.store_dir} ...", flush=True)
    store = RefgetStore.open_local(str(args.store_dir))
    report, provenance = apply_plan(
        store, args.store_dir, plan, seqsets, args.cache_dir, lock,
        ingest_jobs=args.ingest_jobs, filter_jobs=args.filter_jobs,
    )
    print(f"\nn_sequences: {report.sequences_before} -> "
          f"{report.sequences_after_removal} (after removal) -> "
          f"{report.sequences_final}")
    print(f"net change: {report.sequences_final - report.sequences_before:+d}")
    print(f"collections removed: {report.collections_removed}; "
          f"created: {len(report.collections_added)}; "
          f"seqsets reprocessed: {len(report.seqsets_processed)}")
    for digest, n in sorted(report.collections_added.items(),
                            key=lambda kv: -kv[1]):
        print(f"  + {digest}  n_sequences={n}")
    for namespace, (added, removed) in sorted(report.namespaces_reconciled.items()):
        if added or removed:
            print(f"  {namespace}: +{added} -{removed}")

    touched = [
        source for source in sources
        if source.owner in set(report.seqsets_processed)
    ]
    merged = build_lock.merge_into_lock(
        lock, config_path=args.config, download_dir=args.cache_dir,
        touched_sources=touched, collection_by_cachepath=provenance,
        store=store, seqsets=seqsets,
        dropped_scopes={(r["kind"], r["owner"]) for r in plan.removed},
    )
    build_lock.validate_lock(merged)
    if args.no_lock:
        print("\nskipping lock write (--no-lock)")
    else:
        build_lock.write_lock(args.lock, merged)
        print(f"\nwrote {args.lock} "
              f"({len(build_lock.lock_files(merged))} input file(s), "
              f"{merged['outputs']['n_collections']} collection(s))")
    return 0


def _owners_of(plan, seqsets, download_dir: Path) -> set[str]:
    """Seqset names owning any added/reingest source, keyed by cache path."""
    rel_to_owner: dict[str, str] = {}
    for entry in seqsets:
        for _, url in entry.iter_shard_urls():
            path = build_store.mirror_cache_path(download_dir, url)
            try:
                rel = str(path.relative_to(download_dir))
            except ValueError:
                rel = str(path)
            rel_to_owner[rel] = entry.name
    touched = set(plan.added) | {rel for rel, _ in plan.reingest}
    return {rel_to_owner[rel] for rel in touched if rel in rel_to_owner}
