#!/usr/bin/env python
"""Answer file <-> refget-digest provenance from a build.lock.json + the store.

The lock records, per collection, the source files that contributed to it
(``outputs.collections[].from``); the store records each collection's sequence
membership. Composing the two answers both directions on demand, so no giant
provenance table has to be stored.

A collection can have **many** contributing files -- byte-identical sources
produce one collection, and filtering collapsed eight Ensembl collections into
three. 43 collections in the committed lock have two or more contributors and
one has twelve. The inverse holds: one file yields exactly one collection.

Modes:
    --file SUBSTR    list the sequence digests contributed by matching source(s)
    --digest DIGEST  list the source file(s) whose collection contains that digest
    --list           list every collection, its size, and its contributors

Examples:
    uv run python provenance.py --list
    uv run python provenance.py --file human.6.rna
    uv run python provenance.py --digest SQ.abc123…   # SQ. prefix optional
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_lock  # noqa: E402
from store_census import collection_members, strip_sq  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent


def _open_store(store_dir: Path):
    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    return store


def _matching_collections(lock: dict, needle: str) -> list[dict]:
    """Collections with a contributor whose path, URL, or owner matches."""
    by_rel = {record["cache_path"]: record for record in build_lock.lock_files(lock)}
    out = []
    for collection in build_lock.lock_collections(lock):
        hits = [
            rel for rel in collection["from"]
            if needle in rel
            or needle in by_rel.get(rel, {}).get("url", "")
            or needle in by_rel.get(rel, {}).get("owner", "")
        ]
        if hits:
            out.append(collection)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--lock", type=Path, default=REPO_ROOT / "build.lock.json")
    ap.add_argument("--store-dir", type=Path, default=REPO_ROOT / "store")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--file", type=str, help="source URL/cache-path/owner substring")
    g.add_argument("--digest", type=str, help="a sequence sha512t24u (SQ. optional)")
    g.add_argument("--list", action="store_true",
                   help="list all collections + contributors")
    ap.add_argument("--limit", type=int, default=None,
                    help="cap digests printed for --file")
    args = ap.parse_args()

    lock = build_lock.load_lock(args.lock)
    collections = build_lock.lock_collections(lock)

    if args.list:
        for collection in collections:
            sources = collection["from"] or ["(no known contributor)"]
            print(f"{collection['digest']}\t{collection['n_sequences']:>8}\t"
                  f"{len(sources)} file(s)")
            for rel in sources:
                print(f"\t\t{rel}")
        unattributed = sum(1 for c in collections if not c["from"])
        print(f"\n{len(collections)} collection(s), {unattributed} with no known "
              f"contributor", file=sys.stderr)
        return 0

    store = _open_store(args.store_dir)

    if args.file:
        matches = _matching_collections(lock, args.file)
        if not matches:
            print(f"no source matches {args.file!r}", file=sys.stderr)
            return 1
        for collection in matches:
            digests = collection_members(store, collection["digest"])
            print(f"# collection={collection['digest']}  "
                  f"({len(digests)} sequences)  from:")
            for rel in collection["from"]:
                print(f"#   {rel}")
            shown = digests[: args.limit] if args.limit else digests
            for digest in shown:
                print(digest)
            if args.limit and len(digests) > args.limit:
                print(f"# … +{len(digests) - args.limit} more", file=sys.stderr)
        return 0

    # --digest: which collections contain it. An O(all collections) scan, ~56s
    # on the production store. Deliberately not backed by a cached reverse
    # index: that is 1.8M entries to build, persist, and keep honest, to save
    # under a minute on an interactive query.
    want = strip_sq(args.digest)
    found = [
        collection for collection in collections
        if want in set(collection_members(store, collection["digest"]))
    ]
    if not found:
        print(f"digest {want} not found in any collection", file=sys.stderr)
        return 1
    by_rel = {record["cache_path"]: record for record in build_lock.lock_files(lock)}
    print(f"# digest {want} appears in {len(found)} collection(s):")
    for collection in found:
        print(f"# collection={collection['digest']}")
        for rel in collection["from"]:
            record = by_rel.get(rel, {})
            print(f"{record.get('owner', '?')}\t{record.get('url', rel)}")
        if not collection["from"]:
            print("?\t(no known contributor)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
