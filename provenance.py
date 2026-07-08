#!/usr/bin/env python
"""Answer file <-> refget-digest provenance from a build.lock.json + the store.

Each ingested source file became exactly one collection, and the store persists
every collection's sequence membership. So the lock's ``url -> collection_digest``
plus the store is all we need — no giant provenance table is stored; the mapping
is derived on demand here.

Modes:
    --file SUBSTR    list the sequence digests contributed by matching source(s)
    --digest DIGEST  list the source file(s) whose collection contains that digest
    --list           list every source and its collection digest / count

Examples:
    uv run python provenance.py --list
    uv run python provenance.py --file human.6.rna
    uv run python provenance.py --digest SQ.abc123…   # SQ. prefix optional
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gtars.refget import RefgetStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_lock  # noqa: E402
from build_store import SQ_PREFIX  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent


def collection_members(store: RefgetStore, coll_digest: str) -> list[str]:
    store.load_collection(coll_digest)
    level2 = store.get_collection_level2(coll_digest)
    out = []
    for seq_id in level2["sequences"]:
        out.append(seq_id[len(SQ_PREFIX):] if seq_id.startswith(SQ_PREFIX) else seq_id)
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
    g.add_argument("--list", action="store_true", help="list all sources + collections")
    ap.add_argument("--limit", type=int, default=None, help="cap digests printed for --file")
    args = ap.parse_args()

    lock = build_lock.load_lock(args.lock)
    fasta = [s for s in lock["sources"] if s.get("collection_digest")]

    if args.list:
        for s in fasta:
            print(f"{s['collection_digest']}\t{s['n_sequences']:>8}\t"
                  f"{s['owner']}\t{s['cache_path']}")
        print(f"\n{len(fasta)} collection-bearing source(s)", file=sys.stderr)
        return 0

    store = RefgetStore.open_local(str(args.store_dir))

    if args.file:
        matches = [s for s in fasta
                   if args.file in s["url"] or args.file in s["cache_path"]
                   or args.file in s["owner"]]
        if not matches:
            print(f"no source matches {args.file!r}", file=sys.stderr)
            return 1
        for s in matches:
            digests = collection_members(store, s["collection_digest"])
            print(f"# {s['owner']}  {s['cache_path']}  "
                  f"collection={s['collection_digest']}  ({len(digests)} sequences)")
            shown = digests[: args.limit] if args.limit else digests
            for d in shown:
                print(d)
            if args.limit and len(digests) > args.limit:
                print(f"# … +{len(digests) - args.limit} more", file=sys.stderr)
        return 0

    # --digest: which source collections contain it (invert once)
    want = args.digest[len(SQ_PREFIX):] if args.digest.startswith(SQ_PREFIX) else args.digest
    found = []
    for s in fasta:
        if want in set(collection_members(store, s["collection_digest"])):
            found.append(s)
    if not found:
        print(f"digest {want} not found in any source collection", file=sys.stderr)
        return 1
    print(f"# digest {want} appears in {len(found)} source file(s):")
    for s in found:
        print(f"{s['owner']}\t{s['url']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
