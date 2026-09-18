#!/usr/bin/env python
"""Compare a mutated store against an unmutated baseline copy.

Written to check a ``gks-refgetstore sync`` run, but the comparison is generic:
point it at any two stores and it reports what moved.

    uv run python tools/compare_store_mutation.py --baseline store \
        --mutated store.pre-filter

With ``--expect-removed-digests <tsv> --digest-column padded_digest`` it also
asserts that every listed digest is gone and, via ``--survivor-column``, that a
paired digest survived -- the "did we delete exactly the right things" check.

Note on the gtars API: ``get_sequence_metadata`` and
``get_sequence_metadata_by_alias`` return ``None`` for a missing entry rather
than raising, so presence must be tested with ``is not None``. Wrapping them in
``try/except`` silently reports everything as present.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from gks_refgetstore.store_census import collection_census, seq_path


class Checks:
    def __init__(self) -> None:
        self.failed = 0
        self.passed = 0

    def __call__(self, label, actual, expected) -> bool:
        good = actual == expected
        self.passed += good
        self.failed += not good
        print(f"[{'PASS' if good else 'FAIL'}] {label}: {actual!r}"
              + ("" if good else f"  (expected {expected!r})"))
        return good

    def note(self, label, value) -> None:
        print(f"[ -- ] {label}: {value}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--mutated", type=Path, required=True)
    ap.add_argument("--expect-removed-digests", type=Path, default=None,
                    help="TSV listing digests that must no longer resolve")
    ap.add_argument("--digest-column", default="padded_digest")
    ap.add_argument("--survivor-column", default=None,
                    help="column naming a digest that must still resolve")
    ap.add_argument("--expect-sequence-delta", type=int, default=None,
                    help="required n_sequences change (e.g. -445)")
    ap.add_argument("--namespace-sample", nargs="*", default=[],
                    help="namespaces to report alias counts for")
    ap.add_argument("--forbid-alias-prefix", default=None,
                    help="fail if any namespace still holds an alias with "
                         "this prefix (e.g. CHR_)")
    args = ap.parse_args(argv)

    from gtars.refget import RefgetStore

    base = RefgetStore.open_local(str(args.baseline)); base.set_quiet(True)
    mut = RefgetStore.open_local(str(args.mutated)); mut.set_quiet(True)
    check = Checks()
    b, m = base.stats(), mut.stats()

    print("=== counts ===")
    check.note("baseline n_sequences", b["n_sequences"])
    check.note("mutated n_sequences", m["n_sequences"])
    delta = int(m["n_sequences"]) - int(b["n_sequences"])
    if args.expect_sequence_delta is not None:
        check("n_sequences delta", delta, args.expect_sequence_delta)
    else:
        check.note("n_sequences delta", f"{delta:+d}")

    print("\n=== collections ===")
    bc, mc = collection_census(base), collection_census(mut)
    added = {d: n for d, n in mc.items() if d not in bc}
    gone = {d: n for d, n in bc.items() if d not in mc}
    check.note("baseline collections", len(bc))
    check.note("mutated collections", len(mc))
    print(f"  removed ({len(gone)}):")
    for d, n in sorted(gone.items(), key=lambda kv: -kv[1]):
        print(f"    - {d}  n_sequences={n}")
    print(f"  created ({len(added)}):")
    for d, n in sorted(added.items(), key=lambda kv: -kv[1]):
        print(f"    + {d}  n_sequences={n}")

    if args.expect_removed_digests:
        print("\n=== expected removals ===")
        rows = list(csv.DictReader(
            args.expect_removed_digests.open(), delimiter="\t"))
        targets = [r[args.digest_column] for r in rows if r.get(args.digest_column)]
        check.note("digests listed", len(targets))
        on_disk = [d for d in targets if seq_path(args.mutated, d).exists()]
        check("target .seq files remaining", len(on_disk), 0)
        # is not None: gtars returns None for a missing digest, never raises.
        resolving = [d for d in targets
                     if mut.get_sequence_metadata(d) is not None]
        check("target digests still resolving", len(resolving), 0)
        if args.survivor_column:
            survivors = {r[args.survivor_column] for r in rows
                         if r.get(args.survivor_column)}
            missing = [d for d in survivors
                       if mut.get_sequence_metadata(d) is None]
            check.note("paired survivors listed", len(survivors))
            check("survivors wrongly removed", len(missing), 0)

    print("\n=== alias namespaces ===")
    ns_b = set(base.list_sequence_alias_namespaces())
    ns_m = set(mut.list_sequence_alias_namespaces())
    check("namespace set unchanged", ns_m == ns_b, True)
    if ns_m != ns_b:
        print("   only in mutated:", sorted(ns_m - ns_b))
        print("   only in baseline:", sorted(ns_b - ns_m))
    for ns in args.namespace_sample:
        mb = len(base.list_sequence_aliases(ns) or [])
        mm = len(mut.list_sequence_aliases(ns) or [])
        check.note(f"{ns} aliases", f"{mb} -> {mm} ({mm - mb:+d})")

    if args.forbid_alias_prefix:
        print(f"\n=== no {args.forbid_alias_prefix} aliases remain ===")
        leftover = {}
        for ns in sorted(ns_m):
            n = sum(1 for a in (mut.list_sequence_aliases(ns) or [])
                    if a.startswith(args.forbid_alias_prefix))
            if n:
                leftover[ns] = n
        check(f"namespaces holding {args.forbid_alias_prefix} aliases",
              leftover, {})

    print(f"\n{check.passed} passed, {check.failed} failed")
    return 1 if check.failed else 0


if __name__ == "__main__":
    sys.exit(main())
