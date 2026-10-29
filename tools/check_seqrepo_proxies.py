#!/usr/bin/env python3
"""Resolve every RefgetStore alias through both vrs-python data proxies.

For each sequence alias in the store, asks vrs-python's SeqRepo data proxy
(``seqrepo+file://``) and its RefgetStore data proxy (``refget+file://``) for
the alias's refget accession, and requires both to return ``SQ.<digest>`` for
the digest the store files the alias under. Accession namespaces that
vrs-python infers from a bare accession (``refseq``, ``lrg``, ``insdc``) are
also resolved bare, the way HGVS and SPDI expressions name them. Sequence
lengths are compared for every sequence, and residues for a deterministic
sample: both ends of each sampled sequence.

Needs vrs-python with its refget extra plus biocommons.seqrepo, so run it from
that environment, e.g.:

    uv run --project <vrs-python checkout> python tools/check_seqrepo_proxies.py \\
        --seqrepo <root>/<instance> --store <store dir> --out <tsv>

Exits non-zero if any identifier resolves differently.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

BARE_NAMESPACES = ("refseq", "lrg", "insdc")
SLICE = 25


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seqrepo", type=Path, required=True, help="SeqRepo instance directory")
    ap.add_argument("--store", type=Path, required=True, help="RefgetStore directory")
    ap.add_argument("--out", type=Path, required=True,
                    help="TSV of every identifier that resolved differently")
    ap.add_argument("--slice-every", type=int, default=25,
                    help="compare residues for every Nth sequence (default 25)")
    ap.add_argument("--namespace", action="append", default=None,
                    help="only resolve aliases in this store namespace (repeatable); "
                         "lengths and residues are still compared for every sequence")
    args = ap.parse_args()

    from ga4gh.vrs.dataproxy import coerce_accession_namespace, create_dataproxy
    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(str(args.store))
    store.set_quiet(True)
    seqrepo = create_dataproxy(f"seqrepo+file://{args.seqrepo.resolve()}")
    refget = create_dataproxy(f"refget+file://{args.store.resolve()}")

    started = time.monotonic()
    rows = []
    n = bare = 0
    for ns in sorted(store.list_sequence_alias_namespaces()):
        if args.namespace and ns not in args.namespace:
            continue
        for alias in sorted(store.list_sequence_aliases(ns) or []):
            digest = store.get_sequence_metadata_by_alias(ns, alias).sha512t24u
            identifiers = [f"{ns}:{alias}"]
            if ns in BARE_NAMESPACES:
                try:
                    if coerce_accession_namespace(alias) == f"{ns}:{alias}":
                        identifiers.append(alias)
                        bare += 1
                except ValueError:
                    pass
            for ident in identifiers:
                n += 1
                expected = f"SQ.{digest}"
                got_sr = seqrepo.derive_refget_accession(ident)
                got_rg = refget.derive_refget_accession(ident)
                if not (got_sr == got_rg == expected):
                    rows.append(("accession", ident, expected, got_sr, got_rg))
    alias_seconds = time.monotonic() - started

    metas = sorted(store.list_sequences(), key=lambda m: m.sha512t24u)
    sliced = 0
    for i, meta in enumerate(metas):
        ident = f"ga4gh:SQ.{meta.sha512t24u}"
        a, b = seqrepo.get_sequence_length(ident), refget.get_sequence_length(ident)
        if not (a == b == meta.length):
            rows.append(("length", ident, meta.length, a, b))
        if i % args.slice_every == 0:
            sliced += 1
            for start, end in ((0, SLICE), (max(meta.length - SLICE, 0), meta.length)):
                a = seqrepo.get_sequence(ident, start, end)
                b = refget.get_sequence(ident, start, end)
                if a != b or len(a) != end - start:
                    rows.append(("residues", f"{ident}[{start}:{end}]", "", a, b))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(("check", "identifier", "expected", "seqrepo", "refget"))
        writer.writerows(rows)
    print(f"identifiers resolved: {n} ({bare} bare accessions) in {alias_seconds:.0f}s")
    print(f"sequences: {len(metas)} lengths compared, {sliced} sliced at both ends")
    print(f"differences: {len(rows)} -> {args.out}")
    print("PASS" if not rows else "FAIL")
    return 1 if rows else 0


if __name__ == "__main__":
    sys.exit(main())
