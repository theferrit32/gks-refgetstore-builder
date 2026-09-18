#!/usr/bin/env python
"""Scan a store for sequences whose stored payload does not re-digest to its key.

Generates the `--deep` known-bad baseline
(`seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv`). It uses
``store_census.redigest_sequence`` -- the same primitive ``verify --deep`` uses
-- so the baseline cannot disagree with the check it feeds.

    uv run python runs/2026-09-12-deep-verify/scan_roundtrip.py \
        --store store --out roundtrip.tsv [--limit N]

Single-threaded on purpose: whether gtars releases the GIL inside
``stream_sequence`` has not been measured, and a thread pool that does not help
would make the timing numbers in the README meaningless.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from gks_refgetstore import store_census

# Diagnosed by comparing the decoded payload against the source FASTA record for
# each group, not inferred. Keyed by alphabet because that is what separates
# them mechanically; the diagnosis itself is not mechanical, which is why verify
# cannot derive it and this column has to be carried forward by hand.
#
#   protein: ENSP00000473614 (GPX4, release-76 pep, 180 aa). Source has U at
#            index 109; the store decodes A there. Source charset
#            ACDEFGHIKLMNPQRSTUVWY, decoded charset ACDEFGHIKLMNPQRSTVWY.
#   dnaio:   ENSP00000499040.1 (NOTCH2, release-113 pep, 12 aa). Source
#            MCVTYHNGTGYC, decoded MCVTYVNGTGYC -- H -> V at index 5. Every
#            residue is a legal IUPAC nucleotide code, so the alphabet is
#            misdetected and the IUPAC round trip is itself lossy.
CAUSES = {
    "protein": (
        "gtars protein alphabet has no U (selenocysteine); residue replaced on "
        "encode. Upstream defect, not bit rot."
    ),
    "dnaio": (
        "short protein whose residues are all legal IUPAC nucleotide codes; "
        "alphabet misdetected as nucleotide and encoded lossily. Upstream "
        "defect, not bit rot."
    ),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--store", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(str(args.store))
    store.set_quiet(True)
    metadata = store_census.sequence_metadata_by_digest(store)
    digests = sorted(metadata)
    if args.limit:
        digests = digests[: args.limit]

    started = time.monotonic()
    rows = []
    by_alphabet: dict[str, int] = {}
    total_by_alphabet: dict[str, int] = {}
    for i, digest in enumerate(digests, 1):
        meta = metadata[digest]
        # AlphabetType is a gtars enum, not a str; everything downstream (TSV
        # cells, dict keys, the cause lookup) wants its name.
        alphabet = str(meta.alphabet)
        total_by_alphabet[alphabet] = total_by_alphabet.get(alphabet, 0) + 1
        if i % 100_000 == 0:
            rate = i / max(1e-9, time.monotonic() - started)
            print(f"  {i:,}/{len(digests):,} ({rate:,.0f}/s) "
                  f"{len(rows)} mismatch(es)", file=sys.stderr, flush=True)
        actual = store_census.redigest_sequence(store, digest)
        if actual != digest:
            by_alphabet[alphabet] = by_alphabet.get(alphabet, 0) + 1
            rows.append((digest, meta.name, alphabet, str(meta.length),
                         actual, CAUSES.get(alphabet, "undiagnosed")))

    took = time.monotonic() - started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        handle.write("digest\tname\talphabet\tlength\tredigest\tcause\n")
        for row in sorted(rows, key=lambda r: (r[2], r[1])):
            handle.write("\t".join(row) + "\n")

    print(f"\nscanned {len(digests):,} in {took:.0f}s "
          f"({len(digests)/max(1e-9,took):,.0f}/s)")
    print("per-alphabet totals and round-trip failures:")
    for alphabet in sorted(total_by_alphabet):
        print(f"  {alphabet:10} {total_by_alphabet[alphabet]:>9,}  "
              f"failures={by_alphabet.get(alphabet, 0)}")
    print(f"wrote {args.out} ({len(rows)} row(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
