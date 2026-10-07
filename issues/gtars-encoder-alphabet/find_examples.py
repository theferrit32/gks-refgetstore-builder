#!/usr/bin/env python
"""List stored sequences whose returned bytes do not match their digest.

Finds example accessions in an existing store for the cases described in
README.md in the same directory. For a self-contained reproduction that needs
no store, use repro_minimal.py instead.

Requires the gks-refgetstore-builder repository
(https://github.com/theferrit32/gks-refgetstore-builder) and a built store:

    uv run python issues/gtars-encoder-alphabet/find_examples.py \
        --store store --out issues/gtars-encoder-alphabet/examples.tsv

For each sequence in the selected alphabets, it re-digests the bytes returned
by ``stream_sequence`` and compares the result with the stored digest. It stops
once it has ``--limit`` examples per alphabet, so it finishes in seconds rather
than scanning all 1.78M sequences. The accession prefix (``ENSP``, ``NP_``,
``NR_``) is reported in its own column.

The complete list of affected sequences in the gks-refgetstore-builder store is
``seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv``, produced
by a full scan (``runs/2026-09-12-deep-verify/scan_roundtrip.py``, following
``RUNBOOK.md``). This script's output is a sample and does not replace it.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

from gks_refgetstore import store_census

# These are the only two alphabets in which returned bytes have differed from
# the input, and re-digesting a payload is the expensive part of the scan.
# Restricting to them by default is what makes this finish in seconds instead
# of ten minutes; pass --alphabet explicitly to widen or narrow it.
DEFAULT_ALPHABETS = ("protein", "dnaio")

# ENSP00000473614 -> ENSP, NR_103745.1 -> NR_, NP_001034.1 -> NP_. The prefix is
# the source's accession class, which is the axis the README's breakdown uses.
ACCESSION_PREFIX = re.compile(r"^([A-Za-z]+_?)")

# `ensembl-114` and `ensembl-76` are the same namespace at different releases;
# folding the suffix keeps the alias column from restating one accession once
# per release that carried it.
RELEASE_SUFFIX = re.compile(r"-\d+$")

ALIAS_CAP = 6

COLUMNS = ("alphabet", "prefix", "name", "length", "digest", "redigest", "aliases")


def accession_prefix(name: str) -> str:
    """Leading accession class of a record name, or ``?`` if it has none."""
    match = ACCESSION_PREFIX.match(name or "")
    return match.group(1) if match else "?"


def aliases_for(store, digest: str, cap: int = ALIAS_CAP) -> str:
    """``namespace:accession`` aliases for one sequence, condensed.

    A sequence carried by 40 Ensembl releases has 40 aliases that differ only
    in the release suffix, which would make this column hundreds of characters
    wide and the table unreadable. Release suffixes are folded away and the
    result is capped, because the column exists to identify the sequence, not
    to enumerate every place it appears.

    A sequence with no alias is normal -- it means every collection holding it
    named it without a recognized namespace -- so this degrades to an empty
    cell rather than failing the scan.
    """
    try:
        pairs = store.get_aliases_for_sequence(digest)
    except Exception:  # noqa: BLE001 - a sequence without aliases is not an error here
        return ""
    folded = sorted({
        f"{RELEASE_SUFFIX.sub('', namespace)}:{accession}" for namespace, accession in pairs
    })
    if len(folded) <= cap:
        return ",".join(folded)
    return ",".join(folded[:cap]) + f",+{len(folded) - cap} more"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=None,
                        help="TSV to write; without it, results only go to stdout")
    parser.add_argument("--alphabet", action="append", default=None,
                        help=f"alphabet to scan, repeatable (default: {', '.join(DEFAULT_ALPHABETS)})")
    parser.add_argument("--limit", type=int, default=5,
                        help="stop after this many examples per alphabet (0 = no limit)")
    args = parser.parse_args()

    from gtars.refget import RefgetStore

    alphabets = tuple(args.alphabet) if args.alphabet else DEFAULT_ALPHABETS
    limit = args.limit or None

    store = RefgetStore.open_local(str(args.store))
    store.set_quiet(True)
    metadata = store_census.sequence_metadata_by_digest(store)

    # Sorted so a rerun against an unchanged store returns the same examples;
    # an arbitrary five would make the README's accessions unreproducible.
    candidates = sorted(
        digest for digest, meta in metadata.items() if str(meta.alphabet) in alphabets
    )
    print(
        f"{len(metadata):,} sequences in store, {len(candidates):,} in "
        f"{'/'.join(alphabets)}; re-digesting until "
        f"{limit if limit else 'all'} example(s) per alphabet",
        file=sys.stderr,
    )

    started = time.monotonic()
    rows = []
    found: dict[str, int] = dict.fromkeys(alphabets, 0)
    examined = 0
    for digest in candidates:
        if limit and all(found.get(a, 0) >= limit for a in alphabets):
            break
        meta = metadata[digest]
        alphabet = str(meta.alphabet)
        if limit and found.get(alphabet, 0) >= limit:
            continue
        examined += 1
        redigest = store_census.redigest_sequence(store, digest)
        if redigest == digest:
            continue
        found[alphabet] = found.get(alphabet, 0) + 1
        rows.append((alphabet, accession_prefix(meta.name), meta.name, str(meta.length),
                     digest, redigest, aliases_for(store, digest)))

    took = time.monotonic() - started
    widths = [max(len(COLUMNS[i]), *(len(row[i]) for row in rows)) if rows else len(COLUMNS[i])
              for i in range(len(COLUMNS))]
    print("  ".join(name.ljust(width) for name, width in zip(COLUMNS, widths, strict=True)))
    for row in sorted(rows):
        print("  ".join(cell.ljust(width) for cell, width in zip(row, widths, strict=True)))

    by_prefix: dict[tuple[str, str], int] = {}
    for alphabet, prefix, *_ in rows:
        by_prefix[(alphabet, prefix)] = by_prefix.get((alphabet, prefix), 0) + 1
    print(f"\n{len(rows)} example(s) from {examined:,} re-digests in {took:.1f}s",
          file=sys.stderr)
    for (alphabet, prefix), count in sorted(by_prefix.items()):
        print(f"  {alphabet:8} {prefix:6} {count}", file=sys.stderr)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", encoding="utf-8") as handle:
            handle.write(
                "# Example sequences whose returned bytes do not match their digest;\n"
                "# see README.md. A sample produced by find_examples.py. The complete\n"
                "# list is seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv\n"
                "# in https://github.com/theferrit32/gks-refgetstore-builder.\n"
                f"# --store {args.store} --alphabet {'/'.join(alphabets)} "
                f"--limit {args.limit}\n"
            )
            handle.write("\t".join(COLUMNS) + "\n")
            for row in sorted(rows):
                handle.write("\t".join(row) + "\n")
        print(f"wrote {args.out} ({len(rows)} row(s))", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
