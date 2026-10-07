#!/usr/bin/env python
"""Check whether correcting the ``DnaIupac`` decoding table recovers stored data.

gtars' ``DnaIupac`` encoding table gives D, H and V three distinct codes
(0b1101, 0b1110, 0b1111). The decoding table maps 0b1101 to H and 0b1110 to V.
Because the codes are distinct, the stored bits still identify the original
residue, and changing those two decoding entries should return the original
sequence without re-importing it. See README.md in the same directory.

This tests that against an existing store. For every ``dnaio`` sequence it
reads the ``.seq`` payload from disk, decodes it with two tables, and compares
the digest of each result with the stored digest:

  current    -- gtars 0.10.0's table as published. Checked against the bytes
                gtars itself returns, to confirm this decoder matches gtars.
  corrected  -- the same table with 0b1101 -> D and 0b1110 -> H.

Requires the gks-refgetstore-builder repository
(https://github.com/theferrit32/gks-refgetstore-builder) and a built store:

    uv run python issues/gtars-encoder-alphabet/check_decoder_fix.py --store store

Sequences that round-trip with the current table are included, to confirm the
corrected table leaves them unchanged.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from gks_refgetstore import store_census

# DNA_IUPAC_DECODING_ARRAY, gtars-refget/src/digest/alphabet.rs:234-253 at
# gtars-python-v0.10.0 (90141b68), indexed by 4-bit code.
CURRENT = "NACMGRSKTWYDBHVV"
CORRECTED = CURRENT[:0b1101] + "DH" + CURRENT[0b1111:]
assert (CURRENT[0b1101], CURRENT[0b1110]) == ("H", "V")
assert (CORRECTED[0b1101], CORRECTED[0b1110], CORRECTED[0b1111]) == ("D", "H", "V")


def unpack(payload: bytes, length: int, table: str) -> bytes:
    """Decode 4-bit MSB-first codes, as ``SequenceEncoder`` packs them.

    The ``dnaio`` set includes whole chromosomes (NC_000001.11 is 249 Mbp), so
    each nibble is mapped with ``bytes.translate`` and the halves interleaved by
    slice assignment -- both C loops -- rather than one Python call per base.
    """
    high = bytes(ord(table[b >> 4]) for b in range(256))
    low = bytes(ord(table[b & 0x0F]) for b in range(256))
    out = bytearray(2 * len(payload))
    out[0::2] = payload.translate(high)
    out[1::2] = payload.translate(low)
    return bytes(out[:length])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--store", type=Path, required=True)
    args = parser.parse_args()

    from gtars.refget import RefgetStore, sha512t24u_digest

    store = RefgetStore.open_local(str(args.store))
    store.set_quiet(True)
    metadata = store_census.sequence_metadata_by_digest(store)
    dnaio = sorted(d for d, m in metadata.items() if str(m.alphabet) == "dnaio")

    outcome: Counter[tuple[bool, bool]] = Counter()
    differs_from_gtars = []
    for digest in dnaio:
        meta = metadata[digest]
        payload = store_census.seq_path(args.store, digest).read_bytes()
        current = sha512t24u_digest(unpack(payload, int(meta.length), CURRENT))
        corrected = sha512t24u_digest(unpack(payload, int(meta.length), CORRECTED))
        # Compared by digest: gtars' output is streamed and hashed in chunks,
        # so a chromosome is never held in memory twice.
        if current != store_census.redigest_sequence(store, digest):
            differs_from_gtars.append(digest)
        outcome[(current == digest, corrected == digest)] += 1

    total_bp = sum(int(metadata[d].length) for d in dnaio)
    print(f"{len(dnaio)} dnaio sequence(s), {total_bp:,} bp, in {args.store}")
    if differs_from_gtars:
        print(f"stopping: the current table gives different output from gtars for "
              f"{len(differs_from_gtars)} sequence(s), so this decoder does not match "
              "gtars and the comparison below would not be meaningful")
        return 2
    print("the current table reproduces gtars' output for all of them\n")
    print(f"  {'current round-trips':<21}{'corrected round-trips':<23}sequences")
    for (ok_current, ok_corrected), n in sorted(outcome.items()):
        print(f"  {ok_current!s:<21}{ok_corrected!s:<23}{n}")

    changed = sum(n for (cur, cor), n in outcome.items() if cur and not cor)
    not_recovered = sum(n for (_, cor), n in outcome.items() if not cor)
    if changed or not_recovered:
        print(f"\n{not_recovered} sequence(s) do not round-trip with the corrected table; "
              f"{changed} of them round-trip with the current table")
        return 1
    print("\nevery dnaio sequence round-trips with the corrected table")
    return 0


if __name__ == "__main__":
    sys.exit(main())
