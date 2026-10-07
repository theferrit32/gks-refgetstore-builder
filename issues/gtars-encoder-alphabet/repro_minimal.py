#!/usr/bin/env python
"""Round-trip test cases for gtars' RefgetStore encoded storage mode.

Adds a small set of sequences to a RefgetStore, reads them back, and reports
where the returned bytes differ from the input. See README.md in the same
directory for an explanation of each case.

Tested with gtars 0.9.2 and 0.10.0, which behave identically. Depends only on
gtars: it needs no existing store, downloaded data or network access.

    uv run --with gtars==0.10.0 python repro_minimal.py

Three stages, each printed as a table:

  1. detection  -- ``digest_sequence`` on each case, showing the alphabet gtars
                   selects and the sequence digest. No store is involved.
  2. round trip -- adds the cases to ``RefgetStore.in_memory()`` at the default
                   ``StorageMode.Encoded``, reads each back, and re-digests the
                   returned bytes with gtars' own ``sha512t24u_digest``.
  3. control    -- the same round trip at ``StorageMode.Raw``. Every case
                   round-trips here, which places the differences in encoding
                   rather than in FASTA parsing.

A final per-residue check then tests individual residues one at a time. It is
informational and does not affect the exit code.

Exit status is 0 when every result matches the expectations in ``CASES``, and
non-zero when any result differs from them -- including a case that now
round-trips, so the script can also be used to confirm a fix.
"""

from __future__ import annotations

import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

from gtars.refget import RefgetStore, digest_sequence, sha512t24u_digest

# ENSP00000473614 (GPX4, glutathione peroxidase 4), 180 aa, from Ensembl
# release 76 pep. Included inline so the script is self-contained; its digest is
# checked below to confirm the sequence was copied correctly. Selenocysteine (U)
# is at index 109.
GPX4 = (
    "MGRAGAGSPGRRRQRCQSRGRRRPRAPRRRKAPACRRRRARRRRKKPCPRSLRPEIHECPK"
    "SQDPCASRDDWRCARSMHEFSAKDIDGHMVNLDKYRGFVCIVTNVASQUGKTEVNYTQLVD"
    "LHARYAECGLRILAFPCNQFGKQEPGSNEEIKEFAAGYNVKFDMFSKICVNGDDAHPL"
)
GPX4_DIGEST = "z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV"

# (key, sequence, expected alphabet, expected to differ after the round trip,
#  description)
CASES = [
    (
        "selenoprotein",
        GPX4,
        "protein",
        True,
        "GPX4 ENSP00000473614; U (selenocysteine) at index 109",
    ),
    (
        "protein-as-dnaio",
        "MCVTYHNGTGYC",
        "dnaio",
        True,
        "NOTCH2 ENSP00000499040.1; every residue is also an IUPAC nucleotide code",
    ),
    (
        "nucleotide-iupac",
        "ACGTHACGTVACGTNACGT",
        "dnaio",
        True,
        "nucleotide sequence with IUPAC ambiguity codes",
    ),
    (
        "pyrrolysine",
        "MKWVTFISLLOFLFSSAYS",
        "ASCII",
        False,
        "control: O selects the ASCII alphabet",
    ),
    (
        "plain-protein",
        "MKWVTFISLLFLFSSAYS",
        "protein",
        False,
        "control: 20 standard residues only",
    ),
    (
        "plain-dna",
        "ACGTACGTACGT",
        "dna2bit",
        False,
        "control: unambiguous nucleotide",
    ),
]

KEY_WIDTH = max(len(key) for key, *_ in CASES)


def charset(seq: str) -> str:
    """Distinct residues in ``seq``, sorted."""
    return "".join(sorted(set(seq)))


def first_difference(given: str, returned: str) -> tuple[int, str, str] | None:
    """``(index, input residue, returned residue)`` of the first difference."""
    for index, (want, got) in enumerate(zip(given, returned, strict=False)):
        if want != got:
            return index, want, got
    return None


def round_trip(sequences: dict[str, str], *, encoded: bool) -> dict[str, str]:
    """Add ``sequences`` to a store and read each one back.

    ``RefgetStore.in_memory`` writes nothing to disk and the FASTA lives in a
    temporary directory, so running this leaves no files behind.
    """
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "cases.fa"
        fasta.write_text(
            "".join(f">{name}\n{seq}\n" for name, seq in sequences.items()),
            encoding="ascii",
        )
        store = RefgetStore.in_memory()
        store.set_quiet(True)
        if not encoded:
            store.disable_encoding()
        store.add_sequence_collection_from_fasta(str(fasta))
        return {
            name: store.stream_sequence(sha512t24u_digest(seq.encode())).read_all()
            for name, seq in sequences.items()
        }


def stage_1_detection() -> list[str]:
    """Alphabet selection and digests, from ``digest_sequence`` alone."""
    print("stage 1 -- detection (pure functions, no store)")
    print(f"  {'case':<{KEY_WIDTH}}  {'len':>4}  {'alphabet':<9}  {'sha512t24u':<32}  expect")
    discrepancies = []
    for key, seq, want_alphabet, _, _ in CASES:
        record = digest_sequence(seq.encode())
        alphabet = str(record.metadata.alphabet)
        flag = "" if alphabet == want_alphabet else f"  <- expected {want_alphabet}"
        if flag:
            discrepancies.append(
                f"stage 1: {key}: alphabet {alphabet}, expected {want_alphabet}"
            )
        print(
            f"  {key:<{KEY_WIDTH}}  {len(seq):>4}  {alphabet:<9}  "
            f"{record.metadata.sha512t24u:<32}  {want_alphabet}{flag}"
        )
    if digest_sequence(GPX4.encode()).metadata.sha512t24u != GPX4_DIGEST:
        raise SystemExit(
            f"the inline GPX4 sequence does not digest to {GPX4_DIGEST}; "
            "it was not copied correctly, so the results below would not apply"
        )
    print(f"\n  The GPX4 sequence digests to {GPX4_DIGEST}, matching ENSP00000473614.")
    return discrepancies


def stage_round_trip(stage: str, heading: str, *, encoded: bool) -> list[str]:
    """One round trip at a given storage mode, printed as a table.

    At ``StorageMode.Encoded`` each case is compared with its own expectation.
    At ``StorageMode.Raw`` nothing is encoded, so every case is expected to
    round-trip.
    """
    print(f"\n{stage} -- {heading}")
    sequences = {key: seq for key, seq, *_ in CASES}
    returned = round_trip(sequences, encoded=encoded)

    print(f"  {'case':<{KEY_WIDTH}}  {'expect':<7}  {'actual':<7}  substitution")
    discrepancies = []
    changed = []
    for key, seq, _, case_differs, _ in CASES:
        want_differs = case_differs and encoded
        payload = returned[key]
        difference = first_difference(seq, payload)
        differs = payload != seq
        want = "differs" if want_differs else "matches"
        got = "differs" if differs else "matches"
        detail = ""
        if difference:
            index, input_residue, returned_residue = difference
            detail = f"index {index}: {input_residue} -> {returned_residue}"
        elif differs:
            # Same prefix, different length.
            detail = f"length {len(seq)} -> {len(payload)}"
        if detail:
            changed.append((key, seq, payload, detail))
        mark = "" if differs == want_differs else "  <- unexpected"
        if mark:
            if want_differs:
                why = "this case now round-trips; the behaviour described in README.md may have changed"
            elif not encoded:
                why = "this sequence does not round-trip at StorageMode.Raw"
            else:
                why = "a control sequence no longer round-trips"
            discrepancies.append(f"{stage}: {key}: expected {want}, got {got}; {why}")
        print(f"  {key:<{KEY_WIDTH}}  {want:<7}  {got:<7}  {detail}{mark}")

    for key, seq, payload, detail in changed:
        print(f"\n  {key}  ({detail})")
        print(f"    input charset     {charset(seq)}")
        print(f"    returned charset  {charset(payload)}")
        print(f"    input digest      {sha512t24u_digest(seq.encode())}   <- the stored digest")
        print(f"    returned digest   {sha512t24u_digest(payload.encode())}")
    return discrepancies


def per_residue_check() -> None:
    """Which residues round-trip, tested one at a time.

    Informational only; it does not affect the exit code. Each residue is
    placed in a context that keeps alphabet selection on the alphabet being
    tested.
    """
    print("\nper-residue check -- each residue in a fixed context, StorageMode.Encoded")
    probes = {
        "nucleotide": ("ACGTACGT{}ACGTACGT", "RYSWKMBDHVNU"),
        "protein": ("MKWVTFISLLFLFSSAYS{}MKWVTFISLLFLFSSAYS", "BJOUXZ"),
    }
    print(f"  {'context':<10}  {'residue':<7}  {'alphabet':<9}  result")
    for label, (template, residues) in probes.items():
        sequences = {f"{label}-{residue}": template.format(residue) for residue in residues}
        returned = round_trip(sequences, encoded=True)
        for residue in residues:
            name = f"{label}-{residue}"
            given = sequences[name]
            alphabet = str(digest_sequence(given.encode()).metadata.alphabet)
            difference = first_difference(given, returned[name])
            result = f"returned as {difference[2]}" if difference else "round-trips"
            print(f"  {label:<10}  {residue:<7}  {alphabet:<9}  {result}")
    print(
        "\n  D and H are encoded with distinct codes, which the decoding table maps\n"
        "  to H and V. U (in both contexts) and B (protein) share a code with T or\n"
        "  A. Residues that select the ASCII alphabet are stored byte for byte."
    )


def main() -> int:
    print(f"gtars {version('gtars')} -- RefgetStore encoded-storage round trips\n{'=' * 62}\n")
    discrepancies = stage_1_detection()
    discrepancies += stage_round_trip(
        "stage 2", "round trip through RefgetStore (StorageMode.Encoded, the default)", encoded=True
    )
    raw_discrepancies = stage_round_trip(
        "stage 3", "the same round trip at StorageMode.Raw (control)", encoded=False
    )
    if raw_discrepancies:
        print(
            "\n  Some sequences also differ at StorageMode.Raw, so the differences\n"
            "  are not limited to encoding."
        )
    else:
        print(
            "\n  Every case round-trips at StorageMode.Raw, including the three that\n"
            "  differ at StorageMode.Encoded. The same bytes reach the store in both\n"
            "  modes, so the differences arise in encoding, not in FASTA parsing."
        )
    discrepancies += raw_discrepancies

    per_residue_check()

    print(f"\n{'=' * 62}")
    if discrepancies:
        print(f"{len(discrepancies)} result(s) differ from the expectations in README.md:")
        for discrepancy in discrepancies:
            print(f"  - {discrepancy}")
        return 1
    print("results match README.md: 3 cases differ after the round trip, 3 controls match.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
