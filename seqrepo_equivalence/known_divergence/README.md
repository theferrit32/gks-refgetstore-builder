# Known divergence fixtures

Checked-in ground truth for differences that are **explained**, so a verifier can
report them as known rather than re-raising them as unexplained failures — and,
more importantly, so a genuinely new case stands out instead of drowning.

Two independent fixtures live here. They answer different questions, are
consumed by different tools, and share only the convention: a TSV whose header
comments explain what it pins, a `cause` column that is a diagnosis rather than a
measurement, and a generator kept beside it.

| fixture | question | consumed by |
|---|---|---|
| `ensembl_vs_seqrepo_digest_divergence.tsv` | why does this store's digest differ from seqrepo's for the same accession? | `gks-refgetstore-check` |
| `gtars_encoding_roundtrip.tsv` | which stored payloads do not re-digest to the digest they are filed under? | `gks-refgetstore verify --store --deep` |

---

# Ensembl digest divergence vs seqrepo

Why some Ensembl accessions legitimately carry a different sha512t24u in this
store than in the biocommons seqrepo `2024-12-20` snapshot.

## Files

- **`ensembl_vs_seqrepo_digest_divergence.tsv`** — the fixture. 116 rows,
  `accession · seqrepo_seq_id · refget_sha512t24u · cause`. Consumed by
  `gks-refgetstore-check`, which validates it against the pinned `ensembl-113`
  namespace. Read the header comments before changing anything.
- **`diagnose_mismatches.py`** — the generator. Byte-diagnoses every
  shared-alias mismatch emitted by `../full_parity.py` and assigns each a cause.

## The two causes

**`star_normalization` (115 rows).** bioutils' `seq_seqhash` strips `*` before
hashing; gtars hashes the sequence as published. Same sequence, two digests.
Confirmed by re-hashing — strip `*` from the stored sequence, take sha512t24u,
and seqrepo's digest comes back, 115/115.

**`sequence_divergence` (1 row, `ENST00000668831.1`).** Not a hashing artifact:
Ensembl truncated this transcript 1214 → 1120 bases at the 112 → 113 boundary
without bumping `.1`. seqrepo agrees with releases 97–112, this store with
113–116. It is kept here because it was the first observed case of silent
release drift, which is now measured across all 42 releases in
[`../ENSEMBL_RELEASE_DRIFT.md`](../ENSEMBL_RELEASE_DRIFT.md).

## Running the generator

`diagnose_mismatches.py` needs `biocommons.seqrepo`, `bioutils`, and
`ga4gh.core`, which are **deliberately not in `pyproject.toml`** — the build and
verify paths depend on `gtars` alone, and this is a one-off analysis aid, not
part of that contract. Install them into a throwaway environment:

```sh
uv run --with biocommons.seqrepo --with bioutils --with ga4gh.vrs \
    python seqrepo_equivalence/known_divergence/diagnose_mismatches.py --help
```

Regeneration is rarely the right move. The fixture is pinned to one seqrepo
snapshot and one Ensembl release; if either changes, the honest step is a new
fixture alongside this one, not an edit in place.

---

# gtars encoder round-trip failures

Sequences whose stored payload does not re-digest to the digest it is filed
under. Unlike the fixture above, this is not a disagreement between two tools
about how to hash the same sequence. A row here is a **data-correctness
defect**: the store returns bytes that differ from what was imported.

## Files

- **`gtars_encoding_roundtrip.tsv`** — the fixture,
  `digest · name · alphabet · length · redigest · cause`. Consumed by
  `gks-refgetstore verify --store --deep` as its known-bad baseline.
  **Currently empty**: every sequence in the store round-trips.

## Current state

The store is built with gtars from
[databio/gtars#273](https://github.com/databio/gtars/pull/273), pinned as
described in the top-level README's
[gtars revision](../../README.md#gtars-revision) section. A full
`verify --store --deep` against the empty fixture passed with 0 errors over
1,779,052 sequences; see
[`../../runs/2026-10-07-gtars-pr273-rebuild/`](../../runs/2026-10-07-gtars-pr273-rebuild/README.md).

Stores built with gtars 0.10 or earlier return 133 sequences with different
residues than were imported:

| alphabet | sequences | differ |
|---|---:|---:|
| dna2bit | 1,426,421 | 0 |
| dna3bit | 643 | 0 |
| dnaio | 107 | 27 |
| protein | 351,881 | 106 |

- **`protein` (106)** — selenocysteine `U` has no code in the older protein
  alphabet and is stored as alanine's code. `ENSP00000473614` (GPX4,
  release-76 `pep`, 180 aa) has `U` at index 109 and is returned with `A`.
- **`dnaio` (27)** — the older IUPAC decoding table returns `D` as `H` and `H`
  as `V`. 26 are short Ensembl proteins whose residues are all IUPAC nucleotide
  codes, such as `ENSP00000499040.1` (NOTCH2, 12 aa: `MCVTYHNGTGYC` returned as
  `MCVTYVNGTGYC`); one is the RefSeq lncRNA `NR_103745.1`.

That list is preserved as
[`../../runs/2026-10-07-gtars-pr273-rebuild/artifacts/prior-gtars_encoding_roundtrip.tsv`](../../runs/2026-10-07-gtars-pr273-rebuild/artifacts/prior-gtars_encoding_roundtrip.tsv),
and the scan that produced it is
[`../../runs/2026-09-12-deep-verify/`](../../runs/2026-09-12-deep-verify/README.md).
The cause in the gtars source is written up in
[`../../issues/gtars-encoder-alphabet/`](../../issues/gtars-encoder-alphabet/README.md).

## Why it is a baseline and not a suppression flag

A suppression flag has two outcomes. A baseline has three:

| in baseline | re-digests correctly | outcome |
|---|---|---|
| yes | no | `warn` — known defect, rendered with its cause |
| **no** | **no** | **`error`** — a new defect; this is the regression that matters |
| **yes** | **yes** | `info` — gtars was fixed |

The third row is the reason. When an encoder is fixed, the tool says so, rather
than the fix passing unnoticed behind a flag that was still suppressing it. The
affected collections then need re-importing with the fixed gtars: their digests
were always right, so only the payloads change.

With the fixture empty, only the second row can occur: any sequence that does
not round-trip is an `error`.

## Adding rows

Only for a diagnosed case: a new encoder limitation, or a deliberate decision to
publish a store built with a gtars that has one. See the `--deep` baseline
recipe in [`../../RUNBOOK.md`](../../RUNBOOK.md). The `cause` column cannot be
inferred: at verify time an encoder defect and bit rot are indistinguishable, so
causes are written by hand and any `undiagnosed` row is investigated.
