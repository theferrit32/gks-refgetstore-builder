# Known Ensembl digest divergence

Why some Ensembl accessions legitimately carry a different sha512t24u in this
store than in the biocommons seqrepo `2024-12-20` snapshot, and the fixture that
keeps those cases from being re-reported as unexplained failures.

## Files

- **`ensembl_vs_seqrepo_digest_divergence.tsv`** — the fixture. 116 rows,
  `accession · seqrepo_seq_id · refget_sha512t24u · cause`. Consumed by
  `verify_store.py`, which validates it against the pinned `ensembl-113`
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
