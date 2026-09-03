# Ensembl's pre-110 N-padded scaffolds inflate the store by ~21 GiB (43%)

**Labels:** `data-quality` · `storage` · `sources` · `resolved`

> **Resolved by option A, with one change.** Option A below proposed an
> N-fraction threshold. The implemented rule matches on **record name** instead —
> Ensembl marks these records with a `CHR_` prefix, which is both sound and
> complete across releases 76–109, so a name rule is exact where a threshold
> would be a guess and would have needed tuning against legitimately N-rich
> sequences (chrY is 55% `N`). Each affected seqset declares
> `exclude = { file_classes = ["dna.toplevel"], record_prefixes = ["CHR_"] }`.
> See the *Outcome* section of
> [`seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md`](seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md).
> Acceptance criterion 1 (a documented threshold) is therefore not applicable;
> the rest hold. Pinning the ingested outcome in the lock is tracked separately
> in [`ISSUE-lock-should-pin-ingest-outcome.md`](ISSUE-lock-should-pin-ingest-outcome.md).

## Summary

Ensembl releases 76–109 publish alt loci and patch scaffolds in `dna.toplevel`
not at their true length, but padded out to roughly the length of their parent
chromosome, with `N` filling everything outside the placed interval. Release 110
(FASTA staged 2023-04-21) switched to true-length records.

The store ingests releases 75–116, so it holds **445 padded records totalling
60,047,447,030 bases**. Because `N` cannot be represented in 2 bits, these are
stored in the `dna3bit` alphabet at 3 bits per base: **20.97 GiB, or 43% of the
48.9 GiB store**, encoding no sequence the store does not already hold at true
length.

## What a padded record looks like

`HSCHR1_2_CTG3` is a 256 kb alternate haplotype of the PRAME region of
chromosome 1. In release 100 it appears as `CHR_HSCHR1_2_CTG3`:

| | total length | real (non-`N`) bases | |
|---|---:|---:|---|
| `ensembl-100:CHR_HSCHR1_2_CTG3` | 248,975,002 | 256,271 | **0.10%** |
| `ensembl-100:CHR_HG708_PATCH` | 159,270,839 | 589,656 | 0.37% |
| `ensembl-100:1` (a real chromosome, for contrast) | 248,956,422 | 230,481,012 | 92.58% |
| `ensembl-116:HSCHR1_2_CTG3` (same scaffold, modern form) | 256,271 | 256,271 | 100% |

The padding encodes the scaffold's placement: the real bases start at the
chromosome offset NCBI records as `parent_start` in
`alt_scaffold_placement.txt`, so offset *i* in the padded record equals offset
*i* on the parent chromosome. Padding was applied only to sequences that have a
parent-chromosome placement — the 127 unplaced and 42 unlocalized scaffolds stay
at true length in every release.

## Why it is expensive here

Sequence storage is `Encoded`. `dna2bit` holds ACGT at 2 bits per base;
any `N` forces `dna3bit` at 3 bits per base. A padded record is >99% `N`, so it
pays full price for padding.

| alphabet | sequences | bases |
|---|---:|---:|
| dna2bit | 1,426,421 | 4,972,744,406 |
| protein | 351,881 | 184,072,924 |
| **dna3bit** | **1,088** | **98,911,237,632** |
| dnaio | 107 | 1,823,846,320 |

`dna3bit` is 1,088 sequences but **93% of all bases in the store**.

## The padded content is fully redundant

- All **445** distinct `CHR_`-prefixed names have a true-length counterpart
  already in the store. Records whose content exists only in padded form: **0**.
- Extracting the real span from each padded record and re-hashing reproduces the
  true-length digest exactly for all 445 — 387 forward, 58 reverse-complemented
  (exactly those with placement orientation `-`), 0 failures.
- The padded records carry **no alias outside `ensembl-76`–`ensembl-109`**, and
  no `refseq` or `insdc` accession. They are an Ensembl FASTA presentation, not
  accessioned sequence entities. NCBI has published `HSCHR1_2_CTG3` at 256,271
  bases in all 15 assembly reports containing it, from GRCh38 through p14.

## Who consumes alt/patch sequences

Alt and patch scaffolds matter, but consumers address them by accession, which
only ever resolves to the true-length form.

- **ClinVar** (`ClinVarVCVRelease_2026-0621`, full scan of 1,184,400,464 lines):
  828,074 alt/patch references against 47,745,798 chromosome references —
  **1.70%**, across 502 distinct accessions. Every one is ≤ 6,530,008 bases,
  i.e. all true-length; none padded. All 502 resolve in the store today.
- **ClinGen Allele Registry** accepts HGVS on `NT_187517.1` and `KI270766.1` and
  validates against real sequence.
- **seqrepo 2024-12-20** holds all 445 padded digests, but reachable through
  exactly one biological namespace, `Ensembl`, under the `CHR_` names. The
  true-length form additionally carries `NCBI:NT_187517.1`, `GRCh38`,
  `GRCh38.p1`–`p12`, and UCSC-style `chr1_KI270766v1_alt`.

## Impact

- Store is 48.9 GiB / 1,779,811 objects; ~21 GiB and no additional objects are
  padding.
- Publication cost scales with it: ~1.78M PUTs and 48.9 GiB per uploaded
  snapshot, against 9.2 GiB for the 2026-07-22 store.
- Every full rebuild decodes and re-encodes the padding.

## Options

| option | store size | trade-off |
|---|---:|---|
| **A. Filter near-100%-`N` records at ingest** | 48.9 → **27.9 GiB** | Removes only the padding. Keeps `ensembl-75`–`109` chromosome aliases and real scaffolds. Requires a record-level filter in `build_store.py`, a rebuild, and a new lock. Threshold must not catch legitimately `N`-rich sequences — chrY is 55% `N` and some `ENST` records exceed 90%. |
| **B. Drop `dna.toplevel` for releases ≤109** | 48.9 → **15.6 GiB** | `sources.toml` edit only, no new build logic. Also discards **11.19 GiB of Ensembl GRCh37 genomic sequence** (release 75) and 1.10 GiB of superseded GRCh38 chromosome variants, which are not padding. Loses all genomic aliases for releases 75–109. |
| **C. Leave as-is** | 48.9 GiB | Faithful to upstream; every digest is a correct refget digest of published bytes. No rebuild. |

Under A or B the 445 padded digests move from "in both" to "seqrepo-only" in the
parity report: digest coverage 85.715% → 85.677%, gap 163,429 → 163,874.

## Suggested acceptance criteria

If option A is taken:

1. A documented, configurable threshold, with the `N`-fraction distribution
   across all genomic records recorded to justify the cutoff.
2. The filter reports every skipped record by name, length, and `N` fraction, so
   exclusions are auditable rather than silent.
3. `verify_store.py` still passes, including its namespace-set equality check.
4. A regenerated `build.lock.json`, and a parity re-run confirming the expected
   85.677% and no other movement.
5. `seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md` updated with the outcome.

## Evidence

- Report: [`seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md`](seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md)
- Per-record data: `seqrepo_equivalence/ensembl_padded_scaffolds.tsv` (445 rows × 27 columns)
- Per-release data: `seqrepo_equivalence/ensembl_dna_representation_by_release.tsv`
- Regenerate both: `uv run python tools/ensembl_padding_probe.py`

Smallest reproduction:

```python
from gtars.refget import RefgetStore
s = RefgetStore.open_local("store")
for ns, name in (("ensembl-100", "CHR_HSCHR1_2_CTG3"),
                 ("ensembl-116", "HSCHR1_2_CTG3")):
    m = s.get_sequence_metadata_by_alias(ns, name)
    print(ns, name, f"length={m.length:,}", str(m.alphabet))
```
