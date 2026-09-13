# Full round-trip scan: the gtars encoder baseline

Every sequence in the store, streamed out, re-digested from the bytes the store
actually serves, and compared to the digest it is filed under.

This is the check the build lock could never make. The lock records digests,
counts, and collection membership, and all of those can be perfectly correct
while the *payload* is wrong — which is exactly what was found.

## Result

**133 of 1,779,052 sequences (0.0075%) do not re-digest to their own digest.**

| alphabet | sequences | round-trip failures |
|---|---:|---:|
| dna2bit | 1,426,421 | 0 |
| dna3bit | 643 | 0 |
| **dnaio** | 107 | **27** |
| **protein** | 351,881 | **106** |

Scanned at 2,778 sequences/second, single-threaded: **10m40s** for the full
store. That is the number that made `--deep` a routine check rather than an
overnight job.

## Two stores, one baseline

`store.pre-filter/` (1,779,052 sequences) and `store/` (1,779,497) were scanned
independently. The output TSVs are **byte-identical**. The 445-sequence
difference is entirely N-padded `dna3bit` scaffolds — 1,088 in `store/` against
643 — and every one of them round-trips cleanly. One baseline therefore serves
both stores, and the padding exclusion is confirmed to have had no effect on
this defect class either way.

## The defects

Neither is bit rot. Each digest is **correct** — computed from the source bytes
at ingest. Each payload is **wrong**, because gtars' encoder cannot represent
every residue it was given. Re-ingesting from the same source FASTA reproduces
every row byte for byte, which is why `repair` refuses them.

### `protein` — no `U` (selenocysteine), 106 rows

Selenocysteine is a standard proteinogenic amino acid and appears throughout
Ensembl's `pep` distributions. The alphabet lacks it and the residue is replaced.

Confirmed against the source FASTA, not inferred:

```
ENSP00000473614  (GPX4), release-76 pep, 180 aa
  source charset            ACDEFGHIKLMNPQRSTUVWY
  decoded charset           ACDEFGHIKLMNPQRSTVWY     <- no U
  index 109                 source U  ->  decoded A
  sha512t24u(source)        z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV   <- the store's key
  sha512t24u(decoded)       4otSS66T9mEsIeu_jmydFFDB_x-LFqsb
```

The source digest and the store's key agree exactly, which is the whole point:
ingest hashed the right bytes and then stored different ones.

### `dnaio` — misdetected alphabet, then a lossy round trip, 27 rows

A short protein whose residues all happen to be legal IUPAC nucleotide codes is
detected as nucleotide. That alone would be harmless if the encoding were
lossless; it is not.

```
ENSP00000499040.1  (NOTCH2), release-113 pep, 12 aa
  source                    MCVTYHNGTGYC
  decoded                   MCVTYVNGTGYC
                                 ^ index 5: H -> V
  sha512t24u(source)        13xx4yn7JousnAXLxTBX-bbWg4TitdN0   <- the store's key
  sha512t24u(decoded)       fQh3MU3P7sOQ7EihvRecrb08e86fBG0G
```

Every residue (`M C V T Y H N G T G Y C`) is a legal IUPAC nucleotide code, so
the misdetection is at least understandable for a 12-character input. The
corruption is the second half: `H` (A/C/T) comes back as `V` (A/C/G), so `dnaio`
does not round-trip its own alphabet. That is arguably the more serious of the
two defects, because it affects genuine nucleotide sequences carrying ambiguity
codes, not only short proteins that resemble them.

Affected lengths run from 3 to 29 residues.

## What was produced

- **`roundtrip.tsv`** — the scan output, 133 rows,
  `digest · name · alphabet · length · redigest · cause`.
- **`seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv`** — the
  same rows with an explanatory header, checked in as the `verify --deep`
  baseline.
An upstream gtars issue is drafted but **not filed, and deliberately not
tracked** — it is a local working file, not part of this record. Everything it
would say about the defects is above, in "The defects", including both
reproducers diffed against their source FASTA records.

## The `cause` column is a diagnosis

It cannot be inferred at verify time: an encoder defect and bit rot are
indistinguishable from a digest mismatch alone. The two causes here were
established by pulling each reproducer's record out of its source FASTA and
diffing it against the decoded payload, then confirming the pattern held across
the group. `scan_roundtrip.py` assigns causes by alphabet because that is what
separates the groups mechanically — the diagnosis behind that mapping is in the
script's header comment, and any future row landing as `undiagnosed` needs the
same treatment rather than a guess.

## Why a baseline rather than a suppression flag

A flag has two outcomes; a baseline has three, and the third is the one worth
having:

| in baseline | re-digests correctly | outcome |
|---|---|---|
| yes | no | `warn` — known defect, rendered with its cause |
| no | no | **`error`** — a new defect; this is the regression that matters |
| yes | **yes** | `info` — gtars was fixed |

When the encoder is fixed upstream, `verify --deep` says so instead of the fix
passing unnoticed behind a flag that is still suppressing it. At that point the
affected collections need re-ingesting: their digests were always right, so only
the payloads change.

## What this means for the published store

The store contains 133 sequences it cannot serve correctly. They resolve, return
bytes of the right length, and are wrong at one or more residues. Anything
downstream comparing a selenoprotein against Ensembl will see a mismatch, and
the mismatch will be real.

This is not fixed here. Fixing it means either patching gtars or encoding those
sequences outside its alphabet model, both larger than the lock/verify work, and
neither of which that work depends on.

## Reproducing

```sh
uv run python runs/2026-09-12-deep-verify/scan_roundtrip.py \
    --store store.pre-filter --out roundtrip.tsv
```

Or as the verifier sees it:

```sh
gks-refgetstore verify --store --deep
```
