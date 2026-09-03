# seqrepo equivalence

Analysis tooling for one question: **how does this RefgetStore relate to a
biocommons seqrepo snapshot?**

This is an analysis aid, not part of the build contract. The builder loads
authoritative current sources; a seqrepo snapshot holds historical versions and
some sequences it cannot name. The two are expected to differ, and the tools
here describe how.

## Interpreting the numbers

**Comparison is by sha512t24u digest.** The digest is the sequence, so a
sequence counts as shared regardless of what either side calls it. Alias
agreement is a separate question from sequence membership; a sequence held under
a different alias spelling is still held.

**Two coverage figures, and they differ a lot.** Of seqrepo's 1,144,093 digests,
980,664 are in the store — 85.715%. Of the 163,429 that are not, 128,092 are
*digest-only*: seqrepo stores the bytes and their MD5/SEGUID/SHA1/VMC
re-encodings but no biological accession, so no FTP source corresponds to them
and no build can load them. Against the 1,016,001 seqrepo sequences that do
carry an accession, coverage is **96.522%**. The remaining shortfall is 35,286
sequences, 3.1% of seqrepo, mostly NCBI predicted models (`XR_`, `XM_`, `XP_`)
which NCBI renumbers on every annotation release.

**Digest mismatches are measured against the rolling `ensembl` namespace**,
which tracks the newest release. Ensembl sometimes republishes a sequence under
an unchanged `accession.version`, and a change that is later reverted leaves no
trace in that comparison. Mismatch counts therefore set a floor on drift, not a
measure of it — [`ENSEMBL_RELEASE_DRIFT.md`](ENSEMBL_RELEASE_DRIFT.md) measures
it release-to-release instead. Validate pinned artifacts against a pinned
`ensembl-N` namespace, never against `ensembl`.

## Why the namespaces differ

Neither side's namespace list is a subset of the other:

| | namespaces | what they are |
|---|---|---|
| seqrepo only | `MD5`, `SEGUID`, `SHA1`, `VMC` | digest-synthesis namespaces — 1:1 re-encodings of the sequence digest, derivable on demand. The store omits them by design; digest-level coverage already accounts for these sequences. |
| seqrepo only | `NCBI34/35/36`, `hs37d5`, `hs37-1kg`, `JRGv1/2`, `CHM1_1.1`, older GRCh3x patches | assemblies the build does not load. Their sequences are often still present by digest. |
| store only | `ensembl-75` … `ensembl-116` | one immutable namespace per Ensembl release. seqrepo exposes a single `Ensembl` namespace, so there is no counterpart. |
| store only | `lrg` | EBI LRG. The `2024-12-20` snapshot has no LRG namespace, so `NAMESPACE_MAP` omits an `LRG` entry and `lrg` is expected to land entirely in `refget_only`. |

seqrepo's `Ensembl` namespace holds the latest version it has recorded for each
accession, which is not the same as one Ensembl release: accessions loaded at
different times sit at different releases. In `2024-12-20`, chromosome `Y`
matches store namespaces `ensembl-79`–`109` while the proteins match `ensembl-113`.

## Which tool answers which question

| question | tool |
|---|---|
| what does each side hold that the other does not, and why do the namespaces differ | `compare_inventories.py` |
| PASS/FAIL backwards-compatibility verdict, every difference categorized | `verify_seqrepo_equivalence.py` |
| coverage percentages and per-prefix (`NM_`/`NP_`/`XM_`…) rollups | `parity_membership.py` |
| per-alias detail, and which source first contributed each store-only digest | `full_parity.py` |
| how much of the gap a candidate source FASTA would close | `probe_source_coverage.py` |
| how often Ensembl changed a sequence without changing its version | `ENSEMBL_RELEASE_DRIFT.md` |
| why releases 76–109 hold alt/patch scaffolds at chromosome length, and what dropping them would cost | `ENSEMBL_N_PADDED_SCAFFOLDS.md` |

### Contents

- `compare_inventories.py` — stepwise inventory comparison. Names what it is
  comparing at each step and writes complete TSV listings: namespace
  inventories for both sides, a namespace correspondence table with a reason per
  row, the full sequence difference in both directions, and alias-level
  agreement within corresponding namespaces. Nothing is sampled or truncated.
- `verify_seqrepo_equivalence.py` — categorizes every difference (`matched`,
  digest mismatch, `superseded_old_version`, `backfill_candidate`,
  `alias_naming_gap`, `sequence_missing`) into a JSON report plus a
  `build_gaps.tsv` feed, ending with a PASS/FAIL verdict.
- `parity_membership.py` — digest-membership table and Markdown summary; the
  source of the coverage percentages and per-prefix rollups.
- `full_parity.py` — alias-union table, missing and mismatched seqrepo aliases,
  shared-alias mismatches, earliest-source attribution per store-only digest.
- `probe_source_coverage.py` — candidate-source coverage against a gap table.
- `known_divergence/` — regression fixture for Ensembl accessions whose seqrepo
  digest differs from the store's by design, plus its generator. Pinned to the
  `2024-12-20` snapshot and to store namespace `ensembl-113`; see that
  directory's README.
- `ENSEMBL_RELEASE_DRIFT.md` — release-to-release sequence drift measurement.
- `ENSEMBL_N_PADDED_SCAFFOLDS.md` — Ensembl's former N-padded alt/patch scaffold
  representation: 445 records emitted at parent-chromosome length in releases
  76–109, their recoverability, storage cost, and downstream reachability. Backed
  by `ensembl_padded_scaffolds.tsv` (one row per padded record) and
  `ensembl_dna_representation_by_release.tsv` (one row per release 75–116), both
  regenerated by `tools/ensembl_padding_probe.py`.

## Usage

Run from the repo root. Defaults point at `./store`, the `2024-12-20` snapshot,
and `runs/2026-08-26-full-rebuild/parity/`.

    uv run python seqrepo_equivalence/compare_inventories.py \
        --out-dir runs/DATE/inventory
    uv run python seqrepo_equivalence/verify_seqrepo_equivalence.py -v

`compare_inventories.py` takes `--seqrepo`, `--store`, `--out-dir` (required),
and `--include-shared`, and works against any snapshot and any store directory.
For the verifier see `--help`: `--seqrepo`, `--store`, `--known-divergent`,
`--report`, `--gap-list`, `--fail-on {none,gaps,mismatch}`.

Against a different seqrepo snapshot, build a new fixture in `known_divergence/`
for that version rather than editing the existing one. The set of `*`-divergent
accessions changes between snapshots, and Ensembl strips `*` upstream in some
releases, which retires cases entirely.

## Outputs

Bulk TSV outputs are reproducible and gitignored. Retain them locally and record
their schemas, row counts, byte sizes, SHA-256 values, generator commands, and
headline categories in the relevant run manifest.
