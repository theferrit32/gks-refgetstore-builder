# 2026-08-26 full rebuild

Clean end-to-end rebuild: download cache, store, and lock all regenerated from
scratch, then the seqrepo parity analysis run against the result. Supersedes the
removed `2026-07-02-seqrepo-parity` record.

- Store: **1,779,497** sequences, **240** collections, 68 alias namespaces,
  **0** seqset warnings.
- Locked inputs: **335** sources in the repo-root `build.lock.json`. The 23 with
  a null `collection_digest` are `assembly_report` files, which are not
  collections.
- Ensembl: releases **75–116**, each pinned to its own `ensembl-N` namespace,
  with the rolling `ensembl` namespace tracking 116.
- New this build: EBI LRG (`lrg`, 5,905 aliases, 442 store-only digests).
- Verification: `verify_store.py` PASS; cache-vs-lock `matched=335 changed=0
  new=0 missing=0 corrupt=0`, set-diff 0/0; `pytest` 73 passed.

## Commands

Recorded in [`run-manifest.json`](run-manifest.json) with logs under `logs/`.
Two entries have reconstructed rather than captured exit codes, and say so in
`evidence.missing`: `fetch` (serial, superseded by `fetch-parallel`) and `build`
(killed mid-ingest at ensembl-89 to pick up batched ingest; superseded by
`build-resume`, which reused this run's completed collections). `run_record`'s
`finally` block does not run under SIGKILL, so their `exit_code` and
`finished_utc` were derived from log mtimes afterwards.

## Parity against seqrepo 2024-12-20

Outputs in [`parity/`](parity/); the bulk TSVs are gitignored and regenerable.

| metric | value |
|---|---:|
| digests in both | 980,664 |
| seqrepo-only | 163,429 |
| store-only | 798,833 |
| coverage, all digests | 85.715% |
| **coverage of seqrepo sequences that carry an accession** | **96.522%** |

**Read the second coverage figure, not the first.** 128,092 of the 163,429
missing digests — 78% — are *digest-only*: seqrepo holds the bytes and their
MD5/SEGUID/SHA1/VMC re-encodings but no biological accession anywhere. Nothing
on any FTP server corresponds to them, so no build can close that gap; they
would have to be copied byte-for-byte out of seqrepo. The actionable shortfall
is **35,286** sequences (3.1% of seqrepo), dominated by NCBI's predicted-model
tail (`XR_` 25,203, `XM_` 4,136, `XP_` 1,662), which NCBI renumbers on every
annotation release — a fresh seqrepo build from today's FTP would hit the same
wall. Full decomposition in
[`parity/parity-summary.md`](parity/parity-summary.md).

## Inventory comparison

[`inventory/`](inventory/) holds the complete, unsampled listings of what each
side has that the other does not, written by
`seqrepo_equivalence/compare_inventories.py`:

    uv run python seqrepo_equivalence/compare_inventories.py \
        --out-dir runs/2026-08-26-full-rebuild/inventory

| file | contents | rows |
|---|---|---:|
| `00_summary.tsv` | every metric below, as `step · metric · value` | 16 |
| `01_seqrepo_namespaces.tsv` | seqrepo namespaces, alias rows, distinct digests | 35 |
| `02_store_namespaces.tsv` | store namespaces, alias rows, distinct digests | 69 |
| `03_namespace_comparison.tsv` | namespace correspondence, one reason per row | 81 |
| `04_sequences_seqrepo_only.tsv` | every seqrepo-only digest, with length, alphabet, and whether it has any biological accession | 163,429 |
| `05_sequences_store_only.tsv` | every store-only digest, with length, alphabet, aliases | 798,833 |
| `07_alias_comparison_by_namespace.tsv` | per corresponding namespace: shared aliases, digest agreement, one-sided aliases | 23 |
| `08_alias_differences.tsv` | every differing alias, with both digests | 1,003,010 |

Namespaces: 23 correspond, 12 are seqrepo-only (`MD5`, `SEGUID`, `SHA1`, `VMC`
plus `NCBI34/35/36`, `hs37d5`, `hs37-1kg`, `JRGv1/2`, `CHM1_1.1`), and 46 are
store-only (42 pinned `ensembl-N` releases, `insdc`, `lrg`, `GRCh38.p13`,
`GRCh38.p14`). The rolling `ensembl` namespace corresponds to seqrepo's
`Ensembl`. Alias-level agreement within corresponding namespaces is total apart
from Ensembl: `NCBI`→`refseq` agrees on all 826,115 shared aliases, and every
assembly namespace agrees on all of theirs.

Files 04, 05, and 08 are bulk and gitignored; the rest are committed. The whole
set regenerates from the store plus the snapshot in roughly eight minutes.

## Ensembl silent release drift

The parity run reported 317 shared-alias digest mismatches. 115 are the
`*`-stripping artifact; the other 202 led to
[`../../seqrepo_equivalence/ENSEMBL_RELEASE_DRIFT.md`](../../seqrepo_equivalence/ENSEMBL_RELEASE_DRIFT.md),
which measures how often Ensembl republishes the same `accession.version` with
different bytes. Two conclusions bear on how this record should be read:

- The store is faithful to upstream — verified by rebuilding the affected
  accessions from the cached FASTAs alone (11,636/11,636 digests identical).
- The mismatch count is a **weak** drift metric. It compares only against the
  rolling `ensembl` namespace, so it misses anything Ensembl reverted before the
  newest release: it surfaced 200 of 4,407 actual silent changes.
