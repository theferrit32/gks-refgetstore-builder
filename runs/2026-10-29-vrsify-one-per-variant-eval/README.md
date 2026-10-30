# vrsify over one expression per ClinVar variation: RefgetStore vs SeqRepo

clinvar-ingest's `vrsify` step translated ClinVar release 2026-10-24's
`variation_identity` (4,692,217 variations, one selected expression each) to GA4GH
VRS objects with the four sequence backends of
[`../2026-10-29-vrsify-all-genomic-eval/`](../2026-10-29-vrsify-all-genomic-eval/README.md),
three rounds per configuration, on the same machine type: with 8 workers for every
backend, and with 1 worker for the remote RefgetStore and SeqRepo.

**Results**

- All 18 runs report the same 4,692,217 records and the same per-class counts
  (outputs were not compared record by record).
- With 1 worker, SeqRepo translates in 1.64× the remote RefgetStore's time
  (1,340.7 s against 818.6 s), 1.75× including setup; the remote store downloads
  2.53 GB, each sequence it touches once, and needs no setup.
- With 8 workers, the RefgetStore arms are within 9% of each other (A 202.5 s,
  B 216.0 s, D 221.2 s); SeqRepo takes 337.6 s (1.67×).
- 8 workers on 8 vCPUs (4 cores × 2 hyperthreads) are 3.8–4.0× as fast as 1.

Times in this record are UTC.

## Input

`gs://clinvar-gkm/2026-10-24/dev/scratch/full/vi-0*.jsonl.gz`: 17 objects,
111,180,465 bytes, the whole of
`clingen-dev.clinvar_2026_10_24_v2_6_1.variation_identity` (4,692,217 rows, its
row count) as newline-delimited JSON. Each provenance file lists every object's
generation and MD5. The `vi-0*` pattern keeps vrsify output written beside the
export out of the input.

## Backends, software and machines

As in [`../2026-10-29-vrsify-all-genomic-eval/`](../2026-10-29-vrsify-all-genomic-eval/README.md):
A is the `genomic-2026-10-29` RefgetStore archive with whole-sequence loads, B the
same store read remotely from `clingen-public` with a cache per worker, D the
archive with range reads, Cg the SeqRepo built from the store. Image
`gcr.io/clingen-dev/clinvar-ingest:vrsify-exp-12774e5`; Cloud Batch `us-east1`,
**c4-highmem-8** (task request 8,000 milli-CPU, 55,296 MiB), 100 GB
`hyperdisk-balanced`, labels `experiment=vrsify-store-comparison`,
`input=one-per-variant`, `round`, `arm`, `workers`, `machine`, `store`.

## Commands

From the clinvar-ingest worktree, with [`artifacts/submit_arms.sh`](artifacts/submit_arms.sh),
for each round `r1`, `r2`, `r3`:

    IN='gs://clinvar-gkm/2026-10-24/dev/scratch/full/vi-0*.jsonl.gz'
    IMG=gcr.io/clingen-dev/clinvar-ingest:vrsify-exp-12774e5
    MACHINE=c4-highmem-8 ARCHIVES=gcs WORKERS=8 LABELS="round=r1 input=one-per-variant" \
        bash submit_arms.sh $IMG "$IN" gs://clinvar-gkm/2026-10-24/dev/scratch/opv-hm-r1 A B D Cg
    MACHINE=c4-highmem-8 ARCHIVES=gcs WORKERS=1 LABELS="round=r1 input=one-per-variant" \
        bash submit_arms.sh $IMG "$IN" gs://clinvar-gkm/2026-10-24/dev/scratch/opv-hm-1w-r1 B Cg

All 18 jobs were submitted at once and ran concurrently on separate VMs.

## Runs

[`artifacts/batch-jobs.txt`](artifacts/batch-jobs.txt) lists them (workers, round,
arm, Batch job). Outputs and provenance are under
`gs://clinvar-gkm/2026-10-24/dev/scratch/opv-hm-<round>-<arm>/` (8 workers) and
`…/opv-hm-1w-<round>-<arm>/` (1 worker); the provenance files are copied here as
`artifacts/provenance-<8w|1w>-<round>-<arm>.json`.

| workers | round | arm | Batch job | running → succeeded (UTC, 2026-10-30) | job run time |
| --- | --- | --- | --- | --- | ---: |
| 8 | r1 | A | `vrsify-2026-10-24-20261030021356` | 02:14:53 → 02:21:03 | 6.2 min |
| 8 | r1 | B | `vrsify-2026-10-24-20261030021359` | 02:14:52 → 02:18:54 | 4.0 min |
| 8 | r1 | D | `vrsify-2026-10-24-20261030021403` | 02:14:52 → 02:21:21 | 6.5 min |
| 8 | r1 | Cg | `vrsify-2026-10-24-20261030021407` | 02:15:00 → 02:22:39 | 7.7 min |
| 1 | r1 | B | `vrsify-2026-10-24-20261030021410` | 02:14:58 → 02:28:45 | 13.8 min |
| 1 | r1 | Cg | `vrsify-2026-10-24-20261030021414` | 02:15:01 → 02:38:28 | 23.4 min |
| 8 | r2 | A | `vrsify-2026-10-24-20261030021417` | 02:15:06 → 02:21:03 | 6.0 min |
| 8 | r2 | B | `vrsify-2026-10-24-20261030021421` | 02:15:12 → 02:19:05 | 3.9 min |
| 8 | r2 | D | `vrsify-2026-10-24-20261030021424` | 02:15:19 → 02:21:24 | 6.1 min |
| 8 | r2 | Cg | `vrsify-2026-10-24-20261030021427` | 02:15:18 → 02:23:16 | 8.0 min |
| 1 | r2 | B | `vrsify-2026-10-24-20261030021431` | 02:15:37 → 02:29:38 | 14.0 min |
| 1 | r2 | Cg | `vrsify-2026-10-24-20261030021434` | 02:15:34 → 02:39:33 | 24.0 min |
| 8 | r3 | A | `vrsify-2026-10-24-20261030021438` | 02:15:28 → 02:21:04 | 5.6 min |
| 8 | r3 | B | `vrsify-2026-10-24-20261030021441` | 02:15:41 → 02:19:29 | 3.8 min |
| 8 | r3 | D | `vrsify-2026-10-24-20261030021444` | 02:15:37 → 02:21:21 | 5.7 min |
| 8 | r3 | Cg | `vrsify-2026-10-24-20261030021448` | 02:15:38 → 02:22:56 | 7.3 min |
| 1 | r3 | B | `vrsify-2026-10-24-20261030021451` | 02:15:41 → 02:29:37 | 13.9 min |
| 1 | r3 | Cg | `vrsify-2026-10-24-20261030021455` | 02:15:47 → 02:40:58 | 25.2 min |

## Results

Every run: 4,692,217 records; Allele 4,598,094 translated and 2,573 errors,
CopyNumberChange 33,184 and 539, CopyNumberCount 49,034 and 95 (the range-copy
records of clinvar-ingest#327), and 8,698 records of untranslated classes
(Unknown 6,872, Not Available 1,154, Haplotype 672) with a null output.

Per configuration, means of 3 rounds (memory is the sum of each worker's peak
RSS; CPU busy is averaged over workers):

| workers | arm | translate s (r1 / r2 / r3) | mean | vs 8-worker A | setup | memory | CPU busy | voluntary ctx sw | downloaded during run | slowest record |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | A | 202.3 / 202.4 / 202.8 | 202.5 s | 1.00× | 131.7 s | 15.9 GiB | 0.944 | 49,411 | — | 48.1 s |
| 8 | D | 221.0 / 221.1 / 221.4 | 221.2 s | 1.09× | 126.5 s | 2.9 GiB | 0.945 | 43,269 | — | 64.5 s |
| 8 | B | 223.1 / 214.5 / 210.4 | 216.0 s | 1.07× | 0.1 s | 15.5 GiB | 0.890 | 408,267 | 14.29 GB | 47.6 s |
| 8 | Cg | 320.4 / 355.7 / 336.7 | 337.6 s | 1.67× | 104.5 s | 1.9 GiB | 0.960 | 21,681 | — | 137.8 s |
| 1 | B | 810.0 / 825.3 / 820.5 | 818.6 s | 4.04× | 0.1 s | 2.8 GiB | 0.978 | 109,107 | 2.53 GB | 24.4 s |
| 1 | Cg | 1,269.2 / 1,334.6 / 1,418.2 | 1,340.7 s | 6.62× | 95.7 s | 0.7 GiB | 1.000 | 11 | — | 78.7 s |

- **One worker.** SeqRepo takes 1.64× B's translation time, 1.75× with setup
  (1,436 s against 819 s). B's worker downloads each sequence it touches once,
  2.53 GB, and waits on those downloads for about 2% of its time; SeqRepo's is
  CPU-bound throughout.
- **Eight workers.** A, B and D are within 9% of each other on this input; each
  of B's 8 workers keeps its own cache (14.3 GB in all against 2.53 GB for one).
  SeqRepo is 1.67× A and varies most between rounds.
- **Scaling.** 8 workers are 3.79× (B) and 3.97× (Cg) as fast as 1 on 8 vCPUs,
  which are 4 cores with 2 hyperthreads each.
- **Slowest records** are deletions in the run of N at the start of GRCh37 chr15
  (`NC_000015.9:g.5185del` and its neighbours), about 20 million single-base reads
  each during normalization: 24 s on B and 79 s on SeqRepo with one worker; with 8
  workers sharing the cores, 48 s on A and B, 65 s on D, 138 s on SeqRepo.
- Setup for the archive arms (A and D 3.59 GB, Cg 3.61 GB from `clingen-public`)
  took 96–132 s here, with all 18 jobs starting at once.

## Caveats

- This record was written after the runs; `submit_arms.sh` ran from a session
  scratch directory and is copied here unchanged.
- Outputs were compared by per-class counts only. The all-genomic evaluation,
  whose input contains every genomic expression selected here, found all four
  backends' outputs identical record for record.
- Earlier runs of this input on c4-standard-8 with `genomic-2026-10-27` (8-worker A
  184.8 s mean, with more spread between rounds) differ in machine type, store
  version and archive bucket together, so they are not directly comparable.
