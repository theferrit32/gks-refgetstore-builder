# vrsify over every genomic ClinVar expression: RefgetStore vs SeqRepo

clinvar-ingest's `vrsify` step translated every genomic expression in ClinVar
release 2026-10-24 (29,088,406 expressions of 4,681,155 variations) to GA4GH VRS
objects four times, with four sequence backends that hold identical sequences and
aliases. It ran on Cloud Batch, three rounds per backend, all on one machine type.

**Results**

- All 12 outputs are identical: the same 29,088,406 (key, VRS object or error)
  pairs from every backend in every round.
- Translation time with 8 workers on c4-highmem-8, mean of 3 rounds: RefgetStore
  archive 58.3 min; RefgetStore archive with range reads 60.6 min (1.04×, at 9% of
  the memory); remote RefgetStore with per-worker caches 68.4 min (1.17×);
  SeqRepo 81.8 min (1.40×).
- Of 10,646,293 (variation, sequence) groups of Alleles that carry two or more
  expressions, 10,646,285 agree on one VRS id; each of the 8 that do not holds
  expressions describing different changes in ClinVar's input. 72 of 1,105
  copy-number groups disagree; those inspected differ in their input coordinates.
- 99.35% of expressions translate. ClinVar's per-expression `issue` labels
  overstate failures: 98.8% of the expressions labelled "sequence for accession
  not supported by vrs-python release" translate.

Times in this record are UTC.

## Input: every genomic expression of release 2026-10-24

`variation_identity` holds one expression per variation, chosen by precedence.
The input here is every genomic candidate that choice is made from, built by a
one-off procedure modelled on clinvar-gkm's
`src/procedures/variation-identity-proc.sql` (at `e0fe8e9`, the definition
deployed in `clingen-dev.clinvar_ingest`):
[`artifacts/variation-all-genomic-exprs-proc.sql`](artifacts/variation-all-genomic-exprs-proc.sql).
It reads only tables the release already holds (`variation_identity`,
`variation_hgvs`, `variation_loc`) and keeps every candidate of precedence 1–5:

| precedence | fmt | candidates |
| --- | --- | --- |
| 1 | spdi | canonical SPDI |
| 2 | hgvs | top-level genomic HGVS, every assembly, every spelling ClinVar gives on the accession |
| 3 | gnomad | gnomAD-style id from the location, every assembly |
| 4 | hgvs | HGVS built from the location where no gnomAD id exists |
| 5 | hgvs | other genomic HGVS: NG_, LRG, alternate loci, patches |

An expression string two precedences produce is kept once, at the higher one.
Rows are limited to the classes vrsify translates (Allele, CopyNumberChange,
CopyNumberCount). Columns follow `variation_identity`'s, plus `hgvs_type`,
`clinvar_expression` (before the `del<N>` cleanup) and `selected` (the row
`variation_identity` chose: 4,681,155 rows, one per variation with a genomic
selection).

The procedure was not created in BigQuery. It ran once as a script, its
parameters turned into `DECLARE`s by
[`artifacts/proc_as_script.py`](artifacts/proc_as_script.py):

    python3 proc_as_script.py variation-all-genomic-exprs-proc.sql \
        "on_date=DATE '2026-10-24'" \
        "export_uri='gs://clinvar-gkm/2026-10-24/dev/scratch/all-genomic/vi-all-genomic-exprs-*.jsonl.gz'" \
        > variation-all-genomic-exprs-2026-10-24.sql
    bq --project_id=clingen-dev query --use_legacy_sql=false < variation-all-genomic-exprs-2026-10-24.sql

([`artifacts/variation-all-genomic-exprs-2026-10-24.sql`](artifacts/variation-all-genomic-exprs-2026-10-24.sql)
is the script that ran.) It processed 8.3 GB in about 17 s and wrote:

- table `clingen-dev.clinvar_2026_10_24_v2_6_1.variation_all_genomic_exprs`:
  29,088,406 rows
- export `gs://clinvar-gkm/2026-10-24/dev/scratch/all-genomic/vi-all-genomic-exprs-*.jsonl.gz`:
  119 objects, 506,963,347 bytes, newline-delimited JSON with integers as
  strings, as `bq extract` writes the `variation_identity` export

| precedence | fmt | GRCh38 | GRCh37 | NCBI36 | no assembly |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | spdi | 4,592,946 | 13 | | 1 |
| 2 | hgvs | 4,619,524 | 4,651,011 | 14,241 | |
| 3 | gnomad | 4,602,058 | 4,602,041 | 2,841 | |
| 4 | hgvs | 2,998 | 21,301 | 4,889 | |
| 5 | hgvs | | 2 | | 5,974,540 |

## Sequence backends

The backends hold the same 141,875 sequences and 173,732 aliases: the genomic
RefgetStore `genomic-2026-10-29` (built from `sources.genomic.toml`: GRCh38
and GRCh37 with their patches, NCBI36, T2T-CHM13v2.0, RefSeqGene with its
history, LRG, and RefSeq genomic biological regions) and a biocommons SeqRepo
written from it with `gks-refgetstore seqrepo`
([`../2026-10-29-genomic-seqrepo-update/`](../2026-10-29-genomic-seqrepo-update/README.md)).
Both vrs-python data proxies resolve all 319,331 identifiers to the same
digests. Publication:
[`../2026-10-29-genomic-2026-10-29-publish/`](../2026-10-29-genomic-2026-10-29-publish/README.md).

| arm | backend | what the job reads (all in us-east1) | how sequence is read |
| --- | --- | --- | --- |
| A | RefgetStore, local | `gs://clingen-public/refgetstore/genomic-2026-10-29.tar.gz` (3,587,446,072 bytes), extracted at startup | each worker loads a whole sequence on first use and keeps it (vrs-python's default) |
| B | RefgetStore, remote | `https://storage.googleapis.com/clingen-public/refgetstore/genomic-2026-10-29`, no archive | as A; each worker downloads the sequence files it uses into its own cache directory |
| D | RefgetStore, local | as A | range reads: `get_substring` reads only the requested bases from the `.seq` file; no sequence is kept in memory |
| Cg | SeqRepo | `gs://clingen-public/seqrepo/genomic-2026-10-29.tar.gz` (3,611,375,334 bytes), extracted at startup | SeqRepo's bgzip/faidx reads, with its file-descriptor cache on (`SEQREPO_FD_CACHE_MAXSIZE=none`) |

D is clinvar-ingest's `--refget-reads ranges`: in each worker it replaces
vrs-python's `RefgetSequenceDataProxy._ensure_sequence_loaded` with a no-op
(`clinvar_ingest/vrsify/runner.py`), so gtars serves `get_substring` from disk
instead of after `load_sequence`. It is not stock vrs-python behavior.

## Software and machines

- Image `gcr.io/clingen-dev/clinvar-ingest:vrsify-exp-12774e5`
  (`sha256:1cd99f1625626f6c31e2dce62f415027c8ef65ce0f385246b6bfe1113c28aa4d`),
  built from clinvar-ingest `kf/320-vrsify` at `12774e5`
- vrs-python `0.1.dev869+g8ce6d702e` (theferrit32/vrs-python `8ce6d702e3b9`),
  gtars `0.11.0` (theferrit32/gtars `7831e09f5960`), hgvs `2.0.0a0`,
  biocommons.seqrepo `0.6.11`, pysam `0.24.1`, as each provenance file records
- vrsify runs offline: UTA and network sequence fetching are disabled in the
  workers, so an accession missing from the backend is an error
- Cloud Batch, `us-east1`, one task per job on **c4-highmem-8** (8 vCPUs; task
  request 8,000 milli-CPU and 55,296 MiB), 100 GB `hyperdisk-balanced` boot disk,
  `maxRunDuration` 14,400 s, `--workers 8`, record timeout 300 s
- Each job, VM and disk is labelled `app=clinvar-ingest-vrsify`,
  `experiment=vrsify-store-comparison`, `arm`, `round`, `input=all-genomic`,
  `machine=c4-highmem-8`, `store=genomic-2026-10-29`, `trigger=manual`,
  `release-date=2026-10-24` (from clinvar-ingest `vrsify-submit --label`, not yet
  committed when these jobs were submitted)

## Commands

Jobs were submitted from the clinvar-ingest worktree with
[`artifacts/submit_arms.sh`](artifacts/submit_arms.sh), one invocation per round:

    MACHINE=c4-highmem-8 ARCHIVES=gcs WORKERS=8 LABELS="round=r1 input=all-genomic" \
        bash submit_arms.sh gcr.io/clingen-dev/clinvar-ingest:vrsify-exp-12774e5 \
        'gs://clinvar-gkm/2026-10-24/dev/scratch/all-genomic/vi-all-genomic-exprs-*.jsonl.gz' \
        gs://clinvar-gkm/2026-10-24/dev/scratch/allgen-hm-r1 A B D Cg

(`r2` and `r3` likewise). For each arm the script runs, for example for D:

    clinvar-ingest vrsify-submit 2026-10-24 --project clingen-dev --region us-east1 \
        --image gcr.io/clingen-dev/clinvar-ingest:vrsify-exp-12774e5 \
        --boot-disk-size 100 --machine-type c4-highmem-8 \
        --label experiment=vrsify-store-comparison --label arm=d \
        --label store=genomic-2026-10-29 --label workers=8 \
        --label machine=c4-highmem-8 --label round=r1 --label input=all-genomic -- \
        --input-uri 'gs://clinvar-gkm/2026-10-24/dev/scratch/all-genomic/vi-all-genomic-exprs-*.jsonl.gz' \
        --output-uri gs://clinvar-gkm/2026-10-24/dev/scratch/allgen-hm-r1-D/vi-vrs-python.jsonl.gz \
        --refget-store gs://clingen-public/refgetstore/genomic-2026-10-29.tar.gz \
        --refget-reads ranges --workers 8

A takes the same `--refget-store` without `--refget-reads`; B takes
`--refget-store https://storage.googleapis.com/clingen-public/refgetstore/genomic-2026-10-29`;
Cg takes `--seqrepo gs://clingen-public/seqrepo/genomic-2026-10-29.tar.gz`.

The 12 outputs were downloaded to
`/Volumes/RefgetStores/eval-outputs/allgen-hm.TM1y/<round>-<arm>.jsonl.gz`
(each the size of its GCS object) and analyzed with
[`artifacts/allgen_analysis.py`](artifacts/allgen_analysis.py):

    uv run --with duckdb python -I allgen_analysis.py \
        /Volumes/RefgetStores/eval-outputs/analysis.cDNg \
        /Volumes/RefgetStores/eval-outputs/allgen-hm.TM1y

Its output is [`artifacts/analysis.log`](artifacts/analysis.log); its DuckDB
database (tables `a`, every r1-A record with its key fields, and `groups`, the
per-(variation, sequence) VRS id counts) is
`/Volumes/RefgetStores/eval-outputs/analysis.cDNg/allgen.duckdb`.
[`artifacts/compare_arms.py`](artifacts/compare_arms.py) compares outputs record
by record on the same key; it was used on the smoke tests.

## Runs

[`artifacts/batch-jobs.txt`](artifacts/batch-jobs.txt) lists the jobs. Each
wrote `vi-vrs-python.jsonl.gz` and `vi-vrs-python.provenance.json` under
`gs://clinvar-gkm/2026-10-24/dev/scratch/allgen-hm-<round>-<arm>/`; the provenance
files are copied here as `artifacts/provenance-<round>-<arm>.json`.

| round | arm | Batch job | running → succeeded (UTC) | job run time |
| --- | --- | --- | --- | ---: |
| r1 | A | `vrsify-2026-10-24-20261029210252` | 21:03:43 → 22:06:08 | 62.4 min |
| r1 | B | `vrsify-2026-10-24-20261029210257` | 21:03:54 → 22:15:02 | 71.1 min |
| r1 | D | `vrsify-2026-10-24-20261029210302` | 21:03:54 → 22:07:09 | 63.2 min |
| r1 | Cg | `vrsify-2026-10-24-20261029210311` | 21:04:11 → 22:28:42 | 84.5 min |
| r2 | A | `vrsify-2026-10-24-20261029213221` | 21:33:09 → 22:33:49 | 60.7 min |
| r2 | B | `vrsify-2026-10-24-20261029213225` | 21:33:18 → 22:42:32 | 69.2 min |
| r2 | D | `vrsify-2026-10-24-20261029213228` | 21:33:17 → 22:36:34 | 63.3 min |
| r2 | Cg | `vrsify-2026-10-24-20261029213232` | 21:33:19 → 22:58:05 | 84.8 min |
| r3 | A | `vrsify-2026-10-24-20261029213235` | 21:33:30 → 22:33:20 | 59.8 min |
| r3 | B | `vrsify-2026-10-24-20261029213238` | 21:33:28 → 22:42:58 | 69.5 min |
| r3 | D | `vrsify-2026-10-24-20261029213242` | 21:33:37 → 22:36:55 | 63.3 min |
| r3 | Cg | `vrsify-2026-10-24-20261029213245` | 21:33:32 → 22:56:20 | 82.8 min |

Rounds 2 and 3 ran at the same time on separate VMs; round 1 started 30 minutes
earlier.

### Runs that preceded these

- **Smoke tests, one input shard** (`vi-all-genomic-exprs-000000000000.jsonl.gz`,
  244,991 rows), 8 workers on c4-standard-8, all four arms. Against the
  `genomic-2026-10-27` artifacts (jobs `…20261029182631` to `…182642`, outputs
  `scratch/allgen-smoke-<arm>/`): identical outputs, 4,717 errors. Against
  `genomic-2026-10-29` (jobs `…20261029203150` to `…203201`, outputs
  `scratch/allgen-smoke1029-<arm>/`): identical outputs, 1,622 errors; 3,095 NG_
  records that failed against 2026-10-27 translate. Whole-sequence workers
  peaked at up to 2,984 MiB each, 22.1 GiB for 8.
- **Full input on c4-standard-8** (jobs `…20261029203955` A, `…203959` B,
  `…204002` D, `…204005` Cg): A and B failed after 5,657,950 and 5,380,800
  records with `BrokenProcessPool` (a worker killed abruptly), their 8
  whole-sequence workers having outgrown the 27,648 MiB task memory; D and Cg
  were cancelled so that every arm would run on one machine type. Batch does not
  forward kernel messages to Cloud Logging, so the out-of-memory kill is inferred
  from the worker death and the c4-highmem-8 peaks below (31.8–32.1 GiB).

## Results

### Outputs

Every output holds 29,088,406 records with distinct keys, where a key is
(`variation_id`, `fmt`, `assembly_version`, `source`). An order-independent
fingerprint of (key, output) — two 64-bit XOR-folded hashes — is
`3c5e751223f99130/2eea0c7272a6e8e8` for all 12, so no backend or round changed
any output. Per output:

| class | translated | error |
| --- | ---: | ---: |
| Allele | 28,782,039 | 185,835 |
| CopyNumberChange | 41,697 | 3,259 |
| CopyNumberCount | 75,419 | 157 |
| total | 28,899,155 | **189,251 (0.65%)** |

### Time and resources

Per run, from the provenance files. "Setup" is fetching and extracting the
backend; "translate" is the translation itself. Memory is each worker's peak
RSS; context switches and CPU are summed or averaged over the 8 workers.

| arm | round | setup s | translate s | worker peak RSS max (MiB) | sum of peaks (GiB) | downloaded during run (GB) | CPU busy | voluntary ctx sw | involuntary ctx sw | slowest record s |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A | r1 | 86.8 | 3,584.9 | 4,485 | 31.8 | — | 0.981 | 166,605 | 4,287,245 | 46.5 |
| A | r2 | 86.7 | 3,478.6 | 4,485 | 31.9 | — | 0.981 | 174,673 | 4,262,043 | 45.0 |
| A | r3 | 81.9 | 3,427.7 | 4,485 | 31.8 | — | 0.981 | 173,578 | 4,308,644 | 45.2 |
| B | r1 | 0.1 | 4,174.6 | 4,502 | 32.1 | 31.5 | 0.846 | 1,772,693 | 3,983,584 | 45.8 |
| B | r2 | 0.1 | 4,062.8 | 4,501 | 31.9 | 31.4 | 0.857 | 1,700,975 | 3,987,730 | 46.1 |
| B | r3 | 0.1 | 4,078.4 | 4,504 | 32.0 | 31.5 | 0.858 | 1,863,224 | 3,957,186 | 64.2 |
| D | r1 | 85.4 | 3,630.4 | 748 | 2.9 | — | 0.981 | 157,962 | 4,501,082 | 63.0 |
| D | r2 | 86.3 | 3,635.4 | 740 | 2.9 | — | 0.981 | 152,427 | 4,517,879 | 62.9 |
| D | r3 | 85.8 | 3,636.5 | 744 | 2.9 | — | 0.981 | 154,906 | 4,455,224 | 63.4 |
| Cg | r1 | 54.4 | 4,939.4 | 778 | 2.5 | — | 0.985 | 48,638 | 4,593,277 | 140.1 |
| Cg | r2 | 54.2 | 4,956.1 | 777 | 2.5 | — | 0.985 | 49,061 | 4,585,789 | 140.9 |
| Cg | r3 | 55.9 | 4,833.6 | 777 | 2.5 | — | 0.985 | 49,068 | 4,590,196 | 138.5 |

Means over the three rounds:

| arm | translate | sd | vs A | setup | setup + translate | throughput |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A | 3,497.1 s (58.3 min) | 80.2 s | 1.000× | 85.1 s | 59.7 min | 8,318 records/s |
| D | 3,634.1 s (60.6 min) | 3.3 s | 1.039× | 85.8 s | 62.0 min | 8,004 records/s |
| B | 4,105.3 s (68.4 min) | 60.5 s | 1.174× | 0.1 s | 68.4 min | 7,086 records/s |
| Cg | 4,909.7 s (81.8 min) | 66.4 s | 1.404× | 54.8 s | 82.7 min | 5,925 records/s |

- **D** costs 4% over A in time and uses 9% of its memory (2.9 against 31.8 GiB
  for 8 workers). Its times are also the most repeatable.
- **B** spends its extra 17% waiting for sequence downloads: CPU busy falls to
  85%, and voluntary context switches rise tenfold. Each of its 8 workers fills
  its own cache, so it downloads 31.5 GB in all, about 3.9 GB per worker, against
  the 3.59 GB archive A and D fetch once.
- **Cg** (SeqRepo) is 40% slower than A on the same data, at D's memory.
- The slowest records in every run are deletions and duplications inside the run
  of N at the start of GRCh37 chr15 (`NC_000015.9:g.5185del`,
  `g.2897_2953del`, `g.7211_7214del`, `g.2831_2832dup`). VOCA normalization
  rolls over the N run one base at a time, about 20 million single-base reads per
  record, so these records show each backend's per-read cost: about 45 s on A
  and B, 63 s on D, 140 s on SeqRepo.

### One VRS id per variation and sequence

For each (variation, sequence) — the sequence being the translated object's
`location.sequenceReference.refgetAccession`, so aliases of one sequence count
once — do all of a variation's translated expressions give the same VRS id?
From r1-A, identical to every other output:

| class | groups | with 2 or more expressions | agree | disagree |
| --- | ---: | ---: | ---: | ---: |
| Allele | 13,584,997 | 10,646,293 (25,843,335 records) | 10,646,285 | 8 |
| CopyNumberChange | 40,598 | 1,092 | 1,026 | 66 |
| CopyNumberCount | 75,406 | 13 | 7 | 6 |

The disagreements come from ClinVar's input, not from translation:

- All 8 Allele groups hold expressions that describe different changes.
  Variation 1164852 has `NC_000016.9:g.138773dup` (one base duplicated) and
  gnomAD `16-138772-A-AGG` (two bases inserted); 1328344 has
  `NC_000019.9:g.40392316G>A` and `g.40408254G>A`; 37122 has
  `NG_008806.1:g.75301C>T` and `g.75301T>C`; 95104 has `NG_009382.2:g.32608G>A`
  and `g.32608G>T`; 167734 has `NG_032976.1:g.10901=` beside `g.10901dup`;
  1510451, 989413 and 995996 have NG_/LRG expressions at two positions. Where
  the same change is spelled on both an NG_ and its LRG, both give one id.
- The copy-number groups inspected have ClinVar HGVS and location-derived HGVS
  (precedence 4) that differ by one or two bases (variation 1172761:
  `NC_000002.12:g.26398040_26425787del` against `g.26398038_26425785del`), or
  uncertain ranges that differ in an outer bound (100778, 32067). Copy-number
  locations are not normalized, so different coordinates give different ids.

### ClinVar's `issue` against the outcome

`variation_identity_build` labels expressions it expects vrs-python to reject.
vrsify ignores the label and attempts every expression. The labels here are
per expression: `variation_all_genomic_exprs` recomputes the HGVS rules for
every spelling and does not carry the per-variation issues `variation_identity`
holds (such as `range copies are not supported.`).

| issue | records | translated | error |
| --- | ---: | ---: | ---: |
| (none) | 27,444,311 | 27,407,645 | 36,666 |
| sequence for accession not supported by vrs-python release | 1,508,628 | 1,490,879 | 17,749 |
| repeat expressions are not supported. | 129,669 | 0 | 129,669 |
| unsupported hgvs expression. | 5,694 | 631 | 5,063 |
| expression contains unbalaned paretheses. | 97 | 0 | 97 |
| intronic positions are not resolvable in sequence. | 7 | 0 | 7 |

The accession rule predates LRG, GenBank and the other namespaces in this store;
the repeat rule matches vrs-python.

Errors on expressions without an issue:

| error | records | cause |
| --- | ---: | --- |
| `UnresolvableAccessionError: Unable to resolve <accession>` | 29,968 | accession in no source: chiefly NG_ versions absent from both RefSeqGene and the biological-region file |
| `UnresolvableAccessionError: Unable to resolve GRCh37:MT` | 3,115 | GRCh37 gnomAD-style mitochondrial Allele ids; the store's `GRCh37` namespace mirrors the original GRCh37 report, which has no mitochondrion (it entered with GRCh37.p2). clinvar-ingest#326 |
| `ValueError: Unable to parse data as hgvs variation` | 3,053 | gnomAD-style expressions of copy-number variations: vrsify sends every copy-number record to vrs-python's `CnvTranslator` as HGVS, the only format it supports. All 3,053 gnomAD rows of CopyNumberChange and CopyNumberCount, including 21 large mtDNA deletions. clinvar-ingest#328 |
| `DataProxyValidationError: Reference mismatch` | 315 | the expression's reference base is not the sequence's: 266 on NCBI36, 49 on GRCh37/38 |
| `KeyError: 'absolute_copies'` | 127 | CopyNumberCount records with only `range_copies` (95 variations, all labelled `range copies are not supported.` in `variation_identity`); vrsify's `translate.py` reads `absolute_copies` unconditionally, though vrs-python accepts a Range for `copies`. clinvar-ingest#327 |
| `ValueError: Must provide a 'del' or 'dup'` | 53 | copy-number HGVS that is neither |
| `ValidationError` on `Range`, uncertain ranges, interval bounds | 23 | malformed or unsupported ranges |
| `not translated: non-genomic expression` | 9 | `c.` expressions in rows typed genomic |
| other | 3 | |

## Data

| what | where |
| --- | --- |
| input table | `clingen-dev.clinvar_2026_10_24_v2_6_1.variation_all_genomic_exprs` |
| input export | `gs://clinvar-gkm/2026-10-24/dev/scratch/all-genomic/vi-all-genomic-exprs-*.jsonl.gz` |
| outputs | `gs://clinvar-gkm/2026-10-24/dev/scratch/allgen-hm-<r1\|r2\|r3>-<A\|B\|D\|Cg>/vi-vrs-python.jsonl.gz` (each about 2.88 GB) |
| provenance | the same prefixes, `vi-vrs-python.provenance.json`; copies in [`artifacts/`](artifacts/) |
| smoke-test outputs | `gs://clinvar-gkm/2026-10-24/dev/scratch/allgen-smoke-<arm>/`, `…/allgen-smoke1029-<arm>/` |
| RefgetStore | `gs://clingen-public/refgetstore/genomic-2026-10-29/` and `.tar.gz`; also `https://static.ferriter.dev/refgetstore/genomic-2026-10-29/` and `.tar.gz` |
| SeqRepo | `gs://clingen-public/seqrepo/genomic-2026-10-29.tar.gz`, with its `.sha256` and `genomic-2026-10-29.provenance.json` |
| local copies of the outputs | `/Volumes/RefgetStores/eval-outputs/allgen-hm.TM1y/` (32.2 GiB) |
| analysis database | `/Volumes/RefgetStores/eval-outputs/analysis.cDNg/allgen.duckdb` (10 GB) |
| local stores | `/Volumes/RefgetStores/genomic-2026-10-29`, `/Volumes/RefgetStores/seqrepo/genomic-2026-10-29` |

## Caveats

- This record was written after the runs. The commands above are the ones
  submitted; the submission script, BigQuery script and analysis script ran from
  a session scratch directory and are copied here unchanged.
- vrs-python and gtars are forks pinned by commit (see Software); D replaces a
  private vrs-python method at run time.
- B's cost includes duplicated downloads: each worker caches independently,
  because gtars writes cache files non-atomically, so workers cannot share one
  cache directory safely.
- Single-worker timings were measured on the one-expression-per-variation input
  only, not on this input.
- ~1.4% of the 2026-10-24 all-genomic expressions are on NG_ records that are not
  RefSeqGene; `genomic-2026-10-29` added NCBI's `human.biological_region.fna.gz`,
  which supplies 92% of them. The rest remain the main error source above.
