# Add NCBI36 and T2T-CHM13v2.0

The store gained the NCBI36 (hg18) and T2T-CHM13v2.0 assemblies: two
collections, 127 new sequences and their assembly aliases, added by selective
builds from the manifest committed as `9fa66d1`.

| | before | after |
| --- | ---: | ---: |
| sequences | 1,781,533 | 1,781,660 |
| collections | 236 | 238 |
| locked input files | 336 | 340 |

## NCBI36

NCBI publishes the NCBI36 assembly report (`GCF_000001405.12`) but no genomic
FASTA, and no GenBank NCBI36 assembly directory exists. NCBI36 is therefore two
manifest entries:

- `ncbi36_refseq`, an `efetch_accessions` seqset over
  `accession_lists/ncbi36_refseq.txt`: the report's 123 RefSeq accessions (24
  `NC_` chromosomes, 98 `NT_` scaffolds, 1 `NG_` alternate scaffold). `fetch`
  retrieved all 123 in one efetch request; the cached FASTA is 832,545,989
  bytes, SHA-256 `d2612de0…`. `build --seqset ncbi36_refseq` made collection
  `GFxy8xYVWbzxOKJcWiq97r3nGdTQ33yE`: 103 sequences new to the store, 20 already
  present from other sources.
- `NCBI36`, a report-only assembly. `build --assembly NCBI36` resolved all 123
  report rows against those sequences and added 149 aliases.

## T2T-CHM13v2.0

The RefSeq assembly `GCF_009914755.1`. `build --assembly T2T-CHM13v2.0` made
collection `eM0RH5r4wQtxMzLTl1NUkkehbww-cWvM` with the 24 nuclear chromosomes
(`NC_060925.1`–`NC_060948.1`) and added 120 sequence and 2 collection aliases.
The report's 25th row, the mitochondrion, is GenBank-only (`CP068254.1`,
RefSeq `na`): it is absent from the RefSeq FASTA and so from the store, and its
row was skipped. `CP068254.1` is not rCRS (`NC_012920.1`): it starts at rCRS
position 577 and differs at three sites.

## Validation

- `verify`: pass. 340 locked files, 1,781,660 sequences and 238 collections
  reachable, every payload on disk.
- `tests`: 272 passed. The first run failed one test that assumed a single
  efetch seqset; it is now parameterized over both checked-in lists.
- `vrs-python-check` ([`artifacts/check_new_sources.py`](artifacts/check_new_sources.py)):
  for each of the 23 accessions on which the 2026-10-27 ClinVar genomic
  evaluation succeeded only against SeqRepo (22 NCBI36 chromosomes and
  `NC_060946.1`), one such expression
  ([`artifacts/seqrepo_only_samples.json`](artifacts/seqrepo_only_samples.json))
  translated through vrs-python's RefgetStore backend with the store
  sandboxed read-only ([`artifacts/live-readonly.sb`](artifacts/live-readonly.sb)).
  23/23 succeeded with SeqRepo's VRS ID.

## Notes

`sync` was not used: it reprocesses seqsets only, and assemblies need a
selective `build`. `--force-lock` was required because the new files had no
lock baseline.

The record was initialized from a dirty tree; every command ran against the
`sources.toml` committed as `9fa66d1`, whose SHA-256 the lock records.

The published store, `refgetstore/2026-10-07`, has none of these sources.
