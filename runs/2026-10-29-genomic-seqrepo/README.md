# A SeqRepo instance mirroring the genomic store

`/Volumes/RefgetStores/seqrepo/genomic-2026-10-27` is a biocommons SeqRepo
instance holding exactly the sequences and sequence aliases of the genomic-only
RefgetStore `/Volumes/RefgetStores/genomic-2026-10-27`
([`../2026-10-27-genomic-r2-upload/`](../2026-10-27-genomic-r2-upload/README.md)).
It is the SeqRepo side of a RefgetStore-versus-SeqRepo comparison on identical
data.

| | SeqRepo instance | RefgetStore |
| --- | ---: | ---: |
| sequences | 13,615 | 13,615 |
| residues | 13,787,182,066 | 13,787,182,066 |
| sequence aliases (accession and assembly namespaces) | 45,471 | 45,471 |
| digest aliases (`VMC`, `MD5`, `SHA1`, `SEGUID`) | 54,460 | — |
| files | 41 | 13,681 |
| bytes | 3,528,109,422 | 4,886,047,853 |

## How it was built

    gks-refgetstore seqrepo --config sources.genomic.toml --lock build.lock.json \
        --lock-check-mode subset --store-dir /Volumes/RefgetStores/genomic-2026-10-27 \
        --target-dir /Volumes/RefgetStores/seqrepo --instance-name genomic-2026-10-27

with `biocommons.seqrepo` 0.6.11, `bioutils` 0.6.2, `pysam` 0.24.1, htslib
`bgzip` 1.24, and gtars from `theferrit32/gtars@7831e09`. The run took 1,841 s
to store every sequence and 165 s to read them all back.

- **Store and lock.** `sources.genomic.toml` selects 33 collections from the
  full `build.lock.json`, the same subset its store build used. The store holds
  exactly those collections, with the locked sizes, and no other sequence.
- **Writes.** Each store sequence was passed once to the seqrepo library's
  `SeqRepo.store(seq, aliases)` with every alias the store holds for it, then
  `commit()`. Seqrepo computed every `seq_id`, wrote 13 bgzip FASTA files of
  about 1 Gbp each with their `.fai`/`.gzi` indexes, inserted the alias rows,
  and added its four digest aliases per sequence. Every `seq_id` equals the
  store's sha512t24u digest.
- **Namespaces.** In `aliases.sqlite3`, the store's `refseq` is `NCBI` (10,437)
  and `lrg` is `LRG` (5,905); seqrepo answers lookups under either spelling.
  `insdc` (996) and the 25 assembly namespaces (`GRCh38`, `GRCh38.p1`–`.p14`,
  `GRCh37`, its patches, `NCBI36`, `T2T-CHM13v2.0`; 28,133 aliases) are stored
  under the store's names.
- **Assembly names.** `seqrepo add-assembly-names` was not run. The store's
  assembly namespaces already hold every sequence and UCSC name that bioutils'
  assembly data lists for their RefSeq accessions, with the same sequence; that
  command would also add 29 `NCBI34`/`NCBI35` names the store does not have.

[`artifacts/genomic-2026-10-27.provenance.json`](artifacts/genomic-2026-10-27.provenance.json)
records the manifest, lock and store digests, the tool versions, the counts,
the verification result, and the size and SHA-256 of every file in the
instance. The same file sits beside the instance.

## Validation

- **Built-in verification** (`build`): every sequence read back from SeqRepo
  re-hashes to its `seq_id` and matches the store's length and MD5; the two
  sequence sets are equal; the 45,471 aliases are equal in both directions; each
  sequence has one alias in each digest namespace; no alias row is non-current.
- **Data proxies** (`proxy-check`, [`tools/check_seqrepo_proxies.py`](../../tools/check_seqrepo_proxies.py)):
  through vrs-python's `seqrepo+file://` and `refget+file://` data proxies,
  62,809 identifiers resolve to the same `SQ.` digest, the one the store files
  them under: every namespaced alias plus 17,338 bare `refseq`, `lrg` and
  `insdc` accessions. All 13,615 lengths and the first and last 25 residues of
  545 sequences agree. 0 differences
  ([`artifacts/proxy-differences.tsv`](artifacts/proxy-differences.tsv) is
  header-only).
- **VRS IDs** (`samples-refget`, `samples-seqrepo-fdcache`): the 16 expressions
  in [`../2026-10-27-genomic-r2-upload/artifacts/remote_samples.json`](../2026-10-27-genomic-r2-upload/artifacts/remote_samples.json)
  (NG_, LRG, CM, NW_, NT_, chrM, GRCh38 chr22, NCBI36 chr21, T2T chr22)
  translate to the expected VRS IDs on both backends, through the eval
  harness's `WorkerRuntime` with UTA `uta_20241220`. This check is about VRS ID
  equality; its timings are not a performance comparison.

The SeqRepo sample check runs with `SEQREPO_FD_CACHE_MAXSIZE=none`, as
clinvar-ingest's vrsify does. Without it, seqrepo opens the bgzip file and its
indexes for every read. `samples-seqrepo`, run without it, was at full CPU for
24 minutes without finishing and was stopped (exit 143); native stack samples
showed it in htslib's file open and close. Its output was block-buffered, so the
sample it stalled on is not in its log. With the cache enabled the 16 samples
completed in 34 s, and the longest pause was at
`NC_000021.7:g.13636378_28138533dup`, a duplication spanning the NCBI36 chr21
p-arm gap, whose normalization reads many short ranges.

`tests`: the repository's 281 tests pass, including those for the subcommand.
