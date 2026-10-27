# Upload the genomic-only store to refgetstore/genomic-2026-10-27

`theferrit32-public:theferrit32-public/refgetstore/genomic-2026-10-27` holds the
genomic-only store built from [`sources.genomic.toml`](../../sources.genomic.toml):
every assembly (GRCh38 and GRCh37 with their patches, NCBI36, T2T-CHM13v2.0),
NCBI36's efetched sequences, RefSeqGene with its history, and LRG.

Public URL: <https://static.ferriter.dev/refgetstore/genomic-2026-10-27/>

| | genomic store | full store (`store/`) |
| --- | ---: | ---: |
| sequences | 13,615 | 1,781,660 |
| collections | 33 | 238 |
| bytes | 4,886,047,853 | ~24 GB |
| sequence index + aliases | ~4.4 MB | ~1.36 GB |

## The store

It was built from the full build's lock, not by discovery, so it is exactly a
subset of the full store and used the same download cache:

    gks-refgetstore build --config sources.genomic.toml \
        --store-dir <dir> --lock <copy of build.lock.json> \
        --locked-sources --lock-check-mode subset --no-lock

The build re-hashed every cached source against the lock and ingested only the
lock's files owned by the manifest's 29 entries. It exited 0 in 237 s. Its 5
warnings name `NW_003871056.2`, `NW_004166864.1` and `NW_004504303.1`, which
the GRCh37.p11 and p12 reports reference and which no source supplies; the full
store lacks them too.

[`artifacts/check_subset_store.py`](artifacts/check_subset_store.py)
([output](artifacts/subset-check.txt)) passes: every collection the full lock
attributes to those 29 sources is present with the lock's sequence count and no
other is, every indexed sequence has its payload and no payload is unindexed,
and every sequence and collection file is byte-identical to the full store's.

## Upload

`upload_store.sh --remote theferrit32-public --bucket theferrit32-public
--prefix refgetstore/genomic-2026-10-27 --store-dir
/Volumes/RefgetStores/genomic-2026-10-27 --config sources.genomic.toml
--transfers 64 --checkers 64`, with `RCLONE_CHECKSUM=true`, from commit
`b70c3b4`. `--config` makes `manifest.json` list the genomic sources rather than
`sources.toml`'s.

- The prefix did not exist. A dry run first predicted 13,682 copies (the
  store's 13,681 files plus `manifest.json`) and no deletions.
- 4.55 GiB in 116 s.
- `rclone check`: **0 differences, 13,682 matching files.**
- [`artifacts/published-manifest.json`](artifacts/published-manifest.json)
  records `"dirty": true`. The tree was clean at `b70c3b4`, as this record's
  `git` block shows; the flag counts the untracked directory of this record
  itself, which `run_record init` had just created.

## Reading it back

Through vrs-python's `refget+https` backend, the first check
([`logs/public-read-check.log`](logs/public-read-check.log)) resolved none of 16
expressions. That is not the upload: gtars' `open_remote` defers loading the
sequence index, `get_sequence_metadata_by_alias` returns `None` until it is
loaded, and the Python binding has no call that loads it. Any remote store
read through vrs-python is affected.

Calling `load_sequence` on a digest that does not exist loads the index before
it raises. With that ([`artifacts/check_remote_sources.py`](artifacts/check_remote_sources.py)),
all 16 ([`artifacts/remote_samples.json`](artifacts/remote_samples.json): NG_,
LRG, CM, NW_, NT_, chrM, GRCh38 chr22, NCBI36 chr21, T2T chr22) translate
through the public URL with the full store's VRS IDs
([log](logs/public-read-check-index-loaded.log)).

A direct probe also showed `get_substring` without `load_sequence` reading 10
bases of GRCh38 chr22 by HTTP range: the cache held 4.3 MB of metadata and no
`.seq` file. vrs-python calls `load_sequence` before every read, so it downloads
whole payloads instead; this check's cache reached 138 MB.
