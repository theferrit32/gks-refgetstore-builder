# Add RefSeq genomic biological regions to the genomic SeqRepo

`/Volumes/RefgetStores/seqrepo/genomic-2026-10-29` is the genomic SeqRepo
instance extended to mirror the genomic store `genomic-2026-10-29`, which adds
NCBI's RefSeq genomic biological regions
([`../2026-10-29-biological-region/`](../2026-10-29-biological-region/README.md)).
It is an APFS clone of `genomic-2026-10-27`
([`../2026-10-29-genomic-seqrepo/`](../2026-10-29-genomic-seqrepo/README.md))
with the new sequences added by `gks-refgetstore seqrepo --add-missing`.

| | genomic-2026-10-27 | genomic-2026-10-29 |
| --- | ---: | ---: |
| sequences | 13,615 | 141,875 |
| residues | 13,787,182,066 | 13,875,659,588 |
| sequence aliases (accession and assembly namespaces) | 45,471 | 173,732 |
| `NCBI` (`refseq`) aliases | 10,437 | 138,698 |
| alias rows, including seqrepo's digest aliases | 99,931 | 741,232 |
| files | 41 | 74 |
| bytes | 3,528,109,422 | 3,751,758,315 |

## How it was extended

    cp -c -R seqrepo/genomic-2026-10-27 seqrepo/genomic-2026-10-29
    gks-refgetstore seqrepo --config sources.genomic.toml --lock build.lock.json \
        --lock-check-mode subset --store-dir /Volumes/RefgetStores/genomic-2026-10-29 \
        --target-dir /Volumes/RefgetStores/seqrepo --instance-name genomic-2026-10-29 \
        --add-missing

The store matched the 34 collections the lock pins for `sources.genomic.toml`.
`--add-missing` opened the clone writeable and passed to `SeqRepo.store()` the
128,260 store sequences the instance lacked, each with its `refseq` aliases
(128,261, two records sharing a sequence), then committed: 88,477,522
residues in 163 s. No existing sequence lacked an alias. Seqrepo wrote them to
11 new bgzip files and left the existing 13 unchanged.

## Verification

- **Built-in** (`add-missing-verify`): every one of the 141,875 sequences read
  back from SeqRepo re-hashes to its `seq_id` and matches the store's length and
  MD5; the sequence sets are equal; the 173,732 aliases are equal in both
  directions; each sequence has one alias in each digest namespace. The
  sequence pass took 35 s.
- **Data proxies** (`proxy-check`): 319,331 identifiers (every alias plus
  145,599 bare `refseq`, `lrg` and `insdc` accessions) resolve to the same `SQ.`
  digest through vrs-python's `seqrepo+file://` and `refget+file://` proxies;
  all 141,875 lengths and the ends of 5,675 sequences agree. 0 differences.
- **Formerly missing accessions** (`accessions`,
  [`artifacts/check_accessions.py`](artifacts/check_accessions.py)):
  `NG_042865.1`, `NG_046166.1`, `NG_046616.1` and `NG_050744.1` resolve in
  neither 2026-10-27 backend and to the same digest in both 2026-10-29
  backends.
- **VRS IDs** (`samples-seqrepo`, `samples-refget`): the 16 expressions in
  [`../2026-10-27-genomic-r2-upload/artifacts/remote_samples.json`](../2026-10-27-genomic-r2-upload/artifacts/remote_samples.json)
  translate to the expected VRS IDs on both new backends. This check is about
  VRS ID equality, not timing; the SeqRepo runs use
  `SEQREPO_FD_CACHE_MAXSIZE=none`, as clinvar-ingest's vrsify does.

The first run, `add-missing` (exit 143), committed the additions and was
stopped during verification. Its MD5 comparison looked up each sequence's MD5
alias with a query that SQLite answers by scanning every MD5 row, quadratic in
the sequence count. Verification now reads the MD5 aliases in one query and
keeps bgzip files open across reads. `add-missing-verify` reran the same
command: it found nothing to add, verified the instance, and wrote the
provenance record, so that record's `write` block counts zero sequences. The
counts above come from the first run's log.

[`artifacts/genomic-2026-10-29.provenance.json`](artifacts/genomic-2026-10-29.provenance.json)
is the provenance record, also beside the instance. `genomic-2026-10-27` is
unchanged: every file in it hashes as before.
