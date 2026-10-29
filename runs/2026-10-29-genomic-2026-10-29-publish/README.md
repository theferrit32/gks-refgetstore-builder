# Publish genomic-2026-10-29

`genomic-2026-10-29` is `genomic-2026-10-27` plus the RefSeq genomic biological
regions ([`../2026-10-29-biological-region/`](../2026-10-29-biological-region/README.md)),
and `seqrepo/genomic-2026-10-29` is the SeqRepo holding the same sequences and
aliases ([`../2026-10-29-genomic-seqrepo-update/`](../2026-10-29-genomic-seqrepo-update/README.md)).
This record publishes both under new prefixes; every `2026-10-27` object is
unchanged.

| object | where |
| --- | --- |
| store, per file | <https://static.ferriter.dev/refgetstore/genomic-2026-10-29/> (R2) and <https://storage.googleapis.com/clingen-public/refgetstore/genomic-2026-10-29/> (GCS, us-east1) |
| store archive | `refgetstore/genomic-2026-10-29.tar.gz` on both, 3,587,446,072 bytes, SHA-256 `da845a982251f368acd23f221ed8a5b4c173684da1797eb2bc21f0631991b6cd` |
| SeqRepo archive | `gs://clingen-public/seqrepo/genomic-2026-10-29.tar.gz`, 3,611,375,334 bytes, SHA-256 `181e9b794515302364fd0aab9df51fbfe16d802617b7b80d67682cc812e93c9c`, with the build's `genomic-2026-10-29.provenance.json` beside it |

Each archive has a `.sha256` beside it and holds the single directory
`genomic-2026-10-29/`.

## How it was published

- **Store to R2:** `RCLONE_CHECKSUM=true ./upload_store.sh --remote
  theferrit32-public --bucket theferrit32-public --prefix
  refgetstore/genomic-2026-10-29 --store-dir
  /Volumes/RefgetStores/genomic-2026-10-29 --config sources.genomic.toml --copy
  --transfers 64 --checkers 64`. The store directory, cloned from
  `genomic-2026-10-27`, still held that store's `manifest.json`; the script
  regenerated it for this store and destination.
- **Store to GCS:** the same script with an environment-only GCS remote
  (`RCLONE_CONFIG_CLINGENGCS_*`, gcloud application-default credentials) and
  `--no-manifest`, so the mirror carries the R2 `manifest.json`, as the
  `genomic-2026-10-27` mirror does.
- **Archives:** made with the same flags as the earlier archives:

      COPYFILE_DISABLE=1 tar --no-mac-metadata --no-xattrs --no-acls \
          --uid 0 --gid 0 --uname root --gname root -C <parent> -cf - genomic-2026-10-29 | pigz -6

  then `shasum -a 256`, and `rclone copy --checksum` to each destination.

## Validation

- `rclone check` after each upload: 0 differences for the store on R2 and on
  GCS (141,943 files each) and for every archive and its `.sha256`.
- Both archives extract to trees identical to their source directories.
- Every published archive serves its full length over unauthenticated HTTPS,
  and each served `.sha256` holds its archive's digest.
- From both store URLs, [`artifacts/gcs_readback.py`](artifacts/gcs_readback.py)
  translates the 16 sample expressions through clinvar-ingest vrsify's worker
  to the expected VRS IDs, and
  [`artifacts/remote_ng_check.py`](artifacts/remote_ng_check.py) resolves
  NG_042865.1, NG_046166.1 and NG_046616.1, absent from `genomic-2026-10-27`,
  to the same digests with NCBI's lengths.

## Missing manifest entries

Nine commands ran under `run_record exec` while `upload-r2` was still running,
and `upload-r2`'s exit rewrote the manifest from the copy it read when it
started, dropping their entries: the SeqRepo archive's create, verify, upload
and check; the GCS mirror upload; and the store archive's create, verify, GCS
upload and GCS check. Their logs in [`logs/`](logs/) are intact, and the
commands are the ones above. `run_record exec` does not lock the manifest, so
concurrent execs into one record are unsafe.
