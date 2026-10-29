# Mirror the genomic-only store to clingen-public

`gs://clingen-public/refgetstore/genomic-2026-10-27` is a copy of the store
published on R2 at `refgetstore/genomic-2026-10-27`
([`../2026-10-27-genomic-r2-upload/`](../2026-10-27-genomic-r2-upload/README.md)),
in the same us-east1 GCS region as the archives and the Cloud Batch jobs that
read it. `clingen-public` grants `allUsers` object read, which gtars needs: it
reads a remote store over unauthenticated HTTPS.

Public URL: <https://storage.googleapis.com/clingen-public/refgetstore/genomic-2026-10-27/>

This record is reconstructed: the upload was not run under `run_record exec`.
Its rclone log is [`artifacts/upload.log`](artifacts/upload.log), and the
read-back check was rerun under this record.

## Upload

From `/Volumes/RefgetStores/genomic-2026-10-27` with `upload_store.sh` at
`d8d8a8a`, through a GCS remote defined only by environment variables
(gcloud application-default credentials):

    RCLONE_CONFIG_CLINGENGCS_TYPE="google cloud storage" \
    RCLONE_CONFIG_CLINGENGCS_ENV_AUTH=true \
    RCLONE_CONFIG_CLINGENGCS_BUCKET_POLICY_ONLY=true \
    RCLONE_CHECKSUM=true RCLONE_STATS=30s \
    ./upload_store.sh --remote clingengcs --bucket clingen-public \
        --prefix refgetstore/genomic-2026-10-27 \
        --store-dir /Volumes/RefgetStores/genomic-2026-10-27 \
        --no-manifest --copy --transfers 64 --checkers 64

- The prefix did not exist, and a dry run first planned 4.550 GiB.
- `--no-manifest` uploads the store directory's existing `manifest.json`, which
  is byte-identical to R2's, so the mirror matches R2 file for file. Its
  `destination` field still names the R2 bucket.
- `--copy` is additive: nothing at the destination could be deleted.
- 4.550 GiB in about 4.5 minutes.

## Validation

- `rclone check`: **0 differences, 13,682 matching files.**
- `rgstore.json` serves over unauthenticated HTTPS (HTTP 200) with the local
  file's MD5.
- [`artifacts/gcs_readback.py`](artifacts/gcs_readback.py) opens the public URL
  through clinvar-ingest vrsify's own worker setup (offline, per-worker cache)
  and translates the 16 sample expressions from
  [`../2026-10-27-genomic-r2-upload/artifacts/remote_samples.json`](../2026-10-27-genomic-r2-upload/artifacts/remote_samples.json),
  which cover NG_, LRG, CM, NW_, NT_, chrM, GRCh38, NCBI36 and T2T. All 16
  match the expected VRS IDs (`logs/readback.log`), and the store's digests
  match the R2 copy's.
