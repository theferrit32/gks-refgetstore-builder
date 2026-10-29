# Publish the genomic SeqRepo as a single archive

`seqrepo/genomic-2026-10-27.tar.gz` is the SeqRepo instance built from the
genomic store ([`../2026-10-29-genomic-seqrepo/`](../2026-10-29-genomic-seqrepo/README.md))
as one gzipped tar. It holds the same sequences and aliases as
`refgetstore/genomic-2026-10-27`, so the two backends can be compared on
identical data. It adds no data.

- Public URL: <https://storage.googleapis.com/clingen-public/seqrepo/genomic-2026-10-27.tar.gz>
  (`gs://clingen-public/seqrepo/genomic-2026-10-27.tar.gz`)
- SHA-256 `bdb0f01b76538e39b9f3a71ca52b4252c3630f5350c9b7b01dd4db66f35d9a93`,
  3,501,905,945 bytes, published beside it as `genomic-2026-10-27.tar.gz.sha256`
- The build's provenance (verification counts, versions, every file's SHA-256) is
  published beside it as `genomic-2026-10-27.provenance.json` and kept here in
  [`artifacts/`](artifacts/)
- Contents: the single directory `genomic-2026-10-27/` holding the instance's
  41 files, with `aliases.sqlite3` at its top. clinvar-ingest's `vrsify
  --seqrepo` accepts the archive URL: it verifies the `.sha256`, extracts it
  and finds the instance one level down.

## How it was made

With the same flags as the store archive
([`../2026-10-27-genomic-archive/`](../2026-10-27-genomic-archive/README.md)):
no AppleDouble entries or extended attributes, `root:root` ownership. The
instance directory was already mode 0755.

    COPYFILE_DISABLE=1 tar --no-mac-metadata --no-xattrs --no-acls \
        --uid 0 --gid 0 --uname root --gname root \
        -C /Volumes/RefgetStores/seqrepo -cf - genomic-2026-10-27 | pigz -6

It took 16 s and is barely smaller than the instance (3.50 GB from 3.53 GB):
seqrepo already stores its sequences bgzip-compressed.

The three files were uploaded with `rclone copy --checksum` through a GCS
remote defined only by `RCLONE_CONFIG_*` environment variables (gcloud
application-default credentials).

## Validation

- Extracting the archive gives a tree identical to the instance directory
  (`diff -rq`).
- `rclone check`: 0 differences across the archive, its `.sha256` and the
  provenance JSON.
- Unauthenticated HTTPS serves all three at full length, and the served
  `.sha256` holds the archive's digest.
