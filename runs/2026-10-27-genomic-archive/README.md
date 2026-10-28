# Publish the genomic-only store as a single archive

`refgetstore/genomic-2026-10-27.tar.gz` is the store published at
`refgetstore/genomic-2026-10-27/`
([`../2026-10-27-genomic-r2-upload/`](../2026-10-27-genomic-r2-upload/README.md))
as one gzipped tar, for consumers that want the whole store on local disk: one
download and one checksum instead of 13.7k requests. It adds no data. The
per-file prefix stays the layout gtars reads remotely.

- Public URL: <https://static.ferriter.dev/refgetstore/genomic-2026-10-27.tar.gz>
- SHA-256 `9b8903a19d8ff64b204b3b0a4cb9f9cc80026c88516d45a9caae2deb183a65ab`,
  3,535,731,383 bytes, published beside it as `genomic-2026-10-27.tar.gz.sha256`
- Contents: the single directory `genomic-2026-10-27/` holding the store's
  13,681 files and its `manifest.json`

## How it was made

macOS `tar` adds AppleDouble `._*` entries and extended attributes by default,
which would extract as stray files on Linux; the archive omits both and records
`root:root` ownership:

    COPYFILE_DISABLE=1 tar --no-mac-metadata --no-xattrs --no-acls \
        --uid 0 --gid 0 --uname root --gname root \
        -C /Volumes/RefgetStores -cf - genomic-2026-10-27 | pigz -6

The first archive (`create-archive`) recorded the top-level directory as mode
0700, inherited from the `mktemp` directory the store was built in. The
directory was set to 0755 and the archive rebuilt (`create-archive-mode-755`);
only the rebuilt archive was uploaded.

## Validation

- Extracting the uploaded archive gives a tree identical to the store
  directory (`diff -rq`).
- `rclone check`: 0 differences across the archive and its `.sha256`.
- The public URL serves the full content length and the `.sha256` serves the
  archive's digest. clinvar-ingest's `vrsify` step downloaded it by default,
  verified the digest, extracted it, and translated the 2026-10-24
  variation_identity delta to output identical to a run on the local store.
