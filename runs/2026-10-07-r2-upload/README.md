# 2026-10-07 R2 upload: wrong shard layout

An upload of the store from
[`../2026-10-07-gtars-pr273-rebuild/`](../2026-10-07-gtars-pr273-rebuild/README.md)
to `theferrit32-public:theferrit32-public/refgetstore/2026-10-07`. The transfer
completed and `rclone check` passed, but the published layout was wrong: most
payloads were at keys no reader requests. Status **failed**. The prefix was
corrected in place by a later sync from a case-sensitive rebuild; see the
[runs index](../README.md).

## What ran

`upload_store.sh --remote theferrit32-public --bucket theferrit32-public
--prefix refgetstore/2026-10-07 --transfers 64 --checkers 64`, detached, with
rclone stats logged every minute (`RCLONE_LOG_FILE`, `RCLONE_STATS=1m`).

- Manifest: [`artifacts/published-manifest.json`](artifacts/published-manifest.json),
  written into the store and uploaded with it. It records commit `7aad2ae`
  (clean), gtars `0.11.0`, 1,779,361 objects and 23,807,108,941 bytes.
- Duration: 4 h 37 m (17:01–21:38 UTC).
- Seven `PutObject` requests failed with HTTP 502 from R2 in the first attempt.
  rclone's second attempt sent them; its first attempt therefore skipped the
  mirror-mode deletion step, which had nothing to delete in a new prefix.
- `rclone check`: 1,779,362 matching files, 0 differences
  ([`artifacts/rclone-summary.txt`](artifacts/rclone-summary.txt)).

## Why the layout was wrong

The store was built on a case-insensitive APFS volume. gtars stores each payload
at `sequences/<first two digest characters>/<digest>.seq`, and base64url digests
make `zz`, `Zz`, `zZ` and `ZZ` four shards. On that volume they were one
directory, named after whichever prefix was created first, so the store had
1,444 shard directories instead of 4,096. Locally nothing failed, since lookups
ignore case. rclone uploaded the real directory names, and R2 keys are
case-sensitive:

- 1,151,582 of 1,779,052 payloads (64.7%) were uploaded under a shard whose name
  differs in case from their digest prefix
  ([`artifacts/uploaded-layout.txt`](artifacts/uploaded-layout.txt)).
- Reading over the public URL with `RefgetStore.open_remote`, GPX4
  (`z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV`) returned HTTP 404: it was at
  `sequences/Z_/…`, and readers request `sequences/z_/…`.
- `rclone check` passed because it compares the local tree with the bucket, and
  those agreed.

The published `refgetstore/2026-07-22` store has the same problem. An affected
payload there returns 404 at its expected key and 206 at its uploaded key.

## What changed because of it

- Stores are built on a case-sensitive APFS volume, with `./store` a symlink to
  it (README, "Case-sensitive store directory").
- `build`, `sync --apply` and `repair --apply` refuse a store directory on a
  case-insensitive filesystem.
- `upload_store.sh` logs rclone stats periodically when run without a terminal,
  which is what made this run's progress visible.
