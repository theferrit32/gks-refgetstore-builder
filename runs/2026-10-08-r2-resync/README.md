# Corrective sync of refgetstore/2026-10-07

`theferrit32-public:theferrit32-public/refgetstore/2026-10-07` was synced from the
case-sensitive rebuild in
[`../2026-10-08-case-sensitive-rebuild/`](../2026-10-08-case-sensitive-rebuild/README.md),
replacing the wrongly sharded upload recorded in
[`../2026-10-07-r2-upload/`](../2026-10-07-r2-upload/README.md). The prefix now
holds exactly that store, with every payload at the key readers request.

Public URL: <https://static.ferriter.dev/refgetstore/2026-10-07/>

## What ran

`upload_store.sh --remote theferrit32-public --bucket theferrit32-public
--prefix refgetstore/2026-10-07 --store-dir /Volumes/RefgetStores/store
--transfers 64 --checkers 64`, detached, with `RCLONE_CHECKSUM=true` so rclone
compared MD5s from the bucket listing rather than modification times. Without
it, the ~627,000 objects already at their correct keys would have had their
metadata rewritten for new modification times.

- Manifest: [`artifacts/published-manifest.json`](artifacts/published-manifest.json),
  uploaded with the store. Commit `51a1b9d` (clean), gtars `0.11.0`, 1,779,361
  objects, 23,807,108,941 bytes.
- A dry run first predicted 1,151,583 copies (1,151,582 payloads at their
  correct keys plus `rgstore.json`, whose timestamps changed) and 1,151,583
  deletions (the 1,151,582 wrong-key payloads plus the previous
  `manifest.json`), and nothing else.
- `rclone sync` mirrors: after an error-free transfer it deletes every object
  under the prefix that is not in the source. No request failed, so the
  deletions ran.
- Duration: 2 h (22:12–00:12 UTC).

## Verification

[`artifacts/rclone-summary.txt`](artifacts/rclone-summary.txt):

- `rclone check`: **0 differences, 1,779,362 matching files** (the store's
  1,779,361 files plus `manifest.json`). It counts objects present only in the
  bucket as differences, so the prefix holds exactly the store and nothing left
  over.
- No error lines in the rclone log.

[`artifacts/remote-checks.txt`](artifacts/remote-checks.txt):

- 4,096 shard prefixes under `sequences/` in the bucket.
- `RefgetStore.open_remote` on the public URL, with the pinned gtars: GPX4
  (`z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV`) returns `U` at index 109 and
  `NR_103745.1` returns `D` at index 852; both round-trip. GPX4 returned HTTP 404
  from this prefix before the sync.
- 60 sampled payloads with mixed- or upper-case digest prefixes: HTTP 206 at the
  correct key, HTTP 404 at a different-case key.

## Reading this store

Read it with gtars at `theferrit32/gtars@7831e09f` or a later release that
includes databio/gtars#273. Older gtars opens it without error and returns
selenocysteine `U` as `X`.

`refgetstore/2026-07-22`, the previously published store, still has the
case-folded layout and was not changed by this run.
