# Retire refgetstore/2026-06-24 and refgetstore/2026-07-22

The two superseded prefixes were deleted from
`theferrit32-public:theferrit32-public/refgetstore/`, leaving only the published
`refgetstore/2026-10-07`
([`../2026-10-08-r2-resync/`](../2026-10-08-r2-resync/README.md)).

Both were built on a case-insensitive filesystem, so most of their payloads were
published under shard keys that readers do not request. Both were also built
with gtars older than databio/gtars#273. Their records are kept:
[`../2026-06-24-r2-upload/`](../2026-06-24-r2-upload/README.md) and
[`../2026-07-22-published-store/`](../2026-07-22-published-store/README.md).

## Commands

| label | what | result |
| --- | --- | --- |
| `size-before-2026-06-24` | `rclone size` | 783,641 objects, 3,582,897,034 bytes |
| `size-before-2026-07-22` | `rclone size` | 1,213,755 objects, 5,646,144,318 bytes |
| `size-before-2026-10-07` | `rclone size` (baseline, must not change) | 1,779,362 objects, 23,807,123,936 bytes |
| `dryrun-delete-2026-06-24` | `rclone delete --dry-run` | 783,641 objects, all within the store layout |
| `dryrun-delete-2026-07-22` | `rclone delete --dry-run` | 1,213,755 objects, all within the store layout |
| `purge-2026-06-24` | `rclone purge` | 783,641 objects deleted, 3.337 GiB freed, 24 m |
| `purge-2026-07-22` | `rclone purge` | 1,213,755 objects deleted, 5.258 GiB freed |
| `prefixes-after` | `rclone lsf --dirs-only refgetstore/` | `2026-10-07/` only |
| `size-after-2026-10-07` | `rclone size` | 1,779,362 objects, 23,807,123,936 bytes: unchanged |

Each purge logged one error-level line at startup: rclone could not read the
bucket's versioning status (R2 refuses `GetBucketVersioning`) and proceeded as
unversioned. No deletion failed; both purges exited 0.

After the purges, `https://static.ferriter.dev/refgetstore/2026-06-24/rgstore.json`
and `…/2026-07-22/rgstore.json` return HTTP 404, and `…/2026-10-07/rgstore.json`
returns 200.

## Local stores removed at the same time

Only the current store (`store` → `/Volumes/RefgetStores/store`, built from the
current `sources.toml`) and the download cache were kept. Removed:

| path | size | what it was |
| --- | ---: | --- |
| `store.pre-pr273/` | 49 GB | the gtars 0.9.2 build from the previous manifest |
| `store.pre-filter/` | 28 GB | the current digest set, case-folded |
| `store.2026-07-22/` | 9.2 GB | source of `refgetstore/2026-07-22` |
| `store.dev/` | 6.0 GB | dev-manifest store; rebuilds in minutes |
| `store.case-folded-2026-10-07` (Trash) | 28 GB | source of the first `refgetstore/2026-10-07` upload |

Free space in the disk's APFS container went from 131 GiB to 231 GiB.
