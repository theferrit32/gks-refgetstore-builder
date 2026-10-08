# Rebuild onto a case-sensitive volume

The store was rebuilt onto a case-sensitive APFS volume so that every sequence
payload sits in the shard directory named by its own digest prefix. Its contents
are unchanged from
[`../2026-10-07-gtars-pr273-rebuild/`](../2026-10-07-gtars-pr273-rebuild/README.md);
only the directory layout differs.

## Why

gtars stores each payload at `sequences/<first two digest characters>/<digest>.seq`.
The previous store was written to a case-insensitive APFS volume, where
prefixes differing only in case (`zz`, `Zz`, `zZ`, `ZZ`) share one directory.
It had 1,444 shard directories where a full store needs 4,096, and 1,151,582
payloads (64.7%) in a directory named for a different-case prefix. That reads
correctly in place, but not once uploaded: see
[`../2026-10-07-r2-upload/`](../2026-10-07-r2-upload/README.md).

## Setup

A case-sensitive APFS volume in the internal disk's container, sharing its free
space:

    diskutil apfs addVolume disk3 APFSX RefgetStores

It mounts at `/Volumes/RefgetStores`. A smoke build of the LRG seqset there
produced 2,765 shard directories for 4,516 sequences, more than the 1,444 a
case-insensitive filesystem allows, with every payload in its own prefix's
directory.

## Commands

| label | exit | what |
| --- | --- | --- |
| `tests` | 0 | 244 passed |
| `build` | 0 | `build --store-dir /Volumes/RefgetStores/store --locked-sources`, from a copy of the committed lock and the existing `downloads/` cache; 25 minutes |
| `verify-deep` | 0 | `verify --store --deep` against the empty known-bad baseline: all 1,779,052 sequences re-digest to their own digest, 0 errors, 0 warnings |

## Outcome

[`artifacts/comparison.txt`](artifacts/comparison.txt):

- The lock this build wrote ([`artifacts/build.lock.json`](artifacts/build.lock.json))
  has the same inputs, digest roots, counts and collection list as the
  committed `build.lock.json`. The committed lock is kept, since it describes
  this store exactly.
- 4,096 shard directories; no payload is in a directory that differs from its
  digest prefix.
- Every one of the 1,779,052 payloads is byte-identical to the previous
  store's, as are `sequences.rgsi`, `collections.rgci` and the collection
  files. `rgstore.json` differs only in its `created_at` and `modified`
  timestamps.

The written lock records the worktree as dirty (`7aad2ae-dirty`). The lock
takes git state when it is written, at the end of the build, and by then the
case-sensitivity check was being added to the worktree. The build ran clean
`7aad2ae` code: this record was initialized with `--require-clean` immediately
before it.

`./store` is now a symlink to `/Volumes/RefgetStores/store`, so every command's
default `--store-dir` reaches this store. The previous store was moved aside.
