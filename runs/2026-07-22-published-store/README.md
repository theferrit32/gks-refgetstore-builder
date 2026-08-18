# 2026-07-22 published store

This is the durable record for the currently published store. The authoritative
manifest remains unchanged at `store.2026-07-22/manifest.json`; its byte-identical
copy is [`published-manifest.json`](published-manifest.json), SHA-256
`36be27a75988684a738d07527ab0bd5cc2f186a381fbce295d10ca0ac5b0c62a`.

- Store: 1,213,617 sequences and 107 collections.
- Locked inputs: 129 in [`build.lock.json`](build.lock.json).
- Upload: 1,213,754 store objects totaling 5,646,134,737 bytes.
- Destination: `theferrit32-public:theferrit32-public/refgetstore/2026-07-22`.
- Verification: local root manifest/index hashes were verified after preserving
  the store; the repository verifier passed 67/67; the upload was checksum
  verified. The exact verifier and upload logs are missing.
- Upload metadata: clean commit `dc04e5e285a4f8a294298d58fb2d5814adf4fe8a`.

## Dirty build lock versus clean upload

The build lock truthfully records commit `9a213c4d1227cdda8068fb140b54354b9ec8d03c`
as `9a213c4-dirty`: the build ran with later, uncommitted assembly-report alias
changes in the worktree. Those changes were subsequently reviewed and committed
as clean commit `dc04e5e285a4f8a294298d58fb2d5814adf4fe8a`, which is what the upload
manifest records. The differing metadata describes two times in the same
build-to-publication chain; it does not imply two different published stores.

Neither the exact build log nor the exact upload log survives. This record does
not recreate them or imply that reconstructed logs are original evidence.
