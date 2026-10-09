# Run records

Durable records currently retained:

- [`2026-06-24-r2-upload/`](2026-06-24-r2-upload/) — reconstructed first R2 upload.
- [`2026-07-22-published-store/`](2026-07-22-published-store/) — authoritative currently published store.
- [`2026-07-23-ensembl-r116-genomic-experiment/`](2026-07-23-ensembl-r116-genomic-experiment/) — isolated, non-publishable genomic experiment.
- [`2026-08-26-full-rebuild/`](2026-08-26-full-rebuild/) — clean end-to-end rebuild
  (fresh cache, store, and lock) plus the seqrepo parity analysis run against it.
- [`2026-09-12-lock-schema-v4/`](2026-09-12-lock-schema-v4/) — the committed build
  lock converted from schema `/2` to `/4`. Converted, not rebuilt: a rebuild
  would have destroyed the preserved upstream drift.
- [`2026-09-12-deep-verify/`](2026-09-12-deep-verify/) — full round-trip scan of
  the store, producing the `verify --deep` encoder baseline. Found 133 sequences
  whose stored payload does not re-digest to its own digest.
- [`2026-10-07-gtars-pr273-rebuild/`](2026-10-07-gtars-pr273-rebuild/) — the
  current lock: a full rebuild from the locked cache with gtars from
  databio/gtars#273. Every sequence round-trips; the encoder baseline is empty.
- [`2026-10-07-r2-upload/`](2026-10-07-r2-upload/) — **failed**: that store's
  upload to `refgetstore/2026-10-07`. It transferred and checked cleanly, but the
  store had been built on a case-insensitive filesystem, so 64.7% of payloads
  were published under a shard key no reader requests. `refgetstore/2026-07-22`
  has the same problem.
- [`2026-10-08-case-sensitive-rebuild/`](2026-10-08-case-sensitive-rebuild/) —
  the current store: the same contents rebuilt onto a case-sensitive APFS
  volume, which `./store` now links to. Byte-identical payloads, correct shard
  layout.
- [`2026-10-08-r2-resync/`](2026-10-08-r2-resync/) — `refgetstore/2026-10-07`
  corrected in place from that store: `rclone check` reports 0 differences, and
  reads over the public URL return the right residues.

Removed: `2026-07-02-seqrepo-parity/`. Its build and parity analysis were fully
superseded by `2026-08-26-full-rebuild/`, which covers a larger source set and
regenerates the same reports; retaining both invited citing the stale numbers.

See [`RUNBOOK.md`](../RUNBOOK.md) for the lifecycle and retention policy.
