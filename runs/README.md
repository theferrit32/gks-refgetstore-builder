# Run records

Durable records currently retained:

- [`2026-06-24-r2-upload/`](2026-06-24-r2-upload/) — reconstructed first R2 upload.
- [`2026-07-22-published-store/`](2026-07-22-published-store/) — authoritative currently published store.
- [`2026-07-23-ensembl-r116-genomic-experiment/`](2026-07-23-ensembl-r116-genomic-experiment/) — isolated, non-publishable genomic experiment.
- [`2026-08-26-full-rebuild/`](2026-08-26-full-rebuild/) — clean end-to-end rebuild
  (fresh cache, store, and lock) plus the seqrepo parity analysis run against it.

Removed: `2026-07-02-seqrepo-parity/`. Its build and parity analysis were fully
superseded by `2026-08-26-full-rebuild/`, which covers a larger source set and
regenerates the same reports; retaining both invited citing the stale numbers.

See [`RUNBOOK.md`](../RUNBOOK.md) for the lifecycle and retention policy.
