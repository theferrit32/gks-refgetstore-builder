# Run records

Durable records currently retained:

- [`2026-06-24-r2-upload/`](2026-06-24-r2-upload/) — reconstructed first R2 upload.
- [`2026-07-02-seqrepo-parity/`](2026-07-02-seqrepo-parity/) — full build and seqrepo parity analysis.
- [`2026-07-22-published-store/`](2026-07-22-published-store/) — authoritative currently published store.
- [`2026-07-23-ensembl-r116-genomic-experiment/`](2026-07-23-ensembl-r116-genomic-experiment/) — isolated, non-publishable genomic experiment.

See [`RUNBOOK.md`](../RUNBOOK.md) for the lifecycle and retention policy.

## Future commit candidates

After review, suitable commit candidates are the run framework and runbook,
the four run manifests and compact summaries, the preserved small logs and
analysis outputs, the experiment evidence/configuration, and the reusable tools
under `seqrepo_equivalence/`. Generated stores, downloads, bulk parity tables,
and caches remain local and ignored. Nothing in this organization pass is
staged or committed.
