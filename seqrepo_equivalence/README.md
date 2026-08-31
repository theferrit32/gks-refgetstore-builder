# seqrepo equivalence

Tooling to check whether a built RefgetStore is **backwards-compatible with a
biocommons seqrepo snapshot** — i.e. does it contain (by digest) the sequences
seqrepo knows, and do the namespaces/aliases line up. This is an *analysis /
regression* aid, separate from the core `gks-refgetstore` build: the builder
loads authoritative current sources, and seqrepo accumulates historical
versions, so some divergence is expected and documented here rather than treated
as a build defect.

The store's `lrg` namespace is store-only by design and is expected to land
entirely in `refget_only`: the `2024-12-20` seqrepo snapshot has no LRG namespace
at all, so `NAMESPACE_MAP` deliberately omits an `LRG` entry.

## Contents

- `verify_seqrepo_equivalence.py` — exhaustively compares a store against a
  seqrepo snapshot. Reads seqrepo's `aliases.sqlite3` directly via stdlib
  `sqlite3` and the store's on-disk alias TSVs (gtars-only otherwise). Every
  difference is categorized (`matched`, digest mismatch, `superseded_old_version`,
  `backfill_candidate`, `alias_naming_gap`, `sequence_missing`) into a JSON report
  plus a `build_gaps.tsv` feed, ending with a PASS/FAIL verdict.
- `known_divergence/` — the regression fixture for Ensembl accessions whose
  seqrepo digest differs from the gtars/refget digest **by design**, plus the
  script that generated it. Without the fixture those mismatches would be
  reported as *unexplained* and flip the verdict to fail. It is pinned to the
  `2024-12-20` snapshot **and** to store namespace `ensembl-113`; see that
  directory's README.
- `ENSEMBL_RELEASE_DRIFT.md` — measures how often Ensembl republishes the same
  `accession.version` with different bytes, across all 42 release namespaces.
  Read this before drawing conclusions from any digest mismatch: it explains why
  mismatch-vs-seqrepo is a **weak** drift metric, and why the rolling `ensembl`
  namespace must not be used to validate anything pinned.
- `parity_membership.py` — emits the exhaustive digest-membership table and a
  compact Markdown summary. Uniquely produces the seqrepo digest **coverage %**
  headline and the per-prefix (`NM_`/`NP_`/`XM_`…) rollups.
- `full_parity.py` — emits the exhaustive biological alias-union table, every
  missing/mismatched SeqRepo alias, shared-alias mismatches, and earliest-source
  attribution for every store-only digest.
- `probe_source_coverage.py` — measures candidate-source coverage against a
  verifier gap table while caching compact accession lists.

## Usage

Run from the repo root (defaults point at `./store`, the 2024-12-20 snapshot,
and the `runs/2026-08-26-full-rebuild/parity/` output paths):

    uv run python seqrepo_equivalence/verify_seqrepo_equivalence.py -v

Key flags (see `--help`): `--seqrepo PATH` (snapshot dir), `--store PATH`
(default `./store`), `--known-divergent PATH` (default the fixture here),
`--report PATH`, `--gap-list PATH`, `--fail-on {none,gaps,mismatch}`.

If you point this at a different seqrepo snapshot, build a **new** fixture in
`known_divergence/` for that version rather than editing the existing one — the
set of `*`-divergent accessions changes between snapshots, and Ensembl has been
observed stripping `*` upstream too (release 114), which retires cases entirely.

Bulk parity, gap, mismatch, contribution, and inventory TSV outputs are
reproducible and ignored. Retain them locally and record their schemas, row
counts, byte sizes, SHA-256 values, generator commands, and headline categories
in the relevant run manifest.
