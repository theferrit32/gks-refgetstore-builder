# seqrepo equivalence

Tooling to check whether a built RefgetStore is **backwards-compatible with a
biocommons seqrepo snapshot** — i.e. does it contain (by digest) the sequences
seqrepo knows, and do the namespaces/aliases line up. This is an *analysis /
regression* aid, separate from the core `gks-refgetstore` build: the builder
loads authoritative current sources, and seqrepo accumulates historical
versions, so some divergence is expected and documented here rather than treated
as a build defect.

## Contents

- `verify_seqrepo_equivalence.py` — exhaustively compares a store against a
  seqrepo snapshot. Reads seqrepo's `aliases.sqlite3` directly via stdlib
  `sqlite3` and the store's on-disk alias TSVs (gtars-only otherwise). Every
  difference is categorized (`matched`, digest mismatch, `superseded_old_version`,
  `backfill_candidate`, `alias_naming_gap`, `sequence_missing`) into a JSON report
  plus a `build_gaps.tsv` feed, ending with a PASS/FAIL verdict.
- `ensembl_known_divergent.txt` — regression fixture: Ensembl accessions whose
  seqrepo digest differs from the gtars/refget digest **by design** (the ~116
  ENSP sequences with embedded `*` stop codons, which seqrepo strips before
  digesting and gtars does not — see the root README). Without this fixture those
  mismatches would be reported as *unexplained* and flip the verdict to fail.
  **This fixture corresponds to the biocommons seqrepo `2024-12-20` snapshot.**

## Usage

Run from the repo root (defaults point at `../store` and the 2024-12-20 snapshot):

    uv run python seqrepo_equivalence/verify_seqrepo_equivalence.py -v

Key flags (see `--help`): `--seqrepo PATH` (snapshot dir), `--store PATH`
(default `./store`), `--known-divergent PATH` (default the fixture here),
`--report PATH`, `--gap-list PATH`, `--fail-on {none,gaps,mismatch}`.

If you point this at a different seqrepo snapshot, regenerate
`ensembl_known_divergent.txt` for that version — the set of `*`-divergent
accessions can change between snapshots.
