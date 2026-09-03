# The build lock pins inputs but not the ingested outcome

**Labels:** `provenance` · `build-lock` · `verification`

## Summary

`build.lock.json` records, per source, both the upstream `sha256` and the
resulting `collection_digest` (`build_lock.py:129-196`). The drift preflight
compares only `sha256` — see `evaluate_cache_vs_lock` (`build_lock.py:282`) and
`evaluate_sources_vs_lock` (`:308`), driven from `build_store.py:1806-1830`.

Anything that changes *what gets ingested* from *unchanged upstream bytes* is
therefore invisible to the check:

- a `sources.toml` record-exclusion rule edited, added, or removed
- a regression in the exclusion filter — a wrong prefix, an unanticipated header
  format
- a change in a derived-FASTA converter (`gbff_to_fasta`, `lrg_zip_to_fasta`,
  `filter_fasta_records`)
- a gtars upgrade that digests or encodes differently

In each case the lock preflight passes, the build produces a different store,
and a fresh lock is written over the old one with no signal that the outcome
moved.

## Why a record count is not the fix

An earlier design pinned an expected record count per release in `sources.toml`.
It was dropped, for two reasons that are worth recording so it is not
reintroduced:

1. **Upstream changes are already covered, and covered better.** If a provider
   republishes a file such that its record count changes, its `sha256` changes,
   and the preflight flags that exact file *before* ingest.
2. **It cannot see a compensating change.** One record fewer from one file and
   one more from another leaves the total identical. A count is a lossy proxy
   for the outcome.

It also required hand-maintained numbers across 34 near-identical config blocks,
which is a standing maintenance burden and a spurious-failure source.

## Proposed direction

Compare the realized `collection_digest` against the value the lock already
records. It is content-addressed, so it pins the actual ingested result and is
immune to compensating errors. What exists today:

```python
# build_lock.py:194-195
rec["collection_digest"] = coll["collection_digest"] if coll else None
rec["n_sequences"] = coll["n_sequences"] if coll else None
```

Design questions to settle rather than assume:

- **Where the check runs.** A preflight cannot know the digest before ingest, so
  this is either a post-ingest gate or a mode of `gks-refgetstore verify`.
- **How an intentional change is expressed.** Editing an exclusion rule
  legitimately moves the digest. The check needs a strict/informational split in
  the shape of the existing `--lock-check-mode {strict,subset,ignore}`
  (`cli.py:59-63`).
- **Whether the lock should record the exclusion rule itself** — file classes and
  record prefixes — so a deliberate config change is attributable separately from
  a code regression. `build.sources_toml.sha256` currently detects *that* the
  config changed, not *what* changed or *where*.
- **First builds and new sources**, which have no prior digest to compare.
- **Scope.** This is not specific to record exclusion; it applies to every
  derived-FASTA path, and to gtars upgrades. Worth deciding whether the gate
  covers all collections or only sources with a declared transformation.

## Current mitigations

Not nothing, but not a pin:

- A declared exclusion matching **zero** records raises
  (`build_store.filter_fasta_records`), catching the most likely regression — the
  filter silently no-ops and retains everything it was meant to remove.
- Every dropped record is written to `<artifact>.excluded.tsv`, so the outcome is
  inspectable rather than inferred.
- `SeqsetStats.records_excluded` surfaces per-seqset counts in the build summary.

None of these detect a filter that drops the *wrong* records while dropping the
right number.

## Related

- `seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md` — the exclusion this
  surfaced from.
- `ISSUE-ensembl-padded-scaffold-bloat.md` — the originating issue.
