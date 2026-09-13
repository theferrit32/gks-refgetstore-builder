# The build lock pins inputs but not the ingested outcome

**Status: CLOSED** — resolved by `ingest_spec` (schema /3) and by
`outputs` + `verify --store` (schema /4). See [Resolution](#resolution).

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

---

## Resolution

Closed in two steps, because the issue turned out to contain two questions
rather than one.

**Schema /3 pinned the *intent*.** `inputs.files[].ingest_spec` records the
declared transformation per file — derived format, and the record-exclusion
rule when it applies to *that* file class — with `ingest_spec_sha256` as its
canonical digest. `plan_sync` classifies a source as `reingest` when either the
bytes or the spec moved, so editing an exclusion rule over unchanged upstream
bytes is now detected.

**Schema /4 pinned the *outcome*.** `outputs` is a census enumerated from the
store: counts, a canonical root over the sequence digest set, a root over the
collection digest set, and every collection with its contributing files.
`verify --store` checks it. This is what the issue was actually asking for — the
realized result, content-addressed, immune to compensating errors.

### The four design questions, answered

**Where the check runs.** Not in the preflight — the issue was right that a
preflight cannot know a digest before ingest. It runs in `verify --store`, as
its own command, against the lock a previous build wrote. `build` still gates on
*inputs* before opening the store; `verify` judges outputs afterwards. Splitting
them is what lets outcome verification be offline, repeatable, and runnable
without a build.

**How an intentional change is expressed.** Through the existing lock-write
path, not a new mode. An intentional change moves `ingest_spec_sha256`, which
makes `sync` re-ingest and rewrite `outputs` from the resulting store. The
strict/informational split the issue anticipated proved unnecessary: a changed
outcome with an unchanged spec is a regression, and a changed outcome with a
changed spec is already explained by the spec.

**Whether the lock should record the exclusion rule itself.** Yes, and it does —
`ingest_spec` carries the literal `record_prefixes`. That was the point:
`sources_toml_sha256` detects *that* the config changed, not *what*, so a
one-character typo in a prefix and an unrelated comment edit were
indistinguishable. Per-file specs make a deliberate change attributable to the
file it affects, and keep a `dna.toplevel` edit from marking its `cdna` sibling
as changed.

**First builds and new sources.** They have no prior digest, and that is
recorded rather than papered over: a full build writes `outputs` fresh from the
store, and a collection whose contributor is unknown is `"from": []` rather than
omitted. `verify --store` compares the store to *whatever the lock last
recorded*, so a first build establishes the baseline and the second build
onwards is checked against it.

**Scope.** The gate covers every collection, not only transformed sources. The
roots are over the complete digest sets, so a gtars upgrade that encodes or
digests differently moves `sequences_root` regardless of which source it came
from — which is exactly the case the issue named that a per-source check would
have missed.

### What it caught immediately

`verify --store --deep` found 133 sequences whose stored payload does not
re-digest to the digest it is filed under: 106 selenoproteins losing `U` because
gtars' protein alphabet lacks it, and 27 short proteins misdetected as
nucleotide. Not bit rot, not a filter regression, and not fixable by
re-ingesting — a real data-correctness defect that had been invisible for the
entire life of the store. It is pinned as a baseline and reported on every
`--deep` run.

### The premise the issue inherited, and which was wrong

The issue quotes `rec["collection_digest"] = …` approvingly as the thing to
compare against. That field encoded one file → one collection, which was never
true: content-addressing merges byte-identical sources, and 43 collections in
the committed lock have two or more contributors. Schema /4 replaces it with
`outputs.collections[].from`, a list, in the honest direction.
