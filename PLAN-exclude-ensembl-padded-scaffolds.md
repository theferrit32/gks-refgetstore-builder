# Exclude Ensembl's N-padded alt/patch scaffolds at ingest

## Context

Ensembl releases 76–109 publish alt loci and patch scaffolds in `dna.toplevel`
padded to their parent chromosome's length with `N`, named with a `CHR_` prefix.
Release 110 switched to true-length records. The store ingests releases 75–116,
so it holds **445 padded records / 60,047,447,030 bases**. Because `N` forces the
`dna3bit` alphabet (3 bits/base), that is **20.97 GiB — 43% of the 48.9 GiB
store** — encoding nothing the store does not already hold at true length.

All 445 have a true-length counterpart already present; re-hashing the real span
of each padded record reproduces the unpadded digest exactly (387 forward, 58
reverse-complemented, 0 failures). The padded records carry no accession and no
alias outside `ensembl-76`–`ensembl-109`.

Background: `seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md`,
`ISSUE-ensembl-padded-scaffold-bloat.md`, and the per-record data in
`seqrepo_equivalence/ensembl_padded_scaffolds.tsv`.

Outcome: store drops to ~27.9 GiB with no loss of sequence content.

## Decisions

**Declarative, not heuristic.** No N-fraction inspection. Exclusion is driven by
record *name*, declared per seqset. Verified across releases 76, 90, 100 and 109:
zero padded records lack the `CHR_` prefix, zero prefixed records are unpadded,
and releases 110–116 contain no `CHR_` records at all. The marker is both sound
and complete.

**No synthesized aliases.** If release 100 did not publish `HSCHR1_2_CTG3`, we do
not invent `ensembl-100:HSCHR1_2_CTG3`. Excluded records simply vanish from those
release namespaces. The rolling `ensembl:HSCHR1_2_CTG3` alias still exists via
release 116, and the sequence remains reachable through `refseq`/`insdc`.

**Release 75 is untouched** — it is GRCh37 and has zero `CHR_` records, so its
11.19 GiB of genomic sequence is unaffected.

## Config shape

Nested inline table on each of the 34 affected seqsets (`sources.toml`, releases
76–109 only):

```toml
exclude = { file_classes = ["dna.toplevel"], record_prefixes = ["CHR_"] }
```

**No expected-count field.** An earlier draft pinned a per-release record count
in config. It is not worth it: if upstream republishes a file such that the count
changes, its `sha256` changes, and the lock preflight already flags that exact
file *before* ingest — earlier and better attributed than a post-hoc count check.
A hand-maintained count across 34 releases is pure maintenance burden and a
spurious-failure source.

A count would also not catch a *compensating* regression — one record fewer from
one file and one more from another leaves the total unchanged — so it is a weak
proxy for the outcome even in the case it was meant to guard.

Two cheap measures instead, neither needing per-release numbers:

- **Non-empty invariant.** If a file class declares an exclusion and *zero*
  records match, raise. This catches the realistic regression — the filter
  silently no-ops and 21 GiB of padding is quietly retained. No configuration,
  no maintenance.
- **The audit log** (below) records every dropped record, so the outcome is
  always inspectable rather than inferred.

The residual gap — a filter regression that changes *which* records are dropped
while upstream bytes stay identical — is not closed by either, and is tracked as
a follow-up (see **Follow-up issue** below).

- `SeqsetConfig` (`build_store.py:68-121`) gains `exclude: dict | None = None`.
  TOML nested tables arrive as plain dicts through `SeqsetConfig(**entry)`
  (`load_config`, `build_store.py:525-561`), which rejects unknown keys with a
  `TypeError` — so the field must be declared.
- `__post_init__` (`build_store.py:93-170`) converts it to a frozen
  `RecordExclusion` dataclass and validates: unknown keys rejected,
  `record_prefixes` and `file_classes` both non-empty, and every entry in
  `file_classes` present in this seqset's own `file_class`/`file_classes`.
  Follow the existing `raise ValueError(f"seqset {self.name!r}: …")` style used
  by the 13 guards already there.

## Implementation

**1. The filter — `build_store.py`, beside the other derived-FASTA writers.**

```python
def filter_fasta_records(src: Path, out_path: Path, prefixes: tuple[str, ...],
                         excluded_log: Path) -> tuple[int, int]:
    """Copy src to out_path, dropping records whose name matches a prefix.

    Returns (kept, dropped). Writes every dropped record's name and length to
    excluded_log so exclusions are auditable rather than silent.
    """
```

Stream gzip→gzip, copying kept records' lines verbatim. Mirror the established
`.part` + `replace()` idiom from `gbff_to_fasta` (`build_store.py:627-664`) and
`lrg_zip_to_fasta` (`:683-756`).

Write **gzip level 1**, not the `gzip` module default. See the benchmark below:
level 6 costs an extra 347s per file to save 140 MB, and Python's `gzip.open`
defaults to level 9.

**2. The resolver.** `resolve_filtered_fasta(src, exclusion, force) -> Path`,
matching the `(Path, bool) -> Path` contract and cache semantics of
`resolve_gbff_fasta` (`:667-676`): reuse when the output exists, is non-empty,
and `not force`. The output is a **sibling of the source artifact in the download
cache**, exactly like `<artifact>.fasta` for gbff and lrg_zip — same convention,
same reuse rules, differing only in that it is gzipped because it is ~30x larger
than those.

**Naming — this is the trap.** `build_lock.records_from_log`
(`build_lock.py:94-119`) recovers a source cache path by stripping a trailing
`.fasta`, which works only because existing resolvers *append* that suffix. Write
`<artifact>.filtered.fa.gz` and extend `records_from_log` to strip a shared
`DERIVED_SUFFIXES` tuple used by both sides, rather than the current hardcoded
`.fasta`. This path is only exercised by `gks-refgetstore lock --from-log`;
normal builds populate provenance directly at `build_store.py:1449-1456`.
Leaving it unfixed would fail silently.

**3. Filtration is a preflight pass, run in parallel.** Downloads already work
this way — `run_build` calls `ensure_download` for every source at
`build_store.py:1790-1802`, *before* the seqset loop, and `process_seqset` then
consumes cached files. Filtration becomes a second preflight in the same shape:

```python
def prepare_filtered_sources(seqsets, cache_dir, jobs, force) -> dict[Path, Path]:
    """Materialize every declared filtered FASTA. Returns source -> filtered."""
```

Run it over a `ThreadPoolExecutor`. `zlib` releases the GIL during
decompress/compress, so threads scale nearly linearly here — the same reason
`--ingest-jobs` works. This matters: decompression dominates at ~245s per
release, so 34 serial releases would be ~2.3 hours, versus roughly 20 minutes at
8 workers.

Add a separate **`--filter-jobs`** (default `INGEST_JOBS_DEFAULT`, i.e.
`min(8, cpu_count)`) rather than reusing `--ingest-jobs`. The two have different
bottlenecks — filtration is zlib/IO-bound, ingest is Rust-CPU- and memory-bound —
and `--ingest-jobs` is documented as the knob for capping peak memory
(`cli.py:46-51`). Someone lowering it to control memory should not silently
serialize filtration.

**4. The hook — `process_seqset` phase 1 (`build_store.py:1399-1426`).** After
`require_prepared_source` and the `DERIVED_FASTA_RESOLVERS` step, swap in the
filtered path when the current URL's file class is declared in `exclude`. By this
point the preflight has already produced the file, so this is a cache hit;
resolving here rather than in the preflight keeps `process_seqset` usable
standalone (as the tests drive it). `iter_shard_urls()` (`:171-188`) drops
per-index metadata, so enumerate it and add a small
`SeqsetConfig.file_class_for_index(i)` helper returning `file_classes[i]` or
`file_class`.

Keep provenance keyed on the **source cache path** (`target`), not the filtered
file — the comment at `:1451-1453` explains why, and the lock's `sha256` must
keep covering upstream bytes.

**4. Failure behaviour.** A zero-match violation or IO error raises, which the
existing `except Exception` at `:1418-1424` converts to `stats.warnings += 1`,
tripping the zero-warning gates at `:1858-1862` and `:1602-1609`. Refusing
publication on a surprise is the correct default — no new gate needed.

The message must identify the file, not just the seqset: seqset name, file class,
URL, and the configured prefixes.

**6. Reporting.** Add `records_excluded: int` to `SeqsetStats`
(`build_store.py:496-511`) and surface it in `print_summary` (`:1656`), following
the existing `rows_skipped_non_refseq` precedent (`:496`, used at `:1130`).

## Measured cost

Benchmarked on release 100's real `dna.toplevel` — 1.11 GB gz, 64.20 GB
uncompressed, 194 records kept and 445 dropped:

| variant | wall | over decompress | output |
|---|---:|---:|---:|
| decompress + parse only | 244.8s | — | — |
| filter → uncompressed | 243.2s | ~0 (noise) | 3.15 GB |
| **filter → gzip level 1** | **272.4s** | **+27.6s** | **1.02 GB** |
| filter → gzip level 6 | 619.3s | +374.5s | 0.88 GB |

Level 6 spends an extra 347s per file to save 140 MB over level 1 — 3.3 CPU-hours
across the set for 4.8 GB. Level 1 costs ~10% overhead for a 3.1x reduction.
Python's `gzip.open` defaults to level 9, so the default is worse still.

**The compressed file barely shrinks: 1.11 GB → 1.02 GB (−7.9%), while its
uncompressed content drops 64.20 → 3.15 GB (−95.1%).** Runs of `N` cost gzip
almost nothing. The padding is nearly free in the distributed file and expensive
only in the store, which uses fixed-width 3-bit encoding with no run-length
compression. That asymmetry is the entire reason this issue exists, and it is why
the saving appears in the store and not in the cache.

**This is not purely relocated work.** Today gtars decompresses those same
64.20 GB and digests them — a full old release (4 files at `jobs=4`) imports in
134–195s, and only 10 of 42 `dna/` directories have an `.rgsi` sidecar, so these
files really are read in full. After the change we add two costs: a 3.15 GB
compress plus a 3.15 GB decompress by gtars that did not exist before (~5% of the
bytes), and the 64.20 GB decompress moves from Rust running concurrently with
three sibling files to Python running in one thread.

Net: roughly **+2 hours one-time** across 34 releases, cut to ~20 minutes by
`--filter-jobs 8`. Thereafter the cached filtered file carries its own `.rgsi`,
so every subsequent build reads 3.15 GB instead of 64.20 GB and is *faster* than
today, and the 21 GiB store saving is permanent.

## Constraints

- **No git operations.** Do not push, and do not commit unless asked. The 9
  existing unpushed commits and the uncommitted working tree stay as they are.
- **Never delete from `downloads/`.** The cache is 69 GB of verified upstream
  artifacts and is the reason a rebuild is possible at all. Derived files written
  *into* the cache tree are additions, not replacements.
- **`store/` and `build.lock.json` may be modified or deleted.** Neither has been
  published; the published artifact is the separate `store.2026-07-22/`.

## Applying the change: incremental, not a full rebuild

A full rebuild is **not** required. gtars 0.9.2 exposes
`remove_collection(digest, remove_orphan_sequences=True)` — confirmed against the
installed package's own `.pyi`, not the local checkout. Sequences are stored one
file per digest at `store/sequences/<2-char>/<digest>.seq`, so orphan collection
reclaims real bytes.

Only the 34 `dna.toplevel` collections change; the 102 `cdna`/`ncrna`/`pep`
collections, all assemblies, LRG, and RefSeq are untouched, as is the rolling
`ensembl` namespace (it tracks release 116).

Per affected release, the old collection digest is already recorded in
`build.lock.json` keyed by source cache path:

0. Run `prepare_filtered_sources` once for all 34 releases, in parallel
   (`--filter-jobs`). Everything below then works from cached files.
1. Take that release's filtered `dna.toplevel` from the preflight output.
2. `remove_collection(old_digest, remove_orphan_sequences=True)`.
3. `add_sequence_collection_from_fasta(filtered_path)` → new digest.
4. Rewrite `store/aliases/sequences/ensembl-N.tsv` from that release's four
   collections — three unchanged, one new. **Removing a collection does not
   remove aliases**; `remove_sequence_alias` is a separate API, so stale `CHR_`
   aliases would otherwise dangle onto deleted digests. `build_store.py:1622`
   already writes a namespace TSV directly for the rolling namespace, so there is
   in-repo precedent for this.

Prefer two passes — remove all 34, then ingest all 34 — so orphan collection is
reasoned about once. A padded sequence shared across releases 76–109 only becomes
an orphan when the last referencing collection is gone.

### Three behaviours to verify empirically before trusting this

The local `/Users/kferrite/dev/gtars` checkout is **v0.9.0 / python 0.9.1**
(commit `99ca142`), behind the installed **0.9.2**, so its Rust source cannot be
used to settle semantics. Prototype on a small scratch store first:

- **Orphan GC under lazy loading.** Collections load as metadata stubs
  (`is_collection_loaded`). If orphan detection consults only *loaded*
  collections, it could delete sequences still referenced by a stub. Test:
  build a store with two collections sharing a sequence, remove one without
  loading the other, confirm the shared sequence survives. If it does not,
  `load_all_collections()` first — and measure the memory cost.
- **`.seq` files are actually unlinked** and `sequences.rgsi` / `collections.rgci`
  are rebuilt consistently by `write()`.
- **Alias TSV rewrite survives `write()`** rather than being clobbered by
  in-memory state.

Also measure the cost of `remove_orphan_sequences=True`: each call may scan all
remaining collections, so 34 calls could be O(34 × all collections).

### Fallback and safety net

Clone the store before mutating — `cp -Rc store store.pre-filter` is near-free on
APFS and makes a failed incremental attempt recoverable without re-ingesting
48.9 GiB. If any of the three behaviours above does not hold, fall back to a full
rebuild from the (untouched) download cache, roughly 2–3 hours.

## Verification

Nothing in `verify_store.py` hardcodes a count that shifts: the floors
(`:141-155`) are loose, `GROUND_TRUTH`/`ENSEMBL_GROUND_TRUTH` (`:39-49`,
`:61-71`) contain no `CHR_` entries and are pinned to `ensembl-113`, and the
known-divergence fixture has zero `CHR_` rows. The one structurally sensitive
check is namespace-set equality (`:156-162`) — safe here because chromosomes and
real scaffolds remain, so every `ensembl-N` namespace survives.

1. **Smoke test one release first**, into a scratch store: confirm 445 dropped,
   `ensembl-100:CHR_HSCHR1_2_CTG3` absent, `ensembl-100:1` present with its
   digest unchanged.
2. Apply to all 34, then `uv run pytest tests/ -q`, `uv run python
   verify_store.py`, and `gks-refgetstore verify --lock build.lock.json
   --strict-set`.
3. Regenerate `build.lock.json`. Upstream `sha256` values are unchanged; only the
   34 `dna.toplevel` `collection_digest` values move. The exclusion rule is
   captured because `build.sources_toml.sha256` is recorded in the lock.
4. Confirm the 445 padded digests are gone: each `store/sequences/<shard>/
   <digest>.seq` from `ensembl_padded_scaffolds.tsv` should no longer exist, and
   no alias in any namespace should resolve to one.
5. Re-run parity and confirm the predicted movement and nothing else: digest
   coverage **85.715% → 85.677%**, seqrepo-only **163,429 → 163,874**.
6. Confirm store size lands near **27.9 GiB**.
7. Diff against `store.pre-filter/` to confirm nothing outside the 34 collections
   and the `ensembl-76`–`ensembl-109` namespaces moved.

## Tests — `tests/test_source_resolution.py`

Match the file's conventions: hand-rolled `Store` doubles, `monkeypatch` of
`ensure_download` to `pytest.fail` to prove no network, real TOML in `tmp_path`,
full-sentence test names.

1. `exclude` validation cluster, extending the group at `:25-39` — unknown key,
   empty `record_prefixes`, empty `file_classes`, a `file_classes` entry not
   among the seqset's own file classes; plus a valid construction.
2. `filter_fasta_records` keeps and drops the right records, returns
   `(kept, dropped)`, and writes the audit log.
3. A declared exclusion that matches **zero** records raises, and the message
   names the offending file class and URL. This is the regression guard that
   replaces a pinned count, so assert on the message content.
4. Filtered output is cached until `force=True` — parallel to the existing
   derived-FASTA cache test.
5. Exclusion applies **only** to declared file classes: a two-file seqset where
   the undeclared second file also contains `CHR_`-named records leaves that file
   untouched.
6. Real-store round trip via `RefgetStore.in_memory()`, modelled on
   `:671-693` — `ensembl-N:CHR_x` absent, `ensembl-N:1` present.
7. A seqset with no `exclude` is byte-identical to today (no derived file
   written, same collection digest).
8. `prepare_filtered_sources` produces the same files at `jobs=1` and `jobs=4`,
   and `process_seqset` afterwards performs **no** filtering work — assert it
   reuses the preflight output rather than regenerating it (monkeypatch
   `filter_fasta_records` to `pytest.fail`, mirroring the established
   `ensure_download` no-network idiom).
9. `--filter-jobs` reaches `prepare_filtered_sources`. Extend the existing
   `Args` doubles in `tests/test_deferred_assembly_reports.py:39-58`, which
   already carry `ingest_jobs` and `min_free_gb`, so the new attribute does not
   break the six existing doubles.

## Follow-up issue

Write `ISSUE-lock-should-pin-ingest-outcome.md` in the repo root, matching the
style of `ISSUE-ensembl-padded-scaffold-bloat.md`. Do not implement it here.

**Problem.** The lock pins *upstream inputs*, not the *ingested outcome*.
`source_record` (`build_lock.py:129-196`) stores both `sha256` and
`collection_digest`, but the drift preflight (`evaluate_cache_vs_lock`,
`evaluate_sources_vs_lock`, `build_lock.py:282`, `:308`, driven from
`build_store.py:1806-1830`) compares only `sha256`. Anything that changes what
gets ingested from unchanged bytes is therefore invisible: a record-exclusion
rule edit, a filter-logic regression, a derived-FASTA converter change
(`gbff_to_fasta`, `lrg_zip_to_fasta`), or a gtars version that digests
differently.

**Why a record count is not the fix.** It cannot distinguish a compensating
change — one record fewer from one file, one more from another — and it needs
hand-maintained per-release numbers.

**Proposed direction.** Compare the realized `collection_digest` against the
lock's recorded value, since the lock already carries it. That is
content-addressed, so it pins the actual outcome and is immune to compensating
errors. Design questions the issue should raise rather than settle:

- Where to check — preflight cannot know the digest before ingest, so this is
  likely a post-ingest gate or a `verify` subcommand mode.
- How to express an intentional change. A config edit legitimately moves the
  digest; the check needs a strict/informational split like the existing
  `--lock-check-mode {strict,subset,ignore}` (`cli.py:59-63`).
- Whether the lock should also record the exclusion rule itself (prefixes and
  file classes), so a deliberate config change is attributable separately from a
  code regression. `build.sources_toml.sha256` currently detects *that* config
  changed but not *what* or *where*.
- First-build and new-source cases, which have no prior digest to compare.

## Docs

- `sources.toml` field comment block (`:1-53`) — add `exclude` in the existing
  two-space aligned style, plus a prose rationale comment above the Ensembl
  section banner (`:504-513`).
- `README.md` seqset field table (`:247-261`) for `exclude`, and the build-flag
  section for `--filter-jobs` alongside `--ingest-jobs` and `--min-free-gb`.
- `cli.py` `_add_build` (`:34-68`) — `--filter-jobs`, default
  `build_store.INGEST_JOBS_DEFAULT`, with help text explaining it is zlib-bound
  and independent of `--ingest-jobs`.
- `seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md` — record the outcome.
- `ISSUE-ensembl-padded-scaffold-bloat.md` — note the resolution, or delete it if
  it is superseded by the doc. Its option A described an N-fraction threshold;
  record that the implemented approach is name-based instead, and why.

## Repo state

Checked, no action taken. Remote `main` is `dc04e5e2`, an ancestor of our HEAD —
**nothing to pull**. We are 9 commits ahead, unpushed; `git push` fails on SSH
auth (`ssh-add -l` reports no identities), so pushing needs your key loaded.
The working tree has 6 uncommitted paths (the padding report, two TSVs,
`tools/ensembl_padding_probe.py`, the README edit, the issue draft).

Per the constraints above, this work does **not** commit or push anything. The
uncommitted paths stay uncommitted unless you ask otherwise. Note the consequence:
`build.lock.json` records `build.git.dirty`, so a lock regenerated now will be
marked dirty — accurate, and consistent with how the 2026-08-26 build recorded
`a755a69-dirty`.

## Risks

- **Incremental mutation is the main risk.** Three gtars behaviours are assumed
  and must be proven on a scratch store first (above). The `store.pre-filter/`
  clone plus the untouched download cache mean the worst case is a full rebuild,
  not data loss.
- **Disk.** ~35 GB of gzipped derived files (1.02 GB × 34), plus ~49 GB briefly
  for the pre-filter clone (APFS clone: near-zero until blocks diverge). Free
  space is 337 GB. Nothing is removed from `downloads/`.
- **Time.** ~2 hours one-time, dominated by decompressing 64.20 GB per release —
  not by compression, which is ~28s per file at level 1. `--filter-jobs 8` cuts
  the wall time to roughly 20 minutes. Steady-state builds get faster afterwards.
- **Immutable-namespace semantics.** `ensembl-76`–`ensembl-109` change content
  under unchanged names. That is intended, but anything pinned to those
  namespaces sees scaffolds disappear. The store has not been published with
  these namespaces, so no released artifact is invalidated.
- **gtars version skew.** The local checkout is behind the installed 0.9.2; treat
  it as background reading only. Behaviour must be confirmed against the
  installed package.
