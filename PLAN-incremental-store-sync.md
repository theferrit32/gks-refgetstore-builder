# `gks-refgetstore sync` — applying manifest changes to a built store

## Goal

Apply any `sources.toml` change to an existing store as incremental adds and
deletes: sources that did not change keep their lock validation and are not
re-ingested, sources that did change are removed and re-ingested, and the lock
is updated only where the store actually moved.

The Ensembl `CHR_` exclusion committed in `4dfc856` is the first real
application, not a special case. Its classification — 34 sources, all
`reingest` — must fall out of the general rules rather than being hardcoded.

## Why this rather than a one-off script

A purpose-built script would hardcode 8 collection digests and 34 release
numbers pasted out of an analysis. It would work once and be wrong the next time
the manifest changes. `sync` derives the same set from `sources.toml` and the
lock, so it stays correct.

The objection to doing this first — debugging new code during a long destructive
mutation — is answered by `store.pre-filter/`: run `sync` against the clone
first, verify, then against `store/`. A second APFS clone costs 7 minutes and
~1 GiB.

## Measured gtars 0.9.2 behaviour

All against the **installed** 0.9.2. The local `gtars` checkout is `99ca142`
(v0.9.0/0.9.1) and was not consulted.

**1. Orphan collection silently no-ops on a lazily loaded store.** After
`open_local`, `n_collections_loaded` is 0 and
`remove_collection(d, remove_orphan_sequences=True)` removes the collection but
collects **nothing** — a true orphan survived, `n_sequences` unchanged, no
error. After `load_all_collections()` it is correct: the orphan goes, a sequence
still referenced by a remaining collection stays. The lazy path under-collects;
it does not over-delete. This is the highest-value finding — the operation looks
like it worked and reclaims zero bytes.

**2. `load_all_collections()` is cheap.** 14.2 s, RSS 3.68 → 5.33 GiB, on 240
collections / 1,779,497 sequences. Metadata only; `n_sequences_loaded` stays 0.

**3. `write()` is incremental.** Retained `.seq` files kept inode and mtime
across a remove/re-add cycle; only `sequences.rgsi` (447 MB) and
`collections.rgci` were rewritten, and removed sequences were unlinked.

**3b. `remove_collection` persists immediately — there is no in-memory staging.**
A store opened with `open_local` is in persisting mode, and removal unlinks the
orphaned `.seq` files *and* rewrites `sequences.rgsi` / `collections.rgci`
without waiting for an explicit `write()`.

Discovered the hard way: a probe run against the `store.pre-filter/` clone that
deliberately never called `write()` nonetheless left it at 1,779,027 `.seq`
files and 232 collections — the post-removal state, on disk, permanently. The
result is internally consistent, not corrupt, but the clone is no longer a
pre-filter backup.

Two consequences for the design:

- **`--dry-run` must never call `remove_collection`.** Classification and
  prediction only. There is no "remove, inspect, abort" path.
- **The safety clone must be made immediately before the mutation**, and treated
  as spent the moment any removal runs against it. Call
  `disable_persistence()` first when a probe genuinely needs to be
  non-destructive.

**4. Aliases dangle.** Removing a collection leaves aliases listed *and still
resolving*. `remove_sequence_alias` returns `True`, survives `write()` +
reopen, and rewrites the namespace TSV.

**5. Orphan GC costs ~40–48 s per call** regardless of how much it frees — it
scans the remaining collection graph. Eight calls is ~6 min. This makes a sync
that removes many sources O(n²); see *Open questions*.

**6. Ingest cost is per-record, not per-byte — it is dominated by creating one
file per sequence.** Ablated on release-100 cdna (70 MB gz → 376 MB, 190,522
records):

| variant | time | `.seq` written |
|---|---:|---:|
| fresh, full | 28.39 s | 178,516 |
| fresh, `disable_persistence` | **1.46 s** | 0 |
| fresh, `disable_encoding` | 29.99 s | 178,516 |
| fresh, `disable_ancillary_digests` | 31.41 s | 178,516 |
| raw gzip decompress only | 0.23 s | — |

Persistence is **26.93 s of 28.39 s (95%)** — 178,516 file creations at ~151 µs
each. Decompression is negligible; encoding and ancillary digests are free (both
ablations landed slower than baseline — noise, plus larger raw bytes to write
for `disable_encoding`).

**7. Re-ingesting a collection the store already holds is ~8× cheaper**, because
those sequences are already on disk and are not rewritten:

| | time |
|---|---:|
| fresh | 27.95 s |
| re-ingest, `.rgsi` sidecar present | 3.16 s |
| re-ingest, sidecar deleted | 5.80 s |

The sidecar is worth 2.6 s — about 11% of the saving, not the mechanism. It is
written on first ingest only, and is not recreated if deleted.

> This is what makes seqset-granularity sync affordable, and removes the need
> for a bespoke per-source ingest path.

**8. The two source shapes have opposite cost profiles.** Many small records are
file-creation bound; few huge records are decompression bound. Release-100
`dna.toplevel` is ~639 records / 64.20 GB uncompressed, so persistence is
irrelevant and decompression (~245 s) dominates.

Filtering inverts that source completely — 64.20 GB → **2.93 GiB** (−95.1%),
because the padding *was* the content:

| release-100 `dna.toplevel` | uncompressed | decompress | fresh ingest |
|---|---:|---:|---:|
| original | 64.20 GB | ~245 s | — |
| filtered | 2.93 GiB | 4.2 s | **11.5 s** (194 records) |

So ingesting the filtered files is nearly free. The dominant cost of the whole
operation is the **filtering preflight**, which must still decompress the 34
originals once.

**7. `list_collections()` returns a paginated dict**, not a list —
`{"results": [...], "pagination": {"page": 0, "page_size": 100, "total": N}}`.
With 240 collections it must be paged.

## Removal order is required, not preferred

Removing all 8 affected collections in memory:

```
before: n_sequences=1779497  n_collections=240
  zF7K0k1BaCi7  40.34s  n_seq=1779473     <- releases 76-78
  ... five more, all 1779473 ...
  qJ79liNTAD-L  45.58s  n_seq=1779027     <- releases 98-109, last reference
after:  n_sequences=1779027  n_collections=232
```

The first seven freed **24 sequences between them**; the eighth freed **446**. A
sequence published by several release groups only becomes an orphan when the
last collection referencing it is gone. **All removals must precede all
ingests** — interleaving per source reclaims almost nothing.

### The counts are a mid-flight checkpoint

Removing all 8 frees **470**, not 445. The extra 25 are non-`CHR_` records
existing only in the pre-110 toplevel collections; they are not padded, not
matched by the exclusion, and are restored by re-ingest.

| point | `n_sequences` | delta |
|---|---:|---|
| before | 1,779,497 | |
| after removing all 8 | 1,779,027 | −470 |
| after re-ingest | **1,779,052** | **−445 net** |

−445 is exactly the count of distinct `CHR_` names, reached independently of the
padding analysis. `sync` should predict both numbers and abort on disagreement.

## Change 1 — pin the ingest spec in the lock (schema `/3`)

This is the prerequisite. Today the preflight compares only upstream `sha256`,
so a `sources.toml` exclusion edit changes nothing it looks at and those 34
sources would classify as `unchanged` and be skipped — silently wrong.

`build_lock.source_record` gains two fields:

```json
"ingest_spec": {"format": "fasta",
                "exclude": {"record_prefixes": ["CHR_"]}},
"ingest_spec_sha256": "<sha256 of canonical JSON>"
```

Canonical form is `json.dumps(spec, sort_keys=True, separators=(",", ":"))`.

**The exclusion is recorded only when it applies to *this source's* file class**
— `entry.exclusion.applies_to(source.file_class)`. Otherwise editing the
`dna.toplevel` rule would mark every `cdna`/`ncrna`/`pep` source in the same
seqset as changed. `file_class_for_index` already provides the pairing.

**Migration falls out of the rules.** The existing lock is schema `/2` and has
no `ingest_spec_sha256`. Treat absence as *unknown*, and resolve conservatively:

- no spec in lock **and** no spec desired → `unchanged`
- no spec in lock **and** a spec desired → `reingest`

For the current lock that classifies the 34 `dna.toplevel` sources as `reingest`
and all 301 others as `unchanged`, with no special-casing anywhere.

## Change 2 — classification, as a pure function

New in `build_lock.py`, no store access, exhaustively unit-testable:

```python
@dataclass
class SyncPlan:
    unchanged: list[str]                    # cache_path
    added: list[str]
    removed: list[dict]                     # full lock records; need collection_digest
    reingest: list[tuple[str, str]]         # (cache_path, reason)

def plan_sync(sources, lock, download_dir, spec_by_url) -> SyncPlan
```

| class | condition | action |
|---|---|---|
| `unchanged` | in lock, `sha256` matches, spec matches | validate only |
| `added` | not in lock | download + ingest |
| `removed` | in lock, absent from `sources.toml` | remove collection + aliases |
| `reingest` | `sha256` **or** `ingest_spec_sha256` moved | remove, then ingest |

`evaluate_cache_vs_lock` and `evaluate_sources_vs_lock` already compute most of
the inputs; `plan_sync` composes them and adds the spec dimension.

## Change 3 — store mutation operations

New `store_sync.py`:

```python
def load_for_mutation(store) -> None:
    """load_all_collections(). Required before any removal: gtars 0.9.2
    silently skips orphan collection on a lazily loaded store."""

def remove_collections(store, digests) -> RemovalReport:
    """Remove distinct collection digests, orphan GC on. Asserts
    load_for_mutation ran -- the silent no-op is the failure to design against."""

def reconcile_namespace(store, store_dir, namespace, desired) -> tuple[int, int]:
    """Make the namespace hold exactly `desired`; returns (added, removed)."""
```

`reconcile_namespace` is the piece with no existing equivalent.
`load_immutable_aliases` is add-only and raises on a changed digest, which is
the correct *build-time* invariant but is exactly what a reconcile must be
permitted to do.

### It must write the namespace TSV directly, not call the alias APIs

Two measurements decide this.

**`load_sequence_aliases` merges; it does not replace.** Loading a 2-row TSV
into a namespace holding 3 rows left all 3 in place. So it cannot be used to
shrink a namespace, and the "replaced complete namespace" log at
`build_store.py:1857` is a misnomer — that path works because `_write_alias_tsv`
writes the file *directly* to the store's alias directory, and
`load_sequence_aliases` is only called the first time, to register the
namespace.

**`remove_sequence_alias` is O(namespace) per call**, at 16.81 ms each on a
3,000-alias namespace — it rewrites the whole TSV every time, the same
quadratic already documented for `add_sequence_alias` at `build_store.py:1583`.
The `ensembl-N` namespaces hold ~360,000 aliases, ~120× larger, so roughly 2 s
per removal. The `CHR_` reconciliation is ~10,000 removals across 34
namespaces: **on the order of 5 hours.** Not viable.

**A hand-written namespace TSV is durable.** Overwriting
`store/aliases/sequences/<ns>.tsv` and reopening showed the new contents, and a
subsequent `store.write()` did **not** clobber it — the file still held exactly
the written rows after write + reopen.

So `reconcile_namespace` computes the complete desired map and writes it with
`_write_alias_tsv` — one atomic O(namespace) write instead of O(namespace ×
removals). This reuses the existing rolling-namespace mechanism rather than
inventing one.

Ordering constraint that follows: alias TSVs are written **after**
`store.write()`, and the store is not written again afterwards, since an open
store's in-memory alias map does not know about the hand-written file. Sync
reopens read-only at the end to verify.

## Change 4 — reuse the existing ingest path

`sync` operates at **seqset scope**. Any seqset containing an `added`,
`removed`, or `reingest` source is reprocessed in full through the existing
`process_seqset` / `process_release_groups`, feeding `alias_sink`. Finding 6 is
what makes this affordable: the unchanged `cdna`/`ncrna`/`pep` shards re-ingest
at ~3 s each rather than ~26 s.

Seqset scope also matches `merge_into_lock`, which already keys on
`(kind, owner)`.

Release namespaces span several seqsets (`dna.toplevel` + `cdna` + `ncrna` +
`pep` all feed `ensembl-100`), so reconciliation happens per *namespace* after
every contributing seqset is processed — the grouping `process_release_groups`
already performs. The only change there is replacing the add-only
`load_immutable_aliases` call with `reconcile_namespace` when running under
`sync`.

## Change 5 — CLI

```
gks-refgetstore sync [--config sources.toml] [--store store] [--lock build.lock.json]
                     [--dry-run | --apply] [--seqset NAME]
                     [--filter-jobs N] [--ingest-jobs N]
```

`--dry-run` is the default; `--apply` is required to mutate. Dry run prints the
classification with a reason per source and the predicted `n_sequences` delta,
and exits non-zero if the plan is empty when the user expected work.

`merge_into_lock` needs a `dropped_scopes` argument so `removed` sources
disappear from the lock rather than being carried forward.

## Execution order

1. Classify. No store access.
2. Print the plan; stop unless `--apply`.
3. Preflight downloads and derived FASTAs for `added` + `reingest`.
4. Open store, `load_for_mutation`.
5. **Remove every affected collection**, before any ingest.
6. Check `n_sequences` against the predicted post-removal figure; abort on
   mismatch.
7. Ingest the affected seqsets.
8. `write()`.
9. Reconcile the affected namespaces by writing their TSVs directly — after
   `write()`, and with no `write()` afterwards.
10. Reopen read-only; check `n_sequences` against the predicted final figure and
    the namespace set against the pre-sync set. Abort on mismatch.
11. `merge_into_lock`, write lock.

## Tests

Following `tests/test_source_resolution.py` conventions — real TOML in
`tmp_path`, hand-rolled doubles, `monkeypatch` of `ensure_download` to
`pytest.fail` to prove no network.

1. `plan_sync` classification, table-driven, one case per rule.
2. Both schema-`/2` migration cases: absent spec + none desired → `unchanged`;
   absent spec + spec desired → `reingest`.
3. `ingest_spec_for` omits an exclusion that does not apply to that file class,
   so editing the `dna.toplevel` rule leaves `cdna` sources `unchanged`.
4. **Removal without `load_for_mutation` raises** rather than silently
   no-opping. This encodes finding 1 so it cannot regress.
5. `reconcile_namespace` adds missing, removes surplus, leaves correct rows —
   and asserts it did **not** go through `remove_sequence_alias`, which is the
   O(namespace) trap.
6. Integration: 3-collection store in `tmp_path`, add an exclusion, `sync
   --apply`, assert the collection was replaced, the orphan is gone, aliases
   reconciled, and the lock changed only in the touched scope.
7. Removal ordering: a sequence shared by two collections survives removing the
   first and is collected on the second.

## Applying it to `store/`

1. Second APFS clone as the working target, or run against `store.pre-filter/`
   first. ~7 min, ~1 GiB real.
2. `sync --dry-run` — expect exactly 34 `reingest`, 0 `added`, 0 `removed`, 301
   `unchanged`.
3. `sync --apply` against the clone. Verify.
4. `sync --apply` against `store/`.

Estimated **~35–40 min**, from the measured per-step costs:

| step | estimate | basis |
|---|---:|---|
| filter 34 toplevel files | **~20 min** | 34 × ~273 s at `--filter-jobs 8` |
| `load_all_collections` | 14 s | measured |
| remove 8 collections | ~6 min | 8 × ~44 s measured |
| ingest 34 filtered toplevel | ~4 min | 8 fresh × 11.5 s + 26 re-ingest |
| ingest ~102 cdna/ncrna/pep | ~5 min | ~3 s each, short-circuited |
| `write()` + 34 alias TSVs | ~1 min | 447 MB index rewrite |

Filtering is ~55% of the run. It cannot be skipped for the 26 non-representative
releases: sync reprocesses at seqset scope, so every affected seqset needs its
filtered file present.

**Rejected optimization.** Only 8 filtered files are strictly needed for the
*store*, since the 34 releases collapse to 8 collection digests. Filtering just
those and assigning each group's digest to its member releases would save ~15
min — but it makes the outcome depend on an assumed content identity rather than
a measured one, which is precisely the hardcoding this design exists to avoid.
Not worth 15 minutes.

### Acceptance

- `n_sequences` = **1,779,052** (−445); `n_collections` returns to 240 after
  dipping to 232.
- Store ≈ **27.9 GiB**, from 48.9.
- Every `.seq` path in `ensembl_padded_scaffolds.tsv` is gone.
- No alias resolves to a padded digest; the namespace **set** is unchanged — all
  34 `ensembl-N` namespaces survive, since chromosomes and real scaffolds remain.
- `ensembl-100:1`, `X`, `MT` digests unchanged; `ensembl-100:CHR_HSCHR1_2_CTG3`
  absent.
- `pytest tests/ -q`, `verify_store.py`, and `verify --lock --strict-set` pass.
- Lock diff touches only the 34 `dna.toplevel` `collection_digest` /
  `n_sequences` values plus `build.*` and the new spec fields.
- Parity: coverage **85.715% → 85.677%**, seqrepo-only **163,429 → 163,874**.

## Open questions

- **Orphan GC is O(n) per call.** Fine at 8. A sync removing hundreds of sources
  would want `remove_orphan_sequences=False` on all but the last call, but gtars
  0.9.2 exposes no standalone orphan-collection entry point. Worth raising
  upstream — it is a natural API gap, and the workaround (carry the flag only on
  the final removal) depends on undocumented ordering behaviour.
- **First builds**, where every source is `added` and there is no prior lock:
  `sync` should degrade to the existing build path rather than special-case.

## Constraints

- Nothing is deleted from `downloads/`; derived files there are additions.
- `store/` and `build.lock.json` may be modified. Neither is published — the
  published artifact is `store.2026-07-22/`.
- `store.pre-filter/` is the rollback; worst case is a rebuild from the
  untouched cache (~2–3 h), not data loss.
- No pushing. Commits only when asked.
