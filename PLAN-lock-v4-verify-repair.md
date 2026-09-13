# Build lock v4, verify/status/repair, and a fast dev manifest

> **Status: implemented.** Retained as the design record — the reasoning behind
> the shape, and the measurements the decisions rest on. Behaviour is documented
> in `README.md`; this file is why it is that way.

## Context

`build.lock.json` is supposed to serve two roles:

1. **Upstream reproducibility** — re-pulling the sources gives the same bytes.
2. **Local completeness** — the download cache *and the built store* still
   contain everything they should.

Role 2 is entirely unimplemented for the store. `RefgetStore` appears in
`build_lock.py` only inside `run_lock`; **neither verify path has ever opened
the store.** The lock records `collection_digest`, `n_sequences`, and
`build.store.*`, and nothing has ever checked them. A deleted `.seq` file is
invisible.

`verify` also calls `resolve_sources`, which fetches provider manifests, so
offline verification is impossible today and results depend on live upstream
state.

This plan closes both gaps, restructures the lock around `inputs`/`outputs`,
adds `status` and `repair`, and introduces a small dev manifest so the loop is
minutes rather than hours.

The project is unpublished and unused elsewhere. **No backwards compatibility.**
Old schema versions are deleted outright, not shimmed.

## The finding that reorders priorities

Building `verify --deep` as a design exercise turned up a real data-correctness
bug, confirmed directly against the store and the source FASTAs:

**gtars' `protein` alphabet has no `U` (selenocysteine). Every selenoprotein is
silently corrupted at ingest.**

```
ENSP00000473614  GPX4 (glutathione peroxidase 4), 180 aa
  source charset   ACDEFGHIKLMNPQRSTUVWY
  decoded charset  ACDEFGHIKLMNPQRSTVWY      <- no U
  position 109     source U  ->  stored A
```

Six of six `U`-containing records sampled from release-113 `pep` are present in
the store **under the correct digest** — the digest was computed from the
original bytes — but decode without the `U`. So:

- the digest is right,
- the stored payload is wrong,
- the store cannot serve these sequences correctly,
- **re-ingesting reproduces it**, because it is an encoder limitation.

Full scan of `store.pre-filter` (1,779,052 sequences):

| alphabet | sequences | round-trip failures |
|---|---:|---:|
| dna2bit | 1,426,421 | 0 |
| dna3bit | 643 | 0 |
| **protein** | 351,881 | **106** |
| **dnaio** | 107 | **27** |

The 27 `dnaio` cases are a second, distinct bug: short proteins whose residues
happen to all be legal IUPAC nucleotide codes get misdetected as nucleotide and
encoded lossily (`ENSP00000499040.1`, 12 aa).

This is not bit rot and `repair` must not attempt it.

**Decision: leave it in the store and fix it upstream in gtars later.** The
encoder is not ours, no amount of re-ingesting helps, and the affected sequences
are a 0.0075% tail. What this work owes it is *visibility*: these 133 become the
worked example proving `verify --deep` surfaces real data defects rather than
theoretical ones. See "Flagging the encoder mismatch" below.

## Decisions already settled

- **Single lock file, reorganized into `inputs` / `outputs`** (not split, not a
  flat `sources` list).
- **Full scope**: `verify` with four targets, `status`, and `repair`.
- **`status` always exits 0.** It describes; `verify` judges and gates CI.
- **Do not refresh the drifted NCBI files.** They are live test material.

## Hazard to guard first

`ensure_download` re-downloads whenever the cached file's *provider* checksum
mismatches. Upstream has since changed **all 30** remaining `mRNA_Prot` files'
MD5s (not just the two withdrawn shards). So running `build` or `fetch` against
`sources.toml` today would overwrite exactly the bytes `build.lock.json` pins
and destroy the preserved drift.

**Fix as part of this work**: when a lock covers the file and the cached bytes
match the lock's `sha256`, keep the cache and report provider drift rather than
silently re-downloading. `--force-download` remains the override. For a
lock-based tool the lock must win.

---

## Part 1 — `sources.dev.toml`

19 files / **3.19 GB** / **~3 minutes**, versus 335 files / 106 GB / hours.
Every source is already cached and still matches its provider checksum.

| entry | shape covered |
|---|---|
| `GRCh37` | assembly, `load_fasta=true` (smallest, 872 MB) |
| `GRCh37.p13` | **`fasta_path` local override** — zero production usage |
| `GRCh37.p11`, `GRCh37.p12` | report-only assemblies; deferred-report path |
| `refseq_human_refseqgene_dev` | **`url_template` + `shard_range`** — zero production usage; mutable, no checksum manifest, so no download risk |
| `lrg_public` | `format="lrg_zip"`, derived FASTA, no manifest anywhere |
| `refseq_history_grch37_protein_dev` | `md5_manifest_urls` + `file_classes` + **cold-start `exclude`** |
| `ensembl_release_76` | Ensembl bsd-sum, `release`/`rolling_namespace`, production `CHR_` exclude (filtered file already cached) |
| `ensembl_release_113` | Ensembl without exclude; the only mutable release; makes the rolling namespace **advance** |

**The exclusion prefix is verified, not guessed.** Across the three GRCh37
protein files: 88,782 `NP_`, 5,761 `XP_`, 26 `YP_`. `XP_` matches a real subset
of **every** file. `YP_` would be a trap — absent from the first file, and
`filter_fasta_records` raises on a rule that drops nothing.

**Deliberately excluded: shape C** (`url_pattern` + `checksum_manifest_url`).
Building it would re-download 39 drifted files and destroy the test material. It
stays covered by unit tests with a stubbed `fetch_text`, and read-only against
the production config.

**The dev cache must be `downloads/` itself** — `fasta_path` is repo-root
relative, derived artifacts are written as siblings of their source, and
re-deriving release 76's filtered toplevel costs 273s. A separate cache dir is
not viable.

Artifacts: `store.dev/` (already matched by `.gitignore`'s `/store*/`),
`build.dev.lock.json` (add to `.gitignore`; `sources.dev.toml` is tracked).

## Part 2 — lock schema `/4`

```json
{
  "schema": "gks-refgetstore-build-lock/4",
  "build": {"timestamp_utc":…, "git":{…}, "gtars_version":"0.9.2"},
  "inputs": {
    "sources_toml_sha256": "…",
    "files": [{"kind","owner","url","cache_path","mutable","sha256","bytes",
               "present","provider_checksum","provider_checksum_algorithm",
               "provider_checksum_blocks","checksum_url","file_class",
               "ingest_spec","ingest_spec_sha256"}]
  },
  "outputs": {
    "n_sequences": 1779052, "n_collections": 235,
    "sequences_root": "sha256:…", "collections_root": "sha256:…",
    "collections": [{"digest":…, "n_sequences":…, "from":[cache_path, …]}]
  }
}
```

**`outputs.collections[].from` is a list**, which corrects a false premise. The
old `sources[].collection_digest` assumed one file → one collection; filtering
collapsed 8 Ensembl collections into 3, and 43 digests in the current lock
already have ≥2 contributors (one has 12). `provenance.py`'s docstring asserts
the false version.

**`outputs` is a census enumerated from the store** at write time, not
synthesized from provenance. Provenance only supplies the `from` annotation, so
a collection with no known contributor is recorded honestly as `"from": []`.

`upstream_md5` is dropped — losslessly recoverable from `provider_checksum` +
`provider_checksum_algorithm`. `build.store.*` counts were being written as
strings by `store.stats()`; the census uses `len()` and the bug disappears.

**Roots** (canonical, documented so a reimplementation can't drift):

```python
def digest_root(digests):
    h = hashlib.sha256()
    for d in sorted(set(digests)):     # bare digests, no SQ. prefix
        h.update(d.encode("ascii")); h.update(b"\n")
    return "sha256:" + h.hexdigest()
```

Measured on the real store: enumerate 2.5s + sort/hash 1.0s = **~4s added to
every lock write**. Verified value for `store.pre-filter`:
`sha256:3ecf0d7221a60fa038aa12baac6b1da646b4917513b176a1d0661b1541ed921c`.

New `store_census.py` holds the read-only primitives shared by the lock writer,
verify, repair, provenance, and `tools/compare_store_mutation.py`:
`collection_census`, `sequence_digests`, `collection_members`,
`on_disk_sequence_digests`, `digest_root`, `redigest_sequence`, `seq_path`.

Use `stream_sequence` for re-digesting, never `get_sequence()` — O(1) memory
instead of loading 250 Mbp into RSS.

## Part 3 — `verify`

```
gks-refgetstore verify [--cache] [--store] [--remote] [--manifest] [--all] [--deep]
                       [--no-hash] [--jobs N] [--limit N] [--known-bad TSV]
```

Default with no scope flags: `--cache --store`. **Offline and lock-driven.**
Only `--remote` and `--manifest` touch the network. Drop the no-lock integrity
mode — it resolves the manifest over the network, which is the bug being fixed.

Findings carry a stable `code`, which is the contract `repair` dispatches on:
`cache.missing`, `cache.sha256_mismatch`, `cache.gzip_corrupt`,
`store.sequences_root_mismatch`, `store.collection_missing`,
`store.seq_file_missing`, `store.orphan_sequence`, `store.seq_digest_mismatch`,
`manifest.file_added`, `manifest.file_removed`, etc.

### The store ladder — measured, not estimated

| level | proves | does not prove | cost |
|---|---|---|---|
| **L0** roots + counts | the store holds **exactly** the digest sets the lock recorded; any single add/remove flips a root | payloads exist or are readable | **~7s** |
| **L1** membership | every locked collection is present, fully resolvable, backed by `.seq` files; nothing stranded | any byte is correct | **~85s** |
| **L2** `--deep` | the store can actually serve each sequence and return what its digest promises | nothing about the lock — needs no lock | **~18 min** |

L1 uses one `scandir` pass (5.8s, 1,779,052 files) rather than 21M `exists()`
calls. L1 confirmed **0 orphans** in the current store.

L2 is ~18 minutes single-threaded (1,600 seq/s), not the hours I first
estimated. **Do not add `--jobs` to L2 without measuring** whether gtars
releases the GIL in `stream_sequence`.

`--limit` is accepted for `--cache`/`--remote` and **rejected for `--store`** —
a partial root is meaningless.

### Flagging the encoder mismatch

The 133 known defects are not suppressed. They are pinned as a **baseline**, so
the CLI reports them without drowning a real regression, and so a future gtars
fix is detected rather than missed.

`--known-bad` (default `seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv`)
is a checked-in TSV of `digest, name, alphabet, length, redigest, cause`,
following the existing `known_divergence/` convention. Three-way semantics:

| digest in baseline | re-digests correctly | outcome |
|---|---|---|
| yes | no | `warn` — known defect, rendered with its diagnosed cause |
| **no** | **no** | **`error`** — a new defect, which is the regression we care about |
| **yes** | **yes** | `info` — **gtars fixed it**; update the baseline |

That third row is the reason a baseline beats a suppression flag: it turns the
upstream fix into something the tool tells us about.

Findings render with enough to act on, not just a digest pair:

```
$ gks-refgetstore verify --store --deep

[warn ] store.seq_digest_mismatch  z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV
        ENSP00000473614   protein   180 aa
        decoded bytes re-digest to 4otSS66T9mEsIeu_jmydFFDB_x-LFqsb
        known cause: gtars protein alphabet has no U (selenocysteine);
                     residue replaced on encode. Upstream defect, not bit rot.

store   1,779,052 sequences re-digested in 17m42s
        1,778,919 ok
              106 known encoder defects (protein, no U)      [baseline]
               27 known encoder defects (alphabet misdetect) [baseline]
                0 new
PASS (133 known, 0 new)
```

Generate the initial baseline by running `--deep` once and capturing the
findings into a `runs/<date>-deep-verify/` record — the same way every other
ground-truth fixture in this repo was produced. The `cause` column is
hand-curated from this investigation; verify cannot infer it, because at verify
time an encoder defect and bit rot are indistinguishable.

File the gtars issue with the two reproducers: `ENSP00000473614` (`U` → `A`) and
`ENSP00000499040.1` (12 aa protein encoded as IUPAC nucleotide).

Make `remote.*` findings severity `warn`: `Content-Length` is missing on several
endpoints and mutable endpoints legitimately change size. `--manifest` is the
better network signal for the same cost.

## Part 4 — `status`

Four legs, one implementation each, shared with `verify` (no second copy of
"resolve and diff"):

| leg | mechanism | network |
|---|---|---|
| manifest ↔ lock | `sha256_file(config)` vs `inputs.sources_toml_sha256` | no |
| manifest ↔ cache | `resolve_sources` + `plan_sync` | yes (`--offline` skips) |
| lock ↔ cache | `evaluate_cache_vs_lock` | no |
| lock ↔ store | `verify_store(level=0)` | no |

Always exits 0.

## Part 5 — `repair`

```
gks-refgetstore repair [--apply] [--cache] [--store]
```

- Cache findings → unlink and re-fetch, **accepting only bytes matching the
  locked `sha256`**. Repair restores the locked state; it never re-baselines.
  New upstream bytes are a `sync` / `--force-lock` decision.
- Store findings → map digest → collection → `from` → owner scopes, then hand a
  synthesized `SyncPlan` to the existing `store_sync.apply_plan`, which already
  sequences prepare-filter → `load_for_mutation` → remove → re-ingest → write →
  reconcile in the only safe order.
- `store.seq_digest_mismatch` → **unrepairable**; re-ingesting reproduces the
  encoder defect.
- Cache repair always precedes store repair, and aborts before touching the
  store if it fails — removal persists immediately and has no rollback.
- A root mismatch with no attributable cause → refuse and point at `sync`.

Extract `fetch_sources.fetch_many(...)` so `run_fetch` and `repair` share one
downloader.

## Cross-cutting bugs to fix in this work

1. **`ensure_download` clobbers lock-pinned bytes** on provider-checksum drift
   (above). Highest priority — it is destructive and currently live.
2. **`prepare_filtered_sources` filters the raw cache path** while
   `process_seqset` filters the *derived* FASTA. A seqset combining
   `format="lrg_zip"` with `exclude` hands the filter a ZIP, drops nothing, and
   raises. Latent — production never combines them.
3. **Re-ingest scope is narrower than removal scope.** A collection is removed
   because one contributor changed, but carries every contributor's
   contribution. `apply_plan` expands by *namespace*, which accidentally covers
   Ensembl and misses cross-namespace sharing (e.g. a RefSeq file appearing
   under both `GCF/000/001/405/…` and `annotation_releases/9606/…`). The `from`
   list makes this fixable: add `reingest_scope(lock, rels)` and a regression
   test for the cross-namespace case.
4. **`fetch_sources.py:161`** — unguarded `locked_by_url[url]` raises `KeyError`
   mid-loop after some files are already fetched.
5. **`apply_locked_sources` builds `ResolvedSource` positionally** from lock
   keys. Move it into `build_lock.py` (which already imports `ResolvedSource`,
   so no circular import) and use keyword construction. This also removes the
   **duplicated schema allow-list** in `build_store.LOCKED_SOURCE_SCHEMAS`,
   leaving exactly one `SCHEMA` constant.

## Migration surface

Every `lock["sources"]` reader goes through new accessors (`lock_files`,
`lock_collections`, `collection_by_file`, `files_by_collection`) — never raw
subscripting. Sites: `build_lock` (5), `build_store.run_build`,
`fetch_sources`, `store_sync`, plus two readers that currently **bypass
`load_lock` with raw `json.loads`** and would silently produce wrong output:
`seqrepo_equivalence/full_parity.py:76` and
`runs/2026-07-23-…/compare.py:53`.

`provenance.py` is redesigned around `outputs.collections[].from`, and its false
one-file-one-collection premise deleted. `--digest` stays an O(all collections)
scan — 56s measured, fine interactively; do not build a cached reverse index.

**The committed lock is converted once, not rebuilt.** A rebuild would destroy
the preserved drift. Write a throwaway `convert.py` inside a
`runs/<date>-lock-schema-v4/` record that reshapes the `/2` lock and enumerates
`outputs` from `store/` (1,779,497 seqs / 240 collections — `store/` has *not*
had the padding sync applied). **`ingest_spec` must be written `null` for every
file**: that is the truth about `store/`, which was built before the exclusion
rule existed.

`verify_store.py` stays separate — it is semantics-driven (seqrepo ground truth,
assembly reports), not lock-driven. Only delete its two loose floors
(`n_collections >= 4`, `n_sequences >= 900`), which L0 subsumes exactly.

## Order of work

1. `store_census.py` + tests; repoint `tools/compare_store_mutation.py`. Additive.
2. Schema `/4` + accessors + `validate_lock`; move `apply_locked_sources`; fix
   the `ensure_download` clobber and the `fetch_sources` `KeyError`.
3. Readers: `provenance.py`, `full_parity.py`, the run-record `compare.py`, plus
   the one-shot lock conversion.
4. `sources.dev.toml` — land it here so everything after is testable in minutes.
5. `verify.py` + CLI, including the known-bad baseline; capture the 133 failures
   into `known_divergence/gtars_encoding_roundtrip.tsv` via a
   `runs/<date>-deep-verify/` record; file the gtars issue.
6. `status` (pure composition).
7. `repair` + the `reingest_scope` fix.

Docs land with the behaviour change, not in a trailing commit: `README.md`
(verify and build-lock sections), `RUNBOOK.md`. Mark
`ISSUE-lock-should-pin-ingest-outcome.md` closed, answering its four questions.

## Verification

**Unit** (all offline; `tests/conftest.py` already blocks `urlopen`):
`test_lock_schema.py` (record shape, `upstream_md5` reconstruction,
`validate_lock` invariants, `apply_locked_sources` with **shuffled key order** —
the regression test for the positional bug), `test_store_verify.py` (tiny real
store in `tmp_path`: unlink a `.seq` → exactly one `store.seq_file_missing`;
add a sequence → root mismatch; corrupt a `.seq` → L2 fires while L0/L1 stay
clean, proving the levels are independent), `test_verify_cache.py`,
`test_repair.py` (dispatch table; `reingest_scope` expands a shared collection
across namespaces; cache-before-store ordering), `test_status.py`.

**Integration, on the dev manifest** (~3 min per cycle):
`build --config sources.dev.toml --store-dir store.dev --lock build.dev.lock.json`,
then `verify --all`, then mutate (`shard_range` 1–2 → 1–3, or add a prefix to
`record_prefixes`) and re-run `status` and `sync --dry-run`.

**The acceptance test for the whole feature** is the three-way split on the real
preserved drift:

| command | expected |
|---|---|
| `verify --cache --store` | **PASS** — the store matches the lock, which is what it claims |
| `verify --manifest` | **FAIL** — 2 × `manifest.file_removed` (`human.16.*`), 1 × `manifest.file_added` (`refseqgene.9`) |
| `verify --store --deep` | **PASS with 133 known, 0 new** — the encoder baseline |
| `status` | describes all of it, exit 0 |

That split is the point: a legitimately stale lock is not a broken store, and a
known upstream defect is not a regression.

Add one deliberately induced case per level, to prove the levels are
independent and that the baseline does not mask anything:

- unlink a `.seq` → L1 `store.seq_file_missing`, L0 clean
- corrupt a `.seq` → L2 `store.seq_digest_mismatch` at **`error`** (not in the
  baseline), L0/L1 clean
- remove a digest from the baseline → that known defect re-reports as `error`

## Deferred, deliberately

The selenocysteine encoder defect is **not fixed here**. It stays in the store,
pinned as a baseline and reported by `verify --store --deep`, and goes upstream
as a gtars issue. Fixing it would mean either patching gtars or re-encoding
those sequences outside the library's alphabet model, both of which are larger
than this work and neither of which the lock/verify redesign depends on.

Two consequences worth stating plainly rather than discovering later:

- **The published store will contain 133 sequences that cannot be served
  correctly.** They resolve, and return bytes, and those bytes are wrong at one
  residue. Anything downstream comparing a selenoprotein against Ensembl will
  see a mismatch.
- When gtars is fixed, the store must be **rebuilt or those collections
  re-ingested** — the digests are already correct, so only the payloads change.
  `verify --deep` reporting the baseline as `info` is what signals that moment.

---

## Outcome — what actually happened

Everything above landed. Deltas worth recording:

**L2 is faster than estimated.** 2,778 seq/s on `store.pre-filter`, so the full
scan is **10m40s**, not ~18 minutes. `store/` runs slower (~1,950/s) because it
still holds the N-padded scaffolds. `--jobs` was *not* added to L2, as planned.

**The 133 reproduced exactly.** An independent full scan of `store.pre-filter`
returned 106 protein and 27 `dnaio` failures against the same per-alphabet
totals — 1,426,421 / 643 / 107 / 351,881. The estimate in this plan was a
measurement, and it held.

**The dev manifest paid for itself on its first build.** A full build over
`sources.dev.toml` crashed in `validate_lock`: `outputs.collections[].from`
claimed a contributor that `inputs.files` did not contain. Cause —
`ingest_assembly` keys build provenance by `mirror_cache_path(fasta_url)` even
when the FASTA came from a `fasta_path` override, while `resolve_sources`
deliberately emits no remote source for an overridden FASTA. Production has zero
`fasta_path` usage, so nothing else in the repo could have reached it. Fixed by
filtering `from` to recorded input paths, which is what `from` means; the
affected collection is recorded as `"from": []`, and the omission is logged.
A `fasta_path` override is therefore not pinned by the lock — stated in
`sources.dev.toml` rather than left to be rediscovered.

**Two dead functions turned up during the `ensure_download` fix.**
`resolve_fasta_source` and `resolve_report_source` had no callers — superseded
by `require_prepared_source` — and were the only other `ensure_download` call
sites. Deleting them left `ensure_download` with exactly one caller, which made
the lock-guard a parameter rather than plumbing through the assembly path.

**The cross-namespace re-ingest case is hypothetical, not current.** The plan
cited a RefSeq file published under both `GCF/000/001/405/…` and
`annotation_releases/9606/…`; in the committed lock both URLs sit in the same
seqset, so namespace expansion would have covered it. All 31 collections with
contributors in more than one seqset are Ensembl release groups, which the old
expansion also covered. The `from` closure is still the right fix — it is
correct by construction rather than by coincidence, and splitting one seqset's
URL list would break the coincidence — but the comments and the regression test
now say so rather than implying a live bug.

**The run-record `compare.py` was deliberately *not* migrated.** It and its `/2`
lock are hash-pinned artifacts of a sealed run record that validates `PASS`.
Rewriting either would change a hash the manifest pins and falsify the record,
and the script is only ever pointed at its own co-located lock. The rule is now
written down in `RUNBOOK.md`: live code reads locks through `load_lock`; sealed
run-record scripts read their own evidence.

**The dev cycle is slower than the plan's 3-minute estimate, and still worth
it.** Build 6m, `verify --all` 23s, `sync --apply` ~2m,
`status` instant. The `verify --all` run's only findings were two
`remote.size_mismatch` warnings on the mutable RefSeqGene endpoint — a live
confirmation that `remote.*` had to be `warn` rather than `error`.

**Every induced case was run for real, not only in unit tests.** On the dev
store: unlink one `.seq` → exactly one `store.seq_file_missing` at L1 with L0
clean; `repair --store --apply` → collection removed, re-ingested, and the
store's `sequences_root` and `collections_root` back to their pre-damage
values byte for byte. The sync round trip does the same at manifest level —
widening `shard_range` to `[1, 3]` and back returns both roots to exactly what
they were.

**The gtars issue is drafted, not filed, and not tracked.** Filing it upstream
is a separate, outward-facing action, so the draft stays a local working file.
Both reproducers — diffed against their source FASTA records — are recorded in
`runs/2026-09-12-deep-verify/README.md`, which is where they belong regardless
of whether an issue is ever filed.
