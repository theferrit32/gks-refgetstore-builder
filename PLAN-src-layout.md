# Restructure into a `src/` package, and break the import cycle

## Context

Twelve Python modules sit loose in the repo root alongside the manifests, the
lock, the docs and the run records. There is no `__init__.py` anywhere; every
module is top-level and imports its siblings by bare name (`import build_lock`).
That works only because the editable install is a `.pth` file listing the **repo
root itself**, and because eight `sys.path.insert(0, …)` calls — each hardcoding
its own relative depth to the root — patch over the places that mechanism does
not reach.

Two real defects follow, and both are fixed by moving:

* **The wheel cannot import.** Packaging is an explicit
  `[tool.hatch.build.targets.wheel].include` allowlist of six files. It omits
  `store_census.py`, `verify.py`, `repair.py`, `store_sync.py`, `provenance.py`
  and `verify_store.py` — four of which `cli.py` imports at module level. A
  non-editable install's `gks-refgetstore` would `ImportError` on startup.
  Nobody noticed because the editable install bypasses the allowlist entirely.
* **One `sys.path` insert is already wrong.**
  `runs/2026-09-12-deep-verify/artifacts/scan_roundtrip.py` computes
  `parents[2]`, which from `artifacts/` is `runs/`, not the repo root. Harmless
  only because that copy is a frozen snapshot that never runs — but it shows the
  depth-counting is not maintainable by hand.

Intended outcome: one importable package with a working entry point, zero
`sys.path` manipulation anywhere, no module that assumes it lives in the repo
root, and an acyclic module graph.

## Decisions

| | |
|---|---|
| Layout | `src/gks_refgetstore/` — a stale cwd copy cannot shadow the install |
| Imports inside the package | relative (`from .sources import …`) |
| Scope | move **and** break the `build_store` ↔ `build_lock` cycle |
| Standalone scripts | into the package, with console entry points |
| Path defaults | relative to cwd; `fasta_path` relative to its manifest |

Not in scope: decomposing the rest of `build_store.py` beyond the extraction
below, and the run-manifest question (`run_record.py add` copying scripts into
`artifacts/` and pinning them, leaving two copies where one is live) — that is
the agreed follow-up.

---

## Part 1 — extract `sources.py`, breaking the cycle

Done **first, while the files are still flat**, so it is reviewable as a content
change with ordinary diffs before anything moves.

Today there is exactly one cycle in the repo: `build_lock` → `build_store` is
eager (`ResolvedSource`, `SeqsetConfig`, `mirror_cache_path`, `sha256_file`,
plus `build_store.DERIVED_SUFFIXES`), and `build_store` → `build_lock` is
deferred inside `run_build` with a `# lazy: avoids a circular import` comment.
Everything `build_lock` needs is the *source/config model*, not the build
engine — so extracting that model removes the cycle rather than hiding it.

New `sources.py` takes, from `build_store.py`:

- constants `NA_VALUES`, `SEQSET_FORMATS`, `DERIVED_SUFFIXES`
- `RecordExclusion`, `AssemblyConfig`, `SeqsetConfig`, `ResolvedSource`
- `parse_checksum_manifest`, `parse_ncbi_md5checksums`,
  `parse_ensembl_checksum_manifest`, `bsd_sum_file`, `natural_sort_key`,
  `_fetch_text`, `resolve_sources`, `load_config`
- `mirror_cache_path`, `md5_file`, `sha256_file`, `provider_checksum_matches`

`build_store.py` keeps the engine: the gbff/lrg/filter conversion,
`ensure_download`, `check_free_space`, the alias machinery, `ingest_assembly`,
`process_seqset`, `process_release_groups`, `run_build`, and the stats
dataclasses. It drops ~550 lines.

Resulting module-level graph, acyclic:

```
sources       -> ()
store_census  -> ()
build_lock    -> sources, store_census
build_store   -> sources, build_lock      (now EAGER)
store_sync    -> sources, build_store, build_lock
verify        -> sources, build_lock, store_census
fetch_sources -> sources, build_lock
repair        -> sources, build_lock, store_census, fetch_sources, store_sync, verify
provenance    -> build_lock, store_census
verify_store  -> sources
inventory_sources -> sources
cli           -> build_lock, build_store, fetch_sources, repair, store_sync, verify
```

Two consequences to handle deliberately:

- **`build_lock` stops importing `build_store` at all**, so the module handle
  and the contract comment at its top both go. `build_store` imports
  `build_lock` eagerly and the `# lazy` comment goes too.
- **`sha256_file`'s docstring is now wrong.** It currently explains it lives in
  `build_store` "because `ensure_download` needs it … and `build_lock` imports
  this module". After the split both import it from `sources`; rewrite the
  comment to say what is true.
- **`DERIVED_SUFFIXES` loses locality** with the resolvers that produce those
  names (`DERIVED_FASTA_RESOLVERS` stays in `build_store`). Cross-reference both
  ways in comments — a resolver inventing an unregistered suffix would silently
  break `build_lock.records_from_log`.

Keep `build_lock`'s `sha256_file` re-export (it is in `__all__`, and several
tests use `build_lock.sha256_file`).

Then hoist the four nested imports that were never cycle workarounds —
`verify.verify_manifest`, `verify.run_status`, `repair._collections_publishing`,
`repair.repair_store`, and both in `store_sync` — to module level.

**Also flagged:** `build_store.iter_source_urls` has no callers anywhere. Its
docstring claims `fetch_sources` uses it; `fetch_sources` uses `resolve_sources`
instead. Same shape as the two dead resolvers deleted last change. Recommend
deleting it with the extraction rather than carrying a function whose docstring
names a caller that does not exist.

## Part 2 — move into `src/gks_refgetstore/`

A pure `git mv` plus mechanical import rewrites, so the diff reads as a rename.

```
src/gks_refgetstore/
    __init__.py          docstring only
    cli.py  sources.py  build_lock.py  build_store.py  store_census.py
    verify.py  repair.py  store_sync.py  fetch_sources.py  provenance.py
    verify_store.py  inventory_sources.py
    generate_ncbi_source_candidates.py
```

`__init__.py` stays deliberately near-empty: re-exporting submodules there would
make importing the package pull in `gtars` and every module, undoing the lazy
CLI startup. No `__version__` either — it would duplicate the static version in
`pyproject.toml` and drift; nothing currently needs one.

Stays at the repo root: `sources.toml`, `sources.dev.toml`, `build.lock.json`,
`tests/`, `tools/`, `runs/`, `seqrepo_equivalence/`, the docs, `upload_store.sh`.
`seqrepo_equivalence/*.py` compute `REPO_ROOT = HERE.parent`, which stays
correct precisely because that directory does not move.

**Delete seven of the eight `sys.path.insert` calls**, and the ~10
`# noqa: E402` suppressions that exist only to allow imports after them. Nothing
needs them once the package is installed. (The eighth is the frozen
`artifacts/scan_roundtrip.py`, left untouched.)

| file | becomes |
|---|---|
| `tests/conftest.py` | nothing — pytest finds the installed package |
| `tools/compare_store_mutation.py` | `from gks_refgetstore.store_census import …` |
| `seqrepo_equivalence/full_parity.py` | `from gks_refgetstore import build_lock` |
| `runs/2026-09-12-lock-schema-v4/convert.py` | `from gks_refgetstore import build_lock` |
| `runs/2026-09-12-deep-verify/scan_roundtrip.py` | `from gks_refgetstore import store_census` |
| `fetch_sources.py`, `provenance.py` | nothing — same package |

The `runs/` working scripts are **not** hash-pinned (only
`artifacts/scan_roundtrip.py` and `2026-07-23/compare.py` are), so editing them
is safe and `run_record.py validate` keeps passing. `2026-07-23/compare.py` has
no first-party imports and needs no change. The frozen
`artifacts/scan_roundtrip.py` is left exactly as it is — it is the snapshot of
what ran.

`pyproject.toml`:

```toml
[project.scripts]
gks-refgetstore       = "gks_refgetstore.cli:main"
gks-refgetstore-check = "gks_refgetstore.verify_store:main"

[tool.hatch.build.targets.wheel]
packages = ["src/gks_refgetstore"]
```

The `include` allowlist goes entirely — that is the bug. `sources.toml` drops
out of the wheel: with cwd-relative defaults it is *user* data, not package
data, which is the correct model for a build tool. `verify_store.py` already has
a `main() -> int`, so the second entry point needs no extraction.

Keep `[tool.pytest.ini_options] testpaths`; do **not** add `pythonpath`. Tests
resolving through the install is the point of the src layout.

## Part 3 — path defaults

Delete `REPO_ROOT` from all four modules that define it.

- `cli.py`: `DEFAULT_CONFIG = Path("sources.toml")`, and likewise `downloads`,
  `store`, `build.lock.json`. `inventory_sources.py` already does exactly this,
  so it is the in-repo precedent rather than a new convention.
- `verify.py`: `DEFAULT_KNOWN_BAD =
  Path("seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv")`.
- `verify_store.py`: its five `HERE /` constants — `STORE_PATH`, `DOWNLOAD_DIR`,
  `CONFIG_PATH`, `KNOWN_DIVERGENT_PATH` — become cwd-relative. Both fixtures
  stay in `seqrepo_equivalence/known_divergence/`, so their shared README stays
  true.
- `provenance.py`: `--lock` and `--store-dir` defaults.
- `build_store.py`: the single functional use, `(REPO_ROOT / entry.fasta_path)`
  in `ingest_assembly`. Resolve it **once in `load_config`**, against the
  manifest's own directory, onto a new non-init `AssemblyConfig` field —
  following the `resolved_sources` / `exclusion` pattern `SeqsetConfig` already
  uses. `ingest_assembly` then reads the resolved path and `build_store` has no
  repo-root concept left.

Behaviour change to state plainly in the README: the CLI now resolves relative
defaults against the **current directory**, so it must be run from the repo root
or given explicit paths. `fasta_path` in a manifest is now relative to that
manifest rather than to the repo root — for `sources.dev.toml`, which lives at
the root, these are the same path, so no data migration.

## Part 4 — docs

Small surface, because most references are to directories that do not move.

- **`README.md`** — the Layout tree (which also currently omits `tests/`,
  `tools/` and `generate_ncbi_source_candidates.py`); `uv run cli.py <cmd>` →
  `uv run python -m gks_refgetstore.cli`; the `verify_store.py` /
  `provenance.py` invocations → `gks-refgetstore-check` and
  `python -m gks_refgetstore.provenance`; the config-reference table row that
  reads *"`fasta_path` … relative to repo root"*; and a note on cwd-relative
  defaults. The many `gks-refgetstore <cmd>` examples are unaffected — the
  console script name does not change.
- **`sources.dev.toml`** — its header gives "`fasta_path` resolves relative to
  the repo root" as one of three reasons the dev cache must be `downloads/`, and
  repeats it in the `GRCh37.p13` block. The *conclusion* survives (the manifest
  is at the root, so the resolved path is identical), but the reasoning has to
  be restated as manifest-relative.
- **`seqrepo_equivalence/known_divergence/README.md`** — two mentions of
  `verify_store.py` as the fixture's consumer.
- **`upload_store.sh`** — one prose mention of `build_store.py`. The script
  itself must stay at the repo root: its `HERE` drives `uv run --project`,
  `git -C`, the default store dir and `sources.toml`.
- **`RUNBOOK.md`** — no changes. Every path it cites is under `tools/`, `runs/`
  or `seqrepo_equivalence/`.
- **`PLAN-*.md`, `ISSUE-*.md`, `runs/**/README.md`, and the sealed
  `run-manifest.json` files** — deliberately **not** updated. They are records
  of past states: the PLANs carry `build_store.py:1806-1830`-style line
  citations, and two manifests record
  `validations[].command = ["uv","run","python","verify_store.py", …]`.
  Those commands will no longer run, which is correct — they describe what was
  executed then. Rewriting any of it would falsify history, the same principle
  `RUNBOOK.md` already applies to sealed records.

While `.gitignore` is open, two bits of hygiene: `/parity_membership.py` and
`/probe_source_coverage.py` are dead root-anchored rules left over from when
copies sat at the root (the real files are tracked under
`seqrepo_equivalence/`), and `.DS_Store`, `.pytest_cache/` and `.ruff_cache/`
are currently not ignored at all. Also delete the stale root `__pycache__/`,
which still holds `.pyc` files for two modules that no longer have source
(`lock_sources`, `verify_cache`) — harmless once the root leaves `sys.path`,
but it is exactly the confusion the move is meant to end.

## Verified unaffected

Worth stating so they are not re-investigated mid-implementation:

- **All four lock files.** `cache_path` is relative to the *cache directory*,
  never the repo root, and no lock stores a `.py` path.
- **`uv.lock`** — `source = { editable = "." }` is directory-level.
- **`tests/` and `tools/` do not move**, so `tests/conftest.py`'s
  `parent.parent`, `test_store_sync.py`'s `repo / "build.lock.json"`, and
  `test_run_record_streaming.py` invoking `tools/run_record.py` by filesystem
  path all keep working.
- **`seqrepo_equivalence/`'s own `HERE.parent` roots** and its sibling-module
  imports (`from verify_seqrepo_equivalence import …`, which resolve through
  Python's script-directory rule).
- **`.claude/settings.local.json`** contains only a `WebFetch` domain
  allowlist — no paths.
- **`runs/2026-07-23-…/compare.py`** has no first-party imports.

## Order of work

1. Extract `sources.py`; make `build_store` → `build_lock` eager; hoist the
   non-circular nested imports; delete `iter_source_urls`. Flat layout, tests
   green. **Commit.**
2. `git mv` into `src/gks_refgetstore/`; rewrite imports; delete the `sys.path`
   inserts and stale `noqa`s; rewrite `pyproject.toml`; `uv sync`. **Commit.**
3. Path defaults + `fasta_path` resolution + docs. **Commit.**

Splitting 1 from 2 is what keeps the move mechanically verifiable.

## Verification

**Install and imports**

    uv sync
    uv run python -c "import gks_refgetstore.cli as m; print(m.__file__)"   # under src/
    grep -rn "sys.path" --include=*.py . | grep -v .venv                   # empty

**No cycle, as a regression guard.** Add a test that imports each package module
in a fresh subprocess, in isolation — a reintroduced cycle fails in one import
order only, which a normal test run would not catch.

**Suite and lint**

    uv run pytest              # 222 pass, plus the new cycle + fasta_path tests
    uvx ruff check --exclude .venv .

**The packaging bug, which has never been verified.** This is the check that
proves the change did what it claims:

    uv build
    # install the wheel into a throwaway venv, then:
    gks-refgetstore --help          # must not ImportError

**CLI surface**

    gks-refgetstore --help; for each of build verify status repair sync fetch lock
    gks-refgetstore-check --help
    gks-refgetstore status --offline        # roots agree, exit 0
    gks-refgetstore verify --store          # PASS, ~95s on store/
    uv run python -m gks_refgetstore.provenance --list   # 240 collections, 0 unattributed

**Dev manifest, fast path** (no rebuild needed — `store.dev/` and
`build.dev.lock.json` are current):

    gks-refgetstore fetch  --config sources.dev.toml --dry-run   # 19 skip, ~8s
    gks-refgetstore verify --config sources.dev.toml --store-dir store.dev \
        --lock build.dev.lock.json --all                         # PASS, 2 remote warns, ~23s

**The scripts that reach into the package**

    uv run python tools/run_record.py validate runs/*/                 # 6 PASS
    uv run python runs/2026-09-12-deep-verify/scan_roundtrip.py \
        --store store.pre-filter --out /tmp/rt.tsv --limit 5000        # 0 failures
    uv run python runs/2026-09-12-lock-schema-v4/convert.py --help     # imports
    uv run python tools/compare_store_mutation.py --baseline store --mutated store.pre-filter

**`fasta_path` is the one behaviour with no cheap end-to-end check** — a full dev
build is ~6 minutes. Cover it with a unit test asserting `load_config` resolves
`fasta_path` against the manifest's directory for a manifest written to
`tmp_path`, and treat a dev rebuild as optional final confirmation.
