# Reproducible run records

`runs/` keeps compact evidence for builds, publication, parity analysis, and
isolated experiments. A record is a directory named `YYYY-MM-DD-short-name`
with a `run-manifest.json`, a concise `README.md`, and only the small artifacts
needed to audit its result.

## Lifecycle

1. **Initialize.** Start from the repository root. Use `--require-clean` when a
   run is meant to describe committed code, or `--reconstructed` only when
   recording a historical run whose original capture is incomplete.

       uv run python tools/run_record.py init runs/YYYY-MM-DD-name \
         --kind build --title "Locked production build" --require-clean

2. **Execute.** Record each meaningful command separately. Output is streamed
   unchanged to the terminal and `logs/<label>.log`. A nonzero command is
   recorded with its exit code and marks the record failed before the same
   code is returned to the caller.

       uv run python tools/run_record.py exec runs/YYYY-MM-DD-name \
         --label tests -- uv run pytest

3. **Attach small artifacts.** Give each input or output a role. Files are
   copied into `artifacts/` and recorded with their byte size and SHA-256.

       uv run python tools/run_record.py add runs/YYYY-MM-DD-name \
         --role build-lock build.lock.json

   The default limit is 1 MiB. `--allow-large` is an explicit escape hatch,
   not the normal policy: reproducible bulk tables, stores, downloaded FASTAs,
   and caches should be represented by schemas, row counts, hashes, and compact
   summaries rather than retained in Git.

4. **Finalize.** Add concise `result_summary`, `validations`, and `publication`
   facts to the manifest where applicable, then set the outcome.

       uv run python tools/run_record.py finalize runs/YYYY-MM-DD-name \
         --status succeeded

5. **Validate.** Validation checks the schema, directory/run ID relationship,
   safe artifact paths, artifact sizes and hashes, and recorded command logs.

       uv run python tools/run_record.py validate runs/YYYY-MM-DD-name

6. **Review for a future commit.** Inspect the manifest, logs, artifact sizes,
   Markdown links, and `git status`. Run records are not automatically staged
   or committed.

The manifest captures exact argv, working directory, exit code, log identity,
Git commit and dirty state, a diff hash, allowlisted tool versions, artifact
roles and identities, result summaries, validation outcomes, publication
destination, evidence quality, and missing evidence. It deliberately never
captures environment variables or credentials.

## Records are sealed; they are not migrated

Every artifact in a finished record is pinned by SHA-256 in its manifest, and
`validate` checks those hashes. A record therefore **does not move with the
code**. When a schema changes, its scripts and its data stay as they were.

`runs/2026-07-23-ensembl-r116-genomic-experiment/compare.py` reads a build lock
with a raw `json.loads` and expects the `/2` shape. That is correct and stays:
it is only ever pointed at `build.ensembl-dna-r116.lock.json`, the `/2` artifact
sealed beside it, and the two are internally consistent forever. Rewriting
either to the current schema would change a hash the manifest pins and falsify
the record — a worse outcome than a script that only reads its own evidence.
The same applies to `runs/2026-07-22-published-store/build.lock.json`.

The rule: **live code reads locks through `build_lock.load_lock` and the
accessors; sealed run-record scripts read their own co-located artifacts.** If a
historical analysis needs re-running against current data, that is a new record,
not an edit to an old one.

## Workflow recipes

### Locked build

Initialize with `--require-clean`. Record cache verification, the locked build,
the test suite, and both store verifications as separate `exec` commands. Add
the config, build lock, store manifest, and compact verifier summary. Never add
the store or download cache.

Use `sources.dev.toml` for the edit-test loop before pointing anything at the
production manifest; a dev build and a full `verify --all` take minutes rather
than hours. A dev run is not a run record — records describe the real thing.

### Regenerating the `--deep` encoder baseline

`seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv` pins the
sequences whose stored payload does not re-digest to its own digest. It is
regenerated only when the set genuinely moves — a gtars upgrade, or a store
rebuild that adds affected sequences.

    uv run python runs/YYYY-MM-DD-deep-verify/scan_roundtrip.py \
      --store store --out roundtrip.tsv

Record the scan as an `exec` command, add the resulting TSV as an artifact, and
state in the record's README which store was scanned and what the per-alphabet
totals were. The `cause` column is a **diagnosis**, not a measurement: verify
cannot infer it, because at verify time an encoder defect and bit rot are
indistinguishable. Carry causes forward rather than re-deriving them, and
investigate any row that lands as `undiagnosed`.

A baseline entry that starts round-tripping correctly is reported by
`verify --deep` as `info`, not silently dropped. That is the signal that gtars
was fixed, and the point at which the affected collections need re-ingesting —
their digests were always right; only the payloads change.

### Publication

Record upload dry-run and live upload as separate commands and labels. Use a
fresh local manifest for the exact store being uploaded; record destination,
object and byte totals, transfer/checker concurrency, and `rclone check`
results. Do not include remote credentials, rclone configuration, or tokens.

### Seqrepo parity

Run `seqrepo_equivalence/verify_seqrepo_equivalence.py` with explicit store,
seqrepo, report, and gap-list paths. Run
`seqrepo_equivalence/parity_membership.py` with explicit output paths. Retain
the report, compact summary, and a reasonably sized build log. For the
multi-million-row membership/gap tables and probe caches, retain only schema,
row count, byte size, SHA-256, generator, and headline categories.

### Isolated Ensembl experiments

Use a dedicated config, lock, cache, and copy-on-write playground. Pass
`--playground`, `--baseline`, `--lock`, `--out-tsv`, and `--out-json` explicitly
to the comparison script. Verify the preserved baseline before and after the
experiment, and state clearly that the playground is not publishable.
