# gks-refgetstore-builder

Build a custom [`gtars.refget.RefgetStore`](https://github.com/databio/gtars)
archive augmented with the per-sequence and per-collection aliases that GA4GH
tooling (e.g. [vrs-python](https://github.com/ga4gh/vrs-python)) needs to use
refget as its **alias backend** — a drop-in replacement for seqrepo's alias
mapping.

The store maps human-readable accessions to GA4GH sequence digests across:

- Assembly aliases: `GRCh38:chr1`, `GRCh38.p14:chrM`, etc.
- Cross-assembly genomic accessions: `refseq:NC_000001.11`, `insdc:CM000663.2`.
- Collection-level: `refseq:GCF_*`, `insdc:GCA_*`.
- **NCBI RefSeq transcripts + proteins**: `refseq:NM_000551.3`,
  `refseq:NP_000542.1`, etc., loaded from NCBI RefSeq mRNA/Prot shards.
- **Ensembl genomic, transcript, and protein sequences** from releases 75–116.
  Immutable aliases live under `ensembl-N`; `ensembl` is an exact rolling
  snapshot of release 116. Each release includes unmasked `dna.toplevel`, cDNA,
  ncRNA, and peptide FASTAs.

No dependency on `bioutils` or `biocommons.seqrepo` — the loader reads NCBI
assembly reports directly and derives transcript/protein aliases from FASTA
header names. The only runtime dependency is `gtars`.

> Extracted from the `misc/refgetstore/` directory of vrs-python. The
> functional-equivalence harnesses that compare a built store against seqrepo
> through vrs-python's dataproxy remain in vrs-python, since they depend on it.

## Prerequisites

- Python >= 3.11 (uses `tomllib`).
- [`uv`](https://docs.astral.sh/uv/) for environment + dependency management.
- `gtars`, installed via `uv sync`. It is currently built from source at a
  pinned commit (see [gtars revision](#gtars-revision)), so `uv sync` also needs
  a Rust toolchain (`rustup`, stable).
- Network access to `ftp.ncbi.nlm.nih.gov` and `ftp.ensembl.org` (or
  pre-populated cache dirs).
- A **case-sensitive filesystem for the store directory**. The defaults on
  macOS (APFS) and Windows (NTFS) are case-insensitive, and a store written to
  one cannot be published. See
  [Case-sensitive store directory](#case-sensitive-store-directory).
- Disk: the default manifest includes 23 NCBI assembly releases, complete
  published GRCh37/GRCh38 RefSeq release history, and all four canonical FASTA
  classes for Ensembl releases 75–116.
  Its cache and store therefore need substantially more capacity than a
  current-only build. Use `gks-refgetstore fetch --dry-run` to enumerate the
  selected files before provisioning storage.

## Setup

    uv sync

### gtars revision

gtars is pinned to a commit rather than a PyPI release:

| | |
| --- | --- |
| repository | [`theferrit32/gtars`](https://github.com/theferrit32/gtars), a fork of [`databio/gtars`](https://github.com/databio/gtars) |
| branch | [`fix/encodings-buildable`](https://github.com/theferrit32/gtars/tree/fix/encodings-buildable) |
| commit | [`7831e09f596042e3983c7a123c3f6c019ca6c05f`](https://github.com/theferrit32/gtars/commit/7831e09f596042e3983c7a123c3f6c019ca6c05f) |
| package version | `0.11.0` (unreleased; no 0.11.0 wheels were published to PyPI) |

The commit is the head of [databio/gtars#273](https://github.com/databio/gtars/pull/273)
(`bbe7fd97`, "Fix sequences that change when stored in encoded mode") plus one
commit that adds the missing `StorageMode::Zstd` case to the Python and R
bindings, without which the Python package does not compile. #273 corrects the
encoded-storage round trip for `D`/`H` and RNA `U` in the IUPAC nucleotide
alphabet and for `U`/`B`/`Z`/`O`/`J` in the protein alphabet. With an earlier
gtars, 133 sequences in the store are returned with different residues than
were imported; see
[`issues/gtars-encoder-alphabet/`](issues/gtars-encoder-alphabet/README.md).

The pin lives in `pyproject.toml` under `[tool.uv.sources]`, and `uv.lock`
records the resolved commit. Every build lock records it too: alongside
`gtars_version`, the lock's `build.gtars_source` names the repository and commit
gtars was built from, read from the installed package's metadata (`null` for a
PyPI wheel).

Stores written with this revision use residue codes that gtars 0.10 and earlier
do not know about, and those versions decode them as different letters without
raising an error (selenocysteine `U` reads back as `X`). The store format
version is unchanged, so nothing prevents an older gtars from opening such a
store. Read stores built here with this revision or a later release that
includes #273.

Move to a PyPI release once one includes #273: replace the
`[tool.uv.sources]` entry with a version constraint, then `uv lock`.

### Case-sensitive store directory

gtars stores each sequence payload at
`sequences/<first two characters of its digest>/<digest>.seq`. Digests are
base64url, so `zz`, `Zz`, `zZ` and `ZZ` are four different shard directories.

macOS APFS and Windows NTFS are case-insensitive by default (case-preserving,
but lookups ignore case), so on them those four shards become one directory,
named after whichever was created first. Nothing local notices: lookups ignore
case, so the store builds, reads and passes `verify --deep`. The damage shows up
once the store leaves that filesystem. Copied to a case-sensitive filesystem, or
uploaded to object storage such as R2 or S3, whose keys are case-sensitive,
about two thirds of payloads sit at a path no reader requests. A full store
needs 4,096 shard directories; a case-insensitive filesystem holds at most
1,444.

`build`, `sync --apply` and `repair --apply` check this before doing anything
else and refuse a store directory on a case-insensitive filesystem. The check
writes a probe file in the store directory, or in its nearest existing parent,
tests whether the uppercase spelling of its name refers to the same file, and
removes it. `--allow-case-insensitive-fs` skips the refusal for a store that is
only ever read in place, such as a local dev store.

**On macOS**, put the store on a case-sensitive APFS volume. A volume shares
free space with the other volumes in its container, grows and shrinks with its
contents, and needs no size set in advance:

    diskutil apfs list                                  # find the container, e.g. disk3
    diskutil apfs addVolume disk3 APFSX RefgetStores    # mounts at /Volumes/RefgetStores

Creating it asks for an administrator's authorization, because it changes the
internal disk's container. Add `-passprompt` to encrypt it. The volume belongs
to the user who creates it, and needs no special rights to use afterwards.

Then build onto it, and point `./store` at it so every command's default
`--store-dir` keeps working (the checkout this README describes is set up this
way):

    gks-refgetstore build --store-dir /Volumes/RefgetStores/store
    ln -s /Volumes/RefgetStores/store store

`.gitignore` ignores the `store` symlink. The case check follows it, so it
passes for the volume, not for the repository's own filesystem. The download
cache can stay on the default volume, since its file names do not depend on
case.

**On Linux**, ext4 and xfs are case-sensitive by default, and nothing needs
setting up. **On Windows**, build inside WSL2's own Linux filesystem rather than
under `/mnt/c`.

A store already built on a case-insensitive filesystem cannot be fixed where it
is, because each of its shard directories holds several prefixes' files and the
filesystem cannot separate them. Either copy it to a case-sensitive filesystem
and move each payload into the directory named by its own digest prefix, or
rebuild it there: with `--locked-sources` and an intact cache a full rebuild
takes about half an hour and reproduces the same digests and payload bytes.

## Layout

    sources.toml                   # declarative, authoritative source manifest
    sources.dev.toml               # every config shape, ~3 GB: the fast edit-test loop
    build.lock.json                # provenance lock for the most recent build (see below)
    src/gks_refgetstore/
        cli.py                     # gks-refgetstore CLI (build/verify/status/repair/sync/fetch/lock)
        sources.py                 # the source model: manifest schema, resolution, cache paths, digests
        build_store.py             # build engine: download, convert, filter, ingest, alias
        build_lock.py              # the build lock: schema, accessors, validation, sync planning
        store_census.py            # read-only store primitives (counts, roots, membership, re-digest)
        verify.py                  # verify + status: cache, store, manifest, remote
        repair.py                  # restore the locked state after verify finds damage
        store_sync.py              # incremental add/remove against a built store
        fetch_sources.py           # cache pre-fetch (fetch subcommand) + the shared downloader
        inventory_sources.py       # exact remote sizes + pre-download space-budget gate
        provenance.py              # file <-> refget-digest queries from the lock + store
        verify_store.py            # post-build gtars self-checks (gks-refgetstore-check)
        generate_ncbi_source_candidates.py  # draft manifest entries from an NCBI assembly listing
    tests/                         # pytest suite (resolves through the installed package)
    tools/                         # run records, store diffing, upstream drift probes
    seqrepo_equivalence/           # optional backwards-compat check vs a seqrepo snapshot
    RUNBOOK.md                     # reproducible run-record lifecycle and policies
    runs/                          # compact historical manifests, summaries, logs, evidence
    downloads/                     # cached FASTA + assembly_report.txt (gitignored)
    store/                         # local output/playground RefgetStore (gitignored)
    store.2026-07-22/              # preserved published store (local, gitignored)

The package is a `src/` layout, so the code only ever resolves through the
install (`uv sync`) — a stale copy in the working directory cannot shadow it.

## CLI

Everything runs through one CLI, `gks-refgetstore`. Two equivalent invocation
styles:

    gks-refgetstore <cmd> ...                    # installed console script (after `uv sync`)
    uv run python -m gks_refgetstore.cli <cmd> ...   # same code, via the module

| subcommand | does | exits non-zero on |
| --- | --- | --- |
| `build` | ingest the manifest into a RefgetStore, write the lock | source drift, ingest failure |
| `verify` | check the cache and the store against the lock | any error-severity finding |
| `status` | describe manifest/lock/cache/store agreement | **never** |
| `repair` | restore the locked state after `verify` finds damage | unrepairable or failed repair |
| `sync` | apply `sources.toml` changes incrementally | missing sources |
| `fetch` | pre-populate the download cache (no ingest) | download failure |
| `lock` | regenerate a lock for a build that already ran | — |

Global: `-v/--verbose` for debug logging. Every subcommand shares
`--config` (default `./sources.toml`) and `--cache-dir` (default `./downloads`).
`gks-refgetstore <cmd> --help` prints the full flag list.

Every relative default — `sources.toml`, `downloads/`, `store/`,
`build.lock.json` — resolves against the **current directory**, so run the CLI
from the repo root or pass the paths explicitly. The manifest, the cache and the
store are your data; the installed package does not reach back into itself for
them.

### The development manifest

`sources.dev.toml` covers every config *shape* in `sources.toml` with 19 files
and ~3.19 GB, against 335 files and 106 GB. Everything it names is already
cached and still matches its provider checksum, so it downloads nothing.

Measured: build ~6 minutes, `verify --all` 23 seconds, `sync --apply` ~2
minutes, `status` instant.

    gks-refgetstore build  --config sources.dev.toml --store-dir store.dev --lock build.dev.lock.json --allow-case-insensitive-fs
    gks-refgetstore verify --config sources.dev.toml --store-dir store.dev --lock build.dev.lock.json --all
    gks-refgetstore status --config sources.dev.toml --store-dir store.dev --lock build.dev.lock.json

`store.dev/` is only ever read in place, so on a case-insensitive filesystem it
is built with `--allow-case-insensitive-fs`. Drop the flag, or point
`--store-dir` at a case-sensitive volume, if the dev store is to be copied or
uploaded.

It shares `downloads/` with the production manifest rather than using its own
cache directory. That is forced, not incidental: `fasta_path` resolves relative
to the manifest, which sits at the repo root alongside `downloads/`; derived
artifacts are written beside their source; and re-deriving release 76's filtered
toplevel costs 273 seconds. The header of
`sources.dev.toml` lists which entry covers which shape, and which one shape is
deliberately left out.

### `build` — ingest the manifest into a RefgetStore

    gks-refgetstore build                          # build the whole manifest
    gks-refgetstore build --assembly GRCh38        # only one assembly (skips seqsets)
    gks-refgetstore build --seqset refseq_human_rna  # only one seqset (skips assemblies)
    gks-refgetstore build --skip-seqsets           # assemblies only

Reads the cache (never re-downloads unless `--force-download`); a missing file
is fetched on demand. Writes/refreshes the build lock at the end (see
[Build lock](#build-lock--verification)).

    --store-dir PATH       output RefgetStore (default: ./store)
    --assembly NAME        only this assembly namespace; --seqset NAME only this seqset
    --skip-assemblies / --skip-seqsets
    --force-download       re-fetch even if cached
    --ingest-jobs N        FASTAs imported concurrently per seqset (default: min(8, cores))
    --filter-jobs N        source files filtered concurrently in the preflight (default: min(8, cores))
    --min-free-gb N        refuse to start if free disk is below this (default 25)
    --allow-case-insensitive-fs
                           build even if --store-dir is case-insensitive (local-only stores; see Prerequisites)
    --lock PATH            build-lock to check against + write (default ./build.lock.json)
    --lock-check-mode {strict,subset,ignore}   pre-flight check vs the lock (default strict)
    --no-lock              don't write the lock (the pre-flight check still runs)
    --force-lock           accept current inputs and write a new lock baseline
    --locked-sources       use the lock's concrete URLs and cached SHA-256 bytes;
                           perform no live pattern resolution

Building into an existing store **appends** (dedup by digest); to rebuild from
scratch, delete the store dir first.

Each seqset's shards are imported in **one batched call**. gtars persists the
global sequence index once per import call, so importing files one at a time is
O(N²) in store size — measured at 9.9x slower for identical work against a 32 MB
index than an empty one, and a full build's index exceeds 90 MB. Batching
amortizes that persist (3.4x on an eight-shard seqset) and `--ingest-jobs`
parallelizes the decompress/digest work that dominates the large Ensembl inputs
(3.0x on `dna.toplevel` at 3 jobs). Results are identical either way — sequences
are content-addressed, so collection digests do not depend on import order.

The seqset is the batching unit deliberately: batching more widely measured only
~4% faster, while release-scoped seqsets must be ingested and published
oldest-first for the rolling namespace to end up holding exactly the newest
release. Lower `--ingest-jobs` if memory is tight; peak RSS grows by roughly
0.45 GiB per concurrent `dna.toplevel`.

Free space is checked once before ingestion rather than monitored throughout —
a full build runs for hours, and learning at the end that the volume filled is
worse than being told at the start. The projection is a rough multiple of the
compressed source bytes, so it warns rather than blocks; only `--min-free-gb`
aborts.

### `fetch` — pre-populate the download cache (no ingest)

    gks-refgetstore fetch                # download every manifest source
    gks-refgetstore fetch --dry-run      # list what would be pulled
    gks-refgetstore fetch --only ensembl_release_116  # one complete release
    gks-refgetstore fetch --kinds assembly_report     # just the small reports
    gks-refgetstore fetch --locked-sources            # offline lock/cache check

Idempotent (skips files already cached), atomic, disk-guarded
(`--min-free-gb`). Useful to warm the cache before a long build, or offline.

### `verify` — check the cache and the store

Four independent targets. With no scope flag the default is `--cache --store`,
which is **offline and lock-driven**; only `--remote` and `--manifest` touch the
network.

    gks-refgetstore verify                  # cache + store, offline
    gks-refgetstore verify --store --deep   # + re-digest every stored sequence
    gks-refgetstore verify --manifest       # has upstream moved since the lock?
    gks-refgetstore verify --all            # all four

| target | question | network |
| --- | --- | --- |
| `--cache` | are the pinned source bytes still on disk and intact? | no |
| `--store` | does the store still hold what the lock recorded? | no |
| `--manifest` | has upstream moved since the lock was written? | yes |
| `--remote` | are the remote objects still the size we saw? | yes |

The split is the point. On the currently committed lock:

    verify --cache --store   PASS   the store matches the lock, which is what it claims
    verify --manifest        FAIL   2 files withdrawn, 1 added, 30 checksums moved
    status                   describes all of it, exit 0

A legitimately stale lock is not a broken store. The old single `verify` could
not tell them apart, because it resolved the manifest over the network to decide
what to check — which also made offline verification impossible and made the
result depend on what upstream happened to be serving that day.

Every finding carries a stable `code` (`cache.sha256_mismatch`,
`store.seq_file_missing`, `manifest.file_removed`, …). That is the contract
`repair` dispatches on, so the codes are interface, not log text.

`--remote` findings are always warnings: several endpoints omit
`Content-Length`, and mutable endpoints legitimately change size. `--manifest`
answers the same question better at the same network cost.

#### The store ladder

| level | proves | does not prove | cost |
| --- | --- | --- | --- |
| **L0** roots + counts | the store holds **exactly** the digest sets the lock recorded — one add or remove flips a root | that any payload exists or is readable | ~7s |
| **L1** membership | every locked collection resolves, every sequence is backed by a `.seq` file, nothing is stranded | that any byte is correct | ~85s |
| **L2** `--deep` | the store can actually serve each sequence, and returns what its digest promises | nothing about the lock — it needs no lock | ~15 min |

The levels are independent: unlinking a `.seq` fires L1 while L0 stays clean,
and corrupting one fires L2 while L0 and L1 stay clean.

`--limit` is accepted for `--cache` and `--remote` and **rejected for
`--store`** — a root is a digest over the complete set, so a partial one is not
a weaker check, it is a different number that means nothing.

#### Encoder round-trip baseline (`--deep`)

`--deep` re-digests every stored sequence from the bytes the store returns and
compares the result with the digest it is filed under. Sequences known not to
round-trip are listed in a **baseline**
(`seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv`) rather
than suppressed, which gives three outcomes instead of two:

| in baseline | re-digests correctly | outcome |
| --- | --- | --- |
| yes | no | `warn` — known defect, rendered with its diagnosed cause |
| **no** | **no** | **`error`** — a new defect; this is the regression that matters |
| **yes** | **yes** | `info` — the encoder was fixed; update the baseline and re-ingest |

**The baseline is empty**, and a full `--deep` run over the store passes with 0
errors. Stores built with gtars 0.10 or earlier return 133 sequences with
different residues than were imported (106 selenoproteins whose `U` is returned
as `A`, and 27 `dnaio` sequences whose `D`/`H` are returned as `H`/`V`). That is
what the [gtars revision](#gtars-revision) pin fixes. The earlier list is kept in
[`runs/2026-10-07-gtars-pr273-rebuild/`](runs/2026-10-07-gtars-pr273-rebuild/README.md);
see also
[`seqrepo_equivalence/known_divergence/README.md`](seqrepo_equivalence/known_divergence/README.md).

### `status` — describe, don't judge

    gks-refgetstore status             # four legs, one of them online
    gks-refgetstore status --offline   # skip the upstream leg

Four legs: manifest↔lock, manifest↔cache, lock↔cache, lock↔store. Each reuses
the implementation `verify` and `sync` already use — there is no second copy of
"resolve and diff".

**`status` always exits 0.** It is meant to be run casually and read; `verify`
is the one that judges and gates CI. A command that fails on a stale lock is a
command people stop running.

### `repair` — restore the locked state

    gks-refgetstore repair             # verify, then plan; changes nothing
    gks-refgetstore repair --apply     # actually repair

Re-fetches damaged cache files, **accepting only bytes that match the lock's
SHA-256**, and re-ingests collections the store has lost. It restores what the
lock records and never re-baselines — new upstream bytes are a `sync` or
`--force-lock` decision.

Two deliberate refusals:

- **`store.seq_digest_mismatch` is not repaired.** Repair re-ingests from the
  same source bytes, and gtars skips a sequence whose digest the store already
  holds, so the existing payload would be left in place. A mismatch means an
  encoder defect or a damaged payload, and needs diagnosing before anything is
  rewritten. A store built with a gtars that has an encoder defect is fixed by
  rebuilding with a gtars that does not.
- **A root mismatch with no attributable cause is refused.** If the store's
  digest set differs from the lock's but no collection is missing and no payload
  is absent, the store has diverged in a way repair cannot name. That is `sync`.

Cache repair always runs first and aborts the run if it fails: store removal
persists immediately and gtars offers no rollback, so a source file discovered
missing *after* the removal leaves the store worse than it started.

### `lock` — (re)generate a build lock without rebuilding

    gks-refgetstore lock --from-log runs/YYYY-MM-DD-build/logs/build.log

Reconstructs the file→collection mapping from a build log and re-hashes the
cache, producing the same lock a full build would. Use it to backfill a lock
for a build that already ran.

### Post-build store self-check (separate entry point)

    gks-refgetstore-check                     # gtars-only checks on ./store
    gks-refgetstore-check PATH                # explicit store dir

## Build lock & verification

A full build writes **`build.lock.json`** (schema
`gks-refgetstore-build-lock/4`). It is small, committed, and has exactly two
jobs, which is what the shape is organized around:

```json
{
  "schema": "gks-refgetstore-build-lock/4",
  "build":   {"timestamp_utc": …, "git": {…}, "gtars_version": "0.11.0",
              "gtars_source": {"url", "vcs", "commit", "subdirectory"} | null},
  "inputs":  {"sources_toml_sha256": …,
              "files": [{"kind","owner","url","cache_path","mutable",
                         "sha256","bytes","present","provider_checksum",
                         "provider_checksum_algorithm",
                         "provider_checksum_blocks","checksum_url",
                         "file_class","ingest_spec","ingest_spec_sha256"}]},
  "outputs": {"n_sequences": 1779497, "n_collections": 240,
              "sequences_root": "sha256:…", "collections_root": "sha256:…",
              "collections": [{"digest": …, "n_sequences": …,
                               "from": ["cache/path", …]}]}
}
```

- **`inputs`** is upstream reproducibility: every concrete file the build read,
  and `ingest_spec` pinning how those bytes were interpreted. The `sha256` pins
  the bytes; the spec pins what was done with them.
- **`outputs`** is local completeness: a **census enumerated from the store**,
  not synthesized from build provenance. Provenance only supplies the `from`
  annotation, so a collection with no known contributor is recorded honestly as
  `"from": []` rather than omitted.

**`outputs.collections[].from` is a list.** One file yields exactly one
collection, but one collection can have many contributing files —
content-addressing merges byte-identical sources, and filtering collapsed eight
Ensembl collections into three. In the committed lock 43 collections have two or
more contributors and one has twelve. `validate_lock` enforces the direction
that *does* hold by rejecting a cache path attributed to two collections.

**Roots** are canonical, and documented so a reimplementation cannot drift:

```python
def digest_root(digests):
    h = hashlib.sha256()
    for d in sorted(set(digests)):     # bare digests, no SQ. prefix
        h.update(d.encode("ascii")); h.update(b"\n")
    return "sha256:" + h.hexdigest()
```

Order- and duplicate-independent by construction, so a root answers exactly one
question — does the store hold the same *set* the lock recorded — and any single
change flips it. Computing both costs about four seconds on the production
store, which is why it happens on every lock write.

Readers go through the accessors (`lock_files`, `lock_collections`,
`collection_by_file`, `files_by_collection`), never raw subscripting, and
through `load_lock`, which validates. Schemas /1–/3 are **rejected outright**,
not shimmed: a /2 lock's `sources` list reads as an empty `inputs.files` under
the v4 accessors, so a tolerant reader would report zero drift and an empty
store census as a pass. The committed lock was converted once rather than
rebuilt — see
[`runs/2026-09-12-lock-schema-v4/`](runs/2026-09-12-lock-schema-v4/) for why a
rebuild would have destroyed test material.

The committed `build.lock.json` is the lock for the **most recent published
build** — a single, latest snapshot, not a history. A periodic process builds
the dataset, uploads the store to object storage, and commits the regenerated
lock here, so the committed lock always corresponds to the currently-published
store and can be referenced from either place. (It may lag the manifest between
runs — after editing `sources.toml`, the lock is refreshed on the next full
build; `verify` handles the interim gap gracefully.)

**Two uses:**

1. **Reproducibility / drift detection** — a live build resolves dynamic feeds,
   downloads and validates all inputs, then compares URL membership, provider
   checksum values, and cached SHA-256 values before opening the store. Any
   change is fatal unless `--force-lock` explicitly accepts a new baseline.
   NCBI MD5 and Ensembl BSD checksums verify provider delivery; the lock's
   SHA-256 identifies the exact bytes consumed by the build.
2. **Provenance** — the lock records each collection's contributing files and
   the store records each collection's membership, so composing the two answers
   file↔digest questions on demand without storing a table:

       uv run python -m gks_refgetstore.provenance --file human.6.rna    # digests from that file
       uv run python -m gks_refgetstore.provenance --digest <sha512t24u> # source file(s) for a digest
       uv run python -m gks_refgetstore.provenance --list                # all collections + contributors

   `--digest` is an O(all collections) scan, ~56s on the production store.
   That is deliberately not backed by a cached reverse index: 1.8M entries to
   build, persist, and keep honest, to save under a minute on an interactive
   query.

3. **Local completeness** — `verify --cache --store` checks that the download
   cache *and the built store* still contain everything the lock says they do.
   Before schema /4 this was entirely unimplemented for the store: the lock
   recorded collection digests and counts, and nothing ever checked them, so a
   deleted `.seq` file was invisible.

**Recommended workflows:**

    # Live discovery with strict drift detection (default)
    gks-refgetstore build

    # Review a legitimate upstream change, then explicitly rebaseline
    gks-refgetstore build --force-lock

    # Rebuild offline from exactly the locked, already-cached bytes
    gks-refgetstore build --locked-sources

**The lock wins over the provider.** When a cached file's *provider* checksum no
longer matches but its bytes still match the lock's `sha256`, the cache is kept
and the drift is reported rather than silently re-downloaded. For a lock-based
tool that is the only defensible default: upstream has since changed the MD5 of
all 30 remaining `mRNA_Prot` shards, so re-downloading would destroy the only
surviving copy of what the lock describes. `--force-download` remains the
override.

**`--lock-check-mode`** controls whether a live build compares against a lock:

| mode | requires | on discrepancy |
| --- | --- | --- |
| `strict` (default) | URL membership, upstream MD5, and cached SHA-256 match | discrepancy stops before ingestion |
| `subset` | retained for CLI compatibility; source drift is still fatal | discrepancy stops before ingestion |
| `ignore` | — | no check |

`--force-lock` is the explicit re-baseline operation. `--locked-sources` is
stricter and offline: it performs no manifest request, and rejects missing or
SHA-256-mismatched cached files.

**Reproducible rebuild from a validated cache:**

    gks-refgetstore fetch --locked-sources          # confirm every locked cache file
    gks-refgetstore build --locked-sources          # rebuild with no live resolution

## Config reference

Each `[[assembly]]` block:

| field | required | description |
| --- | --- | --- |
| `namespace` | yes | Sequence-alias namespace to populate (e.g. `GRCh38`, `GRCh38.p14`). Becomes the namespace of `GRCh38:chr1`-style aliases. |
| `fasta_url` | when `load_fasta=true` | URL of NCBI `*_genomic.fna.gz`. gtars ingests `.gz` directly. |
| `report_url` | yes | URL of the corresponding `*_assembly_report.txt`. |
| `load_fasta` | no, default `true` | If false, do not ingest a FASTA; the assembly report adds aliases only for digests already in the store. Use this for a report-only release or a deliberate alias-only namespace. |
| `fasta_path` | no | Local path overriding the downloaded FASTA, relative to the directory holding this manifest. |
| `checksum_manifest_url` | no | NCBI directory `md5checksums.txt`; genomic FASTA MD5 is required when configured. Some assembly reports have no provider MD5 and are pinned by lock SHA-256. |

Each `[[seqset]]` block (flat FASTA where the header name is the accession):

| field | required | description |
| --- | --- | --- |
| `name` | yes | Human id for logs + CLI selection (`--seqset NAME`). |
| `namespace` | yes | Alias namespace for every ingested sequence (e.g. `refseq`, `ensembl`). |
| `url_template` | one of | URL with `{shard}` placeholder, or plain URL if unsharded. |
| `urls` | one of | Explicit list of URLs (alternative to `url_template`; not valid with `shard_range`). |
| `url_pattern` | one of | HTTPS URL whose basename contains a glob; requires `checksum_manifest_url`. |
| `checksum_manifest_url` | with `url_pattern` | HTTPS NCBI `*.files.installed` URL in the same directory. |
| `md5_manifest_urls` | no | Per-URL NCBI `md5checksums.txt` list for explicit historical inputs. |
| `checksum_manifest_urls` | no | Per-URL Ensembl `CHECKSUMS` list; BSD checksum and 1 KiB block count are validated. |
| `file_class` | no | One source class applied to every resolved file (useful for dynamic patterns). |
| `file_classes` | no | Per-URL source classes used for validation and contribution analysis; mutually exclusive with `file_class`. |
| `release` | release groups | Integer provider release. Release groups are processed numerically oldest-first. |
| `rolling_namespace` | with `release` | Namespace replaced in full after a release succeeds; the immutable namespace must be `<rolling_namespace>-<release>`. |
| `shard_range` | no | `[min, max]` inclusive substituted into `{shard}` in `url_template`. |
| `format` | no, default `fasta` | `fasta`; `gbff` to convert a GenBank flat file to FASTA before ingest; or `lrg_zip` to concatenate the `LRG_N.fasta` members of an EBI LRG bundle, in natural-sorted member order, into one FASTA. Non-FASTA formats are cached as `<artifact>.fasta` beside the download and regenerate offline. |
| `exclude` | no | Table `{file_classes = [...], record_prefixes = [...]}` dropping records from the named file classes before ingest. Matching is on the record name only, never on sequence content. The filtered FASTA is cached as `<artifact>.filtered.fa.gz` beside the download, with `<artifact>.excluded.tsv` listing every dropped record. A rule matching nothing is an error. |

Exactly one of `url_template`, `urls`, or `url_pattern` is required.
`shard_range` is valid only with `url_template`; `checksum_manifest_url` is
valid and required only with `url_pattern`. Pattern matching is local to safe
manifest basenames. Records must be strict `MD5 filename` lines; malformed,
duplicate, path-like, or empty results are errors. Matches are naturally sorted
(`2` before `10`) and each shared manifest is fetched only once per command.

For example, the mutable current RNA feed uses:

```toml
url_pattern = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz"
checksum_manifest_url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.files.installed"
```

NCBI documents `*.files.installed` as the installed file list with MD5 values
in the [species-specific RefSeq README](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/README).
The [RefSeqGene README](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/README.txt)
describes its numbered subsets and frequent updates. Current RefSeq feeds are
therefore intentionally dynamic; historical RefSeq releases, assemblies, and
Ensembl releases remain explicit URLs.

## What gets written

Per assembly run, the loader adds:

- One sequence collection per `load_fasta=true` assembly, named by its
  collection sha512t24u digest.
- Per-sequence aliases under `<namespace>` (the assembly namespace):
  UCSC-style-name (`chr1`, `chrM`, …) and GenBank-Accn (`CM000663.2`, …).
- Per-sequence aliases under the cross-assembly `refseq` namespace
  (`NC_000001.11` → digest) and `insdc` namespace (`CM000663.2` → digest),
  emitted once per digest.
- Per-collection aliases under `refseq` (`GCF_000001405.40`) and `insdc`
  (`GCA_000001405.29`), parsed from the assembly report header.

Per seqset run, the loader adds:

- One sequence collection per FASTA file ingested.
- Per-sequence aliases under `<namespace>` using each FASTA header's first
  token (e.g. `refseq:NM_000551.3`, `ensembl:ENST00000256474.3`,
  `ensembl:ENSP00000256474.3`).
- For the `lrg_zip` seqset, one extra alias per genomic record: the ingested
  header is EBI's record id `lrg:LRG_1g`, and `lrg:LRG_1` is emitted alongside
  it because HGVS and ClinVar spell the genomic record bare (`LRG_1:g.5596del`)
  and never `LRG_1g`. Transcript and protein records (`LRG_1t1`, `LRG_1p1`) are
  already spelled identically upstream and downstream, so they get one alias.

Release-scoped seqsets accumulate every configured file class before publishing
aliases. Conflicting alias-to-digest mappings within a release fail the build.
Immutable `ensembl-N` namespaces are retained, while the complete `ensembl`
TSV is atomically replaced after each release, removing retired aliases.

Aliases for `sha512t24u:…` and `ga4gh:SQ.…` are intentionally NOT written —
they are the raw digest with a prefix and are synthesized at query time by the
consuming alias proxy.

## Version coverage

The manifest keeps release-specific sources explicit while discovering the
membership of mutable current RefSeq feeds from official checksum manifests. It includes:

- **NCBI assemblies:** the original GRCh38 and GRCh37 releases, their listed
  patch releases, and their matching assembly reports. GRCh37.p11 and p12 are
  report-only because NCBI does not publish genomic FASTAs for them.
- **NCBI RefSeq:** current mRNA, protein, and RefSeqGene shards, plus every RNA
  and protein FASTA discoverable in the official GRCh37/GRCh38 assembly and
  taxid-9606 annotation-release histories.
- **Ensembl:** releases **75–116**, each with unmasked `dna.toplevel`, cDNA,
  ncRNA, and peptide FASTAs. CDS, masked DNA, and redundant chromosome or
  primary-assembly subsets are deliberately excluded.
- **EBI LRG:** the complete public Locus Reference Genomic set — 1,325 genomic,
  1,634 transcript, and 1,624 protein records across 1,325 loci — from the single
  aggregate `LRG_public_fasta_files.zip`. The set is frozen at release 810
  (2021-03-31); `pending/`, `stalled/`, and `suppressed/` records are excluded.
  Most LRG sequences are byte-identical to their RefSeq counterparts
  (`LRG_1g` == `NG_007400.1`) and are deduplicated by digest, so this adds
  ~442 new sequences and 5,905 `lrg` aliases (4,580 distinct record ids — LRG_321
  repeats two protein ids across transcripts upstream — plus 1,325 bare genomic
  ids).

### Verifying backwards-compatibility against a seqrepo snapshot

`seqrepo_equivalence/` holds an optional checker that compares a built store
against a biocommons seqrepo snapshot and categorizes every difference into a
report + a PASS/FAIL verdict. See `seqrepo_equivalence/README.md`.

### Known digest divergence on Ensembl proteins with `*`

For ~0.1% of Ensembl protein sequences (115 ENSPs containing embedded stop
codon `*` characters), this store and seqrepo produce different `ga4gh:SQ.*`
identifiers **despite storing identical sequence bytes**. This is caused by
`bioutils.sequences.normalize_sequence()` stripping `*` before computing the
digest in seqrepo, while gtars computes `sha512t24u` directly from the raw
FASTA bytes. The store's digest is correct per the GA4GH VRS specification.
DNA sequences are not affected.

## Published store, experiments, and run records

The current published artifact is the preserved `store.2026-07-22/` build:
1,213,617 sequences and 107 collections, uploaded at
`theferrit32-public:theferrit32-public/refgetstore/2026-07-22`. Its authoritative
manifest remains inside that store; the compact audit record is
[`runs/2026-07-22-published-store/`](runs/2026-07-22-published-store/).

The local `store/` is a build target and is not published automatically. The
prior isolated release-116 experiment is documented at
[`runs/2026-07-23-ensembl-r116-genomic-experiment/`](runs/2026-07-23-ensembl-r116-genomic-experiment/).
Do not infer publication status from a local store directory.

Use [`RUNBOOK.md`](RUNBOOK.md) for new reproducible runs and browse
[`runs/`](runs/) for retained historical evidence.
