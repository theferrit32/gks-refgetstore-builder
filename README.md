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
- **Ensembl transcripts + proteins**: `ensembl:ENST00000256474.3`,
  `ensembl:ENSP00000256474.3`, etc., loaded from Ensembl release cdna /
  ncrna / pep FASTAs (release is pinned in `sources.toml`).

No dependency on `bioutils` or `biocommons.seqrepo` — the loader reads NCBI
assembly reports directly and derives transcript/protein aliases from FASTA
header names. The only runtime dependency is `gtars`.

> Extracted from the `misc/refgetstore/` directory of vrs-python. The
> functional-equivalence harnesses that compare a built store against seqrepo
> through vrs-python's dataproxy remain in vrs-python, since they depend on it.

## Prerequisites

- Python >= 3.11 (uses `tomllib`).
- [`uv`](https://docs.astral.sh/uv/) for environment + dependency management.
- `gtars` (installed via `uv sync`).
- Network access to `ftp.ncbi.nlm.nih.gov` and `ftp.ensembl.org` (or
  pre-populated cache dirs).
- Disk: the default manifest includes 22 NCBI assembly releases, current
  RefSeq and Ensembl files, and selected historical RefSeq/Ensembl archives.
  Its cache and store therefore need substantially more capacity than a
  current-only build. Use `gks-refgetstore fetch --dry-run` to enumerate the
  selected files before provisioning storage.

## Setup

    uv sync

## Layout

    sources.toml                   # declarative, authoritative source manifest
    cli.py                         # gks-refgetstore CLI (build/verify/fetch/lock)
    build_store.py                 # build engine + shared helpers (config, cache paths, ingest)
    build_lock.py                  # build-lock + cache-verification core (used by build & verify)
    fetch_sources.py               # cache pre-fetch (fetch subcommand)
    provenance.py                  # file <-> refget-digest queries from the lock + store
    build.lock.json                # provenance lock for the most recent build (see below)
    verify_store.py                # post-build gtars self-checks (standalone)
    seqrepo_equivalence/           # optional backwards-compat check vs a seqrepo snapshot
    downloads/                     # cached FASTA + assembly_report.txt (gitignored)
    store/                         # output RefgetStore (gitignored)

## CLI

Everything runs through one CLI, `gks-refgetstore`, with four subcommands.
Two equivalent invocation styles:

    gks-refgetstore <cmd> ...     # installed console script (after `uv sync`)
    uv run cli.py <cmd> ...       # run the source directly, no install

Global: `-v/--verbose` for debug logging. Every subcommand shares
`--config` (default `./sources.toml`) and `--cache-dir` (default `./downloads`).
`gks-refgetstore <cmd> --help` prints the full flag list.

### `build` — ingest the manifest into a RefgetStore

    gks-refgetstore build                          # build the whole manifest
    gks-refgetstore build --assembly GRCh38        # only one assembly (skips seqsets)
    gks-refgetstore build --seqset refseq_human_rna  # only one seqset (skips assemblies)
    gks-refgetstore build --skip-seqsets           # assemblies only

Reads the cache (never re-downloads unless `--force-download`); a missing file
is fetched on demand. Writes/refreshes the build lock at the end (see
[Build lock](#build-lock--cache-verification)).

    --store-dir PATH       output RefgetStore (default: ./store)
    --assembly NAME        only this assembly namespace; --seqset NAME only this seqset
    --skip-assemblies / --skip-seqsets
    --force-download       re-fetch even if cached
    --lock PATH            build-lock to check against + write (default ./build.lock.json)
    --lock-check-mode {strict,subset,ignore}   pre-flight check vs the lock (default subset)
    --no-lock              don't write the lock (the pre-flight check still runs)
    --force-lock           write the lock even if a discrepancy would suppress it

Building into an existing store **appends** (dedup by digest); to rebuild from
scratch, delete the store dir first.

### `fetch` — pre-populate the download cache (no ingest)

    gks-refgetstore fetch                # download every manifest source
    gks-refgetstore fetch --dry-run      # list what would be pulled
    gks-refgetstore fetch --only ensembl_human_cdna   # just one source's files
    gks-refgetstore fetch --kinds assembly_report     # just the small reports

Idempotent (skips files already cached), atomic, disk-guarded
(`--min-free-gb`). Useful to warm the cache before a long build, or offline.

### `verify` — check the cache

Three checks, from cheapest to strictest:

    gks-refgetstore verify                          # (A) integrity: present + gzip -t
    gks-refgetstore verify --check-remote           #     + local size vs server Content-Length
    gks-refgetstore verify --lock build.lock.json   # (B) per-file drift vs the lock
    gks-refgetstore verify --lock build.lock.json --strict-set  # (C) + identical file set

- **(A) integrity** needs no lock: every manifest file is present and, for
  `.gz`, passes `gzip -t` (catches truncation/corruption).
- **(B) per-file drift** compares, for each file the build would ingest that
  *also* appears in the lock, its cached `sha256` against the pinned value —
  keyed per file, so it works even when the manifest has more/fewer files than
  the lock. Files with no lock baseline fall back to a gzip check.
- **(C) set match** additionally fails if the build's file set differs from the
  lock's.

Exit code is non-zero on any failure, so `verify` is CI-friendly.

### `lock` — (re)generate a build lock without rebuilding

    gks-refgetstore lock --from-log build_full.log   # -> ./build.lock.json

Reconstructs the file→collection mapping from a build log and re-hashes the
cache, producing the same lock a full build would. Use it to backfill a lock
for a build that already ran.

### Post-build store self-check (separate script)

    uv run python verify_store.py            # gtars-only checks on ./store
    uv run python verify_store.py PATH        # explicit store dir

## Build lock & cache verification

A full build writes **`build.lock.json`** — a provenance record of exactly what
that build consumed: every source's URL, cached path, byte length, `sha256`, and
the collection digest it produced, plus build metadata (timestamp, git commit,
gtars version, `sources.toml` hash, store counts). It is small and committed.

The committed `build.lock.json` is the lock for the **most recent published
build** — a single, latest snapshot, not a history. A periodic process builds
the dataset, uploads the store to object storage, and commits the regenerated
lock here, so the committed lock always corresponds to the currently-published
store and can be referenced from either place. (It may lag the manifest between
runs — after editing `sources.toml`, the lock is refreshed on the next full
build; `verify` handles the interim gap gracefully.)

**Two uses:**

1. **Reproducibility / drift detection** — `build`'s pre-flight check and the
   `verify --lock` command both compare the cache against the pinned hashes. The
   write is gated on *content drift*: a build never silently overwrites a
   known-good lock hash with a changed one.
2. **Provenance** — because each ingested file becomes exactly one collection and
   the store records collection membership, the lock's `url → collection_digest`
   is enough to answer file↔digest questions on demand:

       uv run python provenance.py --file human.6.rna   # digests from that file
       uv run python provenance.py --digest <sha512t24u> # source file(s) for a digest
       uv run python provenance.py --list                # all sources + collections

**`--lock-check-mode`** (on `build`) controls the pre-flight gate:

| mode | requires | on discrepancy |
| --- | --- | --- |
| `subset` (default) | shared files unchanged (no content drift) | drift → warn + don't write the lock |
| `strict` | shared files unchanged **and** identical file set | drift or set change → warn + don't write |
| `ignore` | — | no check |

`--force-lock` overrides the write suppression (a deliberate re-baseline, e.g.
after a "current" source legitimately updated upstream).

**Reproducible rebuild from a validated cache:**

    gks-refgetstore verify --lock build.lock.json   # confirm cache matches the lock
    rm -rf ./store && gks-refgetstore build          # rebuild from the validated cache

## Config reference

Each `[[assembly]]` block:

| field | required | description |
| --- | --- | --- |
| `namespace` | yes | Sequence-alias namespace to populate (e.g. `GRCh38`, `GRCh38.p14`). Becomes the namespace of `GRCh38:chr1`-style aliases. |
| `fasta_url` | when `load_fasta=true` | URL of NCBI `*_genomic.fna.gz`. gtars ingests `.gz` directly. |
| `report_url` | yes | URL of the corresponding `*_assembly_report.txt`. |
| `load_fasta` | no, default `true` | If false, do not ingest a FASTA; the assembly report adds aliases only for digests already in the store. Use this for a report-only release or a deliberate alias-only namespace. |
| `fasta_path` | no | Local path overriding the downloaded FASTA (relative to repo root). |

Each `[[seqset]]` block (flat FASTA where the header name is the accession):

| field | required | description |
| --- | --- | --- |
| `name` | yes | Human id for logs + CLI selection (`--seqset NAME`). |
| `namespace` | yes | Alias namespace for every ingested sequence (e.g. `refseq`, `ensembl`). |
| `url_template` | one of | URL with `{shard}` placeholder, or plain URL if unsharded. |
| `urls` | one of | Explicit list of URLs (alternative to `url_template`; not valid with `shard_range`). |
| `shard_range` | no | `[min, max]` inclusive substituted into `{shard}` in `url_template`. |
| `format` | no, default `fasta` | `fasta`, or `gbff` to convert a GenBank flat file to FASTA before ingest. |

Exactly one of `url_template` or `urls` is required. The same schema handles
the NCBI RefSeq mRNA/Prot shards (sharded via `{shard}`) and the Ensembl
cdna / ncrna / pep releases (single-file, `shard_range` omitted).

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

Aliases for `sha512t24u:…` and `ga4gh:SQ.…` are intentionally NOT written —
they are the raw digest with a prefix and are synthesized at query time by the
consuming alias proxy.

## Version coverage

The manifest pins every source URL to a specific release. It includes:

- **NCBI assemblies:** the original GRCh38 and GRCh37 releases, their listed
  patch releases, and their matching assembly reports. GRCh37.p11 and p12 are
  report-only because NCBI does not publish genomic FASTAs for them.
- **NCBI RefSeq:** current mRNA, protein, and RefSeqGene shards, plus selected
  official per-patch and annotation-release archives for older transcript and
  protein versions.
- **Ensembl:** current release **113** files and selected official historical
  releases for older ENST/ENSP versions.

This is curated historical coverage, not a complete copy of every version ever
published or loaded by seqrepo. A caller can still receive a `KeyError` for an
older accession.version absent from the selected archives. Backfill such a
version from an authoritative per-accession service such as NCBI E-utilities or
INSDC instead of adding an undocumented bulk file.

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
