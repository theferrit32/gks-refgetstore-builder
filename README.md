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
- `gtars` (installed via `uv sync`).
- Network access to `ftp.ncbi.nlm.nih.gov` and `ftp.ensembl.org` (or
  pre-populated cache dirs).
- Disk: the default manifest includes 23 NCBI assembly releases, complete
  published GRCh37/GRCh38 RefSeq release history, and all four canonical FASTA
  classes for Ensembl releases 75–116.
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
    inventory_sources.py           # exact remote sizes + pre-download space-budget gate
    provenance.py                  # file <-> refget-digest queries from the lock + store
    build.lock.json                # provenance lock for the most recent build (see below)
    verify_store.py                # post-build gtars self-checks (standalone)
    seqrepo_equivalence/           # optional backwards-compat check vs a seqrepo snapshot
    RUNBOOK.md                     # reproducible run-record lifecycle and policies
    runs/                          # compact historical manifests, summaries, logs, evidence
    downloads/                     # cached FASTA + assembly_report.txt (gitignored)
    store/                         # local output/playground RefgetStore (gitignored)
    store.2026-07-22/              # preserved published store (local, gitignored)

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
    --ingest-jobs N        FASTAs imported concurrently per seqset (default: min(8, cores))
    --min-free-gb N        refuse to start if free disk is below this (default 25)
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

    gks-refgetstore lock --from-log runs/YYYY-MM-DD-build/logs/build.log

Reconstructs the file→collection mapping from a build log and re-hashes the
cache, producing the same lock a full build would. Use it to backfill a lock
for a build that already ran.

### Post-build store self-check (separate script)

    uv run python verify_store.py            # gtars-only checks on ./store
    uv run python verify_store.py PATH        # explicit store dir

## Build lock & cache verification

A full build writes **`build.lock.json`** — a provenance record of exactly what
that build consumed: every source's URL, provider checksum and algorithm (when
supplied), checksum source, cached path, byte length, `sha256`, file class, and
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

1. **Reproducibility / drift detection** — a live build resolves dynamic feeds,
   downloads and validates all inputs, then compares URL membership, provider
   checksum values, and cached SHA-256 values before opening the store. Any
   change is fatal unless `--force-lock` explicitly accepts a new baseline.
   NCBI MD5 and Ensembl BSD checksums verify provider delivery; the lock's
   SHA-256 identifies the exact bytes consumed by the build.
2. **Provenance** — because each ingested file becomes exactly one collection and
   the store records collection membership, the lock's `url → collection_digest`
   is enough to answer file↔digest questions on demand:

       uv run python provenance.py --file human.6.rna   # digests from that file
       uv run python provenance.py --digest <sha512t24u> # source file(s) for a digest
       uv run python provenance.py --list                # all sources + collections

**Recommended workflows:**

    # Live discovery with strict drift detection (default)
    gks-refgetstore build

    # Review a legitimate upstream change, then explicitly rebaseline
    gks-refgetstore build --force-lock

    # Rebuild offline from exactly the locked, already-cached bytes
    gks-refgetstore build --locked-sources

The first build after introducing dynamic sources encounters the readable v1
lock and requires `--force-lock` to establish the v2 URL/MD5 baseline.

**`--lock-check-mode`** controls whether a live build compares against a lock:

| mode | requires | on discrepancy |
| --- | --- | --- |
| `strict` (default) | URL membership, upstream MD5, and cached SHA-256 match | discrepancy stops before ingestion |
| `subset` | retained for CLI compatibility; source drift is still fatal | discrepancy stops before ingestion |
| `ignore` | — | no check |

`--force-lock` is the explicit re-baseline operation. `--locked-sources` is
stricter and offline: it performs no manifest request, requires a v2 lock, and
rejects missing or SHA-256-mismatched cached files.

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
| `fasta_path` | no | Local path overriding the downloaded FASTA (relative to repo root). |
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
