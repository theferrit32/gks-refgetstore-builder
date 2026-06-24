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
  ncrna / pep FASTAs (release is pinned in `assemblies.toml`).

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
- Disk:
  - ~5 GB for one GRCh38 major + GRCh38.p14 assemblies (the `.fna.gz` files
    are ~900 MB each; the encoded store is ~2–3 GB).
  - ~3-5 GB additional for the RefSeq mRNA+protein shards (64 files,
    ~50-100 MB each).
  - ~500 MB additional for the Ensembl cdna + ncrna + pep FASTAs.

## Setup

    uv sync

## Layout

    assemblies.toml     # declarative source list (edit to add patches/releases)
    build_store.py      # loader entry point
    verify_store.py     # post-build self-checks (gtars-only, no vrs-python)
    download.sh         # optional pre-fetch of assembly FASTAs/reports into downloads/
    downloads/          # cached FASTA + assembly_report.txt (gitignored)
    store/              # output RefgetStore (gitignored)

## Usage

Build everything in `assemblies.toml`:

    uv run python build_store.py

Build only one assembly:

    uv run python build_store.py --assembly GRCh38

Build only one seqset:

    uv run python build_store.py --seqset refseq_human_rna

Skip categories entirely:

    uv run python build_store.py --skip-seqsets
    uv run python build_store.py --skip-assemblies

All CLI flags:

    --config PATH         assemblies.toml (default: ./assemblies.toml)
    --store-dir PATH      output RefgetStore (default: ./store)
    --download-dir PATH   download cache (default: ./downloads)
    --assembly NAME       only process this assembly namespace
    --seqset NAME         only process this seqset name
    --skip-assemblies     skip all [[assembly]] entries
    --skip-seqsets        skip all [[seqset]] entries
    --force-download      re-fetch even if cached files exist

Verify a built store:

    uv run python verify_store.py            # checks ./store
    uv run python verify_store.py PATH        # checks an explicit store dir

## Config reference

Each `[[assembly]]` block:

| field | required | description |
| --- | --- | --- |
| `namespace` | yes | Sequence-alias namespace to populate (e.g. `GRCh38`, `GRCh38.p14`). Becomes the namespace of `GRCh38:chr1`-style aliases. |
| `fasta_url` | yes | URL of NCBI `*_genomic.fna.gz`. gtars ingests `.gz` directly. |
| `report_url` | yes | URL of the corresponding `*_assembly_report.txt`. |
| `load_fasta` | no, default `true` | If false, aliases are added only for digests already present in the store (cheap patch fanout). |
| `fasta_path` | no | Local path overriding the downloaded FASTA (relative to repo root). |

Each `[[seqset]]` block (flat FASTA where the header name is the accession):

| field | required | description |
| --- | --- | --- |
| `name` | yes | Human id for logs + CLI selection. |
| `namespace` | yes | Alias namespace for every ingested sequence (e.g. `refseq`, `ensembl`). |
| `url_template` | yes | URL with `{shard}` placeholder, or plain URL if unsharded. |
| `shard_range` | no | `[min, max]` inclusive substituted into `{shard}`. Omit for a single-file seqset. |

The same schema handles both the NCBI RefSeq mRNA/Prot shards (sharded
with a `{shard}` placeholder) and the Ensembl cdna / ncrna / pep
releases (single-file, `shard_range` omitted).

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

## Version-drift policy

The store pins **one version** of each source:

- **NCBI assemblies**: pinned by GCF accession (e.g. `GCF_000001405.40` for
  GRCh38.p14). Each patch is a superset of the prior, so only the latest patch
  of each major assembly is loaded.
- **NCBI RefSeq transcripts/proteins**: the FTP shards publish only the
  current version of each accession; older versions are not available.
- **Ensembl**: pinned to a single release (currently **r113**, set in
  `assemblies.toml`). Only the latest version of each ENST/ENSP in that
  release is included.

This differs from seqrepo's accumulative model, which retains every version it
has ever loaded. Callers that depend on older transcript/protein versions will
encounter `KeyError` on accessions this store has not loaded.

### Known digest divergence on Ensembl proteins with `*`

For ~0.1% of Ensembl protein sequences (115 ENSPs containing embedded stop
codon `*` characters), this store and seqrepo produce different `ga4gh:SQ.*`
identifiers **despite storing identical sequence bytes**. This is caused by
`bioutils.sequences.normalize_sequence()` stripping `*` before computing the
digest in seqrepo, while gtars computes `sha512t24u` directly from the raw
FASTA bytes. The store's digest is correct per the GA4GH VRS specification.
DNA sequences are not affected.
