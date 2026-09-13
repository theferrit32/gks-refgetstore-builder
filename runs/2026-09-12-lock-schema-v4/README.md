# Build-lock schema /2 → /4 conversion

The committed `build.lock.json` was **converted in place**, not rebuilt.

## Why conversion rather than a rebuild

A rebuild re-resolves every source against its provider. Upstream has since
changed the MD5 of all 30 remaining `mRNA_Prot` shards and withdrawn two
(`human.16.rna`, `human.16.protein`), and a ninth `refseqgene` shard has
appeared. The committed lock is the only surviving record of the bytes the
store was actually built from, and that preserved drift is live test material:
it is what makes `verify --manifest` fail while `verify --cache --store`
passes. Rebuilding would have overwritten it.

## What changed

| /2 | /4 |
|---|---|
| `build.sources_toml.{path,sha256}` | `inputs.sources_toml_sha256` |
| `build.store.{n_sequences,n_collections}` (strings) | `outputs.{n_sequences,n_collections}` (integers, enumerated from `store/`) |
| `sources[]` | `inputs.files[]` |
| `sources[].upstream_md5` | dropped — reconstructed from `provider_checksum` + `provider_checksum_algorithm` |
| `sources[].collection_digest`, `sources[].n_sequences` | `outputs.collections[].{digest,n_sequences,from}` |
| — | `outputs.sequences_root`, `outputs.collections_root` |

The absolute `sources_toml.path` was machine-specific noise in a committed
file. The counts were being written as strings by `store.stats()`; the census
uses `len()` over the enumerated digest sets, so they are integers and they
agree with the roots by construction.

## Outcome

```
inputs.files          335
outputs.n_sequences   1779497
outputs.n_collections 240
sequences_root        sha256:c32ff976bf498cde06c85a6067a43ca012a2e29815944e7b3c3c020cd3797700
collections_root      sha256:16e1fb84f30ca5e8924a5b7b05547bae0b43f9374d4a8c3fe966f03136b9305b
collections with no known contributor  0
lock claimed but absent from the store 0
collections with >1 contributor        43 (widest 12)
```

Two numbers are worth reading carefully.

**43 collections have two or more contributors, and one has twelve.** The /2
schema's `sources[].collection_digest` encoded one file → one collection, and
`provenance.py`'s docstring asserted that as a property. It was never true:
content-addressing collapses byte-identical sources into a single collection.
The `from` *list* is the correction. The inverse direction — one file yields
exactly one collection — does hold, and `validate_lock` enforces it by
rejecting a cache path attributed to two collections.

**312 of 335 files are attributed.** The 23 that are not are the
`assembly_report` files, which contribute aliases rather than sequences and so
produce no collection.

## `ingest_spec` is null for every file, deliberately

`store/` was built before the Ensembl padded-scaffold exclusion rule existed,
so nothing in it was transformed at ingest. Writing the *current* manifest's
specs would claim an exclusion had been applied when it had not, and `sync`
would then classify those sources as `unchanged` and skip exactly the
re-ingest the exclusion requires.

## This store is not the filtered one

`store/` has **not** had the padding sync applied: 1,779,497 sequences across
240 collections. `store.pre-filter/` — despite the name — is the newer,
filtered store, at 1,779,052 across 235.

## Reproducing

```sh
uv run python runs/2026-09-12-lock-schema-v4/convert.py \
    --in build.lock.json --store store --out build.lock.json
```

`convert.py` is throwaway and is retained as the record of how the committed
lock was produced, not as a supported migration path. Schema /4 has no shim:
`build_lock.load_lock` rejects /1, /2, and /3 outright.
