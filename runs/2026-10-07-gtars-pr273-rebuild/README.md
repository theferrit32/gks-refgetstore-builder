# Full rebuild with gtars from databio/gtars#273

`store/` and `build.lock.json` were rebuilt from the validated download cache
with gtars built from
[databio/gtars#273](https://github.com/databio/gtars/pull/273), so that every
stored sequence is returned with the residues that were imported.

## gtars revision

| | |
| --- | --- |
| repository | [`theferrit32/gtars`](https://github.com/theferrit32/gtars), fork of `databio/gtars` |
| branch | `fix/encodings-buildable` |
| commit | `7831e09f596042e3983c7a123c3f6c019ca6c05f` |
| version | `0.11.0` (unreleased) |

The commit is #273's head (`bbe7fd97`) plus one commit adding the
`StorageMode::Zstd` case to the Python and R bindings, which do not compile
without it. That commit is
[`artifacts/0001-gtars-zstd-bindings.patch`](artifacts/0001-gtars-zstd-bindings.patch).
The written lock records the revision as `build.gtars_source`.

## Why

With gtars 0.10 and earlier, 133 sequences in the store were returned with
different residues than were imported: 106 selenoproteins with `U` returned as
`A`, and 27 `dnaio` sequences with `D`/`H` returned as `H`/`V`. That list, as it
stood in `seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv`, is
[`artifacts/prior-gtars_encoding_roundtrip.tsv`](artifacts/prior-gtars_encoding_roundtrip.tsv).
The cause in gtars is written up in
[`issues/gtars-encoder-alphabet/`](../../issues/gtars-encoder-alphabet/README.md).

The 106 protein payloads cannot be corrected on read, because the encoder had
already stored `U` as alanine's code. A rebuild was chosen over rewriting those
payloads in place because it also exercises the full NCBI/Ensembl/LRG pipeline
against the new gtars.

## Inputs

`gks-refgetstore build --locked-sources` took every URL and SHA-256 from the
committed lock and did no source discovery, so nothing was downloaded. Before
the build, `verify --cache` checked all 335 cached files against their locked
SHA-256 ([`artifacts/verify-cache.txt`](artifacts/verify-cache.txt)): 335 ok.

The input file set and every file's SHA-256 are unchanged from the previous
lock, so the upstream drift that lock preserves (republished `mRNA_Prot` shards
and withdrawn files; see `../2026-09-12-lock-schema-v4/`) is preserved here
too. The only input fields that changed are `ingest_spec` and
`ingest_spec_sha256` on 35 files: the current `sources.toml` declares the
`CHR_` exclusion for 34 Ensembl `dna.toplevel` files and the `lrg_zip` format
for the LRG bundle, and the previous lock predated both.

## Commands

| label | exit | what |
| --- | --- | --- |
| `tests` | 0 | 243 passed |
| `build` | 1 | stopped in the preflight filter; see below |
| `tests-2` | 0 | 244 passed, with the fix below |
| `build-2` | 0 | the rebuild: 27 minutes |
| `verify-deep` | 0 | `verify --store --deep` against an empty baseline |

### The first build attempt

`build` stopped about two minutes in, before any store was created, with
`ValueError: exclusion matched no records in Homo_sapiens.GRCh38.cdna.all.fa.gz
(prefixes 'CHR_')`. Under `--locked-sources` a seqset's files arrive in the
lock's order, sorted by cache path, which puts each Ensembl release's `cdna`
file before `dna`. `SeqsetConfig.file_class_for_index` looked the class up by
position in the manifest's `file_classes`, so the `cdna` file was handed the
`dna.toplevel` exclusion, matched nothing, and tripped the guard against a
declared exclusion that removes nothing. It now prefers the class carried by
each resolved source, with a regression test. The 34 partial `.part` filter
outputs the attempt left in `downloads/` were removed.

## Outcome

```
inputs.files          335
outputs.n_sequences   1779052
outputs.n_collections 235
sequences_root        sha256:3ecf0d7221a60fa038aa12baac6b1da646b4917513b176a1d0661b1541ed921c
collections_root      sha256:829d3700161d1f151efc53cddb7e41d97fabc0e032eec5a6d0b93eeebac41f5c
```

Against the previous store (now `store.pre-pr273/`): 445 sequences removed and
none added. All 445 are padded `CHR_` scaffolds from the 8 Ensembl
`dna.toplevel` collections that the exclusion rule replaces with 3 filtered
ones: `dna3bit`, 60,047,447,030 bp of N-padding in total, each with a
true-length counterpart already in the store. The store is 28 GB on disk,
down from 49 GB. `store.pre-filter/` holds exactly the same digest set as this
store.

Verification of the new store:

- `verify --store --deep` against an empty known-bad baseline
  ([`artifacts/empty-known-bad.tsv`](artifacts/empty-known-bad.tsv)): all
  1,779,052 sequences re-digest to their own digest. PASS, 0 errors,
  0 warnings.
- gtars' own `RefgetStore.verify()`
  ([`artifacts/gtars-verify.txt`](artifacts/gtars-verify.txt)): 1,779,052
  checked, 0 failed, in every alphabet.
- All 133 previously affected sequences are present and round-trip.
- The build's 5 warnings are GRCh37.p11/p12 assembly-report rows naming `NW_`
  scaffolds the store does not hold, identical to the 2026-08-26 build.

`seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv` is now
header-only. The written lock is
[`artifacts/build.lock.json`](artifacts/build.lock.json), identical to the
committed `build.lock.json`.

## Reading this store

Payloads use residue codes that gtars 0.10 and earlier do not know, and those
versions decode them as different letters without raising an error
(selenocysteine `U` reads back as `X`). The store format version is unchanged,
so an older gtars still opens the store. Read it with the revision above or a
later release that includes #273.
