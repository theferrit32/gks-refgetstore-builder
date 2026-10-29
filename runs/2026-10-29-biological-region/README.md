# Add RefSeq genomic biological regions

NCBI's RefSeq Functional Elements file of genomic biological regions,
`refseq/H_sapiens/biological_region/human.biological_region.fna.gz`, is now the
`refseq_human_biological_region` seqset in `sources.toml` and
`sources.genomic.toml`. Its 128,261 `NG_` records (recombination regions,
enhancers, pseudogene regions and other non-RefSeqGene regions) are aliased in
the `refseq` namespace. The lock pins the file; it was added to clones of the
full and genomic stores, `/Volumes/RefgetStores/full-2026-10-29` and
`/Volumes/RefgetStores/genomic-2026-10-29`.

| | genomic 2026-10-27 | genomic 2026-10-29 | full (`store/`) | full-2026-10-29 |
| --- | ---: | ---: | ---: | ---: |
| sequences | 13,615 | 141,875 | 1,781,660 | 1,909,920 |
| collections | 33 | 34 | 238 | 239 |
| sequence aliases | 45,471 | 173,732 | — | — |
| `refseq` aliases | 10,437 | 138,698 | — | +128,261 |
| bytes | 4,886,057,165 | 4,962,291,943 | — | — |

The new collection, `G4BWjptZHKGp0F_1ST4DI_T6DIByO72U`, holds 128,261 records:
128,260 distinct sequences (88,477,522 residues), none already in either
store, and two records with identical sequence. Every record name is a new
`refseq` alias.

## The source

The file is the current set, republished in place, with no checksum manifest
in its directory. The lock's SHA-256 of the downloaded file is its integrity
anchor, as for the LRG bundle, and `is_mutable` marks `/biological_region/`
URLs as current endpoints. The cached copy, fetched with `gks-refgetstore
fetch --only refseq_human_biological_region`, is 29,083,824 bytes, SHA-256
`53529e45…95bb728`, holding 128,261 `NG_` records.

## Steps

| label | what |
| --- | --- |
| `fetch` | populate the download cache |
| `clone-genomic`, `clone-full` | APFS clones (`cp -c -R`) of `genomic-2026-10-27` and `store/` |
| `verify-full-before` | the full clone matches the committed lock (L0 and L1) |
| `build-full` | `build --store-dir full-2026-10-29 --seqset refseq_human_biological_region --force-lock`: ingests the file and merges it into `build.lock.json` |
| `verify-full` | `verify --cache --store` on the full clone against the new lock: 341 files ok, 1,909,920 sequences and 239 collections reachable |
| `build-genomic` | `build --config sources.genomic.toml --locked-sources --lock-check-mode subset --no-lock --seqset refseq_human_biological_region` on the genomic clone, the way the genomic store was built |
| `subset-check` | the genomic clone holds exactly the 34 collections the lock attributes to the genomic manifest's 30 sources, and all 141,875 payloads are byte-identical to the full clone's |
| `deep-verify-genomic` | [`artifacts/deep_verify_store.py`](artifacts/deep_verify_store.py): all 141,875 sequences re-digest from the bytes the store serves |
| `tests` | 283 passed |

`build.lock.json` describes `full-2026-10-29`. A lock's `outputs` is a census
of the store it was written against, so the selective build ran on a clone of
the full store rather than on `store/`. `store/` itself is unchanged and is now
one collection behind the lock; `full-2026-10-29` is the full store the lock
describes. `--force-lock` was required because the new file had no lock
baseline.

The published `genomic-2026-10-27` directory is unchanged: its `rgstore.json`
and `sequences.rgsi` hash as they did before the clone.

The lock records the build's code as `1a197d9-dirty`: the manifest and
`is_mutable` changes were uncommitted when it ran, and are committed with this
record.
