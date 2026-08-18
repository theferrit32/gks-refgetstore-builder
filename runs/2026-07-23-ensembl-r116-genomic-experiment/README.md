# Ensembl genomic DNA experiment

This local-only experiment tests whether Ensembl release 116's human genomic
sequences are identical to the GRCh38.p14 sequences already loaded from NCBI,
and what aliases would be added if the Ensembl FASTA were ingested.

The canonical repository configuration (`sources.toml`) is intentionally not
changed by this experiment. A future, separate task will update the current
Ensembl transcript/protein sources to release 116 and move releases 113–115
into the historical lists.

## Safety boundary

- `store.2026-07-22/` is the preserved copy-on-write clone of the published
  store. It contains 1,213,755 files and its root manifest and indexes were
  verified byte-for-byte against the original immediately after cloning.
- `store/` is the experiment/playground store.
- Only the large FASTA/cache remains under ignored
  `downloads.ensembl-experiment/`; this compact record is retained under `runs/`.
- Nothing from this experiment is uploaded to the published R2 prefix.

## Source

Ensembl release 116 human DNA directory:

<https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/>

Artifact under test:

`Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`

The directory's official README identifies the represented assembly as GenBank
assembly `GCA_000001405.29` (GRCh38.p14). We use the unmasked primary assembly:
chromosomes and nonchromosomal primary-assembly sequences, excluding
haplotypes and patches. Ensembl's `CHECKSUMS` entry for this artifact is:

```text
22450 861294 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Those are the traditional Unix `sum` checksum and 1 KiB block count, not an
MD5 or SHA-256 digest. The completed download passed both that provider
checksum and `gzip -t`. Its local SHA-256 is:

```text
d8c3af0094a7bba6125763bad779ec18a81483c739c6ed122094bdf86c187b92
```

## Questions

1. Does every Ensembl primary-assembly record resolve to a sequence digest
   already present in the preserved NCBI-backed store?
2. Do Ensembl chromosome names such as `1`, `X`, and `MT` resolve to the same
   digests as the corresponding GRCh38/RefSeq/INSDC aliases?
3. Which Ensembl records, if any, differ or are absent?
4. What collection and sequence aliases would ingestion add under the
   `ensembl` namespace?
5. Is `ensembl:1` useful enough to justify its cross-species and
   cross-assembly ambiguity, or should genomic aliases use a more specific
   namespace?

## Procedure

1. Download the exact release-116 primary-assembly FASTA into the isolated
   cache and validate the provider checksum.
2. Ingest it into the playground `store/` using
   `sources.ensembl-dna-r116.toml` and a separate experiment lock.
3. Compare every resulting name/digest pair to `store.2026-07-22/`.
4. Report identical, missing, and differing records, plus aliases and store
   objects added.
5. Recheck the preserved store's root hashes after the experiment.

The experiment configuration contains only this one source and is not intended
to replace the canonical manifest.

## Results

The build ingested one new collection,
`oLfPx0NOBKKXMIngGeQ4YewtU4Ge_wKz`, containing 194 records, and added 194
aliases under the `ensembl` namespace.

- 179 of 194 records are byte-identical to sequences in the preserved store.
- All 169 nonchromosomal records match by both name and digest.
- 10 of 25 chromosome records match by both name and digest.
- The other 15 chromosome records have the same lengths but different digests:
  `1`, `2`, `3`, `6`, `7`, `9`, `10`, `12`, `13`, `16`, `17`, `21`, `22`,
  `X`, and `Y`.
- Across those 15 chromosomes there are only 99 differing characters. Every
  difference is an Ensembl `N` where the preserved NCBI/INSDC sequence has an
  IUPAC ambiguity symbol: 2 `B`, 8 `K`, 8 `M`, 27 `R`, 5 `S`, 14 `W`, and
  35 `Y`.
- Consequently, the playground has 15 new sequence digests and one new
  collection (16 additional object files), not 194 new sequences.
- The preserved store remains at 1,213,617 sequences and 107 collections; the
  playground now has 1,213,632 sequences and 108 collections.
- The repository verifier passes 67/67 against both the preserved store and
  the modified playground store.

For example:

```text
ensembl:1          -> 2YnepKM7OkBoOrKmvHbGqguVfF9amCST
GRCh38:1           -> Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO
insdc:CM000663.2   -> Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO
```

Chromosome 1 has the same length in both representations and differs at only
two positions: Ensembl uses `N` where the accessioned chromosome uses `M` and
`R`.

The Ensembl FASTA headers use coordinate names such as `1`; they do not contain
the INSDC accession `CM000663.2`. Ingestion therefore creates `ensembl:1` but
not `ensembl:CM000663.2`. Assigning the latter to the Ensembl-normalized digest
would be incorrect because `CM000663.2` identifies the sequence retaining the
IUPAC ambiguity symbols.

Machine-readable outputs:

- `comparison.tsv`: one row for each of the 194 Ensembl records.
- `comparison-summary.json`: aggregate identity and substitution counts.
- `build.ensembl-dna-r116.lock.json`: source SHA-256 and collection provenance.

The comparison is now location-independent and requires every path explicitly:

    uv run python runs/2026-07-23-ensembl-r116-genomic-experiment/compare.py \
      --playground store \
      --baseline store.2026-07-22 \
      --lock runs/2026-07-23-ensembl-r116-genomic-experiment/build.ensembl-dna-r116.lock.json \
      --out-tsv runs/2026-07-23-ensembl-r116-genomic-experiment/comparison.tsv \
      --out-json runs/2026-07-23-ensembl-r116-genomic-experiment/comparison-summary.json

## Supporting evidence

[`evidence/grc-human-table.png`](evidence/grc-human-table.png) is a screenshot
of the Genome Reference Consortium **Human Genome Assembly GRCh38.p14**
chromosome-length/accession table, captured from the NCBI GRC human assembly
data page: <https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38.p14>. It
supports the chromosome-name-to-GenBank/RefSeq accession comparison; the
machine-readable comparison and build lock remain the primary evidence.
