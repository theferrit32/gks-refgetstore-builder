# Sequence parity: our RefgetStore vs seqrepo 2024-12-20

_Generated 2026-07-02 15:08:37 by `parity_membership.py`. Numbers are computed, not hand-typed._

Companion to `SEQREPO_PARITY.md` (source strategy) and
`verify_seqrepo_equivalence.py` (categorized backwards-compat verdict). This report
is the **exhaustive per-sequence view**; the full membership table is
`parity_by_digest.tsv` (one row per sequence digest).

## Headline

| metric | value |
|---|---|
| seqrepo distinct sequence digests | 1,144,093 |
| our store distinct sequence digests | 1,252,244 |
| union (rows in the TSV) | 1,422,324 |
| **in both** | **974,013** |
| seqrepo-only (the gap) | 170,080 |
| our-store-only (extensions) | 278,231 |
| **seqrepo digest coverage** | **85.134%** |

"Coverage" = fraction of seqrepo's sequences that are present (by digest) in our
store. A digest counts as covered regardless of which alias spelling either side
uses, because the digest IS the sequence.

## What sequence groups we can load

These are the source groups declared in `sources.toml`, all ingested into the
store this build. Alias counts are what actually landed (per store namespace):

| store namespace | aliases |
|---|---|
| `ensembl` | 945,525 |
| `refseq` | 922,649 |
| `GRCh38.p14` | 2,115 |
| `GRCh38.p12` | 1,638 |
| `GRCh38.p11` | 1,607 |
| `GRCh38.p10` | 1,565 |
| `GRCh38.p9` | 1,553 |
| `GRCh38.p8` | 1,537 |
| `GRCh38.p7` | 1,501 |
| `GRCh38.p6` | 1,493 |
| `GRCh38.p5` | 1,489 |
| `GRCh38.p4` | 1,475 |
| `GRCh38.p3` | 1,443 |
| `GRCh38.p2` | 1,427 |
| `GRCh38.p1` | 1,397 |
| `GRCh38` | 1,365 |
| `insdc` | 941 |
| `GRCh37.p13` | 687 |
| `GRCh37.p12` | 637 |
| `GRCh37.p11` | 623 |
| `GRCh37.p10` | 609 |
| `GRCh37.p9` | 555 |
| `GRCh37.p5` | 463 |
| `GRCh37.p2` | 403 |
| `GRCh37` | 276 |

RefSeq (`refseq`) alias breakdown by accession prefix:

| refseq prefix | count |
|---|---|
| `XM_` | 301,074 |
| `XP_` | 212,399 |
| `NM_` | 165,645 |
| `XR_` | 129,316 |
| `NP_` | 74,219 |
| `NR_` | 32,198 |
| `NG_` | 6,844 |
| `NW_` | 477 |
| `NT_` | 415 |
| `NC_` | 49 |
| `YP_` | 13 |

Ensembl (`ensembl`) alias breakdown by accession prefix:

| ensembl prefix | count |
|---|---|
| `ENST` | 696,082 |
| `ENSP` | 249,443 |

Source groups behind these namespaces:
- **RefSeq current transcripts/proteins** — `human.{1..15}.rna` / `.protein`
  shards (NM/NR/XM/XR, NP/XP).
- **RefSeqGene** — `refseqgene.{1..9}.genomic` (NG_).
- **RefSeq historical predicted+curated** — per-patch assembly dirs
  (orig/p2/p5/p7) + annotation-release archive (AR109 family + 110): the older
  XM/XP/XR/NM/NP/NR that current shards no longer carry.
- **RefSeq curated history (GBFF)** — `…knownrefseq_rna.gbff.gz`, the
  replaced/suppressed NM_/NR_ versions available in no bulk FASTA (converted
  in-flight via `format="gbff"`).
- **Ensembl current + historical** — release 113 cdna/ncrna/pep plus releases
  76–112 and GRCh37 r75.
- **Genome assemblies** — GRCh38, GRCh38.p14, GRCh37, GRCh37.p13 (full FASTA),
  plus the patch-ladder alias fanout (GRCh38 p1–p12, GRCh37 select) via
  `assembly_report.txt` only.

## The residual gap — seqrepo-only sequences

Sequences seqrepo has that our store does not (by digest): **170,080** total.
Not all of these are *sourceable* — many seqrepo entries have no biological
accession we could load from FTP. Split by what seqrepo knows them as:

| bucket | count | can we close it? |
|---|---|---|
| **sourceable** (referenced by a current NCBI/Ensembl/GRCh3x alias) | **37,331** | yes — needs eutils backfill or a not-yet-loaded release |
| non-current-assembly-only (only old patch/assembly namespaces: NCBI3x, GRCh37 patches, hs37d5, JRGv1/2, CHM1) | 52 | out of scope — assemblies we deliberately don't load |
| **digest-only** (no biological accession *anywhere* in seqrepo — only MD5/SEGUID/SHA1/VMC) | **132,697** | not from FTP — seqrepo itself records no accession; only obtainable by copying the raw bytes out of seqrepo |

So the **actionable** transcript/protein gap is ~37,331, not the
headline 170,080. The dominant `(none)` bucket below is the digest-only
sequences — opaque entries seqrepo holds without any accession.

Full breakdown by accession prefix (`(none)` = the digest-only bucket):

| prefix | count |
|---|---|
| `(none)` | 132,697 |
| `XR_` | 26,657 |
| `XM_` | 4,699 |
| `NG_` | 2,138 |
| `XP_` | 2,024 |
| `CHR` | 445 |
| `ENST` | 440 |
| `NT_` | 325 |
| `NW_` | 164 |
| `NP_` | 124 |
| `ENSP` | 104 |
| `chr` | 66 |
| `NC_` | 49 |
| `NR_` | 38 |
| `NM_` | 27 |
| `AC_` | 25 |
| `HG1` | 17 |
| `HSCH` | 7 |
| `chrX` | 3 |
| `HG4` | 3 |
| `chrY` | 2 |
| `HG7` | 2 |
| `Y` | 2 |
| `HG3` | 2 |
| `chrM` | 1 |
| `HG6` | 1 |
| `13` | 1 |
| `1` | 1 |
| `HG9` | 1 |
| `12` | 1 |
| `9` | 1 |
| `22` | 1 |
| `16` | 1 |
| `3` | 1 |
| `21` | 1 |
| `7` | 1 |
| `6` | 1 |
| `10` | 1 |
| `X` | 1 |
| `HG5` | 1 |
| `2` | 1 |
| `17` | 1 |
| `hs3` | 1 |
| `U14` | 1 |

Expected sourceable residual (see `SEQREPO_PARITY.md`): the predicted-model tail
(`XR_`/`XM_`, renumbered per annotation release), a handful of fully-suppressed
`NM_`/`NR_`, predicted proteins (`XP_`/`NP_`), and low-count genomic contigs
(`NT_`/`NW_`). Recoverable only via eutils or accepted as documented drift; a
fresh seqrepo build from current FTP would hit the same wall.

## Our extensions — store-only sequences by prefix

Sequences we hold that seqrepo 2024-12-20 does not (newer releases, `insdc`,
GRCh38.p14, the Ensembl superset).

| prefix | count |
|---|---|
| `ENST` | 188,646 |
| `NM_` | 33,897 |
| `XM_` | 27,414 |
| `ENSP` | 10,791 |
| `XR_` | 9,199 |
| `XP_` | 5,858 |
| `NR_` | 1,669 |
| `NP_` | 744 |
| `(none)` | 13 |

## In both — shared sequences by prefix

| prefix | count |
|---|---|
| `ENST` | 270,095 |
| `XM_` | 249,711 |
| `NM_` | 110,626 |
| `XR_` | 106,724 |
| `ENSP` | 99,914 |
| `XP_` | 75,687 |
| `NR_` | 25,795 |
| `NP_` | 17,887 |
| `(none)` | 9,789 |
| `NG_` | 6,843 |
| `HSCH` | 379 |
| `NW_` | 168 |
| `KI2` | 153 |
| `HG1` | 73 |
| `HG2` | 44 |
| `HG9` | 19 |
| `GL0` | 16 |
| `HG3` | 13 |
| `HG7` | 10 |
| `HG4` | 6 |
| `HG5` | 5 |
| `HG8` | 5 |
| `20` | 2 |
| `10` | 2 |
| `6` | 2 |
| `12` | 2 |
| `8` | 2 |
| `11` | 2 |
| `14` | 2 |
| `21` | 2 |
| `22` | 2 |
| `Y` | 2 |
| `2` | 2 |
| `17` | 2 |
| `15` | 2 |
| `18` | 2 |
| `13` | 2 |
| `7` | 2 |
| `9` | 2 |
| `4` | 2 |
| `19` | 2 |
| `1` | 2 |
| `3` | 2 |
| `16` | 2 |
| `5` | 2 |
| `X` | 2 |
| `HG6` | 1 |
| `MT` | 1 |
| `U43` | 1 |

## Known-divergent note

`ensembl_known_divergent.txt` lists **116** Ensembl accessions whose
digests differ from seqrepo by design (mostly `*`-stop-codon normalization). These
are correct-by-spec differences, not gaps: the sequences are present, only the
digest differs, so they surface as distinct digests on each side rather than as
`both`.
