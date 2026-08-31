# Sequence parity: our RefgetStore vs seqrepo 2024-12-20

_Generated 2026-08-31 11:45:35 by `parity_membership.py`. Numbers are computed, not hand-typed._

Companion to the run record README (source strategy) and
`seqrepo_equivalence/verify_seqrepo_equivalence.py` (categorized backwards-compat
verdict). This report is the **exhaustive per-sequence view**; the full membership
table is `parity_by_digest.tsv` (one row per sequence digest).

## Headline

| metric | value |
|---|---|
| seqrepo distinct sequence digests | 1,144,093 |
| our store distinct sequence digests | 1,779,497 |
| union (rows in the TSV) | 1,942,926 |
| **in both** | **980,664** |
| seqrepo-only (the gap) | 163,429 |
| our-store-only (extensions) | 798,833 |
| seqrepo digest coverage, all digests | 85.715% |
| **coverage of seqrepo sequences that carry an accession** | **96.522%** |

"Coverage" = fraction of seqrepo's sequences that are present (by digest) in our
store. A digest counts as covered regardless of which alias spelling either side
uses, because the digest IS the sequence.

**Quote the second number, not the first.** 128,092 of the
163,429 missing digests (78%)
are *digest-only*: seqrepo records no biological accession for them anywhere, so
there is no FTP source to load them from and no build could close that gap. The
raw 85.715% therefore reads as far worse than the position is. The
actionable shortfall is 35,286 sequences
(3.1% of seqrepo), and it
is dominated by NCBI's predicted-model tail, which is renumbered every annotation
release. Full decomposition in *The residual gap* below.

## What sequence groups we can load

These are the source groups declared in `sources.toml`, all ingested into the
store this build. Alias counts are what actually landed (per store namespace):

| store namespace | aliases |
|---|---|
| `ensembl` | 1,052,681 |
| `ensembl-116` | 1,052,681 |
| `refseq` | 914,128 |
| `ensembl-115` | 778,858 |
| `ensembl-114` | 535,513 |
| `ensembl-113` | 535,460 |
| `ensembl-112` | 401,263 |
| `ensembl-111` | 399,303 |
| `ensembl-110` | 399,061 |
| `ensembl-109` | 395,831 |
| `ensembl-108` | 395,109 |
| `ensembl-107` | 393,093 |
| `ensembl-106` | 386,872 |
| `ensembl-105` | 384,028 |
| `ensembl-104` | 374,860 |
| `ensembl-103` | 370,945 |
| `ensembl-102` | 366,968 |
| `ensembl-101` | 362,807 |
| `ensembl-100` | 360,291 |
| `ensembl-99` | 360,146 |
| `ensembl-98` | 359,438 |
| `ensembl-97` | 357,419 |
| `ensembl-96` | 338,804 |
| `ensembl-95` | 336,053 |
| `ensembl-94` | 336,053 |
| `ensembl-93` | 331,804 |
| `ensembl-92` | 331,804 |
| `ensembl-90` | 323,538 |
| `ensembl-91` | 323,538 |
| `ensembl-88` | 321,932 |
| `ensembl-89` | 321,932 |
| `ensembl-75` | 320,230 |
| `ensembl-84` | 319,297 |
| `ensembl-86` | 318,276 |
| `ensembl-87` | 318,276 |
| `ensembl-85` | 318,276 |
| `ensembl-83` | 316,819 |
| `ensembl-82` | 314,995 |
| `ensembl-81` | 314,995 |
| `ensembl-79` | 311,348 |
| `ensembl-80` | 311,348 |
| `ensembl-78` | 307,058 |
| `ensembl-77` | 307,058 |
| `ensembl-76` | 306,097 |
| `lrg` | 5,905 |
| `GRCh38.p14` | 2,115 |
| `GRCh38.p13` | 1,872 |
| `GRCh38.p12` | 1,642 |
| `GRCh38.p11` | 1,611 |
| `GRCh38.p10` | 1,569 |
| `GRCh38.p9` | 1,557 |
| `GRCh38.p8` | 1,541 |
| `GRCh38.p7` | 1,505 |
| `GRCh38.p6` | 1,497 |
| `GRCh38.p5` | 1,489 |
| `GRCh38.p4` | 1,475 |
| `GRCh38.p3` | 1,443 |
| `GRCh38.p2` | 1,427 |
| `GRCh38.p1` | 1,397 |
| `GRCh38` | 1,365 |
| `insdc` | 972 |
| `GRCh37.p13` | 687 |
| `GRCh37.p12` | 661 |
| `GRCh37.p11` | 649 |
| `GRCh37.p10` | 643 |
| `GRCh37.p9` | 583 |
| `GRCh37.p5` | 489 |
| `GRCh37.p2` | 419 |
| `GRCh37` | 276 |

RefSeq (`refseq`) alias breakdown by accession prefix:

| refseq prefix | count |
|---|---|
| `XM_` | 307,107 |
| `XP_` | 218,125 |
| `NM_` | 138,031 |
| `XR_` | 134,489 |
| `NP_` | 76,667 |
| `NR_` | 31,880 |
| `NG_` | 6,844 |
| `NW_` | 508 |
| `NT_` | 415 |
| `NC_` | 49 |
| `YP_` | 13 |

Ensembl (`ensembl`) alias breakdown by accession prefix:

| ensembl prefix | count |
|---|---|
| `ENST` | 669,547 |
| `ENSP` | 382,428 |
| `HSCH` | 349 |
| `KI2` | 153 |
| `HG2` | 104 |
| `HG1` | 42 |
| `GL0` | 16 |
| `HG4` | 6 |
| `HG7` | 4 |
| `HG5` | 2 |
| `HG6` | 2 |
| `HG9` | 2 |
| `1` | 1 |
| `10` | 1 |
| `11` | 1 |
| `12` | 1 |
| `13` | 1 |
| `14` | 1 |
| `15` | 1 |
| `16` | 1 |
| `17` | 1 |
| `18` | 1 |
| `19` | 1 |
| `2` | 1 |
| `20` | 1 |
| `21` | 1 |
| `22` | 1 |
| `3` | 1 |
| `4` | 1 |
| `5` | 1 |
| `6` | 1 |
| `7` | 1 |
| `8` | 1 |
| `9` | 1 |
| `HG3` | 1 |
| `MT` | 1 |
| `X` | 1 |
| `Y` | 1 |

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

Sequences seqrepo has that our store does not (by digest): **163,429** total.
Not all of these are *sourceable* — many seqrepo entries have no biological
accession we could load from FTP. Split by what seqrepo knows them as:

| bucket | count | can we close it? |
|---|---|---|
| **sourceable** (referenced by a current NCBI/Ensembl/GRCh3x alias) | **35,286** | yes — needs eutils backfill or a not-yet-loaded release |
| non-current-assembly-only (only old patch/assembly namespaces: NCBI3x, GRCh37 patches, hs37d5, JRGv1/2, CHM1) | 51 | out of scope — assemblies we deliberately don't load |
| **digest-only** (no biological accession *anywhere* in seqrepo — only MD5/SEGUID/SHA1/VMC) | **128,092** | not from FTP — seqrepo itself records no accession; only obtainable by copying the raw bytes out of seqrepo |

So the **actionable** transcript/protein gap is ~35,286, not the
headline 163,429. The dominant `(none)` bucket below is the digest-only
sequences — opaque entries seqrepo holds without any accession.

Full breakdown by accession prefix (`(none)` = the digest-only bucket):

| prefix | count |
|---|---|
| `(none)` | 128,092 |
| `XR_` | 25,203 |
| `XM_` | 4,136 |
| `NG_` | 2,137 |
| `XP_` | 1,662 |
| `NM_` | 1,367 |
| `NT_` | 325 |
| `NW_` | 164 |
| `ENSP` | 82 |
| `NR_` | 72 |
| `chr` | 66 |
| `NC_` | 49 |
| `NP_` | 37 |
| `AC_` | 25 |
| `chrX` | 3 |
| `HG1` | 2 |
| `chrY` | 2 |
| `chrM` | 1 |
| `HG3` | 1 |
| `hs3` | 1 |
| `U43` | 1 |
| `U14` | 1 |

Expected sourceable residual (see the run record README): the predicted-model tail
(`XR_`/`XM_`, renumbered per annotation release), a handful of fully-suppressed
`NM_`/`NR_`, predicted proteins (`XP_`/`NP_`), and low-count genomic contigs
(`NT_`/`NW_`). Recoverable only via eutils or accepted as documented drift; a
fresh seqrepo build from current FTP would hit the same wall.

## Our extensions — store-only sequences by prefix

Sequences we hold that seqrepo 2024-12-20 does not (newer releases, `insdc`,
GRCh38.p14, the Ensembl superset).

| prefix | count |
|---|---|
| `ENST` | 588,918 |
| `ENSP` | 145,731 |
| `XM_` | 31,731 |
| `XR_` | 11,718 |
| `XP_` | 8,815 |
| `NM_` | 8,605 |
| `NR_` | 1,472 |
| `NP_` | 1,411 |
| `LRG` | 430 |
| `Y` | 1 |
| `HSCH` | 1 |

## In both — shared sequences by prefix

| prefix | count |
|---|---|
| `ENST` | 282,701 |
| `XM_` | 250,351 |
| `NM_` | 109,137 |
| `XR_` | 108,186 |
| `ENSP` | 101,612 |
| `XP_` | 76,184 |
| `NR_` | 25,759 |
| `NP_` | 18,211 |
| `NG_` | 6,844 |
| `HSCH` | 465 |
| `CHR` | 445 |
| `NW_` | 168 |
| `HG1` | 156 |
| `KI2` | 153 |
| `HG2` | 53 |
| `HG9` | 38 |
| `HG3` | 26 |
| `HG7` | 22 |
| `GL0` | 16 |
| `HG4` | 15 |
| `HG5` | 11 |
| `HG8` | 10 |
| `LRG` | 8 |
| `3` | 5 |
| `Y` | 5 |
| `10` | 4 |
| `6` | 4 |
| `13` | 4 |
| `12` | 4 |
| `1` | 4 |
| `17` | 4 |
| `21` | 4 |
| `22` | 4 |
| `2` | 4 |
| `9` | 4 |
| `16` | 4 |
| `7` | 4 |
| `X` | 4 |
| `20` | 3 |
| `11` | 3 |
| `HG6` | 3 |
| `8` | 3 |
| `14` | 3 |
| `15` | 3 |
| `18` | 3 |
| `4` | 3 |
| `19` | 3 |
| `5` | 3 |
| `MT` | 1 |

## Known-divergent note

`ensembl_vs_seqrepo_digest_divergence.tsv` lists **116** Ensembl accessions whose
digests differ from seqrepo by design (mostly `*`-stop-codon normalization). These
are correct-by-spec differences, not gaps: the sequences are present, only the
digest differs, so they surface as distinct digests on each side rather than as
`both`.
