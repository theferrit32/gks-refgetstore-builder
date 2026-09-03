# Ensembl's N-padded alt/patch scaffolds

**What.** Ensembl releases 76 through 109 emit every GRCh38 alt locus and patch
scaffold in `dna.toplevel` at **parent-chromosome scale**: the scaffold's real
bases sit at their placed coordinates on the parent chromosome and every other
position is `N`. The record is named with a `CHR_` prefix. From release 110
onward the same scaffolds are emitted at true length under the bare name.

**Scale.** 445 distinct padded records exist across releases 76–109. They total
**60,047,447,030 bases** and carry **167,367,270** bases of actual scaffold — an
aggregate inflation of **358.8x**. That is **56.7% of every base the store
holds**, in 445 of its 1,779,497 sequences.

**Recoverability.** All 445 collapse exactly to their true-length counterpart,
and all 445 true-length counterparts are already in the store independently.
Nothing exists only in padded form.

Two machine-readable tables back every count below; see [Data files](#data-files).

## Alt loci, patches, and what padding means

GRCh38 is not 25 sequences. Besides the 25 assembled molecules (chromosomes
1–22, X, Y, MT) the Genome Reference Consortium publishes:

| role | GRCh38.p14 count | what it is |
|---|---:|---|
| assembled-molecule | 25 | the chromosomes |
| alt-scaffold | 261 | an **alternate haplotype** for a polymorphic region — a second, equally valid representation of a stretch of a chromosome (the MHC, the KIR locus, …) |
| fix-patch | 164 | a **correction** to a region of a chromosome, staged for the next major assembly |
| novel-patch | 90 | **added sequence** for a region, likewise staged |
| unlocalized-scaffold | 42 | sequence known to belong to a named chromosome, position unknown |
| unplaced-scaffold | 127 | sequence known to be human, chromosome unknown |

Alt loci and patches share one property that the last two rows lack: each is an
alternative representation of a **specific, known interval of a specific parent
chromosome**. NCBI publishes that interval in `alt_scaffold_placement.txt`.

Ensembl's padded representation makes the placement implicit in the file
format. A padded record is roughly as long as its parent chromosome and its real
bases begin at the placement's start coordinate, so **offset *i* in the padded
record is offset *i* on the parent chromosome** for everything up to the
scaffold. A tool that knows only "this record is called `CHR_HSCHR1_2_CTG3`" can
line its bases up against chromosome 1 without consulting a placement table. The
cost is that the overwhelming majority of every such record is padding.

## Worked example: `HSCHR1_2_CTG3`

An alt scaffold for the PRAME region of chromosome 1, GenBank `KI270766.1`,
RefSeq `NT_187517.1`.

| | value |
|---|---:|
| `ensembl-100:CHR_HSCHR1_2_CTG3` length | 248,975,002 |
| real (non-`N`) bases | 256,271 — **0.1029%** |
| real span, 0-based half-open | 13,075,112 .. 13,331,383 |
| digest | `aDMwTJZkNZfyltbc6wRiMSb_DbTZEaRV` |
| alphabet in this store | `dna3bit` |
| `ensembl-116:HSCHR1_2_CTG3` length | 256,271 — 100% real |
| digest | `5sJPeOQINUJr0syBJgsAX7bqCNmmWXJq` |
| alphabet in this store | `dna2bit` |
| chromosome 1 length, for scale | 248,956,422 |

The padded record is 18,580 bases **longer** than chromosome 1 — the padded
length is close to the parent's but not equal to it. The exact rule is derived
[below](#the-construction-rule).

Extracting the real span from the padded record and hashing it with sha512t24u
yields exactly `5sJPeOQINUJr0syBJgsAX7bqCNmmWXJq`. The two records are the same
bases; one is wrapped in 248.7 million `N`.

For contrast, `N` in a *chromosome* is a different phenomenon — centromeres,
telomeres and assembly gaps:

| record | total | real | % real |
|---|---:|---:|---:|
| `ensembl-100:1` (chromosome 1) | 248,956,422 | 230,481,012 | 92.58% |
| `ensembl-100:CHR_HSCHR1_2_CTG3` | 248,975,002 | 256,271 | 0.10% |
| `ensembl-100:CHR_HG708_PATCH` | 159,270,839 | 589,656 | 0.37% |

Across all 445 padded records: median 0.159% real, minimum 0.0064%
(`CHR_HG1817_1_PATCH`, 7,309 real bases in 114,369,395), maximum 5.36%
(`CHR_HG2365_PATCH`, 5,500,449 in 102,714,182).

## The padding offset is the NCBI placement

NCBI's placement row for this scaffold, from the GRCh38.p14
`ALT_REF_LOCI_1/alt_scaffolds/alt_scaffold_placement.txt`:

```
alt_scaf_name   HSCHR1_2_CTG3      alt_scaf_acc   NT_187517.1
parent_name     1                  parent_acc     NC_000001.11
region_name     PRAME_REGION_1     ori            +
alt_scaf_start  20633              alt_scaf_stop  256271
parent_start    13075113           parent_stop    13312803
alt_start_tail  20632              alt_stop_tail  0
```

`parent_start` is 1-based: **13,075,113**. The padded record's real span begins
at 0-based **13,075,112**. They are the same coordinate. The padding offset *is*
the NCBI placement.

### The construction rule

A padded record is the parent chromosome's coordinate frame with the placed
interval **replaced** by the whole scaffold at the scaffold's own length:

```
padded_length = (parent_start - 1) + scaffold_length + (parent_chromosome_length - parent_stop)
                 └─ leading N ─┘     └─ real span ─┘    └────── trailing N ──────┘
```

For `HSCHR1_2_CTG3`: 13,075,112 + 256,271 + (248,956,422 − 13,312,803) =
248,975,002. The 18,580-base excess over chromosome 1 is exactly
256,271 − (13,312,803 − 13,075,113 + 1) — the scaffold is that much longer than
the interval it stands in for. Where a scaffold is *shorter* than its interval,
the padded record is shorter than the parent chromosome; chromosome 1's padded
records range from 248,926,662 to 248,975,002 around a chromosome length of
248,956,422.

This model — leading `N` count, real span length, and trailing `N` count all
predicted from the placement — holds for **441 of the 442** rows with a
GRCh38.p14 placement. The single exception is the placement revision discussed
below, and even there only the leading term is off; its trailing term is exact.

Two details are worth stating precisely, because "padding" invites a stronger
reading than the data supports.

**The padding is a placement, not a liftover.** The alignment between scaffold
and parent contains indels: NCBI's own columns say the aligned portion of the
scaffold starts at `alt_scaf_start` 20,633, leaving an `alt_start_tail` of
20,632 bases that do not align at all. Ensembl lays down the **whole** 256,271-
base scaffold beginning at `parent_start`, tail included. So position *i* in the
padded record equals position *i* on the parent only up to `parent_start`;
inside the scaffold the correspondence degrades across each indel, and after
`parent_stop` everything is shifted by `scaffold_length − interval_length` (here
+18,580). The representation is good for "roughly where does this scaffold sit";
it is not a coordinate transform.

**Minus-strand scaffolds are laid down reverse-complemented.** Of the 445, 365
have placement orientation `+`, 58 have `-`, 19 have `b` (mixed), and 3 have no
current placement row. For the 387 with `+`, `b`, or no orientation, stripping
the flanking `N` runs reproduces the true-length digest directly. For all 58
with `-`, the **reverse complement** of the stripped span reproduces it. Either
way the recovery is exact — see `recovers_unpadded` in the data file:
`forward` 387, `reverse_complement` 58, `no` **0**.

### The cross-check, run over all 445

Column `placement_matches` compares the measured 0-based `real_span_start`
against `parent_start − 1` for every row that has both:

| result | rows |
|---|---:|
| match | 441 |
| disagree | 1 |
| no placement row in GRCh38.p14 | 3 |

The three without a placement row — `CHR_HG107_PATCH`, `CHR_HG1311_PATCH`,
`CHR_HSCHR10_1_CTG4` — were retired from the assembly before p14, so the p14
placement files do not describe them. Their `parent_*` and `region_name` cells
are blank.

The single disagreement is `CHR_HSCHR5_8_CTG1` (novel-patch, region `GUSBP1`,
chromosome 5). Ensembl's offset is 0-based 21,481,418; GRCh38.p14 gives
`parent_start` 21,448,370. GRCh38.p12 gives `parent_start` **21,481,419** —
0-based 21,481,418, matching Ensembl exactly. NCBI extended the placement
leftward in p13; the padded record preserves the p12-era coordinate. The
disagreement is a placement revision, not a padding error.

Additionally, the length of the real span equals the true-length record's
length for **all 445** rows. 30 records have `real_bases < unpadded_length`
(10,208,740 bases in total): those are `N` runs *inside* the scaffold itself,
which the padded form cannot be blamed for.

## The switchover: 109 → 110

Per-release presence of one scaffold across the store's release-scoped
namespaces, `ensembl-75` … `ensembl-116`:

```
rel  75           absent
rel  76 .. 109    CHR_HSCHR1_2_CTG3   len=248975002   aDMwTJZkNZfyltbc6wRiMSb_DbTZEaRV
rel 110 .. 116    HSCHR1_2_CTG3       len=256271      5sJPeOQINUJr0syBJgsAX7bqCNmmWXJq
```

Never both in one release. The digest is identical in all 34 padded releases and
identical in all 7 unpadded ones. Release 75 is the last GRCh37 release
(`Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz`); this is a GRCh38 scaffold, so it
is absent there rather than unpadded.

Globally, one clean cut:

| | rel 109 | rel 110 |
|---|---:|---:|
| `CHR_`-prefixed records | **445** | **0** |
| bare alt/patch records | 169 | 681 |
| chromosomes | 25 | 25 |

The change lands in the release whose `dna.toplevel` was staged **2023-04-21**.
Ensembl FTP `Last-Modified` for each release's `dna.toplevel` (a staging date,
typically 2–3 months before the public release announcement):

| release | staged | release | staged |
|---:|---|---:|---|
| 108 | 2022-10-04 | 112 | 2024-02-13 |
| 109 | 2022-12-13 | 113 | 2024-08-15 |
| **110** | **2023-04-21** | 114 | 2025-01-30 |
| 111 | 2023-10-04 | 115 | 2025-07-07 |
| | | 116 | 2026-03-23 |

The full 42-release series is in
[`ensembl_dna_representation_by_release.tsv`](ensembl_dna_representation_by_release.tsv).
The padded count is not constant across the era — it grows as GRC patch
releases added scaffolds, and every one of the 445 persists to release 109:

| first appears in release | records added | roles |
|---:|---:|---|
| 76 | 261 | 261 alt-scaffold |
| 79 | 30 | 26 fix, 4 novel |
| 81 | 8 | 6 fix, 2 novel |
| 83 | 23 | 5 fix, 18 novel |
| 85 | 8 | 4 fix, 4 novel |
| 88 | 31 | 10 fix, 21 novel |
| 92 | 38 | 17 fix, 21 novel |
| 98 | 46 | 44 fix, 2 novel |

By role the 445 are 261 alt-scaffold, 112 fix-patch, 72 novel-patch. The 261
alt-scaffold count is exactly GRCh38's alt-scaffold set, present from the first
GRCh38 release. The patch counts are smaller than GRCh38.p14's 164 fix + 90
novel because the padded era ends at the patch level Ensembl 109 tracked.

## Why only alt loci and patches were padded

169 records keep bare names and true lengths right through the padded era. They
are exactly the **127 unplaced-scaffold plus 42 unlocalized-scaffold** rows of
the GRCh38 assembly report — verified by matching the GenBank accessions
Ensembl uses as names (`GL000008.2`, `KI270302.1`, …) against the
`GenBank-Accn` column of `GCF_000001405.40_GRCh38.p14_assembly_report.txt`: the
two sets are equal, not merely equinumerous.

That is the whole rule. **Padding was applied to sequences that have a
parent-chromosome placement, and only to those.** An unplaced or unlocalized
scaffold has no interval to be padded into, so it stays at true length in every
release from 76 to 116.

Release 110 does not merely drop the prefix — the bare count goes 169 → 681,
which is 169 + 445 + 67 newly added patch scaffolds, all at true length.

## What the padding costs to store

This store encodes DNA at 2 bits per base when the alphabet is exactly
`{A,C,G,T}`. `N` does not fit in 2 bits, so any record containing `N` falls back
to a 3-bit `dna3bit` encoding. Every padded record contains `N` by construction.

| alphabet | sequences | bases | share of bases |
|---|---:|---:|---:|
| `dna3bit` | **1,088** | **98,911,237,632** | **93.41%** |
| `dna2bit` | 1,426,421 | 4,972,744,406 | 4.70% |
| `dnaio` | 107 | 1,823,846,320 | 1.72% |
| `protein` | 351,881 | 184,072,924 | 0.17% |
| total | 1,779,497 | 105,891,901,282 | |

1,088 sequences out of 1.78 million hold 93% of the bases.

Not all of that is padding. Splitting `dna3bit` by what reaches it — a sequence
is *pre-110-only* if **every** alias it carries anywhere in the store lies in a
namespace `ensembl-75` … `ensembl-109`:

| bucket | sequences | bases | payload at 3 bits |
|---|---:|---:|---:|
| **padded alt/patch records** | **445** | **60,047,447,030** | **20.97 GiB** |
| release-75-only GRCh37 records | 215 | 29,192,103,640 | 10.20 GiB |
| superseded chromosome variants | 25 | 3,145,497,247 | 1.10 GiB |
| other pre-110-only | 4 | 5,527 | ~0 |
| reachable from 110+ or a non-release namespace | 399 | 6,526,184,188 | 2.28 GiB |
| total `dna3bit` | 1,088 | 98,911,237,632 | 34.54 GiB |

The padding is the single largest item but not the whole of it: release 75 is
GRCh37, and its chromosomes and scaffolds are not reproduced by any GRCh38
release or accession namespace, so they are pre-110-only for an unrelated
reason. The 25 superseded chromosome variants are chromosome records whose bytes
changed after release 109 — the release-110 PAR1 unmasking of chromosome `Y` is
one of them (see [`ENSEMBL_RELEASE_DRIFT.md`](ENSEMBL_RELEASE_DRIFT.md)).

The store as built is 48.9 GiB across 1,779,811 objects, of which roughly 37.6
GiB is encoded sequence payload (34.54 `dna3bit` + 1.16 `dna2bit` + 1.70
`dnaio` + 0.17 protein) and the remainder is indexes, aliases and metadata.
Dropping the 445 padded records alone would reclaim about **21 GiB**, leaving
roughly 28 GiB; dropping the entire pre-110-only `dna3bit` bucket would reclaim
about **32 GiB**, leaving roughly 17 GiB. For scale, the previously published
2026-07-22 store held 1,213,617 sequences in 9.2 GiB.

The per-release `dna3bit_bases` column shows the same thing from the other
side: 63,139,886,602 bases at release 109, 3,150,091,660 at release 110 — a 20x
drop with no loss of sequence content.

## Nothing exists only in padded form

Three checks, all exhaustive over the 445:

1. **Every padded record has a true-length counterpart already in the store.**
   443 resolve as `ensembl-110:<bare name>`; the two retired fix-patches
   `HG107_PATCH` and `HG1311_PATCH` resolve as `refseq:NW_015148966.1` and
   `refseq:NW_015148969.1`. Records reachable *only* in padded form: **0**.
2. **Every true-length counterpart's length agrees with NCBI.** 445/445 match
   the `Sequence-Length` column of the GRCh38 assembly report that last listed
   the scaffold.
3. **The true-length bytes are recoverable from the padded record.**
   `recovers_unpadded` is `forward` for 387 and `reverse_complement` for 58;
   `no` for none. (The reverse direction — rebuilding a padded record from the
   true-length one — needs the parent chromosome length and the placement, and
   is exact for the 441 rows where the construction rule holds.)

The padded records **carry no accession at all**. Across all 445, the number of
padded digests carrying *any* alias in a namespace other than `ensembl-76` …
`ensembl-109` is **0** — no `refseq`, no `insdc`, not even the rolling `ensembl`
namespace. They are an Ensembl FASTA presentation, not accessioned sequence
entities.

NCBI never published a padded form. `HSCHR1_2_CTG3` appears at 256,271 bases in
all 15 GRCh38 assembly reports that contain it, from the original GRCh38 through
p14 — one length, one GenBank accession `KI270766.1`, one RefSeq accession
`NT_187517.1`.

## Who actually references alt and patch sequences

### ClinVar

A complete scan of `ClinVarVCVRelease_2026-0621.xml.gz` — all 1,184,400,464
lines — matching every alt and patch accession drawn from the 23 local NCBI
assembly reports:

| sequence role | occurrences | distinct accessions used |
|---|---:|---|
| assembled-molecule (chromosomes) | 47,745,798 | 71 of 98 |
| alt-scaffold | 424,717 | 201 of 540 |
| fix-patch | 297,368 | 176 of 644 |
| novel-patch | 105,989 | 125 of 340 |
| unlocalized-scaffold | 366 | 5 of 120 |
| unplaced-scaffold | 35 | 4 of 318 |

Alt plus patch is 828,074 occurrences — **1.7048%** of chromosome-plus-alt
references — spread over 502 distinct accessions. Alt and patch sequences are a
small but non-trivial part of clinical variant data, and they are referenced by
accession.

Every one of those 502 accessions is at most 6,530,008 bases long; 452 are under
1 Mb and 50 are 1–10 Mb. All are true-length; none is padded. All 502 resolve in
this store today through the `refseq` and `insdc` namespaces, with **0**
missing.

### ClinGen Allele Registry

The registry validates the reference allele against real sequence, and accepts
alt scaffolds by accession:

| HGVS | result |
|---|---|
| `NC_000001.11:g.100A>T` | HTTP 400 `IncorrectReferenceAllele` — `"actualAllele": "N"` |
| `NT_187517.1:g.100A>T` | HTTP 200, normalized to `NT_187517.1:g.100A>T` |
| `KI270766.1:g.100A>T` | HTTP 200, normalized to `NT_187517.1:g.100A>T` |

Chromosome 1 position 100 is inside the leading telomeric `N` run, so `A` is
wrong there and the registry says so. The store confirms the alt scaffold's
1-based position 100 is `A`, consistent with the acceptance — the registry is
reading true-length `NT_187517.1`, and the GenBank synonym resolves to the same
sequence.

### seqrepo 2024-12-20

The snapshot carries **both** forms, with very different reachability. All 445
padded digests are present. For the worked example:

| digest | namespaces in seqrepo |
|---|---|
| padded `aDMwTJZk…` | `Ensembl:CHR_HSCHR1_2_CTG3` — plus `MD5`, `SHA1`, `SEGUID`, `VMC` |
| unpadded `5sJPeOQI…` | `NCBI:NT_187517.1`; `GRCh38` and `GRCh38.p1`–`p12` under both `HSCHR1_2_CTG3` and `chr1_KI270766v1_alt` — plus `MD5`, `SHA1`, `SEGUID`, `VMC` |

This holds across the whole set: for all 445 padded digests the only biological
namespace is `Ensembl` (column
`seqrepo_biological_namespaces_padded`), the rest being the four
digest-synthesis namespaces.

The distinction that matters: **anything resolving by accession reaches the
unpadded form.** The padded form is reachable only by asking Ensembl for a
`CHR_`-prefixed name. No accession, no assembly namespace, and no UCSC-style
name points at it.

## Consequence for the parity report

Dropping the 445 padded sequences from the store moves them from "in both" to
"seqrepo-only" in the digest-membership comparison:

| | with padding | without |
|---|---:|---:|
| seqrepo digests also in store | 980,664 | 980,219 |
| coverage of seqrepo's 1,144,093 digests | 85.715% | 85.677% |
| gap | 163,429 | 163,874 |

A 0.038 percentage-point cost, entirely in digests whose only seqrepo alias is
an Ensembl `CHR_` name.

## Data files

Both files are regenerated by `tools/ensembl_padding_probe.py`. It reads the
local store, the seqrepo `aliases.sqlite3` (read-only), the cached NCBI GRCh38
assembly reports under `downloads/`, and a directory of NCBI
`alt_scaffold_placement.txt` files (see [Reproducing this](#reproducing-this)
for the fetch).

```
uv run python tools/ensembl_padding_probe.py \
    --placements /tmp/placements \
    --out-dir seqrepo_equivalence \
    --measure all --last-modified
```

Runtime is about 8 minutes: ~4 minutes decoding all 445 padded records
(60 Gbase) plus the release scan and the FTP probe. `--measure none` or
`--measure N` caps the decoding; `--last-modified` is the only step that
requires network access to Ensembl.

### `ensembl_padded_scaffolds.tsv` — 445 rows, one per padded record

Sorted by `bare_name`. **All 445 rows are fully measured**; no column is
sampled or partially populated.

| column | meaning |
|---|---|
| `padded_name`, `bare_name` | `CHR_HSCHR1_2_CTG3` and `HSCHR1_2_CTG3` |
| `padded_digest`, `padded_length` | the padded record |
| `unpadded_digest`, `unpadded_length` | its true-length counterpart in this store |
| `inflation_factor` | `padded_length / unpadded_length` |
| `sequence_role` | alt-scaffold / fix-patch / novel-patch, from the assembly reports |
| `refseq_accn`, `genbank_accn` | from the assembly reports, matched on `Sequence-Name == bare_name` |
| `first_padded_release`, `last_padded_release` | scan of `ensembl-75` … `ensembl-116` |
| `first_unpadded_release`, `unpadded_source` | where the true-length form was found; blank release for the two resolved via `refseq` |
| `parent_chromosome`, `parent_start`, `parent_stop`, `region_name`, `orientation` | NCBI GRCh38.p14 placement; **blank for the 3 scaffolds with no p14 placement row** |
| `real_bases`, `real_span_start`, `real_span_end`, `pct_real` | measured by decoding the full padded record; span is 0-based half-open |
| `recovers_unpadded` | `forward` / `reverse_complement` / `no` — how the stripped span reproduces `unpadded_digest` |
| `placement_matches` | `real_span_start == parent_start − 1`; blank where no placement row |
| `in_seqrepo_padded`, `seqrepo_biological_namespaces_padded` | seqrepo `2024-12-20`, excluding `MD5`/`SEGUID`/`SHA1`/`VMC` |

Rows cited above: `CHR_HSCHR1_2_CTG3` (the worked example),
`CHR_HSCHR5_8_CTG1` (`placement_matches = no`), `CHR_HG107_PATCH` and
`CHR_HG1311_PATCH` (`unpadded_source = refseq:…`), `CHR_HG1817_1_PATCH`
(maximum inflation, 15,647.7x), `CHR_HG2365_PATCH` (minimum, 18.7x).

### `ensembl_dna_representation_by_release.tsv` — 42 rows, releases 75–116

Sorted by `release`. Columns: `release`, `assembly` (GRCh37 for 75, GRCh38
otherwise), `chr_prefixed_count`, `bare_alt_patch_count`, `chromosome_count`,
`dna3bit_bases` (summed over that release's DNA records), `ftp_last_modified`
(date portion of the `Last-Modified` header on that release's `dna.toplevel`).

## Reproducing this

Everything runs from the repo root. `store/` is the local RefgetStore; the
seqrepo snapshot is `~/dev/data/seqrepo/2024-12-20`.

Two conventions used throughout: a record is *padded* iff its Ensembl name
starts with `CHR_`, and the DNA records of an `ensembl-N` namespace are the
aliases that do **not** start with `ENS` (which is what `ENSG`/`ENST`/`ENSP`
cDNA, ncRNA and protein entries use).

### Measure real-vs-padded bases, the real span, and digest identity

```python
# uv run python
import base64, hashlib
from gtars.refget import RefgetStore

s = RefgetStore.open_local("store")

def sha512t24u(seq: str) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(seq.encode()).digest()[:24]).decode()

padded = s.get_sequence_metadata_by_alias("ensembl-100", "CHR_HSCHR1_2_CTG3")
bare   = s.get_sequence_metadata_by_alias("ensembl-116", "HSCHR1_2_CTG3")
print("padded", padded.sha512t24u, padded.length, str(padded.alphabet))
print("bare  ", bare.sha512t24u, bare.length, str(bare.alphabet))

seq   = s.get_substring(padded.sha512t24u, 0, padded.length)
real  = padded.length - seq.count("N") - seq.count("n")
start = padded.length - len(seq.lstrip("Nn"))
end   = len(seq.rstrip("Nn"))
print(f"real={real} ({100 * real / padded.length:.4f}%) span={start}..{end} spanlen={end - start}")
print("extracted digest", sha512t24u(seq[start:end]),
      "match:", sha512t24u(seq[start:end]) == bare.sha512t24u)

for ns, alias in [("ensembl-100", "1"), ("ensembl-100", "CHR_HG708_PATCH")]:
    m = s.get_sequence_metadata_by_alias(ns, alias)
    t = s.get_substring(m.sha512t24u, 0, m.length)
    r = m.length - t.count("N") - t.count("n")
    print(f"{ns}:{alias} total={m.length} real={r} ({100 * r / m.length:.2f}%)")
```

```
padded aDMwTJZkNZfyltbc6wRiMSb_DbTZEaRV 248975002 dna3bit
bare   5sJPeOQINUJr0syBJgsAX7bqCNmmWXJq 256271 dna2bit
real=256271 (0.1029%) span=13075112..13331383 spanlen=256271
extracted digest 5sJPeOQINUJr0syBJgsAX7bqCNmmWXJq match: True
ensembl-100:1 total=248956422 real=230481012 (92.58%)
ensembl-100:CHR_HG708_PATCH total=159270839 real=589656 (0.37%)
```

`s.get_sequence_metadata_by_alias` returns `None` for an absent alias rather
than raising, and `.alphabet` is an enum — wrap it in `str()` before formatting.

### Per-release presence trace

```python
# uv run python
from gtars.refget import RefgetStore
s = RefgetStore.open_local("store")
for rel in range(75, 117):
    hits = []
    for name in ("CHR_HSCHR1_2_CTG3", "HSCHR1_2_CTG3"):
        m = s.get_sequence_metadata_by_alias(f"ensembl-{rel}", name)
        if m is not None:
            hits.append(f"{name} len={m.length} {m.sha512t24u}")
    print(rel, "; ".join(hits) or "ABSENT")
```

Release 75 prints `ABSENT`; 76–109 print the padded record; 110–116 print the
bare one, with one digest per era.

### `CHR_` / bare / chromosome counts per release

```python
# uv run python
from gtars.refget import RefgetStore
s = RefgetStore.open_local("store")
CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
for rel in (75, 109, 110, 116):
    dna = [a for a in s.list_sequence_aliases(f"ensembl-{rel}") if not a.startswith("ENS")]
    chr_pref = [a for a in dna if a.startswith("CHR_")]
    chroms   = [a for a in dna if a in CHROMOSOMES]
    bare     = [a for a in dna if a not in CHROMOSOMES and not a.startswith("CHR_")]
    print(f"rel {rel}: CHR_={len(chr_pref)} bare={len(bare)} chromosomes={len(chroms)}")
```

```
rel 75: CHR_=0 bare=272 chromosomes=25
rel 109: CHR_=445 bare=169 chromosomes=25
rel 110: CHR_=0 bare=681 chromosomes=25
rel 116: CHR_=0 bare=681 chromosomes=25
```

### The 169 bare records are exactly unplaced + unlocalized

```python
# uv run python
from collections import Counter
from gtars.refget import RefgetStore

REPORT = ("downloads/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
          "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt")

s = RefgetStore.open_local("store")
CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
bare = {a for a in s.list_sequence_aliases("ensembl-109")
        if not a.startswith(("ENS", "CHR_")) and a not in CHROMOSOMES}

roles = {}   # Ensembl names these records by GenBank accession
for line in open(REPORT):
    if line.startswith("#") or not line.strip():
        continue
    f = line.split("\t")
    roles[f[4]] = f[1]

print(len(bare), Counter(roles.get(a, "NOT-IN-p14") for a in bare))
unplaced = {a for a, r in roles.items() if r in ("unplaced-scaffold", "unlocalized-scaffold")}
print("sets equal:", unplaced == bare)
```

```
169 Counter({'unplaced-scaffold': 127, 'unlocalized-scaffold': 42})
sets equal: True
```

### The construction rule, over all 445

Runs against the generated TSV, so it needs no sequence decoding.

```python
# uv run python
import csv
from gtars.refget import RefgetStore

s = RefgetStore.open_local("store")
chrlen = {c: s.get_sequence_metadata_by_alias("ensembl-100", c).length
          for c in [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]}

ok = bad = skipped = 0
for r in csv.DictReader(open("seqrepo_equivalence/ensembl_padded_scaffolds.tsv"), delimiter="\t"):
    if not r["parent_start"]:
        skipped += 1
        continue
    L = chrlen[r["parent_chromosome"]]
    ps, pe = int(r["parent_start"]), int(r["parent_stop"])
    predicted = (ps - 1) + int(r["unpadded_length"]) + (L - pe)
    lead  = int(r["real_span_start"]) == ps - 1
    trail = int(r["padded_length"]) - int(r["real_span_end"]) == L - pe
    if predicted == int(r["padded_length"]) and lead and trail:
        ok += 1
    else:
        bad += 1
        print("  exception:", r["padded_name"], int(r["padded_length"]), predicted, lead, trail)
print(f"model holds: {ok}  fails: {bad}  no placement row: {skipped}")
```

```
  exception: CHR_HSCHR5_8_CTG1 181668817 181635768 False True
model holds: 441  fails: 1  no placement row: 3
```

### Alphabet breakdown

```python
# uv run python  (~2 min; iterates all 1.78M sequences)
from collections import Counter
from gtars.refget import RefgetStore

s = RefgetStore.open_local("store")
n, b = Counter(), Counter()
for m in s.list_sequences():
    a = str(m.alphabet)
    n[a] += 1
    b[a] += m.length
total = sum(b.values())
for a in sorted(b, key=lambda x: -b[x]):
    print(f"{a:10s} {n[a]:>10,} {b[a]:>16,}  {100 * b[a] / total:5.2f}%")
print(f"{'total':10s} {sum(n.values()):>10,} {total:>16,}")
```

```
dna3bit         1,088   98,911,237,632  93.41%
dna2bit     1,426,421    4,972,744,406   4.70%
dnaio             107    1,823,846,320   1.72%
protein       351,881      184,072,924   0.17%
total       1,779,497  105,891,901,282
```

### `dna3bit` split by reachability

Reads the 445 padded digests from the generated TSV, then buckets every
`dna3bit` sequence by which namespaces reach it.

```python
# uv run python  (~3 min)
import csv, re
from gtars.refget import RefgetStore

s = RefgetStore.open_local("store")
padded = {r["padded_digest"] for r in
          csv.DictReader(open("seqrepo_equivalence/ensembl_padded_scaffolds.tsv"), delimiter="\t")}
CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

buckets = {}
def add(key, length):
    b = buckets.setdefault(key, [0, 0])
    b[0] += 1
    b[1] += length

for m in s.list_sequences():
    if str(m.alphabet) != "dna3bit":
        continue
    namespaces, aliases = set(), set()
    for a in s.get_aliases_for_sequence(m.sha512t24u):
        namespaces.add(a[0] if isinstance(a, (tuple, list)) else getattr(a, "namespace", str(a)))
        if isinstance(a, (tuple, list)) and len(a) > 1:
            aliases.add(a[1])
    releases = {int(ns.split("-")[1]) for ns in namespaces if re.fullmatch(r"ensembl-\d+", ns)}
    other    = {ns for ns in namespaces if not re.fullmatch(r"ensembl-\d+", ns)}
    if other or not releases or max(releases) >= 110:
        add("reachable from 110+ or a non-release namespace", m.length)
    elif m.sha512t24u in padded:
        add("padded alt/patch record", m.length)
    elif releases == {75}:
        add("release-75-only GRCh37 record", m.length)
    elif aliases & CHROMOSOMES:
        add("superseded chromosome variant", m.length)
    else:
        add("other pre-110-only", m.length)

for k in sorted(buckets, key=lambda k: -buckets[k][1]):
    n, b = buckets[k]
    print(f"{k:48s} {n:>5} seqs {b:>16,} bases {b * 3 / 8 / 2**30:>7.2f} GiB")
```

```
padded alt/patch record                            445 seqs   60,047,447,030 bases   20.97 GiB
release-75-only GRCh37 record                      215 seqs   29,192,103,640 bases   10.20 GiB
reachable from 110+ or a non-release namespace     399 seqs    6,526,184,188 bases    2.28 GiB
superseded chromosome variant                       25 seqs    3,145,497,247 bases    1.10 GiB
other pre-110-only                                   4 seqs            5,527 bases    0.00 GiB
```

The 445 landing in `padded alt/patch record` rather than the first bucket is
also the proof that no padded digest carries an alias outside `ensembl-76` …
`ensembl-109`: any `refseq`, `insdc` or rolling `ensembl` alias would have put
it in `other` and sorted it into the first bucket.

### seqrepo aliases for both digests

```sh
SR=~/dev/data/seqrepo/2024-12-20/aliases.sqlite3
for d in aDMwTJZkNZfyltbc6wRiMSb_DbTZEaRV 5sJPeOQINUJr0syBJgsAX7bqCNmmWXJq; do
  echo "== $d"
  sqlite3 "file:$SR?mode=ro" \
    "select namespace, alias from seqalias where seq_id='$d' order by namespace, alias;"
done
```

The padded digest returns five rows — `Ensembl|CHR_HSCHR1_2_CTG3` plus `MD5`,
`SHA1`, `SEGUID`, `VMC`. The unpadded digest returns `NCBI|NT_187517.1`, the
`GRCh38`/`GRCh38.p1`–`p12` entries under both `HSCHR1_2_CTG3` and
`chr1_KI270766v1_alt`, and the same four digest-synthesis rows.

### NCBI placement files

GRCh38.p14 splits alt loci across 35 `ALT_REF_LOCI_<N>` units, with patches in a
separate `PATCHES` unit. 514 placement rows in total.

```sh
DEST=/tmp/placements; mkdir -p "$DEST"
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure
for i in $(seq 1 40); do
  curl -sfL --retry 2 "$BASE/ALT_REF_LOCI_$i/alt_scaffolds/alt_scaffold_placement.txt" \
       -o "$DEST/ALT_REF_LOCI_$i.txt" || rm -f "$DEST/ALT_REF_LOCI_$i.txt"
done
curl -sfL --retry 2 "$BASE/PATCHES/alt_scaffolds/alt_scaffold_placement.txt" -o "$DEST/PATCHES.txt"

awk -F'\t' '$3 == "HSCHR1_2_CTG3" { print FILENAME; print }' "$DEST"/*.txt
```

The loop runs to 40 and discards misses so it survives a future unit count; 35
succeed today. To reproduce the `CHR_HSCHR5_8_CTG1` placement revision, swap
`GCF_000001405.40_GRCh38.p14` for `GCF_000001405.38_GRCh38.p12` and read the
`PATCHES` file.

### Ensembl FTP `Last-Modified`

```sh
for rel in $(seq 75 116); do
  dir="https://ftp.ensembl.org/pub/release-$rel/fasta/homo_sapiens/dna/"
  f=$(curl -s "$dir" | grep -o 'Homo_sapiens[^"]*dna\.toplevel\.fa\.gz' | head -1)
  printf '%s\t%s\t%s\n' "$rel" "$f" \
    "$(curl -sI "$dir$f" | grep -i '^last-modified' | tr -d '\r')"
done
```

```
109	Homo_sapiens.GRCh38.dna.toplevel.fa.gz	Last-Modified: Tue, 13 Dec 2022 00:02:32 GMT
110	Homo_sapiens.GRCh38.dna.toplevel.fa.gz	Last-Modified: Fri, 21 Apr 2023 16:16:00 GMT
```

Release 75's filename is `Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz`, which the
same pattern matches.

### ClinGen Allele Registry

```sh
for h in "NC_000001.11:g.100A>T" "NT_187517.1:g.100A>T" "KI270766.1:g.100A>T"; do
  printf '%s -> HTTP %s\n' "$h" \
    "$(curl -s -o /dev/null -w '%{http_code}' --get \
        --data-urlencode "hgvs=$h" http://reg.genome.network/allele)"
done
```

## What excluding the pre-110 padding would and would not lose

**Would not lose any sequence content.** All 445 padded records reduce to
true-length sequences the store already holds under `ensembl-110`+ or
`refseq`, exactly, by stripping flanking `N` (387) or by stripping and
reverse-complementing (58). No sequence exists only in padded form.

**Would not lose any accession-addressable lookup.** The padded digests carry no
alias outside `ensembl-76`–`ensembl-109`. Nothing that resolves by RefSeq,
GenBank, UCSC-style, or assembly-scoped name touches them, and the downstream
consumers examined here — ClinVar and the ClinGen Allele Registry — reference
alt and patch sequences exclusively by accession at true length.

**Would lose byte-level fidelity to `Homo_sapiens.GRCh38.dna.toplevel.fa.gz` for
releases 76–109.** A consumer that asks `ensembl-100` for `CHR_HSCHR1_2_CTG3`
would get nothing, and a consumer that reconstructs those FASTA files from the
store would not reproduce them. That is a real property of the current build:
each release-scoped namespace is a faithful record of what Ensembl published.

**Would lose one seqrepo correspondence.** 445 digests that seqrepo holds under
`Ensembl:CHR_*` would move from "in both" to "seqrepo-only", taking digest
coverage from 85.715% to 85.677%.

**Would reclaim about 21 GiB** — roughly 43% of the store's 48.9 GiB, and 56.7%
of every base it holds.

The trade is byte-level fidelity to 34 Ensembl releases' `dna.toplevel`, plus
0.038 points of seqrepo digest coverage, against roughly 43% of the store's
size.

## Outcome

The padded records are excluded at ingest. Each affected seqset declares the rule
in `sources.toml`:

```toml
exclude = { file_classes = ["dna.toplevel"], record_prefixes = ["CHR_"] }
```

Releases 76–109 carry it; release 75 (GRCh37) and 110–116 do not, because they
contain no padded records. Matching is on the record name alone — the `CHR_`
prefix is Ensembl's own marker for this representation and is exact across those
releases, so no sequence-content inspection is involved.

Mechanically, filtration is a preflight pass alongside downloads: every declared
source is rewritten once, concurrently (`--filter-jobs`), to a cached
`<artifact>.filtered.fa.gz` beside the download, and that file is what gtars
ingests. A `<artifact>.excluded.tsv` records every dropped record, and a rule
matching zero records is an error rather than a silent no-op.

Verified on release 100: 445 records dropped, 0 warnings, and the surviving
aliases are the production set minus exactly the 445 `CHR_` entries with none
added. Chromosome digests are unchanged —
`ensembl-100:1` is `2YnepKM7OkBoOrKmvHbGqguVfF9amCST` before and after.

The build lock does not pin this outcome; it compares upstream `sha256` only.
Tracked in `ISSUE-lock-should-pin-ingest-outcome.md`.
