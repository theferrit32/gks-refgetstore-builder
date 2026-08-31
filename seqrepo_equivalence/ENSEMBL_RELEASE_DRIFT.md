# Silent sequence drift across Ensembl releases

**Question.** The 2026-08-28 parity run against seqrepo `2024-12-20` reported 317
shared aliases whose digests disagree. 115 were the known `star_normalization`
cases. The other 202 were dismissed as "release drift, not data loss." This
report tests that dismissal.

**Verdict.** The dismissal was directionally right and materially incomplete.
The store is faithful to upstream — that part holds, and is now proven rather
than asserted. But the parity check sees only a small, biased slice of the
drift, and the 202 turn out to be the *tail* of a much larger upstream event.

**Silent drift** here means: the same Ensembl `accession.version` string carries
different sequence bytes in different releases. Ensembl's versioning contract
implies this should not happen; a version suffix exists precisely so that a
sequence change is announced by a version bump.

## Scope: every consecutive release, 75 → 116

41 transitions, comparing digests for aliases present on both sides.

| era | transitions | silent-change events | character |
|---|---:|---:|---|
| 75 → 85 | 10 | 39,738 | frequent ±1 terminal-residue churn |
| 85 → 116 | 31 | 4,413 | frozen, except three discrete events |

**28 of 41 transitions changed nothing at all.** The modern era is far more
stable than the raw total suggests — 4,413 of its 4,413 events fall in just five
transitions, and 4,407 of those in two adjacent ones:

| transition | events | breakdown |
|---|---:|---|
| 100 → 101 | 1 | one length change |
| 109 → 110 | 1 | chromosome `Y` (see below) |
| 112 → 113 | 4 | length changes |
| **113 → 114** | **2,329** | `pos0 X→M` 2,074; `pos0 L→M` 249; length 6 |
| **114 → 115** | **2,078** | `pos0 M→X` 2,073; length 5 |

## The 113 → 114 → 115 event

Release 114 forced the **initiator residue to Methionine** across the board,
then release 115 partially reverted it. Every one of the 2,323 position-0 cases
is a single-character change at index 0 with no length change.

| class | n | at rel 114 | at rel 115–116 | interpretation |
|---|---:|---|---|---|
| ambiguous start | 2,074 | `X` → `M` | **reverted to `X`** (2,073) | start codon is unresolvable (N in the genomic sequence); asserting `M` claims knowledge that isn't there |
| non-ATG start | 249 | `L` → `M` | **persists as `M`** | initiator `CTG` is charged with Met — the `L` was the bug, and the fix stuck |
| trailing stop | 5 | `*` stripped | persists | e.g. `VLRYFDWLL*` → `VLRYFDWLL` |
| replacement | 1 | 321 aa → 502 aa | persists | `ENSP00000520925.1`, no shared 30-mer with its predecessor |

Observed first-residue trajectories across releases 112–116, derived
independently (see *Reproduction*):

```
X X M X X   x2068     ambiguous  -- one-release blip, reverted
L L M M M   x249      non-ATG    -- persistent correction
V V V V V   x4        trailing '*' cases (position 0 never changed)
. I M M M   x1        ENSP00000520925.1 whole-sequence replacement
```

Release 115 keeping the `L→M` fix while backing out `X→M` is coherent: a `CTG`
initiator *is* translated as Met, but a start codon containing an `N` genuinely
has an unknown product. Release 114 over-applied the normalization; 115 narrowed
it to the cases where the codon is actually known.

### Why parity saw only 200 of 4,407

The parity check compares seqrepo against the **rolling `ensembl` namespace**,
which points at the newest release (116). That makes it structurally blind to
any change that is reverted before the newest release:

- **`X→M` (2,074):** reverted at 115, so release 116 again agrees with seqrepo.
  Zero mismatch rows. Completely invisible.
- **`L→M` (249):** persists, so 116 disagrees with seqrepo. 200 surfaced — the
  other 49 are simply absent from seqrepo.

So the mismatch count is not a measure of drift. It is a measure of *drift that
survived to the newest release and happens to overlap seqrepo's accession set*.
Both filters are severe. **The parity report understated this event by ~22x.**

## Chromosome `Y`: release 110 unmasked PAR1

The one non-protein case. Same length in both releases, differing at 2,778,688
positions:

| release | length | `N` bases | non-`N` |
|---|---:|---:|---:|
| ensembl-109 | 57,227,415 | 33,591,060 | 23,636,355 |
| ensembl-110 | 57,227,415 | 30,812,372 | 26,415,043 |

**2,778,688 bases unmasked**, against a GRCh38 PAR1 of 2,781,479 bp. Release 110
stopped hard-masking the Y pseudoautosomal region.

seqrepo's `Ensembl:Y` digest (`yvJNy7T8_u6k…`) matches our **ensembl-79 through
ensembl-109**. Its protein entries correspond to ~113. This is worth stating
plainly: **seqrepo's `Ensembl` namespace is not a single-release snapshot.** It
accretes over time, so "which Ensembl release is seqrepo" has no single answer,
and any parity narrative that assumes one is unreliable.

## The early era (75 → 85)

Dominated by ±1 terminal-residue churn, and it oscillates rather than
progresses. `ENSP00000448466.1` is representative — verified against the raw
upstream `pep.all.fa.gz` for each release:

```
rel83  len=231  ...IKAYPRLGPPTPGE
rel84  len=232  ...KAYPRLGPPTPGEP     <- one residue added
rel85  len=231  ...IKAYPRLGPPTPGE     <- and removed again
```

Release 84 is a one-release anomaly in the same shape as release 114 — an
Ensembl-wide policy applied for exactly one release and then withdrawn. The
mechanism (terminal partial-codon handling) is characterized here but not
root-caused; it predates the store's practical range of interest.

Two transitions have near-disjoint alias sets (82→83 shares only 493 aliases;
95→96 similar). Those are mass version bumps — the honest kind of change, and
invisible to this analysis by construction, since a bumped version is a
different alias.

## Reproduction: the store is faithful

To rule out an ingest artifact, `tools/ensembl_drift_probe.py` rebuilds a
minimal store containing *only* the 2,329 affected accessions, read straight
from the cached upstream FASTAs, touching neither the production store nor the
network:

```
uv run python tools/ensembl_drift_probe.py \
    --releases 112 113 114 115 116 \
    --accessions <(…) --out /tmp/probe-store
```

Cross-checking all 11,645 (accession, release) pairs against production:

```
digest identical: 11,636      both absent: 9      MISMATCH: 0
```

Every trajectory above was then re-derived from the probe store alone and
matched. The drift is upstream in Ensembl's published files; the store records
it correctly. Spot-checked directly in the gzipped source FASTAs as well:

```
rel113  >ENSP00000219542.3 pep chromosome:GRCh38:16:715338:717390:1 …
        LAPAARAGYSEERCSWRGSGLTQEPGSVGQLALACAEGAV…
rel114  >ENSP00000219542.3 pep chromosome:GRCh38:16:715338:717390:1 …
        MAPAARAGYSEERCSWRGSGLTQEPGSVGQLALACAEGAV…
```

Identical header, identical coordinates, identical version — different sequence.

## Consequences

1. **A versioned Ensembl accession is not a stable identifier.** 4,407 accessions
   changed sequence without a version bump inside two releases. Any consumer
   treating `ENSP…​.N` as immutable — VRS normalization caches, HGVS validators,
   anything keyed on accession rather than digest — can silently disagree with
   itself across an Ensembl upgrade.
2. **The release-scoped namespaces earn their cost.** Because the store pins
   every release, all four variants of these sequences remain addressable and
   the whole event is reconstructable after the fact. Against a store that only
   tracked "current," this would have been undetectable.
3. **The rolling `ensembl` namespace is a convenience, not a citation.** It is
   the right default for "what does Ensembl say now" and the wrong thing to
   record in anything durable. Cite `ensembl-116`.
4. **Digest mismatch counts against seqrepo are a weak drift metric** — they miss
   reverted changes entirely and are clipped to seqrepo's accession set. Measure
   drift release-to-release inside the store instead.

## Standing caveats

- Nothing here is a defect in this repository. No corrective action is implied
  for the store; the 2026-08-26 build is correct as published.
- The `L→M` reading rests on the substitution pattern and Ensembl's revert
  behavior, not on inspecting the underlying codons — the store holds no CDS
  FASTA, so the `CTG` initiator claim is inference, well-supported but not
  directly verified.
- The sweep covers every alias in each release namespace, so `cdna`/`ncrna`
  (`ENST`) and `dna` were included alongside `pep` (`ENSP`). Transcripts are far
  quieter than proteins but not silent — the modern era's changes break down as
  2,329 + 2,078 + 1 `ENSP`, 4 `ENST` (all at 112→113, all length changes), and 1
  `dna` (chromosome `Y`). The four `ENST` cases were not characterized further.
