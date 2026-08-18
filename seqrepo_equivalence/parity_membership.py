#!/usr/bin/env python
"""Exhaustive per-sequence parity between our RefgetStore and a seqrepo snapshot.

Where ``verify_seqrepo_equivalence.py`` produces a *categorized* backwards-compat
verdict, this script produces the *exhaustive* view: one row per sequence digest
across the union of both stores, recording whether each sequence is present in
seqrepo, present in our RefgetStore, or both, plus the identifiers each side
knows it by. It also emits a human-readable summary report.

Two outputs:
  * ``parity_by_digest.tsv`` — one row per sha512t24u digest (the membership table).
  * ``parity-summary.md``    — the written report (totals, per-prefix rollups,
    what source groups we can load, the residual gap).

gtars + stdlib only. seqrepo is read from ``aliases.sqlite3`` via stdlib
``sqlite3``; the store's aliases are read from its on-disk ``aliases/sequences``
TSV sidecars. Reuses the proven readers/helpers in
``verify_seqrepo_equivalence.py`` rather than duplicating them.

Examples:
    uv run python seqrepo_equivalence/parity_membership.py
    uv run python seqrepo_equivalence/parity_membership.py --store ./store \
        --out-tsv ./runs/DATE-parity/parity_by_digest.tsv \
        --report ./runs/DATE-parity/parity-summary.md
"""

from __future__ import annotations

import argparse
import datetime as _dt
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path

from gtars.refget import RefgetStore

from verify_seqrepo_equivalence import (
    DEFAULT_SEQREPO,
    EXPECTED_OMIT,
    NAMESPACE_MAP,
    accession_prefix,
    build_refget_digest_set,
    load_refget_aliases,
    open_seqrepo_db,
)

logger = logging.getLogger("parity_membership")

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DEFAULT_STORE = REPO_ROOT / "store"
DEFAULT_RUN_DIR = REPO_ROOT / "runs" / "2026-07-02-seqrepo-parity"
DEFAULT_OUT_TSV = DEFAULT_RUN_DIR / "parity_by_digest.tsv"
DEFAULT_REPORT = DEFAULT_RUN_DIR / "parity-summary.md"
DEFAULT_KNOWN_DIVERGENT = HERE / "ensembl_known_divergent.txt"


# --------------------------------------------------------------------------- seqrepo


def load_seqrepo_membership(
    conn, alias_cap: int
) -> tuple[set[str], dict[str, set[str]], dict[str, list[str]], dict[str, int]]:
    """Single scan of seqrepo's biological aliases.

    Returns:
      digests          - every distinct current seq_id (includes digest-only seqs)
      ns_by_digest     - digest -> set of biological namespaces referencing it
      aliases_by_digest- digest -> up to ``alias_cap`` example accessions
      alias_count      - digest -> total biological alias count

    The four synthesized digest namespaces (MD5/SEGUID/SHA1/VMC) are excluded from
    the alias/namespace detail: they are 1:1 re-encodings of the sequence, present
    for every seqrepo digest, so they carry no membership information. They ARE
    still counted in ``digests`` via the DISTINCT seq_id scan below.
    """
    cur = conn.execute("SELECT DISTINCT seq_id FROM seqalias WHERE is_current=1")
    digests: set[str] = {row[0] for row in cur}

    placeholders = ",".join("?" for _ in EXPECTED_OMIT)
    ns_by_digest: dict[str, set[str]] = defaultdict(set)
    aliases_by_digest: dict[str, list[str]] = defaultdict(list)
    alias_count: dict[str, int] = defaultdict(int)
    cur = conn.execute(
        f"SELECT seq_id, namespace, alias FROM seqalias "
        f"WHERE is_current=1 AND namespace NOT IN ({placeholders})",
        tuple(EXPECTED_OMIT),
    )
    for seq_id, namespace, alias in cur:
        ns_by_digest[seq_id].add(namespace)
        alias_count[seq_id] += 1
        bucket = aliases_by_digest[seq_id]
        if len(bucket) < alias_cap:
            bucket.append(alias)
    return digests, ns_by_digest, aliases_by_digest, alias_count


# --------------------------------------------------------------------------- refget


def load_refget_membership(
    store: RefgetStore, store_path: Path, alias_cap: int
) -> tuple[set[str], dict[str, set[str]], dict[str, list[str]], dict[str, int]]:
    """Digest set + digest->aliases/namespaces from the store's TSV sidecars."""
    digests = build_refget_digest_set(store)
    ns_by_digest: dict[str, set[str]] = defaultdict(set)
    aliases_by_digest: dict[str, list[str]] = defaultdict(list)
    alias_count: dict[str, int] = defaultdict(int)
    for ns in store.list_sequence_alias_namespaces():
        for alias, digest in load_refget_aliases(store_path, ns).items():
            ns_by_digest[digest].add(ns)
            alias_count[digest] += 1
            bucket = aliases_by_digest[digest]
            if len(bucket) < alias_cap:
                bucket.append(alias)
    return digests, ns_by_digest, aliases_by_digest, alias_count


# --------------------------------------------------------------------------- helpers


def _fmt_aliases(sample: list[str], total: int, cap: int) -> str:
    if not sample:
        return ""
    joined = ";".join(sample)
    if total > len(sample):
        joined += f";+{total - len(sample)}"
    return joined


def _pick_prefix(
    digest: str,
    sr_aliases: dict[str, list[str]],
    rg_aliases: dict[str, list[str]],
) -> str:
    """Coarse accession category from a representative alias (seqrepo first)."""
    for src in (sr_aliases, rg_aliases):
        aliases = src.get(digest)
        if aliases:
            return accession_prefix(aliases[0])
    return "(none)"


def refseq_prefix_breakdown(store_path: Path, namespace: str) -> Counter:
    """Per-prefix alias counts within a single store namespace (e.g. refseq)."""
    counts: Counter = Counter()
    for alias in load_refget_aliases(store_path, namespace):
        counts[accession_prefix(alias)] += 1
    return counts


def count_known_divergent(path: Path) -> int:
    if not path.exists():
        return 0
    n = 0
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                n += 1
    return n


# --------------------------------------------------------------------------- report


def render_report(
    *,
    seqrepo_root: Path,
    store_path: Path,
    sr_digests: set[str],
    rg_digests: set[str],
    both: int,
    sr_only: int,
    rg_only: int,
    sr_only_by_prefix: Counter,
    sr_only_sourceable: int,
    sr_only_noncurrent_assembly: int,
    sr_only_digest_only: int,
    rg_only_by_prefix: Counter,
    both_by_prefix: Counter,
    refget_ns_alias_counts: dict[str, int],
    refseq_breakdown: Counter,
    ensembl_breakdown: Counter,
    known_divergent: int,
) -> str:
    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    union = both + sr_only + rg_only
    coverage = 100.0 * both / max(1, len(sr_digests))

    def table(counter: Counter, header: str) -> str:
        rows = sorted(counter.items(), key=lambda kv: -kv[1])
        out = [f"| {header} | count |", "|---|---|"]
        for k, v in rows:
            out.append(f"| `{k}` | {v:,} |")
        return "\n".join(out)

    ns_rows = sorted(refget_ns_alias_counts.items(), key=lambda kv: -kv[1])
    ns_table = "\n".join(
        [f"| `{ns}` | {n:,} |" for ns, n in ns_rows]
    )

    lines = f"""# Sequence parity: our RefgetStore vs seqrepo {seqrepo_root.name}

_Generated {now} by `parity_membership.py`. Numbers are computed, not hand-typed._

Companion to the run record README (source strategy) and
`seqrepo_equivalence/verify_seqrepo_equivalence.py` (categorized backwards-compat
verdict). This report is the **exhaustive per-sequence view**; the full membership
table is `parity_by_digest.tsv` (one row per sequence digest).

## Headline

| metric | value |
|---|---|
| seqrepo distinct sequence digests | {len(sr_digests):,} |
| our store distinct sequence digests | {len(rg_digests):,} |
| union (rows in the TSV) | {union:,} |
| **in both** | **{both:,}** |
| seqrepo-only (the gap) | {sr_only:,} |
| our-store-only (extensions) | {rg_only:,} |
| **seqrepo digest coverage** | **{coverage:.3f}%** |

"Coverage" = fraction of seqrepo's sequences that are present (by digest) in our
store. A digest counts as covered regardless of which alias spelling either side
uses, because the digest IS the sequence.

## What sequence groups we can load

These are the source groups declared in `sources.toml`, all ingested into the
store this build. Alias counts are what actually landed (per store namespace):

| store namespace | aliases |
|---|---|
{ns_table}

RefSeq (`refseq`) alias breakdown by accession prefix:

{table(refseq_breakdown, "refseq prefix")}

Ensembl (`ensembl`) alias breakdown by accession prefix:

{table(ensembl_breakdown, "ensembl prefix")}

Source groups behind these namespaces:
- **RefSeq current transcripts/proteins** — `human.{{1..15}}.rna` / `.protein`
  shards (NM/NR/XM/XR, NP/XP).
- **RefSeqGene** — `refseqgene.{{1..9}}.genomic` (NG_).
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

Sequences seqrepo has that our store does not (by digest): **{sr_only:,}** total.
Not all of these are *sourceable* — many seqrepo entries have no biological
accession we could load from FTP. Split by what seqrepo knows them as:

| bucket | count | can we close it? |
|---|---|---|
| **sourceable** (referenced by a current NCBI/Ensembl/GRCh3x alias) | **{sr_only_sourceable:,}** | yes — needs eutils backfill or a not-yet-loaded release |
| non-current-assembly-only (only old patch/assembly namespaces: NCBI3x, GRCh37 patches, hs37d5, JRGv1/2, CHM1) | {sr_only_noncurrent_assembly:,} | out of scope — assemblies we deliberately don't load |
| **digest-only** (no biological accession *anywhere* in seqrepo — only MD5/SEGUID/SHA1/VMC) | **{sr_only_digest_only:,}** | not from FTP — seqrepo itself records no accession; only obtainable by copying the raw bytes out of seqrepo |

So the **actionable** transcript/protein gap is ~{sr_only_sourceable:,}, not the
headline {sr_only:,}. The dominant `(none)` bucket below is the digest-only
sequences — opaque entries seqrepo holds without any accession.

Full breakdown by accession prefix (`(none)` = the digest-only bucket):

{table(sr_only_by_prefix, "prefix")}

Expected sourceable residual (see the run record README): the predicted-model tail
(`XR_`/`XM_`, renumbered per annotation release), a handful of fully-suppressed
`NM_`/`NR_`, predicted proteins (`XP_`/`NP_`), and low-count genomic contigs
(`NT_`/`NW_`). Recoverable only via eutils or accepted as documented drift; a
fresh seqrepo build from current FTP would hit the same wall.

## Our extensions — store-only sequences by prefix

Sequences we hold that seqrepo 2024-12-20 does not (newer releases, `insdc`,
GRCh38.p14, the Ensembl superset).

{table(rg_only_by_prefix, "prefix")}

## In both — shared sequences by prefix

{table(both_by_prefix, "prefix")}

## Known-divergent note

`ensembl_known_divergent.txt` lists **{known_divergent}** Ensembl accessions whose
digests differ from seqrepo by design (mostly `*`-stop-codon normalization). These
are correct-by-spec differences, not gaps: the sequences are present, only the
digest differs, so they surface as distinct digests on each side rather than as
`both`.
"""
    return lines


# --------------------------------------------------------------------------- main


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--seqrepo", type=Path, default=DEFAULT_SEQREPO)
    ap.add_argument("--store", type=Path, default=DEFAULT_STORE)
    ap.add_argument("--out-tsv", type=Path, default=DEFAULT_OUT_TSV)
    ap.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    ap.add_argument("--known-divergent", type=Path, default=DEFAULT_KNOWN_DIVERGENT)
    ap.add_argument("--alias-cap", type=int, default=8,
                    help="max example aliases listed per side per digest")
    ap.add_argument("--verbose", "-v", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if not args.store.exists():
        logger.error("store not found at %s", args.store)
        return 2

    logger.info("reading seqrepo aliases from %s", args.seqrepo)
    conn = open_seqrepo_db(args.seqrepo)
    sr_digests, sr_ns, sr_aliases, sr_alias_count = load_seqrepo_membership(
        conn, args.alias_cap
    )
    logger.info("seqrepo: %d distinct digests", len(sr_digests))

    logger.info("opening store at %s", args.store)
    store = RefgetStore.open_local(str(args.store))
    store.pull_aliases()
    rg_digests, rg_ns, rg_aliases, rg_alias_count = load_refget_membership(
        store, args.store, args.alias_cap
    )
    logger.info("store: %d distinct digests", len(rg_digests))

    union = sr_digests | rg_digests
    logger.info("union: %d digests; writing %s", len(union), args.out_tsv)

    both = sr_only = rg_only = 0
    sr_only_by_prefix: Counter = Counter()
    rg_only_by_prefix: Counter = Counter()
    both_by_prefix: Counter = Counter()
    # sourceability split of the seqrepo-only gap
    mapped_ns = set(NAMESPACE_MAP)
    sr_only_sourceable = sr_only_noncurrent_assembly = sr_only_digest_only = 0

    with args.out_tsv.open("w", encoding="utf-8") as out:
        out.write(
            "digest\tmembership\tin_seqrepo\tin_refget\tseqrepo_ns\trefget_ns\t"
            "seqrepo_bio_aliases\trefget_aliases\tprefix\n"
        )
        for digest in sorted(union):
            in_sr = digest in sr_digests
            in_rg = digest in rg_digests
            if in_sr and in_rg:
                membership = "both"
                both += 1
            elif in_sr:
                membership = "seqrepo_only"
                sr_only += 1
            else:
                membership = "refget_only"
                rg_only += 1

            prefix = _pick_prefix(digest, sr_aliases, rg_aliases)
            if membership == "both":
                both_by_prefix[prefix] += 1
            elif membership == "seqrepo_only":
                sr_only_by_prefix[prefix] += 1
                ns_set = sr_ns.get(digest, ())
                if not ns_set:
                    sr_only_digest_only += 1
                elif set(ns_set) & mapped_ns:
                    sr_only_sourceable += 1
                else:
                    sr_only_noncurrent_assembly += 1
            else:
                rg_only_by_prefix[prefix] += 1

            seqrepo_ns = ";".join(sorted(sr_ns.get(digest, ()))) if in_sr else ""
            refget_ns = ";".join(sorted(rg_ns.get(digest, ()))) if in_rg else ""
            sr_al = _fmt_aliases(
                sr_aliases.get(digest, []), sr_alias_count.get(digest, 0), args.alias_cap
            )
            rg_al = _fmt_aliases(
                rg_aliases.get(digest, []), rg_alias_count.get(digest, 0), args.alias_cap
            )
            out.write(
                f"{digest}\t{membership}\t{int(in_sr)}\t{int(in_rg)}\t"
                f"{seqrepo_ns}\t{refget_ns}\t{sr_al}\t{rg_al}\t{prefix}\n"
            )

    logger.info(
        "both=%d seqrepo_only=%d refget_only=%d", both, sr_only, rg_only
    )

    refget_ns_alias_counts = {
        ns: len(load_refget_aliases(args.store, ns))
        for ns in store.list_sequence_alias_namespaces()
    }
    refseq_breakdown = refseq_prefix_breakdown(args.store, "refseq")
    ensembl_breakdown = refseq_prefix_breakdown(args.store, "ensembl")
    known_divergent = count_known_divergent(args.known_divergent)

    report = render_report(
        seqrepo_root=args.seqrepo,
        store_path=args.store,
        sr_digests=sr_digests,
        rg_digests=rg_digests,
        both=both,
        sr_only=sr_only,
        rg_only=rg_only,
        sr_only_by_prefix=sr_only_by_prefix,
        sr_only_sourceable=sr_only_sourceable,
        sr_only_noncurrent_assembly=sr_only_noncurrent_assembly,
        sr_only_digest_only=sr_only_digest_only,
        rg_only_by_prefix=rg_only_by_prefix,
        both_by_prefix=both_by_prefix,
        refget_ns_alias_counts=refget_ns_alias_counts,
        refseq_breakdown=refseq_breakdown,
        ensembl_breakdown=ensembl_breakdown,
        known_divergent=known_divergent,
    )
    args.report.write_text(report, encoding="utf-8")
    logger.info("wrote report %s", args.report)

    # console headline
    coverage = 100.0 * both / max(1, len(sr_digests))
    print()
    print(f"seqrepo digests : {len(sr_digests):,}")
    print(f"store digests   : {len(rg_digests):,}")
    print(f"in both         : {both:,}  ({coverage:.3f}% of seqrepo)")
    print(f"seqrepo-only    : {sr_only:,}")
    print(f"store-only      : {rg_only:,}")
    print(f"TSV   : {args.out_tsv}")
    print(f"report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
