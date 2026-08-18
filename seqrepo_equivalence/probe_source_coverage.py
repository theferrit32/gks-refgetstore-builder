#!/usr/bin/env python
"""Measure how much of the seqrepo gap each candidate source FASTA would cover.

Given a list of candidate source FASTAs (NCBI annotation-release / per-patch /
Ensembl-release URLs) and the gap feed from ``verify_seqrepo_equivalence.py``
(``build_gaps.tsv``), this downloads each source, extracts its
accession.version set from the FASTA headers, and reports per-source *marginal*
coverage of the missing/superseded accessions. Use it to right-size the
``[[seqset]]`` history lists in ``sources.toml`` before committing to a
multi-hour full build — keep the sources that pull their weight, drop the
redundant ones.

Caching is provenance-keyed, not name-keyed: different sources routinely ship
identical FASTA basenames with different contents, so the cache key is the full
source descriptor (the URL), stored under ``<cache>/<label>.acc`` (the tiny
sorted accession.version list). The large ``.gz`` is deleted right after header
extraction unless ``--keep-gz`` is set, so peak disk stays small.

Example:
    uv run python seqrepo_equivalence/probe_source_coverage.py \\
        --sources sources.tsv \\          # tab-separated: url<TAB>label
        --gap runs/DATE-parity/build_gaps.tsv \\
        --prefixes NM_,NR_,XM_,XR_ \\
        --report coverage_report.json

``sources.tsv`` order matters: marginal coverage is reported in file order, so
list the sources you most expect to keep first (or oldest-first to see when each
annotation era stops adding anything).
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DEFAULT_RUN_DIR = REPO_ROOT / "runs" / "2026-07-02-seqrepo-parity"


def read_sources(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        url, _, label = line.partition("\t")
        out.append((url.strip(), label.strip() or url.rsplit("/", 1)[-1]))
    return out


def accession_set(gz_path: Path) -> set[str]:
    """First whitespace token of every FASTA header (the accession.version)."""
    accs: set[str] = set()
    with gzip.open(gz_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                accs.add(line[1:].split(None, 1)[0])
    return accs


def fetch_acc(url: str, label: str, cache: Path, keep_gz: bool) -> set[str]:
    """Return the source's accession.version set, caching the .acc list."""
    acc_file = cache / f"{label}.acc"
    if acc_file.exists():
        return set(acc_file.read_text().split())
    gz = cache / f"{label}.gz"
    tmp = gz.with_suffix(".gz.part")
    urllib.request.urlretrieve(url, tmp)  # noqa: S310
    tmp.replace(gz)
    accs = accession_set(gz)
    acc_file.write_text("\n".join(sorted(accs)) + "\n")
    if not keep_gz:
        gz.unlink(missing_ok=True)
    return accs


def load_gap(path: Path, prefixes: tuple[str, ...], namespace: str,
             categories: tuple[str, ...]) -> set[str]:
    """Missing/superseded accession.versions from the verifier gap feed.

    Columns: seqrepo_ns, refget_ns, accession, category, prefix, seqrepo_seq_id.
    """
    missing: set[str] = set()
    lines = path.read_text().splitlines()
    for line in lines[1:]:  # skip header
        c = line.split("\t")
        if len(c) < 5:
            continue
        if c[1] == namespace and c[3] in categories and c[2].startswith(prefixes):
            missing.add(c[2])
    return missing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sources", type=Path, required=True,
                    help="TSV of url<TAB>label, one per line (order = marginal order)")
    ap.add_argument("--gap", type=Path, default=DEFAULT_RUN_DIR / "build_gaps.tsv")
    ap.add_argument("--cache-dir", type=Path, default=REPO_ROOT / ".probe_cache")
    ap.add_argument("--prefixes", default="NM_,NR_,XM_,XR_",
                    help="comma-separated accession prefixes that count as in-scope")
    ap.add_argument("--namespace", default="refseq",
                    help="refget_ns column value to filter the gap feed on")
    ap.add_argument("--categories", default="sequence_missing,superseded_old_version")
    ap.add_argument("--report", type=Path, default=None)
    ap.add_argument("--keep-gz", action="store_true",
                    help="keep downloaded .gz files (default: delete after extract)")
    ap.add_argument("--limit", type=int, default=None,
                    help="only process the first N sources (smoke test)")
    args = ap.parse_args()

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    prefixes = tuple(p for p in args.prefixes.split(",") if p)
    categories = tuple(c for c in args.categories.split(",") if c)

    sources = read_sources(args.sources)
    if args.limit:
        sources = sources[: args.limit]
    gap = load_gap(args.gap, prefixes, args.namespace, categories)
    print(f"gap ({args.namespace} {'/'.join(prefixes)} "
          f"{'+'.join(categories)}): {len(gap):,}", file=sys.stderr)

    rows = []
    union: set[str] = set()
    seen_cov: set[str] = set()
    print(f"{'label':24} {'own':>10} {'marginal':>10} {'cum_cov':>10} {'cum_%':>7}")
    for url, label in sources:
        accs = fetch_acc(url, label, args.cache_dir, args.keep_gz)
        cov = accs & gap
        marginal = len(cov - seen_cov)
        seen_cov |= cov
        union |= accs
        cum_pct = 100 * len(seen_cov) / max(1, len(gap))
        print(f"{label:24} {len(accs):10,} {marginal:10,} "
              f"{len(seen_cov):10,} {cum_pct:6.1f}%")
        rows.append({"label": label, "url": url, "own": len(accs),
                     "marginal_gap": marginal, "covers_gap": len(cov)})

    residual = gap - union
    from collections import Counter
    by_prefix = dict(Counter(a.split("_", 1)[0] + "_" for a in residual))
    summary = {
        "gap_total": len(gap),
        "covered": len(seen_cov),
        "covered_pct": round(100 * len(seen_cov) / max(1, len(gap)), 1),
        "residual": len(residual),
        "residual_by_prefix": by_prefix,
        "sources": rows,
    }
    print(f"\nCOVERED {summary['covered']:,}/{len(gap):,} "
          f"({summary['covered_pct']}%)  RESIDUAL {len(residual):,} {by_prefix}")
    if args.report:
        args.report.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"wrote {args.report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
