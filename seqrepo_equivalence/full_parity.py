#!/usr/bin/env python
"""Write exhaustive alias parity, gap, mismatch, and source-contribution tables."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

from gtars.refget import RefgetStore

from gks_refgetstore import build_lock
from gks_refgetstore.store_census import collection_members

from verify_seqrepo_equivalence import (
    DEFAULT_SEQREPO, EXPECTED_OMIT, NAMESPACE_MAP, build_refget_digest_set,
    load_refget_aliases, open_seqrepo_db,
)


def url_provider(url: str) -> str:
    """Name the publisher a source URL came from.

    Matched on host, not substring: an "ensembl in url else ncbi" test silently
    files every non-Ensembl host under ncbi, which mislabels the EBI LRG bundle.
    """
    for host, provider in (
        ("ftp.ensembl.org", "ensembl"),
        ("ftp.ncbi.nlm.nih.gov", "ncbi"),
        ("ftp.ebi.ac.uk", "ebi"),
    ):
        if host in url:
            return provider
    return "other"


def source_release(source: dict) -> str:
    owner = source["owner"]
    match = re.fullmatch(r"ensembl_release_(\d+)", owner)
    if match:
        return match.group(1)
    url = source["url"]
    annotation = re.search(r"/annotation_releases/9606/([^/]+)/", url)
    if annotation:
        return annotation.group(1)
    if source["kind"].startswith("assembly_"):
        return source["owner"]
    if source["owner"].startswith("refseq_human_"):
        return "current"
    return ""


def source_chronology(source: dict) -> tuple:
    """Order official releases oldest-first, leaving mutable feeds last."""
    owner = source["owner"]
    ensembl = re.fullmatch(r"ensembl_release_(\d+)", owner)
    if ensembl:
        return (1, int(ensembl.group(1)), 0, source["url"])
    if source["kind"].startswith("assembly_"):
        assembly = re.fullmatch(r"GRCh(37|38)(?:\.p(\d+))?", owner)
        if assembly:
            return (
                0, 0, int(assembly.group(1)), int(assembly.group(2) or 0),
                source["url"],
            )
        return (0, 0, 99, 0, source["url"])
    annotation = re.search(r"/annotation_releases/9606/([^/]+)/", source["url"])
    if annotation:
        numbers = tuple(int(value) for value in re.findall(r"\d+", annotation.group(1)))
        return (0, 1, numbers, source["url"])
    # Species-level current feeds are mutable and therefore never an "earliest"
    # source when a historical NCBI source contains the same digest.
    return (0, 2, source["url"])


def source_contributions(
    store: RefgetStore, lock_path: Path, seqrepo_digests: set[str], out_path: Path,
    summary_path: Path,
) -> tuple[int, Counter]:
    # Through load_lock, not json.loads: this reader attributes every store-only
    # digest to the *earliest* source that introduced it, and a silently
    # mis-shaped lock would produce a complete-looking table attributing
    # nothing.
    lock = build_lock.load_lock(lock_path)
    collection_of = build_lock.collection_by_file(lock)
    first: dict[str, dict] = {}
    for source in sorted(build_lock.lock_files(lock), key=source_chronology):
        collection = collection_of.get(source["cache_path"])
        if not collection:
            continue
        for digest in collection_members(store, collection):
            if digest not in seqrepo_digests and digest not in first:
                first[digest] = source
    grouped: Counter = Counter()
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("digest", "provider", "owner", "release", "file_class", "url"))
        for digest, source in sorted(first.items()):
            provider = url_provider(source["url"])
            release = source_release(source)
            file_class = source.get("file_class") or ""
            grouped[(provider, release, file_class)] += 1
            writer.writerow((
                digest, provider, source["owner"], release,
                file_class, source["url"],
            ))
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("provider", "release", "file_class", "store_only_digests"))
        for (provider, release, file_class), count in sorted(grouped.items()):
            writer.writerow((provider, release, file_class, count))
    return len(first), grouped


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--store", type=Path, default=Path("store"))
    parser.add_argument("--seqrepo", type=Path, default=DEFAULT_SEQREPO)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    args = parser.parse_args()
    args.run_dir.mkdir(parents=True, exist_ok=True)

    store = RefgetStore.open_local(str(args.store))
    store_digests = build_refget_digest_set(store)
    conn = open_seqrepo_db(args.seqrepo)
    all_seqrepo_digests = {
        row[0] for row in conn.execute(
            "SELECT DISTINCT seq_id FROM seqalias WHERE is_current=1"
        )
    }
    placeholders = ",".join("?" for _ in EXPECTED_OMIT)
    rows = conn.execute(
        f"SELECT namespace, alias, seq_id FROM seqalias WHERE is_current=1 "
        f"AND namespace NOT IN ({placeholders})",
        tuple(sorted(EXPECTED_OMIT)),
    )
    seqrepo: dict[tuple[str, str], tuple[str, str]] = {}
    seqrepo_namespace_counts: Counter[str] = Counter()
    for namespace, alias, digest in rows:
        canonical = NAMESPACE_MAP.get(namespace, namespace)
        key = (canonical, alias)
        previous = seqrepo.get(key)
        if previous is not None and previous[1] != digest:
            raise ValueError(f"SeqRepo canonical alias collision {canonical}:{alias}")
        seqrepo[key] = (namespace, digest)
        seqrepo_namespace_counts[canonical] += 1

    alias_path = args.run_dir / "parity_by_alias.tsv"
    gaps_path = args.run_dir / "build_gaps.tsv"
    mismatch_path = args.run_dir / "shared_alias_mismatches.tsv"
    counts: Counter[str] = Counter()
    store_namespace_counts: Counter[str] = Counter()
    with (
        alias_path.open("w", encoding="utf-8", newline="") as alias_handle,
        gaps_path.open("w", encoding="utf-8", newline="") as gap_handle,
        mismatch_path.open("w", encoding="utf-8", newline="") as mismatch_handle,
    ):
        aliases = csv.writer(alias_handle, delimiter="\t", lineterminator="\n")
        gaps = csv.writer(gap_handle, delimiter="\t", lineterminator="\n")
        mismatches = csv.writer(mismatch_handle, delimiter="\t", lineterminator="\n")
        header = (
            "namespace", "alias", "seqrepo_namespace", "seqrepo_digest",
            "store_namespace", "store_digest", "status",
        )
        aliases.writerow(header)
        gaps.writerow(header)
        mismatches.writerow(header)
        for namespace in sorted(store.list_sequence_alias_namespaces()):
            for alias, store_digest in load_refget_aliases(args.store, namespace).items():
                store_namespace_counts[namespace] += 1
                sr = seqrepo.pop((namespace, alias), None)
                if sr is None:
                    status = "store_only_alias"
                    row = (namespace, alias, "", "", namespace, store_digest, status)
                else:
                    sr_namespace, sr_digest = sr
                    status = "match" if sr_digest == store_digest else "digest_mismatch"
                    row = (
                        namespace, alias, sr_namespace, sr_digest,
                        namespace, store_digest, status,
                    )
                    if status == "digest_mismatch":
                        mismatches.writerow(row)
                        gaps.writerow(row)
                counts[status] += 1
                aliases.writerow(row)

        for (namespace, alias), (sr_namespace, sr_digest) in sorted(seqrepo.items()):
            status = (
                "seqrepo_only_alias_digest_present"
                if sr_digest in store_digests else "seqrepo_only_missing_sequence"
            )
            row = (namespace, alias, sr_namespace, sr_digest, "", "", status)
            counts[status] += 1
            aliases.writerow(row)
            gaps.writerow(row)

    contribution_count, contribution_groups = source_contributions(
        store, args.lock, all_seqrepo_digests,
        args.run_dir / "store_only_contributions.tsv",
        args.run_dir / "source_contribution_summary.tsv",
    )
    shared_digests = len(store_digests & all_seqrepo_digests)
    summary = {
        "schema": "gks-refgetstore-full-parity/1",
        "seqrepo_path": str(args.seqrepo),
        "store_path": str(args.store),
        "seqrepo_biological_aliases": sum(
            counts[key] for key in (
                "match", "digest_mismatch", "seqrepo_only_alias_digest_present",
                "seqrepo_only_missing_sequence",
            )
        ),
        "store_alias_rows": counts["match"] + counts["digest_mismatch"] + counts["store_only_alias"],
        "digest_membership": {
            "seqrepo": len(all_seqrepo_digests),
            "store": len(store_digests),
            "both": shared_digests,
            "seqrepo_only": len(all_seqrepo_digests) - shared_digests,
            "store_only": len(store_digests) - shared_digests,
        },
        "seqrepo_namespace_alias_counts": dict(sorted(seqrepo_namespace_counts.items())),
        "store_namespace_alias_counts": dict(sorted(store_namespace_counts.items())),
        "status_counts": dict(sorted(counts.items())),
        "store_only_contributions": contribution_count,
        "source_contribution_groups": {
            "|".join(key): value
            for key, value in sorted(contribution_groups.items())
        },
    }
    (args.run_dir / "full-parity-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Full alias and source parity", "", "## Alias status", "",
        "| status | rows |", "|---|---:|",
        *(f"| `{key}` | {value:,} |" for key, value in sorted(counts.items())),
        "", "## Digest membership", "", "| membership | digests |", "|---|---:|",
        *(f"| `{key}` | {value:,} |" for key, value in summary["digest_membership"].items()),
        "", "## Store namespace aliases", "", "| namespace | aliases |", "|---|---:|",
        *(f"| `{key}` | {value:,} |" for key, value in sorted(store_namespace_counts.items())),
        "", f"Store-only digests attributed to an earliest source: {contribution_count:,}.",
    ]
    (args.run_dir / "full-parity-summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
