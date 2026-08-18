#!/usr/bin/env python
"""Compare the experimental Ensembl collection to the preserved store."""

from __future__ import annotations

import csv
import argparse
import json
from collections import Counter
from pathlib import Path

from gtars.refget import RefgetStore

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]


def metadata_digest(store: RefgetStore, namespace: str, alias: str) -> str | None:
    metadata = store.get_sequence_metadata_by_alias(namespace, alias)
    return metadata.sha512t24u if metadata is not None else None


def compare_sequence_content(
    left_store: RefgetStore,
    left_digest: str,
    right_store: RefgetStore,
    right_digest: str,
) -> tuple[int, dict[str, int]]:
    """Count differing characters without materializing whole chromosomes."""
    differences = 0
    substitutions: Counter[str] = Counter()
    left_chunks = left_store.stream_sequence(left_digest)
    right_chunks = right_store.stream_sequence(right_digest)
    for left, right in zip(left_chunks, right_chunks, strict=True):
        if left == right:
            continue
        for left_base, right_base in zip(left, right, strict=True):
            if left_base != right_base:
                differences += 1
                substitutions[f"{left_base}->{right_base}"] += 1
    return differences, dict(sorted(substitutions.items()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--playground", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()

    lock = json.loads(args.lock.read_text())
    source = lock["sources"][0]
    collection_digest = source["collection_digest"]

    playground = RefgetStore.on_disk(str(args.playground))
    baseline = RefgetStore.on_disk(str(args.baseline))
    playground.load_collection(collection_digest)
    level2 = playground.get_collection_level2(collection_digest)

    rows: list[dict[str, object]] = []
    counts: Counter[str] = Counter()
    aggregate_substitutions: Counter[str] = Counter()
    chromosome_names = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

    for name, seq_id in zip(level2["names"], level2["sequences"], strict=True):
        digest = seq_id.removeprefix("SQ.")
        ensembl_metadata = playground.get_sequence_metadata(digest)
        try:
            baseline_metadata = baseline.get_sequence_metadata(digest)
        except Exception:  # a missing digest may be reported as an exception
            baseline_metadata = None

        ensembl_alias_metadata = playground.get_sequence_metadata_by_alias(
            "ensembl", name
        )
        grch38_metadata = baseline.get_sequence_metadata_by_alias("GRCh38", name)
        ensembl_digest = (
            ensembl_alias_metadata.sha512t24u if ensembl_alias_metadata else None
        )
        grch38_digest = grch38_metadata.sha512t24u if grch38_metadata else None
        baseline_present = baseline_metadata is not None
        ensembl_alias_matches = ensembl_digest == digest
        grch38_name_status = (
            "same" if grch38_digest == digest
            else "missing" if grch38_digest is None
            else "different"
        )
        difference_count = 0
        substitutions: dict[str, int] = {}
        if grch38_name_status == "different":
            difference_count, substitutions = compare_sequence_content(
                playground,
                digest,
                baseline,
                grch38_digest,
            )
            aggregate_substitutions.update(substitutions)

        counts["records"] += 1
        counts["baseline_digest_present" if baseline_present
               else "baseline_digest_missing"] += 1
        counts["ensembl_alias_matches" if ensembl_alias_matches
               else "ensembl_alias_mismatch"] += 1
        counts[f"grch38_name_{grch38_name_status}"] += 1
        if name in chromosome_names:
            counts["chromosome_records"] += 1
            counts[f"chromosome_grch38_name_{grch38_name_status}"] += 1

        rows.append({
            "name": name,
            "digest": digest,
            "ensembl_length": ensembl_metadata.length,
            "ensembl_md5": ensembl_metadata.md5,
            "baseline_digest_present": baseline_present,
            "ensembl_alias_matches": ensembl_alias_matches,
            "grch38_name_status": grch38_name_status,
            "grch38_digest": grch38_digest or "",
            "grch38_length": grch38_metadata.length if grch38_metadata else "",
            "grch38_md5": grch38_metadata.md5 if grch38_metadata else "",
            "same_length": (
                ensembl_metadata.length == grch38_metadata.length
                if grch38_metadata else False
            ),
            "character_differences": difference_count,
            "ensembl_to_grch38_substitutions": json.dumps(
                substitutions, separators=(",", ":")
            ),
        })

    with args.out_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "collection_digest": collection_digest,
        "source_url": source["url"],
        "source_sha256": source["sha256"],
        "counts": dict(sorted(counts.items())),
        "total_character_differences": sum(aggregate_substitutions.values()),
        "ensembl_to_grch38_substitutions": dict(
            sorted(aggregate_substitutions.items())
        ),
        "missing_from_baseline": [
            row["name"] for row in rows if not row["baseline_digest_present"]
        ],
        "grch38_name_differences": [
            row["name"] for row in rows
            if row["grch38_name_status"] == "different"
        ],
        "grch38_name_missing_examples": [
            row["name"] for row in rows
            if row["grch38_name_status"] == "missing"
        ][:20],
    }
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
