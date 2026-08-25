#!/usr/bin/env python
"""Byte-diagnose every shared-alias digest mismatch from ``full_parity.py``."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from biocommons.seqrepo import SeqRepo
from bioutils.digests import seq_seqhash
from ga4gh.core.digests import sha512t24u
from gtars.refget import RefgetStore


def difference_profile(left: str, right: str) -> tuple[int | None, int, Counter]:
    first: int | None = None
    substitutions: Counter = Counter()
    count = 0
    for index, (a, b) in enumerate(zip(left, right)):
        if a == b:
            continue
        if first is None:
            first = index
        count += 1
        substitutions[f"{a}>{b}"] += 1
    if len(left) != len(right):
        if first is None:
            first = min(len(left), len(right))
        count += abs(len(left) - len(right))
        substitutions["length_delta"] += abs(len(left) - len(right))
    return first, count, substitutions


def mismatch_cause(
    row: dict, sr_sequence: str, rg_sequence: str, raw_digest: str,
    normalized_digest: str, substitutions: Counter,
) -> str:
    if (
        sr_sequence == rg_sequence
        and "*" in sr_sequence
        and raw_digest == row["store_digest"]
        and normalized_digest == row["seqrepo_digest"]
    ):
        return "ensembl_stop_codon_digest_normalization"
    ambiguity = set("BDHKMRSVWY")
    known_pairs = {
        f"{base}>N" for base in ambiguity
    } | {
        f"N>{base}" for base in ambiguity
    }
    if (
        row["namespace"] == "ensembl"
        and len(sr_sequence) == len(rg_sequence)
        and substitutions
        and set(substitutions) <= known_pairs
    ):
        return "ensembl_iupac_to_n_normalization"
    if sr_sequence != rg_sequence:
        return "sequence_divergence_unexplained"
    return "unexplained_digest_algorithm"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mismatches", type=Path, required=True)
    parser.add_argument("--store", type=Path, default=Path("store"))
    parser.add_argument("--seqrepo", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()

    seqrepo = SeqRepo(str(args.seqrepo))
    store = RefgetStore.open_local(str(args.store))
    records: list[dict] = []
    with args.mismatches.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for index, row in enumerate(rows, 1):
        record = dict(row)
        try:
            sr_sequence = seqrepo.fetch_uri(
                f"{row['seqrepo_namespace']}:{row['alias']}"
            )
            metadata = store.get_sequence_metadata(row["store_digest"])
            store.load_sequence(row["store_digest"])
            rg_sequence = store.get_substring(
                row["store_digest"], 0, metadata.length
            )
            first, difference_count, substitutions = difference_profile(
                sr_sequence, rg_sequence
            )
            raw_digest = sha512t24u(sr_sequence.encode("ascii"))
            normalized_digest = seq_seqhash(sr_sequence, normalize=True)
            record.update({
                "seqrepo_length": len(sr_sequence),
                "store_length": len(rg_sequence),
                "sequences_identical": sr_sequence == rg_sequence,
                "first_difference": first,
                "differing_positions": difference_count,
                "substitutions": json.dumps(dict(sorted(substitutions.items()))),
                "seqrepo_raw_digest": raw_digest,
                "store_raw_digest": sha512t24u(rg_sequence.encode("ascii")),
                "seqrepo_normalized_digest": normalized_digest,
                "stop_codon_count": sr_sequence.count("*") + rg_sequence.count("*"),
                "cause": mismatch_cause(
                    row, sr_sequence, rg_sequence, raw_digest,
                    normalized_digest, substitutions,
                ),
                "error": "",
            })
        except Exception as exc:  # noqa: BLE001
            record.update({
                "seqrepo_length": "", "store_length": "",
                "sequences_identical": "", "first_difference": "",
                "differing_positions": "", "substitutions": "",
                "seqrepo_raw_digest": "", "store_raw_digest": "",
                "seqrepo_normalized_digest": "", "stop_codon_count": "",
                "cause": "fetch_error", "error": str(exc),
            })
        records.append(record)
        if index % 100 == 0 or index == len(rows):
            print(f"diagnosed {index}/{len(rows)}")

    fields = list(records[0]) if records else [
        "namespace", "alias", "seqrepo_namespace", "seqrepo_digest",
        "store_namespace", "store_digest", "status", "seqrepo_length",
        "store_length", "sequences_identical", "first_difference",
        "differing_positions", "substitutions", "seqrepo_raw_digest",
        "store_raw_digest", "seqrepo_normalized_digest", "stop_codon_count",
        "cause", "error",
    ]
    with args.out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)
    causes = Counter(record["cause"] for record in records)
    summary = {
        "schema": "gks-refgetstore-mismatch-diagnosis/1",
        "mismatch_count": len(records),
        "cause_counts": dict(sorted(causes.items())),
        "unexplained_count": (
            causes["fetch_error"]
            + causes["unexplained_digest_algorithm"]
            + causes["sequence_divergence_unexplained"]
        ),
    }
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return int(summary["unexplained_count"] != 0)


if __name__ == "__main__":
    raise SystemExit(main())
