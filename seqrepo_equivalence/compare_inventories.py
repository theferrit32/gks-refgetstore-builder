#!/usr/bin/env python
"""Stepwise inventory comparison: a seqrepo snapshot vs a gtars RefgetStore.

Answers two questions completely, rather than in summary:

  1. Which *sequences* does each side hold that the other does not?
  2. Which *alias namespaces* does each side expose, which correspond, and why
     do the non-corresponding ones exist?

Each step states what it is comparing before it compares it, and writes the full
listing to TSV -- no sampling, no truncation. Counts printed to stdout are
derived from the same rows written to disk.

Comparison is by **sha512t24u digest**, because the digest is the sequence: a
sequence counts as shared no matter what either side calls it. Alias-level
agreement is a separate question, answered in step 5 for the namespaces the two
sides have in common.

Both sides are read directly and read-only: seqrepo from ``aliases.sqlite3`` and
``sequences/db.sqlite3`` via stdlib ``sqlite3``, the store from its on-disk alias
TSVs plus gtars. Nothing is written outside ``--out-dir``.

Examples:
    uv run python seqrepo_equivalence/compare_inventories.py
    uv run python seqrepo_equivalence/compare_inventories.py \
        --seqrepo ~/dev/data/seqrepo/2024-12-20 --store ./store \
        --out-dir ./runs/DATE/inventory
    uv run python seqrepo_equivalence/compare_inventories.py --include-shared
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

from gtars.refget import RefgetStore

from verify_seqrepo_equivalence import (
    DEFAULT_SEQREPO,
    DEFAULT_STORE,
    EXPECTED_OMIT,
    NAMESPACE_MAP,
    load_refget_aliases,
    open_seqrepo_db,
)

ALIAS_EXAMPLE_CAP = 8
ENSEMBL_RELEASE_RE = re.compile(r"^ensembl-\d+$")


# --------------------------------------------------------------------------- io


def write_tsv(path: Path, header: tuple[str, ...], rows) -> int:
    """Write ``rows`` under ``header``. Returns the row count."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
            n += 1
    return n


def step(number: int, title: str, comparing: str) -> None:
    print(f"\n=== Step {number}: {title} ===")
    print(f"    comparing: {comparing}", flush=True)


def note(path: Path, rows: int, out_dir: Path) -> None:
    print(f"    -> {path.relative_to(out_dir)}  ({rows:,} rows)", flush=True)


# --------------------------------------------------------------------------- seqrepo


def seqrepo_namespace_counts(conn: sqlite3.Connection) -> dict[str, tuple[int, int]]:
    """namespace -> (alias rows, distinct digests), current aliases only."""
    return {
        ns: (rows, digests)
        for ns, rows, digests in conn.execute(
            "SELECT namespace, COUNT(*), COUNT(DISTINCT seq_id) FROM seqalias "
            "WHERE is_current=1 GROUP BY namespace"
        )
    }


def seqrepo_digests(conn: sqlite3.Connection) -> set[str]:
    return {row[0] for row in conn.execute(
        "SELECT DISTINCT seq_id FROM seqalias WHERE is_current=1")}


def seqrepo_seqinfo(seqrepo_root: Path, wanted: set[str]) -> dict[str, tuple[int, str]]:
    """digest -> (length, alphabet) for ``wanted``, from sequences/db.sqlite3.

    Absent in a partial snapshot; the caller degrades to blank columns rather
    than failing, since membership is still answerable without it.
    """
    db = seqrepo_root / "sequences" / "db.sqlite3"
    if not db.exists():
        print(f"    WARNING: {db} not found; length/alphabet will be blank")
        return {}
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    return {
        seq_id: (length, alpha)
        for seq_id, length, alpha in conn.execute("SELECT seq_id, len, alpha FROM seqinfo")
        if seq_id in wanted
    }


def seqrepo_alias_detail(
    conn: sqlite3.Connection, wanted: set[str], biological_only: bool
) -> tuple[dict[str, list[str]], dict[str, set[str]], dict[str, int]]:
    """For digests in ``wanted``: example aliases, namespaces, and alias counts."""
    aliases: dict[str, list[str]] = defaultdict(list)
    namespaces: dict[str, set[str]] = defaultdict(set)
    counts: dict[str, int] = defaultdict(int)
    query = "SELECT seq_id, namespace, alias FROM seqalias WHERE is_current=1"
    if biological_only:
        placeholders = ",".join("?" for _ in EXPECTED_OMIT)
        query += f" AND namespace NOT IN ({placeholders})"
        cursor = conn.execute(query, tuple(sorted(EXPECTED_OMIT)))
    else:
        cursor = conn.execute(query)
    for seq_id, namespace, alias in cursor:
        if seq_id not in wanted:
            continue
        namespaces[seq_id].add(namespace)
        counts[seq_id] += 1
        if len(aliases[seq_id]) < ALIAS_EXAMPLE_CAP:
            aliases[seq_id].append(f"{namespace}:{alias}")
    return aliases, namespaces, counts


# --------------------------------------------------------------------------- store


def store_namespace_aliases(store_path: Path, namespaces: list[str]) -> dict[str, dict[str, str]]:
    """namespace -> {alias: digest}, read from the on-disk TSV sidecars."""
    return {ns: load_refget_aliases(store_path, ns) for ns in namespaces}


def store_alias_detail(
    per_namespace: dict[str, dict[str, str]], wanted: set[str]
) -> tuple[dict[str, list[str]], dict[str, set[str]], dict[str, int]]:
    aliases: dict[str, list[str]] = defaultdict(list)
    namespaces: dict[str, set[str]] = defaultdict(set)
    counts: dict[str, int] = defaultdict(int)
    for ns, mapping in per_namespace.items():
        for alias, digest in mapping.items():
            if digest not in wanted:
                continue
            namespaces[digest].add(ns)
            counts[digest] += 1
            if len(aliases[digest]) < ALIAS_EXAMPLE_CAP:
                aliases[digest].append(f"{ns}:{alias}")
    return aliases, namespaces, counts


# --------------------------------------------------------------------------- namespaces


def classify_namespace(canonical: str, in_seqrepo: bool, in_store: bool,
                       seqrepo_names: list[str]) -> tuple[str, str]:
    """Return (status, reason) for one canonical namespace."""
    if in_seqrepo and in_store:
        return "both", "corresponding namespaces compared alias-by-alias in step 5"
    if in_seqrepo:
        if set(seqrepo_names) & EXPECTED_OMIT:
            return ("seqrepo_only",
                    "digest-synthesis namespace omitted from the store by design: "
                    "these are 1:1 re-encodings of the sequence digest, derivable "
                    "on demand, so storing them would duplicate the primary key")
        return ("seqrepo_only",
                "assembly or source the build deliberately does not load; the "
                "sequences may still be present by digest (see step 4)")
    if ENSEMBL_RELEASE_RE.match(canonical):
        return ("store_only",
                "immutable per-release Ensembl namespace; seqrepo keeps one "
                "rolling Ensembl view, so it has no counterpart")
    if canonical == "lrg":
        return "store_only", "EBI LRG; the seqrepo snapshot has no LRG namespace"
    return "store_only", "source loaded by the build that seqrepo does not carry"


# --------------------------------------------------------------------------- main


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seqrepo", type=Path, default=DEFAULT_SEQREPO,
                        help="seqrepo snapshot directory (holds aliases.sqlite3)")
    parser.add_argument("--store", type=Path, default=DEFAULT_STORE,
                        help="RefgetStore directory")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="directory for the TSV outputs")
    parser.add_argument("--include-shared", action="store_true",
                        help="also write the shared-sequence listing (large)")
    args = parser.parse_args()

    for path, what in ((args.seqrepo, "seqrepo snapshot"), (args.store, "store")):
        if not path.exists():
            raise SystemExit(f"{what} not found: {path}")
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    summary: list[tuple[str, str, str]] = []

    def record(step_name: str, metric: str, value) -> None:
        summary.append((step_name, metric, f"{value:,}" if isinstance(value, int) else str(value)))

    print(f"seqrepo : {args.seqrepo}")
    print(f"store   : {args.store}")
    print(f"out-dir : {out}")

    # ---------------------------------------------------------------- step 1
    step(1, "Inventory the seqrepo snapshot",
         "every current alias row in aliases.sqlite3, grouped by namespace")
    conn = open_seqrepo_db(args.seqrepo)
    sr_ns = seqrepo_namespace_counts(conn)
    sr_digests = seqrepo_digests(conn)
    rows = sorted(
        ((ns, r, d, "digest_synthesis" if ns in EXPECTED_OMIT else "biological")
         for ns, (r, d) in sr_ns.items()), key=lambda x: -x[1])
    n = write_tsv(out / "01_seqrepo_namespaces.tsv",
                  ("namespace", "alias_rows", "distinct_digests", "namespace_kind"), rows)
    note(out / "01_seqrepo_namespaces.tsv", n, out)
    print(f"    {len(sr_ns)} namespaces, {sum(r for r, _ in sr_ns.values()):,} alias rows, "
          f"{len(sr_digests):,} distinct sequence digests")
    record("1_seqrepo_inventory", "namespaces", len(sr_ns))
    record("1_seqrepo_inventory", "alias_rows", sum(r for r, _ in sr_ns.values()))
    record("1_seqrepo_inventory", "distinct_digests", len(sr_digests))

    # ---------------------------------------------------------------- step 2
    step(2, "Inventory the RefgetStore",
         "every alias in the store's aliases/sequences/*.tsv sidecars, "
         "plus its distinct sequence digests")
    store = RefgetStore.open_local(str(args.store))
    store_ns_names = sorted(store.list_sequence_alias_namespaces())
    per_ns = store_namespace_aliases(args.store, store_ns_names)
    st_digests: set[str] = set()
    st_meta: dict[str, tuple[int, str]] = {}
    for md in store.list_sequences():
        st_digests.add(md.sha512t24u)
        if md.sha512t24u not in sr_digests:      # only what step 4 will need
            st_meta[md.sha512t24u] = (md.length, md.alphabet)
    rows = sorted(((ns, len(m), len(set(m.values()))) for ns, m in per_ns.items()),
                  key=lambda x: -x[1])
    n = write_tsv(out / "02_store_namespaces.tsv",
                  ("namespace", "alias_rows", "distinct_digests"), rows)
    note(out / "02_store_namespaces.tsv", n, out)
    total_alias_rows = sum(len(m) for m in per_ns.values())
    print(f"    {len(store_ns_names)} namespaces, {total_alias_rows:,} alias rows, "
          f"{len(st_digests):,} distinct sequence digests")
    record("2_store_inventory", "namespaces", len(store_ns_names))
    record("2_store_inventory", "alias_rows", total_alias_rows)
    record("2_store_inventory", "distinct_digests", len(st_digests))

    # ---------------------------------------------------------------- step 3
    step(3, "Compare alias namespaces",
         "seqrepo namespace names mapped through NAMESPACE_MAP, then set-compared "
         "against the store's namespace names")
    print("    NAMESPACE_MAP (seqrepo -> store):")
    for src, dst in sorted(NAMESPACE_MAP.items()):
        print(f"      {src:12s} -> {dst}")
    canonical_to_seqrepo: dict[str, list[str]] = defaultdict(list)
    for ns in sr_ns:
        canonical_to_seqrepo[NAMESPACE_MAP.get(ns, ns)].append(ns)
    store_set = set(store_ns_names)
    ns_rows = []
    for canonical in sorted(set(canonical_to_seqrepo) | store_set):
        src_names = sorted(canonical_to_seqrepo.get(canonical, []))
        in_sr, in_st = bool(src_names), canonical in store_set
        status, reason = classify_namespace(canonical, in_sr, in_st, src_names)
        ns_rows.append((
            canonical, ";".join(src_names), canonical if in_st else "", status,
            sum(sr_ns[s][0] for s in src_names),
            len(per_ns[canonical]) if in_st else 0, reason,
        ))
    n = write_tsv(out / "03_namespace_comparison.tsv",
                  ("canonical_namespace", "seqrepo_namespace", "store_namespace",
                   "status", "seqrepo_alias_rows", "store_alias_rows", "reason"), ns_rows)
    note(out / "03_namespace_comparison.tsv", n, out)
    for status in ("both", "seqrepo_only", "store_only"):
        hits = [r[0] for r in ns_rows if r[3] == status]
        print(f"    {status:13s} {len(hits):3d}  {', '.join(hits[:6])}"
              f"{' …' if len(hits) > 6 else ''}")
        record("3_namespaces", status, len(hits))

    # ---------------------------------------------------------------- step 4
    step(4, "Compare sequences by digest",
         "the set of sha512t24u digests each side holds, independent of alias "
         "spelling; complete listings written for both directions")
    both = sr_digests & st_digests
    sr_only = sr_digests - st_digests
    st_only = st_digests - sr_digests
    print(f"    in both          {len(both):,}")
    print(f"    seqrepo only     {len(sr_only):,}")
    print(f"    store only       {len(st_only):,}")
    record("4_sequences", "in_both", len(both))
    record("4_sequences", "seqrepo_only", len(sr_only))
    record("4_sequences", "store_only", len(st_only))

    print("    resolving seqrepo-only detail …", flush=True)
    sr_info = seqrepo_seqinfo(args.seqrepo, sr_only)
    sr_al, sr_nsmap, sr_cnt = seqrepo_alias_detail(conn, sr_only, biological_only=False)
    rows = []
    for digest in sorted(sr_only):
        length, alpha = sr_info.get(digest, ("", ""))
        names = sr_nsmap.get(digest, set())
        biological = sorted(names - EXPECTED_OMIT)
        rows.append((
            digest, length, alpha, sr_cnt.get(digest, 0),
            ";".join(sorted(names)), ";".join(biological),
            "yes" if biological else "no",
            ";".join(a for a in sr_al.get(digest, []) if not a.startswith(tuple(f"{o}:" for o in EXPECTED_OMIT))),
        ))
    n = write_tsv(out / "04_sequences_seqrepo_only.tsv",
                  ("digest", "length", "alphabet", "seqrepo_alias_rows",
                   "seqrepo_namespaces", "biological_namespaces",
                   "has_biological_accession", "example_biological_aliases"), rows)
    note(out / "04_sequences_seqrepo_only.tsv", n, out)
    nameable = sum(1 for r in rows if r[6] == "yes")
    print(f"      of these, {nameable:,} carry a biological accession and "
          f"{len(rows) - nameable:,} are digest-only (no accession anywhere in seqrepo)")
    record("4_sequences", "seqrepo_only_with_accession", nameable)
    record("4_sequences", "seqrepo_only_digest_only", len(rows) - nameable)

    print("    resolving store-only detail …", flush=True)
    st_al, st_nsmap, st_cnt = store_alias_detail(per_ns, st_only)
    rows = []
    for digest in sorted(st_only):
        length, alpha = st_meta.get(digest, ("", ""))
        rows.append((
            digest, length, alpha, st_cnt.get(digest, 0),
            ";".join(sorted(st_nsmap.get(digest, set()))),
            ";".join(st_al.get(digest, [])),
        ))
    n = write_tsv(out / "05_sequences_store_only.tsv",
                  ("digest", "length", "alphabet", "store_alias_rows",
                   "store_namespaces", "example_store_aliases"), rows)
    note(out / "05_sequences_store_only.tsv", n, out)

    if args.include_shared:
        print("    writing shared listing …", flush=True)
        n = write_tsv(out / "06_sequences_shared.tsv", ("digest",),
                      ((d,) for d in sorted(both)))
        note(out / "06_sequences_shared.tsv", n, out)

    # ---------------------------------------------------------------- step 5
    step(5, "Compare aliases within namespaces both sides have",
         "for each corresponding namespace pair, the alias strings themselves and "
         "whether a shared alias points at the same digest")
    shared_ns = [(r[1].split(";")[0], r[0]) for r in ns_rows if r[3] == "both"]
    per_ns_rows, diff_rows = [], []
    for sr_name, st_name in sorted(shared_ns):
        sr_map = {alias: seq for seq, alias in conn.execute(
            "SELECT seq_id, alias FROM seqalias WHERE is_current=1 AND namespace=?",
            (sr_name,))}
        st_map = per_ns[st_name]
        shared_aliases = set(sr_map) & set(st_map)
        agree = sum(1 for a in shared_aliases if sr_map[a] == st_map[a])
        only_sr, only_st = set(sr_map) - set(st_map), set(st_map) - set(sr_map)
        per_ns_rows.append((sr_name, st_name, len(sr_map), len(st_map),
                            len(shared_aliases), agree, len(shared_aliases) - agree,
                            len(only_sr), len(only_st)))
        for alias in sorted(shared_aliases):
            if sr_map[alias] != st_map[alias]:
                diff_rows.append((sr_name, st_name, alias, "digest_mismatch",
                                  sr_map[alias], st_map[alias]))
        for alias in sorted(only_sr):
            diff_rows.append((sr_name, st_name, alias, "seqrepo_only_alias",
                              sr_map[alias], ""))
        for alias in sorted(only_st):
            diff_rows.append((sr_name, st_name, alias, "store_only_alias",
                              "", st_map[alias]))
        print(f"    {sr_name:12s} -> {st_name:12s} shared={len(shared_aliases):>7,} "
              f"agree={agree:>7,} mismatch={len(shared_aliases) - agree:>5,} "
              f"seqrepo_only={len(only_sr):>7,} store_only={len(only_st):>7,}")
    n = write_tsv(out / "07_alias_comparison_by_namespace.tsv",
                  ("seqrepo_namespace", "store_namespace", "seqrepo_aliases",
                   "store_aliases", "shared_aliases", "shared_digest_agree",
                   "shared_digest_mismatch", "seqrepo_only_aliases",
                   "store_only_aliases"), per_ns_rows)
    note(out / "07_alias_comparison_by_namespace.tsv", n, out)
    n = write_tsv(out / "08_alias_differences.tsv",
                  ("seqrepo_namespace", "store_namespace", "alias", "status",
                   "seqrepo_digest", "store_digest"), diff_rows)
    note(out / "08_alias_differences.tsv", n, out)
    record("5_aliases", "namespace_pairs_compared", len(per_ns_rows))
    record("5_aliases", "difference_rows", len(diff_rows))

    # ---------------------------------------------------------------- summary
    n = write_tsv(out / "00_summary.tsv", ("step", "metric", "value"), summary)
    note(out / "00_summary.tsv", n, out)
    print("\ndone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
