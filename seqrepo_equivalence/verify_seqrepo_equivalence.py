#!/usr/bin/env python
"""Exhaustively compare a biocommons seqrepo snapshot against a gtars RefgetStore.

Answers: *is the RefgetStore backwards-compatible with seqrepo?* — i.e. are all the
sequences seqrepo knows present in the store, do their digests match, and do the
namespaces / aliases line up. Every difference is categorized, and the seqrepo-only
gaps are exported to a build-gap list that drives extending the build.

Design: the core comparison is **gtars + stdlib only**. seqrepo is read directly from
its ``aliases.sqlite3`` via stdlib ``sqlite3`` (no ``biocommons.seqrepo`` import); the
store's aliases are read from the on-disk ``aliases/sequences/<ns>.tsv`` sidecars
rather than ~1M ``get_sequence_by_alias`` calls. Only the optional ``--deep`` byte-level
diagnosis needs ``biocommons.seqrepo`` (and an env that has it + the sequence data).

Exit 0 if backwards-compatible (per --fail-on), 1 if not, 2 on setup error.

See README.md / the project plan for the namespace-mapping rationale.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sqlite3
import sys
import time
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path

from gtars.refget import RefgetStore

logger = logging.getLogger("verify_seqrepo_equivalence")

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent  # this script lives in seqrepo_equivalence/
DEFAULT_SEQREPO = Path("/Users/kferrite/dev/data/seqrepo/2024-12-20")
DEFAULT_STORE = REPO_ROOT / "store"
DEFAULT_KNOWN_DIVERGENT = HERE / "ensembl_known_divergent.txt"
DEFAULT_RUN_DIR = REPO_ROOT / "runs" / "2026-07-02-seqrepo-parity"
DEFAULT_REPORT = DEFAULT_RUN_DIR / "equivalence-report.json"
DEFAULT_GAP_LIST = DEFAULT_RUN_DIR / "build_gaps.tsv"

# seqrepo namespace -> refgetstore namespace, for namespaces that should be equivalent.
NAMESPACE_MAP: dict[str, str] = {
    "NCBI": "refseq",
    "Ensembl": "ensembl",
    "GRCh38": "GRCh38",
    "GRCh37": "GRCh37",
    "GRCh37.p13": "GRCh37.p13",
}

# Digest-synthesis namespaces the store intentionally omits (synthesized at query time).
# Backwards-compat for these is satisfied by digest-level coverage, not alias presence.
EXPECTED_OMIT: set[str] = {"MD5", "SEGUID", "SHA1", "VMC"}

# Old assemblies / patch levels the build deliberately does not load. Out of scope for
# alias equivalence; their sequences are still counted in digest coverage.
VERSION_DRIFT_NAMESPACES: set[str] = {
    "GRCh38.p1", "GRCh38.p2", "GRCh38.p3", "GRCh38.p4", "GRCh38.p5", "GRCh38.p6",
    "GRCh38.p7", "GRCh38.p8", "GRCh38.p9", "GRCh38.p10", "GRCh38.p11", "GRCh38.p12",
    "GRCh37.p2", "GRCh37.p5", "GRCh37.p9", "GRCh37.p10", "GRCh37.p11", "GRCh37.p12",
    "NCBI36", "NCBI35", "NCBI34", "hs37d5", "hs37-1kg", "JRGv1", "JRGv2", "CHM1_1.1",
}

# Accession-bearing seqrepo namespaces we care about for digest attribution
# (everything except the pure-digest synthesis namespaces).
ACCESSION_NAMESPACES_EXCLUDED = EXPECTED_OMIT


# --------------------------------------------------------------------------- models


@dataclass
class AliasComparison:
    seqrepo_ns: str
    refget_ns: str
    group: str = "core"  # "core" (explicit NAMESPACE_MAP) or "additional" (identity-mapped)
    seqrepo_alias_count: int = 0
    refget_alias_count: int = 0
    common: int = 0
    matched: int = 0
    # digest mismatches on the common set
    mismatch_total: int = 0
    mismatch_explained: int = 0
    mismatch_unexplained: int = 0
    mismatch_by_cause: dict[str, int] = field(default_factory=dict)
    # seqrepo-only aliases, categorized under the latest+best-effort-old policy
    superseded_old_version: int = 0  # store has a higher version; old one optional to backfill
    backfill_candidate: int = 0      # store is behind (has only lower versions of this base)
    # the seqrepo alias is absent from the store; split by whether the underlying
    # sequence (its seq_id digest) is in the store at all:
    alias_naming_gap: int = 0        # digest IS in store under a different alias string
    alias_naming_by_prefix: dict[str, int] = field(default_factory=dict)
    sequence_missing: int = 0        # digest absent from the store entirely
    sequence_missing_by_prefix: dict[str, int] = field(default_factory=dict)
    refget_only: int = 0             # store extensions not in seqrepo
    # examples (capped)
    examples_unexplained: list[dict] = field(default_factory=list)
    examples_alias_naming: list[dict] = field(default_factory=list)
    examples_sequence_missing: list[str] = field(default_factory=list)
    examples_backfill: list[dict] = field(default_factory=list)

    @property
    def backwards_compat_coverage_pct(self) -> float:
        """Fraction of seqrepo aliases that resolve in the store to a matching digest.
        Alias-naming gaps do NOT count (the alias string does not resolve), even though
        the underlying sequence is present."""
        if self.seqrepo_alias_count == 0:
            return 100.0
        ok = self.matched + self.mismatch_explained + self.superseded_old_version
        return round(100.0 * ok / self.seqrepo_alias_count, 4)


@dataclass
class DigestCoverage:
    seqrepo_distinct_digests: int = 0
    present_in_refget: int = 0
    missing: int = 0
    coverage_pct: float = 0.0
    missing_referenced_by_mapped: int = 0       # genuinely missing (a mapped current alias needs it)
    missing_version_drift_only: int = 0         # only old-patch / drift namespaces reference it
    missing_other: int = 0
    examples_missing_mapped: list[dict] = field(default_factory=list)


@dataclass
class NamespaceClassification:
    mapped: list[dict] = field(default_factory=list)
    expected_omit: list[dict] = field(default_factory=list)
    version_drift: list[dict] = field(default_factory=list)
    refget_only: list[str] = field(default_factory=list)
    unrecognized: list[dict] = field(default_factory=list)


# --------------------------------------------------------------------------- seqrepo IO


def open_seqrepo_db(seqrepo_root: Path) -> sqlite3.Connection:
    db = seqrepo_root / "aliases.sqlite3"
    if not db.exists():
        raise FileNotFoundError(f"aliases.sqlite3 not found under {seqrepo_root}")
    # read-only URI connection
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    return conn


def seqrepo_namespace_counts(conn: sqlite3.Connection) -> dict[str, int]:
    cur = conn.execute(
        "SELECT namespace, count(*) FROM seqalias WHERE is_current=1 GROUP BY namespace"
    )
    return {ns: n for ns, n in cur.fetchall()}


def load_seqrepo_aliases(conn: sqlite3.Connection, namespace: str) -> dict[str, str]:
    cur = conn.execute(
        "SELECT alias, seq_id FROM seqalias WHERE namespace=? AND is_current=1",
        (namespace,),
    )
    return {alias: seq_id for alias, seq_id in cur.fetchall()}


# --------------------------------------------------------------------------- refget IO


def load_refget_aliases(store_path: Path, namespace: str) -> dict[str, str]:
    """Read store/aliases/sequences/<ns>.tsv directly (alias\\tdigest)."""
    tsv = store_path / "aliases" / "sequences" / f"{namespace}.tsv"
    out: dict[str, str] = {}
    if not tsv.exists():
        # expected for "additional" namespaces the store does not (yet) expose
        logger.debug("no alias TSV for store namespace %s (%s)", namespace, tsv)
        return out
    with tsv.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            alias, _, digest = line.partition("\t")
            if digest:
                out[alias] = digest
    return out


def build_refget_digest_set(store: RefgetStore) -> set[str]:
    return {md.sha512t24u for md in store.list_sequences()}


# --------------------------------------------------------------------------- helpers


def strip_version(accession: str) -> tuple[str, int | None]:
    """Split a trailing numeric ``.N`` version. ``NM_000551.3`` -> (NM_000551, 3)."""
    base, dot, suffix = accession.rpartition(".")
    if dot and suffix.isdigit():
        return base, int(suffix)
    return accession, None


def accession_prefix(accession: str) -> str:
    """Coarse category key, e.g. ``NM_``, ``XP_``, ``ENST``. Falls back to a short stub."""
    if len(accession) >= 3 and accession[2] == "_":
        return accession[:3]
    # Ensembl-style: ENST/ENSP/ENSG...
    for n in (4, 3):
        if accession[:n].isalpha():
            return accession[:n]
    return accession[:3]


def build_version_index(rg_map: dict[str, str]) -> dict[str, set[int]]:
    idx: dict[str, set[int]] = defaultdict(set)
    for alias in rg_map:
        base, ver = strip_version(alias)
        if ver is not None:
            idx[base].add(ver)
    return idx


def load_known_divergent(path: Path) -> dict[str, str]:
    """accession -> cause. Tolerates absence (returns empty)."""
    out: dict[str, str] = {}
    if not path.exists():
        logger.warning("known-divergent file not found at %s; mismatches will be "
                       "reported as UNEXPLAINED", path)
        return out
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                out[parts[0]] = parts[3]
            elif len(parts) >= 1:
                out[parts[0]] = "listed"
    return out


# --------------------------------------------------------------------------- level A


def classify_namespaces(
    seqrepo_counts: dict[str, int], refget_namespaces: list[str]
) -> NamespaceClassification:
    nc = NamespaceClassification()
    refget_targets = set(NAMESPACE_MAP.values())
    for ns, count in sorted(seqrepo_counts.items(), key=lambda kv: -kv[1]):
        if ns in NAMESPACE_MAP:
            nc.mapped.append(
                {"seqrepo_ns": ns, "refget_ns": NAMESPACE_MAP[ns], "seqrepo_count": count}
            )
        elif ns in EXPECTED_OMIT:
            nc.expected_omit.append({"seqrepo_ns": ns, "count": count})
        elif ns in VERSION_DRIFT_NAMESPACES:
            nc.version_drift.append({"seqrepo_ns": ns, "count": count})
        else:
            nc.unrecognized.append({"seqrepo_ns": ns, "count": count})
    for ns in sorted(refget_namespaces):
        if ns not in refget_targets:
            nc.refget_only.append(ns)
    return nc


# --------------------------------------------------------------------------- level B


def compare_aliases(
    seqrepo_ns: str,
    refget_ns: str,
    group: str,
    sr_map: dict[str, str],
    rg_map: dict[str, str],
    known_divergent: dict[str, str],
    refget_digests: set[str],
    max_examples: int,
) -> AliasComparison:
    """Compare one mapped SeqRepo namespace with its RefgetStore counterpart.

    Args:
        seqrepo_ns: Source SeqRepo namespace.
        refget_ns: Corresponding RefgetStore namespace.
        group: Compatibility group, such as ``core`` or ``additional``.
        sr_map: SeqRepo alias-to-digest mapping.
        rg_map: RefgetStore alias-to-digest mapping.
        known_divergent: Accepted alias-to-cause digest divergences.
        refget_digests: All digests present in the RefgetStore.
        max_examples: Per-category cap for diagnostic examples.

    Shared aliases are classified as matching, explained, or unexplained digest
    mismatches. SeqRepo-only aliases are then attributed to version drift, an
    alias-naming gap, a backfill candidate, or a missing sequence.
    """
    cmp = AliasComparison(
        seqrepo_ns=seqrepo_ns,
        refget_ns=refget_ns,
        group=group,
        seqrepo_alias_count=len(sr_map),
        refget_alias_count=len(rg_map),
    )
    version_index = build_version_index(rg_map)
    by_cause: dict[str, int] = defaultdict(int)

    sr_keys = sr_map.keys()
    for alias in sr_keys:
        rg_digest = rg_map.get(alias)
        if rg_digest is not None:
            cmp.common += 1
            if rg_digest == sr_map[alias]:
                cmp.matched += 1
            else:
                cmp.mismatch_total += 1
                cause = known_divergent.get(alias)
                if cause:
                    cmp.mismatch_explained += 1
                    by_cause[cause] += 1
                else:
                    cmp.mismatch_unexplained += 1
                    by_cause["unexplained"] += 1
                    if len(cmp.examples_unexplained) < max_examples:
                        cmp.examples_unexplained.append(
                            {"alias": alias, "seqrepo_seq_id": sr_map[alias],
                             "refget_digest": rg_digest}
                        )
            continue

        # seqrepo-only alias: classify under latest+best-effort-old policy
        base, ver = strip_version(alias)
        versions = version_index.get(base)
        if versions:
            if ver is None or any(v > ver for v in versions):
                cmp.superseded_old_version += 1
            else:
                # store has the base but only lower (or equal-but-not-this) versions
                cmp.backfill_candidate += 1
                if len(cmp.examples_backfill) < max_examples:
                    cmp.examples_backfill.append(
                        {"alias": alias, "refget_versions": sorted(versions)}
                    )
        elif sr_map[alias] in refget_digests:
            # the sequence exists in the store, just not under this alias spelling
            cmp.alias_naming_gap += 1
            pfx = accession_prefix(alias)
            cmp.alias_naming_by_prefix[pfx] = cmp.alias_naming_by_prefix.get(pfx, 0) + 1
            if len(cmp.examples_alias_naming) < max_examples:
                cmp.examples_alias_naming.append(
                    {"alias": alias, "digest": sr_map[alias]}
                )
        else:
            cmp.sequence_missing += 1
            pfx = accession_prefix(alias)
            cmp.sequence_missing_by_prefix[pfx] = (
                cmp.sequence_missing_by_prefix.get(pfx, 0) + 1
            )
            if len(cmp.examples_sequence_missing) < max_examples:
                cmp.examples_sequence_missing.append(alias)

    cmp.refget_only = len(rg_map.keys() - sr_keys)
    cmp.mismatch_by_cause = dict(by_cause)
    return cmp


# --------------------------------------------------------------------------- level C


def compute_digest_coverage(
    conn: sqlite3.Connection, refget_digests: set[str], max_examples: int
) -> DigestCoverage:
    dc = DigestCoverage()
    # distinct seqrepo sequences (every sequence has MD5/SEGUID/... aliases, so this
    # is the full sequence inventory regardless of accession namespace)
    cur = conn.execute("SELECT DISTINCT seq_id FROM seqalias WHERE is_current=1")
    seqrepo_digests = {row[0] for row in cur}
    dc.seqrepo_distinct_digests = len(seqrepo_digests)

    missing = seqrepo_digests - refget_digests
    dc.present_in_refget = len(seqrepo_digests) - len(missing)
    dc.missing = len(missing)
    dc.coverage_pct = round(
        100.0 * dc.present_in_refget / max(1, dc.seqrepo_distinct_digests), 4
    )
    if not missing:
        return dc

    # attribute each missing digest: which accession-bearing namespaces reference it?
    mapped_ns = set(NAMESPACE_MAP.keys())
    cur = conn.execute(
        "SELECT seq_id, namespace FROM seqalias WHERE is_current=1 "
        "AND namespace NOT IN ('MD5','SEGUID','SHA1','VMC')"
    )
    missing_ns: dict[str, set[str]] = defaultdict(set)
    for seq_id, ns in cur:
        if seq_id in missing:
            missing_ns[seq_id].add(ns)

    for seq_id in missing:
        namespaces = missing_ns.get(seq_id, set())
        if namespaces & mapped_ns:
            dc.missing_referenced_by_mapped += 1
            if len(dc.examples_missing_mapped) < max_examples:
                dc.examples_missing_mapped.append(
                    {"digest": seq_id, "seqrepo_namespaces": sorted(namespaces)}
                )
        elif namespaces & VERSION_DRIFT_NAMESPACES:
            dc.missing_version_drift_only += 1
        else:
            dc.missing_other += 1
    return dc


# --------------------------------------------------------------------------- verdict


def compute_verdict(
    comparisons: list[AliasComparison], coverage: DigestCoverage, fail_on: str
) -> dict:
    """Convert namespace comparisons into the backwards-compatibility verdict.

    Args:
        comparisons: Per-namespace alias comparison results.
        coverage: Whole-store digest coverage result, retained for report context.
        fail_on: ``none``, ``mismatch``, or ``gaps`` failure policy.

    Only core namespaces determine the default verdict. Additional patch and
    legacy namespaces remain visible in the report but are not folded into the
    core pass/fail criteria.
    """
    core = [c for c in comparisons if c.group == "core"]
    extra = [c for c in comparisons if c.group == "additional"]
    unexplained = sum(c.mismatch_unexplained for c in core)
    sequence_missing = sum(c.sequence_missing for c in core)
    alias_naming = sum(c.alias_naming_gap for c in core)
    backfill = sum(c.backfill_candidate for c in comparisons)
    mismatches = sum(c.mismatch_total for c in core)

    # additional (patch-ladder / legacy-assembly) namespaces: reported, and their
    # alias-naming gaps are cheap to close (sequence already in store), so surfaced
    # separately rather than folded into the core pass/fail.
    extra_alias_naming = sum(c.alias_naming_gap for c in extra)
    extra_sequence_missing = sum(c.sequence_missing for c in extra)

    criteria = {
        "no_unexplained_digest_mismatches": unexplained == 0,
        "no_missing_sequences_in_core_namespaces": sequence_missing == 0,
        "no_alias_naming_gaps_in_core_namespaces": alias_naming == 0,
    }
    fail_reasons: list[str] = []
    if not criteria["no_unexplained_digest_mismatches"]:
        fail_reasons.append(f"{unexplained} unexplained digest mismatch(es) (core)")
    if not criteria["no_missing_sequences_in_core_namespaces"]:
        fail_reasons.append(
            f"{sequence_missing} missing sequence(s) in core namespaces (absent from build)"
        )
    if not criteria["no_alias_naming_gaps_in_core_namespaces"]:
        fail_reasons.append(
            f"{alias_naming} alias-naming gap(s) in core namespaces (sequence present, "
            "alias string not resolvable — fix by adding alias spellings)"
        )

    if fail_on == "none":
        backwards_compatible = True
    elif fail_on == "mismatch":
        backwards_compatible = not fail_reasons and mismatches == 0
        if mismatches and fail_on == "mismatch":
            fail_reasons.append(f"{mismatches} digest mismatch(es) (--fail-on mismatch)")
    else:  # "gaps" (default)
        backwards_compatible = not fail_reasons

    return {
        "backwards_compatible": backwards_compatible,
        "criteria": criteria,
        "fail_reasons": fail_reasons,
        "warnings": {"backfill_candidates": backfill},
        "additional_namespaces": {
            "alias_naming_gap": extra_alias_naming,
            "sequence_missing": extra_sequence_missing,
            "note": "patch-ladder / legacy-assembly namespaces seqrepo exposes but the "
                    "store does not; alias_naming_gap entries are cheap to close "
                    "(sequence already present), sequence_missing need new sources.",
        },
    }


# --------------------------------------------------------------------------- output


def write_gap_list(path: Path, comparisons: list[AliasComparison],
                   sr_maps: dict[str, dict[str, str]],
                   rg_maps: dict[str, dict[str, str]],
                   refget_digests: set[str]) -> None:
    """Machine-readable feed for Workstream 2: every seqrepo-only accession with its
    remediation category and prefix.

    Columns: seqrepo_ns<TAB>refget_ns<TAB>accession<TAB>category<TAB>prefix<TAB>seqrepo_seq_id
    Categories: sequence_missing (load a source), alias_naming_gap (add alias spelling),
    backfill_candidate (store behind), superseded_old_version (optional old-version backfill).
    """
    with path.open("w", encoding="utf-8") as out:
        out.write("seqrepo_ns\trefget_ns\taccession\tcategory\tprefix\tseqrepo_seq_id\n")
        for cmp in comparisons:
            sr = sr_maps[cmp.seqrepo_ns]
            rg = rg_maps[cmp.seqrepo_ns]
            version_index = build_version_index(rg)
            for alias, seq_id in sr.items():
                if alias in rg:
                    continue  # in store (matched or mismatch handled elsewhere)
                base, ver = strip_version(alias)
                versions = version_index.get(base)
                if versions:
                    if ver is None or any(v > ver for v in versions):
                        category = "superseded_old_version"
                    else:
                        category = "backfill_candidate"
                elif seq_id in refget_digests:
                    category = "alias_naming_gap"
                else:
                    category = "sequence_missing"
                out.write(
                    f"{cmp.seqrepo_ns}\t{cmp.refget_ns}\t{alias}\t{category}\t"
                    f"{accession_prefix(alias)}\t{seq_id}\n"
                )


def render_summary(report: dict) -> None:
    print()
    print("=" * 78)
    print("seqrepo <-> refgetstore equivalence")
    print("=" * 78)
    m = report["meta"]
    print(f"  seqrepo : {m['seqrepo_path']} (snapshot {m['seqrepo_snapshot']})")
    print(f"  store   : {m['store_path']}")
    print(f"  gtars   : {m['gtars_version']}")

    nc = report["namespace_classification"]
    print("\n[A] Namespace classification")
    print(f"    mapped        : {[d['seqrepo_ns'] + '->' + d['refget_ns'] for d in nc['mapped']]}")
    print(f"    expected-omit : {[d['seqrepo_ns'] for d in nc['expected_omit']]}")
    print(f"    version-drift : {[d['seqrepo_ns'] for d in nc['version_drift']]}")
    print(f"    refget-only   : {nc['refget_only']}")
    if nc["unrecognized"]:
        print(f"    UNRECOGNIZED  : {nc['unrecognized']}")

    def _row(c: dict) -> None:
        print(
            f"  {c['seqrepo_ns']:>10s} -> {c['refget_ns']:<10s} "
            f"sr={c['seqrepo_alias_count']:>7d} rg={c['refget_alias_count']:>7d} "
            f"matched={c['matched']:>7d} "
            f"mism={c['mismatch_total']:>4d}(u={c['mismatch_unexplained']}) "
            f"superseded={c['superseded_old_version']:>6d} "
            f"naming={c['alias_naming_gap']:>5d} "
            f"MISSING={c['sequence_missing']:>6d} ext={c['refget_only']:>6d} "
            f"cov={c['backwards_compat_coverage_pct']:.2f}%"
        )

    core = [c for c in report["alias_comparison"] if c.get("group") == "core"]
    extra = [c for c in report["alias_comparison"] if c.get("group") == "additional"]

    print("\n[B] Core namespace alias comparison")
    for c in core:
        _row(c)
        if c["sequence_missing"]:
            tops = sorted(c["sequence_missing_by_prefix"].items(), key=lambda kv: -kv[1])[:8]
            print(f"             missing-sequence by prefix: {tops}")
        if c["alias_naming_gap"]:
            tops = sorted(c["alias_naming_by_prefix"].items(), key=lambda kv: -kv[1])[:8]
            print(f"             alias-naming by prefix: {tops}")

    if extra:
        en = sum(c["alias_naming_gap"] for c in extra)
        em = sum(c["sequence_missing"] for c in extra)
        print("\n[B2] Additional seqrepo namespaces (patch ladders / legacy assemblies,"
              " not exposed by the store)")
        print(f"     totals: alias_naming_gap={en} (sequence present, cheap to add), "
              f"sequence_missing={em} (need new sources)")
        for c in extra:
            print(f"     {c['seqrepo_ns']:>12s}  sr={c['seqrepo_alias_count']:>5d} "
                  f"naming={c['alias_naming_gap']:>5d} MISSING={c['sequence_missing']:>4d}")

    dc = report["digest_coverage"]
    print("\n[C] Sequence/digest coverage (all seqrepo sequences)")
    print(f"    {dc['present_in_refget']}/{dc['seqrepo_distinct_digests']} present "
          f"({dc['coverage_pct']:.2f}%); missing={dc['missing']}")
    print(f"    missing referenced by MAPPED current alias : {dc['missing_referenced_by_mapped']}")
    print(f"    missing version-drift-only                 : {dc['missing_version_drift_only']}")
    print(f"    missing other                              : {dc['missing_other']}")

    v = report["verdict"]
    print("\n" + "-" * 78)
    an = v.get("additional_namespaces", {})
    if an.get("alias_naming_gap") or an.get("sequence_missing"):
        print(f"  additional namespaces: {an.get('alias_naming_gap', 0)} alias-naming gap(s) "
              f"(cheap to add), {an.get('sequence_missing', 0)} missing sequence(s)")
    if v["warnings"]["backfill_candidates"]:
        print(f"  warning: {v['warnings']['backfill_candidates']} backfill-candidate "
              "accession(s) (store behind seqrepo; best-effort old-version backfill)")
    if v["fail_reasons"]:
        for reason in v["fail_reasons"]:
            print(f"  FAIL: {reason}")
    print(f"\nBACKWARDS-COMPATIBLE (core namespaces): "
          f"{'YES' if v['backwards_compatible'] else 'NO'}")


# --------------------------------------------------------------------------- deep mode


def deep_diagnose(comparisons: list[AliasComparison], seqrepo_root: Path) -> None:
    try:
        from biocommons.seqrepo import SeqRepo  # noqa: F401
    except Exception:
        print(
            "ERROR: --deep requires biocommons.seqrepo (not importable here).\n"
            "       Run under an env that has it (e.g. ~/dev/vrs-python/.venv) or use\n"
            "       vrs-python's misc/refgetstore/diagnose_ensembl_mismatches.py.",
            file=sys.stderr,
        )
        raise SystemExit(2)
    # Byte-level diagnosis of only the unexplained mismatches would go here. The
    # validated implementation lives in diagnose_ensembl_mismatches.py; this hook is
    # intentionally a thin bridge so the gtars-only core stays dependency-light.
    total = sum(len(c.examples_unexplained) for c in comparisons)
    logger.info("deep mode: %d unexplained mismatch example(s) available to diagnose", total)


# --------------------------------------------------------------------------- main


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seqrepo", type=Path,
                        default=Path(os.environ.get("SEQREPO_ROOT_DIR", str(DEFAULT_SEQREPO))))
    parser.add_argument("--store", type=Path, default=DEFAULT_STORE)
    parser.add_argument("--known-divergent", type=Path, default=DEFAULT_KNOWN_DIVERGENT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--gap-list", type=Path, default=DEFAULT_GAP_LIST)
    parser.add_argument("--namespaces", type=str, default=None,
                        help="Comma-separated subset of seqrepo namespaces to compare.")
    parser.add_argument("--max-examples", type=int, default=50)
    parser.add_argument("--fail-on", choices=("none", "gaps", "mismatch"), default="gaps")
    parser.add_argument("--deep", action="store_true",
                        help="Byte-level diagnosis of unexplained mismatches (needs seqrepo).")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if not args.seqrepo.exists():
        print(f"ERROR: seqrepo snapshot not found at {args.seqrepo}", file=sys.stderr)
        return 2
    if not (args.store / "rgstore.json").exists():
        print(f"ERROR: no RefgetStore at {args.store}", file=sys.stderr)
        return 2

    start = time.monotonic()
    try:
        import gtars
        gtars_version = getattr(gtars, "__version__", None)
    except Exception:
        gtars_version = None

    conn = open_seqrepo_db(args.seqrepo)
    store = RefgetStore.open_local(str(args.store))
    store.pull_aliases()

    known_divergent = load_known_divergent(args.known_divergent)
    seqrepo_counts = seqrepo_namespace_counts(conn)
    refget_namespaces = list(store.list_sequence_alias_namespaces())
    nc = classify_namespaces(seqrepo_counts, refget_namespaces)

    logger.info("building store digest set")
    refget_digests = build_refget_digest_set(store)

    # Thorough alias equivalence: compare EVERY accession-bearing seqrepo namespace,
    # not just the explicitly-mapped ones. Namespaces in NAMESPACE_MAP are "core"
    # (their refget target may be spelled differently, e.g. NCBI->refseq); every other
    # non-digest namespace is "additional" and compared against an identically-named
    # store namespace (which usually doesn't exist yet — patch ladders, legacy builds).
    selected = set(args.namespaces.split(",")) if args.namespaces else None
    core_keys = list(NAMESPACE_MAP.keys())
    additional = sorted(
        (ns for ns in seqrepo_counts
         if ns not in EXPECTED_OMIT and ns not in NAMESPACE_MAP),
        key=lambda ns: -seqrepo_counts[ns],
    )
    plan: list[tuple[str, str, str]] = (
        [(ns, NAMESPACE_MAP[ns], "core") for ns in core_keys]
        + [(ns, ns, "additional") for ns in additional]
    )

    comparisons: list[AliasComparison] = []
    sr_maps: dict[str, dict[str, str]] = {}
    rg_maps: dict[str, dict[str, str]] = {}
    for seqrepo_ns, refget_ns, group in plan:
        if selected and seqrepo_ns not in selected:
            continue
        logger.info("comparing %s -> %s (%s)", seqrepo_ns, refget_ns, group)
        sr_map = load_seqrepo_aliases(conn, seqrepo_ns)
        rg_map = load_refget_aliases(args.store, refget_ns)
        sr_maps[seqrepo_ns] = sr_map
        rg_maps[seqrepo_ns] = rg_map
        comparisons.append(
            compare_aliases(seqrepo_ns, refget_ns, group, sr_map, rg_map,
                            known_divergent, refget_digests, args.max_examples)
        )

    logger.info("computing digest coverage")
    coverage = compute_digest_coverage(conn, refget_digests, args.max_examples)

    verdict = compute_verdict(comparisons, coverage, args.fail_on)

    report = {
        "meta": {
            "seqrepo_path": str(args.seqrepo),
            "seqrepo_snapshot": args.seqrepo.name,
            "store_path": str(args.store),
            "gtars_version": gtars_version,
            "known_divergent_file": str(args.known_divergent)
            if args.known_divergent.exists() else None,
            "fail_on": args.fail_on,
            "tool_version": "1.0",
        },
        "verdict": verdict,
        "namespace_classification": asdict(nc),
        "alias_comparison": [
            {**asdict(c), "backwards_compat_coverage_pct": c.backwards_compat_coverage_pct}
            for c in comparisons
        ],
        "digest_coverage": asdict(coverage),
        "metrics": {
            "elapsed_seconds": round(time.monotonic() - start, 2),
        },
    }

    args.report.write_text(json.dumps(report, indent=2) + "\n")
    write_gap_list(args.gap_list, comparisons, sr_maps, rg_maps, refget_digests)
    logger.info("wrote %s and %s", args.report, args.gap_list)

    if args.deep:
        deep_diagnose(comparisons, args.seqrepo)

    render_summary(report)
    print(f"\nReport: {args.report}\nGap list: {args.gap_list}")
    print(f"Elapsed: {report['metrics']['elapsed_seconds']}s")

    return 0 if verdict["backwards_compatible"] else 1


if __name__ == "__main__":
    sys.exit(main())
