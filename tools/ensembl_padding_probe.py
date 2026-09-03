#!/usr/bin/env python3
"""Characterize Ensembl's former N-padded alt/patch scaffold representation.

Emits two tab-separated tables into --out-dir:

  ensembl_padded_scaffolds.tsv            one row per distinct CHR_-prefixed
                                          padded record found in ensembl-75..109
  ensembl_dna_representation_by_release.tsv
                                          one row per Ensembl release 75..116

Inputs are the local RefgetStore, a seqrepo aliases.sqlite3 (read-only), the
cached NCBI GRCh38 assembly reports under downloads/, and NCBI
alt_scaffold_placement.txt files (see --placements; fetch with the curl loop
documented in seqrepo_equivalence/ENSEMBL_N_PADDED_SCAFFOLDS.md).

Usage:
    uv run python tools/ensembl_padding_probe.py \
        --out-dir seqrepo_equivalence \
        --placements /tmp/placements
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
import time
from pathlib import Path

RELEASES = list(range(75, 117))
CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
DIGEST_NAMESPACES = {"MD5", "SEGUID", "SHA1", "VMC"}
FTP_BASE = "https://ftp.ensembl.org/pub/release-{rel}/fasta/homo_sapiens/dna/"

REPORT_GLOB = (
    "downloads/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
    "GCF_000001405.*_GRCh38*/GCF_000001405.*_GRCh38*_assembly_report.txt"
)


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------
# NCBI inputs
# --------------------------------------------------------------------------


def load_assembly_reports(repo_root: Path) -> dict[str, dict[str, str]]:
    """Sequence-Name -> role/accessions, unioned over every GRCh38 report.

    Reports are visited in ascending patch order so that the newest statement
    about a scaffold wins; scaffolds retired before p14 survive from the older
    report that last listed them.
    """
    paths = sorted(
        repo_root.glob(REPORT_GLOB),
        key=lambda p: int(re.search(r"GCF_000001405\.(\d+)_", p.name).group(1)),
    )
    if not paths:
        raise SystemExit(f"no GRCh38 assembly reports under {repo_root}")
    out: dict[str, dict[str, str]] = {}
    for path in paths:
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 10:
                continue
            out[f[0]] = {
                "sequence_role": f[1],
                "assigned_molecule": f[2],
                "genbank_accn": f[4],
                "refseq_accn": f[6],
                "assembly_unit": f[7],
                "report_length": f[8],
                "ucsc_name": f[9],
                "last_report": path.parent.name,
            }
    log(f"assembly reports: {len(paths)} files, {len(out)} distinct Sequence-Name")
    return out


def load_placements(placements_dir: Path) -> dict[str, dict[str, str]]:
    """alt_scaf_name -> parent placement, from alt_scaffold_placement.txt files."""
    out: dict[str, dict[str, str]] = {}
    files = sorted(placements_dir.glob("*.txt"))
    for path in files:
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 15:
                continue
            out[f[2]] = {
                "alt_scaf_acc": f[3],
                "parent_chromosome": f[5],
                "parent_acc": f[6],
                "region_name": f[7],
                "orientation": f[8],
                "alt_scaf_start": f[9],
                "alt_scaf_stop": f[10],
                "parent_start": f[11],
                "parent_stop": f[12],
                "alt_start_tail": f[13],
                "alt_stop_tail": f[14],
            }
    log(f"placements: {len(files)} files, {len(out)} scaffolds")
    return out


# --------------------------------------------------------------------------
# seqrepo
# --------------------------------------------------------------------------


def seqrepo_namespaces(db: Path, digests: list[str]) -> dict[str, list[str]]:
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    out: dict[str, list[str]] = {}
    for i in range(0, len(digests), 400):
        chunk = digests[i : i + 400]
        q = "select seq_id, namespace from seqalias where seq_id in (%s)" % (
            ",".join("?" * len(chunk))
        )
        for seq_id, ns in conn.execute(q, chunk):
            if ns in DIGEST_NAMESPACES:
                out.setdefault(seq_id, [])
            else:
                out.setdefault(seq_id, []).append(ns)
    conn.close()
    return {k: sorted(set(v)) for k, v in out.items()}


# --------------------------------------------------------------------------
# store scan
# --------------------------------------------------------------------------


def scan_releases(store):
    """Per release: the non-ENS (DNA) alias partition and its metadata."""
    per_release = {}
    for rel in RELEASES:
        ns = f"ensembl-{rel}"
        aliases = store.list_sequence_aliases(ns)
        dna = [a for a in aliases if not a.startswith("ENS")]
        chr_pref = sorted(a for a in dna if a.startswith("CHR_"))
        chroms = sorted(a for a in dna if a in CHROMOSOMES)
        bare = sorted(a for a in dna if a not in CHROMOSOMES and not a.startswith("CHR_"))
        meta = {}
        for a in dna:
            m = store.get_sequence_metadata_by_alias(ns, a)
            if m is not None:
                meta[a] = (m.sha512t24u, m.length, str(m.alphabet))
        per_release[rel] = {
            "chr_prefixed": chr_pref,
            "bare": bare,
            "chromosomes": chroms,
            "meta": meta,
        }
        log(
            f"ensembl-{rel}: CHR_={len(chr_pref)} bare={len(bare)} "
            f"chrom={len(chroms)} total_aliases={len(aliases)}"
        )
    return per_release


def sha512t24u(seq: str) -> str:
    import base64
    import hashlib

    return base64.urlsafe_b64encode(hashlib.sha512(seq.encode()).digest()[:24]).decode()


COMPLEMENT = str.maketrans(
    "ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb"
)


def measure(store, digest: str, length: int, unpadded_digest: str):
    """(real_bases, span_start, span_end, recovery); N/n is treated as padding.

    recovery says how the padded record collapses to its true-length form:
    `forward` if stripping the leading and trailing N runs reproduces
    unpadded_digest, `reverse_complement` if reverse-complementing that span
    does, otherwise `no`.
    """
    seq = store.get_substring(digest, 0, length)
    real = length - seq.count("N") - seq.count("n")
    if real == 0:
        return 0, None, None, "no"
    start = length - len(seq.lstrip("Nn"))
    end = len(seq.rstrip("Nn"))
    span = seq[start:end]
    if sha512t24u(span) == unpadded_digest:
        rec = "forward"
    elif sha512t24u(span.translate(COMPLEMENT)[::-1]) == unpadded_digest:
        rec = "reverse_complement"
    else:
        rec = "no"
    return real, start, end, rec


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

PADDED_COLUMNS = [
    "padded_name",
    "bare_name",
    "padded_digest",
    "padded_length",
    "unpadded_digest",
    "unpadded_length",
    "inflation_factor",
    "sequence_role",
    "refseq_accn",
    "genbank_accn",
    "first_padded_release",
    "last_padded_release",
    "first_unpadded_release",
    "unpadded_source",
    "parent_chromosome",
    "parent_start",
    "parent_stop",
    "region_name",
    "orientation",
    "real_bases",
    "real_span_start",
    "real_span_end",
    "pct_real",
    "recovers_unpadded",
    "placement_matches",
    "in_seqrepo_padded",
    "seqrepo_biological_namespaces_padded",
]

RELEASE_COLUMNS = [
    "release",
    "assembly",
    "chr_prefixed_count",
    "bare_alt_patch_count",
    "chromosome_count",
    "dna3bit_bases",
    "ftp_last_modified",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--store", default="store")
    ap.add_argument(
        "--seqrepo", default=str(Path.home() / "dev/data/seqrepo/2024-12-20/aliases.sqlite3")
    )
    ap.add_argument("--placements", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--repo-root", default=".", type=Path)
    ap.add_argument(
        "--measure",
        default="all",
        help="'all', 'none', or an integer cap on how many padded records to decode",
    )
    ap.add_argument(
        "--last-modified",
        action="store_true",
        help="probe Ensembl FTP for each release's dna.toplevel Last-Modified header",
    )
    args = ap.parse_args()

    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(args.store)
    reports = load_assembly_reports(args.repo_root)
    placements = load_placements(args.placements)
    per_release = scan_releases(store)

    # ---- padded records --------------------------------------------------
    padded_names: dict[str, dict] = {}
    for rel in RELEASES:
        for name in per_release[rel]["chr_prefixed"]:
            digest, length, _alpha = per_release[rel]["meta"][name]
            rec = padded_names.setdefault(
                name,
                {"digests": {}, "first": rel, "last": rel},
            )
            rec["digests"][digest] = length
            rec["first"] = min(rec["first"], rel)
            rec["last"] = max(rec["last"], rel)
    log(f"distinct CHR_-prefixed names: {len(padded_names)}")

    # ---- unpadded counterparts ------------------------------------------
    # Search the GRCh38 true-length era (110+) first, then the GRCh38 padded era
    # (a bare name there is an unplaced/unlocalized scaffold), and only then
    # release 75, which is GRCh37 and may hold a same-named different sequence.
    bare_index: dict[str, tuple[str, int, int]] = {}  # bare -> (digest, length, first rel)
    for rel in list(range(110, 117)) + list(range(76, 110)) + [75]:
        for name in per_release[rel]["bare"]:
            if name not in bare_index:
                digest, length, _ = per_release[rel]["meta"][name]
                bare_index[name] = (digest, length, rel)

    rows = []
    measure_cap = None
    if args.measure == "none":
        measure_cap = 0
    elif args.measure != "all":
        measure_cap = int(args.measure)

    ordered = sorted(padded_names.items(), key=lambda kv: kv[0][4:])
    t0 = time.time()
    measured = 0
    for name, rec in ordered:
        bare = name[4:]
        digest = max(rec["digests"], key=lambda d: rec["digests"][d])
        length = rec["digests"][digest]
        rep = reports.get(bare, {})
        plc = placements.get(bare, {})

        unp_digest = unp_len = ""
        first_unpadded = ""
        unp_source = ""
        if bare in bare_index:
            unp_digest, unp_len, first_unpadded = bare_index[bare]
            unp_source = f"ensembl-{first_unpadded}"
        else:
            for ns, accn in (("refseq", rep.get("refseq_accn")), ("insdc", rep.get("genbank_accn"))):
                if not accn or accn in ("na", ""):
                    continue
                m = store.get_sequence_metadata_by_alias(ns, accn)
                if m is not None:
                    unp_digest, unp_len = m.sha512t24u, m.length
                    unp_source = f"{ns}:{accn}"
                    break

        row = {
            "padded_name": name,
            "bare_name": bare,
            "padded_digest": digest,
            "padded_length": length,
            "unpadded_digest": unp_digest,
            "unpadded_length": unp_len,
            "inflation_factor": f"{length / unp_len:.1f}" if unp_len else "",
            "sequence_role": rep.get("sequence_role", ""),
            "refseq_accn": rep.get("refseq_accn", ""),
            "genbank_accn": rep.get("genbank_accn", ""),
            "first_padded_release": rec["first"],
            "last_padded_release": rec["last"],
            "first_unpadded_release": first_unpadded,
            "unpadded_source": unp_source,
            "parent_chromosome": plc.get("parent_chromosome", ""),
            "parent_start": plc.get("parent_start", ""),
            "parent_stop": plc.get("parent_stop", ""),
            "region_name": plc.get("region_name", ""),
            "orientation": plc.get("orientation", ""),
            "real_bases": "",
            "real_span_start": "",
            "real_span_end": "",
            "pct_real": "",
            "recovers_unpadded": "",
            "placement_matches": "",
        }

        if measure_cap is None or measured < measure_cap:
            real, st, en, recovery = measure(store, digest, length, unp_digest)
            measured += 1
            row["real_bases"] = real
            row["real_span_start"] = "" if st is None else st
            row["real_span_end"] = "" if en is None else en
            row["pct_real"] = f"{100 * real / length:.4f}"
            row["recovers_unpadded"] = recovery
            if st is not None and row["parent_start"]:
                row["placement_matches"] = (
                    "yes" if st == int(row["parent_start"]) - 1 else "no"
                )
            if measured % 25 == 0:
                log(f"  measured {measured}/{len(ordered)} ({time.time() - t0:.0f}s)")
        rows.append(row)

    # ---- consistency checks ---------------------------------------------
    bad_len = [
        r["bare_name"]
        for r in rows
        if reports.get(r["bare_name"], {}).get("report_length")
        and str(r["unpadded_length"]) != reports[r["bare_name"]]["report_length"]
    ]
    log(
        f"unpadded length vs NCBI assembly report: {len(rows) - len(bad_len)}/{len(rows)} agree"
        + (f"; disagree: {bad_len[:10]}" if bad_len else "")
    )
    log(f"no unpadded counterpart in store: {sum(1 for r in rows if not r['unpadded_digest'])}")
    outside = 0
    for r in rows:
        nss = set()
        for a in store.get_aliases_for_sequence(r["padded_digest"]):
            nss.add(a[0] if isinstance(a, (tuple, list)) else getattr(a, "namespace", str(a)))
        if any(not re.fullmatch(r"ensembl-\d+", ns) for ns in nss):
            outside += 1
    log(f"padded digests carrying any alias outside ensembl-N: {outside}")
    sd = [r["recovers_unpadded"] for r in rows if r["recovers_unpadded"]]
    if sd:
        log(
            f"recovers_unpadded: forward={sd.count('forward')} "
            f"reverse_complement={sd.count('reverse_complement')} no={sd.count('no')} "
            f"(of {len(sd)} measured rows)"
        )
    matches = [r["placement_matches"] for r in rows if r["placement_matches"]]
    if matches:
        log(
            f"placement_matches: yes={matches.count('yes')} no={matches.count('no')} "
            f"(of {len(matches)} comparable rows)"
        )

    # ---- seqrepo ---------------------------------------------------------
    sr = seqrepo_namespaces(Path(args.seqrepo), [r["padded_digest"] for r in rows])
    for r in rows:
        hit = sr.get(r["padded_digest"])
        r["in_seqrepo_padded"] = "yes" if hit is not None else "no"
        r["seqrepo_biological_namespaces_padded"] = ";".join(hit or [])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    p1 = args.out_dir / "ensembl_padded_scaffolds.tsv"
    with p1.open("w", newline="") as fh:
        w = csv.DictWriter(fh, PADDED_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {p1} ({len(rows)} rows)")

    # ---- per-release table ----------------------------------------------
    last_mod = {}
    if args.last_modified:
        import urllib.request

        for rel in RELEASES:
            url = FTP_BASE.format(rel=rel)
            try:
                with urllib.request.urlopen(url, timeout=60) as resp:
                    listing = resp.read().decode("utf-8", "replace")
            except Exception as exc:  # pragma: no cover
                log(f"  rel {rel}: listing failed: {exc}")
                continue
            m = re.findall(r'href="([^"]*dna\.toplevel\.fa\.gz)"', listing)
            if not m:
                log(f"  rel {rel}: no dna.toplevel in listing")
                continue
            req = urllib.request.Request(url + m[0], method="HEAD")
            try:
                with urllib.request.urlopen(req, timeout=60) as resp:
                    last_mod[rel] = resp.headers.get("Last-Modified", "")
                log(f"  rel {rel}: {last_mod[rel]}")
            except Exception as exc:  # pragma: no cover
                log(f"  rel {rel}: HEAD failed: {exc}")

    rel_rows = []
    for rel in RELEASES:
        pr = per_release[rel]
        dna3bit = sum(
            length
            for (_d, length, alpha) in pr["meta"].values()
            if alpha == "dna3bit"
        )
        stamp = last_mod.get(rel, "")
        if stamp:
            from email.utils import parsedate_to_datetime

            stamp = parsedate_to_datetime(stamp).date().isoformat()
        rel_rows.append(
            {
                "release": rel,
                "assembly": "GRCh37" if rel <= 75 else "GRCh38",
                "chr_prefixed_count": len(pr["chr_prefixed"]),
                "bare_alt_patch_count": len(pr["bare"]),
                "chromosome_count": len(pr["chromosomes"]),
                "dna3bit_bases": dna3bit,
                "ftp_last_modified": stamp,
            }
        )
    p2 = args.out_dir / "ensembl_dna_representation_by_release.tsv"
    with p2.open("w", newline="") as fh:
        w = csv.DictWriter(fh, RELEASE_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rel_rows)
    log(f"wrote {p2} ({len(rel_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
