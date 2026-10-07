#!/usr/bin/env python
"""Re-import previously affected sequences into a small new store and check them.

A full rebuild takes hours. This rebuilds only what matters for the round-trip
cases in README.md: every sequence listed in the known-divergence file, plus one
control sequence per alphabet that round-trips today, so a change that fixed the
listed sequences but broke ordinary ones would show up too.

Requires the gks-refgetstore-builder repository
(https://github.com/theferrit32/gks-refgetstore-builder), a built store, its
build lock, and the download cache it was built from.

Two steps, so the same extracted input can be imported with different gtars
versions:

  extract  Reads the existing store, lock and download cache. For each chosen
           sequence it finds a collection holding it, looks up the record name
           and source file, and copies that FASTA record verbatim. Each copy is
           checked against the stored digest, so the input is byte-for-byte
           what was originally imported. Writes FASTA files and manifest.tsv.

  build    Imports those FASTA files into a new on-disk store with whichever
           gtars is installed, writes it, reopens it, and checks that every
           sequence is present and that its returned bytes hash to its digest.

    uv run python issues/gtars-encoder-alphabet/rebuild_subset.py extract --out WORK
    uv run python issues/gtars-encoder-alphabet/rebuild_subset.py build --work WORK

``extract`` only reads the existing store. ``build`` never opens it, and
refuses to write into an existing directory.
"""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from importlib.metadata import version
from pathlib import Path

from gks_refgetstore import store_census
from gks_refgetstore.verify import DEFAULT_KNOWN_BAD, load_known_bad

CONTROL_ALPHABETS = ("dna2bit", "dna3bit", "dnaio", "protein")
# Prefer controls at least this long, so a control is an ordinary record rather
# than a 3-residue fragment; shorter ones are used only if nothing else exists.
CONTROL_MIN_LENGTH = 100
MANIFEST_COLUMNS = ("digest", "name", "role", "alphabet", "length",
                    "cache_path", "fasta")


def sha512t24u(data: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(data).digest()[:24]).decode("ascii")


def read_records(path: Path, wanted: set[str]) -> dict[str, list[str]]:
    """Raw lines (header included) of every record in ``path`` named in ``wanted``."""
    found: dict[str, list[str]] = {}
    current: list[str] | None = None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii", errors="strict") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                current = [line] if name in wanted and name not in found else None
                if current is not None:
                    found[name] = current
            elif current is not None:
                current.append(line)
    return found


def record_digest(lines: list[str]) -> str:
    """Digest of a record's residues as gtars imports them (uppercased)."""
    residues = "".join(line.strip() for line in lines[1:]).upper()
    return sha512t24u(residues.encode("ascii"))


def cmd_extract(args) -> int:
    from gtars.refget import RefgetStore

    out = args.out
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} exists and is not empty; choose a new directory")
    (out / "fasta").mkdir(parents=True)

    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    file_bytes = {f["cache_path"]: f.get("bytes") or 0 for f in lock["inputs"]["files"]}
    sources = {c["digest"]: c["from"] for c in lock["outputs"]["collections"]}

    store = RefgetStore.open_local(str(args.store))
    store.set_quiet(True)
    metadata = store_census.sequence_metadata_by_digest(store)
    known_bad = load_known_bad(args.known_bad)
    targets = {d for d in known_bad if d in metadata}
    print(f"{len(known_bad)} listed sequence(s); {len(targets)} present in {args.store}")

    # Every place each sequence occurs: (record name, source file).
    started = time.monotonic()
    occurrences: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for collection, paths in sources.items():
        for name, digest in store_census.collection_named_members(store, collection):
            for path in paths:
                occurrences[digest].append((name, path))
    print(f"indexed {len(sources)} collection(s) in {time.monotonic() - started:.0f}s")

    def cheapest(digest: str) -> list[tuple[str, str]]:
        # Smallest source file first, so extraction scans as little as possible.
        return sorted(set(occurrences[digest]), key=lambda o: (file_bytes.get(o[1], 1 << 62), o))

    chosen: dict[str, tuple[str, str, str]] = {}   # digest -> (name, path, role)
    for digest in sorted(targets):
        if not occurrences[digest]:
            sys.exit(f"{digest} is in the store but in no locked collection")
        name, path = cheapest(digest)[0]
        chosen[digest] = (name, path, "listed")

    # One control per alphabet, preferring files already being read for the
    # listed sequences, then the smallest other source file.
    scanned = {path for _, path, _ in chosen.values()}
    for alphabet in CONTROL_ALPHABETS:
        pool = [d for d, m in metadata.items()
                if str(m.alphabet) == alphabet and d not in known_bad and occurrences[d]]
        long_enough = [d for d in pool if int(metadata[d].length) >= CONTROL_MIN_LENGTH]
        pool = long_enough or pool
        if not pool:
            print(f"  no control available for {alphabet}")
            continue

        def rank(d: str) -> tuple:
            _, path = cheapest(d)[0]
            in_scanned = any(p in scanned for _, p in occurrences[d])
            return (not in_scanned, file_bytes.get(path, 1 << 62), int(metadata[d].length), d)

        digest = min(pool, key=rank)
        reuse = [o for o in cheapest(digest) if o[1] in scanned]
        name, path = (reuse or cheapest(digest))[0]
        chosen[digest] = (name, path, "control")
        scanned.add(path)

    by_path: dict[str, dict[str, str]] = defaultdict(dict)   # path -> name -> digest
    for digest, (name, path, _) in chosen.items():
        by_path[path][name] = digest

    rows = []
    mismatched = []
    for index, (path, names) in enumerate(sorted(by_path.items())):
        local = args.cache_dir / path
        print(f"  reading {path} ({file_bytes.get(path, 0) / 1e6:,.0f} MB) for {len(names)} record(s)")
        records = read_records(local, set(names))
        fasta = out / "fasta" / f"{index:03d}-{Path(path).name.removesuffix('.gz')}"
        if not fasta.suffix:
            fasta = fasta.with_suffix(".fa")
        with fasta.open("w", encoding="ascii") as handle:
            for name, digest in sorted(names.items()):
                lines = records.get(name)
                if lines is None or record_digest(lines) != digest:
                    mismatched.append((digest, name, path, "missing" if lines is None else "digest differs"))
                    continue
                handle.writelines(lines)
                meta = metadata[digest]
                rows.append({"digest": digest, "name": name, "role": chosen[digest][2],
                             "alphabet": str(meta.alphabet), "length": str(meta.length),
                             "cache_path": path, "fasta": fasta.name})

    with (out / "manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["role"], r["alphabet"], r["name"])))

    roles = Counter((r["role"], r["alphabet"]) for r in rows)
    print(f"\nwrote {len(rows)} record(s) in {len(by_path)} FASTA file(s) to {out}")
    for (role, alphabet), n in sorted(roles.items()):
        print(f"  {role:8} {alphabet:8} {n}")
    for row in rows:
        if row["role"] == "control":
            print(f"  control {row['alphabet']:8} {row['name']} ({row['length']}) from {row['cache_path']}")
    if mismatched:
        print(f"\n{len(mismatched)} record(s) could not be extracted exactly:")
        for item in mismatched:
            print("  ", *item)
        return 1
    return 0


def cmd_build(args) -> int:
    from gtars.refget import RefgetStore

    gtars_version = version("gtars")
    work = args.work
    rows = list(csv.DictReader((work / "manifest.tsv").open(encoding="utf-8"), delimiter="\t"))
    target = args.out or work / f"store-gtars-{gtars_version}"
    if target.exists():
        sys.exit(f"{target} already exists; builds always go into a new directory")

    print(f"gtars {gtars_version}: importing {len(rows)} record(s) into {target}")
    store = RefgetStore.on_disk(str(target))
    store.set_quiet(True)
    import_errors = {}
    for fasta in sorted({r["fasta"] for r in rows}):
        try:
            store.add_sequence_collection_from_fasta(str(work / "fasta" / fasta))
        except Exception as exc:  # noqa: BLE001 - recorded per file and reported
            import_errors[fasta] = f"{type(exc).__name__}: {exc}"
    store.write()
    del store

    reopened = RefgetStore.open_local(str(target))
    reopened.set_quiet(True)
    metadata = store_census.sequence_metadata_by_digest(reopened)

    results = []
    for row in rows:
        digest = row["digest"]
        meta = metadata.get(digest)
        if meta is None:
            outcome, new_alphabet, returned = "absent", "", ""
        else:
            new_alphabet = str(meta.alphabet)
            returned = store_census.redigest_sequence(reopened, digest)
            outcome = "round-trips" if returned == digest else "differs"
        results.append({**row, "alphabet_now": new_alphabet,
                        "returned_digest": returned, "result": outcome})

    gtars_check = ""
    if hasattr(reopened, "verify"):
        report = reopened.verify(jobs=1)
        gtars_check = f"gtars verify(): checked={report.n_checked} ok={report.n_ok} failed={report.n_failed}"

    out_tsv = work / f"results-gtars-{gtars_version}.tsv"
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    tally = Counter((r["role"], r["alphabet"], r["result"]) for r in results)
    print(f"\n  {'role':8} {'alphabet':8} {'result':12} sequences")
    for (role, alphabet, outcome), n in sorted(tally.items()):
        print(f"  {role:8} {alphabet:8} {outcome:12} {n}")
    moved = Counter((r["alphabet"], r["alphabet_now"]) for r in results
                    if r["alphabet_now"] and r["alphabet_now"] != r["alphabet"])
    for (before, after), n in sorted(moved.items()):
        print(f"  alphabet selected changed {before} -> {after}: {n}")
    if gtars_check:
        print(f"  {gtars_check}")
    for fasta, error in sorted(import_errors.items()):
        print(f"  import error in {fasta}: {error}")
    print(f"wrote {out_tsv}")

    not_ok = [r for r in results if r["result"] != "round-trips"]
    if not_ok or import_errors:
        print(f"\n{len(not_ok)} sequence(s) do not round-trip; {len(import_errors)} import error(s)")
        return 1
    print(f"\nall {len(results)} sequence(s) round-trip")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    extract = sub.add_parser("extract", help="copy source records into a work directory")
    extract.add_argument("--store", type=Path, default=Path("store"))
    extract.add_argument("--lock", type=Path, default=Path("build.lock.json"))
    extract.add_argument("--cache-dir", type=Path, default=Path("downloads"))
    extract.add_argument("--known-bad", type=Path, default=DEFAULT_KNOWN_BAD)
    extract.add_argument("--out", type=Path, required=True)
    extract.set_defaults(func=cmd_extract)

    build = sub.add_parser("build", help="import a work directory into a new store")
    build.add_argument("--work", type=Path, required=True)
    build.add_argument("--out", type=Path, default=None,
                       help="new store directory (default: WORK/store-gtars-<version>)")
    build.set_defaults(func=cmd_build)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
