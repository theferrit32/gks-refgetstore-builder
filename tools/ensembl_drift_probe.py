#!/usr/bin/env python
"""Rebuild a minimal RefgetStore from cached Ensembl pep FASTAs to isolate
silent sequence drift -- the same accession.version carrying different bytes in
different releases.

This exists to answer the drift question *without* consulting the production
store, so a finding can be confirmed from upstream files alone. It reads only
the download cache and writes only to --out.

    uv run python tools/ensembl_drift_probe.py --releases 112 113 114 115 116 \
        --accessions <file> --out /tmp/drift-store
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from gtars.refget import RefgetStore

CACHE = Path("downloads/ftp.ensembl.org/pub")


def pep_path(release: int) -> Path:
    return (CACHE / f"release-{release}/fasta/homo_sapiens/pep"
            / "Homo_sapiens.GRCh38.pep.all.fa.gz")


def extract(release: int, wanted: set[str], out_dir: Path) -> tuple[Path, int]:
    """Write a FASTA holding just ``wanted`` from one release's pep file."""
    source = pep_path(release)
    target = out_dir / f"ensembl-{release}.pep.subset.fa"
    kept = 0
    with gzip.open(source, "rt") as handle, target.open("w") as sink:
        emit = False
        for line in handle:
            if line.startswith(">"):
                emit = line[1:].split()[0] in wanted
                kept += emit
            if emit:
                sink.write(line)
    return target, kept


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--releases", type=int, nargs="+", required=True)
    parser.add_argument("--accessions", type=Path, required=True,
                        help="one accession.version per line")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    wanted = {line.strip() for line in args.accessions.read_text().splitlines()
              if line.strip() and not line.startswith("#")}
    args.out.mkdir(parents=True, exist_ok=True)
    store = RefgetStore.on_disk(str(args.out))

    # Record release -> collection digest as the import happens. Two releases
    # routinely share an identical *name* set while differing in sequence, so a
    # digest cannot be recovered afterwards by matching accession lists.
    manifest: dict[str, str] = {}
    for release in sorted(args.releases):
        fasta, kept = extract(release, wanted, args.out)
        # gtars returns (metadata, created?); only the digest matters here.
        metadata, _ = store.add_sequence_collection_from_fasta(str(fasta))
        digest = metadata.digest
        manifest[str(release)] = digest
        print(f"ensembl-{release}: {kept} records -> {digest}")
    store.write()
    (args.out / "releases.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"store written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
