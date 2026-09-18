#!/usr/bin/env python
"""Inventory every resolved input before download and enforce a space budget."""

from __future__ import annotations

import argparse
import json
import shutil
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .sources import load_config, resolve_sources


def remote_size(url: str) -> int | None:
    request = urllib.request.Request(
        url, method="HEAD", headers={"User-Agent": "gks-refgetstore-builder"}
    )
    for _ in range(3):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
                value = response.headers.get("Content-Length")
                return int(value) if value is not None else None
        except (urllib.error.URLError, TimeoutError, OSError):
            continue
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("sources.toml"))
    parser.add_argument("--out-tsv", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--space-path", type=Path, default=Path("."))
    parser.add_argument("--store-gib", type=float, default=15.0)
    parser.add_argument("--parity-gib", type=float, default=10.0)
    parser.add_argument("--reserve-gib", type=float, default=25.0)
    args = parser.parse_args()

    assemblies, seqsets = load_config(args.config)
    release_by_owner = {entry.name: entry.release for entry in seqsets}
    sources = resolve_sources(assemblies, seqsets)
    with ThreadPoolExecutor(max_workers=8) as executor:
        sizes = list(executor.map(remote_size, (source.url for source in sources)))

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", encoding="utf-8") as out:
        out.write(
            "kind\towner\trelease\tfile_class\tbytes\tprovider_algorithm\t"
            "provider_checksum\tprovider_blocks\tchecksum_url\turl\n"
        )
        for source, size in zip(sources, sizes):
            values = (
                source.kind, source.owner, release_by_owner.get(source.owner),
                source.file_class, size, source.provider_checksum_algorithm
                or ("md5" if source.upstream_md5 else None),
                source.provider_checksum or source.upstream_md5,
                source.provider_checksum_blocks, source.checksum_url, source.url,
            )
            out.write("\t".join("" if value is None else str(value) for value in values) + "\n")

    gib = 1024 ** 3
    input_bytes = sum(size or 0 for size in sizes)
    free_bytes = shutil.disk_usage(args.space_path).free
    required_bytes = input_bytes + int(
        (args.store_gib + args.parity_gib + args.reserve_gib) * gib
    )
    summary = {
        "schema": "gks-refgetstore-source-inventory/1",
        "source_count": len(sources),
        "sized_source_count": sum(size is not None for size in sizes),
        "unsized_source_count": sum(size is None for size in sizes),
        "input_bytes": input_bytes,
        "free_bytes": free_bytes,
        "budget": {
            "store_gib": args.store_gib,
            "parity_gib": args.parity_gib,
            "reserve_gib": args.reserve_gib,
            "required_bytes": required_bytes,
        },
        "space_sufficient": free_bytes >= required_bytes,
    }
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["unsized_source_count"] == 0 and summary["space_sufficient"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
