#!/usr/bin/env python
"""One-shot: reshape the committed build lock from schema /2 to /4.

Throwaway by design. It runs once, its output is committed, and it is retained
only as the record of *how* the committed lock came to be -- not as a migration
path. Schema /4 has no backwards compatibility and no shim; every other reader
in the repo rejects /2 outright.

**Why convert rather than rebuild.** Rebuilding would re-resolve every source
against the providers, and upstream has since changed the MD5 of all 30
remaining ``mRNA_Prot`` shards and withdrawn two. The committed lock preserves
that drift, which is live test material for ``verify --manifest``. A rebuild
would destroy it.

What the conversion does:

* ``build.sources_toml`` -> ``inputs.sources_toml_sha256`` (the absolute path it
  also carried was machine-specific noise in a committed file).
* ``sources[]`` -> ``inputs.files[]``, dropping ``upstream_md5`` (recoverable
  from ``provider_checksum`` + ``provider_checksum_algorithm``) and dropping
  ``collection_digest``/``n_sequences``, which move to ``outputs``.
* ``outputs`` is **enumerated from ``store/``**, not carried over from
  ``build.store``. The old counts were strings, and they describe what the build
  reported rather than what the directory holds.
* ``ingest_spec`` is written ``null`` for every file. That is the truth about
  ``store/``: it was built before the record-exclusion rule existed, so no
  source in it was transformed. Writing the *current* manifest's specs would
  claim an exclusion had been applied when it had not, and ``sync`` would then
  skip exactly the re-ingest that the exclusion requires.

``store/`` has not had the padding sync applied; it holds 1,779,497 sequences
across 240 collections, where ``store.pre-filter/`` holds 1,779,052 across 235.

Usage:
    uv run python runs/2026-09-12-lock-schema-v4/convert.py \
        --in build.lock.json --store store --out build.lock.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

import build_lock  # noqa: E402

V2_SCHEMA = "gks-refgetstore-build-lock/2"


def convert_file_record(source: dict) -> dict:
    """One ``sources[]`` entry as an ``inputs.files[]`` entry."""
    return {
        "kind": source["kind"],
        "owner": source["owner"],
        "url": source["url"],
        "cache_path": source["cache_path"],
        "mutable": source["mutable"],
        "provider_checksum": (
            source.get("provider_checksum") or source.get("upstream_md5")
        ),
        "provider_checksum_algorithm": (
            source.get("provider_checksum_algorithm")
            or ("md5" if source.get("upstream_md5") else None)
        ),
        "provider_checksum_blocks": source.get("provider_checksum_blocks"),
        "checksum_url": source.get("checksum_url"),
        "file_class": source.get("file_class"),
        "bytes": source.get("bytes"),
        "sha256": source.get("sha256"),
        "present": source.get("present", False),
        # Deliberately null: see the module docstring.
        "ingest_spec": None,
        "ingest_spec_sha256": None,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="source", type=Path, required=True)
    ap.add_argument("--store", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--config", type=Path, default=REPO_ROOT / "sources.toml")
    args = ap.parse_args(argv)

    from gtars.refget import RefgetStore

    old = json.loads(args.source.read_text(encoding="utf-8"))
    if old.get("schema") != V2_SCHEMA:
        raise SystemExit(f"expected {V2_SCHEMA}, got {old.get('schema')!r}")

    files = sorted(
        (convert_file_record(source) for source in old["sources"]),
        key=lambda r: (r["kind"], r["cache_path"]),
    )

    # ``from`` inverts the old one-file-one-collection field. The old shape
    # could only express the forward direction; 43 of these digests turn out to
    # have two or more contributors and one has twelve, which is the premise
    # the new schema corrects.
    from_map: dict[str, set[str]] = {}
    for source in old["sources"]:
        digest = source.get("collection_digest")
        if digest:
            from_map.setdefault(digest, set()).add(source["cache_path"])

    store = RefgetStore.open_local(str(args.store))
    store.set_quiet(True)
    outputs = build_lock.store_outputs(store, from_map)

    new = {
        "schema": build_lock.SCHEMA,
        # The original build's metadata is preserved verbatim. This conversion
        # is a reshaping, not a new build, and stamping it with today's git
        # commit would misattribute the store to code that did not produce it.
        "build": {
            "timestamp_utc": old["build"]["timestamp_utc"],
            "git": old["build"]["git"],
            "gtars_version": old["build"]["gtars_version"],
        },
        "inputs": {
            "sources_toml_sha256": old["build"]["sources_toml"]["sha256"],
            "files": files,
        },
        "outputs": outputs,
    }
    build_lock.validate_lock(new)

    unattributed = [c["digest"] for c in outputs["collections"] if not c["from"]]
    claimed = set(from_map) - {c["digest"] for c in outputs["collections"]}
    print(f"inputs.files          {len(files)}")
    print(f"outputs.n_sequences   {outputs['n_sequences']}")
    print(f"outputs.n_collections {outputs['n_collections']}")
    print(f"sequences_root        {outputs['sequences_root']}")
    print(f"collections_root      {outputs['collections_root']}")
    print(f"collections with no known contributor  {len(unattributed)}")
    for digest in unattributed:
        print(f"  {digest}")
    print(f"lock claimed but absent from the store {len(claimed)}")
    for digest in sorted(claimed):
        print(f"  {digest}")
    multi = sum(1 for c in outputs["collections"] if len(c["from"]) > 1)
    widest = max((len(c["from"]) for c in outputs["collections"]), default=0)
    print(f"collections with >1 contributor        {multi} (widest {widest})")

    build_lock.write_lock(args.out, new)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
