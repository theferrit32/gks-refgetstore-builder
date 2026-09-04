#!/usr/bin/env python
"""gks-refgetstore: build and verify a gtars RefgetStore from an authoritative
source manifest (``sources.toml``).

Subcommands:
  build    ingest the manifest's sources into a RefgetStore (+ lock write/check)
  verify   verify the local cache: integrity, or per-file drift vs a build lock
  fetch    pre-populate the download cache from the manifest (no ingest)
  lock     (re)generate a build.lock.json for a build that already ran

All subcommands share one library: config loading + cache-path mapping
(``build_store``), the lock/verify core (``build_lock``), and the fetcher
(``fetch_sources``). Runnable via the ``gks-refgetstore`` console script or
directly with ``uv run cli.py <subcommand> ...``.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import build_lock
import build_store
import fetch_sources
import store_sync

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = REPO_ROOT / "sources.toml"
DEFAULT_STORE = REPO_ROOT / "store"
DEFAULT_CACHE = REPO_ROOT / "downloads"
DEFAULT_LOCK = REPO_ROOT / "build.lock.json"


def _add_build(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("build", help="ingest the manifest into a RefgetStore")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--store-dir", type=Path, default=DEFAULT_STORE)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--assembly", default=None,
                   help="process only this assembly namespace (skips seqsets)")
    p.add_argument("--seqset", default=None,
                   help="process only this seqset name (skips assemblies)")
    p.add_argument("--skip-assemblies", action="store_true")
    p.add_argument("--skip-seqsets", action="store_true")
    p.add_argument("--force-download", action="store_true")
    p.add_argument("--ingest-jobs", type=int,
                   default=build_store.INGEST_JOBS_DEFAULT,
                   help="FASTA files imported concurrently per seqset "
                        f"(default {build_store.INGEST_JOBS_DEFAULT}; 1=serial). "
                        "Lower it if peak memory is a concern; the largest "
                        "Ensembl inputs cost ~0.45 GiB per concurrent file")
    p.add_argument("--filter-jobs", type=int,
                   default=build_store.INGEST_JOBS_DEFAULT,
                   help="source files filtered concurrently in the preflight "
                        f"(default {build_store.INGEST_JOBS_DEFAULT}; 1=serial). "
                        "Separate from --ingest-jobs: filtering is zlib-bound "
                        "and cheap in memory, so capping ingest concurrency to "
                        "limit peak RSS should not serialize it")
    p.add_argument("--min-free-gb", type=float, default=25.0,
                   help="stop before ingesting if free disk is below this "
                        "(default 25)")
    p.add_argument("--lock", type=Path, default=DEFAULT_LOCK,
                   help="build-lock path to check against and write/refresh")
    p.add_argument("--no-lock", action="store_true",
                   help="don't write the lock (the pre-flight check still runs)")
    p.add_argument("--lock-check-mode", choices=["strict", "subset", "ignore"],
                   default="strict",
                   help="pre-flight check of to-be-ingested files vs the lock: "
                        "strict (default)=stop on membership or content drift; "
                        "ignore=skip")
    p.add_argument("--force-lock", action="store_true",
                   help="write the lock even if a discrepancy would suppress it")
    p.add_argument("--locked-sources", action="store_true",
                   help="use concrete URLs and SHA-256 values from --lock; no discovery")
    p.set_defaults(func=build_store.run_build)


def _add_sync(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "sync",
        help="apply sources.toml changes to an existing store incrementally",
        description="Classify every manifest source against the build lock, "
                    "then remove and re-ingest only what changed. Sources whose "
                    "upstream bytes and ingest spec both match the lock are "
                    "validated and left untouched.",
    )
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--store-dir", type=Path, default=DEFAULT_STORE)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    p.add_argument("--apply", action="store_true",
                   help="mutate the store. Without it, sync classifies and "
                        "exits. Removal persists immediately and cannot be "
                        "rolled back, so this is never the default")
    p.add_argument("--no-hash", action="store_true",
                   help="skip re-hashing cached sources; classify on the "
                        "ingest spec alone. Faster, but will not notice "
                        "upstream bytes that changed in place")
    p.add_argument("--no-lock", action="store_true",
                   help="apply to the store but don't rewrite the lock")
    p.add_argument("--ingest-jobs", type=int,
                   default=build_store.INGEST_JOBS_DEFAULT)
    p.add_argument("--filter-jobs", type=int,
                   default=build_store.INGEST_JOBS_DEFAULT)
    p.set_defaults(func=store_sync.run_sync)


def _add_verify(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("verify", help="verify the cache (integrity or vs a lock)")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--lock", type=Path, default=None,
                   help="per-file drift check against this build.lock.json "
                        "(omit for plain integrity checks)")
    p.add_argument("--strict-set", action="store_true",
                   help="with --lock, also fail if the build/lock file sets differ")
    p.add_argument("--check-remote", action="store_true",
                   help="also compare local size to the server's Content-Length")
    p.add_argument("--timeout", type=int, default=30)
    p.add_argument("--limit", type=int, default=None, help="only check first N")
    p.set_defaults(func=build_lock.run_verify)


def _add_fetch(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("fetch", help="pre-populate the download cache")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    p.add_argument("--locked-sources", action="store_true",
                   help="verify and use cached concrete sources from --lock; no network")
    p.add_argument("--only", nargs="*", default=None,
                   help="only fetch sources whose owner matches")
    p.add_argument("--kinds", default=None,
                   help="comma-separated: seqset,assembly_fasta,assembly_report")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--retries", type=int, default=3)
    p.add_argument("--min-free-gb", type=float, default=10.0)
    p.add_argument("--jobs", type=int, default=6,
                   help="concurrent downloads (default 6). Providers throttle "
                        "per connection, so one stream leaves most of the link "
                        "idle; keep this modest to stay polite to public FTP")
    p.set_defaults(func=fetch_sources.run_fetch)


def _add_lock(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("lock", help="(re)generate build.lock.json from a build log")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--store-dir", type=Path, default=DEFAULT_STORE)
    p.add_argument("--from-log", type=Path, required=True,
                   help="build log to reconstruct file->collection mapping from")
    p.add_argument("--out", type=Path, default=DEFAULT_LOCK)
    p.set_defaults(func=build_lock.run_lock)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="gks-refgetstore", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    sub = ap.add_subparsers(dest="command", required=True)
    _add_build(sub)
    _add_verify(sub)
    _add_fetch(sub)
    _add_lock(sub)
    _add_sync(sub)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
