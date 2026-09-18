#!/usr/bin/env python
"""gks-refgetstore: build and verify a gtars RefgetStore from an authoritative
source manifest (``sources.toml``).

Subcommands:
  build    ingest the manifest's sources into a RefgetStore (+ lock write/check)
  verify   check the cache and the store against the build lock (judges; gates CI)
  status   describe manifest/lock/cache/store agreement (describes; always exit 0)
  repair   restore the locked state after verify finds damage
  sync     apply sources.toml changes to an existing store incrementally
  fetch    pre-populate the download cache from the manifest (no ingest)
  lock     (re)generate a build.lock.json for a build that already ran

All subcommands share one library: the source model -- manifest schema, source
resolution, cache-path mapping (``sources``) -- the build engine
(``build_store``), the lock core (``build_lock``), the read-only store
primitives (``store_census``), and the fetcher (``fetch_sources``). Runnable via
the ``gks-refgetstore`` console script or with
``python -m gks_refgetstore.cli <subcommand> ...``.

Relative path defaults resolve against the current directory, so run this from
the repo root or pass ``--config``/``--store-dir``/``--cache-dir``/``--lock``.

For a fast edit-test loop, point everything at the development manifest, which
covers every config shape in minutes rather than hours (measured: build ~6m,
``verify --all`` 23s):

    gks-refgetstore build --config sources.dev.toml \\
        --store-dir store.dev --lock build.dev.lock.json
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from . import build_lock, build_store, fetch_sources, repair, store_sync, verify

# Relative to the current directory, not to the package: the manifest, the
# cache and the store are the *user's* data, and an installed tool has no
# business reaching back into its own install tree for them. Run from the repo
# root, or pass the paths explicitly.
DEFAULT_CONFIG = Path("sources.toml")
DEFAULT_STORE = Path("store")
DEFAULT_CACHE = Path("downloads")
DEFAULT_LOCK = Path("build.lock.json")


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
    p = sub.add_parser(
        "verify", help="verify the cache and the store against the build lock",
        description="Four independent targets. With no scope flag the default "
                    "is --cache --store, which is offline and lock-driven; "
                    "only --remote and --manifest touch the network. A stale "
                    "lock and a broken store are different failures, and the "
                    "split is what tells them apart.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--store-dir", type=Path, default=DEFAULT_STORE)
    p.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    p.add_argument("--cache", action="store_true",
                   help="the pinned source bytes are present and intact (offline)")
    p.add_argument("--store", action="store_true",
                   help="the store holds exactly what the lock recorded (offline)")
    p.add_argument("--manifest", action="store_true",
                   help="upstream still publishes what the lock pins (network)")
    p.add_argument("--remote", action="store_true",
                   help="remote objects are still the size we saw (network). "
                        "Findings are warnings: several endpoints omit "
                        "Content-Length and mutable ones legitimately change")
    p.add_argument("--all", action="store_true", help="all four targets")
    p.add_argument("--deep", action="store_true",
                   help="with --store, re-digest every sequence from the bytes "
                        "the store serves (L2; ~18 min on the full store)")
    p.add_argument("--no-hash", action="store_true",
                   help="with --cache, skip SHA-256 (presence and gzip only)")
    p.add_argument("--known-bad", type=Path, default=verify.DEFAULT_KNOWN_BAD,
                   help="round-trip baseline TSV for --deep. Digests listed "
                        "here report as warnings with their diagnosed cause; "
                        "anything not listed is an error")
    p.add_argument("--jobs", type=int, default=6,
                   help="concurrency for cache hashing and remote probes "
                        "(default 6). L2 is single-threaded by design")
    p.add_argument("--timeout", type=int, default=30)
    p.add_argument("--limit", type=int, default=None,
                   help="only check the first N files; rejected with --store, "
                        "where a partial root is meaningless")
    p.set_defaults(func=verify.run_verify)


def _add_status(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "status", help="describe manifest/lock/cache/store agreement",
        description="Four legs: manifest vs lock, manifest vs cache, lock vs "
                    "cache, lock vs store. Always exits 0 -- status describes, "
                    "verify judges.",
    )
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--store-dir", type=Path, default=DEFAULT_STORE)
    p.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    p.add_argument("--offline", action="store_true",
                   help="skip the one leg that resolves the manifest upstream")
    p.set_defaults(func=verify.run_status)


def _add_repair(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "repair", help="restore the locked state after verify finds damage",
        description="Re-fetches cache files to their locked SHA-256 and "
                    "re-ingests collections the store has lost. Restores what "
                    "the lock records; it never re-baselines. New upstream "
                    "bytes are a `sync` decision, not a repair.",
    )
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--store-dir", type=Path, default=DEFAULT_STORE)
    p.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    p.add_argument("--cache", action="store_true", help="repair cache findings")
    p.add_argument("--store", action="store_true", help="repair store findings")
    p.add_argument("--apply", action="store_true",
                   help="actually repair. Without it, repair plans and exits. "
                        "Store removal persists immediately and cannot be "
                        "rolled back, so this is never the default")
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--retries", type=int, default=3)
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--min-free-gb", type=float, default=10.0)
    p.add_argument("--ingest-jobs", type=int,
                   default=build_store.INGEST_JOBS_DEFAULT)
    p.add_argument("--filter-jobs", type=int,
                   default=build_store.INGEST_JOBS_DEFAULT)
    p.set_defaults(func=repair.run_repair)


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
    _add_status(sub)
    _add_repair(sub)
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
