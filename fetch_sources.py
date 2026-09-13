#!/usr/bin/env python
"""Pre-populate the local FASTA/report cache from sources.toml.

Downloads every remote file the manifest references into a cache that mirrors
the source **host + full path**, e.g.

    downloads/ftp.ncbi.nlm.nih.gov/genomes/all/.../GCF_..._rna.fna.gz
    downloads/ftp.ensembl.org/pub/release-112/.../Homo_sapiens.GRCh38.cdna.all.fa.gz

Mirroring the full path (not just the basename) is required because the NCBI
archive reuses basenames across annotation-release snapshots and per-patch
assembly dirs while the contents differ. ``build_store.py`` reads from the same
mirrored layout (see ``mirror_cache_path``), so a file fetched here is reused by
the build without re-downloading.

Idempotent: existing non-empty files are skipped. Downloads are atomic
(``.part`` then rename) and time-bounded. Failures warn-and-continue so one bad
URL doesn't abort a long run; a summary lists them at the end.

Examples:
    gks-refgetstore fetch --dry-run           # list what would be pulled
    gks-refgetstore fetch --limit 1           # smoke test one file
    gks-refgetstore fetch --only refseq_history_rna
    gks-refgetstore fetch --kinds assembly_report
    gks-refgetstore fetch                     # fetch everything
"""

from __future__ import annotations

import shutil
import socket
import sys
import threading
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_store import (ResolvedSource, load_config, mirror_cache_path,  # noqa: E402
                         provider_checksum_matches, resolve_sources,
                         sha256_file)
from build_lock import apply_locked_sources, load_lock, lock_files  # noqa: E402

CHUNK = 1 << 20  # 1 MiB


def human(nbytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if nbytes < 1024:
            return f"{nbytes:.1f}{unit}"
        nbytes /= 1024
    return f"{nbytes:.1f}PB"


def free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 1e9


def find_partials(cache_dir: Path) -> list[Path]:
    """Return leftover ``.part`` files under the cache, newest last."""
    return sorted(cache_dir.rglob("*.part"), key=lambda p: p.stat().st_mtime)


def sweep_partials(cache_dir: Path, remove: bool) -> tuple[int, int]:
    """Report (and optionally delete) interrupted downloads. Returns (n, bytes).

    Downloads are atomic — data streams to ``<name>.part`` and only becomes the
    real filename via ``replace()`` after a complete read — so a leftover
    ``.part`` can never be mistaken for a finished file. It is dead weight, not
    corruption.

    Cleanup normally happens in ``download``'s exception handler, but that never
    runs if the process is killed, which leaks one file per in-flight download.
    Sweeping at startup keeps an interrupted run from silently accumulating
    dozens of gigabytes of orphans across retries.
    """
    partials = find_partials(cache_dir)
    total = sum(p.stat().st_size for p in partials)
    if remove:
        for path in partials:
            path.unlink(missing_ok=True)
    return len(partials), total


def download(url: str, target: Path, timeout: int, retries: int) -> int:
    """Stream ``url`` to ``target`` atomically. Returns bytes written."""
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + ".part")
    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            written = 0
            req = urllib.request.Request(url, headers={"User-Agent": "gks-refgetstore-builder"})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as fh:  # noqa: S310
                while chunk := resp.read(CHUNK):
                    fh.write(chunk)
                    written += len(chunk)
            tmp.replace(target)
            return written
        except (urllib.error.URLError, socket.timeout, TimeoutError) as exc:
            last_exc = exc
            tmp.unlink(missing_ok=True)
            if attempt < retries:
                print(f"    retry {attempt}/{retries - 1} after {exc}", file=sys.stderr)
    raise RuntimeError(f"failed after {retries} attempts: {last_exc}")


@dataclass
class FetchOutcome:
    """What one concurrent download pass did."""

    fetched: int = 0
    failed: int = 0
    total_bytes: int = 0
    failures: list[tuple[str, str]] = field(default_factory=list)
    aborted_for_disk_space: bool = False

    @property
    def ok(self) -> bool:
        return not self.failed and not self.aborted_for_disk_space


def fetch_many(
    queued: Iterable[tuple[ResolvedSource, Path]],
    cache_dir: Path,
    *,
    timeout: int,
    retries: int,
    jobs: int,
    min_free_gb: float,
    accept: Callable[[Path, ResolvedSource], bool] = provider_checksum_matches,
    label: Callable[[ResolvedSource, Path], str] | None = None,
) -> FetchOutcome:
    """Download ``(source, target)`` pairs concurrently, rejecting bad bytes.

    Providers throttle per connection, not per client -- a second stream
    measured full speed while the first kept its own -- so one serial stream
    leaves most of the link idle. The pool is kept modest to stay polite to
    public FTP endpoints.

    ``accept`` decides whether the downloaded bytes are the ones the caller
    wanted; a rejected file is unlinked and counted as a failure rather than
    left in the cache to be mistaken for a good one. ``fetch`` accepts the
    provider's checksum; ``repair`` accepts only the lock's SHA-256, because
    repair restores the locked state and never re-baselines it.
    """
    outcome = FetchOutcome()
    queued = list(queued)
    if not queued:
        return outcome
    label = label or (lambda source, target: str(target.relative_to(cache_dir)))

    def fetch_one(item: tuple[ResolvedSource, Path]) -> int:
        source, target = item
        print(f"  get   {label(source, target)}", flush=True)
        n = download(source.url, target, timeout, retries)
        if not accept(target, source):
            target.unlink(missing_ok=True)
            raise RuntimeError("downloaded bytes rejected")
        return n

    guard = threading.Lock()
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = {}
        for item in queued:
            if free_gb(cache_dir) < min_free_gb:
                print(f"ABORT: free disk {free_gb(cache_dir):.1f}GB below "
                      f"--min-free-gb {min_free_gb}", file=sys.stderr)
                outcome.aborted_for_disk_space = True
                break
            futures[pool.submit(fetch_one, item)] = item
        for future in as_completed(futures):
            source, target = futures[future]
            try:
                n = future.result()
            except Exception as exc:  # noqa: BLE001
                with guard:
                    outcome.failed += 1
                    outcome.failures.append((source.url, str(exc)))
                print(f"        FAILED {label(source, target)}: {exc}",
                      file=sys.stderr)
                continue
            with guard:
                outcome.fetched += 1
                outcome.total_bytes += n
            print(f"        {human(n)}  {label(source, target)}")
    return outcome


def run_fetch(args) -> int:
    """Pre-populate the mirrored source cache from the manifest.

    Args:
        args: Parsed CLI options for source filters, cache location, retries,
            timeout, dry-run mode, and the minimum free-space threshold.

    Existing nonempty files are retained, while missing files are downloaded
    atomically. Individual download failures are reported after the run and
    produce a nonzero result; insufficient free disk space stops further fetches.
    """
    assemblies, seqsets = load_config(args.config)
    if args.only:
        want = set(args.only)
        assemblies = [entry for entry in assemblies if entry.namespace in want]
        seqsets = [entry for entry in seqsets if entry.name in want]
    lock = load_lock(args.lock) if args.locked_sources else None
    sources = (apply_locked_sources(seqsets, lock) if lock is not None
               else resolve_sources(assemblies, seqsets))
    if args.kinds:
        kinds = set(args.kinds.split(","))
        sources = [s for s in sources if s.kind in kinds]
    if args.limit:
        sources = sources[: args.limit]

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    print(f"{len(sources)} source file(s); cache={args.cache_dir} "
          f"free={free_gb(args.cache_dir):.1f}GB", file=sys.stderr)

    # Sweep before starting: a killed run leaves one .part per in-flight
    # download, and those never get reclaimed otherwise. Dry runs only report,
    # so the check is available without side effects.
    n_partial, partial_bytes = sweep_partials(args.cache_dir, remove=not args.dry_run)
    if n_partial:
        verb = "found" if args.dry_run else "removed"
        print(f"{verb} {n_partial} interrupted download(s), {human(partial_bytes)}",
              file=sys.stderr)

    skipped = 0
    pre_failed: list[tuple[str, str]] = []
    locked_by_url = ({record["url"]: record for record in lock_files(lock)}
                     if lock else {})

    # Pass 1 (serial): classify every source. Checksum verification of cached
    # files is local work, and keeping it ordered keeps the log readable.
    queued: list[tuple[object, Path]] = []
    for i, source in enumerate(sources, 1):
        owner, url = source.owner, source.url
        target = mirror_cache_path(args.cache_dir, url)
        rel = target.relative_to(args.cache_dir)
        if args.locked_sources:
            # A source the lock does not cover is a reportable mismatch, not a
            # crash. Raising a KeyError here aborted mid-loop with some files
            # already fetched and no summary of what had happened.
            locked = locked_by_url.get(url)
            if locked is None:
                reason = "not covered by the lock"
            elif not target.exists() or target.stat().st_size == 0:
                reason = "locked cache file missing"
            elif sha256_file(target) != locked.get("sha256"):
                reason = "locked cache file SHA-256 mismatch"
            else:
                skipped += 1
                print(f"[{i}/{len(sources)}] LOCK  {owner}  {rel}")
                continue
            pre_failed.append((url, reason))
            print(f"[{i}/{len(sources)}] FAIL  {owner}  {rel}  ({reason})")
            continue
        if target.exists() and target.stat().st_size > 0:
            if not provider_checksum_matches(target, source):
                print(f"[{i}/{len(sources)}] stale {owner}  {rel}")
            else:
                skipped += 1
                print(f"[{i}/{len(sources)}] skip  {owner}  {rel}")
                continue
        if args.dry_run:
            print(f"[{i}/{len(sources)}] FETCH {owner}  {rel}")
            continue
        queued.append((source, target))

    # Pass 2 (concurrent), shared with repair.
    outcome = fetch_many(
        queued, args.cache_dir, timeout=args.timeout, retries=args.retries,
        jobs=args.jobs, min_free_gb=args.min_free_gb,
        label=lambda source, target: (
            f"{source.owner}  {target.relative_to(args.cache_dir)}"
        ),
    )
    failures = pre_failed + outcome.failures
    failed = len(pre_failed) + outcome.failed

    print(f"\nfetched={outcome.fetched} skipped={skipped} failed={failed} "
          f"downloaded={human(outcome.total_bytes)}", file=sys.stderr)
    if failures:
        print("failures:", file=sys.stderr)
        for url, exc in failures:
            print(f"  {url}\n    {exc}", file=sys.stderr)
    return 1 if failed or outcome.aborted_for_disk_space else 0
