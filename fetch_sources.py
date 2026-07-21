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
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_store import (apply_locked_sources, load_config, md5_file,
                         mirror_cache_path, resolve_sources)  # noqa: E402
from build_lock import load_lock, sha256_file  # noqa: E402

CHUNK = 1 << 20  # 1 MiB


def human(nbytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if nbytes < 1024:
            return f"{nbytes:.1f}{unit}"
        nbytes /= 1024
    return f"{nbytes:.1f}PB"


def free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 1e9


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
    lock = load_lock(args.lock) if args.locked_sources else None
    sources = (apply_locked_sources(seqsets, lock) if lock is not None
               else resolve_sources(assemblies, seqsets))
    if args.only:
        want = set(args.only)
        sources = [s for s in sources if s.owner in want]
    if args.kinds:
        kinds = set(args.kinds.split(","))
        sources = [s for s in sources if s.kind in kinds]
    if args.limit:
        sources = sources[: args.limit]

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    print(f"{len(sources)} source file(s); cache={args.cache_dir} "
          f"free={free_gb(args.cache_dir):.1f}GB", file=sys.stderr)

    fetched = skipped = failed = 0
    aborted_for_disk_space = False
    total_bytes = 0
    failures: list[tuple[str, str]] = []
    locked_by_url = ({entry["url"]: entry for entry in lock.get("sources", [])}
                     if lock else {})
    for i, source in enumerate(sources, 1):
        kind, owner, url = source.kind, source.owner, source.url
        target = mirror_cache_path(args.cache_dir, url)
        rel = target.relative_to(args.cache_dir)
        if args.locked_sources:
            locked = locked_by_url[url]
            if (not target.exists() or target.stat().st_size == 0
                    or sha256_file(target) != locked.get("sha256")):
                failed += 1
                failures.append((url, "locked cache file missing or SHA-256 mismatch"))
                print(f"[{i}/{len(sources)}] FAIL  {owner}  {rel}")
            else:
                skipped += 1
                print(f"[{i}/{len(sources)}] LOCK  {owner}  {rel}")
            continue
        if target.exists() and target.stat().st_size > 0:
            if source.upstream_md5 and md5_file(target) != source.upstream_md5:
                print(f"[{i}/{len(sources)}] stale {owner}  {rel}")
            else:
                skipped += 1
                print(f"[{i}/{len(sources)}] skip  {owner}  {rel}")
                continue
        if args.dry_run:
            print(f"[{i}/{len(sources)}] FETCH {owner}  {rel}")
            continue
        if free_gb(args.cache_dir) < args.min_free_gb:
            print(f"ABORT: free disk {free_gb(args.cache_dir):.1f}GB below "
                  f"--min-free-gb {args.min_free_gb}", file=sys.stderr)
            aborted_for_disk_space = True
            break
        print(f"[{i}/{len(sources)}] get   {owner}  {rel}", flush=True)
        try:
            n = download(url, target, args.timeout, args.retries)
            if source.upstream_md5 and md5_file(target) != source.upstream_md5:
                target.unlink(missing_ok=True)
                raise RuntimeError("NCBI MD5 mismatch")
            fetched += 1
            total_bytes += n
            print(f"           {human(n)}  (cache free {free_gb(args.cache_dir):.1f}GB)")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            failures.append((url, str(exc)))
            print(f"           FAILED: {exc}", file=sys.stderr)

    print(f"\nfetched={fetched} skipped={skipped} failed={failed} "
          f"downloaded={human(total_bytes)}", file=sys.stderr)
    if failures:
        print("failures:", file=sys.stderr)
        for url, exc in failures:
            print(f"  {url}\n    {exc}", file=sys.stderr)
    return 1 if failed or aborted_for_disk_space else 0
