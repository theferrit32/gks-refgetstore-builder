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
from build_store import iter_source_urls, load_config, mirror_cache_path  # noqa: E402

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
    """Execute the ``fetch`` subcommand: pre-populate the cache from the manifest."""
    assemblies, seqsets = load_config(args.config)
    sources = list(iter_source_urls(assemblies, seqsets))
    if args.only:
        want = set(args.only)
        sources = [s for s in sources if s[1] in want]
    if args.kinds:
        kinds = set(args.kinds.split(","))
        sources = [s for s in sources if s[0] in kinds]
    if args.limit:
        sources = sources[: args.limit]

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    print(f"{len(sources)} source file(s); cache={args.cache_dir} "
          f"free={free_gb(args.cache_dir):.1f}GB", file=sys.stderr)

    fetched = skipped = failed = 0
    total_bytes = 0
    failures: list[tuple[str, str]] = []
    for i, (kind, owner, url) in enumerate(sources, 1):
        target = mirror_cache_path(args.cache_dir, url)
        rel = target.relative_to(args.cache_dir)
        if target.exists() and target.stat().st_size > 0:
            skipped += 1
            print(f"[{i}/{len(sources)}] skip  {owner}  {rel}")
            continue
        if args.dry_run:
            print(f"[{i}/{len(sources)}] FETCH {owner}  {rel}")
            continue
        if free_gb(args.cache_dir) < args.min_free_gb:
            print(f"ABORT: free disk {free_gb(args.cache_dir):.1f}GB below "
                  f"--min-free-gb {args.min_free_gb}", file=sys.stderr)
            break
        print(f"[{i}/{len(sources)}] get   {owner}  {rel}", flush=True)
        try:
            n = download(url, target, args.timeout, args.retries)
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
    return 1 if failed else 0
