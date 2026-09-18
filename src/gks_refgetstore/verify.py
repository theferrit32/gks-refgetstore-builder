#!/usr/bin/env python
"""Verify that the cache and the store still match what the build lock recorded.

Four independent targets, each answering a different question:

| target       | question                                          | network |
|--------------|---------------------------------------------------|---------|
| ``--cache``  | are the pinned source bytes still on disk, intact? | no     |
| ``--store``  | does the store still hold what the lock recorded?  | no     |
| ``--manifest``| has upstream moved since the lock was written?    | yes    |
| ``--remote`` | are the remote objects still the size we saw?      | yes    |

Default with no scope flags is ``--cache --store``: **offline and lock-driven**.
There is deliberately no manifest-driven integrity mode. The previous one called
``resolve_sources``, which fetches provider manifests, so "verify" could not run
offline and its result depended on what upstream happened to be serving that
day -- the exact confusion this separation exists to remove. A legitimately
stale lock is not a broken store, and the three-way split says which is which.

Every finding carries a stable ``code``. That is the contract ``repair``
dispatches on, so codes are part of the interface, not log text.

## The store ladder

| level | proves | cost |
|---|---|---|
| **L0** roots + counts | the store holds **exactly** the digest sets the lock recorded; any single add or remove flips a root | ~7s |
| **L1** membership | every locked collection is present and fully resolvable, every sequence is backed by a ``.seq`` file, nothing is stranded | ~85s |
| **L2** ``--deep`` | the store can actually serve each sequence, and returns what its digest promises | ~18 min |

The levels are independent, which is the point: unlinking a ``.seq`` fires L1
while L0 stays clean (the index is untouched), and corrupting one fires L2 while
L0 and L1 stay clean (the file is still there and still the right length).

L2 needs no lock at all -- it compares the store against itself.

## Known encoder defects

L2 on the real store surfaces 133 sequences whose stored payload does not
re-digest to the digest it is filed under. They are not bit rot: gtars' protein
alphabet has no ``U`` (selenocysteine), and 27 short proteins whose residues are
all legal IUPAC nucleotide codes are misdetected as nucleotide and encoded
lossily. Re-ingesting reproduces both, because they are encoder limitations.

They are pinned as a **baseline** rather than suppressed, which gives three
outcomes instead of two:

| in baseline | re-digests correctly | outcome |
|---|---|---|
| yes | no | ``warn`` -- known defect, rendered with its diagnosed cause |
| **no** | **no** | **``error``** -- a new defect, the regression we care about |
| **yes** | **yes** | ``info`` -- gtars fixed it; update the baseline |

That third row is why a baseline beats a suppression flag: it turns the upstream
fix into something the tool reports rather than something we have to remember to
check.
"""

from __future__ import annotations

import csv
import subprocess
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

from . import build_lock, store_census
from .sources import load_config, resolve_sources, sha256_file

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_KNOWN_BAD = (
    REPO_ROOT / "seqrepo_equivalence" / "known_divergence"
    / "gtars_encoding_roundtrip.tsv"
)

ERROR, WARN, INFO = "error", "warn", "info"


@dataclass(frozen=True)
class Finding:
    """One thing verify noticed, in the form repair can act on."""

    code: str
    severity: str
    subject: str
    detail: str = ""
    lines: tuple[str, ...] = ()

    def render(self) -> str:
        head = f"[{self.severity:5}] {self.code}  {self.subject}"
        if self.detail:
            head += f"\n        {self.detail}"
        return "\n".join([head, *(f"        {line}" for line in self.lines)])


@dataclass
class VerifyReport:
    findings: list[Finding] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def add(self, *findings: Finding) -> None:
        self.findings.extend(findings)

    def note(self, *lines: str) -> None:
        self.notes.extend(lines)

    def count(self, severity: str) -> int:
        return sum(1 for f in self.findings if f.severity == severity)

    @property
    def failed(self) -> bool:
        return self.count(ERROR) > 0

    def merge(self, other: "VerifyReport") -> None:
        self.findings.extend(other.findings)
        self.notes.extend(other.notes)


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}PB"


def elapsed(seconds: float) -> str:
    minutes, secs = divmod(int(seconds), 60)
    return f"{minutes}m{secs:02d}s" if minutes else f"{secs}s"


# ---------------------------------------------------------------------- cache

def gzip_intact(path: Path) -> tuple[bool, str]:
    """``gzip -t`` -- validates the decompressed CRC32 + length trailer."""
    proc = subprocess.run(["gzip", "-t", str(path)], capture_output=True, text=True)
    if proc.returncode == 0:
        return True, "gzip ok"
    return False, (proc.stderr.strip() or f"gzip -t exit {proc.returncode}")


def verify_cache(
    lock: dict, cache_dir: Path, *, hash_files: bool = True,
    limit: int | None = None, jobs: int = 1,
) -> VerifyReport:
    """Check every file the lock pins: present, right bytes, gzip-intact.

    Purely local. The lock's ``sha256`` is the authority, not the provider's
    checksum -- the provider may have republished since, which is
    ``--manifest``'s business, not this one's.
    """
    report = VerifyReport()
    records = [r for r in build_lock.lock_files(lock) if r.get("present")]
    skipped = len(build_lock.lock_files(lock)) - len(records)
    if limit:
        records = records[:limit]

    def check(record: dict) -> list[Finding]:
        rel = record["cache_path"]
        path = cache_dir / rel
        if not path.exists() or path.stat().st_size == 0:
            return [Finding("cache.missing", ERROR, rel,
                            "the lock records this file as present")]
        found: list[Finding] = []
        size = path.stat().st_size
        if record.get("bytes") is not None and size != record["bytes"]:
            found.append(Finding(
                "cache.size_mismatch", ERROR, rel,
                f"lock={human(record['bytes'])} actual={human(size)}",
            ))
        if hash_files and record.get("sha256"):
            actual = sha256_file(path)
            if actual != record["sha256"]:
                found.append(Finding(
                    "cache.sha256_mismatch", ERROR, rel,
                    f"lock={record['sha256'][:16]}… actual={actual[:16]}…",
                ))
        if path.suffix == ".gz":
            good, message = gzip_intact(path)
            if not good:
                found.append(Finding("cache.gzip_corrupt", ERROR, rel, message))
        return found

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for found in pool.map(check, records):
            report.add(*found)

    clean = len(records) - len({f.subject for f in report.findings})
    report.note(
        f"cache   {len(records)} locked file(s) checked"
        f"{' (hashing skipped)' if not hash_files else ''}",
        f"        {clean} ok, {len({f.subject for f in report.findings})} with findings",
    )
    if skipped:
        report.note(f"        {skipped} lock record(s) skipped (present=false)")
    return report


# ---------------------------------------------------------------------- store

def load_known_bad(path: Path) -> dict[str, dict]:
    """Read the round-trip baseline TSV, keyed by digest. Missing file = empty."""
    if not path or not Path(path).exists():
        return {}
    rows: dict[str, dict] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            reader = csv.DictReader([line], delimiter="\t",
                                    fieldnames=KNOWN_BAD_COLUMNS)
            row = next(reader)
            if row["digest"] == "digest":  # the header row
                continue
            rows[row["digest"]] = row
    return rows


KNOWN_BAD_COLUMNS = ("digest", "name", "alphabet", "length", "redigest", "cause")


def verify_store_l0(lock: dict, store) -> VerifyReport:
    """Roots and counts: does the store hold exactly the recorded digest sets?

    A root is a digest over the whole *set*, so it cannot be partially right:
    one sequence added, removed, or re-keyed flips it. It says nothing about
    whether any payload is readable, which is what L1 and L2 are for.
    """
    report = VerifyReport()
    outputs = lock.get("outputs", {})
    census = store_census.collection_census(store)
    sequences = store_census.sequence_digests(store)
    actual = {
        "n_sequences": len(sequences),
        "n_collections": len(census),
        "sequences_root": store_census.digest_root(sequences),
        "collections_root": store_census.digest_root(census),
    }
    for key in ("n_sequences", "n_collections"):
        if outputs.get(key) != actual[key]:
            report.add(Finding(
                f"store.{key}_mismatch", ERROR, key,
                f"lock={outputs.get(key)} actual={actual[key]}",
            ))
    for key in ("sequences_root", "collections_root"):
        if outputs.get(key) != actual[key]:
            report.add(Finding(
                f"store.{key}_mismatch", ERROR, key,
                f"lock={outputs.get(key)}",
                (f"actual={actual[key]}",
                 "the store's digest set is not the one the lock recorded"),
            ))
    report.note(
        f"store   L0 roots+counts: {actual['n_sequences']:,} sequences, "
        f"{actual['n_collections']} collections"
    )
    return report


def verify_store_l1(lock: dict, store, store_dir: Path) -> VerifyReport:
    """Membership: every locked collection resolvable, every sequence backed.

    Two questions that are deliberately kept apart:

    * **Does the store have what the lock recorded?** Per locked collection:
      present, right size, level-2 arrays resolvable.
    * **Is the store internally coherent?** Every indexed sequence backed by a
      payload file, every payload file indexed, nothing stranded outside all
      collections. This half needs no lock, so a lock that under-describes the
      store cannot turn one missing collection into a flood of spurious orphan
      findings for the sequences it published.

    One ``scandir`` pass over ``sequences/`` rather than a ``Path.exists()`` per
    digest -- at 1.8M sequences that is 5.8 seconds against several minutes.
    """
    report = VerifyReport()
    census = store_census.collection_census(store)
    indexed = store_census.sequence_digests(store)
    on_disk = store_census.on_disk_sequence_digests(store_dir)
    locked = build_lock.lock_collections(lock)

    for collection in locked:
        digest = collection["digest"]
        if digest not in census:
            report.add(Finding(
                "store.collection_missing", ERROR, digest,
                f"the lock records {collection['n_sequences']} sequence(s) in "
                "this collection; the store has no such collection",
                tuple(f"from: {rel}" for rel in collection["from"]) or
                ("from: (no known contributor)",),
            ))
            continue
        if census[digest] != collection["n_sequences"]:
            report.add(Finding(
                "store.collection_size_mismatch", ERROR, digest,
                f"lock={collection['n_sequences']} actual={census[digest]}",
            ))

    # Collections the store holds but the lock does not describe. A warning, not
    # an error -- L0's collections_root already failed if it matters -- but it
    # is the finding that *explains* that root mismatch.
    for digest in sorted(set(census) - {c["digest"] for c in locked}):
        report.add(Finding(
            "store.collection_unrecorded", WARN, digest,
            f"{census[digest]} sequence(s) in the store, absent from the lock"))

    reachable: set[str] = set()
    for digest in sorted(census):
        try:
            reachable.update(store_census.collection_members(store, digest))
        except Exception as exc:  # noqa: BLE001
            report.add(Finding(
                "store.collection_unresolvable", ERROR, digest, str(exc)))

    for digest in sorted(indexed - on_disk):
        report.add(Finding(
            "store.seq_file_missing", ERROR, digest,
            f"no payload at {store_census.seq_path(store_dir, digest)}"))
    for digest in sorted(indexed - reachable):
        report.add(Finding(
            "store.orphan_sequence", WARN, digest,
            "indexed by the store but reachable from no collection in it"))
    for digest in sorted(on_disk - indexed):
        report.add(Finding(
            "store.stray_seq_file", WARN, digest,
            f"payload on disk that the store does not index: "
            f"{store_census.seq_path(store_dir, digest)}"))

    report.note(
        f"        L1 membership: {len(locked)} locked collection(s); "
        f"{len(reachable):,} sequence(s) reachable, {len(on_disk):,} payload "
        f"file(s) on disk"
    )
    return report


def verify_store_l2(
    store, known_bad: dict[str, dict], *, progress_every: int = 100_000,
) -> VerifyReport:
    """Deep: re-digest every stored sequence from the bytes the store serves.

    Needs no lock -- it compares the store against itself. The digest a sequence
    is filed under was computed from the source bytes at ingest; this recomputes
    it from what can be read back today. A disagreement means the store cannot
    serve that sequence correctly, whether the cause is bit rot or a lossy
    encoder.

    Single-threaded at ~1,600 sequences/second, about 18 minutes on the
    production store. **Do not parallelize this without first measuring whether
    gtars releases the GIL inside ``stream_sequence``**; if it does not, a
    thread pool buys nothing and a process pool would reopen the store per
    worker.
    """
    report = VerifyReport()
    metadata = store_census.sequence_metadata_by_digest(store)
    started = time.monotonic()
    ok = new_bad = 0
    known_hit: set[str] = set()
    by_cause: dict[str, int] = {}

    for i, (digest, meta) in enumerate(sorted(metadata.items()), 1):
        if i % progress_every == 0:
            rate = i / max(1e-9, time.monotonic() - started)
            print(f"        …{i:,}/{len(metadata):,} ({rate:,.0f}/s)",
                  file=sys.stderr, flush=True)
        try:
            actual = store_census.redigest_sequence(store, digest)
        except Exception as exc:  # noqa: BLE001
            report.add(Finding(
                "store.sequence_unreadable", ERROR, digest,
                f"{meta.name}   {meta.alphabet}   {meta.length} aa/bp",
                (str(exc),)))
            continue
        if actual == digest:
            ok += 1
            if digest in known_bad:
                known_hit.add(digest)
                report.add(Finding(
                    "store.seq_digest_fixed", INFO, digest,
                    f"{meta.name}   {meta.alphabet}   {meta.length}",
                    ("this sequence is in the known-bad baseline but now "
                     "re-digests correctly",
                     "gtars appears to have been fixed; drop it from the "
                     "baseline and re-ingest its collections"),
                ))
            continue
        baseline = known_bad.get(digest)
        if baseline:
            known_hit.add(digest)
            by_cause[baseline["cause"]] = by_cause.get(baseline["cause"], 0) + 1
            report.add(Finding(
                "store.seq_digest_mismatch", WARN, digest,
                f"{meta.name}   {meta.alphabet}   {meta.length}",
                (f"decoded bytes re-digest to {actual}",
                 f"known cause: {baseline['cause']}"),
            ))
        else:
            new_bad += 1
            report.add(Finding(
                "store.seq_digest_mismatch", ERROR, digest,
                f"{meta.name}   {meta.alphabet}   {meta.length}",
                (f"decoded bytes re-digest to {actual}",
                 "not in the known-bad baseline: this is a new defect"),
            ))

    took = time.monotonic() - started
    report.note(
        f"        L2 deep: {len(metadata):,} sequence(s) re-digested in "
        f"{elapsed(took)}",
        f"        {ok:,} ok",
    )
    # Broken down by cause, because "133 known" is a number and
    # "106 selenoproteins + 27 misdetected alphabets" is a description. A count
    # moving between the two groups is a signal the total would hide.
    for cause, count in sorted(by_cause.items(), key=lambda kv: -kv[1]):
        report.note(f"        {count:>5} {cause.split(';')[0]}   [baseline]")
    report.note(f"        {new_bad:>5} new")
    stale = set(known_bad) - set(metadata)
    if stale:
        report.note(
            f"        {len(stale)} baseline entr(y/ies) name digests this store "
            "does not hold"
        )
    return report


def verify_store(
    lock: dict, store_dir: Path, *, deep: bool = False,
    known_bad_path: Path | None = None,
) -> VerifyReport:
    """Run the store ladder. L0 and L1 need the lock; L2 does not."""
    from gtars.refget import RefgetStore

    report = VerifyReport()
    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    report.merge(verify_store_l0(lock, store))
    report.merge(verify_store_l1(lock, store, store_dir))
    if deep:
        known_bad = load_known_bad(known_bad_path or DEFAULT_KNOWN_BAD)
        report.merge(verify_store_l2(store, known_bad))
    return report


# --------------------------------------------------------------------- remote

def remote_size(url: str, timeout: int) -> int | None:
    """Content-Length from an HTTP HEAD, or None if unavailable."""
    request = urllib.request.Request(
        url, method="HEAD", headers={"User-Agent": "gks-refgetstore-builder"}
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
            length = response.headers.get("Content-Length")
            return int(length) if length is not None else None
    except (urllib.error.URLError, ValueError, TimeoutError, OSError):
        return None


def verify_remote(
    lock: dict, *, timeout: int = 30, limit: int | None = None, jobs: int = 6,
) -> VerifyReport:
    """Compare each locked byte count to the server's ``Content-Length``.

    Every finding is ``warn``, never ``error``. Several endpoints omit
    ``Content-Length`` entirely, and the mutable endpoints legitimately change
    size between releases -- so a size difference is information, not a verdict.
    ``--manifest`` answers the same question better, at the same network cost,
    by comparing provider checksums instead of guessing from sizes.
    """
    report = VerifyReport()
    records = [r for r in build_lock.lock_files(lock) if r.get("bytes")]
    if limit:
        records = records[:limit]

    def check(record: dict) -> Finding | None:
        size = remote_size(record["url"], timeout)
        if size is None:
            return Finding("remote.size_unavailable", WARN, record["cache_path"],
                           "the server sent no Content-Length")
        if size != record["bytes"]:
            return Finding(
                "remote.size_mismatch", WARN, record["cache_path"],
                f"lock={human(record['bytes'])} remote={human(size)}"
                + ("  (mutable endpoint)" if record.get("mutable") else ""),
            )
        return None

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for finding in pool.map(check, records):
            if finding is not None:
                report.add(finding)
    report.note(f"remote  {len(records)} URL(s) probed")
    return report


# ------------------------------------------------------------------- manifest

def verify_manifest(lock: dict, config_path: Path) -> VerifyReport:
    """Has upstream moved since the lock was written?

    The only target that resolves the manifest live, which is why it is opt-in.
    A finding here means the *lock* is stale, not that the store is broken --
    the distinction the whole four-target split exists to make.
    """
    report = VerifyReport()
    if sha256_file(config_path) != build_lock.lock_sources_toml_sha256(lock):
        report.add(Finding(
            "manifest.config_changed", WARN, str(config_path),
            "sources.toml has been edited since the lock was written",
        ))
    assemblies, seqsets = load_config(config_path)
    sources = resolve_sources(assemblies, seqsets)
    check = build_lock.evaluate_sources_vs_lock(sources, lock)
    by_url = {r["url"]: r for r in build_lock.lock_files(lock)}
    live_by_url = {s.url: s for s in sources}

    for url in check.new:
        report.add(Finding(
            "manifest.file_added", ERROR, url,
            f"published by {live_by_url[url].owner} since the lock was written",
        ))
    for url in check.only_in_lock:
        report.add(Finding(
            "manifest.file_removed", ERROR, url,
            f"the lock pins it for {by_url[url]['owner']}; upstream no longer "
            "lists it",
        ))
    for url, was, now in check.upstream_changed:
        report.add(Finding(
            "manifest.checksum_changed", ERROR, url,
            f"provider checksum {was} -> {now}"
            + ("  (mutable endpoint)" if by_url[url].get("mutable") else ""),
        ))
    report.note(
        f"manifest {len(sources)} resolved source(s); "
        f"added={len(check.new)} removed={len(check.only_in_lock)} "
        f"checksum_changed={len(check.upstream_changed)}"
    )
    return report


# ------------------------------------------------------------------------ cli

PER_CODE_LIMIT = 20


def render(report: VerifyReport, *, show_info: bool = True,
           per_code: int = PER_CODE_LIMIT) -> None:
    """Print findings worst-first, capped per code.

    The cap matters: one collection missing from the lock strands every sequence
    it published, and printing 200,000 identical-shaped lines buries the single
    finding that explains them. The count is always reported in full.
    """
    order = {ERROR: 0, WARN: 1, INFO: 2}
    findings = [f for f in report.findings if show_info or f.severity != INFO]
    shown: dict[str, int] = {}
    for finding in sorted(findings, key=lambda f: (order[f.severity], f.code,
                                                   f.subject)):
        seen = shown.get(finding.code, 0)
        shown[finding.code] = seen + 1
        if seen < per_code:
            print(finding.render())
    for code, total in sorted(shown.items()):
        if total > per_code:
            print(f"[{'':5}] {code}  … and {total - per_code:,} more "
                  f"({total:,} total)")
    if findings:
        print()
    for note in report.notes:
        print(note)
    errors, warns, infos = (report.count(s) for s in (ERROR, WARN, INFO))
    verdict = "FAIL" if errors else "PASS"
    summary = f"{errors} error(s), {warns} warning(s)"
    if infos:
        summary += f", {infos} informational"
    print(f"{verdict} ({summary})")


def run_verify(args) -> int:
    """Execute the ``verify`` subcommand. Never mutates anything."""
    targets = {t for t in ("cache", "store", "remote", "manifest")
               if getattr(args, t, False)}
    if args.all:
        targets = {"cache", "store", "remote", "manifest"}
    if not targets:
        # Offline and lock-driven by default. Reaching the network has to be
        # asked for.
        targets = {"cache", "store"}
    if args.limit and "store" in targets:
        raise SystemExit(
            "--limit cannot be combined with --store: the store's roots are "
            "digests over the complete set, and a partial root is meaningless"
        )

    lock = build_lock.load_lock(args.lock)
    report = VerifyReport()
    if "cache" in targets:
        report.merge(verify_cache(
            lock, args.cache_dir, hash_files=not args.no_hash,
            limit=args.limit, jobs=args.jobs,
        ))
    if "store" in targets:
        report.merge(verify_store(
            lock, args.store_dir, deep=args.deep,
            known_bad_path=args.known_bad,
        ))
    if "manifest" in targets:
        report.merge(verify_manifest(lock, args.config))
    if "remote" in targets:
        report.merge(verify_remote(
            lock, timeout=args.timeout, limit=args.limit, jobs=args.jobs,
        ))

    render(report)
    return 1 if report.failed else 0


# ---------------------------------------------------------------------- status

def run_status(args) -> int:
    """Describe every relationship between manifest, lock, cache, and store.

    Four legs, each reusing the same implementation ``verify`` and ``sync`` use;
    there is no second copy of "resolve and diff" here.

    **Always exits 0.** ``status`` describes; ``verify`` judges and gates CI.
    Conflating the two produces a command nobody can run casually.
    """
    lock = build_lock.load_lock(args.lock)
    print(f"lock      {args.lock}")
    print(f"          written {lock['build']['timestamp_utc']} "
          f"by gtars {lock['build']['gtars_version']}")
    print(f"          {len(build_lock.lock_files(lock))} input file(s), "
          f"{lock['outputs']['n_collections']} collection(s), "
          f"{lock['outputs']['n_sequences']:,} sequence(s)")

    # Leg 1: manifest <-> lock. Local; just a hash of the config file.
    print(f"\nmanifest  {args.config}")
    if not Path(args.config).exists():
        print("          MISSING")
    elif sha256_file(args.config) == build_lock.lock_sources_toml_sha256(lock):
        print("          unchanged since the lock was written")
    else:
        print("          EDITED since the lock was written")

    # Leg 2: manifest <-> cache. The only leg that needs the network.
    if args.offline:
        print("          (upstream comparison skipped: --offline)")
        sources = None
    else:
        assemblies, seqsets = load_config(args.config)
        sources = resolve_sources(assemblies, seqsets)
        check = build_lock.evaluate_sources_vs_lock(sources, lock)
        print(f"          upstream: {len(sources)} source(s) resolve; "
              f"{len(check.new)} added, {len(check.only_in_lock)} removed, "
              f"{len(check.upstream_changed)} checksum(s) changed")
        plan = build_lock.plan_sync(
            sources, lock, args.cache_dir, seqsets=seqsets, hash_files=False,
        )
        print(f"          sync would: unchanged={len(plan.unchanged)} "
              f"added={len(plan.added)} removed={len(plan.removed)} "
              f"reingest={len(plan.reingest)} missing={len(plan.missing)}")

    # Leg 3: lock <-> cache. Local, no hashing -- status is a glance, not a
    # verdict; `verify --cache` is the hashing one.
    print(f"\ncache     {args.cache_dir}")
    build_files = [
        (record["cache_path"], args.cache_dir / record["cache_path"])
        for record in build_lock.lock_files(lock)
    ]
    present = sum(1 for _, path in build_files
                  if path.exists() and path.stat().st_size > 0)
    print(f"          {present}/{len(build_files)} locked file(s) present "
          "(not hashed; use `verify --cache`)")

    # Leg 4: lock <-> store. L0 only: roots and counts, about seven seconds.
    print(f"\nstore     {args.store_dir}")
    if not Path(args.store_dir).exists():
        print("          MISSING")
    else:
        from gtars.refget import RefgetStore

        store = RefgetStore.open_local(str(args.store_dir))
        store.set_quiet(True)
        l0 = verify_store_l0(lock, store)
        if l0.findings:
            for finding in l0.findings:
                print(f"          {finding.code}: {finding.detail}")
        else:
            print(f"          matches the lock: "
                  f"{lock['outputs']['n_sequences']:,} sequence(s), "
                  f"{lock['outputs']['n_collections']} collection(s), "
                  "roots agree")
    return 0
