"""Per-accession sequence retrieval from NCBI E-utilities ``efetch``.

Some NCBI bulk distributions publish only the *current* version of each
accession. RefSeqGene (``NG_``) is one: its FTP FASTA holds ``NG_007400.1``'s
newest version and nothing else, so every superseded version a client may still
cite -- ClinVar cites many -- is absent from any bulk file. E-utilities serves
every version of every accession, so those are fetched one accession at a time.

An efetch seqset declares a checked-in accession list rather than a URL. Its
source "URL" is the real request endpoint and parameters plus a fragment naming
the list's content digest::

    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text#refseqgene_superseded-1a2b3c4d5e6f7a8b.fna.gz

so everything keyed on a source URL -- the cache path, the lock record, drift
detection, sync's add/remove plan -- follows the list's content: editing the
list is a new source, exactly as a new upstream file would be.

The fetched FASTA is written in list order with a fixed gzip header, so the
same list and the same upstream records always produce the same bytes and the
lock's SHA-256 is meaningful. A requested accession NCBI does not return (a
version that never existed, or a withdrawn record) is left out with a warning
and listed in a ``.missing.txt`` sidecar; a transport failure is retried and
then fails the whole fetch, so a partial result is never cached as complete.

Politeness. NCBI allows 3 requests/s without an API key and 10 with one
(``NCBI_API_KEY``). One process-wide :class:`RateLimiter` paces every request,
however many batches run concurrently. When NCBI answers 429, or its JSON
``API rate limit exceeded`` body, every worker pauses for the server's
``Retry-After`` (or a backoff), and the request interval doubles, recovering
gradually on success.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import random
import re
import sqlite3
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlsplit

logger = logging.getLogger("efetch")

EFETCH_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EFETCH_DBS = ("nuccore", "protein")
TOOL = "gks-refgetstore-builder"
GENE_REFSEQGENE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene"

# NCBI recommends POST and at most a few hundred ids per request.
BATCH_SIZE = 200
DEFAULT_JOBS = 3
REQUESTS_PER_SECOND = 3.0
REQUESTS_PER_SECOND_WITH_KEY = 10.0
MAX_INTERVAL = 10.0

ACCESSION_VERSION_RE = re.compile(r"^[A-Z]{1,2}_?\d+\.\d+$")
_FRAGMENT_RE = re.compile(r"^(?P<name>[A-Za-z0-9_.-]+)-(?P<sha>[0-9a-f]{16})\.fna\.gz$")


# ------------------------------------------------------------------ lists

def accession_sort_key(accession: str) -> tuple:
    base, _, version = accession.rpartition(".")
    prefix = re.match(r"^\D*", base).group(0)
    return (prefix, int(base[len(prefix):] or 0), int(version))


def read_accession_list(path: Path) -> tuple[str, ...]:
    """Read a checked-in accession list: one ``accession.version`` per line.

    ``#`` comments and blank lines are ignored. Entries must be well-formed,
    unique and sorted, so a regenerated list diffs cleanly and its digest does
    not depend on how it was assembled.
    """
    ids: list[str] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if not ACCESSION_VERSION_RE.match(line):
            raise ValueError(f"{path}:{number}: not an accession.version: {line!r}")
        ids.append(line)
    if not ids:
        raise ValueError(f"{path}: accession list is empty")
    if len(set(ids)) != len(ids):
        duplicate = next(a for a in ids if ids.count(a) > 1)
        raise ValueError(f"{path}: duplicate accession {duplicate!r}")
    if ids != sorted(ids, key=accession_sort_key):
        raise ValueError(f"{path}: accessions are not sorted")
    return tuple(ids)


def accession_list_sha256(ids: Iterable[str]) -> str:
    return hashlib.sha256("".join(f"{a}\n" for a in ids).encode()).hexdigest()


def efetch_url(db: str, name: str, ids: Iterable[str]) -> str:
    """Source identity of an efetch seqset: request + list digest (see module doc)."""
    query = urllib.parse.urlencode({"db": db, "rettype": "fasta", "retmode": "text"})
    return f"{EFETCH_ENDPOINT}?{query}#{name}-{accession_list_sha256(ids)[:16]}.fna.gz"


def is_efetch_url(url: str) -> bool:
    parts = urlsplit(url)
    return (f"{parts.scheme}://{parts.netloc}{parts.path}" == EFETCH_ENDPOINT
            and bool(_FRAGMENT_RE.match(parts.fragment)))


def efetch_db_of(url: str) -> str:
    return urllib.parse.parse_qs(urlsplit(url).query)["db"][0]


# ----------------------------------------------------------- rate limiting

class RateLimiter:
    """Pace request starts across threads, and slow down when told to.

    ``acquire`` blocks until the next request may start. ``throttled`` is
    called when the server says to slow down: it pauses *every* caller until
    ``retry_after`` has passed and doubles the interval. ``succeeded`` decays
    the interval back toward the base rate.
    """

    def __init__(self, per_second: float, *, clock: Callable[[], float] = time.monotonic,
                 sleep: Callable[[float], None] = time.sleep) -> None:
        self.base_interval = 1.0 / per_second
        self.interval = self.base_interval
        self._next_start = 0.0
        self._lock = threading.Lock()
        self._clock = clock
        self._sleep = sleep

    def acquire(self) -> None:
        with self._lock:
            now = self._clock()
            start = max(now, self._next_start)
            self._next_start = start + self.interval
        if start > now:
            self._sleep(start - now)

    def throttled(self, retry_after: float | None) -> float:
        with self._lock:
            self.interval = min(self.interval * 2, MAX_INTERVAL)
            pause = retry_after if retry_after is not None else self.interval * 2
            self._next_start = max(self._next_start, self._clock() + pause)
            return pause

    def succeeded(self) -> None:
        with self._lock:
            self.interval = max(self.base_interval, self.interval * 0.9)


_DEFAULT_LIMITER: RateLimiter | None = None
_DEFAULT_LIMITER_LOCK = threading.Lock()


def default_limiter() -> RateLimiter:
    """The process-wide limiter: one budget however many fetches are running."""
    global _DEFAULT_LIMITER
    with _DEFAULT_LIMITER_LOCK:
        if _DEFAULT_LIMITER is None:
            rate = (REQUESTS_PER_SECOND_WITH_KEY if os.environ.get("NCBI_API_KEY")
                    else REQUESTS_PER_SECOND)
            _DEFAULT_LIMITER = RateLimiter(rate)
        return _DEFAULT_LIMITER


# ----------------------------------------------------------------- fetching

class RateLimited(Exception):
    def __init__(self, retry_after: float | None) -> None:
        super().__init__(f"rate limited (retry after {retry_after}s)")
        self.retry_after = retry_after


class TransientError(Exception):
    pass


class NotRetrievable(Exception):
    """HTTP 400: NCBI could retrieve none of the requested ids.

    efetch silently omits an unretrievable id requested alongside good ones,
    but rejects a request in which *no* id is retrievable.
    """


@dataclass(frozen=True)
class Response:
    status: int
    headers: dict[str, str]
    body: bytes


def _urllib_post(url: str, data: bytes, timeout: float) -> Response:
    request = urllib.request.Request(url, data=data, headers={"User-Agent": TOOL})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return Response(response.status, dict(response.headers), response.read())
    except urllib.error.HTTPError as exc:
        return Response(exc.code, dict(exc.headers or {}), exc.read() or b"")


def _retry_after(headers: dict[str, str]) -> float | None:
    value = {k.lower(): v for k, v in headers.items()}.get("retry-after")
    try:
        return max(0.0, float(value)) if value is not None else None
    except ValueError:
        return None


def _is_rate_limit_body(body: bytes) -> bool:
    if not body.lstrip().startswith(b"{"):
        return False
    try:
        return "rate limit" in str(json.loads(body).get("error", "")).lower()
    except (ValueError, AttributeError):
        return False


def request_params(db: str, ids: Iterable[str]) -> dict[str, str]:
    params = {"db": db, "rettype": "fasta", "retmode": "text", "tool": TOOL,
              "id": ",".join(ids)}
    # NCBI asks for a contact address; send one only when the operator gave it.
    if email := os.environ.get("NCBI_EMAIL"):
        params["email"] = email
    if api_key := os.environ.get("NCBI_API_KEY"):
        params["api_key"] = api_key
    return params


def fetch_batch(
    db: str, ids: list[str], *, limiter: RateLimiter, post=_urllib_post,
    timeout: float = 300, retries: int = 6, sleep: Callable[[float], None] = time.sleep,
) -> bytes:
    """POST one efetch batch, pacing and retrying politely. Returns FASTA bytes."""
    data = urllib.parse.urlencode(request_params(db, ids)).encode()
    last: Exception | None = None
    for attempt in range(1, retries + 1):
        limiter.acquire()
        try:
            response = post(EFETCH_ENDPOINT, data, timeout)
            if response.status == 429 or _is_rate_limit_body(response.body):
                raise RateLimited(_retry_after(response.headers))
            if response.status >= 500:
                raise TransientError(f"HTTP {response.status}")
            if response.status == 400:
                limiter.succeeded()
                raise NotRetrievable(urllib.parse.unquote_plus(
                    response.body[:300].decode("utf-8", "replace")).strip())
            if response.status != 200:
                raise RuntimeError(
                    f"efetch HTTP {response.status}: {response.body[:200]!r}")
            body = response.body
            if body.strip() and not body.lstrip().startswith(b">"):
                raise RuntimeError(f"efetch returned non-FASTA: {body[:200]!r}")
            limiter.succeeded()
            return body
        except RateLimited as exc:
            last = exc
            pause = limiter.throttled(exc.retry_after)
            logger.warning("NCBI asked to slow down; pausing all requests %.1fs "
                           "(interval now %.2fs)", pause, limiter.interval)
        except (TransientError, urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            delay = min(60.0, 2.0 ** attempt) * (0.5 + random.random() / 2)
            logger.warning("efetch batch failed (%s); retry %d/%d in %.1fs",
                           exc, attempt, retries - 1, delay)
            sleep(delay)
    raise RuntimeError(f"efetch batch of {len(ids)} failed after {retries} attempts: {last}")


def split_records(fasta: bytes) -> dict[str, bytes]:
    """Map each record's header token to its record bytes (blank lines dropped)."""
    records: dict[str, bytes] = {}
    name = None
    lines: list[bytes] = []
    for line in fasta.splitlines():
        if not line.strip():
            continue
        if line.startswith(b">"):
            if name is not None:
                records[name] = b"\n".join(lines) + b"\n"
            name = line[1:].split()[0].decode()
            if name in records:
                raise ValueError(f"efetch returned {name} twice")
            lines = [line]
        elif name is None:
            raise ValueError("efetch output has sequence before any header")
        else:
            lines.append(line)
    if name is not None:
        records[name] = b"\n".join(lines) + b"\n"
    return records


@dataclass
class FetchResult:
    written: int = 0
    bytes: int = 0
    missing: list[str] = field(default_factory=list)


def missing_list_path(target: Path) -> Path:
    return target.parent / (target.name + ".missing.txt")


def materialize(
    db: str, ids: Iterable[str], target: Path, *, jobs: int = DEFAULT_JOBS,
    batch_size: int = BATCH_SIZE, limiter: RateLimiter | None = None,
    fetch: Callable[[list[str]], bytes] | None = None,
) -> FetchResult:
    """Fetch ``ids`` from ``db`` into a gzipped FASTA at ``target``, atomically.

    Batches run on ``jobs`` threads under one shared rate limiter. Records are
    written in ``ids`` order with a fixed gzip header, so output bytes depend
    only on the list and the upstream records.
    """
    ids = list(ids)
    limiter = limiter or default_limiter()
    fetch = fetch or (lambda batch: fetch_batch(db, batch, limiter=limiter))
    batches = [ids[i:i + batch_size] for i in range(0, len(ids), batch_size)]
    done = 0
    progress = threading.Lock()

    def fetch_isolating(batch: list[str]) -> dict[str, bytes]:
        """Fetch ``batch``; on a whole-request rejection, bisect to the culprits."""
        try:
            records = split_records(fetch(batch))
        except NotRetrievable as exc:
            if len(batch) == 1:
                logger.warning("efetch cannot retrieve %s: %s", batch[0], exc)
                return {}
            middle = len(batch) // 2
            return {**fetch_isolating(batch[:middle]), **fetch_isolating(batch[middle:])}
        unexpected = sorted(set(records) - set(batch))
        if unexpected:
            raise ValueError(f"efetch returned unrequested records: {unexpected[:5]}")
        return records

    def run(batch: list[str]) -> dict[str, bytes]:
        nonlocal done
        records = fetch_isolating(batch)
        with progress:
            done += 1
            logger.info("  efetch %s: batch %d/%d (%d records)",
                        target.name, done, len(batches), len(records))
        return records

    records: dict[str, bytes] = {}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for batch_records in pool.map(run, batches):
            records.update(batch_records)

    result = FetchResult(missing=[a for a in ids if a not in records])
    if ids and not records:
        raise RuntimeError(f"efetch retrieved none of the {len(ids)} accessions")
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + ".part")
    # Every gzip parameter is fixed: the lock pins these bytes, and a later
    # fetch of the same list must reproduce them.
    with open(tmp, "wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9
    ) as out:
        for accession in ids:
            if accession in records:
                out.write(records[accession])
                result.written += 1
    missing_path = missing_list_path(target)
    if result.missing:
        logger.warning("efetch returned no record for %d of %d accessions "
                       "(listed in %s): %s", len(result.missing), len(ids),
                       missing_path.name, ", ".join(result.missing[:10]))
        missing_path.write_text("".join(f"{a}\n" for a in result.missing))
    else:
        missing_path.unlink(missing_ok=True)
    tmp.replace(target)
    result.bytes = target.stat().st_size
    return result


# -------------------------------------------------- RefSeqGene version history

def parse_gene_refseqgene(text: str) -> set[str]:
    """Current ``NG_`` accession.versions from NCBI's ``gene_RefSeqGene`` table."""
    lines = [line for line in text.splitlines() if line.strip()]
    header = lines[0].lstrip("#").split("\t")
    try:
        column = header.index("RSG")
    except ValueError as exc:
        raise ValueError(f"gene_RefSeqGene has no RSG column: {header}") from exc
    current = {line.split("\t")[column] for line in lines[1:]}
    bad = sorted(a for a in current if not ACCESSION_VERSION_RE.match(a))
    if bad:
        raise ValueError(f"gene_RefSeqGene has malformed accessions: {bad[:5]}")
    return current


def all_versions(accessions: Iterable[str]) -> set[str]:
    """Every version up to each given one: RefSeq version numbers are dense."""
    out: set[str] = set()
    for accession in accessions:
        base, _, version = accession.rpartition(".")
        out.update(f"{base}.{v}" for v in range(1, int(version) + 1))
    return out


def superseded_history(current: set[str], known: Iterable[str]) -> list[str]:
    """Versions to fetch: all versions of anything known, minus current ones.

    ``known`` carries accessions no longer in ``current`` -- the previous list
    and any seeds -- so a record NCBI retires keeps its history in the list.
    """
    return sorted(all_versions(current | set(known)) - current, key=accession_sort_key)


def seqrepo_ncbi_aliases(aliases_db: Path, prefix: str) -> set[str]:
    con = sqlite3.connect(f"file:{aliases_db}?mode=ro", uri=True)
    try:
        rows = con.execute(
            "select alias from seqalias where namespace = 'NCBI' "
            "and alias like ? escape '\\'", (prefix.replace("_", "\\_") + "%",))
        return {alias for (alias,) in rows if ACCESSION_VERSION_RE.match(alias)}
    finally:
        con.close()


LIST_HEADER = """\
# Superseded and retired RefSeqGene (NG_) accession versions, fetched one
# accession at a time from NCBI E-utilities efetch. NCBI's RefSeqGene FTP holds
# only current versions; this list is every earlier version of every NG_
# accession known here, plus versions of accessions since retired from
# gene_RefSeqGene. Regenerate with `gks-refgetstore refseqgene-history`; the
# command only ever adds entries.
"""


def run_refseqgene_history(args) -> int:
    """Execute the ``refseqgene-history`` subcommand."""
    from .sources import _fetch_text

    current = parse_gene_refseqgene(_fetch_text(args.gene_refseqgene_url))
    known: set[str] = set()
    previous: tuple[str, ...] = ()
    if args.out.exists():
        previous = read_accession_list(args.out)
        known.update(previous)
    for seed in args.seed or []:
        known.update(read_accession_list(seed))
    if args.seqrepo_aliases:
        known.update(seqrepo_ncbi_aliases(args.seqrepo_aliases, "NG_"))
    history = superseded_history(current, known)
    dropped = sorted(set(previous) - set(history), key=accession_sort_key)
    if dropped:
        # A listed version that is now current would mean NCBI rolled a
        # version back, which RefSeq does not do: stop rather than drop it.
        print(f"refusing to drop {len(dropped)} listed accession(s) now current: "
              f"{', '.join(dropped[:10])}", file=sys.stderr)
        return 1
    body = LIST_HEADER + "".join(f"{a}\n" for a in history)
    if args.dry_run:
        sys.stdout.write(body)
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(body, encoding="utf-8")
    print(f"current NG_ versions: {len(current)}; listed: {len(previous)} -> "
          f"{len(history)} (+{len(set(history) - set(previous))})", file=sys.stderr)
    return 0
