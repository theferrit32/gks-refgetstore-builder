"""efetch seqsets: accession lists, source identity, polite fetching, history."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from conftest import file_record, v4_lock

from gks_refgetstore import efetch
from gks_refgetstore.build_lock import apply_locked_sources
from gks_refgetstore.sources import (
    SeqsetConfig,
    load_config,
    mirror_cache_path,
    resolve_sources,
)

IDS = ("NG_000002.1", "NG_000007.1", "NG_000007.2", "NG_007400.1")


def record(accession: str, seq: str = "ACGTN") -> bytes:
    return f">{accession} Homo sapiens test\n{seq}\n".encode()


def fake_ncbi(available: set[str], calls: list[list[str]] | None = None):
    """An efetch stand-in: omits unknown ids, rejects a batch with none known."""
    def fetch(batch: list[str]) -> bytes:
        if calls is not None:
            calls.append(batch)
        found = [a for a in batch if a in available]
        if not found:
            raise efetch.NotRetrievable("Failed to retrieve sequence")
        # NCBI does not promise request order; reverse to prove we don't rely on it.
        return b"\n".join(record(a) for a in reversed(found))
    return fetch


# ------------------------------------------------------------------ lists

def write_list(path: Path, ids, header: str = "# comment\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + "".join(f"{a}\n" for a in ids))
    return path


def test_accession_list_ignores_comments_and_requires_sorted_unique_versions(
    tmp_path: Path,
) -> None:
    assert efetch.read_accession_list(write_list(tmp_path / "ok.txt", IDS)) == IDS
    for bad, match in [
        (["NG_000007.2", "NG_000007.1"], "not sorted"),
        (["NG_000007.1", "NG_000007.1"], "duplicate"),
        (["NG_000007"], "not an accession.version"),
        ([], "empty"),
    ]:
        with pytest.raises(ValueError, match=match):
            efetch.read_accession_list(write_list(tmp_path / "bad.txt", bad))


def test_versions_sort_numerically_not_lexically() -> None:
    ids = ["NG_000007.10", "NG_10.1", "NG_000007.9", "NG_9.1"]
    assert sorted(ids, key=efetch.accession_sort_key) == [
        "NG_000007.9", "NG_000007.10", "NG_9.1", "NG_10.1",
    ]


# --------------------------------------------------------- source identity

def test_source_url_and_cache_path_follow_the_list_content(tmp_path: Path) -> None:
    url = efetch.efetch_url("nuccore", "ng", IDS)
    assert efetch.is_efetch_url(url)
    assert efetch.efetch_db_of(url) == "nuccore"
    assert url.startswith(efetch.EFETCH_ENDPOINT + "?db=nuccore")
    path = mirror_cache_path(tmp_path, url)
    assert path.parent == tmp_path / "eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    assert path.name == f"ng-{efetch.accession_list_sha256(IDS)[:16]}.fna.gz"
    # Any edit to the list is a different source, and a different cache file.
    other = efetch.efetch_url("nuccore", "ng", IDS[:-1])
    assert other != url and mirror_cache_path(tmp_path, other) != path
    # Ordinary URLs keep their mirrored path; a fragment alone does not count.
    plain = "https://ftp.ncbi.nlm.nih.gov/a/b.fna.gz#frag"
    assert not efetch.is_efetch_url(plain)
    assert mirror_cache_path(tmp_path, plain) == tmp_path / "ftp.ncbi.nlm.nih.gov/a/b.fna.gz"


def test_efetch_mode_is_exclusive_and_validated() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        SeqsetConfig("bad", "refseq", urls=["https://x/a"], efetch_accessions="l.txt")
    with pytest.raises(ValueError, match="efetch_db must be one of"):
        SeqsetConfig("bad", "refseq", efetch_accessions="l.txt", efetch_db="gene")
    with pytest.raises(ValueError, match="requires efetch_accessions"):
        SeqsetConfig("bad", "refseq", urls=["https://x/a"], efetch_db="nuccore")
    with pytest.raises(ValueError, match="one FASTA"):
        SeqsetConfig("bad", "refseq", efetch_accessions="l.txt", format="gbff")
    assert SeqsetConfig("ok", "refseq", efetch_accessions="l.txt").efetch_db == "nuccore"
    with pytest.raises(RuntimeError, match="not been loaded"):
        list(SeqsetConfig("ok", "refseq", efetch_accessions="l.txt").iter_shard_urls())


def test_manifest_list_resolves_relative_to_the_manifest(tmp_path: Path) -> None:
    write_list(tmp_path / "cfg/lists/ng.txt", IDS)
    manifest = tmp_path / "cfg/sources.toml"
    manifest.write_text(
        '[[seqset]]\nname = "ng"\nnamespace = "refseq"\n'
        'efetch_accessions = "lists/ng.txt"\nfile_class = "genomic"\n'
    )
    _, (seqset,) = load_config(manifest)
    assert seqset.efetch_ids == IDS
    (source,) = resolve_sources([], [seqset], fetch_text=lambda url: pytest.fail(url))
    assert source.url == efetch.efetch_url("nuccore", "ng", IDS)
    assert source.efetch_ids == IDS
    assert (source.kind, source.owner, source.file_class) == ("seqset", "ng", "genomic")


def test_locked_sources_reattach_the_list_only_while_it_is_unchanged(
    tmp_path: Path,
) -> None:
    write_list(tmp_path / "ng.txt", IDS)
    manifest = tmp_path / "sources.toml"
    manifest.write_text(
        '[[seqset]]\nname = "ng"\nnamespace = "refseq"\nefetch_accessions = "ng.txt"\n'
    )
    _, (seqset,) = load_config(manifest)
    url = efetch.efetch_url("nuccore", "ng", IDS)
    lock = v4_lock([file_record("x.fna.gz", owner="ng", url=url)])
    (source,) = apply_locked_sources([seqset], lock)
    assert source.efetch_ids == IDS

    stale = efetch.efetch_url("nuccore", "ng", IDS[:2])
    lock = v4_lock([file_record("y.fna.gz", owner="ng", url=stale)])
    (source,) = apply_locked_sources([seqset], lock)
    assert source.efetch_ids is None


# ------------------------------------------------------------- materialize

def test_materialize_writes_list_order_deterministically(tmp_path: Path) -> None:
    available = set(IDS)
    one = tmp_path / "one.fna.gz"
    two = tmp_path / "two.fna.gz"
    efetch.materialize("nuccore", IDS, one, batch_size=1, jobs=4,
                       fetch=fake_ncbi(available))
    efetch.materialize("nuccore", IDS, two, batch_size=3, jobs=1,
                       fetch=fake_ncbi(available))
    assert one.read_bytes() == two.read_bytes()
    headers = [line.split()[0][1:].decode()
               for line in gzip.decompress(one.read_bytes()).splitlines()
               if line.startswith(b">")]
    assert headers == list(IDS)
    assert not efetch.missing_list_path(one).exists()


def test_unretrievable_ids_are_isolated_and_listed(tmp_path: Path) -> None:
    ids = ["NG_1.1", "NG_1.2", "NG_1.3", "NG_1.4"]
    calls: list[list[str]] = []
    target = tmp_path / "x.fna.gz"
    result = efetch.materialize("nuccore", ids, target, batch_size=4, jobs=1,
                                fetch=fake_ncbi({"NG_1.1"}, calls))
    assert result.written == 1
    assert result.missing == ["NG_1.2", "NG_1.3", "NG_1.4"]
    assert efetch.missing_list_path(target).read_text() == "NG_1.2\nNG_1.3\nNG_1.4\n"
    # The silently-omitting first request answered; nothing needed bisecting.
    assert calls == [ids]

    # A whole batch NCBI rejects is bisected down to the unretrievable ids.
    calls.clear()
    efetch.materialize("nuccore", ids, target, batch_size=2, jobs=1,
                       fetch=fake_ncbi({"NG_1.1"}, calls))
    assert calls == [["NG_1.1", "NG_1.2"], ["NG_1.3", "NG_1.4"], ["NG_1.3"], ["NG_1.4"]]

    # Fetched again with everything available, the sidecar goes away.
    efetch.materialize("nuccore", ids, target, fetch=fake_ncbi(set(ids)))
    assert not efetch.missing_list_path(target).exists()


def test_materialize_rejects_unrequested_records_and_total_failure(
    tmp_path: Path,
) -> None:
    target = tmp_path / "x.fna.gz"
    with pytest.raises(ValueError, match="unrequested"):
        efetch.materialize("nuccore", ["NG_1.1"], target,
                           fetch=lambda batch: record("NG_1.2"))
    with pytest.raises(RuntimeError, match="none of the 2"):
        efetch.materialize("nuccore", ["NG_1.1", "NG_1.2"], target,
                           fetch=fake_ncbi(set()))
    assert not target.exists()


# --------------------------------------------------------------- politeness

class FakeClock:
    def __init__(self) -> None:
        self.now = 0.0
        self.slept: list[float] = []

    def __call__(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.slept.append(seconds)
        self.now += seconds


def test_rate_limiter_paces_starts_and_backs_off_when_throttled() -> None:
    clock = FakeClock()
    limiter = efetch.RateLimiter(4.0, clock=clock, sleep=clock.sleep)
    for _ in range(3):
        limiter.acquire()
    assert clock.now == pytest.approx(0.5)  # starts at 0, .25, .5

    pause = limiter.throttled(retry_after=3.0)
    assert pause == 3.0
    assert limiter.interval == pytest.approx(0.5)
    limiter.acquire()  # every caller waits out the server's Retry-After
    assert clock.now == pytest.approx(3.5)

    for _ in range(50):
        limiter.succeeded()
    assert limiter.interval == pytest.approx(limiter.base_interval)
    for _ in range(20):
        limiter.throttled(None)
    assert limiter.interval == efetch.MAX_INTERVAL


def scripted_post(*responses: efetch.Response | Exception):
    calls: list[dict] = []
    queue = list(responses)

    def post(url: str, data: bytes, timeout: float) -> efetch.Response:
        from urllib.parse import parse_qs
        calls.append({k: v[0] for k, v in parse_qs(data.decode()).items()})
        item = queue.pop(0)
        if isinstance(item, Exception):
            raise item
        return item
    post.calls = calls  # type: ignore[attr-defined]
    return post


class RecordingLimiter:
    def __init__(self) -> None:
        self.events: list[tuple[str, float | None]] = []
        self.interval = 1.0

    def acquire(self) -> None:
        self.events.append(("acquire", None))

    def throttled(self, retry_after: float | None) -> float:
        self.events.append(("throttled", retry_after))
        return retry_after or 1.0

    def succeeded(self) -> None:
        self.events.append(("succeeded", None))


def test_fetch_batch_honors_rate_limit_responses(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.setenv("NCBI_EMAIL", "ops@example.test")
    ok = efetch.Response(200, {}, record("NG_1.1"))
    post = scripted_post(
        efetch.Response(429, {"Retry-After": "7"}, b""),
        efetch.Response(200, {}, json.dumps(
            {"error": "API rate limit exceeded", "count": "4", "limit": "3"}).encode()),
        ok,
    )
    limiter = RecordingLimiter()
    sleeps: list[float] = []
    body = efetch.fetch_batch("nuccore", ["NG_1.1"], limiter=limiter, post=post,
                              sleep=sleeps.append)
    assert body == ok.body
    assert limiter.events == [
        ("acquire", None), ("throttled", 7.0),
        ("acquire", None), ("throttled", None),
        ("acquire", None), ("succeeded", None),
    ]
    assert sleeps == []  # the limiter, not a private sleep, enforces the pause
    params = post.calls[0]
    assert params["tool"] == efetch.TOOL and params["email"] == "ops@example.test"
    assert params["id"] == "NG_1.1" and "api_key" not in params


def test_fetch_batch_retries_server_errors_then_gives_up() -> None:
    sleeps: list[float] = []
    post = scripted_post(efetch.Response(503, {}, b""), OSError("reset"),
                         efetch.Response(200, {}, record("NG_1.1")))
    assert efetch.fetch_batch("nuccore", ["NG_1.1"], limiter=RecordingLimiter(),
                              post=post, sleep=sleeps.append) == record("NG_1.1")
    assert len(sleeps) == 2 and sleeps[1] > sleeps[0] / 2

    post = scripted_post(*[efetch.Response(502, {}, b"")] * 3)
    with pytest.raises(RuntimeError, match="failed after 3 attempts"):
        efetch.fetch_batch("nuccore", ["NG_1.1"], limiter=RecordingLimiter(),
                           post=post, retries=3, sleep=lambda s: None)


def test_fetch_batch_classifies_rejections_and_non_fasta() -> None:
    post = scripted_post(efetch.Response(
        400, {}, b"+Error%3A+F+a+i+l+e+d++t+o++r+e+t+r+i+e+v+e"))
    with pytest.raises(efetch.NotRetrievable, match="Error"):
        efetch.fetch_batch("nuccore", ["NG_1.99"], limiter=RecordingLimiter(), post=post)
    post = scripted_post(efetch.Response(200, {}, b"<html>maintenance</html>"))
    with pytest.raises(RuntimeError, match="non-FASTA"):
        efetch.fetch_batch("nuccore", ["NG_1.1"], limiter=RecordingLimiter(), post=post)


# ------------------------------------------------------------------ history

GENE_REFSEQGENE = (
    "#tax_id\tGeneID\tSymbol\tRSG\n"
    "9606\t1\tA\tNG_000001.3\n"
    "9606\t2\tB\tNG_000001.3\n"   # one NG_ can serve several genes
    "9606\t3\tC\tNG_000005.1\n"
)


def test_history_is_every_earlier_version_plus_retired_accessions() -> None:
    current = efetch.parse_gene_refseqgene(GENE_REFSEQGENE)
    assert current == {"NG_000001.3", "NG_000005.1"}
    # NG_000009 was retired from gene_RefSeqGene at .2; its .1 is implied.
    assert efetch.superseded_history(current, ["NG_000009.2"]) == [
        "NG_000001.1", "NG_000001.2", "NG_000009.1", "NG_000009.2",
    ]
    with pytest.raises(ValueError, match="no RSG column"):
        efetch.parse_gene_refseqgene("#tax_id\tGeneID\n9606\t1\n")


def test_history_command_only_adds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                   capsys: pytest.CaptureFixture) -> None:
    monkeypatch.setattr("gks_refgetstore.sources._fetch_text",
                        lambda url: GENE_REFSEQGENE)
    out = tmp_path / "lists/ng.txt"
    write_list(out, ["NG_000009.1"])  # retired since; must survive regeneration
    seed = write_list(tmp_path / "seed.txt", ["NG_000012.2"])
    args = SimpleNamespace(out=out, gene_refseqgene_url="x", seed=[seed],
                           seqrepo_aliases=None, dry_run=False)
    assert efetch.run_refseqgene_history(args) == 0
    assert efetch.read_accession_list(out) == (
        "NG_000001.1", "NG_000001.2", "NG_000009.1", "NG_000012.1", "NG_000012.2",
    )
    assert out.read_text().startswith("# Superseded")

    # A listed version that NCBI now calls current is refused, not dropped.
    write_list(out, ["NG_000005.1"])
    args.seed = None
    assert efetch.run_refseqgene_history(args) == 1
    assert "refusing to drop" in capsys.readouterr().err


@pytest.mark.parametrize("name, path, prefixes", [
    ("refseq_human_refseqgene_superseded", "accession_lists/refseqgene_superseded.txt", ("NG_",)),
    ("ncbi36_refseq", "accession_lists/ncbi36_refseq.txt", ("NC_", "NG_", "NT_")),
])
def test_checked_in_list_is_valid_and_configured(name: str, path: str, prefixes: tuple[str, ...]) -> None:
    root = Path(__file__).parents[1]
    ids = efetch.read_accession_list(root / path)
    assert all(a.startswith(prefixes) for a in ids)
    digest = hashlib.sha256("".join(f"{a}\n" for a in ids).encode()).hexdigest()
    assert efetch.accession_list_sha256(ids) == digest
    _, seqsets = load_config(root / "sources.toml")
    (seqset,) = [s for s in seqsets if s.name == name]
    assert seqset.efetch_ids == ids and seqset.namespace == "refseq"
