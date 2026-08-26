from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

import pytest

import build_lock
import build_store
from build_store import (ResolvedSource, SeqsetConfig, apply_locked_sources,
                         load_config, parse_checksum_manifest,
                         parse_ensembl_checksum_manifest, resolve_sources)


BASE = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/"
MANIFEST = BASE + "human.files.installed"


def pattern(name: str, suffix: str) -> SeqsetConfig:
    return SeqsetConfig(name, "refseq", url_pattern=BASE + "human.*." + suffix,
                        checksum_manifest_url=MANIFEST)


def test_config_source_modes_are_exclusive_and_pattern_pair_is_required() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        SeqsetConfig("bad", "x", url_template="https://x/a", urls=["https://x/b"])
    with pytest.raises(ValueError, match="must be set together"):
        SeqsetConfig("bad", "x", url_pattern="https://x/a/*.gz")
    with pytest.raises(ValueError, match="only valid"):
        SeqsetConfig("bad", "x", urls=["https://x/a"], shard_range=[1, 2])
    with pytest.raises(ValueError, match="require HTTPS"):
        SeqsetConfig("bad", "x", url_pattern="ftp://x/a/*.gz",
                     checksum_manifest_url="ftp://x/a/files.installed")
    with pytest.raises(ValueError, match="format must be one of"):
        SeqsetConfig("bad", "x", url_template="https://x/a", format="zip")
    assert SeqsetConfig(
        "ok", "lrg", urls=["https://x/a.zip"], format="lrg_zip"
    ).format == "lrg_zip"


def test_derived_fasta_resolvers_cover_every_non_fasta_format() -> None:
    assert set(build_store.DERIVED_FASTA_RESOLVERS) | {"fasta"} == set(
        build_store.SEQSET_FORMATS
    )


@pytest.mark.parametrize("text", [
    "not-a-record", "0" * 32 + "  dir/file.gz", "0" * 32 + "  ../file.gz",
    "0" * 32 + "  file.gz\n" + "1" * 32 + "  file.gz",
    "0" * 32 + " *file.gz",
])
def test_manifest_rejects_malformed_duplicate_or_unsafe_records(text: str) -> None:
    with pytest.raises(ValueError):
        parse_checksum_manifest(text)


def test_ensembl_checksum_parser_and_bsd_sum(tmp_path: Path) -> None:
    assert parse_ensembl_checksum_manifest("65 1 one.fa.gz\n") == {
        "one.fa.gz": ("00065", 1)
    }
    path = tmp_path / "one"
    path.write_bytes(b"A")
    assert build_store.bsd_sum_file(path) == ("00065", 1)
    with pytest.raises(ValueError, match="unsafe"):
        parse_ensembl_checksum_manifest("1 1 ../one.fa.gz\n")


def test_shared_manifest_filtering_natural_order_gaps_and_changing_counts() -> None:
    calls: list[str] = []
    text = "\n".join([
        "a" * 32 + "  human.10.rna.fna.gz",
        "b" * 32 + "  human.2.rna.fna.gz",
        "c" * 32 + "  human.7.protein.faa.gz",
        "d" * 32 + "  README",
    ])

    def fetch(url: str) -> str:
        calls.append(url)
        return text

    seqsets = [pattern("rna", "rna.fna.gz"), pattern("protein", "protein.faa.gz")]
    sources = resolve_sources([], seqsets, fetch)
    assert calls == [MANIFEST]
    assert [Path(source.url).name for source in sources] == [
        "human.2.rna.fna.gz", "human.10.rna.fna.gz", "human.7.protein.faa.gz"
    ]
    assert [source.upstream_md5 for source in sources] == ["b" * 32, "a" * 32, "c" * 32]


def test_pattern_file_class_applies_to_every_resolved_source() -> None:
    seqset = SeqsetConfig(
        "rna", "refseq", url_pattern=BASE + "human.*.rna.fna.gz",
        checksum_manifest_url=MANIFEST, file_class="rna",
    )
    manifest = "\n".join(
        f"{value * 32}  human.{index}.rna.fna.gz"
        for index, value in ((1, "a"), (2, "b"))
    )
    sources = resolve_sources([], [seqset], lambda _: manifest)
    assert [source.file_class for source in sources] == ["rna", "rna"]
    with pytest.raises(ValueError, match="only one file class mode"):
        SeqsetConfig(
            "bad", "refseq", urls=["https://x/a"], file_class="rna",
            file_classes=["rna"],
        )


def test_pattern_empty_match_is_an_error() -> None:
    with pytest.raises(ValueError, match="matched no"):
        resolve_sources([], [pattern("rna", "rna.fna.gz")],
                        lambda _: "a" * 32 + "  other.gz")


def test_refseqgene_current_eight_file_result() -> None:
    base = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/"
    seqset = SeqsetConfig("genes", "refseq", url_pattern=base + "refseqgene.*.genomic.fna.gz",
                          checksum_manifest_url=base + "refseqgene.files.installed")
    manifest = "\n".join(f"{i:032x}  refseqgene.{i}.genomic.fna.gz" for i in range(1, 9))
    assert len(resolve_sources([], [seqset], lambda _: manifest)) == 8


def test_lock_classifies_membership_and_same_url_md5_change() -> None:
    live = [ResolvedSource("seqset", "rna", BASE + "human.1.rna.fna.gz", "b" * 32),
            ResolvedSource("seqset", "rna", BASE + "human.2.rna.fna.gz", "c" * 32)]
    lock = {"schema": build_lock.SCHEMA, "sources": [
        {"kind": "seqset", "owner": "rna", "url": live[0].url, "upstream_md5": "a" * 32},
        {"kind": "seqset", "owner": "rna", "url": BASE + "human.9.rna.fna.gz",
         "upstream_md5": "d" * 32},
    ]}
    result = build_lock.evaluate_sources_vs_lock(live, lock)
    assert result.new == [live[1].url]
    assert result.only_in_lock == [BASE + "human.9.rna.fna.gz"]
    assert result.upstream_changed == [(live[0].url, "a" * 32, "b" * 32)]
    assert result.drift


def test_lock_records_ncbi_md5_as_generic_provider_checksum(tmp_path: Path) -> None:
    source = ResolvedSource(
        "seqset", "rna", BASE + "human.1.rna.fna.gz", "a" * 32,
        checksum_url=MANIFEST,
    )
    target = build_store.mirror_cache_path(tmp_path, source.url)
    target.parent.mkdir(parents=True)
    target.write_bytes(b"cached")
    record = build_lock.source_record(source, tmp_path, {}, hash_files=False)
    assert record["upstream_md5"] == "a" * 32
    assert record["provider_checksum"] == "a" * 32
    assert record["provider_checksum_algorithm"] == "md5"


def test_locked_sources_uses_v2_concrete_urls_without_resolution(tmp_path: Path) -> None:
    seqset = pattern("rna", "rna.fna.gz")
    locked = ResolvedSource("seqset", "rna", BASE + "human.3.rna.fna.gz", "a" * 32)
    lock = {"schema": build_lock.SCHEMA, "sources": [locked.__dict__]}
    assert apply_locked_sources([seqset], lock) == [locked]
    assert list(seqset.iter_shard_urls()) == [("1", locked.url)]
    with pytest.raises(ValueError, match="v2"):
        apply_locked_sources([seqset], {"schema": build_lock.V1_SCHEMA, "sources": []})


def test_cache_sha256_classification(tmp_path: Path) -> None:
    cached = tmp_path / "x.gz"
    cached.write_bytes(b"new")
    rel = "x.gz"
    lock = {"sources": [{"cache_path": rel, "sha256": hashlib.sha256(b"old").hexdigest()}]}
    result = build_lock.evaluate_cache_vs_lock([(rel, cached), ("missing", tmp_path / "missing")], lock)
    assert result.changed[0][0] == rel
    assert result.missing == ["missing"]


def test_seqset_ingestion_reuses_prepared_cache_without_downloading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    url = "https://example.test/one.fa"
    seqset = SeqsetConfig("one", "x", url_template=url)
    target = build_store.mirror_cache_path(tmp_path, url)
    target.parent.mkdir(parents=True)
    target.write_bytes(b">one\nA\n")

    monkeypatch.setattr(
        build_store, "ensure_download",
        lambda *_args, **_kwargs: pytest.fail("ingestion attempted a download"),
    )

    class Meta:
        digest = "collection"
        n_sequences = 1

    class Store:
        def add_sequence_collection_from_fasta(self, path):
            assert Path(path) == target
            return Meta(), True

        def load_collection(self, _digest): pass
        def get_collection_level2(self, _digest):
            return {"names": ["one"], "sequences": ["SQ.digest"]}
        def list_sequence_aliases(self, _namespace): return []
        def load_sequence_aliases(self, _namespace, _path): return 1

    stats = build_store.process_seqset(Store(), seqset, tmp_path)
    assert stats.shards_processed == 1


def test_partial_lock_merge_replaces_complete_touched_scopes(tmp_path: Path) -> None:
    config = tmp_path / "sources.toml"
    config.write_text("")
    cache = tmp_path / "cache"

    old_url = "https://example.test/old.fa"
    changed_url = "https://example.test/changed.fa"
    added_url = "https://example.test/added.fa"
    keep_url = "https://example.test/keep.fa"
    for url, content in ((changed_url, b"new changed"), (added_url, b"added")):
        path = build_store.mirror_cache_path(cache, url)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    def old_record(kind: str, owner: str, url: str) -> dict:
        return {
            "kind": kind, "owner": owner, "url": url,
            "cache_path": str(build_store.mirror_cache_path(cache, url)
                              .relative_to(cache)),
            "sha256": "old", "upstream_md5": None,
        }

    existing = {"schema": build_lock.SCHEMA, "sources": [
        old_record("seqset", "rna", old_url),
        old_record("seqset", "rna", changed_url),
        old_record("seqset", "protein", keep_url),
    ]}
    touched = [
        ResolvedSource("seqset", "rna", changed_url),
        ResolvedSource("seqset", "rna", added_url),
    ]

    class Store:
        @staticmethod
        def stats(): return {}

    merged = build_lock.merge_into_lock(
        existing, config_path=config, download_dir=cache,
        touched_sources=touched, collection_by_cachepath={}, store=Store(),
    )
    by_url = {source["url"]: source for source in merged["sources"]}
    assert set(by_url) == {changed_url, added_url, keep_url}
    assert by_url[changed_url]["sha256"] == hashlib.sha256(b"new changed").hexdigest()
    assert by_url[keep_url]["sha256"] == "old"
    assert not build_lock.evaluate_sources_vs_lock(touched, merged).drift


def test_locked_sources_with_force_download_remains_offline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    url = "https://example.test/one.fa"
    config = tmp_path / "sources.toml"
    config.write_text('[[seqset]]\nname="one"\nnamespace="x"\n'
                      f'url_template="{url}"\n')
    cache = tmp_path / "cache"
    target = build_store.mirror_cache_path(cache, url)
    target.parent.mkdir(parents=True)
    target.write_bytes(b">one\nA\n")
    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, {
        "schema": build_lock.SCHEMA,
        "sources": [{
            "kind": "seqset", "owner": "one", "url": url,
            "cache_path": str(target.relative_to(cache)),
            "sha256": build_lock.sha256_file(target), "upstream_md5": None,
        }],
    })

    monkeypatch.setattr(
        build_store, "ensure_download",
        lambda *_args, **_kwargs: pytest.fail("locked build accessed the network"),
    )
    seen: list[bool] = []

    def ingest(_store, entry, _cache, _provenance, refresh_derived=False,
               alias_sink=None, jobs=1):
        seen.append(refresh_derived)
        stats = build_store.SeqsetStats(entry.name, entry.namespace)
        stats.shards_processed = 1
        return stats

    monkeypatch.setattr(build_store, "process_seqset", ingest)

    class Store:
        storage_mode = "test"
        @staticmethod
        def on_disk(_path): return Store()
        @staticmethod
        def write(): pass
        @staticmethod
        def stats(): return {}
        @staticmethod
        def list_sequence_alias_namespaces(): return []
        @staticmethod
        def list_sequence_aliases(_namespace): return []
        @staticmethod
        def list_collection_alias_namespaces(): return []

    monkeypatch.setattr(build_store, "RefgetStore", Store)
    args = type("Args", (), {
        "config": config, "cache_dir": cache, "store_dir": tmp_path / "store",
        "lock": lock_path, "locked_sources": True, "skip_assemblies": False,
        "skip_seqsets": False, "assembly": None, "seqset": None,
        "lock_check_mode": "strict", "force_download": True,
        "force_lock": False, "no_lock": True,
        "ingest_jobs": build_store.INGEST_JOBS_DEFAULT, "min_free_gb": 0.0,
    })()
    assert build_store.run_build(args) == 0
    assert seen == [True]


def test_v1_lock_stops_build_before_store_is_opened(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "sources.toml"
    config.write_text('[[seqset]]\nname="one"\nnamespace="x"\n'
                      'url_template="https://example.test/one.fa"\n')
    cache = tmp_path / "cache"
    target = build_store.mirror_cache_path(cache, "https://example.test/one.fa")
    target.parent.mkdir(parents=True)
    target.write_bytes(b">x\nA\n")
    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, {
        "schema": build_lock.V1_SCHEMA,
        "sources": [{"kind": "seqset", "owner": "one",
                     "url": "https://example.test/one.fa",
                     "cache_path": str(target.relative_to(cache)),
                     "sha256": build_lock.sha256_file(target)}],
    })
    class FailStore:
        @staticmethod
        def on_disk(*_):
            pytest.fail("store opened before drift gate")

    monkeypatch.setattr(build_store, "RefgetStore", FailStore)
    args = type("Args", (), {
        "config": config, "cache_dir": cache, "store_dir": tmp_path / "store",
        "lock": lock_path, "locked_sources": False, "skip_assemblies": False,
        "skip_seqsets": False, "assembly": None, "seqset": None,
        "lock_check_mode": "strict", "force_download": False,
        "force_lock": False, "no_lock": False,
    })()
    with pytest.raises(SystemExit, match="v1"):
        build_store.run_build(args)


def test_checked_in_config_loads_and_resolves_without_duplicates(monkeypatch: pytest.MonkeyPatch) -> None:
    assemblies, seqsets = load_config(Path("sources.toml"))
    manifests = {
        "human.files.installed": "\n".join([
            "a" * 32 + "  human.1.rna.fna.gz",
            "b" * 32 + "  human.1.protein.faa.gz",
        ]),
        "refseqgene.files.installed": "c" * 32 + "  refseqgene.1.genomic.fna.gz",
    }
    md5_files: dict[str, set[str]] = {}
    for assembly in assemblies:
        if assembly.checksum_manifest_url:
            bucket = md5_files.setdefault(assembly.checksum_manifest_url, set())
            bucket.add(Path(assembly.report_url).name)
            if assembly.fasta_url:
                bucket.add(Path(assembly.fasta_url).name)
    for seqset in seqsets:
        if seqset.md5_manifest_urls:
            for url, checksum_url in zip(seqset.urls or [], seqset.md5_manifest_urls):
                md5_files.setdefault(checksum_url, set()).add(Path(url).name)

    def manifest_for(url: str) -> str:
        if Path(url).name == "md5checksums.txt":
            return "\n".join(
                f"{'d' * 32}  ./{name}" for name in sorted(md5_files[url])
            )
        if Path(url).name != "CHECKSUMS":
            return manifests[Path(url).name]
        match = build_store.re.search(r"release-(\d+)", url)
        assert match
        release = int(match.group(1))
        assembly = "GRCh37.75" if release == 75 else "GRCh38"
        directory = Path(build_store.urlsplit(url).path).parent.name
        suffix = {
            "dna": "dna.toplevel.fa.gz", "cdna": "cdna.all.fa.gz",
            "ncrna": "ncrna.fa.gz", "pep": "pep.all.fa.gz",
        }[directory]
        return f"1 1 Homo_sapiens.{assembly}.{suffix}\n"

    sources = resolve_sources(assemblies, seqsets, manifest_for)
    urls = [source.url for source in sources]
    assert len(urls) == len(set(urls))
    assert all("//" not in url.split("://", 1)[1] for url in urls)
    ensembl = [source for source in sources if source.owner.startswith("ensembl_release_")]
    assert len(ensembl) == 42 * 4
    assert {source.file_class for source in ensembl} == {
        "dna.toplevel", "cdna", "ncrna", "pep"
    }


def test_release_groups_are_chronological_atomic_and_replace_rolling_namespace(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = [
        SeqsetConfig("r2", "ensembl-2", urls=["https://x/r2.fa"],
                     file_classes=["x"], release=2, rolling_namespace="ensembl"),
        SeqsetConfig("r1", "ensembl-1", urls=["https://x/r1.fa"],
                     file_classes=["x"], release=1, rolling_namespace="ensembl"),
    ]
    seen: list[str] = []

    def ingest(_store, entry, _cache, _provenance=None, _refresh=False,
               alias_sink=None, jobs=1):
        seen.append(entry.name)
        alias_sink.update({"shared": f"digest-{entry.release}"})
        if entry.release == 1:
            alias_sink["retired"] = "old"
        stats = build_store.SeqsetStats(entry.name, entry.namespace)
        stats.shards_processed = 1
        return stats

    monkeypatch.setattr(build_store, "process_seqset", ingest)

    class Store:
        def __init__(self): self.aliases = {}
        def list_sequence_aliases(self, namespace):
            return list(self.aliases.get(namespace, {}))
        def load_sequence_aliases(self, namespace, path):
            rows = dict(line.rstrip().split("\t") for line in Path(path).read_text().splitlines())
            self.aliases.setdefault(namespace, {}).update(rows)
            return len(rows)

    store = Store()
    build_store.process_release_groups(store, entries, tmp_path, tmp_path / "store")
    assert seen == ["r1", "r2"]
    assert store.aliases["ensembl-1"] == {"retired": "old", "shared": "digest-1"}
    assert store.aliases["ensembl-2"] == {"shared": "digest-2"}
    rolling = (tmp_path / "store/aliases/sequences/ensembl.tsv").read_text()
    assert rolling == "shared\tdigest-2\n"


def test_release_group_rejects_conflicting_aliases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = [
        SeqsetConfig("a", "ensembl-1", urls=["https://x/a"], file_classes=["a"],
                     release=1, rolling_namespace="ensembl"),
        SeqsetConfig("b", "ensembl-1", urls=["https://x/b"], file_classes=["b"],
                     release=1, rolling_namespace="ensembl"),
    ]
    def ingest(_store, entry, _cache, _provenance=None, _refresh=False,
               alias_sink=None, jobs=1):
        digest = "one" if entry.name == "a" else "two"
        previous = alias_sink.get("same")
        if previous is not None and previous != digest:
            raise ValueError("release 1: alias 'same' maps to both digests")
        alias_sink["same"] = digest
        stats = build_store.SeqsetStats(entry.name, entry.namespace)
        stats.shards_processed = 1
        return stats
    monkeypatch.setattr(build_store, "process_seqset", ingest)
    with pytest.raises(ValueError, match="maps to both"):
        build_store.process_release_groups(object(), entries, tmp_path, tmp_path / "store")


def test_release_group_does_not_publish_partial_release(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entry = SeqsetConfig(
        "r1", "ensembl-1", urls=["https://x/a"], file_classes=["x"],
        release=1, rolling_namespace="ensembl",
    )
    def incomplete(*_args, **_kwargs):
        stats = build_store.SeqsetStats(entry.name, entry.namespace)
        stats.warnings = 1
        return stats
    monkeypatch.setattr(build_store, "process_seqset", incomplete)
    class Store:
        @staticmethod
        def list_sequence_alias_namespaces(): return []
    with pytest.raises(RuntimeError, match="aliases not published"):
        build_store.process_release_groups(
            Store(), [entry], tmp_path, tmp_path / "store"
        )
    assert not (tmp_path / "store/aliases/sequences/ensembl.tsv").exists()


def test_load_config_rejects_incomplete_ensembl_release(tmp_path: Path) -> None:
    config = tmp_path / "sources.toml"
    config.write_text(
        '[[seqset]]\nname="r75"\nnamespace="ensembl-75"\nrelease=75\n'
        'rolling_namespace="ensembl"\nurls=["https://x/cdna"]\n'
        'file_classes=["cdna"]\n'
    )
    with pytest.raises(ValueError, match="expected file classes"):
        load_config(config)


def test_immutable_alias_loader_rejects_existing_digest_change() -> None:
    class Store:
        @staticmethod
        def list_sequence_aliases(_namespace): return ["A.1"]
        @staticmethod
        def get_sequence_metadata_by_alias(_namespace, _alias):
            return type("Metadata", (), {"sha512t24u": "old"})()

    with pytest.raises(ValueError, match="immutable alias collision"):
        build_store.load_immutable_aliases(Store(), "refseq", {"A.1": "new"})


def test_real_store_release_aliases_survive_write_and_reopen(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    store_path = tmp_path / "store"
    entries = []
    for release, fasta in (
        (1, b">shared\nA\n>retired\nC\n"),
        (2, b">shared\nG\n"),
    ):
        url = f"https://example.test/r{release}.fa"
        target = build_store.mirror_cache_path(cache, url)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(fasta)
        entries.append(SeqsetConfig(
            f"r{release}", f"ensembl-{release}", urls=[url],
            file_classes=["x"], release=release, rolling_namespace="ensembl",
        ))
    store = build_store.RefgetStore.on_disk(str(store_path))
    build_store.process_release_groups(store, entries, cache, store_path)
    store.write()
    reopened = build_store.RefgetStore.open_local(str(store_path))
    reopened.pull_aliases()
    assert reopened.get_sequence_by_alias("ensembl-1", "retired") is not None
    assert reopened.get_sequence_by_alias("ensembl", "retired") is None
    rolling = reopened.get_sequence_by_alias("ensembl", "shared")
    immutable = reopened.get_sequence_by_alias("ensembl-2", "shared")
    assert rolling is not None and immutable is not None
    assert rolling.metadata.sha512t24u == immutable.metadata.sha512t24u


def _lrg_member(locus: str, extra: str = "", trailing_newline: bool = True) -> bytes:
    text = (
        f">{locus}g (genomic sequence)\nACGT\n"
        f">{locus}t1 (transcript t1 of {locus})\nACG\n"
        f">{locus}p1 (protein translated from transcript t1 of {locus})\nMA\n"
        f"{extra}"
    )
    return text.encode() if trailing_newline else text.rstrip("\n").encode()


def _write_lrg_zip(path: Path, members: dict[str, bytes]) -> Path:
    with zipfile.ZipFile(path, "w") as archive:
        for name, payload in members.items():
            archive.writestr(name, payload)
    return path


def test_lrg_zip_conversion_is_ordered_newline_safe_and_complete(tmp_path: Path) -> None:
    # Deliberately archived in the wrong order, and LRG_2 both lacks a trailing
    # newline and carries a second transcript/protein pair.
    zip_path = _write_lrg_zip(tmp_path / "bundle.zip", {
        "LRG_10.fasta": _lrg_member("LRG_10"),
        "LRG_2.fasta": _lrg_member(
            "LRG_2",
            extra=(
                ">LRG_2t2 (transcript t2 of LRG_2)\nAC\n"
                ">LRG_2p2 (protein translated from transcript t2 of LRG_2)\nM\n"
            ),
            trailing_newline=False,
        ),
        "LRG_1.fasta": _lrg_member("LRG_1"),
    })
    out = tmp_path / "bundle.fasta"

    assert build_store.lrg_zip_to_fasta(zip_path, out) == (3, 11)

    text = out.read_text()
    headers = [line[1:].split()[0] for line in text.splitlines() if line.startswith(">")]
    assert headers == [
        "LRG_1g", "LRG_1t1", "LRG_1p1",
        "LRG_2g", "LRG_2t1", "LRG_2p1", "LRG_2t2", "LRG_2p2",
        "LRG_10g", "LRG_10t1", "LRG_10p1",
    ]
    # The newline guard must keep LRG_2p2's sequence and LRG_10g's header on
    # separate lines; a fused member boundary would leave a mid-line ">".
    assert all(
        ">" not in line[1:] and (line.startswith(">") or ">" not in line)
        for line in text.splitlines()
    )
    assert text.endswith("\n")
    assert not (tmp_path / "bundle.fasta.part").exists()


@pytest.mark.parametrize("members, message", [
    ({}, "no LRG_N.fasta members"),
    ({"README.txt": b"hello\n"}, "unexpected member"),
    ({"LRG_1.fasta": b">NG_007400.1\nACGT\n"}, "genomic LRG header"),
    ({"LRG_1.fasta": b">LRG_1g (genomic sequence)\nACGT\n"}, "expected a genomic"),
    ({"LRG_1.fasta": b""}, "is empty"),
])
def test_lrg_zip_conversion_rejects_malformed_bundles(
    tmp_path: Path, members: dict[str, bytes], message: str
) -> None:
    zip_path = _write_lrg_zip(tmp_path / "bad.zip", members)
    with pytest.raises(ValueError, match=message):
        build_store.lrg_zip_to_fasta(zip_path, tmp_path / "bad.fasta")


@pytest.mark.parametrize("fmt", ["lrg_zip", "gbff"])
def test_derived_fasta_is_cached_until_force(tmp_path: Path, fmt: str) -> None:
    if fmt == "lrg_zip":
        source = _write_lrg_zip(tmp_path / "b.zip", {"LRG_1.fasta": _lrg_member("LRG_1")})
    else:
        source = tmp_path / "b.gbff"
        source.write_text("VERSION     NM_1.1\nORIGIN\n        1 acgt\n//\n")
    resolve = build_store.DERIVED_FASTA_RESOLVERS[fmt]

    out = resolve(source, False)
    assert out == source.parent / (source.name + ".fasta")
    original = out.read_text()

    out.write_text(">sentinel\nA\n")
    assert resolve(source, False).read_text() == ">sentinel\nA\n"
    assert resolve(source, True).read_text() == original


def test_lrg_alias_names_expands_only_the_genomic_record() -> None:
    assert build_store.lrg_alias_names("LRG_1g") == ("LRG_1g", "LRG_1")
    assert build_store.lrg_alias_names("LRG_123g") == ("LRG_123g", "LRG_123")
    assert build_store.lrg_alias_names("LRG_1t1") == ("LRG_1t1",)
    assert build_store.lrg_alias_names("LRG_1p1") == ("LRG_1p1",)
    assert build_store.lrg_alias_names("NG_007400.1") == ("NG_007400.1",)


def test_lrg_seqset_ingest_emits_both_genomic_aliases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    url = "https://ftp.example.test/lrgex/fasta/LRG_public_fasta_files.zip"
    target = build_store.mirror_cache_path(tmp_path, url)
    target.parent.mkdir(parents=True)
    _write_lrg_zip(target, {"LRG_1.fasta": _lrg_member("LRG_1")})
    seqset = SeqsetConfig("lrg_public", "lrg", urls=[url], format="lrg_zip")

    monkeypatch.setattr(
        build_store, "ensure_download",
        lambda *_args, **_kwargs: pytest.fail("ingestion attempted a download"),
    )

    class Meta:
        digest = "collection"
        n_sequences = 3

    class Store:
        def add_sequence_collection_from_fasta(self, path):
            assert Path(path) == target.parent / (target.name + ".fasta")
            return Meta(), True

        def load_collection(self, _digest): pass
        def get_collection_level2(self, _digest):
            return {
                "names": ["LRG_1g", "LRG_1t1", "LRG_1p1"],
                "sequences": ["SQ.g", "SQ.t", "SQ.p"],
            }

    sink: dict[str, str] = {}
    stats = build_store.process_seqset(Store(), seqset, tmp_path, alias_sink=sink)
    assert stats.shards_processed == 1
    # build_name_to_digest_map strips the SQ. prefix.
    assert sink == {"LRG_1g": "g", "LRG_1": "g", "LRG_1t1": "t", "LRG_1p1": "p"}


def test_real_store_ingests_mixed_lrg_nucleotide_and_protein_records(
    tmp_path: Path
) -> None:
    cache = tmp_path / "cache"
    url = "https://ftp.example.test/lrgex/fasta/LRG_public_fasta_files.zip"
    target = build_store.mirror_cache_path(cache, url)
    target.parent.mkdir(parents=True)
    _write_lrg_zip(target, {
        "LRG_1.fasta": _lrg_member("LRG_1"),
        "LRG_2.fasta": _lrg_member("LRG_2"),
    })
    seqset = SeqsetConfig("lrg_public", "lrg", urls=[url], format="lrg_zip")

    store = build_store.RefgetStore.in_memory()
    stats = build_store.process_seqset(store, seqset, cache)
    assert stats.shards_processed == 1
    assert stats.sequences_ingested == 6

    genomic = store.get_sequence_metadata_by_alias("lrg", "LRG_1g")
    assert store.get_sequence_metadata_by_alias(
        "lrg", "LRG_1"
    ).sha512t24u == genomic.sha512t24u
    assert store.get_sequence_metadata_by_alias("lrg", "LRG_2p1") is not None


class BatchStore:
    """Records how ingestion was invoked: batched vs one file at a time."""

    def __init__(self, batch_error: Exception | None = None,
                 short_by: int = 0) -> None:
        self.batch_calls: list[tuple[list[str], int]] = []
        self.single_calls: list[str] = []
        self.batch_error = batch_error
        self.short_by = short_by
        self.aliases: dict[str, str] = {}

    def _meta(self, path: str):
        name = Path(path).stem
        return type("Meta", (), {"digest": f"coll-{name}", "n_sequences": 1})()

    def add_sequence_collections_from_fastas(self, paths, jobs=1):
        self.batch_calls.append((list(paths), jobs))
        if self.batch_error is not None:
            raise self.batch_error
        kept = paths[: len(paths) - self.short_by] if self.short_by else paths
        return [(self._meta(p), True) for p in kept]

    def add_sequence_collection_from_fasta(self, path):
        self.single_calls.append(str(path))
        if str(path).endswith("bad.fa"):
            raise RuntimeError("corrupt shard")
        return self._meta(str(path)), True

    def load_collection(self, _digest): pass

    def get_collection_level2(self, digest):
        return {"names": [f"n-{digest}"], "sequences": [f"SQ.d-{digest}"]}

    def list_sequence_aliases(self, _ns): return []

    def load_sequence_aliases(self, _ns, path):
        rows = Path(path).read_text().splitlines()
        self.aliases.update(dict(r.split("\t") for r in rows if r))
        return len(rows)


def _shard_seqset(tmp_path: Path, names: list[str]) -> SeqsetConfig:
    urls = []
    for name in names:
        url = f"https://example.test/{name}"
        target = build_store.mirror_cache_path(tmp_path, url)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b">x\nA\n")
        urls.append(url)
    return SeqsetConfig("multi", "ns", urls=urls, file_classes=["c"] * len(urls))


def test_seqset_ingests_every_shard_in_one_batched_call(tmp_path: Path) -> None:
    seqset = _shard_seqset(tmp_path, ["a.fa", "b.fa", "c.fa"])
    store = BatchStore()
    provenance: dict[str, dict] = {}

    stats = build_store.process_seqset(
        store, seqset, tmp_path, provenance, jobs=4
    )

    assert len(store.batch_calls) == 1, "expected exactly one batched import"
    paths, jobs = store.batch_calls[0]
    assert [Path(p).name for p in paths] == ["a.fa", "b.fa", "c.fa"]
    assert jobs == 4
    assert store.single_calls == []
    assert stats.shards_processed == 3
    # provenance is keyed by source cache path, one entry per shard, in order
    assert [Path(k).name for k in provenance] == ["a.fa", "b.fa", "c.fa"]
    assert provenance[
        str(build_store.mirror_cache_path(tmp_path, "https://example.test/b.fa"))
    ]["collection_digest"] == "coll-b"


def test_batched_ingest_failure_falls_back_to_serial_and_names_the_file(
    tmp_path: Path,
) -> None:
    seqset = _shard_seqset(tmp_path, ["good.fa", "bad.fa", "other.fa"])
    store = BatchStore(batch_error=RuntimeError("batch blew up"))

    stats = build_store.process_seqset(store, seqset, tmp_path, jobs=4)

    assert len(store.batch_calls) == 1
    assert [Path(p).name for p in store.single_calls] == [
        "good.fa", "bad.fa", "other.fa"
    ]
    # only the one corrupt shard is lost; its siblings still land
    assert stats.warnings == 1
    assert stats.shards_processed == 2


def test_batched_result_length_mismatch_falls_back_rather_than_misaligning(
    tmp_path: Path,
) -> None:
    seqset = _shard_seqset(tmp_path, ["a.fa", "b.fa", "c.fa"])
    store = BatchStore(short_by=1)

    stats = build_store.process_seqset(store, seqset, tmp_path, jobs=2)

    assert len(store.batch_calls) == 1
    assert len(store.single_calls) == 3, "must re-ingest serially, not mis-map"
    assert stats.shards_processed == 3


def test_unreadable_shard_is_excluded_from_the_batch(tmp_path: Path) -> None:
    seqset = _shard_seqset(tmp_path, ["a.fa", "b.fa"])
    missing = "https://example.test/gone.fa"
    seqset.urls.append(missing)
    seqset.file_classes.append("c")
    store = BatchStore()

    stats = build_store.process_seqset(store, seqset, tmp_path, jobs=2)

    paths, _ = store.batch_calls[0]
    assert [Path(p).name for p in paths] == ["a.fa", "b.fa"]
    assert stats.warnings == 1 and stats.shards_processed == 2


def test_batched_and_serial_ingest_produce_the_same_real_store(
    tmp_path: Path,
) -> None:
    """Batching must not change digests -- sequences are content-addressed."""
    names = ["s1.fa", "s2.fa"]
    payloads = [b">one\nACGT\n>two\nGGTT\n", b">three\nTTTT\n"]
    for cache in ("batched", "serial"):
        for name, payload in zip(names, payloads):
            t = build_store.mirror_cache_path(
                tmp_path / cache, f"https://example.test/{name}"
            )
            t.parent.mkdir(parents=True, exist_ok=True)
            t.write_bytes(payload)

    def run(cache: str, jobs: int) -> list[str]:
        seqset = SeqsetConfig(
            "s", "ns",
            urls=[f"https://example.test/{n}" for n in names],
            file_classes=["c", "c"],
        )
        store = build_store.RefgetStore.in_memory()
        prov: dict[str, dict] = {}
        build_store.process_seqset(
            store, seqset, tmp_path / cache, prov, jobs=jobs
        )
        return [v["collection_digest"] for v in prov.values()]

    assert run("batched", 4) == run("serial", 1)


def test_disk_preflight_aborts_below_threshold_and_warns_on_tight_projection(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    url = "https://example.test/big.fa.gz"
    target = build_store.mirror_cache_path(tmp_path, url)
    target.parent.mkdir(parents=True)
    target.write_bytes(b"x" * 4096)
    sources = [ResolvedSource("seqset", "s", url)]
    store_dir = tmp_path / "store"
    store_dir.mkdir()

    # Plenty of room: no warning, no abort.
    monkeypatch.setattr(build_store, "free_gib", lambda _p: 500.0)
    caplog.set_level("WARNING")
    build_store.check_free_space(tmp_path, store_dir, sources, min_free_gb=25.0)
    assert "may run out part-way" not in caplog.text

    # Below the hard threshold: refuse to start a multi-hour build.
    monkeypatch.setattr(build_store, "free_gib", lambda _p: 3.0)
    with pytest.raises(SystemExit, match="below --min-free-gb"):
        build_store.check_free_space(tmp_path, store_dir, sources, min_free_gb=25.0)

    # Above the threshold but under the projection: warn, still proceed.
    monkeypatch.setattr(build_store, "free_gib", lambda _p: 30.0)
    monkeypatch.setattr(build_store, "STORE_SIZE_FACTOR", 1e9)
    caplog.clear()
    build_store.check_free_space(tmp_path, store_dir, sources, min_free_gb=25.0)
    assert "may run out part-way" in caplog.text


def test_disk_preflight_tolerates_sources_missing_from_cache(tmp_path: Path) -> None:
    store_dir = tmp_path / "store"
    store_dir.mkdir()
    build_store.check_free_space(
        tmp_path, store_dir,
        [ResolvedSource("seqset", "s", "https://example.test/absent.fa")],
        min_free_gb=0.0,
    )
