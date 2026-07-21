from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

import build_lock
import build_store
from build_store import (ResolvedSource, SeqsetConfig, apply_locked_sources,
                         load_config, parse_checksum_manifest, resolve_sources)


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


@pytest.mark.parametrize("text", [
    "not-a-record", "0" * 32 + "  dir/file.gz", "0" * 32 + "  ../file.gz",
    "0" * 32 + "  file.gz\n" + "1" * 32 + "  file.gz",
    "0" * 32 + " *file.gz",
])
def test_manifest_rejects_malformed_duplicate_or_unsafe_records(text: str) -> None:
    with pytest.raises(ValueError):
        parse_checksum_manifest(text)


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

    def ingest(_store, entry, _cache, _provenance, refresh_derived=False):
        seen.append(refresh_derived)
        return build_store.SeqsetStats(entry.name, entry.namespace)

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
        def list_collection_alias_namespaces(): return []

    monkeypatch.setattr(build_store, "RefgetStore", Store)
    args = type("Args", (), {
        "config": config, "cache_dir": cache, "store_dir": tmp_path / "store",
        "lock": lock_path, "locked_sources": True, "skip_assemblies": False,
        "skip_seqsets": False, "assembly": None, "seqset": None,
        "lock_check_mode": "strict", "force_download": True,
        "force_lock": False, "no_lock": True,
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
    sources = resolve_sources(assemblies, seqsets,
                              lambda url: manifests[Path(url).name])
    urls = [source.url for source in sources]
    assert len(urls) == len(set(urls))
    assert all("//" not in url.split("://", 1)[1] for url in urls)
