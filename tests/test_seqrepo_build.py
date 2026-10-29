"""``gks-refgetstore seqrepo`` against a real gtars store and real seqrepo writes.

The store is small but carries every alias shape the genomic store has: a
RefSeq accession, assembly-scoped names (``GRCh38:1``, ``GRCh38:chr1``, a
GenBank accession in the assembly namespace), ``insdc``, both LRG spellings of
one genomic record, and an LRG protein with selenocysteine (``U``).
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

from conftest import file_record, lock_for_store, write_fasta
from gks_refgetstore import build_lock, cli, seqrepo_build, store_census

pytest.importorskip("biocommons.seqrepo")

needs_bgzip = pytest.mark.skipif(shutil.which("bgzip") is None,
                                 reason="seqrepo commits through htslib's bgzip")

GENOME = [
    ("NC_000001.11", "ACGTNACGTTTGACCAGTACGATTACAGATTACAGGGCCCAATTN"),
    ("NC_012920.1", "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCC"),
]
LRG = [
    ("LRG_1g", "TTGACCAGTAGGGATTTACCCAGTAGGACATTACAGGACCA"),
    ("LRG_1p1", "MKTAYIAKQRQISFVKSHFSRQUEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ"),
]

MANIFEST = """
[[seqset]]
name = "genome"
namespace = "refseq"
urls = ["https://example.test/genome.fa"]
"""
LRG_ENTRY = """
[[seqset]]
name = "lrg_set"
namespace = "lrg"
urls = ["https://example.test/lrg.fa"]
"""


def _digest(store, ns: str, alias: str) -> str:
    return store.get_sequence_metadata_by_alias(ns, alias).sha512t24u


@pytest.fixture
def aliased(tmp_path: Path):
    """``(store_dir, lock_path, manifest_path, subset_manifest_path)``."""
    from gtars.refget import RefgetStore

    store_dir = tmp_path / "store"
    store = RefgetStore.on_disk(str(store_dir))
    store.set_quiet(True)
    genome, _ = store.add_sequence_collection_from_fasta(
        str(write_fasta(tmp_path / "genome.fa", GENOME)))
    lrg, _ = store.add_sequence_collection_from_fasta(
        str(write_fasta(tmp_path / "lrg.fa", LRG)))
    by_name = {name: store.get_sequence_by_name(genome.digest, name).metadata.sha512t24u
               for name, _ in GENOME}
    by_name |= {name: store.get_sequence_by_name(lrg.digest, name).metadata.sha512t24u
                for name, _ in LRG}
    for ns, alias, name in [
        ("refseq", "NC_000001.11", "NC_000001.11"),
        ("refseq", "NC_012920.1", "NC_012920.1"),
        ("GRCh38", "1", "NC_000001.11"),
        ("GRCh38", "chr1", "NC_000001.11"),
        ("GRCh38", "CM000663.2", "NC_000001.11"),
        ("GRCh38", "MT", "NC_012920.1"),
        ("insdc", "CM000663.2", "NC_000001.11"),
        ("lrg", "LRG_1g", "LRG_1g"),
        ("lrg", "LRG_1", "LRG_1g"),
        ("lrg", "LRG_1p1", "LRG_1p1"),
    ]:
        store.add_sequence_alias(ns, alias, by_name[name])
    store.write()

    files = [file_record("example.test/genome.fa", owner="genome",
                         url="https://example.test/genome.fa"),
             file_record("example.test/lrg.fa", owner="lrg_set",
                         url="https://example.test/lrg.fa")]
    reopened = RefgetStore.open_local(str(store_dir))
    lock = lock_for_store(reopened, store_dir, files, from_map={
        genome.digest: {"example.test/genome.fa"},
        lrg.digest: {"example.test/lrg.fa"},
    })
    lock_path = tmp_path / "build.lock.json"
    build_lock.write_lock(lock_path, lock)
    manifest = tmp_path / "sources.toml"
    manifest.write_text(MANIFEST + LRG_ENTRY)
    subset = tmp_path / "sources.subset.toml"
    subset.write_text(MANIFEST)
    return store_dir, lock_path, manifest, subset


def _run(store_dir, lock_path, manifest, target, *extra) -> int:
    return cli.main(["seqrepo", "--config", str(manifest), "--lock", str(lock_path),
                     "--store-dir", str(store_dir), "--target-dir", str(target),
                     "--instance-name", "tiny", *extra])


def _open(path: Path, **kwargs):
    from biocommons.seqrepo import SeqRepo

    return SeqRepo(str(path), **kwargs)


@needs_bgzip
def test_the_instance_mirrors_the_store(aliased, tmp_path: Path) -> None:
    from gtars.refget import RefgetStore

    store_dir, lock_path, manifest, _ = aliased
    target = tmp_path / "seqrepo"
    assert _run(store_dir, lock_path, manifest, target) == 0

    store = RefgetStore.open_local(str(store_dir))
    sr = _open(target / "tiny")
    assert {r["seq_id"] for r in sr.sequences} == {
        m.sha512t24u for m in store.list_sequences()}
    # seqrepo stores refseq as NCBI and lrg as LRG, and answers both spellings.
    nc1 = _digest(store, "refseq", "NC_000001.11")
    for ns, alias in [("refseq", "NC_000001.11"), ("NCBI", "NC_000001.11"),
                      ("GRCh38", "1"), ("GRCh38", "chr1"), ("GRCh38", "CM000663.2"),
                      ("insdc", "CM000663.2")]:
        assert {r["seq_id"] for r in sr.aliases.find_aliases(namespace=ns, alias=alias)
                } == {nc1}
    lrg1 = _digest(store, "lrg", "LRG_1")
    assert {r["namespace"] for r in sr.aliases.find_aliases(seq_id=lrg1)} >= {"LRG", "lrg"}
    # The selenoprotein is stored with its U.
    p1 = _digest(store, "lrg", "LRG_1p1")
    assert "U" in sr.sequences.fetch(p1)
    assert sr.aliases.find_aliases(namespace="ga4gh", alias=f"SQ.{p1}")

    record = json.loads((target / "tiny.provenance.json").read_text())
    assert record["verification"]["ok"] is True
    assert record["write"]["sequences"] == 4
    assert record["verification"]["counts"]["store_aliases"] == 10
    # 10 store aliases plus seqrepo's four digest aliases per sequence.
    assert record["verification"]["counts"]["seqrepo_alias_rows"] == 10 + 4 * 4
    # Written in a temporary sibling and renamed: nothing else is left behind.
    assert sorted(p.name for p in target.iterdir()) == ["tiny", "tiny.provenance.json"]
    assert (target / "tiny").stat().st_mode & 0o777 == 0o755


@needs_bgzip
def test_an_existing_instance_is_refused(aliased, tmp_path: Path) -> None:
    store_dir, lock_path, manifest, _ = aliased
    (tmp_path / "seqrepo" / "tiny").mkdir(parents=True)
    with pytest.raises(SystemExit, match="already exists"):
        _run(store_dir, lock_path, manifest, tmp_path / "seqrepo")


@needs_bgzip
def test_verify_only_reports_an_alias_the_store_lacks(aliased, tmp_path: Path) -> None:
    store_dir, lock_path, manifest, _ = aliased
    target = tmp_path / "seqrepo"
    assert _run(store_dir, lock_path, manifest, target) == 0
    sr = _open(target / "tiny", writeable=True)
    seq_id = next(iter(sr.sequences))["seq_id"]
    sr.aliases.store_alias(seq_id=seq_id, namespace="NCBI35", alias="21")
    sr.commit()
    del sr
    assert _run(store_dir, lock_path, manifest, target, "--verify-only") == 1


def test_a_subset_manifest_selects_its_own_collections(aliased) -> None:
    from gtars.refget import RefgetStore

    store_dir, lock_path, _, subset = aliased
    lock = build_lock.load_lock(lock_path)
    with pytest.raises(SystemExit, match="lock-check-mode subset"):
        seqrepo_build.expected_collections(subset, lock, "strict")
    expected = seqrepo_build.expected_collections(subset, lock, "subset")
    store = RefgetStore.open_local(str(store_dir))
    assert len(expected) == 1
    (unpinned,) = set(store_census.collection_census(store)) - set(expected)
    # The store also holds the LRG collection, which the subset does not pin;
    # its sequences are then outside the pinned collections too.
    assert seqrepo_build.check_store_scope(store, expected) == [
        f"collection {unpinned} is in the store but not pinned for this manifest",
        ("the store indexes 2 sequence(s) outside the pinned collections and "
         "lacks 0 of their members"),
    ]


def test_the_full_manifest_matches_the_store(aliased) -> None:
    from gtars.refget import RefgetStore

    store_dir, lock_path, manifest, _ = aliased
    lock = build_lock.load_lock(lock_path)
    expected = seqrepo_build.expected_collections(manifest, lock, "strict")
    store = RefgetStore.open_local(str(store_dir))
    assert seqrepo_build.check_store_scope(store, expected) == []


@needs_bgzip
def test_a_sequence_seqrepo_would_renormalize_is_refused(tmp_path: Path) -> None:
    from gtars.refget import RefgetStore

    store = RefgetStore.on_disk(str(tmp_path / "store"))
    store.set_quiet(True)
    store.add_sequence_collection_from_fasta(
        str(write_fasta(tmp_path / "pep.fa", [("ENSP1.1", "MKTAYIAKQ*")])))
    lib = seqrepo_build.require_seqrepo()
    with pytest.raises(RuntimeError, match="did not file sequence"):
        seqrepo_build.write_instance(lib, store, tmp_path / "sr")


@needs_bgzip
def test_add_missing_stores_only_what_the_instance_lacks(tmp_path: Path) -> None:
    from gtars.refget import RefgetStore

    store_dir = tmp_path / "store"
    store = RefgetStore.on_disk(str(store_dir))
    store.set_quiet(True)
    genome, _ = store.add_sequence_collection_from_fasta(
        str(write_fasta(tmp_path / "genome.fa", GENOME)))
    mt = store.get_sequence_by_name(genome.digest, "NC_012920.1").metadata.sha512t24u
    for name, _ in GENOME:
        digest = store.get_sequence_by_name(genome.digest, name).metadata.sha512t24u
        store.add_sequence_alias("refseq", name, digest)
    store.write()
    lib = seqrepo_build.require_seqrepo()
    instance = tmp_path / "sr"
    first = seqrepo_build.write_instance(lib, RefgetStore.open_local(str(store_dir)),
                                         instance)
    assert (first.sequences, first.completed) == (2, 0)
    before = {p: p.read_bytes() for p in instance.rglob("*.bgz")}

    # The store gains a collection, and an existing sequence gains an alias.
    store = RefgetStore.on_disk(str(store_dir))
    store.set_quiet(True)
    lrg, _ = store.add_sequence_collection_from_fasta(
        str(write_fasta(tmp_path / "lrg.fa", LRG)))
    for name, _ in LRG:
        digest = store.get_sequence_by_name(lrg.digest, name).metadata.sha512t24u
        store.add_sequence_alias("lrg", name, digest)
    store.add_sequence_alias("insdc", "J01415.2", mt)
    store.write()
    store = RefgetStore.open_local(str(store_dir))

    added = seqrepo_build.write_instance(lib, store, instance, add_missing=True)
    assert (added.sequences, added.completed, added.aliases_added) == (2, 1, 3)
    # The existing bgzip files are untouched; the new sequences went to a new one.
    assert all(p.read_bytes() == data for p, data in before.items())
    assert len(list(instance.rglob("*.bgz"))) == len(before) + 1
    assert seqrepo_build.verify_instance(lib, store, instance).ok

    again = seqrepo_build.write_instance(lib, store, instance, add_missing=True)
    assert (again.sequences, again.completed, again.aliases_added) == (0, 0, 0)


@needs_bgzip
def test_add_missing_needs_an_existing_instance(aliased, tmp_path: Path) -> None:
    store_dir, lock_path, manifest, _ = aliased
    with pytest.raises(SystemExit, match="no SeqRepo instance"):
        _run(store_dir, lock_path, manifest, tmp_path / "seqrepo", "--add-missing")
    assert _run(store_dir, lock_path, manifest, tmp_path / "seqrepo") == 0
    # A complete instance has nothing to add and still verifies.
    assert _run(store_dir, lock_path, manifest, tmp_path / "seqrepo", "--add-missing") == 0
    record = json.loads((tmp_path / "seqrepo" / "tiny.provenance.json").read_text())
    assert record["write"]["mode"] == "add-missing"
    assert record["write"]["sequences"] == 0
    assert record["verification"]["ok"] is True


def test_a_missing_extra_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "biocommons.seqrepo", None)
    with pytest.raises(SystemExit, match="--extra seqrepo"):
        seqrepo_build.require_seqrepo()


def test_a_missing_bgzip_says_so(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("SEQREPO_BGZIP_PATH", str(tmp_path / "no-bgzip"))
    with pytest.raises(SystemExit, match="bgzip"):
        seqrepo_build.require_seqrepo()
