"""Sync: classifying manifest changes and applying them to a built store.

The classification tests are pure -- no store, no network. The mutation tests
build a real store in ``tmp_path`` via gtars, because the behaviours they guard
(orphan collection needing loaded collections, alias namespaces not shrinking
through the alias API) are properties of gtars, not of this code, and a double
would assert our assumptions rather than the library's behaviour.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import build_lock
import store_sync
from conftest import file_record, v4_lock
from sources import (ResolvedSource, SeqsetConfig, load_config,
                     mirror_cache_path)
from store_sync import MutationOrderError


def _no_network(*args, **kwargs):
    raise AssertionError("test attempted a provider fetch")


def _present(root: Path, url: str) -> Path:
    """Materialize a cache file at the mirrored path a build would use."""
    path = mirror_cache_path(root, url)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_bytes(b"x")
    return path


EN = "https://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/"
TOPLEVEL = EN + "dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
CDNA = EN + "cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"


def seqset(**kw) -> SeqsetConfig:
    base = dict(
        name="ensembl_release_100", namespace="ensembl-100",
        urls=[TOPLEVEL, CDNA], file_classes=["dna.toplevel", "cdna"],
        release=100, rolling_namespace="ensembl",
    )
    base.update(kw)
    return SeqsetConfig(**base)


EXCLUDE = {"file_classes": ["dna.toplevel"], "record_prefixes": ["CHR_"]}


def source(url: str, file_class: str) -> ResolvedSource:
    return ResolvedSource("seqset", "ensembl_release_100", url,
                          file_class=file_class)


def cache(tmp_path: Path, *rels: str) -> Path:
    """Materialize cache files so plan_sync sees them as present."""
    for rel in rels:
        path = tmp_path / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"content of " + rel.encode())
    return tmp_path


REL_TOP = ("ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/"
           "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
REL_CDNA = ("ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/cdna/"
            "Homo_sapiens.GRCh38.cdna.all.fa.gz")


def record(rel: str, *, sha: str | None = None, spec: dict | None = None,
           digest: str = "COLL1", owner: str = "ensembl_release_100") -> dict:
    """A v4 input record, tagged with the collection it fed.

    The tag is consumed by :func:`lock_with`, which inverts it into the
    ``outputs.collections[].from`` lists the real schema stores.
    """
    rec = file_record(rel, owner=owner, url="https://x/" + rel, sha256=sha,
                      ingest_spec=spec)
    rec["_collection"] = digest
    return rec


def lock_with(*records: dict) -> dict:
    files: list[dict] = []
    from_map: dict[str, list[str]] = {}
    for tagged in records:
        rec = dict(tagged)
        digest = rec.pop("_collection", None)
        files.append(rec)
        if digest:
            from_map.setdefault(digest, []).append(rec["cache_path"])
    collections = [
        {"digest": digest, "n_sequences": 10, "from": sorted(paths)}
        for digest, paths in sorted(from_map.items())
    ]
    return v4_lock(files, collections)


# --------------------------------------------------------------------------
# ingest spec
# --------------------------------------------------------------------------

def test_ingest_spec_is_none_when_nothing_transforms_the_source() -> None:
    assert seqset().ingest_spec("dna.toplevel") is None


def test_ingest_spec_records_an_exclusion_only_for_declared_file_classes() -> None:
    entry = seqset(exclude=EXCLUDE)
    assert entry.ingest_spec("dna.toplevel") == {
        "exclude": {"record_prefixes": ["CHR_"]}
    }
    # The sibling cdna source in the same seqset is untransformed. If this
    # leaked, editing the dna.toplevel rule would re-ingest every cdna source.
    assert entry.ingest_spec("cdna") is None


def test_ingest_spec_records_a_derived_format() -> None:
    entry = SeqsetConfig("lrg", "lrg", urls=["https://x/a.zip"], format="lrg_zip")
    assert entry.ingest_spec(None) == {"format": "lrg_zip"}


def test_canonical_spec_digest_is_order_independent() -> None:
    a = build_lock.canonical_spec_sha256({"format": "x", "exclude": {"y": [1]}})
    b = build_lock.canonical_spec_sha256({"exclude": {"y": [1]}, "format": "x"})
    assert a == b
    assert build_lock.canonical_spec_sha256(None) is None


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------

def test_source_matching_bytes_and_spec_is_unchanged(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    sha = build_lock.sha256_file(tmp_path / REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(record(REL_TOP, sha=sha)),
        tmp_path, seqsets=[seqset()],
    )
    assert plan.unchanged == [REL_TOP]
    assert not plan.touched


def test_changed_upstream_bytes_are_reingested(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(record(REL_TOP, sha="0" * 64)),
        tmp_path, seqsets=[seqset()],
    )
    assert plan.reingest == [(REL_TOP, "upstream bytes changed")]
    assert plan.collection_digests_to_remove == ["COLL1"]


def test_changed_ingest_spec_is_reingested_though_bytes_match(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    sha = build_lock.sha256_file(tmp_path / REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(record(REL_TOP, sha=sha)),
        tmp_path, seqsets=[seqset(exclude=EXCLUDE)],
    )
    assert plan.reingest == [(REL_TOP, "ingest spec changed")]


def test_source_absent_from_manifest_is_removed(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    sha = build_lock.sha256_file(tmp_path / REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(record(REL_TOP, sha=sha),
                  record(REL_CDNA, sha=sha, digest="COLL2")),
        tmp_path, seqsets=[seqset()],
    )
    assert [r["cache_path"] for r in plan.removed] == [REL_CDNA]
    assert plan.collection_digests_to_remove == ["COLL2"]


def test_source_absent_from_lock_is_added(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(), tmp_path, seqsets=[seqset()],
    )
    assert plan.added == [REL_TOP]


def test_source_missing_from_cache_is_reported_not_reingested(tmp_path: Path) -> None:
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(record(REL_TOP, sha="0" * 64)),
        tmp_path, seqsets=[seqset()],
    )
    assert plan.missing == [REL_TOP]
    assert not plan.reingest


def test_duplicate_collection_digests_collapse_for_removal(tmp_path: Path) -> None:
    """34 Ensembl releases share 8 collections; each is removed once."""
    cache(tmp_path, REL_TOP, REL_CDNA)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel"), source(CDNA, "cdna")],
        lock_with(
                  record(REL_TOP, sha="0" * 64, digest="SHARED"),
                  record(REL_CDNA, sha="0" * 64, digest="SHARED")),
        tmp_path, seqsets=[seqset()],
    )
    assert len(plan.reingest) == 2
    assert plan.collection_digests_to_remove == ["SHARED"]


# --------------------------------------------------------------------------
# re-ingest scope
# --------------------------------------------------------------------------
#
# Removal is per collection, so a collection dropped because one contributor
# changed also drops every *other* contributor's sequences. Expanding by alias
# namespace -- what sync did before the lock recorded `from` -- covers the
# Ensembl release groups by accident and misses cross-namespace sharing.

def test_reingest_scope_expands_to_every_contributor_of_a_shared_collection(
) -> None:
    lock = lock_with(
        record(REL_TOP, digest="SHARED"),
        record(REL_CDNA, digest="SHARED"),
    )
    assert build_lock.reingest_scope(lock, [REL_TOP]) == {REL_TOP, REL_CDNA}


def test_reingest_scope_crosses_namespaces(tmp_path: Path) -> None:
    """The case namespace expansion cannot reach.

    One RefSeq protein file published under both the per-patch assembly path
    and an annotation-release path, byte-identical, so they yield one
    collection -- but owned by different seqsets under different alias
    namespaces. Expanding by namespace re-ingests only the one that changed and
    silently drops the other's sequences.

    This shape does not occur in the committed lock: its 31 multi-owner
    collections are all Ensembl release groups, which namespace expansion
    happens to cover. That is the point of the test -- the old approximation was
    right by coincidence, and splitting one seqset's URL list would break it.
    """
    patch_rel = ("ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
                 "GCF_000001405.25_GRCh37.p13/GCF_..._protein.faa.gz")
    release_rel = ("ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/"
                   "105.20220307/GCF_..._protein.faa.gz")
    lock = lock_with(
        record(patch_rel, digest="SHARED", owner="refseq_history_grch37_protein"),
        record(release_rel, digest="SHARED", owner="refseq_annotation_105"),
    )
    scope = build_lock.reingest_scope(lock, [patch_rel])
    assert scope == {patch_rel, release_rel}
    assert build_lock.owners_for_paths(lock, scope) == {
        ("seqset", "refseq_history_grch37_protein"),
        ("seqset", "refseq_annotation_105"),
    }


def test_reingest_scope_leaves_unrelated_collections_alone() -> None:
    lock = lock_with(
        record(REL_TOP, digest="COLL1"),
        record(REL_CDNA, digest="COLL2"),
    )
    assert build_lock.reingest_scope(lock, [REL_TOP]) == {REL_TOP}


def test_reingest_scope_of_an_unattributed_file_is_itself() -> None:
    lock = lock_with(record(REL_TOP, digest=None))
    assert build_lock.reingest_scope(lock, [REL_TOP]) == {REL_TOP}


def test_committed_config_against_committed_lock_targets_only_padded_releases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Guards the real migration: releases 76-109 plus the one lrg_zip source.

    Sources come from the lock via ``apply_locked_sources`` rather than from
    live discovery. Resolving them for real would fetch provider manifests, so
    the test would depend on what NCBI and Ensembl are publishing today -- it
    failed exactly that way when a ninth RefSeqGene shard appeared upstream.
    """
    import re

    repo = Path(__file__).resolve().parent.parent
    lock_path = repo / "build.lock.json"
    if not lock_path.exists():  # a build artifact, not always present
        pytest.skip("build.lock.json not present")
    monkeypatch.setattr("sources._fetch_text", _no_network)

    _, seqsets = load_config(repo / "sources.toml")
    lock = build_lock.load_lock(lock_path)
    sources = build_lock.apply_locked_sources(seqsets, lock)
    # Mirror the cache layout under tmp_path so presence checks pass without
    # depending on the real downloads/ tree. mirror_cache_path is a pure
    # URL-to-path mapping, so the relative keys still match the lock's.
    for source in sources:
        _present(tmp_path, source.url)

    plan = build_lock.plan_sync(
        sources, lock, tmp_path, seqsets=seqsets, hash_files=False,
    )
    assert not plan.missing
    releases = {
        int(m.group(1))
        for rel, _ in plan.reingest
        if (m := re.search(r"release-(\d+)", rel))
    }
    assert releases == set(range(76, 110))
    # Everything re-ingested is either a padded Ensembl release or the lrg_zip
    # source, whose conversion the lock cannot vouch for.
    non_ensembl = [rel for rel, _ in plan.reingest if "release-" not in rel]
    assert len(non_ensembl) == 1 and non_ensembl[0].endswith(".zip")


# --------------------------------------------------------------------------
# store mutation
# --------------------------------------------------------------------------

def build_store_fixture(tmp_path: Path):
    """Two collections sharing a sequence, each with one unique record."""
    from gtars.refget import RefgetStore

    fa = tmp_path / "fa"
    fa.mkdir()
    base = "ACGT" * 40
    (fa / "a.fa").write_text(f">SHARED\n{base}\n>ONLY_A\n{'AACC' * 30}\n")
    (fa / "b.fa").write_text(f">SHARED\n{base}\n>ONLY_B\n{'GGTT' * 30}\n")
    store_dir = tmp_path / "store"
    store_dir.mkdir()
    store = RefgetStore.on_disk(str(store_dir))
    store.set_quiet(True)
    digests = [
        store.add_sequence_collection_from_fasta(str(fa / name))[0].digest
        for name in ("a.fa", "b.fa")
    ]
    store.write()
    return store_dir, digests


def test_removal_without_loading_collections_raises(tmp_path: Path) -> None:
    """gtars silently collects nothing on a lazily loaded store, reclaiming
    zero bytes while appearing to succeed. Fail loudly instead."""
    from gtars.refget import RefgetStore

    store_dir, digests = build_store_fixture(tmp_path)
    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    with pytest.raises(MutationOrderError, match="not fully loaded"):
        store_sync.remove_collections(store, [digests[0]])


def test_removal_collects_orphans_but_keeps_shared_sequences(tmp_path: Path) -> None:
    from gtars.refget import RefgetStore

    store_dir, digests = build_store_fixture(tmp_path)
    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    before = store_sync.store_sequence_count(store)
    store_sync.load_for_mutation(store)
    assert store_sync.remove_collections(store, [digests[0]]) == 1
    store.write()
    # SHARED is still referenced by collection B; only ONLY_A is orphaned.
    assert store_sync.store_sequence_count(store) == before - 1


def test_reconcile_namespace_shrinks_a_namespace(tmp_path: Path) -> None:
    """The alias API cannot shrink a namespace: load_sequence_aliases merges.
    reconcile_namespace writes the file, so surplus aliases actually go."""
    from gtars.refget import RefgetStore

    store_dir, digests = build_store_fixture(tmp_path)
    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    names = {
        r.metadata.name: r.metadata.sha512t24u
        for r in store.get_collection(digests[0]).sequences
    }
    for name, digest in names.items():
        store.add_sequence_alias("ns", name, digest)
    store.write()

    keep = {"SHARED": names["SHARED"]}
    added, removed = store_sync.reconcile_namespace(store_dir, "ns", keep)
    assert (added, removed) == (0, 1)
    reopened = RefgetStore.open_local(str(store_dir))
    assert sorted(reopened.list_sequence_aliases("ns")) == ["SHARED"]


def test_reconcile_namespace_rejects_an_unsafe_name(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unsafe namespace"):
        store_sync.reconcile_namespace(tmp_path, "../escape", {})
