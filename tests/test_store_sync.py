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
from build_store import ResolvedSource, SeqsetConfig
from store_sync import MutationOrderError

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


def lock_with(schema: str, *records: dict) -> dict:
    return {"schema": schema, "build": {}, "sources": list(records)}


def record(rel: str, *, sha: str | None = None, spec: dict | None = None,
           digest: str = "COLL1", owner: str = "ensembl_release_100") -> dict:
    rec = {
        "kind": "seqset", "owner": owner, "url": "https://x/" + rel,
        "cache_path": rel, "collection_digest": digest, "n_sequences": 10,
    }
    if sha is not None:
        rec["sha256"] = sha
    if spec is not None or True:
        rec["ingest_spec"] = spec
        rec["ingest_spec_sha256"] = build_lock.canonical_spec_sha256(spec)
    return rec


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
        lock_with(build_lock.SCHEMA, record(REL_TOP, sha=sha)),
        tmp_path, seqsets=[seqset()],
    )
    assert plan.unchanged == [REL_TOP]
    assert not plan.touched


def test_changed_upstream_bytes_are_reingested(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(build_lock.SCHEMA, record(REL_TOP, sha="0" * 64)),
        tmp_path, seqsets=[seqset()],
    )
    assert plan.reingest == [(REL_TOP, "upstream bytes changed")]
    assert plan.collection_digests_to_remove == ["COLL1"]


def test_changed_ingest_spec_is_reingested_though_bytes_match(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    sha = build_lock.sha256_file(tmp_path / REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(build_lock.SCHEMA, record(REL_TOP, sha=sha)),
        tmp_path, seqsets=[seqset(exclude=EXCLUDE)],
    )
    assert plan.reingest == [(REL_TOP, "ingest spec changed")]


def test_source_absent_from_manifest_is_removed(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    sha = build_lock.sha256_file(tmp_path / REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(build_lock.SCHEMA, record(REL_TOP, sha=sha),
                  record(REL_CDNA, sha=sha, digest="COLL2")),
        tmp_path, seqsets=[seqset()],
    )
    assert [r["cache_path"] for r in plan.removed] == [REL_CDNA]
    assert plan.collection_digests_to_remove == ["COLL2"]


def test_source_absent_from_lock_is_added(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(build_lock.SCHEMA), tmp_path, seqsets=[seqset()],
    )
    assert plan.added == [REL_TOP]


def test_source_missing_from_cache_is_reported_not_reingested(tmp_path: Path) -> None:
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")],
        lock_with(build_lock.SCHEMA, record(REL_TOP, sha="0" * 64)),
        tmp_path, seqsets=[seqset()],
    )
    assert plan.missing == [REL_TOP]
    assert not plan.reingest


def test_duplicate_collection_digests_collapse_for_removal(tmp_path: Path) -> None:
    """34 Ensembl releases share 8 collections; each is removed once."""
    cache(tmp_path, REL_TOP, REL_CDNA)
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel"), source(CDNA, "cdna")],
        lock_with(build_lock.SCHEMA,
                  record(REL_TOP, sha="0" * 64, digest="SHARED"),
                  record(REL_CDNA, sha="0" * 64, digest="SHARED")),
        tmp_path, seqsets=[seqset()],
    )
    assert len(plan.reingest) == 2
    assert plan.collection_digests_to_remove == ["SHARED"]


# --------------------------------------------------------------------------
# migrating a lock written before ingest_spec existed
# --------------------------------------------------------------------------

def test_legacy_lock_leaves_untransformed_sources_alone(tmp_path: Path) -> None:
    cache(tmp_path, REL_TOP)
    sha = build_lock.sha256_file(tmp_path / REL_TOP)
    legacy = lock_with(build_lock.V2_SCHEMA, record(REL_TOP, sha=sha))
    for rec in legacy["sources"]:
        rec.pop("ingest_spec"), rec.pop("ingest_spec_sha256")
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")], legacy, tmp_path,
        seqsets=[seqset()],
    )
    assert plan.unchanged == [REL_TOP]
    assert plan.legacy_schema


def test_legacy_lock_reingests_sources_that_now_declare_a_transform(
    tmp_path: Path,
) -> None:
    """The migration case: an absent spec is unknown, so a declared transform
    cannot be assumed to have been applied."""
    cache(tmp_path, REL_TOP)
    sha = build_lock.sha256_file(tmp_path / REL_TOP)
    legacy = lock_with(build_lock.V2_SCHEMA, record(REL_TOP, sha=sha))
    for rec in legacy["sources"]:
        rec.pop("ingest_spec"), rec.pop("ingest_spec_sha256")
    plan = build_lock.plan_sync(
        [source(TOPLEVEL, "dna.toplevel")], legacy, tmp_path,
        seqsets=[seqset(exclude=EXCLUDE)],
    )
    assert plan.reingest == [(REL_TOP, "ingest spec not pinned by this lock")]


def test_live_config_against_the_committed_lock_targets_only_padded_releases() -> None:
    """Guards the real migration: releases 76-109 plus the one lrg_zip source."""
    import re

    import build_store

    repo = Path(__file__).resolve().parent.parent
    lock_path = repo / "build.lock.json"
    if not lock_path.exists():  # lock is a build artifact, not always present
        pytest.skip("build.lock.json not present")
    assemblies, seqsets = build_store.load_config(repo / "sources.toml")
    lock = build_lock.load_lock(lock_path)
    plan = build_lock.plan_sync(
        build_store.resolve_sources(assemblies, seqsets), lock,
        repo / "downloads", seqsets=seqsets, hash_files=False,
    )
    assert not plan.added and not plan.removed
    releases = {
        int(m.group(1))
        for rel, _ in plan.reingest
        if (m := re.search(r"release-(\d+)", rel))
    }
    assert releases == set(range(76, 110))


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
