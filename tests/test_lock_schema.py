"""Schema /4: record shape, the invariants ``validate_lock`` enforces, and the
accessors every reader is required to go through."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from conftest import file_record, v4_lock
from gks_refgetstore import build_lock, store_census
from gks_refgetstore.sources import ResolvedSource

REL_A = "host/a.fa.gz"
REL_B = "host/b.fa.gz"


# --------------------------------------------------------------------- shape

def test_a_file_record_carries_every_declared_field(tmp_path: Path) -> None:
    source = ResolvedSource("seqset", "owner", "https://example.test/a.fa.gz")
    record = build_lock.file_record(source, tmp_path, hash_files=False)
    assert set(record) == set(build_lock.FILE_FIELDS)


def test_an_absent_cache_file_is_recorded_as_absent(tmp_path: Path) -> None:
    source = ResolvedSource("seqset", "owner", "https://example.test/a.fa.gz")
    record = build_lock.file_record(source, tmp_path)
    assert record["present"] is False
    assert record["bytes"] is None and record["sha256"] is None


def test_upstream_md5_is_reconstructed_not_stored() -> None:
    md5 = "a" * 32
    record = file_record(REL_A, provider_checksum=md5,
                         provider_checksum_algorithm="md5")
    assert "upstream_md5" not in record
    assert build_lock.upstream_md5_of(record) == md5


def test_a_bsd_sum_is_not_mistaken_for_an_md5() -> None:
    record = file_record(REL_A, provider_checksum="01234",
                         provider_checksum_algorithm="bsd-sum",
                         provider_checksum_blocks=17)
    assert build_lock.upstream_md5_of(record) is None


def test_a_source_with_no_provider_checksum_reconstructs_to_none() -> None:
    assert build_lock.upstream_md5_of(file_record(REL_A)) is None


# ---------------------------------------------------------------- validation

def test_a_hand_built_lock_validates() -> None:
    lock = v4_lock([file_record(REL_A)], [{"digest": "C1", "from": [REL_A]}])
    assert build_lock.validate_lock(lock) is lock


@pytest.mark.parametrize("schema", [
    "gks-refgetstore-build-lock/1",
    "gks-refgetstore-build-lock/2",
    "gks-refgetstore-build-lock/3",
    None,
])
def test_every_pre_v4_schema_is_rejected(schema) -> None:
    lock = v4_lock()
    lock["schema"] = schema
    with pytest.raises(build_lock.LockError, match="unsupported build lock schema"):
        build_lock.validate_lock(lock)


@pytest.mark.parametrize("section", ["build", "inputs", "outputs"])
def test_a_missing_section_is_rejected(section: str) -> None:
    lock = v4_lock()
    del lock[section]
    with pytest.raises(build_lock.LockError, match=f"missing its {section!r}"):
        build_lock.validate_lock(lock)


def test_a_record_missing_a_field_is_rejected() -> None:
    record = file_record(REL_A)
    del record["mutable"]
    with pytest.raises(build_lock.LockError, match="missing mutable"):
        build_lock.validate_lock(v4_lock([record]))


def test_a_duplicate_cache_path_is_rejected() -> None:
    lock = v4_lock([file_record(REL_A), file_record(REL_A, url="https://x/other")])
    with pytest.raises(build_lock.LockError, match="duplicate cache_path"):
        build_lock.validate_lock(lock)


def test_a_duplicate_url_is_rejected() -> None:
    lock = v4_lock([file_record(REL_A), file_record(REL_B, url=file_record(REL_A)["url"])])
    with pytest.raises(build_lock.LockError, match="duplicate url"):
        build_lock.validate_lock(lock)


def test_an_ingest_spec_digest_that_does_not_match_its_spec_is_rejected() -> None:
    record = file_record(REL_A, ingest_spec={"format": "lrg_zip"},
                         ingest_spec_sha256="0" * 64)
    with pytest.raises(build_lock.LockError, match="does not digest ingest_spec"):
        build_lock.validate_lock(v4_lock([record]))


def test_a_null_ingest_spec_must_carry_a_null_digest() -> None:
    """The conversion writes null/null for every file; a stray digest beside a
    null spec would claim a transformation that never happened."""
    record = file_record(REL_A, ingest_spec=None, ingest_spec_sha256="0" * 64)
    with pytest.raises(build_lock.LockError, match="does not digest ingest_spec"):
        build_lock.validate_lock(v4_lock([record]))


def test_a_collection_count_that_disagrees_with_the_list_is_rejected() -> None:
    lock = v4_lock([], [{"digest": "C1"}])
    lock["outputs"]["n_collections"] = 7
    with pytest.raises(build_lock.LockError, match="disagrees with 1"):
        build_lock.validate_lock(lock)


def test_a_stale_collections_root_is_rejected() -> None:
    lock = v4_lock([], [{"digest": "C1"}])
    lock["outputs"]["collections_root"] = store_census.digest_root(["C2"])
    with pytest.raises(build_lock.LockError, match="collections_root"):
        build_lock.validate_lock(lock)


def test_a_sequences_root_that_is_not_a_root_is_rejected() -> None:
    lock = v4_lock()
    lock["outputs"]["sequences_root"] = "deadbeef"
    with pytest.raises(build_lock.LockError, match="not a sha256 root"):
        build_lock.validate_lock(lock)


def test_a_collection_claiming_an_unknown_contributor_is_rejected() -> None:
    lock = v4_lock([file_record(REL_A)], [{"digest": "C1", "from": [REL_B]}])
    with pytest.raises(build_lock.LockError, match="unknown contributor"):
        build_lock.validate_lock(lock)


def test_a_file_attributed_to_two_collections_is_rejected() -> None:
    """One collection can have many files; one file cannot have two collections.

    The /2 schema could not express the first and the accessors depend on the
    second, so it is enforced rather than assumed.
    """
    lock = v4_lock(
        [file_record(REL_A)],
        [{"digest": "C1", "from": [REL_A]}, {"digest": "C2", "from": [REL_A]}],
    )
    with pytest.raises(build_lock.LockError, match="more than one collection"):
        build_lock.validate_lock(lock)


def test_a_collection_with_no_known_contributor_is_allowed() -> None:
    """Honesty over tidiness: an empty ``from`` records that the store holds a
    collection no source file explains, rather than hiding it."""
    build_lock.validate_lock(v4_lock([], [{"digest": "C1", "from": []}]))


# ---------------------------------------------------------------- accessors

def test_accessors_invert_each_other() -> None:
    lock = v4_lock(
        [file_record(REL_A), file_record(REL_B)],
        [{"digest": "SHARED", "from": [REL_A, REL_B]}],
    )
    assert build_lock.files_by_collection(lock) == {"SHARED": [REL_A, REL_B]}
    assert build_lock.collection_by_file(lock) == {
        REL_A: "SHARED", REL_B: "SHARED",
    }


def test_owners_for_paths_reports_kind_and_owner() -> None:
    lock = v4_lock([
        file_record(REL_A, kind="seqset", owner="one"),
        file_record(REL_B, kind="assembly_fasta", owner="GRCh37"),
    ])
    assert build_lock.owners_for_paths(lock, {REL_A, REL_B}) == {
        ("seqset", "one"), ("assembly_fasta", "GRCh37"),
    }


def test_accessors_on_an_empty_lock_are_empty() -> None:
    lock = v4_lock()
    assert build_lock.lock_files(lock) == []
    assert build_lock.lock_collections(lock) == []
    assert build_lock.collection_by_file(lock) == {}


# ------------------------------------------------------- apply_locked_sources

def test_locked_sources_is_built_by_keyword_not_by_position() -> None:
    """The regression test for the positional-construction bug.

    ``apply_locked_sources`` used to build ``ResolvedSource`` positionally from
    a series of ``dict.get`` calls. ``ResolvedSource`` has nine fields, six of
    them optional strings, so a reordering on either side produced a
    valid-looking source with its checksum in the algorithm slot -- and nothing
    would have noticed. Shuffling the JSON key order must change nothing, since
    JSON objects are unordered by definition.
    """
    record = file_record(
        REL_A, kind="seqset", owner="one", url="https://example.test/one.fa",
        provider_checksum="00065", provider_checksum_algorithm="bsd-sum",
        provider_checksum_blocks=17, checksum_url="https://example.test/CHECKSUMS",
        file_class="dna.toplevel",
    )
    expected = ResolvedSource(
        kind="seqset", owner="one", url="https://example.test/one.fa",
        upstream_md5=None, provider_checksum="00065",
        provider_checksum_algorithm="bsd-sum", provider_checksum_blocks=17,
        checksum_url="https://example.test/CHECKSUMS", file_class="dna.toplevel",
    )
    rng = random.Random(0)
    for _ in range(8):
        keys = list(record)
        rng.shuffle(keys)
        shuffled = {key: record[key] for key in keys}
        # Round-trip through JSON so the shuffle is the real thing a reader
        # would face, not just a dict with a different insertion order.
        reloaded = json.loads(json.dumps({"k": shuffled}))["k"]
        assert build_lock.apply_locked_sources([], v4_lock([reloaded])) == [expected]


def test_locked_sources_validates_the_lock_it_is_handed() -> None:
    with pytest.raises(build_lock.LockError):
        build_lock.apply_locked_sources([], {"schema": "something-else"})


# ----------------------------------------------------------------- round trip

def test_write_then_load_is_the_identity(tmp_path: Path) -> None:
    lock = v4_lock([file_record(REL_A)], [{"digest": "C1", "from": [REL_A]}])
    path = tmp_path / "build.lock.json"
    build_lock.write_lock(path, lock)
    assert build_lock.load_lock(path) == lock


def test_load_lock_rejects_an_invalid_lock_on_disk(tmp_path: Path) -> None:
    path = tmp_path / "build.lock.json"
    path.write_text(json.dumps({"schema": "gks-refgetstore-build-lock/2",
                                "sources": []}))
    with pytest.raises(build_lock.LockError):
        build_lock.load_lock(path)


# -------------------------------------------------------------- store outputs

def test_outputs_are_censused_from_the_store_not_from_provenance(
    tiny_store,
) -> None:
    outputs = build_lock.store_outputs(tiny_store, {})
    assert outputs["n_collections"] == 2
    assert outputs["n_sequences"] == 3
    # No provenance was supplied, so every collection is honestly unattributed
    # rather than omitted.
    assert all(c["from"] == [] for c in outputs["collections"])


def test_outputs_counts_are_integers(tiny_store) -> None:
    """store.stats() returns strings; a lock recording "3" is not the 3 it
    describes, and compares unequal to itself after a round trip."""
    outputs = build_lock.store_outputs(tiny_store, {})
    assert isinstance(outputs["n_sequences"], int)
    assert isinstance(outputs["n_collections"], int)


def test_outputs_roots_agree_with_their_own_contents(tiny_store) -> None:
    outputs = build_lock.store_outputs(tiny_store, {})
    assert outputs["collections_root"] == store_census.digest_root(
        c["digest"] for c in outputs["collections"]
    )
    assert outputs["sequences_root"] == store_census.digest_root(
        store_census.sequence_digests(tiny_store)
    )


def test_provenance_only_annotates_the_census(tiny_store, tmp_path: Path) -> None:
    digests = sorted(store_census.collection_census(tiny_store))
    provenance = {
        str(tmp_path / "downloads" / "alpha.fa"): {
            "collection_digest": digests[0], "n_sequences": 2,
        },
        # A digest the store does not hold: provenance cannot invent a
        # collection, it can only label one the census already found.
        str(tmp_path / "downloads" / "ghost.fa"): {
            "collection_digest": "NOT_IN_STORE", "n_sequences": 1,
        },
    }
    from_map = build_lock.provenance_from_map(tmp_path / "downloads", provenance)
    outputs = build_lock.store_outputs(tiny_store, from_map)
    by_digest = {c["digest"]: c for c in outputs["collections"]}
    assert by_digest.keys() == set(digests)
    assert by_digest[digests[0]]["from"] == ["alpha.fa"]
    assert by_digest[digests[1]]["from"] == []


def test_an_unrecorded_ingested_file_is_not_claimed_as_a_contributor(
    tiny_store, tmp_path: Path
) -> None:
    """An assembly with a ``fasta_path`` override is ingested from a local file,
    and ``resolve_sources`` deliberately emits no remote source for it -- so
    nothing pins it in ``inputs.files``. Listing it under ``from`` would make the
    lock reference an input it does not contain."""
    digests = sorted(store_census.collection_census(tiny_store))
    provenance = {
        str(tmp_path / "downloads" / "override.fa"): {
            "collection_digest": digests[0], "n_sequences": 2,
        },
    }
    from_map = build_lock.provenance_from_map(tmp_path / "downloads", provenance)
    outputs = build_lock.store_outputs(tiny_store, from_map, known_paths=set())
    assert all(c["from"] == [] for c in outputs["collections"])


def test_build_lock_dict_never_emits_an_unvalidatable_lock(
    tiny_store, tmp_path: Path
) -> None:
    """The regression the dev manifest caught: a full build over an assembly
    with a ``fasta_path`` override wrote a lock its own validator rejected."""
    config = tmp_path / "sources.toml"
    config.write_text("")
    digests = sorted(store_census.collection_census(tiny_store))
    lock = build_lock.build_lock_dict(
        config_path=config,
        download_dir=tmp_path / "downloads",
        resolved_sources=[],           # the override emits no remote source
        collection_by_cachepath={
            str(tmp_path / "downloads" / "override.fa"): {
                "collection_digest": digests[0], "n_sequences": 2,
            },
        },
        store=tiny_store,
    )
    build_lock.validate_lock(lock)
    assert build_lock.files_by_collection(lock)[digests[0]] == []


# ------------------------------------------------------------- build metadata

class _FakeDistribution:
    def __init__(self, direct_url: str | None) -> None:
        self._direct_url = direct_url

    def read_text(self, name: str) -> str | None:
        return self._direct_url if name == "direct_url.json" else None


def test_gtars_built_from_a_commit_records_that_commit(monkeypatch) -> None:
    direct_url = json.dumps({
        "url": "https://github.com/example/gtars",
        "vcs_info": {"vcs": "git", "commit_id": "c" * 40,
                     "requested_revision": "c" * 40},
        "subdirectory": "gtars-python",
    })
    monkeypatch.setattr("importlib.metadata.distribution",
                        lambda name: _FakeDistribution(direct_url))
    assert build_lock._gtars_source() == {
        "url": "https://github.com/example/gtars", "vcs": "git",
        "commit": "c" * 40, "subdirectory": "gtars-python",
    }


def test_gtars_from_a_registry_wheel_records_no_source(monkeypatch) -> None:
    monkeypatch.setattr("importlib.metadata.distribution",
                        lambda name: _FakeDistribution(None))
    assert build_lock._gtars_source() is None
