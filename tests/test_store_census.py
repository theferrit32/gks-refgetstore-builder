"""The read-only store primitives every lock/verify/repair path shares."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

import store_census as sc


def test_digest_root_is_order_and_duplicate_independent():
    digests = ["bbb", "aaa", "ccc"]
    root = sc.digest_root(digests)
    assert root == sc.digest_root(reversed(digests))
    assert root == sc.digest_root(digests + digests)


def test_digest_root_matches_the_documented_construction():
    # Spelled out longhand: the lock's roots are a cross-implementation
    # contract, so the definition is pinned by an independent computation
    # rather than by whatever digest_root happens to return.
    hasher = hashlib.sha256()
    for digest in ("aaa", "bbb"):
        hasher.update(digest.encode("ascii"))
        hasher.update(b"\n")
    assert sc.digest_root(["bbb", "aaa"]) == "sha256:" + hasher.hexdigest()


def test_digest_root_ignores_the_sq_prefix():
    assert sc.digest_root(["SQ.aaa"]) == sc.digest_root(["aaa"])


def test_digest_root_of_nothing_is_defined():
    assert sc.digest_root([]) == "sha256:" + hashlib.sha256().hexdigest()


def test_single_digest_change_flips_the_root():
    assert sc.digest_root(["aaa", "bbb"]) != sc.digest_root(["aaa", "bbc"])
    assert sc.digest_root(["aaa", "bbb"]) != sc.digest_root(["aaa"])


def test_census_counts_match_the_store(tiny_store):
    census = sc.collection_census(tiny_store)
    n_sequences, n_collections = sc.store_counts(tiny_store)
    assert len(census) == n_collections == 2
    assert all(n == 2 for n in census.values())
    # Three distinct sequences across two 2-sequence collections: the shared
    # one is stored once.
    assert n_sequences == 3
    assert sum(census.values()) == 4


def test_store_counts_are_integers(tiny_store):
    # store.stats() hands back strings; the lock records numbers.
    assert all(isinstance(v, int) for v in sc.store_counts(tiny_store))


def test_members_union_equals_the_sequence_index(tiny_store):
    indexed = sc.sequence_digests(tiny_store)
    members = {
        digest
        for collection in sc.collection_census(tiny_store)
        for digest in sc.collection_members(tiny_store, collection)
    }
    assert members == indexed


def test_shared_sequence_belongs_to_both_collections(tiny_store):
    census = sc.collection_census(tiny_store)
    memberships = [set(sc.collection_members(tiny_store, c)) for c in census]
    assert len(set.intersection(*memberships)) == 1


def test_named_members_pair_names_with_digests(tiny_store):
    census = sc.collection_census(tiny_store)
    names = {
        name
        for collection in census
        for name, _ in sc.collection_named_members(tiny_store, collection)
    }
    assert names == {"chr_a", "chr_b", "chr_shared"}


def test_on_disk_digests_match_the_index(tiny_store, tiny_store_dir):
    assert sc.on_disk_sequence_digests(tiny_store_dir) == sc.sequence_digests(tiny_store)


def test_on_disk_digests_of_a_missing_store_is_empty(tmp_path: Path):
    assert sc.on_disk_sequence_digests(tmp_path / "nope") == set()


def test_seq_path_shards_on_the_digest_prefix(tiny_store, tiny_store_dir):
    for digest in sc.sequence_digests(tiny_store):
        path = sc.seq_path(tiny_store_dir, digest)
        assert path.exists()
        assert path.parent.name == digest[:2]
        assert sc.seq_path(tiny_store_dir, "SQ." + digest) == path


def test_redigest_reproduces_every_stored_digest(tiny_store):
    # The whole premise of --deep: for an undamaged store with sequences the
    # encoder handles losslessly, the key and the recomputed digest agree.
    for digest in sc.sequence_digests(tiny_store):
        assert sc.redigest_sequence(tiny_store, digest) == digest


def test_redigest_accepts_a_prefixed_digest(tiny_store):
    digest = next(iter(sc.sequence_digests(tiny_store)))
    assert sc.redigest_sequence(tiny_store, "SQ." + digest) == digest


def test_redigest_sees_a_corrupted_payload(tiny_store_dir):
    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(str(tiny_store_dir))
    store.set_quiet(True)
    digest = sorted(sc.sequence_digests(store))[0]
    path = sc.seq_path(tiny_store_dir, digest)
    payload = bytearray(path.read_bytes())
    payload[-1] ^= 0xFF
    path.write_bytes(bytes(payload))

    reopened = RefgetStore.open_local(str(tiny_store_dir))
    reopened.set_quiet(True)
    assert sc.redigest_sequence(reopened, digest) != digest


@pytest.mark.parametrize("value", ["SQ.abc", "abc"])
def test_strip_sq_is_idempotent(value):
    assert sc.strip_sq(sc.strip_sq(value)) == "abc"
