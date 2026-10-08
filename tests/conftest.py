"""Shared test policy.

No test may reach a provider. Source resolution fetches checksum manifests and
file listings from NCBI and Ensembl, so a test that resolves for real depends on
what upstream published today: one failed exactly that way when a ninth
RefSeqGene shard appeared. Drive resolution from the build lock
(``apply_locked_sources``) or from fixtures instead, and inject ``fetch_text``
where a manifest is genuinely under test.

No test may depend on the case sensitivity of the filesystem it runs on.
Commands that write a store refuse a case-insensitive one, and ``tmp_path`` is
case-insensitive on a default macOS install but not on Linux. The probe is
stubbed as case-sensitive everywhere except tests marked ``real_case_probe``.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

import pytest


def _blocked(*args, **kwargs):
    url = args[0] if args else kwargs.get("url", "")
    raise AssertionError(
        f"test attempted a network request ({url!r}). Resolve from the lock or "
        "a fixture; see tests/conftest.py"
    )


@pytest.fixture(autouse=True)
def block_network(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(urllib.request, "urlopen", _blocked)
    monkeypatch.setattr(urllib.request, "urlretrieve", _blocked)


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "real_case_probe: run the filesystem case-sensitivity probe for real",
    )


@pytest.fixture(autouse=True)
def case_sensitive_store_dirs(request: pytest.FixtureRequest,
                              monkeypatch: pytest.MonkeyPatch) -> None:
    if request.node.get_closest_marker("real_case_probe") is None:
        monkeypatch.setattr("gks_refgetstore.fs_checks.is_case_sensitive",
                            lambda directory: True)


# ---------------------------------------------------------------------- locks

# Building a valid /4 lock by hand is fiddly -- validate_lock enforces the
# ingest-spec digest, the collections root, and the count agreement -- so the
# construction lives here once rather than in every test module.


def file_record(cache_path: str, **overrides) -> dict:
    """A complete ``inputs.files`` record with every field present."""
    from gks_refgetstore import build_lock

    record = {
        "kind": "seqset",
        "owner": "owner",
        "url": "https://example.test/" + cache_path,
        "cache_path": cache_path,
        "mutable": False,
        "provider_checksum": None,
        "provider_checksum_algorithm": None,
        "provider_checksum_blocks": None,
        "checksum_url": None,
        "file_class": None,
        "bytes": None,
        "sha256": None,
        "present": True,
        "ingest_spec": None,
        "ingest_spec_sha256": None,
    }
    record.update(overrides)
    # Keep the pair consistent unless a test is deliberately breaking it.
    if "ingest_spec_sha256" not in overrides:
        record["ingest_spec_sha256"] = build_lock.canonical_spec_sha256(
            record["ingest_spec"]
        )
    return record


def v4_lock(
    files: list[dict] | None = None,
    collections: list[dict] | None = None,
    *,
    sources_toml_sha256: str | None = None,
    n_sequences: int = 0,
    sequences_root: str | None = None,
) -> dict:
    """A schema-/4 lock that satisfies ``validate_lock``."""
    from gks_refgetstore import build_lock
    from gks_refgetstore import store_census

    files = files or []
    collections = collections or []
    for collection in collections:
        collection.setdefault("from", [])
        collection.setdefault("n_sequences", 0)
    digests = [c["digest"] for c in collections]
    return {
        "schema": build_lock.SCHEMA,
        "build": {"timestamp_utc": "2026-01-01T00:00:00+00:00",
                  "git": {}, "gtars_version": "0.9.2"},
        "inputs": {
            "sources_toml_sha256": sources_toml_sha256,
            "files": files,
        },
        "outputs": {
            "n_sequences": n_sequences,
            "n_collections": len(collections),
            "sequences_root": sequences_root or store_census.digest_root([]),
            "collections_root": store_census.digest_root(digests),
            "collections": collections,
        },
    }


def lock_for_store(store, store_dir, files=None, from_map=None) -> dict:
    """A lock whose ``outputs`` is a census of a real store."""
    from gks_refgetstore import build_lock

    lock = v4_lock(files or [])
    lock["outputs"] = build_lock.store_outputs(store, from_map or {})
    return lock


# --------------------------------------------------------------------- store

# Two collections sharing one sequence, so removal, orphan detection, and the
# from-list's many-files-to-one-collection case all have something real to act
# on. Deliberately a *real* gtars store rather than a mock: the behaviours
# verify depends on -- lazy collection stubs, on-disk .seq layout, encoder
# round-tripping -- live in gtars, and a mock would assert our beliefs about
# them rather than the behaviour.
TINY_COLLECTIONS = {
    "alpha.fa": [
        ("chr_a", "ACGTACGTAAGGTTCCAA"),
        ("chr_shared", "TTTTGGGGCCCCAAAATT"),
    ],
    "beta.fa": [
        ("chr_shared", "TTTTGGGGCCCCAAAATT"),
        ("chr_b", "GATTACAGATTACAGATT"),
    ],
}


def write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text(
        "".join(f">{name}\n{seq}\n" for name, seq in records), encoding="utf-8"
    )
    return path


@pytest.fixture
def tiny_store_dir(tmp_path: Path) -> Path:
    """A real two-collection store on disk. Returns its directory."""
    from gtars.refget import RefgetStore

    store_dir = tmp_path / "store.tiny"
    store = RefgetStore.on_disk(str(store_dir))
    store.set_quiet(True)
    for filename, records in TINY_COLLECTIONS.items():
        write_fasta(tmp_path / filename, records)
        store.add_sequence_collection_from_fasta(str(tmp_path / filename))
    store.write()
    return store_dir


@pytest.fixture
def tiny_store(tiny_store_dir: Path):
    """The tiny store, reopened read-only the way every consumer opens it."""
    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(str(tiny_store_dir))
    store.set_quiet(True)
    return store
