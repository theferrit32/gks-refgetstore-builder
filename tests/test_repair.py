"""Repair: what it dispatches on, what it refuses, and the order it works in.

Three properties matter more than the mechanics:

1. Repair restores the **locked** state. It accepts re-downloaded bytes only
   when they match the lock's SHA-256, so it can never quietly re-baseline a
   drifted upstream.
2. Cache repair precedes store repair and aborts the run on failure. Store
   removal persists immediately and gtars offers no rollback, so a source file
   discovered missing *after* the removal leaves the store worse than it began.
3. A defect re-ingestion would reproduce -- the gtars encoder round-trip -- is
   refused, not attempted.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from conftest import file_record, v4_lock
from gks_refgetstore import build_lock, repair, verify

REL_A = "host/a.fa.gz"
REL_B = "host/b.fa.gz"


def finding(code: str, subject: str, severity: str = verify.ERROR):
    return verify.Finding(code, severity, subject)


def lock_two_files_one_collection() -> dict:
    return v4_lock(
        [file_record(REL_A, sha256="a" * 64),
         file_record(REL_B, sha256="b" * 64)],
        [{"digest": "SHARED", "n_sequences": 2, "from": [REL_A, REL_B]}],
    )


# ------------------------------------------------------------ dispatch table

def test_cache_findings_become_refetches() -> None:
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(lock, [
        finding("cache.missing", REL_A),
        finding("cache.sha256_mismatch", REL_B),
    ])
    assert plan.cache_files == [REL_A, REL_B]
    assert plan.store_files == []


def test_a_repeated_cache_finding_is_fetched_once() -> None:
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(lock, [
        finding("cache.sha256_mismatch", REL_A),
        finding("cache.gzip_corrupt", REL_A),
    ])
    assert plan.cache_files == [REL_A]


def test_a_lost_collection_re_ingests_every_contributor() -> None:
    """The correction the ``from`` list made possible.

    Removing a collection because one contributor's damage was noticed also
    removes the sequences the *other* contributor published, so both must be
    re-ingested. Re-ingesting only the file that was noticed silently loses the
    rest.
    """
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(lock, [finding("store.collection_missing", "SHARED")])
    assert plan.store_collections == ["SHARED"]
    assert plan.store_files == [REL_A, REL_B]


def test_a_collection_size_mismatch_is_repairable() -> None:
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(
        lock, [finding("store.collection_size_mismatch", "SHARED")]
    )
    assert plan.store_files == [REL_A, REL_B]


def test_assembly_owners_are_reported_rather_than_re_ingested() -> None:
    """``apply_plan`` ingests seqsets. Assembly ingestion carries its own
    deferred-report sequencing, which ``build`` owns, so repair names the
    command instead of half-doing it."""
    lock = v4_lock(
        [file_record(REL_A, kind="assembly_fasta", owner="GRCh37")],
        [{"digest": "ASM", "n_sequences": 1, "from": [REL_A]}],
    )
    plan = repair.plan_repair(lock, [finding("store.collection_missing", "ASM")])
    assert plan.assembly_owners == {"GRCh37"}


def test_an_encoder_round_trip_defect_is_unrepairable() -> None:
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(
        lock, [finding("store.seq_digest_mismatch", "SOMEDIGEST")]
    )
    assert plan.unrepairable == [("store.seq_digest_mismatch", "SOMEDIGEST")]
    assert not plan.touched
    assert "leaves an existing payload in place" in (
        repair.UNREPAIRABLE["store.seq_digest_mismatch"])


def test_an_unattributable_root_mismatch_points_at_sync() -> None:
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(
        lock, [finding("store.sequences_root_mismatch", "sequences_root")]
    )
    assert not plan.touched
    assert "sync" in repair.UNREPAIRABLE["store.sequences_root_mismatch"]


def test_warnings_repair_does_not_handle_are_listed_not_acted_on() -> None:
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(
        lock, [finding("store.orphan_sequence", "D1", verify.WARN)]
    )
    assert plan.ignored == [("store.orphan_sequence", "D1")]
    assert not plan.touched


def test_informational_findings_are_ignored_entirely() -> None:
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(
        lock, [finding("store.seq_digest_fixed", "D1", verify.INFO)]
    )
    assert plan.ignored == [] and not plan.touched


# ------------------------------------------------------------- accept policy

def test_repair_accepts_only_the_locked_bytes(tmp_path: Path, monkeypatch) -> None:
    """The rule that keeps repair from becoming an accidental re-baseline."""
    payload = b">a\nACGT\n"
    locked_sha = hashlib.sha256(payload).hexdigest()
    lock = v4_lock([file_record(REL_A, sha256=locked_sha)])
    (tmp_path / REL_A).parent.mkdir(parents=True)

    served = {"bytes": b"something else entirely"}

    def fake_download(url, target, timeout, retries):
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(served["bytes"])
        return len(served["bytes"])

    from gks_refgetstore import fetch_sources
    monkeypatch.setattr(fetch_sources, "download", fake_download)
    args = type("Args", (), {
        "cache_dir": tmp_path, "timeout": 1, "retries": 1, "jobs": 1,
        "min_free_gb": 0.0,
    })()

    plan = repair.plan_repair(lock, [finding("cache.missing", REL_A)])
    assert repair.repair_cache(plan, lock, args) is False
    # Bytes that were not the locked ones must not be left in the cache to be
    # mistaken for good ones later.
    assert not (tmp_path / REL_A).exists()

    served["bytes"] = payload
    assert repair.repair_cache(plan, lock, args) is True
    assert (tmp_path / REL_A).read_bytes() == payload


def test_a_file_the_lock_pins_no_sha256_for_cannot_be_repaired(
    tmp_path: Path,
) -> None:
    lock = v4_lock([file_record(REL_A, sha256=None)])
    plan = repair.plan_repair(lock, [finding("cache.missing", REL_A)])
    args = type("Args", (), {
        "cache_dir": tmp_path, "timeout": 1, "retries": 1, "jobs": 1,
        "min_free_gb": 0.0,
    })()
    assert repair.repair_cache(plan, lock, args) is False


# ------------------------------------------------------------------ ordering

def test_cache_repair_runs_first_and_a_failure_stops_the_store(
    tmp_path: Path, monkeypatch, tiny_store, tiny_store_dir,
) -> None:
    from conftest import lock_for_store

    order: list[str] = []
    lock = lock_for_store(tiny_store, tiny_store_dir)
    digest = lock["outputs"]["collections"][0]["digest"]
    lock["inputs"]["files"] = [file_record(REL_A, sha256="a" * 64)]
    lock["outputs"]["collections"][0]["from"] = [REL_A]
    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, lock)

    def cache_repair(*_a, **_k):
        order.append("cache")
        return False  # the source is gone upstream

    def store_repair(*_a, **_k):
        order.append("store")
        return True

    monkeypatch.setattr(repair, "repair_cache", cache_repair)
    monkeypatch.setattr(repair, "repair_store", store_repair)
    monkeypatch.setattr(
        repair.verify, "verify_cache",
        lambda *a, **k: _report(finding("cache.missing", REL_A)),
    )
    monkeypatch.setattr(
        repair.verify, "verify_store",
        lambda *a, **k: _report(finding("store.collection_missing", digest)),
    )
    args = type("Args", (), {
        "cache": False, "store": False, "apply": True, "lock": lock_path,
        "cache_dir": tmp_path, "store_dir": tiny_store_dir,
        "config": tmp_path / "sources.toml", "timeout": 1, "retries": 1,
        "jobs": 1, "min_free_gb": 0.0, "ingest_jobs": 1, "filter_jobs": 1,
    })()
    assert repair.run_repair(args) == 1
    assert order == ["cache"], "the store must not be touched after a cache failure"


def _report(*findings) -> verify.VerifyReport:
    report = verify.VerifyReport()
    report.add(*findings)
    return report


def test_dry_run_mutates_nothing(tmp_path: Path, monkeypatch, tiny_store,
                                 tiny_store_dir) -> None:
    from conftest import lock_for_store

    lock = lock_for_store(tiny_store, tiny_store_dir)
    lock["inputs"]["files"] = [file_record(REL_A, sha256="a" * 64)]
    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, lock)
    monkeypatch.setattr(
        repair.verify, "verify_cache",
        lambda *a, **k: _report(finding("cache.missing", REL_A)),
    )
    monkeypatch.setattr(repair.verify, "verify_store",
                        lambda *a, **k: verify.VerifyReport())
    for name in ("repair_cache", "repair_store"):
        monkeypatch.setattr(
            repair, name,
            lambda *a, **k: pytest.fail("dry run attempted a repair"),
        )
    args = type("Args", (), {
        "cache": True, "store": False, "apply": False, "lock": lock_path,
        "cache_dir": tmp_path, "store_dir": tiny_store_dir,
        "config": tmp_path / "sources.toml", "timeout": 1, "retries": 1,
        "jobs": 1, "min_free_gb": 0.0, "ingest_jobs": 1, "filter_jobs": 1,
    })()
    assert repair.run_repair(args) == 0


def test_a_clean_store_and_cache_need_no_repair(tmp_path: Path, tiny_store,
                                                tiny_store_dir) -> None:
    from conftest import lock_for_store

    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, lock_for_store(tiny_store, tiny_store_dir))
    args = type("Args", (), {
        "cache": False, "store": False, "apply": True, "lock": lock_path,
        "cache_dir": tmp_path, "store_dir": tiny_store_dir,
        "config": tmp_path / "sources.toml", "timeout": 1, "retries": 1,
        "jobs": 1, "min_free_gb": 0.0, "ingest_jobs": 1, "filter_jobs": 1,
    })()
    assert repair.run_repair(args) == 0


def test_repair_refuses_to_mutate_when_a_source_is_absent_from_the_cache(
    tmp_path: Path, tiny_store, tiny_store_dir,
) -> None:
    from conftest import lock_for_store

    lock = lock_for_store(tiny_store, tiny_store_dir)
    lock["inputs"]["files"] = [file_record(REL_A, sha256="a" * 64)]
    lock["outputs"]["collections"][0]["from"] = [REL_A]
    digest = lock["outputs"]["collections"][0]["digest"]
    (tmp_path / "sources.toml").write_text("")
    plan = repair.plan_repair(lock, [finding("store.collection_missing", digest)])
    args = type("Args", (), {
        "cache_dir": tmp_path, "store_dir": tiny_store_dir,
        "config": tmp_path / "sources.toml", "ingest_jobs": 1, "filter_jobs": 1,
    })()
    assert repair.repair_store(plan, lock, args) is False


def test_a_missing_payload_without_a_store_is_reported_not_silently_dropped(
) -> None:
    """Planning a seq_file_missing means asking the store which collection
    publishes that digest. With no store, say so rather than plan nothing."""
    lock = lock_two_files_one_collection()
    plan = repair.plan_repair(lock, [finding("store.seq_file_missing", "D1")])
    assert plan.ignored == [("store.seq_file_missing", "D1")]
    assert not plan.touched


def test_a_missing_payload_resolves_to_its_collections_contributors(
    tiny_store, tiny_store_dir,
) -> None:
    from conftest import lock_for_store

    from gks_refgetstore import store_census

    digests = sorted(store_census.collection_census(tiny_store))
    victim = sorted(store_census.collection_members(tiny_store, digests[0]))[0]
    lock = lock_for_store(
        tiny_store, tiny_store_dir,
        files=[file_record(REL_A), file_record(REL_B)],
        from_map={digests[0]: {REL_A, REL_B}},
    )
    plan = repair.plan_repair(
        lock, [finding("store.seq_file_missing", victim)], tiny_store_dir
    )
    assert digests[0] in plan.store_collections
    assert plan.store_files == [REL_A, REL_B]
