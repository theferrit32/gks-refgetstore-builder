"""``status`` describes; ``verify`` judges.

The one property worth locking down hard is the exit code. ``status`` exists to
be run casually and read, so it must be safe in a shell with ``set -e`` and
useless as a CI gate -- if it ever failed on a stale lock, people would stop
running it, and the four-way picture it gives is the thing that makes a stale
lock legible in the first place.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import build_lock
import verify
from conftest import file_record, lock_for_store, v4_lock


def args_for(tmp_path: Path, lock_path: Path, store_dir: Path, **overrides):
    base = {
        "config": tmp_path / "sources.toml", "cache_dir": tmp_path,
        "store_dir": store_dir, "lock": lock_path, "offline": True,
    }
    base.update(overrides)
    return type("Args", (), base)()


@pytest.fixture
def lock_path(tmp_path: Path, tiny_store, tiny_store_dir) -> Path:
    path = tmp_path / "lock.json"
    build_lock.write_lock(path, lock_for_store(tiny_store, tiny_store_dir))
    (tmp_path / "sources.toml").write_text("# empty\n")
    return path


def test_status_exits_zero_when_everything_agrees(
    tmp_path: Path, lock_path: Path, tiny_store_dir, capsys,
) -> None:
    assert verify.run_status(args_for(tmp_path, lock_path, tiny_store_dir)) == 0
    assert "roots agree" in capsys.readouterr().out


def test_status_exits_zero_on_an_edited_manifest(
    tmp_path: Path, lock_path: Path, tiny_store_dir, capsys,
) -> None:
    (tmp_path / "sources.toml").write_text("# changed\n")
    assert verify.run_status(args_for(tmp_path, lock_path, tiny_store_dir)) == 0
    assert "EDITED" in capsys.readouterr().out


def test_status_exits_zero_on_a_damaged_store(
    tmp_path: Path, lock_path: Path, tiny_store_dir, capsys,
) -> None:
    """The case that would tempt a nonzero exit. A store that no longer matches
    the lock is exactly what status is for describing -- and exactly what
    `verify --store` is for failing on."""
    lock = build_lock.load_lock(lock_path)
    lock["outputs"]["n_sequences"] = 999
    lock["outputs"]["sequences_root"] = "sha256:" + "0" * 64
    build_lock.write_lock(lock_path, lock)

    assert verify.run_status(args_for(tmp_path, lock_path, tiny_store_dir)) == 0
    out = capsys.readouterr().out
    assert "store.n_sequences_mismatch" in out
    assert "store.sequences_root_mismatch" in out

    # ... and verify does fail on the same state. The two commands disagreeing
    # about severity is the design, not an inconsistency.
    verify_args = type("Args", (), {
        "cache": False, "store": True, "remote": False, "manifest": False,
        "all": False, "deep": False, "no_hash": False, "limit": None,
        "jobs": 1, "timeout": 1, "known_bad": None, "lock": lock_path,
        "cache_dir": tmp_path, "store_dir": tiny_store_dir,
        "config": tmp_path / "sources.toml",
    })()
    assert verify.run_verify(verify_args) == 1


def test_status_exits_zero_when_the_store_is_absent(
    tmp_path: Path, lock_path: Path, capsys,
) -> None:
    assert verify.run_status(
        args_for(tmp_path, lock_path, tmp_path / "no-such-store")
    ) == 0
    assert "MISSING" in capsys.readouterr().out


def test_status_exits_zero_when_the_manifest_is_absent(
    tmp_path: Path, lock_path: Path, tiny_store_dir, capsys,
) -> None:
    (tmp_path / "sources.toml").unlink()
    assert verify.run_status(args_for(tmp_path, lock_path, tiny_store_dir)) == 0
    assert "MISSING" in capsys.readouterr().out


def test_offline_status_makes_no_network_request(
    tmp_path: Path, lock_path: Path, tiny_store_dir, capsys,
) -> None:
    """The autouse block_network fixture would raise; asserting the skip
    message keeps the promise explicit rather than incidental."""
    verify.run_status(args_for(tmp_path, lock_path, tiny_store_dir))
    assert "skipped: --offline" in capsys.readouterr().out


def test_status_reports_cache_presence_without_hashing(
    tmp_path: Path, tiny_store, tiny_store_dir, capsys,
) -> None:
    lock = lock_for_store(tiny_store, tiny_store_dir)
    lock["inputs"]["files"] = [
        file_record("host/here.gz", sha256="0" * 64),
        file_record("host/gone.gz", sha256="0" * 64),
    ]
    # Present, but with bytes that do not match the locked SHA-256. status
    # counts it as present anyway -- hashing 106 GB is `verify --cache`'s job.
    (tmp_path / "host").mkdir()
    (tmp_path / "host/here.gz").write_bytes(b"wrong bytes")
    path = tmp_path / "lock.json"
    build_lock.write_lock(path, lock)
    (tmp_path / "sources.toml").write_text("")

    assert verify.run_status(args_for(tmp_path, path, tiny_store_dir)) == 0
    assert "1/2 locked file(s) present" in capsys.readouterr().out


def test_status_refuses_a_lock_it_cannot_read(tmp_path: Path,
                                              tiny_store_dir) -> None:
    """The one thing status will not do is describe a lock it does not
    understand: a /2 lock reads as an empty v4 one, and "0 files, 0
    collections" is a confident lie."""
    path = tmp_path / "lock.json"
    build_lock.write_lock(path, {"schema": "gks-refgetstore-build-lock/2",
                                 "sources": []})
    with pytest.raises(build_lock.LockError):
        verify.run_status(args_for(tmp_path, path, tiny_store_dir))


def test_status_summarizes_a_sync_plan_when_online(
    tmp_path: Path, lock_path: Path, tiny_store_dir, monkeypatch, capsys,
) -> None:
    import build_store

    monkeypatch.setattr(build_store, "load_config", lambda _p: ([], []))
    monkeypatch.setattr(build_store, "resolve_sources", lambda *_a: [])
    assert verify.run_status(
        args_for(tmp_path, lock_path, tiny_store_dir, offline=False)
    ) == 0
    out = capsys.readouterr().out
    assert "upstream:" in out and "sync would:" in out


def test_empty_lock_and_empty_store_still_describes(tmp_path: Path,
                                                    capsys) -> None:
    path = tmp_path / "lock.json"
    build_lock.write_lock(path, v4_lock())
    (tmp_path / "sources.toml").write_text("")
    assert verify.run_status(
        args_for(tmp_path, path, tmp_path / "absent-store")
    ) == 0
    assert "0 input file(s)" in capsys.readouterr().out
