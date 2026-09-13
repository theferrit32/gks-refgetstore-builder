"""Cache verification, and the target-selection rules around it.

The cache target is lock-driven and offline: the lock's ``sha256`` is the
authority, never the provider's checksum. A provider that republished is
``--manifest``'s business, and conflating the two is what made the old
``verify`` unable to run without a network.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import build_lock
import verify
from conftest import file_record, v4_lock

REL = "host/path/file.fa.gz"


def cached(cache_dir: Path, rel: str, payload: bytes, *, gzipped: bool = True) -> dict:
    """Write a cache file and return the lock record that pins it."""
    path = cache_dir / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = gzip.compress(payload) if gzipped else payload
    path.write_bytes(blob)
    return file_record(
        rel, sha256=hashlib.sha256(blob).hexdigest(), bytes=len(blob),
        present=True,
    )


def codes(report: verify.VerifyReport) -> list[str]:
    return sorted(f.code for f in report.findings)


def test_an_intact_cache_passes(tmp_path: Path) -> None:
    record = cached(tmp_path, REL, b">a\nACGT\n")
    report = verify.verify_cache(v4_lock([record]), tmp_path)
    assert report.findings == []
    assert not report.failed


def test_a_missing_file_is_an_error(tmp_path: Path) -> None:
    record = cached(tmp_path, REL, b">a\nACGT\n")
    (tmp_path / REL).unlink()
    report = verify.verify_cache(v4_lock([record]), tmp_path)
    assert codes(report) == ["cache.missing"]
    assert report.failed


def test_an_empty_file_counts_as_missing(tmp_path: Path) -> None:
    record = cached(tmp_path, REL, b">a\nACGT\n")
    (tmp_path / REL).write_bytes(b"")
    assert "cache.missing" in codes(verify.verify_cache(v4_lock([record]), tmp_path))


def test_changed_bytes_are_a_sha256_mismatch(tmp_path: Path) -> None:
    record = cached(tmp_path, REL, b">a\nACGT\n")
    (tmp_path / REL).write_bytes(gzip.compress(b">a\nTTTT\n"))
    report = verify.verify_cache(v4_lock([record]), tmp_path)
    assert "cache.sha256_mismatch" in codes(report)
    assert report.failed


def test_a_truncated_file_is_caught_by_size_and_gzip(tmp_path: Path) -> None:
    record = cached(tmp_path, REL, b">a\n" + b"ACGT" * 500 + b"\n")
    blob = (tmp_path / REL).read_bytes()
    (tmp_path / REL).write_bytes(blob[: len(blob) // 2])
    assert codes(verify.verify_cache(v4_lock([record]), tmp_path)) == [
        "cache.gzip_corrupt", "cache.sha256_mismatch", "cache.size_mismatch",
    ]


def test_gzip_corruption_that_preserves_length_is_still_caught(
    tmp_path: Path,
) -> None:
    """The CRC32 trailer is what makes this detectable at all -- the file is
    the right size and the right name."""
    record = cached(tmp_path, REL, b">a\n" + b"ACGT" * 500 + b"\n")
    blob = bytearray((tmp_path / REL).read_bytes())
    blob[len(blob) // 2] ^= 0xFF
    (tmp_path / REL).write_bytes(bytes(blob))
    found = codes(verify.verify_cache(v4_lock([record]), tmp_path))
    assert "cache.gzip_corrupt" in found
    assert "cache.size_mismatch" not in found


def test_no_hash_still_catches_absence_and_gzip_damage(tmp_path: Path) -> None:
    record = cached(tmp_path, REL, b">a\n" + b"ACGT" * 500 + b"\n")
    blob = bytearray((tmp_path / REL).read_bytes())
    blob[len(blob) // 2] ^= 0xFF
    (tmp_path / REL).write_bytes(bytes(blob))
    report = verify.verify_cache(v4_lock([record]), tmp_path, hash_files=False)
    assert codes(report) == ["cache.gzip_corrupt"]


def test_a_record_the_lock_marks_absent_is_skipped_not_failed(
    tmp_path: Path,
) -> None:
    """``present: false`` is a recorded fact about the build, not damage."""
    record = file_record(REL, present=False)
    report = verify.verify_cache(v4_lock([record]), tmp_path)
    assert report.findings == []
    assert any("skipped" in note for note in report.notes)


def test_a_non_gzip_file_is_not_handed_to_gzip(tmp_path: Path) -> None:
    record = cached(tmp_path, "host/report.txt", b"plain text\n", gzipped=False)
    assert verify.verify_cache(v4_lock([record]), tmp_path).findings == []


def test_limit_checks_only_the_first_n(tmp_path: Path) -> None:
    records = [cached(tmp_path, f"host/{i}.fa.gz", b">a\nACGT\n")
               for i in range(4)]
    for i in range(4):
        (tmp_path / f"host/{i}.fa.gz").unlink()
    report = verify.verify_cache(v4_lock(records), tmp_path, limit=2)
    assert len(report.findings) == 2


def test_verify_cache_never_touches_the_network(tmp_path: Path) -> None:
    """Enforced by the autouse block_network fixture; asserted here so the
    offline-by-default promise is a test, not a docstring."""
    record = cached(tmp_path, REL, b">a\nACGT\n")
    verify.verify_cache(v4_lock([record]), tmp_path)


# ------------------------------------------------------- target selection

def args_for(tmp_path: Path, lock_path: Path, **overrides):
    base = {
        "cache": False, "store": False, "remote": False, "manifest": False,
        "all": False, "deep": False, "no_hash": False, "limit": None,
        "jobs": 1, "timeout": 1, "known_bad": None, "lock": lock_path,
        "cache_dir": tmp_path, "store_dir": tmp_path / "store",
        "config": tmp_path / "sources.toml",
    }
    base.update(overrides)
    return type("Args", (), base)()


def test_no_scope_flag_means_cache_and_store(tmp_path: Path, monkeypatch,
                                             tiny_store_dir, tiny_store) -> None:
    from conftest import lock_for_store

    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, lock_for_store(tiny_store, tiny_store_dir))
    called: list[str] = []
    for name in ("verify_cache", "verify_store", "verify_remote",
                 "verify_manifest"):
        monkeypatch.setattr(
            verify, name,
            lambda *a, _n=name, **k: called.append(_n) or verify.VerifyReport(),
        )
    verify.run_verify(args_for(tmp_path, lock_path, store_dir=tiny_store_dir))
    assert sorted(called) == ["verify_cache", "verify_store"]


def test_all_selects_every_target(tmp_path: Path, monkeypatch,
                                  tiny_store_dir, tiny_store) -> None:
    from conftest import lock_for_store

    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, lock_for_store(tiny_store, tiny_store_dir))
    called: list[str] = []
    for name in ("verify_cache", "verify_store", "verify_remote",
                 "verify_manifest"):
        monkeypatch.setattr(
            verify, name,
            lambda *a, _n=name, **k: called.append(_n) or verify.VerifyReport(),
        )
    verify.run_verify(args_for(tmp_path, lock_path, all=True,
                               store_dir=tiny_store_dir))
    assert sorted(called) == ["verify_cache", "verify_manifest",
                              "verify_remote", "verify_store"]


def test_exit_code_is_one_only_when_something_errored(tmp_path: Path) -> None:
    good = cached(tmp_path, REL, b">a\nACGT\n")
    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, v4_lock([good]))
    assert verify.run_verify(args_for(tmp_path, lock_path, cache=True)) == 0
    (tmp_path / REL).unlink()
    assert verify.run_verify(args_for(tmp_path, lock_path, cache=True)) == 1


def test_warnings_alone_do_not_fail(tmp_path: Path) -> None:
    report = verify.VerifyReport()
    report.add(verify.Finding("remote.size_mismatch", verify.WARN, "x"))
    report.add(verify.Finding("store.seq_digest_fixed", verify.INFO, "y"))
    assert not report.failed
