"""The case-sensitivity preflight for commands that write store payloads."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from gks_refgetstore import build_store, fs_checks, repair, store_sync

MISSING_CONFIG = Path("no-such-dir/sources.toml")


def _observed_case_sensitive(directory: Path) -> bool:
    """Independent answer: create a mixed-case name, look up its lowercase."""
    probe = directory / "ObservedCaseProbe"
    probe.write_text("x")
    try:
        return not (directory / "observedcaseprobe").exists()
    finally:
        probe.unlink()


# ------------------------------------------------------------------- the probe

@pytest.mark.real_case_probe
def test_probe_agrees_with_direct_observation(tmp_path: Path) -> None:
    assert fs_checks.is_case_sensitive(tmp_path) == _observed_case_sensitive(tmp_path)


@pytest.mark.real_case_probe
def test_probe_leaves_nothing_behind(tmp_path: Path) -> None:
    fs_checks.is_case_sensitive(tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.real_case_probe
def test_probe_of_a_missing_dir_checks_its_ancestor_without_creating_it(
    tmp_path: Path,
) -> None:
    target = tmp_path / "store" / "not-yet"
    assert fs_checks.is_case_sensitive(target) == _observed_case_sensitive(tmp_path)
    assert not (tmp_path / "store").exists()


# -------------------------------------------------------------- the refusal

def _case_insensitive(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fs_checks, "is_case_sensitive", lambda directory: False)


def test_a_case_insensitive_store_dir_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    _case_insensitive(monkeypatch)
    with pytest.raises(fs_checks.CaseInsensitiveFilesystemError,
                       match=fs_checks.ALLOW_FLAG):
        fs_checks.require_case_sensitive(tmp_path)


def test_the_allow_flag_skips_the_refusal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    _case_insensitive(monkeypatch)
    fs_checks.require_case_sensitive(tmp_path, allow=True)


def test_build_refuses_before_reading_its_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    _case_insensitive(monkeypatch)
    args = argparse.Namespace(store_dir=tmp_path / "store", config=MISSING_CONFIG,
                              allow_case_insensitive_fs=False)
    with pytest.raises(SystemExit, match="case-insensitive filesystem"):
        build_store.run_build(args)


def test_sync_apply_refuses_before_reading_its_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    _case_insensitive(monkeypatch)
    args = argparse.Namespace(store_dir=tmp_path / "store", config=MISSING_CONFIG,
                              apply=True, allow_case_insensitive_fs=False)
    with pytest.raises(SystemExit, match="case-insensitive filesystem"):
        store_sync.run_sync(args)


def test_sync_dry_run_does_not_probe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A dry run only reads, so it reaches the config instead of refusing."""
    _case_insensitive(monkeypatch)
    args = argparse.Namespace(store_dir=tmp_path / "store", config=MISSING_CONFIG,
                              apply=False, allow_case_insensitive_fs=False)
    with pytest.raises(FileNotFoundError):
        store_sync.run_sync(args)


def test_repair_apply_refuses_before_reading_its_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    _case_insensitive(monkeypatch)
    args = argparse.Namespace(store_dir=tmp_path / "store",
                              lock=tmp_path / "missing.lock.json",
                              cache=False, store=True, apply=True,
                              allow_case_insensitive_fs=False)
    with pytest.raises(SystemExit, match="case-insensitive filesystem"):
        repair.run_repair(args)
