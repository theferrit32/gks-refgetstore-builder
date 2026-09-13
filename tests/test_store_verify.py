"""The store ladder, against a real gtars store in ``tmp_path``.

The claim these tests exist to defend is that **the three levels are
independent**. Each induced failure must fire at exactly one level and leave the
others clean, because that is what makes the ladder worth climbing: L0 in seven
seconds is only useful if a clean L0 means something specific.

| induced                | L0    | L1                       | L2                        |
|------------------------|-------|--------------------------|---------------------------|
| unlink a ``.seq``      | clean | ``store.seq_file_missing``| —                        |
| corrupt a ``.seq``     | clean | clean                    | ``store.seq_digest_mismatch`` |
| add a sequence         | root mismatch | —                | —                         |

A real store, not a double: lazy collection stubs, the on-disk shard layout and
the encoder's round-trip behaviour are gtars' properties, and a mock would
assert our beliefs about them instead.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import build_lock
import store_census
import verify
from conftest import lock_for_store, write_fasta


@pytest.fixture
def locked(tiny_store, tiny_store_dir):
    """A lock that exactly describes the tiny store."""
    return lock_for_store(tiny_store, tiny_store_dir)


def codes(report: verify.VerifyReport) -> list[str]:
    return sorted(f.code for f in report.findings)


def reopen(store_dir: Path):
    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(str(store_dir))
    store.set_quiet(True)
    return store


# ------------------------------------------------------------------ baseline

def test_an_undamaged_store_passes_every_level(locked, tiny_store_dir) -> None:
    report = verify.verify_store(locked, tiny_store_dir, deep=True)
    assert report.findings == []
    assert not report.failed


def test_l2_needs_no_lock_at_all(tiny_store) -> None:
    """L2 compares the store against itself, which is why it is the one check
    that still means something when the lock is missing or wrong."""
    report = verify.verify_store_l2(tiny_store, {})
    assert report.findings == []


# ------------------------------------------------- L0: the roots and counts

def test_adding_a_sequence_flips_the_sequences_root(
    locked, tiny_store_dir, tmp_path: Path
) -> None:
    from gtars.refget import RefgetStore

    store = RefgetStore.on_disk(str(tiny_store_dir))
    store.set_quiet(True)
    write_fasta(tmp_path / "gamma.fa", [("chr_c", "CCCCAAAAGGGGTTTTAA")])
    store.add_sequence_collection_from_fasta(str(tmp_path / "gamma.fa"))
    store.write()

    report = verify.verify_store_l0(locked, reopen(tiny_store_dir))
    assert codes(report) == [
        "store.collections_root_mismatch", "store.n_collections_mismatch",
        "store.n_sequences_mismatch", "store.sequences_root_mismatch",
    ]
    assert report.failed


def test_a_root_mismatch_reports_both_sides(locked, tiny_store_dir) -> None:
    locked["outputs"]["sequences_root"] = store_census.digest_root(["nope"])
    report = verify.verify_store_l0(locked, reopen(tiny_store_dir))
    finding = next(f for f in report.findings
                   if f.code == "store.sequences_root_mismatch")
    rendered = finding.render()
    assert "lock=" in rendered and "actual=" in rendered


# --------------------------------------------------------- L1: membership

def test_unlinking_a_payload_fires_l1_and_leaves_l0_clean(
    locked, tiny_store, tiny_store_dir
) -> None:
    victim = sorted(store_census.sequence_digests(tiny_store))[0]
    store_census.seq_path(tiny_store_dir, victim).unlink()

    store = reopen(tiny_store_dir)
    # L0 is clean: the indexes still name the same digests. That is exactly the
    # blind spot L1 exists to cover -- and the one the old lock had permanently,
    # since no verify path ever opened the store.
    assert verify.verify_store_l0(locked, store).findings == []

    l1 = verify.verify_store_l1(locked, store, tiny_store_dir)
    missing = [f for f in l1.findings if f.code == "store.seq_file_missing"]
    assert len(missing) == 1
    assert missing[0].subject == victim
    assert missing[0].severity == verify.ERROR


def test_a_collection_the_store_lost_is_reported_with_its_sources(
    tiny_store, tiny_store_dir
) -> None:
    digests = sorted(store_census.collection_census(tiny_store))
    lock = lock_for_store(
        tiny_store, tiny_store_dir,
        files=[],
        from_map={digests[0]: {"host/alpha.fa"}},
    )
    lock["outputs"]["collections"].append(
        {"digest": "GONE", "n_sequences": 9, "from": []}
    )
    report = verify.verify_store_l1(lock, tiny_store, tiny_store_dir)
    finding = next(f for f in report.findings
                   if f.code == "store.collection_missing")
    assert finding.subject == "GONE"
    assert "no known contributor" in finding.render()


def test_a_collection_the_lock_does_not_describe_is_a_warning(
    tiny_store, tiny_store_dir
) -> None:
    lock = lock_for_store(tiny_store, tiny_store_dir)
    lock["outputs"]["collections"] = lock["outputs"]["collections"][:1]
    report = verify.verify_store_l1(lock, tiny_store, tiny_store_dir)
    unrecorded = [f for f in report.findings
                  if f.code == "store.collection_unrecorded"]
    assert len(unrecorded) == 1
    assert unrecorded[0].severity == verify.WARN


def test_an_unrecorded_collection_does_not_strand_its_sequences(
    tiny_store, tiny_store_dir
) -> None:
    """Orphan detection asks the *store*, not the lock.

    If it walked only locked collections, dropping one collection from the lock
    would report every sequence it published as stranded -- on the real store
    that is a single explanatory finding buried under 200,000 spurious ones.
    """
    lock = lock_for_store(tiny_store, tiny_store_dir)
    lock["outputs"]["collections"] = lock["outputs"]["collections"][:1]
    report = verify.verify_store_l1(lock, tiny_store, tiny_store_dir)
    assert "store.orphan_sequence" not in codes(report)
    assert codes(report) == ["store.collection_unrecorded"]


def test_render_caps_repeated_findings_but_counts_them_all(capsys) -> None:
    report = verify.VerifyReport()
    for i in range(50):
        report.add(verify.Finding("store.orphan_sequence", verify.WARN,
                                  f"digest{i:03d}"))
    verify.render(report, per_code=5)
    out = capsys.readouterr().out
    assert out.count("[warn ] store.orphan_sequence") == 5
    assert "… and 45 more (50 total)" in out
    assert "0 error(s), 50 warning(s)" in out


def test_a_collection_whose_size_moved_is_reported(locked, tiny_store,
                                                   tiny_store_dir) -> None:
    locked["outputs"]["collections"][0]["n_sequences"] = 99
    report = verify.verify_store_l1(locked, tiny_store, tiny_store_dir)
    assert "store.collection_size_mismatch" in codes(report)


def test_a_stray_payload_file_is_reported(locked, tiny_store,
                                          tiny_store_dir) -> None:
    stray = tiny_store_dir / "sequences" / "zz" / ("z" * 32 + ".seq")
    stray.parent.mkdir(parents=True, exist_ok=True)
    stray.write_bytes(b"\x00")
    report = verify.verify_store_l1(locked, tiny_store, tiny_store_dir)
    assert "store.stray_seq_file" in codes(report)


# --------------------------------------------------------------- L2: deep

def test_corrupting_a_payload_fires_l2_only(locked, tiny_store,
                                            tiny_store_dir) -> None:
    victim = sorted(store_census.sequence_digests(tiny_store))[0]
    path = store_census.seq_path(tiny_store_dir, victim)
    payload = bytearray(path.read_bytes())
    payload[-1] ^= 0xFF
    path.write_bytes(bytes(payload))

    store = reopen(tiny_store_dir)
    # The file is still present and still the right length, so neither the
    # indexes nor the directory listing notice. Only reading the bytes does.
    assert verify.verify_store_l0(locked, store).findings == []
    assert verify.verify_store_l1(locked, store, tiny_store_dir).findings == []

    l2 = verify.verify_store_l2(store, {})
    mismatched = [f for f in l2.findings
                  if f.code == "store.seq_digest_mismatch"]
    assert len(mismatched) == 1
    assert mismatched[0].subject == victim
    assert mismatched[0].severity == verify.ERROR


# -------------------------------------------------- the known-bad baseline

def corrupt_one(tiny_store, tiny_store_dir) -> str:
    victim = sorted(store_census.sequence_digests(tiny_store))[0]
    path = store_census.seq_path(tiny_store_dir, victim)
    payload = bytearray(path.read_bytes())
    payload[-1] ^= 0xFF
    path.write_bytes(bytes(payload))
    return victim


def baseline_for(digest: str, cause: str = "a known encoder defect") -> dict:
    return {digest: {"digest": digest, "name": "x", "alphabet": "protein",
                     "length": "1", "redigest": "?", "cause": cause}}


def test_a_baselined_defect_warns_and_renders_its_cause(
    tiny_store, tiny_store_dir
) -> None:
    victim = corrupt_one(tiny_store, tiny_store_dir)
    report = verify.verify_store_l2(
        reopen(tiny_store_dir), baseline_for(victim, "no U in the alphabet")
    )
    finding = next(f for f in report.findings
                   if f.code == "store.seq_digest_mismatch")
    assert finding.severity == verify.WARN
    assert "known cause: no U in the alphabet" in finding.render()
    assert not report.failed


def test_removing_a_digest_from_the_baseline_re_reports_it_as_an_error(
    tiny_store, tiny_store_dir
) -> None:
    """The baseline must not be able to mask anything by accident."""
    victim = corrupt_one(tiny_store, tiny_store_dir)
    store = reopen(tiny_store_dir)
    assert not verify.verify_store_l2(store, baseline_for(victim)).failed
    assert verify.verify_store_l2(store, {}).failed


def test_a_baselined_digest_that_now_round_trips_is_reported_as_fixed(
    tiny_store, tiny_store_dir
) -> None:
    """The third outcome, and the reason this is a baseline rather than a
    suppression flag: when gtars is fixed, the tool says so."""
    healthy = sorted(store_census.sequence_digests(tiny_store))[0]
    report = verify.verify_store_l2(tiny_store, baseline_for(healthy))
    finding = next(f for f in report.findings
                   if f.code == "store.seq_digest_fixed")
    assert finding.severity == verify.INFO
    assert not report.failed
    assert "gtars appears to have been fixed" in finding.render()


def test_a_baseline_naming_an_absent_digest_is_noted_not_fatal(
    tiny_store,
) -> None:
    report = verify.verify_store_l2(tiny_store, baseline_for("NOT_IN_THIS_STORE"))
    assert not report.failed
    assert any("does not hold" in note for note in report.notes)


# ------------------------------------------------------- the baseline file

def test_the_baseline_tsv_parses_with_its_comments_and_header(
    tmp_path: Path,
) -> None:
    path = tmp_path / "known.tsv"
    path.write_text(
        "# a comment\n"
        "#\n"
        "digest\tname\talphabet\tlength\tredigest\tcause\n"
        "D1\tENSP1\tprotein\t180\tR1\tno U\n",
        encoding="utf-8",
    )
    rows = verify.load_known_bad(path)
    assert list(rows) == ["D1"]
    assert rows["D1"]["cause"] == "no U"
    assert rows["D1"]["name"] == "ENSP1"


def test_an_absent_baseline_is_an_empty_baseline(tmp_path: Path) -> None:
    assert verify.load_known_bad(tmp_path / "nope.tsv") == {}


def test_the_checked_in_baseline_parses() -> None:
    if not verify.DEFAULT_KNOWN_BAD.exists():
        pytest.skip("baseline not generated in this checkout")
    rows = verify.load_known_bad(verify.DEFAULT_KNOWN_BAD)
    assert rows
    assert all(row["cause"] for row in rows.values())
    assert all(row["redigest"] != digest for digest, row in rows.items())


# ------------------------------------------------------------ limit policy

def test_limit_is_refused_for_the_store(tmp_path: Path, locked,
                                        tiny_store_dir) -> None:
    """A root is a digest over the complete set. Hashing the first N of it is
    not a weaker check, it is a different number that means nothing."""
    lock_path = tmp_path / "lock.json"
    build_lock.write_lock(lock_path, locked)
    args = type("Args", (), {
        "cache": False, "store": True, "remote": False, "manifest": False,
        "all": False, "limit": 10, "lock": lock_path,
        "store_dir": tiny_store_dir, "cache_dir": tmp_path,
        "deep": False, "no_hash": False, "known_bad": None, "jobs": 1,
        "timeout": 1, "config": tmp_path / "sources.toml",
    })()
    with pytest.raises(SystemExit, match="partial root is meaningless"):
        verify.run_verify(args)
