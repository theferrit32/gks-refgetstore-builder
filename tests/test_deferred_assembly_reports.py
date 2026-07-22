from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

import build_store


REPORT_COLUMNS = (
    "Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
    "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\t"
    "UCSC-style-name"
)


def write_report(path: Path, accession: str, rows: list[tuple[str, str, str, str]]) -> None:
    lines = [
        f"# RefSeq assembly accession: {accession}",
        f"# GenBank assembly accession: GCA_{accession.removeprefix('GCF_')}",
        f"# {REPORT_COLUMNS}",
    ]
    for sequence_name, genbank, refseq, ucsc in rows:
        lines.append(
            f"{sequence_name}\tassembled-molecule\t{sequence_name}\tChromosome\t"
            f"{genbank}\t=\t{refseq}\tPrimary Assembly\t1\t{ucsc}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def cached_path(cache: Path, url: str) -> Path:
    path = build_store.mirror_cache_path(cache, url)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def build_args(tmp_path: Path, config: Path, **overrides):
    values = {
        "config": config,
        "cache_dir": tmp_path / "cache",
        "store_dir": tmp_path / "store",
        "lock": tmp_path / "build.lock.json",
        "locked_sources": False,
        "skip_assemblies": False,
        "skip_seqsets": False,
        "assembly": None,
        "seqset": None,
        "lock_check_mode": "ignore",
        "force_download": False,
        "force_lock": False,
        "no_lock": True,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_report_only_uses_later_fasta_but_normal_report_stays_collection_scoped(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    base = "https://example.test/"
    early_report = base + "early_report.txt"
    scoped_fasta = base + "scoped.fa"
    scoped_report = base + "scoped_report.txt"
    later_fasta = base + "later.fa"
    later_report = base + "later_report.txt"
    config = tmp_path / "sources.toml"
    config.write_text(
        "\n".join([
            "[[assembly]]", 'namespace="early"', f'report_url="{early_report}"',
            "load_fasta=false", "",
            "[[assembly]]", 'namespace="scoped"', f'report_url="{scoped_report}"',
            f'fasta_url="{scoped_fasta}"', "",
            "[[assembly]]", 'namespace="later"', f'report_url="{later_report}"',
            f'fasta_url="{later_fasta}"', "",
        ])
    )
    cache = tmp_path / "cache"
    cached_path(cache, scoped_fasta).write_text(">OWN.1\nA\n")
    cached_path(cache, later_fasta).write_text(">LATE.1\nC\n>BORROW.1\nG\n")
    write_report(cached_path(cache, early_report), "GCF_000000001.1", [
        ("early-name", "GB_LATE.1", "LATE.1", "chrEarly"),
        ("missing-name", "GB_MISSING.1", "MISSING.1", "chrMissing"),
    ])
    write_report(cached_path(cache, scoped_report), "GCF_000000002.1", [
        ("borrowed-name", "GB_BORROW.1", "BORROW.1", "chrBorrowed"),
    ])
    write_report(cached_path(cache, later_report), "GCF_000000003.1", [
        ("later-name", "GB_LATE.1", "LATE.1", "chrLater"),
    ])

    caplog.set_level("WARNING")
    assert build_store.run_build(build_args(tmp_path, config)) == 0

    store = build_store.RefgetStore.on_disk(str(tmp_path / "store"))
    assert store.get_sequence_by_alias("early", "early-name") is not None
    assert store.get_sequence_by_alias("later", "later-name") is not None
    assert store.get_sequence_by_alias("early", "missing-name") is None
    assert store.get_sequence_by_alias("scoped", "borrowed-name") is None
    assert "completed store does not contain it" in caplog.text
    assert "associated FASTA collection does not contain it" in caplog.text


class FakeStore:
    storage_mode = "test"

    @classmethod
    def on_disk(cls, _path):
        return cls()

    def write(self):
        pass

    def stats(self):
        return {}

    def list_sequence_alias_namespaces(self):
        return []

    def list_collection_alias_namespaces(self):
        return []


@pytest.mark.parametrize(
    ("overrides", "expected"),
    [
        ({}, ["ingest:asm", "seq:seq", "parse:asm", "lookup", "apply:asm"]),
        ({"assembly": "asm"}, ["ingest:asm", "parse:asm", "lookup", "apply:asm"]),
        ({"seqset": "seq"}, ["seq:seq", "lookup"]),
        ({"skip_assemblies": True}, ["seq:seq", "lookup"]),
        ({"skip_seqsets": True}, ["ingest:asm", "parse:asm", "lookup", "apply:asm"]),
    ],
)
def test_build_phase_order_and_partial_modes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, overrides: dict, expected: list[str]
) -> None:
    config = tmp_path / "sources.toml"
    config.write_text(
        "[[assembly]]\nnamespace=\"asm\"\n"
        "report_url=\"https://example.test/report\"\n"
        "fasta_url=\"https://example.test/fasta\"\n\n"
        "[[seqset]]\nname=\"seq\"\nnamespace=\"refseq\"\n"
        "url_template=\"https://example.test/seq\"\n"
    )
    events: list[str] = []

    def ingest(_store, entry, _cache, _provenance):
        events.append(f"ingest:{entry.namespace}")
        return build_store.DeferredAssemblyReport(
            entry, Path("report"), build_store.AssemblyStats(entry.namespace), "coll"
        )

    def seq(_store, entry, _cache, _provenance, refresh_derived=False):
        events.append(f"seq:{entry.name}")
        return build_store.SeqsetStats(entry.name, entry.namespace)

    def parse(context):
        events.append(f"parse:{context.entry.namespace}")
        context.report = build_store.AssemblyReport(None, None, [])

    def lookup(_store, required):
        assert required == set()
        events.append("lookup")
        return {}

    def apply(_store, context, _lookup):
        events.append(f"apply:{context.entry.namespace}")
        return context.stats

    monkeypatch.setattr(build_store, "ensure_download", lambda *args, **kwargs: args[1])
    monkeypatch.setattr(build_store, "RefgetStore", FakeStore)
    monkeypatch.setattr(build_store, "ingest_assembly", ingest)
    monkeypatch.setattr(build_store, "process_seqset", seq)
    monkeypatch.setattr(build_store, "parse_deferred_assembly_report", parse)
    monkeypatch.setattr(build_store, "build_targeted_name_to_digest_map", lookup)
    monkeypatch.setattr(build_store, "apply_assembly_report", apply)

    assert build_store.run_build(build_args(tmp_path, config, **overrides)) == 0
    assert events == expected
    if "apply:asm" in events:
        last_ingest = max(
            events.index("ingest:asm"),
            events.index("seq:seq") if "seq:seq" in events else 0,
        )
        assert last_ingest < events.index("apply:asm")
