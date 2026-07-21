from __future__ import annotations

import io
import tempfile
import tomllib
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch

import generate_ncbi_source_candidates as gen
from build_store import load_config


HEADER = ("#assembly_accession\trefseq_category\ttaxid\tversion_status\t"
          "seq_rel_date\tasm_name\tftp_path\n")


def summary(*rows: str) -> list[gen.Assembly]:
    return gen.parse_summary(io.StringIO(HEADER + "".join(rows)))


class SummaryTests(unittest.TestCase):
    def test_parse_summary_and_na(self) -> None:
        rows = summary(
            "GCF_000001405.23\treference genome\t9606\tsuppressed\t2013-06-28\t"
            "GRCh37.p11\tna\n",
            "GCF_009914755.1\treference genome\t9606\tlatest\t2022-01-01\t"
            "T2T-CHM13v2.0\tftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/x\n",
        )
        self.assertIsNone(rows[0].ftp_path)
        self.assertEqual(rows[0].status, "reference genome")
        self.assertEqual(rows[1].name, "T2T-CHM13v2.0")

    def test_header_with_space_is_also_accepted(self) -> None:
        stream = io.StringIO(HEADER.replace("#assembly", "# assembly") +
                             "GCF_1\tna\t9606\tlatest\t2022\tExample\tna\n")
        self.assertEqual(gen.parse_summary(stream)[0].accession, "GCF_1")

    def test_default_and_override_selection(self) -> None:
        rows = summary(
            "GCF_000001405.40\treference genome\t9606\tlatest\t2022\tGRCh38.p14\tna\n",
            "GCF_009914755.1\treference genome\t9606\tlatest\t2022\tT2T\tna\n",
            "GCF_111111111.1\tna\t9606\tlatest\t2022\tOther\tna\n",
        )
        selected = gen.select_assemblies(rows, 9606, False, ["GCF_111111111.1"],
                                         ["GCF_009914755.1"])
        self.assertEqual({r.accession for r in selected},
                         {"GCF_000001405.40", "GCF_111111111.1"})


class ListingTests(unittest.TestCase):
    def test_annotation_layouts_and_reverse_index(self) -> None:
        root = f"{gen.FTP_ROOT}/9606"
        listings = {
            root: ["109", "109.20211119", "GCF_000001405.40-RS_2024_08"],
            f"{root}/109": ["README", "GCF_000001405.38_GRCh38.p12"],
            f"{root}/109.20211119": ["GCF_000001405.39_GRCh38.p13"],
        }
        index = gen.annotation_run_index(9606, listings.__getitem__)
        self.assertEqual(index["GCF_000001405.38"],
                         [f"{root}/109/GCF_000001405.38_GRCh38.p12"])
        self.assertIn(f"{root}/GCF_000001405.40-RS_2024_08",
                      index["GCF_000001405.40"])

    def test_listing_cache_first_and_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cached = gen.FTPListings(Path(tmp))
            path = cached._cache_path("/somewhere")
            path.parent.mkdir(parents=True)
            path.write_text("b\na\n")
            self.assertEqual(cached.list("/somewhere"), ["b", "a"])
            refreshed = gen.FTPListings(Path(tmp), refresh=True)

            class FakeFTP:
                def __init__(self, *args, **kwargs): pass
                def login(self): pass
                def nlst(self, directory): return [directory + "/new"]
                def quit(self): pass

            with patch.object(gen.ftplib, "FTP", FakeFTP):
                self.assertEqual(refreshed.list("/somewhere"), ["new"])

    def test_cwd_summary_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cwd = Path(tmp)
            fallback = cwd / "assembly_summary_refseq.txt"
            fallback.write_text(HEADER)
            with patch.object(gen.Path, "cwd", return_value=cwd):
                result = gen.obtain_summary(gen.SUMMARY_URLS[0], cwd / "cache", False)
            self.assertEqual(result.name, "assembly_summary_refseq.txt")


class CandidateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assembly = gen.Assembly("GCF_000001405.23", "GRCh37.p11", 9606,
                                     "reference", "2013", "suppressed", None)

    def test_report_only_and_valid_config(self) -> None:
        directory = "/genomes/all/GCF/000/001/405/GCF_000001405.23_GRCh37.p11"
        names = ["GCF_000001405.23_GRCh37.p11_assembly_report.txt"]
        candidates = gen.inspect_directory(self.assembly, directory, lambda _: names, "assembly")
        text = gen.emit_toml([self.assembly], candidates, "assemblies")
        parsed = tomllib.loads(text)
        self.assertFalse(parsed["assembly"][0]["load_fasta"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidate.toml"
            path.write_text(text)
            assemblies, seqsets = load_config(path)
        self.assertEqual(assemblies[0].namespace, "GRCh37.p11")
        self.assertEqual(seqsets, [])

    def test_direct_rs_update_uses_assembly_filename_prefix(self) -> None:
        assembly = gen.Assembly("GCF_000001405.40", "GRCh38.p14", 9606,
                                "reference", "2022", "latest", None)
        directory = "/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2025_08"
        names = [
            "GCF_000001405.40_GRCh38.p14_rna.fna.gz",
            "GCF_000001405.40_GRCh38.p14_protein.faa.gz",
        ]
        candidates = gen.inspect_directory(assembly, directory, lambda _: names,
                                           "annotation_run")
        found = {candidate.file_type: candidate for candidate in candidates
                 if candidate.discovery_status == "discovered"}
        self.assertTrue(found["rna_fasta"].url.endswith(names[0]))
        self.assertTrue(found["protein_fasta"].url.endswith(names[1]))

    def test_url_dedup_history_and_valid_config(self) -> None:
        url = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/109/"
               "GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_rna.fna.gz")
        candidate = gen.Candidate("annotation_run", self.assembly.accession,
                                  self.assembly.name, 9606, "/run", "rna_fasta",
                                  url, None, "discovered")
        deduplicated = gen.deduplicate_candidates([candidate, candidate])
        self.assertEqual(deduplicated[1].discovery_status, "skipped_duplicate_url")
        text = gen.emit_toml([], deduplicated, "history")
        self.assertEqual(text.count(url), 1)
        parsed = tomllib.loads(text)
        self.assertEqual(parsed["seqset"][0]["name"], "refseq_history_grch37_rna")
        self.assertEqual(len(parsed["seqset"]), 1)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidate.toml"
            path.write_text(text)
            _, seqsets = load_config(path)
        self.assertEqual(len(seqsets), 1)

    def test_trailing_slashes_are_normalized_everywhere(self) -> None:
        directory = "/genomes/all/GCF/000/001/405/GCF_000001405.23_GRCh37.p11//"
        prefix = "GCF_000001405.23_GRCh37.p11"
        names = [prefix + "_rna.fna.gz", "README.txt"]
        listed = []

        def listing(path):
            listed.append(path)
            return names

        stderr = io.StringIO()
        with redirect_stderr(stderr):
            candidates = gen.inspect_directory(self.assembly, directory, listing, "annotation_run")
            gen.warn_unpaired(candidates)
        self.assertEqual(listed, [directory.rstrip("/").replace("//", "/")])
        for candidate in candidates:
            self.assertNotIn("//", (candidate.annotation_run_directory or "").lstrip("/"))
            self.assertNotIn("//", (candidate.url or "").split("://")[-1])
            self.assertNotIn("//", (candidate.readme_url or "").split("://")[-1])
        self.assertNotIn("//GCF", stderr.getvalue())

        dirty = gen.Candidate("annotation_run", self.assembly.accession,
                              self.assembly.name, 9606, "/run//", "rna_fasta",
                              "https://example.test/run//x.fna.gz",
                              "https://example.test/run//README", "discovered")
        clean = gen.deduplicate_candidates([dirty])[0]
        self.assertEqual(clean.url, "https://example.test/run/x.fna.gz")
        self.assertEqual(clean.annotation_run_directory, "/run")

    def test_history_is_grouped_in_family_pairs_and_ordered(self) -> None:
        def pair(accession, name, stem):
            return [
                gen.Candidate("assembly", accession, name, 9606, None, "rna_fasta",
                              f"https://example.test/{stem}_rna.fna.gz", None, "discovered"),
                gen.Candidate("assembly", accession, name, 9606, None, "protein_fasta",
                              f"https://example.test/{stem}_protein.faa.gz", None, "discovered"),
            ]

        candidates = (pair("GCF_000001405.40", "GRCh38.p14", "38")
                      + pair("GCF_009914755.1", "T2T-CHM13v2.0", "t2t")
                      + pair("GCF_000001405.25", "GRCh37.p13", "37"))
        parsed = tomllib.loads(gen.emit_toml([], candidates, "history"))
        self.assertEqual([entry["name"] for entry in parsed["seqset"]], [
            "refseq_history_grch37_rna", "refseq_history_grch37_protein",
            "refseq_history_grch38_rna", "refseq_history_grch38_protein",
            "refseq_history_t2t_chm13_rna", "refseq_history_t2t_chm13_protein",
        ])

    def test_arbitrary_family_labels_are_stable_and_collision_safe(self) -> None:
        candidates = [
            gen.Candidate("assembly", "GCF_222222222.1", "Foo-Bar.p2", 9606,
                          None, "rna_fasta", "https://x/b", None, "discovered"),
            gen.Candidate("assembly", "GCF_111111111.1", "Foo Bar.p1", 9606,
                          None, "rna_fasta", "https://x/a", None, "discovered"),
        ]
        first = tomllib.loads(gen.emit_toml([], candidates, "history"))
        second = tomllib.loads(gen.emit_toml([], list(reversed(candidates)), "history"))
        names = [entry["name"] for entry in first["seqset"]]
        self.assertEqual(names, [entry["name"] for entry in second["seqset"]])
        self.assertEqual(names, [
            "refseq_history_foo_bar_gcf_111111111_1_rna",
            "refseq_history_foo_bar_gcf_222222222_1_rna",
        ])

    def test_unpaired_warning(self) -> None:
        candidate = gen.Candidate("assembly", self.assembly.accession,
                                  self.assembly.name, 9606, None, "rna_fasta",
                                  "https://example/x_rna.fna.gz", None, "discovered")
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            gen.warn_unpaired([candidate])
        self.assertIn("no protein counterpart", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
