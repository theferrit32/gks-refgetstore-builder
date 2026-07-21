#!/usr/bin/env python3
"""Discover NCBI assembly/annotation files and emit reviewable TOML fragments.

This is intentionally a manifest *candidate* generator: it writes only stdout
(and an optional JSONL catalog), and never changes sources.toml.
"""

from __future__ import annotations

import argparse
import csv
import ftplib
import json
import re
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable
from urllib.parse import urlsplit

HOST = "ftp.ncbi.nlm.nih.gov"
HTTPS = f"https://{HOST}"
FTP_ROOT = "/genomes/all/annotation_releases"
SUMMARY_URLS = (
    f"{HTTPS}/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt",
    f"{HTTPS}/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt",
)
DEFAULT_PREFIX = "GCF_000001405."
DEFAULT_EXACT = "GCF_009914755.1"
ACCESSION_RE = re.compile(r"GCF_\d{9}\.\d+")
REPORT_ONLY = {
    "GCF_000001405.23": "GRCh37.p11",
    "GCF_000001405.24": "GRCh37.p12",
}
FILE_SUFFIXES = {
    "genomic_fasta": "_genomic.fna.gz",
    "assembly_report": "_assembly_report.txt",
    "rna_fasta": "_rna.fna.gz",
    "protein_fasta": "_protein.faa.gz",
}


@dataclass(frozen=True)
class Assembly:
    accession: str
    name: str
    taxid: int
    status: str
    release_date: str
    version_status: str
    ftp_path: str | None


@dataclass(frozen=True)
class Candidate:
    source_kind: str
    accession: str
    assembly_name: str
    taxid: int
    annotation_run_directory: str | None
    file_type: str
    url: str | None
    readme_url: str | None
    discovery_status: str

    def as_dict(self) -> dict[str, object]:
        return self.__dict__.copy()


def mirror_cache_path(cache_dir: Path, url: str) -> Path:
    """Use the same host/full-path cache convention as build_store.py."""
    parsed = urlsplit(url)
    return cache_dir / parsed.netloc / parsed.path.lstrip("/")


def parse_summary(stream: Iterable[str]) -> list[Assembly]:
    """Stream an NCBI assembly_summary file, tolerating added columns."""
    header: list[str] | None = None
    rows: list[Assembly] = []
    for line in stream:
        # NCBI currently writes ``#assembly_accession``; older fixtures and
        # mirrors may include a space after the comment marker.
        possible_header = line.lstrip("# ").rstrip("\r\n")
        if possible_header.startswith("assembly_accession\t"):
            header = possible_header.split("\t")
            break
    if header is None:
        raise ValueError("assembly summary has no '#assembly_accession' header")
    for row in csv.DictReader(stream, fieldnames=header, delimiter="\t"):
        accession = row.get("assembly_accession", "")
        if not accession:
            continue
        raw_ftp = row.get("ftp_path", "na")
        rows.append(Assembly(
            accession=accession,
            name=row.get("asm_name", accession),
            taxid=int(row.get("taxid") or 0),
            status=row.get("refseq_category", ""),
            release_date=row.get("seq_rel_date", ""),
            version_status=row.get("version_status", ""),
            ftp_path=None if raw_ftp.lower() == "na" else raw_ftp,
        ))
    return rows


def _download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "gks-refgetstore-builder"})
    with urllib.request.urlopen(request, timeout=60) as response, tmp.open("wb") as out:  # noqa: S310
        while chunk := response.read(1 << 20):
            out.write(chunk)
    tmp.replace(target)


def obtain_summary(url: str, cache_dir: Path, refresh: bool) -> Path:
    target = mirror_cache_path(cache_dir, url)
    cwd_fallback = Path.cwd() / Path(urlsplit(url).path).name
    if target.exists() and not refresh:
        print(f"summary cache hit: {target}", file=sys.stderr)
        return target
    if cwd_fallback.exists() and not refresh:
        print(f"summary cwd fallback: {cwd_fallback}", file=sys.stderr)
        return cwd_fallback
    print(f"downloading summary: {url}", file=sys.stderr)
    _download(url, target)
    return target


class FTPListings:
    """Cached FTP NLST access; cache files contain one basename per line."""

    def __init__(self, cache_dir: Path, refresh: bool = False) -> None:
        self.cache_dir, self.refresh = cache_dir, refresh
        self.ftp: ftplib.FTP | None = None

    def _cache_path(self, directory: str) -> Path:
        url = f"ftp://{HOST}{directory.rstrip('/')}/.listing"
        return mirror_cache_path(self.cache_dir, url)

    def list(self, directory: str) -> list[str]:
        cache = self._cache_path(directory)
        if cache.exists() and not self.refresh:
            print(f"listing cache hit: {directory}", file=sys.stderr)
            return [line for line in cache.read_text().splitlines() if line]
        if self.ftp is None:
            print(f"connecting to ftp://{HOST}", file=sys.stderr)
            self.ftp = ftplib.FTP(HOST, timeout=60)
            self.ftp.login()
        print(f"listing FTP: {directory}", file=sys.stderr)
        names = sorted(Path(item.rstrip("/")).name for item in self.ftp.nlst(directory))
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text("".join(f"{name}\n" for name in names))
        return names

    def close(self) -> None:
        if self.ftp is not None:
            try:
                self.ftp.quit()
            except OSError:
                self.ftp.close()


def annotation_run_index(taxid: int, listing: Callable[[str], list[str]]) -> dict[str, list[str]]:
    """Map accessions to annotation-run dirs in known NCBI one/two-level layout."""
    root = f"{FTP_ROOT}/{taxid}"
    runs: dict[str, list[str]] = {}
    for entry in listing(root):
        first = f"{root}/{entry}"
        if ACCESSION_RE.search(entry):
            children = [entry]
            parent = root
        elif re.fullmatch(r"\d+(?:\.\d+)?", entry):
            try:
                children = listing(first)
            except ftplib.error_perm as exc:
                print(f"warning: cannot list {first}: {exc}", file=sys.stderr)
                continue
            parent = first
        else:
            continue
        for child in children:
            match = ACCESSION_RE.search(child)
            if match:
                runs.setdefault(match.group(), []).append(f"{parent}/{child}")
    return {key: sorted(set(value), key=run_sort_key) for key, value in runs.items()}


def run_sort_key(path: str) -> tuple[object, ...]:
    release = path.split("/")[-2] if ACCESSION_RE.search(path.split("/")[-1]) else ""
    update = path.split("/")[-1]
    pieces = re.split(r"(\d+)", release + "/" + update)
    return tuple(int(piece) if piece.isdigit() else piece for piece in pieces)


def assembly_sort_key(assembly: Assembly) -> tuple[object, ...]:
    name = assembly.name
    family = 0 if name.startswith("GRCh37") else 1 if name.startswith("GRCh38") else 2
    patch = re.search(r"\.p(\d+)$", name)
    return family, int(patch.group(1)) if patch else 0, assembly.release_date, assembly.accession


def inspect_directory(assembly: Assembly, directory: str, listing: Callable[[str], list[str]],
                      source_kind: str) -> list[Candidate]:
    try:
        names = set(listing(directory))
    except ftplib.error_perm as exc:
        print(f"warning: cannot list {directory}: {exc}", file=sys.stderr)
        return [Candidate(source_kind, assembly.accession, assembly.name, assembly.taxid,
                          directory if source_kind == "annotation_run" else None,
                          kind, None, None, "directory_unavailable") for kind in FILE_SUFFIXES]
    # Dated/release-number directories are named after the assembly and direct
    # update directories are named ``GCF_…-RS_…``. In both cases the payload
    # files retain the canonical ``<accession>_<assembly-name>`` prefix.
    prefix = f"{assembly.accession}_{assembly.name}"
    readme_name = next((n for n in names if n.lower().startswith("readme")), None)
    readme = f"{HTTPS}{directory}/{readme_name}" if readme_name else f"{HTTPS}{directory}/README.txt"
    found: list[Candidate] = []
    for kind, suffix in FILE_SUFFIXES.items():
        filename = prefix + suffix
        url = f"{HTTPS}{directory}/{filename}" if filename in names else None
        status = "discovered" if url else "missing"
        if url is None:
            print(f"missing {kind}: {HTTPS}{directory}/{filename}", file=sys.stderr)
        found.append(Candidate(source_kind, assembly.accession, assembly.name, assembly.taxid,
                               directory if source_kind == "annotation_run" else None,
                               kind, url, readme, status))
    return found


def deduplicate_candidates(candidates: list[Candidate]) -> list[Candidate]:
    """Keep one discoverable candidate per URL and catalog later copies as skipped."""
    seen: set[str] = set()
    result: list[Candidate] = []
    for candidate in candidates:
        if (candidate.discovery_status == "discovered" and candidate.url
                and candidate.url in seen):
            candidate = Candidate(**{
                **candidate.as_dict(), "discovery_status": "skipped_duplicate_url"
            })
        elif candidate.discovery_status == "discovered" and candidate.url:
            seen.add(candidate.url)
        result.append(candidate)
    return result


def select_assemblies(rows: Iterable[Assembly], taxid: int, all_human: bool,
                      includes: list[str], excludes: list[str]) -> list[Assembly]:
    by_accession = {row.accession: row for row in rows if row.taxid == taxid}
    selected = set(by_accession) if all_human else {
        accession for accession in by_accession
        if accession.startswith(DEFAULT_PREFIX) or accession == DEFAULT_EXACT
    }
    selected.update(includes)
    selected.difference_update(excludes)
    missing = sorted(selected - by_accession.keys())
    for accession in missing:
        print(f"warning: requested accession absent from summaries: {accession}", file=sys.stderr)
    return sorted((by_accession[a] for a in selected if a in by_accession), key=assembly_sort_key)


def _q(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def emit_toml(assemblies: list[Assembly], candidates: list[Candidate], section: str) -> str:
    lines = ["# Generated NCBI source candidates; review before copying to sources.toml."]
    discovered = [candidate for candidate in candidates if candidate.discovery_status == "discovered"]
    if section in ("all", "assemblies"):
        by_accession: dict[str, dict[str, Candidate]] = {}
        for candidate in discovered:
            if candidate.source_kind == "assembly":
                by_accession.setdefault(candidate.accession, {})[candidate.file_type] = candidate
        for assembly in assemblies:
            files = by_accession.get(assembly.accession, {})
            report = files.get("assembly_report")
            fasta = files.get("genomic_fasta")
            if report is None:
                continue
            lines.extend(["", "[[assembly]]", f"namespace = {_q(assembly.name)}"])
            if fasta is not None and assembly.accession not in REPORT_ONLY:
                lines.append(f"fasta_url = {_q(fasta.url or '')}")
                load = True
            else:
                load = False
            lines.extend([f"report_url = {_q(report.url or '')}", f"load_fasta = {str(load).lower()}"])
    if section in ("all", "history"):
        for file_type, name in (("rna_fasta", "refseq_history_rna"),
                                ("protein_fasta", "refseq_history_protein")):
            urls = sorted({c.url for c in discovered if c.file_type == file_type and c.url},
                          key=run_sort_key)
            lines.extend(["", "[[seqset]]", f"name = {_q(name)}",
                          'namespace = "refseq"', "urls = ["])
            lines.extend(f"  {_q(url)}," for url in urls)
            lines.append("]")
    return "\n".join(lines) + "\n"


def warn_unpaired(candidates: list[Candidate]) -> None:
    grouped: dict[tuple[str, str], set[str]] = {}
    for candidate in candidates:
        if candidate.discovery_status == "discovered" and candidate.file_type in ("rna_fasta", "protein_fasta"):
            key = (candidate.accession, candidate.annotation_run_directory or candidate.url.rsplit("/", 1)[0])
            grouped.setdefault(key, set()).add(candidate.file_type)
    for (accession, directory), kinds in sorted(grouped.items()):
        if len(kinds) == 1:
            missing = "protein" if "rna_fasta" in kinds else "RNA"
            print(f"warning: {accession} {directory} has no {missing} counterpart", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxid", type=int, default=9606)
    parser.add_argument("--cache-dir", type=Path, default=Path("downloads"))
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--all-human", action="store_true",
                        help="select every RefSeq assembly for --taxid (name retained for compatibility)")
    parser.add_argument("--include-accession", action="append", default=[])
    parser.add_argument("--exclude-accession", action="append", default=[])
    parser.add_argument("--catalog-out", type=Path)
    parser.add_argument("--section", choices=("all", "assemblies", "history"), default="all")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows: list[Assembly] = []
    for url in SUMMARY_URLS:
        path = obtain_summary(url, args.cache_dir, args.refresh)
        with path.open(encoding="utf-8") as stream:
            rows.extend(parse_summary(stream))
    # Current wins over historical when an accession appears in both.
    rows = list({row.accession: row for row in reversed(rows)}.values())
    assemblies = select_assemblies(rows, args.taxid, args.all_human,
                                   args.include_accession, args.exclude_accession)
    print(f"selected {len(assemblies)} assembly row(s)", file=sys.stderr)
    ftp = FTPListings(args.cache_dir, args.refresh)
    candidates: list[Candidate] = []
    try:
        run_index = annotation_run_index(args.taxid, ftp.list)
        for assembly in assemblies:
            directory = None
            if assembly.ftp_path:
                directory = urlsplit(assembly.ftp_path).path
            elif assembly.accession in REPORT_ONLY:
                base = "/genomes/all/GCF/000/001/405"
                directory = f"{base}/{assembly.accession}_{REPORT_ONLY[assembly.accession]}"
            if directory:
                candidates.extend(inspect_directory(assembly, directory, ftp.list, "assembly"))
            else:
                candidates.extend(Candidate("assembly", assembly.accession, assembly.name,
                                            assembly.taxid, None, kind, None, None,
                                            "no_ftp_path") for kind in FILE_SUFFIXES)
            for run in run_index.get(assembly.accession, []):
                candidates.extend(inspect_directory(assembly, run, ftp.list, "annotation_run"))
    finally:
        ftp.close()
    candidates = deduplicate_candidates(candidates)
    warn_unpaired(candidates)
    if args.catalog_out:
        args.catalog_out.parent.mkdir(parents=True, exist_ok=True)
        with args.catalog_out.open("w", encoding="utf-8") as out:
            for candidate in candidates:
                out.write(json.dumps(candidate.as_dict(), sort_keys=True) + "\n")
        print(f"wrote {len(candidates)} catalog records to {args.catalog_out}", file=sys.stderr)
    sys.stdout.write(emit_toml(assemblies, candidates, args.section))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
