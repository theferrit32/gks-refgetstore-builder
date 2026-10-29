"""Resolve NG_ accessions through vrsify's worker data proxy against a remote store.

Usage: remote_ng_check.py STORE_URL CACHE_ROOT ACCESSION...
"""
import sys
from pathlib import Path

from clinvar_ingest.vrsify import runner, store

url, cache_root, *accessions = sys.argv[1:]
location, _ = store.prepare_store(url, Path(cache_root))
runner._init_worker(runner.SequenceSource("refget", location, cache_root=cache_root), 120)
proxy = runner._WORKER["allele_translators"]["38"].data_proxy
for acc in accessions:
    md = proxy.get_metadata(f"refseq:{acc}")
    ga4gh = [a for a in md["aliases"] if a.startswith("ga4gh:")] if md else []
    print(acc, "length", md and md["length"], ga4gh)
