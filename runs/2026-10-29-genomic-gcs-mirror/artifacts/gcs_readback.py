"""Translate the R2 upload's sample expressions through vrsify's worker against a remote store.

Usage: gcs_readback.py STORE_URL CACHE_ROOT SAMPLES_JSON
"""
import json
import sys
from pathlib import Path

from clinvar_ingest.vrsify import runner, store
from clinvar_ingest.vrsify.translate import translate_record

url, cache_root, samples_path = sys.argv[1:]
location, provenance = store.prepare_store(url, Path(cache_root))
print("provenance:", json.dumps(provenance))
runner._init_worker(runner.SequenceSource("refget", location, cache_root=cache_root), 120)
w = runner._WORKER
ok = 0
for s in json.load(open(samples_path)):
    out = translate_record({"vrs_class": "Allele", "fmt": s["format"], "source": s["expression"]},
                           w["allele_translators"], w["cnv_translator"])
    got = out.get("id") or out.get("errors")
    match = got == s["seqrepo_vrs_id"]
    ok += match
    print(f"{'ok  ' if match else 'DIFF'} {s['expression']:<40} {got}")
print(f"{ok}/16 match")
