"""Translate sample expressions on the refget backend through the eval harness's
WorkerRuntime and compare each VRS ID with the one SeqRepo produced.

Usage: check_new_sources.py HARNESS_DIR REFGET_URI UTA_URL SAMPLES_JSON
"""

import json
import sys

harness_dir, uri, uta_url, samples_path = sys.argv[1:5]
sys.path.insert(0, harness_dir)
import harness  # noqa: E402

samples = json.load(open(samples_path))
rt = harness.WorkerRuntime("refget", uri, uta_url, True)
rt.initialize()
same = 0
for sample in samples:
    row = rt.translate({"expression": sample["expression"], "format": sample["format"]})
    ok = row["status"] == "success" and row["vrs_id"] == sample["seqrepo_vrs_id"]
    same += ok
    exc = (row["exception"] or {}).get("type")
    print(f"{'SAME' if ok else 'DIFF'}  {sample['expression']:50}  {row['status']:8} {row['vrs_id'] or exc}")
print(f"{same}/{len(samples)} match SeqRepo")
