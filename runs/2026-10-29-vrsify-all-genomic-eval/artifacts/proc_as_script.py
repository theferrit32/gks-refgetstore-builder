"""Print a procedure file as a one-off BigQuery script: its parameters become DECLAREs.

Usage: proc_as_script.py PROC_FILE NAME=SQL_VALUE...
"""
import re
import sys

text = open(sys.argv[1]).read()
values = dict(arg.split("=", 1) for arg in sys.argv[2:])
header = re.search(r"^CREATE OR REPLACE PROCEDURE `[^`]+`\((.*)\)\n", text, re.M)
params = [p.strip().split() for p in header.group(1).split(",")]
missing = {name for name, _ in params} - values.keys()
if missing:
    sys.exit(f"missing values for {sorted(missing)}")
declares = "".join(f"DECLARE {name} {kind} DEFAULT {values[name]};\n" for name, kind in params)
print(text[: header.start()] + declares + text[header.end():])
