"""Verify and analyze the all-genomic vrsify outputs.

Usage: allgen_analysis.py WORK_DIR OUTPUT_DIR

OUTPUT_DIR holds <round>-<arm>.jsonl.gz. Prints:
1. per output: record count, distinct keys, and an order-independent fingerprint of
   (key, out), so outputs can be compared for identical content without a join
2. for any output whose fingerprint differs from r1-A, a categorized comparison with it
3. within r1-A: whether a variation's expressions on the same sequence give one VRS id
4. within r1-A: the input `issue` against the output `errors`
A record's key is (variation_id, fmt, assembly_version, source).
"""

import sys
from pathlib import Path

import duckdb

work, outputs = Path(sys.argv[1]), Path(sys.argv[2])
con = duckdb.connect(str(work / "allgen.duckdb"))
con.sql("SET memory_limit = '24GB'")
con.sql(f"SET temp_directory = '{work}/duckdb-tmp'")
con.sql("SET preserve_insertion_order = false")

KEY = "concat_ws('|', json ->> '$.in.variation_id', coalesce(json ->> '$.in.fmt', ''), coalesce(json ->> '$.in.assembly_version', ''), coalesce(json ->> '$.in.source', ''))"


def source(path: Path) -> str:
    return f"read_json_objects('{path}', format = 'newline_delimited')"


def show(sql: str, limit: int = 40) -> None:
    rel = con.sql(sql)
    print("  " + " | ".join(rel.columns))
    for row in rel.fetchmany(limit):
        print("  " + " | ".join(str(v) for v in row))


files = sorted(outputs.glob("*.jsonl.gz"))
print("== 1. fingerprints")
prints = {}
for path in files:
    n, keys, h1, h2 = con.sql(f"""
        SELECT count(*), count(DISTINCT k), bit_xor(hash(k || chr(9) || o)), bit_xor(hash(o || chr(9) || k))
        FROM (SELECT {KEY} AS k, CAST(json -> '$.out' AS VARCHAR) AS o FROM {source(path)})
    """).fetchone()
    prints[path.stem.removesuffix(".jsonl")] = (n, h1, h2)
    print(f"  {path.name:20} records {n:>10,} distinct keys {keys:>10,} fingerprint {h1:x}/{h2:x}")
reference = prints["r1-A"]
same = [name for name, fp in prints.items() if fp == reference]
differ = [name for name, fp in prints.items() if fp != reference]
print(f"  identical to r1-A: {len(same)} of {len(prints)}; differ: {differ or 'none'}")

print("\n== 2. differences from r1-A")
con.sql(f"""
    CREATE OR REPLACE TABLE a AS
    SELECT {KEY} AS k,
           json ->> '$.in.variation_id' AS variation_id,
           json ->> '$.in.vrs_class' AS vrs_class,
           json ->> '$.in.fmt' AS fmt,
           json ->> '$.in.assembly_version' AS assembly_version,
           json ->> '$.in.source' AS source,
           CAST(json ->> '$.in.precedence' AS INTEGER) AS precedence,
           json ->> '$.in.hgvs_type' AS hgvs_type,
           json ->> '$.in.issue' AS issue,
           CAST(json ->> '$.in.selected' AS BOOLEAN) AS selected,
           json ->> '$.out.errors' AS errors,
           json ->> '$.out.id' AS vrs_id,
           json ->> '$.out.type' AS vrs_type,
           json ->> '$.out.location.sequenceReference.refgetAccession' AS sequence,
           CAST(json -> '$.out' AS VARCHAR) AS out
    FROM {source(outputs / 'r1-A.jsonl.gz')}
""")
if not differ:
    print("  none: every output holds the same (key, out) pairs as r1-A")
for name in differ:
    con.sql(f"CREATE OR REPLACE TEMP TABLE o AS SELECT {KEY} AS k, CAST(json -> '$.out' AS VARCHAR) AS out, json ->> '$.out.errors' AS errors FROM {source(outputs / (name + '.jsonl.gz'))}")
    print(f"  -- {name}")
    show("""
        SELECT CASE WHEN a.k IS NULL THEN 'only in other' WHEN o.k IS NULL THEN 'only in r1-A'
                    WHEN a.out = o.out THEN 'identical'
                    WHEN a.errors IS NULL AND o.errors IS NOT NULL THEN 'r1-A translated, other error'
                    WHEN a.errors IS NOT NULL AND o.errors IS NULL THEN 'r1-A error, other translated'
                    WHEN a.errors IS NOT NULL THEN 'both errors, different text'
                    ELSE 'both translated, different objects' END AS category, count(*) AS n
        FROM a FULL JOIN o USING (k) GROUP BY ALL ORDER BY n DESC
    """)
    show("SELECT a.source, left(coalesce(a.errors, a.out), 90), left(coalesce(o.errors, o.out), 90) FROM a JOIN o USING (k) WHERE a.out IS DISTINCT FROM o.out LIMIT 8")

print("\n== 3. one VRS id per variation and sequence (r1-A)")
con.sql("""
    CREATE OR REPLACE TABLE groups AS
    SELECT variation_id, sequence, any_value(vrs_class) AS vrs_class, count(*) AS n,
           count(DISTINCT vrs_id) AS ids, list(DISTINCT fmt ORDER BY fmt) AS fmts,
           bool_or(precedence = 4) AS has_loc_hgvs, bool_or(issue IS NOT NULL) AS has_issue
    FROM a WHERE errors IS NULL AND vrs_id IS NOT NULL AND sequence IS NOT NULL
    GROUP BY ALL
""")
show("""
    SELECT vrs_class, count(*) AS groups, count(*) FILTER (WHERE n > 1) AS with_2_or_more,
           count(*) FILTER (WHERE n > 1 AND ids = 1) AS agree, count(*) FILTER (WHERE ids > 1) AS disagree,
           sum(n) FILTER (WHERE n > 1) AS records_compared
    FROM groups GROUP BY ALL ORDER BY groups DESC
""")
print("  disagreeing groups by the formats they hold:")
show("""
    SELECT vrs_class, fmts, has_loc_hgvs, count(*) AS groups FROM groups WHERE ids > 1
    GROUP BY ALL ORDER BY groups DESC LIMIT 15
""")
print("  examples:")
show("""
    WITH g AS (SELECT variation_id, sequence FROM groups WHERE ids > 1 USING SAMPLE 6 ROWS)
    SELECT a.variation_id, a.vrs_class, a.fmt, a.precedence, a.source, a.vrs_id
    FROM a JOIN g USING (variation_id, sequence) WHERE a.errors IS NULL ORDER BY a.variation_id, a.precedence
""", limit=40)

print("\n== 4. input issue against output errors (r1-A)")
show("""
    SELECT coalesce(issue, '(no issue)') AS issue,
           count(*) AS records,
           count(*) FILTER (WHERE errors IS NULL) AS translated,
           count(*) FILTER (WHERE errors IS NOT NULL) AS errors
    FROM a GROUP BY ALL ORDER BY records DESC
""")
print("  errors on records with no issue, by error kind:")
show("""
    SELECT regexp_extract(errors, '^([A-Za-z]+(Error|Exception)?)', 1) AS kind,
           left(regexp_replace(errors, '[A-Z]{2}_[0-9]+\\.[0-9]+|[0-9]+', 'N', 'g'), 80) AS shape,
           count(*) AS n
    FROM a WHERE issue IS NULL AND errors IS NOT NULL GROUP BY ALL ORDER BY n DESC LIMIT 15
""")
