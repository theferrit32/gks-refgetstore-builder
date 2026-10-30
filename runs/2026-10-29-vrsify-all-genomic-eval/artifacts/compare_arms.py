"""Compare vrsify outputs of the comparison arms with arm A, record by record.

Usage: compare_arms.py WORK_DIR A_OUTPUT OTHER_NAME=OTHER_OUTPUT...
Prints, per other arm, how its records compare with A's, and classifies the differences
by the error text on either side. Records are matched on (variation_id, fmt,
assembly_version, source): the all-genomic input has many rows per variation.
"""

import sys

import duckdb

work, a_path, *others = sys.argv[1:]
con = duckdb.connect(f"{work}/compare.duckdb")
con.sql("SET memory_limit = '12GB'")
con.sql(f"SET temp_directory = '{work}/duckdb-tmp'")
KEY = "variation_id, fmt, assembly_version, source"


def load(name: str, path: str) -> None:
    con.sql(f"""
        CREATE OR REPLACE TABLE {name} AS
        SELECT json ->> '$.in.variation_id' AS variation_id,
               json ->> '$.in.vrs_class' AS vrs_class,
               coalesce(json ->> '$.in.fmt', '') AS fmt,
               coalesce(json ->> '$.in.assembly_version', '') AS assembly_version,
               coalesce(json ->> '$.in.source', '') AS source,
               CAST((json -> '$.out') AS VARCHAR) AS out,
               json ->> '$.out.errors' AS errors
        FROM read_json_objects('{path}', format = 'newline_delimited')
    """)


load("a", a_path)
print("A rows:", con.execute("SELECT count(*) FROM a").fetchall()[0][0])
for spec in others:
    name, path = spec.split("=", 1)
    load(name, path)
    rows, only = con.execute(f"""
        SELECT (SELECT count(*) FROM {name}),
               (SELECT count(*) FROM a FULL JOIN {name} o USING ({KEY})
                WHERE a.variation_id IS NULL OR o.variation_id IS NULL)
    """).fetchall()[0]
    print(f"\n== {name}: {rows} rows; records in only one side: {only}")
    for row in con.execute(f"""
        SELECT CASE
                 WHEN a.out IS NOT DISTINCT FROM o.out THEN 'identical'
                 WHEN a.errors IS NULL AND o.errors LIKE 'TimeoutError%' THEN 'A translated, {name} timed out'
                 WHEN a.errors IS NULL AND o.errors IS NOT NULL THEN 'A translated, {name} error'
                 WHEN a.errors IS NOT NULL AND o.errors IS NULL THEN 'A error, {name} translated'
                 WHEN a.errors IS NOT NULL AND o.errors IS NOT NULL THEN 'both errors, different text'
                 ELSE 'both translated, different objects'
               END AS category,
               count(*) AS n
        FROM a JOIN {name} o USING ({KEY})
        GROUP BY ALL ORDER BY n DESC
    """).fetchall():
        print(f"  {row[1]:>9}  {row[0]}")
    for row in con.execute(f"""
        SELECT a.vrs_class, a.source, left(coalesce(a.errors, a.out), 70), left(coalesce(o.errors, o.out), 110)
        FROM a JOIN {name} o USING ({KEY})
        WHERE a.out IS DISTINCT FROM o.out
        ORDER BY a.vrs_class, a.source LIMIT 12
    """).fetchall():
        print("   ", row)
