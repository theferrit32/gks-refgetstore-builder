"""Check that a subset RefgetStore is complete and is exactly a subset of a full store.

Usage: check_subset_store.py SUBSET_DIR FULL_DIR FULL_LOCK SUBSET_MANIFEST
"""

import hashlib
import json
import sys
from pathlib import Path

subset, full, lock_path = (Path(p) for p in sys.argv[1:4])
meta = json.loads((subset / "rgstore.json").read_text())
template = meta["seqdata_path_template"]


def seq_path(digest: str) -> str:
    return template.replace("%s2", digest[:2]).replace("%s4", digest[:4]).replace("%s", digest)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def index_digests(path: Path) -> list[str]:
    # RGSI/RGCI: '#' header lines, then tab-separated rows whose digest column is named in the header
    header, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            header = line.lstrip("#").split("\t")
        elif line.strip():
            rows.append(dict(zip(header, line.split("\t"))))
    key = next(k for k in ("sha512t24u", "digest") if k in header)
    return [r[key] for r in rows]


problems = []
seqs = index_digests(subset / meta["sequence_index"])
colls = index_digests(subset / meta["collection_index"])
print(f"sequence index: {len(seqs)} digests; collection index: {len(colls)} digests")

on_disk = {p.stem for p in (subset / "sequences").rglob("*.seq")}
missing = [d for d in seqs if d not in on_disk]
extra = sorted(on_disk - set(seqs))
print(f".seq files: {len(on_disk)}; indexed but missing: {len(missing)}; on disk but unindexed: {len(extra)}")
problems += missing + extra

for d in colls:
    if not (subset / meta["collections_path_template"].replace("%s", d)).exists():
        problems.append(f"collection file missing: {d}")

# Every collection the full lock attributes to a source owned by the subset manifest must be
# present with the same sequence count, and the subset must hold no other collection.
import tomllib
lock = json.loads(lock_path.read_text())
manifest = tomllib.load(open(sys.argv[4], "rb"))
owners = {a["namespace"] for a in manifest.get("assembly", [])} | {s["name"] for s in manifest.get("seqset", [])}
owner_of = {f["cache_path"]: f["owner"] for f in lock["inputs"]["files"]}
expected = {c["digest"]: c["n_sequences"] for c in lock["outputs"]["collections"]
            if any(owner_of.get(p) in owners for p in c["from"])}
counts = {}
for line in (subset / meta["collection_index"]).read_text().splitlines():
    if line.strip() and not line.startswith("#"):
        digest, n = line.split("\t")[:2]
        counts[digest] = int(n)
print(f"collections the full lock attributes to the manifest's {len(owners)} sources: {len(expected)}")
for d, n in expected.items():
    if d not in counts:
        problems.append(f"expected collection missing: {d}")
    elif counts[d] != n:
        problems.append(f"collection {d} has {counts[d]} sequences, lock says {n}")
problems += [f"collection not attributed to the manifest: {d}" for d in counts if d not in expected]
print(f"expected collections present with matching counts: {sum(1 for d, n in expected.items() if counts.get(d) == n)}/{len(expected)}")

differ = 0
for d in seqs:
    a, b = subset / seq_path(d), full / seq_path(d)
    if not b.exists() or a.stat().st_size != b.stat().st_size or sha256(a) != sha256(b):
        differ += 1
        problems.append(f"sequence differs from full store: {d}")
print(f"sequence files byte-identical to the full store: {len(seqs) - differ}/{len(seqs)}")
for d in colls:
    a = subset / meta["collections_path_template"].replace("%s", d)
    b = full / meta["collections_path_template"].replace("%s", d)
    if not b.exists() or sha256(a) != sha256(b):
        problems.append(f"collection file differs from full store: {d}")

print("PASS" if not problems else f"FAIL: {len(problems)} problem(s), e.g. {problems[:5]}")
sys.exit(1 if problems else 0)
