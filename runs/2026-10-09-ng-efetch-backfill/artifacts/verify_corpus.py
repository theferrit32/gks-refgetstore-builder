"""Check the efetch'd RefSeqGene history corpus before and after ingest.

1. Every record round-trips through a fresh in-memory gtars store: the decoded
   payload re-digests to the digest of the source sequence.
2. Every record seqrepo also has (NCBI namespace) has the same digest there.
3. Every record resolves in the built store under refseq:<accession>.

Usage: uv run python verify_corpus.py FASTA SEQREPO_ALIASES_SQLITE STORE_DIR
"""

import gzip
import sqlite3
import sys

from gtars.refget import RefgetStore, sha512t24u_digest

fasta, aliases_db, store_dir = sys.argv[1:4]

records: dict[str, str] = {}
name, parts = None, []
for line in gzip.open(fasta, "rt"):
    if line.startswith(">"):
        if name:
            records[name] = "".join(parts)
        name, parts = line[1:].split()[0], []
    else:
        parts.append(line.strip())
records[name] = "".join(parts)
digests = {n: sha512t24u_digest(s.upper().encode()) for n, s in records.items()}
alphabet = "".join(sorted({c for s in records.values() for c in set(s)}))
print(f"records={len(records)} bases={sum(map(len, records.values()))} alphabet={alphabet}")

memory = RefgetStore.in_memory()
memory.set_quiet(True)
meta, _ = memory.add_sequence_collection_from_fasta(fasta)
failures = sum(
    (rec := memory.get_sequence(d)) is None or sha512t24u_digest(rec.decode().encode()) != d
    for d in digests.values()
)
print(f"in-memory collection={meta.digest} n_sequences={meta.n_sequences} "
      f"round_trip_failures={failures}")

con = sqlite3.connect(f"file:{aliases_db}?mode=ro", uri=True)
seqrepo = dict(con.execute(
    "select alias, seq_id from seqalias where namespace = 'NCBI' "
    "and alias like 'NG\\_%' escape '\\'"))
shared = [n for n in digests if n in seqrepo]
print(f"in_seqrepo={len(shared)} seqrepo_digest_mismatches="
      f"{sum(digests[n] != seqrepo[n] for n in shared)} "
      f"not_in_seqrepo={len(digests) - len(shared)}")

store = RefgetStore.open_local(store_dir)
store.set_quiet(True)
unresolved = [n for n, d in digests.items()
              if (m := store.get_sequence_metadata_by_alias("refseq", n)) is None
              or m.sha512t24u != d]
print(f"store_alias_unresolved_or_wrong={len(unresolved)} {unresolved[:5]}")
sys.exit(1 if failures or unresolved or any(digests[n] != seqrepo[n] for n in shared) else 0)
