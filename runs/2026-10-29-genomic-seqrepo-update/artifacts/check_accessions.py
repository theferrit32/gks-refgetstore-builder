"""Resolve accessions through the 2026-10-27 and 2026-10-29 genomic SeqRepo and store.

Prints, for each accession, the refget accession vrs-python's SeqRepo and
RefgetStore data proxies derive from each pair, and exits non-zero unless every
accession resolves identically in both 2026-10-29 backends.

Usage: check_accessions.py ACCESSION...
"""

import sys

from ga4gh.vrs.dataproxy import create_dataproxy

ROOT = "/Volumes/RefgetStores"
proxies = {
    f"{kind}-{date}": create_dataproxy(uri)
    for date in ("2026-10-27", "2026-10-29")
    for kind, uri in (("seqrepo", f"seqrepo+file://{ROOT}/seqrepo/genomic-{date}"),
                      ("refget", f"refget+file://{ROOT}/genomic-{date}"))
}
ok = True
print("accession\t" + "\t".join(proxies))
for accession in sys.argv[1:]:
    got = {name: proxy.derive_refget_accession(accession) for name, proxy in proxies.items()}
    print(accession + "\t" + "\t".join(str(v) for v in got.values()))
    new = got["seqrepo-2026-10-29"]
    ok &= new is not None and new == got["refget-2026-10-29"]
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
