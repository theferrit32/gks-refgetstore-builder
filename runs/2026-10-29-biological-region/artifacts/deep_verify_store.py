"""Re-digest every sequence a store serves (verify's L2), with no lock.

`gks-refgetstore verify --store --deep` also runs L0 and L1 against a lock, which
a store built as a subset of that lock cannot pass. L2 compares the store with
itself, so it applies to a subset store unchanged.

Usage: deep_verify_store.py STORE_DIR
"""

import sys

from gtars.refget import RefgetStore

from gks_refgetstore import verify

store = RefgetStore.open_local(sys.argv[1])
store.set_quiet(True)
report = verify.verify_store_l2(store, {})
verify.render(report)
sys.exit(1 if any(f.severity == verify.ERROR for f in report.findings) else 0)
