"""check_new_sources.py for a remote store: gtars' open_remote defers the sequence
index, and alias lookups return None until it is loaded, which the Python
binding cannot request directly. Load it here, through load_sequence on a digest
that does not exist, before vrs-python uses the store.

Usage: check_remote_sources.py HARNESS_DIR REFGET_URI UTA_URL SAMPLES_JSON
"""

import runpy
import sys
from pathlib import Path

import gtars.refget

_RealStore = gtars.refget.RefgetStore


class _IndexLoadingStore:
    open_local = staticmethod(_RealStore.open_local)

    @staticmethod
    def open_remote(cache_dir, url):
        store = _RealStore.open_remote(cache_dir, url)
        try:
            store.load_sequence("A" * 32)
        except OSError:
            pass  # "Sequence not found in store", raised after the index loads
        return store


gtars.refget.RefgetStore = _IndexLoadingStore
runpy.run_path(str(Path(__file__).with_name("check_new_sources.py")), run_name="__main__")
