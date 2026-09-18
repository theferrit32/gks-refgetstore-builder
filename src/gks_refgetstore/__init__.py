"""Build a gtars RefgetStore archive with the aliases GA4GH tools need.

Deliberately empty of re-exports. Importing a submodule here would make
``import gks_refgetstore`` pull in gtars and the whole engine, which is exactly
what the CLI's lazy dispatch avoids -- ``gks-refgetstore --help`` should not pay
for the ingest path. Import the submodule you need:

    from gks_refgetstore import build_lock, store_census
    from gks_refgetstore.sources import load_config, resolve_sources
"""
