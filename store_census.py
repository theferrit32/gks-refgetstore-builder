#!/usr/bin/env python
"""Read-only primitives for describing what a RefgetStore actually holds.

Everything that needs to *look at* a store rather than build one goes through
here: the lock writer's ``outputs`` census, ``verify``, ``repair``,
``provenance.py``, and ``tools/compare_store_mutation.py``. Keeping them in one
module means the answer to "what is in the store" is computed one way, so a
lock written today and a verify run tomorrow cannot disagree by construction.

Nothing here mutates a store, writes a file, or touches the network.

Two cost notes that shape the callers, both measured on the production store
(1,779,052 sequences / 235 collections):

* :func:`sequence_digests` is ~2.5s and :func:`digest_root` ~1.0s over that set,
  so a full root pair is ~4s. Cheap enough to compute on every lock write.
* :func:`on_disk_sequence_digests` is one ``scandir`` pass (~5.8s) instead of
  1.8M ``Path.exists()`` calls, which is the difference between 6 seconds and
  several minutes.
"""

from __future__ import annotations

import base64
import hashlib
import os
from collections.abc import Iterable, Iterator
from pathlib import Path

SQ_PREFIX = "SQ."
SEQ_SUFFIX = ".seq"
# gtars shards sequence payloads into one directory per 2-character digest
# prefix. The layout is gtars', not ours; it is spelled here because verify has
# to find a payload file without asking the store to load it.
SHARD_PREFIX_LEN = 2


def strip_sq(digest: str) -> str:
    """Bare sha512t24u, with or without the ``SQ.`` refget prefix."""
    return digest[len(SQ_PREFIX):] if digest.startswith(SQ_PREFIX) else digest


def digest_root(digests: Iterable[str]) -> str:
    """Canonical digest of a *set* of digests: ``sha256:<hex>``.

    Defined so a reimplementation cannot drift from it: deduplicate, strip any
    ``SQ.`` prefix, sort as ASCII, then feed each digest followed by a single
    ``\\n`` to SHA-256. Order-independent and duplicate-independent by
    construction, so it answers exactly one question -- does this store hold
    the same *set* of digests the lock recorded -- and any single addition or
    removal flips it.
    """
    hasher = hashlib.sha256()
    for digest in sorted({strip_sq(d) for d in digests}):
        hasher.update(digest.encode("ascii"))
        hasher.update(b"\n")
    return "sha256:" + hasher.hexdigest()


def collection_census(store) -> dict[str, int]:
    """Every collection in the store as ``digest -> n_sequences``.

    Enumerated from the store itself, never synthesized from a lock or from
    build provenance: the census is what makes ``outputs`` a statement about
    reality rather than a restatement of intent.
    """
    out: dict[str, int] = {}
    page = 0
    while True:
        result = store.list_collections(page=page, page_size=100)
        for meta in result["results"]:
            out[strip_sq(meta.digest)] = int(meta.n_sequences)
        pagination = result["pagination"]
        if (pagination["page"] + 1) * pagination["page_size"] >= pagination["total"]:
            return out
        page += 1


def collection_members(store, collection_digest: str) -> list[str]:
    """Bare sequence digests belonging to one collection, in level-2 order.

    Loads the collection: a lazily opened store holds collections as metadata
    stubs whose level-2 arrays are not yet resolvable.
    """
    store.load_collection(collection_digest)
    level2 = store.get_collection_level2(collection_digest)
    return [strip_sq(seq_id) for seq_id in level2["sequences"]]


def collection_named_members(store, collection_digest: str) -> list[tuple[str, str]]:
    """``(record name, bare digest)`` pairs for one collection."""
    store.load_collection(collection_digest)
    level2 = store.get_collection_level2(collection_digest)
    return [
        (name, strip_sq(seq_id))
        for name, seq_id in zip(level2["names"], level2["sequences"], strict=True)
    ]


def sequence_digests(store) -> set[str]:
    """Every sequence digest the store indexes.

    This is the store's own index, not a union over collection membership. The
    two are compared in verify: a digest indexed but unreachable from any
    collection is an orphan.
    """
    return {strip_sq(meta.sha512t24u) for meta in store.list_sequences()}


def sequence_metadata_by_digest(store) -> dict[str, object]:
    """``digest -> SequenceMetadata`` for the whole store, in one pass.

    Deep verification needs each sequence's name, alphabet, and length to
    render an actionable finding; asking per digest would be 1.8M lookups.
    """
    return {strip_sq(meta.sha512t24u): meta for meta in store.list_sequences()}


def seq_path(store_dir: Path, digest: str) -> Path:
    """On-disk payload path for one sequence digest."""
    bare = strip_sq(digest)
    return Path(store_dir) / "sequences" / bare[:SHARD_PREFIX_LEN] / f"{bare}{SEQ_SUFFIX}"


def on_disk_sequence_digests(store_dir: Path) -> set[str]:
    """Digests with a payload file under ``<store_dir>/sequences/``.

    One ``scandir`` walk of the shard directories. Deliberately not
    ``rglob``/``exists`` per digest: at 1.8M sequences that is the difference
    between one pass and 1.8M stat calls.
    """
    root = Path(store_dir) / "sequences"
    found: set[str] = set()
    if not root.is_dir():
        return found
    with os.scandir(root) as shards:
        for shard in shards:
            if not shard.is_dir():
                continue
            with os.scandir(shard.path) as entries:
                for entry in entries:
                    if entry.name.endswith(SEQ_SUFFIX):
                        found.add(entry.name[: -len(SEQ_SUFFIX)])
    return found


def stream_sequence_chunks(store, digest: str) -> Iterator[str]:
    """Chunks of a stored sequence, without materializing the whole payload.

    ``SequenceStream.read_all`` would return a 250 Mbp chromosome as one Python
    string; iterating keeps deep verification at O(1) memory per sequence.
    """
    yield from store.stream_sequence(strip_sq(digest))


def redigest_sequence(store, digest: str) -> str:
    """Re-compute the sha512t24u of what the store actually decodes.

    The digest the store is *keyed* by was computed from the source bytes at
    ingest. This recomputes it from the bytes the store can serve today, so the
    two disagreeing means the stored payload no longer represents the sequence
    its digest promises -- whether from bit rot or from a lossy encoder.
    """
    hasher = hashlib.sha512()
    for chunk in stream_sequence_chunks(store, digest):
        hasher.update(chunk.encode("ascii") if isinstance(chunk, str) else chunk)
    return base64.urlsafe_b64encode(hasher.digest()[:24]).decode("ascii")


def store_counts(store) -> tuple[int, int]:
    """``(n_sequences, n_collections)`` as integers.

    ``store.stats()`` returns these as strings; every caller wants numbers, and
    a lock that recorded ``"1779052"`` instead of ``1779052`` compares unequal
    to itself after a round trip.
    """
    stats = store.stats()
    return int(stats["n_sequences"]), int(stats["n_collections"])
