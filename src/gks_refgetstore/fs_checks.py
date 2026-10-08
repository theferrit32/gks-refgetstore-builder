"""Filesystem preconditions for writing a RefgetStore.

gtars shards sequence payloads into ``sequences/<first two digest
characters>/``. Digests are base64url, so ``zz``, ``Zz``, ``zZ`` and ``ZZ`` are
four different shards. On a case-insensitive filesystem -- the default for
macOS APFS and Windows NTFS -- they are one directory, named by whichever was
created first, and payloads land in a directory whose name does not match their
digest prefix.

Nothing local notices: lookups on that filesystem ignore case, so the store
reads and verifies cleanly. The damage appears once the store leaves it. Copied
to a case-sensitive filesystem, or uploaded to object storage (whose keys are
case-sensitive), most payloads sit at a path no reader will request -- 64.7% of
them in a full build. Checking up front is the only point at which this is
cheap.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

ALLOW_FLAG = "--allow-case-insensitive-fs"


class CaseInsensitiveFilesystemError(RuntimeError):
    """The store directory is on a filesystem that ignores filename case."""


def _nearest_existing(path: Path) -> Path:
    """``path`` if it exists, else its closest existing ancestor.

    A directory that does not exist yet will be created on the same filesystem
    as that ancestor, so probing the ancestor answers the same question without
    creating anything the caller did not ask for.
    """
    path = Path(path).absolute()
    while not path.exists():
        if path.parent == path:
            break
        path = path.parent
    return path


def is_case_sensitive(directory: Path) -> bool:
    """Whether the filesystem holding ``directory`` distinguishes names by case.

    Creates a uniquely named probe file with lowercase letters in its name,
    asks whether the uppercase spelling refers to the same file, and removes
    the probe. ``os.path.samefile`` rather than ``exists`` keeps an unrelated
    file that happens to have the uppercase name from being mistaken for a
    case-insensitive match.
    """
    probe_dir = _nearest_existing(directory)
    fd, name = tempfile.mkstemp(prefix=".gks-case-probe-", dir=probe_dir)
    os.close(fd)
    probe = Path(name)
    try:
        variant = probe.with_name(probe.name.upper())
        return not (variant.exists() and os.path.samefile(probe, variant))
    finally:
        probe.unlink()


def require_case_sensitive(store_dir: Path, *, allow: bool = False) -> None:
    """Refuse to write a store to a case-insensitive filesystem.

    ``allow`` skips the refusal for a store that will only ever be read in
    place on the same filesystem; such a store must not be copied to a
    case-sensitive filesystem or uploaded.
    """
    if allow or is_case_sensitive(store_dir):
        return
    raise CaseInsensitiveFilesystemError(
        f"{store_dir} is on a case-insensitive filesystem. gtars shards sequence "
        "payloads by case-sensitive digest prefix, so here prefixes that differ "
        "only in case share one directory, and a copy or upload of the store "
        "puts most payloads at paths readers will not find.\n"
        "Put the store on a case-sensitive filesystem. On macOS, create a "
        "case-sensitive APFS volume, e.g.\n"
        "    diskutil apfs addVolume <container> APFSX RefgetStores\n"
        "and pass --store-dir /Volumes/RefgetStores/store. "
        f"To write here anyway for local-only use, pass {ALLOW_FLAG}."
    )


def preflight_store_dir(args) -> None:
    """CLI preflight for commands that write store payloads.

    Exits with the explanation rather than a traceback. Reads
    ``args.allow_case_insensitive_fs`` when the command defines it.
    """
    try:
        require_case_sensitive(
            args.store_dir,
            allow=getattr(args, "allow_case_insensitive_fs", False),
        )
    except CaseInsensitiveFilesystemError as exc:
        raise SystemExit(f"error: {exc}") from None
