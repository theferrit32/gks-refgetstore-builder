#!/usr/bin/env python3
"""Create and validate compact, reproducible run records."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA = "gks-refgetstore-run/1"
MANIFEST = "run-manifest.json"
MAX_ARTIFACT_BYTES = 1024 * 1024
LABEL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_text(argv: list[str], cwd: Path | None = None) -> tuple[int, str]:
    try:
        result = subprocess.run(
            argv, cwd=cwd, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False,
        )
    except OSError as exc:
        return 127, str(exc)
    return result.returncode, result.stdout.strip()


def repo_root(start: Path) -> Path:
    code, output = run_text(["git", "rev-parse", "--show-toplevel"], start)
    return Path(output).resolve() if code == 0 else start.resolve()


def git_metadata(root: Path) -> dict[str, Any]:
    commit_code, commit = run_text(["git", "rev-parse", "HEAD"], root)
    status_code, status = run_text(["git", "status", "--porcelain=v1"], root)
    diff = subprocess.run(
        ["git", "diff", "--binary", "HEAD", "--"], cwd=root,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
    ).stdout
    return {
        "commit": commit if commit_code == 0 else None,
        "dirty": status_code == 0 and bool(status),
        "status_text": status if status_code == 0 else None,
        "diff_sha256": hashlib.sha256(diff).hexdigest() if commit_code == 0 else None,
    }


def tool_versions() -> dict[str, str]:
    versions = {"python": platform.python_version()}
    for name, argv in (
        ("git", ["git", "--version"]),
        ("uv", ["uv", "--version"]),
        ("gtars", [sys.executable, "-c", "import gtars; print(gtars.__version__)"]),
    ):
        code, output = run_text(argv)
        if code == 0:
            versions[name] = output.splitlines()[0]
    return versions


def manifest_path(run_dir: Path) -> Path:
    return run_dir.resolve() / MANIFEST


def load_manifest(run_dir: Path) -> tuple[Path, dict[str, Any]]:
    path = manifest_path(run_dir)
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot read {path}: {exc}") from exc
    return path, data


def write_manifest(path: Path, data: dict[str, Any]) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    os.replace(temporary, path)


def relative_record(path: Path, run_dir: Path, role: str) -> dict[str, Any]:
    return {
        "role": role,
        "path": path.relative_to(run_dir).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def cmd_init(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    if run_dir.exists() and any(run_dir.iterdir()):
        raise SystemExit(f"run directory is not empty: {run_dir}")
    root = repo_root(Path.cwd())
    git = git_metadata(root)
    if args.require_clean and git["dirty"]:
        raise SystemExit("--require-clean requested, but the Git worktree is dirty")
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "logs").mkdir(exist_ok=True)
    quality = "reconstructed" if args.reconstructed else "captured"
    data: dict[str, Any] = {
        "schema": SCHEMA,
        "run_id": run_dir.name,
        "kind": args.kind,
        "title": args.title,
        "status": "initialized",
        "timestamps": {"initialized_utc": utc_now()},
        "baseline_run_id": args.baseline,
        "evidence": {"quality": quality, "missing": []},
        "git": git,
        "tool_versions": tool_versions(),
        "commands": [],
        "artifacts": [],
        "result_summary": {},
        "validations": [],
        "publication": None,
    }
    write_manifest(run_dir / MANIFEST, data)
    print(run_dir / MANIFEST)
    return 0


def cmd_exec(args: argparse.Namespace) -> int:
    if not LABEL_RE.fullmatch(args.label):
        raise SystemExit("label must contain only letters, digits, dot, underscore, or dash")
    argv = list(args.command)
    if argv and argv[0] == "--":
        argv.pop(0)
    if not argv:
        raise SystemExit("exec requires a command after --")
    run_dir = args.run_dir.resolve()
    path, data = load_manifest(run_dir)
    log_path = run_dir / "logs" / f"{args.label}.log"
    if log_path.exists() or any(c.get("label") == args.label for c in data["commands"]):
        raise SystemExit(f"command label already exists: {args.label}")
    started = utc_now()
    entry = {
        "label": args.label,
        "argv": argv,
        "working_directory": str(Path.cwd().resolve()),
        "started_utc": started,
        "finished_utc": None,
        "exit_code": None,
        "log": f"logs/{args.label}.log",
    }
    data["commands"].append(entry)
    data["status"] = "running"
    write_manifest(path, data)
    exit_code = 127
    try:
        with log_path.open("wb") as log:
            process = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            assert process.stdout is not None
            for chunk in iter(lambda: process.stdout.read(64 * 1024), b""):
                log.write(chunk)
                log.flush()
                sys.stdout.buffer.write(chunk)
                sys.stdout.buffer.flush()
            exit_code = process.wait()
    except OSError as exc:
        message = f"failed to execute {argv[0]}: {exc}\n".encode()
        log_path.write_bytes(message)
        sys.stderr.buffer.write(message)
    finally:
        entry["finished_utc"] = utc_now()
        entry["exit_code"] = exit_code
        entry["log_bytes"] = log_path.stat().st_size
        entry["log_sha256"] = sha256_file(log_path)
        if exit_code != 0:
            data["status"] = "failed"
        write_manifest(path, data)
    return exit_code


def cmd_add(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    path, data = load_manifest(run_dir)
    source = args.path.resolve()
    if not source.is_file():
        raise SystemExit(f"artifact is not a file: {source}")
    size = source.stat().st_size
    if size > MAX_ARTIFACT_BYTES and not args.allow_large:
        raise SystemExit(
            f"artifact is {size} bytes; limit is {MAX_ARTIFACT_BYTES}. "
            "Use --allow-large only after reviewing the summaries-only policy."
        )
    name = args.copy_as or source.name
    if Path(name).name != name or name in {"", ".", ".."}:
        raise SystemExit("--copy-as must be a single file name")
    destination = run_dir / "artifacts" / name
    destination.parent.mkdir(exist_ok=True)
    if destination.exists():
        raise SystemExit(f"artifact already exists: {destination}")
    shutil.copy2(source, destination)
    data["artifacts"].append(relative_record(destination, run_dir, args.role))
    write_manifest(path, data)
    print(destination)
    return 0


def cmd_finalize(args: argparse.Namespace) -> int:
    path, data = load_manifest(args.run_dir)
    data["status"] = args.status
    data.setdefault("timestamps", {})["finalized_utc"] = utc_now()
    write_manifest(path, data)
    return 0


def validate_one(run_dir: Path) -> list[str]:
    errors: list[str] = []
    path = manifest_path(run_dir)
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return [f"{path}: {exc}"]
    for field in (
        "schema", "run_id", "kind", "title", "status", "timestamps",
        "evidence", "git", "tool_versions", "commands", "artifacts",
        "result_summary", "validations",
    ):
        if field not in data:
            errors.append(f"missing field: {field}")
    if data.get("schema") != SCHEMA:
        errors.append(f"unexpected schema: {data.get('schema')!r}")
    if data.get("run_id") != run_dir.resolve().name:
        errors.append("run_id does not match directory name")
    if data.get("evidence", {}).get("quality") not in {"captured", "reconstructed"}:
        errors.append("evidence.quality must be captured or reconstructed")
    for item in data.get("artifacts", []):
        rel = Path(item.get("path", ""))
        target = (run_dir.resolve() / rel).resolve()
        if rel.is_absolute() or run_dir.resolve() not in target.parents:
            errors.append(f"unsafe artifact path: {rel}")
            continue
        if not target.is_file():
            errors.append(f"missing artifact: {rel}")
            continue
        if target.stat().st_size != item.get("bytes"):
            errors.append(f"size mismatch: {rel}")
        if sha256_file(target) != item.get("sha256"):
            errors.append(f"SHA-256 mismatch: {rel}")
    for command in data.get("commands", []):
        log = command.get("log")
        if log:
            log_path = run_dir.resolve() / log
            if not log_path.is_file():
                errors.append(f"missing command log: {log}")
            else:
                if log_path.stat().st_size != command.get("log_bytes"):
                    errors.append(f"command log size mismatch: {log}")
                if sha256_file(log_path) != command.get("log_sha256"):
                    errors.append(f"command log SHA-256 mismatch: {log}")
        if command.get("exit_code") is None and data.get("status") != "running":
            errors.append(f"command has no exit code: {command.get('label')}")
    return errors


def cmd_validate(args: argparse.Namespace) -> int:
    failed = False
    for run_dir in args.run_dirs:
        errors = validate_one(run_dir)
        if errors:
            failed = True
            print(f"FAIL {run_dir}")
            for error in errors:
                print(f"  {error}")
        else:
            print(f"PASS {run_dir}")
    return int(failed)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    sub = result.add_subparsers(dest="subcommand", required=True)
    init = sub.add_parser("init", help="initialize a run directory")
    init.add_argument("run_dir", type=Path)
    init.add_argument("--kind", required=True)
    init.add_argument("--title", required=True)
    init.add_argument("--baseline")
    mode = init.add_mutually_exclusive_group()
    mode.add_argument("--require-clean", action="store_true")
    mode.add_argument("--reconstructed", action="store_true")
    init.set_defaults(func=cmd_init)

    execute = sub.add_parser("exec", help="execute and log a command")
    execute.add_argument("run_dir", type=Path)
    execute.add_argument("--label", required=True)
    execute.add_argument("command", nargs="+")
    execute.set_defaults(func=cmd_exec)

    add = sub.add_parser("add", help="copy a small artifact into a run record")
    add.add_argument("run_dir", type=Path)
    add.add_argument("--role", required=True)
    add.add_argument("path", type=Path)
    add.add_argument("--copy-as")
    add.add_argument("--allow-large", action="store_true")
    add.set_defaults(func=cmd_add)

    finalize = sub.add_parser("finalize", help="finalize a run")
    finalize.add_argument("run_dir", type=Path)
    finalize.add_argument("--status", choices=("succeeded", "failed"), required=True)
    finalize.set_defaults(func=cmd_finalize)

    validate = sub.add_parser("validate", help="validate one or more run records")
    validate.add_argument("run_dirs", type=Path, nargs="+")
    validate.set_defaults(func=cmd_validate)
    return result


if __name__ == "__main__":
    arguments = parser().parse_args()
    raise SystemExit(arguments.func(arguments))
