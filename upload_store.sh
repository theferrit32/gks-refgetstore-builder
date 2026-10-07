#!/usr/bin/env bash
# Upload a built RefgetStore (the store/ directory from gks-refgetstore build) to object
# storage as individual objects, preserving the relative tree.
#
# A gtars RefgetStore is read lazily, object-by-object (range requests against
# sequences/, collections/, aliases/, the .rgsi/.rgci indexes, and rgstore.json).
# So each file must become its own object -- NEVER tar/repackage the store.
#
# Generic across backends: the destination is any rclone remote (Cloudflare R2,
# AWS S3, GCS). Credentials come from rclone's own config, never this repo.
set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# --- defaults (flags override env override these) ---------------------------
REMOTE="${RGS_REMOTE:-}"                 # required, no default
BUCKET="${RGS_BUCKET:-}"                 # required, no default
DATE="$(date -u +%Y-%m-%d)"
PREFIX="${RGS_PREFIX:-refgetstore/${DATE}}"
STORE_DIR="${RGS_STORE_DIR:-${HERE}/store}"
MODE="sync"                              # sync (mirror) by default; --copy for additive
MAKE_MANIFEST=1
RUN_CHECK=1
DRY_RUN=0
TRANSFERS="${RGS_TRANSFERS:-32}"
CHECKERS="${RGS_CHECKERS:-64}"

usage() {
    cat <<'EOF'
Usage: upload_store.sh --bucket NAME [options]

Uploads the RefgetStore directory to REMOTE:BUCKET/PREFIX as individual objects.

Options:
  --bucket NAME       Target bucket (REQUIRED; also RGS_BUCKET). No default, so
                      nothing can upload to the wrong place by accident.
  --remote NAME       rclone remote (REQUIRED; also RGS_REMOTE).
  --prefix PREFIX     Object key prefix (default: refgetstore/<UTC-date>;
                      RGS_PREFIX). Set explicitly to control the full path.
  --store-dir DIR     Store directory to upload (default: ./store; RGS_STORE_DIR).
  --copy              Use 'rclone copy' (additive) instead of the default
                      'rclone sync' (mirror, deletes extra remote objects).
  --no-manifest       Skip writing manifest.json into the store dir.
  --no-check          Skip the post-upload 'rclone check' verification pass.
  --dry-run           Pass --dry-run to rclone (no objects written/deleted).
  --transfers N       Parallel transfers (default: 32; RGS_TRANSFERS).
  --checkers N        Parallel checkers (default: 64; RGS_CHECKERS).
  -h, --help          Show this help.

Credentials: read by rclone from ~/.config/rclone/rclone.conf (or RCLONE_CONFIG_*
env vars / --config). No secrets live in this repo.
EOF
}

# --- arg parsing ------------------------------------------------------------
while [ $# -gt 0 ]; do
    case "$1" in
        --bucket)     BUCKET="$2"; shift 2 ;;
        --remote)     REMOTE="$2"; shift 2 ;;
        --prefix)     PREFIX="$2"; shift 2 ;;
        --store-dir)  STORE_DIR="$2"; shift 2 ;;
        --copy)       MODE="copy"; shift ;;
        --no-manifest) MAKE_MANIFEST=0; shift ;;
        --no-check)   RUN_CHECK=0; shift ;;
        --dry-run)    DRY_RUN=1; shift ;;
        --transfers)  TRANSFERS="$2"; shift 2 ;;
        --checkers)   CHECKERS="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

# --- preflight --------------------------------------------------------------
command -v rclone >/dev/null 2>&1 || die "rclone not found on PATH"
[ -n "$REMOTE" ] || die "--remote is required (or set RGS_REMOTE)"
[ -n "$BUCKET" ] || die "--bucket is required (or set RGS_BUCKET)"
[ -d "$STORE_DIR" ] || die "store dir not found: $STORE_DIR"
[ -f "$STORE_DIR/rgstore.json" ] || \
    die "$STORE_DIR does not look like a RefgetStore (no rgstore.json)"
[ -n "$(ls -A "$STORE_DIR")" ] || die "store dir is empty: $STORE_DIR"
rclone listremotes 2>/dev/null | grep -qx "${REMOTE}:" || \
    die "rclone remote '${REMOTE}:' not configured (see: rclone listremotes)"

DEST="${REMOTE}:${BUCKET}/${PREFIX}"

echo "Store dir : $STORE_DIR"
echo "Dest      : $DEST"
echo "Mode      : rclone ${MODE}$([ "$DRY_RUN" = 1 ] && echo ' (dry-run)')"

# --- manifest / provenance --------------------------------------------------
if [ "$MAKE_MANIFEST" = 1 ]; then
    echo "Generating manifest.json ..."
    STORE_DIR="$STORE_DIR" PREFIX="$PREFIX" REMOTE="$REMOTE" BUCKET="$BUCKET" \
    REPO_DIR="$HERE" UPLOAD_DATE="$DATE" \
    uv run --project "$HERE" python - <<'PY'
import json, os, subprocess, tomllib
from pathlib import Path

store = Path(os.environ["STORE_DIR"])
repo = Path(os.environ["REPO_DIR"])

def git(*args, default=None):
    try:
        return subprocess.check_output(["git", "-C", str(repo), *args],
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return default

# gtars version
try:
    import gtars
    gtars_version = getattr(gtars, "__version__", None)
except Exception:
    gtars_version = None

# source pins from sources.toml
cfg_path = repo / "sources.toml"
assemblies, seqsets = [], []
if cfg_path.exists():
    cfg = tomllib.loads(cfg_path.read_text())
    for a in cfg.get("assembly", []):
        assemblies.append({"namespace": a.get("namespace"),
                           "fasta_url": a.get("fasta_url"),
                           "report_url": a.get("report_url")})
    for s in cfg.get("seqset", []):
        seqsets.append({"name": s.get("name"), "namespace": s.get("namespace"),
                        "url_template": s.get("url_template"),
                        "shard_range": s.get("shard_range")})

# object count + total bytes (exclude the manifest we are about to write)
count, total = 0, 0
manifest_path = store / "manifest.json"
for p in store.rglob("*"):
    if p.is_file() and p != manifest_path:
        count += 1
        total += p.stat().st_size

manifest = {
    "schema": "refgetstore-manifest/1",
    "upload_date": os.environ["UPLOAD_DATE"],
    "destination": f'{os.environ["REMOTE"]}:{os.environ["BUCKET"]}/{os.environ["PREFIX"]}',
    "repo": {
        "name": "gks-refgetstore-builder",
        "commit": git("rev-parse", "HEAD"),
        "describe": git("describe", "--always", "--dirty", "--tags"),
        "dirty": git("status", "--porcelain", default="") != "",
    },
    "gtars_version": gtars_version,
    "object_count": count,
    "total_bytes": total,
    "assemblies": assemblies,
    "seqsets": seqsets,
}
manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2))
PY
fi

# --- upload -----------------------------------------------------------------
RCLONE_OPTS=(--transfers "$TRANSFERS" --checkers "$CHECKERS" --fast-list
             --s3-no-check-bucket --stats-one-line)
# --progress redraws a terminal line and suppresses rclone's periodic logged
# stats, so a detached run (nohup, run_record exec) would log nothing until the
# end. Use it only on a terminal; otherwise log stats every --stats interval
# (default 1m; override with RCLONE_STATS), to RCLONE_LOG_FILE when set.
if [ -t 1 ]; then
    RCLONE_OPTS+=(--progress)
else
    RCLONE_OPTS+=(--stats "${RCLONE_STATS:-1m}" --stats-log-level NOTICE)
fi
[ "$DRY_RUN" = 1 ] && RCLONE_OPTS+=(--dry-run)

echo "Uploading (rclone ${MODE}) ..."
rclone "$MODE" "$STORE_DIR" "$DEST" "${RCLONE_OPTS[@]}"

# --- verify -----------------------------------------------------------------
if [ "$RUN_CHECK" = 1 ] && [ "$DRY_RUN" != 1 ]; then
    echo "Verifying (rclone check) ..."
    rclone check "$STORE_DIR" "$DEST" --fast-list --s3-no-check-bucket
fi

echo
echo "Done. Store available under: $DEST"
