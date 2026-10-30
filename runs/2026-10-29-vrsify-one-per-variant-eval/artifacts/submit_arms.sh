#!/usr/bin/env bash
# Submit one vrsify Batch job per comparison arm.
#
# Usage: submit_arms.sh IMAGE INPUT_URI OUTPUT_PREFIX ARM...
#   ARM is A, B, C, Cfd, Cg or D; each job writes OUTPUT_PREFIX-ARM/vi-vrs-python.jsonl.gz
#   A   RefgetStore archive, extracted; whole-sequence loads
#   B   RefgetStore remote, per-worker caches; whole-sequence loads
#   C   biocommons SeqRepo 2024-12-20 archive, extracted
#   Cg  SeqRepo built from the genomic RefgetStore (same sequences and aliases), extracted
#   D   RefgetStore archive, extracted; range reads from the .seq files
# Prints "ARM JOB_NAME" per submitted job.

set -euo pipefail

image="$1"
input="$2"
prefix="$3"
shift 3

# ARCHIVES=gcs reads everything from us-east1 GCS instead of R2: the archives (arms A, C, Cg, D)
# and, since 2026-10-29, arm B's remote store, from the public clingen-public mirror (gtars reads
# a remote store over unauthenticated HTTPS). Runs before 2026-10-29 read arm B from R2 in both
# modes. The genomic SeqRepo (arm Cg) is published only to GCS.
# STORE=genomic-2026-10-29 (default) adds the RefSeq biological-region NG_ records to
# genomic-2026-10-27; the RefgetStore arms and Cg hold the same sequences in either version.
# In GCS mode, 2026-10-29's store archive comes from clingen-public; 2026-10-27's from the
# internal static bucket, where it was staged.
store_name="${STORE:-genomic-2026-10-29}"
seqrepo_genomic_archive=gs://clingen-public/seqrepo/$store_name.tar.gz
if [ "${ARCHIVES:-r2}" = gcs ]; then
    if [ "$store_name" = genomic-2026-10-27 ]; then
        refget_archive=gs://clingen-dev-gke-internal-static/refgetstore/$store_name.tar.gz
    else
        refget_archive=gs://clingen-public/refgetstore/$store_name.tar.gz
    fi
    seqrepo_archive=gs://clingen-dev-gke-internal-static/seqrepo_2024-12-20.tar.gz
    refget_remote=https://storage.googleapis.com/clingen-public/refgetstore/$store_name
else
    refget_archive=https://static.ferriter.dev/refgetstore/$store_name.tar.gz
    seqrepo_archive=https://static.ferriter.dev/seqrepo/seqrepo_2024-12-20.tar.gz
    refget_remote=https://static.ferriter.dev/refgetstore/$store_name
fi

cd ~/dev/clinvar-ingest/.claude/worktrees/320-vrsify
for arm in "$@"; do
    case "$arm" in
        A) store=(--refget-store "$refget_archive") ;;
        B) store=(--refget-store "$refget_remote") ;;
        # Cfd: images from 12774e5 on keep SeqRepo's files open (SEQREPO_FD_CACHE_MAXSIZE=none)
        C|Cfd) store=(--seqrepo "$seqrepo_archive") ;;
        Cg) store=(--seqrepo "$seqrepo_genomic_archive") ;;
        D) store=(--refget-store "$refget_archive" --refget-reads ranges) ;;
        *) echo "unknown arm $arm" >&2; exit 2 ;;
    esac
    # WORKERS=N sets vrsify's worker processes (default: one per CPU)
    if [ -n "${WORKERS:-}" ]; then
        store+=(--workers "$WORKERS")
    fi
    # GCP labels on the job, VM and disk for billing: experiment, arm, the store version and
    # worker count, plus LABELS="key=value ..." (e.g. "round=r2 input=all-genomic")
    # MACHINE sets the VM type (default c4-standard-8). From 2026-10-29 the all-genomic runs use
    # c4-highmem-8: arms A and B, which keep whole sequences per worker, exceed c4-standard-8's
    # 27 GiB task memory there.
    machine="${MACHINE:-c4-standard-8}"
    labels=(--label experiment=vrsify-store-comparison --label "arm=$(echo "$arm" | tr '[:upper:]' '[:lower:]')"
            --label "store=$store_name" --label "workers=${WORKERS:-all}" --label "machine=$machine")
    for label in ${LABELS:-}; do
        labels+=(--label "$label")
    done
    job=$(GOOGLE_CLOUD_QUOTA_PROJECT=clingen-dev .venv/bin/clinvar-ingest vrsify-submit 2026-10-24 \
        --project clingen-dev --region us-east1 --image "$image" --boot-disk-size 100 \
        --machine-type "$machine" "${labels[@]}" -- \
        --input-uri "$input" --output-uri "$prefix-$arm/vi-vrs-python.jsonl.gz" "${store[@]}" 2>&1 \
        | grep -o "jobs/vrsify-[0-9-]*" | sed 's#jobs/##')
    echo "$arm $job"
    # Job names carry a seconds timestamp
    sleep 2
done
