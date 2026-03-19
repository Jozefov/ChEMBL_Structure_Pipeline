#!/bin/bash
# Submit all 9 Enamine REAL processing jobs as a dependency chain.
#
# Each job starts only after the previous one finishes (success or failure).
# This way you submit once and forget — over weeks/months all files get processed.
#
# Usage:
#   # Submit all 9 files as a chain:
#   bash scripts/submit_enamine_chain.sh
#
#   # Submit starting from a specific file (e.g., skip already-done files):
#   bash scripts/submit_enamine_chain.sh --start 3
#
#   # Resume mode (all jobs use --resume flag):
#   bash scripts/submit_enamine_chain.sh --resume
#
#   # Dry run (show what would be submitted):
#   bash scripts/submit_enamine_chain.sh --dry-run

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PBS_SCRIPT="${SCRIPT_DIR}/run_enamine.sh"

# All 9 file stems in order (smallest to largest)
FILE_STEMS=(
    "HAC_11_21"
    "HAC_22_23"
    "HAC_24"
    "HAC_25"
    "HAC_26"
    "HAC_27"
    "HAC_28"
    "HAC_29_38_P1"
    "HAC_29_38_P2"
)

# Parse arguments
START_IDX=0
RESUME_FLAG=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --start)
            START_IDX="$2"
            shift 2
            ;;
        --resume)
            RESUME_FLAG=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 [--start N] [--resume] [--dry-run]"
            exit 1
            ;;
    esac
done

echo "============================================"
echo "Enamine REAL Chain Submission"
echo "============================================"
echo "PBS script: $PBS_SCRIPT"
echo "Starting from index: $START_IDX"
echo "Resume mode: $( [ $RESUME_FLAG -eq 1 ] && echo 'yes' || echo 'no' )"
echo "Dry run: $( [ $DRY_RUN -eq 1 ] && echo 'yes' || echo 'no' )"
echo "============================================"
echo ""

PREV_JOB_ID=""

for i in $(seq $START_IDX $((${#FILE_STEMS[@]} - 1))); do
    STEM="${FILE_STEMS[$i]}"

    # Build qsub command
    QSUB_VARS="-v FILE_STEM=${STEM}"
    if [ $RESUME_FLAG -eq 1 ]; then
        QSUB_VARS="${QSUB_VARS},RESUME=1"
    fi

    QSUB_CMD="qsub ${QSUB_VARS}"

    # Add dependency on previous job (if any)
    if [ -n "$PREV_JOB_ID" ]; then
        QSUB_CMD="${QSUB_CMD} -W depend=afterany:${PREV_JOB_ID}"
    fi

    QSUB_CMD="${QSUB_CMD} -N enamine_${STEM} ${PBS_SCRIPT}"

    if [ $DRY_RUN -eq 1 ]; then
        echo "[${i}] Would submit: $QSUB_CMD"
        PREV_JOB_ID="DRYRUN_${i}"
    else
        echo -n "[${i}] Submitting ${STEM}... "
        JOB_ID=$(eval "$QSUB_CMD")
        PREV_JOB_ID=$(echo "$JOB_ID" | tr -d '[:space:]')
        echo "Job ID: $PREV_JOB_ID"
    fi
done

echo ""
echo "============================================"
echo "All ${#FILE_STEMS[@]} jobs submitted as a chain."
echo "Jobs will run sequentially, each starting after the previous finishes."
echo "Monitor with: qstat -u jozefov_147"
echo "============================================"
