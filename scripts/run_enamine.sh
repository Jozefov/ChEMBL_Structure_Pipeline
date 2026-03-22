#!/bin/bash
# PBS job script for processing a single Enamine REAL file on MetaCentrum.
#
# Self-resubmitting: if the file isn't fully processed within the walltime
# budget, the job stages out partial results and resubmits itself with
# RESUME=1. This repeats until the file is done — no manual intervention.
#
# Usage:
#   # Submit for a specific file:
#   FILE_STEM="HAC_11_21" qsub scripts/run_enamine.sh
#
#   # Resume an interrupted job (manual):
#   FILE_STEM="HAC_11_21" RESUME=1 qsub scripts/run_enamine.sh
#
# File stems (short form → full name):
#   HAC_11_21  → 2025.02_Enamine_REAL_HAC_11_21_967M_CXSMILES
#   HAC_22_23  → 2025.02_Enamine_REAL_HAC_22_23_1.1B_CXSMILES
#   HAC_24     → 2025.02_Enamine_REAL_HAC_24_849M_CXSMILES
#   HAC_25     → 2025.02_Enamine_REAL_HAC_25_1.1B_CXSMILES
#   HAC_26     → 2025.02_Enamine_REAL_HAC_26_1.2B_CXSMILES
#   HAC_27     → 2025.02_Enamine_REAL_HAC_27_1.4B_CXSMILES
#   HAC_28     → 2025.02_Enamine_REAL_HAC_28_1.3B_CXSMILES
#   HAC_29_38_P1 → 2025.02_Enamine_REAL_HAC_29_38_2.5B_Part_1_CXSMILES
#   HAC_29_38_P2 → 2025.02_Enamine_REAL_HAC_29_38_2.5B_Part_2_CXSMILES

#PBS -N enamine_process
#PBS -l select=1:ncpus=18:mem=32gb:scratch_local=50gb
#PBS -l walltime=48:00:00
#PBS -j oe

# ============================================================================
# Configuration
# ============================================================================

FILE_STEM="${FILE_STEM:?ERROR: Set FILE_STEM before qsub (e.g. FILE_STEM=HAC_11_21)}"

ENAMINE_DIR="/storage/projects-du-praha/dreams/datasets/enamine_real"
PERSISTENT_OUT="${ENAMINE_DIR}/processed/${FILE_STEM}"
CONDA_ENV="/storage/plzen1/home/jozefov_147/.conda/envs/chembl_pipeline"
REPO_DIR="/storage/plzen1/home/jozefov_147/projects/ChEMBL_Structure_Pipeline"
N_WORKERS=16

# Time budget: stop processing 1 hour before walltime to allow stage-out
# and self-resubmit. PBS walltime is 48h = 172800s, budget = 47h = 169200s.
MAX_TIME=169200

# Map short stems to full filenames
declare -A FILE_MAP=(
    ["HAC_11_21"]="2025.02_Enamine_REAL_HAC_11_21_967M_CXSMILES"
    ["HAC_22_23"]="2025.02_Enamine_REAL_HAC_22_23_1.1B_CXSMILES"
    ["HAC_24"]="2025.02_Enamine_REAL_HAC_24_849M_CXSMILES"
    ["HAC_25"]="2025.02_Enamine_REAL_HAC_25_1.1B_CXSMILES"
    ["HAC_26"]="2025.02_Enamine_REAL_HAC_26_1.2B_CXSMILES"
    ["HAC_27"]="2025.02_Enamine_REAL_HAC_27_1.4B_CXSMILES"
    ["HAC_28"]="2025.02_Enamine_REAL_HAC_28_1.3B_CXSMILES"
    ["HAC_29_38_P1"]="2025.02_Enamine_REAL_HAC_29_38_2.5B_Part_1_CXSMILES"
    ["HAC_29_38_P2"]="2025.02_Enamine_REAL_HAC_29_38_2.5B_Part_2_CXSMILES"
)

# Resolve short stem to full filename if applicable
FULL_STEM="${FILE_MAP[$FILE_STEM]:-$FILE_STEM}"
INPUT_FILE="${ENAMINE_DIR}/${FULL_STEM}.cxsmiles.bz2"

# ============================================================================
# Cleanup stale PBS spool files from previous runs on this node
# ============================================================================
# PBS .OU/.ER files accumulate in /var/spool/pbs/undelivered/ and can fill
# the node's ~1GB local quota. Clean up any old files from previous jobs.
if [ -d /var/spool/pbs/undelivered ]; then
    find /var/spool/pbs/undelivered -user "$USER" -type f -delete 2>/dev/null || true
fi

# ============================================================================
# Environment setup
# ============================================================================

module add mambaforge
export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
export TMPDIR="$SCRATCHDIR/tmp"
mkdir -p "$TMPDIR"
mamba activate "$CONDA_ENV" || conda activate "$CONDA_ENV"

# Add repo to PYTHONPATH so scripts can import chembl_structure_pipeline
export PYTHONPATH="${REPO_DIR}:${PYTHONPATH}"

SCRATCH_OUT="$SCRATCHDIR/processed/${FILE_STEM}"
mkdir -p "$SCRATCH_OUT" "$PERSISTENT_OUT"

echo "============================================"
echo "Enamine REAL Processing Job"
echo "============================================"
echo "File:       $INPUT_FILE"
echo "Output:     $PERSISTENT_OUT"
echo "Scratch:    $SCRATCH_OUT"
echo "Workers:    $N_WORKERS"
echo "Max time:   ${MAX_TIME}s ($(( MAX_TIME / 3600 ))h)"
echo "Start time: $(date)"
echo "Node:       $(hostname)"
echo "============================================"

# Verify input exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# ============================================================================
# Resume support
# ============================================================================

RESUME_FLAG=""
if [ "${RESUME:-0}" = "1" ] && [ -f "$PERSISTENT_OUT/checkpoint.json" ]; then
    echo "Resuming from checkpoint..."
    cat "$PERSISTENT_OUT/checkpoint.json"
    echo ""
    # Copy partial outputs from persistent storage to scratch
    for f in molecules.tsv.gz errors.tsv.gz; do
        [ -f "$PERSISTENT_OUT/$f" ] && cp "$PERSISTENT_OUT/$f" "$SCRATCH_OUT/$f"
    done
    RESUME_FLAG="--resume"
fi

# ============================================================================
# Main processing
# ============================================================================

# Use pbzip2 for parallel decompression (2 threads)
# Remaining 16 CPUs for RDKit workers
# Note: lbzip2 not available on MetaCentrum, pbzip2 is
# Redirect all output to log file (not PBS stdout) to avoid filling node quota.
# PBS .OU files accumulate on the compute node's local disk and are quota-limited.
LOG_FILE="$SCRATCH_OUT/processing.log"

pbzip2 -dc -p2 "$INPUT_FILE" \
    | python3 "$REPO_DIR/scripts/enamine_process.py" \
        --output-dir "$PERSISTENT_OUT" \
        --scratch-dir "$SCRATCH_OUT" \
        --workers "$N_WORKERS" \
        --batch-size 100000 \
        --chunk-size 5000 \
        --chunk-timeout 3600 \
        --max-time "$MAX_TIME" \
        --skip-header \
        $RESUME_FLAG \
    > "$LOG_FILE" 2>&1

EXIT_CODE=${PIPESTATUS[1]}

# Show last 20 lines of log in PBS output (small enough for quota)
echo ""
echo "--- Last 20 lines of processing log ---"
tail -20 "$LOG_FILE"
echo "--- End ---"

# ============================================================================
# Stage out from scratch to persistent storage
# ============================================================================

echo ""
echo "Staging out results..."
for f in molecules.tsv.gz errors.tsv.gz processing.log; do
    if [ -f "$SCRATCH_OUT/$f" ]; then
        cp "$SCRATCH_OUT/$f" "$PERSISTENT_OUT/$f"
        echo "  Copied $f"
    fi
done

# ============================================================================
# Self-resubmit if not finished (exit code 42 = time budget exhausted)
# ============================================================================

if [ "$EXIT_CODE" -eq 42 ]; then
    echo ""
    echo "============================================"
    echo "Time budget exhausted — resubmitting with resume..."
    echo "============================================"
    NEXT_JOB=$(qsub -v "FILE_STEM=${FILE_STEM},RESUME=1" \
        -N "enamine_${FILE_STEM}" \
        "$REPO_DIR/scripts/run_enamine.sh")
    echo "Resubmitted as: $NEXT_JOB"
    echo "============================================"
    # Exit 0 so the chain dependency (afterany) considers this job "done"
    # and doesn't prematurely release the next file's job.
    # The resubmitted resume job runs independently of the chain.
    EXIT_CODE=0
fi

echo ""
echo "============================================"
echo "Job finished at $(date)"
echo "Exit code: $EXIT_CODE"
echo "Output: $PERSISTENT_OUT"
ls -lh "$PERSISTENT_OUT/"
echo "============================================"

# Cleanup scratch
rm -rf "$SCRATCHDIR"/*

exit $EXIT_CODE
