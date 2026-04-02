#!/bin/bash
#PBS -N canon_4M
#PBS -l select=1:ncpus=16:mem=16gb:scratch_local=10gb
#PBS -l walltime=8:00:00
#PBS -o /storage/plzen1/home/jozefov_147/projects/msngym/data/candidates_generation/pubchem_canon_molecules/canon_4M.stdout
#PBS -e /storage/plzen1/home/jozefov_147/projects/msngym/data/candidates_generation/pubchem_canon_molecules/canon_4M.stderr

PLZEN_HOME=/storage/plzen1/home/jozefov_147
ENV_PREFIX=$PLZEN_HOME/.conda/envs/chembl_pipeline
REPO_DIR=$PLZEN_HOME/projects/ChEMBL_Structure_Pipeline
DATA_DIR=$PLZEN_HOME/projects/msngym/data/candidates_generation/pubchem_canon_molecules
INPUT=$DATA_DIR/MassSpecGym_retrieval_molecules_4M.tsv
OUTPUT=$DATA_DIR/MassSpecGym_retrieval_molecules_4M_canonicalized.tsv

mkdir -p "$SCRATCHDIR/tmp"
export TMPDIR="$SCRATCHDIR/tmp"

module add mambaforge
export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
mamba activate "$ENV_PREFIX"

# Cleanup stale PBS spool files from previous runs on this node
if [ -d /var/spool/pbs/undelivered ]; then
    find /var/spool/pbs/undelivered -user "$USER" -type f -delete 2>/dev/null || true
fi

cp "$INPUT" "$SCRATCHDIR/input.tsv" || exit 2

# Copy partial output if resuming
if [ -f "$OUTPUT" ]; then
    cp "$OUTPUT" "$SCRATCHDIR/output.tsv"
fi
if [ -f "$OUTPUT.progress.json" ]; then
    cp "$OUTPUT.progress.json" "$SCRATCHDIR/output.tsv.progress.json"
    RESUME_FLAG="--resume"
else
    RESUME_FLAG=""
fi

echo "Starting canonicalization of 4M dataset at $(date)"

python "$REPO_DIR/scripts/canonicalize_pubchem.py" \
    "$SCRATCHDIR/input.tsv" \
    "$SCRATCHDIR/output.tsv" \
    --workers 16 \
    --batch-size 500000 \
    --checkpoint-dir "$DATA_DIR" \
    --checkpoint-name "MassSpecGym_retrieval_molecules_4M_canonicalized.tsv" \
    $RESUME_FLAG \
    > "$SCRATCHDIR/canon_4M.log" 2>&1

EXIT_CODE=$?

# Always copy results back
cp "$SCRATCHDIR/output.tsv" "$OUTPUT" 2>/dev/null
cp "$SCRATCHDIR/output.tsv.progress.json" "$OUTPUT.progress.json" 2>/dev/null
cp "$SCRATCHDIR/canon_4M.log" "$DATA_DIR/canon_4M.log" 2>/dev/null

echo "Done at $(date)"
echo "=== Last 20 lines of processing log ==="
tail -20 "$SCRATCHDIR/canon_4M.log" 2>/dev/null || true
clean_scratch
exit $EXIT_CODE
