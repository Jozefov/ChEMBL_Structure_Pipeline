#!/bin/bash
#PBS -N pubchem_canon
#PBS -l select=1:ncpus=48:mem=96gb:scratch_local=40gb
#PBS -l walltime=48:00:00
#PBS -o /storage/plzen1/home/jozefov_147/projects/msngym/data/candidates_generation/pubchem_canon_molecules/pubchem_canon.stdout
#PBS -e /storage/plzen1/home/jozefov_147/projects/msngym/data/candidates_generation/pubchem_canon_molecules/pubchem_canon.stderr

PLZEN_HOME=/storage/plzen1/home/jozefov_147
ENV_PREFIX=$PLZEN_HOME/.conda/envs/chembl_pipeline
REPO_DIR=$PLZEN_HOME/projects/ChEMBL_Structure_Pipeline
INPUT=$PLZEN_HOME/projects/msngym/data/candidates_generation/pubchem_canon_molecules/pubchem.tsv
OUTPUT_DIR=$PLZEN_HOME/projects/msngym/data/candidates_generation/pubchem_canon_molecules
OUTPUT=$OUTPUT_DIR/pubchem_canonicalized.tsv

mkdir -p "$SCRATCHDIR/tmp"
export TMPDIR="$SCRATCHDIR/tmp"

module add mambaforge
export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
mamba activate "$ENV_PREFIX"

# Copy input to scratch for faster I/O
cp "$INPUT" "$SCRATCHDIR/pubchem.tsv" || exit 2

# Copy partial output and progress if resuming from a previous run
if [ -f "$OUTPUT" ]; then
    cp "$OUTPUT" "$SCRATCHDIR/pubchem_canonicalized.tsv"
fi
if [ -f "$OUTPUT.progress.json" ]; then
    cp "$OUTPUT.progress.json" "$SCRATCHDIR/pubchem_canonicalized.tsv.progress.json"
    RESUME_FLAG="--resume"
else
    RESUME_FLAG=""
fi

echo "Starting canonicalization at $(date)"
echo "Input: $INPUT"
echo "Workers: 48, Batch size: 500000"

# Run canonicalization — write to BOTH scratch (fast) and persistent storage (safe)
# The --checkpoint-dir flag saves progress to persistent storage after every batch
python "$REPO_DIR/scripts/canonicalize_pubchem.py" \
    "$SCRATCHDIR/pubchem.tsv" \
    "$SCRATCHDIR/pubchem_canonicalized.tsv" \
    --workers 48 \
    --batch-size 500000 \
    --checkpoint-dir "$OUTPUT_DIR" \
    $RESUME_FLAG \
    2>&1 | tee "$SCRATCHDIR/pubchem_canon.log"

EXIT_CODE=${PIPESTATUS[0]}

echo "Python exit code: $EXIT_CODE"
echo "Finished at $(date)"

# Always copy whatever we have back to persistent storage
cp "$SCRATCHDIR/pubchem_canonicalized.tsv" "$OUTPUT" 2>/dev/null
cp "$SCRATCHDIR/pubchem_canonicalized.tsv.progress.json" "$OUTPUT.progress.json" 2>/dev/null
cp "$SCRATCHDIR/pubchem_canon.log" "$OUTPUT_DIR/pubchem_canon.log" 2>/dev/null

echo "Output: $OUTPUT"
clean_scratch

exit $EXIT_CODE
