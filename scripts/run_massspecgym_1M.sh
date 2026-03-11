#!/bin/bash
#PBS -N canon_1M
#PBS -l select=1:ncpus=16:mem=16gb:scratch_local=10gb
#PBS -l walltime=4:00:00

set -e

PLZEN_HOME=/storage/plzen1/home/jozefov_147
ENV_PREFIX=$PLZEN_HOME/.conda/envs/chembl_pipeline
REPO_DIR=$PLZEN_HOME/projects/ChEMBL_Structure_Pipeline
DATA_DIR=$PLZEN_HOME/projects/msngym/data/candidates_generation/pubchem_canon_molecules
INPUT=$DATA_DIR/MassSpecGym_retrieval_molecules_1M.tsv
OUTPUT=$DATA_DIR/MassSpecGym_retrieval_molecules_1M_canonicalized.tsv

mkdir -p "$SCRATCHDIR/tmp"
export TMPDIR="$SCRATCHDIR/tmp"

module add mambaforge
export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
mamba activate "$ENV_PREFIX"

cp "$INPUT" "$SCRATCHDIR/input.tsv" || exit 2

echo "Starting canonicalization of 1M dataset at $(date)"

# Run python — tee stdout to log, keep stderr visible to PBS
python "$REPO_DIR/scripts/canonicalize_pubchem.py" \
    "$SCRATCHDIR/input.tsv" \
    "$SCRATCHDIR/output.tsv" \
    --workers 16 \
    --batch-size 500000 \
    2>&1 | tee "$SCRATCHDIR/canon_1M.log" || {
    # Copy log even on failure
    cp "$SCRATCHDIR/canon_1M.log" "$DATA_DIR/canon_1M.log" 2>/dev/null
    exit 3
}

cp "$SCRATCHDIR/output.tsv" "$OUTPUT" || exit 4
cp "$SCRATCHDIR/canon_1M.log" "$DATA_DIR/canon_1M.log" || exit 5
echo "Done at $(date)"
clean_scratch
