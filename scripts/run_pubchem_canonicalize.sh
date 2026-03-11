#!/bin/bash
#PBS -N pubchem_canon
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=20gb
#PBS -l walltime=24:00:00

set -e

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

echo "Starting canonicalization at $(date)"
echo "Input: $INPUT"
echo "Workers: 16"

python "$REPO_DIR/scripts/canonicalize_pubchem.py" \
    "$SCRATCHDIR/pubchem.tsv" \
    "$SCRATCHDIR/pubchem_canonicalized.tsv" \
    --workers 16 \
    --checkpoint \
    > "$SCRATCHDIR/pubchem_canon.log" 2>&1 || exit 3

echo "Finished at $(date)"

# Copy results back to persistent storage
cp "$SCRATCHDIR/pubchem_canonicalized.tsv" "$OUTPUT" || exit 4
cp "$SCRATCHDIR/pubchem_canon.log" "$OUTPUT_DIR/pubchem_canon.log" || exit 5

echo "Output: $OUTPUT"
clean_scratch
