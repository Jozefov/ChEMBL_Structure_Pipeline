#!/bin/bash
# Quick smoke test for the Enamine processing pipeline on MetaCentrum.
#
# Streams the first 5000 lines from the smallest Enamine file,
# processes them, verifies output alignment, and cleans up.
#
# Run interactively on a frontend or inside qsub -I:
#   bash scripts/smoke_test_enamine.sh
#
# Or submit as a short job:
#   qsub -l select=1:ncpus=4:mem=8gb:scratch_local=5gb -l walltime=00:30:00 scripts/smoke_test_enamine.sh

#PBS -N enamine_smoke_test
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=5gb
#PBS -l walltime=00:30:00
#PBS -j oe

# No set -e: conda activate and pipe commands can return non-zero spuriously

# ============================================================================
# Configuration
# ============================================================================

ENAMINE_DIR="/storage/projects-du-praha/dreams/datasets/enamine_real"
FIRST_FILE="2025.02_Enamine_REAL_HAC_11_21_967M_CXSMILES.cxsmiles.bz2"
INPUT="$ENAMINE_DIR/$FIRST_FILE"
N_LINES=5000
N_WORKERS=2

CONDA_ENV="/storage/plzen1/home/jozefov_147/.conda/envs/chembl_pipeline"
REPO_DIR="/storage/plzen1/home/jozefov_147/projects/ChEMBL_Structure_Pipeline"

# Use scratch if available (PBS job), otherwise /tmp
OUTPUT_DIR="${SCRATCHDIR:-/tmp}/enamine_smoke_test"

# ============================================================================
# Setup
# ============================================================================

module add mambaforge 2>/dev/null
export CONDA_PKGS_DIRS="${SCRATCHDIR:-/tmp}/conda_pkgs"
mamba activate "$CONDA_ENV" 2>/dev/null || conda activate "$CONDA_ENV" 2>/dev/null

export PYTHONPATH="${REPO_DIR}:${PYTHONPATH:-}"

mkdir -p "$OUTPUT_DIR"

echo "============================================"
echo "Enamine Smoke Test"
echo "============================================"
echo "Input:    $INPUT"
echo "Lines:    $N_LINES"
echo "Workers:  $N_WORKERS"
echo "Output:   $OUTPUT_DIR"
echo "============================================"

# ============================================================================
# 1. Check input file exists
# ============================================================================

if [ ! -f "$INPUT" ]; then
    echo "FAIL: Input file not found: $INPUT"
    exit 1
fi
echo "OK: Input file exists"

# ============================================================================
# 2. Check first few lines (header detection)
# ============================================================================

echo ""
echo "--- First 3 lines of input ---"
bzcat "$INPUT" 2>/dev/null | head -3 || true
echo "--- End ---"
echo ""

FIRST_LINE=$(bzcat "$INPUT" 2>/dev/null | head -1 || true)
if echo "$FIRST_LINE" | grep -qi "smiles"; then
    echo "DETECTED: File has a header line. Using --skip-header."
    SKIP_HEADER="--skip-header"
else
    echo "DETECTED: No header line. Raw CXSMILES data."
    SKIP_HEADER=""
fi

# ============================================================================
# 3. Check lbzip2 availability
# ============================================================================

if command -v pbzip2 &>/dev/null; then
    DECOMPRESS="pbzip2 -dc -p2"
    echo "OK: pbzip2 available (parallel decompression)"
elif command -v lbzip2 &>/dev/null; then
    DECOMPRESS="lbzip2 -dc -n 2"
    echo "OK: lbzip2 available"
else
    DECOMPRESS="bzcat"
    echo "WARN: no parallel decompressor found, falling back to bzcat"
fi

# ============================================================================
# 4. Process N_LINES molecules
# ============================================================================

echo ""
echo "Extracting $N_LINES lines to temp file..."
SAMPLE_FILE="$OUTPUT_DIR/sample.tsv"
$DECOMPRESS "$INPUT" 2>/dev/null | head -n "$N_LINES" > "$SAMPLE_FILE" || true
ACTUAL_LINES=$(wc -l < "$SAMPLE_FILE")
echo "Extracted $ACTUAL_LINES lines"

echo ""
echo "Processing $ACTUAL_LINES molecules..."
cat "$SAMPLE_FILE" | \
    python3 "$REPO_DIR/scripts/enamine_process.py" \
        --output-dir "$OUTPUT_DIR" \
        --workers "$N_WORKERS" \
        --batch-size 2000 \
        --chunk-size 500 \
        --chunk-timeout 120 \
        $SKIP_HEADER
rm -f "$SAMPLE_FILE"

# ============================================================================
# 5. Verify outputs
# ============================================================================

echo ""
echo "============================================"
echo "Verification"
echo "============================================"

# Check files exist
for f in molecules.tsv.gz fingerprints.bin bitsums.bin errors.tsv.gz meta.json; do
    if [ -f "$OUTPUT_DIR/$f" ]; then
        SIZE=$(ls -lh "$OUTPUT_DIR/$f" | awk '{print $5}')
        echo "OK: $f ($SIZE)"
    else
        echo "FAIL: $f missing!"
    fi
done

# Check alignment
FPS_SIZE=$(stat -c%s "$OUTPUT_DIR/fingerprints.bin" 2>/dev/null || stat -f%z "$OUTPUT_DIR/fingerprints.bin")
SUMS_SIZE=$(stat -c%s "$OUTPUT_DIR/bitsums.bin" 2>/dev/null || stat -f%z "$OUTPUT_DIR/bitsums.bin")
FPS_ROWS=$((FPS_SIZE / 256))
SUMS_ROWS=$((SUMS_SIZE / 2))

echo ""
echo "Fingerprints: $FPS_SIZE bytes → $FPS_ROWS molecules"
echo "Bitsums:      $SUMS_SIZE bytes → $SUMS_ROWS molecules"

if [ "$FPS_ROWS" -eq "$SUMS_ROWS" ]; then
    echo "OK: Alignment check passed (fps rows == bitsums rows)"
else
    echo "FAIL: Alignment mismatch! fps=$FPS_ROWS, sums=$SUMS_ROWS"
fi

# Show meta.json
echo ""
echo "--- meta.json ---"
cat "$OUTPUT_DIR/meta.json"
echo ""
echo "--- End ---"

# Show first few output lines
echo ""
echo "--- First 3 lines of molecules.tsv.gz ---"
zcat "$OUTPUT_DIR/molecules.tsv.gz" | head -3
echo "--- End ---"

# ============================================================================
# 6. Summary
# ============================================================================

echo ""
echo "============================================"
N_OK=$(python3 -c "import json; print(json.load(open('$OUTPUT_DIR/meta.json'))['n_molecules'])")
N_ERR=$(python3 -c "import json; print(json.load(open('$OUTPUT_DIR/meta.json'))['n_errors'])")
TOTAL=$((N_OK + N_ERR))
echo "Smoke test complete: $N_OK ok, $N_ERR errors out of $TOTAL"
if [ "$N_OK" -gt 0 ]; then
    echo "SUCCESS: Pipeline is working!"
else
    echo "FAIL: No molecules processed successfully"
    exit 1
fi
echo "============================================"

# Cleanup
rm -rf "$OUTPUT_DIR"
