"""Enamine REAL dataset processing subpackage.

Provides streaming canonicalization, property computation, and fingerprint
generation for the ~11.4 billion molecule Enamine REAL library.

Output per input file (aligned by row index):
  molecules.tsv.gz   — original_smiles, canonical_smiles, inchikey, inchikey14,
                       formula, mass + passthrough columns from input
  fingerprints.bin   — N × 256 uint8 packed Morgan fingerprints
  bitsums.bin        — N × 2 uint16 bit sums (for fast Tanimoto)
  errors.tsv.gz      — failed molecules with original SMILES
  meta.json          — processing stats and FP parameters
"""

# Fingerprint parameters (frozen for reproducibility across all runs)
MORGAN_RADIUS = 2
MORGAN_NBITS = 2048
PACKED_WIDTH = (MORGAN_NBITS + 7) // 8  # 256 bytes per fingerprint

# Core output TSV column names (always present, computed from canonical mol)
TSV_COLUMNS_CORE = (
    "original_smiles",
    "canonical_smiles",
    "inchikey",
    "inchikey14",
    "molecular_formula",
    "monoisotopic_mass",
)

# Input columns to drop (we recompute these from canonical mol)
# Matched case-insensitively against input header
INPUT_COLUMNS_DROP = {"smiles", "mw", "inchikey"}

ERROR_COLUMNS = (
    "original_smiles",
    "error_reason",
)
