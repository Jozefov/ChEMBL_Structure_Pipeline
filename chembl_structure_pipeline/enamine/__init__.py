"""Enamine REAL dataset processing subpackage.

Provides streaming canonicalization, property computation, and fingerprint
generation for the ~11.4 billion molecule Enamine REAL library.

Output per input file (aligned by row index):
  molecules.tsv.gz   — canonical_smiles, inchikey, inchikey14, formula, mass
  fingerprints.bin   — N × 256 uint8 packed Morgan fingerprints
  bitsums.bin        — N × 2 uint16 bit sums (for fast Tanimoto)
  errors.tsv.gz      — failed molecules with original SMILES
  meta.json          — processing stats and FP parameters
"""

# Fingerprint parameters (frozen for reproducibility across all runs)
MORGAN_RADIUS = 2
MORGAN_NBITS = 2048
PACKED_WIDTH = (MORGAN_NBITS + 7) // 8  # 256 bytes per fingerprint

# Output TSV column names
TSV_COLUMNS = (
    "original_smiles",
    "canonical_smiles",
    "inchikey",
    "inchikey14",
    "molecular_formula",
    "monoisotopic_mass",
)

ERROR_COLUMNS = (
    "original_smiles",
    "error_reason",
)
