"""Worker functions for Enamine REAL molecule processing.

Each function runs in a child process via ProcessPoolExecutor.
The module-level imports of standardizer/RDKit happen inside the
worker function body so that globals are initialised fresh in every
spawned process (required for the ``spawn`` start method on macOS).

A single worker call processes one molecule through the full pipeline:
  CXSMILES → strip → canonicalize → InChIKey → formula → mass → fingerprint
"""

import numpy as np

from . import MORGAN_RADIUS, MORGAN_NBITS, PACKED_WIDTH


def strip_cxsmiles(smi: str) -> str:
    """Remove CXSMILES extended notation from a SMILES string.

    CXSMILES appends extended info after `` |`` (space-pipe), e.g.::

        CCO |c:0,1,SgD:2:amine:...|

    This function returns only the SMILES part before the extension.
    """
    # The canonical separator is " |" (space then pipe)
    idx = smi.find(" |")
    if idx != -1:
        return smi[:idx]
    return smi


def process_single_molecule(smi: str):
    """Process one SMILES through the full Enamine pipeline.

    Args:
        smi: Raw SMILES (may contain CXSMILES extensions).

    Returns:
        Tuple of (canonical_smiles, inchikey, inchikey14, formula, mass,
        packed_fp_bytes) on success, or None on any failure.
    """
    from chembl_structure_pipeline import standardizer
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs, Descriptors, rdMolDescriptors
    from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

    # 1. Strip CXSMILES extended notation
    clean_smi = strip_cxsmiles(smi)
    if not clean_smi:
        return None

    # 2. Canonicalize via ChEMBL pipeline
    canon_smi = None
    try:
        canon_smi = standardizer.standardize_and_canonicalize_smiles(clean_smi)
    except Exception:
        canon_smi = None

    # 3. Parse mol from canonical (fallback to original)
    use_smi = canon_smi if canon_smi is not None else clean_smi
    try:
        mol = Chem.MolFromSmiles(use_smi)
    except Exception:
        mol = None
    if mol is None:
        return None

    # 4. InChIKey (full 27-char + first 14-char connectivity layer)
    inchikey = ""
    inchikey14 = ""
    try:
        inchi = MolToInchi(mol)
        if inchi is not None:
            ik = InchiToInchiKey(inchi)
            if ik:
                inchikey = ik
                inchikey14 = ik[:14]
    except Exception:
        pass

    # 5. Molecular formula
    formula = ""
    try:
        formula = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        pass

    # 6. Monoisotopic mass (ExactMolWt)
    mass = 0.0
    try:
        mass = round(Descriptors.ExactMolWt(mol), 6)
    except Exception:
        pass

    # 7. Morgan fingerprint → packed binary (256 bytes)
    try:
        generator = AllChem.GetMorganGenerator(
            radius=MORGAN_RADIUS, fpSize=MORGAN_NBITS,
        )
        fp = generator.GetFingerprint(mol)
        bits = np.zeros((MORGAN_NBITS,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, bits)
        packed = np.packbits(bits).tobytes()  # 256 bytes
    except Exception:
        return None

    return (
        canon_smi if canon_smi is not None else use_smi,
        inchikey,
        inchikey14,
        formula,
        mass,
        packed,  # bytes, len=256
    )


def worker_process_chunk(chunk):
    """Process a list of (index, smiles) tuples in a worker process.

    This is the function submitted to ProcessPoolExecutor. It imports
    all heavy modules inside the function body so that each worker
    process gets fresh module globals.

    Args:
        chunk: List of (int_index, smiles_string) tuples.

    Returns:
        List of (int_index, result_or_None) tuples, where result is
        the 6-tuple from process_single_molecule or None on failure.
    """
    results = []
    for idx, smi in chunk:
        try:
            result = process_single_molecule(smi)
        except Exception:
            result = None
        results.append((idx, result))
    return results
