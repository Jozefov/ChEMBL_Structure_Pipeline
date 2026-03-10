[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# ChEMBL Structure Pipeline (Extended)

A fork of the [ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline) that adds **salt stripping**, **tautomer canonicalization**, and **isotope preservation** in a single function call.

Given a SMILES string, the pipeline normalizes functional groups, strips salts and solvents, neutralizes charges, and canonicalizes tautomers — so that different representations of the same molecule always produce the same SMILES. Stereochemistry (R/S, E/Z) and isotope labels (²H, ¹³C) are preserved.

## Installation

```bash
pip install git+https://github.com/Jozefov/ChEMBL_Structure_Pipeline.git
```

## Usage

```python
from chembl_structure_pipeline.standardizer import standardize_and_canonicalize_smiles

standardize_and_canonicalize_smiles("CC(O)=N")           # tautomer → 'CC(N)=O'
standardize_and_canonicalize_smiles("CC(=O)N")           # same       'CC(N)=O'
standardize_and_canonicalize_smiles("C[NH3+].[Cl-]")     # salt       'CN'
standardize_and_canonicalize_smiles("N[C@@H](C)C(=O)O")  # stereo     'C[C@H](N)C(=O)O'
standardize_and_canonicalize_smiles("N[C@H](C)C(=O)O")   # distinct   'C[C@@H](N)C(=O)O'
standardize_and_canonicalize_smiles("[2H]c1ccccc1")       # isotope    '[2H]c1ccccc1'
```

For Mol objects use `standardize_and_canonicalize_mol(mol)` instead.

All original ChEMBL functions (`standardize_mol`, `get_parent_mol`, `check_molblock`, etc.) remain available and unchanged.

## References

Bento, A.P., Hersey, A., Félix, E. et al. An open source chemical structure curation pipeline using RDKit. *J Cheminform* 12, 51 (2020). https://doi.org/10.1186/s13321-020-00456-1
