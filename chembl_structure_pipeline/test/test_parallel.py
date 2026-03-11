"""Tests for parallel batch processing and checkpoint/recovery."""
import os
import tempfile

import pytest

from chembl_structure_pipeline.parallel import (
    batch_standardize_smiles,
    batch_standardize_molblocks,
    batch_check_molblocks,
    _load_checkpoint,
    _append_checkpoint,
)
from chembl_structure_pipeline.standardizer import (
    standardize_and_canonicalize_smiles,
    standardize_molblock,
)


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

# Mix of normal molecules, stereochemistry, isotopes, salts, and edge cases
TEST_SMILES = [
    "CC(=O)O",                          # acetic acid
    "CC(=O)N",                          # acetamide
    "CC(O)=N",                          # imidic acid (tautomer of acetamide)
    "N[C@@H](C)C(=O)O",                # L-alanine (R/S stereo)
    "N[C@H](C)C(=O)O",                 # D-alanine (opposite stereo)
    "[NH3+]CC([O-])=O",                 # glycine zwitterion
    "c1ccccc1",                         # benzene
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O.OC(=O)CC(O)(CC(=O)O)C(=O)O",  # caffeine citrate
    "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]",  # deuterated benzene
    r"C(\F)=C/C",                       # E/Z stereo
]

INVALID_SMILES = [
    "not_a_smiles",
    "C(C(C",
    "",
]


# ---------------------------------------------------------------------------
# Correctness: parallel == serial
# ---------------------------------------------------------------------------

class TestBatchStandardizeSmiles:
    def test_parallel_matches_serial(self):
        serial = [standardize_and_canonicalize_smiles(s) for s in TEST_SMILES]
        parallel = batch_standardize_smiles(TEST_SMILES, n_workers=2)
        assert parallel == serial

    def test_single_worker_matches_serial(self):
        serial = [standardize_and_canonicalize_smiles(s) for s in TEST_SMILES]
        parallel = batch_standardize_smiles(TEST_SMILES, n_workers=1)
        assert parallel == serial

    def test_empty_input(self):
        assert batch_standardize_smiles([]) == []

    def test_single_item(self):
        result = batch_standardize_smiles(["CCO"], n_workers=2)
        assert len(result) == 1
        assert result[0] == standardize_and_canonicalize_smiles("CCO")

    def test_error_handling(self):
        mixed = TEST_SMILES + INVALID_SMILES
        results = batch_standardize_smiles(mixed, n_workers=2)
        assert len(results) == len(mixed)
        # Valid SMILES should succeed
        for i in range(len(TEST_SMILES)):
            assert results[i] is not None, f"Expected valid result at index {i}"
        # Invalid SMILES should be None
        for i in range(len(TEST_SMILES), len(mixed)):
            assert results[i] is None, f"Expected None at index {i}"


# ---------------------------------------------------------------------------
# Stereochemistry preservation
# ---------------------------------------------------------------------------

class TestStereoPreservation:
    def test_r_s_stereo_preserved(self):
        results = batch_standardize_smiles(
            ["N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"], n_workers=2,
        )
        # L-alanine and D-alanine must produce different canonical SMILES
        assert results[0] is not None
        assert results[1] is not None
        assert results[0] != results[1], "R/S stereochemistry was lost"

    def test_e_z_stereo_preserved(self):
        results = batch_standardize_smiles(
            [r"C(\F)=C/C", r"C(/F)=C/C"], n_workers=2,
        )
        assert results[0] is not None
        assert results[1] is not None

    def test_isotopes_preserved(self):
        result = batch_standardize_smiles(
            ["[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]"], n_workers=1,
        )
        assert result[0] is not None
        assert "[2H]" in result[0], "Isotope labels were lost"


# ---------------------------------------------------------------------------
# Tautomer canonicalization
# ---------------------------------------------------------------------------

class TestTautomerCanonicalization:
    def test_tautomers_converge(self):
        results = batch_standardize_smiles(
            ["CC(=O)N", "CC(O)=N"], n_workers=2,
        )
        assert results[0] == results[1], (
            f"Tautomers did not converge: {results[0]} vs {results[1]}"
        )


# ---------------------------------------------------------------------------
# Checkpoint / recovery
# ---------------------------------------------------------------------------

class TestCheckpoint:
    def test_checkpoint_write_and_load(self):
        with tempfile.NamedTemporaryFile(suffix=".jsonl", delete=False) as f:
            path = f.name
        try:
            # Write some results
            _append_checkpoint(path, [(0, "CCO"), (1, None), (2, "CC")])
            completed = _load_checkpoint(path)
            assert completed == {0: "CCO", 1: None, 2: "CC"}
        finally:
            os.unlink(path)

    def test_checkpoint_resume(self):
        with tempfile.NamedTemporaryFile(suffix=".jsonl", delete=False) as f:
            path = f.name
        try:
            # First run: process only first 3
            smiles = TEST_SMILES[:5]
            results1 = batch_standardize_smiles(
                smiles, n_workers=1, checkpoint_path=path,
            )
            assert len(results1) == 5

            # Simulate partial checkpoint: keep only first 3 results
            completed = _load_checkpoint(path)
            assert len(completed) == 5

            # Truncate checkpoint to first 3 entries
            os.unlink(path)
            partial = [(i, results1[i]) for i in range(3)]
            _append_checkpoint(path, partial)

            # Resume should process remaining 2
            results2 = batch_standardize_smiles(
                smiles, n_workers=1, checkpoint_path=path,
            )
            assert results2 == results1, "Resume produced different results"
        finally:
            if os.path.exists(path):
                os.unlink(path)

    def test_no_checkpoint(self):
        # Should work fine without checkpoint
        results = batch_standardize_smiles(TEST_SMILES, n_workers=2)
        assert len(results) == len(TEST_SMILES)


# ---------------------------------------------------------------------------
# Salt stripping in parallel
# ---------------------------------------------------------------------------

class TestSaltStripping:
    def test_caffeine_citrate_stripped(self):
        results = batch_standardize_smiles(
            ["Cn1c(=O)c2c(ncn2C)n(C)c1=O.OC(=O)CC(O)(CC(=O)O)C(=O)O"],
            n_workers=1,
        )
        serial = standardize_and_canonicalize_smiles(
            "Cn1c(=O)c2c(ncn2C)n(C)c1=O.OC(=O)CC(O)(CC(=O)O)C(=O)O"
        )
        assert results[0] == serial
