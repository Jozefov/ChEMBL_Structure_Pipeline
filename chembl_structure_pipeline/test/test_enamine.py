"""Tests for the enamine subpackage."""

import gzip
import json
import os
import struct
import tempfile

import numpy as np
import pytest

from chembl_structure_pipeline.enamine import (
    MORGAN_RADIUS, MORGAN_NBITS, PACKED_WIDTH, TSV_COLUMNS,
)
from chembl_structure_pipeline.enamine.worker import (
    strip_cxsmiles,
    process_single_molecule,
    worker_process_chunk,
)


# ---------------------------------------------------------------------------
# strip_cxsmiles
# ---------------------------------------------------------------------------

class TestStripCxsmiles:
    def test_plain_smiles(self):
        assert strip_cxsmiles("CCO") == "CCO"

    def test_strip_extension(self):
        assert strip_cxsmiles("CCO |c:0,1|") == "CCO"

    def test_complex_extension(self):
        smi = "C1=CC=CC=C1 |SgD:0,1:amine:notes|"
        assert strip_cxsmiles(smi) == "C1=CC=CC=C1"

    def test_no_space_before_pipe(self):
        # If there's no space before pipe, it's not CXSMILES extension
        assert strip_cxsmiles("CC|stuff") == "CC|stuff"

    def test_empty_string(self):
        assert strip_cxsmiles("") == ""

    def test_only_extension(self):
        assert strip_cxsmiles(" |ext|") == ""

    def test_multiple_pipes(self):
        smi = "C(=O)O |c:0| |extra|"
        # Should strip at first " |"
        assert strip_cxsmiles(smi) == "C(=O)O"


# ---------------------------------------------------------------------------
# process_single_molecule
# ---------------------------------------------------------------------------

class TestProcessSingleMolecule:
    def test_ethanol(self):
        result = process_single_molecule("CCO")
        assert result is not None
        canon_smi, inchikey, inchikey14, formula, mass, packed_fp = result

        assert isinstance(canon_smi, str)
        assert len(canon_smi) > 0
        assert len(inchikey) == 27  # full InChIKey
        assert len(inchikey14) == 14  # connectivity layer
        assert inchikey.startswith(inchikey14)
        assert formula == "C2H6O"
        assert abs(mass - 46.041865) < 0.001
        assert isinstance(packed_fp, bytes)
        assert len(packed_fp) == PACKED_WIDTH  # 256 bytes

    def test_aspirin(self):
        result = process_single_molecule("CC(=O)Oc1ccccc1C(=O)O")
        assert result is not None
        _, _, _, formula, mass, _ = result
        assert formula == "C9H8O4"
        assert abs(mass - 180.042259) < 0.001

    def test_cxsmiles_input(self):
        result = process_single_molecule("CCO |c:0,1|")
        assert result is not None
        canon_smi, _, _, _, _, _ = result
        assert "|" not in canon_smi  # CXSMILES extension stripped

    def test_invalid_smiles(self):
        result = process_single_molecule("not_a_smiles_XYZ123")
        assert result is None

    def test_empty_smiles(self):
        result = process_single_molecule("")
        assert result is None

    def test_fingerprint_not_all_zeros(self):
        result = process_single_molecule("c1ccccc1")  # benzene
        assert result is not None
        packed_fp = result[5]
        # Benzene should have some bits set
        bits = np.unpackbits(np.frombuffer(packed_fp, dtype=np.uint8))
        assert bits.sum() > 0

    def test_different_molecules_different_fps(self):
        r1 = process_single_molecule("CCO")
        r2 = process_single_molecule("c1ccccc1")
        assert r1 is not None and r2 is not None
        assert r1[5] != r2[5]  # different fingerprints


# ---------------------------------------------------------------------------
# worker_process_chunk
# ---------------------------------------------------------------------------

class TestWorkerProcessChunk:
    def test_basic_chunk(self):
        chunk = [(0, "CCO"), (1, "c1ccccc1"), (2, "CC(=O)O")]
        results = worker_process_chunk(chunk)
        assert len(results) == 3
        for idx, result in results:
            assert result is not None
            assert len(result) == 6

    def test_mixed_valid_invalid(self):
        chunk = [(0, "CCO"), (1, "INVALID_SMILES"), (2, "c1ccccc1")]
        results = worker_process_chunk(chunk)
        assert len(results) == 3
        assert results[0][1] is not None  # CCO ok
        assert results[1][1] is None      # invalid
        assert results[2][1] is not None  # benzene ok

    def test_empty_chunk(self):
        results = worker_process_chunk([])
        assert results == []

    def test_indices_preserved(self):
        chunk = [(42, "CCO"), (99, "c1ccccc1")]
        results = worker_process_chunk(chunk)
        indices = [idx for idx, _ in results]
        assert indices == [42, 99]


# ---------------------------------------------------------------------------
# Integration: process_enamine_stream (small scale)
# ---------------------------------------------------------------------------

class TestProcessorIntegration:
    def test_small_stream(self):
        """Process a small list of SMILES and verify aligned outputs."""
        from chembl_structure_pipeline.enamine.processor import process_enamine_stream

        smiles = [
            "CCO",
            "c1ccccc1",
            "INVALID_XYZ",
            "CC(=O)O",
            "CC(=O)Oc1ccccc1C(=O)O",
        ]
        lines = [s + "\n" for s in smiles]

        with tempfile.TemporaryDirectory() as tmpdir:
            stats = process_enamine_stream(
                input_lines=iter(lines),
                output_dir=tmpdir,
                n_workers=2,
                batch_size=10,
                chunk_size=3,
                chunk_timeout=60,
            )

            # Check stats
            assert stats["lines_processed"] == 5
            assert stats["n_molecules"] == 4  # 1 invalid
            assert stats["n_errors"] == 1

            # Check alignment: TSV rows == fps bytes / 256 == sums bytes / 2
            fps_path = os.path.join(tmpdir, "fingerprints.bin")
            sums_path = os.path.join(tmpdir, "bitsums.bin")
            tsv_path = os.path.join(tmpdir, "molecules.tsv.gz")

            fps_size = os.path.getsize(fps_path)
            sums_size = os.path.getsize(sums_path)

            assert fps_size == 4 * PACKED_WIDTH  # 4 valid molecules
            assert sums_size == 4 * 2  # 4 × uint16

            # Count TSV rows (excluding header)
            with gzip.open(tsv_path, "rt") as f:
                lines_out = f.readlines()
            assert lines_out[0].startswith("original_smiles")  # header
            assert len(lines_out) == 5  # header + 4 data rows
            # Verify header has all 6 columns
            header_cols = lines_out[0].rstrip("\n").split("\t")
            assert header_cols == ["original_smiles", "canonical_smiles", "inchikey",
                                   "inchikey14", "molecular_formula", "monoisotopic_mass"]
            # Verify data rows have original_smiles as first column
            for data_line in lines_out[1:]:
                fields = data_line.rstrip("\n").split("\t")
                assert len(fields) == 6

            # Check bitsums are valid uint16
            with open(sums_path, "rb") as f:
                sums_data = np.frombuffer(f.read(), dtype=np.uint16)
            assert len(sums_data) == 4
            assert all(s > 0 for s in sums_data)  # all valid molecules have some bits set

            # Check meta.json exists
            meta_path = os.path.join(tmpdir, "meta.json")
            assert os.path.exists(meta_path)
            with open(meta_path) as f:
                meta = json.load(f)
            assert meta["n_molecules"] == 4

            # Check errors file
            errors_path = os.path.join(tmpdir, "errors.tsv.gz")
            with gzip.open(errors_path, "rt") as f:
                err_lines = f.readlines()
            assert len(err_lines) == 2  # header + 1 error

    def test_resume(self):
        """Process in two halves and verify resume produces correct output."""
        from chembl_structure_pipeline.enamine.processor import process_enamine_stream

        smiles = ["CCO", "c1ccccc1", "CC(=O)O", "CCCC"]
        lines = [s + "\n" for s in smiles]

        with tempfile.TemporaryDirectory() as tmpdir:
            # Process first 2
            stats1 = process_enamine_stream(
                input_lines=iter(lines[:2]),
                output_dir=tmpdir,
                n_workers=1,
                batch_size=10,
                chunk_size=10,
                chunk_timeout=60,
            )
            assert stats1["n_molecules"] == 2

            # Resume with all 4 (first 2 will be skipped via checkpoint)
            stats2 = process_enamine_stream(
                input_lines=iter(lines),
                output_dir=tmpdir,
                n_workers=1,
                batch_size=10,
                chunk_size=10,
                chunk_timeout=60,
                resume_from=0,  # auto-loads from checkpoint
            )
            assert stats2["n_molecules"] == 4

            # Verify final alignment
            fps_size = os.path.getsize(os.path.join(tmpdir, "fingerprints.bin"))
            assert fps_size == 4 * PACKED_WIDTH


# ---------------------------------------------------------------------------
# Merge utility
# ---------------------------------------------------------------------------

class TestMerge:
    def test_merge_two_dirs(self):
        """Create two small per-file outputs and merge them."""
        from chembl_structure_pipeline.enamine.processor import process_enamine_stream
        from chembl_structure_pipeline.enamine.merge import merge_enamine_outputs

        with tempfile.TemporaryDirectory() as tmpdir:
            dir1 = os.path.join(tmpdir, "HAC_A")
            dir2 = os.path.join(tmpdir, "HAC_B")
            merge_dir = os.path.join(tmpdir, "merged")

            # Process two small datasets
            process_enamine_stream(
                input_lines=iter(["CCO\n", "c1ccccc1\n"]),
                output_dir=dir1,
                n_workers=1, batch_size=10, chunk_size=10, chunk_timeout=60,
            )
            process_enamine_stream(
                input_lines=iter(["CC(=O)O\n", "CCCC\n", "CCCCC\n"]),
                output_dir=dir2,
                n_workers=1, batch_size=10, chunk_size=10, chunk_timeout=60,
            )

            # Merge
            result = merge_enamine_outputs([dir1, dir2], merge_dir)

            assert result["total_molecules"] == 5
            assert result["n_files"] == 2

            # Check file offsets
            offsets_path = os.path.join(merge_dir, "file_offsets.json")
            with open(offsets_path) as f:
                offsets = json.load(f)
            assert "HAC_A" in offsets
            assert "HAC_B" in offsets
            assert offsets["HAC_A"]["start_row"] == 0
            assert offsets["HAC_A"]["end_row"] == 2
            assert offsets["HAC_B"]["start_row"] == 2
            assert offsets["HAC_B"]["end_row"] == 5

            # Check InChIKey14 index exists
            index_path = os.path.join(merge_dir, "inchikey14_index.tsv.gz")
            assert os.path.exists(index_path)
            with gzip.open(index_path, "rt") as f:
                index_lines = f.readlines()
            # header + 5 entries
            assert len(index_lines) == 6
