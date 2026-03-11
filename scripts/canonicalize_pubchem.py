#!/usr/bin/env python
"""Canonicalize PubChem SMILES using ChEMBL Structure Pipeline.

Reads a TSV file with columns: cid, smiles, formula, mass
Standardizes and canonicalizes SMILES in parallel, preserving all columns.
Molecules that fail standardization keep their original SMILES.

Usage:
    python canonicalize_pubchem.py input.tsv output.tsv --workers 16
    python canonicalize_pubchem.py input.tsv output.tsv --checkpoint
    python canonicalize_pubchem.py input.tsv output.tsv --resume
"""
import argparse
import csv
import os
import sys
import time

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

from chembl_structure_pipeline.parallel import batch_standardize_smiles


def main():
    parser = argparse.ArgumentParser(
        description="Canonicalize PubChem SMILES with ChEMBL Structure Pipeline",
    )
    parser.add_argument("input", help="Input TSV (cid, smiles, formula, mass)")
    parser.add_argument("output", help="Output TSV with canonicalized SMILES")
    parser.add_argument(
        "--workers", "-w", type=int, default=None,
        help="Number of worker processes (default: all CPUs)",
    )
    parser.add_argument(
        "--checkpoint", action="store_true",
        help="Enable checkpointing for resume support",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from existing checkpoint file",
    )
    parser.add_argument(
        "--checkpoint-file", default=None,
        help="Custom checkpoint file path (default: <output>.checkpoint.jsonl)",
    )
    args = parser.parse_args()

    checkpoint_path = None
    if args.checkpoint or args.resume:
        checkpoint_path = args.checkpoint_file or (args.output + ".checkpoint.jsonl")

    if args.resume and checkpoint_path and not os.path.exists(checkpoint_path):
        print(f"Error: checkpoint file {checkpoint_path} not found", file=sys.stderr)
        sys.exit(1)

    # Read input TSV
    print(f"Reading {args.input}...")
    rows = []
    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    smiles_list = [row["smiles"] for row in rows]
    print(f"Read {len(smiles_list)} molecules")

    # Canonicalize in parallel
    t0 = time.time()
    results = batch_standardize_smiles(
        smiles_list, n_workers=args.workers, checkpoint_path=checkpoint_path,
    )
    elapsed = time.time() - t0

    # Count successes and failures
    n_ok = sum(1 for r in results if r is not None)
    n_fail = len(results) - n_ok
    print(f"Canonicalized {n_ok}/{len(results)} molecules in {elapsed:.1f}s "
          f"({n_fail} failures kept original SMILES)")

    # Compute InChIKey14 (connectivity layer) for each final SMILES
    print("Computing InChIKey14...")
    inchikey14_list = []
    for row, canon_smi in zip(rows, results):
        smi = canon_smi if canon_smi is not None else row["smiles"]
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                inchi = MolToInchi(mol)
                if inchi is not None:
                    inchikey14_list.append(InchiToInchiKey(inchi)[:14])
                else:
                    inchikey14_list.append("")
            else:
                inchikey14_list.append("")
        except Exception:
            inchikey14_list.append("")

    # Write output TSV — keep original SMILES on failure
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["cid", "smiles", "formula", "mass", "inchikey14"],
            delimiter="\t",
        )
        writer.writeheader()
        for row, canon_smi, ik14 in zip(rows, results, inchikey14_list):
            out_row = dict(row)
            if canon_smi is not None:
                out_row["smiles"] = canon_smi
            out_row["inchikey14"] = ik14
            writer.writerow(out_row)

    print(f"Wrote {len(rows)} rows to {args.output}")

    # Clean up checkpoint on success
    if checkpoint_path and os.path.exists(checkpoint_path):
        os.unlink(checkpoint_path)
        print(f"Removed checkpoint file {checkpoint_path}")


if __name__ == "__main__":
    main()
