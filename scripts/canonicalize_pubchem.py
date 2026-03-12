#!/usr/bin/env python
"""Canonicalize SMILES in a TSV file using ChEMBL Structure Pipeline.

Reads any TSV with a 'smiles' column, standardizes and canonicalizes
SMILES in parallel, adds an inchikey14 column, and preserves all
original columns. Molecules that fail standardization keep their
original SMILES.

Processes data in streaming chunks to keep memory usage constant
regardless of file size.

Usage:
    python canonicalize_pubchem.py input.tsv output.tsv --workers 16
    python canonicalize_pubchem.py input.tsv output.tsv --batch-size 500000
    python canonicalize_pubchem.py input.tsv output.tsv --resume
"""
import argparse
import csv
import json
import os
import sys
import time

from chembl_structure_pipeline.parallel import batch_standardize_smiles_with_inchikey


def _load_progress(progress_path):
    """Load number of rows already written from a progress file."""
    if progress_path and os.path.exists(progress_path):
        with open(progress_path) as f:
            data = json.load(f)
            return data.get("rows_done", 0)
    return 0


def _save_progress(progress_path, rows_done):
    """Save progress to file."""
    if progress_path:
        with open(progress_path, "w") as f:
            json.dump({"rows_done": rows_done}, f)


def main():
    parser = argparse.ArgumentParser(
        description="Canonicalize SMILES with ChEMBL Structure Pipeline",
    )
    parser.add_argument("input", help="Input TSV with a 'smiles' column")
    parser.add_argument("output", help="Output TSV with canonicalized SMILES")
    parser.add_argument(
        "--workers", "-w", type=int, default=None,
        help="Number of worker processes (default: all CPUs)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=500_000,
        help="Rows per batch (default: 500000). Controls memory usage.",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from where a previous run left off",
    )
    args = parser.parse_args()

    progress_path = args.output + ".progress.json"

    # Detect header from input
    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        input_fieldnames = list(reader.fieldnames)

    if "smiles" not in input_fieldnames:
        print("Error: input TSV must have a 'smiles' column", file=sys.stderr)
        sys.exit(1)

    output_fieldnames = input_fieldnames + (
        ["inchikey14"] if "inchikey14" not in input_fieldnames else []
    )

    # Resume support: skip already-processed rows
    rows_done = 0
    if args.resume:
        rows_done = _load_progress(progress_path)
        if rows_done > 0:
            print(f"Resuming from row {rows_done}")

    # Open output: append if resuming, write fresh otherwise
    write_mode = "a" if rows_done > 0 else "w"

    t0 = time.time()
    total_ok = 0
    total_fail = 0
    total_rows = 0

    with open(args.input, newline="") as fin, \
         open(args.output, write_mode, newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=output_fieldnames, delimiter="\t")

        if write_mode == "w":
            writer.writeheader()

        # Skip already-processed rows if resuming
        for _ in range(rows_done):
            try:
                next(reader)
            except StopIteration:
                break

        # Process in batches
        batch_num = 0
        while True:
            # Read a batch of rows
            batch = []
            for _ in range(args.batch_size):
                try:
                    batch.append(next(reader))
                except StopIteration:
                    break
            if not batch:
                break

            batch_num += 1
            smiles_list = [row["smiles"] for row in batch]
            print(f"Batch {batch_num}: processing {len(batch)} molecules "
                  f"(rows {rows_done + 1}–{rows_done + len(batch)})...")

            # Canonicalize and compute InChIKey14 in parallel
            results = batch_standardize_smiles_with_inchikey(
                smiles_list, n_workers=args.workers,
            )

            n_ok = sum(1 for canon, _ in results if canon is not None)
            n_fail = len(results) - n_ok
            total_ok += n_ok
            total_fail += n_fail

            # Write results
            for row, (canon_smi, ik14) in zip(batch, results):
                out_row = dict(row)
                if canon_smi is not None:
                    out_row["smiles"] = canon_smi
                out_row["inchikey14"] = ik14
                writer.writerow(out_row)

            rows_done += len(batch)
            total_rows += len(batch)
            fout.flush()

            # Save progress after each batch
            _save_progress(progress_path, rows_done)

            elapsed = time.time() - t0
            rate = total_rows / elapsed if elapsed > 0 else 0
            print(f"  Done: {n_ok}/{len(batch)} ok, {n_fail} failed | "
                  f"Total: {rows_done} rows, {rate:.0f} mol/s")

    elapsed = time.time() - t0
    print(f"\nFinished: {total_ok + total_fail} molecules in {elapsed:.1f}s "
          f"({total_ok} ok, {total_fail} failures kept original SMILES)")
    print(f"Output: {args.output}")

    # Clean up progress file on success
    if os.path.exists(progress_path):
        os.unlink(progress_path)


if __name__ == "__main__":
    main()
