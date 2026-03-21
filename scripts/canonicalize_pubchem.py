#!/usr/bin/env python
"""Canonicalize SMILES in a TSV file using ChEMBL Structure Pipeline.

Reads any TSV with a 'smiles' column, standardizes and canonicalizes
SMILES in parallel, then computes all molecular properties from the
canonical SMILES: formula, monoisotopic mass, InChIKey, and InChIKey 2D
(first 14 characters). Molecules that fail canonicalization are dropped.

Output columns (fixed schema): smiles, formula, mass, inchikey, inchi_key_2D

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
import shutil
import sys
import time

from chembl_structure_pipeline.parallel import batch_standardize_full

OUTPUT_COLUMNS = ["smiles", "formula", "mass", "inchikey", "inchi_key_2D"]


def _load_progress(progress_path):
    """Load number of rows already written from a progress file."""
    if progress_path and os.path.exists(progress_path):
        with open(progress_path) as f:
            data = json.load(f)
            return data.get("rows_done", 0), data.get("rows_written", 0)
    return 0, 0


def _save_progress(progress_path, rows_done, rows_written):
    """Save progress to file."""
    if progress_path:
        with open(progress_path, "w") as f:
            json.dump({"rows_done": rows_done, "rows_written": rows_written}, f)


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
    parser.add_argument(
        "--checkpoint-dir", type=str, default=None,
        help="Directory on persistent storage to checkpoint output after every batch. "
             "Prevents data loss if the job crashes.",
    )
    parser.add_argument(
        "--checkpoint-name", type=str, default=None,
        help="Filename for checkpoint (default: basename of output). "
             "Use when scratch output name differs from desired persistent name.",
    )
    args = parser.parse_args()

    progress_path = args.output + ".progress.json"

    # Persistent checkpoint paths (on reliable storage, not scratch)
    if args.checkpoint_dir:
        os.makedirs(args.checkpoint_dir, exist_ok=True)
        basename = args.checkpoint_name or os.path.basename(args.output)
        checkpoint_output = os.path.join(args.checkpoint_dir, basename)
        checkpoint_progress = checkpoint_output + ".progress.json"
    else:
        checkpoint_output = None
        checkpoint_progress = None

    # Verify input has a smiles column
    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        input_fieldnames = list(reader.fieldnames)

    if "smiles" not in input_fieldnames:
        print("Error: input TSV must have a 'smiles' column", file=sys.stderr)
        sys.exit(1)

    # Resume support: skip already-processed rows
    rows_done = 0
    rows_written = 0
    if args.resume:
        rows_done, rows_written = _load_progress(progress_path)
        if rows_done > 0:
            print(f"Resuming from row {rows_done} ({rows_written} written)")

    # Open output: append if resuming, write fresh otherwise
    write_mode = "a" if rows_done > 0 else "w"

    t0 = time.time()
    total_ok = 0
    total_fail = 0
    total_rows = 0

    with open(args.input, newline="") as fin, \
         open(args.output, write_mode, newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=OUTPUT_COLUMNS, delimiter="\t")

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
                  f"(rows {rows_done + 1}\u2013{rows_done + len(batch)})...")

            # Canonicalize and compute all properties in parallel
            results = batch_standardize_full(
                smiles_list, n_workers=args.workers,
            )

            n_ok = 0
            n_fail = 0

            # Write results — drop failures
            for result in results:
                if result is None:
                    n_fail += 1
                    continue
                canon_smi, formula, mass, inchikey, inchi_key_2d = result
                writer.writerow({
                    "smiles": canon_smi,
                    "formula": formula,
                    "mass": mass,
                    "inchikey": inchikey,
                    "inchi_key_2D": inchi_key_2d,
                })
                n_ok += 1
                rows_written += 1

            total_ok += n_ok
            total_fail += n_fail

            rows_done += len(batch)
            total_rows += len(batch)
            fout.flush()

            # Save progress after each batch (scratch — fast)
            _save_progress(progress_path, rows_done, rows_written)

            # Checkpoint to persistent storage so we never lose work
            if checkpoint_output:
                shutil.copy2(args.output, checkpoint_output)
                _save_progress(checkpoint_progress, rows_done, rows_written)

            elapsed = time.time() - t0
            rate = total_rows / elapsed if elapsed > 0 else 0
            print(f"  Done: {n_ok}/{len(batch)} ok, {n_fail} dropped | "
                  f"Total: {rows_done} read, {rows_written} written, {rate:.0f} mol/s")

    elapsed = time.time() - t0
    print(f"\nFinished: {total_ok + total_fail} molecules in {elapsed:.1f}s "
          f"({total_ok} ok, {total_fail} dropped)")
    print(f"Output: {args.output} ({rows_written} rows)")

    # Clean up progress file on success
    if os.path.exists(progress_path):
        os.unlink(progress_path)


if __name__ == "__main__":
    main()
