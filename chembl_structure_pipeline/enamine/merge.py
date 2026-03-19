"""Post-processing: merge per-file Enamine outputs into a global index.

After all (or a subset of) Enamine REAL files have been processed
individually, this module creates a unified index that enables
molecule lookup across all files without concatenating the large
binary fingerprint files.

Output files:
  file_offsets.json       — maps source dir name → (start_row, end_row)
  inchikey14_index.tsv.gz — sorted (inchikey14, file_id, local_row)
  global_meta.json        — aggregate statistics

Usage:
    python -m chembl_structure_pipeline.enamine.merge \\
        /path/to/processed/HAC_11_21 \\
        /path/to/processed/HAC_22_23 \\
        ... \\
        --output-dir /path/to/processed/merged/
"""

import argparse
import gzip
import json
import os
import sys

from . import PACKED_WIDTH


def merge_enamine_outputs(input_dirs, output_dir):
    """Create a global index across per-file Enamine outputs.

    Does NOT concatenate binary files (they can be multi-TB). Instead,
    creates an offset-based index so any molecule can be located by
    (file_id, local_row).

    Args:
        input_dirs: List of per-file output directories, each containing
            molecules.tsv.gz, fingerprints.bin, bitsums.bin, meta.json.
        output_dir: Directory for merged index outputs.

    Returns:
        dict with aggregate statistics.
    """
    os.makedirs(output_dir, exist_ok=True)

    file_offsets = {}
    global_row = 0
    total_molecules = 0
    total_errors = 0
    total_lines = 0

    # First pass: collect metadata and compute offsets
    for input_dir in sorted(input_dirs):
        dir_name = os.path.basename(input_dir)
        meta_path = os.path.join(input_dir, "meta.json")

        if not os.path.exists(meta_path):
            print(f"WARNING: skipping {input_dir} (no meta.json)")
            continue

        with open(meta_path) as f:
            meta = json.load(f)

        n_mol = meta.get("n_molecules", 0)
        n_err = meta.get("n_errors", 0)
        n_lines = meta.get("lines_processed", 0)

        # Verify binary file sizes match metadata
        fps_path = os.path.join(input_dir, "fingerprints.bin")
        sums_path = os.path.join(input_dir, "bitsums.bin")

        if os.path.exists(fps_path):
            fps_size = os.path.getsize(fps_path)
            expected_fps = n_mol * PACKED_WIDTH
            if fps_size != expected_fps:
                print(f"WARNING: {dir_name}/fingerprints.bin size mismatch: "
                      f"{fps_size} != {expected_fps} (expected for {n_mol} molecules)")

        if os.path.exists(sums_path):
            sums_size = os.path.getsize(sums_path)
            expected_sums = n_mol * 2
            if sums_size != expected_sums:
                print(f"WARNING: {dir_name}/bitsums.bin size mismatch: "
                      f"{sums_size} != {expected_sums}")

        file_offsets[dir_name] = {
            "dir": input_dir,
            "start_row": global_row,
            "end_row": global_row + n_mol,
            "n_molecules": n_mol,
            "n_errors": n_err,
            "n_lines": n_lines,
        }

        global_row += n_mol
        total_molecules += n_mol
        total_errors += n_err
        total_lines += n_lines

        print(f"  {dir_name}: {n_mol:,} molecules "
              f"(global rows {file_offsets[dir_name]['start_row']:,}"
              f"–{file_offsets[dir_name]['end_row']:,})")

    # Write file offsets
    offsets_path = os.path.join(output_dir, "file_offsets.json")
    with open(offsets_path, "w") as f:
        json.dump(file_offsets, f, indent=2)

    print(f"\nTotal: {total_molecules:,} molecules across "
          f"{len(file_offsets)} files")

    # Second pass: build InChIKey14 index
    print("\nBuilding InChIKey14 index...")
    index_path = os.path.join(output_dir, "inchikey14_index.tsv.gz")
    n_index = _build_inchikey14_index(file_offsets, index_path)

    # Write global metadata
    global_meta = {
        "total_molecules": total_molecules,
        "total_errors": total_errors,
        "total_lines_processed": total_lines,
        "n_files": len(file_offsets),
        "unique_inchikey14": n_index,
        "file_offsets_path": offsets_path,
        "index_path": index_path,
    }
    meta_path = os.path.join(output_dir, "global_meta.json")
    with open(meta_path, "w") as f:
        json.dump(global_meta, f, indent=2)

    print(f"\nMerge complete:")
    print(f"  File offsets: {offsets_path}")
    print(f"  InChIKey14 index: {index_path} ({n_index:,} entries)")
    print(f"  Global meta: {meta_path}")

    return global_meta


def _build_inchikey14_index(file_offsets, output_path):
    """Build sorted InChIKey14 → (file_id, local_row) index.

    Streams through each per-file molecules.tsv.gz, extracts
    inchikey14 column, and writes a sorted index.

    Args:
        file_offsets: dict from merge_enamine_outputs.
        output_path: Path for output .tsv.gz.

    Returns:
        Number of index entries written.
    """
    # Collect (inchikey14, file_id, local_row) tuples
    entries = []

    for dir_name, info in file_offsets.items():
        tsv_path = os.path.join(info["dir"], "molecules.tsv.gz")
        if not os.path.exists(tsv_path):
            print(f"  WARNING: {tsv_path} not found, skipping")
            continue

        local_row = 0
        with gzip.open(tsv_path, "rt") as f:
            header = f.readline()  # skip header
            for line in f:
                parts = line.rstrip("\n").split("\t")
                # inchikey14 is column index 2 (canonical_smiles, inchikey, inchikey14, ...)
                if len(parts) > 2 and parts[2]:
                    entries.append((parts[2], dir_name, local_row))
                local_row += 1

        print(f"  {dir_name}: {local_row:,} rows indexed")

    # Sort by inchikey14 for binary search
    entries.sort(key=lambda x: x[0])

    # Write sorted index
    with gzip.open(output_path, "wt") as f:
        f.write("inchikey14\tfile_id\tlocal_row\n")
        for ik14, file_id, local_row in entries:
            f.write(f"{ik14}\t{file_id}\t{local_row}\n")

    return len(entries)


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-file Enamine outputs into a global index",
    )
    parser.add_argument(
        "input_dirs", nargs="+",
        help="Per-file output directories to merge",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory for merged index outputs",
    )
    args = parser.parse_args()

    merge_enamine_outputs(args.input_dirs, args.output_dir)


if __name__ == "__main__":
    main()
