#!/usr/bin/env python3
"""CLI entry point for Enamine REAL processing.

Reads CXSMILES from stdin (piped from lbzip2/bzcat), canonicalizes
molecules, computes properties and fingerprints, and writes aligned
output files.

Usage:
    lbzip2 -dc -n 2 HAC_11_21.cxsmiles.bz2 | \\
        python enamine_process.py \\
            --output-dir /persistent/processed/HAC_11_21/ \\
            --scratch-dir $SCRATCHDIR/HAC_11_21/ \\
            --workers 16

    # Resume an interrupted run:
    lbzip2 -dc -n 2 HAC_11_21.cxsmiles.bz2 | \\
        python enamine_process.py \\
            --output-dir /persistent/processed/HAC_11_21/ \\
            --scratch-dir $SCRATCHDIR/HAC_11_21/ \\
            --workers 16 --resume
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Process Enamine REAL CXSMILES: canonicalize + fingerprint",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Persistent storage dir for binary outputs (fps, bitsums, meta, checkpoint)",
    )
    parser.add_argument(
        "--scratch-dir", default=None,
        help="Scratch dir for TSV outputs (defaults to --output-dir)",
    )
    parser.add_argument(
        "--workers", "-w", type=int, default=16,
        help="Number of RDKit worker processes (default: 16)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=100_000,
        help="Lines to read per batch from stdin (default: 100000)",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=5_000,
        help="SMILES per worker chunk (default: 5000)",
    )
    parser.add_argument(
        "--chunk-timeout", type=int, default=3600,
        help="Seconds before killing a hung worker chunk (default: 3600)",
    )
    parser.add_argument(
        "--smiles-column", type=int, default=0,
        help="0-based column index for SMILES in input (default: 0)",
    )
    parser.add_argument(
        "--skip-header", action="store_true",
        help="Skip the first line of input (header row)",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from checkpoint (skips already-processed lines)",
    )
    args = parser.parse_args()

    # Import here to avoid slow RDKit import when just checking --help
    from chembl_structure_pipeline.enamine.processor import process_enamine_stream

    # Read from stdin
    stats = process_enamine_stream(
        input_lines=sys.stdin,
        output_dir=args.output_dir,
        scratch_dir=args.scratch_dir,
        n_workers=args.workers,
        batch_size=args.batch_size,
        chunk_size=args.chunk_size,
        chunk_timeout=args.chunk_timeout,
        smiles_column=args.smiles_column,
        skip_header=args.skip_header,
        resume_from=0 if not args.resume else 0,  # processor auto-loads checkpoint
    )

    print(f"\nOutput: {args.output_dir}")
    if args.scratch_dir:
        print(f"Scratch (stage out needed): {args.scratch_dir}")


if __name__ == "__main__":
    main()
