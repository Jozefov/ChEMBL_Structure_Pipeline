"""CLI for parallel batch processing of molecules.

Usage:
    python -m chembl_structure_pipeline standardize input.smi output.smi --workers 16
    python -m chembl_structure_pipeline standardize input.smi output.smi --checkpoint
    python -m chembl_structure_pipeline standardize input.smi output.smi --resume
    python -m chembl_structure_pipeline check input.sdf output.tsv --workers 16
"""
import argparse
import os
import sys
import time


def _read_smiles(path):
    """Read SMILES from a .smi file (one SMILES per line, optional tab-separated name)."""
    smiles = []
    names = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            smiles.append(parts[0].strip())
            names.append(parts[1].strip() if len(parts) > 1 else "")
    return smiles, names


def _read_sdf(path):
    """Read molblocks from an SDF file (split on $$$$)."""
    molblocks = []
    with open(path) as f:
        content = f.read()
    for block in content.split("$$$$"):
        block = block.strip()
        if block:
            molblocks.append(block + "\n")
    return molblocks


def _write_smiles(path, results, names=None):
    """Write SMILES results to a .smi file."""
    with open(path, "w") as f:
        for i, smi in enumerate(results):
            if smi is None:
                smi = ""
            if names and names[i]:
                f.write(f"{smi}\t{names[i]}\n")
            else:
                f.write(f"{smi}\n")


def _write_sdf(path, molblocks):
    """Write molblocks to an SDF file."""
    with open(path, "w") as f:
        for mb in molblocks:
            if mb is None:
                mb = ""
            f.write(mb)
            if not mb.endswith("\n"):
                f.write("\n")
            f.write("$$$$\n")


def cmd_standardize(args):
    from .parallel import batch_standardize_smiles, batch_standardize_molblocks

    checkpoint_path = None
    if args.checkpoint or args.resume:
        checkpoint_path = args.checkpoint_file or (args.output + ".checkpoint.jsonl")

    if args.resume and checkpoint_path and not os.path.exists(checkpoint_path):
        print(f"Error: checkpoint file {checkpoint_path} not found", file=sys.stderr)
        sys.exit(1)

    ext = os.path.splitext(args.input)[1].lower()

    t0 = time.time()

    if ext in (".smi", ".smiles", ".csv", ".txt"):
        smiles, names = _read_smiles(args.input)
        print(f"Read {len(smiles)} SMILES from {args.input}")
        results = batch_standardize_smiles(
            smiles, n_workers=args.workers, checkpoint_path=checkpoint_path,
        )
        _write_smiles(args.output, results, names)
    elif ext in (".sdf", ".mol"):
        molblocks = _read_sdf(args.input)
        print(f"Read {len(molblocks)} molblocks from {args.input}")
        results = batch_standardize_molblocks(
            molblocks, n_workers=args.workers, checkpoint_path=checkpoint_path,
        )
        _write_sdf(args.output, results)
    else:
        print(f"Error: unsupported file format '{ext}'. Use .smi or .sdf",
              file=sys.stderr)
        sys.exit(1)

    elapsed = time.time() - t0
    print(f"Wrote {len(results)} results to {args.output} in {elapsed:.1f}s")

    # Clean up checkpoint on successful completion
    if checkpoint_path and os.path.exists(checkpoint_path):
        os.unlink(checkpoint_path)
        print(f"Removed checkpoint file {checkpoint_path}")


def cmd_check(args):
    from .parallel import batch_check_molblocks

    checkpoint_path = None
    if args.checkpoint or args.resume:
        checkpoint_path = args.checkpoint_file or (args.output + ".checkpoint.jsonl")

    molblocks = _read_sdf(args.input)
    print(f"Read {len(molblocks)} molblocks from {args.input}")

    t0 = time.time()
    results = batch_check_molblocks(
        molblocks, n_workers=args.workers, checkpoint_path=checkpoint_path,
    )
    elapsed = time.time() - t0

    with open(args.output, "w") as f:
        f.write("index\tmax_penalty\tissues\n")
        for i, penalties in enumerate(results):
            if penalties:
                max_pen = max(p[0] for p in penalties)
                issues = "; ".join(f"{p[0]}:{p[1]}" for p in penalties)
            else:
                max_pen = 0
                issues = ""
            f.write(f"{i}\t{max_pen}\t{issues}\n")

    print(f"Wrote {len(results)} check results to {args.output} in {elapsed:.1f}s")

    if checkpoint_path and os.path.exists(checkpoint_path):
        os.unlink(checkpoint_path)


def cmd_get_parent(args):
    from .parallel import batch_get_parent_molblocks

    checkpoint_path = None
    if args.checkpoint or args.resume:
        checkpoint_path = args.checkpoint_file or (args.output + ".checkpoint.jsonl")

    molblocks = _read_sdf(args.input)
    print(f"Read {len(molblocks)} molblocks from {args.input}")

    t0 = time.time()
    results = batch_get_parent_molblocks(
        molblocks, n_workers=args.workers, checkpoint_path=checkpoint_path,
    )
    elapsed = time.time() - t0

    parent_molblocks = [r[0] for r in results]
    _write_sdf(args.output, parent_molblocks)

    print(f"Wrote {len(results)} parent molblocks to {args.output} in {elapsed:.1f}s")

    if checkpoint_path and os.path.exists(checkpoint_path):
        os.unlink(checkpoint_path)


def main():
    parser = argparse.ArgumentParser(
        prog="chembl_structure_pipeline",
        description="ChEMBL Structure Pipeline — parallel batch processing",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Common arguments
    def add_common_args(sub):
        sub.add_argument("input", help="Input file (.smi or .sdf)")
        sub.add_argument("output", help="Output file")
        sub.add_argument(
            "--workers", "-w", type=int, default=None,
            help="Number of worker processes (default: all CPUs)",
        )
        sub.add_argument(
            "--checkpoint", action="store_true",
            help="Enable checkpointing for resume support",
        )
        sub.add_argument(
            "--resume", action="store_true",
            help="Resume from existing checkpoint file",
        )
        sub.add_argument(
            "--checkpoint-file", default=None,
            help="Custom checkpoint file path (default: <output>.checkpoint.jsonl)",
        )

    # standardize
    p_std = subparsers.add_parser(
        "standardize", help="Standardize and canonicalize molecules",
    )
    add_common_args(p_std)
    p_std.set_defaults(func=cmd_standardize)

    # check
    p_chk = subparsers.add_parser(
        "check", help="Check molblocks for structural issues",
    )
    add_common_args(p_chk)
    p_chk.set_defaults(func=cmd_check)

    # get-parent
    p_par = subparsers.add_parser(
        "get-parent", help="Extract parent molblocks (strip salts/solvents)",
    )
    add_common_args(p_par)
    p_par.set_defaults(func=cmd_get_parent)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
