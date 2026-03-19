"""Streaming batch processor for Enamine REAL datasets.

Reads SMILES from an iterator (typically sys.stdin piped from lbzip2/bzcat),
dispatches work to a pool of worker processes, and writes results
incrementally to aligned output files:

  molecules.tsv.gz   — text properties (canonical SMILES, InChIKey, etc.)
  fingerprints.bin   — packed Morgan fingerprints (N × 256 bytes)
  bitsums.bin        — per-molecule bit sums (N × 2 bytes, uint16)
  errors.tsv.gz      — failed molecules
  meta.json          — processing statistics

Cannot reuse ``parallel._batch_process()`` because it requires all items
in memory upfront. Instead, this module replicates the wave-based
ProcessPoolExecutor submission pattern from ``parallel.py`` but feeds
batches from a streaming iterator and writes results incrementally.
"""

import gzip
import json
import os
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

import numpy as np

from . import MORGAN_RADIUS, MORGAN_NBITS, PACKED_WIDTH, TSV_COLUMNS, ERROR_COLUMNS
from .worker import worker_process_chunk


# Popcount lookup table for vectorised bitsum computation
# Same pattern as compute_candidates_pool_fps.py from MSnGym
_POPCOUNT_TABLE = np.unpackbits(
    np.arange(256, dtype=np.uint8)[:, None], axis=1
).sum(axis=1).astype(np.float32)


def _compute_bitsums_batch(packed_fps_list):
    """Compute bit sums from a list of packed fingerprint bytes.

    Args:
        packed_fps_list: List of bytes objects (each 256 bytes).

    Returns:
        np.ndarray of uint16 bit sums, shape (len(packed_fps_list),).
    """
    if not packed_fps_list:
        return np.empty((0,), dtype=np.uint16)
    arr = np.frombuffer(b"".join(packed_fps_list), dtype=np.uint8).reshape(-1, PACKED_WIDTH)
    sums = _POPCOUNT_TABLE[arr].sum(axis=1)
    return sums.astype(np.uint16)


def _save_checkpoint(checkpoint_path, lines_processed, n_ok, n_errors):
    """Atomically save checkpoint to JSON file."""
    if not checkpoint_path:
        return
    dirn = os.path.dirname(checkpoint_path) or "."
    fd, tmp = tempfile.mkstemp(dir=dirn, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump({
                "lines_processed": lines_processed,
                "n_ok": n_ok,
                "n_errors": n_errors,
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            }, f)
        os.replace(tmp, checkpoint_path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _load_checkpoint(checkpoint_path):
    """Load checkpoint if it exists.

    Returns:
        dict with keys lines_processed, n_ok, n_errors, or empty dict.
    """
    if checkpoint_path and os.path.exists(checkpoint_path):
        with open(checkpoint_path) as f:
            return json.load(f)
    return {}


def _chunkify(items, chunk_size):
    """Split a list into chunks of at most chunk_size."""
    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


def process_enamine_stream(
    input_lines,
    output_dir,
    scratch_dir=None,
    n_workers=16,
    batch_size=100_000,
    chunk_size=5_000,
    chunk_timeout=3600,
    smiles_column=0,
    skip_header=False,
    resume_from=0,
):
    """Stream-process Enamine CXSMILES input to aligned output files.

    Args:
        input_lines: Iterator of raw input lines (from stdin/pipe).
        output_dir: Directory for persistent outputs (fingerprints.bin,
            bitsums.bin, meta.json, checkpoint.json). Must exist.
        scratch_dir: Directory for temporary outputs (molecules.tsv.gz,
            errors.tsv.gz). Falls back to output_dir if None.
        n_workers: Number of worker processes for RDKit.
        batch_size: Lines to read per batch from stdin.
        chunk_size: SMILES per worker chunk.
        chunk_timeout: Seconds before killing a hung worker chunk.
        smiles_column: 0-based column index for SMILES in input TSV.
        skip_header: Skip the first input line.
        resume_from: Number of input lines already processed (for resume).

    Returns:
        dict with processing statistics.
    """
    if scratch_dir is None:
        scratch_dir = output_dir

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(scratch_dir, exist_ok=True)

    # Output file paths
    tsv_path = os.path.join(scratch_dir, "molecules.tsv.gz")
    errors_path = os.path.join(scratch_dir, "errors.tsv.gz")
    fps_path = os.path.join(output_dir, "fingerprints.bin")
    sums_path = os.path.join(output_dir, "bitsums.bin")
    checkpoint_path = os.path.join(output_dir, "checkpoint.json")

    # Load checkpoint for resume
    ckpt = _load_checkpoint(checkpoint_path)
    if resume_from == 0 and ckpt:
        resume_from = ckpt.get("lines_processed", 0)

    total_ok = ckpt.get("n_ok", 0) if resume_from > 0 else 0
    total_errors = ckpt.get("n_errors", 0) if resume_from > 0 else 0
    lines_processed = resume_from

    # Open mode depends on resume
    is_resume = resume_from > 0
    tsv_mode = "ab" if is_resume else "wb"
    fps_mode = "ab" if is_resume else "wb"

    line_iter = iter(input_lines)

    # Skip header if requested
    if skip_header and not is_resume:
        try:
            next(line_iter)
        except StopIteration:
            pass

    # Skip already-processed lines on resume
    if resume_from > 0:
        skip_start = time.time()
        skipped = 0
        # If we skipped the header in the original run, skip it again
        if skip_header:
            try:
                next(line_iter)
            except StopIteration:
                pass
        for _ in range(resume_from):
            try:
                next(line_iter)
            except StopIteration:
                break
            skipped += 1
        skip_time = time.time() - skip_start
        print(f"Resume: skipped {skipped} lines in {skip_time:.1f}s "
              f"(continuing from line {resume_from})")

    t0 = time.time()

    # Open all output files
    f_tsv = gzip.open(tsv_path, tsv_mode)
    f_errors = gzip.open(errors_path, tsv_mode)
    f_fps = open(fps_path, fps_mode)
    f_sums = open(sums_path, fps_mode)

    try:
        # Write headers for fresh run
        if not is_resume:
            header = "\t".join(TSV_COLUMNS) + "\n"
            f_tsv.write(header.encode())
            err_header = "\t".join(ERROR_COLUMNS) + "\n"
            f_errors.write(err_header.encode())

        batch_num = 0
        while True:
            # Read a batch of lines from stdin
            batch_lines = []
            for _ in range(batch_size):
                try:
                    line = next(line_iter)
                except StopIteration:
                    break
                line = line.rstrip("\n\r")
                if not line:
                    continue
                batch_lines.append(line)

            if not batch_lines:
                break

            batch_num += 1

            # Extract SMILES from the specified column
            batch_smiles = []
            for line in batch_lines:
                parts = line.split("\t")
                if smiles_column < len(parts):
                    batch_smiles.append(parts[smiles_column])
                else:
                    batch_smiles.append("")

            print(f"Batch {batch_num}: {len(batch_smiles)} molecules "
                  f"(lines {lines_processed + 1}\u2013"
                  f"{lines_processed + len(batch_smiles)})...")

            # Dispatch to workers
            batch_ok, batch_err = _process_batch(
                batch_smiles=batch_smiles,
                batch_lines=batch_lines,
                n_workers=n_workers,
                chunk_size=chunk_size,
                chunk_timeout=chunk_timeout,
                f_tsv=f_tsv,
                f_fps=f_fps,
                f_sums=f_sums,
                f_errors=f_errors,
            )

            total_ok += batch_ok
            total_errors += batch_err
            lines_processed += len(batch_smiles)

            # Flush all files
            f_tsv.flush()
            f_errors.flush()
            f_fps.flush()
            f_sums.flush()

            # Checkpoint
            _save_checkpoint(checkpoint_path, lines_processed, total_ok, total_errors)

            elapsed = time.time() - t0
            rate = (total_ok + total_errors - (ckpt.get("n_ok", 0) + ckpt.get("n_errors", 0)) if is_resume else total_ok + total_errors) / max(elapsed, 1)
            print(f"  Done: {batch_ok}/{len(batch_smiles)} ok, {batch_err} errors | "
                  f"Total: {lines_processed} lines, {total_ok} ok, "
                  f"{total_errors} errors, {rate:.0f} mol/s")

    finally:
        f_tsv.close()
        f_errors.close()
        f_fps.close()
        f_sums.close()

    elapsed = time.time() - t0

    stats = {
        "lines_processed": lines_processed,
        "n_molecules": total_ok,
        "n_errors": total_errors,
        "processing_time_s": round(elapsed, 1),
        "mol_per_s": round((total_ok + total_errors) / max(elapsed, 1), 1),
        "n_workers": n_workers,
        "batch_size": batch_size,
        "chunk_size": chunk_size,
        "chunk_timeout": chunk_timeout,
        "fp_radius": MORGAN_RADIUS,
        "fp_nbits": MORGAN_NBITS,
        "fp_packed_width": PACKED_WIDTH,
    }

    # Write meta.json
    meta_path = os.path.join(output_dir, "meta.json")
    with open(meta_path, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"\nFinished: {lines_processed} lines in {elapsed:.1f}s "
          f"({total_ok} ok, {total_errors} errors)")

    return stats


def _process_batch(
    batch_smiles,
    batch_lines,
    n_workers,
    chunk_size,
    chunk_timeout,
    f_tsv,
    f_fps,
    f_sums,
    f_errors,
):
    """Process one batch of SMILES through the worker pool.

    Uses the wave-based submission pattern from parallel.py: submit
    n_workers chunks at a time, wait for any to complete, refill.

    Returns:
        (n_ok, n_errors) counts for this batch.
    """
    indexed = list(enumerate(batch_smiles))
    chunks = list(_chunkify(indexed, chunk_size))

    if not chunks:
        return 0, 0

    # Collect results keyed by index
    results = {}
    n_ok = 0
    n_errors = 0

    effective_workers = min(n_workers, len(chunks))

    executor = ProcessPoolExecutor(max_workers=effective_workers)
    try:
        future_to_chunk = {}
        submit_times = {}

        # Initial fill
        for chunk in chunks[:effective_workers]:
            fut = executor.submit(worker_process_chunk, chunk)
            future_to_chunk[fut] = chunk
            submit_times[fut] = time.monotonic()
        next_idx = effective_workers

        pending = set(future_to_chunk.keys())

        while pending:
            done_set, still_pending = wait(
                pending, timeout=30, return_when=FIRST_COMPLETED,
            )

            for fut in done_set:
                chunk = future_to_chunk[fut]
                try:
                    chunk_results = fut.result()
                    for idx, result in chunk_results:
                        results[idx] = result
                except Exception as exc:
                    # Entire chunk failed
                    for idx, _ in chunk:
                        results[idx] = None
                    print(f"WARNING: chunk of {len(chunk)} molecules "
                          f"crashed: {exc}")

            pending = still_pending

            # Check for timed-out futures
            if chunk_timeout is not None and pending:
                now = time.monotonic()
                timed_out = {f for f in pending
                             if now - submit_times[f] > chunk_timeout}
                for fut in timed_out:
                    chunk = future_to_chunk[fut]
                    for idx, _ in chunk:
                        results[idx] = None
                    print(f"WARNING: chunk of {len(chunk)} molecules "
                          f"timed out after {chunk_timeout}s")
                    fut.cancel()
                pending -= timed_out

            # Refill
            free_slots = effective_workers - len(pending)
            while free_slots > 0 and next_idx < len(chunks):
                chunk = chunks[next_idx]
                next_idx += 1
                fut = executor.submit(worker_process_chunk, chunk)
                future_to_chunk[fut] = chunk
                submit_times[fut] = time.monotonic()
                pending.add(fut)
                free_slots -= 1

    finally:
        try:
            for proc in executor._processes.values():
                if proc.is_alive():
                    proc.kill()
        except Exception:
            pass
        executor.shutdown(wait=False, cancel_futures=True)

    # Write results in order (preserves alignment)
    ok_fps = []
    for i in range(len(batch_smiles)):
        result = results.get(i)
        if result is not None:
            canon_smi, inchikey, inchikey14, formula, mass, packed_fp = result
            # Write TSV line
            line = f"{canon_smi}\t{inchikey}\t{inchikey14}\t{formula}\t{mass}\n"
            f_tsv.write(line.encode())
            # Collect fingerprint for binary write
            f_fps.write(packed_fp)
            ok_fps.append(packed_fp)
            n_ok += 1
        else:
            # Write to errors file
            original = batch_lines[i] if i < len(batch_lines) else ""
            # Extract just the SMILES for the error log
            parts = original.split("\t")
            smi = parts[0] if parts else original
            err_line = f"{smi}\tprocessing_failed\n"
            f_errors.write(err_line.encode())
            n_errors += 1

    # Compute and write bitsums in bulk
    if ok_fps:
        bitsums = _compute_bitsums_batch(ok_fps)
        f_sums.write(bitsums.tobytes())

    return n_ok, n_errors
