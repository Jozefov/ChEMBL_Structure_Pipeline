"""Parallel batch processing for ChEMBL Structure Pipeline.

All functions accept lists of string inputs (SMILES or molblocks),
distribute work across multiple processes to bypass the GIL, and
return results in input order. RDKit Mol objects never cross process
boundaries — only strings are serialized.

Checkpoint support: pass ``checkpoint_path`` to save progress after
each chunk. If a job is interrupted, re-running with the same
checkpoint file automatically skips already-processed items.
"""
import json
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Checkpoint helpers (JSON Lines format — handles multi-line molblocks
# and tuple results cleanly)
# ---------------------------------------------------------------------------

def _load_checkpoint(path):
    """Load completed indices and their results from a checkpoint file.

    Returns:
        dict mapping int index → result (string, tuple, or None)
    """
    completed = {}
    if path and os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                    idx = rec["i"]
                    completed[idx] = rec["r"] if rec["s"] == "ok" else None
                except (json.JSONDecodeError, KeyError):
                    continue
    return completed


def _append_checkpoint(path, results_with_indices):
    """Atomically append a batch of results to the checkpoint file.

    Writes to a temp file first, then appends to the real file to
    minimise corruption risk if the process is killed mid-write.
    """
    if not path:
        return
    dirn = os.path.dirname(path) or "."
    fd, tmp = tempfile.mkstemp(dir=dirn, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            for idx, result in results_with_indices:
                status = "ok" if result is not None else "err"
                rec = {"i": idx, "s": status, "r": result}
                f.write(json.dumps(rec) + "\n")
        # Append tmp contents to the real checkpoint
        with open(tmp) as src, open(path, "a") as dst:
            dst.write(src.read())
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Worker functions — executed in child processes
# ---------------------------------------------------------------------------
# Each worker imports the module inside the function body so that
# module-level globals (_normalizer, _tautomer_enumerator, cached
# salts/solvents) are initialised fresh in every worker process.
# This is required for the ``spawn`` start method (macOS default).

def _worker_standardize_smiles(chunk):
    """Process a list of (index, smiles) tuples."""
    from . import standardizer
    results = []
    for idx, smi in chunk:
        try:
            results.append((idx, standardizer.standardize_and_canonicalize_smiles(smi)))
        except Exception:
            results.append((idx, None))
    return results


def _worker_standardize_molblocks(chunk):
    from . import standardizer
    results = []
    for idx, ctab in chunk:
        try:
            results.append((idx, standardizer.standardize_molblock(ctab)))
        except Exception:
            results.append((idx, None))
    return results


def _worker_get_parent_molblocks(chunk):
    from . import standardizer
    results = []
    for idx, ctab in chunk:
        try:
            mb, exclude = standardizer.get_parent_molblock(ctab)
            results.append((idx, (mb, exclude)))
        except Exception:
            results.append((idx, (None, None)))
    return results


def _worker_check_molblocks(chunk):
    from .checker import check_molblock
    results = []
    for idx, mb in chunk:
        try:
            results.append((idx, check_molblock(mb)))
        except Exception:
            results.append((idx, ((7, "Processing error"),)))
    return results


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------

def _chunkify(indexed_items, n_chunks):
    """Split a list of (index, item) into n roughly equal chunks."""
    k, remainder = divmod(len(indexed_items), n_chunks)
    chunks = []
    start = 0
    for i in range(n_chunks):
        end = start + k + (1 if i < remainder else 0)
        if start < end:
            chunks.append(indexed_items[start:end])
        start = end
    return chunks


# ---------------------------------------------------------------------------
# Generic batch dispatcher
# ---------------------------------------------------------------------------

def _batch_process(items, worker_fn, n_workers, chunk_size, checkpoint_path,
                   default_result=None):
    """Generic parallel batch processor with checkpoint support.

    Args:
        items: list of input strings
        worker_fn: worker function that accepts [(index, item), ...]
        n_workers: number of worker processes (None = cpu_count)
        chunk_size: items per chunk (None = auto)
        checkpoint_path: path to checkpoint TSV (None = no checkpointing)
        default_result: fallback value for missing results

    Returns:
        list of results in input order
    """
    if not items:
        return []

    n_total = len(items)

    # Load checkpoint if it exists
    completed = _load_checkpoint(checkpoint_path) if checkpoint_path else {}
    if completed:
        print(f"Checkpoint: {len(completed)}/{n_total} items already processed, "
              f"resuming {n_total - len(completed)} remaining")

    # Build work queue: only items not yet completed
    work = [(i, items[i]) for i in range(n_total) if i not in completed]

    if not work:
        # Everything already done — assemble from checkpoint
        return [completed.get(i, default_result) for i in range(n_total)]

    if n_workers is None:
        n_workers = os.cpu_count() or 1
    n_workers = min(n_workers, len(work))

    # Serial fallback for small workloads
    if n_workers <= 1:
        chunk_results = worker_fn(work)
        if checkpoint_path:
            _append_checkpoint(checkpoint_path, chunk_results)
        for idx, result in chunk_results:
            completed[idx] = result
        return [completed.get(i, default_result) for i in range(n_total)]

    # Determine chunking
    if chunk_size is not None:
        n_chunks = max(1, -(-len(work) // chunk_size))
    else:
        # Use more chunks than workers for better progress granularity
        # (each chunk completion triggers a checkpoint write)
        n_chunks = min(len(work), n_workers * 4)

    chunks = _chunkify(work, n_chunks)

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for chunk_result in executor.map(worker_fn, chunks):
            # Save progress after each chunk
            if checkpoint_path:
                _append_checkpoint(checkpoint_path, chunk_result)
            for idx, result in chunk_result:
                completed[idx] = result
            done = len(completed)
            if done < n_total:
                print(f"Progress: {done}/{n_total} "
                      f"({100 * done / n_total:.1f}%)")

    print(f"Done: {n_total}/{n_total} items processed")
    return [completed.get(i, default_result) for i in range(n_total)]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def batch_standardize_smiles(
    smiles_list: List[str],
    n_workers: Optional[int] = None,
    chunk_size: Optional[int] = None,
    checkpoint_path: Optional[str] = None,
) -> List[Optional[str]]:
    """Standardize and canonicalize a list of SMILES in parallel.

    Args:
        smiles_list: Input SMILES strings.
        n_workers: Number of worker processes. Defaults to os.cpu_count().
        chunk_size: Items per chunk. Defaults to auto.
        checkpoint_path: Path to checkpoint file for resume support.
            None disables checkpointing.

    Returns:
        List of canonical SMILES (None for failures), same order as input.
    """
    return _batch_process(
        smiles_list, _worker_standardize_smiles,
        n_workers, chunk_size, checkpoint_path,
        default_result=None,
    )


def batch_standardize_molblocks(
    molblocks: List[str],
    n_workers: Optional[int] = None,
    chunk_size: Optional[int] = None,
    checkpoint_path: Optional[str] = None,
) -> List[Optional[str]]:
    """Standardize a list of molblocks in parallel.

    Args:
        molblocks: Input molblock strings.
        n_workers: Number of worker processes. Defaults to os.cpu_count().
        chunk_size: Items per chunk. Defaults to auto.
        checkpoint_path: Path to checkpoint file for resume support.

    Returns:
        List of standardized molblocks (None for failures), same order as input.
    """
    return _batch_process(
        molblocks, _worker_standardize_molblocks,
        n_workers, chunk_size, checkpoint_path,
        default_result=None,
    )


def batch_get_parent_molblocks(
    molblocks: List[str],
    n_workers: Optional[int] = None,
    chunk_size: Optional[int] = None,
    checkpoint_path: Optional[str] = None,
) -> List[Tuple[Optional[str], Optional[bool]]]:
    """Extract parent molblocks (salt-stripped) in parallel.

    Args:
        molblocks: Input molblock strings.
        n_workers: Number of worker processes. Defaults to os.cpu_count().
        chunk_size: Items per chunk. Defaults to auto.
        checkpoint_path: Path to checkpoint file for resume support.

    Returns:
        List of (molblock, exclude_flag) tuples, same order as input.
    """
    return _batch_process(
        molblocks, _worker_get_parent_molblocks,
        n_workers, chunk_size, checkpoint_path,
        default_result=(None, None),
    )


def batch_check_molblocks(
    molblocks: List[str],
    n_workers: Optional[int] = None,
    chunk_size: Optional[int] = None,
    checkpoint_path: Optional[str] = None,
) -> List[tuple]:
    """Check a list of molblocks for structural issues in parallel.

    Args:
        molblocks: Input molblock strings.
        n_workers: Number of worker processes. Defaults to os.cpu_count().
        chunk_size: Items per chunk. Defaults to auto.
        checkpoint_path: Path to checkpoint file for resume support.

    Returns:
        List of penalty tuples, same order as input.
    """
    return _batch_process(
        molblocks, _worker_check_molblocks,
        n_workers, chunk_size, checkpoint_path,
        default_result=((7, "Processing error"),),
    )
