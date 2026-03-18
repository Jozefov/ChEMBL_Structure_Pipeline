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
import time
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
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


def _worker_standardize_smiles_with_inchikey(chunk):
    """Process a list of (index, smiles) tuples, returning (canon_smi, inchikey14).

    No per-molecule timeout here — hung C-level RDKit code can't be
    interrupted by SIGALRM. Instead, _batch_process() uses per-chunk
    process-level timeouts via future.result(timeout=...) which kills
    the entire worker process if it hangs.
    """
    from . import standardizer
    from rdkit import Chem
    from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

    results = []
    for idx, smi in chunk:
        canon = None
        try:
            canon = standardizer.standardize_and_canonicalize_smiles(smi)
        except Exception:
            canon = None

        # Compute InChIKey14 from canonical (or original) SMILES
        use_smi = canon if canon is not None else smi
        ik14 = ""
        try:
            mol = Chem.MolFromSmiles(use_smi)
            if mol is not None:
                inchi = MolToInchi(mol)
                if inchi is not None:
                    ik14 = InchiToInchiKey(inchi)[:14]
        except Exception:
            pass
        results.append((idx, (canon, ik14)))
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
                   default_result=None, chunk_timeout=1800):
    """Generic parallel batch processor with checkpoint support.

    Uses submit() + as_completed() so that a hung chunk does not block
    progress from other chunks.  Each chunk gets a process-level timeout
    (default 300 s).  If a chunk times out, its worker process is killed
    by the executor and all molecules in that chunk get default_result.

    Args:
        items: list of input strings
        worker_fn: worker function that accepts [(index, item), ...]
        n_workers: number of worker processes (None = cpu_count)
        chunk_size: items per chunk (None = auto)
        checkpoint_path: path to checkpoint TSV (None = no checkpointing)
        default_result: fallback value for missing results
        chunk_timeout: seconds to wait for a single chunk before killing
            the worker (default 1800). None disables timeouts.

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
        # 10K per chunk balances load vs overhead.  If a chunk times out
        # only these 10K molecules are lost, not the whole batch.
        chunk_size = 10_000
        n_chunks = max(1, -(-len(work) // chunk_size))

    chunks = _chunkify(work, n_chunks)

    # Submit chunks in waves of n_workers to avoid the timeout-from-
    # submission-time bug.  Previously all chunks were submitted up front
    # and the timeout was measured from submission, not execution start.
    # Chunks queued behind busy workers would sit in the queue, exceed the
    # timeout before they even started, and get killed as "timed out".
    #
    # Now we submit at most n_workers chunks at a time, wait for any to
    # finish, and refill the pool.  Every in-flight chunk is actually
    # running on a worker, so the timeout is meaningful.
    future_to_chunk = {}
    submit_times = {}
    n_timed_out = 0
    n_ok = 0
    n_fail = 0

    executor = ProcessPoolExecutor(max_workers=n_workers)
    try:
        # Initial fill: submit up to n_workers chunks
        for chunk in chunks[:n_workers]:
            fut = executor.submit(worker_fn, chunk)
            future_to_chunk[fut] = chunk
            submit_times[fut] = time.monotonic()
        next_chunk_idx = n_workers

        pending = set(future_to_chunk.keys())

        while pending:
            # Wait for any future to complete, polling every 30s
            done_set, still_pending = wait(
                pending, timeout=30, return_when=FIRST_COMPLETED,
            )

            # Process completed futures (returned out-of-order)
            for fut in done_set:
                chunk = future_to_chunk[fut]
                try:
                    chunk_result = fut.result()
                    if checkpoint_path:
                        _append_checkpoint(checkpoint_path, chunk_result)
                    for idx, result in chunk_result:
                        completed[idx] = result
                    n_ok += len(chunk)
                except Exception as exc:
                    n_timed_out += len(chunk)
                    n_fail += len(chunk)
                    failed = [(idx, default_result) for idx, _ in chunk]
                    if checkpoint_path:
                        _append_checkpoint(checkpoint_path, failed)
                    for idx, result in failed:
                        completed[idx] = result
                    print(f"WARNING: chunk of {len(chunk)} molecules "
                          f"crashed: {exc} — marked as failed")

            pending = still_pending

            # Check for timed-out futures still pending
            if chunk_timeout is not None and pending:
                now = time.monotonic()
                timed_out = {f for f in pending
                             if now - submit_times[f] > chunk_timeout}
                for fut in timed_out:
                    chunk = future_to_chunk[fut]
                    n_timed_out += len(chunk)
                    n_fail += len(chunk)
                    failed = [(idx, default_result) for idx, _ in chunk]
                    if checkpoint_path:
                        _append_checkpoint(checkpoint_path, failed)
                    for idx, result in failed:
                        completed[idx] = result
                    print(f"WARNING: chunk of {len(chunk)} molecules "
                          f"timed out after {chunk_timeout}s — "
                          f"marked as failed")
                    fut.cancel()
                pending -= timed_out

            # Refill: submit new chunks to replace completed/timed-out ones
            free_slots = n_workers - len(pending)
            while free_slots > 0 and next_chunk_idx < len(chunks):
                chunk = chunks[next_chunk_idx]
                next_chunk_idx += 1
                fut = executor.submit(worker_fn, chunk)
                future_to_chunk[fut] = chunk
                submit_times[fut] = time.monotonic()
                pending.add(fut)
                free_slots -= 1

            # Report progress
            done_count = len(completed)
            if done_count < n_total and (done_set or n_timed_out):
                print(f"  Progress: {done_count}/{n_total} "
                      f"({100 * done_count / n_total:.1f}%) | "
                      f"{n_ok} ok, {n_fail} failed")
    finally:
        # Force-kill any hung worker processes (they're stuck in C code
        # and won't respond to gentle shutdown). This is safe because
        # we've already checkpointed all completed results above.
        try:
            for proc in executor._processes.values():
                if proc.is_alive():
                    proc.kill()
        except Exception:
            pass
        executor.shutdown(wait=False, cancel_futures=True)

    if n_timed_out:
        print(f"WARNING: {n_timed_out} molecules timed out or crashed total")
    print(f"Done: {n_total}/{n_total} items processed "
          f"({n_ok} ok, {n_fail} failed)")
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


def batch_standardize_smiles_with_inchikey(
    smiles_list: List[str],
    n_workers: Optional[int] = None,
    chunk_size: Optional[int] = None,
    checkpoint_path: Optional[str] = None,
) -> List[Tuple[Optional[str], str]]:
    """Standardize SMILES and compute InChIKey14 in parallel.

    Combines canonicalization and InChIKey computation in the same worker
    to avoid a serial InChIKey pass after parallel standardization.

    Returns:
        List of (canonical_smiles_or_None, inchikey14) tuples, same order as input.
    """
    return _batch_process(
        smiles_list, _worker_standardize_smiles_with_inchikey,
        n_workers, chunk_size, checkpoint_path,
        default_result=(None, ""),
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
