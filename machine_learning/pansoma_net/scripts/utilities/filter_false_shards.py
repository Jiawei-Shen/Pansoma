#!/usr/bin/env python3
import argparse
import os
import re
import sys
import numpy as np
from typing import List, Tuple

def log(msg: str) -> None:
    print(msg, flush=True)

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def list_shard_bases(data_dir: str) -> List[str]:
    """
    Return a list of FULL base paths (without _data.npy) for every *_data.npy file.
    Example:
      .../COLO829T_ONT_chr10_00000_data.npy -> base .../COLO829T_ONT_chr10_00000
    """
    pattern = re.compile(r"^(.*\d{5})_data\.npy$")
    bases = []
    for fname in os.listdir(data_dir):
        m = pattern.match(fname)
        if not m:
            continue
        base = os.path.join(data_dir, m.group(1))
        bases.append(base)
    bases.sort()
    return bases

def write_filtered_one(
    base: str,
    out_dir: str,
    suffix: str,
    filter_fraction: float,
    false_label: int,
    unknown_label: int,
    keep_unknown: bool,
    seed: int,
    keep_order: bool,
    save_indices: bool,
    chunk_rows: int = 1024,
) -> Tuple[int, int, int]:
    """
    Filter one shard base:
      base_data.npy, base_labels.npy  -> out_base_data.npy, out_base_labels.npy
    Returns: (N_in, N_out, N_false_removed)
    """
    data_path = base + "_data.npy"
    labels_path = base + "_labels.npy"

    if not os.path.isfile(data_path):
        raise FileNotFoundError(f"Missing data: {data_path}")
    if not os.path.isfile(labels_path):
        raise FileNotFoundError(f"Missing labels: {labels_path}")

    labels = np.load(labels_path, mmap_mode="r")
    if labels.ndim != 1:
        raise ValueError(f"Labels must be 1D: {labels_path} shape={labels.shape}")
    N = labels.shape[0]

    # base_keep: keep all by default; optionally drop unknown
    if keep_unknown:
        keep_mask = np.ones(N, dtype=bool)
    else:
        keep_mask = (labels != unknown_label)

    is_false = (labels == false_label) & keep_mask
    false_idx = np.flatnonzero(is_false)
    n_false = false_idx.size

    n_remove = int(round(filter_fraction * n_false))
    n_remove = max(0, min(n_remove, n_false))

    rng = np.random.default_rng(seed)
    if n_remove > 0:
        remove_idx = rng.choice(false_idx, size=n_remove, replace=False)
        keep_mask[remove_idx] = False

    keep_idx = np.flatnonzero(keep_mask)
    if keep_idx.size == 0:
        raise RuntimeError(f"After filtering, no samples remain for {os.path.basename(base)}")

    if not keep_order:
        rng.shuffle(keep_idx)

    data = np.load(data_path, mmap_mode="r")
    if data.shape[0] != N:
        raise ValueError(f"Data/labels mismatch for {base}: data N={data.shape[0]} vs labels N={N}")

    out_base_name = os.path.basename(base) + f"_{suffix}"
    out_base = os.path.join(out_dir, out_base_name)
    out_data_path = out_base + "_data.npy"
    out_labels_path = out_base + "_labels.npy"
    out_idx_path = out_base + "_kept_indices.npy" if save_indices else None

    out_labels = np.asarray(labels[keep_idx], dtype=np.int8)
    np.save(out_labels_path, out_labels)

    if out_idx_path is not None:
        np.save(out_idx_path, keep_idx.astype(np.int64))

    # Write data with open_memmap to produce a valid .npy without loading whole array
    from numpy.lib.format import open_memmap
    outN = keep_idx.size
    out_shape = (outN,) + data.shape[1:]
    mm = open_memmap(out_data_path, mode="w+", dtype=data.dtype, shape=out_shape)

    for i0 in range(0, outN, chunk_rows):
        i1 = min(outN, i0 + chunk_rows)
        sel = keep_idx[i0:i1]
        mm[i0:i1] = data[sel]

    del mm
    return N, outN, n_remove

def main():
    ap = argparse.ArgumentParser(
        description="Filter (downsample) a portion of false (label=0) for EVERY *_data.npy in a directory."
    )
    ap.add_argument("--data-dir", required=True, help="Directory with *_data.npy + *_labels.npy")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--filter-fraction", type=float, default=0.80,
                    help="Fraction of FALSE examples to REMOVE (0.8 removes 80%% of falses)")
    ap.add_argument("--false-label", type=int, default=0)
    ap.add_argument("--unknown-label", type=int, default=-1)
    ap.add_argument("--drop-unknown", action="store_true", help="Drop unknown samples entirely")
    ap.add_argument("--keep-order", action="store_true", help="Keep original order (default shuffles)")
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--suffix", type=str, default="filtered", help="Suffix added to output base name")
    ap.add_argument("--save-indices", action="store_true",
                    help="Save *_kept_indices.npy per shard for mapping back to original")
    ap.add_argument("--chunk-rows", type=int, default=1024,
                    help="Rows per copy chunk when writing (memory/perf tradeoff)")
    args = ap.parse_args()

    if not (0.0 <= args.filter_fraction <= 1.0):
        sys.exit("--filter-fraction must be in [0, 1].")
    if not os.path.isdir(args.data_dir):
        sys.exit(f"Not a directory: {args.data_dir}")

    ensure_dir(args.out_dir)
    bases = list_shard_bases(args.data_dir)
    if not bases:
        sys.exit(f"No *_data.npy found in {args.data_dir}")

    log(f"Found {len(bases)} data shard file(s)")

    total_in = total_out = total_removed = 0
    for k, base in enumerate(bases):
        bn = os.path.basename(base)
        log(f"\n[{k+1}/{len(bases)}] {bn}")

        # deterministic per-file seed (so reruns are stable but different across files)
        seed = args.seed + (abs(hash(bn)) % 1_000_000)

        N, outN, removed = write_filtered_one(
            base=base,
            out_dir=args.out_dir,
            suffix=args.suffix,
            filter_fraction=args.filter_fraction,
            false_label=args.false_label,
            unknown_label=args.unknown_label,
            keep_unknown=(not args.drop_unknown),
            seed=seed,
            keep_order=args.keep_order,
            save_indices=args.save_indices,
            chunk_rows=args.chunk_rows,
        )
        log(f"  kept {outN:,}/{N:,} (removed false={removed:,})")

        total_in += N
        total_out += outN
        total_removed += removed

    log("\nDone.")
    log(f"Total in : {total_in:,}")
    log(f"Total out: {total_out:,}")
    log(f"Total false removed: {total_removed:,}")

if __name__ == "__main__":
    main()
