#!/usr/bin/env python3
import argparse
import os
import re
import sys
import numpy as np
from typing import Dict, Tuple, List

def log(msg: str) -> None:
    print(msg, flush=True)

def detect_shard_bases(data_dir: str) -> Dict[int, str]:
    """
    Find files like:
      <prefix><00000>_data.npy
    and return:
      shard_idx -> base path without "_data.npy"/"_labels.npy"
    Example:
      COLO829T_ONT_chr10_00000_data.npy -> idx=0, base=.../COLO829T_ONT_chr10_00000
    """
    pattern = re.compile(r"^(.*?)(\d{5})_data\.npy$")
    mapping: Dict[int, str] = {}
    for fname in os.listdir(data_dir):
        m = pattern.match(fname)
        if not m:
            continue
        prefix, idx_str = m.group(1), m.group(2)
        idx = int(idx_str)
        base = os.path.join(data_dir, prefix + idx_str)
        if idx in mapping and mapping[idx] != base:
            log(f"[WARN] Multiple data files for shard {idx:05d}: "
                f"{os.path.basename(mapping[idx])} vs {fname}; using {fname}")
        mapping[idx] = base
    return mapping

def get_data_labels_paths(base: str) -> Tuple[str, str]:
    return base + "_data.npy", base + "_labels.npy"

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def write_filtered_shard(
    data_path: str,
    labels_path: str,
    out_data_path: str,
    out_labels_path: str,
    filter_fraction: float,
    false_label: int,
    keep_unknown: bool,
    unknown_label: int,
    seed: int,
    keep_order: bool,
    save_indices_path: str,
) -> Tuple[int, int, int]:
    """
    Returns: (N_in, N_out, N_false_removed)
    """
    if not os.path.isfile(labels_path):
        raise FileNotFoundError(f"Missing labels: {labels_path}")
    if not os.path.isfile(data_path):
        raise FileNotFoundError(f"Missing data: {data_path}")

    labels = np.load(labels_path, mmap_mode="r")
    if labels.ndim != 1:
        raise ValueError(f"Labels must be 1D. Got shape={labels.shape} from {labels_path}")
    N = labels.shape[0]

    # Figure out which indices are eligible to remove (false only)
    is_false = (labels == false_label)
    if keep_unknown:
        base_keep = np.ones(N, dtype=bool)  # keep everything initially
    else:
        base_keep = (labels != unknown_label)  # drop unknown entirely if requested

    false_idx = np.flatnonzero(is_false & base_keep)  # only consider falses that are otherwise kept
    n_false = false_idx.size

    # Compute how many falses to remove
    n_remove = int(round(filter_fraction * n_false))
    n_remove = max(0, min(n_remove, n_false))

    rng = np.random.default_rng(seed)
    if n_remove > 0:
        remove_idx = rng.choice(false_idx, size=n_remove, replace=False)
        base_keep[remove_idx] = False

    keep_idx = np.flatnonzero(base_keep)
    if keep_idx.size == 0:
        raise RuntimeError("After filtering, no samples remain. Adjust --filter-fraction / flags.")

    if not keep_order:
        rng.shuffle(keep_idx)

    # Load data via mmap and write filtered subset via memmap
    data = np.load(data_path, mmap_mode="r")
    if data.shape[0] != N:
        raise ValueError(f"Data/label length mismatch: data N={data.shape[0]} vs labels N={N}")

    outN = keep_idx.size
    out_shape = (outN,) + data.shape[1:]
    out_dtype = data.dtype

    # Write labels
    out_labels = np.asarray(labels[keep_idx], dtype=np.int8)
    np.save(out_labels_path, out_labels)

    # Write indices if requested (useful for mapping back to original)
    if save_indices_path is not None:
        np.save(save_indices_path, keep_idx.astype(np.int64))

    # Write data (memmap) to avoid huge RAM usage
    # NOTE: np.save() wants an ndarray; for big arrays we create a .npy-compatible file by:
    #  - creating a memmap using np.lib.format.open_memmap (writes valid .npy with header)
    from numpy.lib.format import open_memmap
    mm = open_memmap(out_data_path, mode="w+", dtype=out_dtype, shape=out_shape)

    # Copy rows in chunks to reduce overhead (and keep memory low)
    # Tune chunk size if you want
    chunk = 1024
    for i0 in range(0, outN, chunk):
        i1 = min(outN, i0 + chunk)
        sel = keep_idx[i0:i1]
        mm[i0:i1] = data[sel]

    # Flush memmap to disk
    del mm

    n_removed = n_remove
    return N, outN, n_removed

def main():
    ap = argparse.ArgumentParser(
        description="Filter (downsample) a portion of false (label=0) examples from shard *_data.npy / *_labels.npy."
    )
    ap.add_argument("--data-dir", required=True, help="Directory containing *_data.npy and *_labels.npy")
    ap.add_argument("--out-dir", required=True, help="Output directory for filtered shards")
    ap.add_argument(
        "--filter-fraction",
        type=float,
        default=0.80,
        help="Fraction of FALSE examples (label==--false-label) to REMOVE. Example 0.8 removes 80%% of falses.",
    )
    ap.add_argument("--false-label", type=int, default=0, help="Integer label value treated as 'false'")
    ap.add_argument("--unknown-label", type=int, default=-1, help="Unknown label value")
    ap.add_argument(
        "--keep-unknown",
        action="store_true",
        help="Keep unknown-label samples (default: keep).",
    )
    ap.add_argument(
        "--drop-unknown",
        action="store_true",
        help="Drop unknown-label samples entirely (overrides --keep-unknown).",
    )
    ap.add_argument("--seed", type=int, default=13, help="RNG seed for reproducible filtering")
    ap.add_argument(
        "--keep-order",
        action="store_true",
        help="Keep original order of samples (default: may shuffle kept indices).",
    )
    ap.add_argument(
        "--suffix",
        type=str,
        default="filtered",
        help="Suffix appended to base filename before _data/_labels. Example: base_00000_filtered_data.npy",
    )
    ap.add_argument(
        "--save-indices",
        action="store_true",
        help="Also save kept indices as *_kept_indices.npy per shard (maps filtered->original).",
    )
    args = ap.parse_args()

    if args.filter_fraction < 0.0 or args.filter_fraction > 1.0:
        sys.exit("--filter-fraction must be in [0, 1].")

    if not os.path.isdir(args.data_dir):
        sys.exit(f"Not a directory: {args.data_dir}")
    ensure_dir(args.out_dir)

    keep_unknown = True
    if args.drop_unknown:
        keep_unknown = False
    elif args.keep_unknown:
        keep_unknown = True  # explicit keep

    shard_bases = detect_shard_bases(args.data_dir)
    if not shard_bases:
        sys.exit(f"No *_data.npy shards found in: {args.data_dir}")

    log(f"Found {len(shard_bases)} shard(s) in {args.data_dir}")
    total_in = total_out = total_removed = 0

    for shard_idx in sorted(shard_bases.keys()):
        base = shard_bases[shard_idx]
        data_path, labels_path = get_data_labels_paths(base)

        # Output naming:
        #   <basename>_<suffix>_data.npy / _labels.npy
        base_name = os.path.basename(base)  # e.g., COLO829T_ONT_chr10_00000
        out_base = os.path.join(args.out_dir, f"{base_name}_{args.suffix}")

        out_data_path = out_base + "_data.npy"
        out_labels_path = out_base + "_labels.npy"
        out_idx_path = (out_base + "_kept_indices.npy") if args.save_indices else None

        log(f"\n[shard {shard_idx:05d}]")
        log(f"  in : {os.path.basename(data_path)} / {os.path.basename(labels_path)}")
        log(f"  out: {os.path.basename(out_data_path)} / {os.path.basename(out_labels_path)}")
        log(f"  filter-fraction={args.filter_fraction} on false-label={args.false_label} "
            f"(keep_unknown={keep_unknown}, unknown_label={args.unknown_label})")

        N, outN, removed = write_filtered_shard(
            data_path=data_path,
            labels_path=labels_path,
            out_data_path=out_data_path,
            out_labels_path=out_labels_path,
            filter_fraction=args.filter_fraction,
            false_label=args.false_label,
            keep_unknown=keep_unknown,
            unknown_label=args.unknown_label,
            seed=args.seed + shard_idx,  # vary per shard deterministically
            keep_order=args.keep_order,
            save_indices_path=out_idx_path,
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
