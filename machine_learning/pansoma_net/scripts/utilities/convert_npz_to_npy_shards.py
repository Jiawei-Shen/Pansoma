#!/usr/bin/env python3
import argparse
import os
import numpy as np
from glob import glob
from multiprocessing import Pool, cpu_count


def detect_xy_keys(npz):
    keys = list(npz.files)
    key_set = set(keys)

    if {"data", "labels"} <= key_set:
        return "data", "labels"
    if {"x", "y"} <= key_set:
        return "x", "y"
    if set(keys) == {"arr_0", "arr_1"}:
        return "arr_0", "arr_1"

    raise RuntimeError(f"Unsupported key set {keys}")


def convert_one(args):
    """One conversion task executed in a separate process"""
    npz_path, overwrite, remove_npz = args

    try:
        dir_name = os.path.dirname(npz_path)
        base = os.path.splitext(os.path.basename(npz_path))[0]

        out_x = os.path.join(dir_name, base + "_x.npy")
        out_y = os.path.join(dir_name, base + "_y.npy")

        if os.path.exists(out_x) and os.path.exists(out_y) and not overwrite:
            return f"[SKIP] {npz_path}"

        with np.load(npz_path) as data:
            x_key, y_key = detect_xy_keys(data)
            x_arr = data[x_key]
            y_arr = data[y_key]

        np.save(out_x, x_arr)
        np.save(out_y, y_arr)

        if remove_npz:
            os.remove(npz_path)

        return f"[OK] Converted {npz_path}"

    except Exception as e:
        return f"[ERR] {npz_path}: {e}"


def convert_split(root, split, pattern, workers, overwrite, remove_npz):
    split_dir = os.path.join(root, split)
    paths = sorted(glob(os.path.join(split_dir, pattern)))

    if not paths:
        print(f"[WARN] No NPZ found in {split_dir}")
        return

    print(f"[INFO] Converting {len(paths)} NPZ shards in {split_dir} using {workers} workers")

    tasks = [(p, overwrite, remove_npz) for p in paths]

    with Pool(workers) as pool:
        for msg in pool.imap_unordered(convert_one, tasks):
            print(msg, flush=True)

    print(f"[DONE] Split {split}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=str, help="root folder containing train/ val/")
    ap.add_argument("--splits", nargs="+", default=["train", "val"])
    ap.add_argument("--pattern", type=str, default="shard_*.npz")
    ap.add_argument("--workers", type=int, default=cpu_count())
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--remove-npz", action="store_true")
    args = ap.parse_args()

    root = os.path.abspath(args.root)

    for split in args.splits:
        convert_split(
            root=root,
            split=split,
            pattern=args.pattern,
            workers=args.workers,
            overwrite=args.overwrite,
            remove_npz=args.remove_npz,
        )


if __name__ == "__main__":
    main()
