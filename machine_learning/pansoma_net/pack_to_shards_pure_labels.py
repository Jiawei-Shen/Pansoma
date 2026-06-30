#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np


def list_samples_by_class(split_dir):
    """
    Expecting structure like:

        split_dir/
          false/
            *.npy
          true/
            *.npy

    Returns:
        neg_paths, pos_paths
    """
    neg_paths = []
    pos_paths = []

    if not os.path.isdir(split_dir):
        return neg_paths, pos_paths

    for cls in sorted(os.listdir(split_dir)):
        p = os.path.join(split_dir, cls)
        if not os.path.isdir(p):
            continue

        name = cls.lower()
        # "true" -> positives; everything else -> negatives
        if "true" in name or name in {"1", "pos", "positive"}:
            target = pos_paths
        else:
            target = neg_paths

        for root, _, files in os.walk(p):
            for f in files:
                if f.endswith(".npy"):
                    target.append(os.path.join(root, f))

    neg_paths = sorted(neg_paths)
    pos_paths = sorted(pos_paths)
    return neg_paths, pos_paths


def infer_shape(first_path):
    a = np.load(first_path, mmap_mode="r")
    return a.shape, a.dtype


def write_shards_for_class(
    paths,
    labels_value,
    shard_size,
    split_outdir,
    shard_id_start,
    C,
    H,
    W,
    dtype,
    split_name,
    class_name,
):
    """
    Pack all `paths` (all same class) into shards of at most `shard_size`.

    Each shard is pure-class (all labels_value).

    Returns:
        next_shard_id
    """
    N = len(paths)
    if N == 0:
        print(f"[{split_name}] no {class_name} samples to write.")
        return shard_id_start

    print(
        f"[{split_name}] {class_name}: N={N} "
        f"shape=(6,{H},{W}) dtype={dtype}, shard_size={shard_size}"
    )

    shard_id = shard_id_start
    start = 0
    while start < N:
        end = min(start + shard_size, N)
        n_this = end - start

        base = f"shard_{shard_id:05d}"
        X_path = os.path.join(split_outdir, f"{base}_x.npy")
        y_path = os.path.join(split_outdir, f"{base}_y.npy")

        X = np.memmap(X_path, mode="w+", dtype=dtype, shape=(n_this, C, H, W))
        y = np.memmap(y_path, mode="w+", dtype=np.int64, shape=(n_this,))

        for i, p in enumerate(paths[start:end]):
            arr = np.load(p, mmap_mode="r")
            if arr.shape != (C, H, W):
                raise ValueError(
                    f"Shape mismatch at {p}: got {arr.shape}, expect {(C, H, W)}"
                )
            X[i] = arr

        # all labels in this shard are the same
        y[:] = labels_value

        del X
        del y

        print(
            f"  wrote {X_path} | {y_path} "
            f"[{n_this} {class_name} samples, shard_id={shard_id}]"
        )

        start = end
        shard_id += 1

    return shard_id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="dataset root containing train/ and val/")
    ap.add_argument(
        "-o",
        "--outdir",
        default=None,
        help="Output root for sharded dataset (default: same as root)",
    )
    ap.add_argument(
        "--split",
        choices=["train", "val"],
        default=None,
        help="If set, pack only this split; else pack both train and val.",
    )
    ap.add_argument(
        "--shard_size",
        type=int,
        default=4096,
        help="Max samples per shard file.",
    )
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    outdir = os.path.abspath(args.outdir or root)
    os.makedirs(outdir, exist_ok=True)

    splits = [args.split] if args.split else ["train", "val"]

    for split in splits:
        split_dir = os.path.join(root, split)
        if not os.path.isdir(split_dir):
            print(f"[skip] {split_dir} not found")
            continue

        neg_paths, pos_paths = list_samples_by_class(split_dir)

        if len(neg_paths) == 0 and len(pos_paths) == 0:
            print(f"[skip] no samples in {split_dir}")
            continue

        example_path = (neg_paths or pos_paths)[0]
        shape, dtype = infer_shape(example_path)
        if len(shape) != 3:
            raise ValueError(
                f"Expected sample shape (C,H,W), got {shape} at {example_path}"
            )
        C, H, W = shape
        if C != 6:
            print(
                f"[warn] Expected C=6 channels, but got C={C} in {example_path}",
                file=sys.stderr,
            )

        N_neg = len(neg_paths)
        N_pos = len(pos_paths)
        print(
            f"[{split}] total: neg={N_neg}, pos={N_pos}, "
            f"shape=(6,{H},{W}) dtype={dtype}"
        )

        split_outdir = os.path.join(outdir, split)
        os.makedirs(split_outdir, exist_ok=True)

        shard_id = 0

        # write all FALSE first
        shard_id = write_shards_for_class(
            paths=neg_paths,
            labels_value=0,
            shard_size=args.shard_size,
            split_outdir=split_outdir,
            shard_id_start=shard_id,
            C=C,
            H=H,
            W=W,
            dtype=dtype,
            split_name=split,
            class_name="neg(false)",
        )

        # then all TRUE
        shard_id = write_shards_for_class(
            paths=pos_paths,
            labels_value=1,
            shard_size=args.shard_size,
            split_outdir=split_outdir,
            shard_id_start=shard_id,
            C=C,
            H=H,
            W=W,
            dtype=dtype,
            split_name=split,
            class_name="pos(true)",
        )

        num_neg_shards = (N_neg + args.shard_size - 1) // args.shard_size
        num_pos_shards = (N_pos + args.shard_size - 1) // args.shard_size
        total_shards = num_neg_shards + num_pos_shards

        print(
            f"[{split}] DONE. "
            f"neg_shards={num_neg_shards}, pos_shards={num_pos_shards}, "
            f"total_shards={total_shards}, shard_size={args.shard_size}"
        )
        print(
            f"  (Shard indices: 0..{num_neg_shards-1} = all FALSE, "
            f"{num_neg_shards}..{total_shards-1} = all TRUE)\n"
        )


if __name__ == "__main__":
    main()
