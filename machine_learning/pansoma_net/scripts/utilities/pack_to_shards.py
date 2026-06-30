#!/usr/bin/env python3
import os
import glob
import time
import argparse
import numpy as np
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def load_one(sample):
    """Helper for threaded loading: sample = (path, label)."""
    path, label = sample
    try:
        x = np.load(path)  # (6, 201, 100)
    except Exception as e:
        raise RuntimeError(f"Error loading {path}: {e}")
    return x, label


def parse_args():
    parser = argparse.ArgumentParser(
        description="Shard .npy dataset into mixed-class .npy shards"
    )

    parser.add_argument(
        "--src",
        type=str,
        required=True,
        help="Source root containing class subfolders (true/false)",
    )
    parser.add_argument(
        "--dst",
        type=str,
        required=True,
        help="Destination root for output shards",
    )
    parser.add_argument(
        "--shard_size",
        type=int,
        default=4096,
        help="Number of samples per shard",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads to use when loading .npy files",
    )
    parser.add_argument(
        "--classes",
        type=str,
        nargs="+",
        default=["false", "true"],
        help="Class folder names (default: false true)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for shuffling samples",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    SRC_ROOT = args.src
    DST_ROOT = args.dst
    SHARD_SIZE = args.shard_size
    N_THREADS = args.threads
    CLASSES = args.classes
    SEED = args.seed

    os.makedirs(DST_ROOT, exist_ok=True)

    log("Starting sharding script.")
    log(f"SRC_ROOT = {SRC_ROOT}")
    log(f"DST_ROOT = {DST_ROOT}")
    log(f"SHARD_SIZE = {SHARD_SIZE:,}")
    log(f"N_THREADS = {N_THREADS}")
    log(f"CLASSES = {CLASSES}")
    log(f"SEED = {SEED}")

    # -------------------------------
    # Step 1: Scan dataset
    # -------------------------------
    log("Scanning dataset directories...")

    all_samples = []  # list of (path, label)
    per_class_counts = {}

    for cls in CLASSES:
        cls_dir = os.path.join(SRC_ROOT, cls)
        files = sorted(glob.glob(os.path.join(cls_dir, "*.npy")))

        if len(CLASSES) == 2:
            # binary classification
            label = 1 if cls == "true" else 0
        else:
            # multi-class: label = index
            label = CLASSES.index(cls)

        for p in files:
            all_samples.append((p, label))

        per_class_counts[label] = per_class_counts.get(label, 0) + len(files)
        log(f"  Found {len(files):,} '{cls}' samples in {cls_dir} (label={label})")

    total = len(all_samples)
    log(f"Total samples found: {total:,}")
    if total == 0:
        log("No samples found. Exiting.")
        return

    # Report global label distribution
    log("Global label distribution:")
    for lbl, cnt in sorted(per_class_counts.items()):
        frac = cnt / total if total > 0 else 0.0
        log(f"  label {lbl}: {cnt:,} ({frac:.3%})")

    # -------------------------------
    # Step 2: Shuffle samples
    # -------------------------------
    log("Shuffling all samples to mix classes proportionally...")
    rng = np.random.default_rng(SEED)
    indices = np.arange(total)
    rng.shuffle(indices)
    all_samples = [all_samples[i] for i in indices]

    log("Beginning sharding into .npy files...")

    # -------------------------------
    # Step 3: Create shards
    # -------------------------------
    shard_idx = 0

    for i in range(0, total, SHARD_SIZE):
        chunk = all_samples[i:i + SHARD_SIZE]
        log(
            f"\n--- Creating shard {shard_idx} "
            f"({i:,} → {i + len(chunk):,}) containing {len(chunk)} samples ---"
        )

        xs = []
        ys = []

        # Multithreaded loading
        with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
            futures = [executor.submit(load_one, sample) for sample in chunk]

            for f in tqdm(
                as_completed(futures),
                total=len(futures),
                desc=f"Loading shard {shard_idx}",
                leave=False,
            ):
                try:
                    x, label = f.result()
                except Exception as e:
                    log(str(e))
                    raise
                xs.append(x)
                ys.append(label)

        # Stack
        log("Stacking arrays...")
        xs = np.stack(xs, axis=0)  # shape: (N, 6, 201, 100)
        ys = np.array(ys, dtype=np.int64)  # shape: (N,)

        # Basic label distribution check for this shard
        unique_labels, counts = np.unique(ys, return_counts=True)
        shard_dist_str = ", ".join(
            f"label {int(l)}: {int(c)} ({c / len(ys):.3%})"
            for l, c in zip(unique_labels, counts)
        )
        log(f"Shard {shard_idx} label distribution: {shard_dist_str}")

        # Save to .npy instead of .npz
        data_path = os.path.join(DST_ROOT, f"shard_{shard_idx:05d}_data.npy")
        labels_path = os.path.join(DST_ROOT, f"shard_{shard_idx:05d}_labels.npy")

        log(f"Saving data to {data_path}")
        np.save(data_path, xs)

        log(f"Saving labels to {labels_path}")
        np.save(labels_path, ys)

        log(
            f"Saved shard {shard_idx}: "
            f"data shape={xs.shape}, labels shape={ys.shape}"
        )
        shard_idx += 1

    log("\nSharding completed successfully!")
    log(f"Total shards created: {shard_idx}")


if __name__ == "__main__":
    main()
