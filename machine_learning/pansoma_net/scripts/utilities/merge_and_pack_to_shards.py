#!/usr/bin/env python3
import os
import re
import glob
import time
import argparse
import numpy as np
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def load_one(sample):
    """
    Helper for threaded loading.
    sample = (path, label)
    """
    path, label = sample
    try:
        x = np.load(path)  # expected shape: (6, 201, 100)
    except Exception as e:
        raise RuntimeError(f"Error loading {path}: {e}")
    return x, label


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Merge per-chromosome SNV 6-channel datasets into combined NPY shards "
            "for train/val splits (chr1 -> val, chr2-22 -> train by default)."
        )
    )

    parser.add_argument(
        "--src-base",
        type=str,
        required=True,
        help="Source base (e.g. /scratch/.../6ch_training_data_SNV)",
    )
    parser.add_argument(
        "--dest-base",
        type=str,
        required=True,
        help=(
            "Destination base for merged shards "
            "(e.g. /scratch/.../6ch_training_data_SNV/ALL_chr_merged_REAL_sharded)"
        ),
    )
    parser.add_argument(
        "--shard-size",
        type=int,
        default=4096,
        help="Number of samples per shard (per split).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads to use when loading .npy files.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for shuffling samples within each split.",
    )
    parser.add_argument(
        "--classes",
        type=str,
        nargs="+",
        default=["false", "true"],
        help="Class folder names (default: false true).",
    )
    parser.add_argument(
        "--val-chr",
        type=str,
        nargs="+",
        default=["chr1"],
        help="Chromosomes to route into val split (default: chr1).",
    )
    parser.add_argument(
        "--train-chr",
        type=str,
        nargs="+",
        default=[
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
        ],
        help="Chromosomes to route into train split (default: chr2–chr22).",
    )

    return parser.parse_args()


def discover_chr_dirs(src_base: str):
    """
    Find all SNV_chr*_6chan_tensor_dataset directories under src_base.
    Returns list of (chrom, dir_path).
    """
    pattern = os.path.join(src_base, "SNV_chr*_6chan_tensor_dataset")
    dirs = sorted(glob.glob(pattern))
    chr_dirs = []

    for d in dirs:
        bn = os.path.basename(d)
        m = re.match(r"^SNV_(chr[^_]+)_6chan_tensor_dataset$", bn)
        if not m:
            log(f"Skipping {d}: does not look like SNV_chr*_6chan_tensor_dataset")
            continue
        chrom = m.group(1)
        chr_dirs.append((chrom, d))

    return chr_dirs


def build_route_map(val_chr, train_chr):
    """
    Build chromosome -> split map.
    """
    route = {}
    for c in val_chr:
        route[c] = "val"
    for c in train_chr:
        route[c] = "train"
    return route


def scan_all_samples(
    src_base: str, route_map, classes, manifest_path, merged_ndjson_path
):
    """
    Walk all SNV_chr*_6chan_tensor_dataset directories, route each file into
    train/val based on chromosome, and collect:
      - samples_by_split["train"] = [(path, label_int), ...]
      - samples_by_split["val"]   = [(path, label_int), ...]
    Also builds a manifest and merged NDJSON.
    """
    chr_dirs = discover_chr_dirs(src_base)
    if not chr_dirs:
        raise RuntimeError(f"No SNV_chr*_6chan_tensor_dataset dirs found under {src_base}")

    log("Discovered per-chromosome dirs:")
    for chrom, d in chr_dirs:
        log(f"  {chrom}: {d}")

    # Map class name -> label int (binary or multi-class)
    cls_to_label = {}
    if len(classes) == 2:
        # Binary: treat "true" as 1 and "false" as 0 by convention
        for cls in classes:
            cls_to_label[cls] = 1 if cls == "true" else 0
    else:
        # Multi-class: index in classes list
        for idx, cls in enumerate(classes):
            cls_to_label[cls] = idx

    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    samples_by_split = {"train": [], "val": []}
    per_split_counts = {"train": 0, "val": 0}

    with open(manifest_path, "w") as mf, open(merged_ndjson_path, "w") as nd_all:
        for chrom, d in chr_dirs:
            split = route_map.get(chrom)
            if split is None:
                log(f"↷ Skipping {d} ({chrom} not in route map)")
                continue

            log(f"→ Routing {chrom} from {d} to split = '{split}'")

            # Collect .npy files from origin splits train/val/test and class folders
            for origin_split in ("train", "val", "test"):
                for cls in classes:
                    src_dir = os.path.join(d, origin_split, cls)
                    if not os.path.isdir(src_dir):
                        continue

                    files = sorted(glob.glob(os.path.join(src_dir, "*.npy")))
                    if not files:
                        continue

                    label_int = cls_to_label[cls]
                    for p in files:
                        samples_by_split[split].append((p, label_int))
                        per_split_counts[split] += 1
                        # manifest: original_path, routed_split, class_name, chrom
                        mf.write(f"{p}\t{split}\t{cls}\t{chrom}\n")

            # Append this chromosome's NDJSON (if present) into the merged one
            ndjson_path = os.path.join(d, "classification_summary.ndjson")
            if os.path.isfile(ndjson_path) and os.path.getsize(ndjson_path) > 0:
                with open(ndjson_path, "r") as nd_in:
                    for line in nd_in:
                        nd_all.write(line)

    log("Finished scanning all chromosomes.")
    for split in ("train", "val"):
        log(f"  {split}: {per_split_counts[split]:,} samples")

    return samples_by_split


def shard_one_split(
    split_name: str,
    samples,
    dest_base: str,
    shard_size: int,
    n_threads: int,
    seed: int,
):
    """
    Given a list of (path, label) samples belonging to one split (train or val),
    shuffle and write shard_XXXXX_data.npy and shard_XXXXX_labels.npy under
    dest_base/split_name.
    """
    if not samples:
        log(f"No samples for split '{split_name}', skipping sharding.")
        return 0

    log(f"\n==== Sharding split '{split_name}' ====")
    log(f"Total samples before shuffling: {len(samples):,}")

    rng = np.random.default_rng(seed)
    idx = np.arange(len(samples))
    rng.shuffle(idx)
    samples = [samples[i] for i in idx]
    log(f"Samples for '{split_name}' have been shuffled.")

    split_dest = os.path.join(dest_base, split_name)
    os.makedirs(split_dest, exist_ok=True)

    shard_idx = 0
    total = len(samples)

    for i in range(0, total, shard_size):
        chunk = samples[i : i + shard_size]
        log(
            f"\n--- Creating {split_name} shard {shard_idx} "
            f"({i:,} → {i + len(chunk):,}) containing {len(chunk)} samples ---"
        )

        xs = []
        ys = []

        # Multithreaded loading
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = [executor.submit(load_one, sample) for sample in chunk]

            for f in tqdm(
                as_completed(futures),
                total=len(futures),
                desc=f"Loading {split_name} shard {shard_idx}",
                leave=False,
            ):
                x, label = f.result()
                xs.append(x)
                ys.append(label)

        log("Stacking arrays...")
        xs = np.stack(xs, axis=0)  # (N, 6, 201, 100)
        ys = np.array(ys, dtype=np.int64)  # (N,)

        # Label distribution check for this shard
        unique_labels, counts = np.unique(ys, return_counts=True)
        shard_dist_str = ", ".join(
            f"label {int(l)}: {int(c)} ({c / len(ys):.3%})"
            for l, c in zip(unique_labels, counts)
        )
        log(f"{split_name} shard {shard_idx} label distribution: {shard_dist_str}")

        # >>> PATCHED NAMES TO MATCH NPYShardDataset <<<
        data_path = os.path.join(split_dest, f"shard_{shard_idx:05d}_data.npy")
        labels_path = os.path.join(split_dest, f"shard_{shard_idx:05d}_labels.npy")

        log(f"Saving data to {data_path}")
        np.save(data_path, xs)

        log(f"Saving labels to {labels_path}")
        np.save(labels_path, ys)

        log(
            f"Saved {split_name} shard {shard_idx}: "
            f"data shape={xs.shape}, labels shape={ys.shape}"
        )
        shard_idx += 1

    log(f"\nSharding for split '{split_name}' completed. Total shards: {shard_idx}")
    return shard_idx


def main():
    args = parse_args()

    SRC_BASE = args.src_base
    DEST_BASE = args.dest_base
    SHARD_SIZE = args.shard_size
    N_THREADS = args.threads
    SEED = args.seed
    CLASSES = args.classes
    VAL_CHR = args.val_chr
    TRAIN_CHR = args.train_chr

    os.makedirs(DEST_BASE, exist_ok=True)

    log("Starting merge+shard pipeline.")
    log(f"SRC_BASE   = {SRC_BASE}")
    log(f"DEST_BASE  = {DEST_BASE}")
    log(f"SHARD_SIZE = {SHARD_SIZE:,}")
    log(f"N_THREADS  = {N_THREADS}")
    log(f"CLASSES    = {CLASSES}")
    log(f"VAL_CHR    = {VAL_CHR}")
    log(f"TRAIN_CHR  = {TRAIN_CHR}")
    log(f"SEED       = {SEED}")

    merged_ndjson_path = os.path.join(DEST_BASE, "classification_summary_merged.ndjson")
    manifest_path = os.path.join(DEST_BASE, "placement_manifest.tsv")

    route_map = build_route_map(VAL_CHR, TRAIN_CHR)

    # -------------------------------
    # Step 1: Scan all per-chr dirs & route samples
    # -------------------------------
    log("\nScanning per-chromosome datasets and building sample lists...")
    samples_by_split = scan_all_samples(
        SRC_BASE,
        route_map,
        CLASSES,
        manifest_path,
        merged_ndjson_path,
    )

    # -------------------------------
    # Step 2: Shard each split separately
    # -------------------------------
    total_train_shards = shard_one_split(
        "train",
        samples_by_split["train"],
        DEST_BASE,
        SHARD_SIZE,
        N_THREADS,
        seed=SEED,
    )

    total_val_shards = shard_one_split(
        "val",
        samples_by_split["val"],
        DEST_BASE,
        SHARD_SIZE,
        N_THREADS,
        seed=SEED + 1,  # different shuffle for val if desired
    )

    log("\nAll done!")
    log(f"Total train shards: {total_train_shards}")
    log(f"Total val shards:   {total_val_shards}")
    log(f"Merged NDJSON:      {merged_ndjson_path}")
    log(f"Manifest:           {manifest_path}")


if __name__ == "__main__":
    main()
