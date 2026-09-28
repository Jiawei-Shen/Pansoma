#!/usr/bin/env python3
"""Prepare binary 5-channel train/validation shards for PansomaNet."""

from __future__ import annotations

import argparse
import json
import os
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


DATA_SUFFIX = "_data.npy"
LABEL_SUFFIX = "_labels.npy"
KNOWN_LABELS = (0, 1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter labeled 5-channel tensor shards to binary labels 0/1 and "
            "prepare the train/ and val/ directories required by PansomaNet."
        )
    )
    parser.add_argument(
        "--input-dir",
        nargs="+",
        help=(
            "Tensor directories to split reproducibly at sample level. Use this "
            "mode for a quick test when chromosome-separated inputs are unavailable."
        ),
    )
    parser.add_argument(
        "--train-input-dir",
        nargs="+",
        help="Tensor directories assigned entirely to the training split.",
    )
    parser.add_argument(
        "--val-input-dir",
        nargs="+",
        help="Tensor directories assigned entirely to the validation split.",
    )
    parser.add_argument("--output-dir", required=True, help="Output dataset root.")
    parser.add_argument(
        "--val-fraction",
        type=float,
        default=0.2,
        help="Validation fraction for --input-dir mode (default: 0.2).",
    )
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=256,
        help="Samples copied per chunk to limit memory use (default: 256).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing generated .npy files in OUTPUT_DIR/train and val.",
    )
    return parser.parse_args()


def resolve_input_dirs(values: Optional[Iterable[str]], option: str) -> list[Path]:
    result = []
    for value in values or []:
        path = Path(value).expanduser().resolve()
        if not path.is_dir():
            raise SystemExit(f"ERROR: {option} is not a directory: {path}")
        result.append(path)
    return result


def discover_shards(input_dirs: Iterable[Path]) -> list[tuple[Path, Path, str]]:
    shards: list[tuple[Path, Path, str]] = []
    output_names: set[str] = set()
    for input_dir in input_dirs:
        directory_prefix = re.sub(r"[^A-Za-z0-9_.-]+", "_", input_dir.name)
        found = 0
        for data_path in sorted(input_dir.glob(f"*{DATA_SUFFIX}")):
            labels_path = data_path.with_name(
                data_path.name[: -len(DATA_SUFFIX)] + LABEL_SUFFIX
            )
            if not labels_path.is_file():
                raise SystemExit(f"ERROR: missing label shard for {data_path}")
            shard_name = data_path.name[: -len(DATA_SUFFIX)]
            output_name = f"{directory_prefix}_{shard_name}"
            if output_name in output_names:
                raise SystemExit(
                    "ERROR: duplicate output shard name from input directories: "
                    f"{output_name}"
                )
            output_names.add(output_name)
            shards.append((data_path, labels_path, output_name))
            found += 1
        if found == 0:
            raise SystemExit(f"ERROR: no *{DATA_SUFFIX} shards found in {input_dir}")
    return shards


def validate_shard(data_path: Path, labels_path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.load(data_path, mmap_mode="r")
    labels = np.load(labels_path, mmap_mode="r")
    if data.ndim != 4 or data.shape[1] != 5:
        raise SystemExit(
            f"ERROR: expected a 5-channel (N,5,H,W) shard, got {data.shape}: {data_path}"
        )
    if labels.ndim != 1 or data.shape[0] != labels.shape[0]:
        raise SystemExit(
            "ERROR: data/label shape mismatch: "
            f"{data_path}={data.shape}, {labels_path}={labels.shape}"
        )
    return data, labels


def random_split_indices(
    labels: np.ndarray,
    val_fraction: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    train_parts = []
    val_parts = []
    for label in KNOWN_LABELS:
        indices = np.flatnonzero(labels == label)
        rng.shuffle(indices)
        if len(indices) > 1:
            val_count = min(
                len(indices) - 1,
                max(1, round(len(indices) * val_fraction)),
            )
        else:
            val_count = 0
        val_parts.append(indices[:val_count])
        train_parts.append(indices[val_count:])

    train = np.concatenate(train_parts).astype(np.int64, copy=False)
    val = np.concatenate(val_parts).astype(np.int64, copy=False)
    # Keep selected indices in source order for efficient memory-mapped reads.
    # The trainer shuffles samples independently during every epoch.
    train.sort()
    val.sort()
    return train, val


def known_indices(labels: np.ndarray) -> np.ndarray:
    return np.flatnonzero(np.isin(labels, KNOWN_LABELS)).astype(np.int64, copy=False)


def write_subset(
    data: np.ndarray,
    labels: np.ndarray,
    indices: np.ndarray,
    output_dir: Path,
    output_name: str,
    chunk_size: int,
) -> None:
    if len(indices) == 0:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    data_output = output_dir / f"{output_name}{DATA_SUFFIX}"
    labels_output = output_dir / f"{output_name}{LABEL_SUFFIX}"
    data_tmp = data_output.with_name(f".{data_output.name}.tmp")
    labels_tmp = labels_output.with_name(f".{labels_output.name}.tmp")
    data_tmp.unlink(missing_ok=True)
    labels_tmp.unlink(missing_ok=True)

    out_data = np.lib.format.open_memmap(
        data_tmp,
        mode="w+",
        dtype=data.dtype,
        shape=(len(indices), *data.shape[1:]),
    )
    out_labels = np.lib.format.open_memmap(
        labels_tmp,
        mode="w+",
        dtype=np.int64,
        shape=(len(indices),),
    )
    for start in range(0, len(indices), chunk_size):
        end = min(start + chunk_size, len(indices))
        source_indices = indices[start:end]
        out_data[start:end] = data[source_indices]
        out_labels[start:end] = labels[source_indices]
    out_data.flush()
    out_labels.flush()
    del out_data, out_labels
    os.replace(data_tmp, data_output)
    os.replace(labels_tmp, labels_output)


def prepare_output(output_dir: Path, force: bool) -> None:
    generated = []
    for split in ("train", "val"):
        split_dir = output_dir / split
        if split_dir.is_dir():
            generated.extend(split_dir.glob("*.npy"))
    summary = output_dir / "split_summary.json"
    if summary.exists():
        generated.append(summary)
    if generated and not force:
        raise SystemExit(
            f"ERROR: generated files already exist under {output_dir}; use --force to replace them"
        )
    if force:
        for path in generated:
            path.unlink()


def update_counts(counter: Counter[int], values: np.ndarray) -> None:
    labels, counts = np.unique(values, return_counts=True)
    counter.update({int(label): int(count) for label, count in zip(labels, counts)})


def main() -> int:
    args = parse_args()
    random_dirs = resolve_input_dirs(args.input_dir, "--input-dir")
    train_dirs = resolve_input_dirs(args.train_input_dir, "--train-input-dir")
    val_dirs = resolve_input_dirs(args.val_input_dir, "--val-input-dir")

    random_mode = bool(random_dirs)
    explicit_mode = bool(train_dirs or val_dirs)
    if random_mode == explicit_mode:
        raise SystemExit(
            "ERROR: use either --input-dir, or both --train-input-dir and --val-input-dir"
        )
    if explicit_mode and (not train_dirs or not val_dirs):
        raise SystemExit(
            "ERROR: chromosome-separated mode requires both --train-input-dir and --val-input-dir"
        )
    if not 0.0 < args.val_fraction < 1.0:
        raise SystemExit("ERROR: --val-fraction must be between 0 and 1")
    if args.chunk_size < 1:
        raise SystemExit("ERROR: --chunk-size must be at least 1")
    if set(train_dirs).intersection(val_dirs):
        raise SystemExit("ERROR: the same tensor directory cannot be used for train and validation")

    output_dir = Path(args.output_dir).expanduser().resolve()
    prepare_output(output_dir, args.force)
    rng = np.random.default_rng(args.seed)
    totals = Counter()
    label_counts = {"input": Counter(), "train": Counter(), "val": Counter()}
    shard_summaries = []

    assignments: list[tuple[str, Path, Path, str]] = []
    if random_mode:
        assignments.extend(("random", *shard) for shard in discover_shards(random_dirs))
    else:
        assignments.extend(("train", *shard) for shard in discover_shards(train_dirs))
        assignments.extend(("val", *shard) for shard in discover_shards(val_dirs))

    for assignment, data_path, labels_path, output_name in assignments:
        data, labels = validate_shard(data_path, labels_path)
        update_counts(label_counts["input"], labels)
        known = int(np.isin(labels, KNOWN_LABELS).sum())
        totals["input"] += len(labels)
        totals["ignored"] += len(labels) - known

        if assignment == "random":
            train_indices, val_indices = random_split_indices(labels, args.val_fraction, rng)
        elif assignment == "train":
            train_indices, val_indices = known_indices(labels), np.empty(0, dtype=np.int64)
        else:
            train_indices, val_indices = np.empty(0, dtype=np.int64), known_indices(labels)

        for split, indices in (("train", train_indices), ("val", val_indices)):
            write_subset(
                data,
                labels,
                indices,
                output_dir / split,
                output_name,
                args.chunk_size,
            )
            totals[split] += len(indices)
            update_counts(label_counts[split], labels[indices])

        shard_summaries.append(
            {
                "data": str(data_path),
                "labels": str(labels_path),
                "input_samples": int(len(labels)),
                "ignored_samples": int(len(labels) - known),
                "train_samples": int(len(train_indices)),
                "val_samples": int(len(val_indices)),
            }
        )
        print(
            f"[shard] {data_path.name}: input={len(labels):,} "
            f"train={len(train_indices):,} val={len(val_indices):,} "
            f"ignored={len(labels) - known:,}"
        )

    if totals["train"] == 0 or totals["val"] == 0:
        raise SystemExit(
            "ERROR: both train and validation must contain at least one known label (0 or 1)"
        )

    summary = {
        "schema_version": 1,
        "mode": "random_sample_split" if random_mode else "explicit_input_split",
        "seed": args.seed if random_mode else None,
        "val_fraction": args.val_fraction if random_mode else None,
        "known_labels": list(KNOWN_LABELS),
        "totals": {key: int(value) for key, value in totals.items()},
        "label_counts": {
            split: {str(key): int(value) for key, value in sorted(counts.items())}
            for split, counts in label_counts.items()
        },
        "shards": shard_summaries,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "split_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(
        f"[summary] train={totals['train']:,} val={totals['val']:,} "
        f"ignored_non_linear={totals['ignored']:,} output={output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
