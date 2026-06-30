#!/usr/bin/env python3
"""
Fast batch inference over a directory of .npy files (recursively), supporting:

1) SINGLE mode: each .npy is one tensor of shape (6, H, W) or (H, W, 6)
2) SHARD  mode: each .npy is a shard array of shape (N, 6, H, W) or (N, H, W, 6)

IMPORTANT REVISION:
- ❌ Predictions JSON/CSV output REMOVED
- ❌ No accumulation of prediction results for dumping
- ✅ Designed to be piped directly into downstream VCF generation

Key fixes for shard mode:
- Force mmap for shard arrays
- Per-worker LRU cache of opened memmaps
- Safe DataLoader settings
"""

import argparse
import os
import sys
import time
from typing import Dict, List, Tuple, Optional, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

# -----------------------------------------------------------------------------
# Model import
# -----------------------------------------------------------------------------

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(THIS_DIR, "..")))
from mynet import ConvNeXtCBAMClassifier  # noqa: E402

# -----------------------------------------------------------------------------
# Normalization (MUST MATCH TRAINING)
# -----------------------------------------------------------------------------

VAL_MEAN = torch.tensor([
    18.417816162109375, 12.649129867553711, -0.5452527403831482,
    24.723854064941406, 4.690611362457275, 0.2813551473402196
], dtype=torch.float32)

VAL_STD = torch.tensor([
    25.028322219848633, 14.809632301330566, 0.6181337833404541,
    29.972835540771484, 7.9231791496276855, 0.7659083659074717
], dtype=torch.float32)

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

def _scan_tree(start_dir: str) -> List[str]:
    stack = [start_dir]
    out = []
    while stack:
        d = stack.pop()
        try:
            with os.scandir(d) as it:
                for e in it:
                    if e.is_dir(follow_symlinks=False):
                        stack.append(e.path)
                    elif e.name.lower().endswith(".npy"):
                        out.append(os.path.abspath(e.path))
        except Exception:
            continue
    return out


def list_npy_files_parallel(root: str, workers: int) -> List[str]:
    root = os.path.abspath(root)
    top_dirs, files = [], []

    with os.scandir(root) as it:
        for e in it:
            if e.is_dir(follow_symlinks=False):
                top_dirs.append(e.path)
            elif e.name.lower().endswith(".npy"):
                files.append(os.path.abspath(e.path))

    if workers <= 1:
        for d in top_dirs:
            files.extend(_scan_tree(d))
        return sorted(files)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(_scan_tree, d) for d in top_dirs]
        for f in tqdm(as_completed(futures), total=len(futures), desc="Scan"):
            files.extend(f.result())

    return sorted(files)


def detect_input_mode(files: List[str]) -> str:
    for p in files:
        try:
            arr = np.load(p, mmap_mode="r")
            if arr.ndim == 4:
                return "shard"
            if arr.ndim == 3:
                return "single"
        except Exception:
            continue
    return "single"


def to_chw_float32(x: np.ndarray) -> np.ndarray:
    if x.shape[0] == 6:
        return x.astype(np.float32, copy=False)
    if x.shape[-1] == 6:
        return np.transpose(x, (2, 0, 1)).astype(np.float32, copy=False)
    raise ValueError(f"Unexpected shape: {x.shape}")

# -----------------------------------------------------------------------------
# Dataset
# -----------------------------------------------------------------------------

class NpyDataset(Dataset):
    def __init__(self, files: List[str], input_mode: str, cache_size: int = 2):
        self.files = files
        self.input_mode = input_mode
        self.cache_size = max(1, cache_size)

        self._cache: Dict[str, np.ndarray] = {}
        self._cache_order: List[str] = []

        self._offsets = [0]
        self._counts = []

        for f in tqdm(self.files, desc="Index"):
            try:
                arr = np.load(f, mmap_mode="r")
                n = arr.shape[0] if input_mode == "shard" else 1
            except Exception:
                n = 0
            self._counts.append(n)
            self._offsets.append(self._offsets[-1] + n)

        self.total = self._offsets[-1]
        if self.total == 0:
            raise RuntimeError("No valid samples found")

    def __len__(self):
        return self.total

    def _locate(self, idx):
        lo, hi = 0, len(self._offsets) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if self._offsets[mid + 1] <= idx:
                lo = mid + 1
            else:
                hi = mid
        return lo, idx - self._offsets[lo]

    def _get_arr(self, path: str):
        if path in self._cache:
            return self._cache[path]
        arr = np.load(path, mmap_mode="r")
        self._cache[path] = arr
        self._cache_order.append(path)
        if len(self._cache_order) > self.cache_size:
            ev = self._cache_order.pop(0)
            self._cache.pop(ev, None)
        return arr

    def __getitem__(self, idx):
        file_idx, local_idx = self._locate(idx)
        path = self.files[file_idx]
        arr = self._get_arr(path)

        if self.input_mode == "shard":
            x = arr[local_idx]
            index_in_shard = local_idx
        else:
            x = arr
            index_in_shard = 0

        x = to_chw_float32(x)
        return torch.from_numpy(x), path, index_in_shard

# -----------------------------------------------------------------------------
# Model
# -----------------------------------------------------------------------------

def load_model(ckpt_path: str, device, depths, dims):
    ckpt = torch.load(ckpt_path, map_location=device)
    genotype_map = ckpt.get("genotype_map", {})
    model = ConvNeXtCBAMClassifier(
        in_channels=6,
        class_num=len(genotype_map) or 2,
        depths=depths,
        dims=dims,
    ).to(device)
    model.load_state_dict(ckpt.get("model_state_dict", ckpt))
    model.eval()
    return model, genotype_map

# -----------------------------------------------------------------------------
# Inference (NO result accumulation)
# -----------------------------------------------------------------------------

def run_inference(model, dl, device, class_names):
    with torch.no_grad():
        for imgs, shard_paths, index_in_shard in dl:
            imgs = imgs.to(device, non_blocking=True)
            logits = model(imgs)
            probs = torch.softmax(logits, dim=1)
            top_p, top_i = probs.max(dim=1)

            for i in range(len(imgs)):
                yield {
                    "shard_path": shard_paths[i],
                    "index_in_shard": int(index_in_shard[i]),
                    "pred_class": class_names[int(top_i[i])],
                    "pred_prob": float(top_p[i]),
                }

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True)
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--batch_size", type=int, default=64)
    ap.add_argument("--num_workers", type=int, default=8)
    ap.add_argument("--depths", type=int, nargs="+", default=[3,3,27,3])
    ap.add_argument("--dims", type=int, nargs="+", default=[192,384,768,1536])
    ap.add_argument("--device", choices=["auto","cpu","cuda"], default="auto")
    args = ap.parse_args()

    device = torch.device(
        "cuda" if args.device == "auto" and torch.cuda.is_available() else "cpu"
    )

    files = list_npy_files_parallel(args.input_dir, args.num_workers)
    mode = detect_input_mode(files)
    print(f"Detected input_mode={mode}, files={len(files)}")

    ds = NpyDataset(files, input_mode=mode)
    dl = DataLoader(
        ds,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        pin_memory=(device.type=="cuda"),
        persistent_workers=True,
    )

    model, genotype_map = load_model(args.ckpt, device, args.depths, args.dims)
    class_names = [k for k,_ in sorted(genotype_map.items(), key=lambda x: x[1])] or ["false","true"]

    print("Starting inference (streaming, no JSON/CSV)...")

    n = 0
    for _ in tqdm(run_inference(model, dl, device, class_names), total=len(ds)):
        n += 1

    print(f"Finished inference: {n} samples")

if __name__ == "__main__":
    main()
