#!/usr/bin/env python3
"""
Fast batch inference over a directory of .npy files (recursively), supporting:

1) SINGLE mode: each .npy is one tensor of shape (6, H, W) or (H, W, 6)
2) SHARD  mode: each .npy is a shard array of shape (N, 6, H, W) or (N, H, W, 6)

Key fixes for shard mode:
- Force mmap for shard arrays (prevents loading multi-GB shard per sample).
- Per-worker LRU cache of opened memmaps (prevents reopening files per sample).
- Safer DataLoader settings (persistent_workers, prefetch).

Optional:
- Join predictions with a variant summary NDJSON produced by your data generator
  (expects fields: shard_index, index_within_shard, plus any meta fields).
"""

import argparse
import os
import sys
import json
import time
from typing import Dict, List, Tuple, Optional, Any, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

# Import your model the same way as in training
THIS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(THIS_DIR, "..")))
from mynet import ConvNeXtCBAMClassifier  # noqa: E402

# --- EXACT SAME NORMALIZATION AS VAL (COLO829T Testing) ---
VAL_MEAN = torch.tensor([
    18.417816162109375, 12.649129867553711, -0.5452527403831482,
    24.723854064941406, 4.690611362457275, 0.2813551473402196
], dtype=torch.float32)
VAL_STD = torch.tensor([
    25.028322219848633, 14.809632301330566, 0.6181337833404541,
    29.972835540771484, 7.9231791496276855, 0.7659083659074717
], dtype=torch.float32)

# ----------------------------- Utilities ------------------------------------- #

def _scan_tree(start_dir: str) -> List[str]:
    """Iteratively scan one subtree with os.scandir (fast)."""
    stack = [start_dir]
    out: List[str] = []
    while stack:
        d = stack.pop()
        try:
            with os.scandir(d) as it:
                for e in it:
                    try:
                        if e.is_dir(follow_symlinks=False):
                            stack.append(e.path)
                        else:
                            if e.name.lower().endswith(".npy"):
                                out.append(os.path.abspath(e.path))
                    except OSError:
                        continue
        except (PermissionError, FileNotFoundError):
            continue
    return out


def _format_eta(seconds: float) -> str:
    seconds = int(max(0, seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{h:d}h{m:02d}m{s:02d}s"
    elif m > 0:
        return f"{m:d}m{s:02d}s"
    else:
        return f"{s:d}s"


def list_npy_files_parallel(root: str, workers: int) -> List[str]:
    """Parallel file discovery with tqdm bar."""
    root = os.path.abspath(os.path.expanduser(root))
    top_dirs, files = [], []

    try:
        with os.scandir(root) as it:
            for e in it:
                try:
                    if e.is_dir(follow_symlinks=False):
                        top_dirs.append(e.path)
                    elif e.name.lower().endswith(".npy"):
                        files.append(os.path.abspath(e.path))
                except OSError:
                    continue
    except FileNotFoundError:
        return []

    if workers <= 1 or not top_dirs:
        for d in top_dirs:
            files.extend(_scan_tree(d))
        return sorted(set(files))

    start = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex, \
         tqdm(total=len(top_dirs), desc="Scan", unit="dir", leave=False, dynamic_ncols=True) as bar:
        futures = [ex.submit(_scan_tree, d) for d in top_dirs]
        for fut in as_completed(futures):
            files.extend(fut.result())
            bar.update(1)
            done = bar.n
            elapsed = time.time() - start
            speed = done / max(1e-9, elapsed)
            eta = (len(top_dirs) - done) / max(1e-9, speed)
            bar.set_postfix_str(f"speed={speed:.1f} dir/s ETA={_format_eta(eta)}")

    return sorted(set(files))


def read_file_list(txt_path: str) -> List[str]:
    """Read one path per line."""
    with open(txt_path, "rb") as f:
        lines = f.read().splitlines()
    try:
        return [ln.decode("utf-8").strip() for ln in lines if ln.strip()]
    except UnicodeDecodeError:
        return [ln.strip() for ln in lines if ln.strip()]


def save_file_list(paths: List[str], txt_path: str) -> None:
    """Write one absolute path per line."""
    txt_path = os.path.abspath(os.path.expanduser(txt_path))
    os.makedirs(os.path.dirname(txt_path) or ".", exist_ok=True)
    with open(txt_path, "w") as f:
        for p in paths:
            f.write(p + "\n")
    print(f"Saved {len(paths)} paths to {txt_path}")


def print_data_info(files: List[str], mmap: bool, k: int) -> None:
    print("\n=== Data Information ===")
    print(f"Total .npy files: {len(files)} | mmap={mmap}")
    for i, p in enumerate(files[:k]):
        try:
            arr = np.load(p, mmap_mode="r" if mmap else None)
            print(f"  [{i+1}] {p} | shape={arr.shape}, dtype={arr.dtype}")
        except Exception as e:
            print(f"  [{i+1}] {p} | <error: {e}>")
    print("=" * 26)


def detect_input_mode(files: List[str], mmap: bool) -> str:
    """
    Auto-detect by inspecting the first readable file:
      - if ndim==4 -> shard
      - if ndim==3 -> single
    """
    for p in files:
        try:
            arr = np.load(p, mmap_mode="r" if mmap else None)
            if arr.ndim == 4:
                return "shard"
            if arr.ndim == 3:
                return "single"
        except Exception:
            continue
    return "single"


def to_chw_float32(x: np.ndarray, channels: int) -> np.ndarray:
    """
    Convert numpy array to float32 CHW:
      - CHW if already (C,H,W)
      - HWC if (H,W,C)
    """
    if x.ndim != 3:
        raise ValueError(f"Expected 3D tensor, got shape={x.shape}")
    if x.shape[0] == channels:
        chw = x
    elif x.shape[-1] == channels:
        chw = np.transpose(x, (2, 0, 1))
    else:
        raise ValueError(f"Unexpected shape {x.shape}; expected channels={channels}")
    return chw.astype(np.float32, copy=False)


# ------------------------- Optional meta join ------------------------------- #

def load_summary_ndjson(path: str) -> Dict[Tuple[int, int], Dict[str, Any]]:
    """
    Load NDJSON keyed by (shard_index, index_within_shard).

    WARNING: This can be large. If it's huge, prefer a streaming join:
    write predictions per-shard and stream NDJSON in shard order.
    """
    mp: Dict[Tuple[int, int], Dict[str, Any]] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue
            si = obj.get("shard_index")
            ii = obj.get("index_within_shard")
            if si is None or ii is None:
                continue
            try:
                key = (int(si), int(ii))
            except Exception:
                continue
            mp[key] = obj
    return mp


def parse_shard_index_from_filename(path: str) -> Optional[int]:
    """
    Expect filenames like shard_00012_data.npy.
    Returns 12.
    """
    base = os.path.basename(path)
    parts = base.split("_")
    # shard_00012_data.npy -> ["shard","00012","data.npy"]
    if len(parts) >= 3:
        try:
            return int(parts[1])
        except Exception:
            return None
    return None


# ---------------------------- Dataset ---------------------------------------- #

class NpyShardsDataset(Dataset):
    """
    Dataset over shards, exposing samples inside shards.

    Returns per item:
      - tensor: Float32 CHW (6,H,W)
      - shard_path: str
      - index_in_shard: int
      - shard_index: Optional[int] parsed from filename

    IMPORTANT:
      - In shard mode, we FORCE mmap to avoid loading multi-GB shard arrays per sample.
      - Each worker keeps a small LRU cache of opened memmaps.
    """
    def __init__(
        self,
        files: List[str],
        channels: int,
        mmap: bool = False,
        transform=None,
        input_mode: str = "shard",   # "single" or "shard"
        cache_size: int = 2,         # per-worker open-shard cache
    ):
        self.files = files
        self.channels = channels
        self.mmap = mmap
        self.transform = transform
        self.input_mode = input_mode
        self.cache_size = max(1, int(cache_size))

        if self.channels != 6:
            raise ValueError(f"Expected 6 channels, got {channels}")

        # Per-process cache (each DataLoader worker has its own dataset copy)
        self._arr_cache: Dict[str, np.ndarray] = {}
        self._arr_cache_order: List[str] = []

        # Build an index: global_sample_idx -> (file_idx, local_idx)
        # Use mmap while indexing so we don't load full arrays just to read shapes.
        self._offsets: List[int] = [0]
        self._counts: List[int] = []

        for fp in tqdm(self.files, desc="Index", unit="file", dynamic_ncols=True, leave=False):
            try:
                arr = np.load(fp, mmap_mode="r")
            except Exception:
                self._counts.append(0)
                self._offsets.append(self._offsets[-1])
                continue

            if self.input_mode == "single":
                self._counts.append(1 if arr.ndim == 3 else 0)
            else:
                self._counts.append(int(arr.shape[0]) if arr.ndim == 4 else 0)

            self._offsets.append(self._offsets[-1] + self._counts[-1])

        self.total_samples = self._offsets[-1]
        if self.total_samples == 0:
            raise RuntimeError("No valid samples found in provided .npy files.")

    def __len__(self):
        return self.total_samples

    def _locate(self, global_idx: int) -> Tuple[int, int]:
        # binary search in offsets
        lo, hi = 0, len(self._offsets) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if self._offsets[mid + 1] <= global_idx:
                lo = mid + 1
            else:
                hi = mid
        file_idx = lo
        local_idx = global_idx - self._offsets[file_idx]
        return file_idx, local_idx

    def _get_array(self, path: str) -> np.ndarray:
        """
        Return an ndarray for this path.

        In shard mode: ALWAYS mmap (safety).
        In single mode: honor --mmap.
        Cache a few open arrays per worker process (LRU).
        """
        if path in self._arr_cache:
            return self._arr_cache[path]

        if self.input_mode == "shard":
            arr = np.load(path, mmap_mode="r")
        else:
            arr = np.load(path, mmap_mode="r" if self.mmap else None)

        self._arr_cache[path] = arr
        self._arr_cache_order.append(path)

        if len(self._arr_cache_order) > self.cache_size:
            evict = self._arr_cache_order.pop(0)
            self._arr_cache.pop(evict, None)

        return arr

    def __getitem__(self, idx: int):
        file_idx, local_idx = self._locate(idx)
        shard_path = self.files[file_idx]
        shard_index = parse_shard_index_from_filename(shard_path)

        arr = self._get_array(shard_path)

        if self.input_mode == "single":
            x = arr
            index_in_shard = 0
        else:
            x = arr[local_idx]
            index_in_shard = int(local_idx)

        chw = to_chw_float32(x, self.channels)
        t = torch.from_numpy(chw)

        if self.transform:
            t = self.transform(t)

        return t, shard_path, index_in_shard, shard_index


# --------------------------- Model Loading ---------------------------------- #

def build_and_load_model(ckpt_path: str, device: torch.device,
                         depths: List[int], dims: List[int]):
    if not os.path.isfile(ckpt_path):
        raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")
    checkpoint = torch.load(ckpt_path, map_location=device)

    genotype_map = checkpoint.get("genotype_map", {})
    in_channels = int(checkpoint.get("in_channels", 6))

    model = ConvNeXtCBAMClassifier(
        in_channels=in_channels,
        class_num=len(genotype_map) or 2,
        depths=depths,
        dims=dims
    ).to(device)

    state = checkpoint.get("model_state_dict", checkpoint)
    model.load_state_dict(state, strict=True)
    model.eval()
    return model, genotype_map, in_channels


# --------------------------- Inference Core --------------------------------- #

def run_inference(
    model,
    dl,
    device,
    class_names: List[str],
    total_samples: int,
    no_probs: bool = False,
    summary_map: Optional[Dict[Tuple[int, int], Dict[str, Any]]] = None
):
    results = []
    processed, t0 = 0, time.time()

    with tqdm(total=total_samples, desc="Infer", unit="sample", dynamic_ncols=True, leave=True) as bar:
        with torch.no_grad():
            for images, shard_paths, index_in_shard, shard_index in dl:
                images = images.to(device, non_blocking=True)

                outputs = model(images)
                if isinstance(outputs, tuple):
                    outputs = outputs[0]

                probs = torch.softmax(outputs, dim=1)
                top_prob, top_idx = probs.max(dim=1)

                probs_cpu = None
                if not no_probs:
                    probs_cpu = probs.detach().cpu().numpy()

                bs = images.shape[0]
                for i in range(bs):
                    pred_idx = int(top_idx[i])
                    pred_name = class_names[pred_idx] if pred_idx < len(class_names) else str(pred_idx)

                    rec: Dict[str, Any] = {
                        "shard_path": str(shard_paths[i]),
                        "index_in_shard": int(index_in_shard[i]),
                        "pred_class": pred_name,
                        "pred_prob": float(top_prob[i]),
                    }

                    # Optional join with summary map
                    if summary_map is not None:
                        si_val = shard_index[i]
                        # si_val may be a Python int or torch scalar
                        if hasattr(si_val, "item"):
                            si_val = si_val.item()
                        if si_val is not None:
                            try:
                                key = (int(si_val), int(index_in_shard[i]))
                                meta = summary_map.get(key)
                                if meta:
                                    rec["meta"] = meta
                            except Exception:
                                pass

                    if not no_probs and probs_cpu is not None:
                        rec["probs"] = {class_names[j]: float(probs_cpu[i, j]) for j in range(len(class_names))}

                    results.append(rec)

                processed += bs
                bar.update(bs)
                elapsed = time.time() - t0
                speed = processed / max(1e-9, elapsed)
                eta = (total_samples - processed) / max(1e-9, speed)
                bar.set_postfix_str(f"speed={speed:.1f} samp/s ETA={_format_eta(eta)}")

    return results


def write_outputs(results, csv_path: str, json_path: str):
    if json_path:
        with open(json_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Wrote JSON predictions: {json_path}")

    if csv_path:
        import csv
        any_probs = any("probs" in r for r in results)
        any_meta = any("meta" in r for r in results)

        header = ["shard_path", "index_in_shard", "pred_class", "pred_prob"]
        class_cols: List[str] = []
        if any_probs:
            all_classes = set()
            for r in results:
                if "probs" in r:
                    all_classes.update(r["probs"].keys())
            class_cols = sorted(all_classes)
            header += class_cols

        meta_cols: List[str] = []
        if any_meta:
            preferred = [
                "node_id", "variant_key", "v_pos", "v_type", "v_ref", "v_alt",
                "alt_allele_frequency", "coverage_at_locus"
            ]
            meta_cols = preferred
            header += [f"meta.{k}" for k in meta_cols]

        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            for r in results:
                row = [r["shard_path"], r["index_in_shard"], r["pred_class"], f"{r['pred_prob']:.6f}"]
                if any_probs:
                    row.extend([f"{r.get('probs', {}).get(c, 0.0):.6f}" for c in class_cols])
                if any_meta:
                    m = r.get("meta", {}) if isinstance(r.get("meta"), dict) else {}
                    row.extend([m.get(k, "") for k in meta_cols])
                w.writerow(row)

        print(f"Wrote CSV predictions: {csv_path}")


# ----------------------------- Main ----------------------------------------- #

def main():
    p = argparse.ArgumentParser(description="Infer over .npy files (single or shard) with shard-safe mmap caching.")
    p.add_argument("--input_dir", required=False, help="Directory to scan for .npy files (recursive).")
    p.add_argument("--file_list", required=False, help="Text file with one .npy path per line (skips scanning).")
    p.add_argument("--save_file_list", required=False, help="Write discovered .npy paths to this text file.")
    p.add_argument("--ckpt", required=True)
    p.add_argument("--output_csv", default="predictions.csv")
    p.add_argument("--output_json", default="predictions.json")

    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--num_workers", type=int, default=8, help="DataLoader workers.")
    p.add_argument("--scan_workers", type=int, default=None, help="Directory scanning workers (default: num_workers).")

    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--compile", action="store_true")
    p.add_argument("--no_probs", action="store_true")

    # Still allow mmap flag; shard mode will force mmap regardless (safety).
    p.add_argument("--mmap", action="store_true", help="Use mmap for single mode; shard mode always uses mmap.")
    p.add_argument("--cache_size", type=int, default=2, help="Per-worker LRU cache size for opened shard memmaps.")

    p.add_argument("--depths", type=int, nargs="+", default=[3, 3, 27, 3])
    p.add_argument("--dims", type=int, nargs="+", default=[192, 384, 768, 1536])
    p.add_argument("--sample_info_k", type=int, default=3)

    p.add_argument("--input_mode", choices=["auto", "single", "shard"], default="auto",
                   help="auto: detect. single: each .npy is one sample. shard: each .npy contains (N,6,H,W).")
    p.add_argument("--summary_ndjson", default=None,
                   help="Optional variant_summary.ndjson to join by (shard_index, index_in_shard). "
                        "Requires shard filename like shard_00012_data.npy.")

    args = p.parse_args()

    # Device
    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    if device.type == "cuda":
        torch.backends.cudnn.benchmark = True
        torch.backends.cuda.matmul.allow_tf32 = True

    # Build model
    print(f"Loading Model from {args.ckpt}")
    print(f"depths={args.depths}, dims={args.dims}.")
    model, genotype_map, in_channels = build_and_load_model(args.ckpt, device, args.depths, args.dims)
    if in_channels != 6:
        raise ValueError(f"Expected 6-channel input, got {in_channels}")

    if args.compile:
        try:
            model = torch.compile(model)
        except Exception as e:
            print(f"torch.compile not enabled: {e}")

    # Class names
    if len(genotype_map) == 0:
        num_classes = next((m.out_features for _, m in model.named_modules()
                            if isinstance(m, torch.nn.Linear)), 2)
        class_names = [str(i) for i in range(num_classes)]
    else:
        class_names = [None] * len(genotype_map)
        for cname, cidx in genotype_map.items():
            if 0 <= cidx < len(class_names):
                class_names[cidx] = str(cname)
        class_names = [c if c else str(i) for i, c in enumerate(class_names)]

    # Resolve file list vs directory scan
    if args.file_list:
        files = read_file_list(args.file_list)
        print(f"Loaded {len(files)} paths from --file_list")
    else:
        if not args.input_dir:
            raise SystemExit("Either --file_list or --input_dir must be provided.")
        scan_workers = args.scan_workers if args.scan_workers is not None else args.num_workers
        files = list_npy_files_parallel(args.input_dir, workers=scan_workers)
        if not files:
            raise SystemExit(f"No .npy files found in {args.input_dir}")
        if args.save_file_list:
            save_file_list(files, args.save_file_list)

    print_data_info(files, args.mmap, args.sample_info_k)

    # Decide input mode
    if args.input_mode == "auto":
        mode = detect_input_mode(files, mmap=args.mmap)
        print(f"Auto-detected input_mode={mode}")
    else:
        mode = args.input_mode
        print(f"Using input_mode={mode}")

    # Optional summary join
    summary_map = None
    if args.summary_ndjson:
        print(f"Loading summary NDJSON: {args.summary_ndjson} (may be large)")
        summary_map = load_summary_ndjson(args.summary_ndjson)
        print(f"Loaded {len(summary_map)} summary records keyed by (shard_index, index_in_shard)")

    # Dataset & Loader
    print("Building Data Loader...")
    from torchvision import transforms
    norm_transform = transforms.Normalize(mean=VAL_MEAN.tolist(), std=VAL_STD.tolist())

    ds = NpyShardsDataset(
        files,
        channels=in_channels,
        mmap=args.mmap,
        transform=norm_transform,
        input_mode=mode,
        cache_size=args.cache_size,
    )

    # DataLoader knobs that reduce stalls on large memmaps
    pin_memory = (device.type == "cuda")
    persistent_workers = (args.num_workers > 0)
    prefetch_factor = 2 if args.num_workers > 0 else None

    dl = DataLoader(
        ds,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        pin_memory=pin_memory,
        shuffle=False,
        persistent_workers=persistent_workers,
        prefetch_factor=prefetch_factor,
    )

    # Run
    t0 = time.time()
    results = run_inference(
        model,
        dl,
        device,
        class_names,
        total_samples=len(ds),
        no_probs=args.no_probs,
        summary_map=summary_map,
    )
    elapsed = time.time() - t0
    print(f"\nFinished {len(ds)} samples in {elapsed:.2f}s ({len(ds)/elapsed:.2f} samp/s)")

    # Save
    write_outputs(
        results,
        csv_path=args.output_csv if args.output_csv else None,
        json_path=args.output_json if args.output_json else None,
    )


if __name__ == "__main__":
    main()
