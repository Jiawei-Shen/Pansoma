#!/usr/bin/env python3
import os, sys, json, math, argparse
import numpy as np
import h5py
from typing import List, Tuple, Dict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

# ---------- logging ----------
def log(msg): print(f"[LOG] {msg}", flush=True)

# ---------- scanning ----------
def _scan_one_class(args):
    cls_dir, label, resolve_symlinks = args
    out = []
    try:
        for fn in sorted(os.listdir(cls_dir)):
            if not fn.lower().endswith(".npy"): continue
            fp = os.path.join(cls_dir, fn)
            rp = os.path.realpath(fp) if resolve_symlinks else fp
            # if not os.path.exists(rp):
            #     raise FileNotFoundError(f"Broken symlink or missing file: {fp} -> {rp}")
            out.append((rp, label))
    except PermissionError:
        pass
    return out

def scan_split(root: str, split: str, resolve_symlinks: bool = True, scan_threads: int = 0)\
        -> Tuple[List[Tuple[str, int]], Dict[str, int]]:
    log(f"Scanning split '{split}' in {root}")
    split_dir = os.path.join(os.path.expanduser(root), split)
    if not os.path.isdir(split_dir):
        raise FileNotFoundError(f"Missing split dir: {split_dir}")

    class_names = sorted([d.name for d in os.scandir(split_dir) if d.is_dir()])
    if not class_names:
        raise FileNotFoundError(f"No class subdirectories under {split_dir}")
    class_to_idx = {name: i for i, name in enumerate(class_names)}

    tasks = []
    for cls in class_names:
        cls_dir = os.path.join(split_dir, cls)
        label = class_to_idx[cls]
        tasks.append((cls_dir, label, resolve_symlinks))

    samples = []
    if scan_threads and scan_threads > 0:
        with ThreadPoolExecutor(max_workers=scan_threads) as ex:
            for res in ex.map(_scan_one_class, tasks):
                samples.extend(res)
    else:
        for t in tasks:
            samples.extend(_scan_one_class(t))

    log(f"Completed scan for {split}: total {len(samples)} samples")
    return samples, class_to_idx

# ---------- utilities ----------
def merge_class_maps(maps: List[Dict[str, int]]) -> Dict[str, int]:
    if not maps: return {}
    base = maps[0]
    for m in maps[1:]:
        if m != base:
            raise ValueError(f"class_to_idx mismatch:\n{base}\nvs\n{m}")
    return base

def check_or_infer_shape(path: str) -> Tuple[int, int, int]:
    arr = np.load(path, mmap_mode="r")
    if arr.ndim != 3: raise ValueError(f"{os.path.basename(path)} has {arr.ndim} dims; expected 3 (C,H,W)")
    if arr.shape[0] != 6:
        if arr.shape[-1] == 6:
            raise ValueError(f"{path} is HWC with 6 channels; transpose to CHW offline")
        raise ValueError(f"{path} first dim != 6; got {arr.shape}")
    return arr.shape  # (C,H,W)

def pad_to_shape(x: np.ndarray, target_hw: Tuple[int, int]) -> np.ndarray:
    C, H, W = x.shape
    tH, tW = target_hw
    if (H, W) == (tH, tW): return x
    if H > tH or W > tW:
        raise ValueError(f"Sample {(H,W)} larger than target {(tH,tW)}; increase --pad_to")
    out = np.zeros((C, tH, tW), dtype=x.dtype)
    out[:, :H, :W] = x
    return out

def shard_ranges(n: int, shard_size: int):
    if shard_size <= 0: return [(0, n)]
    s, out = 0, []
    while s < n:
        e = min(s + shard_size, n)
        out.append((s, e))
        s = e
    return out

# ---------- shard writer (safe to run in subprocess) ----------
def write_h5_shard(out_path: str, paths: List[str], labels: np.ndarray,
                   class_to_idx: Dict[str, int], dtype: str,
                   pad_to: Tuple[int, int] = None, compression: str = "lzf",
                   chunk_rows: int = 64) -> str:
    C, H, W = check_or_infer_shape(paths[0])
    if pad_to: H, W = pad_to
    n = len(paths)
    chunk_rows = max(1, min(chunk_rows, n))
    chunks = (chunk_rows, C, H, W)

    if dtype not in ("float32", "float16"):
        raise ValueError("--dtype must be float32 or float16")
    np_dtype = np.float32 if dtype == "float32" else np.float16
    str_dt = h5py.string_dtype(encoding="utf-8")

    with h5py.File(out_path, "w") as h5:
        h5.attrs["class_to_idx"] = json.dumps(class_to_idx, ensure_ascii=False)
        h5.attrs.update(dict(channels=C, height=H, width=W, dtype=dtype, compression=compression))
        dset_x = h5.create_dataset("images", shape=(n, C, H, W), dtype=np_dtype,
                                   chunks=chunks, compression=compression, shuffle=True)
        dset_y = h5.create_dataset("labels", shape=(n,), dtype=np.int32,
                                   chunks=(min(4096, n),), compression=compression, shuffle=True)
        dset_p = h5.create_dataset("paths", shape=(n,), dtype=str_dt)

        for i, p in enumerate(paths):
            arr = np.load(p, mmap_mode="r")
            if arr.shape[0] != 6:
                raise ValueError(f"{p} not CHW-6; got {arr.shape}")
            if pad_to:
                arr = pad_to_shape(arr, (H, W))
            if arr.dtype != np_dtype:
                arr = arr.astype(np_dtype, copy=False)
            dset_x[i, ...] = arr
            dset_y[i] = int(labels[i])
            dset_p[i] = p
            if ((i+1) % 1000) == 0 or (i+1) == n:
                # return strings only (no h5py objects across processes)
                pass
    return out_path

# ---------- pack split (optionally parallel) ----------
def pack_split(out_dir: str, all_samples: List[Tuple[str, int]], class_to_idx: Dict[str, int],
               split: str, shard_size: int, dtype: str, pad_to: str,
               compression: str, chunk_rows: int, workers: int):
    if not all_samples:
        log(f"[{split}] No samples; skipping.")
        return
    paths = [p for p, _ in all_samples]
    labels = np.array([y for _, y in all_samples], dtype=np.int32)
    pad_hw = None
    if pad_to:
        tH, tW = map(int, pad_to.lower().split("x"))
        pad_hw = (tH, tW)

    os.makedirs(out_dir, exist_ok=True)
    ranges = shard_ranges(len(paths), shard_size)
    digits = int(math.ceil(math.log10(max(1, len(ranges)))))

    # prepare tasks
    tasks = []
    for shard_idx, (s, e) in enumerate(ranges):
        shard_paths = paths[s:e]
        shard_labels = labels[s:e]
        shard_tag = f"{split}.shard_{str(shard_idx).zfill(digits)}.h5"
        out_path = os.path.join(out_dir, shard_tag)
        tasks.append((out_path, shard_paths, shard_labels, class_to_idx, dtype,
                      pad_hw, compression, chunk_rows))

    log(f"[{split}] Writing {len(tasks)} shard(s) with workers={workers}")
    if workers and workers > 0:
        # process-level parallelism
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(write_h5_shard, *t): os.path.basename(t[0]) for t in tasks}
            done = 0
            for f in as_completed(futs):
                shard_name = futs[f]
                try:
                    _ = f.result()
                    done += 1
                    log(f"[{split}] Done {done}/{len(tasks)}: {shard_name}")
                except Exception as e:
                    log(f"[{split}] ❌ Shard {shard_name} failed: {e}")
                    raise
    else:
        # single-process
        for i, t in enumerate(tasks, 1):
            log(f"[{split}] (serial) {i}/{len(tasks)} -> {os.path.basename(t[0])}")
            write_h5_shard(*t)

    log(f"[{split}] Completed all shards ✅")

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Pack 6-channel .npy datasets (train/val) into HDF5 shards (parallel).")
    ap.add_argument("roots", nargs="+", help="Each root contains train/ and val/ with class folders.")
    ap.add_argument("-o", "--out_dir", required=True, help="Output directory for HDF5 files.")
    ap.add_argument("--shard_size", type=int, default=20000, help="Max samples per HDF5 file (0 = one file).")
    ap.add_argument("--dtype", choices=["float32", "float16"], default="float32")
    ap.add_argument("--pad_to", type=str, default=None, help="e.g. '128x128'. If omitted, shapes must match.")
    ap.add_argument("--compression", choices=["lzf", "gzip", "None"], default="lzf")
    ap.add_argument("--chunk_rows", type=int, default=64)
    ap.add_argument("--no_resolve_symlinks", action="store_true")
    ap.add_argument("--workers", type=int, default=0, help="Parallel shard writers (processes).")
    ap.add_argument("--scan_threads", type=int, default=0, help="Threads for directory scanning.")
    args = ap.parse_args()

    if args.compression == "None": args.compression = None
    # keep NumPy/BLAS threads from exploding per process
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    log("==== HDF5 PACK START ====")
    log(f"Roots: {args.roots}")
    log(f"Output: {args.out_dir} | workers={args.workers} | scan_threads={args.scan_threads}")

    split_samples = {"train": [], "val": []}
    split_class_maps = {"train": [], "val": []}

    for root in args.roots:
        for split in ("train", "val"):
            samples, cls_map = scan_split(root, split,
                                          resolve_symlinks=not args.no_resolve_symlinks,
                                          scan_threads=args.scan_threads)
            split_samples[split].extend(samples)
            split_class_maps[split].append(cls_map)
            log(f"[{split}] {root}: +{len(samples)} samples")

    cls_map_train = merge_class_maps(split_class_maps["train"])
    cls_map_val   = merge_class_maps(split_class_maps["val"])
    if cls_map_train != cls_map_val:
        raise ValueError(f"class_to_idx mismatch between train and val:\n{cls_map_train}\nvs\n{cls_map_val}")
    class_to_idx = cls_map_train
    log(f"class_to_idx: {class_to_idx}")

    # Stable order: by label then path
    for split in ("train", "val"):
        split_samples[split].sort(key=lambda x: (x[1], x[0]))
        log(f"[{split}] Sorted {len(split_samples[split])} samples")

    # If no pad, verify shapes are identical on a few samples
    if args.pad_to is None:
        ref_path = split_samples["train"][0][0] if split_samples["train"] else split_samples["val"][0][0]
        ref_c, ref_h, ref_w = check_or_infer_shape(ref_path)
        log(f"Reference shape: {(ref_c, ref_h, ref_w)}")

    out_dir = os.path.abspath(args.out_dir)
    log(f"Writing to {out_dir}")

    pack_split(out_dir, split_samples["train"], class_to_idx, "train",
               shard_size=args.shard_size, dtype=args.dtype, pad_to=args.pad_to,
               compression=args.compression, chunk_rows=args.chunk_rows, workers=args.workers)

    pack_split(out_dir, split_samples["val"], class_to_idx, "val",
               shard_size=args.shard_size, dtype=args.dtype, pad_to=args.pad_to,
               compression=args.compression, chunk_rows=args.chunk_rows, workers=args.workers)

    log("==== HDF5 PACK COMPLETE ✅ ====")

if __name__ == "__main__":
    main()
