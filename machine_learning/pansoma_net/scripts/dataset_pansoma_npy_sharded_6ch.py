#!/usr/bin/env python3
import os
import glob
from typing import List, Tuple, Dict, Union

import numpy as np
import torch
import torchvision.transforms as transforms
from torch.utils.data import Dataset, DataLoader, ConcatDataset, Sampler


# ─────────────────────────────────────────────────────────────
#  Label map
# ─────────────────────────────────────────────────────────────

GENOTYPE_MAP: Dict[str, int] = {
    "false": 0,
    "true": 1,
}


# ─────────────────────────────────────────────────────────────
#  Sharded NPY dataset with shard_labels + index_map
# ─────────────────────────────────────────────────────────────

class NPYShardDataset(Dataset):
    """
    Dataset for sharded .npy files written as:

        shard_*_x.npy: (N, 6, H, W)
        shard_*_y.npy: (N,)

    This version:

      • Only stores shard paths + sizes in __init__.
      • Precomputes index_map[global_idx] = (shard_idx, local_idx).
      • Lazily opens shards with np.load(..., mmap_mode="r") via a small LRU cache.
      • Records shard_labels (0 for all-false, 1 for all-true, -1 for mixed/unknown).
    """

    def __init__(
        self,
        root_dir: str,
        transform=None,
        return_paths: bool = False,
        max_cached_shards: int = 2,
    ):
        self.root_dir = os.path.abspath(os.path.expanduser(root_dir))
        self.transform = transform
        self.return_paths = return_paths
        self.max_cached_shards = max_cached_shards

        if not os.path.isdir(self.root_dir):
            raise FileNotFoundError(f"Root directory not found: {self.root_dir}")

        x_paths = sorted(glob.glob(os.path.join(self.root_dir, "*_x.npy")))
        y_paths = sorted(glob.glob(os.path.join(self.root_dir, "*_y.npy")))

        if not x_paths:
            raise ValueError(f"No *_x.npy shard files found in {self.root_dir}")
        if len(x_paths) != len(y_paths):
            raise ValueError(
                f"Shard count mismatch in {self.root_dir}: "
                f"{len(x_paths)} *_x.npy vs {len(y_paths)} *_y.npy"
            )

        self.x_paths: List[str] = x_paths
        self.y_paths: List[str] = y_paths
        self.shard_sizes: List[int] = []
        self.shard_labels: List[int] = []  # 0 (all-false), 1 (all-true), -1 (mixed/unknown)

        self.C = None
        self.H = None
        self.W = None

        total = 0
        # ── Scan each shard once to get N, shape, and label purity ──
        for x_path, y_path in zip(self.x_paths, self.y_paths):
            # x: (N, 6, H, W)
            x_arr = np.load(x_path, mmap_mode="r")
            if x_arr.ndim != 4 or x_arr.shape[1] != 6:
                raise ValueError(
                    f"{x_path}: expected x shape (N, 6, H, W), got {x_arr.shape}"
                )
            n, C, H, W = x_arr.shape
            if self.C is None:
                self.C, self.H, self.W = C, H, W
            else:
                if (C, H, W) != (self.C, self.H, self.W):
                    raise ValueError(
                        f"{x_path}: inconsistent shard shape {x_arr.shape}, "
                        f"expected (*, {self.C}, {self.H}, {self.W})"
                    )
            del x_arr

            # y: (N,)
            y_arr = np.load(y_path, mmap_mode="r")
            if y_arr.shape[0] != n:
                raise ValueError(
                    f"{y_path}: label length mismatch: {y_arr.shape[0]} vs {n}"
                )

            # Determine shard label purity: all 0, all 1, or mixed
            uniq = np.unique(y_arr)
            if uniq.size == 1 and int(uniq[0]) in (0, 1):
                shard_label = int(uniq[0])  # 0 or 1
            else:
                shard_label = -1  # mixed / unexpected
            self.shard_labels.append(shard_label)

            del y_arr

            self.shard_sizes.append(n)
            total += n

        if total == 0:
            raise ValueError(
                f"NPYShardDataset from {self.root_dir} has zero usable samples."
            )

        # ─────────────────────────────────────────────────
        # Precompute global index -> (shard_idx, local_idx)
        # ─────────────────────────────────────────────────
        self.index_map = np.empty((total, 2), dtype=np.int32)
        self.shard_offsets = np.empty(len(self.shard_sizes), dtype=np.int64)

        pos = 0
        for shard_idx, n in enumerate(self.shard_sizes):
            self.shard_offsets[shard_idx] = pos
            self.index_map[pos:pos + n, 0] = shard_idx
            self.index_map[pos:pos + n, 1] = np.arange(n, dtype=np.int32)
            pos += n

        # LRU cache: shard_idx -> (x_memmap, y_memmap)
        from collections import OrderedDict
        self._shard_cache: "OrderedDict[int, Tuple[np.ndarray, np.ndarray]]" = OrderedDict()

        n_false = sum(1 for s in self.shard_labels if s == 0)
        n_true = sum(1 for s in self.shard_labels if s == 1)
        n_mixed = sum(1 for s in self.shard_labels if s < 0)

        print(
            f"Initialized NPYShardDataset from {self.root_dir}: "
            f"{len(self.x_paths)} shards, {total} samples total. "
            f"shape per sample = (6, {self.H}, {self.W}), "
            f"max_cached_shards={self.max_cached_shards}. "
            f"Shard label summary: false={n_false}, true={n_true}, mixed/unknown={n_mixed}"
        )

    # ─────────────────────────────────────────────────
    # Basic Dataset API
    # ─────────────────────────────────────────────────
    def __len__(self) -> int:
        return self.index_map.shape[0]

    def _get_shard_arrays(self, shard_idx: int):
        """LRU-memmap loading of a shard."""
        if shard_idx in self._shard_cache:
            x_arr, y_arr = self._shard_cache.pop(shard_idx)
            self._shard_cache[shard_idx] = (x_arr, y_arr)
            return x_arr, y_arr

        x_path = self.x_paths[shard_idx]
        y_path = self.y_paths[shard_idx]

        x_arr = np.load(x_path, mmap_mode="r")  # (N, 6, H, W)
        y_arr = np.load(y_path, mmap_mode="r")  # (N,)

        self._shard_cache[shard_idx] = (x_arr, y_arr)
        if len(self._shard_cache) > self.max_cached_shards:
            old_idx, (old_x, old_y) = self._shard_cache.popitem(last=False)
            del old_x, old_y

        return x_arr, y_arr

    def __getitem__(self, idx: int):
        if idx < 0 or idx >= len(self):
            raise IndexError(f"Index {idx} out of range 0..{len(self)-1}")

        shard_idx, local_idx = self.index_map[idx]
        shard_idx = int(shard_idx)
        local_idx = int(local_idx)

        x_arr, y_arr = self._get_shard_arrays(shard_idx)

        x = x_arr[local_idx]  # (6, H, W)
        y = y_arr[local_idx]  # scalar

        # Make a contiguous, writable float32 copy
        x_np = np.array(x, dtype=np.float32, copy=True)
        x_tensor = torch.from_numpy(x_np)
        y_tensor = torch.tensor(int(y), dtype=torch.long)

        if self.transform is not None:
            x_tensor = self.transform(x_tensor)

        if self.return_paths:
            sample_id = f"{os.path.basename(self.x_paths[shard_idx])}#{local_idx}"
            return x_tensor, y_tensor, sample_id
        else:
            return x_tensor, y_tensor


# ─────────────────────────────────────────────────────────────
#  ShardWindowSampler
# ─────────────────────────────────────────────────────────────

class ShardWindowSampler(Sampler[int]):
    """
    Sample indices by *windows of shards* to:

      • Keep only a limited number of shard files “hot” at a time.
      • Enforce a class-imbalance-aware schedule over shards.

    Assumptions:
      • NPYShardDataset.shard_labels are:
          0 = all-false shard
          1 = all-true shard
         -1 = mixed/unknown (these are ignored by this sampler).
      • All shards have equal (or similar) sample counts (e.g. 4096).
    """

    def __init__(
        self,
        dataset: NPYShardDataset,
        false_shards_per_window: int,
        true_shards_per_window: int,
        shuffle_within_window: bool = True,
        seed: int = 1234,
    ):
        if not isinstance(dataset, NPYShardDataset):
            raise TypeError("ShardWindowSampler requires an NPYShardDataset.")

        self.dataset = dataset
        self.false_shards_per_window = max(1, int(false_shards_per_window))
        self.true_shards_per_window = max(1, int(true_shards_per_window))
        self.shuffle_within_window = shuffle_within_window
        self.seed = int(seed)
        self._epoch = 0

        # Split shards by label
        self.false_shards = [i for i, l in enumerate(dataset.shard_labels) if l == 0]
        self.true_shards = [i for i, l in enumerate(dataset.shard_labels) if l == 1]

        if len(self.false_shards) == 0 or len(self.true_shards) == 0:
            raise ValueError(
                "ShardWindowSampler needs at least one all-false shard and one all-true shard.\n"
                "Check that shard_labels were correctly inferred."
            )

        # Precompute windows:
        #   window k uses:
        #     • contiguous chunk of false_shards
        #     • tail of true_shards: last K shards, where K grows with k
        nF = len(self.false_shards)
        nT = len(self.true_shards)
        wF = self.false_shards_per_window
        wT = self.true_shards_per_window

        self.windows: List[Tuple[List[int], List[int]]] = []
        num_windows = (nF + wF - 1) // wF  # ceil

        for w in range(num_windows):
            f_start = w * wF
            f_end = min(f_start + wF, nF)
            if f_start >= f_end:
                break
            false_window = self.false_shards[f_start:f_end]

            # Expand tail of true shards with window index (like your example)
            k_true = min((w + 1) * wT, nT)
            true_window = self.true_shards[-k_true:] if k_true > 0 else []

            self.windows.append((false_window, true_window))

        if not self.windows:
            raise ValueError("ShardWindowSampler constructed zero windows; check configuration.")

        # Precompute length in samples (may be > len(dataset) due to positive oversampling)
        self._epoch_length = 0
        for false_window, true_window in self.windows:
            for shard_idx in list(false_window) + list(true_window):
                self._epoch_length += int(self.dataset.shard_sizes[shard_idx])

        # For convenience, keep shard_offsets as Python list
        self.shard_offsets = self.dataset.shard_offsets.tolist()

        print(
            f"[ShardWindowSampler] {len(self.windows)} windows | "
            f"false_shards_per_window={self.false_shards_per_window}, "
            f"true_shards_per_window={self.true_shards_per_window} | "
            f"epoch_length={self._epoch_length} samples"
        )

    def __len__(self) -> int:
        return self._epoch_length

    def __iter__(self):
        # Different seed each epoch for reproducible but changing shuffles
        g = torch.Generator()
        g.manual_seed(self.seed + self._epoch)
        self._epoch += 1

        for (false_window, true_window) in self.windows:
            shard_list = list(false_window) + list(true_window)
            if not shard_list:
                continue

            # Collect all global indices for shards in this window
            all_indices = []
            for shard_idx in shard_list:
                offset = self.shard_offsets[shard_idx]
                n = int(self.dataset.shard_sizes[shard_idx])
                if n <= 0:
                    continue
                shard_indices = torch.arange(offset, offset + n, dtype=torch.long)
                all_indices.append(shard_indices)

            if not all_indices:
                continue

            all_indices = torch.cat(all_indices, dim=0)

            # Shuffle within the window
            if self.shuffle_within_window:
                perm = torch.randperm(all_indices.numel(), generator=g)
                all_indices = all_indices[perm]

            # Yield indices for this window
            for idx in all_indices.tolist():
                yield int(idx)


# ─────────────────────────────────────────────────────────────
#  get_data_loader (unchanged API)
# ─────────────────────────────────────────────────────────────

def _to_list(x):
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        return list(x)
    return [x]


def get_data_loader(
    data_dir: Union[str, List[str], Tuple],
    dataset_type: str,
    batch_size: int = 32,
    num_workers: int = 16,
    shuffle: bool = False,
    return_paths: bool = False,
):
    """
    NPY-sharded version of your original get_data_loader, same interface.

    Layout per root:
        root/
          train/
            shard_train_000_x.npy
            shard_train_000_y.npy
            ...
          val/
            shard_val_000_x.npy
            shard_val_000_y.npy
            ...

    Supports:
      • data_dir = "/path/to/root"
      • data_dir = ["/root1", "/root2", ...]
      • data_dir = (train_roots, val_roots)  # split-mode
    """
    # Decide roots & subfolders
    if (
        isinstance(data_dir, (list, tuple))
        and len(data_dir) == 2
        and (isinstance(data_dir[0], (str, list, tuple)) and isinstance(data_dir[1], (str, list, tuple)))
    ):
        roots = _to_list(data_dir[0] if dataset_type == "train" else data_dir[1])
        subfolders = ["train", "val"]  # original behavior
    else:
        roots = _to_list(data_dir)
        subfolders = [dataset_type]

    dataset_dirs: List[str] = []
    for r in roots:
        r = os.path.abspath(os.path.expanduser(r))
        for sf in subfolders:
            dataset_dirs.append(os.path.join(r, sf))

    missing = [p for p in dataset_dirs if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"Dataset path(s) do not exist: {missing}")

    # Same 6-channel normalization as before
    transform = transforms.Compose([
        transforms.Normalize(
            mean=[
                18.417816162109375,
                12.649129867553711,
                -0.5452527403831482,
                24.723854064941406,
                4.690611362457275,
                0.2813551473402196,
            ],
            std=[
                25.028322219848633,
                14.809632301330566,
                0.6181337833404541,
                29.972835540771484,
                7.9231791496276855,
                0.7659083659074717,
            ],
        )
    ])

    datasets: List[Dataset] = []
    for d in dataset_dirs:
        ds = NPYShardDataset(
            root_dir=d,
            transform=transform,
            return_paths=return_paths,
            max_cached_shards=16,   # adjust if you want more/less hot shards
        )
        datasets.append(ds)

    if len(datasets) == 1:
        dataset = datasets[0]
    else:
        dataset = ConcatDataset(datasets)

    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        pin_memory=True,
    )

    return loader, GENOTYPE_MAP


# Simple smoke test
if __name__ == "__main__":
    data_root = "/path/to/your_6channel_npy_sharded_dataset"  # contains train/ and val/
    batch_size = 16
    num_workers = 0

    if data_root == "/path/to/your_6channel_npy_sharded_dataset":
        print("🛑 Please update 'data_root' to your NPY-sharded dataset root (with train/ and val/).")
    else:
        try:
            print(f"\n--- Loading Training Data (Batch Size: {batch_size}) ---")
            train_loader, class_map = get_data_loader(
                data_dir=data_root,
                dataset_type="train",
                batch_size=batch_size,
                num_workers=num_workers,
                shuffle=True,
            )
            print(f"✅ Loaded {len(train_loader.dataset)} training samples.")
            print(f"Genotype map: {class_map}")

            for i, (data, labels) in enumerate(train_loader):
                print(f"Batch {i + 1}: data.shape={data.shape}, labels[0:5]={labels[:5]}")
                if i >= 1:
                    break

        except Exception as e:
            print(f"Error loading training data: {e}")
