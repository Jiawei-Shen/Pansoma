import torch
import torchvision.transforms as transforms
from torch.utils.data import Dataset, DataLoader, ConcatDataset
import numpy as np
import os, glob, torch



class NpyDataset(Dataset):
    """
    --- REVISED ---
    Custom PyTorch Dataset to load 6-channel .npy files.
    Accepts .npy files with the shape (6, W, H).
    """

    def __init__(self, root_dir, transform=None, return_paths=False, resolve_symlinks=False):
        self.root_dir = os.path.expanduser(root_dir)
        self.transform = transform
        self.return_paths = return_paths
        self.samples = []
        self.classes = []
        self.class_to_idx = {}
        self.resolve_symlinks = resolve_symlinks

        if not os.path.isdir(self.root_dir):
            raise FileNotFoundError(f"Root directory not found: {self.root_dir}")

        class_folders = sorted([d.name for d in os.scandir(self.root_dir) if d.is_dir()])
        if not class_folders:
            raise FileNotFoundError(f"No class subdirectories found in {self.root_dir}")

        self.classes = class_folders
        self.class_to_idx = {cls_name: idx for idx, cls_name in enumerate(self.classes)}

        for cls_name in class_folders:
            class_path = os.path.join(self.root_dir, cls_name)
            label = self.class_to_idx[cls_name]

            for file_name in sorted(os.listdir(class_path)):
                if not file_name.lower().endswith(".npy"):
                    continue

                file_path = os.path.join(class_path, file_name)
                # Resolve real path once (avoid re-resolving every sample)
                # real_path = os.path.realpath(file_path) if self.resolve_symlinks else file_path
                #
                # if not os.path.exists(real_path):
                #     raise FileNotFoundError(f"Broken symlink: {file_path} -> {real_path}")
                # self.samples.append((real_path, label))
                self.samples.append((file_path, label))

        if len(self.samples) == 0:
            raise ValueError(f"No .npy files found in {self.root_dir}")

        print(
            f"Initialized NpyDataset from {self.root_dir}: Found {len(self.samples)} samples in {len(self.classes)} classes.")

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        file_path, label = self.samples[idx]
        try:
            # image_np = np.load(file_path)
            image_np = np.load(file_path, mmap_mode='r')
        except Exception as e:
            raise RuntimeError(f"Failed to load .npy file {file_path}: {e}")

        if not isinstance(image_np, np.ndarray):
            raise TypeError(f"File {file_path} did not load as a NumPy array (loaded type: {type(image_np)}).")

        # --- REVISED: Validate for 6 channels in (C, W, H) format ---
        # This check now assumes the input shape is (6, W, H) as requested.
        if image_np.ndim != 3 or image_np.shape[0] != 6:
            raise ValueError(
                f"Loaded .npy file {file_path} has an unsupported shape {image_np.shape}. "
                f"Expected shape is (6, W, H)."
            )

        image_tensor = torch.from_numpy(image_np.copy()).float()
        if self.transform:
            image_tensor = self.transform(image_tensor)

        if self.return_paths:
            return image_tensor, label, file_path
        else:
            return image_tensor, label


def get_data_loader(data_dir, dataset_type, batch_size=32, num_workers=16, shuffle: bool = False,
                    return_paths: bool = False):
    """
    Load dataset(s) using NpyDataset.

    Behavior:
      • If data_dir is (train_roots, val_roots):
          - Pick roots by dataset_type ("train" -> train_roots, "val" -> val_roots)
          - For EACH root, include BOTH 'train' and 'val' subfolders.
      • Else (str or list/tuple of str): back-compat mode:
          - Include only the requested `dataset_type` subfolder for each root.

    Returns: (DataLoader, class_to_idx)
    """
    import os
    from torch.utils.data import DataLoader, ConcatDataset
    from torchvision import transforms

    def _to_list(x):
        if x is None:
            return []
        if isinstance(x, (list, tuple)):
            return list(x)
        return [x]

    # Decide which roots and which subfolders to include
    if isinstance(data_dir, (list, tuple)) and len(data_dir) == 2 \
       and (isinstance(data_dir[0], (str, list, tuple)) and isinstance(data_dir[1], (str, list, tuple))):
        # New split-mode: (train_roots, val_roots)
        roots = _to_list(data_dir[0] if dataset_type == "train" else data_dir[1])
        subfolders = ["train", "val"]  # ← include BOTH for each root
    else:
        # Back-compat: single root or list of roots; only requested split
        roots = _to_list(data_dir)
        subfolders = [dataset_type]

    # Expand to concrete dataset directories
    dataset_dirs = []
    for r in roots:
        r = os.path.abspath(os.path.expanduser(r))
        for sf in subfolders:
            dataset_dirs.append(os.path.join(r, sf))

    # Existence check (strict, like original)
    missing = [p for p in dataset_dirs if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"Dataset path(s) do not exist: {missing}")

    # 6-channel normalization (unchanged placeholders)
    # transform = transforms.Compose([
    #     transforms.Normalize(
    #         mean=[18.417816162109375, 12.649129867553711, -0.5452527403831482,
    #               24.723854064941406, 4.690611362457275, 10.659969329833984],
    #         std=[25.028322219848633, 14.809632301330566, 0.6181337833404541,
    #              29.972835540771484, 7.9231791496276855, 27.151996612548828]
    #     )
    # ])

    transform = transforms.Compose([
        transforms.Normalize(
            mean=[18.417816162109375, 12.649129867553711, -0.5452527403831482,
                  24.723854064941406, 4.690611362457275, 0.2813551473402196],
            std=[25.028322219848633, 14.809632301330566, 0.6181337833404541,
                 29.972835540771484, 7.9231791496276855, 0.7659083659074717]
        )
    ])

    # Build datasets and verify consistent class_to_idx
    datasets = []
    unified_class_to_idx = None
    for d in dataset_dirs:
        ds = NpyDataset(root_dir=d, transform=transform, return_paths=return_paths)
        if unified_class_to_idx is None:
            unified_class_to_idx = ds.class_to_idx
        else:
            if ds.class_to_idx != unified_class_to_idx:
                raise ValueError("class_to_idx mismatch across dataset roots/subfolders.")
        datasets.append(ds)

    dataset = datasets[0] if len(datasets) == 1 else ConcatDataset(datasets)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle,
                        num_workers=num_workers, pin_memory=True)

    return loader, unified_class_to_idx


def get_inference_data_loader(
    data_dir,
    batch_size: int = 256,
    num_workers: int = 4,
    shuffle: bool = False,
    return_paths: bool = True,
    allow_transpose: bool = True,   # if arrays are [H,W,6], transpose to [6,H,W] (slower)
    enforce_float32: bool = True,   # cast to float32 if needed
    pin_memory: bool = True,
):
    """
    High-throughput *inference* DataLoader for unlabeled 6-channel .npy tensors.
    - Recursively loads ALL *.npy files under the given dir(s), including subfolders.
    - Accepts: str, list[str], or (train_roots, val/test_roots) — uses the second element if a pair.
    - No CPU normalization (do it on GPU for speed).
    - Uses np.load(..., mmap_mode="r") to avoid extra copies.

    Returns:
        (loader, {})  # empty class map signals 'unlabeled' mode
    """
    import os
    import glob
    import numpy as np
    import torch
    from torch.utils.data import Dataset, DataLoader

    def _to_list(x):
        if x is None:
            return []
        if isinstance(x, (list, tuple)):
            return list(x)
        return [x]

    # Resolve roots: str | list[str] | (train_roots, val/test_roots)
    if (
        isinstance(data_dir, (list, tuple))
        and len(data_dir) == 2
        and isinstance(data_dir[1], (str, list, tuple))
    ):
        roots = _to_list(data_dir[1])  # use "val/test" side for inference when a pair is given
    else:
        roots = _to_list(data_dir)

    # Collect *.npy recursively (stable order, de-duplicated)
    files = []
    for r in roots:
        r = os.path.abspath(os.path.expanduser(r))
        if not os.path.isdir(r):
            raise FileNotFoundError(f"Not a directory: {r}")
        found = glob.glob(os.path.join(r, "**", "*.npy"), recursive=True)
        files.extend(found)

    seen = set()
    files = [f for f in sorted(files) if not (f in seen or seen.add(f))]
    if not files:
        raise FileNotFoundError(f"No .npy files found under: {', '.join(map(str, roots))}")

    class _NpyUnlabeledDataset(Dataset):
        """Loads 6-channel arrays; returns (tensor, -1, path) for inference."""
        def __init__(self, paths, return_paths: bool):
            self.paths = paths
            self.return_paths = return_paths

        def __len__(self):
            return len(self.paths)

        def __getitem__(self, idx):
            path = self.paths[idx]
            # Memory-map to reduce CPU copies
            arr = np.load(path, mmap_mode="r")

            # Ensure channel-first [6,H,W]
            if arr.ndim != 3:
                raise RuntimeError(f"{os.path.basename(path)} has {arr.ndim} dims; expected 3")

            if arr.shape[0] == 6:
                pass  # already [6,H,W]
            elif allow_transpose and arr.shape[-1] == 6:
                # NOTE: transpose is slower; prefer saving arrays channel-first for speed
                arr = np.transpose(arr, (2, 0, 1))
            else:
                raise RuntimeError(
                    f"{os.path.basename(path)} shape {arr.shape}; expected [6,H,W]"
                    + (" or [H,W,6] (set allow_transpose=True)" if not allow_transpose else "")
                )

            x = torch.from_numpy(arr)
            if enforce_float32 and x.dtype != torch.float32:
                x = x.float()  # cheap cast when underlying is float32 memmap

            y = -1  # unlabeled
            if self.return_paths:
                return x, y, path
            return x, y

    ds = _NpyUnlabeledDataset(files, return_paths=return_paths)

    # IMPORTANT loader knobs for shared storage / many small files:
    # - persistent_workers=False (avoids stuck workers & reduces contention)
    # - prefetch_factor small (2) if workers > 0
    loader_kwargs = dict(
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        pin_memory=pin_memory,
        persistent_workers=False,
    )
    if num_workers and num_workers > 0:
        loader_kwargs["prefetch_factor"] = 2

    loader = DataLoader(ds, **loader_kwargs)
    return loader, {}  # {} => unlabeled/inference mode


if __name__ == "__main__":
    # --- IMPORTANT ---
    # Update this path to the root directory containing your 'train', 'val', and 'test' folders.
    data_root = "/path/to/your_6channel_dataset"
    # --- IMPORTANT ---

    batch_size = 16
    num_workers = 0

    if data_root == "/path/to/your_6channel_dataset":
        print("🛑 Please update the 'data_root' variable in the script to your dataset's actual path.")
    else:
        try:
            print(f"\n--- Loading Training Data (Batch Size: {batch_size}) ---")
            train_loader, class_map = get_data_loader(
                data_dir=data_root,
                dataset_type="train",
                batch_size=batch_size,
                num_workers=num_workers,
                shuffle=True
            )
            print(f"✅ Loaded {len(train_loader.dataset)} training samples from {os.path.join(data_root, 'train')}.")
            print(f"Class-to-index mapping: {class_map}")

            print("\n--- Checking a few training batches ---")
            for i, (data, labels) in enumerate(train_loader):
                # Now data.shape will be [batch_size, 6, H, W]
                print(f"Batch {i + 1}: Data shape: {data.shape}, Labels: {labels[:5]}...")
                if i >= 2:
                    break
            if not train_loader:
                print("Train loader is empty.")

        except (FileNotFoundError, ValueError, RuntimeError) as e:
            print(f"Error loading training data: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")