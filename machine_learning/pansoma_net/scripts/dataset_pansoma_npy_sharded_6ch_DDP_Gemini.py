import torch
import torchvision.transforms as transforms
from torch.utils.data import Dataset, DataLoader, ConcatDataset
import numpy as np
import os
import glob


class NpyDataset(Dataset):
    def __init__(self, root_dir, transform=None, return_paths=False, load_to_ram=False):
        self.root_path = os.path.abspath(os.path.expanduser(root_dir))
        self.transform = transform
        self.return_paths = return_paths
        self.load_to_ram = load_to_ram  # <--- New Flag
        self.samples = []
        self.mode = None
        self.data_in_ram = None  # Will hold the numpy array if load_to_ram is True

        # 1. Check if input is a specific file (Single Shard Mode)
        if os.path.isfile(self.root_path) and self.root_path.endswith("_data.npy"):
            self.mode = 'sharded'
            self._init_sharded([self.root_path])

        # 2. Check if input is a directory
        elif os.path.isdir(self.root_path):
            shard_files = sorted(glob.glob(os.path.join(self.root_path, "*_data.npy")))
            if len(shard_files) > 0:
                self.mode = 'sharded'
                self._init_sharded(shard_files)
            else:
                self.mode = 'folder'
                self._init_folder()
        else:
            raise FileNotFoundError(f"Path not found: {self.root_path}")

        if len(self.samples) == 0:
            raise ValueError(f"No data found in {self.root_path}")

        # Infer classes (only needed for map generation)
        unique_labels = sorted(list(set(s[2] for s in self.samples)))
        self.class_to_idx = {str(l): l for l in unique_labels}

        # --- RAM LOADING LOGIC ---
        # We only support RAM loading if initialized with a SINGLE shard file
        # (to prevent accidental explosion if a directory is passed)
        if self.load_to_ram and self.mode == 'sharded':
            if len(set(s[0] for s in self.samples)) == 1:
                data_path = self.samples[0][0]
                # print(f"   -> Loading {os.path.basename(data_path)} into RAM...")
                try:
                    self.data_in_ram = np.load(data_path)  # Full read
                except Exception as e:
                    raise RuntimeError(f"Failed to load {data_path} into RAM: {e}")
            else:
                print("Warning: load_to_ram=True ignored because multiple files are present in this Dataset instance.")

    def _init_sharded(self, shard_files):
        for data_path in shard_files:
            label_path = data_path.replace("_data.npy", "_labels.npy")
            if not os.path.exists(label_path):
                # We raise error here to be caught by the trainer's try/except block
                raise FileNotFoundError(f"Missing label file for {data_path}")

            labels = np.load(label_path)
            for local_idx, label in enumerate(labels):
                self.samples.append((data_path, local_idx, int(label)))

    def _init_folder(self):
        class_folders = sorted([d.name for d in os.scandir(self.root_path) if d.is_dir()])
        self.class_to_idx = {cls_name: idx for idx, cls_name in enumerate(class_folders)}
        for cls_name in class_folders:
            class_path = os.path.join(self.root_path, cls_name)
            label = self.class_to_idx[cls_name]
            for file_name in sorted(os.listdir(class_path)):
                if file_name.lower().endswith(".npy"):
                    file_path = os.path.join(class_path, file_name)
                    self.samples.append((file_path, -1, label))

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        path, local_idx, label = self.samples[idx]

        try:
            if self.mode == 'sharded':
                if self.data_in_ram is not None:
                    # --- RAM ACCESS (Fastest) ---
                    image_np = self.data_in_ram[local_idx].copy()
                else:
                    # --- DISK ACCESS (Slower) ---
                    mmap_arr = np.load(path, mmap_mode='r')
                    image_np = mmap_arr[local_idx].copy()
            else:
                image_np = np.load(path)
        except Exception as e:
            raise RuntimeError(f"Failed to load {path} at idx {local_idx}: {e}")

        if image_np.ndim != 3 or image_np.shape[0] != 6:
            raise ValueError(f"Shape mismatch: {image_np.shape}. Expected (6, H, W).")

        image_tensor = torch.from_numpy(image_np).float()

        if self.transform:
            image_tensor = self.transform(image_tensor)

        if self.return_paths:
            return image_tensor, label, path
        else:
            return image_tensor, label


# Helper for Validation loader (Standard Concatenation)
def get_data_loader(data_dir, dataset_type, batch_size=32, num_workers=16, shuffle=False, return_paths=False):
    def _to_list(x):
        if x is None: return []
        if isinstance(x, (list, tuple)): return list(x)
        return [x]

    if isinstance(data_dir, (list, tuple)) and len(data_dir) == 2 and isinstance(data_dir[0], (list, tuple, str)):
        roots = _to_list(data_dir[0] if dataset_type == "train" else data_dir[1])
        subfolders = ["train", "val"]
    else:
        roots = _to_list(data_dir)
        subfolders = [dataset_type]

    final_dirs = []
    for r in roots:
        r = os.path.abspath(os.path.expanduser(r))
        if len(glob.glob(os.path.join(r, "*_data.npy"))) > 0:
            final_dirs.append(r)
        else:
            found = False
            for sf in subfolders:
                t = os.path.join(r, sf)
                if os.path.exists(t): final_dirs.append(t); found = True
            if not found and os.path.exists(r): final_dirs.append(r)

    if not final_dirs: raise FileNotFoundError(f"No '{dataset_type}' folders found.")

    transform = transforms.Compose([
        transforms.Normalize(mean=[18.417, 12.649, -0.545, 24.723, 4.690, 0.281],
                             std=[25.028, 14.809, 0.618, 29.972, 7.923, 0.765])
    ])

    datasets = []
    unified_map = None
    for d in final_dirs:
        try:
            ds = NpyDataset(root_dir=d, transform=transform, return_paths=return_paths)
            if unified_map is None: unified_map = ds.class_to_idx
            datasets.append(ds)
        except:
            continue

    if not datasets: raise RuntimeError("No datasets created.")
    final_ds = ConcatDataset(datasets)
    return DataLoader(final_ds, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers,
                      pin_memory=True), unified_map