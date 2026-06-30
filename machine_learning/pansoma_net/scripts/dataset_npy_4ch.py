import torch
import torchvision.transforms as transforms
from torch.utils.data import Dataset, DataLoader
import os
import numpy as np


class NpyDataset(Dataset):
    """
    Custom PyTorch Dataset to load 4-channel .npy files.
    Accepts .npy files with shapes (4, H, W) or (H, W, 4).
    Files with shape (H, W, 4) will be automatically transposed to (4, H, W).
    """

    def __init__(self, root_dir, transform=None, return_paths=False):
        self.root_dir = os.path.expanduser(root_dir)
        self.transform = transform
        # --- MODIFIED: Added return_paths handling ---
        self.return_paths = return_paths
        self.samples = []
        self.classes = []
        self.class_to_idx = {}

        if not os.path.isdir(self.root_dir):
            raise FileNotFoundError(f"Root directory not found: {self.root_dir}")

        class_folders = sorted([d.name for d in os.scandir(self.root_dir) if d.is_dir()])
        if not class_folders:
            raise FileNotFoundError(f"No class subdirectories found in {self.root_dir}")

        self.classes = class_folders
        self.class_to_idx = {cls_name: idx for idx, cls_name in enumerate(self.classes)}

        for cls_name in class_folders:
            class_path = os.path.join(self.root_dir, cls_name)
            for file_name in sorted(os.listdir(class_path)):
                if file_name.lower().endswith(".npy"):
                    file_path = os.path.join(class_path, file_name)
                    self.samples.append((file_path, self.class_to_idx[cls_name]))

        if len(self.samples) == 0:
            raise ValueError(f"No .npy files found in {self.root_dir}")

        print(
            f"Initialized NpyDataset from {self.root_dir}: Found {len(self.samples)} samples in {len(self.classes)} classes.")

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        file_path, label = self.samples[idx]
        try:
            image_np = np.load(file_path)
        except Exception as e:
            raise RuntimeError(f"Failed to load .npy file {file_path}: {e}")

        if not isinstance(image_np, np.ndarray):
            raise TypeError(f"File {file_path} did not load as a NumPy array (loaded type: {type(image_np)}).")

        if image_np.ndim == 3 and image_np.shape[-1] == 4:
            image_np = image_np.transpose((2, 0, 1))

        if image_np.ndim != 3 or image_np.shape[0] != 4:
            raise ValueError(
                f"Loaded .npy file {file_path} has an unsupported shape {image_np.shape}. "
                f"Expected shape is (4, H, W) or (H, W, 4)."
            )

        image_tensor = torch.from_numpy(image_np.copy()).float()
        if self.transform:
            image_tensor = self.transform(image_tensor)

        # --- MODIFIED: Return path if requested ---
        if self.return_paths:
            return image_tensor, label, file_path
        else:
            return image_tensor, label


# --- NEW: Updated get_data_loader function from user ---
def get_data_loader(data_dir, dataset_type, batch_size=32, num_workers=8, shuffle: bool = False,
                    return_paths: bool = False):
    """
    Load dataset from the given path using NpyDataset.
    Returns a DataLoader and a class-to-index mapping.
    """
    dataset_path = os.path.join(data_dir, dataset_type)

    if not os.path.exists(dataset_path):
        raise FileNotFoundError(f"Dataset path {dataset_path} does not exist.")

    transform = transforms.Compose([
        transforms.Normalize(mean=[1.8503281380292511, -1.487870932531317, -0.5408549062596621, 47.234985618049116],
                             std=[2.4918272597486815, 2.754400126426803, 0.541133137335665, 24.340171339993212])
    ])

    # --- MODIFIED: Pass return_paths flag to the dataset ---
    dataset = NpyDataset(root_dir=dataset_path, transform=transform, return_paths=return_paths)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers, pin_memory=True)

    return loader, dataset.class_to_idx

if __name__ == "__main__":
    # --- IMPORTANT ---
    # Update this path to the root directory containing your 'train', 'val', and 'test' folders.
    data_root = "/path/to/your_4channel_dataset"
    # --- IMPORTANT ---

    batch_size = 16
    num_workers = 0

    if data_root == "/path/to/your_4channel_dataset":
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
                # Now data.shape will be [batch_size, 4, H, W]
                print(f"Batch {i + 1}: Data shape: {data.shape}, Labels: {labels[:5]}...")
                if i >= 2:
                    break
            if not train_loader:
                print("Train loader is empty.")

            print(f"\n--- Loading Validation Data (Batch Size: {batch_size}) ---")
            val_loader, _ = get_data_loader(
                data_dir=data_root,
                dataset_type="val",
                batch_size=batch_size,
                num_workers=num_workers,
                shuffle=False
            )
            print(f"✅ Loaded {len(val_loader.dataset)} validation samples from {os.path.join(data_root, 'val')}.")
            if val_loader and len(val_loader) > 0:
                val_data_sample, _ = next(iter(val_loader))
                print(f"  Sample validation data batch shape: {val_data_sample.shape}")
            else:
                print("Validation loader is empty or could not be iterated.")

            print(f"\n--- Loading Test Data (Batch Size: {batch_size}) ---")
            test_loader, _ = get_data_loader(
                data_dir=data_root,
                dataset_type="test",
                batch_size=batch_size,
                num_workers=num_workers,
                shuffle=False
            )
            print(f"✅ Loaded {len(test_loader.dataset)} test samples from {os.path.join(data_root, 'test')}.")
            if test_loader and len(test_loader) > 0:
                test_data_sample, _ = next(iter(test_loader))
                print(f"  Sample test data batch shape: {test_data_sample.shape}")
            else:
                print("Test loader is empty or could not be iterated.")

        except FileNotFoundError as e:
            print(f"Error: {e}")
        except ValueError as e:
            print(f"Error: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")