import torch
import torchvision.transforms as transforms  # Keep for type hinting if transform is passed
from torch.utils.data import Dataset
import os
import numpy as np


class NpyDataset(Dataset):
    """
    Custom PyTorch Dataset to load 4-channel .npy files.
    Assumes .npy files are loaded with shape (4, Height, Width).
    """

    def __init__(self, root_dir, transform=None):
        self.root_dir = os.path.expanduser(root_dir)
        self.transform = transform
        self.samples = []
        self.classes = []  # Added for consistency if you need it later
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
            if not os.path.isdir(class_path):  # Should not happen if filtered above
                continue

            for file_name in sorted(os.listdir(class_path)):  # Renamed 'file' to 'filename'
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
            # Load .npy file, expecting (4, Height, Width) e.g. (4, 201, 100)
            image_np = np.load(file_path)
        except Exception as e:
            raise RuntimeError(f"Failed to load .npy file {file_path}: {e}")

        if not isinstance(image_np, np.ndarray):
            raise TypeError(f"File {file_path} did not load as a NumPy array (loaded type: {type(image_np)}).")

        # Ensure the shape is (4, H, W)
        if image_np.ndim != 3 or image_np.shape[0] != 4:
            raise ValueError(f"Loaded .npy file {file_path} has unexpected shape {image_np.shape}. Expected (4, H, W).")

        # Convert CHW NumPy array to CHW PyTorch tensor.
        # NO PERMUTATION NEEDED if np.load already gives (4, H, W)
        image_tensor = torch.from_numpy(image_np.copy())  # .copy() for safety

        # The input to self.transform will be a tensor, e.g., torch.int8 with shape (4, H, W)
        if self.transform:
            image_tensor = self.transform(image_tensor)

        return image_tensor, label


def get_data_loader(data_dir, dataset_type, batch_size=32):
    """
    Load dataset from the given path using NpyDataset.
    Returns a DataLoader and a class-to-index mapping.

    :param data_dir: Path to the dataset directory (should contain "train", "val", "test")
    :param dataset_type: One of "train", "val", or "test"
    :param batch_size: Batch size for the DataLoader
    :return: DataLoader and class-to-index mapping
    """
    dataset_path = os.path.join(data_dir, dataset_type)

    if not os.path.exists(dataset_path):
        raise FileNotFoundError(f"Dataset path {dataset_path} does not exist.")

    transform = transforms.Compose([
        transforms.Normalize(mean=[0.5] * 6, std=[0.5] * 6)  # Normalize all 6 channels independently
    ])

    dataset = NpyDataset(root_dir=dataset_path, transform=transform)
    shuffle = dataset_type == "train"  # Only shuffle training data

    loader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, num_workers=4)

    return loader, dataset.class_to_idx  # Return DataLoader and class mapping


if __name__ == "__main__":
    data_root = "/path/to/organized_pileups_dataset_6channels"  # Update this path

    # Load training data
    train_loader, class_map = get_data_loader(data_root, dataset_type="train", batch_size=16)
    print(f"Loaded {len(train_loader.dataset)} training samples.")

    # Load validation data
    val_loader, _ = get_data_loader(data_root, dataset_type="val", batch_size=16)
    print(f"Loaded {len(val_loader.dataset)} validation samples.")

    # Load test data
    test_loader, _ = get_data_loader(data_root, dataset_type="test", batch_size=16)
    print(f"Loaded {len(test_loader.dataset)} test samples.")

    print(f"Class-to-index mapping: {class_map}")
