import torch
import torchvision.transforms as transforms
from torch.utils.data import DataLoader, Dataset
import os
import numpy as np

# from PIL import Image # Not strictly needed if only .npy and ToTensor/Normalize are used directly on numpy arrays

# Define the allowed .npy extension
NPY_EXTENSION = ".npy"


def make_dataset(directory, class_to_idx):
    """
    Scans a directory for .npy files and creates a list of samples.
    Args:
        directory (str): Root directory path.
        class_to_idx (dict): Dictionary mapping class names to class indices.
    Returns:
        list: List of (sample_path, class_idx) tuples
    """
    instances = []
    directory = os.path.expanduser(directory)
    for target_class in sorted(class_to_idx.keys()):
        class_index = class_to_idx[target_class]
        target_dir = os.path.join(directory, target_class)
        if not os.path.isdir(target_dir):
            continue
        for root, _, fnames in sorted(os.walk(target_dir, followlinks=True)):
            for fname in sorted(fnames):
                if fname.lower().endswith(NPY_EXTENSION):
                    path = os.path.join(root, fname)
                    item = path, class_index
                    instances.append(item)
    return instances


class NpyFolder(Dataset):
    """
    A custom dataset loader for .npy files structured like ImageFolder.
    Assumes .npy files contain NumPy arrays that can be converted to tensors.
    """

    def __init__(self, root, transform=None):
        """
        Args:
            root (string): Directory with all the NPY files, organized in subfolders per class.
            transform (callable, optional): Optional transform to be applied on a sample.
        """
        if not os.path.isdir(root):
            raise FileNotFoundError(f"Root directory {root} not found.")

        self.classes, self.class_to_idx = self._find_classes(root)
        if not self.classes:
            raise FileNotFoundError(f"Found no class folders in {root}.")

        self.samples = make_dataset(root, self.class_to_idx)
        if not self.samples:
            raise FileNotFoundError(f"Found no {NPY_EXTENSION} files in {root} subdirectories. "
                                    f"Please check your class folders: {self.classes}")

        self.transform = transform

    def _find_classes(self, dir_path):  # Renamed dir to dir_path to avoid confusion with built-in dir()
        """
        Finds the class folders in a dataset.
        Args:
            dir_path (string): Root directory path.
        Returns:
            tuple: (classes, class_to_idx) where classes are relative names
                for classes and class_to_idx is a dictionary.
        """
        classes = [d.name for d in os.scandir(dir_path) if d.is_dir()]
        classes.sort()
        class_to_idx = {cls_name: i for i, cls_name in enumerate(classes)}
        return classes, class_to_idx

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, index):
        """
        Args:
            index (int): Index
        Returns:
            tuple: (sample, target) where target is class_index of the target class.
        """
        path, target = self.samples[index]
        try:
            sample = np.load(path)
        except Exception as e:
            raise RuntimeError(f"Error loading .npy file {path}: {e}")

        if not isinstance(sample, np.ndarray):
            raise TypeError(f"Loaded sample from {path} is not a NumPy array, but {type(sample)}.")

        if self.transform:
            sample = self.transform(sample)

        if not isinstance(sample, torch.Tensor):
            try:
                temp_transform = transforms.ToTensor()
                sample = temp_transform(sample)
            except Exception as e:
                raise TypeError(f"Sample from {path} is of type {type(sample)} and could not be converted to a Tensor "
                                f"after transforms. Ensure ToTensor() is correctly placed or data is compatible. Error: {e}")
        return sample, target


def get_data_loader(data_dir, train_or_val_subdir, batch_size=32, shuffle=True):
    """
    Load dataset from the given path using NpyFolder.
    Returns a DataLoader and a class-to-index mapping.
    """
    full_data_path = os.path.join(data_dir, train_or_val_subdir)
    if not os.path.exists(full_data_path):
        raise FileNotFoundError(f"Dataset path {full_data_path} does not exist.")

    transform = transforms.Compose([
        transforms.ToTensor(),  # Converts np.array (possibly int8) to torch.Tensor (e.g., torch.int8)
        transforms.ConvertImageDtype(torch.float32),  # Converts to torch.float32 and scales values to [0.0, 1.0]
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])

    dataset = NpyFolder(root=full_data_path, transform=transform)

    if len(dataset) == 0:
        raise ValueError(f"No {NPY_EXTENSION} files found in {full_data_path} or its subdirectories. "
                         f"Ensure data is structured with class subfolders containing .npy files.")

    loader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, num_workers=4,
                        pin_memory=True)  # Added pin_memory for potential speedup

    return loader, dataset.class_to_idx


if __name__ == "__main__":
    # Create dummy .npy data for testing
    dummy_data_path = "./dummy_npy_dataset_test"  # Changed name slightly to avoid conflict if run multiple times
    train_path = os.path.join(dummy_data_path, "train")
    val_path = os.path.join(dummy_data_path, "val")

    os.makedirs(os.path.join(train_path, "class_false"), exist_ok=True)
    os.makedirs(os.path.join(train_path, "class_true"), exist_ok=True)
    os.makedirs(os.path.join(val_path, "class_false"), exist_ok=True)
    os.makedirs(os.path.join(val_path, "class_true"), exist_ok=True)

    # Create dummy .npy files with int8 data type to simulate the problematic scenario
    dummy_image_false_train = np.random.randint(-128, 127, size=(64, 64, 3), dtype=np.int8)
    dummy_image_true_train = np.random.randint(-128, 127, size=(64, 64, 3), dtype=np.int8)
    dummy_image_false_val = np.random.randint(-128, 127, size=(64, 64, 3), dtype=np.int8)
    dummy_image_true_val = np.random.randint(-128, 127, size=(64, 64, 3), dtype=np.int8)

    np.save(os.path.join(train_path, "class_false", "sample1.npy"), dummy_image_false_train)
    np.save(os.path.join(train_path, "class_false", "sample2.npy"), dummy_image_false_train)
    np.save(os.path.join(train_path, "class_true", "sample3.npy"), dummy_image_true_train)

    np.save(os.path.join(val_path, "class_false", "sample_v1.npy"), dummy_image_false_val)
    np.save(os.path.join(val_path, "class_true", "sample_v2.npy"), dummy_image_true_val)

    print(f"Dummy data created at {os.path.abspath(dummy_data_path)}")
    print("-" * 30)

    try:
        train_loader, class_map_train = get_data_loader(dummy_data_path,
                                                        train_or_val_subdir="train",
                                                        batch_size=2,
                                                        shuffle=True)
        print(f"Loaded {len(train_loader.dataset)} training samples.")
        print(f"Training Class-to-index mapping: {class_map_train}")
        for i, (data, target) in enumerate(train_loader):
            print(f"Train Batch {i + 1}: Data shape: {data.shape}, Data dtype: {data.dtype}, Target: {target}")
            if i >= 0:  # Print info for the first batch
                break
        print("-" * 30)

        val_loader, class_map_val = get_data_loader(dummy_data_path,
                                                    train_or_val_subdir="val",
                                                    batch_size=1,
                                                    shuffle=False)
        print(f"Loaded {len(val_loader.dataset)} validation samples.")
        print(f"Validation Class-to-index mapping: {class_map_val}")
        for i, (data, target) in enumerate(val_loader):
            print(f"Val Batch {i + 1}: Data shape: {data.shape}, Data dtype: {data.dtype}, Target: {target}")
            if i >= 0:  # Print info for the first batch
                break
        print("-" * 30)

    except Exception as e:  # Catch any exception for testing
        print(f"Error during test: {e}")
        import traceback

        traceback.print_exc()
    finally:
        import shutil

        if os.path.exists(dummy_data_path):
            shutil.rmtree(dummy_data_path)
        print(f"Cleaned up {dummy_data_path}")