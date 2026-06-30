import torchvision.transforms as transforms
import torchvision.datasets as datasets
from torch.utils.data import DataLoader
import os

def get_data_loader(data_dir, train, batch_size=32):
    """
    Load dataset from the given path using ImageFolder.
    Returns a DataLoader and a class-to-index mapping.
    """
    if not os.path.exists(data_dir):
        raise FileNotFoundError(f"Dataset path {data_dir} does not exist.")

    transform = transforms.Compose([
        # transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])

    # subdir = "train" if train else "val"
    subdir = train
    dataset = datasets.ImageFolder(root=os.path.join(data_dir, subdir), transform=transform)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=train, num_workers=4)

    return loader, dataset.class_to_idx  # Return both DataLoader and class mapping

if __name__ == "__main__":
    train_loader, class_map = get_data_loader("/path/to/imagenet-1k", batch_size=16, train=True)
    print(f"Loaded {len(train_loader.dataset)} training samples.")
    print(f"Class-to-index mapping: {class_map}")
