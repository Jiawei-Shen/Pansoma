import torch
import torch.nn as nn
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from pansoma_net import GoogLeNet
from dataset import get_data_loader

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def test_model(data_path, model_path="/home/jiawei/Documents/Dockers/PansomaNet/saved_models/pansoma_net_pileup.pth", dataloader="default"):
    # Load test dataset and class mapping
    test_loader, genotype_map = get_data_loader(data_path, batch_size=32, train="test")

    # Get class names from dataset
    class_names = list(genotype_map.keys())  # Assuming genotype_map is a dictionary of class labels

    # Initialize GoogLeNet
    model = GoogLeNet(num_classes=len(class_names)).to(device)
    model.load_state_dict(torch.load(model_path))
    model.eval()

    # Overall test accuracy tracking
    correct = 0
    total = 0

    # Per-class accuracy tracking
    class_correct = {cls: 0 for cls in class_names}
    class_total = {cls: 0 for cls in class_names}

    with torch.no_grad():
        for images, labels in test_loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            _, predicted = torch.max(outputs, 1)

            total += labels.size(0)
            correct += (predicted == labels).sum().item()

            # Update per-class accuracy
            for label, pred in zip(labels.cpu().numpy(), predicted.cpu().numpy()):
                class_total[class_names[label]] += 1
                if label == pred:
                    class_correct[class_names[label]] += 1

    # Print overall accuracy
    print(f"\nOverall Test Accuracy: {100 * correct / total:.2f}%")

    # Print per-class accuracy
    print("\nPer-Class Test Accuracy:")
    for cls in class_names:
        accuracy = 100 * class_correct[cls] / class_total[cls] if class_total[cls] > 0 else 0
        print(f"  {cls}: {accuracy:.2f}%")


if __name__ == "__main__":
    test_model("/home/jiawei/Documents/Dockers/PansomaNet/data/organized_pileups_dataset")
