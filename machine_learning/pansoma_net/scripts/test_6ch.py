import torch
import torch.nn as nn
import os
import sys
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from pansoma_net import GoogLeNet
from dataset_own_data_6channels import get_data_loader

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def test_model(data_path, model_path="/home/jiawei/Documents/Dockers/PansomaNet/saved_models/pansoma_net_pileup.pth",
               dataloader="default"):
    # Load test dataset and class mapping
    test_loader, genotype_map = get_data_loader(data_path, dataset_type="test", batch_size=32)

    # Get class names from dataset
    class_names = list(genotype_map.keys())  # Assuming genotype_map is a dictionary of class labels

    # Initialize GoogLeNet
    model = GoogLeNet(num_classes=len(class_names)).to(device)
    # 🔹 Load the model properly
    checkpoint = torch.load(model_path)  # Load full checkpoint
    if "model_state" in checkpoint:  # If it's a full checkpoint with metadata
        model.load_state_dict(checkpoint["model_state"])
        print(f"Loaded model from {model_path}, trained for {checkpoint.get('epoch', 'unknown')} epochs.")
    else:  # If it's just state_dict
        model.load_state_dict(checkpoint)
    model.eval()

    # Overall test accuracy tracking
    correct = 0
    total = 0

    # Per-class accuracy, precision, recall, and F1 tracking
    class_correct = {cls: 0 for cls in class_names}
    class_total = {cls: 0 for cls in class_names}

    # Additional metrics: TP, FP, FN for recall and F1-score
    TP = {cls: 0 for cls in class_names}  # True Positives
    FP = {cls: 0 for cls in class_names}  # False Positives
    FN = {cls: 0 for cls in class_names}  # False Negatives

    with torch.no_grad():
        for images, labels in test_loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            _, predicted = torch.max(outputs, 1)

            total += labels.size(0)
            correct += (predicted == labels).sum().item()

            for label, pred in zip(labels.cpu().numpy(), predicted.cpu().numpy()):
                class_total[class_names[label]] += 1
                if label == pred:
                    class_correct[class_names[label]] += 1
                    TP[class_names[label]] += 1  # True Positive
                else:
                    FP[class_names[pred]] += 1  # False Positive
                    FN[class_names[label]] += 1  # False Negative

    # Print overall accuracy
    overall_accuracy = 100 * correct / total
    print(f"\nOverall Test Accuracy: {overall_accuracy:.2f}%")

    # Compute and print per-class metrics
    print("\nPer-Class Performance Metrics:")
    precision_values = []
    recall_values = []
    f1_values = []
    weighted_f1_sum = 0

    for cls in class_names:
        precision = TP[cls] / (TP[cls] + FP[cls]) if (TP[cls] + FP[cls]) > 0 else 0
        recall = TP[cls] / (TP[cls] + FN[cls]) if (TP[cls] + FN[cls]) > 0 else 0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

        precision_values.append(precision)
        recall_values.append(recall)
        f1_values.append(f1_score)
        weighted_f1_sum += f1_score * class_total[cls]

        print(f"  {cls}:")
        print(f"    Accuracy:  {100 * class_correct[cls] / class_total[cls]:.2f}%" if class_total[
                                                                                          cls] > 0 else "    Accuracy: N/A")
        print(f"    Precision: {precision:.2f}")
        print(f"    Recall:    {recall:.2f}")
        print(f"    F1-score:  {f1_score:.2f}")

    # Compute macro and weighted F1-score
    macro_f1 = sum(f1_values) / len(f1_values) if len(f1_values) > 0 else 0
    weighted_f1 = weighted_f1_sum / total if total > 0 else 0

    print("\nOverall F1 Scores:")
    print(f"  Macro F1-score: {macro_f1:.2f}")
    print(f"  Weighted F1-score: {weighted_f1:.2f}")


if __name__ == "__main__":
    test_model("/home/jiawei/Documents/Dockers/PansomaNet/data/organized_pileups_dataset")
