import argparse
import torch
import torch.optim as optim
import torch.nn as nn
import sys
import os
from tqdm import tqdm
from collections import defaultdict
import subprocess
from test import test_model

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from pansoma_net import GoogLeNet
from dataset import get_data_loader

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def init_weights(m):
    if isinstance(m, nn.Conv2d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)


def train_model(data_path, output_path, num_epochs=100, learning_rate=0.0001, batch_size=32, milestone=50):
    train_loader, genotype_map = get_data_loader(data_path, batch_size=batch_size, train="train")
    val_loader, _ = get_data_loader(data_path, batch_size=batch_size, train="val")
    num_classes = len(genotype_map)

    model = GoogLeNet(num_classes=num_classes).to(device)
    model.apply(init_weights)  # Apply to entire model
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    print(f"Using device: {device}")
    os.makedirs(output_path, exist_ok=True)

    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        correct = 0
        total = 0
        class_correct = defaultdict(int)
        class_total = defaultdict(int)

        progress_bar = tqdm(train_loader, desc=f"Epoch {epoch + 1}/{num_epochs}", leave=True)

        batch_count = 0
        for images, labels in progress_bar:
            images, labels = images.to(device), labels.to(device)

            optimizer.zero_grad()
            outputs = model(images)
            # If the model returns a tuple (i.e., contains auxiliary outputs), extract the main output
            if isinstance(outputs, tuple):
                outputs, aux1, aux2 = outputs  # Extracting main output and auxiliary outputs
            loss = criterion(outputs, labels)  # Compute loss only on the main output
            loss.backward()
            optimizer.step()

            running_loss += loss.item()
            batch_count += 1

            _, predicted = torch.max(outputs, 1)
            correct += (predicted == labels).sum().item()
            total += labels.size(0)

            for label, pred in zip(labels, predicted):
                class_total[label.item()] += 1
                if label == pred:
                    class_correct[label.item()] += 1

            progress_bar.set_postfix(loss=f"{running_loss / batch_count:.4f}", acc=f"{(correct / total) * 100:.2f}%")

        progress_bar.close()

        epoch_loss = running_loss / len(train_loader)
        epoch_acc = (correct / total) * 100

        # Compute validation accuracy
        val_loss, val_acc, class_val_acc = evaluate_model(model, val_loader, criterion, genotype_map)

        print(f"Epoch {epoch + 1}/{num_epochs} - Loss: {epoch_loss:.4f}, Accuracy: {epoch_acc:.2f}%, "
              f"Val Loss: {val_loss:.4f}, Val Accuracy: {val_acc:.2f}%\n")

        if (epoch + 1) % milestone == 0:
            milestone_path = os.path.join(output_path, f"pansoma_net_epoch_{epoch + 1}.pth")
            torch.save(model.state_dict(), milestone_path)
            print(f"Milestone model saved at {milestone_path}")

            # Run test_model function from test.py
            print("Running test_model on milestone model...")
            test_model(data_path, milestone_path)

    final_model_path = os.path.join(output_path, "pansoma_net_final.pth")
    torch.save(model.state_dict(), final_model_path)
    print(f"Training complete. Final model saved at {final_model_path}")


def evaluate_model(model, data_loader, criterion, genotype_map):
    model.eval()
    running_loss = 0.0
    correct = 0
    total = 0
    class_correct = defaultdict(int)
    class_total = defaultdict(int)

    with torch.no_grad():
        for images, labels in data_loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            loss = criterion(outputs, labels)

            running_loss += loss.item()
            _, predicted = torch.max(outputs, 1)
            correct += (predicted == labels).sum().item()
            total += labels.size(0)

            for label, pred in zip(labels, predicted):
                class_total[label.item()] += 1
                if label == pred:
                    class_correct[label.item()] += 1

    avg_loss = running_loss / len(data_loader)
    accuracy = (correct / total) * 100

    class_accuracy = {
        class_name: (class_correct[idx] / class_total[idx]) * 100 if class_total[idx] > 0 else 0
        for class_name, idx in genotype_map.items()
    }

    print("Class-wise Validation Accuracy:")
    for class_name, acc in class_accuracy.items():
        print(f"  {class_name}: {acc:.2f}%")

    return avg_loss, accuracy, class_accuracy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GoogLeNet on custom pileup dataset")
    parser.add_argument("data_path", type=str, help="Path to the pileup dataset")
    parser.add_argument("-o", "--output_path", default="../saved_models", type=str, help="Path to save the model")
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=0.0001, help="Learning rate")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--milestone", type=int, default=30, help="Save model every N epochs")

    args = parser.parse_args()
    train_model(args.data_path, args.output_path, args.epochs, args.lr, args.batch_size, args.milestone)
