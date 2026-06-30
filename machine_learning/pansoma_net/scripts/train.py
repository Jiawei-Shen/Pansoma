import argparse
import torch
import torch.optim as optim
import torch.nn as nn
import sys
import os

# Add the root directory to Python's module search path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from pansoma_net import GoogLeNet  # Now it should work
from dataset import get_data_loader
import os

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def train_model(data_path, output_path, num_epochs=10, learning_rate=0.001, batch_size=64):
    train_loader = get_data_loader(data_path, batch_size=batch_size, train=True)
    model = GoogLeNet(num_classes=len(train_loader.dataset.classes)).to(device)

    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    print(f"Using device: {device}")
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for images, labels in train_loader:
            images, labels = images.to(device), labels.to(device)

            optimizer.zero_grad()
            outputs = model(images)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            running_loss += loss.item()

        print(f"Epoch {epoch + 1}/{num_epochs}, Loss: {running_loss / len(train_loader)}")

    os.makedirs(f"{output_path}", exist_ok=True)
    model_path = f"{output_path}/pansoma_net_imagenet.pth"
    torch.save(model.state_dict(), model_path)
    print(f"Training complete. Model saved at {model_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GoogLeNet on ImageNet-1K")
    parser.add_argument("data_path", type=str, help="Path to the ImageNet-1K dataset")
    parser.add_argument("-o", "--output_path", default="../saved_models", type=str, help="Path to the output model")
    parser.add_argument("--epochs", type=int, default=10, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=0.001, help="Learning rate")
    parser.add_argument("--batch_size", type=int, default=64, help="Batch size")

    args = parser.parse_args()
    train_model(args.data_path, output_path=args.output_path, num_epochs=args.epochs, learning_rate=args.lr, batch_size=args.batch_size)
