import torch
import torchvision.transforms as transforms
from PIL import Image
from pansoma_net import GoogLeNet
import os

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def predict(image_path, class_labels):
    model = GoogLeNet(num_classes=len(class_labels)).to(device)
    model.load_state_dict(torch.load("../saved_models/pansoma_net_imagenet.pth"))
    model.eval()

    transform = transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])

    image = Image.open(image_path).convert("RGB")
    image = transform(image).unsqueeze(0).to(device)

    with torch.no_grad():
        output = model(image)
        _, predicted = torch.max(output, 1)

    print(f"Predicted Class: {class_labels[predicted.item()]}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python inference.py <image_path> <data_path>")
    else:
        from dataset import get_data_loader
        loader = get_data_loader(sys.argv[2], batch_size=1, train=True)
        class_labels = loader.dataset.classes
        predict(sys.argv[1], class_labels)
