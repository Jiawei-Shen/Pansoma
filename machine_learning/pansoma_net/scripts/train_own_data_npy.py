import argparse
import torch
import torch.optim as optim
import torch.nn as nn
import sys
import os
from tqdm import tqdm
from collections import defaultdict

# Import StepLR for learning rate scheduling
from torch.optim.lr_scheduler import StepLR

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from pansoma_net import GoogLeNet  # Assuming pansoma_net.py is in the parent directory
# from dataset_npy import get_data_loader  # Assuming this is your latest data loader for .npy files
from dataset_npy_4ch import get_data_loader

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def init_weights(m):
    if isinstance(m, nn.Conv2d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)


def train_model(data_path, output_path, num_epochs=100, learning_rate=0.0001,
                batch_size=32, model_save_milestone=50,
                lr_scheduler_step_size=50, lr_scheduler_gamma=0.1):  # Added LR scheduler params
    train_loader, genotype_map = get_data_loader(
        data_dir=data_path,
        train_or_val_subdir="train",
        batch_size=batch_size,
        shuffle=True
    )
    val_loader, _ = get_data_loader(
        data_dir=data_path,
        train_or_val_subdir="val",
        batch_size=batch_size,
        shuffle=False
    )
    num_classes = len(genotype_map)

    model = GoogLeNet(num_classes=num_classes).to(device)
    model.apply(init_weights)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    # --- Initialize Learning Rate Scheduler using parameters ---
    scheduler = StepLR(optimizer, step_size=lr_scheduler_step_size, gamma=lr_scheduler_gamma)
    # -------------------------------------------------------

    print(f"Using device: {device}")
    print(f"Initial Learning Rate: {learning_rate}")
    print(f"Learning Rate will decay by a factor of {lr_scheduler_gamma} every {lr_scheduler_step_size} epochs.")
    os.makedirs(output_path, exist_ok=True)

    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        correct = 0
        total = 0

        current_epoch_lr = optimizer.param_groups[0]['lr']
        progress_bar = tqdm(train_loader, desc=f"Epoch {epoch + 1}/{num_epochs} LR: {current_epoch_lr:.1e}", leave=True)

        batch_count = 0
        for images, labels in progress_bar:
            images, labels = images.to(device), labels.to(device)

            optimizer.zero_grad()
            outputs = model(images)
            if isinstance(outputs, tuple):
                main_output, aux1, aux2 = outputs
                loss1 = criterion(main_output, labels)
                loss2 = criterion(aux1, labels)
                loss3 = criterion(aux2, labels)
                loss = loss1 + 0.3 * loss2 + 0.3 * loss3
                outputs_for_acc = main_output
            else:
                loss = criterion(outputs, labels)
                outputs_for_acc = outputs

            loss.backward()
            optimizer.step()

            running_loss += loss.item()
            batch_count += 1

            _, predicted = torch.max(outputs_for_acc, 1)
            correct += (predicted == labels).sum().item()
            total += labels.size(0)

            if total > 0:
                progress_bar.set_postfix(loss=f"{running_loss / batch_count:.4f}",
                                         acc=f"{(correct / total) * 100:.2f}%")
            else:
                progress_bar.set_postfix(loss=f"{running_loss / batch_count:.4f}")

        progress_bar.close()

        epoch_train_loss = 0
        epoch_train_acc = 0
        if batch_count > 0 and total > 0:
            epoch_train_loss = running_loss / batch_count
            epoch_train_acc = (correct / total) * 100
        else:
            print(f"Warning: No data processed in epoch {epoch + 1} training phase.")

        val_loss, val_acc, _ = evaluate_model(model, val_loader, criterion, genotype_map)

        # Update epoch summary print to just show LR at the start of desc
        print(
            f"Epoch {epoch + 1}/{num_epochs} Summary - Train Loss: {epoch_train_loss:.4f}, Train Acc: {epoch_train_acc:.2f}%, "
            f"Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.2f}% (LR for this epoch: {current_epoch_lr:.1e})")

        scheduler.step()

        if (epoch + 1) % model_save_milestone == 0:
            milestone_path = os.path.join(output_path, f"pansoma_net_epoch_{epoch + 1}.pth")
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': epoch_train_loss,
                'accuracy': epoch_train_acc,
                'val_loss': val_loss,
                'val_accuracy': val_acc,
                'genotype_map': genotype_map
            }, milestone_path)
            print(f"Milestone model saved at {milestone_path}")

    final_model_path = os.path.join(output_path, "pansoma_net_final.pth")
    torch.save({
        'epoch': num_epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict(),
        'genotype_map': genotype_map
    }, final_model_path)
    print(f"Training complete. Final model saved at {final_model_path}")


#!/usr/bin/env python3
import argparse
import torch
import torch.optim as optim
import torch.nn as nn
import sys
import os
from tqdm import tqdm
from collections import defaultdict

# Import StepLR for learning rate scheduling
from torch.optim.lr_scheduler import StepLR

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from pansoma_net import GoogLeNet  # Assuming pansoma_net.py is in the parent directory
from dataset_npy import get_data_loader  # Assuming this is your latest data loader for .npy files

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def init_weights(m):
    if isinstance(m, nn.Conv2d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)


def print_and_log(message, log_path):
    """
    Prints a message to console and appends it to the specified log file.
    """
    print(message)
    with open(log_path, 'a', encoding='utf-8') as f:
        f.write(message + '\n')


def train_model(data_path, output_path, num_epochs=100, learning_rate=0.0001,
                batch_size=32, model_save_milestone=50,
                lr_scheduler_step_size=50, lr_scheduler_gamma=0.1):
    # Ensure output directory exists and define log file
    os.makedirs(output_path, exist_ok=True)
    log_file = os.path.join(output_path, "training_log.txt")
    # Clear previous log if exists
    if os.path.exists(log_file):
        os.remove(log_file)

    print_and_log(f"Using device: {device}", log_file)
    print_and_log(f"Initial Learning Rate: {learning_rate:.1e}", log_file)
    print_and_log(
        f"Learning Rate will decay by a factor of {lr_scheduler_gamma} every {lr_scheduler_step_size} epochs.",
        log_file
    )
    print_and_log(f"Checkpoints will be saved every {model_save_milestone} epochs into: {output_path}", log_file)

    # Load data
    train_loader, genotype_map = get_data_loader(
        data_path,  # first argument
        "train",  # second positional argument: train_or_val_subdir
        batch_size=batch_size,
        shuffle=True
    )
    val_loader, _ = get_data_loader(
        data_path,
        "val",
        batch_size=batch_size,
        shuffle=False
    )

    num_classes = len(genotype_map)

    model = GoogLeNet(num_classes=num_classes).to(device)
    model.apply(init_weights)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = StepLR(optimizer, step_size=lr_scheduler_step_size, gamma=lr_scheduler_gamma)

    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        correct = 0
        total = 0

        current_lr = optimizer.param_groups[0]['lr']
        progress_bar = tqdm(
            train_loader,
            desc=f"Epoch {epoch + 1}/{num_epochs} LR: {current_lr:.1e}",
            leave=True
        )

        batch_count = 0
        for images, labels in progress_bar:
            images, labels = images.to(device), labels.to(device)

            optimizer.zero_grad()
            outputs = model(images)
            if isinstance(outputs, tuple):
                main_output, aux1, aux2 = outputs
                loss1 = criterion(main_output, labels)
                loss2 = criterion(aux1, labels)
                loss3 = criterion(aux2, labels)
                loss = loss1 + 0.3 * loss2 + 0.3 * loss3
                outputs_for_acc = main_output
            else:
                loss = criterion(outputs, labels)
                outputs_for_acc = outputs

            loss.backward()
            optimizer.step()

            running_loss += loss.item()
            batch_count += 1

            _, predicted = torch.max(outputs_for_acc, 1)
            correct += (predicted == labels).sum().item()
            total += labels.size(0)

            if total > 0:
                avg_loss = running_loss / batch_count
                avg_acc = (correct / total) * 100
                progress_bar.set_postfix(loss=f"{avg_loss:.4f}", acc=f"{avg_acc:.2f}%")
            else:
                avg_loss = running_loss / batch_count
                progress_bar.set_postfix(loss=f"{avg_loss:.4f}")

        progress_bar.close()

        # Compute epoch statistics
        if batch_count > 0 and total > 0:
            train_loss = running_loss / batch_count
            train_acc = (correct / total) * 100
        else:
            train_loss = 0.0
            train_acc = 0.0
            print_and_log(f"Warning: No data processed in epoch {epoch + 1} training phase.", log_file)

        # Validation step
        val_loss, val_acc, class_accuracy = evaluate_model(
            model, val_loader, criterion, genotype_map, log_file
        )

        # Log summary
        summary_msg = (
            f"Epoch {epoch + 1}/{num_epochs} Summary - "
            f"Train Loss: {train_loss:.4f}, Train Acc: {train_acc:.2f}%, "
            f"Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.2f}% (LR: {current_lr:.1e})"
        )
        print_and_log(summary_msg, log_file)

        # Step the scheduler
        scheduler.step()

        # Save checkpoint if milestone
        if (epoch + 1) % model_save_milestone == 0:
            milestone_path = os.path.join(output_path, f"pansoma_net_epoch_{epoch + 1}.pth")
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': train_loss,
                'accuracy': train_acc,
                'val_loss': val_loss,
                'val_accuracy': val_acc,
                'genotype_map': genotype_map
            }, milestone_path)
            print_and_log(f"Milestone model saved at: {milestone_path}", log_file)

    # Save final model
    final_model_path = os.path.join(output_path, "pansoma_net_final.pth")
    torch.save({
        'epoch': num_epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict(),
        'genotype_map': genotype_map
    }, final_model_path)
    print_and_log(f"Training complete. Final model saved at: {final_model_path}", log_file)


def evaluate_model(model, data_loader, criterion, genotype_map, log_file):
    model.eval()
    running_loss = 0.0
    correct = 0
    total = 0
    class_correct = defaultdict(int)
    class_total = defaultdict(int)
    batch_count_eval = 0

    with torch.no_grad():
        for images, labels in data_loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            if isinstance(outputs, tuple):
                outputs = outputs[0]

            loss = criterion(outputs, labels)
            running_loss += loss.item()
            batch_count_eval += 1

            _, predicted = torch.max(outputs, 1)
            correct += (predicted == labels).sum().item()
            total += labels.size(0)

            for lbl, pred in zip(labels, predicted):
                class_total[lbl.item()] += 1
                if lbl == pred:
                    class_correct[lbl.item()] += 1

    if batch_count_eval > 0 and total > 0:
        avg_loss = running_loss / batch_count_eval
        accuracy = (correct / total) * 100
    else:
        avg_loss = 0.0
        accuracy = 0.0
        print_and_log("Warning: Validation data loader is empty or yielded no samples.", log_file)

    # Log class-wise accuracy
    idx_to_class = {v: k for k, v in genotype_map.items()}
    print_and_log("\nClass-wise Validation Accuracy:", log_file)
    class_accuracy_dict = {}
    for class_idx in range(len(genotype_map)):
        class_name = idx_to_class.get(class_idx, f"Class_{class_idx}")
        total_c = class_total[class_idx]
        correct_c = class_correct[class_idx]
        acc = (correct_c / total_c) * 100 if total_c > 0 else 0.0
        class_accuracy_dict[class_name] = acc
        print_and_log(
            f"  {class_name} (Index {class_idx}): {acc:.2f}% ({correct_c}/{total_c})",
            log_file
        )
    print_and_log("", log_file)

    return avg_loss, accuracy, class_accuracy_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GoogLeNet on custom .npy dataset")
    parser.add_argument("data_path", type=str, help="Path to the dataset (containing train/val subdirectories)")
    parser.add_argument("-o", "--output_path", default="../saved_models", type=str,
                        help="Path to save the model and training log")
    parser.add_argument("--epochs", type=int, default=300, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=0.001, help="Initial learning rate")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--milestone", type=int, default=30,
                        help="Save model checkpoint every N epochs (milestone)")

    # --- New arguments for LR decay ---
    parser.add_argument("--lr_decay_epochs", type=int, default=50,
                        help="Number of epochs after which to decay learning rate")
    parser.add_argument("--lr_decay_factor", type=float, default=0.1,
                        help="Factor by which to decay learning rate")
    # ------------------------------------

    args = parser.parse_args()

    train_model(
        data_path=args.data_path,
        output_path=args.output_path,
        num_epochs=args.epochs,
        learning_rate=args.lr,
        batch_size=args.batch_size,
        model_save_milestone=args.milestone,
        lr_scheduler_step_size=args.lr_decay_epochs,
        lr_scheduler_gamma=args.lr_decay_factor
    )



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GoogLeNet on custom .npy dataset")
    parser.add_argument("data_path", type=str, help="Path to the dataset (containing train/val subdirectories)")
    parser.add_argument("-o", "--output_path", default="../saved_models", type=str, help="Path to save the model")
    parser.add_argument("--epochs", type=int, default=300, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=0.001, help="Initial learning rate")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--milestone", type=int, default=30,
                        help="Save model checkpoint every N epochs (model saving milestone)")

    # --- New arguments for LR decay ---
    parser.add_argument("--lr_decay_epochs", type=int, default=50,
                        help="Number of epochs after which to decay learning rate (default: 50)")
    parser.add_argument("--lr_decay_factor", type=float, default=0.1,
                        help="Factor by which to decay learning rate (default: 0.1)")
    # ------------------------------------

    args = parser.parse_args()

    train_model(
        data_path=args.data_path,
        output_path=args.output_path,
        num_epochs=args.epochs,
        learning_rate=args.lr,
        batch_size=args.batch_size,
        model_save_milestone=args.milestone,  # For saving model
        lr_scheduler_step_size=args.lr_decay_epochs,  # For LR scheduler
        lr_scheduler_gamma=args.lr_decay_factor  # For LR scheduler
    )