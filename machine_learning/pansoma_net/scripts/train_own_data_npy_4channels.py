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
from pansoma_net import GoogLeNet  # Ensure pansoma_net.py is in the parent directory and supports in_channels
from dataset_npy_4ch import get_data_loader  # Using 4-channel dataloader

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
                batch_size=32, num_workers=4, model_save_milestone=50,
                lr_scheduler_step_size=50, lr_scheduler_gamma=0.1):
    os.makedirs(output_path, exist_ok=True)
    log_file = os.path.join(output_path, "training_log_4ch.txt")
    if os.path.exists(log_file):
        os.remove(log_file)

    print_and_log(f"Using device: {device}", log_file)
    print_and_log(f"Initial Learning Rate: {learning_rate:.1e}", log_file)
    print_and_log(
        f"Learning Rate will decay by a factor of {lr_scheduler_gamma} every {lr_scheduler_step_size} epochs.",
        log_file
    )
    print_and_log(f"Checkpoints will be saved every {model_save_milestone} epochs into: {output_path}", log_file)
    print_and_log(f"Using {num_workers} workers for data loading.", log_file)

    train_loader, genotype_map = get_data_loader(
        data_dir=data_path, dataset_type="train", batch_size=batch_size,
        num_workers=num_workers, shuffle=True
    )
    val_loader, _ = get_data_loader(
        data_dir=data_path, dataset_type="val", batch_size=batch_size,
        num_workers=num_workers, shuffle=False
    )

    if not genotype_map:
        print_and_log("Error: genotype_map is empty. Check dataloader and dataset structure.", log_file)
        return
    num_classes = len(genotype_map)
    if num_classes == 0:
        print_and_log("Error: Number of classes is 0. Check dataloader.", log_file)
        return
    print_and_log(f"Number of classes: {num_classes}", log_file)
    # Create an ordered list of class names based on their index in genotype_map
    # genotype_map is expected to be {class_name: index}
    sorted_class_names_from_map = sorted(genotype_map.keys(), key=lambda k: genotype_map[k])

    model = GoogLeNet(num_classes=num_classes, in_channels=4).to(device)
    model.apply(init_weights)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = StepLR(optimizer, step_size=lr_scheduler_step_size, gamma=lr_scheduler_gamma)

    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        correct_train = 0  # Renamed to avoid clash with 'correct' from evaluate_model
        total_train = 0  # Renamed

        current_lr = optimizer.param_groups[0]['lr']
        progress_bar = tqdm(
            train_loader,
            desc=f"Epoch {epoch + 1}/{num_epochs} LR: {current_lr:.1e}",
            leave=True  # Keep the progress bar after completion for the epoch
        )

        batch_count = 0
        for images, labels in progress_bar:
            images, labels = images.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(images)

            if isinstance(outputs, tuple) and len(outputs) == 3:
                main_output, aux1, aux2 = outputs
                loss1 = criterion(main_output, labels)
                loss2 = criterion(aux1, labels)
                loss3 = criterion(aux2, labels)
                loss = loss1 + 0.3 * loss2 + 0.3 * loss3
                outputs_for_acc = main_output
            elif isinstance(outputs, torch.Tensor):
                loss = criterion(outputs, labels)
                outputs_for_acc = outputs
            else:
                progress_bar.close()  # Close bar before raising error
                raise TypeError(f"Model output type not recognized: {type(outputs)}")

            loss.backward()
            optimizer.step()

            running_loss += loss.item()
            batch_count += 1
            _, predicted = torch.max(outputs_for_acc, 1)
            correct_train += (predicted == labels).sum().item()
            total_train += labels.size(0)

            if total_train > 0 and batch_count > 0:
                avg_loss_train = running_loss / batch_count
                avg_acc_train = (correct_train / total_train) * 100
                progress_bar.set_postfix(loss=f"{avg_loss_train:.4f}", acc=f"{avg_acc_train:.2f}%")
            elif batch_count > 0:
                avg_loss_train = running_loss / batch_count
                progress_bar.set_postfix(loss=f"{avg_loss_train:.4f}")

        # Ensure the progress bar for training is explicitly closed if not automatically handled by loop end in some environments
        # tqdm usually handles this well when the iterator is exhausted.
        # If it was set with leave=False, it would disappear. With leave=True, it stays.

        if batch_count > 0 and total_train > 0:
            epoch_train_loss = running_loss / batch_count
            epoch_train_acc = (correct_train / total_train) * 100
        else:
            epoch_train_loss = 0.0
            epoch_train_acc = 0.0
            # No need to print warning to log here, progress bar already shows lack of progress if batches = 0

        # Validation step (occurs after training progress bar has finished for the epoch)
        val_loss, val_acc, class_performance_stats_val = evaluate_model(
            model, val_loader, criterion, genotype_map, log_file  # Pass genotype_map
        )

        # --- Print Class-wise Validation Accuracy (as per user's desired order) ---
        if class_performance_stats_val:
            print_and_log("\nClass-wise Validation Accuracy:", log_file)  # Newline for spacing
            for class_name in sorted_class_names_from_map:  # Iterate in defined order
                stats = class_performance_stats_val.get(class_name)
                if stats:
                    print_and_log(
                        f"  {class_name} (Index {stats['idx']}): {stats['acc']:.2f}% ({stats['correct']}/{stats['total']})",
                        log_file)
                else:
                    # This case should ideally not happen if evaluate_model processes all classes in genotype_map
                    print_and_log(
                        f"  {class_name} (Index {genotype_map.get(class_name, 'N/A')}): Data not found in validation.",
                        log_file)

        # --- Print Epoch Summary ---
        summary_msg = (
            f"Epoch {epoch + 1}/{num_epochs} Summary - "
            f"Train Loss: {epoch_train_loss:.4f}, Train Acc: {epoch_train_acc:.2f}%, "
            f"Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.2f}% (LR: {current_lr:.1e})"
        )
        print_and_log(summary_msg, log_file)
        print_and_log("-" * 30, log_file)  # Separator after full epoch summary

        scheduler.step()

        if (epoch + 1) % model_save_milestone == 0 or (epoch + 1) == num_epochs:
            milestone_path = os.path.join(output_path, f"pansoma_net_4ch_epoch_{epoch + 1}.pth")
            torch.save({
                'epoch': epoch + 1, 'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(), 'scheduler_state_dict': scheduler.state_dict(),
                'train_loss': epoch_train_loss, 'train_accuracy': epoch_train_acc,
                'val_loss': val_loss, 'val_accuracy': val_acc,
                'genotype_map': genotype_map, 'in_channels': 4
            }, milestone_path)
            print_and_log(f"Milestone model saved at: {milestone_path}", log_file)

    final_model_path = os.path.join(output_path, "pansoma_net_4ch_final.pth")
    torch.save({
        'epoch': num_epochs, 'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(), 'scheduler_state_dict': scheduler.state_dict(),
        'genotype_map': genotype_map, 'in_channels': 4
    }, final_model_path)
    print_and_log(f"Training complete. Final model saved at: {final_model_path}", log_file)


def evaluate_model(model, data_loader, criterion, genotype_map, log_file):  # Added genotype_map
    model.eval()
    running_loss_eval = 0.0  # Renamed
    correct_eval = 0  # Renamed
    total_eval = 0  # Renamed

    # These store per-class correct counts and total counts, keyed by class index (integer)
    class_correct_counts = defaultdict(int)
    class_total_counts = defaultdict(int)
    batch_count_eval = 0

    if not data_loader or len(data_loader) == 0:
        # print_and_log("Warning: Validation data loader is empty or not provided.", log_file) # Logged by caller if needed
        return 0.0, 0.0, {}

    with torch.no_grad():
        for images, labels in data_loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            if isinstance(outputs, tuple) and len(outputs) > 0:
                outputs = outputs[0]
            elif not isinstance(outputs, torch.Tensor):
                raise TypeError(f"Model output type not recognized during eval: {type(outputs)}")

            loss = criterion(outputs, labels)
            running_loss_eval += loss.item()
            batch_count_eval += 1
            _, predicted = torch.max(outputs, 1)
            correct_eval += (predicted == labels).sum().item()
            total_eval += labels.size(0)

            for lbl_idx, pred_idx in zip(labels.cpu().numpy(), predicted.cpu().numpy()):
                class_total_counts[lbl_idx] += 1
                if lbl_idx == pred_idx:
                    class_correct_counts[lbl_idx] += 1

    avg_loss_eval = 0.0
    overall_accuracy_eval = 0.0
    if batch_count_eval > 0 and total_eval > 0:
        avg_loss_eval = running_loss_eval / batch_count_eval
        overall_accuracy_eval = (correct_eval / total_eval) * 100
    # else: # Caller can log a warning if needed
    # print_and_log("Warning: Validation data loader yielded no samples or batches.", log_file)

    # MODIFIED: Prepare detailed class performance statistics
    class_performance_stats = {}
    if genotype_map:
        # genotype_map is {class_name: index}
        # Iterate through genotype_map to ensure all known classes are reported, even if not in batch
        for class_name, class_idx in genotype_map.items():
            correct_c = class_correct_counts[class_idx]
            total_c = class_total_counts[class_idx]
            acc_c = (correct_c / total_c) * 100 if total_c > 0 else 0.0
            class_performance_stats[class_name] = {
                'acc': acc_c,
                'correct': correct_c,
                'total': total_c,
                'idx': class_idx  # Store the original index associated with the class name
            }
    else:  # Should not happen if dataloader is correct
        print_and_log("Warning: genotype_map is missing in evaluate_model. Cannot compute class names.", log_file)

    return avg_loss_eval, overall_accuracy_eval, class_performance_stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GoogLeNet on 4-channel custom .npy dataset")
    parser.add_argument("data_path", type=str, help="Path to the dataset (containing train/val subdirectories)")
    parser.add_argument("-o", "--output_path", default="./models_4channel", type=str,
                        help="Path to save the model and training log")
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=0.001, help="Initial learning rate")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--num_workers", type=int, default=4,
                        help="Number of worker processes for data loading (default: 4)")
    parser.add_argument("--milestone", type=int, default=30,
                        help="Save model checkpoint every N epochs (milestone)")
    parser.add_argument("--lr_decay_epochs", type=int, default=20,
                        help="Number of epochs after which to decay learning rate")
    parser.add_argument("--lr_decay_factor", type=float, default=0.1,
                        help="Factor by which to decay learning rate")

    args = parser.parse_args()

    train_model(
        data_path=args.data_path, output_path=args.output_path,
        num_epochs=args.epochs, learning_rate=args.lr,
        batch_size=args.batch_size, num_workers=args.num_workers,
        model_save_milestone=args.milestone,
        lr_scheduler_step_size=args.lr_decay_epochs,
        lr_scheduler_gamma=args.lr_decay_factor
    )