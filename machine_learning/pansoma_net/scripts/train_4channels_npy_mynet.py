#!/usr/bin/env python3
import argparse
import torch
import torch.optim as optim
import torch.nn as nn
import sys
import os
from tqdm import tqdm
from collections import defaultdict
import json

# Import StepLR for learning rate scheduling
from torch.optim.lr_scheduler import StepLR

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
# from pansoma_net import GoogLeNet
from mynet import ConvNeXtCBAMClassifier
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


def train_model(data_path, output_path, save_val_results=False, num_epochs=100, learning_rate=0.0001,
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
        log_file)
    print_and_log(f"Checkpoints will be saved every {model_save_milestone} epochs into: {output_path}", log_file)
    print_and_log(f"Using {num_workers} workers for data loading.", log_file)
    if save_val_results:
        print_and_log("Validation classification results will be saved at each milestone.", log_file)

    train_loader, genotype_map = get_data_loader(
        data_dir=data_path, dataset_type="train", batch_size=batch_size,
        num_workers=num_workers, shuffle=True
    )

    # --- MODIFIED: Request paths from the validation loader ---
    # This requires your custom Dataset to yield (image, label, path) when return_paths=True
    try:
        val_loader, _ = get_data_loader(
            data_dir=data_path, dataset_type="val", batch_size=batch_size,
            num_workers=num_workers, shuffle=False, return_paths=True
        )
    except Exception as e:
        print_and_log(f"\nFATAL: Could not create validation data loader with 'return_paths=True'.", log_file)
        print_and_log("Please ensure your 'dataset_npy_4ch.py' can handle this flag and yield file paths.", log_file)
        print_and_log(f"Error details: {e}", log_file)
        return  # Exit if we can't get paths for the required evaluation

    if not genotype_map:
        print_and_log("Error: genotype_map is empty. Check dataloader and dataset structure.", log_file)
        return
    num_classes = len(genotype_map)
    if num_classes == 0:
        print_and_log("Error: Number of classes is 0. Check dataloader.", log_file)
        return
    print_and_log(f"Number of classes: {num_classes}", log_file)
    sorted_class_names_from_map = sorted(genotype_map.keys(), key=lambda k: genotype_map[k])

    model = ConvNeXtCBAMClassifier(in_channels=4, class_num=num_classes).to(device)

    model.apply(init_weights)
    false_count = 48736
    true_count = 268

    weight = torch.tensor([1.0, false_count / true_count])
    weight = weight.to(device)
    criterion = nn.CrossEntropyLoss(weight=weight)
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = StepLR(optimizer, step_size=lr_scheduler_step_size, gamma=lr_scheduler_gamma)

    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        correct_train = 0
        total_train = 0

        current_lr = optimizer.param_groups[0]['lr']
        progress_bar = tqdm(train_loader, desc=f"Epoch {epoch + 1}/{num_epochs} LR: {current_lr:.1e}", leave=True)

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
                progress_bar.close()
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

        epoch_train_loss = (running_loss / batch_count) if batch_count > 0 else 0.0
        epoch_train_acc = (correct_train / total_train) * 100 if total_train > 0 else 0.0

        # --- MODIFIED: Unpack 4 values, including the inference results ---
        val_loss, val_acc, class_performance_stats_val, val_inference_results = evaluate_model(
            model, val_loader, criterion, genotype_map, log_file
        )

        if class_performance_stats_val:
            print_and_log("\nClass-wise Validation Accuracy:", log_file)
            for class_name in sorted_class_names_from_map:
                stats = class_performance_stats_val.get(class_name, {})
                print_and_log(
                    f"  {class_name} (Index {stats.get('idx', 'N/A')}): {stats.get('acc', 0):.2f}% ({stats.get('correct', 0)}/{stats.get('total', 0)})",
                    log_file)

        summary_msg = (
            f"Epoch {epoch + 1}/{num_epochs} Summary - Train Loss: {epoch_train_loss:.4f}, Train Acc: {epoch_train_acc:.2f}%, "
            f"Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.2f}% (LR: {current_lr:.1e})")
        print_and_log(summary_msg, log_file)

        scheduler.step()

        if (epoch + 1) % model_save_milestone == 0 or (epoch + 1) == num_epochs:
            milestone_path = os.path.join(output_path, f"model_epoch_{epoch + 1}.pth")
            torch.save({
                'epoch': epoch + 1, 'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(), 'scheduler_state_dict': scheduler.state_dict(),
                'genotype_map': genotype_map, 'in_channels': 4
            }, milestone_path)
            print_and_log(f"\nMilestone model saved at: {milestone_path}", log_file)

            # --- MODIFIED: Save the validation inference results if requested ---
            if save_val_results:
                result_path = os.path.join(output_path, f"validation_results_epoch_{epoch + 1}.json")
                try:
                    with open(result_path, 'w') as f:
                        json.dump(val_inference_results, f, indent=4)
                    print_and_log(f"Saved validation results for epoch {epoch + 1} to {result_path}", log_file)
                except Exception as e:
                    print_and_log(f"Error saving validation results: {e}", log_file)

        print_and_log("-" * 30, log_file)

    print_and_log(f"Training complete. Final model located in: {output_path}", log_file)


def evaluate_model(model, data_loader, criterion, genotype_map, log_file):
    """
    Evaluates the model and now also returns detailed classification results.
    Assumes data_loader yields (images, labels, paths).
    """
    model.eval()
    running_loss_eval = 0.0
    correct_eval = 0
    total_eval = 0
    class_correct_counts = defaultdict(int)
    class_total_counts = defaultdict(int)
    batch_count_eval = 0

    # --- MODIFIED: Added dictionary to store inference results ---
    inference_results = defaultdict(list)
    idx_to_class = {v: k for k, v in genotype_map.items()}

    if not data_loader or len(data_loader) == 0:
        return 0.0, 0.0, {}, {}

    with torch.no_grad():
        # --- MODIFIED: Expects data_loader to yield paths as the third item ---
        for images, labels, paths in data_loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            if isinstance(outputs, tuple):
                outputs = outputs[0]

            loss = criterion(outputs, labels)
            running_loss_eval += loss.item()
            batch_count_eval += 1
            _, predicted = torch.max(outputs, 1)
            correct_eval += (predicted == labels).sum().item()
            total_eval += labels.size(0)

            for i, pred_idx_tensor in enumerate(predicted):
                pred_idx = pred_idx_tensor.item()
                true_idx = labels[i].item()
                path = paths[i]

                # For accuracy stats
                class_total_counts[true_idx] += 1
                if pred_idx == true_idx:
                    class_correct_counts[true_idx] += 1

                # For inference results file
                predicted_class_name = idx_to_class[pred_idx]
                inference_results[predicted_class_name].append(os.path.basename(path))

    avg_loss_eval = (running_loss_eval / batch_count_eval) if batch_count_eval > 0 else 0.0
    overall_accuracy_eval = (correct_eval / total_eval) * 100 if total_eval > 0 else 0.0

    class_performance_stats = {}
    if genotype_map:
        for class_name, class_idx in genotype_map.items():
            correct_c = class_correct_counts[class_idx]
            total_c = class_total_counts[class_idx]
            acc_c = (correct_c / total_c) * 100 if total_c > 0 else 0.0
            class_performance_stats[class_name] = {'acc': acc_c, 'correct': correct_c, 'total': total_c,
                                                   'idx': class_idx}
    else:
        print_and_log("Warning: genotype_map is missing in evaluate_model.", log_file)

    # --- MODIFIED: Return the inference results dictionary ---
    return avg_loss_eval, overall_accuracy_eval, class_performance_stats, inference_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a Classifier on 4-channel custom .npy dataset")
    parser.add_argument("data_path", type=str, help="Path to the dataset (containing train/val subdirectories)")
    parser.add_argument("-o", "--output_path", default="./saved_models_4channel", type=str,
                        help="Path to save the model and training log")
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=0.001, help="Initial learning rate")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--num_workers", type=int, default=8,
                        help="Number of worker processes for data loading (default: 8)")
    parser.add_argument("--milestone", type=int, default=50,
                        help="Save model and run test inference every N epochs (milestone)")
    parser.add_argument("--lr_decay_epochs", type=int, default=10,
                        help="Number of epochs after which to decay learning rate")
    parser.add_argument("--lr_decay_factor", type=float, default=0.1,
                        help="Factor by which to decay learning rate")
    # --- MODIFIED: Added flag to control saving of validation results ---
    parser.add_argument("--save_val_results", action='store_true',
                        help="If set, save detailed classification results from the validation set at each milestone.")

    args = parser.parse_args()

    train_model(
        data_path=args.data_path, output_path=args.output_path,
        save_val_results=args.save_val_results,  # Pass the new argument
        num_epochs=args.epochs, learning_rate=args.lr,
        batch_size=args.batch_size, num_workers=args.num_workers,
        model_save_milestone=args.milestone,
        lr_scheduler_step_size=args.lr_decay_epochs,
        lr_scheduler_gamma=args.lr_decay_factor
    )
