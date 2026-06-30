#!/usr/bin/env python3
import argparse
import torch
import torch.optim as optim
import torch.nn as nn
import sys
import os
from tqdm import tqdm
from collections import defaultdict
import torch.nn.functional as F
from torch.utils.data import Subset, DataLoader
import json

# --- MODIFIED: Import new schedulers ---
from torch.optim.lr_scheduler import SequentialLR, LinearLR, CosineAnnealingLR

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from mynet import ConvNeXtCBAMClassifier
from dataset_pansoma_npy_6ch import get_data_loader

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class MultiClassFocalLoss(nn.Module):
    """
    Focal Loss for multi-class classification.
    """
    def __init__(self, gamma=2.0, weight=None, reduction='mean'):
        super(MultiClassFocalLoss, self).__init__()
        self.gamma = gamma
        self.weight = weight
        self.reduction = reduction

    def forward(self, logits, targets):
        log_probs = F.log_softmax(logits, dim=1)
        log_pt = log_probs.gather(1, targets.view(-1, 1)).squeeze(1)
        pt = torch.exp(log_pt)

        if self.weight is not None:
            at = self.weight.gather(0, targets)
            log_pt = log_pt * at

        focal_loss = -1 * (1 - pt)**self.gamma * log_pt

        if self.reduction == 'mean':
            return focal_loss.mean()
        elif self.reduction == 'sum':
            return focal_loss.sum()
        return focal_loss


class CombinedFocalWeightedCELoss(nn.Module):
    def __init__(self, initial_lr, pos_weight=None, gamma=2.0):
        super().__init__()
        self.initial_lr = initial_lr
        self.focal_loss = MultiClassFocalLoss(gamma=gamma, weight=pos_weight)
        self.wce_loss = nn.CrossEntropyLoss(weight=pos_weight)

    def forward(self, logits, targets, current_lr):
        focal_weight = 1.0 - (current_lr / self.initial_lr)
        wce_weight = 1.0 - focal_weight
        loss_focal = self.focal_loss(logits, targets)
        loss_wce = self.wce_loss(logits, targets)
        return focal_weight * loss_focal + wce_weight * loss_wce


def init_weights(m):
    if isinstance(m, nn.Conv2d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)


def print_and_log(message, log_path):
    print(message)
    with open(log_path, 'a', encoding='utf-8') as f:
        f.write(message + '\n')


def train_model(data_path, output_path, save_val_results=False, num_epochs=100, learning_rate=0.0001,
                batch_size=32, num_workers=4, loss_type='weighted_ce',
                warmup_epochs=10, weight_decay=0.05, depths=None, dims=None,
                training_data_ratio=1.0, pos_weight=88.0):  # <--- NEW: pos_weight with default 88.0
    os.makedirs(output_path, exist_ok=True)
    log_file = os.path.join(output_path, "training_log_6ch.txt")
    if os.path.exists(log_file):
        os.remove(log_file)

    MIN_SAVE_EPOCH = 5  # save first checkpoint after 5 epochs

    if not (0 < training_data_ratio <= 1.0):
        raise ValueError(f"--training_data_ratio must be in (0,1], got {training_data_ratio}")

    print_and_log(f"Using device: {device}", log_file)
    print_and_log(f"Initial Learning Rate: {learning_rate:.1e}", log_file)
    print_and_log(f"Using Cosine Annealing scheduler with a {warmup_epochs}-epoch linear warmup.", log_file)
    print_and_log(
        f"Checkpointing: snapshot at epoch {MIN_SAVE_EPOCH}, then save only when validation improves (best so far).",
        log_file)
    print_and_log(f"Using {num_workers} workers for data loading.", log_file)
    if save_val_results:
        print_and_log("Will save validation results when a new best is found.", log_file)

    # Build loaders
    train_loader, genotype_map = get_data_loader(
        data_dir=data_path, dataset_type="train", batch_size=batch_size,
        num_workers=num_workers, shuffle=True
    )

    # Randomly subsample training data if requested
    if training_data_ratio < 1.0:
        full_ds = train_loader.dataset
        n = len(full_ds)
        k = max(1, int(round(n * training_data_ratio)))
        idx = torch.randperm(n)[:k].tolist()
        subset = Subset(full_ds, idx)
        train_loader = DataLoader(subset, batch_size=batch_size, shuffle=True,
                                  num_workers=num_workers, pin_memory=True)
        print_and_log(f"Training subset: using {k}/{n} samples (~{training_data_ratio:.2f} of data).", log_file)

    try:
        val_loader, _ = get_data_loader(
            data_dir=data_path, dataset_type="val", batch_size=batch_size,
            num_workers=num_workers, shuffle=False, return_paths=True
        )
    except Exception as e:
        print_and_log(f"\nFATAL: Could not create validation data loader with 'return_paths=True'.", log_file)
        print_and_log("Please ensure your 'dataset_pansoma_npy_6ch.py' can handle this flag.", log_file)
        print_and_log(f"Error details: {e}", log_file)
        return

    if not genotype_map:
        print_and_log("Error: genotype_map is empty. Check dataloader.", log_file)
        return
    num_classes = len(genotype_map)
    if num_classes == 0:
        print_and_log("Error: Number of classes is 0. Check dataloader.", log_file)
        return
    print_and_log(f"Number of classes: {num_classes}", log_file)
    sorted_class_names_from_map = sorted(genotype_map.keys(), key=lambda k: genotype_map[k])

    model = ConvNeXtCBAMClassifier(in_channels=6, class_num=num_classes,
                                   depths=depths, dims=dims).to(device)

    model.apply(init_weights)
    # --- REVISED: simple pos_weight (default 88.0), no min/false/true counts ---
    pos_weight_value = float(pos_weight)
    class_weights = torch.tensor([1.0, pos_weight_value]).to(device)

    if loss_type == "combined":
        criterion = CombinedFocalWeightedCELoss(initial_lr=learning_rate, pos_weight=class_weights)
        print_and_log(f"Using Combined Focal Loss and Weighted CE Loss.", log_file)
    elif loss_type == "weighted_ce":
        criterion = nn.CrossEntropyLoss(weight=class_weights)
        print_and_log(f"Using Weighted CE Loss.", log_file)
    else:
        raise ValueError(f"Unsupported loss_type: {loss_type}")

    optimizer = optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    print_and_log(f"Using AdamW optimizer with weight decay: {weight_decay}", log_file)

    warmup_scheduler = LinearLR(optimizer, start_factor=0.01, total_iters=warmup_epochs)
    main_scheduler = CosineAnnealingLR(optimizer, T_max=num_epochs - warmup_epochs, eta_min=0)
    scheduler = SequentialLR(optimizer, schedulers=[warmup_scheduler, main_scheduler], milestones=[warmup_epochs])

    best_val_acc = float("-inf")
    best_val_loss = float("inf")
    best_epoch = 0

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

            if loss_type == "combined":
                if isinstance(outputs, tuple) and len(outputs) == 3:
                    main_output, aux1, aux2 = outputs
                    loss1 = criterion(main_output, labels, current_lr)
                    loss2 = criterion(aux1, labels, current_lr)
                    loss3 = criterion(aux2, labels, current_lr)
                    loss = loss1 + 0.3 * loss2 + 0.3 * loss3
                    outputs_for_acc = main_output
                elif isinstance(outputs, torch.Tensor):
                    loss = criterion(outputs, labels, current_lr)
                    outputs_for_acc = outputs
                else:
                    progress_bar.close()
                    raise TypeError(f"Model output type not recognized: {type(outputs)}")
            else:
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

        val_loss, val_acc, class_performance_stats_val, val_inference_results, val_metrics = evaluate_model(
            model, val_loader, criterion, genotype_map, log_file, loss_type, current_lr
        )

        if class_performance_stats_val:
            print_and_log("\nClass-wise Validation Accuracy:", log_file)
            for class_name in sorted_class_names_from_map:
                stats = class_performance_stats_val.get(class_name, {})
                print_and_log(
                    f"  {class_name} (Index {stats.get('idx', 'N/A')}): {stats.get('acc', 0):.2f}% ({stats.get('correct', 0)}/{stats.get('total', 0)})",
                    log_file)

        summary_msg = (
            f"Epoch {epoch + 1}/{num_epochs} Summary - "
            f"Train Loss: {epoch_train_loss:.4f}, Train Acc: {epoch_train_acc:.2f}%, "
            f"Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.2f}% | "
            f"Prec: {val_metrics['precision_macro'] * 100:.2f}%, "
            f"Rec: {val_metrics['recall_macro'] * 100:.2f}%, "
            f"F1: {val_metrics['f1_macro'] * 100:.2f}% "
            f"(LR: {current_lr:.1e})"
        )

        print_and_log(summary_msg, log_file)

        # Snapshot at epoch 5
        if (epoch + 1) == MIN_SAVE_EPOCH:
            snap_path = os.path.join(output_path, f"model_epoch_{MIN_SAVE_EPOCH}.pth")
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'genotype_map': genotype_map,
                'in_channels': 6
            }, snap_path)
            print_and_log(f"\nSnapshot saved at epoch {epoch + 1}: {snap_path}", log_file)

        # Save on validation improvement (primary: higher acc; tie-breaker: lower loss), after epoch 5
        improved = (val_acc > best_val_acc) or (val_acc == best_val_acc and val_loss < best_val_loss)
        if (epoch + 1) >= MIN_SAVE_EPOCH and improved:
            best_val_acc = val_acc
            best_val_loss = val_loss
            best_epoch = epoch + 1

            best_path = os.path.join(output_path, "model_best.pth")
            torch.save({
                'epoch': best_epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'genotype_map': genotype_map,
                'in_channels': 6,
                'best_val_acc': best_val_acc,
                'best_val_loss': best_val_loss,
            }, best_path)
            print_and_log(f"\nNew BEST at epoch {best_epoch}: Val Acc {best_val_acc:.2f}%, Val Loss {best_val_loss:.4f}. Saved: {best_path}", log_file)

            if save_val_results:
                result_path = os.path.join(output_path, "validation_results_best.json")
                try:
                    with open(result_path, 'w') as f:
                        json.dump({
                            'epoch': best_epoch,
                            'val_acc': best_val_acc,
                            'val_loss': best_val_loss,
                            'inference_results': val_inference_results
                        }, f, indent=4)
                    print_and_log(f"Saved best validation results to {result_path}", log_file)
                except Exception as e:
                    print_and_log(f"Error saving best validation results: {e}", log_file)

        scheduler.step()
        print_and_log("-" * 30, log_file)

    print_and_log(
        f"Training complete. Best epoch: {best_epoch} with Val Acc {best_val_acc:.2f}% | Val Loss {best_val_loss:.4f}. "
        f"Best model: {os.path.join(output_path, 'model_best.pth')}",
        log_file
    )



def evaluate_model(model, data_loader, criterion, genotype_map, log_file, loss_type, current_lr):
    model.eval()
    running_loss_eval = 0.0
    correct_eval = 0
    total_eval = 0
    class_correct_counts = defaultdict(int)
    class_total_counts = defaultdict(int)
    batch_count_eval = 0

    # NEW: for precision/recall/F1
    tp = defaultdict(int)
    fp = defaultdict(int)
    fn = defaultdict(int)

    inference_results = defaultdict(list)
    idx_to_class = {v: k for k, v in genotype_map.items()}

    if not data_loader or len(data_loader) == 0:
        # Return zeros + empty metrics
        metrics = {
            'precision_macro': 0.0, 'recall_macro': 0.0, 'f1_macro': 0.0,
            'precision_weighted': 0.0, 'recall_weighted': 0.0, 'f1_weighted': 0.0
        }
        return 0.0, 0.0, {}, {}, metrics

    with torch.no_grad():
        for images, labels, paths in data_loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            if isinstance(outputs, tuple):
                outputs = outputs[0]

            # Compute loss
            if loss_type == "combined":
                loss = criterion(outputs, labels, current_lr)
            else:
                loss = criterion(outputs, labels)

            running_loss_eval += loss.item()
            batch_count_eval += 1

            _, predicted = torch.max(outputs, 1)
            correct_eval += (predicted == labels).sum().item()
            total_eval += labels.size(0)

            # Accumulate per-class stats
            for i, pred_idx_tensor in enumerate(predicted):
                pred_idx = int(pred_idx_tensor.item())
                true_idx = int(labels[i].item())
                path = paths[i]

                # accuracy & class-wise accuracy
                class_total_counts[true_idx] += 1
                if pred_idx == true_idx:
                    class_correct_counts[true_idx] += 1
                    tp[true_idx] += 1
                else:
                    fp[pred_idx] += 1
                    fn[true_idx] += 1

                predicted_class_name = idx_to_class[pred_idx]
                inference_results[predicted_class_name].append(os.path.basename(path))

    avg_loss_eval = (running_loss_eval / batch_count_eval) if batch_count_eval > 0 else 0.0
    overall_accuracy_eval = (correct_eval / total_eval) * 100 if total_eval > 0 else 0.0

    # Build class-wise stats dict (unchanged behavior)
    class_performance_stats = {}
    if genotype_map:
        for class_name, class_idx in genotype_map.items():
            correct_c = class_correct_counts[class_idx]
            total_c = class_total_counts[class_idx]
            acc_c = (correct_c / total_c) * 100 if total_c > 0 else 0.0
            class_performance_stats[class_name] = {
                'acc': acc_c, 'correct': correct_c, 'total': total_c, 'idx': class_idx
            }
    else:
        print_and_log("Warning: genotype_map is missing in evaluate_model.", log_file)

    # ---- NEW: precision/recall/F1 (macro & weighted) ----
    class_indices = list(genotype_map.values()) if genotype_map else list(set(list(tp.keys()) + list(fp.keys()) + list(fn.keys())))
    precisions, recalls, f1s = [], [], []
    supports = []

    for c in class_indices:
        tpc = tp[c]
        fpc = fp[c]
        fnc = fn[c]
        denom_p = tpc + fpc
        denom_r = tpc + fnc

        pc = (tpc / denom_p) if denom_p > 0 else 0.0
        rc = (tpc / denom_r) if denom_r > 0 else 0.0
        fc = (2 * pc * rc / (pc + rc)) if (pc + rc) > 0 else 0.0

        precisions.append(pc)
        recalls.append(rc)
        f1s.append(fc)
        supports.append(tpc + fnc)  # ground-truth count for class

    # macro
    if len(class_indices) > 0:
        precision_macro = sum(precisions) / len(precisions)
        recall_macro = sum(recalls) / len(recalls)
        f1_macro = sum(f1s) / len(f1s)
    else:
        precision_macro = recall_macro = f1_macro = 0.0

    # weighted (by support)
    total_support = sum(supports)
    if total_support > 0:
        precision_weighted = sum(p * s for p, s in zip(precisions, supports)) / total_support
        recall_weighted = sum(r * s for r, s in zip(recalls, supports)) / total_support
        f1_weighted = sum(f * s for f, s in zip(f1s, supports)) / total_support
    else:
        precision_weighted = recall_weighted = f1_weighted = 0.0

    metrics = {
        'precision_macro': precision_macro,
        'recall_macro': recall_macro,
        'f1_macro': f1_macro,
        'precision_weighted': precision_weighted,
        'recall_weighted': recall_weighted,
        'f1_weighted': f1_weighted,
    }

    # Return with metrics as 5th element
    return avg_loss_eval, overall_accuracy_eval, class_performance_stats, inference_results, metrics


# --- MODIFIED: helpers for path resolution ---
def _read_paths_file(file_path):
    paths = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                paths.append(os.path.abspath(os.path.expanduser(s)))
    except Exception:
        pass
    return paths


def _resolve_data_roots(primary_path, extra_paths, paths_file):
    candidates = []
    if primary_path:
        candidates.append(os.path.abspath(os.path.expanduser(primary_path)))
    if extra_paths:
        for p in extra_paths:
            candidates.append(os.path.abspath(os.path.expanduser(p)))
    if paths_file:
        candidates.extend(_read_paths_file(paths_file))
    seen = set()
    deduped = []
    for p in candidates:
        if p not in seen:
            seen.add(p)
            deduped.append(p)
    if len(deduped) == 0:
        return primary_path
    if len(deduped) == 1:
        return deduped[0]
    return deduped


# Helper to read a list of roots from a txt file (one path per line)
def _read_paths_file(file_path):
    paths = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            paths.append(os.path.abspath(os.path.expanduser(s)))
    return paths

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a Classifier on 6-channel custom .npy dataset")

    # Make data_path OPTIONAL; we will enforce XOR with the files mode
    parser.add_argument("data_path", nargs="?", type=str,
                        help="Dataset root containing 'train/' and 'val/' (Mode A).")

    parser.add_argument("-o", "--output_path", default="./saved_models_6channel", type=str, help="Path to save model")
    parser.add_argument("--depths", type=int, nargs='+', default=[3, 3, 27, 3],
                        help="A list of depths for the ConvNeXt stages (e.g., 3 3 27 3)")
    parser.add_argument("--dims", type=int, nargs='+', default=[192, 384, 768, 1536],
                        help="A list of dimensions for the ConvNeXt stages (e.g., 192 384 768 1536)")

    parser.add_argument("--epochs", type=int, default=50, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=0.0001, help="Initial learning rate")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--num_workers", type=int, default=8, help="Number of workers for data loading")
    # parser.add_argument("--milestone", type=int, default=10, help="Save model every N epochs")

    # Optimizer / scheduler (unchanged)
    parser.add_argument("--warmup_epochs", type=int, default=3, help="Number of epochs for linear LR warmup")
    parser.add_argument("--weight_decay", type=float, default=0.01, help="Weight decay for AdamW optimizer")
    parser.add_argument("--save_val_results", action='store_true', help="Save validation results at each milestone.")
    parser.add_argument("--loss_type", type=str, default="weighted_ce", choices=["combined", "weighted_ce"],
                        help="Loss function to use")

    # NEW: Mode B (files) — both must be provided together when using files mode
    parser.add_argument("--train_data_paths_file", type=str, default=None,
                        help="Text file listing TRAIN dataset roots (one per line).")
    parser.add_argument("--val_data_paths_file", type=str, default=None,
                        help="Text file listing VAL dataset roots (one per line).")
    parser.add_argument("--training_data_ratio", type=float, default=1.0,
                        help="Proportion of training data to use (0–1]. Randomly subsamples the training set.")

    # NEW: single positive-class weight knob (default 88)
    parser.add_argument("--pos_weight", type=float, default=88.0,
                        help="Positive class weight applied to class index 1. Default: 88.0")

    args = parser.parse_args()

    # ---- Enforce: exactly one of (data_path) OR (both files) ----
    has_base = args.data_path is not None
    has_both_files = (args.train_data_paths_file is not None) and (args.val_data_paths_file is not None)

    if not (0 < args.training_data_ratio <= 1.0):
        parser.error(f"--training_data_ratio must be in (0,1], got {args.training_data_ratio}")

    if has_base and has_both_files:
        parser.error("Provide either positional data_path (Mode A) OR both --train_data_paths_file and "
                     "--val_data_paths_file (Mode B), not both.")
    if not has_base and not has_both_files:
        parser.error("You must provide exactly one input mode:\n"
                     "  • Mode A: data_path\n"
                     "  • Mode B: --train_data_paths_file and --val_data_paths_file")

    # Build the argument passed into train_model:
    #  - Mode A: a single string root (backward compatible)
    #  - Mode B: a pair (train_roots, val_roots) for the revised get_data_loader
    if has_base:
        data_path_or_pair = os.path.abspath(os.path.expanduser(args.data_path))
    else:
        train_roots = _read_paths_file(args.train_data_paths_file)
        val_roots   = _read_paths_file(args.val_data_paths_file)
        if not train_roots:
            parser.error(f"--train_data_paths_file is empty or unreadable: {args.train_data_paths_file}")
        if not val_roots:
            parser.error(f"--val_data_paths_file is empty or unreadable: {args.val_data_paths_file}")
        # Pair: get_data_loader(dataset_type="train"/"val") will pick the right side and
        # include BOTH 'train' and 'val' subfolders from each root (per your revised dataloader)
        data_path_or_pair = (train_roots, val_roots)

    # Hand off; train_model still discovers loaders via get_data_loader(...)
    train_model(
        data_path=data_path_or_pair, output_path=args.output_path,
        save_val_results=args.save_val_results,
        num_epochs=args.epochs, learning_rate=args.lr,
        batch_size=args.batch_size, num_workers=args.num_workers,
        loss_type=args.loss_type,
        warmup_epochs=args.warmup_epochs,
        weight_decay=args.weight_decay,
        depths=args.depths,
        dims=args.dims,
        training_data_ratio=args.training_data_ratio,
        pos_weight=args.pos_weight,  # <--- pass through
    )
