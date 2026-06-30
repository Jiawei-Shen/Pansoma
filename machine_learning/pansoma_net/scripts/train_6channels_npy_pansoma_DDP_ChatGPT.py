#!/usr/bin/env python3
import argparse
import json
import os
import sys
import gc
import queue
import threading
from collections import defaultdict
import math

import numpy as np
import torch
import torch.distributed as dist
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.nn.parallel import DistributedDataParallel
from torch.optim.lr_scheduler import SequentialLR, LinearLR, CosineAnnealingLR
from tqdm import tqdm

# ---- env + backend knobs (helps speed) ----
os.environ.setdefault("OMP_NUM_THREADS", "1")
torch.backends.cudnn.benchmark = True

# local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from mynet import ConvNeXtCBAMClassifier  # noqa: E402

# Globals updated in __main__ with DDP
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
IS_MAIN_PROCESS = True  # rank-0 logging only

# Normalization tensors (created once after device is set)
NORM_MEAN = None
NORM_STD = None


# ---------------- Losses ----------------
class MultiClassFocalLoss(nn.Module):
    def __init__(self, gamma=2.0, weight=None, reduction='mean'):
        super().__init__()
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

        focal_loss = -1 * (1 - pt) ** self.gamma * log_pt
        if self.reduction == 'mean':
            return focal_loss.mean()
        if self.reduction == 'sum':
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
        return focal_weight * self.focal_loss(logits, targets) + wce_weight * self.wce_loss(logits, targets)


# ---------------- Utils ----------------
def init_weights(m):
    if isinstance(m, (nn.Conv2d, nn.Linear)):
        nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)


def print_and_log(message, log_path):
    if not IS_MAIN_PROCESS:
        return
    print(message, flush=True)
    with open(log_path, 'a', encoding='utf-8') as f:
        f.write(message + '\n')


def _state_dict(m):
    return m.module.state_dict() if hasattr(m, "module") else m.state_dict()


def _load_state_dict(m, state):
    if hasattr(m, "module"):
        m.module.load_state_dict(state)
    else:
        m.load_state_dict(state)


def _unique_path(path: str) -> str:
    if not os.path.exists(path):
        return path
    base, ext = os.path.splitext(path)
    n = 2
    cand = f"{base}_v{n}{ext}"
    while os.path.exists(cand):
        n += 1
        cand = f"{base}_v{n}{ext}"
    return cand


# ---------------- Shard discovery ----------------
def _discover_shards_for_root(root):
    """
    root: .../ALL_chr_merged_REAL_sharded_npy
    Expects:
      root/train/shard_XXXXX_data.npy
      root/train/shard_XXXXX_labels.npy
      root/val/shard_XXXXX_data.npy
      root/val/shard_XXXXX_labels.npy
    Returns:
      train_shards, val_shards (lists of dicts)
    """
    train_dir = os.path.join(root, "train")
    val_dir = os.path.join(root, "val")

    def _discover_in_dir(d):
        shards = []
        if not os.path.isdir(d):
            return shards
        data_paths = sorted(
            p for p in (os.path.join(d, f) for f in os.listdir(d))
            if p.endswith("_data.npy")
        )
        for dp in data_paths:
            lp = dp.replace("_data.npy", "_labels.npy")
            if not os.path.exists(lp):
                continue
            y_tmp = np.load(lp, mmap_mode="r")
            n = int(y_tmp.shape[0])
            del y_tmp
            shards.append({"x_path": dp, "y_path": lp, "num_samples": n})
        return shards

    train_shards = _discover_in_dir(train_dir)
    val_shards = _discover_in_dir(val_dir)
    return train_shards, val_shards


def _discover_all_shards(roots):
    """
    roots: list of sharded_npy roots
    Returns: all_train_shards, all_val_shards, genotype_map
    """
    all_train = []
    all_val = []
    for r in roots:
        tr, va = _discover_shards_for_root(r)
        all_train.extend(tr)
        all_val.extend(va)

    genotype_map = {"false": 0, "true": 1}
    return all_train, all_val, genotype_map


# ---------------- Shard prefetcher & iterators ----------------
class ShardPrefetcher:
    """
    Background prefetcher:
      - Loads each shard into memory (mmap) in a background thread.
      - Main thread consumes (local_idx, shard_dict, x_np, y_np).
    """

    def __init__(self, shard_specs, max_prefetch: int = 2):
        self.shard_specs = list(shard_specs)
        self.max_prefetch = max_prefetch
        self._queue = queue.Queue(maxsize=max_prefetch)
        self._stop = False
        self._thread = threading.Thread(target=self._worker, daemon=True)
        self._thread.start()

    def _worker(self):
        try:
            for local_order_idx, shard in enumerate(self.shard_specs):
                if self._stop:
                    break
                x_path = shard["x_path"]
                y_path = shard["y_path"]
                # IMPORTANT: use mmap to avoid loading full shards into RAM
                x_np = np.load(x_path, mmap_mode="r")
                y_np = np.load(y_path, mmap_mode="r")
                self._queue.put((local_order_idx, shard, x_np, y_np))
            self._queue.put(None)
        except Exception as e:
            self._queue.put(e)

    def __iter__(self):
        while True:
            item = self._queue.get()
            if item is None:
                return
            if isinstance(item, Exception):
                raise item
            yield item

    def close(self):
        self._stop = True
        try:
            self._queue.put_nowait(None)
        except Exception:
            pass


def _maybe_normalize(batch_x, do_norm: bool):
    global NORM_MEAN, NORM_STD
    if not do_norm:
        return batch_x
    # NORM_MEAN / NORM_STD are created once after device is set
    return (batch_x - NORM_MEAN) / NORM_STD


def _iter_sharded_batches(
    shards,
    batch_size: int,
    epoch: int,
    shuffle_shards: bool = True,
    shuffle_within_shard: bool = True,
    drop_last: bool = False,
    training_data_ratio: float = 1.0,
    max_steps: int = None,         # NEW
    do_normalize: bool = True,            # NEW
):
    """
    Non-DDP iterator over shards for this process.
    Uses ShardPrefetcher to overlap 'np.load' of next shard with training on current shard.

    NEW:
      - max_steps: stop yielding after max_steps batches (so generator closes cleanly).
      - do_normalize: apply on-GPU normalization using cached tensors.
    """
    num_shards = len(shards)
    if num_shards == 0:
        return

    shard_indices = np.arange(num_shards)
    if shuffle_shards:
        rng = np.random.default_rng(1234 + epoch)
        rng.shuffle(shard_indices)

    if training_data_ratio < 1.0:
        num_keep = max(1, int(num_shards * training_data_ratio))
        shard_indices = shard_indices[:num_keep]

    ordered_shards = [shards[i] for i in shard_indices]
    prefetcher = ShardPrefetcher(ordered_shards, max_prefetch=1)  # reduced prefetch

    yielded = 0
    try:
        for local_order_idx, shard, x_np, y_np in prefetcher:
            n = int(shard["num_samples"])
            if n <= 0:
                del x_np, y_np
                continue

            idx = np.arange(n)
            if shuffle_within_shard:
                rng = np.random.default_rng(epoch * 1337 + local_order_idx)
                rng.shuffle(idx)

            start = 0
            while start < n:
                if max_steps is not None and yielded >= max_steps:
                    return  # IMPORTANT: ends generator -> finally runs -> prefetch stops

                end = start + batch_size
                if end > n:
                    if drop_last:
                        break
                    end = n
                batch_idx = idx[start:end]

                batch_x_np = x_np[batch_idx]
                batch_y_np = y_np[batch_idx].astype(np.int64, copy=False)

                batch_x = torch.from_numpy(batch_x_np).float().pin_memory()
                batch_y = torch.from_numpy(batch_y_np).pin_memory()

                batch_x = batch_x.to(device, non_blocking=True)
                batch_y = batch_y.to(device, non_blocking=True)

                batch_x = _maybe_normalize(batch_x, do_normalize)

                yield batch_x, batch_y
                yielded += 1

                start = end

            del x_np, y_np
            gc.collect()
    finally:
        prefetcher.close()


def _iter_sharded_batches_ddp(
    shards,
    batch_size: int,
    epoch: int,
    world_size: int,
    rank: int,
    shuffle_shards: bool = True,
    shuffle_within_shard: bool = True,
    drop_last: bool = False,
    training_data_ratio: float = 1.0,
    max_steps: int = None,         # NEW
    do_normalize: bool = True,            # NEW
):
    """
    DDP-aware shard iterator.

    NEW:
      - max_steps: stop yielding after max_steps batches (so generator closes cleanly).
      - do_normalize: apply on-GPU normalization using cached tensors.
    """
    num_shards = len(shards)
    if num_shards == 0:
        return

    shard_indices = np.arange(num_shards)
    if shuffle_shards:
        rng = np.random.default_rng(4321 + epoch)
        rng.shuffle(shard_indices)

    if training_data_ratio < 1.0:
        num_keep = max(world_size, int(round(num_shards * training_data_ratio)))
    else:
        num_keep = num_shards

    num_keep = (num_keep // world_size) * world_size
    if num_keep == 0:
        return

    shard_indices = shard_indices[:num_keep]
    local_indices = shard_indices[rank::world_size]
    if len(local_indices) == 0:
        return

    local_shards = [shards[i] for i in local_indices]
    prefetcher = ShardPrefetcher(local_shards, max_prefetch=1)  # reduced prefetch

    yielded = 0
    try:
        for local_order_idx, shard, x_np, y_np in prefetcher:
            n = int(shard["num_samples"])
            if n <= 0:
                del x_np, y_np
                continue

            idx = np.arange(n)
            if shuffle_within_shard:
                rng = np.random.default_rng(epoch * 7331 + rank * 17 + local_order_idx)
                rng.shuffle(idx)

            start = 0
            while start < n:
                if max_steps is not None and yielded >= max_steps:
                    return  # IMPORTANT: ends generator -> finally runs -> prefetch stops

                end = start + batch_size
                if end > n:
                    if drop_last:
                        break
                    end = n
                batch_idx = idx[start:end]

                batch_x_np = x_np[batch_idx]
                batch_y_np = y_np[batch_idx].astype(np.int64, copy=False)

                batch_x = torch.from_numpy(batch_x_np).float().pin_memory()
                batch_y = torch.from_numpy(batch_y_np).pin_memory()

                batch_x = batch_x.to(device, non_blocking=True)
                batch_y = batch_y.to(device, non_blocking=True)

                batch_x = _maybe_normalize(batch_x, do_normalize)

                yield batch_x, batch_y
                yielded += 1

                start = end

            del x_np, y_np
            gc.collect()
    finally:
        prefetcher.close()


# ---------------- Eval (works with any batch iterator) ----------------
def evaluate_model(model,
                   data_iter,
                   criterion,
                   genotype_map,
                   log_file,
                   loss_type,
                   current_lr,
                   ddp=False,
                   world_size=1,
                   collect_infer=False):
    model.eval()
    num_classes = len(genotype_map) if genotype_map else 0

    correct_eval = torch.zeros(1, device=device, dtype=torch.long)
    total_eval = torch.zeros(1, device=device, dtype=torch.long)
    loss_sum = torch.zeros(1, device=device, dtype=torch.float)

    tp = torch.zeros(num_classes, device=device, dtype=torch.long)
    fp = torch.zeros(num_classes, device=device, dtype=torch.long)
    fn = torch.zeros(num_classes, device=device, dtype=torch.long)
    class_correct_counts = torch.zeros(num_classes, device=device, dtype=torch.long)
    class_total_counts = torch.zeros(num_classes, device=device, dtype=torch.long)

    batch_count_eval = 0
    inference_results = defaultdict(list) if collect_infer else {}
    idx_to_class = {v: k for k, v in genotype_map.items()} if genotype_map else {}

    with torch.no_grad():
        for batch in data_iter:
            images, labels = batch
            images = images.to(device, non_blocking=True)
            labels = labels.to(device, non_blocking=True)

            outputs = model(images)
            if isinstance(outputs, tuple):
                outputs = outputs[0]

            if loss_type == "combined":
                loss = criterion(outputs, labels, current_lr)
            else:
                loss = criterion(outputs, labels)

            loss_sum += loss.detach()
            batch_count_eval += 1

            _, predicted = torch.max(outputs, 1)
            correct_eval += (predicted == labels).sum()
            total_eval += labels.size(0)

            for i in range(labels.size(0)):
                pred_idx = int(predicted[i])
                true_idx = int(labels[i])
                class_total_counts[true_idx] += 1
                if pred_idx == true_idx:
                    class_correct_counts[true_idx] += 1
                    tp[true_idx] += 1
                else:
                    if pred_idx < num_classes:
                        fp[pred_idx] += 1
                    fn[true_idx] += 1

                if collect_infer and genotype_map:
                    # Minimal structure placeholder; keep behavior disabled unless requested
                    # inference_results[idx_to_class[true_idx]].append(...)
                    pass

    if ddp and world_size > 1 and dist.is_initialized():
        dist.all_reduce(correct_eval, op=dist.ReduceOp.SUM)
        dist.all_reduce(total_eval, op=dist.ReduceOp.SUM)
        dist.all_reduce(loss_sum, op=dist.ReduceOp.SUM)
        if num_classes > 0:
            dist.all_reduce(tp, op=dist.ReduceOp.SUM)
            dist.all_reduce(fp, op=dist.ReduceOp.SUM)
            dist.all_reduce(fn, op=dist.ReduceOp.SUM)
            dist.all_reduce(class_correct_counts, op=dist.ReduceOp.SUM)
            dist.all_reduce(class_total_counts, op=dist.ReduceOp.SUM)

    denom_batches = max(1, batch_count_eval * (world_size if (ddp and world_size > 1) else 1))
    avg_loss_eval = loss_sum.item() / denom_batches
    overall_accuracy_eval = (correct_eval.item() / max(1, total_eval.item())) * 100.0

    class_performance_stats = {}
    if genotype_map:
        for class_name, class_idx in genotype_map.items():
            correct_c = int(class_correct_counts[class_idx].item())
            total_c = int(class_total_counts[class_idx].item())
            acc_c = (correct_c / total_c * 100.0) if total_c > 0 else 0.0
            class_performance_stats[class_name] = {
                'acc': acc_c, 'correct': correct_c, 'total': total_c, 'idx': class_idx
            }

    pos_idx = None
    if genotype_map:
        for name, idx in genotype_map.items():
            if str(name).lower() == "true":
                pos_idx = idx
                break
    if pos_idx is None:
        pos_idx = 1 if num_classes > 1 else 0

    tpc = float(tp[pos_idx].item() if pos_idx < num_classes else 0.0)
    fpc = float(fp[pos_idx].item() if pos_idx < num_classes else 0.0)
    fnc = float(fn[pos_idx].item() if pos_idx < num_classes else 0.0)

    precision_true = (tpc / (tpc + fpc)) if (tpc + fpc) > 0 else 0.0
    recall_true = (tpc / (tpc + fnc)) if (tpc + fnc) > 0 else 0.0
    f1_true = (2 * precision_true * recall_true / (precision_true + recall_true)) if (precision_true + recall_true) > 0 else 0.0

    metrics = {
        'precision_true': precision_true,
        'recall_true': recall_true,
        'f1_true': f1_true,
        'pos_class_idx': pos_idx,
    }
    return avg_loss_eval, overall_accuracy_eval, class_performance_stats, inference_results, metrics


# ---------------- Train (sharded NPY only) ----------------
def train_model_sharded(train_shards,
                        val_shards,
                        output_path,
                        save_val_results=False,
                        num_epochs=100,
                        learning_rate=1e-4,
                        batch_size=32,
                        loss_type='weighted_ce',
                        warmup_epochs=10,
                        weight_decay=0.05,
                        depths=None,
                        dims=None,
                        training_data_ratio=1.0,
                        ddp=False,
                        world_size=1,
                        rank=0,
                        resume=None,
                        pos_weight=88.0,
                        genotype_map=None):
    os.makedirs(output_path, exist_ok=True)
    log_file = os.path.join(output_path, "training_log_6ch_sharded.txt")
    if os.path.exists(log_file) and IS_MAIN_PROCESS:
        os.remove(log_file)

    MIN_SAVE_EPOCH = 5

    if not (0 < training_data_ratio <= 1.0):
        raise ValueError(f"--training_data_ratio must be in (0,1], got {training_data_ratio}")

    if genotype_map is None:
        genotype_map = {"false": 0, "true": 1}
    num_classes = len(genotype_map)

    print_and_log(f"Using device: {device}", log_file)
    print_and_log(f"Initial Learning Rate: {learning_rate:.1e}", log_file)
    print_and_log(f"[Sharded NPY mode] #train_shards={len(train_shards)}, #val_shards={len(val_shards)}", log_file)

    approx_used_global_batches = None

    if train_shards:
        approx_train_samples = sum(int(s["num_samples"]) for s in train_shards)
        approx_batches = approx_train_samples // (batch_size * max(1, world_size))
        print_and_log(
            f"  Approx total train samples: {approx_train_samples:,} "
            f"({approx_batches:,} global batches @ batch_size={batch_size})",
            log_file
        )

        approx_used_global_batches = max(1, approx_batches)

        if training_data_ratio < 1.0:
            approx_used_samples = int(approx_train_samples * training_data_ratio)
            approx_used_global_batches = max(
                1, approx_used_samples // (batch_size * max(1, world_size))
            )
            print_and_log(
                f"  training_data_ratio={training_data_ratio:.3f} -> "
                f"~{approx_used_samples:,} samples, ~{approx_used_global_batches:,} global batches per epoch "
                f"(assuming shards are similar size)",
                log_file
            )
    else:
        approx_used_global_batches = None

    print_and_log(f"Number of classes: {num_classes}", log_file)
    print_and_log(f"depths={depths}, dims={dims}. \n", log_file)

    # ---- Model / parallelism ----
    model = ConvNeXtCBAMClassifier(in_channels=6, class_num=num_classes,
                                   depths=depths, dims=dims).to(device)

    if ddp:
        print_and_log(f"Wrapping model in DistributedDataParallel on {device}.", log_file)
        model = DistributedDataParallel(
            model,
            device_ids=[rank] if device.type == "cuda" else None,
            output_device=rank if device.type == "cuda" else None,
            gradient_as_bucket_view=True,
            broadcast_buffers=False,
        )

    model.apply(init_weights)

    # ---- Loss / Optim / Sched ----
    pos_weight_value = float(pos_weight)
    class_weights = torch.ones(num_classes, device=device, dtype=torch.float32)
    if num_classes >= 2:
        class_weights[1] = pos_weight_value

    print_and_log(f"Class weights: {class_weights.tolist()} (from --pos_weight={pos_weight_value})", log_file)

    if loss_type == "combined":
        criterion = CombinedFocalWeightedCELoss(initial_lr=learning_rate, pos_weight=class_weights)
        print_and_log("Using Combined(Focal + Weighted CE) Loss.", log_file)
    elif loss_type == "weighted_ce":
        criterion = nn.CrossEntropyLoss(weight=class_weights)
        print_and_log("Using Weighted CE Loss.", log_file)
    else:
        raise ValueError(f"Unsupported loss_type: {loss_type}")

    optimizer = optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    warmup_scheduler = LinearLR(optimizer, start_factor=0.01, total_iters=warmup_epochs)
    main_scheduler = CosineAnnealingLR(optimizer, T_max=max(1, num_epochs - warmup_epochs), eta_min=0)
    scheduler = SequentialLR(optimizer, schedulers=[warmup_scheduler, main_scheduler],
                             milestones=[warmup_epochs])

    # ---- Resume ----
    start_epoch = 0
    best_epoch = 0
    best_f1_true = float("-inf")
    best_val_acc = float("-inf")
    best_val_loss = float("inf")
    best_rec_true = 0.0
    last_best_ckpt_path = None

    if resume is not None and os.path.isfile(resume):
        try:
            checkpoint = torch.load(resume, map_location=device)
            _load_state_dict(model, checkpoint['model_state_dict'])
            if 'optimizer_state_dict' in checkpoint:
                optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            if 'scheduler_state_dict' in checkpoint:
                scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
            start_epoch = int(checkpoint.get('epoch', 0))
            best_f1_true = float(checkpoint.get('best_f1_true', best_f1_true))
            best_rec_true = float(checkpoint.get('best_rec_true', best_rec_true))
            best_val_acc = float(checkpoint.get('best_val_acc', best_val_acc))
            best_val_loss = float(checkpoint.get('best_val_loss', best_val_loss))
            best_epoch = int(checkpoint.get('epoch', best_epoch))
            print_and_log(f"Resumed from '{resume}' at epoch {start_epoch}.", log_file)
        except Exception as e:
            print_and_log(f"WARNING: Failed to load checkpoint '{resume}': {e}", log_file)

    sorted_class_names_from_map = sorted(genotype_map.keys(), key=lambda k: genotype_map[k])

    # ---- Train loop ----
    for epoch in range(start_epoch, num_epochs):
        model.train()
        running_loss = 0.0
        correct_train = 0
        total_train = 0
        batch_count = 0

        current_lr = optimizer.param_groups[0]['lr']

        if ddp and world_size > 1:
            desc = f"[Rank {rank}] Epoch {epoch + 1}/{num_epochs} LR: {current_lr:.1e}"
        else:
            desc = f"Epoch {epoch + 1}/{num_epochs} LR: {current_lr:.1e}"

        # NEW: per-rank total for tqdm and step-cap
        if approx_used_global_batches is not None:
            steps_per_rank = math.ceil(approx_used_global_batches / max(1, world_size))
            total_local = steps_per_rank
        else:
            steps_per_rank = None
            total_local = None

        # Build iterators (NEW: pass max_steps so iterator ends cleanly)
        if ddp and world_size > 1:
            train_iter = _iter_sharded_batches_ddp(
                shards=train_shards,
                batch_size=batch_size,
                epoch=epoch,
                world_size=world_size,
                rank=rank,
                shuffle_shards=True,
                shuffle_within_shard=True,
                drop_last=False,
                training_data_ratio=training_data_ratio,
                max_steps=steps_per_rank,     # NEW
                do_normalize=True,            # keep behavior
            )
        else:
            train_iter = _iter_sharded_batches(
                shards=train_shards,
                batch_size=batch_size,
                epoch=epoch,
                shuffle_shards=True,
                shuffle_within_shard=True,
                drop_last=False,
                training_data_ratio=training_data_ratio,
                max_steps=steps_per_rank,     # NEW
                do_normalize=True,            # keep behavior
            )

        # Ensure tqdm closes promptly
        with tqdm(
            train_iter,
            desc=desc,
            total=total_local,
            leave=True,
            disable=not IS_MAIN_PROCESS
        ) as progress_bar:
            for images, labels in progress_bar:
                optimizer.zero_grad(set_to_none=True)
                outputs = model(images)

                if loss_type == "combined":
                    if isinstance(outputs, tuple) and len(outputs) == 3:
                        main_output, aux1, aux2 = outputs
                        loss = (criterion(main_output, labels, current_lr)
                                + 0.3 * criterion(aux1, labels, current_lr)
                                + 0.3 * criterion(aux2, labels, current_lr))
                        outputs_for_acc = main_output
                    else:
                        loss = criterion(outputs, labels, current_lr)
                        outputs_for_acc = outputs
                else:
                    if isinstance(outputs, tuple) and len(outputs) == 3:
                        main_output, aux1, aux2 = outputs
                        loss = (criterion(main_output, labels)
                                + 0.3 * criterion(aux1, labels)
                                + 0.3 * criterion(aux2, labels))
                        outputs_for_acc = main_output
                    else:
                        loss = criterion(outputs, labels)
                        outputs_for_acc = outputs

                loss.backward()
                optimizer.step()

                running_loss += loss.item()
                batch_count += 1
                _, predicted = torch.max(outputs_for_acc, 1)
                correct_train += (predicted == labels).sum().item()
                total_train += labels.size(0)

                if IS_MAIN_PROCESS and total_train > 0 and batch_count > 0:
                    avg_loss_train = running_loss / batch_count
                    avg_acc_train = (correct_train / total_train) * 100.0
                    progress_bar.set_postfix(loss=f"{avg_loss_train:.4f}",
                                             acc=f"{avg_acc_train:.2f}%")

        epoch_train_loss = (running_loss / batch_count) if batch_count > 0 else 0.0
        epoch_train_acc = (correct_train / max(1, total_train)) * 100.0

        # ---- Validation on shards ----
        if val_shards:
            if ddp and world_size > 1:
                val_iter = _iter_sharded_batches_ddp(
                    shards=val_shards,
                    batch_size=batch_size,
                    epoch=epoch,
                    world_size=world_size,
                    rank=rank,
                    shuffle_shards=False,
                    shuffle_within_shard=False,
                    drop_last=False,
                    training_data_ratio=1.0,
                    max_steps=None,        # use all val
                    do_normalize=True,
                )
            else:
                val_iter = _iter_sharded_batches(
                    shards=val_shards,
                    batch_size=batch_size,
                    epoch=epoch,
                    shuffle_shards=False,
                    shuffle_within_shard=False,
                    drop_last=False,
                    training_data_ratio=1.0,
                    max_steps=None,        # use all val
                    do_normalize=True,
                )

            val_loss, val_acc, class_stats_val, val_infer_lists, val_metrics = evaluate_model(
                model, val_iter, criterion, genotype_map, log_file, loss_type, current_lr,
                ddp=ddp, world_size=world_size,
                collect_infer=save_val_results
            )
        else:
            val_loss, val_acc, class_stats_val, val_infer_lists, val_metrics = (
                0.0, 0.0, {}, {}, {
                    'precision_true': 0.0,
                    'recall_true': 0.0,
                    'f1_true': 0.0,
                    'pos_class_idx': None
                }
            )

        if IS_MAIN_PROCESS and class_stats_val:
            print_and_log("\nClass-wise Validation Accuracy:", log_file)
            for class_name in sorted_class_names_from_map:
                s = class_stats_val.get(class_name, {})
                print_and_log(
                    f"  {class_name} (idx {s.get('idx','N/A')}): "
                    f"{s.get('acc',0):.2f}% ({s.get('correct',0)}/{s.get('total',0)})",
                    log_file
                )

        val_prec_true = val_metrics.get('precision_true', 0.0)
        val_rec_true = val_metrics.get('recall_true', 0.0)
        val_f1_true = val_metrics.get('f1_true', 0.0)
        pos_idx = val_metrics.get('pos_class_idx', None)

        print_and_log(
            f"Epoch {epoch + 1}/{num_epochs} "
            f"| Train Loss {epoch_train_loss:.4f} Acc {epoch_train_acc:.2f}% "
            f"| Val Loss {val_loss:.4f} Acc {val_acc:.2f}% "
            f"| Prec(true) {val_prec_true*100:.2f}% Rec(true) {val_rec_true*100:.2f}% "
            f"F1(true) {val_f1_true:.4f} "
            f"(LR {current_lr:.1e}{', pos_idx='+str(pos_idx) if pos_idx is not None else ''})",
            log_file
        )

        improved = (
            (val_f1_true > best_f1_true) or
            (val_f1_true == best_f1_true and val_rec_true > best_rec_true) or
            (val_f1_true == best_f1_true and val_rec_true == best_rec_true and val_loss < best_val_loss)
        )

        if (epoch + 1) == MIN_SAVE_EPOCH and IS_MAIN_PROCESS:
            snap_path = _unique_path(os.path.join(output_path, f"model_epoch_{MIN_SAVE_EPOCH}.pth"))
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': _state_dict(model),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'genotype_map': genotype_map,
                'in_channels': 6
            }, snap_path)
            print_and_log(f"Snapshot saved at epoch {epoch + 1}: {snap_path}", log_file)

        if (epoch + 1) >= MIN_SAVE_EPOCH and improved and IS_MAIN_PROCESS:
            best_f1_true = val_f1_true
            best_rec_true = val_rec_true
            best_val_acc = val_acc
            best_val_loss = val_loss
            best_epoch = epoch + 1

            best_path = _unique_path(
                os.path.join(output_path, f"model_e{best_epoch:03d}_f1_{best_f1_true:.4f}.pth")
            )
            payload = {
                'epoch': best_epoch,
                'model_state_dict': _state_dict(model),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'genotype_map': genotype_map,
                'in_channels': 6,
                'best_f1_true': best_f1_true,
                'best_rec_true': best_rec_true,
                'best_val_acc': best_val_acc,
                'best_val_loss': best_val_loss,
            }
            torch.save(payload, best_path)
            last_best_ckpt_path = best_path
            print_and_log(
                f"New BEST @epoch {best_epoch}: F1(true) {best_f1_true:.4f} | "
                f"Rec(true) {best_rec_true*100:.2f}% | Val Acc {best_val_acc:.2f}% | "
                f"Val Loss {best_val_loss:.4f}. Saved: {best_path}",
                log_file
            )

            if save_val_results:
                result_path = _unique_path(
                    os.path.join(output_path, f"validation_results_e{best_epoch:03d}_f1_{best_f1_true:.4f}.json")
                )
                try:
                    with open(result_path, 'w') as f:
                        json.dump({
                            'epoch': best_epoch,
                            'f1_true': best_f1_true,
                            'recall_true': best_rec_true,
                            'val_acc': best_val_acc,
                            'val_loss': best_val_loss,
                            'inference_results': val_infer_lists
                        }, f, indent=4)
                    print_and_log(f"Saved best validation results to {result_path}", log_file)
                except Exception as e:
                    print_and_log(f"Error saving best validation results: {e}", log_file)

        scheduler.step()
        print_and_log("-" * 30, log_file)

        # Explicitly drop iterators (helps release refs promptly)
        del train_iter
        if val_shards:
            del val_iter
        gc.collect()

    final_msg = (
        f"Training complete. Best epoch: {best_epoch} "
        f"| F1(true) {best_f1_true:.4f} | Rec(true) {best_rec_true*100:.2f}% "
        f"| Val Acc {best_val_acc:.2f}% | Val Loss {best_val_loss:.4f}. "
        f"{'Best model: '+last_best_ckpt_path if last_best_ckpt_path else 'No best checkpoint saved.'}"
    )
    print_and_log(final_msg, log_file)


# ---------------- Main ----------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train ConvNeXtCBAM on 6-channel sharded .npy dataset (DDP-friendly, shard prefetching)."
    )

    parser.add_argument("--data_paths", type=str, nargs='+', default=None,
                        help="Sharded NPY roots. Each should contain 'train/' and 'val/' with shard_*_data.npy/labels.npy.")

    parser.add_argument("-o", "--output_path", default="./saved_models_6channel_sharded", type=str,
                        help="Path to save model")
    parser.add_argument("--depths", type=int, nargs='+', default=[3, 3, 27, 3])
    parser.add_argument("--dims", type=int, nargs='+', default=[192, 384, 768, 1536])
    parser.add_argument("--epochs", type=int, default=70)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--batch_size", type=int, default=32)

    parser.add_argument("--num_workers", type=int, default=0,
                        help="Unused in sharded mode; present for compatibility.")
    parser.add_argument("--prefetch_factor", type=int, default=4,
                        help="Unused in sharded mode; present for compatibility.")
    parser.add_argument("--mp_context", type=str, default=None,
                        choices=[None, "fork", "forkserver", "spawn"],
                        help="Unused in sharded mode; present for compatibility.")

    parser.add_argument("--warmup_epochs", type=int, default=3)
    parser.add_argument("--weight_decay", type=float, default=0.01)
    parser.add_argument("--save_val_results", action='store_true')
    parser.add_argument("--loss_type", type=str, default="weighted_ce",
                        choices=["combined", "weighted_ce"])

    parser.add_argument("--training_data_ratio", type=float, default=1.0)

    parser.add_argument("--ddp", action="store_true")
    parser.add_argument("--data_parallel", action="store_true",
                        help="Ignored in shard mode; use --ddp instead.")
    parser.add_argument("--local_rank", type=int, default=None,
                        help="Torch launcher may pass this, but we mostly use LOCAL_RANK env.")

    parser.add_argument("--resume", type=str, default=None)
    parser.add_argument("--pos_weight", type=float, default=88.0)

    args, _unknown = parser.parse_known_args()

    if not args.data_paths:
        parser.error("This script is shard-only. Please provide --data_paths root1 root2 ...")

    roots = [os.path.abspath(os.path.expanduser(p)) for p in args.data_paths]

    # DDP init
    if args.ddp:
        if not torch.cuda.is_available():
            raise RuntimeError("DDP requires CUDA available.")
        local_rank_env = os.environ.get("LOCAL_RANK", None)
        if local_rank_env is None:
            local_rank = 0 if args.local_rank is None else int(args.local_rank)
        else:
            local_rank = int(local_rank_env)

        torch.cuda.set_device(local_rank)
        device = torch.device(f"cuda:{local_rank}")
        dist.init_process_group(backend="nccl", init_method="env://")
        world_size = dist.get_world_size()
        rank = dist.get_rank()
        IS_MAIN_PROCESS = (rank == 0)
        if IS_MAIN_PROCESS:
            print(f"[DDP] World size={world_size} | Local rank={local_rank} | Global rank={rank}")
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        world_size = 1
        rank = 0
        IS_MAIN_PROCESS = True

    # Create normalization tensors ONCE (per process, after device is set)
    NORM_MEAN = torch.tensor(
        [18.41781616, 12.64912987, -0.54525274, 24.72385406, 4.69061136, 0.28135515],
        device=device, dtype=torch.float32
    ).view(1, 6, 1, 1)
    NORM_STD = torch.tensor(
        [25.02832222, 14.80963230, 0.61813378, 29.97283554, 7.92317915, 0.76590837],
        device=device, dtype=torch.float32
    ).view(1, 6, 1, 1)

    # Discover shards
    train_shards, val_shards, genotype_map = _discover_all_shards(roots)
    if IS_MAIN_PROCESS:
        print(f"[Sharded NPY mode] roots={roots}")
        print(f"  Train shards: {len(train_shards)} | Val shards: {len(val_shards)}")
        approx_train_samples = sum(int(s["num_samples"]) for s in train_shards) if train_shards else 0
        approx_batches = (approx_train_samples // (args.batch_size * max(1, world_size))) if approx_train_samples > 0 else 0
        print(
            f"  Approx total train samples: {approx_train_samples:,} "
            f"({approx_batches:,} global batches @ batch_size={args.batch_size})"
        )

    train_model_sharded(
        train_shards=train_shards,
        val_shards=val_shards,
        output_path=args.output_path,
        save_val_results=args.save_val_results,
        num_epochs=args.epochs,
        learning_rate=args.lr,
        batch_size=args.batch_size,
        loss_type=args.loss_type,
        warmup_epochs=args.warmup_epochs,
        weight_decay=args.weight_decay,
        depths=args.depths,
        dims=args.dims,
        training_data_ratio=args.training_data_ratio,
        ddp=args.ddp,
        world_size=world_size,
        rank=rank,
        resume=args.resume,
        pos_weight=args.pos_weight,
        genotype_map=genotype_map,
    )

    if args.ddp and dist.is_initialized():
        dist.destroy_process_group()