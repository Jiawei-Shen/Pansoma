#!/usr/bin/env python3
import argparse
import os
import sys
import glob
import random
import torch
import torch.distributed as dist
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.nn.parallel import DistributedDataParallel
from torch.optim.lr_scheduler import SequentialLR, LinearLR, CosineAnnealingLR
from torch.utils.data import DataLoader
from torch.utils.data.distributed import DistributedSampler
from torchvision import transforms
from tqdm import tqdm

os.environ.setdefault("OMP_NUM_THREADS", "1")
torch.backends.cudnn.benchmark = True

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
try:
    from mynet import ConvNeXtCBAMClassifier
    from dataset_pansoma_npy_sharded_6ch_DDP_Gemini import get_data_loader, NpyDataset
except ImportError:
    pass

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
IS_MAIN_PROCESS = True


# ---------------- Losses & Utils ----------------
class CombinedFocalWeightedCELoss(nn.Module):
    def __init__(self, initial_lr, pos_weight=None, gamma=2.0):
        super().__init__()
        self.ce = nn.CrossEntropyLoss(weight=pos_weight)
        self.gamma = gamma
        self.pos_weight = pos_weight
        self.initial_lr = initial_lr

    def forward(self, logits, targets, current_lr):
        log_pt = F.log_softmax(logits, dim=1).gather(1, targets.view(-1, 1)).squeeze(1)
        pt = torch.exp(log_pt)
        if self.pos_weight is not None: log_pt = log_pt * self.pos_weight.gather(0, targets)
        focal = -1 * (1 - pt) ** self.gamma * log_pt
        fw = 1.0 - (current_lr / self.initial_lr)
        return fw * focal.mean() + (1.0 - fw) * self.ce(logits, targets)


def init_weights(m):
    if isinstance(m, (nn.Conv2d, nn.Linear)):
        nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None: nn.init.constant_(m.bias, 0)


def print_and_log(msg, path):
    if IS_MAIN_PROCESS:
        print(msg, flush=True)
        with open(path, 'a') as f: f.write(msg + '\n')


def _unique_path(path):
    if not os.path.exists(path): return path
    base, ext = os.path.splitext(path)
    n = 2
    while os.path.exists(f"{base}_v{n}{ext}"): n += 1
    return f"{base}_v{n}{ext}"


def _read_paths_file(file_path):
    paths = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                paths.append(os.path.abspath(os.path.expanduser(line.strip())))
    return paths


# ---------------- Core Helper ----------------
def _find_all_shards(roots, subfolder="train"):
    files = []
    for r in roots:
        r = os.path.abspath(os.path.expanduser(r))
        files.extend(glob.glob(os.path.join(r, "*_data.npy")))
        files.extend(glob.glob(os.path.join(r, subfolder, "*_data.npy")))
    return sorted(list(set(files)))


def _get_normalization_transform():
    return transforms.Compose([
        transforms.Normalize(mean=[18.417, 12.649, -0.545, 24.723, 4.690, 0.281],
                             std=[25.028, 14.809, 0.618, 29.972, 7.923, 0.765])
    ])


# ---------------- Train ----------------
def train_model(data_path, output_path, num_epochs=70, learning_rate=1e-4, batch_size=32, num_workers=4,
                ddp=False, local_rank=0, resume=None, pos_weight=88.0, loss_type='weighted_ce',
                depths=None, dims=None, training_data_ratio=1.0, warmup_epochs=3, weight_decay=0.01,
                prefetch_factor=4, mp_context=None, **kwargs):
    os.makedirs(output_path, exist_ok=True)
    log_file = os.path.join(output_path, "training_log.txt")
    if os.path.exists(log_file) and IS_MAIN_PROCESS: os.remove(log_file)

    print_and_log(f"Device: {device} | Epochs: {num_epochs} | Strategy: Just-In-Time RAM Loading", log_file)

    if depths is None: depths = [3, 3, 27, 3]
    if dims is None: dims = [192, 384, 768, 1536]

    # 1. Resolve Data Paths
    if isinstance(data_path, (list, tuple)) and len(data_path) == 2:
        train_roots, val_roots = data_path
    elif isinstance(data_path, list):
        train_roots = val_roots = data_path
    else:
        train_roots = val_roots = [data_path]

    # 2. Find all shard files
    train_shard_paths = _find_all_shards(train_roots, "train")
    if not train_shard_paths: raise ValueError("No shards found!")

    # 3. Build Validation Loader
    val_loader, genotype_map = get_data_loader(val_roots, "val", batch_size, num_workers, False, return_paths=True)
    num_classes = len(genotype_map)

    total_samples = len(train_shard_paths) * 4096
    total_batches = total_samples // batch_size
    print_and_log(f"Found {len(train_shard_paths)} shards. Est. batches: {total_batches}", log_file)

    # 4. Model Setup
    model = ConvNeXtCBAMClassifier(in_channels=6, class_num=num_classes, depths=depths, dims=dims).to(device)
    if ddp: model = DistributedDataParallel(model, device_ids=[local_rank], output_device=local_rank)
    model.apply(init_weights)

    # 5. Optim
    cw = torch.ones(num_classes, device=device);
    if num_classes > 1: cw[1] = float(pos_weight)

    if loss_type == 'combined':
        criterion = CombinedFocalWeightedCELoss(learning_rate, cw)
    else:
        criterion = nn.CrossEntropyLoss(weight=cw)

    optimizer = optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    scheduler = SequentialLR(optimizer,
                             schedulers=[LinearLR(optimizer, 0.01, total_iters=warmup_epochs),
                                         CosineAnnealingLR(optimizer, T_max=max(1, num_epochs - warmup_epochs))],
                             milestones=[warmup_epochs])

    start_epoch = 0
    best_f1 = -1.0
    if resume and os.path.isfile(resume):
        ckpt = torch.load(resume, map_location=device)
        model.module.load_state_dict(ckpt['model_state_dict']) if hasattr(model, "module") else model.load_state_dict(
            ckpt['model_state_dict'])
        start_epoch = int(ckpt.get('epoch', 0))
        best_f1 = float(ckpt.get('best_f1_true', -1.0))
        print_and_log(f"Resumed from epoch {start_epoch}", log_file)

    transform = _get_normalization_transform()

    # 6. Loop
    for epoch in range(start_epoch, num_epochs):
        model.train()
        running_loss = 0.0
        correct = 0
        total = 0
        batches_done = 0
        curr_lr = optimizer.param_groups[0]['lr']

        g = torch.Generator()
        g.manual_seed(epoch + 1000)
        perm = torch.randperm(len(train_shard_paths), generator=g).tolist()

        pbar = tqdm(total=total_batches, desc=f"Ep {epoch + 1}", disable=not IS_MAIN_PROCESS)

        for i, shard_idx in enumerate(perm):
            shard_path = train_shard_paths[shard_idx]

            try:
                ds = NpyDataset(root_dir=shard_path, transform=transform, load_to_ram=True)
                sampler = DistributedSampler(ds, shuffle=True, drop_last=False) if ddp else None
                if ddp: sampler.set_epoch(epoch)

                loader = DataLoader(ds, batch_size=batch_size, shuffle=(sampler is None),
                                    sampler=sampler, num_workers=num_workers,
                                    pin_memory=True, persistent_workers=(num_workers > 0),
                                    prefetch_factor=prefetch_factor if num_workers > 0 else None)

                for images, labels in loader:
                    images, labels = images.to(device, non_blocking=True), labels.to(device, non_blocking=True)
                    optimizer.zero_grad(set_to_none=True)

                    out = model(images)
                    if isinstance(out, tuple): out = out[0]

                    loss = criterion(out, labels, curr_lr) if loss_type == "combined" else criterion(out, labels)
                    loss.backward()
                    optimizer.step()

                    running_loss += loss.item()
                    _, pred = torch.max(out, 1)
                    correct += (pred == labels).sum().item()
                    total += labels.size(0)
                    batches_done += 1

                    if IS_MAIN_PROCESS: pbar.update(1)

                del loader
                del ds
            except Exception as e:
                if IS_MAIN_PROCESS: print(f"Skipping shard {os.path.basename(shard_path)}: {e}")
                continue

        if IS_MAIN_PROCESS: pbar.close()

        # --- Validation & Save ---
        world_size = dist.get_world_size() if ddp else 1
        val_loss, val_acc, _, _, metrics = evaluate_model(model, val_loader, criterion, genotype_map, loss_type,
                                                          curr_lr, ddp, world_size)
        f1 = metrics.get('f1_true', 0.0)

        tr_loss = running_loss / max(1, batches_done)
        tr_acc = correct / max(1, total) * 100
        print_and_log(
            f"Ep {epoch + 1} | Tr Loss {tr_loss:.4f} Acc {tr_acc:.1f}% | Val Loss {val_loss:.4f} Acc {val_acc:.1f}% F1 {f1:.4f}",
            log_file)

        if IS_MAIN_PROCESS:
            state = {'epoch': epoch + 1,
                     'model_state_dict': model.module.state_dict() if hasattr(model, 'module') else model.state_dict(),
                     'best_f1_true': best_f1}
            if (epoch + 1) == 5: torch.save(state, _unique_path(f"{output_path}/model_epoch_5.pth"))
            if f1 > best_f1:
                best_f1 = f1
                torch.save(state, _unique_path(f"{output_path}/model_best.pth"))
                print_and_log(f"Saved Best F1: {best_f1:.4f}", log_file)

        scheduler.step()


def evaluate_model(model, loader, criterion, map, loss_type, lr, ddp, world_size):
    model.eval()
    loss_sum, corr, tot = torch.zeros(3, device=device)
    tp = torch.zeros(len(map), device=device)
    fp = torch.zeros(len(map), device=device)
    fn = torch.zeros(len(map), device=device)

    with torch.no_grad():
        for batch in loader:
            imgs, lbls = batch[0].to(device), batch[1].to(device)
            out = model(imgs)
            if isinstance(out, tuple): out = out[0]
            loss_sum += (criterion(out, lbls, lr) if loss_type == "combined" else criterion(out, lbls))
            _, p = torch.max(out, 1)
            corr += (p == lbls).sum()
            tot += lbls.size(0)

            for i in range(lbls.size(0)):
                pi, ti = int(p[i]), int(lbls[i])
                if pi == ti:
                    tp[ti] += 1
                else:
                    if pi < len(map): fp[pi] += 1
                    fn[ti] += 1

    if ddp and world_size > 1:
        dist.all_reduce(loss_sum);
        dist.all_reduce(corr);
        dist.all_reduce(tot)
        dist.all_reduce(tp);
        dist.all_reduce(fp);
        dist.all_reduce(fn)

    avg_loss = loss_sum.item() / max(1, len(loader) * world_size)
    acc = corr.item() / max(1, tot.item()) * 100

    pos_idx = next((v for k, v in map.items() if str(k).lower() == "true"), 1)
    tpc, fpc, fnc = tp[pos_idx].item(), fp[pos_idx].item(), fn[pos_idx].item()
    prec = tpc / (tpc + fpc) if (tpc + fpc) > 0 else 0
    rec = tpc / (tpc + fnc) if (tpc + fnc) > 0 else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0

    return avg_loss, acc, {}, {}, {'f1_true': f1, 'recall_true': rec}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_paths", type=str, nargs='+')
    parser.add_argument("--train_data_paths_file", type=str)
    parser.add_argument("--val_data_paths_file", type=str)
    parser.add_argument("-o", "--output_path", default="./saved_models")
    parser.add_argument("--ddp", action="store_true")
    parser.add_argument("--local_rank", type=int, default=0)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=70)

    parser.add_argument("--depths", type=int, nargs='+', default=[3, 3, 27, 3])
    parser.add_argument("--dims", type=int, nargs='+', default=[192, 384, 768, 1536])
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--loss_type", default="weighted_ce")
    parser.add_argument("--pos_weight", type=float, default=88.0)
    parser.add_argument("--num_workers", type=int, default=8)
    parser.add_argument("--warmup_epochs", type=int, default=3)
    parser.add_argument("--weight_decay", type=float, default=0.01)
    parser.add_argument("--resume", type=str)
    parser.add_argument("--prefetch_factor", type=int, default=4)
    parser.add_argument("--mp_context", type=str, default=None)

    args, _ = parser.parse_known_args()

    if args.train_data_paths_file:
        data_in = (_read_paths_file(args.train_data_paths_file), _read_paths_file(args.val_data_paths_file))
    else:
        data_in = args.data_paths

    if args.ddp:
        # --- CRITICAL FIX FOR DDP ---
        # Read LOCAL_RANK from Environment Variable provided by torchrun
        local_rank = int(os.environ.get("LOCAL_RANK", 0))
        torch.cuda.set_device(local_rank)
        dist.init_process_group(backend="nccl")
        args.local_rank = local_rank  # Update args so train_model gets the right rank
        IS_MAIN_PROCESS = (dist.get_rank() == 0)
        if IS_MAIN_PROCESS: print(f"DDP Init: Size={dist.get_world_size()}")

    train_model(
        data_in, args.output_path,
        batch_size=args.batch_size,
        epochs=args.epochs,
        ddp=args.ddp,
        local_rank=args.local_rank,
        depths=args.depths,
        dims=args.dims,
        learning_rate=args.lr,
        loss_type=args.loss_type,
        pos_weight=args.pos_weight,
        num_workers=args.num_workers,
        warmup_epochs=args.warmup_epochs,
        weight_decay=args.weight_decay,
        resume=args.resume,
        prefetch_factor=args.prefetch_factor,
        mp_context=args.mp_context
    )

    if args.ddp: dist.destroy_process_group()