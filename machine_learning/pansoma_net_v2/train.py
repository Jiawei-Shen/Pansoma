"""Train PansomaNetV2 on merged Pansoma tensor sets.

    cd machine_learning
    python -m pansoma_net_v2.train --tensors <set> [<set> ...] --output <dir> [--val-chroms chr20] ...
    torchrun --nproc_per_node=N -m pansoma_net_v2.train --ddp ...           # several GPUs

A <set> is a directory holding merged SNV/ and INDEL/ (e.g. .../COLO829T_Illumina/v3_tensors). Classes are
the labels 0 non, 1 somatic, 2 germline; tensors labelled -1 are not used. Chromosomes split train / val /
test (test chromosomes are left out entirely). The encoder's z-score statistics are fitted once on the
training tensors (or taken from --resume) and saved in every checkpoint. The best checkpoint is chosen by
the validation F1 of the somatic class.
"""
import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import torch
import torch.distributed as dist
import torch.nn as nn
from torch.nn.parallel import DistributedDataParallel
from torch.utils.data import DataLoader, Subset

from .data import CLASSES, KINDS, EpochSampler, TensorDataset, load_parts
from .encode import PLANES, compute_stats
from .model import PansomaNetV2

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--tensors", nargs="+", required=True, help="merged tensor sets (directories with SNV/ INDEL/)")
    p.add_argument("--kinds", nargs="+", default=list(KINDS), choices=KINDS)
    p.add_argument("--output", required=True)
    p.add_argument("--cache-dir", help="index cache (default: <output>/index_cache)")
    p.add_argument("--train-chroms", nargs="+", help="default: every autosome not in --val-chroms/--test-chroms")
    p.add_argument("--val-chroms", nargs="+", default=["chr20"])
    p.add_argument("--test-chroms", nargs="+", default=[], help="left out of training and validation")
    p.add_argument("--epochs", type=int, default=30)
    p.add_argument("--epoch-samples", type=int, help="training tensors per epoch (a fresh random subset each epoch)")
    p.add_argument("--val-samples", type=int, help="fixed random subset of the validation tensors")
    p.add_argument("--batch-size", type=int, default=64, help="per GPU")
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--weight-decay", type=float, default=0.05)
    p.add_argument("--warmup-epochs", type=float, default=1.0)
    p.add_argument("--class-weights", default="balanced",
                   help="'balanced' (n / (3 n_c) on the training tensors), 'none', or three numbers w_non,w_som,w_germ")
    p.add_argument("--depths", type=int, nargs=4, default=[3, 3, 27, 3])
    p.add_argument("--dims", type=int, nargs=4, default=[192, 384, 768, 1536])
    p.add_argument("--front", type=int, nargs="+", default=[64, 64], help="widths of the 1x1 layers before the stem")
    p.add_argument("--drop-path", type=float, default=0.1)
    p.add_argument("--stats-samples", type=int, default=20000, help="training tensors used to fit the z-score")
    p.add_argument("--amp", choices=["bf16", "off"], default="bf16")
    p.add_argument("--num-workers", type=int, default=8)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--ddp", action="store_true")
    p.add_argument("--resume", help="checkpoint to continue from (its statistics are kept)")
    p.add_argument("--rebuild-index", action="store_true")
    return p.parse_args(argv)


class Run:
    """Device, DDP rank and rank-0 logging."""

    def __init__(self, args):
        self.ddp = args.ddp
        if self.ddp:
            dist.init_process_group(backend="nccl")
            self.rank, self.world = dist.get_rank(), dist.get_world_size()
            self.local = int(os.environ.get("LOCAL_RANK", 0))
            torch.cuda.set_device(self.local)
            self.device = torch.device("cuda", self.local)
        else:
            self.rank, self.world, self.local = 0, 1, 0
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.main = self.rank == 0
        self.out = Path(args.output)
        if self.main:
            self.out.mkdir(parents=True, exist_ok=True)

    def log(self, message):
        if self.main:
            print(message, flush=True)
            with open(self.out / "train.log", "a") as f:
                f.write(message + "\n")

    def record(self, row):
        if self.main:
            with open(self.out / "metrics.jsonl", "a") as f:
                f.write(json.dumps(row) + "\n")

    def barrier(self):
        if self.ddp:
            dist.barrier()


def split_chroms(args, available):
    val, test = set(args.val_chroms), set(args.test_chroms)
    train = set(args.train_chroms) if args.train_chroms else {c for c in AUTOSOMES if c not in val | test}
    if train & (val | test) or val & test:
        raise SystemExit(f"chromosome sets overlap: train {sorted(train & (val | test))}, val/test {sorted(val & test)}")
    missing = (train | val | test) - set(available)
    return sorted(train), sorted(val), sorted(test), sorted(missing)


def class_weights(spec, counts):
    if spec == "none":
        return torch.ones(len(CLASSES))
    if spec == "balanced":
        counts = np.maximum(counts, 1)
        return torch.tensor(counts.sum() / (len(CLASSES) * counts), dtype=torch.float32)
    w = [float(x) for x in spec.split(",")]
    if len(w) != len(CLASSES):
        raise SystemExit("--class-weights needs three numbers: non,somatic,germline")
    return torch.tensor(w)


def metrics_from_confusion(cm):
    """cm[true, pred] -> accuracy and per-class precision / recall / F1."""
    cm = cm.astype(np.float64)
    out = {"accuracy": float(np.trace(cm) / max(cm.sum(), 1))}
    for k, name in enumerate(CLASSES):
        tp, fp, fn = cm[k, k], cm[:, k].sum() - cm[k, k], cm[k, :].sum() - cm[k, k]
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        out[name] = {"precision": prec, "recall": rec, "f1": 2 * prec * rec / (prec + rec) if prec + rec else 0.0,
                     "support": int(cm[k, :].sum())}
    return out


def evaluate(model, loader, criterion, run, amp):
    model.eval()
    cm = torch.zeros(len(CLASSES), len(CLASSES), dtype=torch.long, device=run.device)
    loss_sum = torch.zeros(2, dtype=torch.float64, device=run.device)  # weighted loss sum, weight sum
    with torch.no_grad():
        for x, blocks, y in loader:
            x, blocks, y = (t.to(run.device, non_blocking=True) for t in (x, blocks, y))
            with torch.autocast(run.device.type, dtype=torch.bfloat16, enabled=amp):
                logits = model(x, blocks)
            w = criterion.weight[y].sum() if criterion.weight is not None else \
                torch.tensor(float(y.numel()), device=run.device)
            loss_sum += torch.stack([criterion(logits.float(), y) * w, w]).double()
            pred = logits.argmax(1)
            cm += torch.bincount(y * len(CLASSES) + pred, minlength=len(CLASSES) ** 2).view(len(CLASSES), -1)
    if run.ddp:
        dist.all_reduce(cm)
        dist.all_reduce(loss_sum)
    m = metrics_from_confusion(cm.cpu().numpy())
    m["loss"] = float(loss_sum[0] / max(loss_sum[1], 1))
    m["confusion"] = cm.cpu().tolist()
    return m


def make_loader(dataset, sampler, num_workers, batch_size, persistent=True):
    """Workers come from a forkserver, not fork: a forked worker's RSS includes every page it shares with the
    main process (CUDA context, pinned memory), and Slurm limits the RSS summed over processes."""
    workers = dict(num_workers=num_workers, persistent_workers=persistent, prefetch_factor=4,
                   multiprocessing_context="forkserver") if num_workers > 0 else {}
    return DataLoader(dataset, batch_size=batch_size, sampler=sampler, pin_memory=torch.cuda.is_available(),
                      drop_last=False, **workers)


def main(argv=None):
    args = parse_args(argv)
    run = Run(args)
    torch.manual_seed(args.seed)
    cache = args.cache_dir or str(Path(args.output) / "index_cache")
    if run.main:  # rank 0 builds / refreshes the index cache, the others read it
        load_parts(args.tensors, args.kinds, cache, rebuild=args.rebuild_index)
    run.barrier()
    everything = load_parts(args.tensors, args.kinds, cache)
    available = sorted({c for index, _ in everything for c in index.meta["chroms"]})
    train_chroms, val_chroms, test_chroms, missing = split_chroms(args, available)
    run.log(f"train {train_chroms}\nval {val_chroms}\ntest (left out) {test_chroms}"
            + (f"\nnot in the data: {missing}" if missing else ""))
    train_parts = [(i, i.select(train_chroms)) for i, _ in everything]
    val_parts = [(i, i.select(val_chroms)) for i, _ in everything]
    train_set, val_set = TensorDataset(train_parts), TensorDataset(val_parts)
    counts = train_set.class_counts()
    for index, pos in train_parts:
        labels = index.arrays["label"][pos]
        run.log(f"{index.dir}: labels {index.meta['labels']['version']} ({index.meta['labels']['created']}), "
                f"train {len(pos):,} " + " ".join(f"{c}={int((labels == k).sum()):,}" for k, c in enumerate(CLASSES)))
    run.log(f"training tensors {len(train_set):,} ({dict(zip(CLASSES, counts.tolist()))}), validation {len(val_set):,}")
    if len(train_set) == 0 or len(val_set) == 0:
        raise SystemExit("empty training or validation set")
    if args.val_samples and args.val_samples < len(val_set):
        keep = np.sort(np.random.default_rng(args.seed).choice(len(val_set), args.val_samples, replace=False))
        val_set = Subset(val_set, keep.tolist())

    model = PansomaNetV2(len(CLASSES), args.depths, args.dims, args.front, args.drop_path)
    checkpoint = torch.load(args.resume, map_location="cpu", weights_only=False) if args.resume else None
    if checkpoint is not None:
        if checkpoint["config"] != model.config:
            raise SystemExit(f"--resume model config {checkpoint['config']} differs from the arguments' {model.config}")
        model.load_state_dict(checkpoint["model_state_dict"])
        stats = checkpoint["stats"]
        run.log(f"resumed from {args.resume} (epoch {checkpoint['epoch']}); statistics kept: {json.dumps(stats)}")
    else:
        stats = [None]
        if run.main:
            order = EpochSampler(len(train_set), args.stats_samples, shuffle=True, seed=args.seed + 12345)
            t0 = time.time()
            stats[0] = compute_stats(make_loader(train_set, order, args.num_workers, 256, persistent=False),
                                     args.stats_samples)
            run.log(f"statistics from {stats[0]['tensors']:,} training tensors ({time.time() - t0:.0f} s): "
                    + json.dumps({k: v for k, v in stats[0].items() if k != "tensors"}))
            (run.out / "stats.json").write_text(json.dumps(stats[0], indent=2) + "\n")
        if run.ddp:
            dist.broadcast_object_list(stats, src=0)
        stats = stats[0]
    model.encoder.set_stats(stats)
    model.to(run.device)
    net = DistributedDataParallel(model, device_ids=[run.local]) if run.ddp else model

    weights = class_weights(args.class_weights, counts).to(run.device)
    run.log(f"class weights {dict(zip(CLASSES, [round(w, 4) for w in weights.tolist()]))}")
    criterion = nn.CrossEntropyLoss(weight=weights, ignore_index=-1)
    train_sampler = EpochSampler(len(train_set), args.epoch_samples, True, run.rank, run.world, args.seed)
    val_sampler = EpochSampler(len(val_set), None, False, run.rank, run.world)
    train_loader = make_loader(train_set, train_sampler, args.num_workers, args.batch_size)
    val_loader = make_loader(val_set, val_sampler, args.num_workers, args.batch_size, persistent=False)
    steps_per_epoch = max(1, len(train_loader))
    total_steps, warmup = args.epochs * steps_per_epoch, int(args.warmup_epochs * steps_per_epoch)
    optimizer = torch.optim.AdamW(net.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    schedule = lambda s: (s + 1) / max(1, warmup) if s < warmup else \
        0.5 * (1 + math.cos(math.pi * (s - warmup) / max(1, total_steps - warmup)))  # noqa: E731
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, schedule)
    start_epoch, best = 0, {"somatic_f1": -1.0, "loss": math.inf, "epoch": None}
    if checkpoint is not None:
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
        start_epoch, best = checkpoint["epoch"], checkpoint.get("best", best)
    amp = args.amp == "bf16" and run.device.type == "cuda"
    n_params = sum(p.numel() for p in model.parameters())
    run.log(f"model {model.config}, {n_params / 1e6:.1f} M parameters, input planes {len(PLANES)}; "
            f"device {run.device} x {run.world}, amp {amp}; {steps_per_epoch} steps per epoch")

    def save(path, epoch, val):
        if run.main:
            payload = dict(format="pansoma_net_v2", config=model.config, model_state_dict=model.state_dict(),
                           optimizer_state_dict=optimizer.state_dict(), scheduler_state_dict=scheduler.state_dict(),
                           epoch=epoch, best=best, val=val, classes=list(CLASSES), planes=list(PLANES),
                           stats=model.encoder.stats(), args=vars(args),
                           data=[dict(directory=str(i.dir), labels=i.meta["labels"], train=int(len(p)))
                                 for i, p in train_parts],
                           chroms=dict(train=train_chroms, val=val_chroms, test=test_chroms))
            tmp = run.out / (path + ".tmp")
            torch.save(payload, tmp)
            tmp.replace(run.out / path)

    for epoch in range(start_epoch, args.epochs):
        net.train()
        train_sampler.set_epoch(epoch)
        t0, loss_sum, seen, correct = time.time(), 0.0, 0, 0
        for step, (x, blocks, y) in enumerate(train_loader):
            x, blocks, y = (t.to(run.device, non_blocking=True) for t in (x, blocks, y))
            with torch.autocast(run.device.type, dtype=torch.bfloat16, enabled=amp):
                logits = net(x, blocks)
            loss = criterion(logits.float(), y)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            scheduler.step()
            loss_sum += loss.item() * y.numel()
            seen += y.numel()
            correct += int((logits.argmax(1) == y).sum())
            if run.main and (step + 1) % 200 == 0:
                print(f"epoch {epoch + 1} step {step + 1}/{steps_per_epoch} loss {loss_sum / seen:.4f} "
                      f"lr {scheduler.get_last_lr()[0]:.2e} {seen * run.world / (time.time() - t0):.0f} tensors/s",
                      flush=True)
        train_time = time.time() - t0
        val = evaluate(net, val_loader, criterion, run, amp)
        som = val["somatic"]
        improved = som["f1"] > best["somatic_f1"] or (som["f1"] == best["somatic_f1"] and val["loss"] < best["loss"])
        if improved:
            best = {"somatic_f1": som["f1"], "loss": val["loss"], "epoch": epoch + 1}
        row = dict(epoch=epoch + 1, train_loss=loss_sum / max(seen, 1), train_accuracy=correct / max(seen, 1),
                   train_seconds=round(train_time, 1), tensors_per_second=round(seen * run.world / train_time, 1),
                   lr=scheduler.get_last_lr()[0], val=val)
        run.record(row)
        run.log(f"epoch {epoch + 1}/{args.epochs}: train loss {row['train_loss']:.4f} acc {row['train_accuracy']:.4f} "
                f"({train_time:.0f} s) | val loss {val['loss']:.4f} acc {val['accuracy']:.4f} | "
                + " | ".join(f"{c} P {val[c]['precision']:.3f} R {val[c]['recall']:.3f} F1 {val[c]['f1']:.3f}"
                             for c in CLASSES) + (" | best" if improved else ""))
        save("last.pth", epoch + 1, val)
        if improved:
            save("best.pth", epoch + 1, val)
    run.log(f"done; best epoch {best['epoch']} somatic F1 {best['somatic_f1']:.4f} -> {run.out / 'best.pth'}")
    if run.ddp:
        dist.destroy_process_group()


if __name__ == "__main__":
    main()
