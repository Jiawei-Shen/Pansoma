"""Train PansomaNetV2 on merged Pansoma tensor sets.

    cd machine_learning
    python -m pansoma_net_v2.train --tensors <set> [<set> ...] --kinds SNV --output <dir> [--test-chroms chr1] ...
    torchrun --nproc_per_node=N -m pansoma_net_v2.train --ddp ...           # several GPUs

A <set> is a directory holding merged SNV/ and INDEL/ (e.g. .../COLO829T_Illumina/v3_tensors). Classes are
the labels 0 non, 1 somatic, 2 germline; tensors labelled -1 are not used. The test chromosomes (default chr1)
are left out; the others train, except the validation: whole ~1 Mb node blocks (--val-fraction of them), or
--val-chroms. The encoder's z-score statistics (and with --scalars the scalars') are fitted once on the
training tensors (or taken from --resume) and saved in every checkpoint. Each epoch the somatic threshold of
the best validation F1 is found (metrics.py); the best checkpoint has the highest thresholded somatic F1 (or
--select ap: somatic average precision) and stores its threshold.
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

from . import metrics
from .env import triton_libcuda
from .data import (CLASSES, KINDS, SCALARS, EpochSampler, TensorDataset, block_split, load_parts, somatic_truth,
                   validation_truth)

AF = SCALARS.index("af")
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
    p.add_argument("--test-chroms", nargs="+", default=["chr1"], help="left out of training and validation")
    p.add_argument("--val-chroms", nargs="+", default=[], help="validation chromosomes (default: node blocks)")
    p.add_argument("--val-fraction", type=float, default=0.05, help="node blocks of the training chromosomes")
    p.add_argument("--val-block-nodes", type=int, default=20000)
    p.add_argument("--select", choices=["f1", "ap", "truth_f1"], default="f1",
                   help="best checkpoint by the tensor-level somatic F1 at its best threshold (default), the somatic AP, "
                        "or the truth-level F1 at its best threshold (metrics.truth_report); the checkpoint stores "
                        "that criterion's threshold")
    p.add_argument("--scalars", action="store_true", help="feed the site scalars (data.SCALARS) to the head")
    p.add_argument("--epochs", type=int, default=30)
    p.add_argument("--epoch-samples", type=int, help="training tensors per epoch (a fresh random subset each epoch)")
    p.add_argument("--val-samples", type=int, help="fixed random subset of the validation tensors")
    p.add_argument("--batch-size", type=int, default=64, help="per GPU")
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--weight-decay", type=float, default=0.05)
    p.add_argument("--warmup-epochs", type=float, default=1.0)
    p.add_argument("--class-weights", default="balanced",
                   help="'balanced' (n / (3 n_c) over an epoch's tensors), 'sqrt' (its square root), 'none', or three "
                        "numbers w_non,w_som,w_germ")
    p.add_argument("--non-fraction", type=float, default=1.0,
                   help="each epoch: every somatic / germline tensor and a fresh random fraction of the non ones")
    p.add_argument("--keep-non-af", type=float,
                   help="non tensors with AF >= this are taken every epoch; the sampled rest weigh 1/--non-fraction "
                        "in the loss and the class weights use all tensors, so the expected loss is that of all non "
                        "(without the weights an AF-dependent sampling would bias p_somatic by AF)")
    p.add_argument("--depths", type=int, nargs=4, default=[3, 3, 27, 3])
    p.add_argument("--dims", type=int, nargs=4, default=[192, 384, 768, 1536])
    p.add_argument("--front", type=int, nargs="+", default=[64, 64], help="widths of the 1x1 layers before the stem")
    p.add_argument("--drop-path", type=float, default=0.1)
    p.add_argument("--stats-samples", type=int, default=20000, help="training tensors used to fit the z-score")
    p.add_argument("--amp", choices=["bf16", "off"], default="bf16")
    p.add_argument("--no-compile", dest="compile", action="store_false",
                   help="skip torch.compile (default on CUDA: 3x faster on an H100, ~2 min to compile)")
    p.add_argument("--no-channels-last", dest="channels_last", action="store_false")
    p.add_argument("--num-workers", type=int, default=8,
                   help="training DataLoader workers: reads from BeeGFS are latency-bound, so more workers read "
                        "faster (HG008 Illumina SNV: 351 / 826 / 2,290 tensors/s with 14 / 28 / 48)")
    p.add_argument("--val-workers", type=int, default=8)
    p.add_argument("--prefetch-factor", type=int, default=2, help="batches in flight per worker")
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
            self.local = int(os.environ.get("LOCAL_RANK", 0))
            torch.cuda.set_device(self.local)
            self.device = torch.device("cuda", self.local)
            dist.init_process_group(backend="nccl", device_id=self.device)
            self.rank, self.world = dist.get_rank(), dist.get_world_size()
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
    """counts: tensors per class in one epoch."""
    if spec == "none":
        return torch.ones(len(CLASSES))
    if spec in ("balanced", "sqrt"):
        counts = np.maximum(counts, 1)
        w = counts.sum() / (len(CLASSES) * counts)
        return torch.tensor(np.sqrt(w) if spec == "sqrt" else w, dtype=torch.float32)
    w = [float(x) for x in spec.split(",")]
    if len(w) != len(CLASSES):
        raise SystemExit("--class-weights needs three numbers: non,somatic,germline")
    return torch.tensor(w)


def predict_probs(model, loader, device, amp):
    """(labels, probabilities) of loader's samples on this rank."""
    model.eval()
    labels, probs = [], []
    with torch.no_grad():
        for x, blocks, scalars, y in loader:
            with torch.autocast(device.type, dtype=torch.bfloat16, enabled=amp):
                logits = model(x.to(device, non_blocking=True), blocks.to(device, non_blocking=True),
                               scalars.to(device, non_blocking=True))
            probs.append(torch.softmax(logits.float(), 1).cpu().numpy())
            labels.append(y.numpy())
    return (np.concatenate(labels) if labels else np.zeros(0, np.int64),
            np.concatenate(probs) if probs else np.zeros((0, len(CLASSES)), np.float32))


def evaluate(model, loader, criterion, run, amp, truth=None):
    """Validation report (metrics.report at the best tensor-level somatic threshold), weighted loss and, with
    truth = (truth keys, per-sample matched keys), metrics.truth_report at the same threshold; all ranks."""
    labels, probs = predict_probs(model, loader, run.device, amp)
    order = np.asarray(list(loader.sampler), np.int64)[:len(labels)]  # dataset indices of this rank's samples
    if run.ddp:
        parts = [None] * run.world
        dist.all_gather_object(parts, (order, labels, probs))
        order, labels, probs = (np.concatenate([p[k] for p in parts]) for k in range(3))
    s = np.argsort(order, kind="stable")
    order, labels, probs = order[s], labels[s], probs[s]
    threshold, _, _, _ = metrics.best_somatic_threshold(labels, probs)
    out = metrics.report(labels, probs, threshold)
    if truth is not None:
        keys, matches = truth
        out["truth"] = metrics.truth_report(keys, [matches[i] for i in order], labels, probs, threshold)
    w = criterion.weight.cpu().numpy() if criterion.weight is not None else np.ones(len(CLASSES))
    nll = -np.log(np.clip(probs[np.arange(len(labels)), labels], 1e-12, None))
    out["loss"] = float((w[labels] * nll).sum() / max(w[labels].sum(), 1e-12))
    return out


def describe(report):
    """One line: somatic AP, F1 / P / R at the report's threshold, argmax F1; germline AP and argmax F1."""
    t, a = report.get("thresholded"), report["argmax"]
    line = f"somatic AP {report['somatic_ap']:.3f}"
    if t:
        s = t["somatic"]
        line += f" F1 {s['f1']:.3f} (P {s['precision']:.3f} R {s['recall']:.3f} @ p>={report['threshold']:.3f})"
    line += f" argmax F1 {a['somatic']['f1']:.3f}"
    if report.get("truth"):
        tr, at = report["truth"], report["truth"]["at_threshold"]
        line += (f" | truth F1 {at['f1']:.3f} (P {at['precision']:.3f} R {at['recall']:.3f}; {tr['truth_alleles']:,} "
                 f"truth, ceiling {tr['ceiling']:.3f}, best F1 {tr['best']['f1']:.3f} @ p>={tr['best']['threshold']:.3f})")
    return line + (f" | germline AP {report['germline_ap']:.3f} F1 {a['germline']['f1']:.3f} | "
                   f"non F1 {a['non']['f1']:.3f}")


def gpu_peak(run):
    """Peak GPU memory (GiB) of this epoch, the maximum over ranks: allocated by tensors, reserved by the caching
    allocator (what nvidia-smi shows, minus the CUDA context)."""
    if run.device.type != "cuda":
        return None
    peak = torch.tensor([torch.cuda.max_memory_allocated(run.device), torch.cuda.max_memory_reserved(run.device)],
                        dtype=torch.float64, device=run.device) / 2 ** 30
    if run.ddp:
        dist.all_reduce(peak, op=dist.ReduceOp.MAX)
    return {"allocated": round(float(peak[0]), 2), "reserved": round(float(peak[1]), 2)}


def make_loader(dataset, sampler, num_workers, batch_size, persistent=True, prefetch_factor=2):
    """Workers come from a forkserver, not fork: a forked worker's RSS includes every page it shares with the
    main process (CUDA context, pinned memory), and Slurm limits the RSS summed over processes."""
    workers = dict(num_workers=num_workers, persistent_workers=persistent, prefetch_factor=prefetch_factor,
                   multiprocessing_context="forkserver") if num_workers > 0 else {}
    return DataLoader(dataset, batch_size=batch_size, sampler=sampler, pin_memory=torch.cuda.is_available(),
                      drop_last=False, **workers)


WORKER_START_FAILURES = (RuntimeError, FileNotFoundError, EOFError, ConnectionError)


def retry_workers(action, log, what, tries=3, wait=5.0):
    """action() with fresh DataLoader workers, again (up to `tries`) when they fail to start. A forkserver
    worker opens the parent's semaphores by name in /dev/shm; on tequila three of four jobs started at once
    lost those names (FileNotFoundError in SemLock._rebuild, "DataLoader worker exited unexpectedly"). Once
    started, workers hold the semaphores, so the loaders are persistent and only their start is retried."""
    for attempt in range(1, tries + 1):
        try:
            return action()
        except WORKER_START_FAILURES as e:
            if attempt == tries:
                raise
            log(f"{what}: DataLoader workers failed to start ({type(e).__name__}: {str(e)[:160]}); retry {attempt}")
            time.sleep(wait)


def started_loader(dataset, sampler, num_workers, batch_size, log, what, prefetch_factor=2):
    """A persistent DataLoader whose workers are running (one batch fetched), retried as in retry_workers."""
    def start():
        loader = make_loader(dataset, sampler, num_workers, batch_size, persistent=True,
                             prefetch_factor=prefetch_factor)
        if num_workers > 0 and len(loader):
            next(iter(loader))
        return loader
    return retry_workers(start, log, what)


def main(argv=None):
    args = parse_args(argv)
    run = Run(args)
    if run.main:
        (run.out / "args.json").write_text(json.dumps(vars(args), indent=2) + "\n")
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
    # training: labels >= 0; validation (like the test): evaluation labels (off-reference tensors without a
    # truth match count as non), since the threshold chosen on it is applied to the test chromosome
    if val_chroms:
        train_parts = [(i, i.select(train_chroms)) for i, _ in everything]
        val_parts = [(i, i.select(val_chroms, evaluation=True)) for i, _ in everything]
    else:
        train_parts, val_parts = [], []
        for i, _ in everything:
            tr, va = block_split(i, i.select(train_chroms, labelled=False), args.val_fraction, args.val_block_nodes,
                                 args.seed)
            train_parts.append((i, tr[i.arrays["label"][tr] >= 0]))
            val_parts.append((i, va[i.arrays["eval_label"][va] >= 0]))
        run.log(f"validation: {args.val_fraction:.0%} of the {args.val_block_nodes}-node blocks of the training chromosomes")
    train_set, val_set = TensorDataset(train_parts), TensorDataset(val_parts, labels="eval_label")
    n_off = sum(int((i.arrays["eval_label"][v] != i.arrays["label"][v]).sum()) for i, v in val_parts)
    run.log(f"validation includes {n_off:,} off-reference tensors without a truth match, as non")
    counts = train_set.class_counts()
    for index, pos in train_parts:
        labels = index.arrays["label"][pos]
        lm = index.meta["labels"]
        run.log(f"{index.dir}: labels {lm.get('version') or (lm.get('format'), (lm.get('rules_sha256') or '')[:12])} "
                f"({lm['created']}), "
                f"train {len(pos):,} " + " ".join(f"{c}={int((labels == k).sum()):,}" for k, c in enumerate(CLASSES)))
    run.log(f"training tensors {len(train_set):,} ({dict(zip(CLASSES, counts.tolist()))}), validation {len(val_set):,} "
            f"({dict(zip(CLASSES, val_set.class_counts().tolist()))})")
    if len(train_set) == 0 or len(val_set) == 0:
        raise SystemExit("empty training or validation set")
    keep = None
    if args.val_samples and args.val_samples < len(val_set):
        keep = np.sort(np.random.default_rng(args.seed).choice(len(val_set), args.val_samples, replace=False))
    # truth-level validation: the truth alleles of the validation region (with or without a tensor) and, per
    # validation sample, the truth alleles its tensor stands for; keys are (part, truth id)
    val_truth, truth_keys = None, set()
    for p, ((index, _), (_, va)) in enumerate(zip(train_parts, val_parts)):
        t_all = somatic_truth(index, chroms=val_chroms or train_chroms)
        if t_all is None:
            run.log(f"{index.dir}: no somatic.recall.tsv, no truth-level validation")
            truth_keys = None
            break
        vt = set(t_all) if val_chroms else validation_truth(
            index, index.select(train_chroms, labelled=False), t_all, args.val_fraction, args.val_block_nodes, args.seed)
        truth_keys |= {(p, t) for t in vt}
    if truth_keys is not None:
        ids = keep.tolist() if keep is not None else range(len(val_set))
        matches = [{(int(val_set.part[i]), int(t)) for t in val_parts[val_set.part[i]][0].truth_of(val_set.position[i])}
                   for i in ids]
        val_truth = (truth_keys, matches)
        run.log(f"truth-level validation: {len(truth_keys):,} somatic truth alleles in the validation region "
                f"({len(set().union(*matches) & truth_keys) if matches else 0:,} with a validation tensor)")
    if keep is not None:
        val_set = Subset(val_set, keep.tolist())

    model = PansomaNetV2(len(CLASSES), args.depths, args.dims, args.front, args.drop_path,
                         scalars=len(SCALARS) if args.scalars else 0)
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
            stats[0] = retry_workers(lambda: compute_stats(make_loader(train_set, order, args.num_workers, 256,
                                                                       persistent=False), args.stats_samples),
                                     run.log, "statistics")
            run.log(f"statistics from {stats[0]['tensors']:,} training tensors ({time.time() - t0:.0f} s): "
                    + json.dumps({k: v for k, v in stats[0].items() if k != "tensors"}))
            (run.out / "stats.json").write_text(json.dumps(stats[0], indent=2) + "\n")
        if run.ddp:
            dist.broadcast_object_list(stats, src=0)
        stats = stats[0]
    model.encoder.set_stats(stats)
    if args.scalars and checkpoint is None:  # all training tensors (the arrays are in memory)
        s = train_set.scalars.astype(np.float64)
        model.set_scalar_stats(s.mean(0), s.std(0))
        run.log("scalar statistics: " + json.dumps({n: [round(float(m), 4), round(float(d), 4)]
                                                   for n, m, d in zip(SCALARS, s.mean(0), s.std(0))}))
    model.to(run.device)
    if args.channels_last and run.device.type == "cuda":
        model.to(memory_format=torch.channels_last)
    net = DistributedDataParallel(model, device_ids=[run.local]) if run.ddp else model
    if args.compile and run.device.type == "cuda":
        run.log(f"torch.compile (Triton libcuda directory: {triton_libcuda()})")
        net = torch.compile(net)  # after DDP, as PyTorch recommends; the checkpoint saves `model`

    # class weights follow one epoch's classes; with --keep-non-af the importance weights restore all non
    importance = args.keep_non_af is not None and args.non_fraction < 1.0
    epoch_counts = counts.astype(np.float64) * np.array([1.0 if importance else args.non_fraction, 1.0, 1.0])
    weights = class_weights(args.class_weights, epoch_counts).to(run.device)
    criterion = nn.CrossEntropyLoss(weight=weights, ignore_index=-1)
    per_sample = nn.CrossEntropyLoss(weight=weights, ignore_index=-1, reduction="none")
    keep = train_set.scalars[:, AF] >= args.keep_non_af if args.keep_non_af is not None else None
    train_sampler = EpochSampler(len(train_set), args.epoch_samples, True, run.rank, run.world, args.seed,
                                 labels=train_set.label, non_fraction=args.non_fraction, keep=keep)
    if importance:
        n_keep = int((keep & (train_set.label == 0)).sum())
        run.log(f"non tensors with AF >= {args.keep_non_af}: {n_keep:,} taken every epoch; the others sampled at "
                f"{args.non_fraction} with importance weight {1 / args.non_fraction:.2f}")
    run.log(f"class weights {dict(zip(CLASSES, [round(w, 4) for w in weights.tolist()]))}; each epoch "
            f"{train_sampler.epoch_size():,} tensors (non fraction {args.non_fraction})")
    val_sampler = EpochSampler(len(val_set), None, False, run.rank, run.world)
    train_loader = started_loader(train_set, train_sampler, args.num_workers, args.batch_size, run.log, "training",
                                  args.prefetch_factor)
    val_loader = started_loader(val_set, val_sampler, args.val_workers if args.num_workers else 0, args.batch_size,
                                run.log, "validation", args.prefetch_factor)
    steps_per_epoch = max(1, len(train_loader))
    total_steps, warmup = args.epochs * steps_per_epoch, int(args.warmup_epochs * steps_per_epoch)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay,
                                  fused=run.device.type == "cuda")
    schedule = lambda s: (s + 1) / max(1, warmup) if s < warmup else \
        0.5 * (1 + math.cos(math.pi * (s - warmup) / max(1, total_steps - warmup)))  # noqa: E731
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, schedule)
    start_epoch, best = 0, {"score": -1.0, "loss": math.inf, "epoch": None, "threshold": None}
    if checkpoint is not None:
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
        start_epoch, best = checkpoint["epoch"], checkpoint.get("best", best)
    amp = args.amp == "bf16" and run.device.type == "cuda"
    n_params = sum(p.numel() for p in model.parameters())
    run.log(f"model {model.config}, {n_params / 1e6:.1f} M parameters, input planes {len(PLANES)}; "
            f"device {run.device} x {run.world}, amp {amp}, compile {args.compile and run.device.type == 'cuda'}, "
            f"channels_last {args.channels_last and run.device.type == 'cuda'}; {steps_per_epoch} steps per epoch")

    def save(path, epoch, val, threshold):
        if run.main:
            payload = dict(format="pansoma_net_v2", config=model.config, model_state_dict=model.state_dict(),
                           optimizer_state_dict=optimizer.state_dict(), scheduler_state_dict=scheduler.state_dict(),
                           epoch=epoch, best=best, val=val, somatic_threshold=threshold, classes=list(CLASSES),
                           planes=list(PLANES), scalars=list(SCALARS) if args.scalars else [],
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
        if run.device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(run.device)
        t0, loss_sum, seen, correct = time.time(), 0.0, 0, 0
        for step, (x, blocks, scalars, y) in enumerate(train_loader):
            x, blocks, scalars, y = (t.to(run.device, non_blocking=True) for t in (x, blocks, scalars, y))
            with torch.autocast(run.device.type, dtype=torch.bfloat16, enabled=amp):
                logits = net(x, blocks, scalars)
            if importance:  # a sampled low-AF non tensor stands for 1 / non_fraction of them
                iw = torch.where((y == 0) & (scalars[:, AF] < args.keep_non_af), 1.0 / args.non_fraction, 1.0)
                loss = (per_sample(logits.float(), y) * iw).sum() / (weights[y] * iw).sum()
            else:
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
        val = evaluate(net, val_loader, criterion, run, amp, val_truth)
        if args.select == "truth_f1" and val.get("truth"):
            score, threshold = val["truth"]["best"]["f1"], val["truth"]["best"]["threshold"]
        else:
            score = val["thresholded"]["somatic"]["f1"] if args.select == "f1" else val["somatic_ap"]
            threshold = val["threshold"]
        improved = score > best["score"] or (score == best["score"] and val["loss"] < best["loss"])
        if improved:
            best = {"score": score, "select": args.select, "loss": val["loss"], "epoch": epoch + 1,
                    "threshold": threshold}
        gpu = gpu_peak(run)
        row = dict(epoch=epoch + 1, train_loss=loss_sum / max(seen, 1), train_accuracy=correct / max(seen, 1),
                   train_seconds=round(train_time, 1), tensors_per_second=round(seen * run.world / train_time, 1),
                   lr=scheduler.get_last_lr()[0], gpu_peak_gib=gpu, val=val)
        run.record(row)
        run.log(f"epoch {epoch + 1}/{args.epochs}: train loss {row['train_loss']:.4f} acc {row['train_accuracy']:.4f} "
                f"({train_time:.0f} s, {row['tensors_per_second']:.0f} tensors/s"
                + (f", GPU peak {gpu['allocated']:.1f} / {gpu['reserved']:.1f} GiB allocated / reserved" if gpu else "")
                + f") | val loss {val['loss']:.4f} | " + describe(val) + (" | best" if improved else ""))
        save("last.pth", epoch + 1, val, threshold)
        if improved:
            save("best.pth", epoch + 1, val, threshold)
    run.log(f"done; best epoch {best['epoch']} ({best.get('select', args.select)} {best['score']:.4f}, "
            f"somatic threshold {best['threshold']}) -> {run.out / 'best.pth'}")
    if run.ddp:
        dist.destroy_process_group()


if __name__ == "__main__":
    main()
