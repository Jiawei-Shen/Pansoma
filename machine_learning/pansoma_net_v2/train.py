"""Train PansomaNetV2 on merged Pansoma tensor sets.

    cd machine_learning
    python -m pansoma_net_v2.train --tensors <set> [<set> ...] --kinds SNV --output <dir> [--test-chroms chr1] ...
    torchrun --nproc_per_node=N -m pansoma_net_v2.train --ddp ...           # several GPUs

A <set> is a directory holding merged SNV/ and INDEL/ (e.g. .../COLO829T_Illumina/tensors). Classes are
the labels 0 non, 1 somatic, 2 germline; tensors labelled -1 are not used. The test chromosomes (default chr1)
are left out; the others train, except the validation: whole ~1 Mb node blocks (--val-fraction of them), or
--val-chroms. The encoder's z-score statistics (and with --scalars the scalars') are fitted once on the
training tensors (or taken from --resume) and saved in every checkpoint; with --norm dataset (the default) BQ is
z-scored per dataset instead, each training dataset (set x kind) with statistics fitted on tensors of all its
chromosomes (labels not used), saved in the checkpoint; predict reuses them for a training dataset and fits them the
same way on any other (encode.py). MAPQ and path count keep the training statistics. Each epoch the somatic threshold of
the best validation F1 is found (metrics.py); the best checkpoint has the highest thresholded somatic F1 and
stores its threshold.

--site-features DIR feeds each tensor's population-database features (<DIR>/<SET>.<KIND>.npz, site_features.py) to
the head as extra scalars, through a fixed map to [0, 1]; the checkpoint records DIR and the map, and predict loads
the same files. With site features the validation also reports the post-PoN somatic threshold
(metrics.postpon_report: tensors the PoN rule tags are never called): val somatic_threshold_postpon, val.postpon and
the report-only keys somatic_at_recall_postpon, somatic_at_recall_postpon_callable and truth_threshold_postpon. The
checkpoint records the sha256 of every feature file read (feature_files), so predict can warn when one has changed.
The best epoch has the highest pre-PoN per-tensor somatic F1 (--best-by f1, the default) or, with --best-by postpon,
the highest post-PoN per-tensor somatic F1 (needs --site-features). The checkpoint records its call rule
(call_rule), the PASS threshold graph_vcf applies by default: with --best-by postpon "postpon-recall" at call_recall
(--call-recall, 0.8 by default): the threshold at which the post-PoN validation recall of the callable positives (those
the PoN does not tag) reaches it, val somatic_at_recall_postpon_callable; else "recall" (val somatic_at_recall at
graph_vcf's --target-recall, 0.9 by default). The post-PoN best-F1 threshold (somatic_threshold_postpon) stays
available as graph_vcf --threshold postpon.

--non-reasons R [R ...] trains and validates the tensors whose label reason is R (e.g. repeat_indel_without_truth,
-1 in the label files) as non where they lie inside the benchmark (the confident region or the somatic truth BED);
--non-reasons-in-bed R [R ...] does the same for reasons such as outside_confident_region where the tensor lies inside
the somatic truth BED. predict scores them the same way. The label files are not changed.

The data steps (indexes, label changes, features, chromosome split, validation selection) are `prepare_data`.
"""
import argparse
import json
import math
import os
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch
import torch.distributed as dist
import torch.nn as nn
from torch.nn.parallel import DistributedDataParallel
from torch.utils.data import DataLoader, Subset

from . import metrics
from . import site_features as sitef
from .env import triton_libcuda
from .chunks import ChunkLoader, build_cache
from .data import (CHECKPOINT_PATH_COUNT, CLASSES, KINDS, SCALARS, EpochSampler, TensorDataset, block_split,
                   common_path_count, load_parts, read_candidate_ids, somatic_truth, validation_truth)
from .encode import PLANES, compute_stats
from .model import PansomaNetV2, full_config, no_decay

AF = SCALARS.index("af")
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--tensors", nargs="+", required=True, help="merged tensor sets (directories with SNV/ INDEL/)")
    p.add_argument("--kinds", nargs="+", default=list(KINDS), choices=KINDS)
    p.add_argument("--output", required=True)
    p.add_argument("--cache-dir", help="index cache (default: <output>/index_cache)")
    p.add_argument("--train-chroms", nargs="+", help="default: every autosome not in --val-chroms/--test-chroms")
    p.add_argument("--test-chroms", nargs="*", default=["chr1"],
                   help="left out of training and validation, in every set (--test-chroms alone: none)")
    p.add_argument("--set-test-chroms", nargs="+", action="append", default=[], metavar=("SET", "CHROM"),
                   help="SET (one of --tensors) CHROM [CHROM ...]: also left out of training and validation, in that set "
                        "only (e.g. a test sample's chr1 while the other sets train on theirs); repeatable")
    p.add_argument("--val-chroms", nargs="+", default=[], help="validation chromosomes (default: node blocks)")
    p.add_argument("--val-fraction", type=float, default=0.05, help="node blocks of the training chromosomes")
    p.add_argument("--val-block-nodes", type=int, default=20000)
    p.add_argument("--scalars", action="store_true", help="feed the site scalars (data.SCALARS) to the head")
    p.add_argument("--site-features", metavar="DIR",
                   help="feed each tensor's site features (<DIR>/<SET>.<KIND>.npz: population AF of the exact allele and "
                        "of the position, PoN tag, homopolymer length; site_features.py) to the head as extra scalars; "
                        "a set without its file, or with a stale one, is an error; the checkpoint records DIR and the "
                        "transform, and their PoN tags also give the post-PoN validation thresholds")
    p.add_argument("--best-by", choices=["f1", "postpon"], default="f1",
                   help="best epoch (best.pth): highest validation somatic F1 (f1), or highest validation post-PoN "
                        "somatic F1, PoN-tagged tensors never called (postpon; needs --site-features); postpon also "
                        "makes the post-PoN recall threshold of --call-recall the checkpoint's call rule (graph_vcf's "
                        "default)")
    p.add_argument("--call-recall", type=float, default=0.8, choices=metrics.RECALLS,
                   help="with --best-by postpon: the call rule's post-PoN validation recall of the callable positives "
                        "(call_rule postpon-recall, call_recall)")
    p.add_argument("--epochs", type=int, default=30)
    p.add_argument("--epoch-samples", type=int, help="training tensors per epoch (a fresh random subset each epoch)")
    p.add_argument("--val-samples", type=int, help="fixed random subset of the validation tensors")
    p.add_argument("--batch-size", type=int, default=64, help="per GPU")
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--weight-decay", type=float, default=0.05)
    p.add_argument("--warmup-epochs", type=float, default=1.0)
    p.add_argument("--class-weights", default="balanced",
                   help="'balanced' (n / (3 n_c) over an epoch's tensors), 'sqrt' (its square root), 'pow:E' (balanced "
                        "to the power E: pow:0.5 = sqrt), 'none', or three numbers w_non,w_som,w_germ")
    p.add_argument("--non-fraction", type=float, default=1.0,
                   help="each epoch: every somatic / germline tensor and a fresh random fraction of the non ones")
    p.add_argument("--ignore-reasons", nargs="+", default=[],
                   help="label reasons left out of training and validation, as if -1 (e.g. "
                        "residual_partial_somatic_truth allele_partial_somatic_truth); the label files are unchanged")
    p.add_argument("--exclude-candidates",
                   help="file of candidate ids (first column; a header line candidate_id is skipped) taken out before "
                        "the model, as if never built: left out of training and validation (-1); predict leaves them "
                        "out of the test and graph_vcf out of the calls (predict's default is this file)")
    p.add_argument("--non-reasons", nargs="+", default=[],
                   help="label reasons trained and validated as class 0 (non) instead of their label where the tensor "
                        "lies inside the benchmark (confident region or somatic truth BED), e.g. "
                        "repeat_indel_without_truth (-1 in the label files); the label files are unchanged; predict "
                        "scores them as 0 too")
    p.add_argument("--non-reasons-in-bed", nargs="+", default=[],
                   help="label reasons trained and validated as class 0 where the tensor's GRCh38 position is inside the "
                        "set's somatic truth BED, e.g. outside_confident_region (somatic BED minus germline BED); "
                        "predict scores them as 0 too")
    p.add_argument("--keep-non-af", type=float,
                   help="non tensors with AF >= this are taken every epoch; the sampled rest weigh 1/--non-fraction "
                        "in the loss and the class weights use all tensors, so the expected loss is that of all non "
                        "(without the weights an AF-dependent sampling would bias p_somatic by AF)")
    p.add_argument("--depths", type=int, nargs=4, default=[3, 3, 27, 3])
    p.add_argument("--dims", type=int, nargs=4, default=[192, 384, 768, 1536])
    p.add_argument("--front", type=int, nargs="+", default=[64, 64], help="widths of the 1x1 layers before the stem")
    p.add_argument("--drop-path", type=float, default=0.1)
    p.add_argument("--drop-planes", nargs="+", default=[], metavar="PLANE",
                   help="encoder planes left out of the input (encode.PLANES), e.g. path_count")
    p.add_argument("--block", choices=["v2", "v1"], default="v2",
                   help="v2: ConvNeXt blocks as published; v1: mynet's (row-wise scalar GRN, GELU after each residual "
                        "sum); see model.py")
    p.add_argument("--decay-all", action="store_true",
                   help="weight decay on every parameter (as runs before 2026-09-29); default: not on biases, norm "
                        "weights and GRN gamma / beta")
    p.add_argument("--stats-samples", type=int, default=20000,
                   help="tensors used to fit the z-score: of the training data, and with --norm dataset of each dataset")
    p.add_argument("--norm", choices=["dataset", "train"], default="dataset",
                   help="dataset: BQ z-scored with the statistics of each sample's own dataset (a set's SNV/ or "
                        "INDEL/), fitted on tensors of all its chromosomes here (labels not used) and by predict on "
                        "each other predicted dataset; train: with the training statistics for every dataset (every "
                        "checkpoint before 2026-10-13). MAPQ and path count always use the training statistics "
                        "(norm dataset checkpoints of 2026-10-13 to 10-19 z-scored MAPQ per dataset too)")
    p.add_argument("--amp", choices=["bf16", "off"], default="bf16")
    p.add_argument("--no-compile", dest="compile", action="store_false",
                   help="skip torch.compile (default on CUDA: 3x faster on an H100, ~2 min to compile)")
    p.add_argument("--no-channels-last", dest="channels_last", action="store_false")
    p.add_argument("--num-workers", type=int, default=8,
                   help="DataLoader workers for the statistics (random reads from BeeGFS are latency-bound: more "
                        "workers read faster); 0 also reads the validation without workers")
    p.add_argument("--val-workers", type=int, default=8)
    p.add_argument("--micro-batches", type=int, default=1,
                   help="forward/backward each batch in this many parts (gradient accumulation: the same update as the "
                        "whole batch, loss normalized over the whole batch, with the GPU memory of one part)")
    p.add_argument("--prefetch-factor", type=int, default=2, help="batches in flight per worker")
    p.add_argument("--read-threads", type=int, default=16,
                   help="training: threads reading chunks of consecutive shard rows (chunks.ChunkLoader)")
    p.add_argument("--chunk-rows", type=int, default=512, help="rows per sequential read (512 = 83 MB)")
    p.add_argument("--window-chunks", type=int, default=32, help="chunks mixed into one shuffled pool of batches")
    p.add_argument("--chunk-cache", help="zstd chunk cache of the training shards (default: <cache-dir>/chunks; "
                   "built where missing, see chunks.py)")
    p.add_argument("--no-chunk-cache", dest="use_chunk_cache", action="store_false",
                   help="read the chunks from the raw shards")
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


def set_test_chroms(args):
    """{resolved --tensors path: sorted chromosomes} of --set-test-chroms."""
    sets = {Path(s).resolve() for s in args.tensors}
    held = {}
    for entry in args.set_test_chroms:
        if len(entry) < 2:
            raise SystemExit(f"--set-test-chroms {' '.join(entry)}: give SET CHROM [CHROM ...]")
        path, chroms = Path(entry[0]).resolve(), entry[1:]
        if path not in sets:
            raise SystemExit(f"--set-test-chroms: {entry[0]} is not one of --tensors")
        if set(chroms) - set(AUTOSOMES):
            raise SystemExit(f"--set-test-chroms: {sorted(set(chroms) - set(AUTOSOMES))} are not autosomes")
        held[path] = sorted(set(held.get(path, [])) | set(chroms))
    return held


def class_weights(spec, counts):
    """counts: tensors per class in one epoch."""
    if spec == "none":
        return torch.ones(len(CLASSES))
    if spec in ("balanced", "sqrt") or spec.startswith("pow:"):
        try:
            power = {"balanced": 1.0, "sqrt": 0.5}.get(spec) or float(spec[4:])
        except ValueError:
            raise SystemExit(f"--class-weights {spec}: pow:E needs a number E")
        counts = np.maximum(counts, 1)
        w = counts.sum() / (len(CLASSES) * counts)
        return torch.tensor(w ** power, dtype=torch.float32)
    w = [float(x) for x in spec.split(",")]
    if len(w) != len(CLASSES):
        raise SystemExit("--class-weights needs three numbers: non,somatic,germline")
    return torch.tensor(w)


def epoch_class_weights(args, counts):
    """(class weights of the loss, importance) for the training class counts: the weights follow one epoch's classes;
    with --keep-non-af (importance) the importance weights restore all non."""
    importance = args.keep_non_af is not None and args.non_fraction < 1.0
    epoch_counts = counts.astype(np.float64) * np.array([1.0 if importance else args.non_fraction, 1.0, 1.0])
    return class_weights(args.class_weights, epoch_counts), importance


def predict_probs(model, loader, device, amp):
    """(labels, probabilities) of loader's samples on this rank."""
    model.eval()
    labels, probs = [], []
    with torch.no_grad():
        for x, blocks, scalars, y, part in loader:
            with torch.autocast(device.type, dtype=torch.bfloat16, enabled=amp):
                logits = model(x.to(device, non_blocking=True), blocks.to(device, non_blocking=True),
                               scalars.to(device, non_blocking=True), part.to(device, non_blocking=True))
            probs.append(torch.softmax(logits.float(), 1).cpu().numpy())
            labels.append(y.numpy())
    return (np.concatenate(labels) if labels else np.zeros(0, np.int64),
            np.concatenate(probs) if probs else np.zeros((0, len(CLASSES)), np.float32))


def evaluate(model, loader, criterion, run, amp, truth=None, pon=None):
    """Validation report (metrics.report at the best tensor-level somatic threshold), weighted loss and, with
    truth = (truth keys, other-kind truth keys, per-sample matched keys), metrics.truth_report at the same
    threshold; with pon (each validation sample's PoN tag, site_features) the post-PoN report
    (metrics.postpon_report) as postpon, somatic_threshold_postpon, somatic_at_recall_postpon,
    somatic_at_recall_postpon_callable and, with truth, truth_threshold_postpon; all ranks."""
    return validation_report(*validation_predictions(model, loader, run, amp), criterion, truth, pon)


def validation_predictions(model, loader, run, amp):
    """(dataset indices, labels, probabilities) of the loader's samples over all ranks, in dataset-index order."""
    labels, probs = predict_probs(model, loader, run.device, amp)
    order = np.asarray(list(loader.sampler), np.int64)[:len(labels)]  # dataset indices of this rank's samples
    if run.ddp:
        parts = [None] * run.world
        dist.all_gather_object(parts, (order, labels, probs))
        order, labels, probs = (np.concatenate([p[k] for p in parts]) for k in range(3))
    s = np.argsort(order, kind="stable")
    return order[s], labels[s], probs[s]


def validation_report(order, labels, probs, criterion, truth=None, pon=None):
    """evaluate's report from validation_predictions (order: the samples' dataset indices, which index truth's
    matches and pon)."""
    threshold, _, _, _ = metrics.best_somatic_threshold(labels, probs)
    out = metrics.report(labels, probs, threshold)
    if truth is not None:
        keys, other, matches = truth
        out["truth"] = metrics.truth_report(keys, [matches[i] for i in order], labels, probs, threshold, other)
    if pon is not None:
        tr = None
        if truth is not None:
            keys, other, matches = truth
            tr = (keys, [matches[i] for i in order], other)
        post = metrics.postpon_report(labels, probs, pon[order], tr)
        out["postpon"] = post
        out["somatic_threshold_postpon"] = post["threshold"]
        out["somatic_at_recall_postpon"] = post["at_recall"]
        out["somatic_at_recall_postpon_callable"] = post["at_recall_callable"]
        if tr is not None:
            out["truth_threshold_postpon"] = post["truth"]["best"]["threshold"]
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
                 f"truth, {at['other_tp']:,} of the other kind found, ceiling {tr['ceiling']:.3f}, best F1 "
                 f"{tr['best']['f1']:.3f} @ p>={tr['best']['threshold']:.3f})")
    if report.get("postpon"):
        q = report["postpon"]
        line += (f" | post-PoN F1 {q['f1']:.3f} (P {q['precision']:.3f} R {q['recall']:.3f} @ p>={q['threshold']:.3f}; "
                 f"{q['tagged']:,} tagged, ceiling {q['ceiling']:.3f})")
        if q.get("truth"):
            b = q["truth"]["best"]
            line += f" truth F1 {b['f1']:.3f} @ p>={b['threshold']:.3f}"
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


def dataset_stats(index, positions, samples, seed, num_workers, log, what):
    """The z-score statistics of one dataset (compute_stats over up to `samples` random tensors of index's positions,
    labels not used): a row of the encoder's per-dataset table under --norm dataset."""
    if not len(positions):
        raise SystemExit(f"{what}: no tensors to fit the statistics on")
    dataset = TensorDataset([(index, np.asarray(positions))])
    order = EpochSampler(len(dataset), samples, shuffle=True, seed=seed)
    return retry_workers(lambda: compute_stats(make_loader(dataset, order, num_workers, 256, persistent=False), samples),
                         log, what)


def started_loader(dataset, sampler, num_workers, batch_size, log, what, prefetch_factor=2):
    """A persistent DataLoader whose workers are running (one batch fetched), retried as in retry_workers."""
    def start():
        loader = make_loader(dataset, sampler, num_workers, batch_size, persistent=True,
                             prefetch_factor=prefetch_factor)
        if num_workers > 0 and len(loader):
            next(iter(loader))
        return loader
    return retry_workers(start, log, what)


def exclude(index, candidate_ids):
    """Set the training (label) and validation (eval_label) labels of the index's tensors whose candidate id is in
    candidate_ids to -1, in place (--exclude-candidates); returns how many."""
    a = index.arrays
    pos = np.array([k for k, c in enumerate(index.candidates()) if c in candidate_ids], dtype=np.int64)
    label, evl = a["label"].copy(), a["eval_label"].copy()
    label[pos], evl[pos] = -1, -1
    a["label"], a["eval_label"] = label, evl
    return len(pos)


def as_non(index, reasons):
    """Set the training and evaluation label of the index's tensors whose label reason is in `reasons` to 0 where they
    lie inside the benchmark (the confident region or the somatic truth BED, where test-time calls are scored), in
    place (--non-reasons); outside it they stay -1. Returns (tensors with the reasons, of them set to 0)."""
    a = index.arrays
    mask = np.isin(a["reason"], [k for k, r in enumerate(index.meta["reasons"]) if r in set(reasons)])
    inside = mask & (a["in_region"] | in_somatic_bed(index))
    label, evl = a["label"].copy(), a["eval_label"].copy()
    label[inside] = 0
    evl[inside] = 0
    a["label"], a["eval_label"] = label, evl
    return int(mask.sum()), int(inside.sum())


def in_somatic_bed(index):
    """Tensors whose pos0 (GRCh38 start; off-reference: the anchor's middle base, as the labeller's region test) lies
    inside the somatic truth BED of the set's labels (labels.manifest.json); the somatic truth is complete there,
    whatever the germline BED says. A point test: agrees with the labeller's span test on 99.3 % of HCC1395T Illumina
    INDEL tensors (the rest cross a BED edge)."""
    from .bed import Bed
    manifest = Path(index.dir) / "labels.manifest.json"
    truth = (json.loads(manifest.read_text()).get("truth") or {}) if manifest.exists() else {}
    path = (truth.get("somatic") or {}).get("bed")
    a, names = index.arrays, index.meta["chroms"]
    inside = np.zeros(len(a["label"]), bool)
    if not path:
        return inside
    bed = Bed.read(path)
    for c, name in enumerate(names):
        if name not in bed.data:
            continue
        sel = np.flatnonzero((a["chrom"] == c) & (a["pos0"] >= 0))
        starts, ends = bed.data[name]
        k = np.searchsorted(starts, a["pos0"][sel], side="right") - 1
        ok = k >= 0
        ok[ok] = ends[k[ok]] > a["pos0"][sel][ok]
        inside[sel[ok]] = True
    return inside


def as_non_in_bed(index, reasons):
    """--non-reasons-in-bed: tensors with a reason in `reasons` (e.g. outside_confident_region) whose pos0 is inside the
    somatic truth BED become 0 for training and evaluation, in place; returns the count."""
    a = index.arrays
    mask = np.isin(a["reason"], [k for k, r in enumerate(index.meta["reasons"]) if r in set(reasons)]) & \
        in_somatic_bed(index)
    label, evl = a["label"].copy(), a["eval_label"].copy()
    label[mask] = 0
    evl[mask] = 0
    a["label"], a["eval_label"] = label, evl
    return int(mask.sum())


def attach_features(index, directory, model_input):
    """site_features.attach as a clean exit: a missing or stale feature file stops the run."""
    try:
        return sitef.attach(index, directory, model_input)
    except ValueError as e:
        raise SystemExit(f"site features: {e}") from None


def prepare_data(args, log, refresh_index=True, barrier=None):
    """The data of a training run from its arguments: the index of every set x kind with the label changes of
    --exclude-candidates / --non-reasons / --non-reasons-in-bed, the site features (--site-features, which also give
    the PoN tags), the chromosome split, the training and validation selections (node-block split or --val-chroms,
    --set-test-chroms, --ignore-reasons), the datasets, the --val-samples subset, the truth-level validation keys
    and the validation samples' PoN tags. log: the log function; refresh_index: build / refresh the index cache first
    (with --rebuild-index: rebuild it), train's rank 0; barrier: called after that (DDP).
    Returns a SimpleNamespace: cache, everything ([(KindIndex, positions)] of every set x kind), reads_path_count,
    path_count, haplotypes, site_record, feature_files ({"site_features": {<SET>.<KIND>: sha256}} or None),
    train_chroms, val_chroms, test_chroms, held (--set-test-chroms), train_parts, val_parts, train_set, val_set (a
    Subset with --val-samples), val_keep, val_truth, val_pon, counts (training class counts)."""
    cache = args.cache_dir or str(Path(args.output) / "index_cache")
    if refresh_index:  # rank 0 builds / refreshes the index cache, the others read it
        load_parts(args.tensors, args.kinds, cache, rebuild=args.rebuild_index)
    if barrier is not None:
        barrier()
    everything = load_parts(args.tensors, args.kinds, cache)
    reads_path_count = "path_count" not in args.drop_planes  # the only plane made from channel 6
    if reads_path_count:
        try:
            path_count, haplotypes = common_path_count(everything)
        except ValueError as e:
            raise SystemExit(str(e)) from None
        log(f"channel 6 (path count): {path_count}" + (f", --haplotypes {haplotypes}" if haplotypes else ""))
    else:  # the model never reads channel 6: sets of any storage
        path_count, haplotypes = None, sorted({i.haplotypes for i, _ in everything if i.haplotypes is not None})
        log("channel 6 (path count) not read (--drop-planes path_count)")
    if args.exclude_candidates:
        excluded = read_candidate_ids(args.exclude_candidates)
        log(f"--exclude-candidates {args.exclude_candidates}: {len(excluded):,} candidate ids")
        for index, _ in everything:
            log(f"{index.dir}: {exclude(index, excluded):,} tensors excluded")
    if args.non_reasons:
        for index, _ in everything:
            n, inside = as_non(index, args.non_reasons)
            log(f"{index.dir}: --non-reasons {' '.join(args.non_reasons)}: {n:,} tensors, {inside:,} inside the "
                f"benchmark trained and validated as non")
    if args.non_reasons_in_bed:
        for index, _ in everything:
            n = as_non_in_bed(index, args.non_reasons_in_bed)
            log(f"{index.dir}: --non-reasons-in-bed {' '.join(args.non_reasons_in_bed)}: {n:,} tensors inside "
                f"the somatic BED trained and validated as non")
    files, loaded = {}, {}  # <SET>.<KIND> -> sha256 of the file read / its arrays' record
    for index, _ in everything:  # before any TensorDataset: they read index.arrays["site_features"]
        if args.site_features:
            data = attach_features(index, args.site_features, True)
            files[sitef.key(index)], loaded[sitef.key(index)] = data["sha256"], data
            log(sitef.describe(index, data) + " (model input)")
            for w in sitef.pon_check(data):
                log(f"warning: {w}")
    site_record = sitef.record(args.site_features, loaded) if args.site_features else None
    feature_files = {"site_features": files} if files else None
    available = sorted({c for index, _ in everything for c in index.meta["chroms"]})
    train_chroms, val_chroms, test_chroms, missing = split_chroms(args, available)
    log(f"train {train_chroms}\nval {val_chroms}\ntest (left out) {test_chroms}"
        + (f"\nnot in the data: {missing}" if missing else ""))
    held = set_test_chroms(args)
    for path, chroms in held.items():
        log(f"test (left out) of {path} only: {chroms}")

    def chroms_of(index, chroms):
        """The chromosomes of `chroms` that index's set trains / validates on (without its --set-test-chroms)."""
        return [c for c in chroms if c not in held.get(index.dir.parent, ())]

    # training: labels >= 0; validation (like the test): evaluation labels (off-reference tensors without a
    # truth match count as non), since the threshold chosen on it is applied to the test chromosome
    if val_chroms:
        train_parts = [(i, i.select(chroms_of(i, train_chroms))) for i, _ in everything]
        val_parts = [(i, i.select(chroms_of(i, val_chroms), evaluation=True)) for i, _ in everything]
    else:
        train_parts, val_parts = [], []
        for i, _ in everything:
            tr, va = block_split(i, i.select(chroms_of(i, train_chroms), labelled=False), args.val_fraction,
                                 args.val_block_nodes, args.seed)
            train_parts.append((i, tr[i.arrays["label"][tr] >= 0]))
            val_parts.append((i, va[i.arrays["eval_label"][va] >= 0]))
        log(f"validation: {args.val_fraction:.0%} of the {args.val_block_nodes}-node blocks of the training chromosomes")
    if args.ignore_reasons:
        def without_ignored(i, pos):
            codes = [k for k, r in enumerate(i.meta["reasons"]) if r in set(args.ignore_reasons)]
            return pos[~np.isin(i.arrays["reason"][pos], codes)]
        before = sum(len(p) for _, p in train_parts), sum(len(p) for _, p in val_parts)
        train_parts = [(i, without_ignored(i, p)) for i, p in train_parts]
        val_parts = [(i, without_ignored(i, p)) for i, p in val_parts]
        log(f"left out (--ignore-reasons {' '.join(args.ignore_reasons)}): "
            f"{before[0] - sum(len(p) for _, p in train_parts):,} training, "
            f"{before[1] - sum(len(p) for _, p in val_parts):,} validation tensors")
    train_set, val_set = TensorDataset(train_parts), TensorDataset(val_parts, labels="eval_label")
    if val_parts and val_set.shape != train_set.shape:
        raise SystemExit(f"validation tensors {val_set.shape} differ from the training tensors {train_set.shape}")
    log(f"tensors of shape {train_set.shape}")
    n_off = sum(int((i.arrays["eval_label"][v] != i.arrays["label"][v]).sum()) for i, v in val_parts)
    log(f"validation includes {n_off:,} off-reference tensors without a truth match, as non")
    counts = train_set.class_counts()
    for index, pos in train_parts:
        labels = index.arrays["label"][pos]
        lm = index.meta["labels"]
        log(f"{index.dir}: labels {lm.get('version') or (lm.get('format'), (lm.get('rules_sha256') or '')[:12])} "
            f"({lm['created']}), "
            f"train {len(pos):,} " + " ".join(f"{c}={int((labels == k).sum()):,}" for k, c in enumerate(CLASSES)))
    log(f"training tensors {len(train_set):,} ({dict(zip(CLASSES, counts.tolist()))}), validation {len(val_set):,} "
        f"({dict(zip(CLASSES, val_set.class_counts().tolist()))})")
    if len(train_set) == 0 or len(val_set) == 0:
        raise SystemExit("empty training or validation set")
    val_keep = None
    if args.val_samples and args.val_samples < len(val_set):
        val_keep = np.sort(np.random.default_rng(args.seed).choice(len(val_set), args.val_samples, replace=False))
    # truth-level validation: the truth alleles of the validation region (with or without a tensor) and, per
    # validation sample, the truth alleles its tensor stands for; keys are (part, truth id)
    # (other_keys: the truth alleles of the other kind, true calls when a validation tensor matches one)
    val_truth, truth_keys, other_keys = None, set(), set()
    for p, ((index, _), (_, va)) in enumerate(zip(train_parts, val_parts)):
        t_all = somatic_truth(index, chroms=chroms_of(index, val_chroms or train_chroms))
        if t_all is None:
            log(f"{index.dir}: no somatic.recall.tsv, no truth-level validation")
            truth_keys = None
            break
        vt = set(t_all) if val_chroms else validation_truth(
            index, index.select(chroms_of(index, train_chroms), labelled=False), t_all, args.val_fraction,
            args.val_block_nodes, args.seed)
        truth_keys |= {(p, t) for t in vt}
        other_keys |= {(p, t) for t in somatic_truth(index, chroms_of(index, val_chroms or train_chroms),
                                                     other_kind=True)}
    if truth_keys is not None:
        ids = val_keep.tolist() if val_keep is not None else range(len(val_set))
        matches = [{(int(val_set.part[i]), int(t)) for t in val_parts[val_set.part[i]][0].truth_of(val_set.position[i])}
                   for i in ids]
        val_truth = (truth_keys, other_keys, matches)
        log(f"truth-level validation: {len(truth_keys):,} somatic truth alleles in the validation region "
            f"({len(set().union(*matches) & truth_keys) if matches else 0:,} with a validation tensor)")
    val_pon = None
    if args.site_features:  # each validation sample's PoN tag, in val_set order (then the --val-samples subset)
        val_pon = np.concatenate([i.arrays["pon_tag"][va] for i, va in val_parts]).astype(np.int8)
        if val_keep is not None:
            val_pon = val_pon[val_keep]
        vl = val_set.label if val_keep is None else val_set.label[val_keep]
        log(f"post-PoN validation: {int((val_pon == 1).sum()):,} of {len(val_pon):,} validation tensors PoN-tagged "
            f"(" + ", ".join(f"{c} {int(((val_pon == 1) & (vl == k)).sum()):,}" for k, c in enumerate(CLASSES))
            + f"), {int((val_pon == -1).sum()):,} without a GRCh38 allele (callable)")
    if val_keep is not None:
        val_set = Subset(val_set, val_keep.tolist())

    return SimpleNamespace(
        cache=cache, everything=everything, reads_path_count=reads_path_count, path_count=path_count,
        haplotypes=haplotypes, site_record=site_record, feature_files=feature_files,
        train_chroms=train_chroms, val_chroms=val_chroms, test_chroms=test_chroms, held=held, train_parts=train_parts,
        val_parts=val_parts, train_set=train_set, val_set=val_set, val_keep=val_keep, val_truth=val_truth,
        val_pon=val_pon, counts=counts)


def main(argv=None):
    args = parse_args(argv)
    if args.best_by == "postpon" and not args.site_features:
        raise SystemExit("--best-by postpon needs the PoN tags of --site-features")
    run = Run(args)
    if run.main:
        (run.out / "args.json").write_text(json.dumps(vars(args), indent=2) + "\n")
    torch.manual_seed(args.seed)
    prep = prepare_data(args, run.log, refresh_index=run.main, barrier=run.barrier)
    train_parts, train_set, val_set, counts = prep.train_parts, prep.train_set, prep.val_set, prep.counts
    path_count, haplotypes, reads_path_count = prep.path_count, prep.haplotypes, prep.reads_path_count
    cache, held = prep.cache, prep.held

    model = PansomaNetV2(len(CLASSES), args.depths, args.dims, args.front, args.drop_path,
                         scalars=len(SCALARS) if args.scalars else 0, block=args.block, drop_planes=args.drop_planes,
                         norm=args.norm, site_features=len(sitef.SITE_FEATURES) if args.site_features else 0)
    checkpoint = torch.load(args.resume, map_location="cpu", weights_only=False) if args.resume else None
    if checkpoint is not None:
        if full_config(checkpoint["config"]) != model.config:
            raise SystemExit(f"--resume model config {checkpoint['config']} differs from the arguments' {model.config}")
        recorded = checkpoint.get("site_features")
        if recorded and recorded.get("version") != sitef.VERSION:
            raise SystemExit(f"--resume checkpoint's site features are {recorded.get('version')}, this code's "
                             f"{sitef.VERSION}")
        trained = checkpoint.get("path_count_storage", CHECKPOINT_PATH_COUNT)
        if reads_path_count and trained != path_count:
            raise SystemExit(f"--resume checkpoint was trained on channel 6 {trained}, the tensors are {path_count}")
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
    part_stats = []
    if args.norm == "dataset":  # one row per training dataset (TensorDataset part order), kept on --resume
        saved = {d["directory"]: d["stats"] for d in (checkpoint or {}).get("dataset_stats", [])}
        box = [None]
        if run.main:
            box[0] = []
            for i, _ in train_parts:
                if str(i.dir) in saved:
                    box[0].append(saved[str(i.dir)])
                    run.log(f"dataset statistics of {i.dir} kept from --resume")
                    continue
                t0 = time.time()
                # every chromosome of the set, the test one too (labels not used): one chromosome is not a typical
                # sample of a set (chr1 holds most pericentromeric tensors), and predict fits other sets the same way
                st = dataset_stats(i, i.select(None, labelled=False), args.stats_samples, args.seed + 12345,
                                   args.num_workers, run.log, f"dataset statistics of {i.dir}")
                box[0].append(st)
                run.log(f"dataset statistics of {i.dir} from {st['tensors']:,} tensors of all its chromosomes "
                        f"({time.time() - t0:.0f} s; z-scored per dataset: {', '.join(model.encoder.dataset_norm)}): "
                        + json.dumps({k: st[k] for k in model.encoder.dataset_norm}))
            (run.out / "dataset_stats.json").write_text(json.dumps(
                [dict(directory=str(i.dir), stats=st) for (i, _), st in zip(train_parts, box[0])], indent=2) + "\n")
        if run.ddp:
            dist.broadcast_object_list(box, src=0)
        part_stats = box[0]
    model.encoder.set_part_stats(part_stats)
    if args.scalars and checkpoint is None:  # all training tensors (the arrays are in memory)
        s = train_set.scalars[:, :len(SCALARS)].astype(np.float64)  # without the site features
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
    weights, importance = epoch_class_weights(args, counts)
    weights = weights.to(run.device)
    criterion = nn.CrossEntropyLoss(weight=weights, ignore_index=-1)
    per_sample = nn.CrossEntropyLoss(weight=weights, ignore_index=-1, reduction="none")
    keep_af = train_set.scalars[:, AF] >= args.keep_non_af if args.keep_non_af is not None else None
    # rank 0 of 1: ChunkLoader splits the epoch over the ranks itself
    train_sampler = EpochSampler(len(train_set), args.epoch_samples, True, 0, 1, args.seed,
                                 labels=train_set.label, non_fraction=args.non_fraction, keep=keep_af)
    if importance:
        n_keep = int((keep_af & (train_set.label == 0)).sum())
        run.log(f"non tensors with AF >= {args.keep_non_af}: {n_keep:,} taken every epoch; the others sampled at "
                f"{args.non_fraction} with importance weight {1 / args.non_fraction:.2f}")
    run.log(f"class weights {dict(zip(CLASSES, [round(w, 4) for w in weights.tolist()]))}; each epoch "
            f"{train_sampler.epoch_size():,} tensors (non fraction {args.non_fraction})")
    val_sampler = EpochSampler(len(val_set), None, False, run.rank, run.world)
    chunk_cache = None
    if args.use_chunk_cache:
        files = [train_set.files[f] for f in np.unique(train_set.shard).tolist()]
        where = args.chunk_cache or str(Path(cache) / "chunks")
        if run.main:  # rank 0 builds what is missing, the others wait
            build_cache(files, where, args.chunk_rows, args.read_threads, run.log)
        run.barrier()
        chunk_cache = build_cache(files, where, args.chunk_rows, args.read_threads, run.log)
    train_loader = ChunkLoader(train_set, train_sampler, args.batch_size, args.read_threads, args.chunk_rows,
                               args.window_chunks, run.rank, run.world, args.seed, cache=chunk_cache)
    run.log(f"training reads chunks of {args.chunk_rows} rows with {args.read_threads} threads, "
            f"{args.window_chunks} chunks per shuffled pool, from "
            + (f"the zstd chunk cache {where}" if chunk_cache else "the raw shards"))
    val_loader = started_loader(val_set, val_sampler, args.val_workers if args.num_workers else 0, args.batch_size,
                                run.log, "validation", args.prefetch_factor)
    steps_per_epoch = max(1, len(train_loader))
    total_steps, warmup = args.epochs * steps_per_epoch, int(args.warmup_epochs * steps_per_epoch)
    if args.decay_all:
        groups = [dict(params=list(model.parameters()))]
    else:
        named = list(model.named_parameters())
        groups = [dict(params=[p for n, p in named if not no_decay(n, p)]),
                  dict(params=[p for n, p in named if no_decay(n, p)], weight_decay=0.0)]
        run.log(f"weight decay {args.weight_decay} on {sum(p.numel() for p in groups[0]['params']) / 1e6:.1f} M parameters, "
                f"none on {sum(p.numel() for p in groups[1]['params']):,} (biases, norms, GRN)")
    optimizer = torch.optim.AdamW(groups, lr=args.lr, weight_decay=args.weight_decay, fused=run.device.type == "cuda")
    schedule = lambda s: (s + 1) / max(1, warmup) if s < warmup else \
        0.5 * (1 + math.cos(math.pi * (s - warmup) / max(1, total_steps - warmup)))  # noqa: E731
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, schedule)
    start_epoch, best = 0, {"score": -1.0, "loss": math.inf, "epoch": None, "threshold": None}
    if checkpoint is not None:
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
        start_epoch, best = checkpoint["epoch"], checkpoint.get("best", best)
        if best.get("select", "f1") != args.best_by:  # scores of another selection are not comparable
            run.log(f"--resume: the checkpoint's best epoch {best.get('epoch')} was selected by "
                    f"{best.get('select', 'f1')}, this run selects by {args.best_by}: best epoch chosen afresh")
            best = {"score": -1.0, "loss": math.inf, "epoch": None, "threshold": None}
    amp = args.amp == "bf16" and run.device.type == "cuda"
    n_params = sum(p.numel() for p in model.parameters())
    run.log(f"model {model.config}, {n_params / 1e6:.1f} M parameters, input planes {len(PLANES)}; "
            f"device {run.device} x {run.world}, amp {amp}, compile {args.compile and run.device.type == 'cuda'}, "
            f"channels_last {args.channels_last and run.device.type == 'cuda'}; {steps_per_epoch} steps per epoch")

    # graph_vcf's default PASS threshold
    call_rule, call_recall = ("postpon-recall", args.call_recall) if args.best_by == "postpon" else ("recall", None)

    def save(path, epoch, val, threshold):
        if run.main:
            payload = dict(format="pansoma_net_v2", config=model.config, model_state_dict=model.state_dict(),
                           optimizer_state_dict=optimizer.state_dict(), scheduler_state_dict=scheduler.state_dict(),
                           epoch=epoch, best=best, val=val, somatic_threshold=threshold, classes=list(CLASSES),
                           planes=list(PLANES), scalars=list(SCALARS) if args.scalars else [],
                           site_features=prep.site_record,
                           feature_files=prep.feature_files, call_rule=call_rule, call_recall=call_recall,
                           stats=model.encoder.stats(), args=vars(args), path_count_storage=path_count,
                           dataset_stats=[dict(directory=str(i.dir), stats=st) for (i, _), st in zip(train_parts, part_stats)],
                           path_count_haplotypes=haplotypes, tensor_shape=list(train_set.shape),
                           data=[dict(directory=str(i.dir), labels=i.meta["labels"], train=int(len(p)))
                                 for i, p in train_parts],
                           chroms=dict(train=prep.train_chroms, val=prep.val_chroms, test=prep.test_chroms,
                                       **({"set_test": {str(k): v for k, v in held.items()}} if held else {})))
            tmp = run.out / (path + ".tmp")
            torch.save(payload, tmp)
            tmp.replace(run.out / path)

    for epoch in range(start_epoch, args.epochs):
        net.train()
        train_sampler.set_epoch(epoch)
        if run.device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(run.device)
        t0, loss_sum, seen, correct = time.time(), 0.0, 0, 0
        for step, (x, blocks, scalars, y, part) in enumerate(train_loader):
            x, blocks, scalars, y, part = (t.to(run.device, non_blocking=True) for t in (x, blocks, scalars, y, part))
            optimizer.zero_grad(set_to_none=True)
            if args.micro_batches == 1:
                with torch.autocast(run.device.type, dtype=torch.bfloat16, enabled=amp):
                    logits = net(x, blocks, scalars, part)
                if importance:  # a sampled low-AF non tensor stands for 1 / non_fraction of them
                    iw = torch.where((y == 0) & (scalars[:, AF] < args.keep_non_af), 1.0 / args.non_fraction, 1.0)
                    loss = (per_sample(logits.float(), y) * iw).sum() / (weights[y] * iw).sum()
                else:
                    loss = criterion(logits.float(), y)
                loss.backward()
                correct += int((logits.argmax(1) == y).sum())
            else:  # the batch's loss, normalized over the whole batch, from parts: the same gradient, less memory
                iw = (torch.where((y == 0) & (scalars[:, AF] < args.keep_non_af), 1.0 / args.non_fraction, 1.0)
                      if importance else torch.ones_like(y, dtype=torch.float32))
                denom = (weights[y] * iw).sum()
                loss = torch.zeros((), device=run.device)
                for piece in torch.arange(len(y), device=run.device).chunk(args.micro_batches):
                    with torch.autocast(run.device.type, dtype=torch.bfloat16, enabled=amp):
                        logits = net(x[piece], blocks[piece], scalars[piece], part[piece])
                    piece_loss = (per_sample(logits.float(), y[piece]) * iw[piece]).sum() / denom
                    piece_loss.backward()
                    loss += piece_loss.detach()
                    correct += int((logits.argmax(1) == y[piece]).sum())
            optimizer.step()
            scheduler.step()
            loss_sum += loss.item() * y.numel()
            seen += y.numel()
            if run.main and (step + 1) % 200 == 0:
                print(f"epoch {epoch + 1} step {step + 1}/{steps_per_epoch} loss {loss_sum / seen:.4f} "
                      f"lr {scheduler.get_last_lr()[0]:.2e} {seen * run.world / (time.time() - t0):.0f} tensors/s",
                      flush=True)
        train_time = time.time() - t0
        val = evaluate(net, val_loader, criterion, run, amp, prep.val_truth, prep.val_pon)
        threshold = val["threshold"]
        score = val["postpon"]["f1"] if args.best_by == "postpon" else val["thresholded"]["somatic"]["f1"]
        improved = score > best["score"] or (score == best["score"] and val["loss"] < best["loss"])
        if improved:
            best = {"score": score, "select": args.best_by, "loss": val["loss"], "epoch": epoch + 1,
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
    run.log(f"done; best epoch {best['epoch']} ({best.get('select', 'f1')} {best['score']:.4f}, "
            f"somatic threshold {best['threshold']}) -> {run.out / 'best.pth'}")
    if run.ddp:
        dist.destroy_process_group()


if __name__ == "__main__":
    main()
