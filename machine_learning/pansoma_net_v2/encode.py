"""Input encoding: int8 tensors (B, 8, 200, 101) -> float planes (B, 36, 200, 104), on the GPU.

The stored channels (tensor format `indexed-gam-candidate`): 0 read base, 1 base quality, 2 site allele,
3 MAPQ, 4 alignment operation, 5 graph base, 6 path count, 7 strand. A cell without evidence is 0 in every
channel, so `covered = ch4 > 0`.

Planes, in order (`PLANES`):
  one-hot, not normalized    read base, site allele, graph base, operation (6 each: codes 1-6), strand (2)
  masked z-score             base quality, MAPQ, path count: (value - mean) / std on valid covered cells,
                             0 elsewhere (valid: BQ >= 0, MAPQ >= 0, path count > 0)
  flags                      covered (ch4 > 0), bq_missing (covered and BQ = -1: a gap or no quality)
  derived                    differs (covered and read base != graph base); the row's site block, like
                             DeepVariant's read_supports_variant filled over the whole read: A1 (the
                             representative allele), ALT (A2..Ak), REF, OTHER, each 1 on the covered cells
                             of the rows in that block (from the summary's row_groups, see data.py)

`drop` leaves planes out of the output (e.g. path_count: the path count says how common the node is in HPRC,
which separates GRCh38 nodes shared by all haplotypes from alternative-branch nodes). Every plane is 0 where there
is no evidence, so padding stays 0 after encoding. The width is zero-padded from 101 to 104 so the 4 x 4, stride 4
stem covers the last column.

Z-score statistics (`compute_stats`: mean and std over the valid covered cells of a sample of tensors):
  - the encoder's buffers `mean` / `std`, fitted on the training tensors and saved in the checkpoint, are used for
    path count and MAPQ always and, without a dataset table, for BQ too (checkpoints with config norm "train", the
    only kind before 2026-10-13). Path count is a property of the graph (stored relative to its haplotypes); MAPQ
    comes from the pipeline's own vg alignment, whose scale is the same for every dataset, while its mean depends on
    which regions the tensors come from (chr1 alone: 0.5 SD below chr2-22 on ONT, mostly the pericentromeric tensors);
  - with norm "dataset" (`set_part_stats`), the encoder's `dataset_norm` channels (new models: `DATASET_NORM`, BQ) are
    z-scored with the statistics of the sample's own dataset (a set's SNV/ or INDEL/ directory, fitted on tensors of
    all its chromosomes): `forward(..., part)` picks each sample's row of the table. Training fits one row per
    training dataset; predict takes the checkpoint's row for a training dataset and fits one on any other. BQ scales
    depend on the instrument, chemistry and basecaller (binned or full-range, e.g. HG008T vs COLO829T ONT 41.2 vs
    29.8), so a model then sees each dataset's BQ relative to that dataset, not to its training data. Norm "dataset"
    checkpoints of 2026-10-13 to 2026-10-19 record no `dataset_norm` and z-scored BQ and MAPQ per dataset
    (`DATASET_NORM_BEFORE_20261019`, model.full_config).
"""
import torch
import torch.nn as nn
import torch.nn.functional as F

CATEGORICAL = (("read_base", 0, 6), ("site_allele", 2, 6), ("graph_base", 5, 6), ("operation", 4, 6),
               ("strand", 7, 2))
CONTINUOUS = (("base_quality", 1, 0), ("mapping_quality", 3, 0), ("path_count", 6, 1))  # (name, channel, min valid)
DATASET_NORM = ("base_quality",)  # z-scored per dataset under norm "dataset" (models from 2026-10-19)
DATASET_NORM_BEFORE_20261019 = ("base_quality", "mapping_quality")  # norm "dataset" checkpoints without dataset_norm
BLOCKS = ("A1", "ALT", "REF", "OTHER")
PLANES = tuple([f"{name}={code}" for name, _, k in CATEGORICAL for code in range(1, k + 1)]
               + [name for name, _, _ in CONTINUOUS]
               + ["covered", "bq_missing", "differs"] + [f"row_{b}" for b in BLOCKS])
N_PLANES = len(PLANES)  # 36
WIDTH = 104


def valid_masks(x):
    """{channel name: bool mask of the cells whose value enters the z-score} for int8/int16 x (B, 8, H, W)."""
    covered = x[:, 4] > 0
    return {name: covered & (x[:, c] >= low) for name, c, low in CONTINUOUS}


class TensorEncoder(nn.Module):
    def __init__(self, stats=None, drop=(), dataset_norm=DATASET_NORM):
        super().__init__()
        unknown = set(drop) - set(PLANES)
        if unknown:
            raise ValueError(f"unknown planes {sorted(unknown)}; planes: {', '.join(PLANES)}")
        names = [name for name, _, _ in CONTINUOUS]
        if set(dataset_norm) - set(names):
            raise ValueError(f"dataset_norm {sorted(set(dataset_norm) - set(names))}: not one of {names}")
        self.dataset_norm = tuple(name for name in names if name in set(dataset_norm))  # channel order
        self.drop = tuple(p for p in PLANES if p in set(drop))
        self.keep = [k for k, p in enumerate(PLANES) if p not in set(drop)]
        self.n_planes = len(self.keep)
        self.register_buffer("mean", torch.zeros(len(CONTINUOUS)))
        self.register_buffer("std", torch.ones(len(CONTINUOUS)))
        self.register_buffer("fitted", torch.zeros((), dtype=torch.bool))
        # per-dataset statistics (norm "dataset"): one row per dataset, set by train / predict for the data at hand;
        # not part of the state dict (they belong to the data, not to the model)
        self.register_buffer("part_mean", torch.zeros(0, len(CONTINUOUS)), persistent=False)
        self.register_buffer("part_std", torch.ones(0, len(CONTINUOUS)), persistent=False)
        if stats is not None:
            self.set_stats(stats)

    def set_stats(self, stats):
        """stats: {name: {"mean": m, "std": s, ...}} as returned by compute_stats."""
        for j, (name, _, _) in enumerate(CONTINUOUS):
            self.mean[j] = float(stats[name]["mean"])
            self.std[j] = max(float(stats[name]["std"]), 1e-6)
        self.fitted.fill_(True)

    def stats(self):
        return {name: {"mean": float(self.mean[j]), "std": float(self.std[j])}
                for j, (name, _, _) in enumerate(CONTINUOUS)}

    def set_part_stats(self, part_stats):
        """part_stats: [compute_stats result] per dataset (row p for samples of part p). The `dataset_norm` channels
        take the dataset's mean and std, the others the encoder's (fitted) buffers. An empty list removes the table."""
        if part_stats and not bool(self.fitted):
            raise RuntimeError("set the encoder's statistics (set_stats) before the per-dataset ones")
        mean = torch.zeros(len(part_stats), len(CONTINUOUS))
        std = torch.ones(len(part_stats), len(CONTINUOUS))
        for p, st in enumerate(part_stats):
            for j, (name, _, _) in enumerate(CONTINUOUS):
                if name in self.dataset_norm:
                    mean[p, j], std[p, j] = float(st[name]["mean"]), max(float(st[name]["std"]), 1e-6)
                else:
                    mean[p, j], std[p, j] = self.mean[j], self.std[j]
        self.part_mean = mean.to(self.mean.device)
        self.part_std = std.to(self.std.device)

    def part_stats(self):
        """The per-dataset table as [{name: {"mean", "std"}}], empty without one."""
        return [{name: {"mean": float(self.part_mean[p, j]), "std": float(self.part_std[p, j])}
                 for j, (name, _, _) in enumerate(CONTINUOUS)} for p in range(self.part_mean.shape[0])]

    def forward(self, x, blocks, dtype=torch.float32, part=None):
        """x: (B, 8, H, W) int8; blocks: (B, 4) row ends of the A1, ALT, REF and OTHER blocks; part: (B,) dataset of
        each sample, a row of the per-dataset table (ignored without a table: the buffers for every sample)."""
        # checked only in eager mode: bool() of a GPU buffer is a graph break, and under torch.compile the
        # code after it was recompiled every step until Dynamo gave up (recompile_limit) and ran it eagerly
        if not torch.compiler.is_compiling() and not bool(self.fitted):
            raise RuntimeError("TensorEncoder has no statistics: fit them with compute_stats or load a checkpoint")
        x = x.to(torch.int16)
        covered = x[:, 4] > 0
        planes = [x[:, c] == code for _, c, k in CATEGORICAL for code in range(1, k + 1)]
        valid = valid_masks(x)
        if part is not None and self.part_mean.shape[0]:
            row = part.to(device=x.device, dtype=torch.long)
            mean, std = self.part_mean[row][:, :, None, None], self.part_std[row][:, :, None, None]  # (B, C, 1, 1)
        else:
            mean, std = self.mean[None, :, None, None], self.std[None, :, None, None]
        for j, (name, c, _) in enumerate(CONTINUOUS):
            z = (x[:, c].to(torch.float32) - mean[:, j]) / std[:, j]
            planes.append(torch.where(valid[name], z, torch.zeros((), device=z.device)))
        planes += [covered, covered & (x[:, 1] < 0), covered & (x[:, 0] != x[:, 5])]
        rows = torch.arange(x.shape[2], device=x.device)[None, :]            # (1, H)
        ends = blocks.to(device=x.device, dtype=torch.long)                  # (B, 4)
        starts = torch.cat([torch.zeros_like(ends[:, :1]), ends[:, :-1]], 1)
        for b in range(len(BLOCKS)):
            in_block = (rows >= starts[:, b:b + 1]) & (rows < ends[:, b:b + 1])   # (B, H)
            planes.append(covered & in_block[:, :, None])
        out = torch.stack([p.to(dtype) for p in planes], 1)                 # (B, 36, H, W)
        if self.drop:
            out = out[:, self.keep]
        return F.pad(out, (0, WIDTH - out.shape[-1])) if out.shape[-1] < WIDTH else out


def compute_stats(loader, max_samples):
    """Mean and std of BQ, MAPQ and path count over the valid covered cells of up to max_samples tensors from
    loader (batches whose first item is x); padding and invalid cells (-1) never enter."""
    total = {name: torch.zeros(3, dtype=torch.float64) for name, _, _ in CONTINUOUS}  # count, sum, sum of squares
    seen = 0
    for batch in loader:
        x = batch[0][:max_samples - seen].to(torch.int16)
        for name, mask in valid_masks(x).items():
            c = next(c for n, c, _ in CONTINUOUS if n == name)
            v = x[:, c][mask].to(torch.float64)
            total[name] += torch.stack([torch.tensor(float(v.numel()), dtype=torch.float64), v.sum(), (v * v).sum()])
        seen += x.shape[0]
        if seen >= max_samples:
            break
    stats = {}
    for name, (n, s, ss) in total.items():
        if n == 0:
            raise ValueError(f"no valid cells for {name} in {seen} tensors")
        mean = s / n
        stats[name] = {"mean": float(mean), "std": float(torch.sqrt(torch.clamp(ss / n - mean * mean, min=0))),
                       "cells": int(n)}
    stats["tensors"] = seen
    return stats
