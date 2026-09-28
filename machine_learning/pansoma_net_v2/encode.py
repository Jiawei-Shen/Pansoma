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

Every plane is 0 where there is no evidence, so padding stays 0 after encoding. The z-score statistics are
fitted once on covered cells of the training set (`compute_stats`) and live in the encoder's buffers, i.e. in
the model checkpoint; they are never refitted on other data. The width is zero-padded from 101 to 104 so the
4 x 4, stride 4 stem covers the last column.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F

CATEGORICAL = (("read_base", 0, 6), ("site_allele", 2, 6), ("graph_base", 5, 6), ("operation", 4, 6),
               ("strand", 7, 2))
CONTINUOUS = (("base_quality", 1, 0), ("mapping_quality", 3, 0), ("path_count", 6, 1))  # (name, channel, min valid)
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
    def __init__(self, stats=None):
        super().__init__()
        self.register_buffer("mean", torch.zeros(len(CONTINUOUS)))
        self.register_buffer("std", torch.ones(len(CONTINUOUS)))
        self.register_buffer("fitted", torch.zeros((), dtype=torch.bool))
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

    def forward(self, x, blocks, dtype=torch.float32):
        """x: (B, 8, H, W) int8; blocks: (B, 4) row ends of the A1, ALT, REF and OTHER blocks."""
        if not bool(self.fitted):
            raise RuntimeError("TensorEncoder has no statistics: fit them with compute_stats or load a checkpoint")
        x = x.to(torch.int16)
        covered = x[:, 4] > 0
        planes = [x[:, c] == code for _, c, k in CATEGORICAL for code in range(1, k + 1)]
        valid = valid_masks(x)
        for j, (name, c, _) in enumerate(CONTINUOUS):
            z = (x[:, c].to(torch.float32) - self.mean[j]) / self.std[j]
            planes.append(torch.where(valid[name], z, torch.zeros((), device=z.device)))
        planes += [covered, covered & (x[:, 1] < 0), covered & (x[:, 0] != x[:, 5])]
        rows = torch.arange(x.shape[2], device=x.device)[None, :]            # (1, H)
        ends = blocks.to(device=x.device, dtype=torch.long)                  # (B, 4)
        starts = torch.cat([torch.zeros_like(ends[:, :1]), ends[:, :-1]], 1)
        for b in range(len(BLOCKS)):
            in_block = (rows >= starts[:, b:b + 1]) & (rows < ends[:, b:b + 1])   # (B, H)
            planes.append(covered & in_block[:, :, None])
        out = torch.stack([p.to(dtype) for p in planes], 1)                 # (B, 36, H, W)
        return F.pad(out, (0, WIDTH - out.shape[-1])) if out.shape[-1] < WIDTH else out


def compute_stats(loader, max_samples):
    """Mean and std of BQ, MAPQ and path count over the valid covered cells of up to max_samples tensors from
    loader (batches of (x, blocks, label)); padding and invalid cells (-1) never enter."""
    total = {name: torch.zeros(3, dtype=torch.float64) for name, _, _ in CONTINUOUS}  # count, sum, sum of squares
    seen = 0
    for x, _, _ in loader:
        x = x[:max_samples - seen].to(torch.int16)
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
