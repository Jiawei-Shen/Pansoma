"""Experiment (2026-10-13): BQ normalized per read, without changing the package. Runs pansoma_net_v2.train or
.predict with the encoder's BQ plane replaced by

    (BQ - median BQ of that read's valid covered cells in the 101-column window) / std

where std is the BQ std of the sample's dataset (--norm dataset, encode.py) or the training one. A run offset of a
whole read set (the tumor and normal reads of an in-silico mixture: Illumina mean BQ 29.43 vs 30.77) cancels; a single
low-quality base inside a read stays; a wholly low-quality read looks like an average one (known limitation). Every
other plane, MAPQ (one value per read, z-scored per dataset) and path count included, is encode.py's.

    python readnorm.py train   <pansoma_net_v2.train options>        (writes <output>/readnorm.json)
    torchrun --standalone --nproc_per_node=N readnorm.py train ... --ddp
    python readnorm.py predict <pansoma_net_v2.predict options>      (needs readnorm.json beside the checkpoint)
    python readnorm.py selftest
Run from machine_learning/ with it on PYTHONPATH (ris/joint_run.sh --readnorm does both).
"""
import json
import os
import sys
from pathlib import Path

import torch

import pansoma_net_v2.encode as E

BQ_PLANE = E.PLANES.index("base_quality")
BQ_CHANNEL, BQ_ROW = dict((n, c) for n, c, _ in E.CONTINUOUS)["base_quality"], [n for n, _, _ in E.CONTINUOUS].index("base_quality")
_forward = E.TensorEncoder.forward


def read_centered_bq(self, x, part=None):
    """(B, H, W) float32: BQ minus the median of its read's valid cells, over the dataset's (or training) BQ std; 0 on
    invalid cells and padding."""
    x = x.to(torch.int16)
    valid = (x[:, 4] > 0) & (x[:, BQ_CHANNEL] >= 0)
    bq = x[:, BQ_CHANNEL].to(torch.float32)
    med = torch.nanmedian(torch.where(valid, bq, torch.full_like(bq, float("nan"))), dim=-1, keepdim=True).values
    if part is not None and self.part_std.shape[0]:
        std = self.part_std[part.to(device=x.device, dtype=torch.long), BQ_ROW][:, None, None]
    else:
        std = self.std[BQ_ROW]
    return torch.where(valid, (bq - med) / std, torch.zeros((), device=bq.device))


def forward(self, x, blocks, dtype=torch.float32, part=None):
    out = _forward(self, x, blocks, dtype, part)
    if BQ_PLANE in self.keep:
        k = self.keep.index(BQ_PLANE)
        z = read_centered_bq(self, x, part)
        out = out.clone()
        out[:, k, :, :z.shape[-1]] = z.to(out.dtype)
    return out


def patch():
    E.TensorEncoder.forward = forward


def selftest():
    """The patched BQ plane against a direct computation; every other plane unchanged."""
    import numpy as np
    from pansoma_net_v2.tests.test_encode import STATS, batch
    x, blocks = batch(8, seed=3)
    enc = E.TensorEncoder(STATS)
    plain = _forward(enc, x, blocks)
    patched = forward(enc, x, blocks)
    other = [k for k in range(plain.shape[1]) if k != BQ_PLANE]
    assert torch.equal(plain[:, other], patched[:, other]), "other planes changed"
    xn = x.numpy().astype(np.float64)
    valid = (xn[:, 4] > 0) & (xn[:, 1] >= 0)
    for b in range(len(xn)):
        for r in range(xn.shape[2]):
            v = valid[b, r]
            want = np.zeros(xn.shape[3])
            if v.any():
                vals = np.sort(xn[b, 1, r, v])
                want[v] = (xn[b, 1, r, v] - vals[(len(vals) - 1) // 2]) / STATS["base_quality"]["std"]  # lower median
            got = patched[b, BQ_PLANE, r, :xn.shape[3]].numpy()
            assert np.allclose(got, want, atol=1e-5), (b, r)
    assert float(patched[:, BQ_PLANE, :, xn.shape[3]:].abs().max()) == 0.0, "padding not zero"
    print("readnorm selftest OK: BQ plane read-centered (lower median), other planes unchanged, padding zero")


def main():
    mode, argv = sys.argv[1], sys.argv[2:]
    if mode == "selftest":
        return selftest()
    patch()
    if mode == "train":
        from pansoma_net_v2 import train
        out = Path(argv[argv.index("--output") + 1])
        if int(os.environ.get("RANK", "0")) == 0:
            out.mkdir(parents=True, exist_ok=True)
            (out / "readnorm.json").write_text(json.dumps(
                {"bq": "read_centered_median", "std": "dataset (--norm dataset) or training", "script": __file__}) + "\n")
        train.main(argv)
    elif mode == "predict":
        from pansoma_net_v2 import predict
        ckpt = Path(argv[argv.index("--checkpoint") + 1])
        if not (ckpt.parent / "readnorm.json").exists():
            raise SystemExit(f"{ckpt.parent}/readnorm.json missing: not a readnorm run")
        predict.main(argv)
    else:
        raise SystemExit(f"mode must be train, predict or selftest, not {mode}")


if __name__ == "__main__":
    main()
