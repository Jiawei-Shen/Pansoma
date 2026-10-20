"""Train exactly as `python -m pansoma_net_v2.train`, and keep every epoch's checkpoint: each time train writes
last.pth (once per epoch, with that epoch's validation report and thresholds), a copy goes to
<output>/epochs/epoch_NN.pth. For checkpoint-selection experiments (2026-10-19: COLO829T Illumina d9, where two
near-tied epochs gave PASS F1 0.778 and 0.624). Usage: the train options, e.g. from machine_learning/:
    python <repo>/ris/tools/save_epochs.py --tensors <set>/tensors --kinds SNV --output <run> ...
"""
import shutil
import sys
from pathlib import Path

import torch

from pansoma_net_v2 import train

_save = torch.save


def save(obj, f, *args, **kwargs):
    _save(obj, f, *args, **kwargs)
    path = Path(f) if isinstance(f, (str, Path)) else None
    if path is not None and path.name == "last.pth.tmp" and isinstance(obj, dict) and "epoch" in obj:
        keep = path.parent / "epochs"
        keep.mkdir(exist_ok=True)
        shutil.copyfile(path, keep / f"epoch_{int(obj['epoch']):02d}.pth")


torch.save = save  # train calls torch.save(payload, <output>/last.pth.tmp) then renames it

if __name__ == "__main__":
    train.main(sys.argv[1:])
