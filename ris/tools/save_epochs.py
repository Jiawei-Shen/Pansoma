"""Train exactly as `python -m pansoma_net_v2.train`, and keep every epoch's checkpoint: each time train writes
last.pth (once per epoch, with that epoch's validation report and thresholds), a copy goes to
<output>/epochs/epoch_NN.pth. For checkpoint-selection experiments (2026-10-19: COLO829T Illumina d9, where two
near-tied epochs gave PASS F1 0.778 and 0.624). Usage: the train options, e.g. from machine_learning/:
    python <repo>/ris/tools/save_epochs.py --tensors <set>/tensors --kinds SNV --output <run> ...
--stop-after-epoch N (taken off before the train options) ends the run once epoch N is saved, on every rank as the
next epoch would start, with --epochs and so the learning-rate schedule unchanged: epochs 1-N are those of the full
run (2026-10-20: epoch 6 of the COLO829T Illumina d9 --norm dataset seeds). Under torchrun with --ddp as run.sh.
"""
import shutil
import sys
from pathlib import Path

import torch

from pansoma_net_v2 import data, train

_save = torch.save


def save(obj, f, *args, **kwargs):
    _save(obj, f, *args, **kwargs)
    path = Path(f) if isinstance(f, (str, Path)) else None
    if path is not None and path.name == "last.pth.tmp" and isinstance(obj, dict) and "epoch" in obj:
        keep = path.parent / "epochs"
        keep.mkdir(exist_ok=True)
        shutil.copyfile(path, keep / f"epoch_{int(obj['epoch']):02d}.pth")


torch.save = save  # train calls torch.save(payload, <output>/last.pth.tmp) then renames it


def stop_after(n):
    set_epoch = data.EpochSampler.set_epoch

    def stop(self, epoch):  # train: train_sampler.set_epoch(epoch) first in every epoch (0-based), on every rank
        if self.shuffle and epoch >= n:
            print(f"stop: epoch {n} saved (--stop-after-epoch {n})", flush=True)
            if torch.distributed.is_initialized():
                torch.distributed.destroy_process_group()
            raise SystemExit(0)
        set_epoch(self, epoch)

    data.EpochSampler.set_epoch = stop


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--stop-after-epoch" in argv:
        i = argv.index("--stop-after-epoch")
        stop_after(int(argv[i + 1]))
        del argv[i:i + 2]
    train.main(argv)
