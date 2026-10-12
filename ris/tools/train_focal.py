"""pansoma_net_v2.train with a focal loss in place of the cross entropy (ris experiment; the package is not changed).

    python ris/tools/train_focal.py [train options]      FOCAL_GAMMA (default 2) sets gamma

Focal loss (Lin et al. 2017): -w_y (1 - p_y)^gamma log p_y, with the class weights w of --class-weights as alpha,
normalized like the weighted cross entropy (sum over the batch / sum of w_y), ignore_index -1 as before. gamma 0 is the
cross entropy. train.py builds its losses as nn.CrossEntropyLoss(weight=..., ignore_index=-1[, reduction="none"]);
this script swaps that class in the train module only, then runs train.main. Validation, checkpoints and predict are
unchanged (the validation loss is still the weighted NLL)."""
import json
import os
import sys
import types

import torch
import torch.nn as tnn

import pansoma_net_v2.train as T

GAMMA = float(os.environ.get("FOCAL_GAMMA", "2"))


class FocalLoss(tnn.Module):
    def __init__(self, weight=None, ignore_index=-100, reduction="mean"):
        super().__init__()
        self.weight, self.ignore_index, self.reduction = weight, ignore_index, reduction

    def forward(self, logits, target):
        valid = target != self.ignore_index
        y = torch.where(valid, target, torch.zeros_like(target))
        logp = torch.log_softmax(logits, dim=1).gather(1, y[:, None]).squeeze(1)
        w = self.weight[y] if self.weight is not None else torch.ones_like(logp)
        w = torch.where(valid, w, torch.zeros_like(w))
        loss = -w * (1.0 - logp.exp()).clamp_min(1e-6).pow(GAMMA) * logp  # clamp: finite gradient for gamma < 1
        if self.reduction == "none":
            return loss
        return loss.sum() / w.sum().clamp_min(1e-12)


shim = types.ModuleType("nn_focal")
shim.__dict__.update(tnn.__dict__)
shim.CrossEntropyLoss = FocalLoss
T.nn = shim

if __name__ == "__main__":
    argv = sys.argv[1:]
    out = argv[argv.index("--output") + 1] if "--output" in argv else None
    if out:
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "focal.json"), "w") as f:
            json.dump({"loss": "focal", "gamma": GAMMA}, f)
    print(f"focal loss, gamma {GAMMA}", flush=True)
    T.main(argv)
