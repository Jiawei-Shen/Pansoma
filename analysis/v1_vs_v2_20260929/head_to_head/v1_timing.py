"""Load the published v1 HG008_WGS_SNV.pth on CPU and time forward passes on real v1 chr1 val tensors (read-only)."""
import os
import sys
import time

import numpy as np
import torch

torch.set_num_threads(4)
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net")
from mynet import ConvNeXtCBAMClassifier  # noqa

CKPT = "/scratch/jshen/Github/Pansoma/pretrained_model/Pansoma_zenodo/pretrained_model_Pansoma/HG008/HG008_WGS_SNV.pth"
VAL = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset"
MEAN = torch.tensor([18.41781616, 12.64912987, -0.54525274, 24.72385406, 4.69061136]).view(1, 5, 1, 1)
STD = torch.tensor([25.02832222, 14.80963230, 0.61813378, 29.97283554, 7.92317915]).view(1, 5, 1, 1)


def load_model():
    ck = torch.load(CKPT, map_location="cpu", mmap=True, weights_only=False)
    m = ConvNeXtCBAMClassifier(in_channels=5, class_num=2, depths=[3, 3, 27, 3], dims=[192, 384, 768, 1536])
    m.load_state_dict(ck["model_state_dict"], strict=True)
    m.eval()
    return m


if __name__ == "__main__":
    t = time.time()
    model = load_model()
    print("model load s", round(time.time() - t, 1), flush=True)
    x = np.load(f"{VAL}/HG008T_Illumina_chr1_00000_data.npy", mmap_mode="r")
    for bs in (16, 32):
        xb = torch.from_numpy(np.ascontiguousarray(x[:bs])).float()
        xb = (xb - MEAN) / STD
        with torch.inference_mode():
            t = time.time()
            out = model(xb)
            dt = time.time() - t
        print("batch", bs, "s", round(dt, 2), "tensors/s", round(bs / dt, 1), flush=True)
    print(torch.softmax(out, 1)[:5])
