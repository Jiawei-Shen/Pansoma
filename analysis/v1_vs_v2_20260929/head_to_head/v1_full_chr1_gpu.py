#!/usr/bin/env python3
"""FULL v1 inference (GPU) on every v1 chr1 SNV val tensor, for the exact truth-level curve. NOT RUN HERE (CPU: ~9 h).

Inputs (read-only):
  checkpoint  /scratch/jshen/Github/Pansoma/pretrained_model/Pansoma_zenodo/pretrained_model_Pansoma/HG008/HG008_WGS_SNV.pth
              (= GoogleNet/HG008T_GIAB_AF-HPRC_CE_Large_Model_V2_weight200/model_e067_f1_0.1784.pth; epoch 67)
  tensors     /scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/
              HG008T_Illumina_chr1_000{00..10}_data.npy  (345,376 x 5 x 201 x 100 int8; row = shard*32768 + index)
Output: <out>/v1_chr1_val_p_true.npy (float32, 345,376, same row order as variant_summary_classified.ndjson)

Run (one GPU, ~10-20 min on an A100/H100; ~1 GB GPU memory for bs 256 fp16):
  cd /scratch/jshen/Github/Pansoma/machine_learning/pansoma_net
  python -u <scratch>/v1_full_chr1_gpu.py --out <scratch>/v1_full --batch 256
then:
  python <scratch>/combine.py --v1-full <scratch>/v1_full/v1_chr1_val_p_true.npy HG008_Illumina_SNV_base ...
  (or load the npy in place of the stratified sample: every row, weight 1)
"""
import argparse
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net")
from mynet import ConvNeXtCBAMClassifier  # noqa

CKPT = "/scratch/jshen/Github/Pansoma/pretrained_model/Pansoma_zenodo/pretrained_model_Pansoma/HG008/HG008_WGS_SNV.pth"
VAL = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset"
MEAN = torch.tensor([18.41781616, 12.64912987, -0.54525274, 24.72385406, 4.69061136]).view(1, 5, 1, 1)
STD = torch.tensor([25.02832222, 14.80963230, 0.61813378, 29.97283554, 7.92317915]).view(1, 5, 1, 1)

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--batch", type=int, default=256)
ap.add_argument("--fp16", action="store_true", help="autocast fp16 (the v1 test script's CUDA default); off = fp32")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
dev = torch.device("cuda")
ck = torch.load(CKPT, map_location="cpu", mmap=True, weights_only=False)
model = ConvNeXtCBAMClassifier(in_channels=5, class_num=2, depths=[3, 3, 27, 3], dims=[192, 384, 768, 1536])
model.load_state_dict(ck["model_state_dict"], strict=True)
model.eval().to(dev)
mean, std = MEAN.to(dev), STD.to(dev)
out = []
t0 = time.time()
with torch.inference_mode():
    for s in range(11):
        x = np.load(f"{VAL}/HG008T_Illumina_chr1_{s:05d}_data.npy", mmap_mode="r")
        for b in range(0, len(x), a.batch):
            xb = (torch.from_numpy(np.ascontiguousarray(x[b:b + a.batch])).to(dev).float() - mean) / std
            with torch.autocast("cuda", dtype=torch.float16, enabled=a.fp16):
                p = torch.softmax(model(xb).float(), 1)[:, 1]
            out.append(p.cpu().numpy())
        print(f"shard {s} done, {sum(map(len, out))} rows, {time.time() - t0:.0f}s", flush=True)
p = np.concatenate(out).astype(np.float32)
assert len(p) == 345376, len(p)
np.save(f"{a.out}/v1_chr1_val_p_true.npy", p)
print("saved", len(p))
