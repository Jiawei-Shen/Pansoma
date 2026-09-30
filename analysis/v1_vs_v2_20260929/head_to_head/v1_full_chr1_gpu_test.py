#!/usr/bin/env python3
"""FULL v1 inference (GPU) on every v1 chr1 SNV *testing* tensor (361,394 = 345,376 GRCh38 + 16,018 off-reference),
for the exact truth-level curve at every threshold. NOT RUN HERE (CPU 4 threads: ~7-9 tensors/s -> ~12 h).

Inputs (read-only):
  checkpoint  default: /scratch/jshen/Github/Pansoma/pretrained_model/Pansoma_zenodo/pretrained_model_Pansoma/HG008/HG008_WGS_SNV.pth
              (byte-identical to /scratch/jshen/Github/GoogleNet/HG008T_GIAB_AF-HPRC_CE_Large_Model_V2_weight200/model_e067_f1_0.1784.pth)
  tensors     /scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV/5ch_testing_data_SNV_chr1/shard_000{00..11}_data.npy
              (int8 (N,5,201,100); row order = variant_summary.ndjson order = shard*32768 + index_within_shard)
Output: <out>/v1_chr1_test_p_true.<tag>.npy (float32, 361,394, variant_summary.ndjson order)

Normalization = test_5channels_npy_pansoma.py VAL_MEAN/VAL_STD. CPU check in this session: the same code on the e053
checkpoint reproduced the PROB of the logged May-14 GPU VCF to within 1.3e-4 (v1_infer_check.out).
"""
import argparse, glob, os, sys, time
import numpy as np
import torch
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net")
from mynet import ConvNeXtCBAMClassifier  # noqa

TEST = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV/5ch_testing_data_SNV_chr1"
CKPT = "/scratch/jshen/Github/Pansoma/pretrained_model/Pansoma_zenodo/pretrained_model_Pansoma/HG008/HG008_WGS_SNV.pth"
MEAN = torch.tensor([18.41781616, 12.64912987, -0.54525274, 24.72385406, 4.69061136]).view(1, 5, 1, 1)
STD = torch.tensor([25.02832222, 14.80963230, 0.61813378, 29.97283554, 7.92317915]).view(1, 5, 1, 1)
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--ckpt", default=CKPT)
ap.add_argument("--tag", default="e067")
ap.add_argument("--batch", type=int, default=256)
ap.add_argument("--fp16", action="store_true", help="autocast fp16 (what test_5channels_npy_pansoma.py does on CUDA)")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
dev = torch.device("cuda")
ck = torch.load(a.ckpt, map_location="cpu", weights_only=False)
model = ConvNeXtCBAMClassifier(in_channels=5, class_num=2, depths=[3, 3, 27, 3], dims=[192, 384, 768, 1536])
model.load_state_dict(ck["model_state_dict"], strict=True)
model.eval().to(dev)
print("checkpoint", a.ckpt, "epoch", ck.get("epoch"), flush=True)
mean, std = MEAN.to(dev), STD.to(dev)
out, t0 = [], time.time()
shards = sorted(glob.glob(f"{TEST}/shard_*_data.npy"))
assert len(shards) == 12, shards
with torch.inference_mode():
    for s, path in enumerate(shards):
        x = np.load(path, mmap_mode="r")
        for b in range(0, len(x), a.batch):
            xb = (torch.from_numpy(np.ascontiguousarray(x[b:b + a.batch])).to(dev).float() - mean) / std
            with torch.autocast("cuda", dtype=torch.float16, enabled=a.fp16):
                p = torch.softmax(model(xb).float(), 1)[:, 1]
            out.append(p.cpu().numpy())
        print(f"{os.path.basename(path)} done, {sum(map(len, out))} rows, {time.time() - t0:.0f}s", flush=True)
p = np.concatenate(out).astype(np.float32)
assert len(p) == 361394, len(p)
np.save(f"{a.out}/v1_chr1_test_p_true.{a.tag}.npy", p)
print("saved", len(p), f"{time.time() - t0:.0f}s")
