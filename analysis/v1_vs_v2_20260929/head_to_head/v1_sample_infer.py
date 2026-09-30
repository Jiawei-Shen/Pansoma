"""Stratified-sample CPU inference of the published v1 HG008_WGS_SNV.pth on v1 chr1 val tensors (read-only inputs).

Strata (inside the confident region = somatic BED ∩ germline BED, v1 AF floor 0.08 is intrinsic):
  T: every v1 tensor labeled true (670)
  G: v1 false tensors at a HG008-N germline SNP truth allele (23,015)  -> random sample
  N: other v1 false tensors (224,882)                                  -> random sample
usage: v1_sample_infer.py <part 1|2> ; part 1 = T + first halves of G/N samples, part 2 = second halves.
"""
import os
import pickle
import sys
import time

import numpy as np
import torch

HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from v1_timing import load_model, MEAN, STD, VAL  # noqa

torch.set_num_threads(4)
N_G, N_N = 3000, 4500
part = int(sys.argv[1])

recs = pickle.load(open(f"{HERE}/v1_val_recs.pkl", "rb"))  # (pos0, alt, lab, inconf, isgerm, sidx, iws, af)
rows = np.arange(len(recs))
lab = np.array([r[2] for r in recs])
inconf = np.array([r[3] for r in recs])
isgerm = np.array([r[4] for r in recs])
T = rows[lab == 1]
G_all = rows[(lab == 0) & inconf & isgerm]
N_all = rows[(lab == 0) & inconf & ~isgerm]
rng = np.random.default_rng(20260929)
G = np.sort(rng.choice(G_all, N_G, replace=False))
Nn = np.sort(rng.choice(N_all, N_N, replace=False))
if part == 1:
    sel = np.concatenate([T, G[: N_G // 2], Nn[: N_N // 2]])
else:
    sel = np.concatenate([G[N_G // 2:], Nn[N_N // 2:]])
strata = np.array(["T" if lab[i] == 1 else ("G" if isgerm[i] else "N") for i in sel])
order = np.argsort(sel)
sel, strata = sel[order], strata[order]
print(f"part {part}: {len(sel)} tensors; strata sizes G_all={len(G_all)} N_all={len(N_all)} T={len(T)}", flush=True)

t0 = time.time()
model = load_model()
print("model loaded", round(time.time() - t0, 1), flush=True)
shards = {}
probs = np.zeros(len(sel), dtype=np.float32)
BS = 16
t0 = time.time()
with torch.inference_mode():
    for b in range(0, len(sel), BS):
        idx = sel[b:b + BS]
        xs = []
        for i in idx:
            s, k = divmod(int(i), 32768)
            if s not in shards:
                shards[s] = np.load(f"{VAL}/HG008T_Illumina_chr1_{s:05d}_data.npy", mmap_mode="r")
            xs.append(np.asarray(shards[s][k]))
        xb = (torch.from_numpy(np.stack(xs)).float() - MEAN) / STD
        probs[b:b + BS] = torch.softmax(model(xb), 1)[:, 1].numpy()
        if (b // BS) % 50 == 0:
            el = time.time() - t0
            print(f"{b + len(idx)}/{len(sel)} {el:.0f}s eta {el / (b + len(idx)) * (len(sel) - b - len(idx)):.0f}s", flush=True)
np.savez(f"{HERE}/v1_sample_part{part}.npz", rows=sel, strata=strata, p_true=probs,
         n_G_all=len(G_all), n_N_all=len(N_all), n_G=N_G, n_N=N_N)
print("done", round(time.time() - t0, 1), flush=True)
