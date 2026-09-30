"""Part 3 of the stratified sample (published HG008_WGS_SNV.pth = e067): 1,000 more G and 3,000 more N rows drawn
at random from the rows not sampled in parts 1-2 (union stays a simple random sample per stratum). Read-only inputs."""
import pickle, sys, time
import numpy as np
import torch
HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from v1_timing import load_model, MEAN, STD, VAL  # noqa
torch.set_num_threads(4)
recs = pickle.load(open(f"{HERE}/v1_val_recs.pkl", "rb"))
rows = np.arange(len(recs))
lab = np.array([r[2] for r in recs]); inconf = np.array([r[3] for r in recs]); isgerm = np.array([r[4] for r in recs])
done = np.concatenate([np.load(f"{HERE}/v1_sample_part{k}.npz")["rows"] for k in (1, 2)])
G_all = rows[(lab == 0) & inconf & isgerm]; N_all = rows[(lab == 0) & inconf & ~isgerm]
rng = np.random.default_rng(20260930)
G = rng.choice(np.setdiff1d(G_all, done), 1000, replace=False)
Nn = rng.choice(np.setdiff1d(N_all, done), 3000, replace=False)
sel = np.sort(np.concatenate([G, Nn]))
strata = np.array(["G" if isgerm[i] else "N" for i in sel])
print("part 3:", len(sel), flush=True)
model = load_model()
shards, probs, BS, t0 = {}, np.zeros(len(sel), np.float32), 16, time.time()
with torch.inference_mode():
    for b in range(0, len(sel), BS):
        xs = []
        for i in sel[b:b + BS]:
            s, k = divmod(int(i), 32768)
            if s not in shards:
                shards[s] = np.load(f"{VAL}/HG008T_Illumina_chr1_{s:05d}_data.npy", mmap_mode="r")
            xs.append(np.asarray(shards[s][k]))
        xb = (torch.from_numpy(np.stack(xs)).float() - MEAN) / STD
        probs[b:b + BS] = torch.softmax(model(xb), 1)[:, 1].numpy()
        if (b // BS) % 50 == 0:
            el = time.time() - t0
            print(f"{b + BS}/{len(sel)} {el:.0f}s", flush=True)
np.savez(f"{HERE}/v1_sample_part3.npz", rows=sel, strata=strata, p_true=probs,
         n_G_all=len(G_all), n_N_all=len(N_all), n_G=1000, n_N=3000)
print("done", round(time.time() - t0, 1), flush=True)
