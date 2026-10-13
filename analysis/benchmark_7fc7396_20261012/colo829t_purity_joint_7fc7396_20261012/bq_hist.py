"""Base-quality value histograms of COLO829T SNV tensors per purity (first 2,000 tensors of chr2 shard 0, covered
cells, BQ >= 0): does the normal (COLO829BL) read set carry a different BQ alphabet than the tumor set?
Run in the container: ris/submit.sh -g 0 -n 2 -m 16 NAME py.sh $PWD/analysis/benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012/bq_hist.py"""
import numpy as np, sys, json
D = "/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors"
for plat in ("Illumina", "fiberseq"):
    for q in ("", "_P75", "_P50", "_P25"):
        s = f"COLO829T_{plat}{q}"
        a = np.load(f"{D}/{s}/tensors/SNV/chr2_shard_00000_data.npy", mmap_mode="r")
        x = np.asarray(a[:2000, 1])  # channel 1 = base quality, 2000 tensors
        cov = np.asarray(a[:2000, 4]) > 0  # covered cells (operation channel)
        v = x[cov & (x >= 0)]
        u, c = np.unique(v, return_counts=True)
        top = sorted(zip(c, u), reverse=True)[:6]
        print(f"{s:24s} distinct BQ values {len(u):3d}; top: " + ", ".join(f"{int(b)}:{n / c.sum():.3f}" for n, b in top))
