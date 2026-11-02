"""Q3: GRCh38 intervals of the validation node blocks (train.py block_split, seed s, fraction 0.05, 20,000-node
blocks) of chr2-22, per tensor set (SNV index cache), and their overlap between graphs and seeds.

Block span (as data.validation_truth assigns truth): each block starts at the smallest pos0 of its tensors; a GRCh38
position belongs to the last block starting at or before it. The last block of a chromosome is capped at the largest
pos0 + 1. Also reported: the 'hull' definition [min pos0, max pos0 + 1) of each validation block's tensors.
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np

CACHE = Path("/scratch/jshen/data/pansoma_net_v2_runs/index_cache")
FILES = {
    ("HG008T_Illumina", "linear"): "HG008T_Illumina_linear.tensors.SNV.df0b6f07d683.npz",
    ("HG008T_Illumina", "d9"): "HG008T_Illumina.tensors.SNV.0062217f8edf.npz",
    ("HG008T_Illumina", "d46"): "HG008T_Illumina_d46.tensors.SNV.91645eef0a6a.npz",
    ("HG008T_PacBio", "linear"): "HG008T_PacBio_linear.tensors.SNV.bdf4b4d4fd76.npz",
    ("HG008T_PacBio", "d9"): "HG008T_PacBio.tensors.SNV.63ee02d2a0bc.npz",
    ("HG008T_PacBio", "d46"): "HG008T_PacBio_d46.tensors.SNV.98441ebf0fb7.npz",
    ("HG008T_ONT", "linear"): "HG008T_ONT_linear.tensors.SNV.21cc6b8fc0f9.npz",
    ("HG008T_ONT", "d9"): "HG008T_ONT.tensors.SNV.91cf0ad9502b.npz",
    ("HG008T_ONT", "d46"): "HG008T_ONT_d46.tensors.SNV.4da62ace82a9.npz",
    ("COLO829T_Illumina", "linear"): "COLO829T_Illumina_linear.tensors.SNV.ecb5c1be7fbe.npz",
    ("COLO829T_Illumina", "d9"): "COLO829T_Illumina.tensors.SNV.0428e5b605a2.npz",
    ("COLO829T_Illumina", "d46"): "COLO829T_Illumina_d46.tensors.SNV.d651e40b0e9a.npz",
    ("COLO829T_fiberseq", "linear"): "COLO829T_fiberseq_linear.tensors.SNV.03675fba739d.npz",
    ("COLO829T_fiberseq", "d9"): "COLO829T_fiberseq.tensors.SNV.3b5c36659e1e.npz",
    ("COLO829T_fiberseq", "d46"): "COLO829T_fiberseq_d46.tensors.SNV.038636730ce3.npz",
    ("COLO829T_ONT", "linear"): "COLO829T_ONT_linear.tensors.SNV.f772e110510c.npz",
    ("COLO829T_ONT", "d9"): "COLO829T_ONT.tensors.SNV.d910a5c87a69.npz",
    ("COLO829T_ONT", "d46"): "COLO829T_ONT_d46.tensors.SNV.46029440f280.npz",
}
SEEDS, FRACTION, BLOCK = (0, 1, 2), 0.05, 20000
CHROMS = [f"chr{i}" for i in range(2, 23)]


def block_u(seed, chrom, block):  # = data._block_u
    return int.from_bytes(hashlib.sha1(f"{seed}:{chrom}:{block}".encode()).digest()[:8], "big") / 2 ** 64


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def total(iv):
    return sum(e - s for s, e in iv)


def inter(a, b):
    i = j = 0
    t = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            t += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t


BIN = 10000
BINS = {}


def spans(npz):
    z = np.load(CACHE / npz)
    meta = json.loads(str(z["meta"]))
    names = meta["chroms"]
    chrom, node, pos0, off = z["chrom"], z["node"], z["pos0"], z["off_reference"]
    BINS[npz] = {}
    for c in CHROMS:
        if c in names:
            m = (chrom == names.index(c)) & (~off) & (pos0 >= 0)
            BINS[npz][c] = (node[m] // BLOCK, pos0[m] // BIN)
    out = {}
    for c in CHROMS:
        if c not in names:
            continue
        m = (chrom == names.index(c)) & (pos0 >= 0)
        b, p = node[m] // BLOCK, pos0[m]
        order = np.lexsort((p, b))
        b, p = b[order], p[order]
        ub, first = np.unique(b, return_index=True)
        last = np.r_[first[1:], len(b)] - 1
        out[c] = dict(block=ub, start=p[first], hull_end=p[last] + 1, end_cap=int(p.max()) + 1,
                      n=np.diff(np.r_[first, len(b)]))
    return out, int(len(chrom))


def val_intervals(sp, seed):
    """{chrom: merged [start, end) list} for the 'tiling' and 'hull' definitions, and block counts."""
    tile, hull, nval, nblk = {}, {}, 0, 0
    for c, d in sp.items():
        isval = np.array([block_u(seed, c, int(k)) < FRACTION for k in d["block"]])
        nval += int(isval.sum()); nblk += len(isval)
        o = np.argsort(d["start"], kind="stable")
        st = d["start"][o]
        en = np.r_[st[1:], d["end_cap"]]
        v = isval[o]
        tile[c] = merge([[int(s), int(e)] for s, e, x in zip(st, en, v) if x and e > s])
        hull[c] = merge([[int(s), int(e)] for s, e, x in zip(d["start"], d["hull_end"], isval) if x])
    return tile, hull, nval, nblk


res = {"sets": {}, "pairs": []}
ivs = {}
for (s, g), f in FILES.items():
    sp, n = spans(f)
    covered = sum(int(d["end_cap"] - d["start"].min()) for d in sp.values())
    for seed in SEEDS:
        tile, hull, nval, nblk = val_intervals(sp, seed)
        # bins: a 10-kb GRCh38 bin is validation when > 50 % of its on-reference tensors are in validation blocks
        binv = {}
        ntens = nvt = 0
        for c, (blk, b) in BINS[f].items():
            isv = np.array([block_u(seed, c, int(k)) < FRACTION for k in np.unique(blk)])
            u, inv = np.unique(blk, return_inverse=True)
            tv = isv[inv]
            ntens += len(tv); nvt += int(tv.sum())
            ub, bi = np.unique(b, return_inverse=True)
            frac = np.bincount(bi, weights=tv) / np.bincount(bi)
            binv[c] = (ub, frac > 0.5, frac)
        tile = {c: v for c, v in tile.items()}
        valbins = {c: [[int(x) * BIN, int(x) * BIN + BIN] for x in ub[m]] for c, (ub, m, _) in binv.items()}
        mixed = sum(int(((fr > 0) & (fr < 1)).sum()) for _, (_, _, fr) in binv.items())
        nb = sum(len(ub) for _, (ub, _, _) in binv.items())
        ivs[(s, g, seed)] = {"tile": tile, "hull": hull, "bins": {c: merge(v) for c, v in valbins.items()}}
        res["sets"][f"{s}|{g}|{seed}"] = dict(blocks=nblk, val_blocks=nval, val_bp_tile=sum(total(v) for v in tile.values()),
                                               val_bp_bins=sum(total(v) for v in ivs[(s, g, seed)]["bins"].values()),
                                               bins_with_tensors=nb, bins_mixed=mixed, onref_tensors=ntens, onref_val_tensors=nvt,
                                               val_bp_hull=sum(total(v) for v in hull.values()), span_bp=covered,
                                               tensors=n, npz=f)
        print(s, g, seed, res["sets"][f"{s}|{g}|{seed}"], flush=True)


def pair(a, b, which):
    A, B = ivs[a][which], ivs[b][which]
    i = sum(inter(A.get(c, []), B.get(c, [])) for c in CHROMS)
    ta, tb = sum(total(v) for v in A.values()), sum(total(v) for v in B.values())
    return dict(a=list(a), b=list(b), def_=which, a_bp=ta, b_bp=tb, inter_bp=i, jaccard=i / (ta + tb - i) if ta + tb - i else None,
                frac_a_in_b=i / ta if ta else None)


sets = sorted({k[0] for k in FILES})
for s in sets:
    for seed in SEEDS:
        for g1, g2 in (("linear", "d9"), ("linear", "d46"), ("d9", "d46")):
            for w in (0, 1, 2):
                res["pairs"].append(pair((s, g1, seed), (s, g2, seed), ("tile", "hull", "bins")[w]) | {"kind": "graphs, same seed"})
    for g in ("linear", "d9", "d46"):
        for s1, s2 in ((0, 1), (0, 2), (1, 2)):
            for w in (0, 1, 2):
                res["pairs"].append(pair((s, g, s1), (s, g, s2), ("tile", "hull", "bins")[w]) | {"kind": "seeds, same graph"})
    # same graph, different dataset (node IDs identical): should be ~identical
for g in ("linear", "d9", "d46"):
    for w in (0, 1, 2):
        res["pairs"].append(pair(("HG008T_Illumina", g, 0), ("COLO829T_ONT", g, 0), ("tile", "hull", "bins")[w]) | {"kind": "datasets, same graph+seed"})
for p in res["pairs"]:
    p["def"] = p.pop("def_")
json.dump(res, open(sys.argv[1], "w"), indent=1)
print("done")
