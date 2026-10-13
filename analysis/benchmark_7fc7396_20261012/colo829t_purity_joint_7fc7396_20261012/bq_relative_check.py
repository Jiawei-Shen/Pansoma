"""Does a per-read BQ normalization remove the read-source cue of the in-silico mixtures? (bq_alt_ref.py's test, three
encodings.) Per COLO829T SNV set, chr2 shard 0 (training data), tensors labelled somatic (1) or germline (2) with >= 2 A1
and >= 2 REF rows: the mean BQ at the anchor column of the A1 rows minus that of the REF rows, with BQ
- raw: as stored (what the model's global z-score keeps, up to scale);
- read_centered: BQ minus the median BQ of that read's covered cells in the 101-column window;
- read_rank: the percentile rank of the cell's BQ among that read's covered cells (0-1).
The AUC of the difference alone, somatic vs germline (lower A1 value = somatic), should stay at its P100 value at every
purity if the encoding hides the read source. Run in the container: ris/submit.sh -g 0 -n 2 -m 16 NAME py.sh $PWD/bq_relative_check.py"""
import json
import re
from pathlib import Path

import numpy as np

D = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors")
HERE = Path(__file__).resolve().parent
RG = re.compile(r'"row_groups": (\[[^\]]*\])')
ANCHOR = re.compile(r'"anchor_column": (\d+)')
SHARD = re.compile(r'"shard_file": "([^"]+)"')
MIN_CELLS = 10  # rows with fewer covered BQ cells in the window are left out of the per-read statistics


def auc(pos, neg):
    x = np.concatenate([pos, neg])
    _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    r = np.argsort(np.argsort(x, kind="mergesort"), kind="mergesort") + 1.0
    s = np.zeros(len(cnt)); np.add.at(s, inv, r); r = (s / cnt)[inv]
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def encodings(t):
    """(H, W) BQ and coverage of one tensor -> {name: (H, W) float with NaN outside valid cells}."""
    bq = t[1].astype(np.float32); valid = (t[4] > 0) & (t[1] >= 0)
    raw = np.where(valid, bq, np.nan)
    n = valid.sum(1)
    med = np.nanmedian(np.where(valid, bq, np.nan), axis=1, keepdims=True) if valid.any() else np.zeros((bq.shape[0], 1))
    centered = np.where(valid & (n[:, None] >= MIN_CELLS), bq - med, np.nan)
    # percentile rank within the row (average rank of ties), over valid cells
    rank = np.full(bq.shape, np.nan, np.float32)
    for i in np.nonzero(n >= MIN_CELLS)[0]:
        v = bq[i, valid[i]]; order = v.argsort(kind="mergesort"); rk = np.empty(len(v)); rk[order] = np.arange(len(v))
        _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True); s = np.zeros(len(cnt)); np.add.at(s, inv, rk)
        rank[i, valid[i]] = (s / cnt)[inv] / max(len(v) - 1, 1)
    return {"raw": raw, "read_centered": centered, "read_rank": rank}


res = {}
for plat in ("Illumina", "fiberseq", "ONT"):
    for q in ("", "_P75", "_P50", "_P25"):
        s = f"COLO829T_{plat}{q}"; t = D / s / "tensors" / "SNV"
        x = np.load(t / "chr2_shard_00000_data.npy", mmap_mode="r"); lab = np.load(t / "chr2_shard_00000_labels.npy")
        diffs = {e: {1: [], 2: []} for e in ("raw", "read_centered", "read_rank")}
        with open(t / "chr2_variant_summary.ndjson") as f:
            for i, line in enumerate(f):
                if SHARD.search(line).group(1) != "chr2_shard_00000_data.npy":
                    break
                if lab[i] not in (1, 2):
                    continue
                groups = json.loads(RG.search(line).group(1)); col = int(ANCHOR.search(line).group(1))
                a1 = [(g["start_row"], g["end_row"]) for g in groups if g["allele"] == "A1"]
                ref = [(g["start_row"], g["end_row"]) for g in groups if g["allele"] == "REF"]
                if not a1 or not ref:
                    continue
                enc = encodings(np.asarray(x[i]))
                for e, m in enc.items():
                    va = m[a1[0][0]:a1[0][1], col]; vr = m[ref[0][0]:ref[0][1], col]
                    va, vr = va[~np.isnan(va)], vr[~np.isnan(vr)]
                    if len(va) >= 2 and len(vr) >= 2:
                        diffs[e][int(lab[i])].append(va.mean() - vr.mean())
        res[s] = {}
        for e, d in diffs.items():
            so, ge = np.array(d[1]), np.array(d[2])
            res[s][e] = {"somatic_n": len(so), "somatic_median": float(np.median(so)), "germline_median": float(np.median(ge)),
                         "auc": auc(-so, -ge)}
        print(f"{s:24s} " + " | ".join(f"{e}: AUC {r['auc']:.3f} (som {r['somatic_median']:+.3f}, germ {r['germline_median']:+.3f})"
                                       for e, r in res[s].items()), flush=True)
(HERE / "bq_relative_check.json").write_text(json.dumps(res, indent=1) + "\n")
