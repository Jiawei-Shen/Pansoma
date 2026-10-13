"""Does the base quality of the ALT reads tell somatic from germline in the in-silico dilutions? Per COLO829T SNV set,
on chr2 shard 0 (training data): for each tensor labelled somatic (1) or germline (2) with >= 2 A1 rows and >= 2 REF
rows, the mean BQ at the anchor column of the A1 (representative-allele) rows minus that of the REF rows (covered
cells, BQ >= 0). In a dilution, somatic ALT reads come from the tumor read set only, germline ALT reads from both, so a
tumor / normal BQ offset shows as a somatic-specific shift; the AUC is how well that difference alone ranks somatic
above germline. Run in the container: ris/submit.sh -g 0 -n 2 -m 16 NAME py.sh $PWD/bq_alt_ref.py"""
import json
import re
from pathlib import Path

import numpy as np

D = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors")
HERE = Path(__file__).resolve().parent
RG = re.compile(r'"row_groups": (\[[^\]]*\])')
ANCHOR = re.compile(r'"anchor_column": (\d+)')
SHARD = re.compile(r'"shard_file": "([^"]+)"')


def auc(pos, neg):
    """P(a random positive > a random negative), ties half."""
    x = np.concatenate([pos, neg]); r = np.argsort(np.argsort(x, kind="mergesort"), kind="mergesort") + 1.0
    _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)  # average ranks of ties
    s = np.zeros(len(cnt)); np.add.at(s, inv, r); r = (s / cnt)[inv]
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


res = {}
for plat in ("Illumina", "fiberseq"):
    for q in ("", "_P75", "_P50", "_P25"):
        s = f"COLO829T_{plat}{q}"; t = D / s / "tensors" / "SNV"
        x = np.load(t / "chr2_shard_00000_data.npy", mmap_mode="r")
        lab = np.load(t / "chr2_shard_00000_labels.npy")
        diffs = {1: [], 2: []}
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
                bq = np.asarray(x[i, 1, :, col], dtype=np.int16); cov = np.asarray(x[i, 4, :, col]) > 0
                ra = np.r_[a1[0][0]:a1[0][1]]; rr = np.r_[ref[0][0]:ref[0][1]]
                va = bq[ra][cov[ra] & (bq[ra] >= 0)]; vr = bq[rr][cov[rr] & (bq[rr] >= 0)]
                if len(va) >= 2 and len(vr) >= 2:
                    diffs[int(lab[i])].append(va.mean() - vr.mean())
        so, ge = np.array(diffs[1]), np.array(diffs[2])
        res[s] = {"somatic": {"n": len(so), "median": float(np.median(so)), "mean": float(so.mean())},
                  "germline": {"n": len(ge), "median": float(np.median(ge)), "mean": float(ge.mean())},
                  "auc_lower_alt_bq_is_somatic": auc(-so, -ge)}
        r = res[s]
        print(f"{s:24s} A1 - REF mean BQ: somatic n {r['somatic']['n']:5d} median {r['somatic']['median']:+.2f} | germline n "
              f"{r['germline']['n']:5d} median {r['germline']['median']:+.2f} | AUC (lower A1 BQ = somatic) "
              f"{r['auc_lower_alt_bq_is_somatic']:.3f}", flush=True)
(HERE / "bq_alt_ref.json").write_text(json.dumps(res, indent=1) + "\n")
