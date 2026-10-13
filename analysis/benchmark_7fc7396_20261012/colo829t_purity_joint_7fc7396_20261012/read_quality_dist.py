"""How often is a whole read of low quality? Per set, the first 3,000 tensors of chr2 shard 0 (SNV): for every row (read)
with >= 10 covered BQ cells in the window, the read's median BQ; quantiles of that median, the share of reads far below
the set's typical read (median more than 5 / 10 Q under the set's median read median) and below fixed BQ values, and
the same for the A1 (ALT-carrying) rows of somatic and of germline tensors. Shows how much a per-read BQ normalization
could hide. Run in the container: ris/submit.sh -g 0 -n 2 -m 16 NAME py.sh $PWD/read_quality_dist.py"""
import json
import re
from pathlib import Path

import numpy as np

D = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors")
HERE = Path(__file__).resolve().parent
RG = re.compile(r'"row_groups": (\[[^\]]*\])')
N = 3000
SETS = ["HG008T_Illumina", "HG008T_PacBio", "HG008T_ONT", "COLO829T_Illumina", "COLO829T_Illumina_P25",
        "COLO829T_fiberseq", "COLO829T_ONT"]
res = {}
for s in SETS:
    t = D / s / "tensors" / "SNV"
    x = np.asarray(np.load(t / "chr2_shard_00000_data.npy", mmap_mode="r")[:N, [1, 4]])  # BQ, operation
    lab = np.load(t / "chr2_shard_00000_labels.npy")[:N]
    a1_end = []
    with open(t / "chr2_variant_summary.ndjson") as f:
        for i, line in enumerate(f):
            if i >= N:
                break
            g = json.loads(RG.search(line).group(1))
            a1_end.append(next((r["end_row"] for r in g if r["allele"] == "A1"), 0))
    bq = x[:, 0].astype(np.float32); valid = (x[:, 1] > 0) & (x[:, 0] >= 0)
    ncell = valid.sum(2)
    med = np.nanmedian(np.where(valid, bq, np.nan), axis=2)  # (N, rows)
    ok = ncell >= 10
    m = med[ok]
    ref = float(np.median(m))
    a1 = np.zeros_like(ok); som = np.zeros_like(ok); ger = np.zeros_like(ok)
    for i, e in enumerate(a1_end):
        a1[i, :e] = True
        som[i] = lab[i] == 1; ger[i] = lab[i] == 2
    def share(sel):
        v = med[ok & sel]
        return {"n": int(len(v)), "below_ref_minus5": round(float((v < ref - 5).mean()), 4) if len(v) else None,
                "below_ref_minus10": round(float((v < ref - 10).mean()), 4) if len(v) else None}
    res[s] = {"rows": int(ok.sum()), "read_median_quantiles": {q: float(np.quantile(m, q)) for q in (0.01, 0.05, 0.25, 0.5, 0.75)},
              "set_reference": ref, "below_ref_minus5": round(float((m < ref - 5).mean()), 4),
              "below_ref_minus10": round(float((m < ref - 10).mean()), 4), "below_20": round(float((m < 20).mean()), 4),
              "somatic_A1_rows": share(a1 & som), "germline_A1_rows": share(a1 & ger), "all_REF_or_other_rows": share(~a1)}
    r = res[s]
    print(f"{s:22s} reads {r['rows']:,} | read-median BQ q01/q05/q25/q50/q75 "
          + "/".join(f"{v:.0f}" for v in r["read_median_quantiles"].values())
          + f" | < median-5: {r['below_ref_minus5']:.4f}, < median-10: {r['below_ref_minus10']:.4f}, < 20: {r['below_20']:.4f}"
          + f" | somatic A1 rows < median-5: {r['somatic_A1_rows']['below_ref_minus5']}, germline A1 rows: {r['germline_A1_rows']['below_ref_minus5']}",
          flush=True)
(HERE / "read_quality_dist.json").write_text(json.dumps(res, indent=1) + "\n")
