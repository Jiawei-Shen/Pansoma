"""E10b add-on (R1.3, 2026-10-20): COLO829T Illumina SNV, the HPRC v1.1 (d9) arm at epoch 6 (the user's decision of
2026-10-20, RIS summary `analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/README.md`, b662f00) instead of
the validation-selected epoch 8 used by E10 / E16 / giab_strata.py, per SMaHT stratum and GIAB stratum, against the
GRCh38-only (linear) graph; paired block bootstrap as E16 (1-Mb blocks of the BED, 2,000 resamples, seed 20261015).

Call sets (chr1, PASS after the PoN, rtg eval_bed pon_calls; positions = VCF POS of tp-baseline / fn / tp / fp):
  linear  reviewer_cpu_20261013/COLO829_Illumina_SNV_linear_7fc7396/vcf_chr1
  d46     reviewer_cpu_20261013/COLO829_Illumina_SNV_d46_7fc7396/vcf_chr1
  d9_ep8  reviewer_cpu_20261013/COLO829_Illumina_SNV_7fc7396_ddp4_b1024/vcf_chr1   (the d9 of E10 / E16, epoch 8)
  d9_ep6  ris_7fc7396_20261012/COLO829_Illumina_SNV_d9_normbq_sweep_e06/vcf_chr1   (seed 0 at epoch 6, --norm dataset)
Writes colo_ill_d9_ep6.md next to this script.
usage: /wanglab/jshen/anaconda3/bin/python colo_ill_d9_ep6.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "E16_bootstrap"))
sys.path.insert(0, str(HERE))
import callsets as C  # noqa: E402
import e16_bootstrap as E  # noqa: E402
import giab_strata as G  # noqa: E402

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
RUNS = {"linear": R / "reviewer_cpu_20261013/COLO829_Illumina_SNV_linear_7fc7396/vcf_chr1",
        "d46": R / "reviewer_cpu_20261013/COLO829_Illumina_SNV_d46_7fc7396/vcf_chr1",
        "d9_ep8": R / "reviewer_cpu_20261013/COLO829_Illumina_SNV_7fc7396_ddp4_b1024/vcf_chr1",
        "d9_ep6": R / "ris_7fc7396_20261012/COLO829_Illumina_SNV_d9_normbq_sweep_e06/vcf_chr1"}
SMAHT = {s: Path(f"/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_{s}_v2.bed") for s in ("easy", "difficult", "extreme")}


def bed_intervals(path):
    iv = sorted((s, e) for s, e in E.read_bed_chr1(str(path)))
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    a = np.array(merged, dtype=np.int64).reshape(-1, 2)
    return a[:, 0], a[:, 1]


def main():
    units = {}
    for arm, vdir in RUNS.items():
        rep = json.loads((vdir / "eval_bed" / "report.json").read_text())
        u = C.Unit(uid=f"{arm}", run=vdir.parent.name, set="COLO829T_Illumina", kind="SNV", role="insample",
                   variant="standard", model=vdir.parent.name)
        u.rtg["bed"] = vdir / "eval_bed" / "vcfeval" / "pon_calls"
        u.truth, u.bed = C.resolved(rep["truth"]), C.resolved(rep["bed"])
        units[arm] = (u, rep["pon"]["at_threshold"])
    beds = {u.bed for u, _ in units.values()}
    assert len(beds) == 1, beds
    strata = [(f"SMaHT {s}", bed_intervals(p)) for s, p in SMAHT.items()]
    strata += [(lab, G.chr1_intervals(G.GIAB / rel)) for _, rel, lab in G.STRATA]
    strata.insert(0, ("all (BED)", None))
    bs = E.BlockSets(C.fasta_path(), 2000, 20261015)
    blocks = bs.get("bed", next(iter(beds)), E.WIDTHS[0])
    pos = {arm: {k: np.array(v, dtype=np.int64) for k, v in E.positions(u, "bed").items()} for arm, (u, _) in units.items()}
    res, tot = {}, {}
    for arm in units:
        for name, iv in strata:
            pp = {k: (v if iv is None else v[G.inside(v, iv)]).tolist() for k, v in pos[arm].items()}
            m = E.block_counts(pp, blocks)
            tot[(arm, name)] = m.sum(axis=0)
            res[(arm, name)] = blocks.mult @ m
    sg = lambda x: "-" if x is None or np.isnan(x) else f"{x:+.3f}"  # noqa: E731
    f3 = lambda x: "-" if x is None or np.isnan(x) else f"{x:.3f}"  # noqa: E731
    L = ["# COLO829T Illumina SNV: HPRC v1.1 (d9) at epoch 6 against epoch 8, per stratum (chr1, BED, after the PoN)", "",
         "Whole-BED check against each run's report.json (P / R / F1): " + "; ".join(
             f"{a} {r['precision']:.3f} / {r['recall']:.3f} / {r['f1']:.3f}" for a, (_, r) in units.items()), "",
         "| stratum | arm | PASS calls | TP / FP / FN | P | R | F1 | dF1 vs linear [95% CI] |", "|---|---|---:|---|---|---|---|---|"]
    out = []
    for name, _ in strata:
        for arm in ("linear", "d9_ep8", "d9_ep6", "d46"):
            t = tot[(arm, name)]
            p, r, f = E.metrics(t)
            d = lo = hi = np.nan
            if arm != "linear":
                diff = E.metrics(res[(arm, name)])[2] - E.metrics(res[("linear", name)])[2]
                d = float(f - E.metrics(tot[("linear", name)])[2])
                lo, hi, _, _ = E.ci(diff)
            out.append(dict(stratum=name, arm=arm, tp_baseline=int(t[0]), fn=int(t[1]), tp=int(t[2]), fp=int(t[3]),
                            P=float(p), R=float(r), F1=float(f), dF1=d, dF1_lo=lo, dF1_hi=hi))
            L.append(f"| {name} | {arm} | {int(t[2] + t[3]):,} | {int(t[0])} / {int(t[3])} / {int(t[1])} | {f3(float(p))} | "
                     f"{f3(float(r))} | {f3(float(f))} | {sg(d)} [{sg(lo)}, {sg(hi)}] |")
    (HERE / "colo_ill_d9_ep6.md").write_text("\n".join(L) + "\n")
    (HERE / "colo_ill_d9_ep6.json").write_text(json.dumps(out, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
