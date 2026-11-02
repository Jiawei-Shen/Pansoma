"""Q1/Q2: per-run chr1 SNV P/R/F1 after PoN (BED, PASS), raw (before PoN) F1, PASS threshold, best epoch."""
import json, re, statistics, sys
from pathlib import Path
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
CPU = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
SUM = json.load(open("/scratch/jshen/Github/Pansoma/analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/summary.json"))
SETS = [("HG008", "Illumina"), ("HG008", "PacBio"), ("HG008", "ONT"), ("COLO829", "Illumina"), ("COLO829", "fiberseq"), ("COLO829", "ONT")]

def report(run):
    for base in (RIS, CPU):
        f = base / run / "vcf_chr1" / "eval_bed" / "report.json"
        if f.exists():
            return base, json.load(open(f))
    return None, None

def epoch(run):
    for base in (RIS, CPU):
        f = base / run / "train.log"
        if f.exists():
            m = re.findall(r"done; best epoch (\d+) \(f1 ([0-9.]+), somatic threshold ([0-9.e-]+)\)", f.read_text())
            if m:
                return int(m[-1][0]), float(m[-1][1]), float(m[-1][2])
    return None, None, None

def row(run):
    base, r = report(run)
    ep, vf1, bt = epoch(run)
    out = {"run": run, "src": None if base is None else base.name, "epoch": ep, "val_f1": vf1}
    if r:
        p, w = r["pon"]["at_threshold"], r["raw"]["at_threshold"]
        out.update(threshold=r["threshold"], P=p["precision"], R=p["recall"], F1=p["f1"], tp=p["tp_call"], fp=p["fp"],
                   rawP=w["precision"], rawR=w["recall"], rawF1=w["f1"], bestF1=r["pon"]["best"]["f1"], unplaced=r.get("unplaced_pass"))
    return out

def runs_train(s, p, arm):
    if arm == "linear":
        b = f"{s}_{p}_SNV_linear_7fc7396"
    elif arm == "d9":
        b = f"{s}_{p}_SNV_7fc7396_ddp4_b1024"
    else:
        b = f"{s}_{p}_SNV_d46_7fc7396"
    return [b, b + "_s1", b + "_s2"]

def runs_normbq(s, p, arm, chosen=False):
    b = f"{s}_{p}_SNV_{arm}_normbq"
    if chosen and (s, p, arm) == ("COLO829", "Illumina", "d9"):
        return [b + "_sweep_e06", b + "_s1_sweep_e05", b + "_s2_sweep_e06"]
    return [b, b + "_s1", b + "_s2"]

res = {"train": {}, "normbq": {}}
for s, p in SETS:
    for arm in ("linear", "d9", "d46"):
        res["train"][f"{s}|{p}|{arm}"] = [row(r) for r in runs_train(s, p, arm)]
    for arm in ("d9", "d46"):
        res["normbq"][f"{s}|{p}|{arm}"] = [row(r) for r in runs_normbq(s, p, arm)]
    if (s, p) == ("COLO829", "Illumina"):
        res["normbq"][f"{s}|{p}|d9_chr1epochs"] = [row(r) for r in runs_normbq(s, p, "d9", True)]
# summary.json values (RIS originals) for cross-check
seeds = {(x["set"], x["arm"]): x["scores"] for x in SUM["seeds"]}
res["summary_train"] = {f"{k[0]}|{k[1]}": v for k, v in seeds.items()}
res["summary_normbq"] = {f"{x['set']}|{x['arm']}": x["scores"] for x in SUM["normbq_seeds"]}
json.dump(res, open(sys.argv[1], "w"), indent=1)
