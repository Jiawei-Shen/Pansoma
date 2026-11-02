import json, glob, statistics as st
from pathlib import Path
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012"); CPU = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
SETS = [("HG008", "Illumina", "HG008T Illumina"), ("HG008", "PacBio", "HG008T PacBio"), ("HG008", "ONT", "HG008T ONT"),
        ("COLO829", "Illumina", "COLO829T Illumina"), ("COLO829", "fiberseq", "COLO829T fiberseq"), ("COLO829", "ONT", "COLO829T ONT")]
def rep(run):
    for b in (RIS, CPU):
        f = b / run / "vcf_chr1/eval_bed/report.json"
        if f.exists():
            r = json.load(open(f)); t = r["pon"]["at_threshold"]
            return dict(P=t["precision"], R=t["recall"], F=t["f1"], best=r["pon"]["best"]["f1"], src=str(f))
    return None
def ap(run):
    for b in (RIS, CPU):
        g = glob.glob(str(b / run / "test_chr1/*.metrics.json"))
        if g:
            return json.load(open(g[0]))["somatic_ap"]
    return None
def fmt(x):
    return "-" if x is None else f"{x:.3f}"
rows = []
print("| dataset | graph | setting | with ch6: run(s) | with ch6 P / R / F1 after PoN (best F1) | with ch6 somatic AP | no ch6 run | no ch6 P / R / F1 (best F1) | no ch6 AP | ΔF1 PASS |")
print("|---|---|---|---|---|---|---|---|---|---|")
for s, p, name in SETS:
    # d9, norm train (no-ch6 on the control's GPU count)
    if (s, p) == ("HG008", "Illumina"):
        w, n = ["HG008_Illumina_SNV_7fc7396_ddp4_b1024"], "HG008_Illumina_SNV_nopc_7fc7396_ddp4_b1024"
    elif s == "HG008":
        w, n = [f"HG008_{p}_SNV_7fc7396"], f"HG008_{p}_SNV_nopc_7fc7396"
    else:
        w, n = [f"COLO829_{p}_SNV_h90"], f"COLO829_{p}_SNV_nopc_7fc7396"
    groups = [("v1.1 d9", "norm train (pre-norm code); 1 GPU except HG008T Illumina 4 GPU", w, n),
              ("v2.1 d46", "norm train, 4 GPU", [f"{s}_{p}_SNV_d46_7fc7396", f"{s}_{p}_SNV_d46_7fc7396_s1", f"{s}_{p}_SNV_d46_7fc7396_s2"], f"{s}_{p}_SNV_nopc_d46_7fc7396"),
              ("v2.1 d46", "normbq (BQ per dataset), 4 GPU", [f"{s}_{p}_SNV_d46_normbq", f"{s}_{p}_SNV_d46_normbq_s1", f"{s}_{p}_SNV_d46_normbq_s2"], f"{s}_{p}_SNV_nopc_d46_normbq")]
    for g, setting, wr, nr in groups:
        W = [rep(r) for r in wr]; A = [ap(r) for r in wr]
        N = rep(nr); An = ap(nr)
        if all(x is None for x in W):
            wtxt = "(not local)"
        else:
            wtxt = "; ".join(f"{fmt(x['P'])}/{fmt(x['R'])}/{fmt(x['F'])} ({fmt(x['best'])})" if x else "-" for x in W)
        atxt = " / ".join(fmt(a) for a in A)
        ntxt = "(not run/synced)" if N is None else f"{fmt(N['P'])}/{fmt(N['R'])}/{fmt(N['F'])} ({fmt(N['best'])})"
        dF = "-" if N is None or not any(W) else f"{N['F'] - st.mean([x['F'] for x in W if x]):+.3f}"
        print(f"| {name} | {g} | {setting} | {', '.join('`'+r+'`' for r in wr)} | {wtxt} | {atxt} | `{nr}` | {ntxt} | {fmt(An)} | {dF} |")
