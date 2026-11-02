import json, statistics as st
d = json.load(open("/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/main_runs.json"))
SETS = [("HG008", "Illumina", "HG008T Illumina"), ("HG008", "PacBio", "HG008T PacBio"), ("HG008", "ONT", "HG008T ONT"),
        ("COLO829", "Illumina", "COLO829T Illumina"), ("COLO829", "fiberseq", "COLO829T fiberseq"), ("COLO829", "ONT", "COLO829T ONT")]
ARMN = {"linear": "GRCh38-only", "d9": "HPRC v1.1 (d9)", "d46": "HPRC v2.1 (d46)", "d9_chr1epochs": "HPRC v1.1 (d9), chr1-chosen epochs 6/5/6"}
def ms(x):
    return f"{st.mean(x):.3f} ± {st.stdev(x):.3f}"
def cell(x):
    return " / ".join(f"{v:.3f}" for v in x) + f" ({ms(x)})"
out = []
for setting, sumkey in (("train", "summary_train"), ("normbq", "summary_normbq")):
    out.append(f"\n### setting {setting}\n")
    out.append("| dataset | reference | P s0/s1/s2 (mean ± SD) | R | F1 after PoN | F1 before PoN (raw, PASS) | PASS threshold s0/s1/s2 | best epoch s0/s1/s2 |")
    out.append("|---|---|---|---|---|---|---|---|")
    means = {}
    for s, p, name in SETS:
        arms = ["linear", "d9", "d46"] if setting == "train" else ["d9", "d46"] + (["d9_chr1epochs"] if (s, p) == ("COLO829", "Illumina") else [])
        for arm in arms:
            rows = d[setting][f"{s}|{p}|{arm}"]
            sk = f"{name}|{arm}"
            if setting == "normbq" and (s, p) == ("COLO829", "Illumina") and arm in ("d9", "d9_chr1epochs"):
                sk = f"{name}|d9 validation-selected (epoch 8)" if arm == "d9" else f"{name}|d9"
            summ = d[sumkey].get(sk)
            if summ and all(summ):
                P = [x["p"] for x in summ]; R = [x["r"] for x in summ]; F = [x["pass"] for x in summ]; src = "RIS summary.json"
            else:
                P = [x["P"] for x in rows]; R = [x["R"] for x in rows]; F = [x["F1"] for x in rows]; src = "local report.json"
            raw = [x["rawF1"] for x in rows]
            th = [x["threshold"] for x in rows]
            ep = [str(x["epoch"]) for x in rows]
            means[(name, arm)] = dict(P=P, R=R, F=F, raw=raw)
            out.append(f"| {name} | {ARMN[arm]} | {cell(P)} | {cell(R)} | {cell(F)} | {cell(raw)} | {' / '.join(f'{t:.3f}' for t in th)} | {' / '.join(ep)} |")
    # derived
    out.append("")
    if setting == "train":
        dl = []
        for s, p, name in SETS:
            for arm in ("d9", "d46"):
                a, b = means[(name, arm)], means[(name, "linear")]
                dF = st.mean(a["F"]) - st.mean(b["F"]); dP = st.mean(a["P"]) - st.mean(b["P"]); dR = st.mean(a["R"]) - st.mean(b["R"])
                draw = st.mean(a["raw"]) - st.mean(b["raw"])
                wins = sum(x > y for x, y in zip(a["F"], b["F"]))
                dl.append((name, arm, dF, dP, dR, draw, wins, [round(x - y, 3) for x, y in zip(a["F"], b["F"])]))
        out.append("| dataset | HPRC arm | ΔF1 after PoN (mean) | ΔP | ΔR | ΔF1 before PoN | seed-matched wins | per-seed ΔF1 |")
        out.append("|---|---|---|---|---|---|---|---|")
        for r in dl:
            out.append(f"| {r[0]} | {ARMN[r[1]]} | {r[2]:+.3f} | {r[3]:+.3f} | {r[4]:+.3f} | {r[5]:+.3f} | {r[6]}/3 | {r[7]} |")
        out.append(f"\nHPRC − GRCh38-only ΔF1 after PoN: {min(r[2] for r in dl):+.3f} to {max(r[2] for r in dl):+.3f}; "
                   f"ΔP {min(r[3] for r in dl):+.3f} to {max(r[3] for r in dl):+.3f}; ΔR {min(r[4] for r in dl):+.3f} to {max(r[4] for r in dl):+.3f}; "
                   f"seed-matched wins {sum(r[6] for r in dl)}/36")
    out.append("\n| dataset | v2.1 − v1.1 ΔF1 after PoN | ΔP after | ΔR after | ΔF1 before PoN | ΔP before | ΔR before |")
    out.append("|---|---|---|---|---|---|---|")
    vv = []
    for s, p, name in SETS:
        d9keys = ["d9"] + (["d9_chr1epochs"] if setting == "normbq" and (s, p) == ("COLO829", "Illumina") else [])
        for k in d9keys:
            a, b = means[(name, "d46")], means[(name, k)]
            ra = [x for x in d[setting][f"{s}|{p}|d46"]]; rb = [x for x in d[setting][f"{s}|{p}|{k}"]]
            rP = st.mean([x["rawP"] for x in ra]) - st.mean([x["rawP"] for x in rb])
            rR = st.mean([x["rawR"] for x in ra]) - st.mean([x["rawR"] for x in rb])
            row = (name + ("" if k == "d9" else " (v1.1 at chr1-chosen epochs 6/5/6)"), st.mean(a["F"]) - st.mean(b["F"]), st.mean(a["P"]) - st.mean(b["P"]),
                   st.mean(a["R"]) - st.mean(b["R"]), st.mean(a["raw"]) - st.mean(b["raw"]), rP, rR)
            vv.append((k, row))
            out.append(f"| {row[0]} | {row[1]:+.3f} | {row[2]:+.3f} | {row[3]:+.3f} | {row[4]:+.3f} | {row[5]:+.3f} | {row[6]:+.3f} |")
    v = [r for k, r in vv if k == "d9"]
    out.append(f"\nv2.1 − v1.1 (validation-selected): after PoN {min(r[1] for r in v):+.3f} to {max(r[1] for r in v):+.3f}; before PoN {min(r[4] for r in v):+.3f} to {max(r[4] for r in v):+.3f}; before-PoN ΔP {min(r[5] for r in v):+.3f} to {max(r[5] for r in v):+.3f}")
    v5 = [r for k, r in vv if k == "d9" and not r[0].startswith("COLO829T Illumina")]
    out.append(f"without COLO829T Illumina: after {min(r[1] for r in v5):+.3f} to {max(r[1] for r in v5):+.3f}; before {min(r[4] for r in v5):+.3f} to {max(r[4] for r in v5):+.3f}")
print("\n".join(out))
