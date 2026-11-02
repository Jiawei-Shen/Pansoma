"""Collect the inference-only off-node mask runs: chr1, BED, after the PoN, PASS at the checkpoint's recall-0.9
validation threshold; best F1 of the chr1 curve; chr1 per-tensor somatic AP; and, for nomask, the existing
reviewer_cpu_20261013 / RIS report as a reproduction check. Writes results.md / results.json next to this script."""
import glob, json
from pathlib import Path
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
HERE = Path(__file__).resolve().parent
tasks = [l.split("\t") for l in (HERE / "tasks.tsv").read_text().split("\n") if l.strip()]
def rep(d):
    f = d / "vcf_chr1/eval_bed/report.json"
    if not f.exists():
        return None
    r = json.load(open(f)); t = r["pon"]["at_threshold"]
    out = dict(P=t["precision"], R=t["recall"], F1=t["f1"], TP=t["tp_call"], FP=t["fp"], best=r["pon"]["best"]["f1"],
               raw_F1=r["raw"]["at_threshold"]["f1"], threshold=r["threshold"])
    g = glob.glob(str(d / "test_chr1/*.metrics.json"))
    if g:
        out["AP"] = json.load(open(g[0]))["somatic_ap"]
    return out
rows = {}
for src, s, mode in tasks:
    rows.setdefault((src, s), {})[mode.strip()] = rep(R / "offnode_mask_20261030" / f"{src}__{mode.strip()}")
for (src, s), v in rows.items():
    for base in ("reviewer_cpu_20261013", "ris_7fc7396_20261012"):
        e = rep(R / base / src)
        if e:
            v["existing"] = e | {"src": base}
            break
f = lambda x, k: "-" if not x or x.get(k) is None else f"{x[k]:.3f}"
L = ["| set | model | nomask P / R / F1 (best; AP) | mask P / R / F1 (best; AP) | ΔF1 PASS | Δbest | existing F1 (repro check) |",
     "|---|---|---|---|---|---|---|"]
for (src, s), v in rows.items():
    a, b, e = v.get("nomask"), v.get("mask"), v.get("existing")
    d = "-" if not (a and b) else f"{b['F1'] - a['F1']:+.3f}"
    db = "-" if not (a and b) else f"{b['best'] - a['best']:+.3f}"
    L.append(f"| {s} | `{src}` | {f(a,'P')} / {f(a,'R')} / {f(a,'F1')} ({f(a,'best')}; {f(a,'AP')}) | "
             f"{f(b,'P')} / {f(b,'R')} / {f(b,'F1')} ({f(b,'best')}; {f(b,'AP')}) | {d} | {db} | {f(e,'F1')} |")
(HERE / "results.md").write_text("\n".join(L) + "\n")
(HERE / "results.json").write_text(json.dumps({f"{k[0]}|{k[1]}": v for k, v in rows.items()}, indent=1))
print("\n".join(L))
