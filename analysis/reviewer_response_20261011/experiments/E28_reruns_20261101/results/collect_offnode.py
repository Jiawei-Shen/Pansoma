"""E28 (1): v2.1 SNV seed 0 retrained with and without the off-node mask (same 1-H100 setup). chr1, BED, after the PoN,
PASS at the checkpoint's validation threshold; best F1 of the chr1 curve; chr1 per-tensor somatic AP; validation-selected
epoch and validation F1. The RIS original (4 GPUs) is shown as a reference. Writes offnode.md / offnode.json here."""
import glob, json, re
from pathlib import Path
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
HERE = Path(__file__).resolve().parent
SETS = [("HG008", "Illumina"), ("HG008", "PacBio"), ("HG008", "ONT"), ("COLO829", "Illumina"), ("COLO829", "fiberseq"), ("COLO829", "ONT")]
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
    lg = d / "train.log"
    if lg.exists():
        m = re.findall(r"done; best epoch (\d+) \(f1 ([0-9.]+)", lg.read_text())
        if m:
            out["epoch"], out["val_f1"] = int(m[-1][0]), float(m[-1][1])
    return out
rows = {}
for s, p in SETS:
    rows[f"{s}_{p}"] = dict(mask=rep(R / "offnode_retrain_20261101" / f"{s}_{p}_SNV_d46_offnode_mask"),
                            nomask=rep(R / "offnode_retrain_20261101" / f"{s}_{p}_SNV_d46_offnode_nomask"),
                            ris=next((x for x in (rep(R / b / f"{s}_{p}_SNV_d46_7fc7396") for b in ("ris_7fc7396_20261012", "reviewer_cpu_20261013")) if x), None))
f = lambda x, k: "-" if not x or x.get(k) is None else (str(x[k]) if isinstance(x[k], int) else f"{x[k]:.3f}")
L = ["| chr1 SNV, v2.1, seed 0 | no mask P / R / F1 (best; AP; epoch) | mask P / R / F1 (best; AP; epoch) | ΔP | ΔR | ΔF1 PASS | Δbest | ΔAP | original (4 GPUs) F1; epoch |",
     "|---|---|---|---|---|---|---|---|---|"]
for k, v in rows.items():
    a, b, e = v["nomask"], v["mask"], v["ris"]
    dd = lambda q: "-" if not (a and b) else f"{b[q] - a[q]:+.3f}"
    L.append(f"| {k} | {f(a,'P')} / {f(a,'R')} / {f(a,'F1')} ({f(a,'best')}; {f(a,'AP')}; {f(a,'epoch')}) | "
             f"{f(b,'P')} / {f(b,'R')} / {f(b,'F1')} ({f(b,'best')}; {f(b,'AP')}; {f(b,'epoch')}) | "
             f"{dd('P')} | {dd('R')} | {dd('F1')} | {dd('best')} | {dd('AP')} | {f(e,'F1')}; {f(e,'epoch')} |")
(HERE / "offnode.md").write_text("\n".join(L) + "\n")
(HERE / "offnode.json").write_text(json.dumps(rows, indent=1))
print("\n".join(L))
