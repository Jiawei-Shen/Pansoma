"""Tables of the HPRC v2.1 d46 runs of 2026-10-12 (repo 7fc7396, the v2 README recipes) against the HPRC v1.1 d9 runs
of the same sample, platform and kind. Standard library only (runs on the login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py

Reads, per run: train.log and metrics.jsonl (training, best epoch), <run>/test_chr1/*.metrics.json (chr1 per tensor and
per truth allele) and <run>/<vcf dir>/eval_{bed,nobed}/report.json (vcf.sh / vcfeval_only.sh). Writes tables.md,
summary.json and runs/<run>/ (train.log, <vcf dir>.eval_bed.txt, <vcf dir>.eval_nobed.txt) next to this script.
"""
import glob
import json
import os
import re
import shutil
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
HERE = Path(__file__).resolve().parent
SAMPLES = {"HG008T": ("HG008", ("Illumina", "PacBio", "ONT")), "COLO829T": ("COLO829", ("Illumina", "fiberseq", "ONT"))}
KINDS = ("SNV", "INDEL")


def d9_run(sample, platform, kind):
    """The HPRC v1.1 d9 run compared with: (run, vcf dir, note)."""
    s = SAMPLES[sample][0]
    if sample == "HG008T":  # 2026-10-12, the same recipes (Illumina SNV on 4 GPUs)
        if (platform, kind) == ("Illumina", "SNV"):
            return f"{s}_Illumina_SNV_7fc7396_ddp4_b1024", "vcf_chr1", "7fc7396, 4 GPUs"
        return f"{s}_{platform}_{kind}_7fc7396", "vcf_chr1", "7fc7396, 1 GPU"
    if kind == "SNV":  # 2026-10-03, the same SNV recipe on one GPU; the SNV PoN rule is unchanged since
        return f"{s}_{platform}_SNV_h90", "vcf_chr1", "9d6bd8d, 1 GPU"
    # 2026-10-06, the INDEL recipe at lr 4e-4 (45d56ad), re-scored with the 7fc7396 PoN
    return f"{s}_{platform}_INDEL_small", "vcf_chr1_pon7fc7396", "45d56ad lr 4e-4, 1 GPU"


def d46_run(sample, platform, kind):
    return f"{SAMPLES[sample][0]}_{platform}_{kind}_d46_7fc7396"


def load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def training(d):
    out = {}
    try:
        log = (d / "train.log").read_text(errors="replace")
    except OSError:
        return out
    m = re.search(r"training tensors ([\d,]+) \((\{[^)]*\})\), validation ([\d,]+) \((\{[^)]*\})\)", log)
    if m:
        out["train_tensors"] = int(m.group(1).replace(",", ""))
        out["train_counts"] = json.loads(m.group(2).replace("'", '"'))
        out["val_tensors"] = int(m.group(3).replace(",", ""))
        out["val_counts"] = json.loads(m.group(4).replace("'", '"'))
    m = re.search(r"done; best epoch (\d+)", log)
    if m:
        out["best_epoch"] = int(m.group(1))
    m = re.search(r"device cuda(?::\d+)? x (\d+)", log)
    if m:
        out["gpus"] = int(m.group(1))
    epochs = []
    try:
        with open(d / "metrics.jsonl") as f:
            epochs = [json.loads(line) for line in f if line.strip()]
    except OSError:
        pass
    out["epochs"] = len(epochs)
    if epochs:
        rest = epochs[1:] or epochs
        out["train_minutes"] = round(sum(e["train_seconds"] for e in epochs) / 60, 1)
        out["tensors_per_second"] = round(sum(e["tensors_per_second"] for e in rest) / len(rest))
    best = next((e for e in epochs if e["epoch"] == out.get("best_epoch")), None)
    if best:
        v = best["val"]
        out["val"] = {"somatic_ap": v["somatic_ap"], "somatic_f1": v["thresholded"]["somatic"]["f1"],
                      "threshold_r09": (v["somatic_at_recall"].get("0.9") or {}).get("threshold")}
    return out


def test(d):
    files = glob.glob(str(d / "test_chr1" / "*.metrics.json"))
    m = load(files[0]) if files else None
    if not m:
        return {}
    t = m["truth"]
    return {"tensors": m["tensors"], "somatic_tensors": m["thresholded"]["somatic"]["support"],
            "somatic_ap": m["somatic_ap"], "truth_alleles": t["truth_alleles"], "with_tensor": t["with_tensor"],
            "ceiling": t["ceiling"], "truth_best": t["best"]}


def rtg(vcf_dir, where):
    r = load(vcf_dir / where / "report.json")
    if not r:
        return None
    out = {"truth_records": r["truth_records"], "threshold": r.get("threshold"), "unplaced_pass": r.get("unplaced_pass"),
           "truth_pon_tagged": (r.get("truth_pon_tags") or {}).get("tagged")}
    for s in ("raw", "pon"):
        if s in r:
            a, b = r[s]["at_threshold"], r[s]["best"]
            out[s] = {"baseline": r[s]["baseline"], "evaluated": a["tp_call"] + a["fp"], "tp": a["tp_baseline"],
                      "fp": a["fp"], "fn": a["fn"], "p": a["precision"], "r": a["recall"], "f1": a["f1"],
                      "ceiling": r[s]["ceiling"],
                      "best": {k: b[k] for k in ("f1", "p_somatic", "precision", "recall", "tp_baseline", "fp", "fn")}}
    return out


def one(run, vcf_dir):
    d = RUNS / run
    out = {"run": run, "vcf_dir": vcf_dir, "training": training(d), "test": test(d),
           "eval_bed": rtg(d / vcf_dir, "eval_bed"), "eval_nobed": rtg(d / vcf_dir, "eval_nobed")}
    o = HERE / "runs" / run
    for src, dst in [(d / "train.log", "train.log")] + [(d / vcf_dir / w / "report.txt", f"{vcf_dir}.{w}.txt")
                                                         for w in ("eval_bed", "eval_nobed")]:
        if src.is_file():
            o.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, o / dst)
    return out


def f3(x):
    return "-" if x is None else f"{x:.3f}"


def pon(x, where="eval_bed"):
    e = x[where]
    return e["pon"] if e and e.get("pon") else None


def tables(rows):
    L = ["## chr1 PASS F1 after the PoN, in the BED: HPRC v2.1 d46 against HPRC v1.1 d9 (the chr1 curve's best F1 in "
         "brackets)\n",
         "| sample | platform | kind | d9 | d46 | Δ PASS F1 | d9 run |", "|---|---|---|---:|---:|---:|---|"]
    for r in rows:
        a, b = pon(r["d9"]), pon(r["d46"])
        L.append(f"| {r['sample']} | {r['platform']} | {r['kind']} | "
                 f"{f3(a and a['f1'])} ({f3(a and a['best']['f1'])}) | **{f3(b and b['f1'])}** ({f3(b and b['best']['f1'])}) | "
                 f"{'-' if not (a and b) else f'{b['f1'] - a['f1']:+.3f}'} | `{r['d9']['run']}` ({r['d9_note']}) |")
    L.append("\n## chr1 details, both graphs (rtg vcfeval in the BED; PASS = validation recall-0.9 threshold)\n")
    L.append("| sample | platform | kind | graph | test tensors (somatic) | truth with a tensor (ceiling) | "
             "PASS evaluated (unplaced) | raw F1 | PoN TP / FP | PoN P / R / F1 | PoN best F1 at p | PoN ceiling | "
             "no BED: PoN F1 |")
    L.append("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        for g in ("d9", "d46"):
            x = r[g]; t, e = x["test"], x["eval_bed"]; p, n = pon(x), pon(x, "eval_nobed")
            tt = f"{t['tensors']:,} ({t['somatic_tensors']:,})" if t else "-"
            tc = f"{t['with_tensor']} / {t['truth_alleles']} ({f3(t['ceiling'])})" if t else "-"
            if not p:
                L.append(f"| {r['sample']} | {r['platform']} | {r['kind']} | {g} | {tt} | {tc} | | | | | | | |"); continue
            b = p["best"]
            L.append(f"| {r['sample']} | {r['platform']} | {r['kind']} | {g} | {tt} | {tc} | "
                     f"{e['raw']['evaluated']:,} ({e['unplaced_pass']:,}) | {f3(e['raw']['f1'])} | {p['tp']} / {p['fp']:,} | "
                     f"{f3(p['p'])} / {f3(p['r'])} / **{f3(p['f1'])}** | {f3(b['f1'])} at {b['p_somatic']:.3g} | "
                     f"{f3(p['ceiling'])} | {f3(n and n['f1'])} |")
    L.append("\n## Training of the d46 runs (chr2-22; 5 % of the node blocks validate; chr1 left out)\n")
    L.append("| run | GPUs | training tensors (somatic) | validation (somatic) | epochs | best epoch | val somatic AP | "
             "val somatic F1 | PASS threshold (val recall 0.9) | train min | tensors/s |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        x = r["d46"]; t = x["training"]; v = t.get("val", {})
        L.append(f"| `{x['run']}` | {t.get('gpus', '-')} | {t.get('train_tensors', 0):,} ({t.get('train_counts', {}).get('somatic', 0):,}) | "
                 f"{t.get('val_tensors', 0):,} ({t.get('val_counts', {}).get('somatic', 0):,}) | {t.get('epochs')} | "
                 f"{t.get('best_epoch')} | {f3(v.get('somatic_ap'))} | {f3(v.get('somatic_f1'))} | {f3(v.get('threshold_r09'))} | "
                 f"{t.get('train_minutes')} | {t.get('tensors_per_second')} |")
    return "\n".join(L) + "\n"


def main():
    rows = []
    for sample, (_, platforms) in SAMPLES.items():
        for k in KINDS:
            for p in platforms:
                r9, v9, note = d9_run(sample, p, k)
                rows.append({"sample": sample, "platform": p, "kind": k, "d9_note": note,
                             "d9": one(r9, v9), "d46": one(d46_run(sample, p, k), "vcf_chr1")})
    (HERE / "summary.json").write_text(json.dumps(rows, indent=1) + "\n")
    (HERE / "tables.md").write_text(tables(rows))
    print(tables(rows))


if __name__ == "__main__":
    main()
