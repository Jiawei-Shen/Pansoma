"""Tables of the COLO829T purity runs of 2026-10-12 (repo 7fc7396, the v2 README recipes): one joint model per platform
and kind, trained on the undiluted set and its P75 / P50 / P25 dilutions together (chr2-22), tested on chr1 of each
purity separately. Standard library only (runs on the login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py

Reads, per run: train.log and metrics.jsonl (training, best epoch), <run>/test_chr1_<SET>/*.metrics.json (chr1 per
tensor and per truth allele) and <run>/vcf_chr1_<SET>/eval_{bed,nobed}/report.json (ris/joint_pipeline.sh: vcf.sh
--set); for the comparison at P100 the earlier single-set runs <run>/test_chr1, <run>/<vcf dir>. Writes tables.md,
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
PLATFORMS = ("Illumina", "fiberseq", "ONT")
KINDS = ("SNV", "INDEL")
PURITIES = ("P100", "P75", "P50", "P25")


def set_name(platform, purity):
    return f"COLO829T_{platform}" + ("" if purity == "P100" else f"_{purity}")


def joint_run(platform, kind):
    return f"COLO829_{platform}_{kind}_joint_7fc7396" + ("_gpu4" if kind == "SNV" else "")  # SNV on 4 GPUs


# the earlier single-set P100 runs (same sets): SNV the same recipe (2026-10-03, the SNV PoN rule unchanged since);
# INDEL lr 4e-4 (2026-10-06), re-scored with the 7fc7396 PoN (ris/vcfeval_only.sh) into vcf_chr1_pon7fc7396
EARLIER = {("SNV", p): (f"COLO829_{p}_SNV_h90", "vcf_chr1") for p in PLATFORMS}
EARLIER.update({("INDEL", p): (f"COLO829_{p}_INDEL_small", "vcf_chr1_pon7fc7396") for p in PLATFORMS})


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
        out["gpu_peak_gib"] = max((e.get("gpu_peak_gib") or {}).get("allocated", 0) for e in epochs)
    best = next((e for e in epochs if e["epoch"] == out.get("best_epoch")), None)
    if best:
        v = best["val"]
        out["val"] = {"somatic_ap": v["somatic_ap"], "somatic_f1": v["thresholded"]["somatic"]["f1"],
                      "threshold_t": v["threshold"],
                      "threshold_r09": (v["somatic_at_recall"].get("0.9") or {}).get("threshold")}
    m = re.search(r"device cuda(?::\d+)? x (\d+)", log)
    if m:
        out["gpus"] = int(m.group(1))
    return out


def test(pred_dir):
    files = glob.glob(str(pred_dir / "*.metrics.json"))
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


def copy_files(run, vcf_dirs):
    d, o = RUNS / run, HERE / "runs" / run
    pairs = [(d / "train.log", "train.log")] + [(d / v / w / "report.txt", f"{v}.{w}.txt")
                                                 for v in vcf_dirs for w in ("eval_bed", "eval_nobed")]
    for src, dst in pairs:
        if src.is_file():
            o.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, o / dst)


def f3(x):
    return "-" if x is None else f"{x:.3f}"


def collect():
    res = {"joint": {}, "earlier": {}}
    for p in PLATFORMS:
        for k in KINDS:
            run = joint_run(p, k)
            d = RUNS / run
            if not d.is_dir():
                continue
            per = {}
            for q in PURITIES:
                s = set_name(p, q)
                per[q] = {"set": s, "test": test(d / f"test_chr1_{s}"),
                          "eval_bed": rtg(d / f"vcf_chr1_{s}", "eval_bed"), "eval_nobed": rtg(d / f"vcf_chr1_{s}", "eval_nobed")}
            res["joint"][f"{p} {k}"] = {"run": run, "training": training(d), "purity": per}
            copy_files(run, [f"vcf_chr1_{set_name(p, q)}" for q in PURITIES])
            er, ev = EARLIER[(k, p)]
            if (RUNS / er).is_dir():
                res["earlier"][f"{p} {k}"] = {"run": er, "vcf_dir": ev, "training": training(RUNS / er),
                                             "test": test(RUNS / er / "test_chr1"), "eval_bed": rtg(RUNS / er / ev, "eval_bed"),
                                             "eval_nobed": rtg(RUNS / er / ev, "eval_nobed")}
                copy_files(er, [ev])
    return res


def tables(res):
    J, E = res["joint"], res["earlier"]
    L = ["## chr1 PASS F1 after the PoN, per purity (rtg vcfeval in the BED; the chr1 curve's best F1 in brackets)\n",
         "| platform | kind | P100 | P75 | P50 | P25 |", "|---|---|---:|---:|---:|---:|"]
    for key, j in J.items():
        p, k = key.split()
        cells = []
        for q in PURITIES:
            e = j["purity"][q]["eval_bed"]
            cells.append(f"**{f3(e['pon']['f1'])}** ({f3(e['pon']['best']['f1'])})" if e and e.get("pon") else "-")
        L.append(f"| {p} | {k} | " + " | ".join(cells) + " |")
    L.append("\n## Training (one joint model per platform and kind: P100 + P75 + P50 + P25, chr2-22; 5 % of the node "
             "blocks validate, the same blocks in every set)\n")
    L.append("| run | GPUs | training tensors (somatic) | validation (somatic) | epochs | best epoch | val somatic AP | "
             "val somatic F1 (t) | PASS threshold (val recall 0.9) | train min | tensors/s |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for key, j in J.items():
        t = j["training"]; v = t.get("val", {})
        L.append(f"| `{j['run']}` | {t.get('gpus', '-')} | {t.get('train_tensors', 0):,} ({t.get('train_counts', {}).get('somatic', 0):,}) | "
                 f"{t.get('val_tensors', 0):,} ({t.get('val_counts', {}).get('somatic', 0):,}) | {t.get('epochs')} | "
                 f"{t.get('best_epoch')} | {f3(v.get('somatic_ap'))} | {f3(v.get('somatic_f1'))} ({f3(v.get('threshold_t'))}) | "
                 f"{f3(v.get('threshold_r09'))} | {t.get('train_minutes')} | {t.get('tensors_per_second')} |")
    L.append("\n## chr1 per purity: tensors, truth alleles, rtg vcfeval in the BED (PASS = validation recall-0.9 threshold)\n")
    L.append("| platform | kind | purity | test tensors (somatic) | somatic AP | truth with a tensor (ceiling) | "
             "PASS evaluated (unplaced) | raw P / R / F1 | PoN TP / FP | PoN P / R / F1 | PoN chr1 best F1 at p (P, R) | "
             "PoN ceiling | no BED: PoN F1 |")
    L.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for key, j in J.items():
        p, k = key.split()
        for q in PURITIES:
            x = j["purity"][q]; t, e, n = x["test"], x["eval_bed"], x["eval_nobed"]
            if not e or not e.get("pon"):
                L.append(f"| {p} | {k} | {q} | {format(t['tensors'], ',') if t else '-'} | | | | | | | | | |"); continue
            raw, pon, b = e["raw"], e["pon"], e["pon"]["best"]
            L.append(f"| {p} | {k} | {q} | {t['tensors']:,} ({t['somatic_tensors']:,}) | {f3(t['somatic_ap'])} | "
                     f"{t['with_tensor']} / {t['truth_alleles']} ({f3(t['ceiling'])}) | {raw['evaluated']:,} ({e['unplaced_pass']:,}) | "
                     f"{f3(raw['p'])} / {f3(raw['r'])} / {f3(raw['f1'])} | {pon['tp']} / {pon['fp']:,} | "
                     f"{f3(pon['p'])} / {f3(pon['r'])} / **{f3(pon['f1'])}** | {f3(b['f1'])} at {b['p_somatic']:.3g} "
                     f"({f3(b['precision'])}, {f3(b['recall'])}) | {f3(pon['ceiling'])} | "
                     f"{f3(n['pon']['f1']) if n and n.get('pon') else '-'} |")
    L.append("\n## P100: the joint model against the earlier single-set model of the same set (chr1, PoN, in the BED)\n")
    L.append("| platform | kind | joint: PASS F1 (best) | single-set run | single-set: PASS F1 (best) | Δ PASS F1 |")
    L.append("|---|---|---:|---|---:|---:|")
    for key, j in J.items():
        p, k = key.split()
        e = j["purity"]["P100"]["eval_bed"]; o = E.get(key)
        oe = o and o["eval_bed"]
        mine = e["pon"]["f1"] if e and e.get("pon") else None
        theirs = oe["pon"]["f1"] if oe and oe.get("pon") else None
        L.append(f"| {p} | {k} | {f3(mine)} ({f3(e['pon']['best']['f1']) if mine is not None else '-'}) | "
                 f"`{o['run']}/{o['vcf_dir']}` | {f3(theirs)} ({f3(oe['pon']['best']['f1']) if theirs is not None else '-'}) | "
                 f"{'-' if mine is None or theirs is None else f'{mine - theirs:+.3f}'} |" if o else
                 f"| {p} | {k} | {f3(mine)} | - | - | - |")
    return "\n".join(L) + "\n"


def main():
    res = collect()
    (HERE / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    (HERE / "tables.md").write_text(tables(res))
    print(tables(res))


if __name__ == "__main__":
    main()
