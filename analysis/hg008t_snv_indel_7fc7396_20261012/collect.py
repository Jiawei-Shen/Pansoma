"""Tables of the HG008T SNV and INDEL runs of 2026-10-12 (repo 7fc7396, the README recipes) and of the earlier runs
they are compared with. Standard library only (runs on the login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py

Reads, per run: train.log (training counts, best epoch), metrics.jsonl (validation of the best epoch), the chr1 test
<run>/test_chr1/*.metrics.json (per tensor and truth-level), and the rtg vcfeval reports <run>/<vcf dir>/eval_bed and
eval_nobed/report.json (vcf.sh / vcfeval_only.sh). Writes tables.md and summary.json next to this script, and copies
each run's train.log and rtg report.txt files to runs/<run>/ (train.log, <vcf dir>.eval_bed.txt, <vcf dir>.eval_nobed.txt).
"""
import glob
import json
import os
import re
import shutil
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
HERE = Path(__file__).resolve().parent
PLATFORMS = ("Illumina", "PacBio", "ONT")
# (label, run name, vcf directory): the runs of this analysis, then the earlier ones
NEW = [(f"{k} {p}", f"HG008_{p}_{k}_7fc7396" + ("_ddp4_b1024" if (k, p) == ("SNV", "Illumina") else ""), "vcf_chr1")
       for k in ("SNV", "INDEL") for p in PLATFORMS]  # HG008T Illumina SNV: 4 H100s (DDP), batch 1024 per GPU
OLD = ([(f"SNV {p}", f"HG008_{p}_SNV_h90", "vcf_chr1") for p in PLATFORMS]
       + [(f"INDEL {p}", f"HG008_{p}_INDEL_small", "vcf_chr1_pon7fc7396") for p in PLATFORMS]
       + [(f"INDEL {p}", f"HG008_{p}_INDEL_small", "vcf_chr1") for p in PLATFORMS])


def load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def training(d):
    """Counts and the best epoch's validation of a run, from train.log and metrics.jsonl."""
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
    m = re.search(r"truth-level validation: (\d+) somatic truth alleles in the validation region \((\d+) with", log)
    if m:
        out["val_truth"], out["val_truth_with_tensor"] = int(m.group(1)), int(m.group(2))
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
        out["train_minutes"] = round(sum(e["train_seconds"] for e in epochs) / 60, 1)
        out["tensors_per_second"] = round(sum(e["tensors_per_second"] for e in epochs[1:] or epochs)
                                          / len(epochs[1:] or epochs))
        out["gpu_peak_gib"] = max((e.get("gpu_peak_gib") or {}).get("allocated", 0) for e in epochs)
    best = next((e for e in epochs if e["epoch"] == out.get("best_epoch")), None)
    if best:
        v = best["val"]
        r9 = v["somatic_at_recall"].get("0.9") or {}
        out["val"] = {"somatic_ap": v["somatic_ap"], "somatic_f1": v["thresholded"]["somatic"]["f1"],
                      "somatic_support": v["thresholded"]["somatic"]["support"], "threshold_t": v["threshold"],
                      "threshold_r09": r9.get("threshold"), "truth_f1": v["truth"]["at_threshold"]["f1"],
                      "truth_best_f1": v["truth"]["best"]["f1"], "truth_ceiling": v["truth"]["ceiling"]}
    return out


def test(d):
    """The chr1 test (predict's .metrics.json): per tensor and against the truth VCF."""
    files = glob.glob(str(d / "test_chr1" / "*.metrics.json"))
    m = load(files[0]) if files else None
    if not m:
        return {}
    t = m["truth"]
    out = {"tensors": m["tensors"], "somatic_tensors": m["thresholded"]["somatic"]["support"],
           "somatic_ap": m["somatic_ap"], "threshold_t": m["threshold"],
           "truth_alleles": t["truth_alleles"], "with_tensor": t["with_tensor"], "ceiling": t["ceiling"],
           "truth_ap": t["ap"], "truth_at_t": t["at_threshold"], "truth_best": t["best"]}
    if "offref_rescue" in m:
        out["offref_rescue"] = m["offref_rescue"]
    return out


def rtg(d, vcf_dir, where):
    r = load(d / vcf_dir / where / "report.json")
    if not r:
        return None
    out = {"truth_records": r["truth_records"], "threshold": r.get("threshold"), "unplaced_pass": r.get("unplaced_pass"),
           "truth_pon_tagged": (r.get("truth_pon_tags") or {}).get("tagged")}
    for s in ("raw", "pon"):
        if s not in r:
            continue
        a, b = r[s]["at_threshold"], r[s]["best"]
        out[s] = {"baseline": r[s]["baseline"], "evaluated": a["tp_call"] + a["fp"], "tp": a["tp_baseline"],
                  "fp": a["fp"], "fn": a["fn"], "p": a["precision"], "r": a["recall"], "f1": a["f1"],
                  "ceiling": r[s]["ceiling"],
                  "best": {k: b[k] for k in ("f1", "p_somatic", "precision", "recall", "tp_baseline", "fp", "fn")}}
    return out


def collect(label, name, vcf_dir):
    d = RUNS / name
    return {"label": label, "run": name, "vcf_dir": vcf_dir, "training": training(d), "test": test(d),
            "eval_bed": rtg(d, vcf_dir, "eval_bed"), "eval_nobed": rtg(d, vcf_dir, "eval_nobed")}


def f3(x):
    return "-" if x is None else f"{x:.3f}"


def prf(s):
    return f"{f3(s['p'])} / {f3(s['r'])} / **{f3(s['f1'])}**" if s else "-"


def best(s):
    if not s:
        return "-"
    b = s["best"]
    return f"{f3(b['f1'])} at {b['p_somatic']:.3g} ({f3(b['precision'])}, {f3(b['recall'])}; {b['tp_baseline']} / {b['fp']:,})"


def tables(new, old):
    L = []
    L.append("## Training (chr2-22, 5 % of the node blocks validate; chr1 left out)\n")
    L.append("| run | training tensors (somatic) | validation (somatic) | epochs | best epoch | val somatic AP | "
             "val somatic F1 (t) | PASS threshold (val recall 0.9) | train min | tensors/s | GPU peak GiB (allocated) |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for x in new:
        t = x["training"]
        v = t.get("val", {})
        L.append(f"| {x['label']} | {t.get('train_tensors', 0):,} ({t.get('train_counts', {}).get('somatic', 0):,}) | "
                 f"{t.get('val_tensors', 0):,} ({t.get('val_counts', {}).get('somatic', 0):,}) | {t.get('epochs')} | "
                 f"{t.get('best_epoch')} | {f3(v.get('somatic_ap'))} | {f3(v.get('somatic_f1'))} ({f3(v.get('threshold_t'))}) | "
                 f"{f3(v.get('threshold_r09'))} | {t.get('train_minutes')} | {t.get('tensors_per_second')} | "
                 f"{t.get('gpu_peak_gib')} |")
    L.append("\n## chr1 test per tensor and against the truth alleles (predict, `metrics.truth_report`)\n")
    L.append("| run | test tensors (somatic) | somatic AP | truth alleles | with a tensor (ceiling) | "
             "truth F1 at predict's t (P, R) | truth best F1 (P, R) |")
    L.append("|---|---:|---:|---:|---:|---:|---:|")
    for x in new:
        t = x["test"]
        if not t:
            L.append(f"| {x['label']} | - | | | | | |"); continue
        a, b = t["truth_at_t"], t["truth_best"]
        L.append(f"| {x['label']} | {t['tensors']:,} ({t['somatic_tensors']:,}) | {f3(t['somatic_ap'])} | {t['truth_alleles']} | "
                 f"{t['with_tensor']} ({f3(t['ceiling'])}) | {f3(a['f1'])} ({f3(a['precision'])}, {f3(a['recall'])}) | "
                 f"{f3(b['f1'])} ({f3(b['precision'])}, {f3(b['recall'])}) |")
    for where, title in (("eval_bed", "in the GIAB benchmark BED"), ("eval_nobed", "whole chr1, no BED")):
        L.append(f"\n## chr1 calls, rtg vcfeval, {title} (PASS = validation recall-0.9 threshold)\n")
        L.append("| run | truth | PASS evaluated (unplaced) | raw P / R / F1 | PoN TP / FP | PoN P / R / F1 | "
                 "PoN chr1 best F1 at p (P, R; TP / FP) | PoN ceiling | truth tagged by the PoN |")
        L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for x in new:
            e = x[where]
            if not e:
                L.append(f"| {x['label']} | - | | | | | | | |"); continue
            raw, pon = e.get("raw"), e.get("pon")
            L.append(f"| {x['label']} | {e['truth_records'] if where == 'eval_nobed' else raw['baseline']} | "
                     f"{raw['evaluated']:,} ({e['unplaced_pass']:,}) | {prf(raw)} | {pon['tp']} / {pon['fp']:,} | {prf(pon)} | "
                     f"{best(pon)} | {f3(pon['ceiling'])} | {e['truth_pon_tagged']} of {e['truth_records']} |")
    L.append("\n## Against the earlier runs (chr1, PoN, in the BED: PASS F1, chr1 best F1)\n")
    L.append("| run | this run (7fc7396) | earlier run | earlier PASS F1 | earlier best F1 | Δ PASS F1 |")
    L.append("|---|---:|---|---:|---:|---:|")
    for x in new:
        e = x["eval_bed"]
        mine = e["pon"]["f1"] if e and e.get("pon") else None
        mine_s = f"{f3(mine)} (best {f3(e['pon']['best']['f1'])})" if mine is not None else "-"
        for o in old:
            if o["label"] != x["label"]:
                continue
            oe = o["eval_bed"]
            of = oe["pon"]["f1"] if oe and oe.get("pon") else None
            L.append(f"| {x['label']} | {mine_s} | `{o['run']}/{o['vcf_dir']}` | {f3(of)} | "
                     f"{f3(oe['pon']['best']['f1']) if of is not None else '-'} | "
                     f"{'-' if mine is None or of is None else f'{mine - of:+.3f}'} |")
    return "\n".join(L) + "\n"


def copy_files(label, name, vcf_dir):
    """The run's train.log and rtg report.txt files, into runs/<run>/."""
    d, o = RUNS / name, HERE / "runs" / name
    for src, dst in [(d / "train.log", "train.log")] + [
            (d / vcf_dir / w / "report.txt", f"{vcf_dir}.{w}.txt") for w in ("eval_bed", "eval_nobed")]:
        if src.is_file():
            o.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, o / dst)


def main():
    new = [collect(*x) for x in NEW]
    old = [collect(*x) for x in OLD]
    for x in NEW + OLD:
        copy_files(*x)
    (HERE / "summary.json").write_text(json.dumps({"runs": new, "earlier": old}, indent=1) + "\n")
    (HERE / "tables.md").write_text(tables(new, old))
    print(tables(new, old))


if __name__ == "__main__":
    main()
