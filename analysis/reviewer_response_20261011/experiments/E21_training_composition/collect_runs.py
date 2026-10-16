"""E21 step 2: one row per trained model (36 single-set runs on the three reference arms + 6 COLO829T purity joint
runs): training / validation tensors per class, tensors left out by --ignore-reasons, epochs, best epoch,
validation somatic AP, PASS threshold (validation recall 0.9), parameter count. Standard library only.

    python collect_runs.py            # writes runs_training.tsv and runs_training.md next to this script

Sources (read-only):
- the 36 RIS runs of 2026-10-12: /scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/<run>/{train.log,
  metrics.jsonl, args.json};
- the six COLO829T HPRC v1.1 d9 runs (COLO829_<platform>_SNV_h90, 2026-10-03; COLO829_<platform>_INDEL_small,
  2026-10-06), whose run directories are on RIS only: the train.log copies in
  analysis/benchmark_7fc7396_20261012/*/runs/<run>/ and the training block of that folder's summary.json (collected
  on RIS from metrics.jsonl).
"""
import csv
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path("/scratch/jshen/Github/Pansoma")
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
BENCH = REPO / "analysis/benchmark_7fc7396_20261012"
PLATFORMS = {"HG008": ("Illumina", "PacBio", "ONT"), "COLO829": ("Illumina", "fiberseq", "ONT")}


def runs():
    """(run, arm, sample, platform, kind) of the 42 models."""
    out = []
    for s, plats in PLATFORMS.items():
        for p in plats:
            for k in ("SNV", "INDEL"):
                if s == "HG008":
                    d9 = "HG008_Illumina_SNV_7fc7396_ddp4_b1024" if (p, k) == ("Illumina", "SNV") else f"HG008_{p}_{k}_7fc7396"
                else:
                    d9 = f"COLO829_{p}_SNV_h90" if k == "SNV" else f"COLO829_{p}_INDEL_small"
                sample = s + "T"
                out.append((d9, "d9", sample, p, k))
                out.append((f"{s}_{p}_{k}_d46_7fc7396", "d46", sample, p, k))
                out.append((f"{s}_{p}_{k}_linear_7fc7396", "linear", sample, p, k))
                if s == "COLO829":
                    out.append((f"COLO829_{p}_{k}_joint_7fc7396" + ("_gpu4" if k == "SNV" else ""), "d9 joint P100+P75+P50+P25",
                                sample, p, k))
    return out


def find_log(run):
    p = RIS / run / "train.log"
    if p.exists():
        return p, "ris_7fc7396_20261012"
    for q in sorted(BENCH.glob(f"*/runs/{run}/train.log")):
        return q, str(q.parent.parent.parent.relative_to(REPO))
    return None, None


def summary_training(run):
    """The training block of the benchmark summary.json files for a run (the RIS-only runs)."""
    for f in sorted(BENCH.glob("*/summary.json")):
        d = json.loads(f.read_text())
        stack = [d]
        while stack:
            x = stack.pop()
            if isinstance(x, dict):
                if x.get("run") == run and "training" in x:
                    return x["training"]
                stack.extend(x.values())
            elif isinstance(x, list):
                stack.extend(x)
    return None


def num(s):
    return int(s.replace(",", ""))


def parse_log(text):
    o = {}
    m = re.search(r"training tensors ([\d,]+) \((\{[^)]*\})\), validation ([\d,]+) \((\{[^)]*\})\)", text)
    if m:
        o["train"] = num(m.group(1))
        o["train_counts"] = json.loads(m.group(2).replace("'", '"'))
        o["val"] = num(m.group(3))
        o["val_counts"] = json.loads(m.group(4).replace("'", '"'))
    m = re.search(r"left out \(--ignore-reasons ([^)]*)\): ([\d,]+) training, ([\d,]+) validation", text)
    if m:
        o["ignore_reasons"], o["ignored_train"], o["ignored_val"] = m.group(1), num(m.group(2)), num(m.group(3))
    m = re.search(r"validation includes ([\d,]+) off-reference tensors", text)
    if m:
        o["val_offref_as_non"] = num(m.group(1))
    m = re.search(r"truth-level validation: ([\d,]+) somatic truth alleles in the validation region \(([\d,]+) with", text)
    if m:
        o["val_truth"], o["val_truth_with_tensor"] = num(m.group(1)), num(m.group(2))
    m = re.search(r"class weights (\{[^}]*\})", text)
    if m:
        o["class_weights"] = m.group(1).replace("'", "")
    m = re.search(r"\}, ([\d.]+) M parameters", text)
    if m:
        o["params_M"] = float(m.group(1))
    m = re.search(r"'depths': (\[[^]]*\]), 'dims': (\[[^]]*\])", text)
    if m:
        o["arch"] = f"depths {m.group(1)} dims {m.group(2)}"
    m = re.search(r"device cuda(?::\d+)? x (\d+)", text)
    o["gpus"] = int(m.group(1)) if m else 1
    m = re.search(r"(\d+) steps per epoch", text)
    if m:
        o["steps_per_epoch"] = int(m.group(1))
    o["epochs_run"] = len(re.findall(r"^epoch \d+/\d+:", text, re.M))
    m = re.search(r"^epoch \d+/(\d+):", text, re.M)
    if m:
        o["epochs"] = int(m.group(1))
    m = re.search(r"done; best epoch (\d+) \(f1 ([\d.]+), somatic threshold ([\d.]+)\)", text)
    if m:
        o["best_epoch"], o["best_val_f1"], o["threshold_best_f1"] = int(m.group(1)), float(m.group(2)), float(m.group(3))
    ep = re.findall(r"^epoch (\d+)/\d+: train loss ([\d.]+) .*?\| val loss ([\d.]+) \| somatic AP ([\d.]+)", text, re.M)
    if ep:
        ep = [(int(e), float(tl), float(vl), float(ap)) for e, tl, vl, ap in ep]
        best = next((x for x in ep if x[0] == o.get("best_epoch")), None)
        vmin = min(ep, key=lambda x: x[2])
        amax = max(ep, key=lambda x: x[3])
        o.update(val_loss_best_epoch=best[2] if best else None, val_loss_min=vmin[2], val_loss_min_epoch=vmin[0],
                 val_loss_last=ep[-1][2], train_loss_best_epoch=best[1] if best else None, train_loss_last=ep[-1][1],
                 val_ap_max=amax[3], val_ap_max_epoch=amax[0], val_ap_last=ep[-1][3])
    o["sets"] = re.findall(r"pansoma_v2_tensors/([^/]+)/tensors/(?:SNV|INDEL): labels", text)
    o["labels_created"] = re.findall(r"labels \('truth-labels', '([0-9a-f]+)'\) \(([^)]+)\)", text)
    return o


def metrics_best(run, best_epoch):
    p = RIS / run / "metrics.jsonl"
    if not p.exists():
        return None
    for line in p.read_text().splitlines():
        if line.strip():
            e = json.loads(line)
            if e["epoch"] == best_epoch:
                v = e["val"]
                return {"somatic_ap": v["somatic_ap"], "threshold_r09": (v["somatic_at_recall"].get("0.9") or {}).get("threshold"),
                        "val_somatic_f1": v["thresholded"]["somatic"]["f1"], "germline_ap": v.get("germline_ap"),
                        "truth": v.get("truth")}
    return None


COLUMNS = ["run", "arm", "sample", "platform", "kind", "sets", "params_M", "gpus", "train", "train_non", "train_somatic",
           "train_germline", "val", "val_non", "val_somatic", "val_germline", "ignored_train", "ignored_val",
           "val_offref_as_non", "val_truth", "val_truth_with_tensor", "epochs", "epochs_run", "best_epoch",
           "val_somatic_ap", "val_somatic_f1", "threshold_best_f1", "threshold_r09", "class_weights",
           "steps_per_epoch", "train_loss_best_epoch", "train_loss_last", "val_loss_best_epoch", "val_loss_min",
           "val_loss_min_epoch", "val_loss_last", "val_ap_max", "val_ap_max_epoch", "val_ap_last", "source",
           "metrics_source"]


def main():
    rows = []
    for run, arm, sample, platform, kind in runs():
        log, source = find_log(run)
        if log is None:
            rows.append(dict(run=run, arm=arm, sample=sample, platform=platform, kind=kind, source="MISSING"))
            continue
        o = parse_log(log.read_text(errors="replace"))
        r = dict(run=run, arm=arm, sample=sample, platform=platform, kind=kind, sets="+".join(o.get("sets", [])),
                 params_M=o.get("params_M"), gpus=o.get("gpus"), train=o.get("train"), val=o.get("val"),
                 ignored_train=o.get("ignored_train"), ignored_val=o.get("ignored_val"),
                 val_offref_as_non=o.get("val_offref_as_non"), val_truth=o.get("val_truth"),
                 val_truth_with_tensor=o.get("val_truth_with_tensor"), epochs=o.get("epochs"),
                 epochs_run=o.get("epochs_run"), best_epoch=o.get("best_epoch"),
                 threshold_best_f1=o.get("threshold_best_f1"), class_weights=o.get("class_weights"), source=source,
                 **{k: o.get(k) for k in ("train_loss_best_epoch", "train_loss_last", "val_loss_best_epoch", "val_loss_min",
                                          "val_loss_min_epoch", "val_loss_last", "val_ap_max", "val_ap_max_epoch",
                                          "val_ap_last", "steps_per_epoch")})
        for c in ("non", "somatic", "germline"):
            r[f"train_{c}"] = o.get("train_counts", {}).get(c)
            r[f"val_{c}"] = o.get("val_counts", {}).get(c)
        m = metrics_best(run, o.get("best_epoch"))
        if m:
            r.update(val_somatic_ap=m["somatic_ap"], val_somatic_f1=m["val_somatic_f1"], threshold_r09=m["threshold_r09"],
                     metrics_source="metrics.jsonl (best epoch)")
        else:
            t = summary_training(run) or {}
            v = t.get("val", {})
            r.update(val_somatic_ap=v.get("somatic_ap"), val_somatic_f1=v.get("somatic_f1"),
                     threshold_r09=v.get("threshold_r09"), metrics_source="benchmark summary.json (RIS metrics.jsonl)")
        r["_labels_created"] = o.get("labels_created")
        rows.append(r)
    with open(HERE / "runs_training.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r[k]) for k in COLUMNS})
    with open(HERE / "runs_training.labels_created.json", "w") as f:
        json.dump({r["run"]: r.get("_labels_created") for r in rows}, f, indent=1)

    def n(x):
        return "-" if x is None else f"{x:,}"

    def f3(x):
        return "-" if x is None else f"{x:.3f}"

    L = ["| run | arm | params (M) | GPUs | train non / somatic / germline | val non / somatic / germline | "
         "left out by --ignore-reasons (train / val) | val off-ref as non | epochs (best) | val somatic AP | "
         "PASS threshold (val R0.9) |", "|---|---|---:|---:|---|---|---|---:|---|---:|---:|"]
    for r in rows:
        L.append(f"| `{r['run']}` | {r['arm']} | {r.get('params_M') or '-'} | {r.get('gpus') or '-'} | "
                 f"{n(r.get('train_non'))} / {n(r.get('train_somatic'))} / {n(r.get('train_germline'))} | "
                 f"{n(r.get('val_non'))} / {n(r.get('val_somatic'))} / {n(r.get('val_germline'))} | "
                 f"{n(r.get('ignored_train'))} / {n(r.get('ignored_val'))} | {n(r.get('val_offref_as_non'))} | "
                 f"{r.get('epochs_run')}/{r.get('epochs')} ({r.get('best_epoch')}) | {f3(r.get('val_somatic_ap'))} | "
                 f"{f3(r.get('threshold_r09'))} |")
    (HERE / "runs_training.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
