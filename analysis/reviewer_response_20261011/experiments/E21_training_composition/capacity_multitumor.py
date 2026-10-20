"""E21 add-on (R1.7, 2026-10-20): (a) training / validation tensors per class of the 12 multi-tumor models, per training
set; (b) training and validation loss per epoch of the model-size runs (E13 on HG008T PacBio d46 SNV, and the HG008T
Illumina INDEL pair of 2026-10-06), with each run's chr1 PASS F1.

(a) reads train.log of /scratch/jshen/data/pansoma_net_v2_runs/multicancer{4,6}_<platform>_<kind>: the per-set
    "train N non= somatic= germline=" lines (training blocks, after the validation split and the partial left-outs; they
    add up to the run's "training tensors" line), the totals of training and validation, the left-out partials, the
    held-out test chromosomes and the best epoch.
(b) reads metrics.jsonl (train_loss, val.loss per epoch) and train.log (best epoch) of the runs in SIZE; PASS F1 from the
    run's vcf_chr1(_pon)/eval_bed/report.json (pon, at_threshold), else from the RIS summary (e13 entry; the seed-0
    50.4 M run has no local vcf_chr1).
Writes multitumor_training.tsv, capacity_loss.tsv and capacity_multitumor.md next to this script.
usage: /wanglab/jshen/anaconda3/bin/python capacity_multitumor.py
"""
import csv
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
RIS_SUMMARY = HERE.parents[2] / "benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/summary.json"
MULTI = [f"multicancer{n}_{p}_{k}" for n in (4, 6) for k in ("SNV", "INDEL") for p in ("Illumina", "PacBio", "ONT")]
SIZE = [  # (kind, dataset, recipe, run dir)
    ("SNV", "HG008T PacBio d46", "E13 (lr 1e-4, 20 epochs)", "ris_7fc7396_20261012/HG008_PacBio_SNV_d46_e13_8p7M"),
    ("SNV", "HG008T PacBio d46", "E13 (lr 1e-4, 20 epochs)", "ris_7fc7396_20261012/HG008_PacBio_SNV_d46_7fc7396"),
    ("SNV", "HG008T PacBio d46", "E13 (lr 1e-4, 20 epochs)", "ris_7fc7396_20261012/HG008_PacBio_SNV_d46_7fc7396_s1"),
    ("SNV", "HG008T PacBio d46", "E13 (lr 1e-4, 20 epochs)", "ris_7fc7396_20261012/HG008_PacBio_SNV_d46_7fc7396_s2"),
    ("SNV", "HG008T PacBio d46", "E13 (lr 1e-4, 20 epochs)", "ris_7fc7396_20261012/HG008_PacBio_SNV_d46_e13_199p6M"),
    ("INDEL", "HG008T Illumina HPRC v1.1 (d9)", "2026-10-06 (lr 4e-4, 12 epochs, sqrt weights)", "HG008_Illumina_INDEL_nopartial_b1024_small"),
    ("INDEL", "HG008T Illumina HPRC v1.1 (d9)", "2026-10-06 (lr 4e-4, 12 epochs, sqrt weights)", "HG008_Illumina_INDEL_nopartial_b1024"),
]
SET_RE = re.compile(r"^(\S+)/tensors/(SNV|INDEL): labels .*?, train ([\d,]+) non=([\d,]+) somatic=([\d,]+) germline=([\d,]+)")
TOT_RE = re.compile(r"^training tensors ([\d,]+) \(\{'non': (\d+), 'somatic': (\d+), 'germline': (\d+)\}\), "
                    r"validation ([\d,]+) \(\{'non': (\d+), 'somatic': (\d+), 'germline': (\d+)\}\)")
num = lambda s: int(s.replace(",", ""))  # noqa: E731


def log_lines(run):
    return (RUNS / run / "train.log").read_text().splitlines()


def best_epoch(lines):
    m = [re.match(r"done; best epoch (\d+)", x) for x in lines]
    return next(int(x.group(1)) for x in m if x)


def params(lines):
    return next(float(re.search(r"([\d.]+) M parameters, input", x).group(1)) for x in lines if "M parameters, input" in x)


def multitumor():
    rows = []
    for run in MULTI:
        L = log_lines(run)
        held = [re.sub(r"^test \(left out\) of (\S+)/tensors only: ", r"\1 ", x) for x in L if x.startswith("test (left out) of")]
        tot = next(TOT_RE.match(x) for x in L if TOT_RE.match(x))
        left = next(re.search(r"([\d,]+) training, ([\d,]+) validation", x) for x in L if x.startswith("left out"))
        common = dict(run=run, kind=run.rsplit("_", 1)[1], params_M=params(L), best_epoch=best_epoch(L))
        sets = []
        for x in L:
            m = SET_RE.match(x)
            if m:
                name = Path(m.group(1)).name
                test = next((h.split(" ", 1)[1] for h in held if Path(h.split(" ", 1)[0]).name == name), "")
                sets.append(dict(common, set=name, test_chroms=test, train=num(m.group(3)), train_non=num(m.group(4)),
                                 train_somatic=num(m.group(5)), train_germline=num(m.group(6)), val="", val_non="",
                                 val_somatic="", val_germline="", left_out_train="", left_out_val=""))
        g = [num(v) for v in tot.groups()]
        for c, i in (("train_non", 1), ("train_somatic", 2), ("train_germline", 3)):
            assert sum(s[c] for s in sets) == g[i], (run, c)
        sets.append(dict(common, set="all", test_chroms="", train=g[0], train_non=g[1], train_somatic=g[2],
                         train_germline=g[3], val=g[4], val_non=g[5], val_somatic=g[6], val_germline=g[7],
                         left_out_train=num(left.group(1)), left_out_val=num(left.group(2))))
        rows += sets
    return rows


def size_runs():
    ris = {e["run"]: e["score"] for e in json.loads(RIS_SUMMARY.read_text())["e13"]}
    out, curves = [], []
    for kind, data, recipe, rel in SIZE:
        d = RUNS / rel
        M = [json.loads(x) for x in (d / "metrics.jsonl").read_text().splitlines() if x.strip()]
        L = log_lines(rel)
        tl, vl = [m["train_loss"] for m in M], [m["val"]["loss"] for m in M]
        i = min(range(len(vl)), key=vl.__getitem__)
        rep = next((p for p in (d / "vcf_chr1/eval_bed/report.json", d / "vcf_chr1_pon/eval_bed/report.json")
                    if p.exists() and "pon" in json.loads(p.read_text())), None)
        if rep is not None:
            a = json.loads(rep.read_text())["pon"]["at_threshold"]
            f1, src = a["f1"], str(rep.relative_to(RUNS))
        else:
            f1, src = ris[d.name]["pass"], "RIS summary.json e13"
        b = best_epoch(L)
        out.append(dict(kind=kind, dataset=data, recipe=recipe, run=d.name, params_M=params(L), epochs=len(M),
                        val_loss_min=vl[i], val_loss_min_epoch=i + 1, val_loss_last=vl[-1], last_over_min=vl[-1] / vl[i],
                        train_loss_last=tl[-1], best_epoch=b, val_loss_best=vl[b - 1], pass_f1=f1, f1_source=src))
        curves += [dict(run=d.name, params_M=params(L), epoch=e + 1, train_loss=tl[e], val_loss=vl[e]) for e in range(len(M))]
    return out, curves


def write_tsv(name, rows):
    with open(HERE / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    mt = multitumor()
    size, curves = size_runs()
    write_tsv("multitumor_training.tsv", mt)
    write_tsv("capacity_loss.tsv", curves)
    n = lambda x: f"{x:,}" if isinstance(x, int) else x  # noqa: E731
    L = ["# Model size and multi-tumor training sets (generated by capacity_multitumor.py)", "",
         "## A. Model size: training and validation loss (cross-entropy with the class weights) and chr1 PASS F1", "",
         "PASS F1 = chr1, benchmark BED, after population filtering, at the validation recall-0.9 threshold. The INDEL pair "
         "was scored with the INDEL population filter of 2026-10-06 (AF >= 0.01), the same for both runs.", "",
         "| kind | dataset | recipe | run | params (M) | val loss min (epoch) | val loss last / min | train loss last | "
         "selected epoch (val loss) | PASS F1 |", "|---|---|---|---|---:|---|---:|---:|---|---:|"]
    for r in size:
        L.append(f"| {r['kind']} | {r['dataset']} | {r['recipe']} | `{r['run']}` | {r['params_M']} | "
                 f"{r['val_loss_min']:.3f} ({r['val_loss_min_epoch']}) | {r['last_over_min']:.2f} | {r['train_loss_last']:.4f} | "
                 f"{r['best_epoch']} ({r['val_loss_best']:.3f}) | {r['pass_f1']:.4f} |")
    L += ["", "Validation loss per epoch:", ""]
    for r in size:
        v = [c["val_loss"] for c in curves if c["run"] == r["run"]]
        L.append(f"- `{r['run']}` ({r['params_M']} M): " + " ".join(f"{x:.2f}" for x in v))
    L += ["", "## B. Multi-tumor models: training tensors per set and class, validation per run", "",
          "Training = chr2-22 (HCC1395T: chr2-19, 21, 22) minus the validation blocks (5 % of the 20,000-node blocks) and minus "
          "the partial somatic matches (left out); validation is reported per run only (train.log does not split it by set).", "",
          "| run | params (M) | set | test (held out) | train non / somatic / germline | val non / somatic / germline | "
          "left out train / val | selected epoch |", "|---|---:|---|---|---|---|---|---:|"]
    for r in mt:
        val = f"{n(r['val_non'])} / {n(r['val_somatic'])} / {n(r['val_germline'])}" if r["set"] == "all" else ""
        lo = f"{n(r['left_out_train'])} / {n(r['left_out_val'])}" if r["set"] == "all" else ""
        L.append(f"| `{r['run']}` | {r['params_M']} | {'**all**' if r['set'] == 'all' else r['set']} | {r['test_chroms']} | "
                 f"{n(r['train_non'])} / {n(r['train_somatic'])} / {n(r['train_germline'])} | {val} | {lo} | "
                 f"{r['best_epoch'] if r['set'] == 'all' else ''} |")
    (HERE / "capacity_multitumor.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[:20]))


if __name__ == "__main__":
    main()
