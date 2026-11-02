#!/usr/bin/env python3
"""R1.8 / Q3: how much do off-reference tensors without a truth match (reason off_reference_no_truth_match; -1 in
training, scored as non = 0 in validation and in predict's chr1 test) weigh in model selection and the PASS threshold?

Read-only, from stored files of the main models:
  SNV   ris_7fc7396_20261012/<S>_<P>_SNV_d46_normbq{,_s1,_s2}   (best epoch = highest pre-PoN per-tensor somatic F1;
        PASS = val somatic_at_recall 0.9)
  INDEL release_s{0,1,2}_<P>_INDEL                                (best epoch = highest post-PoN per-tensor somatic F1;
        PASS = val post-PoN callable recall 0.8, "ppr8")
Per model: train.log (validation off-reference non count, class counts), metrics.jsonl (per-epoch selection score,
margin best vs runner-up, the selection and PASS thresholds of the best epoch). No per-tensor validation predictions
and no per-epoch checkpoints are stored, so the validation effect cannot be recomputed; as a proxy the stored chr1
test predictions of the selected checkpoint (same model, same tensor population) give the per-tensor somatic F1 at the
selection threshold with and without the off-reference non tensors, and how many of them pass each threshold.
"""
import gzip
import json
import re
from pathlib import Path

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
OFF = "off_reference_no_truth_match"


def models():
    for d in sorted((R / "ris_7fc7396_20261012").glob("*_SNV_d46_normbq*")):
        yield "SNV", d
    for d in sorted(R.glob("release_s[012]_*_INDEL")):
        yield "INDEL", d


def log_facts(d):
    text = (d / "train.log").read_text()
    off = int(re.search(r"validation includes ([\d,]+) off-reference", text).group(1).replace(",", ""))
    m = re.search(r"validation ([\d,]+) \(\{'non': (\d+), 'somatic': (\d+), 'germline': (\d+)\}\)", text)
    best = re.search(r"done; best epoch (\d+)", text)
    return dict(val_offref_non=off, val_tensors=int(m.group(1).replace(",", "")), val_non=int(m.group(2)),
                val_somatic=int(m.group(3)), val_germline=int(m.group(4)), best_epoch=int(best.group(1)))


def epoch_facts(d, kind):
    rows = [json.loads(l) for l in open(d / "metrics.jsonl")]
    rows = [r for r in rows if "val" in r]
    score = (lambda v: v["postpon"]["f1"]) if kind == "INDEL" else (lambda v: v["thresholded"]["somatic"]["f1"])
    s = [(score(r["val"]), r["epoch"]) for r in rows]
    s_sorted = sorted(s, reverse=True)
    best = max(rows, key=lambda r: score(r["val"]))
    v = best["val"]
    if kind == "INDEL":
        t_sel = v["somatic_threshold_postpon"]
        t_pass = v["somatic_at_recall_postpon_callable"]["0.8"]["threshold"]
        fp_val = v["postpon"]["calls"] - round(v["postpon"]["precision"] * v["postpon"]["calls"])
    else:
        t_sel = v["threshold"]
        t_pass = v["somatic_at_recall"]["0.9"]["threshold"]
        cm = v["thresholded"]["confusion"]  # rows = true non/somatic/germline, cols = predicted
        fp_val = cm[0][1] + cm[2][1]
    return dict(best_epoch_metrics=best["epoch"], best_score=round(s_sorted[0][0], 4),
                runner_up_epoch=s_sorted[1][1], runner_up_score=round(s_sorted[1][0], 4),
                margin=round(s_sorted[0][0] - s_sorted[1][0], 4), t_sel=t_sel, t_pass=t_pass, val_fp_at_sel=int(fp_val),
                epochs=len(rows))


def prf(tp, fp, npos):
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / npos if npos else 0.0
    return round(p, 4), round(r, 4), round(2 * p * r / (p + r), 4) if p + r else 0.0


def test_facts(path, kind, t_sel, t_pass):
    lab, ps, off, pon = [], [], [], []
    with gzip.open(path, "rt") as f:
        for line in f:
            r = json.loads(line)
            if not r["in_test"] or r["test_label"] is None or r["test_label"] < 0:
                continue
            lab.append(r["test_label"])
            ps.append(r["p_somatic"])
            off.append(r["reason"] == OFF and r["test_label"] == 0)
            pon.append(r.get("pon_tag", 0))
    npos = sum(1 for x in lab if x == 1)
    out = dict(test_tensors=len(lab), test_non=sum(1 for x in lab if x == 0), test_offref_non=sum(off), positives=npos)
    post = kind == "INDEL"  # INDEL selection is post-PoN: a PoN-tagged tensor is never called
    for name, t in (("sel", t_sel), ("pass", t_pass)):
        called = [p >= t and not (post and g == 1) for p, g in zip(ps, pon)]
        tp = sum(1 for c, x in zip(called, lab) if c and x == 1)
        fp = sum(1 for c, x in zip(called, lab) if c and x != 1)
        fp_off = sum(1 for c, o in zip(called, off) if c and o)
        out[f"{name}_tp"], out[f"{name}_fp"], out[f"{name}_fp_offref"] = tp, fp, fp_off
        out[f"{name}_prf_with"] = prf(tp, fp, npos)
        out[f"{name}_prf_without_offref"] = prf(tp, fp - fp_off, npos)
    return out


def main():
    res = {}
    for kind, d in models():
        rec = dict(kind=kind, **log_facts(d), **epoch_facts(d, kind))
        assert rec["best_epoch"] == rec["best_epoch_metrics"], (d, rec)
        tests = {}
        for p in sorted(d.glob("test_chr1*/*.predictions.ndjson.gz")):
            name = p.name.split(".tensors.")[0]
            if not name.startswith(("HG008T", "COLO829T")):
                continue
            tests[name] = test_facts(p, kind, rec["t_sel"], rec["t_pass"])
        rec["tests"] = tests
        res[d.name] = rec
        print(d.name, json.dumps({k: v for k, v in rec.items() if k != "tests"}), flush=True)
        for n, t in tests.items():
            print("   ", n, json.dumps(t), flush=True)
    Path(__file__).with_suffix(".json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
