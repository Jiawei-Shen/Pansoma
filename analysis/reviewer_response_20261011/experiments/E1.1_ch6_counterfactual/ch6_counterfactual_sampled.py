"""E1.1 on a class-stratified sample of the chr1 test tensors (for the GRCh38-only graph's 218-414 k test tensors).

The same modes, edits, thresholds and outputs as ch6_counterfactual.py (imported: the same Counterfactual dataset,
training statistics, LOW rule and per-tensor PoN), on: every somatic tensor (label 1), every off-reference tensor,
and at most --max-germline (30,000) on-reference germline (2) and --max-non (20,000) on-reference non (0) random
tensors (seed 0). Each sampled
tensor of class c carries the weight n_c / kept_c, and every count, average precision and false-call number below is
weighted (with all weights 1 this gives exactly ch6_counterfactual.py's numbers: --check). Truth-level TP are exact
(every somatic tensor is kept); FP, germline and non counts are estimates.

usage: python ch6_counterfactual_sampled.py --run RUN --set SET --out DIR [--max-non 20000] [--max-germline 30000] [--threads 16]"""
import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch6_counterfactual as cf  # noqa: E402
from pansoma_net_v2 import metrics  # noqa: E402
from pansoma_net_v2.data import somatic_truth  # noqa: E402
from pansoma_net_v2.model import PansomaNetV2  # noqa: E402

OWN_CACHE = Path("/scratch/jshen/data/reviewer_experiments_20261013/E1.1_ch6_counterfactual/index_cache")


def wap(positive, score, w):
    """Weighted average precision (sklearn's step sum; ties together); positives must have weight 1."""
    positive = np.asarray(positive, bool)
    n_pos = int(positive.sum())
    if n_pos == 0:
        return 0.0
    order = np.argsort(-score, kind="stable")
    s = score[order]
    tp = np.cumsum(positive[order].astype(np.float64))
    tot = np.cumsum(w[order].astype(np.float64))
    last = np.r_[s[1:] != s[:-1], True]
    tp, tot = tp[last], tot[last]
    return float(np.sum(np.diff(np.r_[0.0, tp / n_pos]) * (tp / tot)))


def wtruth(truth, other, matches, lab, ps, w, thr):
    """metrics.truth_report's at-threshold counts and AP with weighted negatives (false calls)."""
    other = set(other) - truth
    counted = truth | other
    best, neg, negw = {}, [], []
    for k in range(len(lab)):
        if lab[k] < 0:
            continue
        if lab[k] == 1:
            for key in set(matches[k]) & counted:
                best[key] = max(float(ps[k]), best.get(key, -np.inf))
        else:
            neg.append(ps[k])
            negw.append(w[k])
    own = np.array([v for key, v in best.items() if key in truth])
    oth = np.array([v for key, v in best.items() if key in other])
    neg, negw = np.array(neg), np.array(negw)
    n = len(truth)
    tp, otp = int((own >= thr).sum()), int((oth >= thr).sum())
    fp = float(negw[neg >= thr].sum()) if len(neg) else 0.0
    prec, rec, f1 = (float(v) for v in metrics.truth_prf(n, tp, otp, fp))
    score = np.concatenate([own, oth, neg]).astype(np.float64)
    kind = np.r_[np.zeros(len(own)), np.ones(len(oth)), np.full(len(neg), 2)]
    wt = np.r_[np.ones(len(own) + len(oth)), negw]
    order = np.argsort(-score, kind="stable")
    s, kd, ww = score[order], kind[order], wt[order]
    last = np.r_[s[1:] != s[:-1], True] if len(s) else np.zeros(0, bool)
    ctp, cotp, cfp = np.cumsum(kd == 0)[last], np.cumsum(kd == 1)[last], np.cumsum(ww * (kd == 2))[last]
    p_, r_, _ = metrics.truth_prf(n, ctp, cotp, cfp)
    ap = float(np.sum(np.diff(np.r_[0.0, r_]) * p_)) if n else 0.0
    return dict(truth_alleles=n, with_tensor=len(own), ap=ap, tp=tp, other_tp=otp, fp=round(fp), fn=n - tp,
                precision=prec, recall=rec, f1=f1)


def summarize_w(mode, probs, ctx):
    lab, thr, test, w = ctx["labels"], ctx["thr"], ctx["test"], ctx["weights"]
    pr, lt, wt, off = probs[test], lab[test], w[test], ctx["off"][test]
    ps = pr[:, 1]
    am = pr.argmax(1) == 1
    call = ps >= thr
    r = dict(mode=mode, tensors=int(len(lt)), weighted_tensors=float(wt.sum()))
    r["somatic_ap"] = wap(lt == 1, ps, wt)
    r["germline_ap"] = wap(lt == 2, pr[:, 2], np.where(lt == 2, 1.0, wt))  # germline weights: an approximation
    for name, m in (("somatic", lt == 1), ("germline", lt == 2), ("non_onref", (lt == 0) & ~off),
                    ("non_offref", (lt == 0) & off)):
        r[name] = dict(n=int(round(wt[m].sum())), ge_thr=int(round(wt[call & m].sum())),
                       argmax_somatic=int(round(wt[am & m].sum())), sampled=int(m.sum()),
                       p_median=float(np.median(ps[m])) if m.any() else None,
                       p_mean=float(np.average(ps[m], weights=wt[m])) if m.any() else None)
    r["calls_ge_thr"] = int(round(wt[call].sum()))
    tp_t = int((call & (lt == 1)).sum())
    r["per_tensor"] = dict(tp=tp_t, fp=r["calls_ge_thr"] - tp_t, precision=tp_t / max(r["calls_ge_thr"], 1),
                           recall=tp_t / max(int((lt == 1).sum()), 1))
    base = ctx.get("base_probs")
    if base is not None:
        bt = base[test][:, 1]
        r["flips_vs_base"] = dict(new_calls=int(round(wt[call & ~(bt >= thr)].sum())),
                                  lost_calls=int(round(wt[~call & (bt >= thr)].sum())),
                                  somatic_delta_p_median=float(np.median(ps[lt == 1] - bt[lt == 1])) if (lt == 1).any() else None,
                                  germline_delta_p_median=float(np.median(ps[lt == 2] - bt[lt == 2])) if (lt == 2).any() else None)
    r["truth"] = wtruth(ctx["truth"], ctx["other"], ctx["matches"], lt, ps, wt, thr)
    keep = ~off & ~ctx["pon"][test]
    r["truth_after_pon"] = wtruth(ctx["truth"], ctx["other"], ctx["matches"], np.where(keep, lt, -1), ps, wt, thr)
    r["germline_after_pon_ge_thr"] = int(round(wt[call & keep & (lt == 2)].sum()))
    r["non_after_pon_ge_thr"] = int(round(wt[call & keep & (lt == 0)].sum()))
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--modes", nargs="+", default=list(cf.MODES))
    ap.add_argument("--max-non", type=int, default=20000, help="at most this many on-reference non test tensors")
    ap.add_argument("--max-germline", type=int, default=30000, help="at most this many on-reference germline test tensors")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--train-sample", type=int, default=3000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--check", metavar="DIR", help="re-score ch6_counterfactual.py's saved probabilities of DIR with all "
                                                     "weights 1 and compare (no inference)")
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cf.log(f"run {args.run} set {args.set} threads {args.threads} host {os.uname().nodename} job {os.environ.get('SLURM_JOB_ID')}")
    ck = torch.load(cf.RUNS / args.run / "best.pth", map_location="cpu", weights_only=False)
    model = PansomaNetV2.from_checkpoint(ck).to(cf.DEVICE).eval()
    thr = float(ck["val"]["somatic_at_recall"]["0.9"]["threshold"])
    index = cf.open_index(cf.TENSORS / args.set / "tensors" / "SNV", OWN_CACHE)
    if index.path_count != ck.get("path_count_storage"):
        raise SystemExit("ch6 storage differs from the checkpoint's")
    a, reasons = index.arrays, index.meta["reasons"]
    pos = index.select(["chr1"], labelled=False)
    ignore = set(ck.get("args", {}).get("ignore_reasons") or [])
    lab_all = a["eval_label"].astype(np.int64).copy()
    lab_all[np.isin(a["reason"], [k for k, r in enumerate(reasons) if r in ignore])] = -1
    test_pos = pos[lab_all[pos] >= 0]
    rng = np.random.default_rng(0)
    keep, weight_of, sampling = [], {}, {}
    lt, offt = lab_all[test_pos], a["off_reference"][test_pos].astype(bool)
    for name, m in (("somatic", lt == 1), ("offref", offt & (lt != 1)), ("germline", (lt == 2) & ~offt),
                    ("non_onref", (lt == 0) & ~offt)):
        idx = np.flatnonzero(m)
        cap = {"germline": args.max_germline, "non_onref": args.max_non}.get(name)
        if args.check or cap is None or len(idx) <= cap:
            pick, w = idx, 1.0
        else:
            pick, w = np.sort(rng.choice(idx, cap, replace=False)), len(idx) / cap
        keep.append(pick)
        for p in test_pos[pick]:
            weight_of[int(p)] = w
        sampling[name] = dict(n=int(len(idx)), kept=int(len(pick)), weight=w)
    ev_pos = np.sort(test_pos[np.concatenate(keep)])
    if args.check:
        ev_pos = np.load(Path(args.check) / "positions.npy")
    if args.limit:
        ev_pos = ev_pos[:args.limit]
    test = np.isin(ev_pos, test_pos)
    weights = np.array([weight_of.get(int(p), 1.0) for p in ev_pos])
    labels, off = lab_all[ev_pos], a["off_reference"][ev_pos].astype(bool)
    cf.log(f"chr1 tensors {len(pos):,}; test {len(test_pos):,}; evaluated {len(ev_pos):,} ({Counter(labels[test].tolist())}); "
           f"sampling {json.dumps(sampling)}")
    work = Path(args.check) if args.check else out
    stats_file = work / "train_ch6_stats.json"
    if stats_file.exists():
        tstats = json.loads(stats_file.read_text())
    else:
        tstats = cf.train_ch6_stats(index, ck, args.train_sample)
        tstats["checkpoint_mean"] = float(ck["stats"]["path_count"]["mean"])
        tstats["checkpoint_std"] = float(ck["stats"]["path_count"]["std"])
        stats_file.write_text(json.dumps(tstats, indent=1) + "\n")
    params = dict(train_mean=int(round(tstats["mean"])), train_median=int(tstats["median"]))
    pon_file = work / "pon_tags.json"
    if pon_file.exists():
        tags = {int(k): v for k, v in json.loads(pon_file.read_text()).items()}
    else:
        tags = cf.pon_tag(cf.grch38_alleles(index, ev_pos[test & ~off]))
        pon_file.write_text(json.dumps({str(k): v for k, v in tags.items()}) + "\n")
    pon = np.array([bool(tags.get(int(p))) for p in ev_pos])
    truth = somatic_truth(index, chroms={"chr1"})
    other = somatic_truth(index, {"chr1"}, True)
    ctx = dict(labels=labels, thr=thr, test=test, off=off, pon=pon, weights=weights, truth=set(truth),
               other=set(other or {}), matches=[set(index.truth_of(p).tolist()) for p in ev_pos[test]])
    pos_file = out / "positions.npy"
    if not args.check:
        if pos_file.exists() and not np.array_equal(np.load(pos_file), ev_pos):
            raise SystemExit(f"{pos_file} holds other tensors")
        np.save(pos_file, ev_pos)
    rng = np.random.default_rng(1)
    test_idx = np.flatnonzero(test)
    donors = np.empty(len(ev_pos), np.int64)
    donors[test_idx] = test_idx[rng.permutation(len(test_idx))]
    donors[~test] = rng.choice(test_idx, int((~test).sum()))
    results = {}
    for mode in args.modes:
        f = work / f"probs_{mode}.npy"
        if mode == "train_median" and params["train_median"] == cf.ALL:
            cf.log("train_median is 100: skipped (as all_H)")
            continue
        if f.exists() and (work / "profile.npz").exists():
            probs = np.load(f)
            z = np.load(work / "profile.npz")
            prof = (z["site"], z["below"], z["distinct"])
        elif args.check:
            continue
        else:
            cf.log(f"{mode}: start")
            ds = cf.Counterfactual([(index, ev_pos)], mode, params, donors)
            if mode == "base":
                probs, prof = cf.run_mode(model, ds, args, profile=True)
                np.savez(out / "profile.npz", site=prof[0], below=prof[1], distinct=prof[2])
            else:
                probs = cf.run_mode(model, ds, args)
            np.save(f, probs.astype(np.float32))
        if mode == "base":
            ctx["base_probs"] = probs
            site = prof[0]
            offneg = test & off & (labels == 0) & (site > 0)
            if offneg.sum() >= 50:
                params["low"], params["low_rule"] = int(np.median(site[offneg])), f"median site code of {int(offneg.sum())} off-reference test tensors"
            else:
                params["low"], params["low_rule"] = int(np.percentile(site[test & (site > 0)], 5)), "5th percentile of the evaluated test tensors' site codes (fewer than 50 off-reference)"
            groups = {"somatic": test & (labels == 1), "germline": test & (labels == 2),
                      "non_onref": test & (labels == 0) & ~off, "non_offref": test & (labels == 0) & off}
            params["site_code"] = {g: dict(n=int(m.sum()), median=float(np.median(site[m])) if m.any() else None,
                                           p10=float(np.percentile(site[m], 10)) if m.any() else None,
                                           frac_site_100=float((site[m] == cf.ALL).mean()) if m.any() else None,
                                           frac_cells_below_100=float(np.mean(prof[1][m])) if m.any() else None,
                                           distinct_codes_pm10_median=float(np.median(prof[2][m])) if m.any() else None)
                                   for g, m in groups.items()}
            cf.log(f"LOW = {params['low']} ({params['low_rule']})")
        results[mode] = summarize_w(mode, probs, ctx)
        s = results[mode]
        cf.log(f"{mode}: AP {s['somatic_ap']:.4f} | somatic >=thr {s['somatic']['ge_thr']}/{s['somatic']['n']} | germline >=thr "
               f"{s['germline']['ge_thr']}/{s['germline']['n']} | truth TP {s['truth']['tp']} FP {s['truth']['fp']} F1 {s['truth']['f1']:.3f} | "
               f"after PoN F1 {s['truth_after_pon']['f1']:.3f}")
    report = dict(run=args.run, set=args.set, threshold_val_recall_0_9=thr, checkpoint_somatic_threshold=float(ck["somatic_threshold"]),
                  epoch=ck["epoch"], haplotypes=index.haplotypes, train_ch6=tstats, params=params, chr1_tensors=int(len(pos)),
                  test_tensors=int(len(test_pos)), evaluated=int(len(ev_pos)), sampling=sampling, limit=args.limit,
                  truth_alleles=len(truth), pon_tagged_test=int((pon & test).sum()), results=results)
    if args.check:
        ref = json.loads((Path(args.check) / "report.json").read_text())["results"]
        bad = []
        for mode, r in results.items():
            x = ref[mode]
            for k1, k2 in (("somatic_ap", None), ("truth", "tp"), ("truth", "fp"), ("truth", "ap"), ("truth_after_pon", "f1"),
                           ("truth_after_pon", "ap"), ("germline", "ge_thr"), ("non_onref", "ge_thr"), ("calls_ge_thr", None)):
                v1, v2 = (r[k1], x[k1]) if k2 is None else (r[k1][k2], x[k1][k2])
                if abs(v1 - v2) > 1e-9:
                    bad.append((mode, k1, k2, v1, v2))
        print("check against ch6_counterfactual.py:", "identical" if not bad else bad)
        return
    (out / "report.json").write_text(json.dumps(report, indent=1) + "\n")
    cf.log("done")


if __name__ == "__main__":
    main()
