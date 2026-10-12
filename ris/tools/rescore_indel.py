"""Call-level INDEL rescorer (ris experiment): one logistic model over label-free evidence per call, fitted on the
held-out tuning chromosome (chr2) of the training sets, one shared threshold, applied once to the test chromosome.

Features per call of a default-decision linear VCF (graph_vcf -> linear_vcf of exp_calls.sh, FILTER PASS/LowQual,
P_SOMATIC >= 0.01), all available for an unlabelled sample:
  model      logit P_SOMATIC, logit P_GERMLINE, logit P_NON
  population log10 max AF of the exact allele (gnomAD, CoLoRSdb); log10 max AF of any length-changing allele within
             [pos-6, pos+5] (population variability of the locus, Mutect2's STR idea)
  repeat     repeat tract length of the event's unit (FASTA, capped 40) / 10; homopolymer >= 6 flag
  cluster    other calls within 5 kb and 20 kb (log1p; Mutect2's clustered-events idea)
  alignment  LEFTSHIFT > 0 flag (the event moved when left-aligned on GRCh38)
Feature set v2 adds (from the tumor-only tool survey): the exact allele absent from both population databases
(DeepSomatic / Mutect2 treat absence specially); ClairS-TO's reference sequence entropy (k = 5 over the 33 bp around
POS, LowSeqEntropy < 0.9); Mutect2's tandem-repeat terms (repeat unit length, repeat units in the reference, polymerase
slippage eligibility: unit x units >= 8 and a one-unit change). fit --feature-set chooses the columns; a model records them.
The PoN script's rule stays a hard filter (exact allele AF >= floor tagged); the rescorer ranks the untagged calls.

    features RUN_DIR SET CHROM OUT.npz        features (+ labels from RUN_DIR/vcf_SET_CHROM/eval_bed if present)
    fit A.npz [B.npz ...] --out model.json    pooled + leave-one-sample-out fits, shared threshold, chr2 block CV
    apply --model model.json --variant V --features X.npz --vcf LIN.vcf.gz --out OUT.vcf.gz
"""
import argparse
import json
import math
import os
import sys

import numpy as np
import pysam

PON = "/storage2/fs1/epigenome/Active/jiawei.s/Dataset/panel_of_normal_VCFs"
FASTA = "/storage2/fs1/epigenome/Active/jiawei.s/Dataset/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
FEATURES = ["logit_p_somatic", "logit_p_germline", "logit_p_non", "log_pon_af", "log_near_af", "tract", "hp6",
            "dens5k", "dens20k", "leftshift"]
V2_EXTRA = ["pon_absent", "entropy", "low_entropy", "ru_len", "rpa_ref", "slip_elig"]
ALL_FEATURES = FEATURES + V2_EXTRA
FEATURE_SETS = {"v1": FEATURES, "v2": ALL_FEATURES}
MIN_SCORE = 0.01
FLOORS = (0.01, 0.02)
TOL = 1 - 1e-6  # the PoN script's float32 tolerance


def logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def info(field):
    out = {}
    for item in field.split(";"):
        k, _, v = item.partition("=")
        out[k] = v if _ else True
    return out


def pop_af(vcf, chrom, pos, ref, alt):
    """(max AF of the exact allele, max AF of any length-changing allele within [pos-6, pos+5])."""
    exact = near = 0.0
    contig = chrom if chrom in vcf.header.contigs else chrom.replace("chr", "")
    try:
        recs = vcf.fetch(contig, max(0, pos - 7), pos + 5)
    except ValueError:
        return 0.0, 0.0
    for c in recs:
        alts = c.alts or ()
        af = c.info.get("AF")
        af = af if isinstance(af, tuple) else (af,) * len(alts)
        for k, a in enumerate(alts):
            if af[k] is None or not a or a.startswith("<"):
                continue
            if len(a) != len(c.ref) and pos - 6 <= c.pos <= pos + 5:
                near = max(near, af[k])
            if c.pos == pos and c.ref.upper() == ref.upper() and a.upper() == alt.upper():
                exact = max(exact, af[k])
    return exact, near


def entropy(fa, chrom, pos, k=5, span=16):
    """ClairS-TO's calculate_sequence_entropy (src/postfilter_variants.py): the 2*span+1 = 33 bp of the reference
    around POS; one k-mer ends at each of the 33 bases (the first k-1 are left-padded with A, as its rolling 2-bit hash
    starting at 0 does); p = count / 33; Shannon entropy / ln 33. LowSeqEntropy is entropy < 0.9."""
    lo = max(0, pos - 1 - span)
    seq = fa.fetch(chrom, lo, pos + span).upper()
    n = len(seq)
    if n < 2:
        return 1.0
    padded = "A" * (k - 1) + "".join(c if c in "ACGT" else "A" for c in seq)
    counts = {}
    for i in range(n):
        kmer = padded[i:i + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    p = np.array(list(counts.values()), dtype=float) / n
    return float(-(p * np.log(p)).sum() / math.log(n))


def repeat(fa, chrom, pos, ref, alt):
    """(tract length in bases of the event's smallest repeat unit around the event, homopolymer >= 6, unit length,
    units changed by the event)."""
    ev = (alt[1:] if len(alt) > len(ref) else ref[1:]).upper()
    if not ev:
        return 0, 0, 0, 0
    unit = ev
    for k in range(1, len(ev) + 1):
        if len(ev) % k == 0 and ev[:k] * (len(ev) // k) == ev:
            unit = ev[:k]
            break
    lo = max(0, pos - 200)
    seq = fa.fetch(chrom, lo, pos + 200).upper()
    i = pos - lo  # 0-based index of the base after the padding base
    right = 0
    while seq[i + right:i + right + len(unit)] == unit:
        right += len(unit)
    left = 0
    while i - left - len(unit) >= 0 and seq[i - left - len(unit):i - left] == unit:
        left += len(unit)
    tract = left + right
    return tract, int(len(unit) == 1 and tract >= 6), len(unit), len(ev) // len(unit)


def load_labels(eval_dir):
    """{(chrom, pos, ref, alt): 1/0} from the raw curve of a vcfeval.py output; number of truth alleles."""
    d = f"{eval_dir}/vcfeval/raw_curve"
    lab = {}
    for name, y in (("tp", 1), ("fp", 0)):
        with pysam.VariantFile(f"{d}/{name}.vcf.gz") as v:
            for r in v:
                lab[(r.chrom, r.pos, r.ref, r.alts[0])] = y
    n_truth = 0
    for name in ("tp-baseline", "fn"):
        with pysam.VariantFile(f"{d}/{name}.vcf.gz") as v:
            n_truth += sum(1 for _ in v)
    return lab, n_truth


def features(a):
    import glob
    vdir = f"{a.run}/vcf_{a.set}_{a.chrom}"
    lin = [p for p in glob.glob(f"{vdir}/*.INDEL.linear.vcf.gz") if not p.endswith(".pon.vcf.gz")]
    if len(lin) != 1:
        raise SystemExit(f"{vdir}: expected one default-decision linear VCF, found {lin}")
    lin = lin[0]
    gn = pysam.VariantFile(f"{PON}/af-only-gnomad.hg38.vcf.gz")
    co = pysam.VariantFile(f"{PON}/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
    fa = pysam.FastaFile(FASTA)
    rows, keys, ids, kept = [], [], [], []
    with pysam.VariantFile(lin) as v:
        for idx, r in enumerate(v):
            filt = set(r.filter.keys())
            ps = r.info.get("P_SOMATIC")
            # scored: PASS or LowQual only (no LowAF or other filter, whatever the order), P_SOMATIC >= min score
            ok = ps is not None and ps >= MIN_SCORE and (not filt or filt <= {"PASS", "LowQual"})
            kept.append(ok)
            if not ok:
                continue
            alt = r.alts[0]
            e1, n1 = pop_af(gn, r.chrom, r.pos, r.ref, alt)
            e2, n2 = pop_af(co, r.chrom, r.pos, r.ref, alt)
            tract, hp6, ulen, dunits = repeat(fa, r.chrom, r.pos, r.ref, alt)
            rows.append([ps, r.info.get("P_GERMLINE", 0.0), r.info.get("P_NON", 0.0), max(e1, e2), max(n1, n2),
                         tract, hp6, r.pos, int(r.info.get("LEFTSHIFT", 0) or 0), entropy(fa, r.chrom, r.pos),
                         ulen, dunits])
            keys.append((r.chrom, r.pos, r.ref, alt))
            ids.append(r.id or "")
    t = np.array(rows, dtype=float).reshape(-1, 12)
    pos = t[:, 7]
    order = np.argsort(pos)
    sp = pos[order]
    dens = {}
    for w in (5000, 20000):
        c = np.searchsorted(sp, sp + w, side="right") - np.searchsorted(sp, sp - w, side="left") - 1
        dd = np.empty_like(c)
        dd[order] = c
        dens[w] = dd
    ulen = np.maximum(t[:, 10], 1)
    rpa = t[:, 5] / ulen  # repeat units in the reference tract
    X = np.column_stack([logit(t[:, 0]), logit(t[:, 1]), logit(t[:, 2]), np.log10(t[:, 3] + 1e-5),
                         np.log10(t[:, 4] + 1e-5), np.minimum(t[:, 5], 40) / 10, t[:, 6], np.log1p(dens[5000]),
                         np.log1p(dens[20000]), (t[:, 8] > 0).astype(float),
                         # v2
                         (t[:, 3] == 0).astype(float), t[:, 9], (t[:, 9] < 0.9).astype(float),
                         np.minimum(t[:, 10], 6) / 6, np.log1p(rpa), ((t[:, 5] >= 8) & (t[:, 11] == 1)).astype(float)])
    y = np.full(len(X), -1, dtype=np.int8)
    n_truth = -1
    eval_dir = f"{vdir}/eval_bed"
    if os.path.exists(f"{eval_dir}/vcfeval/raw_curve/tp.vcf.gz"):
        lab, n_truth = load_labels(eval_dir)
        y = np.array([lab.get(k, -1) for k in keys], dtype=np.int8)
    np.savez_compressed(a.out, X=X, y=y, pon_af=t[:, 3], p_somatic=t[:, 0], pos=pos, n_truth=n_truth,
                        kept=np.array(kept), ids=np.array(ids), set=a.set, chrom=a.chrom, run=a.run, vcf=lin,
                        features=np.array(ALL_FEATURES))
    print(f"{a.set} {a.chrom}: {len(X):,} calls (of {len(kept):,} records), labelled {(y >= 0).sum():,} "
          f"(TP {(y == 1).sum():,}), truth {n_truth}", flush=True)


# ---------------------------------------------------------------- fitting

def fit_logistic(X, y, w, l2=1.0, iters=50):
    """Weighted L2 logistic regression (IRLS); intercept unpenalized. Returns coef with intercept first."""
    Xb = np.column_stack([np.ones(len(X)), X])
    beta = np.zeros(Xb.shape[1])
    pen = np.full(Xb.shape[1], l2)
    pen[0] = 0.0
    for _ in range(iters):
        z = Xb @ beta
        p = 1 / (1 + np.exp(-z))
        g = Xb.T @ (w * (p - y)) + pen * beta
        H = (Xb * (w * p * (1 - p))[:, None]).T @ Xb + np.diag(pen)
        step = np.linalg.solve(H, g)
        beta -= step
        if np.max(np.abs(step)) < 1e-8:
            break
    return beta


def score(model, X):
    X = X[:, model["cols"]]
    Z = (X - np.array(model["mean"])) / np.array(model["std"])
    return 1 / (1 + np.exp(-(model["coef"][0] + Z @ np.array(model["coef"][1:]))))


def set_f1(s, y, n_truth, thr):
    keep = s >= thr
    tp = float((y[keep] == 1).sum())
    return 2 * tp / (keep.sum() + n_truth) if n_truth > 0 else 0.0


def rows(d, floor):
    """Fit rows of one set: labelled, untagged by the PoN at this floor."""
    return (d["y"] >= 0) & ~(d["pon_af"] >= floor * TOL)


def fit_variant(data, floor, l2=1.0, idx=None, cols=None):
    """Pooled fit over the given sets (each set weight 1000 in total); shared threshold maximizing the mean chr2
    F1 over those sets."""
    Xs, ys, ws = [], [], []
    for k, d in enumerate(data):
        m = rows(d, floor) if idx is None else rows(d, floor) & idx[k]
        Xs.append(d["X"][m][:, cols])
        ys.append(d["y"][m].astype(float))
        ws.append(np.full(m.sum(), 1000.0 / max(m.sum(), 1)))
    X, y, w = np.vstack(Xs), np.concatenate(ys), np.concatenate(ws)
    mean, std = X.mean(0), X.std(0)
    std[std == 0] = 1
    coef = fit_logistic((X - mean) / std, y, w, l2)
    model = dict(mean=mean.tolist(), std=std.tolist(), coef=coef.tolist(), floor=floor, cols=list(cols))
    return model


def best_threshold(model, data, floor, idx=None, n_truth_frac=None):
    grid = 1 / (1 + np.exp(-np.linspace(-6, 2, 321)))
    best = (-1, None)
    for thr in grid:
        f1s = []
        for k, d in enumerate(data):
            m = rows(d, floor) if idx is None else rows(d, floor) & idx[k]
            nt = d["n_truth"] if n_truth_frac is None else d["n_truth"] * n_truth_frac[k]
            f1s.append(set_f1(score(model, d["X"][m]), d["y"][m], nt, thr))
        mf = float(np.mean(f1s))
        if mf > best[0]:
            best = (mf, float(thr), f1s)
    return best


def p_only(data, floor, idx=None, n_truth_frac=None):
    """Baseline: a shared threshold on P_SOMATIC alone (same PoN floor), chosen the same way."""
    grid = np.unique(np.round(np.concatenate([d["p_somatic"] for d in data]), 3))
    best = (-1, None, None)
    for thr in grid:
        f1s = []
        for k, d in enumerate(data):
            m = rows(d, floor) if idx is None else rows(d, floor) & idx[k]
            nt = d["n_truth"] if n_truth_frac is None else d["n_truth"] * n_truth_frac[k]
            f1s.append(set_f1(d["p_somatic"][m], d["y"][m], nt, thr))
        mf = float(np.mean(f1s))
        if mf > best[0]:
            best = (mf, float(thr), f1s)
    return best


def fit(a):
    data = []
    for path in a.npz:
        z = np.load(path, allow_pickle=False)
        d = {k: z[k] for k in z.files}
        d["n_truth"] = int(d["n_truth"])
        d["set"] = str(d["set"])
        if d["n_truth"] <= 0:
            raise SystemExit(f"{path}: no labels (run on the tuning chromosome)")
        if str(d["chrom"]) == "chr1":
            raise SystemExit(f"{path}: chr1 is the test chromosome; fit only on the tuning chromosome")
        data.append(d)
    names = [d["set"] for d in data]
    sample = lambda s: s.split("_")[0]  # noqa: E731
    feats = FEATURE_SETS[a.feature_set]
    stored = [str(x) for x in data[0]["features"]]
    cols = [stored.index(f) for f in feats]
    floors = tuple(a.floors)
    out = dict(features=feats, feature_set=a.feature_set, floors=floors, sets=names, npz=a.npz, l2=a.l2,
               variants={}, note=(
        "Primary: 'pooled' (all sets' chr2, shared threshold and PoN floor chosen on chr2). "
        "'loso_<SAMPLE>' scores SAMPLE's sets with a model fitted only on the other sample's chr2: the estimate for "
        "an unlabelled sample. Chosen before any test-chromosome evaluation."))
    groups = {"pooled": (list(range(len(data))), list(range(len(data))))}
    for s in sorted({sample(n) for n in names}):
        train = [k for k, n in enumerate(names) if sample(n) != s]
        apply_to = [k for k, n in enumerate(names) if sample(n) == s]
        if train and apply_to:
            groups[f"loso_{s}"] = (train, apply_to)
    for name, (train, apply_to) in groups.items():
        tdata = [data[k] for k in train]
        best = None
        for floor in floors:
            model = fit_variant(tdata, floor, a.l2, cols=cols)
            mf, thr, f1s = best_threshold(model, tdata, floor)
            if best is None or mf > best[0]:
                best = (mf, thr, f1s, model, floor)
        mf, thr, f1s, model, floor = best
        model.update(threshold=thr, train_sets=[names[k] for k in train], apply_sets=[names[k] for k in apply_to],
                     chr2_train_mean_f1=mf, chr2_train_f1=dict(zip([names[k] for k in train], f1s)))
        # chr2 F1 of the sets it applies to (for loso: out of sample)
        model["chr2_apply_f1"] = {names[k]: set_f1(score(model, data[k]["X"][rows(data[k], floor)]),
                                                   data[k]["y"][rows(data[k], floor)], data[k]["n_truth"], thr)
                                  for k in apply_to}
        pb = p_only(tdata, floor)
        model["p_only"] = dict(threshold=pb[1], chr2_train_mean_f1=pb[0],
                               chr2_apply_f1={names[k]: set_f1(data[k]["p_somatic"][rows(data[k], floor)],
                                                               data[k]["y"][rows(data[k], floor)],
                                                               data[k]["n_truth"], pb[1]) for k in apply_to})
        out["variants"][name] = model
        print(f"[{name}] floor {floor} threshold {thr:.4f} train mean F1 {mf:.3f} | apply F1 "
              + " ".join(f"{s} {v:.3f}" for s, v in model["chr2_apply_f1"].items())
              + f" | p-only apply " + " ".join(f"{s} {v:.3f}" for s, v in model["p_only"]["chr2_apply_f1"].items()),
              flush=True)
        print("   coef " + " ".join(f"{f}={c:+.2f}" for f, c in zip(["icpt"] + feats, model["coef"])), flush=True)
    # chr2 2-fold block CV of the pooled design: fit on half the 5 Mb blocks, evaluate on the other half
    rng = np.random.default_rng(0)
    cv = {n: [] for n in names}
    cvp = {n: [] for n in names}
    for rep in range(a.cv_repeats):
        blocks = [np.floor(d["pos"] / 5e6).astype(int) for d in data]
        nb = max(int(b.max()) for b in blocks) + 1
        half = rng.permutation(nb) < nb // 2
        for fold in (half, ~half):
            tr = [fold[b] for b in blocks]
            te = [~fold[b] for b in blocks]
            frac_tr = [float(np.mean(t[d["y"] == 1])) if (d["y"] == 1).any() else 0.5 for t, d in zip(tr, data)]
            frac_te = [1 - f for f in frac_tr]
            floor = out["variants"]["pooled"]["floor"]
            model = fit_variant(data, floor, a.l2, idx=tr, cols=cols)
            _, thr, _ = best_threshold(model, data, floor, idx=tr, n_truth_frac=frac_tr)
            _, pthr, _ = p_only(data, floor, idx=tr, n_truth_frac=frac_tr)
            for k, d in enumerate(data):
                m = rows(d, floor) & te[k]
                nt = d["n_truth"] * frac_te[k]
                cv[names[k]].append(set_f1(score(model, d["X"][m]), d["y"][m], nt, thr))
                cvp[names[k]].append(set_f1(d["p_somatic"][m], d["y"][m], nt, pthr))
    out["chr2_block_cv"] = {n: dict(rescorer=float(np.mean(cv[n])), p_only=float(np.mean(cvp[n])),
                                    rescorer_wins=float(np.mean(np.array(cv[n]) > np.array(cvp[n]))))
                            for n in names}
    for n, v in out["chr2_block_cv"].items():
        print(f"   CV {n}: p-only {v['p_only']:.3f} -> rescorer {v['rescorer']:.3f} (wins {v['rescorer_wins']:.0%})")
    json.dump(out, open(a.out, "w"), indent=1)


# ---------------------------------------------------------------- applying

def apply(a):
    import gzip
    import subprocess
    model_all = json.load(open(a.model))
    model = model_all["variants"][a.variant]
    z = np.load(a.features, allow_pickle=False)
    if os.path.realpath(str(z["vcf"])) != os.path.realpath(a.vcf):
        raise SystemExit(f"{a.features} was computed from {z['vcf']}, not {a.vcf}")
    s = score(model, z["X"])
    kept = z["kept"]
    out_score = np.zeros(len(kept))
    out_score[np.flatnonzero(kept)] = s
    thr = model["threshold"]
    plain = a.out[:-3]
    with gzip.open(a.vcf, "rt") as f, open(plain, "w") as o:
        i = 0
        for line in f:
            if line.startswith("##pansoma_threshold="):
                o.write(f"##pansoma_threshold={thr}\n")
                continue
            if line.startswith("#CHROM"):
                o.write('##FILTER=<ID=NoScore,Description="Not scored by the rescorer (P_SOMATIC < 0.01 or a filter other '
                        'than PASS/LowQual): left out of the call set and of the vcfeval curve">\n')
                o.write('##INFO=<ID=P_RESCORE,Number=1,Type=Float,Description="ris rescore_indel.py probability">\n')
                o.write(f"##ris_rescore={json.dumps(dict(model=os.path.abspath(a.model), variant=a.variant, threshold=thr, floor=model['floor']))}\n")
                o.write(line)
                continue
            if line.startswith("#"):
                o.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            p = float(out_score[i])
            if kept[i]:
                c[6] = "PASS" if p >= thr else "LowQual"
            elif c[6] in ("PASS", "LowQual", "."):
                c[6] = "NoScore"
            c[5] = f"{min(50.0, -10 * math.log10(max(1e-12, 1 - p))):.3f}" if kept[i] else "0"
            c[7] = (c[7] + ";" if c[7] != "." else "") + f"P_RESCORE={p:.5f}"
            o.write("\t".join(c) + "\n")
            i += 1
    if i != len(kept):
        raise SystemExit(f"{a.vcf}: {i} records, features have {len(kept)}")
    for old in (a.out, a.out + ".tbi"):
        if os.path.exists(old):
            os.remove(old)
    pysam.tabix_index(plain, preset="vcf", force=True)
    print(f"{a.out}: {int(kept.sum()):,} scored, PASS {int((kept & (out_score >= thr)).sum()):,} at {thr:.4f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("features")
    f.add_argument("run")
    f.add_argument("set")
    f.add_argument("chrom")
    f.add_argument("out")
    g = sub.add_parser("fit")
    g.add_argument("npz", nargs="+")
    g.add_argument("--out", required=True)
    g.add_argument("--l2", type=float, default=1.0)
    g.add_argument("--cv-repeats", type=int, default=20)
    g.add_argument("--feature-set", choices=sorted(FEATURE_SETS), default="v1")
    g.add_argument("--floors", nargs="+", type=float, default=list(FLOORS))
    h = sub.add_parser("apply")
    h.add_argument("--model", required=True)
    h.add_argument("--variant", required=True)
    h.add_argument("--features", required=True)
    h.add_argument("--vcf", required=True)
    h.add_argument("--out", required=True)
    a = ap.parse_args()
    {"features": features, "fit": fit, "apply": apply}[a.cmd](a)


if __name__ == "__main__":
    main()
