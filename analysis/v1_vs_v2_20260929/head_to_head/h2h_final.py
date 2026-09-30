"""Final head-to-head, v1 vs v2, HG008 Illumina SNV chr1, against the same 697 truths (PASS SNP, in somatic BED;
somatic.recall.tsv) and the same region (somatic BED ∩ germline BED). Read-only on all inputs.

Units: TP = distinct truth records with a called tensor (max score per truth). FP = distinct called (pos0, alt) in the
region that are not one of the 697 alleles; sites that are v2 label-1 tensors (partial INDEL matches, other-kind
somatic) are neutral for both models. Off-reference false calls (no GRCh38 position) are reported separately.

Rows:
  v2 <run>                exact, all 256,296 test tensors (p_somatic)
  v1 e067 (published)     = pretrained HG008_WGS_SNV.pth (byte-identical to model_e067_f1_0.1784.pth); stratified
                            random sample of the v1 chr1 GRCh38 candidates in the region (all 670 true, G and N strata
                            weighted by N_stratum / n_sampled)
  v1 e053 (VCF)           exact, the logged GPU run on all v1 chr1 test tensors, but the VCF keeps only PROB >= 0.5
Each with no PoN and with the 4 panels of normals (rule_repo = scripts/filter_panel_of_normals.py matching; rule_af001).
"""
import glob, gzip, json, pickle, sys
from collections import Counter, defaultdict
import numpy as np

HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from common import V2, truth_697, confident_chr1  # noqa
from curve import summary  # noqa
from pon_annot import annotate, rule_repo, rule_af001  # noqa

RUNS = ["HG008_Illumina_SNV_base", "HG008_Illumina_SNV_scalars", "HG008_Illumina_SNV_keephard_w100", "HG008_Illumina_SNV_allnon_w100"]
PRED = "/scratch/jshen/data/pansoma_net_v2_runs/{}/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz"
R1 = "/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB"
VCF_E053 = R1 + "/AF_HPRC/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.{}.vcf.gz"
VCF_N = R1 + "/AF_HPRC_HG008N_added/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.{}.vcf.gz"

truth = truth_697(); NT = len(truth)
tkey = {(t["pos0"], t["alt"]): tid for tid, t in truth.items()}
tref = {(t["pos0"], t["alt"]): t["ref"] for t in truth.values()}
conf = confident_chr1()
germ = set()
for line in open(f"{HERE}/germ_chr1_snp.tsv"):
    p = line.split("\t"); germ.add((int(p[0]), p[2]))

# ---------------- v2 per-tensor ----------------
lin, lab = [], []
for line in open(f"{V2}/SNV/chr1_labels.ndjson"):
    r = json.loads(line); g = r.get("grch38")
    lin.append((g["pos0"], g["ref"], g["alt"]) if g else None); lab.append(r["label"])
lab = np.array(lab)
neutral = {(k[0], k[2]) for k, l in zip(lin, lab) if l == 1 and k is not None and (k[0], k[2]) not in tkey}
af2 = np.load(f"{HERE}/v2_af.npy")
v2 = {}
for run in RUNS:
    rows = []
    with gzip.open(PRED.format(run), "rt") as f:
        for i, line in enumerate(f):
            p = json.loads(line)
            rows.append((p["in_test"], p["test_label"], [t for t in p.get("truth_ids") or [] if t in truth], p["p_somatic"]))
    v2[run] = rows
print("v2 tensors", len(lin), "neutral (v2 label-1 not-697) sites", len(neutral))

# ---------------- v1 e067 sample ----------------
recs = pickle.load(open(f"{HERE}/v1_val_recs.pkl", "rb"))  # (pos0, alt, lab, inconf, isgerm, sidx, iws, af)
parts = sorted(glob.glob(f"{HERE}/v1_sample_part*.npz"))
S_rows, S_str, S_p = [], [], []
for pth in parts:
    z = np.load(pth); S_rows.append(z["rows"]); S_str.append(z["strata"]); S_p.append(z["p_true"])
    nG_all, nN_all = int(z["n_G_all"]), int(z["n_N_all"])
S_rows, S_str, S_p = np.concatenate(S_rows), np.concatenate(S_str), np.concatenate(S_p)
assert len(set(S_rows.tolist())) == len(S_rows)
wG, wN = nG_all / (S_str == "G").sum(), nN_all / (S_str == "N").sum()
print(f"v1 sample parts={len(parts)} T={(S_str == 'T').sum()} G={(S_str == 'G').sum()} (w {wG:.2f}) N={(S_str == 'N').sum()} (w {wN:.2f})")


def v1_ref(pos0, alt):
    return None  # filled from FASTA below


# reference base for v1 sites (needed for allele-level PoN matching)
import pysam
fa = pysam.FastaFile(json.load(open(f"{V2}/SNV/labels.manifest.json"))["fasta"])


def refbase(pos0):
    return fa.fetch("chr1", pos0, pos0 + 1).upper()


# ---------------- v1 VCFs ----------------
def read_vcf(path):
    out = {}
    for line in gzip.open(path, "rt"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        for a in f[4].split(","):
            k = (f[0], int(f[1]) - 1, f[3].upper(), a.upper())
            out[k] = max(out.get(k, 0.0), float(info["PROB"]))
    return out


vcf = {"e053_linear": read_vcf(VCF_E053.format("linear")), "e053_graph": read_vcf(VCF_E053.format("graph")),
       "e053_ownPoN": read_vcf(VCF_E053.format("linear.filtered_PoN")),
       "N_linear": read_vcf(VCF_N.format("linear")), "N_all": read_vcf(VCF_N.format("all")),
       "N_ownPoN": read_vcf(VCF_N.format("linear.filtered_PoN"))}
for k, v in vcf.items():
    print("VCF", k, "records", len(v), "min PROB", round(min(v.values()), 3) if v else None)

# ---------------- PoN annotation for every site any row can call ----------------
sites = set()
for run in RUNS:
    for (it, tl, tids, p), k in zip(v2[run], lin):
        if it and k is not None and (p >= 0.005 or tids or tl == 2):
            sites.add(k)
for r in S_rows:
    pos0, alt = recs[r][0], recs[r][1]
    sites.add((pos0, refbase(pos0), alt))
for name in ("e053_linear", "N_linear"):
    for (c, pos0, ref, alt) in vcf[name]:
        if c == "chr1":
            sites.add((pos0, ref, alt))
ann = annotate(sites)
print("PoN-annotated sites", len(sites))
RULES = {"noPoN": lambda k: False,
         "PoN_repo": lambda k: rule_repo(ann[k]) if k in ann else False,
         "PoN_af001": lambda k: rule_af001(ann[k]) if k in ann else False}


def weighted(tp, fp, fw):
    tp = np.asarray(tp, float); fp = np.asarray(fp, float); fw = np.asarray(fw, float)
    thr = np.unique(np.concatenate([tp, fp]))[::-1]
    t = len(tp) - np.searchsorted(np.sort(tp), thr, side="left")
    o = np.argsort(fp); cw = np.concatenate([[0], np.cumsum(fw[o])])
    f = cw[-1] - cw[np.searchsorted(fp[o], thr, side="left")]
    rec = t / NT; prec = np.where(t + f > 0, t / np.maximum(t + f, 1e-9), 1.0)
    f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(prec + rec, 1e-12), 0)
    return dict(thr=thr, tp=t, fp=np.round(f).astype(int), rec=rec, prec=prec, f1=f1)


res = {}
comp = {}
for rn, drop in RULES.items():
    # v2
    for run in RUNS:
        for af_min in (0.0, 0.08):
            tp, fp, off = {}, {}, {}
            fpk = {}
            for i, ((it, tl, tids, p), k) in enumerate(zip(v2[run], lin)):
                if not it or af2[i] < af_min:
                    continue
                if tids:
                    if k is not None and drop(k):
                        continue
                    for t in tids:
                        tp[t] = max(tp.get(t, 0.0), p)
                elif tl in (0, 2):
                    if k is None:
                        off[i] = p
                    elif (k[0], k[2]) not in tkey and not drop(k):
                        key = (k[0], k[2]); fp[key] = max(fp.get(key, 0.0), p)
                        fpk[key] = "G" if tl == 2 else fpk.get(key, "N")
            tag = f"v2 {run.replace('HG008_Illumina_SNV_', '')} AF>={af_min}"
            c = weighted(list(tp.values()), list(fp.values()), np.ones(len(fp)))
            res[f"{rn} | {tag} | GRCh38 FP"] = summary(c, NT, len(tp))
            c2 = weighted(list(tp.values()), list(fp.values()) + list(off.values()), np.ones(len(fp) + len(off)))
            res[f"{rn} | {tag} | + off-ref FP"] = summary(c2, NT, len(tp))
            if af_min == 0.0:
                comp[f"{rn} | {tag}"] = (tp, fp, fpk, off)
    # v1 e067 sample
    tp, fps, fws, fst = {}, [], [], []
    for r, s, p in zip(S_rows, S_str, S_p):
        pos0, alt = recs[r][0], recs[r][1]
        k = (pos0, refbase(pos0), alt)
        if s == "T":
            if (pos0, alt) in tkey and not drop(k):
                t = tkey[(pos0, alt)]; tp[t] = max(tp.get(t, 0.0), float(p))
            continue
        if (pos0, alt) in tkey or (pos0, alt) in neutral or drop(k):
            continue
        fps.append(float(p)); fws.append(wG if s == "G" else wN); fst.append(s)
    c = weighted(list(tp.values()), fps, fws)
    res[f"{rn} | v1 e067 published (sample est.) | GRCh38 FP"] = summary(c, NT, len(tp))
    fst = np.array(fst); fps = np.array(fps); fws = np.array(fws)
    res[f"{rn} | v1 e067 published (sample est.) | FP non-germline only"] = summary(weighted(list(tp.values()), fps[fst == "N"], fws[fst == "N"]), NT, len(tp))
    comp[f"{rn} | v1 e067"] = (tp, fps, fws, fst)
    # v1 e053 VCF (PROB >= 0.5 only)
    for name in ("e053_linear", "N_linear"):
        tp, fp, cls = {}, {}, Counter()
        for (ch, pos0, ref, alt), p in vcf[name].items():
            if ch != "chr1":
                continue
            k = (pos0, ref, alt)
            if drop(k):
                continue
            if (pos0, alt) in tkey:
                t = tkey[(pos0, alt)]; tp[t] = max(tp.get(t, 0.0), p)
            elif conf.contains0(pos0):
                if (pos0, alt) in neutral:
                    cls["neutral"] += 1; continue
                fp[(pos0, alt)] = p; cls["germline" if (pos0, alt) in germ else "non"] += 1
            else:
                cls["outside_region"] += 1
        lab_ = "v1 e053 (VCF, PROB>=0.5)" if name == "e053_linear" else "v1 HG008N-added graph e064 (VCF, PROB>=0.5)"
        res[f"{rn} | {lab_} | GRCh38 FP"] = summary(weighted(list(tp.values()), list(fp.values()), np.ones(len(fp))), NT, len(tp))
        res[f"{rn} | {lab_} | GRCh38 FP"]["fp_classes_all_thr"] = dict(cls)
        res[f"{rn} | {lab_} | GRCh38 FP"]["at_0.5"] = dict(TP=len(tp), FP=len(fp), P=round(len(tp) / max(1, len(tp) + len(fp)), 4), R=round(len(tp) / NT, 4))

# own-PoN VCFs of the v1 pipeline (as produced), same scoring
for name, lab_ in (("e053_ownPoN", "v1 e053 own PoN VCF"), ("N_ownPoN", "v1 HG008N-added own PoN VCF")):
    tp, fp = {}, {}
    for (ch, pos0, ref, alt), p in vcf[name].items():
        if (pos0, alt) in tkey:
            t = tkey[(pos0, alt)]; tp[t] = max(tp.get(t, 0.0), p)
        elif conf.contains0(pos0) and (pos0, alt) not in neutral:
            fp[(pos0, alt)] = p
    res[f"ownPoN | {lab_} | GRCh38 FP"] = summary(weighted(list(tp.values()), list(fp.values()), np.ones(len(fp))), NT, len(tp))
    res[f"ownPoN | {lab_} | GRCh38 FP"]["at_0.5"] = dict(TP=len(tp), FP=len(fp), P=round(len(tp) / max(1, len(tp) + len(fp)), 4), R=round(len(tp) / NT, 4))

# does rule_repo / rule_af001 reproduce the v1 pipeline's own PoN removal?
lin53 = {k for k in vcf["e053_linear"] if k[0] == "chr1"}
kept_own = {k for k in vcf["e053_ownPoN"]}
agree = Counter()
for k in lin53:
    s = (k[1], k[2], k[3]); own = k in kept_own
    agree[("own_kept" if own else "own_removed", "repo_kept" if not rule_repo(ann[s]) else "repo_removed",
           "af001_kept" if not rule_af001(ann[s]) else "af001_removed")] += 1
res["pon_rule_vs_v1_own_filter(e053 linear calls)"] = {str(k): v for k, v in agree.items()}
print("PoN agreement:", dict(agree))

# off-reference false calls
offv1 = [p for (ch, *_), p in vcf["e053_graph"].items() if p >= 0.5]
res["offref"] = {"v1_e053_graph_calls_ge0.5": len(offv1), "v1_offref_candidates_chr1": 16018}
for run in RUNS:
    offp = np.array([p for (it, tl, tids, p), k in zip(v2[run], lin) if it and k is None and tl in (0, 2)])
    res["offref"][run] = dict(n=len(offp), ge0_5=int((offp >= 0.5).sum()))

# FP composition at matched truth recall (no PoN): germline vs non
def comp_at(recall, key):
    tp, *rest = comp[key]
    ts = np.sort(np.array(list(tp.values())))[::-1]
    if len(ts) < int(np.ceil(recall * NT)):
        return None
    thr = ts[int(np.ceil(recall * NT)) - 1]
    if "v1 e067" in key:
        _, fps, fws, fst = comp[key]
        m = fps >= thr
        return dict(thr=round(float(thr), 4), germline=int(round(fws[m & (fst == "G")].sum())), non=int(round(fws[m & (fst == "N")].sum())))
    _, fp, fpk, off = comp[key]
    g = sum(1 for k, p in fp.items() if p >= thr and fpk[k] == "G"); n = sum(1 for k, p in fp.items() if p >= thr and fpk[k] == "N")
    return dict(thr=round(float(thr), 4), germline=g, non=n, offref=int(sum(p >= thr for p in off.values())))


res["fp_composition_noPoN"] = {f"{k} R={r}": comp_at(r, k) for k in comp if k.startswith("noPoN") for r in (0.5, 0.8, 0.87, 0.9)}
res["fp_composition_PoN_repo"] = {f"{k} R={r}": comp_at(r, k) for k in comp if k.startswith("PoN_repo") for r in (0.5, 0.8, 0.87)}
json.dump(res, open(f"{HERE}/h2h_final.json", "w"), indent=1, default=str)


def fmt(v):
    if v is None:
        return "  n/a"
    return f"{v[0]:.3f}"


print("\nrow | ceiling | R@P.05 R@P.10 R@P.15 R@P.20 | P@R.5 P@R.7 P@R.8 P@R.9 | maxF1 (P, R)")
for k, s in res.items():
    if not isinstance(s, dict) or "maxF1" not in s:
        continue
    m = s["maxF1"]
    print(f"{k} | {s['ceiling']:.3f} | " + " ".join(fmt(s[f'R@P{p:.2f}']) for p in (0.05, 0.10, 0.15, 0.20)) + " | " +
          " ".join(fmt(s[f'P@R{r:.1f}']) for r in (0.5, 0.7, 0.8, 0.9)) + f" | {m[0]:.3f} ({m[1]:.3f}, {m[2]:.3f})" +
          (f" | at0.5 {s['at_0.5']}" if "at_0.5" in s else ""))
for k in ("offref", "fp_composition_noPoN", "fp_composition_PoN_repo", "pon_rule_vs_v1_own_filter(e053 linear calls)"):
    print(k, json.dumps(res[k], indent=0))
