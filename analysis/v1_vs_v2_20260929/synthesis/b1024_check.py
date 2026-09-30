"""HG008_Illumina_SNV_allnon_{sqrt,w1_200_10}_b1024 on chr1 with the my_check.py scorer (697 SNP truths, somatic BED ∩
germline BED, strict alleles, repo-rule PoN): the three decision rules of predict (threshold t, argmax, the checkpoint's
validation recall-0.9 threshold), without and with the PoN, plus the max F1 over all thresholds. Read-only on inputs;
adds PoN lookups for new sites to synthesis/my_pon.pkl."""
import contextlib, gzip, io, json, os, pickle, sys
import pysam
sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/synthesis")
with contextlib.redirect_stdout(io.StringIO()):
    import my_check as M
RUNS = ["HG008_Illumina_SNV_allnon_sqrt_b1024", "HG008_Illumina_SNV_allnon_w1_200_10_b1024"]


def read(run):
    """Per GRCh38 allele: max p_somatic over its tensors and whether any of them is argmax somatic (the pred field is the
    threshold-t call)."""
    d = M.RUNS + run + "/test_chr1/"
    pf = [x for x in os.listdir(d) if x.endswith("SNV.predictions.ndjson.gz")][0]
    p_max, amax = {}, set()
    with gzip.open(d + pf, "rt") as fp, open(M.V2 + "/SNV/chr1_labels.ndjson") as fl:
        for lp, ll in zip(fp, fl):
            p, lab = json.loads(lp), json.loads(ll)
            assert p["candidate_id"] == lab["candidate_id"]
            if not p["in_test"] or lab["grch38"] is None:
                continue
            g = lab["grch38"]; k = (g["pos0"], g["ref"].upper(), g["alt"].upper())
            p_max[k] = max(p_max.get(k, -1), p["p_somatic"])
            if p["p_somatic"] > p["p_non"] and p["p_somatic"] >= p["p_germline"]:  # argmax, ties to the lower class index
                amax.add(k)
    rules = json.load(open(d + pf.replace("predictions.ndjson.gz", "metrics.json")))["rules"]
    return p_max, amax, rules


def at(sites, called, pon=None):
    """P / R / F1 of the called alleles: TP distinct truths, FP in-conf alleles not truth (PoN-flagged dropped)."""
    tids, fp = set(), 0
    for k in called:
        if pon is not None and pon.get(k, False):
            continue
        if k in M.TRUTH:
            tids.add(M.TRUTH[k])
        elif M.inconf(k[0]):
            fp += 1
    tp = len(tids); P = tp / (tp + fp) if tp + fp else 1.0; R = tp / M.NT
    return P, R, 2 * P * R / (P + R) if P + R else 0.0, tp, fp


data = {r: read(r) for r in RUNS}
need = sorted({k for s, _, _ in data.values() for k, v in s.items() if v >= 0.02 and (k in M.TRUTH or M.inconf(k[0]))}
              - set(M.cache))
if need:
    vf = [pysam.VariantFile(M.PON + f) for f in ["af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
                                                   "1000g_pon.hg38.vcf.gz", "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"]]
    for p0, ref, alt in need:
        hit = []
        for v, m in zip(vf, ["allele", "allele", "position", "position"]):
            hit.append(any(r.pos == p0 + 1 and (m == "position" or (r.ref.upper() == ref and alt in [x.upper() for x in (r.alts or ())]))
                           for r in v.fetch("chr1", p0, p0 + 1)))
        M.cache[(p0, ref, alt)] = tuple(hit)
    pickle.dump(M.cache, open(M.CACHE, "wb"))
PONF = {k: any(v) for k, v in M.cache.items()}
print(f"new PoN lookups: {len(need)}; truths {M.NT}")
print("run | rule (threshold) | no PoN: P / R / F1 (TP, FP) | repo PoN: P / R / F1 (TP, FP)")
for run, (s, amax, rules) in data.items():
    name = run.replace("HG008_Illumina_SNV_allnon_", "").replace("_b1024", "")
    for rule, called in [("threshold_t", None), ("argmax", amax), ("val_recall_0.9", None)]:
        thr = rules[rule]["threshold"]
        if called is None:
            called = {k for k, v in s.items() if v >= thr}
        a, b = at(s, called), at(s, called, PONF)
        print(f"{name:<11} | {rule:<14} ({'-' if thr is None else f'{thr:.3f}'}) | {a[0]:.3f} / {a[1]:.3f} / {a[2]:.3f} ({a[3]}, {a[4]})"
              f" | {b[0]:.3f} / {b[1]:.3f} / {b[2]:.3f} ({b[3]}, {b[4]})")
    m0, m1 = M.summ(M.curve(s)), M.summ(M.curve(s, pon=PONF))
    print(f"{name:<11} | max F1 over thresholds: no PoN {m0['maxF1'][0]:.3f} (P {m0['maxF1'][1]:.3f} R {m0['maxF1'][2]:.3f} @ "
          f"{m0['maxF1'][3]:.3f}) | repo PoN {m1['maxF1'][0]:.3f} (P {m1['maxF1'][1]:.3f} R {m1['maxF1'][2]:.3f} @ {m1['maxF1'][3]:.3f})")
