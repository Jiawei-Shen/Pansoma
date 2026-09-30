"""All finished HG008T SNV runs on chr1 with the my_check.py scorer (697 SNP truths, somatic BED ∩ germline BED, strict
alleles) and the PoN rule at a gnomAD / CoLoRSdb floor of 1e-4 (allele match; dbSNP SAO != 2; 1000G), beside the
0.001 floor: the three decision rules of predict and the max F1 over thresholds (chr1-tuned, reference only).
Read-only on inputs; adds lookups to pon_rich.pkl."""
import contextlib, gzip, io, json, os, pickle, sys
import pysam
sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/synthesis")
with contextlib.redirect_stdout(io.StringIO()):
    import pon_variants as V
B, M = V.B, V.M


def rule(floor):
    return lambda c: ((c["gnomad_af"] or 0) >= floor or (c["colors_af"] or 0) >= floor
                      or c["dbsnp_sao"] not in (None, 2) or c["kg_allele"])


RUNS = [("HG008_Illumina_SNV_allnon_sqrt_b1024", "HG008T_Illumina"), ("HG008_Illumina_SNV_allnon_w1_200_10_b1024", "HG008T_Illumina"),
        ("HG008_Illumina_SNV_nopartial_af008", "HG008T_Illumina"), ("HG008_Illumina_SNV_nopartial_af008_nopc", "HG008T_Illumina"),
        ("HG008_PacBio_SNV_allnon_sqrt_b1024", "HG008T_PacBio"), ("HG008_ONT_SNV_allnon_sqrt_b1024", "HG008T_ONT")]


def read(run, sample):
    d = M.RUNS + run + "/test_chr1/"
    pf = [x for x in os.listdir(d) if x.endswith("SNV.predictions.ndjson.gz")][0]
    s, amax = {}, set()
    labels = f"/scratch/jshen/data/pansoma_v2_tensors/{sample}/tensors/SNV/chr1_labels.ndjson"
    with gzip.open(d + pf, "rt") as fp, open(labels) as fl:
        for lp, ll in zip(fp, fl):
            p, lab = json.loads(lp), json.loads(ll)
            assert p["candidate_id"] == lab["candidate_id"]
            if not p["in_test"] or lab["grch38"] is None:
                continue
            g = lab["grch38"]; k = (g["pos0"], g["ref"].upper(), g["alt"].upper())
            s[k] = max(s.get(k, -1), p["p_somatic"])
            if p["p_somatic"] > p["p_non"] and p["p_somatic"] >= p["p_germline"]:  # argmax, ties to the lower class
                amax.add(k)
    return s, amax, json.load(open(d + pf.replace("predictions.ndjson.gz", "metrics.json")))["rules"]


data = {r: read(r, smp) for r, smp in RUNS}
need = sorted({k for s, _, _ in data.values() for k, v in s.items() if v >= 0.02 and (k in M.TRUTH or M.inconf(k[0]))}
              - set(V.cache))
if need:
    vf = [pysam.VariantFile(M.PON + f) for f in V.FILES]
    for k in need:
        V.cache[k] = V.lookup(vf, *k)
    pickle.dump(V.cache, open(V.CACHE, "wb"))
print(f"new lookups {len(need)}; truths {M.NT}")
for floor in (1e-4, 1e-3):
    pon = {k: rule(floor)(c) for k, c in V.cache.items()}
    tagged = len({t for k, t in M.TRUTH.items() if pon.get(k, False)})
    print(f"\n== PoN floor {floor:g}: truths tagged {tagged}")
    print("run | rule (threshold) | no PoN: P / R / F1 (TP, FP) | PoN: P / R / F1 (TP, FP)")
    for run, (s, amax, rules) in data.items():
        name = run.replace("HG008_", "").replace("_SNV", "")
        for r in ("threshold_t", "argmax", "val_recall_0.9"):
            thr = rules[r]["threshold"]
            called = amax if r == "argmax" else {k for k, v in s.items() if v >= thr}
            a, b = B.at(s, called), B.at(s, called, pon)
            print(f"{name:<32} | {r:<14} ({'-' if thr is None else f'{thr:.3f}'}) | {a[0]:.3f} / {a[1]:.3f} / {a[2]:.3f} "
                  f"({a[3]}, {a[4]}) | {b[0]:.3f} / {b[1]:.3f} / {b[2]:.3f} ({b[3]}, {b[4]})")
        m1 = M.summ(M.curve(s, pon=pon))["maxF1"]
        print(f"{name:<32} | max F1 over thresholds with the PoN: {m1[0]:.3f} (P {m1[1]:.3f} R {m1[2]:.3f} @ {m1[3]:.3f})")
    m1 = M.summ(M.curve(M.V1["v1_e053"], pon=pon))["maxF1"]
    print(f"{'v1 e053':<32} | max F1 over thresholds with the PoN: {m1[0]:.3f} (P {m1[1]:.3f} R {m1[2]:.3f})")
