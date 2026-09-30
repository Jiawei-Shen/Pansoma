"""keephard_w100_nopartial (SNV partials left out of training/validation) vs its control keephard_w100 and v1, on chr1 with
the my_check.py scorer (697 SNP truths, somatic BED ∩ germline BED, strict alleles, repo-rule PoN). Read-only on inputs;
adds PoN lookups for new sites to synthesis/my_pon.pkl."""
import collections, contextlib, io, json, pickle, sys
import numpy as np
import pysam
sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/synthesis")
with contextlib.redirect_stdout(io.StringIO()):
    import my_check as M
RUN = "HG008_Illumina_SNV_keephard_w100_nopartial"
M.V2R[RUN] = M.read_v2(RUN)
s_new = M.V2R[RUN][0]
need = sorted({k for k, v in s_new.items() if v >= 0.02 and (k in M.TRUTH or M.inconf(k[0]))} - set(M.cache))
if need:
    vf = [pysam.VariantFile(M.PON + f) for f in ["af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
                                                   "1000g_pon.hg38.vcf.gz", "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"]]
    for p0, ref, alt in need:
        hit = []
        for v, m in zip(vf, ["allele", "allele", "position", "position"]):
            h = any(r.pos == p0 + 1 and (m == "position" or (r.ref.upper() == ref and alt in [x.upper() for x in (r.alts or ())]))
                    for r in v.fetch("chr1", p0, p0 + 1))
            hit.append(h)
        M.cache[(p0, ref, alt)] = tuple(hit)
    pickle.dump(M.cache, open(M.CACHE, "wb"))
PONF = {k: any(v) for k, v in M.cache.items()}
lab = {}
for line in open(M.V2 + "/SNV/chr1_labels.ndjson"):
    r = json.loads(line); g = r["grch38"]
    if g is not None:
        lab.setdefault((g["pos0"], g["ref"].upper(), g["alt"].upper()), r["reason"])
print(f"new PoN lookups: {len(need)}")
print("run | no PoN: P@TP608, max F1 (P, R) | repo PoN: P@TP561, max F1 (P, R) | false calls @TP561: germline / non-germline (near_truth_allele_mismatch)")
for name, sites in [("v1 e053", M.V1["v1_e053"])] + [(r.replace("HG008_Illumina_SNV_", "v2 "), M.V2R[r][0]) for r in
                     ("HG008_Illumina_SNV_base", "HG008_Illumina_SNV_keephard_w100", "HG008_Illumina_SNV_allnon_w100", RUN)]:
    c0, c1 = M.curve(sites), M.curve(sites, pon=PONF)
    m0, m1 = M.summ(c0), M.summ(c1)
    a0, a1 = M.at_tp(c0, 608), M.at_tp(c1, 561)
    p0 = f"{a0[1] / (a0[1] + a0[2]):.3f}" if a0 else "n.a."; p1 = f"{a1[1] / (a1[1] + a1[2]):.3f}" if a1 else "n.a."
    comp = collections.Counter()
    if a1:
        for k, v in sites.items():
            if v >= a1[0] and k not in M.TRUTH and M.inconf(k[0]) and not PONF.get(k, False):
                comp["germ" if k in M.GERM else ("near" if lab.get(k) == "near_truth_allele_mismatch" else "other")] += 1
    print(f"{name:<26} | {p0}, {m0['maxF1'][0]:.3f} ({m0['maxF1'][1]:.3f}, {m0['maxF1'][2]:.3f}) | {p1}, {m1['maxF1'][0]:.3f} "
          f"({m1['maxF1'][1]:.3f}, {m1['maxF1'][2]:.3f}) | {comp['germ']} / {comp['near'] + comp['other']} ({comp['near']})")
