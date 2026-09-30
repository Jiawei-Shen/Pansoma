"""PoN variants on chr1 SNV with the my_check.py scorer (697 SNP truths, somatic BED ∩ germline BED, strict alleles):
for the batch-1024 sqrt / 1:200:10 runs at their validation recall-0.9 threshold and v1 e053, how many truths each PoN
rule removes, the P / R / F1 at that threshold and the max F1 over thresholds (chr1-tuned, reference only).
Per-allele lookups (gnomAD AF, dbSNP SAO, 1000G position / allele, CoLoRSdb position / allele AF) are cached in
synthesis/pon_rich.pkl. Read-only on inputs."""
import contextlib, io, os, pickle, sys
import pysam
sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/synthesis")
with contextlib.redirect_stdout(io.StringIO()):
    import b1024_check as B
M = B.M
CACHE = os.path.join(M.OUT, "pon_rich.pkl")
FILES = ["af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz", "1000g_pon.hg38.vcf.gz",
         "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"]


def lookup(vf, p0, ref, alt):
    """dict: gnomad_af (None = allele absent), dbsnp_sao (None = allele absent), kg_pos, kg_allele, colors_pos,
    colors_af (None = allele absent)."""
    out = dict(gnomad_af=None, dbsnp_sao=None, kg_pos=False, kg_allele=False, colors_pos=False, colors_af=None)
    for name, v in zip(("gnomad", "dbsnp", "kg", "colors"), vf):
        for r in v.fetch("chr1", p0, p0 + 1):
            if r.pos != p0 + 1:
                continue
            alts = [x.upper() for x in (r.alts or ())]
            hit = r.ref.upper() == ref and alt in alts
            if name == "gnomad" and hit:
                af = r.info["AF"][alts.index(alt)]
                out["gnomad_af"] = max(out["gnomad_af"] or 0.0, af)
            elif name == "dbsnp" and hit:
                out["dbsnp_sao"] = r.info.get("SAO", 0)
            elif name == "kg":
                out["kg_pos"] = True; out["kg_allele"] |= hit
            elif name == "colors":
                out["colors_pos"] = True
                if hit:
                    out["colors_af"] = max(out["colors_af"] or 0.0, r.info["AF"][alts.index(alt)])
    return out


runs = {n: B.data[r][0] for n, r in (("sqrt", B.RUNS[0]), ("1:200:10", B.RUNS[1]))}
thr = {n: B.data[r][2]["val_recall_0.9"]["threshold"] for n, r in (("sqrt", B.RUNS[0]), ("1:200:10", B.RUNS[1]))}
v1 = M.V1["v1_e053"]
keys = set(M.TRUTH) | {k for s in list(runs.values()) + [v1] for k, v in s.items() if v >= 0.02 and (k in M.TRUTH or M.inconf(k[0]))}
cache = pickle.load(open(CACHE, "rb")) if os.path.exists(CACHE) else {}
need = sorted(keys - set(cache))
if need:
    vf = [pysam.VariantFile(M.PON + f) for f in FILES]
    for k in need:
        cache[k] = lookup(vf, *k)
    pickle.dump(cache, open(CACHE, "wb"))

RULES = {
    "none": lambda c: False,
    "repo: gnomAD allele + dbSNP allele + 1000G pos + CoLoRSdb pos": lambda c: c["gnomad_af"] is not None or
        c["dbsnp_sao"] is not None or c["kg_pos"] or c["colors_pos"],
    "all four, allele match": lambda c: c["gnomad_af"] is not None or c["dbsnp_sao"] is not None or c["kg_allele"] or
        c["colors_af"] is not None,
    "repo without dbSNP": lambda c: c["gnomad_af"] is not None or c["kg_pos"] or c["colors_pos"],
    "repo, dbSNP non-somatic only (SAO != 2)": lambda c: c["gnomad_af"] is not None or c["dbsnp_sao"] not in (None, 2) or
        c["kg_pos"] or c["colors_pos"],
    "AF>=0.001 (gnomAD, CoLoRSdb) + dbSNP non-somatic + 1000G, allele match": lambda c: (c["gnomad_af"] or 0) >= 1e-3 or
        (c["colors_af"] or 0) >= 1e-3 or c["dbsnp_sao"] not in (None, 2) or c["kg_allele"],
    "AF>=0.001 (gnomAD, CoLoRSdb) + 1000G, allele match": lambda c: (c["gnomad_af"] or 0) >= 1e-3 or
        (c["colors_af"] or 0) >= 1e-3 or c["kg_allele"],
    "gnomAD allele (any AF)": lambda c: c["gnomad_af"] is not None,
    "gnomAD AF>=0.001": lambda c: (c["gnomad_af"] or 0) >= 1e-3,
    "gnomAD + CoLoRSdb allele (any AF)": lambda c: c["gnomad_af"] is not None or c["colors_af"] is not None,
}
tids = set(M.TRUTH.values())
print(f"lookups {len(cache)} (new {len(need)}); truths {len(tids)}; thresholds {thr}")
print("PoN rule | truths tagged | sqrt @ val R0.9: TP FP P / R / F1 | 1:200:10 @ val R0.9: P / R / F1 | max F1 (chr1-tuned): sqrt, 1:200:10, v1 e053")
for name, rule in RULES.items():
    pon = {k: rule(c) for k, c in cache.items()}
    tagged = len({t for k, t in M.TRUTH.items() if pon.get(k, False)})
    cols = []
    for n, s in runs.items():
        P, R, F, tp, fp = B.at(s, {k for k, v in s.items() if v >= thr[n]}, pon)
        cols.append((P, R, F, tp, fp))
    mx = [M.summ(M.curve(s, pon=pon))["maxF1"][0] for s in (runs["sqrt"], runs["1:200:10"], v1)]
    a, b = cols
    print(f"{name} | {tagged} | {a[3]} {a[4]} {a[0]:.3f} / {a[1]:.3f} / {a[2]:.3f} | {b[0]:.3f} / {b[1]:.3f} / {b[2]:.3f} | "
          f"{mx[0]:.3f}, {mx[1]:.3f}, {mx[2]:.3f}")
