"""Raw PoN facts for a set of chr1 SNV sites (pos0, ref, alt); cached. Read-only on the PoN VCFs."""
import os, pickle
import pysam
PON = "/scratch/jshen/data/Panels_of_Normal"
FILES = dict(gnomad="af-only-gnomad.hg38.vcf.gz", dbsnp="Homo_sapiens_assembly38.dbsnp138.vcf.gz",
             pon1000g="1000g_pon.hg38.vcf.gz", colors="CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
CACHE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head/pon_sites.pkl"


def annotate(sites):
    """sites: iterable of (pos0, ref, alt). Returns dict site -> (gnomad_af or -1, dbsnp, pon1000g, colors_pos, colors_maxaf)."""
    cache = pickle.load(open(CACHE, "rb")) if os.path.exists(CACHE) else {}
    todo = sorted(set(sites) - set(cache))
    if todo:
        v = {k: pysam.VariantFile(f"{PON}/{f}") for k, f in FILES.items()}
        for pos0, ref, alt in todo:
            g, d, p, cp, caf = -1.0, False, False, False, -1.0
            for r in v["gnomad"].fetch("chr1", pos0, pos0 + 1):
                if r.pos == pos0 + 1 and r.ref.upper() == ref and alt in [a.upper() for a in (r.alts or ())]:
                    g = max(g, float(r.info["AF"][[a.upper() for a in r.alts].index(alt)]))
            for r in v["dbsnp"].fetch("chr1", pos0, pos0 + 1):
                if r.pos == pos0 + 1 and r.ref.upper() == ref and alt in [a.upper() for a in (r.alts or ())]:
                    d = True
            for r in v["pon1000g"].fetch("chr1", pos0, pos0 + 1):
                if r.pos == pos0 + 1:
                    p = True
            for r in v["colors"].fetch("chr1", pos0, pos0 + 1):
                if r.pos == pos0 + 1:
                    cp = True
                    try:
                        caf = max(caf, max(float(x) for x in r.info["AF"]))
                    except Exception:
                        pass
            cache[(pos0, ref, alt)] = (g, d, p, cp, caf)
        pickle.dump(cache, open(CACHE, "wb"))
    return cache


def rule_repo(a):
    """scripts/filter_panel_of_normals.py defaults: gnomAD allele, dbSNP allele, 1000G position, CoLoRSdb position."""
    return a[0] >= 0 or a[1] or a[2] or a[3]


def rule_af001(a):
    """v2_pipeline chr1_pon.py (ClairS-TO-like): gnomAD allele AF>=0.001, dbSNP allele, 1000G position, CoLoRSdb AF>=0.001."""
    return a[0] >= 0.001 or a[1] or a[2] or (a[3] and a[4] >= 0.001)
