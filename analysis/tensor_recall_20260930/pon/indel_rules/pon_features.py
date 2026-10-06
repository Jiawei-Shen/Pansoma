"""PoN features of one allele and the INDEL PoN rules compared in this folder.

A feature tuple describes every PoN record that starts at the allele's POS (the only records the repo rule and
ClairS-TO look at):
  g_pos, g_af, g_bi         gnomAD: a record at POS; AF bucket of the exact allele (-1 = not there); exact allele in a
                            record with this single ALT (DeepSomatic matches the whole variant key, REF + all ALTs)
  d_pos, d_allele, d_nonsom, d_common, d_bi
                            dbSNP 138: record at POS; exact allele; exact allele not SAO=2; exact allele with COMMON=1
                            (>= 1% in a 1000G population); exact allele in a single-ALT record
  k_pos, k_allele, k_bi     1000G PoN (GATK, Mutect2 calls on 1000G normals): record at POS; exact allele; single ALT
  c_pos_af, c_af            CoLoRSdb: highest AF of any allele at POS (bucket); AF bucket of the exact allele
AF buckets: 0 < 1e-4, 1 [1e-4, 1e-3), 2 [1e-3, 0.01), 3 [0.01, 0.05), 4 [0.05, 0.1), 5 >= 0.1; -1 = no record / allele.
"""
from bisect import bisect_right

PON_DIR = "/scratch/jshen/data/Pansoma/panel_of_normal_VCFs"
PON_FILES = ["af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz", "1000g_pon.hg38.vcf.gz",
             "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"]
EDGES = [1e-4, 1e-3, 1e-2, 0.05, 0.1]
FIELDS = ("g_pos", "g_af", "g_bi", "d_pos", "d_allele", "d_nonsom", "d_common", "d_bi", "k_pos", "k_allele", "k_bi",
          "c_pos_af", "c_af")


def bucket(af):
    return -1 if af is None else bisect_right(EDGES, af)


def pon_record(rec, p):
    """(REF, ALTs, AFs, somatic, common) of a pysam record of PoN p (0 gnomAD, 1 dbSNP, 2 1000G, 3 CoLoRSdb)."""
    alts = tuple((a or "").upper() for a in (rec.alts or ()))
    afs = (None,) * len(alts)
    if p in (0, 3):
        v = rec.info.get("AF")
        afs = v if isinstance(v, tuple) else (v,) * len(alts)
    somatic = common = False
    if p == 1:
        somatic = rec.info.get("SAO") == 2
        common = rec.info.get("COMMON") == 1
    return (rec.ref or "").upper(), alts, afs, somatic, common


def features(records, ref, alt):
    """records: {p: [pon_record(...) of records starting at POS]}; ref/alt upper-case VCF alleles."""
    f = dict.fromkeys(FIELDS, 0)
    f["g_af"] = f["c_af"] = f["c_pos_af"] = -1
    for p, recs in records.items():
        for rref, alts, afs, somatic, common in recs:
            if p == 0:
                f["g_pos"] = 1
            elif p == 1:
                f["d_pos"] = 1
            elif p == 2:
                f["k_pos"] = 1
            else:
                f["c_pos_af"] = max(f["c_pos_af"], max((bucket(a) for a in afs if a is not None), default=-1))
            if rref != ref or alt not in alts:
                continue
            k = alts.index(alt)
            single = len(alts) == 1
            if p == 0:
                f["g_af"] = max(f["g_af"], bucket(afs[k]))
                f["g_bi"] |= single
            elif p == 1:
                f["d_allele"] = 1
                f["d_nonsom"] |= not somatic
                f["d_common"] |= common
                f["d_bi"] |= single
            elif p == 2:
                f["k_allele"] = 1
                f["k_bi"] |= single
            else:
                f["c_af"] = max(f["c_af"], bucket(afs[k]))
    return tuple(int(f[k]) for k in FIELDS)


def as_dict(t):
    return dict(zip(FIELDS, t))


# ---- rules: feature dict -> tagged? ------------------------------------------------------------------------------
B = {"1e-4": 1, "1e-3": 2, "0.01": 3, "0.05": 4, "0.1": 5}  # bucket index of an AF floor
RULES = {
    # what the repo uses now (scripts/filter_panel_of_normals.py)
    "repo (allele; gnomAD/CoLoRSdb AF>=1e-4; dbSNP non-somatic; 1000G)":
        lambda f: f["g_af"] >= B["1e-4"] or f["d_nonsom"] or f["k_allele"] or f["c_af"] >= B["1e-4"],
    # ClairS-TO defaults: sites files gnomAD/CoLoRSdb AF>=0.001; allele match gnomAD+dbSNP, position 1000G+CoLoRSdb
    "ClairS-TO (gnomAD AF>=1e-3 allele; dbSNP non-somatic allele; 1000G position; CoLoRSdb AF>=1e-3 position)":
        lambda f: f["g_af"] >= B["1e-3"] or f["d_nonsom"] or f["k_pos"] or f["c_pos_af"] >= B["1e-3"],
    "ClairS-TO AFs, allele match everywhere (the repo rule of a3365c7)":
        lambda f: f["g_af"] >= B["1e-3"] or f["d_nonsom"] or f["k_allele"] or f["c_af"] >= B["1e-3"],
    # DeepSomatic tumor-only: one PoN of dbSNP + gnomAD + 1000G calls, exact variant key (no CoLoRSdb)
    "DeepSomatic-like (exact variant key in dbSNP / gnomAD / 1000G; no CoLoRSdb)":
        lambda f: f["g_bi"] or f["d_bi"] or f["k_bi"],
    # Mutect2: hard PoN = the 1000G panel; gnomAD only as a germline prior (approximated by a common-allele floor)
    "Mutect2-like (1000G PoN allele + gnomAD AF>=0.01)":
        lambda f: f["k_allele"] or f["g_af"] >= B["0.01"],
    "1000G PoN only (allele)": lambda f: f["k_allele"],
    "repo without CoLoRSdb": lambda f: f["g_af"] >= B["1e-4"] or f["d_nonsom"] or f["k_allele"],
    "repo, dbSNP COMMON only": lambda f: f["g_af"] >= B["1e-4"] or f["d_common"] or f["k_allele"] or f["c_af"] >= B["1e-4"],
}
for floor in ("1e-4", "1e-3"):
    RULES[f"population AF >= {floor} only (gnomAD or CoLoRSdb allele)"] = (
        lambda f, b=B[floor]: f["g_af"] >= b or f["c_af"] >= b)
for floor in ("0.01", "0.05", "0.1"):
    RULES[f"population AF >= {floor} (gnomAD or CoLoRSdb allele) + 1000G + dbSNP COMMON"] = (
        lambda f, b=B[floor]: f["g_af"] >= b or f["c_af"] >= b or f["k_allele"] or f["d_common"])
    RULES[f"population AF >= {floor} only (gnomAD or CoLoRSdb allele)"] = (
        lambda f, b=B[floor]: f["g_af"] >= b or f["c_af"] >= b)
