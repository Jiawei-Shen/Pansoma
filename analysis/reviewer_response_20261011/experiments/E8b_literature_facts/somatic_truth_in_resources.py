#!/usr/bin/env python3
"""How many somatic truth alleles of each benchmark sit in each population / PoN resource (E8b, R2.6, R1.5).

For every ALT of every record of a truth VCF (no BED restriction; all chromosomes in the VCF, or one region),
query the four local resources (/scratch/jshen/data/Pansoma/panel_of_normal_VCFs) and count:
  - exact allele matches (POS, REF, ALT) per resource, with AF floors for gnomAD / CoLoRSdb;
  - dbSNP matches by SAO (variant allele origin: 0 unspecified, 1 germline, 2 somatic, 3 both);
  - position-only matches (any record at POS) for 1000G PoN and CoLoRSdb (ClairS-TO's default rule for these two);
  - the Pansoma rule, computed with the pipeline's own scripts/filter_panel_of_normals.py (imported, not copied);
  - the ClairS-TO default rule (README v0.5.1 / Nat Commun 2025 Methods): gnomAD AF >= 0.001 exact allele,
    dbSNP non-somatic exact allele, 1000G PoN position, CoLoRSdb AF >= 0.001 position.
A somatic truth allele matched by a resource is either a recurrent / coincident population allele, a germline
allele mislabelled somatic in the truth, or a somatic allele submitted to the resource (dbSNP SAO=2).

Usage: somatic_truth_in_resources.py OUT.tsv [--workers N]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from multiprocessing import Pool

import pysam

sys.path.insert(0, "/scratch/jshen/Github/Pansoma/scripts")
import filter_panel_of_normals as pon  # noqa: E402  (pipeline code, imported read-only)

PON_DIR = "/scratch/jshen/data/Pansoma/panel_of_normal_VCFs"
PONS = {
    "gnomAD": f"{PON_DIR}/af-only-gnomad.hg38.vcf.gz",
    "dbSNP": f"{PON_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
    "1000G": f"{PON_DIR}/1000g_pon.hg38.vcf.gz",
    "CoLoRSdb": f"{PON_DIR}/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz",
}
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
TRUTH = {  # name: (VCF, region or None)
    "HG008T_GIAB_v0.2": ("/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/"
                         "HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz", None),
    "COLO829T_SMaHT": ("/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz", None),
    "HCC1395_SEQC2_v1.2.1": ("/scratch/jshen/data/CASTLE/HCC1395/truth/HCC1395_somatic_SEQC2_v1.2.1.vcf.gz", None),
    "HCC1954_CASTLE": ("/scratch/jshen/data/CASTLE/HCC1954/truth/HCC1954_somatic_CASTLE.vcf.gz", None),
    "H1437_CASTLE": ("/scratch/jshen/data/CASTLE/H1437/truth/H1437_somatic_CASTLE.vcf.gz", None),
    "H2009_CASTLE": ("/scratch/jshen/data/CASTLE/H2009/truth/H2009_somatic_CASTLE.vcf.gz", None),
    # HapMap mixture positives are germline alleles of the minor donors: contrast, chr20 only
    "HapMap_positives_chr20": ("/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HapMap_somatic_snv_indel.vcf.gz",
                               "chr20"),
}


def matches(vcf: pysam.VariantFile, contig: str, pos: int, ref: str, alt: str):
    """(exact-allele hits as dicts, any record at POS)."""
    hits, at_pos = [], False
    for rec in vcf.fetch(contig, pos - 1, pos):
        if rec.pos != pos:
            continue
        at_pos = True
        if (rec.ref or "").upper() != ref:
            continue
        alts = [(a or "").upper() for a in (rec.alts or ())]
        if alt not in alts:
            continue
        k = alts.index(alt)
        af = rec.info.get("AF") if "AF" in rec.header.info else None
        if isinstance(af, tuple):
            af = af[k]
        ac = rec.info.get("AC") if "AC" in rec.header.info else None
        if isinstance(ac, tuple):
            ac = ac[k]
        sao = rec.info.get("SAO") if "SAO" in rec.header.info else None
        hits.append(dict(af=af, ac=ac, sao=sao))
    return hits, at_pos


def ge(af, floor):
    return af is not None and af >= floor * (1 - 1e-6)


def run_set(name: str) -> tuple[str, dict]:
    path, region = TRUTH[name]
    vcfs = {k: pysam.VariantFile(v) for k, v in PONS.items()}
    pon_list = [vcfs["gnomAD"], vcfs["dbSNP"], vcfs["1000G"], vcfs["CoLoRSdb"]]
    fasta = pysam.FastaFile(FASTA)
    counts: dict[str, Counter] = {"SNV": Counter(), "INDEL": Counter()}
    truth = pysam.VariantFile(path)
    it = truth.fetch(region) if region else truth.fetch()
    for rec in it:
        contig = rec.contig
        for alt in rec.alts or ():
            if not alt or alt.startswith("<") or alt in ("*", "."):
                continue
            ref, alt_u = rec.ref.upper(), alt.upper()
            kind = pon.record_kind(ref, [alt_u])
            c = counts[kind]
            c["n"] += 1
            g, _ = matches(vcfs["gnomAD"], contig, rec.pos, ref, alt_u)
            d, _ = matches(vcfs["dbSNP"], contig, rec.pos, ref, alt_u)
            t, t_pos = matches(vcfs["1000G"], contig, rec.pos, ref, alt_u)
            r, r_pos = matches(vcfs["CoLoRSdb"], contig, rec.pos, ref, alt_u)
            # one CoLoRSdb record per allele at most, but position-only needs the AF of any record at POS
            r_pos_af001 = False
            if r_pos:
                for prec in vcfs["CoLoRSdb"].fetch(contig, rec.pos - 1, rec.pos):
                    if prec.pos == rec.pos:
                        afs = prec.info.get("AF")
                        afs = afs if isinstance(afs, tuple) else (afs,)
                        if any(ge(a, 0.001) for a in afs):
                            r_pos_af001 = True
            c["gnomAD_any"] += bool(g)
            c["gnomAD_AF>=1e-4"] += any(ge(h["af"], 1e-4) for h in g)
            c["gnomAD_AF>=1e-3"] += any(ge(h["af"], 1e-3) for h in g)
            c["gnomAD_AF>=1e-2"] += any(ge(h["af"], 1e-2) for h in g)
            c["gnomAD_AC<=2"] += any(h["ac"] is not None and h["ac"] <= 2 for h in g)
            c["dbSNP_any"] += bool(d)
            c["dbSNP_SAO=2"] += any(h["sao"] == 2 for h in d)
            c["dbSNP_SAO=3"] += any(h["sao"] == 3 for h in d)
            c["dbSNP_SAO!=2"] += any(h["sao"] != 2 for h in d)
            c["1000G_allele"] += bool(t)
            c["1000G_position"] += t_pos
            c["CoLoRSdb_any"] += bool(r)
            c["CoLoRSdb_AF>=1e-4"] += any(ge(h["af"], 1e-4) for h in r)
            c["CoLoRSdb_AF>=1e-3"] += any(ge(h["af"], 1e-3) for h in r)
            c["CoLoRSdb_AC=1"] += any(h["ac"] == 1 for h in r)
            c["CoLoRSdb_position"] += r_pos
            # Pansoma rule, by the pipeline's own functions
            rules = pon.PON_RULES[kind]
            tagged = False
            for pvcf, rule in zip(pon_list, rules):
                if rule is None:
                    continue
                if rule == "indel_af":
                    kw = dict(min_af=pon.indel_min_af(fasta, contig, rec.pos, ref, alt_u))
                else:
                    kw = rule
                if pon.record_matches(pvcf, contig, rec.pos, ref, [alt_u], **kw):
                    tagged = True
                    break
            c["Pansoma_rule"] += tagged
            # ClairS-TO default rule
            clairs = (any(ge(h["af"], 1e-3) for h in g) or any(h["sao"] != 2 for h in d) or t_pos or r_pos_af001)
            c["ClairS-TO_rule"] += bool(clairs)
    return name, counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--sets", nargs="*", default=list(TRUTH))
    a = ap.parse_args()
    with Pool(a.workers) as pool:
        results = dict(pool.map(run_set, a.sets))
    cols = ["n", "Pansoma_rule", "ClairS-TO_rule", "gnomAD_any", "gnomAD_AF>=1e-4", "gnomAD_AF>=1e-3",
            "gnomAD_AF>=1e-2", "gnomAD_AC<=2", "dbSNP_any", "dbSNP_SAO=2", "dbSNP_SAO=3", "dbSNP_SAO!=2",
            "1000G_allele", "1000G_position", "CoLoRSdb_any", "CoLoRSdb_AF>=1e-4", "CoLoRSdb_AF>=1e-3",
            "CoLoRSdb_AC=1", "CoLoRSdb_position"]
    with open(a.out, "w") as fh:
        fh.write("truth_set\tkind\t" + "\t".join(cols) + "\n")
        for name in a.sets:
            for kind in ("SNV", "INDEL"):
                c = results[name][kind]
                fh.write(f"{name}\t{kind}\t" + "\t".join(str(c[k]) for k in cols) + "\n")


if __name__ == "__main__":
    main()
