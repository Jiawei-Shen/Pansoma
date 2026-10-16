#!/usr/bin/env python3
"""Is an individual's genome inside a population resource?  Private-allele test on chr20 SNVs (E8b, R2.6, R1.5).

The sample lists of the GATK af-only gnomAD file and of CoLoRSdb v1.1.0 are not published per individual.
If an individual contributed to a resource, its private alleles appear there as singletons (AC = 1 for a
heterozygous allele, AC = 2 for a homozygous one); if it did not, its private alleles are mostly absent and a
singleton match needs another carrier inside the resource. So the number of the individual's SNV alleles that
match a resource record with AC <= 2, per 1,000 of the individual's alleles, separates members from
non-members; individuals of the same ancestry serve as controls for each other.

Individuals and their germline SNV alleles on chr20 (GT allele index > 0; SNVs only, so no normalisation issue):
  HapMap donors HG002, HG00438, HG02257, HG02486, HG02622: SMaHT HAPMAP6 v1.6 records (INFO/<donor>_GT),
    i.e. their alleles absent from HG005 (the file holds only those);
  HG005: Pansoma HapMap_truth/HG005_germline_snv_indel.vcf.gz (INFO/HG005_GT);
  HG008-N: dipcall of the HG008-N assembly; COLO829BL: dipcall of the COLO829BL assembly;
  HCC1395BL / HCC1954 / H1437 / H2009 normals: CASTLE (DeepSomatic) germline VCFs, QUAL >= 10.
dbSNP 138 has no AC: reported is the fraction of the individual's gnomAD-absent alleles that dbSNP holds
(a 1000 Genomes phase-1 submission would put a donor's private alleles in dbSNP).

Usage: membership_private_alleles.py OUT.tsv [--contig chr20]
"""
from __future__ import annotations

import argparse
from collections import Counter

import pysam

PON_DIR = "/scratch/jshen/data/Pansoma/panel_of_normal_VCFs"
RESOURCES = {
    "gnomAD": f"{PON_DIR}/af-only-gnomad.hg38.vcf.gz",
    "CoLoRSdb": f"{PON_DIR}/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz",
}
DBSNP = f"{PON_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
HAPMAP6 = "/scratch/jshen/data/HapMap/HapMap_trueset/SMHTHAPMAP6_GRCh38_v1.6_somatic_benchmark_snvs.vcf.gz"
INDIVIDUALS = {  # name: (VCF, mode, sample or INFO key)
    "HG002": (HAPMAP6, "info", "HG002_GT"),
    "HG00438": (HAPMAP6, "info", "HG00438_GT"),
    "HG02257": (HAPMAP6, "info", "HG02257_GT"),
    "HG02486": (HAPMAP6, "info", "HG02486_GT"),
    "HG02622": (HAPMAP6, "info", "HG02622_GT"),
    "HG005": ("/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HG005_germline_snv_indel.vcf.gz", "info", "HG005_GT"),
    "HG008-N": ("/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz", "gt", None),
    "COLO829BL": ("/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.vcf.gz", "gt", None),
    "HCC1395BL": ("/scratch/jshen/data/CASTLE/HCC1395/truth/HCC1395BL_germline_DeepVariant_hybrid_QUAL10.vcf.gz", "gt", None),
    "HCC1954_normal": ("/scratch/jshen/data/CASTLE/HCC1954/truth/HCC1954_germline_DeepVariant_hybrid_QUAL10.vcf.gz", "gt", None),
    "H1437_normal": ("/scratch/jshen/data/CASTLE/H1437/truth/H1437_germline_DeepVariant_hybrid_QUAL10.vcf.gz", "gt", None),
    "H2009_normal": ("/scratch/jshen/data/CASTLE/H2009/truth/H2009_germline_DeepVariant_hybrid_QUAL10.vcf.gz", "gt", None),
}


def snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and alt not in "*." and ref != alt


def records(path: str, contig: str):
    vcf = pysam.VariantFile(path)
    try:
        yield from vcf.fetch(contig)
    except ValueError:  # no index: stream the file
        for rec in vcf:
            if rec.contig == contig:
                yield rec


def alleles(path: str, mode: str, key: str | None, contig: str) -> dict[tuple[int, str, str], int]:
    """(pos, ref, alt) -> copies of the ALT carried (1 het, 2 hom)."""
    out: dict[tuple[int, str, str], int] = {}
    for rec in records(path, contig):
        if mode == "info":
            gt = rec.info.get(key)
            if gt is None:
                continue
            gt = gt[0] if isinstance(gt, tuple) else gt
            idx = [int(x) for x in str(gt).replace("|", "/").split("/") if x.isdigit()]
        else:
            smp = rec.samples[0]
            idx = [x for x in (smp.get("GT") or ()) if x is not None]
        for a in set(i for i in idx if i > 0):
            if a - 1 >= len(rec.alts or ()):
                continue
            ref, alt = rec.ref.upper(), rec.alts[a - 1].upper()
            if snv(ref, alt):
                out[(rec.pos, ref, alt)] = idx.count(a)
    return out


def load_resource(path: str, contig: str) -> dict[tuple[int, str, str], int]:
    """(pos, ref, alt) -> AC of that ALT (SNV ALTs only)."""
    out = {}
    for rec in pysam.VariantFile(path).fetch(contig):
        ref = rec.ref.upper()
        if len(ref) != 1:
            continue
        acs = rec.info.get("AC")
        acs = acs if isinstance(acs, tuple) else (acs,) * len(rec.alts or ())
        for alt, ac in zip(rec.alts or (), acs):
            alt = alt.upper()
            if snv(ref, alt) and ac is not None:
                out[(rec.pos, ref, alt)] = max(out.get((rec.pos, ref, alt), 0), int(ac))
    return out


def load_dbsnp(contig: str) -> set[tuple[int, str, str]]:
    out = set()
    for rec in pysam.VariantFile(DBSNP).fetch(contig):
        ref = rec.ref.upper()
        if len(ref) != 1:
            continue
        for alt in rec.alts or ():
            if snv(ref, alt.upper()):
                out.add((rec.pos, ref, alt.upper()))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--contig", default="chr20")
    a = ap.parse_args()
    res = {name: load_resource(path, a.contig) for name, path in RESOURCES.items()}
    dbsnp = load_dbsnp(a.contig)
    cols = ["individual", "snv_alleles"]
    for r in RESOURCES:
        cols += [f"{r}_absent", f"{r}_AC1", f"{r}_AC2", f"{r}_AC3-10", f"{r}_AC>10", f"{r}_AC<=2_per_1000"]
    cols += ["gnomAD_absent_in_dbSNP", "gnomAD_absent_in_dbSNP_frac"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for name, (path, mode, key) in INDIVIDUALS.items():
            al = alleles(path, mode, key, a.contig)
            row = [name, str(len(al))]
            for r, table in res.items():
                c = Counter()
                for k in al:
                    ac = table.get(k)
                    c["absent" if ac is None else "AC1" if ac == 1 else "AC2" if ac == 2
                      else "AC3-10" if ac <= 10 else "AC>10"] += 1
                low = c["AC1"] + c["AC2"]
                row += [str(c["absent"]), str(c["AC1"]), str(c["AC2"]), str(c["AC3-10"]), str(c["AC>10"]),
                        f"{1000 * low / max(1, len(al)):.2f}"]
            absent = [k for k in al if k not in res["gnomAD"]]
            n_db = sum(k in dbsnp for k in absent)
            row += [str(n_db), f"{n_db / max(1, len(absent)):.4f}"]
            fh.write("\t".join(row) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
