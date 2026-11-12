#!/usr/bin/env python3
"""Truth VAF tables, Illumina and PacBio, for HG008T, COLO829T and HCC1395T (chr1-22, every truth allele of the truth
VCF the labels use; SNV / INS / DEL, other alleles skipped).

  HG008T, HCC1395T: read counts of count_vaf.py (counts/<T>_<Illumina|PacBio>.truth_vaf.tsv).
  COLO829T: the SMaHT truth VCF's own VAF_Ill / VAF_PB (and RGN), the source of the earlier COLO829T violin; no read
            counts (no whole-genome COLO829T Illumina tumor BAM here).
Writes, next to this script:
  <T>_truth_VAF_Illumina_PacBio.csv  one row per truth allele: CHROM, POS, REF, ALT, TYPE, PASSED, IN_BED, [RGN,]
                                     DP_Ill, ALT_COUNT_Ill, VAF_Ill, DP_PB, ALT_COUNT_PB, VAF_PB
                                     (input of COLO829T_VAF_violin.R: columns VAF_Ill / VAF_PB)
  counts/COLO829T_<P>.truth_vaf.tsv   COLO829T in the count_vaf.py layout (DP / counts NA), so every
                                     counts/<T>_<P>.truth_vaf.tsv is an input of HG008T_VAF_violin.R (columns TYPE, VAF)
  summary.md                         n, alleles without reads, quartiles and low-VAF shares (DP > 0) per tumor,
                                     platform and type; the HCC1395T Illumina library check (SRR7890824 vs
                                     WGS_FD_T_1 + 2 on chr1)
"""
import csv
import gzip
from pathlib import Path

import numpy as np

from count_vaf import AUTOSOMES, classify, in_bed, read_bed

HERE = Path(__file__).resolve().parent
C = HERE / "counts"
COLO_VCF = "/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz"
COLO_BED = "/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/SMaHT_v2_easy_difficult_extreme.union.bed"
KEY = ("CHROM", "POS", "REF", "ALT")


def read_counts(path):
    return {tuple(r[k] for k in KEY): r for r in csv.DictReader(open(path), delimiter="\t")}


def colo829t():
    bed = read_bed(COLO_BED)
    rows = []
    with gzip.open(COLO_VCF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, pos, _, ref, alts, _, filt, info = line.rstrip("\n").split("\t")[:8]
            if c not in AUTOSOMES:
                continue
            inf = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
            for alt in alts.split(","):
                kind = classify(ref.upper(), alt.upper())
                if kind == "OTHER":
                    continue
                rows.append(dict(CHROM=c, POS=pos, REF=ref, ALT=alt, TYPE=kind, PASSED=str(filt in ("PASS", ".")),
                                 IN_BED=str(in_bed(bed, c, int(pos) - 1, int(pos) - 1 + len(ref))), RGN=inf.get("RGN", ""),
                                 DP_Ill="NA", ALT_COUNT_Ill="NA", VAF_Ill=inf.get("VAF_Ill", "NA"),
                                 DP_PB="NA", ALT_COUNT_PB="NA", VAF_PB=inf.get("VAF_PB", "NA")))
    return rows


def counted(t):
    ill, pb = read_counts(C / f"{t}_Illumina.truth_vaf.tsv"), read_counts(C / f"{t}_PacBio.truth_vaf.tsv")
    assert ill.keys() == pb.keys(), t
    return [dict(CHROM=a["CHROM"], POS=a["POS"], REF=a["REF"], ALT=a["ALT"], TYPE=a["TYPE"], PASSED=a["PASSED"],
                 IN_BED=a["IN_BED"], DP_Ill=a["DP"], ALT_COUNT_Ill=a["ALT_COUNT"], VAF_Ill=a["VAF"],
                 DP_PB=pb[k]["DP"], ALT_COUNT_PB=pb[k]["ALT_COUNT"], VAF_PB=pb[k]["VAF"]) for k, a in ill.items()]


def stats(rows, s):
    """Alleles with no read (DP 0, VAF written as 0 as in the earlier table) are counted apart, not in the quartiles."""
    dp0 = sum(r[f"DP_{s}"] == "0" for r in rows)
    v = np.array([float(r[f"VAF_{s}"]) for r in rows if r[f"VAF_{s}"] not in ("NA", "") and r[f"DP_{s}"] != "0"])
    if not len(v):
        return f"{len(rows):,} | {dp0} | | | | |"
    q = np.quantile(v, [0.25, 0.5, 0.75])
    return (f"{len(rows):,} | {dp0} | {q[0]:.3f} | {q[1]:.3f} | {q[2]:.3f} | {100 * np.mean(v < 0.05):.1f}% | "
            f"{100 * np.mean(v < 0.1):.1f}%")


def main():
    md = ["# Truth VAF, Illumina and PacBio (chr1-22, all truth alleles of the label truth VCF)", "",
          "Quartiles and low-VAF shares over the alleles with reads (DP > 0); COLO829T = SMaHT VAF_Ill / VAF_PB (no DP).", "",
          "| tumor | platform | type | n | DP = 0 | Q25 | median | Q75 | VAF < 0.05 | VAF < 0.1 |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for t in ("HG008T", "COLO829T", "HCC1395T"):
        rows = colo829t() if t == "COLO829T" else counted(t)
        cols = list(rows[0])
        with open(HERE / f"{t}_truth_VAF_Illumina_PacBio.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        if t == "COLO829T":
            for p, s in (("Illumina", "Ill"), ("PacBio", "PB")):
                with open(C / f"COLO829T_{p}.truth_vaf.tsv", "w") as f:
                    f.write("CHROM\tPOS\tREF\tALT\tTYPE\tPASSED\tIN_BED\tDP\tREF_COUNT\tALT_COUNT\tVAF\n")
                    for r in rows:
                        f.write(f"{r['CHROM']}\t{r['POS']}\t{r['REF']}\t{r['ALT']}\t{r['TYPE']}\t{r['PASSED']}\t"
                                f"{r['IN_BED']}\tNA\tNA\tNA\t{r['VAF_' + s]}\n")
        for p, s in (("Illumina", "Ill"), ("PacBio", "PB")):
            for typ, sel in (("SNV", ("SNV",)), ("INDEL", ("INS", "DEL"))):
                md.append(f"| {t} | {p} | {typ} | " + stats([r for r in rows if r["TYPE"] in sel], s) + " |")
        print(t, len(rows), flush=True)
    chk = C / "HCC1395T_Illumina_FDT_chr1.truth_vaf.tsv"
    if chk.exists():
        a, b = read_counts(C / "HCC1395T_Illumina.truth_vaf.tsv"), read_counts(chk)
        ks = [k for k in b if k in a and int(a[k]["DP"]) and int(b[k]["DP"])]
        x = np.array([float(a[k]["VAF"]) for k in ks])
        y = np.array([float(b[k]["VAF"]) for k in ks])
        dpa = np.median([int(a[k]["DP"]) for k in ks])
        dpb = np.median([int(b[k]["DP"]) for k in ks])
        md += ["", "HCC1395T Illumina library check, chr1 truth alleles with reads in both: SRR7890824 (used above) vs the "
               "WGS_FD_T_1 + 2 reads Pansoma uses (E8 HCC1395T_Illumina.chr1_chr20.bam), same counting:",
               f"n {len(ks):,}, Pearson r {np.corrcoef(x, y)[0, 1]:.3f}, median VAF {np.median(x):.3f} vs {np.median(y):.3f}, "
               f"median DP {dpa:.0f} vs {dpb:.0f}."]
    (HERE / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
