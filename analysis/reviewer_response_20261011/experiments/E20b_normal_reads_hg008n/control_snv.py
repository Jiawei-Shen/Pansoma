"""E20b control: HG008-N normal and HG008-T tumor PacBio HiFi reads (linear GRCh38 BAMs, pileup, MAPQ >= 20) at the 140
absorbed truth SNVs and at 300 random truth SNVs that are not absorbed (PASS, in the benchmark BED, chr1-22, seed 1).
The normal is read remotely from GIAB (one connection, retries), +- 1 bp per site.
Writes $D/control/snv_pileup.tsv and prints the pooled and per-locus ALT fractions.
usage: CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt /wanglab/jshen/anaconda3/bin/python control_snv.py
"""
import collections
import csv
import random
import sys
import time

import pysam

D = "/scratch/jshen/data/reviewer_experiments_20261013/E20b_normal_reads_hg008n"
URL = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Revio_20240313/"
       "HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam")
BAI = f"{D}/idx/HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam.bai"
TUMOR = ("/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/"
         "HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam")
TRUTH = "/scratch/jshen/data/pansoma_v2_tensors/HG008T_PacBio_d46/truth/somatic.graph.tsv"
csv.field_size_limit(sys.maxsize)


def count(bam, c, p, alt):
    a = t = 0
    for col in bam.pileup(c, p, p + 1, truncate=True, min_base_quality=0, stepper="nofilter"):
        for pr in col.pileups:
            x = pr.alignment
            if x.is_secondary or x.is_supplementary or x.mapping_quality < 20 or pr.is_del or pr.is_refskip:
                continue
            t += 1
            a += x.query_sequence[pr.query_position] == alt
    return a, t


def main():
    absorbed = {r["truth_id"]: r for r in csv.DictReader(open(f"{D}/loci_d46.tsv"), delimiter="\t") if r["kind2"] == "SNV"}
    all_abs = {r["truth_id"] for r in csv.DictReader(open(f"{D}/loci_d46.tsv"), delimiter="\t")}
    pool = [r for r in csv.DictReader(open(TRUTH), delimiter="\t")
            if r["kind"] == "SNP" and r["passed"] == "True" and r["in_bed"] == "True" and r["truth_id"] not in all_abs
            and r["chrom"][3:].isdigit()]
    ctrl = random.Random(1).sample(pool, 300)
    sites = [("absorbed", r) for r in absorbed.values()] + [("control", r) for r in ctrl]
    tumor = pysam.AlignmentFile(TUMOR)
    normal = pysam.AlignmentFile(URL, index_filename=BAI)
    out = []
    for i, (grp, r) in enumerate(sites):
        c, p, alt = r["chrom"], int(r["vcf_pos"]) - 1, r["vcf_alt"]
        for attempt in range(8):
            try:
                n = count(normal, c, p, alt)
                break
            except (OSError, ValueError) as e:
                print("retry", i, c, p, attempt, e, file=sys.stderr, flush=True)
                time.sleep(5 * (attempt + 1))
                normal = pysam.AlignmentFile(URL, index_filename=BAI)
        else:
            raise SystemExit(f"failed {c}:{p}")
        t = count(tumor, c, p, alt)
        out.append(dict(group=grp, truth_id=r["truth_id"], chrom=c, pos=p + 1, alt=alt, normal_alt=n[0], normal_depth=n[1],
                        tumor_alt=t[0], tumor_depth=t[1]))
        if i % 100 == 0:
            print(i, flush=True)
    with open(f"{D}/control/snv_pileup.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    for grp in ("absorbed", "control"):
        rr = [x for x in out if x["group"] == grp]
        na, nd = sum(x["normal_alt"] for x in rr), sum(x["normal_depth"] for x in rr)
        ta, td = sum(x["tumor_alt"] for x in rr), sum(x["tumor_depth"] for x in rr)
        fr = [x["normal_alt"] / x["normal_depth"] for x in rr if x["normal_depth"] >= 5]
        b = collections.Counter("<0.1" if f < 0.1 else ">=0.9" if f >= 0.9 else "0.1-0.9" for f in fr)
        print(f"{grp}: {len(rr)} SNVs, normal ALT {na}/{nd} ({100 * na / max(nd, 1):.1f}%), tumor ALT {ta}/{td} "
              f"({100 * ta / max(td, 1):.1f}%), per-locus normal (>= 5 reads, {len(fr)}): {dict(b)}")


if __name__ == "__main__":
    main()
