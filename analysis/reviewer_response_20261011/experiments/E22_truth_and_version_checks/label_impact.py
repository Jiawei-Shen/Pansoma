#!/usr/bin/env python3
"""E22 check 1: how many SNV tensors labelled somatic (1) in the CASTLE tensor sets carry a truth allele that is a
common gnomAD allele (AF >= 1 %) or lies in a normal-LOH segment (make_tables.py / summary_tables.json)?

Reads tensors/SNV/chr*_labels.ndjson of each set (read-only; only lines with "label": 1 are parsed) and the per-allele
tables W/<truth set>.alleles.tsv.gz. A tensor counts as "common" when any of its somatic truth alleles is common.
Writes HERE/label_impact.tsv and prints it.
"""
import csv
import glob
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
W = Path("/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks")
T = Path("/scratch/jshen/data/pansoma_v2_tensors")
SETS = {"H1437_CASTLE": "H1437T", "H2009_CASTLE": "H2009T", "HCC1954_CASTLE": "HCC1954T", "HCC1395_SEQC2": "HCC1395T"}


def main():
    st = json.loads((HERE / "summary_tables.json").read_text())
    rows_out = []
    for truth, prefix in SETS.items():
        common, truth_snvs = set(), set()
        with gzip.open(W / f"{truth}.alleles.tsv.gz", "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["kind"] != "SNV":
                    continue
                k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                truth_snvs.add(k)
                if r["gnomad_af"] and float(r["gnomad_af"]) >= 0.01:
                    common.add(k)
        segs = st[truth]["segments"]

        def in_seg(c, p):
            return any(c == g["chrom"] and g["start_mb"] * 1_000_000 <= p < g["end_mb"] * 1_000_000 for g in segs)

        for platform in ("PacBio", "ONT", "Illumina"):
            d = T / f"{prefix}_{platform}_d46" / "tensors" / "SNV"
            n = n_common = n_seg = n_test = n_test_common = 0
            uniq = set()
            for f in sorted(glob.glob(str(d / "chr*_labels.ndjson"))):
                with open(f) as fh:
                    for line in fh:
                        if '"label": 1,' not in line:
                            continue
                        rec = json.loads(line)
                        if rec["label"] != 1:
                            continue
                        truths = list(rec["somatic"])
                        if rec.get("partial_truth"):  # partial labels: the overlapped truth allele
                            truths.append(rec["partial_truth"])
                        keys = [(t["chrom"], t["vcf_pos"], t["vcf_ref"].upper(), t["vcf_alt"].upper())
                                for t in truths if t.get("vcf_pos") is not None]
                        is_common = any(k in common for k in keys)
                        seg = any(in_seg(k[0], k[1]) for k in keys)
                        n += 1
                        uniq.update(keys)
                        n_common += is_common
                        n_seg += seg
                        if rec["chrom"] == "chr1":
                            n_test += 1
                            n_test_common += is_common
            # truth-level: distinct truth SNV alleles with a somatic-labelled tensor, by class
            cls = {"common": [0, 0], "segment": [0, 0], "rest": [0, 0]}
            for k in truth_snvs:
                c = "segment" if in_seg(k[0], k[1]) else ("common" if k in common else "rest")
                cls[c][0] += 1
                cls[c][1] += k in uniq
            rows_out.append(dict(truth=truth, tensor_set=f"{prefix}_{platform}_d46", somatic_snv_tensors=n,
                                 truth_allele_common=n_common, in_normal_loh_segments=n_seg,
                                 chr1_somatic=n_test, chr1_common=n_test_common,
                                 truth_snvs_with_tensor=len(uniq & truth_snvs), truth_snvs=len(truth_snvs),
                                 segment_truth_with_tensor=f"{cls['segment'][1]}/{cls['segment'][0]}",
                                 common_outside_segments_with_tensor=f"{cls['common'][1]}/{cls['common'][0]}",
                                 rest_with_tensor=f"{cls['rest'][1]}/{cls['rest'][0]}"))
    cols = list(rows_out[0])
    with open(HERE / "label_impact.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows_out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print((HERE / "label_impact.tsv").read_text())


if __name__ == "__main__":
    main()
