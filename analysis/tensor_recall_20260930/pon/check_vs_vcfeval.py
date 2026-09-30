"""Check pon_scan.py's tensor alleles and rule against a real vcfeval PoN output (chr1).

Every record of <run>/vcf_chr1/eval_bed/<set>.<kind>.linear.pon.vcf.gz (linear_vcf + filter_panel_of_normals.py) must be
an allele of scan/chr1.keys.tsv, with the same PoNs (INFO PANSOMA_PON) as the scan's rule bits.
"""
import sys
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
NAMES = ["PoN1_gnomAD", "PoN2_dbSNP", "PoN3_1000G", "PoN4_CoLoRSdb"]
keys = {}
with open(HERE / "scan/chr1.keys.tsv") as f:
    next(f)
    for line in f:
        pos, ref, alt, _, _, rule = line.rstrip("\n").split("\t")
        keys[(int(pos), ref, alt)] = int(rule)
for path in sys.argv[1:]:
    n = found = same = 0
    with pysam.VariantFile(path) as vcf:
        for r in vcf.fetch("chr1"):
            n += 1
            k = (r.pos, r.ref, r.alts[0])
            if k not in keys:
                continue
            found += 1
            tags = r.info.get("PANSOMA_PON") or ()
            bits = sum(1 << NAMES.index(t) for t in tags)
            same += bits == keys[k]
    print(f"{path}: {n:,} records, {found:,} are scan alleles, {same:,} with the same PoN tags")
