"""Comparator VCF with an allele-fraction floor: a PASS record whose allele fraction (ClairS-TO FORMAT AF, DeepSomatic
FORMAT VAF, first sample) is below MIN_AF becomes FILTER LowAF; everything else is unchanged.
usage: rewrite_minaf.py CALLER MIN_AF IN.vcf.gz OUT.vcf (plain text)"""
import collections, gzip, sys
caller, min_af, src, dst = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
key = {"clairs_to": "AF", "deepsomatic": "VAF"}[caller]
n = collections.Counter()
with gzip.open(src, "rt") as i, open(dst, "w") as o:
    for line in i:
        if line[0] == "#":
            if line.startswith("#CHROM"):
                o.write(f'##FILTER=<ID=LowAF,Description="Allele fraction ({key}) below {min_af}">\n')
            o.write(line); continue
        f = line.rstrip("\n").split("\t")
        n["records"] += 1
        if f[6] == "PASS":
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            v = fmt.get(key, ".").split(",")[0]
            if v not in (".", "") and float(v) < min_af:
                f[6] = "LowAF"; n["PASS->LowAF"] += 1
            else:
                n["PASS kept"] += 1
        o.write("\t".join(f) + "\n")
print(src, dict(n), flush=True)
