"""Exact-allele overlap of a somatic truth VCF with the matched normal's dipcall germline VCF, chr1-22.
Alleles are split, trimmed and left-aligned on GRCh38; MNVs split into SNVs (the looser, exact-allele definition).
usage: overlap.py SOM_VCF GERM_VCF SOM_BED OUT_ROWS.tsv"""
import sys, collections, pysam
REF = pysam.FastaFile("/scratch/jshen/REF/GRCh38.primary.fasta")
CHR = [f"chr{i}" for i in range(1, 23)]

def alleles(rec, seq):
    for i, alt in enumerate(rec.alts or (), 1):
        if alt.startswith("<") or alt == "*":
            continue
        r, a, p = rec.ref.upper(), alt.upper(), rec.pos - 1
        while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
            r, a = r[:-1], a[:-1]
        while r and a and r[0] == a[0]:
            r, a, p = r[1:], a[1:], p + 1
        if len(r) == len(a):
            for k, (x, y) in enumerate(zip(r, a)):
                if x != y:
                    yield (rec.chrom, p + k, x, y, "SNV"), i
        elif r and a:
            yield (rec.chrom, p, r, a, "COMPLEX"), i
        else:
            s = r or a
            while p > 0 and seq[p - 1] == s[-1]:
                p, s = p - 1, seq[p - 1] + s[:-1]
            yield (rec.chrom, p, s if r else "", "" if r else s, "INDEL"), i

def load(path, gt=False):
    out = {}
    vf = pysam.VariantFile(path)
    sample = list(vf.header.samples)[0] if vf.header.samples else None
    seq, cur = None, None
    for rec in vf:
            if rec.chrom not in CHR:
                continue
            if rec.chrom != cur:
                seq, cur = REF.fetch(rec.chrom).upper(), rec.chrom
            filt = ";".join(rec.filter.keys()) or "PASS"
            g = rec.samples[sample].get("GT") if sample else None
            for key, i in alleles(rec, seq):
                out[key] = (filt, g, i, rec.pos)
    return out

som, germ = load(sys.argv[1]), load(sys.argv[2])
bed = collections.defaultdict(list)
for line in open(sys.argv[3]):
    c, s, e = line.split()[:3]
    bed[c].append((int(s), int(e)))
import bisect
for c in bed:
    bed[c].sort()
starts = {c: [s for s, _ in v] for c, v in bed.items()}
def inbed(c, p):
    i = bisect.bisect_right(starts.get(c, []), p) - 1
    return i >= 0 and bed[c][i][1] > p
cnt = collections.Counter()
rows = []
for key, (f, g, i, pos) in som.items():
    passed = f in ("PASS", ".")
    kind = key[4]
    cnt[(kind, "all", passed)] += 1
    if key in germ:
        gf, gg, gi, _ = germ[key]
        cnt[(kind, "in_normal", passed)] += 1
        cnt[(kind, "in_normal_bed" if inbed(key[0], key[1]) else "in_normal_outbed", passed)] += 1
        rows.append("\t".join(map(str, key + (f, gf, gg, inbed(key[0], key[1]))))) 
for k in sorted(cnt):
    print(*k, cnt[k], sep="\t")
open(sys.argv[4], "w").write("\n".join(rows) + "\n")
