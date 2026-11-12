#!/usr/bin/env python3
"""Read-count VAF of every somatic truth allele in one tumor BAM (GRCh38), chr1-22.

The counting rules are those of the earlier HG008T truth-VAF table
(/scratch/jshen/tmp/HG008_trueset_germline_overlapping/count_vaf_from_vcf_with_pysam.py), unchanged:
  - pysam pileup at the VCF position (stepper samtools, truncate), primary alignments only, MAPQ >= 20, base quality >= 20
    at the read base over POS, overlapping mates not merged;
  - SNV: DP = reads with any base at POS, REF_COUNT / ALT_COUNT = reads with the REF / ALT base;
  - INS / DEL (REF / ALT share the anchor base): DP = reads with a base at the anchor, ALT_COUNT = reads whose indel
    right after the anchor has the truth length (+len inserted / -len deleted), REF_COUNT = reads with no indel there;
  - VAF = ALT_COUNT / DP (0 when DP = 0); multi-ALT records are split; other alleles (MNP, complex) are skipped.
Added: PASSED (FILTER PASS or .) and IN_BED (the REF span, anchor base included, lies inside one merged interval of the
truth BED), and one worker per chromosome.

usage: count_vaf.py VCF BAM OUT.tsv [--bed BED] [--chroms chr1,...] [--threads 8] [--min-mapq 20] [--min-baseq 20]
"""
import argparse
import gzip
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
import pysam

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
ARGS = None
BED = None


def classify(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt) and alt.startswith(ref) and len(ref) == 1:
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt) and len(alt) == 1:
        return "DEL"
    return "OTHER"


def read_bed(path):
    """Merged half-open intervals per chromosome -> (starts, ends)."""
    iv = defaultdict(list)
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                c, s, e = line.split("\t")[:3]
                iv[c].append((int(s), int(e)))
    out = {}
    for c, items in iv.items():
        merged = []
        for s, e in sorted(items):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = (np.array([m[0] for m in merged]), np.array([m[1] for m in merged]))
    return out


def in_bed(bed, chrom, start, end):
    if bed is None or chrom not in bed:
        return False
    starts, ends = bed[chrom]
    k = np.searchsorted(starts, start, side="right") - 1
    return bool(k >= 0 and ends[k] >= end)


def read_vcf(path, chroms):
    recs = defaultdict(list)
    stats = Counter()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, pos, _, ref, alts, _, filt = line.split("\t", 7)[:7]
            if c not in chroms:
                stats["other_chromosome"] += 1
                continue
            for alt in alts.split(","):
                if alt in (".", "*") or alt.startswith("<"):
                    continue
                kind = classify(ref.upper(), alt.upper())
                stats[kind] += 1
                if kind != "OTHER":
                    recs[c].append((int(pos), ref.upper(), alt.upper(), kind, filt in ("PASS", ".")))
    return recs, stats


def count(bam, chrom, pos, ref, alt, kind):
    pos0 = pos - 1
    dp = rc = ac = 0
    for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper="samtools",
                          min_base_quality=ARGS.min_baseq, min_mapping_quality=ARGS.min_mapq,
                          ignore_overlaps=False, ignore_orphans=False):
        if col.reference_pos != pos0:
            continue
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            aln = pr.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.mapping_quality < ARGS.min_mapq:
                continue
            q = pr.query_position
            if q is None:
                continue
            bq = aln.query_qualities[q] if aln.query_qualities is not None else 0
            if bq < ARGS.min_baseq:
                continue
            dp += 1
            if kind == "SNV":
                b = aln.query_sequence[q].upper()
                rc += b == ref
                ac += b == alt
            else:
                want = len(alt) - len(ref) if kind == "INS" else -(len(ref) - len(alt))
                if pr.indel == want:
                    ac += 1
                elif pr.indel == 0:
                    rc += 1
    return dp, rc, ac


def work(item):
    chrom, recs = item
    bam = pysam.AlignmentFile(ARGS.bam, "rb")
    rows = []
    for pos, ref, alt, kind, passed in recs:
        dp, rc, ac = count(bam, chrom, pos, ref, alt, kind)
        rows.append((chrom, pos, ref, alt, kind, passed, in_bed(BED, chrom, pos - 1, pos - 1 + len(ref)), dp, rc, ac))
    bam.close()
    print(f"[INFO] {chrom}: {len(rows)} alleles", file=sys.stderr, flush=True)
    return rows


def main():
    global ARGS, BED
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf")
    p.add_argument("bam")
    p.add_argument("out")
    p.add_argument("--bed")
    p.add_argument("--chroms", default=",".join(AUTOSOMES))
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=20)
    ARGS = p.parse_args()
    BED = read_bed(ARGS.bed) if ARGS.bed else None
    chroms = ARGS.chroms.split(",")
    recs, stats = read_vcf(ARGS.vcf, set(chroms))
    print(f"[INFO] alleles: {dict(stats)}", file=sys.stderr, flush=True)
    items = sorted(recs.items(), key=lambda kv: -len(kv[1]))
    with Pool(ARGS.threads) as pool:
        results = {rows[0][0]: rows for rows in pool.imap_unordered(work, items) if rows}
    with open(ARGS.out, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tTYPE\tPASSED\tIN_BED\tDP\tREF_COUNT\tALT_COUNT\tVAF\n")
        for c in chroms:
            for chrom, pos, ref, alt, kind, passed, bed, dp, rc, ac in results.get(c, []):
                out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{kind}\t{passed}\t{bed}\t{dp}\t{rc}\t{ac}\t"
                          f"{ac / dp if dp else 0.0:.6f}\n")
    print(f"[INFO] wrote {ARGS.out}; skipped {stats.get('OTHER', 0)} other alleles, "
          f"{stats.get('other_chromosome', 0)} records off chr1-22", file=sys.stderr)


if __name__ == "__main__":
    main()
