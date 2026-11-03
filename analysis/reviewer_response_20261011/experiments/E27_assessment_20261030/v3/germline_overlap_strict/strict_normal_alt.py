"""Strict germline overlap: does a matched-normal haplotype carry GRCh38 + truth over the whole tandem array?
Same rules as analysis/graph_absorbed_somatic_colo829t_20261005/c3_normal_frame.py (array growth, GRCh38 anchors, dipcall
haplotype reconstruction, germline differences outside the array reverted by the a3 aligner), but with the dipcall
haplotype for every locus (no assembly placement), so that both tumors are treated alike. Every PASS chr1-22 truth allele
of tensor_recall_20260930/per_truth_<T>.tsv.
usage: strict_normal_alt.py TUMOR DIP_VCF(indexed) DIP_BED OUT.tsv [LIMIT]
  HG008T: /scratch/jshen/data/reviewer_experiments_20261013/germline_overlap_strict/HG008N.dip.vcf.gz (an indexed copy of
          HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz) and that directory's .dip.bed
  COLO829T: /scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.vcf.gz and .dip.bed (read only)"""
import bisect, collections, csv, sys
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_colo829t_20261005")
import pysam
import c3_normal_frame as c3

T, VCF, BED, OUT = sys.argv[1:5]
LIMIT = int(sys.argv[5]) if len(sys.argv) > 5 else 10 ** 9
n_done = 0
gfa = pysam.FastaFile(c3.GR)
dvcf = pysam.VariantFile(VCF)
bed = collections.defaultdict(list)
for ln in open(BED):
    c, s, e = ln.split()[:3]
    bed[c].append((int(s), int(e)))
for c in bed:
    bed[c].sort()
inbed = lambda c, s, e: (lambda i: i >= 0 and bed[c][i][0] <= s and e <= bed[c][i][1])(bisect.bisect_right(bed[c], (s, 10 ** 12)) - 1)

cols = ["truth_id", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "kind", "in_bed", "arr_lo", "arr_hi", "anchor_a", "anchor_b",
        "h0_dip", "h1_dip", "h0_na", "h1_na", "h0_allele", "h1_allele", "normal_carries_alt"]
with open(OUT, "w") as w:
    w.write("\t".join(cols) + "\n")
    for r in csv.DictReader(open(f"/scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/per_truth_{T}.tsv"), delimiter="\t"):
        if r["passed"] != "True" or n_done >= LIMIT:
            continue
        L = dict(truth_id=r["truth_id"], chrom=r["chrom"], vcf_pos=r["vcf_pos"], vcf_ref=r["vcf_ref"].upper(), vcf_alt=r["vcf_alt"].upper())
        p0, vr, va = int(L["vcf_pos"]) - 1, L["vcf_ref"], L["vcf_alt"]
        g, g0, lo, hi, capped = c3.array(gfa, L)
        a, b, _ = c3.anchors(g, g0, lo, hi, {})
        o = dict(truth_id=r["truth_id"], chrom=r["chrom"], vcf_pos=r["vcf_pos"], vcf_ref=vr, vcf_alt=va, kind=r["kind"],
                 in_bed=r["in_bed"], arr_lo=lo, arr_hi=hi, anchor_a=a if a is not None else "", anchor_b=b if b is not None else "")
        alle = {}
        if a is None or b is None:
            for h in (0, 1):
                o[f"h{h}_dip"], o[f"h{h}_na"], alle[h] = "NA", "no_anchors", "NA"
        else:
            R = g[a - g0:b - g0]; A = c3.edit(R, a, [(p0, vr, va)])
            recs = [(x.pos - 1, x.ref.upper(), tuple((z or "").upper() for z in x.alts or ()), x.samples[0]["GT"],
                     not x.filter.keys() or "PASS" in x.filter.keys()) for x in dvcf.fetch(r["chrom"], a, b)]
            ib = inbed(r["chrom"], a, b)
            for h in (0, 1):
                ds, why = c3.dip_hap(R, a, b, recs, h, ib)
                o[f"h{h}_dip"] = "NA" if ds is None else "ALT" if ds == A else "REF" if ds == R else "other"
                o[f"h{h}_na"] = why
                if ds is None:
                    alle[h] = "NA"; continue
                dr = c3.diff_runs(ds, R)
                if dr is None:
                    core = ds
                else:
                    arr = [d for d in dr if d[0] <= hi - a and d[1] >= lo - a]
                    core = c3.edit(R, 0, [(r0, R[r0:r1], x) for r0, r1, x in arr])
                alle[h] = "ALT" if core == A else "REF" if core == R else "germline_len" if len(core) != len(R) else "germline_seq"
        o["h0_allele"], o["h1_allele"] = alle[0], alle[1]
        o["normal_carries_alt"] = "ALT" in (alle[0], alle[1])
        w.write("\t".join(str(o[c]) for c in cols) + "\n"); n_done += 1
