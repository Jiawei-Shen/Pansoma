"""E9i (R3.3): matched-normal read evidence for COLO829T graph-coordinate SNV calls.

For each call of an E9c run (`<E9c DATA>/runs/<run>/calls.tsv`: unplaced PASS calls and the placed controls ctrl_tp /
ctrl_fp), count reads that contain the call's ALT 21-mer and its REF 21-mer (either strand) among the reads of a GRCh38
BAM that overlap the call's GRCh38 locus (E9c `locus_chrom` / `locus_lo` / `locus_hi`) +- 500 bp:
  normal  COLO829BL ONT (SMaHT D001 + D003, minimap2 to GRCh38), the matched normal;
  tumor   COLO829T ONT (`COLO829T_ONT_merged.bam`), as a check that the ALT is in the tumor's reads.
ALT 21-mer: positions site-10 .. site+10 of the call's ALT context (E9c `alt_read`, the consensus of the call's ALT reads;
`alt_graph` when missing); REF 21-mer: the same window of `ref_graph` (`ref_read` when missing). Primary records only,
MAPQ >= 0 (the locus is often in a region where the graph and GRCh38 differ).
Classes per call: normal ALT reads >= 2 -> 'ALT in normal'; 0 -> 'absent from normal' if normal REF + ALT reads >= 5,
else 'low normal coverage'; 1 -> 'one normal read'. Controls give the assay's error rates (true somatic placed TPs should
be absent from the normal, germline placed FPs present).
Writes <out>/<run>.tsv (or <out>/<run>.<i>.tsv for shard i of n: every n-th call). usage: python normal_reads.py RUN OUT_DIR [I N]
"""
import csv
import sys
from pathlib import Path

import pysam

E9C_RUNS = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation/runs")
NORMAL = [Path("/scratch/jshen/data/COLO829T/COLO829BL_ONT/SMHTCOLO829BL-X-X-M45-D001-uwsc-SMAFIOYE6UJE-sentieon_minimap2_202308.01_GRCh38.aligned.sorted.bam"),
          Path("/scratch/jshen/data/COLO829T/COLO829BL_ONT/SMHTCOLO829BL-X-X-M45-D003-uwsc-SMAFIAMORWKW-sentieon_minimap2_202308.01_GRCh38.aligned.sorted.bam")]
TUMOR = [Path("/scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.bam")]
K, PAD = 21, 500
GROUPS = ("unplaced_pass", "ctrl_tp", "ctrl_fp")
COMP = str.maketrans("ACGT", "TGCA")


def kmer(ctx, site):
    if not ctx or not site or "-" not in site:
        return None
    s = int(site.split("-")[0])
    lo, hi = s - K // 2, s + K // 2 + 1
    if lo < 0 or hi > len(ctx):
        return None
    k = ctx[lo:hi].upper()
    return k if set(k) <= set("ACGT") else None


def count(bams, chrom, lo, hi, alt, ref):
    alt_rc, ref_rc = alt.translate(COMP)[::-1], ref.translate(COMP)[::-1]
    n_alt = n_ref = n = 0
    for bam in bams:
        for r in bam.fetch(chrom, max(0, lo - PAD), hi + PAD):
            if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                continue
            q = r.query_sequence
            a = alt in q or alt_rc in q
            b = ref in q or ref_rc in q
            n += 1
            n_alt += a
            n_ref += b and not a
    return n_alt, n_ref, n


def classify(n_alt, n_ref):
    if n_alt >= 2:
        return "ALT in normal"
    if n_alt == 1:
        return "one normal read"
    return "absent from normal" if n_ref >= 5 else "low normal coverage"


def main():
    run, out = sys.argv[1], Path(sys.argv[2])
    shard, nshards = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (0, 1)
    out.mkdir(parents=True, exist_ok=True)
    normal = [pysam.AlignmentFile(str(p)) for p in NORMAL]
    tumor = [pysam.AlignmentFile(str(p)) for p in TUMOR]
    cache = {}
    rows = []
    with open(E9C_RUNS / run / "calls.tsv") as f:
        for i, r in enumerate(csv.DictReader(f, delimiter="\t")):
            if i % nshards != shard:
                continue
            if r["group"] not in GROUPS or r["kind_event"] != "SNP" or not r["locus_chrom"]:
                continue
            alt = kmer(r["alt_read"], r["alt_read_site"]) or kmer(r["alt_graph"], r["alt_graph_site"])
            ref = kmer(r["ref_graph"], r["ref_graph_site"]) or kmer(r["ref_read"], r["ref_read_site"])
            row = dict(group=r["group"], candidate_id=r["candidate_id"], locus=f"{r['locus_chrom']}:{r['locus_lo']}-{r['locus_hi']}",
                       locus_status=r["locus_status"], homopolymer=r["homopolymer"], str_unit=r["str_unit"])
            if not alt or not ref or alt == ref:
                row.update(status="no k-mer")
                rows.append(row)
                continue
            key = (r["locus_chrom"], int(r["locus_lo"]), int(r["locus_hi"]), alt, ref)
            if key not in cache:
                na, nr, nn = count(normal, *key)
                ta, tr, tn = count(tumor, *key)
                cache[key] = dict(normal_alt=na, normal_ref=nr, normal_reads=nn, tumor_alt=ta, tumor_ref=tr, tumor_reads=tn,
                                  status=classify(na, nr))
            row.update(alt_kmer=alt, ref_kmer=ref, **cache[key])
            rows.append(row)
    cols = ["group", "candidate_id", "locus", "locus_status", "homopolymer", "str_unit", "alt_kmer", "ref_kmer", "normal_alt",
            "normal_ref", "normal_reads", "tumor_alt", "tumor_ref", "tumor_reads", "status"]
    name = f"{run}.tsv" if nshards == 1 else f"{run}.{shard}.tsv"
    with open(out / (name + ".tmp"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    (out / (name + ".tmp")).rename(out / name)
    print(run, len(rows), "calls,", len(cache), "distinct loci")


if __name__ == "__main__":
    main()
