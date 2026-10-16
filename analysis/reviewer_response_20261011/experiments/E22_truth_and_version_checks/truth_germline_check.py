#!/usr/bin/env python3
"""E22 check 1: do the CASTLE somatic truth sets contain germline alleles of the matched normal or common population
alleles? (R2.2, R2.6, R1.m1)

Inputs (read-only): sets.tsv (truth VCF, BED, DeepSomatic full benchmark VCF, normal sample, somatic-only VCF, the
label germline VCF) and the subsets written by subset_resources.sh into W.

Per somatic truth allele (every ALT of every record of the truth VCF; SNV = 1-bp REF and ALT, everything else INDEL):
  in_bed            inside the somatic BED of the set (the labels' confident region; truth alleles are labelled 1
                    inside or outside it)
  normal_exact      the normal's DeepVariant-hybrid call (FILTER GERMLINE record of the full benchmark VCF) carries the
                    same allele (POS, REF, ALT) in its genotype; any QUAL; *_q10 = QUAL >= 10 (the label germline VCF)
  label_germ        the same allele is in the label germline VCF (cross-check of normal_exact_q10)
  normal_pos_other  a normal germline call at the same POS, other allele
  merged_germline   the allele is an ALT of a FILTER GERMLINE record but not in the normal's genotype (bcftools merge
                    of a somatic call into a germline site; the somatic-only VCF keeps it as PASS)
  gnomad_af, colorsdb_af   exact-allele AF (max over records) in the PoN resources
  tumor_vaf         mean of the per-technology tumor VAFs of the DeepSomatic somatic-only record (CASTLE sets; the
                    SEQC2 set takes them from the HCC1395 somatic-only VCF where the allele is there)
Regional: 1-Mb bins of the normal's genome-wide germline SNVs (GERMLINE records, QUAL >= 10, ALT in genotype);
het fraction = het / (het + hom-alt). A bin where the normal has lost one haplotype (LOH in the normal) has a het
fraction near 0. Truth SNVs are counted per het-fraction class of their bin.
Germline baseline: the share of each normal's chr20 germline SNVs (QUAL >= 10, in the set's BED) with gnomAD AF >= 1 %.

Writes: W/<set>.alleles.tsv.gz (one row per truth allele), HERE/summary.json, HERE/tables.md.
"""
from __future__ import annotations

import bisect
import csv
import gzip
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

HERE = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E22_truth_and_version_checks")
W = Path("/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks")
COMMON = 0.01
BIN = 1_000_000
HET_CLASSES = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 1.01)]
MIN_BIN_SNVS = 100  # bins with fewer normal germline SNVs get class "sparse"


def read_bed(path):
    iv = defaultdict(list)
    for line in open(path):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        c, s, e = line.split()[:3]
        iv[c].append((int(s), int(e)))
    out = {}
    for c, v in iv.items():
        v.sort()
        merged = []
        for s, e in v:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = ([m[0] for m in merged], [m[1] for m in merged])
    return out


def in_bed(bed, chrom, pos1):
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    i = bisect.bisect_right(starts, pos1 - 1) - 1
    return i >= 0 and ends[i] > pos1 - 1


def gt_alleles(gt: str):
    return {int(a) for a in gt.replace("|", "/").split("/") if a not in (".", "")}


def af_table(path):
    best = {}
    if not path.exists():
        raise FileNotFoundError(path)
    for line in open(path):
        c, p, r, alts, afs = line.rstrip("\n").split("\t")
        afl = afs.split(",")
        for i, a in enumerate(alts.split(",")):
            try:
                af = float(afl[i] if len(afl) > 1 else afl[0])
            except ValueError:
                continue
            k = (c, int(p), r.upper(), a.upper())
            best[k] = max(best.get(k, 0.0), af)
    return best


def tumor_vafs(path):
    out = {}
    with pysam.VariantFile(path) as v:
        samples = list(v.header.samples)
        for r in v:
            for i, a in enumerate(r.alts or ()):
                vals = []
                for s in samples:
                    smp = r.samples[s]
                    gt = smp.get("GT")
                    if gt is None or all(g is None for g in gt):
                        continue
                    vaf = smp.get("VAF")
                    if vaf is None:
                        continue
                    x = vaf[i] if isinstance(vaf, tuple) and len(vaf) > i else (vaf[0] if isinstance(vaf, tuple) else vaf)
                    if x is not None:
                        vals.append(float(x))
                if vals:
                    out[(r.chrom, r.pos, r.ref.upper(), a.upper())] = (statistics.fmean(vals), len(vals))
    return out


def normal_bins(path):
    """{(chrom, bin): [het, hom]} of the normal's germline SNVs with QUAL >= 10 and the ALT in the genotype;
    also the chr20 SNV keys {(pos, ref, alt): het?}."""
    bins = defaultdict(lambda: [0, 0])
    chr20 = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, p, r, alts, q, gt = line.rstrip("\n").split("\t")
            if len(r) != 1:
                continue
            try:
                if float(q) < 10:
                    continue
            except ValueError:
                continue
            al = gt_alleles(gt)
            alts = alts.split(",")
            snv_alts = [i for i in al if i > 0 and len(alts[i - 1]) == 1]
            if not snv_alts:
                continue
            hom = len(al) == 1 and 0 not in al
            het = 0 in al and len(al) == 2
            b = bins[(c, int(p) // BIN)]
            if hom:
                b[1] += 1
            elif het:
                b[0] += 1
            if c == "chr20":
                for i in snv_alts:
                    chr20[(int(p), r.upper(), alts[i - 1].upper())] = het
    return dict(bins), chr20


def het_class(bins, chrom, pos):
    b = bins.get((chrom, pos // BIN))
    if not b or sum(b) < MIN_BIN_SNVS:
        return "sparse"
    f = b[0] / sum(b)
    for lo, hi in HET_CLASSES:
        if lo <= f < hi:
            return f"het{lo:.1f}-{min(hi, 1.0):.1f}"
    return "sparse"


def analyse(row, line_bins, vaf_cache):
    name = row["name"]
    bed = read_bed(row["bed"])
    gn = af_table(W / f"{name}.gnomad.tsv")
    co = af_table(W / f"{name}.colorsdb.tsv")
    # full benchmark VCF records at the truth positions (normal column only)
    normal_exact, normal_pos, merged = {}, set(), set()
    full_filter = {}
    for line in open(W / f"{name}.full.tsv"):
        c, p, r, alts, q, flt, gt, vaf, ad = line.rstrip("\n").split("\t")
        p = int(p)
        al = gt_alleles(gt)
        for i, a in enumerate(alts.split(","), 1):
            k = (c, p, r.upper(), a.upper())
            full_filter.setdefault(k, set()).add(flt)
            if flt == "GERMLINE":
                if i in al:
                    qq = float(q) if q not in (".", "") else 0.0
                    normal_exact[k] = (qq, gt)
                else:
                    merged.add(k)
        if flt == "GERMLINE" and any(x > 0 for x in al):
            normal_pos.add((c, p))
    label_germ = set()
    lg = W / f"{name}.labelgerm.tsv"
    if lg.exists():
        for line in open(lg):
            c, p, r, alts, gt = line.rstrip("\n").split("\t")
            al = gt_alleles(gt)
            for i, a in enumerate(alts.split(","), 1):
                if i in al:
                    label_germ.add((c, int(p), r.upper(), a.upper()))
    sonly = row["somatic_only_vcf"]
    if sonly not in vaf_cache:
        vaf_cache[sonly] = tumor_vafs(sonly)
    vafs = vaf_cache[sonly]
    line = Path(row["full_vcf"]).name.split("_")[0]
    bins = line_bins[line]

    rows = []
    with pysam.VariantFile(row["truth_vcf"]) as v:
        for r in v:
            flts = list(r.filter.keys())
            if flts not in (["PASS"], []):
                continue
            for a in r.alts or ():
                if not a or a.startswith("<") or a in ("*", "."):
                    continue
                k = (r.chrom, r.pos, r.ref.upper(), a.upper())
                kind = "SNV" if len(k[2]) == 1 and len(k[3]) == 1 else "INDEL"
                ne = normal_exact.get(k)
                tv = vafs.get(k)
                rows.append(dict(
                    chrom=k[0], pos=k[1], ref=k[2], alt=k[3], kind=kind,
                    in_bed=int(in_bed(bed, k[0], k[1])),
                    normal_exact=int(ne is not None),
                    normal_exact_q10=int(ne is not None and ne[0] >= 10),
                    normal_gt=ne[1] if ne else "",
                    label_germ=int(k in label_germ),
                    normal_pos_other=int(ne is None and (k[0], k[1]) in normal_pos),
                    merged_germline=int(k in merged),
                    full_filter=";".join(sorted(full_filter.get(k, {"absent"}))),
                    gnomad_af=gn.get(k, ""), colorsdb_af=co.get(k, ""),
                    tumor_vaf=f"{tv[0]:.3f}" if tv else "", n_tech=tv[1] if tv else "",
                    normal_bin_class=het_class(bins, k[0], k[1]),
                ))
    cols = list(rows[0])
    with gzip.open(W / f"{name}.alleles.tsv.gz", "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in rows:
            fh.write("\t".join(str(x[c]) for c in cols) + "\n")
    return rows


def summarize(rows, kind, bed_only):
    sel = [x for x in rows if x["kind"] == kind and (x["in_bed"] or not bed_only)]
    n = len(sel)
    common = lambda x: x["gnomad_af"] != "" and x["gnomad_af"] >= COMMON
    colors_common = lambda x: x["colorsdb_af"] != "" and x["colorsdb_af"] >= COMMON
    s = Counter()
    for x in sel:
        s["n"] += 1
        s["normal_exact"] += x["normal_exact"]
        s["normal_exact_q10"] += x["normal_exact_q10"]
        s["label_germ"] += x["label_germ"]
        s["normal_exact_hom"] += int(x["normal_exact"] and len(gt_alleles(x["normal_gt"])) == 1)
        s["normal_pos_other"] += x["normal_pos_other"]
        s["merged_germline"] += x["merged_germline"]
        s["full_filter_GERMLINE"] += int("GERMLINE" in x["full_filter"])
        s["full_filter_PASS"] += int("PASS" in x["full_filter"])
        s["full_absent"] += int(x["full_filter"] == "absent")
        s["gnomad_ge1pct"] += int(common(x))
        s["gnomad_any"] += int(x["gnomad_af"] != "")
        s["colorsdb_ge1pct"] += int(colors_common(x))
        s["gnomad_or_colorsdb_ge1pct"] += int(common(x) or colors_common(x))
        s["gnomad_ge1pct_and_normal"] += int(common(x) and x["normal_exact"])
        s["gnomad_ge1pct_not_normal"] += int(common(x) and not x["normal_exact"])
        s["normal_or_gnomad_ge1pct"] += int(common(x) or x["normal_exact"])
    out = dict(s)
    # tumor VAF of SNVs, common vs not
    for lab, f in (("common", common), ("rare", lambda x: not common(x) and not x["normal_exact"])):
        v = [float(x["tumor_vaf"]) for x in sel if f(x) and x["tumor_vaf"] != ""]
        out[f"tumor_vaf_{lab}_n"] = len(v)
        out[f"tumor_vaf_{lab}_median"] = round(statistics.median(v), 3) if v else None
        out[f"tumor_vaf_{lab}_ge0.9"] = sum(1 for y in v if y >= 0.9)
    # normal het-fraction class of the bin
    cls = defaultdict(Counter)
    for x in sel:
        lab = "common_not_normal" if (common(x) and not x["normal_exact"]) else (
            "normal_exact" if x["normal_exact"] else "rare_not_normal")
        cls[x["normal_bin_class"]][lab] += 1
        cls[x["normal_bin_class"]]["all"] += 1
    out["bin_class"] = {k: dict(v) for k, v in sorted(cls.items())}
    # concentration: share of 1-Mb bins holding half of the common (not normal) vs the rare truth SNVs
    for lab, f in (("common_not_normal", lambda x: common(x) and not x["normal_exact"]),
                   ("rare_not_normal", lambda x: not common(x) and not x["normal_exact"])):
        b = Counter((x["chrom"], x["pos"] // BIN) for x in sel if f(x))
        tot, acc, k = sum(b.values()), 0, 0
        for c in sorted(b.values(), reverse=True):
            if acc >= tot / 2:
                break
            acc += c
            k += 1
        out[f"bins_half_{lab}"] = k
        out[f"bins_any_{lab}"] = len(b)
        top = [dict(chrom=c, start_mb=i, n=m) for (c, i), m in b.most_common(10)]
        out[f"top_bins_{lab}"] = top
    return out


def bin_genome_classes(bins):
    c = Counter()
    for (ch, i), b in bins.items():
        if ch in ("chrX", "chrY", "chrM") or "_" in ch:
            continue
        if sum(b) < MIN_BIN_SNVS:
            c["sparse"] += 1
            continue
        f = b[0] / sum(b)
        for lo, hi in HET_CLASSES:
            if lo <= f < hi:
                c[f"het{lo:.1f}-{min(hi, 1.0):.1f}"] += 1
    return dict(c)


def lowhet_segments(bins, thr=0.2):
    """Runs of consecutive autosomal 1-Mb bins with het fraction < thr (>= MIN_BIN_SNVS SNVs), sparse bins bridged."""
    segs = []
    by_chrom = defaultdict(dict)
    for (ch, i), b in bins.items():
        by_chrom[ch][i] = b
    for ch in sorted(by_chrom, key=lambda x: (len(x), x)):
        if ch in ("chrX", "chrY", "chrM") or "_" in ch:
            continue
        cur = None
        for i in range(max(by_chrom[ch]) + 1):
            b = by_chrom[ch].get(i)
            if not b or sum(b) < MIN_BIN_SNVS:
                continue
            low = b[0] / sum(b) < thr
            if low:
                if cur is None:
                    cur = [ch, i, i]
                else:
                    cur[2] = i
            else:
                if cur is not None:
                    segs.append(cur)
                    cur = None
        if cur is not None:
            segs.append(cur)
    return [dict(chrom=c, start_mb=s, end_mb=e + 1, mb=e + 1 - s) for c, s, e in segs]


def main():
    sets = list(csv.DictReader(open(HERE / "sets.tsv"), delimiter="\t"))
    lines = sorted({Path(r["full_vcf"]).name.split("_")[0] for r in sets})
    line_bins, chr20 = {}, {}
    for ln in lines:
        line_bins[ln], chr20[ln] = normal_bins(W / f"{ln}.normal_germline.tsv.gz")
        print(ln, "normal bins", len(line_bins[ln]), "chr20 SNVs", len(chr20[ln]), file=sys.stderr, flush=True)
    # germline baseline: chr20 germline SNVs of each normal with gnomAD AF >= 1 %, per set BED
    want = defaultdict(set)
    for ln in lines:
        for k in chr20[ln]:
            want[k].add(ln)
    gaf = {}
    with gzip.open(W / "gnomad_chr20_snps.tsv.gz", "rt") as fh:
        for line in fh:
            p, r, alts, afs = line.rstrip("\n").split("\t")
            afl = afs.split(",")
            for i, a in enumerate(alts.split(",")):
                k = (int(p), r.upper(), a.upper())
                if k in want:
                    try:
                        af = float(afl[i] if len(afl) > 1 else afl[0])
                    except ValueError:
                        continue
                    gaf[k] = max(gaf.get(k, 0.0), af)
    summary = dict(common_af=COMMON, bin_bp=BIN, min_bin_snvs=MIN_BIN_SNVS, sets={}, normals={})
    for ln in lines:
        summary["normals"][ln] = {"genome_bin_classes": bin_genome_classes(line_bins[ln]),
                                  "lowhet_segments_lt0.2": lowhet_segments(line_bins[ln])}
    vaf_cache = {}
    for row in sets:
        name = row["name"]
        print("set", name, file=sys.stderr, flush=True)
        rows = analyse(row, line_bins, vaf_cache)
        ln = Path(row["full_vcf"]).name.split("_")[0]
        bed = read_bed(row["bed"])
        base = Counter()
        for (p, r, a), het in chr20[ln].items():
            if not in_bed(bed, "chr20", p):
                continue
            base["n"] += 1
            base["gnomad_ge1pct"] += int(gaf.get((p, r, a), 0.0) >= COMMON)
            base["gnomad_any"] += int((p, r, a) in gaf)
            base["het_n"] += int(het)
            base["het_gnomad_ge1pct"] += int(het and gaf.get((p, r, a), 0.0) >= COMMON)
        summary["sets"][name] = dict(
            used_by_pansoma=row["used_by_pansoma"], truth_vcf=row["truth_vcf"], bed=row["bed"],
            normal_line=ln, germline_baseline_chr20_in_bed=dict(base),
            SNV_all=summarize(rows, "SNV", False), SNV_bed=summarize(rows, "SNV", True),
            INDEL_all=summarize(rows, "INDEL", False), INDEL_bed=summarize(rows, "INDEL", True))
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("wrote", HERE / "summary.json", file=sys.stderr)


if __name__ == "__main__":
    main()
