#!/usr/bin/env python3
"""E11b step 4: score one (arm, platform) alignment set against the simulation truth.

Per read (primary record only):
  graph view (GAF from the GAM): mapped, MAPQ, perfect (whole read aligned, no mismatch / indel: cs has only ':'),
      mismatches, off-reference nodes (node not on the GRCh38 path of that graph, from the pipeline's grch38_path index)
  GRCh38 view (vg surject BAM): placed, chrom, start, end (CHM13 for the chm13 control arm)
  correct50 = placed on the true chromosome with start or end within 50 bp of the true GRCh38 interval
  overlap   = placed on the true chromosome and overlapping the true interval (+-50 bp)

Writes <out>/<arm>.<platform>.reads.npz (per-read arrays aligned with truth.npz rows) and <out>/<arm>.<platform>.json
(overall numbers, MAPQ calibration, strata).
"""
import argparse
import gzip
import os
import re
import sys
import time

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e11b_common as C  # noqa: E402

GRCH38_PATH = {
    "linear": "/scratch/jshen/data/Linear_Reference_VG_Indexes/graph_index/hg38_linear.giraffe.grch38_path",
    "d9": "/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path",
    "d46": "/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.grch38_path",
}
NODE = re.compile(r"[<>](\d+)")
CSOP = re.compile(r"([:*+\-~])([0-9]+|[A-Za-z]+)")
MAPQ_BINS = [(0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 60)]


def log(*x):
    print(time.strftime("%H:%M:%S"), *x, file=sys.stderr, flush=True)


def parse_gaf(path, n_rows, n_pairs, paired, chrom_arr):
    mapped = np.zeros(n_rows, bool)
    mapq = np.full(n_rows, -1, np.int16)
    perfect = np.zeros(n_rows, bool)
    mism = np.full(n_rows, -1, np.int32)
    indel = np.full(n_rows, -1, np.int32)
    clip = np.full(n_rows, -1, np.int32)
    offref_nodes = np.zeros(n_rows, np.int32)
    nodes = np.zeros(n_rows, np.int32)
    seen = np.zeros(n_rows, bool)
    dup = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("@"):  # GAF header lines (vg 1.77 writes "@HD VN:Z:1.0")
                continue
            f = line.split("\t", 12)
            name = f[0]
            if paired:
                base, mate = name.rsplit("/", 1)
                r = (int(mate) - 1) * n_pairs + int(base[1:])
            else:
                r = int(name[1:])
            if seen[r]:
                dup += 1
                continue
            seen[r] = True
            if f[5] == "*":
                continue
            mapped[r] = True
            mapq[r] = int(f[11])
            qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
            cs = ""
            for t in f[12].split("\t"):
                if t.startswith("cs:Z:"):
                    cs = t[5:].rstrip("\n")
                    break
            mm = ind = 0
            for op, val in CSOP.findall(cs):
                if op == "*":
                    mm += 1
                elif op in "+-":
                    ind += len(val)
            mism[r] = mm
            indel[r] = ind
            clip[r] = qlen - (qe - qs)
            perfect[r] = (mm == 0 and ind == 0 and qs == 0 and qe == qlen and cs != "")
            ids = np.fromiter((int(x) for x in NODE.findall(f[5])), np.int64)
            nodes[r] = len(ids)
            if chrom_arr is not None and len(ids):
                ok = ids < len(chrom_arr)
                offref_nodes[r] = int(np.count_nonzero(chrom_arr[ids[ok]] < 0)) + int(np.count_nonzero(~ok))
    return dict(g_seen=seen, g_mapped=mapped, g_mapq=mapq, g_perfect=perfect, g_mism=mism, g_indel=indel,
                g_clip=clip, g_nodes=nodes, g_offref_nodes=offref_nodes), dup


def parse_bam(path, n_rows, n_pairs, paired):
    placed = np.zeros(n_rows, bool)
    chrom = np.full(n_rows, -1, np.int16)
    start = np.full(n_rows, -1, np.int64)
    end = np.full(n_rows, -1, np.int64)
    bmapq = np.full(n_rows, -1, np.int16)
    nm = np.full(n_rows, -1, np.int32)
    seen = np.zeros(n_rows, bool)
    names = {}
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        for i, ref in enumerate(bam.references):
            short = ref.split("#")[-1]
            names[i] = C.CHROMS.index(short) if short in C.CHROMS else 100 + i
        for rec in bam.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            qn = rec.query_name
            if paired:
                base = qn.rsplit("/", 1)[0]
                mate = 1 if rec.is_read1 else 2
                r = (mate - 1) * n_pairs + int(base[1:])
            else:
                r = int(qn[1:])
            seen[r] = True
            if rec.is_unmapped:
                continue
            placed[r] = True
            chrom[r] = names[rec.reference_id]
            start[r] = rec.reference_start
            end[r] = rec.reference_end
            bmapq[r] = rec.mapping_quality
            if rec.has_tag("NM"):
                nm[r] = rec.get_tag("NM")
    return dict(b_seen=seen, b_placed=placed, b_chrom=chrom, b_start=start, b_end=end, b_mapq=bmapq, b_nm=nm)


def rate(num, den):
    return None if den == 0 else round(float(num) / float(den), 6)


def summarize(mask, R, T):
    n = int(mask.sum())
    m = mask
    gm = R["g_mapped"][m]
    pl = R["b_placed"][m]
    cor = R["correct50"][m]
    ov = R["overlap"][m]
    wrong = pl & ~cor
    perf = R["g_perfect"][m]
    mq = R["g_mapq"][m]
    off = R["g_offref_nodes"][m] > 0
    return dict(
        reads=n,
        unmapped=int((~gm).sum()), unmapped_frac=rate((~gm).sum(), n),
        graph_mapped_not_placed=int((gm & ~pl).sum()),
        placed=int(pl.sum()),
        correct50=int(cor.sum()), correct50_frac=rate(cor.sum(), n),
        wrong=int(wrong.sum()), wrong_frac_of_placed=rate(wrong.sum(), pl.sum()),
        wrong_not_overlapping=int((pl & ~ov).sum()),
        perfect=int(perf.sum()), perfect_frac=rate(perf.sum(), n),
        perfect_wrong=int((perf & wrong).sum()), perfect_wrong_frac_of_perfect=rate((perf & wrong).sum(), perf.sum()),
        mapq_ge30=int((mq >= 30).sum()), wrong_mapq_ge30=int((wrong & (mq >= 30)).sum()),
        wrong_frac_mapq_ge30=rate((wrong & (mq >= 30)).sum(), (pl & (mq >= 30)).sum()),
        mapq_eq60=int((mq >= 60).sum()), wrong_mapq60=int((wrong & (mq >= 60)).sum()),
        wrong_frac_mapq60=rate((wrong & (mq >= 60)).sum(), (pl & (mq >= 60)).sum()),
        mapq0_frac=rate(((mq == 0) & gm).sum(), n),
        offref_reads=int(off.sum()), offref_frac=rate(off.sum(), n),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=list(GRCH38_PATH) + ["chm13"])
    ap.add_argument("--platform", required=True, choices=["illumina", "hifi", "ont"])
    ap.add_argument("--base", default="/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    T = dict(np.load(f"{a.base}/sim/{a.platform}/truth.npz"))
    n_rows = len(T["read"])
    paired = a.platform == "illumina"
    n_pairs = n_rows // 2 if paired else n_rows
    if paired:
        assert np.all(T["read"][:n_pairs] == np.arange(n_pairs)) and np.all(T["mate"][:n_pairs] == 1)
    else:
        assert np.all(T["read"] == np.arange(n_rows))
    valid = np.ones(n_rows, bool)
    if a.arm == "chm13":
        # CHM13 control: same reads, truth = the template lifted to CHM13 (chm13_truth.py); reads without a one-to-one
        # CHM13 image are excluded from every count
        Z = np.load(f"{a.base}/chm13_truth/{a.platform}.chm13_truth.npz")
        T["chrom"], T["ref_start"], T["ref_end"] = Z["chrom"], Z["ref_start"], Z["ref_end"]
        valid = Z["ok"].astype(bool)
        chrom_arr = None
    else:
        chrom_arr = np.load(f"{GRCH38_PATH[a.arm]}/chrom.npy")
    log("truth", n_rows, "valid", int(valid.sum()), "chrom index", None if chrom_arr is None else len(chrom_arr))
    adir = f"{a.base}/aln/{a.arm}"
    G, dup = parse_gaf(f"{adir}/{a.platform}.gaf.gz", n_rows, n_pairs, paired, chrom_arr)
    log("gaf parsed; duplicate records skipped", dup)
    B = parse_bam(f"{adir}/{a.platform}.bam", n_rows, n_pairs, paired)
    log("bam parsed")
    R = {**G, **B}
    same = R["b_placed"] & (R["b_chrom"] == T["chrom"])
    d = np.minimum(np.abs(R["b_start"] - T["ref_start"]), np.abs(R["b_end"] - T["ref_end"]))
    R["correct50"] = same & (d <= 50)
    R["overlap"] = same & (R["b_start"] < T["ref_end"] + 50) & (R["b_end"] > T["ref_start"] - 50)
    R["dist"] = np.where(same, d, -1)
    np.savez_compressed(f"{a.out}/{a.arm}.{a.platform}.reads.npz", **R)

    res = dict(arm=a.arm, platform=a.platform, reads=n_rows, gaf_records_seen=int(G["g_seen"].sum()),
               bam_records_seen=int(B["b_seen"].sum()), gaf_duplicates=dup,
               mapq_gaf_eq_bam=rate((R["g_mapq"] == R["b_mapq"])[R["b_placed"]].sum(), R["b_placed"].sum()))
    res["valid_reads"] = int(valid.sum())
    res["overall"] = summarize(valid, R, T)
    if paired:
        res["error_free_reads"] = summarize(valid & (T["n_err"] == 0), R, T)
    # MAPQ calibration (graph MAPQ, reads placed on GRCh38)
    cal = []
    pl = R["b_placed"] & valid
    for lo, hi in MAPQ_BINS:
        m = pl & (R["g_mapq"] >= lo) & (R["g_mapq"] <= hi)
        w = m & ~R["correct50"]
        exp = float(np.mean(10.0 ** (-R["g_mapq"][m] / 10.0))) if m.any() else None
        cal.append(dict(mapq=f"{lo}-{hi}" if lo != hi else str(lo), placed=int(m.sum()), wrong=int(w.sum()),
                        error_rate=rate(w.sum(), m.sum()), expected_from_mapq=None if exp is None else round(exp, 6)))
    res["mapq_calibration"] = cal
    # strata
    nv = T["n_snv"] + T["n_indel"] + T["n_sv"]
    strata = {
        "hg008n_variants_in_read": {
            "0": nv == 0, "1": nv == 1, "2-3": (nv >= 2) & (nv <= 3), ">=4": nv >= 4},
        "sv_or_insertion": {
            "read spans an SV (>=50 bp) of its haplotype": T["n_sv"] > 0,
            ">=10 read bases absent from GRCh38": T["inserted_bases"] >= 10,
            "neither": (T["n_sv"] == 0) & (T["inserted_bases"] < 10)},
        "distance_to_nearest_hg008n_variant": {
            "inside read": T["var_dist"] == 0, "1-50 bp": (T["var_dist"] >= 1) & (T["var_dist"] <= 50),
            "51-500 bp": (T["var_dist"] > 50) & (T["var_dist"] <= 500), ">500 bp": T["var_dist"] > 500},
        "smaht_region": {"easy": T["smaht"] == 0, "difficult": T["smaht"] == 1, "extreme": T["smaht"] == 2},
    }
    strf = T["str_bases"] / T["tmpl_len"]
    if paired:
        strata["str_fraction_of_read"] = {"0": strf == 0, "(0, 0.5)": (strf > 0) & (strf < 0.5), ">=0.5": strf >= 0.5}
        strata["homopolymer_ge7"] = {"yes": T["hp7_bases"] > 0, "no": T["hp7_bases"] == 0}
    else:
        strata["str_fraction_of_read"] = {"<1%": strf < 0.01, "1-5%": (strf >= 0.01) & (strf < 0.05),
                                          ">=5%": strf >= 0.05}
    if a.arm not in ("linear", "chm13"):
        strata["graph_alignment_uses_offref_node"] = {"yes": G["g_offref_nodes"] > 0,
                                                      "no": G["g_mapped"] & (G["g_offref_nodes"] == 0)}
    res["strata"] = {k: {s: summarize(m & valid, R, T) for s, m in v.items()} for k, v in strata.items()}
    C.write_json(f"{a.out}/{a.arm}.{a.platform}.json", res)
    o = res["overall"]
    log(a.arm, a.platform, "correct50", o["correct50_frac"], "wrong/placed", o["wrong_frac_of_placed"],
        "perfect", o["perfect_frac"], "perfect_wrong", o["perfect_wrong"], "unmapped", o["unmapped"])


if __name__ == "__main__":
    main()
