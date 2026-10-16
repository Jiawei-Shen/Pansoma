#!/usr/bin/env python3
"""E11b step 1: one-to-one HG008-N haplotype segments on GRCh38 chr1 / chr20, with a base-level hap -> GRCh38 liftover.

Source of truth = dipcall's minimap2 -x asm5 alignments of the curated HG008-N v6.2 haplotypes to GRCh38
(dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.hap{1,2}.paf.gz, cg:Z CIGAR) and dipcall's per-haplotype
"covered by exactly one contig" regions (HG008N_GRCh38_dipcall.hap{1,2}.bed = paftools call R lines).

A segment is a maximal haplotype interval such that
  * it lies in one primary (tp:A:P) alignment of chrN_hapK -> chrN (N in --chroms), MAPQ >= 5, block length >= 50 kb;
  * its GRCh38 image lies inside one interval of hapK.bed (reference side covered once);
  * no haplotype base is covered by a second primary / inversion alignment (query side covered once).

Outputs (in --out):
  segments.tsv         seg, hap, contig, qstart, qend, aln, chrom, tstart, tend, strand
  alignments.npz       per alignment: M-block arrays (query offset along the aligned strand, ref start, length)
  segments_meta.json   inputs, filters, totals
"""
import argparse
import gzip
import json
import re
import sys

import numpy as np

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def read_bed(path, chroms):
    out = {c: [] for c in chroms}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] in out:
                out[f[0]].append((int(f[1]), int(f[2])))
    return {c: np.array(sorted(v), dtype=np.int64).reshape(-1, 2) for c, v in out.items()}


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def multi_covered(iv):
    """Intervals covered by >= 2 of the given intervals."""
    ev = []
    for s, e in iv:
        ev.append((s, 1))
        ev.append((e, -1))
    ev.sort()
    out, depth, start = [], 0, None
    for p, d in ev:
        nd = depth + d
        if depth < 2 <= nd:
            start = p
        elif depth >= 2 > nd:
            if p > start:
                out.append([start, p])
        depth = nd
    return merge(out)


def subtract(seg, bad):
    """seg = (s, e); bad = sorted merged list; returns list of pieces."""
    s, e = seg
    out = []
    cur = s
    for bs, be in bad:
        if be <= cur:
            continue
        if bs >= e:
            break
        if bs > cur:
            out.append((cur, bs))
        cur = max(cur, be)
        if cur >= e:
            break
    if cur < e:
        out.append((cur, e))
    return out


def cigar_blocks(cg, tstart):
    """M blocks: qoff (offset along the aligned query strand), tpos (absolute ref), length."""
    q = 0
    t = tstart
    qo, tp, ln = [], [], []
    for n, op in CIG.findall(cg):
        n = int(n)
        if op in "M=X":
            qo.append(q)
            tp.append(t)
            ln.append(n)
            q += n
            t += n
        elif op == "I":
            q += n
        elif op in "DN":
            t += n
        else:
            raise ValueError("unexpected CIGAR op " + op)
    return np.array(qo, np.int64), np.array(tp, np.int64), np.array(ln, np.int64), q, t


def ref_to_qoff_lo(qo, tp, ln, t):
    """Smallest query offset of an aligned base with ref >= t."""
    k = np.searchsorted(tp + ln, t, side="right")  # first block whose end > t
    if k >= len(tp):
        return None
    if tp[k] >= t:
        return int(qo[k])
    return int(qo[k] + (t - tp[k]))


def ref_to_qoff_hi(qo, tp, ln, t):
    """1 + largest query offset of an aligned base with ref < t."""
    k = np.searchsorted(tp, t, side="left") - 1  # last block starting before t
    if k < 0:
        return None
    if tp[k] + ln[k] <= t:
        return int(qo[k] + ln[k])
    return int(qo[k] + (t - tp[k]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dipcall-prefix", default="/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall")
    ap.add_argument("--chroms", default="chr1,chr20")
    ap.add_argument("--min-aln", type=int, default=50000)
    ap.add_argument("--min-mapq", type=int, default=5)
    ap.add_argument("--min-seg", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    chroms = a.chroms.split(",")

    seg_rows = []
    aln_arrays = {}
    aln_rows = []
    meta = {"inputs": {}, "filters": vars(a), "per_hap": {}}
    for hap in (1, 2):
        paf = f"{a.dipcall_prefix}.hap{hap}.paf.gz"
        bedp = f"{a.dipcall_prefix}.hap{hap}.bed"
        meta["inputs"][f"hap{hap}"] = {"paf": paf, "bed": bedp}
        bed = read_bed(bedp, chroms)
        recs = []
        qcov = {}
        with gzip.open(paf, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[5] == "*":
                    continue
                tags = dict((x[:2], x[5:]) for x in f[12:])
                tp = tags.get("tp", "")
                qname, qs, qe = f[0], int(f[2]), int(f[3])
                if tp in ("P", "I", "i"):
                    qcov.setdefault(qname, []).append((qs, qe))
                if tp != "P" or f[5] not in chroms:
                    continue
                if qname != f"{f[5]}_hap{hap}":
                    continue
                if int(f[11]) < a.min_mapq or int(f[10]) < a.min_aln:
                    continue
                recs.append((qname, int(f[1]), qs, qe, f[4], f[5], int(f[7]), int(f[8]), int(f[11]), tags["cg"]))
        hap_bases = 0
        for r in recs:
            qname, qlen, qs, qe, strand, chrom, ts, te, mapq, cg = r
            qo, tp_, ln, qspan, tend = cigar_blocks(cg, ts)
            assert qspan == qe - qs and tend == te, (r[:9], qspan, tend)
            aid = len(aln_rows)
            aln_rows.append(dict(aln=aid, hap=hap, contig=qname, qstart=qs, qend=qe, strand=strand, chrom=chrom,
                                 tstart=ts, tend=te, mapq=mapq))
            aln_arrays[f"a{aid}_qo"] = qo
            aln_arrays[f"a{aid}_tp"] = tp_
            aln_arrays[f"a{aid}_ln"] = ln
            bad = multi_covered(qcov.get(qname, []))
            b = bed[chrom]
            sel = b[(b[:, 1] > ts) & (b[:, 0] < te)]
            for bs, be in sel:
                bs, be = max(int(bs), ts), min(int(be), te)
                lo = ref_to_qoff_lo(qo, tp_, ln, bs)
                hi = ref_to_qoff_hi(qo, tp_, ln, be)
                if lo is None or hi is None or hi - lo < a.min_seg:
                    continue
                if strand == "+":
                    q0, q1 = qs + lo, qs + hi
                else:
                    q0, q1 = qe - hi, qe - lo
                for p0, p1 in subtract((q0, q1), bad):
                    if p1 - p0 < a.min_seg:
                        continue
                    seg_rows.append((len(seg_rows), hap, qname, p0, p1, aid, chrom, bs, be, strand))
                    hap_bases += p1 - p0
        meta["per_hap"][f"hap{hap}"] = {"alignments": len([x for x in aln_rows if x["hap"] == hap]),
                                       "segment_bases": hap_bases}
        print(f"hap{hap}: {len(recs)} alignments, {hap_bases:,} segment bases", file=sys.stderr)

    with open(f"{a.out}/segments.tsv", "w") as fh:
        fh.write("seg\thap\tcontig\tqstart\tqend\taln\tchrom\tbed_start\tbed_end\tstrand\n")
        for row in seg_rows:
            fh.write("\t".join(map(str, row)) + "\n")
    with open(f"{a.out}/alignments.tsv", "w") as fh:
        keys = list(aln_rows[0].keys())
        fh.write("\t".join(keys) + "\n")
        for r in aln_rows:
            fh.write("\t".join(str(r[k]) for k in keys) + "\n")
    np.savez(f"{a.out}/alignments.npz", **aln_arrays)
    meta["segments"] = len(seg_rows)
    meta["segment_bases"] = int(sum(r[4] - r[3] for r in seg_rows))
    with open(f"{a.out}/segments_meta.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps({k: meta[k] for k in ("segments", "segment_bases", "per_hap")}), file=sys.stderr)


if __name__ == "__main__":
    main()
