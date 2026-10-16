"""Shared helpers for E11b: haplotype segments, hap -> GRCh38 liftover, read-level strata."""
import gzip
import json
import os

import numpy as np

DIPCALL = "/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall"
HAP_FASTA = {1: "/scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_polished6.2_hap1.fasta",
             2: "/scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_polished6.2_hap2.fasta.gz"}
SMAHT = {1: "/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_difficult_v2.bed",
         2: "/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_extreme_v2.bed"}  # easy (0) = the rest
CHROMS = ["chr1", "chr20"]
REF_LEN = {"chr1": 248956422, "chr20": 64444167}

# ACGT -> 0..3, everything else -> 4
CODE = np.full(256, 4, np.uint8)
for i, c in enumerate(b"ACGT"):
    CODE[c] = i
    CODE[c + 32] = i
BASES = np.frombuffer(b"ACGTN", np.uint8)
COMP = np.array([3, 2, 1, 0, 4], np.uint8)


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return [dict(zip(head, r)) for r in rows]


class Alignments:
    """dipcall hap -> GRCh38 alignments (M blocks) and the one-to-one segments built by build_segments.py."""

    def __init__(self, simdir):
        self.alns = read_tsv(f"{simdir}/alignments.tsv")
        for r in self.alns:
            for k in ("aln", "hap", "qstart", "qend", "tstart", "tend", "mapq"):
                r[k] = int(r[k])
        z = np.load(f"{simdir}/alignments.npz")
        self.qo, self.tp, self.ln, self.cum = [], [], [], []
        for r in self.alns:
            a = r["aln"]
            qo, tp, ln = z[f"a{a}_qo"], z[f"a{a}_tp"], z[f"a{a}_ln"]
            self.qo.append(qo)
            self.tp.append(tp)
            self.ln.append(ln)
            self.cum.append(np.concatenate([[0], np.cumsum(ln)[:-1]]))
        self.segs = read_tsv(f"{simdir}/segments.tsv")
        for r in self.segs:
            for k in ("seg", "hap", "qstart", "qend", "aln", "bed_start", "bed_end"):
                r[k] = int(r[k])

    def _qoff(self, a, q):
        r = self.alns[a]
        return (q - r["qstart"]) if r["strand"] == "+" else (r["qend"] - 1 - q)

    def _lift_off(self, a, off):
        qo, tp, ln = self.qo[a], self.tp[a], self.ln[a]
        k = np.searchsorted(qo, off, side="right") - 1
        k = np.clip(k, 0, len(qo) - 1)
        inb = (off >= qo[k]) & (off < qo[k] + ln[k])
        t = np.where(inb, tp[k] + (off - qo[k]), tp[k] + ln[k])
        return t, inb

    def _cum_m(self, a, x):
        """number of aligned (M) bases with query offset < x"""
        qo, ln, cum = self.qo[a], self.ln[a], self.cum[a]
        k = np.searchsorted(qo, x, side="right") - 1
        kk = np.clip(k, 0, len(qo) - 1)
        v = cum[kk] + np.minimum(np.maximum(x - qo[kk], 0), ln[kk])
        return np.where(k < 0, 0, v)

    def lift_reads(self, aln, qa, qb, hap_strand):
        """Template [qa, qb) on the haplotype (hap coords), read strand vs hap (+1/-1) ->
        GRCh38 chrom index, ref_start, ref_end (0-based, half open), ref strand, aligned bases."""
        n = len(aln)
        ref_start = np.zeros(n, np.int64)
        ref_end = np.zeros(n, np.int64)
        ref_strand = np.zeros(n, np.int8)
        aligned = np.zeros(n, np.int64)
        chrom = np.zeros(n, np.int8)
        for a in np.unique(aln):
            m = aln == a
            r = self.alns[a]
            o_a = self._qoff(a, qa[m])
            o_b = self._qoff(a, qb[m] - 1)
            lo, hi = np.minimum(o_a, o_b), np.maximum(o_a, o_b)
            t_lo, in_lo = self._lift_off(a, lo)
            t_hi, in_hi = self._lift_off(a, hi)
            s = t_lo
            e = np.where(in_hi, t_hi + 1, t_hi)
            e = np.maximum(e, s)
            ref_start[m] = s
            ref_end[m] = e
            aligned[m] = self._cum_m(a, hi + 1) - self._cum_m(a, lo)
            sgn = 1 if r["strand"] == "+" else -1
            ref_strand[m] = hap_strand[m] * sgn
            chrom[m] = CHROMS.index(r["chrom"])
        return chrom, ref_start, ref_end, ref_strand, aligned


def load_contig(hap, contig):
    import pysam
    with pysam.FastaFile(HAP_FASTA[hap]) as fa:
        s = fa.fetch(contig)
    return CODE[np.frombuffer(s.encode(), np.uint8)]


def repeat_masks(codes, min_str=12, max_period=6, hp_min=7):
    """cumsum arrays (len n+1, int32) of bases in short tandem repeats (period 1-6, tract >= max(12, 3 x period))
    and of bases in homopolymers >= 7 bp."""
    n = len(codes)
    strm = np.zeros(n, bool)
    hpm = np.zeros(n, bool)
    for p in range(1, max_period + 1):
        eq = (codes[p:] == codes[:-p]) & (codes[p:] < 4)
        d = np.diff(np.concatenate(([False], eq, [False])).view(np.int8))
        st = np.flatnonzero(d == 1)
        en = np.flatnonzero(d == -1)  # eq[st:en] true -> tract [st, en + p)
        tl = en - st + p
        for need, mask in ((max(min_str, 3 * p), strm),) + (((hp_min, hpm),) if p == 1 else ()):
            ok = tl >= need
            cov = (np.bincount(st[ok], minlength=n + 1).astype(np.int32)
                   - np.bincount(np.minimum(en[ok] + p, n), minlength=n + 1).astype(np.int32))
            mask |= np.cumsum(cov[:-1], dtype=np.int32) > 0
    return strm, hpm


def repeat_cumsums(codes, cache):
    """repeat_masks() as cumsum arrays (len n+1), cached as packed bits in `cache` (.npz)"""
    n = len(codes)
    if os.path.exists(cache):
        z = np.load(cache)
        strm = np.unpackbits(z["str"], count=n).astype(bool)
        hpm = np.unpackbits(z["hp"], count=n).astype(bool)
    else:
        strm, hpm = repeat_masks(codes)
        np.savez(cache + f".tmp{os.getpid()}.npz", str=np.packbits(strm), hp=np.packbits(hpm))
        os.replace(cache + f".tmp{os.getpid()}.npz", cache)
    return (np.concatenate([[0], np.cumsum(strm, dtype=np.int32)]),
            np.concatenate([[0], np.cumsum(hpm, dtype=np.int32)]))


SMAHT_CACHE = "/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/sim/smaht_chr1_chr20.npz"


def smaht_cumsums(cache=SMAHT_CACHE):
    """per chrom: (cum_difficult, cum_extreme) int32 arrays of length L+1 on GRCh38 (codes cached as packed bits)"""
    codes = {}
    if os.path.exists(cache):
        z = np.load(cache)
        for c in CHROMS:
            codes[c] = (np.unpackbits(z[f"{c}_d"], count=REF_LEN[c]).astype(bool),
                        np.unpackbits(z[f"{c}_x"], count=REF_LEN[c]).astype(bool))
    else:
        for c in CHROMS:
            code = np.zeros(REF_LEN[c], np.int8)
            for k, path in SMAHT.items():
                with open(path) as fh:
                    for line in fh:
                        f = line.split("\t", 3)
                        if f[0] == c:
                            code[int(f[1]):int(f[2])] = k
            codes[c] = (code == 1, code == 2)
        np.savez(cache + f".tmp{os.getpid()}.npz", **{f"{c}_{t}": np.packbits(codes[c][i]) for c in CHROMS
                                        for i, t in enumerate("dx")})
        os.replace(cache + f".tmp{os.getpid()}.npz", cache)
    return {c: (np.concatenate([[0], np.cumsum(codes[c][0], dtype=np.int32)]),
                np.concatenate([[0], np.cumsum(codes[c][1], dtype=np.int32)])) for c in CHROMS}


def smaht_class(cums, chrom, s, e):
    """0 easy, 1 difficult, 2 extreme by majority of the read's GRCh38 bases (ties -> harder);
    zero-length intervals use the flanking base."""
    cls = np.zeros(len(s), np.int8)
    for ci, c in enumerate(CHROMS):
        m = chrom == ci
        if not m.any():
            continue
        cd, cx = cums[c]
        ss = s[m]
        ee = np.maximum(e[m], ss + 1)
        ee = np.minimum(ee, len(cd) - 1)
        ss = np.minimum(ss, ee - 1)
        L = (ee - ss).astype(np.float64)
        fx = (cx[ee] - cx[ss]) / L
        fd = (cd[ee] - cd[ss]) / L
        cls[m] = np.where(fx >= 0.5, 2, np.where(fx + fd >= 0.5, 1, 0))
    return cls


def load_hap_variants(hap):
    """paftools call V lines of hapN.var.gz (the per-haplotype input of dipcall's VCF) on CHROMS."""
    v = {c: [] for c in CHROMS}
    with gzip.open(f"{DIPCALL}.hap{hap}.var.gz", "rt") as fh:
        for line in fh:
            if line[0] != "V":
                continue
            f = line.split("\t")
            if f[1] not in v:
                continue
            s, e = int(f[2]), int(f[3])
            ref, alt = f[6], f[7]
            lr = 0 if ref == "-" else len(ref)
            la = 0 if alt == "-" else len(alt)
            kind = 0 if (lr == 1 and la == 1) else (2 if max(lr, la) >= 50 else 1)
            v[f[1]].append((s, max(e, s + 1), kind))
    out = {}
    for c in CHROMS:
        a = np.array(sorted(v[c]), np.int64).reshape(-1, 3)
        out[c] = a
    return out


def variant_features(vars_by_hap, hap, chrom, s, e):
    """per read: counts of the read's own haplotype variants (SNV, indel < 50, SV >= 50) inside the true GRCh38
    interval, and the distance to the nearest one (0 = inside)."""
    n = len(s)
    nsnv = np.zeros(n, np.int32)
    nind = np.zeros(n, np.int32)
    nsv = np.zeros(n, np.int32)
    dist = np.full(n, 10 ** 9, np.int64)
    for h in (1, 2):
        for ci, c in enumerate(CHROMS):
            m = (hap == h) & (chrom == ci)
            if not m.any():
                continue
            v = vars_by_hap[h][c]
            vs, ve, vk = v[:, 0], v[:, 1], v[:, 2]
            ss = s[m]
            ee = np.maximum(e[m], ss + 1)
            hi = np.searchsorted(vs, ee, side="left")       # variants starting before read end
            endmax = np.maximum.accumulate(ve)
            lo = np.searchsorted(endmax, ss, side="right")   # first variant whose (running max) end > read start
            idx = np.flatnonzero(hi > lo)
            cs = [np.concatenate([[0], np.cumsum(vk == k)]) for k in (0, 1, 2)]
            out = [np.zeros(len(ss), np.int32) for _ in range(3)]
            for k in range(3):
                out[k] = (cs[k][hi] - cs[k][np.minimum(lo, hi)]).astype(np.int32)
            nsnv[m], nind[m], nsv[m] = out
            d = np.full(len(ss), 10 ** 9, np.int64)
            d[idx] = 0
            # nearest before / after
            prev = lo - 1
            okp = prev >= 0
            dp = np.where(okp, ss - endmax[np.maximum(prev, 0)] + 1, 10 ** 9)
            okn = hi < len(vs)
            dn = np.where(okn, vs[np.minimum(hi, len(vs) - 1)] - ee + 1, 10 ** 9)
            d = np.where(d == 0, 0, np.minimum(dp, dn))
            dist[m] = d
    return nsnv, nind, nsv, dist


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1)
