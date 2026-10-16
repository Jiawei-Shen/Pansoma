#!/usr/bin/env python3
"""E11b step 2: simulate reads of known GRCh38 origin from the HG008-N v6.2 haplotype segments (build_segments.py).

  illumina : 2 x 150 pairs, fragment ~ N(400, 50) clipped to [250, 600], uniform substitutions (--sub, 0.2 %)
  hifi     : single reads, length ~ N(15 kb, 3 kb) clipped to [5, 30] kb, 0.1 % errors (0.04 % sub, 0.03 % ins,
             0.03 % del, 1 bp each)
  ont      : templates as for hifi (length ~ N(15 kb, 4 kb) clipped to [3, 30] kb); errors added by pbsim3
             (--strategy templ, ERRHMM-ONT-HQ model); qualities replaced by a constant Q20

Templates are drawn uniformly over the segment bases of both haplotypes (no template with a non-ACGT base),
strand 50/50. Every read's truth is its template's GRCh38 interval through dipcall's hap -> GRCh38 alignment.

Outputs (in --out/<platform>/): reads (*.fq.gz), truth.npz (one row per read; mate 1/2 for Illumina, 0 otherwise)
"""
import argparse
import gzip
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e11b_common as C  # noqa: E402

PBSIM = "/scratch/wzhang/bin/pbsim"
ONT_MODEL = "/scratch/wzhang/repo/pbsim3/data/ERRHMM-ONT-HQ.model"
PIGZ = "/scratch/jshen/bin/pigz"


def log(*x):
    print(time.strftime("%H:%M:%S"), *x, file=sys.stderr, flush=True)


def draw_templates(rng, segs, n, lengths):
    """lengths: array of n template lengths -> (seg index, start) with the template inside the segment."""
    seg_len = np.array([s["qend"] - s["qstart"] for s in segs], np.int64)
    cum = np.concatenate([[0], np.cumsum(seg_len)])
    out_seg = np.empty(n, np.int64)
    out_start = np.empty(n, np.int64)
    todo = np.arange(n)
    while len(todo):
        g = rng.integers(0, cum[-1], len(todo))
        si = np.searchsorted(cum, g, side="right") - 1
        st = g - cum[si]
        ok = st + lengths[todo] <= seg_len[si]
        out_seg[todo[ok]] = si[ok]
        out_start[todo[ok]] = st[ok]
        todo = todo[~ok]
    return out_seg, out_start


def get_seqs(contigs, segs, seg_idx, qa, L):
    """template code arrays (list) in hap orientation"""
    out = [None] * len(seg_idx)
    for i in range(len(seg_idx)):
        c = contigs[(segs[seg_idx[i]]["hap"], segs[seg_idx[i]]["contig"])]
        out[i] = c[qa[i]:qa[i] + L[i]]
    return out


def revcomp(codes):
    return C.COMP[codes][::-1]


def features(al, cums, vars_by_hap, rep, hap, contig_key, aln, qa, qb, strand):
    chrom, rs, re_, rstr, aligned = al.lift_reads(aln, qa, qb, strand)
    cls = C.smaht_class(cums, chrom, rs, re_)
    nsnv, nind, nsv, dist = C.variant_features(vars_by_hap, hap, chrom, rs, re_)
    strb = np.zeros(len(qa), np.int64)
    hpb = np.zeros(len(qa), np.int64)
    for key in np.unique(contig_key):
        m = contig_key == key
        cs_str, cs_hp = rep[key]
        strb[m] = cs_str[qb[m]] - cs_str[qa[m]]
        hpb[m] = cs_hp[qb[m]] - cs_hp[qa[m]]
    return dict(chrom=chrom, ref_start=rs, ref_end=re_, ref_strand=rstr, aligned_bases=aligned,
                inserted_bases=(qb - qa) - aligned, smaht=cls, n_snv=nsnv, n_indel=nind, n_sv=nsv, var_dist=dist,
                str_bases=strb, hp7_bases=hpb)


def fastq_writer(path):
    p = subprocess.Popen([PIGZ, "-p", "4", "-c"], stdin=subprocess.PIPE, stdout=open(path, "wb"))
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", required=True, help="dir with segments.tsv / alignments.*")
    ap.add_argument("--platform", choices=["illumina", "hifi", "ont"], required=True)
    ap.add_argument("--n", type=int, required=True, help="pairs (illumina) or reads")
    ap.add_argument("--seed", type=int, default=20261013)
    ap.add_argument("--sub", type=float, default=0.002)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.out, exist_ok=True)

    al = C.Alignments(a.sim)
    segs = al.segs
    log("segments", len(segs))
    contigs = {}
    rep = {}
    for s in segs:
        key = (s["hap"], s["contig"])
        if key not in contigs:
            contigs[key] = C.load_contig(*key)
            rep[key] = C.repeat_cumsums(contigs[key], f"{a.sim}/repeats_hap{key[0]}_{key[1]}.npz")
            log("loaded", key, len(contigs[key]))
    keys = list(contigs.keys())
    cums = C.smaht_cumsums()
    vars_by_hap = {h: C.load_hap_variants(h) for h in (1, 2)}
    log("strata loaded")

    n = a.n
    if a.platform == "illumina":
        L = np.clip(np.rint(rng.normal(400, 50, n)), 250, 600).astype(np.int64)
    elif a.platform == "hifi":
        L = np.clip(np.rint(rng.normal(15000, 3000, n)), 5000, 30000).astype(np.int64)
    else:
        L = np.clip(np.rint(rng.normal(15000, 4000, n)), 3000, 30000).astype(np.int64)

    # draw until every template is pure ACGT
    seg_idx, start = draw_templates(rng, segs, n, L)
    for _ in range(50):
        qa = np.array([segs[i]["qstart"] for i in seg_idx]) + start
        bad = np.array([np.any(contigs[(segs[seg_idx[i]]["hap"], segs[seg_idx[i]]["contig"])][qa[i]:qa[i] + L[i]] > 3)
                        for i in range(n)])
        if not bad.any():
            break
        s2, st2 = draw_templates(rng, segs, int(bad.sum()), L[bad])
        seg_idx[bad], start[bad] = s2, st2
    qa = np.array([segs[i]["qstart"] for i in seg_idx]) + start
    frag_strand = np.where(rng.random(n) < 0.5, 1, -1).astype(np.int8)
    hap = np.array([segs[i]["hap"] for i in seg_idx], np.int8)
    aln = np.array([segs[i]["aln"] for i in seg_idx], np.int64)
    ckey = np.array([keys.index((segs[i]["hap"], segs[i]["contig"])) for i in seg_idx], np.int64)
    rep_by_idx = {keys.index(k): v for k, v in rep.items()}
    log("templates drawn", n)

    rows = {}

    def add_rows(d):
        for k, v in d.items():
            rows.setdefault(k, []).append(np.asarray(v))

    if a.platform == "illumina":
        RL = 150
        # read 1 = fragment start in fragment orientation; read 2 = the other end, opposite strand
        r1_qa = np.where(frag_strand == 1, qa, qa + L - RL)
        r2_qa = np.where(frag_strand == 1, qa + L - RL, qa)
        names = np.arange(n)
        f1 = fastq_writer(f"{a.out}/sim_R1.fq.gz")
        f2 = fastq_writer(f"{a.out}/sim_R2.fq.gz")
        qual = b"F" * RL
        nerr = np.zeros((2, n), np.int32)
        for mate, rqa, rstr, fh in ((1, r1_qa, frag_strand, f1), (2, r2_qa, -frag_strand, f2)):
            buf = []
            for i in range(n):
                c = contigs[keys[ckey[i]]][rqa[i]:rqa[i] + RL]
                if rstr[i] < 0:
                    c = C.COMP[c][::-1]
                c = c.copy()
                err = np.flatnonzero(rng.random(RL) < a.sub)
                if len(err):
                    c[err] = (c[err] + rng.integers(1, 4, len(err))) % 4
                nerr[mate - 1, i] = len(err)
                buf.append(b"@i%d/%d\n%s\n+\n%s\n" % (i, mate, C.BASES[c].tobytes(), qual))
                if len(buf) >= 200000:
                    fh.stdin.write(b"".join(buf))
                    buf = []
            fh.stdin.write(b"".join(buf))
            fh.stdin.close()
            fh.wait()
            ft = features(al, cums, vars_by_hap, rep_by_idx, hap, ckey, aln, rqa, rqa + RL, rstr)
            add_rows(dict(read=names, mate=np.full(n, mate, np.int8), hap=hap, seg=seg_idx, frag_len=L,
                          tmpl_len=np.full(n, RL), n_err=nerr[mate - 1], hap_start=rqa, hap_end=rqa + RL,
                          hap_strand=rstr, **ft))
            log("mate", mate, "written")
    else:
        names = np.arange(n)
        tmpl_strand = frag_strand
        seqs = []
        for i in range(n):
            c = contigs[keys[ckey[i]]][qa[i]:qa[i] + L[i]]
            if tmpl_strand[i] < 0:
                c = C.COMP[c][::-1]
            seqs.append(c)
        nerr = np.zeros(n, np.int32)
        if a.platform == "hifi":
            fh = fastq_writer(f"{a.out}/sim.fq.gz")
            buf = []
            for i in range(n):
                c = seqs[i]
                m = len(c)
                ne = rng.binomial(m, 0.001)
                pos = np.sort(rng.choice(m, ne, replace=False))[::-1]
                kinds = rng.choice(3, ne, p=[0.4, 0.3, 0.3])
                c = list(c)
                for p_, k in zip(pos, kinds):
                    if k == 0:
                        c[p_] = (c[p_] + rng.integers(1, 4)) % 4
                    elif k == 1:
                        c.insert(p_, int(rng.integers(0, 4)))
                    else:
                        del c[p_]
                c = np.array(c, np.uint8)
                nerr[i] = ne
                buf.append(b"@h%d\n%s\n+\n%s\n" % (i, C.BASES[c].tobytes(), b"F" * len(c)))
            fh.stdin.write(b"".join(buf))
            fh.stdin.close()
            fh.wait()
        else:
            tdir = f"{a.out}/pbsim"
            os.makedirs(tdir, exist_ok=True)
            with open(f"{tdir}/templates.fa", "wb") as fh:
                for i in range(n):
                    fh.write(b">o%d\n%s\n" % (i, C.BASES[seqs[i]].tobytes()))
            cmd = [PBSIM, "--strategy", "templ", "--template", f"{tdir}/templates.fa", "--method", "errhmm",
                   "--errhmm", ONT_MODEL, "--accuracy-mean", "0.98", "--prefix", f"{tdir}/ont", "--seed",
                   str(a.seed)]
            log(" ".join(cmd))
            subprocess.run(cmd, check=True, stdout=open(f"{tdir}/pbsim.log", "w"), stderr=subprocess.STDOUT)
            # pbsim3 template mode: read S_k = template k (1-based, file order), forward strand of the template
            fh = fastq_writer(f"{a.out}/sim.fq.gz")
            buf = []
            k = 0
            with open(f"{tdir}/ont.fastq", "rb") as fq:
                while True:
                    h = fq.readline()
                    if not h:
                        break
                    s = fq.readline().rstrip(b"\n")
                    fq.readline()
                    fq.readline()
                    assert h.strip() == b"@S_%d" % (k + 1), (h, k)
                    buf.append(b"@o%d\n%s\n+\n%s\n" % (k, s, b"5" * len(s)))
                    k += 1
            assert k == n, (k, n)
            fh.stdin.write(b"".join(buf))
            fh.stdin.close()
            fh.wait()
            os.remove(f"{tdir}/ont.fastq")
            nerr[:] = -1
        ft = features(al, cums, vars_by_hap, rep_by_idx, hap, ckey, aln, qa, qa + L, tmpl_strand)
        add_rows(dict(read=names, mate=np.zeros(n, np.int8), hap=hap, seg=seg_idx, frag_len=L, tmpl_len=L,
                      n_err=nerr, hap_start=qa, hap_end=qa + L, hap_strand=tmpl_strand, **ft))
    truth = {k: np.concatenate(v) for k, v in rows.items()}
    np.savez_compressed(f"{a.out}/truth.npz", **truth)
    C.write_json(f"{a.out}/sim_meta.json", dict(vars(a), reads=int(len(truth["read"]))))
    log("done", len(truth["read"]))


if __name__ == "__main__":
    main()
