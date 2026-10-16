"""E11a: the GRCh38 BAM arm. One region-parallel pass over a coordinate-sorted, indexed BAM; per-read statistics for
the same name-hash sample as scan_gam.py (zlib.crc32(read name) % MOD < KEEP), from the PRIMARY record of each read
(secondary records are counted; supplementary records are counted and their aligned query bases summed per read).

Each record is read once: a region job keeps the records whose reference_start lies in its interval; the unplaced
unmapped records come from fetch('*').
Per read (same columns as scan_gam.py; the graph-only columns off_* stay 0):
  mate 0/1/2 from flags 0x40/0x80; read_len = query length incl. soft and hard clips; mapped = not flag 0x4
  match / mm: '='/'X' CIGAR bases when the CIGAR uses them, else M bases split by the MD tag (or NM - I - D);
  mm_ev = runs of consecutive mismatched reference positions; ins/ins_ev, dl/del_ev from I/D operations;
  sc_l / sc_r: S and H at the alignment ends; chrom (index into the BAM's references), gstart/gend = reference
  start/end of the primary record; nchrom = 1 + supplementary records on other contigs (0 when unmapped)
Extra columns: supp (supplementary records), supp_bp (their aligned query bases), supp_other_chrom.
Output: <out>.npz (rows, contigs) and <out>.summary.json (record counts per flag class over the sampled reads, and
MAPQ histogram of all primary records of the file).

Run: python scan_bam.py --bam B --out PREFIX [--mod 1000 --keep 1] [--procs 16]
"""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
import os
import re
import sys
import time
import zlib

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scan_gam import ROW  # noqa: E402

BROW = np.dtype(ROW.descr + [('supp', 'u2'), ('supp_bp', 'u4'), ('supp_other_chrom', 'u2')])
MD_RE = re.compile(r'(\d+)|(\^[A-Za-z]+)|([A-Za-z])')
CHUNK = 20_000_000


def md_mismatches(md):
    """(mismatch bases, runs of adjacent mismatches) from an MD tag."""
    bases = runs = 0
    last_was_mm = False
    for num, dele, mm in MD_RE.findall(md):
        if mm:
            bases += 1
            if not last_was_mm:
                runs += 1
            last_was_mm = True
        elif num:
            if int(num) > 0:
                last_was_mm = False
        else:
            last_was_mm = False
    return bases, runs


def primary_row(r):
    cig = r.cigartuples or []
    mate = 1 if r.flag & 0x40 else 2 if r.flag & 0x80 else 0
    clip = sum(n for op, n in cig if op in (4, 5))
    qlen = r.query_length if r.query_length else len(r.query_sequence or '')
    hard = sum(n for op, n in cig if op == 5)
    row = dict(name=r.query_name.encode()[:80], mate=mate, read_len=qlen + hard, mapped=0, mapq=r.mapping_quality,
               chrom=-1, gstart=-1, gend=-1, nmap=len(cig))
    if r.is_unmapped:
        return row
    eq = x = m = ins = ins_ev = dl = del_ev = 0
    x_ev = 0
    for op, n in cig:
        if op == 7:
            eq += n
        elif op == 8:
            x += n
            x_ev += 1
        elif op == 0:
            m += n
        elif op == 1:
            ins += n
            ins_ev += 1
        elif op == 2:
            dl += n
            del_ev += 1
    if m:
        if r.has_tag('MD'):
            mm, mm_ev = md_mismatches(r.get_tag('MD'))
        elif r.has_tag('NM'):
            mm = max(0, r.get_tag('NM') - ins - dl)
            mm_ev = mm
        else:
            mm = mm_ev = 0
        match = m - mm
    else:
        # adjacent X runs separated by nothing are one CIGAR op already
        mm, mm_ev, match = x, x_ev, eq
    i = sc_l = 0
    while i < len(cig) and cig[i][0] in (4, 5):
        sc_l += cig[i][1]
        i += 1
    j, sc_r = len(cig) - 1, 0
    while j >= i and cig[j][0] in (4, 5):
        sc_r += cig[j][1]
        j -= 1
    row.update(mapped=1, match=match, mm=mm, mm_ev=mm_ev, ins=ins, ins_ev=ins_ev, dl=dl, del_ev=del_ev,
               sc_l=sc_l, sc_r=sc_r, chrom=r.reference_id, gstart=r.reference_start, gend=r.reference_end,
               nchrom=1, ref_visits=0)
    return row


def scan_region(job):
    bam, contig, start, end, mod, keep = job
    t0 = time.time()
    crc = zlib.crc32
    rows, supp = [], defaultdict(lambda: [0, 0, set()])
    counts = dict(records=0, primary=0, secondary=0, supplementary=0, sel_primary=0, sel_secondary=0,
                  sel_supplementary=0, unmapped_primary=0)
    hist = np.zeros(256, dtype=np.int64)
    with pysam.AlignmentFile(bam, 'rb', threads=1) as f:
        it = f.fetch(until_eof=True, contig='*') if contig == '*' else f.fetch(contig, start, end)
        for r in it:
            if contig != '*' and not (start <= r.reference_start < end):
                continue
            counts['records'] += 1
            flag = r.flag
            if flag & 0x100:
                counts['secondary'] += 1
            elif flag & 0x800:
                counts['supplementary'] += 1
            else:
                counts['primary'] += 1
                hist[r.mapping_quality] += 1
                if flag & 0x4:
                    counts['unmapped_primary'] += 1
            name = r.query_name
            if crc(name.encode()) % mod >= keep:
                continue
            if flag & 0x100:
                counts['sel_secondary'] += 1
            elif flag & 0x800:
                counts['sel_supplementary'] += 1
                mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
                s = supp[(name, mate)]
                s[0] += 1
                s[1] += r.query_alignment_length
                s[2].add(r.reference_id)
            else:
                counts['sel_primary'] += 1
                rows.append(primary_row(r))
    supp = {k: (v[0], v[1], sorted(v[2])) for k, v in supp.items()}
    out = np.zeros(len(rows), dtype=BROW)
    for k in ROW.names:
        if rows:
            out[k] = [r.get(k, 0) for r in rows]
    return dict(region=(contig, start, end), rows=out, supp=supp, counts=counts, hist=hist, seconds=time.time() - t0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bam', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--mod', type=int, default=1000)
    p.add_argument('--keep', type=int, default=1)
    p.add_argument('--procs', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
    p.add_argument('--chunk', type=int, default=CHUNK)
    p.add_argument('--test-region', help='test: contig:start-end only')
    args = p.parse_args()
    t0 = time.time()
    with pysam.AlignmentFile(args.bam, 'rb') as f:
        refs = list(zip(f.references, f.lengths))
    if args.test_region:
        c, se = args.test_region.split(':')
        s, e = (int(v) for v in se.split('-'))
        jobs = [(args.bam, c, s, e, args.mod, args.keep)]
    else:
        jobs = [(args.bam, c, s, min(L, s + args.chunk), args.mod, args.keep)
                for c, L in sorted(refs, key=lambda x: -x[1]) for s in range(0, L, args.chunk)]
        jobs.append((args.bam, '*', 0, 0, args.mod, args.keep))
    with ProcessPoolExecutor(max(1, args.procs), mp_context=multiprocessing.get_context('fork')) as pool:
        results = list(pool.map(scan_region, jobs))
    supp = defaultdict(lambda: [0, 0, set()])
    for res in results:
        for k, (n, bp, chroms) in res['supp'].items():
            s = supp[k]
            s[0] += n
            s[1] += bp
            s[2].update(chroms)
    out = np.concatenate([res['rows'] for res in results])
    keys = set()
    for i in range(len(out)):
        key = (out['name'][i].decode(), int(out['mate'][i]))
        keys.add(key)
        s = supp.get(key)
        if s:
            out['supp'][i] = s[0]
            out['supp_bp'][i] = s[1]
            other = len(set(s[2]) - {int(out['chrom'][i])})
            out['supp_other_chrom'][i] = other
            if out['chrom'][i] >= 0:
                out['nchrom'][i] = 1 + other
    counts = {k: sum(res['counts'][k] for res in results) for k in results[0]['counts']}
    np.savez_compressed(args.out + '.npz', rows=out, contigs=np.array(json.dumps([c for c, _ in refs])))
    summ = dict(bam=args.bam, mod=args.mod, keep=args.keep, regions=len(jobs), complete=not args.test_region,
                selected=int(len(out)), counts=counts,
                supplementary_without_primary=len(set(supp) - keys),
                mapq_hist_primary=np.sum([res['hist'] for res in results], axis=0).tolist(),
                region_seconds_max=round(max(res['seconds'] for res in results), 1),
                wall_seconds=round(time.time() - t0, 1))
    json.dump(summ, open(args.out + '.summary.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != 'mapq_hist_primary'}))


if __name__ == '__main__':
    main()
