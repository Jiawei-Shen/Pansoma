"""s4b: ALT / REF k-mer read counts at the s4a loci in one sample's GRCh38 reads.

usage: s4_count.py TUMOR LABEL BAM[,BAM...] SHARD NSHARDS [SELECT [TAG]]

Loci of DATA/s4/<TUMOR>/loci.tsv (SELECT all (default) / cand / germ: only loci with that purpose flag), split into
NSHARDS contiguous chunks (sorted by position); this task does chunk SHARD. With TAG, loci already in an earlier count
table of the label (counts/<LABEL>.*.tsv not named <LABEL>.<TAG>.*) are left out first (locus ids are content hashes),
and the output is counts/<LABEL>.<TAG>.<SHARD>.tsv.
Window of a locus = [a_start - W, a_end + W), W = longest k-mer + 20 bp.
Records: primary and supplementary (no secondary, unmapped, QC-fail or duplicate flag) overlapping the window.
For each record only the part of the read aligned to the window is searched (CIGAR walk; insertions inside the window
included; a soft clip next to an aligned end inside the window included up to the window edge), so copies of a k-mer
elsewhere in a long read do not count.
A fragment is counted once per locus: paired reads (Illumina) by (pair start, pair end, |TLEN|, orientation), which also
merges PCR / optical duplicates that the BAM does not flag and overlapping mates; other reads by name.
Columns (fragments): frags (overlapping the window), span (a record aligned over [a_start - k, a_end + k), k = longest
k-mer), alt / ref (ALT / REF
k-mer on either strand in the searched part), span_mq10 / alt_mq10 / ref_mq10 (the same with MAPQ >= 10).
Writes DATA/s4/<TUMOR>/counts/<LABEL>.<SHARD>.tsv.
"""
import bisect
import csv
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, rc  # noqa: E402

SKIP = 0x4 | 0x100 | 0x200 | 0x400
EXTRA = 20
GAP = 1000
MAX_BLOCK = 200_000
COLS = ('frags', 'span', 'alt', 'ref', 'span_mq10', 'alt_mq10', 'ref_mq10')


def blocks(loci):
    cur, end = [], 0
    for lc in loci:
        if cur and lc['chrom'] == cur[0]['chrom'] and lc['ws'] <= end + GAP and lc['we'] - cur[0]['ws'] <= MAX_BLOCK:
            cur.append(lc)
            end = max(end, lc['we'])
        else:
            if cur:
                yield cur
            cur, end = [lc], lc['we']
    if cur:
        yield cur


def cigar_index(cig, rstart):
    """op starts (ref, query, is_match) of the reference-consuming ops plus an end sentinel; soft-clip lengths."""
    rs, qs, mt = [], [], []
    r, q, lclip, rclip, seen, q_last = rstart, 0, 0, 0, False, 0
    for op, n in cig:
        if op == 4:
            if seen:
                rclip += n
            else:
                lclip += n
            q += n
        elif op == 1:
            q += n
        elif op in (0, 7, 8):
            seen = True
            rs.append(r)
            qs.append(q)
            mt.append(True)
            r += n
            q += n
            q_last = q
        elif op in (2, 3):
            seen = True
            rs.append(r)
            qs.append(q)
            mt.append(False)
            r += n
            q_last = q
    rs.append(r)
    qs.append(q_last)
    mt.append(False)
    return rs, qs, mt, lclip, rclip


def qat(ix, x):
    """query position aligned to reference position x (inside a deletion: the next aligned query base)."""
    rs, qs, mt = ix[0], ix[1], ix[2]
    i = max(0, bisect.bisect_right(rs, x) - 1)
    return qs[i] + (x - rs[i]) if mt[i] else qs[i]


def frag_key(rd):
    if rd.is_paired and not rd.mate_is_unmapped and not rd.is_supplementary and rd.reference_id == rd.next_reference_id:
        a, b = sorted((rd.reference_start, rd.next_reference_start))
        return ('p', a, b, abs(rd.template_length), rd.is_read1 != rd.is_reverse)
    return rd.query_name


def main():
    tumor, label, bams, shard, nshards = sys.argv[1], sys.argv[2], sys.argv[3].split(','), int(sys.argv[4]), int(sys.argv[5])
    select = sys.argv[6] if len(sys.argv) > 6 else 'all'
    tag = sys.argv[7] if len(sys.argv) > 7 else None
    if select not in ('all', 'cand', 'germ'):
        raise SystemExit(f'SELECT {select}')
    out_dir = f'{DATA}/s4/{tumor}/counts'
    done = set()
    if tag:
        import glob
        for p in glob.glob(f'{out_dir}/{label}.*.tsv'):
            if os.path.basename(p).startswith(f'{label}.{tag}.'):
                continue
            with open(p) as f:
                next(f)
                done.update(line.split('\t', 1)[0] for line in f)
    loci = []
    for r in csv.DictReader(open(f'{DATA}/s4/{tumor}/loci.tsv'), delimiter='\t'):
        if select != 'all' and r.get(select) != '1':
            continue
        if r['locus_id'] in done:
            continue
        w = max(len(r['alt_k']), len(r['ref_k'])) + EXTRA
        r['ws'], r['we'] = max(0, int(r['a_start']) - w), int(r['a_end']) + w
        r['ss'], r['se'] = int(r['a_start']) - (w - EXTRA), int(r['a_end']) + (w - EXTRA)
        loci.append(r)
    lo, hi = len(loci) * shard // nshards, len(loci) * (shard + 1) // nshards
    loci = sorted(loci[lo:hi], key=lambda x: (x['chrom'], x['ws']))
    files = [pysam.AlignmentFile(b, 'rb') for b in bams]
    os.makedirs(out_dir, exist_ok=True)
    out = f'{out_dir}/{label}.{tag}.{shard}.tsv' if tag else f'{out_dir}/{label}.{shard}.tsv'
    counts = {}
    n_rec = 0
    for blk in blocks(loci):
        blk = sorted(blk, key=lambda x: x['ws'])
        sets = {x['locus_id']: tuple(set() for _ in COLS) for x in blk}
        starts = [x['ws'] for x in blk]
        span = max(x['we'] - x['ws'] for x in blk)
        pats = [(x['locus_id'], x['ws'], x['we'], x['ss'], x['se'], x['alt_k'], rc(x['alt_k']), x['ref_k'], rc(x['ref_k']))
                for x in blk]
        c0, s0, e0 = blk[0]['chrom'], blk[0]['ws'], max(x['we'] for x in blk)
        for f in files:
            if c0 not in f.references:
                continue
            for rd in f.fetch(c0, s0, e0):
                if rd.flag & SKIP:
                    continue
                seq = rd.query_sequence
                if not seq or not rd.cigartuples:
                    continue
                rs, re_ = rd.reference_start, rd.reference_end
                j = bisect.bisect_right(starts, re_)
                ix = None
                mq = rd.mapping_quality >= 10
                fk = None
                for i in range(j - 1, -1, -1):
                    lid, ws, we, ss, se, ak, akr, rk, rkr = pats[i]
                    if ws < rs - span:
                        break
                    if we <= rs or ws >= re_:
                        continue
                    if ix is None:
                        ix = cigar_index(rd.cigartuples, rs)
                        fk = frag_key(rd)
                        n_rec += 1
                    q0 = qat(ix, max(ws, rs))
                    if rs > ws:
                        q0 -= min(ix[3], rs - ws)
                    q1 = qat(ix, min(we, re_))
                    if re_ < we:
                        q1 += min(ix[4], we - re_)
                    part = seq[max(0, q0):q1]
                    st = sets[lid]
                    st[0].add(fk)
                    sp = rs <= ss and re_ >= se
                    if sp:
                        st[1].add(fk)
                        if mq:
                            st[4].add(fk)
                    if ak in part or akr in part:
                        st[2].add(fk)
                        if mq:
                            st[5].add(fk)
                    if rk in part or rkr in part:
                        st[3].add(fk)
                        if mq:
                            st[6].add(fk)
        for lid, st in sets.items():
            counts[lid] = [len(x) for x in st]
    with open(out + '.tmp', 'w') as w:
        w.write('locus_id\t' + '\t'.join(COLS) + '\n')
        for lc in loci:
            w.write(lc['locus_id'] + '\t' + '\t'.join(map(str, counts[lc['locus_id']])) + '\n')
    os.replace(out + '.tmp', out)
    print(tumor, label, select, tag, shard, 'loci', len(loci), 'records', n_rec, 'skipped as counted', len(done), 'done')


if __name__ == '__main__':
    main()
