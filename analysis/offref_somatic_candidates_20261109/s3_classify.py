"""s3: assembly class of every s1 row of one tumor, from the s2 alignments to the normal and tumor haplotypes and GRCh38.

usage: s3_classify.py TUMOR

Strict hit (E9d / E9h h50 strict rule): alignment covering >= 90 % of the query, NM <= 1, no inserted or deleted base,
no mismatch within one base of the site interval (a deletion: the two bases at the junction).
Event window (common.event_window): site +- 1, widened to the tandem runs touching it and to every equivalent placement
of an INDEL, so the window does not depend on where minimap2 puts a gap (left-aligned on the target, i.e. on the right
in query coordinates for '-' hits).
Differences of an alignment (diffs): mismatches, inserted query bases, deletion junctions and unaligned query ends; an
unaligned end that overlaps the window is split at the window edge.

Per row (ALT query alt_q, REF-swap query ref_q):
  asm_class   no_query     s1 status not ok, or alt_q shorter than 30 bp
              short_flank  < 10 read-consensus bases on a side of the site, or < 5 beyond the event window
              normal       a strict hit of alt_q in a normal haplotype (the A1 haplotype is the patient's germline)
              tumor_only   a strict hit in a tumor haplotype and none in the normal ones
              neither      no strict hit
  normal_side (tumor_only / neither; what the normal assembly has instead)
              ref_exact       ref_q has a strict normal hit
              site_only       the best normal alignment of alt_q (fewest differences outside the event window, then NM,
                              then AS) differs from it only inside the window
              site_plus_other differences inside the window and outside it (<= 10 outside)
              no_site_diff    no difference inside the window (the normal carries the ALT there)
              no_alignment    no normal alignment covering >= 90 % of alt_q, or > 10 differences outside the window
  t_loci / t_max   strict-hit loci (contig, 50-bp bin; adjacent bins merged) of alt_q in the tumor: sum over the two
                   haplotypes / largest per haplotype; n_loci / n_max the same for ref_q in the normal; n_alt_loci for
                   alt_q in the normal; multi_copy = t_max > 1 or n_max > 1
  tumor_pos / normal_pos  target:contig:pos (0-based, leftmost base of the site on either strand; a zero-length site:
                   the base right of the junction) of a strict hit (the haplotype with the fewest loci) of alt_q in the
                   tumor / of ref_q in the normal (ref_exact), else of the site on the best normal alignment of alt_q
  normal_aln       target:contig:start-end:strand of the best normal alignment of alt_q; normal_out its differences
                   outside the window; normal_best_nm its NM
  norm_q / norm_win  site_only / site_plus_other: the normal sequence of that alignment in alt_q orientation (aligned
                   part only) and the event window mapped onto it (the patient's own allele at the site)
  anchor_mapped    chrom:start-end of the site of ref_q (else alt_q) on its primary GRCh38 alignment with MAPQ >= 20
                   (g38_mapq); alt_in_grch38 = strict GRCh38 loci of alt_q (the ALT context exists in GRCh38)
Writes DATA/s3/<TUMOR>/rows.tsv.gz (s1 columns + the above), candidates.tsv.gz (asm_class tumor_only, or neither with
normal_side ref_exact / site_only) and counts.tsv.
"""
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, MIN_LEN, TARGETS, qid, rc, ref_site, s1_files, site_of, usable, windows  # noqa: E402

CS = re.compile(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})')
MAX_OUT = 10
MIN_MAPQ = 20
BIN = 50


def site_exact(cs, start, lo, hi):
    """E9d e9d_asm.site_exact: cs walk from query position `start`; True when [lo, hi) lies inside the alignment with no
    mismatch / insertion inside and no deletion strictly inside."""
    p = start
    if lo < start:
        return False
    for op in CS.findall(cs):
        t = op[0]
        if t == ':':
            p += int(op[1:])
        elif t == '*':
            if lo <= p < hi:
                return False
            p += 1
        elif t == '+':
            n = len(op) - 1
            if p < hi and p + n > lo:
                return False
            p += n
        elif t == '-':
            if lo < p < hi:
                return False
        elif t == '~':
            return False
        if p >= hi and t == ':':
            return True
    return p >= hi


def diffs(cs, qs, qe, qlen, strand):
    """differences of an alignment as (lo, hi) intervals in original query coordinates; a deletion from the query is
    the zero-length interval at its junction; unaligned query ends count as differences."""
    out = []
    p = qs if strand == '+' else qlen - qe          # walk in alignment (query-forward-on-target) orientation
    start = p
    for op in CS.findall(cs):
        t = op[0]
        if t == ':':
            p += int(op[1:])
        elif t == '*':
            out.append((p, p + 1))
            p += 1
        elif t == '+':
            n = len(op) - 1
            out.append((p, p + n))
            p += n
        elif t in '-~':
            out.append((p, p))
    end = p
    if start > 0:
        out.append((0, start))
    if end < qlen:
        out.append((end, qlen))
    if strand == '-':
        out = [(qlen - hi, qlen - lo) for lo, hi in out]
    return out


def split(d, wl, wh, qlen):
    """(inside, outside) the window [wl, wh): an interval with bases is inside when it overlaps the window; a zero-length
    junction when wl <= x <= wh; an unaligned query end is split at the window edges."""
    ins, out = [], []
    for lo, hi in d:
        if lo == hi:
            (ins if wl <= lo <= wh else out).append((lo, hi))
            continue
        if not (lo < wh and hi > wl):
            out.append((lo, hi))
            continue
        if lo == 0 or hi == qlen:
            if lo < wl:
                out.append((lo, wl))
                lo = wl
            if hi > wh:
                out.append((wh, hi))
                hi = wh
        ins.append((lo, hi))
    return ins, out


def target_map(q, cs, qs, qe, qlen, strand):
    """target bases of the alignment (target-forward) and tmap[j] = target offset of the junction before base j of the
    aligned query (in alignment orientation); a deletion's junction maps after the deleted bases."""
    qq = q[qs:qe] if strand == '+' else rc(q)[qlen - qe:qlen - qs]
    t, tmap, i, x = [], [0] * (len(qq) + 1), 0, 0
    for op in CS.findall(cs):
        c = op[0]
        if c == ':':
            n = int(op[1:])
            t.append(qq[i:i + n])
            for j in range(n):
                tmap[i + j] = x + j
            i += n
            x += n
        elif c == '*':
            t.append(op[1].upper())
            tmap[i] = x
            i += 1
            x += 1
        elif c == '+':
            n = len(op) - 1
            for j in range(n):
                tmap[i + j] = x
            i += n
        elif c == '-':
            t.append(op[1:].upper())
            x += len(op) - 1
    tmap[i] = x
    return ''.join(t), tmap


def tjunc(tmap, qs, qe, strand, p):
    """target offset (from the alignment start, target-forward) of query junction p."""
    p = min(max(p, qs), qe)
    return tmap[p - qs] if strand == '+' else tmap[qe - p]


def site_on_target(q, aln, a, b):
    """(start, end) target-forward interval of query site [a, b) on an alignment; a zero-length site: (x, x + 1) with x
    the base right of the junction."""
    ts, qs, qe, qlen, strand, cs = aln
    _, tmap = target_map(q, cs, qs, qe, qlen, strand)
    x1, x2 = sorted((tjunc(tmap, qs, qe, strand, a), tjunc(tmap, qs, qe, strand, b)))
    return ts + x1, ts + max(x2, x1 + 1)


def norm_context(q, aln, wl, wh):
    """normal sequence of an alignment of q in q's orientation (aligned part) and window [wl, wh) mapped onto it."""
    ts, qs, qe, qlen, strand, cs = aln
    t, tmap = target_map(q, cs, qs, qe, qlen, strand)
    x1, x2 = tjunc(tmap, qs, qe, strand, wl), tjunc(tmap, qs, qe, strand, wh)
    if strand == '+':
        return t, x1, x2
    return rc(t), len(t) - x1, len(t) - x2


def nloci(d):
    """loci of a {(contig, bin): pos} map, adjacent bins merged."""
    return sum(1 for (c, k) in d if (c, k - 1) not in d)


def load_targets(tumor):
    t = []
    for line in open(TARGETS):
        f = line.rstrip('\n').split('\t')
        if f[0] == tumor:
            t.append((f[1], f[2]))
    if sorted(r for _, r in t) != ['grch38', 'normal', 'normal', 'tumor', 'tumor']:
        raise SystemExit(f'{tumor}: targets {t}')
    return t


def parse(tumor, qinfo, targets):
    strict = defaultdict(lambda: defaultdict(dict))  # qid -> target -> {(contig, bin): site pos}
    best = {}                                         # qid -> best normal alignment
    g38 = {}                                          # qid -> primary GRCh38 alignment
    fa_time = os.path.getmtime(f'{DATA}/s2/{tumor}/queries.fa')
    for target, role in targets:
        paf = f'{DATA}/s2/{tumor}/{target}.paf'
        if not os.path.exists(paf):
            raise SystemExit(f'missing {paf}')
        if os.path.getmtime(paf) < fa_time:
            raise SystemExit(f'stale {paf}: older than queries.fa')
        n_line = n_unknown = 0
        for line in open(paf):
            n_line += 1
            f = line.rstrip('\n').split('\t')
            name = f[0]
            info = qinfo.get(name)
            if info is None:
                n_unknown += 1
                continue
            a, b, wl, wh = info
            qlen, qs, qe, strand = int(f[1]), int(f[2]), int(f[3]), f[4]
            ts, te = int(f[7]), int(f[8])
            tags = dict((x[:2], x[5:]) for x in f[12:])
            nm = int(tags.get('NM', 99))
            cs = tags.get('cs', '')
            if role == 'normal':
                _, out = split(diffs(cs, qs, qe, qlen, strand), wl, wh, qlen)
                score = (len(out), nm, -int(tags.get('AS', 0)))
                if name not in best or score < best[name][:3]:
                    best[name] = (*score, target, f[5], ts, te, qs, qe, qlen, strand, cs)
            elif role == 'grch38' and tags.get('tp') == 'P':
                mapq = int(f[11])
                if name not in g38 or mapq > g38[name][0]:
                    g38[name] = (mapq, f[5], ts, qs, qe, qlen, strand, cs)
            if nm > 1 or '+' in cs or '-' in cs or '~' in cs:
                continue
            lo, hi = max(0, a - 1), min(qlen, max(b, a) + 1)
            if strand == '-':
                lo2, hi2, start = qlen - hi, qlen - lo, qlen - qe
            else:
                lo2, hi2, start = lo, hi, qs
            if site_exact(cs, start, lo2, hi2):
                # no indels in a strict hit: query position -> target position is linear
                pos = ts + (a - qs) if strand == '+' else te - (max(a, b) - qs)
                strict[name][target].setdefault((f[5], pos // BIN), pos)
        if n_unknown:
            raise SystemExit(f'{paf}: {n_unknown} of {n_line} lines name queries not in queries.tsv (stale s2)')
        print(f'{target}: {n_line:,} PAF lines', flush=True)
    return strict, best, g38


def main():
    tumor = sys.argv[1]
    targets = load_targets(tumor)
    roles = dict(targets)
    qtsv = f'{DATA}/s2/{tumor}/queries.tsv'
    stale = [p for p in s1_files(tumor) if os.path.getmtime(p) > os.path.getmtime(qtsv)]
    if stale:
        raise SystemExit(f'{tumor}: s1 tables newer than s2 queries.tsv: {stale[:3]}')
    qinfo = {}
    for r in csv.DictReader(open(qtsv), delimiter='\t'):
        a, b = map(int, r['site'].split('-'))
        wl, wh = map(int, r['win'].split('-'))
        qinfo[r['qid']] = (a, b, wl, wh)
    strict, best, g38 = parse(tumor, qinfo, targets)
    print(f'{tumor}: {len(qinfo):,} queries, {len(strict):,} with a strict hit, {len(best):,} with a normal alignment, '
          f'{len(g38):,} with a primary GRCh38 alignment', flush=True)
    out = f'{DATA}/s3/{tumor}'
    os.makedirs(out, exist_ok=True)
    cols = None
    counts = defaultdict(int)
    with gzip.open(f'{out}/rows.tsv.gz.tmp', 'wt') as wr, gzip.open(f'{out}/candidates.tsv.gz.tmp', 'wt') as wc:
        for p in s1_files(tumor):
            for r in csv.DictReader(open(p), delimiter='\t'):
                new = dict(asm_class='', normal_side='', t_loci=0, t_max=0, n_loci=0, n_max=0, n_alt_loci=0,
                           multi_copy='', tumor_pos='', normal_pos='', normal_aln='', normal_out='', normal_best_nm='',
                           norm_q='', norm_win='', anchor_mapped='', g38_mapq='', alt_in_grch38='', q_alt='', q_ref='')
                u = usable(r)
                if u != 'ok':
                    new['asm_class'] = u
                else:
                    alt, ref = r['alt_q'], r['ref_q']
                    a0, b0 = site_of(r)
                    ra, rb = map(int, ref_site(r).split('-'))
                    qa = qid(alt, r['site'])
                    qr = qid(ref, f'{ra}-{rb}') if len(ref) >= MIN_LEN else ''
                    if qa not in qinfo or (qr and qr not in qinfo):
                        raise SystemExit(f'{p}: query {qa}/{qr} not in s2 queries.tsv; rerun s2')
                    new['q_alt'], new['q_ref'] = qa, qr
                    hit = strict.get(qa, {})
                    alt_n = [t for t in hit if roles[t] == 'normal']
                    alt_t = [t for t in hit if roles[t] == 'tumor']
                    new['alt_in_grch38'] = sum(nloci(hit[t]) for t in hit if roles[t] == 'grch38')
                    new['t_loci'] = sum(nloci(hit[t]) for t in alt_t)
                    new['t_max'] = max([nloci(hit[t]) for t in alt_t], default=0)
                    new['n_alt_loci'] = sum(nloci(hit[t]) for t in alt_n)
                    if alt_t:
                        t0 = min(alt_t, key=lambda t: nloci(hit[t]))
                        (c, _), pos = next(iter(hit[t0].items()))
                        new['tumor_pos'] = f'{t0}:{c}:{pos}'
                    rhit = strict.get(qr, {}) if qr else {}
                    ref_n = [t for t in rhit if roles[t] == 'normal']
                    new['n_loci'] = sum(nloci(rhit[t]) for t in ref_n)
                    new['n_max'] = max([nloci(rhit[t]) for t in ref_n], default=0)
                    new['multi_copy'] = int(new['t_max'] > 1 or new['n_max'] > 1)
                    new['asm_class'] = 'normal' if alt_n else 'tumor_only' if alt_t else 'neither'
                    b = best.get(qa)
                    if b:
                        n_out, nm, _, tgt, contig, ts, te, qs, qe, qlen, strand, cs = b
                        new['normal_out'], new['normal_best_nm'] = n_out, nm
                        new['normal_aln'] = f'{tgt}:{contig}:{ts}-{te}:{strand}'
                    if new['asm_class'] != 'normal':
                        if ref_n:
                            new['normal_side'] = 'ref_exact'
                            t0 = min(ref_n, key=lambda t: nloci(rhit[t]))
                            (c, _), pos = next(iter(rhit[t0].items()))
                            new['normal_pos'] = f'{t0}:{c}:{pos}'
                        elif b is None or b[0] > MAX_OUT:
                            new['normal_side'] = 'no_alignment'
                        else:
                            (wl, wh), _ = windows(r)
                            ins, outside = split(diffs(cs, qs, qe, qlen, strand), wl, wh, qlen)
                            side = 'no_site_diff' if not ins else 'site_plus_other' if outside else 'site_only'
                            new['normal_side'] = side
                            aln = (ts, qs, qe, qlen, strand, cs)
                            s, _ = site_on_target(alt, aln, a0, b0)
                            new['normal_pos'] = f'{tgt}:{contig}:{s}'
                            if side != 'no_site_diff':
                                nq, nl, nh = norm_context(alt, aln, wl, wh)
                                new['norm_q'], new['norm_win'] = nq, f'{nl}-{nh}'
                    for q, (sa, sb), name in ((ref, (ra, rb), qr), (alt, (a0, b0), qa)):
                        g = g38.get(name) if name else None
                        if g and g[0] >= MIN_MAPQ and g[3] <= sa and sb <= g[4]:
                            mapq, contig, ts, qs, qe, qlen, strand, cs = g
                            s, e = site_on_target(q, (ts, qs, qe, qlen, strand, cs), sa, sb)
                            new['anchor_mapped'], new['g38_mapq'] = f'{contig}:{s}-{e}', mapq
                            break
                if cols is None:
                    cols = list(r) + list(new)
                    wr.write('\t'.join(cols) + '\n')
                    wc.write('\t'.join(cols) + '\n')
                line = '\t'.join(str({**r, **new}[c]) for c in cols) + '\n'
                wr.write(line)
                counts[(r['set'], r['kind'], r['group'], new['asm_class'], new['normal_side'], new['multi_copy'])] += 1
                if new['asm_class'] == 'tumor_only' or (new['asm_class'] == 'neither' and
                                                        new['normal_side'] in ('ref_exact', 'site_only')):
                    wc.write(line)
    os.replace(f'{out}/rows.tsv.gz.tmp', f'{out}/rows.tsv.gz')
    os.replace(f'{out}/candidates.tsv.gz.tmp', f'{out}/candidates.tsv.gz')
    with open(f'{out}/counts.tsv', 'w') as w:
        w.write('set\tkind\tgroup\tasm_class\tnormal_side\tmulti_copy\tn\n')
        for k, v in sorted(counts.items()):
            w.write('\t'.join(map(str, k)) + f'\t{v}\n')
    print(f'{tumor}: done -> {out}')


if __name__ == '__main__':
    main()
