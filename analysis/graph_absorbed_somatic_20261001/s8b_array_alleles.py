"""Step 8b: whole-tandem-array alleles of the HG008-N haplotypes and the HG008-T contigs at every locus (sequence level).

Why: s6 (minimap2 sr, short context) judges only the s1 event window (period 1-6 stretches touching the event +-5 bp). In
interrupted / compound arrays and VNTRs (period > 6) a length difference can slide out of that window, and a split
representation (45 + 4 bp for a 49-bp deletion) passes the one-or-two-difference rule (audit hg008n_hap_membership).
Here every haplotype is compared over the whole array between common unique anchors, plus a second source (dipcall).
Inputs: loci.tsv (s1); asm_seq_status.tsv (s6 placement: graph_win_t0/t1 on the contig, strand of ref_span / alt_span);
  GRCh38 FASTA; HG008-N v6.2 bothhaps FASTA; HG008-T v3.2 FASTA (s6 copy); HG008-N dipcall dip.vcf.gz + dip.bed (the
  labels-manifest germline truth).
Method
 1 array [arr_lo, arr_hi) (GRCh38, 0-based half-open), grown by exact periodic stretches (seq[i] == seq[i + p] runs)
   touching it with 1-bp tolerance: first round from the s1 event core (SNV base / INDEL bases after the padding base)
   with the s1 length rule max(2p, p + 3) for period 1-6 and >= 2 copies for period 7-60; then, repeated until stable,
   period 1-60 stretches of >= max(2p, 10) bp (adjacent / compound arrays, VNTRs); united with the s1 event window;
   growth stops at 3 kb.
 2 common anchors: GRCh38 24-mers ending at a <= arr_lo / starting at b >= arr_hi, nearest first, up to 400 bp out, not
   low-complexity (>= 3 bases; no period <= 6 matching >= 80 % of the k-mer), once in GRCh38 +-2.5 kb and exactly once in
   each placed contig segment (s6 projected window +-2.5 kb + array length, oriented to GRCh38 by the s6 strand; the right
   anchor after the left one). The nearest anchor every placed contig has wins, else the nearest that both normal haps
   have, else the one most contigs have; a contig without the chosen anchors is NA (anchor).
 3 S = a contig's sequence between the anchors; R = GRCh38[a, b); A = R with the truth allele applied.
   call: ALT (S == A) / ALT_plus_flank (S == A with the hap's own dipcall records that lie outside the array, in [a, arr_lo)
   or [arr_hi, b), applied) / REF (S == R) / other; len_change = len(S) - len(R).
 4 dipcall per normal hap: R with that hap's PASS ('.') records applied by GT (hap1|hap2); NA if [a, b) is not inside one
   dip.bed interval, a non-PASS record overlaps [a, b), a record crosses a or b, a GT is missing or applied records overlap.
 5 per normal hap: seq = the assembly S, else the dipcall sequence; asm_eq_dip where both exist.
 6 tumor contigs: each S against the two normal seqs: hap1 / hap2 / hom (both haps equal) / novel. tumor_verdict:
   novel_ALT (a novel copy == A) / novel_other / no_change_hom (haps equal, every copy equals them) / no_change_both
   (copies equal hap1 and hap2) / all_eq_hap1 / all_eq_hap2 / NA (no tumor S or a normal seq missing).
 7 truth INFO event (asm_event 'chrN_hapK:pos0-REF-ALT', normal-assembly frame, 0-based): when it lies inside the event
   hap's anchored stretch, it is applied to the event hap's S: T_info = the tumor event-hap allele by GIAB's own
   normal-frame call; T_info_class eq_ALT / eq_other_hap / eq_ALT=other_hap / novel / event_outside_stretch / NA.
Output: array_alleles.tsv (loci.tsv order; sequences included, '' when NA).
Assumptions: the s6 placement picks the right copy (a wrong copy usually lacks the anchors -> NA, or gives another
sequence than dipcall -> asm_eq_dip False); exact periodic runs miss imperfect VNTR copies, so an array can end early (both
sides are still compared over the same stretch); dipcall GT order = hap1|hap2 (audit: 4,825 het SNVs, 0 swapped);
'*' alleles are left to the overlapping deletion record.
Run: python s8b_array_alleles.py [truth_id,...] (with ids: a test printout only); login node, 10 min IO-bound FASTA
fetches, 72 MB (2026-10-01 run: $T/s8b_array_alleles.login.log); the dipcall VCF is read in one pass (no index).
"""
import bisect, collections, csv, re, sys
import numpy as np
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
GR = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
ASM = {'normal': '/scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_bothhaps_polished6.2.fasta.gz',
       'tumor': '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/assembly/HG008T_v3.2.fasta'}
DIP = '/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip'
K, CTX, SHIFT, MAXARR = 24, 2500, 400, 3000
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
rc = lambda s: s.translate(RC)[::-1]


def stretches(g, g0):
    """[(start, end, p)] exact periodic stretches of g (absolute coordinates), period 1-60, >= 2 copies and >= p + 3."""
    a, out = np.frombuffer(g.encode(), dtype=np.uint8), []
    for p in range(1, 61):
        d = np.diff(np.concatenate(([0], (a[:-p] == a[p:]).astype(np.int8), [0])))
        s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        m = (e - s) >= max(p, 3)
        out += [(g0 + i, g0 + j + p, p) for i, j in zip(s[m], e[m])]
    return out


def grow(st, lo, hi):
    """round 0 from the event core [lo, hi) (the s1 rule, plus period 7-60 with 2 copies), then adjacent arrays."""
    ok1 = lambda x, y, p: (y - x >= max(2 * p, p + 3)) if p <= 6 else y - x >= 2 * p
    ok2 = lambda x, y, p: y - x >= max(2 * p, 10)
    for rnd in range(100):
        f = ok1 if rnd == 0 else ok2
        n = [(x, y) for x, y, p in st if x <= hi + 1 and y >= lo - 1 and f(x, y, p)]
        nlo, nhi = min([lo] + [x for x, _ in n]), max([hi] + [y for _, y in n])
        if nhi - nlo > MAXARR:
            return lo, hi, True
        if (nlo, nhi) == (lo, hi) and rnd > 0:
            return lo, hi, False
        lo, hi = nlo, nhi
    return lo, hi, True


def lowc(k):
    return len(set(k)) < 3 or any(sum(k[i] == k[i + p] for i in range(K - p)) >= 0.8 * (K - p) for p in range(1, 7))


def occ(s, k, start=0):
    i = s.find(k, start)
    return [] if i < 0 else [i] if s.find(k, i + 1) < 0 else [i, s.find(k, i + 1)]


def anchors(g, g0, lo, hi, segs, normals):
    """(a, b, {contig: (left offset, right offset)}) - see docstring 2."""
    def pick(cands):
        best, first_all, first_norm = None, None, None
        for pos, ok in cands:
            n = len(ok)
            if n == len(segs):
                return pos, ok
            if first_norm is None and normals <= set(ok):
                first_norm = (pos, ok)
            if best is None or n > len(best[1]):
                best = (pos, ok)
        return first_norm or best or (None, {})
    def left():
        for e in range(lo, lo - SHIFT, -1):
            k = g[e - K - g0:e - g0]
            if len(k) < K or 'N' in k or lowc(k) or len(occ(g, k)) != 1:
                continue
            yield e, {c: o[0] for c, s in segs.items() for o in [occ(s, k)] if len(o) == 1}
    a, lo_ok = pick(left())
    if a is None:
        return None, None, {}
    def right():
        for b in range(hi, hi + SHIFT):
            k = g[b - g0:b - g0 + K]
            if len(k) < K or 'N' in k or lowc(k) or len(occ(g, k)) != 1:
                continue
            ok = {}
            for c, i in lo_ok.items():
                o = occ(segs[c], k, i + K)
                if len(o) == 1 and segs[c].find(k) == o[0]:
                    ok[c] = (i + K, o[0])
            yield b, ok
    b, ok = pick(right())
    return (a, b, ok) if b is not None else (a, None, {})


def edit(R, a, eds):
    out, cur = [], a
    for p0, r, x in sorted(eds):
        if p0 < cur or R[p0 - a:p0 - a + len(r)] != r:
            return None
        out += [R[cur - a:p0 - a], x]; cur = p0 + len(r)
    return ''.join(out) + R[cur - a:]


def dip_hap(R, a, b, recs, h, inbed):
    """(sequence, NA reason, the applied records) for hap index h (0 = hap1)."""
    if not inbed:
        return None, 'outside_dip_bed', []
    eds = []
    for p0, r, alts, gt, ok in recs:
        if not ok:
            return None, 'dip_nonpass_record', []
        g = gt[h] if gt and len(gt) > h else None
        if g is None:
            return None, 'dip_gt_missing', []
        if g == 0 or alts[g - 1] == '*':
            continue
        if p0 < a or p0 + len(r) > b:
            return None, 'dip_record_crosses_anchor', []
        eds.append((p0, r, alts[g - 1]))
    s = edit(R, a, eds)
    return (s, '', eds) if s is not None else (None, 'dip_overlapping_records', [])


def main():
    gfa = pysam.FastaFile(GR); fas = {k: pysam.FastaFile(v) for k, v in ASM.items()}
    clen = {k: dict(zip(f.references, f.lengths)) for k, f in fas.items()}
    loci = list(csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t'))
    test = sys.argv[1].split(',') if sys.argv[1:] else None             # test: comma-separated truth_ids -> stdout only
    if test:
        loci = [L for L in loci if L['truth_id'] in test]
    AS = collections.defaultdict(list)
    for r in csv.DictReader(open(f'{D}/asm_seq_status.tsv'), delimiter='\t'):
        AS[r['truth_id']].append(r)
    bed = collections.defaultdict(list)
    for l in open(f'{DIP}.bed'):
        c, s, e = l.split()[:3]; bed[c].append((int(s), int(e)))
    for c in bed:
        bed[c].sort()
    inbed = lambda c, s, e: (lambda i: i >= 0 and bed[c][i][0] <= s and e <= bed[c][i][1])(bisect.bisect_right(bed[c], (s, 10 ** 12)) - 1)
    out = []
    for n, L in enumerate(loci):
        t, ch, p0 = L['truth_id'], L['chrom'], int(L['vcf_pos']) - 1
        vr, va = L['vcf_ref'].upper(), L['vcf_alt'].upper()
        g0 = max(0, p0 - CTX - MAXARR); g = gfa.fetch(ch, g0, p0 + len(vr) + CTX + MAXARR).upper()
        core = (p0, p0 + 1) if len(vr) == 1 and len(va) == 1 else (p0 + 1, max(p0 + len(vr), p0 + 1))   # the s1 event core
        lo, hi, capped = grow(stretches(g, g0), *core)
        lo, hi = min(lo, int(L['ev_lo0'])), max(hi, int(L['ev_hi0']))
        segs, meta = {}, {}
        for r in AS[t]:
            sp = r['ref_span'] or r['alt_span']
            if not sp or not r['graph_win_t0']:
                continue
            c, asm = r['contig'], r['assembly']; strand = sp[-1]
            t0, t1 = sorted((int(r['graph_win_t0']), int(r['graph_win_t1'])))
            s0, s1 = max(0, t0 - CTX - (hi - lo)), min(clen[asm][c], t1 + CTX + (hi - lo))
            if s1 <= s0:
                continue
            s = fas[asm].fetch(c, s0, s1).upper()
            segs[(asm, c)] = s if strand == '+' else rc(s); meta[(asm, c)] = (s0, s1, strand)
        normals = {k for k in segs if k[0] == 'normal'}
        a, b, ok = anchors(g, g0, lo, hi, segs, normals)
        o = dict(truth_id=t, arr_lo0=lo, arr_hi0=hi, arr_len=hi - lo, arr_capped=capped, anchor_a0=a if a is not None else '',
                 anchor_b0=b if b is not None else '', n_contigs_placed=len(segs), n_contigs_anchored=len(ok))
        R = A = None
        if a is not None and b is not None:
            R = g[a - g0:b - g0]; A = edit(R, a, [(p0, vr, va)])
        o['_R'], o['_A'], o['_inbed'] = R, A, R is not None and inbed(ch, a, b)
        o['_S'] = {}
        for k, (i, j) in ok.items():
            s0, s1, strand = meta[k]
            S = segs[k][i:j]
            c0, c1 = (s0 + i, s0 + j) if strand == '+' else (s1 - j, s1 - i)
            o['_S'][k] = (S, c0, c1, strand)
        o['_L'], o['_placed'] = L, set(segs)
        out.append(o)
        if n % 500 == 0:
            print('anchors', n, len(loci), flush=True)
    # dipcall records overlapping [a, b) (one pass)
    iv = collections.defaultdict(list)
    for o in out:
        o['_recs'] = []
        if o['_R'] is not None:
            iv[o['_L']['chrom']].append((o['anchor_a0'], o['anchor_b0'], o))
    for c in iv:
        iv[c].sort(key=lambda z: z[0])
    st = {c: [z[0] for z in v] for c, v in iv.items()}
    spn = max((b - a for v in iv.values() for a, b, _ in v), default=0) + 1
    for rec in pysam.VariantFile(f'{DIP}.vcf.gz'):
        v = iv.get(rec.chrom)
        if not v:
            continue
        q0, q1 = rec.pos - 1, rec.pos - 1 + len(rec.ref)
        for a, b, o in v[bisect.bisect_left(st[rec.chrom], q0 - spn):bisect.bisect_left(st[rec.chrom], q1)]:
            if q0 < b and q1 > a:
                gt = rec.samples[0]['GT']
                o['_recs'].append((q0, rec.ref.upper(), tuple((x or '').upper() for x in rec.alts or ()), gt,
                                   not rec.filter.keys() or 'PASS' in rec.filter.keys()))
    rows = []
    for o in out:
        L, R, A = o['_L'], o['_R'], o['_A']
        ev = L['asm_event'].split(':')[0] if L['asm_event'] else ''
        evh = ev.split('_')[-1] if ev else ''
        seqs = {}
        for h, hk in (('hap1', 0), ('hap2', 1)):
            key = ('normal', f"{L['chrom']}_{h}")
            S = o['_S'].get(key)
            if R is None:
                ds, why, eds = None, 'no_anchors', []
            else:
                ds, why, eds = dip_hap(R, o['anchor_a0'], o['anchor_b0'], o['_recs'], hk, o['_inbed'])
            fl = [e for e in eds if e[0] + len(e[1]) <= o['arr_lo0'] or e[0] >= o['arr_hi0']]
            Af = edit(R, o['anchor_a0'], [(int(L['vcf_pos']) - 1, L['vcf_ref'].upper(), L['vcf_alt'].upper())] + fl) if fl and R else None
            call = lambda s: 'NA' if s is None else 'ALT' if s == A else 'ALT_plus_flank' if Af and s == Af else 'REF' if s == R else 'other'
            sa = S[0] if S else None
            o.update({f'{h}_asm_call': call(sa), f'{h}_asm_len_change': len(sa) - len(R) if sa is not None else '',
                      f'{h}_asm_na': '' if sa is not None else 'no_anchors' if R is None else 'not_placed' if key not in o['_placed'] else 'anchor',
                      f'{h}_dip_call': call(ds), f'{h}_dip_len_change': len(ds) - len(R) if ds is not None else '', f'{h}_dip_na': why,
                      f'{h}_asm_eq_dip': (sa == ds) if sa is not None and ds is not None else '',
                      f'{h}_n_dip_records': len(eds), f'{h}_contig_span': f'{S[1]}-{S[2]}{S[3]}' if S else '',
                      f'{h}_seq': sa if sa is not None else ds if ds is not None else ''})
            seqs[h] = sa if sa is not None else ds
        tum = [(k[1], v[0]) for k, v in sorted(o['_S'].items()) if k[0] == 'tumor']
        cls = []
        for c, s in tum:
            cl = 'NA' if None in seqs.values() else ('hom' if s == seqs['hap1'] == seqs['hap2'] else 'hap1' if s == seqs['hap1'] else
                                                    'hap2' if s == seqs['hap2'] else 'novel')
            cls.append((c, cl, s == A, len(s) - len(R)))
        if not cls or any(x[1] == 'NA' for x in cls):
            ver = 'NA'
        elif any(x[1] == 'novel' for x in cls):
            ver = 'novel_ALT' if any(x[1] == 'novel' and x[2] for x in cls) else 'novel_other'
        else:
            k = {x[1] for x in cls}
            ver = 'no_change_hom' if k == {'hom'} else 'no_change_both' if k == {'hap1', 'hap2'} else f'all_eq_{k.pop()}' if len(k) == 1 else 'no_change_both'
        o.update(n_tumor_contigs=len([k for k in o['_S'] if k[0] == 'tumor']),
                 tumor_copies=';'.join(f'{c}:{cl}:{"ALT" if al else "notALT"}:{dl:+d}' for c, cl, al, dl in cls),
                 tumor_any_ALT_array=any(x[2] for x in cls) if cls else '', tumor_verdict=ver)
        # truth INFO event applied to the event hap's anchored stretch
        tic = 'NA'
        if evh in ('hap1', 'hap2') and ('normal', ev) in o['_S']:
            S, c0, c1, strand = o['_S'][('normal', ev)]
            m = re.match(r'(\d+)-([ACGTN]*)-([ACGTN]*)$', L['asm_event'].split(':', 1)[1].upper())
            q0, er, ea = int(m.group(1)), m.group(2), m.group(3)
            if c0 <= q0 and q0 + len(er) <= c1:
                fw = S if strand == '+' else rc(S); off = q0 - c0
                if fw[off:off + len(er)] == er:
                    T = fw[:off] + ea + fw[off + len(er):]; T = T if strand == '+' else rc(T)
                    oth = seqs['hap2' if evh == 'hap1' else 'hap1']
                    tic = 'eq_ALT=other_hap' if T == A and T == oth else 'eq_ALT' if T == A else 'eq_other_hap' if T == oth else 'novel'
                    o['T_info_len_change'] = len(T) - len(R)
                else:
                    tic = 'event_ref_mismatch'
            else:
                tic = 'event_outside_stretch'
        o.update(event_hap=evh, T_info_class=tic)
        rows.append({k: v for k, v in o.items() if not k.startswith('_')})
    cols = []
    for r in rows:
        cols += [c for c in r if c not in cols]
    if test:
        for r in rows:
            print('\n'.join(f'  {c}: {str(r.get(c, ""))[:400]}' for c in cols)); print()
        return
    with open(f'{D}/array_alleles.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for r in rows:
            w.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    C = collections.Counter
    print('loci', len(rows), 'anchored', sum(r['anchor_b0'] != '' for r in rows), 'capped', sum(r['arr_capped'] for r in rows))
    for h in ('hap1', 'hap2'):
        print(h, 'asm', C(r[f'{h}_asm_call'] for r in rows), 'dip', C(r[f'{h}_dip_call'] for r in rows))
        print(h, 'asm_eq_dip', C(r[f'{h}_asm_eq_dip'] for r in rows), 'dip_na', C(r[f'{h}_dip_na'] for r in rows))
    print('tumor', C(r['tumor_verdict'] for r in rows), 'T_info', C(r['T_info_class'] for r in rows))


if __name__ == '__main__':
    main()
