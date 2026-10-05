"""COLO829T step 1: GRCh38-frame repeat context of EVERY chr1-22 COLO829T truth allele (SNV and INDEL, 43,947 alleles).

Port of the GRCh38-frame part of analysis/graph_absorbed_somatic_20261001/audit/b1_truth_context.py (functions tract, minper,
context copied unchanged) plus the s1_loci.py event window rule (event_window copied unchanged). The HG008 normal-frame /
dipcall / GIAB-BED / PoN parts of b1 are not ported (COLO829T truth has no patient-frame INFO event).
inputs   analysis/tensor_recall_20260930/per_truth_COLO829T.tsv (all truth alleles + per-platform status / class),
         $D/variant_set.tsv (c0: the 497 absorbed = perfect-bypass truths, per-platform *_perfect), GRCh38 FASTA.
method   context in GRCh38 (b1): INDEL: shared prefix k of REF/ALT, inserted/deleted bases S, unit u = minimal period of S,
         indel_len = len(ALT) - len(REF) after the prefix; tract L = longest exact period-u stretch touching the event
         (u <= 6; for a pure insertion the point is widened by 1 bp each side), for u > 6 the copies of S right of the event;
         ctx_class  HP>=7 (u 1, L >= 7) / HP4-6 (u 1, L 4-6) / STR2-6>=3copies (u 2-6, L >= 3u) / VNTR>6 (u > 6, L >= 2u) /
         complex (REF and ALT both change after the prefix) / none; ctx_period / ctx_tract 0 when no stretch is found.
         SNV: the longest period 1-6 stretch (>= 2 periods) covering the base gives ctx_period / ctx_tract and the same
         classes (HP>=7 / HP4-6 / STR2-6>=3copies / none); cpg (C with 3' G or G with 5' C), cpg_ti (C>T / G>A at CpG),
         titv, alt_hp_len (homopolymer of the ALT base the SNV creates), tri (GRCh38 trinucleotide), alt_extends_hp>=5.
         Event window (s1): union of every period 1-6 stretch (length >= max(2p, p + 3)) touching the event, +-5 bp ->
         ev_lo0 / ev_hi0 (0-based half-open GRCh38).
         b5rule_class / _period / _tract: the patient-frame rule of audit/b5_normal_frame.py applied in GRCh38 (extra, for a
         later same-rule comparison with the COLO829BL frame): as ctx_class, except an INDEL with u > 6 is VNTR>6 only if S
         sits right next to the event, else in_STR(pP) when the longest period 1-6 stretch touching it is >= 10 bp.
outputs  $D/grch38_context.tsv (one row per truth allele, per_truth order; absorbed / n_platforms / *_perfect from c0;
         indel_seq = the inserted (INS) / deleted (DEL) bases), printed tables (class x kind, absorbed share).
check    `c1_grch38_context.py --check-hg008`: runs the same functions on 20 HG008T truths (2 per kind x ctx_class, the
         one complex INDEL, 1 more absorbed truth; seed 1) and compares every context column with
         audit/biology/b1_truth_context.tsv and ev_lo0 / ev_hi0 with loci.tsv of the HG008 run where the truth is in it;
         b5_rule on 2 HG008-N v6.2 events per nf_kind x nf_class vs audit/biology/b5_normal_frame.tsv.
assumes  per_truth truth_id == variant_set truth_id (c0 checks the same keys); truth REF/ALT are upper-case ACGT.
"""
import csv, collections, random, re, sys
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
R = '/scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930'
H = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
fa = pysam.FastaFile('/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
PLAT = ('fiberseq', 'ONT', 'Illumina')
CTX = ['indel_unit', 'indel_len', 'ctx_period', 'ctx_tract', 'ctx_class', 'cpg', 'cpg_ti', 'titv', 'alt_hp_len', 'tri',
       'alt_extends_hp>=5']


def tract(seq, lo, hi, periods=range(1, 7)):                          # b1, unchanged
    """longest exact periodic stretch (period p) touching [lo, hi) in seq: (p, start, end)."""
    best = (0, lo, lo)
    for p in periods:
        i = max(0, lo - 200)
        while i < min(len(seq) - p, hi + 200):
            if seq[i] != seq[i + p]: i += 1; continue
            j = i
            while j < len(seq) - p and seq[j] == seq[j + p]: j += 1
            a, b = i, j + p
            if a <= hi and b >= lo and b - a > best[2] - best[1] and b - a >= 2 * p: best = (p, a, b)
            i = j + 1
    return best


def minper(s):                                                        # b1, unchanged
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s: return p
    return len(s)


def context(c, pos0, ref, alt):                                       # b1, unchanged
    F = 300; s0 = max(0, pos0 - F); seq = fa.fetch(c, s0, pos0 + len(ref) + F).upper(); o = pos0 - s0
    out = {}
    if len(ref) == 1 and len(alt) == 1:
        p, a, b = tract(seq, o, o + 1)
        out.update(ctx_period=p, ctx_tract=b - a)
        out['cpg'] = (ref == 'C' and seq[o + 1] == 'G') or (ref == 'G' and seq[o - 1] == 'C')
        out['cpg_ti'] = (ref == 'C' and alt == 'T' and seq[o + 1] == 'G') or (ref == 'G' and alt == 'A' and seq[o - 1] == 'C')
        out['titv'] = 'Ti' if {ref, alt} in ({'A', 'G'}, {'C', 'T'}) else 'Tv'
        l = 0
        while seq[o - 1 - l] == alt: l += 1
        r = 0
        while seq[o + 1 + r] == alt: r += 1
        out['alt_hp_len'] = l + r + 1
        out['tri'] = seq[o - 1:o + 2]
        P, L = p, b - a
        out['ctx_class'] = ('HP>=7' if P == 1 and L >= 7 else 'HP4-6' if P == 1 and L >= 4 else
                            'STR2-6>=3copies' if P >= 2 and L >= 3 * P else 'none')
        out['alt_extends_hp>=5'] = out['alt_hp_len'] >= 5
        return out
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]: k += 1
    r2, a2 = ref[k:], alt[k:]
    if r2 and a2:
        out.update(ctx_class='complex', indel_unit=0, ctx_period=0, ctx_tract=0, indel_len=len(alt) - len(ref)); return out
    S = r2 or a2; u = minper(S); out['indel_unit'] = u; out['indel_len'] = len(a2) - len(r2)
    lo = o + k; hi = o + k + len(r2)
    if u <= 6:
        p, a, b = tract(seq, lo - 1 if not r2 else lo, hi + 1 if not r2 else hi, periods=[u])
    else:
        p, a, b = (u, lo, hi)
        n = 0
        while seq[hi + n * u: hi + (n + 1) * u] == S: n += 1
        b = hi + n * u
    L = b - a
    out.update(ctx_period=p, ctx_tract=L)
    out['ctx_class'] = ('HP>=7' if u == 1 and L >= 7 else 'HP4-6' if u == 1 and L >= 4 else
                        'STR2-6>=3copies' if 2 <= u <= 6 and L >= 3 * u else
                        'VNTR>6' if u > 6 and L >= 2 * u else 'none')
    return out


def b5_rule(seq, o, ref, alt):
    """audit/b5_normal_frame.py classify() rule on a sequence (INDEL part, unchanged): u > 6 -> VNTR>6 if S is next to the
    event, else in_STR(pP) if the longest period 1-6 stretch touching it is >= 10 bp. SNVs: same rule as b1."""
    p, a, b = 0, 0, 0
    if len(ref) == 1 and len(alt) == 1:
        p, a, b = tract(seq, o, o + 1, range(1, 7)); L = b - a
        return ('HP>=7' if p == 1 and L >= 7 else 'HP4-6' if p == 1 and L >= 4 else
                'STR2-6>=3copies' if p >= 2 and L >= 3 * p else 'none'), p, L
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]: k += 1
    r2, a2 = ref[k:], alt[k:]
    if r2 and a2: return 'complex', 0, 0
    S = r2 or a2; u = minper(S); lo = o + k; hi = o + k + len(r2)
    p, a, b = tract(seq, lo - 1 if not r2 else lo, hi + 1 if not r2 else hi, [u] if u <= 6 else range(1, 7))
    L = b - a
    return ('HP>=7' if u == 1 and L >= 7 else 'HP4-6' if u == 1 and L >= 4 else
            'STR2-6>=3copies' if 2 <= u <= 6 and L >= 3 * u else
            'VNTR>6' if u > 6 and (seq[hi:hi + len(S)] == S or seq[lo - len(S):lo] == S) else
            ('in_STR(p%d)' % p if L >= 10 else 'none')), p, L


def event_window(chrom, pos0, vref, valt):                            # s1_loci.py, unchanged
    lref = len(vref)
    if lref == 1 and len(valt) == 1:
        e_lo, e_hi = pos0, pos0 + 1
    else:
        e_lo, e_hi = pos0 + 1, max(pos0 + lref, pos0 + 1)
    s0 = max(0, e_lo - 300)
    seq = fa.fetch(chrom, s0, e_hi + 300).upper()
    lo, hi = e_lo - s0, e_hi - s0
    L, R = lo, hi
    for p in range(1, 7):
        i = 0
        while i < len(seq) - p:
            if seq[i] != seq[i + p]:
                i += 1; continue
            j = i
            while j < len(seq) - p and seq[j] == seq[j + p]:
                j += 1
            a, b = i, j + p
            if b - a >= max(2 * p, p + 3) and a <= hi + 1 and b >= lo - 1:
                L, R = min(L, a), max(R, b)
            i = j + 1
    return s0 + L - 5, s0 + R + 5


def check_hg008(n=20):
    b1 = list(csv.DictReader(open(f'{H}/audit/biology/b1_truth_context.tsv'), delimiter='\t'))
    loci = {r['truth_id']: r for r in csv.DictReader(open(f'{H}/loci.tsv'), delimiter='\t')}
    by = collections.defaultdict(list)
    for r in b1: by[(r['kind2'], r['ctx_class'])].append(r)
    rng = random.Random(1); pick = []
    for k in sorted(by):                     # 2 per kind x class (complex has 1), absorbed ones first (they have s1 windows)
        g = [r for r in by[k] if r['truth_id'] in loci]; g = g if len(g) >= 2 else by[k]
        pick += rng.sample(g, min(2, len(g)))
    ids = {r['truth_id'] for r in pick}
    pick += rng.sample([r for r in b1 if r['truth_id'] in loci and r['truth_id'] not in ids], max(0, n - len(pick)))
    bad = 0
    for r in pick:
        o = context(r['chrom'], int(r['vcf_pos']) - 1, r['vcf_ref'], r['vcf_alt'])
        diff = [k for k in CTX if str(o.get(k, '')) != r[k]]
        ev = event_window(r['chrom'], int(r['vcf_pos']) - 1, r['vcf_ref'], r['vcf_alt']); L = loci.get(r['truth_id'])
        evd = L is not None and (str(ev[0]), str(ev[1])) != (L['ev_lo0'], L['ev_hi0'])
        bad += bool(diff) or evd
        print(r['truth_id'], r['kind2'], r['vcf_ref'][:12], r['vcf_alt'][:12], 'b1', r['ctx_class'], r['indel_unit'],
              r['ctx_tract'], 'port', o['ctx_class'], o.get('indel_unit', ''), o['ctx_tract'], 'diff', diff or '-',
              'ev', ev, 's1', (L['ev_lo0'], L['ev_hi0']) if L else 'n/a')
    print(f'{len(pick)} HG008T truths, {bad} differ;', collections.Counter((r['kind2'], r['ctx_class']) for r in pick))
    # b5_rule: on HG008-N v6.2 patient-frame events, vs b5_normal_frame.tsv (2 per nf_kind x nf_class)
    G = '/scratch/jshen/data/HG008_GIAB'
    FA = {'hap1': pysam.FastaFile(f'{G}/HG008N_curatedv6_250714_polished6.2_hap1.fasta'),
          'hap2': pysam.FastaFile(f'{G}/HG008N_curatedv6_250714_polished6.2_hap2.fasta.gz')}
    b5 = [r for r in csv.DictReader(open(f'{H}/audit/biology/b5_normal_frame.tsv'), delimiter='\t') if r['nf_ref_ok'] == 'True']
    by = collections.defaultdict(list)
    for r in b5: by[(r['nf_kind'], r['nf_class'].split('(')[0])].append(r)
    bad = 0; pick = [x for k in sorted(by) for x in rng.sample(by[k], min(2, len(by[k])))]
    for r in pick:
        m = re.match(r'(.+)_(hap[12]):(\d+)-([^-]*)-(.*)', r['nf_event']); ctg, p0 = m.group(1) + '_' + m.group(2), int(m.group(3))
        ref, alt = m.group(4).upper(), m.group(5).upper(); s0 = max(0, p0 - 300)
        cl, _, L = b5_rule(FA[m.group(2)].fetch(ctg, s0, p0 + len(ref) + 300).upper(), p0 - s0, ref, alt)
        bad += (cl, str(L)) != (r['nf_class'], r['nf_tract']) and r['nf_kind'] in ('SNV', 'INDEL')
        print(r['truth_id'], r['nf_event'][:50], 'b5', r['nf_class'], r['nf_tract'], 'port', cl, L)
    print(f'b5_rule: {len(pick)} HG008-N events, {bad} differ (SNV / INDEL; MNV / complex give complex 0 here)')


if __name__ == '__main__' and '--check-hg008' in sys.argv:
    check_hg008(); sys.exit()
vs = {r['truth_id']: r for r in csv.DictReader(open(f'{D}/variant_set.tsv'), delimiter='\t')}
rows = list(csv.DictReader(open(f'{R}/per_truth_COLO829T.tsv'), delimiter='\t'))
cols = (['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'ins_del', 'indel_seq', 'in_bed', 'RGN', 'VAF_Ill',
         'VAF_PB', 'absorbed', 'n_platforms'] + [f'{p}_perfect' for p in PLAT] + [f'{p}_{x}' for p in PLAT for x in ('status', 'class')]
        + CTX + ['b5rule_class', 'b5rule_period', 'b5rule_tract', 'ev_lo0', 'ev_hi0', 'ev_len'])
out = []
with open(f'{D}/grch38_context.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for t in rows:
        c, pos0, ref, alt = t['chrom'], int(t['vcf_pos']) - 1, t['vcf_ref'], t['vcf_alt']
        assert fa.fetch(c, pos0, pos0 + len(ref)).upper() == ref, t['truth_id']
        snv = len(ref) == 1 and len(alt) == 1; v = vs.get(t['truth_id'])
        o = dict(t, kind2='SNV' if snv else 'INDEL', absorbed=v is not None, n_platforms=v['n_platforms'] if v else 0,
                 ins_del='' if snv else 'INS' if len(alt) > len(ref) else 'DEL' if len(alt) < len(ref) else 'MNV',
                 indel_seq='' if snv else (alt[1:] if len(alt) > len(ref) else ref[1:]) if ref[0] == alt[0] else '')
        for p in PLAT: o[f'{p}_perfect'] = v[f'{p}_perfect'] if v else False
        o.update(context(c, pos0, ref, alt))
        s0 = max(0, pos0 - 300); o['b5rule_class'], o['b5rule_period'], o['b5rule_tract'] = b5_rule(
            fa.fetch(c, s0, pos0 + len(ref) + 300).upper(), pos0 - s0, ref, alt)
        o['ev_lo0'], o['ev_hi0'] = event_window(c, pos0, ref, alt); o['ev_len'] = o['ev_hi0'] - o['ev_lo0']
        w.write('\t'.join(str(o.get(k, '')) for k in cols) + '\n'); out.append(o)
print('rows', len(out), 'absorbed', sum(o['absorbed'] for o in out), 'of', len(vs))
ORDER = ['HP>=7', 'HP4-6', 'STR2-6>=3copies', 'VNTR>6', 'complex', 'none']
for kd in ('SNV', 'INDEL'):
    sel = [o for o in out if o['kind2'] == kd]
    print(f'\n{kd}: class | all chr1-22 | absorbed (union) | % absorbed | ' + ' | '.join(f'{p} perfect' for p in PLAT))
    for k in ORDER + ['total']:
        s = [o for o in sel if k == 'total' or o['ctx_class'] == k]
        if not s: continue
        a = sum(o['absorbed'] for o in s)
        print(f'{k} | {len(s):,} | {a} | {100 * a / len(s):.1f}% | ' + ' | '.join(str(sum(o[f"{p}_perfect"] == 'True' for o in s)) for p in PLAT))
I = [o for o in out if o['kind2'] == 'INDEL']
print('\nINDEL ctx_class x b5rule_class (all / absorbed):', {k: (v, sum(o['absorbed'] for o in I if (o['ctx_class'], o['b5rule_class'].split('(')[0]) == k))
      for k, v in sorted(collections.Counter((o['ctx_class'], o['b5rule_class'].split('(')[0]) for o in I).items())})

# Results (run 2026-10-05, login node: 469 s wall, 152 s CPU, max RSS 0.19 GB)
# --check-hg008: 20 HG008T truths (b1 classes INDEL HP>=7 3, HP4-6 2, STR 2, VNTR>6 2, none 2, complex 1; SNV HP>=7 2,
#   HP4-6 2, STR 2, none 2): 0 differ in the 11 context columns, ev_lo0 / ev_hi0 equal to s1 loci.tsv for all 20;
#   b5_rule on 20 HG008-N v6.2 events (incl. VNTR>6, in_STR(p1), in_STR(p2)): 0 differ in class and tract.
# grch38_context.tsv: 43,947 rows (SNV 42,035, INDEL 1,912 = DEL 1,024 + INS 888; no MNV / complex), absorbed 497 of 497;
#   per-platform perfect SNV / INDEL = fiberseq 64 / 218, ONT 58 / 170, Illumina 53 / 396 (= c0).
# GRCh38 frame (b1 ctx_class): all chr1-22 truth | absorbed union (% of class) | fiberseq / ONT / Illumina perfect
#   SNV   HP>=7            545 |   4 (0.7%)  |  4 /  3 /   2
#         HP4-6          5,914 |   7 (0.1%)  |  7 /  6 /   7
#         STR2-6>=3copies 1,322 |  14 (1.1%)  | 13 /  9 /  11
#         none          34,254 |  42 (0.1%)  | 40 / 40 /  33
#         total         42,035 |  67 (0.2%)  | 64 / 58 /  53
#   INDEL HP>=7          1,167 | 244 (20.9%) | 81 / 47 / 238
#         HP4-6             97 |   3 (3.1%)  |  2 /  3 /   3
#         STR2-6>=3copies   275 | 120 (43.6%) | 86 / 71 /  99
#         VNTR>6              4 |   2 (50.0%) |  2 /  2 /   2
#         none             369 |  61 (16.5%) | 47 / 47 /  54
#         total          1,912 | 430 (22.5%) | 218 / 170 / 396
# INDEL b1 'none' by the b5 rule in GRCh38 (all / absorbed): in_STR 79 / 32 (unit > 6 inside a period 1-6 stretch >= 10 bp,
#   mostly (TA)n), VNTR>6 10 / 0, none 280 / 29.
