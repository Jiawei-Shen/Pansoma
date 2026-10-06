"""Audit b1 (biology lens): context, normal-frame event and HG008-N dipcall identity for ALL HG008T truth alleles.

Independent of s0-s9 (reads only their outputs to label which truths are absorbed / which category):
inputs   per_truth_HG008T.tsv (all 17,186 truth alleles + statuses), truth VCF (all HG008Nv62SOMATICVARIANT events, not
         only the first one that s1 kept), GIAB v0.2 all / nogermlineoverlap BEDs, GRCh38 FASTA, the bcftools-normalised
         indexed copy of the HG008-N dipcall VCF (records with FILTER PASS/'.'), dip.bed, truth_pon_HG008T.tsv,
         $D/per_variant.tsv (category, HPRC numbers of the absorbed set).
method   (1) context in GRCh38: maximal period-1..6 tract covering the event (seq[i]==seq[i+p] runs); INDEL unit = minimal
         period of the inserted/deleted sequence; SNV CpG (C>T with 3' G / G>A with 5' C), Ti/Tv;
         (2) truth INFO events: kind (SNV / INDEL / MNV) of every normal-assembly-frame event, and how many truth alleles
         share the same event string;
         (3) dipcall: window pos0-60 .. end+60 widened to every overlapping record; alt window = truth applied;
         ident_any = one dipcall allele applied alone gives the alt window; hapK_eq_alt = the whole phased hap K over
         the window equals the alt window; hapK_len_minus_alt = length difference (STR length vs the tumor allele).
outputs  b1_truth_context.tsv (one row per truth allele, per_truth order).
assumes  per_truth truth_id == variant_set truth_id; dipcall GT is phased hap1|hap2 (as s8 states).
"""
import csv, collections, re
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
O = f'{D}/audit/biology'
R = '/scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930'
G = '/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark'
fa = pysam.FastaFile('/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
dip = pysam.VariantFile('/scratch/jshen/tmp/HG008N_vs_HPRC_isec/HG008N_dipcall.clean.norm.sorted.vcf.gz')
key = lambda c, p, r, a: (c, int(p), r, a)


def bed(path):
    b = collections.defaultdict(list)
    for l in open(path):
        c, s, e = l.split()[:3]; b[c].append((int(s), int(e)))
    for c in b: b[c].sort()
    return b


import bisect
def inbed(b, c, s, e):                      # [s, e) fully inside one interval
    iv = b.get(c, []); i = bisect.bisect_right(iv, (s, 10**12)) - 1
    return i >= 0 and iv[i][0] <= s and e <= iv[i][1]


ALL, NOG = bed(f'{G}/HG008-T_somatic_smvar_benchmark_v0.2_all.bed'), bed(f'{G}/HG008-T_somatic_smvar_benchmark_v0.2_nogermlineoverlap.bed')
DIPB = bed('/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed')
events = {}
for rec in pysam.VariantFile(f'{G}/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz'):
    v = rec.info.get('HG008Nv62SOMATICVARIANT') or ()
    v = (v,) if isinstance(v, str) else v
    for a in rec.alts or ():
        events[key(rec.chrom, rec.pos, rec.ref, a)] = [e for e in v if e]
pon = {r['truth_id']: r for r in csv.DictReader(open(f'{R}/pon/truth_pon_HG008T.tsv'), delimiter='\t')}
pv = {r['truth_id']: r for r in csv.DictReader(open(f'{D}/per_variant.tsv'), delimiter='\t')}


def tract(seq, lo, hi, periods=range(1, 7)):
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


def minper(s):
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s: return p
    return len(s)


def context(c, pos0, ref, alt):
    F = 300; s0 = max(0, pos0 - F); seq = fa.fetch(c, s0, pos0 + len(ref) + F).upper(); o = pos0 - s0
    out = {}
    if len(ref) == 1 and len(alt) == 1:
        p, a, b = tract(seq, o, o + 1)
        out.update(ctx_period=p, ctx_tract=b - a)
        out['cpg'] = (ref == 'C' and seq[o + 1] == 'G') or (ref == 'G' and seq[o - 1] == 'C')
        out['cpg_ti'] = (ref == 'C' and alt == 'T' and seq[o + 1] == 'G') or (ref == 'G' and alt == 'A' and seq[o - 1] == 'C')
        out['titv'] = 'Ti' if {ref, alt} in ({'A', 'G'}, {'C', 'T'}) else 'Tv'
        # ALT base extends an adjacent homopolymer (left or right run of the ALT base)
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
        # copies of a longer unit next to the event
        n = 0
        while seq[hi + n * u: hi + (n + 1) * u] == S: n += 1
        b = hi + n * u
    L = b - a
    out.update(ctx_period=p, ctx_tract=L)
    out['ctx_class'] = ('HP>=7' if u == 1 and L >= 7 else 'HP4-6' if u == 1 and L >= 4 else
                        'STR2-6>=3copies' if 2 <= u <= 6 and L >= 3 * u else
                        'VNTR>6' if u > 6 and L >= 2 * u else 'none')
    return out


def dipcheck(c, pos0, ref, alt):
    ws, we = pos0 - 60, pos0 + len(ref) + 60
    recs = []
    for _ in range(3):                          # widen to overlapping records
        recs = [r for r in dip.fetch(c, max(0, ws), we) if (not r.filter.keys() or 'PASS' in r.filter.keys())]
        nws = min([ws] + [r.start for r in recs]); nwe = max([we] + [r.stop for r in recs])
        if (nws, nwe) == (ws, we): break
        ws, we = nws, nwe
    W = fa.fetch(c, ws, we).upper()
    altw = W[:pos0 - ws] + alt + W[pos0 - ws + len(ref):]
    ident = False
    for r in recs:
        for a in r.alts or ():
            if W[:r.start - ws] + a + W[r.stop - ws:] == altw: ident = True
    haps = []
    for h in (0, 1):
        sel, bad = [], False
        for r in recs:
            gt = r.samples[0]['GT']
            if len(gt) < 2: bad = True; continue
            g = gt[h]
            if g is None: bad = True; continue
            if g > 0: sel.append((r.start, r.stop, r.alts[g - 1]))
        sel.sort(); s = []; last = ws
        for a0, b0, al in sel:
            if a0 < last: bad = True; break
            s.append(W[last - ws:a0 - ws]); s.append(al); last = b0
        s.append(W[last - ws:]); hs = ''.join(s)
        haps.append(None if bad else hs)
    return dict(dip_n_records=len(recs), dip_ident_any=ident,
                dip_hap1_eq_alt=None if haps[0] is None else haps[0] == altw,
                dip_hap2_eq_alt=None if haps[1] is None else haps[1] == altw,
                dip_hap1_len_minus_alt='' if haps[0] is None else len(haps[0]) - len(altw),
                dip_hap2_len_minus_alt='' if haps[1] is None else len(haps[1]) - len(altw),
                in_dip_bed=inbed(DIPB, c, pos0, pos0 + len(ref)))


rows = list(csv.DictReader(open(f'{R}/per_truth_HG008T.tsv'), delimiter='\t'))
cols = None
with open(f'{O}/b1_truth_context.tsv', 'w') as w:
    for t in rows:
        c, pos0, ref, alt = t['chrom'], int(t['vcf_pos']) - 1, t['vcf_ref'], t['vcf_alt']
        kd = 'SNV' if len(ref) == 1 and len(alt) == 1 else 'INDEL'
        o = dict(truth_id=t['truth_id'], chrom=c, vcf_pos=t['vcf_pos'], vcf_ref=ref, vcf_alt=alt, kind2=kd,
                 in_bed=t['in_bed'], absorbed=t['truth_id'] in pv,
                 category=pv[t['truth_id']]['category'] if t['truth_id'] in pv else '',
                 sub_class=pv[t['truth_id']]['sub_class'] if t['truth_id'] in pv else '')
        ev = events.get(key(c, t['vcf_pos'], ref, alt), [])
        ks = set()
        for e in ev:
            m = re.match(r'(.+):(\d+)-([^-]+)-(.+)', e)
            if m: ks.add('SNV' if len(m.group(3)) == 1 == len(m.group(4)) else 'MNV' if len(m.group(3)) == len(m.group(4)) else 'INDEL')
        o.update(asm_events=','.join(ev), n_asm_events=len(ev), asm_kinds='+'.join(sorted(ks)) or 'none',
                 in_all_bed=inbed(ALL, c, pos0, pos0 + len(ref)), in_nogermline_bed=inbed(NOG, c, pos0, pos0 + len(ref)))
        o.update(context(c, pos0, ref, alt))
        o.update(dipcheck(c, pos0, ref, alt))
        p = pon.get(t['truth_id'], {})
        o.update(pon_rule=p.get('rule', ''), gnomAD_AF=p.get('gnomAD_AF', ''), CoLoRSdb_AF=p.get('CoLoRSdb_AF', ''),
                 dbSNP_AF=p.get('dbSNP_AF', ''), G1000_AF=p.get('1000G_AF', ''))
        for P in ('PacBio', 'ONT', 'Illumina'): o[f'{P}_status'] = t[f'{P}_status']
        if cols is None:
            cols = list(o) + ['ctx_period', 'ctx_tract', 'ctx_class', 'cpg', 'cpg_ti', 'titv', 'alt_hp_len', 'tri',
                              'alt_extends_hp>=5', 'indel_unit', 'indel_len']
            cols = list(dict.fromkeys(cols)); w.write('\t'.join(cols) + '\n')
        w.write('\t'.join(str(o.get(k, '')) for k in cols) + '\n')
# shared events
rows = list(csv.DictReader(open(f'{O}/b1_truth_context.tsv'), delimiter='\t'))
cnt = collections.Counter(e for r in rows for e in set(r['asm_events'].split(',')) if e)
with open(f'{O}/b1_truth_context.tsv', 'w') as w:
    cols = list(rows[0]) + ['max_truths_sharing_event']
    w.write('\t'.join(cols) + '\n')
    for r in rows:
        r['max_truths_sharing_event'] = max([cnt[e] for e in r['asm_events'].split(',') if e] or [0])
        w.write('\t'.join(str(r[k]) for k in cols) + '\n')
print('rows', len(rows))
