"""Audit b5 (biology lens): mutational context of every truth allele in the patient's own genome (HG008-N v6.2 frame).

inputs   b1_truth_context.tsv (asm_events = every HG008Nv62SOMATICVARIANT event of the truth allele: contig chrN_hapK,
         0-based position, REF, ALT), HG008-N v6.2 hap1 / hap2 FASTA (indexed).
method   per event: check that the assembly has REF at the position (ref_ok); kind SNV / MNV / INDEL; INDEL unit = minimal
         period of the inserted / deleted bases after the shared prefix; repeat tract = longest exact period-p stretch
         (p = unit if <= 6, else 1..6) touching the event in the assembly; class HP>=7 / HP4-6 / STR2-6>=3copies /
         VNTR>6 / none, as in b1 but in the patient's sequence; SNV CpG transition (C>T before G, G>A after C).
         The truth row takes its first event with ref_ok (all events are listed).
outputs  b5_normal_frame.tsv: truth_id, n_events, nf_event, nf_ref_ok, nf_kind, nf_len, nf_unit, nf_tract, nf_class,
         nf_cpg_ti, nf_titv.
assumes  the INFO positions are 0-based (memory note; checked here by ref_ok).
"""
import csv, re
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
O = f'{D}/audit/biology'
G = '/scratch/jshen/data/HG008_GIAB'
FA = {'hap1': pysam.FastaFile(f'{G}/HG008N_curatedv6_250714_polished6.2_hap1.fasta'),
      'hap2': pysam.FastaFile(f'{G}/HG008N_curatedv6_250714_polished6.2_hap2.fasta.gz')}


def tract(seq, lo, hi, periods):
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


def classify(e):
    m = re.match(r'(.+)_(hap[12]):(\d+)-([^-]*)-(.*)', e)
    if not m: return None
    ctg, hp, p0, ref, alt = m.group(1) + '_' + m.group(2), m.group(2), int(m.group(3)), m.group(4).upper(), m.group(5).upper()
    fa = FA[hp]
    if ctg not in fa.references: return dict(nf_event=e, nf_ref_ok='no_contig')
    s0 = max(0, p0 - 300); seq = fa.fetch(ctg, s0, p0 + len(ref) + 300).upper(); o = p0 - s0
    out = dict(nf_event=e, nf_ref_ok=seq[o:o + len(ref)] == ref)
    if len(ref) == 1 and len(alt) == 1:
        p, a, b = tract(seq, o, o + 1, range(1, 7)); L = b - a
        out.update(nf_kind='SNV', nf_len=0, nf_unit=0, nf_tract=L,
                   nf_class='HP>=7' if p == 1 and L >= 7 else 'HP4-6' if p == 1 and L >= 4 else
                   'STR2-6>=3copies' if p >= 2 and L >= 3 * p else 'none',
                   nf_cpg_ti=(ref == 'C' and alt == 'T' and seq[o + 1] == 'G') or (ref == 'G' and alt == 'A' and seq[o - 1] == 'C'),
                   nf_titv='Ti' if {ref, alt} in ({'A', 'G'}, {'C', 'T'}) else 'Tv')
        return out
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]: k += 1
    r2, a2 = ref[k:], alt[k:]
    if r2 and a2:
        out.update(nf_kind='MNV' if len(ref) == len(alt) else 'complex', nf_len=len(alt) - len(ref), nf_unit=0, nf_tract=0, nf_class='complex')
        return out
    S = r2 or a2; u = minper(S); lo = o + k; hi = o + k + len(r2)
    p, a, b = tract(seq, lo - 1 if not r2 else lo, hi + 1 if not r2 else hi, [u] if u <= 6 else range(1, 7))
    L = b - a
    out.update(nf_kind='INDEL', nf_len=len(a2) - len(r2), nf_unit=u, nf_tract=L,
               nf_class='HP>=7' if u == 1 and L >= 7 else 'HP4-6' if u == 1 and L >= 4 else
               'STR2-6>=3copies' if 2 <= u <= 6 and L >= 3 * u else
               'VNTR>6' if u > 6 and (seq[hi:hi + len(S)] == S or seq[lo - len(S):lo] == S) else
               ('in_STR(p%d)' % p if L >= 10 else 'none'))
    return out


cols = ['truth_id', 'n_events', 'nf_event', 'nf_ref_ok', 'nf_kind', 'nf_len', 'nf_unit', 'nf_tract', 'nf_class', 'nf_cpg_ti', 'nf_titv']
with open(f'{O}/b5_normal_frame.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for r in csv.DictReader(open(f'{O}/b1_truth_context.tsv'), delimiter='\t'):
        evs = [e for e in r['asm_events'].split(',') if e]
        res = [x for x in (classify(e) for e in evs) if x]
        best = next((x for x in res if x.get('nf_ref_ok') is True), res[0] if res else {})
        o = dict(truth_id=r['truth_id'], n_events=len(evs), **best)
        w.write('\t'.join(str(o.get(k, '')) for k in cols) + '\n')
