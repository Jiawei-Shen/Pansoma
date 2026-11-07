"""E5e step 4: per-locus verdicts and tables (COLO829T, HPRC v2.1 d46, normal-absent absorbed alleles).

Tumor assembly (e2 hits): each tumor copy T (every exact-anchor hit in hifiasm hap1 / hap2) against the two COLO829BL
sequences of c3 (hapX, hapY; None when the hap allele is NA): hom (both haps equal and T equals them) / hapX / hapY /
novel. Verdict = E5d s8b step 6, line for line: NA (no copy or a normal sequence missing) / novel_ALT (a novel copy == A)
/ novel_other / no_change_hom / no_change_both / all_eq_<hap>. The HG008T row of E5d counts 'tumor contig in neither
normal hap' = novel_ALT + novel_other and 'tumor contig = GRCh38 + truth' = novel_ALT over all normal-absent INDELs.
Second view: the R3.7 per-tumor-hap rule of altlen_tumor_asm (truth / normal / novel per tumor hap -> outcome).
Check: the 22 ALT-length loci of E27 v3 colo829t_altlen_tumor_asm (same method) must give the same hits and outcomes.
Tumor reads (e3 b4.tsv): EXP seen / loci with status ok, per platform, absorbed and control. Strict variant (not in
E5d): only loci whose EXP length differs from both normal dipcall lengths (NE_len, NO_len), where a length match cannot
come from a germline allele; computed for HG008T too from E5d's b4_all.tsv and the audit's b2 rows.
Outputs: tables.md, summary.json (this folder); $DATA/per_locus.tsv.
Run: python e4_summary.py   (login node, seconds)
"""
import collections, csv, glob, json
from pathlib import Path

csv.field_size_limit(10 ** 9)
HERE = Path(__file__).resolve().parent
DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E5e_colo829t_tumor_evidence'
NF = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005/normal_frame.tsv'
A22 = '/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E27_assessment_20261030/v3/colo829t_altlen_tumor_asm'
E5D = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46'
BIO = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/audit/biology'
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))  # noqa: E731


def pct(a, b):
    return f'{a:,} / {b:,} ({100 * a / b:.1f}%)' if b else '-'


def asm_verdicts(loci):
    nf = {r['truth_id']: r for r in rd(NF)}
    hits = collections.defaultdict(list)
    for f in sorted(glob.glob(f'{DATA}/tumor_asm/hits.*.tsv')):
        for h in rd(f):
            hits[h['truth_id']].append(h)
    out = {}
    for L in loci:
        t = L['truth_id']
        if L['in_nf'] != 'True' or L['set'] not in ('absent_INDEL', 'absent_SNV'):
            continue
        n = nf[t]
        A = n['A_seq']
        seqs = {h: (None if n[f'{h}_allele'] in ('NA', '') else n[f'{h}_seq']) for h in ('hapX', 'hapY')}
        copies = [h for h in hits[t] if h['tumor_hap'] != 'none']
        cls = []
        for h in copies:
            s = h['T_seq']
            cl = 'NA' if None in seqs.values() else ('hom' if s == seqs['hapX'] == seqs['hapY'] else 'hapX' if s == seqs['hapX']
                                                     else 'hapY' if s == seqs['hapY'] else 'novel')
            cls.append((h['tumor_hap'], cl, s == A))
        if not cls or any(x[1] == 'NA' for x in cls):
            ver = 'NA'
        elif any(x[1] == 'novel' for x in cls):
            ver = 'novel_ALT' if any(x[1] == 'novel' and x[2] for x in cls) else 'novel_other'
        else:
            k = {x[1] for x in cls}
            ver = 'no_change_hom' if k == {'hom'} else 'no_change_both' if k == {'hapX', 'hapY'} else f'all_eq_{k.pop()}' if len(k) == 1 else 'no_change_both'
        per_hap = collections.defaultdict(list)
        for h in copies:
            per_hap[h['tumor_hap']].append(h)
        st = []
        for v in per_hap.values():
            st.append('truth' if any(x['T_eq_A'] == 'True' for x in v) else
                      'normal' if any(x['T_eq_hapX'] == 'True' or x['T_eq_hapY'] == 'True' for x in v) else 'novel')
        oc = ('no hit' if not st else 'tumor carries the truth allele (not in normal)' if 'truth' in st else
              'tumor carries only normal haplotype sequences' if all(s == 'normal' for s in st) else
              'one tumor hap = normal hap, other hap novel non-truth sequence' if 'normal' in st else 'novel non-truth sequence only')
        out[t] = dict(set=L['set'], tumor_verdict=ver, n_copies=len(copies), n_t_hap1=len(per_hap['t_hap1']),
                      n_t_hap2=len(per_hap['t_hap2']), r37_outcome=oc,
                      copies=';'.join(f'{a}:{b}:{"ALT" if c else "notALT"}' for a, b, c in cls))
    return out, hits


def check22(hits, V):
    ref = rd(f'{A22}/altlen_tumor_asm_22.tsv')
    summ = {r['truth_id']: r['outcome'] for r in rd(f'{A22}/altlen_tumor_asm_22_summary.tsv')}
    key = lambda h: (h['truth_id'], h['tumor_hap'], h['contig'], h['pos0'], h['strand'], h['len_T'], h['T_eq_A'],  # noqa: E731
                     h['T_eq_R'], h['T_eq_hapX'], h['T_eq_hapY'])
    a = sorted(key(h) for h in ref)
    b = sorted(key(h) for t in summ for h in hits[t])
    V22, _ = asm_verdicts([dict(truth_id=t, in_nf='True', set='absent_INDEL') for t in summ])   # independent of the E5d class
    oc = sum(V22[t]['r37_outcome'] == summ[t] for t in summ)
    # the 22 by the current E5d class (2026-11-05: a hap with the truth plus other germline records counts as carrying it)
    cat = {r['truth_id']: r['category'] for r in rd(f'{E5D}/per_allele.tsv') if r['sample'] == 'COLO829T' and r['arm'] == 'd46'}
    by = collections.Counter((cat[t], summ[t]) for t in summ)
    return dict(loci=len(summ), hit_rows_ref=len(a), hit_rows_here=len(b), hit_rows_identical=a == b, outcomes_identical=oc,
                by_e5d_class={f'{k[0]}: {k[1]}': v for k, v in sorted(by.items())})


def reads_tables(loci):
    b4 = rd(f'{DATA}/tumor_reads/b4.tsv')
    b2 = {(r['truth_id'], r['platform']): r for r in rd(f'{DATA}/tumor_reads/b2.tsv')}
    ids = {'absorbed': {L['truth_id'] for L in loci if L['set'] == 'absent_INDEL'},
           'control': {L['truth_id'] for L in loci if L['set'] == 'control_INDEL'}}
    res = {}
    for P in ('Illumina', 'Fiberseq', 'ONT'):
        for s in ('absorbed', 'control'):
            rows = [o for o in b4 if o['platform'] == P and o['set'] == s and o['truth_id'] in ids[s]]   # current loci.tsv
            ok = [o for o in rows if o['status'] == 'ok']
            strict = [o for o in ok if int(o['EXP']) not in (int(b2[(o['truth_id'], P)]['NE_len']), int(b2[(o['truth_id'], P)]['NO_len']))]
            res[f'{s}_{P}'] = dict(loci=len(rows), status=dict(collections.Counter(o['status'] for o in rows)), ok=len(ok),
                                   EXP_seen=sum(o['EXP_seen'] == 'True' for o in ok), strict_ok=len(strict),
                                   strict_EXP_seen=sum(o['EXP_seen'] == 'True' for o in strict))
    return res


def hg008t_strict():
    """The same strict variant for HG008T (E5d d46 normal-absent INDELs + the audit controls)."""
    pa = [r for r in rd(f'{E5D}/per_allele.tsv') if r['sample'] == 'HG008T' and r['arm'] == 'd46' and r['kind2'] == 'INDEL' and r['category'] == 'absent']
    ids = {r['truth_id'] for r in pa}
    b2 = {}
    for f in (f'{BIO}/b2_tumor_reads.tsv', f'{E5D}/hg008t/b2_new.tsv'):
        for r in rd(f):
            b2[(r['truth_id'], r['platform'], r['set'])] = r       # b2_new (read second) wins, as in E5d e3
    b4 = rd(f'{E5D}/hg008t/b4_all.tsv')
    res = {}
    for P in ('Illumina', 'PacBio'):
        for s in ('absorbed', 'control'):
            ok = [o for o in b4 if o['platform'] == P and o['set'] == s and o['status'] == 'ok' and o['kind2'] == 'INDEL'
                  and (s == 'control' or o['truth_id'] in ids)]
            strict = [o for o in ok if int(o['EXP']) not in (int(b2[(o['truth_id'], P, s)]['NE_len']), int(b2[(o['truth_id'], P, s)]['NO_len']))]
            res[f'{s}_{P}'] = dict(ok=len(ok), EXP_seen=sum(o['EXP_seen'] == 'True' for o in ok), strict_ok=len(strict),
                                   strict_EXP_seen=sum(o['EXP_seen'] == 'True' for o in strict))
    return res


def main():
    loci = rd(f'{DATA}/loci.tsv')
    V, hits = asm_verdicts(loci)
    S = dict(counts=dict(collections.Counter(L['set'] for L in loci)))
    for s in ('absent_INDEL', 'absent_SNV'):
        v = [x for x in V.values() if x['set'] == s]
        tv = collections.Counter(x['tumor_verdict'] for x in v)
        uniq = [x for x in v if x['n_t_hap1'] <= 1 and x['n_t_hap2'] <= 1]
        tvu = collections.Counter(x['tumor_verdict'] for x in uniq)
        S[f'asm_{s}'] = dict(n=len(v), novel_ALT=tv['novel_ALT'], novel_any=tv['novel_ALT'] + tv['novel_other'], verdicts=dict(tv),
                             r37=dict(collections.Counter(x['r37_outcome'] for x in v)),
                             unique_copies=dict(n=len(uniq), novel_ALT=tvu['novel_ALT'], novel_any=tvu['novel_ALT'] + tvu['novel_other']))
    S['check22'] = check22(hits, V)
    S['reads'] = reads_tables(loci)
    S['hg008t_reads_strict'] = hg008t_strict()
    json.dump(S, open(HERE / 'summary.json', 'w'), indent=1)
    with open(f'{DATA}/per_locus.tsv', 'w') as w:
        cols = ['truth_id', 'set', 'tumor_verdict', 'n_copies', 'n_t_hap1', 'n_t_hap2', 'r37_outcome', 'copies']
        w.write('\t'.join(cols) + '\n')
        for t, x in V.items():
            w.write('\t'.join(str(dict(x, truth_id=t)[k]) for k in cols) + '\n')
    a, sn = S['asm_absent_INDEL'], S['asm_absent_SNV']
    L = ['# E5e tables (COLO829T, HPRC v2.1 d46, normal-absent absorbed alleles; generated by e4_summary.py)', '',
         '## 1. Tumor assembly (hifiasm of the tumor\'s own HiFi + ONT reads) at the normal-absent alleles', '',
         '| alleles | n | tumor copy in neither normal hap (novel_ALT + novel_other) | tumor copy = GRCh38 + truth (novel_ALT) | verdicts |',
         '|---|---:|---|---|---|',
         f"| INDEL (E5d absent) | {a['n']} | {pct(a['novel_any'], a['n'])} | {pct(a['novel_ALT'], a['n'])} | {a['verdicts']} |",
         f"| SNV in normal_frame (49 of 52) | {sn['n']} | {pct(sn['novel_any'], sn['n'])} | {pct(sn['novel_ALT'], sn['n'])} | {sn['verdicts']} |",
         '', f"Loci with at most one copy per tumor hap: INDEL {a['unique_copies']['n']}, novel_any {pct(a['unique_copies']['novel_any'], a['unique_copies']['n'])}, "
         f"novel_ALT {pct(a['unique_copies']['novel_ALT'], a['unique_copies']['n'])}.", '',
         f"R3.7 per-tumor-hap outcome, INDEL: {a['r37']}", '', f"SNV: {sn['r37']}", '',
         f"Check against E27 v3 colo829t_altlen_tumor_asm (22 ALT-length loci, the E5d class of each in by_e5d_class): {S['check22']}", '',
         '## 2. Tumor reads at the normal-absent absorbed INDELs (audit b2 / b4 rule, as E5d for HG008T)', '',
         '| platform | absorbed: EXP seen / loci ok | controls (400 random non-absorbed truth INDELs) | strict (EXP differs from both normal lengths): absorbed | strict: controls |',
         '|---|---|---|---|---|']
    for P in ('Illumina', 'Fiberseq', 'ONT'):
        x, c = S['reads'][f'absorbed_{P}'], S['reads'][f'control_{P}']
        L.append(f"| {P} | {pct(x['EXP_seen'], x['ok'])} | {pct(c['EXP_seen'], c['ok'])} | {pct(x['strict_EXP_seen'], x['strict_ok'])} | {pct(c['strict_EXP_seen'], c['strict_ok'])} |")
    L += ['', 'Status counts (absorbed): ' + '; '.join(f"{P} {S['reads'][f'absorbed_{P}']['status']}" for P in ('Illumina', 'Fiberseq', 'ONT')), '',
          '### HG008T, same strict variant (E5d d46 rows; not recomputed reads)', '',
          '| platform | absorbed (E5d rule) | controls | strict: absorbed | strict: controls |', '|---|---|---|---|---|']
    for P in ('Illumina', 'PacBio'):
        x, c = S['hg008t_reads_strict'][f'absorbed_{P}'], S['hg008t_reads_strict'][f'control_{P}']
        L.append(f"| {P} | {pct(x['EXP_seen'], x['ok'])} | {pct(c['EXP_seen'], c['ok'])} | {pct(x['strict_EXP_seen'], x['strict_ok'])} | {pct(c['strict_EXP_seen'], c['strict_ok'])} |")
    open(HERE / 'tables.md', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
