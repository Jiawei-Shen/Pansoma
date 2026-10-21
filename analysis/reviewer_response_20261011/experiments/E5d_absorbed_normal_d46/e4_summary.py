"""E5d step 4: per absorbed allele, does the matched normal carry it, what does the tumor carry, and the d9 / d46 tables.

Sets (e1, sets.tsv): absorbed on the arm = no candidate on >= 1 platform, a tensor on the GRCh38-only arm, the ALT spelled
by a graph path (E5b; the set of R1.4). Each allele is classified with the 2026-10-01 rules that do not depend on reads or
HPRC walks (every set allele is exactly spelled, so the 'closest' / 'with_germline' / element rules of s9 do not apply):
HG008T (analysis/graph_absorbed_somatic_20261001 s9 section 3-4):
  hap call = s8b array call (assembly; dipcall where the assembly is NA; assembly and dipcall disagreeing on ALT ->
  unresolved); both NA -> the s6 narrow minimap2 call (s9 carries() without the giraffe walk, which is used there only to
  confirm a 'no'). cars = haps calling ALT; pattern by the truth INFO event hap; T_info_class eq_other_hap without a
  carrier = patient-frame present. Category: present (other_hap_only / no_event_hap / patient frame; not when the INFO
  event gives a novel tumor allele), ambiguous (ALT on both haps, on the event hap only, INFO-novel conflict, a hap
  unresolved), absent (both haps resolved, neither carries the ALT).
  Tumor: s8b tumor_verdict (novel_ALT = a tumor contig equals GRCh38 + truth over the array, novel_other = another
  sequence in neither normal hap). Tumor reads: e3 b4_all.tsv (Illumina, the audit's EXP rule) against the audit's
  random non-absorbed truth INDEL controls.
COLO829T (analysis/graph_absorbed_somatic_colo829t_20261005 c3 / c7): hap allele from normal_frame.tsv (every truth
  INDEL and the 67 SNVs of 2026-10-05); SNVs outside it: COLO829BL dipcall at the site (inside dip.bed, PASS, phased GT
  hapY|hapX). present = a hap allele ALT (both -> ambiguous); unresolved = a hap NA and the other not ALT; absent = both
  resolved, neither ALT; germline_like = a hap has exactly the ALT length with other bases (alt_len_haps) and no hap ALT.
On nodes of other HPRC individuals: holds by the AF filter for every set allele (d9 / d46 keep only nodes and edges of
>= 10 % of the 90 / 464 haplotypes; HG008 and COLO829 are in neither graph); alt_min_paths = the fewest haplotypes on a
non-GRCh38 node of the ALT path (E5b; empty when the path uses only skip edges).
Outputs: summary.json and tables.md in this folder; per_allele.tsv in $DATA.
Run: python e4_summary.py   (login node, a minute)
"""
import collections, csv, gzip, itertools, json, statistics
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46'
OLD = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
COLO = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
CDIP = '/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip'
E5B = '/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46'
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))  # noqa: E731


def apply(region, base, edits):
    s, cur = [], 0
    for p0, a, b in sorted(edits):
        off = p0 - base
        if off < cur or off + len(a) > len(region):
            return None
        s += [region[cur:off], b]
        cur = off + len(a)
    return ''.join(s) + region[cur:]


def diffs(L, s):
    b, e = int(L['win_start0']), []
    for d in filter(None, (s or {}).get('hap_diffs', '').split(';')):
        p, v = d.split(':', 1)
        p0 = int(p) - 1
        x = tuple(v[1:].split('>')) if v[0] == 'X' else (v[1:], '') if v[0] == 'D' else ('', v[1:])
        p0 += v[0] == 'I'
        if L['ref_hap'][p0 - b:p0 - b + len(x[0])] == x[0] and b <= p0 <= b + len(L['ref_hap']):
            e.append((p0, *x))
    return e


def truth_in_diffs(L, s):
    b, e = int(L['win_start0']), diffs(L, s)
    return any(apply(L['ref_hap'], b, list(c)) == L['alt_hap'] for k in (1, 2) for c in itertools.combinations(e, k))


def narrow(s, L):
    """s9 carries() with no giraffe walk (gi = '')."""
    mm = (s or {}).get('status', 'none')
    if s and s.get('spans_agree') == 'False':
        return 'unresolved'
    if mm == 'ALT':
        ok = apply(L['ref_hap'], int(L['win_start0']), diffs(L, s)) == L['alt_hap'] or truth_in_diffs(L, s)
        return 'ALT' if ok else 'unresolved'
    if mm == 'other' and truth_in_diffs(L, s):
        return 'ALT_plus_other'
    if mm in ('REF', 'other'):
        return 'no'
    return 'unresolved'


def hap_call(ar, h, nr):
    x, y = ar[f'{h}_asm_call'], ar[f'{h}_dip_call']
    isalt = lambda c: c in ('ALT', 'ALT_plus_flank')  # noqa: E731
    if x != 'NA' and y != 'NA' and isalt(x) != isalt(y):
        return 'unresolved'
    c = x if x != 'NA' else y
    if c != 'NA':
        return 'ALT' if isalt(c) else 'no'
    return 'no' if nr == 'no' else 'unresolved'


def hg008t_alleles():
    """{truth_id: dict(category, pattern, tumor_verdict, ...)} for the old 2,478 and the E5d new loci."""
    out = {}
    for d in (OLD, f'{DATA}/hg008t'):
        loci = {r['truth_id']: r for r in rd(f'{d}/loci.tsv')}
        arr = {r['truth_id']: r for r in rd(f'{d}/array_alleles.tsv')}
        st = {}
        for r in rd(f'{d}/asm_seq_status.tsv'):
            if r['assembly'] == 'normal':
                st[(r['truth_id'], r['contig'].rsplit('_', 1)[-1])] = r
        for t, L in loci.items():
            ar = arr[t]
            ev = L['asm_event']
            evh = ev.split(':')[0].split('_')[-1] if ev else ''
            calls = {h: hap_call(ar, h, narrow(st.get((t, h)), L)) for h in ('hap1', 'hap2')}
            cars = [h for h in ('hap1', 'hap2') if calls[h] == 'ALT']
            unres = [h for h in ('hap1', 'hap2') if calls[h] == 'unresolved']
            tic = ar['T_info_class']
            if cars:
                pat = 'no_event_hap' if not evh else 'both' if len(cars) == 2 else 'event_hap_only' if cars == [evh] else 'other_hap_only'
            elif tic == 'eq_other_hap':
                pat, cars = 'other_hap_only:patient_frame', ['x']
            else:
                pat = 'unresolved' if unres else 'neither'
            if pat in ('both', 'event_hap_only'):
                cat = 'ambiguous'
            elif cars:
                cat = 'ambiguous' if tic == 'novel' else 'present'
            elif pat == 'unresolved':
                cat = 'ambiguous'
            else:
                cat = 'absent'
            out[t] = dict(category=cat, pattern=pat, tumor_verdict=ar['tumor_verdict'], hap1=calls['hap1'], hap2=calls['hap2'])
    return out


def colo_alleles(sets):
    nf = {r['truth_id']: r for r in rd(f'{COLO}/normal_frame.tsv')}
    dip = pysam.VariantFile(f'{CDIP}.vcf.gz')
    bed = collections.defaultdict(list)
    for line in open(f'{CDIP}.bed'):
        c, a, b = line.split()[:3]
        bed[c].append((int(a), int(b)))
    out = {}
    for r in sets:
        if r['sample'] != 'COLO829T' or r['truth_id'] in out:
            continue
        t = r['truth_id']
        if t in nf:
            al = {h: nf[t][f'{h}_allele'] for h in ('hapX', 'hapY')}
            alt_len = nf[t]['alt_len_haps']
            src = 'normal_frame'
        else:                                                        # SNVs outside normal_frame: dipcall at the site
            assert r['kind2'] == 'SNV', t
            c, p0 = r['chrom'], int(r['vcf_pos']) - 1
            inb = any(a <= p0 < b for a, b in bed[c])
            al = {h: 'NA' for h in ('hapX', 'hapY')}
            if inb:
                al = {'hapX': 'REF', 'hapY': 'REF'}
                for rec in dip.fetch(c, p0, p0 + 1):
                    if rec.filter.keys() and 'PASS' not in rec.filter.keys():
                        al = {h: 'NA' for h in al}
                        continue
                    gt = rec.samples[0]['GT']
                    for h, i in (('hapY', 0), ('hapX', 1)):
                        g = gt[i] if len(gt) > i else None
                        if g is None:
                            al[h] = 'NA'
                        elif g > 0 and rec.start == p0 and rec.ref == r['vcf_ref'] and rec.alts[g - 1] == r['vcf_alt']:
                            al[h] = 'ALT'
                        elif g > 0:
                            al[h] = 'germline_other'
            alt_len, src = '', 'dipcall_site'
        cars = [h for h, a in al.items() if a == 'ALT']
        unres = [h for h, a in al.items() if a == 'NA' and h not in cars]
        cat = 'ambiguous' if len(cars) == 2 else 'present' if cars else 'ambiguous' if unres else 'absent'
        out[t] = dict(category=cat, germline_like=bool(alt_len) and not cars, source=src, hapX=al['hapX'], hapY=al['hapY'])
    return out


def pct(a, b):
    return f'{a:,} ({100 * a / b:.1f}%)' if b else '-'


def main():
    sets = rd(f'{DATA}/sets.tsv')
    H, C = hg008t_alleles(), colo_alleles(sets)
    b4 = rd(f'{DATA}/hg008t/b4_all.tsv')
    S, rows = {}, []
    for (sample, arm) in (('HG008T', 'd9'), ('HG008T', 'd46'), ('COLO829T', 'd9'), ('COLO829T', 'd46')):
        x = [r for r in sets if r['sample'] == sample and r['arm'] == arm]
        A = H if sample == 'HG008T' else C
        for r in x:
            rows.append(dict(r, **{k: v for k, v in A[r['truth_id']].items()}))
        s = dict(n=len(x), SNV=sum(r['kind2'] == 'SNV' for r in x), INDEL=sum(r['kind2'] == 'INDEL' for r in x),
                 pon_tagged=sum(r['pon'] not in ('', 'none') for r in x))
        for k in ('SNV', 'INDEL', 'all'):
            y = [r for r in x if k == 'all' or r['kind2'] == k]
            cc = collections.Counter(A[r['truth_id']]['category'] for r in y)
            s[f'{k}_cat'] = dict(cc)
        mp = [int(r['alt_min_paths']) for r in x if r['alt_min_paths'] not in ('', 'None')]
        s['alt_min_paths'] = dict(n=len(mp), min=min(mp), median=statistics.median(mp), skip_only=len(x) - len(mp))
        if sample == 'HG008T':
            ab = [r for r in x if r['kind2'] == 'INDEL' and H[r['truth_id']]['category'] == 'absent']
            tv = collections.Counter(H[r['truth_id']]['tumor_verdict'] for r in ab)
            s['absent_INDEL_tumor'] = dict(n=len(ab), novel_ALT=tv['novel_ALT'], novel_any=tv['novel_ALT'] + tv['novel_other'], verdicts=dict(tv))
            ids = {r['truth_id'] for r in ab}
            for P in ('Illumina', 'PacBio'):
                ok = [o for o in b4 if o['platform'] == P and o['set'] == 'absorbed' and o['truth_id'] in ids and o['status'] == 'ok']
                s[f'absent_INDEL_reads_{P}'] = dict(ok=len(ok), EXP_seen=sum(o['EXP_seen'] == 'True' for o in ok))
        else:
            ind = [r for r in x if r['kind2'] == 'INDEL']
            pres = sum(C[r['truth_id']]['category'] == 'present' for r in ind)
            gl = sum(C[r['truth_id']]['germline_like'] for r in ind)
            s['INDEL_present'], s['INDEL_germline_like'] = pres, gl
            s['INDEL_up_to_germline_like'] = pres + gl
        S[f'{sample}_{arm}'] = s
    for P in ('Illumina', 'PacBio'):
        ok = [o for o in b4 if o['platform'] == P and o['set'] == 'control' and o['kind2'] == 'INDEL' and o['status'] == 'ok']
        S[f'control_INDEL_reads_{P}'] = dict(ok=len(ok), EXP_seen=sum(o['EXP_seen'] == 'True' for o in ok))
    # check: these rules on the 2026-10-01 / 10-05 read-level sets against their published categories
    pv = {r['truth_id']: r for r in rd(f'{OLD}/per_variant.tsv')}
    m = collections.Counter((pv[t]['category'].replace('HG008N_', '').replace('_broad', '').replace('_rare', '').replace('_HPRC_other', '').replace('other_', ''),
                             H[t]['category']) for t in pv)
    S['check_HG008T_old_vs_rules'] = {f'{a}->{b}': n for (a, b), n in sorted(m.items())}
    cpv = {r['truth_id']: r for r in rd(f'{COLO}/hprc_per_variant.tsv')}
    Cold = colo_alleles([dict(r, sample='COLO829T') for r in rd(f'{COLO}/variant_set.tsv')])
    m = collections.Counter((cpv[t]['category'].replace('normal_', '').replace('_broad', '').replace('_rare', '').replace('_HPRC_other', '').replace('other_', ''),
                             Cold[t]['category']) for t in cpv)
    S['check_COLO829T_old_vs_rules'] = {f'{a}->{b}': n for (a, b), n in sorted(m.items())}
    cols = list(rows[0])
    with open(f'{DATA}/per_allele.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for r in rows:
            w.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    (HERE / 'summary.json').write_text(json.dumps(S, indent=1) + '\n')
    tables(S)
    print(json.dumps(S, indent=1))


def tables(S):
    L = ['# E5d tables (generated by e4_summary.py)', '',
         '## 1. Absorbed truth alleles and the matched normal (set = no candidate on >= 1 platform, tensor on the GRCh38-only arm, ALT spelled by the graph)', '',
         '| sample | graph | absorbed SNV / INDEL | PoN-tagged | normal carries ALT (present) | absent from both normal haps | ambiguous |',
         '|---|---|---|---|---|---|---|']
    for k, nm in (('HG008T_d9', 'HG008T | HPRC v1.1 (d9)'), ('HG008T_d46', 'HG008T | HPRC v2.1 (d46)'),
                  ('COLO829T_d9', 'COLO829T | HPRC v1.1 (d9)'), ('COLO829T_d46', 'COLO829T | HPRC v2.1 (d46)')):
        s = S[k]
        c = s['all_cat']
        L.append(f"| {nm} | {s['SNV']:,} / {s['INDEL']:,} | {pct(s['pon_tagged'], s['n'])} | {pct(c.get('present', 0), s['n'])} | "
                 f"{pct(c.get('absent', 0), s['n'])} | {pct(c.get('ambiguous', 0), s['n'])} |")
    L += ['', 'Per kind (present / absent / ambiguous):', '']
    for k in ('HG008T_d9', 'HG008T_d46', 'COLO829T_d9', 'COLO829T_d46'):
        s = S[k]
        L.append(f"- {k}: SNV {s['SNV_cat']}; INDEL {s['INDEL_cat']}; ALT-path non-GRCh38 nodes: min {s['alt_min_paths']['min']}, "
                 f"median {s['alt_min_paths']['median']} haplotypes (n {s['alt_min_paths']['n']}; skip edges only {s['alt_min_paths']['skip_only']})")
    L += ['', '## 2. HG008T: tumor assembly and tumor reads at the absorbed INDELs absent from the normal', '',
          '| graph | absent INDELs | tumor contig in neither normal hap (novel_ALT + novel_other) | tumor contig = GRCh38 + truth (novel_ALT) | Illumina reads show the expected allele | PacBio reads |',
          '|---|---|---|---|---|---|']
    for k, nm in (('HG008T_d9', 'HPRC v1.1 (d9)'), ('HG008T_d46', 'HPRC v2.1 (d46)')):
        s = S[k]
        t = s['absent_INDEL_tumor']
        i, p = s['absent_INDEL_reads_Illumina'], s['absent_INDEL_reads_PacBio']
        L.append(f"| {nm} | {t['n']:,} | {pct(t['novel_any'], t['n'])} | {pct(t['novel_ALT'], t['n'])} | {pct(i['EXP_seen'], i['ok'])} | {pct(p['EXP_seen'], p['ok'])} |")
    ci, cp = S['control_INDEL_reads_Illumina'], S['control_INDEL_reads_PacBio']
    L.append(f"| random non-absorbed truth INDELs (audit controls) | - | - | - | {pct(ci['EXP_seen'], ci['ok'])} | {pct(cp['EXP_seen'], cp['ok'])} |")
    L += ['', '## 3. COLO829T: germline-like absorbed INDELs', '', '| graph | absorbed INDELs | normal carries ALT | + a normal hap has the ALT length, other bases | up to (germline-like) |', '|---|---|---|---|---|']
    for k, nm in (('COLO829T_d9', 'HPRC v1.1 (d9)'), ('COLO829T_d46', 'HPRC v2.1 (d46)')):
        s = S[k]
        L.append(f"| {nm} | {s['INDEL']:,} | {pct(s['INDEL_present'], s['INDEL'])} | {s['INDEL_germline_like']} | {pct(s['INDEL_up_to_germline_like'], s['INDEL'])} |")
    L += ['', '## 4. Check: these rules on the read-level d9 sets of 2026-10-01 / 10-05 (published category -> rule here)', '',
          f"- HG008T (2,478): {S['check_HG008T_old_vs_rules']}", f"- COLO829T (497): {S['check_COLO829T_old_vs_rules']}"]
    (HERE / 'tables.md').write_text('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
