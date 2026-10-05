"""COLO829T step 0: the somatic truth alleles that got no candidate because the ALT reads align perfectly to a graph branch.

Same definition as analysis/graph_absorbed_somatic_20261001/s0_variant_set.py (HG008T), for COLO829T fiberseq / ONT /
Illumina: status no_candidate (analysis/tensor_recall_20260930/per_truth_COLO829T.tsv; tensor sets of 2026-09-30, now in
pansoma_v2_tensors/backup_ch6_linear100_20261001) AND the read-level reason of the COLO829T miss analysis
(analysis/tensor_recall_20260930/colo829t_miss/<platform>/{snv,indel}_no_candidate.tsv, column `reason`, majority of the
ALT reads) is site_bypassed:branch:no_edit or site_bypassed:skip_edge:no_edit.
Output: $D/variant_set.tsv (union of the three platforms; per platform status / class / reason / perfect flag; truth VAF_Ill,
VAF_PB and SMaHT region RGN), printed counts.
"""
import csv, collections
A = '/scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930'
OUT = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
PLAT = ('fiberseq', 'ONT', 'Illumina')
PERFECT = {'site_bypassed:branch:no_edit', 'site_bypassed:skip_edge:no_edit'}
key = lambda r: (r['chrom'], str(r['vcf_pos']), r['vcf_ref'], r['vcf_alt'])
miss = {}
for p in PLAT:
    for f in ('snv_no_candidate.tsv', 'indel_no_candidate.tsv'):
        for r in csv.DictReader(open(f'{A}/colo829t_miss/{p}/{f}'), delimiter='\t'):
            miss[(p, key(r))] = r
truth = list(csv.DictReader(open(f'{A}/per_truth_COLO829T.tsv'), delimiter='\t'))
rows, missing = [], collections.Counter()
for t in truth:
    hit = {}
    for p in PLAT:
        m = miss.get((p, key(t))); st = t[f'{p}_status']
        if st == 'no_candidate' and m is None:
            missing[p] += 1
        hit[p] = st == 'no_candidate' and m is not None and m['reason'] in PERFECT
        t[f'{p}_reason'] = m['reason'] if m and st == 'no_candidate' else ''
        t[f'{p}_alt_reads'] = (int(m['alt_exact']) + int(m['alt_like'])) if m and st == 'no_candidate' else ''
        t[f'{p}_spanning'] = m['spanning'] if m and st == 'no_candidate' else ''
        t[f'{p}_perfect'] = hit[p]
    if any(hit.values()):
        t['kind2'] = 'SNV' if t['kind'] == 'SNP' else 'INDEL'
        t['n_platforms'] = sum(hit.values())
        rows.append(t)
cols = ['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'kind2', 'passed', 'in_bed', 'VAF_Ill', 'VAF_PB', 'RGN', 'n_platforms']
for p in PLAT:
    cols += [f'{p}_perfect', f'{p}_status', f'{p}_class', f'{p}_reason', f'{p}_alt_reads', f'{p}_spanning']
with open(f'{OUT}/variant_set.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for t in rows:
        w.write('\t'.join(str(t.get(c, '')) for c in cols) + '\n')
auto = {f'chr{i}' for i in range(1, 23)}
print('no_candidate without a read-level row', dict(missing))
print('union', len(rows), collections.Counter(t['kind2'] for t in rows))
for p in PLAT:
    for kd in ('SNV', 'INDEL'):
        isk = lambda r: (r['kind'] == 'SNP') == (kd == 'SNV')
        nc = [r for r in truth if r['chrom'] in auto and isk(r) and r[f'{p}_status'] == 'no_candidate']
        s = [t for t in rows if t[f'{p}_perfect'] and t['kind2'] == kd]
        print(p, kd, 'no_candidate', len(nc), 'perfect', len(s), f'{100 * len(s) / max(1, len(nc)):.1f}%',
              dict(collections.Counter(t[f'{p}_class'] for t in s)), 'chr1', sum(t['chrom'] == 'chr1' for t in s))
print('n_platforms', dict(collections.Counter((t['kind2'], t['n_platforms']) for t in rows)))
print('RGN', dict(collections.Counter((t['kind2'], t['RGN']) for t in rows)))
