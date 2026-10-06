"""Step 0: the HG008T somatic truth alleles that got no candidate because the ALT reads align perfectly to a graph branch.

Per platform (PacBio HiFi, ONT-UL, Illumina; current statuses from analysis/tensor_recall_20260930/per_truth_HG008T.tsv):
status no_candidate AND read-level reason (analysis/hg008_somatic_miss_20260926/<platform>/{snv,indel}_no_candidate.tsv,
column `reason`, majority of the ALT reads) is
  site_bypassed:branch:no_edit     - ALT reads leave GRCh38 through non-GRCh38 node(s) and carry no edit at the site +-5 bp
  site_bypassed:skip_edge:no_edit  - ALT reads take an edge that skips GRCh38 bases (or a repeat loop), no edit
INDEL: this is exactly class I1. SNV: S3 / S4b mostly (the class column is kept).
Output: variant_set.tsv, one row per truth allele in the union of the three platforms, with per-platform status / class /
reason and flags.
"""
import csv, collections, os
A = '/scratch/jshen/Github/Pansoma/analysis'
OUT = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
PLAT = {'PacBio': 'pacbio_v6', 'ONT': 'ont', 'Illumina': 'illumina'}
PERFECT = {'site_bypassed:branch:no_edit', 'site_bypassed:skip_edge:no_edit'}
key = lambda r: (r['chrom'], str(r['vcf_pos']), r['vcf_ref'], r['vcf_alt'])
miss = {}
for p, d in PLAT.items():
    for f in ('snv_no_candidate.tsv', 'indel_no_candidate.tsv'):
        for r in csv.DictReader(open(f'{A}/hg008_somatic_miss_20260926/{d}/{f}'), delimiter='\t'):
            miss[(p, key(r))] = r
truth = list(csv.DictReader(open(f'{A}/tensor_recall_20260930/per_truth_HG008T.tsv'), delimiter='\t'))
rows, cnt = [], collections.Counter()
for t in truth:
    k = key(t); hit = {}
    for p in PLAT:
        m = miss.get((p, k))
        st = t[f'{p}_status']
        if st == 'no_candidate' and m is None:
            cnt[('no_candidate without read-level row', p)] += 1
        hit[p] = st == 'no_candidate' and m is not None and m['reason'] in PERFECT
        t[f'{p}_reason'] = m['reason'] if m and st == 'no_candidate' else ''
        t[f'{p}_alt_reads'] = (int(m['alt_exact']) + int(m['alt_like'])) if m and st == 'no_candidate' else ''
        t[f'{p}_spanning'] = m['spanning'] if m and st == 'no_candidate' else ''
        t[f'{p}_perfect'] = hit[p]
    if any(hit.values()):
        t['kind2'] = 'SNV' if t['kind'] == 'SNP' or (len(t['vcf_ref']) == 1 and len(t['vcf_alt']) == 1) else 'INDEL'
        t['n_platforms'] = sum(hit.values())
        rows.append(t)
cols = ['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'kind2', 'passed', 'in_bed', 'n_platforms']
for p in PLAT:
    cols += [f'{p}_perfect', f'{p}_status', f'{p}_class', f'{p}_reason', f'{p}_alt_reads', f'{p}_spanning']
with open(f'{OUT}/variant_set.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for t in rows:
        w.write('\t'.join(str(t.get(c, '')) for c in cols) + '\n')
print('union', len(rows), collections.Counter(t['kind2'] for t in rows))
for p in PLAT:
    for kd in ('SNV', 'INDEL'):
        s = [t for t in rows if t[f'{p}_perfect'] and t['kind2'] == kd]
        print(p, kd, len(s), dict(collections.Counter(t[f'{p}_class'] for t in s)))
print(dict(cnt))
print('n_platforms', collections.Counter((t['kind2'], t['n_platforms']) for t in rows))
print('chr1 BED', collections.Counter(t['kind2'] for t in rows if t['chrom'] == 'chr1' and t['in_bed'] == 'True'))
