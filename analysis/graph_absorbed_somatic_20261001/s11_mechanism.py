"""Step 11: mutational context and local LOH of the absorbed truths, per final category (s9), and the absorption rate by
repeat context over all HG008T truth alleles.

Inputs: $D/per_variant.tsv (s9 categories), $D/variant_set.tsv (absorbed = perfect bypass on >= 1 platform) and three
tables of the verification step (scripts copied to audit/ here):
  audit/biology/b1_truth_context.tsv  every HG008T truth allele (17,186): GIAB BEDs, PoN AFs, statuses per platform
  audit/biology/b5_normal_frame.tsv   the truth's event in the HG008-N v6.2 frame (truth INFO HG008Nv62SOMATICVARIANT,
                                      0-based): kind SNV / INDEL, length, repeat unit, tract length and class
                                      (HP>=7 / HP4-6 / STR2-6>=3copies / VNTR>6 / none) in the patient's own sequence
  audit/biology/b6_local_loh.tsv      Illumina allele balance at <= 8 phased heterozygous dipcall SNVs within +-20 kb:
                                      both_retained / O_lost (other normal hap lost) / E_lost (event hap lost) / ...
Patient-frame class of an event: 'HP>=7 1-unit' = 1-bp INDEL in a homopolymer of >= 7 bp; 'STR 1-unit' = one repeat unit
of a period 2-6 STR; 'repeat multi-unit' = several units; 'repeat substitution' = SNV inside HP>=4 / STR / VNTR;
'non-repeat SNV (CpG Ti)' / 'non-repeat SNV' / 'non-repeat INDEL'.
Output: mechanism.md (tables), printed.
"""
import collections, csv
from pathlib import Path
D = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001')
B = D / 'audit/biology'
A = Path(__file__).resolve().parent
csv.field_size_limit(10**9)
pv = {r['truth_id']: r for r in csv.DictReader(open(D / 'per_variant.tsv'), delimiter='\t')}
b1 = {r['truth_id']: r for r in csv.DictReader(open(B / 'b1_truth_context.tsv'), delimiter='\t')}
b5 = {r['truth_id']: r for r in csv.DictReader(open(B / 'b5_normal_frame.tsv'), delimiter='\t')}
b6 = {r['truth_id']: r for r in csv.DictReader(open(B / 'b6_local_loh.tsv'), delimiter='\t')}
AUTO = {f'chr{i}' for i in range(1, 23)}


def pclass(t):
    n = b5.get(t)
    if not n or not n['nf_kind']:
        return 'no patient-frame event'
    cls, kind = n['nf_class'], n['nf_kind']
    rep = cls not in ('none', '')
    if kind == 'INDEL':
        units = abs(int(n['nf_len'])) // max(1, int(n['nf_unit'] or 1))
        if not rep:
            return 'non-repeat INDEL'
        if cls == 'HP>=7' and units == 1:
            return 'HP>=7 1-unit'
        if cls.startswith('STR') and units == 1:
            return 'STR 1-unit'
        return 'repeat multi-unit / other' if units != 1 else f'{cls} 1-unit'
    if rep:
        return 'repeat substitution'
    return 'non-repeat SNV (CpG Ti)' if n['nf_cpg_ti'] == 'True' else 'non-repeat SNV'


CATS = ['HG008N_present_broad', 'HG008N_present_rare', 'HG008N_absent_HPRC_other', 'other_ambiguous']
out = ['# Mechanism tables (s11; categories from s9)', '']


def table(title, rows, cols, fn):
    out.extend([f'### {title}', '', '| ' + ' | '.join([''] + cols) + ' |', '|' + '---|' * (len(cols) + 1)])
    for r in rows:
        out.append('| ' + ' | '.join([r] + [str(fn(r, c)) for c in cols]) + ' |')
    out.append('')


# 1. patient-frame event class per category (GRCh38 kind = the truth's VCF kind)
cnt = collections.Counter()
for t, r in pv.items():
    cnt[(pclass(t), r['category'], r['kind2'])] += 1
classes = sorted({k[0] for k in cnt}, key=lambda c: -sum(v for k, v in cnt.items() if k[0] == c))
cols = [f'{c} {k}' for c in CATS for k in ('SNV', 'INDEL')]
table('Patient-frame event class x category (GRCh38 kind in the column; union, chr1-22)', classes, cols,
      lambda r, c: cnt[(r, c.rsplit(' ', 1)[0], c.rsplit(' ', 1)[1])])

# 2. absorption rate by patient-frame context over all chr1-22 truth alleles
absorbed = set(pv)
rate = collections.Counter()
for t, r in b1.items():
    if r['chrom'] not in AUTO:
        continue
    n = b5.get(t, {})
    c = pclass(t)
    if c == 'HP>=7 1-unit':
        L = int(n['nf_tract'] or 0)
        c = 'HP 1-unit, tract ' + ('7-9' if L < 10 else '10-14' if L < 15 else '15-19' if L < 20 else '20-29' if L < 30 else '>=30')
    rate[(c, 'all')] += 1
    rate[(c, 'abs')] += t in absorbed
ORDER = ['HP 1-unit, tract 7-9', 'HP 1-unit, tract 10-14', 'HP 1-unit, tract 15-19', 'HP 1-unit, tract 20-29', 'HP 1-unit, tract >=30',
         'HP4-6 1-unit', 'STR 1-unit', 'repeat multi-unit / other', 'repeat substitution', 'non-repeat INDEL',
         'non-repeat SNV (CpG Ti)', 'non-repeat SNV']
rows = [r for r in ORDER if (r, 'all') in rate] + sorted({k[0] for k in rate} - set(ORDER))
table('Absorption rate by patient-frame context, all HG008T truth alleles on chr1-22', rows, ['truth alleles', 'absorbed', '%'],
      lambda r, c: rate[(r, 'all')] if c == 'truth alleles' else rate[(r, 'abs')] if c == 'absorbed'
      else f"{100 * rate[(r, 'abs')] / max(1, rate[(r, 'all')]):.1f}")

# 3. local LOH state per category
loh = collections.Counter((b6[t]['state'] if t in b6 else 'not tested', r['category']) for t, r in pv.items())
states = sorted({k[0] for k in loh})
table('Local LOH around the locus (Illumina allele balance at phased het SNVs +-20 kb) x category (union, chr1-22)',
      states, CATS, lambda r, c: loh[(r, c)])
ctl = collections.Counter(r['state'] for r in b6.values() if r['set'] == 'control')
out.append('Controls (800 random non-absorbed truths): ' + ', '.join(f'{k} {v}' for k, v in ctl.most_common()))
open(A / 'mechanism.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
