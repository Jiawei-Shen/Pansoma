"""Step 14: repeat context of the absorbed (perfect-bypass) INDEL truths, per platform and per category.

Inputs: per_variant.tsv (this folder: kind, *_perfect, category); audit/biology/b1_truth_context.tsv (GRCh38 frame, every
HG008T truth allele) and audit/biology/b5_normal_frame.tsv (patient frame) of $D (scripts in audit/).
Repeat class of an event (b1 in GRCh38, b5 in the HG008-N v6.2 event haplotype = patient frame, the event of the truth INFO
HG008Nv62SOMATICVARIANT):
  unit u   = minimal period of the inserted / deleted bases (after the shared prefix)
  tract L  = longest exact period-u stretch touching the event in the reference sequence of that frame
  HP>=7    = u 1, L >= 7 (homopolymer of >= 7 bp)          HP4-6 = u 1, L 4-6
  STR      = u 2-6, L >= 3u (>= 3 copies)                   VNTR>6 = u > 6, the unit repeated next to the event
  in_STR   = (patient frame only) the event's unit does not repeat, but it sits inside a period 1-6 stretch >= 10 bp
  none     = none of these (not in a repeat)
  SNV / complex = the patient-frame event is a substitution (class from the stretch around the base) or REF and ALT both
                  change; no event = the truth has no HG008Nv62SOMATICVARIANT
Units changed = |length change| / u (1-unit = one base of a homopolymer / one copy of an STR unit).
Outputs: indel_repeat_context.tsv (one row per INDEL truth of the set) and indel_repeat_context.md (tables), printed.
"""
import collections, csv
from pathlib import Path
A = Path(__file__).resolve().parent
B = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/audit/biology')
csv.field_size_limit(10**9)
pv = [r for r in csv.DictReader(open(A / 'per_variant.tsv'), delimiter='\t') if r['kind'] == 'INDEL']
b1 = {r['truth_id']: r for r in csv.DictReader(open(B / 'b1_truth_context.tsv'), delimiter='\t')}
b5 = {r['truth_id']: r for r in csv.DictReader(open(B / 'b5_normal_frame.tsv'), delimiter='\t')}
AUTO = {f'chr{i}' for i in range(1, 23)}
absorbed = {r['truth_id'] for r in pv}


def pclass(t):
    n = b5.get(t, {})
    if not n.get('nf_kind'):
        return 'no event'
    c = n['nf_class']
    if n['nf_kind'] != 'INDEL':
        return f"{n['nf_kind']} in {c}" if c not in ('none', 'complex', '') else n['nf_kind']
    return 'in_STR' if c.startswith('in_STR') else c


def units(t):
    n = b5.get(t, {})
    if n.get('nf_kind') != 'INDEL' or n['nf_class'] not in ('HP>=7', 'HP4-6', 'STR2-6>=3copies', 'VNTR>6'):
        return ''
    u = abs(int(n['nf_len'])) // max(1, int(n['nf_unit']))
    return '1 unit' if u == 1 else '2-3 units' if u <= 3 else '>3 units'


def tract_bin(t):
    n = b5.get(t, {})
    if n.get('nf_kind') != 'INDEL' or n['nf_class'] != 'HP>=7':
        return ''
    L = int(n['nf_tract'])
    return '7-9' if L < 10 else '10-14' if L < 15 else '15-19' if L < 20 else '20-29' if L < 30 else '>=30'


rows = []
for r in pv:
    t = r['truth_id']; n = b5.get(t, {}); g = b1[t]
    rows.append(dict(truth_id=t, chrom=r['chrom'], vcf_pos=r['vcf_pos'], vcf_ref=r['vcf_ref'], vcf_alt=r['vcf_alt'],
                     PacBio_perfect=r['PacBio_perfect'], ONT_perfect=r['ONT_perfect'], Illumina_perfect=r['Illumina_perfect'],
                     category=r['category'], grch38_class=g['ctx_class'], grch38_unit=g['indel_unit'], grch38_tract=g['ctx_tract'],
                     patient_event=n.get('nf_event', ''), patient_kind=n.get('nf_kind', ''), patient_len_change=n.get('nf_len', ''),
                     patient_unit=n.get('nf_unit', ''), patient_tract=n.get('nf_tract', ''), patient_class=pclass(t),
                     units_changed=units(t), hp_tract_bin=tract_bin(t)))
cols = list(rows[0])
with open(A / 'indel_repeat_context.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for x in rows:
        w.write('\t'.join(str(x[c]) for c in cols) + '\n')

SETS = [('PacBio', lambda x: x['PacBio_perfect'] == 'True'), ('ONT', lambda x: x['ONT_perfect'] == 'True'),
        ('Illumina', lambda x: x['Illumina_perfect'] == 'True'), ('union', lambda x: True)]
pct = lambda a, b: f'{a:,} ({100 * a / b:.1f}%)' if b else '0'
out = ['# Repeat context of the absorbed INDEL truths (s14)', '',
       'Patient frame = the truth INFO HG008Nv62SOMATICVARIANT event in the HG008-N v6.2 event haplotype; GRCh38 frame = the truth '
       'VCF record in GRCh38. Classes in the s14 docstring.', '']


def table(title, key, order):
    out.extend([f'### {title}', '', '| class | ' + ' | '.join(s for s, _ in SETS) + ' |', '|---|' + '---:|' * len(SETS)])
    sel = {s: [x for x in rows if f(x)] for s, f in SETS}
    cnt = {s: collections.Counter(key(x) for x in sel[s]) for s, _ in SETS}
    keys = [k for k in order if any(cnt[s][k] for s, _ in SETS)] + sorted({k for s, _ in SETS for k in cnt[s]} - set(order))
    for k in keys:
        out.append(f'| {k} | ' + ' | '.join(pct(cnt[s][k], len(sel[s])) for s, _ in SETS) + ' |')
    out.append('| total | ' + ' | '.join(f'{len(sel[s]):,}' for s, _ in SETS) + ' |')
    out.append('')


ORDER = ['HP>=7', 'STR2-6>=3copies', 'HP4-6', 'VNTR>6', 'in_STR', 'SNV in HP>=7', 'SNV in STR2-6>=3copies', 'SNV in HP4-6', 'SNV',
         'MNV', 'complex', 'none', 'no event']
table('Patient frame: repeat class of the INDEL truths, per perfect set', lambda x: x['patient_class'], ORDER)
table('Patient frame: units changed (repeat INDELs only)', lambda x: x['units_changed'] or '(not a repeat INDEL)',
      ['1 unit', '2-3 units', '>3 units', '(not a repeat INDEL)'])
table('Patient frame: homopolymer length (HP>=7 INDELs only)', lambda x: x['hp_tract_bin'] or '(not HP>=7 INDEL)',
      ['7-9', '10-14', '15-19', '20-29', '>=30', '(not HP>=7 INDEL)'])
table('GRCh38 frame: repeat class of the INDEL truths, per perfect set', lambda x: x['grch38_class'],
      ['HP>=7', 'STR2-6>=3copies', 'HP4-6', 'VNTR>6', 'complex', 'none'])
# GRCh38 'none' explained by the patient frame
gn = collections.Counter(x['patient_class'] for x in rows if x['grch38_class'] == 'none')
out += ['GRCh38-frame "none" INDELs (union) in the patient frame: ' + ', '.join(f'{k} {v}' for k, v in gn.most_common()), '']
# per category (union)
CATS = ['HG008N_present_broad', 'HG008N_present_rare', 'HG008N_absent_HPRC_other', 'other_ambiguous']
cc = collections.Counter((x['category'], x['patient_class']) for x in rows)
ks = [k for k in ORDER if any(cc[(c, k)] for c in CATS)]
out += ['### Patient frame class x category (union INDELs)', '', '| class | ' + ' | '.join(CATS) + ' |', '|---|' + '---:|' * len(CATS)]
for k in ks:
    out.append(f'| {k} | ' + ' | '.join(str(cc[(c, k)]) for c in CATS) + ' |')
out.append('')
# background: all chr1-22 truth INDELs (GRCh38 frame) and absorption rate by class
allI = [t for t, g in b1.items() if g['kind2'] == 'INDEL' and g['chrom'] in AUTO]
bg = collections.Counter(b1[t]['ctx_class'] for t in allI)
ab = collections.Counter(b1[t]['ctx_class'] for t in allI if t in absorbed)
out += ['### All 8,496 chr1-22 truth INDELs (GRCh38 frame) and the share absorbed (union)', '',
        '| class | truth INDELs | absorbed | % absorbed |', '|---|---:|---:|---:|']
for k in ['HP>=7', 'STR2-6>=3copies', 'HP4-6', 'VNTR>6', 'complex', 'none']:
    if bg[k]:
        out.append(f'| {k} | {bg[k]:,} | {ab[k]:,} | {100 * ab[k] / bg[k]:.1f}% |')
out.append(f'| total | {len(allI):,} | {sum(ab.values()):,} | {100 * sum(ab.values()) / len(allI):.1f}% |')
pf = collections.Counter(pclass(t) for t in allI); pa = collections.Counter(pclass(t) for t in allI if t in absorbed)
out += ['', 'Patient frame, same: ' + '; '.join(f'{k} {pa[k]:,}/{pf[k]:,} ({100 * pa[k] / pf[k]:.1f}%)' for k in ORDER if pf[k])]
open(A / 'indel_repeat_context.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
