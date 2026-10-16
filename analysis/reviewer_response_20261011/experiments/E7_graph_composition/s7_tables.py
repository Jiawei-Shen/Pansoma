"""E7 step 7: markdown tables (tables.md) and summary.json from the TSVs of s3 / s4 / s6 (this folder)."""
import collections, csv, json, os

OUT = os.path.dirname(os.path.abspath(__file__))


def read(name):
    return list(csv.DictReader(open(f'{OUT}/{name}'), delimiter='\t'))


def pct(x, nd=1):
    return '' if x in ('', None) else f'{100 * float(x):.{nd}f}'


curves = read('curves.tsv')
C = {(r['sample'], r['role'], r['scope'], r['stratum'], r['pool'], int(r['k']), r['threshold_type'], r['threshold']): r
     for r in curves}
ALLK = {'all': 89, 'AFR': 46, 'nonAFR': 42}
md = []
summary = {}


def cell(s, role, scope, stratum, pool, k, tt, tv, sd=False):
    r = C.get((s, role, scope, stratum, pool, k, tt, tv))
    if r is None:
        return ''
    v = pct(r['mean'])
    if sd and int(r['n_draws']) > 1:
        v += f' ± {100 * float(r["sd"]):.1f}'
    return v


COLS = [('germline', 'SNV', '(a) germline SNV in graph'), ('germline', 'INDEL', '(a) germline INDEL'),
        ('somatic', 'SNV', '(b) somatic SNV absorbed'), ('somatic', 'INDEL', '(b) somatic INDEL'),
        ('somatic', 'INDEL_HP>=7', '(b) somatic INDEL, HP>=7'), ('somatic', 'INDEL_not_HP>=7', '(b) somatic INDEL, not HP>=7')]
n_alleles = {}
for r in curves:
    if r['pool'] == 'all' and r['k'] == '89' and r['threshold_type'] == 'abs' and r['threshold'] == '1':
        n_alleles[(r['sample'], r['role'], r['scope'], r['stratum'])] = int(r['n_alleles'])

md.append('### A. Allele counts (FILTER PASS, inside the truth BED, chr1-22; chr1 in brackets)\n')
md.append('| Sample | ' + ' | '.join(f'{role} {st}' for role, st, _ in COLS) + ' |')
md.append('|---|' + '---|' * len(COLS))
for s in ('HG008', 'COLO829', 'HCC1395'):
    md.append(f'| {s} | ' + ' | '.join(f"{n_alleles.get((s, role, 'bed', st), 0):,} ({n_alleles.get((s, role, 'bed_chr1', st), 0):,})"
                                       for role, st, _ in COLS) + ' |')

for s in ('HG008', 'COLO829', 'HCC1395'):
    md.append(f'\n### B-{s}. Haplotype count k x minimum carrier count d, all 89 HPRC v1.1 haplotypes (% of alleles in the graph; '
              'mean ± sd over 20 draws; chr1-22, PASS, in BED)\n')
    md.append('| k | d | ' + ' | '.join(c[2] for c in COLS) + ' |')
    md.append('|---|---|' + '---|' * len(COLS))
    for k in (2, 4, 8, 16, 32, 64, 89):
        for d in ('1', '2', '5', '9', '18'):
            if int(d) > k:
                continue
            md.append(f'| {k} | {d} | ' + ' | '.join(cell(s, role, 'bed', st, 'all', k, 'abs', d, sd=True) for role, st, _ in COLS) + ' |')

md.append('\n### C. Relative threshold f (d = max(1, ceil(f k))), all haplotypes, HG008 (% in the graph)\n')
md.append('| k | f | d | ' + ' | '.join(c[2] for c in COLS) + ' |')
md.append('|---|---|---|' + '---|' * len(COLS))
for k in (8, 16, 32, 64, 89):
    for f in ('0.0', '0.02', '0.05', '0.1', '0.2'):
        r = C.get(('HG008', 'germline', 'bed', 'SNV', 'all', k, 'rel', f))
        if r is None:
            continue
        md.append(f"| {k} | {'>=1 hap' if f == '0.0' else f} | {r['d']} | " +
                  ' | '.join(cell('HG008', role, 'bed', st, 'all', k, 'rel', f, sd=True) for role, st, _ in COLS) + ' |')

md.append('\n### D. Population of the drawn haplotypes: AFR (46 haplotypes) vs non-AFR (AMR + EAS + SAS, 42), equal k '
          '(% in the graph, mean over 20 draws; k = 42 for non-AFR is its whole pool)\n')
for tt, tv, label in (('abs', '1', 'd = 1'), ('rel', '0.1', 'f = 10%')):
    md.append(f'\n{label}:\n')
    md.append('| Sample | k | germline SNV AFR / non-AFR | germline INDEL AFR / non-AFR | somatic SNV AFR / non-AFR | '
              'somatic INDEL AFR / non-AFR | somatic INDEL HP>=7 AFR / non-AFR |')
    md.append('|---|---|---|---|---|---|---|')
    for s in ('HG008', 'COLO829', 'HCC1395'):
        for k in (2, 4, 8, 16, 32, 42):
            vals = []
            for role, st in (('germline', 'SNV'), ('germline', 'INDEL'), ('somatic', 'SNV'), ('somatic', 'INDEL'),
                             ('somatic', 'INDEL_HP>=7')):
                vals.append(f"{cell(s, role, 'bed', st, 'AFR', k, tt, tv)} / {cell(s, role, 'bed', st, 'nonAFR', k, tt, tv)}")
            md.append(f'| {s} | {k} | ' + ' | '.join(vals) + ' |')

gs = read('graph_size.tsv')
md.append('\n### E. Graph size: HPRC v1.1 alleles kept (chr1-22, millions; all-haplotype pool)\n')
md.append('| k | d | SNV | INDEL <= 50 bp | of which HP>=7 | SV > 50 bp |')
md.append('|---|---|---|---|---|---|')
G = {(r['region'], r['stratum'], r['pool'], int(r['k']), r['threshold_type'], r['threshold']): r for r in gs}
for k in (2, 8, 32, 89):
    for d in ('1', '2', '9', '18'):
        if int(d) > k:
            continue
        md.append(f'| {k} | {d} | ' + ' | '.join(
            f"{float(G[('chr1-22', st, 'all', k, 'abs', d)]['mean_alleles']) / 1e6:.3f}" for st in ('SNV', 'INDEL<=50', 'INDEL_HP>=7', 'SV>50')) + ' |')

cal = read('calibration.tsv')
md.append('\n### F. The real graphs on the curves: predicted (all 89 haplotypes, d = 9) vs measured candidate loss '
          '(PASS, in BED, chr1-22)\n')
md.append('Measured loss = truth alleles with a tensor on the GRCh38-only arm but none on d9 (d46), over alleles with a '
          'GRCh38-only tensor. "if in" / "if out" = d9 loss among alleles the model puts in / out of d9.\n')
md.append('| Sample | Platform | Role | Stratum | n | predicted in d9 | tensor linear / d9 / d46 | measured loss d9 | '
          'loss d46 | d9 loss if in / if out | share of d9 loss predicted |')
md.append('|---|---|---|---|---|---|---|---|---|---|---|')
for r in cal:
    if r['scope'] != 'bed' or r['stratum'] not in ('SNV', 'INDEL', 'INDEL_HP>=7', 'INDEL_not_HP>=7'):
        continue
    md.append(f"| {r['sample']} | {r['platform']} | {r['role']} | {r['stratum']} | {int(r['n']):,} | {pct(r['predicted_in_d9'])} | "
              f"{pct(r.get('tensor_linear'))} / {pct(r.get('tensor_d9'))} / {pct(r.get('tensor_d46'))} | "
              f"{pct(r.get('loss_d9_vs_linear'))} | {pct(r.get('loss_d46_vs_linear'))} | "
              f"{pct(r.get('loss_d9_if_pred_in'))} / {pct(r.get('loss_d9_if_pred_out'))} | {pct(r.get('share_of_d9_loss_pred_in'))} |")

bins = read('calibration_bins.tsv')
md.append('\n### G. Measured loss by HPRC v1.1 carrier count (alleles with a GRCh38-only tensor; PacBio / fiber-seq; '
          'PASS, in BED, chr1-22)\n')
md.append('| Sample | Role | Stratum | carriers (of 89) | n | loss d9 | loss d46 |')
md.append('|---|---|---|---|---|---|---|')
for r in bins:
    if r['platform'] not in ('PacBio', 'PacBio_fiberseq') or r['stratum'] not in ('SNV', 'INDEL'):
        continue
    md.append(f"| {r['sample']} | {r['role']} | {r['stratum']} | {r['hprc_carriers']} | {int(r['n_linear_tensor']):,} | "
              f"{pct(r['loss_d9'])} | {pct(r['loss_d46'])} |")

ret = read('d9_retention.tsv')
md.append('\n### H. How sharp is the d9 filter? HPRC v1.1 alleles found as a simple allele in the d9 deconstruct VCF, by '
          'carrier count (chr1-22)\n')
md.append('| Stratum | ' + ' | '.join(f'{a}-{b}' if a != b else str(a) for a, b in
                                      ((1, 1), (2, 4), (5, 6), (7, 7), (8, 8), (9, 9), (10, 12), (13, 20), (21, 44), (45, 89))) + ' |')
md.append('|---|' + '---|' * 10)
R = collections.defaultdict(lambda: [0, 0])
for r in ret:
    R[(r['stratum'], int(r['hprc_carriers']))][0] += int(r['hprc_alleles'])
    R[(r['stratum'], int(r['hprc_carriers']))][1] += int(r['in_d9_vcf'])
for st in ('SNV', 'INDEL', 'INDEL_HP>=7', 'INDEL_not_HP>=7', 'SV>50'):
    vals = []
    for a, b in ((1, 1), (2, 4), (5, 6), (7, 7), (8, 8), (9, 9), (10, 12), (13, 20), (21, 44), (45, 89)):
        nn = sum(R[(st, c)][0] for c in range(a, b + 1))
        dd = sum(R[(st, c)][1] for c in range(a, b + 1))
        vals.append(f'{100 * dd / nn:.1f}' if nn else '')
    md.append(f'| {st} | ' + ' | '.join(vals) + ' |')

con = read('d9_content.tsv')
md.append('\n### I. Sample alleles in HPRC v1.1: any carrier, predicted d9 (>= 9 carriers), and found as a simple allele in '
          'the d9 deconstruct VCF (PASS, in BED, chr1-22)\n')
md.append('| Sample | Role | Stratum | n | any HPRC carrier | predicted d9 | in d9 VCF | in d9 VCF if predicted / if not |')
md.append('|---|---|---|---|---|---|---|---|')
for r in con:
    if r['scope'] != 'bed' or r['stratum'] not in ('SNV', 'INDEL', 'INDEL_HP>=7', 'INDEL_not_HP>=7'):
        continue
    md.append(f"| {r['sample']} | {r['role']} | {r['stratum']} | {int(r['n']):,} | {pct(r['in_hprc_any'])} | "
              f"{pct(r['predicted_d9'])} | {pct(r['in_d9_vcf'])} | {pct(r['in_d9_vcf_given_pred'])} / {pct(r['in_d9_vcf_given_not_pred'])} |")

rel = read('relatedness.tsv')
md.append('\n### J. Relatedness check (R3.m3): rare germline SNVs of the benchmark normal (1-2 HPRC carriers) carried by '
          'each HPRC v1.1 individual; top three per sample\n')
md.append('| Sample | rare SNVs | HPRC individual | superpopulation | shared | share | ratio to superpopulation median |')
md.append('|---|---|---|---|---|---|---|')
for s in ('HG008', 'COLO829', 'HCC1395'):
    for r in [x for x in rel if x['sample'] == s][:3]:
        md.append(f"| {s} | {int(r['rare_snvs']):,} | {r['hprc_individual']} | {r['superpopulation']} | {int(r['shared']):,} | "
                  f"{pct(r['share'], 2)} | {r['ratio_to_superpop_median']} |")

open(f'{OUT}/tables.md', 'w').write('\n'.join(md) + '\n')

# summary.json: the headline numbers
for s in ('HG008', 'COLO829', 'HCC1395'):
    summary[s] = {}
    for k in (8, 16, 32, 64, 89):
        for d in ('1', '9'):
            if int(d) > k:
                continue
            summary[s][f'k{k}_d{d}'] = {f'{role}_{st}': float(C[(s, role, 'bed', st, 'all', k, 'abs', d)]['mean'])
                                         for role, st, _ in COLS if (s, role, 'bed', st, 'all', k, 'abs', d) in C}
summary['calibration'] = [{k: r[k] for k in ('sample', 'platform', 'role', 'stratum', 'n', 'predicted_in_d9',
                                              'loss_d9_vs_linear', 'loss_d46_vs_linear')}
                          for r in cal if r['scope'] == 'bed' and r['stratum'] in ('SNV', 'INDEL', 'INDEL_HP>=7')]
json.dump(summary, open(f'{OUT}/summary.json', 'w'), indent=1)
print('\n'.join(md[:40]))
