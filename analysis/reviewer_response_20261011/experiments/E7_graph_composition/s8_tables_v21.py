"""E7 step 8 (HPRC v2.1 axis): markdown tables (v2.1/tables.md) and v2.1/summary.json from the v2.1 outputs of s3 / s4 /
s6, plus the v1.1 curves (./curves.tsv) for the comparison at equal haplotype count (section K)."""
import collections, csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = f'{HERE}/v2.1'
NH = 463


def read(path):
    return list(csv.DictReader(open(path), delimiter='\t'))


def pct(x, nd=1):
    return '' if x in ('', None) else f'{100 * float(x):.{nd}f}'


def index(rows):
    return {(r['sample'], r['role'], r['scope'], r['stratum'], r['pool'], int(r['k']), r['threshold_type'],
             r['threshold']): r for r in rows}


C = index(read(f'{OUT}/curves.tsv'))
C11 = index(read(f'{HERE}/curves.tsv'))
SAMPLES = ('HG008', 'COLO829', 'HCC1395')
COLS = [('germline', 'SNV', '(a) germline SNV in graph'), ('germline', 'INDEL', '(a) germline INDEL'),
        ('somatic', 'SNV', '(b) somatic SNV absorbed'), ('somatic', 'INDEL', '(b) somatic INDEL'),
        ('somatic', 'INDEL_HP>=7', '(b) somatic INDEL, HP>=7'), ('somatic', 'INDEL_not_HP>=7', '(b) somatic INDEL, not HP>=7')]
md, summary = [], {}


def cell(c, s, role, scope, stratum, pool, k, tt, tv, sd=False):
    r = c.get((s, role, scope, stratum, pool, k, tt, tv))
    if r is None:
        return ''
    v = pct(r['mean'])
    if sd and int(r['n_draws']) > 1:
        v += f' ± {100 * float(r["sd"]):.1f}'
    return v


def val(c, *k):
    r = c.get(k)
    return None if r is None else round(100 * float(r['mean']), 2)


md.append('# E7, HPRC v2.1 axis: tables\n')
md.append('In-silico graphs from k of the 463 HPRC v2.1 haplotypes (CHM13 + 231 diploid samples; GRCh38 is the backbone). '
          'Scope: FILTER PASS, inside the truth BED, chr1-22, unless stated. Mean ± sd over 20 random draws; one draw '
          'when k is the whole pool. Built by `s8_tables_v21.py` from `curves.tsv`, `graph_size.tsv`, `calibration*.tsv`, '
          '`d46_content.tsv`, `relatedness.tsv` (this folder) and the v1.1 `../curves.tsv`.\n')

KS = (2, 4, 8, 16, 32, 64, 89, 128, 256, NH)
for s in SAMPLES:
    md.append(f'\n### A-{s}. Haplotype count k x minimum carrier count d, all 463 haplotypes (% of alleles in the graph)\n')
    md.append('| k | d | ' + ' | '.join(c[2] for c in COLS) + ' |')
    md.append('|---|---|' + '---|' * len(COLS))
    for k in KS:
        for d in ('1', '2', '5', '9', '18', '46', '92'):
            if int(d) > k or (s, 'germline', 'bed', 'SNV', 'all', k, 'abs', d) not in C:
                continue
            md.append(f'| {k} | {d} | ' + ' | '.join(cell(C, s, role, 'bed', st, 'all', k, 'abs', d, sd=True)
                                                    for role, st, _ in COLS) + ' |')

for f, label in (('0.1', 'B. Fixed relative threshold f = 10 % (d = ceil(0.1 k))'),
                 ('0.0', 'C. No frequency filter (d = 1)')):
    md.append(f'\n### {label}: the haplotype-count axis up to 463 (% in the graph)\n')
    md.append('| Sample | k | d | ' + ' | '.join(c[2] for c in COLS) + ' |')
    md.append('|---|---|---|' + '---|' * len(COLS))
    for s in SAMPLES:
        for k in KS:
            r = C.get((s, 'germline', 'bed', 'SNV', 'all', k, 'rel', f))
            if r is None:
                continue
            md.append(f"| {s} | {k} | {r['d']} | " + ' | '.join(cell(C, s, role, 'bed', st, 'all', k, 'rel', f, sd=True)
                                                              for role, st, _ in COLS) + ' |')

POPS = ('AFR', 'AMR', 'EAS', 'EUR', 'SAS')
md.append('\n### D. Superpopulation of the drawn haplotypes at equal k = 60 (the whole EUR pool; AFR 140, AMR 88, EAS 102, '
          'SAS 72 haplotypes; % in the graph, mean of 20 draws)\n')
for tt, tv, label in (('abs', '1', 'd = 1'), ('rel', '0.1', 'f = 10 % (d = 6)')):
    md.append(f'\n{label}:\n')
    md.append('| Sample | Stratum | ' + ' | '.join(POPS) + ' |')
    md.append('|---|---|' + '---|' * len(POPS))
    for s in SAMPLES:
        for role, st, lab in COLS[:5]:
            md.append(f'| {s} | {lab} | ' + ' | '.join(cell(C, s, role, 'bed', st, p, 60, tt, tv) for p in POPS) + ' |')
            summary.setdefault('superpop_k60', {}).setdefault(f'{s} {role} {st} {label}', {
                p: val(C, s, role, 'bed', st, p, 60, tt, tv) for p in POPS})

gs = read(f'{OUT}/graph_size.tsv')
G = {(r['region'], r['stratum'], r['pool'], int(r['k']), r['threshold_type'], r['threshold']): r for r in gs}
md.append('\n### E. Graph size: HPRC v2.1 alleles kept (chr1-22, millions; all-haplotype pool)\n')
md.append('| k | d | SNV | INDEL <= 50 bp | of which HP>=7 | SV > 50 bp |')
md.append('|---|---|---|---|---|---|')
for k in (8, 32, 89, 256, NH):
    for d in ('1', '2', '9', '46'):
        if int(d) > k or ('chr1-22', 'SNV', 'all', k, 'abs', d) not in G:
            continue
        md.append(f'| {k} | {d} | ' + ' | '.join(
            f"{float(G[('chr1-22', st, 'all', k, 'abs', d)]['mean_alleles']) / 1e6:.3f}"
            for st in ('SNV', 'INDEL<=50', 'INDEL_HP>=7', 'SV>50')) + ' |')
summary['graph_size_M'] = {f'k{k}_d{d}_{st}': round(float(G[('chr1-22', st, 'all', k, 'abs', d)]['mean_alleles']) / 1e6, 3)
                           for k in (89, NH) for d in ('1', '9', '46') for st in ('SNV', 'INDEL<=50')
                           if ('chr1-22', st, 'all', k, 'abs', d) in G}

cal = read(f'{OUT}/calibration.tsv')
md.append('\n### F. The real d46 graph on the v2.1 curves: predicted (>= 46 of 463 carriers) vs measured candidate loss '
          'on identical reads (PASS, in BED, chr1-22)\n')
md.append('Measured loss = truth alleles with a tensor on the GRCh38-only arm but none on d46 (d9), over alleles with a '
          'GRCh38-only tensor. "if in / if out" = d46 loss among alleles the v2.1 model puts in / out of d46.\n')
md.append('| Sample | Platform | Role | Stratum | n | predicted in d46 | tensor linear / d9 / d46 | loss d46 | loss d9 | '
          'd46 loss if in / if out | share of d46 loss predicted |')
md.append('|---|---|---|---|---|---|---|---|---|---|---|')
for r in cal:
    if r['scope'] != 'bed' or r['stratum'] not in ('SNV', 'INDEL', 'INDEL_HP>=7', 'INDEL_not_HP>=7'):
        continue
    md.append(f"| {r['sample']} | {r['platform']} | {r['role']} | {r['stratum']} | {int(r['n']):,} | "
              f"{pct(r['predicted_in_d46'])} | {pct(r.get('tensor_linear'))} / {pct(r.get('tensor_d9'))} / "
              f"{pct(r.get('tensor_d46'))} | {pct(r.get('loss_d46_vs_linear'))} | {pct(r.get('loss_d9_vs_linear'))} | "
              f"{pct(r.get('loss_d46_if_pred_in'))} / {pct(r.get('loss_d46_if_pred_out'))} | "
              f"{pct(r.get('share_of_d46_loss_pred_in'))} |")
    if r['stratum'] in ('SNV', 'INDEL') and r.get('loss_d46_if_pred_in'):
        summary.setdefault('calibration_d46', {})[f"{r['sample']} {r['platform']} {r['role']} {r['stratum']}"] = dict(
            predicted_in_d46=pct(r['predicted_in_d46']), loss_d46=pct(r['loss_d46_vs_linear']),
            loss_if_in=pct(r['loss_d46_if_pred_in']), loss_if_out=pct(r['loss_d46_if_pred_out']),
            share_predicted=pct(r['share_of_d46_loss_pred_in']))

bins = read(f'{OUT}/calibration_bins.tsv')
md.append('\n### G. Measured loss by HPRC v2.1 carrier count (of 463; alleles with a GRCh38-only tensor; PacBio / '
          'fiber-seq; PASS, in BED, chr1-22)\n')
md.append('| Sample | Role | Stratum | carriers (of 463) | n | loss d46 | loss d9 |')
md.append('|---|---|---|---|---|---|---|')
for r in bins:
    if r['platform'] not in ('PacBio', 'PacBio_fiberseq') or r['stratum'] not in ('SNV', 'INDEL'):
        continue
    md.append(f"| {r['sample']} | {r['role']} | {r['stratum']} | {r['hprc_carriers']} | {int(r['n_linear_tensor']):,} | "
              f"{pct(r['loss_d46'])} | {pct(r['loss_d9'])} |")

con = read(f'{OUT}/d46_content.tsv')
md.append('\n### H. Sample alleles in HPRC v2.1: any carrier and predicted d46 (>= 46 of 463 carriers) (PASS, in BED, '
          'chr1-22)\n')
md.append('| Sample | Role | Stratum | n | any v2.1 carrier | predicted d46 |')
md.append('|---|---|---|---|---|---|')
for r in con:
    if r['scope'] != 'bed' or r['stratum'] not in ('SNV', 'INDEL', 'INDEL_HP>=7', 'INDEL_not_HP>=7'):
        continue
    md.append(f"| {r['sample']} | {r['role']} | {r['stratum']} | {int(r['n']):,} | {pct(r['in_hprc_any'])} | "
              f"{pct(r['predicted_d46'])} |")

rel = read(f'{OUT}/relatedness.tsv')
md.append('\n### J. Relatedness (R3.m3): share of each normal\'s rare germline SNVs (1-2 of 463 v2.1 carriers; PASS, in '
          'BED) carried by each HPRC v2.1 individual; top three per sample and the highest ratio per superpopulation\n')
md.append('| Sample | rare SNVs | HPRC individual | superpopulation | shared | share % | ratio to superpopulation median |')
md.append('|---|---|---|---|---|---|---|')
for s in SAMPLES:
    rr = [r for r in rel if r['sample'] == s]
    for r in rr[:3]:
        md.append(f"| {s} | {int(r['rare_snvs']):,} | {r['hprc_individual']} | {r['superpopulation']} | "
                  f"{int(r['shared']):,} | {100 * float(r['share']):.2f} | {r['ratio_to_superpop_median']} |")
    summary.setdefault('relatedness_top', {})[s] = [(r['hprc_individual'], r['superpopulation'], float(r['share']),
                                                     r['ratio_to_superpop_median']) for r in rr[:3]]
md.append('\n| Sample | ' + ' | '.join(f'max ratio {p} (individual, share %)' for p in POPS) + ' |')
md.append('|---|' + '---|' * len(POPS))
for s in SAMPLES:
    vals = []
    for p in POPS:
        rr = [r for r in rel if r['sample'] == s and r['superpopulation'] == p and r['ratio_to_superpop_median'] != '']
        if not rr:
            vals.append('')
            continue
        b = max(rr, key=lambda r: float(r['ratio_to_superpop_median']))
        vals.append(f"{b['ratio_to_superpop_median']} ({b['hprc_individual']}, {100 * float(b['share']):.2f})")
    md.append(f'| {s} | ' + ' | '.join(vals) + ' |')
md.append('\nShare of the rare SNVs carried by the median individual of each superpopulation (%):\n')
md.append('| Sample | ' + ' | '.join(POPS) + ' |')
md.append('|---|' + '---|' * len(POPS))
for s in SAMPLES:
    vals = []
    for p in POPS:
        sh = sorted(float(r['share']) for r in rel if r['sample'] == s and r['superpopulation'] == p)
        vals.append(f'{100 * sh[len(sh) // 2]:.2f}' if sh else '')
    md.append(f'| {s} | ' + ' | '.join(vals) + ' |')

md.append('\n### K. HPRC v1.1 vs v2.1 at equal haplotype count k = 89 (v1.1: all 89 haplotypes, one graph; v2.1: 20 '
          'random subsets of 89) and the whole v2.1 set (k = 463) (% in the graph)\n')
md.append('| Sample | threshold | ' + ' | '.join(c[2] for c in COLS) + ' |')
md.append('|---|---|' + '---|' * len(COLS))
for s in SAMPLES:
    for tt, tv, lab in (('abs', '1', 'd = 1'), ('abs', '9', 'd = 9'), ('rel', '0.1', 'f = 10 %')):
        md.append(f'| {s} | {lab}: v1.1 k 89 / v2.1 k 89 / v2.1 k 463 | ' + ' | '.join(
            f"{cell(C11, s, role, 'bed', st, 'all', 89, tt, tv)} / {cell(C, s, role, 'bed', st, 'all', 89, tt, tv)} / "
            f"{cell(C, s, role, 'bed', st, 'all', NH, tt, tv)}" for role, st, _ in COLS) + ' |')
        summary.setdefault('v11_vs_v21', {})[f'{s} {lab}'] = {
            f'{role} {st}': dict(v11_k89=val(C11, s, role, 'bed', st, 'all', 89, tt, tv),
                                 v21_k89=val(C, s, role, 'bed', st, 'all', 89, tt, tv),
                                 v21_k463=val(C, s, role, 'bed', st, 'all', NH, tt, tv)) for role, st, _ in COLS}

for s in SAMPLES:
    summary.setdefault('k_axis_f10', {})[s] = {
        f'{role} {st}': {k: val(C, s, role, 'bed', st, 'all', k, 'rel', '0.1') for k in KS} for role, st, _ in COLS}
    summary.setdefault('k_axis_d1', {})[s] = {
        f'{role} {st}': {k: val(C, s, role, 'bed', st, 'all', k, 'rel', '0.0') for k in KS} for role, st, _ in COLS}
open(f'{OUT}/tables.md', 'w').write('\n'.join(md) + '\n')
json.dump(summary, open(f'{OUT}/summary.json', 'w'), indent=1)
print('wrote', f'{OUT}/tables.md', len(md), 'lines')
