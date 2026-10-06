"""Step 13: which HPRC populations carry the graph allele of each absorbed truth (per locus and per category).

Inputs: per_variant.tsv (this folder: category, sub_class, kind); $D/hprc_membership.tsv.gz ('allele' rows: every HPRC
haplotype + CHM13 per locus, traversal_status, carries_exact_allele, carries_set_any_path); $D/hprc_sample_metadata.tsv.
Carrier of a locus = a haplotype with a complete traversal that spells the exact ALT allele over the tandem array
(level 'exact'); for the absent loci without any exact carrier (sub_class element_set_other_allele / recombinant_pieces)
also the haplotypes that walk all somatic elements of an ALT path (level 'nodes'). Frequencies = carriers / complete
haplotypes of that population at those loci (pooled). CHM13 is reported apart (reference, not a population).
Outputs: per_locus_populations.tsv (one row per truth allele: carriers / complete per superpopulation and per 1000G
population, the superpopulation pattern, private flags) and populations.md (tables).
"""
import collections, csv, gzip
from pathlib import Path
A = Path(__file__).resolve().parent
D = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001')
SP = ['AFR', 'AMR', 'EAS', 'SAS']
meta = {r['sample']: r for r in csv.DictReader(open(D / 'hprc_sample_metadata.tsv'), delimiter='\t')}
POPS = sorted({(r['superpopulation'], r['population_code']) for r in meta.values() if r['superpopulation'] in SP})
PNAME = {r['population_code']: r['population_name'] for r in meta.values()}
NHAP = collections.Counter(r['population_code'] for r in meta.values() if r['superpopulation'] in SP for _ in (1, 2))
pv = {r['truth_id']: r for r in csv.DictReader(open(A / 'per_variant.tsv'), delimiter='\t')}


def group(r):
    c = r['category']
    if c == 'HG008N_absent_HPRC_other':
        return 'absent: exact allele in HPRC' if r['sub_class'] == 'exact_allele' else 'absent: only the nodes in HPRC'
    return {'HG008N_present_broad': 'HG008-N present, broad', 'HG008N_present_rare': 'HG008-N present, rare'}.get(c, 'ambiguous')


GROUPS = ['HG008-N present, broad', 'HG008-N present, rare', 'absent: exact allele in HPRC', 'absent: only the nodes in HPRC', 'ambiguous']
L = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))   # truth -> level -> Counter
with gzip.open(D / 'hprc_membership.tsv.gz', 'rt') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['element_id'] != 'allele' or r['traversal_status'] != 'complete':
            continue
        t, s = r['truth_id'], r['sample']
        if s == 'CHM13':
            L[t]['chm13']['exact'] += r['carries_exact_allele'] == 'True'
            L[t]['chm13']['nodes'] += r['carries_set_any_path'] == 'True'
            continue
        pop, sp = meta[s]['population_code'], meta[s]['superpopulation']
        for key in (pop, sp):
            L[t]['complete'][key] += 1
            L[t]['exact'][key] += r['carries_exact_allele'] == 'True'
            L[t]['nodes'][key] += r['carries_set_any_path'] == 'True'
        if r['carries_exact_allele'] == 'True':
            L[t]['exact_ind'][s] += 1
        if r['carries_set_any_path'] == 'True':
            L[t]['nodes_ind'][s] += 1


def level(t):
    return 'nodes' if group(pv[t]) == 'absent: only the nodes in HPRC' else 'exact'


def pattern(t):
    lv = level(t)
    sps = [s for s in SP if L[t][lv][s] > 0]
    return '+'.join(sps) if sps else 'none'


rows = []
for t, r in pv.items():
    lv = level(t)
    pops_c = [p for _, p in POPS if L[t][lv][p] > 0]
    ind = L[t][lv + '_ind']
    rows.append(dict(truth_id=t, chrom=r['chrom'], vcf_pos=r['vcf_pos'], kind=r['kind'], category=r['category'], sub_class=r['sub_class'],
                     group=group(r), carrier_level=lv, sp_pattern=pattern(t), n_superpops=len([s for s in SP if L[t][lv][s] > 0]),
                     n_populations=len(pops_c), carrier_populations=','.join(pops_c), n_individuals=len(ind),
                     private=('one individual' if len(ind) == 1 else 'one population' if len(pops_c) == 1 else 'one superpopulation'
                              if pattern(t) in SP else 'none' if not ind else 'shared'),
                     CHM13=L[t]['chm13'][lv] > 0,
                     **{f'{s}_carriers': L[t][lv][s] for s in SP}, **{f'{s}_complete': L[t]['complete'][s] for s in SP},
                     **{f'{p}_carriers': L[t][lv][p] for _, p in POPS}, **{f'{p}_complete': L[t]['complete'][p] for _, p in POPS}))
cols = list(rows[0])
with open(A / 'per_locus_populations.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for x in rows:
        w.write('\t'.join(str(x[c]) for c in cols) + '\n')

out = ['# HPRC populations carrying the graph allele of the absorbed truths (s13)', '',
       'Carrier = complete HPRC haplotype spelling the exact ALT allele over the tandem array; for "absent: only the nodes in '
       'HPRC" = haplotype walking all somatic elements of an ALT path (no haplotype has the exact allele there). Panel: '
       + ', '.join(f'{p} {NHAP[p]}' for _, p in POPS) + ' haplotypes (AFR 46, AMR 32, EAS 8, SAS 2). Frequencies pooled over '
       'the loci (carriers / complete haplotypes). Loci counts with >= 1 carrier grow with the number of haplotypes of the '
       'population, so compare the frequencies.', '']
K = ('SNV', 'INDEL', 'all')
sel = lambda g, k: [x for x in rows if x['group'] == g and (k == 'all' or x['kind'] == k)]
# table 1: superpopulations
out += ['## 1. Superpopulations: loci with >= 1 carrier (share) and pooled carrier frequency', '',
        '| group | kind | loci | ' + ' | '.join(f'{s} loci' for s in SP) + ' | ' + ' | '.join(f'{s} freq' for s in SP) + ' | CHM13 loci |',
        '|' + '---|' * (4 + 2 * len(SP))]
for g in GROUPS:
    for k in K:
        x = sel(g, k)
        if not x:
            continue
        lv = lambda t: level(t)
        cells = [f"{sum(L[r['truth_id']][lv(r['truth_id'])][s] > 0 for r in x)} ({100 * sum(L[r['truth_id']][lv(r['truth_id'])][s] > 0 for r in x) / len(x):.0f}%)" for s in SP]
        fr = []
        for s in SP:
            c = sum(L[r['truth_id']][lv(r['truth_id'])][s] for r in x); n = sum(L[r['truth_id']]['complete'][s] for r in x)
            fr.append(f'{c / n:.3f}' if n else '')
        out.append(f"| {g} | {k} | {len(x)} | " + ' | '.join(cells) + ' | ' + ' | '.join(fr) + f" | {sum(r['CHM13'] for r in x)} |")
out.append('')
# table 2: populations
out += ['## 2. 1000 Genomes populations: loci with >= 1 carrier / pooled carrier frequency', '',
        '| superpop | population | haplotypes | ' + ' | '.join(f'{g} ({k})' for g in GROUPS[:4] for k in ('SNV', 'INDEL')) + ' |',
        '|' + '---|' * (3 + 8)]
for s, p in POPS:
    cells = []
    for g in GROUPS[:4]:
        for k in ('SNV', 'INDEL'):
            x = sel(g, k)
            c = sum(L[r['truth_id']][level(r['truth_id'])][p] for r in x); n = sum(L[r['truth_id']]['complete'][p] for r in x)
            cells.append(f"{sum(L[r['truth_id']][level(r['truth_id'])][p] > 0 for r in x)} / {c / n:.3f}" if n else '')
    out.append(f'| {s} | {p} ({PNAME[p]}) | {NHAP[p]} | ' + ' | '.join(cells) + ' |')
out.append('| | n loci | | ' + ' | '.join(str(len(sel(g, k))) for g in GROUPS[:4] for k in ('SNV', 'INDEL')) + ' |')
out.append('')
# table 3: superpopulation patterns
pats = collections.Counter((x['group'], x['sp_pattern']) for x in rows)
allp = sorted({p for _, p in pats}, key=lambda p: (-sum(v for (g, q), v in pats.items() if q == p), p))
out += ['## 3. Which superpopulations carry the allele (pattern per locus)', '',
        '| pattern | ' + ' | '.join(GROUPS) + ' |', '|' + '---|' * (1 + len(GROUPS))]
for p in allp:
    out.append(f'| {p} | ' + ' | '.join(str(pats[(g, p)]) for g in GROUPS) + ' |')
out.append('')
# table 4: private / shared
pr = collections.Counter((x['group'], x['private']) for x in rows)
out += ['## 4. Shared or private', '',
        '| | ' + ' | '.join(GROUPS) + ' |', '|' + '---|' * (1 + len(GROUPS))]
for p in ('shared', 'one superpopulation', 'one population', 'one individual', 'none'):
    out.append(f'| {p} | ' + ' | '.join(str(pr[(g, p)]) for g in GROUPS) + ' |')
out.append('')
out.append('shared = carriers in >= 2 superpopulations; one superpopulation = >= 2 populations of one superpopulation or ... '
           '(first match: one individual, one population, one superpopulation).')
out.append('')
# table 5: single-superpopulation loci by population
one = [x for x in rows if x['n_superpops'] == 1]
c5 = collections.Counter((x['group'], x['sp_pattern'], x['carrier_populations']) for x in one)
out += ['## 5. Loci carried by one superpopulation only: by superpopulation and populations', '',
        '| group | superpop | carrier populations | loci |', '|---|---|---|---|']
for (g, s, pp), v in sorted(c5.items(), key=lambda kv: (GROUPS.index(kv[0][0]), kv[0][1], -kv[1])):
    out.append(f'| {g} | {s} | {pp} | {v} |')
open(A / 'populations.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
