"""Step 12: every number the README quotes, recomputed from this folder's tables (after s10 / s11) and the s0 inputs,
printed as 'label: value' so the README can be refreshed after a re-run (e.g. when more ONT read paths are merged).
Inputs: per_variant.tsv, per_node.tsv, hprc_individuals.tsv, hprc_populations.tsv, tables.md, mechanism.md (this
folder); analysis/tensor_recall_20260930/per_truth_HG008T.tsv; $D/per_variant.tsv (path_choice).
"""
import collections, csv, re, statistics
from pathlib import Path
A = Path(__file__).resolve().parent
D = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001')
csv.field_size_limit(10**9)
P3 = ('PacBio', 'ONT', 'Illumina')
K2 = ('SNV', 'INDEL')
CATS = ['HG008N_present_broad', 'HG008N_present_rare', 'HG008N_absent_HPRC_other', 'other_ambiguous']
AB, PR = 'HG008N_absent_HPRC_other', ('HG008N_present_broad', 'HG008N_present_rare')
pv = list(csv.DictReader(open(A / 'per_variant.tsv'), delimiter='\t'))
by = {r['truth_id']: r for r in pv}
pn = list(csv.DictReader(open(A / 'per_node.tsv'), delimiter='\t'))
full = {r['truth_id']: r for r in csv.DictReader(open(D / 'per_variant.tsv'), delimiter='\t')}
auto = {f'chr{i}' for i in range(1, 23)}
tr = [r for r in csv.DictReader(open(A.parent / 'tensor_recall_20260930/per_truth_HG008T.tsv'), delimiter='\t') if r['chrom'] in auto]
out = []
say = lambda k, v: out.append(f'{k}: {v}')
pct = lambda a, b: f'{100 * a / b:.1f}%' if b else ''
med = lambda xs: round(statistics.median(xs), 3) if xs else ''
fl = lambda xs: [float(x) for x in xs if x not in ('', 'NA')]

# section 1
for p in P3:
    for kd in K2:
        isk = lambda r: (r['kind'] == 'SNP') == (kd == 'SNV')
        nc = [r for r in tr if isk(r) and r[f'{p}_status'] == 'no_candidate']
        s = [r for r in pv if r[f'{p}_perfect'] == 'True' and r['kind'] == kd]
        rd = [r for r in s if r[f'{p}_read_data'] == 'True']
        say(f'S1 {p} {kd} no_candidate / perfect / share / d9 path spells truth / chr1 BED perfect',
            f"{len(nc)} / {len(s)} / {pct(len(s), len(nc))} / {sum(r['graph_match'] != 'closest' for r in s)} / "
            f"{sum(r['chrom'] == 'chr1' and r['in_bed'] == 'True' for r in s)}")
        say(f'S1 {p} {kd} perfect loci with read data / with own truth reads', f"{len(rd)} / {sum(int(r[f'{p}_n_truth_reads'] or 0) > 0 for r in rd)}")
for kd in K2:
    s = [r for r in pv if r['kind'] == kd]
    say(f'S1 union {kd} / d9 path spells truth / closest', f"{len(s)} / {sum(r['graph_match'] != 'closest' for r in s)} / {sum(r['graph_match'] == 'closest' for r in s)}")
    say(f'S1 union {kd} match exact / with_germline / closest', ' / '.join(str(sum(r['graph_match'] == m for r in s)) for m in ('exact', 'with_germline', 'closest')))
    say(f'S1 union {kd} truth_read_support yes / no_truth_spelling_read / no_data', ' / '.join(
        str(sum(r['truth_read_support'] == m for r in s)) for m in ('yes', 'no_truth_spelling_read', 'no_perfect_platform_read_data')))
say('S1 read data loci PacBio / ONT / Illumina', ' / '.join(str(sum(r[f'{p}_read_data'] == 'True' for r in pv)) for p in P3))
say('S1 all-three-platform INDEL / SNV', f"{sum(r['kind'] == 'INDEL' and all(r[f'{p}_perfect'] == 'True' for p in P3) for r in pv)} / "
    f"{sum(r['kind'] == 'SNV' and all(r[f'{p}_perfect'] == 'True' for p in P3) for r in pv)}")
cls = collections.Counter()
for r in pv:
    if r['kind'] == 'SNV':
        c = {r[f'{p}_miss_class'][:3] for p in P3 if r[f'{p}_perfect'] == 'True'}
        if not c & {'S3_', 'S4b'}:
            cls[r['category']] += 1
say('S1 SNVs with only S1/S2/S4a classes, by category', dict(cls))
say('S1 in_bed False', sum(r['in_bed'] == 'False' for r in pv))
t = open(A / 'tables.md').read()
rv = re.search(r'### Read-vs-graph.*?\n\n(.*?)\n\n', t, re.S).group(1).splitlines()[2:]
for line in rv:
    c = [x.strip() for x in line.strip('|').split('|')]
    say(f'S1/S2 read-vs-graph {c[0]} {c[1]}: perfect / read data / truth reads / germline-only majority / majority = graph',
        f'{c[2]} / {c[3]} / {c[5]} / {c[6]} / {c[8]}')
# section 2
nodes = {x for r in pn if r['type'] == 'branch' for x in re.findall(r'\d+', r['nodes'])}
say('S2 per_node rows / distinct elements / loci / branch nodes / skip edges',
    f"{len(pn)} / {len({r['element_id'] for r in pn})} / {len({r['truth_id'] for r in pn})} / {len(nodes)} / {len({(r['x'], r['y']) for r in pn if r['type'] == 'skip'})}")
npl = collections.Counter(r['truth_id'] for r in pn)
say('S2 loci with 1 / 2+ elements / max / none', f"{sum(v == 1 for v in npl.values())} / {sum(v > 1 for v in npl.values())} / {max(npl.values())} / {len(pv) - len(npl)}")
for lab, f in (('absent', lambda c: c == AB), ('present', lambda c: c in PR)):
    for kd in K2:
        loci = collections.defaultdict(set)
        for r in pn:
            if f(r['category']) and by[r['truth_id']]['kind'] == kd:
                loci[r['subtype']].add(r['truth_id'])
        say(f'S2 {lab} {kd} loci by subtype', {k: len(v) for k, v in sorted(loci.items())})
say('S2 absent elements by alone', dict(collections.Counter(r['alone'] for r in pn if r['category'] == AB)))
say('S2 path choice', dict(collections.Counter(re.sub(r':\d+/\d+$', '', r['path_choice']) for r in full.values())))
# section 3
ab = [r for r in pv if r['category'] == AB]
prs = [r for r in pv if r['category'] in PR]
ev = collections.Counter()
for r in pv:
    n = sum(r[f'HG008N_hap{h}_carries_ALT'] == 'ALT' for h in (1, 2))
    unres = 'unresolved' in (r['HG008N_hap1_carries_ALT'], r['HG008N_hap2_carries_ALT'])
    ev[(r['kind'], 'unresolved' if unres else n, r['category'])] += 1
say('S3 normal haps carrying ALT (kind, n haps, category)', dict(sorted(ev.items(), key=str)))
say('S3 present pattern', dict(collections.Counter(r['event_hap_pattern'] for r in prs)))
say('S3 present dipcall identical / tumor all_eq', f"{sum(r['dipcall_relation'] == 'identical germline allele' for r in prs)} / "
    f"{sum(r['HG008T_verdict'].startswith('all_eq') for r in prs)} of {len(prs)}")
say('S3 absent allele class (kind, class)', dict(collections.Counter((r['kind'], r['HG008N_allele_class']) for r in ab)))
say('S3 absent tumor verdict', dict(collections.Counter(r['HG008T_verdict'] for r in ab)))
say('S3 absent novel_other with T_info novel', sum(r['HG008T_verdict'] == 'novel_other' and r['T_info_class'] == 'novel' for r in ab))
w = lambda r: (r['HG008N_hap1_walks_elements'], r['HG008N_hap2_walks_elements'])
allw = [r for r in ab if 'all' in w(r)]
say('S3 absent loci where a normal hap walks all elements (n, by kind, by allele class) / walks >= 1',
    f"{len(allw)} {dict(collections.Counter(r['kind'] for r in allw))} {dict(collections.Counter(r['HG008N_allele_class'] for r in allw))} / "
    f"{sum(bool({'all', 'some'} & set(w(r))) for r in ab)}")
wn = lambda r: 'True' in (r['HG008N_hap1_giraffe_walk_contains'], r['HG008N_hap2_giraffe_walk_contains'])
say('S3 element rows walked by a normal hap: absent / present', f"{sum(wn(r) for r in pn if r['category'] == AB)} of {sum(r['category'] == AB for r in pn)} / "
    f"{sum(wn(r) for r in pn if r['category'] in PR)} of {sum(r['category'] in PR for r in pn)}")
# section 4
say('S4 median complete HPRC haplotypes', statistics.median(int(r['hprc_n_complete']) for r in pv))
for kd in K2:
    s = [r for r in ab if r['kind'] == kd]
    say(f'S4 absent {kd} median exact / exact (>=1 carrier) / element set / element min',
        f"{med(fl(r['hprc_exact_allele_freq'] for r in s))} / {med(fl(r['hprc_exact_allele_freq'] for r in s if r['hprc_support'] == 'exact_allele'))} / "
        f"{med(fl(r['hprc_element_set_freq'] for r in s))} / {med(fl(r['hprc_element_min_freq'] for r in s))}")
    say(f'S4 absent {kd} sub_class', dict(collections.Counter(r['sub_class'] for r in s)))
    say(f'S4 absent {kd} no exact carrier', sum(r['hprc_exact_allele_carriers'] in ('0', '') for r in s))
say('S4 present median exact', med(fl(r['hprc_exact_allele_freq'] for r in prs)))
say('S4 absent median exact (INDEL+SNV)', med(fl(r['hprc_exact_allele_freq'] for r in ab)))
ex = [r for r in ab if r['hprc_support'] == 'exact_allele']
say('S4 absent exact_allele loci / median individuals / single individual',
    f"{len(ex)} / {statistics.median(int(r['hprc_exact_n_individuals']) for r in ex)} / {sum(r['hprc_exact_n_individuals'] == '1' for r in ex)}")
sp = collections.Counter(sum(int(r[f'{s}_exact_carriers'] or 0) > 0 for s in ('AFR', 'AMR', 'EAS', 'SAS')) for r in ex)
say('S4 absent exact_allele loci by number of superpopulations', dict(sorted(sp.items())))
rare = [r for r in pv if r['hprc_exact_freq_class'] == '<0.2']
say('S4 rare exact carriers range / present_rare with 0', f"{min(int(r['hprc_exact_allele_carriers']) for r in rare)}-"
    f"{max(int(r['hprc_exact_allele_carriers']) for r in rare)} / {sum(r['category'] == 'HG008N_present_rare' and r['hprc_exact_allele_carriers'] == '0' for r in pv)}")
say('S4 CHM13 exact absent / present', f"{sum(r['CHM13_carries_exact'] == 'True' for r in ab)} / {sum(r['CHM13_carries_exact'] == 'True' for r in prs)}")
ind = [r for r in csv.DictReader(open(A / 'hprc_individuals.tsv'), delimiter='\t') if r['sample'] != 'CHM13']
for kd in K2:
    v = sorted((int(r[f'{AB}_{kd}_exact']), r['sample'] + '#' + r['hap'], int(r[f'{AB}_{kd}_complete'])) for r in ind)
    say(f'S4 per haplotype absent {kd} exact: median / min / max / median fraction of complete',
        f"{statistics.median(x[0] for x in v)} / {v[0][:2]} / {v[-1][:2]} / {med([a / c for a, _, c in v if c])}")
pops = list(csv.DictReader(open(A / 'hprc_populations.tsv'), delimiter='\t'))
for s in ('AFR', 'AMR', 'EAS', 'SAS', 'reference'):
    rows = [r for r in pops if r['superpopulation'] == s]
    rng = lambda col: f"{min(float(r[col]) for r in rows):.3f}-{max(float(r[col]) for r in rows):.3f}"
    say(f'S4 population {s}: codes / haplotypes / absent INDEL / absent SNV / present_broad INDEL',
        f"{'/'.join(r['population_code'] for r in rows)} / {'/'.join(r['haplotypes'] for r in rows)} / {rng(AB + '_INDEL_carrier_fraction')} / "
        f"{rng(AB + '_SNV_carrier_fraction')} / {rng('HG008N_present_broad_INDEL_carrier_fraction')}")
for m in re.finditer(r'\| (HG008N_\w+) \| exact allele \| (\d+) \| (\d+) \| (.*?) \| (.*?) \| (.*?) \| (.*?) \| (\d+) \|', t):
    say(f'S4 O/E exact {m.group(1)}', f'AFR {m.group(4)} AMR {m.group(5)} EAS {m.group(6)} SAS {m.group(7)}')
m = re.search(r'\| d9 \| INDEL \| (\d+) \| (\d+) \| (\d+) \|', t); m2 = re.search(r'\| full \| INDEL \| (\d+) \| (\d+) \| (\d+) \|', t)
say('S4 VCF identical carrier sets INDEL d9 / full', f'{m.group(3)} / {m2.group(3)}')
# section 5
for scope, f in (('chr1-22', lambda r: True), ('chr1 BED', lambda r: r['chrom'] == 'chr1' and r['in_bed'] == 'True'),
                 ('chr2-22', lambda r: r['chrom'] != 'chr1'), ('nogermline BED', lambda r: r['in_nogermline_bed'] == 'True')):
    for setn in P3 + ('union',):
        for kd in K2:
            s = [r for r in pv if f(r) and r['kind'] == kd and (setn == 'union' or r[f'{setn}_perfect'] == 'True')]
            say(f'S5 {scope} {setn} {kd}', ' / '.join(str(sum(r['category'] == c for r in s)) for c in CATS) + f' ({len(s)})')
tot = collections.Counter(r['category'] for r in pv)
say('S5 totals absent / present / ambiguous', f"{tot[AB]} ({pct(tot[AB], len(pv))}) / {tot[PR[0]] + tot[PR[1]]} ({pct(tot[PR[0]] + tot[PR[1]], len(pv))}) / "
    f"{tot['other_ambiguous']} ({pct(tot['other_ambiguous'], len(pv))})")
say('S5 ambiguous sub_reason', dict(collections.Counter(r['sub_reason'] for r in pv if r['category'] == 'other_ambiguous')))
pon = re.search(r'### PoN: %.*?\n\n.*?\n\n(.*?)\n\n', t, re.S).group(1).splitlines()[2:]
for line in pon:
    c = [x.strip() for x in line.strip('|').split('|')]
    if c[0] in (AB, 'not absorbed (other HG008T truths)'):
        say(f'S5 PoN {c[0]} {c[1]} n / repo / 1e-3 / 0.01 / 0.05', ' / '.join(c[2:]))
mech = open(A / 'mechanism.md').read()
row = lambda lab: re.search(rf'\| {re.escape(lab)} \| (.*?) \|\n', mech).group(1)
say('S5 mechanism absent INDEL HP>=7 1-unit / STR 1-unit (cols: present_broad SNV, INDEL, present_rare SNV, INDEL, absent SNV, INDEL, ambiguous SNV, INDEL)',
    f"{row('HP>=7 1-unit')} | {row('STR 1-unit')}")
say('S5 mechanism absent SNV rows: repeat substitution / CpG / non-repeat SNV / multi-unit / no event',
    ' | '.join(row(x) for x in ('repeat substitution', 'non-repeat SNV (CpG Ti)', 'non-repeat SNV', 'repeat multi-unit / other', 'no patient-frame event')))
print('\n'.join(out))
