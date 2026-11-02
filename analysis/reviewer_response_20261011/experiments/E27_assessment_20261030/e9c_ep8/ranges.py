"""R3.3 ranges over the twelve in-sample SNV runs (rcpu, d9 + d46), from the isolated E9c / E9i tables (epoch 8)."""
import json, re, sys
A = '/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_ep8'
T = open(f'{A}/EXP/tables.md').read().split('\n## ')
sec = {s.split('.')[0]: s for s in T[1:]}
def rows(s):
    return [[c.strip() for c in l.strip('|').split('|')] for l in s.split('\n') if l.startswith('| rcpu') or l.startswith('| ' + 'rcpu/')]
pct = lambda x: float(re.search(r'\(([\d.]+) %\)', x).group(1))
runs = [r for r in rows(sec['1']) if r[3] == 'SNV' and 'xs__' not in r[0] and 'norescue' not in r[0] and 'joint' not in r[0]]
print(len(runs), 'runs:', [r[0].split('__')[1] for r in runs])
for i, name in [(10, 'population'), (13, 'germline'), (14, 'tumor-only')]:
    v = [pct(r[i]) for r in runs]; print(name, min(v), max(v), 'ep8 COLO Ill:', [pct(r[i]) for r in runs if 'COLO829_Illumina_SNV_7fc7396_ddp4' in r[0]])
names = {r[0] for r in runs}
c3 = [r for r in rows(sec['3']) if r[0] in names]
def cls(s, k):
    m = re.search(rf'{k} (\d+)', s); return int(m.group(1)) if m else 0
tp = [cls(r[4], 'tumor_only') for r in c3 if r[1] == 'ctrl_tp']; fp = [cls(r[4], 'germline') for r in c3 if r[1] == 'ctrl_fp']
print('ctrl_tp tumor_only', min(tp), max(tp), len(tp), ' ctrl_fp germline', min(fp), max(fp), len(fp))
up = [(r[0], cls(r[4], 'tumor_only') / int(r[2])) for r in c3 if r[1] == 'unplaced_pass']
v = [x for _, x in up]; print('assemblies alone tumor_only %', round(100 * min(v), 1), round(100 * max(v), 1))
# section 8: keys PASS on all three platforms per tumor/graph group
g8 = [r for r in rows(sec['8']) if r[1] in names]
for r in g8:
    d = json.loads(r[5]); n_all = {k: v.get('2', 0) for k, v in d.items()}; tot = sum(n_all.values())
    print('  ', r[0], r[1].split('__')[1], 'platforms', r[2], 'all3', tot, 'pop%', round(100 * n_all.get('population', 0) / tot, 1) if tot else None, 'tumor_only', n_all.get('tumor_only', 0))
