#!/usr/bin/env python3
"""Single-tumor INDEL ppr8 checkpoints (RIS, origin analysis/indel_ppr8_9919d4a_20261026) scored on HCC1395T.

Reads every <run>/vcf_<chrom>_<set>/eval_bed/report.json under ris_ppr8_20261026 and writes, next to this script,
results.tsv (one row per checkpoint and test set) and tables.md (seed-0 own-chr1 checks against the RIS scores; 3-seed
means on HCC1395T chr1 / chr20). P / R / F1 = after the population filter, benchmark BED, PASS at the checkpoint's
call rule (post-PoN validation recall 0.8).
"""
import glob
import json
import os
import re
import sys
from collections import defaultdict

R = '/scratch/jshen/data/pansoma_net_v2_runs/ris_ppr8_20261026'
HERE = os.path.dirname(os.path.abspath(__file__))
# RIS seed-0 scores on the model's own chr1 (origin analysis/indel_ppr8_9919d4a_20261026/summary.json, single.<m>.d46.s0)
RIS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '../../../../../indel_ppr8_9919d4a_20261026/summary.json')

rows = []
for f in sorted(glob.glob(f'{R}/*/vcf_chr*_*/eval_bed/report.json')):
    rep = json.load(open(f))
    d = rep['pon']['at_threshold']
    run, v = f.split('/')[-4], f.split('/')[-3]
    m = re.match(r'(.*)_INDEL_d46_ppr8(?:_s(\d))?$', run)
    model, seed = m.group(1), int(m.group(2) or 0)
    chrom, tset = re.match(r'vcf_(chr\d+)_(.*)', v).groups()
    rows.append(dict(model=model, seed=seed, set=tset, chrom=chrom, p=d['precision'], r=d['recall'], f1=d['f1'],
                     tp=d['tp_baseline'], fp=d['fp'], truth=rep['pon']['baseline']))

with open(f'{HERE}/results.tsv', 'w') as o:
    o.write('model\tseed\tset\tchrom\tprecision\trecall\tf1\ttp\tfp\ttruth\n')
    for x in rows:
        o.write('\t'.join(str(x[k]) for k in ('model', 'seed', 'set', 'chrom', 'p', 'r', 'f1', 'tp', 'fp', 'truth')) + '\n')

out = ['# Single-tumor INDEL ppr8 checkpoints (RIS) on HCC1395T', '',
       'chr1 / chr20, benchmark BED, after the population filter, PASS at the checkpoint call rule '
       '(post-PoN validation recall 0.8). Checkpoints: RIS `analysis/indel_ppr8_9919d4a_20261026` (3 seeds), '
       'predicted on a tequila MIG slice from the repository at 3576902.', '']
out += ['## Seed 0 on its own chr1 (check against RIS)', '', '| model | set | P / R / F1 | TP / FP | RIS P / R / F1 | RIS TP / FP |',
        '|---|---|---|---|---|---|']
ris = json.load(open(RIS))['single'] if RIS else {}
for x in rows:
    if x['seed'] == 0 and not x['set'].startswith('HCC1395T'):
        s = ris.get(x['model'], {}).get('d46', {}).get('s0')
        rs = f"{s['p']:.3f} / {s['r']:.3f} / {s['f1']:.3f} | {s['tp']} / {s['fp']}" if s else '- | -'
        out.append(f"| {x['model']} | {x['set']} | {x['p']:.3f} / {x['r']:.3f} / {x['f1']:.3f} | {x['tp']} / {x['fp']} | {rs} |")
out += ['', '## HCC1395T, mean over seeds 0-2 (Supplementary Table R1.7c uses chr1)', '',
        '| model | set | chrom | n | P / R / F1 (mean) | F1 range | FP per seed |', '|---|---|---|---|---|---|---|']
g = defaultdict(list)
for x in rows:
    if x['set'].startswith('HCC1395T'):
        g[(x['set'], x['chrom'], x['model'])].append(x)
for (tset, chrom, model), L in sorted(g.items(), key=lambda k: (k[0][0], int(k[0][1][3:]), k[0][2])):
    n = len(L)
    mean = {k: sum(x[k] for x in L) / n for k in ('p', 'r', 'f1')}
    out.append(f"| {model} | {tset} | {chrom} | {n} | {mean['p']:.3f} / {mean['r']:.3f} / {mean['f1']:.3f} | "
               f"{min(x['f1'] for x in L):.3f}-{max(x['f1'] for x in L):.3f} | "
               f"{', '.join(str(x['fp']) for x in sorted(L, key=lambda y: y['seed']))} |")
open(f'{HERE}/tables.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
