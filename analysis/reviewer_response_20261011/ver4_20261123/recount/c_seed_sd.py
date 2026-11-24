#!/usr/bin/env python3
"""ver4 recount (c): seed standard deviations of chr1 F1 for R3.9 point 1, and the comparator count of R3.10.

    /wanglab/jshen/anaconda3/bin/python -B c_seed_sd.py

Every F1 is `pon.at_threshold.f1` of <run>/<vcf dir>/eval_bed/report.json (rtg vcfeval, chr1, benchmark regions, PASS
calls after the population filter). SD = sample standard deviation (n - 1) over seeds 0 / 1 / 2 (statistics.stdev).
INDEL (R = /scratch/jshen/data/pansoma_net_v2_runs):
  six-tumor   R/release_s{0,1,2}_<P>_INDEL
  four-tumor  R/unseen_s{0,1,2}_<P>_INDEL
  VCF dir vcf_<rule>_chr1_<SET>_d46, rule = ppr8 on Illumina / PacBio HiFi (COLO829T: Fiber-seq) and pp on ONT (ver4);
  the letter-ver3 rule set (ppr8 on all platforms, ONT too) is computed with the same code to reproduce its range.
SNV: COLO829T single-tumor models of the current pipeline, R/snv_newpipe_20261109/
  COLO829_{Illumina,fiberseq}_SNV_{d46,d9}_newpipe[_s1|_s2], VCF dir vcf_chr1_COLO829T_<P>[_d46] (rule 2).
Comparators (R3.10 "outperform both comparator callers"): E8 score, <caller>/<SET>/INDEL/report.json
  eval_bed.pass_pon.at_threshold (the tool's PASS after the same filter), as analysis/latest_results_20261123/collect.py.
Checks (asserted): TP + FN equal over the three seeds, over six-tumor / four-tumor, over ppr8 / pp, equal to
report.json pon.baseline and raw.baseline, and equal to both comparators' truth count; the per-seed F1 strings of the
ver4 rule set equal sections 1b / 3b / 6 of latest_results_20261123.md; an SD computed from F1 recomputed from the
counts differs from the report-f1 SD by < 0.0002.
Writes c_seed_sd.md, c_seed_sd.json, c_seed_sd.tsv next to this script.
"""
import csv
import json
import re
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = Path('/scratch/jshen/data/pansoma_net_v2_runs')
N = R / 'snv_newpipe_20261109'
E8 = Path('/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score')
LR = Path('/wanglab/jshen/scratch_jshen/Pansoma_papers/ver1/reviewer_response/latest_results_20261123.md')
TUMORS = ('HG008T', 'COLO829T', 'HCC1395T')
PLATS = ('Illumina', 'PacBio', 'ONT')
MODELS = (('release', 'six-tumor'), ('unseen', 'four-tumor'))
SEEDS = (0, 1, 2)
VER4 = {'Illumina': 'ppr8', 'PacBio': 'ppr8', 'ONT': 'pp'}
OLD = {p: 'ppr8' for p in PLATS}


def ds(t, p):
    return 'COLO829T_fiberseq' if (t, p) == ('COLO829T', 'PacBio') else f'{t}_{p}'


def pname(t, p):
    return {'Illumina': 'Illumina', 'ONT': 'ONT', 'PacBio': 'Fiber-seq' if t == 'COLO829T' else 'PacBio HiFi'}[p]


def report(d):
    rep = json.loads((Path(d) / 'eval_bed' / 'report.json').read_text())
    a = rep['pon']['at_threshold']
    assert a['tp_baseline'] + a['fn'] == rep['pon']['baseline'] == rep['raw']['baseline'], d
    pr = a['tp_call'] / (a['tp_call'] + a['fp']) if a['tp_call'] + a['fp'] else 0.0
    rc = a['tp_baseline'] / (a['tp_baseline'] + a['fn'])
    f1c = 2 * pr * rc / (pr + rc) if pr + rc else 0.0
    assert abs(f1c - a['f1']) < 6e-5, (d, f1c, a['f1'])          # report f1 is the count F1 rounded to 4 decimals
    return dict(f1=a['f1'], f1_counts=f1c, truth=a['tp_baseline'] + a['fn'], tp=a['tp_baseline'], fp=a['fp'], fn=a['fn'],
                precision=a['precision'], recall=a['recall'], call=(rep.get('steps', {}).get('graph_vcf', {}) or {}).get('call', ''))


def cell(cells):
    f = [c['f1'] for c in cells]
    assert len({c['truth'] for c in cells}) == 1, [c['truth'] for c in cells]
    sd, sdc = statistics.stdev(f), statistics.stdev([c['f1_counts'] for c in cells])
    assert abs(sd - sdc) < 2e-4, (sd, sdc)
    return dict(f1=f, mean=statistics.mean(f), sd=sd, truth=cells[0]['truth'])


def summary(sds):
    return dict(n=len(sds), min=min(sds), max=max(sds), median=statistics.median(sds))


def rng(s):
    return f"{s['min']:.3f}–{s['max']:.3f} (median {s['median']:.3f}; n = {s['n']})"


def lr_table(title):
    """Rows of a section of latest_results_20261123.md, keyed (tumor, platform name) -> list of cells."""
    txt = LR.read_text().split(title, 1)[1]
    txt = re.split(r'\n##+ ', txt, 1)[0]
    out = {}
    for line in txt.splitlines():
        c = [x.strip() for x in line.strip().strip('|').split('|')]
        if len(c) > 3 and c[0] in ('HCC1395T', 'COLO829T', 'HG008T', 'Illumina', 'Fiber-seq', 'ONT'):
            out[(c[0], c[1])] = c
    return out


def main():
    S, rows = {}, []
    # ------------------------------------------------------------------ INDEL
    indel = {}
    for name, rules in (('ver4', VER4), ('ver3', OLD)):
        for m, _ in MODELS:
            for t in TUMORS:
                for p in PLATS:
                    cs = [report(R / f'{m}_s{s}_{p}_INDEL' / f'vcf_{rules[p]}_chr1_{ds(t, p)}_d46') for s in SEEDS]
                    indel[(name, m, t, p)] = dict(cell(cs), rule=rules[p], calls=sorted({c['call'] for c in cs}), per_seed=cs)
    for t in TUMORS:                                   # truth equal over models and rule sets
        for p in PLATS:
            assert len({indel[(n, m, t, p)]['truth'] for n in ('ver4', 'ver3') for m, _ in MODELS}) == 1, (t, p)
    for (name, m, t, p), c in indel.items():
        rows.append(dict(kind='INDEL', rules=name, model=dict(MODELS)[m], tumor=t, platform=pname(t, p), rule=c['rule'], truth=c['truth'],
                         f1_s0=c['f1'][0], f1_s1=c['f1'][1], f1_s2=c['f1'][2], mean=f"{c['mean']:.4f}", sd=f"{c['sd']:.4f}"))

    def sds(name, models, tumors):
        return [indel[(name, m, t, p)]['sd'] for m in models for t in tumors for p in PLATS]
    two, three = ('HG008T', 'COLO829T'), TUMORS
    for name in ('ver4', 'ver3'):
        S[f'INDEL_{name}'] = {
            'six-tumor + four-tumor, HG008T + COLO829T (12)': summary(sds(name, ('release', 'unseen'), two)),
            'six-tumor, HG008T + COLO829T (6)': summary(sds(name, ('release',), two)),
            'four-tumor, HG008T + COLO829T (6)': summary(sds(name, ('unseen',), two)),
            'six-tumor + four-tumor, three tumors (18)': summary(sds(name, ('release', 'unseen'), three)),
            'six-tumor, three tumors (9)': summary(sds(name, ('release',), three)),
            'four-tumor, three tumors (9)': summary(sds(name, ('unseen',), three)),
            'six-tumor, HCC1395T (3)': summary(sds(name, ('release',), ('HCC1395T',))),
            'four-tumor, HCC1395T (3)': summary(sds(name, ('unseen',), ('HCC1395T',)))}
    o = S['INDEL_ver3']['six-tumor + four-tumor, HG008T + COLO829T (12)']
    S['reproduces_letter_ver3_indel_sd'] = (f"{o['min']:.3f}", f"{o['max']:.3f}", f"{o['median']:.3f}") == ('0.003', '0.027', '0.014')
    assert S['reproduces_letter_ver3_indel_sd'], o

    # latest_results_20261123.md: the per-seed F1 strings of sections 1b (six-tumor) and 3b (four-tumor)
    lr1, lr3 = lr_table('### 1b. INDEL'), lr_table('### 3b. INDEL')
    n_lr = 0
    for t in TUMORS:
        for p in PLATS:
            c6, c4 = indel[('ver4', 'release', t, p)], indel[('ver4', 'unseen', t, p)]
            a = lr1[(t, pname(t, p))]
            assert a[6] == ' / '.join(f'{x:.3f}' for x in c6['f1']) and a[7] == f"{c6['mean']:.3f} ± {c6['sd']:.3f}", (t, p, a[6:], c6)
            assert int(a[2]) == c6['truth']
            b = lr3[(t, pname(t, p))]
            assert b[6] == ' / '.join(f'{x:.3f}' for x in c4['f1']) + f" ({c4['mean']:.3f})", (t, p, b[6], c4)
            n_lr += 2
    S['check_latest_results_indel_cells'] = n_lr

    # ------------------------------------------------------------------ SNV: COLO829T single-tumor, three seeds
    snv = {}
    lr6 = lr_table('## 6. COLO829T single-tumor')
    for p in ('Illumina', 'fiberseq'):
        for g, gname in (('d46', 'HPRC v2.1'), ('d9', 'HPRC v1.1')):
            cs = []
            for s in SEEDS:
                run = N / f"COLO829_{p}_SNV_{g}_newpipe{'' if s == 0 else f'_s{s}'}"
                d = sorted(run.glob(f'vcf_chr1_COLO829T_{p}*'))
                assert len(d) == 1, (run, d)
                cs.append(dict(report(d[0]), vcf=str(d[0])))
            c = dict(cell(cs), per_seed=cs)
            snv[(p, g)] = c
            a = lr6[('Fiber-seq' if p == 'fiberseq' else 'Illumina', gname)]
            assert [x.split(' / ')[2] for x in a[2:5]] == [f'{x:.3f}' for x in c['f1']] and a[5] == f"{c['mean']:.3f}", (p, g, a)
            rows.append(dict(kind='SNV', rules='rule 2', model='COLO829T single-tumor', tumor='COLO829T', platform='Fiber-seq' if p == 'fiberseq' else p,
                             rule=gname, truth=c['truth'], f1_s0=c['f1'][0], f1_s1=c['f1'][1], f1_s2=c['f1'][2], mean=f"{c['mean']:.4f}", sd=f"{c['sd']:.4f}"))
    for p in ('Illumina', 'fiberseq'):
        assert snv[(p, 'd46')]['truth'] == snv[(p, 'd9')]['truth'], p
    S['SNV_COLO829T_single_tumor'] = summary([c['sd'] for c in snv.values()])
    S['SNV_COLO829T_single_tumor_v21'] = summary([snv[(p, 'd46')]['sd'] for p in ('Illumina', 'fiberseq')])
    S['SNV_COLO829T_single_tumor_v11'] = summary([snv[(p, 'd9')]['sd'] for p in ('Illumina', 'fiberseq')])
    assert f"{snv[('Illumina', 'd46')]['sd']:.3f}" == '0.048'      # letter_ver3_audit.md: "The seed SD of COLO829T Illumina on v2.1 is now 0.048"
    # which six-tumor SNV seeds exist (for the statement that no six-tumor SNV seed SD can be given)
    S['six_tumor_SNV_seed_runs_with_chr1_vcf'] = {
        str(d.name): len(list(d.glob('vcf_chr1_*/eval_bed/report.json'))) for d in sorted(N.glob('multicancer6_*_SNV_newpipe*'))}

    # ------------------------------------------------------------------ R3.10: six-tumor INDEL against both comparators
    wins = {}
    comp_rows = []
    for t in TUMORS:
        for p in PLATS:
            c = indel[('ver4', 'release', t, p)]
            e = {k: json.loads((E8 / k / ds(t, p) / 'INDEL' / 'report.json').read_text())['eval_bed']['pass_pon']['at_threshold'] for k in ('deepsomatic', 'clairs_to')}
            for a in e.values():
                assert a['tp_baseline'] + a['fn'] == c['truth'], (t, p)
            best = max(a['f1'] for a in e.values())
            wins[(t, p)] = dict(seed0=c['f1'][0] > best, mean=c['mean'] > best, n_seeds=sum(x > best for x in c['f1']))
            comp_rows.append((t, p, c, e['deepsomatic']['f1'], e['clairs_to']['f1']))
    S['six_tumor_INDEL_beats_both_comparators'] = {
        'HG008T + COLO829T (6), seed 0': sum(wins[(t, p)]['seed0'] for t in two for p in PLATS),
        'HG008T + COLO829T (6), seed mean': sum(wins[(t, p)]['mean'] for t in two for p in PLATS),
        'three tumors (9), seed 0': sum(wins[(t, p)]['seed0'] for t in three for p in PLATS),
        'three tumors (9), seed mean': sum(wins[(t, p)]['mean'] for t in three for p in PLATS)}
    lr_wins = re.search(r'### 1b\. INDEL.*?Pansoma has the highest chr1 F1 in (\d+) of 9 cells \(seed 0\)', LR.read_text(), re.S).group(1)
    assert int(lr_wins) == S['six_tumor_INDEL_beats_both_comparators']['three tumors (9), seed 0']

    S['INDEL_cells'] = {f'{n}|{m}|{t}|{p}': {k: v for k, v in c.items() if k != 'per_seed'} for (n, m, t, p), c in indel.items()}
    S['SNV_cells'] = {f'{p}|{g}': dict(f1=c['f1'], mean=c['mean'], sd=c['sd'], truth=c['truth'], vcf=[x['vcf'] for x in c['per_seed']]) for (p, g), c in snv.items()}
    (HERE / 'c_seed_sd.json').write_text(json.dumps(S, indent=1) + '\n')
    with open(HERE / 'c_seed_sd.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # ------------------------------------------------------------------ tables
    f3 = lambda xs: ' / '.join(f'{x:.3f}' for x in xs)  # noqa: E731
    L = ['# ver4 recount (c): seed SD of chr1 F1 (R3.9 point 1) and the comparator count of R3.10', '',
         'Generated by `c_seed_sd.py`. chr1, benchmark regions, PASS after the population filter; F1 = report.json '
         '`pon.at_threshold.f1`; SD = sample SD (n - 1) over seeds 0 / 1 / 2.', '',
         '## 1. INDEL models on HPRC v2.1 (PASS: ppr8 on Illumina / PacBio HiFi / Fiber-seq, pp on ONT)', '',
         '| tumor | platform | chr1 truth | six-tumor chr1 F1 seeds 0 / 1 / 2 | six-tumor chr1 mean ± SD | four-tumor chr1 F1 seeds 0 / 1 / 2 | four-tumor chr1 mean ± SD |',
         '|---|---|---:|---|---|---|---|']
    for t in TUMORS:
        for p in PLATS:
            a, b = indel[('ver4', 'release', t, p)], indel[('ver4', 'unseen', t, p)]
            L.append(f"| {t} | {pname(t, p)} | {a['truth']} | {f3(a['f1'])} | {a['mean']:.3f} ± {a['sd']:.3f} | {f3(b['f1'])} | {b['mean']:.3f} ± {b['sd']:.3f} |")
    L += ['', '## 2. Range and median of the INDEL seed SDs', '',
          '| SDs pooled over | ver4 rules (ONT at pp) | letter-ver3 rules (ONT at ppr8) |', '|---|---|---|']
    for k in S['INDEL_ver4']:
        L.append(f"| {k} | {rng(S['INDEL_ver4'][k])} | {rng(S['INDEL_ver3'][k])} |")
    L += ['', 'The letter-ver3 sentence "0.003–0.027 (median 0.014) for final INDEL models" is the first row of the right column (reproduced, asserted).', '',
          '## 3. SNV: COLO829T single-tumor models of the current pipeline (rule 2), COLO829T chr1', '',
          '| platform | graph | chr1 truth | chr1 F1 seeds 0 / 1 / 2 | chr1 mean ± SD |', '|---|---|---:|---|---|']
    for (p, g), c in snv.items():
        L.append(f"| {'Fiber-seq' if p == 'fiberseq' else p} | {'HPRC v2.1' if g == 'd46' else 'HPRC v1.1'} | {c['truth']:,} | {f3(c['f1'])} | {c['mean']:.3f} ± {c['sd']:.3f} |")
    L += ['', f"SD range over the four model sets: {rng(S['SNV_COLO829T_single_tumor'])}; HPRC v2.1 only: "
          f"{rng(S['SNV_COLO829T_single_tumor_v21'])}; HPRC v1.1 only: {rng(S['SNV_COLO829T_single_tumor_v11'])}.",
          '', 'Six-tumor SNV models: chr1 call sets per run directory (seed 0 = no suffix): ' +
          ', '.join(f'`{k}` {v}' for k, v in S['six_tumor_SNV_seed_runs_with_chr1_vcf'].items()) +
          '. No six-tumor SNV seed SD can be computed on this cluster.', '',
          '## 4. Six-tumor INDEL models against both comparators (R3.10), chr1 F1 after the filter', '',
          '| tumor | platform | six-tumor chr1 F1 seed 0 | six-tumor chr1 seed mean | DeepSomatic-TO chr1 | ClairS-TO chr1 | above both: seed 0 | above both: seed mean | seeds above both |',
          '|---|---|---|---|---|---|---|---|---|']
    for t, p, c, d, e in comp_rows:
        w = wins[(t, p)]
        L.append(f"| {t} | {pname(t, p)} | {c['f1'][0]:.3f} | {c['mean']:.3f} | {d:.3f} | {e:.3f} | {'yes' if w['seed0'] else 'no'} | {'yes' if w['mean'] else 'no'} | {w['n_seeds']} of 3 |")
    L += ['', 'Cells where the six-tumor model is above both comparators: ' + '; '.join(f'{k}: {v}' for k, v in S['six_tumor_INDEL_beats_both_comparators'].items()) + '.', '',
          '## 5. Checks', '',
          '- Every report.json: TP + FN = `pon.baseline` = `raw.baseline`; `f1` equals the F1 recomputed from the counts to 4 decimals; the SD from the recomputed F1 differs by < 0.0002.',
          '- Truth counts equal over the three seeds, six-tumor / four-tumor, ppr8 / pp call sets of a tumor x platform, and equal to both comparators\' truth count; '
          'single-tumor SNV: equal over seeds and over v1.1 / v2.1 per platform.',
          f"- The {n_lr} per-seed INDEL F1 strings (and six-tumor mean ± SD, truth) equal `latest_results_20261123.md` §1b / §3b; the 4 single-tumor SNV rows equal its §6; "
          f"the seed-0 comparator count on nine cells equals its \"{lr_wins} of 9\".",
          '- Old number reproduced with the same code: INDEL SD range 0.003–0.027, median 0.014 (12 SDs, ONT at ppr8) = letter ver3; '
          'COLO829T Illumina v2.1 SNV seed SD 0.048 = letter_ver3_audit.md.',
          '- PASS rule strings of the ver4 INDEL call sets (report.json `steps.graph_vcf.call`, text in brackets): ' +
          '; '.join(f"{r}: {' | '.join(sorted({c.split(' (', 1)[1][:-1] for (n, m, t, p), v in indel.items() if n == 'ver4' and v['rule'] == r for c in v['calls']}))}" for r in ('ppr8', 'pp')) + '.']
    (HERE / 'c_seed_sd.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
