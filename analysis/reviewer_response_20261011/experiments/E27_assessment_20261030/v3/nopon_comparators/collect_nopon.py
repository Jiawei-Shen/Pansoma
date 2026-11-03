"""ClairS-TO / DeepSomatic-TO without their own population filters vs their default (E8 scoring, benchmark BED, PASS).
Without: ClairS-TO rerun with --disable_nonsomatic_tagging (clairs_to_rerun/, E8_comparators/score_nopon_rerun/);
DeepSomatic-TO PON -> PASS (rewrite_nopon.py, E8_comparators/score_nopon/; exact: postprocess_variants applies the PoN
last and only to PASS calls). The ClairS-TO FILTER rewrite (NonSomatic with a PoN flag -> PASS, score_nopon/) is not
equivalent to the rerun (the tagging runs before realignment, postfilter and haplotype filtering, which process PASS
records only); it is listed separately.
usage: python collect_nopon.py > tables.md   (also writes summary.json)"""
import json
E = '/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators'
DS = ('HG008T_Illumina', 'HG008T_PacBio', 'HG008T_ONT', 'COLO829T_Illumina', 'COLO829T_fiberseq', 'COLO829T_ONT',
      'HCC1395T_Illumina', 'HCC1395T_PacBio', 'HCC1395T_ONT',
      'HCC1395T_Illumina_chr20', 'HCC1395T_PacBio_chr20', 'HCC1395T_ONT_chr20')
NOPON = {'clairs_to': 'score_nopon_rerun', 'deepsomatic': 'score_nopon'}
M = ('precision', 'recall', 'f1', 'tp_baseline', 'fp', 'fn')


def cell(r, s='pass'):
    a = r['eval_bed'][s]['at_threshold']
    return f"{a['precision']:.3f} / {a['recall']:.3f} / **{a['f1']:.3f}** ({a['tp_baseline']}/{a['fp']}/{a['fn']})"


def load(src, c, d, k):
    return json.load(open(f'{E}/{src}/{c}/{d}/{k}/report.json'))


print('# Comparators without their own population filters (chr1; HCC1395T also chr20)')
print('\nClairS-TO: rerun with --disable_nonsomatic_tagging. DeepSomatic-TO: PON -> PASS (its GERMLINE class and population-AF '
      'input stay). Default + our PoN = E8 pass_pon.')
out = {}
for k in ('SNV', 'INDEL'):
    print(f'\n## {k}, in the benchmark BED: P / R / F1 (TP/FP/FN)\n')
    print('| dataset | caller | PASS without any population filtering | PASS with the tool\'s own filter (default) | default + our PoN | PASS records without / with own filter |')
    print('|---|---|---|---|---|---|')
    for d in DS:
        for c, name in (('clairs_to', 'ClairS-TO'), ('deepsomatic', 'DeepSomatic-TO')):
            n, o = load(NOPON[c], c, d, k), load('score', c, d, k)
            print(f"| {d} | {name} | {cell(n)} | {cell(o)} | {cell(o, 'pass_pon')} | {n['counts']['pass']:,} / {o['counts']['pass']:,} |")
            out[f'{c}/{d}/{k}'] = {s: {m: (n['eval_bed']['pass'] if s == 'nopon' else o['eval_bed'][s])['at_threshold'][m] for m in M}
                                   for s in ('nopon', 'pass', 'pass_pon')}
            if c == 'clairs_to':
                w = load('score_nopon', c, d, k)
                out[f'{c}/{d}/{k}']['nopon_filter_rewrite'] = {m: w['eval_bed']['pass']['at_threshold'][m] for m in M}
print('\n## ClairS-TO without its PoN step: rerun vs FILTER rewrite (not equivalent), P / R / F1 (TP/FP/FN), PASS records\n')
print('| dataset | kind | rerun --disable_nonsomatic_tagging | rewrite NonSomatic (PoN flag) -> PASS |')
print('|---|---|---|---|')
for k in ('SNV', 'INDEL'):
    for d in DS:
        r, w = load('score_nopon_rerun', 'clairs_to', d, k), load('score_nopon', 'clairs_to', d, k)
        print(f"| {d} | {k} | {cell(r)} {r['counts']['pass']:,} | {cell(w)} {w['counts']['pass']:,} |")
json.dump(out, open('summary.json', 'w'), indent=1)
