import json
E = '/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators'
DS = ('HG008T_Illumina', 'HG008T_PacBio', 'HG008T_ONT', 'HCC1395T_Illumina', 'HCC1395T_PacBio', 'HCC1395T_ONT',
      'HCC1395T_Illumina_chr20', 'HCC1395T_PacBio_chr20', 'HCC1395T_ONT_chr20')
def cell(r, s):
    a = r['eval_bed'][s]['at_threshold']
    return f"{a['precision']:.3f} / {a['recall']:.3f} / **{a['f1']:.3f}** ({a['tp_baseline']}/{a['fp']}/{a['fn']})"
out = {}
for k in ('SNV', 'INDEL'):
    print(f'\n## chr1 / chr20 {k}, in the benchmark BED: P / R / F1 (TP/FP/FN)\n')
    print('| dataset | caller | PASS without any population filtering | PASS with the tool\'s own filter (default) | default + our PoN | PASS records without / with own filter |')
    print('|---|---|---|---|---|---|')
    for d in DS:
        for c, name in (('clairs_to', 'ClairS-TO'), ('deepsomatic', 'DeepSomatic-TO')):
            n = json.load(open(f'{E}/score_nopon/{c}/{d}/{k}/report.json'))
            o = json.load(open(f'{E}/score/{c}/{d}/{k}/report.json'))
            print(f"| {d} | {name} | {cell(n, 'pass')} | {cell(o, 'pass')} | {cell(o, 'pass_pon')} | {n['counts']['pass']:,} / {o['counts']['pass']:,} |")
            out[f'{c}/{d}/{k}'] = {s: {m: n['eval_bed']['pass']['at_threshold'][m] if s == 'nopon' else o['eval_bed'][s]['at_threshold'][m]
                                       for m in ('precision', 'recall', 'f1', 'tp_baseline', 'fp', 'fn')}
                                   for s in ('nopon', 'pass', 'pass_pon')}
json.dump(out, open('summary.json', 'w'), indent=1)
