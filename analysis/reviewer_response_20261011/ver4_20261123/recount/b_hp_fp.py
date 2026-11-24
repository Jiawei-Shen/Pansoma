#!/usr/bin/env python3
"""ver4 recount (b): R3.10 point 4. Share of HG008T chr1 false-positive INDEL calls in homopolymers >= 7 bp for the
six-tumor release INDEL models, and the population-allele counts of the HG008T chr1 truth INDELs.

    PYTHONDONTWRITEBYTECODE=1 /wanglab/jshen/anaconda3/bin/python -B b_hp_fp.py

Source of the letter's "44-70% of HG008T false calls": E17a_indel_strata/final_models.md section 6 (e17a_final_models.py,
2026-10-29): release_s{0,1,2}_<P>_INDEL, PASS = ppr8 on ALL three platforms (ONT too), seeds pooled (FP summed over the
three seeds), chr1, benchmark BED, after the population filter. Here the same counting (E17a's strata() and Classifier,
imported unchanged; call-side FP = eval_bed/vcfeval/pon_calls/fp.vcf.gz, class of the call by
truth_labels.RepeatContext on GRCh38: HP>=7 = the inserted / deleted bases are one repeated base in a run >= 7 bp) is
run
  1. with the old rule set (ppr8 everywhere) -> must give final_models.tsv row by row (asserted);
  2. with the ver4 PASS rules: ppr8 on Illumina / PacBio HiFi, pp on ONT (vcf_pp_chr1_HG008T_ONT_d46).
Truth side (the same for every call set): the chr1 truth INDELs in the BED = tp-baseline + fn records of pon_calls; the
population tag is read from eval_bed/truth.INDEL.pon.vcf.gz (written by vcfeval with scripts/filter_panel_of_normals.py)
and RE-DERIVED here from the two panel VCFs: exact POS / REF / ALT in gnomAD or CoLoRSdb with that ALT's AF >= 0.01 when
the ALT changes a homopolymer of >= 8 bp (filter_panel_of_normals.homopolymer_run on the calls' ##reference FASTA), else
AF >= 1e-4. The AF lookup below is written here (not the filter's record_matches); the filter's own functions are run as
a second check. Both must equal the file's tags (asserted).
Writes b_hp_fp.md, b_hp_fp.json, b_hp_fp_per_seed.tsv, b_truth_indel_pon.tsv next to this script.
"""
import csv
import gzip
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
REPO = Path('/scratch/jshen/Github/Pansoma')
E17A = REPO / 'analysis/reviewer_response_20261011/experiments/E17a_indel_strata'
sys.path.insert(0, str(E17A))
sys.path.insert(0, str(REPO / 'scripts'))
import e17a_indel_strata as E  # noqa: E402
from e17a_indel_strata import C  # noqa: E402
import filter_panel_of_normals as F  # noqa: E402

RUNS = Path('/scratch/jshen/data/pansoma_net_v2_runs')
PLATS = ('Illumina', 'PacBio', 'ONT')
PNAME = {'Illumina': 'Illumina', 'PacBio': 'PacBio HiFi', 'ONT': 'ONT'}
SEEDS = (0, 1, 2)
VER4_RULE = {'Illumina': 'ppr8', 'PacBio': 'ppr8', 'ONT': 'pp'}


def unit(p, seed, rule):
    run = f'release_s{seed}_{p}_INDEL'
    d = RUNS / run / f'vcf_{rule}_chr1_HG008T_{p}_d46' / 'eval_bed'
    w = d / 'vcfeval' / 'pon_calls'
    assert (d / 'report.json').exists() and C._rtg_complete(w), d
    u = C.Unit(uid=f'{run}/vcf_{rule}_chr1_HG008T_{p}_d46', run=run, set=f'HG008T_{p}_d46', kind='INDEL', role='release',
               variant=f's{seed}', model=run)
    u.rtg['bed'], u.report['bed'] = w, json.loads((d / 'report.json').read_text())
    u.truth_pon = d / 'truth.INDEL.pon.vcf.gz'
    return u


def count(u, classify):
    t, skipped, has_pon = E.strata(u, classify)
    rep = u.report['bed']
    at = rep['pon']['at_threshold']
    a = t[('all', 'all')]
    assert (a[0], a[1], a[2], a[3]) == (at['tp_baseline'], at['fn'], at['tp_call'], at['fp']), u.uid   # report.json
    assert a[0] + a[1] == rep['pon']['baseline'] == rep['raw']['baseline'], u.uid
    hp = t.get(('context', 'HP>=7'), [0] * 5)
    return dict(uid=u.uid, call=rep['steps']['graph_vcf']['call'], threshold=rep['threshold'], truth=a[0] + a[1],
                truth_tagged=a[4], truth_hp7=hp[0] + hp[1], truth_hp7_tagged=hp[4], tp=a[0], fn=a[1], tp_call=a[2], fp=a[3],
                fp_hp7=hp[3], tp_call_hp7=hp[2], skipped=skipped,
                fp_by_context={c: t.get(('context', c), [0] * 5)[3] for c in E.CONTEXTS})


def af_of(pon, contig, pos, ref, alt):
    """Highest AF of the exact allele (POS, REF, ALT) in an indexed panel VCF; None when the allele is not there."""
    best = None
    for rec in pon.fetch(contig, max(0, pos - 1), pos):
        if rec.pos != pos or (rec.ref or '').upper() != ref:
            continue
        alts = [(a or '').upper() for a in rec.alts or ()]
        af = rec.info.get('AF')
        af = af if isinstance(af, tuple) else (af,) * len(alts)
        for k, a in enumerate(alts):
            if a == alt and af[k] is not None:
                best = af[k] if best is None else max(best, af[k])
    return best


def truth_side(u, classify, out_tsv):
    """Population tags of the chr1 truth INDELs in the BED: the file's, and re-derived from the panels with the rule."""
    w = u.rtg['bed']
    inbed = [(r[0], r[1], r[2], r[3]) for n in ('tp-baseline', 'fn') for r in C.read_records(w / f'{n}.vcf.gz')]
    assert len(inbed) == len(set(inbed))
    rep = u.report['bed']
    fasta_path = None
    with gzip.open(rep['calls'], 'rt') as f:
        for line in f:
            if not line.startswith('##'):
                break
            if line.startswith('##reference=file://'):
                fasta_path = line[len('##reference=file://'):].strip()
    fasta = pysam.FastaFile(fasta_path)
    pons = dict(zip(F.PON_NAMES, (pysam.VariantFile(p) for p in rep['pon_files'])))
    assert 'gnomad' in rep['pon_files'][0].lower() and 'colors' in rep['pon_files'][3].lower()
    assert (F.INDEL_HOMOPOLYMER, F.INDEL_MIN_AF_HOMOPOLYMER, F.INDEL_MIN_AF) == (8, 0.01, 0.0001)
    assert F.PON_RULES['INDEL'] == ('indel_af', None, None, 'indel_af')
    file_rows = {}
    with gzip.open(u.truth_pon, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            c = line.rstrip('\n').split('\t', 8)
            alts = c[4].upper().split(',')                 # key = first ALT, as E17a / callsets.read_records
            assert (c[0], int(c[1]), c[3].upper(), alts[0]) not in file_rows
            file_rows[(c[0], int(c[1]), c[3].upper(), alts[0])] = (c[6], c[7], alts)
    assert len(file_rows) == rep['truth_records'] and set(inbed) <= set(file_rows)
    n_multi = sum(len(v[2]) > 1 for v in file_rows.values())
    rows = []
    for key in sorted(file_rows, key=lambda k: (k[1], k[2], k[3])):
        chrom, pos, ref, alt = key
        flt, info, alts = file_rows[key]
        tag_file = 'PanelOfNormals' in flt.split(';')
        pon_file = next((x.split('=', 1)[1] for x in info.split(';') if x.startswith('PANSOMA_PON=')), '')
        hp = F.homopolymer_run(fasta, chrom, pos, ref, alt)
        floor = 0.01 if hp >= 8 else 1e-4
        afs = {n: af_of(pons[n], F.resolve_contig(pons[n], chrom), pos, ref, alt) for n in ('PoN1_gnomAD', 'PoN4_CoLoRSdb')}
        hit, hit_filter = [], []
        for n in ('PoN1_gnomAD', 'PoN4_CoLoRSdb'):      # the filter tags a record when any of its ALTs passes its own floor
            ok = False
            for a in alts:
                fl = 0.01 if F.homopolymer_run(fasta, chrom, pos, ref, a) >= 8 else 1e-4
                x = af_of(pons[n], F.resolve_contig(pons[n], chrom), pos, ref, a)
                ok = ok or (x is not None and x >= fl * (1 - 1e-6))   # float32 AF tolerance as the filter
            if ok:
                hit.append(n)
            if any(F.record_matches(pons[n], F.resolve_contig(pons[n], chrom), pos, ref, (a,),
                                    min_af=F.indel_min_af(fasta, chrom, pos, ref, a)) for a in alts):
                hit_filter.append(n)
        assert hit == hit_filter, key
        assert bool(hit) == tag_file and ','.join(hit) == pon_file, (key, hit, flt, pon_file)
        cls = classify(chrom, pos, ref, alt)
        rows.append(dict(chrom=chrom, pos=pos, ref=ref, alt=alt, in_bed=key in set_inbed(inbed), type=cls[0] if cls else '',
                         length_bin=cls[1] if cls else '', context=cls[2] if cls else '', filter_homopolymer_run=hp,
                         af_floor=floor, gnomAD_AF=afs['PoN1_gnomAD'], CoLoRSdb_AF=afs['PoN4_CoLoRSdb'],
                         tagged=bool(hit), panels=','.join(hit), n_alts=len(alts)))
    with open(out_tsv, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        wr.writeheader()
        for r in rows:
            wr.writerow({k: ('' if v is None else f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    b = [r for r in rows if r['in_bed']]
    tg = [r for r in b if r['tagged']]
    return dict(fasta=fasta_path, panels=[rep['pon_files'][0], rep['pon_files'][3]], chr1_truth=len(rows), multi_alt_records=n_multi,
                multi_alt_in_bed=sum(r['n_alts'] > 1 for r in b),
                chr1_tagged=sum(r['tagged'] for r in rows), in_bed=len(b), in_bed_tagged=len(tg),
                in_bed_context=dict(Counter(r['context'] for r in b)), in_bed_tagged_context=dict(Counter(r['context'] for r in tg)),
                in_bed_tagged_by_panel=dict(gnomAD=sum('PoN1_gnomAD' in r['panels'] for r in tg), CoLoRSdb=sum('PoN4_CoLoRSdb' in r['panels'] for r in tg)),
                in_bed_filter_hp_ge8=sum(r['filter_homopolymer_run'] >= 8 for r in b),
                in_bed_tagged_filter_hp_ge8=sum(r['filter_homopolymer_run'] >= 8 for r in tg),
                in_bed_tagged_filter_hp_ge7=sum(r['filter_homopolymer_run'] >= 7 for r in tg),
                in_bed_tagged_floor_001=sum(r['af_floor'] == 0.01 for r in tg), in_bed_tagged_floor_1e4=sum(r['af_floor'] != 0.01 for r in tg),
                report_truth_pon_tags=rep['truth_pon_tags'], recall_bound=1 - len(tg) / len(b))


_INBED = {}


def set_inbed(inbed):
    if id(inbed) not in _INBED:
        _INBED[id(inbed)] = set(inbed)
    return _INBED[id(inbed)]


def main():
    units = {(p, s, r): unit(p, s, r) for p in PLATS for s in SEEDS for r in sorted({'ppr8', VER4_RULE[p]})}
    classify = E.Classifier(C.fasta_path(list(units.values())))
    cnt = {k: count(u, classify) for k, u in units.items()}

    # ---- check 1: the old rule set reproduces E17a final_models.tsv (per seed) and final_models.md section 6 (pooled)
    old = {}
    for r in csv.DictReader(open(E17A / 'final_models.tsv'), delimiter='\t'):
        if r['group'] == 'release' and r['dataset'].startswith('HG008T_') and (r['axis'], r['cls']) in (('all', 'all'), ('context', 'HP>=7')):
            old[(r['dataset'].split('_')[1], int(r['seed'][1]), r['cls'])] = r
    n_old = 0
    for p in PLATS:
        for s in SEEDS:
            c = cnt[(p, s, 'ppr8')]
            a, h = old[(p, s, 'all')], old[(p, s, 'HP>=7')]
            assert (c['tp'], c['fn'], c['fp'], c['truth_tagged']) == (int(a['tp_baseline']), int(a['fn']), int(a['fp']), int(a['truth_pon_tagged'])), (p, s)
            assert (c['fp_hp7'], c['truth_hp7'], c['truth_hp7_tagged']) == (int(h['fp']), int(h['truth']), int(h['truth_pon_tagged'])), (p, s)
            n_old += 1
    md6 = {}
    for line in (E17A / 'final_models.md').read_text().split('## 6. False positives in long homopolymers')[1].splitlines():
        c = [x.strip() for x in line.strip().strip('|').split('|')]
        if len(c) == 6 and c[0].startswith('HG008T_') and c[1] == 'Pansoma release':
            md6[c[0].split('_')[1]] = c[3]
    pooled_old = {p: sum(cnt[(p, s, 'ppr8')]['fp_hp7'] for s in SEEDS) / sum(cnt[(p, s, 'ppr8')]['fp'] for s in SEEDS) for p in PLATS}
    for p in PLATS:
        assert f'{pooled_old[p]:.3f}' == md6[p], (p, pooled_old[p], md6[p])

    # ---- truth side (one call set; then asserted identical for all)
    T = truth_side(units[('Illumina', 0, 'ppr8')], classify, HERE / 'b_truth_indel_pon.tsv')
    for k, c in cnt.items():
        assert (c['truth'], c['truth_tagged'], c['truth_hp7'], c['truth_hp7_tagged']) == \
            (T['in_bed'], T['in_bed_tagged'], T['in_bed_context'].get('HP>=7', 0), T['in_bed_tagged_context'].get('HP>=7', 0)), k

    # ---- tables
    res = {}
    rows = []
    for label, rules in (('letter ver3 (E17a final_models: ppr8 on all platforms)', {p: 'ppr8' for p in PLATS}),
                         ('ver4 (ppr8 Illumina / PacBio HiFi, pp ONT)', VER4_RULE)):
        res[label] = {}
        for p in PLATS:
            cs = [cnt[(p, s, rules[p])] for s in SEEDS]
            shares = [c['fp_hp7'] / c['fp'] for c in cs]
            res[label][p] = dict(rule=rules[p], fp_hp7=[c['fp_hp7'] for c in cs], fp=[c['fp'] for c in cs], share=shares,
                                 share_seed_mean=statistics.mean(shares), share_pooled=sum(c['fp_hp7'] for c in cs) / sum(c['fp'] for c in cs),
                                 tp=[c['tp'] for c in cs], calls=[c['tp_call'] + c['fp'] for c in cs],
                                 precision_hp7_pooled=sum(c['tp_call_hp7'] for c in cs) / max(1, sum(c['tp_call_hp7'] + c['fp_hp7'] for c in cs)),
                                 call=[c['call'] for c in cs], fp_by_context=[c['fp_by_context'] for c in cs])
            for s, c in zip(SEEDS, cs):
                rows.append(dict(rules=label, platform=p, seed=s, rule=rules[p], uid=c['uid'], truth=c['truth'], tp=c['tp'], fn=c['fn'], fp=c['fp'],
                                 fp_hp7=c['fp_hp7'], share_fp_hp7=f"{c['fp_hp7'] / c['fp']:.4f}", **{f'fp_{k}': v for k, v in c['fp_by_context'].items()}))
    with open(HERE / 'b_hp_fp_per_seed.tsv', 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        wr.writeheader()
        wr.writerows(rows)
    v4 = res['ver4 (ppr8 Illumina / PacBio HiFi, pp ONT)']
    S = dict(truth=T, fp_hp7=res, check_old_per_seed_rows=n_old, check_old_pooled_md=md6,
             ver4_range_seed_mean=[min(v4[p]['share_seed_mean'] for p in PLATS), max(v4[p]['share_seed_mean'] for p in PLATS)],
             ver4_range_pooled=[min(v4[p]['share_pooled'] for p in PLATS), max(v4[p]['share_pooled'] for p in PLATS)],
             ver4_range_per_seed=[min(x for p in PLATS for x in v4[p]['share']), max(x for p in PLATS for x in v4[p]['share'])])
    (HERE / 'b_hp_fp.json').write_text(json.dumps(S, indent=1) + '\n')

    pc = lambda x: f'{100 * x:.1f}%'  # noqa: E731
    L = ['# ver4 recount (b): HG008T chr1 INDEL false calls in long homopolymers, and population alleles in the truth', '',
         'Generated by `b_hp_fp.py`. Six-tumor release INDEL models `release_s{0,1,2}_<P>_INDEL` on HPRC v2.1, HG008T chr1 '
         '(held out), benchmark regions, PASS calls after the population filter (`eval_bed/vcfeval/pon_calls`). A false call is '
         '"in a homopolymer >= 7 bp" when its own inserted / deleted bases are one repeated base in a GRCh38 run of >= 7 bp '
         '(`truth_labels.RepeatContext` class HP>=7, as E17a).', '',
         '## 1. Share of false-positive calls in homopolymers >= 7 bp, HG008T chr1 after the filter', '',
         '| PASS rules | platform | rule | FP in HP>=7 / all FP, seed 0 | seed 1 | seed 2 | share seed 0 / 1 / 2 | seed mean of the share | pooled share (FP summed over seeds) |',
         '|---|---|---|---|---|---|---|---|---|']
    for label, d in res.items():
        for p in PLATS:
            x = d[p]
            L.append(f"| {label} | {PNAME[p]} | {x['rule']} | " + ' | '.join(f"{a} / {b}" for a, b in zip(x['fp_hp7'], x['fp'])) +
                     f" | {' / '.join(pc(s) for s in x['share'])} | {pc(x['share_seed_mean'])} | {pc(x['share_pooled'])} |")
    o = res['letter ver3 (E17a final_models: ppr8 on all platforms)']
    L += ['', f"Letter ver3's \"44–70%\" = the pooled shares of the first three rows ({' / '.join(pc(o[p]['share_pooled']) for p in PLATS)}; "
          'E17a `final_models.md` §6, reproduced). ver4: seed means ' + ' / '.join(pc(v4[p]['share_seed_mean']) for p in PLATS) +
          ' (Illumina / PacBio HiFi / ONT), pooled ' + ' / '.join(pc(v4[p]['share_pooled']) for p in PLATS) +
          f"; range of the seed means {pc(S['ver4_range_seed_mean'][0])}–{pc(S['ver4_range_seed_mean'][1])}, of the single seeds "
          f"{pc(S['ver4_range_per_seed'][0])}–{pc(S['ver4_range_per_seed'][1])}.", '',
          '### False calls by repeat context of the call, ver4 rules (FP per seed 0 / 1 / 2)', '',
          '| platform | rule | ' + ' | '.join(E.CONTEXTS) + ' | all |', '|---|---|' + '---|' * (len(E.CONTEXTS) + 1)]
    for p in PLATS:
        x = v4[p]
        L.append(f"| {PNAME[p]} | {x['rule']} | " + ' | '.join(' / '.join(str(c[k]) for c in x['fp_by_context']) for k in E.CONTEXTS) +
                 f" | {' / '.join(str(v) for v in x['fp'])} |")
    hp7, hp7t = T['in_bed_context'].get('HP>=7', 0), T['in_bed_tagged_context'].get('HP>=7', 0)
    L += ['', '## 2. Population alleles among the HG008T chr1 truth INDELs (current INDEL filter rule)', '',
          '| quantity | value | source |', '|---|---|---|',
          f"| chr1 truth INDELs / in the benchmark regions | {T['chr1_truth']} / {T['in_bed']} | `truth.INDEL.pon.vcf.gz` records / `pon_calls` tp-baseline + fn (= report.json `baseline`) |",
          f"| population alleles (tagged), chr1 / in the benchmark regions | {T['chr1_tagged']} / **{T['in_bed_tagged']} of {T['in_bed']}** ({pc(T['in_bed_tagged'] / T['in_bed'])}; recall bound {T['recall_bound']:.3f}) | "
          'FILTER PanelOfNormals in `truth.INDEL.pon.vcf.gz`; re-derived here from gnomAD + CoLoRSdb with the rule (identical for every record) |',
          f"| of the tagged in-region truth INDELs, in homopolymers >= 7 bp | **{hp7t}** (of {hp7} truth INDELs in HP>=7) | `RepeatContext` class HP>=7 of the truth record |",
          f"| tagged in-region truth INDELs by context | {T['in_bed_tagged_context']} | same |",
          f"| in-region truth INDELs by context | {T['in_bed_context']} | same |",
          f"| tagged by gnomAD / by CoLoRSdb (in regions; an allele can be in both) | {T['in_bed_tagged_by_panel']['gnomAD']} / {T['in_bed_tagged_by_panel']['CoLoRSdb']} | `PANSOMA_PON` |",
          f"| tagged at the AF >= 0.01 floor (the ALT changes a homopolymer >= 8 bp by the filter's `homopolymer_run`) / at the AF >= 1e-4 floor | {T['in_bed_tagged_floor_001']} / {T['in_bed_tagged_floor_1e4']} | `b_truth_indel_pon.tsv` `af_floor` |",
          '', f"Rule constants read from `scripts/filter_panel_of_normals.py` (asserted): INDEL_HOMOPOLYMER = {F.INDEL_HOMOPOLYMER}, "
          f"INDEL_MIN_AF_HOMOPOLYMER = {F.INDEL_MIN_AF_HOMOPOLYMER}, INDEL_MIN_AF = {F.INDEL_MIN_AF}; INDEL panels gnomAD and CoLoRSdb only. "
          f"Panels: `{T['panels'][0]}`, `{T['panels'][1]}`; FASTA `{T['fasta']}`.", '',
          '## 3. Checks', '',
          f"- Old rule set (ppr8 on all three platforms): TP / FN / FP, tagged truth, FP in HP>=7, truth in HP>=7 and its tagged count equal E17a `final_models.tsv` for all {n_old} release HG008T call sets; "
          f"the pooled shares equal `final_models.md` §6 ({', '.join(f'{p} {v}' for p, v in md6.items())}).",
          f"- Every call set used ({len(cnt)}): TP / FN / TP-call / FP equal report.json `pon.at_threshold`; truth = {T['in_bed']} = `pon.baseline` = `raw.baseline`; "
          f"truth, tagged truth ({T['in_bed_tagged']}), truth in HP>=7 ({hp7}) and its tagged count ({hp7t}) are identical across platforms, seeds and rules.",
          f"- Population tags: for each of the {T['chr1_truth']} chr1 truth INDELs the tag re-derived here (own exact-allele AF lookup in the two panels + the rule's floors) equals "
          'the tag from the filter script\'s own `record_matches` / `indel_min_af` and the FILTER / PANSOMA_PON of `truth.INDEL.pon.vcf.gz` '
          f"(report.json `truth_pon_tags`: tagged {T['report_truth_pon_tags']['tagged']} of {T['report_truth_pon_tags']['records']}).",
          '- PASS rule of each ver4 call set (report.json `steps.graph_vcf.call`): ' +
          '; '.join(f"{PNAME[p]} ({v4[p]['rule']}): {' | '.join(sorted({c.split(' (', 1)[1][:-1] for c in v4[p]['call']}))}" for p in PLATS) + '.',
          f"- Multi-allelic truth records: {T['multi_alt_records']} on chr1, {T['multi_alt_in_bed']} in the regions; a record is keyed by its first ALT (as E17a) "
          'and tagged when any of its ALTs passes its own floor (as the filter).']
    (HERE / 'b_hp_fp.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
