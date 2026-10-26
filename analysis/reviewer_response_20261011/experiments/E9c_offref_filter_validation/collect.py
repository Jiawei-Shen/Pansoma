"""E9c / E9d tables from DATA/runs/*/calls.tsv, DATA/gold/*/calls.tsv and DATA/asm/<sample>/hits.tsv.

usage: collect.py            -> EXP/tables.md, EXP/summary.json, EXP/per_call.tsv.gz

Graph PoN (population tag) of a call, rule R(flank, min_ac): a walk of the run's graph near the call spells the ALT
context (in-graph), or GRCh38 does, or >= min_ac of the 88 HPRC v1.1 haplotypes (+ CHM13) at the call's GRCh38 bubble
interval do; the ALT context is the graph-path context or the read-consensus context with `flank` bases each side
(A = adaptive). Main rule (fixed before looking at chr1): R(A, 1).

Assembly class (E9d) of the ALT context (read consensus when >= 40 bp, else graph path), minimap2 -x sr hits that cover
>= 90 % of the context with NM <= 1 and match the event exactly:
  HG008T   germline (in HG008-N v6.2) / tumor_only (in HG008-T v3.2, not in HG008-N) / unsupported (in neither, but the
           REF context is found) / unverifiable (neither the ALT nor the REF context is found); multi-copy flag when a
           context hits > 2 loci in one assembly (both haplotypes = 2)
  COLO829T the same with the COLO829BL donor-specific assembly (hapX / hapY) as normal and the COLO829T hifiasm haplotypes
           as tumor (2026-10-23; before: germline / not_in_normal / unverifiable, no tumor assembly)
Exclusive categories of an unplaced PASS call, in this order: population (graph PoN), germline, tumor_only,
unverifiable, unsupported (a sample without a tumor assembly: population, germline, not_in_normal, unverifiable).
"""
import csv
import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e9c_lib as L  # noqa: E402
from e9d_asm import qid  # noqa: E402

csv.field_size_limit(sys.maxsize)
E9A = '/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection'
FLANKS = ('5', '10', '15', 'A')
MIN_AC = (1, 2, 5, 9)
CATS = {'HG008T': ('population', 'germline', 'tumor_only', 'unverifiable', 'unsupported'),
        'COLO829T': ('population', 'germline', 'tumor_only', 'unverifiable', 'unsupported')}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def pop_tag(r, flank='A', min_ac=1, srcs='gr', ingraph=True, hprc=True):
    for s in srcs:
        k = f'{s}{flank}'
        if ingraph and num(r.get(f'ingraph_{k}')) > 0:
            return True
        if hprc and (num(r.get(f'hprc_{k}')) >= min_ac or num(r.get(f'hprc_{k}_grch38')) > 0):
            return True
    return False


def load_asm(sample):
    d = f'{L.DATA}/asm/{sample}'
    if not os.path.exists(f'{d}/hits.tsv'):
        return None
    hits = defaultdict(lambda: defaultdict(int))   # qid -> role -> max loci in one target
    for r in csv.DictReader(open(f'{d}/hits.tsv'), delimiter='\t'):
        hits[r['qid']][r['role']] = max(hits[r['qid']][r['role']], int(r['n_loci']))
    mapped = {r['qid'] for r in csv.DictReader(open(f'{d}/queries.tsv'), delimiter='\t')}
    return hits, mapped


def asm_class(r, sample, asm):
    if asm is None or sample not in CATS:
        return 'not_checked', ''
    hits, mapped = asm
    alt = 'alt_read' if len(r.get('alt_read', '')) >= 40 else 'alt_graph'
    if not r.get(alt) or not r.get(f'{alt}_site'):
        return 'not_checked', ''
    qa = qid(r[alt], r[f'{alt}_site'])
    if qa not in mapped:
        return 'not_checked', ''
    ha = hits.get(qa, {})
    refs = [qid(r[c], r[f'{c}_site']) for c in ('ref_graph', 'ref_read') if r.get(c) and r.get(f'{c}_site')]
    href = [hits.get(q, {}) for q in refs]
    ref_any = any(h.get('normal', 0) + h.get('tumor', 0) > 0 for h in href)
    multi = max([ha.get('normal', 0), ha.get('tumor', 0)] + [h.get('normal', 0) for h in href] + [0]) > 2
    if ha.get('normal', 0) > 0:
        c = 'germline'
    elif sample not in L.TUMOR_ASM:
        c = 'not_in_normal' if ref_any else 'unverifiable'
    elif ha.get('tumor', 0) > 0:
        c = 'tumor_only'
    else:
        c = 'unsupported' if ref_any else 'unverifiable'
    return c, 'multi_copy' if multi else ''


def category(r, sample):
    if r['pop']:
        return 'population'
    return r['asm'] if r['asm'] != 'not_checked' else 'not_checked'


def run_meta(name, rj):
    parts = name.split('__')
    source = parts[0]
    model = parts[1] if source != 'rcpu' or parts[1] != 'xs' else '__'.join(parts[1:4])
    s = rj.get('tensors', '').rstrip('/').split('/')[-2] if rj.get('tensors') else ''
    platform = next((p for p in ('Illumina', 'PacBio', 'ONT', 'fiberseq') if p in s), '')
    purity = next((p for p in ('P75', 'P50', 'P25') if s.endswith(p)), 'P100')
    rtype = 'legacy' if source == 'legacy' else ('cross-sample' if '__xs__' in f'__{name}' or name.startswith('rcpu__xs')
                                                 else 'in-sample')
    return dict(source=source, set=s, platform=platform, purity=purity, rtype=rtype, model=model)


def bed_intervals(path):
    out = defaultdict(list)
    if not path or not os.path.exists(path):
        return None
    op = gzip.open if path.endswith('.gz') else open
    for line in op(path, 'rt'):
        f = line.split('\t')
        if len(f) >= 3 and not line.startswith(('#', 'track')):
            out[f[0]].append((int(f[1]), int(f[2])))
    for c in out:
        out[c].sort()
    return out


def inside(bed, chrom, lo, hi):
    import bisect
    iv = bed.get(chrom, [])
    k = bisect.bisect_right(iv, (lo, float('inf'))) - 1
    return k >= 0 and iv[k][0] <= lo and hi <= iv[k][1]


def prf(tp, fp, rec):
    p = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * p * rec / (p + rec) if p + rec else 0.0
    return round(p, 4), round(f1, 4)


def e9a_projection(name):
    parts = name.split('__')
    if parts[0] != 'rcpu':
        return None
    run, vdir = '__'.join(parts[1:-1]), parts[-1]
    for p in (f'{E9A}/runs/{run}/{vdir}/projected.tsv', f'{E9A}/{run}.projected.tsv', f'{E9A}/{run}/{vdir}/projected.tsv'):
        if os.path.exists(p):
            return {r['candidate_id'].split(';')[0]: r for r in csv.DictReader(open(p), delimiter='\t')}
    return None


def main():
    runs = []
    for p in sorted(glob.glob(f'{L.DATA}/runs/*/calls.tsv')):
        d = os.path.dirname(p)
        rj = json.load(open(f'{d}/run.json'))
        runs.append((os.path.basename(d), rj, list(csv.DictReader(open(p), delimiter='\t'))))
    golds = []
    for p in sorted(glob.glob(f'{L.DATA}/gold/*/calls.tsv')):
        d = os.path.dirname(p)
        golds.append((os.path.basename(d), json.load(open(f'{d}/run.json')), list(csv.DictReader(open(p), delimiter='\t'))))
    asm = {s: load_asm(s) for s in CATS}
    summary = dict(runs=[], gold=[], controls=[], grid=[], strata=[], cross_platform=[], gold_in_runs=[])
    per_call = []
    for name, rj, rows in runs:
        meta = run_meta(name, rj)
        sample = rj['sample']
        proj = e9a_projection(name)
        for r in rows:
            r['pop'] = pop_tag(r) if r.get('tensor') == 'True' else False
            r['asm'], r['asm_flag'] = asm_class(r, sample, asm.get(sample)) if r.get('tensor') == 'True' else ('no_tensor', '')
            r['cat'] = category(r, sample) if r.get('tensor') == 'True' else 'no_tensor'
            r['proj_class'] = (proj.get(r['candidate_id']) or {}).get('proj_class', '') if proj else ''
            per_call.append(dict(run=name, group=r['group'], candidate_id=r['candidate_id'], kind=r.get('kind_event', ''),
                                 af=r.get('af', ''), homopolymer=r.get('homopolymer', ''), str_unit=r.get('str_unit', ''),
                                 str_copies=r.get('str_copies', ''), locus_status=r.get('locus_status', ''),
                                 locus=f"{r.get('locus_chrom', '')}:{r.get('locus_lo', '')}-{r.get('locus_hi', '')}",
                                 pop=int(r['pop']), hprc_gA=r.get('hprc_gA', ''), hprc_rA=r.get('hprc_rA', ''),
                                 ingraph_gA=r.get('ingraph_gA', ''), grch38_gA=r.get('hprc_gA_grch38', ''),
                                 asm=r['asm'], asm_flag=r['asm_flag'], category=r['cat'], proj_class=r['proj_class']))
        up = [r for r in rows if r['group'] == 'unplaced_pass' and r.get('tensor') == 'True']
        U = rj.get('unplaced_pass', len(up))
        n_proc = len([r for r in rows if r['group'] == 'unplaced_pass'])
        cats = Counter(r['cat'] for r in up)
        frac = {c: cats[c] / len(up) if up else 0.0 for c in set(cats) | set(CATS.get(sample, ()))}
        rec = dict(name=name, **meta, graph=rj['graph'], kind=rj['kind'], sample=sample, unplaced_pass=U,
                   processed=n_proc, with_context=len(up), categories=dict(cats),
                   pop_ingraph=sum(pop_tag(r, hprc=False) for r in up), pop_hprc=sum(pop_tag(r, ingraph=False) for r in up),
                   multi_copy=sum(r['asm_flag'] == 'multi_copy' for r in up),
                   tumor_only_multi_copy=sum(r['asm_flag'] == 'multi_copy' and r['asm'] == 'tumor_only' for r in up))
        if proj:
            rec['e9a_proj_class_of_unplaced_pass'] = dict(Counter(r['proj_class'] for r in up))
            rec['e9a_category_by_class'] = {k: dict(Counter(r['cat'] for r in up if r['proj_class'] == k))
                                            for k in ('i', 'ii', 'iii')}
        # corrected precision (scaled to all unplaced PASS calls when they were sampled)
        for ev in ('eval_nobed', 'eval_bed'):
            rep = rj.get('report', {}).get(ev)
            if not rep:
                continue
            tp, fp, base = rep['pon']['tp_call'], rep['pon']['fp'], rep.get('baseline') or 0
            tpb = rep['pon']['tp_baseline']
            recall = tpb / base if base else 0.0
            if ev == 'eval_nobed':
                Ue = U
                f_ = frac
            else:
                bed = bed_intervals(rep.get('bed'))
                if bed is None:
                    continue
                inb = [r for r in up if r.get('locus_lo') in ('', None) or r.get('locus_status') not in ('bubble', 'one_side')
                       or inside(bed, r['locus_chrom'], int(r['locus_lo']), int(r['locus_hi']))]
                Ue = U * (len(inb) / len(up)) if up else 0
                c2 = Counter(r['cat'] for r in inb)
                f_ = {c: c2[c] / len(inb) if inb else 0.0 for c in set(c2) | set(frac)}
            rem = Ue * (1 - f_.get('population', 0))
            som = Ue * f_.get('tumor_only', 0)
            out = dict(tp=tp, fp=fp, unplaced_pass=round(Ue, 1), recall=round(recall, 4))
            out['placed_only'] = prf(tp, fp, recall)
            out['all_unplaced_fp'] = prf(tp, fp + Ue, recall)
            out['after_graph_pon'] = prf(tp, fp + rem, recall)
            out['after_graph_pon_tumor_only_tp'] = prf(tp + som, fp + rem - som, recall)
            out['remaining_unplaced_after_graph_pon'] = round(rem, 1)
            out['unplaced_share_of_pass'] = round(Ue / (tp + fp + Ue), 4) if tp + fp + Ue else 0.0
            out['unplaced_share_after_graph_pon'] = round(rem / (tp + fp + rem), 4) if tp + fp + rem else 0.0
            rec[ev] = out
        summary['runs'].append(rec)
        # controls and grid
        for g in ('ctrl_tp', 'ctrl_fp', 'ctrl_pon', 'unplaced_pass', 'gold'):
            rs = [r for r in rows if r['group'] == g and r.get('tensor') == 'True']
            if not rs:
                continue
            summary['controls'].append(dict(name=name, group=g, n=len(rs), pop=sum(r['pop'] for r in rs),
                                            asm=dict(Counter(r['asm'] for r in rs))))
            for fl in FLANKS:
                for ac in MIN_AC:
                    summary['grid'].append(dict(name=name, group=g, flank=fl, min_ac=ac, n=len(rs),
                                                tagged=sum(pop_tag(r, fl, ac) for r in rs)))
        # repeat strata (R1.3)
        for r in up:
            hp, sc = int(num(r.get('homopolymer'))), int(num(r.get('str_copies')))
            r['stratum'] = 'homopolymer>=7' if hp >= 7 else ('STR' if sc >= 3 else 'other')
        for s in ('homopolymer>=7', 'STR', 'other'):
            rs = [r for r in up if r['stratum'] == s]
            summary['strata'].append(dict(name=name, stratum=s, n=len(rs), categories=dict(Counter(r['cat'] for r in rs))))
        # gold keys of this run's chromosomes
        gs = [r for r in rows if r['group'] == 'gold']
        if gs:
            summary['gold_in_runs'].append(dict(name=name, gold_with_tensor=len(gs),
                                                pass_unplaced=sum(r['filter'] == 'PASS' for r in gs),
                                                pass_and_tagged=sum(r['filter'] == 'PASS' and r['pop'] for r in gs),
                                                tagged=sum(r['pop'] for r in gs),
                                                asm=dict(Counter(r['asm'] for r in gs))))
    # cross-platform reproducibility of unplaced PASS keys (in-sample and legacy runs, P100 only)
    keysets = {}
    for name, rj, rows in runs:
        meta = run_meta(name, rj)
        if meta['rtype'] == 'cross-sample' or meta['purity'] != 'P100':
            continue
        keys = {f[2].split(';')[0] for f in L.vcf_records(rj['unplaced_vcf']) if f[6] == 'PASS'}
        keysets[name] = (meta, rj, keys)
    groups = defaultdict(list)
    for name, (meta, rj, keys) in keysets.items():
        groups[(meta['source'], rj['graph'], rj['sample'], rj['kind'])].append(name)
    by_name = {name: rows for name, _, rows in runs}
    for gk, names in groups.items():
        if len(names) < 2:
            continue
        for name in names:
            others = [keysets[o][2] for o in names if o != name]
            up = [r for r in by_name[name] if r['group'] == 'unplaced_pass' and r.get('tensor') == 'True']
            tab = defaultdict(Counter)
            for r in up:
                tab[r['cat']][sum(r['candidate_id'] in o for o in others)] += 1
            allk = keysets[name][2]
            summary['cross_platform'].append(dict(group='/'.join(gk), name=name, platforms=len(names),
                                                  unplaced_pass=len(allk),
                                                  in_another=sum(any(k in o for o in others) for k in allk),
                                                  by_category={c: dict(v) for c, v in tab.items()}))
    # gold calibration (all chromosomes; chr2-22 = calibration, chr1 = held out)
    gpon = {}
    gp = f'{L.DATA}/gold_grch38_pon/events.tsv'
    if os.path.exists(gp):
        gpon = {r['event']: r['grch38_pon_tagged'] == '1' for r in csv.DictReader(open(gp), delimiter='\t')}
    for name, rj, rows in golds:
        rs = [r for r in rows if r.get('tensor') == 'True']
        for r in rs:
            r['asm'], r['asm_flag'] = asm_class(r, 'HG008T', asm.get('HG008T'))
            per_call.append(dict(run=f'gold__{name}', group='gold', candidate_id=r['candidate_id'], kind=r.get('kind_event', ''),
                                 af=r.get('af', ''), homopolymer=r.get('homopolymer', ''), str_unit=r.get('str_unit', ''),
                                 str_copies=r.get('str_copies', ''), locus_status=r.get('locus_status', ''),
                                 locus=f"{r.get('locus_chrom', '')}:{r.get('locus_lo', '')}-{r.get('locus_hi', '')}",
                                 pop=int(pop_tag(r)), hprc_gA=r.get('hprc_gA', ''), hprc_rA=r.get('hprc_rA', ''),
                                 ingraph_gA=r.get('ingraph_gA', ''), grch38_gA=r.get('hprc_gA_grch38', ''),
                                 asm=r['asm'], asm_flag=r['asm_flag'], category='', proj_class=''))
        for part, sel in (('chr2-22', lambda r: r['chrom'] != 'chr1'), ('chr1', lambda r: r['chrom'] == 'chr1'),
                          ('all', lambda r: True)):
            ss = [r for r in rs if sel(r)]
            for fl in FLANKS:
                for ac in MIN_AC:
                    summary['gold'].append(dict(set=name, part=part, flank=fl, min_ac=ac, n=len(ss),
                                                tagged=sum(pop_tag(r, fl, ac) for r in ss),
                                                tagged_ingraph=sum(pop_tag(r, fl, ac, hprc=False) for r in ss)))
            both = Counter((int(pop_tag(r)), int(gpon.get(r['gold_event'], False))) for r in ss)
            summary['gold'].append(dict(set=name, part=part, flank='pon_vs_graph', n=len(ss),
                                        graph_pon=sum(v for (g, _), v in both.items() if g),
                                        grch38_pon=sum(v for (_, p), v in both.items() if p),
                                        both=both.get((1, 1), 0), neither=both.get((0, 0), 0)))
            tonly = [r for r in ss if r['asm'] == 'tumor_only']
            summary['gold'].append(dict(set=name, part=part, flank='tumor_only_gold', n=len(tonly),
                                        graph_pon=sum(pop_tag(r) for r in tonly),
                                        grch38_pon=sum(gpon.get(r['gold_event'], False) for r in tonly)))
            summary['gold'].append(dict(set=name, part=part, flank='asm', n=len(ss), asm=dict(Counter(r['asm'] for r in ss)),
                                        repeat=dict(Counter('homopolymer>=7' if num(r.get('homopolymer')) >= 7 else
                                                            ('STR' if num(r.get('str_copies')) >= 3 else 'other') for r in ss))))
    json.dump(summary, open(f'{L.EXP}/summary.json', 'w'), indent=1)
    cols = list(per_call[0]) if per_call else []
    with gzip.open(f'{L.EXP}/per_call.tsv.gz', 'wt') as w:
        w.write('\t'.join(cols) + '\n')
        for r in per_call:
            w.write('\t'.join(str(r[c]) for c in cols) + '\n')
    write_tables(summary)
    print(f'{len(runs)} runs, {len(golds)} gold sets -> {L.EXP}/tables.md, summary.json, per_call.tsv.gz')


def pct(a, b):
    return f'{a} ({100 * a / b:.1f} %)' if b else '-'


def write_tables(s):
    out = ['# E9c / E9d tables (generated by collect.py; do not edit)', '']
    out += ['## 1. Unplaced PASS calls per run: graph PoN and assembly classes (chr1)', '',
            'Fractions of the unplaced PASS calls with a tensor context (`with ctx`); `processed` < `unplaced PASS` '
            'when a random sample was processed. Categories are exclusive in the order population > germline > '
            'tumor-only > unverifiable > unsupported.', '',
            '| run | type | graph | kind | sample | platform | purity | unplaced PASS | processed | with ctx | population '
            '(graph PoN) | in-graph | HPRC hap | germline (normal asm) | tumor-only | not in normal | unverifiable | '
            'unsupported | not checked |',
            '|' + '---|' * 19]
    for r in s['runs']:
        c, n = r['categories'], r['with_context']
        out.append(f"| {r['name']} | {r['rtype']} | {r['graph']} | {r['kind']} | {r['sample']} | {r['platform']} | "
                   f"{r['purity']} | {r['unplaced_pass']} | {r['processed']} | {n} | {pct(c.get('population', 0), n)} | "
                   f"{pct(r['pop_ingraph'], n)} | {pct(r['pop_hprc'], n)} | {pct(c.get('germline', 0), n)} | "
                   f"{pct(c.get('tumor_only', 0), n)} | {pct(c.get('not_in_normal', 0), n)} | "
                   f"{pct(c.get('unverifiable', 0), n)} | {pct(c.get('unsupported', 0), n)} | {pct(c.get('not_checked', 0), n)} |")
    out += ['', '## 2. Corrected precision, chr1, after the GRCh38 PoN, at the run\'s PASS threshold', '',
            '`placed`: rtg on placed calls (the run\'s report); `all FP`: every unplaced PASS call as a false call; '
            '`graph PoN`: unplaced calls tagged by the graph PoN removed, the rest false; `+tumor-only TP`: the '
            'assembly-supported (tumor-only) remainder counted as true. Recall is rtg\'s (truth is on GRCh38). In the '
            'BED, an unplaced call counts when its GRCh38 bubble interval lies inside the BED (calls without an '
            'interval count as inside).', '',
            '| run | eval | TP | FP | unplaced PASS | unplaced share of PASS | recall | P placed | F1 placed | P all FP | '
            'F1 all FP | P graph PoN | F1 graph PoN | P +tumor-only | F1 +tumor-only | unplaced left | unplaced share left |',
            '|' + '---|' * 17]
    for r in s['runs']:
        for ev in ('eval_bed', 'eval_nobed'):
            e = r.get(ev)
            if not e:
                continue
            out.append(f"| {r['name']} | {ev[5:]} | {e['tp']} | {e['fp']} | {e['unplaced_pass']} | "
                       f"{e['unplaced_share_of_pass']} | {e['recall']} | "
                       f"{e['placed_only'][0]} | {e['placed_only'][1]} | {e['all_unplaced_fp'][0]} | {e['all_unplaced_fp'][1]} | "
                       f"{e['after_graph_pon'][0]} | {e['after_graph_pon'][1]} | {e['after_graph_pon_tumor_only_tp'][0]} | "
                       f"{e['after_graph_pon_tumor_only_tp'][1]} | {e['remaining_unplaced_after_graph_pon']} | "
                       f"{e['unplaced_share_after_graph_pon']} |")
    out += ['', '## 3. Controls: graph PoN tag rate (main rule) and assembly classes', '',
            '`ctrl_tp` placed PASS calls that rtg matched to the somatic truth (after the PoN); `ctrl_fp` placed false '
            'calls after the PoN; `ctrl_pon` placed PASS calls the GRCh38 PoN tags (population alleles); random '
            'samples of up to 150 each.', '',
            '| run | group | n | graph PoN tagged | assembly classes |', '|---|---|---|---|---|']
    for r in s['controls']:
        out.append(f"| {r['name']} | {r['group']} | {r['n']} | {pct(r['pop'], r['n'])} | "
                   + ', '.join(f'{k} {v}' for k, v in sorted(r['asm'].items())) + ' |')
    out += ['', '## 4. Gold set: assembly-verified off-reference somatic keys (all chromosomes)', '',
            'False-tag rate of the graph PoN by flank and minimum number of HPRC haplotypes (GRCh38 and in-graph '
            'walks always count). chr2-22 is the calibration part; chr1 is held out.', '',
            '| set | part | flank | min AC | n | tagged | of which in-graph |', '|---|---|---|---|---|---|---|']
    for r in s['gold']:
        if r['flank'] in ('asm', 'pon_vs_graph', 'tumor_only_gold'):
            continue
        out.append(f"| {r['set']} | {r['part']} | {r['flank']} | {r['min_ac']} | {r['n']} | {pct(r['tagged'], r['n'])} | "
                   f"{r['tagged_ingraph']} |")
    out += ['', '| set | part | n | assembly classes | repeat context |', '|---|---|---|---|---|']
    for r in s['gold']:
        if r['flank'] == 'asm':
            out.append(f"| {r['set']} | {r['part']} | {r['n']} | " + ', '.join(f'{k} {v}' for k, v in sorted(r['asm'].items()))
                       + ' | ' + ', '.join(f'{k} {v}' for k, v in sorted(r['repeat'].items())) + ' |')
    out += ['', 'Graph PoN (main rule) against the GRCh38 PoN (repository rule) applied to the same events\' GRCh38 truth '
            'records; and on the gold keys whose ALT context is tumor-only in the assemblies (clean positives).', '',
            '| set | part | n | graph PoN | GRCh38 PoN | both | neither | tumor-only gold n | graph PoN | GRCh38 PoN |',
            '|---|---|---|---|---|---|---|---|---|---|']
    pv = {(r['set'], r['part']): r for r in s['gold'] if r['flank'] == 'pon_vs_graph'}
    tg = {(r['set'], r['part']): r for r in s['gold'] if r['flank'] == 'tumor_only_gold'}
    for k, r in pv.items():
        t = tg.get(k, dict(n=0, graph_pon=0, grch38_pon=0))
        out.append(f"| {k[0]} | {k[1]} | {r['n']} | {pct(r['graph_pon'], r['n'])} | {pct(r['grch38_pon'], r['n'])} | "
                   f"{r['both']} | {r['neither']} | {t['n']} | {pct(t['graph_pon'], t['n'])} | {pct(t['grch38_pon'], t['n'])} |")
    out += ['', '## 5. Gold keys on chr1 in each run', '',
            '| run | gold keys with a tensor | PASS (unplaced) | tagged | PASS and tagged | assembly classes |',
            '|---|---|---|---|---|---|']
    for r in s['gold_in_runs']:
        out.append(f"| {r['name']} | {r['gold_with_tensor']} | {r['pass_unplaced']} | {r['tagged']} | {r['pass_and_tagged']} | "
                   + ', '.join(f'{k} {v}' for k, v in sorted(r['asm'].items())) + ' |')
    out += ['', '## 6. Rule grid on unplaced PASS calls and controls (tagged / n)', '',
            '| run | group | ' + ' | '.join(f'f{f} AC>={a}' for f in FLANKS for a in MIN_AC) + ' |',
            '|---|---|' + '---|' * (len(FLANKS) * len(MIN_AC))]
    grid = defaultdict(dict)
    for r in s['grid']:
        grid[(r['name'], r['group'])][(r['flank'], r['min_ac'])] = (r['tagged'], r['n'])
    for (name, g), v in grid.items():
        out.append(f'| {name} | {g} | ' + ' | '.join(f"{v[(f, a)][0]}/{v[(f, a)][1]}" for f in FLANKS for a in MIN_AC) + ' |')
    out += ['', '## 7. Repeat context of unplaced PASS calls (R1.3)', '',
            '| run | stratum | n | categories |', '|---|---|---|---|']
    for r in s['strata']:
        out.append(f"| {r['name']} | {r['stratum']} | {r['n']} | " + ', '.join(f'{k} {v}' for k, v in sorted(r['categories'].items())) + ' |')
    out += ['', '## 8. Cross-platform reproducibility of unplaced PASS keys (same node:offset:allele, same graph)', '',
            '| group | run | platforms | unplaced PASS keys | also PASS on another platform | by category: {category: {n other platforms: calls}} |',
            '|---|---|---|---|---|---|']
    for r in s['cross_platform']:
        out.append(f"| {r['group']} | {r['name']} | {r['platforms']} | {r['unplaced_pass']} | "
                   f"{pct(r['in_another'], r['unplaced_pass'])} | {json.dumps(r['by_category'], sort_keys=True)} |")
    open(f'{L.EXP}/tables.md', 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
