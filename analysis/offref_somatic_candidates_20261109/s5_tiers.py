"""s5: tiers of off-reference somatic candidates of one tumor, from the s3 assembly classes and the s4 read counts.

usage: s5_tiers.py TUMOR

Unit = key (candidate_id; the d46 node-frame allele, the same in the three platform sets of a tumor). Rows of a key are
its tensors in the three sets (SNV / INDEL tables).

Read rules per row and sample (counts of DATA/s4/<TUMOR>/counts/<LABEL>.*.tsv at the row's locus):
  context   hp (homopolymer >= 6 at the site, in alt_q or the REF context), str (period 2-6 run with >= 3 copies), other
  normal    other context: present  alt >= 2 and alt / (alt + ref) >= 0.05
                           absent   ref >= MIN_REF[platform], alt <= 1 and alt / (alt + ref) < 0.03
            hp / str (SNV or INDEL; a graph-frame SNV in a run is often a length change in the patient frame): only
            Illumina reads (long-read normals show the error background), present alt >= 3 and fraction >= 0.15,
            absent ref >= 10 and fraction < 0.10 (calibrated on the controls: Illumina normal fraction >= 0.15 in 98-99 %
            of germline and 0.6-3.6 % of somatic controls in runs; somatic p95 0.06-0.10)
            (ref = the patient-REF k-mer: graph REF for ref_exact rows, the normal assembly's allele for site_only)
            row normal = present if an informative normal is present, absent if one is absent and none present,
            else no_depth; no locus = no_locus
  tumor     support on a platform: alt >= 3 and alt / (alt + ref) >= 0.05 (hp / str: Illumina only, fraction >= 0.15)
Key features (over its rows): sets (platforms with a tensor), best assembly evidence (tumor_only with ref_exact /
site_only > tumor_only other > neither with ref_exact / site_only), normal (present over any row wins, then absent),
evidence platforms = platforms with read support at a row locus plus platforms with a tensor of the key (a tensor
needs >= 8 % of that platform's reads with the ALT); for hp / str rows only Illumina counts, so they reach T2 at most;
clean = some row with n_a1 >= 4, n_split == 0, multi_copy 0 and a specific k-mer.
Tiers (first match; T1-T3 need tumor read support at a row locus on >= 1 platform, which shows the GRCh38 locus is
where the ALT reads are, so that the normal check there means something):
  rejected_normal  normal present
  (no tumor read support at any row locus: asm_only for tumor_only keys, else other)
  T1               tumor_only + ref_exact / site_only, normal absent, evidence on >= 2 platforms, clean
  T2               tumor_only + ref_exact / site_only, normal absent (one platform, or not clean)
  T3_reads_only    neither + ref_exact / site_only, normal absent, evidence on >= 2 platforms, clean
  asm_only         tumor_only, normal no_depth / no_locus or no read support (no read check possible)
  other
Writes DATA/s5/<TUMOR>/keys.tsv.gz, candidates.tsv (T1 / T2 / T3 off-reference keys) and summary.md.
"""
import csv
import glob
import gzip
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, context, platform  # noqa: E402

BAMS = '/scratch/jshen/Github/Pansoma/tmp/offref_somatic_candidates_20261109/bams.tsv'
MIN_REF = {'Illumina': 10, 'PacBio': 8, 'HiFi': 8, 'fiberseq': 8, 'ONT': 10}
GOOD = ('ref_exact', 'site_only')
TIERS = ('T1', 'T2', 'T3_reads_only', 'asm_only', 'rejected_normal', 'other')
PLATS = {'HG008T': ('Illumina', 'ONT', 'PacBio'), 'COLO829T': ('Illumina', 'ONT', 'fiberseq')}
PATTERNS = []   # filled in main() for the tumor


def normal_call(c, rep=False):
    alt, ref = c['alt'], c['ref']
    frac = alt / (alt + ref) if alt + ref else (1.0 if alt else 0.0)
    if rep:
        if alt >= 3 and frac >= 0.15:
            return 'present'
        if ref >= 10 and frac < 0.10:
            return 'absent'
        return 'low_depth'
    if alt >= 2 and frac >= 0.05:
        return 'present'
    if ref >= c['min_ref'] and alt <= 1 and frac < 0.03:
        return 'absent'
    return 'low_depth'


def tumor_support(c, rep=False):
    alt, ref = c['alt'], c['ref']
    return alt >= 3 and alt / (alt + ref) >= (0.15 if rep else 0.05)


def load_counts(tumor):
    labels = [line.rstrip('\n').split('\t')[1] for line in open(BAMS) if line.split('\t')[0] == tumor]
    counts = {}
    for lab in labels:
        files = sorted(glob.glob(f'{DATA}/s4/{tumor}/counts/{lab}.*.tsv'))
        if not files:
            print(f'warning: no counts for {lab}', file=sys.stderr)
            continue
        d = {}
        for p in files:
            for r in csv.DictReader(open(p), delimiter='\t'):
                d[r['locus_id']] = {k: int(v) for k, v in r.items() if k != 'locus_id'}
        counts[lab] = d
    return counts


def main():
    tumor = sys.argv[1]
    pl = PLATS[tumor]
    PATTERNS[:] = [','.join(sorted(c)) for c in ([pl[0]], [pl[1]], [pl[2]], [pl[0], pl[1]], [pl[0], pl[2]], [pl[1], pl[2]], pl)]
    counts = load_counts(tumor)
    loci = {r['locus_id']: r for r in csv.DictReader(open(f'{DATA}/s4/{tumor}/loci.tsv'), delimiter='\t')}
    rl = {}
    for r in csv.DictReader(gzip.open(f'{DATA}/s4/{tumor}/row_loci.tsv.gz', 'rt'), delimiter='\t'):
        rl[(r['set'], r['kind'], r['candidate_id'])] = r
    keys = defaultdict(list)
    funnel = Counter()
    for r in csv.DictReader(gzip.open(f'{DATA}/s3/{tumor}/rows.tsv.gz', 'rt'), delimiter='\t'):
        funnel[(r['set'], r['kind'], r['group'], r['asm_class'], r['normal_side'])] += 1
        x = rl.get((r['set'], r['kind'], r['candidate_id']))
        if x is None:
            continue
        row = {k: r[k] for k in ('set', 'kind', 'chrom', 'candidate_id', 'group', 'label', 'reason', 'af', 'n_a1',
                                 'n_split', 'asm_class', 'normal_side', 'multi_copy', 'anchor', 'anchor_mapped')}
        row['known'], row['locus'], row['anchor_src'] = x['known'], x['locus'], x['anchor_src']
        row['ctx'] = context(r) if r['asm_class'] not in ('no_query', 'short_flank') else ''
        lc = loci.get(x['locus'])
        row['specific'] = lc['specific'] if lc else ''
        row['ref_src'] = lc['ref_src'] if lc else ''
        rep = row['ctx'] in ('hp', 'str')
        ncalls, sup, rc_ = [], [], {}
        for lab, d in counts.items():
            c = d.get(x['locus']) if lc else None
            if c is None:
                continue
            p = platform(lab)
            rc_[lab] = f"{c['alt']}/{c['ref']}/{c['span']}"
            informative = not rep or p == 'Illumina'
            if lab.startswith('N_') and informative:
                ncalls.append(normal_call({**c, 'min_ref': MIN_REF[p]}, rep))
            if lab.startswith('T_') and informative and tumor_support(c, rep):
                sup.append(p)
        row['normal'] = ('no_locus' if not lc else 'present' if 'present' in ncalls else
                         'absent' if 'absent' in ncalls else 'no_depth')
        row['support'] = ','.join(sorted(set(sup)))
        row['inf_set'] = platform(r['set']) if (not rep or platform(r['set']) == 'Illumina') else ''
        row['reads'] = ';'.join(f'{k}={v}' for k, v in sorted(rc_.items()))
        keys[(row['group'], r['candidate_id'])].append(row)
    out = f'{DATA}/s5/{tumor}'
    os.makedirs(out, exist_ok=True)
    recs = []
    for (group, key), rows in keys.items():
        sets = sorted({platform(r['set']) for r in rows})
        tum_good = [r for r in rows if r['asm_class'] == 'tumor_only' and r['normal_side'] in GOOD]
        tum_any = [r for r in rows if r['asm_class'] == 'tumor_only']
        nei_good = [r for r in rows if r['asm_class'] == 'neither' and r['normal_side'] in GOOD]
        nc = [r['normal'] for r in rows]
        normal = 'present' if 'present' in nc else 'absent' if 'absent' in nc else (
            'no_depth' if 'no_depth' in nc else 'no_locus')
        support = set()
        for r in rows:
            support.update(x for x in r['support'].split(',') if x)
        evidence = support | {r['inf_set'] for r in rows if r['inf_set']}
        clean = any(r['n_a1'] and int(r['n_a1']) >= 4 and r['n_split'] == '0' and r['multi_copy'] == '0' and
                    r['specific'] == '1' for r in rows)
        if normal == 'present':
            tier = 'rejected_normal'
        elif not support:
            tier = 'asm_only' if tum_any else 'other'
        elif tum_good and normal == 'absent' and len(evidence) >= 2 and clean:
            tier = 'T1'
        elif tum_good and normal == 'absent':
            tier = 'T2'
        elif nei_good and normal == 'absent' and len(evidence) >= 2 and clean:
            tier = 'T3_reads_only'
        elif tum_any and normal in ('no_depth', 'no_locus'):
            tier = 'asm_only'
        else:
            tier = 'other'
        r0 = (tum_good or tum_any or nei_good or rows)[0]
        recs.append(dict(tumor=tumor, group=group, key=key, kind=r0['kind'], chrom=r0['chrom'], tier=tier,
                         sets=','.join(sets), support=','.join(sorted(support)), evidence=','.join(sorted(evidence)), normal=normal,
                         asm=';'.join(f"{platform(r['set'])}:{r['asm_class']}/{r['normal_side']}" for r in rows),
                         ctx=r0['ctx'], known=r0['known'], reason=r0['reason'], label=r0['label'], af=r0['af'],
                         anchor=r0['anchor'] or r0['anchor_mapped'], anchor_src=r0['anchor_src'],
                         specific=r0['specific'], ref_src=r0['ref_src'], clean=int(clean), reads=r0['reads']))
    cols = list(recs[0])
    with gzip.open(f'{out}/keys.tsv.gz', 'wt') as w:
        w.write('\t'.join(cols) + '\n')
        for x in recs:
            w.write('\t'.join(str(x[c]) for c in cols) + '\n')
    with open(f'{out}/candidates.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for x in sorted(recs, key=lambda x: (x['tier'], x['chrom'], x['key'])):
            if x['group'] == 'offref' and x['tier'] in ('T1', 'T2', 'T3_reads_only'):
                w.write('\t'.join(str(x[c]) for c in cols) + '\n')
    # events: off-reference keys of a tier and kind merged when their anchors lie within 10 bp (one event, several
    # node frames); keys without an anchor count as one event each
    ev = {}
    for t in ('T1', 'T2', 'T3_reads_only'):
        for k in ('SNV', 'INDEL'):
            xs = [x for x in recs if x['group'] == 'offref' and x['tier'] == t and x['kind'] == k]
            pos = sorted((x['anchor'].rsplit(':', 1)[0], int(x['anchor'].rsplit(':', 1)[1].split('-')[0]))
                         for x in xs if x['anchor'])
            n_ev = sum(1 for x in xs if not x['anchor'])
            for i, (c, p) in enumerate(pos):
                if i == 0 or pos[i - 1][0] != c or p - pos[i - 1][1] > 10:
                    n_ev += 1
            ev[(t, k)] = n_ev
    # summary
    with open(f'{out}/summary.md', 'w') as w:
        w.write(f'# {tumor} off-reference somatic candidates (keys)\n\n')
        w.write('| group | kind | ' + ' | '.join(TIERS) + ' | chr1 T1 | chr1 T2 |\n|' + '---|' * (len(TIERS) + 4) + '\n')
        tab = Counter((x['group'], x['kind'], x['tier']) for x in recs)
        c1 = Counter((x['group'], x['kind'], x['tier']) for x in recs if x['chrom'] == 'chr1')
        for g in sorted({x['group'] for x in recs}):
            for k in ('SNV', 'INDEL'):
                w.write(f'| {g} | {k} | ' + ' | '.join(str(tab[(g, k, t)]) for t in TIERS) +
                        f' | {c1[(g, k, "T1")]} | {c1[(g, k, "T2")]} |\n')
        w.write('\n## Off-reference candidates by platform sets with a tensor of the key\n\n| tier | kind | keys | events (anchors '
                'within 10 bp merged) | ' + ' | '.join(PATTERNS) + ' |\n|' + '---|' * (len(PATTERNS) + 4) + '\n')
        for t in ('T1', 'T2', 'T3_reads_only'):
            for k in ('SNV', 'INDEL'):
                xs = [x for x in recs if x['group'] == 'offref' and x['tier'] == t and x['kind'] == k]
                pc = Counter(x['sets'] for x in xs)
                w.write(f'| {t} | {k} | {len(xs)} | {ev[(t, k)]} | ' +
                        ' | '.join(str(pc[p]) for p in PATTERNS) + ' |\n')
        w.write('\n## Off-reference keys by label reason (in-BED keys without truth = FDR proxy; residual_partial = next '
                'to GRCh38 somatic truth)\n\n| reason | kind | keys | ' + ' | '.join(TIERS) + ' |\n|' +
                '---|' * (len(TIERS) + 3) + '\n')
        rt = Counter((x['reason'], x['kind'], x['tier']) for x in recs if x['group'] == 'offref')
        for reason in sorted({x['reason'] for x in recs if x['group'] == 'offref'}):
            for k in ('SNV', 'INDEL'):
                n = sum(rt[(reason, k, t)] for t in TIERS)
                if n:
                    w.write(f'| {reason} | {k} | {n} | ' + ' | '.join(str(rt[(reason, k, t)]) for t in TIERS) + ' |\n')
        w.write('\n## Known off-reference somatic keys\n\n| known | ' + ' | '.join(TIERS) + ' |\n|' +
                '---|' * (len(TIERS) + 1) + '\n')
        kt = Counter((x['known'], x['tier']) for x in recs if x['known'])
        for kn in sorted({x['known'] for x in recs if x['known']}):
            w.write(f'| {kn} | ' + ' | '.join(str(kt[(kn, t)]) for t in TIERS) + ' |\n')
        w.write('\n## Funnel (rows): set, kind, group, asm_class, normal_side\n\n| set | kind | group | asm_class | '
                'normal_side | rows |\n|---|---|---|---|---|---|\n')
        for k, v in sorted(funnel.items()):
            w.write('| ' + ' | '.join(k) + f' | {v} |\n')
    print(tumor, Counter((x['group'], x['tier']) for x in recs))


if __name__ == '__main__':
    main()
