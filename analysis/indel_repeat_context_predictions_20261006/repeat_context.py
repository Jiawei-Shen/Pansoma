"""Repeat context (homopolymer / STR / other) of the HG008T truth INDELs and of Pansoma's INDEL predictions (chr1 test).

Inputs
  truth   analysis/graph_absorbed_somatic_20261001 audit tables (data dir): b1_truth_context.tsv (every HG008T truth allele
          classed in GRCh38 with the b1 rule), b5_normal_frame.tsv (the same in the patient frame: the truth INFO
          HG008Nv62SOMATICVARIANT event on the HG008-N v6.2 event haplotype).
  calls   /scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_INDEL_*/ (the seven current HG008T Illumina INDEL runs; no
          PacBio / ONT INDEL model exists): vcf_chr1/HG008T_Illumina.INDEL.linear.vcf.gz (GRCh38 calls, FILTER PASS = p_somatic
          at the run's threshold), vcf_chr1/...linear.unplaced.vcf.gz (off-reference calls, no GRCh38 allele), and the rtg
          vcfeval outputs on the chr1 GIAB BED: vcf_chr1/eval_bed/vcfeval/raw_calls/{tp,fp}.vcf.gz (PASS calls before the
          PoN) and vcf_chr1_pon/eval_bed/vcfeval/pon_calls/{tp,fp}.vcf.gz (after the 2026-10-06 INDEL PoN: gnomAD / CoLoRSdb
          exact allele AF >= 0.01).
Repeat class (the b1 rule, GRCh38): unit = minimal period of the inserted / deleted bases after the shared prefix; tract =
  the longest exact period-unit stretch touching the event; HP>=7 (unit 1, tract >= 7) / HP4-6 / STR (unit 2-6, >= 3 copies)
  / VNTR>6 (unit > 6, >= 2 copies) / complex (REF and ALT both change) / none. Predictions are classed with the same
  function; truth classes are read from b1 (same function) and b5.
Outputs: tables.md (printed), calls_context.tsv (one row per evaluated PASS call: run, set raw/pon, tp/fp, class, unit,
  tract, length, p_somatic).
"""
import collections, csv, gzip
from pathlib import Path
import pysam
A = Path(__file__).resolve().parent
B = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/audit/biology')
R = Path('/scratch/jshen/data/pansoma_net_v2_runs')
RUNS = ['HG008_Illumina_INDEL_base', 'HG008_Illumina_INDEL_base_b1024', 'HG008_Illumina_INDEL_nopartial_b1024',
        'HG008_Illumina_INDEL_nopartial_b1024_nopc', 'HG008_Illumina_INDEL_nopartial_b1024_small',
        'HG008_Illumina_INDEL_nopartial_b1024_nopc_small', 'HG008_Illumina_INDEL_nopartial_b1024_nopc_small_lr1e4']
PON_DIR = {'HG008_Illumina_INDEL_base': 'vcf_chr1_r09_pon'}       # base: evaluated at its R0.9 threshold
fa = pysam.FastaFile('/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
CLASSES = ['HP>=7', 'STR2-6>=3copies', 'HP4-6', 'VNTR>6', 'complex', 'none']


def tract(seq, lo, hi, periods=range(1, 7)):               # b1_truth_context.py, verbatim
    best = (0, lo, lo)
    for p in periods:
        i = max(0, lo - 200)
        while i < min(len(seq) - p, hi + 200):
            if seq[i] != seq[i + p]: i += 1; continue
            j = i
            while j < len(seq) - p and seq[j] == seq[j + p]: j += 1
            a, b = i, j + p
            if a <= hi and b >= lo and b - a > best[2] - best[1] and b - a >= 2 * p: best = (p, a, b)
            i = j + 1
    return best


def minper(s):
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s: return p
    return len(s)


def context(c, pos0, ref, alt):                             # b1 context(), INDEL branch
    F = 300; s0 = max(0, pos0 - F); seq = fa.fetch(c, s0, pos0 + len(ref) + F).upper(); o = pos0 - s0
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]: k += 1
    r2, a2 = ref[k:], alt[k:]
    if r2 and a2:
        return dict(ctx_class='complex', indel_unit=0, ctx_tract=0, indel_len=len(alt) - len(ref))
    S = r2 or a2; u = minper(S); lo = o + k; hi = o + k + len(r2)
    if u <= 6:
        p, a, b = tract(seq, lo - 1 if not r2 else lo, hi + 1 if not r2 else hi, periods=[u])
    else:
        n = 0
        while seq[hi + n * u: hi + (n + 1) * u] == S: n += 1
        a, b = lo, hi + n * u
    L = b - a
    cls = ('HP>=7' if u == 1 and L >= 7 else 'HP4-6' if u == 1 and L >= 4 else 'STR2-6>=3copies' if 2 <= u <= 6 and L >= 3 * u
           else 'VNTR>6' if u > 6 and L >= 2 * u else 'none')
    return dict(ctx_class=cls, indel_unit=u, ctx_tract=L, indel_len=len(a2) - len(r2))


def vcf_records(path):
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '#': continue
            c, pos, _id, ref, alt, qual, filt, info = line.split('\t')[:8]
            yield c, int(pos) - 1, ref.upper(), alt.upper(), filt, dict(x.split('=', 1) for x in info.split(';') if '=' in x)


# truth
b1 = {r['truth_id']: r for r in csv.DictReader(open(B / 'b1_truth_context.tsv'), delimiter='\t')}
b5 = {r['truth_id']: r for r in csv.DictReader(open(B / 'b5_normal_frame.tsv'), delimiter='\t')}
auto = {f'chr{i}' for i in range(1, 23)}
tr = [t for t, r in b1.items() if r['kind2'] == 'INDEL' and r['chrom'] in auto]
tr1 = [t for t in tr if b1[t]['chrom'] == 'chr1' and b1[t]['in_bed'] == 'True']
# self-check: the copied classifier reproduces b1 on the truth
bad = sum(context(b1[t]['chrom'], int(b1[t]['vcf_pos']) - 1, b1[t]['vcf_ref'].upper(), b1[t]['vcf_alt'].upper())['ctx_class'] != b1[t]['ctx_class'] for t in tr1)
assert bad == 0, bad


def pclass(t):
    n = b5.get(t)
    if not n or not n['nf_kind']: return 'no INFO event'
    c, k = n['nf_class'], n['nf_kind']
    if k == 'INDEL': return 'in_STR' if c.startswith('in_STR') else c
    return 'substitution in a repeat' if c not in ('none', 'complex', '') else 'substitution'


pct = lambda a, b: f'{a:,} ({100 * a / b:.1f}%)' if b else '0'
out = ['# Repeat context: HG008T truth INDELs vs Pansoma INDEL predictions (chr1 test)', '']
out += ['## 1. Truth INDELs', '', '| class | chr1-22 (GRCh38 frame) | chr1 BED (GRCh38 frame) | chr1-22 (patient frame) | chr1 BED (patient frame) |', '|---|---:|---:|---:|---:|']
g_all, g_1 = collections.Counter(b1[t]['ctx_class'] for t in tr), collections.Counter(b1[t]['ctx_class'] for t in tr1)
p_all, p_1 = collections.Counter(pclass(t) for t in tr), collections.Counter(pclass(t) for t in tr1)
PCL = ['HP>=7', 'STR2-6>=3copies', 'HP4-6', 'VNTR>6', 'in_STR', 'substitution in a repeat', 'substitution', 'none', 'no INFO event']
for k in CLASSES + ['in_STR', 'substitution in a repeat', 'substitution', 'no INFO event']:
    if g_all[k] or p_all[k]:
        out.append(f'| {k} | {pct(g_all[k], len(tr)) if k in CLASSES else ""} | {pct(g_1[k], len(tr1)) if k in CLASSES else ""} | {pct(p_all[k], len(tr))} | {pct(p_1[k], len(tr1))} |')
out.append(f'| total | {len(tr):,} | {len(tr1):,} | {len(tr):,} | {len(tr1):,} |')
out.append('')

# predictions
rows, T = [], {}
for run in RUNS:
    d = R / run
    for s, sub in (('raw', 'vcf_chr1/eval_bed/vcfeval/raw_calls'), ('pon', f'{PON_DIR.get(run, "vcf_chr1_pon")}/eval_bed/vcfeval/pon_calls')):
        for res in ('tp', 'fp'):
            p = d / sub / f'{res}.vcf.gz'
            if not p.exists(): continue
            for c, pos0, ref, alt, filt, info in vcf_records(p):
                x = context(c, pos0, ref, alt)
                rows.append(dict(run=run, set=s, result=res, chrom=c, pos=pos0 + 1, ref=ref, alt=alt, p_somatic=info.get('P_SOMATIC', ''), **x))
    T[run] = dict(raw_pass=sum(1 for *_, f, _ in vcf_records(d / 'vcf_chr1/HG008T_Illumina.INDEL.linear.vcf.gz') if f == 'PASS'),
                  unplaced_pass=sum(1 for *_, f, _ in vcf_records(d / 'vcf_chr1/HG008T_Illumina.INDEL.linear.unplaced.vcf.gz') if f == 'PASS'))
with open(A / 'calls_context.tsv', 'w') as w:
    cols = list(rows[0]); w.write('\t'.join(cols) + '\n')
    for r in rows: w.write('\t'.join(str(r[c]) for c in cols) + '\n')

out += ['## 2. Pansoma PASS INDEL calls on chr1 (inside the GIAB BED, GRCh38-placed; rtg vcfeval TP / FP), by class', '',
        'raw = PASS at the run threshold before the PoN; pon = after the INDEL PoN (gnomAD / CoLoRSdb AF >= 0.01). '
        'Truth chr1 BED for comparison: ' + ', '.join(f'{k} {g_1[k]}' for k in CLASSES if g_1[k]) + f' (total {len(tr1)}).', '']
for run in RUNS:
    rr = [r for r in rows if r['run'] == run]
    if not rr: continue
    out += [f'### {run}', '', f"PASS calls on chr1: {T[run]['raw_pass']:,} placed + {T[run]['unplaced_pass']:,} off-reference (no GRCh38 allele, not classed).", '',
            '| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |',
            '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for k in CLASSES:
        c = {(s, res): sum(1 for r in rr if r['set'] == s and r['result'] == res and r['ctx_class'] == k) for s in ('raw', 'pon') for res in ('tp', 'fp')}
        n_raw, n_pon = c[('raw', 'tp')] + c[('raw', 'fp')], c[('pon', 'tp')] + c[('pon', 'fp')]
        if not (n_raw or n_pon or g_1[k]): continue
        pr = lambda tp, n: f'{tp / n:.3f}' if n else ''
        rc = lambda tp: f'{tp / g_1[k]:.3f}' if g_1[k] else ''
        out.append(f"| {k} | {pct(n_raw, sum(1 for r in rr if r['set'] == 'raw'))} | {c[('raw', 'tp')]} | {c[('raw', 'fp')]} | {pr(c[('raw', 'tp')], n_raw)} | "
                   f"{pct(n_pon, sum(1 for r in rr if r['set'] == 'pon'))} | {c[('pon', 'tp')]} | {c[('pon', 'fp')]} | {pr(c[('pon', 'tp')], n_pon)} | {g_1[k]} | {rc(c[('raw', 'tp')])} | {rc(c[('pon', 'tp')])} |")
    n_raw, n_pon = sum(r['set'] == 'raw' for r in rr), sum(r['set'] == 'pon' for r in rr)
    tp_raw, tp_pon = sum(r['set'] == 'raw' and r['result'] == 'tp' for r in rr), sum(r['set'] == 'pon' and r['result'] == 'tp' for r in rr)
    out.append(f'| total | {n_raw:,} | {tp_raw} | {n_raw - tp_raw} | {tp_raw / n_raw:.3f} | {n_pon:,} | {tp_pon} | {n_pon - tp_pon} | {tp_pon / max(1, n_pon):.3f} | {len(tr1)} | {tp_raw / len(tr1):.3f} | {tp_pon / len(tr1):.3f} |')
    # HP tract bins for FP vs TP (raw)
    bins = lambda L: '7-9' if L < 10 else '10-14' if L < 15 else '15-19' if L < 20 else '20-29' if L < 30 else '>=30'
    hp = [r for r in rr if r['set'] == 'raw' and r['ctx_class'] == 'HP>=7']
    tb = collections.Counter((r['result'], bins(int(r['ctx_tract']))) for r in hp)
    tt = collections.Counter(bins(int(b1[t]['ctx_tract'])) for t in tr1 if b1[t]['ctx_class'] == 'HP>=7')
    out += ['', 'HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): ' +
            ', '.join(f"{b} {tb[('tp', b)]} / {tb[('fp', b)]} (truth {tt[b]})" for b in ('7-9', '10-14', '15-19', '20-29', '>=30')), '']
open(A / 'tables.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
