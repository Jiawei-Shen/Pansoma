"""Independent re-check (synthesis agent): v1 e053 VCF vs v2 runs on the 697 chr1 SNP truths.
Read-only on all inputs; writes only to this scratch dir (pon cache)."""
import gzip, json, bisect, collections, os, pickle, sys, time
import numpy as np
import pysam

T0 = time.time()
OUT = os.path.dirname(os.path.abspath(__file__))
V2 = '/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors'
DV = '/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/'
TRUTH_ALL = DV + 'HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz'
SBED = DV + 'HG008-T_somatic_smvar_benchmark_v0.2_all.bed'
GBED = '/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed'
GVCF = '/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz'
V1D = '/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/'
V1N = '/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC_HG008N_added/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/'
RUNS = '/scratch/jshen/data/pansoma_net_v2_runs/'
PON = '/scratch/jshen/data/Panels_of_Normal/'


def bed(path, chrom='chr1'):
    iv = []
    for l in open(path):
        if l.startswith(('#', 'track', 'browser')):
            continue
        f = l.split('\t')
        if f[0] == chrom:
            iv.append((int(f[1]), int(f[2])))
    iv.sort(); m = []
    for s, e in iv:
        if m and s <= m[-1][1]:
            m[-1] = (m[-1][0], max(m[-1][1], e))
        else:
            m.append((s, e))
    return m


def inter(a, b):
    i = j = 0; o = []
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if s < e:
            o.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return o


SB = bed(SBED); CONF = inter(SB, bed(GBED)); CST = [s for s, _ in CONF]


def inconf(p0):
    k = bisect.bisect_right(CST, p0) - 1
    return k >= 0 and CONF[k][0] <= p0 < CONF[k][1]


print('conf bp chr1', sum(e - s for s, e in CONF))
# 697 truths
pos697 = {}
with open(V2 + '/somatic.recall.tsv') as f:
    next(f)
    for l in f:
        t, c, p, k, ps, ib, st = l.rstrip('\n').split('\t')[:7]
        if c == 'chr1' and k == 'SNP' and ps == 'True' and ib == 'True':
            pos697[int(p)] = (t, st)
TRUTH = {}  # (pos0, ref, alt) -> tid
for r in pysam.VariantFile(TRUTH_ALL).fetch('chr1'):
    if r.pos in pos697 and len(r.ref) == 1 and all(len(a) == 1 for a in r.alts):
        for a in r.alts:
            TRUTH[(r.pos - 1, r.ref.upper(), a.upper())] = pos697[r.pos][0]
NT = len(set(TRUTH.values()))
print('697 check:', len(pos697), 'allele keys', len(TRUTH), 'distinct tids', NT,
      collections.Counter(s for _, s in pos697.values()))
# germline chr1 SNP alleles
GERM = set()
_seen1 = False
for r in pysam.VariantFile(GVCF):
    if r.chrom != 'chr1':
        if _seen1:
            break
        continue
    _seen1 = True
    if len(r.ref) == 1:
        for a in r.alts or ():
            if len(a) == 1:
                GERM.add((r.pos - 1, r.ref.upper(), a.upper()))
print('germline chr1 SNP alleles', len(GERM), 't=%.0f' % (time.time() - T0))


def read_vcf(path):
    d = {}
    for r in pysam.VariantFile(path).fetch('chr1'):
        if len(r.ref) != 1:
            continue
        for a in r.alts:
            if len(a) != 1:
                continue
            k = (r.pos - 1, r.ref.upper(), a.upper())
            d[k] = max(d.get(k, -1), float(r.info['PROB']))
    return d


def read_v2(run):
    pf = [x for x in os.listdir(RUNS + run + '/test_chr1') if x.endswith('SNV.predictions.ndjson.gz')][0]
    sites, offref, nt, neutral = {}, [], 0, set()
    with gzip.open(RUNS + run + '/test_chr1/' + pf, 'rt') as fp, open(V2 + '/SNV/chr1_labels.ndjson') as fl:
        for lp, ll in zip(fp, fl):
            p = json.loads(lp); lab = json.loads(ll)
            assert p['candidate_id'] == lab['candidate_id']
            if not p['in_test']:
                continue
            nt += 1
            g = lab['grch38']
            if g is None:
                offref.append(p['p_somatic']); continue
            k = (g['pos0'], g['ref'].upper(), g['alt'].upper())
            sites[k] = max(sites.get(k, -1), p['p_somatic'])
            if p['test_label'] == 1 and k not in TRUTH:
                neutral.add(k)
    return sites, np.array(offref), nt, neutral


def curve(sites, pon=None, offref=None, neutral=None):
    """TP distinct truths; FP distinct in-conf alleles not truth (PoN-flagged sites dropped). Returns arrays."""
    items = []
    for k, s in sites.items():
        if pon is not None and pon.get(k, False):
            continue
        if k in TRUTH:
            items.append((s, 1, TRUTH[k]))
        elif inconf(k[0]) and not (neutral and k in neutral):
            items.append((s, 0, None))
    if offref is not None:
        items += [(s, 0, None) for s in offref]
    items.sort(key=lambda x: -x[0])
    seen = set(); tp = fp = 0; pts = []
    for i, (s, t, tid) in enumerate(items):
        if t:
            if tid not in seen:
                seen.add(tid); tp += 1
        else:
            fp += 1
        if i + 1 == len(items) or items[i + 1][0] != s:
            pts.append((s, tp, fp))
    return pts


def summ(pts):
    P = np.array([tp / (tp + fp) if tp + fp else 1 for _, tp, fp in pts]); R = np.array([tp / NT for _, tp, _ in pts])
    F = np.where(P + R > 0, 2 * P * R / (P + R + 1e-12), 0)
    out = {}
    for x in (0.05, 0.10, 0.15, 0.20):
        m = P >= x; out['R@P%.2f' % x] = round(float(R[m].max()), 3) if m.any() else None
    for y in (0.5, 0.7, 0.8, 0.9):
        m = R >= y; out['P@R%.1f' % y] = round(float(P[m].max()), 3) if m.any() else None
    i = int(F.argmax()); out['maxF1'] = (round(float(F[i]), 3), round(float(P[i]), 3), round(float(R[i]), 3), round(float(pts[i][0]), 4))
    out['ceiling'] = round(float(R.max()), 3)
    return out


def at_tp(pts, tp_target):
    for s, tp, fp in pts:
        if tp >= tp_target:
            return s, tp, fp
    return None


V1 = {'v1_e053': read_vcf(V1D + 'pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.vcf.gz'),
      'v1_e053_ownPoN': read_vcf(V1D + 'pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.filtered_PoN.vcf.gz')}
V2R = {}
for run in ['HG008_Illumina_SNV_base', 'HG008_Illumina_SNV_keephard_w100', 'HG008_Illumina_SNV_scalars', 'HG008_Illumina_SNV_allnon_w100']:
    V2R[run] = read_v2(run)
    s, o, nt, neu = V2R[run]
    print(run, 'in_test', nt, 'on-ref sites', len(s), 'offref', len(o), 'label1-not-697 sites', len(neu),
          'truth sites present', len(set(TRUTH) & set(s)), 't=%.0f' % (time.time() - T0))
print('v1 truth present (PROB>=0.5)', len(set(TRUTH) & set(V1['v1_e053'])))

# ---------------- PoN lookup (repo rule: gnomAD allele, dbSNP allele, 1000G pos, CoLoRSdb pos) -------------
CACHE = OUT + '/my_pon.pkl'
cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
need = set(V1['v1_e053'])
for run, (s, o, nt, neu) in V2R.items():
    need |= {k for k, v in s.items() if v >= 0.02 and (k in TRUTH or inconf(k[0]))}
todo = sorted(need - set(cache))
print('PoN lookups needed', len(todo))
if todo:
    vf = [pysam.VariantFile(PON + f) for f in ['af-only-gnomad.hg38.vcf.gz', 'Homo_sapiens_assembly38.dbsnp138.vcf.gz',
                                                 '1000g_pon.hg38.vcf.gz', 'CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz']]
    mode = ['allele', 'allele', 'position', 'position']
    for n, (p0, ref, alt) in enumerate(todo):
        hit = []
        for v, m in zip(vf, mode):
            h = False
            for r in v.fetch('chr1', p0, p0 + 1):
                if r.pos != p0 + 1:
                    continue
                if m == 'position' or (r.ref.upper() == ref and alt in [a.upper() for a in (r.alts or ())]):
                    h = True; break
            hit.append(h)
        cache[(p0, ref, alt)] = tuple(hit)
    pickle.dump(cache, open(CACHE, 'wb'))
print('PoN done t=%.0f' % (time.time() - T0))
PONF = {k: any(v) for k, v in cache.items()}
# sites with p<0.02 not looked up -> treated as not PoN (they are never above the thresholds used for R<=0.85 anyway)

# agreement of my PoN rule vs v1's own filtered VCF
own = set(V1['v1_e053_ownPoN'])
agree = collections.Counter((k in own, not PONF[k]) for k in V1['v1_e053'])
print('v1 e053 calls: (kept by v1 own PoN, kept by my repo-rule PoN) ->', dict(agree))

print('\n==== curves (FP = distinct in-conf alleles, strict allele match; v2 label-1 non-697 SNV sites counted as FP) ====')
rows = {}
rows['v1_e053 noPoN'] = curve(V1['v1_e053'])
rows['v1_e053 ownPoN VCF'] = curve(V1['v1_e053_ownPoN'])
rows['v1_e053 myPoN'] = curve(V1['v1_e053'], pon=PONF)
for run, (s, o, nt, neu) in V2R.items():
    r = run.replace('HG008_Illumina_SNV_', 'v2_')
    rows[r + ' noPoN'] = curve(s)
    rows[r + ' noPoN +offref'] = curve(s, offref=o)
    rows[r + ' noPoN neutral-partial'] = curve(s, neutral=neu)
    rows[r + ' myPoN'] = curve(s, pon=PONF)
    rows[r + ' myPoN +offref'] = curve(s, pon=PONF, offref=o)
for k, pts in rows.items():
    print(k, summ(pts))

print('\n==== matched TP comparison ====')
for tgt_name, tgt in [('v1 e053 PROB>=0.5 noPoN', 'v1_e053 noPoN'), ('v1 e053 own PoN', 'v1_e053 ownPoN VCF')]:
    s, tp, fp = rows[tgt][-1]
    print(f'{tgt_name}: TP {tp} FP {fp} P {tp/(tp+fp):.4f} R {tp/NT:.4f}')
    for k in rows:
        if k.startswith('v2_') and (('noPoN' in k) == ('noPoN' in tgt) or ('myPoN' in k and 'PoN' in tgt)):
            if 'neutral' in k:
                continue
            a = at_tp(rows[k], tp)
            if a:
                print(f'   {k}: thr {a[0]:.4f} TP {a[1]} FP {a[2]} P {a[1]/(a[1]+a[2]):.4f}')

# FP composition at matched TP 608 (no PoN, on-ref) and 561 (PoN)
print('\n==== FP composition (germline = exact HG008-N dipcall SNP allele) ====')


def comp(sites, thr, pon=None, af=None):
    c = collections.Counter()
    for k, s in sites.items():
        if s < thr or k in TRUTH or not inconf(k[0]):
            continue
        if pon is not None and pon.get(k, False):
            continue
        c['germ' if k in GERM else 'nongerm'] += 1
    return dict(c)


print('v1 e053 noPoN (>=0.5)', comp(V1['v1_e053'], 0.5))
print('v1 e053 myPoN (>=0.5)', comp(V1['v1_e053'], 0.5, PONF))
print('v1 e053 ownPoN VCF', comp(V1['v1_e053_ownPoN'], 0.5))
for run, (s, o, nt, neu) in V2R.items():
    r = run.replace('HG008_Illumina_SNV_', 'v2_')
    a = at_tp(rows[r + ' noPoN'], 608); b = at_tp(rows[r + ' myPoN'], 555)
    print(r, 'noPoN@TP608 thr %.4f' % a[0], comp(s, a[0]), 'offref>=thr', int((o >= a[0]).sum()),
          '| myPoN@TP555 thr %.4f' % b[0], comp(s, b[0], PONF), 'offref>=thr', int((o >= b[0]).sum()))
json.dump({k: summ(v) for k, v in rows.items()}, open(OUT + '/my_check.json', 'w'), indent=1)
print('done t=%.0f' % (time.time() - T0))
