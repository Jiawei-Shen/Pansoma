"""Read-only verification of v1 (.dat/.idx) HG008 Illumina SNV artifacts against the v2 697-truth definition.
Outputs to stdout. No writes outside this scratch dir."""
import json, gzip, bisect, collections, sys, os, glob, time
import pysam

T0 = time.time()
V2 = '/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors'
DV = '/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/'
TRUTH_ALL = DV + 'HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz'
TRUTH_SNV = DV + 'HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.snv.vcf.gz'   # the file v1 labels used
SBED = DV + 'HG008-T_somatic_smvar_benchmark_v0.2_all.bed'
GBED = '/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed'
GVCF = '/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz'
P1 = '/scratch/jshen/data/Pansoma/HG008_GIAB/'
R2 = '/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/'
CHROMS = [f'chr{i}' for i in range(1, 23)]


def bed(path):
    iv = collections.defaultdict(list)
    for l in open(path):
        if l.startswith(('#', 'track', 'browser')):
            continue
        f = l.split('\t')
        iv[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, v in iv.items():
        v.sort(); m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        out[c] = m
    return out


def intersect(a, b):
    i = j = 0; out = []
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


sb, gb = bed(SBED), bed(GBED)
CONF = {c: intersect(sb.get(c, []), gb.get(c, [])) for c in CHROMS}
STARTS = {c: [s for s, e in CONF[c]] for c in CHROMS}


def inbed(c, pos1):
    p = pos1 - 1; st = STARTS.get(c)
    if not st:
        return False
    k = bisect.bisect_right(st, p) - 1
    return k >= 0 and CONF[c][k][0] <= p < CONF[c][k][1]


print('confident bp chr1', sum(e - s for s, e in CONF['chr1']))

# 697 truths (v2 definition)
rows = [l.rstrip('\n').split('\t') for l in open(V2 + '/somatic.recall.tsv')]
hdr = rows[0]; rows = [dict(zip(hdr, r)) for r in rows[1:]]
t697 = [r for r in rows if r['chrom'] == 'chr1' and r['kind'] == 'SNP' and r['passed'] == 'True' and r['in_bed'] == 'True']
print('n697', len(t697), collections.Counter(r['status'] for r in t697))
pos697 = {int(r['vcf_pos']): r['truth_id'] for r in t697}
truth = {}   # (pos, ref, alt) -> truth_id  (chr1, 697 only)
for r in pysam.VariantFile(TRUTH_ALL).fetch('chr1'):
    if r.pos in pos697 and len(r.ref) == 1 and all(len(a) == 1 for a in r.alts):
        for a in r.alts:
            truth[(r.pos, r.ref.upper(), a.upper())] = pos697[r.pos]
print('697 truth allele keys', len(truth), 'distinct', len(set(truth.values())))

# all somatic truth records, per chrom (to classify v1 'true' labels)
snv_keys = collections.defaultdict(set); snv_filters = collections.Counter(); allsom_pos = collections.defaultdict(set)
for r in pysam.VariantFile(TRUTH_SNV):
    snv_filters[tuple(r.filter.keys()) or ('.',)] += 1
    for a in r.alts:
        snv_keys[r.chrom].add((r.pos, r.ref.upper(), a.upper()))
for r in pysam.VariantFile(TRUTH_ALL):
    allsom_pos[r.chrom].add(r.pos)
print('snv truth VCF records', sum(snv_filters.values()), 'FILTER values', dict(snv_filters), 'chr1 SNV alleles', len(snv_keys['chr1']))

# germline SNP alleles (HG008-N dipcall), all autosomes
germ = collections.defaultdict(set)
for r in pysam.VariantFile(GVCF):
    if r.chrom not in CONF:
        continue
    if len(r.ref) == 1:
        for a in r.alts or []:
            if len(a) == 1 and a != '*':
                germ[r.chrom].add((r.pos, r.ref.upper(), a.upper()))
print('germline SNP alleles chr1', len(germ['chr1']), 'autosomes', sum(len(v) for v in germ.values()), f't={time.time()-T0:.0f}s')


def classify_set(path, chrom):
    """v1 variant_summary_classified.ndjson -> Counter[(label, bed, class)], truth ids found."""
    c = collections.Counter(); found = set(); n = 0
    for l in open(path):
        r = json.loads(l); n += 1
        gp = r.get('genomic_position')
        if gp is None:
            c[(r['classification'], 'noPos', '-')] += 1; continue
        k = (gp, r['v_ref'].upper(), r['v_alt'].upper())
        ib = 'inBED' if inbed(chrom, gp) else 'outBED'
        if chrom == 'chr1' and k in truth:
            cls = 'truth697'; found.add(truth[k])
        elif k in snv_keys[chrom]:
            cls = 'other_somatic_SNV_truth'
        elif k in germ[chrom]:
            cls = 'HG008N_germline'
        elif gp in allsom_pos[chrom]:
            cls = 'pos_of_other_somatic'
        else:
            cls = 'other'
        c[(r['classification'], ib, cls)] += 1
        # sanity: v1 'true' should equal exact match to the snv truth VCF
        if (r['classification'] == 'true') != (k in snv_keys[chrom]):
            c[('LABEL_RULE_MISMATCH', ib, cls)] += 1
    return c, found, n


for design in ['AF_HPRC', 'AF_HPRC_HG008N_added']:
    vpath = P1 + design + '/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson'
    c, found, n = classify_set(vpath, 'chr1')
    print(f'\n[{design}] chr1 val tensors {n}; by (label, BED, class):')
    for k, v in sorted(c.items()):
        print('   ', k, v)
    print(f'[{design}] ceiling: {len(found)} of 697 = {len(found)/697:.4f}')
    tot = collections.Counter()
    for ch in CHROMS[1:]:
        tp = P1 + design + f'/5ch_training_data_SNV/train/SNV_{ch}_5chan_tensor_dataset/variant_summary_classified.ndjson'
        if not os.path.exists(tp):
            print('   missing', tp); continue
        cc, _, nn = classify_set(tp, ch)
        for (lab, ib, cls), v in cc.items():
            tot[(lab, ib, cls)] += v
    print(f'[{design}] train chr2-22 by (label, BED, class):')
    for k, v in sorted(tot.items()):
        print('   ', k, v)
    t_true = sum(v for (l, b, c_), v in tot.items() if l == 'true')
    t_germ = sum(v for (l, b, c_), v in tot.items() if l == 'false' and c_ == 'HG008N_germline')
    print(f'[{design}] train true {t_true}, false-that-are-HG008N-germline {t_germ} (ratio {t_germ/max(1,t_true):.1f}:1)', f't={time.time()-T0:.0f}s')

# v1 test set chr1 (AF_HPRC) vs val: map GRCh38 node starts from val
val_start = {}
for l in open(P1 + 'AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson'):
    r = json.loads(l)
    if r.get('genomic_position') is not None:
        val_start[r['node_id']] = r['genomic_position'] - r['v_pos']
tc = collections.Counter(); tfound = set()
for l in open(P1 + 'AF_HPRC/5ch_testing_data_SNV/5ch_testing_data_SNV_chr1/variant_summary.ndjson'):
    r = json.loads(l)
    s = val_start.get(r['node_id'])
    if s is None:
        tc['node_not_in_val_(off-GRCh38 or new)'] += 1; continue
    k = (s + r['v_pos'], r['v_ref'].upper(), r['v_alt'].upper())
    tc['on_val_node'] += 1
    if k in truth:
        tfound.add(truth[k])
print('\n[AF_HPRC] chr1 testing set:', dict(tc), 'truth697 found', len(tfound))


# ---- existing v1 prediction VCFs on chr1 ----
VCFS = {
    'AF_HPRC e053 linear (May14, ALL_chr)': R2 + 'AF_HPRC/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.vcf.gz',
    'AF_HPRC e053 linear+4PoN (Apr9, ALL_chr)': R2 + 'AF_HPRC/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.filtered_PoN.vcf.gz',
    'AF_HPRC e053 linear (Mar20, Chr1)': R2 + 'AF_HPRC/pansoma_HG008T_WGS_Chr1/pansoma-to_SNV_pansoma_HG008T_WGS_Chr1.linear.vcf.gz',
    'AF_HPRC e053 linear+4PoN (Mar19, Chr1)': R2 + 'AF_HPRC/pansoma_HG008T_WGS_Chr1/pansoma-to_SNV_pansoma_HG008T_WGS_Chr1.linear.filtered_PoN.vcf.gz',
    'AF_HPRC run1 linear (Mar9, tmp)': '/scratch/jshen/tmp/HG008T_WGS_AF-HPRC_pansoma.test.linear.vcf.gz',
    'AF_HPRC run1 linear+4PoN (Mar9, tmp)': '/scratch/jshen/tmp/HG008T_WGS_AF-HPRC_pansoma.test.linear.filtered_PoN.vcf.gz',
    'HG008N_added linear (Apr5)': R2 + 'AF_HPRC_HG008N_added/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.vcf.gz',
    'HG008N_added linear+4PoN (Apr6)': R2 + 'AF_HPRC_HG008N_added/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.filtered_PoN.vcf.gz',
    'HG008N_added e064 linear (Apr7, gthreshold02)': R2 + 'AF_HPRC_HG008N_added/pansoma_HG008T_WGS_ALL_chr_gthreshold02/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.vcf.gz',
}


def curve(items):
    items.sort(key=lambda x: -x[0])
    tps = set(); fp = 0; pts = []; last = {}
    for p, x, t in items:
        if t == 'TP':
            tps.add(x)
        else:
            fp += 1
        last[p] = (len(tps), fp)
    for th in sorted(last, reverse=True):
        tp, f = last[th]; P = tp / (tp + f); R = tp / 697; F = 2 * P * R / (P + R) if P + R else 0
        pts.append((th, tp, f, P, R, F))
    out = {'maxF1': max(pts, key=lambda x: x[5])}
    for pr in [0.05, 0.10, 0.15, 0.20]:
        c = [x for x in pts if x[3] >= pr]; out[f'R@P{pr}'] = round(max(c, key=lambda x: x[4])[4], 4) if c else None
    for rc in [0.5, 0.7, 0.8, 0.9]:
        c = [x for x in pts if x[4] >= rc]; out[f'P@R{rc}'] = round(max(c, key=lambda x: x[3])[3], 4) if c else None
    out['lowest_threshold'] = pts[-1]
    return out


print('\nv1 VCFs, strict counting (TP = distinct of the 697; FP = distinct called chr1 SNV alleles in somatic BED ∩ germline BED not among the 697):')
for name, path in VCFS.items():
    if not os.path.exists(path):
        print(name, 'MISSING', path); continue
    calls = {}
    for r in pysam.VariantFile(path).fetch('chr1'):
        for a in r.alts:
            if len(r.ref) != 1 or len(a) != 1:
                continue
            k = (r.pos, r.ref.upper(), a.upper())
            calls[k] = max(calls.get(k, 0.0), float(r.info.get('PROB', 1.0)))
    items = []; fpc = collections.Counter(); n_all_fp_nobed = 0
    for k, p in calls.items():
        if k in truth:
            items.append((p, truth[k], 'TP'))
        else:
            if k not in snv_keys['chr1']:
                n_all_fp_nobed += 1
            if inbed('chr1', k[0]):
                cls = 'germ' if k in germ['chr1'] else ('somother' if k[0] in allsom_pos['chr1'] else 'other')
                items.append((p, k, cls)); fpc[cls] += 1
    o = curve(items)
    th, tp, f, P, R, F = o['maxF1']; lt = o['lowest_threshold']
    print(f'* {name}: calls {len(calls)} minPROB {min(calls.values()):.3f}; at lowest threshold TP {lt[1]} FP {lt[2]} (HG008N germline {fpc["germ"]}, other-somatic-pos {fpc["somother"]}, other {fpc["other"]}) P {lt[3]:.4f} R {lt[4]:.4f} F1 {lt[5]:.4f}; '
          f'maxF1 {F:.4f} (P {P:.4f} R {R:.4f} thr {th:.3f}); ' + ' '.join(f'{k}={v}' for k, v in o.items() if k.startswith(('R@', 'P@')))
          + f'; [no-BED FP vs snv VCF: {n_all_fp_nobed}]')
print(f'done t={time.time()-T0:.0f}s')
