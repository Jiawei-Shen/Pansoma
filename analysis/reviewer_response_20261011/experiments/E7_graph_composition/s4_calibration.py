"""E7 step 4: place the real graphs on the in-silico curves, allele by allele.

For every truth allele (scope: FILTER PASS, inside its BED, chr1-22; and chr1 only) it joins
  - the in-silico prediction for the full HPRC v1.1 haplotype set: carrier count c among the 89 non-GRCh38 haplotypes
    (s3 joined.npz); 'predicted in d9' = c >= 9;
  - whether the real d9 graph's own deconstruct VCF has the allele as a simple key (s5; lower bound, complex ALT
    alleles of d9 are not decomposed);
  - the measured status in the tensor sets (tensors/{somatic,germline}.recall.tsv, the source of
    analysis/graph_ablation_tensors_20261005/tables.md section B): 'tensor' = status tensor_representative.
Measured graph loss = alleles with a tensor on the GRCh38-only (linear) arm but none on d9 (or d46), over alleles with a
linear tensor. Outputs: calibration.tsv, calibration_bins.tsv, d9_retention.tsv, d9_content.tsv (this folder).
v2.1 (E7_GRAPH=v2.1, added 2026-10-16): carrier count c among the 463 v2.1 haplotypes; 'predicted in d46' = c >= 46
(the d46 graph keeps nodes of >= 46 of its 464 haplotypes; GRCh38 never carries an ALT). HPRC publishes no deconstruct
VCF of the d46 graph, so there is no graph-content check (steps 2-3); the calibration is the measured candidate loss.
Outputs in ./v2.1: calibration.tsv, calibration_bins.tsv (v2.1 carrier bins), d46_content.tsv.
"""
import collections, csv, gzip, hashlib, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import AUTOSOMES, FI, FINE, G, GOUT, GRAPH, GROUPS, GWORK, ROLES, SAMPLES, WORK, fine_code

OUT = GOUT
JOINED = f'{GWORK}/joined.npz'
R, DR = G['real'], G['d_real']  # the real graph of this HPRC release and its carrier threshold (d9: 9, d46: 46)
CHROMS = AUTOSOMES
if os.environ.get('E7_TEST'):  # quick test with s3's test output (chr21 only)
    OUT = os.environ['E7_TEST']
    JOINED = f'{OUT}/joined.npz'
    CHROMS = ['chr21']
SETS = '/scratch/jshen/data/pansoma_v2_tensors'
PLATFORMS = {'HG008': [('PacBio', 'HG008T_PacBio'), ('ONT', 'HG008T_ONT'), ('Illumina', 'HG008T_Illumina')],
             'COLO829': [('PacBio_fiberseq', 'COLO829T_fiberseq'), ('ONT', 'COLO829T_ONT'), ('Illumina', 'COLO829T_Illumina')],
             'HCC1395': [('PacBio', 'HCC1395T_PacBio'), ('ONT', 'HCC1395T_ONT'), ('Illumina', 'HCC1395T_Illumina')]}
ARMS = {'linear': '_linear', 'd9': '', 'd46': '_d46'}
if os.environ.get('E7_TEST'):
    PLATFORMS = {k: v[:1] for k, v in PLATFORMS.items()}
CBINS = [(0, 0), (1, 1), (2, 4), (5, 8), (9, 17), (18, 44), (45, 89)] if GRAPH == 'v1.1' else \
    [(0, 0), (1, 1), (2, 9), (10, 22), (23, 34), (35, 45), (46, 57), (58, 92), (93, 231), (232, 463)]

J = np.load(JOINED)
W = J['masks'] if 'masks' in J.files else np.stack([J['lo'], J['hi']], axis=1)
n = len(W)
cnt = sum(np.bitwise_count(W[:, w]).astype(np.int64) for w in range(W.shape[1]))
cnt[~J['found']] = 0
sample, role, tid, chrom, fine = J['sample'], J['role'], J['truth_id'], J['chrom'], J['fine']
inbed = J['passed'] & J['in_bed']
index = {}
for si in range(len(SAMPLES)):
    for ri in range(len(ROLES)):
        sel = np.flatnonzero((sample == si) & (role == ri))
        index[(si, ri)] = (sel, dict(zip(tid[sel].tolist(), range(len(sel)))))

# 1. measured statuses
ids = collections.defaultdict(set)
status = {}  # (sample, platform, arm) -> int8 array over all n (-1 = no table): 1 tensor, 0 no tensor
for si, s in enumerate(SAMPLES):
    for plat, base in PLATFORMS[s]:
        for arm, suf in ARMS.items():
            d = f'{SETS}/{base}{suf}/tensors'
            if not os.path.exists(f'{d}/somatic.recall.tsv'):
                continue
            st = np.full(n, -1, np.int8)
            for ri, r in enumerate(ROLES):
                sel, pos = index[(si, ri)]
                with open(f'{d}/{r}.recall.tsv') as fh:
                    hdr = fh.readline().rstrip('\n').split('\t')
                    it, ist = hdr.index('truth_id'), hdr.index('status')
                    miss = 0
                    hs = hashlib.md5()
                    for line in fh:
                        f = line.split('\t', ist + 1)
                        hs.update('\t'.join(f[:4]).encode())
                        j = pos.get(int(f[it]))
                        if j is None:
                            miss += 1
                            continue
                        st[sel[j]] = 1 if f[ist].rstrip('\n') == 'tensor_representative' else 0
                ids[(s, r)].add(hs.hexdigest())
            status[(s, plat, arm)] = st
            print(s, plat, arm, 'rows not in joined (non-autosomal / complex):', miss, flush=True)
for k, v in ids.items():
    print('truth table identity (truth_id, chrom, vcf_pos, kind) across sets', k, 'distinct hashes:', len(v), flush=True)
    if len(v) != 1:
        print('WARNING: truth tables differ between sets for', k, flush=True)

# 2. d9 deconstruct presence (simple keys); v1.1 only
print('d9 keys', flush=True)
in_d9vcf = np.zeros(n, bool)
keys_by_chrom = collections.defaultdict(list)
for si, s in (enumerate(SAMPLES) if GRAPH == 'v1.1' else []):
    for ri, r in enumerate(ROLES):
        sel, pos = index[(si, ri)]
        fn = f'{WORK}/sample_alleles/{s}.{r}.tsv.gz'
        if CHROMS != AUTOSOMES and not os.path.exists(fn):
            continue
        with gzip.open(fn, 'rt') as fh:
            fh.readline()
            for line in fh:
                f = line.split('\t', 6)
                j = pos.get(int(f[0]))
                if f[1] not in CHROMS or j is None:
                    continue
                keys_by_chrom[f[1]].append((sel[j], f[5].split(';') if f[5] else []))
ret = collections.Counter()  # (fine, c) -> [n, in d9]
for c in (CHROMS if GRAPH == 'v1.1' else []):
    d9 = set()
    with gzip.open(f'{WORK}/d9_alleles/{c}.tsv.gz', 'rt') as fh:
        fh.readline()
        for line in fh:
            d9.add(line.split('\t', 1)[0])
    for i, ks in keys_by_chrom.pop(c, []):
        in_d9vcf[i] = bool(ks) and all(k in d9 for k in ks)
    with gzip.open(f'{WORK}/hprc_alleles/{c}.tsv.gz', 'rt') as fh:
        fh.readline()
        for line in fh:
            k, kind, rc, u, tr, il, nc, m = line.rstrip('\n').split('\t')
            fc = fine_code(kind, rc, il)
            ret[(fc, int(nc), 'n')] += 1
            ret[(fc, int(nc), 'd9')] += k in d9
    print(c, len(d9), flush=True)

# 3. d9 retention per HPRC carrier count
rows = []
for g in ['SNV', 'INDEL', 'INDEL_HP>=7', 'INDEL_not_HP>=7', 'SV>50']:
    fis = [FI[f] for f in GROUPS[g]]
    for c in range(90):
        nn = sum(ret[(f, c, 'n')] for f in fis)
        if nn:
            rows.append(dict(stratum=g, hprc_carriers=c, hprc_alleles=nn, in_d9_vcf=sum(ret[(f, c, 'd9')] for f in fis),
                             fraction=round(sum(ret[(f, c, 'd9')] for f in fis) / nn, 4)))
if rows:
    with open(f'{OUT}/d9_retention.tsv', 'w') as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)

# 4. calibration per (sample, platform, role, stratum, scope)
rows, brows, crows = [], [], []
for si, s in enumerate(SAMPLES):
    for ri, r in enumerate(ROLES):
        for scope in ('bed', 'bed_chr1'):
            base = (sample == si) & (role == ri) & inbed & ((chrom == 0) if scope == 'bed_chr1' else True)
            for g, members in GROUPS.items():
                gm = base & np.isin(fine, [FI[f] for f in members])
                N = int(gm.sum())
                if not N:
                    continue
                pred9 = cnt[gm] >= DR
                crow = dict(sample=s, role=r, scope=scope, stratum=g, n=N,
                            in_hprc_any=round(float((cnt[gm] >= 1).mean()), 4),
                            **{f'predicted_{R}': round(float(pred9.mean()), 4)})
                if GRAPH == 'v1.1':
                    crow.update(in_d9_vcf=round(float(in_d9vcf[gm].mean()), 4),
                                in_d9_vcf_given_pred=round(float(in_d9vcf[gm][pred9].mean()), 4) if pred9.any() else '',
                                in_d9_vcf_given_not_pred=round(float(in_d9vcf[gm][~pred9].mean()), 4) if (~pred9).any() else '')
                crows.append(crow)
                for plat, _ in PLATFORMS[s]:
                    st = {a: status.get((s, plat, a)) for a in ARMS}
                    rec = dict(sample=s, platform=plat, role=r, scope=scope, stratum=g, n=N,
                               **{f'predicted_in_{R}': round(float(pred9.mean()), 4)})
                    for a, v in st.items():
                        rec[f'tensor_{a}'] = round(float((v[gm] == 1).mean()), 4) if v is not None else ''
                    lin = st['linear']
                    for a in ('d9', 'd46'):
                        v = st[a]
                        if lin is None or v is None:
                            rec[f'loss_{a}_vs_linear'] = ''
                            continue
                        has = gm & (lin == 1)
                        lost = has & (v != 1)
                        rec[f'loss_{a}_vs_linear'] = round(float(lost.sum() / max(1, has.sum())), 4)
                        if a == R:
                            p = has & (cnt >= DR)
                            q = has & (cnt < DR)
                            rec[f'loss_{R}_if_pred_in'] = round(float((lost & p).sum() / max(1, p.sum())), 4)
                            rec[f'loss_{R}_if_pred_out'] = round(float((lost & q).sum() / max(1, q.sum())), 4)
                            rec['n_linear_tensor'] = int(has.sum())
                            rec[f'share_of_{R}_loss_pred_in'] = round(float((lost & p).sum() / max(1, lost.sum())), 4)
                    rows.append(rec)
                    if scope == 'bed' and lin is not None:
                        for lo_, hi_ in CBINS:
                            b = gm & (lin == 1) & (cnt >= lo_) & (cnt <= hi_)
                            nb = int(b.sum())
                            br = dict(sample=s, platform=plat, role=r, stratum=g, hprc_carriers=f'{lo_}-{hi_}' if lo_ != hi_ else str(lo_),
                                      n_linear_tensor=nb)
                            for a in ('d9', 'd46'):
                                v = st[a]
                                br[f'loss_{a}'] = round(float((b & (v != 1)).sum() / nb), 4) if (v is not None and nb) else ''
                            brows.append(br)
for name, rr in (('calibration.tsv', rows), ('calibration_bins.tsv', brows), (f'{R}_content.tsv', crows)):
    keys = []
    for x in rr:
        for k in x:
            if k not in keys:
                keys.append(k)
    with open(f'{OUT}/{name}', 'w') as fo:
        w = csv.DictWriter(fo, fieldnames=keys, delimiter='\t', restval=''); w.writeheader(); w.writerows(rr)
print('done', flush=True)
