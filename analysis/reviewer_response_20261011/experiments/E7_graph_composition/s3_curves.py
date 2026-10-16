"""E7 step 3: in-silico graphs from k random HPRC v1.1 (E7_GRAPH=v2.1: v2.1) haplotypes with a minimum carrier count d.

An allele is 'in the graph' when at least d of the k drawn haplotypes carry it (the GRCh38 backbone is always in). For the
sample's germline alleles this is (a), germline that is no longer a candidate; for its somatic truth alleles it is
(b), somatic alleles the graph absorbs. A truth allele with several primitives (MNP-like records) counts as in the graph
only when one haplotype set carries all of them (AND of the masks).
Inputs: WORK/hprc_alleles/chr*.tsv.gz (s1), WORK/sample_alleles/<sample>.<role>.tsv.gz (s2),
        analysis/graph_absorbed_somatic_20261001/hprc_individuals.tsv (superpopulation per haplotype).
Pools: all 89 non-GRCh38 haplotypes (CHM13 + 44 x 2); AFR (23 x 2 = 46); non-AFR (AMR 16 + EAS 4 + SAS 1 = 42 haplotypes;
CHM13 left out of both population pools). k in KS, 20 random draws per k (seeded), one draw when k is the whole pool.
For every draw the full carrier-count histogram per (sample, role, scope, stratum) is kept, so any threshold can be read:
absolute d in D_ABS and relative f in F_REL (d = max(1, ceil(f k))).
Scopes: bed = FILTER PASS and inside the truth's BED (somatic benchmark BED; dipcall BED for germline), chr1-22;
        bed_chr1 = the same on chr1; all = every chr1-22 allele.
Outputs (OUT = this folder): curves.tsv (figure-ready), graph_size.tsv (HPRC alleles kept per graph);
        WORK/joined.npz (per sample allele: sample, role, truth_id, chrom, stratum, passed, in_bed, carrier mask) for s4.
v2.1 (E7_GRAPH=v2.1, added 2026-10-16): 463 haplotypes (CHM13 + 231 x 2), superpopulations from the HPRC release 2
metadata (common.haplotypes). Pools: all 463; AFR 140, AMR 88, EAS 102, EUR 60, SAS 72 (CHM13 in none of them).
k in 2-256 and 89 (the size of v1.1) for the whole pool, d up to 92; for each superpopulation k in 2-32 plus k = 60 (the EUR pool, so
all five are compared at equal k) plus the whole pool. Graph size only for the all-haplotype pool. Outputs in ./v2.1 and
WORK/v2.1. Carrier masks are NW = ceil(NHAP / 64) uint64 words (joined.npz 'masks'; the v1.1 file of 2026-10-13 has the
same two words as 'lo' / 'hi').
"""
import array, collections, csv, gzip, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import AUTOSOMES, FI, FINE, GOUT, GRAPH, GROUPS, GWORK, NHAP, NW, ROLES, SAMPLES, WORK, fine_code, \
    haplotypes, mask_words

OUT = GOUT
os.makedirs(OUT, exist_ok=True)
JOINED = f'{GWORK}/joined.npz'
CHROMS = AUTOSOMES
if os.environ.get('E7_TEST'):  # quick test: one chromosome, outputs in E7_TEST
    OUT = JOINED = os.environ['E7_TEST']
    JOINED += '/joined.npz'
    CHROMS = ['chr21']
KS = [2, 4, 8, 16, 32, 64] if GRAPH == 'v1.1' else [2, 4, 8, 16, 32, 64, 89, 128, 256]
D_ABS = [1, 2, 5, 9, 18] if GRAPH == 'v1.1' else [1, 2, 5, 9, 18, 46, 92]
F_REL = [0.0, 0.02, 0.05, 0.10, 0.20]  # 0.0 = at least one haplotype
NDRAW = 20
SCOPES = ['all', 'bed', 'bed_chr1']
NB = NHAP + 1  # carrier counts 0..NHAP


# haplotype bit order = VCF column order (CHM13, then hap 1 / hap 2 per sample); same order as hprc_individuals.tsv
haps = haplotypes()
SUPER = [h['superpopulation'] for h in haps]
if GRAPH == 'v1.1':
    POOLS = {'all': list(range(NHAP)),
             'AFR': [i for i, s in enumerate(SUPER) if s == 'AFR'],
             'nonAFR': [i for i, s in enumerate(SUPER) if s in ('AMR', 'EAS', 'SAS')]}
else:
    POOLS = {'all': list(range(NHAP)),
             **{sp: [i for i, s in enumerate(SUPER) if s == sp] for sp in ('AFR', 'AMR', 'EAS', 'EUR', 'SAS')}}
KEQ = min(len(v) for p, v in POOLS.items() if p != 'all')  # v2.1: equal-k comparison of the superpopulations
print({p: len(v) for p, v in POOLS.items()}, 'KEQ', KEQ, flush=True)

# 1. sample alleles (keys per chromosome)
print('loading sample alleles', flush=True)
meta = collections.defaultdict(list)  # per field
keys_by_chrom = collections.defaultdict(list)  # chrom -> [(row index, [keys])]
n = 0
for si, s in enumerate(SAMPLES):
    for ri, role in enumerate(ROLES):
        fn = f'{WORK}/sample_alleles/{s}.{role}.tsv.gz'
        if CHROMS != AUTOSOMES and not os.path.exists(fn):
            continue
        with gzip.open(fn, 'rt') as fh:
            fh.readline()
            for line in fh:
                tid, c, kind, p, b, keys, nk, rc, u, tr, il = line.rstrip('\n').split('\t')
                if c not in CHROMS:
                    continue
                meta['sample'].append(si); meta['role'].append(ri); meta['truth_id'].append(int(tid))
                meta['chrom'].append(AUTOSOMES.index(c)); meta['fine'].append(fine_code(nk, rc, il))
                meta['passed'].append(p == 'T'); meta['in_bed'].append(b == 'T')
                meta['tract'].append(int(tr)); meta['ilen'].append(int(il))
                keys_by_chrom[c].append((n, keys.split(';') if keys else []))
                n += 1
        print(s, role, n, flush=True)
masks = np.zeros((n, NW), np.uint64)
found = np.zeros(n, bool)

# 2. HPRC alleles: join the masks, keep every allele for graph size
print('loading HPRC alleles', flush=True)
h_w, h_fine, h_chr1 = [array.array('Q') for _ in range(NW)], array.array('b'), array.array('b')
for c in CHROMS:
    d = {}
    with gzip.open(f'{GWORK}/hprc_alleles/{c}.tsv.gz', 'rt') as fh:
        fh.readline()
        for line in fh:
            k, kind, rc, u, tr, il, nc, m = line.rstrip('\n').split('\t')
            m = int(m, 16)
            d[k] = m
            for w, x in zip(h_w, mask_words(m)):
                w.append(x)
            h_fine.append(fine_code(kind, rc, il)); h_chr1.append(c == 'chr1')
    for i, ks in keys_by_chrom.pop(c, []):
        ms = [d.get(k) for k in ks]
        if ks and all(m is not None for m in ms):
            m = ms[0]
            for x in ms[1:]:
                m &= x
            masks[i] = mask_words(m)
            found[i] = True
    print(c, len(d), flush=True)
    del d
h_masks = np.stack([np.frombuffer(w, np.uint64) for w in h_w], axis=1) if NW > 1 else \
    np.frombuffer(h_w[0], np.uint64)[:, None]
del h_w
h_fine = np.frombuffer(h_fine, np.int8).astype(np.int64); h_chr1 = np.frombuffer(h_chr1, np.int8).astype(bool)
M = {k: np.array(v) for k, v in meta.items()}
np.savez_compressed(JOINED, masks=masks, found=found, **M)
print('joined', n, 'found in HPRC', int(found.sum()), flush=True)

# 3. draws
rng = np.random.default_rng(20261013)


def count(mw, sel):
    """carriers among the drawn haplotypes, per row of the (n, NW) mask array"""
    s = 0
    for h in sel:
        s |= 1 << int(h)
    out = np.zeros(len(mw), np.int64)
    for w, x in enumerate(mask_words(s)):
        if x:
            out += np.bitwise_count(mw[:, w] & np.uint64(x))
    return out


ng = len(SAMPLES) * len(ROLES) * len(FINE)
cell = (M['sample'] * len(ROLES) + M['role']) * len(FINE) + M['fine']
inbed = M['passed'] & M['in_bed']
chr1 = M['chrom'] == 0
scope_idx = {'all': np.arange(n), 'bed': np.flatnonzero(inbed), 'bed_chr1': np.flatnonzero(inbed & chr1)}
draws = []  # (pool, k, draw, haplotype list)
for pool, members in POOLS.items():
    ks = [k for k in KS if k < len(members)] + [len(members)]
    if GRAPH == 'v1.1' and pool == 'AFR':
        ks = sorted(set(ks) | {len(POOLS['nonAFR'])})
    if GRAPH != 'v1.1' and pool != 'all':
        ks = sorted({k for k in KS if k <= 32} | {KEQ, len(members)})
    for k in ks:
        for r in range(1 if k == len(members) else NDRAW):
            draws.append((pool, k, r, sorted(rng.choice(members, size=k, replace=False).tolist())))
print('draws', len(draws), flush=True)
hist = {}  # (pool, k, draw, scope) -> (ng, NB) counts
ghist = {}  # (pool, k, draw, region) -> (len(FINE), NB) for HPRC alleles
for pool, k, r, sel in draws:
    cnt = count(masks, sel)
    for sc, idx in scope_idx.items():
        hist[(pool, k, r, sc)] = np.bincount(cell[idx] * NB + cnt[idx], minlength=ng * NB).reshape(ng, NB)
    if GRAPH != 'v1.1' and pool != 'all':
        continue
    gc = count(h_masks, sel)
    for reg, idx in (('chr1-22', slice(None)), ('chr1', h_chr1)):
        ghist[(pool, k, r, reg)] = np.bincount(h_fine[idx] * NB + gc[idx], minlength=len(FINE) * NB).reshape(len(FINE), NB)
np.save(f'{OUT}/hist_draws.npy' if CHROMS != AUTOSOMES else f'{GWORK}/hist_draws.npy', {'hist': hist, 'ghist': ghist, 'draws': draws}, allow_pickle=True)


def thresholds(k):
    out = [('abs', d, d) for d in D_ABS if d <= k]
    out += [('rel', f, max(1, math.ceil(f * k - 1e-9))) for f in F_REL]
    return out


# 4. curves.tsv
rows = []
by = collections.defaultdict(list)
for (pool, k, r, sc), h in hist.items():
    by[(pool, k, sc)].append(h)
for (pool, k, sc), hs in sorted(by.items()):
    H = np.stack(hs)  # draws x ng x NB
    suffix = H[:, :, ::-1].cumsum(axis=2)[:, :, ::-1]  # n with count >= j
    for si, s in enumerate(SAMPLES):
        for ri, role in enumerate(ROLES):
            base = (si * len(ROLES) + ri) * len(FINE)
            for g, members in GROUPS.items():
                fi = [base + FI[f] for f in members]
                tot = H[0, fi, :].sum()
                if tot == 0:
                    continue
                for tt, tv, d in thresholds(k):
                    frac = suffix[:, fi, d].sum(axis=1) / tot
                    rows.append(dict(sample=s, role=role, scope=sc, stratum=g, pool=pool, k=k, threshold_type=tt,
                                     threshold=tv, d=d, n_alleles=int(tot), n_draws=len(hs),
                                     mean=round(float(frac.mean()), 6), sd=round(float(frac.std(ddof=1)) if len(hs) > 1 else 0.0, 6),
                                     min=round(float(frac.min()), 6), max=round(float(frac.max()), 6)))
with open(f'{OUT}/curves.tsv', 'w') as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print('curves rows', len(rows), flush=True)

# 5. graph_size.tsv: HPRC alleles kept by each in-silico graph
rows = []
gby = collections.defaultdict(list)
for (pool, k, r, reg), h in ghist.items():
    gby[(pool, k, reg)].append(h)
GG = {'SNV': ['SNV'], 'INDEL<=50': GROUPS['INDEL'], 'INDEL_HP>=7': GROUPS['INDEL_HP>=7'], 'SV>50': ['SV>50'],
      'other': ['other']}
for (pool, k, reg), hs in sorted(gby.items()):
    H = np.stack(hs)
    suffix = H[:, :, ::-1].cumsum(axis=2)[:, :, ::-1]
    for g, members in GG.items():
        fi = [FI[f] for f in members]
        for tt, tv, d in thresholds(k):
            v = suffix[:, fi, d].sum(axis=1)
            rows.append(dict(region=reg, stratum=g, pool=pool, k=k, threshold_type=tt, threshold=tv, d=d,
                             n_draws=len(hs), mean_alleles=round(float(v.mean()), 1),
                             sd=round(float(v.std(ddof=1)) if len(hs) > 1 else 0.0, 1),
                             hprc_alleles_total=int(H[0, fi, :].sum())))
with open(f'{OUT}/graph_size.tsv', 'w') as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print('done', flush=True)
