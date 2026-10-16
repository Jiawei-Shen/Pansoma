"""E7 step 6 (R3.m3): is any HPRC v1.1 individual a close relative of a benchmark donor?

A first-degree relative shares about half of a person's rare alleles; unrelated people share only what their ancestry
makes common. For each benchmark normal (HG008-N dipcall, COLO829BL dipcall, HCC1395BL DeepVariant; FILTER PASS, inside
the BED, chr1-22) take its germline SNVs that are rare in HPRC v1.1 (carried by 1-2 of the 89 non-GRCh38 haplotypes),
and count for every HPRC individual how many of them it carries (either haplotype). Output: relatedness.tsv (this folder):
per sample and HPRC individual, the shared count and its ratio to the median over individuals of the same superpopulation.
Input: WORK/joined.npz (s3).
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import AUTOSOMES, FI, SAMPLES, WORK

OUT = os.path.dirname(os.path.abspath(__file__))
JOINED = f'{WORK}/joined.npz'
CHROMS = AUTOSOMES
if os.environ.get('E7_TEST'):  # quick test with s3's test output (chr21 only)
    OUT = os.environ['E7_TEST']
    JOINED = f'{OUT}/joined.npz'
    CHROMS = ['chr21']
haps = list(csv.DictReader(open('/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/'
                                'hprc_individuals.tsv'), delimiter='\t'))
J = np.load(JOINED)
cnt = (np.bitwise_count(J['lo']) + np.bitwise_count(J['hi'])).astype(np.int64)
rows = []
for si, s in enumerate(SAMPLES):
    sel = ((J['sample'] == si) & (J['role'] == 0) & J['passed'] & J['in_bed'] & (J['fine'] == FI['SNV']) &
           J['found'] & (cnt >= 1) & (cnt <= 2))
    lo, hi = J['lo'][sel], J['hi'][sel]
    per_hap = []
    for b in range(89):
        w, bit = (lo, b) if b < 64 else (hi, b - 64)
        per_hap.append(int(((w >> np.uint64(bit)) & np.uint64(1)).sum()))
    ind = {}
    for h, v in zip(haps, per_hap):
        ind.setdefault((h['sample'], h['superpopulation']), 0)
        ind[(h['sample'], h['superpopulation'])] += v
    med = {}
    for (name, sp), v in ind.items():
        med.setdefault(sp, []).append(v)
    med = {sp: float(np.median(v)) for sp, v in med.items()}
    for (name, sp), v in sorted(ind.items(), key=lambda x: -x[1]):
        rows.append(dict(sample=s, hprc_individual=name, superpopulation=sp, rare_snvs=int(sel.sum()), shared=v,
                         share=round(v / max(1, int(sel.sum())), 4), ratio_to_superpop_median=round(v / med[sp], 2) if med[sp] else ''))
with open(f'{OUT}/relatedness.tsv', 'w') as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
for s in SAMPLES:
    top = [r for r in rows if r['sample'] == s][:3]
    print(s, [(r['hprc_individual'], r['superpopulation'], r['shared'], r['share'], r['ratio_to_superpop_median']) for r in top])
