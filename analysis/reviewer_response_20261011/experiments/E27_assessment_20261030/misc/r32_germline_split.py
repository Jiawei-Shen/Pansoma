"""R3.2 fact check (2026-10-30, read-only): candidate retention of truth alleles split by whether the allele is
represented in the graph, per sample x platform x reference arm.

Same inputs and definitions as E7 s4_calibration.py (analysis/reviewer_response_20261011/experiments/E7_graph_composition):
  - alleles: E7 joined.npz (v1.1 and v2.1; same row order, checked), scope = FILTER PASS & in BED, chr1-22 ('bed')
    and chr1 only ('bed_chr1');
  - 'tensor' = status tensor_representative in <set>/tensors/{germline,somatic}.recall.tsv;
  - represented, three definitions:
      v11_pred: >= 9 of the 89 HPRC v1.1 haplotypes carry the allele (E7 'predicted in d9');
      v11_vcf:  allele (all primitive keys) is in the real d9 graph's deconstruct VCF (E7 s5; lower bound);
      v21_pred: >= 46 of the 463 HPRC v2.1 haplotypes carry it (E7 'predicted in d46'; no d46 deconstruct VCF exists).
Output: TSV to argv[1].
"""
import collections, gzip, sys
import numpy as np

sys.path.insert(0, '/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E7_graph_composition')
from common import AUTOSOMES, FI, GROUPS, SAMPLES, ROLES, WORK  # noqa: E402

SETS = '/scratch/jshen/data/pansoma_v2_tensors'
PLATFORMS = {'HG008': [('PacBio', 'HG008T_PacBio'), ('ONT', 'HG008T_ONT'), ('Illumina', 'HG008T_Illumina')],
             'COLO829': [('fiberseq', 'COLO829T_fiberseq'), ('ONT', 'COLO829T_ONT'), ('Illumina', 'COLO829T_Illumina')]}
ARMS = {'linear': '_linear', 'd9': '', 'd46': '_d46'}


def counts(W):
    return sum(np.bitwise_count(W[:, w]).astype(np.int64) for w in range(W.shape[1]))


J1 = np.load(f'{WORK}/joined.npz')
J2 = np.load(f'{WORK}/v2.1/joined.npz')
for k in ('sample', 'role', 'truth_id', 'chrom', 'fine', 'passed', 'in_bed'):
    assert np.array_equal(J1[k], J2[k]), k
c89 = counts(J1['masks']); c89[~J1['found']] = 0
c463 = counts(J2['masks']); c463[~J2['found']] = 0
sample, role, tid, chrom, fine = J1['sample'], J1['role'], J1['truth_id'], J1['chrom'], J1['fine']
inbed = J1['passed'] & J1['in_bed']
n = len(sample)
index = {}
for si in range(len(SAMPLES)):
    for ri in range(len(ROLES)):
        sel = np.flatnonzero((sample == si) & (role == ri))
        index[(si, ri)] = (sel, dict(zip(tid[sel].tolist(), range(len(sel)))))
print('joined loaded', n, flush=True)

# real d9 deconstruct VCF presence (as E7 s4 step 2)
in_d9vcf = np.zeros(n, bool)
keys_by_chrom = collections.defaultdict(list)
for si, s in enumerate(SAMPLES):
    if s not in PLATFORMS:
        continue
    for ri, r in enumerate(ROLES):
        sel, pos = index[(si, ri)]
        with gzip.open(f'{WORK}/sample_alleles/{s}.{r}.tsv.gz', 'rt') as fh:
            fh.readline()
            for line in fh:
                f = line.split('\t', 6)
                j = pos.get(int(f[0]))
                if f[1] not in AUTOSOMES or j is None:
                    continue
                keys_by_chrom[f[1]].append((sel[j], f[5].split(';') if f[5] else []))
for c in AUTOSOMES:
    d9 = set()
    with gzip.open(f'{WORK}/d9_alleles/{c}.tsv.gz', 'rt') as fh:
        fh.readline()
        for line in fh:
            d9.add(line.split('\t', 1)[0])
    for i, ks in keys_by_chrom.pop(c, []):
        in_d9vcf[i] = bool(ks) and all(k in d9 for k in ks)
    print(c, 'd9 keys', len(d9), flush=True)

REP = {'v11_pred': c89 >= 9, 'v11_vcf': in_d9vcf, 'v21_pred': c463 >= 46}

out = open(sys.argv[1], 'w')
out.write('sample\tplatform\trole\tscope\tstratum\tdefinition\trepresented\tn\ttensor_linear\ttensor_d9\ttensor_d46\n')
for si, s in enumerate(SAMPLES):
    if s not in PLATFORMS:
        continue
    for plat, base in PLATFORMS[s]:
        st = {}
        for arm, suf in ARMS.items():
            v = np.full(n, -1, np.int8)
            for ri, r in enumerate(ROLES):
                sel, pos = index[(si, ri)]
                with open(f'{SETS}/{base}{suf}/tensors/{r}.recall.tsv') as fh:
                    hdr = fh.readline().rstrip('\n').split('\t')
                    it, ist = hdr.index('truth_id'), hdr.index('status')
                    for line in fh:
                        f = line.split('\t', ist + 1)
                        j = pos.get(int(f[it]))
                        if j is not None:
                            v[sel[j]] = 1 if f[ist].rstrip('\n') == 'tensor_representative' else 0
            st[arm] = v
            print(s, plat, arm, 'loaded', flush=True)
        for ri, r in enumerate(ROLES):
            for scope in ('bed', 'bed_chr1'):
                base_m = (sample == si) & (role == ri) & inbed & ((chrom == 0) if scope == 'bed_chr1' else True)
                for g in ('SNV', 'INDEL'):
                    gm = base_m & np.isin(fine, [FI[f] for f in GROUPS[g]])
                    for dname, rep in REP.items():
                        for flag, m in (('yes', gm & rep), ('no', gm & ~rep), ('all', gm)):
                            if dname != 'v11_pred' and flag == 'all':
                                continue
                            row = [s, plat, r, scope, g, dname if flag != 'all' else 'all', flag, int(m.sum())]
                            row += [int((m & (st[a] == 1)).sum()) for a in ARMS]
                            out.write('\t'.join(map(str, row)) + '\n')
        out.flush()
out.close()
print('done', flush=True)
