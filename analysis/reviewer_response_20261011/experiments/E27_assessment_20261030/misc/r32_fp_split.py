"""R3.2 fact check (2026-10-30, read-only): E4b's chr1 PASS false positives after population filtering (in the BED,
in-sample, main threshold) re-split by whether the FP allele is represented in the graph.

Input: E4b per-FP tables DATA/E4b/<call set>.eval_bed.fp.v2.tsv.gz (columns npos/nref/nalt = allele normalized by
E4b; hprc_ac = carriers among the 89 HPRC v1.1 haplotypes; class = E4b exclusive class). Added here: n_carriers among
the 463 HPRC v2.1 haplotypes from E7 WORK/v2.1/hprc_alleles/chr1.tsv.gz (same left-aligned primitive keys).
Represented: v1.1 any = hprc_ac > 0 (E4b's flag); v1.1 d9 = hprc_ac >= 9; v2.1 d46 = v2.1 carriers >= 46.
Output TSV to argv[1].
"""
import csv, gzip, glob, os, sys

E4B = '/scratch/jshen/data/reviewer_experiments_20261013/call_level/E4b'
V21 = '/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/v2.1/hprc_alleles/chr1.tsv.gz'
files = sorted(f for f in glob.glob(f'{E4B}/*__vcf_chr1.eval_bed.fp.v2.tsv.gz')
               if '__norescue' not in f and 'joint' not in f and not os.path.basename(f).startswith('xs__'))
data, keys = {}, set()
for f in files:
    rows = list(csv.DictReader(gzip.open(f, 'rt'), delimiter='\t'))
    data[f] = rows
    for r in rows:
        keys.add(f"chr1:{r['npos']}:{r['nref']}:{r['nalt']}")
v21 = {}
with gzip.open(V21, 'rt') as fh:
    fh.readline()
    for line in fh:
        k, rest = line.split('\t', 1)
        if k in keys:
            v21[k] = int(rest.split('\t')[5])
print('keys', len(keys), 'found in v2.1 table', len(v21), file=sys.stderr)
out = csv.writer(open(sys.argv[1], 'w'), delimiter='\t')
out.writerow(['call_set', 'fp', 'germline_same', 'gs_v11_any', 'gs_v11_ge9', 'gs_v21_any', 'gs_v21_ge46',
              'gs_not_v11_any', 'gs_not_v21_ge46', 'germline_overlap_or_near', 'hprc_class', 'population_classes',
              'repeat', 'other', 'all_fp_v11_any', 'all_fp_v21_ge46', 'other_v21_ge46'])
for f, rows in data.items():
    name = os.path.basename(f).replace('__vcf_chr1.eval_bed.fp.v2.tsv.gz', '')
    def c(pred):
        return sum(1 for r in rows if pred(r))
    ac = lambda r: int(r['hprc_ac'] or 0)
    nc = lambda r: v21.get(f"chr1:{r['npos']}:{r['nref']}:{r['nalt']}", 0)
    gs = lambda r: r['class'] == 'germline_same'
    out.writerow([name, len(rows), c(gs), c(lambda r: gs(r) and ac(r) > 0), c(lambda r: gs(r) and ac(r) >= 9),
                  c(lambda r: gs(r) and nc(r) > 0), c(lambda r: gs(r) and nc(r) >= 46),
                  c(lambda r: gs(r) and ac(r) == 0), c(lambda r: gs(r) and nc(r) < 46),
                  c(lambda r: r['class'] in ('germline_overlap', 'near_germline')), c(lambda r: r['class'] == 'hprc'),
                  c(lambda r: r['class'] in ('pop_below_floor', 'pop_other_repr')), c(lambda r: r['class'] == 'repeat'),
                  c(lambda r: r['class'] == 'other'), c(lambda r: ac(r) > 0), c(lambda r: nc(r) >= 46),
                  c(lambda r: r['class'] == 'other' and nc(r) >= 46)])
