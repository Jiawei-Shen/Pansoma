"""E5b step 1: one row per somatic truth allele, its candidate status on every tensor set (platform x reference arm).

Inputs (read-only): /scratch/jshen/data/pansoma_v2_tensors/<SET>/tensors/somatic.recall.tsv (per truth allele: status
tensor_representative / tensor_non_representative_allele / filtered / no_candidate / ...) and <SET>/truth/somatic.graph.tsv
(the truth allele itself: chrom, vcf_pos, vcf_ref, vcf_alt). truth_id is the row of the truth table, identical on every
arm of a sample (checked below: same chrom / vcf_pos / kind on every set).
Arms: d9 = <SET> (HPRC v1.1 d9), d46 = <SET>_d46 (HPRC v2.1 d46), linear = <SET>_linear (GRCh38 only).
Output: status_matrix.tsv.gz in the data dir (columns sample, truth_id, chrom, vcf_pos, vcf_ref, vcf_alt, kind, passed, in_bed,
then <platform>_<arm> = status code). Status codes: T tensor_representative, A tensor_non_representative_allele,
F filtered, N no_candidate, other = first letters.
"""
import csv, gzip, os, sys
csv.field_size_limit(1 << 30)
T = '/scratch/jshen/data/pansoma_v2_tensors'
OUT = '/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46'
SAMPLES = {'HG008T': (['Illumina', 'PacBio', 'ONT'], ['d9', 'd46', 'linear']),
           'COLO829T': (['Illumina', 'fiberseq', 'ONT'], ['d9', 'd46', 'linear']),
           'HCC1395T': (['Illumina', 'PacBio', 'ONT'], ['d46'])}
CODE = {'tensor_representative': 'T', 'tensor_non_representative_allele': 'A', 'filtered': 'F', 'no_candidate': 'N',
        'complex_allele': 'C', 'no_unique_grch38_node_key': 'K'}


def set_name(sample, plat, arm):
    return f'{sample}_{plat}' + ('' if arm == 'd9' else f'_{arm}')


os.makedirs(OUT, exist_ok=True)
cols = ['sample', 'truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'passed', 'in_bed']
allcols = list(cols)
rows_out = []
for sample, (plats, arms) in SAMPLES.items():
    truth = None
    status = {}
    for plat in plats:
        for arm in arms:
            s = set_name(sample, plat, arm)
            g = list(csv.DictReader(open(f'{T}/{s}/truth/somatic.graph.tsv'), delimiter='\t'))
            r = list(csv.DictReader(open(f'{T}/{s}/tensors/somatic.recall.tsv'), delimiter='\t'))
            assert len(g) == len(r), s
            for a, b in zip(g, r):
                assert a['truth_id'] == b['truth_id'] and a['chrom'] == b['chrom'] and a['vcf_pos'] == b['vcf_pos'] and a['kind'] == b['kind'], s
            if truth is None:
                truth = [{c: a[c] for c in ('truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'passed', 'in_bed')} for a in g]
            else:
                assert [t['truth_id'] for t in truth] == [a['truth_id'] for a in g], s
                for t, a in zip(truth, g):
                    assert (t['chrom'], t['vcf_pos'], t['vcf_ref'], t['vcf_alt']) == (a['chrom'], a['vcf_pos'], a['vcf_ref'], a['vcf_alt']), s
            status[f'{plat}_{arm}'] = [CODE.get(b['status'], b['status'][:2]) for b in r]
            print(s, len(r), file=sys.stderr)
    for k in status:
        if k not in allcols:
            allcols.append(k)
    for i, t in enumerate(truth):
        d = dict(t, sample=sample)
        for k, v in status.items():
            d[k] = v[i]
        rows_out.append(d)
plat_cols = [c for c in allcols if c not in cols]
with gzip.open(f'{OUT}/status_matrix.tsv.gz', 'wt') as w:
    w.write('\t'.join(cols + plat_cols) + '\n')
    for d in rows_out:
        w.write('\t'.join(str(d.get(c, '')) for c in cols + plat_cols) + '\n')
print(len(rows_out), 'rows', file=sys.stderr)
