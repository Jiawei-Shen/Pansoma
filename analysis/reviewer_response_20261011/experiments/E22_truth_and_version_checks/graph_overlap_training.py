"""E22 check 3, training effect: do the truth SNVs that are HPRC v2.1 d46 alleles become somatic training examples (R1.m1)?

Joins the per-allele graph overlap (WORK graph_overlap.alleles.tsv.gz, graph_overlap.py) with each tensor set's per-allele
candidate status (tensors/somatic.recall.tsv: tensor_representative = a tensor whose representative allele is the truth
allele, i.e. a somatic training example; no_candidate; filtered; tensor_non_representative_allele) by truth_id.
Scope: truth SNVs and INDELs (separately), FILTER PASS, in the benchmark BED, chr1-22 (the CASTLE lines are trained on every autosome).
Output: graph_overlap_training.tsv (this folder).
"""
import csv, gzip, os

HERE = os.path.dirname(os.path.abspath(__file__))
W = '/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks/graph_overlap.alleles.tsv.gz'
SETS = '/scratch/jshen/data/pansoma_v2_tensors'
LINES = {'HCC1954': 'HCC1954T', 'H1437': 'H1437T', 'H2009': 'H2009T', 'HCC1395': 'HCC1395T'}
allcar = {'SNV': {}, 'INDEL': {}}
with gzip.open(W, 'rt') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['sample'] in LINES and r['stratum'] in allcar and r['bed'] == '1':
            allcar[r['stratum']][(r['sample'], r['truth_id'])] = int(r['carriers'])
out = []
for (stratum, car), (s, t) in ((x, y) for x in allcar.items() for y in LINES.items()):
    for plat in ('Illumina', 'PacBio', 'ONT'):
        st = {}
        with open(f'{SETS}/{t}_{plat}_d46/tensors/somatic.recall.tsv') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                st[r['truth_id']] = r['status']
        for grp in ('in_d46', 'not_in_d46'):
            ids = [i for (ss, i), c in car.items() if ss == s and (c >= 46) == (grp == 'in_d46')]
            n = len(ids)
            tens = sum(st.get(i) == 'tensor_representative' for i in ids)
            nocand = sum(st.get(i) == 'no_candidate' for i in ids)
            out.append(dict(stratum=stratum, sample=s, platform=plat, group=grp, n=n, training_example=tens,
                            pct_training_example=round(100 * tens / n, 2) if n else '',
                            no_candidate=nocand, pct_no_candidate=round(100 * nocand / n, 2) if n else ''))
        tot = sum(st.get(i) == 'tensor_representative' for (ss, i), c in car.items() if ss == s)
        d46 = sum(st.get(i) == 'tensor_representative' for (ss, i), c in car.items() if ss == s and c >= 46)
        out.append(dict(stratum=stratum, sample=s, platform=plat, group='share_of_somatic_examples_that_are_d46_alleles', n=tot,
                        training_example=d46, pct_training_example=round(100 * d46 / tot, 2) if tot else '',
                        no_candidate='', pct_no_candidate=''))
with open(f'{HERE}/graph_overlap_training.tsv', 'w') as fo:
    w = csv.DictWriter(fo, list(out[0]), delimiter='\t')
    w.writeheader()
    w.writerows(out)
for r in out:
    print('\t'.join(map(str, r.values())))
