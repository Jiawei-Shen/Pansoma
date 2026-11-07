"""E5e step 1: COLO829T loci for the tumor-assembly and tumor-read checks.

Absorbed set = E5d's COLO829T HPRC v2.1 (d46) set (401 = 54 SNV + 347 INDEL) with E5d's per-allele normal category
(E5d per_allele.tsv: absent 319 INDEL + 52 SNV). Patient frame (anchors, R / A / normal hap sequences, event hap, the
normal-frame event length) = analysis/graph_absorbed_somatic_colo829t_20261005 c3 normal_frame.tsv (every COLO829T
chr1-22 truth INDEL and the 67 SNVs of 2026-10-05; 3 of the 52 absent SNVs are outside it).
Controls = 400 random COLO829T truth INDELs of normal_frame.tsv that are in neither the d9 absorbed set (normal_frame
'absorbed') nor the d46 set (random.Random(9).sample over the truth_id-sorted pool; the HG008 audit used seed 9 too).
Event hap for the audit's b2 window code: dipcall of COLO829BL was run as 'run-dipcall GRCh38 hapY.fa hapX.fa', so
GT = hapY|hapX: hapY -> 'hap1' (GT index 0), hapX -> 'hap2' (GT index 1).
Output: $DATA/loci.tsv
Run: python e1_loci.py   (login node, seconds)
"""
import csv, random

csv.field_size_limit(10 ** 9)
DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E5e_colo829t_tumor_evidence'
E5D = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46/per_allele.tsv'
NF = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005/normal_frame.tsv'
GT_HAP = {'hapY': 'hap1', 'hapX': 'hap2'}


def main():
    nf = {r['truth_id']: r for r in csv.DictReader(open(NF), delimiter='\t')}
    pa = [r for r in csv.DictReader(open(E5D), delimiter='\t') if r['sample'] == 'COLO829T' and r['arm'] == 'd46']
    assert len(pa) == 401, len(pa)
    d46 = {r['truth_id'] for r in pa}
    rows = []
    for r in pa:
        if r['category'] != 'absent':
            continue
        n = nf.get(r['truth_id'])
        rows.append(dict(set=f"absent_{r['kind2']}", truth_id=r['truth_id'], chrom=r['chrom'], vcf_pos=r['vcf_pos'],
                         vcf_ref=r['vcf_ref'], vcf_alt=r['vcf_alt'], kind2=r['kind2'], category=r['category'],
                         in_nf=n is not None, event_hap_nf=n['event_hap'] if n else '',
                         event_hap=GT_HAP.get(n['event_hap'], '') if n else '',
                         pf_len_change=n['pf_len_change'] if n else '', pf_class=n['pf_class'] if n else ''))
    pool = sorted((r for t, r in nf.items() if r['kind2'] == 'INDEL' and r['absorbed'] == 'False' and t not in d46),
                  key=lambda r: int(r['truth_id']))
    for n in random.Random(9).sample(pool, 400):
        rows.append(dict(set='control_INDEL', truth_id=n['truth_id'], chrom=n['chrom'], vcf_pos=n['vcf_pos'],
                         vcf_ref=n['vcf_ref'], vcf_alt=n['vcf_alt'], kind2='INDEL', category='control', in_nf=True,
                         event_hap_nf=n['event_hap'], event_hap=GT_HAP.get(n['event_hap'], ''),
                         pf_len_change=n['pf_len_change'], pf_class=n['pf_class']))
    cols = list(rows[0])
    with open(f'{DATA}/loci.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for r in rows:
            w.write('\t'.join(str(r[k]) for k in cols) + '\n')
    import collections
    print(collections.Counter((r['set'], r['in_nf']) for r in rows), 'control pool', len(pool))


if __name__ == '__main__':
    main()
