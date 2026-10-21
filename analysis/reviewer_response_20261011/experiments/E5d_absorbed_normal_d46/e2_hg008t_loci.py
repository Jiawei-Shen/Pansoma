"""E5d step 2 (HG008T): loci.tsv for the set alleles (sets.tsv, d9 or d46) that the 2026-10-01 analysis
(analysis/graph_absorbed_somatic_20261001, 2,478 alleles) did not cover, in that analysis' s1 format, so its s6 / s8b
steps (copied here as e2_s6_assembly_seq.py / e2_s8b_array_alleles.py, only the data directory changed) run on them.

Window = E5b's d9 window (s2_windows.py), which equals s1's window at all 2,478 old loci (first / last node, window and
event window checked identical); ref_hap / alt_hap from the GRCh38 FASTA; asm_event = the truth INFO
HG008Nv62SOMATICVARIANT (first value, as s1). Outputs ($DATA/hg008t/): variant_set.tsv, loci.tsv.
Run: python e2_hg008t_loci.py   (login node, seconds)
"""
import csv, gzip
from pathlib import Path

import pysam

E5B = '/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46'
OLD = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
OUT = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46'
FASTA = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
TRUTH = '/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz'
COLS = ['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'status', 'first_node', 'last_node', 'n_ref_nodes',
        'win_start0', 'win_end0', 'ev_lo0', 'ev_hi0', 'asm_event', 'ref_hap', 'alt_hap']


def main():
    old = {r['truth_id'] for r in csv.DictReader(open(f'{OLD}/variant_set.tsv'), delimiter='\t')}
    sets = [r for r in csv.DictReader(open(f'{OUT}/sets.tsv'), delimiter='\t') if r['sample'] == 'HG008T']
    new = {}
    for r in sets:
        if r['truth_id'] not in old:
            new[r['truth_id']] = r
    loci = {r['truth_id']: r for r in csv.DictReader(gzip.open(f'{E5B}/loci.tsv.gz', 'rt'), delimiter='\t') if r['sample'] == 'HG008T'}
    keys = {(r['chrom'], int(r['vcf_pos']), r['vcf_ref'], r['vcf_alt']) for r in new.values()}
    asm = {}
    for rec in pysam.VariantFile(TRUTH):
        for a in rec.alts or ():
            if (rec.chrom, rec.pos, rec.ref, a) in keys:
                v = rec.info.get('HG008Nv62SOMATICVARIANT')
                asm[(rec.chrom, rec.pos, rec.ref, a)] = (v[0] if isinstance(v, tuple) else v) or ''
    fa = pysam.FastaFile(FASTA)
    d = Path(OUT) / 'hg008t'
    d.mkdir(exist_ok=True)
    ids = sorted(new, key=int)
    with open(d / 'variant_set.tsv', 'w') as w:
        w.write('\t'.join(['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'kind2', 'arms']) + '\n')
        for t in ids:
            r = new[t]
            arms = ','.join(sorted({s['arm'] for s in sets if s['truth_id'] == t}))
            w.write('\t'.join([t, r['chrom'], r['vcf_pos'], r['vcf_ref'], r['vcf_alt'], r['kind'], r['kind2'], arms]) + '\n')
    n_bad = 0
    with open(d / 'loci.tsv', 'w') as w:
        w.write('\t'.join(COLS) + '\n')
        for t in ids:
            r, L = new[t], loci[t]
            chrom, pos1, vref, valt = r['chrom'], int(r['vcf_pos']), r['vcf_ref'].upper(), r['vcf_alt'].upper()
            g = 'd9' if L['d9_status'] == 'ok' else 'd46'      # GRCh38 window; s6 / s8b use only its coordinates
            out = dict(truth_id=t, chrom=chrom, vcf_pos=pos1, vcf_ref=r['vcf_ref'], vcf_alt=r['vcf_alt'], kind2=r['kind2'],
                       status=L[f'{g}_status'], ev_lo0=L['ev_lo0'], ev_hi0=L['ev_hi0'],
                       asm_event=asm.get((chrom, pos1, r['vcf_ref'], r['vcf_alt']), ''))
            if L[f'{g}_status'] == 'ok':
                lo, hi = int(L[f'{g}_win_start0']), int(L[f'{g}_win_end0'])
                region = fa.fetch(chrom, lo, hi + 1).upper()
                off = pos1 - 1 - lo
                assert region[off:off + len(vref)] == vref
                out.update(first_node=L[f'{g}_first'], last_node=L[f'{g}_last'], n_ref_nodes=L[f'{g}_n_ref_nodes'],
                           win_start0=lo, win_end0=hi, ref_hap=region, alt_hap=region[:off] + valt + region[off + len(vref):])
            else:
                n_bad += 1
            w.write('\t'.join(str(out.get(c, '')) for c in COLS) + '\n')
    print(len(ids), 'new HG008T loci', n_bad, 'without a d9 window')


if __name__ == '__main__':
    main()
