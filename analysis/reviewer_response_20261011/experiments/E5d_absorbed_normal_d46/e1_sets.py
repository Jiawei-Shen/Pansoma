"""E5d step 1: the absorbed truth alleles of HG008T and COLO829T on HPRC v2.1 (d46) and, by the same rule, on HPRC v1.1
(d9), with the graph path that spells each ALT.

Set rule (E5b s5 'absorbed', the set R1.4 quotes): on at least one platform the truth allele has no candidate on the arm
(status N), a tensor on the GRCh38-only arm (status T), and a path of the arm's graph spells GRCh38 + the truth allele
between the window anchors (E5b s4 'spelled'). Statuses and spelling flags are read from E5b's per_truth.tsv.gz; nothing
is re-decided here. Platforms: HG008T Illumina / PacBio / ONT, COLO829T Illumina / fiberseq / ONT.
The node path is re-derived with E5b's own spells() on its subgraph (same inputs, so the same answer; the script asserts
that every set allele is spelled again).
Outputs ($DATA): sets.tsv (sample, truth_id, kind, arm, platforms absorbed, PoN tag, chr1 / BED flags, element summary,
alt_min_paths, node path); and in this folder colo829t_d46_absorbed_loci.tsv (for E20) and hg008t_d46_absorbed_loci.tsv.
Run: python e1_sets.py   (Slurm: tmp/reviewer_cpu_20261021/E5d/e1_sets.sbatch; loads both reference-path indexes)
"""
import csv, gzip, sys
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'E5b_absorbed_d46'))
from e5b_common import DATA as E5B, FASTA, RefPath, SubGraph, spells, elements  # noqa: E402

OUT = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46'
PLATFORMS = {'HG008T': ('Illumina', 'PacBio', 'ONT'), 'COLO829T': ('Illumina', 'fiberseq', 'ONT')}


def main():
    rows = [r for r in csv.DictReader(gzip.open(f'{E5B}/per_truth.tsv.gz', 'rt'), delimiter='\t') if r['sample'] in PLATFORMS]
    loci = {(r['sample'], r['truth_id']): r for r in csv.DictReader(gzip.open(f'{E5B}/loci.tsv.gz', 'rt'), delimiter='\t')}
    fa = pysam.FastaFile(FASTA)
    out = []
    for arm in ('d9', 'd46'):
        R = RefPath(arm)
        print(arm, 'reference path loaded', file=sys.stderr, flush=True)
        G = SubGraph(f'{E5B}/subgraph_{arm}.gfa')
        for r in rows:
            plats = [p for p in PLATFORMS[r['sample']]
                     if r[f'{p}_{arm}'] == 'N' and r[f'{p}_linear'] == 'T' and r[f'{arm}_spelled'] == '1']
            if not plats:
                continue
            L = loci[(r['sample'], r['truth_id'])]
            chrom, pos0, vref, valt = r['chrom'], int(r['vcf_pos']) - 1, r['vcf_ref'].upper(), r['vcf_alt'].upper()
            lo, hi = int(L[f'{arm}_win_start0']), int(L[f'{arm}_win_end0'])
            region = fa.fetch(chrom, lo, hi + 1).upper()
            off = pos0 - lo
            assert region[off:off + len(vref)] == vref
            alt = region[:off] + valt + region[off + len(vref):]
            ok, path, _ = spells(G, R, chrom, int(L[f'{arm}_first']), int(L[f'{arm}_last']), lo, hi, alt)
            assert ok, (r['sample'], r['truth_id'], arm)
            out.append(dict(sample=r['sample'], truth_id=r['truth_id'], chrom=chrom, vcf_pos=r['vcf_pos'], vcf_ref=r['vcf_ref'],
                            vcf_alt=r['vcf_alt'], kind=r['kind'], kind2='SNV' if r['kind'] == 'SNP' else 'INDEL', arm=arm,
                            platforms=','.join(plats), pon=r['pon'], in_bed=r['in_bed'], passed=r['passed'],
                            absorbed_d9_reads=r['absorbed_d9_reads'],
                            elements=';'.join(':'.join(map(str, e)) for e in elements(R, path)),
                            alt_min_paths=r[f'{arm}_alt_min_paths'], win_start0=lo, win_end0=hi,
                            node_path=','.join(f'{n}{o}' for n, o in path)))
        print(arm, 'done', sum(o['arm'] == arm for o in out), file=sys.stderr, flush=True)
    cols = list(out[0])
    with open(f'{OUT}/sets.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for o in out:
            w.write('\t'.join(str(o[c]) for c in cols) + '\n')
    for sample, name in (('COLO829T', 'colo829t'), ('HG008T', 'hg008t')):
        keep = ['chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'truth_id', 'platforms', 'pon', 'elements', 'node_path']
        with open(HERE / f'{name}_d46_absorbed_loci.tsv', 'w') as w:
            w.write('\t'.join(keep) + '\n')
            for o in out:
                if o['sample'] == sample and o['arm'] == 'd46':
                    w.write('\t'.join(str(o[c]) for c in keep) + '\n')


if __name__ == '__main__':
    main()
