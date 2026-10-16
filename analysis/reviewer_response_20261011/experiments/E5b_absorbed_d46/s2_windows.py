"""E5b step 2: one GRCh38 window per somatic truth allele, and its anchor nodes on the d9 and d46 graphs.

Window (graph_absorbed_somatic_20261001 s1_loci.py): the GRCh38 nodes overlapping [pos0 - 60, pos0 + len(REF) + 60],
widened to the event window +-20 bp (the periodic stretch around the event); first / last of them = anchors.
ref_hap / alt_hap are the GRCh38 sequence of the window and the same with the truth allele applied; they depend on the
graph only through the anchors' extent. Writes loci.tsv.gz (one row per sample x truth_id, the d9 and d46 anchors and
spans) and window_nodes_<graph>.txt (every GRCh38 window node, the seeds of s3_extract.sh's `vg find -N`).
Run: python s2_windows.py
"""
import csv, gzip, sys
import pysam
from e5b_common import DATA, FASTA, RefPath, event_window

FLANK = 60
fa = pysam.FastaFile(FASTA)
refs = {g: RefPath(g) for g in ('d9', 'd46')}
seeds = {g: set() for g in refs}
cols = ['sample', 'truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'ev_lo0', 'ev_hi0']
for g in refs:
    cols += [f'{g}_status', f'{g}_first', f'{g}_last', f'{g}_win_start0', f'{g}_win_end0', f'{g}_n_ref_nodes']
n = 0
with gzip.open(f'{DATA}/status_matrix.tsv.gz', 'rt') as f, gzip.open(f'{DATA}/loci.tsv.gz', 'wt') as w:
    w.write('\t'.join(cols) + '\n')
    for r in csv.DictReader(f, delimiter='\t'):
        chrom, pos1, vref, valt = r['chrom'], int(r['vcf_pos']), r['vcf_ref'].upper(), r['vcf_alt'].upper()
        pos0 = pos1 - 1
        assert fa.fetch(chrom, pos0, pos0 + len(vref)).upper() == vref, (r['sample'], r['truth_id'])
        ev = event_window(fa, chrom, pos0, vref, valt)
        lo, hi = min(pos0 - FLANK, ev[0] - 20), max(pos0 + len(vref) + FLANK, ev[1] + 20)
        out = {c: r[c] for c in ('sample', 'truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind')}
        out.update(ev_lo0=ev[0], ev_hi0=ev[1])
        for g, R in refs.items():
            nodes = R.window_nodes(chrom, lo, hi)
            if not nodes:
                out[f'{g}_status'] = 'no_window'
                continue
            out.update({f'{g}_status': 'ok', f'{g}_first': nodes[0][0], f'{g}_last': nodes[-1][0],
                        f'{g}_win_start0': nodes[0][1], f'{g}_win_end0': nodes[-1][2], f'{g}_n_ref_nodes': len(nodes)})
            seeds[g].update(x[0] for x in nodes)
        w.write('\t'.join(str(out.get(c, '')) for c in cols) + '\n')
        n += 1
        if n % 10000 == 0:
            print(n, file=sys.stderr, flush=True)
for g, s in seeds.items():
    with open(f'{DATA}/window_nodes_{g}.txt', 'w') as w:
        w.write('\n'.join(map(str, sorted(s))) + '\n')
    print(g, len(s), 'window nodes', file=sys.stderr)
print(n, 'loci', file=sys.stderr)
