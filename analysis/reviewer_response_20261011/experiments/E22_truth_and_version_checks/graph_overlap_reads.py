"""E22 check 3, examples: tumor and normal PacBio HiFi read counts at the graph-overlap examples (R1.m1).

For every row of graph_overlap_examples.tsv: REF / ALT / other base counts at the SNV in the CASTLE tumor BAM and in the
matched-normal BAM, both read from the public CASTLE bucket over https with their .bai (only the reads over the site are
fetched; the local tumor copies in /scratch/jshen/data/CASTLE have no index). Filters as in DeepSomatic's pileups: primary, non-duplicate, MAPQ >= 5, base quality >= 10.
Output: graph_overlap_examples_reads.tsv (this folder).
Run from a scratch directory: samtools caches the remote .bai in the current directory.
"""
import csv, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SAMTOOLS = '/wanglab/jshen/anaconda3/envs/sambcfenv/bin/samtools'
NUM = {'HCC1954': '1954', 'H1437': '1437', 'H2009': '2009'}
BUCKET = 'https://storage.googleapis.com/brain-genomics-public/publications/park2024_deepsomatic/bams/HiFi/pbmm2_grch38_bams'
TUMOR = BUCKET + '/{n}_Tumor_HiFi.GRCh38.sorted.bam'
NORMAL = BUCKET + '/{n}_Normal_HiFi.GRCh38.sorted.bam'


def counts(bam, chrom, pos, ref, alt):
    out = subprocess.run([SAMTOOLS, 'mpileup', '-r', f'{chrom}:{pos}-{pos}', '-B', '-q', '5', '-Q', '10',
                          '--ff', 'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY', bam],
                         capture_output=True, text=True, check=True).stdout.strip()
    if not out:
        return 0, 0, 0
    bases = out.split('\t')[4]
    r = a = o = 0
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == '^':
            i += 2
            continue
        if c in '+-':
            j = i + 1
            while j < len(bases) and bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j])
            continue
        if c in '.,' or c.upper() == ref:  # no -f reference: matches are printed as bases
            r += 1
        elif c.upper() == alt:
            a += 1
        elif c.upper() in 'ACGTN*':
            o += 1
        i += 1
    return r, a, o


rows = list(csv.DictReader(open(f'{HERE}/graph_overlap_examples.tsv'), delimiter='\t'))
out = []
for x in rows:
    n = NUM[x['sample']]
    tr, ta, to = counts(TUMOR.format(n=n), x['chrom'], x['pos'], x['ref'], x['alt'])
    nr, na, no = counts(NORMAL.format(n=n), x['chrom'], x['pos'], x['ref'], x['alt'])
    out.append({**x, 'tumor_ref': tr, 'tumor_alt': ta, 'tumor_other': to, 'tumor_alt_af': round(ta / max(1, tr + ta + to), 3),
                'normal_ref': nr, 'normal_alt': na, 'normal_other': no,
                'normal_alt_af': round(na / max(1, nr + na + no), 3)})
    print(x['sample'], x['chrom'], x['pos'], 'tumor ref/alt/other', tr, ta, to, 'normal', nr, na, no, flush=True)
with open(f'{HERE}/graph_overlap_examples_reads.tsv', 'w') as fo:
    w = csv.DictWriter(fo, list(out[0]), delimiter='\t')
    w.writeheader()
    w.writerows(out)
