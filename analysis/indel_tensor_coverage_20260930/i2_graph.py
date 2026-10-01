"""Graph view of the I2 somatic INDEL truths of HG008T Illumina (no germline information used).

For each I2 truth: the graph alleles at the site (sequences spelled by the d9 graph between the GRCh38 nodes ~60 bp
around it; each as a length change vs GRCh38), the truth's length change Lt, and the reads' main residual INDEL Lres
(analysis/hg008_somatic_miss_20260926/illumina/indel_no_candidate.tsv). The graph allele the reads used is taken as
Lg = Lt - Lres when the graph has an allele of that length. Classes:
  graph allele shorter, same direction - |Lg| < |Lt|: the graph holds part of the length change, the residual adds the rest
  graph allele longer, same direction  - |Lg| > |Lt|: the reads take a longer graph allele and the residual goes back
  graph allele of the other direction  - e.g. truth DEL, reads take a graph INS branch
  no graph allele of length Lt - Lres  - the residual does not complete any graph allele to the truth
  no residual INDEL                    - reads show no INDEL edit; nearest graph allele given
Output: i2_classes[_<chroms>].tsv
"""
import csv, re, sys
from collections import Counter
import pysam
src = open('/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/classify.py').read().split("fa = pysam.FastaFile(FASTA)\nchroms")[0]
ns = {}; exec(src, ns)
db, graph_haplotypes = ns['db'], ns['graph_haplotypes']
fa = pysam.FastaFile(ns['FASTA'])
MISS = '/scratch/jshen/Github/Pansoma/analysis/hg008_somatic_miss_20260926/illumina/indel_no_candidate.tsv'


def top_indel(s):
    best = None
    for e in s.split('; '):
        m = re.match(r'(\S+?)\((\d+)\)=(.*)', e)
        if not m:
            continue
        p = m.group(1).split(':')
        if p[2] not in ('INS', 'DEL'):
            continue
        ra = p[3].split('>')
        L = len(ra[1].split('@')[0]) if p[2] == 'INS' else -len(ra[0])
        if best is None or int(m.group(2)) > best[1]:
            best = (L, int(m.group(2)), m.group(1))
    return best


chroms = set(sys.argv[1:])
rows = [r for r in csv.DictReader(open(MISS), delimiter='\t') if r['final_class'].startswith('I2') and (not chroms or r['chrom'] in chroms)]
out = open(f'/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/i2_classes{"_" + "_".join(sorted(chroms)) if chroms else ""}.tsv', 'w')
out.write('chrom\tvcf_pos\tvcf_ref\tvcf_alt\tin_bed\treads_reason\ttruth_len\tgraph_alleles\tresidual_len\tresidual_reads\tspanning\tgraph_len_used\tclass\n')
for k, r in enumerate(rows):
    chrom, pos1, vref, valt = r['chrom'], int(r['vcf_pos']), r['vcf_ref'], r['vcf_alt']
    lt = len(valt) - len(vref)
    nodes = ns['window_nodes'](chrom, pos1 - 61, pos1 - 1 + len(vref) + 60)   # GRCh38 coordinates from the pipeline's grch38_path
    if not nodes:
        out.write(f"{chrom}\t{pos1}\t{vref}\t{valt}\t{r['in_bed']}\t{r['reason']}\t{lt}\t\t\t\t{r['spanning']}\t\tno GRCh38 window (reverse node)\n")
        continue
    first, last = nodes[0], nodes[-1]
    region = fa.fetch(chrom, first[1], last[2] + 1).upper()
    off = pos1 - 1 - first[1]
    alt_hap = region[:off] + valt.upper() + region[off + len(vref):]
    haps = graph_haplotypes(chrom, first, last, limit=3000, max_len=len(region) + 300)
    single = sorted({len(h) - len(region) for h in haps if h != region and ns['lev'](h, region) == abs(len(h) - len(region))}) \
        if len(haps) <= 300 else []                                   # graph alleles that are one pure INDEL vs GRCh38
    lens = single
    t = top_indel(r['residual_edits(reads)=tensor_label'])
    lg, cls = None, None
    if alt_hap in haps:
        cls = 'truth haplotype IS a graph path'
    elif t is None:
        near = sorted(haps, key=lambda h: abs(len(h) - len(alt_hap)))[:300]
        d, h = min((ns['lev'](h, alt_hap), h) for h in near)
        lg = len(h) - len(region)
        cls = 'no residual INDEL: reads follow a graph allele that is not the truth' if h != region else 'no residual INDEL: closest is GRCh38'
    else:
        done = [h for h in haps if len(alt_hap) - len(h) == t[0] and ns['lev'](h, alt_hap) == abs(t[0])]
        if not done:
            cls = 'residual does not turn any graph allele into the truth'
        else:
            h = min(done, key=lambda h: ns['lev'](h, region))
            lg = len(h) - len(region)
            if h == region:
                cls = 'reads on GRCh38 + residual = truth (graph allele not used)'
            elif lg == 0:
                cls = 'graph allele same length (SNV-type branch) + residual = truth'
            elif (lg > 0) != (lt > 0):
                cls = 'graph allele of the other direction + residual = truth'
            elif abs(lg) < abs(lt):
                cls = 'graph allele shorter, same direction + residual = truth'
            else:
                cls = 'graph allele longer, same direction + residual = truth'
    out.write(f"{chrom}\t{pos1}\t{vref}\t{valt}\t{r['in_bed']}\t{r['reason']}\t{lt}\t{','.join(map(str, lens))}\t"
              f"{t[0] if t else ''}\t{t[1] if t else ''}\t{r['spanning']}\t{'' if lg is None else lg}\t{cls}\n")
    if k % 100 == 0:
        print(k, len(rows), flush=True)
out.close()
