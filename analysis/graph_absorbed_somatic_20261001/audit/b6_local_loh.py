"""Audit b6 (biology lens): is the normal event haplotype still present in the tumor around each locus (local LOH)?

inputs   b2_tumor_reads.tsv (loci + event hap: absorbed set and controls), b1_truth_context.tsv, the HG008-N dipcall VCF
         (normalised indexed copy), HG008-T Illumina 195x GRCh38 BAM, FASTA.
method   up to 8 nearest phased heterozygous dipcall SNVs (GT 0|1 / 1|0, PASS, no other dipcall record within 10 bp,
         not inside a homopolymer >= 4) within +-20 kb of the truth; tumor Illumina reads (MAPQ >= 20, base quality
         >= 20, primary, not dup, one read per pair) counted for the hap1 and hap2 base; fE = event-hap base reads /
         (event + other). States: E_lost fE < 0.1, O_lost fE > 0.9, both_retained 0.25-0.75, imbalanced otherwise;
         no_hets when < 2 informative SNVs or < 20 reads.
outputs  b6_local_loh.tsv: truth_id, set, category, event_hap, n_snv, reads_E, reads_O, fE, state, fE_per_snv.
assumes  dipcall GT order = v6.2 hap1|hap2 and phase holds over +-20 kb (same assembly contig).
"""
import csv, sys, multiprocessing as mp
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
O = f'{D}/audit/biology'
DIP = '/scratch/jshen/tmp/HG008N_vs_HPRC_isec/HG008N_dipcall.clean.norm.sorted.vcf.gz'
BAM = '/scratch/jshen/data/HG008_GIAB/HG008-T_Illumina_195x_GRCh38-GIABv3.bam'
FASTA = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
FL = 20000


def base_at(rd, pos):
    r = rd.reference_start; q = 0
    for op, l in rd.cigartuples:
        if op in (0, 7, 8):
            if r <= pos < r + l: return q + pos - r
            r += l; q += l
        elif op in (1, 4): q += l
        elif op in (2, 3):
            if r <= pos < r + l: return None
            r += l
        if r > pos: return None
    return None


def work(t):
    dip = pysam.VariantFile(DIP); bam = pysam.AlignmentFile(BAM); fa = pysam.FastaFile(FASTA)
    c, p0 = t['chrom'], int(t['vcf_pos']) - 1
    eh = {'hap1': 0, 'hap2': 1}.get(t['event_hap'])
    o = dict(truth_id=t['truth_id'], set=t['set'], category=t['category'], event_hap=t['event_hap'])
    if eh is None: return dict(o, state='no_event_hap')
    recs = [r for r in dip.fetch(c, max(0, p0 - FL), p0 + FL) if not r.filter.keys() or 'PASS' in r.filter.keys()]
    pos = [r.start for r in recs]
    cand = []
    for i, r in enumerate(recs):
        if len(r.ref) != 1 or len(r.alts[0]) != 1: continue
        gt = r.samples[0]['GT']
        if len(gt) != 2 or None in gt or sorted(gt) != [0, 1]: continue
        if (i > 0 and r.start - recs[i - 1].stop < 10) or (i + 1 < len(recs) and recs[i + 1].start - r.stop < 10): continue
        if abs(r.start - p0) < 50: continue
        ctx = fa.fetch(c, r.start - 3, r.start + 4).upper()
        if any(ctx[k:k + 4] == ctx[k] * 4 for k in range(4)): continue
        cand.append((abs(r.start - p0), r.start, r.ref, r.alts[0], gt))
    cand.sort()
    tot = [0, 0]; per = []
    for _, s, ref, alt, gt in cand[:8]:
        al = [ref, alt]; hb = (al[gt[0]], al[gt[1]]); cnt = [0, 0]; seen = set()
        for rd in bam.fetch(c, s, s + 1):
            if rd.is_secondary or rd.is_supplementary or rd.is_duplicate or rd.is_qcfail or rd.mapping_quality < 20: continue
            if rd.query_name in seen: continue
            q = base_at(rd, s)
            if q is None or rd.query_qualities[q] < 20: continue
            seen.add(rd.query_name); b = rd.query_sequence[q]
            if b == hb[0]: cnt[0] += 1
            elif b == hb[1]: cnt[1] += 1
        if sum(cnt) >= 10:
            e, oo = cnt[eh], cnt[1 - eh]; tot[0] += e; tot[1] += oo; per.append(round(e / (e + oo), 2))
    o.update(n_snv=len(per), reads_E=tot[0], reads_O=tot[1], fE_per_snv=','.join(map(str, per)))
    if len(per) < 2 or sum(tot) < 20: o['state'] = 'no_hets'; return o
    f = tot[0] / sum(tot); o['fE'] = round(f, 3)
    o['state'] = 'E_lost' if f < 0.1 else 'O_lost' if f > 0.9 else 'both_retained' if 0.25 <= f <= 0.75 else 'imbalanced'
    return o


if __name__ == '__main__':
    b1 = {r['truth_id']: r for r in csv.DictReader(open(f'{O}/b1_truth_context.tsv'), delimiter='\t')}
    todo = {}
    for r in csv.DictReader(open(f'{O}/b2_tumor_reads.tsv'), delimiter='\t'):
        if r['platform'] != 'PacBio': continue
        todo[r['truth_id']] = dict(b1[r['truth_id']], set=r['set'], category=r['category'], event_hap=r['event_hap'])
    todo = list(todo.values())
    if len(sys.argv) > 1: todo = todo[:int(sys.argv[1])]
    with mp.Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 16) as pool:
        res = pool.map(work, todo, chunksize=8)
    cols = ['truth_id', 'set', 'category', 'event_hap', 'n_snv', 'reads_E', 'reads_O', 'fE', 'state', 'fE_per_snv']
    with open(f'{O}/b6_local_loh.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for o in res: w.write('\t'.join(str(o.get(k, '')) for k in cols) + '\n')
    print('loci', len(res))
