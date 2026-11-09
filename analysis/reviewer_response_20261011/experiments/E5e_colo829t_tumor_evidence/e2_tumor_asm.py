"""E5e step 2: COLO829T tumor assembly at the normal-absent absorbed alleles (d46).

The COLO829T tumor assembly is hifiasm of the tumor's own HiFi + ONT reads (hap1 / hap2, E9c asm_targets), i.e. built
from the same tumor reads the truth set was called from: agreement with the truth shows consistency, not independent
support (R3.7). Method = E27 v3 colo829t_altlen_tumor_asm/altlen_tumor_asm.py (R3.7 method) on every locus instead of 22:
the c3 normal-frame anchors (24-mers ending at anchor_a0 / starting at anchor_b0; whole tandem array between common
unique anchors, chosen on the COLO829BL haplotypes) are searched exactly in every tumor contig, both strands; each hit
T = the tumor bases between them (oriented to GRCh38, at most len(R) + 2 kb). Compared with R = GRCh38[a, b),
A = R with the truth allele, and the COLO829BL hap sequences of c3 (assembly, else dipcall; '' is a valid sequence when
the hap deletes the whole stretch).
Input: $DATA/loci.tsv (e1), sets absent_INDEL and absent_SNV inside normal_frame.tsv.
Output: $DATA/tumor_asm/hits.<chunk>.tsv (one row per hit, a 'none' row when a locus has no hit in a tumor hap set).
Run: python e2_tumor_asm.py CHUNK NCHUNKS   (Slurm array e2_tumor_asm.sbatch; the 22-locus run: 7.5 min, 2.8 GB)
"""
import csv, gzip, sys

csv.field_size_limit(10 ** 9)
DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E5e_colo829t_tumor_evidence'
NF = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005/normal_frame.tsv'
FA = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
TUM = {'t_hap1': '/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation/asm_targets/COLO829T_hifiasm_hap1.fa.gz',
       't_hap2': '/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation/asm_targets/COLO829T_hifiasm_hap2.fa.gz'}
K = 24
COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path, want=None):
    seqs, name, buf = {}, None, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None and (want is None or name in want):
                    seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (want is None or name in want):
        seqs[name] = ''.join(buf).upper()
    return seqs


def main(chunk, nchunks):
    L = [r for r in csv.DictReader(open(f'{DATA}/loci.tsv'), delimiter='\t')
         if r['set'] in ('absent_INDEL', 'absent_SNV') and r['in_nf'] == 'True']
    L = L[chunk::nchunks]
    want = {r['truth_id'] for r in L}
    nf = {r['truth_id']: r for r in csv.DictReader(open(NF), delimiter='\t') if r['truth_id'] in want}
    g = read_fasta(FA, {r['chrom'] for r in L})
    loci = []
    for r in L:
        n = nf[r['truth_id']]
        a, b = int(n['anchor_a0']), int(n['anchor_b0'])
        s = g[n['chrom']]
        assert s[a:b] == n['R_seq'], (r['truth_id'], 'R_seq mismatch')
        loci.append(dict(set=r['set'], n=n, L=s[a - K:a], Rt=s[b:b + K], maxlen=len(n['R_seq']) + 2000))
    del g
    print('chunk', chunk, 'loci', len(loci), flush=True)
    found = {i: [] for i in range(len(loci))}
    for hap, path in TUM.items():
        contigs = read_fasta(path)
        print(hap, 'contigs', len(contigs), sum(map(len, contigs.values())), flush=True)
        for name, seq in contigs.items():
            for i, lo in enumerate(loci):
                for strand, left, right in (('+', lo['L'], lo['Rt']), ('-', rc(lo['Rt']), rc(lo['L']))):
                    p = seq.find(left)
                    while p != -1:
                        q = seq.find(right, p + K, p + K + lo['maxlen'] + K)
                        if q != -1:
                            mid = seq[p + K:q]
                            found[i].append((hap, name, p, strand, mid if strand == '+' else rc(mid)))
                        p = seq.find(left, p + 1)
        del contigs
    with open(f'{DATA}/tumor_asm/hits.{chunk}.tsv', 'w') as fo:
        fo.write('set\ttruth_id\tchrom\tvcf_pos\tvcf_ref\tvcf_alt\tkind2\talt_len_haps\tpf_class\tlen_R\tlen_A\tlen_hapX\t'
                 'len_hapY\thapX_allele\thapY_allele\ttumor_hap\tcontig\tpos0\tstrand\tlen_T\tT_eq_A\tT_eq_R\tT_eq_hapX\t'
                 'T_eq_hapY\tT_seq\n')
        for i, lo in enumerate(loci):
            r = lo['n']
            A, R, X, Y = r['A_seq'], r['R_seq'], r['hapX_seq'], r['hapY_seq']
            base = [lo['set'], r['truth_id'], r['chrom'], r['vcf_pos'], r['vcf_ref'], r['vcf_alt'], r['kind2'],
                    r['alt_len_haps'], r['pf_class'], len(R), len(A), len(X), len(Y), r['hapX_allele'], r['hapY_allele']]
            hits = found[i] or [('none', '', '', '', None)]
            for hap, name, p, strand, T in hits:
                if T is None:
                    fo.write('\t'.join(map(str, base + [hap] + [''] * 9)) + '\n')
                    continue
                fo.write('\t'.join(map(str, base + [hap, name, p, strand, len(T), T == A, T == R, T == X, T == Y, T])) + '\n')
    print('done', flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
