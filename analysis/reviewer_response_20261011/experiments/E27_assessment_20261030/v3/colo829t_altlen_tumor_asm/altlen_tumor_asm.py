"""R2.3 v2.1 (d46) check (2026-11-03): the R3.7 tumor-assembly comparison (tmp/assessment_20261030/misc/r37_tumor_asm_35.py,
unchanged method) at the 22 absorbed d46 COLO829T truth INDELs where a COLO829BL haplotype has exactly the ALT length
with other bases (E5d section 4.4); loci and sequences from c3 normal_frame.tsv (every truth INDEL), ids in argv[2].
usage: altlen_tumor_asm.py OUT.tsv IDS.txt"""
import csv, gzip, sys

PV = '/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_colo829t_20261005/per_variant.tsv'
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


def lev(a, b):
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


csv.field_size_limit(10 ** 9)
IDS = [x.strip() for x in open(sys.argv[2]) if x.strip()]
NF = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005/normal_frame.tsv'
rows = [r for r in csv.DictReader(open(NF), delimiter='\t') if r['truth_id'] in IDS]
assert len(rows) == len(IDS), (len(rows), len(IDS))
chroms = {r['chrom'] for r in rows}
g = read_fasta(FA, chroms)
print('GRCh38 loaded', len(g), flush=True)
loci = []
for r in rows:
    a, b = int(r['anchor_a0']), int(r['anchor_b0'])
    s = g[r['chrom']]
    assert s[a:b] == r['R_seq'], (r['truth_id'], 'R_seq mismatch')
    loci.append(dict(r=r, L=s[a - K:a], Rt=s[b:b + K], maxlen=len(r['R_seq']) + 2000))
del g

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

with open(sys.argv[1], 'w') as fo:
    fo.write('truth_id\tchrom\tvcf_pos\tvcf_ref\tvcf_alt\talt_len_haps\tpf_class\tlen_R\tlen_A\tlen_hapX\tlen_hapY\t'
             'hapX_allele\thapY_allele\ttumor_hap\tcontig\tpos0\tstrand\tlen_T\tT_eq_A\tT_eq_R\tT_eq_hapX\tT_eq_hapY\t'
             'lev_T_A\tlev_T_hapX\tlev_T_hapY\n')
    for i, lo in enumerate(loci):
        r = lo['r']
        A, R, X, Y = r['A_seq'], r['R_seq'], r['hapX_seq'], r['hapY_seq']
        base = [r['truth_id'], r['chrom'], r['vcf_pos'], r['vcf_ref'], r['vcf_alt'], r['alt_len_haps'], r['pf_class'],
                len(R), len(A), len(X), len(Y), r['hapX_allele'], r['hapY_allele']]
        hits = found[i] or [('none', '', '', '', None)]
        for hap, name, p, strand, T in hits:
            if T is None:
                fo.write('\t'.join(map(str, base + [hap, '', '', '', '', '', '', '', '', '', '', ''])) + '\n')
                continue
            fo.write('\t'.join(map(str, base + [hap, name, p, strand, len(T), T == A, T == R, T == X, T == Y,
                                                 lev(T, A), lev(T, X) if X else '', lev(T, Y) if Y else ''])) + '\n')
print('done', flush=True)
