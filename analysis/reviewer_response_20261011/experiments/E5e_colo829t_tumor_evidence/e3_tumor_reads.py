"""E5e step 3: COLO829T tumor read alleles at the normal-absent absorbed INDELs (d46) and at 400 random controls.

The HG008T numbers (E5d e3) use the 2026-10-01 biology audit's b2_tumor_reads.py work() and the b4 rule. The same code is
used here, imported unchanged, with three module paths set to the COLO829 inputs:
  DIP  = COLO829BL dipcall vs GRCh38 (/scratch/qfu/COLO829BL_DSA/dipcall_hg38), cleaned as the HG008-N copy was
         (ACGTN alleles only; bcftools 1.18 norm -f GRCh38 no-alt -m -any; sort; index): $DATA/COLO829BL_dipcall.clean.norm.sorted.vcf.gz
  BAMS = the linear GRCh38 COLO829T tumor BAMs: Illumina (SMaHT sentieon bwa-mem), Fiber-seq (PacBio Revio, SMaHT pbmm2,
         merged), ONT (SMaHT minimap2 2.26 map-ont, merged).
  event hap = the c3 normal-frame event hap (hapY -> dipcall GT index 0 'hap1', hapX -> index 1 'hap2'; e1).
work(): window = the longest period-1..6 tract touching the event +-5 bp, widened to every overlapping PASS dipcall record
and the truth REF; reads with MAPQ >= 20, primary, not dup / qcfail, one per Illumina pair, <= 300 per locus; read allele =
bases between the window anchors; NE / NO lengths from the dipcall haplotypes.
b4 (E5d e3 b4_row, imported): EXP = NE_len + nf_len, nf_len = c3 pf_len_change (the normal-frame event length, as b5's
nf_len for HG008T); EXP 'seen' at >= max(3, 15 %) of the reads; loci with >= 10 reads and nf_len != 0 only.
Outputs ($DATA/tumor_reads/): b2.tsv (one row per locus x platform), b4.tsv.
Run: python e3_tumor_reads.py [workers]   (Slurm: e3_tumor_reads.sbatch)
"""
import csv, importlib.util, multiprocessing as mp, sys

DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E5e_colo829t_tumor_evidence'
BIO = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/audit/biology'
E5D = '/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E5d_absorbed_normal_d46'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


b2 = load('b2', f'{BIO}/b2_tumor_reads.py')
b2.DIP = f'{DATA}/COLO829BL_dipcall.clean.norm.sorted.vcf.gz'
b2.BAMS = {'Illumina': '/scratch/qfu/COLO829T_WGS/SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIZO17Q9I-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.bam',
           'Fiberseq': '/scratch/jshen/data/COLO829T/COLO829T_fiberseq/SMHTCOLO829T-X-X-M45-B004-uwsc-MERGED-pbmm2_1.13.0_GRCh38.aligned.sorted.bam',
           'ONT': '/scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.bam'}
e3 = load('e3_hg008t', f'{E5D}/e3_hg008t_tumor_reads.py')   # for b4_row only (its b2 import is the same audit file)


def work(t):
    return b2.work(t)


def main(workers):
    loci = [r for r in csv.DictReader(open(f'{DATA}/loci.tsv'), delimiter='\t')
            if r['set'] in ('absent_INDEL', 'control_INDEL')]
    todo = [dict(r, set='absorbed' if r['set'] == 'absent_INDEL' else 'control') for r in loci]
    with mp.Pool(workers) as pool:
        res = pool.map(work, todo, chunksize=2)
    rows = [o for rr in res for o in rr]
    cols = list(rows[0])
    with open(f'{DATA}/tumor_reads/b2.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for o in rows:
            w.write('\t'.join(str(o[k]) for k in cols) + '\n')
    b5 = {r['truth_id']: dict(nf_len=r['pf_len_change']) for r in loci}
    out = [e3.b4_row({k: str(v) for k, v in o.items()}, b5) for o in rows]
    cols = ['truth_id', 'platform', 'set', 'kind2', 'n', 'status', 'EXP', 'f_EXP', 'EXP_seen']
    with open(f'{DATA}/tumor_reads/b4.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for o in out:
            w.write('\t'.join(str(o.get(k, '')) for k in cols) + '\n')
    print('loci', len(loci), 'rows', len(out))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
