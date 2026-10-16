"""List the vcf directories E9c can process: finished (eval_nobed/report.json), one *.linear.unplaced.vcf.gz, on the
d9 or d46 graph (the GRCh38-only graph has no off-reference nodes).

Sources:
  rcpu    /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1[_<SET>] (the RIS 7fc7396 models,
          CPU re-predictions of the main session; read only)
  legacy  the local pre-RIS runs of the same recipes (SNV pow075_S_lr1e4_e20 of 2026-09-30/10-01, INDEL small_nofilt
          of 2026-10-10), available now
usage: discover.py [--all]   prints name<TAB>vcfdir, only those without DATA/runs/<name>/calls.tsv unless --all
"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e9c_lib as L  # noqa: E402

LEGACY = [f'{p}_{k}' for p in ('HG008_Illumina', 'HG008_PacBio', 'HG008_ONT', 'COLO829T_Illumina', 'COLO829T_fiberseq',
                                'COLO829T_ONT') for k in ('SNV_pow075_S_lr1e4_e20', 'INDEL_small_nofilt')]


def candidates():
    for d in sorted(glob.glob(f'{L.RUNS}/reviewer_cpu_20261013/*/vcf_chr1*')):
        if os.path.isdir(d):
            run = d.split('/')[-2]
            yield f'rcpu__{run}__{os.path.basename(d)}', d
    for run in LEGACY:
        d = f'{L.RUNS}/{run}/vcf_chr1'
        if os.path.isdir(d):
            yield f'legacy__{run}__vcf_chr1', d


def ready(d):
    if not os.path.exists(f'{d}/eval_nobed/report.json'):
        return False
    u = glob.glob(f'{d}/*.linear.unplaced.vcf.gz')
    if len(u) != 1:
        return False
    return L.graph_of(L.vcf_header(u[0]).get('pansoma_graph_index', '')) in L.GRAPHS


def main():
    everything = '--all' in sys.argv
    for name, d in candidates():
        if ready(d) and (everything or not os.path.exists(f'{L.DATA}/runs/{name}/calls.tsv')):
            print(f'{name}\t{d}')


if __name__ == '__main__':
    main()
