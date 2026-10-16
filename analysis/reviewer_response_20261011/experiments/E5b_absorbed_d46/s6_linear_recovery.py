"""E5b step 6 (collect, needs the CPU re-predictions): which truth INDELs lost on the pangenome does the linear arm call
after the PoN, and which do d9 / d46 call?

Reads, for every HG008T / COLO829T INDEL run of the three arms, rtg vcfeval's tp-baseline.vcf.gz of the after-PoN
evaluation in the benchmark BED: <run>/vcf_chr1/eval_bed/vcfeval/pon_curve (every call, any score: the curve ceiling)
and .../pon_calls (PASS calls). The run directories are those the main session writes under
/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/ (read-only here); a run without the files is
reported as pending. Joins them with per_truth.tsv.gz (s5): per platform, the chr1 in-BED truth INDELs that are
  lost_d9 / lost_d46  - a tensor on linear, none on the arm (and the absorbed subset: no candidate, ALT spelled)
and counts how many of them are a baseline TP of the linear run (curve / PASS) and of the arm's own run.
Writes linear_recovery.md (experiment dir). Run: python s6_linear_recovery.py [runs_root]
"""
import csv, gzip, os, sys
from e5b_common import DATA

EXP = os.path.dirname(os.path.abspath(__file__))
ROOT = sys.argv[1] if len(sys.argv) > 1 else '/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013'
RUNS = {  # (sample, platform) -> {arm: run name}
    ('HG008T', 'Illumina'): dict(linear='HG008_Illumina_INDEL_linear_7fc7396', d9='HG008_Illumina_INDEL_7fc7396', d46='HG008_Illumina_INDEL_d46_7fc7396'),
    ('HG008T', 'PacBio'): dict(linear='HG008_PacBio_INDEL_linear_7fc7396', d9='HG008_PacBio_INDEL_7fc7396', d46='HG008_PacBio_INDEL_d46_7fc7396'),
    ('HG008T', 'ONT'): dict(linear='HG008_ONT_INDEL_linear_7fc7396', d9='HG008_ONT_INDEL_7fc7396', d46='HG008_ONT_INDEL_d46_7fc7396'),
    # COLO829T d9: the RIS d9 INDEL model is the purity-joint one (P100 + P75 + P50 + P25), tested on the P100 set
    ('COLO829T', 'Illumina'): dict(linear='COLO829_Illumina_INDEL_linear_7fc7396', d9='COLO829_Illumina_INDEL_joint_7fc7396', d46='COLO829_Illumina_INDEL_d46_7fc7396'),
    ('COLO829T', 'fiberseq'): dict(linear='COLO829_fiberseq_INDEL_linear_7fc7396', d9='COLO829_fiberseq_INDEL_joint_7fc7396', d46='COLO829_fiberseq_INDEL_d46_7fc7396'),
    ('COLO829T', 'ONT'): dict(linear='COLO829_ONT_INDEL_linear_7fc7396', d9='COLO829_ONT_INDEL_joint_7fc7396', d46='COLO829_ONT_INDEL_d46_7fc7396'),
}


def tp_keys(run, which, set_name):
    """Baseline TP alleles (chrom, pos, REF, ALT) of an rtg output directory, or None if missing. The VCF directory is
    vcf_chr1 (one-set runs) or vcf_chr1_<set> (joint runs: the P100 set)."""
    hits = []
    for d in (f'{ROOT}/{run}/vcf_chr1/eval_bed/vcfeval/{which}', f'{ROOT}/{run}/vcf_chr1_{set_name}/eval_bed/vcfeval/{which}'):
        p = f'{d}/tp-baseline.vcf.gz'
        if os.path.exists(p) and os.path.exists(f'{d}/done'):
            hits.append(p)
    if not hits:
        return None
    keys = set()
    with gzip.open(hits[0], 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            x = line.split('\t', 6)
            for a in x[4].split(','):
                keys.add((x[0], int(x[1]), x[3].upper(), a.upper()))
    return keys


rows = [r for r in csv.DictReader(gzip.open(f'{DATA}/per_truth.tsv.gz', 'rt'), delimiter='\t')
        if r['kind'] in ('INS', 'DEL') and r['chrom'] == 'chr1' and r['in_bed'] == 'True' and r['sample'] in ('HG008T', 'COLO829T')]
L = ['# Lost / absorbed truth INDELs called after the PoN (chr1, in the BED; rtg baseline TPs)', '',
     'curve = any score (the ceiling), PASS = at the run\'s recall-0.9 threshold. "pending" = the run\'s rtg output is not there yet.', '',
     '| sample | platform | set | n | not PoN-tagged | linear TP curve / PASS | d9 TP curve / PASS | d46 TP curve / PASS |',
     '|---|---|---|---:|---:|---|---|---|']
pending = []
for (s, p), arms in RUNS.items():
    tp = {}
    for arm, run in arms.items():
        for which in ('pon_curve', 'pon_calls'):
            k = tp_keys(run, which, f'{s}_{p}')
            tp[(arm, which)] = k
            if k is None:
                pending.append(f'{run}/{which}')
    sub = [r for r in rows if r['sample'] == s]
    key = lambda r: (r['chrom'], int(r['vcf_pos']), r['vcf_ref'].upper(), r['vcf_alt'].upper())
    sets = {'all truth INDEL': sub}
    for g in ('d9', 'd46'):
        lost = [r for r in sub if r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] != 'T']
        sets[f'lost on {g}'] = lost
        sets[f'absorbed on {g}'] = [r for r in lost if r[f'{p}_{g}'] == 'N' and r[f'{g}_spelled'] == '1']
    for name, x in sets.items():
        cells = []
        for arm in ('linear', 'd9', 'd46'):
            a, b = tp.get((arm, 'pon_curve')), tp.get((arm, 'pon_calls'))
            if arm not in arms:
                cells.append('no run')
            elif a is None or b is None:
                cells.append('pending')
            else:
                cells.append(f'{sum(key(r) in a for r in x):,} / {sum(key(r) in b for r in x):,}')
        L.append(f'| {s} | {p} | {name} | {len(x):,} | {sum(r["pon"] == "" for r in x):,} | ' + ' | '.join(cells) + ' |')
L += ['', f'pending rtg outputs: {len(pending)}', ''] + [f'- {x}' for x in pending]
open(f'{EXP}/linear_recovery.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
