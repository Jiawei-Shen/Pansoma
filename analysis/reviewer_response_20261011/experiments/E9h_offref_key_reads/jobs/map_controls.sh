#!/bin/bash
# E9h control queries (build_controls.py): minimap2 2.28 -ax sr with the settings of map.sh, one assembly haplotype per
# array task: 0 HG008-N v6.2 hap1, 1 HG008-N v6.2 hap2, 2 HG008-T v3.2 hap1, 3 HG008-T v3.2 hap2.
#SBATCH -p general
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 1:00:00
set -euo pipefail
G=/scratch/jshen/data/HG008_GIAB
D=/scratch/jshen/data/reviewer_experiments_20261013/E9h_offref_key_reads
NAMES=(N_hap1 N_hap2 T_hap1 T_hap2)
FASTA=($G/HG008N_curatedv6_250714_polished6.2_hap1.fasta $G/HG008N_curatedv6_250714_polished6.2_hap2.fasta.gz \
       $G/HG008T_v3.2.hap1.fasta.gz $G/HG008T_v3.2.hap2.fasta.gz)
i=${SLURM_ARRAY_TASK_ID}
T=${NAMES[$i]}
F=${FASTA[$i]}
module load minimap2/2.28
echo "task $i $T $F host $(hostname) job ${SLURM_JOB_ID:-} start $(date)"
/usr/bin/time -f "$T controls minimap2 %es %MKB" minimap2 -ax sr --secondary=yes -N 20 --cs -t ${SLURM_CPUS_PER_TASK:-8} \
  "$F" $D/controls.fa > $D/$T.controls.sam.tmp 2> $D/$T.controls.minimap2.log
mv $D/$T.controls.sam.tmp $D/$T.controls.sam
tail -3 $D/$T.controls.minimap2.log
echo "end $(date)"
