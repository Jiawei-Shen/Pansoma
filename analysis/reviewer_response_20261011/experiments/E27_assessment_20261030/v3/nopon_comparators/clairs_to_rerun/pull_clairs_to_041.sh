#!/bin/bash
#SBATCH -p general
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -J pull_clairs041
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/clairs_to_nopon_rerun_20261102/logs/pull_%j.out
# ClairS-TO v0.4.1 (the version of the user's COLO829T runs) as a SIF; build temp on node-local /tmp.
set -euo pipefail
module load singularity/4.1.0
SIF=/scratch/jshen/bin/clairs-to_v0.4.1.sif
[ -s $SIF ] && { echo "$SIF exists"; exit 0; }
T=/tmp/jshen_singularity_$SLURM_JOB_ID; mkdir -p $T/cache $T/tmp
export SINGULARITY_CACHEDIR=$T/cache SINGULARITY_TMPDIR=$T/tmp
(cd $T && singularity pull clairs-to_v0.4.1.sif docker://hkubal/clairs-to:v0.4.1)
cp $T/clairs-to_v0.4.1.sif $SIF.partial && mv $SIF.partial $SIF
echo "node-local temp left for a literal-path cleanup: $(hostname):$T"
singularity exec $SIF /opt/bin/run_clairs_to --version || true
