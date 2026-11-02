#!/bin/bash
#SBATCH -p general
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 6:00:00
set -euo pipefail
S=COLO829T
E=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_ep8/code
D=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_ep8/DATA/asm/$S
module load minimap2/2.28
echo "sample $S host $(hostname) job ${SLURM_JOB_ID:-}"
/wanglab/jshen/anaconda3/bin/python $E/e9d_asm.py queries $S
/wanglab/jshen/anaconda3/bin/python -c "import sys; sys.path.insert(0, '$E'); import e9c_lib as L; [print(t, f) for t, r, f in L.ASSEMBLIES['$S']]" | while read T F; do
  /usr/bin/time -f "$T minimap2 %es %MKB" minimap2 -cx sr --cs --secondary=yes -N 50 -t ${SLURM_CPUS_PER_TASK:-8} "$F" $D/queries.fa > $D/$T.paf.tmp 2> $D/$T.minimap2.log
  mv $D/$T.paf.tmp $D/$T.paf
  tail -2 $D/$T.minimap2.log
done
/wanglab/jshen/anaconda3/bin/python $E/e9d_asm.py classify $S
