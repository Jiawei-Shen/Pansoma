#!/bin/bash
# s2: minimap2 -x sr of one tumor's unique queries to one target (assembly haplotype or GRCh38; PAF with cs;
# secondaries, -N 50). Only alignments covering >= 90 % of the query with NM <= 60 are kept (NM counts the site's own
# INDEL bases, so the cap is loose; s3 ranks normal alignments by their differences outside the event window).
# $1 = TUMOR; task = SLURM_ARRAY_TASK_ID (0-4) = the tumor's line in targets.tsv.
#SBATCH -p general
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 3:00:00
set -euo pipefail
J=/scratch/jshen/Github/Pansoma/tmp/offref_somatic_candidates_20261109
D=/scratch/jshen/data/pansoma_net_v2_runs/offref_somatic_candidates_20261109/s2/$1
read T NAME ROLE FA < <(awk -v t=$1 '$1 == t' $J/targets.tsv | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
echo "tumor $T target $NAME ($ROLE) $FA host $(hostname) job ${SLURM_JOB_ID:-}"
module load minimap2/2.28
/usr/bin/time -f "minimap2 %es %MKB" minimap2 -cx sr --cs --secondary=yes -N 50 -t ${SLURM_CPUS_PER_TASK:-8} "$FA" $D/queries.fa 2> $D/$NAME.minimap2.log \
  | awk -F'\t' '{nm = 99; for (i = 13; i <= NF; i++) if (substr($i, 1, 5) == "NM:i:") nm = substr($i, 6) + 0;
                 if ($4 - $3 >= 0.9 * $2 && nm <= 60) print}' > $D/$NAME.paf.tmp
mv $D/$NAME.paf.tmp $D/$NAME.paf
tail -3 $D/$NAME.minimap2.log
wc -l $D/$NAME.paf
