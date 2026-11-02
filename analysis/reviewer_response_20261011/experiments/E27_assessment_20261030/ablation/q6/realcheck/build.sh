#!/bin/bash
#SBATCH -p general
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH -J q6_realcheck
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/realcheck/%x_%j.out
# args: SET GAM GRAPH_INDEX HAPLOTYPES CHR_INDEX SNV_MIN_AF INDEL_MIN_AF
set -euo pipefail
S=$1; GAM=$2; GI=$3; H=$4; CHR=$5; SAF=$6; IAF=$7
Q=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/realcheck
P=/wanglab/jshen/anaconda3/bin/python
cd /scratch/jshen/Github/Pansoma
echo "host $(hostname) $(date)"
[ -s $Q/$S.nodes.txt ] || $P $Q/pick_nodes.py $S 400 $Q/$S.nodes.txt
O=$Q/$S.build
/usr/bin/time -f "build %es %MKB" $P -m indexed_gam_pipeline_v4.run build --gam $GAM --output $O/shared --nodes $Q/$S.nodes.txt \
  --graph-index $GI --haplotypes $H --snv-min-af $SAF --indel-min-af $IAF --snv-output $O/SNV --indel-output $O/INDEL \
  --debug-rows --rows 200 --width 101 --min-af 0.05 --max-node-span 10000 --min-mapq 10 --min-variants 3 --min-allele-bq 10 --max-indel-len 50 \
  --max-node-reads 800 --chr-index $CHR --early-af-filter
/wanglab/jshen/anaconda3/envs/sambcfenv/bin/python $Q/check_mask.py $O/SNV $GI
/wanglab/jshen/anaconda3/envs/sambcfenv/bin/python $Q/check_mask.py $O/INDEL $GI
echo "end $(date)"
