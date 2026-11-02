#!/bin/bash
#SBATCH -J assess_m2_fulltruth
#SBATCH -p general
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators/logs/hapmap_m2_fulltruth-%j.out
# Mutect2 (E23 chr1 rerun, tumor-only, no germline resource / PoN) scored on the FULL chr1 HapMap benchmark BED, the
# protocol of Pansoma's own E23 evaluation (regions.bed = chr1, -e SMHTHAPMAP6 v1.6 autosome BED), so graph-carried
# truth counts as FN; then call_vaf_by_value.py with the atomized graph alleles (both columns: not graph-carried / all).
set -euo pipefail
D=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators/hapmap_m2_full; mkdir -p $D
O=/scratch/jshen/data/reviewer_experiments_20261013/E23_mutect2; R=/scratch/jshen/data/pansoma_net_v2_runs
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E23_hapmap_lowvaf
G=/scratch/jshen/data/reviewer_experiments_20261013/E23_graph_alleles.mgi_non_HG005.d9.atomized.tsv.gz
BED=/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/SMHTHAPMAP6_GRCh38_v1.6_somatic_benchmark_autosomes.bed
RTG=/scratch/jshen/bin/rtg-tools-3.12.1/rtg; SDF=/scratch/jshen/REF/GRCh38.primary.sdf; P=/wanglab/jshen/anaconda3/bin/python
if ! command -v java > /dev/null; then module load java; fi
declare -A TRUTH=([SNV]=$R/HapMap_Illumina_SNV_full/vcf_chr1/eval_bed/truth.SNV.vcf.gz [INDEL]=$R/HapMap_Illumina_INDEL/vcf_chr1/eval_bed/truth.INDEL.vcf.gz)
declare -A REG=([SNV]=$R/HapMap_Illumina_SNV_full/vcf_chr1/eval_bed/regions.bed [INDEL]=$R/HapMap_Illumina_INDEL/vcf_chr1/eval_bed/regions.bed)
for m in pass pass_cluster; do for k in SNV INDEL; do
  out=$D/eval.$m.$k
  [ -d $out ] || $RTG RTG_MEM=8g vcfeval -b ${TRUTH[$k]} -c $O/sets/m2.$m.$k.vcf.gz -t $SDF -o $out --squash-ploidy \
    --sample ALT,ALT --bed-regions ${REG[$k]} -e $BED --vcf-score-field INFO.TLOD -T 4
  (cd $E && PYTHONDONTWRITEBYTECODE=1 $P call_vaf_by_value.py --eval-dir $out --kind $k --graph-alleles $G \
    --title "Mutect2 ($m), full chr1 BED (graph-carried truth counted)" --output $D/call_vaf_by_value.Mutect2_$m.$k.md)
done; done
echo done
