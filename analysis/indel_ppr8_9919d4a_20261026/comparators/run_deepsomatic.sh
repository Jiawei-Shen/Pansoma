#!/bin/bash
# DeepSomatic 1.9.0 tumor-only on chr1 of one BAM, in google/deepsomatic:1.9.0-gpu (call_variants on the GPU), with the
# settings of the user's P100 runs (outputs/deepsomatic_COLO829T_*): the platform's *_TUMOR_ONLY model,
# GRCh38.primary.fasta, --regions chr1, the default PoN filtering; num_shards = the job's slots.
#   usage: run_deepsomatic.sh NAME BAM MODEL_TYPE OUT   (submit_deepsomatic.sh)
set -euo pipefail
NAME=$1 BAM=$2 MODEL=$3 OUT=$4
REF=/storage2/fs1/epigenome/Active/jiawei.s/Dataset/REF/GRCh38.primary.fasta
RD=/opt/deepvariant/bin/deepsomatic/run_deepsomatic
N=${LSB_DJOB_NUMPROC:-16}
mkdir -p "$OUT/logs" "$OUT/intermediate_results_dir"
echo "host $(hostname) job ${LSB_JOBID:-}; $(date); bam $BAM; model $MODEL; shards $N"
nvidia-smi -L || true
HELP=$($RD --helpfull 2>&1 || true)
grep -q -- "$MODEL" <<< "$HELP" || { echo "model type $MODEL is not in run_deepsomatic --helpfull"; exit 2; }
# postprocess_variants starts one worker per visible core, each loading the PoN VCF (E8: killed at 98 GB): 4 workers
extra=(); grep -q postprocess_variants_extra_args <<< "$HELP" && extra=(--postprocess_variants_extra_args=cpus=4)
CMD=($RD --model_type="$MODEL" --ref="$REF" --reads_tumor="$BAM" --sample_name_tumor="$NAME"
     --output_vcf="$OUT/$NAME.chr1.vcf.gz" --output_gvcf="$OUT/$NAME.chr1.g.vcf.gz" --num_shards="$N"
     --logging_dir="$OUT/logs" --intermediate_results_dir="$OUT/intermediate_results_dir" --regions=chr1
     --use_default_pon_filtering=true "${extra[@]}")
echo "${CMD[*]}" | tee "$OUT/command.txt"
time "${CMD[@]}"
ls -la "$OUT"
date
