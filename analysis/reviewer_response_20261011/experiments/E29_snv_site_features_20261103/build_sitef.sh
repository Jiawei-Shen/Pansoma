#!/bin/bash
# SNV site features (code/: repository ca50ab0 + the SNV builder of site_features.py) of one tensor set.
# usage: sbatch build_sitef.sh SET    -> /scratch/jshen/data/snv_sitefeat_20261103/site_features/<SET>.SNV.npz
#SBATCH -p general
#SBATCH -x maotai,guinness,meow,miao
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 4:00:00
set -euo pipefail
SET=$1
module load bcftools/1.18
P=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs
T=/scratch/jshen/data/pansoma_v2_tensors/$SET/tensors
FASTA=$(/wanglab/jshen/anaconda3/bin/python -c "import json; print(json.load(open('$T/SNV/labels.manifest.json'))['fasta'])")
echo "host $(hostname) job ${SLURM_JOB_ID:-} $SET fasta $FASTA $(date)"
cd /scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103/code/machine_learning
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python -m pansoma_net_v2.site_features build --kind SNV \
  --tensors $T --fasta $FASTA --gnomad $P/af-only-gnomad.hg38.vcf.gz \
  --colors $P/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz --dbsnp $P/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --g1000 $P/1000g_pon.hg38.vcf.gz --out /scratch/jshen/data/snv_sitefeat_20261103/site_features \
  --cache-dir /scratch/jshen/data/pansoma_net_v2_runs/index_cache --workers 4 --bcftools $(which bcftools)
