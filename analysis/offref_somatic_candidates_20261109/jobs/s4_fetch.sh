#!/bin/bash
# s4c: HG008-N reads ($1 = HiFi | Illumina) over the HG008T regions, one remote connection; resumable (rerun = continue).
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 16:00:00
set -euo pipefail
export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
echo "host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "s4_fetch %es %MKB" /wanglab/jshen/anaconda3/bin/python -I -B $A/s4_fetch_hg008n.py $1
