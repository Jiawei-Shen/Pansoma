#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J q6_side_colo_ill_d9
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/side/logs/%x_%j.out
# COLO829T_Illumina (d9): the local index cache (2026-10-01) predates the 2026-10-02 rebuild of this set, so no cache
# check here; offnode.side_arrays checks the node order against the index predict builds.
/wanglab/jshen/anaconda3/envs/sambcfenv/bin/python -u /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/side/build_side.py \
  /scratch/jshen/data/pansoma_v2_tensors/COLO829T_Illumina/tensors/SNV /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/side/COLO829T_Illumina.SNV.side.npz
