#!/bin/bash
#SBATCH -J s7_giraffe
#SBATCH -p general
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s7_giraffe-%j.out
# s7_assembly_paths.py: assembly windows (HG008-N v6.2 + HG008-T v3.2) on the d9 graph, vg 1.65 giraffe (short-read index)
set -euo pipefail
D=/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001
G=/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9
VG=/scratch/wzhang/bin/vg_v1.65.0
/usr/bin/time -v $VG giraffe -Z $G.gbz -m $G.shortread.withzip.min -z $G.shortread.zipcodes -d $G.dist -f $D/asm_windows.fq -t 16 > $D/asm_windows.gam 2> $D/asm_windows.giraffe.log
$VG view -a $D/asm_windows.gam > $D/asm_windows.gam.json
tail -n 25 $D/asm_windows.giraffe.log | grep -E 'Elapsed|Maximum resident'
wc -l $D/asm_windows.gam.json
