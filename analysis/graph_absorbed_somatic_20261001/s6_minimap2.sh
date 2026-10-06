#!/bin/bash
#SBATCH -J s6_minimap2
#SBATCH -p general
#SBATCH -c 16
#SBATCH --mem=26G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s6_minimap2-%j.out
# s6_assembly_seq.py: REF/ALT context queries vs the HG008-N v6.2 (both haps) and HG008-T v3.2 assemblies
set -euo pipefail
D=/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001
MM=/opt/apps/minimap2/2.28/minimap2
$MM -ax sr --secondary=yes -N 20 --cs -t 16 /scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_bothhaps_polished6.2.fasta.gz $D/asm_queries.fa > $D/asm_normal.sam 2> $D/asm_normal.minimap2.log
$MM -ax sr --secondary=yes -N 20 --cs -t 16 /scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/assembly/HG008T_v3.2.fasta $D/asm_queries.fa > $D/asm_tumor.sam 2> $D/asm_tumor.minimap2.log
tail -n 3 $D/asm_normal.minimap2.log $D/asm_tumor.minimap2.log
