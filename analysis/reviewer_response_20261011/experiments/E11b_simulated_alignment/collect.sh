#!/bin/bash
# E11b collect: score every finished (arm, platform) set and rebuild the tables.
#   bash collect.sh            (login node is fine: ~2-5 min per Illumina set, < 8 GB)
set -euo pipefail
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E11b_simulated_alignment
D=/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment
PY=/wanglab/jshen/anaconda3/bin/python
mkdir -p "$D/eval"
# CHM13 control truth (needs both haplotype -> CHM13 PAFs from chm13_paf.sbatch)
if [[ -s $D/chm13_truth/HG008N_CHM13_dipcall.hap1.bed && -s $D/chm13_truth/HG008N_CHM13_dipcall.hap2.bed ]]; then
  [[ -s $D/chm13_truth/segments.tsv ]] || \
    $PY "$E/build_segments.py" --dipcall-prefix "$D/chm13_truth/HG008N_CHM13_dipcall" --out "$D/chm13_truth"
  for P in illumina hifi ont; do
    if [[ -s $D/sim/$P/truth.npz && ! -s $D/chm13_truth/$P.chm13_truth.npz ]]; then
      $PY "$E/chm13_truth.py" --platform "$P"
    fi
  done
fi
for P in illumina hifi ont; do
  for A in linear d9 d46 chm13; do
    [[ $A != chm13 || -s $D/chm13_truth/$P.chm13_truth.npz ]] || continue
    if [[ -s $D/aln/$A/$P.bam && -s $D/aln/$A/$P.gaf.gz && ! -s $D/eval/$A.$P.json ]]; then
      $PY "$E/evaluate.py" --arm "$A" --platform "$P" --out "$D/eval"
    fi
  done
done
mkdir -p "$E/results"
# truth check (first 20,000 Illumina rows) and ONT error profile of the simulated reads (2,000 reads vs their templates)
if [[ -s $D/sim/illumina/truth.npz && ! -s $E/results/truth_check.txt ]]; then
  $PY "$E/check_truth.py" "$D/sim/illumina" /scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta \
    > "$E/results/truth_check.txt"
fi
if [[ -s $D/sim/ont/sim.fq.gz && ! -s $E/results/ont_error_profile.txt ]]; then
  gzip -dc "$D/sim/ont/sim.fq.gz" | head -n 8000 > "$D/sim/ont/first2000.fq"
  /opt/apps/minimap2/2.28/minimap2 -c -t 2 "$D/sim/ont/pbsim/templates.fa" "$D/sim/ont/first2000.fq" 2>/dev/null \
    | awk '$1==$6 && $0 ~ /tp:A:P/ {for(i=13;i<=NF;i++) if($i ~ /^NM:i:/){split($i,x,":"); nm[$1]+=x[3]; bl[$1]+=$11}}
           END{for(r in nm) print nm[r]/bl[r]}' | sort -g \
    | awk '{v[NR]=$1} END{printf "ONT-like reads vs their own templates (minimap2 -c, primary blocks): n=%d, NM/block median %.4f, q25 %.4f, q75 %.4f, q90 %.4f\n", NR, v[int(NR/2)], v[int(NR/4)], v[int(3*NR/4)], v[int(0.9*NR)]}' \
    > "$E/results/ont_error_profile.txt"
fi
[[ -s $E/results/segment_stats.json ]] || \
  $PY "$E/segment_stats.py" "$D/sim" /scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta "$E/results"
$PY "$E/summarize.py" --evals "$D/eval" --out "$E/results"
