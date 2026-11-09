#!/bin/bash
# E11b collect: score every finished (arm, platform) set and rebuild the tables.
#   [E11B_SAMPLE=COLO829BL] bash collect.sh   (login node is fine: ~2-5 min per Illumina set, < 8 GB)
set -euo pipefail
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E11b_simulated_alignment
PY=/wanglab/jshen/anaconda3/bin/python
export E11B_SAMPLE=${E11B_SAMPLE:-HG008N}
# sample settings from e11b_common.SAMPLES: data dir, CHM13 dipcall prefix, contig rule, results dir
read -r D CHM13_PREFIX RULE RES < <($PY -c "import sys; sys.path.insert(0, '$E'); import e11b_common as C; \
print(C.BASE, C.CFG['chm13'], C.CFG['contig_rule'], C.CFG['results'])")
RES=$E/$RES
echo "[INFO] sample $E11B_SAMPLE data $D results $RES"
mkdir -p "$D/eval" "$D/chm13_truth"
# CHM13 control truth (needs both haplotype -> CHM13 PAFs; HG008-N: chm13_paf.sbatch, COLO829BL: qfu's dipcall_chm13)
if [[ -s $CHM13_PREFIX.hap1.bed && -s $CHM13_PREFIX.hap2.bed ]]; then
  [[ -s $D/chm13_truth/segments.tsv ]] || \
    $PY "$E/build_segments.py" --dipcall-prefix "$CHM13_PREFIX" --contig-rule "$RULE" --out "$D/chm13_truth"
  for P in illumina hifi ont; do
    if [[ -s $D/sim/$P/truth.npz && ! -s $D/chm13_truth/$P.chm13_truth.npz ]]; then
      $PY "$E/chm13_truth.py" --platform "$P"
    fi
  done
fi
for P in illumina hifi ont; do
  for A in linref linear d9 d46 chm13; do
    [[ $A != chm13 || -s $D/chm13_truth/$P.chm13_truth.npz ]] || continue
    if [[ -s $D/aln/$A/$P.bam && ( $A == linref || -s $D/aln/$A/$P.gaf.gz ) && ! -s $D/eval/$A.$P.json ]]; then
      $PY "$E/evaluate.py" --arm "$A" --platform "$P" --out "$D/eval"
    fi
  done
done
mkdir -p "$RES"
# truth check (first 20,000 Illumina rows) and ONT error profile of the simulated reads (2,000 reads vs their templates)
if [[ -s $D/sim/illumina/truth.npz && ! -s $RES/truth_check.txt ]]; then
  $PY "$E/check_truth.py" "$D/sim/illumina" /scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta \
    > "$RES/truth_check.txt"
fi
if [[ -s $D/sim/ont/sim.fq.gz && ! -s $RES/ont_error_profile.txt ]]; then
  { gzip -dc "$D/sim/ont/sim.fq.gz" || true; } | head -n 8000 > "$D/sim/ont/first2000.fq"  # gzip gets SIGPIPE under pipefail
  /opt/apps/minimap2/2.28/minimap2 -c -t 2 "$D/sim/ont/pbsim/templates.fa" "$D/sim/ont/first2000.fq" 2>/dev/null \
    | awk '$1==$6 && $0 ~ /tp:A:P/ {for(i=13;i<=NF;i++) if($i ~ /^NM:i:/){split($i,x,":"); nm[$1]+=x[3]; bl[$1]+=$11}}
           END{for(r in nm) print nm[r]/bl[r]}' | sort -g \
    | awk '{v[NR]=$1} END{printf "ONT-like reads vs their own templates (minimap2 -c, primary blocks): n=%d, NM/block median %.4f, q25 %.4f, q75 %.4f, q90 %.4f\n", NR, v[int(NR/2)], v[int(NR/4)], v[int(3*NR/4)], v[int(0.9*NR)]}' \
    > "$RES/ont_error_profile.txt"
fi
[[ -s $RES/segment_stats.json ]] || \
  $PY "$E/segment_stats.py" "$D/sim" /scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta "$RES"
$PY "$E/summarize.py" --evals "$D/eval" --out "$RES"
