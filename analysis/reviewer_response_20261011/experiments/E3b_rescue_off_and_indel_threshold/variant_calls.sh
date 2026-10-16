#!/bin/bash
# E3b: two re-calls of an existing chr1 prediction set (no new inference), each through the pipeline's own
# jobs/vcf.sh (graph_vcf -> linear_vcf -> PoN -> rtg vcfeval, in the BED and without it):
#   norescue  SNV runs: the predictions with predict's off-reference rescue fields removed (offref_call,
#             p_offref_*), so graph_vcf calls every tensor by the checkpoint's validation recall-0.9 threshold
#             (R3.4: complete output without the rescue). Output: <run>__norescue/vcf_chr1
#   valbest   INDEL runs: PASS threshold = the checkpoint's validation-chosen best truth-level F1 threshold
#             (metrics.jsonl val.truth.best.threshold at the best epoch; chosen on chr2-22 validation blocks, never on
#             chr1) instead of the recall-0.9 threshold (R1.6 / R3.10 / R3.m1: a held-out INDEL operating point).
#             Output: <run>/vcf_chr1_valbest
# usage: variant_calls.sh MODE RUN     (RUN under pansoma_net_v2_runs/reviewer_cpu_20261013; read-only inputs)
set -euo pipefail
MODE=$1; RUN=$2
R=/scratch/jshen/data/pansoma_net_v2_runs; S=$R/reviewer_cpu_20261013; P=/wanglab/jshen/anaconda3/bin/python
VCF=/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2/jobs/vcf.sh
SRC=$(readlink -f $S/$RUN/best.pth); SRCRUN=$(basename $(dirname $SRC))
echo "$MODE $RUN (checkpoint $SRC) host $(hostname) job ${SLURM_JOB_ID:-}"
case $MODE in
  norescue)
    O=$S/${RUN}__norescue
    [ -s $O/vcf_chr1/eval_nobed/report.json ] && { echo "done already"; exit 0; }
    mkdir -p $O/test_chr1
    ln -sfn $SRC $O/best.pth
    for f in $S/$RUN/test_chr1/*.SNV.predictions.ndjson.gz; do
      b=$(basename $f)
      $P -c "
import gzip, json, sys
drop = ('offref_call', 'p_offref_somatic', 'p_offref_germline', 'p_offref_non')
n = m = 0
with gzip.open('$f', 'rt') as i, gzip.open('$O/test_chr1/$b.tmp', 'wt') as o:
    for line in i:
        r = json.loads(line); n += 1
        if any(k in r for k in drop):
            m += 1
            for k in drop: r.pop(k, None)
        o.write(json.dumps(r, separators=(',', ':')) + '\n')
print(f'{n} records, rescue fields removed from {m}')
"
      mv $O/test_chr1/$b.tmp $O/test_chr1/$b
      cp $S/$RUN/test_chr1/${b%.predictions.ndjson.gz}.metrics.json $O/test_chr1/
    done
    bash $VCF reviewer_cpu_20261013/${RUN}__norescue SNV ;;
  valbest)
    O=$S/$RUN/vcf_chr1_valbest
    [ -s $O/eval_nobed/report.json ] && { echo "done already"; exit 0; }
    T=$($P -c "
import json, torch
ck = torch.load('$SRC', map_location='cpu', weights_only=False, mmap=True)
ep = ck['epoch']
rows = [json.loads(l) for l in open('$R/ris_7fc7396_20261012/$SRCRUN/metrics.jsonl')]
r = [x for x in rows if x.get('epoch') == ep][0]
print(r['val']['truth']['best']['threshold'])
")
    echo "best epoch validation truth-best threshold: $T"
    echo "$T" > $S/$RUN/valbest_threshold.txt
    THRESHOLD=$T bash $VCF reviewer_cpu_20261013/$RUN INDEL $O 2>&1 | sed 's/^/  /' ;;
esac
