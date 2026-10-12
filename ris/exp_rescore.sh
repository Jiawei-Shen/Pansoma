#!/bin/bash
# CPU job: the call-level INDEL rescorer (tools/rescore_indel.py) over exp_submit.sh runs.
#   exp_rescore.sh [--chr2-only] [--fit "FIT OPTIONS"] [--feats DIR] NAME RUN:SET [RUN:SET ...]
# --fit passes options to rescore_indel.py fit (e.g. "--feature-set v2 --floors 0.01 0.02 0.05"); --feats reuses the
# feature cache of another rescore run directory (features do not depend on the fit options).
# 1. features of each RUN:SET on chr2 (labelled by the run's chr2 rtg curve) -> NAME/feats/
# 2. fit on chr2 only: 'pooled' (all given sets) and 'loso_<SAMPLE>' (fitted on the other sample's sets) -> model.json
# 3. unless --chr2-only: features on chr1, then each set's chr1 calls rescored with 'pooled' and with its loso
#    variant, PoN-tagged at the floor chosen on chr2 (tools/pon_tag.py), rtg vcfeval in the BED and without it.
# Output: $PANSOMA_RUNS/NAME/
set -euo pipefail
shopt -s inherit_errexit  # a failing command inside $(...) stops the job
source "$(dirname "$0")/env.sh"
CHR2_ONLY=0 FIT="" FEATS=""
while [ $# -gt 0 ]; do
  case $1 in --chr2-only) CHR2_ONLY=1; shift;; --fit) FIT=$2; shift 2;; --feats) FEATS=$2; shift 2;; *) break;; esac
done
NAME=$1; shift
D=$PANSOMA_RUNS/$NAME; mkdir -p "$D"
if [ -n "$FEATS" ]; then ln -sfn "$FEATS" "$D/feats"; else mkdir -p "$D/feats"; fi
SDF=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf
P=$PANSOMA_DATA/panel_of_normal_VCFs
PONS=("$P/af-only-gnomad.hg38.vcf.gz" "$P/Homo_sapiens_assembly38.dbsnp138.vcf.gz" "$P/1000g_pon.hg38.vcf.gz"
      "$P/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
cd "$REPO/machine_learning"
export PYTHONPATH="$REPO/machine_learning:$REPO/scripts:$REPO/ris/tools${PYTHONPATH:+:$PYTHONPATH}"
RS="$REPO/ris/tools/rescore_indel.py"
chroms=(chr2); [ $CHR2_ONLY = 1 ] || chroms=(chr2 chr1)
npz2=()
for rs in "$@"; do
  run=${rs%%:*} set=${rs#*:}
  for c in "${chroms[@]}"; do
    f=$D/feats/${run}__${set}__$c.npz
    [ -s "$f" ] || python -u "$RS" features "$PANSOMA_RUNS/$run" "$set" "$c" "$f" &
  done
  npz2+=("$D/feats/${run}__${set}__chr2.npz")
  while [ "$(jobs -r | wc -l)" -ge "${LSB_DJOB_NUMPROC:-4}" ]; do sleep 2; done
done
wait
python -u "$RS" fit "${npz2[@]}" --out "$D/model.json" $FIT | tee "$D/fit.txt"
[ $CHR2_ONLY = 1 ] && exit 0
for rs in "$@"; do
  run=${rs%%:*} set=${rs#*:}
  read -r TRUTH BED < <(python -c "import json, sys; t = json.load(open(sys.argv[1]))['truth']['somatic']; print(t['vcf'], t['bed'])" \
      "$PANSOMA_DATA/pansoma_v2_tensors/$set/tensors/INDEL/labels.manifest.json")
  for v in pooled "loso_${set%%_*}"; do
    floor=$(python -c "import json, sys; print(json.load(open(sys.argv[1]))['variants'][sys.argv[2]]['floor'])" "$D/model.json" "$v")
    o=$D/vcf_${run}__${set}_chr1_$v; rm -rf "$o"; mkdir -p "$o"
    lin=$(python -c "import numpy as np, sys; print(str(np.load(sys.argv[1])['vcf']))" "$D/feats/${run}__${set}__chr1.npz")
    out=$o/$(basename "$lin" .linear.vcf.gz).rescored.vcf.gz
    python -u "$RS" apply --model "$D/model.json" --variant "$v" --features "$D/feats/${run}__${set}__chr1.npz" --vcf "$lin" --out "$out"
    tagged=${out%.vcf.gz}.floor$floor.pon.vcf.gz
    python -u "$REPO/ris/tools/pon_tag.py" "$floor" "$out" "$tagged" "${PONS[@]}"
    python -u -m pansoma_net_v2.vcfeval --calls "$out" --pon-vcf "$tagged" --truth "$TRUTH" --bed "$BED" --sdf "$SDF" --rtg "$RTG" --output "$o/eval_bed"
    python -u -m pansoma_net_v2.vcfeval --calls "$out" --pon-vcf "$tagged" --truth "$TRUTH" --sdf "$SDF" --rtg "$RTG" --output "$o/eval_nobed"
  done
done
for x in "$D"/vcf_*; do echo "-- $(basename "$x")"; grep -E "^(raw|pon) " "$x/eval_bed/report.txt"; done
