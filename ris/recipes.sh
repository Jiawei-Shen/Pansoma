# Sourced by snv_pipeline.sh and joint_pipeline.sh: RECIPE = the train options of the recipe of $KIND in
# machine_learning/pansoma_net_v2/README.md ("SNV recipe", "INDEL recipe"; 7fc7396: INDEL lr 1e-4, was 4e-4 in 45d56ad;
# 9919d4a: the INDEL pipeline, i.e. labels that follow the benchmark, site features, best epoch by the post-PoN
# validation F1 and PASS at the post-PoN validation recall 0.8 of the callable positives, "ppr8").
# SNV, by platform (PLATFORM, else from the set names SET / SETS: a name with "_ONT" is ONT, any other is PacBio HiFi /
# Fiber-seq / Illumina; a mix of the two is an error):
#   PacBio HiFi, Fiber-seq, Illumina: sets built at SNV AF 0.05, the SNV site features, non tensors with AF <= 0.1
#     sampled 25% per epoch at loss weight 1, best epoch by the post-PoN validation F1, PASS at the post-PoN validation
#     recall 0.95 of the callable positives (E30:
#     analysis/reviewer_response_20261011/experiments/E30_snv_af005_20261103);
#   ONT: the original SNV recipe (sets at SNV AF 0.08, no site features, every negative, PASS at validation recall 0.9).
# The site-feature recipes read the files of every training and predicted set from SITE_FEATURES (site_features.sh,
# KIND=SNV for SNV; env.sh's default).
SF_DIR="${SITE_FEATURES:-${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}/site_features}"
case ${KIND:-SNV} in
  SNV)   # the 50.4 M model
    names=",${SET:-}${SETS:+,$SETS},"
    platform=${PLATFORM:-}
    if [ -z "$platform" ]; then
      case $names in
        *_ONT*) platform=ONT; [[ $(tr ',' '\n' <<< "$names" | grep -v '^$' | grep -cv _ONT) -gt 0 ]] &&
                  { echo "recipes.sh: ONT and other sets in one SNV model ($names); set PLATFORM" >&2; exit 2; } ;;
        *) platform=PacBio_Illumina;;
      esac
    fi
    SNV_BASE=(--depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75
              --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth)
    case $platform in
      ONT) RECIPE=("${SNV_BASE[@]}" --non-fraction 1.0);;
      PacBio|Illumina|PacBio_Illumina|fiberseq)
           RECIPE=("${SNV_BASE[@]}" --non-fraction 0.25 --keep-non-af 0.100001 --no-importance-weights
                   --site-features "$SF_DIR" --best-by postpon --call-recall 0.95);;
      *) echo "recipes.sh: PLATFORM must be PacBio, Illumina, fiberseq or ONT, not $platform" >&2; exit 2;;
    esac;;
  INDEL) # the 8.7 M model, every candidate kept
    RECIPE=(--depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --batch-size 1024 --lr 1e-4 --weight-decay 0.05
            --warmup-epochs 1 --class-weights sqrt --non-fraction 1.0 --block v2
            --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
            --non-reasons repeat_indel_without_truth --non-reasons-in-bed outside_confident_region
            --site-features "$SF_DIR" --best-by postpon);;
  *) echo "KIND must be SNV or INDEL, not ${KIND}" >&2; exit 2;;
esac
