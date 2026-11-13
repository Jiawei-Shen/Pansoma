# Sourced by snv_pipeline.sh and joint_pipeline.sh: RECIPE = the train options of the recipe of $KIND in
# machine_learning/pansoma_net_v2/README.md ("SNV recipe", "INDEL recipe"; 7fc7396: INDEL lr 1e-4, was 4e-4 in 45d56ad;
# 9919d4a: the INDEL pipeline, i.e. labels that follow the benchmark, site features, best epoch by the post-PoN
# validation F1 and PASS at the post-PoN validation recall 0.8 of the callable positives, "ppr8"). INDEL on ONT (user
# 2026-11-13): PASS at the post-PoN validation best-F1 threshold instead ("pp", train --call-rule postpon), the rule
# the final batch pre-registered; ppr8's lower threshold let in hundreds of ONT HP-error calls (HCC1395T chr1 six-tumor
# model 0.304 at ppr8 vs 0.564 at pp, three-seed means). The platform is PLATFORM (Illumina, PacBio, ONT) or, when unset,
# ONT if a set name (SET of snv_pipeline.sh, SETS of joint_pipeline.sh) contains _ONT.
# SNV, every platform (PacBio HiFi, Fiber-seq, Illumina, ONT): the SNV site features, non tensors with AF <= 0.1
#   sampled 25% per epoch at loss weight 1, best epoch by the post-PoN validation F1, PASS at the post-PoN validation
#   recall 0.95 of the callable positives (E30: analysis/reviewer_response_20261011/experiments/E30_snv_af005_20261103).
#   Only the build's SNV AF floor differs: 0.05 for PacBio HiFi / Fiber-seq / Illumina, 0.08 for ONT.
# The site-feature recipes read the files of every training and predicted set from SITE_FEATURES (site_features.sh,
# KIND=SNV for SNV; env.sh's default).
SF_DIR="${SITE_FEATURES:-${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}/site_features}"
case ${KIND:-SNV} in
  SNV)   # the 50.4 M model
    RECIPE=(--depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75
            --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
            --non-fraction 0.25 --keep-non-af 0.100001 --no-importance-weights
            --site-features "$SF_DIR" --best-by postpon --call-recall 0.95);;
  INDEL) # the 8.7 M model, every candidate kept
    RECIPE=(--depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --batch-size 1024 --lr 1e-4 --weight-decay 0.05
            --warmup-epochs 1 --class-weights sqrt --non-fraction 1.0 --block v2
            --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
            --non-reasons repeat_indel_without_truth --non-reasons-in-bed outside_confident_region
            --site-features "$SF_DIR" --best-by postpon)
         plat=${PLATFORM:-}
         if [ -z "$plat" ]; then case ",${SET:-}${SETS:-}," in *_ONT[,_]*) plat=ONT;; esac; fi
         if [ "$plat" = ONT ]; then RECIPE+=(--call-rule postpon); fi;;
  *) echo "KIND must be SNV or INDEL, not ${KIND}" >&2; exit 2;;
esac
