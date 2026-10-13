# Sourced by snv_pipeline.sh and joint_pipeline.sh: RECIPE = the train options of the recipe of $KIND in
# machine_learning/pansoma_net_v2/README.md ("SNV recipe", "INDEL recipe"; 7fc7396: INDEL lr 1e-4, was 4e-4 in 45d56ad).
case ${KIND:-SNV} in
  SNV)   # the 50.4 M model
    RECIPE=(--depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75
            --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth);;
  INDEL) # the 8.7 M model, every candidate kept
    RECIPE=(--depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --batch-size 1024 --lr 1e-4 --class-weights sqrt
            --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth);;
  *) echo "KIND must be SNV or INDEL, not ${KIND}" >&2; exit 2;;
esac
