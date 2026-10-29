# E14b: ONT purity dilutions without the read-source base-quality cue (R1.5, R2.4, R3.10, R1.7)

In the in-silico dilutions of COLO829T with COLO829BL, tumor- and normal-derived reads come from different runs and
differ in base quality, a cue a real tumor lacks. For Illumina it inflated the dilution scores (R1.5; E26 §5 and the
RIS SNV run of 2026-10-13), but those readnorm runs also changed the normalisation and the SNV batch, so they are not a
clean comparison. This experiment is a clean pair for ONT: two models that differ only in the BQ plane.

## Setup

- `run_joint.sh NAME control|readnorm` (copied from `tmp/ont_readnorm_20261026/`, with `submit.sh`): COLO829T ONT
  P100 + P75 + P50 + P25 (HPRC v1.1 d9 sets of 2026-10-01), SNV, train chr2-22 (5 % node blocks for validation in every
  set), test chr1 of each set. Recipe of the RIS run `COLO829_ONT_SNV_joint_7fc7396_gpu4` (50.4 M network, pow:0.75,
  lr 1e-4, 20 epochs) with its global batch of 4,096 on one H100 (`--batch-size 4096 --micro-batches 4`),
  `--norm train`, seed 0, same code (repo HEAD of 2026-10-26) and hardware (tequila H100) for both.
- `readnorm` trains and predicts through `ris/tools/readnorm.py` (BQ centred on each read's own median); `control` runs
  the package unchanged.
- Prediction with the off-reference rescue on (`--offref-site90 10`), as the RIS predictions; scoring `vcf.sh`
  (rtg vcfeval, chr1, benchmark BED, PoN), PASS = the model's recall-0.9 validation threshold, shared by all purities.
- Jobs (tequila, 2026-10-27): control 452323 and readnorm 452328 (9 h 18 min each), vcf 452324-452327 / 452329-452332.
  SNV chunk caches 452409-452412 (the first try on meow, 452319-452322, was BeeGFS-stalled and cancelled).

## Results (`tables.md`, `collect.sh`)

PASS F1 after population filtering (chr1, BED):

| purity | control | readnorm | difference |
|---|---|---|---|
| 100 % | 0.733 | 0.724 | -0.009 |
| 75 % | 0.721 | 0.693 | -0.028 |
| 50 % | 0.691 | 0.666 | -0.025 |
| 25 % | 0.714 | 0.700 | -0.014 |

1. **No dilution scores above the undiluted sample, with or without the cue** (control 0.691-0.721 against 0.733;
   readnorm 0.666-0.700 against 0.724). Illumina before the fix had P75-P25 at 0.804-0.840 against 0.675 at P100.
2. **Centring costs 0.009 at 100 % purity, where there is no read-source difference, and 0.014-0.028 in the
   dilutions.** The excess over P100 (0.005-0.019) bounds what the cue added to the ONT dilution scores: at most about
   0.02 F1 [one seed; seed-to-seed spread of this recipe not measured here].
3. Before population filtering the same order holds (control 0.578 / 0.464 / 0.402 / 0.467, readnorm 0.537 / 0.403 /
   0.364 / 0.441). Recall at PASS is the same in both (0.639-0.662 after filtering); the difference is in precision.
4. Both runs score above the RIS run of 2026-10-12 (0.677 / 0.654 / 0.639 / 0.666; `analysis/benchmark_7fc7396_20261012/
   colo829t_purity_joint_7fc7396_20261012/README.md`), which used the code of 7fc7396 and four GPUs [not
   investigated]; the comparison here is within the pair only.

## Caveats

- SNVs only: ONT INDEL dilutions were at best F1 0.051-0.068 at every purity, without the Illumina pattern.
- One seed per model; one tumor.
