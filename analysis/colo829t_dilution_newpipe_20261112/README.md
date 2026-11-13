# COLO829T purity dilutions with the current SNV pipeline (2026-11-12)

The user (2026-11-12): rerun the COLO829T dilution data with the latest pipeline, on the original dilution sets (not
AF 0.05 rebuilds); ONT with the same recipe as the other platforms (f0171af, pulled during the runs), and no ONT
negative with AF < 0.08 in training. `machine_learning/pansoma_net_v2` Python code last changed in 1497adb.

## What runs (`submit.sh`; LSF job ids in `$PANSOMA_RUNS/logs/submit_dilution_newpipe_20261112.txt`)

Sets: `COLO829T_<P>`, `COLO829T_<P>_P75`, `_P50`, `_P25` for Illumina, fiberseq and ONT (HPRC v1.1 d9, SNV build AF
0.08, labels `5430533aa100`), the sets of the earlier dilution runs. The design of those runs
(`../benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012`) is kept: one model per platform on the four
purities, chr1 of every set left out, 5 % of the chr2-22 node blocks validate, chr1 of each purity tested with one
threshold per model; 4 H100s per model under DDP (`--batch-size 1024` per GPU, global batch 4,096), seed 0.

| step | what | jobs |
|---|---|---|
| S | SNV site features of the 8 Illumina / fiberseq sets (`KIND=SNV ris/site_features.sh`) into `$PANSOMA_RUNS/site_features/<SET>.SNV.npz` | 771046-771053 |
| A | `COLO829_Illumina_SNV_joint_newpipe` (4 H100, 16 slots, 320 GB) | 771059, vcf 771060, 771062-771064 |
| A | `COLO829_fiberseq_SNV_joint_newpipe` (4 H100, 16, 256 GB) | 771065, vcf 771066-771069 |
| B | SNV site features of `COLO829T_ONT_P75` / `_P50` / `_P25` (`COLO829T_ONT.SNV.npz` existed) | 772331-772333 |
| B | `COLO829_ONT_SNV_joint_newpipe` (4 H100, 16, 384 GB; after the three site-feature jobs) | 772334, vcf 772335-772338 |
| C | best of epochs 1-5 (the user: test the best epoch once 5 epochs are done): `snapshot_best5.sh` copies each run's `best.pth` after epoch 5 into `COLO829_<P>_SNV_joint_newpipe_e5/` (fiberseq by hand at 17:34, epoch 3; Illumina 772808, ONT 772817), then `predict_set.sh` per purity (one H100 each, moved ahead of the ONT training with `btop`) and `vcf.sh --set` | fiberseq 772800-772807, Illumina 772809-772816, ONT 772818-772825 (predict, vcf alternating) |
| D | final PASS thresholds (the user, after the tests below): `set_call_threshold.sh` writes the validation call threshold - 0.10 into `COLO829_<P>_SNV_joint_newpipe_e5/best.pth` of Illumina and fiberseq, the calls at the recipe's threshold move to `vcf_chr1_<SET>_cr095`, `vcf.sh --set` without a threshold fills `vcf_chr1_<SET>` again; `check_checkpoint.sh` compares the changed checkpoints with the unchanged copies | 776584 (Illumina), 776589 (fiberseq), vcf 776585-776588, 776590-776593, check 776689, 776690 |

The SNV recipe (`ris/recipes.sh`, one for every platform since f0171af): the 50.4 M model, `--epochs 20 --batch-size
1024 --lr 1e-4 --class-weights pow:0.75 --block v2`, partial somatic labels ignored, `--norm dataset`,
`--site-features`, `--non-fraction 0.25 --keep-non-af 0.100001 --no-importance-weights` (non tensors with AF <= 0.1
sampled 25 % per epoch at loss weight 1), `--best-by postpon --call-recall 0.95` (PASS at the post-PoN validation
recall 0.95 of the callable positives).

- Illumina, fiberseq: the recipe's AF 0.05 tensor floor is not used here: the sets stay at 0.08 (the user), so these
  runs differ from the earlier ones in the recipe only.
- ONT: the sets are built at SNV AF 0.08, the pipeline's ONT floor (manifest `min_af` 0.08 in all four; lowest tensor
  AF on chr22 of P100 and P25 exactly 0.08), so no negative with AF < 0.08 exists to train on. chr22 tensors at AF
  0.08-0.1 / above: P100 non 1,771 / 2,249, germline 107 / 7,293, somatic 14 / 476; P25 non 1,564 / 2,453, germline
  86 / 7,370, somatic 43 / 373.

Not kept: the first ONT submission (771054, repository 0c32edf, when `recipes.sh` still gave ONT the original SNV
recipe: no site features, every negative, PASS at validation recall 0.9) was killed in epoch 2 after the recipe change
(epoch 1: 668 s, 4,594 tensors/s); its partial directory is `$PANSOMA_RUNS/COLO829_ONT_SNV_joint_ontrecipe_killed_20261112`.
An ONT site-feature run under another name (771257-771260, 771262-771265) was submitted and killed before it started.

INDEL is not rerun: `COLO829_<P>_INDEL_joint_ppr8` (`../indel_ppr8_9919d4a_20261026`, step C) already ran the current
INDEL recipe on these sets (the same `args.json` options as `ris/recipes.sh` gives today).

## 10 epochs (the user, 2026-11-12: enough for these dilution models; the recipes keep `--epochs 20`)

- Illumina and fiberseq: the trainings (771059, 771065) were killed after epochs 10 and 15. Post-PoN validation F1 by
  epoch: Illumina 0.897, 0.925, 0.932, **0.936** (epoch 4), then 0.924-0.935 to epoch 10; fiberseq 0.836, 0.859,
  **0.861** (epoch 3), then falling to 0.826 at epoch 15. Their best of epochs 1-10 is the step-C snapshot
  (`COLO829_<P>_SNV_joint_newpipe_e5/best.pth`, byte-identical to the run's `best.pth` when it was killed), so the
  step-C tests are their results; the runs' own vcf jobs were killed.
- ONT: the 20-epoch job (772334) was killed while pending; `joint_run_stop.sh 10` (773816, vcf 773819, 773823, 773827,
  773829) trains with the 20-epoch learning-rate schedule and ends after epoch 10 (`ris/tools/save_epochs.py
  --stop-after-epoch 10`), the same first 10 epochs as the other two platforms; its best of epochs 1-5 is tested too
  (step C).

## Results

`tables.md` (`collect.py`): chr1, benchmark BED, after the PoN, rtg vcfeval, against ClairS-TO and DeepSomatic-TO on
the same BAMs (`../indel_ppr8_9919d4a_20261026/comparators`) and the earlier joint models (all runs finished
2026-11-12 20:00; final thresholds set 21:30).

### Final PASS thresholds (the user, 2026-11-12)

Illumina and fiberseq call at the validation call threshold - 0.10, ONT at the recipe's threshold (the post-PoN
validation recall 0.95 of the callable positives):

| platform | model (checkpoint) | best epoch | recipe's threshold | final PASS threshold |
|---|---|---|---|---|
| Illumina | `COLO829_Illumina_SNV_joint_newpipe_e5/best.pth` | 4 | 0.916548 | **0.816548** |
| fiberseq | `COLO829_fiberseq_SNV_joint_newpipe_e5/best.pth` | 3 | 0.603872 | **0.503872** |
| ONT | `COLO829_ONT_SNV_joint_newpipe/best.pth` | 3 | 0.632432 | **0.632432** (unchanged) |

The Illumina and fiberseq thresholds are written into the checkpoints (`set_call_threshold.sh`, step D):
`val.somatic_at_recall_postpon_callable["0.95"]` is now `{threshold, manual: true, offset: -0.1, note, validation: the
record train wrote}` and `call_threshold_manual` repeats it, so `graph_vcf` / `vcf.sh` without a threshold call at the
final threshold. `call_rule` and `call_recall` are unchanged, so the VCF header still names the rule "validation
post-PoN recall 0.95 of the callable positives" beside the manual threshold; the validation recall at the manual
threshold was not measured. The unchanged checkpoints stay as `COLO829_<P>_SNV_joint_newpipe/best.pth` (byte-identical
to the `_e5` files before the change); `check_checkpoint.sh`: the 603 weight tensors and every other entry of the
changed checkpoints equal the unchanged ones. The chr1 predictions (`test_chr1_<SET>`) are those made before the
change (the probabilities do not depend on the threshold; their `.metrics.json` keep the recipe's threshold).
Directories of `COLO829_<P>_SNV_joint_newpipe_e5`: `vcf_chr1_<SET>` the final threshold (identical TP / FP / FN to
`vcf_chr1_<SET>_tm010`, the test by `vcf.sh THRESHOLD`; `collect.py` asserts it), `vcf_chr1_<SET>_cr095` the recipe's
threshold, `vcf_chr1_<SET>_tm015` - 0.15.

PASS P / R / **F1** at the final thresholds, and F1 of the other callers:

| platform | purity | new pipeline, final | TP / FP / FN | ClairS-TO F1 | DeepSomatic-TO F1 | earlier joint (7fc7396) F1 | new pipeline at the recipe's threshold F1 |
|---|---|---|---|---|---|---|---|
| Illumina (p >= 0.8165) | P100 | 0.588 / 0.733 / **0.653** | 1,790 / 1,253 / 653 | 0.551 | 0.565 | 0.675 | 0.693 |
| | P75 | 0.644 / 0.763 / **0.698** | 1,864 / 1,031 / 579 | 0.661 | 0.500 | 0.840 | 0.781 |
| | P50 | 0.624 / 0.737 / **0.675** | 1,800 / 1,086 / 643 | 0.606 | 0.479 | 0.833 | 0.762 |
| | P25 | 0.591 / 0.689 / **0.636** | 1,684 / 1,165 / 759 | 0.550 | 0.450 | 0.804 | 0.728 |
| fiberseq (p >= 0.5039) | P100 | 0.814 / 0.768 / **0.790** | 1,877 / 430 / 566 | 0.762 | 0.632 | 0.779 | 0.795 |
| | P75 | 0.682 / 0.735 / **0.708** | 1,795 / 836 / 648 | 0.688 | 0.701 | 0.687 | 0.712 |
| | P50 | 0.662 / 0.699 / **0.680** | 1,708 / 872 / 735 | 0.659 | 0.676 | 0.665 | 0.679 |
| | P25 | 0.697 / 0.674 / **0.685** | 1,647 / 716 / 796 | 0.614 | 0.623 | 0.720 | 0.701 |
| ONT (p >= 0.6324) | P100 | 0.795 / 0.720 / **0.756** | 1,759 / 454 / 684 | 0.715 | 0.633 | 0.677 | the final threshold |
| | P75 | 0.728 / 0.701 / **0.714** | 1,712 / 639 / 731 | 0.664 | 0.621 | 0.654 | |
| | P50 | 0.688 / 0.676 / **0.682** | 1,652 / 750 / 791 | 0.630 | 0.609 | 0.639 | |
| | P25 | 0.714 / 0.664 / **0.688** | 1,623 / 649 / 820 | 0.565 | 0.572 | 0.666 | |

- At the final thresholds the new pipeline is above both tumor-only callers in all 12 cells (smallest margins:
  fiberseq P50 / P75 against DeepSomatic-TO, +0.004 / +0.007; Illumina P75 against ClairS-TO, +0.037) and above the
  earlier joint models in 7 of 12 (below: Illumina, -0.022 / -0.142 / -0.158 / -0.168, and fiberseq P25, -0.035).
- Against the recipe's threshold the final one gives Illumina +219 / +14 / +2 / +2 TP for +732 / +584 / +611 / +670
  false calls (F1 -0.040 / -0.082 / -0.087 / -0.092) and fiberseq +30 / +55 / +64 / +12 TP for +75 / +128 / +120 /
  +127 false calls (F1 -0.005 / -0.004 / +0.001 / -0.015).

### At the recipe's threshold (before the change)

PASS F1, P100 / P75 / P50 / P25:

| platform | new pipeline | ClairS-TO | DeepSomatic-TO | earlier joint (7fc7396) |
|---|---|---|---|---|
| Illumina (epoch 4, p >= 0.917) | 0.693 / 0.781 / 0.762 / 0.728 | 0.551 / 0.661 / 0.606 / 0.550 | 0.565 / 0.500 / 0.479 / 0.450 | 0.675 / 0.840 / 0.833 / 0.804 |
| fiberseq (epoch 3, p >= 0.604) | 0.795 / 0.712 / 0.679 / 0.701 | 0.762 / 0.688 / 0.659 / 0.614 | 0.632 / 0.701 / 0.676 / 0.623 | 0.779 / 0.687 / 0.665 / 0.720 |
| ONT (epoch 3, p >= 0.632) | 0.756 / 0.714 / 0.682 / 0.688 | 0.715 / 0.664 / 0.630 / 0.565 | 0.633 / 0.621 / 0.609 / 0.572 | 0.677 / 0.654 / 0.639 / 0.666 |

- At the recipe's threshold the new pipeline is above both tumor-only callers in all 12 cells (smallest margins: fiberseq P50 / P75 against
  DeepSomatic-TO, +0.003 / +0.011) and above the earlier joint models in 8 of 12 (below: Illumina P75 / P50 / P25,
  -0.059 / -0.071 / -0.076, and fiberseq P25, -0.019).
- ONT: post-PoN validation F1 by epoch 0.815, 0.847, **0.861** (epoch 3), 0.854, 0.848, 0.847, 0.816, 0.833, 0.837,
  0.813; the run ended after epoch 10 (78 min in all), its `best.pth` is the epoch-5 snapshot, so the two ONT tests
  are the same checkpoint and give the same numbers.
- Recall is bound by the tensors: of the 2,443 chr1 truth SNVs 82 % (P100) to 71 % (P25) have a tensor in every
  platform (no tensor: AF < 0.08 209-328, fewer than 3 ALT reads 69-321, no candidate 40-217; the sets'
  `somatic.recall.tsv`); the PoN tags 54-62 of those with a tensor. Below the threshold with a tensor: Illumina
  375 / 28 / 8 / 9, fiberseq 99 / 142 / 169 / 48, ONT 177 / 148 / 126 / 58. The 375 of Illumina P100: 358 of the 383
  tensors have somatic as the most probable class (248 at p 0.8-0.917, 110 at 0.5-0.8), 265 at AF 0.85-1.0 and 92 at
  0.4-0.6 (germline-like AF in the pure tumor); the one threshold of the four purities cuts them.

- Illumina against the earlier joint model: recall equal or higher (0.643 / 0.757 / 0.736 / 0.689 against 0.530 /
  0.742 / 0.733 / 0.689), precision lower (0.75-0.81 against 0.93-0.97: 447-521 false calls per purity against
  57-96). The chr1 curves' best F1 at P75 / P50 / P25 are equal (0.830 / 0.832 / 0.804 against 0.841 / 0.834 / 0.805),
  at P100 lower (0.694 against 0.784).

### Lower PASS thresholds (the user, 2026-11-12: Illumina -0.15, fiberseq -0.10, then Illumina -0.10; the tests before
the final thresholds were set)

The same checkpoints and chr1 predictions scored again with `vcf.sh --set SET RUN SNV THRESHOLD SUFFIX`
(`vcf_chr1_<SET>_tm015` / `_tm010`; jobs 774105-774108, 774132-774135; fiberseq -0.15 finished for P100 and P25 only,
774109 / 774112). After the PoN, TP / FP, P / R / F1:

| set | call rule | Illumina p >= 0.767 or fiberseq p >= 0.504 | fiberseq p >= 0.454 |
|---|---|---|---|
| Illumina P100 | 1,571 / 521, 0.751 / 0.643 / 0.693 | 1,853 / 1,609, 0.535 / 0.758 / 0.628 | |
| Illumina P75 | 1,850 / 447, 0.805 / 0.757 / 0.781 | 1,870 / 1,318, 0.587 / 0.765 / 0.664 | |
| Illumina P50 | 1,798 / 475, 0.791 / 0.736 / 0.762 | 1,802 / 1,375, 0.567 / 0.738 / 0.641 | |
| Illumina P25 | 1,682 / 495, 0.773 / 0.689 / 0.728 | 1,687 / 1,484, 0.532 / 0.691 / 0.601 | |
| fiberseq P100 | 1,847 / 355, 0.839 / 0.756 / 0.795 | 1,877 / 430, 0.814 / 0.768 / 0.790 | 1,889 / 472, 0.800 / 0.773 / 0.786 |
| fiberseq P75 | 1,740 / 708, 0.711 / 0.712 / 0.712 | 1,795 / 836, 0.682 / 0.735 / 0.708 | |
| fiberseq P50 | 1,644 / 752, 0.686 / 0.673 / 0.679 | 1,708 / 872, 0.662 / 0.699 / 0.680 | |
| fiberseq P25 | 1,635 / 589, 0.735 / 0.669 / 0.701 | 1,647 / 716, 0.697 / 0.674 / 0.685 | 1,655 / 789, 0.677 / 0.677 / 0.677 |

Illumina: only P100 gains recall (+282 TP; the truth tensors at p 0.77-0.92 with germline-like AF), at +1,088 false
calls; the dilutions gain 4-20 TP for 871-989 false calls. fiberseq: +12 to +64 TP for +75 to +128 false calls, F1
within -0.015 to +0.001.

Illumina at -0.10 (p >= 0.817; `vcf_chr1_<SET>_tm010`, jobs 774218-774221), TP / FP, P / R / F1: P100 1,790 / 1,253,
0.588 / 0.733 / 0.653; P75 1,864 / 1,031, 0.644 / 0.763 / 0.698; P50 1,800 / 1,086, 0.624 / 0.737 / 0.675; P25 1,684 /
1,165, 0.591 / 0.689 / 0.636. Between the call rule and -0.10 the false calls more than double in every purity
(+584 to +732) while only P100 gains truth (+219; dilutions +2 to +14).
