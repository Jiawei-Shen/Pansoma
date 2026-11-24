# Response Figure 1.2 (ver4): graph context outside the candidate node masked at inference

Six-tumor current SNV models (`snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe`, seed 0, HPRC v2.1), chr1 SNV, benchmark regions, PASS calls (rule 2: the checkpoint's threshold at validation callable recall 0.95 after the population filter). Full tensor = the existing call sets; own node only = the same checkpoint and threshold with every tensor cell outside the candidate's own node visit set to 0 at inference (the site columns are kept). Δ = full tensor − own node only (difference of the reports' 4-decimal values, rounded to 3 decimals).

## chr1, after the population filter

| Tumor (chr1) | Platform | Truth (chr1) | Full tensor P / R / F1 (chr1) | Own node only P / R / F1 (chr1) | ΔP | ΔR | ΔF1 | Full TP / FP / FN | Own node TP / FP / FN |
|---|---|---|---|---|---|---|---|---|---|
| HG008T | Illumina | 697 | 0.540 / 0.776 / 0.637 | 0.470 / 0.547 / 0.506 | +0.070 | +0.230 | +0.132 | 541 / 460 / 156 | 381 / 429 / 316 |
| HG008T | PacBio HiFi | 697 | 0.542 / 0.814 / 0.651 | 0.577 / 0.610 / 0.593 | −0.035 | +0.204 | +0.057 | 567 / 479 / 130 | 425 / 311 / 272 |
| HG008T | ONT | 697 | 0.626 / 0.656 / 0.640 | 0.638 / 0.311 / 0.418 | −0.012 | +0.344 | +0.222 | 457 / 273 / 240 | 217 / 123 / 480 |
| COLO829T | Illumina | 2,443 | 0.775 / 0.756 / 0.766 | 0.708 / 0.524 / 0.603 | +0.067 | +0.232 | +0.163 | 1,848 / 536 / 595 | 1,281 / 528 / 1,162 |
| COLO829T | Fiber-seq | 2,443 | 0.864 / 0.749 / 0.803 | 0.885 / 0.507 / 0.645 | −0.021 | +0.242 | +0.158 | 1,830 / 287 / 613 | 1,239 / 161 / 1,204 |
| COLO829T | ONT | 2,443 | 0.787 / 0.708 / 0.745 | 0.575 / 0.447 / 0.502 | +0.213 | +0.261 | +0.243 | 1,729 / 467 / 714 | 1,091 / 808 / 1,352 |
| HCC1395T | Illumina | 3,440 | 0.863 / 0.905 / 0.883 | 0.849 / 0.693 / 0.763 | +0.013 | +0.212 | +0.120 | 3,113 / 495 / 327 | 2,383 / 423 / 1,057 |
| HCC1395T | PacBio HiFi | 3,440 | 0.877 / 0.858 / 0.867 | 0.890 / 0.595 / 0.714 | −0.013 | +0.262 | +0.154 | 2,950 / 413 / 490 | 2,048 / 253 / 1,392 |
| HCC1395T | ONT | 3,440 | 0.760 / 0.676 / 0.716 | 0.621 / 0.467 / 0.533 | +0.139 | +0.209 | +0.183 | 2,326 / 733 / 1,114 | 1,606 / 979 / 1,834 |

## chr1, before the population filter

| Tumor (chr1) | Platform | Full tensor P / R / F1 (chr1, before) | Own node only P / R / F1 (chr1, before) | ΔF1 (before) | Full TP / FP / FN (before) | Own node TP / FP / FN (before) |
|---|---|---|---|---|---|---|
| HG008T | Illumina | 0.537 / 0.778 / 0.635 | 0.223 / 0.549 / 0.317 | +0.318 | 542 / 468 / 155 | 383 / 1,338 / 314 |
| HG008T | PacBio HiFi | 0.542 / 0.816 / 0.651 | 0.577 / 0.611 / 0.594 | +0.058 | 569 / 481 / 128 | 426 / 312 / 271 |
| HG008T | ONT | 0.626 / 0.656 / 0.640 | 0.638 / 0.311 / 0.418 | +0.222 | 457 / 273 / 240 | 217 / 123 / 480 |
| COLO829T | Illumina | 0.771 / 0.760 / 0.765 | 0.477 / 0.528 / 0.501 | +0.264 | 1,856 / 550 / 587 | 1,291 / 1,418 / 1,152 |
| COLO829T | Fiber-seq | 0.863 / 0.750 / 0.802 | 0.883 / 0.508 / 0.645 | +0.158 | 1,832 / 291 / 611 | 1,241 / 164 / 1,202 |
| COLO829T | ONT | 0.784 / 0.709 / 0.744 | 0.574 / 0.447 / 0.502 | +0.242 | 1,731 / 477 / 712 | 1,092 / 812 / 1,351 |
| HCC1395T | Illumina | 0.862 / 0.908 / 0.884 | 0.647 / 0.697 / 0.671 | +0.213 | 3,124 / 502 / 316 | 2,397 / 1,309 / 1,043 |
| HCC1395T | PacBio HiFi | 0.877 / 0.860 / 0.868 | 0.890 / 0.598 / 0.715 | +0.153 | 2,959 / 416 / 481 | 2,056 / 255 / 1,384 |
| HCC1395T | ONT | 0.749 / 0.678 / 0.712 | 0.616 / 0.468 / 0.532 | +0.179 | 2,331 / 780 / 1,109 | 1,610 / 1,002 / 1,830 |

## Direction of the change (9 comparisons, after the population filter)

| Metric | Lower with own node only (Δ > 0) | Higher with own node only (Δ < 0) | Equal | Range of Δ | Mean Δ |
|---|---|---|---|---|---|
| Recall | 9 | 0 | 0 | +0.204 to +0.344 | +0.244 |
| Precision | 5 | 4 | 0 | −0.035 to +0.213 | +0.047 |
| F1 | 9 | 0 | 0 | +0.057 to +0.243 | +0.159 |
| F1 before the filter | 9 | 0 | 0 | +0.058 to +0.318 | +0.201 |

Truth alleles lost by the mask per comparison (full TP − own-node TP): 160, 142, 240, 567, 591, 638, 730, 902, 720 (sum 4,690 of 15,361 full-tensor TP). False calls, own node − full: -31, -168, -150, -8, -126, +341, -72, -160, +246 (sum -128 on 4,143).

## Share of the tensor cells the mask zeroes (chr1 samples, mask_fraction.py)

| Set (chr1) | Tensors sampled (all / somatic label) | Tensors with a masked cell (all / somatic) | Non-empty cells zeroed (all / somatic) |
|---|---|---|---|
| HG008T_Illumina_d46 | 3,000 / 825 | 71.3% / 78.1% | 45.0% / 34.3% |
| HG008T_PacBio_d46 | 3,000 / 706 | 39.9% / 75.5% | 65.5% / 32.6% |
| HG008T_ONT_d46 | 3,000 / 704 | 34.0% / 75.1% | 26.6% / 32.9% |
| COLO829T_Illumina_d46_snvaf05 | 3,000 / 2,181 | 75.2% / 72.3% | 42.8% / 28.1% |
| COLO829T_fiberseq_d46_snvaf05 | 3,000 / 2,126 | 32.4% / 71.7% | 57.5% / 29.4% |
| COLO829T_ONT_d46 | 3,000 / 2,019 | 39.9% / 72.0% | 29.3% / 29.6% |
| HCC1395T_Illumina_d46_snvaf05 | 3,000 / 3,000 | 79.7% / 71.3% | 41.4% / 25.7% |
| HCC1395T_PacBio_d46_snvaf05 | 3,000 / 3,000 | 49.2% / 70.9% | 68.4% / 27.4% |
| HCC1395T_ONT_d46 | 3,000 / 2,658 | 46.6% / 70.6% | 19.3% / 27.5% |

## Checks

- Patched code without the mask flag, HG008T PacBio HiFi chr1 (model PacBio): after the filter TP / FP / FN 567 / 481 / 130 vs the existing report 567 / 479 / 130; before the filter 569 / 483 / 128 vs 569 / 481 / 128; threshold 0.8157289028167725 vs 0.8157289028167725: DIFFERENT.
- Truth counts, truth VCF, BED, checkpoint, call rule, threshold and number of predictions equal between the two arms in 9 of 9 finished comparisons.
- PROBLEM: the patched code without the mask does NOT reproduce the existing HG008T PacBio report
