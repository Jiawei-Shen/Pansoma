# COLO829T chr1 by expected truth VAF (computed 2026-10-30, dilution_vaf.py)

Expected VAF = truth VCF VAF (VAF_Ill Illumina, VAF_PB fiberseq / ONT) x measured purity factor (purity.md: 0.66 / 0.40 / 0.20; ONT 0.70 / 0.44 / 0.22). chr1, PASS truth in the BED.

## 1. Candidate stage per tensor set (cells: truth n; tensor / filtered / no candidate)

| set | platform | purity | kind | <0.05 | 0.05-0.10 | 0.10-0.20 | >=0.20 | all | VAF unknown |
|---|---|---|---|---|---|---|---|---|---|
| COLO829T_Illumina (default build (d9, AF 0.08)) | Illumina | P100 | SNV | 303; 10 / 267 / 26 | 146; 48 / 95 / 3 | 142; 120 / 19 / 3 | 1852; 1830 / 3 / 19 | 2443; 2008 / 384 / 51 | 0 |
| COLO829T_Illumina (default build (d9, AF 0.08)) | Illumina | P100 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 131; 57 / 3 / 71 | 131; 57 / 3 / 71 | 0 |
| COLO829T_Illumina_P75 (default build (d9, AF 0.08)) | Illumina | P75 | SNV | 391; 6 / 340 / 45 | 163; 75 / 83 / 5 | 121; 107 / 13 / 1 | 1768; 1748 / 2 / 18 | 2443; 1936 / 438 / 69 | 0 |
| COLO829T_Illumina_P75 (default build (d9, AF 0.08)) | Illumina | P75 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 21; 7 / 1 / 13 | 110; 48 / 3 / 59 | 131; 55 / 4 / 72 | 0 |
| COLO829T_Illumina_P50 (default build (d9, AF 0.08)) | Illumina | P50 | SNV | 506; 19 / 398 / 89 | 130; 68 / 59 / 3 | 266; 253 / 10 / 3 | 1541; 1524 / 1 / 16 | 2443; 1864 / 468 / 111 | 0 |
| COLO829T_Illumina_P50 (default build (d9, AF 0.08)) | Illumina | P50 | INDEL | 0; 0 / 0 / 0 | 12; 4 / 1 / 7 | 60; 25 / 3 / 32 | 59; 24 / 2 / 33 | 131; 53 / 6 / 72 | 0 |
| COLO829T_Illumina_P25 (default build (d9, AF 0.08)) | Illumina | P25 | SNV | 636; 18 / 427 / 191 | 266; 216 / 44 / 6 | 899; 879 / 9 / 11 | 642; 635 / 2 / 5 | 2443; 1748 / 482 / 213 | 0 |
| COLO829T_Illumina_P25 (default build (d9, AF 0.08)) | Illumina | P25 | INDEL | 12; 0 / 4 / 8 | 60; 14 / 14 / 32 | 58; 21 / 4 / 33 | 1; 1 / 0 / 0 | 131; 36 / 22 / 73 | 0 |
| COLO829T_fiberseq (default build (d9, AF 0.08)) | fiberseq | P100 | SNV | 307; 5 / 287 / 15 | 135; 29 / 101 / 5 | 157; 147 / 5 / 5 | 1844; 1825 / 4 / 15 | 2443; 2006 / 397 / 40 | 0 |
| COLO829T_fiberseq (default build (d9, AF 0.08)) | fiberseq | P100 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 131; 67 / 11 / 53 | 131; 67 / 11 / 53 | 0 |
| COLO829T_fiberseq_P75 (default build (d9, AF 0.08)) | fiberseq | P75 | SNV | 411; 6 / 340 / 65 | 160; 86 / 68 / 6 | 122; 115 / 3 / 4 | 1750; 1733 / 1 / 16 | 2443; 1940 / 412 / 91 | 0 |
| COLO829T_fiberseq_P75 (default build (d9, AF 0.08)) | fiberseq | P75 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 5; 4 / 0 / 1 | 126; 61 / 14 / 51 | 131; 65 / 14 / 52 | 0 |
| COLO829T_fiberseq_P50 (default build (d9, AF 0.08)) | fiberseq | P50 | SNV | 509; 16 / 389 / 104 | 142; 86 / 50 / 6 | 313; 304 / 3 / 6 | 1479; 1464 / 2 / 13 | 2443; 1870 / 444 / 129 | 0 |
| COLO829T_fiberseq_P50 (default build (d9, AF 0.08)) | fiberseq | P50 | INDEL | 0; 0 / 0 / 0 | 1; 0 / 0 / 1 | 43; 23 / 8 / 12 | 87; 38 / 10 / 39 | 131; 61 / 18 / 52 | 0 |
| COLO829T_fiberseq_P25 (default build (d9, AF 0.08)) | fiberseq | P25 | SNV | 651; 24 / 430 / 197 | 313; 263 / 43 / 7 | 783; 762 / 11 / 10 | 696; 688 / 5 / 3 | 2443; 1737 / 489 / 217 | 0 |
| COLO829T_fiberseq_P25 (default build (d9, AF 0.08)) | fiberseq | P25 | INDEL | 1; 0 / 0 / 1 | 43; 22 / 8 / 13 | 83; 34 / 10 / 39 | 4; 2 / 0 / 2 | 131; 58 / 18 / 55 | 0 |
| COLO829T_ONT (default build (d9, AF 0.08)) | ONT | P100 | SNV | 307; 20 / 250 / 37 | 135; 30 / 103 / 2 | 157; 120 / 33 / 4 | 1844; 1828 / 7 / 9 | 2443; 1998 / 393 / 52 | 0 |
| COLO829T_ONT (default build (d9, AF 0.08)) | ONT | P100 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 131; 69 / 23 / 39 | 131; 69 / 23 / 39 | 0 |
| COLO829T_ONT_P75 (default build (d9, AF 0.08)) | ONT | P75 | SNV | 400; 9 / 344 / 47 | 153; 62 / 89 / 2 | 133; 107 / 23 / 3 | 1757; 1744 / 6 / 7 | 2443; 1922 / 462 / 59 | 0 |
| COLO829T_ONT_P75 (default build (d9, AF 0.08)) | ONT | P75 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 5; 4 / 0 / 1 | 126; 65 / 23 / 38 | 131; 69 / 23 / 39 | 0 |
| COLO829T_ONT_P50 (default build (d9, AF 0.08)) | ONT | P50 | SNV | 472; 9 / 381 / 82 | 149; 42 / 102 / 5 | 109; 86 / 21 / 2 | 1713; 1701 / 6 / 6 | 2443; 1838 / 510 / 95 | 0 |
| COLO829T_ONT_P50 (default build (d9, AF 0.08)) | ONT | P50 | INDEL | 0; 0 / 0 / 0 | 1; 0 / 0 / 1 | 29; 14 / 8 / 7 | 101; 50 / 19 / 32 | 131; 64 / 27 / 40 | 0 |
| COLO829T_ONT_P25 (default build (d9, AF 0.08)) | ONT | P25 | SNV | 621; 10 / 446 / 165 | 109; 55 / 50 / 4 | 441; 413 / 23 / 5 | 1272; 1258 / 11 / 3 | 2443; 1736 / 530 / 177 | 0 |
| COLO829T_ONT_P25 (default build (d9, AF 0.08)) | ONT | P25 | INDEL | 1; 0 / 0 / 1 | 29; 14 / 7 / 8 | 79; 40 / 12 / 27 | 22; 9 / 3 / 10 | 131; 63 / 22 / 46 | 0 |
| COLO829T_Illumina_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P100 | SNV | 303; 211 / 77 / 15 | 146; 133 / 10 / 3 | 142; 138 / 1 / 3 | 1852; 1830 / 3 / 19 | 2443; 2312 / 91 / 40 | 0 |
| COLO829T_Illumina_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P100 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 131; 57 / 4 / 70 | 131; 57 / 4 / 70 | 0 |
| COLO829T_Illumina_P75_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P75 | SNV | 391; 224 / 139 / 28 | 163; 152 / 6 / 5 | 121; 119 / 1 / 1 | 1768; 1748 / 2 / 18 | 2443; 2243 / 148 / 52 | 0 |
| COLO829T_Illumina_P75_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P75 | INDEL | 0; 0 / 0 / 0 | 0; 0 / 0 / 0 | 21; 8 / 0 / 13 | 110; 49 / 3 / 58 | 131; 57 / 3 / 71 | 0 |
| COLO829T_Illumina_P50_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P50 | SNV | 506; 247 / 182 / 77 | 130; 127 / 1 / 2 | 266; 262 / 1 / 3 | 1541; 1524 / 1 / 16 | 2443; 2160 / 185 / 98 | 0 |
| COLO829T_Illumina_P50_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P50 | INDEL | 0; 0 / 0 / 0 | 12; 5 / 0 / 7 | 60; 28 / 0 / 32 | 59; 24 / 3 / 32 | 131; 57 / 3 / 71 | 0 |
| COLO829T_Illumina_P25_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P25 | SNV | 636; 233 / 234 / 169 | 266; 258 / 4 / 4 | 899; 888 / 0 / 11 | 642; 637 / 0 / 5 | 2443; 2016 / 238 / 189 | 0 |
| COLO829T_Illumina_P25_lowaf_chr1 (low-AF rebuild (d9, AF 0.015, node-alt 0.01)) | Illumina | P25 | INDEL | 12; 5 / 0 / 7 | 60; 26 / 2 / 32 | 58; 23 / 3 / 32 | 1; 1 / 0 / 0 | 131; 55 / 5 / 71 | 0 |

## 2. Final PASS recall after the PoN per expected-VAF bin (rtg pon_calls; cells: R (TP/truth))

| call set | model | platform | purity | <0.05 | 0.05-0.10 | 0.10-0.20 | >=0.20 | P / R / F1 (report.json pon) | unplaced PASS |
|---|---|---|---|---|---|---|---|---|---|
| `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_fiberseq_SNV_joint_7fc7396_gpu4/vcf_chr1_COLO829T_fiberseq/eval_bed` | joint 4-purity SNV, d9 (COLO829T chr2-22 of all purities) | fiberseq | P100 | 0.013 (4/307) | 0.185 (25/135) | 0.879 (138/157) | 0.835 (1540/1844) | 0.881 / 0.699 / 0.779 | 705 |
| `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_fiberseq_SNV_joint_7fc7396_gpu4/vcf_chr1_COLO829T_fiberseq_P75/eval_bed` | joint 4-purity SNV, d9 (COLO829T chr2-22 of all purities) | fiberseq | P75 | 0.012 (5/411) | 0.500 (80/160) | 0.828 (101/122) | 0.787 (1377/1750) | 0.742 / 0.640 / 0.687 | 547 |
| `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_fiberseq_SNV_joint_7fc7396_gpu4/vcf_chr1_COLO829T_fiberseq_P50/eval_bed` | joint 4-purity SNV, d9 (COLO829T chr2-22 of all purities) | fiberseq | P50 | 0.026 (13/509) | 0.542 (77/142) | 0.898 (281/313) | 0.776 (1147/1479) | 0.715 / 0.621 / 0.665 | 626 |
| `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_fiberseq_SNV_joint_7fc7396_gpu4/vcf_chr1_COLO829T_fiberseq_P25/eval_bed` | joint 4-purity SNV, d9 (COLO829T chr2-22 of all purities) | fiberseq | P25 | 0.029 (19/651) | 0.789 (247/313) | 0.917 (718/783) | 0.912 (635/696) | 0.786 / 0.663 / 0.719 | 688 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_control_tq/vcf_chr1_P100/eval_bed` | joint 4-purity SNV, d9, E14b control | ONT | P100 | 0.033 (10/307) | 0.156 (21/135) | 0.567 (89/157) | 0.782 (1442/1844) | 0.858 / 0.639 / 0.733 | 1047 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_control_tq/vcf_chr1_P75/eval_bed` | joint 4-purity SNV, d9, E14b control | ONT | P75 | 0.018 (7/400) | 0.314 (48/153) | 0.654 (87/133) | 0.834 (1465/1757) | 0.798 / 0.658 / 0.721 | 1076 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_control_tq/vcf_chr1_P50/eval_bed` | joint 4-purity SNV, d9, E14b control | ONT | P50 | 0.013 (6/472) | 0.208 (31/149) | 0.633 (69/109) | 0.862 (1477/1713) | 0.740 / 0.648 / 0.691 | 1160 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_control_tq/vcf_chr1_P25/eval_bed` | joint 4-purity SNV, d9, E14b control | ONT | P25 | 0.011 (7/621) | 0.404 (44/109) | 0.873 (385/441) | 0.924 (1175/1272) | 0.779 / 0.659 / 0.714 | 1061 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_readnorm_tq/vcf_chr1_P100/eval_bed` | joint 4-purity SNV, d9, E14b per-read BQ | ONT | P100 | 0.039 (12/307) | 0.185 (25/135) | 0.592 (93/157) | 0.803 (1480/1844) | 0.804 / 0.659 / 0.724 | 3501 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_readnorm_tq/vcf_chr1_P75/eval_bed` | joint 4-purity SNV, d9, E14b per-read BQ | ONT | P75 | 0.018 (7/400) | 0.359 (55/153) | 0.684 (91/133) | 0.829 (1456/1757) | 0.731 / 0.659 / 0.693 | 3288 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_readnorm_tq/vcf_chr1_P50/eval_bed` | joint 4-purity SNV, d9, E14b per-read BQ | ONT | P50 | 0.013 (6/472) | 0.235 (35/149) | 0.651 (71/109) | 0.851 (1458/1713) | 0.690 / 0.643 / 0.666 | 3195 |
| `/scratch/jshen/data/pansoma_net_v2_runs/COLO829T_ONT_SNV_joint_readnorm_tq/vcf_chr1_P25/eval_bed` | joint 4-purity SNV, d9, E14b per-read BQ | ONT | P25 | 0.013 (8/621) | 0.394 (43/109) | 0.875 (386/441) | 0.928 (1181/1272) | 0.742 / 0.662 / 0.700 | 2909 |
| `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_Illumina_SNV_d46_7fc7396/vcf_chr1/eval_bed` | single-tumor COLO829T SNV, d46 (in-sample) | Illumina | P100 | 0.026 (8/303) | 0.301 (44/146) | 0.810 (115/142) | 0.852 (1578/1852) | 0.867 / 0.714 / 0.783 | 516 |
| `/scratch/jshen/data/pansoma_net_v2_runs/multicancer6_Illumina_SNV/vcf_chr1_COLO829T_Illumina_d46_normfix/eval_bed` | multicancer6 SNV, d46 (COLO829T chr2-22 in training) | Illumina | P100 | 0.026 (8/303) | 0.308 (45/146) | 0.817 (116/142) | 0.781 (1447/1852) | 0.880 / 0.661 / 0.755 | 604 |
| `/scratch/jshen/data/pansoma_net_v2_runs/multicancer4_Illumina_SNV/vcf_chr1_COLO829T_Illumina_d46_normfix/eval_bed` | multicancer4 SNV, d46 (COLO829T unseen) | Illumina | P100 | 0.026 (8/303) | 0.308 (45/146) | 0.803 (114/142) | 0.654 (1211/1852) | 0.827 / 0.564 / 0.671 | 922 |
| `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_fiberseq_SNV_d46_7fc7396/vcf_chr1/eval_bed` | single-tumor COLO829T SNV, d46 (in-sample) | fiberseq | P100 | 0.013 (4/307) | 0.193 (26/135) | 0.904 (142/157) | 0.844 (1556/1844) | 0.866 / 0.707 / 0.779 | 257 |
| `/scratch/jshen/data/pansoma_net_v2_runs/multicancer6_PacBio_SNV/vcf_chr1_COLO829T_fiberseq_d46_normfix/eval_bed` | multicancer6 SNV, d46 (COLO829T chr2-22 in training) | fiberseq | P100 | 0.013 (4/307) | 0.193 (26/135) | 0.879 (138/157) | 0.824 (1519/1844) | 0.880 / 0.691 / 0.774 | 801 |
| `/scratch/jshen/data/pansoma_net_v2_runs/multicancer4_PacBio_SNV/vcf_chr1_COLO829T_fiberseq_d46_normfix/eval_bed` | multicancer4 SNV, d46 (COLO829T unseen) | fiberseq | P100 | 0.013 (4/307) | 0.193 (26/135) | 0.911 (143/157) | 0.497 (917/1844) | 0.892 / 0.446 / 0.595 | 1146 |
| `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_ONT_SNV_d46_7fc7396/vcf_chr1/eval_bed` | single-tumor COLO829T SNV, d46 (in-sample) | ONT | P100 | 0.046 (14/307) | 0.185 (25/135) | 0.675 (106/157) | 0.838 (1545/1844) | 0.787 / 0.692 / 0.736 | 4569 |
| `/scratch/jshen/data/pansoma_net_v2_runs/multicancer6_ONT_SNV/vcf_chr1_COLO829T_ONT_d46_normfix/eval_bed` | multicancer6 SNV, d46 (COLO829T chr2-22 in training) | ONT | P100 | 0.046 (14/307) | 0.170 (23/135) | 0.688 (108/157) | 0.810 (1493/1844) | 0.817 / 0.670 / 0.737 | 1318 |
| `/scratch/jshen/data/pansoma_net_v2_runs/multicancer4_ONT_SNV/vcf_chr1_COLO829T_ONT_d46_normfix/eval_bed` | multicancer4 SNV, d46 (COLO829T unseen) | ONT | P100 | 0.046 (14/307) | 0.170 (23/135) | 0.669 (105/157) | 0.331 (610/1844) | 0.497 / 0.308 / 0.380 | 2211 |

## 3. E15b: HapMap 200-row SNV model (HapMap_Illumina_SNV_full) on COLO829T Illumina P100, low-AF rebuild, tensor level, no PoN (cells: R (called/truth))

| threshold | value | <0.05 | 0.05-0.10 | 0.10-0.20 | >=0.20 | all | truth-level P / R / F1 (metrics.json) |
|---|---|---|---|---|---|---|---|
| pipeline (val recall 0.9) | 0.9876 | 0.432 (131/303) | 0.747 (109/146) | 0.521 (74/142) | 0.005 (9/1852) | 0.132 (323/2443) | 0.325 / 0.132 / 0.188 |
| threshold_t | 0.8714 | 0.647 (196/303) | 0.863 (126/146) | 0.852 (121/142) | 0.010 (19/1852) | 0.189 (462/2443) | 0.202 / 0.189 / 0.195 |
