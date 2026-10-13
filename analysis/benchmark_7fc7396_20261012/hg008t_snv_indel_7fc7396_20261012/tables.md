## Training (chr2-22, 5 % of the node blocks validate; chr1 left out)

| run | training tensors (somatic) | validation (somatic) | epochs | best epoch | val somatic AP | val somatic F1 (t) | PASS threshold (val recall 0.9) | train min | tensors/s | GPU peak GiB (allocated) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SNV Illumina | 1,903,802 (7,134) | 118,218 (469) | 20 | 4 | 0.190 | 0.285 (0.672) | 0.251 | 125.0 | 5280 | 49.42 |
| SNV PacBio | 381,661 (7,098) | 18,851 (465) | 20 | 7 | 0.227 | 0.330 (0.727) | 0.314 | 35.8 | 3614 | 48.83 |
| SNV ONT | 509,038 (7,124) | 24,767 (470) | 20 | 8 | 0.210 | 0.301 (0.798) | 0.449 | 47.5 | 3611 | 48.83 |
| INDEL Illumina | 92,108 (3,121) | 6,009 (175) | 12 | 3 | 0.110 | 0.214 (0.392) | 0.119 | 2.6 | 8719 | 23.04 |
| INDEL PacBio | 1,477,288 (3,588) | 98,097 (194) | 12 | 4 | 0.093 | 0.193 (0.392) | 0.040 | 33.3 | 8992 | 23.04 |
| INDEL ONT | 1,890,765 (3,203) | 124,585 (181) | 12 | 5 | 0.055 | 0.168 (0.282) | 0.018 | 43.5 | 8771 | 23.04 |

## chr1 test per tensor and against the truth alleles (predict, `metrics.truth_report`)

| run | test tensors (somatic) | somatic AP | truth alleles | with a tensor (ceiling) | truth F1 at predict's t (P, R) | truth best F1 (P, R) |
|---|---:|---:|---:|---:|---:|---:|
| SNV Illumina | 190,697 (659) | 0.199 | 697 | 659 (0.945) | 0.279 (0.210, 0.418) | 0.285 (0.213, 0.433) |
| SNV PacBio | 32,880 (658) | 0.233 | 697 | 658 (0.944) | 0.307 (0.261, 0.373) | 0.310 (0.258, 0.389) |
| SNV ONT | 44,116 (663) | 0.232 | 697 | 663 (0.951) | 0.297 (0.229, 0.423) | 0.298 (0.249, 0.370) |
| INDEL Illumina | 9,577 (364) | 0.094 | 588 | 317 (0.539) | 0.148 (0.113, 0.218) | 0.153 (0.108, 0.265) |
| INDEL PacBio | 151,524 (481) | 0.076 | 588 | 413 (0.702) | 0.158 (0.113, 0.265) | 0.159 (0.113, 0.269) |
| INDEL ONT | 204,120 (471) | 0.031 | 588 | 408 (0.694) | 0.080 (0.094, 0.070) | 0.091 (0.087, 0.095) |

## chr1 calls, rtg vcfeval, in the GIAB benchmark BED (PASS = validation recall-0.9 threshold)

| run | truth | PASS evaluated (unplaced) | raw P / R / F1 | PoN TP / FP | PoN P / R / F1 | PoN chr1 best F1 at p (P, R; TP / FP) | PoN ceiling | truth tagged by the PoN |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SNV Illumina | 697 | 5,682 (116) | 0.106 / 0.862 / **0.188** | 579 / 703 | 0.452 / 0.831 / **0.585** | 0.607 at 0.388 (0.506, 0.759; 529 / 516) | 0.888 | 61 of 702 |
| SNV PacBio | 697 | 5,540 (211) | 0.108 / 0.857 / **0.191** | 577 / 584 | 0.497 / 0.828 / **0.621** | 0.641 at 0.395 (0.536, 0.798; 556 / 482) | 0.888 | 61 of 702 |
| SNV ONT | 697 | 5,843 (1,409) | 0.102 / 0.855 / **0.182** | 576 / 661 | 0.466 / 0.826 / **0.596** | 0.620 at 0.612 (0.538, 0.730; 509 / 437) | 0.897 | 61 of 702 |
| INDEL Illumina | 587 | 3,596 (563) | 0.059 / 0.359 / **0.101** | 140 / 434 | 0.244 / 0.238 / **0.241** | 0.264 at 0.209 (0.346, 0.213; 125 / 236) | 0.257 | 413 of 635 |
| INDEL PacBio | 587 | 11,759 (2,075) | 0.021 / 0.431 / **0.041** | 150 / 2,630 | 0.054 / 0.256 / **0.089** | 0.260 at 0.259 (0.356, 0.204; 120 / 217) | 0.259 | 413 of 635 |
| INDEL ONT | 587 | 41,052 (7,288) | 0.006 / 0.400 / **0.011** | 134 / 17,268 | 0.008 / 0.228 / **0.015** | 0.140 at 0.168 (0.196, 0.109; 64 / 262) | 0.238 | 413 of 635 |

## chr1 calls, rtg vcfeval, whole chr1, no BED (PASS = validation recall-0.9 threshold)

| run | truth | PASS evaluated (unplaced) | raw P / R / F1 | PoN TP / FP | PoN P / R / F1 | PoN chr1 best F1 at p (P, R; TP / FP) | PoN ceiling | truth tagged by the PoN |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SNV Illumina | 702 | 5,786 (116) | 0.105 / 0.862 / **0.186** | 583 / 735 | 0.442 / 0.831 / **0.577** | 0.602 at 0.394 (0.501, 0.754; 529 / 526) | 0.887 | 61 of 702 |
| SNV PacBio | 702 | 5,624 (211) | 0.107 / 0.856 / **0.190** | 581 / 609 | 0.488 / 0.828 / **0.614** | 0.636 at 0.395 (0.528, 0.798; 560 / 500) | 0.887 | 61 of 702 |
| SNV ONT | 702 | 5,971 (1,409) | 0.101 / 0.855 / **0.180** | 580 / 698 | 0.454 / 0.826 / **0.586** | 0.608 at 0.612 (0.522, 0.728; 511 / 468) | 0.896 | 61 of 702 |
| INDEL Illumina | 635 | 3,889 (563) | 0.058 / 0.357 / **0.100** | 148 / 486 | 0.233 / 0.233 / **0.233** | 0.252 at 0.209 (0.325, 0.206; 131 / 272) | 0.254 | 413 of 635 |
| INDEL PacBio | 635 | 13,300 (2,075) | 0.021 / 0.433 / **0.040** | 161 / 3,129 | 0.049 / 0.254 / **0.082** | 0.247 at 0.259 (0.328, 0.198; 126 / 258) | 0.257 | 413 of 635 |
| INDEL ONT | 635 | 44,712 (7,288) | 0.006 / 0.402 / **0.011** | 143 / 18,825 | 0.007 / 0.225 / **0.015** | 0.134 at 0.175 (0.201, 0.101; 64 / 254) | 0.235 | 413 of 635 |

## Against the earlier runs (chr1, PoN, in the BED: PASS F1, chr1 best F1)

| run | this run (7fc7396) | earlier run | earlier PASS F1 | earlier best F1 | Δ PASS F1 |
|---|---:|---|---:|---:|---:|
| SNV Illumina | 0.585 (best 0.607) | `HG008_Illumina_SNV_h90/vcf_chr1` | 0.592 | 0.612 | -0.007 |
| SNV PacBio | 0.621 (best 0.641) | `HG008_PacBio_SNV_h90/vcf_chr1` | 0.621 | 0.641 | +0.000 |
| SNV ONT | 0.596 (best 0.620) | `HG008_ONT_SNV_h90/vcf_chr1` | 0.596 | 0.620 | +0.000 |
| INDEL Illumina | 0.241 (best 0.264) | `HG008_Illumina_INDEL_small/vcf_chr1_pon7fc7396` | 0.246 | 0.263 | -0.004 |
| INDEL Illumina | 0.241 (best 0.264) | `HG008_Illumina_INDEL_small/vcf_chr1` | 0.218 | 0.243 | +0.023 |
| INDEL PacBio | 0.089 (best 0.260) | `HG008_PacBio_INDEL_small/vcf_chr1_pon7fc7396` | 0.097 | 0.231 | -0.008 |
| INDEL PacBio | 0.089 (best 0.260) | `HG008_PacBio_INDEL_small/vcf_chr1` | 0.089 | 0.220 | -0.000 |
| INDEL ONT | 0.015 (best 0.140) | `HG008_ONT_INDEL_small/vcf_chr1_pon7fc7396` | 0.016 | 0.147 | -0.002 |
| INDEL ONT | 0.015 (best 0.140) | `HG008_ONT_INDEL_small/vcf_chr1` | 0.016 | 0.138 | -0.001 |
