Cosine similarity of SBS96 spectra to the chr1 somatic truth (chr1, benchmark regions; PASS SNV calls after the population filter, TP + FP; Pansoma = six-tumor model on HPRC v2.1, seed 0).

| Tumor | Platform | chr1 truth SNVs | Pansoma, n | Pansoma, cosine to truth | DeepSomatic-TO, n | DeepSomatic-TO, cosine to truth | ClairS-TO, n | ClairS-TO, cosine to truth | Truth vs matched-normal germline, cosine |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HG008T | Illumina | 697 | 1,001 | 0.975 | 1,649 | 0.921 | 1,435 | 0.912 | 0.585 |
| HG008T | PacBio HiFi | 697 | 1,046 | 0.974 | 1,321 | 0.943 | 1,292 | 0.918 | 0.585 |
| HG008T | ONT | 697 | 730 | 0.973 | 1,954 | 0.878 | 1,399 | 0.899 | 0.585 |
| COLO829T | Illumina | 2,443 | 2,384 | 0.994 | 4,779 | 0.913 | 4,335 | 0.913 | 0.428 |
| COLO829T | Fiber-seq | 2,443 | 2,117 | 0.997 | 1,853 | 0.982 | 2,334 | 0.995 | 0.428 |
| COLO829T | ONT | 2,443 | 2,196 | 0.996 | 2,915 | 0.984 | 2,233 | 0.993 | 0.428 |
| HCC1395T | Illumina | 3,440 | 3,608 | 0.998 | 4,083 | 0.993 | 3,667 | 0.993 | 0.524 |
| HCC1395T | PacBio HiFi | 3,440 | 3,363 | 0.997 | 4,165 | 0.994 | 3,512 | 0.994 | 0.524 |
| HCC1395T | ONT | 3,440 | 3,059 | 0.995 | 3,488 | 0.994 | 3,174 | 0.993 | 0.524 |

Cosine similarity of the same call sets to the matched normal's germline spectrum (chr1, same regions) and the share of T>C substitutions.

| Tumor | Platform | Germline SNVs | Pansoma, cosine to germline | DeepSomatic-TO, cosine to germline | ClairS-TO, cosine to germline | T>C share: truth | germline | Pansoma | DeepSomatic-TO | ClairS-TO |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HG008T | Illumina | 296,786 | 0.651 | 0.815 | 0.821 | 0.089 | 0.329 | 0.115 | 0.178 | 0.166 |
| HG008T | PacBio HiFi | 296,786 | 0.674 | 0.763 | 0.827 | 0.089 | 0.329 | 0.121 | 0.160 | 0.186 |
| HG008T | ONT | 296,786 | 0.518 | 0.858 | 0.847 | 0.089 | 0.329 | 0.045 | 0.207 | 0.199 |
| COLO829T | Illumina | 298,057 | 0.439 | 0.596 | 0.620 | 0.075 | 0.327 | 0.072 | 0.125 | 0.073 |
| COLO829T | Fiber-seq | 298,057 | 0.400 | 0.457 | 0.490 | 0.075 | 0.327 | 0.034 | 0.108 | 0.109 |
| COLO829T | ONT | 298,057 | 0.421 | 0.549 | 0.494 | 0.075 | 0.327 | 0.060 | 0.144 | 0.111 |
| HCC1395T | Illumina | 273,569 | 0.537 | 0.597 | 0.595 | 0.114 | 0.332 | 0.124 | 0.145 | 0.127 |
| HCC1395T | PacBio HiFi | 273,569 | 0.509 | 0.592 | 0.589 | 0.114 | 0.332 | 0.111 | 0.137 | 0.137 |
| HCC1395T | ONT | 273,569 | 0.527 | 0.576 | 0.575 | 0.114 | 0.332 | 0.120 | 0.134 | 0.131 |
