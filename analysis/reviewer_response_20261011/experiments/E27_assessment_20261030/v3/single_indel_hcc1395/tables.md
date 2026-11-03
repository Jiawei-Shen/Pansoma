# Single-tumor INDEL ppr8 checkpoints (RIS) on HCC1395T

chr1 / chr20, benchmark BED, after the population filter, PASS at the checkpoint call rule (post-PoN validation recall 0.8). Checkpoints: RIS `analysis/indel_ppr8_9919d4a_20261026` (3 seeds), predicted on a tequila MIG slice from the repository at 3576902.

## Seed 0 on its own chr1 (check against RIS)

| model | set | P / R / F1 | TP / FP | RIS P / R / F1 | RIS TP / FP |
|---|---|---|---|---|---|
| COLO829_Illumina | COLO829T_Illumina_d46 | 0.107 / 0.221 / 0.144 | 29 / 243 | 0.107 / 0.221 / 0.144 | 29 / 243 |
| COLO829_ONT | COLO829T_ONT_d46 | 0.004 / 0.221 / 0.007 | 29 / 8116 | 0.004 / 0.221 / 0.007 | 29 / 8121 |
| COLO829_fiberseq | COLO829T_fiberseq_d46 | 0.108 / 0.191 / 0.138 | 25 / 207 | 0.108 / 0.191 / 0.138 | 25 / 207 |
| HG008_Illumina | HG008T_Illumina_d46 | 0.343 / 0.213 / 0.263 | 125 / 240 | 0.343 / 0.213 / 0.263 | 125 / 240 |
| HG008_ONT | HG008T_ONT_d46 | 0.016 / 0.198 / 0.029 | 116 / 7274 | 0.016 / 0.198 / 0.029 | 116 / 7280 |
| HG008_PacBio | HG008T_PacBio_d46 | 0.410 / 0.206 / 0.274 | 121 / 174 | 0.410 / 0.206 / 0.274 | 121 / 174 |

## HCC1395T, mean over seeds 0-2 (Supplementary Table R1.7c uses chr1)

| model | set | chrom | n | P / R / F1 (mean) | F1 range | FP per seed |
|---|---|---|---|---|---|---|
| COLO829_Illumina | HCC1395T_Illumina_d46 | chr1 | 3 | 0.229 / 0.744 / 0.349 | 0.320-0.392 | 261, 363, 397 |
| HG008_Illumina | HCC1395T_Illumina_d46 | chr1 | 3 | 0.166 / 0.298 / 0.212 | 0.172-0.246 | 200, 223, 165 |
| COLO829_Illumina | HCC1395T_Illumina_d46 | chr20 | 3 | 0.213 / 0.865 / 0.341 | 0.309-0.383 | 110, 145, 155 |
| HG008_Illumina | HCC1395T_Illumina_d46 | chr20 | 3 | 0.189 / 0.429 / 0.259 | 0.222-0.286 | 74, 102, 54 |
| COLO829_ONT | HCC1395T_ONT_d46 | chr1 | 3 | 0.006 / 0.649 / 0.012 | 0.009-0.015 | 11429, 16003, 18073 |
| HG008_ONT | HCC1395T_ONT_d46 | chr1 | 3 | 0.003 / 0.659 / 0.006 | 0.005-0.007 | 22103, 34827, 26571 |
| COLO829_ONT | HCC1395T_ONT_d46 | chr20 | 3 | 0.006 / 0.722 / 0.013 | 0.010-0.015 | 3985, 4653, 5972 |
| HG008_ONT | HCC1395T_ONT_d46 | chr20 | 3 | 0.004 / 0.746 / 0.007 | 0.006-0.008 | 6836, 11284, 8875 |
| COLO829_fiberseq | HCC1395T_PacBio_d46 | chr1 | 3 | 0.247 / 0.233 / 0.227 | 0.081-0.347 | 124, 34, 88 |
| HG008_PacBio | HCC1395T_PacBio_d46 | chr1 | 3 | 0.241 / 0.283 / 0.259 | 0.235-0.281 | 115, 95, 148 |
| COLO829_fiberseq | HCC1395T_PacBio_d46 | chr20 | 3 | 0.341 / 0.405 / 0.334 | 0.237-0.409 | 59, 10, 41 |
| HG008_PacBio | HCC1395T_PacBio_d46 | chr20 | 3 | 0.351 / 0.603 / 0.440 | 0.411-0.462 | 43, 38, 62 |
