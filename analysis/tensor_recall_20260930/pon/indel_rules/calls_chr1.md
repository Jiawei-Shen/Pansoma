# chr1 INDEL calls (HG008T Illumina, chr1 in the BED, 587 truth records) under each INDEL PoN rule

## HG008_Illumina_INDEL_nopartial_b1024

All scored calls: TP 233, FP 5,369; PASS = the run's validation recall-0.9 threshold.

| rule | tagged TP / FP (all calls) | chr1 PASS TP | chr1 PASS FP | P | R | F1 | chr1 best F1 over thresholds |
|---|---:|---:|---:|---:|---:|---:|---:|
| no PoN | 0 / 0 | 217 | 3,835 | 0.054 | 0.370 | **0.094** | 0.150 (P 0.110 R 0.235) |
| repo (allele; gnomAD/CoLoRSdb AF>=1e-4; dbSNP non-somatic; 1000G) | 171 / 4,772 | 52 | 271 | 0.161 | 0.089 | **0.114** | 0.126 (P 0.240 R 0.085) |
| ClairS-TO (gnomAD AF>=1e-3 allele; dbSNP non-somatic allele; 1000G position; CoLoRSdb AF>=1e-3 position) | 174 / 5,005 | 49 | 160 | 0.234 | 0.083 | **0.123** | 0.131 (P 0.196 R 0.099) |
| ClairS-TO AFs, allele match everywhere (the repo rule of a3365c7) | 164 / 4,689 | 59 | 323 | 0.154 | 0.101 | **0.122** | 0.136 (P 0.224 R 0.097) |
| DeepSomatic-like (exact variant key in dbSNP / gnomAD / 1000G; no CoLoRSdb) | 36 / 2,537 | 181 | 1,930 | 0.086 | 0.308 | **0.134** | 0.178 (P 0.151 R 0.216) |
| Mutect2-like (1000G PoN allele + gnomAD AF>=0.01) | 7 / 1,474 | 210 | 2,762 | 0.071 | 0.358 | **0.118** | 0.164 (P 0.128 R 0.228) |
| 1000G PoN only (allele) | 0 / 97 | 217 | 3,770 | 0.054 | 0.370 | **0.095** | 0.151 (P 0.111 R 0.235) |
| repo without CoLoRSdb | 58 / 3,116 | 160 | 1,508 | 0.096 | 0.273 | **0.142** | 0.180 (P 0.164 R 0.199) |
| repo, dbSNP COMMON only | 168 / 4,747 | 55 | 279 | 0.165 | 0.094 | **0.119** | 0.132 (P 0.244 R 0.090) |
| population AF >= 1e-4 only (gnomAD or CoLoRSdb allele) | 168 / 4,736 | 55 | 284 | 0.162 | 0.094 | **0.119** | 0.132 (P 0.242 R 0.090) |
| population AF >= 1e-3 only (gnomAD or CoLoRSdb allele) | 160 / 4,642 | 63 | 343 | 0.155 | 0.107 | **0.127** | 0.142 (P 0.243 R 0.101) |
| population AF >= 0.01 (gnomAD or CoLoRSdb allele) + 1000G + dbSNP COMMON | 83 / 4,035 | 137 | 749 | 0.155 | 0.233 | **0.186** | 0.232 (P 0.309 R 0.186) |
| population AF >= 0.01 only (gnomAD or CoLoRSdb allele) | 80 / 3,996 | 140 | 777 | 0.153 | 0.239 | **0.186** | 0.235 (P 0.306 R 0.191) |
| population AF >= 0.05 (gnomAD or CoLoRSdb allele) + 1000G + dbSNP COMMON | 18 / 2,921 | 199 | 1,602 | 0.110 | 0.339 | **0.167** | 0.220 (P 0.199 R 0.247) |
| population AF >= 0.05 only (gnomAD or CoLoRSdb allele) | 12 / 2,546 | 205 | 1,879 | 0.098 | 0.349 | **0.154** | 0.215 (P 0.187 R 0.254) |
| population AF >= 0.1 (gnomAD or CoLoRSdb allele) + 1000G + dbSNP COMMON | 8 / 2,109 | 209 | 2,253 | 0.085 | 0.356 | **0.137** | 0.186 (P 0.146 R 0.254) |
| population AF >= 0.1 only (gnomAD or CoLoRSdb allele) | 2 / 1,396 | 215 | 2,780 | 0.072 | 0.366 | **0.120** | 0.177 (P 0.142 R 0.235) |

## HG008_Illumina_INDEL_nopartial_b1024_nopc

All scored calls: TP 235, FP 6,433; PASS = the run's validation recall-0.9 threshold.

| rule | tagged TP / FP (all calls) | chr1 PASS TP | chr1 PASS FP | P | R | F1 | chr1 best F1 over thresholds |
|---|---:|---:|---:|---:|---:|---:|---:|
| no PoN | 0 / 0 | 214 | 3,684 | 0.055 | 0.365 | **0.095** | 0.143 (P 0.097 R 0.269) |
| repo (allele; gnomAD/CoLoRSdb AF>=1e-4; dbSNP non-somatic; 1000G) | 173 / 5,633 | 52 | 232 | 0.183 | 0.089 | **0.119** | 0.136 (P 0.333 R 0.085) |
| ClairS-TO (gnomAD AF>=1e-3 allele; dbSNP non-somatic allele; 1000G position; CoLoRSdb AF>=1e-3 position) | 176 / 6,001 | 49 | 135 | 0.266 | 0.083 | **0.127** | 0.134 (P 0.412 R 0.080) |
| ClairS-TO AFs, allele match everywhere (the repo rule of a3365c7) | 166 / 5,540 | 59 | 280 | 0.174 | 0.101 | **0.127** | 0.148 (P 0.308 R 0.097) |
| DeepSomatic-like (exact variant key in dbSNP / gnomAD / 1000G; no CoLoRSdb) | 36 / 2,960 | 179 | 1,694 | 0.096 | 0.305 | **0.146** | 0.189 (P 0.160 R 0.230) |
| Mutect2-like (1000G PoN allele + gnomAD AF>=0.01) | 7 / 1,734 | 207 | 2,499 | 0.076 | 0.353 | **0.126** | 0.173 (P 0.129 R 0.262) |
| 1000G PoN only (allele) | 0 / 124 | 214 | 3,602 | 0.056 | 0.365 | **0.097** | 0.144 (P 0.099 R 0.269) |
| repo without CoLoRSdb | 59 / 3,657 | 158 | 1,310 | 0.108 | 0.269 | **0.154** | 0.200 (P 0.189 R 0.213) |
| repo, dbSNP COMMON only | 170 / 5,597 | 55 | 247 | 0.182 | 0.094 | **0.124** | 0.142 (P 0.333 R 0.090) |
| population AF >= 1e-4 only (gnomAD or CoLoRSdb allele) | 170 / 5,582 | 55 | 251 | 0.180 | 0.094 | **0.123** | 0.142 (P 0.329 R 0.090) |
| population AF >= 1e-3 only (gnomAD or CoLoRSdb allele) | 162 / 5,475 | 63 | 309 | 0.169 | 0.107 | **0.131** | 0.155 (P 0.302 R 0.104) |
| population AF >= 0.01 (gnomAD or CoLoRSdb allele) + 1000G + dbSNP COMMON | 85 / 4,735 | 138 | 693 | 0.166 | 0.235 | **0.195** | 0.246 (P 0.291 R 0.213) |
| population AF >= 0.01 only (gnomAD or CoLoRSdb allele) | 82 / 4,688 | 141 | 724 | 0.163 | 0.240 | **0.194** | 0.246 (P 0.283 R 0.218) |
| population AF >= 0.05 (gnomAD or CoLoRSdb allele) + 1000G + dbSNP COMMON | 19 / 3,410 | 197 | 1,448 | 0.120 | 0.336 | **0.177** | 0.231 (P 0.203 R 0.267) |
| population AF >= 0.05 only (gnomAD or CoLoRSdb allele) | 13 / 2,977 | 203 | 1,764 | 0.103 | 0.346 | **0.159** | 0.218 (P 0.179 R 0.278) |
| population AF >= 0.1 (gnomAD or CoLoRSdb allele) + 1000G + dbSNP COMMON | 8 / 2,488 | 206 | 2,043 | 0.092 | 0.351 | **0.145** | 0.192 (P 0.148 R 0.276) |
| population AF >= 0.1 only (gnomAD or CoLoRSdb allele) | 2 / 1,672 | 212 | 2,651 | 0.074 | 0.361 | **0.123** | 0.173 (P 0.128 R 0.269) |

