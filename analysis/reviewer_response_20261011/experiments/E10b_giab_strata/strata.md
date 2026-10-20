# E10b: GIAB v3.5 stratifications (chr1, PASS after the PoN, in the benchmark BED; 1-Mb blocks, 2,000 resamples)

## 1. Summary: HPRC graph - GRCh38-only (linear) graph, paired F1 difference

Per stratum and variant kind: the 12 pairs (HG008T / COLO829T x Illumina / PacBio-fiberseq / ONT x v1.1 d9 / v2.1 d46, where both arms exist). Truth = truth records of the stratum in the BED (range over the six sets).

| stratum | chr1 Mb | kind | pairs | truth per set | dF1 range | dF1 > 0 | CI above 0 | CI below 0 |
|---|---:|---|---:|---|---|---:|---:|---:|
| not in any difficult region | 183.7 | SNV | 12 | 517-2000 | -0.001 to +0.176 | 11 | 7 | 0 |
| all difficult regions | 46.8 | SNV | 12 | 180-443 | +0.054 to +0.483 | 12 | 10 | 0 |
| low mappability or segdup | 23.2 | SNV | 12 | 67-176 | +0.033 to +0.528 | 12 | 10 | 0 |
| low mappability | 19.8 | SNV | 12 | 50-144 | +0.034 to +0.548 | 12 | 10 | 0 |
| segmental duplications | 14.8 | SNV | 12 | 48-95 | -0.063 to +0.331 | 11 | 9 | 0 |
| tandem repeats and homopolymers | 17.4 | SNV | 12 | 87-166 | +0.031 to +0.513 | 12 | 11 | 0 |
| homopolymers >= 7 bp | 6.7 | SNV | 12 | 42-70 | +0.025 to +0.160 | 12 | 5 | 0 |
| tandem repeats <= 50 bp | 2.3 | SNV | 12 | 14-20 | -0.080 to +0.342 | 9 | 4 | 0 |
| satellites | 3.4 | SNV | 8 | 2-11 | +0.002 to +0.490 | 8 | 5 | 0 |
| L1H > 500 bp | 0.3 | SNV | 7 | 1-3 | -0.250 to +0.214 | 3 | 0 | 0 |
| GC < 25% or > 65% | 14.2 | SNV | 12 | 59-157 | +0.008 to +0.232 | 12 | 8 | 0 |
| not in any difficult region | 183.7 | INDEL | 12 | 7-35 | -0.001 to +0.280 | 11 | 9 | 0 |
| all difficult regions | 46.8 | INDEL | 12 | 124-552 | -0.053 to +0.090 | 4 | 4 | 7 |
| low mappability or segdup | 23.2 | INDEL | 12 | 29-41 | -0.078 to +0.125 | 4 | 2 | 3 |
| low mappability | 19.8 | INDEL | 12 | 26-28 | -0.080 to +0.123 | 4 | 2 | 3 |
| segmental duplications | 14.8 | INDEL | 12 | 22-30 | -0.077 to +0.056 | 2 | 1 | 3 |
| tandem repeats and homopolymers | 17.4 | INDEL | 12 | 108-547 | -0.053 to +0.089 | 4 | 4 | 7 |
| homopolymers >= 7 bp | 6.7 | INDEL | 12 | 75-465 | -0.065 to +0.083 | 4 | 4 | 7 |
| tandem repeats <= 50 bp | 2.3 | INDEL | 12 | 21-36 | -0.033 to +0.130 | 5 | 2 | 0 |
| L1H > 500 bp | 0.3 | INDEL | 5 | 1-1 | +0.011 to +0.714 | 5 | 1 | 0 |
| GC < 25% or > 65% | 14.2 | INDEL | 12 | 14-50 | -0.040 to +0.088 | 5 | 2 | 7 |

## 2. Paired differences (HPRC - linear graph)

| set | kind | arm | stratum | truth | dP [95% CI] | dR [95% CI] | dF1 [95% CI] |
|---|---|---|---|---:|---|---|---|
| COLO829T_Illumina | INDEL | d9 | not in any difficult region | 7 | +0.067 [+0.020, +0.130] | +0.429 [+0.000, +0.857] | +0.119 [+0.039, +0.217] |
| COLO829T_Illumina | INDEL | d9 | all difficult regions | 124 | +0.042 [+0.023, +0.061] | +0.032 [-0.027, +0.090] | +0.063 [+0.034, +0.090] |
| COLO829T_Illumina | INDEL | d9 | low mappability or segdup | 29 | +0.113 [+0.024, +0.220] | -0.034 [-0.088, +0.000] | +0.125 [+0.028, +0.242] |
| COLO829T_Illumina | INDEL | d9 | low mappability | 26 | +0.115 [+0.023, +0.232] | +0.000 [+0.000, +0.000] | +0.123 [+0.031, +0.242] |
| COLO829T_Illumina | INDEL | d9 | segmental duplications | 22 | +0.066 [-0.012, +0.195] | -0.045 [-0.111, +0.000] | +0.056 [+0.001, +0.188] |
| COLO829T_Illumina | INDEL | d9 | tandem repeats and homopolymers | 108 | +0.038 [+0.021, +0.056] | +0.046 [-0.019, +0.108] | +0.059 [+0.032, +0.084] |
| COLO829T_Illumina | INDEL | d9 | homopolymers >= 7 bp | 75 | +0.026 [+0.011, +0.042] | +0.067 [-0.013, +0.140] | +0.041 [+0.016, +0.066] |
| COLO829T_Illumina | INDEL | d9 | tandem repeats <= 50 bp | 21 | +0.026 [-0.036, +0.113] | -0.048 [-0.158, +0.000] | +0.011 [-0.048, +0.082] |
| COLO829T_Illumina | INDEL | d9 | satellites | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| COLO829T_Illumina | INDEL | d9 | L1H > 500 bp | 1 | +0.833 [+0.250, +0.909] | +0.000 [+0.000, +0.000] | +0.714 [+0.143, +0.833] |
| COLO829T_Illumina | INDEL | d9 | GC < 25% or > 65% | 14 | +0.040 [+0.000, +0.087] | +0.000 [+0.000, +0.000] | +0.062 [+0.016, +0.129] |
| COLO829T_Illumina | INDEL | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | INDEL | d46 | not in any difficult region | 7 | +0.108 [+0.037, +0.197] | +0.429 [+0.000, +0.857] | +0.185 [+0.072, +0.311] |
| COLO829T_Illumina | INDEL | d46 | all difficult regions | 124 | +0.051 [+0.029, +0.076] | +0.000 [-0.064, +0.065] | +0.070 [+0.037, +0.103] |
| COLO829T_Illumina | INDEL | d46 | low mappability or segdup | 29 | +0.054 [-0.005, +0.134] | -0.103 [-0.238, +0.000] | +0.061 [+0.005, +0.156] |
| COLO829T_Illumina | INDEL | d46 | low mappability | 26 | +0.046 [-0.006, +0.130] | -0.077 [-0.238, +0.000] | +0.050 [+0.010, +0.146] |
| COLO829T_Illumina | INDEL | d46 | segmental duplications | 22 | +0.026 [-0.027, +0.130] | -0.091 [-0.222, +0.000] | +0.013 [-0.018, +0.125] |
| COLO829T_Illumina | INDEL | d46 | tandem repeats and homopolymers | 108 | +0.049 [+0.027, +0.074] | +0.028 [-0.044, +0.096] | +0.071 [+0.038, +0.104] |
| COLO829T_Illumina | INDEL | d46 | homopolymers >= 7 bp | 75 | +0.039 [+0.017, +0.065] | +0.040 [-0.048, +0.123] | +0.056 [+0.022, +0.093] |
| COLO829T_Illumina | INDEL | d46 | tandem repeats <= 50 bp | 21 | +0.047 [-0.033, +0.170] | -0.048 [-0.158, +0.000] | +0.023 [-0.043, +0.105] |
| COLO829T_Illumina | INDEL | d46 | satellites | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| COLO829T_Illumina | INDEL | d46 | L1H > 500 bp | 1 | - [-, -] | -1.000 [-1.000, -1.000] | - [-, -] |
| COLO829T_Illumina | INDEL | d46 | GC < 25% or > 65% | 14 | +0.059 [+0.000, +0.135] | +0.000 [+0.000, +0.000] | +0.088 [+0.023, +0.186] |
| COLO829T_Illumina | INDEL | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | SNV | d9 | not in any difficult region | 2000 | -0.010 [-0.022, +0.003] | +0.127 [+0.105, +0.148] | +0.082 [+0.066, +0.098] |
| COLO829T_Illumina | SNV | d9 | all difficult regions | 443 | +0.036 [-0.111, +0.102] | +0.063 [+0.023, +0.102] | +0.054 [-0.060, +0.120] |
| COLO829T_Illumina | SNV | d9 | low mappability or segdup | 176 | +0.019 [-0.056, +0.077] | +0.034 [-0.040, +0.106] | +0.033 [-0.051, +0.098] |
| COLO829T_Illumina | SNV | d9 | low mappability | 144 | +0.020 [-0.037, +0.080] | +0.056 [-0.021, +0.129] | +0.034 [-0.032, +0.107] |
| COLO829T_Illumina | SNV | d9 | segmental duplications | 95 | -0.047 [-0.121, +0.198] | -0.042 [-0.156, +0.065] | -0.063 [-0.140, +0.142] |
| COLO829T_Illumina | SNV | d9 | tandem repeats and homopolymers | 166 | +0.018 [-0.091, +0.058] | +0.078 [+0.023, +0.132] | +0.031 [-0.087, +0.085] |
| COLO829T_Illumina | SNV | d9 | homopolymers >= 7 bp | 70 | +0.136 [-0.136, +0.255] | +0.143 [+0.039, +0.244] | +0.140 [-0.005, +0.245] |
| COLO829T_Illumina | SNV | d9 | tandem repeats <= 50 bp | 20 | -0.137 [-0.400, +0.118] | -0.050 [-0.214, +0.125] | -0.080 [-0.243, +0.069] |
| COLO829T_Illumina | SNV | d9 | satellites | 11 | +0.001 [+0.000, +0.011] | +0.000 [+0.000, +0.000] | +0.002 [+0.000, +0.020] |
| COLO829T_Illumina | SNV | d9 | L1H > 500 bp | 3 | +0.333 [+0.000, +1.000] | +0.333 [+0.000, +1.000] | - [-, -] |
| COLO829T_Illumina | SNV | d9 | GC < 25% or > 65% | 157 | -0.061 [-0.161, +0.030] | +0.051 [-0.013, +0.115] | +0.008 [-0.054, +0.068] |
| COLO829T_Illumina | SNV | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | SNV | d46 | not in any difficult region | 2000 | -0.027 [-0.041, -0.014] | +0.154 [+0.129, +0.176] | +0.091 [+0.072, +0.110] |
| COLO829T_Illumina | SNV | d46 | all difficult regions | 443 | +0.665 [+0.310, +0.738] | +0.034 [-0.014, +0.080] | +0.483 [+0.157, +0.567] |
| COLO829T_Illumina | SNV | d46 | low mappability or segdup | 176 | +0.659 [+0.440, +0.760] | +0.017 [-0.073, +0.110] | +0.528 [+0.253, +0.606] |
| COLO829T_Illumina | SNV | d46 | low mappability | 144 | +0.674 [+0.457, +0.784] | +0.021 [-0.074, +0.118] | +0.548 [+0.287, +0.631] |
| COLO829T_Illumina | SNV | d46 | segmental duplications | 95 | +0.520 [+0.326, +0.641] | -0.011 [-0.148, +0.124] | +0.331 [+0.116, +0.449] |
| COLO829T_Illumina | SNV | d46 | tandem repeats and homopolymers | 166 | +0.646 [+0.405, +0.753] | +0.012 [-0.054, +0.077] | +0.513 [+0.232, +0.587] |
| COLO829T_Illumina | SNV | d46 | homopolymers >= 7 bp | 70 | +0.291 [-0.031, +0.476] | +0.043 [-0.076, +0.157] | +0.160 [-0.016, +0.306] |
| COLO829T_Illumina | SNV | d46 | tandem repeats <= 50 bp | 20 | +0.000 [-0.200, +0.222] | +0.000 [-0.143, +0.150] | +0.000 [-0.136, +0.132] |
| COLO829T_Illumina | SNV | d46 | satellites | 11 | +0.427 [+0.249, +0.999] | +0.000 [+0.000, +0.000] | +0.330 [+0.117, +0.498] |
| COLO829T_Illumina | SNV | d46 | L1H > 500 bp | 3 | +0.200 [+0.000, +0.667] | +0.333 [+0.000, +1.000] | - [-, -] |
| COLO829T_Illumina | SNV | d46 | GC < 25% or > 65% | 157 | +0.005 [-0.060, +0.070] | +0.025 [-0.038, +0.092] | +0.019 [-0.035, +0.075] |
| COLO829T_Illumina | SNV | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | INDEL | d9 | not in any difficult region | 7 | -0.000 [-0.005, +0.004] | +0.000 [-0.429, +0.444] | -0.001 [-0.009, +0.008] |
| COLO829T_ONT | INDEL | d9 | all difficult regions | 124 | -0.000 [-0.001, +0.000] | -0.040 [-0.100, +0.017] | -0.001 [-0.002, +0.000] |
| COLO829T_ONT | INDEL | d9 | low mappability or segdup | 29 | -0.001 [-0.002, +0.000] | -0.069 [-0.176, +0.000] | -0.001 [-0.004, +0.000] |
| COLO829T_ONT | INDEL | d9 | low mappability | 26 | -0.000 [-0.001, +0.000] | -0.038 [-0.100, +0.000] | -0.001 [-0.003, +0.000] |
| COLO829T_ONT | INDEL | d9 | segmental duplications | 22 | -0.001 [-0.004, +0.000] | -0.091 [-0.222, +0.000] | -0.002 [-0.007, +0.001] |
| COLO829T_ONT | INDEL | d9 | tandem repeats and homopolymers | 108 | -0.000 [-0.001, +0.000] | -0.028 [-0.090, +0.030] | -0.001 [-0.001, +0.000] |
| COLO829T_ONT | INDEL | d9 | homopolymers >= 7 bp | 75 | -0.000 [-0.001, +0.000] | -0.013 [-0.083, +0.056] | -0.000 [-0.001, +0.000] |
| COLO829T_ONT | INDEL | d9 | tandem repeats <= 50 bp | 21 | -0.005 [-0.015, +0.001] | -0.095 [-0.250, +0.000] | -0.010 [-0.028, +0.002] |
| COLO829T_ONT | INDEL | d9 | satellites | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| COLO829T_ONT | INDEL | d9 | L1H > 500 bp | 1 | +0.066 [+0.000, +0.262] | +0.000 [+0.000, +0.000] | +0.107 [-0.010, +0.361] |
| COLO829T_ONT | INDEL | d9 | GC < 25% or > 65% | 14 | +0.001 [-0.000, +0.004] | +0.071 [+0.000, +0.250] | +0.002 [-0.001, +0.009] |
| COLO829T_ONT | INDEL | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | INDEL | d46 | not in any difficult region | 7 | +0.001 [-0.003, +0.005] | -0.143 [-0.500, +0.000] | +0.002 [-0.006, +0.009] |
| COLO829T_ONT | INDEL | d46 | all difficult regions | 124 | -0.001 [-0.001, -0.000] | -0.048 [-0.106, +0.008] | -0.002 [-0.003, -0.001] |
| COLO829T_ONT | INDEL | d46 | low mappability or segdup | 29 | -0.001 [-0.003, -0.000] | -0.103 [-0.226, +0.000] | -0.003 [-0.006, -0.000] |
| COLO829T_ONT | INDEL | d46 | low mappability | 26 | -0.001 [-0.003, -0.000] | -0.077 [-0.185, +0.000] | -0.002 [-0.006, -0.000] |
| COLO829T_ONT | INDEL | d46 | segmental duplications | 22 | -0.002 [-0.005, +0.000] | -0.136 [-0.286, +0.000] | -0.004 [-0.009, -0.000] |
| COLO829T_ONT | INDEL | d46 | tandem repeats and homopolymers | 108 | -0.001 [-0.001, -0.000] | -0.028 [-0.083, +0.022] | -0.001 [-0.002, -0.001] |
| COLO829T_ONT | INDEL | d46 | homopolymers >= 7 bp | 75 | -0.001 [-0.001, -0.000] | -0.027 [-0.092, +0.034] | -0.001 [-0.002, -0.000] |
| COLO829T_ONT | INDEL | d46 | tandem repeats <= 50 bp | 21 | -0.003 [-0.010, +0.000] | -0.048 [-0.158, +0.000] | -0.006 [-0.019, +0.001] |
| COLO829T_ONT | INDEL | d46 | satellites | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| COLO829T_ONT | INDEL | d46 | L1H > 500 bp | 1 | +0.006 [-0.050, +0.089] | +0.000 [+0.000, +0.000] | +0.011 [-0.089, +0.143] |
| COLO829T_ONT | INDEL | d46 | GC < 25% or > 65% | 14 | -0.000 [-0.001, +0.000] | +0.000 [+0.000, +0.000] | -0.001 [-0.002, -0.000] |
| COLO829T_ONT | INDEL | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | SNV | d9 | not in any difficult region | 2000 | +0.179 [+0.156, +0.205] | +0.172 [+0.145, +0.198] | +0.176 [+0.157, +0.197] |
| COLO829T_ONT | SNV | d9 | all difficult regions | 443 | +0.483 [+0.344, +0.565] | +0.097 [+0.045, +0.146] | +0.427 [+0.277, +0.503] |
| COLO829T_ONT | SNV | d9 | low mappability or segdup | 176 | +0.439 [+0.320, +0.572] | +0.085 [-0.006, +0.175] | +0.438 [+0.323, +0.532] |
| COLO829T_ONT | SNV | d9 | low mappability | 144 | +0.435 [+0.311, +0.576] | +0.104 [+0.006, +0.200] | +0.451 [+0.330, +0.553] |
| COLO829T_ONT | SNV | d9 | segmental duplications | 95 | +0.278 [+0.172, +0.435] | +0.011 [-0.134, +0.144] | +0.284 [+0.189, +0.388] |
| COLO829T_ONT | SNV | d9 | tandem repeats and homopolymers | 166 | +0.442 [+0.233, +0.512] | +0.066 [-0.006, +0.137] | +0.427 [+0.183, +0.497] |
| COLO829T_ONT | SNV | d9 | homopolymers >= 7 bp | 70 | +0.022 [-0.098, +0.133] | +0.129 [+0.014, +0.259] | +0.064 [-0.033, +0.173] |
| COLO829T_ONT | SNV | d9 | tandem repeats <= 50 bp | 20 | +0.469 [+0.247, +0.694] | +0.100 [-0.095, +0.312] | +0.299 [+0.119, +0.484] |
| COLO829T_ONT | SNV | d9 | satellites | 11 | +0.496 [-0.004, +0.997] | -0.273 [-0.625, +0.000] | +0.258 [+0.068, +0.562] |
| COLO829T_ONT | SNV | d9 | L1H > 500 bp | 3 | -0.200 [-0.750, +0.000] | -0.333 [-1.000, +0.000] | -0.250 [-0.600, +0.000] |
| COLO829T_ONT | SNV | d9 | GC < 25% or > 65% | 157 | +0.326 [+0.138, +0.468] | +0.121 [+0.060, +0.187] | +0.232 [+0.122, +0.340] |
| COLO829T_ONT | SNV | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | SNV | d46 | not in any difficult region | 2000 | +0.181 [+0.156, +0.207] | +0.168 [+0.139, +0.196] | +0.175 [+0.154, +0.196] |
| COLO829T_ONT | SNV | d46 | all difficult regions | 443 | +0.450 [+0.270, +0.619] | +0.065 [+0.012, +0.114] | +0.395 [+0.235, +0.495] |
| COLO829T_ONT | SNV | d46 | low mappability or segdup | 176 | +0.323 [+0.177, +0.576] | +0.057 [-0.045, +0.158] | +0.354 [+0.207, +0.495] |
| COLO829T_ONT | SNV | d46 | low mappability | 144 | +0.319 [+0.172, +0.591] | +0.076 [-0.036, +0.181] | +0.364 [+0.213, +0.518] |
| COLO829T_ONT | SNV | d46 | segmental duplications | 95 | +0.175 [+0.065, +0.460] | +0.021 [-0.143, +0.180] | +0.209 [+0.081, +0.378] |
| COLO829T_ONT | SNV | d46 | tandem repeats and homopolymers | 166 | +0.596 [+0.377, +0.687] | +0.036 [-0.034, +0.108] | +0.485 [+0.243, +0.556] |
| COLO829T_ONT | SNV | d46 | homopolymers >= 7 bp | 70 | +0.155 [+0.029, +0.286] | +0.057 [-0.049, +0.176] | +0.106 [+0.004, +0.212] |
| COLO829T_ONT | SNV | d46 | tandem repeats <= 50 bp | 20 | +0.585 [+0.347, +0.813] | +0.100 [-0.095, +0.312] | +0.342 [+0.159, +0.533] |
| COLO829T_ONT | SNV | d46 | satellites | 11 | +0.662 [-0.004, +0.998] | -0.273 [-0.625, +0.000] | +0.277 [+0.023, +0.563] |
| COLO829T_ONT | SNV | d46 | L1H > 500 bp | 3 | -0.200 [-0.750, +0.000] | -0.333 [-1.000, +0.000] | -0.250 [-0.600, +0.000] |
| COLO829T_ONT | SNV | d46 | GC < 25% or > 65% | 157 | +0.340 [+0.166, +0.469] | +0.089 [+0.029, +0.151] | +0.221 [+0.112, +0.321] |
| COLO829T_ONT | SNV | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | INDEL | d9 | not in any difficult region | 7 | +0.036 [+0.010, +0.075] | +0.143 [+0.000, +0.500] | +0.064 [+0.019, +0.126] |
| COLO829T_fiberseq | INDEL | d9 | all difficult regions | 124 | -0.011 [-0.018, -0.006] | -0.040 [-0.113, +0.028] | -0.021 [-0.032, -0.010] |
| COLO829T_fiberseq | INDEL | d9 | low mappability or segdup | 29 | +0.002 [-0.005, +0.010] | -0.034 [-0.088, +0.000] | +0.003 [-0.010, +0.018] |
| COLO829T_fiberseq | INDEL | d9 | low mappability | 26 | +0.004 [-0.005, +0.016] | -0.038 [-0.100, +0.000] | +0.006 [-0.010, +0.026] |
| COLO829T_fiberseq | INDEL | d9 | segmental duplications | 22 | -0.001 [-0.011, +0.007] | -0.045 [-0.111, +0.000] | -0.002 [-0.019, +0.013] |
| COLO829T_fiberseq | INDEL | d9 | tandem repeats and homopolymers | 108 | -0.010 [-0.017, -0.005] | -0.037 [-0.120, +0.042] | -0.019 [-0.031, -0.008] |
| COLO829T_fiberseq | INDEL | d9 | homopolymers >= 7 bp | 75 | -0.016 [-0.026, -0.008] | -0.013 [-0.119, +0.086] | -0.029 [-0.046, -0.014] |
| COLO829T_fiberseq | INDEL | d9 | tandem repeats <= 50 bp | 21 | -0.020 [-0.057, +0.006] | -0.095 [-0.250, +0.000] | -0.033 [-0.085, +0.009] |
| COLO829T_fiberseq | INDEL | d9 | satellites | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| COLO829T_fiberseq | INDEL | d9 | L1H > 500 bp | 1 | +0.190 [-0.050, +0.857] | +0.000 [+0.000, +0.000] | +0.250 [-0.100, +0.778] |
| COLO829T_fiberseq | INDEL | d9 | GC < 25% or > 65% | 14 | -0.008 [-0.019, -0.000] | -0.071 [-0.250, +0.000] | -0.014 [-0.036, -0.001] |
| COLO829T_fiberseq | INDEL | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | INDEL | d46 | not in any difficult region | 7 | +0.033 [-0.003, +0.079] | +0.000 [-0.429, +0.429] | +0.058 [-0.006, +0.133] |
| COLO829T_fiberseq | INDEL | d46 | all difficult regions | 124 | -0.013 [-0.020, -0.007] | -0.040 [-0.113, +0.031] | -0.023 [-0.035, -0.012] |
| COLO829T_fiberseq | INDEL | d46 | low mappability or segdup | 29 | +0.010 [-0.005, +0.032] | -0.069 [-0.160, +0.000] | +0.015 [-0.010, +0.049] |
| COLO829T_fiberseq | INDEL | d46 | low mappability | 26 | +0.011 [-0.008, +0.039] | -0.077 [-0.185, +0.000] | +0.015 [-0.013, +0.056] |
| COLO829T_fiberseq | INDEL | d46 | segmental duplications | 22 | +0.001 [-0.016, +0.027] | -0.091 [-0.222, +0.000] | -0.001 [-0.020, +0.042] |
| COLO829T_fiberseq | INDEL | d46 | tandem repeats and homopolymers | 108 | -0.011 [-0.018, -0.005] | -0.028 [-0.108, +0.050] | -0.020 [-0.032, -0.010] |
| COLO829T_fiberseq | INDEL | d46 | homopolymers >= 7 bp | 75 | -0.018 [-0.027, -0.009] | -0.013 [-0.120, +0.087] | -0.032 [-0.049, -0.016] |
| COLO829T_fiberseq | INDEL | d46 | tandem repeats <= 50 bp | 21 | -0.009 [-0.035, +0.010] | -0.048 [-0.158, +0.000] | -0.016 [-0.056, +0.013] |
| COLO829T_fiberseq | INDEL | d46 | satellites | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| COLO829T_fiberseq | INDEL | d46 | L1H > 500 bp | 1 | +0.357 [+0.000, +0.889] | +0.000 [+0.000, +0.000] | +0.417 [+0.000, +0.818] |
| COLO829T_fiberseq | INDEL | d46 | GC < 25% or > 65% | 14 | -0.008 [-0.020, -0.001] | -0.071 [-0.250, +0.000] | -0.015 [-0.036, -0.002] |
| COLO829T_fiberseq | INDEL | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | SNV | d9 | not in any difficult region | 2000 | -0.026 [-0.038, -0.015] | +0.106 [+0.083, +0.126] | +0.057 [+0.042, +0.072] |
| COLO829T_fiberseq | SNV | d9 | all difficult regions | 443 | +0.102 [-0.175, +0.397] | +0.065 [+0.017, +0.113] | +0.091 [-0.101, +0.266] |
| COLO829T_fiberseq | SNV | d9 | low mappability or segdup | 176 | +0.100 [-0.137, +0.508] | +0.051 [-0.031, +0.135] | +0.105 [-0.123, +0.377] |
| COLO829T_fiberseq | SNV | d9 | low mappability | 144 | +0.091 [-0.141, +0.529] | +0.076 [-0.013, +0.161] | +0.106 [-0.128, +0.410] |
| COLO829T_fiberseq | SNV | d9 | segmental duplications | 95 | +0.019 [-0.177, +0.439] | -0.032 [-0.144, +0.082] | +0.018 [-0.189, +0.319] |
| COLO829T_fiberseq | SNV | d9 | tandem repeats and homopolymers | 166 | +0.329 [+0.048, +0.476] | +0.042 [-0.042, +0.122] | +0.216 [+0.048, +0.330] |
| COLO829T_fiberseq | SNV | d9 | homopolymers >= 7 bp | 70 | -0.112 [-0.250, +0.027] | +0.129 [+0.021, +0.243] | +0.025 [-0.079, +0.132] |
| COLO829T_fiberseq | SNV | d9 | tandem repeats <= 50 bp | 20 | +0.305 [-0.036, +0.578] | +0.050 [-0.111, +0.235] | +0.160 [-0.050, +0.351] |
| COLO829T_fiberseq | SNV | d9 | satellites | 11 | +0.964 [+0.000, +0.990] | -0.273 [-0.625, +0.000] | +0.241 [-0.300, +0.496] |
| COLO829T_fiberseq | SNV | d9 | L1H > 500 bp | 3 | -0.100 [-0.500, +0.000] | +0.000 [+0.000, +0.000] | -0.071 [-0.333, +0.000] |
| COLO829T_fiberseq | SNV | d9 | GC < 25% or > 65% | 157 | -0.069 [-0.265, +0.153] | +0.096 [+0.044, +0.157] | +0.016 [-0.112, +0.132] |
| COLO829T_fiberseq | SNV | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | SNV | d46 | not in any difficult region | 2000 | -0.018 [-0.032, -0.005] | +0.091 [+0.070, +0.109] | +0.052 [+0.037, +0.065] |
| COLO829T_fiberseq | SNV | d46 | all difficult regions | 443 | +0.356 [+0.169, +0.468] | +0.050 [+0.007, +0.090] | +0.211 [+0.097, +0.299] |
| COLO829T_fiberseq | SNV | d46 | low mappability or segdup | 176 | +0.408 [+0.219, +0.578] | +0.040 [-0.038, +0.123] | +0.295 [+0.144, +0.398] |
| COLO829T_fiberseq | SNV | d46 | low mappability | 144 | +0.414 [+0.214, +0.609] | +0.069 [-0.013, +0.156] | +0.324 [+0.161, +0.436] |
| COLO829T_fiberseq | SNV | d46 | segmental duplications | 95 | +0.291 [+0.106, +0.584] | -0.042 [-0.155, +0.076] | +0.214 [+0.064, +0.343] |
| COLO829T_fiberseq | SNV | d46 | tandem repeats and homopolymers | 166 | +0.292 [+0.054, +0.513] | +0.036 [-0.026, +0.094] | +0.196 [+0.041, +0.322] |
| COLO829T_fiberseq | SNV | d46 | homopolymers >= 7 bp | 70 | +0.002 [-0.139, +0.149] | +0.043 [-0.070, +0.152] | +0.029 [-0.076, +0.131] |
| COLO829T_fiberseq | SNV | d46 | tandem repeats <= 50 bp | 20 | +0.341 [+0.013, +0.639] | +0.000 [-0.200, +0.208] | +0.134 [-0.090, +0.347] |
| COLO829T_fiberseq | SNV | d46 | satellites | 11 | +0.679 [+0.000, +0.961] | +0.000 [+0.000, +0.000] | +0.489 [+0.000, +0.632] |
| COLO829T_fiberseq | SNV | d46 | L1H > 500 bp | 3 | -0.100 [-0.500, +0.000] | +0.000 [+0.000, +0.000] | -0.071 [-0.333, +0.000] |
| COLO829T_fiberseq | SNV | d46 | GC < 25% or > 65% | 157 | +0.177 [+0.050, +0.305] | +0.076 [+0.024, +0.135] | +0.117 [+0.051, +0.189] |
| COLO829T_fiberseq | SNV | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d9 | not in any difficult region | 35 | +0.246 [+0.163, +0.337] | +0.114 [-0.040, +0.280] | +0.280 [+0.186, +0.375] |
| HG008T_Illumina | INDEL | d9 | all difficult regions | 552 | +0.147 [+0.117, +0.178] | -0.029 [-0.049, -0.009] | +0.085 [+0.066, +0.106] |
| HG008T_Illumina | INDEL | d9 | low mappability or segdup | 41 | +0.022 [-0.106, +0.203] | -0.146 [-0.250, -0.054] | -0.064 [-0.144, +0.030] |
| HG008T_Illumina | INDEL | d9 | low mappability | 28 | +0.120 [-0.093, +0.435] | -0.107 [-0.231, +0.000] | -0.022 [-0.116, +0.094] |
| HG008T_Illumina | INDEL | d9 | segmental duplications | 30 | -0.009 [-0.137, +0.181] | -0.133 [-0.250, -0.034] | -0.076 [-0.153, +0.033] |
| HG008T_Illumina | INDEL | d9 | tandem repeats and homopolymers | 547 | +0.147 [+0.117, +0.179] | -0.029 [-0.049, -0.009] | +0.085 [+0.064, +0.106] |
| HG008T_Illumina | INDEL | d9 | homopolymers >= 7 bp | 465 | +0.136 [+0.106, +0.170] | -0.032 [-0.053, -0.013] | +0.080 [+0.058, +0.102] |
| HG008T_Illumina | INDEL | d9 | tandem repeats <= 50 bp | 36 | +0.183 [+0.000, +0.400] | +0.056 [+0.000, +0.143] | +0.097 [+0.019, +0.212] |
| HG008T_Illumina | INDEL | d9 | satellites | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d9 | L1H > 500 bp | 1 | - [-, -] | +0.000 [+0.000, +0.000] | - [-, -] |
| HG008T_Illumina | INDEL | d9 | GC < 25% or > 65% | 50 | +0.077 [+0.007, +0.144] | -0.080 [-0.167, -0.018] | +0.036 [-0.033, +0.091] |
| HG008T_Illumina | INDEL | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d46 | not in any difficult region | 35 | +0.179 [+0.109, +0.256] | +0.143 [+0.000, +0.296] | +0.222 [+0.139, +0.305] |
| HG008T_Illumina | INDEL | d46 | all difficult regions | 552 | +0.157 [+0.126, +0.193] | -0.029 [-0.050, -0.009] | +0.090 [+0.068, +0.112] |
| HG008T_Illumina | INDEL | d46 | low mappability or segdup | 41 | +0.003 [-0.107, +0.135] | -0.146 [-0.250, -0.054] | -0.068 [-0.146, +0.021] |
| HG008T_Illumina | INDEL | d46 | low mappability | 28 | +0.080 [-0.097, +0.299] | -0.107 [-0.231, +0.000] | -0.027 [-0.118, +0.085] |
| HG008T_Illumina | INDEL | d46 | segmental duplications | 30 | -0.017 [-0.137, +0.146] | -0.133 [-0.250, -0.034] | -0.077 [-0.154, +0.029] |
| HG008T_Illumina | INDEL | d46 | tandem repeats and homopolymers | 547 | +0.159 [+0.126, +0.195] | -0.031 [-0.051, -0.011] | +0.089 [+0.067, +0.111] |
| HG008T_Illumina | INDEL | d46 | homopolymers >= 7 bp | 465 | +0.149 [+0.115, +0.186] | -0.037 [-0.058, -0.017] | +0.083 [+0.060, +0.107] |
| HG008T_Illumina | INDEL | d46 | tandem repeats <= 50 bp | 36 | +0.218 [+0.040, +0.437] | +0.083 [+0.000, +0.186] | +0.130 [+0.026, +0.255] |
| HG008T_Illumina | INDEL | d46 | satellites | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d46 | L1H > 500 bp | 1 | - [-, -] | +0.000 [+0.000, +0.000] | - [-, -] |
| HG008T_Illumina | INDEL | d46 | GC < 25% or > 65% | 50 | +0.086 [+0.010, +0.163] | -0.060 [-0.159, +0.022] | +0.051 [-0.023, +0.116] |
| HG008T_Illumina | INDEL | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | SNV | d9 | not in any difficult region | 517 | +0.001 [-0.014, +0.016] | +0.056 [+0.026, +0.087] | +0.015 [-0.003, +0.031] |
| HG008T_Illumina | SNV | d9 | all difficult regions | 180 | +0.196 [+0.121, +0.262] | +0.011 [-0.060, +0.081] | +0.197 [+0.115, +0.268] |
| HG008T_Illumina | SNV | d9 | low mappability or segdup | 67 | +0.293 [+0.162, +0.394] | -0.149 [-0.237, -0.070] | +0.231 [+0.091, +0.335] |
| HG008T_Illumina | SNV | d9 | low mappability | 50 | +0.320 [+0.183, +0.430] | -0.180 [-0.281, -0.079] | +0.254 [+0.109, +0.365] |
| HG008T_Illumina | SNV | d9 | segmental duplications | 48 | +0.297 [+0.123, +0.438] | -0.188 [-0.311, -0.081] | +0.173 [+0.016, +0.308] |
| HG008T_Illumina | SNV | d9 | tandem repeats and homopolymers | 87 | +0.134 [+0.067, +0.195] | +0.011 [-0.088, +0.122] | +0.152 [+0.071, +0.223] |
| HG008T_Illumina | SNV | d9 | homopolymers >= 7 bp | 42 | +0.070 [+0.015, +0.127] | -0.024 [-0.105, +0.054] | +0.068 [+0.003, +0.129] |
| HG008T_Illumina | SNV | d9 | tandem repeats <= 50 bp | 14 | -0.030 [-0.205, +0.142] | +0.286 [+0.000, +0.583] | +0.024 [-0.188, +0.250] |
| HG008T_Illumina | SNV | d9 | satellites | 2 | - [-, -] | -0.500 [-1.000, +0.000] | - [-, -] |
| HG008T_Illumina | SNV | d9 | L1H > 500 bp | 1 | +0.167 [+0.000, +0.667] | +0.000 [+0.000, +0.000] | +0.167 [+0.000, +0.500] |
| HG008T_Illumina | SNV | d9 | GC < 25% or > 65% | 59 | +0.197 [+0.112, +0.289] | +0.136 [+0.016, +0.255] | +0.195 [+0.108, +0.286] |
| HG008T_Illumina | SNV | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | SNV | d46 | not in any difficult region | 517 | -0.013 [-0.029, +0.002] | +0.060 [+0.033, +0.090] | +0.003 [-0.014, +0.020] |
| HG008T_Illumina | SNV | d46 | all difficult regions | 180 | +0.167 [+0.100, +0.217] | +0.083 [+0.022, +0.144] | +0.190 [+0.114, +0.249] |
| HG008T_Illumina | SNV | d46 | low mappability or segdup | 67 | +0.232 [+0.129, +0.315] | +0.000 [-0.070, +0.071] | +0.248 [+0.123, +0.337] |
| HG008T_Illumina | SNV | d46 | low mappability | 50 | +0.246 [+0.140, +0.334] | +0.020 [-0.067, +0.102] | +0.281 [+0.148, +0.372] |
| HG008T_Illumina | SNV | d46 | segmental duplications | 48 | +0.222 [+0.118, +0.336] | -0.042 [-0.128, +0.034] | +0.210 [+0.088, +0.321] |
| HG008T_Illumina | SNV | d46 | tandem repeats and homopolymers | 87 | +0.113 [+0.050, +0.168] | +0.034 [-0.067, +0.142] | +0.137 [+0.056, +0.208] |
| HG008T_Illumina | SNV | d46 | homopolymers >= 7 bp | 42 | +0.062 [+0.002, +0.130] | -0.024 [-0.105, +0.054] | +0.060 [-0.009, +0.130] |
| HG008T_Illumina | SNV | d46 | tandem repeats <= 50 bp | 14 | -0.030 [-0.188, +0.121] | +0.286 [+0.000, +0.583] | +0.024 [-0.172, +0.233] |
| HG008T_Illumina | SNV | d46 | satellites | 2 | -0.006 [-0.333, +0.000] | -0.500 [-1.000, +0.000] | - [-, -] |
| HG008T_Illumina | SNV | d46 | L1H > 500 bp | 1 | +0.167 [+0.000, +0.667] | +0.000 [+0.000, +0.000] | +0.167 [+0.000, +0.500] |
| HG008T_Illumina | SNV | d46 | GC < 25% or > 65% | 59 | +0.163 [+0.088, +0.246] | +0.186 [+0.077, +0.298] | +0.184 [+0.108, +0.265] |
| HG008T_Illumina | SNV | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d9 | not in any difficult region | 35 | +0.089 [+0.059, +0.122] | +0.086 [+0.000, +0.188] | +0.140 [+0.098, +0.185] |
| HG008T_ONT | INDEL | d9 | all difficult regions | 552 | -0.004 [-0.005, -0.003] | -0.013 [-0.032, +0.006] | -0.007 [-0.009, -0.005] |
| HG008T_ONT | INDEL | d9 | low mappability or segdup | 41 | -0.002 [-0.005, +0.000] | -0.049 [-0.125, +0.000] | -0.004 [-0.009, +0.001] |
| HG008T_ONT | INDEL | d9 | low mappability | 28 | -0.002 [-0.006, +0.001] | -0.071 [-0.174, +0.000] | -0.004 [-0.012, +0.002] |
| HG008T_ONT | INDEL | d9 | segmental duplications | 30 | -0.001 [-0.004, +0.001] | -0.033 [-0.115, +0.000] | -0.002 [-0.008, +0.003] |
| HG008T_ONT | INDEL | d9 | tandem repeats and homopolymers | 547 | -0.004 [-0.005, -0.003] | -0.011 [-0.030, +0.008] | -0.007 [-0.009, -0.005] |
| HG008T_ONT | INDEL | d9 | homopolymers >= 7 bp | 465 | -0.004 [-0.005, -0.003] | +0.000 [-0.018, +0.018] | -0.007 [-0.009, -0.005] |
| HG008T_ONT | INDEL | d9 | tandem repeats <= 50 bp | 36 | -0.005 [-0.025, +0.012] | -0.028 [-0.128, +0.063] | -0.009 [-0.041, +0.020] |
| HG008T_ONT | INDEL | d9 | satellites | 0 | +0.000 [+0.000, +0.000] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d9 | L1H > 500 bp | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| HG008T_ONT | INDEL | d9 | GC < 25% or > 65% | 50 | -0.004 [-0.008, -0.001] | -0.060 [-0.132, +0.000] | -0.008 [-0.015, -0.002] |
| HG008T_ONT | INDEL | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d46 | not in any difficult region | 35 | +0.054 [+0.034, +0.077] | +0.086 [+0.000, +0.188] | +0.089 [+0.057, +0.121] |
| HG008T_ONT | INDEL | d46 | all difficult regions | 552 | -0.004 [-0.005, -0.003] | -0.004 [-0.023, +0.016] | -0.007 [-0.009, -0.005] |
| HG008T_ONT | INDEL | d46 | low mappability or segdup | 41 | -0.001 [-0.005, +0.002] | +0.000 [-0.109, +0.086] | -0.002 [-0.010, +0.004] |
| HG008T_ONT | INDEL | d46 | low mappability | 28 | -0.001 [-0.005, +0.002] | +0.000 [-0.097, +0.103] | -0.002 [-0.010, +0.004] |
| HG008T_ONT | INDEL | d46 | segmental duplications | 30 | -0.000 [-0.005, +0.003] | +0.033 [-0.091, +0.129] | -0.001 [-0.009, +0.006] |
| HG008T_ONT | INDEL | d46 | tandem repeats and homopolymers | 547 | -0.004 [-0.005, -0.003] | -0.004 [-0.024, +0.016] | -0.007 [-0.009, -0.005] |
| HG008T_ONT | INDEL | d46 | homopolymers >= 7 bp | 465 | -0.004 [-0.005, -0.003] | +0.009 [-0.012, +0.028] | -0.007 [-0.009, -0.005] |
| HG008T_ONT | INDEL | d46 | tandem repeats <= 50 bp | 36 | -0.008 [-0.027, +0.009] | -0.028 [-0.129, +0.065] | -0.013 [-0.043, +0.015] |
| HG008T_ONT | INDEL | d46 | satellites | 0 | +0.000 [+0.000, +0.000] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d46 | L1H > 500 bp | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| HG008T_ONT | INDEL | d46 | GC < 25% or > 65% | 50 | -0.004 [-0.008, -0.001] | -0.060 [-0.132, +0.000] | -0.008 [-0.016, -0.003] |
| HG008T_ONT | INDEL | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | SNV | d9 | not in any difficult region | 517 | +0.000 [-0.012, +0.013] | -0.004 [-0.028, +0.019] | -0.001 [-0.015, +0.014] |
| HG008T_ONT | SNV | d9 | all difficult regions | 180 | +0.183 [+0.126, +0.247] | -0.050 [-0.102, +0.000] | +0.148 [+0.092, +0.209] |
| HG008T_ONT | SNV | d9 | low mappability or segdup | 67 | +0.209 [+0.118, +0.307] | -0.090 [-0.190, +0.000] | +0.178 [+0.087, +0.275] |
| HG008T_ONT | SNV | d9 | low mappability | 50 | +0.226 [+0.123, +0.336] | -0.060 [-0.158, +0.048] | +0.219 [+0.113, +0.323] |
| HG008T_ONT | SNV | d9 | segmental duplications | 48 | +0.241 [+0.102, +0.412] | -0.146 [-0.275, -0.024] | +0.157 [+0.030, +0.290] |
| HG008T_ONT | SNV | d9 | tandem repeats and homopolymers | 87 | +0.196 [+0.118, +0.281] | -0.046 [-0.110, +0.012] | +0.156 [+0.079, +0.227] |
| HG008T_ONT | SNV | d9 | homopolymers >= 7 bp | 42 | +0.140 [+0.063, +0.238] | +0.000 [-0.065, +0.071] | +0.116 [+0.048, +0.187] |
| HG008T_ONT | SNV | d9 | tandem repeats <= 50 bp | 14 | +0.126 [+0.000, +0.270] | +0.000 [+0.000, +0.000] | +0.115 [+0.000, +0.238] |
| HG008T_ONT | SNV | d9 | satellites | 2 | - [-, -] | -0.500 [-1.000, +0.000] | - [-, -] |
| HG008T_ONT | SNV | d9 | L1H > 500 bp | 1 | +0.167 [+0.000, +0.750] | +0.000 [+0.000, +0.000] | +0.214 [+0.000, +0.667] |
| HG008T_ONT | SNV | d9 | GC < 25% or > 65% | 59 | +0.138 [+0.066, +0.218] | -0.017 [-0.119, +0.083] | +0.101 [+0.031, +0.173] |
| HG008T_ONT | SNV | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | SNV | d46 | not in any difficult region | 517 | -0.003 [-0.018, +0.011] | +0.014 [-0.011, +0.039] | +0.000 [-0.016, +0.016] |
| HG008T_ONT | SNV | d46 | all difficult regions | 180 | +0.086 [+0.031, +0.142] | +0.028 [-0.019, +0.076] | +0.091 [+0.035, +0.150] |
| HG008T_ONT | SNV | d46 | low mappability or segdup | 67 | +0.122 [+0.040, +0.216] | +0.060 [-0.028, +0.154] | +0.144 [+0.049, +0.248] |
| HG008T_ONT | SNV | d46 | low mappability | 50 | +0.127 [+0.040, +0.237] | +0.080 [-0.035, +0.205] | +0.161 [+0.053, +0.279] |
| HG008T_ONT | SNV | d46 | segmental duplications | 48 | +0.114 [+0.016, +0.252] | +0.042 [-0.079, +0.159] | +0.130 [+0.019, +0.265] |
| HG008T_ONT | SNV | d46 | tandem repeats and homopolymers | 87 | +0.090 [+0.027, +0.147] | +0.000 [-0.044, +0.045] | +0.091 [+0.022, +0.154] |
| HG008T_ONT | SNV | d46 | homopolymers >= 7 bp | 42 | +0.060 [-0.002, +0.126] | +0.000 [-0.093, +0.091] | +0.056 [-0.014, +0.127] |
| HG008T_ONT | SNV | d46 | tandem repeats <= 50 bp | 14 | -0.009 [-0.106, +0.097] | +0.000 [+0.000, +0.000] | -0.009 [-0.106, +0.093] |
| HG008T_ONT | SNV | d46 | satellites | 2 | +0.186 [+0.000, +0.311] | +0.000 [+0.000, +0.000] | +0.259 [+0.000, +0.357] |
| HG008T_ONT | SNV | d46 | L1H > 500 bp | 1 | -0.167 [-0.500, +0.000] | -1.000 [-1.000, -1.000] | - [-, -] |
| HG008T_ONT | SNV | d46 | GC < 25% or > 65% | 59 | +0.054 [-0.007, +0.122] | +0.034 [-0.060, +0.131] | +0.053 [-0.009, +0.123] |
| HG008T_ONT | SNV | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d9 | not in any difficult region | 35 | +0.167 [+0.106, +0.236] | +0.486 [+0.308, +0.656] | +0.250 [+0.166, +0.338] |
| HG008T_PacBio | INDEL | d9 | all difficult regions | 552 | -0.042 [-0.053, -0.032] | -0.031 [-0.052, -0.012] | -0.053 [-0.067, -0.040] |
| HG008T_PacBio | INDEL | d9 | low mappability or segdup | 41 | -0.057 [-0.096, -0.018] | -0.024 [-0.111, +0.058] | -0.078 [-0.126, -0.023] |
| HG008T_PacBio | INDEL | d9 | low mappability | 28 | -0.058 [-0.107, -0.010] | +0.000 [-0.097, +0.103] | -0.080 [-0.140, -0.011] |
| HG008T_PacBio | INDEL | d9 | segmental duplications | 30 | -0.061 [-0.112, -0.013] | +0.000 [-0.095, +0.091] | -0.076 [-0.132, -0.018] |
| HG008T_PacBio | INDEL | d9 | tandem repeats and homopolymers | 547 | -0.042 [-0.052, -0.032] | -0.033 [-0.053, -0.014] | -0.053 [-0.067, -0.041] |
| HG008T_PacBio | INDEL | d9 | homopolymers >= 7 bp | 465 | -0.052 [-0.064, -0.041] | -0.030 [-0.052, -0.011] | -0.065 [-0.078, -0.052] |
| HG008T_PacBio | INDEL | d9 | tandem repeats <= 50 bp | 36 | +0.001 [-0.034, +0.034] | +0.000 [-0.083, +0.081] | +0.001 [-0.047, +0.046] |
| HG008T_PacBio | INDEL | d9 | satellites | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d9 | L1H > 500 bp | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| HG008T_PacBio | INDEL | d9 | GC < 25% or > 65% | 50 | -0.023 [-0.047, -0.002] | -0.060 [-0.153, +0.021] | -0.035 [-0.071, -0.003] |
| HG008T_PacBio | INDEL | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d46 | not in any difficult region | 35 | +0.143 [+0.083, +0.204] | +0.486 [+0.308, +0.656] | +0.221 [+0.136, +0.303] |
| HG008T_PacBio | INDEL | d46 | all difficult regions | 552 | -0.037 [-0.046, -0.029] | -0.031 [-0.050, -0.012] | -0.048 [-0.058, -0.037] |
| HG008T_PacBio | INDEL | d46 | low mappability or segdup | 41 | -0.045 [-0.083, -0.017] | -0.049 [-0.133, +0.000] | -0.060 [-0.107, -0.023] |
| HG008T_PacBio | INDEL | d46 | low mappability | 28 | -0.047 [-0.091, -0.012] | -0.036 [-0.115, +0.000] | -0.062 [-0.113, -0.017] |
| HG008T_PacBio | INDEL | d46 | segmental duplications | 30 | -0.051 [-0.105, -0.013] | -0.033 [-0.121, +0.000] | -0.063 [-0.119, -0.017] |
| HG008T_PacBio | INDEL | d46 | tandem repeats and homopolymers | 547 | -0.038 [-0.046, -0.029] | -0.033 [-0.052, -0.015] | -0.048 [-0.059, -0.038] |
| HG008T_PacBio | INDEL | d46 | homopolymers >= 7 bp | 465 | -0.046 [-0.057, -0.036] | -0.030 [-0.051, -0.012] | -0.056 [-0.068, -0.044] |
| HG008T_PacBio | INDEL | d46 | tandem repeats <= 50 bp | 36 | -0.012 [-0.048, +0.016] | +0.000 [-0.083, +0.081] | -0.013 [-0.061, +0.027] |
| HG008T_PacBio | INDEL | d46 | satellites | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d46 | L1H > 500 bp | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| HG008T_PacBio | INDEL | d46 | GC < 25% or > 65% | 50 | -0.027 [-0.050, -0.009] | -0.060 [-0.153, +0.021] | -0.040 [-0.073, -0.011] |
| HG008T_PacBio | INDEL | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | SNV | d9 | not in any difficult region | 517 | +0.015 [-0.000, +0.028] | +0.015 [-0.008, +0.039] | +0.016 [+0.000, +0.031] |
| HG008T_PacBio | SNV | d9 | all difficult regions | 180 | +0.274 [+0.156, +0.367] | -0.039 [-0.102, +0.027] | +0.228 [+0.103, +0.326] |
| HG008T_PacBio | SNV | d9 | low mappability or segdup | 67 | +0.289 [+0.143, +0.394] | -0.104 [-0.222, +0.000] | +0.260 [+0.087, +0.383] |
| HG008T_PacBio | SNV | d9 | low mappability | 50 | +0.300 [+0.151, +0.412] | -0.080 [-0.194, +0.044] | +0.303 [+0.123, +0.423] |
| HG008T_PacBio | SNV | d9 | segmental duplications | 48 | +0.378 [+0.191, +0.532] | -0.146 [-0.306, +0.000] | +0.271 [+0.053, +0.422] |
| HG008T_PacBio | SNV | d9 | tandem repeats and homopolymers | 87 | +0.309 [+0.159, +0.415] | -0.069 [-0.147, +0.009] | +0.238 [+0.082, +0.343] |
| HG008T_PacBio | SNV | d9 | homopolymers >= 7 bp | 42 | +0.196 [+0.098, +0.300] | -0.024 [-0.106, +0.059] | +0.118 [+0.031, +0.198] |
| HG008T_PacBio | SNV | d9 | tandem repeats <= 50 bp | 14 | +0.203 [+0.071, +0.347] | -0.071 [-0.250, +0.000] | +0.157 [+0.027, +0.272] |
| HG008T_PacBio | SNV | d9 | satellites | 2 | - [-, -] | -0.500 [-1.000, +0.000] | - [-, -] |
| HG008T_PacBio | SNV | d9 | L1H > 500 bp | 1 | +0.000 [+0.000, +0.000] | +0.000 [+0.000, +0.000] | - [-, -] |
| HG008T_PacBio | SNV | d9 | GC < 25% or > 65% | 59 | +0.178 [+0.095, +0.272] | +0.085 [-0.020, +0.193] | +0.158 [+0.075, +0.240] |
| HG008T_PacBio | SNV | d9 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | SNV | d46 | not in any difficult region | 517 | -0.005 [-0.020, +0.010] | +0.019 [-0.005, +0.044] | +0.001 [-0.015, +0.015] |
| HG008T_PacBio | SNV | d46 | all difficult regions | 180 | +0.232 [+0.109, +0.323] | -0.022 [-0.094, +0.049] | +0.206 [+0.079, +0.311] |
| HG008T_PacBio | SNV | d46 | low mappability or segdup | 67 | +0.242 [+0.075, +0.350] | -0.179 [-0.319, -0.043] | +0.201 [+0.010, +0.325] |
| HG008T_PacBio | SNV | d46 | low mappability | 50 | +0.270 [+0.102, +0.387] | -0.160 [-0.298, +0.000] | +0.255 [+0.054, +0.383] |
| HG008T_PacBio | SNV | d46 | segmental duplications | 48 | +0.270 [+0.046, +0.422] | -0.271 [-0.471, -0.094] | +0.149 [-0.094, +0.305] |
| HG008T_PacBio | SNV | d46 | tandem repeats and homopolymers | 87 | +0.239 [+0.103, +0.335] | -0.023 [-0.081, +0.033] | +0.219 [+0.072, +0.318] |
| HG008T_PacBio | SNV | d46 | homopolymers >= 7 bp | 42 | +0.125 [+0.036, +0.226] | +0.000 [-0.067, +0.071] | +0.092 [+0.019, +0.169] |
| HG008T_PacBio | SNV | d46 | tandem repeats <= 50 bp | 14 | +0.143 [+0.034, +0.274] | +0.000 [+0.000, +0.000] | +0.133 [+0.032, +0.232] |
| HG008T_PacBio | SNV | d46 | satellites | 2 | +0.495 [+0.000, +0.997] | +0.000 [+0.000, +0.000] | +0.490 [+0.321, +0.778] |
| HG008T_PacBio | SNV | d46 | L1H > 500 bp | 1 | +0.333 [+0.000, +0.750] | +1.000 [+1.000, +1.000] | - [-, -] |
| HG008T_PacBio | SNV | d46 | GC < 25% or > 65% | 59 | +0.150 [+0.071, +0.234] | +0.153 [+0.022, +0.280] | +0.160 [+0.070, +0.253] |
| HG008T_PacBio | SNV | d46 | CMRG falsely duplicated genes | 0 | - [-, -] | - [-, -] | - [-, -] |

## 3. Per call set

| set | kind | arm | stratum | PASS calls | TP / FP / FN | P [95% CI] | R [95% CI] | F1 [95% CI] |
|---|---|---|---|---:|---|---|---|---|
| COLO829T_Illumina | SNV | linear | not in any difficult region | 1,298 | 1186 / 112 / 814 | 0.914 [0.895, 0.931] | 0.593 [0.567, 0.620] | 0.719 [0.697, 0.740] |
| COLO829T_Illumina | SNV | linear | all difficult regions | 2,346 | 237 / 2109 / 206 | 0.101 [0.048, 0.452] | 0.535 [0.494, 0.580] | 0.170 [0.088, 0.492] |
| COLO829T_Illumina | SNV | linear | low mappability or segdup | 2,146 | 91 / 2055 / 85 | 0.042 [0.019, 0.263] | 0.517 [0.448, 0.593] | 0.078 [0.036, 0.350] |
| COLO829T_Illumina | SNV | linear | low mappability | 2,109 | 75 / 2034 / 69 | 0.036 [0.015, 0.250] | 0.521 [0.446, 0.605] | 0.067 [0.030, 0.341] |
| COLO829T_Illumina | SNV | linear | segmental duplications | 330 | 45 / 285 / 50 | 0.136 [0.058, 0.420] | 0.474 [0.374, 0.574] | 0.212 [0.102, 0.457] |
| COLO829T_Illumina | SNV | linear | tandem repeats and homopolymers | 2,115 | 83 / 2032 / 83 | 0.039 [0.017, 0.285] | 0.500 [0.423, 0.572] | 0.073 [0.033, 0.360] |
| COLO829T_Illumina | SNV | linear | homopolymers >= 7 bp | 85 | 33 / 52 / 37 | 0.388 [0.208, 0.694] | 0.471 [0.357, 0.585] | 0.426 [0.273, 0.595] |
| COLO829T_Illumina | SNV | linear | tandem repeats <= 50 bp | 13 | 12 / 1 / 8 | 0.923 [0.750, 1.000] | 0.600 [0.391, 0.800] | 0.727 [0.538, 0.867] |
| COLO829T_Illumina | SNV | linear | satellites | 1,856 | 3 / 1853 / 8 | 0.002 [0.000, 0.043] | 0.273 [0.000, 0.556] | 0.003 [0.001, 0.077] |
| COLO829T_Illumina | SNV | linear | L1H > 500 bp | 3 | 0 / 3 / 3 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_Illumina | SNV | linear | GC < 25% or > 65% | 114 | 97 / 17 / 60 | 0.851 [0.779, 0.908] | 0.618 [0.541, 0.697] | 0.716 [0.650, 0.775] |
| COLO829T_Illumina | SNV | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | SNV | d9 | not in any difficult region | 1,593 | 1440 / 153 / 560 | 0.904 [0.887, 0.920] | 0.720 [0.700, 0.742] | 0.802 [0.786, 0.817] |
| COLO829T_Illumina | SNV | d9 | all difficult regions | 1,927 | 265 / 1662 / 178 | 0.138 [0.071, 0.388] | 0.598 [0.554, 0.642] | 0.224 [0.127, 0.472] |
| COLO829T_Illumina | SNV | d9 | low mappability or segdup | 1,572 | 97 / 1475 / 79 | 0.062 [0.029, 0.242] | 0.551 [0.480, 0.630] | 0.111 [0.054, 0.345] |
| COLO829T_Illumina | SNV | d9 | low mappability | 1,500 | 83 / 1417 / 61 | 0.055 [0.025, 0.226] | 0.576 [0.497, 0.669] | 0.101 [0.049, 0.334] |
| COLO829T_Illumina | SNV | d9 | segmental duplications | 458 | 41 / 417 / 54 | 0.090 [0.031, 0.494] | 0.432 [0.346, 0.528] | 0.148 [0.058, 0.462] |
| COLO829T_Illumina | SNV | d9 | tandem repeats and homopolymers | 1,690 | 96 / 1594 / 70 | 0.057 [0.027, 0.204] | 0.578 [0.497, 0.654] | 0.103 [0.052, 0.306] |
| COLO829T_Illumina | SNV | d9 | homopolymers >= 7 bp | 82 | 43 / 39 / 27 | 0.524 [0.394, 0.643] | 0.614 [0.493, 0.725] | 0.566 [0.455, 0.657] |
| COLO829T_Illumina | SNV | d9 | tandem repeats <= 50 bp | 14 | 11 / 3 / 9 | 0.786 [0.545, 1.000] | 0.550 [0.312, 0.778] | 0.647 [0.417, 0.824] |
| COLO829T_Illumina | SNV | d9 | satellites | 1,118 | 3 / 1115 / 8 | 0.003 [0.000, 0.083] | 0.273 [0.000, 0.556] | 0.005 [0.001, 0.133] |
| COLO829T_Illumina | SNV | d9 | L1H > 500 bp | 3 | 1 / 2 / 2 | 0.333 [0.000, 1.000] | 0.333 [0.000, 1.000] | 0.333 [0.200, 0.857] |
| COLO829T_Illumina | SNV | d9 | GC < 25% or > 65% | 133 | 105 / 28 / 52 | 0.789 [0.687, 0.882] | 0.669 [0.597, 0.736] | 0.724 [0.655, 0.784] |
| COLO829T_Illumina | SNV | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | SNV | d46 | not in any difficult region | 1,684 | 1493 / 191 / 507 | 0.887 [0.868, 0.904] | 0.747 [0.726, 0.768] | 0.811 [0.794, 0.827] |
| COLO829T_Illumina | SNV | d46 | all difficult regions | 329 | 252 / 77 / 191 | 0.766 [0.696, 0.830] | 0.569 [0.523, 0.616] | 0.653 [0.605, 0.697] |
| COLO829T_Illumina | SNV | d46 | low mappability or segdup | 134 | 94 / 40 / 82 | 0.701 [0.586, 0.811] | 0.534 [0.455, 0.617] | 0.606 [0.528, 0.681] |
| COLO829T_Illumina | SNV | d46 | low mappability | 110 | 78 / 32 / 66 | 0.709 [0.577, 0.830] | 0.542 [0.453, 0.634] | 0.614 [0.524, 0.697] |
| COLO829T_Illumina | SNV | d46 | segmental duplications | 67 | 44 / 23 / 51 | 0.657 [0.479, 0.827] | 0.463 [0.352, 0.576] | 0.543 [0.424, 0.652] |
| COLO829T_Illumina | SNV | d46 | tandem repeats and homopolymers | 124 | 85 / 39 / 81 | 0.685 [0.569, 0.807] | 0.512 [0.436, 0.588] | 0.586 [0.508, 0.662] |
| COLO829T_Illumina | SNV | d46 | homopolymers >= 7 bp | 53 | 36 / 17 / 34 | 0.679 [0.553, 0.793] | 0.514 [0.392, 0.629] | 0.585 [0.471, 0.678] |
| COLO829T_Illumina | SNV | d46 | tandem repeats <= 50 bp | 13 | 12 / 1 / 8 | 0.923 [0.750, 1.000] | 0.600 [0.385, 0.800] | 0.727 [0.526, 0.875] |
| COLO829T_Illumina | SNV | d46 | satellites | 7 | 3 / 4 / 8 | 0.429 [0.250, 1.000] | 0.273 [0.000, 0.556] | 0.333 [0.149, 0.545] |
| COLO829T_Illumina | SNV | d46 | L1H > 500 bp | 5 | 1 / 4 / 2 | 0.200 [0.000, 1.000] | 0.333 [0.000, 1.000] | 0.250 [0.143, 0.750] |
| COLO829T_Illumina | SNV | d46 | GC < 25% or > 65% | 118 | 101 / 17 / 56 | 0.856 [0.780, 0.922] | 0.643 [0.575, 0.713] | 0.735 [0.676, 0.791] |
| COLO829T_Illumina | SNV | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | SNV | linear | not in any difficult region | 1,639 | 1100 / 539 / 900 | 0.671 [0.641, 0.700] | 0.550 [0.522, 0.580] | 0.605 [0.582, 0.626] |
| COLO829T_ONT | SNV | linear | all difficult regions | 2,261 | 225 / 2036 / 218 | 0.100 [0.053, 0.238] | 0.508 [0.462, 0.561] | 0.166 [0.096, 0.325] |
| COLO829T_ONT | SNV | linear | low mappability or segdup | 2,001 | 82 / 1919 / 94 | 0.041 [0.020, 0.116] | 0.466 [0.387, 0.555] | 0.075 [0.037, 0.190] |
| COLO829T_ONT | SNV | linear | low mappability | 1,867 | 68 / 1799 / 76 | 0.036 [0.017, 0.114] | 0.472 [0.391, 0.569] | 0.068 [0.033, 0.187] |
| COLO829T_ONT | SNV | linear | segmental duplications | 777 | 41 / 736 / 54 | 0.053 [0.028, 0.104] | 0.432 [0.315, 0.550] | 0.094 [0.052, 0.171] |
| COLO829T_ONT | SNV | linear | tandem repeats and homopolymers | 1,472 | 82 / 1390 / 84 | 0.056 [0.024, 0.271] | 0.494 [0.417, 0.571] | 0.100 [0.046, 0.351] |
| COLO829T_ONT | SNV | linear | homopolymers >= 7 bp | 78 | 30 / 48 / 40 | 0.385 [0.262, 0.516] | 0.429 [0.302, 0.551] | 0.405 [0.288, 0.510] |
| COLO829T_ONT | SNV | linear | tandem repeats <= 50 bp | 32 | 11 / 21 / 9 | 0.344 [0.156, 0.577] | 0.550 [0.304, 0.765] | 0.423 [0.216, 0.620] |
| COLO829T_ONT | SNV | linear | satellites | 1,191 | 5 / 1186 / 6 | 0.004 [0.001, 0.200] | 0.455 [0.125, 0.778] | 0.008 [0.002, 0.305] |
| COLO829T_ONT | SNV | linear | L1H > 500 bp | 5 | 3 / 2 / 0 | 0.600 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.750 [0.444, 1.000] |
| COLO829T_ONT | SNV | linear | GC < 25% or > 65% | 209 | 90 / 119 / 67 | 0.431 [0.294, 0.626] | 0.573 [0.500, 0.649] | 0.492 [0.381, 0.609] |
| COLO829T_ONT | SNV | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | SNV | d9 | not in any difficult region | 1,698 | 1444 / 254 / 556 | 0.850 [0.828, 0.870] | 0.722 [0.702, 0.744] | 0.781 [0.765, 0.797] |
| COLO829T_ONT | SNV | d9 | all difficult regions | 460 | 268 / 192 / 175 | 0.583 [0.507, 0.657] | 0.605 [0.559, 0.651] | 0.594 [0.541, 0.644] |
| COLO829T_ONT | SNV | d9 | low mappability or segdup | 202 | 97 / 105 / 79 | 0.480 [0.365, 0.623] | 0.551 [0.475, 0.634] | 0.513 [0.424, 0.608] |
| COLO829T_ONT | SNV | d9 | low mappability | 176 | 83 / 93 / 61 | 0.472 [0.354, 0.619] | 0.576 [0.493, 0.664] | 0.519 [0.421, 0.621] |
| COLO829T_ONT | SNV | d9 | segmental duplications | 127 | 42 / 85 / 53 | 0.331 [0.213, 0.505] | 0.442 [0.338, 0.556] | 0.378 [0.273, 0.500] |
| COLO829T_ONT | SNV | d9 | tandem repeats and homopolymers | 187 | 93 / 94 / 73 | 0.497 [0.431, 0.566] | 0.560 [0.489, 0.643] | 0.527 [0.467, 0.592] |
| COLO829T_ONT | SNV | d9 | homopolymers >= 7 bp | 96 | 39 / 57 / 31 | 0.406 [0.308, 0.506] | 0.557 [0.443, 0.681] | 0.470 [0.372, 0.563] |
| COLO829T_ONT | SNV | d9 | tandem repeats <= 50 bp | 16 | 13 / 3 / 7 | 0.812 [0.583, 1.000] | 0.650 [0.412, 0.862] | 0.722 [0.500, 0.880] |
| COLO829T_ONT | SNV | d9 | satellites | 4 | 2 / 2 / 9 | 0.500 [0.000, 1.000] | 0.182 [0.000, 0.500] | 0.267 [0.100, 0.667] |
| COLO829T_ONT | SNV | d9 | L1H > 500 bp | 5 | 2 / 3 / 1 | 0.400 [0.000, 1.000] | 0.667 [0.000, 1.000] | 0.500 [0.250, 0.833] |
| COLO829T_ONT | SNV | d9 | GC < 25% or > 65% | 144 | 109 / 35 / 48 | 0.757 [0.675, 0.833] | 0.694 [0.629, 0.763] | 0.724 [0.664, 0.779] |
| COLO829T_ONT | SNV | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | SNV | d46 | not in any difficult region | 1,685 | 1436 / 249 / 564 | 0.852 [0.830, 0.872] | 0.718 [0.696, 0.739] | 0.779 [0.762, 0.796] |
| COLO829T_ONT | SNV | d46 | all difficult regions | 462 | 254 / 208 / 189 | 0.550 [0.408, 0.724] | 0.573 [0.528, 0.619] | 0.561 [0.472, 0.645] |
| COLO829T_ONT | SNV | d46 | low mappability or segdup | 253 | 92 / 161 / 84 | 0.364 [0.234, 0.640] | 0.523 [0.446, 0.602] | 0.429 [0.319, 0.578] |
| COLO829T_ONT | SNV | d46 | low mappability | 222 | 79 / 143 / 65 | 0.356 [0.224, 0.650] | 0.549 [0.461, 0.642] | 0.432 [0.314, 0.596] |
| COLO829T_ONT | SNV | d46 | segmental duplications | 189 | 43 / 146 / 52 | 0.228 [0.130, 0.521] | 0.453 [0.343, 0.565] | 0.303 [0.202, 0.482] |
| COLO829T_ONT | SNV | d46 | tandem repeats and homopolymers | 135 | 88 / 47 / 78 | 0.652 [0.557, 0.746] | 0.530 [0.457, 0.610] | 0.585 [0.511, 0.651] |
| COLO829T_ONT | SNV | d46 | homopolymers >= 7 bp | 63 | 34 / 29 / 36 | 0.540 [0.409, 0.683] | 0.486 [0.360, 0.615] | 0.511 [0.391, 0.620] |
| COLO829T_ONT | SNV | d46 | tandem repeats <= 50 bp | 14 | 13 / 1 / 7 | 0.929 [0.750, 1.000] | 0.650 [0.412, 0.862] | 0.765 [0.556, 0.909] |
| COLO829T_ONT | SNV | d46 | satellites | 3 | 2 / 1 / 9 | 0.667 [0.000, 1.000] | 0.182 [0.000, 0.500] | 0.286 [0.111, 0.615] |
| COLO829T_ONT | SNV | d46 | L1H > 500 bp | 5 | 2 / 3 / 1 | 0.400 [0.000, 1.000] | 0.667 [0.000, 1.000] | 0.500 [0.250, 0.833] |
| COLO829T_ONT | SNV | d46 | GC < 25% or > 65% | 135 | 104 / 31 / 53 | 0.770 [0.681, 0.854] | 0.662 [0.596, 0.731] | 0.712 [0.649, 0.770] |
| COLO829T_ONT | SNV | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | SNV | linear | not in any difficult region | 1,405 | 1283 / 122 / 717 | 0.913 [0.894, 0.931] | 0.641 [0.619, 0.666] | 0.754 [0.735, 0.771] |
| COLO829T_fiberseq | SNV | linear | all difficult regions | 641 | 242 / 399 / 201 | 0.378 [0.256, 0.608] | 0.546 [0.502, 0.594] | 0.446 [0.346, 0.583] |
| COLO829T_fiberseq | SNV | linear | low mappability or segdup | 462 | 91 / 371 / 85 | 0.197 [0.114, 0.409] | 0.517 [0.445, 0.595] | 0.285 [0.185, 0.466] |
| COLO829T_fiberseq | SNV | linear | low mappability | 415 | 75 / 340 / 69 | 0.181 [0.102, 0.401] | 0.521 [0.442, 0.606] | 0.268 [0.169, 0.459] |
| COLO829T_fiberseq | SNV | linear | segmental duplications | 326 | 44 / 282 / 51 | 0.135 [0.063, 0.333] | 0.463 [0.366, 0.573] | 0.209 [0.110, 0.400] |
| COLO829T_fiberseq | SNV | linear | tandem repeats and homopolymers | 296 | 86 / 210 / 80 | 0.291 [0.175, 0.577] | 0.518 [0.439, 0.594] | 0.372 [0.255, 0.550] |
| COLO829T_fiberseq | SNV | linear | homopolymers >= 7 bp | 50 | 36 / 14 / 34 | 0.720 [0.575, 0.845] | 0.514 [0.392, 0.635] | 0.600 [0.484, 0.699] |
| COLO829T_fiberseq | SNV | linear | tandem repeats <= 50 bp | 22 | 9 / 13 / 11 | 0.409 [0.174, 0.750] | 0.450 [0.222, 0.667] | 0.429 [0.205, 0.642] |
| COLO829T_fiberseq | SNV | linear | satellites | 140 | 5 / 135 / 6 | 0.036 [0.006, 1.000] | 0.455 [0.125, 0.778] | 0.066 [0.016, 0.800] |
| COLO829T_fiberseq | SNV | linear | L1H > 500 bp | 4 | 2 / 2 / 1 | 0.500 [0.000, 1.000] | 0.667 [0.000, 1.000] | 0.571 [0.286, 1.000] |
| COLO829T_fiberseq | SNV | linear | GC < 25% or > 65% | 133 | 94 / 39 / 63 | 0.707 [0.575, 0.832] | 0.599 [0.518, 0.678] | 0.648 [0.563, 0.723] |
| COLO829T_fiberseq | SNV | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | SNV | d9 | not in any difficult region | 1,684 | 1494 / 190 / 506 | 0.887 [0.869, 0.905] | 0.747 [0.727, 0.769] | 0.811 [0.795, 0.827] |
| COLO829T_fiberseq | SNV | d9 | all difficult regions | 565 | 271 / 294 / 172 | 0.480 [0.265, 0.760] | 0.612 [0.560, 0.660] | 0.538 [0.365, 0.696] |
| COLO829T_fiberseq | SNV | d9 | low mappability or segdup | 337 | 100 / 237 / 76 | 0.297 [0.123, 0.759] | 0.568 [0.481, 0.658] | 0.390 [0.199, 0.684] |
| COLO829T_fiberseq | SNV | d9 | low mappability | 316 | 86 / 230 / 58 | 0.272 [0.109, 0.762] | 0.597 [0.503, 0.689] | 0.374 [0.183, 0.701] |
| COLO829T_fiberseq | SNV | d9 | segmental duplications | 266 | 41 / 225 / 54 | 0.154 [0.050, 0.662] | 0.432 [0.318, 0.556] | 0.227 [0.088, 0.570] |
| COLO829T_fiberseq | SNV | d9 | tandem repeats and homopolymers | 150 | 93 / 57 / 73 | 0.620 [0.513, 0.717] | 0.560 [0.476, 0.642] | 0.589 [0.515, 0.656] |
| COLO829T_fiberseq | SNV | d9 | homopolymers >= 7 bp | 74 | 45 / 29 / 25 | 0.608 [0.479, 0.730] | 0.643 [0.531, 0.754] | 0.625 [0.516, 0.713] |
| COLO829T_fiberseq | SNV | d9 | tandem repeats <= 50 bp | 14 | 10 / 4 / 10 | 0.714 [0.444, 0.933] | 0.500 [0.250, 0.737] | 0.588 [0.345, 0.776] |
| COLO829T_fiberseq | SNV | d9 | satellites | 2 | 2 / 0 / 9 | 1.000 [1.000, 1.000] | 0.182 [0.000, 0.500] | 0.308 [0.118, 0.667] |
| COLO829T_fiberseq | SNV | d9 | L1H > 500 bp | 5 | 2 / 3 / 1 | 0.400 [0.000, 1.000] | 0.667 [0.000, 1.000] | 0.500 [0.250, 0.833] |
| COLO829T_fiberseq | SNV | d9 | GC < 25% or > 65% | 171 | 109 / 62 / 48 | 0.637 [0.423, 0.846] | 0.694 [0.628, 0.761] | 0.665 [0.515, 0.776] |
| COLO829T_fiberseq | SNV | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | SNV | d46 | not in any difficult region | 1,636 | 1464 / 172 / 536 | 0.895 [0.878, 0.911] | 0.732 [0.711, 0.754] | 0.805 [0.789, 0.822] |
| COLO829T_fiberseq | SNV | d46 | all difficult regions | 360 | 264 / 96 / 179 | 0.733 [0.608, 0.849] | 0.596 [0.547, 0.641] | 0.658 [0.593, 0.715] |
| COLO829T_fiberseq | SNV | d46 | low mappability or segdup | 162 | 98 / 64 / 78 | 0.605 [0.408, 0.848] | 0.557 [0.477, 0.646] | 0.580 [0.462, 0.696] |
| COLO829T_fiberseq | SNV | d46 | low mappability | 143 | 85 / 58 / 59 | 0.594 [0.389, 0.856] | 0.590 [0.507, 0.686] | 0.592 [0.464, 0.723] |
| COLO829T_fiberseq | SNV | d46 | segmental duplications | 94 | 40 / 54 / 55 | 0.426 [0.212, 0.821] | 0.421 [0.308, 0.546] | 0.423 [0.263, 0.592] |
| COLO829T_fiberseq | SNV | d46 | tandem repeats and homopolymers | 158 | 92 / 66 / 74 | 0.582 [0.396, 0.811] | 0.554 [0.479, 0.623] | 0.568 [0.452, 0.675] |
| COLO829T_fiberseq | SNV | d46 | homopolymers >= 7 bp | 54 | 39 / 15 / 31 | 0.722 [0.596, 0.837] | 0.557 [0.446, 0.667] | 0.629 [0.521, 0.716] |
| COLO829T_fiberseq | SNV | d46 | tandem repeats <= 50 bp | 12 | 9 / 3 / 11 | 0.750 [0.470, 1.000] | 0.450 [0.208, 0.688] | 0.563 [0.300, 0.757] |
| COLO829T_fiberseq | SNV | d46 | satellites | 7 | 5 / 2 / 6 | 0.714 [0.333, 1.000] | 0.455 [0.125, 0.778] | 0.556 [0.250, 0.857] |
| COLO829T_fiberseq | SNV | d46 | L1H > 500 bp | 5 | 2 / 3 / 1 | 0.400 [0.000, 1.000] | 0.667 [0.000, 1.000] | 0.500 [0.250, 0.833] |
| COLO829T_fiberseq | SNV | d46 | GC < 25% or > 65% | 120 | 106 / 14 / 51 | 0.883 [0.816, 0.943] | 0.675 [0.608, 0.744] | 0.765 [0.712, 0.815] |
| COLO829T_fiberseq | SNV | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | SNV | linear | not in any difficult region | 923 | 450 / 473 / 67 | 0.488 [0.447, 0.529] | 0.870 [0.839, 0.899] | 0.625 [0.589, 0.661] |
| HG008T_Illumina | SNV | linear | all difficult regions | 714 | 98 / 616 / 82 | 0.137 [0.093, 0.207] | 0.544 [0.461, 0.628] | 0.219 [0.158, 0.304] |
| HG008T_Illumina | SNV | linear | low mappability or segdup | 401 | 36 / 365 / 31 | 0.090 [0.047, 0.195] | 0.537 [0.391, 0.695] | 0.154 [0.085, 0.292] |
| HG008T_Illumina | SNV | linear | low mappability | 369 | 28 / 341 / 22 | 0.076 [0.037, 0.180] | 0.560 [0.405, 0.712] | 0.134 [0.069, 0.280] |
| HG008T_Illumina | SNV | linear | segmental duplications | 257 | 22 / 235 / 26 | 0.086 [0.042, 0.197] | 0.458 [0.295, 0.644] | 0.144 [0.074, 0.288] |
| HG008T_Illumina | SNV | linear | tandem repeats and homopolymers | 488 | 34 / 454 / 53 | 0.070 [0.040, 0.127] | 0.391 [0.284, 0.500] | 0.118 [0.071, 0.196] |
| HG008T_Illumina | SNV | linear | homopolymers >= 7 bp | 142 | 16 / 126 / 26 | 0.113 [0.060, 0.182] | 0.381 [0.238, 0.523] | 0.174 [0.099, 0.258] |
| HG008T_Illumina | SNV | linear | tandem repeats <= 50 bp | 24 | 7 / 17 / 7 | 0.292 [0.125, 0.500] | 0.500 [0.222, 0.857] | 0.368 [0.162, 0.563] |
| HG008T_Illumina | SNV | linear | satellites | 172 | 1 / 171 / 1 | 0.006 [0.000, 0.333] | 0.500 [0.000, 1.000] | 0.011 [0.005, 0.500] |
| HG008T_Illumina | SNV | linear | L1H > 500 bp | 3 | 1 / 2 / 0 | 0.333 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.500 [0.286, 1.000] |
| HG008T_Illumina | SNV | linear | GC < 25% or > 65% | 122 | 37 / 85 / 22 | 0.303 [0.208, 0.402] | 0.627 [0.489, 0.746] | 0.409 [0.298, 0.509] |
| HG008T_Illumina | SNV | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | SNV | d9 | not in any difficult region | 981 | 479 / 502 / 38 | 0.488 [0.447, 0.529] | 0.926 [0.901, 0.950] | 0.640 [0.602, 0.675] |
| HG008T_Illumina | SNV | d9 | all difficult regions | 300 | 100 / 200 / 80 | 0.333 [0.270, 0.399] | 0.556 [0.458, 0.647] | 0.417 [0.346, 0.485] |
| HG008T_Illumina | SNV | d9 | low mappability or segdup | 68 | 26 / 42 / 41 | 0.382 [0.269, 0.491] | 0.388 [0.258, 0.545] | 0.385 [0.269, 0.497] |
| HG008T_Illumina | SNV | d9 | low mappability | 48 | 19 / 29 / 31 | 0.396 [0.261, 0.521] | 0.380 [0.241, 0.538] | 0.388 [0.260, 0.510] |
| HG008T_Illumina | SNV | d9 | segmental duplications | 34 | 13 / 21 / 35 | 0.382 [0.208, 0.536] | 0.271 [0.125, 0.458] | 0.317 [0.162, 0.468] |
| HG008T_Illumina | SNV | d9 | tandem repeats and homopolymers | 172 | 35 / 137 / 52 | 0.203 [0.142, 0.274] | 0.402 [0.290, 0.519] | 0.270 [0.194, 0.348] |
| HG008T_Illumina | SNV | d9 | homopolymers >= 7 bp | 82 | 15 / 67 / 27 | 0.183 [0.104, 0.273] | 0.357 [0.214, 0.500] | 0.242 [0.145, 0.333] |
| HG008T_Illumina | SNV | d9 | tandem repeats <= 50 bp | 42 | 11 / 31 / 3 | 0.262 [0.118, 0.429] | 0.786 [0.500, 1.000] | 0.393 [0.194, 0.575] |
| HG008T_Illumina | SNV | d9 | satellites | 0 | 0 / 0 / 2 | - [-, -] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_Illumina | SNV | d9 | L1H > 500 bp | 2 | 1 / 1 / 0 | 0.500 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.667 [0.400, 1.000] |
| HG008T_Illumina | SNV | d9 | GC < 25% or > 65% | 90 | 45 / 45 / 14 | 0.500 [0.374, 0.624] | 0.763 [0.635, 0.872] | 0.604 [0.483, 0.707] |
| HG008T_Illumina | SNV | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | SNV | d46 | not in any difficult region | 1,014 | 481 / 533 / 36 | 0.474 [0.433, 0.516] | 0.930 [0.908, 0.951] | 0.628 [0.591, 0.664] |
| HG008T_Illumina | SNV | d46 | all difficult regions | 372 | 113 / 259 / 67 | 0.304 [0.251, 0.363] | 0.628 [0.543, 0.706] | 0.409 [0.347, 0.471] |
| HG008T_Illumina | SNV | d46 | low mappability or segdup | 112 | 36 / 76 / 31 | 0.321 [0.232, 0.419] | 0.537 [0.397, 0.686] | 0.402 [0.302, 0.506] |
| HG008T_Illumina | SNV | d46 | low mappability | 90 | 29 / 61 / 21 | 0.322 [0.225, 0.427] | 0.580 [0.439, 0.714] | 0.414 [0.300, 0.527] |
| HG008T_Illumina | SNV | d46 | segmental duplications | 65 | 20 / 45 / 28 | 0.308 [0.203, 0.442] | 0.417 [0.267, 0.596] | 0.354 [0.235, 0.486] |
| HG008T_Illumina | SNV | d46 | tandem repeats and homopolymers | 203 | 37 / 166 / 50 | 0.182 [0.125, 0.250] | 0.425 [0.306, 0.544] | 0.255 [0.181, 0.333] |
| HG008T_Illumina | SNV | d46 | homopolymers >= 7 bp | 86 | 15 / 71 / 27 | 0.174 [0.098, 0.267] | 0.357 [0.214, 0.500] | 0.234 [0.137, 0.331] |
| HG008T_Illumina | SNV | d46 | tandem repeats <= 50 bp | 42 | 11 / 31 / 3 | 0.262 [0.122, 0.405] | 0.786 [0.500, 1.000] | 0.393 [0.196, 0.552] |
| HG008T_Illumina | SNV | d46 | satellites | 3 | 0 / 3 / 2 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_Illumina | SNV | d46 | L1H > 500 bp | 2 | 1 / 1 / 0 | 0.500 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.667 [0.400, 1.000] |
| HG008T_Illumina | SNV | d46 | GC < 25% or > 65% | 103 | 48 / 55 / 11 | 0.466 [0.353, 0.583] | 0.814 [0.696, 0.910] | 0.593 [0.478, 0.696] |
| HG008T_Illumina | SNV | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | SNV | linear | not in any difficult region | 985 | 474 / 511 / 43 | 0.481 [0.439, 0.523] | 0.917 [0.893, 0.940] | 0.631 [0.593, 0.667] |
| HG008T_ONT | SNV | linear | all difficult regions | 506 | 113 / 393 / 67 | 0.223 [0.169, 0.286] | 0.628 [0.533, 0.711] | 0.329 [0.260, 0.401] |
| HG008T_ONT | SNV | linear | low mappability or segdup | 267 | 38 / 229 / 29 | 0.142 [0.092, 0.226] | 0.567 [0.420, 0.713] | 0.228 [0.154, 0.330] |
| HG008T_ONT | SNV | linear | low mappability | 241 | 30 / 211 / 20 | 0.124 [0.074, 0.207] | 0.600 [0.459, 0.738] | 0.206 [0.128, 0.317] |
| HG008T_ONT | SNV | linear | segmental duplications | 186 | 24 / 162 / 24 | 0.129 [0.072, 0.240] | 0.500 [0.327, 0.691] | 0.205 [0.121, 0.336] |
| HG008T_ONT | SNV | linear | tandem repeats and homopolymers | 285 | 40 / 245 / 47 | 0.140 [0.088, 0.212] | 0.460 [0.339, 0.576] | 0.215 [0.143, 0.301] |
| HG008T_ONT | SNV | linear | homopolymers >= 7 bp | 113 | 15 / 98 / 27 | 0.133 [0.076, 0.198] | 0.357 [0.217, 0.489] | 0.194 [0.114, 0.275] |
| HG008T_ONT | SNV | linear | tandem repeats <= 50 bp | 35 | 11 / 24 / 3 | 0.314 [0.150, 0.483] | 0.786 [0.500, 1.000] | 0.449 [0.238, 0.621] |
| HG008T_ONT | SNV | linear | satellites | 72 | 1 / 71 / 1 | 0.014 [0.000, 0.333] | 0.500 [0.000, 1.000] | 0.027 [0.012, 0.500] |
| HG008T_ONT | SNV | linear | L1H > 500 bp | 6 | 1 / 5 / 0 | 0.167 [0.000, 0.600] | 1.000 [1.000, 1.000] | 0.286 [0.182, 0.800] |
| HG008T_ONT | SNV | linear | GC < 25% or > 65% | 113 | 45 / 68 / 14 | 0.398 [0.297, 0.505] | 0.763 [0.625, 0.878] | 0.523 [0.411, 0.621] |
| HG008T_ONT | SNV | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | SNV | d9 | not in any difficult region | 980 | 472 / 508 / 45 | 0.482 [0.439, 0.525] | 0.913 [0.886, 0.938] | 0.631 [0.591, 0.667] |
| HG008T_ONT | SNV | d9 | all difficult regions | 256 | 104 / 152 / 76 | 0.406 [0.336, 0.479] | 0.578 [0.484, 0.661] | 0.477 [0.406, 0.546] |
| HG008T_ONT | SNV | d9 | low mappability or segdup | 91 | 32 / 59 / 35 | 0.352 [0.250, 0.464] | 0.478 [0.333, 0.641] | 0.405 [0.301, 0.520] |
| HG008T_ONT | SNV | d9 | low mappability | 77 | 27 / 50 / 23 | 0.351 [0.239, 0.474] | 0.540 [0.387, 0.694] | 0.425 [0.301, 0.545] |
| HG008T_ONT | SNV | d9 | segmental duplications | 46 | 17 / 29 / 31 | 0.370 [0.222, 0.562] | 0.354 [0.200, 0.548] | 0.362 [0.225, 0.525] |
| HG008T_ONT | SNV | d9 | tandem repeats and homopolymers | 107 | 36 / 71 / 51 | 0.336 [0.241, 0.440] | 0.414 [0.299, 0.532] | 0.371 [0.274, 0.466] |
| HG008T_ONT | SNV | d9 | homopolymers >= 7 bp | 55 | 15 / 40 / 27 | 0.273 [0.159, 0.404] | 0.357 [0.214, 0.500] | 0.309 [0.189, 0.423] |
| HG008T_ONT | SNV | d9 | tandem repeats <= 50 bp | 25 | 11 / 14 / 3 | 0.440 [0.222, 0.652] | 0.786 [0.500, 1.000] | 0.564 [0.323, 0.740] |
| HG008T_ONT | SNV | d9 | satellites | 0 | 0 / 0 / 2 | - [-, -] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_ONT | SNV | d9 | L1H > 500 bp | 3 | 1 / 2 / 0 | 0.333 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.500 [0.286, 1.000] |
| HG008T_ONT | SNV | d9 | GC < 25% or > 65% | 82 | 44 / 38 / 15 | 0.537 [0.409, 0.662] | 0.746 [0.625, 0.851] | 0.624 [0.504, 0.722] |
| HG008T_ONT | SNV | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | SNV | d46 | not in any difficult region | 1,006 | 481 / 525 / 36 | 0.478 [0.436, 0.519] | 0.930 [0.908, 0.950] | 0.632 [0.594, 0.666] |
| HG008T_ONT | SNV | d46 | all difficult regions | 381 | 118 / 263 / 62 | 0.310 [0.249, 0.380] | 0.656 [0.570, 0.736] | 0.421 [0.351, 0.495] |
| HG008T_ONT | SNV | d46 | low mappability or segdup | 159 | 42 / 117 / 25 | 0.264 [0.179, 0.381] | 0.627 [0.471, 0.779] | 0.372 [0.263, 0.498] |
| HG008T_ONT | SNV | d46 | low mappability | 135 | 34 / 101 / 16 | 0.252 [0.159, 0.380] | 0.680 [0.536, 0.814] | 0.368 [0.245, 0.512] |
| HG008T_ONT | SNV | d46 | segmental duplications | 107 | 26 / 81 / 22 | 0.243 [0.145, 0.406] | 0.542 [0.364, 0.733] | 0.335 [0.214, 0.496] |
| HG008T_ONT | SNV | d46 | tandem repeats and homopolymers | 174 | 40 / 134 / 47 | 0.230 [0.164, 0.305] | 0.460 [0.340, 0.582] | 0.307 [0.223, 0.393] |
| HG008T_ONT | SNV | d46 | homopolymers >= 7 bp | 78 | 15 / 63 / 27 | 0.192 [0.108, 0.292] | 0.357 [0.216, 0.500] | 0.250 [0.148, 0.354] |
| HG008T_ONT | SNV | d46 | tandem repeats <= 50 bp | 36 | 11 / 25 / 3 | 0.306 [0.143, 0.484] | 0.786 [0.500, 1.000] | 0.440 [0.229, 0.622] |
| HG008T_ONT | SNV | d46 | satellites | 5 | 1 / 4 / 1 | 0.200 [0.000, 0.333] | 0.500 [0.000, 1.000] | 0.286 [0.200, 0.500] |
| HG008T_ONT | SNV | d46 | L1H > 500 bp | 2 | 0 / 2 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_ONT | SNV | d46 | GC < 25% or > 65% | 104 | 47 / 57 / 12 | 0.452 [0.342, 0.558] | 0.797 [0.667, 0.908] | 0.577 [0.460, 0.676] |
| HG008T_ONT | SNV | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | SNV | linear | not in any difficult region | 944 | 464 / 480 / 53 | 0.492 [0.451, 0.532] | 0.897 [0.873, 0.923] | 0.635 [0.601, 0.669] |
| HG008T_PacBio | SNV | linear | all difficult regions | 601 | 112 / 489 / 68 | 0.186 [0.105, 0.314] | 0.622 [0.535, 0.705] | 0.287 [0.179, 0.425] |
| HG008T_PacBio | SNV | linear | low mappability or segdup | 385 | 40 / 345 / 27 | 0.104 [0.046, 0.262] | 0.597 [0.439, 0.746] | 0.177 [0.085, 0.373] |
| HG008T_PacBio | SNV | linear | low mappability | 361 | 32 / 329 / 18 | 0.089 [0.037, 0.245] | 0.640 [0.490, 0.782] | 0.156 [0.069, 0.361] |
| HG008T_PacBio | SNV | linear | segmental duplications | 301 | 25 / 276 / 23 | 0.083 [0.030, 0.275] | 0.521 [0.333, 0.714] | 0.143 [0.057, 0.384] |
| HG008T_PacBio | SNV | linear | tandem repeats and homopolymers | 399 | 41 / 358 / 46 | 0.103 [0.047, 0.242] | 0.471 [0.351, 0.595] | 0.169 [0.084, 0.330] |
| HG008T_PacBio | SNV | linear | homopolymers >= 7 bp | 85 | 16 / 69 / 26 | 0.188 [0.104, 0.295] | 0.381 [0.235, 0.522] | 0.252 [0.152, 0.361] |
| HG008T_PacBio | SNV | linear | tandem repeats <= 50 bp | 37 | 11 / 26 / 3 | 0.297 [0.138, 0.457] | 0.786 [0.500, 1.000] | 0.431 [0.222, 0.600] |
| HG008T_PacBio | SNV | linear | satellites | 205 | 1 / 204 / 1 | 0.005 [0.000, 0.333] | 0.500 [0.000, 1.000] | 0.010 [0.004, 0.500] |
| HG008T_PacBio | SNV | linear | L1H > 500 bp | 2 | 0 / 2 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_PacBio | SNV | linear | GC < 25% or > 65% | 106 | 40 / 66 / 19 | 0.377 [0.268, 0.482] | 0.678 [0.532, 0.804] | 0.485 [0.361, 0.589] |
| HG008T_PacBio | SNV | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | SNV | d9 | not in any difficult region | 932 | 472 / 460 / 45 | 0.506 [0.466, 0.547] | 0.913 [0.889, 0.935] | 0.651 [0.614, 0.685] |
| HG008T_PacBio | SNV | d9 | all difficult regions | 228 | 105 / 123 / 75 | 0.461 [0.385, 0.536] | 0.583 [0.500, 0.663] | 0.515 [0.444, 0.583] |
| HG008T_PacBio | SNV | d9 | low mappability or segdup | 84 | 33 / 51 / 34 | 0.393 [0.290, 0.500] | 0.493 [0.349, 0.661] | 0.437 [0.326, 0.550] |
| HG008T_PacBio | SNV | d9 | low mappability | 72 | 28 / 44 / 22 | 0.389 [0.272, 0.506] | 0.560 [0.406, 0.721] | 0.459 [0.338, 0.574] |
| HG008T_PacBio | SNV | d9 | segmental duplications | 39 | 18 / 21 / 30 | 0.462 [0.308, 0.622] | 0.375 [0.220, 0.575] | 0.414 [0.264, 0.567] |
| HG008T_PacBio | SNV | d9 | tandem repeats and homopolymers | 85 | 35 / 50 / 52 | 0.412 [0.303, 0.528] | 0.402 [0.293, 0.517] | 0.407 [0.303, 0.508] |
| HG008T_PacBio | SNV | d9 | homopolymers >= 7 bp | 39 | 15 / 24 / 27 | 0.385 [0.233, 0.542] | 0.357 [0.217, 0.500] | 0.370 [0.232, 0.494] |
| HG008T_PacBio | SNV | d9 | tandem repeats <= 50 bp | 20 | 10 / 10 / 4 | 0.500 [0.250, 0.714] | 0.714 [0.400, 0.938] | 0.588 [0.333, 0.765] |
| HG008T_PacBio | SNV | d9 | satellites | 0 | 0 / 0 / 2 | - [-, -] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_PacBio | SNV | d9 | L1H > 500 bp | 2 | 0 / 2 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_PacBio | SNV | d9 | GC < 25% or > 65% | 81 | 45 / 36 / 14 | 0.556 [0.426, 0.683] | 0.763 [0.636, 0.871] | 0.643 [0.519, 0.745] |
| HG008T_PacBio | SNV | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | SNV | d46 | not in any difficult region | 974 | 474 / 500 / 43 | 0.487 [0.446, 0.527] | 0.917 [0.894, 0.940] | 0.636 [0.600, 0.670] |
| HG008T_PacBio | SNV | d46 | all difficult regions | 258 | 108 / 150 / 72 | 0.419 [0.349, 0.490] | 0.600 [0.513, 0.681] | 0.493 [0.424, 0.561] |
| HG008T_PacBio | SNV | d46 | low mappability or segdup | 81 | 28 / 53 / 39 | 0.346 [0.244, 0.455] | 0.418 [0.288, 0.577] | 0.378 [0.272, 0.491] |
| HG008T_PacBio | SNV | d46 | low mappability | 67 | 24 / 43 / 26 | 0.358 [0.239, 0.478] | 0.480 [0.333, 0.640] | 0.410 [0.286, 0.529] |
| HG008T_PacBio | SNV | d46 | segmental duplications | 34 | 12 / 22 / 36 | 0.353 [0.200, 0.514] | 0.250 [0.120, 0.429] | 0.293 [0.150, 0.446] |
| HG008T_PacBio | SNV | d46 | tandem repeats and homopolymers | 114 | 39 / 75 / 48 | 0.342 [0.252, 0.439] | 0.448 [0.328, 0.561] | 0.388 [0.295, 0.478] |
| HG008T_PacBio | SNV | d46 | homopolymers >= 7 bp | 51 | 16 / 35 / 26 | 0.314 [0.188, 0.449] | 0.381 [0.234, 0.519] | 0.344 [0.218, 0.457] |
| HG008T_PacBio | SNV | d46 | tandem repeats <= 50 bp | 25 | 11 / 14 / 3 | 0.440 [0.227, 0.643] | 0.786 [0.500, 1.000] | 0.564 [0.323, 0.737] |
| HG008T_PacBio | SNV | d46 | satellites | 2 | 1 / 1 / 1 | 0.500 [0.000, 1.000] | 0.500 [0.000, 1.000] | 0.500 [0.333, 1.000] |
| HG008T_PacBio | SNV | d46 | L1H > 500 bp | 3 | 1 / 2 / 0 | 0.333 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.500 [0.286, 1.000] |
| HG008T_PacBio | SNV | d46 | GC < 25% or > 65% | 93 | 49 / 44 / 10 | 0.527 [0.406, 0.643] | 0.831 [0.714, 0.925] | 0.645 [0.529, 0.743] |
| HG008T_PacBio | SNV | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | INDEL | linear | not in any difficult region | 126 | 4 / 122 / 3 | 0.032 [0.008, 0.064] | 0.571 [0.143, 1.000] | 0.060 [0.016, 0.119] |
| COLO829T_Illumina | INDEL | linear | all difficult regions | 1,430 | 25 / 1405 / 99 | 0.017 [0.009, 0.029] | 0.202 [0.125, 0.279] | 0.032 [0.017, 0.051] |
| COLO829T_Illumina | INDEL | linear | low mappability or segdup | 638 | 6 / 632 / 23 | 0.009 [0.002, 0.028] | 0.207 [0.069, 0.421] | 0.018 [0.004, 0.052] |
| COLO829T_Illumina | INDEL | linear | low mappability | 597 | 4 / 593 / 22 | 0.007 [0.001, 0.025] | 0.154 [0.028, 0.400] | 0.013 [0.003, 0.044] |
| COLO829T_Illumina | INDEL | linear | segmental duplications | 170 | 3 / 167 / 19 | 0.018 [0.000, 0.043] | 0.136 [0.000, 0.333] | 0.031 [0.009, 0.074] |
| COLO829T_Illumina | INDEL | linear | tandem repeats and homopolymers | 1,369 | 19 / 1350 / 89 | 0.014 [0.007, 0.024] | 0.176 [0.103, 0.252] | 0.026 [0.013, 0.043] |
| COLO829T_Illumina | INDEL | linear | homopolymers >= 7 bp | 546 | 15 / 531 / 60 | 0.027 [0.014, 0.043] | 0.200 [0.107, 0.299] | 0.048 [0.025, 0.074] |
| COLO829T_Illumina | INDEL | linear | tandem repeats <= 50 bp | 62 | 3 / 59 / 18 | 0.048 [0.000, 0.109] | 0.143 [0.000, 0.300] | 0.072 [0.023, 0.154] |
| COLO829T_Illumina | INDEL | linear | satellites | 431 | 0 / 431 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_Illumina | INDEL | linear | L1H > 500 bp | 6 | 1 / 5 / 0 | 0.167 [0.000, 0.667] | 1.000 [1.000, 1.000] | 0.286 [0.167, 0.857] |
| COLO829T_Illumina | INDEL | linear | GC < 25% or > 65% | 227 | 4 / 223 / 10 | 0.018 [0.000, 0.038] | 0.286 [0.000, 0.572] | 0.033 [0.008, 0.071] |
| COLO829T_Illumina | INDEL | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | INDEL | d9 | not in any difficult region | 71 | 7 / 64 / 0 | 0.099 [0.036, 0.173] | 1.000 [1.000, 1.000] | 0.179 [0.069, 0.295] |
| COLO829T_Illumina | INDEL | d9 | all difficult regions | 487 | 29 / 458 / 95 | 0.060 [0.037, 0.082] | 0.234 [0.152, 0.318] | 0.095 [0.060, 0.129] |
| COLO829T_Illumina | INDEL | d9 | low mappability or segdup | 41 | 5 / 36 / 24 | 0.122 [0.029, 0.241] | 0.172 [0.042, 0.412] | 0.143 [0.037, 0.281] |
| COLO829T_Illumina | INDEL | d9 | low mappability | 33 | 4 / 29 / 22 | 0.121 [0.024, 0.250] | 0.154 [0.028, 0.400] | 0.136 [0.035, 0.273] |
| COLO829T_Illumina | INDEL | d9 | segmental duplications | 24 | 2 / 22 / 20 | 0.083 [0.000, 0.231] | 0.091 [0.000, 0.319] | 0.087 [0.032, 0.253] |
| COLO829T_Illumina | INDEL | d9 | tandem repeats and homopolymers | 460 | 24 / 436 / 84 | 0.052 [0.032, 0.072] | 0.222 [0.143, 0.300] | 0.085 [0.053, 0.116] |
| COLO829T_Illumina | INDEL | d9 | homopolymers >= 7 bp | 374 | 20 / 354 / 55 | 0.053 [0.031, 0.077] | 0.267 [0.162, 0.370] | 0.089 [0.052, 0.126] |
| COLO829T_Illumina | INDEL | d9 | tandem repeats <= 50 bp | 27 | 2 / 25 / 19 | 0.074 [0.000, 0.200] | 0.095 [0.000, 0.231] | 0.083 [0.036, 0.206] |
| COLO829T_Illumina | INDEL | d9 | satellites | 3 | 0 / 3 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_Illumina | INDEL | d9 | L1H > 500 bp | 1 | 1 / 0 / 0 | 1.000 [1.000, 1.000] | 1.000 [1.000, 1.000] | 1.000 [1.000, 1.000] |
| COLO829T_Illumina | INDEL | d9 | GC < 25% or > 65% | 70 | 4 / 66 / 10 | 0.057 [0.000, 0.123] | 0.286 [0.000, 0.572] | 0.095 [0.024, 0.198] |
| COLO829T_Illumina | INDEL | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_Illumina | INDEL | d46 | not in any difficult region | 50 | 7 / 43 / 0 | 0.140 [0.051, 0.241] | 1.000 [1.000, 1.000] | 0.246 [0.098, 0.388] |
| COLO829T_Illumina | INDEL | d46 | all difficult regions | 366 | 25 / 341 / 99 | 0.068 [0.042, 0.098] | 0.202 [0.128, 0.279] | 0.102 [0.064, 0.143] |
| COLO829T_Illumina | INDEL | d46 | low mappability or segdup | 47 | 3 / 44 / 26 | 0.064 [0.000, 0.156] | 0.103 [0.000, 0.286] | 0.079 [0.022, 0.186] |
| COLO829T_Illumina | INDEL | d46 | low mappability | 38 | 2 / 36 / 24 | 0.053 [0.000, 0.150] | 0.077 [0.000, 0.250] | 0.062 [0.024, 0.174] |
| COLO829T_Illumina | INDEL | d46 | segmental duplications | 23 | 1 / 22 / 21 | 0.043 [0.000, 0.167] | 0.045 [0.000, 0.200] | 0.044 [0.029, 0.182] |
| COLO829T_Illumina | INDEL | d46 | tandem repeats and homopolymers | 348 | 22 / 326 / 86 | 0.063 [0.037, 0.091] | 0.204 [0.129, 0.279] | 0.096 [0.058, 0.137] |
| COLO829T_Illumina | INDEL | d46 | homopolymers >= 7 bp | 270 | 18 / 252 / 57 | 0.067 [0.037, 0.099] | 0.240 [0.140, 0.343] | 0.104 [0.059, 0.150] |
| COLO829T_Illumina | INDEL | d46 | tandem repeats <= 50 bp | 21 | 2 / 19 / 19 | 0.095 [0.000, 0.250] | 0.095 [0.000, 0.231] | 0.095 [0.042, 0.230] |
| COLO829T_Illumina | INDEL | d46 | satellites | 12 | 0 / 12 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_Illumina | INDEL | d46 | L1H > 500 bp | 0 | 0 / 0 / 1 | - [-, -] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_Illumina | INDEL | d46 | GC < 25% or > 65% | 52 | 4 / 48 / 10 | 0.077 [0.000, 0.170] | 0.286 [0.000, 0.572] | 0.121 [0.032, 0.255] |
| COLO829T_Illumina | INDEL | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | INDEL | linear | not in any difficult region | 678 | 6 / 672 / 1 | 0.009 [0.003, 0.017] | 0.857 [0.500, 1.000] | 0.018 [0.006, 0.033] |
| COLO829T_ONT | INDEL | linear | all difficult regions | 17,739 | 31 / 17708 / 93 | 0.002 [0.001, 0.002] | 0.250 [0.174, 0.328] | 0.003 [0.002, 0.005] |
| COLO829T_ONT | INDEL | linear | low mappability or segdup | 3,599 | 7 / 3592 / 22 | 0.002 [0.000, 0.005] | 0.241 [0.077, 0.444] | 0.004 [0.001, 0.009] |
| COLO829T_ONT | INDEL | linear | low mappability | 3,084 | 5 / 3079 / 21 | 0.002 [0.000, 0.004] | 0.192 [0.056, 0.400] | 0.003 [0.001, 0.009] |
| COLO829T_ONT | INDEL | linear | segmental duplications | 1,715 | 4 / 1711 / 18 | 0.002 [0.000, 0.005] | 0.182 [0.000, 0.375] | 0.005 [0.001, 0.010] |
| COLO829T_ONT | INDEL | linear | tandem repeats and homopolymers | 17,367 | 24 / 17343 / 84 | 0.001 [0.001, 0.002] | 0.222 [0.140, 0.308] | 0.003 [0.002, 0.004] |
| COLO829T_ONT | INDEL | linear | homopolymers >= 7 bp | 14,428 | 19 / 14409 / 56 | 0.001 [0.001, 0.002] | 0.253 [0.147, 0.359] | 0.003 [0.001, 0.004] |
| COLO829T_ONT | INDEL | linear | tandem repeats <= 50 bp | 341 | 3 / 338 / 18 | 0.009 [0.000, 0.020] | 0.143 [0.000, 0.300] | 0.017 [0.005, 0.037] |
| COLO829T_ONT | INDEL | linear | satellites | 1,667 | 0 / 1667 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_ONT | INDEL | linear | L1H > 500 bp | 13 | 1 / 12 / 0 | 0.077 [0.000, 0.267] | 1.000 [1.000, 1.000] | 0.143 [0.104, 0.444] |
| COLO829T_ONT | INDEL | linear | GC < 25% or > 65% | 1,314 | 3 / 1311 / 11 | 0.002 [0.000, 0.005] | 0.214 [0.000, 0.444] | 0.005 [0.001, 0.010] |
| COLO829T_ONT | INDEL | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | INDEL | d9 | not in any difficult region | 710 | 6 / 704 / 1 | 0.008 [0.003, 0.016] | 0.857 [0.500, 1.000] | 0.017 [0.005, 0.032] |
| COLO829T_ONT | INDEL | d9 | all difficult regions | 20,397 | 26 / 20370 / 98 | 0.001 [0.001, 0.002] | 0.210 [0.136, 0.289] | 0.003 [0.002, 0.004] |
| COLO829T_ONT | INDEL | d9 | low mappability or segdup | 3,661 | 5 / 3656 / 24 | 0.001 [0.000, 0.003] | 0.172 [0.042, 0.412] | 0.003 [0.001, 0.007] |
| COLO829T_ONT | INDEL | d9 | low mappability | 3,125 | 4 / 3121 / 22 | 0.001 [0.000, 0.004] | 0.154 [0.028, 0.400] | 0.003 [0.001, 0.007] |
| COLO829T_ONT | INDEL | d9 | segmental duplications | 1,781 | 2 / 1779 / 20 | 0.001 [0.000, 0.004] | 0.091 [0.000, 0.319] | 0.002 [0.001, 0.007] |
| COLO829T_ONT | INDEL | d9 | tandem repeats and homopolymers | 20,065 | 21 / 20043 / 87 | 0.001 [0.001, 0.002] | 0.194 [0.124, 0.273] | 0.002 [0.001, 0.003] |
| COLO829T_ONT | INDEL | d9 | homopolymers >= 7 bp | 17,135 | 18 / 17116 / 57 | 0.001 [0.001, 0.002] | 0.240 [0.141, 0.342] | 0.002 [0.001, 0.003] |
| COLO829T_ONT | INDEL | d9 | tandem repeats <= 50 bp | 296 | 1 / 295 / 20 | 0.003 [0.000, 0.011] | 0.048 [0.000, 0.143] | 0.006 [0.006, 0.023] |
| COLO829T_ONT | INDEL | d9 | satellites | 1,654 | 0 / 1654 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_ONT | INDEL | d9 | L1H > 500 bp | 7 | 1 / 6 / 0 | 0.143 [0.000, 0.429] | 1.000 [1.000, 1.000] | 0.250 [0.167, 0.667] |
| COLO829T_ONT | INDEL | d9 | GC < 25% or > 65% | 1,411 | 4 / 1406 / 10 | 0.004 [0.001, 0.008] | 0.286 [0.077, 0.556] | 0.007 [0.001, 0.015] |
| COLO829T_ONT | INDEL | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_ONT | INDEL | d46 | not in any difficult region | 501 | 5 / 496 / 2 | 0.010 [0.002, 0.020] | 0.714 [0.333, 1.000] | 0.020 [0.004, 0.039] |
| COLO829T_ONT | INDEL | d46 | all difficult regions | 27,386 | 25 / 27361 / 99 | 0.001 [0.001, 0.001] | 0.202 [0.127, 0.279] | 0.002 [0.001, 0.003] |
| COLO829T_ONT | INDEL | d46 | low mappability or segdup | 6,601 | 4 / 6597 / 25 | 0.001 [0.000, 0.002] | 0.138 [0.026, 0.350] | 0.001 [0.000, 0.004] |
| COLO829T_ONT | INDEL | d46 | low mappability | 5,675 | 3 / 5672 / 23 | 0.001 [0.000, 0.002] | 0.115 [0.000, 0.320] | 0.001 [0.000, 0.004] |
| COLO829T_ONT | INDEL | d46 | segmental duplications | 3,734 | 1 / 3733 / 21 | 0.000 [0.000, 0.002] | 0.045 [0.000, 0.200] | 0.001 [0.000, 0.004] |
| COLO829T_ONT | INDEL | d46 | tandem repeats and homopolymers | 26,941 | 21 / 26920 / 87 | 0.001 [0.000, 0.001] | 0.194 [0.122, 0.271] | 0.002 [0.001, 0.002] |
| COLO829T_ONT | INDEL | d46 | homopolymers >= 7 bp | 21,618 | 17 / 21601 / 58 | 0.001 [0.000, 0.001] | 0.227 [0.130, 0.325] | 0.002 [0.001, 0.002] |
| COLO829T_ONT | INDEL | d46 | tandem repeats <= 50 bp | 355 | 2 / 353 / 19 | 0.006 [0.000, 0.015] | 0.095 [0.000, 0.231] | 0.011 [0.005, 0.028] |
| COLO829T_ONT | INDEL | d46 | satellites | 2,550 | 0 / 2550 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_ONT | INDEL | d46 | L1H > 500 bp | 12 | 1 / 11 / 0 | 0.083 [0.000, 0.286] | 1.000 [1.000, 1.000] | 0.154 [0.105, 0.466] |
| COLO829T_ONT | INDEL | d46 | GC < 25% or > 65% | 1,623 | 3 / 1620 / 11 | 0.002 [0.000, 0.004] | 0.214 [0.000, 0.444] | 0.004 [0.001, 0.008] |
| COLO829T_ONT | INDEL | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | INDEL | linear | not in any difficult region | 135 | 6 / 129 / 1 | 0.044 [0.014, 0.083] | 0.857 [0.500, 1.000] | 0.085 [0.028, 0.153] |
| COLO829T_fiberseq | INDEL | linear | all difficult regions | 1,505 | 33 / 1472 / 91 | 0.022 [0.014, 0.031] | 0.266 [0.188, 0.350] | 0.041 [0.026, 0.056] |
| COLO829T_fiberseq | INDEL | linear | low mappability or segdup | 412 | 6 / 406 / 23 | 0.015 [0.004, 0.032] | 0.207 [0.069, 0.421] | 0.027 [0.008, 0.058] |
| COLO829T_fiberseq | INDEL | linear | low mappability | 357 | 5 / 352 / 21 | 0.014 [0.003, 0.033] | 0.192 [0.056, 0.400] | 0.026 [0.007, 0.059] |
| COLO829T_fiberseq | INDEL | linear | segmental duplications | 278 | 3 / 275 / 19 | 0.011 [0.000, 0.028] | 0.136 [0.000, 0.333] | 0.020 [0.006, 0.049] |
| COLO829T_fiberseq | INDEL | linear | tandem repeats and homopolymers | 1,414 | 27 / 1387 / 81 | 0.019 [0.012, 0.027] | 0.250 [0.172, 0.339] | 0.035 [0.023, 0.050] |
| COLO829T_fiberseq | INDEL | linear | homopolymers >= 7 bp | 829 | 21 / 808 / 54 | 0.025 [0.015, 0.037] | 0.280 [0.176, 0.390] | 0.046 [0.027, 0.066] |
| COLO829T_fiberseq | INDEL | linear | tandem repeats <= 50 bp | 90 | 3 / 87 / 18 | 0.033 [0.000, 0.073] | 0.143 [0.000, 0.300] | 0.054 [0.017, 0.115] |
| COLO829T_fiberseq | INDEL | linear | satellites | 122 | 0 / 122 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_fiberseq | INDEL | linear | L1H > 500 bp | 7 | 1 / 6 / 0 | 0.143 [0.000, 0.500] | 1.000 [1.000, 1.000] | 0.250 [0.167, 0.719] |
| COLO829T_fiberseq | INDEL | linear | GC < 25% or > 65% | 221 | 5 / 216 / 9 | 0.023 [0.004, 0.047] | 0.357 [0.083, 0.667] | 0.043 [0.009, 0.086] |
| COLO829T_fiberseq | INDEL | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | INDEL | d9 | not in any difficult region | 87 | 7 / 80 / 0 | 0.080 [0.029, 0.143] | 1.000 [1.000, 1.000] | 0.149 [0.056, 0.250] |
| COLO829T_fiberseq | INDEL | d9 | all difficult regions | 2,685 | 28 / 2657 / 96 | 0.010 [0.006, 0.015] | 0.226 [0.147, 0.305] | 0.020 [0.012, 0.028] |
| COLO829T_fiberseq | INDEL | d9 | low mappability or segdup | 297 | 5 / 292 / 24 | 0.017 [0.003, 0.038] | 0.172 [0.042, 0.412] | 0.031 [0.007, 0.067] |
| COLO829T_fiberseq | INDEL | d9 | low mappability | 221 | 4 / 217 / 22 | 0.018 [0.003, 0.044] | 0.154 [0.028, 0.400] | 0.032 [0.007, 0.077] |
| COLO829T_fiberseq | INDEL | d9 | segmental duplications | 205 | 2 / 203 / 20 | 0.010 [0.000, 0.027] | 0.091 [0.000, 0.319] | 0.018 [0.007, 0.050] |
| COLO829T_fiberseq | INDEL | d9 | tandem repeats and homopolymers | 2,631 | 23 / 2608 / 85 | 0.009 [0.005, 0.012] | 0.213 [0.140, 0.285] | 0.017 [0.010, 0.024] |
| COLO829T_fiberseq | INDEL | d9 | homopolymers >= 7 bp | 2,241 | 20 / 2221 / 55 | 0.009 [0.005, 0.013] | 0.267 [0.171, 0.361] | 0.017 [0.010, 0.025] |
| COLO829T_fiberseq | INDEL | d9 | tandem repeats <= 50 bp | 74 | 1 / 73 / 20 | 0.014 [0.000, 0.047] | 0.048 [0.000, 0.143] | 0.021 [0.018, 0.076] |
| COLO829T_fiberseq | INDEL | d9 | satellites | 53 | 0 / 53 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_fiberseq | INDEL | d9 | L1H > 500 bp | 3 | 1 / 2 / 0 | 0.333 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.500 [0.286, 1.000] |
| COLO829T_fiberseq | INDEL | d9 | GC < 25% or > 65% | 271 | 4 / 267 / 10 | 0.015 [0.000, 0.033] | 0.286 [0.000, 0.572] | 0.028 [0.007, 0.062] |
| COLO829T_fiberseq | INDEL | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| COLO829T_fiberseq | INDEL | d46 | not in any difficult region | 77 | 6 / 71 / 1 | 0.078 [0.023, 0.145] | 0.857 [0.500, 1.000] | 0.143 [0.046, 0.250] |
| COLO829T_fiberseq | INDEL | d46 | all difficult regions | 3,026 | 28 / 2998 / 96 | 0.009 [0.006, 0.013] | 0.226 [0.146, 0.308] | 0.018 [0.011, 0.025] |
| COLO829T_fiberseq | INDEL | d46 | low mappability or segdup | 162 | 4 / 158 / 25 | 0.025 [0.005, 0.058] | 0.138 [0.026, 0.350] | 0.042 [0.010, 0.094] |
| COLO829T_fiberseq | INDEL | d46 | low mappability | 119 | 3 / 116 / 23 | 0.025 [0.000, 0.066] | 0.115 [0.000, 0.320] | 0.041 [0.011, 0.105] |
| COLO829T_fiberseq | INDEL | d46 | segmental duplications | 86 | 1 / 85 / 21 | 0.012 [0.000, 0.045] | 0.045 [0.000, 0.200] | 0.019 [0.013, 0.077] |
| COLO829T_fiberseq | INDEL | d46 | tandem repeats and homopolymers | 3,004 | 24 / 2980 / 84 | 0.008 [0.005, 0.011] | 0.222 [0.143, 0.300] | 0.015 [0.009, 0.022] |
| COLO829T_fiberseq | INDEL | d46 | homopolymers >= 7 bp | 2,687 | 20 / 2667 / 55 | 0.007 [0.004, 0.011] | 0.267 [0.162, 0.370] | 0.014 [0.008, 0.021] |
| COLO829T_fiberseq | INDEL | d46 | tandem repeats <= 50 bp | 83 | 2 / 81 / 19 | 0.024 [0.000, 0.065] | 0.095 [0.000, 0.231] | 0.038 [0.017, 0.102] |
| COLO829T_fiberseq | INDEL | d46 | satellites | 34 | 0 / 34 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| COLO829T_fiberseq | INDEL | d46 | L1H > 500 bp | 2 | 1 / 1 / 0 | 0.500 [0.000, 1.000] | 1.000 [1.000, 1.000] | 0.667 [0.400, 1.000] |
| COLO829T_fiberseq | INDEL | d46 | GC < 25% or > 65% | 275 | 4 / 271 / 10 | 0.015 [0.000, 0.032] | 0.286 [0.000, 0.572] | 0.028 [0.007, 0.060] |
| COLO829T_fiberseq | INDEL | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | linear | not in any difficult region | 147 | 26 / 121 / 9 | 0.177 [0.115, 0.243] | 0.743 [0.588, 0.895] | 0.286 [0.194, 0.376] |
| HG008T_Illumina | INDEL | linear | all difficult regions | 1,353 | 126 / 1227 / 426 | 0.093 [0.076, 0.110] | 0.228 [0.197, 0.260] | 0.132 [0.111, 0.153] |
| HG008T_Illumina | INDEL | linear | low mappability or segdup | 78 | 8 / 70 / 33 | 0.103 [0.043, 0.175] | 0.195 [0.089, 0.308] | 0.134 [0.058, 0.217] |
| HG008T_Illumina | INDEL | linear | low mappability | 49 | 5 / 44 / 23 | 0.102 [0.024, 0.190] | 0.179 [0.043, 0.333] | 0.130 [0.033, 0.233] |
| HG008T_Illumina | INDEL | linear | segmental duplications | 50 | 5 / 45 / 25 | 0.100 [0.024, 0.188] | 0.167 [0.050, 0.292] | 0.125 [0.037, 0.217] |
| HG008T_Illumina | INDEL | linear | tandem repeats and homopolymers | 1,339 | 124 / 1215 / 423 | 0.093 [0.076, 0.110] | 0.227 [0.194, 0.259] | 0.131 [0.110, 0.153] |
| HG008T_Illumina | INDEL | linear | homopolymers >= 7 bp | 1,070 | 119 / 951 / 346 | 0.111 [0.091, 0.132] | 0.256 [0.220, 0.292] | 0.155 [0.130, 0.181] |
| HG008T_Illumina | INDEL | linear | tandem repeats <= 50 bp | 59 | 1 / 58 / 35 | 0.017 [0.000, 0.057] | 0.028 [0.000, 0.088] | 0.021 [0.018, 0.073] |
| HG008T_Illumina | INDEL | linear | satellites | 1 | 0 / 1 / 0 | 0.000 [0.000, 0.000] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | linear | L1H > 500 bp | 3 | 0 / 3 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_Illumina | INDEL | linear | GC < 25% or > 65% | 168 | 10 / 158 / 40 | 0.060 [0.028, 0.098] | 0.200 [0.098, 0.306] | 0.092 [0.044, 0.145] |
| HG008T_Illumina | INDEL | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d9 | not in any difficult region | 71 | 30 / 41 / 5 | 0.423 [0.310, 0.535] | 0.857 [0.731, 0.967] | 0.566 [0.447, 0.672] |
| HG008T_Illumina | INDEL | d9 | all difficult regions | 459 | 110 / 349 / 442 | 0.240 [0.199, 0.284] | 0.199 [0.167, 0.232] | 0.218 [0.183, 0.252] |
| HG008T_Illumina | INDEL | d9 | low mappability or segdup | 16 | 2 / 14 / 39 | 0.125 [0.000, 0.333] | 0.049 [0.000, 0.129] | 0.070 [0.029, 0.179] |
| HG008T_Illumina | INDEL | d9 | low mappability | 9 | 2 / 7 / 26 | 0.222 [0.000, 0.600] | 0.071 [0.000, 0.188] | 0.108 [0.043, 0.261] |
| HG008T_Illumina | INDEL | d9 | segmental duplications | 11 | 1 / 10 / 29 | 0.091 [0.000, 0.333] | 0.033 [0.000, 0.120] | 0.049 [0.036, 0.174] |
| HG008T_Illumina | INDEL | d9 | tandem repeats and homopolymers | 451 | 108 / 343 / 439 | 0.239 [0.198, 0.284] | 0.197 [0.165, 0.230] | 0.216 [0.182, 0.251] |
| HG008T_Illumina | INDEL | d9 | homopolymers >= 7 bp | 421 | 104 / 317 / 361 | 0.247 [0.202, 0.295] | 0.224 [0.188, 0.263] | 0.235 [0.197, 0.274] |
| HG008T_Illumina | INDEL | d9 | tandem repeats <= 50 bp | 15 | 3 / 12 / 33 | 0.200 [0.000, 0.438] | 0.083 [0.000, 0.182] | 0.118 [0.038, 0.250] |
| HG008T_Illumina | INDEL | d9 | satellites | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d9 | L1H > 500 bp | 0 | 0 / 0 / 1 | - [-, -] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_Illumina | INDEL | d9 | GC < 25% or > 65% | 44 | 6 / 38 / 44 | 0.136 [0.050, 0.220] | 0.120 [0.041, 0.207] | 0.128 [0.048, 0.206] |
| HG008T_Illumina | INDEL | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d46 | not in any difficult region | 87 | 31 / 56 / 4 | 0.356 [0.256, 0.462] | 0.886 [0.767, 0.975] | 0.508 [0.390, 0.613] |
| HG008T_Illumina | INDEL | d46 | all difficult regions | 439 | 110 / 329 / 442 | 0.251 [0.208, 0.297] | 0.199 [0.168, 0.233] | 0.222 [0.187, 0.257] |
| HG008T_Illumina | INDEL | d46 | low mappability or segdup | 19 | 2 / 17 / 39 | 0.105 [0.000, 0.263] | 0.049 [0.000, 0.129] | 0.067 [0.027, 0.169] |
| HG008T_Illumina | INDEL | d46 | low mappability | 11 | 2 / 9 / 26 | 0.182 [0.000, 0.444] | 0.071 [0.000, 0.188] | 0.103 [0.042, 0.250] |
| HG008T_Illumina | INDEL | d46 | segmental duplications | 12 | 1 / 11 / 29 | 0.083 [0.000, 0.286] | 0.033 [0.000, 0.120] | 0.048 [0.034, 0.167] |
| HG008T_Illumina | INDEL | d46 | tandem repeats and homopolymers | 425 | 107 / 318 / 440 | 0.252 [0.208, 0.299] | 0.196 [0.164, 0.230] | 0.220 [0.184, 0.256] |
| HG008T_Illumina | INDEL | d46 | homopolymers >= 7 bp | 392 | 102 / 290 / 363 | 0.260 [0.214, 0.310] | 0.219 [0.184, 0.257] | 0.238 [0.200, 0.278] |
| HG008T_Illumina | INDEL | d46 | tandem repeats <= 50 bp | 17 | 4 / 13 / 32 | 0.235 [0.050, 0.467] | 0.111 [0.024, 0.220] | 0.151 [0.040, 0.286] |
| HG008T_Illumina | INDEL | d46 | satellites | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_Illumina | INDEL | d46 | L1H > 500 bp | 0 | 0 / 0 / 1 | - [-, -] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_Illumina | INDEL | d46 | GC < 25% or > 65% | 48 | 7 / 41 / 43 | 0.146 [0.054, 0.239] | 0.140 [0.051, 0.234] | 0.143 [0.054, 0.227] |
| HG008T_Illumina | INDEL | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | linear | not in any difficult region | 366 | 30 / 336 / 5 | 0.082 [0.055, 0.112] | 0.857 [0.733, 0.969] | 0.150 [0.102, 0.198] |
| HG008T_ONT | INDEL | linear | all difficult regions | 11,145 | 108 / 11037 / 444 | 0.010 [0.008, 0.012] | 0.196 [0.164, 0.228] | 0.018 [0.015, 0.022] |
| HG008T_ONT | INDEL | linear | low mappability or segdup | 1,081 | 8 / 1073 / 33 | 0.007 [0.003, 0.013] | 0.195 [0.094, 0.310] | 0.014 [0.006, 0.024] |
| HG008T_ONT | INDEL | linear | low mappability | 811 | 6 / 805 / 22 | 0.007 [0.003, 0.013] | 0.214 [0.080, 0.357] | 0.014 [0.005, 0.026] |
| HG008T_ONT | INDEL | linear | segmental duplications | 813 | 6 / 807 / 24 | 0.007 [0.003, 0.013] | 0.200 [0.077, 0.333] | 0.014 [0.005, 0.025] |
| HG008T_ONT | INDEL | linear | tandem repeats and homopolymers | 11,038 | 104 / 10934 / 443 | 0.009 [0.008, 0.012] | 0.190 [0.158, 0.222] | 0.018 [0.014, 0.022] |
| HG008T_ONT | INDEL | linear | homopolymers >= 7 bp | 10,038 | 93 / 9945 / 372 | 0.009 [0.007, 0.012] | 0.200 [0.165, 0.236] | 0.018 [0.014, 0.022] |
| HG008T_ONT | INDEL | linear | tandem repeats <= 50 bp | 186 | 5 / 181 / 31 | 0.027 [0.006, 0.052] | 0.139 [0.032, 0.258] | 0.045 [0.010, 0.087] |
| HG008T_ONT | INDEL | linear | satellites | 25 | 0 / 25 / 0 | 0.000 [0.000, 0.000] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | linear | L1H > 500 bp | 7 | 0 / 7 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_ONT | INDEL | linear | GC < 25% or > 65% | 1,011 | 10 / 1001 / 40 | 0.010 [0.004, 0.016] | 0.200 [0.100, 0.304] | 0.019 [0.008, 0.031] |
| HG008T_ONT | INDEL | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d9 | not in any difficult region | 193 | 33 / 160 / 2 | 0.171 [0.121, 0.225] | 0.943 [0.854, 1.000] | 0.289 [0.213, 0.364] |
| HG008T_ONT | INDEL | d9 | all difficult regions | 16,964 | 101 / 16863 / 451 | 0.006 [0.005, 0.007] | 0.183 [0.154, 0.214] | 0.012 [0.009, 0.014] |
| HG008T_ONT | INDEL | d9 | low mappability or segdup | 1,087 | 6 / 1081 / 35 | 0.006 [0.002, 0.010] | 0.146 [0.057, 0.244] | 0.011 [0.004, 0.019] |
| HG008T_ONT | INDEL | d9 | low mappability | 753 | 4 / 749 / 24 | 0.005 [0.001, 0.011] | 0.143 [0.033, 0.273] | 0.010 [0.003, 0.021] |
| HG008T_ONT | INDEL | d9 | segmental duplications | 805 | 5 / 800 / 25 | 0.006 [0.001, 0.013] | 0.167 [0.050, 0.292] | 0.012 [0.003, 0.024] |
| HG008T_ONT | INDEL | d9 | tandem repeats and homopolymers | 16,904 | 98 / 16806 / 449 | 0.006 [0.005, 0.007] | 0.179 [0.149, 0.210] | 0.011 [0.009, 0.014] |
| HG008T_ONT | INDEL | d9 | homopolymers >= 7 bp | 16,316 | 93 / 16223 / 372 | 0.006 [0.004, 0.007] | 0.200 [0.167, 0.235] | 0.011 [0.009, 0.014] |
| HG008T_ONT | INDEL | d9 | tandem repeats <= 50 bp | 187 | 4 / 183 / 32 | 0.021 [0.005, 0.043] | 0.111 [0.026, 0.222] | 0.036 [0.009, 0.072] |
| HG008T_ONT | INDEL | d9 | satellites | 7 | 0 / 7 / 0 | 0.000 [0.000, 0.000] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d9 | L1H > 500 bp | 7 | 0 / 7 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_ONT | INDEL | d9 | GC < 25% or > 65% | 1,232 | 7 / 1225 / 43 | 0.006 [0.002, 0.010] | 0.140 [0.051, 0.234] | 0.011 [0.003, 0.019] |
| HG008T_ONT | INDEL | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d46 | not in any difficult region | 242 | 33 / 209 / 2 | 0.136 [0.096, 0.179] | 0.943 [0.854, 1.000] | 0.238 [0.172, 0.303] |
| HG008T_ONT | INDEL | d46 | all difficult regions | 18,075 | 106 / 17969 / 446 | 0.006 [0.005, 0.007] | 0.192 [0.161, 0.226] | 0.011 [0.009, 0.014] |
| HG008T_ONT | INDEL | d46 | low mappability or segdup | 1,283 | 8 / 1275 / 33 | 0.006 [0.002, 0.011] | 0.195 [0.067, 0.317] | 0.012 [0.004, 0.022] |
| HG008T_ONT | INDEL | d46 | low mappability | 943 | 6 / 937 / 22 | 0.006 [0.002, 0.012] | 0.214 [0.057, 0.400] | 0.012 [0.004, 0.023] |
| HG008T_ONT | INDEL | d46 | segmental duplications | 993 | 7 / 986 / 23 | 0.007 [0.002, 0.013] | 0.233 [0.065, 0.385] | 0.014 [0.004, 0.026] |
| HG008T_ONT | INDEL | d46 | tandem repeats and homopolymers | 17,973 | 102 / 17871 / 445 | 0.006 [0.004, 0.007] | 0.186 [0.156, 0.220] | 0.011 [0.009, 0.014] |
| HG008T_ONT | INDEL | d46 | homopolymers >= 7 bp | 17,314 | 97 / 17217 / 368 | 0.006 [0.004, 0.007] | 0.209 [0.175, 0.246] | 0.011 [0.009, 0.014] |
| HG008T_ONT | INDEL | d46 | tandem repeats <= 50 bp | 214 | 4 / 210 / 32 | 0.019 [0.004, 0.039] | 0.111 [0.024, 0.220] | 0.032 [0.008, 0.067] |
| HG008T_ONT | INDEL | d46 | satellites | 37 | 0 / 37 / 0 | 0.000 [0.000, 0.000] | - [-, -] | - [-, -] |
| HG008T_ONT | INDEL | d46 | L1H > 500 bp | 7 | 0 / 7 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_ONT | INDEL | d46 | GC < 25% or > 65% | 1,286 | 7 / 1279 / 43 | 0.005 [0.002, 0.010] | 0.140 [0.051, 0.234] | 0.010 [0.003, 0.019] |
| HG008T_ONT | INDEL | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | linear | not in any difficult region | 139 | 16 / 123 / 19 | 0.115 [0.062, 0.176] | 0.457 [0.296, 0.636] | 0.184 [0.105, 0.268] |
| HG008T_PacBio | INDEL | linear | all difficult regions | 1,534 | 134 / 1400 / 418 | 0.087 [0.073, 0.103] | 0.243 [0.209, 0.276] | 0.128 [0.108, 0.149] |
| HG008T_PacBio | INDEL | linear | low mappability or segdup | 116 | 9 / 107 / 32 | 0.078 [0.035, 0.129] | 0.220 [0.111, 0.348] | 0.115 [0.053, 0.182] |
| HG008T_PacBio | INDEL | linear | low mappability | 79 | 6 / 73 / 22 | 0.076 [0.024, 0.138] | 0.214 [0.074, 0.385] | 0.112 [0.038, 0.196] |
| HG008T_PacBio | INDEL | linear | segmental duplications | 66 | 5 / 61 / 25 | 0.076 [0.019, 0.140] | 0.167 [0.050, 0.292] | 0.104 [0.032, 0.182] |
| HG008T_PacBio | INDEL | linear | tandem repeats and homopolymers | 1,516 | 132 / 1384 / 415 | 0.087 [0.072, 0.103] | 0.241 [0.208, 0.274] | 0.128 [0.108, 0.148] |
| HG008T_PacBio | INDEL | linear | homopolymers >= 7 bp | 1,180 | 123 / 1057 / 342 | 0.104 [0.086, 0.124] | 0.265 [0.227, 0.304] | 0.150 [0.125, 0.174] |
| HG008T_PacBio | INDEL | linear | tandem repeats <= 50 bp | 88 | 4 / 84 / 32 | 0.045 [0.010, 0.096] | 0.111 [0.023, 0.220] | 0.065 [0.016, 0.130] |
| HG008T_PacBio | INDEL | linear | satellites | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | linear | L1H > 500 bp | 4 | 0 / 4 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_PacBio | INDEL | linear | GC < 25% or > 65% | 211 | 12 / 199 / 38 | 0.057 [0.030, 0.088] | 0.240 [0.130, 0.350] | 0.092 [0.049, 0.138] |
| HG008T_PacBio | INDEL | linear | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d9 | not in any difficult region | 117 | 33 / 84 / 2 | 0.282 [0.200, 0.370] | 0.943 [0.854, 1.000] | 0.434 [0.331, 0.535] |
| HG008T_PacBio | INDEL | d9 | all difficult regions | 2,567 | 117 / 2450 / 435 | 0.046 [0.035, 0.057] | 0.212 [0.179, 0.245] | 0.075 [0.060, 0.091] |
| HG008T_PacBio | INDEL | d9 | low mappability or segdup | 391 | 8 / 383 / 33 | 0.020 [0.006, 0.061] | 0.195 [0.074, 0.333] | 0.037 [0.011, 0.102] |
| HG008T_PacBio | INDEL | d9 | low mappability | 343 | 6 / 337 / 22 | 0.017 [0.005, 0.062] | 0.214 [0.077, 0.385] | 0.032 [0.009, 0.106] |
| HG008T_PacBio | INDEL | d9 | segmental duplications | 330 | 5 / 325 / 25 | 0.015 [0.002, 0.062] | 0.167 [0.033, 0.321] | 0.028 [0.006, 0.101] |
| HG008T_PacBio | INDEL | d9 | tandem repeats and homopolymers | 2,514 | 114 / 2400 / 433 | 0.045 [0.035, 0.056] | 0.208 [0.175, 0.242] | 0.074 [0.059, 0.091] |
| HG008T_PacBio | INDEL | d9 | homopolymers >= 7 bp | 2,101 | 109 / 1992 / 356 | 0.052 [0.042, 0.063] | 0.234 [0.199, 0.273] | 0.085 [0.069, 0.102] |
| HG008T_PacBio | INDEL | d9 | tandem repeats <= 50 bp | 87 | 4 / 83 / 32 | 0.046 [0.010, 0.095] | 0.111 [0.024, 0.220] | 0.065 [0.017, 0.129] |
| HG008T_PacBio | INDEL | d9 | satellites | 10 | 0 / 10 / 0 | 0.000 [0.000, 0.000] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d9 | L1H > 500 bp | 4 | 0 / 4 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_PacBio | INDEL | d9 | GC < 25% or > 65% | 266 | 9 / 257 / 41 | 0.034 [0.015, 0.059] | 0.180 [0.085, 0.286] | 0.057 [0.026, 0.096] |
| HG008T_PacBio | INDEL | d9 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d46 | not in any difficult region | 128 | 33 / 95 / 2 | 0.258 [0.180, 0.338] | 0.943 [0.854, 1.000] | 0.405 [0.300, 0.500] |
| HG008T_PacBio | INDEL | d46 | all difficult regions | 2,346 | 117 / 2229 / 435 | 0.050 [0.040, 0.060] | 0.212 [0.180, 0.245] | 0.081 [0.066, 0.096] |
| HG008T_PacBio | INDEL | d46 | low mappability or segdup | 217 | 7 / 210 / 34 | 0.032 [0.012, 0.059] | 0.171 [0.074, 0.288] | 0.054 [0.021, 0.096] |
| HG008T_PacBio | INDEL | d46 | low mappability | 172 | 5 / 167 / 23 | 0.029 [0.007, 0.062] | 0.179 [0.048, 0.345] | 0.050 [0.013, 0.103] |
| HG008T_PacBio | INDEL | d46 | segmental duplications | 163 | 4 / 159 / 26 | 0.025 [0.005, 0.051] | 0.133 [0.034, 0.240] | 0.041 [0.011, 0.083] |
| HG008T_PacBio | INDEL | d46 | tandem repeats and homopolymers | 2,306 | 114 / 2192 / 433 | 0.049 [0.040, 0.060] | 0.208 [0.176, 0.241] | 0.080 [0.065, 0.096] |
| HG008T_PacBio | INDEL | d46 | homopolymers >= 7 bp | 1,864 | 109 / 1755 / 356 | 0.058 [0.047, 0.071] | 0.234 [0.199, 0.274] | 0.094 [0.076, 0.112] |
| HG008T_PacBio | INDEL | d46 | tandem repeats <= 50 bp | 120 | 4 / 116 / 32 | 0.033 [0.007, 0.071] | 0.111 [0.024, 0.220] | 0.051 [0.013, 0.105] |
| HG008T_PacBio | INDEL | d46 | satellites | 23 | 0 / 23 / 0 | 0.000 [0.000, 0.000] | - [-, -] | - [-, -] |
| HG008T_PacBio | INDEL | d46 | L1H > 500 bp | 2 | 0 / 2 / 1 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | - [-, -] |
| HG008T_PacBio | INDEL | d46 | GC < 25% or > 65% | 298 | 9 / 289 / 41 | 0.030 [0.014, 0.050] | 0.180 [0.085, 0.286] | 0.052 [0.024, 0.085] |
| HG008T_PacBio | INDEL | d46 | CMRG falsely duplicated genes | 0 | 0 / 0 / 0 | - [-, -] | - [-, -] | - [-, -] |
