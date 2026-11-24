# E19 on the current models (strata_latest.py)

chr1, benchmark regions, PASS calls after the population filter. Stratum = germline variants of the normal within 50 bp of the variant (0 / 1 / ≥2); recall by the truth variant's stratum, precision by the call's stratum. Models, PASS rules and directories: the docstring of strata_latest.py. One seed (seed 0), except INDEL on HPRC v2.1: three seeds, shown as the best seed (highest whole-chr1 PASS F1 on that tumor and platform) and as the mean of the three. ΔF1 = HPRC − GRCh38-only of the same platform. pending = the call set does not exist yet.

## 1. chr1 INDEL, HG008T, six-tumor model (truth INDELs 0 / 1 / ≥2: 452 / 84 / 51)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 139 / 34 / 17 | 0.423 / 0.105 / 0.105 | 0.226 / 0.048 / 0.039 | 0.294 / 0.066 / 0.057 | – |
|  | HPRC v1.1 | 85 / 13 / 7 | 0.448 / 0.235 / 0.222 | 0.153 / 0.048 / 0.039 | 0.228 / 0.079 / 0.067 | -0.067 / +0.014 / +0.010 |
|  | HPRC v2.1, best seed (seed 2) | 75 / 17 / 10 | 0.516 / 0.227 / 0.167 | 0.177 / 0.060 / 0.039 | 0.264 / 0.094 / 0.063 | -0.031 / +0.029 / +0.006 |
|  | HPRC v2.1, mean of 3 seeds | 78.3 / 23.3 / 12.7 | 0.472 / 0.207 / 0.119 | 0.155 / 0.071 / 0.033 | 0.233 / 0.106 / 0.051 | -0.062 / +0.040 / -0.006 |
| PacBio HiFi | GRCh38-only | 325 / 137 / 79 | 0.246 / 0.062 / 0.060 | 0.235 / 0.107 / 0.098 | 0.240 / 0.078 / 0.074 | – |
|  | HPRC v1.1 | 107 / 23 / 11 | 0.473 / 0.258 / 0.083 | 0.212 / 0.095 / 0.020 | 0.293 / 0.139 / 0.032 | +0.053 / +0.061 / -0.042 |
|  | HPRC v2.1, best seed (seed 2) | 85 / 9 / 6 | 0.525 / 0.500 / 0.143 | 0.208 / 0.107 / 0.020 | 0.298 / 0.176 / 0.034 | +0.058 / +0.098 / -0.040 |
|  | HPRC v2.1, mean of 3 seeds | 91.0 / 13.0 / 6.0 | 0.496 / 0.419 / 0.145 | 0.198 / 0.107 / 0.020 | 0.283 / 0.170 / 0.034 | +0.042 / +0.091 / -0.040 |
| ONT | GRCh38-only | 197 / 28 / 9 | 0.140 / 0.034 / 0.100 | 0.071 / 0.012 / 0.020 | 0.094 / 0.018 / 0.033 | – |
|  | HPRC v1.1 | 36 / 6 / 1 | 0.526 / 0.250 / 0.500 | 0.088 / 0.024 / 0.020 | 0.152 / 0.043 / 0.038 | +0.058 / +0.026 / +0.005 |
|  | HPRC v2.1, best seed (seed 0) | 34 / 7 / 2 | 0.564 / 0.222 / 0.333 | 0.097 / 0.024 / 0.020 | 0.166 / 0.043 / 0.037 | +0.072 / +0.025 / +0.004 |
|  | HPRC v2.1, mean of 3 seeds | 32.7 / 5.7 / 2.0 | 0.566 / 0.265 / 0.333 | 0.094 / 0.024 / 0.020 | 0.162 / 0.044 / 0.037 | +0.068 / +0.026 / +0.004 |

## 2. chr1 INDEL, COLO829T, six-tumor model (truth INDELs 0 / 1 / ≥2: 92 / 22 / 17)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 56 / 16 / 7 | 0.211 / 0.111 / 0.125 | 0.163 / 0.091 / 0.059 | 0.184 / 0.100 / 0.080 | – |
|  | HPRC v1.1 | 31 / 8 / 2 | 0.340 / 0.111 / 0.333 | 0.174 / 0.045 / 0.059 | 0.230 / 0.065 / 0.100 | +0.046 / -0.035 / +0.020 |
|  | HPRC v2.1, best seed (seed 2) | 20 / 9 / 2 | 0.412 / 0.182 / 0.333 | 0.152 / 0.091 / 0.059 | 0.222 / 0.121 / 0.100 | +0.038 / +0.021 / +0.020 |
|  | HPRC v2.1, mean of 3 seeds | 30.7 / 19.0 / 6.3 | 0.358 / 0.106 / 0.111 | 0.170 / 0.076 / 0.020 | 0.226 / 0.085 / 0.033 | +0.042 / -0.015 / -0.047 |
| Fiber-seq | GRCh38-only | 144 / 70 / 60 | 0.065 / 0.041 / 0.016 | 0.109 / 0.136 / 0.059 | 0.081 / 0.063 / 0.026 | – |
|  | HPRC v1.1 | 39 / 12 / 7 | 0.235 / 0.077 / 0.125 | 0.130 / 0.045 / 0.059 | 0.168 / 0.057 / 0.080 | +0.087 / -0.006 / +0.054 |
|  | HPRC v2.1, best seed (seed 1) | 26 / 6 / 5 | 0.381 / 0.143 / 0.167 | 0.174 / 0.045 / 0.059 | 0.239 / 0.069 / 0.087 | +0.158 / +0.006 / +0.061 |
|  | HPRC v2.1, mean of 3 seeds | 25.7 / 4.3 / 3.3 | 0.343 / 0.198 / 0.250 | 0.145 / 0.045 / 0.059 | 0.203 / 0.073 / 0.094 | +0.122 / +0.010 / +0.068 |
| ONT | GRCh38-only | 48 / 6 / 1 | 0.020 / 0.000 / 0.000 | 0.011 / 0.000 / 0.000 | 0.014 / 0.000 / 0.000 | – |
|  | HPRC v1.1 | 22 / 6 / 1 | 0.214 / 0.000 / 0.500 | 0.065 / 0.000 / 0.059 | 0.100 / 0.000 / 0.105 | +0.086 / +0.000 / +0.105 |
|  | HPRC v2.1, best seed (seed 0) | 33 / 9 / 3 | 0.233 / 0.000 / 0.250 | 0.109 / 0.000 / 0.059 | 0.148 / 0.000 / 0.095 | +0.134 / +0.000 / +0.095 |
|  | HPRC v2.1, mean of 3 seeds | 26.7 / 6.7 / 2.0 | 0.257 / 0.048 / 0.361 | 0.094 / 0.015 / 0.059 | 0.136 / 0.023 / 0.100 | +0.122 / +0.023 / +0.100 |

## 3. chr1 INDEL, HCC1395T, six-tumor model (truth INDELs 0 / 1 / ≥2: 119 / 12 / 2)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 57 / 28 / 7 | 0.593 / 0.222 / 0.222 | 0.697 / 0.667 / 1.000 | 0.641 / 0.333 / 0.364 | – |
|  | HPRC v1.1 | 65 / 20 / 6 | 0.567 / 0.310 / 0.250 | 0.714 / 0.750 / 1.000 | 0.632 / 0.439 / 0.400 | -0.009 / +0.106 / +0.036 |
|  | HPRC v2.1, best seed (seed 2) | 51 / 8 / 5 | 0.625 / 0.556 / 0.286 | 0.714 / 0.833 / 1.000 | 0.667 / 0.667 / 0.444 | +0.026 / +0.333 / +0.081 |
|  | HPRC v2.1, mean of 3 seeds | 69.0 / 15.7 / 9.0 | 0.559 / 0.417 / 0.198 | 0.714 / 0.833 / 1.000 | 0.626 / 0.547 / 0.327 | -0.015 / +0.214 / -0.037 |
| PacBio HiFi | GRCh38-only | 123 / 51 / 24 | 0.403 / 0.136 / 0.077 | 0.697 / 0.667 / 1.000 | 0.511 / 0.225 / 0.143 | – |
|  | HPRC v1.1 | 65 / 15 / 3 | 0.564 / 0.348 / 0.400 | 0.706 / 0.667 / 1.000 | 0.627 / 0.457 / 0.571 | +0.116 / +0.232 / +0.429 |
|  | HPRC v2.1, best seed (seed 2) | 46 / 8 / 2 | 0.641 / 0.429 / 0.500 | 0.689 / 0.500 / 1.000 | 0.664 / 0.462 / 0.667 | +0.153 / +0.236 / +0.524 |
|  | HPRC v2.1, mean of 3 seeds | 61.0 / 9.7 / 4.0 | 0.578 / 0.410 / 0.361 | 0.692 / 0.556 / 1.000 | 0.629 / 0.470 / 0.522 | +0.118 / +0.245 / +0.379 |
| ONT | GRCh38-only | 135 / 26 / 7 | 0.318 / 0.161 / 0.222 | 0.529 / 0.417 / 1.000 | 0.397 / 0.233 / 0.364 | – |
|  | HPRC v1.1 | 44 / 5 / 1 | 0.600 / 0.500 / 0.667 | 0.555 / 0.417 / 1.000 | 0.576 / 0.455 / 0.800 | +0.179 / +0.222 / +0.436 |
|  | HPRC v2.1, best seed (seed 2) | 50 / 5 / 1 | 0.587 / 0.545 / 0.667 | 0.597 / 0.500 / 1.000 | 0.592 / 0.522 / 0.800 | +0.194 / +0.289 / +0.436 |
|  | HPRC v2.1, mean of 3 seeds | 61.0 / 6.3 / 2.0 | 0.540 / 0.509 / 0.556 | 0.597 / 0.528 / 1.000 | 0.566 / 0.514 / 0.700 | +0.169 / +0.282 / +0.336 |

## 4. chr1 INDEL, HG008T, four-cell-line model (truth INDELs 0 / 1 / ≥2: 452 / 84 / 51)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 41 / 24 / 6 | 0.431 / 0.077 / 0.143 | 0.069 / 0.024 / 0.020 | 0.118 / 0.036 / 0.034 | – |
|  | HPRC v1.1 | 25 / 5 / 2 | 0.603 / 0.167 / 0.333 | 0.084 / 0.012 / 0.020 | 0.148 / 0.022 / 0.037 | +0.029 / -0.014 / +0.003 |
|  | HPRC v2.1, best seed (seed 0) | 42 / 7 / 1 | 0.506 / 0.300 / 0.500 | 0.095 / 0.036 / 0.020 | 0.160 / 0.064 / 0.038 | +0.042 / +0.027 / +0.003 |
|  | HPRC v2.1, mean of 3 seeds | 32.0 / 4.7 / 1.3 | 0.545 / 0.378 / 0.333 | 0.083 / 0.032 / 0.013 | 0.144 / 0.058 / 0.025 | +0.026 / +0.022 / -0.009 |
| PacBio HiFi | GRCh38-only | 133 / 51 / 28 | 0.281 / 0.073 / 0.000 | 0.115 / 0.048 / 0.000 | 0.163 / 0.058 / 0.000 | – |
|  | HPRC v1.1 | 48 / 16 / 7 | 0.515 / 0.200 / 0.125 | 0.113 / 0.048 / 0.020 | 0.185 / 0.077 / 0.034 | +0.022 / +0.019 / +0.034 |
|  | HPRC v2.1, best seed (seed 2) | 40 / 6 / 3 | 0.545 / 0.400 / 0.250 | 0.106 / 0.048 / 0.020 | 0.178 / 0.085 / 0.036 | +0.015 / +0.028 / +0.036 |
|  | HPRC v2.1, mean of 3 seeds | 41.3 / 8.0 / 4.7 | 0.515 / 0.438 / 0.187 | 0.097 / 0.052 / 0.020 | 0.164 / 0.091 / 0.035 | +0.000 / +0.033 / +0.035 |
| ONT | GRCh38-only | 18 / 17 / 5 | 0.379 / 0.000 / 0.000 | 0.024 / 0.000 / 0.000 | 0.046 / 0.000 / 0.000 | – |
|  | HPRC v1.1 | 34 / 14 / 5 | 0.521 / 0.222 / 0.167 | 0.082 / 0.048 / 0.020 | 0.141 / 0.078 / 0.035 | +0.096 / +0.078 / +0.035 |
|  | HPRC v2.1, best seed (seed 2) | 22 / 3 / 3 | 0.600 / 0.400 / 0.250 | 0.073 / 0.024 / 0.020 | 0.130 / 0.045 / 0.036 | +0.084 / +0.045 / +0.036 |
|  | HPRC v2.1, mean of 3 seeds | 20.3 / 6.0 / 1.3 | 0.604 / 0.270 / 0.583 | 0.069 / 0.024 / 0.020 | 0.123 / 0.043 / 0.038 | +0.077 / +0.043 / +0.038 |

## 5. chr1 INDEL, COLO829T, four-cell-line model (truth INDELs 0 / 1 / ≥2: 92 / 22 / 17)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 38 / 5 / 2 | 0.224 / 0.167 / 0.333 | 0.120 / 0.045 / 0.059 | 0.156 / 0.071 / 0.100 | – |
|  | HPRC v1.1 | 11 / 4 / 1 | 0.522 / 0.000 / 0.500 | 0.130 / 0.000 / 0.059 | 0.209 / 0.000 / 0.105 | +0.053 / -0.071 / +0.005 |
|  | HPRC v2.1, best seed (seed 0) | 18 / 3 / 1 | 0.419 / 0.250 / 0.500 | 0.141 / 0.045 / 0.059 | 0.211 / 0.077 / 0.105 | +0.055 / +0.005 / +0.005 |
|  | HPRC v2.1, mean of 3 seeds | 14.0 / 3.3 / 0.3 | 0.459 / 0.233 / 0.833 | 0.127 / 0.045 / 0.059 | 0.198 / 0.076 / 0.109 | +0.042 / +0.005 / +0.009 |
| Fiber-seq | GRCh38-only | 74 / 27 / 5 | 0.140 / 0.036 / 0.000 | 0.130 / 0.045 / 0.000 | 0.135 / 0.040 / 0.000 | – |
|  | HPRC v1.1 | 30 / 9 / 7 | 0.286 / 0.100 / 0.125 | 0.130 / 0.045 / 0.059 | 0.179 / 0.063 / 0.080 | +0.044 / +0.023 / +0.080 |
|  | HPRC v2.1, best seed (seed 0) | 28 / 13 / 9 | 0.317 / 0.071 / 0.100 | 0.141 / 0.045 / 0.059 | 0.195 / 0.056 / 0.074 | +0.061 / +0.016 / +0.074 |
|  | HPRC v2.1, mean of 3 seeds | 44.3 / 10.3 / 5.7 | 0.221 / 0.091 / 0.167 | 0.127 / 0.045 / 0.059 | 0.159 / 0.060 / 0.085 | +0.025 / +0.020 / +0.085 |
| ONT | GRCh38-only | 466 / 73 / 10 | 0.011 / 0.000 / 0.091 | 0.054 / 0.000 / 0.059 | 0.018 / 0.000 / 0.071 | – |
|  | HPRC v1.1 | 134 / 23 / 6 | 0.069 / 0.042 / 0.143 | 0.109 / 0.045 / 0.059 | 0.085 / 0.043 / 0.083 | +0.067 / +0.043 / +0.012 |
|  | HPRC v2.1, best seed (seed 2) | 55 / 10 / 3 | 0.113 / 0.091 / 0.000 | 0.076 / 0.045 / 0.000 | 0.091 / 0.061 / 0.000 | +0.073 / +0.061 / -0.071 |
|  | HPRC v2.1, mean of 3 seeds | 100.3 / 19.0 / 2.3 | 0.078 / 0.058 / 0.222 | 0.083 / 0.045 / 0.039 | 0.078 / 0.049 / 0.067 | +0.060 / +0.049 / -0.005 |

## 6. chr1 INDEL, HCC1395T, four-cell-line model (truth INDELs 0 / 1 / ≥2: 119 / 12 / 2)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 51 / 21 / 4 | 0.622 / 0.300 / 0.333 | 0.706 / 0.750 / 1.000 | 0.661 / 0.429 / 0.500 | – |
|  | HPRC v1.1 | 35 / 9 / 3 | 0.703 / 0.471 / 0.400 | 0.697 / 0.667 / 1.000 | 0.700 / 0.552 / 0.571 | +0.039 / +0.123 / +0.071 |
|  | HPRC v2.1, best seed (seed 1) | 30 / 5 / 1 | 0.732 / 0.643 / 0.667 | 0.689 / 0.750 / 1.000 | 0.710 / 0.692 / 0.800 | +0.049 / +0.264 / +0.300 |
|  | HPRC v2.1, mean of 3 seeds | 32.3 / 5.0 / 2.0 | 0.713 / 0.623 / 0.522 | 0.675 / 0.694 / 1.000 | 0.694 / 0.656 / 0.679 | +0.032 / +0.227 / +0.179 |
| PacBio HiFi | GRCh38-only | 93 / 30 / 13 | 0.459 / 0.189 / 0.071 | 0.664 / 0.583 / 0.500 | 0.543 / 0.286 / 0.125 | – |
|  | HPRC v1.1 | 49 / 13 / 1 | 0.632 / 0.409 / 0.667 | 0.706 / 0.750 / 1.000 | 0.667 / 0.529 / 0.800 | +0.124 / +0.244 / +0.675 |
|  | HPRC v2.1, best seed (seed 1) | 44 / 6 / 3 | 0.642 / 0.500 / 0.400 | 0.664 / 0.500 / 1.000 | 0.653 / 0.500 / 0.571 | +0.110 / +0.214 / +0.446 |
|  | HPRC v2.1, mean of 3 seeds | 44.3 / 8.3 / 2.7 | 0.641 / 0.439 / 0.433 | 0.664 / 0.528 / 1.000 | 0.652 / 0.476 / 0.603 | +0.109 / +0.190 / +0.478 |
| ONT | GRCh38-only | 128 / 29 / 10 | 0.312 / 0.094 / 0.091 | 0.487 / 0.250 / 0.500 | 0.380 / 0.136 / 0.154 | – |
|  | HPRC v1.1 | 76 / 10 / 1 | 0.483 / 0.333 / 0.667 | 0.597 / 0.417 / 1.000 | 0.534 / 0.370 / 0.800 | +0.154 / +0.234 / +0.646 |
|  | HPRC v2.1, best seed (seed 2) | 39 / 6 / 1 | 0.639 / 0.500 / 0.667 | 0.580 / 0.500 / 1.000 | 0.608 / 0.500 / 0.800 | +0.228 / +0.364 / +0.646 |
|  | HPRC v2.1, mean of 3 seeds | 53.0 / 8.3 / 1.3 | 0.563 / 0.397 / 0.556 | 0.557 / 0.444 / 0.833 | 0.559 / 0.418 / 0.656 | +0.179 / +0.282 / +0.502 |

## 7. chr1 SNV, HG008T, six-tumor model (truth SNVs 0 / 1 / ≥2: 561 / 97 / 39)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 354 / 80 / 83 | 0.482 / 0.304 / 0.057 | 0.588 / 0.361 / 0.128 | 0.530 / 0.330 / 0.079 | – |
|  | HPRC v1.1 | 386 / 65 / 27 | 0.549 / 0.511 / 0.289 | 0.836 / 0.701 / 0.282 | 0.662 / 0.591 / 0.286 | +0.132 / +0.261 / +0.207 |
|  | HPRC v2.1 | 369 / 60 / 31 | 0.555 / 0.535 / 0.279 | 0.820 / 0.711 / 0.308 | 0.662 / 0.611 / 0.293 | +0.132 / +0.280 / +0.214 |
| PacBio HiFi | GRCh38-only | 520 / 109 / 51 | 0.465 / 0.377 / 0.164 | 0.806 / 0.680 / 0.256 | 0.590 / 0.485 / 0.200 | – |
|  | HPRC v1.1 | 392 / 68 / 10 | 0.558 / 0.507 / 0.500 | 0.881 / 0.722 / 0.256 | 0.683 / 0.596 / 0.339 | +0.093 / +0.110 / +0.139 |
|  | HPRC v2.1 | 398 / 61 / 20 | 0.551 / 0.523 / 0.355 | 0.872 / 0.691 / 0.282 | 0.675 / 0.596 / 0.314 | +0.086 / +0.110 / +0.114 |
| ONT | GRCh38-only | 353 / 57 / 19 | 0.526 / 0.452 / 0.240 | 0.697 / 0.485 / 0.154 | 0.599 / 0.468 / 0.188 | – |
|  | HPRC v1.1 | 305 / 45 / 9 | 0.587 / 0.563 / 0.526 | 0.774 / 0.598 / 0.256 | 0.668 / 0.580 / 0.345 | +0.068 / +0.112 / +0.157 |
|  | HPRC v2.1 | 217 / 45 / 11 | 0.645 / 0.545 / 0.421 | 0.704 / 0.557 / 0.205 | 0.673 / 0.551 / 0.276 | +0.074 / +0.083 / +0.088 |

## 8. chr1 SNV, COLO829T, six-tumor model (truth SNVs 0 / 1 / ≥2: 2091 / 282 / 70)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 1418 / 45 / 17 | 0.481 / 0.783 / 0.673 | 0.627 / 0.574 / 0.500 | 0.544 / 0.663 / 0.574 | – |
|  | HPRC v1.1 | 543 / 60 / 22 | 0.746 / 0.778 / 0.672 | 0.761 / 0.745 / 0.643 | 0.753 / 0.761 / 0.657 | +0.209 / +0.098 / +0.083 |
|  | HPRC v2.1 | 475 / 44 / 17 | 0.771 / 0.825 / 0.717 | 0.764 / 0.734 / 0.614 | 0.768 / 0.777 / 0.662 | +0.223 / +0.114 / +0.088 |
| Fiber-seq | GRCh38-only | 1961 / 339 / 337 | 0.436 / 0.349 / 0.097 | 0.724 / 0.645 / 0.514 | 0.544 / 0.453 / 0.163 | – |
|  | HPRC v1.1 | 343 / 56 / 12 | 0.833 / 0.786 / 0.793 | 0.815 / 0.730 / 0.657 | 0.824 / 0.757 / 0.719 | +0.280 / +0.304 / +0.556 |
|  | HPRC v2.1 | 249 / 31 / 7 | 0.865 / 0.865 / 0.857 | 0.760 / 0.706 / 0.600 | 0.809 / 0.777 / 0.706 | +0.265 / +0.324 / +0.543 |
| ONT | GRCh38-only | 532 / 49 / 21 | 0.687 / 0.756 / 0.580 | 0.559 / 0.539 / 0.414 | 0.616 / 0.629 / 0.483 | – |
|  | HPRC v1.1 | 519 / 64 / 16 | 0.747 / 0.757 / 0.724 | 0.735 / 0.706 / 0.600 | 0.741 / 0.730 / 0.656 | +0.125 / +0.101 / +0.173 |
|  | HPRC v2.1 | 392 / 58 / 17 | 0.792 / 0.771 / 0.707 | 0.714 / 0.691 / 0.586 | 0.751 / 0.729 / 0.641 | +0.135 / +0.100 / +0.157 |

## 9. chr1 SNV, HCC1395T, six-tumor model (truth SNVs 0 / 1 / ≥2: 3032 / 356 / 52)

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 601 / 57 / 5 | 0.820 / 0.842 / 0.872 | 0.902 / 0.854 / 0.654 | 0.859 / 0.848 / 0.747 | – |
|  | HPRC v1.1 | 455 / 62 / 13 | 0.858 / 0.842 / 0.759 | 0.904 / 0.927 / 0.788 | 0.880 / 0.882 / 0.774 | +0.021 / +0.034 / +0.026 |
|  | HPRC v2.1 | 434 / 52 / 9 | 0.864 / 0.862 / 0.820 | 0.906 / 0.916 / 0.788 | 0.884 / 0.888 / 0.804 | +0.025 / +0.040 / +0.057 |
| PacBio HiFi | GRCh38-only | 539 / 78 / 10 | 0.831 / 0.802 / 0.783 | 0.877 / 0.888 / 0.692 | 0.853 / 0.843 / 0.735 | – |
|  | HPRC v1.1 | 401 / 63 / 4 | 0.865 / 0.829 / 0.895 | 0.845 / 0.857 / 0.654 | 0.855 / 0.843 / 0.756 | +0.001 / -0.000 / +0.021 |
|  | HPRC v2.1 | 353 / 54 / 6 | 0.881 / 0.851 / 0.857 | 0.859 / 0.868 / 0.692 | 0.870 / 0.860 / 0.766 | +0.016 / +0.017 / +0.031 |
| ONT | GRCh38-only | 537 / 55 / 18 | 0.805 / 0.823 / 0.640 | 0.729 / 0.716 / 0.615 | 0.765 / 0.766 / 0.627 | – |
|  | HPRC v1.1 | 356 / 50 / 7 | 0.854 / 0.832 / 0.816 | 0.689 / 0.697 / 0.596 | 0.763 / 0.758 / 0.689 | -0.002 / -0.007 / +0.061 |
|  | HPRC v2.1 | 632 / 85 / 16 | 0.764 / 0.744 / 0.660 | 0.675 / 0.694 / 0.596 | 0.717 / 0.718 / 0.626 | -0.048 / -0.048 / -0.001 |

## 10. Truth variants with an exact candidate per stratum, and the whole-chr1 result

| chr1 call set | Exact candidate (0 / 1 / ≥2) | chr1 all: P / R / F1 |
|---|---|---|
| INDEL, six-tumor model, HG008T Illumina, GRCh38-only | 90.5% / 77.4% / 62.7% | 0.362 / 0.184 / 0.244 |
| INDEL, six-tumor model, HG008T Illumina, HPRC v1.1 | 45.8% / 27.4% / 17.6% | 0.417 / 0.128 / 0.196 |
| INDEL, six-tumor model, HG008T Illumina, HPRC v2.1, seed 0 | 50.0% / 31.0% / 17.6% | 0.376 / 0.111 / 0.171 |
| INDEL, six-tumor model, HG008T Illumina, HPRC v2.1, seed 1 | 50.0% / 31.0% / 17.6% | 0.379 / 0.138 / 0.202 |
| INDEL, six-tumor model, HG008T Illumina, HPRC v2.1, seed 2 (best) | 50.0% / 31.0% / 17.6% | 0.460 / 0.148 / 0.224 |
| INDEL, six-tumor model, HG008T PacBio HiFi, GRCh38-only | 98.2% / 92.9% / 86.3% | 0.182 / 0.204 / 0.192 |
| INDEL, six-tumor model, HG008T PacBio HiFi, HPRC v1.1 | 56.2% / 33.3% / 21.6% | 0.427 / 0.179 / 0.252 |
| INDEL, six-tumor model, HG008T PacBio HiFi, HPRC v2.1, seed 0 | 57.5% / 35.7% / 21.6% | 0.468 / 0.160 / 0.239 |
| INDEL, six-tumor model, HG008T PacBio HiFi, HPRC v2.1, seed 1 | 57.5% / 35.7% / 21.6% | 0.448 / 0.170 / 0.247 |
| INDEL, six-tumor model, HG008T PacBio HiFi, HPRC v2.1, seed 2 (best) | 57.5% / 35.7% / 21.6% | 0.510 / 0.177 / 0.263 |
| INDEL, six-tumor model, HG008T ONT, GRCh38-only | 88.7% / 73.8% / 72.5% | 0.127 / 0.058 / 0.080 |
| INDEL, six-tumor model, HG008T ONT, HPRC v1.1 | 58.4% / 35.7% / 19.6% | 0.500 / 0.073 / 0.128 |
| INDEL, six-tumor model, HG008T ONT, HPRC v2.1, seed 0 (best) | 58.2% / 38.1% / 21.6% | 0.522 / 0.080 / 0.139 |
| INDEL, six-tumor model, HG008T ONT, HPRC v2.1, seed 1 | 58.2% / 38.1% / 21.6% | 0.524 / 0.075 / 0.131 |
| INDEL, six-tumor model, HG008T ONT, HPRC v2.1, seed 2 | 58.2% / 38.1% / 21.6% | 0.548 / 0.078 / 0.137 |
| INDEL, six-tumor model, COLO829T Illumina, GRCh38-only | 90.2% / 81.8% / 70.6% | 0.186 / 0.137 / 0.158 |
| INDEL, six-tumor model, COLO829T Illumina, HPRC v1.1 | 54.3% / 27.3% / 5.9% | 0.305 / 0.137 / 0.189 |
| INDEL, six-tumor model, COLO829T Illumina, HPRC v2.1, seed 0 | 58.7% / 27.3% / 5.9% | 0.262 / 0.130 / 0.173 |
| INDEL, six-tumor model, COLO829T Illumina, HPRC v2.1, seed 1 | 58.7% / 27.3% / 5.9% | 0.176 / 0.145 / 0.159 |
| INDEL, six-tumor model, COLO829T Illumina, HPRC v2.1, seed 2 (best) | 58.7% / 27.3% / 5.9% | 0.354 / 0.130 / 0.190 |
| INDEL, six-tumor model, COLO829T Fiber-seq, GRCh38-only | 94.6% / 100.0% / 94.1% | 0.049 / 0.107 / 0.067 |
| INDEL, six-tumor model, COLO829T Fiber-seq, HPRC v1.1 | 63.0% / 31.8% / 11.8% | 0.194 / 0.107 / 0.138 |
| INDEL, six-tumor model, COLO829T Fiber-seq, HPRC v2.1, seed 0 | 68.5% / 31.8% / 11.8% | 0.350 / 0.107 / 0.164 |
| INDEL, six-tumor model, COLO829T Fiber-seq, HPRC v2.1, seed 1 (best) | 68.5% / 31.8% / 11.8% | 0.327 / 0.137 / 0.194 |
| INDEL, six-tumor model, COLO829T Fiber-seq, HPRC v2.1, seed 2 | 68.5% / 31.8% / 11.8% | 0.275 / 0.107 / 0.154 |
| INDEL, six-tumor model, COLO829T ONT, GRCh38-only | 98.9% / 86.4% / 88.2% | 0.018 / 0.008 / 0.011 |
| INDEL, six-tumor model, COLO829T ONT, HPRC v1.1 | 65.2% / 27.3% / 17.6% | 0.194 / 0.053 / 0.084 |
| INDEL, six-tumor model, COLO829T ONT, HPRC v2.1, seed 0 (best) | 67.4% / 27.3% / 17.6% | 0.196 / 0.084 / 0.118 |
| INDEL, six-tumor model, COLO829T ONT, HPRC v2.1, seed 1 | 67.4% / 27.3% / 17.6% | 0.204 / 0.076 / 0.111 |
| INDEL, six-tumor model, COLO829T ONT, HPRC v2.1, seed 2 | 67.4% / 27.3% / 17.6% | 0.290 / 0.069 / 0.111 |
| INDEL, six-tumor model, HCC1395T Illumina, GRCh38-only | 90.8% / 91.7% / 100.0% | 0.503 / 0.699 / 0.585 |
| INDEL, six-tumor model, HCC1395T Illumina, HPRC v1.1 | 84.9% / 83.3% / 100.0% | 0.513 / 0.722 / 0.600 |
| INDEL, six-tumor model, HCC1395T Illumina, HPRC v2.1, seed 0 | 87.4% / 83.3% / 100.0% | 0.527 / 0.744 / 0.617 |
| INDEL, six-tumor model, HCC1395T Illumina, HPRC v2.1, seed 1 | 87.4% / 83.3% / 100.0% | 0.426 / 0.714 / 0.534 |
| INDEL, six-tumor model, HCC1395T Illumina, HPRC v2.1, seed 2 (best) | 87.4% / 83.3% / 100.0% | 0.602 / 0.729 / 0.660 |
| INDEL, six-tumor model, HCC1395T PacBio HiFi, GRCh38-only | 89.9% / 91.7% / 100.0% | 0.320 / 0.699 / 0.439 |
| INDEL, six-tumor model, HCC1395T PacBio HiFi, HPRC v1.1 | 84.0% / 83.3% / 100.0% | 0.531 / 0.707 / 0.606 |
| INDEL, six-tumor model, HCC1395T PacBio HiFi, HPRC v2.1, seed 0 | 86.6% / 75.0% / 100.0% | 0.523 / 0.684 / 0.593 |
| INDEL, six-tumor model, HCC1395T PacBio HiFi, HPRC v2.1, seed 1 | 86.6% / 75.0% / 100.0% | 0.520 / 0.692 / 0.594 |
| INDEL, six-tumor model, HCC1395T PacBio HiFi, HPRC v2.1, seed 2 (best) | 86.6% / 75.0% / 100.0% | 0.616 / 0.677 / 0.645 |
| INDEL, six-tumor model, HCC1395T ONT, GRCh38-only | 81.5% / 83.3% / 100.0% | 0.294 / 0.526 / 0.377 |
| INDEL, six-tumor model, HCC1395T ONT, HPRC v1.1 | 78.2% / 75.0% / 100.0% | 0.593 / 0.549 / 0.570 |
| INDEL, six-tumor model, HCC1395T ONT, HPRC v2.1, seed 0 | 79.0% / 66.7% / 100.0% | 0.491 / 0.602 / 0.541 |
| INDEL, six-tumor model, HCC1395T ONT, HPRC v2.1, seed 1 | 79.0% / 66.7% / 100.0% | 0.534 / 0.594 / 0.562 |
| INDEL, six-tumor model, HCC1395T ONT, HPRC v2.1, seed 2 (best) | 79.0% / 66.7% / 100.0% | 0.585 / 0.594 / 0.590 |
| INDEL, four-cell-line model, HG008T Illumina, GRCh38-only | 90.5% / 77.4% / 62.7% | 0.324 / 0.058 / 0.098 |
| INDEL, four-cell-line model, HG008T Illumina, HPRC v1.1 | 45.8% / 27.4% / 17.6% | 0.556 / 0.068 / 0.121 |
| INDEL, four-cell-line model, HG008T Illumina, HPRC v2.1, seed 0 (best) | 50.0% / 31.0% / 17.6% | 0.485 / 0.080 / 0.137 |
| INDEL, four-cell-line model, HG008T Illumina, HPRC v2.1, seed 1 | 50.0% / 31.0% / 17.6% | 0.552 / 0.063 / 0.113 |
| INDEL, four-cell-line model, HG008T Illumina, HPRC v2.1, seed 2 | 50.0% / 31.0% / 17.6% | 0.534 / 0.066 / 0.118 |
| INDEL, four-cell-line model, HG008T PacBio HiFi, GRCh38-only | 98.2% / 92.9% / 86.3% | 0.209 / 0.095 / 0.131 |
| INDEL, four-cell-line model, HG008T PacBio HiFi, HPRC v1.1 | 56.2% / 33.3% / 21.6% | 0.441 / 0.095 / 0.157 |
| INDEL, four-cell-line model, HG008T PacBio HiFi, HPRC v2.1, seed 0 | 57.5% / 35.7% / 21.6% | 0.431 / 0.080 / 0.135 |
| INDEL, four-cell-line model, HG008T PacBio HiFi, HPRC v2.1, seed 1 | 57.5% / 35.7% / 21.6% | 0.485 / 0.082 / 0.140 |
| INDEL, four-cell-line model, HG008T PacBio HiFi, HPRC v2.1, seed 2 (best) | 57.5% / 35.7% / 21.6% | 0.520 / 0.090 / 0.154 |
| INDEL, four-cell-line model, HG008T ONT, GRCh38-only | 88.7% / 73.8% / 72.5% | 0.216 / 0.019 / 0.034 |
| INDEL, four-cell-line model, HG008T ONT, HPRC v1.1 | 58.4% / 35.7% / 19.6% | 0.442 / 0.072 / 0.123 |
| INDEL, four-cell-line model, HG008T ONT, HPRC v2.1, seed 0 | 58.2% / 38.1% / 21.6% | 0.574 / 0.060 / 0.108 |
| INDEL, four-cell-line model, HG008T ONT, HPRC v2.1, seed 1 | 58.2% / 38.1% / 21.6% | 0.517 / 0.053 / 0.096 |
| INDEL, four-cell-line model, HG008T ONT, HPRC v2.1, seed 2 (best) | 58.2% / 38.1% / 21.6% | 0.562 / 0.061 / 0.111 |
| INDEL, four-cell-line model, COLO829T Illumina, GRCh38-only | 90.2% / 81.8% / 70.6% | 0.224 / 0.099 / 0.138 |
| INDEL, four-cell-line model, COLO829T Illumina, HPRC v1.1 | 54.3% / 27.3% / 5.9% | 0.448 / 0.099 / 0.162 |
| INDEL, four-cell-line model, COLO829T Illumina, HPRC v2.1, seed 0 (best) | 58.7% / 27.3% / 5.9% | 0.405 / 0.115 / 0.179 |
| INDEL, four-cell-line model, COLO829T Illumina, HPRC v2.1, seed 1 | 58.7% / 27.3% / 5.9% | 0.464 / 0.099 / 0.164 |
| INDEL, four-cell-line model, COLO829T Illumina, HPRC v2.1, seed 2 | 58.7% / 27.3% / 5.9% | 0.448 / 0.099 / 0.162 |
| INDEL, four-cell-line model, COLO829T Fiber-seq, GRCh38-only | 94.6% / 100.0% / 94.1% | 0.109 / 0.099 / 0.104 |
| INDEL, four-cell-line model, COLO829T Fiber-seq, HPRC v1.1 | 63.0% / 31.8% / 11.8% | 0.233 / 0.107 / 0.147 |
| INDEL, four-cell-line model, COLO829T Fiber-seq, HPRC v2.1, seed 0 (best) | 68.5% / 31.8% / 11.8% | 0.231 / 0.115 / 0.153 |
| INDEL, four-cell-line model, COLO829T Fiber-seq, HPRC v2.1, seed 1 | 68.5% / 31.8% / 11.8% | 0.154 / 0.092 / 0.115 |
| INDEL, four-cell-line model, COLO829T Fiber-seq, HPRC v2.1, seed 2 | 68.5% / 31.8% / 11.8% | 0.177 / 0.107 / 0.133 |
| INDEL, four-cell-line model, COLO829T ONT, GRCh38-only | 98.9% / 86.4% / 88.2% | 0.011 / 0.046 / 0.017 |
| INDEL, four-cell-line model, COLO829T ONT, HPRC v1.1 | 65.2% / 27.3% / 17.6% | 0.069 / 0.092 / 0.078 |
| INDEL, four-cell-line model, COLO829T ONT, HPRC v2.1, seed 0 | 67.4% / 27.3% / 17.6% | 0.071 / 0.084 / 0.077 |
| INDEL, four-cell-line model, COLO829T ONT, HPRC v2.1, seed 1 | 67.4% / 27.3% / 17.6% | 0.056 / 0.069 / 0.061 |
| INDEL, four-cell-line model, COLO829T ONT, HPRC v2.1, seed 2 (best) | 67.4% / 27.3% / 17.6% | 0.105 / 0.061 / 0.077 |
| INDEL, four-cell-line model, HCC1395T Illumina, GRCh38-only | 90.8% / 91.7% / 100.0% | 0.556 / 0.714 / 0.625 |
| INDEL, four-cell-line model, HCC1395T Illumina, HPRC v1.1 | 84.9% / 83.3% / 100.0% | 0.664 / 0.699 / 0.681 |
| INDEL, four-cell-line model, HCC1395T Illumina, HPRC v2.1, seed 0 | 87.4% / 83.3% / 100.0% | 0.659 / 0.654 / 0.657 |
| INDEL, four-cell-line model, HCC1395T Illumina, HPRC v2.1, seed 1 (best) | 87.4% / 83.3% / 100.0% | 0.721 / 0.699 / 0.710 |
| INDEL, four-cell-line model, HCC1395T Illumina, HPRC v2.1, seed 2 | 87.4% / 83.3% / 100.0% | 0.713 / 0.692 / 0.702 |
| INDEL, four-cell-line model, HCC1395T PacBio HiFi, GRCh38-only | 89.9% / 91.7% / 100.0% | 0.390 / 0.654 / 0.489 |
| INDEL, four-cell-line model, HCC1395T PacBio HiFi, HPRC v1.1 | 84.0% / 83.3% / 100.0% | 0.601 / 0.714 / 0.653 |
| INDEL, four-cell-line model, HCC1395T PacBio HiFi, HPRC v2.1, seed 0 | 86.6% / 75.0% / 100.0% | 0.601 / 0.669 / 0.633 |
| INDEL, four-cell-line model, HCC1395T PacBio HiFi, HPRC v2.1, seed 1 (best) | 86.6% / 75.0% / 100.0% | 0.621 / 0.654 / 0.637 |
| INDEL, four-cell-line model, HCC1395T PacBio HiFi, HPRC v2.1, seed 2 | 86.6% / 75.0% / 100.0% | 0.614 / 0.647 / 0.630 |
| INDEL, four-cell-line model, HCC1395T ONT, GRCh38-only | 81.5% / 83.3% / 100.0% | 0.271 / 0.466 / 0.343 |
| INDEL, four-cell-line model, HCC1395T ONT, HPRC v1.1 | 78.2% / 75.0% / 100.0% | 0.473 / 0.586 / 0.523 |
| INDEL, four-cell-line model, HCC1395T ONT, HPRC v2.1, seed 0 | 79.0% / 66.7% / 100.0% | 0.558 / 0.541 / 0.550 |
| INDEL, four-cell-line model, HCC1395T ONT, HPRC v2.1, seed 1 | 79.0% / 66.7% / 100.0% | 0.455 / 0.534 / 0.491 |
| INDEL, four-cell-line model, HCC1395T ONT, HPRC v2.1, seed 2 (best) | 79.0% / 66.7% / 100.0% | 0.626 / 0.579 / 0.602 |
| SNV, six-tumor model, HG008T Illumina, GRCh38-only | 98.9% / 97.9% / 87.2% | 0.417 / 0.531 / 0.467 |
| SNV, six-tumor model, HG008T Illumina, HPRC v1.1 | 98.6% / 92.8% / 43.6% | 0.534 / 0.786 / 0.636 |
| SNV, six-tumor model, HG008T Illumina, HPRC v2.1 | 98.2% / 94.8% / 59.0% | 0.540 / 0.776 / 0.637 |
| SNV, six-tumor model, HG008T PacBio HiFi, GRCh38-only | 99.5% / 97.9% / 76.9% | 0.437 / 0.758 / 0.554 |
| SNV, six-tumor model, HG008T PacBio HiFi, HPRC v1.1 | 98.4% / 92.8% / 43.6% | 0.550 / 0.824 / 0.659 |
| SNV, six-tumor model, HG008T PacBio HiFi, HPRC v2.1 | 98.0% / 94.8% / 51.3% | 0.542 / 0.813 / 0.651 |
| SNV, six-tumor model, HG008T ONT, GRCh38-only | 99.6% / 97.9% / 69.2% | 0.509 / 0.637 / 0.566 |
| SNV, six-tumor model, HG008T ONT, HPRC v1.1 | 98.9% / 92.8% / 48.7% | 0.583 / 0.720 / 0.644 |
| SNV, six-tumor model, HG008T ONT, HPRC v2.1 | 98.9% / 94.8% / 48.7% | 0.626 / 0.656 / 0.641 |
| SNV, six-tumor model, COLO829T Illumina, GRCh38-only | 88.4% / 90.1% / 91.4% | 0.505 / 0.618 / 0.556 |
| SNV, six-tumor model, COLO829T Illumina, HPRC v1.1 | 87.5% / 89.0% / 81.4% | 0.747 / 0.756 / 0.751 |
| SNV, six-tumor model, COLO829T Illumina, HPRC v2.1 | 87.3% / 88.3% / 82.9% | 0.775 / 0.756 / 0.766 |
| SNV, six-tumor model, COLO829T Fiber-seq, GRCh38-only | 86.9% / 85.8% / 84.3% | 0.396 / 0.709 / 0.509 |
| SNV, six-tumor model, COLO829T Fiber-seq, HPRC v1.1 | 86.5% / 85.8% / 81.4% | 0.826 / 0.801 / 0.814 |
| SNV, six-tumor model, COLO829T Fiber-seq, HPRC v2.1 | 86.4% / 85.5% / 84.3% | 0.864 / 0.749 / 0.803 |
| SNV, six-tumor model, COLO829T ONT, GRCh38-only | 82.4% / 82.3% / 77.1% | 0.691 / 0.552 / 0.614 |
| SNV, six-tumor model, COLO829T ONT, HPRC v1.1 | 82.0% / 82.3% / 74.3% | 0.748 / 0.727 / 0.737 |
| SNV, six-tumor model, COLO829T ONT, HPRC v2.1 | 82.1% / 81.6% / 78.6% | 0.787 / 0.708 / 0.745 |
| SNV, six-tumor model, HCC1395T Illumina, GRCh38-only | 95.8% / 98.9% / 94.2% | 0.823 / 0.894 / 0.857 |
| SNV, six-tumor model, HCC1395T Illumina, HPRC v1.1 | 95.5% / 98.0% / 94.2% | 0.855 / 0.905 / 0.879 |
| SNV, six-tumor model, HCC1395T Illumina, HPRC v2.1 | 95.5% / 98.0% / 92.3% | 0.863 / 0.905 / 0.883 |
| SNV, six-tumor model, HCC1395T PacBio HiFi, GRCh38-only | 92.9% / 96.3% / 90.4% | 0.828 / 0.875 / 0.851 |
| SNV, six-tumor model, HCC1395T PacBio HiFi, HPRC v1.1 | 92.6% / 95.5% / 90.4% | 0.861 / 0.843 / 0.852 |
| SNV, six-tumor model, HCC1395T PacBio HiFi, HPRC v2.1 | 92.5% / 95.2% / 88.5% | 0.877 / 0.858 / 0.867 |
| SNV, six-tumor model, HCC1395T ONT, GRCh38-only | 77.0% / 79.5% / 80.8% | 0.804 / 0.726 / 0.763 |
| SNV, six-tumor model, HCC1395T ONT, HPRC v1.1 | 76.9% / 78.7% / 80.8% | 0.851 / 0.688 / 0.761 |
| SNV, six-tumor model, HCC1395T ONT, HPRC v2.1 | 76.8% / 78.9% / 78.8% | 0.760 / 0.676 / 0.716 |

## 11. Paired recall in stratum ≥2, HPRC vs GRCh38-only (same truth variants)

| chr1 comparison | truth ≥2 | found by HPRC only / by GRCh38-only only | exact McNemar p |
|---|---|---|---|
| INDEL, six-tumor model, HG008T Illumina, HPRC v1.1 | 51 | 1 / 1 | 1.000 |
| INDEL, six-tumor model, HG008T Illumina, HPRC v2.1 (best seed, seed 2) | 51 | 0 / 0 | - |
| INDEL, six-tumor model, HG008T PacBio HiFi, HPRC v1.1 | 51 | 0 / 4 | 0.125 |
| INDEL, six-tumor model, HG008T PacBio HiFi, HPRC v2.1 (best seed, seed 2) | 51 | 0 / 4 | 0.125 |
| INDEL, six-tumor model, HG008T ONT, HPRC v1.1 | 51 | 0 / 0 | - |
| INDEL, six-tumor model, HG008T ONT, HPRC v2.1 (best seed, seed 0) | 51 | 0 / 0 | - |
| INDEL, six-tumor model, COLO829T Illumina, HPRC v1.1 | 17 | 0 / 0 | - |
| INDEL, six-tumor model, COLO829T Illumina, HPRC v2.1 (best seed, seed 2) | 17 | 0 / 0 | - |
| INDEL, six-tumor model, COLO829T Fiber-seq, HPRC v1.1 | 17 | 1 / 1 | 1.000 |
| INDEL, six-tumor model, COLO829T Fiber-seq, HPRC v2.1 (best seed, seed 1) | 17 | 1 / 1 | 1.000 |
| INDEL, six-tumor model, COLO829T ONT, HPRC v1.1 | 17 | 1 / 0 | 1.000 |
| INDEL, six-tumor model, COLO829T ONT, HPRC v2.1 (best seed, seed 0) | 17 | 1 / 0 | 1.000 |
| INDEL, six-tumor model, HCC1395T Illumina, HPRC v1.1 | 2 | 0 / 0 | - |
| INDEL, six-tumor model, HCC1395T Illumina, HPRC v2.1 (best seed, seed 2) | 2 | 0 / 0 | - |
| INDEL, six-tumor model, HCC1395T PacBio HiFi, HPRC v1.1 | 2 | 0 / 0 | - |
| INDEL, six-tumor model, HCC1395T PacBio HiFi, HPRC v2.1 (best seed, seed 2) | 2 | 0 / 0 | - |
| INDEL, six-tumor model, HCC1395T ONT, HPRC v1.1 | 2 | 0 / 0 | - |
| INDEL, six-tumor model, HCC1395T ONT, HPRC v2.1 (best seed, seed 2) | 2 | 0 / 0 | - |
| INDEL, four-cell-line model, HG008T Illumina, HPRC v1.1 | 51 | 0 / 0 | - |
| INDEL, four-cell-line model, HG008T Illumina, HPRC v2.1 (best seed, seed 0) | 51 | 0 / 0 | - |
| INDEL, four-cell-line model, HG008T PacBio HiFi, HPRC v1.1 | 51 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, HG008T PacBio HiFi, HPRC v2.1 (best seed, seed 2) | 51 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, HG008T ONT, HPRC v1.1 | 51 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, HG008T ONT, HPRC v2.1 (best seed, seed 2) | 51 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, COLO829T Illumina, HPRC v1.1 | 17 | 0 / 0 | - |
| INDEL, four-cell-line model, COLO829T Illumina, HPRC v2.1 (best seed, seed 0) | 17 | 0 / 0 | - |
| INDEL, four-cell-line model, COLO829T Fiber-seq, HPRC v1.1 | 17 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, COLO829T Fiber-seq, HPRC v2.1 (best seed, seed 0) | 17 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, COLO829T ONT, HPRC v1.1 | 17 | 1 / 1 | 1.000 |
| INDEL, four-cell-line model, COLO829T ONT, HPRC v2.1 (best seed, seed 2) | 17 | 0 / 1 | 1.000 |
| INDEL, four-cell-line model, HCC1395T Illumina, HPRC v1.1 | 2 | 0 / 0 | - |
| INDEL, four-cell-line model, HCC1395T Illumina, HPRC v2.1 (best seed, seed 1) | 2 | 0 / 0 | - |
| INDEL, four-cell-line model, HCC1395T PacBio HiFi, HPRC v1.1 | 2 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, HCC1395T PacBio HiFi, HPRC v2.1 (best seed, seed 1) | 2 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, HCC1395T ONT, HPRC v1.1 | 2 | 1 / 0 | 1.000 |
| INDEL, four-cell-line model, HCC1395T ONT, HPRC v2.1 (best seed, seed 2) | 2 | 1 / 0 | 1.000 |
| SNV, six-tumor model, HG008T Illumina, HPRC v1.1 | 39 | 6 / 0 | 0.031 |
| SNV, six-tumor model, HG008T Illumina, HPRC v2.1 | 39 | 7 / 0 | 0.016 |
| SNV, six-tumor model, HG008T PacBio HiFi, HPRC v1.1 | 39 | 1 / 1 | 1.000 |
| SNV, six-tumor model, HG008T PacBio HiFi, HPRC v2.1 | 39 | 1 / 0 | 1.000 |
| SNV, six-tumor model, HG008T ONT, HPRC v1.1 | 39 | 4 / 0 | 0.125 |
| SNV, six-tumor model, HG008T ONT, HPRC v2.1 | 39 | 3 / 1 | 0.625 |
| SNV, six-tumor model, COLO829T Illumina, HPRC v1.1 | 70 | 12 / 2 | 0.013 |
| SNV, six-tumor model, COLO829T Illumina, HPRC v2.1 | 70 | 13 / 5 | 0.096 |
| SNV, six-tumor model, COLO829T Fiber-seq, HPRC v1.1 | 70 | 14 / 4 | 0.031 |
| SNV, six-tumor model, COLO829T Fiber-seq, HPRC v2.1 | 70 | 8 / 2 | 0.109 |
| SNV, six-tumor model, COLO829T ONT, HPRC v1.1 | 70 | 14 / 1 | <0.001 |
| SNV, six-tumor model, COLO829T ONT, HPRC v2.1 | 70 | 14 / 2 | 0.004 |
| SNV, six-tumor model, HCC1395T Illumina, HPRC v1.1 | 52 | 9 / 2 | 0.065 |
| SNV, six-tumor model, HCC1395T Illumina, HPRC v2.1 | 52 | 10 / 3 | 0.092 |
| SNV, six-tumor model, HCC1395T PacBio HiFi, HPRC v1.1 | 52 | 2 / 4 | 0.688 |
| SNV, six-tumor model, HCC1395T PacBio HiFi, HPRC v2.1 | 52 | 2 / 2 | 1.000 |
| SNV, six-tumor model, HCC1395T ONT, HPRC v1.1 | 52 | 5 / 6 | 1.000 |
| SNV, six-tumor model, HCC1395T ONT, HPRC v2.1 | 52 | 3 / 4 | 1.000 |

## Pending call sets

none
