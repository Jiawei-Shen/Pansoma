# Letter ver4: tables from existing reports of the current models (make_tables.py)

## nofilter_snv: chr1 SNVs, no population filter, six-tumor models on three references (P / R / F1)

Seed 0. Comparator columns: the same call sets as in ver3, P / R / F1 recomputed from the counts.

| chr1 SNVs, no population filter (P / R / F1) | Pansoma GRCh38-only | Pansoma HPRC v1.1 | Pansoma HPRC v2.1 | ClairS-TO | DeepSomatic-TO |
|---|---|---|---|---|---|
| HG008T Illumina | 0.013 / 0.547 / 0.026 | 0.525 / 0.789 / 0.630 | 0.537 / 0.778 / 0.635 | 0.003 / 0.902 / 0.006 | 0.029 / 0.935 / 0.056 |
| HG008T PacBio HiFi | 0.006 / 0.786 / 0.013 | 0.550 / 0.825 / 0.660 | 0.542 / 0.816 / 0.651 | 0.003 / 0.894 / 0.006 | 0.059 / 0.829 / 0.109 |
| HG008T ONT | 0.009 / 0.663 / 0.018 | 0.583 / 0.720 / 0.644 | 0.626 / 0.656 / 0.641 | 0.003 / 0.884 / 0.006 | 0.047 / 0.910 / 0.090 |
| COLO829T Illumina | 0.050 / 0.636 / 0.093 | 0.744 / 0.759 / 0.752 | 0.771 / 0.760 / 0.766 | 0.010 / 0.796 / 0.020 | 0.097 / 0.858 / 0.174 |
| COLO829T Fiber-seq | 0.044 / 0.730 / 0.083 | 0.826 / 0.802 / 0.814 | 0.863 / 0.750 / 0.802 | 0.011 / 0.771 / 0.022 | 0.255 / 0.569 / 0.352 |
| COLO829T ONT | 0.044 / 0.567 / 0.082 | 0.747 / 0.728 / 0.737 | 0.784 / 0.709 / 0.744 | 0.010 / 0.709 / 0.020 | 0.152 / 0.714 / 0.251 |
| HCC1395T Illumina | 0.039 / 0.909 / 0.074 | 0.851 / 0.908 / 0.879 | 0.862 / 0.908 / 0.884 | 0.014 / 0.884 / 0.028 | 0.139 / 0.944 / 0.242 |
| HCC1395T PacBio HiFi | 0.034 / 0.890 / 0.066 | 0.861 / 0.844 / 0.852 | 0.877 / 0.860 / 0.868 | 0.014 / 0.845 / 0.028 | 0.257 / 0.944 / 0.404 |
| HCC1395T ONT | 0.022 / 0.740 / 0.042 | 0.850 / 0.690 / 0.762 | 0.749 / 0.678 / 0.712 | 0.013 / 0.752 / 0.025 | 0.217 / 0.781 / 0.340 |

## nofilter_indel: chr1 INDELs, no population filter, six-tumor models on three references (P / R / F1)

GRCh38-only and HPRC v1.1: seed 0; HPRC v2.1: mean of three seeds. PASS = ppr8 (ONT: pp).

| chr1 INDELs, no population filter (P / R / F1) | Pansoma GRCh38-only | Pansoma HPRC v1.1 | Pansoma HPRC v2.1 | ClairS-TO | DeepSomatic-TO |
|---|---|---|---|---|---|
| HG008T Illumina | 0.016 / 0.445 / 0.031 | 0.412 / 0.128 / 0.195 | 0.322 / 0.137 / 0.179 | 0.004 / 0.101 / 0.007 | 0.017 / 0.206 / 0.032 |
| HG008T PacBio HiFi | 0.015 / 0.475 / 0.029 | 0.398 / 0.189 / 0.256 | 0.464 / 0.177 / 0.256 | 0.004 / 0.092 / 0.007 | 0.080 / 0.037 / 0.051 |
| HG008T ONT | 0.015 / 0.097 / 0.027 | 0.478 / 0.073 / 0.127 | 0.525 / 0.078 / 0.136 | 0.004 / 0.087 / 0.007 | 0.043 / 0.019 / 0.026 |
| COLO829T Illumina | 0.006 / 0.298 / 0.013 | 0.300 / 0.137 / 0.188 | 0.227 / 0.145 / 0.154 | 0.002 / 0.191 / 0.003 | 0.010 / 0.244 / 0.019 |
| COLO829T Fiber-seq | 0.015 / 0.198 / 0.027 | 0.195 / 0.115 / 0.144 | 0.311 / 0.117 / 0.169 | 0.002 / 0.160 / 0.004 | 0.124 / 0.092 / 0.105 |
| COLO829T ONT | 0.006 / 0.008 / 0.006 | 0.189 / 0.053 / 0.083 | 0.228 / 0.076 / 0.113 | 0.002 / 0.137 / 0.003 | 0.016 / 0.015 / 0.016 |
| HCC1395T Illumina | 0.015 / 0.744 / 0.029 | 0.466 / 0.729 / 0.569 | 0.417 / 0.744 / 0.498 | 0.006 / 0.699 / 0.012 | 0.022 / 0.835 / 0.042 |
| HCC1395T PacBio HiFi | 0.010 / 0.797 / 0.020 | 0.508 / 0.714 / 0.594 | 0.515 / 0.699 / 0.589 | 0.005 / 0.526 / 0.011 | 0.271 / 0.647 / 0.382 |
| HCC1395T ONT | 0.019 / 0.534 / 0.036 | 0.570 / 0.549 / 0.559 | 0.526 / 0.596 / 0.558 | 0.006 / 0.564 / 0.012 | 0.249 / 0.481 / 0.328 |

## t31_three_refs: chr1 SNV after the population filter, six-tumor models on three references

Seed 0. GRCh38-only models take no site features. 'no site features' HPRC rows: RIS models, results.md. The paired 95 % intervals of ΔF1 (block bootstrap) are added by the build step from the regions item for the rows whose call sets are on this cluster.

| tumor | platform | reference | precision | recall | F1 | ΔF1, HPRC − GRCh38-only |
|---|---|---|---|---|---|---|
| HG008T | Illumina | GRCh38-only | 0.417 | 0.531 | 0.467 | – |
| HG008T | Illumina | HPRC v1.1 | 0.534 | 0.786 | 0.636 | +0.169 |
| HG008T | Illumina | HPRC v1.1, no site features | 0.523 | 0.706 | 0.601 | +0.134 |
| HG008T | Illumina | HPRC v2.1 | 0.540 | 0.776 | 0.637 | +0.170 |
| HG008T | Illumina | HPRC v2.1, no site features | 0.552 | 0.723 | 0.626 | +0.159 |
| HG008T | PacBio HiFi | GRCh38-only | 0.437 | 0.758 | 0.554 | – |
| HG008T | PacBio HiFi | HPRC v1.1 | 0.550 | 0.824 | 0.659 | +0.105 |
| HG008T | PacBio HiFi | HPRC v1.1, no site features | 0.530 | 0.762 | 0.625 | +0.071 |
| HG008T | PacBio HiFi | HPRC v2.1 | 0.542 | 0.813 | 0.651 | +0.096 |
| HG008T | PacBio HiFi | HPRC v2.1, no site features | 0.521 | 0.769 | 0.621 | +0.067 |
| HG008T | ONT | GRCh38-only | 0.509 | 0.637 | 0.566 | – |
| HG008T | ONT | HPRC v1.1 | 0.583 | 0.720 | 0.644 | +0.079 |
| HG008T | ONT | HPRC v1.1, no site features | 0.594 | 0.613 | 0.603 | +0.038 |
| HG008T | ONT | HPRC v2.1 | 0.626 | 0.656 | 0.641 | +0.075 |
| HG008T | ONT | HPRC v2.1, no site features | 0.621 | 0.620 | 0.620 | +0.055 |
| COLO829T | Illumina | GRCh38-only | 0.505 | 0.618 | 0.556 | – |
| COLO829T | Illumina | HPRC v1.1 | 0.747 | 0.756 | 0.751 | +0.196 |
| COLO829T | Illumina | HPRC v1.1, no site features | 0.766 | 0.752 | 0.759 | +0.203 |
| COLO829T | Illumina | HPRC v2.1 | 0.775 | 0.756 | 0.766 | +0.210 |
| COLO829T | Illumina | HPRC v2.1, no site features | 0.825 | 0.756 | 0.789 | +0.233 |
| COLO829T | Fiber-seq | GRCh38-only | 0.396 | 0.709 | 0.509 | – |
| COLO829T | Fiber-seq | HPRC v1.1 | 0.826 | 0.801 | 0.814 | +0.305 |
| COLO829T | Fiber-seq | HPRC v1.1, no site features | 0.847 | 0.765 | 0.804 | +0.295 |
| COLO829T | Fiber-seq | HPRC v2.1 | 0.864 | 0.749 | 0.803 | +0.294 |
| COLO829T | Fiber-seq | HPRC v2.1, no site features | 0.849 | 0.758 | 0.801 | +0.293 |
| COLO829T | ONT | GRCh38-only | 0.691 | 0.552 | 0.614 | – |
| COLO829T | ONT | HPRC v1.1 | 0.748 | 0.727 | 0.737 | +0.123 |
| COLO829T | ONT | HPRC v1.1, no site features | 0.754 | 0.723 | 0.738 | +0.124 |
| COLO829T | ONT | HPRC v2.1 | 0.787 | 0.708 | 0.745 | +0.131 |
| COLO829T | ONT | HPRC v2.1, no site features | 0.742 | 0.740 | 0.741 | +0.127 |
| HCC1395T | Illumina | GRCh38-only | 0.823 | 0.894 | 0.857 | – |
| HCC1395T | Illumina | HPRC v1.1 | 0.855 | 0.905 | 0.879 | +0.022 |
| HCC1395T | Illumina | HPRC v1.1, no site features | 0.841 | 0.889 | 0.864 | +0.008 |
| HCC1395T | Illumina | HPRC v2.1 | 0.863 | 0.905 | 0.883 | +0.027 |
| HCC1395T | Illumina | HPRC v2.1, no site features | 0.828 | 0.894 | 0.859 | +0.003 |
| HCC1395T | PacBio HiFi | GRCh38-only | 0.828 | 0.875 | 0.851 | – |
| HCC1395T | PacBio HiFi | HPRC v1.1 | 0.861 | 0.843 | 0.852 | +0.001 |
| HCC1395T | PacBio HiFi | HPRC v1.1, no site features | 0.839 | 0.864 | 0.851 | +0.001 |
| HCC1395T | PacBio HiFi | HPRC v2.1 | 0.877 | 0.858 | 0.867 | +0.017 |
| HCC1395T | PacBio HiFi | HPRC v2.1, no site features | 0.831 | 0.871 | 0.851 | 0.000 |
| HCC1395T | ONT | GRCh38-only | 0.804 | 0.726 | 0.763 | – |
| HCC1395T | ONT | HPRC v1.1 | 0.851 | 0.688 | 0.761 | −0.002 |
| HCC1395T | ONT | HPRC v1.1, no site features | 0.815 | 0.710 | 0.759 | −0.004 |
| HCC1395T | ONT | HPRC v2.1 | 0.760 | 0.676 | 0.716 | −0.047 |
| HCC1395T | ONT | HPRC v2.1, no site features | 0.663 | 0.708 | 0.684 | −0.079 |

## t32_indel_seeds: Per-seed chr1 INDEL results of the final models (ONT at pp)

PASS = ppr8 on Illumina / PacBio HiFi / Fiber-seq, pp on ONT. Seeds 0 / 1 / 2 (mean ± sample SD).

| tumor | platform | training set | precision | recall | F1 |
|---|---|---|---|---|---|
| HG008T | Illumina | 6 tumors | 0.376 / 0.379 / 0.460 (0.405 ± 0.048) | 0.111 / 0.138 / 0.148 (0.132 ± 0.019) | 0.171 / 0.202 / 0.224 (0.199 ± 0.027) |
| HG008T | Illumina | 4 tumors | 0.485 / 0.552 / 0.534 (0.524 ± 0.035) | 0.080 / 0.063 / 0.066 (0.070 ± 0.009) | 0.137 / 0.113 / 0.118 (0.123 ± 0.013) |
| HG008T | PacBio HiFi | 6 tumors | 0.468 / 0.448 / 0.510 (0.475 ± 0.031) | 0.160 / 0.170 / 0.177 (0.169 ± 0.009) | 0.239 / 0.247 / 0.263 (0.249 ± 0.012) |
| HG008T | PacBio HiFi | 4 tumors | 0.431 / 0.485 / 0.520 (0.479 ± 0.045) | 0.080 / 0.082 / 0.090 (0.084 ± 0.005) | 0.135 / 0.140 / 0.154 (0.143 ± 0.010) |
| HG008T | ONT | 6 tumors | 0.522 / 0.524 / 0.548 (0.531 ± 0.014) | 0.080 / 0.075 / 0.078 (0.078 ± 0.003) | 0.139 / 0.131 / 0.137 (0.136 ± 0.004) |
| HG008T | ONT | 4 tumors | 0.574 / 0.517 / 0.562 (0.551 ± 0.030) | 0.060 / 0.053 / 0.061 (0.058 ± 0.005) | 0.108 / 0.096 / 0.111 (0.105 ± 0.008) |
| COLO829T | Illumina | 6 tumors | 0.262 / 0.176 / 0.354 (0.264 ± 0.089) | 0.130 / 0.145 / 0.130 (0.135 ± 0.009) | 0.173 / 0.159 / 0.190 (0.174 ± 0.015) |
| COLO829T | Illumina | 4 tumors | 0.405 / 0.464 / 0.448 (0.439 ± 0.030) | 0.115 / 0.099 / 0.099 (0.104 ± 0.009) | 0.179 / 0.164 / 0.162 (0.168 ± 0.009) |
| COLO829T | Fiber-seq | 6 tumors | 0.350 / 0.327 / 0.275 (0.317 ± 0.039) | 0.107 / 0.137 / 0.107 (0.117 ± 0.018) | 0.164 / 0.194 / 0.154 (0.170 ± 0.021) |
| COLO829T | Fiber-seq | 4 tumors | 0.231 / 0.154 / 0.177 (0.187 ± 0.039) | 0.115 / 0.092 / 0.107 (0.104 ± 0.012) | 0.153 / 0.115 / 0.133 (0.134 ± 0.019) |
| COLO829T | ONT | 6 tumors | 0.196 / 0.204 / 0.290 (0.230 ± 0.052) | 0.084 / 0.076 / 0.069 (0.076 ± 0.008) | 0.118 / 0.111 / 0.111 (0.113 ± 0.004) |
| COLO829T | ONT | 4 tumors | 0.071 / 0.056 / 0.105 (0.077 ± 0.025) | 0.084 / 0.069 / 0.061 (0.071 ± 0.012) | 0.077 / 0.061 / 0.077 (0.072 ± 0.009) |
| HCC1395T | Illumina | 6 tumors | 0.527 / 0.426 / 0.602 (0.518 ± 0.089) | 0.744 / 0.714 / 0.729 (0.729 ± 0.015) | 0.617 / 0.534 / 0.660 (0.603 ± 0.064) |
| HCC1395T | Illumina | 4 tumors | 0.659 / 0.721 / 0.713 (0.698 ± 0.034) | 0.654 / 0.699 / 0.692 (0.682 ± 0.024) | 0.657 / 0.710 / 0.702 (0.690 ± 0.029) |
| HCC1395T | PacBio HiFi | 6 tumors | 0.523 / 0.520 / 0.616 (0.553 ± 0.055) | 0.684 / 0.692 / 0.677 (0.684 ± 0.008) | 0.593 / 0.594 / 0.645 (0.611 ± 0.030) |
| HCC1395T | PacBio HiFi | 4 tumors | 0.601 / 0.621 / 0.614 (0.612 ± 0.010) | 0.669 / 0.654 / 0.647 (0.657 ± 0.011) | 0.633 / 0.637 / 0.630 (0.634 ± 0.004) |
| HCC1395T | ONT | 6 tumors | 0.491 / 0.534 / 0.585 (0.537 ± 0.047) | 0.602 / 0.594 / 0.594 (0.596 ± 0.004) | 0.541 / 0.562 / 0.590 (0.564 ± 0.025) |
| HCC1395T | ONT | 4 tumors | 0.558 / 0.455 / 0.626 (0.546 ± 0.086) | 0.541 / 0.534 / 0.579 (0.551 ± 0.024) | 0.550 / 0.491 / 0.602 (0.548 ± 0.055) |

## t32_snv_seeds: Per-seed chr1 SNV results of the COLO829T single-tumor models of the current pipeline

PASS = rule 2. Seeds 0 / 1 / 2 (mean ± sample SD).

| tumor | platform | training set | precision | recall | F1 |
|---|---|---|---|---|---|
| COLO829T | Illumina | COLO829T only, SNV, HPRC v2.1 | 0.802 / 0.665 / 0.828 (0.765 ± 0.088) | 0.791 / 0.799 / 0.805 (0.798 ± 0.007) | 0.797 / 0.725 / 0.816 (0.779 ± 0.048) |
| COLO829T | Illumina | COLO829T only, SNV, HPRC v1.1 | 0.749 / 0.648 / 0.677 (0.691 ± 0.052) | 0.800 / 0.794 / 0.773 (0.789 ± 0.014) | 0.774 / 0.714 / 0.722 (0.736 ± 0.032) |
| COLO829T | Fiber-seq | COLO829T only, SNV, HPRC v2.1 | 0.852 / 0.856 / 0.853 (0.854 ± 0.002) | 0.786 / 0.786 / 0.794 (0.789 ± 0.005) | 0.818 / 0.820 / 0.823 (0.820 ± 0.003) |
| COLO829T | Fiber-seq | COLO829T only, SNV, HPRC v1.1 | 0.851 / 0.857 / 0.872 (0.860 ± 0.011) | 0.788 / 0.784 / 0.782 (0.785 ± 0.003) | 0.819 / 0.819 / 0.825 (0.821 ± 0.003) |

## t27_snv_nearby: chr1 SNV, HG008T, six-tumor model, by nearby germline variants

Copied from E19_germline_density/strata_latest.md section 7.

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

## t28_indel_nearby: chr1 INDEL, HG008T, six-tumor model, by nearby germline variants

Copied from E19_germline_density/strata_latest.md section 1.

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 139 / 34 / 17 | 0.423 / 0.105 / 0.105 | 0.226 / 0.048 / 0.039 | 0.294 / 0.066 / 0.057 | – |
|  | HPRC v1.1 | 85 / 13 / 7 | 0.448 / 0.235 / 0.222 | 0.153 / 0.048 / 0.039 | 0.228 / 0.079 / 0.067 | −0.067 / +0.014 / +0.010 |
|  | HPRC v2.1, best seed (seed 2) | 75 / 17 / 10 | 0.516 / 0.227 / 0.167 | 0.177 / 0.060 / 0.039 | 0.264 / 0.094 / 0.063 | −0.031 / +0.029 / +0.006 |
|  | HPRC v2.1, mean of 3 seeds | 78.3 / 23.3 / 12.7 | 0.472 / 0.207 / 0.119 | 0.155 / 0.071 / 0.033 | 0.233 / 0.106 / 0.051 | −0.062 / +0.040 / −0.006 |
| PacBio HiFi | GRCh38-only | 325 / 137 / 79 | 0.246 / 0.062 / 0.060 | 0.235 / 0.107 / 0.098 | 0.240 / 0.078 / 0.074 | – |
|  | HPRC v1.1 | 107 / 23 / 11 | 0.473 / 0.258 / 0.083 | 0.212 / 0.095 / 0.020 | 0.293 / 0.139 / 0.032 | +0.053 / +0.061 / −0.042 |
|  | HPRC v2.1, best seed (seed 2) | 85 / 9 / 6 | 0.525 / 0.500 / 0.143 | 0.208 / 0.107 / 0.020 | 0.298 / 0.176 / 0.034 | +0.058 / +0.098 / −0.040 |
|  | HPRC v2.1, mean of 3 seeds | 91.0 / 13.0 / 6.0 | 0.496 / 0.419 / 0.145 | 0.198 / 0.107 / 0.020 | 0.283 / 0.170 / 0.034 | +0.042 / +0.091 / −0.040 |
| ONT | GRCh38-only | 197 / 28 / 9 | 0.140 / 0.034 / 0.100 | 0.071 / 0.012 / 0.020 | 0.094 / 0.018 / 0.033 | – |
|  | HPRC v1.1 | 36 / 6 / 1 | 0.526 / 0.250 / 0.500 | 0.088 / 0.024 / 0.020 | 0.152 / 0.043 / 0.038 | +0.058 / +0.026 / +0.005 |
|  | HPRC v2.1, best seed (seed 0) | 34 / 7 / 2 | 0.564 / 0.222 / 0.333 | 0.097 / 0.024 / 0.020 | 0.166 / 0.043 / 0.037 | +0.072 / +0.025 / +0.004 |
|  | HPRC v2.1, mean of 3 seeds | 32.7 / 5.7 / 2.0 | 0.566 / 0.265 / 0.333 | 0.094 / 0.024 / 0.020 | 0.162 / 0.044 / 0.037 | +0.068 / +0.026 / +0.004 |

## absorbed_v21: Sequence context of the somatic truth INDELs absorbed by the HPRC v2.1 graph (chr1–22)

Copied from graph_absorbed_somatic_hcc1395t_20261112/tables_v21_patient_frame.md section 1.

| Truth INDELs (chr1–22) | HG008T | COLO829T | HCC1395T |
|---|---|---|---|
| Absorbed by the graph (union across platforms) | 1,904 (22.4%) | 394 (20.6%) | 46 (2.8%) |
| Homopolymer ≥ 7 bp | 1,508 (79.2%) | 223 (56.6%) | 24 (52.2%) |
| Short tandem repeat (unit 2–6 bp, ≥ 3 copies) | 268 (14.1%) | 65 (16.5%) | 3 (6.5%) |
| Homopolymer + short tandem repeat, total | 1,776 (93.3%) | 288 (73.1%) | 27 (58.7%) |
| Other repeats | 7 (0.4%) | 18 (4.6%) | 1 (2.2%) |
| Substitution within repetitive sequence | 100 (5.3%) | 10 (2.5%) | 2 (4.3%) |
| Others | 19 (1.0%) | 68 (17.3%) | 13 (28.3%) |
| Non-repetitive sequence | 2 (0.1%) | 10 (2.5%) | 3 (6.5%) |

## facts (numbers quoted in sentences)

```
{
 "snv_raw_f1_range": {
  "linear_nositef": "0.013–0.093",
  "d9": "0.630–0.879",
  "d46": "0.635–0.884"
 },
 "snv_raw_f1_v21_by_tumor": {
  "HG008T": "0.635–0.651",
  "COLO829T": "0.744–0.802",
  "HCC1395T": "0.712–0.884"
 },
 "t31_gain_hg_colo_with_sitef": "0.075–0.305",
 "t31_gain_hg_colo_no_sitef": "0.038–0.295",
 "t31_hcc_ont_v21_no_sitef": -0.079,
 "t31_hcc_delta": {
  "Illumina|d9": 0.022,
  "Illumina|d46": 0.027,
  "PacBio|d9": 0.001,
  "PacBio|d46": 0.017,
  "ONT|d9": -0.002,
  "ONT|d46": -0.047
 },
 "snv_v21_vs_comparators": {
  "HG008T Illumina": {
   "four": 0.623,
   "six": 0.637,
   "deepsomatic": 0.529,
   "clairs": 0.553
  },
  "HG008T PacBio": {
   "four": 0.626,
   "six": 0.651,
   "deepsomatic": 0.557,
   "clairs": 0.591
  },
  "HG008T ONT": {
   "four": 0.548,
   "six": 0.641,
   "deepsomatic": 0.457,
   "clairs": 0.552
  },
  "COLO829T Illumina": {
   "four": 0.736,
   "six": 0.766,
   "deepsomatic": 0.565,
   "clairs": 0.551
  },
  "COLO829T PacBio": {
   "four": 0.639,
   "six": 0.803,
   "deepsomatic": 0.632,
   "clairs": 0.762
  },
  "COLO829T ONT": {
   "four": 0.365,
   "six": 0.745,
   "deepsomatic": 0.633,
   "clairs": 0.715
  },
  "HCC1395T Illumina": {
   "four": 0.879,
   "six": 0.883,
   "deepsomatic": 0.848,
   "clairs": 0.832
  },
  "HCC1395T PacBio": {
   "four": 0.862,
   "six": 0.867,
   "deepsomatic": 0.839,
   "clairs": 0.814
  },
  "HCC1395T ONT": {
   "four": 0.443,
   "six": 0.716,
   "deepsomatic": 0.76,
   "clairs": 0.76
  }
 },
 "indel_six_above_both_hg_colo": [
  "HG008T PacBio HiFi",
  "HG008T ONT",
  "COLO829T Fiber-seq",
  "COLO829T ONT"
 ],
 "indel_sd_18": {
  "range": "0.004–0.064",
  "median": 0.014
 },
 "snv_sd_4": {
  "range": "0.003–0.048",
  "median": 0.018
 },
 "t27_delta": {
  "none": "0.068–0.132",
  "one_or_more": "0.083–0.280"
 }
}
```

