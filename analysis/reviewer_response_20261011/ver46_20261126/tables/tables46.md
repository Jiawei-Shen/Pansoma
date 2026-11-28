## t12. Response Tables 1.2 = 2.10 = 3.4: chr1 INDELs before the population filter, P / R / F1 (seed 0, no site features in any model)

| chr1 INDELs, no population filter (P / R / F1) | Pansoma GRCh38-only | Pansoma HPRC v1.1 | Pansoma HPRC v2.1 |
|---|---|---|---|
| HG008T Illumina | 0.016 / 0.445 / 0.031 | 0.125 / 0.208 / 0.156 | 0.140 / 0.182 / 0.158 |
| HG008T PacBio HiFi | 0.015 / 0.475 / 0.029 | 0.111 / 0.278 / 0.159 | 0.139 / 0.239 / 0.176 |
| HG008T ONT | 0.015 / 0.097 / 0.027 | 0.101 / 0.082 / 0.090 | 0.121 / 0.083 / 0.099 |
| COLO829T Illumina | 0.006 / 0.298 / 0.013 | 0.082 / 0.153 / 0.106 | 0.064 / 0.160 / 0.091 |
| COLO829T Fiber-seq | 0.015 / 0.198 / 0.027 | 0.051 / 0.122 / 0.072 | 0.070 / 0.153 / 0.096 |
| COLO829T ONT | 0.006 / 0.008 / 0.006 | 0.039 / 0.046 / 0.042 | 0.026 / 0.038 / 0.031 |
| HCC1395T Illumina | 0.015 / 0.744 / 0.029 | 0.136 / 0.744 / 0.230 | 0.184 / 0.767 / 0.297 |
| HCC1395T PacBio HiFi | 0.010 / 0.797 / 0.020 | 0.117 / 0.729 / 0.201 | 0.212 / 0.714 / 0.326 |
| HCC1395T ONT | 0.019 / 0.534 / 0.036 | 0.186 / 0.534 / 0.276 | 0.201 / 0.571 / 0.297 |

- HPRC v1.1 F1 > GRCh38-only: 9 of 9
- HPRC v2.1 F1 > GRCh38-only: 9 of 9
- HPRC v2.1 F1 > HPRC v1.1: 7 of 9
- HPRC v2.1 F1 > both comparators (3 decimals): 6 of 9
- HPRC v1.1 F1 > both comparators (3 decimals): 6 of 9
- HG008T F1 range GRCh38-only: 0.027-0.031
- HG008T F1 range HPRC v1.1: 0.090-0.159
- HG008T F1 range HPRC v2.1: 0.099-0.176
- COLO829T F1 range GRCh38-only: 0.006-0.027
- COLO829T F1 range HPRC v1.1: 0.042-0.106
- COLO829T F1 range HPRC v2.1: 0.031-0.096
- HCC1395T F1 range GRCh38-only: 0.020-0.036
- HCC1395T F1 range HPRC v1.1: 0.201-0.276
- HCC1395T F1 range HPRC v2.1: 0.297-0.326

## t32. Response Table 3.2, SNV rows: six-tumor SNV models on HPRC v2.1, seeds 0 / 1 / 2 (mean ± SD), chr1 after the population filter

| tumor | platform | training set | precision | recall | F1 |
|---|---|---|---|---|---|
| HG008T | Illumina | 6 tumors, SNV | 0.540 / 0.576 / 0.567 (0.561 ± 0.018) | 0.776 / 0.783 / 0.782 (0.780 ± 0.004) | 0.637 / 0.664 / 0.657 (0.653 ± 0.014) |
| HG008T | PacBio HiFi | 6 tumors, SNV | 0.542 / 0.536 / 0.534 (0.537 ± 0.004) | 0.813 / 0.829 / 0.799 (0.814 ± 0.015) | 0.651 / 0.651 / 0.640 (0.647 ± 0.006) |
| HG008T | ONT | 6 tumors, SNV | 0.626 / 0.617 / 0.609 (0.617 ± 0.009) | 0.656 / 0.736 / 0.730 (0.707 ± 0.045) | 0.641 / 0.671 / 0.664 (0.659 ± 0.016) |
| COLO829T | Illumina | 6 tumors, SNV | 0.775 / 0.840 / 0.836 (0.817 ± 0.036) | 0.756 / 0.718 / 0.747 (0.740 ± 0.020) | 0.766 / 0.774 / 0.789 (0.776 ± 0.012) |
| COLO829T | Fiber-seq | 6 tumors, SNV | 0.864 / 0.845 / 0.840 (0.850 ± 0.013) | 0.749 / 0.770 / 0.772 (0.764 ± 0.013) | 0.803 / 0.806 / 0.805 (0.804 ± 0.002) |
| COLO829T | ONT | 6 tumors, SNV | 0.787 / 0.790 / 0.754 (0.777 ± 0.020) | 0.708 / 0.722 / 0.726 (0.719 ± 0.010) | 0.745 / 0.755 / 0.740 (0.747 ± 0.008) |
| HCC1395T | Illumina | 6 tumors, SNV | 0.863 / 0.874 / 0.870 (0.869 ± 0.006) | 0.905 / 0.897 / 0.891 (0.897 ± 0.007) | 0.883 / 0.885 / 0.880 (0.883 ± 0.003) |
| HCC1395T | PacBio HiFi | 6 tumors, SNV | 0.877 / 0.873 / 0.870 (0.874 ± 0.003) | 0.858 / 0.856 / 0.863 (0.859 ± 0.004) | 0.867 / 0.864 / 0.866 (0.866 ± 0.002) |
| HCC1395T | ONT | 6 tumors, SNV | 0.760 / 0.767 / 0.611 (0.713 ± 0.088) | 0.676 / 0.684 / 0.674 (0.678 ± 0.005) | 0.716 / 0.723 / 0.641 (0.693 ± 0.045) |

- SD of F1 over the three seeds: 0.002–0.045 (median 0.008); largest: HCC1395T ONT 0.045 (0.716 / 0.723 / 0.641); HG008T ONT 0.016 (0.641 / 0.671 / 0.664); HG008T Illumina 0.014 (0.637 / 0.664 / 0.657)

