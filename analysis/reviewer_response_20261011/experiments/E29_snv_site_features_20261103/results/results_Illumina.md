## Illumina: SNV F1 (P / R) after the PoN, benchmark BED, PASS at each caller's own rule

| test | Pansoma + site features (ppr9) | Pansoma (normfix) | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|
| HCC1395T chr20 | 0.821 (0.862 / 0.783) | 0.809 (0.839 / 0.781) | 0.845 (0.769 / 0.937) | 0.730 (0.649 / 0.835) |
| HCC1395T chr1 | 0.852 (0.901 / 0.807) | 0.850 (0.879 / 0.822) | 0.848 (0.781 / 0.927) | 0.832 (0.806 / 0.859) |
| HG008T chr1 | 0.644 (0.604 / 0.690) | 0.601 (0.592 / 0.611) | 0.529 (0.376 / 0.889) | 0.552 (0.410 / 0.845) |
| COLO829T chr1 | 0.732 (0.839 / 0.650) | 0.755 (0.880 / 0.661) | 0.565 (0.427 / 0.835) | 0.551 (0.431 / 0.764) |

PASS threshold HCC1395T chr20: site features 0.867, normfix 0.735
PASS threshold HCC1395T chr1: site features 0.867, normfix 0.735
PASS threshold HG008T chr1: site features 0.867, normfix 0.735
PASS threshold COLO829T chr1: site features 0.867, normfix 0.735

## Illumina: HCC1395T chr20 recall by read AF (truth SNVs with a Pansoma tensor in the VCF; after the PoN, BED)

| AF | truth | site features | normfix | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|---|
| [0.08, 0.25] | 231 | 0.939 | 0.965 | 0.987 | 0.835 |
| [0.25, 0.4] | 222 | 0.923 | 0.937 | 0.986 | 0.914 |
| [0.4, 0.6] | 433 | 0.836 | 0.811 | 0.970 | 0.940 |
| [0.6, 0.85] | 98 | 0.776 | 0.796 | 0.959 | 0.918 |
| [0.85, 1.0] | 9 | 0.556 | 0.444 | 0.889 | 0.556 |
| no tensor in the VCF | 113 | 0.009 | 0.000 | 0.593 | 0.221 |

## Illumina: HCC1395T chr1 recall by read AF (truth SNVs with a Pansoma tensor in the VCF; after the PoN, BED)

| AF | truth | site features | normfix | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|---|
| [0.08, 0.25] | 796 | 0.931 | 0.979 | 0.985 | 0.835 |
| [0.25, 0.4] | 896 | 0.922 | 0.949 | 0.978 | 0.953 |
| [0.4, 0.6] | 749 | 0.861 | 0.896 | 0.960 | 0.940 |
| [0.6, 0.85] | 453 | 0.848 | 0.819 | 0.951 | 0.967 |
| [0.85, 1.0] | 242 | 0.740 | 0.645 | 0.930 | 0.971 |
| no tensor in the VCF | 304 | 0.007 | 0.000 | 0.510 | 0.197 |

## Rule 2 (registered 2026-11-03 21:57, before this result) and the validation best-F1 cutoff: at_cutoffs.py

Illumina epoch 12: validation cutoffs recall 0.9 (ppr9) p >= 0.867, recall 0.95 p >= 0.777, best F1 p >= 0.765
| test | recall 0.9 (ppr9) | recall 0.95 | best F1 |
|---|---|---|---|
| HCC1395T chr20 | 0.820 (0.862 / 0.783) | 0.830 (0.831 / 0.829) | 0.830 (0.827 / 0.833) |
| HCC1395T chr1 | 0.852 (0.901 / 0.807) | 0.864 (0.878 / 0.851) | 0.865 (0.876 / 0.855) |
| HG008T chr1 | 0.644 (0.604 / 0.690) | 0.652 (0.557 / 0.786) | 0.650 (0.550 / 0.795) |
| COLO829T chr1 | 0.732 (0.839 / 0.650) | 0.760 (0.814 / 0.713) | 0.759 (0.808 / 0.716) |

HCC1395T chr20 recall by read AF

| AF | truth | recall 0.9 (ppr9) | recall 0.95 | best F1 | normfix PASS |
|---|---|---|---|---|---|
| [0.08, 0.25] | 231 | 0.939 | 0.952 | 0.952 | 0.965 |
| [0.25, 0.4] | 222 | 0.923 | 0.959 | 0.959 | 0.937 |
| [0.4, 0.6] | 433 | 0.836 | 0.905 | 0.915 | 0.811 |
| [0.6, 0.85] | 98 | 0.776 | 0.867 | 0.867 | 0.796 |
| [0.85, 1.0] | 9 | 0.556 | 0.667 | 0.667 | 0.444 |

HCC1395T chr1 recall by read AF

| AF | truth | recall 0.9 (ppr9) | recall 0.95 | best F1 | normfix PASS |
|---|---|---|---|---|---|
| [0.08, 0.25] | 796 | 0.931 | 0.956 | 0.961 | 0.979 |
| [0.25, 0.4] | 896 | 0.922 | 0.953 | 0.955 | 0.949 |
| [0.4, 0.6] | 749 | 0.861 | 0.923 | 0.928 | 0.896 |
| [0.6, 0.85] | 453 | 0.848 | 0.912 | 0.918 | 0.819 |
| [0.85, 1.0] | 242 | 0.740 | 0.847 | 0.860 | 0.645 |

## Diagnostic: matched recall (matched_recall.py)

## Illumina HCC1395T chr20: recall 0.781 (normfix PASS, p >= 0.735) vs site features at p >= 0.870 (recall 0.781)
precision there: normfix 0.839, site features 0.864

| AF | truth | site features, matched recall | normfix PASS |
|---|---|---|---|
| [0.08, 0.25] | 231 | 0.939 | 0.965 |
| [0.25, 0.4] | 222 | 0.923 | 0.937 |
| [0.4, 0.6] | 433 | 0.831 | 0.811 |
| [0.6, 0.85] | 98 | 0.776 | 0.796 |
| [0.85, 1.0] | 9 | 0.556 | 0.444 |

## Illumina HCC1395T chr1: recall 0.822 (normfix PASS, p >= 0.735) vs site features at p >= 0.842 (recall 0.822)
precision there: normfix 0.879, site features 0.892

| AF | truth | site features, matched recall | normfix PASS |
|---|---|---|---|
| [0.08, 0.25] | 796 | 0.942 | 0.979 |
| [0.25, 0.4] | 896 | 0.932 | 0.949 |
| [0.4, 0.6] | 749 | 0.893 | 0.896 |
| [0.6, 0.85] | 453 | 0.854 | 0.819 |
| [0.85, 1.0] | 242 | 0.764 | 0.645 |

