## PacBio: SNV F1 (P / R) after the PoN, benchmark BED, PASS at each caller's own rule

| test | Pansoma + site features (ppr9) | Pansoma (normfix) | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|
| HCC1395T chr20 | 0.794 (0.847 / 0.747) | 0.788 (0.804 / 0.773) | 0.835 (0.740 / 0.957) | 0.795 (0.778 / 0.812) |
| HCC1395T chr1 | 0.828 (0.895 / 0.769) | 0.825 (0.843 / 0.809) | 0.839 (0.766 / 0.927) | 0.814 (0.805 / 0.822) |
| HG008T chr1 | 0.659 (0.614 / 0.710) | 0.625 (0.586 / 0.670) | 0.557 (0.425 / 0.806) | 0.591 (0.455 / 0.844) |
| COLO829T chr1 | 0.777 (0.879 / 0.697) | 0.774 (0.880 / 0.691) | 0.632 (0.732 / 0.555) | 0.762 (0.779 / 0.745) |

PASS threshold HCC1395T chr20: site features 0.875, normfix 0.622
PASS threshold HCC1395T chr1: site features 0.875, normfix 0.622
PASS threshold HG008T chr1: site features 0.875, normfix 0.622
PASS threshold COLO829T chr1: site features 0.875, normfix 0.622

## PacBio: HCC1395T chr20 recall by read AF (truth SNVs with a Pansoma tensor in the VCF; after the PoN, BED)

| AF | truth | site features | normfix | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|---|
| [0.08, 0.25] | 193 | 0.927 | 0.969 | 0.990 | 0.912 |
| [0.25, 0.4] | 226 | 0.903 | 0.925 | 0.973 | 0.907 |
| [0.4, 0.6] | 404 | 0.879 | 0.894 | 0.978 | 0.943 |
| [0.6, 0.85] | 125 | 0.672 | 0.752 | 0.952 | 0.912 |
| [0.85, 1.0] | 7 | 0.571 | 0.571 | 0.857 | 0.571 |
| no tensor in the VCF | 151 | 0.000 | 0.000 | 0.841 | 0.119 |

## PacBio: HCC1395T chr1 recall by read AF (truth SNVs with a Pansoma tensor in the VCF; after the PoN, BED)

| AF | truth | site features | normfix | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|---|
| [0.08, 0.25] | 673 | 0.944 | 0.961 | 0.979 | 0.854 |
| [0.25, 0.4] | 775 | 0.920 | 0.952 | 0.970 | 0.957 |
| [0.4, 0.6] | 881 | 0.877 | 0.952 | 0.966 | 0.949 |
| [0.6, 0.85] | 427 | 0.763 | 0.836 | 0.963 | 0.939 |
| [0.85, 1.0] | 251 | 0.789 | 0.801 | 0.888 | 0.944 |
| no tensor in the VCF | 433 | 0.005 | 0.000 | 0.679 | 0.088 |

## Diagnostic (after the results; not pre-specified): the same comparison at matched recall

## PacBio HCC1395T chr20: recall 0.773 (normfix PASS, p >= 0.622) vs site features at p >= 0.859 (recall 0.773)
precision there: normfix 0.804, site features 0.842

| AF | truth | site features, matched recall | normfix PASS |
|---|---|---|---|
| [0.08, 0.25] | 193 | 0.943 | 0.969 |
| [0.25, 0.4] | 226 | 0.925 | 0.925 |
| [0.4, 0.6] | 404 | 0.908 | 0.894 |
| [0.6, 0.85] | 125 | 0.744 | 0.752 |
| [0.85, 1.0] | 7 | 0.571 | 0.571 |

## PacBio HCC1395T chr1: recall 0.809 (normfix PASS, p >= 0.622) vs site features at p >= 0.840 (recall 0.809)
precision there: normfix 0.843, site features 0.879

| AF | truth | site features, matched recall | normfix PASS |
|---|---|---|---|
| [0.08, 0.25] | 673 | 0.961 | 0.961 |
| [0.25, 0.4] | 775 | 0.937 | 0.952 |
| [0.4, 0.6] | 881 | 0.935 | 0.952 |
| [0.6, 0.85] | 427 | 0.845 | 0.836 |
| [0.85, 1.0] | 251 | 0.884 | 0.801 |


## Cutoffs from the checkpoint's validation (`at_cutoffs.py`; after the results)

PacBio epoch 6: validation cutoffs recall 0.9 (ppr9) p >= 0.875, recall 0.95 p >= 0.825, best F1 p >= 0.795
| test | recall 0.9 (ppr9) | recall 0.95 | best F1 |
|---|---|---|---|
| HCC1395T chr20 | 0.794 (0.847 / 0.747) | 0.813 (0.832 / 0.795) | 0.815 (0.823 / 0.807) |
| HCC1395T chr1 | 0.828 (0.895 / 0.769) | 0.845 (0.872 / 0.820) | 0.847 (0.863 / 0.831) |
| HG008T chr1 | 0.659 (0.614 / 0.710) | 0.645 (0.549 / 0.780) | 0.650 (0.538 / 0.821) |
| COLO829T chr1 | 0.777 (0.879 / 0.697) | 0.792 (0.859 / 0.734) | 0.793 (0.840 / 0.752) |

HCC1395T chr20 recall by read AF

| AF | truth | recall 0.9 (ppr9) | recall 0.95 | best F1 | normfix PASS |
|---|---|---|---|---|---|
| [0.08, 0.25] | 193 | 0.927 | 0.953 | 0.964 | 0.969 |
| [0.25, 0.4] | 226 | 0.903 | 0.929 | 0.938 | 0.925 |
| [0.4, 0.6] | 404 | 0.879 | 0.938 | 0.943 | 0.894 |
| [0.6, 0.85] | 125 | 0.672 | 0.816 | 0.880 | 0.752 |
| [0.85, 1.0] | 7 | 0.571 | 0.571 | 0.571 | 0.571 |

HCC1395T chr1 recall by read AF

| AF | truth | recall 0.9 (ppr9) | recall 0.95 | best F1 | normfix PASS |
|---|---|---|---|---|---|
| [0.08, 0.25] | 673 | 0.944 | 0.964 | 0.969 | 0.961 |
| [0.25, 0.4] | 775 | 0.920 | 0.946 | 0.954 | 0.952 |
| [0.4, 0.6] | 881 | 0.877 | 0.948 | 0.959 | 0.952 |
| [0.6, 0.85] | 427 | 0.763 | 0.874 | 0.906 | 0.836 |
| [0.85, 1.0] | 251 | 0.789 | 0.908 | 0.932 | 0.801 |
