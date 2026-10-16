# Lost / absorbed truth INDELs called after the PoN (chr1, in the BED; rtg baseline TPs)

curve = any score (the ceiling), PASS = at the run's recall-0.9 threshold. "pending" = the run's rtg output is not there yet.

| sample | platform | set | n | not PoN-tagged | linear TP curve / PASS | d9 TP curve / PASS | d46 TP curve / PASS |
|---|---|---|---:|---:|---|---|---|
| HG008T | Illumina | all truth INDEL | 588 | 209 | 167 / 154 | 151 / 140 | 151 / 141 |
| HG008T | Illumina | lost on d9 | 250 | 27 | 26 / 23 | 0 / 0 | 3 / 2 |
| HG008T | Illumina | absorbed on d9 | 114 | 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| HG008T | Illumina | lost on d46 | 235 | 25 | 24 / 20 | 2 / 1 | 0 / 0 |
| HG008T | Illumina | absorbed on d46 | 93 | 2 | 2 / 2 | 0 / 0 | 0 / 0 |
| HG008T | PacBio | all truth INDEL | 588 | 209 | 177 / 153 | 152 / 150 | 152 / 151 |
| HG008T | PacBio | lost on d9 | 228 | 30 | 29 / 27 | 0 / 0 | 2 / 2 |
| HG008T | PacBio | absorbed on d9 | 64 | 1 | 1 / 1 | 0 / 0 | 0 / 0 |
| HG008T | PacBio | lost on d46 | 222 | 31 | 30 / 28 | 3 / 2 | 0 / 0 |
| HG008T | PacBio | absorbed on d46 | 66 | 2 | 2 / 2 | 0 / 0 | 0 / 0 |
| HG008T | ONT | all truth INDEL | 588 | 209 | 140 / 140 | 140 / 134 | 139 / 139 |
| HG008T | ONT | lost on d9 | 141 | 18 | 18 / 18 | 0 / 0 | 1 / 1 |
| HG008T | ONT | absorbed on d9 | 38 | 1 | 1 / 1 | 0 / 0 | 0 / 0 |
| HG008T | ONT | lost on d46 | 136 | 20 | 20 / 20 | 3 / 3 | 0 / 0 |
| HG008T | ONT | absorbed on d46 | 39 | 2 | 2 / 2 | 0 / 0 | 0 / 0 |
| COLO829T | Illumina | all truth INDEL | 131 | 50 | 32 / 29 | pending | 35 / 32 |
| COLO829T | Illumina | lost on d9 | 55 | 11 | 6 / 4 | pending | 0 / 0 |
| COLO829T | Illumina | absorbed on d9 | 26 | 1 | 0 / 0 | pending | 0 / 0 |
| COLO829T | Illumina | lost on d46 | 51 | 11 | 6 / 4 | pending | 0 / 0 |
| COLO829T | Illumina | absorbed on d46 | 23 | 1 | 0 / 0 | pending | 0 / 0 |
| COLO829T | fiberseq | all truth INDEL | 131 | 50 | 39 / 39 | 36 / 36 | 34 / 34 |
| COLO829T | fiberseq | lost on d9 | 60 | 13 | 12 / 12 | 0 / 0 | 0 / 0 |
| COLO829T | fiberseq | absorbed on d9 | 21 | 1 | 1 / 1 | 0 / 0 | 0 / 0 |
| COLO829T | fiberseq | lost on d46 | 55 | 13 | 12 / 12 | 0 / 0 | 0 / 0 |
| COLO829T | fiberseq | absorbed on d46 | 21 | 1 | 1 / 1 | 0 / 0 | 0 / 0 |
| COLO829T | ONT | all truth INDEL | 131 | 50 | 37 / 37 | pending | 32 / 30 |
| COLO829T | ONT | lost on d9 | 54 | 10 | 7 / 7 | pending | 0 / 0 |
| COLO829T | ONT | absorbed on d9 | 16 | 0 | 0 / 0 | pending | 0 / 0 |
| COLO829T | ONT | lost on d46 | 53 | 10 | 7 / 7 | pending | 0 / 0 |
| COLO829T | ONT | absorbed on d46 | 18 | 1 | 1 / 1 | pending | 0 / 0 |

pending rtg outputs: 4

- COLO829_Illumina_INDEL_joint_7fc7396/pon_curve
- COLO829_Illumina_INDEL_joint_7fc7396/pon_calls
- COLO829_ONT_INDEL_joint_7fc7396/pon_curve
- COLO829_ONT_INDEL_joint_7fc7396/pon_calls
