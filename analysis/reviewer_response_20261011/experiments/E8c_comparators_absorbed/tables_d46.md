# E8c: linear tumor-only callers at the alleles the HPRC v2.1 (d46) graph absorbs (chr1, benchmark BED)

Absorbed = no candidate on HPRC v2.1 (d46) on that platform and the ALT spelled by a d46 path (E5b strict test, a lower bound).
Kept = the truth allele is a true positive of the caller at that stage (rtg tp-baseline). In brackets: share of the
absorbed alleles, and for comparison the share of all chr1 truth alleles of that kind in the BED.

| tumor | platform | kind | caller | absorbed / all truth | any record: absorbed (all truth) | PASS | PASS + our PoN |
|---|---|---|---|---|---|---|---|
| HG008T | Illumina | SNV | DeepSomatic-TO | 10 / 697 | 5 (50%) (96%) | 0 (0%) (91%) | 0 (0%) (89%) |
| HG008T | Illumina | SNV | ClairS-TO | 10 / 697 | 3 (30%) (94%) | 0 (0%) (86%) | 0 (0%) (85%) |
| HG008T | Illumina | INDEL | DeepSomatic-TO | 108 / 588 | 88 (81%) (79%) | 5 (5%) (17%) | 0 (0%) (14%) |
| HG008T | Illumina | INDEL | ClairS-TO | 108 / 588 | 49 (45%) (50%) | 0 (0%) (8%) | 0 (0%) (8%) |
| HG008T | PacBio HiFi | SNV | DeepSomatic-TO | 10 / 697 | 7 (70%) (96%) | 0 (0%) (82%) | 0 (0%) (81%) |
| HG008T | PacBio HiFi | SNV | ClairS-TO | 10 / 697 | 4 (40%) (95%) | 0 (0%) (86%) | 0 (0%) (84%) |
| HG008T | PacBio HiFi | INDEL | DeepSomatic-TO | 78 / 588 | 33 (42%) (47%) | 0 (0%) (4%) | 0 (0%) (3%) |
| HG008T | PacBio HiFi | INDEL | ClairS-TO | 78 / 588 | 2 (3%) (17%) | 0 (0%) (7%) | 0 (0%) (7%) |
| HG008T | ONT | SNV | DeepSomatic-TO | 10 / 697 | 5 (50%) (95%) | 0 (0%) (89%) | 0 (0%) (87%) |
| HG008T | ONT | SNV | ClairS-TO | 10 / 697 | 3 (30%) (95%) | 0 (0%) (85%) | 0 (0%) (83%) |
| HG008T | ONT | INDEL | DeepSomatic-TO | 64 / 588 | 23 (36%) (24%) | 1 (2%) (2%) | 1 (2%) (2%) |
| HG008T | ONT | INDEL | ClairS-TO | 64 / 588 | 4 (6%) (14%) | 0 (0%) (7%) | 0 (0%) (7%) |
| COLO829T | Illumina | SNV | DeepSomatic-TO | 3 / 2443 | 3 (100%) (88%) | 0 (0%) (85%) | 0 (0%) (83%) |
| COLO829T | Illumina | SNV | ClairS-TO | 3 / 2443 | 3 (100%) (82%) | 0 (0%) (77%) | 0 (0%) (76%) |
| COLO829T | Illumina | INDEL | DeepSomatic-TO | 34 / 131 | 23 (68%) (82%) | 1 (3%) (21%) | 0 (0%) (18%) |
| COLO829T | Illumina | INDEL | ClairS-TO | 34 / 131 | 22 (65%) (89%) | 0 (0%) (15%) | 0 (0%) (15%) |
| COLO829T | Fiber-seq | SNV | DeepSomatic-TO | 3 / 2443 | 3 (100%) (92%) | 0 (0%) (56%) | 0 (0%) (56%) |
| COLO829T | Fiber-seq | SNV | ClairS-TO | 3 / 2443 | 2 (67%) (83%) | 0 (0%) (75%) | 0 (0%) (74%) |
| COLO829T | Fiber-seq | INDEL | DeepSomatic-TO | 26 / 131 | 16 (62%) (64%) | 0 (0%) (9%) | 0 (0%) (8%) |
| COLO829T | Fiber-seq | INDEL | ClairS-TO | 26 / 131 | 7 (27%) (24%) | 0 (0%) (11%) | 0 (0%) (11%) |
| COLO829T | ONT | SNV | DeepSomatic-TO | 3 / 2443 | 3 (100%) (83%) | 0 (0%) (70%) | 0 (0%) (69%) |
| COLO829T | ONT | SNV | ClairS-TO | 3 / 2443 | 2 (67%) (79%) | 0 (0%) (69%) | 0 (0%) (68%) |
| COLO829T | ONT | INDEL | DeepSomatic-TO | 21 / 131 | 12 (57%) (34%) | 0 (0%) (2%) | 0 (0%) (2%) |
| COLO829T | ONT | INDEL | ClairS-TO | 21 / 131 | 6 (29%) (20%) | 0 (0%) (9%) | 0 (0%) (9%) |
| HCC1395T | Illumina | SNV | DeepSomatic-TO | 4 / 3440 | 4 (100%) (96%) | 0 (0%) (93%) | 0 (0%) (93%) |
| HCC1395T | Illumina | SNV | ClairS-TO | 4 / 3440 | 4 (100%) (92%) | 0 (0%) (86%) | 0 (0%) (86%) |
| HCC1395T | Illumina | INDEL | DeepSomatic-TO | 4 / 133 | 4 (100%) (92%) | 2 (50%) (80%) | 0 (0%) (77%) |
| HCC1395T | Illumina | INDEL | ClairS-TO | 4 / 133 | 3 (75%) (86%) | 0 (0%) (69%) | 0 (0%) (69%) |
| HCC1395T | PacBio HiFi | SNV | DeepSomatic-TO | 3 / 3440 | 3 (100%) (96%) | 0 (0%) (93%) | 0 (0%) (93%) |
| HCC1395T | PacBio HiFi | SNV | ClairS-TO | 3 / 3440 | 2 (67%) (88%) | 0 (0%) (83%) | 0 (0%) (82%) |
| HCC1395T | PacBio HiFi | INDEL | DeepSomatic-TO | 3 / 133 | 3 (100%) (82%) | 0 (0%) (63%) | 0 (0%) (62%) |
| HCC1395T | PacBio HiFi | INDEL | ClairS-TO | 3 / 133 | 1 (33%) (62%) | 0 (0%) (52%) | 0 (0%) (52%) |
| HCC1395T | ONT | SNV | DeepSomatic-TO | 4 / 3440 | 3 (75%) (80%) | 0 (0%) (77%) | 0 (0%) (77%) |
| HCC1395T | ONT | SNV | ClairS-TO | 4 / 3440 | 3 (75%) (77%) | 0 (0%) (73%) | 0 (0%) (73%) |
| HCC1395T | ONT | INDEL | DeepSomatic-TO | 2 / 133 | 2 (100%) (63%) | 0 (0%) (48%) | 0 (0%) (48%) |
| HCC1395T | ONT | INDEL | ClairS-TO | 2 / 133 | 1 (50%) (62%) | 0 (0%) (56%) | 0 (0%) (56%) |
