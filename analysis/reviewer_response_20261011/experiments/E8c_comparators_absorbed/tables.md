# E8c: linear tumor-only callers at the alleles the HPRC v1.1 graph absorbs (chr1, benchmark BED)

Absorbed = ALT reads align without edit to an existing HPRC v1.1 (d9) path on that platform (no Pansoma candidate).
Kept = the truth allele is a true positive of the caller at that stage (rtg tp-baseline). In brackets: share of the
absorbed alleles, and for comparison the share of all chr1 truth alleles of that kind in the BED.

| tumor | platform | kind | caller | absorbed / all truth | any record: absorbed (all truth) | PASS | PASS + our PoN |
|---|---|---|---|---|---|---|---|
| HG008T | Illumina | SNV | DeepSomatic-TO | 12 / 697 | 4 (33%) (96%) | 0 (0%) (91%) | 0 (0%) (89%) |
| HG008T | Illumina | SNV | ClairS-TO | 12 / 697 | 3 (25%) (94%) | 0 (0%) (86%) | 0 (0%) (85%) |
| HG008T | Illumina | INDEL | DeepSomatic-TO | 136 / 588 | 115 (85%) (79%) | 6 (4%) (17%) | 0 (0%) (14%) |
| HG008T | Illumina | INDEL | ClairS-TO | 136 / 588 | 60 (44%) (50%) | 0 (0%) (8%) | 0 (0%) (8%) |
| HG008T | PacBio HiFi | SNV | DeepSomatic-TO | 14 / 697 | 7 (50%) (96%) | 0 (0%) (82%) | 0 (0%) (81%) |
| HG008T | PacBio HiFi | SNV | ClairS-TO | 14 / 697 | 4 (29%) (95%) | 0 (0%) (86%) | 0 (0%) (84%) |
| HG008T | PacBio HiFi | INDEL | DeepSomatic-TO | 102 / 588 | 42 (41%) (47%) | 0 (0%) (4%) | 0 (0%) (3%) |
| HG008T | PacBio HiFi | INDEL | ClairS-TO | 102 / 588 | 4 (4%) (17%) | 0 (0%) (7%) | 0 (0%) (7%) |
| HG008T | ONT | SNV | DeepSomatic-TO | 14 / 697 | 5 (36%) (95%) | 0 (0%) (89%) | 0 (0%) (87%) |
| HG008T | ONT | SNV | ClairS-TO | 14 / 697 | 4 (29%) (95%) | 0 (0%) (85%) | 0 (0%) (83%) |
| HG008T | ONT | INDEL | DeepSomatic-TO | 74 / 588 | 21 (28%) (24%) | 0 (0%) (2%) | 0 (0%) (2%) |
| HG008T | ONT | INDEL | ClairS-TO | 74 / 588 | 5 (7%) (14%) | 0 (0%) (7%) | 0 (0%) (7%) |
| COLO829T | Illumina | SNV | DeepSomatic-TO | 4 / 2443 | 3 (75%) (88%) | 0 (0%) (85%) | 0 (0%) (83%) |
| COLO829T | Illumina | SNV | ClairS-TO | 4 / 2443 | 3 (75%) (82%) | 0 (0%) (77%) | 0 (0%) (76%) |
| COLO829T | Illumina | INDEL | DeepSomatic-TO | 34 / 131 | 26 (76%) (82%) | 1 (3%) (21%) | 0 (0%) (18%) |
| COLO829T | Illumina | INDEL | ClairS-TO | 34 / 131 | 25 (74%) (89%) | 0 (0%) (15%) | 0 (0%) (15%) |
| COLO829T | Fiber-seq | SNV | DeepSomatic-TO | 5 / 2443 | 3 (60%) (92%) | 0 (0%) (56%) | 0 (0%) (56%) |
| COLO829T | Fiber-seq | SNV | ClairS-TO | 5 / 2443 | 2 (40%) (83%) | 0 (0%) (75%) | 0 (0%) (74%) |
| COLO829T | Fiber-seq | INDEL | DeepSomatic-TO | 23 / 131 | 19 (83%) (64%) | 0 (0%) (9%) | 0 (0%) (8%) |
| COLO829T | Fiber-seq | INDEL | ClairS-TO | 23 / 131 | 5 (22%) (24%) | 0 (0%) (11%) | 0 (0%) (11%) |
| COLO829T | ONT | SNV | DeepSomatic-TO | 4 / 2443 | 3 (75%) (83%) | 0 (0%) (70%) | 0 (0%) (69%) |
| COLO829T | ONT | SNV | ClairS-TO | 4 / 2443 | 2 (50%) (79%) | 0 (0%) (69%) | 0 (0%) (68%) |
| COLO829T | ONT | INDEL | DeepSomatic-TO | 19 / 131 | 10 (53%) (34%) | 0 (0%) (2%) | 0 (0%) (2%) |
| COLO829T | ONT | INDEL | ClairS-TO | 19 / 131 | 6 (32%) (20%) | 0 (0%) (9%) | 0 (0%) (9%) |
