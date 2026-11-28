# Ranges for the R3.2 / R2.5 sentences, recomputed on the ver4.6 tables (no site features on either reference)

Range = min-max over the rows of the group, from the unrounded values, in the cell format of the tables. "smallest n" = the smallest denominator among the rows of the range. Retained percentages of classes with fewer than 20 alleles are left out of a range (stated in the cell).

### FP before population filtering

| quantity (sentence it replaces) | HG008T + COLO829T (six rows) | all nine rows | HCC1395T (three rows) |
|---|---|---|---|
| SNV: graph-represented germline, share of GRCh38-graph FP (R3.2 / Table 3.8 caption "71-86%"; R2.5 "over 70%") | 71-85% (smallest n 28,178) | 71-87% (smallest n 28,178) | 86-87% (smallest n 77,905) |
| SNV: graph-represented germline alleles that keep a candidate on v2.1 ("0.2-0.6%") | 0.1-0.7% (smallest n 22,744) | 0.1-0.7% (smallest n 22,744) | 0.1-0.7% (smallest n 66,643) |
| SNV: graph-absent germline alleles that keep a candidate on v2.1 ("83-88%") | 83-90% (smallest n 3,106) | 83-90% (smallest n 3,106) | 86-87% (smallest n 10,718) |
| SNV: graph-absent germline, share of GRCh38-graph FP | 11-16% (smallest n 28,178) | 11-16% (smallest n 28,178) | 12-14% (smallest n 77,905) |
| SNV: non-germline, share of GRCh38-graph FP | 1-17% (smallest n 28,178) | 0.7-17% (smallest n 28,178) | 0.7-1% (smallest n 77,905) |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, Illumina ("0.1-0.3%") | <0.1-0.2% (smallest n 5,055) | <0.1-0.2% (smallest n 5,055) | <0.1% (smallest n 5,550) |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, long reads ("15-22%") | 15-29% (smallest n 79) | 0.8-29% (smallest n 79) | 0.8-2.6% (smallest n 2,925) |
| INDEL: graph-represented germline, share of GRCh38-graph FP | 43-86% (smallest n 180) | 43-86% (smallest n 180) | 79-84% (smallest n 3,713) |
| INDEL: graph-absent germline alleles that keep a candidate on v2.1 | 25-66% (smallest n 378; 1 of 6 rows left out, class under 20 alleles) | 25-75% (smallest n 378; 1 of 9 rows left out, class under 20 alleles) | 62-75% (smallest n 568) |
| INDEL: non-germline, share of GRCh38-graph FP, ONT ("43-57%") | 13-49% (smallest n 180) | 6-49% (smallest n 180) | 6% (smallest n 3,713) |
| INDEL: non-germline, share of GRCh38-graph FP, Illumina + PacBio HiFi / Fiber-seq | 0.8-36% (smallest n 1,764) | 0.8-36% (smallest n 1,764) | 2-5% (smallest n 6,594) |
| SNV: HPRC v2.1 FP as % of GRCh38-graph FP (last column / first count column) | 2-6% (smallest n 28,178) | 2-6% (smallest n 28,178) | 2-3% (smallest n 77,905) |
| INDEL: HPRC v2.1 FP as % of GRCh38-graph FP | 4-104% (smallest n 180) | 3-104% (smallest n 180) | 3-8% (smallest n 3,713) |

### FP after population filtering

| quantity (sentence it replaces) | HG008T + COLO829T (six rows) | all nine rows | HCC1395T (three rows) |
|---|---|---|---|
| SNV: graph-represented germline, share of GRCh38-graph FP (R3.2 / Table 3.8 caption "71-86%"; R2.5 "over 70%") | 0-0.9% (smallest n 429) | 0-0.9% (smallest n 429) | 0% (smallest n 610) |
| SNV: graph-represented germline alleles that keep a candidate on v2.1 ("0.2-0.6%") | - (6 rows under 20 alleles) | - (9 rows under 20 alleles) | - (3 rows under 20 alleles) |
| SNV: graph-absent germline alleles that keep a candidate on v2.1 ("83-88%") | 75-99% (smallest n 106) | 75-100% (smallest n 106) | 99-100% (smallest n 425) |
| SNV: graph-absent germline, share of GRCh38-graph FP | 7-79% (smallest n 429) | 7-85% (smallest n 429) | 64-85% (smallest n 610) |
| SNV: non-germline, share of GRCh38-graph FP | 20-93% (smallest n 429) | 15-93% (smallest n 429) | 15-36% (smallest n 610) |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, Illumina ("0.1-0.3%") | - (2 rows under 20 alleles) | - (3 rows under 20 alleles) | - (1 rows under 20 alleles) |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, long reads ("15-22%") | 6.4% (smallest n 47; 3 of 4 rows left out, class under 20 alleles) | 6.4% (smallest n 47; 5 of 6 rows left out, class under 20 alleles) | - (2 rows under 20 alleles) |
| INDEL: graph-represented germline, share of GRCh38-graph FP | 0-9% (smallest n 55) | 0-9% (smallest n 55) | 0-0.5% (smallest n 92) |
| INDEL: graph-absent germline alleles that keep a candidate on v2.1 | 25-91% (smallest n 44; 1 of 6 rows left out, class under 20 alleles) | 25-100% (smallest n 28; 1 of 9 rows left out, class under 20 alleles) | 92-100% (smallest n 28) |
| INDEL: non-germline, share of GRCh38-graph FP, ONT ("43-57%") | 80-96% (smallest n 55) | 80-96% (smallest n 55) | 83% (smallest n 168) |
| INDEL: non-germline, share of GRCh38-graph FP, Illumina + PacBio HiFi / Fiber-seq | 14-68% (smallest n 79) | 14-68% (smallest n 79) | 47-63% (smallest n 92) |
| SNV: HPRC v2.1 FP as % of GRCh38-graph FP (last column / first count column) | 13-104% (smallest n 429) | 13-203% (smallest n 429) | 97-203% (smallest n 610) |
| INDEL: HPRC v2.1 FP as % of GRCh38-graph FP | 26-211% (smallest n 55) | 26-211% (smallest n 55) | 44-82% (smallest n 92) |

