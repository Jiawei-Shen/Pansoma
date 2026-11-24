# Ranges for the R3.2 / R2.5 sentences, recomputed on the ver4 tables

Range = min-max over the rows of the group, from the unrounded values, in the cell format of the tables. "smallest n" = the smallest denominator among the rows of the range. Retained percentages of classes with fewer than 20 alleles are left out of a range (stated in the cell).

### FP before population filtering

| quantity (sentence it replaces) | HG008T + COLO829T (six rows) | all nine rows | HCC1395T (three rows) |
|---|---|---|---|
| SNV: graph-represented germline, share of GRCh38-graph FP (R3.2 / Table 3.8 caption "71-86%"; R2.5 "over 70%") | 84% (smallest n 121,139) | 84% (smallest n 121,139) | - |
| SNV: graph-represented germline alleles that keep a candidate on v2.1 ("0.2-0.6%") | 0.3% (smallest n 101,866) | 0.3% (smallest n 101,866) | - |
| SNV: graph-absent germline alleles that keep a candidate on v2.1 ("83-88%") | 85% (smallest n 16,200) | 85% (smallest n 16,200) | - |
| SNV: graph-absent germline, share of GRCh38-graph FP | 13% (smallest n 121,139) | 13% (smallest n 121,139) | - |
| SNV: non-germline, share of GRCh38-graph FP | 3% (smallest n 121,139) | 3% (smallest n 121,139) | - |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, Illumina ("0.1-0.3%") | 0.3% (smallest n 31,450) | 0.3% (smallest n 31,450) | - |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, long reads ("15-22%") | - | - | - |
| INDEL: graph-represented germline, share of GRCh38-graph FP | 79% (smallest n 39,930) | 79% (smallest n 39,930) | - |
| INDEL: graph-absent germline alleles that keep a candidate on v2.1 | 44% (smallest n 6,316) | 44% (smallest n 6,316) | - |
| INDEL: non-germline, share of GRCh38-graph FP, ONT ("43-57%") | - | - | - |
| INDEL: non-germline, share of GRCh38-graph FP, Illumina + PacBio HiFi / Fiber-seq | 5% (smallest n 39,930) | 5% (smallest n 39,930) | - |
| SNV: HPRC v2.1 FP as % of GRCh38-graph FP (last column / first count column) | 2% (smallest n 121,139) | 2% (smallest n 121,139) | - |
| INDEL: HPRC v2.1 FP as % of GRCh38-graph FP | 7% (smallest n 39,930) | 7% (smallest n 39,930) | - |

### FP after population filtering

| quantity (sentence it replaces) | HG008T + COLO829T (six rows) | all nine rows | HCC1395T (three rows) |
|---|---|---|---|
| SNV: graph-represented germline, share of GRCh38-graph FP (R3.2 / Table 3.8 caption "71-86%"; R2.5 "over 70%") | 2% (smallest n 1,089) | 2% (smallest n 1,089) | - |
| SNV: graph-represented germline alleles that keep a candidate on v2.1 ("0.2-0.6%") | - (1 rows under 20 alleles) | - (1 rows under 20 alleles) | - |
| SNV: graph-absent germline alleles that keep a candidate on v2.1 ("83-88%") | 93% (smallest n 625) | 93% (smallest n 625) | - |
| SNV: graph-absent germline, share of GRCh38-graph FP | 57% (smallest n 1,089) | 57% (smallest n 1,089) | - |
| SNV: non-germline, share of GRCh38-graph FP | 41% (smallest n 1,089) | 41% (smallest n 1,089) | - |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, Illumina ("0.1-0.3%") | 0% (smallest n 73) | 0% (smallest n 73) | - |
| INDEL: graph-represented germline alleles that keep a candidate on v2.1, long reads ("15-22%") | - | - | - |
| INDEL: graph-represented germline, share of GRCh38-graph FP | 5% (smallest n 1,348) | 5% (smallest n 1,348) | - |
| INDEL: graph-absent germline alleles that keep a candidate on v2.1 | 56% (smallest n 346) | 56% (smallest n 346) | - |
| INDEL: non-germline, share of GRCh38-graph FP, ONT ("43-57%") | - | - | - |
| INDEL: non-germline, share of GRCh38-graph FP, Illumina + PacBio HiFi / Fiber-seq | 69% (smallest n 1,348) | 69% (smallest n 1,348) | - |
| SNV: HPRC v2.1 FP as % of GRCh38-graph FP (last column / first count column) | 73% (smallest n 1,089) | 73% (smallest n 1,089) | - |
| INDEL: HPRC v2.1 FP as % of GRCh38-graph FP | 29% (smallest n 1,348) | 29% (smallest n 1,348) | - |

