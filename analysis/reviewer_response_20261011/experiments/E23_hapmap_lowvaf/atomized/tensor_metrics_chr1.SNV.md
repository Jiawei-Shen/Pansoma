# E23 tensor level, chr1 test: HapMap_Illumina_SNV_full, SNV, threshold 0.98756 (validation recall 0.9)

| tensors | somatic | called | TP | precision | recall | F1 |
|---|---|---|---|---|---|---|
| all (3,354,082) | 132,899 | 120,142 | 119,148 | 0.992 | 0.897 | 0.942 |
| graph-carried left out (3,351,624) | 130,542 | 120,047 | 119,053 | 0.992 | 0.912 | 0.950 |

Graph-carried left out; recall by truth VAF, precision by tensor allele fraction:

| VAF | recall | precision | F1 |
|---|---|---|---|
| 1-2 % | 0.858 (37,616 / 43,841) | 0.991 (24,104 / 24,320) | 0.920 |
| 2-5 % | 0.948 (12,379 / 13,061) | 0.987 (32,414 / 32,851) | 0.967 |
| 5-10 % | 0.951 (64,913 / 68,254) | 0.995 (57,474 / 57,772) | 0.972 |
| 10-16.5 % | 0.901 (3,400 / 3,775) | 0.991 (4,873 / 4,915) | 0.944 |
| 1-16.5 % | 0.918 (118,308 / 128,931) | 0.992 (118,865 / 119,858) | 0.953 |
