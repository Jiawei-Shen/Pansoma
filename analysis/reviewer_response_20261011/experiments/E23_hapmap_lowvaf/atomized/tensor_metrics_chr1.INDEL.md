# E23 tensor level, chr1 test: HapMap_Illumina_INDEL, INDEL, threshold 0.33516 (validation recall 0.9)

| tensors | somatic | called | TP | precision | recall | F1 |
|---|---|---|---|---|---|---|
| all (52,847) | 30,709 | 28,334 | 23,851 | 0.842 | 0.777 | 0.808 |
| graph-carried left out (51,236) | 29,120 | 27,477 | 23,004 | 0.837 | 0.790 | 0.813 |

Graph-carried left out; recall by truth VAF, precision by tensor allele fraction:

| VAF | recall | precision | F1 |
|---|---|---|---|
| 1-2 % | 0.720 (5,852 / 8,133) | 0.782 (3,725 / 4,763) | 0.750 |
| 2-5 % | 0.777 (2,513 / 3,236) | 0.829 (9,797 / 11,811) | 0.802 |
| 5-10 % | 0.859 (13,104 / 15,251) | 0.901 (7,484 / 8,310) | 0.879 |
| 10-16.5 % | 0.741 (969 / 1,308) | 0.790 (1,193 / 1,511) | 0.764 |
| 1-16.5 % | 0.803 (22,438 / 27,928) | 0.841 (22,199 / 26,395) | 0.822 |
