# E9k part 5: tensors on the HPRC graphs for all COLO829BL-DSA somatic SNVs (dsa_all_tensors.py)

44,810 DSA somatic SNVs. A tensor has the SNV when its representative allele A1 is the SNV (by the candidate key
of the DSA window's edit on the graph, or by the SNV's GRCh38 lift); "site" = a tensor there with another A1.

## SNVs with a tensor (exact A1)

| DSA class | SNVs | COLO829T_Illumina | COLO829T_fiberseq | COLO829T_ONT | COLO829T_Illumina_d46 | COLO829T_fiberseq_d46 | COLO829T_ONT_d46 | any v1.1 set | any v2.1 set | any of the six |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| lifts, same base, in the SMaHT GRCh38 truth | 35,150 | 34,924 (99.4%) | 34,906 (99.3%) | 34,983 (99.5%) | 34,926 (99.4%) | 34,897 (99.3%) | 34,995 (99.6%) | 35,053 (99.7%) | 35,071 (99.8%) | 35,091 (99.8%) |
| lifts, same base, not in the truth | 4,160 | 1,216 (29.2%) | 1,229 (29.5%) | 1,405 (33.8%) | 1,210 (29.1%) | 1,229 (29.5%) | 1,401 (33.7%) | 1,500 (36.1%) | 1,520 (36.5%) | 1,574 (37.8%) |
| lifts, on a germline allele | 148 | 39 (26.4%) | 31 (20.9%) | 32 (21.6%) | 33 (22.3%) | 28 (18.9%) | 29 (19.6%) | 41 (27.7%) | 35 (23.6%) | 42 (28.4%) |
| no GRCh38 position in the chain (DSA only) | 5,352 | 392 (7.3%) | 280 (5.2%) | 275 (5.1%) | 392 (7.3%) | 322 (6.0%) | 287 (5.4%) | 428 (8.0%) | 467 (8.7%) | 516 (9.6%) |
| all | 44,810 | 36,571 (81.6%) | 36,446 (81.3%) | 36,695 (81.9%) | 36,561 (81.6%) | 36,476 (81.4%) | 36,712 (81.9%) | 37,022 (82.6%) | 37,093 (82.8%) | 37,223 (83.1%) |

## Where each SNV lands on the graph (DSA window, E9j classes) and whether a set has its tensor

### d9

| DSA class | graph class | SNVs | COLO829T_Illumina: exact / site | COLO829T_fiberseq: exact / site | COLO829T_ONT: exact / site |
|---|---|---:|---:|---:|---:|
| lifts, same base, in the SMaHT GRCh38 truth | A | 34,653 | 34,600 / 0 | 34,573 / 0 | 34,573 / 0 |
| lifts, same base, in the SMaHT GRCh38 truth | B | 24 | 23 / 0 | 21 / 0 | 22 / 0 |
| lifts, same base, in the SMaHT GRCh38 truth | C | 46 | 0 / 0 | 1 / 0 | 1 / 1 |
| lifts, same base, in the SMaHT GRCh38 truth | D | 427 | 301 / 0 | 311 / 0 | 387 / 0 |
| lifts, same base, not in the truth | A | 2,858 | 1,072 / 3 | 994 / 3 | 1,011 / 3 |
| lifts, same base, not in the truth | B | 34 | 28 / 0 | 28 / 0 | 29 / 0 |
| lifts, same base, not in the truth | C | 38 | 2 / 0 | 2 / 0 | 3 / 0 |
| lifts, same base, not in the truth | D | 1,230 | 114 / 1 | 205 / 2 | 362 / 4 |
| lifts, on a germline allele | A | 18 | 12 / 0 | 9 / 1 | 10 / 1 |
| lifts, on a germline allele | B | 30 | 27 / 0 | 22 / 0 | 22 / 0 |
| lifts, on a germline allele | C | 48 | 0 / 0 | 0 / 0 | 0 / 1 |
| lifts, on a germline allele | D | 52 | 0 / 0 | 0 / 5 | 0 / 8 |
| no GRCh38 position in the chain (DSA only) | A | 2,134 | 266 / 23 | 148 / 21 | 144 / 30 |
| no GRCh38 position in the chain (DSA only) | B | 183 | 126 / 3 | 132 / 1 | 131 / 1 |
| no GRCh38 position in the chain (DSA only) | C | 43 | 0 / 0 | 0 / 0 | 0 / 0 |
| no GRCh38 position in the chain (DSA only) | D | 2,992 | 0 / 0 | 0 / 0 | 0 / 0 |

### d46

| DSA class | graph class | SNVs | COLO829T_Illumina_d46: exact / site | COLO829T_fiberseq_d46: exact / site | COLO829T_ONT_d46: exact / site |
|---|---|---:|---:|---:|---:|
| lifts, same base, in the SMaHT GRCh38 truth | A | 34,638 | 34,585 / 0 | 34,551 / 0 | 34,555 / 0 |
| lifts, same base, in the SMaHT GRCh38 truth | B | 29 | 28 / 0 | 25 / 0 | 27 / 0 |
| lifts, same base, in the SMaHT GRCh38 truth | C | 39 | 1 / 0 | 1 / 0 | 2 / 0 |
| lifts, same base, in the SMaHT GRCh38 truth | D | 444 | 312 / 0 | 320 / 0 | 411 / 0 |
| lifts, same base, not in the truth | A | 2,840 | 1,063 / 3 | 983 / 8 | 996 / 2 |
| lifts, same base, not in the truth | B | 34 | 28 / 0 | 28 / 0 | 27 / 0 |
| lifts, same base, not in the truth | C | 40 | 4 / 0 | 5 / 0 | 3 / 0 |
| lifts, same base, not in the truth | D | 1,246 | 115 / 1 | 213 / 7 | 375 / 3 |
| lifts, on a germline allele | A | 13 | 6 / 0 | 3 / 3 | 5 / 1 |
| lifts, on a germline allele | B | 29 | 27 / 0 | 25 / 1 | 24 / 1 |
| lifts, on a germline allele | C | 52 | 0 / 0 | 0 / 0 | 0 / 0 |
| lifts, on a germline allele | D | 54 | 0 / 0 | 0 / 4 | 0 / 7 |
| no GRCh38 position in the chain (DSA only) | A | 1,834 | 253 / 16 | 167 / 39 | 145 / 32 |
| no GRCh38 position in the chain (DSA only) | B | 422 | 139 / 5 | 155 / 6 | 142 / 4 |
| no GRCh38 position in the chain (DSA only) | C | 41 | 0 / 0 | 0 / 0 | 0 / 0 |
| no GRCh38 position in the chain (DSA only) | D | 3,055 | 0 / 0 | 0 / 0 | 0 / 0 |

## Labels of the exact tensors (all SNVs)

- COLO829T_Illumina: somatic / representative_allele_is_somatic_truth 35,033; non / confident_no_truth_allele 778; ignore / outside_confident_region 321; non / near_truth_allele_mismatch 154; ignore / off_reference_no_truth_match 122; somatic / residual_partial_somatic_truth 108; ignore / not_on_unique_grch38_node 34; germline / representative_allele_is_germline_truth 20; somatic / allele_partial_somatic_truth 1 (found by key+grch38 35,544, key 598, grch38 429)
- COLO829T_fiberseq: somatic / representative_allele_is_somatic_truth 34,996; non / confident_no_truth_allele 819; ignore / outside_confident_region 198; non / near_truth_allele_mismatch 166; ignore / off_reference_no_truth_match 116; somatic / residual_partial_somatic_truth 104; ignore / not_on_unique_grch38_node 33; germline / representative_allele_is_germline_truth 14 (found by key+grch38 35,497, grch38 551, key 398)
- COLO829T_ONT: somatic / representative_allele_is_somatic_truth 35,076; non / confident_no_truth_allele 884; ignore / outside_confident_region 302; non / near_truth_allele_mismatch 174; ignore / off_reference_no_truth_match 115; somatic / residual_partial_somatic_truth 96; ignore / not_on_unique_grch38_node 33; germline / representative_allele_is_germline_truth 14; somatic / allele_partial_somatic_truth 1 (found by key+grch38 35,510, grch38 792, key 393)
- COLO829T_Illumina_d46: somatic / representative_allele_is_somatic_truth 35,031; non / confident_no_truth_allele 774; ignore / outside_confident_region 305; non / near_truth_allele_mismatch 149; ignore / off_reference_no_truth_match 117; somatic / residual_partial_somatic_truth 105; ignore / not_on_unique_grch38_node 61; germline / representative_allele_is_germline_truth 18; somatic / allele_partial_somatic_truth 1 (found by key+grch38 35,512, key 606, grch38 443)
- COLO829T_fiberseq_d46: somatic / representative_allele_is_somatic_truth 34,977; non / confident_no_truth_allele 815; ignore / outside_confident_region 221; non / near_truth_allele_mismatch 164; ignore / off_reference_no_truth_match 114; somatic / residual_partial_somatic_truth 101; ignore / not_on_unique_grch38_node 70; germline / representative_allele_is_germline_truth 14 (found by key+grch38 35,468, grch38 566, key 442)
- COLO829T_ONT_d46: somatic / representative_allele_is_somatic_truth 35,080; non / confident_no_truth_allele 890; ignore / outside_confident_region 289; non / near_truth_allele_mismatch 171; ignore / off_reference_no_truth_match 110; somatic / residual_partial_somatic_truth 95; ignore / not_on_unique_grch38_node 63; germline / representative_allele_is_germline_truth 13; somatic / allele_partial_somatic_truth 1 (found by key+grch38 35,494, grch38 823, key 395)
