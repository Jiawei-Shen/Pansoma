## 1. 召回（全部 somatic allele / PASS 且在 somatic BED 内）

| 状态 | 种类 | PacBio | ONT | Illumina |
|---|---|---:|---:|---:|
| 全部 代表 allele | SNP (8,690) | 8,234 (94.8%) | 8,275 (95.2%) | 8,279 (95.3%) |
| 全部 非代表 allele | SNP (8,690) | 2 (0.0%) | 5 (0.1%) | 4 (0.0%) |
| 全部 filtered | SNP (8,690) | 54 (0.6%) | 68 (0.8%) | 68 (0.8%) |
| 全部 没有候选 | SNP (8,690) | 400 (4.6%) | 342 (3.9%) | 339 (3.9%) |
| 全部 代表 allele | DEL (3,545) | 2,102 (59.3%) | 1,864 (52.6%) | 1,785 (50.4%) |
| 全部 非代表 allele | DEL (3,545) | 216 (6.1%) | 446 (12.6%) | 49 (1.4%) |
| 全部 filtered | DEL (3,545) | 337 (9.5%) | 510 (14.4%) | 390 (11.0%) |
| 全部 没有候选 | DEL (3,545) | 890 (25.1%) | 725 (20.5%) | 1,321 (37.3%) |
| 全部 代表 allele | INS (4,951) | 1,959 (39.6%) | 1,764 (35.6%) | 1,770 (35.8%) |
| 全部 非代表 allele | INS (4,951) | 167 (3.4%) | 416 (8.4%) | 29 (0.6%) |
| 全部 filtered | INS (4,951) | 721 (14.6%) | 1,042 (21.0%) | 519 (10.5%) |
| 全部 没有候选 | INS (4,951) | 2,104 (42.5%) | 1,729 (34.9%) | 2,633 (53.2%) |
| PASS∩BED 代表 allele | SNP (8,543) | 8,115 (95.0%) | 8,152 (95.4%) | 8,159 (95.5%) |
| PASS∩BED 非代表 allele | SNP (8,543) | 2 (0.0%) | 5 (0.1%) | 4 (0.0%) |
| PASS∩BED filtered | SNP (8,543) | 51 (0.6%) | 67 (0.8%) | 61 (0.7%) |
| PASS∩BED 没有候选 | SNP (8,543) | 375 (4.4%) | 319 (3.7%) | 319 (3.7%) |
| PASS∩BED 代表 allele | DEL (3,357) | 2,019 (60.1%) | 1,794 (53.4%) | 1,727 (51.4%) |
| PASS∩BED 非代表 allele | DEL (3,357) | 196 (5.8%) | 412 (12.3%) | 47 (1.4%) |
| PASS∩BED filtered | DEL (3,357) | 310 (9.2%) | 475 (14.1%) | 353 (10.5%) |
| PASS∩BED 没有候选 | DEL (3,357) | 832 (24.8%) | 676 (20.1%) | 1,230 (36.6%) |
| PASS∩BED 代表 allele | INS (4,625) | 1,850 (40.0%) | 1,667 (36.0%) | 1,673 (36.2%) |
| PASS∩BED 非代表 allele | INS (4,625) | 156 (3.4%) | 381 (8.2%) | 29 (0.6%) |
| PASS∩BED filtered | INS (4,625) | 653 (14.1%) | 966 (20.9%) | 476 (10.3%) |
| PASS∩BED 没有候选 | INS (4,625) | 1,966 (42.5%) | 1,611 (34.8%) | 2,447 (52.9%) |

## 2. 同一个 truth 在三个平台上的状态（全部 somatic allele）

「T」= 有 tensor（代表或非代表），「F」= filtered，「N」= 没有候选；顺序 PacBio / ONT / Illumina。

| PacBio/ONT/Illumina | SNP | DEL | INS |
|---|---:|---:|---:|
| T/T/T | 8168 | 1643 | 1660 |
| N/N/N | 241 | 560 | 1523 |
| N/F/N | 19 | 115 | 369 |
| T/T/N | 36 | 278 | 104 |
| F/F/N | 2 | 103 | 262 |
| T/T/F | 11 | 176 | 102 |
| F/T/N | 5 | 109 | 158 |
| N/N/F | 39 | 60 | 80 |
| T/F/T | 8 | 85 | 81 |
| T/F/N | 4 | 74 | 85 |
| F/F/F | 4 | 21 | 122 |
| N/T/N | 24 | 50 | 55 |
| F/T/F | 1 | 34 | 81 |
| T/F/F | 2 | 50 | 60 |
| N/N/T | 38 | 46 | 14 |
| F/N/N | 5 | 24 | 62 |
| N/F/F | 2 | 34 | 41 |
| F/F/T | 14 | 13 | 14 |
| N/F/T | 13 | 15 | 8 |
| N/T/T | 22 | 6 | 5 |
| F/N/F | 6 | 11 | 15 |
| F/T/T | 11 | 10 | 6 |
| T/N/N | 3 | 8 | 15 |
| F/N/T | 6 | 12 | 1 |
| T/N/T | 3 | 4 | 10 |
| N/T/F | 2 | 4 | 9 |
| T/N/F | 1 | 0 | 9 |
| **三个平台都有 tensor** | **8168** (94.0%) | **1643** (46.3%) | **1660** (33.5%) |
| **至少一个平台有 tensor** | **8372** (96.3%) | **2617** (73.8%) | **2477** (50.0%) |
| **三个平台都没有候选** | **241** (2.8%) | **560** (15.8%) | **1523** (30.8%) |

## 3. 三个平台都没有候选的 truth：按 PacBio 的类别，看另两个平台是否同类

| 种类 | PacBio 类别 | 个数 | 三个平台同类 |
|---|---|---:|---:|
| DEL | I1 | 276 | 208 |
| DEL | I2 | 242 | 209 |
| DEL | I4 | 14 | 11 |
| DEL | I3 | 12 | 5 |
| DEL | I7 | 6 | 3 |
| DEL | I5 | 6 | 0 |
| DEL | I6 | 4 | 1 |
| INS | I2 | 706 | 514 |
| INS | I1 | 667 | 496 |
| INS | I3 | 71 | 44 |
| INS | I4 | 44 | 26 |
| INS | I7 | 17 | 10 |
| INS | I6 | 11 | 0 |
| INS | I5 | 7 | 1 |
| SNP | S4b | 86 | 37 |
| SNP | S3 | 46 | 29 |
| SNP | S4a | 43 | 21 |
| SNP | S4c | 32 | 24 |
| SNP | S2 | 30 | 14 |
| SNP | S1 | 4 | 1 |

## 4. 没有候选的原因分类（每个平台自己的未命中）

| 种类 | 类别 | PacBio | ONT | Illumina |
|---|---|---:|---:|---:|
| INS | I1_allele_fully_in_graph | 962 | 683 | 1136 |
| INS | I2_partly_in_graph_residual_edit | 961 | 879 | 1056 |
| INS | I3_grch38_path_different_edit | 83 | 95 | 81 |
| INS | I4_no_or_few_alt_reads | 51 | 40 | 264 |
| INS | I5_low_mapq_or_no_reads | 17 | 1 | 84 |
| INS | I6_longer_than_50bp | 11 | 11 | 0 |
| INS | I7_unclear | 19 | 20 | 12 |
| DEL | I1_allele_fully_in_graph | 415 | 306 | 780 |
| DEL | I2_partly_in_graph_residual_edit | 419 | 369 | 448 |
| DEL | I3_grch38_path_different_edit | 17 | 18 | 25 |
| DEL | I4_no_or_few_alt_reads | 16 | 20 | 42 |
| DEL | I5_low_mapq_or_no_reads | 12 | 0 | 17 |
| DEL | I6_longer_than_50bp | 4 | 4 | 1 |
| DEL | I7_unclear | 7 | 7 | 8 |
| DEL | I9_truth_edit_on_unbuilt_node | 0 | 1 | 0 |
| SNP | S1_low_mapq_or_no_reads | 48 | 1 | 55 |
| SNP | S2_alt_absent_in_reads | 46 | 36 | 40 |
| SNP | S3_alt_is_existing_graph_allele | 65 | 46 | 76 |
| SNP | S4a_repeat_snv_edit_on_branch | 70 | 61 | 57 |
| SNP | S4b_repeat_absorbed_by_graph_path | 127 | 154 | 65 |
| SNP | S4c_repeat_no_clear_alt | 44 | 44 | 46 |
| INS | **合计** | **2104** | **1729** | **2633** |
| DEL | **合计** | **890** | **725** | **1321** |
| SNP | **合计** | **400** | **342** | **339** |

## 5. 标签（truth-labels-v2）


**SNV**

| 标签 | reason | PacBio | ONT | Illumina |
|---:|---|---:|---:|---:|
| 1 | representative_allele_is_somatic_truth | 8,234 (0.7%) | 8,275 (0.3%) | 8,279 (0.2%) |
| 2 | representative_allele_is_germline_truth | 332,598 (28.2%) | 343,284 (14.5%) | 311,650 (6.2%) |
| 0 | confident_no_truth_allele | 51,228 (4.3%) | 256,169 (10.8%) | 3,052,652 (61.0%) |
| -1 | outside_confident_region | 589,635 (50.0%) | 1,511,377 (63.7%) | 605,382 (12.1%) |
| -1 | not_on_unique_grch38_node | 106,031 (9.0%) | 109,773 (4.6%) | 208,089 (4.2%) |
| -1 | near_truth_allele_mismatch | 35,590 (3.0%) | 57,732 (2.4%) | 789,784 (15.8%) |
| -1 | truth_matches_non_representative_allele | 209 (0.0%) | 269 (0.0%) | 393 (0.0%) |
| -1 | germline_truth_filtered | 56,814 (4.8%) | 85,627 (3.6%) | 24,659 (0.5%) |
| -1 | somatic_truth_filtered | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) |
| | 合计 | 1,180,339 | 2,372,506 | 5,000,888 |
| | 0 : 1 | 6 : 1 | 31 : 1 | 369 : 1 |

**INDEL**

| 标签 | reason | PacBio | ONT | Illumina |
|---:|---|---:|---:|---:|
| 1 | representative_allele_is_somatic_truth | 4,097 (0.2%) | 3,673 (0.1%) | 3,562 (2.4%) |
| 2 | representative_allele_is_germline_truth | 147,508 (7.7%) | 157,914 (6.3%) | 55,800 (36.8%) |
| 0 | confident_no_truth_allele | 1,307,535 (68.2%) | 1,681,012 (66.9%) | 6,552 (4.3%) |
| -1 | outside_confident_region | 79,965 (4.2%) | 134,661 (5.4%) | 20,628 (13.6%) |
| -1 | not_on_unique_grch38_node | 123,317 (6.4%) | 169,499 (6.8%) | 22,886 (15.1%) |
| -1 | near_truth_allele_mismatch | 241,902 (12.6%) | 321,046 (12.8%) | 39,098 (25.8%) |
| -1 | truth_matches_non_representative_allele | 8,300 (0.4%) | 37,176 (1.5%) | 1,166 (0.8%) |
| -1 | germline_truth_filtered | 4,066 (0.2%) | 6,032 (0.2%) | 1,773 (1.2%) |
| -1 | somatic_truth_filtered | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) |
| | 合计 | 1,916,690 | 2,511,013 | 151,465 |
| | 0 : 1 | 319 : 1 | 458 : 1 | 2 : 1 |

## 6. 标签质量问题

| | PacBio | ONT | Illumina |
|---|---:|---:|---:|
| 残余 edit 的 tensor 被标成 germline（不同候选数） | 301 | 194 | 317 |
| SNV 未命中里 germline 真值在同一位置、同一 ALT | 38 | 35 | 30 |
| site 的代表 allele 是 germline、另一个 allele 是 somatic truth → 标 2 | 254 | 232 | 40 |
