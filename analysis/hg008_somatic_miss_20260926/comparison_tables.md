## 1. 召回（全部 somatic allele / PASS 且在 somatic BED 内）

| 状态 | 种类 | PacBio | ONT | Illumina |
|---|---|---:|---:|---:|
| 全部 代表 allele | SNP (8,690) | 8,234 (94.8%) | 8,275 (95.2%) | 8,279 (95.3%) |
| 全部 非代表 allele | SNP (8,690) | 2 (0.0%) | 5 (0.1%) | 4 (0.0%) |
| 全部 filtered | SNP (8,690) | 56 (0.6%) | 68 (0.8%) | 68 (0.8%) |
| 全部 没有候选 | SNP (8,690) | 398 (4.6%) | 342 (3.9%) | 339 (3.9%) |
| 全部 代表 allele | DEL (3,545) | 2,102 (59.3%) | 1,864 (52.6%) | 1,785 (50.4%) |
| 全部 非代表 allele | DEL (3,545) | 216 (6.1%) | 446 (12.6%) | 49 (1.4%) |
| 全部 filtered | DEL (3,545) | 344 (9.7%) | 510 (14.4%) | 390 (11.0%) |
| 全部 没有候选 | DEL (3,545) | 883 (24.9%) | 725 (20.5%) | 1,321 (37.3%) |
| 全部 代表 allele | INS (4,951) | 1,960 (39.6%) | 1,764 (35.6%) | 1,770 (35.8%) |
| 全部 非代表 allele | INS (4,951) | 167 (3.4%) | 416 (8.4%) | 29 (0.6%) |
| 全部 filtered | INS (4,951) | 731 (14.8%) | 1,042 (21.0%) | 519 (10.5%) |
| 全部 没有候选 | INS (4,951) | 2,093 (42.3%) | 1,729 (34.9%) | 2,633 (53.2%) |
| PASS∩BED 代表 allele | SNP (8,543) | 8,115 (95.0%) | 8,152 (95.4%) | 8,159 (95.5%) |
| PASS∩BED 非代表 allele | SNP (8,543) | 2 (0.0%) | 5 (0.1%) | 4 (0.0%) |
| PASS∩BED filtered | SNP (8,543) | 53 (0.6%) | 67 (0.8%) | 61 (0.7%) |
| PASS∩BED 没有候选 | SNP (8,543) | 373 (4.4%) | 319 (3.7%) | 319 (3.7%) |
| PASS∩BED 代表 allele | DEL (3,357) | 2,019 (60.1%) | 1,794 (53.4%) | 1,727 (51.4%) |
| PASS∩BED 非代表 allele | DEL (3,357) | 196 (5.8%) | 412 (12.3%) | 47 (1.4%) |
| PASS∩BED filtered | DEL (3,357) | 317 (9.4%) | 475 (14.1%) | 353 (10.5%) |
| PASS∩BED 没有候选 | DEL (3,357) | 825 (24.6%) | 676 (20.1%) | 1,230 (36.6%) |
| PASS∩BED 代表 allele | INS (4,625) | 1,851 (40.0%) | 1,667 (36.0%) | 1,673 (36.2%) |
| PASS∩BED 非代表 allele | INS (4,625) | 156 (3.4%) | 381 (8.2%) | 29 (0.6%) |
| PASS∩BED filtered | INS (4,625) | 662 (14.3%) | 966 (20.9%) | 476 (10.3%) |
| PASS∩BED 没有候选 | INS (4,625) | 1,956 (42.3%) | 1,611 (34.8%) | 2,447 (52.9%) |

## 2. 同一个 truth 在三个平台上的状态（全部 somatic allele）

「T」= 有 tensor（代表或非代表），「F」= filtered，「N」= 没有候选；顺序 PacBio / ONT / Illumina。

| PacBio/ONT/Illumina | SNP | DEL | INS |
|---|---:|---:|---:|
| T/T/T | 8168 | 1643 | 1660 |
| N/N/N | 241 | 560 | 1522 |
| N/F/N | 19 | 114 | 363 |
| T/T/N | 36 | 278 | 105 |
| F/F/N | 2 | 104 | 268 |
| T/T/F | 11 | 176 | 102 |
| F/T/N | 5 | 111 | 160 |
| N/N/F | 39 | 60 | 80 |
| T/F/T | 8 | 85 | 81 |
| T/F/N | 4 | 74 | 85 |
| F/F/F | 4 | 22 | 123 |
| N/T/N | 24 | 48 | 52 |
| F/T/F | 1 | 34 | 81 |
| T/F/F | 2 | 50 | 60 |
| N/N/T | 38 | 44 | 14 |
| F/N/N | 5 | 24 | 63 |
| N/F/F | 2 | 33 | 40 |
| F/F/T | 16 | 13 | 14 |
| N/F/T | 11 | 15 | 8 |
| F/N/F | 6 | 11 | 15 |
| N/T/T | 22 | 5 | 5 |
| F/T/T | 11 | 11 | 6 |
| T/N/N | 3 | 8 | 15 |
| F/N/T | 6 | 14 | 1 |
| T/N/T | 3 | 4 | 10 |
| N/T/F | 2 | 4 | 9 |
| T/N/F | 1 | 0 | 9 |
| **三个平台都有 tensor** | **8168** (94.0%) | **1643** (46.3%) | **1660** (33.5%) |
| **至少一个平台有 tensor** | **8372** (96.3%) | **2617** (73.8%) | **2477** (50.0%) |
| **三个平台都没有候选** | **241** (2.8%) | **560** (15.8%) | **1522** (30.7%) |

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
| INS | I1 | 666 | 496 |
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

## 5. 标签（PacBio truth-labels rules 197b5d25、ONT truth-labels rules 197b5d25、Illumina truth-labels rules 197b5d25）


**SNV**

| 标签 | reason | PacBio | ONT | Illumina |
|---:|---|---:|---:|---:|
| 1 | allele_partial_somatic_truth | 1 (0.0%) | 3 (0.0%) | 2 (0.0%) |
| 1 | representative_allele_is_somatic_truth | 8,234 (0.7%) | 8,275 (0.3%) | 8,277 (0.2%) |
| 1 | residual_partial_somatic_truth | 601 (0.1%) | 570 (0.0%) | 2,347 (0.0%) |
| 2 | representative_allele_is_germline_truth | 332,598 (28.2%) | 343,284 (14.5%) | 310,349 (6.2%) |
| 2 | representative_allele_is_germline_truth_gap_filtered | 18,229 (1.5%) | 26,891 (1.1%) | 5,089 (0.1%) |
| 0 | confident_no_truth_allele | 51,229 (4.3%) | 256,169 (10.8%) | 2,041,479 (40.8%) |
| 0 | near_truth_allele_mismatch | 35,311 (3.0%) | 57,455 (2.4%) | 594,093 (11.9%) |
| 0 | truth_matches_non_representative_allele | 137 (0.0%) | 160 (0.0%) | 325 (0.0%) |
| -1 | below_snv_min_af | 0 (0.0%) | 0 (0.0%) | 1,312,554 (26.2%) |
| -1 | not_on_unique_grch38_node | 33,043 (2.8%) | 37,329 (1.6%) | 20,923 (0.4%) |
| -1 | off_reference_no_truth_match | 48,953 (4.1%) | 48,284 (2.0%) | 126,818 (2.5%) |
| -1 | outside_confident_region | 652,007 (55.2%) | 1,594,086 (67.2%) | 578,632 (11.6%) |
| | **1 合计** | **8,836** | **8,848** | **10,626** |
| | **2 合计** | **350,827** | **370,175** | **315,438** |
| | **0 合计** | **86,677** | **313,784** | **2,635,897** |
| | **−1 合计** | **734,003** | **1,679,699** | **2,038,927** |
| | 总数 | 1,180,343 | 2,372,506 | 5,000,888 |
| | 0 : 1 | 10 : 1 | 35 : 1 | 248 : 1 |

**INDEL**

| 标签 | reason | PacBio | ONT | Illumina |
|---:|---|---:|---:|---:|
| 1 | allele_partial_somatic_truth | 93 (0.0%) | 391 (0.0%) | 32 (0.0%) |
| 1 | representative_allele_is_somatic_truth | 4,098 (0.2%) | 3,673 (0.1%) | 3,562 (2.4%) |
| 1 | residual_partial_somatic_truth | 2,649 (0.1%) | 2,319 (0.1%) | 1,645 (1.1%) |
| 2 | allele_partial_germline_truth | 4,540 (0.2%) | 16,407 (0.7%) | 609 (0.4%) |
| 2 | representative_allele_is_germline_truth | 147,786 (7.7%) | 157,914 (6.3%) | 55,800 (36.8%) |
| 2 | representative_allele_is_germline_truth_gap_filtered | 1,544 (0.1%) | 2,177 (0.1%) | 373 (0.2%) |
| 2 | residual_partial_germline_truth | 37,554 (2.0%) | 44,265 (1.8%) | 13,768 (9.1%) |
| 0 | confident_no_truth_allele | 1,309,897 (68.1%) | 1,681,012 (66.9%) | 6,552 (4.3%) |
| 0 | near_truth_allele_mismatch | 205,183 (10.7%) | 275,572 (11.0%) | 24,343 (16.1%) |
| 0 | truth_matches_non_representative_allele | 3,577 (0.2%) | 20,053 (0.8%) | 497 (0.3%) |
| -1 | not_on_unique_grch38_node | 8,684 (0.5%) | 12,075 (0.5%) | 973 (0.6%) |
| -1 | off_reference_no_truth_match | 99,441 (5.2%) | 138,291 (5.5%) | 18,218 (12.0%) |
| -1 | outside_confident_region | 98,698 (5.1%) | 156,864 (6.2%) | 25,093 (16.6%) |
| | **1 合计** | **6,840** | **6,383** | **5,239** |
| | **2 合计** | **191,424** | **220,763** | **70,550** |
| | **0 合计** | **1,518,657** | **1,976,637** | **31,392** |
| | **−1 合计** | **206,823** | **307,230** | **44,284** |
| | 总数 | 1,923,744 | 2,511,013 | 151,465 |
| | 0 : 1 | 222 : 1 | 310 : 1 | 6 : 1 |

## 6. 标签质量问题

| | PacBio | ONT | Illumina |
|---|---:|---:|---:|
| 残余 edit 的 tensor 被标成 germline（不同候选数） | 432 | 271 | 330 |
| SNV 未命中里 germline 真值在同一位置、同一 ALT | 38 | 35 | 30 |
| site 的代表 allele 是 germline、另一个 allele 是 somatic truth → 标 2 | 254 | 234 | 40 |
