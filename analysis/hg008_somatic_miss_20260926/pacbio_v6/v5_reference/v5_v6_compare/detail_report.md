
### I2 部分在图里 + 残余 edit  (DEL 429, INS 817)

| 最好的相关 tensor 标签 | DEL v5 | DEL v6 | DEL v6new | INS v5 | INS v6 | INS v6new |
|---|---:|---:|---:|---:|---:|---:|
| somatic (1) | 2 | 4 | 4 | 2 | 2 | 2 |
| germline (2) | 71 | 66 | 66 | 74 | 71 | 71 |
| non (0) | 4 | 4 | 3 | 4 | 48 | 2 |
| −1 near_truth_allele_mismatch | 294 | 253 | 254 | 510 | 374 | 420 |
| −1 truth_matches_non_representative_allele | 4 | 36 | 36 | 10 | 25 | 25 |
| −1 outside_confident_region | 20 | 14 | 14 | 38 | 30 | 30 |
| −1 not_on_unique_grch38_node | 7 | 18 | 18 | 119 | 175 | 175 |
| −1 other | 0 | 0 | 0 | 0 | 0 | 0 |
| only SNV tensor | 22 | 25 | 25 | 49 | 42 | 42 |
| no related tensor | 5 | 9 | 9 | 11 | 50 | 50 |
| **有相关 INDEL tensor** | **402** | **395** | **395** | **757** | **725** | **725** |

### I3 走 GRCh38，edits 写法不同  (DEL 55, INS 80)

| 最好的相关 tensor 标签 | DEL v5 | DEL v6 | DEL v6new | INS v5 | INS v6 | INS v6new |
|---|---:|---:|---:|---:|---:|---:|
| somatic (1) | 0 | 31 | 31 | 0 | 0 | 0 |
| germline (2) | 4 | 3 | 3 | 4 | 4 | 4 |
| non (0) | 0 | 0 | 0 | 0 | 0 | 0 |
| −1 near_truth_allele_mismatch | 41 | 13 | 13 | 34 | 41 | 41 |
| −1 truth_matches_non_representative_allele | 0 | 1 | 1 | 0 | 0 | 0 |
| −1 outside_confident_region | 6 | 1 | 1 | 0 | 1 | 1 |
| −1 not_on_unique_grch38_node | 1 | 1 | 1 | 9 | 5 | 5 |
| −1 other | 0 | 0 | 0 | 0 | 0 | 0 |
| only SNV tensor | 3 | 3 | 3 | 31 | 27 | 27 |
| no related tensor | 0 | 2 | 2 | 2 | 2 | 2 |
| **有相关 INDEL tensor** | **52** | **50** | **50** | **47** | **51** | **51** |

### I1 完全在图里（没有 edits）  (DEL 435, INS 827)

| 最好的相关 tensor 标签 | DEL v5 | DEL v6 | DEL v6new | INS v5 | INS v6 | INS v6new |
|---|---:|---:|---:|---:|---:|---:|
| somatic (1) | 1 | 0 | 0 | 3 | 6 | 6 |
| germline (2) | 32 | 41 | 41 | 42 | 45 | 45 |
| non (0) | 2 | 1 | 1 | 0 | 12 | 2 |
| −1 near_truth_allele_mismatch | 128 | 141 | 141 | 182 | 215 | 225 |
| −1 truth_matches_non_representative_allele | 0 | 3 | 3 | 1 | 5 | 5 |
| −1 outside_confident_region | 7 | 11 | 11 | 20 | 24 | 24 |
| −1 not_on_unique_grch38_node | 21 | 16 | 16 | 92 | 78 | 78 |
| −1 other | 0 | 0 | 0 | 1 | 1 | 1 |
| only SNV tensor | 40 | 39 | 39 | 65 | 53 | 53 |
| no related tensor | 204 | 183 | 183 | 421 | 388 | 388 |
| **有相关 INDEL tensor** | **191** | **213** | **213** | **341** | **386** | **386** |

### I4–I7  (DEL 39, INS 102)

| 最好的相关 tensor 标签 | DEL v5 | DEL v6 | DEL v6new | INS v5 | INS v6 | INS v6new |
|---|---:|---:|---:|---:|---:|---:|
| somatic (1) | 0 | 0 | 0 | 1 | 0 | 0 |
| germline (2) | 0 | 2 | 2 | 4 | 5 | 5 |
| non (0) | 1 | 1 | 1 | 0 | 0 | 0 |
| −1 near_truth_allele_mismatch | 10 | 8 | 8 | 24 | 25 | 25 |
| −1 truth_matches_non_representative_allele | 0 | 0 | 0 | 0 | 0 | 0 |
| −1 outside_confident_region | 2 | 1 | 1 | 4 | 3 | 3 |
| −1 not_on_unique_grch38_node | 1 | 1 | 1 | 3 | 3 | 3 |
| −1 other | 0 | 0 | 0 | 0 | 0 | 0 |
| only SNV tensor | 5 | 5 | 5 | 12 | 14 | 14 |
| no related tensor | 20 | 21 | 21 | 54 | 52 | 52 |
| **有相关 INDEL tensor** | **14** | **13** | **13** | **36** | **36** | **36** |

### 全部 v5 no-candidate INDEL  (DEL 958, INS 1826)

| 最好的相关 tensor 标签 | DEL v5 | DEL v6 | DEL v6new | INS v5 | INS v6 | INS v6new |
|---|---:|---:|---:|---:|---:|---:|
| somatic (1) | 3 | 35 | 35 | 6 | 8 | 8 |
| germline (2) | 107 | 112 | 112 | 124 | 125 | 125 |
| non (0) | 7 | 6 | 5 | 4 | 60 | 4 |
| −1 near_truth_allele_mismatch | 473 | 415 | 416 | 750 | 655 | 711 |
| −1 truth_matches_non_representative_allele | 4 | 40 | 40 | 11 | 30 | 30 |
| −1 outside_confident_region | 35 | 27 | 27 | 62 | 58 | 58 |
| −1 not_on_unique_grch38_node | 30 | 36 | 36 | 223 | 261 | 261 |
| −1 other | 0 | 0 | 0 | 1 | 1 | 1 |
| only SNV tensor | 70 | 72 | 72 | 157 | 136 | 136 |
| no related tensor | 229 | 215 | 215 | 488 | 492 | 492 |
| **有相关 INDEL tensor** | **659** | **671** | **671** | **1181** | **1198** | **1198** |

### I2 + I3：相关 INDEL tensors（去重）按标签

| 标签 / reason | v5 | v6 | v6new |
|---|---:|---:|---:|
| 2 representative_allele_is_germline_truth | 156 | 147 | 147 |
| 1 representative_allele_is_somatic_truth | 4 | 37 | 37 |
| 0 confident_no_truth_allele | 8 | 56 | 6 |
| -1 near_truth_allele_mismatch | 1375 | 877 | 902 |
| -1 not_on_unique_grch38_node | 617 | 416 | 416 |
| -1 outside_confident_region | 111 | 60 | 60 |
| -1 truth_matches_non_representative_allele | 30 | 79 | 79 |

### 全部 tensors：0 和 −1 按 reason


**SNV**

| 标签 | reason | v5 | v6 | v6new |
|---|---|---:|---:|---:|
| 0 | confident_no_truth_allele | 51,228 | 52,347 | 51,228 |
| -1 | germline_truth_filtered_or_outside_bed | 74,076 | 74,076 | 74,076 |
| -1 | near_truth_allele_mismatch | 35,590 | 34,471 | 35,590 |
| -1 | not_on_unique_grch38_node | 106,031 | 106,031 | 106,031 |
| -1 | outside_confident_region | 589,635 | 589,635 | 589,635 |
| -1 | truth_matches_non_representative_allele | 209 | 209 | 209 |
| 1 | representative_allele_is_somatic_truth | 8,234 | 8,234 | 8,234 |
| 2 | representative_allele_is_germline_truth | 315,336 | 315,336 | 315,336 |
| | **合计** | **1,180,339** | **1,180,339** | **1,180,339** |

**INDEL**

| 标签 | reason | v5 | v6 | v6new |
|---|---|---:|---:|---:|
| 0 | confident_no_truth_allele | 562,499 | 1,314,677 | 1,307,535 |
| -1 | germline_truth_filtered_or_outside_bed | 6,462 | 5,522 | 5,522 |
| -1 | near_truth_allele_mismatch | 281,657 | 234,760 | 241,902 |
| -1 | not_on_unique_grch38_node | 149,755 | 123,317 | 123,317 |
| -1 | outside_confident_region | 152,498 | 79,965 | 79,965 |
| -1 | truth_matches_non_representative_allele | 3,031 | 8,300 | 8,300 |
| 1 | representative_allele_is_somatic_truth | 6,777 | 4,097 | 4,097 |
| 2 | representative_allele_is_germline_truth | 165,639 | 146,052 | 146,052 |
| | **合计** | **1,328,318** | **1,916,690** | **1,916,690** |

**INDEL 0 / −1 按类型和长度**

| 标签 | DEL1 v5 | DEL1 v6 | DEL1 v6new | DEL>1 v5 | DEL>1 v6 | DEL>1 v6new | INS1 v5 | INS1 v6 | INS1 v6new | INS>1 v5 | INS>1 v6 | INS>1 v6new |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 492,626 | 1,159,544 | 1,155,344 | 8,698 | 4,211 | 3,671 | 56,636 | 146,473 | 145,449 | 4,539 | 4,449 | 3,071 |
| -1 | 365,642 | 267,287 | 271,487 | 119,595 | 51,840 | 52,380 | 56,067 | 90,681 | 91,705 | 52,099 | 42,056 | 43,434 |
