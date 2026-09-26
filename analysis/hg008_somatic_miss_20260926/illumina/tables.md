
### Illumina somatic 召回


**全部 somatic 等位基因**

| 状态 | SNP | DEL | INS |
|---|---:|---:|---:|
| tensor_representative | 8,279 (95.3%) | 1,785 (50.4%) | 1,770 (35.8%) |
| tensor_non_representative_allele | 4 (0.0%) | 49 (1.4%) | 29 (0.6%) |
| filtered | 68 (0.8%) | 390 (11.0%) | 519 (10.5%) |
| no_candidate | 339 (3.9%) | 1,321 (37.3%) | 2,633 (53.2%) |
| 合计 | 8,690 | 3,545 | 4,951 |

**PASS 且在 somatic BED 内**

| 状态 | SNP | DEL | INS |
|---|---:|---:|---:|
| tensor_representative | 8,159 (95.5%) | 1,727 (51.4%) | 1,673 (36.2%) |
| tensor_non_representative_allele | 4 (0.0%) | 47 (1.4%) | 29 (0.6%) |
| filtered | 61 (0.7%) | 353 (10.5%) | 476 (10.3%) |
| no_candidate | 319 (3.7%) | 1,230 (36.6%) | 2,447 (52.9%) |
| 合计 | 8,543 | 3,357 | 4,625 |

### ONT / PacBio v6 → Illumina 状态变化

**SNP：ONT 状态（行）→ Illumina 状态（列）**

| ONT \ Illumina | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 8197 | 0 | 14 | 64 | 8275 |
| non-rep | 2 | 2 | 0 | 1 | 5 |
| filtered | 35 | 0 | 8 | 25 | 68 |
| no_cand | 45 | 2 | 46 | 249 | 342 |

**DEL：ONT 状态（行）→ Illumina 状态（列）**

| ONT \ Illumina | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1417 | 13 | 113 | 321 | 1864 |
| non-rep | 217 | 12 | 101 | 116 | 446 |
| filtered | 99 | 14 | 105 | 292 | 510 |
| no_cand | 52 | 10 | 71 | 592 | 725 |

**INS：ONT 状态（行）→ Illumina 状态（列）**

| ONT \ Illumina | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1427 | 10 | 111 | 216 | 1764 |
| non-rep | 223 | 11 | 81 | 101 | 416 |
| filtered | 97 | 6 | 223 | 716 | 1042 |
| no_cand | 23 | 2 | 104 | 1600 | 1729 |

**SNP：PacBio v6 状态（行）→ Illumina 状态（列）**

| PacBio v6 \ Illumina | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 8177 | 1 | 13 | 43 | 8234 |
| non-rep | 0 | 1 | 1 | 0 | 2 |
| filtered | 31 | 0 | 11 | 12 | 54 |
| no_cand | 71 | 2 | 43 | 284 | 400 |

**DEL：PacBio v6 状态（行）→ Illumina 状态（列）**

| PacBio v6 \ Illumina | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1649 | 11 | 135 | 307 | 2102 |
| non-rep | 50 | 22 | 91 | 53 | 216 |
| filtered | 29 | 6 | 66 | 236 | 337 |
| no_cand | 57 | 10 | 98 | 725 | 890 |

**INS：PacBio v6 状态（行）→ Illumina 状态（列）**

| PacBio v6 \ Illumina | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1712 | 12 | 96 | 139 | 1959 |
| non-rep | 17 | 10 | 75 | 65 | 167 |
| filtered | 17 | 4 | 218 | 482 | 721 |
| no_cand | 24 | 3 | 130 | 1947 | 2104 |

**ONT 没有候选的 truth，按 ONT 类别看 Illumina 状态**

| ONT 类别 | 种类 | rep | non-rep | filtered | no_cand |
|---|---|---:|---:|---:|---:|
| S1 | SNP | 0 | 0 | 0 | 1 |
| S2 | SNP | 2 | 0 | 3 | 31 |
| S3 | SNP | 8 | 0 | 6 | 32 |
| S4a | SNP | 11 | 0 | 8 | 42 |
| S4b | SNP | 20 | 2 | 25 | 107 |
| S4c | SNP | 4 | 0 | 4 | 36 |
| I1 | DEL | 7 | 0 | 14 | 285 |
| I2 | DEL | 40 | 10 | 53 | 266 |
| I3 | DEL | 3 | 0 | 3 | 12 |
| I4 | DEL | 2 | 0 | 0 | 18 |
| I6 | DEL | 0 | 0 | 0 | 4 |
| I7 | DEL | 0 | 0 | 1 | 6 |
| I9 | DEL | 0 | 0 | 0 | 1 |
| I1 | INS | 4 | 0 | 19 | 660 |
| I2 | INS | 9 | 1 | 71 | 798 |
| I3 | INS | 9 | 1 | 11 | 74 |
| I4 | INS | 0 | 0 | 2 | 38 |
| I5 | INS | 0 | 0 | 0 | 1 |
| I6 | INS | 0 | 0 | 0 | 11 |
| I7 | INS | 1 | 0 | 1 | 18 |

**PacBio v6 没有候选的 truth，按 PacBio v6 类别看 Illumina 状态**

| PacBio v6 类别 | 种类 | rep | non-rep | filtered | no_cand |
|---|---|---:|---:|---:|---:|
| S1 | SNP | 22 | 0 | 0 | 26 |
| S2 | SNP | 4 | 0 | 6 | 36 |
| S3 | SNP | 12 | 1 | 5 | 47 |
| S4a | SNP | 12 | 0 | 10 | 48 |
| S4b | SNP | 18 | 0 | 15 | 94 |
| S4c | SNP | 3 | 1 | 7 | 33 |
| I1 | DEL | 5 | 0 | 17 | 393 |
| I2 | DEL | 48 | 10 | 75 | 286 |
| I3 | DEL | 0 | 0 | 3 | 14 |
| I4 | DEL | 2 | 0 | 0 | 14 |
| I5 | DEL | 2 | 0 | 2 | 8 |
| I6 | DEL | 0 | 0 | 0 | 4 |
| I7 | DEL | 0 | 0 | 1 | 6 |
| I1 | INS | 4 | 0 | 37 | 921 |
| I2 | INS | 14 | 2 | 85 | 860 |
| I3 | INS | 0 | 1 | 5 | 77 |
| I4 | INS | 2 | 0 | 1 | 48 |
| I5 | INS | 3 | 0 | 1 | 13 |
| I6 | INS | 0 | 0 | 0 | 11 |
| I7 | INS | 1 | 0 | 1 | 17 |

### Illumina 没有候选的 INDEL：原因分类（Illumina INS 2,633、DEL 1,321；ONT INS 1,729、DEL 725；PacBio v6 INS 2,104、DEL 890）

| 类别 | INS Illumina | INS ONT | INS PacBio | DEL Illumina | DEL ONT | DEL PacBio | ALT 比例中位数 (Illumina) | 按无误 ALT reads 判定 (Illumina) | PASS∩BED (Illumina) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| I1_allele_fully_in_graph | 1136 (43.1%) | 683 | 962 | 780 (59.0%) | 306 | 415 | 0.889 | 1610 | 1790 |
| I2_partly_in_graph_residual_edit | 1056 (40.1%) | 879 | 961 | 448 (33.9%) | 369 | 419 | 0.953 | 803 | 1397 |
| I3_grch38_path_different_edit | 81 (3.1%) | 95 | 83 | 25 (1.9%) | 18 | 17 | 0.628 | 43 | 102 |
| I4_no_or_few_alt_reads | 264 (10.0%) | 40 | 51 | 42 (3.2%) | 20 | 16 | 0.017 | 0 | 277 |
| I5_low_mapq_or_no_reads | 84 (3.2%) | 1 | 17 | 17 (1.3%) | 0 | 12 | - | 0 | 91 |
| I6_longer_than_50bp | 0 (0.0%) | 11 | 11 | 1 (0.1%) | 4 | 4 | 0.667 | 1 | 0 |
| I7_unclear | 12 (0.5%) | 20 | 19 | 8 (0.6%) | 7 | 7 | 1.0 | 8 | 20 |
| 合计 | 2633 | 1729 | 2104 | 1321 | 725 | 890 | | | |

**更细的原因（reads 的多数表示）**

| reason | INS | DEL |
|---|---:|---:|
| site_bypassed:branch:no_edit | 1129 | 104 |
| site_bypassed:branch:edit_on_path_taken | 1051 | 121 |
| site_bypassed:skip_edge:no_edit | 7 | 676 |
| site_bypassed:skip_edge:edit_on_path_taken | 5 | 327 |
| fewer_than_3_ALT_reads | 135 | 24 |
| no_ALT_read | 129 | 18 |
| site_on_grch38:other_edit | 81 | 25 |
| no_read_spans_site | 84 | 17 |
| site_on_grch38:no_local_edit | 12 | 8 |
| indel_longer_than_50 | 0 | 1 |

**长度（Illumina）**

| | 1 | 2-5 | 6-20 | 21-50 | >50 |
|---|---:|---:|---:|---:|---:|
| INS 有 tensor | 1604 | 169 | 26 | 0 | 0 |
| INS filtered | 204 | 264 | 44 | 7 | 0 |
| INS 没有候选 | 829 | 923 | 670 | 167 | 44 |
| DEL 有 tensor | 1354 | 306 | 147 | 27 | 0 |
| DEL filtered | 153 | 200 | 35 | 2 | 0 |
| DEL 没有候选 | 529 | 459 | 259 | 59 | 15 |

### 残余 edits（≥3 条 ALT reads）的 tensor 标签，按位点计（一个位点可有多个标签）

| 类别 | somatic | germline | non | ignore | no tensor | (无 ≥3 reads 的残余 edit) |
|---|---:|---:|---:|---:|---:|---:|
| DEL I1 | 1 | 3 | 0 | 242 | 179 | 439 |
| DEL I2 | 7 | 102 | 0 | 343 | 165 | 29 |
| DEL I3 | 1 | 6 | 0 | 22 | 7 | 1 |
| DEL I4 | 0 | 0 | 0 | 0 | 0 | 42 |
| DEL I5 | 0 | 0 | 0 | 0 | 0 | 17 |
| DEL I6 | 0 | 0 | 0 | 0 | 0 | 1 |
| DEL I7 | 0 | 0 | 0 | 1 | 1 | 6 |
| INS I1 | 0 | 3 | 0 | 579 | 364 | 422 |
| INS I2 | 3 | 122 | 0 | 862 | 399 | 92 |
| INS I3 | 13 | 17 | 0 | 57 | 21 | 5 |
| INS I4 | 0 | 0 | 0 | 0 | 0 | 264 |
| INS I5 | 0 | 0 | 0 | 0 | 0 | 84 |
| INS I7 | 0 | 0 | 0 | 4 | 3 | 6 |

### 对照（有 tensor 的 truth，同样方法分析）

| reason | SNP | INS | DEL |
|---|---:|---:|---:|
| key_observed_on_target_node | 50 | 54 | 109 |
| site_bypassed:skip_edge:edit_on_path_taken | 0 | 2 | 16 |
| site_bypassed:branch:edit_on_path_taken | 0 | 4 | 5 |
| site_on_grch38:other_edit | 0 | 1 | 4 |
| fewer_than_3_ALT_reads | 0 | 2 | 0 |
| site_bypassed:branch:no_edit | 0 | 1 | 1 |
| no_read_spans_site | 0 | 1 | 0 |

### SNV 没有候选（Illumina 339；ONT 342；PacBio v6 400）

| 类别 | Illumina | Illumina 占比 | Illumina ALT 比例中位数 | ONT | PacBio v6 |
|---|---:|---:|---:|---:|---:|
| S1_low_mapq_or_no_reads | 55 | 16.2% | 0.0 | 1 | 48 |
| S2_alt_absent_in_reads | 40 | 11.8% | 0.0 | 36 | 46 |
| S3_alt_is_existing_graph_allele | 76 | 22.4% | 0.913 | 46 | 65 |
| S4a_repeat_snv_edit_on_branch | 57 | 16.8% | 0.818 | 61 | 70 |
| S4b_repeat_absorbed_by_graph_path | 65 | 19.2% | 0.621 | 154 | 127 |
| S4c_repeat_no_clear_alt | 46 | 13.6% | 0.0 | 44 | 44 |
