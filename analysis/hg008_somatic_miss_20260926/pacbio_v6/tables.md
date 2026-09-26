
### v6 somatic 召回（新标签规则；召回与规则无关）


**全部 somatic 等位基因**

| 状态 | SNP | DEL | INS |
|---|---:|---:|---:|
| tensor_representative | 8,234 (94.8%) | 2,102 (59.3%) | 1,959 (39.6%) |
| tensor_non_representative_allele | 2 (0.0%) | 216 (6.1%) | 167 (3.4%) |
| filtered | 54 (0.6%) | 337 (9.5%) | 721 (14.6%) |
| no_candidate | 400 (4.6%) | 890 (25.1%) | 2,104 (42.5%) |
| 合计 | 8,690 | 3,545 | 4,951 |

**PASS 且在 somatic BED 内**

| 状态 | SNP | DEL | INS |
|---|---:|---:|---:|
| tensor_representative | 8,115 (95.0%) | 2,019 (60.1%) | 1,850 (40.0%) |
| tensor_non_representative_allele | 2 (0.0%) | 196 (5.8%) | 156 (3.4%) |
| filtered | 51 (0.6%) | 310 (9.2%) | 653 (14.1%) |
| no_candidate | 375 (4.4%) | 832 (24.8%) | 1,966 (42.5%) |
| 合计 | 8,543 | 3,357 | 4,625 |

### v5 → v6 状态变化（INDEL）

**DEL：v5 状态（行）→ v6 状态（列）**

| v5 \ v6 | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1936 | 88 | 10 | 9 | 2043 |
| non-rep | 9 | 26 | 2 | 0 | 37 |
| filtered | 123 | 85 | 269 | 30 | 507 |
| no_cand | 34 | 17 | 56 | 851 | 958 |

**INS：v5 状态（行）→ v6 状态（列）**

| v5 \ v6 | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1875 | 29 | 4 | 15 | 1923 |
| non-rep | 4 | 92 | 0 | 12 | 108 |
| filtered | 78 | 46 | 695 | 275 | 1094 |
| no_cand | 2 | 0 | 22 | 1802 | 1826 |

**v5 没有候选的 INDEL，按 v5 类别看 v6 状态**

| v5 类别 | DEL rep | DEL non-rep | DEL filtered | DEL no_cand | INS rep | INS non-rep | INS filtered | INS no_cand |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| I1 | 0 | 0 | 22 | 413 | 2 | 0 | 12 | 813 |
| I2 | 3 | 14 | 31 | 381 | 0 | 0 | 10 | 807 |
| I3 | 31 | 3 | 2 | 19 | 0 | 0 | 0 | 80 |
| I4 | 0 | 0 | 1 | 16 | 0 | 0 | 0 | 51 |
| I5 | 0 | 0 | 0 | 12 | 0 | 0 | 0 | 17 |
| I6 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 10 |
| I7 | 0 | 0 | 0 | 8 | 0 | 0 | 0 | 24 |
| 合计 | 34 | 17 | 56 | 851 | 2 | 0 | 22 | 1802 |

### v6 没有候选的 INDEL：原因分类（v6 INS 2,104、DEL 890；v5 INS 1,826、DEL 958）

| 类别 | INS v6 | INS v5 | DEL v6 | DEL v5 | ALT 比例中位数 (v6) | 按无误 ALT reads 判定 (v6) | PASS∩BED (v6) |
|---|---:|---:|---:|---:|---:|---:|---:|
| I1_allele_fully_in_graph | 962 (45.7%) | 827 | 415 (46.6%) | 435 | 0.778 | 1089 | 1288 |
| I2_partly_in_graph_residual_edit | 961 (45.7%) | 817 | 419 (47.1%) | 429 | 0.824 | 1040 | 1293 |
| I3_grch38_path_different_edit | 83 (3.9%) | 80 | 17 (1.9%) | 55 | 0.589 | 48 | 93 |
| I4_no_or_few_alt_reads | 51 (2.4%) | 51 | 16 (1.8%) | 17 | 0.0 | 0 | 60 |
| I5_low_mapq_or_no_reads | 17 (0.8%) | 17 | 12 (1.3%) | 12 | - | 0 | 28 |
| I6_longer_than_50bp | 11 (0.5%) | 10 | 4 (0.4%) | 2 | 0.53 | 12 | 10 |
| I7_unclear | 19 (0.9%) | 24 | 7 (0.8%) | 8 | 0.689 | 8 | 26 |
| 合计 | 2104 | 1826 | 890 | 958 | | | |

**更细的原因（reads 的多数表示）**

| reason | INS | DEL |
|---|---:|---:|
| site_bypassed:branch:no_edit | 957 | 113 |
| site_bypassed:branch:edit_on_path_taken | 954 | 114 |
| site_bypassed:skip_edge:edit_on_path_taken | 7 | 305 |
| site_bypassed:skip_edge:no_edit | 5 | 302 |
| site_on_grch38:other_edit | 83 | 17 |
| no_ALT_read | 35 | 11 |
| no_read_spans_site | 17 | 12 |
| site_on_grch38:no_local_edit | 19 | 7 |
| fewer_than_3_ALT_reads | 16 | 5 |
| indel_longer_than_50 | 11 | 4 |

**长度（v6）**

| | 1 | 2-5 | 6-20 | 21-50 | >50 |
|---|---:|---:|---:|---:|---:|
| INS 有 tensor | 1787 | 262 | 67 | 10 | 0 |
| INS filtered | 372 | 301 | 40 | 8 | 0 |
| INS 没有候选 | 478 | 793 | 633 | 156 | 44 |
| DEL 有 tensor | 1768 | 382 | 141 | 27 | 0 |
| DEL filtered | 143 | 183 | 11 | 0 | 0 |
| DEL 没有候选 | 125 | 400 | 289 | 61 | 15 |

### 残余 edits（≥3 条 ALT reads）的 tensor 标签，按位点计（一个位点可有多个标签）

| 类别 | somatic | germline | non | ignore | no tensor | (无 ≥3 reads 的残余 edit) |
|---|---:|---:|---:|---:|---:|---:|
| DEL I1 | 2 | 28 | 0 | 137 | 132 | 188 |
| DEL I2 | 8 | 62 | 0 | 369 | 185 | 6 |
| DEL I3 | 0 | 5 | 0 | 17 | 5 | 0 |
| DEL I4 | 0 | 0 | 0 | 0 | 0 | 16 |
| DEL I5 | 0 | 0 | 0 | 0 | 0 | 12 |
| DEL I6 | 0 | 0 | 0 | 0 | 0 | 4 |
| DEL I7 | 0 | 0 | 0 | 6 | 4 | 1 |
| INS I1 | 3 | 34 | 0 | 416 | 299 | 352 |
| INS I2 | 4 | 73 | 0 | 918 | 534 | 6 |
| INS I3 | 11 | 14 | 0 | 67 | 33 | 0 |
| INS I4 | 0 | 0 | 0 | 0 | 0 | 51 |
| INS I5 | 0 | 0 | 0 | 0 | 0 | 17 |
| INS I6 | 0 | 0 | 0 | 3 | 0 | 8 |
| INS I7 | 0 | 2 | 0 | 8 | 2 | 10 |

### 对照（有 tensor 的 INDEL truth，同样方法分析）

| reason | INS | DEL |
|---|---:|---:|
| key_observed_on_target_node | 67 | 114 |
| site_bypassed:skip_edge:edit_on_path_taken | 1 | 5 |
| site_bypassed:branch:no_edit | 3 | 2 |
| site_bypassed:skip_edge:no_edit | 0 | 4 |
| site_bypassed:branch:edit_on_path_taken | 2 | 0 |
| site_on_grch38:other_edit | 2 | 0 |

### SNV 没有候选（400；分类沿用 v5 分析，SNV tensors 两版相同）

| 类别 | 个数 | 占比 |
|---|---:|---:|
| S1_low_mapq_or_no_reads | 48 | 12.0% |
| S2_alt_absent_in_reads | 46 | 11.5% |
| S3_alt_is_existing_graph_allele | 65 | 16.2% |
| S4a_repeat_snv_edit_on_branch | 70 | 17.5% |
| S4b_repeat_absorbed_by_graph_path | 127 | 31.8% |
| S4c_repeat_no_clear_alt | 44 | 11.0% |
