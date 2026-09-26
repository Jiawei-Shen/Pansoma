
### ONT somatic 召回


**全部 somatic 等位基因**

| 状态 | SNP | DEL | INS |
|---|---:|---:|---:|
| tensor_representative | 8,275 (95.2%) | 1,864 (52.6%) | 1,764 (35.6%) |
| tensor_non_representative_allele | 5 (0.1%) | 446 (12.6%) | 416 (8.4%) |
| filtered | 68 (0.8%) | 510 (14.4%) | 1,042 (21.0%) |
| no_candidate | 342 (3.9%) | 725 (20.5%) | 1,729 (34.9%) |
| 合计 | 8,690 | 3,545 | 4,951 |

**PASS 且在 somatic BED 内**

| 状态 | SNP | DEL | INS |
|---|---:|---:|---:|
| tensor_representative | 8,152 (95.4%) | 1,794 (53.4%) | 1,667 (36.0%) |
| tensor_non_representative_allele | 5 (0.1%) | 412 (12.3%) | 381 (8.2%) |
| filtered | 67 (0.8%) | 475 (14.1%) | 966 (20.9%) |
| no_candidate | 319 (3.7%) | 676 (20.1%) | 1,611 (34.8%) |
| 合计 | 8,543 | 3,357 | 4,625 |

### PacBio v6 → ONT 状态变化

**SNP：PacBio v6 状态（行）→ ONT 状态（列）**

| PacBio \ ONT | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 8209 | 4 | 14 | 7 | 8234 |
| non-rep | 1 | 1 | 0 | 0 | 2 |
| filtered | 17 | 0 | 20 | 17 | 54 |
| no_cand | 48 | 0 | 34 | 318 | 400 |

**DEL：PacBio v6 状态（行）→ ONT 状态（列）**

| PacBio \ ONT | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1704 | 248 | 140 | 10 | 2102 |
| non-rep | 45 | 100 | 69 | 2 | 216 |
| filtered | 79 | 74 | 137 | 47 | 337 |
| no_cand | 36 | 24 | 164 | 666 | 890 |

**INS：PacBio v6 状态（行）→ ONT 状态（列）**

| PacBio \ ONT | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 1539 | 243 | 154 | 23 | 1959 |
| non-rep | 36 | 48 | 72 | 11 | 167 |
| filtered | 146 | 99 | 398 | 78 | 721 |
| no_cand | 43 | 26 | 418 | 1617 | 2104 |

**PacBio 没有候选的 truth，按 PacBio 类别看 ONT 状态**

| PacBio 类别 | 种类 | rep | non-rep | filtered | no_cand |
|---|---|---:|---:|---:|---:|
| S1 | SNP | 40 | 0 | 4 | 4 |
| S2 | SNP | 2 | 0 | 7 | 37 |
| S3 | SNP | 1 | 0 | 2 | 62 |
| S4a | SNP | 0 | 0 | 11 | 59 |
| S4b | SNP | 4 | 0 | 9 | 114 |
| S4c | SNP | 1 | 0 | 1 | 42 |
| I1 | DEL | 29 | 11 | 83 | 292 |
| I2 | DEL | 4 | 12 | 75 | 328 |
| I3 | DEL | 0 | 0 | 3 | 14 |
| I4 | DEL | 0 | 1 | 1 | 14 |
| I5 | DEL | 3 | 0 | 2 | 7 |
| I6 | DEL | 0 | 0 | 0 | 4 |
| I7 | DEL | 0 | 0 | 0 | 7 |
| I1 | INS | 26 | 13 | 233 | 690 |
| I2 | INS | 9 | 10 | 175 | 767 |
| I3 | INS | 1 | 0 | 5 | 77 |
| I4 | INS | 3 | 2 | 1 | 45 |
| I5 | INS | 4 | 1 | 4 | 8 |
| I6 | INS | 0 | 0 | 0 | 11 |
| I7 | INS | 0 | 0 | 0 | 19 |

### ONT 没有候选的 INDEL：原因分类（ONT INS 1,729、DEL 725；PacBio v6 INS 2,104、DEL 890）

| 类别 | INS ONT | INS PacBio | DEL ONT | DEL PacBio | ALT 比例中位数 (ONT) | 按无误 ALT reads 判定 (ONT) | PASS∩BED (ONT) |
|---|---:|---:|---:|---:|---:|---:|---:|
| I1_allele_fully_in_graph | 683 (39.5%) | 962 | 306 (42.2%) | 415 | 0.647 | 747 | 923 |
| I2_partly_in_graph_residual_edit | 879 (50.8%) | 961 | 369 (50.9%) | 419 | 0.656 | 831 | 1166 |
| I3_grch38_path_different_edit | 95 (5.5%) | 83 | 18 (2.5%) | 17 | 0.52 | 39 | 106 |
| I4_no_or_few_alt_reads | 40 (2.3%) | 51 | 20 (2.8%) | 16 | 0.016 | 0 | 54 |
| I5_low_mapq_or_no_reads | 1 (0.1%) | 17 | 0 (0.0%) | 12 | - | 0 | 1 |
| I6_longer_than_50bp | 11 (0.6%) | 11 | 4 (0.6%) | 4 | 0.481 | 10 | 10 |
| I7_unclear | 20 (1.2%) | 19 | 7 (1.0%) | 7 | 0.667 | 8 | 26 |
| I9_truth_edit_on_unbuilt_node | 0 (0.0%) | 0 | 1 (0.1%) | 0 | 0.549 | 1 | 1 |
| 合计 | 1729 | 2104 | 725 | 890 | | | |

**更细的原因（reads 的多数表示）**

| reason | INS | DEL |
|---|---:|---:|
| site_bypassed:branch:edit_on_path_taken | 873 | 116 |
| site_bypassed:branch:no_edit | 679 | 112 |
| site_bypassed:skip_edge:edit_on_path_taken | 6 | 253 |
| site_bypassed:skip_edge:no_edit | 4 | 194 |
| site_on_grch38:other_edit | 95 | 18 |
| fewer_than_3_ALT_reads | 22 | 9 |
| no_ALT_read | 18 | 11 |
| site_on_grch38:no_local_edit | 20 | 7 |
| indel_longer_than_50 | 11 | 4 |
| no_read_spans_site | 1 | 0 |
| key_observed_on_non_target_node | 0 | 1 |

**长度（ONT）**

| | 1 | 2-5 | 6-20 | 21-50 | >50 |
|---|---:|---:|---:|---:|---:|
| INS 有 tensor | 1812 | 323 | 38 | 7 | 0 |
| INS filtered | 484 | 438 | 109 | 11 | 0 |
| INS 没有候选 | 341 | 595 | 593 | 156 | 44 |
| DEL 有 tensor | 1733 | 422 | 128 | 27 | 0 |
| DEL filtered | 235 | 244 | 31 | 0 | 0 |
| DEL 没有候选 | 68 | 299 | 282 | 61 | 15 |

### 残余 edits（≥3 条 ALT reads）的 tensor 标签，按位点计（一个位点可有多个标签）

| 类别 | somatic | germline | non | ignore | no tensor | (无 ≥3 reads 的残余 edit) |
|---|---:|---:|---:|---:|---:|---:|
| DEL I1 | 3 | 20 | 0 | 114 | 35 | 159 |
| DEL I2 | 5 | 34 | 0 | 332 | 101 | 11 |
| DEL I3 | 1 | 2 | 0 | 17 | 5 | 0 |
| DEL I4 | 0 | 0 | 0 | 0 | 0 | 20 |
| DEL I6 | 0 | 0 | 0 | 0 | 0 | 4 |
| DEL I7 | 0 | 0 | 0 | 3 | 0 | 4 |
| DEL I9 | 0 | 0 | 0 | 1 | 0 | 0 |
| INS I1 | 2 | 16 | 0 | 376 | 74 | 272 |
| INS I2 | 3 | 44 | 0 | 726 | 173 | 106 |
| INS I3 | 11 | 18 | 0 | 72 | 20 | 6 |
| INS I4 | 0 | 0 | 0 | 0 | 0 | 40 |
| INS I5 | 0 | 0 | 0 | 0 | 0 | 1 |
| INS I6 | 0 | 0 | 0 | 1 | 1 | 9 |
| INS I7 | 0 | 1 | 0 | 7 | 2 | 11 |

### 对照（有 tensor 的 truth，同样方法分析）

| reason | SNP | INS | DEL |
|---|---:|---:|---:|
| key_observed_on_target_node | 50 | 52 | 100 |
| site_bypassed:branch:edit_on_path_taken | 0 | 15 | 2 |
| site_bypassed:branch:no_edit | 0 | 10 | 1 |
| site_bypassed:skip_edge:no_edit | 0 | 0 | 11 |
| site_bypassed:skip_edge:edit_on_path_taken | 0 | 0 | 8 |
| site_on_grch38:other_edit | 0 | 1 | 0 |

### SNV 没有候选（ONT 342；PacBio v6 400）

| 类别 | ONT | ONT 占比 | ONT ALT 比例中位数 | PacBio v6 |
|---|---:|---:|---:|---:|
| S1_low_mapq_or_no_reads | 1 | 0.3% | - | 48 |
| S2_alt_absent_in_reads | 36 | 10.5% | 0.0 | 46 |
| S3_alt_is_existing_graph_allele | 46 | 13.5% | 0.762 | 65 |
| S4a_repeat_snv_edit_on_branch | 61 | 17.8% | 0.615 | 70 |
| S4b_repeat_absorbed_by_graph_path | 154 | 45.0% | 0.459 | 127 |
| S4c_repeat_no_clear_alt | 44 | 12.9% | 0.0 | 44 |
