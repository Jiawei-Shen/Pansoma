# HG008 ONT-UL：没有候选的 somatic 真值，原因分析（和 PacBio v6 对比）

- 数据：HG008-T ONT-UL R10.4.1 dorado 0.8.1 sup，54×。tensors 来自 `indexed_gam_pipeline_v3`，job 367607/367608（常染色体，format v6：indel 左对齐、site 单元；**没有 supplement 轮**）。目录是 `pansoma_v2_tensors/Liss_lab_Northeastern-ONT-UL-20241216/v3_tensors`。
- 标签：v3 冻结源码的 `truth_labels`，已包含 INS 等价范围修复，和 PacBio v6 现在的规则相同。
- 真值：somatic `HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz` + `_all.bed`；germline dipcall `dip.vcf.gz` + `dip.bed`。
- 对象：ONT 召回报告里"没有候选"的 somatic 等位基因，共 **2,796** 个（SNV 342，INS 1,729，DEL 725）。SNV 和 INDEL 都按 reads 逐个分析了。
- 对照：随机抽 250 个有 ONT tensor 的 truth（INDEL ≥2 bp 150 个、1 bp 50 个，SNV 50 个）。
- 方法和 PacBio v6 分析相同（`tmp/somatic_miss_analysis_v6/REPORT.md` 第 8 节）。解码用 v3 冻结源码的 Python 解码器（build 用的是结果相同的 native 版本），每条 record 只保留位点附近 ±600 个 read 碱基。在 24 个 ONT 位点上和整条 record 解码比较过，结果完全相同，0 差异。SNV 的 S1–S4c 分类规则是从 v5 SNV 分析的中间结果反推出来的，在 v5 的 454 个 SNV 位点上和原分类 100% 一致。

## 1. 结论

1. **主要原因和 PacBio 一样：somatic allele 已经全部或部分在图里。** INDEL 未命中里 I1 + I2 占 INS 的 **90.3%**（1,562/1,729）、DEL 的 **93.1%**（675/725）。SNV 里 S3 + S4a + S4b 占 76.3%（261/342）。
2. **ONT 的 tensor 覆盖和 PacBio 差不多，但代表 allele 更少。**
   - truth 有 tensor（代表或非代表 allele）的比例：DEL 65.2%（PacBio 65.4%），INS 44.0%（PacBio 42.9%）。
   - truth 是代表 allele 的比例：DEL 52.6%（PacBio 59.3%），INS 35.6%（PacBio 39.6%）。
   - 原因是 ONT 的噪声 allele 更多：PacBio 里是代表 allele 的，有 248 个 DEL、243 个 INS 在 ONT 里变成了同一 site 的非代表 allele。
3. **ONT 的"没有候选"更少，"filtered"更多。** PacBio 没有候选的 truth 里，ONT 有 418 个 INS、164 个 DEL 变成了 filtered，reads 写出了 truth 的 edit，但数量不够。ONT filtered 的原因大多是 min_variants（INS 831，DEL 397）：ONT 的 reads 把同一个 allele 分散成多种写法，每种都不到 3 条。
4. **ONT-UL 解决了低 MAPQ 区域。** I5 只有 1 个（PacBio 29），SNV 的 S1 只有 1 个（PacBio 48）。PacBio S1 的 48 个 SNV，在 ONT 里有 40 个成为代表 allele：ultra-long reads 能唯一比对到片段重复区。
5. **I2 的残余 edit 大多有 tensor，但标签不对。** INS I2 879 个位点里，529 个的残余 edit tensor 在分支节点上（−1）。残余 edit tensor 被标成 **germline** 的有 **194 个不同候选**（INS 位点 79，DEL 56，SNV 27）。
6. **没有 supplement 的代价很小**：只有 1 个 DEL（chr11:134897366）是 reads 写出了 truth edit，但它所在节点不是 target（I9）。
7. SNV：35 个 somatic SNV 的 germline 真值在同一位置有同一个 ALT（PacBio 分析里是 38 个），不适合当 somatic 正例。

## 2. 召回总览

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

## 3. PacBio v6 → ONT 的变化

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

## 4. INDEL 原因分类（ONT）

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
| key_observed_on_non_target_node | 0 | 1 |
| no_read_spans_site | 1 | 0 |

**长度（ONT）**

| | 1 | 2-5 | 6-20 | 21-50 | >50 |
|---|---:|---:|---:|---:|---:|
| INS 有 tensor | 1812 | 323 | 38 | 7 | 0 |
| INS filtered | 484 | 438 | 109 | 11 | 0 |
| INS 没有候选 | 341 | 595 | 593 | 156 | 44 |
| DEL 有 tensor | 1733 | 422 | 128 | 27 | 0 |
| DEL filtered | 235 | 244 | 31 | 0 | 0 |
| DEL 没有候选 | 68 | 299 | 282 | 61 | 15 |

## 5. 残余 edits 的 tensor 标签

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

逐条列表：`residual_edits_labelled_germline.tsv`（204 行，194 个不同候选）。

## 6. SNV（342）

### SNV 没有候选（ONT 342；PacBio v6 400）

| 类别 | ONT | ONT 占比 | ONT ALT 比例中位数 | PacBio v6 |
|---|---:|---:|---:|---:|
| S1_low_mapq_or_no_reads | 1 | 0.3% | - | 48 |
| S2_alt_absent_in_reads | 36 | 10.5% | 0.0 | 46 |
| S3_alt_is_existing_graph_allele | 46 | 13.5% | 0.762 | 65 |
| S4a_repeat_snv_edit_on_branch | 61 | 17.8% | 0.615 | 70 |
| S4b_repeat_absorbed_by_graph_path | 154 | 45.0% | 0.459 | 127 |
| S4c_repeat_no_clear_alt | 44 | 12.9% | 0.0 | 44 |

**和命中 SNV 的对比**（`snv_compare.py`；命中 8,275 个，未命中 342 个）

| 特征 | 命中 | 没有候选 |
|---|---:|---:|
| SNV 所在的 GRCh38 节点只有 1 bp | 2.6% | **70.8%** |
| 节点 2–32 bp | 4.3% | 21.6% |
| 节点 >32 bp | 93.1% | 7.6% |
| 50 bp 内有 ≥2 个 germline 变异 | 3.7% | **76.9%** |
| 50 bp 内没有 germline 变异 | 83.7% | 7.6% |
| germline 真值在同一位置、同一个 ALT | 0.1% | **10.2%（35 个）** |

## 7. 例子

每行都用 `show_locus.py` 把 reads 逐条打出来看过，完整输出在 `show_examples.txt`。

| 位点 | 真值 | 类别 | reads 实际的样子 |
|---|---|---|---|
| chr15:83610975 | C>CAAAA | I1 | 48 条 ALT reads 走非 GRCh38 节点 16745890→88→87→85，全程没有 edits |
| chr1:237212078 | CT>C | I1 | ALT reads 走跳过节点 4377655 的边（4377654→4377656），没有 edits |
| chr13:65370940 | C>CATATAT | I2 | ALT reads 走 AT 重复分支 14478371/14478370，再加一个 2 bp 插入（残余 edit） |
| chr3:121747150 | 26 bp TA 重复缺失 | I2 | ALT reads 走跳过边 35323394→35323385，再 `DEL` 12 bp（残余 edit） |
| chr11:105058573 | +TATCTATCTATCTA | I3 | ALT reads 走一个 TATC 重复分支，再在 GRCh38 节点上加 2 bp 插入，写法和 truth 不同 |
| chr11:134897366 | TTTG>T | I9 | 22 条 ALT reads 写出 truth 的 `DEL`，但所在节点不是 target（ONT 没有 supplement 轮） |
| chr13:76153195 | G>A | S3 | 52 条 ALT reads 走 SNP 分叉的 A 节点 14096196，没有错配 |
| chr12:44973851 | T>G | S4b | TA/GA 重复交界：ALT reads 走现成的非 GRCh38 节点 10669552，没有 edits |
| chr20:52733809 | A>C | S4a | CT/AT 重复：reads 的 CT 长度各不相同，走多条分支，SNV 成了别处的错配（68 条 reads 里只有 1 条完全等于 ALT） |
| chr8:55438589 | T>A | S2 | polyA 前：reads 都是 `CA` 加长短不一的 polyA，没有 read 带 ALT |
| chr17:26640943 | T>A | S1 | 没有 MAPQ>10 的 read 跨过这个位点（ONT 唯一的 S1） |

## 8. 对照

### 对照（有 tensor 的 truth，同样方法分析）

| reason | SNP | INS | DEL |
|---|---:|---:|---:|
| key_observed_on_target_node | 50 | 52 | 100 |
| site_bypassed:branch:edit_on_path_taken | 0 | 15 | 2 |
| site_bypassed:skip_edge:no_edit | 0 | 0 | 11 |
| site_bypassed:branch:no_edit | 0 | 10 | 1 |
| site_bypassed:skip_edge:edit_on_path_taken | 0 | 0 | 8 |
| site_on_grch38:other_edit | 0 | 1 | 0 |

对照组里，SNV 50/50、INS 52/78、DEL 100/122 是"truth 键出现在 target 节点上"。其余对照是多数 reads 绕开 GRCh38，但仍有足够的 reads 写出 truth 的 edit。

## 9. 建议

1. **标签**：和 PacBio 一样。tensor 落在 somatic 真值的等价范围内时，不标 germline 而标 −1（ONT 194 个）。用单倍型投影来匹配真值，可以同时解决 I2/I3 和分支节点上的残余 edit。
2. **ONT 的 min_variants / AF**：ONT 有 1,552 个 somatic INDEL 是 filtered，比 PacBio 的 1,058 多很多，主要卡在 min_variants。如果单倍型投影把同一 allele 的不同写法合并计数，这部分会大幅恢复。
3. **代表 allele**：ONT 的噪声 allele 常常抢走代表位置（非代表 allele DEL 446、INS 416）。在训练里可以考虑用"site 中任一 allele 匹配 truth"来打标签。
4. **supplement**：ONT 只丢了 1 个 truth，可以继续不开。

## 10. 文件

都在 `tmp/somatic_miss_analysis_ont/`：

| 文件 | 内容 |
|---|---|
| `REPORT.md` | 本报告 |
| `indel_no_candidate.tsv` | 2,454 个 INDEL 未命中，每个一行（列同 PacBio v6） |
| `snv_no_candidate.tsv` | 342 个 SNV 未命中，另有 SNV 专用列：等长 reads 数、SNV 位置的碱基分布、ALT reads 走的节点、节点长度、是否 target、discovery 统计 |
| `controls.tsv` | 250 个对照 |
| `transitions.tsv` | 每个 truth 的 PacBio v6 状态/类别 → ONT 状态/类别 |
| `residual_edits_labelled_germline.tsv` | 残余 edit tensor 被标成 germline 的 204 行 |
| `show_examples.txt`、`examples_pick.tsv`、`show/` | 例子的 reads 逐条显示 |
| `tables.md` | 全部表格（`tables.py` 生成） |
| `reads.py`、`snv_detail.py`、`snv_nodes.py`、`snv_compare.py`、`report.py`、`tables.py`、`show_locus.py`、`validate_trim.py` | 分析脚本 |
| `reads_part*.jsonl`、`snv_detail.jsonl`、`snv_nodes.json` | 原始结果 |
