# HG008 PacBio v6：没有候选的 somatic 真值，原因分析

- 数据：HG008-T PacBio Revio 116×，**v6 tensors**（job 364673/364850，常染色体，format v6：indel 左对齐、site 单元、supplement runs）。
- 标签：2026-09-25 用当前规则（含 INS 等价范围修复，commit `9b8b1b3`）重打过的 `v6_tensors` 标签。召回状态和标签规则无关。
- 真值：somatic `HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz` + `_all.bed`；germline dipcall `dip.vcf.gz` + `dip.bed`。
- 对象：v6 召回报告里状态为"没有候选"的 somatic 等位基因，共 **3,394** 个（SNV 400，INS 2,104，DEL 890）。"没有候选"指：它的任何一种节点坐标写法，既不是任何 tensor 的 allele，也不在 `filtered_candidates` 里。
- 对照：随机抽 200 个 v6 有 tensor 的 somatic INDEL（≥2 bp 150 个，1 bp 50 个），用同样的方法分析。
- **SNV 部分沿用 v5 分析**：v5 和 v6 的 SNV tensors 完全相同。v5 的 404 个 SNV 未命中里，有 4 个在 v6 变成了 filtered（有候选但没过滤），其余 400 个的分类不变。**INDEL 部分全部用 v6 重新做了**。

## 1. 结论

1. **主要原因没变：somatic INDEL 在泛基因组图里已经全部或部分存在。** 候选只来自 reads 相对它所走节点的 edits。图里已经有这个 allele 时，ALT reads 走现成的分支或跳过边，要么完全没有 edits（I1），要么只剩一个写法不同的残余 edit（I2）。v6 里 **I1 + I2 占 INS 的 91.4%（1,923/2,104）、DEL 的 93.7%（834/890）**，大多数是按完全无测序错误的 ALT reads 判定的。
2. **INS 未命中比 v5 多了 278 个（1,826 → 2,104），DEL 少了 68 个（958 → 890）。**
   - DEL 的改善主要来自 v6 合并跨节点缺失：v5 的 I3 DEL 55 个里 31 个在 v6 成为代表 allele，I3 DEL 降到 17 个。
   - INS 变多，是因为 v6 的左对齐把一部分 edit 推进了分支节点：
     - v5 里有 tensor、v6 里没有候选的 truth 共 **36 个**（INS 27，DEL 9）；
     - v5 里 filtered、v6 里没有候选的 **305 个**（INS 275，DEL 30）。
     这些全部是 I1/I2，其中 330/341 的序列背景里有 ≥8 bp 的 homopolymer（全部 INDEL 未命中是 2,221/2,994）。ALT reads 走一条人群插入分支（非 GRCh38 节点），剩下的 +A/+T 被左对齐到重复区起点，正好落在这个分支节点上。于是候选没有 GRCh38 坐标（标签 −1 not_on_unique_grch38_node），也不等于任何 truth 键。例子见第 6 节的 chr1:15398784、chr1:81119123。
3. **I2 的残余 edit 大多有 tensor，但标签不对。** I2 INS 961 个位点里，658 个的残余 edit tensor 在分支节点上（−1）。按 ≥3 条 ALT reads 统计，残余 edit tensor 被标成 **germline** 的有 **301 个不同候选**（INS 位点 123，DEL 位点 95），这些是标签错误；大部分是 −1，不参与训练。
4. **次要原因**（和 v5 基本一样）：
   - ALT reads 很少或没有（I4：INS 51，DEL 16）；
   - MAPQ≤10 或没有 reads（I5：29）；
   - reads 里的 indel 超过 50 bp（I6：15）；
   - 不确定（I7：26）。
   v5 的 I3 INS（80 个，edits 写法不同）在 v6 还是 83 个，一个都没恢复。
5. **没有"truth edit 出现了却没建"的情况**：3,394 个未命中里，没有一个位点的 ALT reads 产生了 truth 键的候选。v6 的 target 集合（常染色体 parts + 3 轮 supplement）覆盖了所有 reads 写出 truth 键的节点。

## 2. 召回总览

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

## 3. v5 → v6 的变化（INDEL）

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

要点：
- DEL：v5 没有候选的 958 个里，34 个变成代表 allele、17 个变成非代表 allele、56 个变成 filtered。主要来自 I3（跨节点缺失合并）和 I2。
- INS：v5 没有候选的 1,826 个里，只有 2 个变成代表 allele。反方向有 275 个从 filtered 变成没有候选，27 个从有 tensor 变成没有候选，原因见第 1 节第 2 条。

## 4. INDEL 原因分类（v6）

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

## 5. 残余 edits 的 tensor 标签

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

- 标签的含义：somatic 基本是附近另一个 somatic truth 的 tensor；germline 是错标；ignore 是 −1（near_truth_allele_mismatch、not_on_unique_grch38_node 等）；no tensor 是残余 edit 没有建成 tensor（AF 或 reads 数不够）。
- 残余 edit tensor 被标成 germline 的逐条列表：`residual_edits_labelled_germline.tsv`（301 行）。

## 6. 例子

每行都用 `show_locus.py` 把 reads 的局部序列和图路径逐条打出来看过，完整输出在 `show_examples.txt`。

| 位点 | 真值 | 类别 | reads 实际的样子 |
|---|---|---|---|
| chr16:52827238 | C>CA | I1 | 103 条 ALT reads 走多一个 A 的非 GRCh38 节点 19238194，全程没有 edits |
| chr12:25045848 | TAA>T | I1 | 99 条无误 ALT reads 走跳过节点 10524294 的边（10524295→10524293），没有 edits |
| chr20:58872348 | T>TAA | I2 | ALT reads 走 +A 分支 31293240，再在重复区起点 +A；残余 edit `31293240:0:INS:>A`（85 reads）在分支节点上，标签 −1 |
| chr2:80669853 | ATT>A | I2 | ALT reads 走跳过一个 T 的边，再 `DEL T`；残余 edit `26920655:0:DEL:T>`（80 reads）被标成 **germline** |
| chr13:107514434 | T>TATCA | I3 | ATCT 重复：ALT reads 走多一个重复单元的分支，再带一个 T>A 错配；候选是 SNP `13612951:0:SNP:T>A`（72 reads），标签 −1 |
| chr15:34847140 | CAAA>C | I3 | reads 跳过 AAA，剩下的差异写成一个 A>C 错配（56 reads 的 SNP 候选），标签 −1 |
| chr1:15398784 | C>CAA | I2（v5 有 tensor） | ALT reads 走 +A 分支 2659415，再 +A；v6 把这个 +A 左对齐到分支节点上，候选 `2659415:0:INS:>A`（71 reads）没有 GRCh38 坐标，标签 −1 |
| chr1:81119123 | G>GTT | I2（v5 有 tensor） | 同上：+T 分支 1207140，再 +T，左对齐后落在分支节点上 |
| chr1:21780738 | A>AT | I1（v5 是 filtered） | ALT reads 走 +T 分支 2824685，没有 edits；v5 里少数 reads 在 GRCh38 节点上写出的 truth 键，v6 左对齐后到了分支节点上 |

## 7. SNV（400；沿用 v5 分析）

### 7.1 和命中的 SNV 对比（`snv_compare.py`）

| 特征 | 命中（8,234） | 没有候选（v5 的 404） |
|---|---:|---:|
| SNV 所在的 GRCh38 节点只有 1 bp（图里在这个位置本来就有 SNP 分叉） | 2.7% | **61.6%** |
| 节点 2–32 bp | 4.2% | 21.0% |
| 节点 >32 bp | 93.1% | 17.3% |
| 50 bp 内有 ≥2 个 germline 变异 | 3.7% | **70.5%** |
| 50 bp 内没有 germline 变异 | 83.8% | 15.1% |
| 50 bp 内还有别的 somatic 真值 | 7% | 58% |
| germline 真值在同一位置、同一个 ALT | 0.1% | **9.4%（38 个）** |

### 7.2 原因分类

| 类别 | 个数 | 占比 | ALT 比例中位数 | 含义 |
|---|---:|---:|---:|---|
| **S3** ALT 是图里现成的等位基因 | 65 | 16.2% | 0.84 | 序列背景简单（reads 与 REF 等长），带 ALT 的 reads 在这个位置走非 GRCh38 的节点（1 bp SNP 分叉的另一个碱基），没有错配。ALT 比例很高，很多接近 1（例如 chr13:87814362，121/121）。 |
| **S4a** 重复区，SNV 以 edit 形式出现在别处 | 70 | 17.5% | 0.70 | reads 走了分支，或因附近的 indel 位置错开，SNV 成了另一个节点或另一个位置上的错配。候选存在，但 ID 和真值键不同；其中 20 个位点的这个tensor被标成 **germline**。 |
| **S4b** 重复区，被图的路径吸收 | 127 | 31.8% | 0.50 | 重复单元的组成或长度改变（例如交界处 TA→GA），reads 走现成的跳过边或分支，没有 SNV edit。 |
| **S4c** 重复区，没有清楚的 ALT reads | 44 | 11.0% | 0 | 局部单倍型和 GRCh38 不同（indel、STR 长度不同），但没有 read 更接近"GRCh38 + SNV"。 |
| **S2** 序列背景简单，reads 里没有 ALT | 46 | 11.5% | 0 | 大多数 reads 与 REF 等长，SNV 位置上是 REF 碱基，ALT 读数为 0 或 1–2 条。例如一个 VNTR 扩增的等位基因被投影成了一串 SNV（chr10:132826472 簇）；或者这个 ALT 在这份 PacBio 数据里本来就不存在。 |
| **S1** 所有 reads MAPQ≤10，或没有 reads | 48 | 12.0% | — | 片段重复区（chr1:143–145 Mb、chr15:20–30 Mb、chr7:62/74–75 Mb 等），被 `min_mapq 10` 全部过滤掉。 |

- S4a + S4b + S4c 共 241 个（60.3%），都是局部单倍型和 GRCh38 长度不同的重复区。v5 的统计：按"大多数 reads 附近有 indel"筛出的 239 个里，85%（204）在 20 bp 内有正常样本的 germline indel。
- **38 个 SNV 的 germline 真值在同一位置有同一个 ALT**（somatic benchmark 与 dipcall 冲突），例如 chr11:43270499 C>T（germline `C>T(2/1)`）。其中 33 个落在 S3/S4，ALT 比例中位数 0.86，可能是小范围的杂合性丢失（基因转换）。这些位点不适合当作 somatic 正例。

v5 的 404 个里，下面 4 个在 v6 有了候选，但没通过过滤：chr15:37166386 G>C、chr7:49588809 T>A（S4a，min_af），chr10:18220238 G>T、chr11:93286605 T>A（S4b，min_variants）。它们的节点不在 v6 的主 parts 里，只在 supplement 轮里被建（已核对）。

## 8. 方法

读取层面的分析和 v5 一样（见 `v5_reference/v5_REPORT.md` 第 2 节），改动如下：
1. 解码用 **v6 冻结源码**（`v6_run/source`）的 `decode_alignment`（左对齐 indel）。v6 会跨节点左对齐，所以每条 record 只保留位点附近的 mappings，再加两侧各 ≥600 个 read 碱基（`reads.py` 的 `trim`）。在 60 个位点上和整条 record 解码比较过，结果完全相同，0 差异（`validate_trim.py`）。
2. target 节点 = v6 的常染色体 parts + supplement 01–03。
3. 残余 edits 的标签取自 v6 tensors（每个 tensor 的全部 site alleles）。
4. 类别的对应关系和 v5 相同，另外加了 I8–I10（复杂替换 / truth edit 出现在未建或已建的节点上），但未命中里都没有出现。
5. 对照：200 个 INDEL 里 181 个是"truth 键出现在 target 节点上"（INS 67/75，DEL 114/125），其余是多数 reads 绕开 GRCh38、但仍有足够的 reads 走 GRCh38。

### 对照（有 tensor 的 INDEL truth，同样方法分析）

| reason | INS | DEL |
|---|---:|---:|
| key_observed_on_target_node | 67 | 114 |
| site_bypassed:skip_edge:edit_on_path_taken | 1 | 5 |
| site_bypassed:branch:no_edit | 3 | 2 |
| site_bypassed:skip_edge:no_edit | 0 | 4 |
| site_on_grch38:other_edit | 2 | 0 |
| site_bypassed:branch:edit_on_path_taken | 2 | 0 |

## 9. 建议

1. **标签（可以马上做）**：tensor 如果匹配 germline 真值，但落在某个 somatic 真值的等价范围内，改标 −1。这样可以去掉 301 个 INDEL 残余 edit 的 germline 错标（SNV 另有约 30 个）。
2. **标签（需要设计）**：按单倍型比对真值。把 ALT reads 走的分支和残余 edit 合起来投影到 GRCh38，再和"GRCh38 + 真值"比较。这样 I2/I3 的残余 edit tensor，包括分支节点上没有 GRCh38 坐标的那些，都能对上真值。
3. **左对齐（需要讨论）**：v6 的左对齐让 36 个 truth 丢了 tensor、305 个从 filtered 变成没有候选，都是 edit 被推进了人群插入分支节点。单倍型投影（第 2 条）能同时解决。也可以考虑不让 indel 左对齐越过"GRCh38 → 非 GRCh38"的节点边界，但这会改变 tensor，需要评估。
4. **候选生成（需要讨论）**：I1 这类 allele 完全在图里的变异，靠 edits 永远看不到，需要基于图路径的候选。对 tumor-only 来说，germline 的分支使用会带来大量干扰。

## 10. 文件

都在 `tmp/somatic_miss_analysis_v6/`：

| 文件 | 内容 |
|---|---|
| `REPORT.md` | 本报告 |
| `indel_no_candidate.tsv` | 2,994 个 INDEL 未命中，每个一行：`final_class`、reads 数、REF/ALT exact/like、ALT 比例、附近 somatic/germline 数、残余 edits 及其 v6 tensor 标签、ALT reads 的表示、序列背景 |
| `snv_no_candidate.tsv` | 400 个 SNV 未命中（v5 分析的行，去掉 v6 里变成 filtered 的 4 个） |
| `controls.tsv` | 200 个 INDEL 对照 |
| `transitions.tsv` | 每个 INDEL truth 的 v5 状态/类别 → v6 状态/类别 |
| `residual_edits_labelled_germline.tsv` | 残余 edit tensor 被标成 germline 的 301 行 |
| `show_examples.txt`、`examples_pick.tsv`、`show/` | 例子的 reads 逐条显示 |
| `tables.md` | 本报告的全部表格（`tables.py` 生成） |
| `reads.py`、`report.py`、`tables.py`、`show_locus.py`、`validate_trim.py`、`pick_examples.py` | 分析脚本 |
| `reads_part*.jsonl` | 每个位点的原始结果 |
| `snv_from_v5/` | v5 的 SNV 分析文件（`snv_detail.*`、`snv_nodes.*`、`snv_compare.py` 等） |
| `v5_reference/` | v5 的召回表、INDEL 未命中表、报告、标签报告，以及 v5/v6 对比表（v5 数据本身已删除） |

列说明：`ref_exact`/`alt_exact` 是局部序列完全等于 REF/ALT 单倍型的 read 数；`ref_like`/`alt_like` 是按编辑距离更接近哪一边的 read 数；`alt_fraction` = (alt_exact + alt_like) / spanning；`residual_edits(reads)=tensor_label` 是 ≥3 条 ALT reads 在这里产生的候选，以及它们的 v6 tensor 标签。
