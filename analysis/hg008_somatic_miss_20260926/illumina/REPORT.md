# HG008 Illumina WGS：没有候选的 somatic 真值，原因分析（和 ONT、PacBio v6 对比）

- 数据：HG008-T p23 BCM Illumina WGS。tensors 来自 `indexed_gam_pipeline_v3`：
  - 代码 `6ef7c69`，续跑前打了 `12c29d3` 的 linkage 补丁；
  - 常染色体，format v6；
  - 深度超过 1 万的 683 个 node 下采样到 10,000 条 read；
  - 没有 supplement 轮。
  - 目录是 `pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors`。
- 标签：用 `truth-labels-v2` 原地重打过，和 PacBio v6、ONT 现在的规则相同；旧标签清单备份在 `tmp/relabel_20260925/Illumina/`。召回状态和重打前完全一样。
- 真值：
  - somatic：`HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz` + `_all.bed`；
  - germline：dipcall `dip.vcf.gz` + `dip.bed`。
- 对象：Illumina 召回报告里「没有候选」的 somatic 等位基因，共 **4,293** 个：SNV 339，INS 2,633，DEL 1,321。每个都按 reads 分析了。
- 对照：随机抽 250 个有 Illumina tensor 的真值（INDEL ≥2 bp 150 个、1 bp 50 个，SNV 50 个）。
- 方法：和 PacBio v6、ONT 两份分析相同，脚本从 `tmp/somatic_miss_analysis_ont/` 复制过来，只改了路径。
  - 解码用这次运行冻结源码里的 Python 解码器。
  - 每条 record 只保留位点两侧各 600 个 read 碱基。Illumina 的 read 只有 150 bp，这一步等于保留整条 record，所以没有另做截断验证。
  - SNV 的 S1–S4c 分类规则和 ONT 分析相同。

## 1. 结论

1. **主要原因和长读段一样：somatic allele 已经全部或部分在图里。**
   - INDEL 未命中里 I1 + I2 占 INS 的 **83.2%**（2,192 / 2,633），占 DEL 的 **93.0%**（1,228 / 1,321）。
   - DEL 里以 I1 为主（59%），最常见的原因是 ALT reads 走图里现成的跳过边、没有任何 edits（676 个）。
2. **Illumina 的 INDEL 召回低于长读段，多出来的漏检主要也是 I1。**
   - 真值有 tensor（代表或非代表 allele）的比例：DEL 51.7%（ONT 65.2%，PacBio 65.4%）；INS 36.3%（ONT 44.0%，PacBio 42.9%）。
   - ONT 有代表 tensor、Illumina 没有候选的：DEL 321 个，其中 279 个是 I1；INS 216 个，其中 141 个是 I1。
   - 例子 chr4:32683412 GA>G：两个平台的 ALT reads 都走跳过边、没有 edits。ONT 有 tensor，是因为部分 ONT read 又多丢了一个 A（同聚物噪声），左对齐后正好落在真值位置。Illumina 不产生这种噪声 edit。
   - 也就是说，**长读段对一部分 I1 类真值的召回，来自噪声 edits 碰巧落在真值位置**。这个判断目前来自例子，还没有逐个统计。
3. **Illumina 特有的漏检：长插入、长重复和片段重复区。**
   - **I4（ALT reads 太少）：** INS 264 个（ONT 40，PacBio 51）。长度多在 6–20 bp（125）、21–50 bp（57）、>50 bp（23），中位只有 15 条 read 跨过位点。长的 STR 插入在短读段里几乎看不到 ALT，例子 chr11:14293124（+16 bp 的 AT 重复）50 条 read 里只有 1 条 ALT。
   - **I5（没有 read 跨过位点或低 MAPQ）：** INS 84 个（ONT 1），中位 0 条 read 跨过，低 MAPQ 比例只有 1.6%，说明主要是 150 bp 跨不过长重复，例子 chr1:34802227 是在 220 bp 的 GGTTT 重复里插入 60 bp。DEL 17 个，低 MAPQ 比例中位 0.99，是片段重复区，例子 chr1:144665794。
   - **I6（>50 bp）：** Illumina 只有 1 个（ONT 15）。短读段写不出 >50 bp 的 edit，这类真值落进了 I4 / I5。
   - **SNV 的 S1（低 MAPQ）：** 55 个（ONT 1，PacBio 48），在片段重复区没有 MAPQ>10 的短读段（例如 chr1p36 PRAMEF 区），分布在 chr7（12）、chr1（9）、chr2（7）等。ONT 有 34 个有 tensor（33 个代表），PacBio 17 个。
4. **Illumina 的噪声 allele 少，filtered 也少。**
   - 非代表 allele：DEL 49、INS 29（ONT 446、416）。
   - filtered：DEL 390、INS 519（ONT 510、1,042）。
   - 短读段很少把同一个 allele 写成多种样子，所以卡在「同一写法不到 3 条」的情况少。
5. **I1 / I2 位点上 Illumina 几乎看不到 REF read。**
   - I1 / I2 的「ALT 比例」中位数 0.89 / 0.95（ONT 0.65 / 0.66）；对照组只有 0.59。
   - 在两个平台都漏掉的 1,752 个 I1 / I2 位点上，像 REF 的 read 中位数 Illumina 0 条、ONT 6 条。
   - 逐条看的 4 个位点里，ONT 的「像 REF」多是同聚物长度噪声；Illumina 的 read 长度都等于 ALT，被归为「other」的是长 T 串之后零散的替换错误（Illumina 在长同聚物之后的典型错误）。
   - 这里的「ALT 比例」是「比 REF 更接近 ALT」的 read 比例，不等于 VAF。
6. **残余 edits 的 tensor 被标成 germline 的更多：317 个不同候选**（ONT 194），集中在 INS I2（122 个位点）和 DEL I2（102 个位点）。
7. **下采样没有影响真值：** 683 个深 node 上没有任何 somatic 真值的键。
8. SNV：30 个 somatic SNV 的 germline 真值在同一位置、同一个 ALT（ONT 35），不适合当 somatic 正例。

## 2. 召回总览

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

**三个平台对比（全部 somatic 等位基因）**

| | Illumina | ONT | PacBio v6 |
|---|---:|---:|---:|
| SNV 代表 allele | 95.3% | 95.2% | 94.8% |
| DEL 有 tensor（代表 + 非代表） | 51.7% | 65.2% | 65.4% |
| DEL 代表 allele | 50.4% | 52.6% | 59.3% |
| INS 有 tensor（代表 + 非代表） | 36.3% | 44.0% | 42.9% |
| INS 代表 allele | 35.8% | 35.6% | 39.6% |
| 没有候选 SNV / DEL / INS | 339 / 1,321 / 2,633 | 342 / 725 / 1,729 | 400 / 890 / 2,104 |

## 3. ONT / PacBio v6 → Illumina 的变化

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

**SNP：ONT 状态（行）→ Illumina 状态（列）**

| ONT \ Illumina | rep | non-rep | filtered | no_cand | 合计 |
|---|---:|---:|---:|---:|---:|
| rep | 8197 | 0 | 14 | 64 | 8275 |
| non-rep | 2 | 2 | 0 | 1 | 5 |
| filtered | 35 | 0 | 8 | 25 | 68 |
| no_cand | 45 | 2 | 46 | 249 | 342 |

**ONT 有 tensor 或 filtered、Illumina 没有候选时，Illumina 的类别**

| ONT 状态 | 种类 | 个数 | I1 | I2 | I3 | I4 | I5 |
|---|---|---:|---:|---:|---:|---:|---:|
| rep | DEL | 321 | 279 | 33 | 2 | 5 | 2 |
| rep | INS | 216 | 141 | 58 | 3 | 8 | 6 |
| filtered | DEL | 292 | 181 | 95 | 6 | 5 | 5 |
| filtered | INS | 716 | 364 | 290 | 7 | 41 | 14 |

**ONT / PacBio v6 没有候选的真值，在 Illumina 的状态（按对方的类别）**

| 类别 | 种类 | ONT → Illumina rep / non-rep / filtered / no_cand | PacBio v6 → Illumina rep / non-rep / filtered / no_cand |
|---|---|---|---|
| S1 | SNP | 0 / 0 / 0 / 1 | 22 / 0 / 0 / 26 |
| S3 | SNP | 8 / 0 / 6 / 32 | 12 / 1 / 5 / 47 |
| S4b | SNP | 20 / 2 / 25 / 107 | 18 / 0 / 15 / 94 |
| I1 | DEL | 7 / 0 / 14 / 285 | 5 / 0 / 17 / 393 |
| I2 | DEL | 40 / 10 / 53 / 266 | 48 / 10 / 75 / 286 |
| I1 | INS | 4 / 0 / 19 / 660 | 4 / 0 / 37 / 921 |
| I2 | INS | 9 / 1 / 71 / 798 | 14 / 2 / 85 / 860 |
| I5 | DEL | — | 2 / 0 / 2 / 8 |
| I5 | INS | 0 / 0 / 0 / 1 | 3 / 0 / 1 / 13 |

长读段的 I1 / I2 在 Illumina 里几乎都还是没有候选：图里已有的 allele，三种平台都「看不见」。完整的状态转移表在 `tables.md`。

## 4. INDEL 原因分类（Illumina）

### Illumina 没有候选的 INDEL：原因分类

| 类别 | INS Illumina | INS ONT | INS PacBio | DEL Illumina | DEL ONT | DEL PacBio | ALT 比例中位数 (Illumina) | 按无误 ALT reads 判定 (Illumina) | PASS∩BED (Illumina) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| I1_allele_fully_in_graph | 1136 (43.1%) | 683 | 962 | 780 (59.0%) | 306 | 415 | 0.889 | 1610 | 1790 |
| I2_partly_in_graph_residual_edit | 1056 (40.1%) | 879 | 961 | 448 (33.9%) | 369 | 419 | 0.953 | 803 | 1397 |
| I3_grch38_path_different_edit | 81 (3.1%) | 95 | 83 | 25 (1.9%) | 18 | 17 | 0.628 | 43 | 102 |
| I4_no_or_few_alt_reads | 264 (10.0%) | 40 | 51 | 42 (3.2%) | 20 | 16 | 0.017 | 0 | 277 |
| I5_low_mapq_or_no_reads | 84 (3.2%) | 1 | 17 | 17 (1.3%) | 0 | 12 | - | 0 | 91 |
| I6_longer_than_50bp | 0 (0.0%) | 11 | 11 | 1 (0.1%) | 4 | 4 | 0.667 | 1 | 0 |
| I7_unclear | 12 (0.5%) | 20 | 19 | 8 (0.6%) | 7 | 7 | 1.0 | 8 | 20 |
| I9_truth_edit_on_unbuilt_node | 0 | 0 | 0 | 0 | 1 | 0 | | | |
| 合计 | 2633 | 1729 | 2104 | 1321 | 725 | 890 | | | |

「ALT 比例」= 比 REF 更接近 ALT 的 read 占跨过位点的 read 的比例（第 1 节第 5 点）。

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

**I1 / I2 / I4 / I5 的长度和 reads**

| 类别 | 种类 | 个数 | 1 bp | 2-5 | 6-20 | 21-50 | >50 | 跨过位点的 reads（中位数） | 低 MAPQ 比例（中位数） |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| I1 | INS | 1136 | 620 | 343 | 150 | 21 | 2 | 42 | 0.007 |
| I1 | DEL | 780 | 461 | 215 | 82 | 21 | 1 | 48 | 0.006 |
| I2 | INS | 1056 | 161 | 487 | 343 | 63 | 2 | 30 | 0.007 |
| I2 | DEL | 448 | 43 | 219 | 151 | 28 | 7 | 38 | 0.008 |
| I4 | INS | 264 | 22 | 37 | 125 | 57 | 23 | 15 | 0.011 |
| I4 | DEL | 42 | 9 | 9 | 12 | 10 | 2 | 6 | 0.020 |
| I5 | INS | 84 | 11 | 14 | 23 | 19 | 17 | 0 | 0.016 |
| I5 | DEL | 17 | 8 | 2 | 4 | 0 | 3 | 0 | 0.99 |

**长度（Illumina，全部 INDEL 真值）**

| | 1 | 2-5 | 6-20 | 21-50 | >50 |
|---|---:|---:|---:|---:|---:|
| INS 有 tensor | 1604 | 169 | 26 | 0 | 0 |
| INS filtered | 204 | 264 | 44 | 7 | 0 |
| INS 没有候选 | 829 | 923 | 670 | 167 | 44 |
| DEL 有 tensor | 1354 | 306 | 147 | 27 | 0 |
| DEL filtered | 153 | 200 | 35 | 2 | 0 |
| DEL 没有候选 | 529 | 459 | 259 | 59 | 15 |

21 bp 以上的 INS，Illumina 一个 tensor 都没有（ONT 7 个）。

## 5. 残余 edits 的 tensor 标签

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

逐条列表：`residual_edits_labelled_germline.tsv`，334 行，317 个不同候选（INS 189 行、DEL 121 行、SNV 24 行）。

## 6. SNV（339）

### SNV 没有候选（Illumina 339；ONT 342；PacBio v6 400）

| 类别 | Illumina | Illumina 占比 | Illumina ALT 比例中位数 | ONT | PacBio v6 |
|---|---:|---:|---:|---:|---:|
| S1_low_mapq_or_no_reads | 55 | 16.2% | - | 1 | 48 |
| S2_alt_absent_in_reads | 40 | 11.8% | 0.0 | 36 | 46 |
| S3_alt_is_existing_graph_allele | 76 | 22.4% | 0.913 | 46 | 65 |
| S4a_repeat_snv_edit_on_branch | 57 | 16.8% | 0.818 | 61 | 70 |
| S4b_repeat_absorbed_by_graph_path | 65 | 19.2% | 0.621 | 154 | 127 |
| S4c_repeat_no_clear_alt | 46 | 13.6% | 0.0 | 44 | 44 |

PacBio 的 48 个 S1，Illumina 有 22 个成了代表 allele（ONT 40 个）。反过来，Illumina 的 55 个 S1 里，ONT 有 34 个有 tensor（33 个代表），PacBio 17 个；另有 14 个在 ONT 里也没有候选。

**和命中 SNV 的对比**（`snv_compare.py`；命中 8,279 个，未命中 339 个）

| 特征 | 命中 | 没有候选 |
|---|---:|---:|
| SNV 所在的 GRCh38 节点只有 1 bp | 2.9% | **64.3%** |
| 节点 2–32 bp | 4.4% | 18.3% |
| 节点 >32 bp | 92.7% | 17.4% |
| 50 bp 内有 ≥2 个 germline 变异 | 4.1% | **67.3%** |
| 50 bp 内没有 germline 变异 | 83.4% | 18.0% |
| germline 真值在同一位置、同一个 ALT | 0.1% | **8.8%（30 个）** |

## 7. 例子

每行都用 `show_locus.py` 把 reads 逐条打出来看过，完整输出在 `show_examples.txt`（36 个按类别挑的位点）。ONT 的对比显示在 `ont_show_compare/`。

| 位点 | 真值 | 类别 | reads 实际的样子 |
|---|---|---|---|
| chr13:54018939 | GT>G | I1 | 119 条 ALT reads 走跳过边 14386371→14386373，全程没有 edits；只有 1 条 REF |
| chr1:32405697 | A>AT | I1 | 所有 read 都是 17 个 T，走多一个 T 的非 GRCh38 节点 362336，没有 edits；「other」是 T 串之后零散的替换错误。ONT 同一位点 REF、ALT 都有，还有 18–23 个 T 的噪声 read |
| chr4:32683412 | GA>G | I1（ONT 有 tensor） | 两个平台的 ALT reads 都走跳过边 37858837→37858839、没有 edits；ONT 的 tensor 来自又多丢一个 A 的噪声 read |
| chr20:58872348 | T>TAA | I2 | reads 走「+A」的非 GRCh38 节点 31293240；92 条再多写一个 1 bp 插入（残余 edit），合起来等于 ALT；其余 81 条多数只走分支（+A），少数带替换错误 |
| chr1:198137519 | +TATAGAGAGAG | I3 | ALT reads 在 GRCh38 节点上写出 8 bp 插入，位置和写法与真值的键不同；REF 52、ALT 35 条 |
| chr11:14293124 | +16 bp AT 重复 | I4 | 50 条 read 里只有 1 条带 ALT，大多数是短 2 bp 的另一个 allele（DD） |
| chr1:34802227 | +60 bp GGTTT 重复 | I5 | 重复有 220 bp，150 bp 的 read 跨不过去，没有任何 read 跨过位点 |
| chr1:144665794 | AT>A | I5 | chr1q21 片段重复区，没有 MAPQ>10 的 read |
| chr1:13086788 | C>T | S1 | chr1p36 片段重复区，没有 MAPQ>10 的 read |
| chr13:76153195 | G>A | S3 | ALT reads 走 SNP 分叉上的 A 节点，没有错配（ONT 例子里同一个位点） |

## 8. 对照

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

对照组里 SNV 50/50、INS 54/65、DEL 109/135 是「真值的键出现在 target 节点上」，和 ONT 对照相近。对照的 ALT 比例中位数 INDEL 0.585、SNV 0.544，正常。

## 9. 建议

1. **标签**：和 PacBio、ONT 的建议相同。tensor 落在 somatic 真值的等价范围内时，不标 germline 而标 −1（Illumina 317 个候选）；用单倍型投影匹配真值。
2. **I1 是三个平台共同的主要问题，而且长读段的部分召回来自噪声。** 只靠 edits 找候选，看不见图里已有的 somatic allele。建议优先做「基于路径的候选」：read 走了非 GRCh38 的分支或跳过边、且这条路径对应的 allele 在该位点有足够支持时，也生成候选。这对 Illumina 的收益最大（I1 INS 1,136、DEL 780）。
3. **Illumina 的长插入和长重复（I4 / I5，约 400 个）是短读段本身的限制**，靠调参解决不了，适合交给长读段数据。
4. **INDEL 的 AF 门槛**：之前的测试显示，降到 0.05 只多召回 92 个真值 indel，而没有候选的有一半以上，瓶颈不在 AF。
5. **supplement 和下采样**：Illumina 没有丢任何 I9；683 个深 node 上没有 somatic 真值。两者都可以维持现状。

## 10. 文件

都在 `tmp/somatic_miss_analysis_illumina/`：

| 文件 | 内容 |
|---|---|
| `REPORT.md` | 本报告 |
| `indel_no_candidate.tsv` | 3,954 个 INDEL 未命中，每个一行（列同 ONT / PacBio v6） |
| `snv_no_candidate.tsv` | 339 个 SNV 未命中，另有 SNV 专用列 |
| `controls.tsv` | 250 个对照 |
| `transitions.tsv` | 每个真值的 PacBio v6、ONT、Illumina 状态和类别 |
| `residual_edits_labelled_germline.tsv` | 残余 edit tensor 被标成 germline 的 334 行 |
| `show_examples.txt`、`examples_pick.tsv`、`show/` | 例子的 reads 逐条显示 |
| `ont_show_compare/`、`refgap_loci.txt` | 同一批位点在 ONT 里的 reads（第 1 节第 2、5 点） |
| `tables.md` | 全部表格（`tables.py` 生成） |
| `reads.py`、`snv_detail.py`、`snv_nodes.py`、`snv_compare.py`、`report.py`、`tables.py`、`show_locus.py`、`pick_examples.py` | 分析脚本 |
| `reads_part*.jsonl`、`snv_detail.jsonl`、`snv_nodes.json`、`summary.json` | 原始结果 |
| `somatic.graph.tsv`、`somatic.recall.v1.tsv` | 分析时用的真值表和召回表快照 |
