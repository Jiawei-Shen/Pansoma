# HG008 PacBio v5：没有候选的 somatic 真值，原因分析

- 数据：HG008-T PacBio Revio 116×，v5 张量（job 363651），合并后按 chr 打的标签（job 365671，已含下文的插入范围修复）。
- 真值：somatic `HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz` + `_all.bed`；germline dipcall `dip.vcf.gz` + `dip.bed`。
- 对象：召回报告里状态为"没有候选"的 somatic 等位基因，共 **3,188** 个（SNV 404，INS 1,826，DEL 958）。"没有候选"指：它的任何一种节点坐标写法，既不是任何张量的等位基因，也没有出现在 `filtered_candidates` 里。
- 对照：随机抽 250 个"有张量"的 somatic 真值（SNV 50，INS 78，DEL 122），用同样的方法分析。

## 1. 结论

1. **主要原因：这些 somatic 等位基因在泛基因组图里已经全部或部分存在。** pipeline 只在 reads 相对它所走的图节点有编辑（错配、插入、缺失）的地方产生候选。HPRC 图里已经有这个等位基因时，带 ALT 的 reads 会走那条现成的分支，或走一条跳过参考碱基的边，于是要么完全没有编辑，要么只剩一个写法不同的小编辑。
   - INDEL：**90%**（INS 1,644/1,826，DEL 864/958）属于这种情况，而且大多数（INS 1,241、DEL 648）是按完全无测序错误的 ALT reads 判定的。
   - SNV：约 **66%**（S3 + S4a + S4b，266/404）。
2. **SNV 的未命中集中在重复区和多态区域**，和命中的 SNV 差别非常明显（第 3.1 节）：61.6% 落在 1 bp 的 GRCh38 节点上（命中的只有 2.7%），70.5% 在 50 bp 内有 ≥2 个 germline 变异（命中的 3.7%）。
3. **标签的影响**：在"部分被图吸收"的位点上，残余编辑本身有张量，但 reads 实际带的是 somatic 单倍型。其中 **370 个不同的候选被标成了 germline**（DEL 169，INS 170，SNP 31），这是标签错误；另有约 4,500 个标成了 −1（不参与训练），见 `residual_edits_labelled_germline.tsv`。
4. **次要原因**：
   - 所有 reads 都 MAPQ≤10 或没有 reads（SNV 48，INDEL 29）；
   - PacBio reads 里根本没有这个 ALT（SNV 46 个在简单序列背景下也没有 ALT reads，INDEL 68 个）；
   - indel 超过 50 bp（按 reads 判定 12 个；按真值长度 59 个）；
   - 跨节点缺失在 v5 里被拆成两个候选（v6 已修复）。
5. **分析中发现并已修复的标签 bug**：判断"附近 10 bp 内有没有真值"时，插入只用了最左边的等价位置。修复后重新打标签，**37,065 个 INDEL 张量和 1,119 个 SNV 张量**从 non 改为 −1，somatic 和 germline 的数量不变。

## 2. 方法

### 2.1 读取层面的分析（`reads.py`）

对每个位点：
1. 取出这个位点附近 ±300 bp 的 GRCh38 节点上的全部 GAM 记录，保留 MAPQ>10 的（和 build 一致）。
2. 用 **v5 冻结源码**的 `decode_alignment` 解码，只解码位点附近的 mapping。v5 按 mapping 逐个解码，这样得到的候选和整条解码完全一样；在 15 个位点上验证过，0 差异。
3. 取变异等价范围外侧、距离 ≥30 bp 的最近 GRCh38 匹配碱基作为锚点（最远 2 kb），把两个锚点之间 read 的碱基拼成局部序列，和 REF 单倍型以及"GRCh38 + 真值"的 ALT 单倍型比较：
   - 完全相等 → exact；
   - 否则比编辑距离，离 ALT 更近 → ALT-like，离 REF 更近 → REF-like，一样远 → unclear。
4. 对每条 ALT reads 判断它在变异区间里走的是不是 GRCh38 的碱基，以及线性位置是否连续：
   - 有非 GRCh38 节点 → **走了分支**；
   - 只有 GRCh38 节点但位置有跳跃 → **走了跳过参考碱基的边**；
   - 都不是 → **走 GRCh38**。
   同时记录它在这里产生的候选、是否等于真值键，以及这些候选的张量被打了什么标签。
5. 位点的主要原因 = 完全无误的 ALT reads 的多数表示；没有这种 read 时，才用带测序错误的 ALT reads。

### 2.2 SNV 逐个细查（`snv_detail.py`）

在重复区里，"局部序列等于 GRCh38 + SNV"几乎不可能成立，所以对每个 SNV 额外统计：
- 覆盖位点的 reads 里，有多少条局部长度和 REF 相同（两侧 ±20 bp 内没有 indel）；
- 这些等长 reads 在 SNV 位置上是 REF、ALT 还是其他碱基；
- ALT reads 在这个位置走的是 GRCh38 节点还是别的节点。

### 2.3 核对

- 对照组全部按预期分类：SNV 50/50、INS 69/78、DEL 93/122 是"ALT reads 在目标节点上产生了真值键"。其余对照是多数 reads 绕开了 GRCh38，但仍有足够的 reads 走 GRCh38，所以也有张量。
- 每个类别都挑了例子，把 reads 的局部序列和图路径逐条打印出来看过（`show_locus.py`，见第 5 节）。

## 3. SNV（404）

### 3.1 和命中的 SNV 对比（`snv_compare.py`）

| 特征 | 命中（8,234） | 没有候选（404） |
|---|---:|---:|
| SNV 所在的 GRCh38 节点只有 1 bp（图里在这个位置本来就有 SNP 分叉） | 2.7% | **61.6%** |
| 节点 2–32 bp | 4.2% | 21.0% |
| 节点 >32 bp | 93.1% | 17.3% |
| 50 bp 内有 ≥2 个 germline 变异 | 3.7% | **70.5%** |
| 50 bp 内没有 germline 变异 | 83.8% | 15.1% |
| 50 bp 内还有别的 somatic 真值 | 7% | 58% |
| germline 真值在同一位置、同一个 ALT | 0.1% | **9.4%（38 个）** |

### 3.2 原因分类

| 类别 | 个数 | 占比 | ALT 比例中位数 | 含义 |
|---|---:|---:|---:|---|
| **S3** ALT 是图里现成的等位基因 | 65 | 16.1% | 0.84 | 序列背景简单（reads 与 REF 等长），带 ALT 的 reads 在这个位置走非 GRCh38 的节点（1 bp SNP 分叉的另一个碱基），没有错配。ALT 比例很高，很多接近 1（例如 chr13:87814362，121/121）。 |
| **S4a** 重复区，SNV 以编辑形式出现在别处 | 72 | 17.8% | 0.70 | reads 走了分支，或因附近的 indel 位置错开，SNV 成了另一个节点或另一个位置上的错配。候选存在，但 ID 和真值键不同；其中 20 个位点的这个张量被标成 **germline**。 |
| **S4b** 重复区，被图的路径吸收 | 129 | 31.9% | 0.50 | 重复单元的组成或长度改变（例如交界处 TA→GA），reads 走现成的跳过边或分支，没有 SNV 编辑。 |
| **S4c** 重复区，没有清楚的 ALT reads | 44 | 10.9% | 0 | 局部单倍型和 GRCh38 不同（indel、STR 长度不同），但没有 read 更接近"GRCh38 + SNV"。 |
| **S2** 序列背景简单，reads 里没有 ALT | 46 | 11.4% | 0 | 大多数 reads 与 REF 等长，SNV 位置上是 REF 碱基，ALT 读数为 0 或 1–2 条。例如一个 VNTR 扩增的等位基因被投影成了一串 SNV（chr10:132826472 簇）；或者这个 ALT 在这份 PacBio 数据里本来就不存在。 |
| **S1** 所有 reads MAPQ≤10，或没有 reads | 48 | 11.9% | — | 片段重复区（chr1:143–145 Mb、chr15:20–30 Mb、chr7:62/74–75 Mb 等），被 `min_mapq 10` 全部过滤掉。 |

- S4a + S4b + S4c 共 245 个（60.6%），都是局部单倍型和 GRCh38 长度不同的重复区。按"大多数 reads 附近有 indel"筛出的 239 个里，85%（204）在 20 bp 内有正常样本的 germline indel。
- **38 个 SNV 的 germline 真值在同一位置有同一个 ALT**（somatic benchmark 与 dipcall 冲突），例如 chr11:43270499 C>T（germline `C>T(2/1)`）。其中 33 个落在 S3/S4，ALT 比例中位数 0.86，可能是小范围的杂合性丢失（基因转换）。这些位点不适合当作 somatic 正例。

## 4. INDEL（INS 1,826，DEL 958）

### 4.1 长度：没有候选的集中在 2–50 bp

| | 1 bp | 2–5 bp | 6–20 bp | 21–50 bp | >50 bp |
|---|---:|---:|---:|---:|---:|
| INS 有张量 | 1,711 | 162 | 45 | 5 | 0 |
| INS 没有候选 | 315 | 690 | 623 | 154 | 44 |
| DEL 有张量 | 1,663 | 254 | 109 | 17 | 0 |
| DEL 没有候选 | 107 | 452 | 315 | 69 | 15 |

### 4.2 原因分类

| 类别 | INS | DEL | 含义 |
|---|---:|---:|---|
| **I1** 等位基因完全在图里（无编辑） | 827（45.3%） | 435（45.4%） | ALT reads 走现成分支（多出来的 T、TT、A… 节点），或走跳过参考碱基的边，全程无编辑。INS 613 个、DEL 325 个按完全无误的 ALT reads 判定。 |
| **I2** 部分在图里 + 残余编辑 | 817（44.7%） | 429（44.8%） | ALT reads 走一条接近的分支或跳过边，剩下的差异成为一个小编辑（例如 somatic +TT 变成"分支 +3T、再 DEL T"）。残余编辑的张量被标成 germline 的位点：INS 77、DEL 82。 |
| **I3** 走 GRCh38，但编辑写法不同 | 80（4.4%） | 55（5.7%） | 例如 13 bp 的 polyT 缺失被 v5 拆成两个相邻节点上的 `DEL T` + `DEL 12×T`（v6 已合并）；或写成错配加 indel。 |
| **I4** 没有或只有 1–2 条 ALT reads | 51（2.8%） | 17（1.8%） | 多为 STR：germline 的 STR 长度已经和 GRCh38 不同，没有 read 等于"GRCh38 + ALT"。 |
| **I5** MAPQ≤10 或没有 reads | 17（0.9%） | 12（1.3%） | |
| **I6** reads 里的 indel 超过 50 bp | 10（0.5%） | 2（0.2%） | 例如 GGTTT VNTR +60 bp、CTTT 重复 +145 bp。 |
| **I7** 不确定 | 24（1.3%） | 8（0.8%） | ALT-like reads 走 GRCh38，但附近没有编辑（多在长重复区，编辑距离判断不稳定）。 |

## 5. 例子

每个类别都有例子，reads 的局部序列和图路径都逐条看过，见 **`examples.tsv`**（34 个）。几个典型的：

| 位点 | 真值 | 类别 | reads 实际的样子 |
|---|---|---|---|
| chr5:24851177 | A>AT | I1 | 66 条无误 ALT reads 走多一个 T 的非 GRCh38 节点 41817187，没有编辑 |
| chr3:151219891 | AC>A | I1 | 42 条 ALT reads 走跳过一个 C 的边（36110549→36110552），没有编辑 |
| chr1:210849365 | C>CTT | I2 | 35 条 ALT reads 走多 3 个 T 的分支，再 `DEL T`：somatic 插入变成了一个缺失候选（标签 −1） |
| chr4:108843148 | CAAAA>C | I2 | ALT reads 跳过 3 个 A，再 `DEL A`；这个 `−A` 张量被标成 **germline** |
| chr2:11519449 | −13 个 T | I3 | 拆成相邻节点上的 `DEL T` + `DEL 12×T` 两个候选（v6 已修复） |
| chr13:87814362 | A>T | S3 | 121/121 条 reads 走 SNP 分叉的另一个节点，没有错配；germline 同位点 A>T(2/1) |
| chr10:10988294 | T>G | S4b | TA\|GA 交界一个单元改变，47 条无误 ALT reads 走跳过边 + 非 GRCh38 的 GA 节点，全程无编辑 |
| chr3:163930271 | T>G | S4a | reads 把它写成 GRCh38 节点 36344263 上的 T>G（57 reads，位置因 2 bp 缺失错开），张量被标成 germline |
| chr10:132826472 | T>C（簇内 8 个） | S2 | 29 条 reads 完全等于 GRCh38，另 26 条带一段长得多的 GAGCC VNTR；这些 SNV 像是 VNTR 等位基因的投影 |
| chr1:103652711 | G>A | S1 | 84 条 reads 全部 MAPQ≤10（1q21 片段重复） |

## 6. 建议

1. **标签（可以马上做）**：一个张量如果匹配 germline 真值，但位于某个 somatic 真值的等价范围内（或 10 bp 内），改标 **−1**，而不是 germline。这样可以去掉上述 370 个错误的 germline 标签。
2. **标签（需要设计）**：按单倍型比对真值，而不是只比较编辑 ID。把张量里 ALT reads 走的图路径和编辑拼成局部序列，和"GRCh38 + 真值"的单倍型比较。I2、S4a 里残余编辑的张量就能被正确标成 somatic。
3. **候选生成（需要讨论）**：I1、S3、S4b 这类"等位基因已经在图里"的 somatic 变异，只靠编辑永远看不到。要检测它们，需要基于图路径的候选：在每个分叉处统计 reads 走各分支的比例，把 reads 走了少见分支（通道 6 的路径数很小）的位置也作为候选。这对 tumor-only 是个设计问题，germline 的分支使用会造成大量干扰。
4. **v6**：跨节点缺失已经合并成一个候选，真值键也已支持 `@` 写法，I3 里这部分会自然恢复。
5. **真值冲突**：38 个 somatic/germline 同位点同等位基因的 SNV，建议从 somatic 正例里排除或单独标记。

## 7. 文件

都在 `tmp/somatic_miss_analysis_20260923/`：

| 文件 | 内容 |
|---|---|
| `REPORT.md` | 本报告 |
| `snv_no_candidate.tsv` | 404 个 SNV 未命中，每个一行：`final_class`、等长 reads 数、SNV 位置的碱基分布、ALT reads 走的节点、ALT/REF reads 数、附近的 somatic/germline 数量、节点长度和 discovery 统计、残余编辑及其张量标签 |
| `indel_no_candidate.tsv` | 2,784 个 INDEL 未命中，列同上（没有 SNV 专用列，另有 `context` 序列） |
| `controls.tsv` | 250 个对照 |
| `examples.tsv` | 34 个逐条看过 reads 的例子及说明 |
| `residual_edits_labelled_germline.tsv` | 残余编辑张量被标成 germline 的 380 行（370 个不同候选） |
| `reads.py`、`snv_detail.py`、`show_locus.py`、`report.py`、`nearby.py`、`snv_compare.py` | 分析脚本（可复现） |
| `reads*.jsonl`、`snv_detail.jsonl`、`nearby.json` | 原始结果 |

列说明：`ref_exact`/`alt_exact` 是局部序列完全等于 REF/ALT 单倍型的 read 数；`ref_like`/`alt_like` 是按编辑距离更接近的 read 数；`unclear` 是到两边一样远的 read 数；`alt_fraction` = (alt_exact + alt_like) / spanning；`residual_edits(reads)=tensor_label` 是 ≥3 条 ALT reads 在这里产生的候选及其张量标签。
