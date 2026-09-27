# truth-labels-v5：tensor 标签算法详解

- 代码：`indexed_gam_pipeline_v3/tensor_postprocessing/truth_labels.py`（`VERSION = "truth-labels-v5"`，最终规则 commit `6fa641f`，golden `ec7b3d2`）。
- 坐标换算：`indexed_gam_pipeline_v3/tensor_postprocessing/reference_path.py`（`ReferencePath.linear`）。
- 命令：`python -m indexed_gam_pipeline_v3.tensor_postprocessing label --tensors … --somatic-vcf … --somatic-bed … --germline-vcf … --germline-bed … --fasta … --reference-path … --truth-dir …`
- 本文的例子全部取自 **HG008 PacBio `v6_tensors` 的 v5 标签**（2026-09-27 重打），计数来自三个平台 v5 的 `labels.manifest.json`。

## 0. 名词

| 名词 | 意思 |
|---|---|
| tensor / site | 一个 site 是同一节点、同一起点、同一类型（SNV 或 INDEL）的一组候选 allele，只画成一个 tensor |
| A1 / 代表 allele | site 里 ALT reads 最多的 allele。tensor 的 `candidate_id` 就是它，画在最上面那组行 |
| 节点坐标 / `candidate_id` | builder 给候选的身份：`节点:起点:类型:REF>ALT`，坐标在节点的正向链上。跨几个节点的缺失写成 `节点:起点:DEL:REF>@n2+n3` |
| GRCh38 节点 | GRCh38 参考路径恰好经过一次的节点，它的每个碱基都有唯一的 GRCh38 坐标 |
| 分支节点 | 不在 GRCh38 路径上的节点（人群里其它单倍型的 allele），或者 GRCh38 经过不止一次的节点。这类节点没有唯一的 GRCh38 坐标 |
| truth 键 | 把一个 truth allele 的每个等价位置换算成的节点坐标写法，格式和 `candidate_id` 一样 |
| span（lo, hi） | 一个 truth allele 所有等价位置覆盖的 GRCh38 区间 |
| confident 区 | somatic BED ∩ germline BED（HG008：`HG008-T_somatic_smvar_benchmark_v0.2_all.bed` ∩ dipcall `dip.bed`） |
| PASS | VCF 的 FILTER 是 `PASS` 或 `.` |

标签值：**1** somatic，**2** germline，**0** non（确定不是真实变异的负样本），**−1** ignore（训练时丢掉）。

## 1. 整体流程

```
truth VCF ──(1) 拆分 allele、枚举等价位置、换算成节点键──▶ truth 表（键 → truth allele；span；PASS；in_bed）
                                                                │
tensor（variant_summary 的一行）──(2) 求位置：GRCh38 坐标，或分支节点的 anchor 区间 ──┐   │
                                                                                   ▼   ▼
                                                              (3) 按固定顺序判断，第一条符合的规则决定标签
```

召回表（`somatic.recall.tsv`）只看 truth 键有没有出现在某个 tensor 的 allele 或 filtered 候选里，**和标签规则无关**。

## 2. 第一步：准备 truth allele

### 2.1 拆分（`split_alleles`）

- 只取 chr1–22。每个 ALT 单独成为一个 allele；`<DEL>`、`*`、`.`、带 N 的跳过。
- 先去掉 REF/ALT 共同的尾部，再去掉共同的头部（记下起点）。
- 去完以后：
  - REF 和 ALT 等长：每个不同的碱基拆成一个 SNV。MNV `AC>GT` 就变成两个 SNV；
  - REF 为空：INS；ALT 为空：DEL；
  - 都不为空、长度不同：COMPLEX。COMPLEX 不可能等于任何候选，只用来判断"附近有 truth"。

例子：VCF `chr20:151663 TA>T` → 去掉共同头部 `T` → DEL，pos0 = 151663，REF = `A`。

### 2.2 等价位置（`placements`）

重复序列里的 indel 有很多等价写法（VCF 左对齐，vg 可能写在另一端）。代码先把 allele 左移到最左，再一步步右移，把每个位置都列出来（最多 5,000 个）：

```
GRCh38    C T T T T T G
truth     在 C 之后插入 T（VCF: C>CT）
等价位置  C^TTTTTG, CT^TTTTG, CTT^TTTG, …, CTTTTT^G   共 6 个插入点
```

DEL 同理：在 `TTTTT` 里删一个 T，5 个位置都算。SNV 只有一个位置。

### 2.3 节点键（`Locator.keys`）

每个等价位置都通过 GRCh38 路径换算成节点坐标：
- 节点在路径上是反向的，就把位置翻过来、碱基取反向互补；
- INS 落在两个节点的交界处，两个节点上各给一个键；
- DEL 跨过几个连续、同方向的 GRCh38 节点，给跨节点键 `…@n2+n3`；
- 路径经过不止一次的节点不给键。

`chr20:151663 TA>T` 得到的键是 `30192875:0:DEL:A>`，正好是 PacBio 上某个 tensor 的 `candidate_id`，所以那个 tensor 是"完全匹配"。

### 2.4 truth 表里记下的其它东西

- **span** `lo..hi`：所有等价位置覆盖的区间。INS 的 span 从最左的插入点到最右的插入点。
- **in_bed**：span 是否在这个 truth 自己的 BED 里（只记录，1/2 的判断不看它）。
- **passed**：FILTER 是 PASS 或 `.`。germline 另有一类：FILTER 只有 `GAP1`/`GAP2`（dipcall 在一条组装单倍型上调用了、另一条没组装出来），这类也算 germline truth。

## 3. 第二步：tensor 的位置

### 3.1 GRCh38 坐标（`ReferencePath.linear`）

- 候选所在节点是 GRCh38 节点，就换算出 `grch38 = {chrom, pos0, ref, alt, node_reverse}`（碱基转回 GRCh38 正链）。
- 跨节点的缺失要求这些节点都是 GRCh38 节点、方向一致、在参考上首尾相接，否则没有坐标。
- 判断时用的区间（`linear_interval`）：SNV/DEL 是 `[pos0, pos0 + len(REF))`；INS 是插入点两侧各 1 bp，即 `[pos0 − 1, pos0 + 1)`。

### 3.2 分支节点：anchor（`anchor`）

Minigraph-Cactus 按拓扑顺序给节点编号，所以分支节点的 ID 夹在它两侧参考节点的 ID 中间。代码据此给分支节点定位：

1. 在节点 ID 两侧各找 200 个 ID 以内、最近的 GRCh38 节点；
2. 两侧必须在同一条染色体上；
3. **v5 新增**：两个参考节点**之间**的 GRCh38 碱基（不算这两个节点本身）不能超过 1,024 bp。着丝粒附近按 ID 相邻的节点可能相距几 Mb，v4 在那里产生了错误的匹配；
4. 满足以上条件，anchor 区间就是"左侧节点起点 → 右侧节点终点"，写进 `details["anchor"]`。

例子：`31293240:0:INS:>A` 在分支节点上，anchor = `chr20:58872268–58872534`（266 bp）。

找不到或不满足条件的，没有位置（→ −1，见 R4）。

## 4. 第三步：分类规则（按顺序，第一条符合的决定标签）

| 顺序 | 条件 | 标签 | reason |
|---|---|---:|---|
| R1 | A1 的键等于一个 **PASS 的 somatic truth** | 1 | `representative_allele_is_somatic_truth` |
| R1' | A1 等于 somatic truth，但它们都没 PASS | −1 | `somatic_truth_filtered` |
| R2 | A1 等于一个 **PASS 的 germline truth** | 2 | `representative_allele_is_germline_truth` |
| R2' | A1 等于一个 FILTER 只有 GAP1/GAP2 的 germline truth | 2 | `representative_allele_is_germline_truth_gap_filtered` |
| R2'' | A1 只等于其它 FILTER 的 germline truth（如 HET1/HET2） | 先记下，继续往下 | 见 R6 |
| R3 | A1 不是 truth，但 **site 里别的 allele 是 PASS truth** | 重叠 >60% 取 truth 的标签，否则 0 | `allele_partial_somatic_truth` / `allele_partial_germline_truth` / `truth_matches_non_representative_allele` |
| R3' | site 里别的 allele 是没 PASS 的 somatic truth | −1 | `somatic_truth_filtered` |
| R4 | 没有 GRCh38 坐标，也没有 anchor | −1 | `not_on_unique_grch38_node` |
| R5 | 位置不在 confident 区里（anchor 用区间中点判断） | −1 | `outside_confident_region` |
| R6 | R2'' 的情况（在 confident 区里） | 0 | `germline_truth_filtered` |
| R7 | **同一位置有 PASS somatic truth**，A1 和它重叠 >60% | 1 | `residual_partial_somatic_truth` |
| R8 | （只限有 GRCh38 坐标的）同一位置有 PASS germline truth，A1 和它重叠 >60% | 2 | `residual_partial_germline_truth` |
| R9 | 10 bp 内有任何 truth allele（somatic 或 germline） | 0 | `near_truth_allele_mismatch` |
| R10 | 其它所有 | 0 | `confident_no_truth_allele` |

说明：
- **1 和 2 不看 BED**（R1–R3 都在 R5 之前）。0 必须在 confident 区里；−1 只剩 BED 外、没有位置、somatic truth 没 PASS。
- **somatic 优先**：A1 同时等于 somatic 和 germline truth 时，R1 先命中，标 1。
- R3 里先试 somatic，再试 germline。
- R7/R8 的"同一位置"：truth 的 span 和 tensor 的区间相交或紧挨着（有坐标时两边各放宽 1 bp，anchor 放宽 10 bp）。
- R8 只对有 GRCh38 坐标的 tensor 做：**分支节点上的 tensor 永远拿不到 2**。
- R9 的"10 bp 内"用的是所有 truth allele，包括没 PASS 的和 COMPLEX。

### 各规则的真实例子（PacBio v6）

**R1 完全匹配 somatic。** `30217844:21:SNP:G>A` = somatic `chr20:964493 G>A`（0/1，PASS）。67 条 reads 全是 ALT → **1**。

**R2 完全匹配 germline。** `30192018:0:SNP:C>T` = germline `chr20:101582 C>T`（0/1，PASS）→ **2**。INDEL 例子：`30192875:0:DEL:A>` = `chr20:151663 TA>T`（1/1）→ **2**。

**R2' GAP 过滤的 germline。** `30749403:652:DEL:A>` = germline `chr20:26689873 CA>C`，FILTER `GAP2`，GT `1/.`，在 dip.bed 外 → **2**（不看 BED）。

**R3 部分匹配 site 里的另一个 allele。**
- `30246963:0:INS:>TT`：site 里有 `+TT`（13 条 reads，A1）、`+T`（10）、`+TTT`（10）。germline truth `chr20:2077646 A>ATTT`（2/1）是 `+TTT`，不是 A1。`+TT` 和 `+TTT` 的重叠 = 2/3 = 0.667 > 0.6 → **2**（`allele_partial_germline_truth`）。
- `30941026:157:INS:>AA`：A1 是 `+AA`（46 条），somatic truth `chr20:42232813 C>CA` 是 site 里的 `+A`（26 条）。allele 重叠 = 1/2 = 0.5 不够；somatic 再算 A1 reads 的单倍型重叠，得到 1.0 → **1**（`allele_partial_somatic_truth`）。
- `30200048:135:INS:>T`：A1 `+T`（28 条），germline truth `chr20:385860 A>ATT`（0/1）是 site 里的 `+TT`（5 条）。重叠 = 1/2 = 0.5；germline 不算单倍型重叠 → **0**（`truth_matches_non_representative_allele`）。

**R4 没有位置。** `30211411:0:INS:>CC` 在分支节点上，两侧 200 个 ID 内没有合格的参考节点 → **−1**。

**R5 confident 区外。** `30191551:9:DEL:T>` 在 `chr20:60300`，不在 somatic BED ∩ dip.bed 里 → **−1**。分支节点也一样：`30191605:101:DEL:C>` 的 anchor `chr20:61453–61534` 在区外 → **−1**。

**R7 部分匹配同一位置的 somatic truth（残余 edit）。**
- `31293240:0:INS:>A`（miss 分析里的例子）：somatic truth `chr20:58872348 T>TAA`。ALT reads 走一个多一个 A 的分支节点，再带一个 `+A` 残余 edit。这个 tensor 在分支节点上，anchor `chr20:58872268–58872534` 覆盖 truth 位置。A1 reads 的单倍型重叠 = 0.75 → **1**。
- `30602891:0:SNP:T>C`：somatic truth 是 `chr20:19677686 CTTTT>C`（删 4 个 T），但 reads 在旁边写出了一个 SNV。单倍型重叠 1.0 → **1**。所以 SNV tensor 也可以部分匹配一个 DEL truth。

**R8 部分匹配同一位置的 germline truth。** `30209021:0:DEL:TTTTTTT>@30209022`（删 7 个 T，跨两个节点），germline truth `chr20:775268 CTTTTTTTTTTT>CTTT`（删 8 个 T，1/2）。共同删掉的碱基 = 7，除以较长的 8 = 0.875 → **2**。

**R9 附近有 truth 但对不上。** `30191903:0:INS:>T`（`chr20:92969`，14 条 reads，AF 0.25），10 bp 内有 truth allele，但既不相同也重叠不够 → **0**。

**R10 附近什么都没有。** `30191732:122:INS:>A`（`chr20:73326`，5 条 reads，AF 0.09）→ **0**。分支节点也会走到这里：`30194624:0:SNP:T>C` 在 anchor `chr20:214823–214878`，7 条 reads 全是 ALT（AF 1.0）→ **0**。

## 5. 重叠怎么算

`MIN_OVERLAP = 0.6`，必须**大于** 0.6。有两种算法，取较大的那个。

### 5.1 allele 重叠（`allele_overlap`，somatic 和 germline 都用）

只比较同类型的 INDEL；SNV 对 SNV（ALT 不同）或类型不同，重叠都是 0。

- **DEL 对 DEL**：两者都删掉的 GRCh38 碱基数 ÷ 较长的那个。truth 的位置取整个 span（任何等价位置都行）。
  - 例子（R8）：删 7 个 T 对删 8 个 T，7/8 = 0.875 ✔。
  - 删 1 个 T 对删 2 个 T，1/2 = 0.5 ✘。
- **INS 对 INS**：插入点必须在 truth 的 span 里。truth 的插入序列按插入点在重复里的位置轮转后，和 tensor 的插入序列求最长公共子序列（LCS），÷ 较长的那个。
  - `+TT` 对 `+TTT`：LCS 2 ÷ 3 = 0.667 ✔。
  - `+T` 对 `+TT`：1 ÷ 2 = 0.5 ✘。
  - `+AAATA` 对 germline `chr20:2284361 A>AAAAT`（插入 `AAAT`，插入点轮转后对齐）：LCS 4 ÷ 5 = 0.8 ✔（PacBio `30251303:0:INS:>AAATA`，标 2）。

### 5.2 单倍型重叠（`haplotype_overlap`，只对 somatic truth）

用 tensor 里 A1 那组 reads 的实际序列来判断它们带的是不是 truth 的变化，所以"走分支 + 残余 edit""跳过节点 + 错配"这类写法也能认出来。

1. 从 tensor 的 channel 0 取最多 20 条 A1 行、10 条 REF 行，每行至少 30 个碱基（tensor 宽 101 列，所以每条大约 100 bp 的局部序列）。
2. 取 truth 两侧各 90 bp 的 GRCh38 序列作 REF 单倍型，在其中换入 truth allele 得到 truth 单倍型。
3. 每条 read 分别算它嵌进 REF 单倍型、嵌进 truth 单倍型所需的编辑数：`d_ref`、`d_truth`。分支节点的 read 不知道方向，正反两个方向取较好的。
4. 从 `d_ref`、`d_truth` 里减去 REF reads 编辑数的中位数 `b`，算作测序背景噪声。
5. `event` = truth 的长度（SNV 为 1，INDEL 为 max(len REF, len ALT)）。
6. 每条 read：`shared = (d_ref + event − d_truth) / 2`，`overlap = shared / max(d_ref, event)`。
7. 所有 A1 reads 取中位数；少于 3 条 A1 reads 时不算。

用 truth `+AA`（event = 2）、背景 b = 0 来算几种 read：

| read 带的变化 | d_ref | d_truth | shared | overlap |
|---|---:|---:|---:|---:|
| 正好 `+AA` | 2 | 0 | 2 | **1.0** |
| `+AAA` | 3 | 1 | 2 | 0.667 |
| `+A` | 1 | 1 | 1 | 0.5 |
| 旁边一个无关的错配 | 1 | 3 | 0 | 0 |

A1 reads 的中位数超过 0.6，就算部分匹配这个 somatic truth。

## 6. 输出

每个 `v*_tensors/<SNV|INDEL>/` 里：
- `<chrom>_shard_NNNNN_labels.npy`（int8，和 tensor 一一对应）。
- `<chrom>_labels.ndjson`，每个 tensor 一行：`candidate_id`、`site_id`、`label`、`label_name`、`reason`、`somatic`/`germline`（键完全相同的 truth，含 `representative`，表示是不是 A1）、`grch38`、`anchor`、`partial`（`allele` / `residual` / null）、`overlap`、`partial_truth`。
- `labels.manifest.json`：版本、每条染色体和每个 reason 的计数、`partial` 计数、truth 和 BED 的 SHA-256。

## 7. 三个平台的计数（v5）

| 标签 | reason | PacBio SNV | PacBio INDEL | ONT SNV | ONT INDEL | Illumina SNV | Illumina INDEL |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | representative_allele_is_somatic_truth | 8,234 | 4,097 | 8,275 | 3,673 | 8,279 | 3,562 |
| 1 | allele_partial_somatic_truth | 0 | 66 | 0 | 133 | 0 | 14 |
| 1 | residual_partial_somatic_truth | 374 | 2,008 | 305 | 1,505 | 1,424 | 962 |
| 2 | representative_allele_is_germline_truth | 332,598 | 147,508 | 343,284 | 157,914 | 311,650 | 55,800 |
| 2 | representative_allele_is_germline_truth_gap_filtered | 18,229 | 1,543 | 26,891 | 2,177 | 5,145 | 373 |
| 2 | allele_partial_germline_truth | 0 | 1,666 | 0 | 2,235 | 0 | 151 |
| 2 | residual_partial_germline_truth | 0 | 5,124 | 0 | 6,715 | 0 | 5,608 |
| 0 | confident_no_truth_allele | 69,253 | 1,325,307 | 272,721 | 1,711,480 | 3,089,145 | 7,737 |
| 0 | near_truth_allele_mismatch | 66,466 | 315,924 | 89,452 | 421,759 | 898,940 | 50,219 |
| 0 | truth_matches_non_representative_allele | 138 | 6,477 | 163 | 34,483 | 334 | 973 |
| −1 | outside_confident_region | 652,004 | 98,339 | 1,594,086 | 156,864 | 661,461 | 25,093 |
| −1 | not_on_unique_grch38_node | 33,043 | 8,631 | 37,329 | 12,075 | 24,510 | 973 |
| | 总数 | 1,180,339 | 1,916,690 | 2,372,506 | 2,511,013 | 5,000,888 | 151,465 |

`somatic_truth_filtered` 和 `germline_truth_filtered` 在 HG008 三个平台上都是 0（somatic truth 全部 PASS；没有 tensor 的 A1 只匹配到 HET 等其它过滤的 germline allele），表里没列。

## 8. 版本历史

| 版本 | 主要变化 |
|---|---|
| v1 | 只有完全匹配：A1 = somatic → 1；A1 = germline（PASS 且在 dip.bed 内）→ 2；confident 区内、10 bp 内没有 truth → 0；其它 −1。INS 的"附近"范围后来修正为整个等价区间（commit `9b8b1b3`） |
| v2 | 1/2 只要求 PASS，不看 BED；没 PASS → −1 |
| v3 | "附近有 truth 但对不上" −1 → 0；GAP1/GAP2 的 germline → 2；其它过滤的 germline → 0 |
| v4 | 部分匹配（重叠 >60%）拿 truth 的标签（allele / residual 两种；somatic 另算单倍型重叠）；分支节点按参考邻居定位；−1 只剩 BED 外、没位置、somatic 没 PASS |
| v5 | 分支节点两侧参考节点之间的 GRCh38 碱基超过 1,024 bp 就不定位（→ −1） |

## 9. 读代码时注意到、值得讨论的点

以下都是按现在代码的实际行为整理的，数字来自 PacBio `v6_tensors`（v5）。

1. **分支节点上的 tensor 永远拿不到 2。** R8（germline 部分匹配）只对有 GRCh38 坐标的 tensor 做。分支节点上 AF ≥ 0.3、却标 0 的：INDEL 41,614 个，SNV 33,976 个。这里面应该有不少是人群分支上的真实 germline allele（例如上面 `30194624:0:SNP:T>C`，7/7 条 reads 是 ALT）。标 1 的分支节点 tensor 只有 INDEL 950、SNV 221。
2. **truth 在 site 的次要 allele 上、重叠不到 60% → 0。** INDEL 6,477 个，其中 3,998 个的 truth allele reads 不少于 A1 的一半。例子 `chr20:385860`：A1 `+T` 28 条，germline truth `+TT` 5 条 → 0。site 里确实有真实 germline allele 的 reads，却当成负样本。
3. **重叠阈值对短 indel 很严。** `+T` 对 `+TT`、删 1 个对删 2 个的重叠都是 0.5，永远过不了 0.6；`+TT` 对 `+TTT` 是 0.667 就能过。同聚物里差一个碱基，短的算"不同"，长的算"相同"。
4. **SNV 对 SNV 只看完全相同。** 同一位置 ALT 不同（例如 truth G>A、A1 G>C）的 allele 重叠为 0；somatic 还有单倍型重叠兜底，germline 没有。
5. **germline 不做单倍型重叠。** 所以"走分支 + 残余 edit"写法的 germline 变异，只能靠 R8 的 allele 重叠（且只在 GRCh38 节点上），多数落到 R9 的 0。
6. **somatic 优先。** A1 同时等于 somatic 和 germline truth 时标 1。HG008 有 209 个 somatic allele 和 germline 完全相同（多为正常样本杂合）。
7. **"附近有 truth 但对不上"一律 0。** 其中包括 somatic 单倍型的残余 edit、重叠不到 60% 的：miss 分析里 I1–I3 未命中的残余 edit，标 0 的有 PacBio 1,031、ONT 1,254、Illumina 2,729 个（v4 统计，v5 只改了分支节点定位）。
8. **A1 正好等于 germline truth、但 reads 其实是 somatic 单倍型的残余。** 例如 germline `TA>T`、somatic `TAA>T`，残余 edit `DEL A` 等于 germline 键 → 2。PacBio 的 330 个"残余 edit 标 2"里，这类约 140 个；另有约 87 个是 ALT reads 上附带的真实 germline 变异（这些标 2 是对的）。
9. **1/2 不看 BED，0 要看。** confident 区外的 truth 匹配照样是 1/2，而同样位置的非 truth tensor 是 −1。分支节点判断是否在 confident 区只看 anchor 区间的中点。
10. **单倍型重叠只用 tensor 窗口里的序列。** 每条 read 大约 100 bp、最多 20 条，取中位数；A1 少于 3 条时不算。长插入（>50 bp）或跨出窗口的事件算不准。
