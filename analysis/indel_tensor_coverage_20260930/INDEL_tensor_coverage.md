# INDEL truth VCF 中的 tensor 覆盖情况（HG008T Illumina，按 graph 重新统计）

更新：2026-09-30。"未生成 candidate"的 I1 / I2 不再只按 reads 判定，而是在 d9 graph 上按序列重新检查（方法见最后一节）。
"有 tensor"和"被 builder 阈值过滤"按 allele key 和阈值判定，不受影响。

## 总览

INDEL truth VCF 共包含 **9,343 条 INDEL 记录**。其中 **862 条位于 chrX/Y**；当前 tensor 仅针对 chr1–22 生成，因此这些记录没有 tensor。

chr1–22 上的记录按多等位位点拆分后，共得到 **8,496 个 INDEL allele**。下表及后续原因统计均以 allele 为单位。

| 状态 | 数量 | 占 chr1–22 allele 比例 |
| --- | ---: | ---: |
| 有 tensor | 3,633 | 42.8% |
| ├─ 作为 A1 | 3,555 | — |
| └─ 作为同一位点的其他 allele | 78 | — |
| 无 tensor：已生成 candidate，但被 builder 阈值过滤 | 909 | 10.7% |
| ├─ ALT reads 少于 3 条 | 637 | — |
| ├─ AF 低于 0.08 | 239 | — |
| └─ 两项阈值均不满足 | 33 | — |
| 无 tensor：未生成 candidate | 3,954 | 46.5% |
| **合计** | **8,496** | **100%** |

> 注：9,343 和 862 的单位是 VCF 记录；8,496 及表中各项的单位是拆分后的 allele，因此不能直接用前两项相减来核对表格合计。

## chr1 BED 测试区域

测试用的 chr1 BED 区域内共有 **588 个 INDEL allele**。

| 状态 | 数量 | 占 chr1 BED allele 比例 |
| --- | ---: | ---: |
| 有 tensor | 239 | 40.6% |
| ├─ 作为 A1 | 236 | — |
| └─ 作为同一位点的其他 allele | 3 | — |
| 无 tensor：已生成 candidate，但被 builder 阈值过滤 | 66 | 11.2% |
| ├─ ALT reads 少于 3 条 | 46 | — |
| ├─ AF 低于 0.08 | 19 | — |
| └─ 两项阈值均不满足 | 1 | — |
| 无 tensor：未生成 candidate | 283 | 48.1% |
| **合计** | **588** | **100%** |

因此，在**仅以现有 tensor 为输入**、并以这 588 个 allele 为评估集合的条件下，chr1 BED 区域 INDEL 模型的理论 recall 上限为 **239 / 588 = 40.6%**。

## 未生成 candidate 的原因（按 graph 重新分组）

"未生成 candidate"指 reads 在 graph 上比对时，没有产生与 truth allele 一致的 edit。原来的 I1 / I2 是按 reads 判定的；现在在 d9 graph
上按序列重新判定，分成 a–d 四组。I3–I7 是 read 层面的原因，没有重新分析。

| 组 | 含义 | chr1–22 | chr2–22 | chr1 BED |
| --- | --- | ---: | ---: | ---: |
| **a** | **truth allele 就是 graph 中的一条 path**：reads 沿这条 path 走就没有 edit，结构上不可能产生 candidate | **2,095** | **1,936** | **146** |
| **b** | **graph 有长度相近的 allele，reads 走它，残余 edit 补齐后正好等于 truth**（真正的 I2，会产生 partial tensor） | **944** | **869** | **72** |
| c | graph 只有相近的 allele（或只有 GRCh38），reads 跟着它走，reads 上看不到 truth allele | 299 | 279 | 19 |
| d | 残余 edit 不能把任何 graph allele 变成 truth，或 reads 走 GRCh38 且单靠残余 edit 就等于 truth | 82 | 76 | 6 |
| I4 | ALT reads 少于 3 条 | 306 | 285 | 19 |
| I3 | reads 留在 GRCh38 上，但拼出的 edit 对应另一个 allele | 106 | 94 | 12 |
| I5 | 该位置没有 MAPQ > 10 的 read 覆盖 | 101 | 91 | 8 |
| I6 / I7 | INDEL 超过 50 bp，或原因不明 | 21 | 20 | 1 |
| **合计** | | **3,954** | **3,650** | **283** |

和原来 I1 / I2 的对应关系：

| 原分类 | → a | → b | → c | → d | 合计 |
| --- | ---: | ---: | ---: | ---: | ---: |
| I1（chr1–22） | 1,803 | — | 113 | — | 1,916 |
| I2（chr1–22） | 292 | 944 | 186 | 82 | 1,504 |
| I1（chr1 BED） | 127 | — | 9 | — | 136 |
| I2（chr1 BED） | 19 | 72 | 10 | 6 | 107 |

- a + b 合计 **3,039 / 3,954（76.9%）**，chr2–22 **2,805 / 3,650（76.8%）**，chr1 BED **218 / 283（77.0%）**：truth allele 全部或部分已由 graph 路径表示，这是未生成
  candidate 的主要原因。
- 原 I2 中只有 **944 / 1,504（63%）** 是真正的"graph allele + 残余 edit = truth"；292 个的 truth 其实已在 graph 中（应属 a），
  186 + 82 个的 reads 上看不到 truth，或残余 edit 和 graph 对不上。
- 原 I1 中 113 个（5.9%）graph 里只有相近的 allele，不是 truth 本身。

## I1：reads 的路径（read 层面，原 I1 全部 1,916 个，chr1–22）

在原 I1 类别中，INS 主要通过 graph branch 表示，DEL 主要通过 skip edge 表示。

| INDEL 类型 | I1 总数 | 走 branch | 走 skip edge |
| --- | ---: | ---: | ---: |
| INS | 1,136 | 1,129 | 7 |
| DEL | 780 | 104 | 676 |
| **合计** | **1,916** | **1,233** | **683** |

Insertion：插入的序列如果正好是旁边 GRCh38 序列的重复（tandem duplication / STR 多一个拷贝），graph 可以不另建 branch node，而是加一条往回的
edge，形成一个环。haplotype 把同一段 GRCh38 node 再走一遍，就多出一个拷贝。read 照着走，同样没有 edit；但它的 GRCh38 坐标出现回退、重复，
分类器就把它也记成了 "skip_edge"。

## I1：graph 拼出 truth ALT 序列的方式

### chr1 BED 内

| graph 路径形式 | INS | DEL |
| --- | ---: | ---: |
| 单个 bubble node：一个 branch node 夹在两个相邻的 GRCh38 node 之间，其序列就是插入序列 | 50 | — |
| 多个 branch node：纯插入，但插入序列被拆成多个 branch node（例如 +AA 被拆成两个 "A" node） | 24 | — |
| 单条 skip edge：只经过 GRCh38 node，一条 edge 正好跳过被删除的碱基 | — | 39 |
| 多处偏离：需要组合多个 bubble 或 skip edge 才能拼出 ALT，多见于较长的 STR | 4 | 7 |
| 加上附近的 germline SNV 后才是 graph path | 3 | — |
| 最接近的 graph allele 差 1 bp | 2 | 3 |
| 最接近的 graph allele 差 ≥ 2 bp | 4 | — |
| **合计** | **87** | **49** |

BED 内合计 **136 个 allele**。前五行（127 个）的 truth allele 就是 graph path（a 组）；后两行（9 个）graph 里只有相近的 allele（c 组）。

### 为什么多处 deletion 也能拼出同一个 ALT？

#### Homopolymer 例子

GRCh38 为 `A TTTTTT G`：共有 6 个 T，每个 T 是一个 node。Truth 删除 2 个 T，ALT 为 `A TTTT G`。

1. 一次删：跳过相邻的两个 T，得到 `A TTTT G`。
2. 分两处删：`A → 跳过 T1 → T2 → T3 → 跳过 T4 → T5 → T6 → G`，同样得到 `A TTTT G`。

第二条路径在两次跳过之间回到 GRCh38，经过 T2、T3。由于这些 T 与被跳过的 T1、T4 在序列上相同，跳过哪两个 T 都会得到同一个 ALT。

#### STR 例子

GRCh38 为 `C TA TA TA TA TA G`：共有 5 个 TA。Truth 删除 2 个 TA，ALT 为 `C TA TA TA G`。

1. 一次删：连续跳过两个 TA，得到 `C TA TA TA G`。
2. 分两处删：跳过第 1 个 TA，经过第 2、3 个，再跳过第 4 个，经过第 5 个；同样得到 `C TA TA TA G`。

重复序列允许 graph 路径在不同位置跳过相同的碱基或重复单元。即使路径多次离开、回到 GRCh38，最终拼出的 ALT 序列也可能完全相同。

### chr1–22（表中的 chr1 是整条染色体，和上一小节的 chr1 BED 内子集不同）

INS：

| 类别 | chr1 | chr2–22 | 合计 |
| --- | ---: | ---: | ---: |
| 单个 bubble node：插入序列正好是一个 variant node | 53 | 579 | 632（55.6%） |
| 多个 branch node：在一条 branch path 上 perfect align | 25 | 272 | 297（26.1%） |
| Branch 替换了一段 GRCh38 | 0 | 4 | 4（0.4%） |
| 多处偏离：需要组合多个 bubble 或 skip edge | 5 | 93 | 98（8.6%） |
| 加上附近的 germline SNV 后才是 graph path | 3 | 39 | 42（3.7%） |
| 最接近的 graph allele 差 1 bp | 2 | 24 | 26（2.3%） |
| 最接近的 graph allele 差 ≥ 2 bp | 5 | 32 | 37（3.3%） |
| **合计** | **93** | **1,043** | **1,136** |

DEL：

| 类别 | chr1 | chr2–22 | 合计 |
| --- | ---: | ---: | ---: |
| 单条 skip edge | 44 | 599 | 643（82.4%） |
| Branch 替换了一段 GRCh38 | 0 | 2 | 2（0.3%） |
| 多处偏离 | 8 | 69 | 77（9.9%） |
| 加上附近的 germline SNV 后才是 graph path | 0 | 8 | 8（1.0%） |
| 最接近的 graph allele 差 1 bp | 3 | 30 | 33（4.2%） |
| 最接近的 graph allele 差 ≥ 2 bp | 0 | 17 | 17（2.2%） |
| **合计** | **55** | **725** | **780** |

## I2：graph 对残余 edit 的影响

### 1. 原 I2 是怎么判定的（read 层面）

Miss 分析中的 `reads.py` 根据 reads 判定 I2；这个分类本身不是对拼接后序列的证明。判定分三步：

1. 挑 ALT reads：比较两端 anchor 之间的 read 序列与 ALT、REF haplotype 的编辑距离。如果与 ALT 的距离严格小于与 REF 的距离，就记为 ALT
   read（`alt_like`）；与 ALT 完全一致的另记为 `alt_exact`。
2. 查看 reads 在 site 的路径和 edit：判断 read 是绕开 GRCh38（走 branch 或 skip edge），还是在 site 范围内（包括两侧各 5 bp）带有任何
   edit；SNV edit 也计入。
3. 按多数 reads 定类：绕开 GRCh38 且没有 edit 的归为 I1；绕开 GRCh38 且有 edit 的归为 I2。

因此，这套判定没有检查"graph allele + 残余 edit"拼出的序列是否等于 truth。可能出现两类误判：

- 在较长的 homopolymer / STR 中，同一位置长度相近的其他 allele 的 reads 也可能因为"离 ALT 更近"而被算作 ALT reads。
- 测序错误产生的 SNV edit 也可能使一个 truth allele 被归入 I2。

### 2. 按 graph 重新分类

对每个原 I2 truth：枚举 graph 在附近能拼出的所有 haplotype，找出哪一条加上 reads 上最主要的残余 INDEL 之后，序列正好等于 truth ALT。

原 I2 在 chr1–22 上共有 **1,504 个 allele**（INS 1,056，DEL 448），其中 **107 个**位于 chr1 BED 内（INS 75，DEL 32）。

| 组 | graph 的情况 | chr1–22 | chr2–22 | chr1 BED |
| --- | --- | ---: | ---: | ---: |
| b | graph 有同方向、更短的 allele；reads 走它，残余 edit 补上剩下的部分 | 700（46.5%） | 644 | 54 |
| b | graph 有同方向、更长的 allele；reads 走它，残余 edit 往回补 | 240（16.0%） | 221 | 18 |
| b | graph allele 方向相反（例如 truth 是 DEL，reads 走 graph 的 INS 分支）+ 残余 edit = truth | 4（0.3%） | 4 | 0 |
| a | graph 中已有 truth 这条 path，但 reads 没有干净地走它 | 292（19.4%） | 271 | 19 |
| c | 没有残余 INDEL：reads 跟着 graph 里另一个 allele 走 | 174（11.6%） | 164 | 10 |
| c | 没有残余 INDEL：最接近的就是 GRCh38 | 12（0.8%） | 12 | 0 |
| d | 残余 edit 不能把任何 graph allele 变成 truth | 65（4.3%） | 60 | 5 |
| d | reads 走 GRCh38，单靠残余 edit 就等于 truth（graph allele 没有用到） | 17（1.1%） | 16 | 1 |
| | **合计** | **1,504** | **1,392** | **107** |

chr1 BED 内的例子：

| graph 的情况 | 例子 |
| --- | --- |
| graph allele 更短 | chr1:13389005 truth −3T：graph 只有 −1T（skip edge），reads 走 −1T，再加 −2T 残余（81 条 spanning reads 中 76 条） |
| graph allele 更短 | chr1:72386968 truth −13T：graph 只有 −4T、−5T 两条 skip edge，reads 走 −5T，再加 −8T 残余（43 条 spanning reads 中 40 条） |
| graph allele 更长 | chr1:18232209 truth +12T：graph 有 +11T 和 +13T，reads 走 +13T，再加 −1T（47 条中 36 条） |
| graph 已有 truth | chr1:64132760 truth −3A：graph 有 −3A，但 29 条 reads 中没有一条和 truth 完全一致 |
| reads 跟着另一个 allele | chr1:167357369 truth −10C：graph 有 −10C 和 −9C，reads 停在 −9C 上，没有残余 INDEL |
| 残余 edit 对不上 | chr1:210311120 truth −8C：graph 有 −6、−5、−3、−2 等，reads 只带 −1C |

graph 的作用：这些位置几乎都在 homopolymer 或 STR 中，HPRC 在那里有一串长度不同的 allele（skip edge 或 branch）。aligner 选最接近
reads 实际长度的现成路径，把剩下的差额写成更短的 edit（gap 更短、得分更高）。于是 truth 被拆成"graph allele + 残余 edit"；builder
只按残余 edit 生成 candidate，它的 allele 与 truth 对不上，只能成为 partial tensor（b 组）。

## 方法

- graph：HPRC v1.1 d9（`hprc-v1.1-mc-grch38.d9`）。node 序列和 edge 取自
  `AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13.path_index.sqlite`；GRCh38 坐标取自 pipeline 的
  `pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path`（sqlite 的 `path_coords` 缺 chr3–9 和部分 chr22 的 GRCh38）。
- 每个 truth：取 truth 两侧各约 60 bp 范围内的 GRCh38 node，枚举从第一个到最后一个 GRCh38 node 之间 graph 能拼出的所有序列（I1 上限
  20,000 条 path，I2 上限 3,000 条 haplotype），与 truth ALT（GRCh38 + truth allele）比较。
- I1：先找与 truth ALT 完全一致的 path 并按离开 GRCh38 的方式分类；找不到时加上附近最多 6 个 PASS germline allele 再找；仍找不到则记最接近
  graph haplotype 的编辑距离。
- I2：reads 的残余 INDEL 取 miss 表（`analysis/hg008_somatic_miss_20260926/illumina/indel_no_candidate.tsv`）中支持 reads 最多的
  INDEL edit；找出加上它后等于 truth ALT 的 graph haplotype，比较该 graph allele 与 truth 的长度和方向。
- I3–I7、"有 tensor"、"被阈值过滤"沿用原统计，未重新分析。
- 脚本和逐 truth 结果在本目录：`classify.py` → `classes_all.tsv`（原 I1，1,916 个），`i2_graph.py` → `i2_classes.tsv`（原 I2，1,504 个）。
  脚本在 `/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/` 运行（输出写到那里），这里是副本。
