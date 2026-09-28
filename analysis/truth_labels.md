# tensor 标签算法详解

- 代码：`indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py`。规则就是这个文件的模块 docstring 和常量；`labels.manifest.json` 里的 `rules_sha256` 是这个文件的 SHA-256，用来标识打标签时用的规则。
- 坐标换算：`indexed_gam_pipeline_v4/tensor_postprocessing/reference_path.py`（`ReferencePath.linear`）。
- 命令（在仓库根目录运行）：`python -m indexed_gam_pipeline_v4.tensor_postprocessing label --tensors … --reference-path … --fasta … --somatic-vcf … --somatic-bed … --germline-vcf … --germline-bed … --truth-dir … [--recall-dir …] [--snv-min-af 0.07] [--indel-min-af …]`。`orchestrate finalize` 在 merge 之后用同样的规则打标签，AF 下限（floor）用 `orchestrate prepare` 时给的 `--label-snv-min-af` / `--label-indel-min-af`（存进 `config.json` 的 `postprocess.labels.snv_min_af` / `indel_min_af`，不设是 null，即不按 AF 去掉），等于 `label` 的 `--snv-min-af` / `--indel-min-af`。短读长数据集在 prepare 时加 `--label-snv-min-af 0.07`，跑完不用再单独 relabel。这两个和 build 的 `--snv-min-af` / `--indel-min-af`（决定哪些 allele 生成 tensor）无关。
- 例子：除非另外注明，都取自 HG008 PacBio `/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/v3_tensors`（2026-09-28 用现行规则打的标签），按 `candidate_id` 在 `{SNV,INDEL}/chr*_labels.ndjson` 里查到的标签、reason、`partial`/`overlap`；reads 数来自同一目录的 `chr*_variant_summary.ndjson`；单倍型重叠例子里每行 read 的 `(d_ref, d_truth)` 和 `b` 是用 `truth_labels.py` 的函数在同一个 tensor 上重算的（结果和 labels.ndjson 的 `overlap` 一致）。第 7 节的计数来自六个数据集的 `labels.manifest.json`。

## 0. 名词

| 名词 | 意思 |
|---|---|
| tensor / site | 一个 site 是同一节点、同一起点、同一类型（SNV 或 INDEL）的一组候选 allele，只画成一个 tensor |
| A1 / 代表 allele | site 里 ALT reads 最多的 allele。tensor 的 `candidate_id` 就是它，画在最上面那组行（A1 的 row group） |
| 节点坐标 / `candidate_id` | builder 给候选的身份：`节点:起点:类型:REF>ALT`，坐标在节点的正向链上。跨几个节点的缺失写成 `节点:起点:DEL:REF>@n2+n3` |
| edit / 残余 edit | edit 是一条 read 相对它所走节点的一处差异（错配、插入、缺失），builder 从 edits 生成候选。残余 edit 是 ALT reads 走了图里已有的分支或跳过边之后，还剩下的、写法和 truth 不同的那部分差异 |
| GRCh38 节点 | GRCh38 参考路径恰好经过一次的节点，它的每个碱基都有唯一的 GRCh38 坐标 |
| 分支节点 | 不在 GRCh38 路径上的节点（人群里其它单倍型的 allele），或者 GRCh38 经过不止一次的节点（HPRC v1.1 d9 图里没有后一种：49,092,514 个 GRCh38 节点都只经过一次）。这类节点没有唯一的 GRCh38 坐标 |
| 等价位置（placement） | 重复序列里同一个 indel 的每一种写法（第 2.2 节） |
| truth 键 | 把一个 truth allele 的每个等价位置换算成的节点坐标写法，格式和 `candidate_id` 一样 |
| span（lo, hi） | 一个 truth allele 所有等价位置覆盖的 GRCh38 区间 |
| PASS | VCF 的 FILTER 是 `PASS` 或 `.` |
| 合格 truth | somatic：FILTER 是 PASS 或 `.`；germline：FILTER 是 PASS 或 `.`，或者只有 `GAP1`/`GAP2`（第 2.4 节） |
| confident 区 | somatic BED ∩ germline BED（HG008：`HG008-T_somatic_smvar_benchmark_v0.2_all.bed` ∩ dipcall `HG008N_GRCh38_dipcall.dip.bed`；COLO829T：`SMaHT_v2_easy_difficult_extreme.union.bed` ∩ COLO829BL dipcall `dipcall_hg38.dip.bed`） |
| anchor | 分支节点的估计位置：按节点 ID 找到的两侧参考节点（连同这两个节点本身）覆盖的 GRCh38 区间（第 3.2 节），写在 labels.ndjson 的 `anchor` 里 |
| partial | 部分匹配：A1 和 truth 不完全相同，但重叠超过 `MIN_OVERLAP`，拿 truth 的标签。`allele` = truth 是同一 site 里的另一个 allele；`residual` = truth 在同一位置、但不是 site 里的 allele（多是残余 edit） |

标签值：**1** somatic，**2** germline，**0** non（确定不是真实变异的负样本），**−1** ignore（训练时丢掉）。

规则里用到的常量（都在 `truth_labels.py` 开头）：

| 常量 | 值 | 用在哪 |
|---|---:|---|
| `NEAR_BP` | 10 | "附近有 truth"的距离（R11）；分支节点找同一位置 truth 时 anchor 两边放宽的距离（R8） |
| `MIN_OVERLAP` | 0.45 | 部分匹配要求重叠**大于**这个值（R4、R8、R10） |
| `ANCHOR_REACH` | 200 | 分支节点两侧各找多少个节点 ID（R5） |
| `ANCHOR_GAP` | 1,024 | 两侧参考节点之间最多隔多少 GRCh38 碱基（R5） |
| `EVIDENCE_ROWS` | 20 | 单倍型重叠用的 A1 行数上限（REF 行取一半，10） |
| `MIN_EVIDENCE_BASES` | 30 | 一行 read 至少要有的碱基数 |
| `HAPLOTYPE_WINDOW` | 90 | 单倍型重叠时 truth 两侧取的 GRCh38 碱基数 |
| `MAX_SHIFTS` | 5,000 | 一个 truth allele 最多枚举的等价位置数 |
| `--snv-min-af` | 不设（短读长数据集 0.07） | 命令行参数，不是常量（`orchestrate prepare` 里是 `--label-snv-min-af`）；设了以后 AF 更低的 SNV tensor 是 −1（R1） |
| `--indel-min-af` | 不设（六个数据集都没用） | 命令行参数，不是常量（`orchestrate prepare` 里是 `--label-indel-min-af`）；设了以后 AF 更低的 INDEL tensor（A1 是 INS 或 DEL）是 −1（R1），和 SNV 的下限分开设、互不影响 |

## 1. 整体流程

```
truth VCF ──(1) 拆分 allele、枚举等价位置、换算成节点键──▶ truth 表（键 → truth allele；span；passed；in_bed）
                                                                │
tensor（variant_summary 的一行）──(2) 求位置：GRCh38 坐标，或分支节点的 anchor 区间 ──┐   │
                                                                                   ▼   ▼
                                  (3) 按固定顺序判断（R1–R12），第一条符合的规则决定标签
```

实际顺序上，GRCh38 坐标一开始就算（R4 的 allele 重叠要用），anchor 只在走到 R5 时才算：R1–R4（AF、A1 是不是 truth、site 里别的 allele 是不是 truth）不需要 anchor。

召回表（`somatic.recall.tsv`、`germline.recall.tsv`）只看 truth 键有没有出现在某个 tensor 的 allele 或 filtered 候选里，**和标签规则无关**。

## 2. 第一步：准备 truth allele

### 2.1 拆分（`split_alleles`）

- 只取 chr1–22。每个 ALT 单独成为一个 allele；`<DEL>` 这类符号 allele、`*`、`.`、带 N 的跳过。
- 先去掉 REF/ALT 共同的尾部，再去掉共同的头部（记下起点）。
- 去完以后：
  - REF 和 ALT 等长：每个不同的碱基拆成一个 SNV。MNV `AC>GT` 就变成两个 SNV；
  - REF 为空：INS；ALT 为空：DEL；
  - 都不为空、长度不同：COMPLEX。COMPLEX 不可能等于任何候选，只用来判断"附近有 truth"。
- SNV 和 DEL 的 REF 和 FASTA 对不上的跳过（记在 `reference_mismatch`）。

例子：VCF `chr20:151663 TA>T` → 去掉共同头部 `T` → DEL，pos0 = 151663，REF = `A`。

### 2.2 等价位置（`placements`）

重复序列里的 indel 有很多等价写法（VCF 左对齐；builder 在节点的正向链上左对齐，在反向节点上就是重复的另一端）。代码先把 allele 左移到最左，再一步步右移，把每个位置都列出来（最多 5,000 个）：

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
- DEL 跨过几个连续、同方向的 GRCh38 节点，给跨节点键 `…@n2+n3`（和 builder 的写法一样）；方向不一致的节点上没有键（builder 只在同方向的 mapping 上合并缺失）；
- 路径经过不止一次的节点不给键。

`chr20:151663 TA>T` 在一串 A 里，有 11 个等价位置，给出 11 个键：`30192874:0:DEL:A>` 和 `30192875:0:DEL:A>` … `30192875:9:DEL:A>`。其中 `30192875:0:DEL:A>`（GRCh38 pos0 151664）正好是 PacBio 上一个 tensor 的 `candidate_id`，所以那个 tensor 是"完全匹配"。

### 2.4 truth 表里记下的其它东西

- **span** `lo..hi`：所有等价位置覆盖的区间。INS 的 span 从最左的插入点到最右的插入点。
- **in_bed**：span 是否在这个 truth 自己的 BED 里（只记录，1/2 的判断不看它）。
- **passed**：FILTER 是 PASS 或 `.`。
- **GAP**：germline 另有一类 FILTER 只有 `GAP1`/`GAP2` 的 allele（dipcall 在一条组装单倍型上调用了、另一条没组装出来）。它在正常样本的基因组里（只是合子型不知道），所以也算 germline truth。
- 合在一起就是第 0 节的"合格 truth"：somatic 要 passed；germline 要 passed 或只有 GAP。其它 FILTER 的 germline（如 dipcall `HET1`/`HET2`：组装的单倍型本身有歧义）不是合格 truth。

## 3. 第二步：tensor 的位置

### 3.1 GRCh38 坐标（`ReferencePath.linear`）

- 候选所在节点是 GRCh38 节点，就换算出 `grch38 = {chrom, pos0, ref, alt, node_reverse}`（碱基转回 GRCh38 正链）。
- 跨节点的缺失要求这些节点都是 GRCh38 节点、方向一致、在参考上首尾相接，否则没有坐标。
- 判断时用的区间（`linear_interval`）：SNV/DEL 是 `[pos0, pos0 + len(REF))`；INS 是插入点两侧各 1 bp，即 `[pos0 − 1, pos0 + 1)`。

### 3.2 分支节点：anchor（`anchor`）

Minigraph-Cactus 按拓扑顺序给节点编号，所以分支节点的 ID 夹在它两侧参考节点的 ID 中间。代码据此给分支节点定位：

1. 在节点 ID 下方、上方各找 200 个 ID 以内、最近的 GRCh38 节点（只找到一侧也可以，区间就是那一个节点）；
2. 找到两侧时，两侧必须在同一条染色体上；
3. 两个参考节点**之间**的 GRCh38 碱基（不算这两个节点本身；节点最长 1,024 bp）不能超过 1,024 bp。按 ID 相邻的参考节点有时在 GRCh38 上相距很远，最远上百 Mb，这样的区间里任何 truth 都能匹配上。不算节点本身的长度，是为了两个长节点之间的 SNV 气泡不会因此丢掉位置；
4. 满足以上条件，anchor 区间就是"GRCh38 上靠左的参考节点起点 → 两个节点里最靠右的终点"，写进 labels.ndjson 的 `anchor`。

例子：`31293240:0:INS:>A` 在分支节点上，anchor = `chr20:58872268–58872534`（266 bp）。

找不到或不满足条件的，没有位置（→ −1，见 R5）：
- `30211411:0:INS:>CC`：节点 ID 两侧 200 个以内都没有 GRCh38 节点；
- chr17 上的 `21708846:0:SNP:A>G`（16 条 reads 里 8 条 ALT）：按 ID 最近的两个参考节点 `21708658`、`21709029` 在 `chr17:43,322,822` 和 `chr17:55,501,408`，中间隔 12,178,585 bp，超过 1,024。

PacBio 上没有位置的 tensor（SNV 33,043、INDEL 8,684）绝大多数是第一种（SNV 32,519、INDEL 8,514）；因为间隔超过 1,024 bp 而没有位置的有 SNV 524、INDEL 170 个，间隔从 2,368 bp 到 173,436,931 bp（chr2 的 `29914995:0:SNP:T>A`），中位数 14,336 bp，超过 1 Mb 的 198 个，397 个在 chr9；两侧不在同一条染色体上的没有。

## 4. 第三步：分类规则（按顺序，第一条符合的决定标签）

| 顺序 | 条件 | 标签 | reason |
|---|---|---:|---|
| R1 | AF 低于这类 tensor 的下限，不管 truth：SNV tensor 低于 `--snv-min-af`，INDEL tensor（A1 是 INS 或 DEL）低于 `--indel-min-af`（各自只在给了时） | −1 | `below_snv_min_af` / `below_indel_min_af` |
| R2 | A1 的键等于一个 **PASS 的 somatic truth** | 1 | `representative_allele_is_somatic_truth` |
| R2' | A1 等于 somatic truth，但它们都没 PASS | −1 | `somatic_truth_filtered` |
| R3 | A1 等于一个 **PASS 的 germline truth** | 2 | `representative_allele_is_germline_truth` |
| R3' | A1 等于一个 FILTER 只有 GAP1/GAP2 的 germline truth | 2 | `representative_allele_is_germline_truth_gap_filtered` |
| R3'' | A1 只等于其它 FILTER 的 germline truth（如 HET1/HET2） | 先记下，继续往下 | 见 R7 |
| R4 | A1 没在 R2/R3 定下来，但 **site 里别的 allele 是合格 truth** | A1 和它重叠 >0.45 取 truth 的标签，否则 0 | `allele_partial_somatic_truth` / `allele_partial_germline_truth` / `truth_matches_non_representative_allele` |
| R4' | site 里别的 allele 只是没 PASS 的 somatic truth | −1 | `somatic_truth_filtered` |
| R5 | 没有 GRCh38 坐标，也没有 anchor | −1 | `not_on_unique_grch38_node` |
| R6 | 位置不在 confident 区里（anchor 用区间中点判断） | −1 | `outside_confident_region` |
| R7 | R3'' 的情况（在 confident 区里） | 0 | `germline_truth_filtered` |
| R8 | **同一位置有 PASS somatic truth**，A1 和它重叠 >0.45 | 1 | `residual_partial_somatic_truth` |
| R9 | 分支节点（有 anchor）、R8 没对上 | −1 | `off_reference_no_truth_match` |
| R10 | （只剩有 GRCh38 坐标的）同一位置有合格 germline truth，A1 和它重叠 >0.45 | 2 | `residual_partial_germline_truth` |
| R11 | 10 bp 内有任何 truth allele（somatic 或 germline） | 0 | `near_truth_allele_mismatch` |
| R12 | 其它所有 | 0 | `confident_no_truth_allele` |

说明：
- **R2–R4 的 1 和 2 不看 BED**（它们都在 R6 之前）；R8、R10 的部分匹配在 R6 之后，只在 confident 区里。R4 的 0（`truth_matches_non_representative_allele`）也不看 BED；其它的 0 都在 confident 区里。
- **−1 的来源**：R1 AF 太低、R2'/R4' somatic truth 没 PASS、R5 没位置、R6 BED 外、R9 分支节点没对上 somatic。R1、R5、R6 是测试时不看 truth 也能做的判断（AF、有没有位置、在不在调用用的 BED 里）：tumor-only caller 在测试时会遇到每一个 tensor，所以 −1 只留给测试时同样能去掉的 tensor，BED 内 GRCh38 节点上其它没有 truth 的 tensor——包括真实变异旁边的错误和假象——都是 0。R9 另有原因：germline truth 只有 GRCh38 上的键，分支节点上"没对上 truth"说明不了什么，多数分支节点上根本没有 truth，所以不当负样本。
- **somatic 优先**：A1 同时等于 somatic 和 germline truth 时，R2 先命中，标 1。
- R4 里先试 somatic，再试 germline；同一类有几个 truth 时取重叠最大的那个。
- R8/R10 的"同一位置"：truth 的 span 和 tensor 的区间相交（有 GRCh38 坐标时两边各放宽 1 bp，所以紧挨着也算；anchor 两边各放宽 10 bp）。
- **分支节点上的 tensor 永远拿不到 2**：R3/R4 要 truth 键，分支节点上没有键；R10 在 R9 之后，只剩有 GRCh38 坐标的 tensor。分支节点只能是 1（R8）或 −1。
- R11 的"10 bp 内"用的是所有 truth allele，包括没 PASS 的和 COMPLEX，并且按整个 span 算：INS 的 span 从最左到最右的插入点，因为 reads 可能把它写在重复里的任何位置。
- R1 在最前面，但 labels.ndjson 里照样记下 `somatic`/`germline` 匹配，所以能看出被它去掉的 truth。

### 各规则的真实例子（HG008 PacBio，另外注明的除外）

**R1 AF 太低（HG008 Illumina，`--snv-min-af 0.07`）。** `39113762:0:SNP:C>A` = somatic `chr4:108442081 C>A`（0/1，PASS），82 条 reads 里 5 条 ALT（AF 0.061）→ **−1**。HG008 Illumina 里 A1 是 somatic truth、被 R1 去掉的只有它和 `51692207:1:SNP:T>A`（`chr7:12908284 T>A`，7/103，AF 0.068）两个。germline 也一样：`30223463:10:SNP:A>T` = germline `chr20:1258971 A>T`（1/0），3/47 条（AF 0.064）→ **−1**。

**R2 完全匹配 somatic。** `30217844:21:SNP:G>A` = somatic `chr20:964493 G>A`（0/1，PASS）。67 条 reads 全是 ALT → **1**。

**R3 完全匹配 germline。** `30192018:0:SNP:C>T` = germline `chr20:101582 C>T`（0/1，PASS）→ **2**。INDEL 例子：`30192875:0:DEL:A>` = `chr20:151663 TA>T`（1/1）→ **2**。

**R3' GAP 过滤的 germline。** `30749403:652:DEL:A>` = germline `chr20:26689873 CA>C`，FILTER `GAP2`，GT `1/.`，在 dip.bed 外 → **2**（不看 BED）。

**R4 site 里的另一个 allele 是 truth。**
- `30246963:0:INS:>TT`：site 里有 `+TT`（13 条 reads，A1）、`+T`（10）、`+TTT`（10）。germline truth `chr20:2077646 A>ATTT`（2/1）是 `+TTT`，不是 A1。`+TT` 和 `+TTT` 的重叠 = 2/3 = 0.667 > 0.45 → **2**（`allele_partial_germline_truth`，overlap 0.667）。
- `30200048:135:INS:>T`：A1 `+T`（28 条），germline truth `chr20:385860 A>ATT`（0/1）是 site 里的 `+TT`（5 条）。重叠 = 1/2 = 0.5 > 0.45 → **2**（`allele_partial_germline_truth`，overlap 0.5）。
- `30941026:157:INS:>AA`：A1 是 `+AA`（46 条），somatic truth `chr20:42232813 C>CA`（0/1）是 site 里的 `+A`（26 条）。allele 重叠 = 1/2 = 0.5 > 0.45 → **1**（`allele_partial_somatic_truth`，overlap 0.5；allele 重叠已经够了，不再算单倍型重叠）。
- `31379887:70:DEL:T>`：A1 `−T`（31 条），germline truth `chr20:62992291 CTTT>C`（1/0）是 site 里的 `−TTT`（24 条；site 里还有 `−TTTT` 14 条、`−TT` 12 条）。重叠 = 1/3 = 0.333 → **0**（`truth_matches_non_representative_allele`）。
- `31372027:40:DEL:A>`：A1 `−A`（49 条），germline truth `chr20:62188128 C>CA`（1/1）是 site 里的 `+A`（46 条）。类型不同，allele 重叠为 0；germline 不算单倍型重叠 → **0**（`truth_matches_non_representative_allele`）。

**R5 没有位置。** `30211411:0:INS:>CC` 在分支节点上，两侧 200 个 ID 内没有 GRCh38 节点 → **−1**。chr17 的 `21708846:0:SNP:A>G` 两侧参考节点相隔 12,178,585 bp → **−1**（见第 3.2 节）。

**R6 confident 区外。** `30191551:9:DEL:T>` 在 `chr20:60300`，不在 somatic BED ∩ dip.bed 里 → **−1**。分支节点也一样：`30191605:101:DEL:C>` 的 anchor `chr20:61453–61534` 在区外 → **−1**。

**R7 其它 FILTER 的 germline（COLO829T fiberseq）。** `21701893:995:DEL:T>` = germline `chr17:22914140 CT>C`，FILTER `HET1`，GT `./2`，在 confident 区里，28 条 reads 里 22 条 ALT → **0**。六个数据集里只有这一个 site（COLO829T fiberseq 和 ONT 各一个 tensor）；HG008 三个平台都是 0 个。

**R8 部分匹配同一位置的 somatic truth（残余 edit）。**
- `31293240:0:INS:>A`（miss 分析里的例子）：somatic truth `chr20:58872348 T>TAA`。ALT reads 走一个多一个 A 的分支节点，再带一个 `+A` 残余 edit。这个 tensor 在分支节点上，anchor `chr20:58872268–58872534` 覆盖 truth 位置。A1 row group 有 86 行，均匀取 20 行，其中 19 行正好带着 `+AA`，单倍型重叠 = 0.75（怎么算出来的见第 5.2 节）→ **1**（overlap 0.75）。
- `30602891:0:SNP:T>C`：somatic truth 是 `chr20:19677686 CTTTT>C`（删 4 个 T），但 reads 在旁边写出了一个 SNV（`chr20:19677690`，7 条 A1 reads）。SNV 对 DEL 的 allele 重叠为 0。7 行 A1 的原始 `(d_ref, d_truth)`（第 5.2 节）是 4 行 (5, 1)、(6, 2)、(4, 1)、(5, 2)，都离 truth 单倍型更近；b = 1，每行 0.8–1.0，单倍型重叠 1.0 → **1**（overlap 1.0）。所以 SNV tensor 也可以部分匹配一个 DEL truth。

**R9 分支节点，没对上 somatic truth。** `30194624:0:SNP:T>C` 在分支节点上，anchor `chr20:214823–214878`，在 confident 区里，7 条 reads 全是 ALT（AF 1.0），附近没有能部分匹配的 somatic truth → **−1**。

**R10 部分匹配同一位置的 germline truth。**
- `30209021:0:DEL:TTTTTTT>@30209022`（删 7 个 T，跨两个节点，6 条 reads 全是 ALT），germline truth `chr20:775268 CTTTTTTTTTTT>CTTT`（删 8 个 T，1/2）。共同删掉的碱基 = 7，除以较长的 8 = 0.875 → **2**。
- `30251303:0:INS:>AAATA`（18 条 A1 reads），germline truth `chr20:2284361 A>AAAAT`（插入 `AAAT`，2/1）。这个 truth 的键里有 `30251303:0:INS:>AATA`，但它不是这个 site 的 allele，所以走 R10：重叠 0.8（第 5.1 节）→ **2**。

**R11 附近有 truth 但对不上。** `30191903:0:INS:>T`（`chr20:92969`，55 条 reads 里 14 条 ALT，AF 0.255），10 bp 内有 truth allele，但既不相同也重叠不够 → **0**。同一个碱基上就有 somatic truth、但 A1 reads 离它不比离 GRCh38 近的，也是这里的 0：`30419351:0:DEL:A>`（第 5.2 节）。

**R12 附近什么都没有。** `30191732:122:INS:>A`（`chr20:73326`，56 条 reads 里 5 条 ALT，AF 0.089）→ **0**。

## 5. 重叠怎么算

`MIN_OVERLAP = 0.45`，必须**大于** 0.45。allele 重叠对 somatic 和 germline 都算；somatic truth 在 allele 重叠不超过 0.45 时再算单倍型重叠，取较大的那个。记在 labels.ndjson 的 `overlap`（三位小数）和 `partial_truth`（对上的那个 truth）里。

### 5.1 allele 重叠（`allele_overlap`，somatic 和 germline 都用）

只比较同类型的 INDEL，而且要有 GRCh38 坐标；SNV 对 SNV（ALT 不同）、类型不同（INS 对 DEL、SNV 对 INDEL）或在分支节点上，重叠都是 0。

- **DEL 对 DEL**：两者都删掉的 GRCh38 碱基数 ÷ 较长的那个。truth 的位置取整个 span（任何等价位置都行）。
  - 例子（R10）：删 7 个 T 对删 8 个 T，7/8 = 0.875 ✔。
  - 删 1 个对删 2 个，1/2 = 0.5 ✔。
  - 删 1 个 T 对删 3 个 T，1/3 = 0.333 ✘（`31379887:70:DEL:T>`，标 0）。
- **INS 对 INS**：插入点必须在 truth 的 span 里。truth 的插入序列按插入点在重复里的位置轮转后，和 tensor 的插入序列求最长公共子序列（LCS），÷ 较长的那个。
  - `+TT` 对 `+TTT`：LCS 2 ÷ 3 = 0.667 ✔（`30246963:0:INS:>TT`，标 2）。
  - `+T` 对 `+TT`：1 ÷ 2 = 0.5 ✔（`30200048:135:INS:>T`，标 2）。
  - `+T` 对 `+TTT`：1 ÷ 3 = 0.333 ✘。
  - `+AAATA` 对 germline `chr20:2284361 A>AAAAT`（插入 `AAAT`，插入点轮转后对齐）：LCS 4 ÷ 5 = 0.8 ✔（`30251303:0:INS:>AAATA`，标 2）。

### 5.2 单倍型重叠（`haplotype_overlap`，只对 somatic truth）

用 tensor 里 A1 那组 reads 的实际序列来判断它们带的是不是 truth 的变化，所以"走分支 + 残余 edit""跳过节点 + 错配"这类写法、以及和 truth 类型不同的 edit（如上面 SNV 对 DEL 的 `30602891:0:SNP:T>C`）也能认出来。分支节点上的 tensor 只能靠它匹配 somatic truth。它只在 allele 重叠不超过 0.45 时才算（`truth_overlap`）。

1. **取 reads。** 从 tensor 的 channel 0（每行 read 的碱基，编码 1–5 = A/C/G/T/N，其它值是空位）在 A1 和 REF 的 row group 里各均匀取最多 20 行、10 行（`EVIDENCE_ROWS`），每行去掉空位后至少要有 30 个碱基（`MIN_EVIDENCE_BASES`），不够的丢掉。tensor 宽 101 列，所以每行是一条 read 大约 100 bp 的局部序列。节点在 GRCh38 上是反向的，就先把 reads 取反向互补。
2. **两条单倍型。** GRCh38 上 truth 两侧各 90 bp（`[pos0 − 90, pos0 + len(REF) + 90)`，pos0、REF 是 VCF 的 allele 去掉共同头尾之后的）作 REF 单倍型，在其中换入 truth allele 得到 truth 单倍型。
3. **原始距离。** 每条 read 分别算它整条放进 REF 单倍型、放进 truth 单倍型所需的最少编辑数（read 从头到尾都要对上，单倍型两端可以空出；`fit_distance`），记作原始的 `d_ref`、`d_truth`。分支节点上的 read 不知道方向，正反两个方向都算，用 min(d_ref, d_truth) 较小的那个方向的两个距离。
4. **背景 `b`。** 每条 REF read 取 min(d_ref, d_truth)，排序后取第 ⌊n/2⌋ 个（从 0 数；偶数条时是中间两个里较大的那个），就是这个位置的测序背景噪声；没有 REF read 时 b = 0。
5. **`event`** = truth 的长度：SNV 为 1，INDEL 为 max(len REF, len ALT)。
6. **每条 A1 read 的 overlap。**
   - 原始 `d_truth ≥ d_ref`（离 truth 单倍型不比离 GRCh38 近）：overlap = 0。
   - 否则 `d_ref`、`d_truth` 各减去 `b`（最小为 0），`shared = max(0, (d_ref + event − d_truth) / 2)`，`overlap = shared / max(d_ref, event)`。
7. **tensor 的单倍型重叠** = A1 reads 的 overlap 排序后第 ⌊n/2⌋ 个（同第 4 步）。合格的 A1 行（第 1 步）少于 3 条时没有单倍型重叠，只看 allele 重叠。

第 6 步的第一条是因为：离两条单倍型一样近（或离 GRCh38 更近）的 read——没覆盖到事件、只带测序错误、或者在同一位置带着别的变化——不是这个 truth 的证据；只做减法的话，b 不小于这两个距离时它会正好得到 0.5，超过 0.45。

用 truth `+AA`（event = 2）来算几种 read：

| read 带的变化 | b | 原始 d_ref | 原始 d_truth | 减去 b 后 | shared | overlap |
|---|---:|---:|---:|---|---:|---:|
| 正好 `+AA` | 0 | 2 | 0 | 2, 0 | 2 | **1.0** |
| `+AAA` | 0 | 3 | 1 | 3, 1 | 2 | 0.667 |
| `+A` | 0 | 1 | 1 | — | — | 0（d_truth ≥ d_ref） |
| 旁边一个无关的错配 | 0 | 1 | 3 | — | — | 0（d_truth ≥ d_ref） |
| 正好 `+AA` | 1 | 2 | 0 | 1, 0 | 1.5 | 0.75 |
| `+AA`，另有 1 个测序错误 | 1 | 3 | 1 | 2, 0 | 2 | **1.0** |
| 正好 `+AA` | 2 | 2 | 0 | 0, 0 | 1 | 0.5 |
| 没覆盖到插入点，另有 2 个测序错误 | 2 | 2 | 2 | — | — | 0（d_truth ≥ d_ref） |

单倍型重叠（和 allele 重叠取较大的）超过 0.45，就算部分匹配这个 somatic truth。由此：
- b = 0 时，只带 truth 一部分的 read 要带超过一半才算：INS truth 长 event、read 带其中 k 个碱基时 d_ref = k、d_truth = event − k。所以 truth `+AA`、read `+A` 记 0。同聚物里差一个碱基的 indel 靠的是 allele 重叠（第 5.1 节），那只在有 GRCh38 坐标时有。
- b 也会压低带着完整事件的 read：原始 (2, 0) 的 read，b = 1 时 0.75，b ≥ 2 时 0.5。

真实例子（HG008 PacBio）：
- `31293240:0:INS:>A`（R8，分支节点，somatic truth `chr20:58872348 T>TAA`）：取的 20 行 A1 里 19 行原始 `(d_ref, d_truth)` = (2, 0)，正好带着 `+AA`，1 行 (3, 1)；10 行 REF 的 min 是 9 个 1、1 个 3，b = 1。于是 19 行各 0.75、1 行 1.0，单倍型重叠 0.75 > 0.45 → **1**（`residual_partial_somatic_truth`，overlap 0.75）。
- `30419351:0:DEL:A>`（`chr20:10302796` 删一个 A，72 条 reads 里 10 条 ALT，AF 0.139）：同一个碱基上是 somatic truth `chr20:10302797 A>T`（0/1，PASS）。10 行 A1 的原始 `(d_ref, d_truth)` 都是 (1, 1)：删掉的 A 离 GRCh38 的 A 和离 truth 的 T 一样远，这些 reads 不带这个 SNV。REF 只有 2 行（min 是 1 和 2，b = 2）。每行记 0，单倍型重叠 0；DEL 对 SNV 的 allele 重叠也是 0 → 不是部分匹配，10 bp 内有 truth → **0**（R11 `near_truth_allele_mismatch`）。
- 分支节点上一样：`30887064:0:DEL:G>`（anchor `chr20:38906296–38906302`，somatic truth `chr20:38906297 T>G`，100 条 reads 里 14 条 ALT）的 14 行 A1 里 13 行 (1, 1)、1 行 (2, 2)，b = 1，单倍型重叠 0 → **−1**（R9 `off_reference_no_truth_match`）。

## 6. 输出

每个 `<tensors>/<SNV|INDEL>/` 里（merge 之后的目录，和 `<chrom>_variant_summary.ndjson` 同序）：
- `<chrom>_shard_NNNNN_labels.npy`（int8，和 tensor 一一对应）。
- `<chrom>_labels.ndjson`，每个 tensor 一行：`candidate_id`、`site_id`、`chrom`、`shard_file`、`index_within_shard`、`label`、`label_name`、`reason`、`somatic`/`germline`（键等于 site 里某个 allele 的 truth，含 `matched_candidate_id` 和 `representative`，表示是不是 A1）、`grch38`、`anchor`（分支节点定了位时）、`partial`（`allele` / `residual` / null），部分匹配时再加 `overlap` 和 `partial_truth`。
- `labels.manifest.json`：`format`（`truth-labels`）、`rules_sha256`（`truth_labels.py` 的 SHA-256）、`created`、标签值、`near_bp`、`min_overlap`、`anchor_reach`、`anchor_gap`、`snv_min_af`、`indel_min_af`（没给时是 null）、`tensors`（tensor 总数）、每条染色体和总的标签计数、每个 reason 的计数、`partial` 计数、truth VCF 和 BED 的路径与 SHA-256、truth 统计、confident 区大小、FASTA 和 reference-path 目录。

另外：`--truth-dir` 下的 `<set>.graph.tsv`（每个 truth allele 和它的键，只取决于图和 VCF）；`--recall-dir`（默认 `--tensors`）下的 `<set>.recall.tsv` 和 `truth_recall.json`（每个 truth allele：是某个 tensor 的 A1、是 site 里的其它 allele、COMPLEX、没有唯一节点键、被 filtered（带原因）、或者没有候选）。

## 7. 六个数据集的计数

目录都在 `/scratch/jshen/data/pansoma_v2_tensors/<数据集>/v3_tensors/{SNV,INDEL}/labels.manifest.json`，2026-09-28 按现行规则打的标签（十二个 manifest 的 `rules_sha256` 都是 `197b5d25bbf4…`）。部分匹配的几行（`allele_partial_*`、`residual_partial_*`）和 manifest 的 `partial` 计数相同。两个 Illumina 数据集用了 `--snv-min-af 0.07`，其它四个没有（六个数据集 build 时 SNV 的 AF 阈值都是 0.06）；六个都没有用 `--indel-min-af`。

**HG008**（`Liss_lab_PacBio_Revio_20240125`、`Liss_lab_Northeastern-ONT-UL-20241216`、`Liss_lab_BCM_Illumina-WGS_20240313`；truth：GIAB `HG008-T_somatic_smvar_benchmark_v0.2` + HG008-N dipcall）

| 标签 | reason | PacBio SNV | PacBio INDEL | ONT SNV | ONT INDEL | Illumina SNV | Illumina INDEL |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | representative_allele_is_somatic_truth | 8,234 | 4,098 | 8,275 | 3,673 | 8,277 | 3,562 |
| 1 | allele_partial_somatic_truth | 1 | 93 | 3 | 391 | 2 | 32 |
| 1 | residual_partial_somatic_truth | 601 | 2,649 | 570 | 2,319 | 2,347 | 1,645 |
| 2 | representative_allele_is_germline_truth | 332,598 | 147,786 | 343,284 | 157,914 | 310,349 | 55,800 |
| 2 | representative_allele_is_germline_truth_gap_filtered | 18,229 | 1,544 | 26,891 | 2,177 | 5,089 | 373 |
| 2 | allele_partial_germline_truth | 0 | 4,540 | 0 | 16,407 | 0 | 609 |
| 2 | residual_partial_germline_truth | 0 | 37,554 | 0 | 44,265 | 0 | 13,768 |
| 0 | confident_no_truth_allele | 51,229 | 1,309,897 | 256,169 | 1,681,012 | 2,041,479 | 6,552 |
| 0 | near_truth_allele_mismatch | 35,311 | 205,183 | 57,455 | 275,572 | 594,093 | 24,343 |
| 0 | truth_matches_non_representative_allele | 137 | 3,577 | 160 | 20,053 | 325 | 497 |
| 0 | germline_truth_filtered | 0 | 0 | 0 | 0 | 0 | 0 |
| −1 | outside_confident_region | 652,007 | 98,698 | 1,594,086 | 156,864 | 578,632 | 25,093 |
| −1 | off_reference_no_truth_match | 48,953 | 99,441 | 48,284 | 138,291 | 126,818 | 18,218 |
| −1 | not_on_unique_grch38_node | 33,043 | 8,684 | 37,329 | 12,075 | 20,923 | 973 |
| −1 | below_snv_min_af | 0 | 0 | 0 | 0 | 1,312,554 | 0 |
| −1 | somatic_truth_filtered | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | **合计** | **8,836** | **6,840** | **8,848** | **6,383** | **10,626** | **5,239** |
| 2 | **合计** | **350,827** | **191,424** | **370,175** | **220,763** | **315,438** | **70,550** |
| 0 | **合计** | **86,677** | **1,518,657** | **313,784** | **1,976,637** | **2,635,897** | **31,392** |
| −1 | **合计** | **734,003** | **206,823** | **1,679,699** | **307,230** | **2,038,927** | **44,284** |
| | 总数 | 1,180,343 | 1,923,744 | 2,372,506 | 2,511,013 | 5,000,888 | 151,465 |

**COLO829T**（`COLO829T_Illumina`、`COLO829T_fiberseq`、`COLO829T_ONT`；truth：`COLO829T_somatic_snv_indel.vcf.gz` + SMaHT BED，germline 是 COLO829BL dipcall）

| 标签 | reason | Illumina SNV | Illumina INDEL | fiberseq SNV | fiberseq INDEL | ONT SNV | ONT INDEL |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | representative_allele_is_somatic_truth | 38,239 | 1,096 | 38,550 | 1,217 | 38,355 | 1,121 |
| 1 | allele_partial_somatic_truth | 1 | 5 | 0 | 9 | 1 | 88 |
| 1 | residual_partial_somatic_truth | 487 | 248 | 246 | 413 | 291 | 546 |
| 2 | representative_allele_is_germline_truth | 331,542 | 59,271 | 370,509 | 158,478 | 397,907 | 182,304 |
| 2 | representative_allele_is_germline_truth_gap_filtered | 5,954 | 544 | 23,521 | 2,062 | 34,055 | 2,839 |
| 2 | allele_partial_germline_truth | 0 | 573 | 0 | 4,781 | 0 | 24,440 |
| 2 | residual_partial_germline_truth | 0 | 15,100 | 0 | 42,833 | 0 | 44,274 |
| 0 | confident_no_truth_allele | 902,500 | 8,606 | 90,005 | 1,061,525 | 715,458 | 3,360,091 |
| 0 | near_truth_allele_mismatch | 332,182 | 25,496 | 51,419 | 245,149 | 118,872 | 439,356 |
| 0 | truth_matches_non_representative_allele | 190 | 493 | 246 | 4,007 | 241 | 23,512 |
| 0 | germline_truth_filtered | 0 | 0 | 0 | 1 | 0 | 1 |
| −1 | outside_confident_region | 336,146 | 15,587 | 1,024,172 | 98,542 | 1,706,426 | 144,439 |
| −1 | off_reference_no_truth_match | 66,201 | 18,442 | 78,703 | 148,906 | 101,246 | 249,781 |
| −1 | not_on_unique_grch38_node | 11,280 | 791 | 38,458 | 14,634 | 51,452 | 20,945 |
| −1 | below_snv_min_af | 480,878 | 0 | 0 | 0 | 0 | 0 |
| −1 | somatic_truth_filtered | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | **合计** | **38,727** | **1,349** | **38,796** | **1,639** | **38,647** | **1,755** |
| 2 | **合计** | **337,496** | **75,488** | **394,030** | **208,154** | **431,962** | **253,857** |
| 0 | **合计** | **1,234,872** | **34,595** | **141,670** | **1,310,682** | **834,571** | **3,822,960** |
| −1 | **合计** | **894,505** | **34,820** | **1,141,333** | **262,082** | **1,859,124** | **415,165** |
| | 总数 | 2,505,600 | 146,252 | 1,715,829 | 1,782,557 | 3,164,304 | 4,493,737 |

- `somatic_truth_filtered` 在六个数据集里都是 0：没有 tensor 的 allele 只匹配到没 PASS 的 somatic truth（HG008 的 somatic truth 全部 PASS）。
- `germline_truth_filtered` 只有 COLO829T fiberseq 和 ONT 各 1 个，是同一个 site（R7 的例子）。
- SNV 的 `allele_partial_germline_truth`、`residual_partial_germline_truth` 都是 0：germline 只算 allele 重叠，而 SNV 对 SNV 只有完全相同才算（第 5.1 节）。SNV 的部分匹配只有 somatic，靠单倍型重叠。
- `below_snv_min_af` 去掉的 SNV 里，A1 是 somatic truth 的：HG008 Illumina 2 个，COLO829T Illumina 389 个；A1 是 germline truth 的：1,825 和 1,660 个。

## 8. 读代码时注意到、值得讨论的点

以下都是按现在代码的实际行为整理的；没有另外注明的数字来自 HG008 PacBio。

### 8.1 还没定的点

1. **分支节点上的 tensor 拿不到 2，没对上 somatic 的也不当 0。** 分支节点只能是 1（R8）或 −1（R5、R6、R9；短读长的 SNV 还有 R1）。PacBio 上定了 anchor 的 tensor：SNV 72,989 个、INDEL 116,695 个；其中标 1 的 SNV 322、INDEL 1,353，R9 的 −1 有 SNV 48,953、INDEL 99,441（AF ≥ 0.3 的 SNV 33,925、INDEL 41,741），其余（SNV 23,714、INDEL 15,901）在 confident 区外。AF 高的那些里应该有不少是人群分支上的真实 germline allele（例如 `30194624:0:SNP:T>C`，7/7 条 reads 是 ALT）：训练时它们既不是 2 也不是 0，而 tumor-only caller 在测试时会遇到它们。
2. **truth 在 site 的另一个 allele 上、A1 和它重叠不够 → 0，BED 外也是 0。** PacBio INDEL 有 3,577 个，其中 3,223 个是同一位置 INS 对 DEL（A1 是缺失、truth 是插入，或反过来；类型不同，allele 重叠是 0；这里面 3,195 个只有 germline truth，germline 又不算单倍型重叠），这里面 2,731 个是同一个碱基的删一个对插一个（如 A1 `−A`、truth `+A`）；2,804 个的 truth allele reads 不少于 A1 的一半，300 个不少于 A1。例子 `31372027:40:DEL:A>`：A1 `−A` 49 条，germline truth（1/1）`+A` 46 条 → 0。site 里确实有真实 germline allele 的 reads，却当成负样本。
3. **SNV 对 SNV 只看完全相同。** 同一位置 ALT 不同（例如 truth G>A、A1 G>C）的 allele 重叠为 0；somatic 还有单倍型重叠兜底，germline 没有。PacBio SNV 的 137 个 `truth_matches_non_representative_allele` 全是这种（136 个 germline，1 个 somatic）。
4. **germline 不做单倍型重叠。** 所以"走分支 + 残余 edit""跳过节点 + 错配"这类写法的 germline 变异，只能靠 allele 重叠（R4、R10，且只在 GRCh38 节点上）：在 GRCh38 节点上对不上的落到 R11/R12 的 0，在分支节点上的是 R9 的 −1。
5. **somatic 优先；但 A1 是 germline 时，site 里别的 somatic allele 不管用。** A1 同时等于 somatic 和 germline truth 时 R2 先命中，标 1：HG008 有 209 个 somatic allele 和某个 germline allele 完全相同（多为正常样本杂合，可能是肿瘤 LOH 或两个 benchmark 冲突），PacBio 上这样的 tensor 有 61 个（INDEL 51、SNV 10）。反过来，A1 等于 germline truth、site 里另一个 allele 是 somatic truth 时 R3 先命中，整个 tensor 标 2：PacBio 254、ONT 234、Illumina 40 个。
6. **somatic 单倍型的残余 edit，重叠不够的标 0。** miss 分析（`analysis/hg008_somatic_miss_20260926`）列出了未命中的 somatic INDEL（I1–I3 类：allele 全部在图里、部分在图里留下残余 edit、GRCh38 路径上写成别的 edit）的 ALT reads 上的残余 edits（至少 3 条 ALT reads 带它的）。按不同候选计，它们所在 tensor 现在的标签：

   | | PacBio | ONT | Illumina |
   |---|---:|---:|---:|
   | 1 | 3,322 | 2,668 | 1,709 |
   | 2 | 430 | 239 | 323 |
   | 0 | 273 | 324 | 1,684 |
   | −1 | 406 | 424 | 729 |
   | 没有 tensor | 1,971 | 679 | 1,762 |

   0 大多是 R11 `near_truth_allele_mismatch`（231 / 237 / 1,679），其余是 R4 的 `truth_matches_non_representative_allele`（42 / 87 / 5）。这些 tensor 带的是 somatic 单倍型的 reads，却是负样本；Illumina 上标 0 的和标 1 的差不多一样多。
7. **残余 edit 正好等于（或部分匹配）一个 germline truth → 2。** 例如 germline `TA>T`、somatic `TAA>T`，残余 edit `DEL A` 等于 germline 键，R3（和 R4）在 R8 之前 → 2。上表里的 2：PacBio 430 个（A1 正好是 germline truth 300 个、R4 allele 部分匹配 108、R10 residual 部分匹配 21、GAP 1），ONT 239，Illumina 323。像上面这个例子，reads 其实是 somatic 单倍型的残余，却标 2；另一部分是 ALT reads 上附带的真实 germline 变异，这些标 2 是对的。
8. **1/2 不看 BED，大部分 0 要看。** confident 区外的 truth 匹配（R2–R4）照样是 1/2（R4 的 0 也不看 BED），而同样位置的非 truth tensor 是 −1（R8、R10 的部分匹配也只在区内）。分支节点判断是否在 confident 区只看 anchor 区间的中点。
9. **单倍型重叠只用 tensor 窗口里的序列。** 每条 read 大约 100 bp、最多 20 条，取中位数；A1 少于 3 条时不算。长插入（>50 bp）或跨出窗口的事件算不准。

### 8.2 现行规则已经定下的点

- **重叠阈值是 0.45（`MIN_OVERLAP`）。** 这样同聚物里差一个碱基的短 indel（`+T` 对 `+TT`、删 1 个对删 2 个，allele 重叠 0.5）和较长的（`+TT` 对 `+TTT`，0.667）一样算同一个事件；差两个碱基的（`+T` 对 `+TTT`、删 1 个对删 3 个，0.333）不算。
- **单倍型重叠只认离 truth 比离 GRCh38 更近的 read。** 原始 `d_truth ≥ d_ref` 的 A1 read 记 0，减去背景 b 不改变这一点（第 5.2 节）：离两条单倍型一样近的 read 不是这个 truth 的证据。所以同一位置有 somatic truth、A1 reads 却不带它的 tensor 不是部分匹配，例如 `30419351:0:DEL:A>`（R11 的 0）、`30887064:0:DEL:G>`（R9 的 −1）。
- **分支节点上没对上 somatic truth 的 tensor 是 −1（R9），不是 0。** germline truth 只有 GRCh38 上的键，分支节点上"没有 truth"不代表没有变异，多数分支节点上根本没有 truth；部分匹配 somatic 的 1 照样保留（PacBio SNV 322、INDEL 1,353）。
- **短读长数据集的 SNV：AF < 0.07 → −1（R1，`--snv-min-af 0.07`，只用于 HG008 Illumina 和 COLO829T Illumina）。** 效果等于用 SNV AF 阈值 0.07 建这两个数据集（build 时是 0.06），但不用重建 tensor；AF 是测试时不看 truth 也能用的过滤条件。HG008 Illumina SNV 有 1,312,554 个（26.2%）、COLO829T Illumina 480,878 个（19.2%）因此是 −1，其中 A1 是 somatic truth 的只有 2 个和 389 个。用 `orchestrate` 跑短读长数据时，在 prepare 加 `--label-snv-min-af 0.07`，finalize 就按这个下限打标签。INDEL 的下限（`--indel-min-af`，prepare 里是 `--label-indel-min-af`）和 SNV 的分开设，六个数据集都没用。
