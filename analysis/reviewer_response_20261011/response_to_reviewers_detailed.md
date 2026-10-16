# Pansoma: point-by-point response to the NBT reviews (NBT-RA72852A) — 详细版

2026-10-12 草稿。简版 `response_to_reviewers.md`（每条 < 350 词）仍是回复信的底稿。本文件对每条意见多给三样东西：我们的判断、
更完整的 response、以及 response 需要的实验或功能的设计。**Pipeline 代码本身没有改**；下面凡是要写代码的，都只是设计。

怎么读：

- **Comment**：审稿意见里的每一个具体要求，逐条列出，回复信里一条都不能漏。
- **评估（内部）**：意见是否成立、我们怎么答、风险在哪。
- **Response**：回复信英文正文，比简版长，最后定稿时再压缩。
- **现有 pipeline 能否解决**：能 / 部分 / 不能，以及已有的证据和文件。
- **实验 / 功能设计**：每个实验有编号（E1–E25，总表在第 6 节），写目的、数据、步骤、输出、判定规则和代价。标注
  "新脚本" 的只需要分析脚本；标注 "改 pipeline" 的要改 `indexed_gam_pipeline_v4` 或 `pansoma_net_v2`，这次都不动。
- `[待补: …]` 是等实验的数字；`[待核实: …]` 是要核实的事实。两者发出前都必须处理掉。
- 数字默认是 chr1（test chromosome）、benchmark BED 内。
  - **SNV 的 call 层面数字**来自 2026-10-01 的六个 SNV 模型，这些模型训练用的是 10-03 重建之前的 tensor，必须重训（E1）。
  - **INDEL 的 call 层面数字**来自 10-10 的 `small_nofilt` runs（lr 4e-4）；10-11 定的 recipe（lr 1e-4）还没有跑（E2）。
  - **candidate 层面数字**（truth allele 有没有 tensor）来自现在的 tensor 集，可以直接用。

---

## 1. 总体策略（内部）

### 1.1 回复的主线

三位审稿人的核心问题是同一个：**性能提升有多少来自 pangenome 本身**（R1.1–R1.3、R3.1、R3.2、R3.6）。其次是：图会吞掉 somatic
（R1.4）、HapMap 不合适（R1.5、R2.1、R3.10）、recall 低（R1.6、R2.3、R2.4）、只在单个肿瘤内训练测试（R1.7、R3.9）、
graph 坐标的 calls 没进指标（R3.3、R3.4）。

建议的主线（R1 自己的说法就是 "the pangenome effectively serving as a population-level generic normal prior"）：

> The pangenome acts as a population prior applied at alignment. On identical reads and identical code it removes about nine
> in ten matched-normal germline candidates before classification, at a cost of 0.3–4 points of somatic SNV retention. The
> classifier also reads the per-node haplotype count. We quantify what this gains end to end (reference-level ablation), what
> it costs (graph-absorbed somatic alleles, which population filtering removes anyway), and how far it generalises (six tumor
> lines, three platforms, held-out lines).

各项主张的现状：

| 主张 | 现状 |
|---|---|
| 图在分类前去掉约九成 germline candidate，somatic SNV 保留只少 0.3–4 个点 | **已有**（六套 × 三个参考臂，candidate 层面，见 R1.2 表） |
| 分类器在用 haplotype count（ch6） | 旧 build 上已有 counterfactual，新 build 上要重做（E1.1） |
| 同一代码下 HPRC 比 GRCh38-only 好多少（end to end） | **没有**，决定性实验（E3） |
| 跨肿瘤泛化 | **没有**（E12） |
| graph 坐标 calls 是真的 | **没有**，而且现在多数可能是 FP（E9） |
| INDEL | 弱，只能诚实写限制（R1.6、R3.10） |

### 1.2 风险和应对

- **风险 a：E3 在 PoN 之后打平。** 之前去掉 ch6 的 nopc ablation 在 PoN 之后就是打平的；linear 臂也可能这样。
  - 准备：PoN 前后都报。
  - 准备：报 candidate 数和算力（图上 tensor 少 43–69 %）。
  - 准备：按区域分层（E10），看困难区域里有没有差别。
  - 准备：报完整输出的指标（E9e）。
  - 如果真打平，主张改成 "the graph moves the population prior upstream and reaches the same accuracy with fewer candidates"，
    不再说 graph 带来精度优势。这一点要在 E3 出结果后和用户一起定。
- **风险 b：unplaced（graph 坐标）calls 大多是 FP。**
  - PoN 后它们占 SNV PASS 的 8–51 %；全算 FP 时，HG008T precision 只有 0.28–0.46。
  - INDEL 也有：HG008T Illumina INDEL raw PASS 里有 534 个 unplaced。
  - 投稿前必须先做 E9（anchor 投影、graph 原生 population filter、组装验证）。
- **风险 c：H1437 PacBio 的异常。**
  - 它的 somatic SNV 只有 70.9 % 有 tensor，其他数据集的 PacBio 是 87–95 %（d46：HCC1395 87.1 %、COLO829T 89.8 %、HG008T 94.9 %）。
  - 抽查的漏检位点 VAF ≈ 0.34、深度 68，不是低 VAF；原因还没查清。
  - 查清之前，不要把 H1437 放进联合训练或测试集。
- **风险 d：INDEL 在长读长上很差。**
  - chr1 PoN 后 F1：HG008T Illumina 0.239、PacBio 0.132、ONT 0.015；COLO829T Illumina 0.116、fiberseq 0.055、ONT 0.004。
  - INDEL 只作为次要结果，不主张优势。
  - 审稿人指出 ClairS-TO 在部分 ONT INDEL 设置上更好，要明说（E8 后核实）。
- **风险 e：三个参考臂的 vg 版本不一致。**
  - GRCh38-only 臂用 vg 1.66（vg 1.77 读不了 1.66 建的 linear minimizer index）。
  - v2.1 臂用 vg 1.77。
  - v1.1 臂的 GAM 是 2025 年比对的 [待核实: 版本]。
  - 简版 R3.1 原来写的 "same vg giraffe version" 不准确，已改（2026-10-12）。
  - 应对：E3 里加一个版本对照（见 R3.1），或者把 linear 臂用 vg 1.77 重建索引后重比对。

### 1.3 两个需要先定的设计决定

1. **主图用哪个。** 建议在 E3 之后把主图换成 HPRC v2.1 d46，理由有三：
   - 四个新细胞系都在 d46 上。
   - v2.1 有 EUR 单倍型，而 HG008 donor 是欧洲人。
   - 六套数据上，d46 的 somatic 保留都 ≥ d9，germline 保留都 < d9（R1.2 表）。

   代价是 COLO829T 的 P25/P50/P75 mixture 现在只有 d9 版，要用 d46 重新比对 9 套。另一个办法是 mixture 继续用 d9，作为单独的实验。
2. **主评估用哪种。** 建议主结果用跨样本：联合训练其他细胞系，在留出细胞系的 chr1 上测试（E12）。同一样本内（chr2–22 训练、
   chr1 测试）的模型只作为上界和 ablation。原稿是 COLO829T 训练、HG008T 测试，R3.9 说的 "independent HG008T benchmark"
   就是这个；现在 12 个模型全是同一样本内训练测试，退了一步。

### 1.4 顺序和依赖

```
E1 (SNV 重训, d9) ─┬─> E3 (三臂 + nopc) ─┬─> E4b / E19 (FP 分解)
E2 (INDEL recipe) ─┘                     ├─> E10 (区域分层)
                                          └─> E16 (seed, bootstrap)
E9a–e (unplaced 投影 / 过滤 / 组装验证 / 完整输出指标) ── 在 E1 之后就能做，结果进所有表
E8 (比较器同协议重跑) ── 和 E1 并行；E5c、E10、E16 都要用
CASTLE finalize + H1437 查清 ─> E12 (跨样本 / leave-one-line-out) ─> E16
E14a (mixture 异常) ─> E14b/c (mixture 模型, HCC1395 dilution)
文字项 (R2.m, R3.m1, R3.m3, R3.m4, R3.8) ── 随时
```

---

## 2. Summary of the main changes (for the letter)

We rebuilt the Pansoma pipeline and revised the manuscript around the questions the reviewers raised.

1. **What does the pangenome contribute?** A reference-level ablation processes the same reads with identical code on a
   GRCh38-only graph, HPRC v1.1 and HPRC v2.1, with and without the haplotype-count channel (R1.1, R3.1).
2. **What does the graph remove, and at what cost?** For somatic truth, matched-normal germline and background candidates we
   report how many each reference keeps (R3.2). We traced every somatic allele hidden by a graph path to the patient's own
   assemblies and to population databases (R1.4).
3. **Is the complete output correct?** Graph-coordinate calls are now projected to GRCh38 where the graph allows it, filtered
   against population haplotypes, validated against assemblies, and included in the headline metrics (R3.3, R3.4).
4. **Does the model generalise?** Four breast and lung cancer lines were added. The headline results come from models tested
   on cell lines held out from training, with repeated training and block-bootstrap intervals (R1.7, R3.9).
5. **Low VAF.** The HapMap mixture was removed from the main benchmark and replaced by tumor/normal dilutions with true somatic
   mutations (R1.5, R2.4, R3.10).

The claims are restricted to tumor-only calling in cancer samples, with SNV and INDEL performance reported separately (R3.8,
R3.10).

（内部：第 2 条已完成；第 1、3、4、5 条要等 E3、E9、E12、E14 做完才能这样写。）

---

## 3. Reviewer 1

### R1.1 "Pangenome-native" is not supported

**Comment.**

- The pangenome is used mainly during vg Giraffe alignment; the downstream representation is a local pileup.
- Neither the input nor the network models graph edges, bubbles/superbubbles, cross-node context or long-range haplotype
  structure, so the method resembles DeepVariant, DeepSomatic or ClairS.
- Replacing linear coordinates with node IDs is not enough; haplotype and path information is discarded.
- Request: substantial additional evidence, **or a re-design**, showing measurable benefit from pangenome-specific
  information beyond alignment.

**评估（内部）。** 部分成立。

- 已经不成立的部分：第一版 tensor 截断在节点边界，也不带 path 信息。现在的 tensor 每行沿着 read 自己的 walk 跨节点，并带两个
  graph 专有通道：ch5 是 read 当前所在节点的碱基，ch6 是经过该节点的 HPRC haplotype 数。旧 build 上的 counterfactual 说明
  模型确实在用 ch6。
- 仍然成立的部分：
  - 去掉 ch6 的 nopc ablation 在 PoN 之后打平，说明 ch6 的信息和 PoN 大量重叠。
  - 不建模 bubble 和 competing branches 是事实。
- 策略：
  - 删掉 "fully exploit"，改称 pangenome-aware。
  - 列出 graph 信息进入的三个环节，每个都给数字。
  - 用 E3（三个参考臂 + 无 ch6 臂）把贡献拆开。
  - 对 "re-design" 的要求：可选做 E24（locus / bubble-aware tensor）。不做的话，在文中明确承认网络不建模拓扑。

**Response.**

We agree that "fully exploit the richness of pangenome data" overstated the design, and we removed it. We also agree that
graph node IDs alone would not make a method pangenome-native. In the revised pipeline, graph information enters at three
stages after alignment, and we now quantify each.

1. *Candidate definition.* Candidates are edits relative to the nodes each read traverses, so alleles carried by population
   haplotypes in the graph produce no candidate. On identical reads, 76–82% of matched-normal germline alleles produce a
   candidate on a GRCh38-only graph built with the same aligner, against 7–13% on HPRC v1.1 and v2.1 (six datasets, three
   platforms).
2. *Path-following tensors.* Each tensor row follows its read's own walk through the graph across node boundaries. Two
   channels are graph-specific: the base of the node the read traverses (channel 5) and the number of HPRC haplotypes through
   that node (channel 6). Rows on a rare branch are therefore marked as such.
3. *The classifier uses channel 6.* Lowering the site haplotype count of true somatic tensors to the typical germline value
   lowers their median somatic probability from 0.89 to 0.10; raising every count to the maximum more than triples germline
   false calls (2,129 to 7,505). [待补: repeated on the final models, six datasets]

To separate these contributions we added a reference-level ablation (R3.1): the same reads, aligned and processed with the
same code on a GRCh38-only graph, HPRC v1.1 and HPRC v2.1, plus an HPRC arm without channel 6. The step from GRCh38-only to
HPRC-without-channel-6 measures alignment and candidate definition; the step to full HPRC measures the haplotype channel.
[待补: chr1 precision, recall and F1 per arm, before and after population filtering, with intervals]

Calls are written in graph coordinates and projected to GRCh38 where a unique projection exists, or through the bubble
anchors of the graph (R3.4).

Pansoma does not model bubbles or long-range haplotypes explicitly, and the revised text says so. We describe it as
pangenome-aware: graph-native in its coordinates, candidates and read representation.
[可选, 若做 E24: We also tested a locus mode in which each tensor includes the reads on the competing branches of the
enclosing bubble together with each branch's haplotype frequency. [待补: result]]

**现有 pipeline 能否解决：部分。**

- 已有：
  - 跨节点的 per-read walk tensor、ch5 / ch6（v4 README §8）。
  - 三个参考臂的 candidate 表（R1.2，现有 tensor 集的 `truth_recall.json`）。
  - 旧 build 上的 counterfactual（`pansoma_net_v2_runs/pathcount_counterfactual_20260930`）和 nopc ablation（`analysis/v1_vs_v2_20260929`）。
- 缺：
  - 新 build 上的 counterfactual（E1.1，要先有 E1）。
  - 三臂模型和无 ch6 臂（E3）。
  - bubble 建模（E24，改 pipeline，可选）。

**实验 / 功能设计。**

- **E1.1 最终模型上的 counterfactual**（新脚本；复用 `pathcount_counterfactual_20260930` 的脚本）。
  1. 对 E1 的六个 SNV 模型，取 chr1 的全部 tensor。
  2. 跑三种改动并推理：
     - somatic truth tensor 的 site 列（±10）ch6 改成 germline tensor 的 site ch6 中位数；
     - germline tensor 的 ch6 改成 H；
     - 在 tensor 之间打乱 ch6。
  3. 输出 p_somatic 分布、PASS 数、PoN 前后 F1 的变化（表 + 一张图）。
  4. 代价：每套 GPU 推理 < 1 h。
- **E24 locus / bubble-aware tensor**（改 pipeline，可选，P2）。
  1. snarl 分解：d9 和 d46 的 `.snarls` 都已有（`AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.snarls`、
     `AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.snarls`）。为每个 target node 找包含它的最小 snarl（≤ 1 kb）。
  2. build 时，除了经过 target node 的 reads，再取经过同一 snarl 其他 branch 的 reads。它们在同一 GAM 区域，batch 的 fetch
     基本已经读到。
  3. 新加一个 plane：每行标出它走的是哪条 branch（target / GRCh38 branch / 其他），以及该 branch 的 haplotype frequency。
     AF 改为在整个 snarl 的 reads 上计算。
  4. 先在 HG008T Illumina、PacBio 和 COLO829T Illumina 上建 chr1 + 训练染色体，和标准 tensor 比较 chr1 F1（PoN 前后，E16 的 CI）。
  5. 判定：PoN 后 F1 提升超出 CI 才进主流程；否则作为补充结果写 "explicit bubble modelling did not add accuracy"，这本身也
     回答了审稿人。
  6. 代价：builder 改动 + goldens 更新 + 3 次 build + 6 个模型，约 2–3 周。

### R1.2 Composition and size of the pangenome

**Comment.**

- No systematic analysis of how the number and diversity of haplotypes, haplotype walks and alternative paths affect node
  pileups and germline identification.
- The graph decides which alleles are represented, and so which variants are suppressed before the classifier.
- With HPRC2 (many more haplotypes): better population representation, but possibly more mapping ambiguity and more somatic
  mutations that overlap population alleles absorbed.

**评估（内部）。** 成立，而且已经有很好的数据。

- d9 和 d46 的 candidate 层面对比现在六套都有了：d46 在每一套上 somatic 保留都略高或持平、germline 保留都略低。所以 "larger
  graph absorbs more somatic" 在同样的频率过滤下不成立。
- 但两张图都按约 10 % 的 haplotype 频率过滤（d9 = 9/90，d46 = 46/464）。决定吸收的主要是频率阈值，而不是 haplotype 数。
  "size" 这条轴被过滤吸收掉了，要说清楚。
- 建议加一个便宜的 in-silico 分析（E7）：用 HPRC 每条单倍型的 allele，同时扫 haplotype 数和频率阈值两条轴，再用 d9 / d46
  的实际吸收做校准。
- mapping ambiguity 用 E11（MAPQ、模拟 reads）回答，见 R1.3 / R3.6。

**Response.**

We added HPRC v2.1 as a third reference. All six datasets were realigned and rebuilt with identical code on a GRCh38-only
graph, HPRC v1.1 (90 haplotypes; nodes kept when on at least 9) and HPRC v2.1 (464 haplotypes, including European
individuals absent from v1.1; nodes kept when on at least 46). The table gives the fraction of truth alleles that reach a
candidate tensor:

| Dataset | Somatic SNV (GRCh38 / v1.1 / v2.1) | Somatic insertions | Somatic deletions | Matched-normal germline |
|---|---|---|---|---|
| HG008T Illumina | 0.984 / 0.952 / 0.956 | 0.744 / 0.358 / 0.396 | 0.856 / 0.504 / 0.525 | 0.757 / 0.077 / 0.072 |
| HG008T PacBio | 0.983 / 0.947 / 0.949 | 0.796 / 0.396 / 0.424 | 0.884 / 0.593 / 0.601 | 0.779 / 0.106 / 0.094 |
| HG008T ONT | 0.986 / 0.951 / 0.954 | 0.529 / 0.356 / 0.361 | 0.629 / 0.526 / 0.524 | 0.775 / 0.116 / 0.107 |
| COLO829T Illumina | 0.909 / 0.901 / 0.901 | 0.898 / 0.502 / 0.544 | 0.928 / 0.631 / 0.656 | 0.791 / 0.082 / 0.076 |
| COLO829T PacBio (Fiber-seq) | 0.903 / 0.898 / 0.898 | 0.962 / 0.535 / 0.578 | 0.965 / 0.718 / 0.729 | 0.822 / 0.117 / 0.108 |
| COLO829T ONT | 0.896 / 0.893 / 0.892 | 0.736 / 0.465 / 0.479 | 0.769 / 0.679 / 0.675 | 0.816 / 0.131 / 0.124 |

On every dataset, v2.1 kept as many or more somatic alleles than v1.1 and fewer germline alleles. At the frequency filter
both graphs use (about 10% of haplotypes), the five-fold larger haplotype set improved germline suppression slightly and
did not increase somatic absorption. [待补: trained models on each arm (R3.1)]

What decides absorption is which alleles the graph contains, which depends on the frequency filter more than on the
haplotype count. We quantified both axes directly. For random subsets of 8–88 HPRC v1.1 haplotypes [and up to 462 v2.1
haplotypes] and frequency thresholds from a single haplotype to 10%, we computed the fraction of matched-normal germline
alleles and of somatic truth alleles the graph would contain, and checked the prediction against the absorption observed on
v1.1 and v2.1. [待补: curves; expected: germline coverage saturates early, while somatic absorption grows mainly as the
frequency threshold is lowered, concentrated in homopolymer and STR indels]

Absorption itself is characterised in R1.4, and mapping ambiguity with more haplotypes in R1.3 and R3.6 (MAPQ distributions
on all three references; simulated reads with known origin).

**现有 pipeline 能否解决：部分。**

- 已有：上表（2026-10-12 从六套 linear / d9 / d46 的 `tensors/truth_recall.json` 算出，status = `tensor_representative`
  的比例）。HCC1395 三套和 H1437 PacBio 只有 d46。
- 缺：三臂模型（E3）、E7、MAPQ（E11a）。

**实验 / 功能设计。**

- **E7 in-silico 图组成分析**（新脚本，不改 pipeline，P1）。
  1. 输入：
     - HPRC v1.1 每条单倍型的 GT：`/scratch/jshen/data/hprc-v1.1-mc-grch38.raw.vcf.gz`（未过滤图的 deconstruct VCF），以及
       `hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz`（已 wave 归一，优先用它）。
     - v2.1 只有 d46 图，没有全图 VCF，要下载 v2.1 的 vcfbub / wave VCF（用户决定）。
  2. 把每个位点拆成 per-haplotype 的 allele 集合（bcftools 在这里坏了，用 pysam）。INDEL 在重复区按 equivalent placement
     归一（复用 labels 的 placement 枚举）。
  3. 查询集合：
     - germline：HG008-N dipcall、COLO829BL dipcall；
     - somatic：HG008T、COLO829T、HCC1395 truth（exact allele）。
  4. 对每个 k ∈ {8, 16, 32, 64, 88}（v2.1 再加 128、256、462），每个 k 抽 20 次，可以按 superpopulation 分层；频率阈值
     f ∈ {≥ 1 条, 2 %, 5 %, 10 %}。allele 算 "represented"，当子集里携带它的单倍型数 ≥ max(1, f·k)。
  5. 输出：germline represented 比例和 somatic represented 比例随 k、f 变化的曲线；按 SNV、HP ≥ 7 INDEL、STR INDEL、
     其他 INDEL 分层。
  6. 校准：k = 88、f = 10 % 的预测对比 d9 的实测（germline no_candidate 比例；perfect-bypass 2,478 / 497）；k = 462、f = 10 %
     对比 d46。
  7. 注意：represented 不等于 absorbed（还取决于比对），所以只作为趋势和上界，用实测校准。
  8. 代价：CPU，几小时；VCF 解析要几十 GB 内存。
- **E7b（可选，P2）一个实测点。** 从 v2.1 全图用 `vg clip -d 23`（约 5 %）做一张新图，建 dist / min / zipcodes，只比对
  HG008T PacBio 的 chr1 reads（从 d46 GAM 取经过 chr1 节点的 read 名，回 FASTQ 取 reads），建 chr1 tensor，看吸收和 germline
  抑制是否落在 E7 的曲线上。需要下载 v2.1 全图 GBZ。

### R1.3 Difficult regions and ONT ambiguity

**Comment.**

- Pangenome alignment can stay ambiguous in repetitive, paralogous, structurally complex and low-mappability regions,
  especially for noisy ONT reads; more haplotypes may add competing mappings.
- Aggregate statistics (more perfect reads, fewer mismatches/indels/soft clips) do not show that the improvement happens
  where it should.
- ONT has more graph-coordinate predictions: improved resolution or extra ambiguity? Region-specific benchmarking and error
  analysis are needed.

**评估（内部）。** 成立。

- 现在只有 node 层面的 branch 统计（`analysis/branch_node_alignments_20261002`）、COLO829T 的 SMaHT 区域，以及漏检分析里的
  S1 类（没有 MAPQ > 10 read）。
- 需要：
  - HG008T 的 GIAB stratification（服务器上没有，要下载）；
  - 三个参考臂 × 区域的 alignment 指标和 calling 指标；
  - ONT unplaced calls 的处理（E9）。
- 预期：ONT 的 unplaced calls 多数是 FP（branch 上的 rare germline 和测序错误），要诚实承认 "ambiguity, now filtered and
  validated"，不能说成 resolution。

**Response.**

We added region-stratified analyses for alignment, calling and graph-coordinate calls.

1. *Alignment.* For the same reads on the GRCh38-only graph, HPRC v1.1 and v2.1 (and the GRCh38 BAM used by the comparators),
   we report per-read MAPQ, mismatches, indels and soft clipping per GIAB v3.5 stratification for HG008T and per SMaHT easy,
   difficult and extreme region for COLO829T. [待补] On HPRC v1.1, 4.9–8.7% of node visits by MAPQ > 5 alignments lie on
   non-reference nodes, and 94.3–98.1% of those visits are edit-free. Within each node-length bin, branch and GRCh38 nodes
   are about equally often perfect, so the gain is not an artefact of short branch nodes.
2. *Correct versus perfect alignment.* Reads simulated from the HG008-N assembly, which is in none of the graphs, measure
   misplacement and MAPQ calibration per stratum on each reference (R3.6). [待补]
3. *Calling.* Precision and recall per stratum, for each reference arm and each comparator, with block-bootstrap intervals.
   [待补] ONT reads reach segmental duplications that short reads miss: among missed HG008T truth SNVs, one has no
   MAPQ > 10 read on ONT, against 48–55 on the other platforms.
4. *ONT graph-coordinate calls.* We tested resolution against ambiguity directly (R3.3, R3.4).
   - The three platforms share one graph, so a call can be compared across platforms in node coordinates.
   - Supporting reads are checked against the tumor and matched-normal assemblies.
   - Each call's ALT is checked against all HPRC haplotypes, including alleles removed by the frequency filter.
   - Each call is placed in a stratum through the GRCh38 anchors of its bubble.

   [待补: fraction reproduced across platforms, confirmed by assembly, and found in population haplotypes; strata] We report
   the outcome as it is: [预期] most ONT graph-coordinate calls before filtering were population alleles on rare branches or
   ONT errors, not newly resolved somatic variants.

**现有 pipeline 能否解决：部分。**

- 已有：
  - d9 上的 branch-node 统计（`analysis/branch_node_alignments_20261002`）；
  - COLO829T 的 SMaHT easy / difficult / extreme 区域；
  - 漏检分析（`analysis/hg008_somatic_miss_20260926`、`analysis/tensor_recall_20260930`）。
- 缺：
  - GIAB stratification（要下载）；
  - 三个参考臂的 per-read 指标（E11a）；
  - 分层评估（E10）；
  - ONT unplaced 的分析（E9）。

**实验 / 功能设计。**

- **E10 区域分层评估**（新脚本 + 下载，P1）。
  1. 下载 GIAB v3.5 GRCh38 stratifications（一个 tarball，几 GB；用户决定）。选这些 stratum：
     - easy（`notinalldifficultregions`）、all difficult；
     - low mappability、segmental duplications；
     - tandem repeats 和 homopolymer（≥ 7）；
     - MHC。

     COLO829T 用已有的 SMaHT easy / difficult / extreme。
  2. calling：对每个参考臂 × caller × 数据集的 PASS calls，在 rtg 输出的 tp / fp / fn VCF 上按 stratum BED 计数（比每层重跑
     vcfeval 快）。PoN 前后都算。
  3. alignment：把 E11a 的 per-read 统计按 read 的 GRCh38 位置分层（graph 臂用 `vg surject` 的位置，或 linear BAM 的位置）。
  4. unplaced calls：用 E9a 找到的 bubble anchor，给每个 unplaced call 一个近似 GRCh38 区间，统计它们的 stratum 分布。预期
     集中在 segdup 和 TR。
  5. 输出：分层 P / R / F1，加 bootstrap CI（E16），以及每层 HPRC 和 linear 的差值。
  6. 代价：CPU，主要是计数脚本。
- E9（unplaced calls 的整套处理）见 R3.3 / R3.4；E11（alignment 正确性）见 R3.6。

### R1.4 Candidate generation removes somatic variants already encoded in the graph

**Comment.**

- Candidates come only from mismatches/indels against the traversed node. If a somatic allele is already an alternative
  path and tumor reads align perfectly to it, no candidate, tensor or call is produced.
- Called a "fatal methodological issue".
- Supplementary Fig. 3 wrongly describes these false negatives as GRCh38 reference-bias artefacts.
- An allele on an alternative path is not necessarily germline: the same allele can be germline in one person and somatic in
  another.

**评估（内部）。** 机制和概念上完全成立，必须接受，并撤回 reference-bias 的说法。但定量上，图特有的损失很小：

- 被吸收的 truth 里，84–87 % 是其他 HPRC 个体的 allele（群体 allele），只有 4.4–6.8 % 是病人自己的 germline。
- 用现在的 PoN，被吸收的 truth 也有 94.7–99 % 会被 PoN 标掉（`tmp/graph_absorbed_pon_20261012`）。所以任何用 population
  filter 的 tumor-only 流程也会丢掉它们。
- 图特有的损失：HG008T 2 个 SNV + 121 个 INDEL（INDEL truth 的 1.3 %），COLO829T 1 个 SNV + 13 个 INDEL（0.6 %）。
- 吸收集中在 homopolymer / STR INDEL。HCC1395 的 SEQC2 truth 里长 homopolymer INDEL 很少，d46 上 truth insertion 有 tensor
  的比例是 0.71–0.80，HG008T 只有 0.36–0.42。

策略：

- 接受机制，给出定量，说明图只是把同一个 population prior 提前到比对阶段，而且不可逆。
- 可选：给不做 population filter 的用户提供 branch candidate 模式。
- 建议先做 E6（不改 pipeline 的可行性分析），看找回多少、代价多大，再决定要不要做 E25。

**Response.**

The reviewer is correct. An allele present on a graph path is not necessarily germline, and reads that align without edit
to that path produce no candidate. We withdrew the reference-bias interpretation of Supplementary Fig. 3 and quantified the
loss.

Using the matched-normal assemblies (HG008-N v6.2; COLO829BL), the HG008-T tumor assembly and the HPRC haplotype walks, we
classified every truth allele whose supporting reads align without edit to an existing path of HPRC v1.1:

- **HG008T:** 2,478 alleles (201 SNVs, 2,277 INDELs). 110 (4.4%) are the patient's germline allele. 2,085 (84.1%) are absent
  from both normal haplotypes but carried by other HPRC individuals (median haplotype frequency 0.17).
- **COLO829T:** 497 alleles. 34 (6.8%) are germline and 430 (86.5%) are alleles of other individuals.

Most are one-unit homopolymer or STR changes, and absorption rises with tract length.

These alleles are population alleles, and tumor-only pipelines remove them at the filtering stage instead. With the
population filtering we apply to all callers (gnomAD and CoLoRSdb allele matches; R2.6), 99% of the absorbed SNVs and 95–97%
of the absorbed INDELs are removed, against 4–10% of all truth SNVs and 52–62% of all truth INDELs. The loss specific to the
graph is therefore 2 SNVs and 121 INDELs for HG008T (1.3% of its truth INDELs) and 1 SNV and 13 INDELs for COLO829T (0.6%).
The graph applies the same population prior as the filter, but at alignment and irreversibly. We now present this as a
limitation of population-based tumor-only calling that the graph makes explicit.

Absorption depends on the repeat content of the truth set. On HCC1395, whose SEQC2 truth contains few long-homopolymer
INDELs, 71–80% of truth insertions reach a candidate on HPRC v2.1, against 36–42% for HG008T.

[待补: comparator recall at the absorbed loci, raw and after the same filtering (E8)]
[可选, 若做 E25: An optional branch-candidate mode treats reads on rare graph branches as candidates; it recovers [待补] of
the absorbed alleles at a cost of [待补] additional candidates.]

**现有 pipeline 能否解决：部分（分析完成，功能没有）。**

- 已有：
  - `analysis/graph_absorbed_somatic_20261001`、`analysis/graph_absorbed_somatic_colo829t_20261005`（分类、read path、组装）；
  - `tmp/graph_absorbed_pon_20261012`（和现行 PoN 的交叉表，2026-10-12）；
  - HCC1395 的 d46 truth_recall。
- 缺：
  - d46 和 HCC1395 的吸收分类（脚本可复用，E5b）；
  - 比较器在这些位点上的表现（E5c，要先有 E8）；
  - branch candidate（E6 / E25）。

**实验 / 功能设计。**

- **E5b d46 / HCC1395 的吸收分类**（新脚本运行，P1）：复用 graph_absorbed 的 s0–s14（HCC1395 没有组装，只做 "HPRC 其他个体
  allele vs 无" 的分类和 PoN 交叉表）。
- **E5c 比较器在吸收位点上**（新脚本，P1，依赖 E8）：在 2,478 / 497 个位点上数 ClairS-TO / DeepSomatic-TO 的 PASS calls，
  分三档：raw、各自默认过滤后、同一 PoN 后。DeepSomatic-TO 在 COLO829T fiberseq chr1 上把 55 % 的 truth INDEL 标成
  germline，这是现成的一个数据点。
- **E6 branch candidate 可行性**（新脚本，不改 pipeline，P1）。
  1. 用 d9 `.snarls` 和 graph index 的 path count，算每个 snarl 每条 branch 的 haplotype frequency。
  2. 扫 chr1 的 sorted GAM（GAI fetch），对每个 snarl 统计每条 branch 上无编辑的 reads 数和 AF（分母是整个 snarl 的 reads）。
  3. branch candidate 规则：branch frequency ≤ f（f = 5 / 10 / 20 %），≥ 3 reads，AF ≥ 0.08，且 branch 不是 GRCh38 branch。
  4. 计数：
     - 规则捕获多少 absorbed truth（HG008T chr1 BED 内 17 SNV + 169 INDEL）；
     - chr1 上一共产生多少 branch candidate（预期多数是病人的 germline）；
     - PoN 后还剩多少 truth。
  5. 判定：PoN 后能多找回的 truth ≤ 20 个，而新增候选上万，就不做 E25，把这个 trade-off 写进文中。
  6. 代价：1–2 天（python 扫 chr1 区域）。
- **E25 branch candidate 模式**（改 pipeline，可选，P2，取决于 E6）。
  1. builder 在 snarl 层面产生 candidate：allele = branch 相对 GRCh38 branch 的差异。
  2. tensor 以 snarl 为窗口，GRCh38 branch 的 reads 一起放进来（和 E24 共用代码）。
  3. 标签：truth placement 和 branch 序列比对。
  4. 模型用同一 recipe。
  5. 报 PoN 前后 recall 的提升和 FP 代价。

### R1.5 HapMap mixture as a somatic benchmark, and data leakage

**Comment.**

- The HapMap "somatic-like" positives are germline alleles of the minor individuals. Pansoma's rationale is that population
  germline alleles in the graph are suppressed, so the two assumptions conflict; the benchmark may systematically
  disadvantage a pangenome caller and is not a biologically appropriate test of acquired low-VAF mutations.
- Leakage: HPRC v1.1 contains four of the six donors (HG00438, HG02257, HG02486, HG02622). The companion paper (Fu, Xin, Miao
  et al., bioRxiv 2026.01.04.697580) raises this, yet the donors were not excluded.

**评估（内部）。** 两点都成立。

- 建议把 HapMap 从主结果里去掉，用含真实获得性突变的 tumor/normal dilution 替代：
  - COLO829T 三个平台的 P25 / P50 / P75 已经建好（d9）；
  - HCC1395 Illumina 用本地数据就能做（E14c）。
- 不建议花大成本做去掉 donor 的图（E23）。除非用户想在 Supplement 里保留 HapMap 的低 VAF 结果，那就要做，否则泄漏问题无解。
- v2.1 也包含这四个 donor（R3.m3），换成 d46 并不能解决泄漏。
- 必须先查 COLO829T fiberseq mixture 的异常：P75 / P50 比 P25 还差，per-tensor AP 只有 0.31 / 0.33，P25 是 0.84。不查清楚，
  dilution 结果不能拿来替代 HapMap。

**Response.**

We agree with both points. Four of the six HapMap donors (HG00438, HG02257, HG02486, HG02622) are in HPRC v1.1 and also in
v2.1, and the positives of the mixture are germline alleles that the graph is designed to absorb, as our companion study
also noted. We removed the HapMap mixture from the main benchmark and withdrew the model trained on it.

Low-VAF performance is now evaluated on tumor/normal dilutions that contain true acquired mutations and were processed with
the same pipeline:

- COLO829T diluted with COLO829BL to nominal purities of 25, 50 and 75% on all three platforms
  [待核实: measured purity per set];
- HCC1395 diluted with HCC1395BL (Illumina, nominal purities 5–50%), scored against the SEQC2 truth. This tumor can be held
  out from training, so the dilutions also test generalisation.

We report recall per truth-VAF bin and per purity (R2.4). [待补]

[二选一]
(a) The HapMap analysis is no longer part of the manuscript.
(b) The Supplement keeps a HapMap analysis on an HPRC graph from which the four donors were removed, and states the limitation
that its positives are germline alleles. [待补]

**现有 pipeline 能否解决：部分。**

- 已有：
  - COLO829T P25 / P50 / P75 九套 tensor（d9）。
  - fiberseq 联合模型的结果：PoN 后 F1 在 P100 / P75 / P50 / P25 分别为 0.785 / 0.699 / 0.663 / 0.722，P75 / P50 的异常没查。
  - SEQC2 数据在本地（`/scratch/jshen/data/SEQC2/`）：tumor WGS_FD_T_1/2/3、normal SRR7890827（HCC1395BL，2×150，约 58x）、
    一个湿实验 titration run `Tumor_Purity_Titration/SPP_GT_1-1_1.sra`（约 148x，purity 档位不明）。
  - 2025-10 建过去掉一个样本的图 `non_HG005_VG_Indexes/mgi_non_HG005.d9.*`，流程可以照搬（E23）。
- 缺：异常排查（E14a）、Illumina / ONT mixture 模型（E14b）、HCC1395 dilution（E14c）。

**实验 / 功能设计。**

- **E14a mixture 异常排查**（新脚本，P1，先做）。
  1. 每个 purity，取 somatic truth tensor 的 A1 AF 分布，和期望（truth VAF × 实测 purity）比较。truth VAF 用 SMaHT 的
     VAF_Ill / VAF_PB。
  2. 比较每套的深度和 read 数（`discovery_report.json`、`node_stats.json`），确认 mixture 的 GAM 来源对（Illumina 旧 mix GAM 曾经
     把 R1 用了两次）。
  3. 核对各套 `labels.manifest.json`：truth、BED、AF floor、规则 SHA 一致。
  4. 按 AF bin 算每个 purity 的 per-tensor AP。看 P75 / P50 的 somatic 是不是被判成 germline：这两档 somatic AF 接近
     germline het 的 0.3–0.5，而 P25 远离这个区间。
  5. 判定：数据问题就重建；模型问题（AF 分布和 germline 重叠）就把 purity 作为训练多样性（E14b），并在文中报告。
- **E14b Illumina / ONT mixture 模型**（不改代码，P1）：照 `jobs/submit_colo829t_fiberseq_mix.sh` 的做法，P100 + P75 + P50 +
  P25 联合训练，chr1 按 purity 分别测试，SNV 和 INDEL 都做。
- **E14c HCC1395 in-silico dilution**（数据都在本地，不用下载，P1；用户决定）。
  1. tumor = WGS_FD_T_1 + _2（约 108x，和现有 HCC1395T_Illumina_d46 同一批 reads），normal = HCC1395BL `SRR7890827`（约 58x）。
  2. 总深度约 60x，混成 P5 / P10 / P20 / P50（P5 = 3x tumor + 57x normal，normal 够用）。用 `seqtk sample` 按 pair 抽样，
     固定 seed，tumor 和 normal 的 read 名不冲突。
  3. normal reads 只用来稀释，不作为 matched normal 输入，和 COLO829T mixture 的做法一样。
  4. giraffe d46（短读长每个 chunk 16 线程，参照 `tmp/d46_illumina_chunks_20261003`）→ gamsort → discovery → prepare（和
     HCC1395T_Illumina_d46 同样的选项：rule (a)、germline QUAL ≥ 10）→ finalize。
  5. 用 HCC1395 同一样本内的模型，以及 E12 里留出 HCC1395 的跨样本模型，各自预测 → vcf → vcfeval，按 purity 和 truth VAF
     bin 报 recall / precision。
  6. 代价：每个 purity 一次 Illumina giraffe（d46 短读长约 12 h / chunk × 8 chunks）+ sort（11–37 h）+ build；4 档约一周
     集群时间，磁盘每档约 1–2 TB 中间文件（做完删 FASTQ 和未排序 GAM）。
- **E14d SEQC2 湿实验 titration**（可选，P2）：先查 `SPP_GT_1-1` 的 SRA metadata，确定 purity（要上网查 NCBI，用户同意后做），
  再按 E14c 第 4–5 步处理。
- **E23 去掉 donor 的 HapMap 图**（可选，P2，只在保留 HapMap 时做）。
  1. 按 `mgi_non_HG005.d9` 的流程，从 v1.1 全图去掉 4 个 donor（`vg gbwt` 去样本，或用 MC 重建时去掉）。
  2. 重新 `vg clip -d 9`，建 dist / min / zipcodes 和 graph index。
  3. 重新比对 HapMap mixture，建 tensor 和标签。
  4. 和原 d9 上的结果对比。
  5. 代价：建图和索引约 1–2 天，比对和 build 约 3–5 天。

### R1.6 Low sensitivity

**Comment.**

- F1 gains come from precision; Supplementary Fig. 1 and Tables 2–3 show lower SNV recall than ClairS-TO and DeepSomatic-TO on
  COLO829T and HG008T.
- INDEL recall is ~4–22% (COLO829T) and ~13–19% (HG008T).
- In tumor-only calling no matched normal or second stage can rescue missed variants, so sensitivity matters; F1 alone does
  not compensate.

**评估（内部）。** 成立。

- 现在的规则：PASS 阈值在 validation 上按 recall 0.9 选（偏 recall），PoN 后 HG008T SNV R 0.83–0.84（旧模型）。
- 但比较器还没在同一协议下重跑（E8）。没有 E8，就不能说 recall 和它们持平。
- INDEL 是真问题，分三段：
  - candidate 阶段：d9 上 truth insertion 有 tensor 的只有 36–54 %。
  - PoN 阶段：现行 PoN 把 HG008T chr1 truth INDEL 的 413 / 635 标掉，任何用它的 caller PoN 后 recall 都 ≤ 0.35。
  - 我们的曲线上限：PoN 后 0.257。
- 回答策略：
  1. 用一张 "漏检漏斗" 定量 recall 在哪里丢（candidate → 分类器 → population filter）。
  2. PASS 规则偏 sensitivity；所有记录都带 FILTER 标签保留下来，没有不可逆的删除。
  3. 和比较器比 matched precision 下的 recall，并给完整 PR 曲线。
  4. INDEL 限制诚实说明。
  5. 可选的改进：低 AF build（E15）、INDEL error panel（E17b）。

**Response.**

We agree that sensitivity is essential for tumor-only calling, and we changed both the operating point and the reporting.

*Operating point.* The PASS threshold is now chosen on validation blocks to reach 90% validation recall, and population
filtering raises precision afterwards. No call is removed irreversibly after candidate generation: records below the
threshold (down to p = 0.01) are kept as LowQual and population matches are tagged, not deleted, so users can move along the
precision–recall curve. We report full curves for every caller, recall at matched precision, and a high-sensitivity
operating point (validation recall 0.95). [待补: final models; comparators under the identical protocol (R3.m2)]

*Where recall is lost.* For every truth allele we now report the stage at which it is lost (Supplementary Table [待补]):

| Stage | SNV (HG008T, chr1) | INDEL (HG008T, chr1) |
|---|---|---|
| Truth allele reaches a candidate tensor | 95% | 36–59% |
| Called (PASS, before population filtering) | 86–87% | [待补] |
| Kept after population filtering | 83–84% | 23–25% |

[待补: final numbers; the same funnel for COLO829T and the held-out lines; comparators' PASS and post-filter recall]

The SNV misses at the candidate stage are mostly graph SNP alleles and repeat absorption on HG008T (R1.4), and low-VAF
subclonal mutations on COLO829T (R2.4).

*INDELs.* The INDEL limitation is real, and we state it. Two causes lie outside the classifier. First, only 36–54% of truth
insertions reach a candidate on HPRC v1.1, against 53–96% on the GRCh38-only graph; the difference is absorption into
population paths of the graph (R1.4). Second, 413 of the 635
HG008T chr1 truth INDELs (65%) are exact population alleles at the frequencies used for filtering, so after population
filtering no tumor-only caller can exceed a recall of 0.35 on this chromosome. The claims of the manuscript are now
restricted accordingly (R3.10).

**现有 pipeline 能否解决：部分。**

- 已有：
  - R0.9 规则；graph_vcf 保留 LowQual（p ≥ 0.01），PoN 只打标签（`filter_panel_of_normals.py` 默认保留记录）；
  - vcfeval 输出完整曲线和 chr1 best point；
  - 漏检分解（`analysis/tensor_recall_20260930`）；
  - PoN 对 truth 的上限（每个 vcfeval report 都有）。
- 缺：E1 / E2 的最终模型、E8 比较器、漏斗脚本（E1.6）。

**实验 / 功能设计。**

- **E1.6 recall 漏斗**（新脚本，P0，随 E1 / E2 / E8 一起出）。
  1. 每个最终模型、每个 truth allele 取四个状态：
     - candidate 状态：`truth_recall.json` / per-truth 表的 status；
     - 是否 PASS：rtg 在 raw PASS 上的 tp-baseline；
     - 是否被 PoN 标记：truth 的 PoN tags；
     - 最终结果：PoN 后 PASS 的 tp-baseline。
  2. 按 SNV / INS / DEL、按 truth VAF bin、按重复区类别分层。
  3. 比较器：truth → PASS → PoN（它们看不到 candidate 阶段）。
  4. 输出：Supplementary Table + 一张 Sankey 或堆叠图。
- E8（比较器）见 R3.m2；E15（低 AF build）见 R2.4；E17（INDEL 分层、error panel、PoN 阈值重选）见 R3.10 / R2.6；E16（CI）见 R3.9。

### R1.7 Training scale and tumor diversity

**Comment.**

- A ~200 M-parameter ConvNeXt trained mainly on COLO829 (40–50 k truth variants): is this enough data?
- Report the exact parameter count and the numbers of positive and negative training examples per platform and variant type,
  and show that the data volume supports the capacity.
- COLO829 is one UV-dominated melanoma; the model may learn tumor-specific features. DeepSomatic trained across several lines
  with extensive augmentation.
- Extend training to tumors with distinct mutational processes and show directly that the model is not tumor-specific.

**评估（内部）。** 成立。

- 现在 SNV 模型 50.4 M、INDEL 8.7 M，已经比 199.6 M 小，而且有 overfit 的证据（199.6 M 在 lr 1e-4 下 validation loss
  0.51 → 2.16）。
- 但所有模型都是同一样本内训练测试。审稿人要的是跨肿瘤证据。
- 需要做：
  1. 参数量和每类样本数的表（从 checkpoint 导出，E21）。
  2. learning curve 和模型大小对比（E13），说明数据量撑得住这个容量。
  3. 6 个细胞系的 leave-one-line-out（E12）。
  4. 按突变签名分层的 recall（UV、APOBEC、吸烟等）。
- 先决条件：CASTLE 三株全部 finalize 并做完 QUAL10 relabel；H1437 的异常查清。
- augmentation：purity dilution 已有（COLO829T）。coverage downsampling 可以在训练时随机丢行实现，属于训练代码的改动
  （pansoma_net_v2），这次只设计。

**Response.**

*Model size.* The SNV network now has 50.4 M parameters (ConvNeXt, depths 3-3-27-3, widths 96–768) and the INDEL network
8.7 M (depths 2-2-6-2, widths 64–512). The earlier 199.6 M model overfitted on our data (validation loss rising from 0.51 to
2.16 by epoch 12) and gave no better test accuracy. Learning curves on 10–100% of the training blocks, and the three model
sizes, show where accuracy saturates. [待补: E13]

*Example counts.* Supplementary Table [待补] lists somatic, germline and non-variant training examples per dataset, platform
and variant type. For example, HG008T SNV sets hold 8.8–10.3 k somatic, 0.31–0.37 M germline and 0.07–1.9 M non-variant
tensors on chr1–22; COLO829T SNV sets hold about 38 k somatic tensors.

*Tumor diversity.* We added two breast cancer lines (HCC1395, HCC1954) and two lung adenocarcinoma lines (H1437, H2009), each
with published multi-technology truth sets (SEQC2, CASTLE) and sequenced on Illumina, PacBio HiFi and ONT. Together with
melanoma (COLO829T, UV signature) and pancreatic cancer (HG008T), the six lines cover distinct mutational processes [待补:
dominant signatures per line] and add about 323 k somatic truth alleles, against about 46 k for COLO829T.

*Tumor specificity, tested directly.*

1. Leave-one-line-out: for each platform, a model trained on five lines is tested on the sixth (chr1 and, for lines with
   few truth variants, all chromosomes). [待补]
2. A model trained on COLO829T alone is tested on the other lines and compared with the leave-one-line-out model on the
   same test line. [待补]
3. Recall is stratified by mutational context, for example UV-type C>T at dipyrimidines and APOBEC-type TCW mutations, to
   show whether a model trained without a signature still detects it. [待补]

*Augmentation.* Purity dilutions (25–75%) are part of training (R2.4). [待补, 若做: coverage downsampling during training]

**现有 pipeline 能否解决：部分。**

- 已有：
  - 模型大小和 overfit 证据（`pansoma_net_v2/README`，SNV recipe）；
  - 标签计数（v4 README §5，每套的 `labels.manifest.json`）；
  - d46 上 HCC1395 三套、H1437 PacBio 已 finalize；HCC1954 / H2009 PacBio 在 finalize，ONT / Illumina 在 sort / build
    （2026-10-12）。
  - 训练代码支持多套 `--tensors` 联合训练（fiberseq mixture 用过）。
- 缺：E12、E13、E21；H1437 排查。

**实验 / 功能设计。**

- **E12 跨样本 / leave-one-line-out**（不改 pipeline；P0，主结果）。
  1. 全部用 d46：6 个细胞系 × 3 个平台（HG008T、COLO829T、HCC1395、HCC1954、H1437、H2009）。
  2. 主模型：每个平台 × SNV / INDEL，train = 五个细胞系的 chr2–22（validation = 5 % node blocks），test = 留出细胞系的
     chr1。truth 少的细胞系另报全基因组。
     - 共 6 folds × 3 平台 × 2 种 = 36 个 run，先做 SNV 的 18 个。
     - 也可以先只做两个最关键的 fold（留出 HG008T、留出 COLO829T），共 6 个 SNV run。
  3. 对照：
     - 同一样本内的模型（E3 的 d46 臂）作为上界；
     - 只用 COLO829T 训练的模型测其他细胞系，量化 tumor specificity。
  4. 签名分层：truth SNV 按 SBS96 context 分组，比较两类模型的 recall：
     - UV：C>T at TCC / CCC，以及 CC>TT；
     - APOBEC：TCW > T / G；
     - 吸烟：C>A。
  5. 样本均衡：CASTLE 的 somatic 多（H2009 168 k，HCC1954 27 k），要按细胞系加权采样，避免一个细胞系主导。这要在训练代码里
     加 per-set 权重（pansoma_net_v2 小改，不是 v4 pipeline）。
  6. 比较器在 CASTLE 细胞系上的比较，只在 DeepSomatic 训练时留出的染色体上做 [待核实: 哪几条]。
  7. 代价：d46 联合训练的数据量约为单套的 5 倍。一个 SNV run 估计在一块 H100 上要 1–2 天；18 个 SNV run 要 3–5 周单卡时间，
     可以多卡并行。
- **E13 learning curve 和容量**（不改 pipeline，P1）。
  1. 数据集：HG008T PacBio SNV（d9 或 d46，和主图一致）。
  2. 训练块按 node block 抽 10 / 25 / 50 / 100 %，validation 和 test 不变。
  3. 模型 8.7 M / 50.4 M / 199.6 M，各在 100 % 数据上训练。
  4. 输出：chr1 PoN 后 F1 和 per-tensor AP 随数据比例、模型大小的变化。
  5. 代价：7 个 run。
- **E21 训练样本计数表** 见 R1.8。

### R1.8 Construction of the training set

**Comment.**

- Positives: truth variants on GRCh38 are projected to graph coordinates and labelled by exact allele match. How is a
  positive assigned when the truth allele is already a graph path, when the projection is ambiguous, or when graph and linear
  representations differ (especially INDELs)?
- Graph-only positions are assigned GRCh38 coordinates indirectly through adjacent 1-bp SNV nodes.
- Absorbed sites never become positives (Supp. Fig. 3).
- Describe labelling and projection explicitly, including whether normalisation is applied on the graph side, and justify
  why the positive set is not systematically distorted.
- Negatives: what makes a candidate "determinable"? Are germline variants explicitly included, how are artefacts defined,
  and how were VAF and purity sampled?

**评估（内部）。** 成立，旧稿确实写得不清楚。

- 现在的规则是完整的（`analysis/truth_labels.md`，v4 tensor_postprocessing README），可以逐条回答。
- 新增的主要是表，加上一个 "projection 没有扭曲正样本" 的定量：
  - truth allele 的状态分布（representative / non-representative allele / filtered / no candidate），按种类分；
  - no candidate 的原因（已有）。
- "通过相邻 1-bp SNV 节点给 graph-only 位置 GRCh38 坐标" 是旧版做法。现在 off-reference tensor 没有 truth 匹配时记 −1，
  不参与训练。
- 可以提 E9f（用病人组装给 off-reference tensor 打标签）作为改进，它同时帮 R3.3 / R3.4。

**Response.**

We now describe labelling in full (Methods and Supplementary Note); the rule set is versioned by checksum in every label
manifest.

*Positives.*

- Each truth allele is expanded to all equivalent placements in repeats (up to 5,000).
- Each placement is converted to node coordinates along the GRCh38 path of the graph. An insertion at a node junction gets a
  key on both nodes; a deletion over consecutive nodes gets a multi-node key.
- Candidates are left-normalised on the node's forward strand by the same convention, so both sides are normalised.
- A tensor is somatic when its representative allele equals a truth key.
- Partial matches (allele or read-haplotype overlap > 0.45) are recorded separately. SNV partial matches are excluded from
  training because they taught low-AF artefacts; excluding them raised chr1 F1 from 0.623 to 0.665.
- Off-reference tensors (on nodes outside the GRCh38 path) are no longer given coordinates through neighbouring nodes. Without
  an exact truth match they are excluded from training (label "ignore").
  [可选, 若做 E9f: they are labelled from the patient's tumor and normal assemblies instead.]
- Alleles absorbed by graph paths cannot become positives; R1.4 quantifies them.

*Is the positive set distorted?* For every truth allele we report whether it reaches a tensor as the representative allele,
appears as a non-representative allele, is removed by a candidate filter, or has no candidate, with the reason (Supplementary
Table [待补]). On HG008T, 95% of truth SNVs reach a tensor on HPRC v1.1 against 98% on GRCh38-only, and the shortfall is
explained by graph SNP alleles and repeat absorption (R1.4). INDELs are affected more, and we report them separately.

*Negatives.* Labels have three classes.

- *Germline* comes from assembly-based matched-normal truth (HG008-N and COLO829BL dipcall; DeepVariant normal calls with
  QUAL ≥ 10 on the matched normal for the SEQC2 and CASTLE lines), giving 0.31–0.42 M germline SNV tensors per dataset.
- *Non-variant* marks a tensor on a GRCh38 node inside the confident region (somatic ∩ germline benchmark BED) without any
  truth allele. These are sequencing and alignment artefacts, including errors next to true variants.
- "Determinable" means exactly this: inside the confident region, with a GRCh38 position, and either matching a truth allele
  or not. All other tensors are ignored in training, because the germline truth cannot confirm them: those outside the
  region, those without a GRCh38 position, and branch-node tensors without a truth match.
- For truth sets built from short-read call ensembles (SEQC2, CASTLE), INDEL tensors in long homopolymers or STRs without a
  truth allele are also ignored, because such truth sets omit most of these INDELs.

*VAF and purity.* Training uses each line's natural VAF distribution, plus purity dilutions (25–75%) for COLO829T (R2.4).
[待补: VAF histogram of training positives per dataset]

**现有 pipeline 能否解决：能（已有，要整理成 Methods 和表）。**

- 规则：`analysis/truth_labels.md`、`indexed_gam_pipeline_v4/tensor_postprocessing/README.md`；标签计数在 v4 README §5。
- 检查：`analysis/label_checks_20260928`（部分匹配的证据、overlap 0.5 的问题和修正）。
- 缺：每个 run 的训练集计数表（E21）；truth 状态表（E1.6 的一部分）。

**实验 / 功能设计。**

- **E21 训练样本计数表**（新脚本，P1）。对每个最终 run，按 checkpoint 的 args 重放 block split 和 `--ignore-reasons`，统计
  train / validation 中每类的数量（per set、per kind）。另附训练正样本的 VAF 直方图。
- **E9f 用组装给 off-reference tensor 打标签**（改 pipeline 的 labels，可选，P2）。
  1. HG008T：对每个 off-reference tensor 的 A1，取 ALT reads 的 consensus（±50 bp），比对到 HG008-T v3.2 两个 haplotype 和
     HG008-N v6.2 两个 haplotype（复用 `pansoma_net_v2_runs/offref_negatives_20260929` 的脚本）。
  2. 规则：
     - 只在 tumor 组装里出现 → 1；
     - normal 组装里出现 → 2；
     - 两边都没有、而且 ALT reads 带测序错误特征 → 0；
     - 其他 → −1。
  3. COLO829T 只有 normal 组装（COLO829BL dipcall haps），只能标 2；其余仍为 −1。
  4. 用新标签训练，去掉 ch6 的 off-ref rescue，用 E9d 评估 off-reference calls。
  5. 代价：标签脚本 + 重新训练 HG008T 三套。

### R1.m1 Additional benchmarks

**Comment.** Benchmark on the SEQC2 HCC1395 lines and the five additional cell lines of the DeepSomatic study.

**评估（内部）。**

- 已加 HCC1395（SEQC2 truth）和 CASTLE 的三株（HCC1954、H1437、H2009）。HCC1937 没做：truth 已在
  `CASTLE/benchmarks/` 里，reads 要下载（用户决定）。
- 关键的公平性问题：
  - DeepSomatic 用这五株（包括 HCC1395）训练过。
  - HCC1954 / H1437 / H2009 的 truth 就是 DeepSomatic 的训练 truth。
  - 所以在这些细胞系上和 DeepSomatic 比，只能在它训练时留出的染色体上做 [待核实]。
  - HCC1395 的 SEQC2 truth 和 DeepSomatic 无关，但 DeepSomatic 见过 HCC1395 的 reads。
- SEQC2 / CASTLE 的 INDEL truth 来自短读长 caller，长 homopolymer / STR INDEL 很少（SEQC2 HP ≥ 7 只占 15.5 %，HG008T 是
  74.9 %），所以用了 rule (a)。这些细胞系的 INDEL 数字不能直接和 HG008T 比。

**Response.**

We added HCC1395 (SEQC2 truth v1.2.1) and three CASTLE lines (HCC1954, H1437, H2009), each sequenced with Illumina, PacBio
HiFi and ONT, and processed them with the identical pipeline on HPRC v2.1. Of the 41,072 HCC1395 truth alleles, 32–37 k reach a
tensor on each platform. [待补: results of the held-out-line models (R1.7) and of the comparators]

Two caveats shape how we use these lines.

- DeepSomatic was trained on these lines, and the CASTLE truth sets of HCC1954, H1437 and H2009 are its training truth.
  Comparisons with DeepSomatic there use only the chromosomes held out from its training [待核实: chromosomes]. HCC1395 has
  an independent truth set (SEQC2) and is our main external comparison.
- These truth sets were built from short-read call ensembles and omit most INDELs in long homopolymers and STRs (15.5% of
  SEQC2 INDELs lie in homopolymers of at least 7 bp, against 74.9% of HG008T INDELs). An INDEL candidate there without a
  truth allele is therefore not treated as a negative, and INDEL results on these lines are reported separately.

[待补, 若做: HCC1937 completes the five DeepSomatic lines.]

**现有 pipeline 能否解决：部分。**

- 已有：
  - HCC1395 三套 tensor（d46，rule (a)，germline QUAL ≥ 10 relabel 在跑）；
  - H1437 PacBio 已 finalize；HCC1954 / H2009 PacBio 在 finalize；其余在 sort / build（2026-10-12 squeue）。
- 缺：模型和评估（E12）；比较器（E8）；DeepSomatic held-out 染色体的核实；HCC1937。

**实验 / 功能设计。**

- E12（跨样本）和 E8（比较器）就是这条的实验。
- **E8b DeepSomatic held-out 核实**（文献，P0）：查 DeepSomatic 论文和 model card，确定每个训练细胞系留出的染色体。如果
  留出的是 chr1，就和我们的测试染色体一致。
- **E12b HCC1937**（可选，P2）：下载 HiFi / ONT / Illumina BAMs（约 0.6 TB），按 `tmp/castle_hcc1395_20261006` 的流程处理。

### R1.m2 Tensor windows truncated at node boundaries

**Comment.**

- Windows are restricted to the current node and truncated at node boundaries, so adjacent-node and path-spanning evidence is
  discarded.
- Short or single-base nodes are hard to handle, and many alternative nodes are shorter than the window.
- A credible graph-native method needs path-aware tensors spanning neighbouring nodes and explicit modelling of edges,
  bubbles, alternative paths and competing path support.

**评估（内部）。** 前半部分描述的是第一版，现在已经解决：窗口以位点为中心，每行沿 read 的 walk 跨节点。后半部分（显式
bubble / competing path support）没有，可选做 E24。

**Response.**

This described the first version, and the revised tensors address it.

- The 101-column window is centred on the site, and each row follows its own read's walk across node boundaries, on
  reference and branch nodes alike.
- The graph base (channel 5) and haplotype count (channel 6) are filled per column from the node the read traverses there.
  Single-base nodes therefore keep their full flanking context, and node transitions inside the window are visible to the
  network.
- These transitions carry signal: true somatic tensors have a median of 0 transitions within ±10 columns of the site and
  germline tensors 2 (AUC 0.72).
- Rows are ordered by the similarity of their edits and node transitions, so reads taking the same alternative path sit
  together. Competing path support is thus visible through the rows' walks, although not as an explicit bubble model (R1.1).

A schematic of a tensor crossing node boundaries is added as Supplementary Fig. [待补].
[可选, 若做 E24: a locus mode adds the reads of competing branches with their haplotype frequency (R1.1).]

**现有 pipeline 能否解决：能（已实现）。**

- 需要：一张跨节点 tensor 的示意图（`scripts/visualize_tensor.py`）。
- 可选：E24。

---

## 4. Reviewer 2

### R2.1 Why include the HapMap benchmark

**Comment.** The HapMap positives are germline variants whose donors are in the HPRC graph, so in theory the pangenome
alignment removes them. Discuss why the benchmark is included and its limits for modelling somatic mutations.

**评估（内部）。** 同 R1.5。R2 的语气温和（"discussion"），但 R1 和 R3 都提了，统一处理：主结果去掉 HapMap，用 dilution 替代。

**Response.**

We agree. The positives of the HapMap mixture are germline alleles, and four of its six donors are in both HPRC graphs, so the
benchmark tests the graph against its own design. We removed it from the main benchmark and replaced it with purity
dilutions that contain true acquired mutations: COLO829T on three platforms and HCC1395 on Illumina (R1.5). [二选一: HapMap no
longer appears / appears only in the Supplement on a donor-removed graph, with its limitation stated.]

**现有 pipeline 能否解决：部分。** 同 R1.5（E14、E23）。

### R2.2 How the truth sets were generated

**Comment.** Briefly describe how each benchmark's truth was generated, in particular to support the claim that truth sets
may contain germline variants and sequencing artefacts.

**评估（内部）。**

- 写 Methods 就行。
- 但 "truth 里有 germline 和 artefact" 这个说法要收敛：我们有直接证据的只是少数位点。改成 "benchmark ambiguities"，并给出
  组装证据。
- 每个 truth 的生成方式按原文献写，下面列的要逐一核实。

**Response.**

We added a description of each truth set (Supplementary Note) [待核实: each derivation against the original publication]:

- *HG008T:* the GIAB HG008-T somatic small-variant draft benchmark v0.2, derived from multi-technology tumor data against the
  curated diploid HG008-N assembly. Our confident region is its benchmark BED intersected with the normal's dipcall regions.
- *COLO829T:* the SMaHT validated SNV and INDEL union, with Illumina and PacBio VAFs, evaluated over the SMaHT easy,
  difficult and extreme regions.
- *HCC1395:* the SEQC2 high- and medium-confidence somatic SNVs and INDELs (v1.2.1), built from multi-site, multi-platform
  WGS and multiple callers, inside the SEQC2 high-confidence regions (v1.2).
- *HCC1954, H1437, H2009:* the CASTLE somatic truth sets and high-confidence regions released with DeepSomatic.

The germline statement is now supported directly and stated more carefully.

- 209 HG008T somatic truth alleles equal an HG008-N germline allele, through loss of heterozygosity or conflicting benchmarks.
- Of the 110 graph-absorbed HG008T alleles present in the normal, 109 equal the germline allele of the other normal haplotype.
- At 33 COLO829T truth INDELs, the COLO829BL assembly carries the ALT over the whole repeat.

We describe these as benchmark ambiguities at a small number of sites, not as errors of the truth sets, and we no longer use
them to explain the lower recall.

**现有 pipeline 能否解决：能。**

- 数字来自 `analysis/truth_labels.md` §8.1 第 5 点和两个 graph_absorbed 分析。
- 需要：按原文献写 Methods。

### R2.3 Explaining the lower recall; normal sample at missed sites

**Comment.**

- Recall is ~20% (COLO829T) and ~10% (HG008T) lower than DeepSomatic; the truth-error explanation should be expanded.
- Explain the false-negative categories "variant node not ranked #1" and "first aligned node mismatch" (Supp. Fig. 2).
- For false negatives that overlap existing graph nodes, check whether reads of the normal sample map preferentially to the
  alternative paths.

**评估（内部）。**

- 第一点：旧的类别来自旧 pipeline，已经被 read 层面的漏检分类取代，可以直接回答。
- 第二点：审稿人要 "normal reads"，我们有更强的证据：单倍型分辨的 matched-normal 组装。
- 可选补一个 normal reads 的实验（E20）：COLO829BL 的 ONT reads 在本地；HG008-N 的 reads 不在本地，要下载。
- recall 差距本身要等 E8 在同一协议下重新量。

**Response.**

Those categories came from an earlier pipeline. We replaced them with read-level classes, computed for every missed truth
allele on six datasets:

- no read with MAPQ > 10;
- fewer than three ALT reads, or AF below the candidate threshold;
- the ALT is an existing graph allele (reads align without edit to a graph path);
- repeat absorption or respelling (the reads carry the event, but the aligner places it as a different edit);
- an indel longer than 50 bp.

For example, of 398 HG008T PacBio SNVs without a candidate, 65 are graph SNP alleles and 239 are repeat absorptions; of its
2,976 INDELs without a candidate, 92% are fully or partly graph-encoded.

For the alleles on graph paths, we checked the normal with haplotype-resolved matched-normal assemblies, which show the
normal's sequence directly. Only 4.4% (HG008T) and 6.8% (COLO829T) of these alleles are present in the normal; the rest are
alleles of other HPRC individuals (R1.4). The HG008-T tumor assembly carries the ALT at 95% of graph-encoded HG008T INDEL
truths. These are real somatic events hidden by the graph, and we now say so.
[可选, 若做 E20: Normal reads agree: COLO829BL ONT reads aligned to the same graph take the alternative path at [待补] of the 34
germline sites and at [待补] of the other sites.]

The recall difference to DeepSomatic-TO is re-measured under the identical protocol (R3.m2) and decomposed by stage (R1.6).
[待补]

**现有 pipeline 能否解决：能（已完成；E20 可选）。**

**实验 / 功能设计。**

- **E20 normal reads 在吸收位点上的走向**（新脚本，可选，P2）。
  1. 位点：COLO829T 的 497 个吸收位点（GRCh38 坐标），±5 kb。
  2. 从 COLO829BL ONT 的 GRCh38 BAM（`/scratch/jshen/data/COLO829T/COLO829BL_ONT/*.bam`，本地）取这些区间的 reads（`samtools
     view -L`），转回 FASTQ。
  3. giraffe（long-read 参数，用和 d9 GAM 相同的 vg 版本 [待核实]）比对到 d9。
  4. 每个位点：normal reads 走 truth allele 那条 branch 的比例。
  5. 和组装分类对照（34 个 germline、430 个其他个体 allele）。
  6. HG008-N 的 reads 不在本地（只有 HG008-T），要下载 GIAB HG008-N 的数据才能做，用户决定。
  7. 代价：半天到一天。

### R2.4 Low-VAF detection and the low-VAF model

**Comment.**

- The 10% read-support threshold limits low-VAF detection in real samples with low purity or high heterogeneity. How should
  true low-VAF SNVs be detected?
- The specialised low-VAF model is poorly described: how does it differ from the main model?
- Does it recover low-VAF variants? At what VAF does the main model stop performing well? Run sensitivity analyses over
  decreasing VAF.
- When should each model be used?

**评估（内部）。**

- 成立。
- 旧的 low-VAF 模型是用 HapMap 训练的，随 HapMap 撤回。
- 要做三件事：
  1. 按 truth VAF bin 报 sensitivity：现有 per-truth 表里有 VAF，加上 predictions 就能出。
  2. dilution：COLO829T 三平台，HCC1395 Illumina（E14）。
  3. 一个低 AF 配置：candidate AF 0.03（E15）。参数在 builder 里都已经有，不用改代码，但 tensor 会多很多。
- 使用建议按 E14 / E15 的结果给出（VAF 低于多少时默认模型 recall 掉到 50 % 以下）。

**Response.**

*Thresholds.* Candidates now require AF ≥ 0.08 and at least three ALT reads with base quality ≥ 10, for SNVs and INDELs alike.
We report sensitivity per truth-VAF bin for every dataset. On COLO829T, the candidate threshold removes 74–82% of truth SNVs
with VAF < 0.06, while 89–99.5% of those with VAF ≥ 0.10 reach a tensor. [待补: recall per VAF bin after calling, all datasets]

*Dilutions.* For low purity, we trained on COLO829T dilutions (25, 50 and 75% purity) together with the undiluted tumor and
report results per purity. The HCC1395 dilutions (5–50%) test the same models on a tumor held out from training. [待补]

*Low-AF configuration.* The earlier HapMap-trained low-VAF model is withdrawn with the HapMap benchmark. It is replaced by a
low-AF configuration that differs from the default in one respect: candidates are built at AF ≥ 0.03 [待定: and two ALT
reads]. It is trained on the dilutions with the same network and recipe. [待补: recall and precision per VAF bin, default
against low-AF]

*Usage.* The default model is recommended when the expected somatic VAF is at least [待补: the VAF at which its recall falls
below 50%], and the low-AF configuration below that, at the cost of [待补]-fold more candidates and lower precision.

**现有 pipeline 能否解决：部分。**

- 已有：
  - VAF 分析（`analysis/tensor_recall_20260930` 的 per-truth 表有 VAF_Ill / VAF_PB）；
  - 九套 mixture tensor；fiberseq 联合模型（P75 / P50 异常没查）；
  - builder 的 `--snv-min-af`、`--indel-min-af`、`--min-variants` 选项。
- 缺：E14a–c、E15、VAF bin 的召回曲线（E1.6 的一部分）。

**实验 / 功能设计。**

- **E15 低 AF 配置**（不改 pipeline，用已有选项；P1）。
  1. 数据：
     - COLO829T Illumina 的 P100 / P75 / P50 / P25（d9；如果主图换 d46 就用 d46 版）；
     - HG008T Illumina（深度最高）。
  2. build：`--snv-min-af 0.03 --indel-min-af 0.05`，`--min-variants` 3 和 2 各试一次。先只建 chr1 + chr2–8 估算 tensor
     增量；增量可接受再建全基因组。
  3. 训练：同一 SNV recipe，P100 + 三个 dilution 联合训练；validation 同样按 recall 0.9 选阈值。
  4. 比较：在 chr1 上按 truth VAF bin（< 0.05、0.05–0.1、0.1–0.2、≥ 0.2）比较默认模型和低 AF 配置的 recall / precision，
     PoN 前后都看。
  5. 输出：使用建议里的 VAF 分界点。
  6. 代价：tensor 预计多 3–5 倍，build 和训练时间相应增加。
- E14（dilution）见 R1.5。

### R2.5 What drives the higher precision

**Comment.**

- Is the precision gain mostly reduced germline contamination, or are there other factors?
- Beyond modelling germline variation at the site, does graph-aware calling help when several germline variants are nearby?

**评估（内部）。**

- 第一点已有一半：PoN 前 HG008T Illumina 的 SNV FP 99 % 是 matched-normal germline。缺的是 linear 臂的对照（E3），以及把
  FP 按 "germline allele 是否在图里" 分开（E4b，同时回答 R3.2）。
- 第二点要一个新分层：按 ±50 bp 内的 germline 数目分层（E19）。预期 graph 在 germline 密集区帮助更大：read 走对 haplotype，
  错配少。但这只是预期，要看数据。

**Response.**

We decomposed the false positives.

- Before population filtering, 99% of HG008T Illumina SNV false calls are matched-normal germline variants.
- The graph removes most germline alleles before classification: on identical reads, 76–82% of germline truth alleles produce
  a candidate on a GRCh38-only graph, against 7–13% on HPRC v1.1 and v2.1 (R1.2).
- Population filtering then removes 95.5–96.6% of the remaining germline SNV tensors. The false positives left are mostly
  rare germline alleles (gnomAD AF < 0.001).

To link this to the graph, we split the false positives of the GRCh38-only and HPRC arms by whether the germline allele is
encoded in the HPRC graph (R3.2). [待补: expected that the precision gain is concentrated in germline alleles the graph
encodes]

For clustered germline variation, we stratified precision and recall by the number of matched-normal germline variants within
±50 bp of each call or truth allele (0, 1, 2–3, ≥ 4), on the HPRC and GRCh38-only arms. [待补]

**现有 pipeline 能否解决：部分。**

- 已有：FP 组成（`analysis/v1_vs_v2_20260929`、`analysis/tensor_recall_20260930` 的 PoN 表）。
- 缺：linear 臂的 calls（E3），E4b，E19。

**实验 / 功能设计。**

- **E4b FP 按 "germline allele 是否在图里" 分解**（新脚本，P1，依赖 E3）。
  1. 取 linear 臂和 HPRC 臂各自 PoN 前后的 PASS FP（rtg 的 fp.vcf）。
  2. 和 matched-normal germline 匹配：HG008-N dipcall、COLO829BL dipcall；CASTLE / HCC1395 用 DeepVariant normal（QUAL ≥ 10）。
  3. 对每个 germline FP 判断 allele 是否在 HPRC 图里：exact allele 在 d9 / d46 的 deconstruct VCF 里
     （`AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz`；d46 没有现成 VCF，用 `vg deconstruct` 生成，或用 labels 的 truth
     placement 在图上找 key）。INDEL 做 placement 归一。
  4. 表：FP 分成 germline-in-graph / germline-not-in-graph / non-germline 三类，每个参考臂一列。HPRC 的精度提升应集中在第一类。
  5. 同时报 non-germline FP（artefact）在两臂的差别，回答 "other factors"。
- **E19 按附近 germline 密度分层**（新脚本，P1，依赖 E3）。
  1. 对每个 truth allele 和每个 FP，数 ±50 bp 内的 matched-normal germline 变异数（dipcall）。
  2. 分 0 / 1 / 2–3 / ≥ 4 四档，算各档的 recall 和 precision，HPRC 臂和 linear 臂对比，加 bootstrap CI（E16）。

### R2.6 Panels of normals

**Comment.**

- Some PoNs were built from common benchmark datasets and include variants of the tested samples (HapMap in particular), which
  may inflate the reduction in germline contamination.
- Fig. 2e shows that > 20% of germline variants are not in the AF-HPRC graph: were they removed by the PoNs?
- Only COLO829 was evaluated without the PoN; add HG008T without the PoN.

**评估（内部）。**

- 第三点已经有数字。
- 第二点要一个小分析（E4c）。
- 第一点要逐个核实这四个资源的样本：
  - HapMap donors 是 1000G 样本，1000G PoN 很可能包含它们，dbSNP 也会有它们的变异。
  - HG008（2023 年后采集）和 COLO829 应该不在这些资源里，但 CoLoRSdb v1.1 收了哪些长读长样本要查。
- INDEL 规则是在 chr1（测试染色体）上选的，这是一个审稿人可能抓的问题。要在 validation 染色体上重新确认（E17c）。

**Response.**

1. *Overlap with the benchmark samples.* HG008 and COLO829 are not among the individuals of gnomAD, dbSNP, the 1000 Genomes
   panel or CoLoRSdb [待核实: sample lists of each resource and version]. The HapMap donors are 1000 Genomes samples and are
   represented in these resources, which is another reason we removed that benchmark (R1.5).
2. *Germline variants not in the graph.* They remain candidates and are labelled germline in training, so the classifier
   learns to reject them. Population filtering removes 95.5–96.6% of germline SNV tensors on all six datasets. Split by graph
   membership, it removes [待补]% of germline tensors whose allele is not in the graph, and the germline false positives left
   after filtering are [待补]% graph-absent.
3. *HG008T without filtering.* Chr1 SNV F1 is 0.18–0.20 without filtering and 0.59–0.62 with it; precision is 0.10–0.11
   against 0.45–0.50, and recall 0.86–0.87 against 0.83–0.84. [待补: final models] All results are now reported both with and
   without population filtering.

Population filtering cannot reach calls without a GRCh38 position; these are handled separately (R3.4). We renamed these
resources "population and background filtering resources" (R3.m4). INDELs use a separate, documented rule: a gnomAD or
CoLoRSdb exact-allele match at population AF ≥ 0.01 for length changes of homopolymers of at least 7 bp, and ≥ 0.001
otherwise. The rule was selected on validation chromosomes and applied unchanged to chr1. [待补: E17c；如果在 validation 上选出
来的规则不同，就改这里]

**现有 pipeline 能否解决：能（E4c、E17c 是小分析）。**

**实验 / 功能设计。**

- **E4c 不在图里的 germline 和 PoN**（新脚本，P1）。
  1. 每套数据的 germline truth tensor（label 2），按 allele 是否在图里分两组（同 E4b 第 3 步）。
  2. 两组分别统计被 PoN 标记的比例。
  3. 最终 PASS 里的 germline FP 也按同样两组统计。
- **E17c 在 validation 染色体上重选 INDEL PoN 规则**（新脚本，P1）。
  1. 对每个 INDEL 模型，在 validation blocks（chr2–22 的 5 % node blocks；`data.block_split` 是确定性的，可以重放）上跑
     predict，写 graph VCF → linear VCF。
  2. 在这些 calls 上扫 PoN 规则：
     - 全局 AF floor：0.001 / 0.005 / 0.01 / 0.02；
     - homopolymer 分档规则：HP ≥ 7 用 0.01、其他用 0.001，以及几个邻近组合。
  3. 选 validation F1 最高的规则，再在 chr1 上只评估一次。
  4. 注意 HCC1395 / CASTLE 的 truth 里 HP ≥ 7 INDEL 很少，只能确认非 homopolymer 那一档；HP 那一档要用 HG008T 和 COLO829T
     的 validation blocks。
  5. 代价：每个模型一次 validation 推理加 rtg。
- **E8c 资源样本核实**（文献，P0）：gnomAD 用的是哪一版（GATK bundle 的 af-only-gnomad.hg38），dbSNP 138，GATK 1000G PoN 的
  样本表，CoLoRSdb v1.1 的样本表。

### R2.7 Alignment statistics overstated

**Comment.** Lines 292–294 (more perfect reads, less soft clipping on the graph) look overstated given Supplementary Fig. 4;
give read counts or proportions.

**评估（内部）。** 成立。node 层面的计数已有，per-read 的三参考对比要扫 GAM（E11a）。

**Response.**

We replaced the figure with counts and proportions and toned down the statement (R3.6). For example, HG008T PacBio has 3.93
billion node visits by MAPQ > 5 alignments on HPRC v1.1; 8.0% are on branch nodes, and 98.0% of those are edit-free, against
92.4% on GRCh38 nodes. For the same reads on the GRCh38 BAM, the GRCh38-only graph, HPRC v1.1 and v2.1, we now report the
number and proportion of reads that are unmapped, aligned without edit, with at least one mismatch, with an indel, and with
soft clipping, together with MAPQ distributions, per platform. [待补]

**现有 pipeline 能否解决：部分。**

- node 层面已有（`analysis/branch_node_alignments_20261002`）。
- 同一批 reads 在三个参考上的 per-read 统计要扫三个 GAM（E11a）。

**实验 / 功能设计。**

- **E11a 同一批 reads 的 per-read 统计**（新脚本，P1）。
  1. 每个平台从 HG008T（以及 COLO829T）抽 2 M 个 read 名，要求在 GRCh38 BAM 和三个 graph GAM 里都有（输入 FASTQ 相同，
     read 数已核对一致）。
  2. GAM 用 v4 的 `gam_reader` 解析（MAPQ、edits、soft clip），BAM 用 pysam。
  3. 每条 read 记录：是否比上、MAPQ、mismatch 数、indel 数、soft-clip 碱基数，以及是否完全匹配。
  4. 输出每个参考 × 平台的计数和比例表，MAPQ 分布图；按 GIAB stratum 分层（E10）。
  5. 代价：几小时 CPU（按 read 名从 sorted GAM 取，或顺序扫一遍取样）。

### R2.m Minor comments 1–8

**Response.**

1. Figure 2d and 2f show truth-set variants, not Pansoma output; they moved to the Supplement.
2. Figures are now cited in order (2e before 2f; Figure 4 panels follow the text, minor 7).
3. The COLO829T panel of Supplementary Fig. 2 is labelled.
4. Point labels in Supplementary Fig. 1 were moved off the points.
5. "Truth overlapping" was replaced by defined categories, for example a call within 10 bp of a truth allele with a different
   allele. [待核实: 第 13 页原文的意思]
6. The truth-set spectrum is shown beside each caller in Figure 4 (R3.5).
7. See item 2.
8. A *candidate* is one allele (node, start, REF, ALT). A *candidate site* is all alleles at one node position and variant
   class, which are encoded together in one tensor.

**现有 pipeline 能否解决：能（改图和文字）。**

---

## 5. Reviewer 3

### R3.1 Reference-level ablation

**Comment.**

- End-to-end benchmarking does not isolate the pangenome: Pansoma differs from the comparators in reference, alignment,
  candidate generation, tensorisation, architecture and post-processing.
- Run the Pansoma pipeline with a GRCh38 linear reference (or a single-path graph) against the full HPRC graph, keeping
  candidate extraction, tensors, architecture, training and filtering thresholds identical, to show how much of the gain the
  pangenome drives.

**评估（内部）。** 这是三位审稿人最核心的要求，也是决定性的实验（E3）。

- tensor 已全部完成：六套 × linear / d9 / d46，标签齐全。
- 缺：
  - 模型。d9 的 SNV 要先重训（E1），其余用同一 recipe。
  - 无 ch6 的臂，用来把 "比对 / candidate" 和 "haplotype 通道" 分开。
- 两个要处理的细节：
  - vg 版本不同（linear 1.66，d46 1.77，d9 [待核实]）。简版原来写的 "same vg giraffe version" 已改（2026-10-12）。
  - linear 图的节点是 32 bp，节点数和 tensor 的 node block 划分不同，所以 validation 块和三臂不完全一样。chr1 测试集是
    同一批 reads、同一个 truth，可比。
- 结果的三种可能和对应写法见下面的判定表，E3 出来后和用户一起定。

**Response.**

We performed this ablation. The reads of all six datasets were aligned with vg Giraffe to three references and processed by
identical code:

- a GRCh38-only graph built by vg autoindex from the same GRCh38 FASTA (a single-path-equivalent graph);
- HPRC v1.1 (frequency-filtered Minigraph-Cactus graph, 90 haplotypes);
- HPRC v2.1 (464 haplotypes).

Candidate extraction (thresholds, read cap, normalisation), tensor construction, labelling rules, network, training recipe,
threshold rule and population filtering are identical across arms. Only the graph and its haplotype count differ. Giraffe
versions are 1.66 for the GRCh38-only graph and [待核实] / 1.77 for the HPRC graphs; a control realignment with one version
shows that this does not affect the comparison [待补, E3 步骤 6].

Candidate-level results are in R1.2 and R3.2.

[待补: Table — chr1 precision, recall and F1 per arm, dataset and variant type, before and after population filtering, with
block-bootstrap intervals for the differences between arms (R3.9)]

On the GRCh38-only graph the haplotype-count channel is constant, so the comparison GRCh38-only against HPRC includes both the
alignment and candidate effect and the population channel. An HPRC arm trained without channel 6 separates the two. [待补]

[待补: a cross-sample check — the COLO829T model of each arm applied to HG008T chr1]

**现有 pipeline 能否解决：部分（tensor 全部完成，模型没训）。**

- 已有：18 套 tensor（`pansoma_v2_tensors/*_linear`、`*`（d9）、`*_d46`），标签齐全；training 支持 `--drop-planes path_count`。
- 缺：模型、评估、版本对照。

**实验 / 功能设计。**

- **E3 三个参考臂的 ablation**（不改 pipeline，P0）。
  1. 臂：
     - A0 GRCh38-only（vg 1.66 autoindex，H = 17，path cover）；
     - A1 HPRC v1.1 d9（H = 90）；
     - A2 HPRC v2.1 d46（H = 464）；
     - A1−ch6、A2−ch6（`--drop-planes path_count`）。
  2. 模型：SNV recipe（50.4 M，pow:0.75，lr 1e-4，20 epochs，val R0.9，off-ref rescue）和 INDEL recipe（8.7 M，sqrt，
     lr 1e-4，12 epochs）。每个臂 × 6 套 × 2 种，A1 就是 E1 / E2。
     - 共 A0 12 + A1 12 + A2 12 + −ch6 臂（先只做 A1 的 SNV 6 个）= 42 个 run。
     - off-ref rescue 只在 HPRC 臂有意义（linear 臂没有 off-reference 节点）。
  3. 评估：
     - chr1 BED，raw 和 PoN 后，PASS 和曲线；
     - 完整输出（E9e）；
     - 按区域分层（E10）；
     - 各臂之间的配对 bootstrap 差值（E16）。
  4. 跨样本检查（不用新训练）：每个臂的 COLO829T 模型直接预测同平台 HG008T 的 chr1，看两臂的差在跨样本时是否保持。
  5. ch6 编码按各图的 H，三臂各自一致（v4 README §8）。
  6. 版本对照：
     - 用 vg 1.77 重建 linear 图索引（autoindex），重新比对 HG008T PacBio 和 Illumina 的 linear 臂；
     - 比较 candidate 层面的 truth 保留率和 chr1 结果。
     - 差 < 0.5 个百分点，就写 "version effect negligible"；否则把六套 linear 都换成 1.77 的比对。linear 比对很快，Illumina 约为
       d46 的 1/14。
  7. 代价：SNV 一个 run 约 5 h（一块 H100，linear 的 tensor 约为 d9 的两倍，时间更长），INDEL 小模型更短，可以用 MIG。合计
     约 200–250 GPU·h。
  8. 判定表：

     | E3 结果（PoN 后，CI） | 写法 |
     |---|---|
     | HPRC > GRCh38-only | 主张成立："the pangenome drives X points of F1"，再用 −ch6 臂拆成两部分 |
     | 两者持平，但 PoN 前 HPRC 明显好 | 图和 population filter 提供的是同一类信息。主张改为：图把这个先验提前到比对阶段，candidate 少 43–69 %，而且覆盖 PoN 到不了的 graph 坐标（前提是 E9 做好）；困难区域看 E10 |
     | HPRC < GRCh38-only | 先查原因（INDEL 吸收、off-ref FP），再决定；不能把这种结果藏起来 |

### R3.2 Composition of the removed candidates

**Comment.**

- Fewer candidates help only if germline variants and artefacts are removed preferentially while true somatic variants are
  kept.
- Compare linear and graph candidate generation for (1) somatic truth, (2) matched-normal germline, (3) background /
  non-truth candidates.
- Stratify the reduction of false positives by whether the germline allele is encoded in the HPRC graph, to connect the 79% /
  73% overlap of normal-sample SNVs / INDELs with the graph to the precision gain.

**评估（内部）。**

- 前三类的对比已经完成，现在是六套 × 三个参考臂（R1.2 的表），比简版更全。
- 背景 candidate 的数字在 `analysis/graph_ablation_tensors_20261005`（d46 和 COLO829T ONT linear 要重跑 `collect.py` 补进表）。
- "按 germline allele 是否在图里分层 FP" 要 E3 的 calls（E4b）。

**Response.**

We did this on identical reads for all six datasets and three references. The fraction of truth alleles that produce a
candidate tensor (full table in R1.2):

| | GRCh38-only | HPRC v1.1 | HPRC v2.1 |
|---|---|---|---|
| Somatic SNVs | 0.896–0.986 | 0.893–0.952 | 0.892–0.956 |
| Somatic insertions | 0.529–0.962 | 0.356–0.535 | 0.361–0.578 |
| Somatic deletions | 0.629–0.965 | 0.504–0.718 | 0.524–0.729 |
| Matched-normal germline variants | 0.757–0.822 | 0.077–0.131 | 0.072–0.124 |

On HPRC v1.1, total tensors fall by 43–69% and non-variant SNV candidates in the confident region by 11–61%
[待补: v2.1]. The graph therefore removes about nine in ten germline candidates while losing 0.3–4 points of somatic SNV
sensitivity. Somatic INDELs are the exception, because repeat-length somatic changes coincide with population alleles in the
graph (R1.4).

To connect graph representation to precision, we split the false positives of the GRCh38-only and HPRC arms into germline
alleles encoded in the HPRC graph, germline alleles absent from it, and non-germline artefacts. [待补: expected that the
reduction is concentrated in germline alleles the graph encodes, which links the 79% / 73% overlap to the precision gain]

**现有 pipeline 能否解决：部分。**

- candidate 层面完成（`truth_recall.json`、`analysis/graph_ablation_tensors_20261005`）。
- 缺：E4a 补表；E4b（依赖 E3）。

**实验 / 功能设计。**

- **E4a candidate 表补全**（脚本已有，P1）：重跑 `analysis/graph_ablation_tensors_20261005/collect.py`，加入 d46 六套和
  COLO829T ONT linear，并加 HCC1395 的 d46 列。
- **E4b** 见 R2.5。

### R3.3 Validation of graph-coordinate variants

**Comment.**

- A subset of final predictions exists only in graph coordinates; the evidence (signatures, read examples) is indirect.
- Validate a representative sample orthogonally: matched-normal alignment evidence, phased haplotype assemblies, or
  cross-platform data.

**评估（内部）。** 成立，而且对我们不利。

- 旧稿把 graph 坐标的 calls 当作亮点，但现在的数据显示它们大部分可能是 FP：
  - PoN 后 unplaced 占 SNV PASS 的 8–51 %；
  - 用同样的组装检验，没有 truth 匹配的 off-reference candidate 里只有 1.2–1.6 % 是 tumor-specific。
- 策略：把这些 calls 分三步处理，诚实报告每一步剩下多少。
  1. 能放回 GRCh38 的先放回去，交给 rtg 评估（E9a）。
  2. 剩下的做 graph 原生的 population filter（E9c）。
  3. 再剩下的做组装验证和跨平台检验（E9d）。
- HG008T 有 tumor 和 normal 组装，可以完整验证；COLO829T 只有 normal 组装和 COLO829BL ONT reads，只能排除 germline。

**Response.**

We added orthogonal validation and changed how graph-coordinate calls are handled.

1. *Placement first.* Many calls without a direct GRCh38 position lie on short branch nodes inside bubbles whose flanking
   nodes are on GRCh38. These are now placed on GRCh38 through the bubble anchors and evaluated against the truth like any
   other call (R3.4). [待补: fraction placed]
2. *Population haplotypes.* For the remaining calls, the ALT and its flanking sequence are compared with all HPRC haplotypes at
   the locus, including alleles removed by the frequency filter. A match marks a population allele, the graph-coordinate
   counterpart of population filtering (R3.4). [待补]
3. *Assemblies.* For each remaining call, the consensus of its ALT reads is aligned to the HG008-T tumor assembly and both
   HG008-N haplotypes (for COLO829T, both COLO829BL haplotypes and the COLO829BL ONT reads). A call is confirmed somatic when
   the ALT is in the tumor and absent from both normal haplotypes. [待补]
4. *Cross-platform.* The three platforms share one graph, so the same call (node, offset, allele) can be compared directly
   across platforms. [待补: reproducibility of confirmed versus unconfirmed calls]
5. *A gold set.* In the patient's own assembly frame, 450 HG008T somatic events lie on non-reference nodes, and 233 of them
   have an exact candidate tensor. We report how many are called. [待补]

The same assembly test, applied to off-reference candidates without a truth match, found only 1.2–1.6% tumor-specific. The
headline metrics therefore count every unconfirmed graph-coordinate call as a false positive (R3.4).

**现有 pipeline 能否解决：部分。**

- 已有：
  - `pansoma_net_v2_runs/assembly_offref_20260929`（450 个事件，233 个 exact key）；
  - `pansoma_net_v2_runs/offref_negatives_20260929`（ALT read consensus 比对到组装的脚本）；
  - 组装：HG008-T v3.2 hap1 / hap2、HG008-N v6.2 hap1 / hap2、COLO829BL dipcall；COLO829BL ONT BAM。
- 缺：E9a、E9c、E9d 在最终 calls 上的运行。

**实验 / 功能设计。**

- **E9d 组装验证和跨平台**（新脚本，P0）。
  1. 输入：E9a、E9c 之后仍然没有 GRCh38 位置、也没被 population 标记的 PASS calls（HG008T 三平台，COLO829T 三平台）。
  2. 对每个 call 取 ALT reads 的 consensus（±50 bp），用 minimap2 比对到 HG008-T v3.2 hap1 / hap2 和 HG008-N v6.2 hap1 / hap2
     （复用 offref_negatives 脚本）。
  3. 分类：
     - tumor-only（confirmed somatic）；
     - 在 normal 里（germline）；
     - 两边都没有（artefact）；
     - ambiguous（多处比对、identity 低）。
  4. COLO829T：
     - 用 COLO829BL 两个 haplotype 排除 germline；
     - COLO829BL ONT reads 在同一节点的 ALT 支持作为 normal 证据（reads 用 E20 的方法比到 d9）；
     - 没有 tumor 组装，所以只能给 "not germline"。
     - 可选：下载 SMaHT 的 COLO829T 组装（用户决定）。
  5. 跨平台：同一 node:offset:ALT 在 ≥ 2 个平台都是 PASS，按上面的类别统计复现率。
  6. gold set：233 个 exact key 里被最终模型 call 到的比例（graph-only somatic 的 sensitivity）。
  7. 输出：Supplementary Table；完整输出指标（E9e）里只把 confirmed 算 TP。

### R3.4 Metrics exclude graph-coordinate calls; "about 2%"; CHM13

**Comment.**

- TP/FP/FN are computed only for calls that project to GRCh38, so precision, recall and F1 describe a subset, not the complete
  output.
- The graph-coordinate fraction is platform dependent and large (> 10%, sometimes > 20% for ONT in Supplementary Table 6),
  which contradicts the Discussion's "about 2%".
- Are graph-coordinate calls truly graph-dependent, or just not projectable to GRCh38? Determine what fraction can be
  represented on T2T-CHM13 and what remains graph-specific.

**评估（内部）。** 成立，是现在最大的短板。

- PoN 后 unplaced PASS 占 SNV PASS 的比例：HG008T PacBio 8 %、Illumina 32 %、ONT 41 %；COLO829T Illumina 12 %、fiberseq
  13 %、ONT 51 %。全部算 FP 时 precision：HG008T 0.28–0.46，COLO829T ONT 0.36。
- INDEL 也有：HG008T Illumina INDEL raw PASS 里有 534 个 unplaced。
- 原因有两个：
  - off-ref rescue 把 site 的 ch6 设成 90，去掉了 population 先验；
  - PoN 只认 GRCh38 坐标，过滤不到这些 calls。
- 要做的（E9 整套）：
  - a. anchor 投影：很多 unplaced 在简单 bubble 的 branch 节点上，其实能放回 GRCh38。
  - b. CHM13 投影。
  - c. graph 原生的 population filter。
  - d. 组装验证。
  - e. 完整输出指标。
  - f.（可选）用组装给 off-ref tensor 打标签，替代 rescue。
- 已有原型：`hprc_offref_filter_20260930`，旧 run 上的结果：
  - 严格规则（21-mer 在 v1.1 全图里找到）标掉 76–90 % 的 unplaced PASS calls，只碰到 9–10 % 的 on-reference TP（只用在
    off-ref calls 上时不影响 on-reference TP）。
  - 但那次 4 个组装验证过的 off-ref somatic 里有 3 个也被标掉了，样本太少，要用 233 个的 gold set 重新校准。
- CHM13：d9 的 node → CHM13 坐标库是完整的，25 条染色体都在（`AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13.path_index.sqlite`，
  2026-10-12 查过）。它的 GRCh38 部分缺 chr3–9、X、Y、M，但 GRCh38 用 pipeline 自己的 `grch38_path`。d46 要另建。

**Response.**

We agree, and the revised evaluation covers the complete output.

1. *More calls are placed on GRCh38.* A call without a direct GRCh38 position on a branch node is now placed through the
   GRCh38 nodes that flank its bubble, by aligning the branch to the GRCh38 interval between them; it is then scored by rtg
   vcfeval like any other call. [待补: fraction of graph-coordinate calls placed this way, per platform]
2. *Population filtering of the rest.* Calls that still have no GRCh38 position are compared with all HPRC haplotypes at the
   locus, including alleles below the graph's frequency filter. [待补: fraction tagged]
3. *Complete-output metrics.* Every run reports two results: the GRCh38-evaluable calls, and the complete output, in which
   every remaining graph-coordinate PASS call counts as a false positive unless the assembly test confirms it (R3.3).
   [待补: complete-output precision, recall and F1, per platform]
4. *"About 2%".* That statement held for one dataset only. Before these steps, graph-coordinate calls were 8–51% of PASS SNV
   calls after population filtering, highest for ONT. [待补: after]
5. *CHM13.* CHM13 is a path in both HPRC graphs, so each remaining call is projected onto CHM13 directly or through its bubble
   anchors. We report the fraction representable on GRCh38 (directly or through anchors), on CHM13 only, and on neither, which
   is the part intrinsic to the graph. [待补]

**现有 pipeline 能否解决：不能（关键部分要新功能）。**

- 已有：
  - `linear_vcf` 把没有 GRCh38 位置的写进 `*.linear.unplaced.vcf.gz`（CHROM = node，INFO START = offset）；
  - vcfeval 报告里有 unplaced 数和 "全算 FP" 的 precision；
  - 原型 `hprc_offref_filter_20260930`（s1–s5）；
  - d9 的 CHM13 坐标库。
- 缺：E9a、E9b、E9c 的正式版，E9e 的报告。

**实验 / 功能设计。**（这一组先写成独立的后处理脚本，验证好后再考虑并进 `linear_vcf` / `vcfeval`，那一步才算 "改 pipeline"。）

- **E9a bubble anchor 投影**（新脚本，P0）。
  1. 输入：`<run>/vcf_chr1/*.linear.unplaced.vcf.gz` 和 graph VCF 里每个 call 的 INFO PATH。
  2. 对每个 call，沿图的边找左右两侧最近的 GRCh38 节点（复用 s5_local 的 `anchors()`：BFS，不穿过 GRCh38 节点）。d9 用
     path_index 里的 edges 表；d46 要先从 GFA 的 L 行建一个边表（GFA 375 GB，流式扫一遍，只留 edges）。
  3. 取 ALT reads 支持最多的 walk 在两个 anchor 之间的序列，和 GRCh38 的 anchor 区间做全局比对（edlib / parasail）。
  4. call 的 offset 落在 match / mismatch 列上，就映射成 GRCh38 POS，REF 取 FASTA 的碱基。落在插入段里，或者该列的 branch
     碱基和 GRCh38 不同（call 正好在 germline 变异的碱基上），就写成复合 allele 并打标记。
  5. 输出 `<name>.anchored.vcf.gz`，和 linear VCF 合并后交给 vcfeval；剩下的写 `<name>.unplaceable.vcf.gz`，并记原因：插入
     序列在 GRCh38 上没有对应、嵌套 / 复杂 bubble、N 个节点内找不到 anchor。
  6. 输出：每个 run 的 unplaced → anchored / unplaceable 比例。anchored 的部分 PoN 也能用了。
- **E9b CHM13 投影**（新脚本，P1）。
  1. 对 unplaceable 的 call：节点在 CHM13 path 上就直接取 CHM13 坐标（d9 path_index）；否则按 E9a 的办法用 CHM13 anchor。
  2. d46：用 `tools/graph_prep ref-path-scan` 扫 CHM13 path 建坐标表 [待核实: `--reference-sample CHM13` 是否可用]。
  3. 输出四类比例：GRCh38 直接、GRCh38 via anchor、只在 CHM13、都不在（graph-specific）。
- **E9c graph 原生的 population filter**（新脚本，P0；只用于 unplaceable calls）。
  1. 每个 call 取 ALT 及两侧 10 bp（21-mer），来自 ALT reads 的 consensus。
  2. 在该 locus（anchor 区间 ±50 bp）的 HPRC v1.1 全图单倍型里找这个序列：用 raw VCF（88 条单倍型，包括 d9 过滤掉的 allele），
     原型是 s5_local。v2.1 需要下载全图 VCF。
  3. 规则：至少一条单倍型拼出 → 标 `PopulationHaplotype`（相当于 SNV PoN 的 "任何群体 allele"）。
  4. 校准：
     - 233 个 gold set 的 off-ref somatic 被标掉的比例要很小；
     - 在 validation 染色体上选 context 长度（11 / 21 / 31 bp）。
  5. 输出：每个 run 的 unplaced PASS 被标掉的比例，以及 gold set 的损失。
- **E9e 完整输出指标**（新脚本，P0；以后并进 vcfeval 报告）。
  - precision_complete = TP_rtg / (TP_rtg + FP_rtg + 剩余 unplaced PASS − confirmed)。
  - recall 仍按 rtg（truth 都在 GRCh38 上），另报 confirmed graph-only somatic 的个数。
  - COLO829T 没有 tumor 组装，剩余的 unplaced 全部算 FP。
- **E9f 用组装给 off-ref tensor 打标签**（改 pipeline，可选，P2）：见 R1.8。做了它就可以去掉 ch6 rescue，从根上减少这类 FP。
- **E9g off-ref 单独阈值**（新脚本，P1，E9f 的轻量版）：只在 validation blocks 上用组装标签（HG008T）评估 off-ref calls，选一个让
  off-ref precision 达到目标的阈值，不重新训练。

### R3.5 Mutational signatures with the truth set

**Comment.**

- In Figure 4a, show the truth-set spectrum side by side with Pansoma, DeepSomatic-TO and ClairS-TO, instead of using it only
  as the reference for cosine similarity; this matters most for HG008T, whose profile is heterogeneous and whose cosine
  similarities are modest.
- Provide the SBS96 profiles in the Supplement, to support the claim that graph-coordinate variants are biological.

**评估（内部）。** 小工作量，要等最终 calls（E1、E8、E9）。graph 坐标的 calls 的三核苷酸 context 要从 graph 上取。

**Response.** Figure 4 now shows the truth-set spectrum beside Pansoma, DeepSomatic-TO and ClairS-TO for each dataset. The
SBS96 profiles of every call set are in the Supplement: true and false positives separately, and graph-coordinate calls
separately, both before and after the steps of R3.4. [待补] For HG008T, the truth spectrum is a mixture of [待补: signatures],
which explains the modest cosine similarity of every caller.

**现有 pipeline 能否解决：能。**

**实验 / 功能设计。**

- **E18 SBS96**（新脚本，P1）。
  1. GRCh38 calls：取 FASTA 的三核苷酸 context。
  2. graph 坐标的 calls：用 ALT reads 的 consensus 或图上的 walk（节点加相邻节点）取 context。
  3. 用 SigProfilerAssignment（COSMIC v3.4）分解。
  4. 输出：truth / Pansoma / 两个比较器的 SBS96 和 cosine，分成 TP、FP、graph-only 三组。

### R3.6 More perfect alignment is not more correct alignment

**Comment.**

- More perfect alignments and fewer mismatches/clips do not prove more correct mapping, because the graph has more sequence
  to match.
- Compare alignment quality across GRCh38, T2T-CHM13 and the full HPRC graph, to separate reference completeness from
  multi-haplotype modelling.
- If feasible, use MAPQ distributions, phased-haplotype concordance, assembly-based truth, segmental duplications and other
  difficult regions.

**评估（内部）。** 成立，而且是很好的建议。

- 最干净的做法是模拟：从 HG008-N 组装（不在任何图里）模拟 reads，真实来源已知，在每个参考上测错位率和 MAPQ 校准（E11b）。
- CHM13 只做 alignment 对照，不做 calling：truth 都在 GRCh38 上，做 calling 要 liftover，得不偿失。
- 真实数据有一个现成的个人基因组参照：HG008-T reads 比到 HG008-N hap1 的 BAM（E11c）。

**Response.**

We agree, revised the wording, and added three analyses.

1. *Same reads, four references.* Per-read MAPQ, mismatch, indel and soft-clip distributions for identical reads on GRCh38
   (the comparators' BAMs), the GRCh38-only graph, HPRC v1.1 and v2.1, per GIAB stratification (R1.3, R2.7). [待补]
2. *Correctness with known origin.* We simulated Illumina, HiFi and ONT reads from the curated HG008-N diploid assembly, which
   is in none of the references, so the true origin of each read is known through the assembly-to-GRCh38 alignment. On
   GRCh38, a CHM13-only graph, the GRCh38-only graph, HPRC v1.1 and v2.1, we measured the misplacement rate and MAPQ
   calibration, overall and in segmental duplications, low-mappability regions and tandem repeats. This separates reference
   completeness (CHM13) from multi-haplotype modelling (HPRC). [待补]
3. *The patient's own genome.* Real HG008-T reads were also aligned to the HG008-N assembly. For reads in difficult regions we
   compare their graph placement with their placement on the personal genome. [待补]

**现有 pipeline 能否解决：部分。**

- 已有：
  - HG008-T reads 比到 HG008-N hap1 的 BAM（`HG008_GIAB/HG008T_align_to_HG008N_hap1`、`HG008T_PacBio_vs_HG008N_hap1`），只有
    hap1；
  - dipcall 的 hap → GRCh38 PAF（`HG008_GIAB/dipcall_HG008N_GRCh38/*.hap1.paf.gz` 等）；
  - 三个图的索引。
- 缺：CHM13-only 图；read simulator；E11b、E11c 的脚本。

**实验 / 功能设计。**

- **E11b 模拟 reads 的比对正确性**（新脚本 + 新索引，P1）。
  1. 来源：HG008-N v6.2 hap1、hap2（本地）。
  2. 模拟：
     - Illumina 2×150（wgsim / ART / mason，任选一个装得上的）；
     - HiFi 和 ONT R10：pbsim3（HiFi 用 CCS 模式，ONT 用 R10 误差模型）。

     每种约 5x 全基因组，或者 30x 只取 chr1 + 几个 segdup 多的染色体。
  3. 真实位置：用 dipcall 的 hap → GRCh38 PAF 做 liftover（paftools liftover）；CHM13 的位置用 minimap2 asm5 把 haps 比到
     CHM13 后同样 liftover。本地如果没有 CHM13 FASTA 要下载（用户决定）。
  4. 比对：
     - GRCh38：bwa-mem2 / minimap2 / pbmm2，和比较器一致；
     - GRCh38-only 图、d9、d46；
     - CHM13-only 图：`vg autoindex` 新建。

     graph 比对 `vg surject` 到 GRCh38（CHM13 图 surject 到 CHM13）。
  5. 指标：
     - 错位率（surject 位置和真实位置不重叠）；
     - 未比上率；
     - MAPQ 校准（各 MAPQ 档的实际错位率）。

     按 GIAB stratum 分层（E10）。
  6. 输出：平台 × 参考的表和 MAPQ 校准图。
  7. 代价：建 CHM13 图索引约 1 天；模拟和比对每个参考几小时到一天。
- **E11c 个人基因组参照**（新脚本，P2）。
  1. 取困难区域（segdup、low mappability）里的 HG008-T reads。
  2. 比较两种位置：
     - graph 位置（surject 到 GRCh38）；
     - hap1 BAM 里的位置，经 dipcall PAF 投到 GRCh38。
  3. 统计一致率。只有 hap1，所以只看落在 hap1 上的 reads。

### R3.7 Caution on reclassifying benchmark false negatives

**Comment.** Graph alignment alone cannot overturn established truth annotations; reclassifying the Supp. Fig. 3 false
negatives needs orthogonal evidence (matched-normal alignments or haplotype-resolved assemblies), otherwise soften the claims.

**评估（内部）。** 成立。我们现在有组装证据，结论也反过来了：这些多数是真的 somatic，被图吞掉了，而不是 truth 的错误。

**Response.** We agree and revised the text. We now provide orthogonal evidence from haplotype-resolved assemblies (R1.4): only
4.4% (HG008T) and 6.8% (COLO829T) of these alleles are present in the matched normal, and the HG008-T tumor assembly carries
95% of the graph-encoded HG008T INDELs. The evidence points the other way from our original interpretation: most of these
alleles are genuine somatic events hidden by population alleles in the graph, not benchmark errors. Supplementary Fig. 3 is
reframed accordingly, and the claim that they are reference-bias artefacts is withdrawn.

**现有 pipeline 能否解决：能（已完成）。**

### R3.8 Scope: tumor-only versus somatic mosaicism

**Comment.** Validation is tumor-centric (two tumors and a synthetic admixture), with no normal-tissue mosaicism dataset such
as SMaHT tissues. Either keep the scope to tumor-only calling consistently, or validate on an independently characterised
mosaicism dataset.

**评估（内部）。** 选 tumor-only。不建议做 SMaHT 组织：VAF 范围不同（常 < 1 %），需要新的训练数据和模型，工作量大，而且偏离
主题。

**Response.** Pansoma is designed and validated for tumor-only calling in cancer samples. We now state this scope consistently
in the title, abstract, results and Discussion, and removed broader claims about somatic variation. Somatic mosaicism in
normal tissue is discussed as future work, because it involves much lower allele fractions and needs its own training and
validation data.

**现有 pipeline 能否解决：能（改文字）。**

### R3.9 Model stability and statistical uncertainty

**Comment.**

- Results are point estimates; assess stability by multi-seed retraining, block-bootstrap confidence intervals over genomic
  intervals, or chromosome-wide variance.
- In particular, retrain across seeds and report the distribution of precision, recall and F1 on the held-out chromosome and
  the independent HG008T benchmark.
- Use chromosome- or block-level resampling for confidence intervals of the difference between Pansoma and the comparators.
- Report sensitivity to the fixed 0.5 threshold.

**评估（内部）。** 成立，都要做。

- 注意 "independent HG008T benchmark"：原稿是 COLO829T 训练、HG008T 测试；现在 12 个模型都是同一样本内的。跨样本测试（E12）
  就是对这一点的回答。
- 阈值：现在不用 0.5，用 validation recall 0.9，另外三种规则都报。

**Response.**

We now report uncertainty for every headline number.

1. *Repeated training.* Each headline model was retrained with [待定: 3 / 5] random seeds; we report the mean and range of
   precision, recall and F1 on the held-out chromosome and on held-out cell lines, including HG008T when it is held out
   (R1.7). [待补]
2. *Block bootstrap.* 95% intervals from 1,000 resamplings of 1-Mb blocks of the test chromosome, for each metric and for the
   paired difference between Pansoma and each comparator (the same blocks for both). [待补]
3. *Operating point.* The fixed 0.5 threshold is no longer used. The PASS threshold is chosen on validation blocks (recall
   0.9), never on the test chromosome. We show full precision–recall curves and the metrics at four operating points:
   validation recall 0.9 and 0.95, the best validation F1 threshold, and argmax. [待补]
4. *Chromosome-level variance.* For [待定: one dataset], a three-fold cross-validation over chromosome groups estimates the
   variance across test chromosomes. [待补, 可选]

**现有 pipeline 能否解决：部分。**

- 已有：train 有 `--seed`；vcfeval 已经输出曲线、chr1 best point，predict 报三种规则。
- 缺：多 seed 的 run；bootstrap 脚本；R0.95 点。

**实验 / 功能设计。**

- **E16 seed、bootstrap 和阈值敏感性**（新脚本 + GPU，P0）。
  1. seed：只对主结果的模型做，即 E12 的主 fold 和 E3 / E1 里进主表的模型。每个 3 个 seed（时间够就 5 个），约 18–30 个额外
     run。
  2. block bootstrap：
     - 把 rtg 输出（tp-baseline、tp、fp、fn VCF）按 chr1 的 1-Mb 块（约 248 块）分组；
     - 有放回抽 1,000 次，算 P / R / F1；
     - 和比较器比时，同一组块配对，算差值的 95 % percentile CI；
     - unplaced 的 FP 按 anchor 位置（E9a）放进块，放不进的按比例分摊。
  3. 阈值：curves 已有；加 val R0.95（graph_vcf 有 `--target-recall`）和 argmax、t 的结果表。
  4. 可选：一个数据集做三折染色体组交叉验证（3 个 run）。

### R3.10 Low-VAF and INDEL claims

**Comment.**

- Qualify the HapMap low-VAF results (germline lines, not somatic lineages).
- Separate SNV from INDEL performance; Pansoma's edge is in SNVs, while INDEL accuracy is modest and below ClairS-TO in some
  ONT settings. Refine overall claims.
- Complement HapMap with a bona fide low-VAF somatic benchmark or an in-silico dilution of validated somatic variants (COLO829T
  or HG008T).
- Stratify INDEL performance by insertion/deletion, length, sequence context (homopolymer, tandem repeat) and platform, to tell
  low-VAF biology, platform error modes and INDEL classes apart.

**评估（内部）。** 全部成立。

- dilution：E14。
- INDEL 分层已有初步结果（`tmp/indel_hp_strata_20261010`）。例如 Illumina 上 HP ≥ 7 和其他 INDEL 的 F1 相近（0.249 / 0.213），
  ONT 上 HP ≥ 7 的 precision 只有 0.006。要在最终模型上出完整表（E17a）。
- INDEL 的主要发现可以作为一个结果写：
  - INDEL 的限制主要来自 population allele 的重叠和长读长 homopolymer 错误，不是分类器；
  - 跨个体的 error panel 可以减少长读长的假 INDEL（E17b，可选）。

**Response.**

*Low VAF.* HapMap is out of the main text, replaced by true-somatic dilutions of COLO829T (three platforms) and HCC1395
(Illumina), with VAF-binned sensitivity (R1.5, R2.4).

*SNV and INDEL are separated throughout.* The abstract and Discussion now state the SNV and INDEL results separately, and INDEL
claims are limited to what the data support. Where ClairS-TO outperforms Pansoma on INDELs (for example ONT [待补]), we say so.

*INDEL stratification.* INDEL performance is stratified by insertion or deletion, length (1, 2–4, 5–15, 16–50 bp), sequence
context (homopolymer ≥ 7 bp, homopolymer 4–6 bp, STR with unit 2–6 bp and ≥ 3 copies, other) and platform, for Pansoma and
the comparators. [待补: final table]

The main finding is already clear.

- 74% of HG008T chr1 truth INDELs lie in homopolymers of at least 5 bp, and most are exact population alleles. Population
  filtering alone caps any tumor-only caller's INDEL recall at 0.35 on HG008T chr1 (R1.6).
- On Illumina, accuracy is similar inside and outside long homopolymers (F1 0.25 against 0.21 on HG008T chr1).
- On long reads, false INDEL calls are dominated by one-base homopolymer errors at allele fractions just above the candidate
  threshold. On ONT, precision in homopolymers of at least 7 bp is below 1% (HG008T chr1).

INDEL accuracy is therefore limited by the overlap with population alleles and by long-read homopolymer errors, more than by
the classifier. [可选, 若做 E17b: An error panel built from unrelated samples of the same platform removes [待补]% of these
false calls with [待补]% loss of truth.]

**现有 pipeline 能否解决：部分。**

- 已有：
  - `analysis/indel_repeat_context_predictions_20261006`；
  - `tmp/indel_hp_strata_20261010`（HP 分层）；
  - `tmp/indel_non_20261006`（长读长 INDEL non 的来源和 error panel 规则）；
  - `tmp/indel_pon_recheck_20261010`（PoN 规则对比）。
- 缺：E2 的最终模型、E17a 的完整表、比较器的同一分层（要 E8）。

**实验 / 功能设计。**

- **E17a INDEL 分层表**（新脚本，扩展 `strata.py`，P1）：
  - 维度：类别 × 长度 × context × 平台；
  - 指标：recall / precision / F1，PoN 前后；
  - 对象：Pansoma 和两个比较器；
  - context 用 truth_labels 的 `RepeatContext` 规则，和 rule (a) 一致。
- **E17b 跨个体 error panel**（新脚本，可选，P2）。
  1. panel：同平台、不相关个体的 non-variant INDEL allele（例如给 COLO829T 用 HG008T 和 CASTLE 的，给 HG008T 用 COLO829T 和
     CASTLE 的），只取 AF < 0.15 的。绝不用 matched normal。
  2. 作为后处理 FILTER 标签，在 validation blocks 上选规则，再在 chr1 上评估一次。
  3. 之前的分析：PacBio / fiberseq 的 non 少 25–33 %，truth 几乎不丢；ONT 只少约 11 %。
- **E17c** 见 R2.6。

### R3.m1 Parameter discrepancies between Methods and Supplement

**Comment.** Main Methods (AF > 10%, 2 supporting reads, MAPQ ≤ 5 excluded, 80 epochs) and Supplementary Methods (AF ≥ 0.08,
≥ 3 ALT observations, mean ALT BQ ≥ 10, MAPQ ≥ 10, INDEL ≤ 50 bp, 80–120 epochs) disagree; upstream and downstream steps are
condensed compared with the network.

**评估（内部）。** 成立。现在只有一套参数，都记录在 manifest 和 checkpoint 里，统一成一张表。Methods 里上下游步骤要写得和网络
一样详细（candidate、标签、VCF、PoN、vcfeval）。

**Response.** One parameter table (Supplementary Table [待补]) now serves the main and Supplementary Methods, and every value
is recorded in each run's manifest. The values are:

- node selection (discovery): alignments with MAPQ > 5; a node is a target when more than 5% of its alignments carry an edit;
- reads with MAPQ ≤ 10 are dropped from candidate counting;
- at least 3 ALT reads, each ALT observation with base quality ≥ 10;
- AF ≥ 0.08 for SNVs and INDELs;
- INDELs up to 50 bp;
- at most 800 reads per node (deterministic subsampling);
- tensors of 8 channels × 200 rows × 101 columns;
- SNV model 50.4 M parameters, 20 epochs; INDEL model 8.7 M parameters, 12 epochs; batch 1,024; learning rate 1e-4; the
  best epoch is selected on validation blocks.

The Methods now describe candidate generation, labelling, VCF conversion, population filtering and evaluation in the same
detail as the network.

**现有 pipeline 能否解决：能。** 参数都在 manifest 和 checkpoint 里（v4 README §4、§5；pansoma_net_v2 README）。

### R3.m2 Comparator setup

**Comment.** Document the comparators as precisely as Pansoma: software versions, aligners and versions, pretrained
checkpoints, exact tumor-only configurations, PASS/filter criteria, reference resources, population/PoN resources, command
lines. Precision, Pansoma's main advantage, is sensitive to filtering configuration.

**评估（内部）。** 成立，而且必须重跑。

- 本地只有 COLO829T fiberseq 的 DeepSomatic 输出和早期 ClairS 的东西，原稿的比较器结果不在这里。
- 要在所有数据集上按同一协议重跑（E8），这同时是 R1.6、R2.3、R3.9 的前提。
- 精度对过滤配置敏感，所以两种过滤都报：
  - 各 caller 自己的默认过滤；
  - 同一个 population filter。

**Response.** A new Supplementary Table lists, for ClairS-TO and DeepSomatic-TO: software version, aligner and version,
model checkpoint, tumor-only configuration, reference, population and panel-of-normals resources, PASS criteria and the
exact command lines. [待补] All callers were run on the same reads (GRCh38 BAMs from the same FASTQ) and scored with the same
truth, benchmark regions, chromosomes and evaluation tool. Each caller is reported twice: with its own default filters, and
with the same population filtering rule applied to Pansoma. Precision–recall curves use each caller's own quality score.
[待补]

**现有 pipeline 能否解决：不能（要重跑）。**

**实验 / 功能设计。**

- **E8 比较器同协议重跑**（新运行，不改 pipeline，P0）。
  1. caller 和版本：ClairS-TO 最新 release 和 DeepSomatic 最新 1.x 的 tumor-only 模型（WGS / PacBio / ONT）[待核实: 版本号、
     模型名]。用 Singularity 镜像，记录镜像 digest。
  2. 输入：每个数据集的 GRCh38 BAM。
     - HG008T Illumina：本地有 GIAB v3 BAM；
     - HG008T PacBio / ONT：GIAB 的 GRCh38 BAM [待核实: 本地有没有]，没有就用 pbmm2 / minimap2 从同一 FASTQ 比；
     - COLO829T：SMaHT BAM；
     - HCC1395 / CASTLE：已下载的 DeepSomatic bucket BAM。

     记录 aligner 和版本。
  3. 运行：全基因组（或 chr1 加 E10 需要的区域），命令行全部写进表。
  4. 过滤，三种都报：
     - 各自默认 PASS（ClairS-TO 自带 PoN 资源；DeepSomatic-TO 自己的 germline / PoN 处理 [待核实]）；
     - 对它们的 PASS 再加我们的 population filter；
     - 用 QUAL 画曲线。
  5. 评估：同一个 vcfeval wrapper（truth、BED、chr1、`--squash-ploidy --sample ALT,ALT`）。
  6. DeepSomatic 训练过的细胞系，只在它留出的染色体上比较（E8b）。
  7. 代价：DeepSomatic 每个样本 chr1 几小时（GPU 更快），ClairS-TO 用 CPU；全部数据集约 1–2 周集群时间。
- **E8c** 资源样本核实：见 R2.6。

### R3.m3 Composition of the pangenome reference

**Comment.** Document the HPRC release and version, included haplotypes, allele-frequency filtering, treatment of structural
variants and graph construction; state whether benchmark individuals or close relatives contributed to the graph.

**评估（内部）。** 改文字就行；v2.1 的人群组成可以从 `tmp/hprc_v21_samples/v21_samples_populations.tsv` 统计（2026-10-12 已算）。

**Response.** The Methods now specify both graphs.

- *HPRC v1.1:* the Minigraph-Cactus GRCh38-based graph, frequency-filtered by keeping nodes traversed by at least 9 of 90
  haplotypes (44 individuals × 2, CHM13 and GRCh38). The individuals are African 23, American 16, East Asian 4, South Asian 1
  and European 0 [待核实].
- *HPRC v2.1:* the Minigraph-Cactus GRCh38-based graph filtered at 46 of 464 haplotypes (462 HPRC haplotypes, GRCh38 and
  CHM13). Of the 227 individuals with population metadata, 68 are African, 50 East Asian, 44 American, 36 South Asian and 29
  European [待核实: the remaining individuals]. Its distance and minimizer indexes were rebuilt with vg 1.77.
- *GRCh38-only graph:* vg autoindex from the GRCh38 no-alt analysis set.

Structural variants remain in the graphs, but candidates are limited to indels of at most 50 bp; longer events are logged and
not called. HG008, COLO829 and the SEQC2 and CASTLE cell lines are in neither graph, and none of their donors are HPRC
individuals [待核实]. The four HapMap donors are in both graphs (R1.5).

**现有 pipeline 能否解决：能。**

### R3.m4 Terminology and language

**Comment.** Inconsistent terminology and several grammatical problems; "panel of normals" is too broad for gnomAD, dbSNP,
1000G and CoLoRSdb; specific sentences; NPU expansion; Pansoma vs Pansoma-TO; "Supplement Figure" vs "Supplementary Figure";
final language edit.

**Response.** We standardised the terms: Pansoma (not Pansoma-TO), "Supplementary Figure", and the node pileup (NPU) expansion
in the main text and Supplementary Methods. "Panel of normals" became "population and background filtering resources". The
quoted sentences were corrected ("One key challenge in calling variants on graph alignments is that they differ fundamentally
…", "Pansoma converts the variants …", "with the technical details described in Methods", "read pileup" in Figure 5), and the
manuscript was language-edited.

**现有 pipeline 能否解决：能（改文字）。**

---

## 6. 实验和功能总表

类型：**无代码** = 现有 job 脚本直接跑；**新脚本** = 只写分析脚本；**改 pipeline** = 改 `indexed_gam_pipeline_v4` 或
`pansoma_net_v2`（这次不动）；**下载** = 要下载数据（用户决定）。

| ID | 内容 | 回应 | 优先级 | 类型 | 依赖 | 估计代价 |
|---|---|---|---|---|---|---|
| E1 | d9 六套 SNV 在新 build 上重训（现行 SNV recipe） | R1.6、R2.6、R3.1，所有 SNV 数字 | P0 | 无代码（`jobs/run.sh`） | — | 6 run × 约 5 h H100 |
| E1.1 | 最终模型上的 ch6 counterfactual | R1.1 | P1 | 新脚本（复用） | E1 | < 6 GPU·h |
| E1.6 | recall 漏斗（candidate → PASS → PoN），按种类、VAF、重复区 | R1.6、R1.8、R2.3 | P0 | 新脚本 | E1、E2、E8 | CPU |
| E2 | d9 六套 INDEL，10-11 recipe（lr 1e-4，不 filter） | R1.6、R3.10 | P0 | 无代码（`jobs/submit_indel.sh`） | — | 6 run |
| E3 | 三个参考臂 + 无 ch6 臂，同 recipe；跨样本检查；vg 版本对照 | R3.1、R1.1、R1.2、R2.5、R3.2 | P0 | 无代码（版本对照要重建 linear 索引） | E1、E2 | 约 42 run，200–250 GPU·h |
| E4a | candidate 表补全（d46、COLO829T ONT linear、HCC1395） | R3.2、R1.2 | P1 | 脚本已有 | — | CPU |
| E4b | FP 按 germline allele 是否在图里分解 | R3.2、R2.5 | P1 | 新脚本 | E3 | CPU |
| E4c | 不在图里的 germline 被 PoN 去掉多少 | R2.6 | P1 | 新脚本 | E1 | CPU |
| E5b | d46、HCC1395 的吸收分类 | R1.4、R1.2 | P1 | 复用 graph_absorbed 脚本 | — | 1–2 天 CPU |
| E5c | 比较器在吸收位点上的 recall | R1.4 | P1 | 新脚本 | E8 | CPU |
| E6 | branch candidate 可行性（找回多少、代价多大） | R1.4 | P1 | 新脚本 | — | 1–2 天 |
| E7 | in-silico 图组成（haplotype 数 × 频率阈值） | R1.2 | P1 | 新脚本（v2.1 部分要下载 VCF） | — | 几小时 CPU |
| E7b | 一个实测点（v2.1 d23 图） | R1.2 | P2 | 下载 + 新索引 | E7 | 3–5 天 |
| E8 | ClairS-TO / DeepSomatic-TO 同协议重跑，全部数据集 | R3.m2、R1.6、R2.3、R3.9、R3.5 | P0 | 新运行 | GRCh38 BAMs | 1–2 周集群 |
| E8b | DeepSomatic 训练时留出的染色体（文献核实） | R1.m1、R3.m2 | P0 | 文献 | — | — |
| E8c | 四个 population 资源的样本表核实 | R2.6 | P0 | 文献 | — | — |
| E9a | unplaced calls 的 bubble anchor 投影 | R3.4、R3.3、R1.3 | P0 | 新脚本 | E1 | 2–4 天开发 |
| E9b | CHM13 投影 | R3.4 | P1 | 新脚本 | E9a | 1–2 天 |
| E9c | graph 原生 population filter（原型 `hprc_offref_filter_20260930`） | R3.4、R1.3 | P0 | 新脚本 | E9a | 2–3 天 |
| E9d | 组装验证 + 跨平台复现 + gold set | R3.3 | P0 | 新脚本（复用 offref_negatives） | E9c | 2–3 天 |
| E9e | 完整输出指标 | R3.4 | P0 | 新脚本 | E9a–d | 1 天 |
| E9f | 用组装给 off-ref tensor 打标签，去掉 rescue | R1.8、R3.4 | P2 | 改 pipeline | — | 1–2 周 |
| E9g | off-ref 单独阈值（validation + 组装标签） | R3.4 | P1 | 新脚本 | E9d | 1–2 天 |
| E10 | 按 GIAB / SMaHT 区域分层（calling + alignment + unplaced） | R1.3、R3.6、R2.7 | P1 | 下载 + 新脚本 | E3、E8 | CPU |
| E11a | 同一批 reads 在四个参考上的 per-read 统计 | R2.7、R3.6、R1.3 | P1 | 新脚本 | — | 几小时 |
| E11b | 从 HG008-N 组装模拟 reads，测错位率和 MAPQ 校准（含 CHM13-only 图） | R3.6、R1.3 | P1 | 新脚本 + 新索引（装 simulator） | — | 3–5 天 |
| E11c | 个人基因组参照（HG008-T reads on HG008-N hap1） | R3.6 | P2 | 新脚本 | — | 1–2 天 |
| E12 | 跨样本 / leave-one-line-out（d46，6 个细胞系 × 3 平台） | R1.7、R1.m1、R3.9 | P0 | 训练代码小改（per-set 权重） | CASTLE finalize、H1437 排查 | 18–36 run，单卡 3–5 周 |
| E12b | HCC1937 | R1.m1 | P2 | 下载 | — | 约 1 周集群 |
| E13 | learning curve 和模型大小 | R1.7 | P1 | 无代码 | E1 | 7 run |
| E14a | COLO829T mixture 异常排查 | R1.5、R2.4 | P1（先做） | 新脚本 | — | 1–2 天 |
| E14b | Illumina / ONT mixture 模型 | R2.4、R3.10 | P1 | 无代码 | E14a | 4 run |
| E14c | HCC1395 in-silico dilution（本地数据） | R1.5、R2.4、R3.10 | P1 | 无代码（新比对 + build） | — | 约 1 周集群 |
| E14d | SEQC2 湿实验 titration | R2.4 | P2 | 查 metadata | — | 同 E14c 一档 |
| E15 | 低 AF 配置（AF 0.03） | R2.4 | P1 | 无代码（选项已有） | E14a | build + 训练，tensor 多 3–5 倍 |
| E16 | 多 seed、block bootstrap、阈值敏感性 | R3.9，所有主表 | P0 | 新脚本 + GPU | E3、E8、E12 | 18–30 run |
| E17a | INDEL 分层表 | R3.10 | P1 | 新脚本 | E2、E8 | CPU |
| E17b | 跨个体 error panel（长读长 INDEL） | R3.10、R1.6 | P2 | 新脚本 | E2 | 1–2 天 |
| E17c | 在 validation 上重选 INDEL PoN 规则 | R2.6 | P1 | 新脚本 | E2 | 每个模型一次 validation 推理 |
| E18 | SBS96 和签名分解 | R3.5、R2.m6 | P1 | 新脚本 | E1、E8、E9 | CPU |
| E19 | 按附近 germline 密度分层 | R2.5 | P1 | 新脚本 | E3 | CPU |
| E20 | normal reads 在吸收位点上的走向（COLO829BL ONT） | R2.3 | P2 | 新脚本 | — | 1 天 |
| E21 | 训练样本计数表、VAF 直方图 | R1.7、R1.8 | P1 | 新脚本 | 最终 run | CPU |
| E22 | 文字和图 | R2.2、R2.m、R3.7、R3.8、R3.m1、R3.m3、R3.m4 | P0（随时） | — | — | — |
| E23 | 去掉 donor 的 HapMap 图 | R1.5 | P2 | 新图（照 mgi_non_HG005） | — | 约 1 周 |
| E24 | locus / bubble-aware tensor | R1.1、R1.m2 | P2 | 改 pipeline | — | 2–3 周 |
| E25 | branch candidate 模式 | R1.4 | P2 | 改 pipeline | E6 | 2–3 周 |

**算力粗估。**

- P0 的 GPU 部分约 80 个 run：E1 6、E2 6、E3 约 30、E12 先做 SNV 18、E16 约 20。
- 单卡算约 6–8 周，多卡并行可以压到 2–3 周。
- E12 的 d46 联合训练最慢，单个 run 1–2 天。

**本地已有、设计里用到的资源（2026-10-12 查过）。**

- 18 套参考臂 tensor（linear / d9 / d46 × 6），9 套 COLO829T mixture（d9），HCC1395 三套和 H1437 PacBio（d46）。
- 组装：HG008-T v3.2 hap1 / hap2，HG008-N v6.2 hap1 / hap2，HG008-N dipcall（含 hap → GRCh38 PAF），COLO829BL dipcall。
- reads：
  - COLO829BL ONT 的 GRCh38 BAM；
  - SEQC2 的 HCC1395 WGS_FD_T_1/2/3、HCC1395BL SRR7890827、湿实验 titration 一个 run；
  - HG008-T 各平台；HG008-N 的 reads 不在本地。
- 图：
  - d9 / d46 的 `.snarls`；
  - d9 的 GRCh38 / CHM13 node 坐标库（CHM13 全部 25 条染色体）；
  - v1.1 全图 raw VCF 和 vcfbub wave VCF；
  - 去掉 HG005 的 d9 图（流程可照搬）。
- 没有、要下载或要装的：
  - GIAB v3.5 stratifications；
  - HPRC v2.1 全图 VCF / GBZ（E7、E7b）；
  - CHM13 FASTA（也可以用 `vg paths` 从图里导出 CHM13 path，不必下载）；
  - read simulator（wgsim / pbsim3）；
  - HG008-N reads（E20）；
  - HCC1937 reads（E12b）；
  - COLO829T 组装（可选）。

---

## 7. 要用户决定的事

1. **主图**：E3 之后是否把主图换成 HPRC v2.1 d46。换的话，COLO829T mixture 要在 d46 上重新比对 9 套，或者 mixture 继续用 d9
   作为单独的实验。
2. **主评估**：主结果是否用跨样本（E12，留出细胞系），同一样本内的模型只作为上界。
3. **HapMap**：直接去掉，还是做去掉 donor 的图（E23）后放进 Supplement。
4. **下载**：GIAB stratifications（E10）、HPRC v2.1 全图 VCF（E7）、HG008-N reads（E20）、HCC1937（E12b）、COLO829T 组装（E9d，
   可选）。
5. **HCC1395 dilution**（E14c）：不用下载，但每档要重新比对，约一周集群时间。
6. **新功能**：E24（bubble-aware tensor）、E25（branch candidate）、E9f（off-ref 组装标签）做不做。建议先做 E6 和 E9a–e，再决定。
7. **rule (a)**：是否也用到 HG008T 和 COLO829T（现在只用于短读长 ensemble 的 truth）。
8. **E3 的判定**：如果 PoN 后打平，主张怎么改（见 R3.1 的判定表）。

---

## 8. 简版要同步改的地方

写详版时发现简版有几处要改。两处事实错误已在简版里直接改掉，其余是可以用新数字加强的地方：

- **R3.1**（已改）：原来写的是 "aligned with the same vg giraffe version"，不对。linear 用 vg 1.66，d46 用 vg 1.77，d9 [待核实]。
- **R1.6**（已改）：原来写的是 "Population filtering alone caps INDEL recall at 0.30 on HG008T chr1"。vcfeval 报告给的是 0.350（chr1 全染色体，
  635 个 truth 里 413 个被标掉），已改成 0.35。
- **R1.4**：要换成 10-12 的新数字：
  - 被吸收的 truth 99 % SNV、95–97 % INDEL 被现行 PoN 标掉；
  - 图特有的损失：HG008T INDEL truth 的 1.3 %，COLO829T 的 0.6 %。
- **R1.2 / R3.2**：candidate 表可以换成六套 × 三个参考臂的完整版（本文件 R1.2）。
- **R2.4 / R1.5**：HCC1395 dilution 的数据是本地的，简版已改（2026-10-12）。
