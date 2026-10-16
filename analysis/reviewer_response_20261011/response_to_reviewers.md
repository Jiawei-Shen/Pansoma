# Pansoma: point-by-point response to the Nature Biotechnology reviews (NBT-RA72852A)

Draft of 2026-10-11, written against the current pipeline: `indexed_gam_pipeline_v4` (candidates, tensors, labels) and
`machine_learning/pansoma_net_v2` (model, VCFs, evaluation).

How to read this draft:

- **Response** is the English text for the response letter.
- **现有 pipeline 能否解决** 是内部评估（中文，发出前删除）：
  - **能** = 已有结果，或现有代码直接能做。
  - **部分** = 已有部分证据，还要补实验。
  - **不能** = 需要新实验或新功能。
- `[待补: …]` 是要等新实验的数字或要核实的事实，发出前必须填上或删掉。
- 数字都是 chr1（test chromosome）、benchmark BED 内，另有说明的除外。
- SNV 的 chr1 数字来自 2026-10-01 的六个 SNV 模型。这些模型用的是 10-03 重建之前的 tensor（那一版已删除），最终数字要等重训（见文末 P0-1）。

## 总表

| # | 意见（简述） | 能否解决 | 主要缺口 |
|---|---|---|---|
| R1.1 | 不是 pangenome-native | 部分 | linear/d46 模型；可选 bubble plane |
| R1.2 | 图的组成和大小 | 部分 | d46 模型；单倍型数梯度实验 |
| R1.3 | 困难区域、ONT 歧义 | 部分 | GIAB stratification；三参考 alignment 指标 |
| R1.4 | 图里已有的 somatic allele 被删 | 部分 | 分析已完成；branch candidate 功能 |
| R1.5 | HapMap 不合适、有泄漏 | 部分 | 去掉 HapMap 或去 donor 的图；dilution 模型；HCC1395 dilution（数据已有） |
| R1.6 | recall 低 | 部分 | SNV 重训；比较器同协议；INDEL 现行 recipe |
| R1.7 | 训练数据规模和多样性 | 部分 | 跨样本 / leave-one-line-out |
| R1.8 | 标签过程不清楚 | 能 | 整理成 Methods 和计数表 |
| R1.m1 | 更多数据集 | 部分 | HCC1395 等的模型和评估；CASTLE truth 是 DeepSomatic 的；HCC1937 没做 |
| R1.m2 | 窗口截断在节点边界 | 能 | 图示 |
| R2.1 | HapMap | 部分 | 同 R1.5 |
| R2.2 | truth 怎么来的 | 能 | 写 Methods |
| R2.3 | recall 原因、正常样本 | 能 | 已完成 |
| R2.4 | 低 VAF | 部分 | mixture 异常；低 AF build |
| R2.5 | precision 从哪来 | 部分 | linear 模型 FP 分层 |
| R2.6 | PoN 问题 | 能 | 核实样本是否在资源里 |
| R2.7 | alignment 统计 | 部分 | 逐 read 三参考统计 |
| R2.m | 图和文字 | 能 | 改图 |
| R3.1 | 参考层面的 ablation | 部分 | 三臂模型 |
| R3.2 | candidate 减少的组成 | 部分 | FP 按 allele 是否在图里分层 |
| R3.3 | graph 坐标 calls 的验证 | 部分 | 组装验证、跨平台 |
| R3.4 | 指标没包括全部输出 | 不能 | off-ref 人群过滤（新功能）；CHM13 投影（先查现有坐标库） |
| R3.5 | SBS96 和 truth 谱 | 能 | 小脚本 |
| R3.6 | perfect 不等于 correct | 部分 | MAPQ、个人组装、CHM13 臂 |
| R3.7 | FN 结论太肯定 | 能 | 已有正交证据 |
| R3.8 | 范围 | 能 | 改文字 |
| R3.9 | 稳定性和不确定性 | 部分 | 多 seed、block bootstrap |
| R3.10 | 低 VAF 和 INDEL 分层 | 部分 | 最终 INDEL 模型后出表 |
| R3.m1–m4 | 参数、比较器、图组成、术语 | m1/m3/m4 能，m2 不能 | 比较器重跑 |

---

## Summary of the main changes (for the letter)

We rebuilt the Pansoma pipeline and revised the manuscript around four questions the reviewers raised.

1. **What does the graph contribute?** We added a reference-level ablation in which the same reads are processed by
   identical code on a GRCh38-only graph, HPRC v1.1 and HPRC v2.1.
2. **What does the graph remove?** We quantified which truth, germline and background candidates each reference
   produces. We also traced every somatic allele hidden by a graph path to the patient's own assemblies.
3. **Is the complete output correct?** We now report metrics with graph-coordinate calls included, and we validated
   those calls orthogonally.
4. **Does the model generalise?** We added four breast and lung cancer lines, cross-sample tests, purity dilutions,
   repeated training and confidence intervals.

The HapMap mixture was removed from the main benchmark. The claims are now restricted to tumor-only calling of SNVs,
with INDEL performance reported separately.

（内部：第 2 条已完成；第 1、3、4 条和 HapMap 的处理都要等对应实验做完才能这样写，见文末 P0 / P1。）

---

## Reviewer 1

### R1.1 "Pangenome-native" is not supported

**Comment.** The graph is used only for alignment. Graph topology, bubbles and haplotype paths are discarded, so the
model is a standard pileup CNN.

**Response.** We agree that "fully exploit the richness of pangenome data" overstated the design, and we removed it.
The Methods now state where graph information enters:

1. *Candidate definition.* Candidates are edits relative to the nodes each read traverses. Alleles already encoded as
   population haplotypes therefore produce no candidate. On identical reads, 76–82% of matched-normal germline truth
   alleles give a candidate on a GRCh38-only graph, against 8–12% on HPRC v1.1.
2. *Path-following tensors.* Each tensor row follows its read's own walk across node boundaries. Channel 5 holds the
   graph base on that walk. Channel 6 holds the number of HPRC haplotypes through each column's node.
3. *Inference uses channel 6.* Setting the site haplotype count of true somatic tensors to 38 of 90 lowers their
   median somatic probability from 0.89 to 0.10. Raising all counts to 90 more than triples germline false calls
   (2,129 to 7,505). After population filtering, a model without channel 6 performs similarly, because the filter
   carries overlapping information.
4. *Output.* Calls are written in graph coordinates and projected to GRCh38 only where a unique projection exists.

We now call Pansoma pangenome-aware and do not claim explicit modelling of bubbles or long-range haplotypes. The
reference-level ablation (R3.1) tests the graph's contribution directly.
[待补: counterfactual and channel-6 ablation repeated on the final models]

**现有 pipeline 能否解决：部分。**

- 已有：linear / d9 tensor 对比（`analysis/graph_ablation_tensors_20261005`）；ch6 counterfactual 和 nopc ablation（旧
  build 的模型）。
- 需要：
  - 在新 build 上用最终 recipe 重跑 counterfactual 和 nopc ablation。
  - linear / d46 模型（R3.1）。
  - 可选新功能：bubble plane（每个 site 所在 bubble 各 branch 的 read 数和 path count），回应 "competing path
    support"。要改 builder 和 tensor 格式；不做就在文字上承认 CNN 不建模拓扑。

### R1.2 Composition and size of the pangenome

**Comment.** How do the number and diversity of haplotypes change candidate generation and germline identification?
Could a larger graph such as HPRC v2 absorb more somatic mutations?

**Response.** We added HPRC v2.1 as a third reference arm. All datasets were realigned and rebuilt with identical code
on three references:

- a GRCh38-only graph (vg autoindex);
- HPRC v1.1 (90 haplotypes, nodes kept when on at least 9);
- HPRC v2.1 (464 haplotypes including European samples absent from v1.1, nodes kept when on at least 46).

The larger graph did not absorb more somatic alleles. On HG008T PacBio, 72.6% of truth alleles reach a tensor on v2.1
against 71.5% on v1.1, while germline candidates fall from 10.6% to 9.4%. Both graphs keep alleles carried by about
10% of haplotypes, and at that filter the haplotype count changed candidate generation little.
[待补: the other datasets; trained models on each arm]

We also characterised absorption itself (R1.4). On HPRC v1.1, 2,478 HG008T and 497 COLO829T truth alleles align
without edit to an existing graph path. 84–87% of these are alleles of other HPRC individuals, mostly one-unit
homopolymer or STR changes. Absorption rises with repeat length.
[待补: absorption per dataset on v2.1; haplotype-count series]

**现有 pipeline 能否解决：部分。**

- 已有：d46 tensor 六套全部完成，另有 HCC1395 三套；HG008T PacBio 的 recall 对比。
- 需要：
  - d46 六套的模型和 chr1 结果（同 recipe）。
  - 每套 d46 的 absorbed truth 统计。`tensor_recall` 的脚本可以直接跑；graph_absorbed 的 s1–s14 可选。
  - 新实验（建议做，回应 "number and diversity"）：只做 HG008T PacBio chr1，用 HPRC 单倍型子集（8 / 16 / 32 / 64 /
    88 条，同一频率过滤）建图，画 absorption 和 germline 抑制随单倍型数的曲线。需要建图、giraffe 索引、重新比对。

### R1.3 Difficult regions and ONT ambiguity

**Comment.** Aggregate alignment statistics do not show gains where pangenomes should help. More ONT graph-coordinate
calls may reflect ambiguity rather than resolution.

**Response.** We added region-stratified analyses.

1. *Alignment.* 4.9–8.7% of node visits by MAPQ > 5 alignments lie on non-reference nodes, and 94.3–98.1% of those
   visits are edit-free. Within each node-length bin, branch and GRCh38 nodes are about equally often perfect, so the
   gain is not an artefact of short branch nodes. ONT reads also reach segmental duplications that short reads miss:
   among missed HG008T truth SNVs, one on ONT has no MAPQ > 10 read, against 48–55 on the other platforms.
   [待补: the same statistics on the GRCh38-only graph, and MAPQ distributions, per GIAB stratification]
2. *Calling.* We report precision and recall per GIAB stratification for HG008T and per SMaHT easy, difficult and
   extreme region for COLO829T, for each reference arm and each comparator. [待补]
3. *ONT graph-coordinate calls.* On chr1, PASS SNV calls without a GRCh38 position number 105 (PacBio), 595 (Illumina)
   and 855 (ONT) for HG008T, and 2,515 for COLO829T ONT. We test resolution against ambiguity in two ways. First, the
   platforms share one graph, so calls can be compared across platforms in node coordinates. Second, supporting reads
   are checked against the tumor and matched-normal assemblies (R3.3). [待补]

**现有 pipeline 能否解决：部分。**

- 已有：d9 上的 branch node 统计（`analysis/branch_node_alignments_20261002`）；COLO829T 的 SMaHT 区域；miss 分析里的
  S1 区域。
- 需要：
  - GIAB GRCh38 stratification BED。服务器上没有，要下载，用户决定。
  - linear / d46 的同样统计：discovery 的 `node_stats.json` 已有，s1 脚本直接跑。MAPQ 分布要扫 GAM，用现有 reader
    写个小脚本。
  - 按区域分层评估 calls：vcfeval 加 BED 交集，小改动。
  - ONT unplaced calls 的跨平台一致性和组装验证（R3.3）。

### R1.4 Candidate generation removes somatic variants already encoded in the graph

**Comment.** If tumor reads align perfectly to an alternative path, no candidate is made. Supplementary Fig. 3
misattributes these false negatives to GRCh38 reference bias.

**Response.** The reviewer is correct, and we have withdrawn the reference-bias interpretation. Using the
matched-normal assemblies (HG008-N v6.2; COLO829BL), the HG008-T tumor assembly and HPRC haplotype walks, we classified
every truth allele whose reads align without edit to an existing path.

- **HG008T:** 2,478 alleles (201 SNVs, 2,277 INDELs). Only 110 (4.4%) are the patient's germline allele. 2,085 (84.1%)
  are absent from both normal haplotypes but carried by other HPRC individuals, with a median frequency of 0.17.
- **COLO829T:** 497 alleles. 34 (6.8%) are germline and 430 (86.5%) belong to other individuals.

Most are one-unit homopolymer or STR changes.

Linear tumor-only pipelines lose the same alleles. Among those absent from the normal, our standard SNV population
filter removes every SNV (8.4% of other truth SNVs), and 93.5% of INDELs match a population allele at AF ≥ 0.01
(45.3% of other truth INDELs). In PacBio fiberseq data,
DeepSomatic-TO labels 55% of COLO829T chr1 truth INDELs as germline. We now present this as a limitation of
tumor-only calling that the graph makes explicit.
[待补: branch-candidate mode and the recall it recovers]

**现有 pipeline 能否解决：部分（分析已完成，功能没有）。**

- 已有：`analysis/graph_absorbed_somatic_20261001`、`graph_absorbed_somatic_colo829t_20261005`、PoN 对比，以及
  DeepSomatic 的对比（COLO829T fiberseq chr1，v1.9）。
- 需要：
  - 新功能（建议做）：branch candidate。在 bubble 处把 "tumor reads 走低频 branch（path count ≤ X）" 也作为候选，
    由 ch6 和模型判断。评估能找回 2,478 / 497 里的多少，FP 代价多少。
  - 比较器在这些位点上的表现（需要同一 BAM 上的比较器 VCF）。
  - HCC1395 等新细胞系的同样统计（脚本可复用）。

### R1.5 HapMap mixture as a somatic benchmark, and data leakage

**Comment.** The HapMap positives are germline alleles that the graph is designed to absorb. Four of the six donors
are in HPRC v1.1.

**Response.** We agree on both points. Four donors (HG00438, HG02257, HG02486, HG02622) are in HPRC v1.1 and also in
v2.1. We removed the HapMap mixture from the main benchmark.

Low-VAF performance is now evaluated with in-silico tumor/normal dilutions of COLO829T: purity 25, 50 and 75% on all
three platforms. These contain true acquired mutations and were built with the same pipeline. The same design on
SEQC2 HCC1395 (Illumina, 5–50% purity, SEQC2 truth) tests low-VAF calling on a second tumor.
[待补: HCC1395 dilution results]

The Supplement keeps a HapMap analysis on a graph with the four donors removed, and states the limitation explicitly.
[待补]

**现有 pipeline 能否解决：部分。**

- 已有：COLO829T P25 / P50 / P75 九套 tensor。fiberseq 联合模型 PoN 后 F1：P100 / P75 / P50 / P25 为
  0.785 / 0.699 / 0.663 / 0.722（实测 purity 0.65 / 0.38 / 0.17），P75、P50 比 P25 差，异常没查。
- 需要：
  - 先查异常，再训练 Illumina 和 ONT 的混合模型。
  - leave-donor-out 图：GBZ 去掉 4 个 donor，重做频率过滤、索引、比对和 tensor。工作量大，也可直接去掉 HapMap。
    用户决定。
  - HCC1395 dilution，数据都在 `/scratch/jshen/data/SEQC2/`，不用下载：tumor = WGS_FD_T_1 + _2（~108x），normal =
    HCC1395BL `SRR7890827`（~58x）。按 ~60x 混 P5 / P10 / P20 / P50，normal 够用；只作稀释，不作 matched normal。
    每档要 giraffe + sort + tensor，用户决定。
  - SEQC2 湿实验 titration `Tumor_Purity_Titration/SPP_GT_1-1_1.sra`（~148x）也在盘上，只有一个 run，purity 待查。

### R1.6 Low sensitivity

**Comment.** F1 gains come from precision. SNV recall is below ClairS-TO and DeepSomatic-TO, and INDEL recall is 4–22%.

**Response.** We changed the calling rule to favour sensitivity and now report recall at every stage. The PASS
threshold is chosen on validation blocks to reach 90% validation recall, and population filtering raises precision
afterwards. On chr1 after filtering, SNV recall is 0.83–0.84 on HG008T (precision 0.45–0.50) and 0.72–0.73 on
COLO829T (precision 0.72–0.88). [待补: retrained models; comparators under the identical protocol]

We decomposed every missed truth allele. 91–95% of truth SNVs reach a tensor. On HG008T the misses are graph SNP alleles
and repeat absorption. On COLO829T they are low-VAF subclonal mutations and regions without confident mapping.

The INDEL limitation is real, and we state it. 43–68% of truth INDELs reach a tensor, and 81–92% of those without
any candidate are graph-encoded (R1.4). Moreover, 72–84% of all truth INDELs are exact population alleles. Population filtering alone
therefore caps INDEL recall at 0.35 on HG008T chr1, for any tumor-only caller that uses it. We now show full
precision–recall curves for every caller so readers can choose their operating point.

**现有 pipeline 能否解决：部分。**

- 已有：SNV 召回分解；R0.9 规则；PoN 前后的上限（`analysis/tensor_recall_20260930`）。
- 需要：
  - 六个 SNV 模型在新 build 上重训。现有 SNV 模型用的是已删除的旧 build；linear / d46 用的是新代码，比较必须在
    同一 build 上。
  - 比较器同协议重跑（R3.m2）。
  - INDEL 按 10-11 定的 recipe（no filter、sqrt、lr 1e-4）跑全部六套，目前还没有完全一致的 run。
  - 想提高 recall 的话：branch candidate（R1.4）和低 AF build（R2.4）。

### R1.7 Training scale and tumor diversity

**Comment.** A ~200 M-parameter model trained mainly on COLO829 (40–50 k variants, one UV-driven melanoma) may be
under-trained and tumor-specific.

**Response.** *Model size.* The SNV network now has 50.4 M parameters and the INDEL network 8.7 M. On our data the
199.6 M model overfitted, with validation loss rising from 0.51 to 2.16 by epoch 12, and gave no better test accuracy.

*Example counts.* We report somatic, germline and non-variant examples per platform and variant type. For example,
HG008T SNV sets hold 8.8–10.3 k somatic, 0.31–0.37 M germline and 0.07–1.9 M non-variant tensors on chr1–22. COLO829T
SNV sets hold about 38 k somatic tensors.

*Diversity.* We added two breast (HCC1395, HCC1954) and two lung (H1437, H2009) cancer lines with published
multi-technology truth sets (SEQC2, CASTLE). They add about 323 k somatic truth alleles, against about 46 k for COLO829T.
Together with melanoma and pancreatic cancer, this covers distinct mutational processes. HCC1395 is complete on all
three platforms. [待补: other lines; signature per line]

*Tumor specificity.* We test it directly with leave-one-cell-line-out training and with cross-sample tests, for example
training on COLO829T and testing on HG008T. [待补]

**现有 pipeline 能否解决：部分。**

- 已有：模型大小；标签计数；HCC1395 三平台 tensor（d46，rule (a)）。另外三株在比对和排序（2026-10-11）。
- 需要：
  - 注意：现在 12 个模型都是同一样本内 chr2–22 训练、chr1 测试。审稿人会认为这是 within-tumor 评估，所以必须加
    跨样本和 leave-one-line-out。
  - 新细胞系都在 d46 上，HG008T / COLO829T 的主结果在 d9。多肿瘤联合训练要统一用 d46（HG008T / COLO829T 的
    d46 六套已完成）。

### R1.8 Construction of the training set

**Comment.** How are positives projected and labelled, including graph-path alleles and normalisation? What makes a
negative "determinable", and are germline variants and artefacts represented?

**Response.** We now describe labelling in full (Methods, Supplementary Note). The rule set is versioned by checksum.

*Positives.* Each truth allele is expanded to all equivalent placements in repeats (up to 5,000). Each placement is
converted to node coordinates along the GRCh38 path. Insertions at node junctions get a key on both nodes, and
deletions over consecutive nodes get a multi-node key. Candidates are left-normalised on the node's forward strand by
the same convention, so both sides are normalised. A tensor is somatic when its representative allele equals a truth
key. Partial matches (allele or read-haplotype overlap > 0.45) are recorded separately. SNV partial matches are
excluded from training because they taught low-AF artefacts: excluding them raised chr1 F1 from 0.623 to 0.665.
Alleles absorbed by graph paths cannot become positives (R1.4).

*Negatives.* Labels have three classes. Germline labels come from assembly-based matched-normal truth (HG008-N and
COLO829BL dipcall), giving 0.31–0.42 M germline SNV tensors per dataset. Non-variant labels mark tensors on a unique
GRCh38 node inside the confident region (somatic ∩ germline benchmark BED) without any truth allele. These are
sequencing and alignment artefacts, including errors next to true variants. "Determinable" means exactly this. Other
tensors are ignored in training, because the germline truth cannot confirm them: those outside the region, without a
GRCh38 position, or on branch nodes without a truth match. Purity varies only through the COLO829T dilutions
(25–75%).

**现有 pipeline 能否解决：能（已有，要整理成 Methods）。**

- 规则全部在 `analysis/truth_labels.md` 和 `tensor_postprocessing/README.md`；标签计数在 v4 README §5。
- 需要：
  - 把每个 run 的训练集正负样本数（chr2–22 去掉 validation block）从 checkpoint 导出成 Supplementary Table。
  - rule (a)（重复区里没有 truth 的 INDEL 记 −1）从 HCC1395 开始用，HG008T / COLO829T 要不要统一，用户决定。

### R1.m1 Additional benchmarks

**Comment.** Add SEQC2 HCC1395 and the other DeepSomatic cell lines.

**Response.** We added HCC1395 (SEQC2 truth; PacBio HiFi, ONT R10 and Illumina) and three CASTLE lines (HCC1954, H1437,
H2009). 32–37 k of the 41,072 HCC1395 truth alleles reach a tensor on each platform. [待补: results]

DeepSomatic was trained on these lines, so comparisons there use only the chromosomes held out from its training.
[待补: 核实 DeepSomatic 的 held-out 染色体] The HCC1954, H1437 and H2009 truth sets are DeepSomatic's own training
truth. We therefore use these three lines for training diversity and leave-one-line-out tests, and compare callers on
HCC1395, whose SEQC2 truth is independent of DeepSomatic.

**现有 pipeline 能否解决：部分。**

- HCC1395 三平台 tensor 已完成（d46）；HCC1954 / H1437 / H2009 正在比对和排序。
- 需要：模型和评估。HCC1937 没包括；它的 truth 已在 `CASTLE/benchmarks/` 里，要覆盖 DeepSomatic 全部五株的话需要加
  reads 和 tensor。
- 注意：HCC1954 / H1437 / H2009 的 truth 是 DeepSomatic 的训练 truth（`CASTLE/README.txt`），在这三株上和
  DeepSomatic 比较对它有利。HCC1395 的 SEQC2 truth 是独立的，但 DeepSomatic 也用 HCC1395 训练过（`benchmarks/`
  里有它的 HCC1395 训练 truth），所以 HCC1395 上的比较也只能在它的 held-out 染色体上做。
- SEQC2 / CASTLE 的 INDEL truth 来自短读长 caller，少了大部分 homopolymer / STR INDEL（SEQC2 HP≥7 15.5 %，HG008T
  74.9 %），所以用了 rule (a)。HCC1395 的 INDEL 数字和 HG008T 不能直接比，文中要说明。

### R1.m2 Tensor windows truncated at node boundaries

**Comment.** Windows stop at node boundaries, so short nodes lose context and path-spanning evidence is discarded.

**Response.** This described the first version. In the revised tensors the 101-column window is centred on the site.
Each row follows its own read's walk across node boundaries, on reference and branch nodes alike. The graph base
(channel 5) and haplotype count (channel 6) come per column from the node the read traverses there. Single-base nodes
therefore keep full flanking context, and node transitions inside the window are visible to the network. True somatic
tensors have a median of 0 transitions within ±10 columns, germline tensors 2 (AUC 0.72). Rows are ordered by the
similarity of their edits and node transitions, so reads taking the same alternative path sit together. Competing path
support is visible through the rows' walks, although not as an explicit bubble model (R1.1).

**现有 pipeline 能否解决：能（已实现）。**

- 需要：一张跨节点 tensor 的示意图（`scripts/visualize_tensor.py`）。
- 可选：bubble plane（R1.1）。

---

## Reviewer 2

### R2.1 Why include the HapMap benchmark

**Response.** We agree. Its positives are germline alleles of donors who are in the HPRC graphs. We removed it from
the main benchmark and replaced it with true-somatic purity dilutions (R1.5).

**现有 pipeline 能否解决：部分。** 同 R1.5。

### R2.2 How the truth sets were generated

**Comment.** Describe each truth set, and support the claim that truth sets contain germline variants and artefacts.

**Response.** We added a description of each truth set (Supplementary Note).

- *HG008T:* the GIAB HG008-T somatic small-variant draft benchmark v0.2. It was derived from multi-technology tumor data
  against the curated diploid HG008-N assembly. Confident regions are its BED intersected with the normal's dipcall
  regions.
- *COLO829T:* the SMaHT validated SNV and INDEL union, with Illumina and PacBio VAFs, over SMaHT easy, difficult and
  extreme regions.
- *HCC1395 and CASTLE lines:* the SEQC2 and CASTLE truth sets.

The germline statement is now supported directly. 209 HG008T somatic truth alleles equal an HG008-N germline allele,
from LOH or conflicting benchmarks. Of the 110 graph-absorbed HG008T alleles present in the normal, 109 equal the
germline allele of the other normal haplotype. At 33 COLO829T truth INDELs, the COLO829BL assembly carries the ALT over the whole repeat. We describe these
as benchmark ambiguities, not errors.

**现有 pipeline 能否解决：能。**

- 数字来自 `analysis/truth_labels.md` §8.1 第 5 点和两个 graph_absorbed 分析。
- 需要：写 Methods。SEQC2 / CASTLE 按原文献写。

### R2.3 Explaining the lower recall; normal sample at missed sites

**Comment.** "Variant node not ranked #1" and "first aligned node mismatch" need explanation. Check whether normal
reads map to the alternative paths at those sites.

**Response.** Those categories came from an earlier pipeline. We replaced them with read-level classes for every missed
truth allele on six datasets:

- no MAPQ > 10 read;
- fewer than three ALT reads;
- the ALT is an existing graph allele;
- repeat absorption or respelling;
- an indel longer than 50 bp.

For example, of 398 HG008T PacBio SNVs without a candidate, 65 are graph SNP alleles and 239 are repeat absorptions. Of
its 2,976 INDELs without a candidate, 92% are fully or partly graph-encoded.

We checked the normal at those sites with haplotype-resolved matched-normal assemblies, which are more direct than
normal-read alignments. Only 4.4% (HG008T) and 6.8% (COLO829T) of graph-absorbed truth alleles are present in the
normal; the rest are alleles of other HPRC individuals (R1.4). The HG008-T tumor assembly carries the ALT at 95% of
graph-encoded HG008T INDEL truths. These are real somatic events hidden by the graph, and we now say so.

**现有 pipeline 能否解决：能（已完成）。**

- 可选：如果审稿人坚持要 "normal reads"，可以把 HG008-N 的 reads 比到 d9 图上再看。需要 HG008-N FASTQ 和 giraffe，
  还没做。

### R2.4 Low-VAF detection and the low-VAF model

**Comment.** The 10% threshold limits low-VAF detection. How does the low-VAF model differ, where does the main model
stop working, and when should each be used?

**Response.** Candidates now require AF ≥ 0.08 and at least three ALT reads, for SNVs and INDELs alike. We report
sensitivity by truth VAF. On COLO829T, this threshold removes 74–82% of truth SNVs with VAF < 0.06, while 89–99.5% of
those with VAF ≥ 0.10 reach a tensor. [待补: recall per VAF bin after calling]

For low purity, we trained on COLO829T dilutions (25, 50 and 75% purity) together with the undiluted tumor, and we
report results per purity. [待补: per-purity F1 on three platforms] The HCC1395 dilutions (R1.5) test the same models
on an unseen tumor.

The earlier HapMap-trained low-VAF model is withdrawn with the HapMap benchmark. A low-AF configuration, built at
AF ≥ 0.03 and trained on the dilutions, replaces it. [待补] We give explicit usage guidance: the default model for
expected somatic VAF ≥ 0.1, and the low-AF model below that. [待补: the VAF at which the default model degrades]

**现有 pipeline 能否解决：部分。**

- 已有：VAF 分析（`tensor_recall` 的 per-truth 表有 VAF_Ill / VAF_PB）；九套 mixture tensor；fiberseq 联合模型（见
  R1.5，P75 / P50 的异常没查）。
- 需要：
  - 查 mixture 异常；训练 Illumina / ONT 混合模型。
  - 低 AF build：`--snv-min-af 0.03`，需要新 tensor 集，候选会多很多。
  - 按 VAF bin 的召回曲线：per-truth 表加 predictions，小脚本。

### R2.5 What drives the higher precision

**Comment.** Is the precision gain explained by reduced germline contamination? Does the graph help when several
germline variants are nearby?

**Response.** We decomposed false positives. Before population filtering, 99% of HG008T Illumina SNV false calls are
matched-normal germline variants. The graph removes most germline alleles before classification: on identical reads,
76–82% of germline truth alleles produce a candidate on a GRCh38-only graph, against 8–12% on HPRC v1.1. Population
filtering then removes 95.5–96.6% of the remaining germline SNV tensors. The false positives left are mostly rare
germline alleles (gnomAD AF < 0.001).
[待补: the same decomposition for the GRCh38-only arm, split by whether the germline allele is in the graph (R3.2)]

For clustered germline variation, we stratify precision and recall by the number of matched-normal germline variants
within ±50 bp, on the HPRC and GRCh38-only arms. [待补]

**现有 pipeline 能否解决：部分。**

- 已有：FP 组成（`analysis/v1_vs_v2_20260929`、`tensor_recall` 的 PoN 表）。
- 需要：
  - linear 模型的 calls，然后按 "allele 是否在图里" 分层 FP。
  - 按附近 germline 数分层：dipcall 加一个小脚本。

### R2.6 Panels of normals

**Comment.** Some panels contain the benchmark samples. Were the >20% germline variants missing from the graph removed
by the panels? Report HG008T without panels.

**Response.**

1. *Leakage.* HG008 and COLO829 are not among the individuals of these resources. [待补: 核实 gnomAD、dbSNP、1000G、
   CoLoRSdb 的样本表] The HapMap donors are, which is another reason we removed that benchmark (R1.5).
2. *Germline variants not in the graph.* They remain candidates and are labelled germline in training. Population
   filtering removes 95.5–96.6% of germline SNV tensors on all six datasets.
3. *HG008T without filtering.* Chr1 SNV F1 is 0.18–0.20 without filtering and 0.59–0.62 with it. Precision is
   0.10–0.11 against 0.45–0.50, and recall 0.86–0.87 against 0.83–0.84.

Filtering cannot reach calls without a GRCh38 position (R3.4). We renamed these "population and background filtering
resources". INDELs now use a separate documented rule: an exact allele with AF ≥ 0.01 in homopolymers of at least
7 bp, otherwise ≥ 0.001.

**现有 pipeline 能否解决：能。**

- 需要：
  - 核实 HG008 / COLO829 是否在 CoLoRSdb、gnomAD 的样本里。
  - INDEL 规则是在 chr1（test）上选的，要在 validation 或独立样本（HCC1395）上重新确认。HCC1395 的 SEQC2 truth
    里 HP≥7 的 INDEL 很少，只能确认 0.001 那一档，0.01 那一档要用 HG008T / COLO829T 的 validation 染色体。

### R2.7 Alignment statistics overstated

**Comment.** Report read counts or proportions for perfect alignments and soft clipping.

**Response.** We replaced the figure with counts and proportions. For example, HG008T PacBio has 3.93 billion node
visits by MAPQ > 5 alignments. 8.0% of them are on branch nodes, and 98.0% of those are edit-free, against 92.4% on
GRCh38 nodes. [待补: per-read perfect, mismatch, indel and soft-clip counts for the same reads on GRCh38-only, HPRC
v1.1 and v2.1] We also toned down the statement (R3.6).

**现有 pipeline 能否解决：部分。**

- node 层面已有。同一批 reads 在三个参考上的逐 read 统计要扫三个 GAM，用现有 reader 写小脚本。

### R2.m Minor comments 1–8

**Response.**

1. Figure 2d and 2f show truth-set variants and have moved to the Supplement.
2. Figures are now cited in order (2e before 2f; Figure 4 panels follow the text, minor 7).
3. The COLO829T panel of Supplementary Fig. 2 is labelled.
4. Point labels in Supplementary Fig. 1 were moved off the points.
5. "Truth overlapping" was replaced by defined categories, for example a call within 10 bp of a truth allele with a
   different allele. [待补: 核实第 13 页原文的意思]
6. The truth-set spectrum is shown beside each caller (R3.5).
7. See item 2.
8. A *candidate* is one allele (node, start, REF, ALT). A *candidate site* is all alleles at one node position and
   variant class, encoded together in one tensor.

**现有 pipeline 能否解决：能（改图和文字）。**

---

## Reviewer 3

### R3.1 Reference-level ablation

**Comment.** Run the identical Pansoma pipeline on GRCh38 (or a single-path graph) and on the HPRC graph to isolate the
pangenome's contribution.

**Response.** We performed this ablation. Reads of all six datasets were aligned with vg Giraffe (1.66 for the GRCh38-only graph, 1.77 for HPRC v2.1; [待核实: v1.1]) to
three references:

- a GRCh38-only graph built by vg autoindex from the same FASTA;
- HPRC v1.1;
- HPRC v2.1.

Tensors were built with identical code, thresholds and labelling; only the graph and its haplotype count differ.
Candidate-level results are in R3.2.
[待补: models trained with the identical recipe on each arm; chr1 precision, recall and F1 before and after population
filtering, with bootstrap intervals (R3.9)]

On the GRCh38-only graph the haplotype-count channel is constant. That comparison therefore includes both the
alignment and candidate effect and the population channel. An HPRC arm without channel 6 separates the two. [待补]

**现有 pipeline 能否解决：部分（tensor 全部完成，模型没训）。**

- 需要：每个参考 × 平台 × 样本 × SNV / INDEL 用同一 recipe 训练。linear 6 + 6、d46 6 + 6，d9 在新 build 上重训
  6 + 6，约 36 个模型。
- 一个 SNV 50.4 M 模型约 5 h（一块 H100）；INDEL 小模型可以在 MIG 上跑。
- linear 图是 32 bp 节点，tensor 数约为 d9 的两倍，训练更久。

### R3.2 Composition of the removed candidates

**Comment.** Compare linear and graph candidate generation for somatic truth, matched-normal germline variants and
background candidates. Stratify false-positive reductions by whether the germline allele is in the graph.

**Response.** We did this on identical reads. The table gives the fraction of truth alleles with a candidate tensor:

| | GRCh38-only | HPRC v1.1 |
|---|---:|---:|
| HG008T somatic SNVs | 98.3–98.6% | 94.7–95.2% |
| COLO829T somatic SNVs | 90.3–90.9% | 89.8–90.1% |
| matched-normal germline variants | 75.7–82.2% | 7.7–11.7% |

On the graph, total tensors fall by 43–69%, and non-variant SNV candidates in the confident region by 11–61%. The graph
thus removes about nine in ten germline candidates and loses 0.5–4 points of somatic SNV sensitivity.

Somatic INDELs are the exception: HG008T insertions drop from 53–80% to 36–40%, because repeat-length somatic changes
coincide with population alleles in the graph (R1.4).
[待补: final false positives split by whether the germline allele is in the HPRC graph, linking the 79% / 73% overlap
to the precision gain]

**现有 pipeline 能否解决：部分。**

- 三类 candidate 的对比已完成（`analysis/graph_ablation_tensors_20261005`）。COLO829T ONT linear 和后来完成的 d46 几套
  还没进表，重跑 `collect.py` 即可。
- FP 按 "allele 是否在图里" 分层需要 linear 模型的 calls，再查图里的 allele（d9 VCF 已有）。

### R3.3 Validation of graph-coordinate variants

**Comment.** Validate a sample of graph-exclusive calls with matched-normal alignments, phased assemblies or
cross-platform data.

**Response.** We added three orthogonal checks. [待补: results]

1. *Assemblies.* For each PASS call without a GRCh38 position, the consensus of its ALT reads is aligned to the HG008-T
   tumor assembly and both HG008-N haplotypes (for COLO829T, the COLO829BL haplotypes). A call is confirmed somatic when
   the ALT is in the tumor and absent from both normal haplotypes.
2. *Cross-platform.* The three platforms share one graph, so the same call is compared directly in node coordinates.
3. *A gold set.* In the patient's own assembly frame, 450 HG008T somatic events lie on non-reference nodes, and 233 have
   an exact tensor. We report how many are called.

The same assembly test, applied to off-reference candidates without a truth match, found only 1.2–1.6% tumor-specific.
Unconfirmed graph-coordinate calls therefore count as false positives in the headline metrics (R3.4).

**现有 pipeline 能否解决：部分。**

- 已有：`assembly_offref_20260929`（450 个 B 类事件，233 个 exact key）；`offref_negatives_20260929`（A1 read consensus
  比对到组装的脚本）。
- 需要：对现有 PASS unplaced calls 跑组装验证和跨平台一致性（同一 d9 节点坐标，小脚本）。
- 风险：unplaced PASS 很多（COLO829T ONT 2,515），大多数可能是假阳性，需要先改进 off-ref 部分（R3.4）。

### R3.4 Metrics exclude graph-coordinate calls; "about 2%"; CHM13

**Comment.** Precision, recall and F1 cover only GRCh38-projectable calls. The graph-coordinate fraction exceeds 10–20%
on ONT, contradicting "about 2%". How many are representable on CHM13?

**Response.** We now report two metrics for every run. The first covers GRCh38-evaluable calls. The second covers the
complete output, where every graph-coordinate PASS call counts as a false positive unless the assembly test confirms
it (R3.3). [待补: graph-coordinate fraction and complete-output precision after the revised off-reference calling]

The "about 2%" statement held for one dataset only and has been corrected.

To separate GRCh38 limitations from graph-specific variation, we project each graph-coordinate call onto the CHM13 path.
CHM13 is part of both HPRC graphs. We report the fraction representable on CHM13. [待补]

**现有 pipeline 能否解决：不能（关键部分需新功能）。**

- 现状严重。PoN 后 unplaced PASS 占 SNV PASS 的比例：HG008T PacBio 8%、Illumina 32%、ONT 41%；COLO829T Illumina 12%、
  fiberseq 13%、ONT 51%。全部算 FP 时 precision 是 HG008T 0.28–0.46、COLO829T ONT 0.36。
- 原因：off-ref rescue 把 ch6 设成 100，去掉了人群先验，而 PoN 到不了 graph 坐标。
- 需要：
  1. graph 原生的人群过滤：检查 off-ref call 的 ALT 序列是否被 HPRC 单倍型 walk 拼出。`hprc_offref_filter` 因
     sqlite BFS 太慢停了，可改用 GFA W 行（walks.ndjson）的局部索引。
  2. 或者给 off-ref calls 一个单独的阈值，在 validation 上选。
  3. CHM13 投影：已有一个 d9 的 node → GRCh38 / CHM13 坐标库
     （`AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13.path_index.sqlite`，2026-06 用旧流程建）。它的 GRCh38
     部分缺 chr3–9，CHM13 部分是否完整要先查。不完整的话，从 GFA W 行取 CHM13 walk 建表（类似 ref-path-scan 的小工具）。
     d46 要另建。

### R3.5 Mutational signatures with the truth set

**Comment.** Show the truth spectrum beside each caller and provide SBS96 profiles, including graph-coordinate calls.

**Response.** Figure 4 now shows the truth-set spectrum beside Pansoma, DeepSomatic-TO and ClairS-TO. The SBS96
profiles of every call set, including the graph-coordinate calls, are in the Supplement. [待补]

**现有 pipeline 能否解决：能。**

- 一个小脚本：VCF 加 FASTA 取三核苷酸上下文；graph 坐标的 call 用节点序列取上下文。分解用 SigProfilerAssignment。
- 要等最终 calls。

### R3.6 More perfect alignment is not more correct alignment

**Comment.** Compare GRCh38, T2T-CHM13 and HPRC. Use MAPQ, phased-haplotype or assembly truth, and difficult regions.

**Response.** We agree and revised the wording. We added three analyses. [待补]

1. Per-read MAPQ and edit distributions for the same reads on the GRCh38-only, HPRC v1.1 and v2.1 graphs.
2. Correctness against the patient's own genome: HG008-T reads are also aligned to the curated HG008-N assembly, so each
   read's graph placement can be checked against its placement on the personal genome.
3. Stratification by GIAB difficult regions (R1.3).

A CHM13-only arm for one dataset separates reference completeness from multi-haplotype modelling. [待补]

**现有 pipeline 能否解决：部分。**

- 已有：HG008T reads 比到 HG008-N hap1 的 BAM（`HG008_GIAB/HG008T_align_to_HG008N_hap1`、
  `HG008T_PacBio_vs_HG008N_hap1`，2026-06），只有 hap1，可以作为 placement 的参照；GAM 统计用小脚本。
- 需要：
  - CHM13-only 图（vg autoindex）和 giraffe，建议只做 HG008T PacBio。
  - GIAB stratification（要下载）。

### R3.7 Caution on reclassifying benchmark false negatives

**Comment.** Graph alignment alone cannot overturn truth annotations. Use orthogonal evidence or soften the claims.

**Response.** We agree and revised the text. We now give orthogonal evidence (R1.4). Haplotype-resolved matched-normal
assemblies show that only 4.4% (HG008T) and 6.8% (COLO829T) of these alleles are germline. The HG008-T tumor assembly
carries 95% of the graph-encoded HG008T INDELs. Most are therefore genuine somatic events hidden by population alleles
in the graph, not benchmark errors. Supplementary Fig. 3 is reframed accordingly.

**现有 pipeline 能否解决：能（已完成）。**

### R3.8 Scope: tumor-only versus somatic mosaicism

**Response.** Pansoma is designed and validated for tumor-only calling in cancer samples. We now state this scope
throughout, in the title, abstract and Discussion. Somatic mosaicism in normal tissue is listed as future work, because
it needs lower allele fractions and its own training data.

**现有 pipeline 能否解决：能（改文字）。** 不建议做 SMaHT 组织：VAF 范围不同，需要新模型。

### R3.9 Model stability and statistical uncertainty

**Comment.** Report multi-seed retraining, block-bootstrap intervals for differences to comparators, and sensitivity
to the fixed 0.5 threshold.

**Response.** We now report four things. [待补]

1. Three to five retrainings with different seeds for each main model.
2. 95% intervals from a block bootstrap over 1-Mb blocks of the test chromosome, for each metric and for the difference
   between Pansoma and each comparator.
3. Full precision–recall curves, and metrics at three operating points: the validation recall-0.9 threshold, the best
   validation F1 threshold and argmax.
4. Cross-sample tests (R1.7).

The fixed 0.5 threshold is no longer used. The PASS threshold is chosen on validation blocks, never on the test
chromosome.

**现有 pipeline 能否解决：部分。**

- 已有：train 有 `--seed`；vcfeval 已输出曲线和三种规则。
- 需要：
  - 多 seed 训练（GPU 时间）。
  - block bootstrap 脚本：rtg 输出的 TP / FP / FN 按位置分块重抽样。
  - 可考虑第二条测试染色体。
  - 可选：SEQC2 的技术重复。`WGS_FD_T_3` 是同一株、另一个文库（盘上已有，~56x），可以看同一模型在重复测序上的
    calls 一致性。

### R3.10 Low-VAF and INDEL claims

**Comment.** Qualify the HapMap low-VAF results, separate SNV from INDEL claims, and stratify INDELs by class, length,
context and platform.

**Response.** *Low VAF:* HapMap is out of the main text, replaced by COLO829T dilutions and VAF-binned sensitivity
(R2.4).

*INDEL:* SNV and INDEL claims are now separated throughout. INDEL performance is stratified by insertion or deletion,
length, sequence context (homopolymer ≥ 7 bp, STR, other) and platform. [待补: final table]

The main finding is already clear. 74% of HG008T chr1 truth INDELs lie in homopolymers of at least 5 bp, and most are
exact population alleles. Pansoma's false INDEL calls are germline homopolymer INDELs at the same allele fractions.
INDEL accuracy is therefore limited by this overlap more than by platform error.

**现有 pipeline 能否解决：部分。**

- 已有：`analysis/indel_repeat_context_predictions_20261006`、`tmp/indel_hp_strata_20261010`。
- 需要：最终 INDEL 模型（10-11 recipe）在六套上的结果出来后，再出分层表。
- HCC1395 / CASTLE 的 INDEL 分层要单独列：它们的 truth 少了 HP / STR INDEL（rule (a)），重复区的 recall 和
  precision 在这些株上没有意义。

### R3.m1 Parameter discrepancies between Methods and Supplement

**Response.** One parameter table now serves both. The values are:

- reads with MAPQ ≤ 10 are dropped;
- at least 3 ALT reads with base quality ≥ 10;
- AF ≥ 0.08 for SNVs and INDELs;
- INDELs up to 50 bp;
- at most 800 reads per node;
- tensors of 200 rows × 101 columns;
- SNV model 50.4 M parameters, 20 epochs;
- INDEL model 8.7 M parameters, 12 epochs.

**现有 pipeline 能否解决：能。** 参数都在 manifest 和 checkpoint 里。

### R3.m2 Comparator setup

**Response.** A new table lists, for ClairS-TO and DeepSomatic-TO, the software and aligner versions, model checkpoints,
command lines, reference and filtering resources, and PASS criteria. [待补] All callers are scored with the same truth,
BED, chromosome and population-filtering rule, and also with each caller's own default filters.

**现有 pipeline 能否解决：不能（需重跑）。**

- 现有只有 COLO829T fiberseq chr1 的 DeepSomatic v1.9 TO 和早期 ClairS 输出。
- 要在全部数据集（包括 HCC1395 等）上重跑并记录版本和参数。

### R3.m3 Composition of the pangenome reference

**Response.** Methods now specify the graphs.

- *HPRC v1.1:* the Minigraph-Cactus GRCh38 graph, keeping nodes on at least 9 of 90 haplotypes (44 individuals × 2,
  CHM13, GRCh38). The individuals are African 23, American 16, East Asian 4, South Asian 1 and European 0.
- *HPRC v2.1:* nodes on at least 46 of 464 haplotypes (462 HPRC haplotypes, GRCh38, CHM13).
  [待补: individuals per population]

Structural variants stay in the graph, but candidates are limited to indels of at most 50 bp. HG008, COLO829 and the
SEQC2 / CASTLE lines are in neither graph; the four HapMap donors are in both.

**现有 pipeline 能否解决：能。** v2.1 的个体和人群可以从 `tmp/hprc_v21_samples` 统计。

### R3.m4 Terminology and language

**Response.** We standardised the terms: Pansoma (not Pansoma-TO), "Supplementary Figure", and the node pileup (NPU)
expansion. "Panel of normals" became "population and background filtering resources". The quoted sentences were
corrected, and the manuscript was language-edited.

**现有 pipeline 能否解决：能（改文字）。**

---

## 需要新增的实验和功能（按优先级）

**P0：挡住多条回复，先做。**

1. 六个 d9 SNV 模型在新 build 上重训。现有 SNV 模型用的是已删除的旧 build。（R1.6、R2.6、R3.1）
2. 三个参考臂 × 六套 × SNV / INDEL 同 recipe 训练和 chr1 评估。（R3.1、R1.2、R3.2、R2.5）
3. 比较器 ClairS-TO / DeepSomatic-TO 同协议重跑，记录版本。（R3.m2、R1.6、R3.9、R3.5）
4. graph 坐标 calls：组装验证、跨平台一致性、off-ref 人群过滤（新功能）、CHM13 投影（先查现有 d9 坐标库）。
   （R3.3、R3.4、R1.3）

**P1：主要回复的补充。**

5. 跨样本和 leave-one-line-out，统一在 d46 上。（R1.7、R1.m1）
6. 多 seed 和 block bootstrap。（R3.9）
7. VAF 和 purity：查 mixture 异常、Illumina / ONT 混合模型、低 AF build；HCC1395 Illumina in-silico dilution
   （SEQC2 tumor + HCC1395BL reads，都在盘上）。（R1.5、R2.4、R3.10）
8. 区域分层：GIAB stratification、三参考的 alignment 指标和 MAPQ。（R1.3、R3.6、R2.7）

**P2：可选的新功能和实验。**

9. branch candidate 模式。（R1.4、R1.1）
10. bubble plane。（R1.1、R1.m2）
11. HPRC 单倍型子集图。（R1.2）
12. 去掉 donor 的 HapMap 图；SEQC2 湿实验 titration（盘上只有 SPP_GT_1-1 一个 run，purity 待查）；HCC1937
    （truth 已有，reads 要下载）。（R1.5、R1.m1）
13. CHM13-only 图（HG008T PacBio）。（R3.6）
14. SBS96、INDEL 分层表、图和文字修改。（R3.5、R3.10、R2.m、R3.m）

**要用户决定的事：**

- HapMap 是直接去掉，还是做去 donor 的图。
- 是否下载 GIAB stratification。
- 是否做 HCC1395 dilution（不用下载），以及 HCC1937 要不要下载 reads。
- 是否做 branch candidate 和 bubble plane。
- rule (a) 是否用到 HG008T / COLO829T。
- INDEL 人群过滤规则要在非 chr1 的数据上重新确认。
