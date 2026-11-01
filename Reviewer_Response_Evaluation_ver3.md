# Pansoma Reviewer Response — Evaluation ver3

更新日期：2026-11-01。评估对象：当前 `Reviewer_Letter_Response_ver3.md` 及六张附图。评估文件名保持不变，response 和原图未作修改。

**总体评价：4/5。40 条实质性意见中，38 条的回复内容达到 OK；2 条仍有必要修改：R1.8 的模型选择/标签表述，以及 R2.code 的可运行版本和复现材料交付。** 后者是多条回复共享的交付问题，统一记录一次；38 条 OK 不等于相应附件均已取得或独立核验。

从面向 Nature Biotechnology 的方法学审稿角度，我对本版的判断是：**主要科学质疑已得到有证据、有边界的回应，适合在针对性修正和材料交付后继续审理；不需要把整套研究重新推倒，也尚不能仅凭这份回信建议直接接收。** 4 分评的是回复充分性，不是文章创新性、期刊录用概率或软件发布完成度的另一个代称。

本次改变评分来自重新校准“当前主张需要什么证据”，不是 ver3 新增了实验，也不是以作者是否愿意补实验决定分数。尤其不能把严格机制分解、VAF 分箱的潜在细节、comparator 模型选择的未证实疑问都变成新的必做门槛；也不能仅凭回信写得具体就假定软件和附件已交付。

## 评估标准、范围与核验边界

采用与此前相同的尺度：**5 分**表示该条回复完整、直接、证据充分；**4–4.5 分＝OK**，足以回应当前质疑，可以有合理且明确的限制；**3–3.5 分＝部分解决**，仍有需要直接处理的解释、分析或交付缺口；**1–2.5 分＝重要缺口**；**0 分＝缺失**。总体分按关键科学问题判断，不取逐项算术平均。一个术语澄清拿 5 分，不表示它比核心方法学问题更重要。

本评估参考 Nature 公开审稿指南中对结果有效性、结论支持度、原创性与重要性的区分，但下面的 5 分制和具体判定均为本次独立评估，并非期刊官方评分。研究还可更强，与当前结论缺少必要证据，应分开判断。[Nature：How to Write a Report](https://www.nature.com/nature/for-referees/how-to-write-a-report)；[Nature：Peer Review](https://www.nature.com/nature/editorial-policies/peer-review)。

- **科学依据来自当前 ver3。** 旧评估只用于检查要求和评分是否一致，不将旧对话或仓库中的其他实验当作 ver3 已提供的证据。
- **不把承认局限本身当成充分回应。** 只有实际证据与收窄后的结论相匹配时，限制说明才可以关闭问题；影响当前结论的错误仍须修正。
- **不把增强研究的实验自动升级为必做实验。** 比如共享 genomic split 会增强 reference attribution，但当前系统比较加配对 channel 消融已能支持较窄主张；如果恢复更强主张，证据要求也随之提高。
- **不重复扣共同交付缺口。** R3.m1 的参数说明、R3.m2 的比较规则可以内容充分，同时其配置/命令附件仍待交付。共同缺口集中在 R2.code，未被视为已经完成。
- **逐条评分只列必要修改。** 每条保留 reviewer 原文、本条需要回应什么、ver3 实际回应、评分与理由。只有未达到 OK 的条目出现“还需做”。“需要回应什么”是判定标准，不是给 OK 项追加任务。

完整评估已覆盖当前 Markdown 的 1,120 行及六张 PNG，并核对主要跨表关系、HapMap 分母、54 组 per-seed P/R/F1、96 个 alignment count/percentage 和 chr2 confirmation 均值。本次更新重新审读争议条目的原问题、完整回复和补充定义；原文件哈希与此前全面核查时一致，因此保留仍然成立的数字/图像核查结果。未重新训练、执行原始 VCF 评价或重跑 bootstrap。

当前包没有完整 revised manuscript、完整 Supplementary Data、逐运行配置及可访问 release。对“全文已修改”“完整表已提供”的说法，只能评价其具体回应是否合理，不能认证实际附件或全文已经一致。这是本次核验范围的限制，不能反推实验没有做。Reviewer 2 的八条 minor comments 分开评分，共 40 条实质意见；三位 reviewer 的开头总评、Reviewer 3 的最终总评和 R1.code 的肯定性意见另列，共 45 个回复单元，后五项不重复计数。

## 当前证据支持的研究主张

| 核心主张 | 当前证据与专业判断 | 尚不能据此推出的结论 |
|---|---|---|
| Pangenome 信息在 alignment 之后仍有预测价值 | 固定比对和候选的 channel 6 配对重训中，12 组过滤前 AP 均下降、10 组过滤后 F1 下降；支持 haplotype-count 信息有用。 | 未证明显式学习 bubble topology、长程 haplotype，或 path-spanning context 本身的独立性能增益。 |
| HPRC 作为参考提高当前 Pansoma 系统的 SNV 表现 | 同流程的三参考比较、三 seeds、候选组成和困难区域结果方向一致；比仅做跨工具 benchmark 更能回答 reference 贡献。 | validation split 未共享，不能把所有增益归为纯 reference 因果效应，更不能精确分解各机制贡献比例。 |
| 方法有潜在的 tumor-only 应用价值 | 两个主要肿瘤各三平台，并增加 HCC1395 及其他 cell-line 分析；多种数据下的 SNV precision/F1 结果有实质价值。 | 两肿瘤×三平台不是六个独立肿瘤；HCC1395 的 within-tumor benchmark 不是独立 tumor-held-out 泛化。held-out SNV 多为单次探索，COLO829T UV transfer 失败限制了普遍适用性。 |
| INDEL 和低 VAF 是当前方法的明确限制 | absorption、完整 truth recall、分层 INDEL 与 dilution/transfer 数据直接揭示局限；当前回复不再包装成普遍优势。 | 不能宣称全面优于 comparator 的 INDEL sensitivity、可靠支持低于约 0.08 的 somatic calling，或已适用于正常组织 mosaicism。 |
| Off-reference calls 已被审查 | normal/assembly/cross-platform 证据与完整输出的保守评价揭示多数额外 calls 缺少充分 somatic 支持，rescue 默认关闭与证据一致。 | 未证明广泛发现新的 reference-external somatic mutations；germline-like spectrum、CHM13 表示性或跨平台共现都不能单独证明 somatic。 |

这是一项范围清楚的技术贡献，可信之处在于用对照和不利结果界定收益与代价。是否达到 Nature Biotechnology 所需的重要性与广泛用途，还取决于完整论文的核心定位和证据整体；回复逐条通过不能替代这一编辑判断。本次不因收窄主张而自动给予高分，也不因方法不是通用解法而自动判为失败。

## 评分速览

| Reviewer | 实质意见数 | 内容达到 OK | 仍需直接修改 |
|---|---:|---:|---:|
| R1 | 10 | 9 | 1 |
| R2 | 16 | 15 | 1 |
| R3 | 14 | 14 | 0 |

| 条目 | 主题 | 评分 / 5 | 判断 |
|---|---|---:|---|
| [R1.1](#r1-1) | 比对之后是否真正利用 pangenome 信息 | 4 | OK |
| [R1.2](#r1-2) | 图的组成与规模 | 4.5 | OK |
| [R1.3](#r1-3) | 困难区域及 ONT mapping ambiguity | 4.5 | OK |
| [R1.4](#r1-4) | 图中已有 somatic allele 被候选生成漏掉 | 4 | OK |
| [R1.5](#r1-5) | HapMap 生物学适用性与 donor overlap | 4.5 | OK |
| [R1.6](#r1-6) | 低 recall，尤其 INDEL sensitivity | 4.5 | OK |
| [R1.7](#r1-7) | 模型容量、训练样本量与跨肿瘤泛化 | 4 | OK |
| [R1.8](#r1-8) | 标签规则、unknown negatives 与模型选择 | 3.5 | 需补充 |
| [R1.m1](#r1-m1) | 新增 HCC1395 与其他 cancer cell lines | 4.5 | OK |
| [R1.m2](#r1-m2) | 节点边界截断及 path-spanning tensor | 4.5 | OK |
| [R2.1](#r2-1) | 为何保留 HapMap benchmark | 4.5 | OK |
| [R2.2](#r2-2) | truth sets 的来源及 benchmark error 主张 | 4 | OK |
| [R2.3](#r2-3) | recall 损失原因及 matched-normal 证据 | 4.5 | OK |
| [R2.4](#r2-4) | 低 VAF 模型、敏感性与适用范围 | 4 | OK |
| [R2.5](#r2-5) | precision 提高的机制及 clustered germline variants | 4.5 | OK |
| [R2.6](#r2-6) | population filtering 的泄漏与独立贡献 | 4 | OK |
| [R2.7](#r2-7) | alignment statistics 是否夸大 | 4.5 | OK |
| [R2.m1](#r2-m1) | 原 Figure 2d/f 的数据来源 | 4.5 | OK |
| [R2.m2](#r2-m2) | Figure 2 panel 引用顺序 | 4.5 | OK |
| [R2.m3](#r2-m3) | Supplementary Figure 2 的 COLO829T 标识 | 4.5 | OK |
| [R2.m4](#r2-m4) | purity/VAF 标签遮挡数据点 | 4.5 | OK |
| [R2.m5](#r2-m5) | truth overlapping 的含义 | 5 | OK |
| [R2.m6](#r2-m6) | Figure 4 加入 truth spectrum | 4.5 | OK |
| [R2.m7](#r2-m7) | Figure 4 panel 顺序 | 4.5 | OK |
| [R2.m8](#r2-m8) | candidate sites 与 candidates | 5 | OK |
| [R2.code](#r2-code) | 依赖、安装失败与可运行发布 | 3 | 需补充 |
| [R3.1](#r3-1) | 三参考比较是否隔离 pangenome 的贡献 | 4 | OK |
| [R3.2](#r3-2) | 候选减少是否选择性保留 somatic | 4 | OK |
| [R3.3](#r3-3) | graph-coordinate calls 的独立验证 | 4.5 | OK |
| [R3.4](#r3-4) | 完整输出指标、graph-only 比例及 CHM13 | 4.5 | OK |
| [R3.5](#r3-5) | truth 与 graph-coordinate SBS96 谱 | 4.5 | OK |
| [R3.6](#r3-6) | perfect alignment 与 correct placement | 4.5 | OK |
| [R3.7](#r3-7) | 不能凭图比对改判 false negatives | 4.5 | OK |
| [R3.8](#r3-8) | 肿瘤应用与正常组织 mosaicism 的范围 | 5 | OK |
| [R3.9](#r3-9) | 训练稳定性、bootstrap 与 operating point | 4 | OK |
| [R3.10](#r3-10) | 低 VAF 和 INDEL 的结论边界 | 4 | OK |
| [R3.m1](#r3-m1) | 参数统一与逐模型可复现材料 | 4 | OK |
| [R3.m2](#r3-m2) | comparator 配置、过滤和模型选择公平性 | 4 | OK |
| [R3.m3](#r3-m3) | 图版本、组成、过滤、SV 与 donor 关系 | 4.5 | OK |
| [R3.m4](#r3-m4) | 术语和语言一致性 | 4.5 | OK |

## 本次评分校准的依据

| 条目 | 原 evaluation → 本次 | 改变或保留判断的原因 |
|---|---|---|
| R3.1 | 3.5 → 4 | 实际三参考比较与 channel 消融有证据价值；已明确 system-level 的边界。此前将严格共享 split 仍列为强制门槛，未充分尊重当前较窄主张。 |
| R2.4 | 3.5 → 4 | 原问题的四个子项均回答；expected VAF 定义存在核验边界，但没有证据证明计算错误，不应由潜在疑问推出必须重算。 |
| R3.m1 | 3.5 → 4 | 具体参数与五个 features 已讲清。未见逐模型附件是真实的交付问题，集中列入 R2.code，不重复算作方法解释失败。 |
| R3.m2 | 3.5 → 4 | PR/filter 规则已明确；ilmn/ilmn_ssrs 差异公开但尚无不公平证据。actual commands 仍纳入共同交付核对。 |
| R1.8 | 3.5 → 3.5 | fixed-checkpoint 与整个选择流程仍被混淆，final INDEL 的 other/排除规则需说明；将原先默认要求定向敏感性实验的做法，改为按作者最终保留的主张决定实验是否必要。 |
| R2.code | 3 → 3 | 占位符、将来时与外部可运行性仍未闭合；不会因其他条目改善而自动提高该分。 |

Reviewer 3 的开头和最终总评同步由 3.5 调至 4，不参与 40 条计数。其余逐条评分维持原判断。此前要求中的 R1.1 配对条件、R1.5 分母、R3.2 proxy 定义、R3.9 表间差异及 chr2 汇总、R3.m2 PR 规则已经得到回应；R1.8 部分未闭合，R3.m1 的内容已解释但附件仍待核验，R2.code 仍需实际交付。不能把已完成的解释继续列为“作者没改”。

## Reviewer 1：逐条评价

<a id="r1-general"></a>

### R1.general｜Reviewer 1 总评

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:25)

**Reviewer 原文**

> This manuscript presents Pansoma, a machine-learning framework for tumor-only somatic variant calling directly from pangenome graph alignments. Pansoma aligns sequencing reads to a human pangenome graph, summarizes the resulting alignments into node-level pileups, and uses image-like tensors derived from these pileups for somatic variant prediction. The idea of incorporating population-scale pangenome information into tumor-only somatic variant calling is potentially interesting because pangenome references may help reduce reference bias and provide additional information for distinguishing germline from somatic variation, with the pangenome effectively serving as a population-level “generic normal” prior.
>
> However, in its current form, I do not think the method demonstrates what it claims. I have serious and fundamental concerns about the current study design and how the pangenome graph is utilized. Pansoma claims to be a pangenome-native somatic variant caller but does not exploit graph-specific information beyond alignment stage. Graph topology, edges, bubbles, cross-node connectivity, and broader haplotype/path information, is all discarded during training and variant inference. As a result, the downstream model resembles a conventional pileup-based deep-learning architecture, such as DeepVariant or DeepSomatic (Pansoma is not very much different from what DeepVariant and DeepSomatic have done), rather than a method specifically designed to model the pangenome itself. The benchmark results alone do not demonstrate that Pansoma meaningfully exploits graph-native information during variant inference. The manuscript therefore overstates its main methodological novelty.
>
> Another critical concern is that Pansoma does not convincingly demonstrate the specific advantages of using a pangenome graph. A major motivation for pangenome references is that population haplotypes and alternative paths can improve read alignment in repetitive, low-mappability, and structurally complex regions where a linear reference is limited. However, the manuscript provides no region-stratified analysis showing that Pansoma actually improves somatic variant calling in these regions, nor does it provide an ablation that separates the contribution of pangenome-based alignment from that of the downstream node-pileup model.
>
> Also, several aspects of the methodological design and evaluation, including the questionable choice of the HapMap mixture as a somatic benchmark and the relatively low recall. More seriously, the candidate-generation strategy may completely miss true somatic variants already represented on alternative graph paths, preventing them from ever reaching model inference. This exposes a fundamental conceptual problem in the current design. All these require substantial methodological reconsideration and raise concerns about whether the central claims of the manuscript are sufficiently supported.
>
> See below for my specific comments:

**本条需要回应什么：** 总体质疑集中在图信息是否真正进入推断、图带来的净收益是否被隔离，以及 absorption、HapMap、低 recall 与训练设计是否削弱核心主张。总回复应准确说明研究现在能支持什么，并把各质疑引到对应证据。

**ver3 实际回应：** 用 R1.1–R1.8 对应上述问题，说明新增 channel 消融、图版本比较、absorption 计数、肿瘤扩展，并移走 HapMap 主 benchmark。

**评分：4/5 — OK。**

作为总览，这段有清楚的对应关系，没有声称所有方法局限已经消失。其作用是导航，具体科学充分性应看逐条评分；不能因为这段组织良好就推断 Reviewer 1 已接受方法。

<a id="r1-1"></a>

### R1.1｜比对之后是否真正利用 pangenome 信息

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:39)

**Reviewer 原文**

> 1. The central claim that Pansoma is “pangenome-native” is not convincingly supported by its method design. The authors state that Pansoma can “fully exploit the richness of pangenome data” (line 89), yet in practice the pangenome is used primarily during vg Giraffe alignment, while the downstream representation remains essentially a local pileup. The released implementation does not explicitly model graph edges, bubbles/superbubbles, cross-node context, or long-range haplotype structure in either the input representation or its neural network. As a result, the method largely amounts to pangenome alignment followed by a standard pileup-based neural network, similar to existing DL variant callers such as DeepVariant, DeepSomatic, and ClairS. Replacing linear genomic coordinates with graph node IDs alone is insufficient to establish a method as pangenome-native. Moreover, haplotype- and path-level information, which represents one of the key advantages of a pangenome graph, is largely discarded. The authors should therefore provide substantial additional evidence or re-design their architecture to demonstrate that Pansoma derives measurable benefit from pangenome-specific information beyond the alignment stage itself.

**本条需要回应什么：** 证明收益并非仅来自换了比对参考或节点坐标；需要在可比条件下测试图特有输入的增益，并说明模型究竟利用了哪些 path/topology 信息、哪些没有利用。

**ver3 实际回应：** 撤回“fully exploit”，改称 pangenome-aware；给出 12 组 channel 6 有/无重训配对，配对内固定 reads、candidates、labels、recipe、seed 和 GPU 数。12 组过滤前 somatic AP 均下降，10 组过滤后 F1 下降，median ΔF1 = −0.023。跨节点 context 被明确列为表征修正，没有宣称已经单独证实其增益。

**评分：4/5 — OK。**

这已经直接回答“beyond alignment”的最低核心要求：固定候选与比对，去掉 haplotype-count 后排名能力下降。v1.1 的五组是一 GPU 对一 GPU，不能把它们误判成与四 GPU control 混配；正文也披露了 COLO829T Illumina 的反例。单 seed 消融、channel 与 population filter 的重叠，以及未显式建模 bubble/长程 haplotype，使其不足以证明“全面 graph-native”，但当前收窄后的主张可以接受。去掉 channel 后仍优于 GRCh38 是系统层面的支持，不构成各机制贡献比例的精确分解。

<a id="r1-2"></a>

### R1.2｜图的组成与规模

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:62)

**Reviewer 原文**

> 2. Follow the first comment, a fundamental issue is that the authors do not provide a systematic analysis of how the composition and size of the pangenome reference affect somatic variant calling. It remains unclear how the number and diversity of reference haplotypes, haplotype walks, and alternative paths influence node-pileup construction and germline identification. This is crucial in understanding pangenome, because the graph itself determines which alleles are already represented as population variation and therefore which variants may be suppressed before reaching the classifier. This becomes important as newer references such as HPRC2 incorporate many more haplotypes. While a larger graph may improve population representation, it may also increase mapping ambiguity or cause true somatic mutations that overlap population alleles to be absorbed into existing graph paths and missed.

**本条需要回应什么：** 考察不同参考的 haplotype 数量、组成、频率过滤、mapping ambiguity、germline suppression 和 somatic retention，避免把 release 差异全部归因于图大小。

**ver3 实际回应：** 比较 GRCh38-only、HPRC v1.1 的 90 haplotypes/d9 与 v2.1 的 464/d46；明确两个 HPRC release 同时改变人数、组成和构图。报告 germline candidate 减少、过滤前后 F1、ONT 错位率及额外 INDEL 损失。

**评分：4.5/5 — OK。**

回应同时覆盖收益和代价，且没有把 v2.1 的 ONT 错位率 1.45% 对 v1.1 的 0.60% 隐去。“release comparison，不是 graph-size causal experiment”的边界正确。审稿人能够据此判断图依赖性；本条不需要为了完整回答而再强制增加一个 haplotype-number 梯度实验。

<a id="r1-3"></a>

### R1.3｜困难区域及 ONT mapping ambiguity

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:74)

**Reviewer 原文**

> 3. The manuscript fails to demonstrate one of the primary practical motivations for using a pangenome reference. Pangenome alignment can remain highly ambiguous in repetitive, paralogous, structurally complex, and low-mappability regions, especially for noisy ONT reads. This concern is also related to graph composition, since adding more haplotypes and alternative paths may introduce additional competing mappings in difficult regions. The authors report that graph alignment increases perfectly aligned reads and reduces mismatches, indels, and soft clipping relative to GRCh38. However, these aggregate statistics do not demonstrate that the improvement occurs in the regions where pangenome references are expected to provide the greatest advantage.
>
> This is important for ONT, for which the authors observe more graph-coordinate predictions and attribute this partly to higher sequencing error and more complex alignment profiles. Without region-specific benchmarking and error analysis, it remains unclear whether these graph-based signals represent improved resolution or simply additional alignment ambiguity. The current results therefore do not establish that Pansoma meaningfully addresses the challenging regions that motivate the use of a pangenome in the first place.

**本条需要回应什么：** 用区域分层的 calling 和 alignment 证据回答图在哪里有益，尤其区分 ONT 的更干净比对与更多歧义，不能仅报告全局 mismatch 数下降。

**ver3 实际回应：** 给出 GIAB difficult、low-mappability/segmental duplication、repeat 以及 SMaHT extreme 分层，报告 F1 差和 paired intervals；同时给出 MAPQ 转移、已知起源模拟 reads，以及 ONT graph-coordinate calls 的验证结果。

**评分：4.5/5 — OK。**

证据与问题匹配，包含长读 INDEL 的负面结果及 v2.1 在 extreme regions 的 MAPQ 损失。12 个 difficult-union 比较均为正、10 个区间排除零的陈述，比只报 precision 更完整。不能据这些结果声称所有困难区域都被解决，但 ver3 已明确限制。区域表的运行可追溯性作为共同交付事项，统一在 R2.code 处理。

<a id="r1-4"></a>

### R1.4｜图中已有 somatic allele 被候选生成漏掉

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:88)

**Reviewer 原文**

> 4.Candidate generation in Pansoma can delete true somatic variants that are encoded in the graph. Pansoma defines candidates only through mismatches and indels against the traversed node. If a true somatic allele is already encoded in an alternative path and tumor reads align perfectly to it, Pansoma treats it as the local graph reference and eliminates it before inference: no candidate, tensor, or call is produced. This is a fatal methodological issue. Supplementary Figure 3 exposes this failure, yet the authors incorrectly describe these false negatives as GRCh38 reference-bias artifacts. The presence of an allele on an alternative path in the pangenome does not mean that it is a germline variant. The same allele can be germline in one individual and somatic in another.

**本条需要回应什么：** 承认 population allele 在某患者身上仍可为 somatic；在完整 truth 上量化候选阶段损失、strict/partial absorption、population-filter 重叠和不可恢复部分，不可借图中的存在性改判 truth。

**ver3 实际回应：** 明确承认 no mismatch 就无 candidate，撤回 reference-bias 解读。给出完整 autosomal truth、线性图有而 v2.1 无候选的 INDEL 比例 24–42%、过滤后额外损失 3.5–7.0%，并区分严格吸收与 residual allele 不再匹配 truth 的部分吸收。说明当前方法不能恢复这些事件，SNV 与 INDEL 分开结论。

**评分：4/5 — OK。**

这是对真实算法盲区的可审查取舍说明；不再用“normal/pangenome 中存在”直接抹去 false negatives，也没有把 comparator 同样过滤掉它们当作无害证明。Reviewer 1 仍可能不接受这种取舍，但“回复充分”不要求作者在这一轮发明完整的 alternative-path caller。4 分对应充分解释与限制，不代表 INDEL sensitivity 问题已解决。

<a id="r1-5"></a>

### R1.5｜HapMap 生物学适用性与 donor overlap

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:114)

**Reviewer 原文**

> 5. The HapMap low-VAF experiment raises a fundamental concern about whether the proposed pangenome-based mechanism is being evaluated with an appropriate somatic benchmark. The manuscript explicitly states that the “somatic-like” variants in the HapMap mixture are actually germline variants from the minor HapMap individuals, which are used to low VAFs by mixing multiple cell lines. However, a central rationale of Pansoma is that germline alleles already represented in the pangenome can align cleanly to existing graph paths and therefore be suppressed as somatic candidates. The authors themselves report that a large fraction of the germline variants in the benchmark normal samples are already present in the HPRC graph, and they argue that this representation helps suppress germline‑like signals. These two assumptions are potentially contradictory: the positive variants in the HapMap benchmark are themselves population germline variants, whereas true acquired somatic mutations would generally not be expected to follow the same population representation. This raises a major concern about the biological validity of the training and evaluation strategy using pangenome graph using synthetic HapMap results. As a result, the HapMap mixture may systematically disadvantage a pangenome-based caller and does not provide a biologically appropriate setting for evaluating its ability to detect true acquired somatic mutations at low VAF.
>
> The HapMap experiment can also suffer from data leakage. HPRC v1.1 contains four of the six HapMap donors used in the mixture HG00438, HG02257, HG02486, and HG02622. The authors themselves also raise this point in a companion paper (Fu, Xin, Miao et al., "Leveraging Human Pangenome for Improved Somatic Variant Detection," bioRxiv 2026.01.04.697580). Yet these donors are not excluded when benchmarking Pansoma.

**本条需要回应什么：** 分别回答合成 germline mixture 是否能支持 somatic 结论、六个 donor 是否进入所用图、共享 allele 怎样影响分母，并对 comparator 使用相同评价条件。

**ver3 实际回应：** 说明实际用的是去掉六 donor 的 82-haplotype 专用图；HapMap 移至 Supplement，仅衡量混合物中的低比例 allele 检测。列出 graph-carried positives、conditional/full-truth recall 和相同分母的 Mutect2，详细解释 span ±1 bp mask 的连带排除。

**评分：4.5/5 — OK。**

主要质疑已闭合。SNV 分母确实可重建为 392,567 − 175,773 − 1,701 + 1 = 215,094，INDEL 为 124,499 − 25,938 − 6,482 = 92,079；最后的 +1 是脚注中未被 mask 的 N-reference SNV，不是算术遗漏。full-truth SNV recall 0.302 对 Mutect2 0.630 被直接承认，因此没有用较好看的 conditional recall 冒充总体 sensitivity。

<a id="r1-6"></a>

### R1.6｜低 recall，尤其 INDEL sensitivity

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:140)

**Reviewer 原文**

> 6. Although Pansoma achieves higher F1 scores in several benchmark settings, these gains appear to be driven primarily by improved precision rather than improved sensitivity. Pansoma does not demonstrate a convincing sensitivity advantage for tumor-only somatic variant calling, which limits its practical value for tumor-only downstream interpretation. Supplementary Figure 1 and Tables 2–3 show that Pansoma generally has lower SNV recall than ClairS-TO and DeepSomatic-TO across COLO829T and HG008T. The situation is more concerning for INDELs, where Pansoma recall is only \~4–22% in COLO829T and \~13–19% in HG008T. In my understanding, sensitivity is of great importance because there is no matched normal sample or second calling stage that can readily rescue somatic variants missed during candidate generation or classification. Missing such a large fraction of true variants is a serious limitation for tumor-only analysis, where false negatives may be permanently excluded from downstream interpretation in many discovery and clinical-prioritization settings. A higher F1 score alone does not compensate for this loss of sensitivity in many practical tumor-only applications.

**本条需要回应什么：** 对完整 truth 给出 P/R/F1 和 loss-stage funnel，区分 operating threshold 与 test-set 最优点，说明是否真正改善 sensitivity，而不是以 F1 抵消漏检。

**ver3 实际回应：** 分开 exact candidate 与不计 credit 的 partial match，报告 filtering 前后及 comparator native recall，给出 matched-precision 分析和候选/过滤上限；明确 INDEL operating recall 0.11–0.18 没有优于投稿版，过滤上限不是 tumor-only 的固有限制。

**评分：4.5/5 — OK。**

回复正面承认最不利的结果，同时交代 SNV recall 的改善与 INDEL 的不同瓶颈。matched-precision 使用 Pansoma 曲线对 comparator 固定 operating point，R3.m2 已解释其不对称性；它只能说明该工作点的取舍，不能证明 Pansoma 的整条 PR 曲线优于 comparator。当前文本没有把它作为所有场景下 sensitivity 更强的结论，故本条可接受。

<a id="r1-7"></a>

### R1.7｜模型容量、训练样本量与跨肿瘤泛化

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:167)

**Reviewer 原文**

> 7. The training data may be insufficient in scale for a model of this size, and insufficient in tumor diversity for a general-purpose tumor-only somatic variant caller. Pansoma uses a ConvNeXt backbone with around 200M parameters, but Pansoma seems to rely primarily on COLO829 for model training, which provides only on 40k-50k truth somatic variants. It is therefore unclear whether this is sufficient to train such a large neural network. The authors should report the exact parameter count and the number of positive and negative training examples per platform and variant type and provide evidence that the data volume supports the chosen capacity. Moreover, COLO829 is also a single melanoma line dominated by UV signatures, so the model may be learning tumor-specific rather than general somatic features. By comparison, DeepSomatic was trained across multiple cancer cell lines with extensive augmentation to build a broader training distribution. Authors should consider extending training to more types of tumors with a distinct mutational process and demonstrate directly that the current model is not tumor specific.

**本条需要回应什么：** 报告参数量及按类别/平台/variant type 的训练量；检验容量是否适当，并用训练时未见过的肿瘤测试不同 mutational processes 下的泛化。

**ver3 实际回应：** 网络降为 SNV 50.4 M、INDEL 8.7 M，提供容量和训练量曲线、训练计数摘要；增加四条 breast/lung cell lines，明确 W/H/T+。跨肿瘤 SNV 是单次 exploratory，INDEL 有三 seed 与预先固定的 chr2 confirmation。直接报告 COLO829T UV 转移失败及 held-out INDEL 落后 ClairS-TO。

**评分：4/5 — OK。**

没有用“把目标肿瘤加入训练后的恢复”冒充外部泛化，也没有隐去不同 cancer signature 的失败。数据支持的是训练分布覆盖条件下的泛化，尚非 general-purpose 保证；ver3 已如此表述。8.7 M 在一个容量例子里与 50.4 M 相当，并不自动证明 50.4 M 不合理。完整 per-model counts 当前仅被引用，归入 R2.code 的共同交付核对，不在本条重复开任务。

<a id="r1-8"></a>

### R1.8｜标签规则、unknown negatives 与模型选择

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:181)

**Reviewer 原文**

> 8. The construction of the training set is not also described in enough detail to be evaluated, and the choices that were made are also not justified. For positive examples, truth variants defined on GRCh38 are projected onto graph coordinates and labeled by exact allele match. However, the same biological variant can be represented differently across graph paths, especially for indels, and graph-only positions are assigned GRCh38 coordinates indirectly through adjacent 1-bp SNV nodes. It is therefore unclear how a positive label is assigned when the truth allele is already represented as a graph path, when a projection is ambiguous, or when the graph and linear representations of the same variant differ. The authors' own false-negative analysis (Supplementary Fig. 3) shows that some truth variants have their apparent allele fraction reduced to nearly zero once reads align cleanly to an alternative path, so these sites never reach candidate generation and cannot become positives at all. The authors should describe the labeling and projection procedure explicitly, including whether normalization was applied on the graph side, and justify why this procedure does not systematically distort the positive set.
>
> For negative examples, the only definition given is that determinable non-matches were labeled negative, with no explanation of what makes a candidate determinable. This matters more for a tumor-only caller than for a paired one, because the model must learn to reject germline alleles and alignment artifacts without a matched normal, and it can only do so if these classes are actually represented among the negatives. The authors should explain how the negative set was assembled, including whether germline variants were explicitly included, how artifacts were defined, and how VAF and tumor-purity conditions were sampled.

**本条需要回应什么：** 明确 graph/linear allele normalization、投影失败与 graph-absorbed positives 的处理；说明 somatic、germline、other、partial/unknown 的定义，以及 train/validation/benchmark 各自如何使用。对于不可完全确定的标签，交代假设、适用边界及其对模型选择的可能影响。

**ver3 实际回应：** 提供有序的九条规则：exact somatic → exact germline → partial 等，明确 exact germline 优先于另一 allele 的 partial somatic；列出 train/validation/benchmark 的差异。承认所有 somatic training examples 位于 GRCh38-path nodes，rule 7 unknowns 只在 validation 作 negatives，且最终 INDEL 加入部分 rule 5 样本为 other。报告所选 checkpoint 的 chr1 影响及无法重建逐 epoch validation 的限制。

**评分：3.5/5 — 部分解决，需补充。**

大部分原始质疑已经回答。exact somatic、exact germline、partial 的优先级明确；未产生 candidate 的 truth 被量化为正例覆盖损失；83–717 个 off-path negatives 的来源和三种评价阶段的规则也已解释。不能再把这些列成未完成事项。

仍有一处实质性逻辑问题：作者称把 rule 7 unknowns 作为 validation negatives “cannot move the PASS threshold”。这只在 **checkpoint 固定** 时成立。实际流程先按 validation F1 选 checkpoint，unknown negatives 可能改变不同 epoch 的 F1 排序，再间接改变被选模型的 positive scores 与阈值。rule 7 占 validation negatives 的 4.7–42.5%（SNV）和 6.6–29.4%（INDEL）；这些比例说明该假设值得明确讨论，但比例大本身并不能证明选择发生了改变或性能被高估。

作者已经承认未检验 checkpoint selection 的敏感性，这是一项有效披露。所选模型在 chr1 上去除相应 calls 后，SNV F1 变化 ≤0.0005、INDEL ≤0.05，能够说明固定模型的输出影响，不能证明换一种 validation 标记仍会选出同一模型。尤其 INDEL 的最大变化 0.05 不宜概括为普遍可忽略。另一方面，chr1 的独立 allele-based benchmark 不依赖这些 candidate labels，因此这项不确定性也不能自动推翻现有 chr1 结果。

对 final INDEL 的 rule 5，需要纠正此前评估中过强的判断：样本位于 somatic benchmark 区域、但不在 germline confident region，并不自动意味着它不能作为 somatic-negative。缺少 germline confidence 主要使其 **germline 与 artifact 的分类身份** 不确定；在所采用的 somatic truth 完整性假设下，将未匹配候选作为 non-somatic 可以是合理选择。问题是正文一处笼统说 confident region 外不参与训练，而补充说明有这一例外；同时 “other = sequencing/alignment errors” 对此例外过窄。本条应修正说明和结论边界，而非先验要求删除全部这类样本。

因此保留 **3.5 分**：有需要修正的流程表述与标签解释，但没有足够证据把它定性为已证实的训练污染、test leakage 或足以使整套 benchmark 失效的偏差。

**还需做：**

1. 将阈值表述限定为固定 checkpoint：说明 threshold 的计算只用 positives，但 checkpoint selection 使用包含 unknown negatives 的 validation F1，因而整个选择流程并非对 unknown labels 不敏感。清楚保留“尚未检验 alternative validation labeling 是否改变所选模型”这一限制。
2. 解释将 off-path unknowns 纳入 validation negatives 的目的，例如保守地惩罚缺少可确认 somatic 支持的预测；说明这是模型选择假设，不能把它等同于这些候选已被 normal 证实为阴性。现有结果支持报告该规则下的独立测试性能，不支持声称换规则后结果必然不变。
3. 在 R1.8 正文直接列出 final INDEL 的 rule 5 例外，将这里的 other 定义为未获 germline 身份确认的 non-somatic/background，而不是全部已证实的测序或比对错误；明确其阴性依据依赖 somatic benchmark 的适用区域与完整性。

**补实验的必要性边界：** 如果撤回“全流程不受影响”的说法，并接受当前方法选择下的性能结论，上述修正可以关闭本条的必要修改；不要求重训全部模型。如果仍要证明标签假设不影响 checkpoint/最终指标，则需要定向 sensitivity analysis。若各 epoch checkpoints 尚存，可重新推断 validation；若不存，才考虑代表性 INDEL 设置的重训。仅重算当前 checkpoint 的 test F1 不能支持这一更强主张。

<a id="r1-m1"></a>

### R1.m1｜新增 HCC1395 与其他 cancer cell lines

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:207)

**Reviewer 原文**

> 1. Beyond COLO829 and HG008, the authors should benchmark Pansoma on additional well-established datasets such as the SEQC2 HCC1395 cell lines, which provides a high-confidence somatic variant truth set, as well as the five additional cell lines evaluated in the DeepSomatic study.

**本条需要回应什么：** 扩展两个原有肿瘤之外的证据，尤其 HCC1395/DeepSomatic datasets；交代 benchmark truth 与训练暴露的独立性。

**ver3 实际回应：** 加入 HCC1395、HCC1954、H1437、H2009，每条三个平台。HCC1395 作为额外 benchmark；CASTLE 仅作训练与 exploratory tests，并解释 truth 构建偏向和 normal-LOH 问题。列出 SNV/INDEL F1 与训练暴露。

**评分：4.5/5 — OK。**

有实质扩展，不把同源训练/真值构建的数据包装成完全外部验证。HCC1395 ONT 落后、six-tumor INDEL 比 four-line 下降也保留在表内。审稿人列举数据集的本意是扩大验证范围；当前回应已满足这一目的，不必为凑足名单而强制新增所有可用 cell lines。

<a id="r1-m2"></a>

### R1.m2｜节点边界截断及 path-spanning tensor

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:226)

**Reviewer 原文**

> 2. Pansoma’s tensor windows are restricted to the single current node and truncated at node boundaries. Adjacent-node and path-spanning evidence has already been discarded. This contradicts the authors’ claim that Pansoma can “fully exploit the richness of pangenome data.” This design also makes it difficult to deal with short nodes or even single-base node. In fact, a large fraction of alternative nodes has the length shorter than the length of Pansoma’s tensor window. A credible graph-native method requires path-aware tensors spanning neighboring nodes rather than arbitrarily terminating the context at node boundaries and explicitly model graph edges, bubbles, alternative-path relationships, cross-node context, or competing path support.

**本条需要回应什么：** 说明短节点/单碱基节点如何保留邻接上下文，展示实际 tensor 中 competing branches/path 信息，而不只是口头更名。

**ver3 实际回应：** 8×200×101 tensor 的每行沿 read walk 跨节点；报告 57–76% SNV、92–97% INDEL 窗口超出当前节点。提供真实 chr1 SNV 的八通道图，38 ALT 与 40 REF reads 在 +21 bp 的 bubble 分别走 18/90 与 72/90 haplotypes 的分支。

**评分：4.5/5 — OK。**

已检查实际图：channel 5 的 graph bases 和 channel 6 的 path count 可见分支差异，candidate 位点、读数和图注一致。它证明了截断已从表征中消除；没有假称模型显式处理 bubble topology，也没有把该实例当成独立性能消融。

<a id="r1-code"></a>

### R1.code｜Reviewer 1 对 README 的肯定

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:237)

**Reviewer 原文**

> README is pretty comprehensive.

**本条需要回应什么：** 这是肯定性意见，没有新的科学问题；礼貌确认并指出新版使用说明即可。

**ver3 实际回应：** 感谢 reviewer，并说明 README 含安装和最小示例，指向 R2.code。

**评分：5/5 — OK。**

就这一句肯定意见而言，答复充分。此分数只评回应是否适当，不认证新版软件已经公开可运行；发布状态单独在 R2.code 评分。

## Reviewer 2：逐条评价

<a id="r2-general"></a>

### R2.general｜Reviewer 2 总评

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:245)

**Reviewer 原文**

> Shen et al. developed Pansoma, a novel pangenome-aware, machine learning-based algorithm for detecting somatic single nucleotide variants (SNVs) and short insertions/deletions (indels) from tumor-only sequencing data. As tumor sequencing becomes increasingly common in clinical practice, matched normal samples are not always available, making it challenging to distinguish somatic mutations from germline variants. Pangenome and other graph-based reference genomes improve the representation of human germline variation, facilitating more accurate identification of true somatic mutations. Pansoma is the first somatic variant caller designed to operate directly on variation graph-aligned sequencing data. The method identifies candidate mismatch sites from graph-aligned reads, stores them in a novel node pileup (NPU) format, generates feature tensors containing sequence and alignment information, and uses a deep convolutional neural network to classify candidate variants as somatic SNVs or indels.
>
> The authors evaluated Pansoma using three benchmark datasets: COLO829T (a melanoma cell line), GIEB HG008T (a pancreatic tumor model), and the SMaHT Consortium HapMap Mixture benchmark, which was designed to assess performance on low variant allele frequency (VAF) variants. Across these datasets, Pansoma outperformed two commonly used tumor-only somatic variant callers, DeepSomatic-TO and ClairS-TO, primarily through improved precision, although this was accompanied by reduced recall. Additionally, a low-VAF version of Pansoma outperformed Mutect2 on the HapMap Mixture benchmark, highlighting its potential for improved detection of low-frequency somatic variants in tumor-only sequencing data.
>
> As the field continues towards adoption of graph-based references, the field needs graph-aware tools for somatic mutation detection. Pansoma is the first to help fill this gap. The manuscript is well written and clear. The manuscript would be strengthened by more discussion on some benchmarking choices as well as an expanded investigation on the reasons for the increased precision and decreased recall of Pansoma. More specific comments included below:

**本条需要回应什么：** 回应总体肯定，并把重点放在 benchmark 选择、precision 增益来源及 recall 损失原因上。

**ver3 实际回应：** 感谢总体评价，将机制与 sensitivity 的解释交给后续逐条回复。

**评分：4.5/5 — OK。**

简洁且聚焦，没有重复 reviewer 对原稿低 VAF 优势的较强判断。总回复本身已完成导航作用。

<a id="r2-1"></a>

### R2.1｜为何保留 HapMap benchmark

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:255)

**Reviewer 原文**

> 1. The HapMap benchmark is a mixture of multiple (germline) samples at various frequencies, allowing for the detection of low-VAF events (as low as 0.5%). However, these are all germline variants, which are in theory removed by the pangenome alignments given that all of the HapMap components are included in the HPRC reference genome. The manuscript would benefit from more discussion as to why the HapMap benchmark is being included in the study as well as the limitations to this benchmarking approach for modelling somatic mutations.

**本条需要回应什么：** 说明 germline mixture 的用途、donor overlap 的真实情况及它不能代表的 somatic 生物学。

**ver3 实际回应：** 说明 donor-free 图、共享 allele 仍被吸收及双分母；明确其不能测试 acquired mutation 或 somatic/germline discrimination，不用 population filtering，因为会标记 95.7% SNV positives。主文低纯度证据改用 COLO829T dilutions。

**评分：4.5/5 — OK。**

用途与限制说得准确，而且与 R1.5 一致。不是用“已移到 Supplement”回避问题，而是清楚解释为何只能把该实验当成混合物中的 allele detection。

<a id="r2-2"></a>

### R2.2｜truth sets 的来源及 benchmark error 主张

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:263)

**Reviewer 原文**

> 2. While it is appreciated that the benchmarks used as quite standard, the manuscript would be strengthened by a brief discussion on how the ground truth datasets are generated for each benchmark, particularly, to strengthen the authors claim that the ground truth may contain germline variants and sequencing artifacts.

**本条需要回应什么：** 逐一交代 truth 的构建、验证、版本和评估区域；只有具备对应证据时才提出 germline contamination 或 artifact，不能凭 graph alignment 推翻真值。

**ver3 实际回应：** 表列 HG008 draft、COLO829、HCC1395、CASTLE、HapMap 的来源与局限，补充 matched-normal germline truth。撤回无直接证据的 artifact 主张；把 normal allele overlap 称为 genotype/LOH ambiguity，H1437 的 population/LOH enrichment 限定为 likely germline contribution。

**评分：4/5 — OK。**

总体充分而谨慎。HG008-T draft 与 assembly-based truth 的局限、CASTLE 对 DeepSomatic 的偏向和 H1437 的 chr1/全基因组区别都已说清。R1.m1 的 24.5% 与本条 gnomAD-only 27.8%涉及不同表述/统计范围，单靠两个百分比不能断定算错；本条也没有据此删除不利 truth records。其证据支持局部歧义，不支持把整套 benchmark 判为不可信，而 ver3 保留了这个边界。

<a id="r2-3"></a>

### R2.3｜recall 损失原因及 matched-normal 证据

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:281)

**Reviewer 原文**

> 3. The authors show that Pansoma out-performs both ClairS-TO and DeepSomatic-TO by increasing precision which leads to improved F1-score overall. However, Pansoma has worse recall. The decrease in recall is \~20% in COLO829T and \~10% in HG008T from DeepSomatic. The authors suggest that some of this reduced recall is due to errors in the ground truth datasets. However, the explanation of the reduced recall could be expanded.
>
> - The authors investigate the reasons behind false negatives for both HG008T and COLO829T in Supplementary Figure 2. They indicate some variants can be explained by limited allele frequency support in either the linear or the graph alignments, however, a large proportion of false negatives were labeled as “variant node not ranked #1” and “first aligned node mismatch”. These labels need more explanation.
>
> - Similarly, the authors propose that many false negatives overlap existing nodes in the pangenome reference suggesting they may correspond to germline alleles rather than somatic SNVs. To further show this, the authors could assess the alignment of the normal samples at the same sites and confirm that reads are mapping to the alternative paths preferentially.

**本条需要回应什么：** 解释原 false-negative 分类，并实际检查 normal 在被图吸收位点的 reads/alleles，验证“可能是 germline”是否成立。

**ver3 实际回应：** 撤掉难解释的旧类别，换成按 allele/read evidence 与首个丢失阶段分类；给出 normal reads 的 REF/ALT path 分析和两套 normal haplotype assemblies。v2.1 strict-absorption 子集只有约 6.7–6.8%由 normal 携带，约 92%不在两条 normal haplotypes；结论转为真正漏检的 somatic events。

**评分：4.5/5 — OK。**

真正执行了 reviewer 建议的 normal 检查，并接受其结果与原解释相反。repeat INDEL read-level evidence 不确定、COLO829T tumor assembly 与被分析 reads 不独立均已披露；v1.1/v2.1 的不同定义在 Supplement 明示，不能直接把它们的计数差当成版本效应。

<a id="r2-4"></a>

### R2.4｜低 VAF 模型、敏感性与适用范围

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:302)

**Reviewer 原文**

> 4. In a similar vein, the authors suggest the majority of false-negatives are a result of intentionally excluding low-VAF variants (as shown in Supplementary Figure 2). Thresholding candidates at 10% read support may limit the detection of true somatic SNV with low VAF as many real-world somatic-mutation analyses will deal with low-VAF given low tumor purity / high heterogeneity.
>
> - How do the authors propose to handle the detection of true low VAF somatic SNVs?
>
> - The authors also developed a specialized low VAF model, however, there is limited description of how this model differs from their main model.
>
> - Does the specialized low VAF model help recover low VAF variants? Is there a VAF threshold that the main model stops performing well on? It would be helpful if the authors conducted sensitivity analyses on decreasing VAF to demonstrate how sensitive their model is to VAF.
>
> - The manuscript would also benefit from an expanded discussion on how they propose to use each model and when.

**本条需要回应什么：** 回答低 VAF 真 somatic 如何处理、专用模型与主模型的差别、降低 VAF 时的敏感性、以及各模型的适用场景；区分 observed-AF 候选门槛和 true/expected VAF 的检测表现。

**ver3 实际回应：** 四个子问题都作了直接回答：主模型需 ≥3 ALT reads 和 observed AF≥0.08；HapMap 模型使用 0.015 floor、200-read tensors，但转到 COLO829T 后 precision 低并漏掉 clonal SNVs；因此只保留为 mixture 实验选项。用 COLO829T 三平台 dilution 报告按 expected VAF 的 recall，并处理 tumor/normal base-quality cue。

**评分：4/5 — OK。**

这条的四个原始问题均得到直接且有数据的回答。作者没有只把 0.10 改为 0.08 就宣称解决低 VAF，而是给出真实肿瘤/normal 的 dilution、分 VAF recall、HapMap 模型向 somatic 数据转移时的低 precision 和 clonal-SNV failure，并据此限制模型用途。原来最值得担心的“用 germline mixture 支撑低 VAF somatic 能力”已经撤回。

Supplementary Note R2.4 已将横轴定义为 expected VAF = truth VAF × measured purity；它是期望分箱，不能与每个位点实际观测到的混合后 VAF 完全等同。当前回复对 measured purity 的估计方法仍不够细，单凭这份文件不能独立核验局部覆盖或 copy-number 对分箱的影响。但是，这既没有证明该缩放错误，也没有推翻文中实际采用的保守结论——当前系统对低于约 0.08 的 somatic alleles 敏感性不足。此前把这种尚未证实的计算疑问直接列为必须重算的缺陷，超出了现有证据。

0.08 是 observed-AF 候选门槛，真实/expected VAF 低于此值时仍可能因采样波动通过；表中非零 recall 与此并不矛盾。本条可接受，4 分反映低 VAF 分析的解释边界以及有限的验证范围，而非认可一个精确、平台无关的生物学 detection limit。

<a id="r2-5"></a>

### R2.5｜precision 提高的机制及 clustered germline variants

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:336)

**Reviewer 原文**

> 5. The authors investigate reasons for the increase in false negatives in Pansoma; however, it would be equally interesting to better understand the decrease in false positives. What is leading to the increased precision of Pansoma? Can the majority be explained by reduced germline contamination or are their other factors?
>
> - A primary motivator for graph-aware somatic mutation detection is distinguishing germline variants from somatic mutations. However, beyond more accurate modeling of germline variation at the site itself, does graph-aware somatic mutation detection improve somatic mutation calling if there are multiple germline variants in the region?

**本条需要回应什么：** 对 false-positive 减少作来源分解，检验是否仅靠同位点 germline 抑制，还是也改善附近有多枚 germline variants 的区域。

**ver3 实际回应：** 给出 germline candidates、过滤前 germline misclassification、过滤后 excess FP 的 germline association、SMaHT extreme-region 贡献，以及 ±50 bp germline density 分层。HG008T 的 gain 随附近 germline 密度增大，COLO829T 则主要来自困难区域；明确这些类别重叠。

**评分：4.5/5 — OK。**

针对 reviewer 的两个问题都有直接量化证据。实际 R2.5 图与正文范围一致，并明确横线是观察范围而非 CI。INDEL 在 clustered germline stratum 中候选 retention 变差也被报告，因此机制没有被写成对所有 variant types 都成立。

<a id="r2-6"></a>

### R2.6｜population filtering 的泄漏与独立贡献

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:352)

**Reviewer 原文**

> 6. The authors call somatic variants using tumour-only models, supplemented with filtering (using 4 PoNs) to remove known germline variants and sequencing artefacts. This filtering approach is standard to many tumor-only variant-calling approaches. However, the authors state numerous times that alignment to the pangenome (HPRC) should improve read mapping (reduce sequence artefacts) and reduce germline-variant contamination.
>
> - Some of the PoNs used were developed using common benchmarking datasets and include variants from the tested samples [the HapMap benchmark in particular]. This may have inflated the reduced germline-variant contamination observed here as compared to real-world samples
>
> - Fig2e shows that >20% of all germline variants were not captured by the AF-HPRC graph. Were these variants successfully removed from the somatic variant calls using the PoNs?
>
> - The authors only evaluate variant calling without the PoN for COLO829. Could the authors also include the performance of the Pansoma without the PoN for HG008T as well?

**本条需要回应什么：** 核查资源是否包含测试样本；追踪图未涵盖的 germline 是否被模型/过滤去掉；补 HG008T 的无过滤评价，区分 sample overlap 与 population allele sharing。

**ver3 实际回应：** 分别做 sample-list 与 rare-allele-sharing 检查，并承认无法完全排除 gnomAD/dbSNP overlap。量化未被图携带的 germline、模型拒绝、filter 移除和残留 FP；列 HG008T 三图三平台有/无 filter 的 F1。HapMap 不用过滤；最终 INDEL 模型本身使用 population features。

**评分：4/5 — OK。**

三个明确问题均有回应。“无证据存在直接样本泄漏”被保持为有限结论，没有升级为完全排除；模型输入中的 AF 与事后 filter 也已分清。剩余 FP 主要是 rare germline 的解释有计数支持。HPRC 与 GRCh38 模型的 normalization 口径在 R3.9a 中区分，数值范围不同不应自动当作表内错误；运行映射完整性在 R2.code 作为共同交付事项统一处理。

<a id="r2-7"></a>

### R2.7｜alignment statistics 是否夸大

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:378)

**Reviewer 原文**

> 7. Lines 292-294 state that graph-based alignment increased perfect aligned reads and reduce soft-clipping. Based on the referenced figure (Supplementary Figure 4), this appears overstated. Please provide actual read counts or show proportion of reads to improve clarity.

**本条需要回应什么：** 用绝对 read 数、分母和比例替代“最多提高几倍”等容易夸大的说法；如实报告 soft clipping 与长读低基数。

**ver3 实际回应：** 提供两个肿瘤×三平台×四参考的 count/percentage 表与图，解释各参考 chr1–22 分母不同、read/mate 计数方式和 primary/supplementary 口径。撤回统一降低 clipping 的说法，区分 reference agreement 与 placement accuracy。

**评分：4.5/5 — OK。**

已逐行核对表中 96 个 count/denominator 比例，均与显示值在四舍五入范围内相符；图中趋势也一致。PacBio perfect reads 从约 0.15–0.18% 到 1.4–1.5% 的小基数被明确揭示，避免以倍数掩盖绝对效应。

<a id="r2-m1"></a>

### R2.m1｜原 Figure 2d/f 的数据来源

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 1. Are the variants shown in Figure 2 (panels d and f) generated from Pansoma or from previous linear-reference analyses? If these are based off the ground-truth variant calls, it is not clear they are added new information and could be moved to supplementary.

**本条需要回应什么：** 明确来自 Pansoma 还是已发表 truth；如果仅描述 benchmark，说明图的定位。

**ver3 实际回应：** 说明两图都来自 published truth sets，分别为 AF distribution 和 SBS96，已移至 Supplement，并修改 legends。

**评分：4.5/5 — OK。**

回答直接且合理，没有把 truth-set 描述包装成新算法结果。当前没有完整改稿可逐页核查移动位置，但回复内容已充分。

<a id="r2-m2"></a>

### R2.m2｜Figure 2 panel 引用顺序

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 2. Figure 2f is referred to before Figure 2e in the main text

**本条需要回应什么：** 消除正文先提 f 后提 e 的顺序混乱。

**ver3 实际回应：** 移走 d/f 后按首次引用顺序重新编号其余 panels。

**评分：4.5/5 — OK。**

处理方法明确。完整 manuscript 未附在此次 ver3 包中，不能声称已独立核验全稿引用顺序；这属于已具体报告的编辑修正，本条可接受。

<a id="r2-m3"></a>

### R2.m3｜Supplementary Figure 2 的 COLO829T 标识

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 3. In Supplementary Figure 2, the graph corresponding to COLO829T is not clearly labeled

**本条需要回应什么：** 让读者清楚识别 tumor、platform 和 variant type。

**ver3 实际回应：** 说明已替换成新版 missed-allele 分类图，且每张 chart 标明三项信息。

**评分：4.5/5 — OK。**

修正针对原问题，且新版图的目的与 R2.3 相符。此次包没有该完整替换图，评价依据是具体回复，并不把“已声明修改”写成已目视确认。

<a id="r2-m4"></a>

### R2.m4｜purity/VAF 标签遮挡数据点

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 4. The purity and VAF labels in Supplementary Figure 1 are often overlapping the points which make it hard to read

**本条需要回应什么：** 清除标签和 points 的重叠，保持 purity 与 VAF 信息可读。

**ver3 实际回应：** 将 purity/AF 标签移到 legend。

**评分：4.5/5 — OK。**

是足以解决该排版问题的具体动作；没有必要再要求额外分析。与 R2.m3 一样，完整新图不在六张附图中，本次不认证最终版面的执行状态。

<a id="r2-m5"></a>

### R2.m5｜truth overlapping 的含义

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 5. What do the authors mean by “truth overlapping” on page 13

**本条需要回应什么：** 明确是否指 exact truth-matching calls，并与 caller 间共同检出的 shared core 区分。

**ver3 实际回应：** 解释原意是匹配 benchmark truth 的 PASS SNVs，即 true positives；另行定义三个 caller 共同检出的 true-positive shared core。

**评分：5/5 — OK。**

定义清楚且避免把“位置相交”和“真正匹配 allele”混为一谈。

<a id="r2-m6"></a>

### R2.m6｜Figure 4 加入 truth spectrum

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 6. In Figure 4a, the authors show the cosine similarity between the signature decomposition of each tool with the signature decomposition of the benchmark truth set. It would be helpful to also include the signature composition of the truth set in the middle barplots as well.

**本条需要回应什么：** 把 truth 的突变谱与三个 caller 直接并列，使差异可逐 channel 比较。

**ver3 实际回应：** 以 SBS96 展示 truth、matched-normal germline 和三个 caller，并附 cosine similarity，见 R3.5a。

**评分：4.5/5 — OK。**

已经检查 HG008T 与 COLO829T 两张实际图：truth 行确实存在，三个平台分别显示；用原始 SBS96 直接比较是对 reviewer 科学目的的合理满足，不必为了保留旧式 signature stacked bars 再补一套重复图。

<a id="r2-m7"></a>

### R2.m7｜Figure 4 panel 顺序

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 7. The order of the figures in Figure 4 does not match how they are referenced in the main text which at times makes it confusing to read

**本条需要回应什么：** 让 panels 与正文首次引用顺序一致。

**ver3 实际回应：** 说明已按 first citation 顺序重新 lettering。

**评分：4.5/5 — OK。**

回复针对性充分。与 Figure 2 的顺序修正采用同一规则；当前没有完整稿件，不能据此独立认证所有交叉引用。

<a id="r2-m8"></a>

### R2.m8｜candidate sites 与 candidates

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:390)

**Reviewer 原文**

> 8. In Supplementary Figure 5, what is the difference between "candidate sites" and "candidates"?

**本条需要回应什么：** 区分 allele 与 site，交代 multiallelic site、variant class 和 tensor 的对应关系。

**ver3 实际回应：** candidate = 一个 (node, offset, REF, ALT) allele；candidate site = 同一 node position、同一 SNV/INDEL class 的 allele 集合；每个 site 产生一个 tensor。

**评分：5/5 — OK。**

定义足够明确，与 R1.8 representative-allele labeling 及 R3.m1 units 一致，解决了原术语歧义。

<a id="r2-code"></a>

### R2.code｜依赖、安装失败与可运行发布

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:419)

**Reviewer 原文**

> The model has a series of dependencies that have prevented full testing, however, we are continuing to troubleshoot to get the model running on our platform.

**本条需要回应什么：** 让 reviewer 实际拿到当前 pipeline 并能在干净环境安装、运行最小例子和核对输出，而不仅是作者本地测试成功。

**ver3 实际回应：** 指出 default branch 仍为旧版及 compiled-extension fallback 两个问题；说明修复、CPU/GPU lockfiles、最小 chr1 测试、运行时间/内存及 tensor 一致性。但 release tag 与 Zenodo DOI 仍是占位符，fresh-clone 复测被写成提交前将要完成。

**评分：3/5 — 部分解决，需补充。**

这是当前最明确的未完成交付。作者已说明旧 default branch、compiled-extension fallback、CPU/GPU 环境和内部 smoke test 的结果，技术诊断有内容。但原审稿人遇到的是实际无法安装/运行；当前 release tag、Zenodo DOI 仍为占位符，fresh-clone 复测仍用将来时，因此外部可运行性尚未在回信中闭合。

这不是因“我在当前文件夹里没找到代码”而推断实验没有做。判断依据还包括作者自己的未完成表述。完整运行配置、逐模型表、实际 comparator 命令和已引用的补充结果没有进入本次可见材料，属于同一交付与可追溯性问题：在本条统一记录，不在 R3.m1、R3.m2 等处反复计为独立科学缺陷。后两条的 OK 评价针对方法解释充分性，并不表示其附件已经过核验。

保留 **3 分**：故障原因和修复方案已回应，但 reviewer 能拿到并运行的确切版本尚未确认。这个问题不要求新的生物学实验，却是正式提交前必须完成的实际工作。

**还需做：**

1. 给出 reviewer 实际可访问的固定 release/commit、模型、最小输入及输出示例，将占位符换为真实链接或有效审稿访问入口。固定 release、审稿专用归档均可；关键是可取用、可定位，不要求一定合入 main，也不把提交前公开 DOI 当成唯一途径。
2. 从该确切版本在干净环境执行 README 的最小例子，记录命令、依赖、期望输出与通过结果。既然回复声明支持并验证 CPU/GPU，就核验相应两种路径；无需重新运行全基因组 benchmark。明确示例从现成 GAM 开始，上游 reads alignment/graph preparation 另有 vg 和参考数据要求。
3. 将回复已引用/承诺的复现材料一并提供，并使主要表图能够追溯到模型和运行：逐模型配置及 split/seed/checkpoint/threshold，涉及不同 normalization、GPU 或 inference precision 的运行对应关系，实际 comparator 命令与模型/资源版本，以及所称 full training counts、strata/comparator tables、PR curves 等附件。可复用现有配置和日志；不要求另造一套复杂归档系统。未交付时不能把 “provided” 当作已核验事实。
4. 用最终独立提交包检查链接和附件可达性。ver3 的六张图目前仍引用 ver2_figures，而两个目录的对应图片完全相同；使提交包包含正确路径即可，不需重画。完成后再把回信中的 will release / will repeat 改为确已完成的事实。

## Reviewer 3：逐条评价

<a id="r3-general"></a>

### R3.general｜Reviewer 3 总评

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:427)

**Reviewer 原文**

> This manuscript presents Pansoma, a machine learning-based framework for tumor-only somatic variant calling that operates directly on pangenome graph alignments. The authors introduce a node-centric pileup representation, construct graph-derived tensors for machine-learning classification, and evaluate the framework using COLO829T, HG008T and SMaHT HapMap mixture datasets across short- and long-read sequencing platforms.
>
> Addressing the primary limitation of tumor-unmatched somatic variant detection, specifically, the discrimination of true somatic events from co-occurring germline alleles and alignment noise, this investigation evaluates the utility of a pangenome graph reference architecture. Incorporating population genetics directly into the alignment workflow rather than relying solely on downstream filtering rules provides a principled computational strategy. The authors substantiate this approach with concrete performance metrics: a sizable subset of benchmark germline variation is directly represented in the AF-HPRC graph, graph-guided alignment minimizes read-mapping errors and suppresses false-positive candidate calls, and the Pansoma workflow consistently outperforms existing methods in SNV calling precision and overall F1 performance across multiple validation benchmarks.
>
> Overall, I find Pansoma technically interesting and potentially important. To my knowledge, while graph-based alignment methods are well established, an end-to-end graph-native framework specifically developed for somatic variant calling represents a meaningful methodological advance.
>
> My main concern is not whether Pansoma performs well; the benchmarks indicate that it does. It is whether the manuscript demonstrates that this advantage is specifically attributable to the pangenome-native representation, the study's central conceptual advance. For this reason, I recommend major revision.

**本条需要回应什么：** 回答核心增益能否归因于 pangenome，而不是把 pipeline 整体优势直接当成 reference representation 的因果证据。

**ver3 实际回应：** 以 reference-level ablation、candidate composition、off-reference validation、complete-output metrics 和 scope narrowing 回应。

**评分：4/5 — OK。**

总回复准确指向三参考比较、候选组成、off-reference 验证和完整输出评价，后文也明确限定三参考分析为 system-level comparison。这里的 “address it directly” 可理解为直接处理该疑问，并不等于声称已完全隔离所有机制。结合 R3.1 和固定候选/比对的 channel 6 消融，已足以概括当前较窄的 pangenome-aware 主张。本段不能替代完整因果分解，但没有必要因沿用 “reference-level ablation” 这一名称而额外判为失败。

<a id="r3-1"></a>

### R3.1｜三参考比较是否隔离 pangenome 的贡献

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:439)

**Reviewer 原文**

> 1. The specific contribution of the pangenome warrants more direct isolation
>
> The study articulates a logical mechanistic framework wherein the pangenome structure accommodates baseline germline variants via alternative haplotypes. This graph-guided alignment minimizes non-biological discrepancies and mapping errors, thereby suppressing upstream candidate noise and translating directly into improved variant precision. Several components of this model are well supported by the current data-graph-based alignment increases the proportion of perfectly aligned reads and reduces substitutions, indels, and soft-clipping, while graph-informed candidate generation substantially lowers the number of putative variants.
>
> Nevertheless, the end-to-end benchmarking does not cleanly disentangle the contribution of the pangenome per se. Pansoma diverges from ClairS-TO and DeepSomatic-TO not only in its reference representation but also in alignment strategy, candidate generation, tensorization, model architecture, and post-processing. As the Methods note, the comparator callers employ their standard GRCh38 workflows, whereas Pansoma relies on vg graph alignment coupled to a bespoke graph-native pipeline. To definitively validate the manuscript's primary claim, a rigorous reference-level ablation study is strongly recommended. Specifically, the authors should compare the Pansoma pipeline's performance using a GRCh38 linear reference (or single-path graph equivalent) directly against the full HPRC pangenome graph, ensuring that candidate extraction, tensor generation, deep learning architecture, training setup, and filtering thresholds remain identical across both arms. Such an analysis is critical to elucidate graph-specific benefits from algorithmic optimizations, thereby demonstrating how much of the performance advance is directly driven by the pangenome architecture.

**本条需要回应什么：** 用同一 Pansoma 流程比较 GRCh38-only 与 HPRC，避免只靠跨工具 benchmark 将优势归因于图。交代实际控制了哪些变量、哪些随参考改变，并使机制结论与可隔离的证据一致。

**ver3 实际回应：** 同 reads、Giraffe、candidate/tensor/label/network/training/filter rules 比较三图，三 seeds，并报告 per-seed 和 bootstrap 结果。明确 20,000-node validation blocks 在不同图上落在不同 GRCh38 区域，彼此 validation 仅有 3–8%重叠；没有共享坐标划分重训，故定义为 system-level comparison。

**评分：4/5 — OK。**

本条现在可以接受，但其证据级别必须准确理解。作者已按 reviewer 的核心方向，在同一 caller 中实际加入 GRCh38-only 对照，而不是仅把 Pansoma 与不同架构的外部工具比较。reads、candidate/tensor/label/network 规则、training recipe、checkpoint/threshold selection rules 和 population filter 保持一致；两种 HPRC 图、三 seeds、六个肿瘤×平台设置给出相当一致的 SNV 方向：35/36 seed-matched 比较支持 HPRC，所有十二个 first-seed block-bootstrap 区间为正。这是有实质价值的 reference-level 系统证据。

其限制也是真实的：validation blocks 根据 node IDs 抽取，三图不共享相同 genomic split；initialisation 与 split 都随 seed 改变。约 90% 的训练区域共同、validation overlap 3–8%，能解释设计，却不等于这些因素完全控制。chr1 bootstrap 也不能消除训练/选择规则所产生的混杂。另一方面，独立抽取的 5% validation 并不自动构成偏向 HPRC 的选择设计；现有一致方向不能仅因未共享 split 就被判为没有意义。相同阈值选择规则产生不同数值阈值，本身不是不公平。

决定本条达到 4 分的，是上述实测证据与明确的 scope narrowing **共同** 起作用。作者已经写明“system-level comparison ... not an isolation ... with every other factor fixed”，未做共享 GRCh38 split，也不据此得出最终 INDEL 优势；再结合 R1.1 固定 alignments/candidates 后重训去掉 channel 6 的配对结果，可以支持：HPRC 在当前统一流程规则下提高这些设置的 SNV 表现，且 haplotype 信息具有额外预测价值。

仍不能从中精确计算 alignment、candidate suppression、context 和 haplotype channel 各自贡献了多少。R3.1 中“mainly alignment and candidate generation”有候选构成及去掉 channel 后仍有增益作为定性依据，不能读成已测得机制贡献比例。共享 split 会增强因果归因，但对已明确收窄的主张，不再将它列为本轮必须完成的实验。4 分表示核心疑问已有合理且有边界的回应，不是将不完全控制的实验评为完美因果证明。

<a id="r3-2"></a>

### R3.2｜候选减少是否选择性保留 somatic

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:465)

**Reviewer 原文**

> 2. Candidate reduction should be evaluated together with retention of bona fide somatic variants
>
> The reduction in candidate burden is one of the strongest mechanistic results. However, fewer candidates are beneficial only if germline variants and technical artifacts are preferentially removed while true somatic variants are retained.
>
> The authors currently emphasize overall candidate reduction, but the composition of removed candidates is not sufficiently resolved.
>
> I suggest comparing linear-reference and graph-based candidate generation for at least:
>
> 1) benchmark somatic truth variants,
>
> 2) matched-normal germline variants,
>
> 3) background/non-truth candidates.
>
> To mechanistic link graph representation with downstream accuracy, the authors should consider stratifying false-positive reductions based on whether the underlying germline allele is encoded within the HPRC pangenome graph. Such an analysis would directly connect the observation that over 79% of normal-sample SNVs and 73% of INDELs overlap with the AF-HPRC graph to the observed gains in precision, thereby demonstrating that graph-level allele representation is the primary driver of germline variant suppression.

**本条需要回应什么：** 分别比较 somatic truth、normal germline、background candidates；把 germline suppression 与 graph content 连接起来，同时量化真正 somatic 被移除的代价。

**ver3 实际回应：** 列三类 SNV/INDEL 的 retention/reduction，结合过滤后 FP；明确 ≥46/463 carrier threshold 是频率 proxy，而不是实际路径逐 allele 验证，披露 CHM13/GRCh38 unfiltered paths 引起的“no”组例外。报告 somatic INDEL 大量损失和长读 INDEL FP 增加。

**评分：4/5 — OK。**

核心三分组和 downstream FP linkage 已提供，且 graph-content proxy 的误分类来源及规模没有隐瞒。它支持“主要移除常见 germline SNVs”，不能等同于每一个 allele 都经过 topology membership 确证；ver3 现有定义已经把这一证据边界交代出来。SNV 与 INDEL 的不同结果也直接回答选择性问题。

<a id="r3-3"></a>

### R3.3｜graph-coordinate calls 的独立验证

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:506)

**Reviewer 原文**

> 3. Variants residing in graph coordinates require more rigorous validation
>
> Among Pansoma's most consequential capabilities is its capacity to retain variants that cannot be unambiguously encoded on GRCh38. The authors report that a subset of final predictions exists exclusively in graph coordinates and offer mutational-signature analyses together with representative read-level examples as evidence of biological plausibility. This is a compelling and potentially important result; however, the supporting evidence remains largely indirect.
>
> The authors are encouraged to perform orthogonal validation on a representative sample of graph-exclusive variants by integrating complementary data sources, including matched-normal alignment evidence, phased haplotype assemblies, or cross-platform validation datasets. Demonstrating the authenticity of these unique calls through independent confirmation would significantly strengthen the manuscript's thesis, proving that pangenome reference architectures reliably extend variant detection into genomic regions inaccessible to conventional linear alignment workflows.

**本条需要回应什么：** 用 normal、assembly、cross-platform 等直接证据验证代表性 off-reference calls，不能仅靠 signatures/个例声称扩大了 somatic 检测范围。

**ver3 实际回应：** 对 12 runs 的 chr1 calls 用 assemblies、normal reads 和其他平台验证，较大集合随机抽 3,000，并设 placed TP/FP controls。区分 population-first classification 与 assembly-only tumor support，承认 COLO829T 自有 reads 组装非独立。多数 calls 是 population/normal supported，只有小子集可能 somatic。

**评分：4.5/5 — OK。**

已达到 reviewer 建议的验证层次，并让结论随不利结果收窄。normal depth≥5 时没有 ALT 只是 suggestive，文中没有当作绝对证伪；cross-platform 三者共有也没有直接等同 somatic。默认关闭 rescue 与这些结果一致。Reviewer 原先期待的广泛 off-reference 优势并未被证明，但回复已充分说明实际可支持的范围。

<a id="r3-4"></a>

### R3.4｜完整输出指标、graph-only 比例及 CHM13

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:522)

**Reviewer 原文**

> 4. The reported benchmark metrics do not evaluate the complete Pansoma output
>
> Conventional TP/FP/FN evaluation is performed only on Pansoma calls that can be projected back to GRCh38. Calls that remain exclusively in graph coordinates are analyzed separately. Thus, the reported precision, recall and F1 values represent the GRCh38-projectable subset of Pansoma calls rather than the complete output of the method.
>
> This distinction is important because the graph-coordinate fraction is platform dependent and can be significant, particularly for ONT data. Supplementary Table 6 reports graph-coordinate fractions exceeding 10% and, in some settings, 20%. In addition, the statement in the Discussion that only “about 2%” of variants remain on alternative graph paths appears inconsistent with the higher fractions reported in Supplementary Table 6 and should be re-evaluated.
>
> An additional consideration is whether calls classified as graph-coordinate are truly dependent on graph representation or are simply not projectable to GRCh38. The authors note that the projection framework can be extended to CHM13. It would therefore be informative to determine what fraction of GRCh38-unprojectable calls can be represented on the T2T-CHM13 reference, and what fraction remains genuinely graph-specific even when a more complete linear reference is considered. Such an analysis would help distinguish limitations of GRCh38 from advantages that are intrinsic to the pangenome graph representation.

**本条需要回应什么：** 让 GRCh38-evaluable 指标与完整输出负担分开；纠正约 2%的概括，并检查不可投影至 GRCh38 是否只是参考不完整、可转到 CHM13。

**ver3 实际回应：** 明确固定 GRCh38 truth 下将全部 graph-coordinate PASS calls 计 FP 的 conservative evaluation，不能测未知 truth 的 recall。列 rescue on/off 的 graph-coordinate share 与 F1、CHM13 representation，以及 flanking bubbles 的 GRCh38 recovery。

**评分：4.5/5 — OK。**

三部分均有定量回答。默认 graph-coordinate share 0.5–10.1%、rescue on 8.0–68.0%，不利的完整输出比较也保留；具有 CHM13 allele 不再被等同 graph-intrinsic novelty。把所有 unknown calls 计 FP 是惩罚性边界而不是真实全空间 accuracy，ver3 第一段已明确这个限制，本条可接受。

<a id="r3-5"></a>

### R3.5｜truth 与 graph-coordinate SBS96 谱

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:546)

**Reviewer 原文**

> 5. The Figure 4 mutational-signature analysis should incorporate the truth set directly
>
> In Figure 4a, the overlap analysis includes the benchmark truth set, yet the mutational-signature stacked bars presented immediately below encompass only Pansoma, ClairS-TO, and DeepSomatic-TO, with the ground-truth signature entering only indirectly as the reference for cosine-similarity calculation. I recommend displaying the spectra side by side, ground truth alongside Pansoma, DeepSomatic-TO, and ClairS-TO, so that the specific signatures underlying inter-caller agreement or disagreement can be readily identified. This is particularly relevant for HG008T, whose mutational profile is heterogeneous and for which the reported cosine similarities are comparatively modest.
>
> To strengthen the claim that graph-coordinate variants reflect true biological events through mutational-signature concordance, the authors should provide the underlying 96-channel single base substitution (SBS96) profiles in the Supplementary Information.

**本条需要回应什么：** 在图中直接并列 truth 与各 caller 的谱，并公开 graph-coordinate calls 的 SBS96；避免只报告 cosine 数字或通过 signature 相似性推断 somatic。

**ver3 实际回应：** 提供 HG008T、COLO829T 两张 truth/normal/三 caller SBS96 图，另有六个 graph-coordinate spectra；明确后者五个更像 germline，撤回以谱吻合支持 somatic 的论证。

**评分：4.5/5 — OK。**

已检查三张实际图，panel 的 tumor/platform、truth counts、calls counts 和正文表一致。graph-coordinate 六格 counts 相加为 6,018，恰等于抽样后 6,025 减 7 个无可用 context 的 calls；COLO829T ONT 的 3,000 抽样已说明。用谱揭示原解释不成立是有说服力的回应，而不是失败。

<a id="r3-6"></a>

### R3.6｜perfect alignment 与 correct placement

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:559)

**Reviewer 原文**

> 6. More perfect alignment does not necessarily imply more correct alignment
>
> The manuscript interprets increased perfect alignment and reduced mismatches/soft clipping as evidence of improved alignment to the pangenome. This interpretation is plausible.
>
> However, because the pangenome contains more alternative sequence paths than GRCh38, an increase in perfect alignments does not by itself prove that reads are mapped more correctly.
>
> To isolate the true mechanism behind the observed performance improvements, a control comparison evaluating alignment quality across GRCh38, T2T-CHM13, and the full HPRC graph is strongly recommended. Such an analysis is crucial to delineate improvements stemming simply from reference completeness (i.e., linear T2T assembly) versus those uniquely conferred by graph-based multi-haplotype modeling.
>
> If feasible, the authors could provide additional evidence using:
>
> MAPQ distributions,
>
> phased haplotype concordance,
>
> assembly-based truth,
>
> repetitive or segmentally duplicated regions,
>
> other difficult genomic regions.
>
> This would help distinguish “more perfectly aligned” from “more accurately aligned.”

**本条需要回应什么：** 用已知来源或其他独立证据测量定位正确性，加入 CHM13 控制以区分 reference completeness 与 multi-haplotype 的作用。

**ver3 实际回应：** 从不在图中的 HG008-N diploid assembly 模拟 chr1/20 reads，在 GRCh38-only、CHM13-only、HPRC v1.1/v2.1 上评估 misplaced、MAPQ 和 difficult regions。报告 v2.1 ONT 更差、SV-spanning reads 多数不可评分及双参考 one-to-one 区域的限制。

**评分：4.5/5 — OK。**

符合 reviewer 指定的对照逻辑。它支持的是共同可评估区域内 Illumina/HiFi placement 改善，不能覆盖 CHM13 新增序列空间；作者明确承认这个外推边界。calling 层没有 CHM13 arm 也未被掩盖，因此不能以“没有覆盖所有参考完整性”否定本条已做的直接测试。

<a id="r3-7"></a>

### R3.7｜不能凭图比对改判 false negatives

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:599)

**Reviewer 原文**

> 7. The conclusions drawn from specific benchmark false negatives should be stated with greater caution, as the current interpretation may be overly definitive.
>
> Supplementary Figure 3 characterizes several benchmark false-negative variants and concludes that they reflect GRCh38 reference-bias artifacts rather than genuinely missed variants, on the grounds that the corresponding reads align perfectly to alternative graph paths. While this interpretation is plausible, graph-alignment behavior alone may be insufficient to overturn established benchmark truth annotations. Reclassifying these events requires orthogonal validation, such as matched-normal alignment evidence or haplotype-resolved assembly data. In the absence of such supporting evidence, the authors are advised to adopt more cautious phrasing and qualify their claims accordingly.

**本条需要回应什么：** 撤回或严格限定 reference-bias artifact 的判断，并用 matched-normal、haplotype assembly 等证据审视原例子。

**ver3 实际回应：** 完全撤回原先概括，给出 v1.1 原定义下的 normal/tumor assemblies、normal reads 和 haplotype walks；多数吸收 allele 不在 normal 中，存在真正 somatic。说明 COLO829T tumor assembly 是自有 long reads 组装及 homopolymer assembly 的局限。

**评分：4.5/5 — OK。**

该回复接受证据所指向的相反结果，并保留 truth 分母，避免通过重新分类抬高性能。v1.1 原定义与当前 v2.1 严格定义不同已在补充说明，不能用其不同分母指控前后矛盾。“有 tumor assembly，但缺独立于被分析 reads 的 assembly”也已区分。

<a id="r3-8"></a>

### R3.8｜肿瘤应用与正常组织 mosaicism 的范围

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:612)

**Reviewer 原文**

> 8. The validation is predominantly tumor-centric; applicability to somatic mosaicism in normal tissue remains unresolved
>
> The manuscript introduces Pansoma broadly as a somatic variant caller, yet its biological validation rests on tumor-derived samples and a synthetic low-VAF admixture. COLO829T is a melanoma-derived benchmark and HG008T a pancreatic cancer sample, while the low-VAF HapMap benchmark comprises germline variants from minor cell-line constituents employed as controlled proxies for somatic-like events. No authentic normal-tissue somatic-mosaicism dataset such as SMaHT tissues appears to be evaluated. This distinction is material, as somatic variants arising in normal tissues occupy a substantially different biological context and frequently present at low variant allele fractions.
>
> If Pansoma is intended specifically for tumor-only somatic variant calling, the manuscript should articulate this scope consistently and refrain from broader claims. If, conversely, the authors intend general applicability to somatic variation beyond malignancy, validation on an independently characterized normal-tissue somatic-mosaicism dataset would considerably strengthen the work.

**本条需要回应什么：** 明确只做 tumor-only cancer calling，或者提供真正 normal-tissue mosaicism 验证；不能从 mixture/低纯度肿瘤直接外推。

**ver3 实际回应：** 选择收窄范围：title/Abstract/Results/Discussion 只限 cancer，说明未训练或验证 normal tissue，0.08 floor 不适合支持低 VAF mosaicism；将 dedicated data/models 列为未来工作。

**评分：5/5 — OK。**

这是 reviewer 明确允许的回应路径，边界清晰。无需为了拿高分强制增加 SMaHT normal-tissue 实验。

<a id="r3-9"></a>

### R3.9｜训练稳定性、bootstrap 与 operating point

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:622)

**Reviewer 原文**

> 9. Assessment of Model Stability and Statistical Uncertainty
>
> The selected performance metrics are appropriate; however, the presentation relies almost exclusively on point estimates. Given the inherent susceptibility of machine learning models to variation during selection and training, a formal assessment of model stability is necessary. The authors are encouraged to evaluate uncertainty using approaches such as multi-seed retrainings, block-bootstrap confidence intervals across genomic intervals, or chromosome-wide variance estimates to ensure the robustness of the reported performance metrics.
>
> A particularly informative analysis would be to retrain Pansoma across multiple random seeds and report the distribution of precision, recall and F1 on both the held-out chromosome and the independent HG008T benchmark. In addition, chromosome-level or genomic-block re-sampling could be used to provide confidence intervals for the performance difference between Pansoma and the comparator callers, rather than reporting only point estimates. Sensitivity of the results to the fixed 0.5 classification threshold would also help establish that the reported advantage is not dependent on a single operating point.

**本条需要回应什么：** 报告 seed variation 与 genomic sampling uncertainty，给出与 comparator 的差值区间，说明阈值敏感性和 held-out tumors 的证据强度。

**ver3 实际回应：** 三参考 SNV 与最终 INDEL 均有三 seeds；SNV 以首 seed 做 paired 1-Mb block bootstrap，其余 seeds 补充方向一致性。明确 seed 同时改变初始化与 validation split，跨肿瘤 SNV 为单次 exploratory，chr2 INDEL CI 仅覆盖 genomic sampling。R3.9b 解释第二集群 full-precision 预测与 R3.9a GPU call set 的微小差异。

**评分：4/5 — OK。**

本次不能再把 HG008T PacBio 的 0.839−0.826=0.013 与表内 +0.009 当成未解释错配：R3.9b 已明确该 run recall 为 0.835、与原集最多差 3 TP/5 FP。主表三 seed mean 与 ΔF1 首 seed 的区分也标明了。54 组 seed-level P/R/F1 与展示精度一致；chr2 六格 F1 差的均值重算为 +0.0345 和 −0.0205，与 +0.034/−0.021 在输入取整范围内相容。未联合传播训练和 genomic uncertainty、只评 chr1 等限制均有声明，所以当前可接受但不是强泛化保证；其 run identity 的归档在 R2.code 作为共同交付事项处理。

<a id="r3-10"></a>

### R3.10｜低 VAF 和 INDEL 的结论边界

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:637)

**Reviewer 原文**

> 10. Need for Qualified Interpretation of Low-VAF and Insertion/Deletion Benchmarks
>
> The utility of the HapMap benchmark for low-VAF assessment is acknowledged; nevertheless, because these synthetic constructs rely on mixed germline lines rather than genuine somatic lineages, the generalizability of the low-VAF results should be stated with appropriate caution. Additionally, the manuscript should clearly distinguish SNV performance from INDEL accuracy. Pansoma's competitive edge is substantially more pronounced in SNV calling, whereas INDEL accuracy remains modest and is surpassed by ClairS-TO in certain ONT settings. Overall performance claims should be refined to accurately reflect these boundaries.
>
> To better define these limits, the authors could complement the HapMap mixture with either a bona fide low-VAF somatic benchmark or an in silico dilution series based on validated somatic variants from COLO829T or HG008T. In addition, INDEL performance should be stratified by insertion/deletion class, event length, sequence context such as homopolymers or tandem repeats, and sequencing platform. These analyses would help determine whether the observed limitations reflect low-VAF biology, platform-specific error modes, or particular classes of indels.

**本条需要回应什么：** 分开 SNV/INDEL 主张，补真实 somatic dilution 与按 insertion/deletion、长度、repeat、平台的分层，并承认低 VAF mixture 的外推限制。

**ver3 实际回应：** 不再宣称 INDEL 优势；给出 final models 的低 F1、目标肿瘤加入/未加入训练的区别、chr2 held-out confirmation 和逐 stratum TP/FP/FN。低纯度改由 COLO829T dilution 支持，HapMap 限为补充。

**评分：4/5 — OK。**

在“证据是否支持当前结论”上已经充分回应：作者撤回普遍 INDEL 优势，披露 held-out INDEL 对 ClairS-TO 六个肿瘤×平台设置全部较差，并列出 candidate loss、population filtering、homopolymer errors 等瓶颈。这些局限会限制方法适用范围，但不等于该回复没有回答质疑。dilution 证据支持低 VAF 能力不足的结论，不能被解释成精确的通用 detection limit，见 R2.4。

表 R3.10 的三 seed 平均计数被取整，不能要求取整后的分项严格加和为 all，也不能要求平均 P/R 的调和平均等于平均 F1。完整 comparator strata 的附件未在本次包中出现，交付核对集中在 R2.code，不额外重复扣分。

<a id="r3-m1"></a>

### R3.m1｜参数统一与逐模型可复现材料

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:652)

**Reviewer 原文**

> 1) Additional methodological details and explicit reconciliation are required
>
> The methodological description appears disproportionately focused on the neural-network architecture, whereas several upstream and downstream steps that directly influence calling performance are comparatively condensed.
>
> More importantly, several key parameters differ between the main Methods and Supplementary Methods.
>
> The main Methods describe:
>
> non-reference AF >10%,
>
> 2 supporting reads,
>
> exclusion of reads with MAPQ ≤5.
>
> The Supplementary Methods instead describe:
>
> AF ≥0.08,
>
> ≥3 alternative observations,
>
> mean alternative BQ ≥10,
>
> MAPQ ≥10,
>
> INDEL length ≤50 bp.
>
> Similarly, the main Methods describe 80 training epochs, whereas the Supplementary Methods state 80–120 epochs.
>
> Resolving these discrepancies as well as more detailed description are essential, as they directly impact the study’s reproducibility and the overall interpretation of model performance.

**本条需要回应什么：** 统一主 Methods 与 Supplement 的真实参数，明确上游候选、张量、标签、网络、选择规则及新增特征，并说明参数与所报告模型的对应关系。可复现附件的实际可取用状态与方法文字充分性分别判断。

**ver3 实际回应：** 明确 node selection MAPQ>5、candidate MAPQ>10、≥3 ALT、每条 BQ≥10、AF≥0.08、1–50 bp，SNV 20 epochs/INDEL 12 epochs；详细定义最终 INDEL 五个 site features、缩放、缺失值与融合层。承诺 per-model Table R3.m1、logs、validation lists 与 checkpoint 归档。

**评分：4/5 — OK。**

原审稿人列出的关键不一致——AF、ALT 支持数、BQ、MAPQ、INDEL 长度和 epochs——已有逐项统一值。node selection 的 MAPQ>5 与 candidate reads 的 MAPQ>10 属于不同步骤，strict inequalities 也已明确；SNV 20 epochs、INDEL 12 epochs 取代旧稿口径，不能因数值与旧稿不同就当作仍未统一。最终 INDEL 五个 features 的输入来源、计算/缩放、缺失值、5→64 融合层，以及 population prior 的含义，都已清楚说明。

所以，本条最核心的“参数究竟是什么、方法能否理解”已经充分回应，评分为 4。完整 revised manuscript、逐模型 Supplementary Table R3.m1、logs/configs 不在当前可见包中，本次不能认证它们已与回复完全一致；这属于实际交付的核验边界，集中列入 R2.code，而不把它重复计作本条未回答。这个处理也与图号/legend 等其他已具体说明、但完整改稿不可见的条目保持同一标准。

同样，旧口径的 61/702 与 current benchmark 60/697 已由 Supplementary Note R2.3 的 benchmark-region 范围差异解释，不再留作未说明的数字冲突。这里的 OK 并不表示仅凭一张全局参数表就已独立复现所有训练。

<a id="r3-m2"></a>

### R3.m2｜comparator 配置、过滤和模型选择公平性

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:702)

**Reviewer 原文**

> 2) The setup and parameters for comparator workflows should be documented in greater detail to ensure fair benchmarking.
>
> Pansoma is described in considerable technical detail, including the vg giraffe version, but the comparator workflows are largely described as their “standard GRCh38-based alignment and calling workflows.”
>
> For reproducibility and fair benchmarking, the authors should provide:
>
> software versions,
>
> aligners and versions,
>
> pretrained checkpoints/models,
>
> exact tumor-only configurations,
>
> PASS/filtering criteria,
>
> reference resources,
>
> population/PoN resources,
>
> command-line parameters.
>
> This is particularly important because the main advantage reported for Pansoma is increased precision, which can be sensitive to filtering configuration.

**本条需要回应什么：** 给出 caller/aligner/checkpoint/数据/资源版本、确切 tumor-only 命令、PASS/filter 规则及训练暴露；确认所选模型与输入技术匹配，性能优势不由不透明的配置选择产生。

**ver3 实际回应：** 给 DeepSomatic 1.9.0、ClairS-TO 0.4.1/0.4.2 与平台 model 名、provider aligners、资源和 allele/position filter 差异。说明 native PASS 叠加 common filter 不是总过滤完全相同，PR curves 另取非 RefCall/NoCall records；纠正 HG008T ONT DeepSomatic 模型，交代 within/held-out exposure。

**评分：4/5 — OK。**

原始意见要求的核心配置已相当具体：软件与 aligner 版本、平台 checkpoint 名、tumor-only 模式、population/PoN 资源、native PASS 与附加 common filter、normalization/evaluation 以及 training exposure 均有说明。HG008T ONT 的 DeepSomatic 模型选择已纠正。尤其重要的是，ver3 已消除之前的 PR 规则矛盾：曲线纳入除 RefCall/NoCall 外的 records，忽略 native FILTER，应用共同 population filter，再按 QUAL 排序；没有一边说解除 native tags、一边又保留 native technical filters。

COLO829T Illumina 使用 ilmn，而另两组使用 ilmn_ssrs，是已公开的配置差异。它可以影响结果，但差异本身不足以证明作者挑选了弱 comparator。公开模型建议也不能替代对特定训练暴露、数据匹配与实际运行的判断：[ClairS-TO 官方模型说明](https://github.com/HKU-BAL/ClairS-TO#pre-trained-models)。现有材料没有支持“已确认不公平”这一结论，因而不再把新增的 ilmn/ilmn_ssrs 定向实验设为本条通过的必要条件。

matched-precision 结果确实不对称：Pansoma 曲线与 comparator 固定工作点相比。但它回答的是“能否在该 comparator 所达到的 precision 下取得更多 recall”，这是有效的工作点分析；它不能说明两工具完整 PR 曲线孰优，也不是各自最优阈值下的能力上限。作者已明确这一计算方式，主证据也不是仅靠该指标。实际每次运行的 exact commands 仍需作为已承诺材料交付，统一在 R2.code 记录。本条 4 分认可方法透明度与原问题的解决，不认证未见的 flags 或模型文件。

<a id="r3-m3"></a>

### R3.m3｜图版本、组成、过滤、SV 与 donor 关系

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:742)

**Reviewer 原文**

> 3) The composition of the pangenome reference should be specified in greater detail
>
> The manuscript refers throughout to an allele-frequency-filtered HPRC graph, and Figure 5 specifies the HPRC v1.1 AF-Graph. Because the pangenome constitutes the central methodological element, the Methods should document the HPRC release and version, the haplotypes included, the allele-frequency filtering criteria, the treatment of structural variants, and the graph-construction procedure. It should further be stated whether any benchmark individuals or haplotypes closely related to them, contributed to the graph used for evaluation. Because the manuscript notes that diagnostic efficacy hinges on graph composition and haplotype representation, disclosing this information is essential to assessing the external validity and translational utility of the proposed pipeline.

**本条需要回应什么：** 说明 release、construction、haplotypes、AF filter、SV treatment，以及测试样本/近亲是否进入图。

**ver3 实际回应：** 列 v1.1 44 individuals+两条参考=90 haplotypes，v2.1 231 individuals+两条参考=464；d9/d46 与参考路径不滤、SV 同规则、Pansoma 只调用1–50 bp。六条 cancer lines 无直接 sample-list overlap；三个 normals 另做 rare-allele sharing，并承认不能正式排除亲缘。

**评分：4.5/5 — OK。**

已覆盖 reviewer 列出的要素，且对亲缘证据强度没有过度承诺。R3.2 的 463 是用于 carrier-frequency proxy 的计数口径，不能直接与包含 GRCh38 的464条 graph paths 当成同一个计数对象而宣布矛盾；proxy 与实际图内容的区分已在 R3.2 明示。

<a id="r3-m4"></a>

### R3.m4｜术语和语言一致性

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:758)

**Reviewer 原文**

> 4) Terminology and Language Consistency
>
> Throughout the manuscript and Supplementary Information, terminology is not always used consistently, and several minor grammatical and stylistic issues are present. The authors are encouraged to perform a careful language and consistency check, particularly for method names, abbreviations, figure/supplementary nomenclature, and technical terminology.
>
> The term “panel of normals” is used broadly for resources that include gnomAD, dbSNP, a 1000 Genomes PoN and CoLoRSdb. These are not all conventional PoNs; “population/background filtering resources” may be more precise.
>
> “One key challenge in performing variant calling on graph alignments is it differs fundamentally…” is grammatically strange and incorrect.
>
> “Pansoma convert the variants…”
>
> “with the technique details described in Methods”
>
> “read pile up” in Figure 5
>
> The expansion and terminology for NPU should be standardized between the main text and Supplementary Methods.
>
> Pansoma versus Pansoma-TO, and “Supplement Figure” versus “Supplementary Figure”, should be standardized throughout.
>
> The manuscript would benefit from a careful final language edit

**本条需要回应什么：** 区分 population databases 与传统 PoN，统一 Pansoma、Supplementary nomenclature、NPU 和具体语法问题。

**ver3 实际回应：** 统称 population filtering resources，仅对真正 PoN 保留该名；新版不产生 NPU，删除过时格式说明；统一方法名/图表名并逐一说明语法修改。

**评分：4.5/5 — OK。**

回应覆盖具体意见，且 NPU 的删除有方法变更原因。当前 response 的主要术语已与此一致；完整稿语言修改不可在未见稿件时认证，但这不妨碍判定本条回复充分。ver3 图资源目录的交付问题在 R2.code 处理。

<a id="r3-overall"></a>

### R3.overall｜Reviewer 3 的最终总体意见

[查看 ver3 对应原文](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md:782)

**Reviewer 原文**

> Overall assessment
>
> In summary, I consider Pansoma a novel and potentially valuable framework for graph-native somatic variant calling. The improvement in tumor-only SNV precision is encouraging, and several analyses support the proposed benefit of pangenome-aware read representation.
>
> The primary issue is whether the manuscript currently demonstrates sufficiently and directly that this advantage is specifically attributable to the pangenome itself, rather than to the broader Pansoma architecture and workflow.
>
> A direct reference-level ablation, a clearer accounting of which candidates are removed by the graph, stronger validation of graph-coordinate calls, and reconciliation of the methodological inconsistencies would substantially strengthen the study.
>
> In addition, the current validation is primarily cancer-centered. Either the scope should remain tumor-only, or broader applicability to somatic mosaicism should be supported by additional data.
>
> I therefore recommend major revision.

**本条需要回应什么：** 综合证明图的贡献、候选的去留、off-reference validation、方法可复现性，并将 biological scope 控制在癌症。

**ver3 实际回应：** 逐项引用 R3.1–R3.4、R3.m1、R3.8 概括已做工作。

**评分：4/5 — OK。**

候选去留、off-reference 验证、SBS96、完整输出负担、参数和 cancer-only 范围均已有针对性证据；三参考系统比较加上 haplotype-channel 消融足以支撑目前收窄的 SNV 方法学主张。原来更强的 graph-native 机制分解、广泛 off-reference somatic discovery 和泛化/低 VAF/INDEL 优势没有获得全面支持，但作者也没有继续将它们作为已证实结论。

因此，总体回复质量达到 4 分。R1.8 的选择流程表述仍需针对性修正；软件和补充材料的可核验交付见 R2.code。本条不再另开一组实验或因同一交付问题重复计分。这个判断表示支持在有限修订与材料落实后继续审理，不等于仅凭 response letter 可以建议期刊直接接收。

## 剩余必要修改及其完成标准

这里只汇总两条未达到 OK 的事项，不给已经 OK 的条目追加任务。

| 对应条目 | 必须处理的内容 | 怎样才算完成 | 不自动要求的工作 |
|---|---|---|---|
| R1.8，方法解释与结论 | 区分 fixed-checkpoint threshold 与包含 checkpoint selection 的完整流程；解释 unknown-negative 假设及 final INDEL rule 5 的 non-somatic/other 身份。 | 回复的正文与补充规则一致，承认未测 selection sensitivity，当前性能结论限定为实际采用的方法；不再声称整个选择流程对这些标签不受影响。 | 不要求全套重训。若保留“模型选择/指标不受影响”的更强结论，才需要针对性分析支持。 |
| R2.code，共同交付 | 可访问版本、模型/最小例子、fresh-clone 测试；兑现已引用的配置、命令、表图和可追溯材料。 | reviewer 可按实际命令取得并运行对应版本，关键结果能找到相应配置与附件，链接不再是占位符。 | 不要求合入 main、全基因组重跑或重画相同图片；公开 DOI 不是唯一审稿访问方式。 |

这两条的重要性不同：R1.8 需要使方法描述与推断逻辑准确，但现有证据没有证明它会逆转主要 SNV 结论；R2.code 直接决定 reviewer 是否能检查和使用软件，是明确的提交条件。它们都不能靠删去一段不利说明来“解决”。

## 已核查、不再列为问题的事项

| 事项 | 当前证据与判断 |
|---|---|
| channel 6 消融是否混配 GPU controls | 每对 with/without 的 GPU 数匹配；部分 v1.1 例子是一 GPU 配对，不能拿它们与四 GPU 主表混算。 |
| path-spanning context 是否已有独立收益证明 | 实际 tensor 支持跨节点表征；当前已明确未隔离其性能效应，因此不当作虚假的消融结论。 |
| HapMap conditional 分母 | graph-carried allele、span ±1 bp mask 的连带排除与一个 N-base 例外可重建分母。 |
| exact germline / partial somatic 优先级 | 有序规则已明确，不能继续说顺序未交代。 |
| germline confidence 外的候选是否一定不能当 somatic-negative | 不能这样推断；需区分 somatic truth 的适用性与 germline/artifact 的类别确定性。R1.8 剩余修改已按此限制。 |
| ≥46/463 是否等同实际 path membership | 已定义为 carrier-frequency proxy，并披露 unfiltered reference-path 例外。 |
| R3.9a/b 的 +0.013 与 +0.009 | full-precision 第二集群的 call set、0.835 recall 及小幅 TP/FP 差异已有说明，非未解释算术错误。 |
| 三 seed 平均、取整与 F1 | 非线性平均和取整会导致 mean P/R 与 mean F1、分项整数与总体之间的表面不等，不能直接据此判错。 |
| graph-coordinate SBS96 counts | 6,018 = 6,025 − 7；六图计数、抽样和 context 排除一致。 |
| COLO829T tumor assembly 是否存在 | 有自有 long reads 的组装；缺的是独立于被分析 reads 的验证来源，原文已区分。 |
| ClairS-TO 默认 INDEL 区域限制 | 官方版本记录说明 v0.3.0 后取消该默认限制，不能用旧规则指控当前 0.4.1/0.4.2 配置。[官方版本记录](https://github.com/HKU-BAL/ClairS-TO#latest-updates) |
| 低 VAF 或 ilmn/ilmn_ssrs 的新增疑问 | 均不作为已证实错误；对证据边界的说明保留在 R2.4、R3.m2，不追加必做对照。 |
| 正常组织 mosaicism | reviewer 允许收窄范围；cancer-only 的回应路径已完成，不要求另起正常组织项目。 |

## 专业审稿总体意见

当前回复支持的核心研究结论是：**Pansoma 将 pangenome alignments、候选生成与 haplotype-aware 输入结合，在已测试的 tumor-only SNV 设置中呈现有意义的 precision/F1 收益；这些收益伴随真实的 allele absorption、训练分布依赖和低 VAF/INDEL 局限。** 证据来自实际对照、候选与错误构成、独立 truth 的测试、已知来源的 alignment 模拟以及对 off-reference calls 的重新验证，不只是作者自述或更换措辞。

从严格但合理的审稿角度，这已经超过“主要靠解释避开问题”的程度。尤其正常支持/assembly 检查纠正了原先将 absorbed alleles 当作 benchmark error 的解释；HapMap 与 off-reference 结果也不再承担它们不能支撑的生物学主张。这些改变使较窄结论更可信，但不会反过来证明更大的技术或临床适用范围。

我仍保留三项对科学贡献的实质判断：三参考实验是系统证据而非完全因果分解；跨肿瘤泛化的独立性与数量有限，平台重复不能替代生物样本重复；方法的净价值主要落在当前支持的 SNV 场景，不能由 INDEL、低 VAF 或 reference-external discovery 补足尚未验证的广度。这些都是阅读完整论文时评估重要性的边界，并非在本次 OK 条目后暗中新增的待办。

**因此维持本次总体 4/5：回复质量已达到总体可接受、需要有限修正的水平。** 不给 4.5–5，是因为 R1.8 还有准确性问题、实际交付未闭合，而且当前材料不足以认证全文与所有附件；不再给 3.5，是因为没有证据把未完成的严格机制分解或尚未证实的配置疑问视为对当前较窄结论的主要反证。若以我作为技术审稿人的判断，我会支持完成这两类事项后继续审理，而不会仅因尚有可增强的实验就要求无限追加。期刊的最终接收决定仍须结合完整改稿的重要性、可用性和编辑判断。

## 文件与证据记录

- 主评估输入：[Reviewer_Letter_Response_ver3.md](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md)，SHA-256：`3c75162663cdb4d5cd830a4657e379485c67c873fa8ab615bdafbf79cdb885f2`。
- 图目录：[Reviewer_Letter_Response_ver3_figures](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver3_figures)。六张图此前已逐张检查；本次确认原始文件未变。Markdown 引用的 ver2 图与对应 ver3 图逐字节相同。
- 本次仅更新此 evaluation，并保留更新前的 evaluation 备份；未修改 ver3/ver2 response、DOCX、图或旧 assessment。
- 当前结论没有将未见的逐运行命令、完整补充表或主文当成已核验事实；没有重新运行训练和 benchmark。
