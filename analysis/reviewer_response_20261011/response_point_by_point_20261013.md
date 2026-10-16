# Pansoma: NBT 审稿意见逐条回应（2026-10-13 版）

这份文件把每条意见的 **状态**（已解决 / 部分解决 / 待解决）、**更新后的 Response** 和 **证据** 放在一起。

- **数据来源**：2026-10-12 的 RIS 结果（`analysis/benchmark_7fc7396_20261012`），以及今天在本地做的实验（`experiments/`）。
- **和其他文件的关系**：
  - `response_to_reviewers_detailed.md`（10-12 详版）仍是实验设计的出处；
  - 每条的 Response 以本文件为准；
  - 简版 `response_to_reviewers.md` 还没同步。
- **单条原文**：`points/<R-ID>.md`。本文件第 6 节是把它们按顺序拼起来的。

状态的定义：
- **已解决**：数字都来自已完成的分析，只剩改稿。
- **部分解决**：主要论点已有数据支撑，部分数字或实验还缺。
- **待解决**：关键证据还没有。

数字默认指 chr1（test chromosome），PoN 后，benchmark BED 内；例外的地方都标明了。

## 1. 状态总表

已解决 5 条，部分解决 24 条，待解决 3 条。

| Point | 问题 | 状态 | 已有证据（关键数字） | 还缺什么 |
|---|---|---|---|---|
| R1.1 | "pangenome-native" 证据不足 | 部分解决 | **E1.1**：18 个模型上做了 ch6 counterfactual；把 ch6 设成 H 后，HG008T 的 somatic AP 从 0.20–0.35 降到 0.06–0.09，germline 误判增加 3.4–7.6 倍。**E3**：HPRC 在 SNV 上比 GRCh38-only 高 0.04–0.32 | 去掉 ch6 重训（GPU RIS）；bootstrap CI |
| R1.2 | 泛基因组组成、单倍型数、HPRC v2.1 | 部分解决 | d46 ≈ d9（PoN 后 −0.021 到 +0.031）。**E7**：决定吸收多少的是频率阈值，不是单倍型数 | matched d9 重训（GPU）；v2.1 VCF 下载（决定） |
| R1.3 | 困难区域；ONT 的歧义 | 部分解决 | **E9a/c/d**（3 个 d46 run）：ONT 的 unplaced calls 主要是人群 allele 和 germline，tumor-only ≤ 1.7%（只看组装） | 按区域分层的 calling 指标（E10）；E11b 比对（被 hold） |
| R1.4 | 图吸收已有的 somatic | 部分解决 | **E5b**：被吸收的 truth 有 95.5–100% 本来就会被 PoN 标掉；吸收只占图臂额外 INDEL 损失的小部分 | linear 臂能否叫回这些（数据同步）；E6 / E25（决定） |
| R1.5 | HapMap mixture、数据泄漏 | 部分解决 | 用 COLO829T dilution 替代 HapMap（fiberseq 干净）；HapMap donor 在图里的情况已核实（v1.1 中 4/6，v2.1 中 6/6） | HapMap 去留（决定）；Illumina / ONT 的 BQ 伪影修复（你在做）；HCC1395 dilution |
| R1.6 | sensitivity 低 | 部分解决 | recall 漏斗；PoN 上限（HG008T INDEL ≤ 0.35） | 比对器按同一协议重打分（等路径）；INDEL held-out 阈值 |
| R1.7 | 训练规模、肿瘤多样性 | 部分解决 | **E21** 训练集组成。**E12a** 跨样本：HG008T 模型 → COLO829T 0.708 / 0.754 / 0.630，→ HCC1395 0.790 / 0.698 / 0.689 | 用 `--norm dataset` 训练、LOLO（leave-one-line-out）、learning curve（GPU） |
| R1.8 | 训练集怎么构建 | 部分解决 | 标签规则和代码出处；VAF 分布（d9 和 HG008T d46） | 其余 set 的行（E21 job 被 hold）；Supplementary Table |
| R1.m1 | 补充 benchmark（HCC1395、CASTLE） | 部分解决 | 新 set 已 finalize；HCC1395 的跨样本数字（E12a） | HCC1395 in-sample 模型（GPU）；CASTLE 的结果；比对器 |
| R1.m2 | tensor 窗口在 node 边界截断 | **已解决** | 跨 node 的 tensor 已实现；E1.1 显示模型在用 path-aware 的 ch6 | 只剩示意图和措辞（E24 可选） |
| R2.1 | 为什么纳入 HapMap | 部分解决 | 事实已由 E8b 核实；替代方案已有 | HapMap 去留（决定） |
| R2.2 | truth set 怎么来的 | 部分解决 | 每个 truth 的构建方式按原始文献核实（E8b） | HG008 用 v0.2 还是 v0.3（决定）；COLO829 文件和公开版是否一致 |
| R2.3 | recall 偏低；漏检位点上的 normal | 部分解决 | 漏检的 read-level 分类，以及 normal 组装证据 | 和 DeepSomatic 的差距（等比对器路径） |
| R2.4 | 低 VAF | 部分解决 | purity joint 模型的结果；P75 / P50 的低谷由 LOH 解释 | call 层面按 VAF 分档；低 AF 配置 E15（决定）；BQ 伪影修复 |
| R2.5 | 精度提升从哪来 | 部分解决 | 两步：图减少 germline candidate（E3），ch6 再挡掉一部分（E1.1）。PoN 后剩下的 placed FP 里，135–146 / 150 是 germline（E9d） | 附近 germline 密度分层（E19）；去掉 ch6 重训 |
| R2.6 | PoN | 部分解决 | 资源组成、private allele 测试、truth-in-resources（E8b）；HG008T 不加 PoN 的数字 | 在 validation 上重选 INDEL PoN（E17c，数据同步） |
| R2.7 | 比对统计夸大 | 部分解决 | **E11a**：同一批 Illumina reads 上 perfect 比例 69.12% → 79.36%（d9）；新增的 perfect reads 里 99.74% 位置不变 | d46 和 COLO829T 的 job 被 hold |
| R2.m | Minor 1–8 | 部分解决 | 大部分只是改图改字 | "truth overlapping" 的原意；SBS96（E18） |
| R3.1 | 参考层面的 ablation | 部分解决 | **E3** 三臂 × 12：SNV 上 HPRC 6/6 高于 linear；HG008T 的 linear 和 d9 都用 vg 1.66.0，所以这组对比没有版本混杂 | matched d9；seed / CI；去掉 ch6；complete-output 扩到其余 run |
| R3.2 | 被去掉的 candidate 由什么组成 | 部分解决 | 三类 candidate 对比 × 三臂；E5b；E9d 对 FP 的拆分 | linear 臂 FP 的同样拆分（E4b，要 call sets） |
| R3.3 | graph 坐标 variant 的验证 | 部分解决 | 组装验证有效：placed TP 有 145–149 / 150 是 tumor-only。unplaced calls 里 population 50–89%，germline 8–45% | d9 / Illumina / fiberseq 的 run（数据同步）；rescue 默认值（决定）；COLO829T 没有 tumor 组装 |
| R3.4 | 指标不含 graph 坐标 call；"about 2%"；CHM13 | 部分解决 | "about 2%" 是错的：实际 8–68%。unplaced SNV 全部来自 off-ref rescue；关掉 rescue 后 complete-output F1 ≥ 0.596 / 0.524 / 0.721，高于 linear。CHM13 不做 | rescue 默认值（决定）；其余 9 个 run |
| R3.5 | mutational signature | **待解决** | — | SBS96：要 call sets（数据同步）和比对器 VCF |
| R3.6 | perfect 不等于 correct | **待解决** | 只有 MAPQ 分布和 E11a | E11b 的模拟比对 job 被 hold（96 CPU，约 2 h） |
| R3.7 | 重新归类 FN 要谨慎 | **已解决** | 单倍型组装证据（HG008T、COLO829T） | 只剩改稿 |
| R3.8 | tumor-only 和 mosaicism | **已解决** | 统一写成 tumor-only cancer calling | 只剩改稿 |
| R3.9 | 稳定性、统计不确定性 | **待解决** | 只有间接证据：CPU / GPU 重跑的差 ≤ 0.005；E12a | 多 seed（GPU）；block bootstrap（数据同步） |
| R3.10 | 低 VAF 和 INDEL 的结论 | 部分解决 | SNV / INDEL 分开报告；dilution 结果 | INDEL held-out 阈值；E17a 分层；BQ 伪影修复 |
| R3.m1 | Methods / Supplement 参数不一致 | 部分解决 | 参数只剩一套；从 checkpoint 能复现 RIS 结果 | 等最终模型定了再出参数表 |
| R3.m2 | 比对器设置 | 部分解决 | 版本、资源、训练数据都已核实（E8b）；打分脚本已就绪 | 你已有结果的路径和元数据 |
| R3.m3 | 泛基因组参考的组成 | **已解决** | E8b 核实了组成；E7 的近亲检查：共享的罕见 SNV ≤ 8.2% | 只剩改稿 |
| R3.m4 | 术语和语言 | **已解决** | — | 只剩改稿 |

## 2. 今天的主要结果

1. **三个参考臂（E3；RIS 结果）。**
   - SNV，PoN 后：d9 / d46 在 6 套数据上都比 GRCh38-only 高，幅度 +0.04 到 +0.32。d46 和 d9 基本一样，差在 −0.021 到 +0.031。
   - INDEL：比 test-tuned best F1（上界），HPRC 在 6/6 套上更高。
   - 机制：GRCh38-only 上，样本 75.7–82.2% 的 germline allele 都成了 candidate；HG008T 上 per-tensor somatic AP 只有 0.01–0.02。
2. **unplaced calls（E3 §2、E9a、E9c、E9d；`points/R3.4.md`）。**
   - 最坏情况（unplaced 全算 FP）下，HG008T ONT 和 COLO829T ONT 的 d46 低于 linear。
   - 所有 unplaced SNV PASS 都来自 predict 的 off-reference rescue。
   - 它们能精确放回 GRCh38 的只有 2.8–7.9%。按优先顺序分类：
     - 人群 allele 50–89%；
     - 病人 germline 8–45%；
     - tumor-only 0.2–0.4%；只看组装的话是 1.7–3.8%。
   - 关掉 rescue，complete-output F1 至少是 0.596 / 0.524 / 0.721，linear 是 0.359 / 0.365 / 0.418。
   - 代价：组装推出的 off-ref somatic key 在 tensor 层面本来能叫出 32/112，关掉后是 0/112。这 112 个 key 没有逐个用 reads 核实过。
3. **ch6 counterfactual（E1.1，18 个模型）。**
   - 模型确实在用 path count。把它设成 H 或抹掉信息后，somatic AP 大幅下降，germline 误判增加 3–7 倍。
   - linear 模型的 ch6 是常数，所以不受影响。
   - 同时 ch6 压制 off-reference somatic：没有 rescue 时，组装 key 叫出 0/112、1/89。
4. **跨样本（E12a，初步）。**
   - HG008T 的 d46 模型 → COLO829T：0.708 / 0.754 / 0.630（COLO829T 自己的模型是 0.783 / 0.779 / 0.737）；→ HCC1395：0.790 / 0.698 / 0.689。
   - COLO829T 的模型主要输在阈值迁移上。例如 Illumina → HCC1395 的 recall 只有 0.157，但 best F1 有 0.810。这正是按 dataset 归一化（f2a0e4e）要解决的。
5. **图吸收（E5b，chr1–22）。**
   - 被吸收（ALT 在图里）的 truth 有 95.5–100% 本来就会被 PoN 标掉。
   - 吸收只占图臂比 linear 额外丢掉的 INDEL 的一小部分。例如 HG008T v1.1 丢掉、PoN 又不标的 481 / 581 / 299 个里，被吸收的只有 35 / 34 / 20 个。
6. **图组成（E7，HPRC v1.1）。**
   - 吸收多少由频率阈值决定，单倍型数影响很小：阈值 10% 时，从 16 个单倍型增加到 89 个，HG008 germline SNV 在图里的比例只从 84.5% 升到 86.9%。
   - 和 benchmark donor 共享罕见 SNV 最多的 HPRC 个体，共享比例也 ≤ 8.2%，说明没有近亲。
7. **比对统计（E11a，HG008T）。**
   - 同一批 Illumina reads，GRCh38-only 上 69.12% 是 perfect，d9 上是 79.36%。
   - 新增的 perfect reads 里，99.74% 的 GRCh38 位置不变。
8. **purity mixture（RIS，你的分析）。**
   - 长读长在 P75 / P50 上的低谷由 LOH 区的 germline 解释。
   - Illumina 和 ONT 的 dilution 带 read-source BQ 伪影，分数偏乐观；fiberseq 没有这个问题。
9. **事实核查（E8b）。**
   - ClairS-TO 的训练排除了 chr20；它的 ssrs 模型用了 CASTLE 细胞系。
   - HPRC 的组成；v1.1 里没有 EUR。
   - HapMap donor 在两张图里。
   - 六株 benchmark 细胞系都不在图里。

## 3. 要你做的

### 3.1 GPU（RIS）

| 优先级 | job | 量 | 回应 |
|---|---|---|---|
| P0 | **matched d9 重训**：HG008T PacBio / ONT SNV 用 4 GPU；COLO829T 三套 SNV 用 7fc7396 recipe、4 GPU；COLO829T 三套 INDEL 用 lr 1e-4 | 8 run | R3.1、R1.2、R2.5、R3.2 |
| P0 | **去掉 ch6 重训**：`experiments/E1.1_ch6_counterfactual/gpu_nopc_commands.sh`（`--drop-planes path_count`） | 27 run | R1.1、R2.5、R3.1 |
| P0 | **`--norm dataset` 的跨样本模型**：HG008T 和 COLO829T 的 d46 × 三平台 SNV（加 Illumina INDEL）；HCC1395T in-sample × 3；CASTLE set 完成后做 LOLO | 6–8 + 3 run，之后加 LOLO | R1.7、R1.m1、R3.9 |
| P0 | **多 seed**：三个臂的主 SNV 模型各加 2 个 seed | 约 36 run | R3.9、R3.1 |
| P1 | 修好 BQ 伪影后重训 purity joint 模型，每个 purity 用自己的 validation 阈值（你在做） | 6 run | R1.5、R2.4、R3.10 |
| P1 | learning curve 和模型大小（E13，命令在 `experiments/E21_training_composition/README.md` §6） | 10 run | R1.7 |
| P2 | E9f：用组装给 off-ref tensor 打标签后重训（如果决定做） | 3 run | R3.3、R3.4 |

推理（predict）都在 tequila GPU 上跑，不占 RIS。

### 3.2 数据同步（RIS → 本地）

把除 checkpoint 以外的 run 目录传过来，包括 `test_chr1*/` 的 predictions 和 `vcf_chr1*/`。需要的 run：
- benchmark 里的 36 + 6 个；
- `HG008_*_SNV_h90`、`HG008_*_INDEL_small`、`COLO829_*_SNV_h90`、`COLO829_*_INDEL_small`。

有了这些，CPU 上就能做：
- 其余 run 的 E9a / E9c / E9d；
- INDEL held-out 阈值（validation 的 best-F1 阈值已存在 `metrics.jsonl` 里）；
- SBS96（E18）；
- bootstrap CI（E16）；
- linear 臂 FP 的拆分（E4b）；
- E5b 的 `linear_recovery.md`。

如果不传，就要在 tequila GPU 上重新 predict。

### 3.3 比对器

给我已有 ClairS-TO / DeepSomatic 结果的路径，以及每套结果的：
- 版本；
- 模型；
- 是否 tumor-only；
- BAM / aligner；
- 命令行。

`experiments/E8_comparators/score_comparator.py` 会按 Pansoma 的协议重新打分（同一 truth、BED、chr1、PoN、rtg），只用 CPU，几分钟。

## 4. 要你决定的事

1. **集群上限。** 被 hold 的 CPU job 还在等你定上限，主要是 E11a、E11b 和 E21；IO 重的一次放 1–2 个。我建议审稿相关 job 同时最多 48 CPU、200 GB。
2. **SNV off-ref rescue 要不要保留为默认。** 有三种选择：
   - 保留：unplaced 会拉低 complete-output；
   - 关掉：F1 更高，但 off-ref somatic 叫不出来；
   - 两种都报。

   graph PoN 会标掉 74–88% 的 gold，所以只适合用来描述，不适合当过滤器。
3. **主图选 d9 还是 d46；主评估用跨样本还是 in-sample。** 跨样本的话，in-sample 结果只作上界。
4. **INDEL 的 operating point。** 用 validation 选的阈值（建议），test-tuned best F1 只作上界。
5. **HapMap。** 直接去掉，还是只在 Supplement 里用 donor-free 的 `mgi_non_HG005.d9`。另外要确认原稿用的是哪张图。
6. **下载。**
   - GIAB stratification（R1.3 / E10）；
   - HPRC v2.1 全图 VCF（R1.2、R3.m3）；
   - COLO829T tumor 组装（R3.3）；
   - SEQC2 titration（R2.4）。
7. **vg 版本对照。** 用 vg 1.77 重新比对 linear，还是只在文字里说明。HG008T 的 linear 和 d9 都是 1.66.0。
8. **新功能做不做。** E24（bubble-aware tensor）、E25 / E6（branch candidate）、E9f（off-ref 组装标签）。
9. **低 VAF。** 做 HCC1395 dilution（E14c）还是 SEQC2 titration（E14d）；E15 低 AF 配置做不做。
10. **truth。**
    - HG008 用 v0.2 还是 v0.3；
    - rule (a) 是否也用到 HG008T / COLO829T；
    - H1437 的异常怎么处理：PacBio somatic SNV 只有 70.9% 有 tensor，truth 里约 30% 被 PoN 标掉。
11. **术语。** "panel of normals" 是否改成 "population filtering resources"（R3.m4）。

## 5. 实验文件夹（`experiments/`）

| ID | 内容 | 状态 | 回应 |
|---|---|---|---|
| E3_reference_arms | 三个参考臂的表；complete-output 最坏情况 | 完成 | R3.1、R1.1、R1.2、R2.5 |
| E1.1_ch6_counterfactual | ch6 counterfactual，18 个模型（tequila GPU） | 完成 | R1.1、R1.m2、R2.5、R3.3 |
| E12a_cross_sample | 16 个跨样本 predict（tequila GPU），初步 | 完成 | R1.7、R1.m1、R3.9 |
| E9a_unplaced_projection | unplaced → GRCh38 投影；complete-output | 3 个 d46 run 完成 | R3.4、R1.3 |
| E9c_offref_filter_validation | graph PoN 和组装验证 | 3 个 d46 run，加 4 个 legacy run | R3.3、R3.4、R1.3、R2.5 |
| E5b_absorbed_d46 | d9 / d46 / linear 的吸收分类 | 完成（`linear_recovery.md` 要数据同步） | R1.4、R3.2 |
| E7_graph_composition | in-silico 图组成；近亲检查 | v1.1 完成（v2.1 要下载） | R1.2、R3.m3 |
| E11a_alignment_stats | 同一批 reads 在四个参考上的统计 | HG008T 部分完成；其余被 hold | R2.7、R3.6、R1.3 |
| E11b_simulated_alignment | 模拟 reads 的比对正确性 | 模拟完成；比对被 hold；CHM13 不做 | R3.6、R1.3 |
| E21_training_composition | 训练集组成；E13 的准备 | 主表完成；其余被 hold | R1.7、R1.8 |
| E8b_literature_facts | 文献和资源事实 | 完成（HapMap 那项被 hold） | 多条 |
| E8_comparators | 已取消（你已有结果）；打分脚本保留 | 取消 | R3.m2 |

要清理的数据：E8 中途留下的 bwa / chr1 BAM 在 `/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators`，约 5.7 GB，删不删由你定。

## 6. 逐条回应

每条的格式：状态 / Comment（要点）/ Response（英文）/ 证据（内部）/ 还缺什么。

---

## 6.1 Reviewer 1

### R1.1 "Pangenome-native" is not supported

**状态：部分解决**: E1.1 已在全部 18 个 SNV 模型上完成（12 个 HPRC 图模型 + 6 个 GRCh38-only 对照），证明分类器在用 ch6（这是 "在用"，不是 "有益"）；E3 给出同 reads、同 tensor 代码下 HPRC 相对 GRCh38-only 的 end-to-end 增益（BED 内、unplaced 不计分）。还缺：去 ch6 重训臂（GPU RIS），用它把 "candidate 定义" 和 "ch6 通道" 的贡献分开；臂间差值的 bootstrap CI；以及组装 off-ref key 的 read-level 核实。

**Comment（要点）。**

- The pangenome is used mainly at vg Giraffe alignment; the downstream representation is a local pileup, so the method resembles DeepVariant, DeepSomatic or ClairS.
- Neither the input nor the network models graph edges, bubbles / superbubbles, cross-node context or long-range haplotype structure.
- Replacing linear coordinates with node IDs is not enough; haplotype and path information is discarded.
- Request: substantial additional evidence, or a re-design, showing measurable benefit from pangenome-specific information beyond alignment.

**Response.**

We agree that "fully exploit the richness of pangenome data" overstated the design. It has been removed, and Pansoma is now described as pangenome-aware. We also agree that node IDs alone would not make a method pangenome-native. Graph information enters at three stages after alignment, and the revision quantifies each.

1. *Candidate definition.* Candidates are edits relative to the nodes each read traverses, so alleles on the graph's haplotypes produce no candidate. On identical reads, 76–82% of matched-normal germline alleles give a candidate on a GRCh38-only graph, against 7–13% on HPRC v1.1 and v2.1 (six datasets; R1.2).
2. *Path-following tensors.* Each tensor row follows its read's walk across node boundaries. Two channels come from the graph: the base of the node (channel 5) and the number of haplotype paths through the node (channel 6), coded relative to the graph's haplotype count.
3. *The classifier uses channel 6.* We edited channel 6 at inference, without retraining, in the SNV models of all six datasets on both HPRC graphs (12 models, chr1; COLO829T on v1.1: the four-purity joint models; new Supplementary Table). Setting it to an uninformative constant (the training mean) or to the all-haplotype count lowers the somatic average precision from 0.199–0.352 to 0.044–0.087 (HG008T) and from 0.639–0.829 to 0.329–0.520 (COLO829T). With the all-haplotype count, 3.4–7.6 times as many germline tensors pass the calling threshold. Channel 6 alone separates somatic from germline tensors (AUC 0.84–0.89, three HG008T sets). On the GRCh38-only graph, where it is constant, the edits change nothing.

The population filter removes most of these extra calls (all-haplotype count: chr1 F1 after filtering, estimated per tensor, changes by −0.063 to +0.007), so the two carry overlapping information. Channel 6 also works against somatic SNVs on off-reference nodes: of the HG008T truth SNVs whose allele lies on an off-reference v1.1 node (placed by aligning the donor's assembly) and has a tensor, the model alone calls 0 of 112 (PacBio) and 1 of 89 (ONT). The pipeline therefore re-scores off-reference tensors with the site counts raised to the all-haplotype code (R3.4), which calls 32 and 46 of them.

End to end, the reference-level ablation (R3.1) measures the benefit. With identical reads and tensor code, both HPRC graphs beat the GRCh38-only graph on all six datasets: by +0.038 to +0.319 in chr1 SNV F1 after population filtering (in the truth BED, unplaced calls not scored; R3.1) and in indel F1 at each curve's test-tuned best point. This contrast combines alignment, candidate definition and channel 6; a model retrained without channel 6 separates the last. [待补: no-channel-6 arm on HPRC v1.1 and v2.1]

Pansoma does not model bubbles, competing branches or long-range haplotypes explicitly, and the revised text says so. [待补, optional: E24 bubble-aware tensor]

**证据（内部）。**

- 76–82 % / 7–13 %：`response_to_reviewers_detailed.md` R1.2 表（2026-10-12，各套 `tensors/truth_recall.json` 的 `tensor_representative` 比例）。Matched-normal germline：linear 0.757–0.822，v1.1 d9 0.077–0.131，v2.1 d46 0.072–0.124。简版写的 "8–12% on HPRC v1.1" 和这张表不一致（v1.1 是 0.077–0.131），以详版表为准。
- ch5 / ch6 定义：`indexed_gam_pipeline_v4/README.md` channel 表（ch5 graph base；ch6 path count，按 `--haplotypes` H 编码，H → 100）。各套 H：linear 17、d9 90、d46 464（`pansoma_v2_tensors/HG008T_PacBio_d46/tensors/SNV/manifest.json` `"haplotypes": 464`）。
- E1.1：`experiments/E1.1_ch6_counterfactual/results/tables.md` 和 `summary.json`（collect.py 2026-10-13 15:32；`summary.json` 的 `done` 有 18 个 run，`pending` 为空）。
  - 18 个 run：HG008T d9 × 3、HG008T d46 × 3、COLO829T d46 × 3、COLO829T d9 × 3，加 6 个 linear。COLO829T 的 d9 是 purity-joint 模型 `COLO829_*_SNV_joint_7fc7396_gpu4`（在 P100 set 上测，PASS 阈值是四个 purity 合起来定的），不是 E3 §1 里的 `COLO829_*_SNV_h90`。所以 E1.1 里 COLO829T d9 的数字和 E3 的 d9 数字不能直接对照。Response 原来写 "the final SNV models"，对 COLO829T d9 不准确，已改成注明 joint 模型。
  - joint 模型的 Illumina / ONT 训练数据包括 P75 / P50 / P25 dilution，这两个平台的 dilution 有 read-source BQ artefact（用户 2026-10-13 的判断；fiberseq 干净）。E1.1 是在 P100 上、同一个模型的各 mode 之间比较，所以结论（方向和倍数）不受影响；但 `COLO829_Illumina_joint` / `COLO829_ONT_joint` 的绝对 AP / F1 不要拿去和别的模型比。
  - 各 mode 的 somatic AP（§4）：
    - HG008T base 0.1986 / 0.2333 / 0.2316（d9 Illumina / PacBio / ONT）、0.3121 / 0.3523 / 0.3264（d46）。`all_H` 0.0712 / 0.0870 / 0.0706、0.0621 / 0.0723 / 0.0612。`train_mean` 0.0601 / 0.0525 / 0.0437、0.0527 / 0.0661 / 0.0513。合起来 0.044–0.087。
    - COLO829T base 0.8199 / 0.8291 / 0.7546（d46 Illumina / fiberseq / ONT）、0.6926 / 0.7135 / 0.6392（d9 joint）。`all_H` 0.5158 / 0.5149 / 0.4049、0.4866 / 0.5149 / 0.3472。`train_mean` 0.5046 / 0.5197 / 0.4213、0.4817 / 0.5116 / 0.3293。合起来 0.329–0.520。
    - `shuffle`：HG008T 0.0277–0.0428，COLO829T 0.1921–0.3079。
  - germline ≥ t 相对 base 的倍数（`all_H`）：HG008T d9 3.58 / 3.50 / 3.79，d46 6.75 / 6.83 / 6.64；COLO829T d46 7.57 / 6.72 / 7.39，d9 joint 3.45 / 3.41 / 4.19。范围 3.41–7.57。
    - `train_mean` 的倍数在 HG008T ONT d9 上是例外（0.27）：这个模型在 `train_mean` 下几乎什么都不叫（somatic 65 / 663 ≥ t，PoN 后 F1 0.145）。所以 Response 里的倍数只引 `all_H`。
  - `site0_low`（只改 site 列，设为 off-reference test tensor 的 site code 中位数 LOW = 22–36）：somatic ≥ t 在 12 个图模型上是 0–56 个（56 是 COLO829T ONT d9 joint，1,984 个里；HG008T 六个是 0 / 0 / 0 / 1 / 0 / 0）。
  - per-tensor PoN 后 F1，base → `all_H`：HG008T d9 0.589 → 0.573、0.621 → 0.566、0.596 → 0.533；d46 0.560 → 0.511、0.605 → 0.561、0.579 → 0.516；COLO829T d46 0.784 → 0.791、0.785 → 0.760、0.743 → 0.705；d9 joint 0.675 → 0.682、0.782 → 0.778、0.706 → 0.664。差值 −0.063 到 +0.007。`shuffle` 损失大得多（PoN 后 0.212–0.551）。
  - linear 6 个 run：8 个 mode 的所有数字和 base 完全一样（ch6 恒为 100，LOW = 100）。
  - §5 组装 off-ref somatic SNV（HG008T，d9，全部常染色体，都不是训练标签；key 来自 `pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv`，即 HG008-N 组装上的 truth somatic 事件，giraffe 放到 d9 的 off-reference 节点）。格式是 "≥ t / argmax somatic / n"：
    - base：PacBio 0 / 0 / 112，ONT 1 / 1 / 89。
    - `all_H`：61 / 44 / 112，49 / 49 / 89。
    - `site10_high`：45 / 32 / 112，45 / 46 / 89。按 E1.1 README，`site10_high` 就是 predict 的 off-reference rescue 变换；predict 对 off-ref tensor 用 argmax（`predict.py` 第 31–39 行，`graph_vcf.py` 的 `offref_call`），所以 rescue 后是 32 / 112 和 46 / 89。
    - d46 的组装 key（`asm_d46_classify.py`）没有进 §5，表里只有两个 d9 run。
    - 这些 tensor 不一定带 somatic reads：旧的 read-level 检查（`pansoma_net_v2_runs/offref_negatives_20260929/overrides.tsv`，HG008T Illumina d9）里，exact key 的组装 SNV tensor 只有 30 个（17 train + 13 held out）的 A1 reads 只出现在 HG008-T 组装里，123 个在 HG008-N 组装里也有（"normal"），33 个两边都没有（"neither"）。PacBio / ONT 的 112 / 89 个没做这个检查。所以 Response 里不再写 "true somatic" 和 "assembly-validated"，改成 "truth SNVs whose allele lies on an off-reference node ... and has a tensor"；rescue 后 "calls 32 and 46" 里有多少是真的 somatic read 证据，还不知道。
  - §6 只看 ch6（`profile_auc.json`，不用模型，只跑了 3 个 run）：site code 的 AUC 0.837（HG008T PacBio d9）、0.890（PacBio d46）、0.892（ONT d46）。site code 中位数 somatic 100，germline 96（d9）/ 79（d46）。
  - §2：E1.1 的 base mode 对 RIS 的 GPU bf16 `predict`（同 checkpoint、同 test tensor、同阈值），somatic AP 差 ≤ 0.0047（最大是 COLO829T fiberseq linear 0.2832 / 0.2879，其次 COLO829T Illumina joint 0.6926 / 0.6971）。注意 E1.1 的 base 只有 4 个 run 是 CPU fp32（`HG008_PacBio`、`HG008_ONT`、`HG008_PacBio_d46`、`HG008_ONT_d46`，AP 0.2333 / 0.2331、0.2316 / 0.2323、0.3523 / 0.3528、0.3264 / 0.3283），其余 14 个是 tequila GPU bf16（E1.1 README "How the runs were done" 第 3 条；`collect.py` 的 §2 标题 "CPU fp32 here" 不准，README 已注明）。
  - §2 per-tensor PoN 和 rtg：图模型上 TP 差 ≤ 8，FP 在 12 个里有 11 个差 ≤ 52（例外：COLO829T ONT d9 joint 399 对 595）。linear 上 rtg 的 FP 多很多：HG008T 837 / 713 / 692 对 1,091 / 973 / 904，COLO829T Illumina 211 对 2,216，ONT 1,063 对 2,575。
    - 可能的原因（未核实）：per-tensor 评估只用 confident region = somatic BED ∩ germline BED（`machine_learning/pansoma_net_v2/data.py` 的 `eval_label`；`COLO829T_Illumina_linear/tensors/SNV/labels.manifest.json` 的 `confident_region`）。rtg 只用 somatic BED。linear 上 germline BED 外的 germline candidate 在 per-tensor 里被排除，在 rtg 里算 FP。
    - 所以 E1.1 的 "after PoN" F1 只能在同一个模型的各 mode 之间比较。臂间比较（R3.1）只用 rtg 的数字（E3）。
  - 早上从 log 里取的数字（base AP 0.2333 / 0.3523，`site0_low` 中位 p 0.016 / 0.003，`all_H` germline 17,327 / 17,256，`shuffle` AP 0.0321 / 0.0277）和 `tables.md` 一致。
- End-to-end：`experiments/E3_reference_arms/tables.md` §1（PoN 后、BED 内、PASS；rtg 不算 unplaced call）。SNV d9 − linear 是 +0.059 到 +0.311，d46 − linear 是 +0.038 到 +0.319。no BED、unplaced 全算 FP 的 worst-case bound（E3 §2）下，HG008T ONT 和 COLO829T fiberseq（d9）/ ONT（d46）不再高于 linear，见 R3.1；所以 Response 里加了 "in the truth BED, unplaced calls not scored"。"identical recipe" 只对 linear 和 d46 完全成立；12 个 d9 模型里有 8 个训练设置不同（E3 README Caveats 1），Response 改成 "identical reads and tensor code"，细节放 R3.1。INDEL best F1：linear 0.199 / 0.173 / 0.082 / 0.069 / 0.092 / 0.019，d9 0.264 / 0.260 / 0.140 / 0.195 / 0.120 / 0.071，d46 0.249 / 0.269 / 0.155 / 0.163 / 0.129 / 0.085。
- 两份来源不一致：E3 README（10-13 10:17）按 PASS F1 写 "INDEL: no pangenome advantage after the PoN, except on Illumina"；benchmark README（commit 0049e36，history 改写前是 d4f898c）改为按 best F1 比，linear 在全部 6 套 INDEL 上都低。这里按较新的 benchmark README 写，PASS F1 的差异放在 R3.1 里说明。INDEL best F1 是 test-tuned 上界，held-out 规则还没定。
- "ch6 和 PoN 信息重叠" 还有一条旧证据：详版 R1.1 评估里的 nopc ablation（`analysis/v1_vs_v2_20260929`），旧 build 去掉 ch6 后 PoN 后打平。这次 Response 改用 E1.1 的新数字，不再引旧 build。

**还缺什么。**

- [GPU RIS] 去 ch6 重训臂（E3 / E1.1）：`experiments/E1.1_ch6_counterfactual/gpu_nopc_commands.sh`，共 27 个 run：d9 和 d46 各 12 个 no-ch6，加 3 个 COLO829T d9 INDEL 的 with-ch6 对照。出了以后替换 Response 里的 [待补]。
- [数据同步] + [CPU 本地] E16：各臂 F1 差值的 chr1 block bootstrap CI。要等 RIS 的 call sets 同步回本地。
- [CPU 本地]（可选）`profile_auc.py` 只跑了 3 个 HG008T run。如果 Response 要写 "all models"，补跑其余 9 个图 run（只读 tensor，不用 GPU）。
- [用户决定]（可选）d46 的组装 off-ref key 打分（`asm_scores.py`，tequila 上很短的 GPU inference，E1.1 README "Remaining"；没有提交 job），以及 HG008T Illumina d9 的组装 key（`--assembly`）。
- [CPU 本地] §5 那 112 / 89 个 PacBio / ONT 组装 key 的 read-level 检查（A1 reads 只在 HG008-T 组装里，还是 HG008-N 里也有；方法同 `offref_negatives_20260929`）。做完才能在 Response 里说 rescue 叫回来的是真的 somatic。
- [CPU 本地] 核实 linear 臂上 per-tensor PoN 和 rtg 的 FP 差异（上面的 germline BED 假设）。这影响的是 E1.1 内部 "after PoN" 数字的解读，不影响 Response 里的 rtg 数字。
- [文字] E1.1 README 现在已有 Results / Status 段（15:45 版）。新的 Supplementary Table 用 `tables.md` §4（AP、germline ≥ t、PoN 后 F1，注明 PoN 后 F1 是 per-tensor 估计、不是 rtg）。`collect.py` §2 的标题 "CPU fp32 here" 要改。
- [用户决定] 是否做 E24（locus / bubble-aware tensor，改 pipeline，约 2–3 周）。不做的话，Response 最后一句的 [待补] 删掉，只保留 "does not model bubbles"。
- [文字] Methods / Abstract：删 "fully exploit"，"pangenome-native" 改为 "pangenome-aware"，写清 ch5 / ch6 的定义和 H 的编码，以及 off-reference rescue 的规则（和 R3.4 一致）。

### R1.2 Composition and size of the pangenome

**状态：部分解决**: d9 / d46 的对比已经齐了（candidate 层面六套，训练后 chr1 十二套），E1.1 也给出了 haplotype 数对 ch6 和 germline 识别的影响，E7 的 in-silico 曲线和用 v1.1 实测吸收做的校准（s4）都已出。还缺：mapping ambiguity（E11a 没出）、E5b graph spelling 已完成（absorbed 只占 graph 臂比 linear 多丢的 INDEL 的一小部分，见 R1.4），剩下的损失原因（not spelled / filtered）还要解释、v2.1 单倍型轴（要下载 VCF，用户决定），以及 matched d9 重训（d46 ≈ d9 的结论目前只有 HG008T Illumina SNV 和 HG008T 三套 INDEL 是只差图）。

**Comment（要点）。**

- No systematic analysis of how the number and diversity of haplotypes, haplotype walks and alternative paths affect node pileups and germline identification.
- The graph decides which alleles are represented, and so which variants are suppressed before the classifier.
- With HPRC2 (many more haplotypes): better population representation, but possibly more mapping ambiguity and more absorbed somatic mutations that overlap population alleles.

**Response.**

We added HPRC v2.1 as a third reference. The reads of all six datasets (HG008T and COLO829T, three platforms each) were aligned to three references and rebuilt with identical tensor code (aligner versions in R3.1):

- a GRCh38-only graph;
- HPRC v1.1 (90 haplotypes, nodes kept when on at least 9);
- HPRC v2.1 (464 haplotypes, including 29 European individuals; v1.1 has none; nodes kept when on at least 46).

Supplementary Table X gives the fraction of truth alleles that reach a candidate tensor. Against v1.1, v2.1 changed it by −0.001 to +0.004 for somatic SNVs, +0.005 to +0.043 for somatic insertions, −0.004 to +0.025 for somatic deletions and −0.012 to −0.005 for matched-normal germline alleles. The five-fold larger set thus suppressed germline candidates slightly better without measurably absorbing more somatic alleles. The haplotype-count channel, which codes how many haplotypes miss a node (100 = none), gives germline tensors lower codes on v2.1 (median 74–80 against 95–96 on v1.1; about 4–6% of the haplotypes in both) and alone separates them from somatic tensors better (AUC 0.890 against 0.837, HG008T PacBio; R1.1). After population filtering, models trained on each graph are equivalent on chr1 (SNV F1 −0.021 to +0.031; indel F1 at the test-tuned best point −0.033 to +0.015). [待补: matched v1.1 retraining]

Absorption depends on which alleles the graph contains, and so on the frequency filter more than on the haplotype count. We measured both axes in silico, with random subsets of the HPRC v1.1 haplotypes:

- *One carrier enough:* the share of HG008 germline SNVs in the graph saturates (0.851 with 8 haplotypes, 0.940 with 32, 0.960 with 89), while the share of HG008T somatic indels keeps growing (0.306, 0.490, 0.590; 0.709 for homopolymers of 7 bp or more).
- *At the 10% filter of both HPRC graphs:* from 16 to 89 haplotypes, somatic indel representation stays at 0.231–0.259 and germline SNV representation at 0.845–0.869.
- *Composition:* 42 non-African haplotypes contain more of HG008's germline SNVs than 42 African ones (0.950 against 0.917), and fewer of its somatic indels (0.462 against 0.521).

On the real v1.1 graph the prediction holds for germline SNVs: 0.87 are predicted in the graph, and 0.89–0.90 lose their candidate relative to GRCh38-only. For somatic indels it is a lower bound (HG008T: 0.26 predicted, 0.38–0.49 lost), because 0.29–0.36 of the indels the graph does not contain also lose their candidate. [待补: cause of this extra loss (R1.4); v2.1 haplotype subsets]

Counting only truth alleles that population filtering would keep, the net loss of somatic candidates on either HPRC graph relative to GRCh38-only is 0.2–0.6% of SNV truth and 2.0–8.2% of indel truth (chr1–22; R1.4). Mapping ambiguity is addressed in R1.3 and R3.6. [待补: MAPQ of the same reads on the three references, E11a]

**证据（内部）。**

- 四个变化范围是用详版 R1.2 表（`response_to_reviewers_detailed.md`，2026-10-12，各套 `tensors/truth_recall.json`）逐行相减得到的：
  - SNV：+0.004 / +0.002 / +0.003 / 0 / 0 / −0.001。
  - INS：+0.038 / +0.028 / +0.005 / +0.042 / +0.043 / +0.014。
  - DEL：+0.021 / +0.008 / −0.002 / +0.025 / +0.011 / −0.004。
  - germline：−0.005 / −0.012 / −0.009 / −0.006 / −0.009 / −0.007。
  - germline 原来写成 "−0.005 to −0.012"，顺序改成 −0.012 to −0.005。
  - 详版原文说 "v2.1 kept as many or more somatic alleles than v1.1"，这不准确：COLO829T ONT 的 SNV 0.893 → 0.892，HG008T ONT 的 DEL 0.526 → 0.524，COLO829T ONT 的 DEL 0.679 → 0.675 都略降。所以这里改写成范围。
- 29 个 EUR、v1.1 EUR 0：`experiments/E8b_literature_facts/README.md` H2、H5。
- ch6 随图变大更有信息（新，E1.1）：
  - `experiments/E1.1_ch6_counterfactual/results/tables.md` §1 和 §6（`profile_auc.json`）。germline test tensor 的 site code 中位数：HG008T d9 96 / 96 / 96，d46 80 / 79 / 79；COLO829T d9 joint 96 / 95 / 95，d46 80 / 76 / 74。somatic 全是 100。
  - 编码（`indexed_gam_pipeline_v4/candidates.py` `encode_count`）：H−49 ≤ count < H 时 code = 100 − (H − count)，即 "少了几条 haplotype"。所以 d9 的 95–96 = 90 条里少 4–5 条（4.4–5.6 %），d46 的 74–80 = 464 条里少 20–26 条（4.3–5.6 %）：按比例两张图一样，code 更低是因为 H 大了五倍、少的条数多了。早上的 Response 写 "germline tensors sit on nodes carried by fewer haplotypes"，按比例不对，已改成 "codes how many haplotypes miss a node ... about 4–6% of the haplotypes in both"；AUC 只有 HG008T PacBio 这一对 d9 / d46。
  - 只用 site code 区分 somatic / germline 的 AUC：HG008T PacBio d9 0.837，d46 0.890；HG008T ONT d46 0.892（d9 ONT 没算，`profile_auc.py` 只跑了 3 个 run）。
  - 模型对 ch6 的依赖（没写进 Response，因为 d9 / d46 的训练设置不完全对齐）：`all_H` 时 germline ≥ t 的倍数，d9 是 3.41–4.19（HG008T 3.58 / 3.50 / 3.79，COLO829T joint 3.45 / 3.41 / 4.19），d46 是 6.64–7.57。base 的 germline ≥ t，HG008T d9 4,864 / 4,951 / 5,147，d46 2,542 / 2,528 / 2,810。
  - 注意：E1.1 的 COLO829T d9 是 purity-joint 模型（P100 上测），HG008T PacBio / ONT d9 是 1 GPU；只有 HG008T Illumina 这一对是只差图。
- 训练后 d46 − d9：
  - SNV（PASS F1，PoN 后，BED 内）：−0.015 / −0.016 / −0.021 / +0.005 / −0.002 / +0.031（`E3_reference_arms/tables.md` §1）。
  - INDEL best F1：−0.015 / +0.009 / +0.015 / −0.033 / +0.009 / +0.013（`analysis/benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/README.md`）。INDEL best F1 是 test-tuned 上界，held-out 规则还没定。
  - raw SNV PASS F1 提高 0.10–0.13：同一 README，例如 HG008T PacBio 0.191 → 0.310。chr1 SNV tensor 数 HG008T PacBio 32,880 → 27,501。这一句从 Response 里删了（字数），放在 R3.1。
  - 注意：d46 − d9 只有在 HG008T Illumina SNV 和 HG008T 三套 INDEL 上是 "只差图"。其余几套的 d9 run 在 global batch（1 GPU 对 4 GPU）或 lr（COLO829T INDEL 用的 4e-4）上也不同。仅 1 → 4 GPU 这一项，就让 HG008T Illumina d9 降了 0.007。
- E7 曲线：`experiments/E7_graph_composition/curves.tsv`，scope `bed`，pool `all`（89 条非 GRCh38 单倍型）、`AFR` 和 `nonAFR`，`threshold_type rel`，每个 k 抽 20 次。源文件是 6 位小数，这里四舍五入到 3 位。
  - HG008 germline SNV，阈值 ≥ 1 条（rel 0.0）：k8 0.850796，k32 0.940044，k89 0.959909。
  - HG008 somatic INDEL，同阈值：k8 0.305994，k32 0.489827，k89 0.590319。INDEL_HP>=7 在 k89 是 0.709247。
  - rel 0.1（k89 时 d = 9，即 d9）：somatic INDEL k16 0.25924、k32 0.231243、k64 0.245172、k89 0.255893；germline SNV k16 0.844506、k89 0.869412。
  - k = 42、rel 0.0：germline SNV nonAFR 0.949734 对 AFR 0.916904；somatic INDEL nonAFR 0.462498 对 AFR 0.520572。AFR 是 46 条里抽 42 条、20 次抽样的均值；nonAFR 一共只有 42 条，所以是全部，`n_draws` = 1（没有抽样误差估计）。
  - somatic SNV 即使阈值为 1 条也很低：HG008 0.046003，COLO829 0.007208（k89）。
  - 图大小（`graph_size.tsv`，chr1-22，k89）：SNV 在阈值 1 条时 21,076,950 个，在 d9 时 6,869,974 个。
- E7 校准（s4，`experiments/E7_graph_composition/calibration.tsv`，2026-10-13 11:56 写出；早上的版本说 "还没产出"，不对）。scope `bed`，chr1-22。`predicted_in_d9` = 89 条单倍型里携带数 ≥ 9；`loss_d9_vs_linear` = linear 有 tensor、d9 没有的比例（分母是 linear 有 tensor 的 allele）：
  - germline SNV：predicted 0.8694（HG008）/ 0.8724（COLO829）；loss 0.8979–0.9019（HG008 三平台）/ 0.8908–0.9000（COLO829）；predicted 的那些里 loss 0.9929–0.9949 / 0.9881–0.9952。d46 的 loss 0.9059–0.9097（HG008）。
  - somatic INDEL：predicted 0.2559（HG008）/ 0.2228（COLO829）；loss HG008 0.4921 / 0.4519 / 0.3754（Illumina / PacBio / ONT），COLO829 0.3812 / 0.3485 / 0.3081；不在预测里的（`loss_d9_if_pred_out`）HG008 0.3195 / 0.3603 / 0.286，COLO829 0.2235 / 0.2406 / 0.2151。
  - 分层：HG008 somatic INDEL_STR 的 loss 0.7263–0.7638（predicted 0.2733），INDEL_HP>=7 0.3353–0.5006（predicted 0.311）。
  - somatic SNV：predicted 0.024（HG008）/ 0.0019（COLO829），loss 0.0333–0.0372 / 0.0054–0.0103。
  - E7 还没有 README，`s7_tables.py` 还没跑。Response 里的数字直接取自 `calibration.tsv`。"lose their candidate" 指没有 `tensor_representative`，partial 匹配不算。
- 净损失 0.2–0.6 % / 2.0–8.2 %：`experiments/E5b_absorbed_d46/tables.md` §1 最后一列（"tensor and not tagged"，linear − arm，按 truth 的百分比，chr1-22）。这是 tensor 层面的上限（没被 PoN tag、又有 tensor 的 truth），不是 call 层面的损失；Response 已改成 "net loss of somatic candidates ... that population filtering would keep"。
  - d9：SNV −0.2 到 −0.6 %，INDEL −2.0 到 −8.2 %。d46：SNV −0.2 到 −0.6 %，INDEL −2.0 到 −8.1 %。
  - 这张表标的是 PRELIMINARY，但 §1 不依赖 graph spelling，只有 "absorbed" 那几列依赖。

**还缺什么。**

- [CPU 本地] E7：跑 `s7_tables.py`、写 README（s4 / s6 的输出已经在了）。之后核对 Response 里的校准数字。
- [CPU 本地] E5b graph spelling（`E5b_s3_extract_d46` 387428 在跑；`E5b_s4_spell_d9` / `_d46` 387430 / 387429 等 dependency，2026-10-13 下午 squeue）。用来解释 "图里没有的 somatic INDEL 也有 0.29–0.36 丢了 candidate"，替换 Response 第一个 [待补]；完成后 R1.4 的吸收分类才是最终版。
- [用户决定] E7 的 v2.1 轴：是否下载 HPRC v2.1 全图 VCF（vcfbub / wave）。不下载的话，删掉 "v2.1 haplotype subsets"，只报 v1.1 子集，加上 d9 / d46 两个实测点。
- [CPU 本地] E11a MAPQ / mapping ambiguity（`E11a_win` 387287、`E11a_scan` 387288 还被 hold 着，2026-10-13 下午 squeue），用来替换第二个 [待补]。
- [GPU RIS] matched d9 重训：HG008T PacBio / ONT SNV 用 4 GPU，COLO829T d9 SNV 在 7fc7396 重训，COLO829T d9 INDEL 用 lr 1e-4。这样 d46 − d9 在 12 套上都只差图，E1.1 的 d9 / d46 "ch6 依赖" 对比也才干净。
- [CPU 本地]（可选）`profile_auc.py` 补跑 d9 ONT / Illumina 和 COLO829T，让 AUC 的 d9 / d46 对比不只靠 HG008T PacBio 一对。
- [用户决定] 主图选 d9 还是 d46（详版 §1.3）。结果打平时，这决定正文报哪张图。
- [文字] 把详版 R1.2 的六套 × 三臂 candidate 表做成 Supplementary Table，并改掉 "as many or more" 的说法。

### R1.3 Difficult regions and ONT ambiguity

**状态：部分解决**: "ONT graph-coordinate calls 是 resolution 还是 ambiguity" 这一半现在有最终数据，但只来自三个 d46 in-sample SNV run 的 chr1 re-prediction（HG008T PacBio 在 CPU 上跑，fp32；HG008T ONT / COLO829T ONT 在 tequila GPU 上跑，bf16）。E9a：unplaced PASS 只有 7.9 / 3.4 / 2.8 % 能写成 GRCh38 variant（class i），PoN 后一个 TP 都没多（PoN 前两个 HG008T run 各多 1 个，是同一个 truth SNV，PoN 会标掉）；class ii（只有 interval）的 221 / 984 / 1,334 个没有匹配到任何漏掉的 truth（FN）。E9c / E9d 的类别互斥，按 population > germline > tumor-only > unverifiable > unsupported 的优先顺序：population（graph PoN）89.2 / 74.1 / 49.8 %，病人 germline 7.9 / 19.2 / 45.4 %，tumor-only 只有 0.4 / 0.2 %（COLO829T 只处理了随机 3,000 / 4,569 个，not in normal 2.4 %）。不看优先顺序时，组装判 tumor-only 的有 9 / 22 个（3.8 / 1.7 %），其中 8 / 19 个被 graph PoN 归成 population。对照组说明组装分类在 placed call 上可靠。但 graph PoN 对 somatic 不特异：d46 gold（全常染色体）被它标掉 90 / 119、72 / 97。chr1 上的 gold 也只有 11 / 9 / 2 个。所以第 5 条按 "added ambiguity" 写，限定在这三个 run 上。还缺：审稿人明确要的按区域分层 calling 指标（E10）还没有数字，也还没有脚本（三个 d46 re-prediction 的 rtg 输出在本地，linear / d9 臂的在 RIS 上）；E11b（模拟 reads 的 MAPQ 校准）的 alignment job 仍 held；rescue 是否保留为默认要用户定。

**Comment（要点）。**

- 在重复、paralog、结构复杂和低 mappability 区域，pangenome 比对仍可能有歧义，ONT 这类高噪声 reads 尤其如此；haplotype 越多，竞争比对位置可能越多。
- 汇总统计（更多 perfect reads、更少 mismatch / indel / soft clip）不能说明改善发生在 pangenome 应该起作用的地方。
- ONT 的 graph-coordinate predictions 更多：是分辨率提高，还是额外的歧义？
- 需要按区域的 benchmark 和错误分析。

**Response.**

We agree and added node-, region- and call-level analyses.

1. *Where alignments change.* On HPRC v1.1, 4.9–8.7% of node visits by MAPQ > 5 alignments are on non-reference nodes, 94.3–98.1% of them edit-free. Per node-length bin, long reads are about as often edit-free there as on GRCh38 nodes (HG008T PacBio, 1 bp: 98.5% vs 97.9%), Illumina less so (6–32 bp: 88.0% vs 94.5%). The gain is edit-free alignment to population alleles (R3.6).

2. *Competing mappings.* Of 14,014,144 HG008T PacBio reads, 95.88% of mapped reads have MAPQ 60 on the GRCh38-only graph and on HPRC v1.1, 83.32% on v2.1 (464 haplotypes, newer giraffe), while reads with MAPQ > 5 barely change (13.34–13.37 million) [待补: simulated-read MAPQ calibration (E11b)]. Chr1 SNV F1 after filtering differs by −0.021 to +0.031 between v1.1 and v2.1.

3. *ONT in difficult regions.* Benchmark SNVs with neither a candidate nor any MAPQ > 10 read number 1 on HG008T ONT against 48 (PacBio) and 55 (Illumina), and 45 on COLO829T ONT against 227 (Fiber-seq) and 237 (Illumina), mostly in SMaHT "extreme" regions; there ONT leaves 4.7% of COLO829T truth SNVs without a candidate, against 7.9–8.8%.

4. *Region-stratified calling.* [待补: P / R / F1 per SMaHT / difficult-region stratum and reference arm, with bootstrap intervals (E10)]

5. *ONT graph-coordinate calls.* ONT makes the most PASS SNVs without a GRCh38 position (chr1, HG008T v1.1 / v2.1: 1,409 / 1,261, against Illumina 116 / 754 and PacBio 211 / 244; COLO829T: 1,675 / 4,566), all from the off-reference rescue (R3.4). In re-predictions of three v2.1 runs (HG008T PacBio, HG008T ONT, COLO829T ONT; 240 / 1,258 / 4,569 calls):
   - 7.9 / 3.4 / 2.8% can be written as a GRCh38 variant, adding no true positive after filtering; 221 / 984 / 1,334 have only a GRCh38 interval, matching no missed truth SNV; 0 / 231 / 3,107 have no anchor.
   - In exclusive priority order, 89.2 / 74.1 / 49.8% match an HPRC allele, 7.9 / 19.2 / 45.4% are patient germline by assembly, and 0.4 / 0.2% tumor-only (3.8 / 1.7% without the priority; COLO829T, 3,000 sampled: 2.4% absent from the normal). Controls: 145–149 of 150 placed HG008T true positives are tumor-only, 135–146 of 150 false positives germline; but the HPRC-allele test also tags 72 of 97 assembly-verified off-reference somatic SNVs (HG008T ONT, v2.1, autosomes).

   Complete-output F1 (chr1, no BED) is 0.533 / 0.355 / 0.369 with the rescue and at least 0.596 / 0.524 / 0.721 without it (GRCh38-only graph 0.359 / 0.365 / 0.418) [待补: rescue default (R3.4)]. The rescue recovers 46 of 89 such SNVs at tensor level (HG008T ONT, v1.1; 1 without); on v2.1 chr1 it calls 5 of 9, all HPRC-tagged. In these three runs the extra ONT calls add ambiguity (mostly population and germline alleles), not resolution.

**证据（内部）。**

- 第 1 条（node 层面，d9，2026-10-02）：`analysis/branch_node_alignments_20261002/README.md` s1 表：branch 上的 visit 占比 4.9 %（COLO829T_Illumina）– 8.7 %（HG008T_ONT），branch visit perfect 94.3 %（COLO829T_ONT）– 98.1 %（COLO829T_Illumina）；长度分档表：HG008T_PacBio 1 bp branch 98.5 % vs GRCh38 97.9 %；HG008T_Illumina 6–32 bp branch 88.0 % vs GRCh38 94.5 %，33–256 bp 79.5 % vs 89.7 %（COLO829T_Illumina 同样：89.5 % vs 95.3 %，81.8 % vs 90.9 %）。所以 detailed doc 里 "branch and GRCh38 nodes are about equally often perfect" 只对长读长成立，Illumina 的长节点不成立，Response 已改。
- 第 2 条（MAPQ，E11a 完整扫描，`complete: true`，job 387166_0–2 已 COMPLETED）：`/scratch/jshen/data/reviewer_experiments_20261013/E11a_alignment_stats/scans/HG008T_PacBio.{linear,d9,d46}.summary.json`：records 均为 14,014,144；with_path 13,624,706 / 13,644,710 / 14,012,552；`mapq_hist_with_path[60]` = 13,064,032 / 13,081,945 / 11,675,897；MAPQ > 5 之和 13,371,096 / 13,350,742 / 13,340,113。95.88 % / 95.88 % / 83.32 % 和 unmapped 2.78 % / 2.64 % / 0.01 % 是用这些计数算的（四舍五入到两位）。交叉核对：d9 的 MAPQ > 5 数 13,350,742 与 `branch_node_alignments_20261002` 的 "MAPQ>5 alignments" 完全一致，d9 的 MAPQ 60 数与 `HG008_PacBio_gam_stats.txt`（`max 60 (13081945 reads)`）一致。
- 第 2 条（vg 版本混杂）：`experiments/E3_reference_arms/README.md` 的 caveat 4：linear vg 1.66、d9 "unverified"、d46 vg 1.77。v1.1 vs v2.1 SNV PoN 后 F1 差 −0.021 到 +0.031：`experiments/E3_reference_arms/tables.md` §1 "d46 - d9" 列。
- 第 3 条：`analysis/tensor_recall_20260930/tables.md` §3 SNV 表 S1 行：HG008T PacBio 48、ONT 1、Illumina 55；COLO829T fiberseq 227、ONT 45、Illumina 237；README 第 144 行：S1 在 SMaHT Extreme（fiberseq 202/227、Illumina 208/237），COLO829T ONT S1 的 Extreme / Difficult / Easy = 38 / 3 / 4。Extreme 区 no candidate：fiberseq 367 (8.8 %)、ONT 194 (4.7 %)、Illumina 328 (7.9 %)，Easy 区 0.4 % / 0.9 % / 0.7 %（"COLO829T: status by SMaHT region" 表）。注意：这些是 2026-09-30 的 d9 tensor（10-03 ch6 重建之前），candidate 层面数字不受 ch6 编码影响，但 build 不同。HG008T 的同样结论另见 `analysis/hg008_somatic_miss_20260926/comparison_tables.md` §4（S1 48 / 1 / 55，全基因组常染色体）。
- 第 5 条计数（RIS）：`experiments/E3_reference_arms/tables.md` §2（7fc7396 RIS runs，chr1，无 BED）：HG008T SNV unplaced PASS Illumina 116 (d9) / 754 (d46)，PacBio 211 / 244，ONT 1,409 / 1,261；COLO829T ONT 1,675 / 4,566。
- 第 5 条的分类都在 re-prediction 上做（`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/`）：`vcf_chr1/*.SNV.linear.stats.json` 的 `unplaced_pass` = 240 / 1,258 / 4,569（RIS 的 E3 §2 是 244 / 1,261 / 4,566，re-prediction 的小差别）。Response 里两组数分开写：第一句是 RIS 计数，冒号后面的分类注明 "in re-predictions ... 240 / 1,258 / 4,569"。re-prediction 不都在 GPU 上（E9a README Data 表）：HG008T PacBio 是 job 387135（CPU，fp32），HG008T ONT 387379_1、COLO829T ONT 387379_2（tequila GPU，bf16）。状态行原来写 "都在 tequila GPU 上"，已改。
- 第 5 条投影（E9a，最终，job 387435_0–2 COMPLETED 15:39 / 15:45 / 15:50）：`experiments/E9a_unplaced_projection/results/tables.md`、`projection_by_run.tsv`、`metrics_by_run.tsv`、`summary.json`（16:00）。顺序 PacBio / ONT / COLO829T ONT：
  - class i（"written as a GRCh38 variant"，`project_unplaced.py` docstring；= detail exact 2 / 19 / 94 + branch_base_differs 17 / 24 / 34，后者是 call 落在 germline 差异的 branch 碱基上）19 (7.9 %) / 43 (3.4 %) / 128 (2.8 %)；class ii（只有 interval）221 (92.1 %) / 984 (78.2 %) / 1,334 (29.2 %)；class iii（没有 anchor）0 (0.0 %) / 231 (18.4 %) / 3,107 (68.0 %)，COLO829T 的 iii 里 no left anchor 2,836。
  - class ii 细分（`summary.json` 的 `details_pass`）：`alt_is_grch38_base` 89 / 517 / 608（HG008T ONT 的 517 / 1,258 是 GRCh38 碱基本身落在 branch node 上；为了字数没进 Response，R3.4 第 4 点可以用），`in_insertion` 108 / 383 / 589，`near_branch_indel` 24 / 84 / 116，`too_long` 0 / 0 / 21。注意：`tables.md` 第二张表的 "ii INDEL context" 列三行都是 0，和 `summary.json` 的 `ii:near_branch_indel` 对不上（四个 ii 子类加起来才等于 class ii 总数），E9a README Caveat 3 确认是 `collect.py` 的 key bug（读 `ii:indel_context_gap`，`project_unplaced.py` 写的是 `ii:near_branch_indel`），class ii 总数和所有指标不受影响；以 `summary.json` 为准。
  - "adding no true positive after filtering"：`metrics_by_run.tsv` 的 `bed_pon` 行，placed TP 582 / 599 / 1,690，`projected_rtg`（加上 class i，用 rtg 打分）TP 不变，FP 650 → 651、788 → 793、457 → 496。PoN 前不是 0：`bed_raw` 行 TP 607 → 608（PacBio）、627 → 628（ONT）、1,730 → 1,730；E9a README §3：这个 PoN 前的 class-i TP 在两个 HG008T run 里是同一个 truth SNV chr1:183,643,262 T>A（node 7923590，`branch_base_differs`），PoN 会标掉它。所以 Response 加了 "after filtering"。
  - "matching no missed truth SNV"：`class_ii_tp` / `class_ii_calls` 全 chr1（`nobed_*` 行）是 0 / 221、0 / 984、0 / 1,334（no BED，与本条的 F1 口径一致）；BED 内是 0 / 208、0 / 940、0 / 1,334。class ii TP 的定义（`evaluate_complete.py` docstring）：rtg 判为 FN 的 truth 落在 interval ±1 bp（`SLACK = 1`）内才算；`alt_is_grch38_base`（89 / 517 / 608）一律算 FP、不参与匹配。interval 里已经被 placed call 叫到的 truth 不算，所以不能写成 "contain no truth variant"。原来的 "the other 221 / 984 / 1,334" 对两个 ONT run 也不对：class i 以外的还有 class iii（0 / 231 / 3,107），现在单独写出来。
  - 投影后的 complete-output F1 约等于 worst case（BED 内，PoN 后）：placed 0.603 / 0.575 / 0.736；unplaced 全算 FP 0.537 / 0.358 / 0.369；投影（i 用 rtg、ii 按 interval、iii 算 FP）0.544 / 0.367 / 0.373。这组数没进 Response（R3.4 的范围），这里只用来说明投影救不回 F1。
  - CHM13 列全是 0：`chm13_pass` = `not_built`（用户放弃了 CHM13）。
- 第 5 条 graph PoN 和组装分类（E9c / E9d，最终）：`experiments/E9c_offref_filter_validation/tables.md`（16:07，`rcpu__*` 行）；job 387431_0–2、387432 / 387433、387434 都已 COMPLETED（16:07:50）。顺序 PacBio / ONT / COLO829T ONT：
  - §1（互斥，顺序 population > germline > tumor-only > unverifiable > unsupported）：population 214 (89.2 %) / 932 (74.1 %) / 1,495 (49.8 %)；germline 19 (7.9 %) / 242 (19.2 %) / 1,362 (45.4 %)；tumor-only 1 (0.4 %) / 3 (0.2 %) / 0 (0.0 %)，COLO829T not in normal 72 (2.4 %)；unverifiable 5 (2.1 %) / 70 (5.6 %) / 71 (2.4 %)；unsupported 1 (0.4 %) / 11 (0.9 %) / 0 (0.0 %)。COLO829T ONT 是 4,569 个里随机 3,000 个（`processed` 列）。COLO829T 没有 tumor assembly，只有 germline / not in normal / unverifiable（`README_methods.part` Method 第 4 条）。
  - 只看组装、不看优先顺序（E9c README §1 第二张表；`tables.md` §3 的 `unplaced_pass` 行）：germline 158 (65.8 %) / 724 (57.6 %) / 2,493 (83.1 %)；tumor-only 9 (3.8 %) / 22 (1.7 %)，其中被 graph PoN 标的 8 / 19；COLO829T not in normal 197 / 3,000，其中被标的 125。Response 的 "3.8 / 1.7% without the priority" 来自这里。COLO829T ONT 的 multi-copy context 有 1,594 / 3,000，germline 类里可能有 paralog 拷贝（E9c README Caveats）。
  - §3 对照（各 150 个随机）：ctrl_tp（placed、rtg 判 TP）组装分类 tumor_only 149 / 145，COLO829T not_in_normal 150；graph PoN 标了 0 / 1 / 0。ctrl_fp（PoN 后的 placed FP）germline 146 / 135 / 110（COLO829T 其余 40 是 not_in_normal）。Response 的 "145–149" / "135–146" 只指 HG008T 两个 run；COLO829T 没有 tumor assembly，对照只能说 "absent from the normal"，FP 里 germline 也只有 110 / 150。
  - §5 gold（chr1，d46）：有 tensor 的 gold key 11 / 9 / 2，unplaced PASS 3 / 5 / 1，graph PoN 标了 9 / 7 / 2，"PASS and tagged" 3 / 5 / 1，也就是 rescue 叫出的 gold 全部会被 graph PoN 去掉。Response 写 "on v2.1 chr1 it calls 5 of 9, all HPRC-tagged"。chr1 上的 gold 很少（每个 run 2–11 个，E9c Caveats），只能当旁证。
  - graph PoN 对 off-reference somatic 不特异（E9c `tables.md` §4、README §4，main rule，全常染色体）：d46 gold PacBio 90 / 119 (75.6 %)、ONT 72 / 97 (74.2 %) 被标；chr2-22 上整个 flank × AC grid 的最低值也有 45.7 % / 42.5 %（adaptive flank，AC ≥ 9）。README Interpretation 3："not specific for off-reference somatic variants"。所以 "population" 不能读成 "一定不是 somatic"，Response 写 "match an HPRC allele"，并给出 72 of 97。组装检查在 gold（按构造是 somatic）上也把 d46 PacBio 11 / 119、ONT 9 / 97 判成 germline，所以在 off-reference 位点上组装分类不如 placed 对照那么干净（Response 只说对照在 placed call 上成立）。
  - §2 校正后 F1（no BED，GRCh38 PoN 后）：placed 0.598 / 0.5633 / 0.7364；unplaced 全算 FP 0.5327 / 0.3548 / 0.369；去掉 graph PoN 标的 0.5901 / 0.4889 / 0.4911；再把 tumor-only 算 TP 0.5908 / 0.4906 / 0.4911。BED 版用 bubble interval 判断是否在 BED 内（计入的 unplaced 226 / 1,182 / 4,569），和 E9a 的 BED 口径不同，所以 BED 的 "all FP" 是 0.5401 / 0.3668 / 0.369，E9a 是 0.537 / 0.358 / 0.369。Response 用 no BED，与 R3.4 一致。
  - §7 重复背景（homopolymer ≥ 7 / STR / other）：PacBio 97 / 60 / 83；HG008T ONT 387 / 219 / 652，其中 population 322 / 185 / 425；COLO829T ONT 303 / 289 / 2,408（other 里 germline 1,276、population 1,075）。没进 Response（字数）。
  - §8 跨平台（d46，HG008T，同一 node:offset:allele）：共享 190 个 key，占 ONT 1,258 的 15.1 %、PacBio 240 的 79.2 %；ONT 这 190 个里 population 167、germline 16、tumor_only 1、unverifiable 6。为了字数没进 Response（R3.3 第 5 点用）。
  - graph PoN 用的是 HPRC v1.1 的 89 条 haplotype 加 in-graph walk（v2.1 的 full-graph VCF 不在本地，`README_methods.part` Data 一节）。
  - legacy COLO829T ONT d9（387265_9）也已 collect：2,515 个里 population 75.5 %、germline 17.8 %、not in normal 2.4 %、unverifiable 4.4 %。开发模型，没进 Response。
- 第 5 条 F1：with rescue 0.533 / 0.355 / 0.369 = E9c §2 no BED "F1 all FP" 0.5327 / 0.3548 / 0.369（re-prediction；RIS 的 E3 §2 是 0.533 / 0.354 / 0.369）；without rescue ≥ 0.596 / 0.524 / 0.721 和 GRCh38-only 0.359 / 0.365 / 0.418 来自 R3.4 证据（"不用 rescue：0.5961 / 0.5237 / 0.7210"，off-ref tensor 没算 AF floor，所以是下界）和 E3 §2 linear 行。原稿第 5 条的 "F1 0.342 / 0.354 vs 0.365" 和 "COLO829T ONT v2.1 0.369 vs 0.418" 删了，第 2 点在 R3.4 里有完整表。2026-10-13 晚复核：从 re-prediction 的 `test_chr1/*.predictions.ndjson.gz` 重数 `off_reference` 且 `p_somatic` ≥ 阈值（0.29122 / 0.23086 / 0.56769）的 tensor = 6 / 162 / 98（off-ref tensor 12,632 / 26,665 / 31,721，excluded 0），用 2TP / (2TP + FP + FN + U) 和 `metrics_by_run.tsv` `nobed_pon` placed 的 586 / 672 / 116、603 / 836 / 99、1,690 / 457 / 753 算出 0.5961 / 0.5237 / 0.7210，和 R3.4 一致。re-prediction 的 `metrics.json` `offref_rescue`：calls 240 / 1,258 / 4,569（= `unplaced_pass`），calls_before 3 / 48 / 262。
- 第 5 条 rescue 和 gold（E1.1）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §5（单元格 = p_somatic >= t / argmax somatic / n (chr1: >= t / n)）：`HG008_ONT`（run `HG008_ONT_SNV_7fc7396`，d9，t = 0.449）base "1 / 1 / 89 (chr1 0 / 8)"，`site10_high`（= predict 的 rescue 变换）"45 / 46 / 89 (chr1 3 / 8)"；E1.1 README："argmax 32 / 112 and 46 / 89, i.e. the pipeline's own off-reference calls"。keys 来自 `pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv`，所有常染色体，"None of them is a training label"。pipeline 对 off-ref SNV tensor 用 rescue 后的 argmax（`machine_learning/pansoma_net_v2/predict.py` docstring "Off-reference rescue"；`graph_vcf.py`："called by offref_call ... instead of the threshold"）。这是 tensor 层面的 argmax（没算 AF floor），Response 现在写 "at tensor level"，和 R3.4 一致。Response 原来的 "only 8 of them on chr1 (4 rescued)" 和 "162 off-reference tensors instead of 1,258"（`HG008T_ONT_d46.tensors.SNV.predictions.ndjson.gz` 自己数的）为了字数删了，后者在 R3.4 Response 里。
- Rescue = unplaced：`/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/HG008_ONT_SNV_7fc7396/test_chr1/HG008T_ONT.tensors.SNV.metrics.json` 的 `offref_rescue`：calls 1,409（= E3 §2 的 unplaced PASS），calls_before 30；HG008T ONT d46 calls 1,261、calls_before 48；COLO829T ONT d46 calls 4,566、calls_before 262。COLO829T ONT d9 的 1,675（`COLO829_ONT_SNV_h90`）metrics 不在本地，没核对。
- 数字不一致（用较新的）：简版 `response_to_reviewers.md` R1.3 写的是 105 (PacBio) / 595 (Illumina) / 855 (ONT) / 2,515 (COLO829T ONT)，来自 09-30 / 10-01 的 legacy `*_SNV_pow075_S_lr1e4_e20` d9 runs。7fc7396 的数字（211 / 116 / 1,409）见 `analysis/benchmark_7fc7396_20261012/hg008t_snv_indel_7fc7396_20261012/README.md`。E9c `README_methods.part` 写 "SNV runs reproduce the PacBio / ONT ones exactly"，和 legacy PacBio 105 vs 7fc7396 211 矛盾；现在 Response 只用 `rcpu__*` 行，不再引用 legacy 的 78.8 % / 83.8 %。
- Job 状态（2026-10-13 晚 sacct）：E9a 387435_0–2、E9c 387431_0–2、E9d 387432 / 387433、collect 387434 都 COMPLETED；legacy E9c 387265_2–7 仍 PENDING（held；其中 _3 / _4 是 HG008 ONT），_8 / _10 FAILED（11:49–11:50），_9（legacy COLO829T ONT）COMPLETED；E11b alignment 387342–387347 和 collect 387349 PENDING；E11a 387166_[5-6]、387287_[2,6-10]、387288_[12-…] PENDING。

**还缺什么。**

- [CPU 本地] E10（三个 d46 re-prediction 先做）：rtg 输出在本地（`reviewer_cpu_20261013/<run>/vcf_chr1/eval_{bed,nobed}`）；按 stratum 数 tp / fp / fn（PoN 前后；COLO829T 用本地 SMaHT v2 BED，HG008T 用 SMaHT v2 + `GRCh38_segdups.bed` 或 GIAB），加 block bootstrap CI（E16）；顺带把 E9a 的 class i / ii interval 分到 stratum 上，给 unplaced calls 一个区域背景。脚本还没写。
- [数据同步] E10 的 linear / d9 臂和其余平台：本地 `ris_7fc7396_20261012/<run>/` 只有 `test_chr1/*.metrics.json`，rtg 输出和 predictions 在 RIS 上；要同步，或者用本地的 `best.pth` 在 tequila GPU 上 re-predict。
- [用户决定] HG008T 的分层用不用 GIAB stratification：本地没有，v3.5 GRCh38 tarball 1.3 GB，另有 v3.6（`experiments/E8b_literature_facts/README.md` G7）；不下载就把稿子里的 "GIAB v3.5 stratification" 改成 SMaHT v2 + segdup。
- [用户决定] SNV off-ref rescue 是否保留为默认（与 R3.3 / R3.4 同一个决定）。现在三个 d46 run 的两边都有数：保留时 complete-output F1 0.533 / 0.355 / 0.369，unplaced 里 tumor-only 只有 1 / 3 个（互斥顺序；只看组装是 9 / 22 个，其中 8 / 19 个被 graph PoN 标；COLO829T 抽样 3,000 个里 not in normal 72，只看组装 197）；关掉时至少 0.596 / 0.524 / 0.721，代价是 d9 gold 46/89 → 1/89（ONT，tensor 层面），d46 chr1 上 PASS 的 gold 3 / 5 / 1 个（全被 graph PoN 标）也会没有。复核：这 3 / 5 / 1 个 key 在 predictions 里的 `p_somatic` 都低于阈值；反过来，ONT 有 1 个 gold key（10160480:0:SNP:A>T，p 0.24051 ≥ 0.23086）不用 rescue 时过阈值，用 rescue 时 `offref_call` 反而是 false（AF floor 没算）。定了以后改第 5 条最后一段的 [待补]。
- [用户决定] HG008T ONT 的 v1.1（d9，1,409 个）要不要也分类：本地没有 7fc7396 d9 re-prediction，只有 legacy HG008 ONT（387265_3 / _4，held）；不跑就只报 v2.1 的分类，v1.1 只给计数（现在 Response 就是这样写的）。
- [CPU 本地]（可选）E9c / E9d 的 COLO829T ONT 只处理了随机 3,000 / 4,569 个 unplaced PASS（seed 0；E9c §2 的 precision 按样本比例外推到 4,569，`unplaced left` 2292.1 是估计）。要全量数就补跑剩下的 1,569 个；现在 Response 已注明 "3,000 sampled"。
- [数据同步] E9a / E9c 只覆盖三个 d46 SNV run（chr1）。HG008T / COLO829T Illumina、COLO829T fiberseq 和 d9 的 graph SNV run 要 predictions（E9a README "Status details"），之后才能说 "ONT 比别的平台多出来的是 population / germline"，而不只是这三个 run 上如此。
- [CPU 本地] E11a：已完成的 window job 387287_0、_1、_3、_4、_5（HG008T Illumina linear / d9 / GRCh38 BAM，HG008T PacBio / ONT GRCh38 BAM）可以先用 `compare.py` 出按 SMaHT / segdup 分层的表；387166_[5-6]、387287_[2,6-10]、387288_[12-…] 仍 held。
- [CPU 本地] E11b：三个平台的模拟都齐了；alignment job 387342–387347 和 collect 387349 仍 held，需要用户放出；CHM13 臂按用户决定已放弃。
- [文字] 简版 R1.3 的 105 / 595 / 855 / 2,515 换成 7fc7396 数字，删掉 "resolution" 的说法，改为 "added ambiguity, mostly population and germline alleles"，并和 Response 一样写明：只有三个 run；类别互斥、有优先顺序；graph PoN 是 "match an HPRC allele"，对 somatic 不特异；"equally often perfect" 限定为长读长；rescue 句子和 R3.3 / R3.4 / R1.1 保持一致（同一组 89 / 112 个 keys；rescue 叫回 46 / 32，不 rescue 1 / 0）。
- [文字] Response 现在约 435 词（不算 "/"），比 400 多：第 5 条为了写明 "互斥优先顺序"、class iii、"after filtering" 和 graph PoN 不特异（72 of 97）多了约 25 词。如需再压，第 1、3 条可以移一部分到 Supplementary。

### R1.4 Candidate generation removes somatic variants already encoded in the graph

**状态：部分解决**: 机制已接受。下面几项都已有数字：read 层面的吸收分类，和现行 PoN 的交叉表，和 GRCh38-only 臂在同一批 reads 上的 tensor 层面对比，ch6 在 classifier 里对 off-reference somatic 的压制（E1.1 §5）。新的是 E5b 最终结果（graph spelling 已跑完，chr1–22）。相对 linear 丢掉、又没被 PoN 标的 truth INDEL 里，"被图吸收"（没有 candidate，ALT 被图 spell 出来）只占小部分：HG008T v1.1 是 481 / 581 / 299 个里的 35 / 34 / 20 个。graph spelling 是严格的精确匹配，这个数是下界；按 read 层面的吸收集合算是 49 / 74 / 40 个（E5b 表 2b），也还是小部分。多数是 "没有 candidate、图也不 spell ALT" 和 "candidate 被 AF 过滤"。吸收的 truth 在每套数据、两个图上都有 95.5–100 % 被 PoN 标掉，"not spelled" 的 INDEL 只有 52.8–78.7 %。还缺几样数字：E5b 的 `linear_recovery.md`（丢掉 / 吸收的 truth INDEL 在三个臂上叫没叫出来）要 RIS 的 36 个 rtg 输出，或者在 tequila GPU 上重新 predict 18 个 INDEL run；在那之前，"linear 保住了这些 INDEL" 只是 tensor 层面的说法；比较器在吸收位点上的表现（E5c）要用户给出已有 ClairS-TO / DeepSomatic 结果的路径（不重跑）；INDEL 的 held-out 阈值规则要从 RIS 同步。branch candidate（E6 / E25）做不做要用户定。

**Comment（要点）。**

- Candidates come only from mismatches / indels against the traversed node. A somatic allele that is already an alternative path, with tumor reads aligning perfectly to it, gives no candidate, tensor or call.
- Called a "fatal methodological issue".
- Supplementary Fig. 3 wrongly describes these false negatives as GRCh38 reference-bias artefacts.
- An allele on an alternative path is not necessarily germline: the same allele can be germline in one person and somatic in another.

**Response.**

The reviewer is correct: a graph-path allele is not necessarily germline, and reads aligning to it without edit give no candidate. We withdrew the reference-bias reading of Supplementary Fig. 3 and quantified the loss.

*Read level (HPRC v1.1, patient assemblies).* ALT reads of 2,478 HG008T truth alleles (201 SNVs, 2,277 INDELs; 26.8% of truth INDELs) and 497 COLO829T alleles (67 SNVs, 430 INDELs; 22.5%) align without edit to an existing path. Only 110 (4.4%) and 34 (6.8%) are the patient's germline allele, 2,085 (84.1%) and 430 (86.5%) are carried only by other HPRC individuals, and 1,853 HG008T INDELs change a homopolymer ≥ 7 bp in the patient's sequence.

*Population filtering.* Tumor-only pipelines remove such alleles at filtering. Our filter (R2.6) tags 98.5–99.0% of these SNVs and 94.7–97.0% of these INDELs, against 3.6–10.2% and 52.0–61.8% of all truth SNVs and INDELs. Absorbed and untagged: HG008T 2 SNVs and 121 INDELs (1.3% of truth INDELs), COLO829T 1 SNV and 13 INDELs (0.6%); on HPRC v2.1, by graph spelling per platform, 0.2–0.7%.

*Against a GRCh38-only graph* (same reads and code, chr1–22; R3.1). Net of gains, HPRC graphs keep 0.2–0.6% fewer truth SNVs and 2.0–8.2% fewer truth INDELs with a tensor after filtering. Absorption is a small part: of the untagged INDELs lost on HPRC v1.1 (HG008T Illumina / PacBio / ONT: 481 / 581 / 299), 35 / 34 / 20 are absorbed by a strict spelling test (read-level: 49 / 74 / 40); most get no candidate without such a path (274 / 339 / 145) or fall below the candidate allele-frequency cut (153 / 125 / 84). The filter tags 95.5–100% of absorbed losses on every dataset and graph, 52.8–78.7% of the other no-candidate INDEL losses. The GRCh38-only graph gives these INDELs a tensor but also makes the patient's germline INDELs candidates: its best chr1 INDEL F1 after filtering is lower on all six datasets (HG008T 0.199 / 0.173 / 0.082 against 0.264 / 0.260 / 0.140 on v1.1), but higher at the validation-chosen threshold on the four long-read ones (HG008T PacBio 0.125 against 0.089) [待补: held-out threshold rule; E5b linear calls].

On HCC1395, whose SEQC2 truth INDELs are rarely population alleles (7.3% tagged) or spelled by HPRC v2.1 (2.9%; HG008T 21.3%), 76–85% reach a candidate on v2.1, against 43–50% for HG008T.

Linear tumor-only callers also drop population-like alleles: an existing DeepSomatic tumor-only run (COLO829T Fiber-seq, chr1) labels 72 of 131 truth INDELs germline, 55 of them filter-tagged [待补: ClairS-TO / DeepSomatic recall at the absorbed loci (E5c)].

The classifier applies the population prior a second time through the haplotype-count channel: of the assembly-verified HG008T somatic SNVs on off-reference v1.1 nodes, 0 of 112 (PacBio) and 1 of 89 (ONT) pass the threshold, against 61 and 49 with the channel at every haplotype (R3.3). The Discussion states both as limitations of population-based tumor-only calling. [可选, 用户决定: branch-candidate mode (E6 / E25)]

**证据（内部）。**

- read 层面吸收（d9，三平台 union）：
  - `analysis/graph_absorbed_somatic_20261001/README.md`："201 SNVs + 2,277 INDELs = 2,478 truth alleles"；"2,085 of the 2,478 (84.1%)"；"median HPRC frequency 0.17"；"110 (4.4%) are an HG008-N germline allele"。
  - `analysis/graph_absorbed_somatic_colo829t_20261005/README.md`："67 SNVs + 430 INDELs = 497 truth alleles (22.5% of all 1,912 chr1-22 truth"；"430 (86.5%)"；"Germline filtering: 34 (6.8%)"；"Overall INDEL absorption: 22.5% vs 26.8%"（HG008T 2,277 / 8,496 = 26.8 %）。
  - "1,853 … homopolymer ≥ 7 bp"：`analysis/graph_absorbed_somatic_20261001/README.md` "Repeat context" 表的 union 列（patient frame，即 HG008-N v6.2 event haplotype 上的序列）；GRCh38 frame 下是 1,688。
- 和现行 PoN（7fc7396 规则）的交叉表：`tmp/graph_absorbed_pon_20261012/crosstab.md`、`crosstab_colo.md`、`README.md`。
  - HG008T：absorbed SNV 199 / 201（99.0 %），INDEL 2,156 / 2,277（94.7 %）；all truth SNV 982 / 9,609（10.2 %），INDEL 5,783 / 9,362（61.8 %）；absorbed INDEL 里 patient_class HP>=7 = 1,853。
  - COLO829T：absorbed SNV 66 / 67（98.5 %），INDEL 417 / 430（97.0 %）；all truth SNV 3.6 %，INDEL 52.0 %。
  - 图特有损失："HG008T 2 SNV + 121 INDEL (1.3 % of 9,362 INDELs), COLO829T 1 SNV + 13 INDEL (0.6 %)"（分母是全部染色体、ALT 拆开后的 truth INDEL）。E5b 表 3 的分母不拆 ALT（HG008T 8,496，tagged 61.7 %），所以和这里的 61.8 % 差 0.1，Response 用交叉表的数。
- E5b 最终（graph spelling 已完成）：`experiments/E5b_absorbed_d46/tables.md`、`summary.json`、`t2_compare_d9.txt`。22:53 第一次生成；23:09 修了 `s5_tables.py` `pon_tags()` 的换行 bug（README Caveat 6）后重新生成，README 说 tables.md / summary.json 的数不受影响。2026-10-13 晚又对表 1、2、2b、3、5、7、8 逐项核对了一遍，和下面引用的一致；job 387428（s3 extract）/ 387429（s4 spell d46）/ 387430（s4 spell d9）COMPLETED 22:47 / 22:50 / 22:52。`s5_prelim.stderr.txt`（11:49，"no spell_d9.tsv.gz yet"）是早上 preliminary 那次的 log，已过时。类别定义见 `s5_tables.py` docstring：absorbed = 该臂 no_candidate 且 ALT 被该臂的图 path spell 出来；nocand_other = no_candidate 但不 spell；filtered = "a candidate below the AF / min-variant cut"；other_allele = tensor_non_representative_allele。
  - 表 1（chr1–22，全部 truth 记录，不限 BED；"tensor and not tagged" 的净差 = 丢的减去多出的，占全部 truth）：SNV −0.2 % 到 −0.6 %；HG008T INDEL −4.8 % / −4.4 %（Illumina，d9 / d46）、−6.0 % / −5.8 %（PacBio）、−2.0 % / −2.0 %（ONT）；COLO829T INDEL −6.8 % / −6.2 %、−8.2 % / −8.1 %、−4.6 % / −4.9 %。和 preliminary 版相同。
  - 表 2（每个平台，相对 linear 丢掉的 truth，按原因拆；"after PoN" = 没被 PoN 标的）：HG008T INDEL d9 lost after PoN 481 / 581 / 299（Illumina / PacBio / ONT），其中 absorbed 35 / 34 / 20，no candidate not spelled 274 / 339 / 145，filtered 153 / 125 / 84，other allele 19 / 83 / 50（四类加起来正好 481 / 581 / 299）；d46 是 464 / 571 / 290，absorbed 39 / 32 / 20，not spelled 280 / 361 / 146，filtered 123 / 95 / 60，other 22 / 83 / 64。COLO829T INDEL d9 133 / 159 / 114，absorbed 4 / 3 / 1；d46 124 / 157 / 116，absorbed 6 / 5 / 2。gained after PoN（HG008T INDEL d9）74 / 68 / 132，481 − 74 = 407 = 表 1 的 −407。
  - 表 2 的 PoN 比例：absorbed 一列 24 行全部在 95.5 %（HG008T ONT INDEL d9，441 个里 421）到 100.0 % 之间；"no candidate, ALT not spelled" 的 INDEL 行 HG008T 71.8–78.7 %、COLO829T 52.8–57.5 %（Response 写 52.8–78.7 %）；SNV 行 HG008T 73.2–87.4 %、COLO829T 只有 7.8–17.5 %（COLO829T 的 SNV 损失主要不是 population allele，Response 没展开）。
  - 表 2b（按 read 层面吸收集合拆，口径不同：union of platforms，09-30 tensor）：HG008T INDEL d9 lost after PoN 里在 read 层面集合中 49 / 74 / 40，其余 432 / 507 / 259。比表 2 的 graph-spelling 口径（35 / 34 / 20）大，因为 read 层面集合是三平台 union、不要求该平台 no_candidate。Response 用表 2，括号里给出表 2b 的 "read-level: 49 / 74 / 40"。
  - spelling 是严格测试，absorbed 是下界（E5b README Caveat 1–2）：要求 window anchor 之间精确的 GRCh38 + truth haplotype，漏掉 "with_germline" 和 "closest"（t1 的 300 个 read 层面吸收位点里 44 个算成了 "not spelled"）；vg find 的 250 bp subgraph 不完整，complete d9 edge table 能 spell 的 256 个里有 22 个（8.6 %）在 subgraph 上 spell 不出来（`t2_compare_d9.txt`：278 / 300 一致），所以 d9 的 absorbed 是下界，d46 没有同样的检查。部分 "not spelled" 其实是吸收。两种口径下吸收都只占 PoN 后损失的小部分（35 / 34 / 20 或 49 / 74 / 40，对 481 / 581 / 299），所以 Response 写 "absorbed by a strict spelling test (read-level: 49 / 74 / 40)"，结论不变。另外三个臂是同一批 reads 的不同 GAM（Caveat 4），"not spelled" 混了比对和 repeat 表示的差别。
  - 表 3（三平台 union）：HG008T INDEL absorbed on ≥ 1 platform d9 1,818（96.8 % tagged，59 after PoN）、d46 1,467（96.5 %，51）；COLO829T d9 382（98.4 %，6）、d46 347（98.0 %，7）。
  - 表 5（与 reads 无关，ALT 是否被图 path spell）：HG008T INDEL d9 2,319 (27.3 %)、d46 1,809 (21.3 %)，d46 spell 的 96.2 % tagged；COLO829T INDEL d46 398 (20.8 %)，97.7 %；HCC1395T INDEL d46 47 (2.9 %)，100.0 %。和 read 层面分析对照：HG008T read 层面吸收的 2,277 INDEL 里 1,839 个被这次的 d9 spelling 确认，COLO829T 430 里 382 个。`t2_compare_d9.txt`：sqlite 和 vg find 两种 spelling 方法在 300 个位点上 278 个一致。
  - 表 7 第二张（d46，no candidate 且 ALT 被 d46 spell，不和 linear 比）：INDEL absorbed after PoN 占 truth INDEL 的 HG008T 0.7 / 0.5 / 0.4 %、COLO829T 0.4 / 0.3 / 0.2 %、HCC1395T 0.0 %（Response 的 "on HPRC v2.1, by graph spelling, 0.2–0.7%"；HCC1395 没写进 Response）。
  - 表 8（chr1 BED，vcfeval 的范围）：相对 linear 丢掉且没被 PoN 标的 truth INDEL（括号里是其中 absorbed）：HG008T d9 27 (0) / 30 (1) / 18 (1)，d46 25 (2) / 31 (2) / 20 (2)；COLO829T d9 11 (1) / 13 (1) / 10 (0)，d46 11 (1) / 13 (1) / 10 (1)。即 chr1 测试范围里吸收造成的、PoN 也挡不住的损失每套 0–2 个。
  - `linear_recovery.md`（23:09）：所有格子都是 "pending"，"pending rtg outputs: 36"（三个臂 × 六套的 `pon_curve` / `pon_calls`，在 RIS 上）。只列出了 chr1 BED 内丢掉 / 吸收的 truth 数（例 HG008T Illumina lost on d9 250 个、其中 not PoN-tagged 27；absorbed on d9 114 个、0），叫没叫出来还没有。修 bug 后 COLO829T 的 absorbed n 是 26 / 23、21 / 21、16 / 18（Illumina / fiberseq / ONT，d9 / d46；22:53 版是 1 / 1、1 / 1、0 / 1，漏掉了所有被 PoN 标的），not PoN-tagged 一列（1 / 1、1 / 1、0 / 1）不变。
- d46 vs d9：E5b 表 1，truth INDEL 有 tensor 的数目 d46 在六套上都多（例 HG008T Illumina 3,555 → 3,821，COLO829T fiberseq 1,210 → 1,260）；表 4，去掉 PoN 后净差 HG008T INDEL +30 / +19 / +1，COLO829T +12 / +2 / −6。Response 为了字数删了 "v2.1 gives a candidate to more truth INDELs than v1.1" 这句。
- HCC1395：E5b 表 1，d46 truth INDEL 有 tensor 0.828 / 0.846 / 0.762（Illumina / PacBio / ONT），HG008T d46 0.450 / 0.497 / 0.429；PoN 标的 HCC1395T INDEL 119 / 1,625（7.3 %，表 7 第一张）。旁证：`experiments/E7_graph_composition/d9_content.tsv`（bed scope，`in_d9_vcf`）somatic INDEL 是 d9 allele 的比例 HG008 0.252、COLO829 0.2108、HCC1395 0.0327（E7 没有 README，算初步）。
- end-to-end：`experiments/E3_reference_arms/tables.md` §1，INDEL best F1（PoN 后，BED 内）linear 0.199 / 0.173 / 0.082 vs d9 0.264 / 0.260 / 0.140 vs d46 0.249 / 0.269 / 0.155（HG008T）；COLO829T 0.069 / 0.092 / 0.019 vs 0.195 / 0.120 / 0.071 vs 0.163 / 0.129 / 0.085。PASS F1（同表）四套长读长都是 linear 高：HG008T PacBio 0.125 vs 0.089、ONT 0.022 vs 0.015，COLO829T fiberseq 0.040 vs 0.016、ONT 0.004 vs 0.003；Illumina 是 d9 高（0.241 vs 0.133、0.115 vs 0.032）。这些都是 BED 内、unplaced calls 不计分（E3 §2 的 complete-output bound 另算）。linear README：PoN 后 ceiling 0.284 / 0.301 / 0.238 vs d9 0.257 / 0.259 / 0.238。Response 原来写 "keeps these INDELs but makes every germline INDEL a candidate"，改成 "gives these INDELs a tensor but also makes the patient's germline INDELs candidates"：E5b 的 lost 是 tensor 层面的 ceiling（README Caveat 3），linear 臂真叫没叫出来要等 `linear_recovery.md`（Response 加了 [待补: ... E5b linear calls]）；germline INDEL 在 linear 上也要过 candidate 的 AF cut，"every" 说过头了。PASS F1 在四套长读长上 linear 高，但绝对值都很低（≤ 0.125；COLO829T ONT 0.004 vs 0.003 / 0.002）。
- ch6 对 off-reference somatic 的压制：`experiments/E1.1_ch6_counterfactual/results/tables.md` §5 和 `results/summary.json`（`assembly_keys`、`results.<mode>.assembly`）。
  - 对象：`pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv` 里有 tensor 的 SNP key（HG008-N assembly 事件，somatic allele 被 giraffe 放在 d9 的非 GRCh38 节点上），全部常染色体；summary.json 里 112 / 89 个 key 的 `off_reference` 全是 true，都不是训练标签。
  - "p_somatic >= t / argmax somatic / n"：base `0 / 0 / 112`（HG008_PacBio，d9）、`1 / 1 / 89`（HG008_ONT，d9）；all_H `61 / 44 / 112`、`49 / 49 / 89`；train_mean `30 / 26 / 112`、`5 / 4 / 89`。
  - site10_high（README："`predict`'s off-reference rescue transform applied to every tensor"）：`45 / 32 / 112`、`45 / 46 / 89`。pipeline 对 off-reference SNV tensor 用 rescue 后的 argmax 叫（`machine_learning/pansoma_net_v2/predict.py` 的 "Off-reference rescue" 说明），所以 Response 写 rescue 叫出 32 / 46。
  - §5 只有 d9 的两个 run。和 R3.3 / R1.1 共用，这里只写一句；为了字数，"all autosomes, none a training label" 和 rescue 的 32 / 46 只留在 R3.3。E9c 的 d46 结果（graph PoN 标掉 chr1 gold 9/11、7/9、2/2）在 R1.3 / R3.3，不放这里。
- DeepSomatic：`analysis/tensor_recall_20260930/pon/indel_rules/README.md` §2b（v1.9.0 `PACBIO_TUMOR_ONLY`，默认 PoN，2025-08 的旧 run）：COLO829T chr1 INDEL 131 个里 GERMLINE 72（55.0 %），PASS 12（9.2 %）；"55 of those 72 are alleles our PoN rule tags"。README 的 caveat：fiberseq 不是 PacBio 模型的训练数据，绝对数要谨慎。按用户决定不重跑比较器，这里只引用已有 run，所以 Response 写 "an existing DeepSomatic tumor-only run"，不写成一般结论。
- DeepSomatic 一句原来写 "Linear tumor-only callers lose the same alleles"：这个 run 只说明 131 个 truth INDEL 里 72 个被标 GERMLINE（55 个我们的 PoN 也标），没说明它们就是图吸收的那批，所以改成 "also drop population-like alleles"；吸收位点上的比较是 E5c。
- 和详版的出入（以新来源为准）：
  - 详版 "99% of the absorbed SNVs"：实际 99.0 %（HG008T）/ 98.5 %（COLO829T），改成 98.5–99.0 %。
  - 详版 "4–10% of all truth SNVs"：实际 3.6 % / 10.2 %。
  - 详版 "71–80% of truth insertions … against 36–42%"：口径是 insertions；这里换成 E5b 的全部 INDEL（76–85 % vs 43–50 %）。
  - 详版只有 read 层面的 1.3 % / 0.6 %；E5b 的端到端对比（INDEL 2.0–8.2 %）更大，因为它包括所有原因。现在 E5b 也给了拆分：端到端损失里吸收只占小部分，主要是 "not spelled" 和 "filtered"。三个都写进 Response，没有隐藏端到端的数。

**还缺什么。**

- [数据同步] E5b `linear_recovery.md`：要 RIS 上 36 个 rtg 输出（linear / d9 / d46 × 六套 INDEL run 的 `pon_curve` / `pon_calls`）；另一条路是 E5b README Status 第 1 条写的，在 tequila GPU 上重新 predict 18 个 in-sample INDEL run（要等用户定 lab cluster 的 cap）。拿到后跑 `collect.sh`（会调 `s6_linear_recovery.py`；`pon_tags()` 的 bug 已修），回答 "丢掉 / 吸收的 truth INDEL 在 linear 臂上是否真被叫出"（chr1 BED 内每套 0–2 个吸收且未被 PoN 标的，见表 8）。
- [数据同步] E5c：用户给出已有 ClairS-TO / DeepSomatic 结果的路径（不重跑）。
- [CPU 本地] E5c：在 2,478 / 497 个吸收位点上数 PASS，分 raw、各自默认过滤、同一 PoN 三档。
- [数据同步] INDEL 的 held-out 阈值规则（validation 上选的阈值）要从 RIS 同步 predictions。之后 linear vs HPRC 的 INDEL 比较才能不用 test 上的 best F1。
- [用户决定] graph-specific loss 主文报哪个数：read 层面 1.3 % / 0.6 %、d46 graph spelling 的 0.2–0.7 %，还是相对 linear 的端到端 INDEL 2.0–8.2 %（其中吸收只占小部分）。建议都报，Response 现在就是都报。
- [用户决定] E6 / E25（branch candidate）做不做；不做就只写成 limitation。按 E5b，吸收的 truth 95.5–100 % 本来就被 PoN 标掉，branch candidate 能多叫回来的、PoN 后还在的 truth 很少（chr1 BED 内每套 0–2 个；strict spelling 口径，是下界，按 read 层面的定义会多一些），可以作为不做的理由。
- [文字] 撤回 Supplementary Fig. 3 的 reference-bias 解释；Discussion 加 limitation（alignment 时的吸收，加上 ch6 在 classifier 里压制 off-reference somatic，和 R3.3 的写法对齐）；简版 R1.4 换成上面的数字（98.5–99.0 %、3.6–10.2 %、端到端对比和 E5b 的原因拆分）。
- [文字] Response 现在约 460 词（不算 "/"），比 400 多。这次为了写明 strict spelling 的下界（read-level 49 / 74 / 40）、"gives a tensor" 和 [待补: E5b linear calls] 多了几个词。如需压到 400 左右，HCC1395 一句和 DeepSomatic 一句（约 60 词）可以移到 Supplementary。

### R1.5 HapMap mixture as a somatic benchmark, and data leakage

**状态：部分解决**: "HapMap 退出主结果、换成含真实 somatic 突变的 dilution" 这条主线已有 2026-10-12 RIS 上 COLO829T 三平台 purity joint model 的结果（fiberseq 没有 artefact，可以直接用），donor 在图里的事实也已核实（E8b），E1.1 又表明训练好的分类器确实用 ch6（path count）压 germline allele（过 PoN 以后影响小得多）；但 Illumina / ONT dilution 带 read-source BQ artefact，分数偏乐观（要在 RIS 上修），HCC1395 dilution（E14c）没开始，donor-free 图上的 HapMap 重跑还在 sort，HapMap 去留和原稿用的是哪张图要用户定。

**Comment（要点）。**

- HapMap 的 "somatic-like" positives 是少数个体（minor donors）的 germline allele；Pansoma 的前提是图会压掉 population germline allele，两者冲突，这个 benchmark 可能系统性地对 pangenome caller 不利，也不是获得性低 VAF 突变的合适检验。
- 数据泄漏：HPRC v1.1 含 6 个 donor 中的 4 个（HG00438、HG02257、HG02486、HG02622）；companion paper（Fu, Xin, Miao et al., bioRxiv 2026.01.04.697580）提过这个问题，但没有排除这些 donor。

**Response.**

We agree with both points. The HAPMAP6 positives are germline alleles of the five minor donors (VAF 0.25–16.5%), not acquired mutations. Four of the six donors (HG00438, HG02257, HG02486, HG02622) are in HPRC v1.1; all six, including HG002 and the 83.5% background HG005, are in HPRC v2.1. Our graphs keep only nodes carried by about 10% or more of their haplotypes, so they absorb the positives that are common among the graph's haplotypes; a donor's private alleles are filtered out either way. The donors' rare alleles are also found in CoLoRSdb, one of our population filters, far more often than those of our cancer benchmarks' normals (chr20 SNVs: 14.5–25.9% against 4.0–7.9%). The benchmark therefore pits the population prior, which the trained model relies on (R2.1), against its own design. [待核实（用户）: the submitted HapMap results appear to come from `mgi_non_HG005.d9`, a year-1 graph without the six donors. If so, we say so and correct the Methods; the germline-positive problem remains.]

We removed the HapMap mixture from the main benchmark and withdrew the HapMap-trained low-VAF model. Low-VAF performance is now measured, through the same pipeline, on in-silico dilutions of COLO829T with COLO829BL on three platforms: nominal purity 75, 50 and 25%, measured 0.66 / 0.40 / 0.20 (Illumina, fiber-seq) and 0.70 / 0.44 / 0.22 (ONT), biased upward at low purity. One model per platform and variant type is trained on chr2–22 of all four purities and tested on chr1 of each; the tumor itself is not held out.

- Fiber-seq SNV F1 after population filtering, at one validation threshold for all purities, is 0.779 / 0.687 / 0.665 / 0.720 at 100 / 75 / 50 / 25% purity (547–705 calls per purity without a GRCh38 position are not scored; R3.4). The share of truth SNVs that any model can reach falls from 0.821 to 0.711.
- The dip at 75 and 50% comes from the tumor's loss of heterozygosity: where the tumor lost the ALT copy of a germline heterozygous SNV, the admixed normal reads put that allele at (1 − purity)/2, inside the somatic allele-fraction range of those mixtures.
- In the Illumina and ONT mixtures, the normal reads have higher base quality than the tumor reads, a cue that a real tumor lacks, so those scores are likely optimistic. [待补: results after correcting this cue.]

[待补 / 用户决定: a second tumor held out from training: the SEQC2 HCC1395 wet-lab titration (75–5% purity, DNA mixed before sequencing) or an in-silico HCC1395/HCC1395BL dilution.]

[用户决定] (a) HapMap leaves the manuscript; or (b) the Supplement reports it only on the donor-free graph (40 HPRC individuals + GRCh38 + CHM13), stating that its positives are germline alleles. [待补]

**证据（内部）。**

- Donor 在图里：`experiments/E8b_literature_facts/README.md` H8（v1.1 有 4/6；v2.1 d46 有 6/6，含 HG002 和 HG005），并在 §5 第 1 条更正了草稿里 "four of its six donors are in both HPRC graphs" 的说法；H1：v1.1 = 44 individuals + GRCh38 + CHM13 = 90 haplotypes，HG002、HG005、NA19240 被排除。
- 图只保留常见 allele（核对时加进 Response 第一段）：E8b H3（v1.1 d9 = `--filter 9`，去掉少于 9 / 90 个 haplotype 覆盖的 node，reference path 不过滤）、H6（v2.1 d46 = 46 / 464，约 10 %）。donor 的 private allele（1–2 个 haplotype）不够 9 / 46，所以 "泄漏" 只影响接近这个门槛的 allele；常见的 positives 不管 donor 在不在图里都会被吸收。这一句是从 filter 定义推出来的，没有 HapMap positives 落在图 node 上的实际计数（见还缺什么）。
- 混合比例和 VAF 范围：E8b M1（83.5 % HG005、10 % HG02622、HG002 / HG02257 / HG02486 各 2 %、HG00438 0.5 %；positives VAF 0.25–16.5 %）；truth 来源 M2（Kong et al. graph truth v1.6，somatic 7,870,692 records，`pansoma_v2_tensors/HapMap_truth/README.txt`）；companion paper M4。
- CoLoRSdb：E8b §4.4（chr20 SNV，"rare matched %" = AC 1–2 的 allele / (不在库里的 + AC 1–2 的)）。HapMap 六个 donor（含 HG005）14.5–25.9 %，HG008-N、COLO829BL 和四个 CASTLE normal 4.0–7.9 %。
  - 核对时改过：原稿 Response 用的是 P11 的 "two to five times"。这是约数，按 §4.4 两端相除并不都落在 2–5 倍里，所以 Response 改成直接写 §4.4 的百分比。
  - 原文标的是 [inference]（consistent with membership），所以 Response 写 "found more often"，没有写 "contains"。
- Donor-free 图：E8b H10（`mgi_non_HG005.seqfile`，40 individuals + GRCh38 + CHM13 = 82 haplotypes）。我又用 `bcftools query -l /scratch/jshen/data/non_HG005_VG_Indexes/mgi_non_HG005.d9.raw.vcf.gz` 核了一遍：41 个 sample（CHM13 + 40 个），HG002 / HG00438 / HG005 / HG02257 / HG02486 / HG02622 / NA19240 都不在里面（0 个）。
- 原稿 HapMap 可能本来就跑在 donor-free 图上（新发现，E8b 没提，要用户确认）：
  - `tmp/hapmap_nonHG005_20261012/submit_prep.sh` 第 2–4 行写着 "the graph the old .dat/.idx run used"；
  - 旧训练数据准备脚本 `/scratch/jshen/Scripts/HapMap_train_data_preparation.sh` 第 22–33 行读的是 `*_mgi_non_HG005.dat/.idx`，第 67 行用的是 `hprc-v1.1-mc-grch38.nonHG005.d9.genomeAD.${chr}.json`；
  - `.dat` 文件的日期是 2025-10-23；
  - 同一目录下还有一个 `hprc-v1.1-mc-grch38.d9.genomeAD.chr1.HapMap_Merged.json`，所以原稿里每张 HapMap 图到底用的哪张图还不能确定。
- HapMap 重跑的状态（2026-10-13 的 squeue / sacct）：
  - graph_prep 387100 COMPLETED；
  - 8 个 sort 里 387086、387089 COMPLETED，387087 / 387088 / 387090 / 387091 / 387092 RUNNING，387093 PENDING；
  - tensor chain 387101–387105 都在 PENDING（Dependency）；
  - `pansoma_v2_tensors/HapMap_Illumina_nonHG005` 还不存在。
- Dilution 结果：`analysis/benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012/tables.md`。
  - fiberseq SNV PoN F1 0.779 / 0.687 / 0.665 / 0.720。
  - ONT 0.677 / 0.654 / 0.639 / 0.666。
  - Illumina 0.675 / 0.840 / 0.833 / 0.804。Illumina P100 偏低是 pooled threshold 0.804 造成的，P100 的 chr1 best 是 0.784。
  - Ceiling 和 measured purity 来自同目录的 `purity.md`：fiberseq SNV with a tensor 0.821 → 0.711；measured purity 0.66 / 0.40 / 0.20（Illumina、fiberseq）、0.70 / 0.44 / 0.22（ONT）。`purity.md` 注明这个估计在低 purity 时偏高（tensor 要求 AF ≥ 0.08），Response 照抄写 "biased upward at low purity"。
  - 这些 F1 都在 PASS（joint checkpoint 的 validation recall-0.9 阈值，四个 purity 共用一个：fiberseq 0.549）。COLO829T 的 BED 覆盖整条染色体，in-BED 和 no-BED 一样。
  - Unplaced（核对时补的）：fiberseq 每个 purity 的 PASS calls 里没有 GRCh38 位置、rtg 不评的有 705 / 547 / 627 / 688 个（P100 / P75 / P50 / P25，`runs/COLO829_fiberseq_SNV_joint_7fc7396_gpu4/vcf_chr1_*.eval_bed.txt`）。把它们全算 FP 的话，PoN 后 precision 是 0.646 / 0.589 / 0.552 / 0.590（报的是 0.881 / 0.741 / 0.715 / 0.786）。所以 Response 加了 "547–705 calls per purity ... are not scored"。
  - In-sample：训练和测试是同一个肿瘤（chr2–22 训练，chr1 测试），而且 dilution 和 P100 共用 reads（README "Data" 段），所以这不是跨肿瘤的结果；Response 加了 "the tumor itself is not held out"。
- LOH dip：同目录 README §3（标题 "reproduced and explained"）。原稿 Response 写 "is biological"，核对时改成 "comes from the tumor's loss of heterozygosity"，和 README 的说法一致（README 的依据是 per-tensor 计数和 AF 峰的位置，不是另做的实验）。
  - fiberseq 在 PASS 阈值 0.549 下被判成 somatic 的 germline tensor：P75 4,538 个（AF 峰在 0.1–0.2，2,107 个），P50 4,470 个（峰在 0.2–0.3，1,776 个）。
  - (1 − purity)/2 = 0.125 / 0.25。
- BQ artefact：同目录 README §2 和 `bq_alt_ref.json` / `bq_relative_check.json`。
  - Illumina 的两个 read group 平均 BQ：tumor.sub 29.43，normal.sub 30.77（COLO829BL normal reads 的 BQ 更高）。
  - P25 时只用这一项区分 somatic 的 AUC：Illumina 0.706，ONT 0.689；fiberseq 只有 0.439–0.459，没有这个 cue。
  - README §2 的结论是 Illumina（和 ONT）dilution 的分数 "likely optimistic for a real low-purity tumor"；benchmark 总 README 同样写 "their scores are likely optimistic. fiberseq has no such cue"。Response 用的是 "likely optimistic"。
- "the trained model relies on" 这个 prior（Response 第一段，细节在 R2.1）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §4。12 个图上的 SNV 模型（HG008T d9 ×3；HG008T 和 COLO829T 的 d46 各 ×3；COLO829T d9 purity joint ×3，在 P100 上测）在 all_H 模式（推理时把 ch6 每个有效 cell 改成 100）下，被判成 somatic 的 germline test tensor 是 base 的 3.41–7.57 倍（chr1，各模型自己的 PASS 阈值）。linear 模型的 ch6 是常数（H = 17，每个 cell 都是 100），没有这个效应。过 per-tensor PoN 以后 F1 的变化小得多（例如 COLO829_Illumina_d46 0.784 → 0.791），所以 "relies on" 说的是分类器本身，不是 PoN 之后的最终 call set（见 R2.1 证据）。这些模型都是 in-sample（同一肿瘤 chr2–22 训练、chr1 测试），不是 HapMap 上的结果。
- 两处来源不一致，都用了较新的数字：
  - 详细版 R1.5 的 fiberseq 数字 0.785 / 0.699 / 0.663 / 0.722 来自 10-06 的 tequila run；这里用 10-12 RIS run 的 0.779 / 0.687 / 0.665 / 0.720。
  - 详细版的 measured purity 0.65 / 0.38 / 0.17 是按 read 来源数的（benchmark README 里写作 "the reviewer notes"）；这里用 `purity.md` 按 AF 算的值。
- E14a（mixture 异常）已经有解释（LOH）；E14b（Illumina / ONT mixture 模型）在 RIS 上已经训完（6 个 joint run，tables.md 的 Training 表）。
- SEQC2 titration：E8b Q4–Q6（75 / 50 / 20 / 10 / 5 %，HiSeq 4000 一式三份）。本地只有 50 % 那一个 run（SRR7890936，约 148x）。
- 撤回的 HapMap 模型：详细版 R2.4 评估写 "旧的 low-VAF 模型是用 HapMap 训练的"；旧训练脚本是 `/scratch/jshen/Scripts/HapMap_SNV_train.sh`。

**还缺什么。**

- [用户决定] HapMap 是 (a) 直接去掉，还是 (b) 在 donor-free 图 `mgi_non_HG005.d9` 上放进 Supplement。重跑（chain 387101–387105）已经在排队，E23 只剩 tensor、标签和模型。
- [用户决定] 确认原稿 HapMap 的图和表是在 `mgi_non_HG005.d9` 还是 `hprc-v1.1-mc-grch38.d9` 上跑的，Methods 当时是怎么写的。这决定 Response 里的泄漏一段怎么措辞。
- [GPU RIS] 修掉 read-source BQ cue 以后，重训 Illumina / ONT 的 purity joint model（SNV + INDEL），重新报 per-purity 数字（E14b）。同时按每个 purity 自己的 validation blocks 定阈值，解决 Illumina P100 被 pooled threshold 压低的问题。
- [用户决定] 第二个肿瘤用哪种数据：E14c in-silico HCC1395 dilution（不用下载，约一周集群时间），还是 E14d SEQC2 湿实验 titration（DNA 先混再测序，按构造就没有 read-source cue [inference]；本地只有 50 % 一档，其他档要下载）。
- [GPU RIS] 如果选 (b)：HapMap_Illumina_nonHG005 的 tensor 做完以后再训练和测试（E23；用哪个模型也要用户定）。
- [CPU 本地] E8b §4.5 里 HapMap chr20 positives 被 PoN 标掉的比例，量化 "PoN 也会去掉 HapMap positives"。Slurm array 387253：sacct（2026-10-13）显示 task 0–5 COMPLETED，HapMap 是 task 6（`HapMap_positives_chr20`），还在 PENDING；做完以后用 E8b 的 `collect.sh` 汇总。
- [CPU 本地]（可选，同 R2.1）数一下 HapMap positives 落在 d9 / d46 / `mgi_non_HG005.d9` 非 reference node 上的个数（例如 chr20），量化 donor 在图里到底多吸收了多少 positives。
- [文字] Methods 补 dilution 的构建方法和 measured purity；Response 里 donor 和图的对应要按 E8b H8 写（v2.1 是 6/6）。

### R1.6 Low sensitivity

**状态：部分解决**: Pansoma 自己的 recall 漏斗（truth 有 tensor → PASS → PoN 后）、PoN 对任何 caller 的 recall 上限和 INDEL 的限制都已有 7fc7396 benchmark 的最终数字，跨样本 recall 也有了 E12a 的 preliminary 数字（旧 checkpoint）；比较器在同一协议下的 recall（user 已有的 ClairS-TO / DeepSomatic 结果待 re-score）、held-out 的 INDEL 阈值规则和用 per-dataset normalization 重训的留出细胞系模型还没有。

**Comment（要点）。**

- F1 的提升来自 precision；Supplementary Fig. 1 和 Tables 2–3 里 COLO829T、HG008T 的 SNV recall 低于 ClairS-TO 和 DeepSomatic-TO。
- INDEL recall 只有约 4–22 %（COLO829T）和 13–19 %（HG008T）。
- tumor-only 没有 matched normal 或第二阶段来补回漏掉的 variant，所以 sensitivity 很重要，F1 不能替代它。

**Response.**

We agree that sensitivity is essential in tumor-only calling.

*Operating point.* The PASS threshold is set on held-out validation blocks of the training chromosomes for a validation recall of 0.9; population filtering follows. After candidate generation only records with p < 0.01 are left out: records below the threshold are kept as LowQual and population matches are tagged in FILTER, not deleted, so users can move along the precision–recall curve, now shown for every caller.

*Where recall is lost* (chr1, benchmark regions, HPRC v1.1 graph; ranges over the three platforms; HG008T INDEL candidates include partial allele matches; Supplementary Table [待补]):

| Recall | SNV HG008T | SNV COLO829T | INDEL HG008T | INDEL COLO829T |
|---|---|---|---|---|
| truth allele reaches a candidate | 0.944–0.951 | 0.812–0.816 | 0.539–0.702 | 0.397–0.458 |
| PASS | 0.855–0.862 | 0.705–0.736 | 0.359–0.431 | 0.321–0.420 |
| PASS, after population filtering | 0.826–0.831 | 0.686–0.717 | 0.228–0.256 | 0.237–0.267 |
| any threshold, after population filtering | 0.888–0.897 | 0.790–0.795 | 0.238–0.259 | 0.237–0.275 |
| bound for any caller using the same filter (whole chr1 truth) | 0.913 | 0.966 | 0.350 | 0.382 |

After the candidate stage, SNV recall depends mainly on the operating point, whereas INDEL recall is bounded by the population filter (below). Candidate-stage SNV losses are graph SNP alleles and repeat absorption on HG008T, most of which population filtering would tag anyway (R1.4, R2.3), and low-VAF subclonal mutations below the candidate allele-fraction or read-support thresholds and regions without confident mapping on COLO829T (R2.3, R2.4).

*INDELs.* The Abstract and Discussion now state this limitation (R3.10). Two limits lie outside the classifier. First, 413 of the 635 HG008T chr1 truth INDELs (81 of 131 on COLO829T) are exact population alleles at the filtering frequencies, so no caller filtered this way exceeds a recall of 0.350 (0.382); our post-filter recall at PASS is close to the curve's ceiling, and long-read INDEL accuracy is limited by precision. Second, the graph absorbs repeat-length INDELs (R1.4): a GRCh38-only graph gives more HG008T truth INDELs a candidate (0.75–0.87 against 0.54–0.70), but these absorbed alleles are almost all population alleles, so its post-filter ceiling is only 0.284 / 0.301 / 0.238 against 0.257 / 0.259 / 0.238. Both limits depend on the truth set: on HCC1395 (SEQC2; chr1–22, all records), only 119 of 1,625 truth INDELs are population-tagged, and 0.76–0.85 reach a candidate on the HPRC v2.1 graph (HG008T: 0.43–0.50).

*Other tumors (preliminary; earlier checkpoints, HPRC v2.1 graph, chr1).* Without retraining, at the threshold chosen on the training tumor, the HG008T SNV models reach a post-filter PASS recall of 0.640 / 0.701 / 0.589 on COLO829T (Illumina / fiberseq / ONT), against 0.715 / 0.707 / 0.692 for COLO829T's own models. The COLO829T models lose more at PASS, mostly through threshold calibration (HCC1395T Illumina: recall 0.157 at precision 0.934; test-tuned best F1 0.810), which per-dataset quality normalization targets (R1.7). [待补: retrained models]

*Comparators.* ClairS-TO and DeepSomatic-TO are scored under the identical protocol (R3.m2): recall at their own PASS, after the shared population filter, and at Pansoma's precision along their QUAL curves. [待补: comparator recall at each stage; INDEL numbers at a validation-chosen threshold; the funnel with the retrained held-out models (R1.7)]

**证据（内部）。**

- HG008T SNV / INDEL（d9，7fc7396）：`analysis/benchmark_7fc7396_20261012/hg008t_snv_indel_7fc7396_20261012/tables.md`。
  - SNV truth 有 tensor 659 / 658 / 663 of 697（0.945 / 0.944 / 0.951，Illumina / PacBio / ONT）；raw PASS R 0.862 / 0.857 / 0.855；PoN 后 PASS R 0.831 / 0.828 / 0.826（P 0.452 / 0.497 / 0.466）；PoN ceiling 0.888 / 0.888 / 0.897；PoN 标掉 61 of 702 chr1 truth SNV。
  - INDEL truth 有 tensor 317 / 413 / 408 of 588（0.539 / 0.702 / 0.694，predict 的 truth report，partial match 算在内）；raw PASS R 0.359 / 0.431 / 0.400；PoN 后 PASS R 0.238 / 0.256 / 0.228；PoN ceiling 0.257 / 0.259 / 0.238；PoN 标掉 413 of 635；chr1 best F1 0.264 / 0.260 / 0.140，对应 R 0.213 / 0.204 / 0.109。
- COLO829T（d9）：`analysis/reviewer_response_20261011/experiments/E8_comparators/tables.md` 的 Pansoma 行（来自 `benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/summary.json`）。
  - SNV（`COLO829_*_SNV_h90`，10-03，1 GPU，同 SNV recipe）：truth 有 tensor 0.816 / 0.816 / 0.812（Illumina / fiberseq / ONT，`hprc_v21_d46_7fc7396_20261012/tables.md`）；raw PASS R 0.729 / 0.736 / 0.705；PoN 后 R 0.711 / 0.717 / 0.686（P 0.860 / 0.857 / 0.727）；PoN ceiling 0.795 / 0.792 / 0.790；PoN 标掉 84 of 2,443。
  - INDEL（`COLO829_*_INDEL_small`，10-06，lr 4e-4，用 7fc7396 PoN 重新打分）：truth 有 tensor 0.427 / 0.458 / 0.397（10-06 的 predict 不算 partial match，和 HG008T 不完全可比）；raw PASS R 0.382 / 0.420 / 0.321；PoN 后 0.267 / 0.244 / 0.237；PoN ceiling 0.275 / 0.244 / 0.237；PoN 标掉 81 of 131。
- PoN bound（"bound for any caller using the same filter"）：E8 `tables.md` 的 "PoN bound on recall" 列：HG008T SNV 0.913、COLO829T SNV 0.966、HG008T INDEL 0.350、COLO829T INDEL 0.382。
  - 计算方式（`E8_comparators/collect.py` L249，`linear_grch38_7fc7396_20261012/summary.json` 的 `truth_records` / `truth_pon_tagged`）：1 − tagged / truth records，truth 用整条 chr1（HG008T 702 SNV / 635 INDEL，eval_bed 和 eval_nobed 相同），不是 BED 内的 697 / 587；COLO829T 的 2,443 / 131 在两种 eval 里也相同。所以表里这一行标 "whole chr1 truth"，其余行是 BED 内。
- GRCh38-only 臂：`benchmark_7fc7396_20261012/README.md`（"truth with a tensor 0.75-0.87 vs 0.54-0.70 on HG008T"；`linear_grch38_7fc7396_20261012/tables.md`：0.827 / 0.872 / 0.748 vs d9 0.539 / 0.702 / 0.694）；PoN ceiling linear 0.284 / 0.301 / 0.238 vs d9 0.257 / 0.259 / 0.238（E8 `tables.md` 的 curve ceiling after PoN 列；`experiments/E3_reference_arms/README.md` 同）。
- "almost all population alleles" / "most of which population filtering would tag anyway"：`tmp/graph_absorbed_pon_20261012/README.md`，被图吸收（perfect bypass）的 HG008T truth 被 PoN 标掉 SNV 199 / 201、INDEL 2,156 / 2,277（94.7 %），三个平台的 union，chr1–22；`E5b_absorbed_d46/tables.md` §2：d9 比 linear 少了 tensor 的 HG008T SNV 293 / 333 / 317（Illumina / PacBio / ONT）里 PoN 标掉 87.7 % / 90.4 % / 90.2 %（chr1–22）。所以是 "most"，不是全部。
- HCC1395：`experiments/E5b_absorbed_d46/tables.md` §1，HCC1395T INDEL n 1,625（not tagged 1,506，即标掉 119）；同表 HG008T INDEL 8,496（3,251），COLO829T 1,912（929），chr1-22。有 tensor 的比例（d46，chr1–22，全部记录）：HCC1395T INDEL 1,346（0.828）/ 1,375（0.846）/ 1,239（0.762）（Illumina / PacBio / ONT），HG008T INDEL 3,821（0.450）/ 4,226（0.497）/ 3,642（0.429）→ Response 的 "0.76–0.85 ... (HG008T: 0.43–0.50)"。该文件标 PRELIMINARY，但只影响 absorbed 列（graph spelling 还在跑），§1 的 tensor / PoN 列不受影响。E8b 的 Slurm 387253 原始输出 `/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts/somatic_truth_in_resources.HCC1395_SEQC2_v1.2.1.tsv` 也是 Pansoma 规则标 119 of 1,625（全部记录，没有 BED），一致。
- 跨样本（E12a，preliminary）：`experiments/E12a_cross_sample/tables.md` 和 `README.md`（d46，chr1，在 target 的 BED 内，PoN 后，阈值 = 模型自己在训练 tumor 上选的 validation R0.9 阈值；输出 `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/xs__<model>__<set>/`）。
  - HG008T 模型 → COLO829T（Illumina / PacBio→fiberseq / ONT）：P 0.792 / 0.817 / 0.678，R 0.640 / 0.701 / 0.589，F1 0.708 / 0.754 / 0.630（best 0.717 / 0.764 / 0.639）；→ HCC1395T：R 0.778 / 0.624 / 0.628，F1 0.790 / 0.698 / 0.689。
  - COLO829T 模型 → HG008T：R 0.255 / 0.548 / 0.420，F1 0.359 / 0.557 / 0.487（best 0.554 / 0.597 / 0.579）；→ HCC1395T：R 0.157 / 0.646 / 0.515，Illumina 那行 thr 0.640、P 0.934、best F1 0.810。best F1 是 chr1 曲线上 test-tuned 的点（E12a `tables.md` 表头 "test-tuned, reference only"），Response 已这样标。
  - unplaced PASS calls（不进 rtg，不在 F1 里）：HG008T 模型 → COLO829T 225 / 550 / 1,388。
  - in-sample d46（`benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/tables.md`）：COLO829T PoN 后 R 0.715 / 0.707 / 0.692（F1 0.783 / 0.779 / 0.737）；HG008T R 0.852 / 0.839 / 0.859（F1 0.570 / 0.605 / 0.575）。
  - 限制（E12a README Caveats）：旧 checkpoint，BQ / MAPQ 用训练 tumor 的统计量 z-score（COLO829T Illumina BQ 中位数 30，HG008T 是 4 级 binned BQ）；一个 seed、一条测试染色体；HCC1395T 的 truth 是 short-read 来源，INDEL 不跨 target 比较；还没有 HCC1395T in-sample 模型。
- Pansoma 的 chr1 call VCF：tequila GPU 上重新 predict 的三个 d46 SNV run（`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/vcf_chr1/`，含 `eval_bed/vcfeval/pon_calls/{tp,fp,fn}.vcf.gz`）复现 RIS：`eval_bed/report.txt` PoN F1 0.603 / 0.575 / 0.736，RIS 0.605 / 0.575 / 0.737。
- 记录保留：`machine_learning/pansoma_net_v2/graph_vcf.py` docstring L33–46（LowQual；"Tensors with p_somatic < --min-score (default 0.01) that are not called are not written"，predict 标 excluded 的 tensor 也不写）；`scripts/filter_panel_of_normals.py`（默认保留记录加 FILTER，`--drop-matched` 才删）。旧稿 "No record is deleted after candidate generation" 太强（p < 0.01 的不写出），Response 已改成 "only records with p < 0.01 are left out"。
- 漏检原因：`analysis/tensor_recall_20260930/README.md` §2（HG008T SNV 的 S3 + S4 占 no-candidate 的 72–89 %；COLO829T 是 S1 Extreme 区域和 S2 低 VAF）；`tensor_recall_20260930/tables.md` §7：COLO829T fiberseq chr1 2,443 个 truth 里缺 380，其中 AF 过滤 271、no candidate 40（09-30 的 build AF 0.06，现在是 0.08，所以只能定性引用）。
- INDEL 的 PoN 规则差异：`experiments/E8b_literature_facts/README.md` §4.5，HG008T 全部 somatic INDEL 里 Pansoma 规则标 5,783（61.8 %），ClairS-TO 默认规则标 8,155（87.1 %）。
- 和旧稿的出入（用新数字）：详版 R1.6 表写 INDEL 有 candidate 36–59 %、PoN 后 23–25 %；简版写 COLO829T SNV R 0.72–0.73、"72–84 % of truth INDELs are exact population alleles"（09-30 的旧 PoN 规则）。7fc7396 规则下全染色体是 HG008T 61.8 %、COLO829T 52.0 %（`tmp/graph_absorbed_pon_20261012/README.md`）。

**还缺什么。**

- [数据同步] user 提供已有 ClairS-TO / DeepSomatic 结果的路径，然后用 `E8_comparators/score_comparator.py` 按同一协议 re-score（E8 re-score；不重跑）：比较器的 PASS / PASS+PoN / 曲线 recall 和 matched-precision recall。
- [CPU 本地] E1.6 漏斗脚本：每个 truth allele 的 tensor / PASS / PoN / 最终状态，按 SNV / INS / DEL、VAF bin、repeat 类别分层，比较器只有 truth → PASS → PoN；出 Supplementary Table。三个 d46 SNV run（HG008T PacBio / ONT、COLO829T ONT）的 call VCF 和 rtg tp / fp / fn 已在 `reviewer_cpu_20261013/`，可以先做；E12a 的 `xs__` 输出也有 VCF，可以顺便给跨样本漏斗。
- [数据同步] 表里 d9 run 和所有 INDEL run 的 chr1 call VCF 要从 RIS 同步（本地 `ris_7fc7396_20261012/<run>/test_chr1/` 只有 metrics）。
- [数据同步] INDEL 的 held-out 阈值规则（在 validation 上选阈值）需要从 RIS 同步 predictions；在此之前 INDEL 只能报 PASS（R0.9）和 test-tuned best F1（上界）。
- [GPU RIS] E12 跨样本 / 留出细胞系模型（新默认 `--norm dataset`，commit f2a0e4e，原 cbc94a7）、leave-one-line-out、HCC1395T in-sample 模型，给留出细胞系的最终漏斗；E12a（旧 checkpoint）只算 preliminary。
- [GPU RIS] COLO829T d9 INDEL 按 7fc7396 recipe（lr 1e-4）重训，让 COLO829T INDEL 漏斗和 HG008T 同一 recipe、同一 partial-match 计数。
- [用户决定] 漏斗的主图（d9 还是 d46）和主评估（同一样本 chr1 还是留出细胞系）；E12a 的 preliminary 跨样本 recall 是现在就写进 Response（如上），还是等 `--norm dataset` 重训的模型。
- [文字] Abstract / Discussion 写明 SNV 和 INDEL 分开、INDEL 的限制和 PoN 上限；Methods 写 operating point 和 LowQual / FILTER 保留。

### R1.7 Training scale and tumor diversity

**状态：部分解决**: 参数量、每个模型每类的训练样本数、overfit 和 early stopping、两株训练细胞系的 mutation spectrum 都已有最终数字（E21）；"模型不是 tumor-specific" 现在有了直接但 preliminary 的证据（E12a cross-sample，一个 seed：HG008T 模型在 COLO829T 上的 PASS F1 比 COLO829T 自己的 in-sample 模型低 0.025–0.107，在 HCC1395T 上 0.689–0.790（没有 in-sample 参照）；COLO829T 模型主要掉在阈值上，best F1 是 test-tuned upper bound），足以支撑主论点；但 E12a 用的是旧 normalization 的 checkpoint、一个 seed，最终的 `--norm dataset` 模型、leave-one-line-out、HCC1395T in-sample 参照和 learning curve（E13）都是 RIS 上的 GPU job，还没有。

**Comment（要点）。**

- 约 200 M 参数的 ConvNeXt 主要用 COLO829（40–50 k truth variants）训练，数据量够不够？
- 报告确切参数量，以及每个平台、每种 variant type 的正 / 负训练样本数，并说明数据量撑得住模型容量。
- COLO829 是单一的 UV 主导的 melanoma，模型可能学到 tumor-specific 特征；DeepSomatic 用了多个细胞系并做了大量 augmentation。
- 把训练扩展到不同 mutational process 的肿瘤，并直接证明模型不是 tumor-specific。

**Response.**

*Model size.* The revised SNV network has 50.4 M parameters (ConvNeXt, depths 3-3-27-3, widths 96–768), the INDEL network 8.7 M (2-2-6-2, 64–512). The original 199.6 M configuration overfitted (HG008T PacBio validation loss 0.51 to 2.16) and was not more accurate (F1 after population filtering 0.604, against 0.616).

*Training examples.* Supplementary Table [待补: number] lists the examples per class for each of the 42 models (one tumor line each; chr1 held out). On the pangenome graphs, SNV models see 7,098–7,274 (HG008T) or 33,472–34,110 (COLO829T) somatic, 256,583–381,159 germline and 67,735–1,620,829 non-variant examples; INDEL models see 3,121–3,761 (HG008T) and 971–1,141 (COLO829T) somatic examples. Every positive is an exact allele match, and each SNV positive is a different truth allele.

*Capacity.* With about 7,000 parameters per HG008T somatic SNV example, the network can fit its training set (final training loss 0.002–0.139), so generalisation rests on early stopping on held-out node blocks (checkpoints from epochs 4–12 of 20). A learning curve and a model-size series test whether accuracy is still data-limited [待补: E13].

*Tumor diversity.* COLO829T positives are 71 % C>T, 65 % at dipyrimidines (UV); HG008T positives are 41 % C>T, 19 % C>A and 8.4 % C>G at TCW. We added two breast (HCC1395, HCC1954) and two lung adenocarcinoma lines (H1437, H2009), each sequenced on Illumina, PacBio HiFi and ONT and processed identically on HPRC v2.1 (311,921 somatic truth alleles on chr1–22) [待补: spectra of the added lines; the CASTLE Illumina sets and H1437 ONT are still being processed].

*Tumor specificity, a first direct test.* We applied each HPRC v2.1 SNV model to chr1 of tumors it was not trained on, at the threshold chosen on its own training tumor (PASS F1 after population filtering, calls with a GRCh38 position, one seed; Supplementary Table [待补]):

- Models trained on HG008T (pancreatic; about 7,200 somatic examples) reach 0.708, 0.754 and 0.630 on COLO829T (Illumina, PacBio HiFi / Fiber-seq, ONT), against 0.783, 0.779 and 0.737 for COLO829T's own models, and 0.790, 0.698 and 0.689 on HCC1395.
- Models trained on COLO829T reach 0.359, 0.557 and 0.487 on HG008T (its own models: 0.570, 0.605, 0.575). Their ranking transfers better than their threshold: on HCC1395 Illumina, recall is 0.157 at precision 0.934, while the best point of the chr1 curve (chosen on the test data, an upper bound) reaches F1 0.810. All these models standardised base and mapping qualities with the training tumor's statistics (mean base quality 33.2 vs 28.5 in the HG008T and COLO829T Illumina training data); the final models standardise them per dataset [待补: E12 numbers with these models].
- Leave-one-line-out models, an HCC1395 in-sample reference and recall by mutational context complete the test [待补: E12].

*Augmentation.* Joint COLO829T models add in-silico purity dilutions (75, 50, 25 %); a base-quality cue of the normal-derived reads in the Illumina and ONT dilutions is removed before reporting them [待补].

**证据（内部）。**

- 参数量：`experiments/E21_training_composition/README.md` Summary（50.4 M，depths 3-3-27-3，widths 96-768；8.7 M，2-2-6-2，64-512；`train.log` 打印）。199.6 M 的 overfit：`machine_learning/pansoma_net_v2/README.md` 237–239 行（HG008T PacBio chr1，12 epochs，validation loss 0.51 -> 2.16；PoN 后 BED 内 P/R/F1 0.473/0.832/0.604 vs 50.4 M 0.493/0.822/0.616）。
- 每个模型的训练样本数：E21 README §1 表（42 个 run）和 "Reading"（graph arms：HG008T SNV somatic 7,098–7,274、germline 256,583–326,341、non 67,735–1,620,829；COLO829T SNV somatic 33,472–34,110；INDEL HG008T 3,121–3,761、COLO829T 971–1,141）；germline 范围 256,583–381,159 来自 E21 Summary（两株的 pangenome 模型）。replay check：18 / 18 一致（d9 全部和 HG008T d46，`runs_replay_check.tsv`），其余 24 个等 collect job。
- 与 detailed doc 的差异：detailed doc 写的 "HG008T SNV 8.8–10.3 k somatic、0.31–0.37 M germline" 是 label manifest 在 chr1–22 的全部 tensor；E21（较新）给的是模型真正训练用的数（chr2–22、去掉 validation block 和 −1），Response 用 E21 的数。
- 每个 SNV 正样本是不同的 truth allele：E21 Summary（"Each SNV positive is a different truth allele"；HG008T PacBio d9：7,098 tensors = 7,098 alleles）。INDEL 不是：E21 §1 / §2 训练 somatic tensor 比背后的 truth allele 略多（HG008T 3,121–3,761 vs 3,115–3,728，COLO829T d9 971–1,078 vs 967–1,072），所以 Response 改成 "each SNV positive is a different truth allele"（原稿写 "Every positive is ... a distinct truth allele"，对 INDEL 不成立）；跨平台 union：HG008T SNV 7,203（d9）/ 7,337（d46），COLO829T 34,508（d9）（E21 §2）。
- Capacity：E21 Summary（约 7,000 parameters / HG008T SNV positive；training loss 第 20 epoch 0.002–0.139，0.002 是 `COLO829_Illumina_SNV_joint_7fc7396_gpu4`，tables.md T1；validation loss 在 epoch 2–11 最低，最后高 1.4–9.5 倍；checkpoint 在 SNV epoch 4–12 / 20，INDEL 2–7 / 12）。
- Spectrum：E21 §3b、`somatic_spectrum.tsv`（COLO829T 71 % C>T、65 % C>T at dipyrimidine、0.8 % C>G at TCW；HG008T 41 % C>T、29 % at dipyrimidine、19 % C>A、8.4 % C>G at TCW、9.4 % CpG C>T）。是 context 计数，不是 signature fit（E18）。
- 新细胞系 truth：E21 §5（HCC1395T 41,072、HCC1954T 25,014、H1437T 84,744、H2009T 161,091 on chr1–22，合计 311,921，含其他 contig 323,516；COLO829T 43,947（46,064）、HG008T 17,186）。
- 新细胞系 tensor（`/scratch/jshen/data/pansoma_v2_tensors/<SET>/tensors/truth_recall.json`，somatic SNV 有 representative tensor 的数，百分比是我按文件数字算的）：HCC1395T Illumina 35,266 / PacBio 34,374 / ONT 30,869（of 39,447）；HCC1954T PacBio 20,762 / 22,622（91.8 %）；H1437T PacBio 56,703 / 79,974（70.9 %，异常未查清）；H2009T PacBio 143,877 / 154,851（92.9 %）；HCC1954T ONT 20,695、H2009T ONT 145,007（QUAL10 relabel 之后的 truth_recall.json，2026-10-13 12:27 / 12:55，somatic 数没变）。
- 新细胞系的进度（squeue 2026-10-13 15:39）：HCC1395T 三平台、CASTLE 三株 PacBio、HCC1954T / H2009T ONT 已 finalize 并 QUAL10 relabel；H1437T ONT finalize（385470）在跑；H2009T Illumina run（385782）在跑；HCC1954T / H1437T Illumina 在 gamsort（385759、385769），后面的 discover / run / finalize 和 relabel（387131–387134）pending。
- **Cross-sample（E12a，preliminary）**：`experiments/E12a_cross_sample/README.md`、`tables.md`、`summary.json`。d46 的 2026-10-12 RIS checkpoint，预测另一个肿瘤 d46 tensor 的 chr1（PacBio ↔ fiberseq，HCC1395T 用同平台），16 个 cross-sample 预测都在 tequila GPU bf16 上（array 387379 task 3–18；CPU array 387323 的 xs task 3–4 在 11:46 取消，`rcpu_predict_vcf-387323_{3,4}.out`），然后 `jobs/vcf.sh`：graph_vcf 用模型自己训练肿瘤上的 validation recall-0.9 阈值（不在 target 上调），PoN 后、target 的 BED 内 rtg vcfeval。输出 `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/xs__<model>__<target>/`。
  - SNV PASS F1（best F1）：HG008T -> COLO829T 0.708（0.717）/ 0.754（0.764）/ 0.630（0.639），COLO829T in-sample 0.783（0.786）/ 0.779（0.783）/ 0.737（0.739）；PASS 差值 0.075 / 0.025 / 0.107 是我按表相减（E12a README 写成 "0.02–0.11"）。HG008T -> HCC1395T 0.790（0.809）/ 0.698（0.775）/ 0.689（0.699）。COLO829T -> HG008T 0.359（0.554）/ 0.557（0.597）/ 0.487（0.579），HG008T in-sample 0.570（0.613）/ 0.605（0.622）/ 0.575（0.604）。COLO829T -> HCC1395T 0.269（0.810）/ 0.741（0.819）/ 0.621（0.691）。顺序都是 Illumina / PacBio-fiberseq / ONT。
  - COLO829T Illumina -> HCC1395T：阈值 0.640，PoN 后 TP / FP / FN 540 / 38 / 2900，P / R 0.934 / 0.157（`tables.md`）。best F1 是 chr1 曲线上选的点（test-tuned），Response 只用来说明 ranking 能迁移。
  - Illumina INDEL（E12a README）：HG008T -> COLO829T 0.076（0.135），in-sample 0.109（0.163）；COLO829T -> HG008T 0.254（0.262），in-sample 0.244（0.249）；-> HCC1395T 0.174（0.247）/ 0.307（0.331）。HCC1395T 的 INDEL truth 缺很多 homopolymer INDEL（E12a Caveat 4 "misses many"；rule (a)），跨 target 不可比；Response 只写 SNV。
  - Preliminary 的原因：E12a Caveats 1–2、5（旧 checkpoint，BQ / MAPQ 用训练肿瘤的统计 z-score；一个 seed、一个测试染色体、没有 CI；还没有 HCC1395T in-sample 模型）。Unplaced PASS calls 不在 F1 里（Caveat 3）。
  - "about 7,200 somatic examples"：E21 graph arms 的 7,098–7,274；E12a README 写 "about 7,200 ... against about 34,000 for COLO829T"。
  - 同一套本地预测 + `vcf.sh` 流程可信：三个 d46 run 的 in-sample 重预测（HG008T PacBio 在 CPU fp32 上，tsingtao job 387135，`rcpu_test_HG008_PacBio_SNV_d46-387135.out`；HG008T ONT 和 COLO829T ONT 在 tequila GPU bf16 上，387379 task 1–2），rtg PoN F1 0.603 / 0.575 / 0.736（`reviewer_cpu_20261013/{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/vcf_chr1/eval_bed/report.txt`），RIS 0.605 / 0.575 / 0.737（`analysis/benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/tables.md`）。
- Normalization：commit f2a0e4e（history 重写前是 cbc94a7；`train --norm dataset`，现为默认）；`pansoma_net_v2/README.md` 23–33 行（BQ / MAPQ per dataset，path count 用训练统计）。BQ 尺度差异：checkpoint 统计 HG008_Illumina_SNV_d46 BQ mean 33.161 / std 8.400（`tmp/reviewer_cpu_20261013/logs/gpu_predict_vcf-387379_3.out`），COLO829_Illumina_SNV_d46 BQ mean 28.545 / std 5.075（`rcpu_predict_vcf-387323_4.out`）；Response 里四舍五入写成 33.2 vs 28.5（checkpoint `stats` 是在训练 tensor 上拟合的）；E12a README：COLO829T Illumina median BQ 30，HG008T 4-level binned BQ。MAPQ std 2.3（HG008T PacBio）vs 8.9（COLO829T Illumina）不在 E12a README 里，出处是 `machine_learning/pansoma_net_v2/README.md` 27–29 行和 f2a0e4e 的 commit message；这个 2.3 和 `HG008_PacBio_SNV_d46` checkpoint 的训练 MAPQ std 11.747（`gpu_predict_vcf-387379_*` 日志）对不上，口径没查，Response 不用 MAPQ 数。
- DeepSomatic 训练：`experiments/E8b_literature_facts/README.md` D4（HCC1395、HCC1937、HCC1954、H1437、H2009 五株）、D9（purity / tumor-in-normal titration 进训练）。ClairS-TO：C2（SS = GIAB 合成数据）、C3（SSRS 在 HCC1937 / HCC1954 / H1437 / H2009 上 fine-tune）。
- Dilution 的 BQ artefact：`analysis/benchmark_7fc7396_20261012/README.md`（Illumina read groups mean BQ 30.77 vs 29.43；P25 时这一差值单独的 AUC 最高 0.71 Illumina、0.69 ONT；fiberseq 没有）。

**还缺什么。**

- [GPU RIS] 用 `--norm dataset`（f2a0e4e 起的默认）重训 d46 单细胞系 SNV 模型（HG008T、COLO829T 各三平台），再在 tequila GPU 上重跑 E12a 的 16 个 cross-sample 预测，替换 Response 里的 preliminary 数字（用户已定：最终 cross-sample 模型在 RIS 上用新默认训练）。
- [GPU RIS] E12 leave-one-line-out（d46）：先做每个平台 "留出 HG008T" 和 "留出 COLO829T" 两个 fold 的 SNV（6 run），再补全六株 LOLO 和 INDEL。
- [GPU RIS] HCC1395T 三平台的 in-sample 模型（chr2–22 训练），作为 HCC1395T chr1 cross-sample 数字的参照（E12a Caveat 5）。
- [数据同步] 上面这些新 checkpoint 从 RIS 拷到本地，供 tequila 上的 cross-sample 预测用。
- [CPU 本地] 用 E12a 已有的 vcfeval 输出（`xs__*/vcf_chr1/eval_bed/vcfeval/pon_calls/` 的 tp / fn）按 mutational context 分层（C>T at dipyrimidine、TCW、C>A）：比如 HG008T 模型在 COLO829T UV 突变上的 recall 对 COLO829T 自己的模型。不用新 job，先在 preliminary 结果上做，最终模型出来后重跑。
- [GPU RIS] E13：8 个 learning-curve run + 2 个 model-size run（命令在 E21 README §6）。
- [数据同步] 把 E13 的 exclude 文件（`/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/e13/`，16 个，70 MB）拷到 RIS。
- [CPU 本地] E21 collect（387193，连同 index array 387137 的后续 task，现为 held）放出后补 HCC1395T / CASTLE 的 spectrum 和 VAF；签名分层要 E18 的 SBS96 脚本。
- [CPU 本地] CASTLE ONT / Illumina 集合 finalize 和 QUAL10 relabel（H1437T ONT 和三株 Illumina 的 pipeline job 在队列里）。
- [用户决定] cross-sample 的阈值怎么报：E12a 显示 COLO829T 模型的阈值迁移不过去（Illumina -> HCC1395T recall 0.157，best F1 0.810）；`--norm dataset` 若修不好，LOLO 是否用多株的 validation block 选阈值，主结果报迁移阈值下的 PASS、best F1 只作 upper bound。
- [用户决定] H1437 PacBio 异常（70.9 %）查清前是否从 E12 排除；联合训练要不要加 per-set 权重（`train.py` 现在没有这个选项，要改 pansoma_net_v2）；主结果是否用 cross-sample（E12）、同一样本内只作上界。
- [GPU RIS] Illumina / ONT dilution 的 BQ artefact 修好后重训 joint mixture 模型（用户在做），否则 Response 里 augmentation 一句要改成只提 fiberseq。
- [文字] 用 E21 tables.md T1 做 Supplementary Table；参数量和 early stopping 写进 Methods；最终模型的 cross-sample 表（E12a 的格式）做另一张 Supplementary Table。

### R1.8 Construction of the training set

**状态：部分解决**: 标签规则（正样本投影、partial、off-reference、"determinable"、germline、artefact）和 d9 臂 + HG008T d46 的训练集覆盖率、VAF 分布都已有最终的代码出处和数字（E21 §1–4、`truth_labels.py`），off-reference tensor 不进训练的后果现在也量化了（E1.1 §5：没有 rescue 时 assembly 验证的 off-reference somatic SNV 只有 0 / 112、1 / 89 过阈值），Response 可以定稿；但 E21 的 collect job（387193）仍被 hold，COLO829T d46、linear 臂、dilution 和 HCC1395T / CASTLE 的覆盖率与 VAF 行还没出，per-truth 状态表也还没整理成 Supplementary Table。

**Comment（要点）。**

- 正样本：GRCh38 上的 truth 投影到 graph 坐标并按 exact allele 匹配；truth allele 本身就是 graph path、投影有歧义、graph 和 linear 表示不同（尤其 INDEL）时，正样本怎么定？
- graph-only 位置通过相邻的 1-bp SNV node 间接得到 GRCh38 坐标。
- 被吸收的位点永远成不了正样本（Supp. Fig. 3）。
- 要明确写出标签和投影过程，包括 graph 一侧是否做 normalisation，并说明正样本集合没有被系统性扭曲。
- 负样本：什么叫 "determinable"？germline 是否显式包含？artefact 怎么定义？VAF 和 purity 怎么取样？

**Response.**

Labelling is now described in full (Methods, Supplementary Note; every label manifest records the rule file's checksum).

*Positives and projection.* Each truth ALT is split, trimmed and, in repeats, expanded to every equivalent placement (up to 5,000). Each placement is converted along the graph's GRCh38 path into the builder's node key. Because all placements are enumerated, no separate normalisation is needed on either side. A tensor is somatic only when its representative allele equals a truth key. Where graph and linear representations differ, matches with allele or read-haplotype overlap above 0.45 are recorded as partial and excluded from training, so every training positive is an exact allele match (in an earlier model version, this raised HG008T Illumina's best chr1 F1 after population filtering from 0.623 to 0.665).

*Graph-only positions* are no longer placed through neighbouring nodes; without a truth match, branch-node tensors are excluded from training. The model thus never sees a positive off GRCh38 and scores a low haplotype path count as not somatic: of the HG008T somatic SNVs on off-reference nodes of HPRC v1.1 (from the HG008-N assembly), 0 of 112 (PacBio) and 1 of 89 (ONT) pass the threshold. Prediction re-scores off-reference SNV tensors with the full haplotype count at the site (32 of 112 and 46 of 89 called), reported separately as calls without a GRCh38 position (R3.3, R3.4).

*Is the positive set distorted?* Alleles that are already graph paths produce no candidate (R1.4), but on HPRC v1.1 94.7–99.0 % of absorbed truth alleles are population-tagged anyway, so the graph-specific loss is 2 SNVs and 121 INDELs for HG008T (1.3 % of its INDELs) and 1 SNV and 13 INDELs for COLO829T (0.6 %). In the training region, 94.8–95.9 % (HG008T, both graphs) and 89.5–90.4 % (COLO829T, HPRC v1.1) of truth SNVs stand behind a training example; the COLO829T shortfall is mostly subclonal SNVs below the candidate threshold. On chr1, 94.4–95.4 % of HG008T truth SNVs have a tensor on the pangenomes, against 97.6–98.0 % on GRCh38 only. INDELs are affected more (42.6–63.6 % reach training). Supplementary Table [待补] gives every truth allele's status.

*Negatives.* Germline is an explicit class from the matched normal's truth (HG008-N and COLO829BL assembly calls; DeepVariant normal calls with QUAL ≥ 10 for SEQC2 / CASTLE): 256,583–381,159 germline SNV examples per pangenome model. Non-variant examples are tensors on a GRCh38 node inside the confident region without a matching truth allele, including errors next to true variants. "Determinable" means exactly this: inside the confident region, at a unique GRCh38 position, and either an exact truth allele or none; all other tensors are ignored.

*VAF and purity* are not resampled: each line keeps its natural distribution (median somatic SNV allele fraction 0.53–0.54 for HG008T, 0.46–0.51 for COLO829T), plus COLO829T purity dilutions in joint models.

**证据（内部）。**

- 规则：`indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py` docstring 1–58 行（展开、node key、label 1 / 2 / 0 / −1、partial）和 `classify` 655–690 行（off-reference node 用参考邻居之间的区间找 partial；没有 → `off_reference_no_truth_match` −1）；中文解释 `analysis/truth_labels.md`；E21 README §4（`MAX_SHIFTS` 5,000、`MIN_OVERLAP` 0.45、规则行号）。
- Partial 不进训练：E21 §4 第 2 条（`ris/recipes.sh` 6、9 行 `--ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth`；`train.py` 337–346 行）。每个模型去掉的 partial 数在 E21 §1 表 "left out" 列（如 HG008T PacBio SNV d9 464 / 28）。
- 0.623 -> 0.665：`analysis/v1_vs_v2_20260929/synthesis/ablation_check.out`（`v2 keephard_w100` vs `keephard_w100_nopartial`，HG008T Illumina SNV，repo PoN 后 max F1，2026-09-29 的旧 build）。这是 chr1 曲线上的最好点（test-tuned），detailed doc 写成 "chr1 F1"，Response 已改成 "best chr1 F1"，并注明是 "an earlier model version"（不是 7fc7396 的模型）。
- Off-reference 不进训练的后果（新，E1.1）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §5（HG008T assembly off-reference somatic SNV：HG008-N assembly 事件里 somatic allele 被 giraffe 放到 d9 的 off-GRCh38 node 上、且有 tensor 的，全部常染色体；E1.1 README："None of them is a training label (their labels are -1 or an ignored partial match)"）。格子是 p ≥ t / argmax somatic / n：d9 `HG008_PacBio_SNV_7fc7396` base 0 / 0 / 112，`HG008_ONT_SNV_7fc7396` base 1 / 1 / 89；`site10_high`（E1.1 README："`predict`'s off-reference rescue transform"，site ±10 列 ch6 低于 100 的 cell 提到 100）45 / 32 / 112、45 / 46 / 89；pipeline 对 off-reference SNV tensor 用 rescue 后的 argmax（`offref_call`，`machine_learning/pansoma_net_v2/README.md` 221–230 行、`predict.py` 31–38 行），所以 Response 写 32 和 46。`all_H` 是 61 / 44 / 112、49 / 49 / 89。只有 d9：d46 的 assembly key 已分好类（`asm_d46_classify.py`），但用 d46 run 打分（`asm_scores.py`）还没跑（E1.1 README）。
- "learned that a low haplotype path count is not somatic"：`pansoma_net_v2/README.md` 221–223 行原话；E1.1 §1 site code median somatic 100、germline 95–96（d9）/ 74–80（d46）；§6 只用 ch6 区分 somatic / germline 的 AUC 0.837–0.892。训练里 `off_reference_no_truth_match` 是 −1（同一 README 147 行）。
- 吸收位点：`tmp/graph_absorbed_pon_20261012/README.md`（HG008T SNV 199 / 201、INDEL 2,156 / 2,277 = 94.7 %；COLO829T SNV 66 / 67（98.5 %）、INDEL 417 / 430 = 97.0 %；graph-specific loss HG008T 2 SNV + 121 INDEL（9,362 个 INDEL 的 1.3 %），COLO829T 1 SNV + 13 INDEL（0.6 %））。94.7–99.0 % 是这四个比例（HG008T SNV 99.0 %）的范围。这是 d9（HPRC v1.1）上的 absorbed 定义（ALT reads 无 edit 地走 d9 的分支 / skip edge），所以 Response 写 "on HPRC v1.1"；d46 版（E5b）还是 preliminary，graph spelling 在跑。
- 训练区覆盖率的范围：HG008T 94.8–95.9 % 包括 d9 和 d46，COLO829T 89.5–90.4 % 只有 d9（E21 §2 表），Response 已注明。
- 训练区覆盖率：E21 §2 表（SNV 94.8–95.9 % HG008T、89.5–90.4 % COLO829T；COLO829T filtered 2,992–3,364；INDEL 42.6–63.6 %）。只覆盖 d9 臂和 HG008T d46，其余等 collect。
- chr1 覆盖率：`analysis/benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/tables.md` "chr1 tensors per graph"（HG008T SNV truth with a tensor：d9 0.945 / 0.944 / 0.951，d46 0.954 / 0.948 / 0.954，linear 0.980 / 0.978 / 0.976）。detailed doc 的 "95 % vs 98 %" 与此一致。
- Germline：E21 §4 第 3 条（HG008-N / COLO829BL dipcall；HCC1395T 和 CASTLE 用 DeepVariant QUAL ≥ 10，`/scratch/jshen/data/CASTLE/README.txt` 2026-10-12 段）；训练用 germline SNV 数 256,583–381,159（E21 Summary）。detailed doc 的 "0.31–0.42 M germline SNV tensors per dataset" 是 label 总数，以 E21（较新、是训练实际用的数）为准。
- Non 和 determinable：E21 §4 第 4、5 条（R11 `near_truth_allele_mismatch`、R12 `confident_no_truth_allele`；−1 的各种情况）。
- Rule (a)：`analysis/truth_labels.md` 160 行（somatic INDEL 里 HP ≥ 7 的比例 SEQC2 HCC1395 15.5 %、HG008T 74.9 %；CASTLE H1437 / H2009 / HCC1954 20.3 % / 23.8 % / 50.8 %），E21 §4 第 5 条。Response 里因字数没写，放在 R1.m1。
- Build 阈值：E21 §4 第 6 条（AF ≥ 0.08、≥ 3 ALT reads、MAPQ > 10、ALT BQ ≥ 10；label 时不设 AF floor）。
- VAF：E21 §3 表（somatic SNV median 0.53–0.54 HG008T、0.46–0.51 COLO829T；AF < 0.20 占 0.4–1.2 % / 6.5–7.3 %；germline SNV median 0.50–0.60，与 somatic 几乎一样）。
- Purity：dilution 的 Illumina / ONT BQ artefact 见 `analysis/benchmark_7fc7396_20261012/README.md`（mixture 数字是 provisional）。

**还缺什么。**

- [CPU 本地] E21 collect（387193，连同 index array 387137 的后续 task，现为 held）：补 COLO829T d46、linear 臂、dilution、HCC1395T、CASTLE PacBio 的 T3 覆盖率、T5 / T6 VAF 和剩下 24 个 replay check。
- [文字] 把 `truth_labels.md` 和 E21 §4 整理成 Methods / Supplementary Note；用 E21 T3 + 每套的 `truth_recall.json` 做 per-truth 状态的 Supplementary Table，附 VAF 直方图（E21 T6）。
- [用户决定] rule (a)（重复区里没有 truth 的 INDEL 记 −1）现在只用于 SEQC2 / CASTLE；联合训练（E12）里 HG008T / COLO829T 是否统一。
- [GPU RIS] Illumina / ONT dilution 的 BQ artefact 修好后（用户在做），"purity dilutions in joint models" 一句才能不加限定地写。
- [用户决定] E9f（用组装给 off-reference tensor 打标签，P2，可选）做不做。E1.1 §5 量化了不给 off-reference 正样本的代价（没有 rescue 时 0 / 112、1 / 89；rescue 后 32 / 112、46 / 89），这让 E9f 比原来更有价值；不做的话 Response 维持 "excluded from training" + rescue 的说法。
- [CPU 本地] E1.1 §5 的 d46 版（`asm_scores.py` 按每个 d46 run 的 LOW 和训练统计给已分好类的 d46 keys 打分，R3.3 同一条）。这一步是推理，按规定放本地 tequila GPU 跑；出来后把 d46 数字补进证据（E9d 的组装验证在跑，可一起核对）。

### R1.m1 Additional benchmarks

**状态：部分解决**: 新细胞系的数据已经就位（HCC1395T 三平台、CASTLE 三株 PacBio、HCC1954T / H2009T ONT 的 d46 tensor 已 finalize 并 relabel），比较器训练用了哪些细胞系和染色体也已核实（E8b）；HCC1395T chr1 现在有了 Pansoma 的第一批数字（E12a：从没见过 HCC1395 的 HG008T / COLO829T 模型，SNV PASS F1 0.269–0.790，unplaced calls 不计），但这些是 preliminary（旧 normalization、一个 seed、没有 HCC1395T in-sample 参照；COLO829T 模型 "损失在阈值" 靠的是 test-tuned best F1），比较器在 HCC1395T 上的同协议重打分要等用户给已有结果的路径，CASTLE 三株上还没有任何结果。

**Comment（要点）。**

- 在 SEQC2 的 HCC1395 tumor / normal 细胞系上做 benchmark。
- 在 DeepSomatic 研究里另外五个细胞系（CASTLE）上做 benchmark。

**Response.**

We added HCC1395 (SEQC2 truth v1.2.1, 41,072 somatic alleles) and three CASTLE lines, HCC1954, H1437 and H2009 (25,014, 84,744 and 161,091 somatic alleles on chr1–22), each sequenced with Illumina, PacBio HiFi and ONT and processed with the identical pipeline on HPRC v2.1 [待补: the CASTLE Illumina sets and H1437 ONT are still being processed]. On HCC1395, 32,108–36,612 of the 41,072 truth alleles reach a tensor, depending on the platform.

*HCC1395 chr1 with models that never saw HCC1395* (a first, single-seed test; 3,440 somatic SNVs in the benchmark region; PASS F1 after population filtering, threshold chosen on the training tumor, calls with a GRCh38 position): models trained on HG008T reach 0.790 (Illumina), 0.698 (PacBio) and 0.689 (ONT); models trained on COLO829T (Fiber-seq for PacBio) reach 0.269, 0.741 and 0.621. The best points of their chr1 curves, chosen on the test data and thus an upper bound, are 0.810, 0.819 and 0.691, so most of their shortfall comes from the transferred threshold (R1.7). These models standardised base and mapping qualities with the training tumor's statistics; the final models do so per dataset. [待补: the retrained and held-out-line models (E12); an HCC1395 in-sample reference; ClairS-TO and DeepSomatic under the same protocol; the CASTLE lines]

How we use these lines follows from how the comparators were trained.

- DeepSomatic's multi-cancer models were trained on HCC1395, HCC1937, HCC1954, H1437 and H2009, with chr1 held out; our test chromosome is also chr1. The CASTLE truth sets of HCC1954, H1437 and H2009 were derived from DeepSomatic's own calls, and its authors note that this evaluation is biased toward DeepSomatic. ClairS-TO's real-tumor model was fine-tuned on HCC1937, HCC1954, H1437 and H2009 (held-out chromosome chr20 in its documentation, chr1 in its paper). We therefore use these three lines mainly for training diversity and held-out-line tests of Pansoma, and report comparator results on them with this caveat.
- HCC1395 has a truth set independent of both tools (SEQC2). ClairS-TO excluded HCC1395 from training; DeepSomatic held out its chr1 but trained on its other chromosomes, whereas our models above never saw HCC1395. HCC1395 chr1 is therefore our main external comparison, together with HG008T, on which neither tool was trained.
- The SEQC2 and CASTLE truth sets likely omit many INDELs in long homopolymers and STRs (15.5 % of SEQC2 INDELs lie in homopolymers of at least 7 bp, against 74.9 % for HG008T). On these lines, INDEL candidates in such repeats without a truth allele are not treated as negatives, and INDEL results are reported separately, not compared with HG008T.

[待定: HCC1937 and Hs578T, the remaining CASTLE lines, are not included / are added.]

**证据（内部）。**

- DeepSomatic：`experiments/E8b_literature_facts/README.md` D2（chr2–20 train、chr21–22 tune、chr1 held out，"In all the models we describe in this paper"）、D4（五株 multi-cancer；Hs578 完全留出）、D6（HG008 是 DeepSomatic 的 held-out 评估细胞系，全基因组）、D8（我们的 CASTLE truth = `DeepSomatic_HCC1395-model_benchmark`，原文 "strongly biased toward the single-cancer trained models ... biased toward DeepSomatic in general"）、D14（r1.10 文档：tumor-only 模型 chr1 held out）、D3 / D10（HCC1395 SEQC2 41,072 用于 single-cancer 和 tumor-only 模型）。所以 DeepSomatic 在 HCC1395 chr1 上是 "同一细胞系、留出染色体"，Pansoma 的 E12a 数字是 "没见过这个细胞系"。注意 Response 第一条说的是 multi-cancer 模型（D4）；我们比的是 tumor-only 模型：论文里只用 HCC1395 训练、chr1 留出（D10），发布的 WGS_TUMOR_ONLY 在 1.9.0 用全部 manuscript 数据重训（D15），长读长 tumor-only checkpoint 用了哪些细胞系没有写（D16，unverified）。chr1 留出这一点对它们都成立（D2、D14），Response 的结论不变。
- ClairS-TO：E8b C2（SS = GIAB HG002 / HG001 / HG003 / HG004 合成数据，所以没训练过 HG008）、C3（SSRS 用 HCC1937、HCC1954、H1437、H2009 fine-tune；"left out the entire HCC1395 cell line"）、C4（held-out 染色体：README chr20，论文 chr1，未定）。
- E8b §5 第 3 条和 §6：R1.m1 这部分 "Resolved"，剩 ClairS-TO 的 chr20 / chr1 问题和 HCC1937。
- Truth 数：E21 README §5（HCC1395T 41,072、HCC1954T 25,014、H1437T 84,744、H2009T 161,091 on chr1–22）。`/scratch/jshen/data/CASTLE/README.txt` 写的是 VCF record 数（HCC1954 26,578 / H1437 87,410 / H2009 168,397，含其他 contig），和 E21 的 chr1–22 allele 数口径不同。
- HCC1395T 进 tensor 的数：`/scratch/jshen/data/pansoma_v2_tensors/HCC1395T_{Illumina,PacBio,ONT}_d46/tensors/truth_recall.json`（somatic representative 36,612 / 35,749 / 32,108 of 41,072），与 detailed doc 的 "32–37 k" 一致。
- **HCC1395T chr1 cross-sample（E12a，preliminary）**：`experiments/E12a_cross_sample/tables.md`、`summary.json`、README。d46 的 2026-10-12 RIS checkpoint，tequila GPU bf16 预测（array 387379；387323 的 xs task 已取消、在 387379 重跑），阈值 = 模型自己训练肿瘤上的 validation recall-0.9（不在 HCC1395T 上调），PoN 后、HCC1395T BED 内 rtg vcfeval；HCC1395T labels 用 rule (a)。输出 `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/xs__<model>__HCC1395T_<platform>_d46/`。
  - SNV PASS F1（best F1）：HG008T 模型 0.790（0.809）/ 0.698（0.775）/ 0.689（0.699）；COLO829T 模型 0.269（0.810）/ 0.741（0.819）/ 0.621（0.691）（Illumina / PacBio（COLO829T 用 fiberseq 模型）/ ONT）。best F1 是 chr1 曲线上选的点（test-tuned），Response 写明 "chosen on the test data and thus an upper bound"，只用来说明 ranking 能迁移、损失主要在阈值；没有 HCC1395T in-sample 模型，所以原稿 "their loss lies in the transferred threshold" 改成 "most of their shortfall comes from the transferred threshold"（相对它们自己的 best F1，而 COLO829T ONT 的 best 0.691 与 HG008T ONT 的 PASS 0.689 / best 0.699 接近）。
  - Unplaced PASS calls 不在 F1 里（E12a Caveat 3），HCC1395T 上数量不小：HG008T 模型 672 / 126 / 941，COLO829T 模型 444 / 114 / 3,559（Illumina / PacBio / ONT，`tables.md` "unplaced PASS" 列；COLO829T ONT -> HCC1395T 的 3,559 多于它 PoN 后的 2,261 个 placed records）。Response 加了 "calls with a GRCh38 position"。
  - COLO829T Illumina -> HCC1395T：阈值 0.640，PoN 后 TP / FP / FN 540 / 38 / 2900，P / R 0.934 / 0.157。
  - 3,440 个 somatic SNV：每行 PoN TP + FN = 3,440（如 540 + 2900、2677 + 763），`xs__COLO829_Illumina_SNV_d46_7fc7396__HCC1395T_Illumina_d46/vcf_chr1/eval_bed/report.txt` 写 "PoN tags 66 of 3,440 truth records"。
  - Illumina INDEL -> HCC1395T：HG008T 模型 0.174（0.247），COLO829T 模型 0.307（0.331）；HCC1395T chr1 INDEL truth 133（35 + 98）。E12a Caveat 4：SEQC2 INDEL truth 缺 homopolymer INDEL，"not comparable across targets"；同一个 target 上两个模型的 PoN ceiling 也不同（COLO829T 模型 0.429、HG008T 模型 0.789；两者写进 VCF 的 raw record 数 538 vs 1,779，`report.txt`；原因没查），所以 Response 没给 INDEL 数字。
  - Preliminary 的原因：E12a Caveats 1、2、5（旧 checkpoint，BQ / MAPQ 用训练肿瘤的统计 z-score；一个 seed；"No in-sample HCC1395T model yet"）。三个 d46 run 的本地 in-sample 重预测复现了 RIS（PoN F1 0.603 / 0.575 / 0.736 vs 0.605 / 0.575 / 0.737；HG008T PacBio 是 CPU fp32（387135），另两个是 tequila GPU bf16（387379 task 1–2），见 R1.7 证据）。
- CASTLE tensor（同一类文件）：HCC1954T PacBio representative 22,876 / 25,014，H1437T PacBio 59,359 / 84,744，H2009T PacBio 148,873 / 161,091，HCC1954T ONT 22,806，H2009T ONT 150,071（ONT 两个是 QUAL10 relabel 后的 truth_recall.json，2026-10-13 12:27 / 12:55）。H1437T PacBio somatic SNV 只有 70.9 %（56,703 / 79,974）有 tensor，异常未查清（detailed doc 风险 c）。
- 进度（squeue 2026-10-13 15:39）：HCC1954T / H2009T ONT 的 QUAL10 relabel 已完成；H1437T ONT finalize（385470）在跑；H2009T Illumina run（385782）在跑；HCC1954T / H1437T Illumina 在 gamsort（385759、385769）；H1437T ONT 和三株 Illumina 的 relabel（387131–387134）pending。
- HP ≥ 7 比例：`analysis/truth_labels.md` 160 行（SEQC2 HCC1395 15.5 %、HG008T 74.9 %，来自 `analysis/indel_repeat_context_predictions_20261006`；CASTLE H1437 / H2009 / HCC1954 20.3 % / 23.8 % / 50.8 %）；CASTLE INDEL truth 要 Illumina 和至少一个长读长 DeepSomatic 模型都 call 到（同一行）。Response 原稿写 "omit most INDELs"；15.5 % vs 74.9 % 是不同肿瘤之间的间接比较（`truth_labels.md` 原话是 "漏掉很多"），所以改成 "likely omit many"。
- 比较器：`experiments/E8_comparators/README.md`（PENDING，没有 call set）；用户 2026-10-13 决定不重跑，只把已有结果按同一协议重打分（路径待给）。
- HCC1937 / Hs578T：本地 `CASTLE/benchmarks/` 有 HCC1937 的 truth 和 BED，没有 reads（detailed doc E12b：约 0.6 TB 下载）；Hs578 的 truth 和 reads 本地都没有。

**还缺什么。**

- [用户决定] 给出已有 ClairS-TO / DeepSomatic 结果的路径；CASTLE 三株上 ClairS-TO 用的是 ssrs 还是 ss 模型（chr20 / chr1 问题，E8b C4）。这是这一条最主要的缺口。
- [CPU 本地] 拿到路径后，用 E8 的 `score_comparator.py` 在 HCC1395T chr1 上按同一协议（BED、PoN、rtg 选项）重打分，和 E12a 的 Pansoma 数字放进同一张表。
- [GPU RIS] 用 `--norm dataset` 重训的 d46 单细胞系模型和 E12 的 held-out-line 模型（包含 HCC1395T / CASTLE），最终数字以它们为准；之后在 tequila GPU 上重跑 -> HCC1395T 的预测。
- [GPU RIS] HCC1395T 三平台的 in-sample 模型（chr2–22 训练），作为 HCC1395T chr1 的参照（E12a Caveat 5）；它和 DeepSomatic 的 "同一细胞系、留出 chr1" 设定可比。
- [CPU 本地] CASTLE ONT / Illumina 集合 finalize 和 QUAL10 relabel（H1437T ONT、三株 Illumina 的 pipeline job 在队列里）；之后对 CASTLE 三株做同样的 cross-sample 预测（推理在 tequila GPU）。
- [用户决定] H1437 PacBio 异常查清前，H1437 是否进测试；HCC1937（下载）和 Hs578T 加不加，决定后改掉 Response 最后一句。
- [文字] Methods 写清 SEQC2 / CASTLE truth 的来历和 rule (a)；INDEL 结果单列。

### R1.m2 Tensor windows truncated at node boundaries

**状态：已解决**: 跨节点的 tensor 已经实现（v4 README §8）；E1.1 的 ch6 counterfactual 已在全部 18 个 SNV 模型上跑完，12 个图模型都显示 classifier 在用 path-aware 的 ch6（linear 模型作对照，没有变化）；剩下的只是示意图和措辞。显式的 bubble 建模（E24）没做，Response 按 limitation 写，做不做是可选的用户决定。

**Comment（要点）。**

- 窗口限制在当前节点内、在节点边界截断，相邻节点和跨路径的证据丢掉了。
- 短节点、单碱基节点很难处理；很多 alternative 节点比窗口还短。
- 可信的 graph-native 方法需要跨相邻节点的 path-aware tensor，并显式建模 edges、bubbles、alternative paths 和 competing path support。

**Response.**

This described the first version; the revised tensors are not truncated at node boundaries.

- Each tensor (8 channels × 200 rows × 101 columns) is centred on the candidate site, and each row follows its read's
  alignment walk across node boundaries, on reference and alternative nodes alike. A short or single-base node
  therefore keeps the flanking context of its neighbouring nodes.
- The graph base (channel 5) and the number of haplotype paths through the node (channel 6, relative to the graph's
  haplotype count) are filled per column from the node the read traverses there.
- Rows are grouped by allele (A1, A2, …, REF, OTHER) and clustered within a group by their events, which include every
  node transition, so reads taking the same alternative path sit together.

The classifier uses this path information.

- Alone, the site column's path count separates somatic from germline chr1 SNV test tensors with an AUC of 0.837–0.892
  (three HG008T sets), and the number of distinct codes within ±10 columns, a proxy for node changes, with 0.805–0.831.
- We edited channel 6 at inference, without retraining, in all 12 graph-based SNV models (two HPRC graphs, HG008T and
  COLO829T, three platforms; chr1 test tensors). With the plane set to its training mean, per-tensor somatic average
  precision falls from 0.199–0.352 to 0.044–0.066 (HG008T) and from 0.639–0.829 to 0.329–0.520 (COLO829T). Setting
  every node to "all haplotypes" has a similar effect and raises the number of germline tensors called somatic 3.4- to
  7.6-fold. On the GRCh38-only graph, where the plane is constant, nothing changes.
- After population filtering (a per-tensor estimate; off-reference tensors left out) these two edits change the
  truth-level F1 less (−0.066 to +0.007, except HG008T ONT on HPRC v1.1 with the training mean, 0.596 to 0.145);
  another tensor's plane lowers it to 0.212–0.551. Whether the plane improves the final calls is the retraining
  ablation of R1.1.
- The plane also suppresses somatic SNVs on nodes outside GRCh38: the HPRC v1.1 PacBio and ONT models score 0 of 112 and
  1 of 89 such HG008T somatic SNVs (assembly-derived) above the threshold, against 61 and 49 with every node set to all
  haplotypes. Off-reference SNVs therefore have a separate calling rule (R3.m1, R3.4), which recovers 32 and 46.

A schematic of a tensor crossing node boundaries is added as Supplementary Fig. [待补].

We agree that the tensors do not model bubbles explicitly: competing path support is represented through the reads'
walks and the per-node haplotype counts, not as a bubble-level object. We state this as a limitation (R1.1).
[可选, 若做 E24: a locus mode adds the reads of competing branches with their haplotype frequency.]

**证据（内部）。**

- 实现：`indexed_gam_pipeline_v4/README.md`
  - §8 "Site layout"：site 从第 50 列开始（101 列居中）；shape `(8, 200, 101)`（README 开头）。
  - §8 "Rows"：blocks `A1, A2, ..., REF, OTHER`；block 内按 average-linkage 排序，events 包括 "every node transition of its
    visible path"。
  - §8 Channels 表：ch5 = "the graph reference base under that row/column"；ch6 = "distinct GBWT paths of the column's
    node, relative to the graph's haplotype count H"。
  - §4 `--chromosomes`："Reads are not filtered: a read on a chr1 target keeps its columns on neighbouring nodes of any
    block"。
- E1.1 全部跑完：`experiments/E1.1_ch6_counterfactual/results/summary.json` 的 `done` 有 18 个 run，`pending` 为空；表在
  `results/tables.md`（`collect.py` 生成）。
  - 12 个图模型：HG008T d9 × 3、HG008T d46 × 3、COLO829T d46 × 3、COLO829T d9 × 3。COLO829T 的 d9 模型是 purity-joint
    模型（P100+P75+P50+P25 联合训练），在 P100 set 上测，PASS 阈值是四个 purity 合起来定的（README "Models and data"）。
  - ch6 单独（§6，无模型）：site code AUC 0.837（HG008_PacBio，中位 100 / 96）、0.890（HG008_PacBio_d46，100 / 79）、
    0.892（HG008_ONT_d46，100 / 79）；distinct codes ±10 AUC 0.805 / 0.829 / 0.831，中位都是 1 / 2。
  - somatic AP（§4），base → train_mean：HG008T d9 0.199 / 0.233 / 0.232 → 0.060 / 0.053 / 0.044；HG008T d46 0.312 /
    0.352 / 0.326 → 0.053 / 0.066 / 0.051；COLO829T d46 0.820 / 0.829 / 0.755 → 0.505 / 0.520 / 0.421；COLO829T d9 joint
    0.693 / 0.714 / 0.639 → 0.482 / 0.512 / 0.329。all_H：HG008T 0.061–0.087，COLO829T 0.347–0.516。shuffle 更低：HG008T
    0.028–0.043，COLO829T 0.192–0.308。
  - germline ≥ t 相对 base 的倍数（all_H，§3）：HG008T d9 3.58 / 3.50 / 3.79，HG008T d46 6.75 / 6.83 / 6.64，COLO829T
    d46 7.57 / 6.72 / 7.39，COLO829T d9 joint 3.45 / 3.41 / 4.19 → 3.41–7.57。
  - PoN 后 truth-level F1（§4 最后一张表），base → all_H / train_mean：变化 −0.066（COLO829_ONT_joint train_mean 0.706 →
    0.640）到 +0.007（COLO829_Illumina_d46 all_H 0.784 → 0.791；COLO829_Illumina_joint all_H 0.675 → 0.682）。例外：
    HG008_ONT（d9）train_mean 0.596 → 0.145（somatic ≥ t 只剩 65 / 663，germline ≥ t 也降到 1,365，0.27×）。shuffle
    PoN 后 F1 0.212–0.421（HG008T）、0.274–0.551（COLO829T）→ Response 写 0.212–0.551。
  - 这个 "PoN 后 F1" 是 E1.1 自己的 per-tensor 估计，不是 rtg F1：只算 on-reference tensor（off-reference / unplaced 不算），
    README "Measures per mode" 明说 "the after-PoN F1 here is a per-tensor estimate for comparing modes, not the rtg F1"；
    base 模式下 TP 和 rtg 差 ≤ 12，FP 在 11 / 12 个图模型里差 ≤ 52（`COLO829_ONT_joint` 399 / 595 例外）。所以 Response
    写成 "a per-tensor estimate; off-reference tensors left out"，并且只说这两种 edit 的变化小，shuffle 的大降也写上。
  - linear 对照：§1 training ch6 mean 100.0、cells at 100 = 1.00（H = 17，每个 cell 都是 100）；§3 六个 linear run 的
    所有模式数字完全相同。
  - 组装 off-reference somatic SNV（§5，HG008T，全部常染色体，都不是训练标签）：HG008_PacBio（d9）base 0 / 112 ≥ t，all_H
    61 / 112；HG008_ONT（d9）base 1 / 89，all_H 49 / 89。site10_high（= `predict` 的 off-reference rescue 变换）下 argmax
    somatic 32 / 112 和 46 / 89（≥ t 45 / 45）；Response 的 "recovers 32 and 46" 用的是 argmax（`graph_vcf` 对
    off-reference SNV 用 `offref_call` = argmax，见 E1.1 README §5 和 `pansoma_net_v2/README.md` "Off-reference rescue"）。
    key 来自 `pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv`（HG008-N 组装窗口经 giraffe 落在 d9
    off-reference 节点上的 truth 事件）。d46 的组装 key 已经分好类（`asm_d46_classify.py`，job 387185：345 个 d46
    off-reference somatic 事件，154 SNP），但还没有模型打分（`results/assembly_keys.json` 只有两个 d9 run；README
    Remaining："Optional, user's decision ... a short GPU inference on tequila"）。
  - §2（README "How the runs were done" 第 3 条更正了 tables.md 的表头）：只有 4 个 run 的 base 模式是 CPU fp32
    （`HG008_PacBio`、`HG008_ONT`、`HG008_PacBio_d46`、`HG008_ONT_d46`，AP 0.2333 / 0.2331、0.2316 / 0.2323、0.3523 /
    0.3528、0.3264 / 0.3283），其余 14 个是 tequila GPU bf16；和 RIS GPU bf16 `predict` 比，18 个 run 的 AP 差都 ≤ 0.0047
    （最大 COLO829_fiberseq_linear 0.2832 / 0.2879）。所以这里的 counterfactual 和 RIS 的预测可比。（原稿写 "CPU fp32 和
    RIS GPU bf16"，按 README 改。）
- 和旧稿不一致：详版 / 简版写的是 "true somatic tensors have a median of 0 transitions within ±10 columns, germline 2
  (AUC 0.72)"。`pansoma_net_v2_runs/pathcount_kind_auc_20260930/` 里只有 `pc_auc.py`，没有保存的输出，数字核不了（而且是
  10-01 ch6 重建之前的 tensor）；改用 E1.1 §6 的数字（定义不同：distinct codes 是 node change 的代理，path count 相同的相邻
  节点数不到）。
- 示意图脚本已存在：`scripts/visualize_tensor.py`。

**还缺什么。**

- [文字] Supplementary Fig.：一个跨节点 tensor 的示意图（`scripts/visualize_tensor.py`，本地出图，挑一个单碱基 branch
  节点上的 site）。
- [GPU RIS] E3 的无 ch6 臂（`--drop-planes path_count`，d9 / d46 重训）：回答 "path 信息是否带来收益"，属于 R1.1；本条
  Response 只引用，不依赖它的数字。
- [用户决定]（可选，主要给 R1.1 / R3.3）E1.1 §5 补 d46 的组装 off-reference key 的模型打分：key 已有
  （`asm_d46/offref_somatic_d46.tsv`），只差 `asm_scores.py` 在三个 d46 run 上的推理（tequila GPU，短；E1.1 README
  Remaining 列为用户决定，没提交）。
- [用户决定] E24（locus / bubble-aware tensor，改 pipeline，2–3 周）做不做；不做就保留 Response 最后一段的 limitation。

---

## 6.2 Reviewer 2

### R2.1 Why include the HapMap benchmark

**状态：部分解决**: 为什么纳入、有哪些局限，这部分文字已经有 E8b 核实过的事实撑着，"分类器靠 graph prior（ch6）压 germline allele" 现在也有 E1.1 的 ch6 counterfactual 直接支持（过 PoN 以后影响小得多；图只保留约 10 % 以上 haplotype 的 node，所以只吸收常见的 positives），替代用的 COLO829T dilution 结果也有了；但 HapMap 最后去掉还是放进 Supplement、原稿用的是哪张图，要用户定，Illumina / ONT dilution 带 read-source BQ artefact、分数偏乐观（同 R1.5）。

**Comment（要点）。**

- HapMap 的 positives 是 germline variant，它们的 donor 在 HPRC 图里，所以理论上 pangenome 比对会把它们去掉。
- 请讨论为什么纳入这个 benchmark，以及它在模拟 somatic mutation 上的局限。

**Response.**

We included the SMaHT HAPMAP6 mixture because it offers a large truth set at very low allele fractions (0.25–16.5%) [待核实（用户）: the rationale given in the original manuscript]. We agree that it has three limits for a pangenome-based tumor-only caller:

1. Its positives are germline alleles of the five minor donors that the HG005 background lacks. They are not acquired mutations, so they carry no somatic mutational process, clonal structure or tumor copy-number context.
2. Four of the donors are in HPRC v1.1, and all six are in HPRC v2.1. Our graphs keep only nodes carried by about 10% or more of their haplotypes, so they absorb the positives that are common among the graph's haplotypes (a donor's private alleles are filtered out even when the donor is in the graph). The trained classifier relies on this prior. On our tumor benchmarks, when the haplotype path-count channel is overwritten at inference so that every read appears to pass through nodes shared by all haplotypes, 3.41 to 7.57 times as many germline test tensors are called somatic (12 graph-based SNV models, chr1). After population filtering the effect on F1 is much smaller (R1.1).
3. The donors are public reference samples. Their rare alleles are found in CoLoRSdb, one of our population filters, far more often than those of our cancer benchmarks' normals (chr20 SNVs: 14.5–25.9% against 4.0–7.9%), so population filtering also works against these positives.

The mixture therefore measures sensitivity to low-fraction germline alleles. It does not measure tumor-only somatic calling.

We removed it from the main benchmark. In its place we use tumor/normal dilutions that contain true somatic mutations: COLO829T at 75, 50 and 25% purity on three platforms, with one model per platform trained across purities (chr2–22) and tested per purity on chr1 (R1.5, R2.4). [待补: Illumina and ONT after removing a read-source base-quality cue (R1.5); a second tumor held out from training (SEQC2 HCC1395).]

[用户决定] Either HapMap no longer appears, or it appears only in the Supplement, computed on a graph without any of the six donors (40 HPRC individuals + GRCh38 + CHM13) and with these limits stated. [待核实（用户）: if the submitted HapMap results already used this donor-free graph, we say so.]

**证据（内部）。**

- 混合比例、五个 minor donors、VAF 0.25–16.5 %：`experiments/E8b_literature_facts/README.md` M1。truth 是 Kong et al. 的 graph truth v1.6（M2）。
- Donor 在图里：E8b H8，v1.1 有 4/6，v2.1 有 6/6。E8b §5 第 1 条更正了详细版 R2.1 的 "four of its six donors are in both HPRC graphs"，这里用的是更正后的说法。
- CoLoRSdb：E8b §4.4（chr20 SNV）。rare matched %：HapMap 六个 donor（含 HG005）14.5–25.9 %，HG008-N、COLO829BL 和四个 CASTLE normal 4.0–7.9 %。核对时把原稿的 "two to five times"（P11 的约数，按 §4.4 两端相除并不都落在 2–5 倍里）改成直接写这两组百分比。原文标为 [inference]，所以 Response 写 "found more often"，不写 "contains"。
- 图只保留常见 allele（核对时加进 Response 第 2 条）：E8b H3（v1.1 d9 = `--filter 9`，去掉少于 9 / 90 个 haplotype 覆盖的 node，reference path 不过滤）和 H6（v2.1 d46 = `--filter 46`，46 / 464，约 10 %）。所以 donor 的 private allele（1–2 个 haplotype）即使 donor 在图里也会被滤掉；图能吸收的是在 haplotype 里常见的 positives。原稿 "The graph is built to absorb exactly such alleles" 比证据强，已改。HapMap positives 里有多少落在图的非 reference node 上还没算（见还缺什么）。
- Donor-free 图：
  - 出处是 E8b H10；
  - 我用 `bcftools query -l mgi_non_HG005.d9.raw.vcf.gz` 核过：40 个个体加 CHM13，六个 donor 和 NA19240 都不在里面（核对时又直接读了 VCF header：41 个 sample，0 个 donor）；
  - 原稿可能就是在这张图上跑的，线索见 R1.5 证据：`tmp/hapmap_nonHG005_20261012/submit_prep.sh` 第 2–4 行和 `/scratch/jshen/Scripts/HapMap_train_data_preparation.sh` 第 22–33 行。
- 模型靠 graph prior 压 germline（Response 第 2 条）：`experiments/E1.1_ch6_counterfactual/results/tables.md`。
  - §4 "Germline tensors with p_somatic >= t"：all_H 模式（推理时把 ch6 每个有效 cell 改成 100，即 "all H haplotypes run through this node"，不重训）相对 base 的倍数：HG008T d9 3.58 / 3.50 / 3.79；HG008T d46 6.75 / 6.83 / 6.64；COLO829T d46 7.57 / 6.72 / 7.39；COLO829T d9 purity joint（在 P100 上测）3.45 / 3.41 / 4.19（Illumina / PacBio-fiberseq / ONT）。范围 3.41–7.57。
  - §6：只用 ch6 的 site code 区分 somatic 和 germline test tensor，AUC 0.837（HG008T PacBio d9）、0.890（PacBio d46）、0.892（ONT d46）；site code 中位数 somatic 100，germline 96（d9）/ 79（d46）。
  - linear 模型（GRCh38-only，H = 17）的 ch6 每个 cell 都是 100（§1 "cells at 100" 1.00），all_H 和 base 完全一样（§3，例如 HG008T PacBio linear 的 germline >= t 都是 127,888），没有这个效应。
  - 过 per-tensor PoN 以后 F1 的变化小得多（§4 "Truth-level F1 at t after the per-tensor PoN"，base → all_H：HG008T Illumina d46 0.560 → 0.511，COLO829T Illumina d46 0.784 → 0.791），因为 PoN 也去掉大部分这类 germline。这个 F1 是 per-tensor 的估计（只算 on-reference tensor），不是 rtg 的 F1（E1.1 README "Method"）。核对时在 Response 第 2 条加了 "After population filtering the effect on F1 is much smaller"，免得读成最终 call set 也差 3–7 倍；细节在 R1.1。
  - 12 个模型都是 in-sample（同一肿瘤 chr2–22 训练、chr1 测试），阈值是各自的 PASS（validation recall-0.9）。
  - 这些都是 cancer benchmark（HG008T、COLO829T），不是 HapMap。"donor-free 图上 HapMap positives 不再被 ch6 压" 是推论 [inference]，Response 里没写。
- 替代用的 dilution 结果：`analysis/benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012/README.md`、`tables.md`（具体数字见 R1.5）。

**还缺什么。**

- [用户决定] HapMap 是直接去掉，还是只在 donor-free 图上放进 Supplement（E23；重跑 chain 387101–387105 在排队）。这一项和 R1.5 一起定。
- [用户决定] 原稿纳入 HapMap 的理由和当时用的图（`mgi_non_HG005.d9` 还是 `hprc-v1.1 d9`），用来定 Response 第一句和最后一句的措辞。
- [GPU RIS] Illumina / ONT dilution 在修掉 read-source BQ cue 以后重训（E14b，同 R1.5）；在这之前它们的分数偏乐观，Response 不引用 dilution 的数字（只指向 R1.5），并在 [待补] 里写明。
- [用户决定] 第二个肿瘤的 dilution 用 E14c（in-silico）还是 E14d（SEQC2 湿实验 titration）。
- [CPU 本地]（可选）数一下 HapMap positives 里有多少落在 d9 / d46 / `mgi_non_HG005.d9` 的非 reference node 上（例如 chr20），量化 "图只吸收常见的 positives"、donor 进不进图差多少。没有这个数之前，Response 第 2 条只写定性的说法。
- [文字] Discussion 里加一段 HapMap 局限，三条和 Response 一致。
- [文字] Response 第 2 条的 E1.1 句子（3.41–7.57 倍）要和 R1.1 / R1.m2 里 E1.1 的写法用同一组数字、同一种 mode 描述（all_H）。

### R2.2 How the truth sets were generated

**状态：部分解决**: 每个 truth 的构建方式已按原始文献 / README 核实（E8b G1–G4、S1–S5、Q1–Q3、D8、M1–M2），"truth 里有 germline" 的说法也有组装证据，SEQC2 / CASTLE 的 repeat INDEL 怎么进标签（rule (a)）也补上了；但 HG008 用 v0.2（未公开的 draft）还是换 v0.3 要用户定，v0.2 自己的构建说明和我们的 COLO829 文件是否等于公开版还没核实，HapMap 那一条要看 R1.5 的去留。

**Comment（要点）。**

- 简要说明每个 benchmark 的 truth 是怎么生成的。
- 特别是要支撑文中 "truth set 可能含 germline variants 和测序 artefacts" 的说法。

**Response.**

We added a Supplementary Note on each truth set. Calls are scored with rtg vcfeval inside each benchmark's regions;
training negatives come from regions confident in both the somatic and the germline truth. Germline truth is the
dipcall of the matched-normal assembly (HG008-N, COLO829BL) or, for the other four lines, DeepSomatic's DeepVariant
calls of the normal (QUAL ≥ 10).

- *HG008T:* GIAB HG008-T somatic small-variant draft benchmark v0.2 (18,952 records: 9,609 SNVs, 9,343 INDELs), a
  pre-release draft provided by GIAB [待定: or the public v0.3 of 2026-04-25]. It is derived from comparing the HG008-T
  tumor assembly with the HG008-N normal assembly (v6.2). [待核实: whether v0.2, like v0.3, curated candidates with
  PacBio HiFi, ONT and Illumina/Element reads.]
- *COLO829T:* the SMaHT COLO829 reference set (Ha et al., 2025): candidates from four Illumina callers (Mutect2, Strelka2,
  RUFUS, VarNet; COLO829 200× vs COLO829BL 300×), validated with PacBio reads (ALT present in the tumor, absent in the
  normal): 44,005 SNVs, 961 insertions, 1,098 deletions. It is scored genome-wide: the SMaHT easy, difficult and
  extreme regions together tile the genome.
- *HCC1395:* SEQC2 v1.2.1 high- and medium-confidence somatic SNVs and INDELs (41,072 records) in the SEQC2
  high-confidence regions v1.2, built from 378 call sets (21 replicates × 3 aligners × 6 callers) scored by SomaticSeq and
  NeuSomatic; it covers only regions analysable with Illumina data.
- *HCC1954, H1437, H2009:* the CASTLE somatic truth and high-confidence regions released with DeepSomatic (26,578 /
  87,410 / 168,397 records). Made with DeepSomatic models, they favour DeepSomatic, as its authors note; we say so
  wherever these lines are compared with it.
- [若保留, R1.5: *HapMap mixture:* SMaHT HAPMAP6; positives are minor donors' germline alleles absent from HG005, expected
  VAF 0.25–16.5%, from the pangenome-graph truth of Kong et al. (v1.6).]

For HCC1395 and the CASTLE lines, unmatched INDEL candidates in homopolymers of ≥ 7 bp or short tandem repeats are not
used as negatives: the short-read SEQC2 truth has few such INDELs (15.5% of its INDELs in homopolymers of ≥ 7 bp,
against 74.9% for HG008T), and the CASTLE lines follow the same rule.

We removed the claim about sequencing artefacts (no direct evidence) and state the germline claim narrowly:

- 209 HG008T somatic truth alleles are identical to an HG008-N germline allele.
- Of the 110 graph-absorbed HG008T truth alleles carried by the normal, 109 are the germline allele of the other normal
  haplotype; flanking phased SNPs indicate local changes, not regional LOH: somatic by GIAB's definition, but
  germline-like to a tumor-only caller.
- At 55 COLO829T truth INDELs a COLO829BL haplotype already carries the truth allele over the whole repeat; 33 of them
  are absorbed by the graph, 32 of these in SMaHT extreme regions.

We describe these as benchmark ambiguities at a small number of sites, not as errors of the truth sets, and no longer
use them to explain the lower recall (R2.3, R3.7).

**证据（内部）。**

- HG008T：`experiments/E8b_literature_facts/README.md` §3.5
  - G1：v0.2 `_tumorvariants.vcf.gz` 18,952 records（9,609 SNV、9,343 INDEL）；`_all.bed` 41,835 regions、2,809,425,965 bp。
  - G2：v0.3（2026-04-25）是 "initial public draft"，18,946 variants；v0.2 不在公开 FTP。
  - G3：v0.3 构建方法（tumor v3.2 vs normal v6.3 组装比较；PacBio HiFi、ONT、Illumina/Element reads）。
  - v0.2 本身：VCF header 有 `INFO/HG008Nv62SOMATICVARIANT`（"somatic variant calls on the HG008N v6.2 assembly"）和
    `vcfeval ... dipHG8NvGRCh38 ... dipHG8Tchr19fixvGRCh38`，说明 v0.2 也是 tumor / normal 组装比较、normal 用 v6.2
    （`/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/`，本次重查 header 确认）。header 里 `bcftools annotate` 的路径含
    `.../repun/HG008/`，提示 v0.2 也走了 repun，但 header 证明不了 reads 的 curation（PacBio HiFi / ONT / Illumina-Element
    是 v0.3 README 的说法）。所以 Response 里 "with candidates checked against ... reads" 改成了 [待核实]，"curated" 和
    "primarily" 也去掉（只写 header 能证明的部分）。
- 用法：`labels.manifest.json`（HG008T_PacBio、COLO829T_ONT、HCC1395T_Illumina_d46）的 `truth.somatic.bed` 就是 `ris/vcf.sh`
  的 vcfeval `--bed`；label 0 的区域 = somatic BED ∩ germline BED（`analysis/truth_labels.md` §0 "confident 区"）。
- COLO829T：E8b §3.6 S2（四个 caller + PacBio 验证）、S3（44,005 / 961 / 1,098 = 我们的 44,005 SNV + 2,059 INDEL）、S4（easy
  / difficult / extreme v2）；`pansoma_v2_tensors/COLO829T_truth/README.txt`："the three are disjoint and tile chr1-22,X,Y,
  so the union is every base"。也就是说 COLO829T 的 "benchmark BED" 等于全基因组：旧稿写 "over the union of SMaHT easy,
  difficult and extreme regions" 会让人以为有一个受限区域，已改成 "scored genome-wide"。COLO829T 的 label 0 区域 =
  COLO829BL dip.bed（同一 README）。
- HCC1395：E8b §3.7 Q1–Q3；`/scratch/jshen/data/CASTLE/README.txt`（41,072 records，BED 2,485,467,324 bp；germline =
  DeepSomatic 的 HCC1395BL DeepVariant hybrid，QUAL ≥ 10）。
- CASTLE：`/scratch/jshen/data/CASTLE/README.txt`（HCC1954 26,578 / H1437 87,410 / H2009 168,397，somatic_only PASS）；E8b D8
  （DeepSomatic 原文 "strongly biased toward the single-cancer trained models ... and is biased toward DeepSomatic in general"）。
- HapMap：E8b M1–M2；`pansoma_v2_tensors/HapMap_truth/README.txt`（v1.6，7,870,692 records）。
- 新加的一段 "For HCC1395 and the CASTLE lines ..."（repeat INDEL 不当负样本，即 rule (a)）：
  - `/scratch/jshen/data/CASTLE/README.txt` "Labels: rule (a) ... an INDEL tensor in a GRCh38 homopolymer >= 7 bp or STR
    (unit 2-6, >= 3 units) without a truth allele is -1 instead of 0"；`analysis/truth_labels.md` 常量表
    `IGNORED_REPEAT_CLASSES`（`HP>=7`、`STR`）和规则 R13；`--ignore-repeat-indels` 那一行写明 `HCC1395T_*`、`HCC1954T_*`、
    `H1437T_*`、`H2009T_*` 都用了，HG008T / COLO829T 的六个数据集没用。
  - 理由：`experiments/E12a_cross_sample/README.md` Caveats 第 4 条 "HCC1395T truth (SEQC2) is short-read derived. Its
    INDEL truth misses many homopolymer INDELs, which is why the labels use rule (a)"；同一条还给了 HCC1395T 的 PoN-limited
    INDEL ceiling 0.43–0.79（HG008T / COLO829T 0.25–0.27），所以各 truth 之间的 INDEL 数字不可比，Supplementary Note 要说。
  - Response 新加的 15.5 % / 74.9 %：`analysis/indel_repeat_context_predictions_20261006/README.md`（SEQC2 HCC1395 v1.2.1
    sINDEL 1,625 个，HP ≥ 7 占 15.5 %；HG008T v0.2 chr1-22 8,496 个，74.9 %）；`analysis/truth_labels.md` R13 段引用同样的数。
    原稿写 "misses many such INDELs"，这是推断（没有独立的漏检测量），改成只陈述比例（"has few such INDELs"）。
  - CASTLE 为什么也用：CASTLE README 只写 "for these truth sets"（用户 2026-10-06 的决定）；`truth_labels.md` R13 段给了
    CASTLE truth 的 HP ≥ 7 比例 H1437 / H2009 / HCC1954 20.3 % / 23.8 % / 50.8 %，并说 CASTLE INDEL 要 Illumina 和至少一个
    长读长 DeepSomatic 模型都 call 到（这一条 E8b 没核过原文）。HCC1954 的 50.8 % 不算低，所以 Response 仍写成 "the CASTLE
    lines follow the same rule"，不给 CASTLE 单独的理由。
- HCC1395T 现在有 chr1 结果了（E12a 跨样本，preliminary）：HG008T 的 d46 模型 SNV PASS F1 0.790 / 0.698 / 0.689
  （Illumina / PacBio / ONT，`experiments/E12a_cross_sample/tables.md`）。这些数字属于 R1.7 / R1.m1，这里只说明 HCC1395 的
  truth 已经在评估里用上了。
- 209：`analysis/truth_labels.md` §8.1 第 5 点。
- 110 / 109 和 flanking SNP（76 both retained、29 other hap lost、0 event hap lost）：
  `analysis/graph_absorbed_somatic_20261001/README.md`（Short answer 和 §3 "Present loci"）。
- 55 / 33 / 32：`analysis/graph_absorbed_somatic_colo829t_20261005/README.md`。Short answer "33 (60.0%) of the 55 truth INDELs
  whose allele the normal carries are absorbed"，"32 of them are in SMaHT Extreme regions"；patient-frame 表 "normal carries
  ALT | 55 | 33 | 60.0%"（分母是 chr1-22 全部 1,912 个 truth INDEL）；定义 "a COLO829BL haplotype already has GRCh38 + truth
  over the whole array"。旧稿只写了被图吸收的 33 个；和 HG008T 的 209（全部 truth）对应的应该是 55，已改。
- 和旧稿不一致：
  - 详版写 HG008T "derived from multi-technology tumor data against the curated diploid HG008-N assembly"；E8b G3 是
    "derived primarily from a tumor-normal assembly comparison"，已改。
  - `truth_labels.md` 把 209 解释为 "可能是肿瘤 LOH 或两个 benchmark 冲突"；但在 110 个 absorbed present 位点上 flanking SNP
    排除了 regional LOH，所以 Response 不写 LOH，209 只说 "identical"，不说原因。
  - SMaHT GitHub README 说 "three variant callers"，论文说四个（E8b S2）；用论文。

**还缺什么。**

- [用户决定] HG008 benchmark：v0.2 写成 "pre-release draft provided by GIAB"，还是换 v0.3（要下载约 0.3 MB，所有 HG008T set
  relabel、重新 vcfeval；counts 几乎一样）（E8b 第 5 节第 7 条）。
- [文字] 核实 v0.2 的构建方法是否与 v0.3 README 描述一致（问 GIAB，或找 v0.2 附带说明）；不一致就只写 header 能证明的部分。
- [用户决定] 是否下载公开的 `Truthset_{SNV,Indel}_COLO829.vcf.gz`（parklab repo）核对是否等于 qfu 的 `Step1_Validated_union_*`
  文件（E8b S5，未核实）。
- [用户决定] HapMap 去留（R1.5）：决定 HapMap 那一条是否保留。
- [用户决定] rule (a) 是否也用到 HG008T 和 COLO829T（详版第 7 节第 7 条）；用的话，新加那一段的 "For HCC1395 and the CASTLE lines" 要改成
  全部数据集，并且这六个数据集要 relabel。
- [文字] Supplementary Note：每个 truth 一段（来源、构建、版本、区域、在 label 和评估里怎么用），引用 E8b 里的原文。

### R2.3 Explaining the lower recall; normal sample at missed sites

**状态：部分解决**: 漏检的 read-level 分类和 matched-normal 组装证据已经完成，可以直接回答后两问；但和 DeepSomatic 的 recall 差距要等 user 已有的比较器结果按同一协议 re-score 才能量化，旧的两个 FN 类别的原始定义也还要核实。

**Comment（要点）。**

- recall 比 DeepSomatic 低约 20 %（COLO829T）和 10 %（HG008T），用 truth 错误来解释需要展开。
- 解释 Supp. Fig. 2 的 false-negative 类别 "variant node not ranked #1" 和 "first aligned node mismatch"。
- 对落在已有 graph node 上的 false negative，检查 normal 样本的 reads 是否优先比对到 alternative path。

**Response.**

*Recall difference to DeepSomatic.* We withdrew the truth-error explanation: the assembly evidence below shows that most graph-absorbed truth alleles are genuine somatic events (R1.4, R3.7). The gap is now measured under one protocol and decomposed by stage (R1.6, R3.m2). [待补: ClairS-TO and DeepSomatic-TO recall at PASS and after the shared population filter, HG008T and COLO829T chr1] Pansoma's chr1 PASS SNV recall after population filtering is 0.826–0.831 (HG008T) and 0.686–0.717 (COLO829T). On HG008T, 61 of the 702 truth SNVs are population alleles that any caller using the same filter loses; on COLO829T, only 0.812–0.816 of the truth reaches a candidate, mainly because of low-VAF subclonal SNVs below the candidate allele-fraction or read-support thresholds (R2.4).

*The two categories.* "Variant node not ranked #1" and "first aligned node mismatch" belonged to the earlier pipeline, which we replaced. [待核实: one-sentence definitions from the original Supplementary Fig. 2] Every truth allele that does not reach a candidate now gets one read-level class, on six datasets (later losses: R1.6): no read with MAPQ > 10; fewer than three ALT reads or an allele fraction below the candidate threshold; the ALT is an existing graph allele; repeat absorption or respelling; an indel longer than 50 bp. For example, of 398 HG008T PacBio truth SNVs without a candidate (chr1–22), 65 are graph SNP alleles and 239 lie in repeats where the reads absorb or respell the event; of its 2,976 INDELs without a candidate, 2,739 (92.0%) are fully or partly encoded in the graph.

*The normal at missed sites on graph nodes.* Instead of realigning normal reads, we used haplotype-resolved matched-normal assemblies (HG008-N v6.2; COLO829BL), which show the normal's sequence directly. Of the truth alleles whose tumor reads align without edit to an existing graph path (HG008T 2,478; COLO829T 497), only 110 (4.4%) and 34 (6.8%) are carried by the matched normal. 2,085 (84.1%) and 430 (86.5%) are absent from both normal haplotypes and lie on nodes walked by other HPRC individuals, mostly as one-unit homopolymer or STR changes. At the 1,935 HG008T INDEL loci absent from the normal, the HG008-T tumor assembly carries a sequence found in neither normal haplotype at 1,931 and exactly the truth allele at 1,652, and the tumor's Illumina reads show the expected allele at 96.2% of them (95.7% at randomly chosen non-absorbed truth INDELs). On HG008T these are real somatic events hidden by the graph, and the manuscript now says so. COLO829T has no tumor assembly; there, up to 15.8% of the absorbed INDELs may be germline-like (a normal haplotype has the ALT length with other bases). [可选 E20: share of these sites where COLO829BL ONT reads take the alternative path]

**证据（内部）。**

- recall 数字：同 R1.6（`analysis/benchmark_7fc7396_20261012/hg008t_snv_indel_7fc7396_20261012/tables.md`；COLO829T 用 `experiments/E8_comparators/tables.md` 的 Pansoma d9 行）：HG008T SNV PoN 后 PASS R 0.831 / 0.828 / 0.826，truth 有 tensor 0.945 / 0.944 / 0.951，PoN 标掉 61 of 702；COLO829T SNV PoN 后 R 0.711 / 0.717 / 0.686，truth 有 tensor 0.816 / 0.816 / 0.812。
- COLO829T 低 VAF：`analysis/tensor_recall_20260930/tables.md` §5（truth VAF < 0.06 的 SNV 74–82 % 被 AF floor 过滤）和 §7（chr1 2,443 个 truth SNV 缺 380 / 380 / 366（fiberseq / ONT / Illumina），其中 filtered AF 271 / 157 / 232、<3 reads 69 / 172 / 83、no candidate 40 / 51 / 51；ONT 以 <3 reads 为主，所以 Response 写 "allele-fraction or read-support thresholds"）；这是 09-30 的 build（SNV AF 0.06），现在是 0.08，只作定性引用。
- read-level 类别：`analysis/tensor_recall_20260930/README.md` §2：HG008T PacBio SNV no candidate 398，S3 65、S4a+b+c 239；INDEL no candidate 2,976，I1+I2 2,739（92.0 %）；六套 I1+I2 = 81–92 %。（no-candidate 不依赖 AF floor；`graph_absorbed_somatic_20261001/README.md` 说明 10-01 rebuild 不会改变这个集合。）
- 组装证据 HG008T：`analysis/graph_absorbed_somatic_20261001/README.md`：2,478（201 SNV + 2,277 INDEL）；110（4.4 %）是 HG008-N germline；2,085（84.1 %）两个 normal hap 都没有，其中 1,854 有 HPRC haplotype 携带完全相同的 allele；1,612 of 1,935 absent INDEL 是 HP ≥ 7 的一碱基变化，209 是 STR 一个 unit；Illumina tumor reads 在 96.2 % 的 absent INDEL 位点显示 tumor allele（随机非吸收 truth INDEL 95.7 %）。
- 组装证据 HG008-T：`graph_absorbed_somatic_20261001/tables.md` "HG008-T shows no change at the array"：HG008N_absent_HPRC_other INDEL False 1,932 / True 3；"HG008-T contigs vs normal haps"：novel_ALT INDEL 1,652、novel_other 279（另 all_eq_hap1 1、no_change_both 1、no_change_hom 1、NA 1，共 1,935）。"neither normal haplotype" = novel_ALT + novel_other = 1,931（和 R3.7 一致）；旧稿写 "differs from both normal haplotypes at 1,932" 多算了 all_eq_hap1 那 1 个（1,932 是 "shows a change"，不是 "两个 normal hap 都没有"），已改成 1,931。
- read-level 类别的措辞：S4a / b / c 的定义是 "the SNV becomes an edit on another branch / is absorbed by a graph path / no clear ALT"（`tensor_recall_20260930/README.md` §2），不全是 absorption，Response 改成 "lie in repeats where the reads absorb or respell the event"。这些数字是 09-30 的 tensor set、chr1–22。
- 组装证据 COLO829T：`analysis/graph_absorbed_somatic_colo829t_20261005/README.md`：497（67 SNV + 430 INDEL）；germline filtering 34（6.8 %）；pangenome-induced 430（86.5 %），其中 390 有 HPRC haplotype 携带完全相同的 allele；吸收的 INDEL 里 325（75.6 %）是整 unit 的 tandem repeat 长度变化。同一 README："the tumor allele was not assembled ... up to 68 (15.8%) of the absorbed INDELs may be germline-like, not only the 33 (7.7%) where the normal carries the ALT exactly"（另 35 个 INDEL 上 COLO829BL 某个 hap 有相同 ALT 长度但碱基不同）。所以 "real somatic events" 只对 HG008T 有 tumor 组装支持，Response 已把 COLO829T 的这个限制写进去。
- 和 PoN 的关系（R1.4 的数字，可在本条引用）：`tmp/graph_absorbed_pon_20261012/README.md`：被吸收的 truth 被 PoN 标掉 HG008T SNV 199/201、INDEL 2,156/2,277（94.7 %），COLO829T 66/67、417/430（97.0 %）；图特有损失 HG008T 2 SNV + 121 INDEL（9,362 个 INDEL 的 1.3 %），COLO829T 1 + 13（0.6 %）。
- 比较器预览（不是同协议的正式结果，不写进 Response）：`experiments/E8_comparators/README.md` §6.1，COLO829T fiberseq chr1 SNV：本地 DeepSomatic 1.9.0 tumor-only + PON：PASS P / R / F1 0.592 / 0.564 / 0.578，PASS + 我们的 PoN 0.733 / 0.557 / 0.633，全部记录 + PoN 曲线 best F1 0.792（P 0.724，R 0.873），ceiling 0.894；SMaHT DeepSomatic 1.6.0（header 没记录 mode，README 推测像 tumor-normal [unverified]）PASS 0.972 / 0.745 / 0.844。Pansoma d9 fiberseq PoN 后 0.857 / 0.717 / 0.781。
- 出入：详版 / 简版 R2.3 的 "HG008-T tumor assembly carries the ALT at 95 % of graph-encoded HG008T INDEL truths" 在 `graph_absorbed_somatic_20261001` 里找不到出处；按 `tables.md` "HG008-T any contig == ALT over the array"，所有吸收的 INDEL 里 True = 75 + 23 + 1,652 + 82 = 1,832 of 2,277（80.5 %，本文件计算），所以 Response 改用 1,931 / 1,652 of 1,935。
- 旧类别的定义：repo 里 grep "ranked #1" / "first aligned node" 只在两份回复稿里出现，没有原始定义。
- 另一类不在 GRCh38 truth 里的漏检（R3.3 / R3.4 / R1.1 主讲，不写进本条 Response）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §5，HG008T 的 assembly 验证 off-reference somatic SNV（somatic allele 落在 GRCh38 不经过的 d9 node 上，全常染色体，都不是训练 label）：不改 ch6（base）时 p_somatic ≥ t 的只有 0 / 112（PacBio d9）和 1 / 89（ONT d9）；ch6 全设成 H（all_H）后 61 / 112、49 / 89；predict 实际用的 off-reference rescue（site ±10 列抬到 H，按 argmax 判，= `site10_high`）argmax somatic 32 / 112、46 / 89。即 graph node 上的 somatic 也会因为 path count 低被压掉，rescue 只救回一部分。这些不计入 benchmark 的 FN，所以不影响和 DeepSomatic 的 recall 差距。

**还缺什么。**

- [数据同步] user 提供已有 DeepSomatic-TO / ClairS-TO 结果路径，按同一协议 re-score HG008T 和 COLO829T chr1，得出 recall 差距（E8 re-score；不重跑）。注意确认原稿用的 DeepSomatic 是 tumor-only 还是 tumor-normal（SMaHT 1.6.0 VCF 像 tumor-normal）。
- [CPU 本地] E1.6 漏斗里加比较器的 truth → PASS → PoN，按 stage 分解差距；E5c：比较器在 2,478 / 497 个吸收位点上的 recall。
- [文字] 从原稿 Supplementary Fig. 2 找回 "variant node not ranked #1" / "first aligned node mismatch" 的定义并各写一句；撤回 truth 错误的解释；Supp. Fig. 2 换成新的 read-level 类别。
- [用户决定] 是否做 E20（COLO829BL ONT reads 本地，约 1 天；HG008-N reads 要下载），还是只用组装回答 normal 的问题。
- [CPU 本地] 可选：在现行 build（SNV AF 0.08）上重算 COLO829T chr1 缺失 truth 的 VAF 分布，替换 09-30 的定性引用；三个 d46 SNV run 的 chr1 call VCF 和 rtg fn 已在 `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/`（COLO829T ONT d46 可直接用）。
- [CPU 本地] 可选：E5b graph spelling（正在跑，今天出结果）给出 d46 和 HCC1395T 上的 absorbed 分类；出来后可在本条补一句 d46 上的吸收情况（现在 `E5b_absorbed_d46/tables.md` 标 PRELIMINARY，absorbed 列全是 0）。

### R2.4 Low-VAF detection and the low-VAF model

**状态：部分解决**: 已经有三类结果：candidate 阈值、COLO829T 按 truth VAF 分 bin 的 candidate 保留率（来自 09-30 的旧 build，SNV AF 下限 0.06，现在是 0.08，要在现在的 build 上重算）、COLO829T 三平台 purity joint model 按 purity 的结果（fiberseq 干净；同一肿瘤 in-sample，unplaced calls 不评）。但 reviewer 最直接问的两件还没有：call 层面按 truth VAF bin 的 sensitivity（"主模型在多低的 VAF 下变差"），和替代旧 low-VAF model 的低 AF 配置（E15）。HCC1395 dilution 也没开始，Illumina / ONT dilution 带 read-source BQ artefact、分数偏乐观（要在 RIS 上修）。call 层面 VAF bin 的输入现在本地已有一套（COLO829T ONT d46 的 in-sample GPU 重新预测，已复现 RIS），可以先做这一套。

**Comment（要点）。**

- 10 % 的 read-support 阈值会限制低 purity、高异质性真实样本里的低 VAF 检测：真正的低 VAF SNV 应该怎么测？
- 专门的 low-VAF model 描述不清：它和主模型有什么不同？
- 它能不能找回低 VAF 变异？主模型在多低的 VAF 开始变差？要求做 VAF 递减的 sensitivity 分析。
- 两个模型各自在什么情况下用？

**Response.**

*Thresholds.* A candidate now needs an allele fraction of at least 0.08 and at least three ALT-supporting reads (base quality ≥ 10), for SNVs and INDELs alike. This floor sets the lower VAF limit: below it, a variant never reaches the classifier. On undiluted COLO829T (SMaHT truth VAF, chr1–22, three platforms), the share of truth SNVs that become a classifier input depends on VAF. These shares come from an earlier build with an SNV floor of 0.06; at the current 0.08 floor the shares below 0.10 can only be lower [待补: the same table on the current build]:

| truth VAF | share that becomes a classifier input |
|---|---|
| < 0.06 | 13.5–18.3% |
| 0.06–0.10 | 54.8–73.9% |
| 0.10–0.20 | 89.3–95.0% |
| 0.20–0.30 | 95.8–97.4% |
| ≥ 0.30 | 99.2–99.5% |

[待补: recall per truth-VAF bin after calling, for every dataset.]

*Decreasing purity.* For each platform, one model was trained on chr2–22 of COLO829T and its in-silico dilutions with COLO829BL (nominal purity 75, 50, 25%) and tested on chr1 of each purity at one validation threshold; training and test share the tumor. On fiber-seq, after population filtering (calls without a GRCh38 position are not scored; R3.4):

- SNV recall is 0.698 / 0.639 / 0.621 / 0.663 from 100% to 25% purity, at precision 0.881 / 0.741 / 0.715 / 0.786 (F1 0.779 / 0.687 / 0.665 / 0.720).
- The share of truth SNVs with a candidate falls from 0.821 to 0.711.
- At 25% purity (median allele fraction of the somatic candidates 0.188), recall (0.663) is close to the highest recall the curve reaches after filtering (0.688), so recall there is limited mainly by candidate generation.
- At 75% and 50% purity the classifier loses more (recall 0.639 / 0.621 against curve maxima of 0.767 / 0.740) and precision falls; the extra false calls are germline SNVs in the tumor's regions of loss of heterozygosity (R1.5).

[待补: Illumina and ONT after removing a read-source base-quality cue that makes their current scores likely optimistic (R1.5); the same models on HCC1395 dilutions, a tumor not seen in training.]

*The low-VAF model.* The original low-VAF model was trained on the HapMap mixture and is withdrawn with it (R1.5). [用户决定: either] it is replaced by a low-AF configuration (candidate floor AF ≥ 0.03; otherwise the default network and recipe, trained on the dilutions) [待补: recall and precision per VAF bin against the default, and the increase in candidates]; [or] the dilution-trained models are our only low-purity models.

*Usage.* We recommend the default model for expected somatic VAFs of about 0.10 and above, where at least 89.3% of truth SNVs became candidates in the earlier build [待补: confirm on the current build and with call-level recall per VAF bin]. Below about 0.06, most variants never become candidates, and [待补: the low-AF configuration or a dedicated ultra-low-VAF method] is needed. [待补: the VAF at which call-level recall drops below 50%.]

**证据（内部）。**

- 候选阈值：`indexed_gam_pipeline_v4/README.md` 第 54、74 行写 `--snv-min-af 0.08 --indel-min-af 0.08`，第 238 行写 `--min-variants 3`（"minimum ALT-supporting records"），第 239 行写 `--min-allele-bq 10`（"minimum base quality of an ALT observation"，Response 的 "base quality ≥ 10"）。COLO829T 的 set 是 build AF 0.08（benchmark colo829t README 的 "Data" 段）。
- 按 truth VAF 分 bin 的 candidate 保留率：`analysis/tensor_recall_20260930/tables.md` §5（COLO829T chr1–22，Illumina 用 VAF_Ill，fiberseq / ONT 用 VAF_PB）。
  - 三个平台 built 比例（< 0.06 / 0.06–0.10 / 0.10–0.20 / 0.20–0.30 / ≥ 0.30）：fiberseq 13.5 / 73.9 / 95.0 / 95.8 / 99.2 %；ONT 15.3 / 54.8 / 89.3 / 97.1 / 99.5 %；Illumina 18.3 / 66.4 / 93.3 / 97.4 / 99.2 %。核对时把漏掉的 0.20–0.30 一行（95.8–97.4 %）补进了 Response 的表。
  - "built" = 有 tensor（过了 candidate filter），"filtered" = 有 candidate 但没过 AF / read 数下限。原稿表头 "share reaching a candidate" 会把 filtered 也算进去，核对时改成 "share that becomes a classifier input"。这是未稀释的 COLO829T（P100），低 VAF 的 truth 是 subclonal 的。
  - 这是 2026-09-30 的 tensor 集，早于 10-01 的 ch6 rebuild。核对时发现原稿这里写错了：`analysis/tensor_recall_20260930/README.md` 第 18 行写的是 "Build thresholds on every set: SNV AF 0.06, INDEL AF 0.08"，不是 0.08 / 0.08；现在的 set 是 build AF 0.08（SNV 和 INDEL，benchmark colo829t README "Data" 段）。AF 下限从 0.06 提到 0.08 只会多滤掉 candidate，所以 0.06–0.10（以及 < 0.06、0.10–0.20 里 read AF 抽样落到 0.08 以下的）比例在现在的 build 上只会更低。Response 的表前已经写明 "earlier build with an SNV floor of 0.06"，Usage 一句也改成 "became candidates in the earlier build"。
  - 详细版写的 "74–82 % filtered（VAF < 0.06）" 和这里一致（filtered 73.9–82.3 %）。
- 按 purity 的结果：`analysis/benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012/tables.md` 的 per-purity 表。
  - fiberseq SNV PoN P / R：0.881 / 0.698、0.741 / 0.639、0.715 / 0.621、0.786 / 0.663。
  - "PoN ceiling" 列（PoN 后 curve 能到的最高 recall，rtg report 的 "ceiling"）：0.795、0.767、0.740、0.688。P25 的 recall 0.663 离 0.688 只差一点；P75 / P50 的 0.639 / 0.621 离 0.767 / 0.740 差得多。原稿 Response 的 "So down to this allele fraction, recall is limited by candidate generation, not by the classifier" 只对 P25 成立，核对时改成只说 P25，并加了 P75 / P50 那一条（原因见 README §3 的 LOH，细节在 R1.5）。
  - PASS 阈值是 joint checkpoint 的 validation recall-0.9，四个 purity 共用一个（fiberseq 0.549）；训练和测试是同一个肿瘤（chr2–22 / chr1），dilution 和 P100 共用 reads，所以不是跨肿瘤的结果（Response 加了 "training and test share the tumor"）。
  - Unplaced：fiberseq 每个 purity 有 705 / 547 / 627 / 688 个 PASS call 没有 GRCh38 位置，rtg 不评（`runs/COLO829_fiberseq_SNV_joint_7fc7396_gpu4/vcf_chr1_*.eval_bed.txt`）；全算 FP 的话 PoN 后 precision 是 0.646 / 0.589 / 0.552 / 0.590。Response 加了 "calls without a GRCh38 position are not scored"。
  - PoN F1：0.779 / 0.687 / 0.665 / 0.720。
- 同目录 `purity.md`：
  - fiberseq SNV with a tensor 0.821（P100）到 0.711（P25）；
  - P25 的 somatic SNV tensor AF 中位数 0.188（read AF，只算有 tensor 的，所以偏高；Response 写 "median allele fraction of the somatic candidates"）；
  - measured purity 0.66 / 0.40 / 0.20。
- Illumina / ONT dilution 是临时数字：同目录 README §2（read-source BQ artefact：COLO829BL normal reads 的 BQ 比 tumor reads 高，Illumina 平均 BQ 29.43 对 30.77；P25 时只用这一项的 AUC Illumina 0.706、ONT 0.689；README 的结论是分数 "likely optimistic for a real low-purity tumor"，fiberseq 没有这个 cue）。Illumina P100 因为 pooled threshold 只有 0.675，chr1 best 是 0.784（README §1）。
- 本地重新预测可以代替从 RIS 同步：
  - `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/vcf_chr1/eval_bed/report.txt` 的 pon 行 F1 是 0.603 / 0.575 / 0.736；RIS 同一模型（`experiments/E12a_cross_sample/tables.md` 的 in-sample 列）是 0.605 / 0.575 / 0.737。
  - `experiments/E1.1_ch6_counterfactual/results/tables.md` §2：18 个 SNV 模型的 base mode，CPU fp32 和 GPU bf16 的 somatic AP 相差不超过 0.005（例如 COLO829_ONT_joint 0.6392 / 0.6392）。
  - COLO829T ONT d46 的 `eval_bed/truth.SNV.vcf.gz` 有 VAF_Ill / VAF_PB，`vcfeval/pon_calls` 和 `raw_calls` 里有 `tp-baseline.vcf.gz` / `fn.vcf.gz`，所以这一套的 VAF bin 不用再跑任何 job。HG008T 的 truth 没有 VAF 字段。
  - 注意这一套有 4,569 个 PASS call 没有 GRCh38 位置（同一 report.txt；全算 FP 的话 PoN 后 precision 0.252，报的是 0.787）。按 VAF bin 算 recall 不受影响，但 precision 不能按 bin 报。
- 跨肿瘤（只在 100 % purity）：`experiments/E12a_cross_sample/README.md`（preliminary）。HG008T 的 d46 模型在 HCC1395T chr1 上 PoN F1 0.790 / 0.698 / 0.689（Illumina / PacBio / ONT）；COLO829T 模型的 PASS 阈值不能直接搬（Illumina → HCC1395T recall 0.157，best F1 0.810）。这说明换肿瘤以后阈值要单独校准，和 Response 里 "same models on HCC1395 dilutions" 的 [待补] 相关；dilution 层面的跨肿瘤结果还没有。E12a 用的是旧 checkpoint（BQ / MAPQ 按训练肿瘤的统计量 z-score），per-dataset normalization（f2a0e4e，原 id cbc94a7）只对新训练的模型有效，所以这些跨肿瘤数字是 preliminary，Response 没有引用。
- 训练正样本里低 VAF 很少：`experiments/E21_training_composition/README.md` §3。somatic SNV 训练 tensor 里 AF < 0.20 的，HG008T 占 0.4–1.2 %，COLO829T 占 6.5–7.3 %；README 原话是 "the COLO829T dilution sets are the only source of them"。
- 撤回 HapMap-trained low-VAF model：依据是详细版 R2.4 的评估；旧训练脚本 `/scratch/jshen/Scripts/HapMap_SNV_train.sh`。
- SEQC2 titration 的 purity 档位：E8b Q4–Q6（75 / 50 / 20 / 10 / 5 %）。本地那个 run 是 50 %（SRR7890936），解决了详细版写的 "purity 档位不明"。
- 低 AF build 的可行性：HapMap 重跑用的就是 builder AF 0.015、`--min-variants 3`（`tmp/hapmap_nonHG005_20261012/make_tensors.sh` 第 7–8 行），说明选项现成，不用改代码。

**还缺什么。**

- [CPU 本地] 在现在的 build（AF 0.08）上重算 COLO829T 按 truth VAF 分 bin 的 candidate 保留率（`tensor_recall_20260930` 的 §5 是 SNV AF 0.06 的旧 build），替换 Response 的表；输入是各 set 的 `tensors/somatic.recall.tsv` 和 truth 的 VAF_Ill / VAF_PB。
- [CPU 本地] call 层面按 truth VAF bin 的 recall（E1.6 的 VAF 部分），COLO829T ONT d46 这一套现在就能做：`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_ONT_SNV_d46_7fc7396/vcf_chr1/eval_bed/vcfeval/pon_calls/` 的 `tp-baseline.vcf.gz` / `fn.vcf.gz`，按 truth 的 VAF_PB 分 bin（< 0.05、0.05–0.1、0.1–0.2、≥ 0.2），PoN 前用 `raw_calls`。
- [数据同步] 其余各套（COLO829T Illumina / fiberseq，以及 P100 和三个 dilution 的 joint model）的 chr1 预测。本地 `ris_7fc7396_20261012/<run>/test_chr1_*` 只有 `.metrics.json`，`best.pth` 已经在本地。
  - 两条路：从 RIS 同步 predictions / vcfeval 输出；或者在 tequila GPU 上重新 predict（和 `reviewer_cpu_20261013` 一样，已证明能复现 RIS；按用户 10-13 的规则 inference 不用 CPU；要用户同意再提交）。
  - HG008T 的 truth 没有 VAF 字段（只有 AC / AN / HG008Nv62SOMATICVARIANT），要先定 VAF 来源（例如 tensor AF），再分 bin。
- [GPU RIS] 修掉 read-source BQ cue 以后重训 Illumina / ONT 的 purity joint model，按每个 purity 的 validation 定阈值（E14b）。
- [用户决定] 做不做 E15 低 AF 配置。build 用 `--snv-min-af 0.03`，tensor 预计多 3–5 倍；如果做，训练要 [GPU RIS]。不做的话，Response 写 "dilution-trained models only"。
- [用户决定] 第二个肿瘤用 E14c（in-silico HCC1395，约一周集群时间），还是 E14d（SEQC2 湿实验 titration，需要下载 75 / 20 / 10 / 5 % 各档，没有 read-source cue）。做完以后的预测在 tequila GPU 上跑（用户 10-13 的规则：inference 不用 CPU），训练是 [GPU RIS]。
- [文字] Methods 里写清楚候选阈值、dilution 的构建方法，以及用法建议的 VAF 分界（等 call 层面的 VAF bin 出来再定）。

### R2.5 What drives the higher precision

**状态：部分解决**: 第一问（精度提升是不是主要来自 germline contamination 减少）有三个参考臂在同一批 reads 上的结果，加上 18 个 SNV 模型的 ch6 counterfactual（E1.1），两步都有数字：图先减少 germline candidate，ch6 再挡掉剩下的一部分。E9a / E9c / E9d 已在三个 d46 SNV run（HG008T PacBio、HG008T ONT、COLO829T ONT，chr1；COLO829T ONT 只处理了随机抽的 3,000 / 4,569）上完成：unplaced PASS calls 大多也是 germline（只看组装，57.6–83.1 % 在病人 normal assembly 里；按 population > germline > tumor-only 的互斥优先顺序是 population 49.8–89.2 %、germline 7.9–45.4 %、tumor-only 0.2–0.4 %，COLO829T 没有 tumor assembly），投影到 GRCh38 几乎不改变 F1。worst case 下 "linear 更高" 的 4 个比较里 d46 的 2 个来自 off-ref rescue：不用 rescue 时 complete-output F1 至少 0.596 / 0.524 / 0.721，高于 linear 的 0.359 / 0.365 / 0.418；用 rescue、扣掉 graph PoN 标的是 0.590 / 0.489 / 0.491，但 graph PoN 也标掉 chr1 上大部分 gold off-ref somatic key（9/11、7/9、2/2，gold 很少），所以只当参考。PoN 后剩下的 placed FP 按 assembly 大多是 germline、而且多数不在 HPRC haplotype 里（d46 臂，每组抽 150）。只有三个 run，都是长读长 SNV。还缺：rescue 默认值和 graph PoN 的用法（用户决定）；另外 2 个反转（d9，要先 re-predict）；第二问（附近有多个 germline 时 graph 是否更有帮助，E19）没有任何数据；E4b 的 linear 臂拆分；去掉 ch6 的重训（GPU RIS）。

**Comment（要点）。**

- Is the precision gain mostly reduced germline contamination, or are there other factors?
- Beyond modelling germline variation at the site, does graph-aware calling help when several germline variants are nearby?

**Response.**

Mostly reduced germline contamination, in two steps: the graph removes most germline candidates at alignment, and the classifier rejects many of the rest using the number of HPRC haplotypes through each node (ch6). We compared a GRCh38-only graph, HPRC v1.1 and v2.1 on the same reads, code and recipe (R3.1), and edited ch6 at inference in all 18 SNV models.

*Fewer germline candidates.*

- 75.7–82.2% of matched-normal germline variants become candidates on the GRCh38-only graph, against 7.7–13.1% on v1.1 and 7.2–12.4% on v2.1 (R1.2). Before filtering, the HG008T chr1 models call 117,205–144,251 germline tensors somatic at PASS on GRCh38-only, 4,864–5,147 on v1.1 and 2,528–2,810 on v2.1 (somatic average precision 0.012–0.016, 0.199–0.233 and 0.310–0.353).
- After filtering (benchmark regions), v1.1 still leads by 0.06–0.12 SNV F1 on HG008T and 0.10–0.31 on COLO829T (v2.1: 0.04–0.10, 0.10–0.32). On HG008T recall is similar (0.79–0.84 against 0.83 for v1.1) and 1.4–1.7 times more false calls survive filtering on GRCh38-only (1.8–7.8 times on COLO829T).
- Graph-coordinate calls, which rtg cannot score, are mostly germline too. Counted as false, they put GRCh38-only ahead in 4 of 12 SNV comparisons, all long-read (R3.4). In three v2.1 SNV runs (HG008T PacBio and ONT, COLO829T ONT; chr1; COLO829T: 3,000 of 4,569 calls sampled), 2.8–7.9% can be written as a GRCh38 allele, the patient's normal assembly carries 57.6–83.1%, and HPRC haplotypes or the graph spell 49.8–89.2%; tumor-only calls are 1.7–3.8% (HG008T), most of them also population-tagged, leaving 0.2–0.4%. All come from the off-reference rescue (R3.4): without it, complete-output F1 is at least 0.596 / 0.524 / 0.721 against 0.359 / 0.365 / 0.418 on GRCh38-only, so both v2.1 reversals disappear. Discounting the HPRC-spelled calls instead gives 0.590 / 0.489 / 0.491, a reference only, as this test also tags 9 of 11, 7 of 9 and 2 of 2 assembly-verified off-reference somatic keys on chr1. [待补: the rescue default; the two v1.1 reversals]

*Haplotype counts.* Germline tensors left on the graph sit on less-shared nodes (median site code 95–96 on v1.1 and 74–80 on v2.1, 100 for somatic; AUC 0.837–0.892 for this code alone, three runs). Setting ch6 to "every haplotype" multiplies the germline tensors called by 3.41–7.57 and lowers somatic AP from 0.199–0.352 to 0.061–0.087 on HG008T and from 0.639–0.829 to 0.347–0.516 on COLO829T (other ch6 edits lower it about as much or more). After filtering the effect is smaller (F1 −0.063 to +0.007): the filter tags nearly all the extra germline calls, and most extra false calls left are not germline. This tests use, not benefit [待补: retraining without ch6 (R1.1)].

*Other factors.* Matched-normal germline variants that the resources do not tag are 79.2–97.9% of the false calls left after filtering on HG008T and 35.1–77.8% on COLO829T (per-tensor estimates). In the same three v2.1 runs the normal assemblies carry the ALT of 146, 135 and 110 of 150 sampled false calls, and only 2.0–10.0% are spelled by HPRC v1.1 haplotypes or the v2.1 graph: they are mostly germline alleles the graph lacks. Before filtering, GRCh38-only also passes more tensors that are neither somatic nor germline (HG008T 936–1,993 against 24–229 on v1.1). [待补: the same split for the GRCh38-only and v1.1 arms (E4b, with R3.2)]

*Clustered germline variation.* [待补: precision and recall by the number of matched-normal germline variants within ±50 bp, HPRC against GRCh38-only, with block-bootstrap intervals (E19)]

*Limitations.* The GRCh38-only arm is our network, not a dedicated linear caller (R3.m2); each arm has one seed, and some v1.1 models differ in batch size or learning rate [待补: matched reruns, seeds (R3.9)].

**证据（内部）。**

- 三臂 candidate 层面（germline truth 有 tensor 的比例）：`response_to_reviewers_detailed.md` R1.2 表（2026-10-12，`tensors/truth_recall.json`）：linear 0.757–0.822、d9 0.077–0.131、d46 0.072–0.124；peer 的 `points/R3.2.md` 已重跑核对一致。
- 三臂 call 层面：`experiments/E3_reference_arms/tables.md`（来源 `analysis/benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/summary.json`）。
  - §3：HG008T somatic AP linear 0.012 / 0.012 / 0.016，d9 0.199 / 0.233 / 0.232，d46 0.310 / 0.353 / 0.328；chr1 SNV test tensors 例 HG008T PacBio 251,303 vs 32,880。
  - §1：HG008T SNV PASS F1 linear 0.469 / 0.513 / 0.537，d9 0.585 / 0.621 / 0.596，d46 0.570 / 0.605 / 0.575；COLO829T linear 0.468 / 0.680 / 0.418，d9 0.778 / 0.781 / 0.706，d46 0.783 / 0.779 / 0.737。linear README："HG008T -0.06 to -0.12, COLO829T -0.10 to -0.31 against d9"。
  - summary.json（eval_bed，pon）：HG008T P linear 0.3344 / 0.3719 / 0.3937 vs d9 0.4516 / 0.497 / 0.4656；R linear 0.7862 / 0.8264 / 0.8422 vs d9 0.8307 / 0.8278 / 0.8264。FP linear 1,091 / 973 / 904 vs d9 703 / 584 / 661；COLO829T 2,216 / 522 / 2,575 vs 284 / 293 / 630。linear README："1.4-1.7x d9's on HG008T, 1.8-7.8x on COLO829T"。
  - d46 − linear（E3 §1）：HG008T +0.101 / +0.092 / +0.038，COLO829T +0.316 / +0.099 / +0.319，所以 Response 写 "v2.1: 0.04–0.10 and 0.10–0.32"；原稿 "HPRC graphs lead by 0.06–0.12" 只对 d9 成立（linear README 也写 "against d9"）。
  - BED 和 unplaced：§1 是 BED 内、unplaced PASS calls 不计分。E3 §2 complete-output bound（不限 BED、unplaced 全算 FP）的 SNV F1：HG008T linear 0.197 / 0.359 / 0.365，d9 0.546 / 0.553 / 0.342，d46 0.414 / 0.533 / 0.354；COLO829T linear 0.468 / 0.680 / 0.418，d9 0.711 / 0.636 / 0.522，d46 0.702 / 0.736 / 0.369。linear 高的是 HG008T ONT d9、HG008T ONT d46、COLO829T fiberseq d9、COLO829T ONT d46，12 个比较里 4 个。这是 worst case；三个 d46 run 上的 E9a / E9c / E9d 结果见下一条。
  - d46 vs d9：`hprc_v21_d46_7fc7396_20261012/README.md`："Raw PASS F1 rises by 0.10-0.13 on every SNV set"；E3 §1 d46 − d9 SNV −0.015 / −0.016 / −0.021 / +0.005 / −0.002 / +0.031。
- Unplaced calls 是什么（新，E9a / E9c / E9d 已完成）。三个 d46 SNV run 在 tequila GPU 上 in-sample re-predict（`reviewer_cpu_20261013`），chr1。下面三个数的顺序都是 HG008T PacBio / HG008T ONT / COLO829T ONT。
  - E9a 投影：`experiments/E9a_unplaced_projection/results/tables.md`（还有 `projection_by_run.tsv`、`metrics_by_run.tsv`、`summary.json`）。
    - unplaced PASS 240 / 1,258 / 4,569。class i（能写成 GRCh38 allele：exact 2 / 19 / 94，branch base differs 17 / 24 / 34）19（7.9 %）/ 43（3.4 %）/ 128（2.8 %）；投影后 PoN 标掉其中 18 / 38 / 88 个，PoN 后没有一个是 truth（E9a README Interpretation 2）。class ii（只有区间）221 / 984 / 1,334。class iii（没有 anchor）0 / 231 / 3,107。
    - class ii 在区间里没有任何 truth TP：BED 内 0 / 208、0 / 940、0 / 1334。
    - PoN 后、no BED 的 complete-output F1：unplaced 全算 FP 是 0.533 / 0.355 / 0.369；投影后（class i 交给 rtg，class ii 按区间，class iii 算 FP）是 0.537 / 0.359 / 0.373。BED 内是 0.537 / 0.358 / 0.369 → 0.544 / 0.367 / 0.373。
    - 所以投影后的 F1 基本就是 worst case。HG008T ONT d46 和 COLO829T ONT d46 投影后仍低于 linear（0.365 / 0.418）。Response 只写了 "only 2.8–7.9% can be written as a GRCh38 allele"（旧稿的 "exact GRCh38 allele" 不准：E9a 的 "exact" 只是 class i 里的一个子类，2 / 19 / 94 个，其余是 branch base differs）。
    - CHM13 列全是 0，因为用户已取消 CHM13。
  - E9c / E9d：`experiments/E9c_offref_filter_validation/tables.md` 的 `rcpu__*` 行（collect job 387434，16:07）。同目录 `README.md` 的 Results / Interpretation / Caveats 已经写好（23:04），下面的数和它一致。
    - §1 unplaced PASS 的分类是互斥的，顺序是 population > germline > tumor-only > unverifiable > unsupported：
      - population（graph PoN）214（89.2 %）/ 932（74.1 %）/ 1,495（49.8 %）；
      - germline（normal assembly）19（7.9 %）/ 242（19.2 %）/ 1,362（45.4 %）；
      - tumor-only 1（0.4 %）/ 3（0.2 %）/ 0（0.0 %）。COLO829T 没有 tumor assembly，这一类结构上就是 0，not in normal 是 72（2.4 %）；所以 Response 写 HG008T 的 0.2–0.4 %，不写 "0.0–0.4 %"；
      - unverifiable 2.1 / 5.6 / 2.4 %，unsupported 0.4 / 0.9 / 0.0 %。
      - COLO829T ONT 只处理了随机抽的 3,000 / 4,569（seed 0），它的百分比是这 3,000 个里的。
    - 只看组装、不按优先顺序（README §1 第二张表；`tables.md` §3 的 `unplaced_pass` 行）：
      - germline 158（65.8 %）/ 724（57.6 %）/ 2,493（83.1 %）；
      - tumor-only 9（3.8 %）/ 22（1.7 %），其中 graph PoN 标掉 8 / 19；COLO829T not in normal 197 / 3,000，其中标掉 125（我从 `per_call.tsv.gz` 数过，和 README 一致）。
      - 所以按优先顺序的 tumor-only 0.4 / 0.2 % 偏低：大部分 tumor-only 被算进了 population。Response 写只看组装的 germline 57.6–83.1 % 和 population 49.8–89.2 %，tumor-only 两种都写（1.7–3.8 %，扣掉 population 后 0.2–0.4 %）；按优先顺序的 germline 7.9–45.4 % 只写在状态行和这里。
    - graph PoN 用主规则（adaptive flank、min AC 1、graph 或 read context）。haplotype 用的是 HPRC v1.1 的 89 条（v2.1 的全图 VCF 不在本地），再加 d46 自己的 graph walk。所以 Response 写 "spelled by HPRC haplotypes or the graph"，FP 那句写明 "HPRC v1.1 haplotypes or the v2.1 graph"。用 v2.1 的 haplotype 时标掉的比例可能更高。
    - §2 graph PoN 之后的 F1（tagged 的 unplaced 去掉，剩下的算 FP）：
      - no BED：0.5901 / 0.4889 / 0.4911；+tumor-only 是 0.5908 / 0.4906 / 0.4911。
      - BED：0.5954 / 0.4988 / 0.4911。
      - 和 E3 §2 linear（no BED，RIS）的 0.359 / 0.365 / 0.418 比，三个 d46 run 都更高。原来 worst case 下 linear 更高的两个 d46 比较（HG008T ONT 0.354、COLO829T ONT 0.369）按这种算法反过来。
      - 但这只是一种计数方式，不是验证过的过滤器：graph PoN 也标掉大部分 gold off-ref somatic key（下面 §5），E9c README "Needs" 把 "graph PoN 只做分类还是当过滤器" 列为用户决定，R3.4 也只把这组数当参考。所以 Response 先写不用 rescue 的下界（0.596 / 0.524 / 0.721，见下面 off-ref rescue 一条），再把 0.590 / 0.489 / 0.491 写成 "a reference only"。
      - COLO829T ONT 扣掉之后仍有 51.6 % 的输出是 unplaced（unplaced share left 0.5163），F1 0.4911 远低于 placed 的 0.7364。
      - COLO829T ONT 的 "unplaced left" 2,292.1 是从 3,000 个抽样外推的。
      - 本次 re-predict 的 unplaced 数和 RIS 略有不同：240 vs 244、1,258 vs 1,261、4,569 vs 4,566。
    - §5：graph PoN 也标掉了 chr1 的 gold off-ref somatic keys，9/11、7/9、2/2（有 tensor 的 gold key 只有 11 / 9 / 2 个，其中 unplaced PASS 3 / 5 / 1 个，全被标掉）。这对 F1 没影响，因为 unplaced call 本来就不会被 rtg 算 TP。但它说明 graph PoN 不能当成保留 somatic 的过滤器来用，这里只拿它给 unplaced calls 分类。
      - chr1 的 gold 太少，比例只看方向。全染色体（README §4，主规则）标掉 74.2–88.3 %；chr2-22 校准的整个 grid 里最低也有 42.5–67.0 %（adaptive flank、AC >= 9），预先定的 "gold 上损失小" 标准没有一个规则满足。
    - §3 controls，每组随机抽 150：
      - placed rtg TP（PoN 后）：tumor-only 149 / 145；COLO829T ONT not_in_normal 150。graph PoN 标掉 0 / 1 / 0。说明组装验证在 placed 位点上分得开 somatic。在 off-ref 的 gold 上它没这么准：每个 gold 集有 9–60 个 somatic key 被判成 germline（例 d9 Illumina 60 / 231，README Caveats）。
      - placed FP（PoN 后）：germline 146 / 135 / 110（COLO829T 另有 40 not_in_normal）。graph PoN 标掉 4（2.7 %）/ 15（10.0 %）/ 3（2.0 %）。Response "Other factors" 段用的就是这组数。
      - GRCh38 PoN 标掉的 placed PASS（ctrl_pon）：germline 147 / 140 / 144，graph PoN 只标 17（11.3 %）/ 13（8.7 %）/ 11（7.3 %）。
    - assembly 的 "germline" 允许 secondary hit（`minimap2 -N 50`），所以 ALT context 出现在别的位置（paralog）也算。`per_call.tsv.gz` 里 ctrl_fp 的 germline 带 `multi_copy` 标记的有 2 / 146、16 / 135、19 / 110（我从 per_call.tsv 数的）。
    - 和 E1.1 的出入：COLO829T ONT d46 PoN 后 FP 里的 germline 比例，E1.1（dipcall 标签，per-tensor PoN）是 191 / 402（47.5 %），E9d（normal assembly，rtg FP 抽 150）是 110 / 150。推断差别来自 dipcall 标签覆盖不全，加上 paralog hit，没有验证。所以 Response 删掉了旧稿的 "lowest for ONT, where most are not germline"，HG008T 的范围合并写成 79.2–97.9 %。
    - §7（repeat context）和 §8（跨平台：HG008T ONT 的 unplaced keys 有 190 个（15.1 %）也在 PacBio 上 PASS）属于 R1.3 / R3.3，这里没用。
  - Off-ref rescue（R3.4 的证据）：三个 run 的 unplaced SNV PASS 全部来自 off-ref rescue。不用 rescue 时 complete-output F1 至少 0.596 / 0.524 / 0.721（R3.4 自己算的下界），也高于 linear 的 0.359 / 0.365 / 0.418。
    - 我按 R3.4 证据里的数重算过：F1 = 2TP / (2TP + FP + FN + U)，TP / FP / FN 586 / 672 / 116、603 / 836 / 99、1,690 / 457 / 753，U（不 rescue 时 p ≥ t 的 off-ref tensor）6 / 162 / 98，得到 0.5961 / 0.5237 / 0.7210，一致。U 没算 AF floor，是上界，所以 F1 是下界。
    - Response 现在把它写在 graph-coordinate 那一条里（"All come from the off-reference rescue ... at least 0.596 / 0.524 / 0.721"），Haplotype counts 段原来的最后一句删掉了。rescue 是否保留是用户决定。
- E1.1 ch6 counterfactual（已完成，18 个 SNV 模型：d9 HG008T 三个、d9 COLO829T purity-joint 三个、d46 六个、linear 六个）：`experiments/E1.1_ch6_counterfactual/results/tables.md`（`collect.py` 生成）和 `results/summary.json`（`pending` 为空）。
  - 设置：chr1 test tensors（和 `predict` 相同），阈值 t = 各模型的 validation recall-0.9（PASS 规则）。推理精度（README "How the runs were done"）：`HG008_PacBio`、`HG008_PacBio_d46` 的 base / all_H / train_mean，`HG008_ONT` 的 base，`HG008_ONT_d46` 的 base / all_H 是 CPU fp32，其余（另外 14 个 run 的全部模式和这四个 run 的其他模式）是 tequila GPU bf16。§2：base 模式和 RIS GPU bf16 `predict` 的 somatic AP 差 ≤ 0.005（最大 `HG008_ONT_linear` 0.0202 / 0.0160 和 `COLO829_fiberseq_linear` 0.2832 / 0.2879）。§2 表头写的 "CPU fp32 here" 是 collect.py 的旧说明，README 已注明只对那四个 run 成立。
  - Response 第一段的 AP（0.012–0.016 / 0.199–0.233 / 0.310–0.353）是 E3 §3 的 RIS 数；ch6 段落的 AP 是 E1.1 的 base（例 HG008_PacBio_d46 CPU fp32 0.3523 vs RIS 0.3528；HG008_Illumina_d46 是 E1.1 自己的 GPU bf16 加权估计 0.3121 vs RIS 0.3099），所以两处上限差 0.001。
  - HG008T / COLO829T Illumina 和六个 linear 套是 class-stratified sample（`ch6_counterfactual_sampled.py`），germline / non 的数是加权估计，truth TP 是精确的。
  - Response 里 linear 的 germline ≥ t 用 §4 的 E1.1 base（117,205 / 127,888 / 144,251，GPU bf16，class-stratified 加权估计）；§2 的 RIS `predict` 全量数是 117,879 / 128,084 / 144,625。
  - "Other factors" 段的 PASS 非 somatic 非 germline tensor（PoN 前）：§3 各 run base 行的 "non on-ref >= t"：HG008T linear 1,993（Illumina，加权）/ 936 / 1,561，d9 229（加权）/ 24 / 89（d9 另有 off-reference non 8 / 6 / 10；d46 是 305 / 44 / 197）；阈值 §1 "PASS threshold"：linear 0.064 / 0.040 / 0.042，d9 0.251 / 0.314 / 0.449。阈值低是因为 linear 的排序差（val somatic AP 约 0.01），这是推断，没单独验证。为了控制字数，Response 不再写阈值，只写 non 的个数。
  - §4 "Germline tensors with p_somatic >= t"，base → all_H（§3 的 "x base"）：HG008T d9 4,864 → 17,422（3.58）、4,951 → 17,327（3.50）、5,147 → 19,496（3.79）；HG008T d46 2,542 → 17,149（6.75）、2,528 → 17,256（6.83）、2,810 → 18,668（6.64）；COLO829T d46 827 → 6,262（7.57）、808 → 5,427（6.72）、942 → 6,961（7.39）；COLO829T d9 joint 666 → 2,299（3.45）、1,799 → 6,136（3.41）、1,396 → 5,855（4.19）。
  - §4 somatic AP，base / all_H / train_mean / shuffle：HG008T d9 0.199 / 0.071 / 0.060 / 0.037、0.233 / 0.087 / 0.053 / 0.032、0.232 / 0.071 / 0.044 / 0.043；HG008T d46 0.312 / 0.062 / 0.053 / 0.033、0.352 / 0.072 / 0.066 / 0.028、0.326 / 0.061 / 0.051 / 0.035；COLO829T d46 0.820 / 0.516 / 0.505 / 0.299、0.829 / 0.515 / 0.520 / 0.192、0.755 / 0.405 / 0.421 / 0.246；joint 0.693 / 0.487 / 0.482 / 0.308、0.714 / 0.515 / 0.512 / 0.220、0.639 / 0.347 / 0.329 / 0.192。COLO829T 的 train_mean 有两个 run 比 all_H 略高（0.520 vs 0.515、0.421 vs 0.405），所以 Response 写 "other ch6 edits lower it about as much or more"。
  - §4 "Truth-level F1 at t after the per-tensor PoN"，base → all_H：HG008T d9 0.589 → 0.573、0.621 → 0.566、0.596 → 0.533；d46 0.560 → 0.511、0.605 → 0.561、0.579 → 0.516；COLO829T d46 0.784 → 0.791、0.785 → 0.760、0.743 → 0.705；joint 0.675 → 0.682、0.782 → 0.778、0.706 → 0.664。这些三位小数相减是 HG008T −0.016 到 −0.063、COLO829T −0.042 到 +0.007，Response 合并写成 −0.063 to +0.007。§3 的倍数 3.41–4.19（d9）和 6.64–7.57（d46）也合并写成 3.41–7.57。
  - train_mean 不是单调的：`HG008_ONT`（d9）train_mean 把所有分数压低（somatic >= t 65 / 663，germline 1,365，PoN 后 F1 0.145），所以 Response 不写 "train_mean 增加 germline call"。
  - linear 六个 run：§1 训练 ch6 "100.0 / 100 / 1.00"（每个 cell 都是 100），所有模式和 base 完全一样（§3）。
  - §1 site code 中位数（somatic / germline）：d9 HG008T 100 / 96（三个），joint 100 / 96、100 / 95、100 / 95；d46 HG008T 100 / 80、100 / 79、100 / 79，COLO829T d46 100 / 80、100 / 76、100 / 74。§6 ch6 单独的 site code AUC：0.837（HG008_PacBio）、0.890（HG008_PacBio_d46）、0.892（HG008_ONT_d46），只有这三个 run。
  - PoN 后多出来的 FP 是什么：`summary.json` 的 `germline_after_pon_ge_thr` / `non_after_pon_ge_thr`，base → all_H：HG008_PacBio（d9）germline 571 → 597、non 12 → 194；12 个 graph run 里 non 的增量都大于 germline 的增量（例 HG008_ONT d9 618 → 678 vs 40 → 266，COLO829_ONT_d46 191 → 202 vs 211 → 506）。
  - §5 off-reference somatic（R1.4 / R3.3）：base 0 / 112、1 / 89；all_H 61 / 112、49 / 89。
- PoN 前 FP 的组成（E1.1 base，`germline.ge_thr` / `truth.fp`，百分比是我按计数算的）：HG008T d9 4,864 / 5,101（95.4 %，加权）、4,951 / 4,981（99.4 %）、5,147 / 5,246（98.1 %）；d46 2,542 / 2,869（88.6 %，加权）、2,528 / 2,577（98.1 %）、2,810 / 3,027（92.8 %）；linear 117,205 / 119,198、127,888 / 128,824、144,251 / 145,811（98.3 / 99.3 / 98.9 %，加权）。
- PoN 后 FP 的组成（E1.1 base，per-tensor PoN，只算 on-reference tensor；`germline_after_pon_ge_thr` / `truth_after_pon.fp`，百分比我算的）：
  - HG008T d9 609 / 688（88.5 %，加权）、571 / 583（97.9 %）、618 / 658（93.9 %）；d46 654 / 826（79.2 %，加权）、618 / 642（96.3 %）、672 / 771（87.2 %）。
  - COLO829T d46 191 / 250（76.4 %，加权）、172 / 221（77.8 %）、191 / 402（47.5 %）；joint 60 / 92（65.2 %，加权）、158 / 210（75.2 %）、140 / 399（35.1 %）。
  - per-tensor PoN 在 graph run 上和 rtg 很接近（§2 例 HG008_PacBio 576 / 583 vs rtg 577 / 584，HG008_ONT 575 / 658 vs 576 / 661）；COLO829_ONT_joint 差得多一些（399 vs 595）。在 linear run 上差很多（HG008_Illumina_linear 837 vs 1091，COLO829_Illumina_linear 211 vs 2216，COLO829_ONT_linear 1,063 vs 2575），所以三臂 PoN 后的对比只用 rtg（E3），不用 E1.1 的 per-tensor PoN F1（例 linear 0.526 / 0.579 / 0.593 和 rtg 的 0.469 / 0.513 / 0.537 不一样）。
- 和详版的出入：
  - 详版 "Before population filtering, 99% of HG008T Illumina SNV false calls are matched-normal germline"：最终 d9 Illumina 模型是 4,864 / 5,101 = 95.4 %（加权估计），HG008T 六个 graph 模型是 88.6–99.4 %；Response 不再写这个百分比，改写各臂 germline call 的绝对数（117,205–144,251 / 4,864–5,147 / 2,528–2,810）。
  - 详版 "the false positives left are mostly rare germline alleles (gnomAD AF < 0.001)"：还是找不到最终模型上 gnomAD AF 分层的文件；Response 改成 "matched-normal germline that the resources do not tag"（上面 E1.1 的数），AF 分层留给 E4b。
- 可比性：E3 README Caveats 1–3（d9 的 HG008T PacBio / ONT SNV 用 1 GPU，COLO829T d9 是旧 run；每臂一个 seed；linear 臂 validation somatic AP 约 0.01）。E1.1 的 COLO829T d9 是 purity-joint 模型（阈值在四个 purity 上合并选），不是 E3 里的 `COLO829_*_SNV_h90`。

**还缺什么。**

- [CPU 本地] E4b：
  - d46 臂已经有抽样的答案（E9c §3 controls，三个 d46 run）：PoN 后的 placed FP 按 normal assembly 146 / 135 / 110 of 150 是 germline，graph PoN 只标 2.0–10.0 %。还缺 linear 和 d9 臂的同样拆分，这是 "精度提升来自图里编码的 germline" 的直接证据。
  - 需要每个臂 chr1 的 call set。三个 d46 run（`HG008_PacBio_SNV_d46`、`HG008_ONT_SNV_d46`、`COLO829_ONT_SNV_d46`）已经在 tequila GPU 上 re-predict 过（`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1/`，`eval_bed/vcfeval/pon_calls/summary.txt` 的 None 行 F1 0.6034 / 0.5749 / 0.7364，RIS 0.605 / 0.575 / 0.737），可以直接用；linear 和 d9 臂还要 re-predict（按 user 10-13 的规则 predict 在 tequila GPU 上跑，不用 CPU），再跑 `vcf.sh`。
  - 把三臂 PoN 前后的 FP 和 HG008-N / COLO829BL dipcall 匹配，再按 allele 是否在 d9 deconstruct VCF 里分三类。d46 要先 `vg deconstruct`。E1.1 已经给了 "germline vs 非 germline" 的拆分，E4b 只差 "germline 是否在图里" 和 gnomAD AF 分层。也可以试 E9c 的 ctrl_fp 流程（`e9c_run.py`：graph PoN + 组装比对），但它在 linear 图上能不能跑还没确认。
- [数据同步] E4b 的另一个办法：从 RIS 同步三臂 run 的 `vcf_chr1/`（本地 `ris_7fc7396_20261012/` 只有 `best.pth` 和 metrics）。
- [CPU 本地] E9a / E9c / E9d 在三个 d46 run 上已完成，4 个反转里 d46 的 2 个（HG008T ONT、COLO829T ONT）已经解释。剩下两个是 d9：HG008T ONT（0.342 vs linear 0.365）和 COLO829T fiberseq `_h90`（0.636 vs 0.680）。要先在 tequila GPU 上 re-predict 出 calls，再跑 E9a / E9c / E9d，补 Response 里的 [待补: the two v1.1 reversals]（和 R3.4 共用）。
- [用户决定] off-ref rescue 是否保留为默认（和 R3.4 共用）。它决定 Response 里 complete-output 的主数：rescue 关（至少 0.596 / 0.524 / 0.721），还是 rescue 开（worst case 0.533 / 0.355 / 0.369；扣掉 graph PoN 标的 0.590 / 0.489 / 0.491 现在只当参考）。
- [用户决定] graph PoN 只做 unplaced calls 的分类，还是也当过滤器（E9c README "Needs"）。主规则在全染色体的 gold off-ref somatic key 上标掉 74.2–88.3 %，grid 上没有一个规则满足 "gold 上损失小"。在用户定之前，Response 把 0.590 / 0.489 / 0.491 写成 "a reference only"。
- [用户决定] / [数据同步] HPRC v2.1 全图 VCF（下载要用户同意，E7 的 v2.1 轴也要）。现在 d46 run 的 graph PoN 只用 HPRC v1.1 的 89 条 haplotype（加 d46 自己的图），"only 2.0–10.0% spelled" 和 population 49.8–89.2 % 都可能偏低。
- [CPU 本地]（可选）COLO829T ONT 剩下的 1,569 个 unplaced PASS 也跑 E9c / E9d，把 3,000 的抽样换成全量（等用户定 cluster cap）。
- [CPU 本地] E19：对每个 truth allele 和 FP 数 ±50 bp 内的 dipcall germline，分四档算 P / R，HPRC 臂 vs linear 臂，加 chr1 block bootstrap。依赖 E4b 的 call sets。第二问完全靠它。
- [GPU RIS] 去掉 ch6 重训（`--drop-planes path_count`，d9 / d46 SNV；命令在 `experiments/E1.1_ch6_counterfactual/gpu_nopc_commands.sh`），回答 "ch6 带来多少"（和 R1.1 共用）。E1.1 只证明模型在用 ch6。
- [GPU RIS] d9 的配对重跑（HG008T PacBio / ONT SNV 用 4 GPU；COLO829T d9 SNV / INDEL 用 7fc7396 recipe）和多 seed，让三臂对比同一 recipe。
- [数据同步] 用户给已有 ClairS-TO / DeepSomatic 结果的路径，按同一协议重新评分，回答 "专门的 linear caller" 的对照（R3.m2）。
- [文字] 详版和简版 R2.5 换成上面两步的写法（candidate 减少 + ch6），删掉 "99% … Illumina" 和 "rare germline (gnomAD AF < 0.001)"；ch6 的段落和 R1.1 / R3.3 的 E1.1 数字对齐。
- [文字] unplaced calls 的分类和 graph PoN 后的 F1 要和 R3.3 / R3.4 用同一组 E9c 数，写法也要一致（互斥优先顺序要写明；graph PoN 后的 F1 只当参考，和 R3.4 一样）。E9c README 的结果段已写好，本次已按 `tables.md` 和 README 核过 R2.5 的数。

### R2.6 Panels of normals

**状态：部分解决**: 第三问（HG008T 不加 PoN）已有最终数字；第一问由 E8b 的文献核实、私有 allele 测试和六套 somatic truth 的 truth-in-resources 计数回答（只差 HapMap）；第二问现在有最终模型上的数（E1.1：PoN 去掉多少 PASS 的 germline call、剩下的 FP 里 germline 占多少），但还没按 "allele 是否在图里" 拆（E4c）；INDEL PoN 规则是在 chr1 上选的，validation 上重选（E17c）要等 RIS 的 predictions。

**Comment（要点）。**

- Some PoNs were built from common benchmark datasets and include variants of the tested samples (HapMap in particular), which may inflate the reduction in germline contamination.
- Fig. 2e shows that > 20% of germline variants are not in the AF-HPRC graph: were they removed by the PoNs?
- Only COLO829 was evaluated without the PoN; add HG008T without the PoN.

**Response.**

*1. Overlap with the benchmark samples.* Our four resources are the files ClairS-TO uses by default: the GATK gnomAD allele-frequency file (a 2017 GRCh38 liftover of gnomAD v2 genome frequencies), dbSNP build 138 (2013), the GATK 1000 Genomes panel (2017) and CoLoRSdb v1.1.0 (2024; 1,381 long-read genomes). HG008 data collection began in 2023, so HG008 cannot be in the first three, and CoLoRSdb's listed cohorts contain no GIAB cancer sample. As the gnomAD and 1000 Genomes panel sample lists are not public, we matched each individual's rare germline SNVs to the resources: the normals of HG008, COLO829, HCC1395 and three CASTLE lines match least of the 12 individuals tested (gnomAD 0.6–2.5%, CoLoRSdb 4.0–7.9%), against 3.4–15.4% and 14.5–25.9% for the HapMap donors. The HapMap donors are therefore likely represented, notably in CoLoRSdb, one more reason to remove that benchmark (R1.5). The somatic truth alleles found in the resources (tagged: HG008T 10.2% of SNVs and 61.8% of INDELs; COLO829T 3.6% / 52.0%; HCC1395 2.8% / 7.3%) are population alleles that recur by chance; no allele of the six somatic truth sets we use (HG008T, COLO829T, HCC1395 and three CASTLE lines) is flagged somatic in dbSNP [待补: HapMap (E8b)]. They cost recall and do not inflate precision.

*2. Germline variants not in the graph.* In the benchmark regions 13% of the matched-normal germline SNVs and 25% of germline INDELs are not HPRC v1.1 graph alleles. Such variants remain candidates and are labelled germline in training. Population filtering tags 95.5–96.6% of the germline SNV candidates on all six datasets. With the final HG008T chr1 SNV models, it removes 87.5–88.5% (v1.1) and 74.3–76.1% (v2.1) of the germline candidates called somatic at PASS (per-tensor estimate, placed calls); the 571–618 and 618–672 left are 88.5–97.9% and 79.2–96.3% of the false calls after filtering. So untagged germline variants, not artefacts, limit HG008T precision (R2.5). [待补: the same split by graph membership (E4c)]

*3. HG008T without filtering*, now reported for every dataset and reference. HG008T SNV, HPRC v1.1, chr1 in the benchmark regions (calls without a GRCh38 position not scored; R3.4): F1 0.182–0.191 without filtering and 0.585–0.621 with it (precision 0.102–0.108 against 0.452–0.497; recall 0.855–0.862 against 0.826–0.831). INDEL F1 at PASS 0.011–0.101 against 0.015–0.241; at the chr1 curve's best point (test-tuned, an upper bound) 0.090–0.150 against 0.140–0.264. Filtering tags 61 of 702 chr1 truth SNVs and 413 of 635 truth INDELs.

*Rules* (renamed "population filtering resources", R3.m4). SNVs: exact allele in gnomAD or CoLoRSdb at population AF ≥ 1e-4, in dbSNP unless flagged somatic, or in the 1000 Genomes panel. INDELs: exact allele in gnomAD or CoLoRSdb at AF ≥ 0.01 when it changes a homopolymer of ≥ 7 bp, ≥ 0.001 otherwise. ClairS-TO's default rule tags more HG008T truth INDELs (87.1% against 61.8%). The INDEL rule was chosen on chromosome 1 calls; [待补: re-selection on validation chromosomes and the resulting chr1 numbers (E17c)]. Calls without a GRCh38 position cannot be filtered this way (R3.4) [待补: graph-native filter for them (E9c)].

**证据（内部）。**

- 资源和样本：`experiments/E8b_literature_facts/README.md`。
  - C9：四个文件和 ClairS-TO 的相同（1000G PoN 2,609,566、dbSNP 60,691,395、CoLoRSdb 49,550,902；gnomAD 总数没本地重数）。
  - P1 / P2：GATK `af-only-gnomad.hg38` 是 2017-08 从 b37 lift 过来的；implied AN 约 40,600，"about 20.3 k genomes"。"gnomAD r2.1" 这个标签不对。
  - P4：dbSNP 138（fileDate 2013-08-06）；P5：1000G PoN 是 Mutect2 在 1000G exome 上的结果（推断），样本表不公开；P6：CoLoRSdb v1.1.0，2024-07-31，1,381 samples；P7：cohort 表，含 127 个 HPRC genomes。
  - P9："HG008 data collection began Fall 2023"；P10：COLO829 等不在 gnomAD v2 里不能用样本表排除（[unverified]）。
  - 4.4 私有 allele 测试（chr20 SNV）：HG008-N 2.5 / 7.9，COLO829BL 1.5 / 6.7，CASTLE normals 0.6–1.8 / 4.0–7.2；HapMap 六个 donor（含 HG005）gnomAD 3.4–15.4，CoLoRSdb 14.5–25.9。
  - 4.5（HG008T v0.2，全部记录）：Pansoma rule SNV 982 / 9,609（10.2 %），INDEL 5,783 / 9,362（61.8 %）；ClairS-TO rule 906（9.4 %）/ 8,155（87.1 %）；dbSNP SAO=2 = 0。
  - 4.5 其余五套（README 11:25 写完之后才跑完）：`/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts/somatic_truth_in_resources.<set>.tsv`（11:54–12:13 写出；`collect.sh` 还没重跑，所以仓库里的 `somatic_truth_in_resources.tsv` 只有 HG008T）。全部记录、不限 BED。`Pansoma_rule` / n，`dbSNP_SAO=2`：
    - COLO829T_SMaHT：SNV 1,606 / 44,005，INDEL 1,070 / 2,059；SAO=2 0 / 0。
    - HCC1395_SEQC2_v1.2.1：SNV 1,102 / 39,447，INDEL 119 / 1,625；0 / 0。
    - HCC1954_CASTLE：SNV 1,956 / 24,033，INDEL 528 / 2,556；0 / 0。
    - H1437_CASTLE：SNV 25,323 / 82,529，INDEL 2,508 / 4,896；0 / 0。
    - H2009_CASTLE：SNV 11,445 / 161,918，INDEL 1,422 / 6,519；0 / 0。
    - 我按计数算的百分比（不在任何表里，Response 没用）：HCC1954 8.1 % / 20.7 %，H1437 30.7 % / 51.2 %，H2009 7.1 % / 21.8 %。COLO829T 和 HCC1395 的计数和下面 crosstab_colo.md / E5b 表 7 的完全一样。
    - H1437 不正常：SNV 里 `gnomAD_AF>=1e-2` 22,982、`1000G_allele` 7,664，远多于其他五套（例 H2009 8,141 / 2,676，n 是 H1437 的两倍）。CASTLE H1437 的 truth 是 DeepSomatic single-cancer 模型做的（E8b D8），可能含 germline。H1437 的数在核实前不进 Response；Response 只用 "六套里没有 dbSNP somatic 标记" 这一条。
    - HapMap_positives_chr20（task 6，387253_6）还 held，没出。
- 其他 truth 被 PoN 标的比例：
  - COLO829T 3.6 % / 52.0 %：`tmp/graph_absorbed_pon_20261012/crosstab_colo.md`（SNV 1,606 / 44,005，INDEL 1,070 / 2,059）。
  - HCC1395T 2.8 % / 7.3 %：`experiments/E5b_absorbed_d46/tables.md` 表 7（SEQC2 v1.2.1，SNV 1,102 / 39,447，INDEL 119 / 1,625）。
- germline 不在图里：`experiments/E7_graph_composition/d9_content.tsv`（bed scope，`in_d9_vcf`）HG008 germline SNV 0.8697、INDEL 0.7513；COLO829 0.8721、0.7457。Response 里的 "13% / 25%" 是 1 − 这些值后取整（HG008 13.0 % / 24.9 %，COLO829 12.8 % / 25.4 %）。E7 没有 README，算初步。
- PoN 标 germline tensor：`analysis/tensor_recall_20260930/README.md` §3（chr1–22，SNV label 2）96.4 / 95.8 / 96.4 / 96.1 / 95.5 / 96.6 %。注意这是 2026-09-30 的 tensor 集（10-03 ch6 重建之前），SNV PoN 规则和现在一样（AF ≥ 1e-4）。graph 上的 germline tensor 就是图没吸收掉的 germline，所以这是 "不在图里的 germline" 的代理，精确拆分要 E4c。
- 最终模型上 PoN 对 germline call 的作用（新）：`experiments/E1.1_ch6_counterfactual/results/summary.json`，base 模式，chr1 test tensors，阈值 = validation recall-0.9，per-tensor PoN（每个 on-reference tensor 的 GRCh38 allele 用 SNV 规则匹配四个 PoN；off-reference tensor 不算，和 rtg 一样）。字段 `germline.ge_thr` → `germline_after_pon_ge_thr`，再除以 `truth_after_pon.fp`：
  - HG008T d9：Illumina 4,864 → 609（FP 688）、PacBio 4,951 → 571（583）、ONT 5,147 → 618（658）。去掉的比例（我算的）87.5 / 88.5 / 88.0 %；剩下占 PoN 后 FP 88.5 / 97.9 / 93.9 %。
  - HG008T d46：2,542 → 654（826）、2,528 → 618（642）、2,810 → 672（771）。去掉 74.3 / 75.6 / 76.1 %；占 FP 79.2 / 96.3 / 87.2 %。
  - COLO829T（没写进 Response）：d46 827 → 191（250）、808 → 172（221）、942 → 191（402）；d9 joint 666 → 60（92）、1,799 → 158（210）、1,396 → 140（399）。ONT 的 PoN 后 FP 里非 germline 占一半以上。
  - Illumina 的数是 class-stratified sample 的加权估计（E1.1 README）。per-tensor PoN 在这些 graph run 上和 rtg 很接近（tables.md §2，例 HG008_PacBio 576 / 583 vs rtg 577 / 584）。
  - 这里的 germline = tensor 标签 2（HG008-N / COLO829BL dipcall germline truth，`analysis/truth_labels.md`），graph 上的 germline tensor 是图没吸收掉的 germline，所以仍然只是 "不在图里的 germline" 的代理；按 allele 是否在 d9 / d46 图里拆要 E4c。
  - d46 去掉的比例低（74–76 % vs 88 %）：d46 已经吸收了更多常见 germline，剩下被叫的 germline 更稀有（推断，没单独验证）。
- HG008T 不加 PoN：`analysis/benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/summary.json`（d9，eval_bed，PASS = validation recall-0.9）：
  - raw P / R / F1：0.1058 / 0.8623 / 0.1884（Illumina）、0.1078 / 0.8565 / 0.1914（PacBio）、0.102 / 0.8551 / 0.1823（ONT）。
  - pon：0.4516 / 0.8307 / 0.5851、0.497 / 0.8278 / 0.6211、0.4656 / 0.8264 / 0.5957。
  - INDEL PASS F1 raw 0.1009 / 0.041 / 0.0114，pon 0.2412 / 0.0891 / 0.0149（PASS = validation recall-0.9 阈值）；best F1（chr1 曲线上选的点，test-tuned 上限，held-out 规则还没定）raw 0.1467 / 0.1503 / 0.0899，pon 0.2637 / 0.2597 / 0.1402。原稿只写 best F1，现在 PASS 和 best 都写。
  - 都是 BED 内、unplaced PASS calls 不计分（d9 SNV unplaced 116 / 211 / 1,409，INDEL 563 / 2,075 / 7,288，`eval_bed.unplaced_pass`；E3 §2）。
  - truth_pon_tagged：SNV 61 / 702，INDEL 413 / 635。`hg008t_snv_indel_7fc7396_20261012/README.md`："after it no caller recalls more than 0.350 of them"。
  - d46 和 linear 的 raw / pon 也在同一个 summary.json 里（例 linear raw F1 0.0095 / 0.0093 / 0.0084）。
- 规则：`hg008t_snv_indel_7fc7396_20261012/README.md` "PoN rules"；dbSNP 只去掉 SAO=2（E8b 4.1）。
- INDEL 规则是在 chr1 上选的：`analysis/tensor_recall_20260930/pon/indel_rules/README.md` §3 Caveats："The call-level columns come from chr1, the test chromosome. The INDEL runs keep no validation predictions (`val_fraction` 0.05, nothing saved), so the floor was compared on chr1."；10-10 的 homopolymer 规则也是在 chr1 calls 上定的（memory `indel-pon-rules`："Caveat: still chosen on chr1"）。
- 和详版的出入：
  - 详版 R2.6 写 "The rule was selected on validation chromosomes and applied unchanged to chr1"，和上面的来源矛盾，Response 已改成 "chosen on chromosome 1 calls"。
  - 详版第 3 点是旧 run 的数（0.18–0.20 / 0.59–0.62），这里换成 7fc7396 benchmark 的数。
  - 详版第 2 点的 "[待补]% of germline tensors whose allele is not in the graph" 还是空的；Response 先写最终模型上 PASS germline call 被 PoN 去掉的比例和剩下 FP 的组成（上面 E1.1 的数）。

**还缺什么。**

- [CPU 本地] E8b：task 1–5 已写出（见证据）；task 6（HapMap，387253_6）还 held，要放出来。跑完后 `bash experiments/E8b_literature_facts/collect.sh`，把七套写进 `somatic_truth_in_resources.tsv` / `.pct.tsv`，Response 再引用表里的百分比。
- [CPU 本地] 核实 H1437 CASTLE truth 为什么有 30 % 左右的 SNV 被 PoN 标（`gnomAD_AF>=1e-2` 22,982 / 82,529）：是不是 truth 里有 germline。这也影响 R1.m1 / R3.m2 里 H1437 的数。
- [CPU 本地] E4c：在最终 tensor 集上把 germline tensor 按 allele 是否在 d9 / d46 图里分两组，各算 PoN 标的比例；PASS germline call 和 PoN 后剩下的 germline FP（E1.1 已有总数）也这样分。按 tensor 做就不需要 call set：`experiments/E1.1_ch6_counterfactual/ch6_counterfactual.py` 已经会算每个 on-reference tensor 的 GRCh38 allele 和 PoN 标记（只是没存进 summary.json），再加 d9 deconstruct VCF / d46 `vg deconstruct` 的 allele 表就能拆。
- [CPU 本地] E9c（在跑：`e9c_runs_d46` 387431、`e9c_collect` 387434）：给没有 GRCh38 坐标的 calls 一个 graph-native PoN，填 Rules 段最后的 [待补]（和 R3.4 共用）。
- [GPU RIS] E17c：每个 INDEL 模型在 validation blocks 上跑 predict 并保存 predictions（现在的 INDEL run 没存），再写 graph / linear VCF。
- [数据同步] E17c 的 validation predictions 和 INDEL 的 held-out 阈值都要从 RIS 同步。
- [CPU 本地] E17c：同步后在 validation calls 上扫 INDEL PoN 规则，选定后在 chr1 只评估一次。
- [用户决定] 如果 E17c 选出的 INDEL 规则和现行规则不同，是否接受换规则并更新所有 INDEL 数字。
- [用户决定] HapMap：从主 benchmark 去掉，还是用已有的去 donor 图 `mgi_non_HG005.d9`（E8b H10）放进 Supplement（和 R1.5 一起定）。
- [数据同步] 用户给已有 ClairS-TO / DeepSomatic 结果的路径，按同一 PoN 重新评分（"the filtering we apply to all callers"）。
- [文字] 资源改名（R3.m4）；"gnomAD r2.1" 改成 GATK af-only 文件的真实来源；详版的 "selected on validation chromosomes" 改掉；详版第 2 点换成上面 E1.1 的数。

### R2.7 Alignment statistics overstated

**状态：部分解决**: 同一批 reads 在不同参考上的计数和比例已经有了：HG008T Illumina 的整份 GAM `vg stats`（d9 vs GRCh38-only）、HG008T PacBio 四个参考、HG008T ONT 三个参考（GRCh38 BAM / GRCh38-only / d9）的 MAPQ 和 unmapped（E11a 完整扫描、`vg stats` + GRCh38 BAM 的 `samtools stats`），以及 node 层面的计数（含按节点长度分档）；GRCh38 BAM 的比对软件已从 BAM header 核实；但承诺的完整表（E11a：两个肿瘤 × 三个平台 × 四个参考，perfect / mismatch / indel / soft clip / MAPQ，按区域分层）还没出：HG008T 的 window scan 已完成、`compare.py` 还没跑，d46 和 COLO829T 的 job 仍 held。

**Comment（要点）。**

- 第 292–294 行说 graph 上 perfect reads 更多、soft clipping 更少，对照 Supplementary Fig. 4 看有些夸大。
- 要给出 read 数或比例。

**Response.**

We agree. Lines 292–294 are rewritten, and Supplementary Fig. 4 is replaced by a table of counts and proportions for the same reads on each reference. Examples:

- *HG008T Illumina* (2,974,694,440 read records, mates counted separately): 2,286,270,813 (76.86%) align full-length without any edit on HPRC v1.1, against 2,009,031,268 (67.54%) on the GRCh38-only graph; 1.42% vs 1.51% are unmapped. Insertion events fall from 49,909,068 to 17,762,998, deletion events from 52,666,034 to 22,692,983, and soft-clip events from 93,004,309 to 85,507,213 (3.17 vs 2.92 per 100 mapped reads). The number of reads at MAPQ 60 is almost the same (2,626,084,493 on HPRC v1.1, 2,628,588,329 on the GRCh38-only graph).
- *Soft clipping on ONT* does not change: COLO829T ONT has 173.45 (GRCh38-only) and 173.44 (HPRC v1.1) soft-clip events per 100 mapped reads.
- *Long reads* are rarely edit-free end to end (HG008T PacBio: 1.39% of reads on HPRC v1.1), so we also count node visits. HG008T PacBio has 3.93 billion node visits by MAPQ > 5 alignments on HPRC v1.1; 8.0% are on branch (non-reference) nodes, and 98.0% of these are edit-free, against 92.4% on GRCh38 nodes. Most branch nodes are 1 bp, and within the same node length the two are about equally often edit-free (1 bp: 98.5% vs 97.9%).
- *MAPQ*, for the same 14,014,144 HG008T PacBio reads: unmapped 0.01% (GRCh38 BAM), 2.78% (GRCh38-only graph), 2.64% (HPRC v1.1) and 0.01% (HPRC v2.1); MAPQ 60 for 94.84%, 95.88%, 95.88% and 83.32% of mapped reads. For the same 3,642,157 HG008T ONT reads, giraffe leaves more reads unmapped than minimap2: 3.89% (GRCh38 BAM), 7.88% (GRCh38-only graph) and 7.77% (HPRC v1.1); MAPQ 60 for 92.44%, 95.18% and 95.02% of mapped reads.

The table states the aligner of each arm and that MAPQ scales are aligner-specific: the GRCh38 BAMs come from BWA-MEM 0.7.17 (Illumina), pbmm2 1.10.0 (PacBio HiFi) and minimap2 2.26 (ONT) for HG008T, and from pbmm2 1.13.0 (fiberseq) and minimap2 2.26 (ONT) for COLO829T; the graph arms use vg giraffe 1.66 (GRCh38-only) and 1.77 (HPRC v2.1) [待核实: giraffe version of the v1.1 alignments]. [待补: the complete same-read table for HG008T and COLO829T × Illumina / PacBio / ONT on all four references (unmapped, perfect, mismatch, indel, soft clip, MAPQ), also per SMaHT region class (E11a)]

We also toned down the statement. The graph mainly removes edits at alleles that the sample shares with the HPRC haplotypes; soft clipping changes little, and more edit-free alignments alone do not show that placements are more correct (R3.6).

**证据（内部）。**

- HG008T Illumina d9：`/scratch/jshen/data/HG008_GIAB/vg_giraffe_GAM/Liss_lab_BCM_Illumina-WGS_20240313/HG008-T_p23_BCM_IlluminaWGS_20240313.vg.stats`：Total alignments 2974694440，Total aligned 2932326192，Total perfect 2286270813，`max 60 (2626084493 reads)`，Insertions 17762998 read events，Deletions 22692983 read events，Softclips 85507213 read events。
- HG008T Illumina GRCh38-only：`/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_BCM_Illumina-WGS_20240313/vg_linear_hg38/HG008T_BCM_IlluminaWGS_hg38_linear_vg.gam.stats.txt`：Total alignments 2974694440，aligned 2929783717，perfect 2009031268，`max 60 (2628588329 reads)`，Insertions 49909068、Deletions 52666034、Softclips 93004309 read events。
- `vg stats` 的 "read events" 是 edit 事件数，不是 read 数（例：COLO829T ONT d9 "Substitutions: 2107426682 bp in 1674490568 read events"，Illumina 的 substitution bp 与 events 相等），所以 Response 写 "insertion / deletion / soft-clip events"，不写 "reads with an insertion"（checker 改）。
- 百分比（76.86 / 67.54 perfect，1.42 / 1.51 unmapped，2.92 / 3.17 和 COLO829T ONT 173.44 / 173.45 soft-clip events per 100 mapped reads，PacBio d9 perfect 1.39）：`experiments/E11a_alignment_stats/results/vg_stats_crosscheck.md`。
- Node 层面：`analysis/branch_node_alignments_20261002/README.md` s1 表，HG008T_PacBio：node visits 3,931,468,159，branch 313,140,918 (8.0 %)，perfect on branch 98.0 %，GRCh38 nodes 92.4 %。按长度分档（同 README "Perfect fraction of visits by node length" 表）：HG008T_PacBio 1 bp branch 98.5 %（267 M visits）vs GRCh38 97.9 %；README 原话 "a visit's chance of an edit grows with node length, so compare within a length bin"。98.0 vs 92.4 主要是节点长度造成的（branch visit 大多在 1 bp 节点上），只引这一对数字正是审稿人说的那种夸大，Response 已加分档的对比（checker 改）。
- PacBio MAPQ，三个 graph 臂：`/scratch/jshen/data/reviewer_experiments_20261013/E11a_alignment_stats/scans/HG008T_PacBio.{linear,d9,d46}.summary.json`（`complete: true`，job 387166_0–2 COMPLETED）：with_path 13,624,706 / 13,644,710 / 14,012,552（records 14,014,144），`mapq_hist_with_path[60]` 13,064,032 / 13,081,945 / 11,675,897。
- PacBio MAPQ，GRCh38 BAM：`/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008_PacBio_bam_stats.txt`（samtools stats 1.18）：raw total sequences 14014144，reads unmapped 1055，`MAPQ 60` 行 13,289,473，primary mapped 14,013,089。
- 94.84 / 95.88 / 95.88 / 83.32 % 和 0.01 / 2.78 / 2.64 / 0.01 % 是我用上面两条的计数算的（两位小数）。
- HG008T ONT（checker 加进 Response 的 MAPQ 条）：3.89 % = 141,830 / 3,642,157；7.88 % = (3,642,157 − 3,355,051) / 3,642,157；7.77 % = (3,642,157 − 3,359,067) / 3,642,157（与 `vg_stats_crosscheck.md` 的 7.77 一致）；MAPQ 60 的 95.02 % = 3,191,911 / 3,359,067。来源同下一行：
- 原始计数：GRCh38 BAM `HG008_ONT_bam_stats.txt` unmapped 141,830 / 3,642,157，MAPQ 60 3,235,858 / 3,500,327（92.44 %）；GRCh38-only graph `scans/HG008T_ONT.linear.summary.json` with_path 3,355,051，MAPQ 60 3,193,314（95.18 %）；d9 `HG008_ONT_gam_stats.txt` aligned 3,359,067，`max 60 (3191911 reads)`。d46 ONT scan 还没跑完。
- vg 版本：`experiments/E3_reference_arms/README.md` 表和 caveat 4（linear vg giraffe 1.66，d9 unverified，d46 vg giraffe 1.77）；d9 GAM 是 2026-03 比对的（目录时间戳），版本仍未查到：唯一的 log `/scratch/jshen/data/HG008_GIAB/vg_giraffe_GAM/Liss_lab_BCM_Illumina-WGS_20240313/log/HG008-T_p23_BCM_IlluminaWGS_20240313.log` 只有四行 "Guessing that ... is ..."（index 类型），没有版本号。
- GRCh38 BAM 的比对软件（新，`samtools view -H` 的 `@PG` 行，samtools 用 `/wanglab/jshen/anaconda3/envs/sambcfenv/bin/samtools`）：
  - HG008T Illumina `HG008-T_Illumina_195x_GRCh38-GIABv3.bam`：`PN:bwa VN:0.7.17-r1188 CL:bwa mem -t 20 ...`；
  - HG008T PacBio `HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam`：`PN:pbmm2 VN:1.10.0 (commit v1.10.0) CL:pbmm2 align ... --preset HIFI ... --unmapped`；
  - HG008T ONT `HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.bam`：`PN:minimap2 VN:2.26-r1175 CL:minimap2 -t 25 -y -x map-ont ...`；
  - COLO829T fiberseq `SMHTCOLO829T-X-X-M45-B004-uwsc-MERGED-pbmm2_1.13.0_GRCh38.aligned.sorted.bam`：`PN:pbmm2 VN:1.13.0 (commit v1.13.0-2-gbcd99f5)`；
  - COLO829T ONT `COLO829T_ONT_merged.bam`：`PN:minimap2 VN:2.26-r1175`（三条 `ID:minimap2*` 都是 2.26；header 里另有一条 `ID:aligner PN:minimap2 VN:2.24-r1122`，和 basecaller `dorado 0.4.2+bb2a133` / `guppy 6.5.7+ca6d6af` 并列，来源不明，可能是上游流程；Response 暂写 2.26，见还缺什么）。
  - 三个 HG008T BAM 的参考都是 GIAB v3 的 `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18`（bwa / pbmm2 的 CL 里可见）。COLO829T Illumina 不在 E11a 的 GRCh38 BAM 臂里，没查。
- E11a 的新输出（12:00 之后核对，都还没进表）：`/scratch/jshen/data/reviewer_experiments_20261013/E11a_alignment_stats/windows/` 有 HG008T_Illumina.{linear,d9,grch38_bam}、HG008T_PacBio.grch38_bam、HG008T_ONT.grch38_bam 的 `.npz` + `.summary.json`（`complete: true`；sacct：387287_0、_1、_3、_4、_5 COMPLETED）。`summary.json` 只有 scan 元数据（窗口数 48、walked reads、耗时），比例要 `compare.py` 算；`results/` 里仍只有 `vg_stats_crosscheck.md`，README 仍是占位符。

**还缺什么。**

- [CPU 本地] E11a：`compare.py` 对不全的 set 也能出表（`analysis_g` / `analysis_w` 各需 ≥ 2 个完整 arm），所以现在就能先跑 HG008T_PacBio（G：三个 graph scan；W：再加 GRCh38 BAM windows）、HG008T_Illumina（W：GRCh38 BAM / linear / d9）和 HG008T_ONT（W：GRCh38 BAM / linear），得到同一批 reads 的 unmapped / perfect / mismatch / indel / soft clip / MAPQ 表和 SMaHT easy / difficult / extreme、segdup 分层；仍 held 的是 387166_[5-6]（HG008T ONT d9 / d46 scan）、387287_[2,6-10]（HG008T Illumina d46 windows、COLO829T windows）、387288_[12-…]（COLO829T scans），放出后再跑一次补全 d46 和 COLO829T。
- [文字] 改写第 292–294 行；把 Supplementary Fig. 4 换成计数表；[待核实] Supplementary Fig. 4 原来画的是哪些数据集、哪个参考。
- [文字] 表注写上各臂的比对软件和版本（GRCh38 BAM 已核实，见上）；d9 GAM 的 giraffe 版本仍未核实（E11b README 认为是 `vg_v1.65.0`，E3 标为 unverified），查不到就在表注写 "version not recorded"；COLO829T ONT BAM 的 header 同时有 minimap2 2.24 和 2.26，需向数据来源（SMaHT）确认最终比对用的是哪个。

### R2.m Minor comments 1–8

**状态：部分解决**: 第 1–4、7、8 条只是改图和文字，第 8 条的定义已和 v4 README 核对；但第 5 条 "truth overlapping" 在原稿第 13 页的意思还没核实（原稿不在 repo 里），第 6 条要 E18 的 SBS96 / truth 谱，E18 还没开始：比较器不重跑、等用户给已有 VCF 的路径，Pansoma 的 chr1 call sets 本地只有三个 d46 SNV run。

**Comment（要点）。**（原审稿文字不在 repo 里；下面按详版的逐条 Response 反推。）

1. Figure 2d、2f 画的是 truth-set variants 还是 Pansoma 的输出，要说清楚。
2. 图的引用顺序（2e 在 2f 之前引用）。
3. Supplementary Fig. 2 的 COLO829T panel 没有标注。
4. Supplementary Fig. 1 的点标签压在点上。
5. 第 13 页的 "truth overlapping" 没有定义。
6. Figure 4 里把 truth-set 的 spectrum 和各 caller 放在一起（同 R3.5）。
7. Figure 4 的 panel 引用顺序和正文不一致。
8. 定义 candidate 和 candidate site。

**Response.**

We thank the reviewer and made the following changes.

1. Figures 2d and 2f show truth-set variants, not Pansoma output. The legends now say so, and both panels moved to the
   Supplement.
2. Figures are now cited in order (2e before 2f).
3. The COLO829T panel of Supplementary Fig. 2 is labelled.
4. Point labels in Supplementary Fig. 1 were moved off the points.
5. "Truth overlapping" was replaced by defined categories, for example a call within 10 bp of a truth allele but with a
   different allele. [待核实: the meaning of the original sentence on page 13]
6. The truth-set spectrum is now shown beside each caller's calls in Figure 4 (R3.5). [待补: E18 panels from the final
   call sets]
7. Figure 4 panels are cited in the order of the text (see item 2).
8. A *candidate* is one allele at a graph position: node, start on the node's forward strand, REF, ALT and type (SNV,
   insertion or deletion). A *candidate site* is all candidates that pass the filters at one node position and variant
   class (SNV or INDEL); they are encoded together in one tensor, each read labelled with the allele it carries (or as
   REF or OTHER).

**证据（内部）。**

- 第 8 条：`indexed_gam_pipeline_v4/README.md` §8 "Candidate identity is `(node, forward start, REF, ALT, kind)` with
  `kind ∈ {SNP, INS, DEL}`"；"Sites — surviving alleles are grouped by `(node, start, SNV|INDEL)` ... all of them are in the
  site's one tensor: every record covering any passing allele is a row, labeled with the allele it carries"。
- 第 5 条用的 10 bp：`analysis/truth_labels.md` 常量表 `NEAR_BP` = 10（"附近有 truth" 的距离，R11）。原稿第 13 页的上下文看不到，
  所以这个对应关系是推测。
- 第 6 条：
  - `points/R3.5.md` 状态还是 "待解决"；`experiments/` 下还没有 E18 目录（2026-10-13 15:40 查过）。
  - 比较器：用户决定不重跑 ClairS-TO / DeepSomatic（用户已有结果，re-score 等路径）；`experiments/E8_comparators/README.md`
    还写着 "PENDING"，已过时。
  - Pansoma 这边，本地现在有 chr1 call sets 的只有三个 d46 SNV run 的 GPU re-prediction（`HG008_PacBio_SNV_d46_7fc7396`、
    `HG008_ONT_SNV_d46_7fc7396`、`COLO829_ONT_SNV_d46_7fc7396`，在
    `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1/`，有 `*.linear.vcf.gz`、
    `*.linear.unplaced.vcf.gz` 和 rtg 输出）；它们复现了 RIS（`eval_bed/report.txt` 的 pon F1 0.603 / 0.575 / 0.736，RIS
    0.605 / 0.575 / 0.737，`experiments/E3_reference_arms/README.md` 的 d46 列）。同目录的 `xs__*` 是跨样本预测（E12a），不是 Figure 4 要的 in-sample call sets。
- 第 8 条措辞：README §8 "Sites" 只把 "surviving alleles"（过了 ALT / AF 等过滤的）放进 site，读不带任何 allele 的标成 `REF` /
  `OTHER`，所以 Response 补了 "that pass the filters" 和 "(or as REF or OTHER)"。
- 第 1 条（2d / 2f 画的是 truth-set variants）出自详版的判断；原稿不在 repo 里，本地核不了，改图时顺手确认。
- 其余各条只涉及图和文字，没有数字。

**还缺什么。**

- [用户决定] 第 5 条：用户确认原稿第 13 页 "truth overlapping" 指的是什么（labels 的 "near truth"，还是评估里 FP 的某一类），
  再定替换成哪个定义。
- [数据同步] 第 6 条：用户已有的 ClairS-TO / DeepSomatic VCF 路径（比较器不重跑）；Pansoma 其余 run 的最终 chr1 call sets
  从 RIS 同步（本地只有上面三个 d46 SNV run）。
- [CPU 本地] 第 6 条：E18 SBS96 / truth 谱（和 R3.5 共用）。
- [文字] 第 1–4、7 条：改 Figure 2、Figure 4、Supplementary Fig. 1 和 2；第 8 条的定义写进 Methods 的术语段（和 R3.m4 统一）。

---

## 6.3 Reviewer 3

### R3.1 Reference-level ablation

**状态：部分解决**: 三个参考臂（GRCh38-only / HPRC v1.1 d9 / HPRC v2.1 d46）× 六套数据 × SNV / INDEL 的模型和 chr1 结果都已出（E3）：SNV 在 PoN 后（BED 内、unplaced 不计分）、INDEL 按 test-tuned best F1，HPRC 在全部 6 套上都高于 linear；no BED + unplaced 全算 FP 的 worst case 下，HG008T ONT（两张图）、COLO829T fiberseq（d9）和 COLO829T ONT（d46）的 SNV 不再高于 linear。新加的 E12a 给了 d46 臂的跨肿瘤初步结果，E1.1 确认 linear 臂的 ch6 不起作用。还缺：去 ch6 臂、matched d9 重训、seed / bootstrap CI、linear 和 d9 臂的跨肿瘤检查、INDEL 的 held-out 阈值、E9 实测的 complete-output 数字扩到其余 run（三个 d46 run 已完成，见 R3.4），以及 COLO829T d9 的 vg 版本核实。

**Comment（要点）。**

- End-to-end benchmarking does not isolate the pangenome: Pansoma differs from the comparators in reference, alignment, candidate generation, tensorisation, architecture and post-processing.
- Run the Pansoma pipeline with a GRCh38 linear reference (or a single-path graph) and with the full HPRC graph. Keep candidate extraction, tensors, architecture, training and filtering thresholds identical, and show how much of the gain the pangenome drives.

**Response.**

We performed this ablation. The same reads of six datasets (HG008T and COLO829T, three platforms each) were aligned with vg Giraffe to three references:

- a GRCh38-only graph built by vg autoindex from the same GRCh38 FASTA (its haplotype-count channel is constant, R1.1);
- HPRC v1.1 (90 haplotypes);
- HPRC v2.1 (464 haplotypes).

Candidate extraction, tensors, labels, truth, network, training recipe, threshold rule (validation recall 0.9) and population filtering use the same code. One SNV and one indel model per dataset and reference were trained on chr2–22 and tested on chr1. The GRCh38-only and v2.1 runs are fully matched; eight of the twelve v1.1 models differ in training details (GPU count, code version or indel learning rate). [待补: matched v1.1 retraining]

*SNVs.* After population filtering (chr1, in the truth BED), both HPRC graphs beat GRCh38-only on all six datasets: F1 0.418–0.680, against 0.585–0.781 (v1.1) and 0.570–0.783 (v2.1), a gain of +0.038 to +0.319. Without population alleles in the graph, most germline variants become candidates (2.2–7.6 times as many chr1 tensors as on v1.1) and raw F1 falls to 0.008–0.075; 1.4–7.8 times as many false positives survive filtering as on v1.1. The two HPRC graphs are equivalent after filtering (−0.021 to +0.031; R1.2).

*Indels.* The recall-0.9 threshold lands far below the best operating point, so we compare each chr1 curve's best F1, a test-tuned upper bound. [待补: validation-chosen threshold] GRCh38-only is lower on all six datasets: 0.019–0.199, against 0.071–0.264 (v1.1) and 0.085–0.269 (v2.1); it keeps more somatic indels as candidates but ranks them worse. At the recall-0.9 threshold it scores higher on the four long-read sets (HG008T PacBio 0.125 against 0.089; all arms are poor there), because its threshold lands at fewer false positives.

*Calls without a GRCh38 position.* If calls outside the truth BED are scored and every unplaced call counts as a false positive (a worst case), the SNV advantage disappears on HG008T ONT (0.342 and 0.354 against 0.365), COLO829T Fiber-seq with v1.1 (0.636 against 0.680) and COLO829T ONT with v2.1 (0.369 against 0.418), and holds elsewhere. R3.4 describes how these calls are resolved. [待补: measured numbers, E9]

*Controls.* Giraffe 1.66.0 aligned the HG008T reads to both GRCh38-only and v1.1, so that contrast (+0.059 to +0.116) has no aligner-version effect; v2.1 used 1.77. [待核实: version of the COLO829T v1.1 alignments] Without retraining (preliminary, v2.1, in the BED, unplaced calls not scored), the HG008T SNV models score 0.708, 0.754 and 0.630 on COLO829T, above COLO829T's own GRCh38-only models (0.468, 0.680, 0.418). The COLO829T models transfer worse at their own threshold, falling below HG008T's GRCh38-only models on Illumina and ONT (0.359 and 0.487 against 0.469 and 0.537), though not at the best point (R1.7). The GRCh38-only arm is Pansoma on a linear graph; R3.m2 compares dedicated linear callers. [待补: cross-tumor check of the GRCh38-only and v1.1 arms; GRCh38-only realigned with 1.77; HPRC arm without channel 6; repeated seeds and block-bootstrap intervals (R3.9)]

**证据（内部）。**

- 三臂表：`experiments/E3_reference_arms/tables.md` §1–§3，数据来自 `analysis/benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/summary.json`。
  - SNV PoN 后 PASS F1（BED 内）：linear 0.469 / 0.513 / 0.537 / 0.468 / 0.680 / 0.418，d9 0.585 / 0.621 / 0.596 / 0.778 / 0.781 / 0.706，d46 0.570 / 0.605 / 0.575 / 0.783 / 0.779 / 0.737。
  - d9 − linear 是 +0.059 到 +0.311，d46 − linear 是 +0.038 到 +0.319。按 test-tuned best F1 比，linear 也在 6 套上都低。
- raw PASS F1：linear 0.008–0.075，d9 0.182–0.571，d46 0.287–0.685（linear README 结果表）。
- chr1 tensor 是 d9 的 2.2–7.6 倍（HG008T PacBio 251,303 对 32,880；E3 §3 逐行算：414,453 / 190,697 = 2.17，251,303 / 32,880 = 7.64）。PoN 后的 FP：HG008T 是 d9 的 1.4–1.7 倍，COLO829T 1.8–7.8 倍（Illumina 2,216 对 284）。见 `linear_grch38_7fc7396_20261012/README.md`。
  - linear README 和 E3 README 写 "every germline variant of the sample becomes a candidate"，说过头了：candidate 表（R1.2）里 linear 上 matched-normal germline 有 candidate 的是 0.757–0.822。Response 改成 "most germline variants"。
- 训练设置（E3 `tables.md` §4、E3 README Caveats 1、d46 README "The d9 runs compared"）：linear 和 d46 完全一致（SNV 4 GPU，INDEL 1 GPU，lr 1e-4）。d9 里只有 HG008T Illumina SNV 和 HG008T 三套 INDEL 一致；HG008T PacBio / ONT SNV 是 1 GPU，COLO829T SNV × 3 是 9d6bd8d 的 `*_SNV_h90`（1 GPU），COLO829T INDEL × 3 是 45d56ad lr 4e-4。共 8 / 12 不一致。早上的 Response 写 "Everything else is identical ... training recipe"，对 d9 臂不成立，已改。
- INDEL best F1（PoN 后）：linear 0.199 / 0.173 / 0.082 / 0.069 / 0.092 / 0.019，d9 0.264 / 0.260 / 0.140 / 0.195 / 0.120 / 0.071，d46 0.249 / 0.269 / 0.155 / 0.163 / 0.129 / 0.085。
  - PASS F1：HG008T PacBio 0.125 / 0.089 / 0.095，ONT 0.022 / 0.015 / 0.015，COLO829T fiberseq 0.040 / 0.016 / 0.020，ONT 0.004 / 0.003 / 0.002。
  - 在 best point 上，linear 的 FP 是 d9 的 1.6–2.9 倍（HG008T），但 TP 不多（PacBio 97 TP / 437 FP，d9 是 120 / 217）。
  - HG008T 有 tensor 的 truth INDEL：linear 0.75–0.87，HPRC 0.54–0.76（E3 README）。
  - 来源不一致：E3 README 按 PASS F1 写的是 "no pangenome advantage after the PoN, except on Illumina"。按较新的 benchmark README（commit 0049e36，history 改写前是 d4f898c）写：主比较用 best F1，PASS 的差异是阈值落点造成的。best F1 是 test-tuned 上界，held-out 规则还没定。
- Complete-output bound（no BED，unplaced 全算 FP）：E3 `tables.md` §2。早上的 Response 只写了 "every unplaced call counts as a false positive"，没说 no BED；已补 "calls outside the truth BED are scored"。HG008T 上 no BED 本身就让 linear 掉很多（placed-only 0.197 / 0.359 / 0.365，BED 内 0.469 / 0.513 / 0.537）；COLO829T 上 no BED 和 BED 内的 placed-only F1 一样（linear README 结果表）。
  - HG008T ONT：d9 0.342、d46 0.354，linear 0.365。
  - COLO829T ONT：d9 0.522、d46 0.369，linear 0.418（d46 有 4,566 个 unplaced PASS，TP 是 1,690）。
  - COLO829T fiberseq：d9 0.636（1,018 个 unplaced）、d46 0.736，linear 0.680。
  - HG008T Illumina（0.546 / 0.414 对 0.197）、HG008T PacBio（0.553 / 0.533 对 0.359）、COLO829T Illumina（0.711 / 0.702 对 0.468）上，两个 HPRC 图都在 linear 之上。
  - 来源不一致：E3 README 写 "SNV d9 / d46 stay above linear on Illumina, PacBio and fiberseq"，和它自己的 `tables.md` §2 对不上（COLO829T fiberseq 的 d9 是 0.636 < 0.680）。这里按 `tables.md` 写。
  - E9a / E9c / E9d 正在三个 d46 run 上跑（HG008T PacBio、HG008T ONT、COLO829T ONT）。它们的输入是今天在 tequila GPU 上重新 predict 的结果（`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1/eval_bed/vcfeval/pon_calls/summary.txt`，"None" 行）：rtg PoN 后 F1 0.6034 / 0.5749 / 0.7364，RIS 是 0.605 / 0.575 / 0.737，复现了。
- linear 臂的 ch6（E1.1，新）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §3–§4。6 个 linear SNV 模型在 8 个 mode 下的 AP、germline ≥ t、truth TP、PoN 后 F1 全部和 base 一样（ch6 恒为 100）。所以 linear 臂确实是 "没有 haplotype 信息" 的对照。HPRC 模型则明显依赖 ch6（R1.1）。
- 不要用 E1.1 的 per-tensor "after PoN" F1 做臂间比较：
  - E1.1 §2 里，图模型的 per-tensor PoN 和 rtg 基本一致；linear 上 rtg 的 FP 多很多（COLO829T Illumina 211 对 2,216，ONT 1,063 对 2,575，HG008T 837 / 713 / 692 对 1,091 / 973 / 904）。
  - 可能的原因（未核实）：per-tensor 评估区域是 somatic BED ∩ germline BED（`pansoma_net_v2/data.py` 的 `eval_label`），rtg 只用 somatic BED。linear 上 germline BED 外的 germline candidate 在 per-tensor 里被排除，在 rtg 里算 FP。
  - 臂间比较一律用 rtg（E3）。如果审稿人问 "linear 是否被 VCF 转换拖累"，需要先把这一点核实。
- 跨肿瘤（E12a，新，初步）：`experiments/E12a_cross_sample/README.md` 和 `tables.md`。只有 d46 臂，旧 checkpoint（BQ / MAPQ 用训练肿瘤的统计量归一化）。阈值是模型自己的 validation recall-0.9 阈值，不在目标上重调。
  - HG008T 模型 → COLO829T：0.708 / 0.754 / 0.630（Illumina / PacBio→fiberseq / ONT；best 0.717 / 0.764 / 0.639）。COLO829T 自己的 d46 模型 0.783 / 0.779 / 0.737，COLO829T 的 linear in-sample（E3 §1）0.468 / 0.680 / 0.418。
  - COLO829T 模型 → HG008T：0.359 / 0.557 / 0.487（best 0.554 / 0.597 / 0.579）。HG008T 的 linear in-sample 是 0.469 / 0.513 / 0.537（best 0.544 / 0.596 / 0.570）。所以这个方向在 PASS 上 Illumina 和 ONT 低于 linear in-sample，best F1 则都略高（PacBio 0.597 对 0.596 基本打平）。E12a README 的解释是阈值校准（ranking 能迁移，p_somatic 的校准不能），per-dataset normalization（commit f2a0e4e，改写前 cbc94a7）就是为这个改的。早上的 Response 只写 "transfer less well"，只报对我们有利的方向；现在两个方向都写了。
  - 两个方向都是 BED 内、unplaced 不计分。HG008T 模型 → COLO829T 的 unplaced PASS 是 225 / 550 / 1,388（TP 1,563 / 1,713 / 1,440），ONT 这一格和 E3 §2 里 complete-output bound 反转的情况同类，所以 Response 里加了 "in the BED, unplaced calls not scored"。
  - 目标都是 P100 set（COLO829T_*_d46、HG008T_*_d46、HCC1395T_*_d46），不涉及有 BQ artefact 的 dilution。
- vg 版本，这次核实了：
  - HG008T d9 三个 GAM 都是 vg 1.66.0 "Navetta"：
    - PacBio：`/scratch/jshen/Log_info/vg_longread_gam_AFHPRC.279410.out`，GAM `...AF-HPRC.gam`，2026-04-16。
    - ONT：`hg008t_ont_align_AFHPRC.279553.out`，2026-04-20。
    - Illumina BCM：`hg008t_giraffe_Liss_lab_BCM_Illumina-WGS_20240313_AF-HPRC.268637.out`，GAM `HG008-T_p23_BCM_IlluminaWGS_20240313.gam`，2026-03-04。
  - linear：vg 1.66.0（`/scratch/jshen/Scripts/vg_autoindex_hg38_linear.sh`、`hg008t_pacbio_ont_vg_linear_hg38_gam.sh`；memory 里记着 vg 1.77 读不了 linear 的 longread minimizer）。HG008T PacBio / ONT linear 的日志 `/scratch/jshen/Log_info/HG008T_longread_vg_linear_gam_379046_{0,1}.out` 写 `vg version v1.66.0 "Navetta"`。HG008T Illumina linear GAM（2026-05-08）没找到 giraffe 日志，版本来自脚本 `hg008t_illumina_wgs_vg_giraffe_surject_bam.sh`（`VG=.../vg_v1.66.0`，输出目录就是这个 GAM 的目录）；之后的 gamsort 用的是 vg 1.77（`gamsort_HG008T_Illumina_linear-379040.out`），排序不改比对。
  - d46：vg 1.77.0（`/scratch/jshen/Scripts/hprc_v21_d46_vg177_giraffe.sh`，dist / min 也是 1.77 重建的）。
  - COLO829T d9 GAM 是 2025 年做的（`COLO829T.gam` 2025-05-09，`COLO829T_fiberseq.gam` 2025-06-11，`COLO829T_ONT_merged.gam` 2025-09-04），日志里找不到 vg 版本，仍待核实。
  - 所以 HG008T 上 linear 和 d9 是同一个 aligner 版本，这一对的差值是 +0.116 / +0.108 / +0.059。其中只有 Illumina 是训练设置完全对齐的一对（都是 4 GPU）；PacBio / ONT 的 d9 是 1 GPU 训的。
- d9 可比性：E3 `tables.md` §4、E3 README Caveats 1。只有 HG008T Illumina SNV 和 HG008T 三套 INDEL 是同 commit、同设置。linear 和 d46 之间完全对齐。
- 小错误：`hprc_v21_d46_7fc7396_20261012/README.md` 的 Data 段写 d46 是 "`--haplotypes 90`"，实际是 464（`pansoma_v2_tensors/HG008T_PacBio_d46/tensors/SNV/manifest.json`、`run/config.json`）。E3 README 写的 464 是对的。
- 相同 reads：E3 README（"read counts identical across arms"）。linear 每个节点 H = 17（GRCh38 + 16 条 path cover），ch6 恒为 100（pipeline README channel 表）。

**还缺什么。**

- [用户决定] vg 版本对照要不要真的做（E3 step 6：vg 1.77 重建 linear 索引，重新比对 HG008T PacBio + Illumina，再 build tensor 和训练）。HG008T 的 linear 和 d9 已经核实是同一版本（1.66.0），所以也可以只在文字里说明，不再重比对。
- [CPU 本地] 核实 COLO829T 三个 d9 GAM 的 vg 版本（2025 年的 giraffe 日志或脚本）。
- [GPU RIS] 去 ch6 臂：`E1.1_ch6_counterfactual/gpu_nopc_commands.sh`，27 runs。
- [GPU RIS] matched d9 重训：HG008T PacBio / ONT SNV 用 4 GPU，COLO829T d9 SNV × 3，COLO829T d9 INDEL × 3 用 lr 1e-4（这 3 个和 nopc 脚本里的对照是同一批）。
- [GPU RIS] seed 重复：每臂至少 2 个额外 seed（E16 / R3.9）。
- [数据同步] + [CPU 本地] E16 chr1 block bootstrap CI（臂间配对差值）。要等 RIS 的 call sets 同步回本地。
- [GPU RIS] 跨肿瘤检查的正式版（E3 step 4 / E12）：
  - 用新的默认 `--norm dataset`（commit f2a0e4e）在 RIS 上训练，三个臂都做（d46 的初步结果见上；linear 和 d9 还没有，`tasks_full_plan.tsv` 里也没有 linear 臂的 xs）。
  - 跨肿瘤 predict 本身可以在 tequila GPU 上跑，但要等 RIS 的新 checkpoint。
- [CPU 本地] E9a / E9c / E9d（在跑：E9a array 387435 的三个 task 15:39–15:50 已结束，结果还没 collect；`e9c_runs_d46` 387431_2 在跑；`e9c_collect` 387434、`e9d_asm_*` 387432 / 387433 等 dependency；2026-10-13 下午 sacct / squeue）。完成后用实测数字替换 complete-output 的 worst-case bound 和 Response 里的 [待补: measured numbers, E9]。
- [CPU 本地] 核实 linear 臂 per-tensor PoN 和 rtg 的 FP 差异（germline BED 假设）。要 linear 臂的 rtg FP VCF（目前只在 RIS 上，[数据同步]）。
- [数据同步] INDEL 的 held-out 阈值（validation 上选阈值）要把 RIS 的 predictions 同步回来。之后替换 Response 里 "upper bound" 的 [待补]。
- [用户决定] INDEL 主比较用 best F1（上界）还是等 held-out 阈值。ONT 的 complete-output bound 怎么写（详版 R3.1 判定表）。
- [文字] 简版 R3.1 里的 "[待核实: v1.1]"：HG008T 改成 vg 1.66.0，COLO829T 仍待核实。

### R3.2 Composition of the removed candidates

**状态：部分解决**: candidate 层面的三类对比（somatic truth、matched-normal germline、background）在六套数据 × 三个参考臂上都做完了，E1.1 补了 classifier 层面（ch6 也在压 germline）。新完成的 E5b 回答了 "somatic 有没有被保住"（tensor 层面，chr1-22）：比 linear 少掉的 somatic allele 里，ALT 被图拼出来的（absorbed）95.5–100 % 本来就会被 PoN 标掉；只算 PoN 不标的 allele，图比 linear 少的 candidate 是 SNV 0.2–0.6 %、INDEL 2.0–8.2 %（占全部 truth，净差），其中 absorbed 只占一小部分（chr1 测试区域里每个 dataset × 图丢 10–31 个 PoN 不标的 INDEL candidate，absorbed 0–2 个）。E9c / E9d 给了 d46 臂 FP 的拆分（只有三个 d46 长读长 SNV run，chr1，每组抽 150）：PoN 后剩下的 placed FP 按 normal assembly 146 / 135 / 110 of 150 是 germline，graph PoN（HPRC v1.1 haplotype 加 d46 图）只标 2.0–10.0 %，也就是图里没有的 germline；unplaced PASS 也大多是 population / germline。还缺：linear 臂 FP 的同样拆分（E4b，要 call sets），这是把 79 % / 73 % 和精度提升连起来的直接证据；丢掉的 INDEL 在 linear 上 PASS 能不能叫出来（E5b `linear_recovery.md`，要 RIS 的 rtg 输出或 tequila re-predict）；去掉 ch6 的重训；rescue 默认值（用户决定）。

**Comment（要点）。**

- Fewer candidates help only if germline variants and artefacts are removed preferentially while true somatic variants are kept.
- Compare linear and graph candidate generation for (1) somatic truth, (2) matched-normal germline variants, (3) background / non-truth candidates.
- Stratify the reduction of false positives by whether the germline allele is encoded in the HPRC graph, to connect the 79 % / 73 % overlap of normal-sample SNVs / INDELs with the graph to the precision gain.

**Response.**

We did this on identical reads and code for all six datasets (HG008T and COLO829T, three platforms each) with a GRCh38-only graph, HPRC v1.1 and HPRC v2.1. Fraction of truth alleles that produce a candidate tensor (range over the six datasets; per-dataset values in R1.2):

| | GRCh38-only | HPRC v1.1 | HPRC v2.1 |
|---|---|---|---|
| Somatic SNVs | 0.896–0.986 | 0.893–0.952 | 0.892–0.956 |
| Somatic insertions | 0.529–0.962 | 0.356–0.535 | 0.361–0.578 |
| Somatic deletions | 0.629–0.965 | 0.504–0.718 | 0.524–0.729 |
| Matched-normal germline variants | 0.757–0.822 | 0.077–0.131 | 0.072–0.124 |

Background candidates:

- Non-truth SNV candidates at GRCh38 positions inside the benchmark regions fall by 11–61 % on HPRC v1.1 and 5–59 % on v2.1; total candidates by 33–69 % and 34–70 %.
- Non-truth INDEL candidates fall by 75–84 % for Illumina but change by −0.8 % to +3.9 % for long reads; many look like low-AF sequencing errors, which a population graph does not contain (COLO829T Fiber-seq: 56.5 % are 1-bp changes in homopolymers ≥ 5 bp at AF < 0.12).

So the HPRC graphs remove 84–90 % (v1.1) and 85–91 % (v2.1) of the GRCh38-only graph's germline candidates, at a cost of 0.3–3.6 points of candidate-level somatic SNV sensitivity. Somatic INDELs lose more, but 73.0–85.6 % of the lost INDELs are population alleles that the filter removes on any arm (R1.4). Lost alleles whose ALT the graph spells are 95.5–100 % tagged; those it does not spell less often (71.8–87.4 % on HG008T, 7.8–57.5 % on COLO829T). Among alleles the filter keeps, GRCh38-only has a candidate for 0.2–0.6 % more of all somatic SNVs and 2.0–8.2 % more of all somatic INDELs (chr1–22, tensor level). Absorption is a small part of this: in the chr1 test regions the graph arms lose 10–31 untagged INDEL candidates per dataset and graph, 0–2 of them absorbed [待补: whether GRCh38-only calls these INDELs].

This matches the graph's content: in the benchmark regions, 87 % of matched-normal germline SNVs and 75 % of germline INDELs are HPRC v1.1 graph alleles, against 0.2–2.5 % of somatic SNVs and 21–25 % of somatic INDELs. [待补: reconcile with the original 79 % / 73 %, which used another definition]

Inside the classifier, the number of HPRC haplotypes on each node (ch6) suppresses the remaining germline candidates: setting it to "every haplotype" at inference multiplies those called somatic 3.4–7.6-fold in all 12 HPRC-graph models, though filtering removes most again (R1.1, R2.5) [待补: retraining without the channel].

At the call level (R3.1), PASS SNV F1 before filtering is 0.008–0.075 on GRCh38-only against 0.18–0.69 on the HPRC graphs; after filtering, 1.4–7.8× more false calls remain on GRCh38-only than on v1.1, and the HPRC graphs keep +0.04 to +0.32 SNV F1 (GRCh38-placed calls, benchmark regions). Graph-coordinate calls are mostly population or patient-germline alleles (three v2.1 SNV runs, HG008T PacBio / HG008T ONT / COLO829T ONT, chr1: the normal assembly carries 57.6–83.1 % of them; R2.5). Counted all false, they put GRCh38-only ahead in 4 of 12 SNV comparisons, all long-read. All come from the off-reference rescue: without it, complete-output F1 (no BED) of the three v2.1 runs is at least 0.596 / 0.524 / 0.721 against 0.359 / 0.365 / 0.418, which reverses both v2.1 cases [待补: the rescue default; the two v1.1 cases] (R2.5, R3.4).

False calls left after filtering on v2.1 are mostly germline alleles the graph does not encode: in the same three runs the normal assemblies carry the ALT of 146, 135 and 110 of 150 sampled false calls, and only 2.0–10.0 % are spelled by HPRC v1.1 haplotypes or the v2.1 graph. [待补: E4b — the GRCh38-only arm's false calls split the same way, to test whether its excess is germline alleles the graph encodes]

**证据（内部）。**

- 三臂 truth 保留率的表来自 `response_to_reviewers_detailed.md` R1.2（2026-10-12，从 18 个 set 的 `tensors/truth_recall.json` 算出，status = `tensor_representative`）。
  - 本次重跑 `analysis/graph_ablation_tensors_20261005/collect.py`，输出写到 session scratchpad：`/tmp/claude-10008/-scratch-jshen-Github-Pansoma/832836ae-84c0-4629-b961-f3e655be6f4d/scratchpad/ablation_summary_20261013.json`。这是临时文件，没进仓库；复现命令是 `python collect.py <out.json>`。
  - 18 个 set 都已 finalized，重跑结果和 R1.2 一致。例：HG008T Illumina germline 3,832,002 / 389,473 / 363,126 of 5,061,892 = 0.757 / 0.077 / 0.072。
  - d46 deletions 是 72.95 %，doc 写的是 0.729，照抄。
- Background（label 0：BED 内、有 GRCh38 位置、没有 truth 的 tensor；branch 节点上没有 truth 的 tensor 是 −1（R9），不在这里），linear → d9 → d46：
  - SNV non：HG008T Illumina 2,114,925 → 1,876,969 → 1,852,495（−11.3 % / −12.4 %）；PacBio 188,910 → 73,454 → 77,478（−61.1 % / −59.0 %）；ONT 339,711 → 197,575 → 206,573（−41.8 % / −39.2 %）。
  - SNV non：COLO829T Illumina 1,285,779 → 934,450 → 920,461（−27.3 % / −28.4 %）；fiberseq 250,913 → 120,489 → 162,510（−52.0 % / −35.2 %）；ONT 577,417 → 410,794 → 546,170（−28.9 % / −5.4 %）。
  - INDEL non：Illumina −77.7 % / −83.5 %（HG008T）、−74.9 % / −81.1 %（COLO829T）；长读长 −0.8 % 到 +3.9 %（例：COLO829T fiberseq 1,272,815 → 1,310,694 → 1,321,934）。
  - Total tensors 下降：d9 33.1–68.8 %，d46 34.4–69.8 %。
  - 长读长 INDEL non 的组成：`tmp/indel_non_20261006/rules.txt`（chr1-22，COLO829T_fiberseq，non 1,310,694，就是 d9 set 的 INDEL non 数）。规则 "1-bp in HP>=5: AF >= 0.12" 的 non_removed 0.565，也就是 56.5 % 的 non 是 HP ≥ 5 里的 1-bp INDEL，而且 AF < 0.12。
- 和旧数字的出入：简版 / 详版写的 "43–69 %" 不含 COLO829T ONT linear（10-05 的 `tables.md` 里还没有这一套）。加上它（10,399,773 → 6,958,282，−33.1 %）后是 33–69 %，以新数字为准。
- somatic SNV 损失 0.3–3.6 个点：
  - d9：linear − d9 = 3.21 / 3.61 / 3.44（HG008T）、0.85 / 0.45 / 0.33（COLO829T）个点；
  - d46：2.82 / 3.40 / 3.18、0.86 / 0.50 / 0.36。
  - 详版写的是 "0.3–4"，这里写成实际范围。
- germline candidate 被去掉的比例（本次核对，同一个 scratchpad JSON，`germline_rep`，1 − graph / linear）：
  - d9：HG008T 0.8984 / 0.8643 / 0.8506，COLO829T 0.8966 / 0.8576 / 0.8392（Illumina / PacBio-fiberseq / ONT），范围 84–90 %。
  - d46：HG008T 0.9052 / 0.8789 / 0.8619，COLO829T 0.9041 / 0.8686 / 0.8487，范围 85–91 %。
  - 12:00 版写的 "about nine in ten" 偏高（最低 0.839），改成实际范围。
- 图的 germline 内容：`experiments/E7_graph_composition/d9_content.tsv`（scope = bed，`in_d9_vcf`）。
  - HG008 germline SNV 0.8697、INDEL 0.7513；COLO829 0.8721、0.7457。
  - somatic SNV 0.0245（HG008）、0.0022（COLO829）；somatic INDEL 0.252、0.2108。
  - E7 的 README 已经写好（2026-10-13，HPRC v1.1 部分 DONE，v2.1 轴要 v2.1 全图 VCF）。README 说明 s5 没有拆 d9 的 COMPLEX allele，所以 `in_d9_vcf` 是下界。
- Classifier 层面（新，E1.1）：`experiments/E1.1_ch6_counterfactual/results/tables.md`，以及 `results/summary.json`、`results/profile_auc.json`。用的是 18 个 7fc7396 SNV 模型，chr1 test tensors。
  - Germline 的 site code 中位数（§1）：d9 的 HG008T 是 96 / 96 / 96，COLO829T joint 是 96 / 95 / 95。d46 的 HG008T 是 80 / 79 / 79，COLO829T 是 80 / 76 / 74。somatic 全部是 100。
  - ch6 单独的 AUC（§6，somatic vs germline，按 site code）：`HG008_PacBio` 0.837、`HG008_PacBio_d46` 0.890、`HG008_ONT_d46` 0.892。只算了这三个 run。
  - all_H 下 germline ≥ t 相对 base 的倍数（§3 的 "x base"）：
    - HG008T d9：3.58 / 3.50 / 3.79；HG008T d46：6.75 / 6.83 / 6.64。
    - COLO829T d46：7.57 / 6.72 / 7.39；COLO829T d9 joint：3.45 / 3.41 / 4.19。
    - 范围 3.41–7.57。
  - all_H 下 PoN 后 F1 的变化（§4，base → all_H）：
    - HG008T d9：0.589 → 0.573、0.621 → 0.566、0.596 → 0.533。
    - HG008T d46：0.560 → 0.511、0.605 → 0.561、0.579 → 0.516。
    - COLO829T d46：0.784 → 0.791、0.785 → 0.760、0.743 → 0.705。
    - COLO829T joint：0.675 → 0.682、0.782 → 0.778、0.706 → 0.664。
    - 范围 −0.063 到 +0.007。
    - E1.1 README 说明了这个 F1 是 per-tensor 的估计，只用来比较模式，不是 rtg F1。
  - Linear 模型的 ch6 是常数（H = 17，每个格子都是 100），所有模式结果都一样。所以 linear 臂没有这一层。
  - 这是推理时的 "有没有用" 测试，不是 "有没有好处"。好处要靠去掉 ch6 重训来看（`E1.1_ch6_counterfactual/gpu_nopc_commands.sh`）。
  - Response 的 classifier 段缩短了：只保留 3.4–7.6 倍；AUC、per-tensor F1 的 −0.063 到 +0.007 和 site code 中位数只在 R2.5 里写。
- call 层面：
  - `experiments/E3_reference_arms/README.md`："Raw PASS F1 is 0.008–0.075 against 0.18–0.69"，"+0.04 to +0.32"。
  - `benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/README.md`："false calls after it 1.4-1.7x d9's on HG008T, 1.8-7.8x on COLO829T (Illumina 2,216 vs 284)"（在 BED 内）。
  - 这些都只算 placed calls。unplaced 全算 FP 时（`E3_reference_arms/tables.md` §2，no BED）：HG008T ONT d9 / d46 0.342 / 0.354 < linear 0.365；COLO829T ONT d46 0.369 < 0.418；COLO829T fiberseq d9 0.636 < 0.680。所以 Response 里给 +0.04 到 +0.32 加了这个限定（细节见 R3.4）。
  - 三个 d46 run 上 E9c 之后的数（见下面 E9c / E9d 一条）：graph PoN 后 complete-output F1（no BED）0.5901 / 0.4889 / 0.4911，都高于 linear。但 graph PoN 也标掉大部分 gold off-ref somatic key，R3.4 只把这组数当参考，R2.5 的 Response 也写成 "a reference only"。所以本点 Response 把原来的 "the advantage disappears in some long-read runs" 换成 "4 of 12 SNV comparisons" 加不用 rescue 的下界 0.596 / 0.524 / 0.721（R3.4 证据，F1 = 2TP / (2TP + FP + FN + U)，U 不含 AF floor，所以是下界；我按 R3.4 的 TP / FP / FN / U 重算一致）。另外两个反转（HG008T ONT d9、COLO829T fiberseq d9）还没做，Response 里是 [待补]。
  - "4 of 12" 只算 SNV（E3 §2）：INDEL 的 complete-output 在长读长上 linear 本来就常常更高（例 HG008T PacBio INDEL linear 0.121 vs d9 0.054、d46 0.035），所以 Response 写明 "SNV comparisons"。
- 丢掉的 somatic allele 是什么（新，E5b final）：`experiments/E5b_absorbed_d46/tables.md`（`s5_tables.py` 生成），还有 `summary.json`、`t2_compare_d9.txt`。范围是 chr1-22 全部 truth record，tensor 层面（是不是 PASS call 不在这里）。"PoN-tagged" 是 truth allele 被四个 PoN 资源标掉。
  - §2，lost vs linear（linear 上有 tensor、图上没有）里的 absorbed（没有 candidate，图能拼出 ALT）的 tagged 比例，d9 / d46：
    - HG008T SNV 三个平台都是 100.0 %；HG008T INDEL 97.8 / 97.0（Illumina）、96.4 / 96.2（PacBio）、95.5 / 95.8（ONT）。
    - COLO829T SNV 98.2 / 100.0（Illumina），其余 100.0；COLO829T INDEL 98.9 / 98.1、98.4 / 97.5、99.1 / 98.4。
    - 范围 95.5–100.0 %，Response 写 "95.5–100 %"。
  - §2 "no candidate, ALT not spelled" 的 tagged 比例：HG008T SNV 73.2–87.4 %，HG008T INDEL 71.8–78.7 %，COLO829T SNV 7.8–17.5 %，COLO829T INDEL 52.8–57.5 %。Response 合并成 HG008T 71.8–87.4 %、COLO829T 7.8–57.5 %。
  - 注意：lost 的全部（不分子类）tagged 比例在 COLO829T SNV 只有 23.4–33.8 %（例 Illumina d9 395 个里 98 个），所以 Response 不写 "most lost alleles are tagged"，只写 absorbed 那一类。
  - §2 INDEL 的 lost（全部）tagged 比例：HG008T 85.6 / 85.1、82.0 / 81.4、83.8 / 83.4 %，COLO829T 80.0 / 79.6、75.2 / 73.4、74.3 / 73.0 %（d9 / d46），范围 73.0–85.6 %。Response 写 "73.0–85.6 % of the lost INDELs are population alleles that the filter removes on any arm"。
  - §2 PoN 后 lost INDEL 的子类（E5b README Interpretation 3）：absorbed HG008T 35 / 34 / 20（d9）、39 / 32 / 20（d46），COLO829T 4 / 3 / 1、6 / 5 / 2；"no candidate, ALT not spelled" HG008T 145–361，filtered 60–153，other allele 19–83。所以 PoN 后的 INDEL 损失主要不是 absorption，Response 写 "Absorption is a small part of this"，并用 §8 的 chr1 BED 数（10–31，absorbed 0–2）。
  - 旧稿的 "because repeat-length somatic changes coincide with population alleles in the graph" 删掉了：E5b README Interpretation 3 说 absorption 只解释一小部分损失，§6 也看不出丢掉的 INDEL 偏向 HP ≥ 7（见下面 §6 一条）。
  - §1 "tensor and not tagged" 的 d9 − linear / d46 − linear（占全部 truth）：
    - SNV −0.3 / −0.3、−0.3 / −0.3、−0.3 / −0.2（HG008T Illumina / PacBio / ONT），−0.6 / −0.6、−0.3 / −0.3、−0.2 / −0.2（COLO829T），范围 0.2–0.6 %。
    - INDEL d9 −4.8 / −6.0 / −2.0（HG008T）、−6.8 / −8.2 / −4.6（COLO829T）；d46 −4.4 / −5.8 / −2.0、−6.2 / −8.1 / −4.9。范围 2.0–8.2 %。
    - 这是净差（图上多出来的 "gained" 也算进去了）。
  - §5 / §7：图能拼出 ALT 的 truth INDEL 里 96.2 %（HG008T）和 97.7 %（COLO829T）被 PoN 标掉（d46）；absorbed 在 PoN 后只占 truth INDEL 的 0.2–0.7 %（HG008T / COLO829T，d46，§7 第二张表；HCC1395T 是 0.0 %）。
  - §6：absorbed INDEL 大多在 HP ≥ 7 里（HG008T d9 1,473 / 1,818）。但 HG008T 的 truth INDEL 本身就有 6,229 / 8,496 在 HP ≥ 7 里；lost on ≥ 1 platform 是 3,155 / 4,275（d9）。COLO829T d9 lost 428 / 799，truth 1,110 / 1,912。所以丢掉的 INDEL 相对 truth 没有明显偏向 HP ≥ 7，不能拿来支持 "repeat-length changes" 的说法。
  - §8（chr1、BED 内，就是 vcfeval 的范围）：PoN 后丢的 INDEL HG008T 18–31、COLO829T 10–13，其中 absorbed 0–2。
  - `t2_compare_d9.txt`：E5b 用 vg find 子图判断 "ALT spelled"。在 300 个 HG008T read-level absorbed locus 上和完整 d9 edge table（sqlite）对比，278 个一致（234 个两边都 spelled，44 个两边都不 spelled）；另外 22 个 sqlite 说 spelled、vg find 子图说没有（sqlite spelled 的 256 个里 8.6 %，E5b README Caveat 2）。所以 d9 的 absorbed 是下界，"not spelled" 里混了一部分 absorbed（只在这 300 个 locus 上测过）；d46 没有同样的检查。
  - E5b 的 README 已经写好（23:00），`tables.md` / `summary.json` 是 final（22:53）。README Interpretation 1 的 "95.5-97.8 % of HG008T INDELs, 97.5-99.1 % of COLO829T INDELs and 98.2-100 % of SNVs" 和上面的 95.5–100 % 一致。
  - `linear_recovery.md`（丢掉的 INDEL 在 linear / d9 / d46 上 PASS 能不能叫出来）还是全部 pending：文件末尾列了 36 个缺的 rtg 输出（18 个 INDEL run × `pon_curve` / `pon_calls`）。E5b README 写的另一条路是 18 个 INDEL in-sample re-predict（tequila GPU，`tasks_full_plan.tsv` 的 `insample_linear` / `insample_d9` / `joint_d9` / `insample_d46` 行），在等用户定 cluster cap。
  - E5b README Caveat 6：`s5_tables.py` 的 `pon_tags()` 留了 INFO 末尾的换行，COLO829T 的 PoN-tagged 记录在 `per_truth.tsv.gz` 里被拆成两行。`tables.md` / `summary.json` 不受影响，但 `linear_recovery.md` 里 COLO829T 的 "absorbed" n 偏少（例 Illumina d9 文件 1、重新拼回 26）。跑 s6 之前要先修。
- Unplaced calls 和 d46 臂 FP 的拆分（新，E9c / E9d）：`experiments/E9c_offref_filter_validation/tables.md` 的 `rcpu__*` 行（collect job 387434，16:07）；方法、结果和 caveat 见同目录 `README.md`（23:04 已写好）。三个 d46 SNV run 在 tequila GPU 上 in-sample re-predict，chr1，顺序是 HG008T PacBio / HG008T ONT / COLO829T ONT。只有这三个 run，都是长读长 SNV；细节和 caveat 也写在 R2.5 的证据里。
  - §1 unplaced PASS（互斥分类，优先顺序 population > germline > tumor-only > unverifiable > unsupported）：population（graph PoN）89.2 / 74.1 / 49.8 %；germline（normal assembly）7.9 / 19.2 / 45.4 %；tumor-only 0.4 / 0.2 / 0.0 %（COLO829T 没有 tumor assembly，not in normal 2.4 %）。COLO829T ONT 处理的是随机抽的 3,000 / 4,569。
  - 只看组装、不按优先顺序（README §1 第二张表）：germline 65.8 / 57.6 / 83.1 %，Response 用的 "57.6–83.1 %" 就是这组；tumor-only 9 / 22 个（3.8 / 1.7 %），其中 graph PoN 标掉 8 / 19。
  - §2 graph PoN 后 F1（no BED）0.5901 / 0.4889 / 0.4911；linear 来自 E3 §2（RIS）0.359 / 0.365 / 0.418。本点 Response 不再写这组数（只当参考，见上面 call 层面一条）。
  - §5 gold：chr1 上有 tensor 的 gold key 只有 11 / 9 / 2 个，graph PoN 标掉 9 / 7 / 2；全染色体主规则标掉 74.2–88.3 %（README §4）。graph PoN 只能拿来分类，不能当保留 somatic 的过滤器；怎么用由用户决定（E9c README "Needs"）。
  - §3 controls（各抽 150）：
    - placed FP（PoN 后）germline 146 / 135 / 110，graph PoN 标 4（2.7 %）/ 15（10.0 %）/ 3（2.0 %）。
    - placed TP 是 tumor-only 149 / 145，COLO829T ONT not_in_normal 150。说明组装验证在 placed 位点上分得开 somatic（在 off-ref 的 gold 上没这么准：每个 gold 集有 9–60 个 somatic key 被判成 germline，README Caveats）。
    - ctrl_fp 的 germline 里带 `multi_copy`（context 在一个组装里 > 2 个位置）的有 2 / 146、16 / 135、19 / 110（我从 `per_call.tsv.gz` 数的），这些可能是 paralog。
    - GRCh38 PoN 标掉的 placed PASS（ctrl_pon）germline 147 / 140 / 144，graph PoN 只标 17（11.3 %）/ 13（8.7 %）/ 11（7.3 %）。也就是说，还能成为 candidate 的 population allele 多数不在 HPRC haplotype 里，和 "图里编码的 germline 被 alignment 吸收" 一致（推断）。这一点没放进 Response。
  - graph PoN 的 haplotype 用的是 HPRC v1.1 的 89 条（v2.1 全图 VCF 不在本地），加上 d46 自己的 graph walk，所以 Response 写 "spelled by HPRC v1.1 haplotypes or the v2.1 graph"。用 v2.1 的 haplotype，2.0–10.0 % 可能会高一些。
  - 这只是 d46 臂的拆分。R3.2 要的是 "linear 臂多出来的 FP 是不是图里编码的 germline"，这要 linear 臂的 FP（E4b）。
  - E9a（投影，`experiments/E9a_unplaced_projection/results/tables.md`）：三个 run 只有 2.8–7.9 %（19 / 43 / 128 个）的 unplaced PASS 能写成 GRCh38 allele（class i，其中 E9a 的 "exact" 子类只有 2 / 19 / 94 个），class ii 的区间里 0 个 truth TP（BED 内 0 / 208、0 / 940、0 / 1,334），所以投影后的 complete-output F1 基本等于 worst case（no BED，PoN 后 0.533 / 0.355 / 0.369 → 0.537 / 0.359 / 0.373；R2.5 / R3.4 有细节）。CHM13 列是 0，因为用户取消了 CHM13。
- 可比性：d9 和 linear 的配对只在 HG008T Illumina SNV 和 HG008T INDEL 上完全一致。
  - HG008T PacBio / ONT SNV 的 d9 用 1 GPU（global batch 1,024），linear 用 4 GPU。
  - COLO829T d9 INDEL 是 lr 4e-4。
  - 见 E3 README Caveats 1。

**还缺什么。**

- [CPU 本地] E4b：
  - d46 臂已经有抽样的拆分（E9c §3 controls，三个 d46 run）：PoN 后的 placed FP 大多是图里没有的 germline。还缺 linear 臂（和 d9 臂）的同样拆分，这是这一点的核心。
  - 需要 linear 臂 12 个 run 的 chr1 call sets。做法是本地 re-predict：用 `tmp/reviewer_cpu_20261013/tasks_full_plan.tsv` 里的 12 行 `insample_linear`，它们还没进当前 `tasks.tsv`。predict 按规定在 tequila GPU 上跑，和 387379 那 18 个任务一样走 `task_gpu.sh`，然后接 `vcf.sh`。
  - 然后把 linear / d9 / d46 PoN 前后的 FP 和 dipcall germline 匹配，再按 allele 是否在 d9 deconstruct VCF（`AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz`）里分三类。
  - d46 没有现成的 VCF：用 `vg deconstruct` 生成，或者用 tensor key 判断。
  - 也可以试 E9c 的 ctrl_fp 流程（`e9c_run.py`：graph PoN + 组装比对），这样三臂用同一个定义；它在 linear 图的 run 上能不能跑还没确认。
- [数据同步] 另一个办法：从 RIS 同步 24 个 linear / d9 / d46 run 的 `vcf_chr1/`（rtg 的 fp.vcf）。本地的 `ris_7fc7396_20261012/` 只有 `best.pth` 和 `test_chr1/*.metrics.json`。
- [数据同步] E5b `linear_recovery.md`：丢掉的 INDEL 在 linear / d9 / d46 上 PASS 能不能叫出来，要从 RIS 同步 36 个 rtg 输出（每个 INDEL run 的 `pon_curve` 和 `pon_calls`，清单在文件末尾），然后重跑 `s6_linear_recovery.py`。另一条路（E5b README Status）：18 个 INDEL in-sample re-predict 在 tequila GPU 上跑，再接 `vcf.sh`，这要等用户定 cluster cap。补 Response 里的 [待补: whether GRCh38-only calls these INDELs]。
- [CPU 本地] 跑 s6 之前先修 E5b `s5_tables.py` 的 `pon_tags()` 换行 bug（`x[7].rstrip('\n')`，README Caveat 6），再 `./collect.sh`。不修的话 `linear_recovery.md` 里 COLO829T 的 absorbed n 偏少（例 Illumina d9 1 vs 26）。
- [CPU 本地] E9a / E9c / E9d 的另外两个反转（d9：HG008T ONT、COLO829T fiberseq `_h90`）：先在 tequila GPU 上 re-predict，再跑（和 R2.5 / R3.4 共用）。
- [用户决定] off-ref rescue 是否保留为默认（R3.4）：决定 call 层面 complete-output 的主数是 rescue 关（至少 0.596 / 0.524 / 0.721），还是 rescue 开（worst case 0.533 / 0.355 / 0.369；扣掉 graph PoN 标的 0.590 / 0.489 / 0.491 现在只当参考）。
- [用户决定] graph PoN 只做分类还是也当过滤器（E9c README "Needs"，和 R2.5 / R3.4 共用）；还有要不要下载 HPRC v2.1 全图 VCF（E7 也要），让 d46 的 "spelled by HPRC haplotypes" 用 v2.1 的 haplotype。
- [CPU 本地] E4a：把本次重跑的 `collect.py` 结果写回 `analysis/graph_ablation_tensors_20261005/`（`tables.md`、`summary.json`），并加 HCC1395 的 d46 列。目前只在 scratchpad 里。
- [GPU RIS] 去掉 ch6 重训：`experiments/E1.1_ch6_counterfactual/gpu_nopc_commands.sh`（`train --drop-planes path_count`），共 27 个 run，其中 24 个 no-ch6、3 个 COLO829T d9 INDEL 的 with-ch6 对照。它回答 "classifier 层面的 germline 压制在 PoN 之后还剩多少好处"。
- [GPU RIS] d9 的配对重跑：HG008T PacBio / ONT 和 COLO829T 三套 SNV 用 4 GPU；COLO829T INDEL 用 lr 1e-4。这样 E4b 的 linear vs d9 是同一 recipe（E3 Caveats 1）。
- [文字] 简版 R3.2 要改：
  - "43–69 %" 改成 33–69 %（含 COLO829T ONT），"0.5–4 points" 改成 0.3–3.6；
  - 表换成三臂版；
  - 说明原稿 79 % / 73 % 的定义，和 E7 的 87 % / 75 % 对齐。
- [文字] E7、E5b、E9c 的 README 都已写好；本次已按各自的 `tables.md` / README 核过 Response 里的 95.5–100 %、73.0–85.6 %、0.2–0.6 % / 2.0–8.2 %、10–31 / 0–2、57.6–83.1 % 和 146 / 135 / 110，都一致。剩下的是简版 R3.2 要跟着改：删掉 "because repeat-length somatic changes coincide with population alleles in the graph"，"4 of 12" 写明是 SNV，complete-output 的主数等 rescue 的用户决定。

### R3.3 Validation of graph-coordinate variants

**状态：部分解决**: 关键证据已经有了。三个 d46 最终 run（HG008T PacBio、HG008T ONT、COLO829T ONT；chr1，本地 in-sample re-predict：PacBio 是 CPU fp32，两个 ONT 是 tequila GPU）的 unplaced SNV PASS calls 做完了三件事：组装验证（E9d）、graph-native PoN（E9c）、投影到 GRCh38（E9a）。结论是 HG008T 上按 population 优先的互斥分类，tumor-only 只有 0.4 / 0.2 %（1 个和 3 个 call）；只看组装是 9 / 22 个（3.8 / 1.7 %），其中 8 / 19 个也被 graph PoN 标成 population。其余主要是 population allele（89.2 / 74.1 / 49.8 %）或 germline（7.9 / 19.2 / 45.4 %；COLO829T 只随机处理了 3,000 / 4,569）。对照说明 assay 本身是有效的：HG008T 上 placed TP 145–149 / 150 是 tumor-only，PoN 后的 placed FP 135–146 / 150 在 normal 组装里（COLO829T 是 150/150 不在 normal、110/150 在 normal）。Response 第 2–5 点已换成最终模型的数。证据的边界：只有三个 run、只有 chr1，chr1 上的 gold 每个 run 只有 2–11 个。还缺三样：一是 HPRC v1.1（d9）的最终 run 和 Illumina / fiberseq 没有做 E9a / E9c / E9d（predictions 在 RIS 上，要数据同步）；二是 COLO829T 只能做 germline 排除（没有 tumor 组装）；三是 rescue 是否保留要用户定，这决定 Response 结尾怎么写。

**Comment（要点）。**

- A subset of the final predictions exists only in graph coordinates. The evidence for them (signatures, read examples) is indirect.
- Validate a representative sample orthogonally: matched-normal alignment evidence, phased haplotype assemblies, or cross-platform data.

**Response.**

We agree. We no longer present graph-coordinate calls as a strength, and every headline metric also counts them as false positives (R3.4). So far we validated the final models' graph-coordinate SNV PASS calls on chr1 in three HPRC v2.1 runs: HG008T PacBio (240 calls), HG008T ONT (1,258) and COLO829T ONT (4,569; a random 3,000 analysed).

1. *Origin.* All come from the off-reference rescue. The classifier treats a low haplotype path count as evidence against a somatic call (R1.1), and reads on nodes off GRCh38 have low counts: the final HPRC v1.1 models call 0 of 112 (PacBio) and 1 of 89 (ONT) assembly-verified off-reference somatic events of HG008T at the PASS threshold. The rescue raises the count around the site and calls 32 and 46 (tensor level, all autosomes).
2. *Placement.* Through the GRCh38 nodes flanking their bubble, only 7.9, 3.4 and 2.8 % become a GRCh38 allele, none a true positive after population filtering; none of the 221, 984 and 1,334 placed only to an interval matches a missed truth variant (R3.4).
3. *Assemblies.* We aligned each call's ALT context to the HG008-T tumor and both HG008-N haplotype assemblies (COLO829T: the COLO829BL normal only). Controls (150 each): 149 and 145 placed true positives are tumor-only (COLO829T: 150 absent from the normal); 146 and 135 placed false positives after population filtering are in the normal (COLO829T: 110). In exclusive classes, in this order, the calls are population alleles (ALT context spelled by GRCh38, an HPRC haplotype or a graph walk) in 89.2, 74.1 and 49.8 %, germline in 7.9, 19.2 and 45.4 %, tumor-only in 0.4 and 0.2 % (COLO829T: 2.4 % absent from the normal), and unverifiable or in no assembly in the rest. By the assemblies alone, 9 (3.8 %) and 22 (1.7 %) are tumor-only, 8 and 19 of them also population alleles.
4. *Population filter.* The graph-native filter tags 0 or 1 of 150 placed true positives per run, but also 9 of 11 and 7 of 9 assembly-verified off-reference somatic events on HG008T chr1, including every PASS one (all chromosomes: 90 of 119 and 72 of 97). It describes calls and cannot certify them.
5. *Cross-platform.* 190 calls (same node, offset and allele) are PASS on both HG008T PacBio and ONT (79.2 % of PacBio's, 15.1 % of ONT's), 167–169 of them population alleles.

On HG008T at most 3.8 % and 1.7 % of graph-coordinate calls are assembly-supported somatic candidates; the rest are mostly population or germline alleles beyond a GRCh38 population filter, so in these runs counting them as false positives (R3.4) is close to their true status. The rescue is a stated limitation; training with assembly-derived off-reference labels is the proper fix. [待补: the rescue default (user decision, R3.4)]

**证据（内部）。**

- 最终模型的数据：`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/`，是 3 个 d46 SNV run 的本地 in-sample re-predict，都带 `vcf_chr1/eval_bed`、`eval_nobed`。HG008T PacBio 是 CPU fp32（job 387135，host tsingtao，log `tmp/reviewer_cpu_20261013/logs/rcpu_test_HG008_PacBio_SNV_d46-387135.out`），两个 ONT 是 tequila GPU bf16（387379_1 / _2）；见 E9a README 的 Data 表和 Caveat 2。之前写的 "都在 tequila GPU 上" 不对。
  - eval_bed PoN F1 是 0.603 / 0.575 / 0.736。RIS 是 0.605 / 0.575 / 0.737（`benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/tables.md`）。
  - unplaced PASS 是 240 / 1,258 / 4,569。RIS 是 244 / 1,261 / 4,566。
  - Response 里的顺序都是 HG008T PacBio、HG008T ONT、COLO829T ONT。
- Origin（E1.1 §5）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §5 和 `results/assembly_keys.json`。
  - 用的是 7fc7396 的 d9 模型 `HG008_PacBio_SNV_7fc7396` 和 `HG008_ONT_SNV_7fc7396`。keys 来自 `assembly_offref_20260929/offref_somatic.tsv`，范围是全部常染色体。
  - key 数：PacBio 112 个（chr1 8 个），ONT 89 个（chr1 8 个），全部 `off_reference`。
  - 标签：112 里有 79 个是 label 1 / `residual_partial_somatic_truth`（被 ignore），33 个是 −1；89 里是 61 / 28。README 原话："None of them is a training label"。
  - 格子的格式是 "p ≥ t / argmax somatic / n"：

    | 模式 | PacBio | ONT |
    |---|---|---|
    | base | 0 / 0 / 112 | 1 / 1 / 89 |
    | all_H | 61 / 44 / 112 | 49 / 49 / 89 |
    | site10_high | 45 / 32 / 112 | 45 / 46 / 89 |
    | chr1 上的 base | 0 / 8 | 0 / 8 |

  - `site10_high` 就是 predict 的 rescue 变换（E1.1 README）。pipeline 对 off-ref tensor 用 rescue 之后的 argmax（`predict.py` docstring、`graph_vcf` 的 `offref_call`），所以 rescue 在 tensor 层面叫出 32/112、46/89。graph_vcf 的 PASS 还要 AF ≥ floor（LowAF），这一步没算，所以 32 / 46 是 PASS 的上界。
  - 每个 key 都有事件的 GRCh38 表示（`assembly_keys.json` 的 `grch38` 字段）。但 off-ref tensor 的 call 是 unplaced，rtg 不给它分，complete-output bound 还把它算 FP。
  - 精度：按 E1.1 README，`HG008_PacBio`（d9）的 base / all_H / train_mean 和 `HG008_ONT`（d9）的 base 是 CPU fp32，其余模式（包括两者的 `site10_high`、ONT 的 all_H）是 tequila GPU bf16。§2 里这两个 run 的 base AP CPU / GPU 是 0.2333 / 0.2331、0.2316 / 0.2323，18 个 run 都在 0.0047 以内。
  - d46 上 rescue 在 PASS 层面的效果（新，E9c §5）：chr1 上有 tensor 的 gold key，PacBio 11 个里 3 个 PASS（unplaced），ONT 9 个里 5 个，COLO829T 2 个里 1 个。这比 E1.1 §5 的 d46 版窄（只有 chr1、只有 rescue 之后），但已经是 PASS 层面（含 AF floor）的数。
- Rescue 就是 unplaced（本次之前已核对）：`ris_7fc7396_20261012/<run>/test_chr1/*.metrics.json` 的 `offref_rescue.calls`，和 `E3_reference_arms/tables.md` §2 的 unplaced PASS 完全相等。
  - HG008T d9：116 / 211 / 1,409。HG008T d46：754 / 244 / 1,261。COLO829T d46：516 / 258 / 4,566。
  - rescue 之前 argmax 判为 somatic 的（`calls_before`）：3 / 7 / 30、9 / 3 / 48、43 / 106 / 262。
  - E3 里 COLO829T 的 d9 SNV 是 `*_SNV_h90`（2026-10-03），metrics.json 不在本地，没核对。按设计 predict 默认 `--offref-site90 10`，每个 SNV off-ref tensor 都会被 rescue。
- Placement（E9a 最终，新）：`experiments/E9a_unplaced_projection/results/tables.md`、`projection_by_run.tsv`、`metrics_by_run.tsv`、`summary.json`。
  - unplaced PASS 240 / 1,258 / 4,569。class i（写成 GRCh38 variant）19 (7.9 %) / 43 (3.4 %) / 128 (2.8 %)；class ii（只有 GRCh38 位置或区间）221 / 984 / 1,334；class iii（没有一致的 anchor 对）0 / 231 / 3,107。
  - class i 里 exact 2 / 19 / 94，branch_base_differs（call 落在一个 germline 差异上）17 / 24 / 34。GRCh38 PoN 标掉 class i 里的 18 / 38 / 88 个。
  - class ii 细分（`summary.json` 的 `details_pass`）：in_insertion 108 / 383 / 589，alt_is_grch38_base 89 / 517 / 608，near_branch_indel 24 / 84 / 116，too_long 0 / 0 / 21。class iii：no_left_anchor 0 / 222 / 2,836，no_right_anchor 0 / 9 / 271。
  - rtg：no BED、PoN 后，加上 class i 后 TP 不变（586 / 603 / 1,690）；PoN 前 HG008T 两个 run 各多 1 个 TP（611 → 612、631 → 632），COLO829T 不变。
  - class ii 的区间里 truth TP 为 0：no BED 0/221、0/984、0/1,334；BED 内 0/208、0/940、0/1,334。
    - 规则（`evaluate_complete.py` docstring 和 `match_intervals`）：只有 rtg 留下的 FN truth 能被认领（区间 ±1 bp，每个 truth 只认领一次）；`alt_is_grch38_base`（89 / 517 / 608）直接算 FP，不查区间。真正查过区间的是 132 / 467 / 726 个（no BED）。
    - 所以 Response 写的是 "matches a missed truth variant"，不能写成 "区间里没有任何 truth"（已改）。
  - `tables.md` 第二张表的 "ii INDEL context" 列原来是 0（`collect.py` 读错 key）。23:08 已修好并重跑（E9a README Caveat 3），现在显示 24 / 84 / 116；`metrics_by_run.tsv` 的 F1 和 TP / FP / FN 重跑前后一致（23:10 核对）。
  - CHM13 列全是 0（`chm13_pass` = `not_built`），因为用户取消了 CHM13。
  - 旧的 dev 结果（legacy d9 `HG008_ONT_SNV_pow075_S_lr1e4_e20`，`tmp/reviewer_cpu_20261013/E9a_unplaced_projection/logs/dev.387367_0.log`）：855 个里 class i 24、alt_is_grch38_base 562。最终模型上 alt_is_grch38_base 的比例低一些，in_insertion 和 class iii 更多。
- 组装分类和 graph PoN（E9c / E9d 最终，新）：`experiments/E9c_offref_filter_validation/tables.md` 里的 `rcpu__*` 行；方法见 `README_methods.part`。
  - §1（互斥，顺序 population > germline > tumor-only > unverifiable > unsupported）：
    - HG008T PacBio 240：population 214 (89.2 %)、germline 19 (7.9 %)、tumor-only 1 (0.4 %)、unverifiable 5 (2.1 %)、unsupported 1 (0.4 %)。
    - HG008T ONT 1,258：932 (74.1 %)、242 (19.2 %)、3 (0.2 %)、70 (5.6 %)、11 (0.9 %)。
    - COLO829T ONT，随机处理 3,000 / 4,569：population 1,495 (49.8 %)、germline 1,362 (45.4 %)、not in normal 72 (2.4 %)、unverifiable 71 (2.4 %)、unsupported 0。
    - 定义：unverifiable = ALT 和 REF context 在哪个组装里都找不到；unsupported = ALT 不在任何组装里、REF 找得到。population = graph PoN 标的（GRCh38、HPRC v1.1 haplotype 或本图的 walk 拼得出 ALT context，`README_methods.part` 第 3 步）。
    - E9d 的 germline 类不是零错误：gold（按构造是 somatic）在 d46 上被判 germline 11/119（PacBio）、9/97（ONT）、33/175（Illumina）；ctrl_fp 里也有 2/150 是 tumor_only。COLO829T ONT 的 context 有 1,594/3,000 是 multi-copy，germline 可能含 paralog（E9c README Caveats）。所以 germline 7.9 / 19.2 / 45.4 % 是 "组装里找得到"，不是确证。
  - §3 只看组装、不先扣 population：unplaced PASS 里 tumor_only 9/240 (3.8 %)、22/1,258 (1.7 %)；COLO829T not_in_normal 197/3,000。即使不先扣 population，tumor-only 也很少。
    - 其中也被 graph PoN 标的：8/9、19/22、125/197（E9c README §1 第二张表；`per_call.tsv.gz` 里 `pop` 列重新数过，一致）。§1 的 1 / 3 / 72 就是减掉这些以后剩下的。
    - population 优先会把真的 somatic 也算成 population：被标的 tumor-only 里有 node 7923590 G>A（HG008T 两个 run 都有，E9a class i），就是 E9a README 里那个 truth SNV chr1:183,643,262 T>A（GRCh38 PoN 也标它）。所以 Response 结尾用只看组装的 3.8 / 1.7 % 作上界，不再用 0.4 %。
    - 跨平台共有的那 1 个 tumor-only 是 8299380:10:SNP:C>T（两个 run 都没被 graph PoN 标，class ii）。
  - §3 对照（每组随机 150 个）：
    - ctrl_tp：tumor_only 149 / 145（HG008T PacBio / ONT），COLO829T not_in_normal 150；graph PoN 标掉 0 / 1 / 0。
    - ctrl_fp（PoN 后的 FP）：germline 146 / 135，COLO829T 110；graph PoN 标掉 4 / 15 / 3。
    - ctrl_pon（GRCh38 PoN 标的 PASS calls）：germline 147 / 140 / 144。
  - §4 gold 校准（全染色体，d46 gold key，主规则 R(A,1)）：graph PoN 标掉 HG008T_PacBio_d46 90/119 (75.6 %)、HG008T_ONT_d46 72/97 (74.2 %)、HG008T_Illumina_d46 139/175 (79.4 %)。组装分类 tumor_only 103/119、81/97、137/175。
    - 注意：gold 本来就是从同一批组装推出来的，所以 103/119 这类数是一致性检查，不是独立的灵敏度。
    - 09-29 的旧 assay（`offref_negatives_20260929`，`v1_queries.py`）在 gold 上只确认了 26/177 tumor。E9d 换了查询方式（tensor 里的 ALT read context，minimap2 sr），所以 Response 里不再写 "lower bound"。
  - §5 chr1 gold（每个 run 里有 tensor 的）：PacBio 11 个，PASS 3、被标 9、PASS 且被标 3，组装 tumor_only 9 / germline 2；ONT 9 个，5 / 7 / 5，tumor_only 8 / germline 1；COLO829T 2 个，1 / 2 / 1。所以 rescue 叫出的 gold 全部会被 graph PoN 扣掉。
    - 数很小（3 + 5 个 PASS gold），只能当方向；全染色体的标记率（§4，90/119、72/97）更稳，Response 第 4 点两个都给了。
    - COLO829T 那 2 个是 HG008-T 的 gold key 碰巧也是 COLO829T 的 candidate，不是 COLO829T 的 somatic；d46 上其中 1 个在 COLO829BL 里是 germline（E9c README Caveats "Gold-set purity"）。
  - §8 跨平台：190 个 key 两个平台都是 PASS，占 HG008T ONT 的 1,258 里 15.1 %、HG008T PacBio 的 240 里 79.2 %。ONT 这边按类别是 population 167、germline 16、unverifiable 6、tumor_only 1；PacBio 这边是 population 169、germline 16、unverifiable 4、tumor_only 1。两个 run 的 tensor context 不同，所以分类差 2 个，Response 写 167–169。只有两个平台（没有 Illumina d46），所以只是 PacBio 对 ONT 的复现率。
  - §7 重复区（R1.3 用）：homopolymer ≥ 7 / STR / other = 97 / 60 / 83（PacBio）、387 / 219 / 652（ONT）、303 / 289 / 2,408（COLO829T，3,000 个）。
  - §2 的更正 F1 在 R3.4。
  - 限制：E9c 的 HPRC haplotype 检查用的是 HPRC v1.1 的 89 个 haplotype（v2.1 的全图 VCF 不在本地，README_methods 写明），d46 run 也一样；in-graph walk 用的是 run 自己的图。
- 旧模型的 E9c 行（`legacy__*`，09-30 / 10-01 的 pow075 模型）：unplaced PASS 被 graph PoN 标 Illumina 469/595（78.8 %）、PacBio 88/105（83.8 %），ctrl_tp 都是 0/150，跨平台 8.4 % / 47.6 %。结论和最终模型一致。
  - INDEL（`legacy__HG008_Illumina_INDEL_small_nofilt`）：unplaced PASS 478/534 被标，但 ctrl_tp 也有 35/145 被标。graph PoN 对 INDEL 更不特异。
- Gold set：`assembly_offref_20260929/offref_somatic.tsv`（450 个事件，233 个 exact tensor key）；d46 的 key 由 `gold_graph.py` 重新对齐。同一批 gold 的 GRCh38 PoN：329/450 被标（73.1 %）。
- "旧模型"：E9c README_methods 说 "RIS 7fc7396 SNV runs reproduce the PacBio / ONT ones exactly"，这不对。benchmark README 写的是 "SNV reproduces the 2026-10-03 runs exactly (PacBio, ONT)"，被复现的是 10-03 那批（`*_SNV_h90`），不是 09-30 / 10-01 的 pow075 runs。

**还缺什么。**

- [数据同步] HPRC v1.1（d9）的 7fc7396 最终 run，以及 HG008T Illumina、COLO829T Illumina / fiberseq 的 d46 run，都还没做 E9a / E9c / E9d。它们的 `vcf_chr1/` 和 predictions 在 RIS 上，要同步，或者在 tequila GPU 上 re-predict。之后在本地 CPU 重跑 `e9c_run.py` + `e9d_asm.py` + `collect.py`。d9 的 unplaced PASS（116 / 211 / 1,409，E3 §2）现在还没有 E9a / E9c / E9d 验证；Response 第 1 点的 d9 数是 E1.1 的 tensor 层面数，不是这三项。
- [用户决定] off-ref rescue 是否保留为默认（代价和收益见 R3.4）。组装验证已经完成（三个 run、chr1）：rescue 产出的 call 在 HG008T 上 tumor-only 只有 1 / 3 个（population 优先；只看组装是 9 / 22 个，其中 8 / 19 个也被 graph PoN 标），chr1 上叫出的 gold（3 / 5 个）也全部会被 graph PoN 扣掉。
- [用户决定] graph PoN 的定位。它在最终模型上 ctrl_tp 只标 0–1/150，但也标掉 chr1 gold 9/11、7/9（全染色体 90/119、72/97）。建议只用来描述 "population allele"，在 complete-output 指标里一律按 FP 计。
- [用户决定] COLO829T 只能做 germline 排除。要不要下载 SMaHT 的 COLO829T 组装（详版 E9d 第 4 步）？
- [用户决定] E1.1 §5 的 d46 版（E1.1 README 列为 "Optional, user's decision"）。E9c §5 已经给出 chr1 上 PASS 层面的数，全常染色体的打分还没做；要做的话放本地 tequila GPU。
- [GPU RIS] E9f：用组装给 off-ref tensor 打标签（详版 R1.8 E9f），重训 HG008T 三套，并去掉 rescue。E9f 的优先级（P2 → P1？）也要用户定。
- [文字] 删掉原稿里 "signatures / read examples" 那类间接论证，改成上面五点。Supplementary Table 用 E9a 每个 run 的 `projected.tsv`（`reviewer_experiments_20261013/E9a_unplaced_projection/runs/<run>/vcf_chr1/`）加 E9c 的 `per_call.tsv.gz` 拼出来：node、offset、class、GRCh38 区间、population / 组装类别。Discussion 里加一句 limitation："the classifier suppresses off-reference somatic calls; the rescue is a stop-gap"。
- [文字] E9c README_methods 里 "RIS 7fc7396 SNV runs reproduce the PacBio / ONT ones exactly" 要改（见上）。

### R3.4 Metrics exclude graph-coordinate calls; "about 2%"; CHM13

**状态：部分解决**: 三个问题都有最终 run 上的证据了。(1) complete-output 的 worst-case 指标（unplaced 全算 FP，36 个 run）和更正后的 unplaced 比例回答了 "about 2%" 和 "metrics cover a subset"。(2) 三个 d46 最终 run（HG008T PacBio、HG008T ONT、COLO829T ONT；chr1）的投影（E9a）：只有 2.8–7.9 % 的 unplaced PASS 能写成 GRCh38 allele，且 PoN 后不多出 TP；class ii（只投到区间的）没有一个能对上 rtg 漏掉的 truth；投影后的 complete-output F1 是 0.537 / 0.359 / 0.373，和 worst case 0.533 / 0.355 / 0.369 几乎一样。组装验证（E9d，R3.3）也说明这些 call 基本不是 somatic：HG008T 上只看组装 tumor-only 是 3.8 / 1.7 %，其中多数也是 population allele（COLO829T 只随机处理了 3,000 / 4,569，且没有 tumor 组装）。这些都只有三个 run、只有 chr1。(3) "是否 graph-specific" 由 E9a 的分类回答，CHM13 已按用户决定不做。还缺：rescue 关掉的正式数（现在 0.596 / 0.524 / 0.721 是不含 AF floor 的下界）和其余 9 个 graph SNV run 的 rescue on / off（要数据同步）；rescue 是否保留要用户定。

**Comment（要点）。**

- TP / FP / FN are computed only for calls that project to GRCh38, so precision, recall and F1 describe a subset, not the complete output.
- The graph-coordinate fraction depends on the platform and is large: above 10 %, and sometimes above 20 % for ONT (Supplementary Table 6). This contradicts the Discussion's "about 2%".
- Are graph-coordinate calls truly graph-dependent, or just not projectable to GRCh38? What fraction can be represented on T2T-CHM13, and what remains graph-specific?

**Response.**

We agree, and the revised evaluation covers the complete output.

1. *"About 2%" was wrong and has been removed.* On chr1, after population filtering and without the benchmark-region BED, graph-coordinate calls are 8.1–52.4 % of the PASS SNV output on HG008T and 17.4–42.1 % on COLO829T (HPRC v1.1), and 11.4–68.0 % on HPRC v2.1; ONT is highest everywhere. Every benchmark table now reports this fraction.
2. *Complete-output metrics.* Every run reports F1 for the GRCh38-evaluable calls and for the complete output, which counts every graph-coordinate PASS call as a false positive (a worst case). HG008T SNV, chr1, same filtering, no BED, HPRC v1.1 / v2.1 / GRCh38-only graph:

   | | Illumina | PacBio | ONT |
   |---|---|---|---|
   | evaluable calls | 0.577 / 0.561 / 0.197 | 0.614 / 0.600 / 0.359 | 0.586 / 0.563 / 0.365 |
   | complete output | 0.546 / 0.414 / 0.197 | 0.553 / 0.533 / 0.359 | 0.342 / 0.354 / 0.365 |

   Under this worst case the pangenome arms stay ahead on Illumina (both samples) and HG008T PacBio. Both fall below the GRCh38-only graph on HG008T ONT, as do v2.1 on COLO829T ONT (0.369 against 0.418) and v1.1 on COLO829T Fiber-seq (0.636 against 0.680).
3. *Projection.* In three v2.1 runs (HG008T PacBio, HG008T ONT, COLO829T ONT), the GRCh38 nodes flanking each call's bubble give a GRCh38 allele for only 19 of 240, 43 of 1,258 and 128 of 4,569, adding no true positive after population filtering; none of the 221, 984 and 1,334 placed only to an interval matches a missed truth variant. On these chr1 re-predictions, complete-output F1 moves only from 0.533 / 0.355 / 0.369 (worst case) to 0.537 / 0.359 / 0.373. On HG008T, only 3.8 and 1.7 % have tumor-only support in the assemblies (R3.3).
4. *Origin.* Every graph-coordinate SNV PASS call comes from the off-reference rescue (R3.3), which raises the haplotype path count on nodes off GRCh38. It recovers assembly-verified off-reference somatic SNVs that rtg cannot credit (HPRC v1.1, HG008T PacBio and ONT: 32 of 112 and 46 of 89 at tensor level, against 0 and 1), and it causes the ONT loss: without it, complete-output F1 of the three v2.1 runs is at least 0.596 / 0.524 / 0.721, above the GRCh38-only graph (0.359 / 0.365 / 0.418). [待补: the default (user decision)]
5. *CHM13.* We did not add a CHM13 projection: CHM13 is a path in both graphs, but our truth sets are on GRCh38, so CHM13 coordinates would not change any score. On the underlying question, only 2.8–7.9 % are GRCh38 variants that failed to project; the rest mostly lie in sequence absent from GRCh38, carry the GRCh38 base on a branch node, or lack anchors, most of them population or germline alleles (R3.3). A Supplementary Table lists each call's node, offset, class, GRCh38 interval and assembly class.

**证据（内部）。**

- Unplaced 比例是 12:00 那版用 `experiments/E3_reference_arms/tables.md` §2 算的（no BED，PoN 后）：U / (TP + FP + U)。PoN 碰不到 unplaced，所以 U 全算进输出。
  - d9 SNV：HG008T 116 / 1,434 = 8.1 %、211 / 1,401 = 15.1 %、1,409 / 2,687 = 52.4 %；COLO829T 425 / 2,447 = 17.4 %、1,018 / 3,063 = 33.2 %、1,675 / 3,981 = 42.1 %。
  - d46 SNV：34.5 / 16.2 / 46.7 %（HG008T）、20.4 / 11.4 / 68.0 %（COLO829T）。
  - INDEL d9：47.0 / 38.7 / 27.8 %（HG008T）、35.3 / 40.6 / 14.5 %（COLO829T）。
  - 详版写的 "HG008T PacBio 8 %、Illumina 32 %、ONT 41 %" 是 legacy 模型（`*_pow075_S_lr1e4_e20`）的数，以 7fc7396 的为准。
- Complete-output F1：`E3_reference_arms/tables.md` §2 的 "F1 placed only" / "F1 unplaced = FP" 两列。
  - HG008T SNV d9 0.577 → 0.546、0.614 → 0.553、0.586 → 0.342；d46 0.561 → 0.414、0.600 → 0.533、0.563 → 0.354；linear 0.197 / 0.359 / 0.365。
  - COLO829T ONT d9 0.522、d46 0.369、linear 0.418。
  - COLO829T fiberseq d9 0.781 → 0.636、d46 0.779 → 0.736、linear 0.680；COLO829T Illumina d9 0.711、d46 0.702、linear 0.468。
  - E3 README 写 "SNV d9 / d46 stay above linear on Illumina, PacBio and fiberseq"，但它自己的 tables.md §2 里 COLO829T fiberseq d9 是 0.636 < linear 0.680。以 tables.md 为准。
  - 这些都是 no BED 的数。BED 内的 PASS F1（E3 §1，0.585 等）是另一组数，不要混用。
- 投影（E9a 最终，新）：`experiments/E9a_unplaced_projection/results/tables.md`、`projection_by_run.tsv`、`metrics_by_run.tsv`、`summary.json`。三个 d46 run 的本地 in-sample re-predict（HG008T PacBio 是 CPU fp32，job 387135；两个 ONT 是 tequila GPU bf16，387379_1 / _2；`pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1`），顺序 HG008T PacBio / HG008T ONT / COLO829T ONT。
  - unplaced PASS 240 / 1,258 / 4,569。class i（写成 GRCh38 variant）19 (7.9 %) / 43 (3.4 %) / 128 (2.8 %)；class ii（只有 GRCh38 位置或区间）221 (92.1 %) / 984 (78.2 %) / 1,334 (29.2 %)；class iii（没有一致的 anchor 对）0 / 231 (18.4 %) / 3,107 (68.0 %)。
  - class i：exact 2 / 19 / 94，branch_base_differs 17 / 24 / 34；GRCh38 PoN 标掉其中 18 / 38 / 88 个；和已有的 placed call 合并 7 / 11 / 11 个。
  - class ii（`summary.json` 的 `details_pass`）：in_insertion 108 / 383 / 589，alt_is_grch38_base 89 / 517 / 608，near_branch_indel 24 / 84 / 116，too_long 0 / 0 / 21。`tables.md` 第二张表的 "ii INDEL context" 列原来显示 0（`collect.py` 读错 key）；23:08 已修好并重跑，现在是 24 / 84 / 116，F1 和 TP / FP / FN 重跑前后一致（23:10 核对）。
  - class iii：no_left_anchor 0 / 222 / 2,836，no_right_anchor 0 / 9 / 271。
  - F1（chr1，no BED，PoN 后，`metrics_by_run.tsv` 的 `nobed_pon`）：

    | 变体 | HG008T PacBio | HG008T ONT | COLO829T ONT |
    |---|---|---|---|
    | published（RIS） | 0.5995 | 0.5633 | 0.7369 |
    | placed（本地） | 0.5980 | 0.5633 | 0.7364 |
    | unplaced 全算 FP | 0.5327 | 0.3548 | 0.3690 |
    | + class i 交给 rtg（ii / iii 不计） | 0.5977 | 0.5620 | 0.7302 |
    | + class ii 按区间（0 TP，全 FP） | 0.5371 | 0.3853 | 0.5668 |
    | + class iii 算 FP（投影后的完整输出） | 0.5371 | 0.3588 | 0.3727 |

  - class i 交给 rtg 后 TP 不变（586 / 603 / 1,690）；PoN 前 HG008T 两个 run 各多 1 个 TP（611 → 612、631 → 632），COLO829T 不变（1,730）。
  - class ii 区间里的 truth TP：no BED 0/221、0/984、0/1,334；BED 内 0/208、0/940、0/1,334。
    - 规则（`evaluate_complete.py`）：只有 rtg 留下的 FN truth 能被认领（±1 bp）；`alt_is_grch38_base`（89 / 517 / 608）直接算 FP，不查区间；真正查过的是 132 / 467 / 726 个。所以 Response 写 "matches a missed truth variant"，不写 "区间里没有 truth"（已改）。
  - 第 3 点的 worst case 用的是本地 re-predict（0.5327 / 0.3548 / 0.3690，四舍五入 0.533 / 0.355 / 0.369）；第 2 点表里 d46 ONT 是 RIS 的 0.354。差 0.001 来自 re-predict，不是错。
  - 投影带来的那一点 F1 不是新 TP：`projected_worst` 的 FP 比 `placed_unplaced_fp` 少 18 / 38 / 89（no BED，PoN 后），主要是被 GRCh38 PoN 标掉的 class i（18 / 38 / 88）。
  - BED 内、PoN 后：unplaced 全算 FP 0.537 / 0.358 / 0.369，投影后 0.544 / 0.367 / 0.373（`tables.md` 第三张表）。COLO829T 的 BED 内和 no BED 数完全一样。
    - E9a 在 BED 内也把全部 240 / 1,258 / 4,569 个 unplaced 算 FP；E9c §2 的 BED 行只算区间落在 BED 内的（226 / 1,182 / 4,569），所以 E9c BED 的 "F1 all FP" 是 0.5401 / 0.3668 / 0.369。两者定义不同，不要混用；Response 只用 no BED。
  - unplaced 占 PoN 后 PASS 的比例，投影前 → 后：0.160 → 0.149、0.466 → 0.457、0.680 → 0.670。
  - CHM13 列全是 0（`chm13_pass` = `not_built`），因为用户取消了 CHM13。
  - 旧的 dev 结果（legacy d9，`tmp/reviewer_cpu_20261013/E9a_unplaced_projection/logs/dev.387367_{0,1}.log`）：HG008T ONT SNV 855 个里 class i 24、alt_is_grch38_base 562；HG008T Illumina INDEL 534 个里 class i 49、ii 485，放回 GRCh38 后 SNV 多 1 个 TP（582 → 583，BED 内 PoN 后），INDEL 不多 TP、多 28 个 FP。最终模型以上面为准。
- Graph PoN 之后的更正 F1（E9c §2，新，`experiments/E9c_offref_filter_validation/tables.md` 的 `rcpu__*` 行，no BED，PoN 后）：
  - placed 0.598 / 0.5633 / 0.7364；all FP 0.5327 / 0.3548 / 0.369；去掉 graph PoN 标的 unplaced、其余算 FP：0.5901 / 0.4889 / 0.4911；再把 tumor-only 算 TP（上界，没和 truth 匹配）：0.5908 / 0.4906 / 0.4911。
  - 剩下的 unplaced：26 / 326 / 2,292.1（COLO829T 只处理了随机 3,000 / 4,569，按比例外推，所以是小数）。
  - graph PoN 也标掉 chr1 gold 9/11、7/9，包括所有被 rescue 叫成 PASS 的 gold（R3.3）。所以这些数只是 "扣掉 population allele 之后" 的参考，不当作主结果。即使扣掉，F1 仍低于不用 rescue 的下界 0.596 / 0.524 / 0.721。
  - 组装分类（population 89.2 / 74.1 / 49.8 %，germline 7.9 / 19.2 / 45.4 %，tumor-only 0.4 / 0.2 %，COLO829T not in normal 2.4 %）的细节见 R3.3。
    - 这是互斥分类，按 population > germline > tumor-only > unverifiable > unsupported 的顺序；COLO829T 只处理了随机 3,000 / 4,569（seed 0）。
    - 只看组装：tumor-only 9/240 (3.8 %)、22/1,258 (1.7 %)，其中 8 / 19 个也被 graph PoN 标（E9c README §1 第二张表，`per_call.tsv.gz` 重新数过）。被标的里有 E9a 那个 truth SNV chr1:183,643,262 T>A（node 7923590）。所以 Response 第 3 点用 3.8 / 1.7 %，不用 "at most 0.4 %"（已改）。
- Rescue = unplaced（本次之前已核对）：`/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/<run>/test_chr1/*.metrics.json` 的 `offref_rescue.calls`，和 E3 §2 的 unplaced PASS 九个 run 全部相等。
  - HG008T d9：116 / 211 / 1,409。HG008T d46：754 / 244 / 1,261。COLO829T d46：516 / 258 / 4,566。
  - rescue 之前 argmax 判为 somatic 的（`calls_before`）：3 / 7 / 30、9 / 3 / 48、43 / 106 / 262。
  - 机制见 `machine_learning/pansoma_net_v2/predict.py` 的 docstring "Off-reference rescue"（site ±10 列的 ch6 抬到 floor 100，然后取 argmax）。graph_vcf 对这些 tensor 用 `offref_call`。
  - COLO829T 的 d9 joint 模型（`*_joint_7fc7396_gpu4`）不是 E3 里的 d9 run，没放进来比较。
  - E3 里 COLO829T 的 d9 SNV 是 `*_SNV_h90`（2026-10-03），metrics.json 不在本地，没核对。按设计 predict 默认 `--offref-site90 10`，每个 SNV off-ref tensor 都会被 rescue。
- Rescue 关掉的估计（12:00 之后自己算的），数据是 `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/`：
  - 读 `test_chr1/*.predictions.ndjson.gz`。数的是 `off_reference` 并且 `p_somatic` ≥ `rules.pipeline.threshold`（0.29122 / 0.23086 / 0.56769）的 tensor。`p_somatic` 是没改过 ch6 的分数。
  - 结果是 6 / 162 / 98，分母是 off-ref tensor 12,632 / 26,665 / 31,721，没有 excluded。用 rescue 的是 240 / 1,258 / 4,569，等于 `vcf_chr1/*.linear.stats.json` 的 `unplaced_pass`。
    - 23:10 独立重数了一遍（同样读 predictions 和 `rules.pipeline.threshold`），6 / 162 / 98 和 off-ref tensor 数都一致；metrics.json 的 `offref_rescue.calls` 240 / 1,258 / 4,569，`calls_before` 3 / 48 / 262。
    - 为什么是下界：没 rescue 时这些 off-ref call 全按 unplaced FP 算；其中如果有被 linear_vcf 放上 GRCh38、被 PoN 标掉或被 LowAF 去掉的，F1 只会更高。placed call 来自 on-ref tensor，不受 rescue 影响。
  - AF floor（LowAF）没有算进去，所以 6 / 162 / 98 是上界。
  - F1 = 2TP / (2TP + FP + FN + U)。TP / FP / FN 取 `vcf_chr1/eval_nobed/report.txt` 的 pon 行：586 / 672 / 116、603 / 836 / 99、1,690 / 457 / 753（和 E9a `metrics_by_run.tsv` 的 placed `nobed_pon` 一致）。
  - 用 rescue：0.5327 / 0.3548 / 0.3690。RIS 是 0.533 / 0.354 / 0.369，差别来自 re-predict 的 TP / FP 有细微不同。
  - 不用 rescue：0.5961 / 0.5237 / 0.7210。只算 placed：0.598 / 0.563 / 0.736。
  - unplaced 比例：用 rescue 时 16.0 / 46.6 / 68.0 %，不用时 0.5 / 10.1 / 4.4 %。
  - re-predict 和 RIS 对得上：eval_bed PoN F1 0.603 / 0.575 / 0.736，RIS 是 0.605 / 0.575 / 0.737（`benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/tables.md`）。
- Rescue 的收益（gold）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §5，d9 模型 `HG008_PacBio_SNV_7fc7396` / `HG008_ONT_SNV_7fc7396`，全部常染色体。
  - base（不 rescue，p ≥ t）：0/112、1/89。`site10_high`（就是 rescue）的 argmax：32/112、46/89。32 / 46 没算 AF floor，是 PASS 的上界。
  - d46 上 PASS 层面的数（新，E9c §5，只有 chr1）：有 tensor 的 gold key 里 PASS（unplaced）的是 PacBio 3/11、ONT 5/9、COLO829T 1/2，全部会被 graph PoN 标掉。
    - 数很小（HG008T 一共 8 个 PASS gold），只能当方向。COLO829T 的 2 个是 HG008-T 的 gold key，不是 COLO829T 的 somatic（其中 1 个在 COLO829BL 里是 germline）。
  - 这些 call 是 unplaced：虽然每个 key 都有 GRCh38 表示（`assembly_keys.json` 的 `grch38`），rtg 不给它们分，complete-output bound 还把它们算 FP。所以 Response 第 4 点写了 "rtg cannot credit"。
- E5b（新，`experiments/E5b_absorbed_d46/tables.md` §2）讲的是 recall 一侧：图上没有 candidate、linear 上有的 truth allele。HG008T SNV d46 上（Illumina / PacBio / ONT）"absorbed"（没有 candidate、ALT 被图拼出）的 109 / 124 / 121 个 100 % 被 PoN 标掉；"not spelled" 的标得少一些（73.2 / 82.8 / 87.4 %）。全表（两个样本、SNV / INDEL、d9 / d46）absorbed 的 PoN 标记率是 95.5–100 %，not spelled 低得多（COLO829T SNV 只有 7.8–17.5 %）。它不改变本点的结论，放在比较 linear 的那几个点里用；`linear_recovery.md` 还在等 RIS 的 rtg 输出。
- CHM13 在两张图里都是 path：`.../logs/prep_graph_d9.387183.log` 和 `prep_graph_d46.387184.log` 里都有 "vg paths -A -S CHM13"。用户 2026-10-13 取消了所有 CHM13 工作，第 5 点就是不做新 CHM13 结果的写法。

**还缺什么。**

- [CPU 本地] rescue 关掉的正式数。三个 d46 run 不用重新 predict：把 predictions 的 `offref_*` 字段去掉，重跑 graph_vcf → linear_vcf → vcfeval，这样 AF floor 也算进去了。出来后替换 Response 第 4 点的 "at least"。
- [数据同步] E3 §2 里其余 9 个 graph SNV run 的 predictions 还在 RIS 上：HG008T d9 3 个、HG008T Illumina d46、COLO829T Illumina / fiberseq d46、COLO829T d9 3 个。要从 RIS 同步，或者在 tequila GPU 上 re-predict。之后对全部 12 个 graph SNV run 都给出 "rescue on / off" 两种 complete-output F1，并跑 E9a 投影（第 3 点现在只有 3 个 d46 run）。
- [用户决定] SNV off-ref rescue 是否保留为默认。现在两边的证据都有了：
  - 代价：三个 d46 run 上 complete-output F1 0.533 / 0.355 / 0.369，关掉时至少 0.596 / 0.524 / 0.721；扣掉 graph PoN 标的也只有 0.590 / 0.489 / 0.491。
  - 收益：d9 上 gold 32/112、46/89（tensor 层面；关掉时 0/112、1/89），d46 chr1 上 PASS 的 gold 3/11、5/9。但这些 call 是 unplaced，rtg 的数字里看不到；组装验证里 rescue 的 call 在 HG008T 上只有 1 / 3 个是 tumor-only（population 优先；只看组装是 9 / 22 个，其中 8 / 19 个也被 graph PoN 标，R3.3）。证据只来自三个 run 的 chr1。
  - 也可以在报告里两种数都给。
- [GPU RIS] E9f：用组装给 off-ref tensor 打标签，重训 HG008T 三套，不用 rescue（同 R3.3）。
- [文字] E3 README "Complete-output bound" 一段要改：COLO829T fiberseq d9 0.636 < linear 0.680；另外 ONT 的劣势要写明来自 rescue。
- [文字] 原稿 Discussion 的 "about 2%" 删除；Supplementary Table 6 换成每个 run 的 unplaced 比例和 F1（evaluable / complete / complete projected / complete without rescue）。逐个 call 的 Supplementary Table 用 E9a 的 `projected.tsv` 加 E9c 的 `per_call.tsv.gz`（同 R3.3）。

### R3.5 Mutational signatures with the truth set

**状态：待解决**: SBS96 还一张都没算：Pansoma 目前只有三个 d46 SNV run（HG008T PacBio / ONT、COLO829T ONT）的 in-sample chr1 call set（tequila GPU 重新 predict，复现 RIS，含 rtg tp / fp；另有 E12a 的 16 个跨样本 call set），其余 run 要从 RIS 同步，比较器 VCF 等 user 给路径，truth 只有训练正样本的 SBS6 谱；graph 坐标 calls 的直接验证已改由 E9c / E9d 完成（三个 d46 run，见 R3.3），不再依赖 SBS96。

**Comment（要点）。**

- Figure 4a 里把 truth 的 spectrum 和 Pansoma、DeepSomatic-TO、ClairS-TO 并排画出来，而不只是当作 cosine similarity 的参照；HG008T 最需要，因为它的 profile 不单一，各 caller 的 cosine 都不高。
- 在 Supplement 里给出 SBS96 profiles，用来支持 "graph 坐标的 variants 是 biological" 的说法。

**Response.**

We agree, and Figure 4 is revised accordingly.

*Truth spectrum beside each caller.* For each dataset, Figure 4a now shows the SBS96 spectrum of the truth set next to those of Pansoma, DeepSomatic-TO and ClairS-TO (chr1, PASS after the same population filter), and the Supplement gives every call set split into true and false positives, with the cosine similarity of each to the truth spectrum. [待补: SBS96 panels and cosine similarities]

*HG008T.* The HG008T truth spectrum is indeed heterogeneous. Among the HG008T somatic SNVs behind our training sets, 40.6–40.8% are C>T (29.2–29.4% of all SNVs are C>T at dipyrimidines), 18.6–18.8% C>A and 12.8–13.1% C>G, whereas COLO829T is dominated by UV-type changes (71.0–71.5% C>T, 64.6–65.1% at dipyrimidines). [待补: COSMIC v3.4 decomposition of each truth spectrum.] With a mixed truth spectrum and a post-filter SNV precision of 0.45–0.50 on HG008T chr1, the cosine similarity of a whole call set to the truth can depend strongly on the spectrum of its false positives; showing true and false positives separately makes this visible. [待补: confirm with the false-positive spectra]

*Training spectrum.* Our main models are trained on the same tumor, so their calls could follow its spectrum. A preliminary cross-sample test (earlier checkpoints, HPRC v2.1 graph, chr1) suggests the models are not tied to it, without excluding a partial effect: applied without retraining to the UV-dominated COLO829T, the HG008T SNV models reach a post-filter PASS F1 of 0.708 / 0.754 / 0.630 (Illumina / fiberseq / ONT), against 0.783 / 0.779 / 0.737 for COLO829T's own models; in the other direction the COLO829T models lose more on HG008T (0.359 / 0.557 / 0.487 against 0.570 / 0.605 / 0.575), mostly through threshold calibration (R1.7). [待补: SBS96 of these cross-sample calls, true and false positives, against the target's truth spectrum]

*Graph-coordinate calls.* Their SBS96 profiles are shown separately, with the trinucleotide context taken from the graph nodes. We no longer use the spectrum as evidence that these calls are biological [用户决定: confirm this withdrawal], because a spectrum cannot distinguish somatic mutations from population alleles that GRCh38 lacks. Instead, these calls are tested directly: placement on GRCh38 through the bubble anchors, a graph-native population filter against HPRC haplotypes, and validation against the tumor and normal assemblies (R3.3, R3.4). The Supplement gives their SBS96 before and after these steps. [待补: SBS96 of graph-coordinate calls before and after each step; share confirmed somatic]

**证据（内部）。**

- truth 谱（只有 SBS6 的 context 计数，没有 signature fit）：`analysis/reviewer_response_20261011/experiments/E21_training_composition/tables.md` T7 和 `somatic_spectrum.tsv`：HG008T 训练正样本 SNV 7,098–7,274 个，C>T 40.6–40.8 %，UV-type C>T 29.2–29.4 %，C>A 18.6–18.8 %，C>G 12.8–13.1 %，TCW C>G 8.3–8.5 %，CpG C>T 9.3–9.6 %；COLO829T 33,472–33,782 个，C>T 71.0–71.5 %，UV-type 64.6–65.1 %，TCW C>G 0.8 %。注意这是训练集（chr2-22 的 positive tensors 背后的 truth allele），不是 chr1 全部 truth。
- HG008T chr1 SNV PoN 后 precision 0.452 / 0.497 / 0.466：`analysis/benchmark_7fc7396_20261012/hg008t_snv_indel_7fc7396_20261012/tables.md`。
- 比较器的 SBS96：`experiments/E8_comparators/tables.md` 的 "SBS96 cosine similarity" 节是空表（E8 已按 user 10-13 的决定取消；README 还写着 PENDING，已过时）。
- call sets 的位置：RIS run 的本地副本 `/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/<run>/test_chr1/` 只有 metrics，没有 VCF。tequila GPU 重新 predict 的三个 d46 SNV run 有完整的 chr1 call set：`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/vcf_chr1/`（`*.linear.vcf.gz`、`*.linear.unplaced.vcf.gz`、`eval_bed/vcfeval/pon_calls/{tp,fp,fn}.vcf.gz`）。
  - `eval_bed/report.txt`（PoN 后 PASS）：TP / FP 582 / 650、599 / 788、1690 / 457；P / R / F1 0.472 / 0.835 / 0.603、0.432 / 0.859 / 0.575、0.787 / 0.692 / 0.736（RIS 0.605 / 0.575 / 0.737，`benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/tables.md`）。
  - 没有 GRCh38 位置、不进 rtg 的 PASS calls（graph 坐标 calls）：240、1,258、4,569。
  - E12a 的 16 个跨样本 run（`reviewer_cpu_20261013/xs__<model>__<set>/vcf_chr1/`）也有 VCF 和 rtg tp / fp。
- 跨样本：`experiments/E12a_cross_sample/tables.md`：HG008T d46 SNV 模型 → COLO829T chr1（BED 内，PoN 后，模型自己的 R0.9 阈值）PASS F1 0.708 / 0.754 / 0.630（Illumina / PacBio→fiberseq / ONT），COLO829T 自己的 d46 模型 0.783 / 0.779 / 0.737；E12a README 的解读："a model trained on one tumor is not tied to that tumor's signature"（HG008T d46 约 7,200 个 somatic SNV 训练 tensor，COLO829T 约 34,000）。反方向 COLO829T d46 模型 → HG008T chr1：PASS F1 0.359 / 0.557 / 0.487，HG008T 自己 0.570 / 0.605 / 0.575；best F1（test-tuned）0.554 / 0.597 / 0.579，也低于 in-sample best 0.613 / 0.622 / 0.604，所以排序也有一部分损失，不全是阈值。旧 checkpoint（BQ / MAPQ 用训练 tumor 的统计量），preliminary。F1 只是间接证据：它不能区分 "不依赖谱" 和 "其他原因造成的损失"，所以 Response 写成 "suggests ... without excluding a partial effect"，直接证据要等跨样本 calls 的 SBS96（TP / FP 分开）。
- graph 坐标 SNV calls 的来源（背景，不写进 Response）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §5，HG008T 的 assembly 验证 off-reference somatic SNV：不改 ch6 时 p_somatic ≥ t 只有 0 / 112（PacBio d9）和 1 / 89（ONT d9），predict 的 off-reference rescue（site ±10 抬到 H，argmax，= `site10_high`）argmax somatic 32 / 112、46 / 89；`E9a_unplaced_projection/README.md`：unplaced SNV calls 大多是 rescue 打出来的。所以 graph 坐标 calls 是不是 somatic 要靠 E9c / E9d 判断，不能靠谱。
- 工具：`/wanglab/jshen/anaconda3/bin/python` 里装了 SigProfilerAssignment 0.2.5（带 `COSMIC_v3.4_SBS_GRCh38.txt`）和 SigProfilerMatrixGenerator 1.3.5；MatrixGenerator 没装 GRCh38 genome（`references/chromosomes/` 下只有 `tsb_BED`、`transcripts`、`context_distributions`、`exome`，没有 `tsb`；10-13 本文件核对），所以三核苷酸 context 要自己从 FASTA 取，再把 96 维矩阵交给 SigProfilerAssignment。
- graph 坐标 calls 多数是 population allele 的初步证据（旧模型，不写进 Response）：`experiments/E9c_offref_filter_validation/tables.md` §1，legacy run `HG008_Illumina_SNV_pow075_S_lr1e4_e20`：595 个 unplaced PASS 里 graph PoN 标 469（78.8 %）；`HG008_PacBio_SNV_pow075_S_lr1e4_e20`：105 里 88（83.8 %）；组装分类 not checked（E9d 未完成）。详版 R3.3 另记录：没有 truth match 的 off-reference candidate 里只有 1.2–1.6 % 是 tumor-specific。
- 设计：详版 R3.5 / E18（graph 坐标 calls 的 context 用 ALT reads consensus 或图上的 walk；SigProfilerAssignment，COSMIC v3.4）。

**还缺什么。**

- [数据同步] 其余 run 的 chr1 call sets（d9 主表的 run，d46 的 HG008T Illumina、COLO829T Illumina / fiberseq；`vcf_chr1/*.linear.vcf.gz`、`*.linear.unplaced.vcf.gz` 和 rtg 的 tp / fp）从 RIS 同步。
- [数据同步] user 已有的 ClairS-TO / DeepSomatic VCF 路径（同时用于 R3.m2 的 re-score）。
- [CPU 本地] E18 SBS96 脚本：truth 和 GRCh38 calls 用 FASTA 取三核苷酸 context，graph 坐标 calls 从 graph index 的 node 序列取 context；按 TP / FP / graph-only 分组；对 truth 算 cosine；SigProfilerAssignment 0.2.5（COSMIC v3.4，已安装）分解 truth 谱。现在就可以在三个 d46 run 和 E12a 的跨样本 call sets 上先做（Response 里 "Training spectrum" 的 [待补]）。
- [CPU 本地] E9c / E9d 在三个 d46 run 上正在跑（graph PoN + 组装验证），结果出来后才能给 graph 坐标 calls 在每一步前后的 SBS96 和 confirmed somatic 的比例。
- [用户决定] 是否撤回 "SBS96 说明 graph 坐标 variants 是 biological" 的说法，改为直接验证（E9c legacy 结果显示多数是 population allele）。
- [文字] Figure 4 改版（truth 并排、TP / FP 分开）；Supplementary SBS96 图。

### R3.6 More perfect alignment is not more correct alignment

**状态：待解决**: 关键证据——已知来源的比对正确性（E11b 模拟 reads）——还没有任何数字：三个平台的模拟都已完成（ONT 11:57 完成），但 alignment job 387342–387347 仍 held；个人基因组对照（E11c）没开始也没定；CHM13 按用户决定不做（写成 limitation）。目前只有 MAPQ 分布（HG008T PacBio 四个参考）、同一批 reads 的 `vg stats` 和 E3 的 end-to-end 结果，够改措辞，不够回答 "correct"。今天的新结果（E1.1、E12a、d46 re-predictions）都不涉及比对正确性。

**Comment（要点）。**

- 更多 perfect alignment、更少 mismatch / soft clip 不等于比对更正确：图里可匹配的序列本来就更多。
- 比较 GRCh38、T2T-CHM13 和完整 HPRC graph 上的比对质量，把 "参考更完整" 和 "多 haplotype 建模" 两个因素分开。
- 如果可行：MAPQ 分布、phased-haplotype 一致性、基于组装的 truth、segmental duplications 和其他困难区域。

**Response.**

We agree and removed the inference from "more perfect" to "more correct".

1. *What the gain is.* For the same HG008T Illumina reads, 76.86% align without any edit on HPRC v1.1 against 67.54% on the GRCh38-only graph, while the number of MAPQ 60 reads is almost unchanged (2,626,084,493 vs 2,628,588,329). On long reads, visits to branch nodes are about as often edit-free as visits to GRCh38 nodes of the same length (HG008T PacBio, 1 bp: 98.5% vs 97.9%). The graph removes edits at alleles that the sample shares with HPRC haplotypes; this is not evidence of more correct placement.

2. *MAPQ on four references, same reads.* For 14,014,144 HG008T PacBio reads, MAPQ is 60 for 94.84% (GRCh38 BAM, pbmm2), 95.88% (GRCh38-only graph), 95.88% (HPRC v1.1) and 83.32% (HPRC v2.1) of mapped reads; MAPQ scales are aligner-specific. HPRC v2.1 aligns nearly every read (0.01% unmapped against 2.64% on v1.1) but with lower MAPQ; it also uses a newer giraffe, so the causes are not separated. [待补: all platforms, per SMaHT region class and segmental duplications (E11a)]

3. *Correctness with known origin.* We simulated Illumina (2.5 million pairs), HiFi-like and ONT-like reads (40,000 each) from the HG008-N v6.2 diploid assembly, which is in neither HPRC graph, over about 285 Mb of each haplotype on chr1 and chr20 where the assembly aligns one-to-one to GRCh38. On the GRCh38-only graph, HPRC v1.1 and v2.1 we measure misplacement, "perfect but wrong" alignments and MAPQ calibration, by region class, repeat content and off-reference node use. [待补: E11b results] The truth excludes regions without a one-to-one alignment, so the rates are a lower bound in the hardest regions.

4. *Personal genome and phased haplotypes.* [待补 / 待定: graph placement of HG008-T reads in difficult regions against their placement on the HG008-N assembly, and phased-haplotype concordance (E11c)]

5. *Completeness versus multi-haplotype modelling.* Our controlled comparison is the GRCh38-only graph against the HPRC graphs with identical reads, pipeline and labels (giraffe versions differ): after population filtering, chr1 SNV F1 in the benchmark BED improves by +0.038 to +0.319 on all six data sets (R3.1). The HPRC graphs add both sequence and haplotypes, so this does not separate the two. We did not add a T2T-CHM13 arm. CHM13 is one of the 90 / 464 haplotypes of both HPRC graphs, but the graphs we use keep only nodes carried by about 10% of haplotypes, so sequence present only in CHM13 is not in them. A CHM13 calling arm would need the GRCh38 truth sets lifted over, and our simulation covers only sequence that GRCh38 resolves, where CHM13's extra sequence is not tested. We state this as a limitation.

**证据（内部）。**

- 第 1 条：`vg stats` 文件（路径见 points/R2.7.md）：d9 perfect 2286270813 / 2974694440，`max 60 (2626084493 reads)`；GRCh38-only perfect 2009031268，`max 60 (2628588329 reads)`；百分比 76.86 / 67.54 来自 `experiments/E11a_alignment_stats/results/vg_stats_crosscheck.md`。Node 层面：`analysis/branch_node_alignments_20261002/README.md` s1 表（HG008T_PacBio）整体 branch 98.0 % vs GRCh38 92.4 %，但同 README 按长度分档的表：1 bp 98.5 % vs 97.9 %，2–5 bp 96.4 % vs 97.0 %，6–32 bp 93.8 % vs 96.1 %；README 原话 "compare within a length bin: branch and GRCh38 nodes are about equally often perfect"。原稿引 98.0 vs 92.4 是节点长度造成的差别，正是审稿人批评的那种夸大，Response 已改成同长度的对比（checker 改）。
- 第 2 条：`/scratch/jshen/data/reviewer_experiments_20261013/E11a_alignment_stats/scans/HG008T_PacBio.{linear,d9,d46}.summary.json`（`complete: true`）：with_path 13,624,706 / 13,644,710 / 14,012,552，MAPQ 60 13,064,032 / 13,081,945 / 11,675,897；GRCh38 BAM `/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008_PacBio_bam_stats.txt`：unmapped 1055，MAPQ 60 13,289,473 / 14,013,089 mapped。百分比是我由计数算的（两位小数）。MAPQ 1–9 档：d46 765,525（5.46 %）vs d9 333,973（2.45 %）；30–59 档 1,427,180（10.19 %）vs 121,596（0.89 %），即 d46 把 MAPQ 60 的一部分移到了中间档，不是变成 MAPQ 0（MAPQ 0：d9 9,926（0.07 %）、d46 9,589（0.07 %）；linear 5,574（0.04 %），checker 核对）。GRCh38 BAM 是 pbmm2 的 MAPQ，和 giraffe 的 MAPQ 不是同一尺度，Response 已注明 "MAPQ scales are aligner-specific"。
- 第 2 条的版本混杂：`experiments/E3_reference_arms/README.md` caveat 4（1.66 / d9 unverified / 1.77）；E11b README 写 d9 的生产版本是 `vg_v1.65.0`，但 E3 标为 unverified，未核实。
- 第 3 条：`experiments/E11b_simulated_alignment/README.md`（Method 1–5，Caveats）；`/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/sim/segments_meta.json`：358 segments，segment_bases 570,972,967（chr1、chr20，min_aln 50,000，min_mapq 5）；`per_hap`：hap1 285,835,588、hap2 285,137,379。571 Mb 是两个 haplotype 之和（chr1 + chr20 的 GRCh38 只有约 313 Mb），原稿 "over 571 Mb of chr1 and chr20" 会被读成 GRCh38 长度，Response 已改成 "about 285 Mb of each haplotype"（checker 改）；`sim/illumina/sim_meta.json` n 2,500,000（reads 5,000,000），`sim/hifi/sim_meta.json` n 40,000（reads 40,000）；ONT（新）：`sim/ont/sim_meta.json` n 40,000，reads 40,000，`truth.npz` 和 `sim.fq.gz` 已在；log `tmp/reviewer_cpu_20261013/E11b_simulated_alignment/logs/sim_E11b_sim_387341.out` 最后是 "11:57:04 done 40000"，Exit status 0（pbsim3 `--method errhmm --errhmm .../ERRHMM-ONT-HQ.model --accuracy-mean 0.98`）。E11b README 的 Results 部分仍是 "Pending"，没有 `results/` 目录（`segment_stats.log` 是空文件）。HG008 不在两个 HPRC graph 里：E11b README Data 段，`experiments/E8b_literature_facts/README.md` H9。
- 第 5 条：E3 `tables.md` §1：d9 − linear +0.059 到 +0.311，d46 − linear +0.038 到 +0.319（六个 SNV 集，chr1，PASS，PoN 后，BED 内）。E3 README caveat 1 / 4：d9 的 HG008T PacBio / ONT SNV 是 1 GPU（linear / d46 是 4 GPU），COLO829T d9 是 10-03 的 `_h90` runs；vg 版本 1.66 / d9 unverified / 1.77。所以 Response 写 "identical reads, pipeline and labels (giraffe versions differ)"，不写 "identical code"（checker 改）。原稿第 5 条没说这个对比本身分不开 completeness 和 multi-haplotype（HPRC 同时加了序列和 haplotype），已补一句。CHM13 在两个图里：E8b H1（v1.1 = 44 个体 + GRCh38 + CHM13 = 90 haplotypes）、H4（v2.1 = 231 个样本 + GRCh38 + CHM13 = 464）；过滤：H3（`--filter 9`，"filter does not apply to the reference"）、H6（`--filter 46`，"dropping nodes and edges supported by fewer than 10% of haplotypes"）。"simulation truth covers only GRCh38-resolvable sequence" 来自 E11b README Caveats 第一条（"only where the assembly-to-GRCh38 truth is one-to-one"）和 Method 1；Caveats 最后一条也这么说，但那条是 CHM13 control 的，删 CHM13 时会一起删，所以改引第一条。
- 第 2 条的 GRCh38 BAM 比对软件（新）：`HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam` 的 `@PG`：`PN:pbmm2 VN:1.10.0 (commit v1.10.0)`，`--preset HIFI`（其余 BAM 见 points/R2.7.md）。
- CHM13：用户 2026-10-13 决定所有 CHM13 工作放弃，第 5 条按 limitation 写，不再标 [用户决定]。CHM13 工作已停：`tmp/reviewer_cpu_20261013/E11b_simulated_alignment/logs/chm13_index_387281.out`（CANCELLED 2026-10-13T11:51:35，跑了 21 min，autoindex 还在建 VG graph）；sacct：387348 `E11b_align_chm13` CANCELLED。hap → CHM13 的 PAF 已完成（387282、387283，`chm13_truth/HG008N_CHM13_dipcall.hap{1,2}.paf.gz`）。detailed doc 的 R3.6 Response 第 2 条和 E11b README 第 6 节仍写着 CHM13-only graph，需要改。

**还缺什么。**

- [CPU 本地] E11b：三个平台的模拟都齐了，由用户放出 alignment job 387342–387347（各 arm 用生产版本的 vg）和 collect 387349，`collect.sh` 填 `results/tables.md`（misplacement、"perfect but wrong"、MAPQ 校准，按 SMaHT class / 重复 / off-reference node 分层），替换第 3 条的 [待补]。
- [文字] E11b README 删掉 CHM13 臂（由 E11b 负责人改）：Reviewer points 段的 "including the CHM13-only control"、Data 表的 hap -> CHM13 truth / T2T-CHM13 / CHM13-only graph 三行、Method 第 6 节、Exact commands 里 387281 / 387282–387283 / 387348 三行、Caveats 最后一条；`collect.sh` 里的 chm13 arm 也要去掉。detailed doc 的 R3.6 Response 第 2 条同样删掉 CHM13-only graph。
- [CPU 本地] E11a：HG008T 的 window scan 已完成（387287_0、_1、_3、_4、_5），`compare.py` 对不全的 set 也能出表，现在可以先跑 HG008T_PacBio / HG008T_Illumina / HG008T_ONT，得到 MAPQ 分布、SMaHT / segdup 分层和 same-read placement concordance；held 的 387166_[5-6]、387287_[2,6-10]、387288_[12-…] 放出后补 d46 和 COLO829T。
- [用户决定] E11c（P2，HG008-T reads 在 HG008-N hap1 BAM 上的位置 vs graph 位置，只有 hap1）做不做；审稿人还点名了 "phased-haplotype concordance"，第 4 条已把它列进 [待定]。不做就把第 4 条改成一句 limitation（写明没有做 personal-genome / phased-haplotype 对照），不要整条删掉，否则这项要求没有回应。
- [用户决定] vg 版本对照（linear graph 用 vg 1.77 重建索引并重比对，E3 caveat 4）做不做；不做就在 Response 里明确写 v2.1 的 MAPQ 变化可能部分来自 giraffe 版本（第 2 条已这样写）。
- [文字] 改 detailed doc / 简版 R3.6：删掉 CHM13-only arm 的描述；"same vg giraffe version" 类的说法不要再出现。

### R3.7 Caution on reclassifying benchmark false negatives

**状态：已解决**: 审稿人要的正交证据（单倍型组装）已经有了，HG008T 和 COLO829T 的分析都已完成，结论反过来（多数是被图吞掉的真 somatic，不是 truth 错误）；剩下的是改 Supplementary Fig. 3 和正文措辞。旧稿里 "tumor assembly carries 95%" 这个数字找不到出处，已换成核实过的数字。数字都是 HPRC v1.1 d9 上的（Response 已写明）；只有主图换成 d46 时，Supp. Fig. 3 才要换成 E5b 的 d46 分类（E5b 已完成）。

**Comment（要点）。**

- 单靠 graph alignment 不能推翻已建立的 truth 注释。
- 把 Supp. Fig. 3 的 false negatives 重新归类，需要正交证据：matched-normal 的 alignments 或 haplotype-resolved assemblies。
- 否则就要把说法收弱。

**Response.**

We agree, and we withdrew the interpretation that these false negatives are reference-bias artefacts or benchmark
errors. We tested them with orthogonal evidence: haplotype-resolved assemblies of the matched normals (HG008-N v6.2;
COLO829BL) and, for HG008T, of the tumor (HG008-T v3.2), together with the HPRC haplotype walks (R1.4).

- *HG008T.* On the HPRC v1.1 graph, 2,478 truth alleles (201 SNVs, 2,277 INDELs) are lost because their reads align
  without edits to an existing graph path. Only 110 (4.4%) are carried by the matched normal, and in 109 of them the
  tumor allele is the germline allele of the other normal haplotype. For 2,085 (84.1%), neither normal haplotype carries
  the allele, but other HPRC haplotypes walk the same nodes (1,854 with exactly this allele). At 1,652 of the 1,935 such
  INDELs the tumor assembly carries exactly the truth allele over the repeat, and at 1,931 it carries a sequence found
  in neither normal haplotype. Illumina tumor reads show the expected allele at 96.2% of these INDEL loci, as at 95.7%
  of randomly chosen truth INDELs that were not absorbed.
- *COLO829T.* Of 497 such alleles, 34 (6.8%) are the normal's germline allele, and 430 (86.5%) are absent from both
  normal haplotypes but lie on nodes walked by other HPRC haplotypes (390 with exactly this allele). Without a tumor
  assembly we cannot resolve 35 further INDELs at which a normal haplotype has the ALT length with different bases, so
  the germline-like share of the absorbed INDELs lies between 7.7% and 15.8%.

The evidence therefore points the other way from our original interpretation: most of these alleles are genuine
somatic events hidden by population alleles in the graph, not truth-set errors. Supplementary Fig. 3 is reframed as
"somatic alleles absorbed by the graph", and the cost of this absorption is quantified in R1.4. Normal reads were not
realigned at these sites (no HG008-N reads were available to us); the matched-normal evidence comes from the
assemblies, whose errors in long homopolymers we can bound only indirectly (HG008-N: about 6–8% of loci, equally at
absorbed and other truth alleles).
[可选, 若做 E20: COLO829BL ONT reads aligned to the same graph take the alternative path at [待补] of the 34 germline
sites and at [待补] of the other sites.]

**证据（内部）。**

- HG008T：`analysis/graph_absorbed_somatic_20261001/README.md`
  - Short answer：union 201 SNV + 2,277 INDEL = 2,478；110（4.4 %）；109 of 110；2,085（84.1 %）。
  - §"Is the graph helping or hurting?"："Illumina reads show the expected tumor allele at 96.2% of absent INDEL loci vs 95.7%
    of random non-absorbed truth INDELs (audit table `b4`)"。
  - "Assumptions and caveats"："No HG008-N reads are on disk"；normal 组装在长 homopolymer 上的错误只能间接估计（~6-8 %）。
- 肿瘤组装：`analysis/graph_absorbed_somatic_20261001/tables.md` "HG008-T contigs vs normal haps over the array"：
  `HG008N_absent_HPRC_other` INDEL 共 1,935（NA 1、all_eq_hap1 1、no_change_both 1、no_change_hom 1、novel_ALT 1,652、
  novel_other 279）。1,652 / 1,935 = 85.4 %，1,652 + 279 = 1,931 / 1,935 = 99.8 %（百分比是我按表里的计数算的）。
- 和旧稿不一致：详版和简版写 "The HG008-T tumor assembly carries 95% of the graph-encoded HG008T INDELs"。在
  graph_absorbed 的 README / tables / report_numbers / audit 里都没有这个数；`populations.md` 里的 "1698 (95%)" 是 "AFR 有
  carrier 的位点比例"，和肿瘤组装无关。最接近的：所有 absorbed INDEL 里 `tumor_any_ALT_array` = True 1,832 / 2,267（有结果的，
  80.8 %，我从 `per_variant.tsv` 算的），absent INDEL 里 novel_ALT 85.4 %。已换成 tables.md 的计数。
  - 另一个可能的出处：`pansoma_net_v2_runs/indel_graph_paths_20260930/assembly/`（10-01，I1 = ALT 是 d9 path 的 truth
    INDEL 对 HG008-T v3.2 组装）。那里没有书面报告；按 `haplotypes_tumor.tsv` 里 "任一 tumor contig status = ALT" 数，I1 是
    1,669 / 1,803（92.6 %，我算的），也不是 95 %，所以仍不用这个数。
- COLO829T：`analysis/graph_absorbed_somatic_colo829t_20261005/README.md` Short answer：497；"Germline filtering: 34 (6.8%)"；
  "Pangenome-induced false negative: 430 (86.5%)"；"between 7.7% (33) and 15.8% (68)"；"No tumor assembly was used"。
- 本次核对改的三处（数字都在上面两个 README 里）：
  - 旧稿 "neither normal haplotype carries the allele but other HPRC haplotypes do"（2,085）和 "430 (86.5%) are alleles of
    other individuals" 说过头了：HG008T README "HPRC haplotypes walk the same nodes. ... At 1,854 of them, at least one HPRC
    haplotype carries exactly the same allele ... At 231 only the nodes are population nodes"；COLO829T README "Neither
    COLO829BL haplotype carries the ALT, but HPRC haplotypes walk the same nodes. At 390 of them ... at least one HPRC
    haplotype carries exactly the allele"。已改成 "walk the same nodes (1,854 / 390 with exactly this allele)"。
  - 加了 "On the HPRC v1.1 graph"：两个分析都在 d9 上（README 开头 "already a path of the HPRC v1.1 d9 graph"），主图可能换
    d46（见下面还缺什么）。
  - 加了 normal 组装误差的 caveat（HG008T README "Assumptions and caveats"："the tumor's retained other haplotype disagrees
    with the assembly at ~6-8% of loci, the same at absorbed and other truths"）。这是 HG008T 的估计，COLO829BL 没有对应的数。
- 核过的旁证（不进 Response）：`tumor_any_ALT_array` INDEL True = 75 + 23 + 1,652 + 82 = 1,832，有结果的 1,832 + 435 = 2,267
  （tables.md "HG008-T any contig == ALT over the array"），80.8 % 对。
- E5b（d46 吸收分类）：`experiments/E5b_absorbed_d46/tables.md` 开头 "PRELIMINARY: graph spelling missing for ['d9', 'd46'];
  "absorbed" columns are 0"，`s5_prelim.stderr.txt` 也是 "no spell_d9.tsv.gz yet" / "no spell_d46.tsv.gz yet"；graph spelling
  （`s3_extract.sh` / `s4_spell.py`）正在跑。本条 Response 的数字都来自 d9 上已完成的两个 graph_absorbed 分析，不受影响。
- PoN 背景（R1.4 用）：`tmp/graph_absorbed_pon_20261012/README.md`：absorbed truth 被现行 PoN 标掉 HG008T SNV 199/201、INDEL
  2,156/2,277 = 94.7 %，COLO829T INDEL 417/430 = 97.0 %。

**还缺什么。**

- [文字] Supplementary Fig. 3 改版（按 present / absent-other-HPRC / ambiguous 分组），正文删掉 reference-bias 的解释；把详版、
  简版 R3.7 和 R2.3 里的 "95%" 改成上面的数字。
- [用户决定] 可选 E20：COLO829BL ONT reads 在 497 个位点重新比对到 d9（约 1 天 CPU）；HG008-N reads 要下载才能做 HG008T 那一半。
- [CPU 本地] 如果主图换成 HPRC v2.1 d46（详版第 7 节用户决定 1），Supp. Fig. 3 要用 d46 的吸收分类：E5b 的 graph spelling 正在跑，
  今天晚些出结果；现在的 `experiments/E5b_absorbed_d46/tables.md` 仍标着 PRELIMINARY（"graph spelling missing for ['d9',
  'd46']"，absorbed 列全是 0），不能引用。

### R3.8 Scope: tumor-only versus somatic mosaicism

**状态：已解决**: 选择审稿人给的第一条路（全文统一为 tumor-only cancer calling），不做 SMaHT 组织的 mosaicism 验证；理由（candidate 门槛 AF ≥ 0.08、≥ 3 ALT reads，训练正样本几乎都是克隆性变异）都有现成数字，剩下的只是改标题、摘要、结果和讨论的措辞；最后一段提到的新增细胞系只是交叉引用（R1.7 / R1.m1），名单等那两条定了再填。

**Comment（要点）。**

- 验证以肿瘤为中心（两个肿瘤加一个合成 admixture），没有正常组织 mosaicism 数据集（如 SMaHT tissues）。
- 二选一：要么全文一致地把范围限定在 tumor-only calling，要么在独立表征过的 mosaicism 数据集上验证。

**Response.**

We keep the scope to tumor-only somatic calling in cancer samples. The title, abstract, Results and Discussion now state
this consistently, and broader claims about "somatic variation" in general were removed.

The method is not built for mosaicism in normal tissue, and the revised Methods say why:

- A candidate needs at least three ALT reads and an allele fraction of at least 0.08 for SNVs and INDELs, so variants at
  the much lower allele fractions typical of mosaicism in normal tissue are not evaluated.
- The training positives are tumor variants, mostly clonal: the somatic SNV training tensors of the single-sample
  models have a median allele fraction of 0.53–0.54 (HG008T) and 0.46–0.51 (COLO829T), and only 0.4–1.2% (HG008T) and
  6.5–7.3% (COLO829T) lie below 0.20.
- The classifier is trained to separate a tumor's somatic variants from the germline variants and sequencing errors of
  the same sample.

Validating on an independently characterised mosaicism dataset such as the SMaHT tissues would need a low-allele-fraction
candidate configuration, new training labels and new models. We describe this in the Discussion as future work rather
than claim it.

Within the cancer scope we broadened the validation beyond HG008T and COLO829T with further cancer cell lines (R1.7,
R1.m1) [待补: the lines with final results; so far only the breast line HCC1395, in a preliminary cross-sample test],
and the HapMap mixture, whose positives are germline alleles of other individuals, is no longer used as evidence of
somatic performance (R1.5).

**证据（内部）。**

- candidate 门槛：`indexed_gam_pipeline_v4/README.md` §4 `run build`（`--min-variants 3`、`--min-allele-bq 10`、
  `--snv-min-af` / `--indel-min-af`）；`/scratch/jshen/data/pansoma_v2_tensors/README.txt`："All sets: build AF 0.08 for SNV
  and INDEL"。
- 训练正样本的 AF：`experiments/E21_training_composition/README.md` Summary "VAF"：median AF 0.53-0.54（HG008T）、0.46-0.51
  （COLO829T）；AF < 0.20 为 0.4-1.2 %（HG008T）、6.5-7.3 %（COLO829T）。同一段还说 germline SNV tensor 的 AF 分布几乎一样，
  "AF alone does not separate the two classes"，所以 Response 不说模型靠 AF 区分。这些数是单样本模型（HG008T、COLO829T
  P100 set）的训练 tensor（E21 §3 表）；COLO829T 的 purity-joint 模型还用了 P75 / P50 / P25 稀释 set（E21："the COLO829T
  dilution sets are the only source of them [low-VAF positives] (joint models)"；Illumina / ONT 稀释有 read-source BQ
  artefact），所以 Response 加了 "of the single-sample models"。候选门槛 AF ≥ 0.08 对所有 set 都一样，不受影响。
- 新增的四个细胞系：`/scratch/jshen/data/CASTLE/README.txt`（HCC1395 乳腺 TNBC；HCC1954 乳腺；H1437、H2009 肺腺癌）。
  - 已 finalize（`pansoma_v2_tensors/<SET>/tensors/SNV/labels.manifest.json` 存在，2026-10-13 15:40 查）：HCC1395T
    Illumina / PacBio / ONT，HCC1954T PacBio / ONT，H1437T PacBio，H2009T PacBio / ONT。还没有：HCC1954T 和 H1437T 的
    Illumina（set 目录都还没有）、H1437T ONT、H2009T Illumina（有 set 目录，没有 manifest）。
  - 有 chr1 结果的只有 HCC1395T，而且是跨样本的 preliminary 结果（`experiments/E12a_cross_sample/tables.md`：HG008T 的 d46
    SNV 模型 PASS F1 0.790 / 0.698 / 0.689，Illumina / PacBio / ONT）；没有 HCC1395T 的 in-sample 模型，其余三个细胞系还没有
    模型。数字属于 R1.7 / R1.m1，这里不给。
  - 旧稿写 "further cancer cell lines (breast and lung adenocarcinoma) were added ... the validation is broader than
    before"：两个肺腺癌细胞系（H1437、H2009）现在既没有模型也没有任何 chr1 结果，HCC1395T 也只有跨样本、旧 normalization
    checkpoint 的 preliminary 结果（E12a Caveats 1、5），这样写会被读成已经验证过。已改成 "broadened ... with further cancer
    cell lines" 加 [待补]（写明目前只有 HCC1395 的 preliminary 跨样本测试），等 R1.7 / R1.m1 定了再填细胞系名单。
- HapMap positives 是其他 donor 的 germline allele：E8b M1。
- "mosaicism 的 AF 常 < 1 %" 只是详版的内部判断，没有文件出处，所以 Response 只写 "much lower allele fractions"。

**还缺什么。**

- [文字] 标题、摘要、Results、Discussion 统一为 tumor-only cancer calling；删掉泛指 "somatic variation" 的句子；Discussion 加
  mosaicism 作为 future work。
- [用户决定] HapMap 去留（R1.5）：Response 最后一句两种情况都成立，但要和 R1.5 / R2.1 的最终措辞一致。
- [文字] 最后一段的细胞系那句要和 R1.7 / R1.m1 的最终状态一致：如果最后只评估了部分细胞系，就写出名字（现在只有 HCC1395T 有
  chr1 结果）。

### R3.9 Model stability and statistical uncertainty

**状态：待解决**: 审稿人要的核心证据仍然没有：主表的每个模型都只有一个 seed，没有 bootstrap 区间，比较器差值的区间还要等比较器结果路径。今天新增的都是间接或 preliminary 的：推理在硬件 / 精度之间稳定（E1.1 §2：18 个 checkpoint 本地重跑（4 个 CPU fp32、14 个 tequila GPU bf16）对 RIS GPU bf16 的 somatic AP 差 ≤ 0.005；本地重预测三个 d46 run（一个 CPU、两个 GPU）的 PoN F1 与 RIS 差 ≤ 0.002）；"独立 HG008T 测试" 有了一个 seed 的 preliminary 数字（E12a，COLO829T 模型 -> HG008T，旧 normalization）；阈值敏感性多了跨肿瘤的一面（E12a：COLO829T 模型的阈值迁移不过去）。另有同 seed 完全可复现、改 global batch 让 F1 变 0.007、早先 INDEL 三 seed 实验和 5-Mb block bootstrap 工具。

**Comment（要点）。**

- 结果都是点估计；用多 seed 重训、基因组区间的 block bootstrap 置信区间或染色体间方差评估稳定性。
- 特别是：多 seed 重训，报告 held-out 染色体和独立 HG008T benchmark 上 precision / recall / F1 的分布。
- 用染色体或 block 级重抽样给 Pansoma 与比较器之间差值的置信区间。
- 报告对固定 0.5 阈值的敏感性。

**Response.**

We now report uncertainty for every headline number.

1. *Repeated training.* Training is deterministic for a given seed and data set (retraining two models with the same seed and settings reproduced their chr1 calls exactly), and inference is stable across hardware: repeating the inference of 18 checkpoints on a second cluster (4 on CPUs in fp32, 14 on GPUs in bf16) reproduced their somatic average precision within 0.005, and repeating three models end to end (one on CPU) gave chr1 F1 after population filtering of 0.603, 0.575 and 0.736, against 0.605, 0.575 and 0.737. The variance we report therefore comes from the seed, which also redraws the validation blocks. Each headline model is retrained with [待定: 3 / 5] seeds, and we report the mean and range of precision, recall and F1 on chr1 and on held-out cell lines, including HG008T when it is held out from training. This replaces the original COLO829T-to-HG008T benchmark; in a first single-seed test, with models that standardised base and mapping qualities with the training tumor's statistics, COLO829T models reach PASS F1 0.359, 0.557 and 0.487 on HG008T chr1 (Illumina, Fiber-seq model on PacBio HiFi, ONT), against 0.570, 0.605 and 0.575 for HG008T's own models (R1.7). [待补: E16, E12]
2. *Block bootstrap.* 95 % intervals come from resampling 1-Mb blocks of chr1 (1,000 replicates) over rtg's true-positive, false-positive and false-negative records, for each metric and for the paired difference to each comparator and between reference arms, using the same blocks for both members of a pair. [待补] The test sets are small (697 somatic SNVs on HG008T chr1; 131 somatic INDELs on COLO829T chr1), and we state where differences, for example between the two HPRC graphs, fall within these intervals.
3. *Operating point.* The fixed 0.5 threshold is no longer used. PASS is the threshold that reaches 90 % recall on held-out validation blocks, never tuned on the test chromosome. For SNVs this rule lands within 0.002–0.043 F1 of the best point of the chr1 precision–recall curve in all twelve pangenome models tested on their own tumor; across tumors it transfers less well (the COLO829T Illumina model reaches recall 0.157 on HCC1395 at its own threshold, while its curve peaks at F1 0.810; R1.7). For INDELs it lands far below the best point on long reads (HPRC v1.1, HG008T PacBio 0.089 against 0.260; ONT 0.015 against 0.140). We therefore show full curves for every model, report the test-tuned best point only as a labelled upper bound, and add a validation-chosen INDEL rule [待补]. Metrics are also given at validation recall 0.95, at the best validation F1 threshold and at argmax [待补].
4. *Chromosome-level variance.* [待定, 可选: a three-fold cross-validation over chromosome groups for one dataset.]

**证据（内部）。**

- 可复现（两个模型，所以 Response 写 "retraining two models"）：`analysis/benchmark_7fc7396_20261012/hg008t_snv_indel_7fc7396_20261012/README.md`（PacBio / ONT SNV 重跑与 2026-10-03 的 run 完全一样：best epoch 7 / 8、PASS 阈值 0.3141991198 / 0.4492301941、rtg TP / FP、unplaced 211 / 1,409；"training is deterministic for a given seed and set"）。
- 推理在硬件 / 精度之间稳定（新）：`experiments/E1.1_ch6_counterfactual/results/tables.md` §2（18 个 checkpoint：d9 HG008T 3、d46 6、COLO829T d9 joint 3、linear 6；同一批 chr1 test tensors、同一阈值，本地 base mode vs RIS GPU bf16 `predict`）。§2 的标题写 "CPU fp32 here"，但 E1.1 README（"How the runs were done" 第 3 条）说明只有 4 个 run 的 base 是 CPU fp32（`HG008_PacBio`、`HG008_ONT`、`HG008_PacBio_d46`、`HG008_ONT_d46`），其余 14 个是 tequila GPU bf16；原稿 Response 的 "CPU (fp32) and GPU (bf16) predictions of 18 checkpoints" 已改成 "4 on CPUs in fp32, 14 on GPUs in bf16"。somatic AP 差我按表相减：最大 0.0047（`COLO829_fiberseq_linear` 0.2832 / 0.2879），其次 0.0045（`COLO829_Illumina_joint` 0.6926 / 0.6971）、0.0042（`HG008_ONT_linear` 0.0202 / 0.0160），这三个都是 GPU vs GPU、而且是 class-stratified 抽样的加权估计；4 个 CPU fp32 run 的差 0.0002–0.0019（最大 `HG008_ONT_d46` 0.3264 / 0.3283）；truth TP 最多差 3（`HG008_PacBio_d46` 606 / 609）。E1.1 README 另说 per-tensor PoN 的 F1 只是估计，不等于 rtg F1（GRCh38-only 套 FP 对不上），所以 Response 只用 AP。
- 本地重预测复现 RIS（新；HG008T PacBio 是 CPU fp32，tsingtao job 387135，`tmp/reviewer_cpu_20261013/logs/rcpu_test_HG008_PacBio_SNV_d46-387135.out`；HG008T ONT、COLO829T ONT 是 tequila GPU bf16，387379 task 1–2；原稿写 "on a second GPU cluster"，已改）：`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/vcf_chr1/eval_bed/report.txt`：PoN TP / FP、F1 582 / 650 0.603、599 / 788 0.575、1690 / 457 0.736；RIS（`analysis/benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/tables.md` chr1 details）585 / 652 0.605、599 / 788 0.575、1690 / 454 0.737。
- 改 seed 会移动 validation blocks：E21 README §6（"Changing `--seed` would also move the validation blocks (`data.block_split` uses it)"）。
- 训练细节带来的变化：`analysis/benchmark_7fc7396_20261012/hg008t_snv_indel_7fc7396_20261012/README.md`（"SNV Illumina on 4 GPUs" 段），HG008T Illumina SNV d9 1 GPU（global batch 1,024）PoN F1 0.592（best 0.612）vs 4 GPU（4,096）0.585（0.607）；E3 README caveat 1："The batch change alone cost HG008T Illumina d9 0.007"。
- SNV PASS 与最好点的差：`analysis/benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012/tables.md` 第一张表（PASS F1 (best F1)），我按表里数字相减：HG008T d9 0.022 / 0.020 / 0.024，d46 0.043 / 0.017 / 0.029；COLO829T d9 0.002 / 0.003 / 0.011，d46 0.003 / 0.004 / 0.002 -> 0.002–0.043。
- INDEL PASS 与最好点：hg008t README（HG008T d9 PacBio 0.089 vs 0.260，ONT 0.015 vs 0.140）；benchmark README：INDEL 现在按 chr1 曲线最好点比较，是 upper bound。用户今天确认：held-out 规则（validation 选阈值）要等 RIS 预测同步。
- 早先的 validation 选阈值规则：`machine_learning/pansoma_net_v2/README.md` INDEL recipe 段（validation best-F1 threshold t：HG008T PacBio 0.214、ONT 0.101 vs recall-0.9 的 0.132、0.015；Illumina 0.194 vs 0.239），来自更早的 `*_INDEL_small_*` runs，不是 7fc7396 的模型。
- 四种规则的 per-tensor 结果（chr1，PoN 前）：本地 CPU 重预测 `HG008_PacBio_SNV_d46_7fc7396`（`tmp/reviewer_cpu_20261013/logs/rcpu_test_HG008_PacBio_SNV_d46-387135.out`、`.../test_chr1/HG008T_PacBio_d46.tensors.SNV.metrics.json` `rules`）：threshold_t（p ≥ 0.782）F1 0.420（vs truth 0.411），argmax 0.364（0.360），val_recall_0.9（p ≥ 0.291）0.315（0.312），pipeline 0.298（0.296）。说明 predict 已经报多种规则，rtg 层面的表还要做。
- 早先的 seed / bootstrap 证据：`ris/README.md` "INDEL: a general method"（45d56ad recipe，三 seed 0 / 1 / 2，ensemble + rescorer，flat PoN；5 Mb block bootstrap 95 % CI，例如 HG008T Illumina 0.275（0.236–0.311），COLO829T ONT 0.111（0.047–0.178）；single seeds 列 HG008T Illumina 0.223–0.285、COLO829T ONT 0.096–0.113）。只作内部依据（不是现在的模型和 PoN），说明 seed 间差可达约 0.06，COLO829T INDEL 区间很宽。
- Bootstrap 工具：`ris/tools/final_summary.py`（5 Mb blocks、默认 2,000 次、`random.Random(0)`，读 vcfeval 的 tp / fp / tp-baseline / fn VCF），可以改成 1 Mb、配对差值。
- chr1 truth 数：hprc_v21 tables.md "chr1 details"（HG008T SNV 697，COLO829T SNV 2,443、INDEL 131）。
- 一个 seed：E3 README caveat 2（"One seed per arm"）；d46 vs d9 "equal within noise"（benchmark README）目前没有区间支撑。
- 独立 HG008T 测试（新，E12a，preliminary）：`experiments/E12a_cross_sample/tables.md`：COLO829T d46 SNV 模型 -> HG008T chr1，PASS F1（best）0.359（0.554）/ 0.557（0.597）/ 0.487（0.579）（Illumina / PacBio（COLO829T 用 fiberseq 模型）/ ONT），HG008T in-sample 0.570（0.613）/ 0.605（0.622）/ 0.575（0.604）。旧 checkpoint（BQ / MAPQ 用训练肿瘤统计），一个 seed、没有 CI（E12a Caveats 1–2）；最终用 RIS 上 `--norm dataset` 训练的模型（用户决定）。
- 跨肿瘤的阈值敏感性（新，E12a）：COLO829T Illumina -> HCC1395T，阈值 0.640，PoN 后 P / R 0.934 / 0.157、F1 0.269，曲线 best F1 0.810；-> HG008T 0.359 vs best 0.554。E12a README："The ranking transfers; the calibration of p_somatic does not." HG008T 模型 PASS 与 best 差小（如 -> COLO829T 0.708 / 0.717）。
- 本地已有的 vcfeval 输出：`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/xs__*/vcf_chr1/eval_bed/vcfeval/pon_calls/{tp,fp,tp-baseline,fn}.vcf.gz`（E12a 16 个 cross-sample）和同目录下 `{HG008_PacBio,HG008_ONT,COLO829_ONT}_SNV_d46_7fc7396/vcf_chr1/eval_bed/vcfeval/pon_calls/`（3 个 in-sample），可以先拿来开发 bootstrap 脚本。
- 比较器：用户 2026-10-13 决定不重跑，已有结果的路径待给；`experiments/E8_comparators/` 有 `score_comparator.py`。

**还缺什么。**

- [GPU RIS] E16 多 seed：每个主模型再训 2 个（或 4 个）seed（E3 三臂 SNV、E12 主 fold），并把 E3 的 d9 matched rerun 一起做，否则臂间差值混入 batch / lr 差异。
- [数据同步] 把 RIS 上各 run 的 vcfeval tp / fp / tp-baseline / fn VCF 和 INDEL 的 validation 预测拷到本地（或在 RIS 上直接跑 bootstrap 和 held-out INDEL 规则）。
- [CPU 本地] 写 E16 bootstrap 脚本（在 `ris/tools/final_summary.py` 基础上：1-Mb blocks、配对差值、unplaced FP 按 anchor 位置分块，依赖 E9a）。本地已有 19 套 vcfeval 输出（E12a 16 个 + 3 个 in-sample 重预测），可以先在上面开发和试跑（比如 HG008T in-sample vs COLO829T -> HG008T 的配对差值），不用等 RIS 同步。
- [CPU 本地] 每个 run 补 val R0.95（graph_vcf `--target-recall`）、best validation F1 阈值、argmax 的 rtg 结果表。
- [用户决定] seed 个数（3 或 5）和哪些模型做；bootstrap block 大小（1 Mb 或沿用 5 Mb）；已有比较器结果的路径（配对 CI 必需）。
- [GPU RIS] E12 留出 HG008T 的模型（`--norm dataset`），回答 "independent HG008T benchmark"；最好每个 fold 也做多 seed，E12a 的单 seed preliminary 数字只能先放 Supplementary 或不用。
- [GPU RIS] 可选：一个数据集的三折染色体组交叉验证（3 run），做不做由用户定。

### R3.10 Low-VAF and INDEL claims

**状态：部分解决**: 以下几项已有数据支撑：HapMap 降级、用 COLO829T 真 somatic dilution 替代、SNV / INDEL 分开报告（7fc7396 benchmark 的数字；COLO829T 用的是 10-03 / 10-06 的单样本模型，INDEL 用 7fc7396 的 PoN 重新打分；INDEL 同时给 best 上界和 PASS）、INDEL 按 homopolymer 的初步分层（10-06 / 10-10 的旧 run）。还缺这几项：最终的 E17a 分层表（INS/DEL × 长度 × context × 平台，含两个比较器）、ClairS-TO / DeepSomatic 在同一协议下的 SNV / INDEL 数字（等用户给结果路径；没有它们之前，"SNV 有优势" 也只能写成有条件的）、INDEL 的 held-out 阈值规则，以及修掉 BQ cue 后的 Illumina / ONT dilution。

**Comment（要点）。**

- HapMap 的低 VAF 结果要加限定：那是 germline 细胞系，不是 somatic lineage。
- SNV 和 INDEL 分开报告。Pansoma 的优势在 SNV；INDEL 准确度一般，在部分 ONT 设置下还低于 ClairS-TO。整体结论要收窄。
- 用真正的低 VAF somatic benchmark，或者对已验证 somatic 变异做 in-silico dilution（COLO829T 或 HG008T），来补充 HapMap。
- INDEL 按 insertion / deletion、长度、序列 context（homopolymer、tandem repeat）和平台分层，好把低 VAF 生物学、平台错误模式和 INDEL 类别区分开。

**Response.**

*Low VAF.* We agree. The HapMap positives are germline alleles, so the mixture is no longer a main benchmark (R1.5). We replaced it with in-silico dilutions of COLO829T (validated somatic variants) at 75, 50 and 25% purity on three platforms (R2.4); on fiber-seq, SNV F1 after population filtering is 0.665–0.779 across purities. [待补: Illumina / ONT dilutions after removing a read-source base-quality cue; HCC1395.]

*SNV and INDEL are reported separately* in the Abstract, Results and Discussion. Any claim of an advantage is now limited to SNVs, and to settings where the re-scored comparators support it. All numbers below are chr1, in the benchmark BED, after population filtering, from models trained on chr2–22 of the same tumor; calls without a GRCh38 position are not scored (R3.4).

- SNV F1, at a threshold chosen on validation data, is 0.585 / 0.621 / 0.596 on HG008T (Illumina / PacBio / ONT) and 0.778 / 0.781 / 0.706 on COLO829T (Illumina / fiber-seq / ONT).
- INDEL F1 is far lower. At the best threshold of the chr1 curve, an upper bound, it is 0.264 / 0.260 / 0.140 on HG008T and 0.195 / 0.120 / 0.071 on COLO829T; at the validation threshold it is 0.241 / 0.089 / 0.015 and 0.115 / 0.016 / 0.003. [待补: a held-out INDEL threshold rule.]
- [待补: ClairS-TO and DeepSomatic-TO under the same protocol. We state every setting where ClairS-TO is better (the reviewer points to ONT INDELs).]

*INDEL stratification.* [待补: final table by insertion/deletion, length (1, 2–4, 5–15, 16–50 bp), repeat context and platform, for Pansoma and both comparators.] Preliminary analyses already separate the factors:

- *INDEL classes.* On HG008T chr1, 74.1% of the truth INDELs in the BED lie in homopolymers of at least 7 bp and 11.1% in STRs. Of all 635 chr1 truth INDELs, 413 match population alleles and are removed by the filter, so no caller using it can exceed a recall of 0.350 (COLO829T: 0.382). Insertions become candidates less often than deletions (HG008T 36–44% against 52–65%).
- *Platform error modes* (HG008T chr1, earlier models at the validation threshold). On Illumina, F1 is similar inside and outside homopolymers ≥ 7 bp (0.249 and 0.213). On ONT, 15,638 of 17,305 false INDEL calls lie in such homopolymers, and only 95 of the 15,733 calls there are true.
- *Low VAF.* On the long reads, the best-threshold INDEL F1 is 0.047–0.073 (fiber-seq) and 0.051–0.068 (ONT) at every purity, from 4–17 true calls, too few to resolve a purity trend.
- *Truth sets.* Short-read truth sets contain far fewer homopolymer (≥ 7 bp) INDELs than the long-read HG008T truth (SEQC2 HCC1395 15.5% against 74.9%), so INDEL results cannot be compared across truth sets.

Within these small numbers, INDEL accuracy is limited mainly by population-allele overlap, long-read homopolymer errors and, on long reads, the threshold, rather than by allele fraction.

**证据（内部）。**

- SNV / INDEL 主数字：
  - HG008T 的数字来自 `analysis/benchmark_7fc7396_20261012/hg008t_snv_indel_7fc7396_20261012/README.md` 的 Results 表：SNV PASS F1 0.585 / 0.621 / 0.596；INDEL chr1 best 0.264 / 0.260 / 0.140（PASS 0.241 / 0.089 / 0.015）。
  - COLO829T d9 单样本模型的数字来自 `analysis/benchmark_7fc7396_20261012/README.md` 的 linear 表：SNV 0.778 / 0.781 / 0.706；INDEL best 0.195 / 0.120 / 0.071。INDEL 的 PASS F1 0.115 / 0.016 / 0.003 来自 `colo829t_purity_joint_7fc7396_20261012/tables.md` 的 "P100: the joint model against the earlier single-set model" 表（和 `experiments/E8_comparators/tables.md` 一致）。核对时把两组 PASS 值加进了 Response，免得只给 test-tuned 的上界。
  - COLO829T 的 INDEL 是 `*_INDEL_small`（lr 4e-4），用 7fc7396 的 PoN 重新打分；SNV 是 `*_SNV_h90`（10-03）。
- 都是 in-sample（同一肿瘤 chr2–22 训练、chr1 测试）、in BED、PoN 后；没有 GRCh38 位置的 PASS call 不评。SNV 的 unplaced：HG008T 116 / 211 / 1,409（hg008t README Results 表），COLO829T 单样本模型 425 / 1,018 / 1,675（`colo829t_purity_joint_7fc7396_20261012/runs/COLO829_*_SNV_h90/vcf_chr1.eval_bed.txt`；全算 FP 的话 PoN 后 precision 0.710 / 0.572 / 0.421，报的是 0.860 / 0.857 / 0.727）。核对时在 Response 里加了这句范围说明（"in the benchmark BED ... from models trained on chr2–22 of the same tumor; calls without a GRCh38 position are not scored"）。
- INDEL 用 best F1 是 test-tuned 的上界：benchmark README 的 "Main results" 段和 hg008t README 都写了 "an upper bound for a threshold rule, not a held-out result"。这句在 commit 0049e36 里（history 改写前的 id 是 d4f898c，引用时用 0049e36）。
- PoN 上限：hg008t README 写着 "tags 413 of the 635 chr1 truth INDELs ... no caller recalls more than 0.350"；COLO829T 的 0.382 来自 `experiments/E8_comparators/tables.md` INDEL 表的 "PoN bound on recall" 列。
- Repeat 组成：`analysis/indel_repeat_context_predictions_20261006/README.md` 和 `tables.md`。
  - HG008T 的 chr1 BED 里，HP≥7 占 436 / 588（74.1 %），STR 占 65（11.1 %）；chr1–22 的 HP≥7 是 6,360 / 8,496（74.9 %）。
  - 其他 truth set：SEQC2 HCC1395 15.5 %，NYGC COLO829 16.3 %，SMaHT COLO829T 61.0 %。
  - 不一致：详细版 / 简版写的是 "74 % ... homopolymers of at least 5 bp"，源文件给的是 ≥ 7 bp 下 74.1 %（repeat_context 的分类）。`tmp/indel_hp_strata_20261010/counts.txt` 用 PoN 脚本的 `homopolymer_run` 算，HP≥7 是 426 / 587，HP≥5 是 432 / 587。这里改用 ≥ 7 bp 的 74.1 %。
- 按 homopolymer 分的 call 层面数字来自 `tmp/indel_hp_strata_20261010/strata.txt` 和 `counts.txt`。这些是 10-10 的 tequila runs（sqrt，lr 4e-4，阈值 val R0.9，PoN 用的是 hp 规则），不是 7fc7396 的 RIS 模型，所以 Response 里写的是 "preliminary"。
  - Illumina：HG008_Illumina_INDEL_small_nofilt 的 HP≥7 F1 是 0.249，non-HP 是 0.213。
  - ONT：HG008_ONT_INDEL_small_nofilt 过 PoN 后，HP≥7 有 15,733 个 call（TP 95，FP 15,638），non-HP 有 1,709 个（FP 1,667）。
  - FP 总数 17,305 来自 `tmp/indel_non_20261006/results_detail_20261011.md`（HG008T ONT "不加过滤, sqrt, 4e-4" 的 val R0.9 行：137 / 17305 / 450）。
- INS / DEL：`analysis/tensor_recall_20260930/README.md` 写着 "HG008T built INS 36–44 % vs DEL 52–65 %; COLO829T INS 50–56 % vs DEL 64–79 %"。这是 candidate 层面的数字，用的是 09-30 的 tensor 集。
- 长读长假 INDEL 的来源，作为补充（Response 没有展开）：`tmp/indel_non_20261006/rules.txt`，chr1–22 的 non tensor。
  - 规则 "1-bp in HP>=5: AF >= 0.12" 去掉的 non：HG008T PacBio 56.3 %，COLO829T fiberseq 56.5 %，HG008T ONT 23.3 %，COLO829T ONT 41.6 %。
  - 跨个体 error panel（AF < 0.15）：PacBio / fiberseq 的 non 少 24.0–32.5 %，truth 只丢 0–19 个（E17b 的依据）。
- Dilution 的 INDEL：`analysis/benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012/tables.md` 的 best-F1 表。
  - fiberseq 0.068 / 0.047 / 0.053 / 0.073，ONT 0.054 / 0.051 / 0.068 / 0.060。
  - README §4 写着 "the best points hold 4-17 true calls, so the purity trend is within noise"。
  - Illumina 是 0.163 / 0.308 / 0.330 / 0.233。Illumina / ONT dilution 带 read-source BQ artefact（normal reads 的 BQ 比 tumor reads 高），README §4 写的是 "the read-source BQ cue (2.) applies to INDEL reads as well, so the dilution numbers are likely optimistic"，所以 Response 没写 Illumina。fiberseq 没有这个 cue（README §2）。
- E12a 跨样本（preliminary，`experiments/E12a_cross_sample/README.md` Caveats 4）也说明 INDEL 不能跨 truth set 比：HCC1395T（SEQC2，short-read truth）的 PoN-limited INDEL ceiling 是 0.43–0.79，HG008T / COLO829T 是 0.25–0.27，README 原话 "The INDEL numbers are not comparable across targets"。Response 没有引用 E12a 的 INDEL 数字。
- 比较器：`experiments/E8_comparators/tables.md` 只有 Pansoma 的行，没有 comparator 的数字；按用户今天的决定，E8 已取消，只做 re-score。所以原稿 Response 的 "We now claim an advantage for SNVs only" 在 re-score 之前没有证据，核对时改成 "Any claim of an advantage is now limited to SNVs, and to settings where the re-scored comparators support it"；"including ONT INDELs" 改成 "(the reviewer points to ONT INDELs)"，因为我们自己还没有这组比较。
- 核对时 Response 的其他改动：
  - *INDEL classes*：74.1 % 的分母是 chr1 BED 里的 588 个，413 / 635 是整条 chr1（hg008t README），原稿把两个分母放在一起没说明，已分开写。"reach a candidate" 改成 "become candidates"（tensor_recall 的 "built"）。
  - *Platform error modes*：标明是 "HG008T chr1, earlier models at the validation threshold"（10-10 的 tequila runs，见下）。
  - *Low VAF*：标明是 best-threshold 的 F1，并写出 4–17 个 TP（README §4）。原稿 "stays at ... at every purity" 读起来像 PASS 结果。
  - *Truth sets*：原稿 "Short-read truth sets ... (SEQC2 HCC1395 15.5%, HG008T 74.9%)" 把 HG008T 放进了 short-read 的括号里；HG008T 是 HiFi + ONT + assembly 的 truth（repeat_context README），已改。
  - 结论句：原稿 "more than by allele fraction" 没提长读长 INDEL 的阈值问题；hg008t README 写 "What limits INDEL on the long reads is the decision rule"（PASS 0.089 / 0.015 对 best 0.260 / 0.140），所以加了 "the threshold"，并加了 "Within these small numbers"。
  - *Low VAF* 段加了 fiberseq SNV 的 per-purity 范围 0.665–0.779（purity tables，PoN 后 PASS），回答 reviewer 第 3 问。

**还缺什么。**

- [数据同步] 7fc7396 各个 INDEL run 的 chr1 vcfeval 输出：6 个 d9、6 个 d46、6 个 linear，再加 3 个 purity joint 的 4 个 purity；要 `vcf_chr1*/eval_bed/vcfeval/pon_calls` 下的 tp / fp / fn / tp-baseline。本地 `ris_7fc7396_20261012` 只有 `.metrics.json` 和 `best.pth`。这些是 E17a 的输入。
- [CPU 本地] E17a 分层表：在 `tmp/indel_hp_strata_20261010/strata.py` 的基础上加 INS/DEL、长度（1、2–4、5–15、16–50 bp）、context（RepeatContext：HP≥7、HP4–6、STR、other）和平台，PoN 前后都报，用最终模型。
- [用户决定] 给出 ClairS-TO / DeepSomatic 已有结果的路径。之后 [CPU 本地] 用 `E8_comparators/score_comparator.py` 在同一协议下 re-score，同样做 E17a 分层，核实 "ClairS-TO 在 ONT INDEL 上更好" 是哪些设置。
- [数据同步] 同步 INDEL 的 validation predictions，定 held-out 阈值规则，替换现在 test-tuned 的 best F1（用户已在跟进）。
- [GPU RIS] 修掉 BQ cue 以后重训 Illumina / ONT 的 dilution INDEL 模型（E14b，同 R1.5、R2.4）。
- [用户决定] 是否做 E17b 跨个体 error panel（长读长 INDEL，可选）。规则已有（`tmp/indel_non_20261006`），可以在最终模型上 [CPU 本地] 评估。
- [文字] Abstract / Discussion 里把 SNV 和 INDEL 分开写；INDEL 只说限制，不主张优势；写明 HCC1395 / CASTLE 的 truth 少了 HP / STR INDEL（rule (a)）。

### R3.m1 Parameter discrepancies between Methods and Supplement

**状态：部分解决**: candidate、tensor、训练 recipe、PoN、评估的参数现在只有一套，都能从 v4 README、pansoma_net_v2 README 和 manifest / checkpoint 核对，从 checkpoint 重跑也复现了 RIS 的 chr1 F1（三个 d46 SNV run，差 ≤ 0.002）；但最终报告用的模型还没定（INDEL 的 held-out 阈值规则待定、跨样本模型要在 RIS 上用新的 `--norm dataset` 重训、COLO829T d9 模型不完全按 recipe（1 GPU、INDEL lr 4e-4）要重跑、d9 的 vg 版本未核实），参数表要等这些定了才能定稿。

**Comment（要点）。**

- 主 Methods（AF > 10 %、2 条 supporting reads、排除 MAPQ ≤ 5、80 epochs）和 Supplementary Methods（AF ≥ 0.08、≥ 3 个 ALT
  observation、mean ALT BQ ≥ 10、MAPQ ≥ 10、INDEL ≤ 50 bp、80–120 epochs）不一致。
- 上游（candidate）和下游（VCF、过滤、评估）步骤写得比网络部分简略得多。

**Response.**

We agree. One parameter table (Supplementary Table [待补]) now serves the main and Supplementary Methods, and every value
is recorded in each run's manifest and checkpoint; re-running prediction and evaluation from three SNV checkpoints on a
second cluster reproduces the reported chr1 F1 within 0.002. The earlier values (including 80–120 epochs) belong to the
first version and are superseded. "MAPQ ≤ 5" and "MAPQ ≥ 10" were two different steps, now both stated:

- *Alignment:* vg giraffe (GRCh38-only graph: vg 1.66; HPRC v2.1: vg 1.77; HPRC v1.1: [待补: version]).
- *Node selection:* alignments with MAPQ > 5; a node is a target when more than 5% of its alignments carry an edit.
- *Candidates:* reads with MAPQ ≤ 10 are dropped; at least 3 ALT reads, each ALT observation with base quality ≥ 10
  (per observation, not a mean); allele fraction ≥ 0.08 for SNVs and INDELs; INDELs up to 50 bp; at most 800 reads per
  node (deterministic subsampling).
- *Tensors:* 8 channels × 200 rows × 101 columns.
- *Training:* chromosomes 2–22, 5% of their ~1-Mb blocks held out for validation; chromosome 1 is the test chromosome.
  SNV: 50.4 M-parameter ConvNeXt, 20 epochs, class weights pow 0.75. INDEL: 8.7 M parameters, 12 epochs, square-root
  class weights. Both: learning rate 1e-4, batch 1,024 per GPU (models trained on four GPUs: global batch 4,096; listed
  per model) [待补: the HPRC v1.1 COLO829T models of the current benchmark deviate (one GPU; INDEL lr 4e-4); reruns
  pending]. The checkpoint is the epoch with the best thresholded somatic F1 on the validation blocks (epochs 4–12 of
  20 for SNV, 2–7 of 12 for INDEL in the reported models).
- *Normalisation:* base and mapping quality are z-scored with the training data's statistics in the models reported so
  far, and per dataset in newer models (now the default, for use across samples, R1.7); the haplotype-count channel
  always with the training statistics. [待补: the setting of each final model]
- *Inference:* bf16 on GPU; another GPU, or fp32 on CPU (4 of 18 SNV models), gives per-tensor somatic average
  precision within 0.005.
- *Calls:* PASS at the validation threshold for somatic recall 0.9 [待补: INDEL operating point]; SNVs on off-reference
  nodes are called by argmax after their path counts in the 21 central columns are raised to the all-haplotype code.
- *Population filtering (R2.6):* exact-allele matches. SNVs: gnomAD or CoLoRSdb at AF ≥ 1e-4, dbSNP records not flagged
  somatic, the 1000 Genomes panel. INDELs: gnomAD or CoLoRSdb at AF ≥ 0.01 in homopolymers of ≥ 7 bp, else ≥ 0.001.
- *Evaluation:* rtg vcfeval against each truth set in its benchmark regions, and on the whole chromosome. Off-reference
  calls have no GRCh38 position; they go to a separate graph-coordinate VCF that vcfeval does not score (R3.4).

The Methods now describe candidate generation, labelling, VCF conversion, population filtering and evaluation at the
same level of detail as the network.

**证据（内部）。**

- 节点选择和 candidate：`indexed_gam_pipeline_v4/README.md` §4 `run discover`（`--min-mapq 5` exclusive；`--node-alt 0.05`
  "imperfect / (perfect + imperfect) > this"）和 `run build`（`--min-mapq 10` exclusive、`--min-variants 3`、`--min-allele-bq 10`、
  `--max-indel-len 50`、`--max-node-reads 800` 按 record SHA-256、`--rows 200 --width 101`）；build AF 0.08：
  `/scratch/jshen/data/pansoma_v2_tensors/README.txt`。
- 训练：`machine_learning/pansoma_net_v2/README.md`
  - "train on chr2-22 (5 % of their ~1 Mb node blocks validate), chr1 left out"；`--val-fraction 0.05`。
  - SNV recipe：50.4 M、`pow:0.75`、batch 1024、lr 1e-4、20 epochs；INDEL recipe：8.7 M、`sqrt`、batch 1024、lr 1e-4、12 epochs。
  - Selection："best.pth is the epoch with the highest thresholded per-tensor somatic F1 on the validation tensors"。
  - Normalization：`--norm dataset` 是 2026-10-13 起的默认；之前的 checkpoint（config 没有 `norm`）BQ / MAPQ 用训练统计量。
    对应 commit f2a0e4e（history 改写前的 id 是 cbc94a7）。旧稿 Response 写 "base and mapping quality are z-scored per
    dataset"（现在时），但现在报告的全部模型（RIS 7fc7396 checkpoint，E1.1 / E12a / E3 用的）都是旧 checkpoint，用训练统计量
    （E12a tables.md 表头 "Old checkpoints: BQ / MAPQ z-scored with the training tumor's statistics"）；已改成两种都写明，
    哪个最终模型用哪种留 [待补]。
  - 为什么参数表必须写明每个模型用哪种 normalisation：`experiments/E12a_cross_sample/README.md` "Reading"——旧 checkpoint
    （训练肿瘤的 BQ / MAPQ 统计量）跨肿瘤时 ranking 能迁移、calibration 不能：COLO829T Illumina SNV d46 模型在 HCC1395T 上
    PASS（阈值 0.640）recall 0.157、precision 0.934，而 chr1 曲线的 best F1 是 0.810（`tables.md`）。
  - Off-reference rescue：`predict --offref-site90 W`（默认 10，中间列 ± 10 = 21 列），argmax。
- DDP：`analysis/benchmark_7fc7396_20261012/README.md`："Multi-GPU runs use DDP with `--batch-size 1024` per GPU (global batch
  N x 1024)"；E3 README caveat 1：d9 HG008T PacBio / ONT SNV 1 GPU，linear / d46 4 GPU；"COLO829T d9 SNV are the 2026-10-03
  one-GPU runs (`*_SNV_h90`)"、"COLO829T d9 INDEL were trained at lr 4e-4"、"The matched d9 reruns are GPU jobs, in the RIS
  list"。所以 "Both: learning rate 1e-4" 对现在 benchmark 里的 COLO829T d9 INDEL 不成立，Response 加了 [待补]。
- best epoch：`experiments/E21_training_composition/README.md` Summary "Capacity"：SNV 4-12 of 20，INDEL 2-7 of 12。
- vg 版本：`experiments/E3_reference_arms/README.md`（linear 1.66、d46 1.77、d9 "unverified"）。
  `experiments/E11b_simulated_alignment/README.md` 和 `align_arm.sh` 对 d9 用 `/scratch/wzhang/bin/vg_v1.65.0`，说是 "the vg
  version of the production GAMs"，但没给出处（GAM header 或当时的 job 脚本），所以仍算未核实。
- PoN 规则：`tmp/graph_absorbed_pon_20261012/README.md`（SNV gnomAD/CoLoRSdb AF ≥ 1e-4、dbSNP non-somatic、1000G；INDEL
  gnomAD/CoLoRSdb exact allele，HP ≥ 7 bp 时 AF ≥ 0.01，否则 ≥ 0.001）；E8b §5 第 5 条（全部 exact allele；dbSNP 只去掉 SAO=2）。
- 评估：`ris/vcf.sh`（`eval_bed` 用 labels manifest 的 somatic BED，`eval_nobed` 全染色体）。unplaced calls：
  `pansoma_net_v2/README.md` "Events with no GRCh38 position go unchanged to `<name>.linear.unplaced.vcf.gz`, in graph
  coordinates"、"they are not evaluated"；三个复现 run 的 report.txt 也写 "240 / 1,258 / 4,569 PASS calls without a GRCh38
  position are not evaluated"。这是 R3.4 的核心问题，Methods 的评估段必须写明，Response 已加一句。
- 复现（Response 第一段新加的一句）：tequila GPU 上从 checkpoint 重跑 predict + `jobs/vcf.sh` 的三个 d46 SNV run
  （`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1/eval_bed/report.txt`，pon 行）：
  `HG008_PacBio_SNV_d46_7fc7396` 0.603、`HG008_ONT_SNV_d46_7fc7396` 0.575、`COLO829_ONT_SNV_d46_7fc7396` 0.736；RIS 的同一
  三个 run 是 0.605 / 0.575 / 0.737（`experiments/E3_reference_arms/README.md` 的 d46 列）。最大差 0.002。
- Inference 精度：`experiments/E1.1_ch6_counterfactual/results/tables.md` §2 + README "How the runs were done" 第 3 条。
  旧稿写 "fp32 on CPU ... within 0.005 (18 SNV models)"，不对：README 说 tables.md §2 的表头 "CPU fp32 here" 过时了，18 个 run
  里只有 4 个的 base 模式是 CPU fp32（`HG008_PacBio` 0.2333 / 0.2331、`HG008_ONT` 0.2316 / 0.2323、`HG008_PacBio_d46`
  0.3523 / 0.3528、`HG008_ONT_d46` 0.3264 / 0.3283，最大差 0.0019），其余 14 个是 tequila H100 NVL 的 GPU bf16 对 RIS GPU
  bf16；18 个合起来最大差 0.0047（COLO829_fiberseq_linear 0.2832 / 0.2879）。已改成 "another GPU, or fp32 on CPU (4 of 18
  SNV models) ... within 0.005"。GPU 用 bf16：E1.1 README "bf16 autocast on a GPU (as `predict` does on a GPU)"，E12a README
  "bf16 as on RIS"（选项是 `pansoma_net_v2/README.md` 的 `--amp bf16|off`）。
- 和旧稿不一致：详版写 "batch 1,024" 不分 GPU 数，实际多卡 run 是全局 4,096，已改。
- 复现那句改成 "three SNV checkpoints"：三个都是 d46 SNV run（PASS、BED 内、PoN 后、unplaced 不算），没有 INDEL 和 d9。
- Off-reference rescue 措辞："set to" 改成 "raised to"（`pansoma_net_v2/README.md`："its ch6 cells below the every-haplotype
  code raised to it in the site columns"；高于 100 的 repeat-visit code 不动，未覆盖的 cell 仍是 0）。

**还缺什么。**

- [用户决定] INDEL 的 operating point：validation 上选阈值的 held-out 规则 vs recall-0.9；test-tuned best F1 只作上界（需要从 RIS
  同步 predictions 后定）。
- [GPU RIS] 最终（跨样本，E12）模型用新默认 `--norm dataset` 训练；参数表按最终模型填 normalisation、GPU 数和 best epoch。
- [GPU RIS] E3 caveat 1 的 matched d9 reruns（COLO829T d9 SNV 4 GPU、INDEL lr 1e-4），否则参数表要逐模型写出偏离。
- [文字] 核实 d9 GAM 的 vg giraffe 版本（GAM 的 header / 当时的 job 脚本；E11b 假定是 1.65.0，没有出处），填进表里。
- [文字] Supplementary Table（从 manifest / checkpoint 自动生成）和 Methods 的上下游段落。

### R3.m2 Comparator setup

**状态：部分解决**: 两个比较器的版本、tumor-only 模型、population / PoN 资源和默认规则、训练数据和 held-out 染色体都已按一手来源核实（E8b），两种 PoN 规则的差异已有六个 truth set 的数字，评估协议和打分脚本也已就绪，同一样本训练的不对称有 E12a 的 preliminary 跨样本数字；但 user 已有那批结果实际用的版本、checkpoint、BAM / aligner 和命令行要等 user 给路径，按同一协议的 re-score 数字也还没有。

**Comment（要点）。**

- 像描述 Pansoma 一样详细地记录比较器：软件版本、aligner 和版本、预训练 checkpoint、确切的 tumor-only 配置、PASS / filter 标准、参考资源、population / PoN 资源、命令行。
- Pansoma 的主要优势是 precision，而 precision 对过滤配置很敏感。

**Response.**

We agree. A new Supplementary Table documents both comparators at the level of detail given for Pansoma. [待补: the exact versions, checkpoints, BAMs, aligners, reference and command lines of the runs reported]

*Software and models.* DeepSomatic [待补: version; latest release 1.10.0] in tumor-only mode (`--model_type` WGS_TUMOR_ONLY, PACBIO_TUMOR_ONLY or ONT_TUMOR_ONLY), and ClairS-TO [待补: version; latest release 0.5.1] with its platform models (Illumina `ilmn`, PacBio `hifi_revio`, ONT `ont_r10_dorado_sup_5khz`; [待补: `ss` or `ssrs`]). ClairS-TO has no model for the dorado sup v5 basecaller of the HG008T ONT data; [待补: the model used].

*Population resources and PASS criteria.* DeepSomatic tumor-only uses a population allele-frequency input channel (DeepVariant 1000 Genomes calls for Illumina, CoLoRSdb for long reads, AF > 0.05); its panel-of-normals filter is applied only with `--use_default_pon_filtering=true` [待补: whether set]. ClairS-TO tags non-somatic variants with four resources by default: gnomAD (AF ≥ 0.001), dbSNP 138 non-somatic sites, the GATK 1000 Genomes PoN and CoLoRSdb (AF ≥ 0.001), matching the allele for gnomAD and dbSNP and the position only for the other two. Our population filter uses the same four source files with exact-allele matching throughout (R2.6). These choices matter for INDELs: applied to the whole truth sets (all records, no benchmark regions), ClairS-TO's default rule tags 87.1% of the HG008T somatic truth INDEL alleles, 74.0% on COLO829T and 13.8% on HCC1395, against 61.8%, 52.0% and 7.3% for ours; for SNVs the two rules agree within about one percentage point on six truth sets.

*Training data.* DeepSomatic holds out chr1 in every model; the tumor-only models of its paper were trained on HCC1395, and the WGS_TUMOR_ONLY model was retrained with all of the paper's datasets from release 1.9.0. ClairS-TO's README states that all models exclude chr20, whereas its paper states chr1; its `ssrs` models were fine-tuned on HCC1937, HCC1954, H1437 and H2009. Neither tool was trained on HG008; COLO829 is a ClairS-TO benchmark sample, not a training sample, and is not named in the DeepSomatic paper. Pansoma's models in the main comparison are trained on chr2–22 of the same tumor; this asymmetry favours Pansoma and is stated beside every table. Without retraining, the HG008T SNV models reach a post-filter PASS F1 of 0.708 / 0.754 / 0.630 on COLO829T chr1 (Illumina / fiberseq / ONT; in-sample 0.783 / 0.779 / 0.737), and the COLO829T models 0.359 / 0.557 / 0.487 on HG008T (in-sample 0.570 / 0.605 / 0.575), mostly through threshold calibration (preliminary: earlier checkpoints, HPRC v2.1 graph; R1.7).

*Scoring.* Every caller is scored with the same truth, benchmark regions, chromosome and evaluation (rtg vcfeval 3.12.1, `--squash-ploidy --sample ALT,ALT`), in three settings: (1) the tool's own PASS; (2) its PASS after our population filter; (3) all its records after our population filter, ranked by its QUAL, as a precision–recall curve. Because precision depends on filtering, all three are reported. [待补: numbers]

**证据（内部）。**

- 事实表：`analysis/reviewer_response_20261011/experiments/E8b_literature_facts/README.md` §3。
  - DeepSomatic：D2 / D14（所有模型 chr2–20 训练、chr21–22 tune、chr1 held out）；D10（论文的 Illumina / PacBio / ONT tumor-only 模型用 HCC1395 SEQC2 训练，chr1 held out）；D12（tumor-only make_examples 的 population VCF：Illumina 用 DeepVariant 1KG AF > 0.05，PacBio / ONT 用 CoLoRSdb AF > 0.05）；D13（PoN filter 默认关，`--use_default_pon_filtering=true` 才开）；D15（1.10.0 是 2026-03-05 的最新 release；1.9.0 的 WGS_TUMOR_ONLY 用全部论文数据重训）；D16（1.10.0 长读长 tumor-only checkpoint 的训练数据 [unverified]）；D6 / D7（HG008 只用作 held-out 评估；COLO829 未在论文里出现）。
  - ClairS-TO：C1（v0.5.1，2026-09-09）；C3（ssrs 在 HCC1937、HCC1954、H1437、H2009 上 fine-tune，HCC1395 完全留出，COLO829 没训练）；C4（README 说全部模型排除 chr20，论文说排除 chr1，冲突未解决）；C5（COLO829 和 HCC1395 是 benchmark 样本，全文没有 HG008）；C6 / C7（四个 PoN 及规则：gnomAD AF ≥ 0.001、dbSNP 138 non-somatic、1000G PoN 全部位点、CoLoRSdb AF ≥ 0.001；gnomAD / dbSNP 要求 allele，1000G / CoLoRSdb 只看位置）；C9（我们的四个 PoN 文件就是 ClairS-TO 用的那四个源文件）；C10（没有 dorado SUP v5 模型）。
- 规则差异的后果：E8b README §4.5，HG008T 全部 somatic INDEL 9,362 个：Pansoma 规则标 5,783（61.8 %），ClairS-TO 规则标 8,155（87.1 %）；SNV 9,609 个：982（10.2 %）vs 906（9.4 %）。
  - 其余五个 truth set：Slurm 387253_1-5 已 COMPLETED（sacct），原始输出在 `/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts/somatic_truth_in_resources.<set>.tsv`（全部记录，ALT 拆开，没有 BED），还没 `bash collect.sh` 汇进 E8b 的表。百分比是本文件按 n 算的。
  - INDEL，Pansoma 规则 vs ClairS-TO 规则：COLO829T_SMaHT 2,059 个，1,070（52.0 %）vs 1,524（74.0 %）；HCC1395_SEQC2_v1.2.1 1,625 个，119（7.3 %）vs 224（13.8 %）；HCC1954_CASTLE 2,556 个，528（20.7 %）vs 1,027（40.2 %）；H1437_CASTLE 4,896 个，2,508（51.2 %）vs 2,767（56.5 %）；H2009_CASTLE 6,519 个，1,422（21.8 %）vs 1,960（30.1 %）。
  - SNV：COLO829T 44,005 个，1,606（3.6 %）vs 1,378（3.1 %）；HCC1395 39,447 个，1,102（2.8 %）vs 1,127（2.9 %）；HCC1954 24,033 个，1,956（8.1 %）vs 1,699（7.1 %）；H1437 82,529 个，25,323（30.7 %）vs 25,591（31.0 %）；H2009 161,918 个，11,445（7.1 %）vs 12,188（7.5 %）。六套里 SNV 两个规则最多差 1.1 个百分点（HCC1954）。
  - 旁注（不写进本条 Response，和 R1.m1 / R2.2 有关）：H1437 CASTLE truth（DeepSomatic 做的）有 30.7 % 的 SNV 被 Pansoma 的 PoN 规则标掉（ClairS-TO 规则 31.0 %），远高于其他 truth set 的 2.8–10.2 %。
  - 387253_6（HapMap chr20 germline 对照）是 JobHeldUser，和本条无关。
- 跨样本（Response 的 asymmetry 句子）：`experiments/E12a_cross_sample/tables.md`：HG008T d46 SNV 模型 → COLO829T chr1（在 BED 内，PoN 后，模型自己在训练 tumor 上选的 R0.9 阈值）PASS F1 0.708 / 0.754 / 0.630（Illumina / PacBio→fiberseq / ONT）；COLO829T 自己的 d46 模型 0.783 / 0.779 / 0.737。反方向 COLO829T 模型 → HG008T 0.359 / 0.557 / 0.487（HG008T 自己 0.570 / 0.605 / 0.575），主要掉在阈值校准上（best F1 0.554 / 0.597 / 0.579，chr1 曲线上 test-tuned 的点，只作参考）。旧 checkpoint（BQ / MAPQ 用训练 tumor 的统计量），preliminary。旧稿 Response 只写了对 Pansoma 有利的 HG008T → COLO829T 方向，现在两个方向都写；数字是 d46 模型，主比较表是 d9，Response 已注明 HPRC v2.1。
- 评估协议和脚本：`experiments/E8_comparators/README.md` §4：`score_comparator.py` 直接调用 `pansoma_net_v2/vcfeval.py` 的 `write_truth` / `run_pon` / `vcfeval`；rtg 3.12.1，`--squash-ploidy --sample ALT,ALT -f QUAL`；四个 call set（pass、pass_pon、all、all_pon），eval_bed 和 eval_nobed；§6.1 已在两个现成的 DeepSomatic VCF 上端到端测过。
- HG008T ONT 是 dorado 0.8.1 sup v5.0.0，最接近的 ClairS-TO 模型是 `ont_r10_dorado_sup_5khz_ssrs`（sup v4.2.0）；COLO829T ONT 是 dorado 0.4.2 sup v4.2.0（E8 README §3）。
- "applied to the whole truth sets (all records, no benchmark regions)"：E8b §4.5 和 387253 的原始输出都是全部 truth 记录、没有 BED（E8b README §7 Caveats "4.5 counts all truth records (no BED)"），而且 ClairS-TO 规则是本项目按 C6 / C7 重新实现的，不是 ClairS-TO 自己的输出。
- "not named in the DeepSomatic paper"：D7（Nature HTML 全文 "COLO" 0 次；DeepSomatic 是否用过 COLO829 训练 [unverified beyond the text search]）。
- 出入：E8 README 写 "Status PENDING, everything is submitted"，已过时：user 10-13 决定不重跑比较器（E8 取消），只 re-score 已有结果；E8 `summary.json` 的 `datasets` 为空。

**还缺什么。**

- [数据同步] user 给出已有 ClairS-TO / DeepSomatic 结果的路径，以及每个 run 的版本、image、模型名、tumor-only 参数、BAM（aligner 和版本）、参考、命令行或日志。
- [CPU 本地] 用 `score_comparator.py` 对这些 VCF 做四种 call set 的 re-score（eval_bed / eval_nobed），`collect.py` 汇成 Supplementary Table 和 matched-precision 表。
- [用户决定] 原稿 / 新稿报告哪一版、哪种模式：DeepSomatic 是否开了 `--use_default_pon_filtering`；是否混用了 tumor-normal 结果（SMaHT 1.6.0 VCF 像 tumor-normal）；CASTLE 细胞系上 ClairS-TO 用 `ss` 还是 `ssrs`（chr20 / chr1 held-out 冲突）。
- [CPU 本地] `bash E8b_literature_facts/collect.sh`，把已完成的五个 truth set（387253_1-5）汇进 E8b §4.5 的表（387253_6 HapMap 对照 JobHeldUser，可以缺）。
- [GPU RIS] 用 `--norm dataset`（commit f2a0e4e，原 cbc94a7）重训的跨样本 / 留出细胞系模型，替换 Training data 段里 E12a 的 preliminary 数字（R1.7）。
- [文字] Supplementary Table（两个比较器 + Pansoma 同格式）；每张比较表旁边写明同一样本训练的不对称。

### R3.m3 Composition of the pangenome reference

**状态：已解决**: E8b 已经按 HPRC 原始文档和本地 GBZ 的样本表核实了这些事实：两张图的 release / 版本、个体和人群组成、频率过滤方式、benchmark donor 在不在图里；SV 的处理也和 pipeline README 对过。近亲的遗传学检查（E7 s6，v1.1）也已出：没有任何 v1.1 个体带了 benchmark normal 超过 8.2 % 的罕见 germline SNV（一级亲属约一半）。剩下的只是改 Methods 文字；v2.1 的近亲检查是可选的（要下载 VCF，用户决定），因为两张图都按 9 / 90、46 / 464 过滤，单个个体或小家庭特有的 allele 本来就进不了图（Response 已加这一句）。

**Comment（要点）。**

- Document the HPRC release and version, the included haplotypes, the allele-frequency filtering, the treatment of structural variants and the graph construction.
- State whether benchmark individuals or their close relatives contributed to the graph.

**Response.**

The Methods now specify each reference.

- *HPRC v1.1 (release 1).* The Minigraph-Cactus graph with GRCh38 as the reference, in the frequency-filtered form distributed by the HPRC. Nodes and edges carried by fewer than 9 of the 90 haplotypes (about 10%) are removed; the GRCh38 path is never filtered.
  - The 90 haplotypes are 44 individuals × 2, plus GRCh38 and CHM13. The HPRC held out HG002, HG005 and NA19240.
  - The 44 individuals are African 23, American 16, East Asian 4, South Asian 1 and European 0.
- *HPRC v2.1 (release 2).* The same construction, filtered at 46 of 464 haplotypes.
  - The 464 haplotypes are 231 individuals × 2, plus GRCh38 and CHM13. The HPRC held out HG00272.
  - The individuals are African 70, East Asian 50, American 44, South Asian 36, European 29, and two without a population label (HG002, HG005).
  - We used the distributed graph and rebuilt its distance and minimizer indexes with vg 1.77.
- *GRCh38-only graph.* Built by vg autoindex (vg 1.66) from the GRCh38 no-alt analysis set that all truth sets use. Every node carries the same 17 paths, so its haplotype-count channel is constant.

*Structural variants.* The graphs keep the structural variants that pass the frequency filter, and reads align along them. Pansoma's candidates, however, are limited to indels of 1–50 bp; longer events are logged and not called.

*Benchmark individuals.* None of the cancer cell lines used here, nor their matched normals, is in either graph: HG008, COLO829, HCC1395, HCC1954, H1437 and H2009 are absent from both sample lists. These lines are anonymised, so we also tested for close relatives genetically. Of the germline SNVs of the HG008, COLO829 and HCC1395 normals that are rare in HPRC v1.1 (on one or two of its 89 non-reference haplotypes), no v1.1 individual carries more than 8.2%, whereas a first-degree relative would carry about half. In addition, because both graphs keep only nodes on at least 9 of 90 or 46 of 464 haplotypes, alleles private to one individual or a small family cannot enter either graph. [optional 待补: the same test on v2.1]

The HapMap-mixture donors are in the graphs: four of the six in v1.1 and all six in v2.1. R1.5 explains how we handle this.

The HG008 donor is of European ancestry, a group absent from v1.1 and present in v2.1. After population filtering, the results on the two graphs are equivalent (R1.2).

**证据（内部）。**

- `experiments/E8b_literature_facts/README.md` §3.4：
  - H1（v1.1：44 个个体 + GRCh38 + CHM13 = 90；HG002 / HG005 / NA19240 被排除；来源 hpp_pangenome_resources README-v1）。
  - H2（AFR 23, AMR 16, EAS 4, SAS 1, EUR 0）。
  - H3（`--filter 9`，参考 path 不过滤，"fewer than 10% of haplotypes"）。
  - H4（v2.1：232 个组装、231 个进图、HG00272 留出；+ GRCh38 + CHM13 = 464；本地 d46 GBZ 有 462 条样本单倍型）。
  - H5（AFR 70，其中 68 个来自 1000G 人群，外加 NA21309 MKK 和 HG06807 ASL；AMR 44, EAS 50, SAS 36, EUR 29；无人群标签 2 个：HG002, HG005）。
  - H6（`--filter 46`）。
  - H8（HapMap donor：v1.1 有 4 个，v2.1 6 个都有）。
  - H9（HG008, COLO829, HCC1395, HCC1954, H1437, H2009 都不在两张图的样本表里；近亲无法从元数据核实）。
  - E8b §6 对 R3.m3 的状态是 "Resolved (H1-H9)"。
- 近亲检查（E7 s6，新）：`experiments/E7_graph_composition/relatedness.tsv`（2026-10-13 11:57 写出；早上的版本说 "还没产出"，不对），脚本 `s6_relatedness.py`。
  - 方法（脚本 docstring）：每个 benchmark normal（HG008-N dipcall、COLO829BL dipcall、HCC1395BL DeepVariant；FILTER PASS，BED 内，chr1-22）在 HPRC v1.1 里罕见的 germline SNV（89 条非 GRCh38 单倍型里只有 1–2 条携带），数每个 HPRC 个体带了多少。docstring 写 "A first-degree relative shares about half of a person's rare alleles"。
  - 罕见 SNV 数：HG008 95,232，COLO829 92,076，HCC1395 74,604。
  - 每个样本共享最多的个体：HG008–HG00733（AMR）7,674 个，0.0806；COLO829–HG00733（AMR）7,126 个，0.0774；HCC1395–HG01175（AMR）6,112 个，0.0819。都远低于一半。
  - 相对同 superpopulation 中位数的倍数最高是 2.3（HCC1395–NA20129，AFR，share 0.0284），其次 2.2（HG008–NA21309，AFR，0.031）、2.01（COLO829–HG02145，AFR，0.0254）。这些都是 AFR 个体：绝对 share 只有 2.5–3.1 %，倍数高是因为 AFR 的中位数低。没有核实这些个体的混合血统，所以 Response 里只写绝对 share。
  - 只覆盖 v1.1；HCC1954、H1437、H2009 和 HapMap donor 没有做。E7 还没有 README（`s7_tables.py` 没跑）。
  - 8.2 % 是三个样本里的最大 share（HCC1395–HG01175 0.0819）；HG008 最大 0.0806，COLO829 最大 0.0774（`relatedness.tsv` 按 share 排序）。
  - 注意分母只含 HPRC 里有 1–2 条携带的 allele（`s6_relatedness.py`：`found & cnt 1–2`）。如果某个 HPRC 个体是一级亲属，这个分母里来自他的比例会高于一半，所以 "about half" 是保守的说法。
- 频率过滤这一句（Response 新加的）：E8b H3（`--filter 9`，"nodes and edges supported by fewer than 10% of haplotypes" 被删）和 H6（`--filter 46`）。一个个体最多贡献 2 条 haplotype，几个人的小家庭一般也少于 9 / 46 条，所以只在这个个体或家庭里的 allele 进不了 d9 / d46。这也部分回答了 v2.1 近亲问题（v2.1 有 29 个 EUR，HG008 donor 是欧洲人，而 v1.1 的遗传检查里没有 EUR 个体可比）：即使 v2.1 里有近亲，被过滤掉的正是近亲才共享的罕见 allele。GRCh38 path 不过滤，但 GRCh38 不是任何 benchmark donor。
- 来源不一致，已更正：
  - 详版 R3.m3 和简版都写 "The four HapMap donors are in both"。E8b H8 / §5 第 1 条指出 v2.1 是六个全在，所以 Response 改成 "four of the six in v1.1 and all six in v2.1"。
  - 详版写 v2.1 "227 with population metadata … [待核实]"。按 E8b H5 已补全为 231 个个体。
- HG008 donor 是欧洲人：E8b G6（GIAB v0.3 README，"derived from a female of European genetic ancestry"）。d46 ≈ d9（PoN 后）：见 R1.2 的证据。
- vg 版本：
  - d46 用 vg 1.77 重建 dist / min：`/scratch/jshen/Scripts/hprc_v21_d46_vg177_giraffe.sh`，注释写 "the d46 GBZ with the vg 1.77 dist and minimizers"，`VG=/scratch/jshen/bin/vg_v1.77.0`。
  - linear 图：`/scratch/jshen/Scripts/vg_autoindex_hg38_linear.sh`（`VG=vg_v1.66.0`，`REF_FASTA=.../GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`）。
  - H = 17：`pansoma_v2_tensors/HG008T_PacBio_linear/README.txt` 写 `--haplotypes 17`；pipeline README 的说法是 "GRCh38 + 16 path-cover paths on every node"。
- INDEL 1–50 bp：`indexed_gam_pipeline_v4/README.md` 第 240 行（`--max-indel-len 50`，"1-50 and <= --width; longer indels (total over mappings) are logged, not candidates"）。
- 同一个 FASTA 用在所有 truth / labels 上：`experiments/E7_graph_composition/common.py` 的注释（"GRCh38 no-alt analysis set, the fasta of every labels manifest"）。

**还缺什么。**

- [文字] Methods 的参考图一节按 Response 改写。简版和详版里 "four HapMap donors are in both" 要改成 v1.1 4/6、v2.1 6/6。
- [CPU 本地] E7 的 README / `s7_tables.py`：把近亲检查的方法和 `relatedness.tsv` 的结论写进去，Response 的 8.2 % 以它为准。
- [用户决定]（可选）v2.1 的近亲检查要 v2.1 全图 VCF（下载；d46 GBZ 本身不行，罕见 allele 已被 46 / 464 过滤掉）。不做的话，删掉 Response 里的 "[optional 待补: the same test on v2.1]"，写明遗传检查只查了 v1.1，v2.1 靠样本表 + 频率过滤那一句。
- [用户决定] HapMap 的处理（R1.5：直接去掉，还是用已有的 donor-free 图 `mgi_non_HG005.d9`），这一点决定 Response 里 HapMap 那句话怎么写。

### R3.m4 Terminology and language

**状态：已解决**: 只涉及术语统一和语言修改，没有要等的数字；剩下的是改稿本身（NPU 的全称要对照原稿确认，"panel of normals" 的改名要在回复信和正文里一起执行）。

**Comment（要点）。**

- 术语不一致，有若干语法问题。
- "panel of normals" 对 gnomAD、dbSNP、1000G、CoLoRSdb 来说范围太宽。
- 指出的几个具体句子要改。
- NPU 要给出全称。
- Pansoma 和 Pansoma-TO 混用；"Supplement Figure" 和 "Supplementary Figure" 混用。
- 最后做一次语言编辑。

**Response.**

We thank the reviewer and standardised the terminology throughout.

- The method is called Pansoma everywhere (not Pansoma-TO), and figures are cited as "Supplementary Figure".
- NPU is expanded at first use in the main text and in the Supplementary Methods [待核实: the expansion used in the
  manuscript].
- "Panel of normals" is no longer used as the collective name of gnomAD, dbSNP, CoLoRSdb and the GATK 1000 Genomes
  panel (only the last is a panel of normals in the usual sense). We call the four population filtering resources, and
  the step population filtering; each resource and the exact matching rule are listed in Methods (R2.6, R3.m1).
- Terms that recur across the paper (candidate, candidate site, off-reference node, unplaced call, graph-absorbed allele)
  are defined once in Methods and used consistently (R2.m, item 8).
- The quoted sentences were corrected, for example "One key challenge in calling variants on graph alignments is that
  they differ fundamentally …", "Pansoma converts the variants …", "with the technical details described in Methods",
  and "read pileup" in Figure 5.
- The manuscript was edited for grammar and language by [待补: who did the language edit].

**证据（内部）。**

- 修改后的句子来自 `response_to_reviewers_detailed.md` R3.m4 的 Response；原审稿意见和原稿都不在 repo 里，无法逐句核对原文。
- NPU 的全称：详版写 "node pileup (NPU)"，没有文件出处可核，所以 Response 里标了 [待核实]。
- 四个资源是什么：`experiments/E8b_literature_facts/README.md` §3.3 P1–P7（gnomAD 是 GATK af-only 文件、dbSNP 138、GATK 1000G
  PoN 由 Mutect2 在 1000 Genomes exome 上的 calls 构成 [E8b 标为 inference]、CoLoRSdb v1.1.0）。只有 GATK 1000G 那个文件本身
  叫 panel of normals，其余三个是人群数据库，支持改名。旧稿 Response 把四个一起说成 "no longer [called] panel of normals"，
  对 GATK 1000G 那个文件不准，已改成 "no longer used as the collective name ... (only the last is a panel of normals in the
  usual sense)"。
- candidate / candidate site 的定义：`indexed_gam_pipeline_v4/README.md` §8（见 `points/R2.m.md` 第 8 条）。

**还缺什么。**

- [文字] 全文替换：Pansoma-TO → Pansoma；Supplement Figure → Supplementary Figure；"panel of normals" / PoN → population
  filtering（回复信里别的条目，如 R2.6、R1.4、R3.m1，也要一起改，否则回复信自己前后不一致）。
- [文字] 对照原稿确认 NPU 的全称，加在首次出现处。
- [文字] 语言编辑（谁来做、是否用专业编辑服务）。
- [用户决定] 是否接受 "population filtering resources" 这个说法替代 PoN（影响所有条目和 Methods）。
