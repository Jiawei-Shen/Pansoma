# Pansoma 主文稿修改清单

对照日期：2026-11-09

这次修订需要更新论文的研究定位、方法实现、训练与评估设计、主要结果和图表。Response 中已经描述了与投稿稿不同的 pipeline，因此不能只在原稿上补几段解释。下面将三个 reviewer 的 **39 条编号意见合并为 21 个修改主题**，按主文稿的实际位置组织，并保留 reviewer 对照关系。

最优先的五项是：撤回把 graph absorption 漏检解释为 truth-set error 的结论；替换已过时的 NPU、五通道和单节点窗口方法；区分同肿瘤留出染色体与未见肿瘤测试；把 HapMap 移至补充材料并收窄 low-VAF 结论；以 reference ablation 和 candidate composition 支撑 SNV 优势。

## 对照文件和定位方式

| 代号 | 文件 | 用途 |
|---|---|---|
| MS | [Pansoma.pdf](/Users/jiaweishen/Library/CloudStorage/Box-Box/Jiawei-Smaht-Paper/Submit_Folder/Pansoma.pdf) | 主文稿基准，23 个 PDF 页面；下文页码和行号均指此文件 |
| RESP | [Reviewer_Letter_Response.docx](/Users/jiaweishen/Library/CloudStorage/Box-Box/Jiawei-Smaht-Paper/Pansoma_Resubmit/Reviewer_Letter_Response.docx) | Reviewer 意见及作者拟回应；读取当前正文和修订插入内容，不把删除文字作为最终承诺 |
| SUP | [Pansoma_Supplementary_Tables_and_Figures.docx](/Users/jiaweishen/Library/CloudStorage/Box-Box/Jiawei-Smaht-Paper/Pansoma_Resubmit/Pansoma_Supplementary_Tables_and_Figures.docx) | 核对回应对应的表、图、方法说明和未完成事项 |
| WORD | [Pansoma_Ver5.docx](/Users/jiaweishen/Library/CloudStorage/Box-Box/Jiawei-Smaht-Paper/Submit_Folder/Pansoma_Ver5.docx) | 同目录可编辑主稿，辅助核对文字 |

**使用说明：** `R1.1` 表示 Reviewer 1 major comment 1；`R1.m1` 表示 Reviewer 1 minor comment 1，其余类推。表图的 `R...` 编号是目前 response/SUP 的工作编号，提交前需要转换成正式编号。每项复选框表示主文稿尚待落实的工作，不表示实验尚未完成。Response 中的 “we revised / we added” 也不视为这份 MS 已经完成修改。

本清单以文件内容之间的对应关系为依据；没有重跑分析或核验实验代码。来源相互矛盾或尚缺对应材料的内容标为 **待核实**，不能直接作为主稿定稿数字。文档中的审稿要求和作者批注均作为待分析材料处理。

## 按主文稿位置整理的修改任务

### M01 收窄标题、摘要和引言的研究定位

**位置：** p1，L1–2、L11–27；p2–3，L66–104。  
**合并意见：** R1.1、R1.5、R1.6、R3.3、R3.4、R3.8、R3.10。

- [ ] 标题明确 **tumor-only somatic variant calling**。可考虑 `Pansoma enables tumor-only somatic variant calling from pangenome graph alignments`，最终题名由作者确定。
- [ ] 将作为方法定位的 `pangenome-native` 改为 `pangenome-aware`；删除 L89 的 `fully exploit the richness of pangenome data`。技术上确实在 graph coordinates 中生成 candidates 和输出 calls 的事实仍可保留，但不能据此声称完整建模 graph topology、bubbles 或长程 haplotypes。
- [ ] 摘要明确证据主要来自 cancer datasets，主要优势针对 SNV；INDEL 表现依赖平台和模型，不能概括为所有 small variants 一致提升。
- [ ] 摘要加入主要模型的 observed candidate AF ≥0.08 限制，并说明低 VAF 和 INDEL sensitivity 是主要局限。
- [ ] L25–27 关于 graph-specific variants 的结论改为“保留 GRCh38-unprojectable candidates 供后续验证”；不能把不能投射至 GRCh38 等同于已验证的额外 somatic discovery，或不能在任何 linear reference 上表达。
- [ ] 引言末段概括新版的跨节点 tensors、haplotype-support feature、多肿瘤训练与直接 reference ablation；删除 NPU 作为现行核心创新的介绍。
- [ ] L73–79、L89–91 加入关键限定：population graph 中存在某个 allele，并不证明它是该患者的 germline allele；这种 prior 同时可能压低真正 somatic alleles 的候选生成。

**正文应传达：** 在所测肿瘤和平台中，pangenome information 有助于降低 germline candidate burden 并改善 SNV calling；它不能替代 patient-specific matched normal，也尚未验证 normal-tissue mosaicism。

### M02 重写 Overview 和现行 pipeline

**位置：** p3–5，Fig. 1、L110–159；p16，L434–451。  
**合并意见：** R1.1、R1.m2、R3.m1、R3.m4。

- [ ] 将 `GAM → NPU files → tensors` 改为从 **sorted, indexed GAM alignments** 直接构建候选和 tensors 的现行流程；删除 L442–448 的 NPU 文件结构说明和对旧 Supplementary Table 1 的引用。
- [ ] Fig. 1a 同步删除 NPU 文件、索引格式等旧流程元素，保留按节点检索/组织 read evidence 的实际作用。
- [ ] Fig. 1b 和 Overview 明确 tensor rows 随 read walk 跨越相邻节点，不能仍描述为局限在单个节点内的 pileup。
- [ ] L143–148 从简单 somatic/non-somatic 分类描述更新为 somatic、germline、other 三类训练标签；具体的 label eligibility 放到 Methods。
- [ ] Fig. 1c 标明 training、validation、held-out chromosome 与 held-out tumor 的关系；Fig. 1d 标明 GRCh38-evaluable output、off-reference output 和过滤适用范围。
- [ ] 更新 mapper 的实际版本、参数和 graph resource 标识。原稿的 vg 1.66.0 不能未经核对就当作所有修订实验统一版本；response 承诺的是可复现的实际配置。

**依据：** RESP R1.m2、R3.m4；SUP Fig. R1.m2、Table R3.m1。

### M03 替换候选生成规则和五通道 tensor 描述

**位置：** p6–7，Fig. 2a–c、L161–200；p10–11，L297–304；p17，L453–490。  
**合并意见：** R1.1、R1.8、R1.m2、R2.m8、R3.m1。

- [ ] 区分 **node selection** 与 **candidate generation** 两阶段，不能再用同一个 MAPQ 阈值描述它们。
- [ ] Node selection：MAPQ >5 的 alignments；超过 5% 在该 node 上含 mismatch/INDEL 时选中该节点。
- [ ] Candidate generation：MAPQ >10；至少 3 个 ALT-supporting observations；每个 observation 的 allele quality ≥10；AF ≥0.08；INDEL 长度 1–50 bp；每 node 最多 800 reads。确认 allele quality 的实现和名称，不沿用“平均 ALT BQ ≥10”的旧说法。
- [ ] 说明 candidate AF 的计算分母、quality/read cap 的先后顺序以及 full pileup 与最多 200 行 tensor 抽样的关系；这影响低 VAF 解释，不能把 tensor 中的 ALT 比例自动当作完整 pileup AF。
- [ ] 删除 L463–464 的 node-boundary truncation；改为 **8 channels × up to 200 reads × 101 columns**，每行按 read walk 保留跨节点 flanking context。
- [ ] 重新列出全部八通道，不是简单在旧五通道后增加三项。SUP Fig. R1.m2 显示的通道为：read bases、base qualities、site allele carried by each read、mapping qualities、alignment operations、per-read graph-reference bases、node haplotype/path support、read strand。
- [ ] 通道编号统一：图中第 7 个 panel 是零起始编号 `channel 6`；解释原始 path count、相对支持率和网络实际编码/归一化之间的关系。当前图例也有编码说明，不能仅称“haplotype fraction”而不解释。
- [ ] 用按 candidate allele 分组、按 alignment events/node transitions 排序的现行规则替代 L469 的笼统“assigned haplotype path and start position”。
- [ ] 明确定义计数单位：candidate = 一个 `(node, offset, REF, ALT)` allele；candidate site = 同一 node position、同一 variant class 的 passing alleles 分组，并生成一个 tensor。SNV 与 INDEL 即便坐标相同，也计为不同 candidate sites。

**依据：** RESP R1.m2、R2.m8、R3.m1；SUP Table R3.m1、Fig. R1.m2。

### M04 在 Methods 新增 pangenome reference 的完整说明

**位置：** p16，L449–451 附近新增 subsection；同步 p8，L237–248，以及各图 graph 名称。  
**合并意见：** R1.2、R1.5、R3.1、R3.m3。

- [ ] 分别写明 HPRC v1.1 d9：44 individuals 的 88 sample haplotypes，加 GRCh38、CHM13，共 90；HPRC v2.1 d46：231 individuals 的 462 sample haplotypes，加两个 reference paths，共 464。
- [ ] 写明 Minigraph-Cactus release、资源文件、下载标识/校验信息、d9/d46 的过滤对象与规则，以及 GRCh38/CHM13 路径例外保留。
- [ ] 解释两个阈值约为 10% haplotype support。v1.1 与 v2.1 不仅样本数不同，composition 和 construction 也不同，因此不能将结果归因于 haplotype 数量的单一变化。
- [ ] 记录 GRCh38-only graph 的 no-alt analysis reference 和 vg autoindex 构建方法；CHM13-only graph 用于 alignment control，勿与三参考 calling ablation 混同。
- [ ] 说明 SV sequence 可保留用于 alignment，但 Pansoma 的 calling scope 是 SNV 和 1–50 bp INDEL。
- [ ] 说明六个 cancer lines 不在 graph sample list；rare-allele sharing 分析不能正式排除亲缘关系。完整 individuals/populations 清单放补充表。
- [ ] HapMap 单独说明 donor-excluded graph：40 个 HPRC v1 individuals，加 GRCh38/CHM13，共 82 haplotypes，d9；不把它写成标准 v1.1/v2.1 graph。

**正文篇幅建议：** Methods 保留 reference 定义和关键排除规则；完整 manifest 放 Supplement。Personalized pangenome 只作为未测试的未来方向。

### M05 重写 truth projection、标签生成及 INDEL 等价表达规则

**位置：** p4，L143–148；p18，L503–508、L517–534；Methods 新增独立 labeling subsection。  
**合并意见：** R1.4、R1.8、R2.2、R2.3、R3.m1。

- [ ] 区分 truth loci 在 GRCh38 embedded path 上的定位、独立发现 read-derived candidates、以及推理后 call projection 三个步骤。
- [ ] 删除 L529–532 从 adjacent 1-bp SNV nodes 借用坐标来推断 graph-only positions 的旧描述。写明 unique path placement、歧义处理和 left-alignment/normalization。
- [ ] 明确 site 的代表 allele 由最多 ALT reads 决定；仅 site 内另一个 allele 与 truth 相符，不自动将该 tensor 标为 somatic。
- [ ] 正文说明 somatic=1、germline=2、other=0 与排除/未标注情况；Supplement 列出规则优先级、somatic ∩ germline confident regions、off-reference candidates 和 final INDEL exceptions。
- [ ] 原有 >45% overlap 的 partial somatic annotations 不等于当前 positive training examples；当前 training 排除 partial somatic matches。不要把 residual INDEL 当成完整 GRCh38 truth allele 的 TP。
- [ ] 说明 final INDEL recipe 将 somatic benchmark 内、germline coverage 外的部分 candidates 当作 other 的假设；这些不是已确认的 sequencing artifacts。
- [ ] 说明未标注 off-path candidates 在 validation 中作为 negatives 的 model-selection 假设；不能把这一假设当作它们确实 non-somatic 的生物学证据。
- [ ] 交代训练用 exact/partial label 规则与评估时 rtg vcfeval haplotype-aware matching 的差别。训练 labels 不能直接用来决定 chr1 TP/FP/FN。
- [ ] 明确 ALT 被 graph path 完整吸收时无 candidate、无 tensor；这种 truth allele 仍保留在 benchmark 分母，漏检计 FN。

**依据：** RESP R1.8；SUP Note R1.8 是更详细的实现说明。相关 validation-negative sensitivity 检查尚有 pending 项，见后文 D08。

### M06 更新网络、训练配置和决策阈值

**位置：** p17–18，L491–516；Fig. 1c；旧 Supplementary Fig. 6 的正文引用。  
**合并意见：** R1.6、R1.7、R1.8、R3.9、R3.m1。

- [ ] 用新版 SNV 50.4 M、INDEL 8.7 M 参数网络替代原 199.6 M 级网络；各自的 depths、widths、额外模块及实际 parameter count 由最终配置表确定，不能保留旧 `[3,3,27,3] / [192,384,768,1536]` 就只改参数数字。
- [ ] 更新训练：SNV 20 epochs、INDEL 12 epochs；AdamW，learning rate 10⁻⁴，weight decay 0.05，1 warm-up epoch 后 cosine decay，batch size 1,024/GPU。删除旧 80 epochs、3 warm-up epochs、weight decay 0.01 及固定 positive weights 50/200 的描述。
- [ ] 交代按 inverse class frequency 的 0.75 次幂（SNV）和 0.5 次幂（INDEL）生成 class weights；GPU count、gradient accumulation、normalization、seeds 和 checkpoint 列入逐模型配置表。
- [ ] SNV 的 checkpoint 和 PASS threshold 分开解释：best validation somatic F1 checkpoint；validation recall 0.9 operating threshold。这不是 test recall 的保证，也不是 candidate AF threshold。
- [ ] **待核实：** RESP R3.m1 将 PASS 概括为 recall 0.9，并提 INDEL best-F1；SUP Table R3.m1/Note R1.8 的 final INDEL 是 **post-filter validation recall 0.8**，checkpoint 选 post-filter validation F1。先统一协议再写入主稿。
- [ ] 增补 final INDEL 的五个 site features 和与 ConvNeXt 的融合方式；原 Methods 和 response 的简表都不足以说明这一变化。SUP 记录：exact-allele population AF、同位点 INDEL 最大 AF、filter tag、homopolymer length、是否有 GRCh38 allele；经 5→64 linear layer/GELU 与 pooled features 合并。
- [ ] 因为这些 INDEL features 使用 GRCh38 allele 和 population priors，原稿“完全不依赖 linear projection”的表述应限定到实际适用步骤，不能概括整个最终 INDEL 模型。
- [ ] 报告各 class 的 training counts、模型容量对照和 learning-curve 结果；正文一句说明选择较小模型的依据，详细 counts/curves 放 Supplement。
- [ ] 区分 main models 无 on-the-fly tensor augmentation，与专门 dilution models 使用 in-silico mixtures。

**依据：** RESP R1.7、R3.m1；SUP Tables R1.7a、R3.m1。Final INDEL 方法变化不能仅藏在补充表脚注里。

### M07 重写 benchmark 数据集和 training exposure

**位置：** p3，L97–104；p7–8，L201–236；p18，L503–508；p19，L550–565。  
**合并意见：** R1.7、R1.m1、R2.2、R3.8、R3.m2。

- [ ] 主 benchmark 改为 HG008T、COLO829T、HCC1395；新增 HCC1954、H1437、H2009 的训练用途和数据来源。不要把六个 training tumors 全部称为独立验证集。
- [ ] 增加 truth-set overview：HG008-T draft v0.2 的 assembly-based tumor/normal comparison；SMaHT COLO829 的 multi-caller + PacBio read support；SEQC2 HCC1395 v1.2.1 high/medium-confidence；CASTLE truth 的 DeepSomatic 来源。版本、variant class、PASS/region 范围和样本量分母须配套。
- [ ] 说明 matched-normal germline truth 的来源：HG008N/COLO829BL 的 assembly/dipcall，SEQC2/CASTLE 的 DeepVariant calls；matched normals 用于 labels、validation evidence 和 benchmark construction，与 tumor-only inference 的使用范围分开。
- [ ] 明确 single-tumor、four-tumor 和 six-tumor 三类模型的作用。Four-tumor models 排除 HG008T 和 COLO829T；six-tumor models 加入这两者的训练染色体。
- [ ] 删除 L506 的统一表述 `HG008T was excluded from model fitting`，并重写 L274–276、L563–565 的 independent/generalizability 结论：它仅适用于确实未包含 HG008T 的模型；six-tumor 的 HG008T chr1 是 within-tumor held-out-chromosome test。
- [ ] 给出按 tumor × model family 的训练/验证/测试染色体表。HCC1395 排除 chr1 和 chr20；**final INDEL 与 SNV 对 CASTLE 染色体使用范围不同，见 D02**。
- [ ] Validation 使用 20,000-node blocks 中的 5%，而非简单逐 candidate random split；注明 seeded hash 以及 graph/seed 间 validation blocks 不完全相同。
- [ ] 解释 CASTLE 不作为主要独立 benchmark 的原因：truth construction 与 comparator 相关；部分 alleles 与 matched-normal LOH 情形相关。保留其训练用途，但避免泛化为其他 truth sets 普遍不可靠。

**依据：** RESP R1.7、R1.m1、R2.2；SUP Tables R1.7a–c、R2.2、R3.m1、R3.m2。

### M08 补全 comparator 配置和公平评估口径

**位置：** p19–20，L550–590，尤其 L575–590；所有性能图图注。  
**合并意见：** R1.6、R2.3、R3.4、R3.9、R3.m2。

- [ ] 写明评估范围：chr1、对应 benchmark regions、GRCh38-projectable calls；不是 complete graph output 的 precision/recall/F1。另有 chr2 held-out-tumor INDEL confirmation 时单独说明。
- [ ] 将“standard workflows”展开为软件版本、aligner/version、checkpoint、tumor-only flags、resource versions 和 PASS 定义。主要版本：DeepSomatic 1.9.0；ClairS-TO 0.4.2（HG008T/HCC1395）与 0.4.1（COLO829T）；详细 per-platform 表用 SUP Table R3.m2。
- [ ] 明确 HG008T ONT DeepSomatic 已改用 `ONT_TUMOR_ONLY`；旧 `PACBIO_TUMOR_ONLY` 的结果不能继续用于 Fig. 3、SBS96、Venn、unfiltered comparisons 或正文结论。
- [ ] 说明 HG008T ONT ClairS-TO SUP v4.2 model 与数据的适配限制，以及 COLO829T Illumina 使用 `ilmn` 而不是其他样本的 `ilmn_ssrs`。
- [ ] 说明训练暴露差异：HCC1395 对 Pansoma multi-tumor 和 DeepSomatic 是 chr1 留出，对 ClairS-TO SSRS 是未见 tumor；不能称它对所有方法都 fully independent。
- [ ] 所有 callers 使用相同 truth/regions、bcftools normalization、rtg vcfeval 3.12.1 和 `--squash-ploidy` 的 allele/haplotype matching；基因型不作为评分条件。
- [ ] 分别定义 native PASS、native PASS 再加 common population filter，以及 PR curve 的 all variant records（排除 RefCall/NoCall，再 common filter，按 QUAL 排序）。这些不是同一调用集合。
- [ ] 说明 harmonized PASS 是在各工具原生过滤之上追加共同过滤，确保共同的最低过滤条件，**并不代表所有工具经历完全相同的总过滤**。
- [ ] 若报告 matched-precision recall，说明它是在 test PR curve 上取得的描述性最大 recall，不能作为预先在 validation 选定的 operating point。

**依据：** RESP R3.m2；SUP Table R3.m2 的脚注比 response 更详细，应同步进入 Methods/可复现 Supplementary Note。

### M09 重做主性能结果并分别评价 SNV 和 INDEL

**位置：** p9–11，Fig. 3、L250–320；p20，L592–598。  
**合并意见：** R1.6、R1.7、R1.m1、R2.3、R3.9、R3.10。

- [ ] Fig. 3 和对应文字用最终模型、正确 comparator checkpoint、统一 filtering/evaluation 后的结果重做，不能将旧图作为主证据而只在正文补新版数字。
- [ ] 主图至少同时展示 precision、recall、F1；正文明确更高 F1 不等于每种使用情形下 sensitivity 足够。完整 PR curves 和 candidate→classification→filter 的 recall funnel 应可追溯。
- [ ] 主要 release models 为 six-tumor SNV/INDEL models；cross-tumor evidence 独立显示 four-tumor models。不要跨平台挑选最优 Pansoma model 后把结果称为同一个 release model 的性能。
- [ ] 加入 HCC1395 benchmark；对 HG008T/COLO829T six-tumor results、four-tumor results、within-tumor single models 分开作图/标注。
- [ ] 原稿“HG008T WGS/PacBio INDEL 均最好”和“ONT 始终较差”的总结需按 final models 重写。例如 SUP Table R1.7c 中，six-tumor HG008T Illumina INDEL F1 为 0.199，DeepSomatic 为 0.208；HG008T ONT 则为 0.168，高于两 comparator。不能保留旧版平台排名。
- [ ] 可报告 response 中 six-tumor ONT INDEL 相对相应 single-tumor model 的提升：HG008T 0.030→0.168、COLO829T 0.006→0.090；同时注明 COLO829T ONT 仍低于 ClairS-TO 的 0.103。这些是特定模型比较，不是“换 graph 就改善 INDEL”的证据。
- [ ] INDEL 按 insertion/deletion、event length、homopolymer/tandem repeat、platform 分层；正文用一段概括，TP/FP/FN 和 truth counts 放 SUP Table R3.10。
- [ ] 不将某些 INDEL classifier/operating-threshold 的局限简化为一个普适原因。需要区分 candidate availability、population filtering、homopolymer errors 和具体 operating threshold 对不同模型的影响。

**依据：** SUP Tables R1.6a、R1.7b/c、R3.9a–c、R3.10；RESP R1.6 的 recall ranges 是其特定分析的结果，不能直接套到 six-tumor release model。

### M10 新增直接证明 pangenome 贡献的 ablation 结果

**位置：** 建议在原 p10，L278–308 的机制解释前新增 Results subsection；Methods 补实验设计。  
**合并意见：** R1.1、R1.2、R1.m2、R3.1、R3.9。

- [ ] 主文展示 **GRCh38-only graph / HPRC v1.1 d9 / HPRC v2.1 d46** 的三参考对照。相同 reads、candidate/tensor pipeline、训练方案与过滤，分别重新训练；说明 haplotype-count normalization 随 reference 调整，GRCh38-only channel 为常数。
- [ ] 清楚报告三种 reference 的 SNV P/R/F1，而不是只说 Pansoma 胜过其他 callers。RESP 报告在六个 tumor–platform datasets × 两个 HPRC graphs 的 12 个配对比较中，mean post-filter SNV F1 增加 0.028–0.292；最终稿需与 Table R3.9a 的模型和 seeds 对齐。
- [ ] 同时报告 pre-filter 与 post-filter 的 SNV reference benefit，建立这种改善不仅来自后置 population filter 的证据。
- [ ] 单独报告 fixed-alignments/fixed-candidates 的 channel-6 retraining ablation：去除后 pre-filter F1 在 12/12 比较下降，post-filter 在 11/12 下降；明确例外为 HG008T ONT v2.1。
- [ ] 区分 cross-node context 的两种测试：inference-only masking 会降低 recall；mask 后重新训练的 v2.1 模型 PASS F1 与 controls 差值均在 0.016 内。因此可以证明模型使用了 context，但不能宣称 cross-node context 带来普遍且很大的独立 F1 增益。
- [ ] 明确 INDEL 在三参考对照中没有一致优势；不能把 SNV 的 12/12 改善外推到 INDEL。
- [ ] Methods 说明 seeds、normalization、GPU/global batch、validation block selection 的实际差异。SUP 各表不是完全同一套模型，跨表拼接数据时须核对。

**主文必须保留：** 三参考 SNV 对照结论和 channel-6 ablation 的关键结论。完整逐 seed、逐平台值放 Supplement，但不能让主文只引用回复信来支撑核心 novelty。

### M11 用候选组成和 false-positive 分析解释 precision 增益

**位置：** p8–11，L237–248、L278–308；p13–14，L368–408。  
**合并意见：** R1.2、R2.5、R2.6、R3.2。

- [ ] 用 somatic truth、matched-normal germline、background/non-truth 三类 retention 替代旧的“总 candidates 变少，因此 precision 更好”的单一步推断。
- [ ] 更新 L298–304 的旧阈值和旧 candidate 总数比例。SUP Table R1.2 的 germline candidate reduction 为 v1.1 84–90%、v2.1 85–91%；注明该表计数范围与 representative-allele 定义。
- [ ] 在 benchmark regions 内，Table R3.2a 报告 v2.1 保留约 89–96% somatic SNV candidates，germline SNV candidate fraction 从 GRCh38-only 的 82–88% 降至 8–9%；这些分母不同于 Table R1.2，不能混用。
- [ ] 增补 graph-represented vs graph-absent germline FP 的 selective reduction。Table R3.2b 中，GRCh38-only pre-filter SNV FP 有 71–86% 对应 graph-represented germline alleles；转到 v2.1 后这批 alleles 仅 0.2–0.6% 保留 tensors。
- [ ] 增补 candidate 周围 ±50 bp 内 germline density 分层，排除 candidate 自身；HG008T 在 ≥2 nearby germline variants 区域的 SNV F1 增益最大，response 给出 +0.107–0.204。
- [ ] 区分 HG008T 的 germline-associated FP reduction 与 COLO829T extreme-region FP reduction；这些类别重叠，不能相加为独立解释比例。
- [ ] 补充 graph-absent rare germline alleles 的残留：模型和过滤都不能完全剔除，HG008T post-filter FP 仍以这类为主。说明 graph 的 population prior 并未解决 patient-specific germline identification。

**依据：** SUP Tables R1.2、R3.2a–c；Fig. R2.5；RESP R2.5/R2.6。原 Fig. 2e 的 graph overlap 可作为背景，但不能独自证明 calling precision 增益。

### M12 撤回 truth-error 推断并重写 false-negative 和 graph absorption 分析

**位置：** p10，L279–288；原 Supplementary Figs. 2–3 的正文引用；p20–21 Discussion。  
**合并意见：** R1.2、R1.4、R1.6、R1.8、R2.2、R2.3、R3.2、R3.7。

- [ ] 删除/重写 L281–283 从“truth variant overlaps graph node”推断它可能是 germline/错误 truth 的论证。ALT 已被 population graph 编码与患者是否 germline 是两个问题。
- [ ] 明确 complete absorption 和 partial absorption：前者 ALT reads 无 edit，后者留下不同于完整 truth allele 的 residual edit；两者都可能导致 genuine somatic alleles 无 exact candidate。
- [ ] 原 `variant node not ranked #1`、`first aligned node mismatch` 分类和旧统计撤回，改为基于 read evidence 与首个 loss stage 的分类；不要新旧分类并列而不给对应关系。
- [ ] 加入完整 truth denominator 下的分阶段漏检：candidate generation、classifier threshold、population filtering；分别展示 SNV 和 INDEL。
- [ ] 同时给出**未过滤时的候选损失**和**共同 population filter 后剩余的额外损失**。例如 Table R1.4：v2.1 相对 GRCh38-only 的 INDEL candidate loss 为 HG008T 29.4–42.1%、COLO829T 23.9–32.0%；after-filter graph-specific loss 分别为 3.5–7.0%、5.4–6.9%。不能只写较小的后一个数字，使读者误以为全部 graph absorption 很小。
- [ ] 说明 graph-specific loss after filtering 0.2–0.7% SNV、3.5–7.0% INDEL 是所列 chr1–22 truth sets 的范围；不是总体 FN rate、chr1 benchmark recall 或 conditional percentage。严格 exact-path absorption subset 又是另一种定义。
- [ ] 将 matched-normal reads、phased normal assemblies、tumor assembly 的证据引入正文，支持“相当多被吸收 alleles 是 genuine somatic events”的解释；**具体 counts 先解决 D03 的 response/SUP 冲突**。
- [ ] 对 normal heterozygous allele 与 haplotype-specific somatic event 的重合，解释 phase ambiguity，不能仅凭 allele overlap 宣判 truth 错误。
- [ ] COLO829T tumor assembly 由正在分析的长读段构建，应称 consistency evidence；不能称完全独立的 orthogonal validation。Assembly 没找回 truth allele 也不能反推它必定是 germline。
- [ ] 这些 alleles 仍保留在 truth set，并在 evaluation regions 内按漏检计 FN。原 Supplementary Fig. 3 标题改为 `Benchmark somatic truth alleles missed through graph absorption`，正文同步。

**主文重点：** graph suppresses population-carried germline variation，也会 suppress population-carried true somatic variation；这是可测量的 sensitivity–precision trade-off，而非“graph 纠正了错误 benchmark”。

### M13 调整 low-VAF、tumor purity 和 HapMap 的位置及结论

**位置：** p7–8，L201–228；p9–11，Fig. 3c、L265–276、L309–320；p19–20，L551–575；p21，L615–618。  
**合并意见：** R1.5、R1.6、R2.1、R2.4、R3.8、R3.10。

- [ ] 将 Fig. 3c 和 HapMap 主结果 L309–320 移到 Supplement；主文保留简短定位和引用。删去其证明 tumor-only somatic accuracy、正常组织 mosaicism 或普遍低 VAF 优势的结论。
- [ ] HapMap 定义为约 500×、truth VAF 1–16.5% 的 mixture allele detection；positives 来自 minor donors 的 germline alleles。
- [ ] 在补充方法中交代所有六个 donors 从专用 graph 排除，仍存在其他 graph individuals 共享的 alleles。共同 evaluation mask 排除 graph-carried alleles 及 reference span ±1 bp，两 caller 使用同一 mask。
- [ ] 明确 HapMap 不应用 population filtering；因此它测试的是低 fraction allele detection，不能测试完整 somatic–germline discrimination。不能沿用“同样 four resources 过滤所有数据”的泛化描述。
- [ ] 更新 low-VAF rerun：candidate floor 0.015、最多 200-read tensor；删除旧最大 500-read tensor 作为现行方法的说法，并核实其他 AF 上限规则。
- [ ] **原稿已包含 COLO829 tumor/normal dilution，所需动作是重做口径和解释，而非仅新增一个原稿没有的实验。** 用 revised joint models、actual purity、quality normalization 与 expected-VAF stratification 替换原“25–100% purity 均稳健”的笼统结论。
- [ ] 说明主模型 AF ≥0.08、≥3 ALT reads 的 detection floor；expected VAF 不等于 observed graph-pileup AF，低于 expected 0.08 的 bins 出现少量 recall 不与 candidate rule 矛盾。
- [ ] RESP 给出的低于 expected VAF 0.08 时 recall ≤0.40、≥0.20 时 0.78–0.93，应限定到 SUP Table R2.4b 的 **Fiber-seq 和 ONT joint dilution models**；该表未提供 Illumina 同类 VAF 分层。
- [ ] 加入 Supplement Note R2.4 的关键控制：Illumina/ONT 做 per-read base-quality centring，避免模型用 tumor/normal read quality 差异识别混合来源；说明 nominal purity 与 measured purity，并披露后者受 candidate floor 的选择偏倚。
- [ ] 说明 joint dilution models 在同一肿瘤 chr2–22 的多个 purity 上训练；这是 purity sensitivity，而非独立 tumor generalization。Illumina comparators 与 Pansoma 使用同源 reads 的不同随机子集，长读段 mixtures 相同。
- [ ] 把 L617–618 “specialized low-VAF model partly addresses this limitation”改为：HapMap model 仅是 exploratory mixture option，不能作为标准 tumor low-VAF calling 的已验证解决方案。

**依据：** SUP Tables R1.5、R2.4a/b、Note R2.4；RESP R1.5、R2.1、R2.4。

### M14 重写 alignment 改善并新增 difficult-region 结果

**位置：** p10，L289–296；建议新增 Results 小节；Methods 新增 stratification/simulation 方法。  
**合并意见：** R1.2、R1.3、R2.5、R2.7、R3.6。

- [ ] 删除把 perfect alignment 增加直接解释为 placement 更准确的表述；区分 reference agreement、read placement 和 downstream calling accuracy。
- [ ] 用 absolute counts + proportions 更新“up to 11% / more than four-fold”。例如 HG008T Illumina perfect-read proportion 为 68.3%→78.5–78.7%，约 10 percentage points；长读段原始比例非常低，大倍数增长不能替代绝对比例。
- [ ] 分开报告 linear GRCh38、GRCh38-only graph、HPRC v1.1、v2.1；指出 soft-clipping 降低主要发生在前两者之间，涉及 mapper/workflow 变化，不能全部归因于 HPRC allelic sequence。
- [ ] 主文补充 SMaHT easy/difficult/extreme、GIAB v3.5 low-mappability/segmental duplication、tandem repeats/homopolymers 等区域的 calling analysis；Methods 给 BED/version/overlap handling，避免只有类别名称没有定义。
- [ ] 报告 HPRC 在 GIAB difficult union 中的 SNV F1 改善，RESP 给出的配对范围为 +0.054–0.483；与 1-Mb block-bootstrap CI 配套，而非只给全基因组均值。
- [ ] 同时报告 ≥7 bp homopolymers 中 long-read INDEL precision 的问题；注明 Table R1.3a 的 INDEL 是 earlier three-reference models，不能冒充 final release INDEL 分层结果。
- [ ] 补充以 phased assemblies 为 truth 的 read simulation，含 linear GRCh38、GRCh38 graph、CHM13 graph、两个 HPRC graphs；限定 chr1/chr20 和在两个 linear references 上都有 one-to-one truth 的可评估区域。
- [ ] 明确 Illumina/HiFi-like reads placement 改善，但 ONT-like reads 在 v2.1 misplacement 约翻倍、集中于 extreme regions。Real-read high-MAPQ loss 也应呈现，不能只保留正向结果。
- [ ] Expanded repeat diversity 导致 noisy-read ambiguity 仅作为解释性假说；没有单独完成的 read-level haplotype-concordance 实验不能写成已完成。

**依据：** SUP Tables R1.3a/b、R2.7、R3.6/R3.6a；Fig. R2.7。

### M15 重新界定 graph-coordinate output 和验证结果

**位置：** p12–16，Fig. 4b/c/e、Fig. 5、L368–431；p18，L517–534；p20–21，L599–607、L615–616。  
**合并意见：** R1.3、R1.8、R3.3、R3.4、R3.8。

- [ ] 用 `GRCh38-unprojectable` 描述投射状态；不要自动称 `graph-exclusive` 或 `graph-specific somatic variants`。大量这类 calls 可在 CHM13 上表示。
- [ ] 删除 L603 的 `about 2%`；用按 sample/platform/reference/output mode 的比例替代。RESP 的 default output 给出 HG008T/COLO829T 0.5–10.1%，HCC1395T 0.4–2.2%；具体表格支撑尚需核实，见 D04。
- [ ] 明确 optional off-reference rescue 默认关闭；default 与 rescue-enabled 输出的数量、验证率、SBS96、CHM13 projection 比例不能混用。
- [ ] 正文报告 matched-normal assembly/read、tumor assembly 和 cross-platform evidence 的主要结论：大多数 off-reference calls 不能被当作已经验证的新增 somatic variants。
- [ ] RESP 报告 default output 的 tumor-only assembly support 为 v2.1 34/541、v1.1 17/287；这是潜在 somatic 支持，不是已建立完整 truth 下的 precision。**先补齐与当前 SUP 对应的 default-mode 明细再使用这些数字。**
- [ ] 删除/改写 L386–398、L427–431 关于保留可靠 tumor-associated signal、超越 linear workflow 的过强结论；SBS96 concordance 不能代替 allele-level validation。
- [ ] 在 Methods 和 benchmark 图注中声明：主 P/R/F1 仅衡量 GRCh38-evaluable calls；图外 calls 不能因缺少 truth 就默认 TP 或从 full-output accuracy 解释中忽略。
- [ ] 若保留 SUP Table R3.4 中“所有 off-reference calls 当 FP”的分析，标为 conservative sensitivity analysis；它不等于真实 full-output accuracy。还要区分 unrestricted chr1 和 benchmark-restricted 结果。
- [ ] Fig. 5 若保留，只能作为 representation 示例，或换成有完整 tumor/normal/CHM13 evidence 的候选。存在于 HG008N alternative sequence 的**背景序列**不等于 ALT allele 已在正常中，也不证明该 ALT 是 somatic；图注必须把 sequence context 与 allele evidence 分清。
- [ ] Discussion 交代 reference-path positive training distribution，以及未标注 off-path validation negatives 可能影响 off-reference behavior；这是假说与模型局限，不是已证明的机制。

**依据：** RESP R3.3/R3.4；SUP Tables R3.3、R3.4。两个来源的 output mode 与表编号尚未完全一致。

### M16 用直接 SBS96 比较替换当前 signature 展示

**位置：** p12–14，Fig. 4a/c、L324–367、L388–398；Fig. 2f 的引用。  
**合并意见：** R2.m1、R2.m5、R2.m6、R2.m7、R3.5。

- [ ] Fig. 4a 按每个 tumor/platform 并排显示 chr1 truth、Pansoma、DeepSomatic-TO、ClairS-TO 的 SBS96，并加入 matched-normal germline spectrum。
- [ ] 将 COSMIC decomposition 的比较文字和 cosine 定义同步改为直接 SBS96 profiles；不能换图后仍用旧 signature-weight vectors 的 cosine 值。
- [ ] 更新 call sets、Venn counts 和 cosine，尤其 HG008T ONT comparator 必须用修正后的模型。
- [ ] L349 的 `truth-overlapping calls` 改为 exact/haplotype-equivalent truth-matching true positives；`shared core` 单独定义为三个 callers 都检出的 TP 集合，不是只在同一位置重叠。
- [ ] Graph-coordinate SBS96 放 Supplement，并交代 graph context extraction、缺失 context 排除数、小样本与 ONT subsampling、不能套用 GRCh38 population filter/benchmark masks 的限制。
- [ ] 标明当前 SUP Fig. R3.5b 使用 rescue-enabled calls；不能把其 6,025-call 总量或谱形当成 default output 的证据。
- [ ] 删除基于 spectrum 相似就判断单个 off-reference call 真实性的推断；正文用其描述 spectrum，不作为 somatic authenticity 的决定性证据。

**依据：** SUP Figs. R3.5a/b、Table R3.5；RESP R2.m6、R3.5。

### M17 增加模型稳定性、区间估计和 operating-point 说明

**位置：** p9–11 性能结果；p18 training Methods；p20 evaluation Methods。  
**合并意见：** R1.6、R1.7、R3.9。

- [ ] 报告三 seed 的 mean ± sample SD，注明哪些实验是三次训练、哪些 multi-tumor SNV/cross-tumor 结果只有一次训练；不能统一声称全部结果均经过 multi-seed validation。
- [ ] 增补 paired 1-Mb block bootstrap、2,000 resamples、95% intervals 的方法及适用对照。只做 HPRC-vs-GRCh38 interval 的分析不能写成已覆盖所有 caller comparisons。
- [ ] PR curves 和 matched-precision recall 补充 F1 point estimates；阈值从 validation 选择，test curves 仅用于事后性能描述。
- [ ] 明确 COLO829T Illumina v1.1 的较早 checkpoints 是查看 chr1 后的 exploratory analysis；正式结果只能用预先定义的 validation selection。SUP Table R3.9a 明确写了这些 early-epoch values 不用于任何 comparison。
- [ ] 不把三 seed dispersion、same-seed rerun reproducibility、跨 cluster inference reproduction 混为一种稳定性证据。SUP 记录部分 final INDEL same-epoch rerun F1 差可达 0.041，不能写完全 bit-reproducible。
- [ ] 将 held-out-tumor chr2 INDEL confirmation 独立报告；SUP Table R1.7c 的平均差为对 DeepSomatic +0.034、对 ClairS-TO −0.021，不能据此声称两个 comparator 均被超过。

**依据：** SUP Tables R3.9a–c、R1.7c、R1.8。RESP 对 R3.9b 的用途有误引用，见 D07。

### M18 重写 population filtering 方法和 unfiltered 结果

**位置：** p4–5，L157–159；p12–14，Fig. 4d、L399–408；p19，L535–549。  
**合并意见：** R1.2、R1.4、R2.6、R3.1、R3.m1、R3.m2、R3.m4。

- [ ] Subsection 改为 `Population filtering`；四类资源统称 population filtering resources。`Panel of normals` 留给真正的 GATK 1000 Genomes PoN 或 comparator 对应功能。
- [ ] 更新 filtering rules：SNV 使用 gnomAD/CoLoRSdb AF ≥10⁻⁴、non-somatic dbSNP 138、GATK 1000 Genomes PoN；INDEL 仅 gnomAD/CoLoRSdb，homopolymer ≥8 bp 时 AF ≥0.01，其余 AF ≥10⁻⁴；按 exact allele match。
- [ ] 删除旧 CoLoRSdb AF ≥0.001 作为 Pansoma 统一阈值的说法；ClairS-TO 的原生 0.001 设置是另一流程，不应一起替换。
- [ ] Fig. 4d 扩展到 HG008T，并按最终分析需要纳入 HCC1395；SNV 与 INDEL 分开报告，明确 reference/model。
- [ ] 将 `PoN-independent` 改为更准确的 `before post-calling population filtering`。DeepSomatic 仍有 population AF input；final Pansoma INDEL 也有 population AF features，不能说它们不使用 population information。
- [ ] 不把 SUP Table R2.6a 与 Table R3.9a 的任意数值直接相减作为纯 filtering effect；两表 HPRC runs 的 normalization 设置并不完全相同，需同 checkpoint/同 prediction 配对比较。
- [ ] 说明 HapMap 未做 population filter；对 gnomAD/dbSNP 样本重叠只能说未发现证据但不能完全排除，因没有完整公开 sample lists。详细 provenance 和 rare-allele overlap 放 Supplement。
- [ ] 主文解释为何 graph 和后置过滤能相互补充，又为何 population tagging 不等于患者 germline status。

**依据：** SUP Tables R2.6a/b、R3.m1、R3.m2；RESP R2.6。

### M19 按修订证据重写 Discussion

**位置：** p20–21，L591–623。  
**合并意见：** R1.1–R1.8、R2.1–R2.7、R3.1–R3.10 中的解释与范围要求。

- [ ] 第一段将优势限定为有证据支持的 SNV settings，并依次连接三参考 ablation、germline suppression 和 graph-derived feature 的证据；不再仅凭与其他工具比较将提升全部归于 graph-native design。
- [ ] 第二段用 cross-node tensors 和 haplotype-support features 替换 NPU 核心创新；说明单一 path-count feature 和局部 read walk 不等于全面建模图拓扑。
- [ ] 增补 graph absorption 的真实 sensitivity cost、INDEL partial representation 与 population-filter overlap；保留原 truth alleles，避免再推断 benchmark errors。
- [ ] 增补 tumor-specific learning 的实证边界：COLO829-only transfer bias、four-tumor 对 unseen COLO829 的不足、six-tumor within-tumor 性质；不能宣称已证明所有 tumor types/mutational processes 的泛化。
- [ ] 改写 L608–613 “larger graph 将带来更多准确性和泛化收益”的预测：增加 haplotypes 同时可能增加 noisy-read ambiguity 和 absorption；个体化 graph 尚未在 Pansoma 测试。
- [ ] 删除 `about 2%` 和 low-VAF model 已解决主模型低频缺陷的说法；写清 off-reference exploratory status、CHM13 representability 和 AF 0.08 floor。
- [ ] 将 normal-tissue mosaicism 明确列为未验证的未来方向；保留 tumor-only cancer scope，避免把研究方法表现写成已验证临床效能。

**写作建议：** Discussion 可组织为“证据支持的贡献 → precision/sensitivity 权衡 → 泛化与 graph-only 边界 → 未来工作”。不需要逐 reviewer 再复述一次全部实验。

### M20 更新 Data availability、Code availability 和 References

**位置：** p21–23，L624–713；p23 Acknowledgements 视新增资源需要更新。  
**合并意见：** R1.m1、R2.2、R3.m2、R3.m3，以及已承诺的可复现性要求。

- [ ] Data availability 增加 SEQC2 HCC1395、CASTLE HCC1954/H1437/H2009、相关 matched-normal/assembly/truth resources；区分训练和测试用途。
- [ ] 增加两个标准 HPRC releases、HapMap donor-excluded graph、GRCh38/CHM13 control references、stratification BEDs、population resource versions 和文件 manifest。
- [ ] Code availability 删除“node pileup generation”指 NPU 文件生成的过时表述，更新为实际 sorted/indexed GAM tensor pipeline、projection、filtering、evaluation scripts。
- [ ] 指向对应 revised release 的 SNV/INDEL six-tumor checkpoints、必要的 four-tumor models、配置/normalization/thresholds、可复现命令和归档版本；不默认原 Zenodo DOI 已包含新版内容。
- [ ] 补充新增 truth-set、HPRC release、mapper、population-filter trade-off 等对应引用，并核对 response 和 References 中作者/年份/版本的一致性。此项是提交前 bibliographic 核验任务，不以回复信中的暂写引文直接定稿。
- [ ] 对 Reviewer 2 code-availability 的 dependency 安装困难，检查修订 release 的安装说明与可复现环境；主文无需展开排错，但 code availability 不应超出实际提供内容。

### M21 同步主图、图注、引用顺序和术语

**位置：** 全文，重点为 Figs. 1–5、L161–178、L324–340、L410–417。  
**合并意见：** R1.m2、R2.m1–R2.m8、R3.5、R3.m4。

- [ ] 对照下方主图表逐一更新图像和图注，不能只更新其中之一。
- [ ] Fig. 2d、2f 移至 Supplement，明确它们来自 benchmark truth sets，不是 Pansoma predictions；按正文首次引用顺序重排剩余 panels。
- [ ] Fig. 4 panels 按首次引用顺序重新编号；删除正文先引用后面 panel 的跳序问题。
- [ ] 更新 Supplementary Fig. 1 的 purity/VAF labels，把重叠标签移到 legend；更新 Supplementary Fig. 2 的 tumor/platform/variant type labels。
- [ ] 全文统一 Pansoma、Supplementary Fig.、Supplementary Table、read pileup、platform names、graph release/filter names，避免 Pansoma-TO、Supplement Figure 等混用。
- [ ] 修正文法并替换已失效概念：L435 的长句重写；`technical details`（L266）；`read pileup`（Fig. 5）；检查 `Pansoma converts`。删除旧 NPU expansion 的争论，因为现行 pipeline 已不输出 NPU。
- [ ] 对“已存在且正确”的部分保留框架而更新细节。例如主稿已说明 SNV/INDEL 分别评估、已存在 purity dilution、已把 unprojectable calls 单独处理；修订重点是其模型/计数/解释边界，不是重复新增相同内容。

## 主图逐一处理建议

下表区分 response 已承诺的动作与本清单为承载新增结果提出的排版建议；最终主图数量不必固定为五张。

| 当前主图 | 必须落实的变化 | 建议承载的内容 |
|---|---|---|
| Fig. 1，p3–4 | 删除 NPU 文件路径；改 cross-node tensor、三类 labels、当前训练/输出/过滤流程 | 保留 overview，但使图与最终 SNV/INDEL pipeline 一致 |
| Fig. 2，p6–7 | 五通道改八通道；去单节点截断；d/f 移 Supplement；e 更新 reference/denominator | 可整合表示方式与 channel ablation；参考 SUP Fig. R1.m2 的真实 tensor |
| Fig. 3，p9–10 | HapMap c 移 Supplement；更正 comparator 后全量更新；区分训练暴露并加 HCC1395；展示 recall | 主要 release performance 与 held-out-tumor evidence；dilution/VAF 可另列 panel/figure |
| Fig. 4，p12–13 | a 改 truth+normal+callers SBS96；d 增加 HG008T；e 区分 default/rescue；c 降低生物学结论并移 Supplement | 建议拆成“reference benefit/FP机制”和“graph-coordinate characterization”，避免一张图承载过多主题 |
| Fig. 5，p15–16 | 不把单个 pileup 例子当作新增 somatic accuracy 证明；补 allele-specific normal/tumor/CHM13 evidence 或降为表示例 | **建议**移到 Supplement，主文版面优先放三参考 ablation、candidate retention 或 FN funnel；不是 response 明确要求必须移图 |
| 当前主图未覆盖的关键结果 | 三参考 ablation；haplotype-channel ablation；FN absorption cost；difficult-region calling | 至少让三参考和 candidate-retention 核心证据进入主图或主文表；不宜全部放 Supplement |

主图和表的每一组数值均应附可追溯标识：tumor、platform、graph release/filter、model/training exposure、seed aggregation、normalization、PASS definition、population filtering、rescue mode、chromosome、benchmark region。可以把完整 manifest 放 Supplement，图注保留影响解读的字段。

## 写入主文前需要解决的来源差异

这些问题直接影响 manuscript 应该写什么。不能因为 response 使用了完成时，就默认所有参数、图表和数字已经统一。

### D01 Final INDEL 的 PASS 和 checkpoint 规则冲突

**RESP：** R3.m1 的 Calling 写 `PASS at the recall-0.9 validation threshold`，并说 INDEL 也用 validation best-F1 threshold；Checkpoint 统一写 highest validation somatic F1。  
**SUP：** Tables R1.7c、R3.m1 和 Note R1.8 对 final INDEL 一致写 **recall 0.8 after population filtering**，checkpoint 选 post-filter validation F1。

**主稿处理：** 把 SNV、earlier INDEL 和 final INDEL 拆成不同 recipe，按实际出图模型确定 operating point，再同步 response。不可在主稿采用 0.9、结果表使用 0.8 而不交代。

### D02 SNV 和 final INDEL 的训练染色体范围不同

**RESP：** R3.m1 简表对三个 CASTLE lines 统一写 all autosomes。  
**SUP：** Table R3.m1 明确 multi-tumor SNV 使用三个 CASTLE lines 的 all autosomes，而 final INDEL 使用 **every tumor 的 chr2–22**，HCC1395 另去掉 chr20。

**主稿处理：** 按 model family 列 split。既不能说全部模型在所有训练 tumors 都 held out chr1，也不能把 SNV 的 CASTLE 全染色体规则套到 final INDEL。

### D03 Graph absorption 的 normal status 和 tumor-assembly 计数未统一

以下比较均指 RESP/SUP 标识为严格 v2.1 absorption cohort 的部分，不是把 v1.1 与 v2.1 混比。

| 指标 | RESP R2.3 / R3.7 | SUP Table R2.3 / Note R2.3 |
|---|---|---|
| HG008T / COLO829T absorbed truth alleles | 1,607 / 401 | 1,607 / 401 |
| Matched normal carries ALT | 110 (6.8%) / 27 (6.7%) | 120 (7.5%) / 35 (8.7%) |
| ALT absent from both normal haplotypes | 92.2% / 92.5%（R3.7） | 1,461 (90.9%) / 363 (90.5%) |
| Normal-absent HG008T INDEL exact tumor-assembly support | 1,153 / 1,351 (85.3%) | 1,148 / 1,346 (85.3%) |
| Normal-absent COLO829T INDEL exact tumor-assembly support | 226 / 319 (70.8%) | 226 / 311 (72.7%) |

SUP Note R2.3 还说明，在 22 个 COLO829T length-matched/sequence-distinct loci 中有 7 个已重新归入 normal-carried，这提示分类曾更新，但不足以在没有最终逐 allele 记录的情况下自行解释全部差值。

**主稿处理：** 统一 cohort definition、normal-status 分类、repeat-context 匹配规则和最终 allele-level table。两个来源都支持“多数吸收 alleles 不能据此判为 truth errors”的方向；具体 percentages 不能挑一版直接抄入。Whole-truth matched-normal overlap 分析与 absorbed subset 分析也必须分开。

### D04 Default output 与 rescue-enabled output 混杂

**RESP：** R3.3/R3.4 报告 default-mode 34/541、17/287 的 tumor-only assembly support，以及 CHM13 projection 396/541、234/287，引用 Tables R3.3a/b、R3.4a/b。  
**SUP：** 目前目录和表题为 Table R3.3、R3.4；Table R3.3 明确是 **with optional off-reference rescue**，Table R3.4 同时汇总 with/without rescue。SUP Fig. R3.5b 也是 rescue-enabled calls，且含 COLO829T ONT 随机抽样。

**主稿处理：** 补齐 default/rescue 的分样本验证和投射明细，并统一工作编号。不能把 rescue-enabled 的 spectrum、validation counts 与 default 的 graph-coordinate fraction 合成一个结果；HCC1395 是否计入 pooled denominator 也需写明。

### D05 Final INDEL 的 population-feature 分支在 response 简表中遗漏

**RESP：** R3.m1 网络简表主要描述 ConvNeXt；R2.6 已承认 INDEL 模型接收 population AF input，但没有展开完整结构。  
**SUP：** Table R3.m1 及脚注列出五个 site features、转换和融合层。

**主稿处理：** 将这条输入分支写入 Methods/模型图；限定“before population filtering”和“不投射到 linear reference”的表述。这是现行模型的重要变化，不是单纯增加几条 filtering rules。

### D06 不同结果表使用的 normalization 和模型不完全一致

**SUP 的明确区别：** Table R3.9a 的 reference SNV ablation 使用 training-set statistics；Table R1.7b 的 caller/generalization comparison 使用 per-dataset base-quality standardization；Table R2.6a 的 HPRC pre-filter results 也采用 per-dataset base-quality standardization，GRCh38-only 则用 training-set statistics。Context ablation 还区分 one-GPU 与 four-GPU/global-batch controls。

**主稿处理：** reference ablation 的 matched comparison 应从对应控制实验取数，纯 filter effect 要用同一预测集。不要跨表挑数字后声称“除 graph/filter 外完全相同”。完整 protocol differences 放 Supplement，主文指出关键分析对应哪一套 protocol。

### D07 表编号与承诺的证据位置不一致

- RESP R3.1 把 paired confidence intervals 指向 Table R3.9b；SUP 的 R3.9b 实际是 final INDEL per-seed results，R3.9c 是 matched-precision recall。
- RESP 多处引用 Table R3.1，SUP 目录中没有独立同名表；目前三参考 SNV 主要数据在 Table R3.9a。需要明确对应关系，不能凭名称补出一张表。
- RESP 的 R3.3a/b、R3.4a/b 与 SUP 的 R3.3、R3.4 不完全对应。
- RESP/SUP 仍有 `Supplementary Table [X]`、`Supplementary Fig. [X]` 等占位；当前 SUP Table R1.7a 是 training-count summary，不等于已提供承诺的完整逐模型 manifest。

**主稿处理：** 先做最终 table/figure crosswalk，逐个检查正文引用能找到相同 output mode、数据集和统计定义的证据。不能只把 `[X]` 替换为任意相近编号。

### D08 部分承诺尚不能认定材料齐全

- RESP R1.m1 仍有 `[待补: multi-tumor INDEL results.]`；SUP 已有 Table R1.7c，但 response 对应文字未同步。
- RESP R3.9 仍有 `[waiting for the figure]`；当前 SUP figure inventory 没列完整 PR-curve figure。已有 matched-precision 表不等于完整 PR curves 已提供。
- SUP Table R1.8 写 `[Pending: the three ONT models and PacBio HiFi 3, retraining.]`；因此不能在主文概括全部平台的 off-path-validation-negative sensitivity 已完成。
- RESP R1.1 的 `(both?)`、R1.4 的 `v2.1?` 仍在现行文字中；需要确定 graph scope 后删除编辑残留。
- 原 Fig. 2d/f 移 Supplement 的新成图、完整 command note、model/config manifests、reference-individual/population table 等，并未都能在当前这个 SUP 文件中明确对应到成品。可能另存于其他文件，但不能据当前材料宣称全部已落地。

**主稿处理：** 不影响结论的文字占位可按已存在表补齐；缺失图表或尚待实验完成的内容保留待办，完成后再写“we show/validated”。

### D09 Dilution 的范围与“同一 reads”描述需要更精确

**RESP：** 低 expected-VAF recall 范围容易被读成三平台通用；comparator 方法总体写 same reads。  
**SUP Note R2.4：** VAF-specific ranges 来自 Fiber-seq/ONT；Illumina 在该 dilution 表中只给总体 P/R/F1；Illumina comparator mixtures 是相同来源 reads 的不同随机子集；quality centring 与 measured purity 另有明确说明。

**主稿处理：** 报告这些细节并使用对应范围，避免把 experimental sampling、platform 与 normalization 差异隐藏在“same protocol”之下。Nominal purity 的 75/50/25% 与 measured purity 不能当作同一字段。

### D10 相似百分比所用分母不同

| 容易混淆的内容 | 应明确的区别 |
|---|---|
| Truth resource 原始大小 vs 分析 truth 数 | Table R2.2 的 resource-level counts 不等于 Table R1.4 的 chr1–22 analysis counts，更不等于 chr1 benchmark-region counts |
| Candidate retention | Table R1.2 要求 representative allele 为 truth；Table R1.4 允许 any candidate allele exact match，包括非代表 allele |
| Germline graph overlap | 原 Fig. 2e 的 whole-autosome overlap 与 R2.6 benchmark-region graph-absent fractions 不能直接做“旧版→新版”的差值 |
| Absorbed alleles | Complete/partial candidate loss、strict exact-path absorption、原 v1.1 read-level cohort 是不同定义 |
| Homopolymer threshold | ≥7 bp 是部分 evaluation strata；≥8 bp 是 INDEL filtering exception，不能为了表面一致而统一成一个阈值 |
| Tumor-only assembly support | 需要说明是否已优先排除 population-represented / normal-supported classes；不是只要 tumor assembly 中存在 ALT 就算 tumor-specific |

**主稿处理：** 图注、主文和 Supplement 给出 numerator/denominator、regions、variant classes 和 grouping rule；这类区别应明确解释，不应自动当作分析错误。

## Reviewer 意见到修改主题的完整对照

同一主题可以回应多个 reviewer；正文应按科学问题组织一次，回复信再分别引用相同正文位置。

| Reviewer 条目 | 内容 | 合并到的主题 |
|---|---|---|
| R1.1 | Pangenome-native claim 和 beyond-alignment evidence | M01、M02、M03、M10 |
| R1.2 | Graph composition、size、population prior 和 sensitivity trade-off | M04、M10、M11、M12、M14 |
| R1.3 | Difficult regions 和 ONT ambiguity | M14、M15 |
| R1.4 | Graph 吸收 true somatic candidates | M05、M12 |
| R1.5 | HapMap 作为 somatic benchmark 的适用性 | M01、M04、M13 |
| R1.6 | Recall、低敏感性和 operating point | M06、M09、M12、M13、M17 |
| R1.7 | Model capacity、训练量、tumor diversity | M06、M07、M09、M17 |
| R1.8 | Projection、positive/negative label 构建 | M03、M05、M06 |
| R1.m1 | 新增 cancer lines 和独立 benchmark | M07、M09、M20 |
| R1.m2 | Single-node truncation 与 path-aware tensors | M02、M03、M10、M21 |
| R2.1 | HapMap graph donors、scope 和 limitations | M04、M13 |
| R2.2 | Truth-set construction 与错误解释 | M05、M07、M12、M20 |
| R2.3 | False-negative categories 和 matched-normal evidence | M09、M12 |
| R2.4 | Low VAF、purity、模型使用边界 | M13 |
| R2.5 | False-positive 来源、nearby germline variants | M11、M14 |
| R2.6 | Population resources、graph-absent germline 和 unfiltered HG008T | M11、M18 |
| R2.7 | Perfect-read counts/proportions 和过强 alignment 解释 | M14 |
| R2.m1 | Fig. 2d/f 的 truth 来源与移图 | M21 |
| R2.m2 | Fig. 2 panel 引用顺序 | M21 |
| R2.m3 | Supplementary Fig. 2 缺少 COLO829T 标签 | M12、M21 |
| R2.m4 | Supplementary Fig. 1 purity/VAF 标签重叠 | M21 |
| R2.m5 | Truth-overlapping 的定义 | M16 |
| R2.m6 | Fig. 4a 直接展示 truth spectrum | M16 |
| R2.m7 | Fig. 4 panel 顺序 | M21 |
| R2.m8 | Candidate site 与 candidate 定义 | M03 |
| R3.1 | Reference-level ablation | M10 |
| R3.2 | Candidate reduction 的组成和 somatic retention | M11、M12 |
| R3.3 | Graph-coordinate calls 验证 | M15 |
| R3.4 | 全输出指标范围、graph fraction、CHM13 projection | M08、M15 |
| R3.5 | Truth/SBS96 直接比较 | M16、M21 |
| R3.6 | Perfect alignment 不等于 correct placement | M14 |
| R3.7 | 撤回 benchmark-error 推断 | M12 |
| R3.8 | Tumor-only vs normal-tissue mosaicism scope | M01、M07、M13、M19 |
| R3.9 | Seeds、bootstrap、threshold robustness | M06、M08、M17 |
| R3.10 | Low-VAF 与 INDEL 结论边界 | M09、M13、M19 |
| R3.m1 | Main/Supplement 参数一致性 | M03、M05、M06 |
| R3.m2 | Comparator workflows 与公平评估 | M08、M18、M20 |
| R3.m3 | HPRC composition、filtering、sample overlap | M04 |
| R3.m4 | Terminology、NPU、语言统一 | M02、M18、M21 |

另有未编号的 code-availability remarks：Reviewer 1 肯定 README，无新增主文实验要求；Reviewer 2 提到 dependencies 阻碍运行，归入 M20 的 release/reproducibility 检查，不额外制造一个科学结果小节。

## 建议的新版正文组织

以下是组织建议，用于把重复意见合并成连贯论文；不是要求逐字采用标题。

1. **Overview of the pangenome-aware Pansoma workflow**：current GAM-based pipeline、cross-node tensor、three-class learning、output scope。
2. **Tumor datasets and evaluation design**：三个主要 benchmark、六肿瘤训练、truth construction、within-tumor 与 held-out-tumor 区别。
3. **Tumor-only SNV and INDEL performance**：release models、完整 P/R/F1、unseen-tumor results；SNV/INDEL 分开解释。
4. **Contribution of the pangenome reference and graph-derived features**：三参考 ablation、channel-6/context ablation、pre/post-filter SNV performance。
5. **Germline suppression and the sensitivity cost of graph absorption**：candidate composition、FP mechanisms、loss-stage funnel、matched-normal/tumor evidence。
6. **Performance in difficult regions and at reduced tumor purity**：region-stratified calling、alignment controls、expected-VAF sensitivity；HapMap 只简短引出 Supplement。
7. **Characterization of GRCh38-unprojectable calls**：default/rescue、CHM13 representability、assembly/read support、SBS96 描述与局限。

Methods 建议按执行顺序组织：datasets/truth → reference construction/resources → alignment/indexing → candidate/tensor generation → labels/splits → SNV/INDEL networks/training → projection/output → population filtering → comparator/evaluation → ablations/statistics/validation。

## 修改顺序和完成检查

1. **先固定结果版本。** 解决 D01–D06 的 recipe、cohort、output mode 和 normalization；整理 model/run manifest，避免写作过程中更换分母或模型。
2. **先写 Methods，再更新 Results。** M02–M08、M18 的方法事实决定各图应该采用哪组结果。把重要新方法从 Supplement 提回主文。
3. **更新主图和结果结论。** 完成 M09–M17，移走 HapMap/背景描述图，确保核心 reference ablation 和 sensitivity trade-off 在主文可见。
4. **最后改标题、摘要、引言末段和 Discussion。** 按最终实际结果落实 M01、M19，避免使用旧版“independent HG008T”“NPU innovation”“truth errors”等叙事。
5. **完成引用和可复现资料。** M20–M21；补齐 pending figures/tables，并将 response 的完成时陈述与最终稿的确切位置对应。

- [ ] 全文没有把 graph-carried allele 自动判为 patient germline 的推断。
- [ ] 全文没有把 six-tumor HG008T/COLO829T chr1 结果称为 held-out-tumor validation。
- [ ] 全文方法与 Fig. 1/2、模型配置、实际出图结果一致。
- [ ] 主结果清楚区分 SNV/INDEL、native/harmonized、pre/post-filter、default/rescue。
- [ ] 每项新增结果都有可找到的正式 table/figure/Methods 引用，已清除 `[X]` 和 pending 标记。
- [ ] 已保留真正漏检的 truth alleles，报告 recall 与 candidate floor，不用 F1 或 population-filter overlap 掩盖 sensitivity 损失。
- [ ] Graph-coordinate calls 的独立验证程度、CHM13 representability 和 benchmark coverage 边界在正文可见。
- [ ] 三位 reviewer 的 39 条编号意见均能通过上表找到正文修改位置；重复实验只解释一次，回复信分别引用。
