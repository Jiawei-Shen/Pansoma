# Pansoma 审稿回复逐点评估

**总体评价：3/5，已有实质进展，但当前版本仍需修订后提交。** 最有说服力的是承认并验证 graph absorption、区分序列匹配与定位准确性、收窄应用范围，以及公开 graph-coordinate 输出的代价。主要障碍是关键结果口径未统一、部分结论仍超出证据、少数问题尚未真正答完。现在更需要整理证据和修正论证，而不是继续增加解释性文字。

来源：[Reviewer_Letter_Response.docx](/Users/jiaweishen/Library/CloudStorage/Box-Box/Jiawei-Smaht-Paper/Pansoma_Resubmit/Reviewer_Letter_Response.docx)。以下保留原编号，压缩转述 reviewer 的要求，覆盖 **39 个编号意见及 1 个代码可用性问题**。蓝字状态和内部备注不视为已完成证据，也不作为本次任务指令。

评分针对当前回复的完整性、证据力度、自洽性和措辞边界，并非实验结果的复算认证。**4–5 分＝OK；3–3.5 分＝部分解决；1–2.5 分＝重要缺口；0 分＝缺失回复。** 总评按核心问题权重判断，不是简单平均。OK 项不列待办；文中所有图表占位符、待补标记和内部备注，统一在提交前清理，不在每点重复扣分。

## Reviewer 1

核心关切：方法是否真正利用泛基因组、是否以不可逆漏检换取精确率，以及训练与评估是否足以支持泛化。

### R1.1 泛基因组方法创新性

- **意见与应答要点：** 节点编号和图比对不足以支持 “pangenome-native”；需证明比对之后的图信息带来独立收益，或收窄主张。
- **现有回复与评分：3/5。** 改称 pangenome-aware、跨节点 tensor 和 channel 6 都切中问题。但这里仍写消融正在运行，R3.1 却已报告去掉 channel 6 的结果；跨参考比较也不能单独证明下游图特征有效。
- **还需做：** 在相同 alignment 和 candidate set 下，分别报告跨节点上下文、channel 6 的增益，一张 P/R/F1 表即可。统一已完成状态；若增益不稳定，保留较窄的方法定位，无需为了名称强行重设计网络。

### R1.2 图大小与组成的影响

- **意见与应答要点：** 解释 haplotype 数量与多样性如何影响 germline 抑制、mapping ambiguity 和 somatic absorption，尤其是更大的 HPRC 图。
- **现有回复与评分：3/5。** 三参考比较及吸收分析方向正确，但 R1.2 称过滤后 v1.1/v2.1 的 SNV F1 差异仅 −0.021 至 +0.031，R3.1 却出现 COLO829T Illumina 相差 0.163。两个版本同时改变规模、组成与构建细节，不能把差异全归因于规模。
- **还需做：** 用同一最终结果表统一 seed、阈值、rescue、区域和过滤口径；引用 R3.6 的实测 ambiguity，替代以既有文献概括当前图的表现。将结论限定为两个 release 的比较；只有坚持“规模的因果作用”才需受控子图实验。

### R1.3 困难区域的真实收益

- **意见与应答要点：** 在 repeats、paralogs、低可比对区等分层验证收益，尤其区分 ONT 新信号与错误比对。
- **现有回复与评分：4.5/5，OK。** 提供区域分层、已知来源模拟和 calling precision；明确承认 v2.1 的 ONT ambiguity 及 rescue 代价。结论与证据边界匹配，没有把更少 mismatch 当作更正确的定位。

### R1.4 已存在于图中的真实体细胞变异被漏掉

- **意见与应答要点：** 正面处理 edit-based candidate generation 的盲区；群体中存在的 ALT 仍可能是该患者的 somatic allele，不能改判为 benchmark 错误。
- **现有回复与评分：2.5/5。** 已承认机制并给出定量损失，是实质进展。但“后续过滤本来也会删除”只能证明两类损失重叠，不能证明候选阶段漏检无害。R2.3 的 normal/tumor 证据反而支持多数被吸收位点是真实漏检；“every tumor-only caller shares”过强。
- **还需做：** 并列给出完整 truth 分母下的吸收损失、过滤重叠和过滤后额外损失，区分全基因组与 chr1。删去普遍化辩护，保留 comparator 原生过滤结果。若不改候选策略，就明确这是当前方法不可恢复的敏感性代价，不能称机制已解决。

### R1.5 HapMap 的生物学有效性与 donor overlap

- **意见与应答要点：** 解释为何 germline 混合物可用于本研究，核实 donor 排除，并避免将其当作真实低 VAF somatic benchmark。
- **现有回复与评分：3/5。** 专用 donor-excluded 图、移至补充材料、撤回 somatic model 定位都合理。但额外排除 graph-represented truth 后的 recall 是条件召回率，可能使图方法看起来更好。
- **还需做：** 给出原始 truth、排除数量及剩余 truth；若比较工具，必须使用同一分母，并明确全量 truth 的敏感性。把实验定位压缩为“限定条件下的低频等位基因检测”，真实肿瘤低纯度结果集中回答在 R2.4。

### R1.6 召回率尤其 INDEL 召回率过低

- **意见与应答要点：** 独立评价 sensitivity，给出漏检环节、匹配 precision 的比较和实际使用边界，不能仅用 F1 抵消低 recall。
- **现有回复与评分：3.5/5。** PR 曲线、first-loss funnel、修订前后结果很有价值。但 INDEL 的 0.342/0.382 上限来自所选过滤规则，不是 tumor-only calling 的不可突破上限；部分 candidate retention 还计入 partial matches。
- **还需做：** 对完整 truth 同时呈现过滤前后 recall，并保留 comparator 原生结果；按统一 exact-allele 标准给出可实现的候选召回上限，partial 单列。把“达到过滤策略的上限”与“已解决敏感性不足”明确分开。

### R1.7 网络容量与肿瘤多样性

- **意见与应答要点：** 报参数量、正负样本数和过拟合证据，并通过不同突变谱、独立肿瘤验证 generalization。
- **现有回复与评分：3.5/5。** 缩小模型、数据量曲线及四肿瘤外部测试回答了主要质疑，也坦诚报告 UV 偏倚。但 multi-tumor INDEL 尚未补齐；加入目标肿瘤 chr2–22 后的 chr1 提升不属于独立肿瘤泛化。
- **还需做：** 补齐各类别训练数量及 multi-tumor INDEL 结果，差结果也应报告；将 within-tumor、held-out-tumor、加入目标肿瘤后的结果分开。默认模型的适用结论必须以真正未见肿瘤结果为依据。

### R1.8 标签与训练集构建

- **意见与应答要点：** 说明投射、INDEL 等价表示、graph-only/ambiguous 位点、正负例定义，以及 VAF/purity 分布是否导致选择偏倚。
- **现有回复与评分：3/5。** 精确匹配、重复区等价位置、somatic/germline/other 三类和排除规则相当具体。但此处称 off-path tensor 全部排除，R3.4 却称仅 93.6–99.7% somatic 例子在 GRCh38 path；R2.2 也保留“off-reference 无 somatic match 才排除”的说法。
- **还需做：** 以一个标签规则表统一这些例外及类别优先级，列 exact、partial、未生成候选、不可标注的数量。说明 matched-normal truth 中采用 partial germline 匹配、somatic 却要求 exact 的影响，并明确结果主要适用于被纳入训练的表示范围。

### R1.m1 增加 HCC1395 等 benchmark

- **意见与应答要点：** 用更多可靠 tumor truth 检验结果，不局限于 COLO829T/HG008T。
- **现有回复与评分：3.5/5。** 新增 HCC1395 和三种 CASTLE 细胞系，并解释 CASTLE 不适合作为独立 benchmark，理由合理；INDEL 仍为空缺。
- **还需做：** 补 HCC1395 的 SNV/INDEL comparator 表和训练暴露关系；将其余三种明确称为训练扩展或探索性分析，避免概称四个独立 benchmark。无需为凑数量强行纳入不可靠 truth。

### R1.m2 tensor 在节点边界被截断

- **意见与应答要点：** 让短节点、单碱基节点也保留相邻节点和 read path 的上下文。
- **现有回复与评分：4/5，OK。** 明确说明 101 列窗口沿 read walk 跨节点，并给出 ALT/REF reads 经过不同分支的具体例子，已回答表示层面的缺陷。性能归因由 R1.1 统一处理。

Reviewer 1 对 README 的肯定不需要补充实质答复。

## Reviewer 2

核心关切：benchmark 的合理性、precision/recall 改变的原因，以及低 VAF 与过滤策略的实际限制。

### R2.1 HapMap 是否适合模拟 somatic mutations

- **意见与应答要点：** 解释混合物用途、图中 donor/allele 重叠及生物学局限。
- **现有回复与评分：3.5/5。** 已区分 donor overlap 与共享 allele，并撤回对 somatic discrimination 的外推。仍欠排除 graph-represented truth 后的明确分母。
- **还需做：** 与 R1.5 共用一份排除统计及评估定义；本条保留用途与局限即可，不必再次长篇解释建图。

### R2.2 truth set 来源及所谓 truth 错误

- **意见与应答要点：** 简述各 truth 的生成和验证方式；只有独立证据才能支持 germline 污染或 artifact 的判断。
- **现有回复与评分：3/5。** 给出了各来源，并撤回无证据的 artifact 判断，但大量训练流程偏离问题。H1437 的 24.5%、27.8%、约 30% 在不同回复中缺少分母解释；“97.5%…chr1:9–12 Mb”和“chr1 仅 19/4,266”也需核对适用范围。
- **还需做：** 改为简短 truth 来源表，列版本、生成/验证方法、区域和已知限制；核实 H1437 数字与区间。将“正常样本也携带 ALT”限定为待解释的 genotype/LOH 歧义，不直接等同整套 truth 错误。

### R2.3 漏检分类及 matched-normal 验证

- **意见与应答要点：** 解释旧的两个不清楚的漏检标签，并用正常样本检验 graph-overlapping alleles 是否真是 germline。
- **现有回复与评分：4.5/5，OK。** 撤回旧分类，采用 first-loss funnel，补 normal reads、双单倍型组装及 tumor assembly；区分 exact ALT 与仅长度相同，并承认多数吸收位点不是 truth 错误。这是全文最有说服力的回应之一。

### R2.4 真正的低 VAF 检测能力与模型选择

- **意见与应答要点：** 回答标准模型何时失效、低 VAF 模型如何不同、是否能救回真实 somatic variants，以及何时使用哪个模型。
- **现有回复与评分：2.5/5。** VAF 分层和 0.08 门槛解释清楚，但主要恢复证据仍来自 HapMap；推荐 multi-tumor 模型并不能绕过候选门槛。先描述 500-read tensor，随后用 200-read 重跑结果，也容易混淆模型身份。
- **还需做：** 用真实 COLO829T dilution 汇总 purity/VAF、truth 数、candidate recall 和最终 P/R/F1，标清是否同肿瘤训练；统一模型与 tensor 设置。若无有效的 tumor low-VAF 模型，就直接说明当前不支持该范围，避免以 HapMap 高 precision 代答。

### R2.5 precision 改善来自哪里

- **意见与应答要点：** 区分 germline contamination 和其他假阳性；检查邻近多个 germline variants 的区域是否受益。
- **现有回复与评分：4.5/5，OK。** 同参考对照、germline-associated FP、困难区分层和 ±50 bp 聚集分析均直接回答问题；也说明不同类别重叠，未冒充独立因果贡献。图与正文一致。

### R2.6 population filtering 是否泄漏或掩盖图的局限

- **意见与应答要点：** 检查资源是否包含测试样本；说明图未覆盖的 germline 如何被清除，并给 HG008T 无过滤结果。
- **现有回复与评分：2.5/5。** 后两问已有数字，但 rare-allele overlap 较低不能证明 gnomAD/dbSNP 完全没有样本贡献，更不足以推出“benchmark does not overstate its effect”。HapMap 段仍待补。
- **还需做：** 将结论改为“未发现明显重叠证据，不能完全排除”；把 sample-list 核实与遗传间接检验分开。补齐 HapMap 段，区分 truth 预筛选与 caller 后过滤，并用一张表固定 HG008T 过滤前后结果。

### R2.7 perfect alignment 与 clipping 的表述夸大

- **意见与应答要点：** 提供实际 read 数或比例，不仅报倍数改善。
- **现有回复与评分：4.5/5，OK。** 已给绝对计数、分母和百分点变化，说明各参考并非完全相同的 read 集，且不同平台的 clipping 并非统一改善；对“序列更吻合”的定位也准确。

### Reviewer 2 的 minor comments

| 编号 | Reviewer 意见与应答要求 | 现有回复评价及评分 |
|---|---|---|
| R2.m1 | 说明 Fig. 2d/f 的变异来自 Pansoma 还是 truth；若只是数据集描述，可移到补充。 | 明确来自已发表 truth，并移到补充、注明来源。**4.5/5，OK。** |
| R2.m2 | 修正 Fig. 2f 先于 2e 引用的问题。 | 已按首次引用重排正文与图例。**5/5，OK。** |
| R2.m3 | Supplementary Fig. 2 的 COLO829T 标识不清。 | 改为显式标注 tumor、platform、variant type。**5/5，OK。** |
| R2.m4 | Supplementary Fig. 1 的 purity/VAF 标签遮挡数据点。 | 标签移至图例，回应直接。**5/5，OK。** |
| R2.m5 | 解释 “truth overlapping”。 | 改为 exact-allele true positives，并独立定义三工具 shared core。**5/5，OK。** |
| R2.m6 | Figure 4 直接展示 truth 的 signature composition。 | 用 truth 与各 caller 并列 SBS96 代替分解柱图，说明了改法，已满足直接比较目的。**4.5/5，OK。** |
| R2.m7 | Figure 4 面板顺序与正文引用不一致。 | 已声明按首次引用顺序统一。**5/5，OK。** |
| R2.m8 | 解释 candidate sites 与 candidates 的差别。 | 区分 allele 与按位置/变异类别分组的 tensor，并解释新版计数单位。**4.5/5，OK。** |

### R2.code 安装依赖阻碍运行

- **意见与应答要点：** Reviewer 因 dependencies 未能完成测试；需要可实际运行的安装路径。
- **现有回复与评分：0/5。** 原文该意见后直接进入 Reviewer 3，没有独立回应。
- **还需做：** 增补简短回复，提供固定版本环境或容器、最小输入与一条运行命令、预期输出，以及干净环境中验证成功的记录。不能用“README 完整”替代解决安装失败。

## Reviewer 3

核心关切：用受控分析证明泛基因组本身的贡献，并完整评估方法输出、统计稳定性及可复现性。

### R3.1 隔离泛基因组贡献的消融

- **意见与应答要点：** 在同一 pipeline 下比较 GRCh38-only 与 HPRC，保持模型、训练、候选和过滤条件可比。
- **现有回复与评分：3.5/5。** 三参考、三平台、三 seed 和无 channel 6 重训是正确设计；但结果与 R1.2、R3.9 的摘要口径不统一。“全部下游条件相同”也需要更精确。
- **还需做：** 固定一张逐 dataset 的主结果表，标明 mean 或单 seed、CI 和各阈值。核实不同图的 node-block 划分是否对应相同基因组训练/验证区间；区分相同阈值选择规则与相同数值阈值。不能把图 node ID 相同或 seed 相同自动视为相同数据划分。

### R3.2 candidate reduction 是否选择性保留 somatic

- **意见与应答要点：** 同时统计 somatic、germline、background，并按 germline allele 是否被图表示，关联候选减少与最终 FP 减少。
- **现有回复与评分：3/5。** SNV 三类候选分析较充分，但主要报告总体图覆盖率，尚未完全给出 reviewer 要的 represented/unrepresented 分层链条；INDEL 的真实候选保留率也没有对应呈现。
- **还需做：** 在同一表中补 SNV/INDEL 的三类 retention，并将 germline 按图内/图外分层，列出 candidate 与最终 FP 数。不能用 SNV 的高保留率概括所有 somatic variants。

### R3.3 graph-coordinate calls 的正交验证

- **意见与应答要点：** 用 normal、phased assembly 或跨平台证据验证代表性 graph-exclusive calls，不能仅凭 signatures 和个例。
- **现有回复与评分：4.5/5，OK。** 提供组装、normal reads、跨平台及阳性/阴性对照；区分自有 reads 组装的一致性支持与独立验证，报告随机抽样和未证实类别，也收窄了生物学结论。无需为本条再扩展一套实验。

### R3.4 只评估可投射子集导致整体指标不完整

- **意见与应答要点：** 报完整输出中的 graph-coordinate 比例及其代价，修正“约 2%”，并检验 CHM13 可表示比例。
- **现有回复与评分：4.5/5，OK。** 已报告 rescue 开关、保守 FP 假设、性能反转和 CHM13 投射，并说明保守 F1 不是真实全输出性能、不测量参考外 recall。也没有把不能投射直接等同“只能由图表示”。

### R3.5 truth 与 graph-coordinate SBS96 展示

- **意见与应答要点：** 将 truth 与 callers 并排，提供 graph-coordinate 的 SBS96，并检验其是否真的支持 somatic interpretation。
- **现有回复与评分：3/5。** 撤回以谱相似性证明 somatic 的推论非常正确，但有可直接核实的图文错配：**R3.5b 图注说 graph-coordinate calls，实际图是 COLO829T 的 truth/normal/三工具常规 benchmark SBS96 比较**，图题还写着 population filtering 和 benchmark regions。
- **还需做：** 将两张现有 tumor 对比图归入 R3.5a，补入真正的 graph-coordinate 图至 R3.5b。统一“22–66”究竟是全基因组 tensor-level detections 还是去重后的 variants，避免与 chr1 的 2–8 个支持位点混用。

### R3.6 perfect alignment 不等于 correct alignment

- **意见与应答要点：** 用 GRCh38、CHM13、HPRC 及定位真值/MAPQ/困难区分析，区分参考完整性与多单倍型贡献。
- **现有回复与评分：4.5/5，OK。** 四参考模拟、MAPQ 校准、ONT 反例及 SV reads 未评分比例都透明；明确共享可评估区不能充分验证 CHM13 完整性，结论没有越界。回复偏长，但科学上成立。

### R3.7 不应把 false negatives 直接改判为 truth 错误

- **意见与应答要点：** 补独立证据，或撤回 Supplementary Fig. 3 中过强的解释。
- **现有回复与评分：3.5/5。** 撤回原解释、承认多数为真正漏检，核心回答充分。但此处说 COLO829T “no tumor assembly is available”，R3.3 明确使用了其 hifiasm tumor assembly。
- **还需做：** 明确是“本项分析未使用”或“没有独立 tumor assembly”，依据实际情况统一，不能同时称不存在和已使用；标明此处计数采用 v1.1 的原始 read-level 定义。

### R3.8 正常组织 mosaicism 的适用性

- **意见与应答要点：** 若主张广义 somatic calling，需真实正常组织验证；否则全文明确限定 tumor-only cancer calling。
- **现有回复与评分：5/5，OK。** 选择收窄范围，明确未经正常组织训练/验证、0.08 门槛及 dilution 的局限，已满足 reviewer 给出的替代路径。无需为当前范围新增正常组织实验。

### R3.9 seed 稳定性与统计不确定性

- **意见与应答要点：** 提供重复训练、基因组区块 CI、与 comparator 的差值 CI，以及阈值敏感性；区分同肿瘤留出与独立 HG008T。
- **现有回复与评分：3.5/5。** 三 seed、paired bootstrap、PR 曲线及未显著的 ClairS-TO 比较均有价值。但独立跨肿瘤测试只有一 seed，不能承接同肿瘤分析的稳定性结论；与 R3.1 的 ΔF1 范围也缺口径说明。
- **还需做：** 明确 bootstrap 使用哪个 seed 或如何汇总，并与主表一致；为核心独立 HG008T 模型补多 seed，或将跨肿瘤稳健性主张降为单次探索结果。seed 同时改变验证划分时，应称联合训练/划分波动，不能解释为仅初始化波动。

### R3.10 低 VAF 与 INDEL 的结论边界

- **意见与应答要点：** 用真实 somatic dilution 支撑低 VAF，按 insertion/deletion、长度、重复背景和平台分层 INDEL，并收窄性能主张。
- **现有回复与评分：3.5/5。** 分层结论及低 TP 数限制合理，但关键 dilution 数据主要通过交叉引用带过；称 R1.5 已给 tumor/normal base-quality controls，实际 R1.5 没有这些内容。
- **还需做：** 给出可定位的真实 dilution 与 INDEL 分层表，至少有各层 truth/TP/FP/FN。修正错误引用；若 quality control 尚未完成，不应写成已完成。低 VAF 主结果与 R2.4 共用，不再重复叙述。

### R3.m1 方法参数不一致

- **意见与应答要点：** 统一 AF、ALT reads、BQ、MAPQ、INDEL 长度与 epochs，并充分说明上下游步骤。
- **现有回复与评分：3/5。** 配置表很有用，也诚实披露看过 chr1 后检查早期 checkpoints。但这与其他回复笼统宣称“从未依据 chr1 选择”冲突；R1.7 的统一 chr2–22 训练描述与此处 CASTLE 使用全部 autosomes 也需对齐。
- **还需做：** 用每个实验的配置与 split manifest 作为唯一依据；主结果只使用预定验证规则，chr1 后验 epoch 结果独立标作探索，不能替换主消融。同步 MAPQ 的 > 与 ≥、模型 epochs、候选与 tensor 单位。

### R3.m2 comparator 配置与公平性

- **意见与应答要点：** 列版本、aligner、checkpoint、tumor-only 参数、过滤、资源和命令，确认比较公平。
- **现有回复与评分：3.5/5。** 信息已相当完整，且修正了 HG008T ONT 曾误用 PACBIO_TUMOR_ONLY 的问题。但在各工具原生 PASS 后再叠加相同过滤，不等于各工具拥有相同总过滤规则。
- **还需做：** 明确分开 native performance 与 harmonized-filter 分析；后者若用于公平归因，应说明原生 population tags 如何处理，并保留技术质量过滤。补 ClairS-TO 实际 INDEL calling BED、位点/allele 匹配规则及训练暴露；确认所有表和 PR 曲线均来自修正后的 ONT 模型。

  官方文档支持这一核查重点：ClairS-TO 的默认资源包含不同的 allele-matching 规则，并有 INDEL 区域设置；DeepSomatic 1.9 的 population-frequency 输入与可开关的后过滤是两个环节。这里只据此提出配置核查，不推定实际运行参数。[ClairS-TO 文档](https://github.com/HKU-BAL/ClairS-TO)；[DeepSomatic 1.9 ONT 文档](https://github.com/google/deepsomatic/blob/r1.9/docs/deepsomatic-case-study-ont-tumor-only.md)。

### R3.m3 图版本与组成

- **意见与应答要点：** 明确 release、haplotypes、AF 过滤、SV 处理、构建和 benchmark donor/亲缘重叠。
- **现有回复与评分：4.5/5，OK。** v1.1/v2.1、d9/d46、参考路径保留、SV 范围及专用 HapMap 图均有交代；对 rare-allele sharing 不能正式排除亲缘的限制表述恰当。

### R3.m4 术语与语言统一

- **意见与应答要点：** 区分 PoN 与 population resources，统一 NPU、方法名、图表名称并修正语言。
- **现有回复与评分：4/5，OK。** 改用 population filtering resources，说明新版不再生成 NPU，并逐项回应命名和语法。方法大幅更新后删除过时 NPU 描述是必要且合理的处理。

## 控制回复长度的建议

**同意需要压缩。** 主要冗余来自同一分析在多位 reviewer 之间反复完整展开，而不是审稿问题本身太多。

- 每个主要回复保留“直接回答、最关键结果、限制与修改位置”；普通点以约 100–180 个英文词为目标，核心方法争议允许更长，minor 一两句即可。
- absorption 证据统一放 R2.3，reference ablation 放 R3.1，低 VAF 放 R2.4，graph-coordinate 验证/完整输出放 R3.3–R3.4。其他点保留针对性结论及交叉引用，避免再次复述整套数据。
- 长表、算法细节和多组数值移至正式 Supplementary Note/Table；回复只留足以改变 reviewer 判断的数字。保留不利结果和适用限制，优先删除背景介绍、反复感谢及防御性论述。
- 先统一最终实验版本，再压缩文字。尤其是 R1.2/R3.1 的图比较、R1.8/R3.4 的 off-path 标签、R3.3/R3.7 的组装状态和 R3.5b 图文对应；仅润色会把这些实质问题藏得更深。
