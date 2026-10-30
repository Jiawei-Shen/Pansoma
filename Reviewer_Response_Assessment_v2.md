# Pansoma 第二版审稿回复复评

**总体评分：3 → 3.5/5。新版明显进步，但仍不建议原样提交。** 40 项中，31 项的回应已足够，9 项仍需处理；19 项上调、2 项下调。上调主要来自新增实质证据和更准确的结论边界。剩余问题集中在消融的可比性、标签规则、数字对应关系和代码发布状态，不需要再全面扩展实验。

评估对象：[第二版回复](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md)及随附的 6 张图；对照[第一次评估](/Users/jiaweishen/Documents/GitHub/Pansoma/Reviewer_Response_Assessment.md)。已通读正文和补充表，检查全部图片，并复核部分表格算术与跨表对应。目录未包含文中另指的完整 Supplementary Data、配置 manifests 或原始运行结果，因此本次仍评估回复的说服力与自洽性，不作实验复现认证。

沿用原标准：**4–5＝OK；3–3.5＝部分解决；1–2.5＝重要缺口；0＝缺失。** 各项显示“上次 → 本次”；总评按核心问题权重判断，不作简单平均。OK 不等于方法没有局限，而是已经充分回答该条质疑。只有未达标项列“还需做”。内部待办不视为完成证据。

## 本轮最重要的变化

已经解决：R3.5b 图文错配、COLO829T 是否有 tumor assembly 的矛盾、H1437 的 chr1/chr9q 区间错误、INDEL 结果缺失，以及将过滤上限误当作方法上限的论证。真实 somatic dilution 和 HapMap 模型转移失败的结果也补得很好。

优先处理以下四组即可：**R1.8 的标签规则；R1.1/R3.1 的消融可比性；R1.5/R3.9 的分母与模型对应；R2.code 的实际可获取版本。** 其中代码回复末尾的待确认是发布状态问题，不能只删掉备注就视为完成。

## Reviewer 1

### R1.1 泛基因组在比对之后的贡献

- **Reviewer 要求／回应应覆盖：** 证明下游图信息带来可测量收益，或收窄 pangenome-native 的创新性主张。
- **评分：3 → 3.5/5。** Channel 6 重训和 AP/P/R/F1 是实质补充，已经超越“换成节点坐标”的解释。但表中混用单次、三 seed 均值及部分单 GPU controls；跨节点 context 消融仍待补。三 seed 均值不是单次消融的严格配对 control。
- **还需做：** 给每对 with/without 标明相同 seed、checkpoint 规则、归一化、GPU/全局 batch；已有匹配运行就直接整理，无需一律重跑。说明 AP 是过滤前 tensor-level 还是过滤后 call-level。跨节点实验完成后再称其收益“已测量”；若不做，限定为表示方式改进。特别核对 v1.1 COLO829T Illumina 的 control F1 **0.779** 与主消融 **0.572** 的实验身份，不能让读者自行猜测。

### R1.2 图大小和组成

- **Reviewer 要求／回应应覆盖：** 比较图 release、germline 抑制、somatic 损失与 mapping ambiguity，区分规模和组成。
- **评分：3 → 4.5/5，OK。** 已明确是 release 比较，不能单独归因于 size；v2.1 的 0.163 例外与主表对齐，新增实测 ambiguity，并报告更大的 INDEL 损失。上次主要问题已解决。

### R1.3 困难区域与 ONT ambiguity

- **Reviewer 要求／回应应覆盖：** 在困难区分层评估实际 calling 收益，排除“更干净但定位错误”的解释。
- **评分：4.5 → 4.5/5，OK。** 新版补充困难区 F1、区间和 INDEL 反例，比只报 precision 更完整；保留 ONT v2.1 和 rescue 的不利结果。已足够回答本条。

### R1.4 graph absorption 导致真实 somatic 漏检

- **Reviewer 要求／回应应覆盖：** 承认已有图路径不代表患者 germline；量化不可恢复的候选损失，不能把它改判为 truth 错误。
- **评分：2.5 → 4/5，OK。** 现在明确区分完整 truth、过滤重叠、strict absorption 和更大的 partial absorption 损失，并承认 **3.5–7.0% 的 INDEL 额外损失不能恢复**。Comparator 原生结果也保留了。“不致命”仍可能不被 Reviewer 1 接受，但回应已形成可审查的取舍论证；不要求为了拿高分必须解决整个算法盲区。

### R1.5 HapMap 适用性与 donor overlap

- **Reviewer 要求／回应应覆盖：** 解释 donor 排除、生物学限制，以及 graph-represented truth 被移除后真实的评估分母。
- **评分：3 → 3.5/5。** 专用图、full/conditional recall 和 Mutect2 同分母比较解决了主要解释问题，但表内排除流程尚未闭合：SNV 的 `392,567 − 175,773 = 216,794`，表中剩 **215,094**；INDEL 应剩 **98,561**，实际为 **92,079**。
- **还需做：** 列出额外排除的 **1,700 SNV、6,482 INDEL** 及原因。可能是按 site masking 连带移除其他 alleles，而非算错，但需要明确。用最终分母核对两工具的 recall，说明 full truth 的 TP 是否重新计算。[定位表格](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:126)

### R1.6 低 recall 和 INDEL sensitivity

- **Reviewer 要求／回应应覆盖：** 报绝对 recall、漏检环节、matched-precision 比较和适用限制，不能用 F1 替代 sensitivity。
- **评分：3.5 → 4.5/5，OK。** Exact 与 partial candidate 已分开，提供完整 truth 分母及 comparator 原生结果；也坦承最终 INDEL operating recall **0.11–0.18 并未优于投稿版**。过滤规则上限不再被当作 tumor-only 的固有限制，回应充分而诚实。

### R1.7 模型容量和跨肿瘤泛化

- **Reviewer 要求／回应应覆盖：** 参数和样本数量、容量合理性、不同突变谱及真正未见肿瘤的测试。
- **评分：3.5 → 4/5，OK。** 多肿瘤 INDEL 已补齐，W/H/T+ 明确区分训练暴露；单 seed 的跨肿瘤 SNV 被限定为探索性，UV 转移失败及 INDEL 不敌 ClairS-TO 均直接报告。不再把加入目标肿瘤后的提升当作外部泛化。

### R1.8 标签和训练集构建

- **Reviewer 要求／回应应覆盖：** 给出可唯一执行的投射、等价表示、标签优先级和排除规则，并说明选择偏倚。
- **评分：3 → 2.5/5。** 原先“somatic 正例是否 off-path”的矛盾已纠正，但新增规则暴露了更实质的问题，所以本项下调：表中 **exact germline 在 partial somatic 前**，后文却说 partial somatic 优先排除；此外，无 truth 的 off-path tensors 被直接当作 validation/test negatives。无标签不等于已证实非 somatic，可能影响 checkpoint 和阈值选择。
- **还需做：** 明确同时 exact-germline/partial-somatic 时的唯一处理顺序；分别列训练、验证、benchmark 的规则。无证据的 off-path 位点应作为 unknown 排除于主要标签评分，或明确作为保守假设并核查是否改变模型选择；不要与经 normal 证实的 negatives 合并。解释 83–717 个 off-path INDEL negatives 的例外条件。[定位规则表](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:188)

### R1.m1 增加独立 benchmark

- **Reviewer 要求／回应应覆盖：** 增加可靠癌症数据，交代 truth 来源和各工具训练暴露。
- **评分：3.5 → 4.5/5，OK。** HCC1395 的 SNV/INDEL 表完整，CASTLE 明确降为训练与探索用途；没有回避 ONT 上落后。既扩展了证据，也避免了把非独立测试包装成新 benchmark。

### R1.m2 节点边界截断

- **Reviewer 要求／回应应覆盖：** 展示短节点也能保留 read path 和跨节点上下文。
- **评分：4 → 4.5/5，OK。** 实际 tensor 图与 38 ALT／40 REF reads、21 bp 外分支的描述对应，表示层面的修复有直接图示。性能增益是否已测量仍由 R1.1 处理，不在此重复要求实验。

R1.code 是对 README 的肯定；新增简短答谢恰当，不额外计入 40 项评分。

## Reviewer 2

### R2.1 为什么保留 HapMap

- **Reviewer 要求／回应应覆盖：** 解释低频检测用途与真正 somatic discrimination 的区别。
- **评分：3.5 → 4/5，OK。** donor 排除、共享 alleles、无 population filtering 和限定到 Supplement 的理由均清楚。本条的科学定位已解决；具体分母表的缺口集中记在 R1.5。

### R2.2 truth set 如何建立及其局限

- **Reviewer 要求／回应应覆盖：** 简明列出 truth 的来源、生成和验证方法，对疑似 germline/artifact 使用证据约束。
- **评分：3 → 4.5/5，OK。** 来源表切题；H1437 的 24.5% 图频率、27.8% gnomAD 频率、约 30% LOH 区间比例已有不同定义，区间改为 chr9q，与 chr1 的少量异常不再冲突。未逐个验证的事件也不再称为已证实错误。

### R2.3 漏检原因及正常样本证据

- **Reviewer 要求／回应应覆盖：** 替换不清楚的漏检标签，并用正常样本检验吸收位点的真实状态。
- **评分：4.5 → 4.5/5，OK。** Stage funnel、normal reads、正常双单倍型和 tumor assembly 相互补充；明确两版本的吸收定义不同，也区分同长度与同 allele。继续保留这段作为 absorption 证据的主要出处即可。

### R2.4 低 VAF 能力与模型选择

- **Reviewer 要求／回应应覆盖：** 标准模型的失效范围、低 VAF 模型的区别和实际效果、何时用哪个模型。
- **评分：2.5 → 4.5/5，OK。** 真正补上了 dilution 的 candidate recall 和最终 P/R/F1，并直接报告 HapMap 模型转到 tumor 后的低 precision 和 clonal recall 失败。500/200 reads 的历史与当前设置也清楚。Illumina comparator 使用不同抽样，故这张表适合描述趋势，不宜解读为严格配对的微小差异优势；本条关于低 VAF 限制的回答已成立。

### R2.5 precision 的来源

- **Reviewer 要求／回应应覆盖：** 区分 germline 抑制、困难区和邻近 germline 聚集的贡献。
- **评分：4.5 → 4.5/5，OK。** 候选、最终 FP 和局部 read 差异都有对应分析；新增 clustered germline 附近 INDEL retention 很差的反例，使机制解释更完整。图中区间也明确不是 CI。

### R2.6 population resources 的重叠与作用

- **Reviewer 要求／回应应覆盖：** 样本泄漏风险、图外 germline 的去除过程、HG008T 无过滤结果。
- **评分：2.5 → 4.5/5，OK。** 已分开 sample-list 核查和遗传间接证据，不再声称完全排除 gnomAD/dbSNP 重叠；HapMap 段补齐，明确不筛 truth 或 calls。HG008T 分平台过滤前后表也直接回答了原问题。

### R2.7 alignment 统计是否夸大

- **Reviewer 要求／回应应覆盖：** 报实际 read 数与比例，解释不同平台的效果和 denominator。
- **评分：4.5 → 4.5/5，OK。** 四参考的完整计数表与图已随附；抽核比例与计数相符，保留不同 read 集及 split alignment 的口径说明。没有继续混淆序列吻合与定位准确性。

### Reviewer 2 的 minor comments

| 编号 | Reviewer 要求／应答内容 | 本轮评价与评分 |
|---|---|---|
| R2.m1 | 说明 Fig. 2d/f 的数据来源，必要时移至 Supplement。 | 明确为 truth 描述并移走。**4.5 → 4.5/5，OK。** |
| R2.m2 | Fig. 2 面板按正文引用排序。 | 回应明确。**5 → 5/5，OK。** |
| R2.m3 | Supplementary Fig. 2 明确 COLO829T 等标签。 | 声明按 tumor/platform/type 标注，足以回应此类修改。**5 → 5/5，OK。** |
| R2.m4 | purity/VAF 标签不要遮挡数据点。 | 移至图例。**5 → 5/5，OK。** |
| R2.m5 | 定义 truth overlapping。 | 改为 matching truth SNV 的 true positives，另定义 shared core。**5 → 5/5，OK。** |
| R2.m6 | 把 truth 谱直接放在 caller 旁。 | 两种 tumor 的实际图均已核对。**4.5 → 4.5/5，OK。** |
| R2.m7 | Figure 4 面板引用顺序一致。 | 回应明确。**5 → 5/5，OK。** |
| R2.m8 | 区分 candidate allele、site 和 tensor。 | 给出新版计数定义。**4.5 → 4.5/5，OK。** |

### R2.code 安装依赖与可运行性

- **Reviewer 要求／回应应覆盖：** 提供审稿人现在可以获取并成功运行的安装方案。
- **评分：0 → 2.5/5。** 新回复给出了具体失败原因、CPU/GPU 环境、最小例子和运行记录，方向很好。但最后明确注明：须等 develop 合并、文件提交、example/checkpoints 发布之后才成立。它目前证明了本地验证的报告，未证明 reviewer 能复现。
- **还需做：** 发布可获取的固定 release/commit 和数据、模型入口，附一条实际命令及预期输出，再从该发布版本复测。无需一定合入 main，固定的可访问 release 也可以；关键是所述资源真实可取。发布前删除或改写“now provides／both are fixed”等已完成表述。[定位待确认项](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:434)

## Reviewer 3

### R3.1 隔离泛基因组贡献

- **Reviewer 要求／回应应覆盖：** 在可比训练和评估条件下，只改变参考，区分 graph 效果与其他 pipeline 改变。
- **评分：3.5 → 3.5/5。** 现在说明 mean 与 first-seed CI 的区别，后验 epoch 也独立标注；旧的大范围数值矛盾基本解决。但各图 validation 区间只有 **2–6% overlap**，会同时改变训练数据、checkpoint 和阈值；“necessarily differ／exactly this ablation”仍把未消除的混杂说得太轻。
- **还需做：** 若保留严格隔离贡献的主张，用相同 GRCh38 区间投射得到共享训练/验证划分，并作关键配对核查；若不补实验，将结论限定为“各参考按同一规则独立训练的系统级比较”，明确划分混杂。主表与 CI 的具体模型对应另见 R3.9。[定位实验设计](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:458)

### R3.2 哪些候选被图去除

- **Reviewer 要求／回应应覆盖：** 同时跟踪 somatic/germline/background，按图是否表示该 allele 分层并连接最终 FP。
- **评分：3 → 3.5/5。** 新表已补 SNV/INDEL 三类及 FP，且承认 somatic INDEL 并未被选择性保留，实质缺口已明显缩小。剩余问题是把“≥46/463 haplotypes 携带”直接定义为图内表示；R3.m3 又说明 GRCh38/CHM13 路径不受该过滤约束，两者不完全等价。
- **还需做：** 若用实际 graph-path 检查，写清判定；若用 carrier-frequency proxy，将列名改为相应频率类别，并说明保留参考路径的例外。机制结论不能把“高频群体 allele”与“图确实能够拼出该 allele”当作同一事件。[定位分层定义](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:496)

### R3.3 graph-coordinate calls 的验证

- **Reviewer 要求／回应应覆盖：** 用 assemblies、normal 或跨平台证据验证代表性 calls，而非只靠 signature 和个例。
- **评分：4.5 → 4.5/5，OK。** 三类证据、随机抽样、对照和自有 reads 组装的非独立性均保留；结论仍限定为小部分 potentially somatic。对“多数并未证实”的报告充分，不必再新增一套验证实验。

### R3.4 完整输出与 CHM13

- **Reviewer 要求／回应应覆盖：** 量化 graph-coordinate 输出的代价，修正约 2%，并区分无法投射 GRCh38 与真正图特有。
- **评分：4.5 → 4/5，OK。** Rescue 开关、CHM13 比例及全部 graph-coordinate calls 计 FP 的 worst case 仍回答了输出负担。下调半分是因为新版“evaluation now covers the complete output”的概括比之前更强；这里成立的是**固定已知 truth 下的保守惩罚评估**，并没有测得参考外真实 recall。给定正文已写明 FP 假设，仍判为足够回应。

### R3.5 SBS96 图与 truth 比较

- **Reviewer 要求／回应应覆盖：** 展示 truth 与 callers，以及真正的 graph-coordinate SBS96，避免据谱相似性确认 somatic。
- **评分：3 → 4.5/5，OK。** 已逐张核对：两种 tumor 的普通比较图均归入 R3.5a；R3.5b 现在确实是六个 graph-coordinate 谱。表中 6,018 个有 context calls 与 6,025 减 7 对应，图上数量一致。旧图文错配已解决。

### R3.6 perfect alignment 与正确定位

- **Reviewer 要求／回应应覆盖：** 使用已知来源 reads 和 CHM13 控制，区分完整性、多单倍型与 ambiguity。
- **评分：4.5 → 4.5/5，OK。** 四参考模拟、MAPQ 过度自信、ONT 反例和 SV reads 不可评分的限制均清楚；没有把 shared evaluable regions 的结论外推到完整参考空间。

### R3.7 false negatives 不应轻率改判

- **Reviewer 要求／回应应覆盖：** 给独立证据，或撤回 benchmark 错误的原解释。
- **评分：3.5 → 4.5/5，OK。** 现在准确称 COLO829T 有自有长 reads 组装、但缺独立组装；v1.1 原始定义也已标明。原来“没有组装／已使用组装”的矛盾已消除，且保留了 normal assembly 本身的局限。

### R3.8 正常组织 mosaicism 的范围

- **Reviewer 要求／回应应覆盖：** 验证正常组织，或明确限定 cancer tumor-only。
- **评分：5 → 5/5，OK。** 范围和 VAF 限制明确，采用 reviewer 允许的收窄路径。无需补正常组织实验。

### R3.9 稳定性与统计不确定性

- **Reviewer 要求／回应应覆盖：** 多 seed、paired genomic-block CI、阈值敏感性及独立肿瘤结果的适当边界。
- **评分：3.5 → 3.5/5。** 单 seed 跨肿瘤 SNV 已明确降为探索性，初始化与划分波动也区分正确。但仍有可定位的模型对应不一致：R3.9a 中 HG008T PacBio v2.1 的 first-seed recall 为 **0.839**，GRCh38-only 为 **0.826**，差为 **0.013**；R3.9b 却为 **0.009**。后者恰与另一处 per-dataset-normalized 模型的 recall **0.835** 相符，提示可能混用了两套运行。
- **还需做：** 按 checkpoint、normalization、seed、阈值建立每张表的 run ID，对照 CI 是否确来自其宣称的主表模型；若 point estimate 是 bootstrap 均值而非原始差，也需说明。另补 chr2 confirmation 的总体 ΔF1 区间如何汇总 tumor/platform/seed，不能只写“95% CI”。不再要求已明确降格为探索性的 SNV 必须新增多 seed。[定位 per-seed 行](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:1040)；[定位 CI 行](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:1068)

### R3.10 低 VAF 和 INDEL 结论

- **Reviewer 要求／回应应覆盖：** 真实 somatic dilution，INDEL 的类型、长度、重复背景和平台分层，明确性能边界。
- **评分：3.5 → 4/5，OK。** Dilution、最终模型分层计数及跨肿瘤 chr2 confirmation 已补上；不再宣传 INDEL 优势，明确对 ClairS-TO 的不利结果。具体 CI 方法集中在 R3.9 处理。三 seed 平均 P/R 的调和平均不一定等于平均 F1，分层平均计数取整也可能相差 1，不应把这些正常现象误报为算术错误。

### R3.m1 参数与复现描述

- **Reviewer 要求／回应应覆盖：** 统一候选、训练、阈值、过滤及每项实验的例外，使方法可复现。
- **评分：3 → 3.5/5。** MAPQ 严格不等式、最终 INDEL 的新阈值、CASTLE chromosome splits 和后验 epoch 已讲清。新问题是最终 INDEL 又新增“五个 site features”，但该包内没有定义；同时多处称完整配置/manifest 已提供，当前目录只有正文和六张图。
- **还需做：** 定义五个特征、输入来源、计算/缩放及缺失值处理，注明是否使用 population prior；把实际 model/config/split manifest 随包或以固定链接提供。统一遗留旧口径，例如 Note R2.3 的 **61/702** 与当前 HG008T SNV truth **697** 的关系。这里要求补齐可查材料，不据缺附件推断实验未做。

### R3.m2 comparator 的配置与公平性

- **Reviewer 要求／回应应覆盖：** 版本、模型、资源、过滤、命令和训练暴露明确，区分原生与共同过滤评估。
- **评分：3.5 → 3.5/5。** 原生过滤叠加共同过滤不等于完全相同过滤，已说明正确；误用 ONT 模型已明确全表替换。ClairS-TO 自 v0.3.0 取消默认 INDEL 区域限制，也与[官方版本记录](https://github.com/HKU-BAL/ClairS-TO#whats-new)相符，不再列作缺口。
- **还需做：** 消除 PR 曲线的规则歧义：一处说除 RefCall/NoCall 外全纳入、忽略 native tags，下一句又说保留 technical quality filters。列清哪些 FILTER tags 保留、哪些解除、QUAL 如何扫描，以及 matched-precision recall 用哪条曲线；用实际命令/筛选表达式固定即可，无需重复工具背景。[定位过滤规则](/Users/jiaweishen/Downloads/Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md:749)

### R3.m3 图版本与组成

- **Reviewer 要求／回应应覆盖：** release、haplotypes、AF filter、参考路径、SV 和 donor/亲缘重叠。
- **评分：4.5 → 4.5/5，OK。** 所需信息完整，对 rare-allele sharing 不能正式排除亲缘的限制恰当。实际图表示与频率 proxy 的区分集中记在 R3.2，不重复增加任务。

### R3.m4 术语与语言

- **Reviewer 要求／回应应覆盖：** population resources/PoN、NPU、方法名称、图表和语法统一。
- **评分：4 → 4.5/5，OK。** 新版主回复的术语更一致，NPU 的淘汰原因明确，具体语言问题均有回应。

## 长度与提交前整理

**仍应压缩主回复，但保留这次新增的不利结果和证据表。** 新版总量约 2.68 万英文分词，含 reviewer 原话和补充表；即使剔除 reviewer 原话及文末补充部分，主回复仍约 1.45 万。相比上次，证据更完整，篇幅却没有缩下来。

- 将完整配置、分层计数、逐 seed 表固定为正式补充材料；主回复只留关键结果和对应位置。不要用“见完整表”代替实际可取得的表。
- 每个主要点保留直接回答、1–2 个决定性结果和适用限制；R1.4/R2.3/R3.7、R1.5/R2.1、R1.6/R2.4/R3.10 继续减少重复。短回复足够的 minor 不再扩写。
- 先完成上述 9 项的核对，再统一提交状态和摘要措辞。尤其不能让 cover letter 的“贡献已隔离／软件已可复现”等概括超过正文实际完成程度。

**本轮判断：值得保留的实质修订很多。下一步最有效的是固定一套可追溯的结果和标签规则，再缩短回复；不是继续增加大范围 benchmark。**
