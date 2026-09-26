# HG008-T：somatic 真值为什么没有候选（PacBio / ONT-UL / Illumina 三平台对比）

整理日期 2026-09-26。三个平台各做了一份逐个位点的分析：把每个"没有候选"的 somatic 真值附近的 reads 取出来，看 ALT reads 在图里怎么走、产生了哪些 edits、这些 edits 的 tensor 被打了什么标签。

## 数据

| 平台 | 测序 | tensors | 目录（`/scratch/jshen/data/pansoma_v2_tensors/`） |
|---|---|---|---|
| PacBio | HG008-T PacBio HiFi Revio 116× | format v6，v2 pipeline（job 364673/364850），有 supplement 轮 | `Liss_lab_PacBio_Revio_20240125/v6_tensors` |
| ONT | HG008-T ONT-UL R10.4.1 dorado 0.8.1 sup 54× | format v6，v3 pipeline（job 367607/367608），没有 supplement 轮 | `Liss_lab_Northeastern-ONT-UL-20241216/v3_tensors` |
| Illumina | HG008-T p23 BCM Illumina WGS | format v6，v3 pipeline（job 367751），深 node 下采样到 1 万条 read | `Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors` |

- 三个平台都只含常染色体，AF 阈值相同（SNV 0.06、INDEL 0.08，至少 3 条 reads）。
- 真值：
  - somatic：GIAB `HG008-T_somatic_smvar_benchmark_v0.2`，17,186 个 allele，全部 PASS；
  - germline：dipcall `HG008N` 的 `dip.vcf.gz` + `dip.bed`。
- 标签：三个平台现在都是 `truth-labels-v2`：
  - 1 和 2 要求 truth 是 PASS，不看 BED；
  - confident 区 = somatic BED ∩ germline BED。
- PacBio 和 ONT 的分析最初在旧标签下生成。标签换成 v2 后，两份表格重新生成并逐行核对过，报告里的数字没有变化。Illumina 的分析是在 v2 之后做的。

## 主要结论

1. **三个平台的主要原因相同：somatic allele 已经全部或部分在 HPRC 图里。** 候选只来自 reads 相对它所走节点的 edits。图里已经有这个 allele 时，ALT reads 走现成的分支或跳过边，要么没有 edits（I1），要么只剩一个写法不同的残余 edit（I2）。
   - INDEL 未命中里 I1 + I2 的占比：INS 91.4% / 90.3% / 83.2%，DEL 93.7% / 93.1% / 93.0%（依次为 PacBio / ONT / Illumina）。
   - 三个平台都没有候选的 truth：INS 1,523 个（30.8%），DEL 560 个（15.8%），SNV 241 个（2.8%）。其中大部分在三个平台被归为同一类，主要是 I1 / I2。
2. **SNV 三个平台都在 95% 左右**；漏掉的集中在图里本来就有 SNP 分叉（S3）或重复区（S4）的位置。
3. **各平台的特点：**
   - **PacBio**：DEL 的代表 allele 比例最高（59.3%）。v6 的 indel 左对齐让 36 个 truth 丢了 tensor、305 个从 filtered 变成没有候选：残余的 +A/+T 被推进了人群插入分支节点，没有 GRCh38 坐标。
   - **ONT**：ultra-long reads 能唯一比对到片段重复区（I5 只有 1 个，S1 只有 1 个）。噪声 allele 多：非代表 allele 多（DEL 12.6%），filtered 多，同一个 allele 被分散成多种 edit 写法。
   - **Illumina**：INDEL 召回最低（有 tensor 的比例 DEL 51.7%、INS 36.3%）。长插入和长重复里 150 bp 的 reads 看不到或跨不过（I4 INS 264，I5 INS 84）；片段重复区的 SNV 没有 MAPQ>10 的 reads（S1 55）。噪声 allele 很少（非代表 allele DEL 1.4%）。
   - 长读段对一部分 I1 truth 的召回，来自同聚物噪声 edits 碰巧左对齐到 truth 位置；Illumina 没有这种噪声，所以 I1 更多（见 Illumina 报告第 1 节第 2 条）。
4. **标签质量问题（三个平台都有，还没处理）：**
   - 残余 edit 的 tensor 被标成 germline：301 / 194 / 317 个不同候选。
   - site 的代表 allele 是 germline、另一个 allele 是 somatic truth，整个 tensor 被标成 2：254 / 232 / 40 个。
   - 209 个 somatic allele 和 germline 完全相同（多为正常样本杂合，可能是肿瘤 LOH 或两个 benchmark 冲突）；SNV 未命中里有 38 / 35 / 30 个属于这种。
5. **负样本的构成差别很大**：
   - PacBio 和 ONT 的 0 主要是 INDEL，大多是 1 bp 同聚物噪声；0 : 1 分别为 319 : 1 和 458 : 1。
   - Illumina 的 0 主要是 SNV：305 万个，0 : 1 为 369 : 1。它的 INDEL 负样本反而很少，0 : 1 只有 2 : 1。

## 目录

| 路径 | 内容 |
|---|---|
| `comparison_tables.md` | 三平台对比表：召回、同一个 truth 在三个平台的状态、三个平台都漏掉的类别、原因分类、标签 0/1/2/−1、标签质量问题（`compare_platforms.py` 生成） |
| `compare_platforms.py` | 生成上面的表 |
| `pacbio_v6/REPORT.md` | PacBio 报告：INDEL 用 v6 重新分析，SNV 沿用 v5（SNV tensors 两版相同）；含 v5 → v6 变化 |
| `ont/REPORT.md` | ONT 报告（SNV 和 INDEL 都分析），和 PacBio 对比 |
| `illumina/REPORT.md` | Illumina 报告（另一个对话用同样方法完成），和 ONT、PacBio 对比 |
| `*/indel_no_candidate.tsv`、`*/snv_no_candidate.tsv` | 每个未命中 truth 一行：类别、reads 数、REF/ALT reads、ALT 比例、附近 germline、残余 edits 及其 tensor 标签、ALT reads 的表示、序列背景 |
| `*/controls.tsv` | 有 tensor 的对照（同样方法分析） |
| `*/transitions.tsv` | 每个 truth 的状态和类别在两个版本或平台之间的变化：PacBio 是 v5 → v6，ONT 是 PacBio → ONT，Illumina 是 ONT / PacBio → Illumina |
| `*/residual_edits_labelled_germline.tsv` | 残余 edit 的 tensor 被标成 germline 的逐条列表 |
| `*/tables.md`、`*/summary.json` | 各报告里的全部表格和计数 |
| `*/examples_pick.tsv`、`*/show_examples.txt`、`*/show/` | 例子位点的 reads 逐条显示（局部序列、edits、图路径） |
| `*/reads.jsonl.gz` | 每个位点的原始结果（`reads.py` 输出） |
| `*/scripts/` | 分析脚本 |
| `ont/snv_detail.jsonl` 等 | SNV 细查结果；PacBio 的在 `pacbio_v6/snv_from_v5/` |
| `pacbio_v6/v5_reference/` | v5 的报告、v5 标签报告、v5 与 v6 的标签对比表（v5 数据已删除） |

## 类别说明

- **INDEL**
  - I1：allele 完全在图里，ALT reads 没有 edits；
  - I2：部分在图里，只剩残余 edit；
  - I3：走 GRCh38 节点，但 edits 写法和 truth 不同；
  - I4：没有或只有 1–2 条 ALT reads；
  - I5：没有 MAPQ>10 的 reads 跨过位点；
  - I6：reads 里的 indel 超过 50 bp；
  - I7：不确定；
  - I9：reads 写出了 truth 的 edit，但节点没建。
- **SNV**
  - S1：MAPQ 低或没有 reads；
  - S2：reads 里没有 ALT；
  - S3：ALT 是图里现成的 SNP allele；
  - S4a：重复区，SNV 成了别处的 edit；
  - S4b：重复区，被图路径吸收；
  - S4c：重复区，没有清楚的 ALT。
  - 规则见 `ont/scripts/report.py` 的 `snv_class`。

## 复现

- 这些分析是在 `tmp/somatic_miss_analysis_{v6,ont,illumina}/` 里跑的，脚本里的路径都指向那里和数据目录，这里存的是副本。
- 流程：`reads.py`（Slurm 分段，读 GAM、按位点解码）→ SNV 另跑 `snv_detail.py` / `snv_nodes.py` → `report.py` → `tables.py` → 例子用 `show_locus.py`。
- 解码用各个 run 冻结的源码。PacBio 在 `v6_run/source`（v2 pipeline 已从仓库删除）；ONT 和 Illumina 在各自的 `v3_run/source`。
- 每条 read 只取位点附近、两侧各 600 bp 的部分。PacBio（60 个位点）和 ONT（24 个位点）上都和整条 read 解码对比过，结果完全相同。
