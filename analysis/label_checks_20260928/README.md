# truth-labels-v6 标签检查（2026-09-28）

v6 规则（重叠阈值 0.45、分支节点没对上 truth 的标 −1、Illumina SNV AF < 0.07 标 −1）已于 2026-09-28 给 HG008 和 COLO829T 各三组数据重打。本目录是之后做的检查：数字、脚本和原始输出都在这里。

## 1. v5 → v6 各组的标签变化

`v5_v6_counts.py`（v5 取自各样本目录的 `labels_backup_v3_tensors_truth-labels-v5_*`）。主要变化：

| 数据 | 类型 | 1 somatic | 2 germline | 0 |
|---|---|---|---|---|
| HG008 PacBio | INDEL | 6,181 → 7,050 | 156,180 → 191,411 | 1,654,001 → 1,518,564 |
| HG008 Illumina | SNV | 9,703 → 11,089 | 316,795 → 315,438 | 3,988,419 → 2,635,565 |
| COLO829T Illumina | SNV | 38,965 → 38,781 | 338,743 → 337,496 | 1,760,336 → 1,234,838 |
| COLO829T ONT | INDEL | 1,393 → 2,400 | 195,066 → 253,853 | 4,131,894 → 3,822,575 |

- AF 0.07 门槛去掉的 SNV 里，A1 本来是完全匹配 somatic truth 的：HG008 Illumina 2 个，**COLO829T Illumina 389 个**（COLO829T 低 AF 的 somatic 更多）。
- 各门槛下 somatic / germline / 0 保留多少，见 `af_offref.txt`；按链过滤的效果见 `strand_filter.txt`（v5 标签下算的）。

## 2. Illumina 的"部分匹配 somatic" SNV（`partial_snv_check.txt`）

HG008 Illumina 有 2,812 个 SNV tensor 靠部分匹配 somatic truth 标成 1。其中 96% 对上的是 INDEL truth（DEL 1,073、INS 1,624），只有 115 个对上 SNV truth。

- ALT reads 只在一条链上的比例：部分匹配的为 21–46%，完全匹配的 somatic SNV 只有 0.4%；PacBio 的部分匹配只有 3–6%。
- 抽 12 个逐条 read 回放（`replay_examples.txt`）后看，偏链多半**不代表是假信号**。短读段跨过长插入时，往往只有一个方向的 reads 能覆盖到，reads 实际上带着 truth 的 indel，只是比对器把它写成了一个 SNV。判断依据是每条 read 到 GRCh38 的编辑数比到 truth 单倍型的多出差不多一个事件长度，例如 chr8 那个 +12A，每条 read 的 `d_ref − d_truth` 是 11–12。

## 3. 发现的问题：单倍型重叠在"没有证据"时也会给出 0.5

`haplotype_overlap` 对每条 A1 read 算 `d_ref`、`d_truth`，先各减去 REF reads 编辑数的中位数 `b`（最低到 0），再算

    shared = (d_ref + event − d_truth) / 2,   overlap = shared / max(d_ref, event)

**当 `b ≥ d_ref` 且 `b ≥ d_truth` 时，两个距离都变成 0，overlap 正好等于 0.5，不管这条 read 带不带 truth 的变化。**

例子：`41365234:104:SNP:C>A`（HG008 Illumina，chr5:2461228，truth CA>C）。每条 A1 read 的 `(d_ref, d_truth)` 是 (1,1)、(2,2)、(3,3)……，也就是到 REF 和到 truth 一样远，reads 不带这个 deletion；b = 2，结果算出 overlap 0.5 > 0.45，被标成 1。

阈值是 0.6 时 0.5 过不了，这个问题被挡住了；降到 0.45 以后就放进来了。`partial_evidence.py` 对所有部分匹配 somatic 的 tensor 重算了原始证据（A1 reads 的 `d_ref − d_truth` 中位数，不减背景）。"没有证据"指这个中位数 ≤ 0，这类 tensor 几乎全部在 0.45–0.6 这一档：

| | INDEL | SNV | 合计 | 占该组 somatic 标签 |
|---|---:|---:|---:|---:|
| HG008 Illumina | 77 | 349 | 426 | 2.6% |
| HG008 PacBio | 221 | 92 | 313 | 2.0% |
| HG008 ONT | 477 | 162 | 639 | 4.0% |
| COLO829T Illumina | 17 | 42 | 59 | 0.1% |
| COLO829T ONT | 655 | 121 | 776 | 1.9%（INDEL 1 的 27%） |
| COLO829T fiberseq | 199 | 50 | 249 | 0.6% |

另外还有一批"弱证据"的（中位数大于 0 但小于事件长度的一半），数量见 `partial_evidence.txt`。

**建议的修复**（需要用户决定）：单倍型重叠只认"比 GRCh38 更接近 truth"的 read。每条 read 若 `d_truth ≥ d_ref`，重叠记 0；或者改成用原始差值 `(d_ref − d_truth) / event` 衡量，不再对两个距离做对称的背景扣除。修复后上表这些 tensor 会从 1 变成 0（在分支节点上的则变成 −1）。

## 文件

| 文件 | 内容 |
|---|---|
| `v5_v6_counts.py` | v5 → v6 每组各 reason 的计数 |
| `af_offref.py` / `af_offref.txt` | SNV AF 门槛对 1/2/0 的影响；分支节点 tensor 的标签分布 |
| `strand_filter.py` / `strand_filter.txt` | Illumina SNV：按链过滤 vs 提高 AF 对 somatic 和 0 的影响（v5 标签） |
| `partial_snv_check.py` / `partial_snv_check.txt` | 部分匹配 somatic 的 SNV：对上的 truth 类型、重叠档、偏链、AF；抽样例子 |
| `replay_examples.py` / `replay_examples.txt` / `partial_snv_examples.json` | 12 个例子逐条 read 的 `d_ref`、`d_truth`、`b` |
| `partial_evidence.py` / `partial_evidence.txt` | 六组所有部分匹配 somatic 的 tensor：原始 read 证据分档 |
