# HG008T Illumina INDEL 为什么差（2026-09-30）

模型 `HG008_Illumina_INDEL_base`，chr1，somatic truth PASS 且在 BED 内的 INDEL 586 个。calls 用 predict 的 best-F1
threshold（`<run>/vcf_chr1/`，rtg `eval_bed`）。脚本只读，数字和原始输出都在本目录。

结论：INDEL 差主要不是模型的问题，是这个样本和 graph 上找 candidate 的方式决定的。ClairS-TO 在同一数据上 F1 0.123，
v2 是 0.128–0.150。

## 1. truth 的构成（`truth_indels.py` → `truth_indels.txt`，逐条见 `truth_indels.tsv`）

- 74% 在 ≥5 bp 的 homopolymer 里（431），STR 51，其他 104。
- 523 个（89%）的 allele 在 CoLoRSdb 里（repo PoN 的 truth 标记）。
- 254 个（43%）是 composite：GRCh38 上的长度变化和 HG008-N assembly 坐标（INFO `HG008Nv62SOMATICVARIANT`）上的不同，
  即 GRCh38 allele = germline 差异 + somatic 事件，例如 normal 上 +10A、GRCh38 上 +17A。

## 2. recall：candidate 阶段就丢了一半

- exact tensor（`tensor_representative`）236 个，另约 78 个只靠 partial tensor 对上；约 270 个（46%）没有任何 tensor。
- 157 个 truth 的 allele 就是 d9 graph 的一个 branch（与 `hprc-v1.1-mc-grch38.d9.vcf.gz` 同一条单倍型序列）：reads 走 branch，
  没有 edit，就没有 candidate。这 157 个里 exact tensor 为 0。
- `no_candidate` 282 个里 195 个是 composite；builder 门槛（≥3 reads、AF 0.08）挡掉 65 个。
- 对比全基因组：somatic INS 53%、DEL 37% 是 `no_candidate`，SNP 只有 4%（`tensors/truth_recall.json`）。

## 3. precision：FP 是真的 germline INDEL（`calls_indels.py` → `calls_indels.txt`）

- BED 内 731 个 FP：672 个（92%）是 germline truth，705 个在 homopolymer 里，AF 基本在 0.3–1。
- 92 个 TP 里 90 个也在 homopolymer 里，AF 基本 ≥ 0.3，分布和 FP 一样。HG008-T 是 cell line，somatic het 和 germline het
  的 VAF 都在 0.5 左右，tumor-only 只剩人群频率能分。
- 不在 GRCh38 上的 PASS 324 个：268 个是 `off_reference_no_truth_match`，31 个是 residual partial somatic。
- 训练集（chr2–22）somatic INDEL tensor 只有 4,856 个：exact 3,307、residual partial 1,518、allele partial 31。

## 4. PoN 子集（`pon_subsets.py` → `pon_subsets.txt`）

在 rtg 的 curve（PASS + LowQual）上去掉各 PoN 标记的 calls（近似：不重算 rtg）。没有一个子集能把 germline 和 somatic 分开：

| PoN | recall 上限 | best F1 |
|---|---|---|
| 无 | 1.00 | 0.150 |
| gnomAD | 0.88 | 0.152 |
| gnomAD + dbSNP + 1000G | 0.64 | 0.178 |
| 四个全用 | 0.09 | 0.113 |

## 可做的方向

1. 分层报告（homopolymer / composite），并训练 COLO829T 的 INDEL 模型，看是不是 HG008-T 特有。
2. INDEL 去掉 partial label 的 ablation；用 SNV 权重初始化、小模型或加 COLO829T 数据（best epoch 是 2，过拟合快）。
3. INDEL 用 gnomAD v3/v4 genomes（按 AF，allele normalize 后匹配），不用 CoLoRSdb。
4. off-reference 事件沿 graph 投影回 GRCh38（composite 的写法），curve 里最多多找回约 78 个 truth。
5. graph branch 当 candidate：recall 上限最多 +157，但 tumor-only 下先验几乎都是 germline，不推荐。
