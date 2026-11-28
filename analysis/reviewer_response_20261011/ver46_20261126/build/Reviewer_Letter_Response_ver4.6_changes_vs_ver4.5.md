# Reviewer_Letter_Response_ver4.6.docx：相对 ver4.5 的全部改动

记录日期：2026-11-26。基准：`Reviewer_Letter_Response_ver4.5.docx`（md5 713dbb63…，没有改动，仍在原处）。ver4.6 是在 ver4.5 上按你 11-26 的清单（T2、T3、T5、T7、T9、T10、T11、Figure 1.3）用最新结果改出来的。本文件里的表是从两个 Word 文件里直接读出来的，所以就是信里现在的内容。编号沿用 `Reviewer_Letter_Response_ver4.5_changes_vs_ver3.md` 的 T / F 编号。

ver4.6 相对 ver4.5 一共改了：9 张表（5 组）、1 张图、14 段文字（11 段图注 / 表注，3 段正文）。其余逐段、逐格相同；ver3 的 6 条批注和已有的修订标记都在；所有改动都是直接改字，没有新增修订标记。

用到的最新结果：

| 数据 | 位置（repo `analysis/`） | 内容 |
|---|---|---|
| 不带 site features 的 six-tumor INDEL 模型，HPRC v1.1 / v2.1 × 三个平台，seed 0 | `m6_nositef_graphs_20261118/results_indel.md` | 11-24 / 25 在 lab cluster 训练完 |
| 不带 site features 的 six-tumor SNV 模型，HPRC v1.1 / v2.1 | `m6_nositef_graphs_20261118/results.md` | RIS，call sets 11-26 已拉回本地 |
| six-tumor SNV 模型 seed 1、2（HPRC v2.1） | `m6_seeds_20261123/results.md` | RIS，11-26 拉回 |
| 不带 haplotype-count channel 的 six-tumor SNV 模型（HPRC v2.1，seed 0） | `m6_nopathcount_20261124/results.md` | RIS，11-26 拉回 |

这次新算的东西都在 `analysis/reviewer_response_20261011/ver46_20261126/`（`tables/`、`regions/`、`strata/`、`fpdecomp/`、`pathcount/`、`figs/`、`build/`），没有训练任何模型，没有改流程代码，没有 commit。

| 编号 | 表 / 图 | 结果 |
|---|---|---|
| T2 | Response Table 1.2 = 2.10 = 3.4 | HPRC v1.1、v2.1 两列换成不带 site features 的 INDEL 模型（seed 0） |
| T3 | Response Table 1.5 = 2.5 = 3.15 | SNV 各列和 INDEL 一列的 HPRC 一侧换成不带 site features 的模型，区间重算 |
| T5 | Response Table 2.6 = 3.8 | 没有改：用本地的 call sets 重算，54 格全部相同 |
| T7 | Response Table 2.8 | HPRC 行换成不带 site features 的 INDEL 模型；每个 reference 一行（去掉 best seed / mean of 3 seeds） |
| T9 | Response Table 3.2 | 4 行 COLO829T single-tumor SNV 换成 9 行 six-tumor SNV（seeds 0 / 1 / 2） |
| T10 | Response Table 3.9 | 最后一列换成不带 site features 的 HPRC v2.1 INDEL 模型 |
| T11 | Response Table 3.10 | 没有改：最新结果不影响这张表；和来源文件核对一致。需要你说明想改什么 |
| F2 | Response Figure 1.3 | 换成 six-tumor 模型有 / 无 haplotype-count channel 的 9 行图，图注和正文一句随之改 |

## T2. Response Table 1.2 = 2.10 = 3.4（chr1 INDEL，无群体过滤）

- 改了什么：HPRC v1.1、HPRC v2.1 两列（每张表 18 格，三张表内容相同）。ver4.5 里这两列是带 site features 的模型，而且 v2.1 是三个 seed 的均值；现在是不带 site features 的 six-tumor INDEL 模型，seed 0（`m6_nositef_graphs_20261118/multicancer6_<P>_INDEL_{d9,d46}_nositef`，过滤前，PASS = ppr8，ONT = pp）。GRCh38-only 一列和两个对照工具的列没有变；GRCh38-only 的 9 格用现在的 report 重算过，和信里相同。
- 数值由计数算出（不是 report 里四位小数再取整），所以有 4 格和 `results_indel.md` 差 0.001。
- 来源：`ver46_20261126/tables/make_tables46.py` → `tables46.md`。

**ver4.5：**

| chr1 INDELs, no population filter (P / R / F1) | Pansoma GRCh38-only | Pansoma HPRC v1.1 | Pansoma HPRC v2.1 | ClairS-TO | DeepSomatic-TO |
|---|---|---|---|---|---|
| HG008T Illumina | 0.016 / 0.445 / 0.031 | 0.412 / 0.128 / 0.195 | 0.322 / 0.137 / 0.179 | 0.004 / 0.101 / 0.007 | 0.017 / 0.206 / 0.032 |
| HG008T PacBio HiFi | 0.015 / 0.475 / 0.029 | 0.398 / 0.189 / 0.256 | 0.464 / 0.177 / 0.256 | 0.004 / 0.092 / 0.007 | 0.080 / 0.037 / 0.051 |
| HG008T ONT | 0.015 / 0.097 / 0.027 | 0.478 / 0.073 / 0.127 | 0.525 / 0.078 / 0.136 | 0.004 / 0.087 / 0.007 | 0.043 / 0.019 / 0.026 |
| COLO829T Illumina | 0.006 / 0.298 / 0.013 | 0.300 / 0.137 / 0.188 | 0.227 / 0.145 / 0.154 | 0.002 / 0.191 / 0.003 | 0.010 / 0.244 / 0.019 |
| COLO829T Fiber-seq | 0.015 / 0.198 / 0.027 | 0.195 / 0.115 / 0.144 | 0.311 / 0.117 / 0.169 | 0.002 / 0.160 / 0.004 | 0.124 / 0.092 / 0.105 |
| COLO829T ONT | 0.006 / 0.008 / 0.006 | 0.189 / 0.053 / 0.083 | 0.228 / 0.076 / 0.113 | 0.002 / 0.137 / 0.003 | 0.016 / 0.015 / 0.016 |
| HCC1395T Illumina | 0.015 / 0.744 / 0.029 | 0.466 / 0.729 / 0.569 | 0.417 / 0.744 / 0.498 | 0.006 / 0.699 / 0.012 | 0.022 / 0.835 / 0.042 |
| HCC1395T PacBio HiFi | 0.010 / 0.797 / 0.020 | 0.508 / 0.714 / 0.594 | 0.515 / 0.699 / 0.589 | 0.005 / 0.526 / 0.011 | 0.271 / 0.647 / 0.382 |
| HCC1395T ONT | 0.019 / 0.534 / 0.036 | 0.570 / 0.549 / 0.559 | 0.526 / 0.596 / 0.558 | 0.006 / 0.564 / 0.012 | 0.249 / 0.481 / 0.328 |

**ver4.6：**

| chr1 INDELs, no population filter (P / R / F1) | Pansoma GRCh38-only | Pansoma HPRC v1.1 | Pansoma HPRC v2.1 | ClairS-TO | DeepSomatic-TO |
|---|---|---|---|---|---|
| HG008T Illumina | 0.016 / 0.445 / 0.031 | 0.125 / 0.208 / 0.156 | 0.140 / 0.182 / 0.158 | 0.004 / 0.101 / 0.007 | 0.017 / 0.206 / 0.032 |
| HG008T PacBio HiFi | 0.015 / 0.475 / 0.029 | 0.111 / 0.278 / 0.159 | 0.139 / 0.239 / 0.176 | 0.004 / 0.092 / 0.007 | 0.080 / 0.037 / 0.051 |
| HG008T ONT | 0.015 / 0.097 / 0.027 | 0.101 / 0.082 / 0.090 | 0.121 / 0.083 / 0.099 | 0.004 / 0.087 / 0.007 | 0.043 / 0.019 / 0.026 |
| COLO829T Illumina | 0.006 / 0.298 / 0.013 | 0.082 / 0.153 / 0.106 | 0.064 / 0.160 / 0.091 | 0.002 / 0.191 / 0.003 | 0.010 / 0.244 / 0.019 |
| COLO829T Fiber-seq | 0.015 / 0.198 / 0.027 | 0.051 / 0.122 / 0.072 | 0.070 / 0.153 / 0.096 | 0.002 / 0.160 / 0.004 | 0.124 / 0.092 / 0.105 |
| COLO829T ONT | 0.006 / 0.008 / 0.006 | 0.039 / 0.046 / 0.042 | 0.026 / 0.038 / 0.031 | 0.002 / 0.137 / 0.003 | 0.016 / 0.015 / 0.016 |
| HCC1395T Illumina | 0.015 / 0.744 / 0.029 | 0.136 / 0.744 / 0.230 | 0.184 / 0.767 / 0.297 | 0.006 / 0.699 / 0.012 | 0.022 / 0.835 / 0.042 |
| HCC1395T PacBio HiFi | 0.010 / 0.797 / 0.020 | 0.117 / 0.729 / 0.201 | 0.212 / 0.714 / 0.326 | 0.005 / 0.526 / 0.011 | 0.271 / 0.647 / 0.382 |
| HCC1395T ONT | 0.019 / 0.534 / 0.036 | 0.186 / 0.534 / 0.276 | 0.201 / 0.571 / 0.297 | 0.006 / 0.564 / 0.012 | 0.249 / 0.481 / 0.328 |


新表的情况：

- HPRC 两个图的 F1 在 9 个数据集上都高于 GRCh38-only（18 of 18）；v2.1 高于 v1.1 的是 7 of 9（COLO829T Illumina、COLO829T ONT 不是）。
- 数值比 ver4.5 低很多：v2.1 的 F1 从 0.113–0.589 变成 0.031–0.326。原因是这些模型不再拿群体频率做输入，过滤前的假阳性多（precision 0.03–0.21）。
- 和对照工具比（过滤前）：v2.1 高于两个对照工具的是 6 of 9。不是的三个：COLO829T Fiber-seq（0.096，DeepSomatic-TO 0.105）、HCC1395T PacBio HiFi（0.326 对 0.382）、HCC1395T ONT（0.297 对 0.328）。ver4.5 的表里是 9 of 9。信里没有句子说这一点，所以没有要改的正文。

连带改的文字：

| 位置 | ver4.5 | ver4.6 |
|---|---|---|
| Table 1.2 / 2.10 / 3.4 表注 | "GRCh38-only and HPRC v1.1, one run; HPRC v2.1, mean of three seeds;" | "one run (seed 0) on each reference;" |
| 同上 | "The HPRC models take population site features as an input (gnomAD and CoLoRSdb allele frequencies and the population-filter tag of each candidate), so their values without the filter are not free of population information; the GRCh38-only models do not." | "None of these models takes population site features as an input."（和 Table 1.1 的写法相同） |
| Reviewer 2 第 6 条第 3 点末句 | "INDELs show the same pattern before filtering (Response Table 2.10); their HPRC models take population allele frequencies as input, so their unfiltered results retain population information." | "INDELs show the same pattern before filtering (Response Table 2.10)." |

## T3. Response Table 1.5 = 2.5 = 3.15（区域分层，过滤前）

- 改了什么：18 行的第 3–9 列（每张表 126 格，三张表内容相同）。SNV 各列的 HPRC 一侧换成不带 site features 的 six-tumor SNV 模型（`multicancer6_<P>_SNV_{d9,d46}_nositef`），paired bootstrap 区间重新算了（2,000 次、1-Mb block，和原来同一套 resample）；最后一列 INDEL F1 的 HPRC 一侧换成不带 site features 的 INDEL 模型。GRCh38-only 一侧（箭头左边的数）没有变。
- 核对：108 个 rtg 目录的合计都等于各自的 report.json；三种 reference 的 truth、BED、测试集版本相同；SNV 的点估计和 RIS 上算的 `m6_nositef_graphs_20261118/regions.md` 逐格相同（216 格，最大差 0.000）。
- 来源：`ver46_20261126/regions/region_nositef.py` → `table_before_filter.md`（过滤后的版本在 `table_after_filter.md`，没有放进信里）。

**ver4.5：**

| chr1 dataset | HPRC graph | SMaHT extreme PASS calls (precision) | SMaHT extreme ΔF1 | GIAB all difficult ΔF1 | GIAB low map or segdup ΔF1 | GIAB tandem repeats and homopolymers ΔF1 | GIAB outside difficult ΔF1 | INDEL F1 homopolymers ≥ 7 bp |
|---|---|---|---|---|---|---|---|---|
| HG008T Illumina | v1.1 | 3,659 (0.007) → 163 (0.215) | +0.261 [0.182, 0.345] | +0.461 [0.391, 0.533] | +0.352 [0.252, 0.470] | +0.327 [0.236, 0.426] | +0.656 [0.621, 0.686] | 0.036 → 0.131 |
| HG008T Illumina | v2.1 | 3,659 (0.007) → 158 (0.209) | +0.251 [0.168, 0.340] | +0.440 [0.364, 0.519] | +0.324 [0.225, 0.449] | +0.330 [0.232, 0.438] | +0.669 [0.632, 0.699] | 0.036 → 0.102 |
| HG008T PacBio HiFi | v1.1 | 10,638 (0.004) → 88 (0.375) | +0.361 [0.266, 0.450] | +0.551 [0.481, 0.615] | +0.473 [0.370, 0.574] | +0.416 [0.307, 0.516] | +0.674 [0.637, 0.707] | 0.037 → 0.214 |
| HG008T PacBio HiFi | v2.1 | 10,638 (0.004) → 120 (0.283) | +0.314 [0.218, 0.408] | +0.508 [0.432, 0.582] | +0.390 [0.285, 0.509] | +0.428 [0.320, 0.533] | +0.680 [0.643, 0.714] | 0.037 → 0.190 |
| HG008T ONT | v1.1 | 6,484 (0.006) → 127 (0.252) | +0.282 [0.177, 0.410] | +0.447 [0.348, 0.548] | +0.353 [0.225, 0.518] | +0.336 [0.234, 0.437] | +0.689 [0.651, 0.723] | 0.025 → 0.052 |
| HG008T ONT | v2.1 | 6,484 (0.006) → 109 (0.294) | +0.309 [0.212, 0.410] | +0.458 [0.374, 0.538] | +0.374 [0.262, 0.518] | +0.370 [0.265, 0.466] | +0.682 [0.642, 0.717] | 0.025 → 0.064 |
| COLO829T Illumina | v1.1 | 5,273 (0.026) → 459 (0.312) | +0.364 [0.259, 0.504] | +0.477 [0.383, 0.581] | +0.319 [0.217, 0.473] | +0.356 [0.239, 0.539] | +0.711 [0.694, 0.726] | 0.012 → 0.111 |
| COLO829T Illumina | v2.1 | 5,273 (0.026) → 399 (0.368) | +0.415 [0.306, 0.546] | +0.512 [0.425, 0.602] | +0.366 [0.258, 0.509] | +0.386 [0.271, 0.546] | +0.715 [0.700, 0.730] | 0.012 → 0.091 |
| COLO829T Fiber-seq | v1.1 | 10,303 (0.014) → 194 (0.727) | +0.632 [0.574, 0.684] | +0.656 [0.621, 0.688] | +0.638 [0.583, 0.689] | +0.585 [0.521, 0.641] | +0.735 [0.720, 0.748] | 0.025 → 0.056 |
| COLO829T Fiber-seq | v2.1 | 10,303 (0.014) → 220 (0.609) | +0.563 [0.478, 0.644] | +0.606 [0.552, 0.651] | +0.528 [0.438, 0.619] | +0.541 [0.467, 0.612] | +0.734 [0.720, 0.748] | 0.025 → 0.024 |
| COLO829T ONT | v1.1 | 5,334 (0.022) → 420 (0.319) | +0.367 [0.269, 0.488] | +0.463 [0.382, 0.551] | +0.324 [0.231, 0.454] | +0.361 [0.258, 0.501] | +0.708 [0.694, 0.723] | 0.000 → 0.049 |
| COLO829T ONT | v2.1 | 5,334 (0.022) → 302 (0.440) | +0.453 [0.319, 0.583] | +0.525 [0.414, 0.613] | +0.404 [0.263, 0.556] | +0.518 [0.458, 0.579] | +0.696 [0.679, 0.711] | 0.000 → 0.050 |
| HCC1395T Illumina | v1.1 | 4,780 (0.030) → 191 (0.743) | +0.732 [0.682, 0.776] | +0.767 [0.744, 0.790] | +0.743 [0.696, 0.788] | +0.737 [0.694, 0.777] | +0.811 [0.801, 0.821] | 0.011 → 0.242 |
| HCC1395T Illumina | v2.1 | 4,780 (0.030) → 181 (0.751) | +0.720 [0.664, 0.775] | +0.771 [0.746, 0.796] | +0.757 [0.702, 0.807] | +0.721 [0.675, 0.766] | +0.817 [0.807, 0.826] | 0.011 → 0.321 |
| HCC1395T PacBio HiFi | v1.1 | 7,522 (0.019) → 142 (0.838) | +0.728 [0.671, 0.782] | +0.728 [0.698, 0.757] | +0.725 [0.665, 0.785] | +0.692 [0.642, 0.741] | +0.795 [0.784, 0.806] | 0.009 → 0.300 |
| HCC1395T PacBio HiFi | v2.1 | 7,522 (0.019) → 151 (0.854) | +0.769 [0.713, 0.820] | +0.770 [0.742, 0.795] | +0.740 [0.677, 0.799] | +0.760 [0.718, 0.799] | +0.807 [0.796, 0.818] | 0.009 → 0.381 |
| HCC1395T ONT | v1.1 | 9,206 (0.014) → 142 (0.725) | +0.635 [0.563, 0.703] | +0.658 [0.620, 0.692] | +0.630 [0.551, 0.711] | +0.602 [0.546, 0.658] | +0.730 [0.717, 0.743] | 0.011 → 0.300 |
| HCC1395T ONT | v2.1 | 9,206 (0.014) → 501 (0.224) | +0.307 [0.166, 0.551] | +0.458 [0.307, 0.617] | +0.260 [0.130, 0.517] | +0.594 [0.530, 0.659] | +0.722 [0.708, 0.735] | 0.011 → 0.286 |

**ver4.6：**

| chr1 dataset | HPRC graph | SMaHT extreme PASS calls (precision) | SMaHT extreme ΔF1 | GIAB all difficult ΔF1 | GIAB low map or segdup ΔF1 | GIAB tandem repeats and homopolymers ΔF1 | GIAB outside difficult ΔF1 | INDEL F1 homopolymers ≥ 7 bp |
|---|---|---|---|---|---|---|---|---|
| HG008T Illumina | v1.1 | 3,659 (0.007) → 382 (0.073) | +0.104 [0.063, 0.149] | +0.181 [0.145, 0.219] | +0.146 [0.097, 0.202] | +0.154 [0.105, 0.205] | +0.224 [0.199, 0.251] | 0.036 → 0.168 |
| HG008T Illumina | v2.1 | 3,659 (0.007) → 256 (0.137) | +0.187 [0.131, 0.251] | +0.313 [0.260, 0.373] | +0.248 [0.183, 0.330] | +0.287 [0.207, 0.362] | +0.420 [0.383, 0.455] | 0.036 → 0.159 |
| HG008T PacBio HiFi | v1.1 | 10,638 (0.004) → 412 (0.085) | +0.131 [0.087, 0.175] | +0.207 [0.170, 0.245] | +0.177 [0.123, 0.236] | +0.176 [0.127, 0.227] | +0.224 [0.199, 0.250] | 0.037 → 0.169 |
| HG008T PacBio HiFi | v2.1 | 10,638 (0.004) → 249 (0.149) | +0.210 [0.141, 0.281] | +0.341 [0.282, 0.406] | +0.269 [0.189, 0.365] | +0.299 [0.224, 0.376] | +0.401 [0.365, 0.436] | 0.037 → 0.178 |
| HG008T ONT | v1.1 | 6,484 (0.006) → 322 (0.090) | +0.129 [0.081, 0.183] | +0.199 [0.161, 0.241] | +0.195 [0.132, 0.270] | +0.150 [0.103, 0.200] | +0.260 [0.228, 0.291] | 0.025 → 0.086 |
| HG008T ONT | v2.1 | 6,484 (0.006) → 142 (0.183) | +0.212 [0.121, 0.304] | +0.348 [0.283, 0.416] | +0.308 [0.207, 0.419] | +0.294 [0.202, 0.386] | +0.468 [0.429, 0.503] | 0.025 → 0.083 |
| COLO829T Illumina | v1.1 | 5,273 (0.026) → 784 (0.179) | +0.226 [0.161, 0.309] | +0.303 [0.244, 0.370] | +0.187 [0.126, 0.280] | +0.259 [0.178, 0.365] | +0.479 [0.453, 0.503] | 0.012 → 0.090 |
| COLO829T Illumina | v2.1 | 5,273 (0.026) → 384 (0.391) | +0.436 [0.367, 0.503] | +0.478 [0.434, 0.521] | +0.399 [0.321, 0.483] | +0.500 [0.442, 0.552] | +0.601 [0.585, 0.616] | 0.012 → 0.107 |
| COLO829T Fiber-seq | v1.1 | 10,303 (0.014) → 545 (0.279) | +0.363 [0.282, 0.455] | +0.401 [0.342, 0.460] | +0.322 [0.235, 0.435] | +0.386 [0.300, 0.466] | +0.509 [0.477, 0.537] | 0.025 → 0.049 |
| COLO829T Fiber-seq | v2.1 | 10,303 (0.014) → 328 (0.451) | +0.500 [0.414, 0.587] | +0.528 [0.477, 0.578] | +0.450 [0.356, 0.553] | +0.558 [0.500, 0.609] | +0.633 [0.617, 0.647] | 0.025 → 0.083 |
| COLO829T ONT | v1.1 | 5,334 (0.022) → 1,110 (0.135) | +0.180 [0.122, 0.266] | +0.241 [0.185, 0.315] | +0.155 [0.099, 0.244] | +0.269 [0.216, 0.330] | +0.421 [0.392, 0.446] | 0.000 → 0.016 |
| COLO829T ONT | v2.1 | 5,334 (0.022) → 624 (0.244) | +0.311 [0.186, 0.477] | +0.375 [0.261, 0.488] | +0.264 [0.148, 0.451] | +0.466 [0.414, 0.514] | +0.583 [0.565, 0.599] | 0.000 → 0.000 |
| HCC1395T Illumina | v1.1 | 4,780 (0.030) → 386 (0.345) | +0.422 [0.370, 0.474] | +0.486 [0.449, 0.524] | +0.422 [0.343, 0.501] | +0.461 [0.414, 0.507] | +0.562 [0.545, 0.580] | 0.011 → 0.112 |
| HCC1395T Illumina | v2.1 | 4,780 (0.030) → 271 (0.517) | +0.579 [0.533, 0.626] | +0.628 [0.600, 0.656] | +0.596 [0.534, 0.658] | +0.589 [0.546, 0.633] | +0.679 [0.664, 0.693] | 0.011 → 0.137 |
| HCC1395T PacBio HiFi | v1.1 | 7,522 (0.019) → 395 (0.344) | +0.445 [0.400, 0.491] | +0.465 [0.433, 0.500] | +0.415 [0.348, 0.480] | +0.448 [0.399, 0.497] | +0.546 [0.529, 0.563] | 0.009 → 0.055 |
| HCC1395T PacBio HiFi | v2.1 | 7,522 (0.019) → 262 (0.519) | +0.594 [0.545, 0.644] | +0.610 [0.577, 0.644] | +0.578 [0.514, 0.644] | +0.573 [0.525, 0.621] | +0.677 [0.664, 0.690] | 0.009 → 0.139 |
| HCC1395T ONT | v1.1 | 9,206 (0.014) → 557 (0.214) | +0.301 [0.234, 0.373] | +0.375 [0.330, 0.419] | +0.305 [0.233, 0.395] | +0.362 [0.316, 0.411] | +0.471 [0.454, 0.486] | 0.011 → 0.079 |
| HCC1395T ONT | v2.1 | 9,206 (0.014) → 1,345 (0.084) | +0.122 [0.063, 0.244] | +0.250 [0.158, 0.381] | +0.103 [0.049, 0.214] | +0.424 [0.360, 0.489] | +0.612 [0.598, 0.627] | 0.011 → 0.111 |


新表的情况（18 个 dataset–graph 对）：

| ΔF1 列 | ver4.5（带 site features） | ver4.6（不带） | 区间在 0 以上 |
|---|---|---|---|
| SMaHT extreme | +0.251 到 +0.769 | +0.104 到 +0.594 | 18 of 18 |
| GIAB all difficult | +0.440 到 +0.771 | +0.181 到 +0.628 | 18 of 18 |
| GIAB low map or segdup | +0.260 到 +0.757 | +0.103 到 +0.596 | 18 of 18 |
| GIAB tandem repeats and homopolymers | +0.327 到 +0.760 | +0.150 到 +0.589 | 18 of 18 |
| GIAB outside difficult | +0.656 到 +0.817 | +0.224 到 +0.679 | 18 of 18 |

- 增益比 ver4.5 小，但每一列 18 个区间都在 0 以上。
- 现在 v2.1 在 9 个数据集里的 8 个上每一列都高于 v1.1（HCC1395T ONT 除外）；ver4.5 的表里两个图大致相同。
- HCC1395T ONT v2.1 是最弱的一行：SMaHT extreme 里还剩 1,345 个 PASS calls（precision 0.084），v1.1 是 557（0.214）。
- INDEL 一列：17 of 18 上升（唯一不升的是 COLO829T ONT v2.1，0.000 → 0.000，两边都没有 TP）。HPRC 一侧的数值多数比 ver4.5 低（例如 HCC1395T Illumina 0.242 / 0.321 → 0.112 / 0.137）。这一列里 COLO829T 和 HCC1395T 的 truth 只有 75 和 26 个。

正文没有需要改的，逐句核对过：Reviewer 1 第 3 条 "in all eighteen dataset–graph comparisons"、"no significantly lower SNV F1 in any difficult stratum"、"HPRC v2.1 reduced false-positive calls and improved F1 across all datasets"、"higher on the HPRC graphs in 17 of the 18 comparisons"，Reviewer 2 第 5 条 "HPRC v2.1 reduces false-positive calls and improves precision in all three tumors on all three platforms"，在新表下都成立。

表注改了一处（三张表相同）："the HPRC models take population site features as an input (population allele frequencies and the population-filter tag of each candidate) and the GRCh38-only models do not." → "none of these models takes population site features as an input."

## T5. Response Table 2.6 = 3.8（SNV 假阳性按胚系 / 图分类）：没有改

这张表在 ver4.5 里已经是两侧都不带 site features 的模型（最后一列是 11-24 你贴过来的 RIS 数）。这次 RIS 的 call sets 已经在本地，用原脚本对着它们把整张表重算了一遍：6 列 × 9 行全部和信里相同（`ver46_20261126/fpdecomp/checks.txt`）。所以表、图注、正文都没有动。

你在清单里给 T5 写的是 "No INDEL no site feature model so far"。这张表本身是 SNV 的，不依赖 INDEL 模型；它的 INDEL 版本是 Table 3.9（T10，已改）。如果你的意思是 Reviewer 2 的部分也要一张 INDEL 的假阳性分类表（现在只有 Reviewer 3 有），数据已经有了（就是新的 Table 3.9），告诉我就加。

## T7. Response Table 2.8（HG008T INDEL，按邻近胚系变异数）

- 改了什么：表从 12 行数据变成 9 行：每个平台 GRCh38-only、HPRC v1.1、HPRC v2.1 各一行。去掉了 "HPRC v2.1, best seed" 和 "HPRC v2.1, mean of 3 seeds" 两种行；HPRC 的 6 行换成不带 site features 的 six-tumor INDEL 模型，seed 0（过滤后，PASS = ppr8，ONT = pp）。GRCh38-only 的 3 行没有变（重算后和信里相同）。
- 核对：每个 call set 的合计等于它的 report.json；truth 三组是 452 / 84 / 51。
- 来源：`ver46_20261126/strata/strata_indel_nositef.py` → `strata_indel_nositef.md`（用的是 E19 原来的分层方法）。

**ver4.5：**

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 139 / 34 / 17 | 0.423 / 0.105 / 0.105 | 0.226 / 0.048 / 0.039 | 0.294 / 0.066 / 0.057 | – |
|  | HPRC v1.1 | 85 / 13 / 7 | 0.448 / 0.235 / 0.222 | 0.153 / 0.048 / 0.039 | 0.228 / 0.079 / 0.067 | −0.067 / +0.014 / +0.010 |
|  | HPRC v2.1, best seed (seed 2) | 75 / 17 / 10 | 0.516 / 0.227 / 0.167 | 0.177 / 0.060 / 0.039 | 0.264 / 0.094 / 0.063 | −0.031 / +0.029 / +0.006 |
|  | HPRC v2.1, mean of 3 seeds | 78.3 / 23.3 / 12.7 | 0.472 / 0.207 / 0.119 | 0.155 / 0.071 / 0.033 | 0.233 / 0.106 / 0.051 | −0.062 / +0.040 / −0.006 |
| PacBio HiFi | GRCh38-only | 325 / 137 / 79 | 0.246 / 0.062 / 0.060 | 0.235 / 0.107 / 0.098 | 0.240 / 0.078 / 0.074 | – |
|  | HPRC v1.1 | 107 / 23 / 11 | 0.473 / 0.258 / 0.083 | 0.212 / 0.095 / 0.020 | 0.293 / 0.139 / 0.032 | +0.053 / +0.061 / −0.042 |
|  | HPRC v2.1, best seed (seed 2) | 85 / 9 / 6 | 0.525 / 0.500 / 0.143 | 0.208 / 0.107 / 0.020 | 0.298 / 0.176 / 0.034 | +0.058 / +0.098 / −0.040 |
|  | HPRC v2.1, mean of 3 seeds | 91.0 / 13.0 / 6.0 | 0.496 / 0.419 / 0.145 | 0.198 / 0.107 / 0.020 | 0.283 / 0.170 / 0.034 | +0.042 / +0.091 / −0.040 |
| ONT | GRCh38-only | 197 / 28 / 9 | 0.140 / 0.034 / 0.100 | 0.071 / 0.012 / 0.020 | 0.094 / 0.018 / 0.033 | – |
|  | HPRC v1.1 | 36 / 6 / 1 | 0.526 / 0.250 / 0.500 | 0.088 / 0.024 / 0.020 | 0.152 / 0.043 / 0.038 | +0.058 / +0.026 / +0.005 |
|  | HPRC v2.1, best seed (seed 0) | 34 / 7 / 2 | 0.564 / 0.222 / 0.333 | 0.097 / 0.024 / 0.020 | 0.166 / 0.043 / 0.037 | +0.072 / +0.025 / +0.004 |
|  | HPRC v2.1, mean of 3 seeds | 32.7 / 5.7 / 2.0 | 0.566 / 0.265 / 0.333 | 0.094 / 0.024 / 0.020 | 0.162 / 0.044 / 0.037 | +0.068 / +0.026 / +0.004 |

**ver4.6：**

| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |
|---|---|---|---|---|---|---|
| Illumina | GRCh38-only | 139 / 34 / 17 | 0.423 / 0.105 / 0.105 | 0.226 / 0.048 / 0.039 | 0.294 / 0.066 / 0.057 | – |
|  | HPRC v1.1 | 103 / 19 / 10 | 0.464 / 0.136 / 0.167 | 0.197 / 0.036 / 0.039 | 0.276 / 0.057 / 0.063 | −0.018 / −0.009 / +0.006 |
|  | HPRC v2.1 | 97 / 16 / 8 | 0.446 / 0.360 / 0.273 | 0.173 / 0.107 / 0.059 | 0.249 / 0.165 / 0.097 | −0.046 / +0.100 / +0.040 |
| PacBio HiFi | GRCh38-only | 325 / 137 / 79 | 0.246 / 0.062 / 0.060 | 0.235 / 0.107 / 0.098 | 0.240 / 0.078 / 0.074 | – |
|  | HPRC v1.1 | 139 / 21 / 7 | 0.444 / 0.323 / 0.222 | 0.246 / 0.119 / 0.039 | 0.316 / 0.174 / 0.067 | +0.076 / +0.096 / −0.007 |
|  | HPRC v2.1 | 124 / 11 / 5 | 0.434 / 0.421 / 0.286 | 0.210 / 0.095 / 0.039 | 0.283 / 0.155 / 0.069 | +0.043 / +0.077 / −0.005 |
| ONT | GRCh38-only | 197 / 28 / 9 | 0.140 / 0.034 / 0.100 | 0.071 / 0.012 / 0.020 | 0.094 / 0.018 / 0.033 | – |
|  | HPRC v1.1 | 70 / 13 / 1 | 0.381 / 0.133 / 0.000 | 0.095 / 0.024 / 0.000 | 0.152 / 0.040 / 0.000 | +0.058 / +0.023 / −0.033 |
|  | HPRC v2.1 | 103 / 6 / 0 | 0.270 / 0.000 / 1.000 | 0.084 / 0.000 / 0.039 | 0.128 / 0.000 / 0.075 | +0.034 / −0.018 / +0.043 |


新表的情况（6 个 HPRC 行）：

| 组 | ΔF1 范围 | ΔF1 > 0 | precision 高于 GRCh38-only | FP 减少 |
|---|---|---|---|---|
| 没有邻近胚系变异 | −0.046 到 +0.076 | 4 of 6 | 6 of 6 | 26–64% |
| 1 个 | −0.018 到 +0.100 | 4 of 6 | 5 of 6 | 44–92% |
| ≥ 2 个 | −0.033 到 +0.043 | 3 of 6 | 5 of 6 | 41–100% |

- 假阳性在 18 格里全部低于 GRCh38-only；ΔF1 有正有负，没有 "邻近胚系变异越多增益越大" 的规律。
- ONT 的 0.000 和 1.000 是很小的计数：v1.1 的 ≥2 组 0 TP、1 FP；v2.1 的 1 个组 0 TP、6 FP；v2.1 的 ≥2 组 2 TP、0 FP。≥2 组每个 call set 只找到 0–5 个 truth（共 51 个）。
- 正文 Reviewer 2 第 5 条第 4 点的 INDEL 一句（"For INDELs, exact-allele candidate availability is lower on the graphs where germline variants cluster, and the results are more variable (Response Table 2.8)."）在新表下仍然成立，没有动。

表注改了两处：

| ver4.5 | ver4.6 |
|---|---|
| "Six-tumor INDEL models after population filtering; GRCh38-only and HPRC v1.1, one training seed; HPRC v2.1, the best of three seeds (highest chr1 F1 on that platform) and the mean of the three seeds;" | "Six-tumor INDEL models after population filtering, one training seed (seed 0) per model;" |
| "The HPRC models take population site features as an input and the GRCh38-only models do not." | "None of these models takes population site features as an input." |

## T9. Response Table 3.2（per-seed chr1 结果）

- 改了什么：最后 4 行（COLO829T single-tumor SNV，HPRC v2.1 / v1.1）换成 9 行 six-tumor SNV 模型（HPRC v2.1，seeds 0 / 1 / 2，过滤后），9 个数据集各一行。前 18 行 INDEL 的数值没有变；它们的 "training set" 一格从 "6 tumors" / "4 tumors" 改成 "6 tumors, INDEL" / "4 tumors, INDEL"，和新的 "6 tumors, SNV" 对应。表从 22 行数据变成 27 行。
- seed 0 = `snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe`（信里所有 six-tumor SNV 结果用的模型），seed 1、2 = `m6_seeds_20261123`。数值由计数算出。
- 去掉 COLO829T single-tumor 的 4 行是我的判断：它们是当时没有 six-tumor seeds 时的替代，信里其他地方都是 six-tumor 模型。想保留的话可以加回去。
- 来源：`ver46_20261126/tables/make_tables46.py` → `tables46.md`。

**ver4.5（最后 4 行，已删）：**

| tumor | platform | training set | precision | recall | F1 |
|---|---|---|---|---|---|
| COLO829T | Illumina | COLO829T only, SNV, HPRC v2.1 | 0.802 / 0.665 / 0.828 (0.765 ± 0.088) | 0.791 / 0.799 / 0.805 (0.798 ± 0.007) | 0.797 / 0.725 / 0.816 (0.779 ± 0.048) |
| COLO829T | Illumina | COLO829T only, SNV, HPRC v1.1 | 0.749 / 0.648 / 0.677 (0.691 ± 0.052) | 0.800 / 0.794 / 0.773 (0.789 ± 0.014) | 0.774 / 0.714 / 0.722 (0.736 ± 0.032) |
| COLO829T | Fiber-seq | COLO829T only, SNV, HPRC v2.1 | 0.852 / 0.856 / 0.853 (0.854 ± 0.002) | 0.786 / 0.786 / 0.794 (0.789 ± 0.005) | 0.818 / 0.820 / 0.823 (0.820 ± 0.003) |
| COLO829T | Fiber-seq | COLO829T only, SNV, HPRC v1.1 | 0.851 / 0.857 / 0.872 (0.860 ± 0.011) | 0.788 / 0.784 / 0.782 (0.785 ± 0.003) | 0.819 / 0.819 / 0.825 (0.821 ± 0.003) |

**ver4.6（新增 9 行）：**

| tumor | platform | training set | precision | recall | F1 |
|---|---|---|---|---|---|
| HG008T | Illumina | 6 tumors, SNV | 0.540 / 0.576 / 0.567 (0.561 ± 0.018) | 0.776 / 0.783 / 0.782 (0.780 ± 0.004) | 0.637 / 0.664 / 0.657 (0.653 ± 0.014) |
| HG008T | PacBio HiFi | 6 tumors, SNV | 0.542 / 0.536 / 0.534 (0.537 ± 0.004) | 0.813 / 0.829 / 0.799 (0.814 ± 0.015) | 0.651 / 0.651 / 0.640 (0.647 ± 0.006) |
| HG008T | ONT | 6 tumors, SNV | 0.626 / 0.617 / 0.609 (0.617 ± 0.009) | 0.656 / 0.736 / 0.730 (0.707 ± 0.045) | 0.641 / 0.671 / 0.664 (0.659 ± 0.016) |
| COLO829T | Illumina | 6 tumors, SNV | 0.775 / 0.840 / 0.836 (0.817 ± 0.036) | 0.756 / 0.718 / 0.747 (0.740 ± 0.020) | 0.766 / 0.774 / 0.789 (0.776 ± 0.012) |
| COLO829T | Fiber-seq | 6 tumors, SNV | 0.864 / 0.845 / 0.840 (0.850 ± 0.013) | 0.749 / 0.770 / 0.772 (0.764 ± 0.013) | 0.803 / 0.806 / 0.805 (0.804 ± 0.002) |
| COLO829T | ONT | 6 tumors, SNV | 0.787 / 0.790 / 0.754 (0.777 ± 0.020) | 0.708 / 0.722 / 0.726 (0.719 ± 0.010) | 0.745 / 0.755 / 0.740 (0.747 ± 0.008) |
| HCC1395T | Illumina | 6 tumors, SNV | 0.863 / 0.874 / 0.870 (0.869 ± 0.006) | 0.905 / 0.897 / 0.891 (0.897 ± 0.007) | 0.883 / 0.885 / 0.880 (0.883 ± 0.003) |
| HCC1395T | PacBio HiFi | 6 tumors, SNV | 0.877 / 0.873 / 0.870 (0.874 ± 0.003) | 0.858 / 0.856 / 0.863 (0.859 ± 0.004) | 0.867 / 0.864 / 0.866 (0.866 ± 0.002) |
| HCC1395T | ONT | 6 tumors, SNV | 0.760 / 0.767 / 0.611 (0.713 ± 0.088) | 0.676 / 0.684 / 0.674 (0.678 ± 0.005) | 0.716 / 0.723 / 0.641 (0.693 ± 0.045) |


SNV 的 F1 在三个 seed 之间的 SD：0.002–0.045（中位数 0.008）。PacBio HiFi / Fiber-seq ≤ 0.006，Illumina 0.003–0.014，ONT 0.008–0.045；最大的是 HCC1395T ONT（0.716 / 0.723 / 0.641，seed 2 的假阳性是 1,474，另两个是 733 和 716）。

连带改的文字：

| 位置 | ver4.5 | ver4.6 |
|---|---|---|
| Table 3.2 表注 | "Per-seed chr1 results of the final INDEL models and of the single-tumor COLO829T SNV models." | "Per-seed chr1 results (seeds 0 / 1 / 2, mean ± s.d.; after population filtering) of the final INDEL models (six-tumor and four-tumor) and of the six-tumor SNV models on HPRC v2.1." |
| Reviewer 3 第 9 条 Training seeds | "… and the single-tumor COLO829T SNV models of the current pipeline (Illumina and Fiber-seq, HPRC v1.1 and v2.1) were trained with three seeds …" | "… and the six-tumor SNV models on HPRC v2.1 were trained with three seeds …" |
| 同上 | "was 0.003–0.048 (median 0.018) for these SNV models" | "was 0.002–0.045 (median 0.008) for these SNV models" |
| 同上 | "the largest spreads were COLO829T Illumina SNV on HPRC v2.1 (F1 0.797 / 0.725 / 0.816)," | "the largest spreads were HCC1395T ONT SNV (F1 0.716 / 0.723 / 0.641)," |
| Table 1.8 表注 | "SNV models: one run (seed 0) each;" | "SNV models: seed 0 (the six-tumor SNV models on HPRC v2.1 were also trained with seeds 1 and 2);" |

Table 1.8 那一处不在你的清单里：原句说 SNV 模型都只有一个 seed，和新的 Table 3.2 矛盾，所以一起改了。INDEL 的 SD（0.004–0.064，中位数 0.014）和另两个最大波动的例子没有变。

## T10. Response Table 3.9（INDEL 假阳性按胚系 / 图分类）

- 改了什么：只有最后一列（HPRC v2.1 FP before population filtering，9 格），换成不带 site features 的 HPRC v2.1 six-tumor INDEL 模型（seed 0）过滤前的 FP 数。前 5 列没有变。
- 为什么只有这一列：前 5 列只用 GRCh38-only 模型的假阳性、matched-normal 胚系 VCF 和 HPRC v2.1 测试集里有没有这个候选（tensor），不用 HPRC 模型的 call。重算后前 5 列 9 行都和信里相同。
- 来源：`ver46_20261126/fpdecomp/fp_decomp_nositef.py` → `table_INDEL_before_filter.md`、`checks.txt`。

下面只列第 1、2、6 列（中间三列没有变）：

**ver4.5：**

| Dataset | GRCh38-graph FP before population filtering | HPRC v2.1 FP before population filtering |
|---|---|---|
| HG008T Illumina | 16,016 | 109 |
| HG008T PacBio HiFi | 18,590 | 115 |
| HG008T ONT | 3,642 | 43 |
| COLO829T Illumina | 5,977 | 49 |
| COLO829T Fiber-seq | 1,764 | 27 |
| COLO829T ONT | 180 | 46 |
| HCC1395T Illumina | 6,594 | 92 |
| HCC1395T PacBio HiFi | 10,207 | 126 |
| HCC1395T ONT | 3,713 | 86 |

**ver4.6：**

| Dataset | GRCh38-graph FP before population filtering | HPRC v2.1 FP before population filtering |
|---|---|---|
| HG008T Illumina | 16,016 | 659 |
| HG008T PacBio HiFi | 18,590 | 868 |
| HG008T ONT | 3,642 | 357 |
| COLO829T Illumina | 5,977 | 308 |
| COLO829T Fiber-seq | 1,764 | 267 |
| COLO829T ONT | 180 | 188 |
| HCC1395T Illumina | 6,594 | 451 |
| HCC1395T PacBio HiFi | 10,207 | 354 |
| HCC1395T ONT | 3,713 | 302 |


需要你知道的一点：**COLO829T ONT 这一行，HPRC v2.1 的 FP（188）现在比 GRCh38-only（180）多。** GRCh38-only 的 ONT INDEL 模型在这个数据集上几乎不出 call（131 个 truth 只找到 1 个，recall 0.008；v2.1 是 5 个）。其余 8 行 v2.1 的 FP 是 GRCh38-only 的 3–15%。信里没有句子说 "INDEL FP 在每个数据集上都减少"（Reviewer 3 第 2 条引用这张表的句子只用了候选保留比例，那些数没有变），所以正文没有动。

表注改了一处："Six-tumor INDEL models (HPRC v2.1: seed 0); the GRCh38-graph models take no population site features as an input and the HPRC v2.1 models do." → "Six-tumor INDEL models without population site features as an input on both references (seed 0)." 表注里候选保留比例的一句（at most 0.2% on Illumina；15–29%；0.8–2.6%）不变。

## T11. Response Table 3.10（graph-coordinate calls 的验证）：没有改

你的清单里 T11 后面没有写要改什么。我核对的结果：

- 这张表是 release 的 six-tumor SNV 模型（HPRC v2.1，seed 0，带 site features）默认输出里的 935 个 graph-coordinate PASS calls。这次的新结果（no site features、seeds 1 / 2、no path-count、不带 site features 的 INDEL 模型）都不改变这 935 个 call。
- 信里的 Table 3.10 和 Table 3.11 与来源文件（`ver4_20261123/offref/table_3.10.md`、`table_3.11.md`）逐格相同。

如果你想要的是下面某一种，告诉我：(a) 把 935 个 call 换成不带 site features 的模型的 graph-coordinate calls 重新验证（需要重跑 assembly / reads 的验证）；(b) 加上 off-reference-only 模型的结果（11-25 的结论是它还不能区分 somatic 和 germline）；(c) 别的。

## F2. Response Figure 1.3（有 / 无 haplotype-count channel）

![Figure 1.3](Reviewer_Letter_Response_ver4.6_figure_changes/Figure_1.3_path_count_channel.png)

上半部分是 ver4.5（= ver3）的图，下半部分是 ver4.6 的图。单独的新图：`Reviewer_Letter_Response_ver4.6_figure_changes/Response_Figure_1.3_ver4.6.png`（和 `.pdf`）。

| | ver4.5（= ver3） | ver4.6 |
|---|---|---|
| 模型 | 10 月的 single-tumor 模型 | six-tumor 模型，当前流程，HPRC v2.1，seed 0 |
| 行 | 6 行：HG008T、COLO829T × 三个平台 | 9 行：加了 HCC1395T |
| 有 − 无，过滤后 F1 | −0.009 到 +0.311 | +0.004 到 +0.063，9 行全为正 |
| 有 − 无，过滤前 F1 | +0.220 到 +0.420 | +0.004 到 +0.068，9 行全为正 |
| 有 − 无，precision | −0.017 到 +0.505 | +0.001 到 +0.126，9 行全为正 |
| 有 − 无，recall | +0.018 到 +0.069 | −0.009 到 +0.022；5 行为正，COLO829T Fiber-seq 相差不到 0.001，3 行 ONT 为负 |
| b 的色标范围 | ±0.5 | ±0.15 |

- 每行的两个模型只差这个 channel：同一批训练集、同一 recipe、同一 seed、同样留出的染色体、同一 PASS 规则；"无" = `m6_nopathcount_20261124/multicancer6_<P>_SNV_newpipe_nopc`（`--drop-planes path_count`）。两边测试集相同（核对了测试集名和 chr1 的 prediction 数）。
- 版式、字体、颜色和 ver3 的图相同；没有加 seed 0 / 1 / 2 的范围线（你 11-24 的意见）。图比原来高（9 行），Word 里宽度不变，高度 2.02 → 2.64 英寸。
- 差别最大的是 COLO829T Illumina（过滤后 FP 536 → 871）和 HCC1395T ONT（733 → 1,353）；HCC1395T Illumina / PacBio HiFi 和 COLO829T ONT 的 F1 差只有 0.004–0.008。和 seed 之间的波动相比：带 channel 的三个 seed 的 F1 在 9 个数据集里有 6 个全部高于不带 channel 的模型，另 3 个是 2 of 3（HG008T PacBio HiFi、COLO829T ONT、HCC1395T ONT；见 `m6_nopathcount_20261124/results.md`）。这一点没有写进信里。
- 来源：`ver46_20261126/pathcount/fig13_six_tumor.py` → `fig_R1.3_six_tumor.png`、`table.md`（9 行的全部数值）。

连带改的文字：

| 位置 | ver4.5 | ver4.6 |
|---|---|---|
| Reviewer 1 第 1 条第 2 点 | "we retrained the SNV models with and without this channel using the same dataset and pipeline. Removing channel 6 reduced SNV F1 before population filtering in all twelve comparisons and after filtering in eleven of twelve; the only exception was HG008T ONT on v2.1. (Response Figure 1.3)" | "we retrained the six-tumor SNV models with and without this channel using the same dataset and pipeline. Removing the channel reduced SNV F1 in all nine comparisons, before population filtering (by 0.004–0.068) and after it (by 0.004–0.063), with lower precision in all nine (Response Figure 1.3)." |
| Figure 1.3 图注 | "The evaluation was performed separately using single tumor models trained on COLO829T and HG008T, with chromosome 1 SNVs restricted to benchmark regions. For each tumor dataset, we trained two models …" | "The evaluation uses the six-tumor models on HPRC v2.1 (HG008T, COLO829T and HCC1395T), with chromosome 1 SNVs restricted to benchmark regions. For each sequencing platform, we trained two models …" |

这样也去掉了 ver4.5 里记录过的两处不一致：正文的 "all twelve … eleven of twelve" 对 6 行的图，和 "channel 6" 对 Figure 1.1 图注里的第 3 个 channel。

另外，11-24 排队的 HCC1395T single-tumor 有 / 无模型（当前流程，11-25 跑完：Illumina 0.890 / 0.881，PacBio HiFi 0.869 / 0.862，ONT 0.595 / 0.649）没有用在这张图里：six-tumor 的对比一次覆盖了九个数据集，而且和信里其他表是同一批模型。

## 表格对应的柱状图

你的 `fig_*/` 文件夹里每张表有一张柱状图。这次改了数据的五张重画成 `_ver4.6`（原有文件都没有动）：

| 表 | 新图（`.png` + `.pdf`） | 和 `_ver4.5` 的差别 |
|---|---|---|
| T2，Table 1.2 = 2.10 = 3.4 | `fig_R1.1_SNV_noPoN/fig_R1.1_INDEL_noPoN_{f1,precision,recall}_ver4.6` | 只换数据 |
| T3，Table 1.5 = 2.5 = 3.15 | `fig_R1.5_region_dF1/fig_R1.5_region_dF1_ver4.6` | 只换数据 |
| T7，Table 2.8 | `fig_R2.8_HG008T_INDEL_nearby_germline/fig_R2.8_HG008T_INDEL_nearby_germline_ver4.6` | 每组 3 根柱（GRCh38-only、v1.1、v2.1），没有 best seed / mean 两根；precision 轴到 1.0（ONT v2.1 的 ≥2 组是 1.000，2 TP、0 FP） |
| T9，Table 3.2 | `fig_R3.2_per_seed/fig_R3.2_per_seed_ver4.6` | 3 列（三个肿瘤），每个平台 3 根柱：INDEL six-tumor、INDEL four-tumor、SNV six-tumor；原来第 4 列（COLO829T-only SNV）没有了 |
| T10，Table 3.9 | `fig_R3.9_INDEL_FP_classes/fig_R3.9_INDEL_FP_classes_ver4.6` | c 的图例是 "HPRC v2.1, no site features"；COLO829T ONT 的 v2.1 柱（188）比 GRCh38-only（180）长 |

原来的 ggplot 脚本不在集群上，集群上也没有 R，所以这五张是按 `_ver4.5` 的版式用 matplotlib 重画的（脚本在 `ver46_20261126/figs/`）。

## 没有改、但和这次的数据有关的地方

- **Reviewer 2 第 6 条第 2 点**（"… the six-tumor models on HPRC v1.1 and v2.1, which take population site features as an input, reject 97.8–98.9% of these SNV candidates …"）：这一句说的是带 site features 的 release 模型，句子本身写明了，所以没有动。不带 site features 的对应数字已经有了（`m6_nositef_graphs_20261118/r26p2_nositef.md`：模型拒掉 85.7–96.4%，过滤再去掉被叫成 somatic 的 67.2–87.4%，最后剩 1.1–2.3%，过滤后的假阳性里 89.5–96.5% 是胚系等位基因）。要不要换，你定。
- **Table 3.1（过滤后）** 仍然只有带 site features 的 HPRC 行（你 11-24 的决定），所以 Table 3.3 / 3.4（过滤前，不带）和 Table 3.1（过滤后，带）的 HPRC 结果不是同一批模型。两张表的图注各自写明了。不带 site features 的 INDEL 模型过滤后的结果（`results_indel.md` 第 2 节）不在信里：chr1 F1 在 v1.1 和 v2.1 上各有 8 of 9 高于 GRCh38-only（HG008T Illumina 不是：0.231 / 0.226 对 0.244）。
- **Reviewer 3 第 10 条** "four of the six … (mean INDEL F1 of three seeds)" 说的是 release 的 INDEL 模型（过滤后），不受这次改动影响。
- `Reviewer_Letter_Response_ver4.5_changes_vs_ver3.md` 第 10 节里其余的遗留问题（旧编号、Table 1.8 "all other models use every example"、候选层面的图表是 AF 0.08 等）都没有动。
