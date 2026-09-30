# HG008T Illumina chr1：靠 partial 匹配拿到 label 1 的 INDEL tensor（2026-09-30）

tensor 集 `/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors/INDEL`，label `truth-labels 197b5d25bbf4`
（规则见 `analysis/truth_labels.md`：R4 `allele_partial_somatic_truth`，R8 `residual_partial_somatic_truth`）。
chr1 上这类 tensor 共 128 个，每个一张图：

| 目录 | reason | n |
|---|---|---:|
| `png/residual_on_ref/` | residual_partial_somatic_truth，GRCh38 节点 | 66 |
| `png/residual_off_ref/` | residual_partial_somatic_truth，分支节点（off-reference，位置用 anchor） | 61 |
| `png/allele_partial/` | allele_partial_somatic_truth | 1 |

图用 `scripts/visualize_tensor.py` 画（8 个 channel，和其它 tensor 图一样），存成 256 色 PNG。标题前三行是本目录加的：
1. `candidate_id` | label 和 reason | overlap（allele 重叠或单倍型重叠，> 0.45 才算 partial）
2. 对上的 truth（VCF 写法和 GT）| truth 在 GRCh38 上的重复 context（`HP<n>` = n 个同一碱基；`STR <unit>x<copies>`）|
   tensor 的位置（GRCh38 节点给坐标，分支节点给 anchor 区间）
3. 两个 INDEL 模型的 chr1 `p_somatic` 和各自的阈值 t，过 t 的写 `called`：`HG008_Illumina_INDEL_base_b1024`（partial 当 1
   训练）、`HG008_Illumina_INDEL_nopartial_b1024`（partial 用 `--ignore-reasons` 去掉）

后三行是 visualizer 原有的：格式、A1 allele 和 AF、site coverage。

## 汇总

| 组 | n | overlap = 0.5 | 0.5–0.99 | 1.0 | HP ≥ 6 | AF 中位数 | A1 reads 中位数 | base 过 t | nopartial 过 t |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| residual_partial_somatic_truth，GRCh38 节点 | 66 | 26 | 33 | 7 | 51 | 0.17 | 24 | 3 | 2 |
| residual_partial_somatic_truth，分支节点 | 61 | 21 | 31 | 9 | 45 | 0.33 | 20 | 19 | 10 |
| allele_partial_somatic_truth | 1 | 1 | 0 | 0 | 1 | 0.14 | 28 | 0 | 0 |

- 128 个里 97 个的 truth 在 ≥ 6 bp 的 homopolymer 里；overlap 正好 0.5（规则的下限附近）的 48 个，1.0 的 16 个。
- GRCh38 节点上的 66 个 AF 中位数只有 0.17，把它们当 1 训练的 base 模型也只把 3 个打过 t；分支节点上的 61 个 AF 高一些（0.33），
  base 过 t 19 个。
- 例子：`png/residual_on_ref/001_170412_0_DEL_T_.png`：truth `chr1:5175534 CTTTTT>C`（删 5 个 T，HP26），A1 是删 1 个 T，
  162 条 reads 里 27 条（AF 0.167），overlap 0.5 来自单倍型重叠。
  `png/residual_off_ref/005_2752798_0_DEL_T_.png`：truth `chr1:18232209 C>C+12T`（HP8），reads 走一个 T 更长的分支节点
  （site 处 path count 很低），在上面再删 1 个 T，114 条里 66 条（AF 0.579），overlap 0.917。

逐个 tensor 的表（带图的链接）：`tables.md`；全部字段（shard、行号、truth、context、AF、两个模型的分数）：`index.tsv`。

## 重画

```
/wanglab/jshen/anaconda3/envs/polymarket-btc-5m-bot/bin/python analysis/indel_partial_tensors_20260930/render.py
```

（在仓库根目录跑；base env 的 matplotlib 在 NumPy 2 下坏了，所以用这个 env。约 8 分钟。）
