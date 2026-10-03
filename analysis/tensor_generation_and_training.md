# 从 GAM 到训练好的模型：tensor 生成与 pansoma_net_v2 训练

在这台服务器上从一个样本的 GAM 生成训练 tensor，再用它们训练和测试模型的完整步骤。命令都可以直接运行。
代码和数据的术语保留英文。现状部分是 2026-09-28 的状态。

- tensor pipeline：`indexed_gam_pipeline_v4/`（所有命令和选项见它的 `README.md`，标签规则见
  `tensor_postprocessing/README.md` 和 `analysis/truth_labels.md`）
- 模型：`machine_learning/pansoma_net_v2/`（`README.md`）

---

## 一、项目概览

### 1. 组成

| 目录 | 作用 |
|---|---|
| `indexed_gam_pipeline_v4/` | 把排好序、建好索引的 GAM 变成 tensor，然后按染色体合并、打标签 |
| `machine_learning/pansoma_net_v2/` | 直接在合并好的 tensor 上训练和预测 |
| `machine_learning/pansoma_net/` | 旧模型，只能读旧的 `.dat/.idx` tensor，新数据不用它 |
| `scripts/visualize_tensor.py` | 把一个 tensor 画成 8 个 panel 的 PNG |

**tensor 的格式**：

- 一个 **site**（同一个 node、同一起点、同一类型 SNV 或 INDEL）对应一个 tensor，大小 `(8, 200, 101)`，int8。
- 8 个 channel：read base、base quality、site allele、MAPQ、alignment operation、graph base、path count、strand。
- 200 行是 reads，101 列是以 site 为中心的窗口。

**标签**：1 = somatic，2 = germline，0 = non（不是真实变异），−1 = ignore（训练时不用）。

### 2. 现有数据（6 套 HPRC v1.1 d9 tensor，都可以直接用来训练）

都在 `/scratch/jshen/data/pansoma_v2_tensors/<样本>/tensors/`，样本目录是 `HG008T_PacBio`、`HG008T_ONT`、
`HG008T_Illumina`、`COLO829T_Illumina`、`COLO829T_fiberseq`、`COLO829T_ONT`。每个样本目录的 `README.txt` 写了输入、
参数、tensor 和标签数量、relabel 命令；数据根目录的 `README.txt` 说明整体布局（还有 linear / HPRC v2.1 d46 graph 和
P25/P50/P75 mixture 的 set）。数量和资源见 `indexed_gam_pipeline_v4/README.md` 第 5 节。

- **参数**：`--haplotypes 90`（channel 6 的 storage 是 `int8-count-haplotypes100`），build AF 的 SNV 和 INDEL 都是
  0.08，标签不设 AF floor（floor 不超过 build AF 时没有作用）。
- **标签都是最新规则**：12 个 `labels.manifest.json`（6 套 × SNV/INDEL）的 `rules_sha256` 都是 `5430533aa100`，和
  repo 里的 `truth_labels.py` 一致。
- 之前的一版（channel 6 是 `int8-count-linear100-log2`，SNV AF 0.06，两套 Illumina 的 SNV label floor 0.07）在
  `pansoma_v2_tensors/backup_ch6_linear100_20261001/`，只有之前训练的 checkpoint 读它。
- **检查一套数据的标签是否最新**：比较 `<set>/SNV/labels.manifest.json` 的 `rules_sha256` 和
  `sha256sum indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py`。
- relabel 正在进行时不要用那套数据训练：先确认对应的 relabel job 已经 COMPLETED。

### 3. 待决定的问题

1. **训练的 job 脚本不在 repo 里。** `run.sh` 和 `build_index.sh` 放在 `/scratch/jshen/data/pansoma_net_v2_runs/jobs/`，
   没有进 git。`run.sh` 的 `FP32CHECK=1` 选项还调用了 `tmp/` 下的一个脚本。可以像 v4 那样挪到
   `machine_learning/pansoma_net_v2/jobs/`。
2. **新样本没有现成的 v4 job 脚本。** `tools/jobs/` 只有 `graph_prep.sh`、`gam_sort.sh`、`relabel.sh`，没有
   discovery → prepare → 提交这一串。重建和 linear / d46 用的脚本在 repo 的 `tmp/`（不进 git）：
   `tmp/rebuild_ch6_20261001/<样本>/`、`tmp/linear_hg38_20261001/`、`tmp/hprc_v21_d46_tensors_20261003/`。下面第二部分
   是完整命令，可以整理成 `indexed_gam_pipeline_v4/tools/jobs/` 下的脚本。
3. **未被 git 跟踪的遗留目录。** `src/`（旧的编译产物 `.so`）、`build/`、根目录的 `__pycache__/` 已被 git 忽略，可以删。
4. **现有的 checkpoint 要重新训练。** 它们都是用 backup 里的旧 channel 6 storage 训练的，不能预测现在的 tensors
   （第三部分第 6 步最后一段）。

---

## 二、生成训练 tensors（从一个新样本的 GAM 开始，CPU，`general` partition）

### 第 0 步：前提和一次性准备

```bash
cd /scratch/jshen/Github/Pansoma
PY=/wanglab/jshen/anaconda3/bin/python
P=indexed_gam_pipeline_v4
G=/scratch/jshen/data/pansoma_v2_tensors/graph_index          # HPRC v1.1 d9 的 graph 文件，已经做好
GAM=/path/to/sample.sorted.gam                                 # 必须排过序，并且旁边有 sample.sorted.gam.gai
S=/scratch/jshen/data/pansoma_v2_tensors/MY_SAMPLE; mkdir -p $S

$PY -m $P.native compile      # 编译 C++ record decoder：比 Python 快约 28 倍，输出完全相同
$PY -m $P.native check        # 确认 builder 会用它；不能用时会说明原因
```

- **GAM 从哪来**：还没有 GAM 时，先按根目录 `README.md` 第 1 节用 vg giraffe 比对，再用 `vg gamsort -i` 排序并生成
  GAI。在这台机器上先 `source scripts/use_vg.sh`。
- **graph 文件**：`$G` 下已有 graph index（SQLite）、reference path 和 chr 表，同一个 graph 不用重做。换 graph 时用
  `indexed_gam_pipeline_v4/tools/jobs/graph_prep.sh`。
- **decoder 编译**：每个 checkout 做一次。prepare 会把编译好的 `.so` 冻结进 run，所以要先 compile 再 prepare。

### 第 1 步：discovery（选 target nodes）

**做什么**：扫一遍整个 GAM。一个 node 在 MAPQ > 5 的 mappings 里，如果超过 5% 在 indel left-normalization 之后带有
edits（和参考不一致的碱基或 indel），就成为 **target node**。只有 target node 上会生成 tensor。

```bash
sbatch -p general -c 48 --mem=20G -t 12:00:00 -J discover -o $S/discover-%j.out --wrap \
  "cd /scratch/jshen/Github/Pansoma && $PY -m $P.run discover --gam $GAM --output $S/discovery \
   --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite --processes 48"
```

- **耗时和内存**（实测）：PacBio 18.6 分钟 / 10.5 GB；ONT-UL 31 分钟 / 16 GB；Illumina 2 小时 25 分钟 / 15 GB；
  COLO829T 三个平台同时跑时各约 3.4 小时。
- **输出**：
  - `target_nodes.txt`：排好序的 node ID。
  - `node_stats.json`：每个 node 的 perfect 和 not_perfect mapping 数，prepare 用它估计每个 task 的耗时。
  - `discovery_report.json`：扫描了多少条记录、选中多少 node 等。
- **注意**：`--output` 必须是新目录或空目录。
- **检查**：`wc -l $S/discovery/target_nodes.txt`，现有样本在 1,600 万到 3,000 万个 node 之间。

### 第 2 步：prepare（准备 run root）

**做什么**：

- 把当前 v4 代码复制一份冻结进 `$S/run/source/`，以后改 repo 不影响这个 run。
- 把 node 列表切成若干 **task**。一个 task 是一段连续的约 15,500 个 target node，由一个 builder 进程构建。
- 记录所有输入文件的指纹，写出 `config.json` 和 `run.sh`。

```bash
TASKS=$(( ($(wc -l < $S/discovery/target_nodes.txt) + 15499) / 15500 ))
H=/scratch/jshen/data/HG008_GIAB
sbatch -p general -c 1 --mem=8G -t 1:00:00 -J prepare -o $S/prepare-%j.out --wrap \
  "cd /scratch/jshen/Github/Pansoma && $PY -m $P.orchestrate prepare --root $S/run --tensors $S/tensors \
   --gam $GAM --nodes $S/discovery/target_nodes.txt --node-stats $S/discovery/node_stats.json \
   --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite --haplotypes 90 \
   --tasks $TASKS --processes 48 --gam-cache-mb 8192 \
   --snv-min-af 0.08 --indel-min-af 0.08 \
   --chromosomes autosome --chr-index $G/hprc-v1.1-mc-grch38.d9.chr_node_ranges.tsv \
   --merge-shard-size 32768 --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path \
   --somatic-vcf  $H/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz \
   --somatic-bed  $H/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed \
   --germline-vcf $H/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz \
   --germline-bed $H/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed \
   --reference-fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   --truth-dir $S/truth"
```

**参数说明**：

| 参数 | 含义 |
|---|---|
| `--haplotypes 90` | 这个 graph 的 haplotype path 数 H，包括 GRCh38 和 CHM13：HPRC v1.1 d9 是 90，v2.1 d46 是 464，只有 GRCh38 的 vg autoindex linear graph 是 17（GRCh38 + 16 条 path cover）。channel 6（path count）按 H 编码：H 条 path → 100，少 1–49 条 → 99–51，更少的 → 1–50，多于 H → 101–127（见 v4 README 第 8 节）。必填，换 graph 时一定要跟着改 |
| `--tasks` | task 个数；按每个 task 约 15,500 个 node 算 |
| `--processes` | 同时跑几个 task，也就是几个 builder 进程。决定内存，见第 3 步的表 |
| `--gam-cache-mb 8192` | 每个进程的 GAM group cache，8 GiB 足够 |
| `--snv-min-af 0.08 --indel-min-af 0.08` | **build 阶段**的 AF 阈值：低于它的 allele 根本不生成 tensor |
| `--chromosomes autosome` | 只保留 chr1–22 上的 target node；chrX/Y/M 和 unplaced 在构建前就去掉 |
| `--merge-shard-size 32768` | 跑完后按染色体合并，每个 shard 文件最多 32,768 个 tensor |
| （不加）`--keep-sources` | 加了会在合并后保留各 task 自己的目录，占用约多一倍磁盘；不加时校验通过后删掉 |
| 4 个 truth 文件 + `--reference-fasta` + `--truth-dir` | 合并后自动打标签。truth 表写到 `--truth-dir` |

**label 阶段的 AF floor**（`--label-snv-min-af`、`--label-indel-min-af`）：AF 更低的 tensor 保留在数据里，但标成 −1
（reason `below_snv_min_af` / `below_indel_min_af`）。build AF 已经是 0.08，floor 不超过 0.08 就没有作用，所以现在都
不加。加了的话会冻结进 `config.json`，finalize 时直接生效，跑完不用再单独 relabel。

**COLO829T 的 truth 文件**：

- somatic VCF：`/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz`
- somatic BED：同一目录下的 `SMaHT_v2_easy_difficult_extreme.union.bed`
- germline：`/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.vcf.gz` 和同目录的 `dipcall_hg38.dip.bed`

**注意**：`--root` 必须是不存在的新目录。prepare 读 6 GB 的 `node_stats.json` 约需 1 分钟、7 GB 内存，所以给了 8G。

### 第 3 步：run（构建所有 task，然后自动合并和打标签）

```bash
sbatch -p general --cpus-per-task=48 --mem=420G --time=14-00:00:00 \
       --output=$S/run/slurm-%j.out $S/run/run.sh
```

**做什么**：

1. 最多同时跑 `--processes` 个 task，预计最慢的先跑。每个 task 是一个独立的 builder 进程，结束后内存全部释放。
2. 每个 task 跑完立刻校验输出（shape、shard 长度、summary 与 manifest 是否一致）。
3. 全部完成后执行 **finalize**：按染色体合并并逐字节校验，再对所有 tensor 打标签。

**各平台设置**（现有 6 个 run 的实测，HG008 / COLO829T）：

| 平台 | `--processes` | 峰值内存 | 用过的 `--mem` | 墙钟时间：tasks + finalize |
|---|---|---|---|---|
| PacBio HiFi | 48 | 305 GiB | 335G | 4.7 + 1.5 小时 |
| Illumina | 48 | 247 / 194 GiB | 277G / 231G | 8.0 + 4.0 / 5.8 + 2.2 小时 |
| ONT-UL | 48 | 403 / 430 GiB | 433G / 448G | 7.2 + 3.1 / 24.2 + 4.6 小时 |
| fiberseq | 36 | 405 GiB | 454G | 4.8 + 2.5 小时 |

- fiberseq 用 48 进程 13 分钟就超过 454G 被杀，所以用 36。
- `--mem` 按峰值加约 30 GB 给。Slurm 不能降低运行中 job 的内存：要降就取消 run job，用新的 `--mem` 重新提交
  `run.sh --resume`，已完成的 task 会保留。
- 6 个 run 和别的 GAM job 同时在同一个存储上跑，墙钟时间是上限。
- `--cpus-per-task` 必须不小于 `--processes`，否则 run 会拒绝启动。
- finalize 里标签那一步单进程，17–36 分钟，内存 12–16 GiB；整个 finalize（合并、标签、删除 task 输出）1.5–4.6 小时。

**监控**：

```bash
cat $S/run/status.json            # status, tensors, tensors_by_type, merged, peak_sampled_rss_kib ...
$PY -c "import json, collections; q = json.load(open('$S/run/queue_status.json')); \
print(q['completed_tasks'], '/', q['total_tasks'], collections.Counter(t['status'] for t in q['tasks'].values()))"
tail $S/run/logs/task_0000.log    # 单个 task 的日志；task_NNNN.resources.txt 是它的时间和内存
```

**失败了怎么办**：

- **重跑没完成的 task**：用同样的 sbatch 命令，把最后的脚本改成 `$S/run/run.sh --resume`。已完成且校验通过的 task
  会跳过，半成品会被移到 `incomplete/`。
- **内存不够**：先在 `config.json` 里把 `processes` 改小（例如 48 改成 36），再 `--resume`。
- **提前挂一个 fallback job**：用 `--dependency=afternotok:<run 的 job id> --kill-on-invalid-dep=yes`，run 失败时
  自动降低进程数并 resume。run job 因为它自己的依赖（例如 prepare）失败而被取消时，`afternotok` 也会触发，所以
  fallback 脚本要先检查 `$S/run/config.json` 是否存在。`sbatch` 提交时就复制了脚本，提交后再改脚本不起作用。
- **只剩合并或标签失败**：`$PY -m $P.orchestrate finalize --root $S/run`。已完成的步骤会跳过，可以重复执行。

### 第 4 步：输出和检查

```
$S/tensors/SNV/  和  $S/tensors/INDEL/
  chrN_shard_NNNNN_data.npy      (n, 8, 200, 101) int8，每个文件最多 32,768 个 tensor
  chrN_variant_summary.ndjson    每行一个 tensor：site_id、candidate_id、AF、各类 read 数、row_groups、alleles[] ...
  chrN_shard_NNNNN_labels.npy    int8 标签，顺序与 shard 相同
  chrN_labels.ndjson             每个 tensor 的标签和 reason（为什么是这个标签）
  labels.manifest.json           各标签和各 reason 的计数、规则常数、rules_sha256、truth 文件的 SHA-256
  manifest.json                  每条染色体有哪些 shard（带 SHA-256）、tensor 格式
$S/tensors/non_autosomal/        chrX/Y/M 等（不用于训练）
$S/tensors/*.recall.tsv, truth_recall.json    每个 truth allele 有没有对应的 tensor
```

**检查**：

- 看 `labels.manifest.json` 里的 `totals`，和 `indexed_gam_pipeline_v4/README.md` 第 5 节的表对比，大致应在同一量级。
- 画一个 tensor 看看（base 环境的 matplotlib 和 NumPy 不兼容，所以用 hunyuanvideo15 环境的 python）：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/pansoma_matplotlib /wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python \
  scripts/visualize_tensor.py $S/tensors/SNV/chr1_shard_00000_data.npy -i 0 -o chr1_0.png
```

### 第 5 步（可选）：以后换规则或 floor 时重新打标签

```bash
sbatch -J relabel -o relabel-%j.out indexed_gam_pipeline_v4/tools/jobs/relabel.sh \
  $S/tensors SOMATIC_VCF SOMATIC_BED GERMLINE_VCF GERMLINE_BED $S/truth [SNV_MIN_AF [INDEL_MIN_AF]]
```

- 不留备份：`label` 只在一个 kind 的标签文件全部写完后才替换旧的，失败时旧标签保持不变。
- `''` 表示跳过其中一个 floor。
- 用 `/scratch/jshen/data/pansoma_v2_tensors/pipeline_code/indexed_gam_pipeline_v4/tools/jobs/relabel.sh` 时，跑的是
  `pipeline_code/` 里冻结的代码（commit 在 `pipeline_code/git_head.txt`）；用 repo 里的那个，跑的是当前 checkout。

---

## 三、训练模型（GPU，`gpu` partition，节点 tequila）

### 基本约定

- **一个模型**对应一个样本 × 一个平台 × 一个 kind（SNV 或 INDEL）。
- **数据划分**：
  - chr1 是测试集，训练和选模型都不用它。
  - chr2–22 按每 20,000 个 node 切成块，随机取 5% 的块做 validation。整块划分，保证 validation 和 training 不共享 reads。
- **训练数据**：只用标签 0/1/2。
- **评估数据**：validation 和 chr1 测试除了 0/1/2，还把 `off_reference_no_truth_match` 这类 −1 当作 non 计分，因为
  实际调用时会遇到它们。其它 −1（BED 外、低于 AF floor、没有 GRCh38 位置）不计分。
- **输入编码**：在 GPU 上把 8 个 channel 变成 36 个 plane。
  - read base、site allele、graph base、operation、strand 做 one-hot。
  - BQ、MAPQ、path count 做 masked z-score。统计量只从训练数据估计一次，存进模型。
  - 另外加几个标记 plane：`covered`、`bq_missing`、`differs`，以及每一行属于 A1/ALT/REF/OTHER 哪一组。
  - 所有 padding 位置严格为 0。
- **例子：HG008 Illumina SNV**（2026-09-28 的训练 job 376757）：
  - 训练集 2,560,638 个 tensor：non 2,274,667、somatic 9,178、germline 276,793。
  - validation 160,113 个：non 144,501、somatic 590、germline 15,022；另加 5,617 个 off-reference 当作 non。

### GPU 资源

tequila 上有：

- 2 块 H100（94 GB）：`--gres=gpu:h100:1`
- 5 个 24 GB MIG slice：`--gres=gpu:24gb:1`
- 6 个 12 GB slice

提交前先 `squeue -p gpu` 看有没有别的 job 占着 H100。

### 第 1 步：建 index cache（CPU，可选但推荐）

```bash
J=/scratch/jshen/data/pansoma_net_v2_runs/jobs
T=/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors     # 或者你自己的 $S/tensors
sbatch -J index -o $J/index-%j.out $J/build_index.sh $T
```

- **做什么**：读所有 summary 和 label，提取训练需要的字段（node、shard 位置、row_groups、scalars、reason），存到
  `/scratch/jshen/data/pansoma_net_v2_runs/index_cache/`。
- **好处**：GPU job 不用再花时间读这些 GB 级的 ndjson。
- **耗时**：约 5 分钟，12G 内存。
- **自动失效**：labels.manifest 变了（例如 relabel 之后），cache 会自动重建。

### 第 2 步：训练 + 在 chr1 上测试（一个 job 做完）

```bash
sbatch -p gpu --gres=gpu:h100:1 -c 16 --mem=40G -t 8:00:00 -J HG008_Illumina_SNV -o $J/%x-%j.out \
  $J/run.sh HG008_Illumina_SNV $T SNV \
  --epochs 12 --non-fraction 0.25 --class-weights sqrt --batch-size 256 --lr 2e-4 --num-workers 14
```

`run.sh` 的三个位置参数：NAME（输出放到 `/scratch/jshen/data/pansoma_net_v2_runs/NAME/`）、tensor set、kind。之后的
参数原样传给 `train`。

**参数说明**（这组设置和 2026-09-28 的 job 376757 一样）：

| 参数 | 含义 | 默认值 |
|---|---|---|
| `--epochs 12` | 训练轮数 | 30 |
| `--non-fraction 0.25` | 每个 epoch 用全部 somatic 和 germline，再加随机抽的 25% non，每轮重新抽 | 1.0 |
| `--class-weights sqrt` | loss 的类别权重：balanced 是 n / (3 · n_c)，sqrt 是它的平方根，比 balanced 温和 | balanced |
| `--batch-size 256` | 每块 GPU 的 batch 大小 | 64 |
| `--lr 2e-4` | 学习率 | 1e-4 |
| `--num-workers 14` | 读数据的进程数 | 8 |

**其它常用参数**：

- `--scalars`：把 site 的数值特征（log coverage、AF、ALT/REF/OTHER 数、allele 数等）也喂给模型。
- `--select ap`：按 somatic PR-AUC 选最佳 checkpoint。默认 f1 是按阈值化后的 somatic F1。
- `--val-chroms chr2`：改为整条染色体做 validation。
- `--resume last.pth`：接着训练。

**速度和耗时**：

- 默认开 bf16 和 torch.compile（编译约 2 分钟），单块 H100 约 1,400 tensors/s，显存峰值约 27 GiB。
- HG008 Illumina SNV 按上面的设置，每个 epoch 约 85 万个 tensor（9,178 + 276,793 + 25% × 2,274,667），算下来约
  10 分钟一轮，12 轮约 2–2.5 小时。这是按速度估算的。
- job 的时间上限给的是 8 小时。

### 第 3 步：两块 GPU 怎么用

**推荐：两个模型同时跑。** 把第 2 步的命令提交两次（例如一个 `SNV`，一个 `INDEL`），每个 job 占一块 H100。原因：

- 单卡开 compile 实测 1,404 tensors/s。
- 两卡 DDP 只测过不开 compile 的情况，才 926 tensors/s。

**一个模型用两块卡（DDP）：**

```bash
sbatch -p gpu --gres=gpu:h100:2 -c 32 --mem=60G -t 8:00:00 -J ddp -o ddp-%j.out --wrap \
  "cd /scratch/jshen/Github/Pansoma/machine_learning && torchrun --nproc_per_node=2 -m pansoma_net_v2.train --ddp \
   --tensors $T --kinds SNV --output /scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_SNV_ddp \
   --cache-dir /scratch/jshen/data/pansoma_net_v2_runs/index_cache \
   --epochs 12 --non-fraction 0.25 --class-weights sqrt --batch-size 256 --lr 2e-4"
```

- DDP 实测所有进程合计内存 47–50 GiB，所以给 60G。
- DDP 加 compile 还没测过。
- 这条命令只训练、不做 chr1 测试；训练完用第 6 步的 `predict` 命令测。

**24 GB MIG slice**：batch 只能到 64–96，约 125 tensors/s，比整块 H100 慢约 10 倍，只适合试跑。

### 第 4 步：监控

```bash
R=/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_SNV
tail -f $R/train.log                        # 每个 epoch 的 loss、速度、validation 结果
tail -3 $R/gpu_usage.csv                    # 每 10 秒一次的显存和 GPU 利用率
tail -3 $R/rss_all_user_processes.txt       # 每 30 秒一次，你账号所有进程的内存合计
```

### 第 5 步：输出怎么看

| 文件 | 内容 |
|---|---|
| `train.log`、`metrics.jsonl` | 每个 epoch 的训练 loss 和 accuracy、速度、GPU 显存峰值，以及 validation 报告：somatic 和 germline 的 average precision（PR-AUC），每类的 precision / recall / F1，somatic F1 最高时的阈值 t |
| `stats.json` | z-score 用的统计量 |
| `best.pth`、`last.pth` | checkpoint，含模型权重、统计量、`somatic_threshold` t、数据划分、训练参数；用 `PansomaNetV2.from_checkpoint(path)` 重建 |
| `test_chr1/<样本>.<set>.SNV.predictions.ndjson.gz` | chr1 每个 tensor 一行：chrom、candidate_id、label（truth）、test_label（参与计分的标签，不计分时为 null）、reason、off_reference、p_non、p_somatic、p_germline、pred |
| `test_chr1/*.metrics.json` | chr1 测试集上的同一套指标 |

`pred` 的规则：p_somatic ≥ t 就判为 somatic，否则取 non 和 germline 中概率较大的那个。

### 第 6 步：不用 job 脚本，直接运行

在 GPU 节点上，或者自己 `srun --pty` 进去以后：

```bash
cd /scratch/jshen/Github/Pansoma/machine_learning
P=/wanglab/jshen/anaconda3/bin/python          # torch 2.8 + timm
$P -m pansoma_net_v2.train --tensors $T --kinds SNV --output runs/X \
   --epochs 12 --non-fraction 0.25 --class-weights sqrt --batch-size 256 --lr 2e-4
$P -m pansoma_net_v2.predict --checkpoint runs/X/best.pth --tensors $T --kinds SNV --chroms chr1 --output runs/X/test_chr1
```

`--tensors` 可以一次给多个 set（例如同一平台的 HG008 和 COLO829T），训练一个跨样本的模型。README 的约定是一个样本一个
模型，所以跨样本的效果需要自己比较。

**channel 6 的 storage 必须一致**：一个模型只读一种 path count 编码（manifest 的 `tensor_storage`）。`int8-count-linear100-log2` 是加 `--haplotypes` 之前在 HPRC v1.1 d9 上建的 tensors（现在在 `backup_ch6_linear100_20261001/`），2026-10-01 以前的 checkpoint 都是用它训练的；`int8-count-haplotypes100` 是现在的 builder 写的（按 H 编码），现有的 6 套 tensor 都是这种。`train` 拒绝把两种混在一起，checkpoint 记下 `path_count_storage`；`predict` 和 `train --resume` 遇到另一种 storage 会在预测之前退出（没有记录的旧 checkpoint 算 `int8-count-linear100-log2`）。所以旧 checkpoint 不能预测新建的 tensors，要用新 tensors 重新训练。用 `--drop-planes path_count` 训练的模型不读 channel 6，不受这个限制。

### 第 7 步：改代码前后的测试

```bash
cd /scratch/jshen/Github/Pansoma && $PY -m unittest discover -s indexed_gam_pipeline_v4/tests -t .      # pipeline：168 个测试，约 70 秒
cd machine_learning && $PY -m unittest discover -s pansoma_net_v2/tests -t .                             # 模型
```

模型测试里的 `test_gpu` 只在 GPU 节点上运行，`test_real_data` 会读 COLO829T Illumina 的真实数据，约 70 秒。
