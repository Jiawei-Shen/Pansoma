# 从 GAM 到训练好的模型：tensor 生成与 pansoma_net_v2 训练

在这台服务器上从一个样本的 GAM 生成训练 tensor，再用它们训练和测试模型的完整步骤。命令都可以直接运行。
代码和数据的术语保留英文。现状部分是 2026-10-12 的状态。

- tensor pipeline：`indexed_gam_pipeline_v4/`（所有命令和选项见它的 `README.md`，标签规则见
  `tensor_postprocessing/README.md` 和 `analysis/truth_labels.md`）
- 模型：`machine_learning/pansoma_net_v2/`（`README.md`）；这台集群上的 Slurm job 脚本在
  `machine_learning/pansoma_net_v2/jobs/`（不进 git，说明见其中的 `README.txt`）

---

## 一、项目概览

### 1. 组成

| 目录 | 作用 |
|---|---|
| `indexed_gam_pipeline_v4/` | 把排好序、建好索引的 GAM 变成 tensor，然后按染色体合并、打标签 |
| `machine_learning/pansoma_net_v2/` | 直接在合并好的 tensor 上训练和预测，再出 VCF、和 truth 比 |
| `machine_learning/pansoma_net_v2/jobs/` | 训练 + chr1 测试、出 VCF、INDEL recipe 提交等 Slurm 脚本（本地副本，git 忽略；日志和 `submitted.txt` 在 `/scratch/jshen/data/pansoma_net_v2_runs/jobs/`） |
| `machine_learning/pansoma_net/` | 旧模型，只能读旧的 `.dat/.idx` tensor，留给旧数据和旧 checkpoint，新数据不用它 |
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
P25/P50/P75 mixture 的 set）。HPRC v2.1 d46 上另有多癌种的 set（`HCC1395T_*`、`HCC1954T_*`、`H1437T_*`、`H2009T_*`，
truth 是 SEQC2 / CASTLE，标签加了 `--ignore-repeat-indels`；还没有 `README.txt`，有的还没做完）。数量和资源见
`indexed_gam_pipeline_v4/README.md` 第 5 节。

- **参数**：`--haplotypes 90`（channel 6 的 storage 是 `int8-count-haplotypes100`），build AF 的 SNV 和 INDEL 都是
  0.08，标签不设 AF floor（floor 不超过 build AF 时没有作用）。
- **标签**：12 个 `labels.manifest.json`（6 套 × SNV/INDEL）是 2026-10-02/03 finalize 打的，`rules_sha256` 都是
  `5430533aa100`，即 `pansoma_v2_tensors/pipeline_code/`（82d0a1d）里那份 `truth_labels.py` 的 SHA-256。repo 里的
  `truth_labels.py` 之后加了 `--ignore-repeat-indels`、改了 docstring 和 `label_run` 的写法，所以 SHA-256 不同；
  不加这个选项时，分类的代码和那份相同。
- **检查标签是否和当前规则一致**：`rules_sha256` 是 `truth_labels.py` 整个文件的 SHA-256，任何改动（注释、docstring
  也算）都会改变它，所以它和 `sha256sum indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py` 不同并不说明标签
  过时，它只说明是哪一份文件打的标签。要知道规则改动会不会改变标签，比较标签本身：relabel 之前记下
  `<set>/tensors/<KIND>/labels.manifest.json` 的 `totals`、`reasons` 和 `sha256sum <set>/tensors/<KIND>/*_labels.npy`，
  relabel 之后再比一遍，都相同就是标签没变。六个数据集的计数也在 `analysis/truth_labels.md` 第 7 节和各样本的
  `README.txt` 里。
- relabel 正在进行时不要用那套数据训练：先确认对应的 relabel job 已经 COMPLETED。

### 3. 待决定的问题

1. **新样本没有现成的 v4 job 脚本串起整个流程。** `tools/jobs/` 有 `graph_prep.sh`、`gam_sort.sh`、
   `build_record_index.sh`、`relabel.sh`，没有 discovery → prepare → run → finalize 这一串。重建、linear / d46 和多癌种
   set 用的脚本在 repo 的 `tmp/`（不进 git）：`tmp/rebuild_ch6_20261001/<样本>/`、`tmp/linear_hg38_20261001/`、
   `tmp/hprc_v21_d46_tensors_20261003/`、`tmp/castle_hcc1395_20261006/`（finalize 用 `tmp/split_finalize_20261003/`）。
   下面第二部分是完整命令，可以整理成 `indexed_gam_pipeline_v4/tools/jobs/` 下的脚本。
2. **未被 git 跟踪的遗留目录。** `src/`（旧的编译产物 `.so`）、`build/`、根目录的 `__pycache__/` 已被 git 忽略，可以删。

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

- **GAM 从哪来**：还没有 GAM 时，先按根目录 `README.md` 第 1 节用 vg giraffe 比对（在这台机器上先
  `source scripts/use_vg.sh`），再用 `sbatch -J gam_sort -o LOG $P/tools/jobs/gam_sort.sh sample.gam` 排序并生成 GAI
  （`vg gamsort -i`，检查通过才发布 `sample.sorted.gam` 和 `.gai`）。长读长的 GAM 再用
  `sbatch -J gri -o LOG $P/tools/jobs/build_record_index.sh sample.sorted.gam` 建 `.gri`，fetch 时只读一个 batch 用到的记录。
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
- **几个 GAM**：一个样本的 GAM 太大、按文库分别排序时，`--gam` 给全部排好序的 GAM（`--gam A B ...`，prepare 也一样），
  结果和一个合起来的 GAM 相同。
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
| 4 个 truth 文件 + `--reference-fasta` + `--truth-dir` | finalize 时合并后打标签。truth 表写到 `--truth-dir` |

**label 阶段的选项**（冻结进 `config.json`，finalize 时直接生效，跑完不用再单独 relabel）：

- AF floor（`--label-snv-min-af`、`--label-indel-min-af`）：AF 更低的 tensor 保留在数据里，但标成 −1（reason
  `below_snv_min_af` / `below_indel_min_af`）。build AF 已经是 0.08，floor 不超过 0.08 就没有作用，所以现在都不加。
- `--label-ignore-repeat-indels`：truth 漏掉大部分同聚物 / STR INDEL 时（SEQC2、CASTLE）加，GRCh38 上同聚物 ≥ 7 bp
  或 STR 里没有 truth 的 INDEL tensor 标 −1 而不是 0（`analysis/truth_labels.md` 的 R13）。HG008 和 COLO829T 不加。

**COLO829T 的 truth 文件**：

- somatic VCF：`/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz`
- somatic BED：同一目录下的 `SMaHT_v2_easy_difficult_extreme.union.bed`
- germline：`/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.vcf.gz` 和同目录的 `dipcall_hg38.dip.bed`

**注意**：`--root` 必须是不存在的新目录。prepare 读 6 GB 的 `node_stats.json` 约需 1 分钟、7 GB 内存，所以给了 8G。

### 第 3 步：run（构建所有 task），然后单独的 finalize job（合并和打标签）

```bash
sbatch -p general --cpus-per-task=48 --mem=420G --time=14-00:00:00 \
       --output=$S/run/slurm-%j.out $S/run/run.sh --no-finalize
# run 的 status.json 是 complete 之后（或者 --dependency=afterok:<run 的 job id>）：
sbatch -p general -c 9 --mem=22G -t 2-00:00:00 -J finalize -o $S/finalize-%j.out --wrap \
  "cd $S/run/source && $PY -m $P.orchestrate finalize --root $S/run"
```

**做什么**：

1. 最多同时跑 `--processes` 个 task，预计最慢的先跑。每个 task 是一个独立的 builder 进程，结束后内存全部释放。
2. 每个 task 跑完立刻校验输出（shape、shard 长度、summary 与 manifest 是否一致）。
3. 全部完成并校验后，`--no-finalize` 让 run 停在 `complete`。**finalize** 是单独的 job，用 run 冻结的代码
   （`$S/run/source`）：按染色体合并（`min(8, processes)` 个 worker）并逐字节校验，删掉 task 输出，再对所有 tensor 打
   标签。run job 的峰值是 194–430 GiB（下表），finalize 只要约 15 GiB，加上召回表扫描的 8 个 worker 约 2.4 GiB，所以给 22G。
   不加 `--no-finalize` 时 run job 最后自己做 finalize，这几个小时里一直占着大内存。

**各平台设置**（现有 6 个 run 的实测，HG008 / COLO829T）：

| 平台 | `--processes` | 峰值内存 | 用过的 `--mem` | 墙钟时间：tasks + finalize |
|---|---|---|---|---|
| PacBio HiFi | 48 | 305 GiB | 335G | 4.7 + 1.5 小时 |
| Illumina | 48 | 247 / 194 GiB | 277G / 231G | 8.0 + 4.0 / 5.8 + 2.2 小时 |
| ONT-UL | 48 | 403 / 430 GiB | 433G / 448G | 7.2 + 3.1 / 24.2 + 4.6 小时 |
| fiberseq | 36 | 405 GiB | 454G | 4.8 + 2.5 小时 |

- fiberseq 用 48 进程 13 分钟就超过 454G 被杀，所以用 36。
- `--mem` 按峰值加约 30 GB 给。Slurm 不能降低运行中 job 的内存：要降就取消 run job，用新的 `--mem` 重新提交
  `run.sh --no-finalize --resume`，已完成的 task 会保留。
- 6 个 run 和别的 GAM job 同时在同一个存储上跑，墙钟时间是上限。
- `--cpus-per-task` 必须不小于 `--processes`，否则 run 会拒绝启动。
- 这 6 个 run 的 finalize（它们冻结的代码里召回表扫描只有一个进程）：标签那一步 17–36 分钟、12–16 GiB；整个
  finalize（合并、标签、删除 task 输出）1.5–4.6 小时。

**监控**：

```bash
cat $S/run/status.json            # status, tensors, tensors_by_type, merged, labeled, peak_sampled_rss_kib ...
$PY -c "import json, collections; q = json.load(open('$S/run/queue_status.json')); \
print(q['completed_tasks'], '/', q['total_tasks'], collections.Counter(t['status'] for t in q['tasks'].values()))"
tail $S/run/logs/task_0000.log    # 单个 task 的日志；task_NNNN.resources.txt 是它的时间和内存
```

**失败了怎么办**：

- **重跑没完成的 task**：用同样的 sbatch 命令，把最后的脚本改成 `$S/run/run.sh --no-finalize --resume`。已完成且校验
  通过的 task 会跳过，半成品会被移到 `incomplete/`。
- **内存不够**：先在 `config.json` 里把 `processes` 改小（例如 48 改成 36），再 `--resume`。
- **提前挂一个 fallback job**：用 `--dependency=afternotok:<run 的 job id> --kill-on-invalid-dep=yes`，run 失败时
  自动降低进程数并 resume。run job 因为它自己的依赖（例如 prepare）失败而被取消时，`afternotok` 也会触发，所以
  fallback 脚本要先检查 `$S/run/config.json` 是否存在。`sbatch` 提交时就复制了脚本，提交后再改脚本不起作用。
- **合并或标签失败**：重新提交上面的 finalize job。已完成的步骤会跳过，可以重复执行。

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
- 画一个 tensor 看看（base 环境的 matplotlib 和 NumPy 不兼容，所以用 hunyuanvideo15 环境的 python；不给 `-o` 时画到
  repo 的 `tensor_figures/<shard 名>/`，git 忽略）：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/pansoma_matplotlib /wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python \
  scripts/visualize_tensor.py $S/tensors/SNV/chr1_shard_00000_data.npy -i 0 -o chr1_0.png
```

### 第 5 步（可选）：以后换规则、truth 或 floor 时重新打标签

```bash
sbatch -J relabel -o relabel-%j.out indexed_gam_pipeline_v4/tools/jobs/relabel.sh \
  $S/tensors SOMATIC_VCF SOMATIC_BED GERMLINE_VCF GERMLINE_BED $S/truth [SNV_MIN_AF [INDEL_MIN_AF]]
```

- 资源写在 repo 这份脚本的 `#SBATCH` 里：`-c 9`（召回表扫描的 8 个进程加主进程）、`--mem=22G`、6 小时。
- 不留备份：`label` 只在一个 kind 的标签文件全部写完后才替换旧的，失败时旧标签保持不变。
- `''` 表示跳过其中一个 floor。`IGNORE_REPEAT_INDELS=1 sbatch ...` 加上 `--ignore-repeat-indels`；reference path 默认用
  当前标签的（`labels.manifest.json` 的 `reference_path`），也可以用环境变量 `REFERENCE_PATH` 指定。
- 用 `/scratch/jshen/data/pansoma_v2_tensors/pipeline_code/indexed_gam_pipeline_v4/tools/jobs/relabel.sh` 时，跑的是
  `pipeline_code/` 里冻结的代码（commit 在 `pipeline_code/git_head.txt`，没有 `IGNORE_REPEAT_INDELS`）；用 repo 里的那个，
  跑的是当前 checkout。
- relabel 会改写 `labels.manifest.json` 的 `created`（和 `rules_sha256`），pansoma_net_v2 的 index cache 因此会自动重建。

---

## 三、训练模型（GPU，`gpu` partition，节点 tequila）

### 基本约定

- **一个模型**对应一个样本 × 一个平台 × 一个 kind（SNV 或 INDEL）。
- **数据划分**：
  - chr1 是测试集，训练和选模型都不用它。
  - chr2–22 按每 20,000 个 node 切成块，随机取 5% 的块做 validation。整块划分，保证 validation 和 training 不共享 reads。
- **训练数据**：只用标签 0/1/2。两个 recipe 还用 `--ignore-reasons residual_partial_somatic_truth
  allele_partial_somatic_truth` 把部分匹配 somatic 的 tensor 在训练和 validation 里当 −1（标签文件不变）。
- **评估数据**：validation 和 chr1 测试除了 0/1/2，还把 `off_reference_no_truth_match` 这类 −1 当作 non 计分，因为
  实际调用时会遇到它们。其它 −1（BED 外、低于 AF floor、没有 GRCh38 位置等）不计分。`--ignore-reasons` 的 tensor 在测试
  时默认也不计分，只有 INDEL 测试把部分匹配 somatic 的 tensor 当 1 计分。
- **输入编码**：在 GPU 上把 8 个 channel 变成 36 个 plane。
  - read base、site allele、graph base、operation、strand 做 one-hot。
  - BQ、MAPQ、path count 做 masked z-score。统计量只从训练数据估计一次，存进模型。
  - 另外加几个标记 plane：`covered`、`bq_missing`、`differs`，以及每一行属于 A1/ALT/REF/OTHER 哪一组。
  - 所有 padding 位置严格为 0。
- **例子：HG008T Illumina SNV**（现在的标签，SNV recipe；`analysis/hg008t_snv_indel_7fc7396_20261012` 里
  `HG008_Illumina_SNV_7fc7396_ddp4_b1024` 的 `train.log`）：
  - 训练集 1,903,802 个 tensor：non 1,620,829、somatic 7,134、germline 275,839。
  - validation 118,218 个：non 103,012、somatic 469、germline 14,737；另加 5,044 个 off-reference 当作 non。
  - `--ignore-reasons` 去掉训练 1,720 个、validation 97 个 tensor。

### GPU 资源

tequila 上有：

- 2 块 H100（94 GB）：`--gres=gpu:h100:1`
- 5 个 24 GB MIG slice：`--gres=gpu:24gb:1`
- 6 个 12 GB slice

提交前先 `squeue -p gpu` 看有没有别的 job 占着 H100。

### job 脚本

```bash
JS=/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2/jobs   # 脚本（git 忽略，说明见 $JS/README.txt）
J=/scratch/jshen/data/pansoma_net_v2_runs/jobs                          # 日志和 submitted.txt
T=/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors       # 或者你自己的 $S/tensors
```

脚本里写死了这台集群的路径：代码 `/scratch/jshen/Github/Pansoma/machine_learning`，run 目录
`/scratch/jshen/data/pansoma_net_v2_runs/<NAME>`，index 和 chunk cache 在 `pansoma_net_v2_runs/index_cache/`。

### 第 1 步：建 index cache 和 chunk cache（CPU，可选但推荐）

```bash
sbatch -J index -o $J/index-%j.out $JS/build_index.sh $T
sbatch -p general -c 4 --mem=8G -t 6:00:00 -J chunks -o $J/chunks-%j.out --wrap \
  "/wanglab/jshen/anaconda3/bin/python $JS/build_chunks_from_manifest.py $T /scratch/jshen/data/pansoma_net_v2_runs/index_cache/chunks"
```

- **index cache**：读所有 summary 和 label，提取训练需要的字段（node、shard 位置、row_groups、scalars、reason），存到
  `pansoma_net_v2_runs/index_cache/`，GPU job 不用再花时间读这些 GB 级的 ndjson。HG008 Illumina 大小的 set 约 5 分钟，
  脚本给 12G 内存。`labels.manifest.json` 的 `format`、`rules_sha256`、`created`、`tensors` 或 AF floor 变了（例如
  relabel 之后），cache 会自动重建。
- **chunk cache**：训练从 zstd 压缩的 chunk 读 tensor（`--chunk-cache`，默认 `<cache-dir>/chunks`）。缺的 chunk
  `train` 会在训练前自己建；提前建可以和 relabel 同时跑（只读 `manifest.json`，默认 chr2–22，`KINDS='SNV INDEL'`）。

### 第 2 步：训练 + 在 chr1 上测试

**SNV**（一个 job 做完训练和 chr1 测试）：

```bash
sbatch -p gpu --gres=gpu:h100:1 -c 8 --mem=40G -t 8:00:00 -J HG008_Illumina_SNV -o $J/%x-%j.out \
  $JS/run.sh HG008_Illumina_SNV $T SNV \
  --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75 \
  --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
sbatch --dependency=afterok:<上面的 job id> -J vcf_HG008_Illumina_SNV -o $J/vcf_HG008_Illumina_SNV-%j.out \
  $JS/vcf.sh HG008_Illumina_SNV SNV                                       # chr1 的 VCF、PoN、rtg vcfeval
```

`run.sh` 的三个位置参数：NAME（输出放到 `/scratch/jshen/data/pansoma_net_v2_runs/NAME/`）、tensor set、kind。之后的
参数原样传给 `train`；训练完用 `best.pth` 预测 chr1。`vcf.sh RUN KIND` 从预测文件名找到 set，truth、BED、FASTA 和
reference path 读它的 `labels.manifest.json`，graph index 读它的 `manifest.json`。

**INDEL**：`submit_indel.sh` 按 recipe 给每个 set 提交三个连起来的 job（chunk cache → `run.sh` → `vcf.sh`），job id
追加到 `$J/submitted.txt`。直接运行它（它自己调用 `sbatch`），不要用 `sbatch` 提交：

```bash
$JS/submit_indel.sh HG008T_Illumina                       # run 名 HG008_Illumina_INDEL_small_nofilt_lr1e4
GPU=h100 TAG=mytag $JS/submit_indel.sh HG008T_PacBio COLO829T_ONT
```

默认在 24 GB MIG slice 上跑（`--micro-batches 4`），`GPU=h100` 用整块 H100；训练 job `-c 8`，`--mem` 短读长 36G、
长读长 48G。run 目录已经存在时它会停下，换一个 `TAG`。

**参数说明**（两个 recipe，`machine_learning/pansoma_net_v2/README.md` 的 "SNV recipe" 和 "INDEL recipe"）：

| 参数 | SNV recipe | INDEL recipe | 含义 | 默认值 |
|---|---|---|---|---|
| `--depths` / `--dims` | `3 3 27 3` / `96 192 384 768`（50.4 M 参数） | `2 2 6 2` / `64 128 256 512`（8.7 M） | 模型大小 | `3 3 27 3` / `192 384 768 1536`（199.6 M） |
| `--epochs` | 20 | 12 | 训练轮数 | 30 |
| `--batch-size` | 1024 | 1024 | 每块 GPU 的 batch 大小 | 64 |
| `--lr` | 1e-4 | 1e-4 | 学习率 | 1e-4 |
| `--class-weights` | `pow:0.75` | `sqrt` | loss 的类别权重：balanced 是 n / (3 · n_c)，sqrt 是它的平方根，`pow:E` 是它的 E 次方 | balanced |
| `--non-fraction` | 1.0 | 1.0 | 每个 epoch 用全部 tensor | 1.0 |
| `--block` | v2 | v2 | ConvNeXt block | v2 |
| `--ignore-reasons` | 两个部分匹配 somatic 的 reason | 同左 | 这些 reason 的 tensor 在训练和 validation 里当 −1 | 不设 |

`--weight-decay 0.05`、`--warmup-epochs 1` 是默认值。出 VCF 时 PASS 用 checkpoint 的 validation recall 0.9 阈值
（`graph_vcf` 的默认），SNV 的 off-reference tensor 用 `predict` 的 rescue；然后 PoN 和 rtg vcfeval。

**其它常用参数**：

- `--micro-batches K`：每个 batch 分 K 份做 forward 和 backward（梯度累积），更新和整批一样，显存只要一份的。
- `--val-chroms chr2`：改为整条染色体做 validation。
- `--resume last.pth`：接着训练（保留 checkpoint 的统计量）。
- `--drop-planes path_count`：不读 channel 6（见第 6 步最后一段）。
- 选模型：`best.pth` 是 validation 上阈值化后 somatic F1 最高的 epoch（同分取 validation loss 低的），存下那个阈值。
- `--scalars`、`--block v1`、`--decay-all`、`--keep-non-af` 等是留给旧 checkpoint 和实验的选项，recipe 不用（README）。

**速度和耗时**（tequila 实测）：

- 默认开 bf16 和 torch.compile（编译约 2 分钟）。
- SNV recipe 一块 H100 约 3,100 tensors/s：HG008T Illumina（约 190 万个训练 tensor）每个 epoch 约 10 分钟，20 轮加
  chr1 测试约 3.7 小时，RSS 峰值约 33 GiB，所以给 40G、8 小时。
- INDEL recipe 大小的模型（8.7 M，batch 1024，12 轮）在 24 GB MIG slice 上，训练加 chr1 测试：HG008T Illumina 33 分钟
  （`--micro-batches 16`，RSS 峰值 21 GiB）、PacBio 2 小时 43 分（`--micro-batches 4`，28 GiB）。

### 第 3 步：两块 GPU 怎么用

**推荐：两个模型同时跑。** 每个 job 占一块 H100（例如一个 SNV，一个别的平台的 SNV）。原因（默认大小的模型，
`machine_learning/pansoma_net_v2/README.md` 的 GPU runs 表）：

- 单卡 batch 256 开 compile 实测 1,404 tensors/s。
- 两卡 DDP 每块 batch 384、不开 compile 只有 926 tensors/s。两卡 DDP 加 compile（`train` 在 DDP 下也默认 compile）在
  tequila 上没测过。

**一个模型用两块卡（DDP）：**

```bash
sbatch -p gpu --gres=gpu:h100:2 -c 32 --mem=60G -t 8:00:00 -J ddp -o $J/ddp-%j.out --wrap \
  "cd /scratch/jshen/Github/Pansoma/machine_learning && torchrun --nproc_per_node=2 -m pansoma_net_v2.train --ddp \
   --tensors $T --kinds SNV --output /scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_SNV_ddp \
   --cache-dir /scratch/jshen/data/pansoma_net_v2_runs/index_cache \
   --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 512 --lr 1e-4 --class-weights pow:0.75 \
   --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth"
```

- `--batch-size` 是每块 GPU 的：两块卡给 512，总 batch 才是 recipe 的 1024。
- 默认大小的模型两卡 DDP（不开 compile）实测所有进程合计内存 47–50 GiB，所以给 60G。
- 这条命令只训练、不做 chr1 测试；训练完用第 6 步的 `predict` 命令测。

**24 GB MIG slice**：默认大小的模型 batch 只能到 64–96，约 125 tensors/s，比整块 H100 慢约 10 倍；INDEL
recipe 的 8.7 M 模型用 `--micro-batches 4` 在上面跑 batch 1024（`submit_indel.sh` 的默认）。

### 第 4 步：监控

```bash
R=/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_SNV
tail -f $R/train.log                        # 每个 epoch 的 loss、速度、validation 结果
tail -3 $R/gpu_usage.csv                    # 每 10 秒一次的显存和 GPU 利用率
tail -3 $R/rss_all_user_processes.txt       # 每 30 秒一次，你账号所有进程的内存合计
cd /scratch/jshen/Github/Pansoma/machine_learning && $PY -m pansoma_net_v2.summary HG008_Illumina_SNV
                                            # 每个 epoch 一行的表、现在在哪一步、chr1 测试和 vcfeval 结果
```

### 第 5 步：输出怎么看

| 文件 | 内容 |
|---|---|
| `train.log`、`metrics.jsonl` | 每个 epoch 的训练 loss 和 accuracy、速度、GPU 显存峰值，以及 validation 报告：somatic 和 germline 的 average precision（PR-AUC），每类的 precision / recall / F1，somatic F1 最高时的阈值 t，和 truth VCF 比的结果 |
| `stats.json` | z-score 用的统计量 |
| `best.pth`、`last.pth` | checkpoint，含模型权重、统计量、`somatic_threshold` t、数据划分、训练参数、channel 6 的 storage；用 `PansomaNetV2.from_checkpoint(path)` 重建 |
| `test_chr1/<样本>.<set>.SNV.predictions.ndjson.gz` | chr1 每个 tensor 一行：chrom、candidate_id、label（truth）、test_label（参与计分的标签，不计分时为 null）、in_test、reason、ignored（`--ignore-reasons` 去掉的）、excluded（`--exclude-candidates` 去掉的）、off_reference、truth_ids、p_non、p_somatic、p_germline、pred；SNV 另有 off-reference rescue 的 p_offref_* 和 offref_call |
| `test_chr1/*.metrics.json` | chr1 测试集上的同一套指标 |
| `vcf_chr1/` | `vcf.sh` 的 graph VCF、GRCh38 VCF（放不到 GRCh38 上的在 `.linear.unplaced.vcf.gz`），以及 `eval_bed`（benchmark BED 内）和 `eval_nobed`（整条 chr1）：PoN 标记过的 VCF 和 rtg vcfeval 的 `report.txt` / `report.json` |

`pred` 的规则：p_somatic ≥ t 就判为 somatic，否则取 non 和 germline 中概率较大的那个。

### 第 6 步：不用 job 脚本，直接运行

在 GPU 节点上，或者自己 `srun --pty` 进去以后：

```bash
cd /scratch/jshen/Github/Pansoma/machine_learning
P=/wanglab/jshen/anaconda3/bin/python          # torch 2.8 + timm
$P -m pansoma_net_v2.train --tensors $T --kinds SNV --output runs/X \
   --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75 \
   --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
$P -m pansoma_net_v2.predict --checkpoint runs/X/best.pth --tensors $T --kinds SNV --chroms chr1 --output runs/X/test_chr1
```

出 VCF 和 vcfeval 的三条命令（`graph_vcf`、`linear_vcf`、`vcfeval`）见 `machine_learning/pansoma_net_v2/README.md` 的
"Calls"。

`--tensors` 可以一次给多个 set（例如同一平台的 HG008 和 COLO829T），训练一个跨样本的模型。README 的约定是一个样本一个
模型，所以跨样本的效果需要自己比较。

**channel 6 的 storage 必须一致**：一个模型只读一种 path count 编码（manifest 的 `tensor_storage`）。现在所有 set 都是
`int8-count-haplotypes100`（按 H 编码）；`int8-count-linear100-log2` 是加 `--haplotypes` 之前在 HPRC v1.1 d9 上建的
tensors（那批数据已经删了），只留给 2026-10-01 以前训练的 checkpoint。`train` 拒绝把两种混在一起，checkpoint 记下
`path_count_storage`；`predict` 和 `train --resume` 遇到另一种 storage 会在预测之前退出（没有记录的旧 checkpoint 算
`int8-count-linear100-log2`）。所以 2026-10-01 以前的 checkpoint 不能预测现在的 tensors。用
`--drop-planes path_count` 训练的模型不读 channel 6，不受这个限制。

### 第 7 步：改代码前后的测试

```bash
cd /scratch/jshen/Github/Pansoma && $PY -m unittest discover -s indexed_gam_pipeline_v4/tests -t .      # pipeline：206 个测试，约 2 分钟
cd machine_learning && $PY -m unittest discover -s pansoma_net_v2/tests -t .                             # 模型：76 个测试
```

模型测试里的 `test_gpu` 只在 GPU 节点上运行，`test_real_data` 会读 COLO829T Illumina 的真实数据，约 70 秒；
`test_gpu` 也检查了 compile + channels_last 在 fp32 下和 eager 的结果一致。
