# Follow-up: legacy count audit and image visualization — 2026-09-18

Rebuilt legacy tensors for the same three nodes using `--variant-type all`:
`tmp/indexed_gam_legacy_comparison_all/`. Both formats produce five candidates.
The original bundled three-tensor legacy example was SNP-only.

Four candidates match exactly in coverage/ALT/REF/other counts. Node 2189 C>A
has identical coverage 72 but legacy ALT/other 71/1 versus v2 70/2. Independent
raw GAM edit traversal identifies one ALT observation with BQ 3, MAPQ 60:
`m84039_240114_012401_s1/24708040/ccs`, reverse mapping 202, substitution T
(forward A). Legacy filters mean ALT BQ (37.14); v2 requires each ALT observation
to pass BQ 10. This accounts for the entire difference. All three independently
recounted SNPs match both summaries under their respective quality policies.
Insertion boundary 138 corresponds to legacy preceding-base position 137;
both insertion allele counts match exactly. AF rounding is also accounted for.

The existing `scripts/visualize_tensor.py` now delegates to the shared module
`src/pangenome_ml_data_generation/tensors/visualization.py`, which supports both
formats. Generated and inspected five v2 and five matching legacy PNGs. Portable
images, compact summary and comparison evidence are in
[samples/README.md](samples/README.md). The source tensors and large debug metadata
remain server-local. The existing 2.3 KB compressed legacy regression fixture is
preserved; no new tensor shards or bulk GAM/graph data are included in synchronization.

Validation: **21 pipeline tests passed**, including the added BQ-3 regression;
**5 visualization tests passed**, checking version selection, complete candidate
ranges, padding trimming, first-row BQ zero, gap-quality and MAPQ masks, distinct
per-row references, PNG validity, metadata matching, and legacy rendering.
The image tests used NumPy 1.26.4 / Matplotlib 3.10.8 in the existing
`/wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python` environment. The default
Python's pre-existing NumPy/Matplotlib ABI incompatibility was not modified.

```bash
python indexed_gam_pipeline/compare_formats.py --legacy tmp/indexed_gam_legacy_comparison_all/variant_summary.ndjson --v2 tmp/indexed_gam_candidate_v2_verified/variant_summary.ndjson --gam tmp/HG008_pacbio_test.sorted.gam --index tmp/indexed_gam_current/rebuilt.gai --node-sqlite /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13_node_index.sqlite --output tmp/indexed_gam_candidate_v2_verified/legacy_comparison.json
python -m unittest discover -s indexed_gam_pipeline/tests -v
MPLCONFIGDIR=/tmp/pansoma_matplotlib /wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python -m unittest discover -s tests -p test_tensor_visualization.py -v
```

No full genome-scale run, model training or inference was performed.

---

# Candidate-centered v2 server validation — 2026-09-18

Implemented default `--format candidate-v2`, with explicit `--format legacy`.
The existing uncommitted server implementation and its ten regression tests were
preserved. The historical results below refer to the five-channel format.

## Checks completed

`python -m unittest discover -s indexed_gam_pipeline/tests -v`: **20 tests passed**.
These include all ten pre-existing regressions plus:

- Original-edit SNP extraction, inclusive 50 bp insertion/deletion support,
  explicit 51 bp and complex-replacement rejection, and merging adjacent I/D
  edits before enforcing the limit.
- Multi-node consumed slices with no boundary padding, different upstream and
  downstream branches, row-specific references, and row-local context insertions.
- Reverse SNP, insertion and deletion identity and tensor orientation, including
  reversal of base-quality order.
- Shared insertion blocks, left-aligned different alleles, full candidate flags,
  gap versus missing-coverage padding, node-edge and partial-interior boundaries.
- Position-specific coverage, partial interval overlap, deletion spans,
  insufficient insertion evidence, repeated node visits counted once per record,
  and retention of separate duplicate records.
- Full counts and AF invariant to row caps; deterministic selection under input
  reordering; overlapping indexed runs merged without multiplying records.
- Synthetic shards, summary indices, manifest shape/dtype, legacy tensor parity,
  and the existing bundled HG008 five-channel sample checksums.

Real indexed retrieval was independently checked against a sequential scan of
**only the 10,000-record HG008 test GAM**, using the previously rebuilt matching
GAI. Result: **290 expected and retrieved records; zero missing/extra records**;
node-segment multisets also matched. Query: 2 merged runs, 7 groups, 7,000 decoded
records. Report: `tmp/indexed_gam_candidate_v2_validation/validation_report.json`.

The final v2 build queried nodes **577, 1841, 2189**, retaining 290 complete
MAPQ-qualified alignments, 5,154,410 read bases, and 4,373 graph context nodes.
It found 24 candidate alleles, filtered 19, and saved **five tensors** in one
`(5, 6, 200, 100)` int16 shard. There were no unsupported events in this subset.

| Candidate (forward, zero-based) | Coverage | ALT | REF | Other |
|---|---:|---:|---:|---:|
| 577:138 insertion A | 115 | 14 | 92 | 9 |
| 577:138 insertion AT | 115 | 8 | 92 | 15 |
| 577:138 T>A | 115 | 92 | 23 | 0 |
| 1841:0 T>C | 117 | 117 | 0 | 0 |
| 2189:0 C>A | 72 | 70 | 0 | 2 |

Containing-node segment counts are not used as candidate coverage: node 577 has
117 retrieved segments but only 115 records reach these candidate positions.
The real subset exercises SNPs and insertions; deletion and 50/51 bp behavior
were tested synthetically, not established from real deletion examples here.

## Inspectable outputs

Final artifacts are local, ignored, and **not committed or uploaded**:

- `tmp/indexed_gam_candidate_v2_verified/shard_00000_data.npy`
- `tmp/indexed_gam_candidate_v2_verified/variant_summary.ndjson` (full row paths
  and column-to-graph mappings enabled)
- `tmp/indexed_gam_candidate_v2_verified/manifest.json`
- `tmp/indexed_gam_candidate_v2_verified/paired_rows.txt`
- `tmp/indexed_gam_candidate_v2_verified/inspection_report.json`

The inspection audit checked all 534 occupied rows, shard/summary/manifest
consistency, full counts/AF, candidate flags, insertion reference gaps, and zero
unused rows. **All 52,626 non-gap reference columns matched the corresponding
oriented graph base**, across 24 nodes represented inside the windows. The
sample manifest in `samples/manifest.json` keeps the legacy sample contract and
adds a `candidate_v2` section with provenance, parameters, per-candidate counts,
source GAM/GAI hashes and final artifact hashes. The text preview displays actual
read/reference pairs, with `.` padding and `-` alignment gaps.

Reproduction (choose a new output directory):

```bash
python indexed_gam_pipeline/run.py validate --gam tmp/HG008_pacbio_test.sorted.gam --index tmp/indexed_gam_current/rebuilt.gai --nodes tmp/indexed_gam_current/known_nodes.txt --output tmp/my_v2_validation
python indexed_gam_pipeline/run.py build --gam tmp/HG008_pacbio_test.sorted.gam --index tmp/indexed_gam_current/rebuilt.gai --nodes tmp/indexed_gam_current/known_nodes.txt --node-sqlite /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13_node_index.sqlite --output tmp/my_v2_samples --batch-nodes 3 --min-variants 3 --variant-type all --shard-size 32 --debug-rows
python indexed_gam_pipeline/inspect_tensor.py tmp/my_v2_samples tmp/my_v2_samples/paired_rows.txt
```

An initial slower trial was interrupted; counting was then changed to inspect
only the anchor mapping and immediate neighbors, and bounded tensor windows
retain multi-node context. Subsequent small builds completed. Intermediate
trial directories are not the final sample provenance.

**No full genome-scale run was performed.** No throughput/memory scaling study,
truth labeling, model training/inference, or graph-path-only variant discovery
was performed. Normalization establishes exact forward-node allele identity,
including reverse complements and adjacent I/D edits; repeat-shift and cross-node
allele equivalence are not implemented. Complex replacements remain unsupported
candidates. A central deletion interval plus its internal context insertions that
cannot fit the requested width fails explicitly instead of truncating it.

---

# 本地测试记录

日期：2026-09-17。

## 当前工作区：HG008 真实 GAM 和图索引

使用 `tmp/HG008_pacbio_test.sorted.gam`（10,000 alignments）和
`/scratch/jshen/data/AF-Filtered_VG_Indexes/` 中匹配 GFA 的节点 SQLite
索引测试。原有 `.gai` 与当前 GAM 不匹配：随机读取时索引区间落在 GAM group
内部。运行 `run.py index` 生成 `tmp/indexed_gam_current/rebuilt.gai` 后：

- 完整 `discover` 扫描 10,000 alignments，观察到 93,602 个 nodes，选中
  20,171 个候选。`node_stats.json` 对全部 93,602 个 nodes 的 `perfect`、
  `not_perfect`、`max_read_length`、`max_cigar_length` 与旧
  `tmp/HG008_pacbio_test.unperfect_nodes.pkl` 逐项一致。另取前 30 个做过
  探索性测试。
- `validate` 对已知的 577、1841、2189 三个候选 node 做独立全 GAM 扫描；
  命中 290 alignments，缺失和多余均为 0，per-node segments 完全相同。
- `build` 从索引 GAM 和图 SQLite 直接生成 3 个 `(5,201,100)` int8 tensors，
  保存在 `tmp/indexed_gam_current/tensors_rebuilt/`。三个 variant key、ALT/REF
  计数、覆盖度和 AF 与现有 NPU 生成的 HG008 tensor summary 一致。
  tensor 数组未按行顺序逐元素一致；旧 NPU 使用原始 GAM，read 行顺序不同。
  按完整五通道 read 行作为多重集合比较，三个 tensors 与旧 NPU 输出全部
  逐字节一致。新生成的 `.npy` 和 summary 与本目录 `samples/` 中压缩的
  sample tensors、summary 也逐字节一致。
- 从原始 GFA 读取这三个 node 的序列，与 SQLite 中的序列逐一一致。
- 另对发现的前 30 个 nodes 完成直接构建；这些低 ID nodes 在默认筛选下
  没有达到 tensor 阈值，因此输出 0 个 tensor。这个现象不代表管线失败。
- 10 项回归测试通过，包含重建索引、随机读取、SQLite 序列冲突检测、
  发现阶段旧统计 schema 对照、NPU tensor 核心对照及真实样本五通道校验。

当前重现命令：

```bash
python indexed_gam_pipeline/run.py index --gam tmp/HG008_pacbio_test.sorted.gam --output tmp/indexed_gam_current/rebuilt.gai
python indexed_gam_pipeline/run.py discover --gam tmp/HG008_pacbio_test.sorted.gam --output tmp/indexed_gam_current/discovery_schema
python indexed_gam_pipeline/run.py validate --gam tmp/HG008_pacbio_test.sorted.gam --index tmp/indexed_gam_current/rebuilt.gai --nodes tmp/indexed_gam_current/known_nodes.txt --output tmp/indexed_gam_current/validation_rebuilt
python indexed_gam_pipeline/run.py build --gam tmp/HG008_pacbio_test.sorted.gam --index tmp/indexed_gam_current/rebuilt.gai --nodes tmp/indexed_gam_current/known_nodes.txt --node-sqlite /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13_node_index.sqlite --node-json tmp/HG008_pacbio_test_candidate_nodes.json --output tmp/indexed_gam_current/tensors_synced --batch-nodes 3 --min-variants 3 --variant-type snp --shard-size 32
python -m unittest discover -s indexed_gam_pipeline/tests -v
```

上述命令的输出目录和 `.gai` 文件需使用新的路径，已有文件不会覆盖。

## 历史 COLO829T 索引检验记录

输入为仓库内 `tmp/COLO829T_3M.sorted.gam`，索引为
`tmp/COLO829T_3M.sorted.gam.gai`，GAI version 1。
使用 `pangenome-ml-data-generation` 环境，protobuf 3.20.3、pysam 0.24.0、
NumPy 2.5.1。未调用 `vg`，未生成 NPU 文件。

先扫描前 100,000 条 alignments，统计到 199,748 个 nodes，按原 node 筛选
阈值 0.05 选取 node ID 最小的 12 个候选，用于小规模索引查询验证。
这份候选清单是测试子集，不能代替完整生产候选清单。

```text
233 329 457 596 600 1196 1198 1219 1236 1283 1522 1587
```

然后对上述 nodes 读取 `.gai` 指定的 GAM 区间，独立顺序扫描整份 GAM 作为基准。

| 检查项 | 结果 |
|---|---:|
| 顺序扫描 alignment 总数 | 3,000,000 |
| 目标 nodes | 12 |
| 索引查询 runs / groups | 21 / 22 |
| 索引查询实际解码 alignments | 22,000 |
| 精确命中的 alignments | 18 |
| 顺序扫描期望命中的 alignments | 18 |
| 缺失 / 多余 alignments | 0 / 0 |
| 每个 node 的 segment 多重集合 | 全部一致 |
| MAPQ 过滤后的 segments 合计 | 18 |
| 索引查询耗时 | 2.306 秒 |
| 顺序扫描核对耗时 | 409.021 秒 |

这里比较完整 alignment 内容的 SHA-256 多重集合，同时比较 node segment 内容的
多重集合；不按 read name 去重。耗时仅针对这次小规模查询及本机环境，
不是全基因组 tensor pipeline 的性能结论。

原始输出保存在仓库内：

- `tmp/indexed_gam_smoke/discovery/discovery_report.json`
- `tmp/indexed_gam_smoke/discovery/target_nodes.txt`
- `tmp/indexed_gam_smoke/validation/validation_report.json`

复现命令见本目录 README 的“本地简单测试”。

## 早期合成回归检查记录

1. 索引查询与顺序扫描一致，覆盖重复记录、同名 reads 和跨 node mappings。
2. 不支持的索引版本和截断 varint 明确失败。
3. Read cursor、MAPQ 边界、末尾 BQ=0 和负链转换。
4. 批次 node 归属，避免跨批次重复计数。
5. 共享 tensor 核心与未修改的旧 `.dat` reader 输出逐元素一致，覆盖 SNP、
   insertion、deletion 和负链输入；`.dat` 测试数据仅在内存中构造。
6. 使用合成图序列与合成 GAM 检查 shards、summary 索引和报告格式。

```bash
python -m unittest discover -s indexed_gam_pipeline/tests -v
```

## 未测试范围

按用户要求，本地不运行 COLO829T 的完整 tensor 生成，也不执行 truth labeling、
模型推理或线性坐标 VCF 输出。匹配的图序列仅在 server 上可用。
HG008 的三个真实 tensor 已与旧 NPU 结果按完整五通道 read 行核对。
COLO829T 的完整 tensor 生成、全基因组吞吐量与内存评估仍未执行。
