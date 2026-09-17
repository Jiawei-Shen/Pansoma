# 本地测试记录

日期：2026-09-17。

## 真实数据：通过

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

## 合成回归检查：6 项全部通过

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
真实 tensor 的新旧结果比较、全基因组吞吐量与内存评估留到 server 进行；
README 已提供正式命令及资源要求。
