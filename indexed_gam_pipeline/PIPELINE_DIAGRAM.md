# Current tensor-building pipeline — diagrams

Drawn from the working-tree code on 2026-09-22 (`run.py`, `build_v2.py`,
`gam_reader.py`, `graph_index.py`, `candidates.py`, `candidate_work.py`,
`edit_columns.py`, `split_outputs.py`, `tensor_storage.py`,
`dynamic_tensor_run.py`). Prose reference: [CURRENT_PIPELINE.md](CURRENT_PIPELINE.md),
[UNIFIED_GRAPH_INDEX.md](UNIFIED_GRAPH_INDEX.md), [DYNAMIC_TENSOR_RUN.md](DYNAMIC_TENSOR_RUN.md).

Default format `indexed-gam-candidate-v4`, schema 4, shape `(7, 200, 101)`,
`int8`, storage encoding `int8-count-div4-v1`.

---

## 1. Offline / one-time inputs

```mermaid
flowchart LR
    subgraph Upstream["Upstream, outside this pipeline"]
        FQ["Reads"] --> MAP["vg giraffe / map"]
        MAP --> GAM0["GAM"]
        GAM0 --> SORT["vg gamsort"]
    end

    SORT --> GAM["Sorted BGZF GAM"]
    SORT --> GAI[".gam.gai  (GAI v0/v1)"]
    GAM -.->|"run.py index, if .gai absent"| GAI

    GBZ["GBZ graph, e.g. hprc-v1.1-mc-grch38.d9.gbz"] --> CPP["gbz_graph_index.cpp<br/>compiled by build_gbz_query.py --graph-index"]
    CPP --> BGI["build_graph_index.py"]
    BGI --> GIDX[("all_graph_nodes.gbz.sqlite<br/>node_id, seq, distinct_path_count<br/>+ graph_metadata")]

    GAM --> DISC["run.py discover<br/>sequential scan, MAPQ &gt; 5<br/>imperfect-mapping fraction &gt; 0.05"]
    DISC --> NODES["target_nodes.txt<br/>sorted positive node IDs"]

    GAM --> BUILD(["run.py build --format candidate-v4"])
    GAI --> BUILD
    GIDX --> BUILD
    NODES --> BUILD
```

`--graph-index` is the current path: one SQLite file replaces the old
`--node-sqlite` sequence source **and** the per-process `--gbz/--gbz-query/--occurrence-cache`
GBWT occurrence cache. Mixing the two input styles is rejected.

---

## 2. Production orchestration (`dynamic_tensor_run.py`)

```mermaid
flowchart TD
    PRIOR[("Prior verified run:<br/>graph index, discovery node list,<br/>passed graph audit + preflight")]
    PRIOR --> PREP["dynamic_tensor_run prepare<br/>--tasks 512 --processes 32 --cache-mb 8192<br/>--snv-min-af .06 --indel-min-af .08"]
    PREP --> FROZEN["run root:<br/>source/ frozen copy + per-file SHA-256<br/>config.json, parts/ node files, run.sh"]
    FROZEN --> SBATCH["sbatch run.sh<br/>32 CPUs, 420 GiB, OMP/BLAS threads = 1"]
    SBATCH --> CTRL["dynamic_tensor_run run<br/>verify inputs + graph-index stamp"]
    CTRL --> QUEUE{{"execute_queue: 512 contiguous node tasks,<br/>at most 32 concurrent, fail-fast"}}
    QUEUE --> T0["task 0000"]
    QUEUE --> T1["task 0001"]
    QUEUE --> TN["... task 0511"]
    T0 --> B0["fresh builder process<br/>run.py build (Section 3)"]
    B0 --> EXIT["builder exits<br/>caches + allocator memory released"]
    EXIT --> VAL["validate_shards per output type<br/>shape, dtype int8, allowed event types"]
    CTRL --> MEM["memory.ndjson: process-tree sample / 30 s<br/>python/status.json: PIDs, timings, state"]
    VAL --> OUT[("SNV/task_NNNN/ and INDEL/task_NNNN/<br/>plus python/task_NNNN/ shared coordinator")]
```

One builder per task, one candidate worker per builder; shards are never merged.
Any task failure stops the queue. There is no automatic resume — recovery is a
new prepared run over the unfinished node intervals (`runs/recover_*.py`).

---

## 3. The builder: per-batch dataflow (the core)

`run.py build` → `build_v2.build` → `_dispatch` → `_build` → `_build_impl`.

```mermaid
flowchart TD
    START(["build(args)"]) --> VALID["split_enabled: validate SNV/INDEL dirs + AF thresholds<br/>check max_indel_len, width, cache limits"]
    VALID --> MODE["select_decode_mode: mean length of first 10 GAM reads<br/>&lt; 1000 bp -&gt; 'full'  or  &gt;= 1000 bp -&gt; 'window'"]
    MODE --> OPEN["open readers once for the whole build:<br/>IndexedGam (LRU group cache, 8 GiB) or VgGam<br/>GraphIndex read-only SQLite, single BEGIN snapshot"]
    OPEN --> MANI["write manifest.json / run_report.json (status=running)"]
    MANI --> BATCH{{"batches(nodes, --batch-nodes 512, --max-node-span 10000)"}}

    subgraph PerBatch["for each target-node batch"]
        direction TB
        S1["<b>1. GAM fetch</b><br/>GAI bins intersecting the batch -&gt; merge virtual-offset runs<br/>seek BGZF, decode groups, node-&gt;record postings<br/>keep MAPQ &gt; 10 and optional --chr<br/>collect context_nodes = targets + every visited node"]
        S1 --> S2["<b>2. Graph lookup</b><br/>GraphIndex.get_nodes(context_nodes), chunks of 900 IDs<br/>-&gt; forward sequence + distinct GBWT path count<br/>missing node = hard error"]
        S2 --> S3["<b>3. Decode edits</b> (decode_alignment per record)<br/>walk read + graph cursors, merge adjacent I / adjacent D<br/>reverse mappings mapped to forward-node coordinates<br/>emit Columns / Visits / Observations, release the protobuf"]
        S3 --> S3b["candidates = observations on <b>target</b> nodes only<br/>by_node index: node -&gt; reads<br/>unsupported_events.ndjson: indel &gt; 50 bp, complex replacement"]
        S3b --> S4["<b>4. Early ALT prefilter</b> (--early-alt-filter)<br/>alt_support_bounds: one vote per record per candidate, BQ &gt;= 10<br/>bound &lt; --min-variants 3 -&gt; reject with alt_support_upper_bound"]
        S4 --> S5["<b>5. Support counting</b> per surviving candidate<br/>tasks of &lt;= 16 candidates per node -&gt; candidate_work<br/>overlap(read, candidate): one vote per record,<br/>ALT &gt; REF &gt; other, earliest mapping breaks ties"]
        S5 --> S6{"coverage = ALT+REF+other<br/>AF = ALT / coverage<br/>ALT &gt;= 3 and AF &gt;= threshold<br/>and variant_type allowed?"}
        S6 -->|no| FILT["filtered_candidates.ndjson<br/>with reasons + counts"]
        S6 -->|yes| S7["<b>6. make_tensor</b> (Section 4)"]
        S7 --> S8["<b>7. Buffer and shard</b><br/>2048 tensors -&gt; shard_XXXXX_data.npy<br/>+ variant_summary.ndjson rows"]
        S8 --> S9["<b>8. Release batch state</b><br/>clear by_node / reads / alignments before the next fetch<br/>append batch_timing.ndjson"]
    end

    BATCH --> S1
    S9 --> BATCH
    BATCH -->|"all batches done or --max-tensors reached"| FIN["final flush, timings,<br/>gam_group_cache + graph-index counters,<br/>manifest status=complete"]

    subgraph Split["--snv-output / --indel-output (single decode, two writers)"]
        SP["SplitOutputs.append: route by event_type"]
        SP --> SNV[("SNV/  AF &gt;= 0.06")]
        SP --> IND[("INDEL/  AF &gt;= 0.08")]
    end
    S8 -.->|split mode| SP
```

Key invariants visible in the code:

- **One GAM pass per batch, shared by both variant types.** `SplitOutputs` only
  fans out at write time; SNV and INDEL never re-fetch or re-decode.
- **Only the parent process touches GAM / SQLite / GBZ.** With `--workers > 1`,
  forked candidate workers see the decoded batch copy-on-write and return
  `(tensor, metadata)` in original candidate order.
- **Counts use all eligible records**, before the 200-row cap.
- Timings accumulate into `gam_fetch_seconds`, `graph_index_seconds`,
  `decode_edits_seconds`, `candidate_tensors_and_shard_writes_seconds`.

### Decode mode

```mermaid
flowchart LR
    A["decode_alignment(mode=...)"] --> F["full: one Column object per alignment column<br/>~104 B each; simple, memory-heavy on long reads"]
    A --> W["window: EditColumns stores one EditRun per edit<br/>columns materialized lazily through ColumnView<br/>split / reference_support computed from edit intervals"]
```

Production HG008 runs currently pin `full` (`config['params']['decode_policy']`),
which is the dominant memory term — see
[runs/memory_diagnosis_362962.md](runs/memory_diagnosis_362962.md).

---

## 4. One candidate → one `(7, 200, 101)` tensor

```mermaid
flowchart TD
    C["Candidate = (node, forward start, REF, ALT, kind)<br/>kind in SNP / INS / DEL"] --> E["eligible = [(read, support, anchor visit)]<br/>from overlap(), all records, no cap yet"]
    E --> W["<b>anchor_window</b> per record — 'anchor-centered-columns-v1'<br/>candidate start pinned at column width//2 = 50<br/>50 preceding columns, anchor, 50 following<br/>long I/D cropped at the window edge -&gt; omitted_context<br/>INS: G slots only with boundary evidence, else all-zero padding<br/>DEL: one D column per deleted graph base<br/>reverse anchor -&gt; row reversed/complemented into candidate orientation"]
    W --> M["mismatch bp = substituted + inserted + deleted bases visible in the window"]
    M --> SORT["<b>sort</b> by  -mismatch bp,  -MAPQ,  record SHA-256,  anchor mapping index"]
    SORT --> GRP["<b>group</b> stably by the visible, candidate-oriented node path<br/>groups in first-occurrence order; order inside a group preserved"]
    GRP --> SAMP["<b>sample</b> up to 200 rows: index floor(i*(N-1)/(K-1)); K=1 -&gt; N//2<br/>row-selection revision 'window-edit-bp-group-uniform-v1'"]
    SAMP --> FILL["fill 7 channels, int8"]
    FILL --> OUTT[("tensor + metadata:<br/>coverage, alt/ref/other, AF, candidate_columns,<br/>row_groups, omitted_context, versions")]
```

| Channel (1-based) | Contents |
|---|---|
| 1 | Read base — A=1, C=2, G=3, T=4, N=5, gap=6; padding 0 |
| 2 | Base quality, clipped to −1…127; −1 = missing quality / no read base |
| 3 | Flags — 1 = difference, 2 = candidate region, 3 = both |
| 4 | Alignment MAPQ, clipped to −1…127; present on aligned gaps |
| 5 | Operation — M=1, X=2, I=3, D=4, complex=5, aligned no-insertion gap=6 |
| 6 | Graph reference base for that row/column, same base encoding |
| 7 | `min(distinct GBWT path count of the column's node // 4, 127)` |

Insertions and insertion-gap slots take their **anchor** node's count;
deletions keep their mapped node's count; missing evidence is zero in all
seven channels.

---

## 5. Outputs and verification

```mermaid
flowchart LR
    BLD["builder"] --> NPY["shard_XXXXX_data.npy — (&lt;=2048, 7, 200, 101) int8"]
    BLD --> SUM["variant_summary.ndjson — identity, coverage, ALT/REF/other, AF,<br/>candidate_columns, row_groups, shard location, versions"]
    BLD --> FC["filtered_candidates.ndjson — rejection reasons"]
    BLD --> UE["unsupported_events.ndjson — oversized indels, complex replacements"]
    BLD --> MF["manifest.json / run_report.json — args, channels, encodings,<br/>index provenance, cache stats, stage timings"]
    BLD --> BT["batch_timing.ndjson — per-batch stage seconds, counters"]
    BLD --> TN["target_nodes.txt"]

    NPY --> VS["validate_shards — shape / dtype / event types (in-run)"]
    NPY --> VE["validate_examples.py — independent GAM + graph re-query<br/>needs --debug-rows"]
    NPY --> VZ["scripts/visualize_tensor.py -&gt; 7-panel PNGs"]
```

`--debug-rows` adds per-row record hashes, full read paths and per-column graph
coordinates; it is required by the independent auditor and is disabled on full
production runs because of its size.
