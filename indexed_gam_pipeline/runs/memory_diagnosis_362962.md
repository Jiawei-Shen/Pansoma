# Actual full-decode memory attribution (job 362962)

Read-only diagnostic, completed in 4m10s on tsingtao. Production job 362961 was not modified. Script: `tmp/profile_heavy_batch.py`; raw results: `tmp/memory_probe_362962.log`.

Reproduced first batch of task 22 from recovery 362305, using its frozen source, full decoding, 512 target nodes, MAPQ > 10 and 8 GiB GAM cache. This measures fetching, graph loading, full decoding and reference release, not candidate counting/tensor writing.

| Stage | RSS GiB |
|---|---:|
| Startup | 0.046 |
| GAM fetch complete | 8.297 |
| Graph loaded | 8.299 |
| Full decoding complete | 25.368 |
| Reads/by_node/alignments released, GC, GAM cache kept | 7.901 |
| Additional diagnostic malloc_trim, cache kept | 7.865 |
| Cache released, GC and trim | 2.945 |

3441 selected alignments contain 73,918,852 read bases. Graph lookup retrieves 8644 distinct nodes and only 61,553 sequence bases. The GAM query scans 78 groups, reports 82,402 decode operations and returns 4402 alignments before MAPQ filtering. Cached groups include records other than the selected targets; accounted cache bytes are 6,881,223,997 (~6.41 GiB). Cache accounting is not whole-process RSS.

Full decoding produces 96,533,414 Column objects, 15,110,694 Visit objects, and 7,309,616 Observations across complete alignments. Only 13,265 unique candidates lie in the target nodes. Column objects alone occupy 9.35 GiB by shallow-size accounting (104 bytes each); their lists add 0.764 GiB. Referenced coordinates, mapping metadata, observations and allele objects add further memory. These counts include context across the entire reads, not only the 512 target nodes. The added ~17.07 GiB during decoding is measured; a fully exclusive object-type heap attribution was not performed.

This batch's references release ~17.47 GiB without evicting the GAM cache. Additional trim frees only ~36 MiB, so persistent decoded-object references explain the dominant decode allocation here. This is not proof that every batch is free of leaks or fragmentation. About 2.945 GiB remains after cache release/trim in this diagnostic; its remaining allocations were not fully attributed.

## Why recovery exceeds the original job's peak

The original 17 unfinished partitions were repartitioned into 512 smaller tasks and still scheduled with 32 simultaneous builders. Total remaining nodes measures total work, not simultaneous full-read footprint. Target node batches do not cap complete-alignment bases or mapping/observation object counts.

Original 48-process job 362293: sampled at ~10m RSS sum 261.60 GiB, largest process 6.58 GiB, no process >15 GiB. At ~60m: 321.55 GiB, one process >15 GiB. At ~100m: 279.81 GiB, none >15 GiB. Later sampled points have only 1–2 processes >15 GiB.

Recovery 32-process job 362305: ~5m RSS sum 161.19 GiB; ~10m 325.49 GiB with 9 processes >15 GiB; ~20m 342.32 GiB with 10; ~60m 336.75 GiB with 10; ~80m 348.47 GiB with 10; ~100m 399.88 GiB; final sample 410.11 GiB. Slurm killed at charged memory 422.31 GiB versus 420 GiB limit.

Long-lived heavy tasks hold their active decoded batch while processing candidates; faster tasks finish and are replaced. Several heavy processes plateau (e.g. PID 4057170 27.41 -> 27.51 GiB from minute 10 to final sample), inconsistent with an assertion that all memory grows indefinitely in every process. Other tasks grow as their batches change. These samples support concentration of heavy live work plus cache and batch peaks, not simply cumulative output size.

The prior batch references survived until after the following GAM fetch and graph lookup. Production 362961 now explicitly releases them at each completed batch. This removes cross-batch overlap but cannot reduce the memory required while a single heavy batch remains active. At the diagnostic's completion, 362961 was RUNNING, 2/121 tasks complete, latest summed RSS ~396.64 GiB. No further production changes made.

Possible follow-up approaches compatible with full decoding: measure and avoid storing irrelevant candidate observations outside target nodes after validating support semantics; replace per-base Python object storage with compact full-length arrays; use a byte/base-aware working-set budget. None is implemented by this investigation.
