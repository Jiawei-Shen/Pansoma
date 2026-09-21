# HG008 bounded candidate optimization benchmark — 2026-09-21

## Status and provenance

- Full job **362252 remains stopped** (`scancel --signal=STOP --full`); no automatic resume or replacement full build was submitted. `scontrol suspend` was denied. STOP preserves the process and its files but may continue consuming allocation/wall time.
- Debug comparison and independent audit: job **362255**, source **dd3863d**.
- Production metadata comparison: job **362258**, source **1915e71** (same builder implementation; additional benchmark mode).
- Inputs: the supplied HG008-T PacBio HiFi AF-HPRC sorted GAM/GAI, matching d9 GBZ, and the existing complete graph SQLite. Discovery and graph indexing were reused.
- Exactly **1,024 discovery nodes**, two batches of 512, **11,144 raw candidates**, no tensor-count cutoff. Early filtering rejects **10,556** candidates; **588** proceed to the original full support/AF filters.
- Each run produces **182 tensors** (89 DEL, 17 INS, 76 SNP), `int32 (182,7,200,101)` in one permitted partial final shard; configured shard capacity remains 2,048.
- 44 regression tests passed, including duplicate records/repeated visits, reverse orientation, I/D/BQ, FIFO bounds, and deterministic one/two-worker outputs. The expected missing-node-99 fixture error is not a production error.

## Timings

Times below are builder wall seconds. Candidate time includes filtering, index/pool overhead, tensor calculation, and shard/summary writing. Each subprocess has a separate GNU time resource report; raw JSON records wrapper time as well. Memory is the maximum sampled process-tree PSS (five-second samples), not a guaranteed peak. Summed RSS double-counts pages shared after fork.

| Mode | GAM fetch | Graph | GBWT | Decode | Candidate + write | Total | Sampled PSS GiB | Sampler seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline_before | 98.34 | 1.39 | 5.90 | 19.92 | 158.57 | 285.98 | 6.19 | 4.89 |
| early_alt_only | 83.95 | 0.89 | 4.57 | 19.79 | 20.18 | 130.82 | 6.15 | 1.67 |
| indexed_one_worker | 80.89 | 1.07 | 4.30 | 20.63 | 19.31 | 127.13 | 6.16 | 1.71 |
| indexed_two_workers | 95.75 | 0.72 | 4.47 | 20.15 | 22.94 | 145.29 | 7.23 | 2.54 |
| baseline_after | 82.39 | 0.90 | 5.00 | 19.94 | 156.38 | 265.57 | 6.18 | 5.01 |
| production_one_worker | 51.48 | 1.03 | 4.08 | 19.32 | 10.32 | 87.30 | 5.42 | 1.06 |
| production_two_workers | 67.47 | 0.14 | 3.85 | 19.77 | 10.33 | 102.74 | 6.08 | 1.63 |

All first five modes use `--debug-rows` for independent auditing. The last two modes disable it, matching the intended full-run configuration. The sampler is part of the benchmark controller, not the production builder; its recorded wall cost runs concurrently and must not simply be subtracted from build time.

## Findings and worker recommendation

1. Early ALT filtering provides the main gain. With debug output held constant, original candidate/write time is 156–159 seconds versus 20.18 seconds for early filtering alone, and 19.31 seconds with a FIFO index. Whole-build speedup for indexed single-worker mode is about **2.09–2.25x** in this bounded comparison. This is not a genome-wide runtime forecast.
2. FIFO indexing has a small additional effect on this sample. Its peak retained estimate is 8,977,424 bytes (~8.56 MiB), with 240 node entries in the busiest batch. The configured cap is 1,000 entries and 64 MiB estimated containers across all workers, not per worker. Batch boundaries clear indices. The 1,000-entry eviction behavior is also tested synthetically.
3. Two workers really process two candidate tasks concurrently. Each task contains at most 16 candidates from one node, with at most two submissions in flight. Inputs are initially shared via Linux fork/COW; the parent alone fetches GAM, reads graph/GBZ, decodes, and writes shards. The per-batch pools and process-to-parent results have overhead; task sizes after prefiltering are small.
4. Debug metadata magnifies IPC: the debug summary is **517,078,603 bytes**, compared with **351,928 bytes** in production mode; tensor NPY size is **102,939,328 bytes** in either mode. Debug metadata was therefore not used as the final basis for choosing a production worker count.
5. In production mode, candidate/write time is **10.32 seconds (one worker)** versus **10.33 seconds (two)**: no measurable gain in this run. GAM fetch varies **51.48 versus 67.47 seconds**, explaining most of the total-time difference. Decode remains about 19–20 seconds. Worker summed timers are not wall time and are not a direct parallel speedup measure. Scheduling/COW/IPC overhead can absorb the small available candidate-stage gain; this benchmark does not individually profile those overhead components.
6. Use **one worker** for now. Two workers increase sampled production PSS from **5.42 to 6.08 GiB** without reducing candidate wall time. Larger/more complex regions or other machines could behave differently. Short-read performance and full-genome timing are not established by this long-read sample.

Explicit recommended build options (in addition to the existing HG008 arguments):

```text
--early-alt-filter --node-index-cache-nodes 1000 --node-index-cache-mb 64 --workers 1
```

CLI defaults remain opt-in. This report does not restart the stopped process: its existing code snapshot cannot acquire these changes via SIGCONT. A future optimized full-run launch needs a new tested source snapshot and separate staging/recovery plan; existing discovery/index files can be reused.

## Correctness

- All debug configurations have identical SHA-256 hashes for every tensor NPY and the entire `variant_summary.ndjson`.
- Independent raw-GAM/graph audit passed for **182 candidates, 10,707 selected rows, 1,076,175 reference columns**, including GBWT count cells, support/AF, padding, row selection/grouping and candidate flags. Audit wall time: **134.00 seconds**.
- Both production configurations have identical tensor and production-summary hashes; their tensors also match the independently audited debug tensors.
- Anchor lookup/window construction, mismatch bp sorting, stable node-path grouping, uniform sampling, seven channels, width 101, and filter thresholds are unchanged. The declined anchor/window optimization was not added.
- Early-rejected log records contain `alt_support_upper_bound` and `coverage_not_evaluated: true`; they intentionally do not claim complete REF/other/AF evaluation. Accepted candidate metadata is unchanged.

## Reproducibility records

- [Debug timing/commands/checksums](hg008_candidate_optimization/debug_status.json)
- [Production timing/commands/checksums](hg008_candidate_optimization/production_status.json)
- [Independent source audit](hg008_candidate_optimization/debug_source_audit.json)
- [Fixed node list](hg008_candidate_optimization/debug_nodes.txt)
- Per-batch timings and snapshot provenance are alongside those files.
- External debug output: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/candidate_optimization_20260921T182037Z`.
- External production output: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/candidate_production_bench_20260921T183601Z`.
