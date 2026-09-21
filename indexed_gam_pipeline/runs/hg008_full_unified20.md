# HG008 full unified-index run: 20 Python builders

Submitted 2026-09-21; job **362287**, node `meow`, allocation **20 CPUs / 400G**,
14-day wall-time limit. Source snapshot: **2ee56f2**. Jobs 362277 and 362278 were
cancelled at the user's request (the batch steps report failure due to cancellation).

Run directory:
`/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/full_unified20_20260921T214210Z`

- GAM: original HG008-T PacBio HiFi 116x AF-HPRC sorted GAM and matching GAI.
- GBZ: `hprc-v1.1-mc-grch38.d9.gbz`.
- Reuse all 16,560,350 discovered candidate nodes; no repeated discovery scan.
- 20 contiguous, disjoint partitions of 828,017 or 828,018 nodes.
- Shared read-only unified sequence/count SQLite index, built as stage 1.
- Python reader, ALT-support upper-bound filter, 8 GiB GAM cache per process.
- No optional node-index cache or nested candidate workers.
- Seven channels, 200 rows, 101 columns, int32; 2048 tensors per full shard.
- Outputs: `python/part_0/` ... `python/part_19/`, no merge or full-content comparison.

55 regression tests passed in 7.520 seconds, including an actual 20-process build
against one SQLite index and full/disjoint partition checks. Production includes
an independent graph-count audit followed by a 20-process bounded real-GAM
preflight before starting the full build.

Observed after submission: input/partition verification completed in **1.933 s**;
GBZ index creation running. A memory sample during indexing was **11.39 GiB**
summed process RSS. This is an initial observation, not the full-run peak or a
completed real-data preflight. Preparation/snapshot/partitioning took **57.575 s**.

Tracking files inside the run directory:

- `status.json`, `memory.ndjson` (per-process/tree RSS every 30 seconds).
- `slurm-362287.out`, `slurm-362287.err`.
- `graph_index_report.json` (written on successful index completion).
- `graph_audit.json` (written after independent count verification).
- `preflight/status.json`, `python/status.json`.
- `python/part_N/batch_timing.ndjson`, `manifest.json`, `run_report.json`.
- `python/part_N.resources.txt`, `job.resources.txt` (GNU time final resource stats).
- `config.json`, `submission.json`, frozen `source/`, and `run.sh`.

`run/active_run.json` points to this job. RSS samples may double-count shared pages
and miss short peaks; GNU time and Slurm accounting complement these samples.
