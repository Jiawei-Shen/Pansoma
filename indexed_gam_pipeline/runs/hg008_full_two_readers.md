# Full HG008, two readers with two independent processes each

Launched 2026-09-21 after the user approved the full runs and requested a simpler
workflow: each process generates its own tensors; **no full tensor comparison
and no output merge/rewrite**. Both jobs were confirmed RUNNING on `meow`.
The original job 362252 remains STOPPED.

| Backend | Job | Processes | Nodes per process |
|---|---|---:|---:|
| Python reader | 362277 | 2 | 8,280,175 |
| VG v1.77.0 | 362278 | 2 | 8,280,175 |

Each process handles its own contiguous half of the 16,560,350-node discovery
result. The first half has 16,193 batches, the second 16,191, with the existing
512-node / 10,000-ID-span limits. Neither full run has a read/node/tensor count
cutoff. Full GAM discovery (14,014,144 records) and graph sequence indexing are
reused, with fingerprints checked.

Only ALT-support upper-bound filtering is enabled. Candidate node-index cache
is off; internal candidate workers = 1. Each independent Python builder has an
8-GiB GAM group cache (16 GiB per backend allocation total), plus other memory.
The VG backend uses `/scratch/jshen/bin/vg_v1.77.0`, exact-node filtering and
complete alignment decoding. The two processes own separate GBWT caches and
output directories. Debug row metadata and the external memory sampler are off.
Each Slurm job requests 4 CPUs, 128 GiB and a 14-day limit; that limit is not a
runtime estimate.

## Preflight

49 regression tests passed, including peer-process cleanup on failure and
partition-boundary output handling. Preflight 362275 completed successfully on
128 nodes from each actual half (256 nodes total), producing 46 tensors. Python
and VG two-process output matched the validated serial reference in tensor bytes
and normalized metadata. The source implementation is commit `64c9746`.

The first 4-GiB-per-process preflight (362272) was stopped deliberately after
identifying a cache working-set concern, not because tensor correctness failed.
Its completed serial reference was retained and checked. The 8-GiB parallel
preflight reused it with unchanged builder implementation and input nodes.

The preflight exercised merging as a correctness test, but merging and global
comparison were subsequently removed from the full-run launch at the user's
request. The submitted 12-line wrapper uses the same tested `run_builders`
function; it does not invoke the helper's full merge/compare modes.

## Outputs and progress

Root:
`/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/full_reader_comparison_20260921T202021Z`

- Python tensors: `python/part_0/`, `python/part_1/`.
- VG tensors: `vg/part_0/`, `vg/part_1/`.
- Each part writes `shard_XXXXX_data.npy`, `variant_summary.ndjson`, its manifest,
  rejection logs and per-batch stage times in `batch_timing.ndjson`.
- NPY shape is `(N,7,200,101)`, int32, N=2048 except each part's final partial shard.
  Consumers should enumerate both part directories for a backend; filenames are
  intentionally local to each directory.
- Backend controllers: `python/status.json` and `vg/status.json`; per-process
  logs and GNU time resource reports are in those backend directories.
- Controller status updates every 30 seconds, per-batch times at each completed
  batch. A failing process stops its peer and preserves outputs/logs.
- End-of-build checks validate each part's shard shape/count and summary
  consistency, without comparing the two backends' tensor contents.

[Launch commands](hg008_full_two_readers/full_submission.json),
[preflight result](hg008_full_two_readers/preflight_status.json),
[configuration/provenance](hg008_full_two_readers/comparison_config.json), and
[actual submitted wrapper](hg008_full_two_readers/run_partition_build.py).
