# Full and edit-window decoding comparison, 2026-09-21

Both implementations remain available through `--decode-mode full|window|auto`.
The default auto policy uses the arithmetic mean length of the first ten
unfiltered records at the start of the input GAM: below 1,000 bp selects full,
otherwise window. The decision is fixed per build, independent of node query
order and GAM reader backend, and written to the manifest. Short inputs sample
all available records. Explicit modes skip sampling.

Window mode retains edit intervals and all candidate/unsupported-event discovery
semantics. It computes reference support from intervals and crops boundary edits
before materializing candidate windows. Complete-read hashes, mapping paths,
repeated visits, row selection and output schemas are preserved. It does not
skip compressed GAM decoding or protobuf parsing. Quality slices, edit strings
and observations still cover the complete alignment.

## Correctness tests

- Existing candidate semantics tests also run with window decoding.
- 100 deterministic randomized alignment cases, each with four mappings including
  repeated node visits, both orientations, M/X/I/D/complex edits, varying or
  absent quality, and widths 1, 2, 11, 101 and 201. Compare columns, visits,
  observations, unsupported events, support, windows, tensor values and metadata.
- A 100,000-bp alignment materializes fewer than 250 columns for a 101-column
  candidate window, with identical output.
- A 10,000-bp background insertion on both orientations is not expanded during
  reference support checking; bounded windows retain identical overflow metadata.
- Build integration at 999, 1,000 and 10,000 bp, full/window/auto, two candidate
  workers: tensor files and summary bytes agree. Auto threshold, first-ten-only,
  mixed lengths, fewer than ten records and empty input are tested separately.

Run with:

```sh
/wanglab/jshen/anaconda3/bin/python -m unittest discover -s indexed_gam_pipeline/tests
```

## Real HG008 data benchmark (Slurm 362288)

Four regions from the two node partitions: 438 alignment records in total.
The first regions reuse the earlier 128-read samples. Additional regions use
node-list chunks 57 and 49 (512 nodes per chunk; these are not necessarily the
production batch numbers when max-node-span splits batches).

Each region runs in separate processes: full/window, then window/full, plus auto.
The table averages the two explicit-mode wall times and uses their maximum RSS.
Runs enable debug row audit, shard size 16, and stop at 32 tensors. The original
GAM retrieval is performed during preparation; builds query the miniature GAM
and read precomputed graph sequences/counts. Original GAM I/O and GBZ startup
are excluded. These are bounded correctness/performance samples, not production
wall-time or memory forecasts.

| Region | Reads | Tensors | Full wall s | Window wall s | Full peak MiB | Window peak MiB |
|---|---:|---:|---:|---:|---:|---:|
| Part 0, initial | 128 | 24 | 4.929 | 2.352 | 433.7 | 207.1 |
| Part 1, initial | 128 | 23 | 3.892 | 1.496 | 390.6 | 154.9 |
| Part 0, chunk 57 | 128 | 32 | 30.851 | 29.847 | 2531.3 | 2458.3 |
| Part 1, chunk 49 | 54 | 32 | 14.266 | 14.411 | 912.9 | 834.5 |

All five runs per region produced identical SHA256 hashes for every tensor shard,
variant summary including row audit, filtered-candidate log, and unsupported-event
log. Auto selected window for all four long-read regions. All shard validations
passed (111 reference tensors across four regions).

| Region | Full decode s | Window decode s | Full candidate/write s | Window candidate/write s |
|---|---:|---:|---:|---:|
| Part 0, initial | 3.199 | 0.693 | 1.166 | 1.353 |
| Part 1, initial | 2.871 | 0.437 | 0.700 | 0.834 |
| Part 0, chunk 57 | 8.225 | 7.101 | 20.527 | 20.644 |
| Part 1, chunk 49 | 3.097 | 2.010 | 10.450 | 11.773 |

Savings vary: total build time falls about 52%, 62%, and 3% in the first three
regions, while the fourth is about 1% slower. Delayed materialization shifts
some work to the candidate stage; reduced decode time does not guarantee faster
total execution. Repeated overlapping windows are not cached in this version.

Reproduction script: `tmp/no_digest_benchmark_20260921/benchmark_windows.py`.
Full JSON and logs: `tmp/no_digest_benchmark_20260921/window_results/report.json`.
The tmp directory is ignored by git; preserve these local artifacts if needed.
Existing production Slurm jobs were not restarted or modified by this test.
