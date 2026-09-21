# HG008 independent partition and VG reader tests — 2026-09-21

## Scope

Slurm job **362266**, source **4e2a943**, completed 0:0 in 4m08s. Exactly two
new treatments were run, as requested; there was no combined VG+partition test
and no new baseline execution. Both treatments enable only the early ALT upper
bound filter. The optional candidate node-index cache is **disabled** and each
builder uses **one candidate worker**. Source discovery/indexing are reused.
The full job **362252 remains STOPPED** and was not resumed.

To compare directly to the previously recorded early-ALT-only run, the test uses
its exact **1,024 nodes**, two batches of 512, `(7,200,101)` seven-channel int32
tensors, 2,048 tensor shard capacity, and **debug row metadata enabled**. This
last setting is for comparison/auditing; the intended full run does not require
it. The baseline is `candidate_optimization_20260921T182037Z/early_alt_only`.
All modes produce 182 tensors: 89 DEL, 17 INS, 76 SNP.

## Treatments and results

| Mode | Build subprocess wall | Additional merge + validation | Sampled process-tree peak PSS |
|---|---:|---:|---:|
| Previous Python single-process, early ALT only | 131.54 s | — | 6.15 GiB |
| Python, two independent contiguous partitions | 141.46 s | 19.52 s | 7.88 GiB |
| VG reader, single builder | 63.39 s | — | 2.35 GiB |

Wall time includes process startup/shutdown and all build stages through NPY and
summary writing. For the partition treatment it is elapsed time while both
builders run, not their summed time. The 19.52 s includes repacking and validating
the merged result, not just copying arrays. Per-part validation and SHA checking
are outside build time. The overall Slurm duration also includes preparation,
validation and controller work. The old baseline was measured in an earlier job;
cluster load and storage/cache state were not controlled, so do not attribute
every timing difference to implementation alone.

Memory is externally sampled every five seconds; the maxima can miss short
peaks. The sampler excludes its own benchmark controller and merge phase, but
includes VG child processes. Its cumulative sampling wall cost was 2.18 seconds
for partitioned Python and 0.36 seconds for VG; the baseline's sampler was 1.67
seconds (see original JSON for precise value). Sampler time is concurrent and
cannot simply be subtracted from build time.

### Python contiguous partitions

- Builder 0 owns nodes 128–2430 (512 nodes); builder 1 owns 2431–4960 (512).
- Each builder independently fetches reads, loads graph/GBZ context, decodes,
  filters candidates, builds tensors, and writes its own output directory.
- GAM cache budget is **4 GiB per builder, 8 GiB total**. These limits are cache
  accounting limits, not whole-process RAM limits. Separate GBWT SQLite copies
  avoid concurrent writes to the same database.
- Builder 0 took 141.02 s internally, builder 1 took 120.99 s. Their GAM fetch
  times were 105.35 s and 105.27 s; both initialized indexes for the same 66 GAM
  groups. Output was 68 and 114 tensors respectively.
- Parts are merged in node/candidate order. Tensor rows are unchanged; only
  summary shard/index locations are renumbered. The resulting final shard has
  182 tensors and is a valid partial final shard under the 2,048 capacity.

No speed improvement appeared in this small workload. The two intervals reuse
the same groups, which benefits a persistent single reader but duplicates cold
loading for independent readers. This is not equivalent to testing distant
halves of the entire 16.56-million-node list, and does not disprove the utility
of partitioning the full dataset. A 4-GiB-per-process cache also differs from the
old single 8-GiB cache, deliberately preserving the combined cache budget.

### VG single builder

| Stage | Previous early-ALT-only Python | VG builder |
|---|---:|---:|
| GAM fetch, filtering and context collection | 83.95 s | 18.57 s |
| Graph sequences | 0.89 s | 0.71 s |
| GBWT occurrence | 4.57 s | 3.95 s |
| Decode edits | 19.79 s | 19.59 s |
| Candidate handling and shard/summary writing | 20.18 s | 19.99 s |
| Builder internal total | 130.82 s | 63.06 s |

The two VG commands took 17.96 seconds combined and returned 16,145,610 bytes of
compressed output. Unlike the earlier reader-only test, this test includes all
subsequent tensor stages. VG speedup over the saved baseline is about **2.08x**
for this workload, primarily from GAM retrieval. No full-genome speedup or wall
time forecast is justified from these two batches.

## Correctness and code changes

47 regression tests passed. New adapter tests cover exact membership filtering
within a broad range, ordering, duplicate records, empty requests, command
failure, alternate-index rejection and changed inputs.

Both treatments' final NPY and complete `variant_summary.ndjson` SHA-256 values
match the saved baseline **exactly**. That baseline's tensors also matched the
previous independently source-audited 182-tensor set. A new independent source
audit was not rerun here because the complete outputs were byte-identical.

The optional backend is available with:

```text
--gam-reader vg --vg /scratch/jshen/bin/vg_v1.77.0
--early-alt-filter --node-index-cache-nodes 0 --workers 1
```

Python remains the default reader. The VG adapter uses a temporary output GAM,
parses it to complete alignment objects, applies exact target-node membership,
and passes those objects into the unchanged candidate/tensor stages. It checks
input timestamps/sizes and requires the adjacent GAM.gai because VG resolves
that index automatically. Nonzero VG exit is fatal. Transport files are cleaned
up and placed under TMPDIR/system temporary storage; no 8-GiB Python group cache
is maintained for the VG backend. This does not replace the stopped full job's
source snapshot or launch a new full build.

## Records

- [Full commands, per-part timings, validation and checksums](hg008_partition_reader/status.json)
- [Source provenance](hg008_partition_reader/provenance.json)
- [Slurm submission](hg008_partition_reader/submission.json)
- Per-batch timings, validation reports and GNU time resources are alongside.
- All tensors/logs: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/partition_reader_bench_20260921T194232Z`.
