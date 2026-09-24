# Target-observation and protobuf-release tests

Only two optimizations were implemented:

1. `decode_alignment(..., target_nodes=...)` retains candidate observations only on the current batch's target nodes. Full columns, visits, read hashes and unsupported-event auditing remain unchanged. Omitting the argument retains the previous API behavior.
2. The builder clears each original protobuf alignment's list entry and loop reference immediately after decoding that record. Graph context has already been loaded. The existing end-of-batch reference release remains in place.

No compact-array representation, window decoding, filtering threshold changes, production job replacement, or resource-configuration changes were made.

## Regression tests

`python -m unittest discover -s indexed_gam_pipeline/tests`: 91 tests, 2 skipped, no failures. New tests cover empty/subset/all target-node sets, repeated mappings, both orientations, SNP/INS/DEL, unsupported events outside target nodes, identical full context, support/tensors/debug metadata, and byte-identical multi-batch outputs with unrestricted versus restricted observations.

## Complete heavy-batch decode experiment

Same first batch of prior recovery task 22: 512 target nodes, 3441 selected complete reads, 73,918,852 read bases, full decoding and 8 GiB GAM cache. Baseline diagnostic job 362962 and optimized diagnostic phase of job 362965 ran on tsingtao.

| Metric | Baseline | Optimized |
|---|---:|---:|
| Process high-water RSS | 25.368 GiB | 20.415 GiB |
| Decode phase wall time | 188.110 s | 155.580 s |
| Full columns | 96,533,414 | 96,533,414 |
| Mapping visits | 15,110,694 | 15,110,694 |
| Retained observations | 7,309,616 | 336,252 |
| Unique target candidates | 13,265 | 13,265 |

Peak memory decreased 19.5%; measured decode time decreased 17.3%. These are single full-batch diagnostic runs, not repeated whole-pipeline measurements. Input fetch took 54.9 s versus 100.9 s, so total diagnostic wall time must not be used to attribute CPU speed to the code change. Diagnostics used GC/trim only after the measured decode stage; production code does not add per-read GC/trim.

An initial full-batch end-to-end ABBA benchmark (362964) was stopped while the baseline was performing candidate counting because even an 8-tensor output cap did not bound that work sufficiently. No complete full-batch end-to-end comparison is claimed.

## Actual-data bounded end-to-end comparison

Fixed first 128 MAPQ-passing reads from the same real batch, saved as a miniature indexed GAM; original graph index and 512-node target list retained. Each run uses full decoding, 8 GiB cache, one worker, original ALT/BQ/AF thresholds, and an 8-tensor cap. All runs produce 4 SNV and 4 INDEL tensors.

Jobs 362965 and 362966 each run baseline/optimized/optimized/baseline on tsingtao. Tensor shard, variant summary, filtered-candidate and unsupported-event files have identical SHA256 hashes across all versions within each sequence.

| Sequence | Version | Wall seconds | Peak RSS KiB |
|---|---|---:|---:|
| First 1 | Baseline | 14.521 | 1,023,280 |
| First 2 | Optimized | 7.834 | 916,024 |
| First 3 | Optimized | 9.240 | 916,064 |
| First 4 | Baseline | 8.638 | 1,023,352 |
| Follow-up 1 | Baseline | 11.447 | 1,022,656 |
| Follow-up 2 | Optimized | 7.436 | 916,008 |
| Follow-up 3 | Optimized | 7.384 | 915,948 |
| Follow-up 4 | Baseline | 8.667 | 1,022,736 |

Cold-start/I/O and CPU variability affect first baseline runs; do not advertise a speedup from the inflated first baseline alone. In the follow-up, comparable warm fetch times (~0.875–0.877 s) accompany optimized wall times 7.384–7.436 s versus baseline 8.667 s, approximately 14–15% lower elapsed time. No consistent slowdown was observed; these small samples do not guarantee the same speedup or memory safety for 32 concurrent production builders.

Artifacts (local, ignored tmp): `tmp/target_release_benchmark/` contains frozen baseline/optimized source, sample GAM, benchmark scripts, raw logs, `sample_results.json`, `warm_results.json`, and `full_optimized_profile.log`. Original complete-batch baseline: `tmp/memory_probe_362962.log`.
