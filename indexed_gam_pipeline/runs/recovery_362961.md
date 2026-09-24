# Batch-release recovery: job 362961

Started 2026-09-22 17:45:13; verified RUNNING on meow, 32 CPUs / 420 GiB.

- Continues only 121 unfinished recovery intervals (129,916 nodes).
- Preserves original 495 tasks plus 391 completed tasks from job 362305. All 1,772 preserved SNV/INDEL manifests were verified; canceled partial outputs are excluded.
- User-requested configuration: full decoding, 32 concurrent processes, 8 GiB GAM cache per process, 512 nodes per batch. No window decoding or reduced concurrency/cache/batch size.
- Frozen source differs from 362305 only in `build_v2.py`: after batch processing, clear `by_node`, `reads`, and `alignments`, and remove loop-variable references before the next GAM fetch/graph lookup.
- Previous implementation retained the prior batch's decoded reads during the next batch's fetch and graph lookup. This can increase peak memory; it is not evidence that all batches accumulate indefinitely.
- Full decoding creates one Column per alignment column (104 bytes for the object alone in this environment, plus references, coordinates and other allocations). GAM cache limits do not cover these decoded objects, protobuf alignments, graph data, candidates, or tensor buffers. The largest sampled process in 362305 was 27.51 GiB. Batch release does not guarantee that RSS immediately returns to the OS or that 420 GiB cannot be exceeded.
- Existing candidate and N/EOF behavior unchanged. Tests: 89 tests, 2 skipped, no failures.
- Output root: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/recovery_remaining121_p32_full_cache8g_release_20260922`.
- On successful completion, `finalize.py` verifies preserved manifest hashes and publishes `combined_outputs.json`, including the original 495, previous recovery's 391, and new 121 tasks exactly once. Original 495 retain earlier N-filter behavior.
- Preparation/finalization helper: `indexed_gam_pipeline/runs/recover_memory_362305.py`.
