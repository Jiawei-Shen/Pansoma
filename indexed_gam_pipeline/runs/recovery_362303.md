# N-filter recovery: job 362303

Submitted 2026-09-22; verified RUNNING on meow with 48 CPUs / 420 GiB.

- Replaces job 362302. `scancel 362302` succeeded; accounting records its batch as CANCELLED and parent as FAILED after termination.
- Only the original 17 unfinished task intervals are rebuilt: 549,850 nodes, repartitioned into 512 tasks, 48 concurrent disposable builders, 8 GiB GAM cache per builder.
- Original 495 completed tasks remain untouched. Their manifest hashes were verified during preparation and will be checked again before the combined catalog is published.
- New source snapshot includes the N candidate filter and BGZF EOF reader fix. SNV AF >= 0.06; INDEL AF >= 0.08; int8 storage unchanged.
- N filtering applies only to the recovered nodes. Preserved original results retain their earlier filtering behavior, as requested.
- Reuses node discovery and GBZ index. Each task exits after building, releasing its process memory.
- Smoke tests: 18 candidate tests and 15 reader/pipeline tests passed, including N filtering and EOF offsets. Earlier full suite: 89 tests, 2 skipped, no failures.
- Output root: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/recovery17_p48_512_nfilter_eoffix_20260922`.
- Finalization publishes `combined_outputs.json` with the original 495 tasks plus only this new recovery's outputs. Outputs from canceled 362302 and `excluded_partial_backup/` are excluded.
- An earlier full-dataset preparation directory (`split512_p48_cache8g_nfilter_eoffix_20260922`) was created before the user corrected scope, but no job was submitted for it and it is excluded from this recovery.

Preparation helper: `indexed_gam_pipeline/runs/prepare_nfilter_recovery.py`.
