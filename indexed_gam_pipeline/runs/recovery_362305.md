# 32-process recovery: job 362305

Started 2026-09-22 14:52:40; verified RUNNING on meow, 32 CPUs / 420 GiB.

- Repeats only the original 17 unfinished tasks' 549,850 nodes, partitioned into 512 tasks.
- Original 495 completed tasks remain untouched; preserved manifest hashes verified during preparation.
- Same frozen pipeline as job 362303: N filter, EOF fix, int8 tensors, SNV AF >= 0.06, INDEL AF >= 0.08, 8 GiB GAM cache per builder. Only concurrency changes from 48 to 32.
- All 512 recovery tasks restart in a fresh output directory; failed job 362303's complete and partial recovery outputs are excluded from the combined catalog.
- Reuses discovery and GBZ index. Each small task exits and releases its process memory.
- Output: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/recovery17_p32_512_nfilter_eoffix_20260922`.
- Finalization publishes `combined_outputs.json` including only original 495 task outputs and this recovery. Original results retain their prior N-filter behavior.
- Preparation: `indexed_gam_pipeline/runs/prepare_recovery32.py`.
