# Latest optimized recovery: job 362967

Started 2026-09-22 19:38:10; verified RUNNING on meow, 32 CPUs / 420 GiB, with no startup stderr.

- Only 74 unfinished recovery tasks / 79,438 nodes are rebuilt.
- Preserves 495 original completed tasks plus 438 completed recovery tasks (391 from 362305, 47 from 362961). All preserved manifest hashes were verified. Previous incomplete outputs are excluded.
- Frozen source includes both tested optimizations: target-node-only candidate observations and immediate release of each decoded protobuf alignment. Prior end-of-batch release and N/EOF fixes remain.
- Candidate/build source hashes match the optimized benchmark snapshot from `target_observations_release_test_20260922.md`.
- Configuration: full decoding, 32 concurrent builders, 8 GiB GAM cache per builder, 512 target nodes per batch, int8 tensors, SNV AF >= 0.06, INDEL AF >= 0.08. Node discovery and graph index are reused.
- Output root: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/recovery_remaining74_p32_full_targetobs_release_20260922`.
- After successful completion, `finalize.py` publishes a verified `combined_outputs.json` containing the preserved 933 task outputs and these 74 new tasks, with separate SNV/INDEL entries. Original 495 retain their earlier N-filter behavior.
- Preparation/finalization: `indexed_gam_pipeline/runs/recover_latest_362961.py`.
