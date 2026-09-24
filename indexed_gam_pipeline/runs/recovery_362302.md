# Recovery of job 362293

Submitted job **362302**, initially RUNNING on **meow**, 48 CPUs / 420 GiB.

Recovery directory:
`/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/recovery_362293_p48_512_eoffix_20260922`

- Original unfinished tasks: 38, 175, 258, 263, 270, 271, 276, 285, 393, 440, 441, 445, 486, 489, 490, 498, 511.
- Exactly 549,850 nodes repartitioned into 512 disjoint tasks (1,073–1,074 nodes each), with 48 concurrent builders.
- Original 495 completed tasks retained in place. Their SNV/INDEL manifest hashes are recorded in `recovery_inventory.json`.
- Only the 17 unfinished tasks' python/SNV/INDEL directories (51 total) were moved out of the original output trees into `excluded_partial_backup/` under the recovery directory. This backup is recoverable and must never be included in result discovery.
- Source snapshot copied from the failed run, with only `indexed_gam_pipeline/gam_reader.py` replaced by the EOF-fixed repository version. Cache remains 8 GiB per process; tensor and filtering settings are unchanged.
- `run.sh` runs the existing build/validation queue, then `recovery.py finalize` publishes `combined_outputs.json` containing only the 495 preserved and 512 recovered output directories. Archives are excluded. The catalog is written only after all recovery tasks pass validation and preserved manifest hashes still match.
- Monitor Slurm job 362302 and recovery `python/status.json`, `status.json`, `memory.ndjson`, and `slurm-362302.err`.

Preparation helper: `indexed_gam_pipeline/runs/recover_362293.py`.
