# HG008 PacBio Revio 116x full tensor run

Input: `/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x.AF-HPRC.sorted.gam`

Output: `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125`

User-confirmed format: candidate-v4, seven channels, `int32`. Full NPY shards
contain `(2048, 7, 200, 100)`; only the last shard can be smaller. No exploratory
read/node/tensor limits apply to the production discovery/build. Candidate
filters: MAPQ > 10, ALT >= 3, AF >= .05, allele BQ >= 10, SNP/INS/DEL <= 50 bp.
Discovery uses MAPQ > 5 and imperfect-mapping fraction > .05.

`full_run.py` runs these timed stages:

1. Verify selected `/scratch/jshen/bin/vg_v1.77.0` version.
2. Build a complete node-sequence SQLite index from the matching d9 GFA (all S
   records, including alternate nodes), recording source identity.
3. Build up to eight HG008 preflight tensors with debug rows; independently
   audit these against GAM records, graph bases and GBWT counts. These are
   temporary verification examples, not production shards.
4. Discover candidate nodes from the entire 165 GB GAM.
5. Build every passing candidate with 2,048 tensors per shard; no `--max-tensors`.
6. Validate every shard header, shape, dtype, byte length and summary index/count.
7. Move completed outputs from `.building_tensors/` into the requested directory.

The job fails explicitly if graph identity, GAM edits, source-file stability,
preflight auditing, or final structural validation fails. It does not silently
skip batches. Failed jobs retain logs and intermediate files for diagnosis;
automatic restart/resume is not implemented.

`run/status.json` records each stage's state, command, wall time and timestamps.
Subprocess stages also have `run/<stage>.resources.txt` from `/usr/bin/time -v`
and `run/<stage>.log`. The graph-index stage logs progress to the Slurm output.
The build saves `batch_timing.ndjson` with per-batch elapsed times, cumulative
stage timing, node ranges, fetch metrics and tensor counts. Manifests are updated
after each completed shard. The final manifest contains total build and stage
times, and `validation_report.json` reconciles all NPY/summary records.

A source snapshot, helper binary and configuration are retained under `run/`.
The 1 GiB GAM cache is additional to tensor, long-read and GBZ memory. Production
uses 32-node batches and a 20,000-alignment batch limit; an oversized batch fails
for investigation rather than truncating evidence. Shards use about 1.15 GB each.
No PNGs or truth labels are generated for the full dataset. Seven-channel-aware
model architecture and preprocessing are required; older five/six-channel
checkpoints do not become compatible merely by using sharded NPY files.
