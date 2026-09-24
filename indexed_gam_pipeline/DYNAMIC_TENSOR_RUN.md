# Dynamic tensor-only rebuild

`dynamic_tensor_run` reuses a prior unified run whose graph build, graph audit,
input verification and preflight completed. Preparation freezes the current working-tree pipeline and records copied file
hashes, including the int8 storage encoding. The prior preflight is provenance
for the reused inputs/index; it is not a preflight of the new source snapshot. The graph index remains in the prior run directory:
do not delete that directory or its index.

```bash
python -m indexed_gam_pipeline.dynamic_tensor_run prepare \
  --prior /path/to/prior/run --root /path/to/new/run \
  --tasks 512 --processes 32 --cache-mb 8192 \
  --snv-min-af 0.06 --indel-min-af 0.08 \
  --max-node-reads 800 --candidate-unit site
sbatch --cpus-per-task=32 --mem=420G --time=14-00:00:00 \
  --output=/path/to/new/run/slurm-%j.out \
  --error=/path/to/new/run/slurm-%j.err /path/to/new/run/run.sh
```

The queue dispatches contiguous node tasks as slots become available. Each task
uses a fresh builder process, then validates its shards in the task supervisor.
The builder exits before validation; all task processes exit before the slot is
reused. This releases caches and allocator-held memory between tasks, but does
not impose a per-task peak memory limit. Large individual nodes can still be slow.
The configured process count limits concurrent tasks, not the number of lightweight
supervisor processes. Each builder uses one candidate worker and single-threaded
numerical libraries through run.sh.

`python/status.json` records pending/running/completed tasks, PIDs and timings.
`python/task_NNNN/` contains each task's independent shards and validation report.
`memory.ndjson` samples the controller's process tree every 30 seconds. Any task
failure stops the queue and terminates remaining task process groups. Automatic
resume is not implemented; an existing output directory is rejected to prevent
accidental overwrites or duplicate append operations. Shards are not merged.

Before deleting an old run's tensor outputs, cancel its job and confirm that it
has stopped. Keep the graph index, frozen source, configuration and audit reports.

With both AF options, each of the 512 node tasks uses one disposable builder.
Each batch fetches GAM records and decodes edits once; both candidate types share
those decoded reads. SNP uses AF >=0.06 and INS/DEL use AF >=0.08 before tensor
construction. Outputs are `SNV/task_NNNN/` and `INDEL/task_NNNN/`, each with its own
bounded shard buffer, manifest, filtered-candidate audit and summary. ALT count
>=3 and the remaining candidate filters are unchanged.

`python/task_NNNN/` holds the shared coordinator manifest, batch timings and audit,
with no mixed tensor shards. Shared read/decode timing and cache statistics are
referenced in the type manifests; do not add those shared values across types.
At most 32 builders run concurrently. The builder exits before both output
validations and the slot is reused. Full decoding remains selected, matching the
previous production job. Both outputs use int8 with missing quality -1 and
channel 7 `min(count // 4, 127)`.

Direct CLI builds use `--variant-type all --snv-output /new/SNV --indel-output
/new/INDEL --snv-min-af .06 --indel-min-af .08`, with `--output /new/shared` for
the shared coordinator. The three directories must be new, distinct and
non-nested. `--max-tensors`, if used for a test, caps the combined output count.
