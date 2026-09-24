# Unified GBZ graph index

The preferred v4 graph input is `--graph-index`. One SQLite file contains every
node's **forward sequence** and **distinct logical GBWT path count**. It replaces
both the GFA-derived `all_graph_nodes.sqlite` sequence source and the per-process
GBZ occurrence cache. It does not alter the seven tensor channels or candidate
filtering/row selection.

## Build once per graph

Run from the repository root, using an environment with the pipeline dependencies.
The native builder additionally needs SQLite development headers/library.

```bash
python indexed_gam_pipeline/build_gbz_query.py \
  --deps /scratch/jshen/Github/gbz-tool/dependency \
  --graph-index --output tmp/gbz_graph_index

/usr/bin/time -v python -m indexed_gam_pipeline.build_graph_index \
  --gbz /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz \
  --builder tmp/gbz_graph_index \
  --output /path/to/run/all_graph_nodes.gbz.sqlite
```

Use a **new filename**. Existing sequence-only databases are not upgraded in
place; neither GFA walk counts nor old occurrence caches are imported.

The C++ builder loads the GBZ once, walks one forward sequence per logical GBWT
path, and deduplicates node visits with a last-seen path token. Reverse visits
are included, while the reverse-complement copy of each entire path is skipped.
Reference paths count too. These are path counts, not unique sample counts.
It then exports all graph nodes, including nodes with zero path count if present.

For compact node IDs, two uint32 arrays cost `8 * (max_node_id + 1)` bytes during
counting (about 459 MiB for 60.1 million IDs), **in addition to the loaded GBZ**.
Sparse IDs use a hash map. The builder does not hold all node sequences or all
extracted paths in extra Python memory. Count-time complexity is proportional to
all path visits; preprocessing can be substantial and is amortized over runs.

Output schema:

```sql
CREATE TABLE nodes (
  node_id INTEGER PRIMARY KEY,
  seq TEXT NOT NULL,
  distinct_path_count INTEGER NOT NULL
);
CREATE TABLE graph_metadata (value TEXT NOT NULL);
```

Metadata records schema/metric, source path/size/mtime/SHA-256, node/path totals,
build completion, and `source_hash_seconds`, `gbz_load_seconds`,
`path_count_seconds`, `sqlite_write_seconds`, `validation_seconds`, and
`total_build_seconds`. Native progress goes to stderr. The wrapper checks source
stability, node count, and SQLite integrity before atomic no-overwrite publication.
A failed build never publishes a usable final filename.

## Tensor building

Keep existing GAM discovery and target-node partitioning. Replace the old
`--node-sqlite`, `--gbz`, `--gbz-query`, `--occurrence-cache` arguments with:

```bash
python indexed_gam_pipeline/run.py build --format candidate-v4 \
  --gam /path/to/input.sorted.gam --nodes /path/to/target_nodes.txt \
  --graph-index /path/to/run/all_graph_nodes.gbz.sqlite \
  --gam-reader python --vg /scratch/jshen/bin/vg_v1.77.0 \
  --early-alt-filter --node-index-cache-nodes 0 --workers 1 \
  --rows 200 --width 101 --shard-size 2048 --output /path/to/new_tensors
```

`--gam-reader vg` also works. Two independent processes receive their own
contiguous candidate-node partitions and output directories, and share this index
through **read-only connections**. No occurrence caches need copying. Each batch
collects all nodes in its retrieved reads, queries sequences and counts together
(in chunks of 900 IDs), decodes edits, applies the ALT upper-bound filter, builds
candidate windows, orders/selects reads, then writes `(2048,7,200,101)` int8 full
shards (the last shard may be shorter).

Missing graph nodes fail explicitly; they do not silently become count zero.
Choose the GBZ matching the GAM's graph: a GAM does not supply an independently
verifiable GBZ fingerprint, and overlapping node IDs alone do not prove a match.
The complete index is self-contained; runtime reads do not rehash/reload GBZ.

`manifest.json`/`batch_timing.ndjson` record `graph_index_seconds` for the combined
sequence/count query. Do not compare it with occurrence time alone: compare it
with the previous **sequence + occurrence** stages. Manifests include index
provenance and query counters.

For debug tensors, `validate_examples.py` accepts `--graph-index` instead of
`--node-sqlite`. Optionally also supply legacy GBZ/helper/cache arguments to audit
counts using the independent search/locate method.

## Orchestrators

- `full_run.py`: set config `unified_graph_index` to the new SQLite path and
  `gbz_index_builder` to the compiled executable. A missing index is built as a
  separately timed stage; an existing complete index is checked against `gbz`.
  GFA and `gbz_query` are unnecessary in this mode.
- `full_reader_comparison.py`: set nested `config.unified_graph_index` and set the
  top-level `graph_index` fingerprint to **that same file**. Prepare a new source
  snapshot/config; retain the input and partition guards. Its two builders share
  the read-only database and skip cache copies.
- Legacy arguments remain supported for existing snapshots/configurations.

## Verification and rollout scope

Native tests compare every node on small dense and sparse GBZ fixtures against
independent GBWT search/locate queries, including revisits and mixed orientations.
End-to-end GAM fixtures compare tensor contents and candidate metadata across old
and unified inputs, then audit the unified tensors against original reads.
Failure tests cover missing nodes, incomplete metadata, read-only enforcement,
conflicting inputs, invalid GBZ, and refusing overwrite.

```bash
GBZ_QUERY=tmp/gbz_node_counts \
GBZ_TOOL=/scratch/jshen/Github/gbz-tool/gbztool \
GBZ_GRAPH_INDEX=tmp/gbz_graph_index \
python -m unittest discover -s indexed_gam_pipeline/tests -q
```

Validation result: **52 tests passed** (2.145 seconds on the fixture suite),
including the configured native integration tests. This time excludes compilation.

The full HPRC unified index has **not** been built or benchmarked as part of these
fixture checks. No full-data speedup is claimed yet. Previously submitted HG008
jobs 362277/362278 keep their frozen legacy pipeline; this change does not migrate
or restart those jobs.

## Full run with 20 Python processes and reused discovery

`prepare_unified_run.py` freezes the source and executables, verifies the original
GAM/GAI/GBZ fingerprints, and streams the existing sorted candidate-node list into
20 disjoint contiguous partitions. It never rescans the GAM for discovery.

```bash
python -m indexed_gam_pipeline.prepare_unified_run \
  --prior-config /path/to/previous/comparison_config.json \
  --output /path/to/new_run --processes 20 \
  --builder tmp/gbz_graph_index --query tmp/gbz_node_counts
sbatch --partition=general --nodes=1 --ntasks=1 --cpus-per-task=20 \
  --mem=400G --time=14-00:00:00 \
  --output=/path/to/new_run/slurm-%j.out \
  --error=/path/to/new_run/slurm-%j.err /path/to/new_run/run.sh
```

The job performs, in order: source/partition verification; all-node unified GBZ
index creation; independent GBWT search/locate count checks at known nodes and
partition starts; 20 concurrent real-GAM preflight builds (32 target nodes per
partition, at most 8 tensors each); full 20-process build and shard validation.
A failed stage prevents the full build. No full-output comparison or shard merge
is performed. Full outputs are `python/part_0/` through `python/part_19/`.
All builders enable the ALT-support upper-bound filter and disable the optional
node-index cache. Each has an 8 GiB GAM cache, one candidate worker, and native
math/OpenMP thread counts set to one to avoid nested oversubscription.

Timing/memory records:

- `status.json`: job stages, wall time, terminal state, sampled peak tree RSS.
- `graph_index_report.json`: offline index phase timings and source fingerprint.
- `memory.ndjson`: every 30 seconds, current stage, per-process RSS/high-water RSS,
  and summed process-tree RSS. This is not PSS; shared pages can be counted twice,
  and short peaks between samples may be missed. No expensive heap/PSS scan.
- `python/status.json`: each partition's command, status, validation and timings.
- `python/part_N/batch_timing.ndjson`: batch/stage times, progress, filtering counts.
- `python/part_N.resources.txt`: GNU time wall/CPU time, peak RSS and I/O counters.
- `job.resources.txt`: GNU time resources for the whole controller/children run.
- Slurm accounting supplies allocation-level memory/time statistics after exit.

Controller tests include an actual 20-process build on a synthetic GAM using one
shared SQLite index; partition tests ensure full coverage without duplicates.
