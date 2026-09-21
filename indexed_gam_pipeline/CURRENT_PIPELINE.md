# Current tensor-building pipeline

Recorded against implementation commit **dd8b41f**. The default format is
**`indexed-gam-candidate-v4`**, schema 4. The latest rerun results are in
[samples/colo829t_100/retest.json](samples/colo829t_100/retest.json), and previous
optimization measurements are in [performance.json](samples/colo829t_100/performance.json).

## Overview

```mermaid
flowchart TD
    GAM[Sorted BGZF GAM alignments] --> Discover[Optional exploratory or full node discovery]
    Discover --> Nodes[Target node list]
    Fixed[Existing target list] --> Nodes
    Nodes --> Fetch[Fetch complete alignments using GAI]
    GAM --> Fetch
    GAI[Matching .gam.gai index] --> Fetch
    Fetch --> Graph[Load sequences for targets and all visited context nodes]
    SQL[Graph sequence SQLite] --> Graph
    GBZ[Matching d9 GBZ and occurrence cache] --> Counts[Verify sequences and obtain distinct GBWT path counts]
    Graph --> Counts
    Counts --> Edits[Decode original GAM edits and identify candidates]
    Edits --> Support[Count ALT / REF / other and filter candidates]
    Support --> Rows[Select up to 200 rows and center the candidate]
    Rows --> Tensor[Fill 7 channels and group rows by visible node path]
    Tensor --> Save[Save NPY shards, NDJSON summaries, manifests and timing]
    Save --> Audit[Independent GAM / graph audit]
    Save --> PNG[Render inspection PNGs]
```

## Selected vg executable

As of 2026-09-19, use **`/scratch/jshen/bin/vg_v1.77.0`** (verified as
`vg version v1.77.0 "Ruby"`) for vg commands, replacing the previous local vg.
On this machine it is a symlink to `/scratch/jshen/bin/vg`.

```bash
source scripts/use_vg.sh
vg version
# Explicit invocation is also supported:
/scratch/jshen/bin/vg_v1.77.0 version
```

Source this after activating Conda so bare `vg` commands in upstream alignment
and graph-extraction scripts resolve to the selected binary. The activation
script verifies that resolution and sets `PANSOMA_VG` to the explicit path.
The v4 tensor builder reads GAM/GAI in Python and queries GBZ through its compiled
helper; neither stage invokes the vg executable. Earlier benchmark records
retain their original provenance.

## 1. Inputs and graph identity

The starting input is an already graph-aligned, sorted BGZF GAM. Mapping and
sorting are upstream prerequisites; truth labeling and model training are
downstream of these candidate tensors.

The COLO829T test uses:

- GAM: `/scratch/jshen/data/COLO829T/illumina/GAM/COLO829T_3M.sorted.gam`.
- GAI: the GAM filename plus `.gai`, unless `--index` overrides it.
- GBZ: `/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz`.
- Graph sequences: the same directory's
  `hprc-v1.1-mc-grch38.d9.GRCh38_CHM13_node_index.sqlite`.
- Target node IDs: one positive numeric ID per line.
- Occurrence cache: a reusable SQLite database made by the v4 pipeline.

The GAM, sequence database, and GBZ must use the same node identities. The GBZ
forward sequence is compared with the sequence database for every requested
context node. Missing nodes or disagreement fail explicitly. For this d9 graph,
node 2753 has 114 bases and 86 GBWT paths; node 51182 is `C` and has 24 paths.

The pipeline reads the existing GAI. If it is absent, `run.py index` can build
one from the sorted GAM. This rerun reused the real-data GAI; index construction
and retrieval against an independent full scan are also covered by synthetic tests.

## 2. Choose target nodes

`run.py discover` scans GAM records sequentially. With its defaults, it uses
alignments with **MAPQ > 5** and counts perfect/imperfect mappings per node.
A mapping is imperfect when an edit changes length or has replacement sequence.
A node is selected when it has at least one imperfect mapping and its imperfect
mapping fraction is **strictly greater than 0.05**. This is a mapping-based
prefilter, not the eventual candidate allele frequency.

`--max-alignments 300000` is an exploratory scan, not a whole-genome discovery
run. Omit that option for a full scan. Optional chromosome and node-whitelist
filters are available. The chromosome filter uses `Alignment.refpos.name`.

The reproducibility benchmark uses the frozen 3,000-node
[benchmark_nodes.txt](samples/colo829t_100/benchmark_nodes.txt), copied from the
historical benchmark input. Do not substitute the first 3,000 nodes of a fresh
discovery run: that is a different selection. The rerun tests discovery separately
and uses the frozen input for byte-for-byte tensor comparisons. For a new analysis,
use the discovery output directly or supply your own target list.

## 3. Retrieve complete alignments in bounded batches

Target IDs are deduplicated and sorted. Batches are bounded by both node count
and ID span; the benchmark uses 32 nodes and a maximum span of 10,000.

The GAI identifies candidate BGZF ranges. Overlapping ranges are merged, GAM
groups are read, and complete alignment records are filtered by actual node
membership. Distinct records with identical names or identical bytes remain
separate observations. A record matching multiple requested nodes is yielded
once within that batch. It can be fetched again in another batch; each candidate
belongs to only its target-node batch.

Build filtering uses **MAPQ > 10** in this test. The optional chromosome filter
is applied here too. All visited nodes are collected, including context nodes
outside the target batch. A maximum-alignment limit fails rather than silently
truncating an oversized batch.

The 64 MiB LRU cache retains serialized groups and node-to-record positions,
avoiding repeated full-group decoding. `--gam-cache-mb 0` disables reuse.
The budget covers retained group data, not total process memory, tensor shards,
the GBZ, or a transient group being decoded.

## 4. Read graph bases and occurrence counts

Graph bases are loaded with bounded SQLite `IN` queries through one read
connection/transaction for the entire build. JSON and GFA sequence inputs also
exist, but SQLite is used here.

For each context node, the occurrence cache supplies its forward sequence and
**number of distinct physical GBWT paths**. On a cache miss, a persistent C++
helper loads the matching GBZ once and queries its GBWT. Because the GBWT is
bidirectional, one forward-node query includes paths traversing either node
orientation. Sequence IDs are converted to physical path IDs and deduplicated.
Repeated visits and reverse-complement copies count once; reference, haplotype,
and generic paths are included. Fragments of one biological haplotype can have
different path IDs, so this is not a unique sample/haplotype frequency.

The cache stores graph path/size/mtime and an initial SHA-256. Source size/mtime
are checked on reuse, and node sequences are checked during lookup. A fully
populated cache avoids loading the GBZ helper. The old GFA W-record count is used
only by the explicitly selected legacy v3 format.

## 5. Decode GAM edits and form candidates

The builder consumes original edits with both read and graph cursors, checks
bounds and replacement sequences, and retains all mappings for row context.
Reverse mappings are converted to forward-node coordinates for candidate identity.

- Substitution edits produce SNP candidates for differing bases.
- Insertions and deletions produce candidates up to 50 bases in this test.
- Adjacent insertions or deletions within one mapping are merged before applying
  the length limit.
- Longer indels and complex replacements are logged as unsupported; complex
  replacement columns can remain as row context.

Candidate identity is `(node, forward position/interval, REF, ALT, event type)`.
There is no repeat left-normalization or equivalence merging across different
nodes. Ordinary GAM match edits are not rescanned to invent SNP candidates.

## 6. Count support and filter

Every eligible overlapping alignment record contributes once per candidate.
Repeated visits in one record resolve to ALT, then REF, then other; the earliest
mapping wins ties.

ALT requires an exact candidate observation meeting the allele-BQ threshold.
SNPs use their base quality, insertions use mean inserted-base quality, and
deletions use the minimum available neighboring-base quality. REF requires
matching graph bases over the affected interval; insertion REF additionally
requires evidence on both sides of the boundary. Remaining overlaps are `other`.
The allele-BQ threshold is an ALT-observation filter, not a universal row-BQ filter.

Coverage is `ALT + REF + other`; AF is `ALT / coverage`. These counts use all
eligible records before the 200-row cap. The inspection benchmark accepts
ALT count >= 1, AF >= 0.05, allele BQ >= 10, and all supported event types.
The normal CLI default for minimum ALT count is **3**, so it must be explicitly
set to 1 to reproduce this low-support inspection set.

## 7. Construct one tensor per passing candidate

Each tensor is **(7, 200, 100), int32**: channels × alignment rows × columns.
There is no dedicated reference row.

Row selection revision `window-edit-bp-group-uniform-v1` first prepares the
visible window for every eligible record. It counts differing read/reference
bases (substitutions, inserted bases and deleted bases); padding and aligned
absence-of-insertion gaps do not count. Edits outside the rendered window and
already omitted insertion columns do not count. Sort by decreasing edit bp,
then decreasing MAPQ, record SHA-256 and anchor mapping index. Stably group by
the visible, candidate-oriented node path; groups follow their first occurrence
in this sorted list and each group retains the edit-bp order. Only then sample
up to 200 rows, keeping order, at indices `floor(i*(N-1)/(K-1))` for K > 1;
for K = 1 take index `N//2`. There is no ALT-first selection or minimum group
quota. Counts and AF use all eligible records before sampling; duplicate records
remain separate. The seven-channel format stays v4; the row-selection revision
is recorded separately and must match between model training and inference.

The candidate interval
is centered as a complete block. Insertions reserve a shared block and use gap
slots for supported absence; missing evidence stays padding. Context follows
each alignment's own graph traversal, so rows can have different neighboring
nodes. A reverse anchor reverses/complements the row into candidate orientation.

| Channel (1-based) | Meaning |
|---|---|
| 1 | Read base: A=1, C=2, G=3, T=4, N=5, gap=6; padding=0 |
| 2 | Base quality; -1 for missing quality or no read base |
| 3 | Flags: value 1 = difference; 2 = candidate region; 3 = both |
| 4 | Alignment MAPQ, present on aligned gaps too; capped at 32767 |
| 5 | Operation: M=1, X=2, I=3, D=4, complex=5, aligned no-insertion gap=6 |
| 6 | Graph-reference base for that row and column, with the same base encoding |
| 7 | Raw distinct GBWT path count of that column's node |

Deletions keep their mapped node's count. Insertions and aligned insertion-gap
slots use their anchor node's count. Missing coverage and unused rows have zero
occurrence. Candidate-region flags can still mark a central slot without coverage;
use the operation/metadata to distinguish padding from aligned evidence.

Grouping takes place before sampling. All seven channels and row metadata
follow the same selected record order; no second sorting occurs after sampling.

## 8. Save, audit, and render

The builder writes:

- `shard_XXXXX_data.npy`: stacked tensor arrays.
- `variant_summary.ndjson`: candidate identity, coverage, support, AF, candidate
  columns, row groups, shard/index, format, and parameters. `--debug-rows` adds
  source record hashes and per-column node coordinates.
- `filtered_candidates.ndjson` and `unsupported_events.ndjson`: rejection evidence.
- `target_nodes.txt`, `manifest.json`, `run_report.json`: inputs, channel definitions,
  provenance, cache statistics, and stage/total wall times.

`--max-tensors 1000` makes this a bounded test. It does not process all 3,000
requested nodes once 1,000 candidates have been written.

`validate_examples.py` independently queries original GAM records and verifies
coverage, selected records, graph bases/orientation, groups, padding, candidate
flags, and occurrence cells. It requires debug row metadata. Regression tests
cover duplicates, reverse/repeated visits, cache eviction, invalid graphs, and
older supported formats.

The visualizer loads each tensor with its metadata and draws seven panels.
It marks the complete candidate block and labels channel 7 `Distinct GBWT paths`.
Rendering is separate from tensor generation and optional for bulk runs. The
published inspection set remains **100 examples** (70 SNP, 20 INS, 10 DEL).

## Latest post-optimization rerun

Implementation **dd8b41f** passed all **40 regression tests**. The full independent
audit passed for **1,000 tensors**, **1,427 alignment rows**, **104,620 graph-base
columns**, and **109,800 occupied occurrence cells**. All tensor bytes match the
previous output. Parsed summaries match too; the CLI writes the BQ threshold as
`10.0` instead of the earlier programmatic benchmark's `10`.

All **100 inspection PNGs** were regenerated, decoded successfully, and matched
the published images pixel for pixel. The published set remains 100 examples.

This rerun took **120.84 s to build** and **53.20 s to render/check the 100 images**.
The build reused occurrence counts but spent **107.63 s reading graph sequences**.
The earlier **8.16 s** warmed-storage result is not a runtime upper bound; stage
measurements show substantial run-to-run variation in graph reads. Discovery,
auditing, regression tests, and rendering are excluded from build timing.
See [the recorded results](samples/colo829t_100/retest.json) for exact versions,
arguments, input fingerprints, checksums, counters, and image checks.

## Reproduce the bounded test

Run from the repository root, using a new output directory. The current machine's
base Python has NumPy, protobuf, and pysam. Rendering uses the separate interpreter
shown below because the base environment's matplotlib has a NumPy ABI conflict.

```bash
python indexed_gam_pipeline/build_gbz_query.py \
  --deps /scratch/jshen/Github/gbz-tool/dependency --output tmp/gbz_node_counts

# Optional fresh discovery; its output is a different input from the frozen benchmark.
python indexed_gam_pipeline/run.py discover \
  --gam /scratch/jshen/data/COLO829T/illumina/GAM/COLO829T_3M.sorted.gam \
  --output tmp/new_discovery --max-alignments 300000 --min-mapq 5 --node-alt 0.05

python indexed_gam_pipeline/run.py build --format candidate-v4 \
  --gam /scratch/jshen/data/COLO829T/illumina/GAM/COLO829T_3M.sorted.gam \
  --nodes indexed_gam_pipeline/samples/colo829t_100/benchmark_nodes.txt \
  --node-sqlite /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13_node_index.sqlite \
  --gbz /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz \
  --gbz-query tmp/gbz_node_counts --occurrence-cache tmp/hprc_gbwt_counts.sqlite \
  --output tmp/new_build --batch-nodes 32 --max-node-span 10000 --gam-cache-mb 64 \
  --shard-size 1000 --max-tensors 1000 --rows 200 --width 100 \
  --min-mapq 10 --min-af 0.05 --min-variants 1 --min-allele-bq 10 \
  --max-indel-len 50 --variant-type all --debug-rows

python indexed_gam_pipeline/validate_examples.py tmp/new_build \
  --gam /scratch/jshen/data/COLO829T/illumina/GAM/COLO829T_3M.sorted.gam \
  --node-sqlite /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13_node_index.sqlite \
  --gbz /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz \
  --gbz-query tmp/gbz_node_counts --occurrence-cache tmp/hprc_gbwt_counts.sqlite \
  --output tmp/new_build/audit.json

# Preview the first 100 tensors; the published gallery uses its saved source_indices.
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/pansoma_matplotlib \
  /wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python scripts/visualize_tensor.py \
  tmp/new_build/shard_00000_data.npy --all-samples --max-samples 100 \
  --output-dir tmp/new_build/images
```

Implementation entry points: [run.py](run.py), [gam_reader.py](gam_reader.py),
[build_v2.py](build_v2.py) (also implements v4), [candidates.py](candidates.py),
[gbz_counts.py](gbz_counts.py), [gbz_node_counts.cpp](gbz_node_counts.cpp), and
[the visualizer](../src/pangenome_ml_data_generation/tensors/visualization.py).
