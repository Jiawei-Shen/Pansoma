# tensor_postprocessing

Steps that run **after** every tensor task of a run has finished. Nothing here imports the
tensor builder (`candidates.py`, `build.py`), so the same code serves every tensor format
(v5, v6, ...); the format is read from each directory's `manifest.json`.

```
graph (once, tools/graph_prep.py)   ref-path-scan ─▶ ref-path-check ─▶ chr-index
run (this package)                  merge ─▶ label
```

| Command | Module | What it does |
|---|---|---|
| `ref-path-scan` | `tools/graph_prep.py` | one pass over the GFA (awk prefilter): node lengths, GRCh38 walk coordinates of every node, node range of every walk of every sample |
| `ref-path-check` | `tools/graph_prep.py` | GFA lengths vs the GBZ graph index, reference node sequences vs the GRCh38 FASTA (100k random nodes each) |
| `chr-index` | `tools/graph_prep.py` | node ID → chromosome block table (`*.chr_node_ranges.tsv/.json`) |
| `merge` | `merge_shards.py` | `task_*/shard_*` (2,048 each) → `<chrom>_shard_*` (32,768 each), byte-verified |
| `label` | `truth_labels.py` | germline/somatic VCF → node candidate keys → one label per tensor, recall report |

Run from the repository root:

```bash
python -m indexed_gam_pipeline_v3.tensor_postprocessing merge|label --help      # run time, frozen with every run
python -m indexed_gam_pipeline_v3.tools.graph_prep ref-path-scan|ref-path-check|chr-index --help   # once per graph
```

The graph-prep code is offline tooling: it lives in `tools/`, which `orchestrate prepare` does not
freeze. The run-time readers of its outputs stay here: `reference_path.ReferencePath`
(gfa-reference-path-v1 directory) and `chr_index.ChrIndex` / `select_nodes` (chr-node-ranges-v1
TSV plus JSON with `tsv_sha256`). The formats are unchanged, so directories and tables made by
earlier versions can be used as they are.

## Chromosome of a node (`chr_index.py`)

Minigraph-Cactus builds one graph per reference chromosome and then numbers the nodes
chromosome by chromosome, so each chromosome's graph is one contiguous node-ID interval.
For HPRC v1.1 d9 the order is the contig-name order: chr1, chr10–chr19, chr2, chr20–chr22,
chr3–chr9, chrEBV, chrM, unplaced, chrX, chrY.

* **chr1–22** blocks are exactly the node sets of `vg chunk -C -p GRCh38#0#chrN` (the whole
  connected component, including every node that is not on GRCh38 — insertions and
  alternative branches, 17–22 % of each component), from `scripts/build_chr_node_filters.sh`.
  Each list must be one gap-free interval.
* **Non-autosomal groups** (chrEBV, chrM, chrX, chrY; `*_random`, `chrUn_*` and HPRC contigs
  assigned to no chromosome → `unplaced`) take the remaining IDs, split at each group's
  smallest GRCh38 node.
* **Walk check:** every walk of every sample in the GFA (one walk = one assembly contig)
  must stay inside one block. A walk leaving an autosome block is an error.

HPRC v1.1 d9 result: 27 blocks cover all 60,118,570 nodes; 8,688 walks of 46 samples,
0 crossing a block.

The same table drives the builder's `--chromosomes all|autosome|chr1,chr2,...`
(`select_nodes`): target nodes outside the chosen blocks are dropped before partitioning
(`orchestrate prepare`) and before the first batch (`run build`).

## In a whole-genome run

`orchestrate prepare` freezes the options below into `config.json`; `orchestrate run` calls
`orchestrate finalize` (merge, then labels) when every task validated.
`finalize` can also be run on its own and skips finished steps.

```
--merge-shard-size 32768      0 = keep the task layout (no merge, no labels)
--keep-sources                keep task_* after the verified merge (default: delete them)
--chr-index TSV               needed for the merge and for --chromosomes other than all
--reference-path DIR          GRCh38 coordinates in the merged summaries; needed for labels
--somatic-vcf --somatic-bed --germline-vcf --germline-bed --reference-fasta --truth-dir   (labels; all or none)
```

## Merge (`merge_shards.py`)

Input: `<root>/outputs.json` of a completed `orchestrate run`. Output per kind (SNV, INDEL):

```
<tensors>/<kind>/<chrom>_shard_NNNNN_data.npy     chr1..chr22, --shard-size tensors each, last one shorter
<tensors>/<kind>/<chrom>_variant_summary.ndjson   source record, plus chrom, shard_file, new shard_index /
                                                  index_within_shard, source_task / source_shard_index /
                                                  source_index_within_shard, grch38 {chrom, pos0, ref, alt, node_reverse}
<tensors>/<kind>/filtered_candidates.ndjson       all task audit streams, task order
<tensors>/<kind>/unsupported_events.ndjson
<tensors>/<kind>/manifest.json, validation_report.json
<tensors>/non_autosomal/<kind>/...                chrX, chrY, chrM, chrEBV, unplaced (not training data)
<root>/batch_timing.ndjson                        every task's batch timings
```

Records keep source order: task index, then shard, then index within the shard. The tasks of one
run cover contiguous node ranges in order, so its merged shards are node-sorted; a merge over tasks
of two node lists (e.g. a run followed by extra tasks for nodes added later) holds the first list's
tensors followed by the second's. Use `node_id` (or `grch38`) when position order matters.
Merges of runs made with v2 (PacBio v6: main tasks, then three supplement rounds) are ordered the
same way, round after round.
`grch38` is null when the node has no unique GRCh38 visit (alternative branches); `pos0` is the
0-based REF start (INS: the boundary between `pos0 - 1` and `pos0`), bases on the contig's forward
strand, no anchor base.

Everything is written to a hidden `.merging/` directory and published only when:
every output shard's data SHA-256 (re-read from disk) equals the hash of the source tensors
copied into it; every summary record equals its source record apart from the position
fields; totals match `outputs.json`; shard headers, lengths and file sizes agree; and random
tensors reloaded from the source task shards equal the merged ones. The task directories are
deleted afterwards unless `--keep-sources` is given; `outputs.json` gets a `merge` section
(original kept as `outputs.pre_merge.json`) and `status.json` becomes `finalized` once the
sources are gone. A second merge of the same run is refused.

## Labels (`truth_labels.py`)

1. **Truth alleles**: every ALT on chr1–22, trimmed of shared prefix/suffix bases; MNVs become
   one SNV per differing base; complex replacements are kept only as "nearby truth".
2. **Equivalent positions**: an indel in a repeat is enumerated at every equivalent placement
   on GRCh38 (shifted left and right). VCFs are left-aligned, and the builder left-normalizes
   on the node's forward strand, which on a reverse-oriented node is the other end of the repeat.
3. **Node keys**: each placement becomes the builder's candidate identity
   `node:start:KIND:REF>ALT` in node-forward coordinates (reverse-oriented nodes flip the
   position and reverse-complement the bases; an insertion at a node junction gets a key on
   both nodes; only nodes with one GRCh38 visit). A deletion over several consecutive reference
   nodes of one orientation gets the builder's multi-node key: `start` on the forward-first
   node, REF over all nodes, the further nodes as an `@n2+n3` suffix
   (`node:start:DEL:REF>@n2+n3`). A span over nodes of mixed orientation has no key (the builder
   only joins a deletion over mappings of one orientation). The `truth_labels.py` module
   docstring still says "a deletion must fit inside one node"; the code is the rule above.
4. **Labels**, from the tensor's representative allele (`candidate_id`):

| value | name | rule |
|---:|---|---|
| 1 | somatic | representative allele is a somatic truth allele with FILTER PASS/`.` (inside or outside the somatic BED) |
| 2 | germline | representative allele is a germline truth allele with FILTER PASS/`.` (inside or outside the germline BED) |
| 0 | non | unique GRCh38 node, inside somatic BED ∩ germline BED, no truth allele at the site, none within 10 bp |
| −1 | ignore | anything else; `reason` says why (`outside_confident_region`, `not_on_unique_grch38_node`, `near_truth_allele_mismatch`, `somatic_truth_filtered`, `germline_truth_filtered`, `truth_matches_non_representative_allele`) |

"Within 10 bp" (`NEAR_BP`) is measured against each truth allele's whole span of equivalent
placements: an insertion in a repeat spans from its leftmost to its rightmost equivalent boundary,
because reads may write it anywhere in the repeat. Any truth set and any FILTER count as nearby
truth.

Outputs, next to the merged shards and in summary order:
`<chrom>_shard_NNNNN_labels.npy` (int8), `<chrom>_labels.ndjson` (label, reason, matched truth
alleles with GT/filter, GRCh38 coordinates), `labels.manifest.json` (truth/BED SHA-256, counts
per chromosome and reason). Truth tables `<set>.graph.tsv` (every truth allele with its keys)
go to `--truth-dir`; `<set>.recall.tsv` and `truth_recall.json` (per truth allele: tensor as
representative / other allele, filtered with reasons, no unique-node key, or no candidate)
to `--recall-dir` (default: the tensor directory).

## HG008 PacBio commands

```bash
G=/scratch/jshen/data/pansoma_v2_tensors/graph_index
P="python -m indexed_gam_pipeline_v3"
# once per graph (the existing files under $G were made by the same code and are still valid)
$P.tools.graph_prep ref-path-scan \
    --gfa /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa --output $G/hprc-v1.1-mc-grch38.d9.grch38_path
$P.tools.graph_prep ref-path-check --path $G/hprc-v1.1-mc-grch38.d9.grch38_path \
    --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite \
    --fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
$P.tools.graph_prep chr-index \
    --components-dir /scratch/jshen/data/AF-Filtered_VG_Indexes/chr_component_vs_GRCh38_summary \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path --output $G/hprc-v1.1-mc-grch38.d9.chr_node_ranges \
    --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite
# manual re-merge / re-label of a root prepared by indexed_gam_pipeline_v3 (orchestrate finalize does both)
$P.tensor_postprocessing merge --root <run root> \
    --chr-index $G/hprc-v1.1-mc-grch38.d9.chr_node_ranges.tsv --shard-size 32768 \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path [--keep-sources]
$P.tensor_postprocessing label --tensors <tensors dir> \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path \
    --fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --somatic-vcf  /scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz \
    --somatic-bed  /scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed \
    --germline-vcf /scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz \
    --germline-bed /scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed \
    --truth-dir /scratch/jshen/data/pansoma_v2_tensors/truth
```

Measured on HG008 (2026-09-23): ref-path-scan 7.4 min / 4.2 GB (49,092,514 GRCh38 nodes, none
visited twice); ref-path-check 0 length and 0 sequence mismatches; chr-index 29 s.

## History

This directory was copied from `indexed_gam_pipeline_v2/tensor_postprocessing` at git commit
`d0d25d6` (`git archive`; `git show d0d25d6:indexed_gam_pipeline_v2/tensor_postprocessing/<file>`
shows the originals). v2 was retired on 2026-09-26 and removed from the repository in `9d61016`;
this copy is now the only one. Changes since the copy:

* `chr_index.py`, `reference_path.py`: the graph-prep code moved verbatim to `tools/graph_prep.py`
  (reference_path `SEPARATORS`, `AWK`, `parse_walk`, `scan`, `check`; chr_index `VERSION` (there
  `CHR_INDEX_VERSION`), `NAMED_GROUPS`, `UNPLACED`, `group_of`, `component_block`, `build`); the
  unused `SELECTIONS`, `ChrIndex.is_autosome` and `ReferencePath.unique` were dropped.
* `__main__.py`: `merge` and `label` only, no module-search-path block.
* `merge_shards.py` (`8367f73`, 2026-09-25): the parallel copy: one job per (kind, chromosome)
  group, rows read and shards written with plain sequential file I/O instead of memory maps, audit
  streams copied in slices to their offsets. Same bytes (goldens, a real 20-task merge); HG008
  Illumina chr22 tasks: copy 227 → 101 s with 8 workers.
* `truth_labels.py` (`ba1dec2`, 2026-09-26, from v2 `5e82150`): truth-labels-v2 — labels 1 and 2
  need a PASS truth allele, and the BEDs only bound the confident region for label 0.

A change to label rules or merged bytes changes the goldens (`tests/golden_hashes.json`); record
them again from the committed change (main README, section 9).
