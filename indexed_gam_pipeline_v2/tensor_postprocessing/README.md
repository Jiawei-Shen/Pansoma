# tensor_postprocessing

Steps that run **after** every tensor task of a run has finished. Nothing here imports the
tensor builder (`candidates.py`, `build.py`), so the same code serves every tensor format
(v5, v6, ...); the format is read from each directory's `manifest.json`.

```
graph (once)   ref-path-scan ─▶ ref-path-check ─▶ chr-index
run            merge ─▶ label
```

| Command | Module | What it does |
|---|---|---|
| `ref-path-scan` | `reference_path.py` | one pass over the GFA (awk prefilter): node lengths, GRCh38 walk coordinates of every node, node range of every walk of every sample |
| `ref-path-check` | `reference_path.py` | GFA lengths vs the GBZ graph index, reference node sequences vs the GRCh38 FASTA (100k random nodes each) |
| `chr-index` | `chr_index.py` | node ID → chromosome block table (`*.chr_node_ranges.tsv/.json`) |
| `merge` | `merge_shards.py` | `task_*/shard_*` (2,048 each) → `<chrom>_shard_*` (32,768 each), byte-verified |
| `label` | `truth_labels.py` | germline/somatic VCF → node candidate keys → one label per tensor, recall report |

Run from the repository root: `python -m indexed_gam_pipeline_v2.tensor_postprocessing <command> --help`.

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
(`orchestrate prepare`) and before the first batch (`run.py build`).

## In a whole-genome run

`orchestrate prepare` freezes the options below into `config.json`; `orchestrate run` calls
`orchestrate finalize` (merge, then labels) when every task validated. `finalize` can also be
run on its own and skips finished steps.

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

Records keep source order, so shards stay node-sorted. `grch38` is null when the node has no
unique GRCh38 visit (alternative branches); `pos0` is the 0-based REF start (INS: the boundary
between `pos0 - 1` and `pos0`), bases on the contig's forward strand, no anchor base.

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
   on GRCh38 (shifted left and right), because vg places it at one end of the repeat in read
   orientation (v5) and VCFs are left-aligned.
3. **Node keys**: each placement becomes the builder's candidate identity
   `node:start:KIND:REF>ALT` in node-forward coordinates (reverse-oriented nodes flip the
   position and reverse-complement the bases; an insertion at a node junction gets a key on
   both nodes; a deletion must fit inside one node; only nodes with one GRCh38 visit).
4. **Labels**, from the tensor's representative allele (`candidate_id`):

| value | name | rule |
|---:|---|---|
| 1 | somatic | representative allele is a somatic truth allele with FILTER PASS/`.` (inside or outside the somatic BED) |
| 2 | germline | representative allele is a germline truth allele with FILTER PASS/`.` (inside or outside the germline BED) |
| 0 | non | unique GRCh38 node, inside somatic BED ∩ germline BED, no truth allele at the site, none within 10 bp |
| −1 | ignore | anything else; `reason` says why (`outside_confident_region`, `not_on_unique_grch38_node`, `near_truth_allele_mismatch`, `somatic_truth_filtered`, `germline_truth_filtered`, `truth_matches_non_representative_allele`) |

Outputs, next to the merged shards and in summary order:
`<chrom>_shard_NNNNN_labels.npy` (int8), `<chrom>_labels.ndjson` (label, reason, matched truth
alleles with GT/filter, GRCh38 coordinates), `labels.manifest.json` (truth/BED SHA-256, counts
per chromosome and reason). Truth tables `<set>.graph.tsv` (every truth allele with its keys)
go to `--truth-dir`; `<set>.recall.tsv` and `truth_recall.json` (per truth allele: tensor as
representative / other allele, filtered with reasons, no unique-node key, or no candidate)
to `--recall-dir`.

## HG008 PacBio v5 commands

```bash
G=/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/graph_index
python -m indexed_gam_pipeline_v2.tensor_postprocessing ref-path-scan \
    --gfa /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa --output $G/hprc-v1.1-mc-grch38.d9.grch38_path
python -m indexed_gam_pipeline_v2.tensor_postprocessing ref-path-check --path $G/hprc-v1.1-mc-grch38.d9.grch38_path \
    --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite \
    --fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
python -m indexed_gam_pipeline_v2.tensor_postprocessing chr-index \
    --components-dir /scratch/jshen/data/AF-Filtered_VG_Indexes/chr_component_vs_GRCh38_summary \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path --output $G/hprc-v1.1-mc-grch38.d9.chr_node_ranges \
    --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite
python -m indexed_gam_pipeline_v2.tensor_postprocessing merge --root <run root> \
    --chr-index $G/hprc-v1.1-mc-grch38.d9.chr_node_ranges.tsv --shard-size 32768 \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path [--keep-sources]
python -m indexed_gam_pipeline_v2.tensor_postprocessing label --tensors <tensors dir> \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path \
    --fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --somatic-vcf  /scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz \
    --somatic-bed  /scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed \
    --germline-vcf /scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz \
    --germline-bed /scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed \
    --truth-dir /scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/truth
```

Measured on HG008 (2026-09-23): ref-path-scan 7.4 min / 4.2 GB (49,092,514 GRCh38 nodes, none
visited twice); ref-path-check 0 length and 0 sequence mismatches; chr-index 29 s.
