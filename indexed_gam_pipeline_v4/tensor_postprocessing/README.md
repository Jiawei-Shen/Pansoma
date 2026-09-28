# tensor_postprocessing

Steps that run **after** every tensor task of a run has finished. Nothing here imports the
tensor builder (`candidates.py`, `build.py`): tensor shapes and candidate fields are read from
each directory's `manifest.json` and variant summaries.

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
python -m indexed_gam_pipeline_v4.tensor_postprocessing merge|label --help      # run time, frozen with every run
python -m indexed_gam_pipeline_v4.tools.graph_prep ref-path-scan|ref-path-check|chr-index --help   # once per graph
```

The graph-prep code is offline tooling. It lives in `tools/`, which `orchestrate prepare` does not
freeze. The run-time readers of its outputs stay here:

* `reference_path.ReferencePath` reads the **gfa-reference-path** directory: arrays by node ID
  (`lengths`, `chrom`, `start0`, `reverse`, `visits`), arrays in reference-walk order (`path_nodes`,
  `path_starts`, `path_reverse`), `walks.ndjson`, and `meta.json`, whose key `format` names the
  format. The HPRC directory under `graph_index/` records `version: gfa-reference-path-v1` and is
  read as it is (`common.EARLIER_FORMAT_NAMES`).
* `chr_index.ChrIndex` / `select_nodes` read the **chr-node-ranges** table: a TSV plus a JSON (key
  `format`) whose `tsv_sha256` must match the TSV. The reader checks only that checksum, so the
  HPRC table (JSON `version: chr-node-ranges-v1`) is used as it is.

## Chromosome of a node (`chr_index.py`)

Minigraph-Cactus builds one graph per reference chromosome and then numbers the nodes
chromosome by chromosome, so each chromosome's graph is one contiguous node-ID interval.
For HPRC v1.1 d9 the order is the contig-name order: chr1, chr10–chr19, chr2, chr20–chr22,
chr3–chr9, chrEBV, chrM, unplaced, chrX, chrY.

* **chr1–22** blocks are exactly the node sets of `vg chunk -C -p GRCh38#0#chrN`: the whole
  connected component, including every node that is not on GRCh38 (insertions and
  alternative branches, 17–22 % of each component). They come from `scripts/build_chr_node_filters.sh`.
  Each list must be one gap-free interval.
* **Non-autosomal groups** (chrEBV, chrM, chrX, chrY; `*_random`, `chrUn_*` and HPRC contigs
  assigned to no chromosome → `unplaced`) take the remaining IDs, split at each group's
  smallest GRCh38 node.
* **Walk check:** every walk of every sample in the GFA (one walk = one assembly contig)
  must stay inside one block. A walk leaving an autosome block is an error.

HPRC v1.1 d9 result: 27 blocks cover all 60,118,570 nodes. None of the 8,688 walks of 46 samples
crosses a block.

The same table drives the builder's `--chromosomes all|autosome|chr1,chr2,...`
(`select_nodes`). Target nodes outside the chosen blocks are dropped before partitioning
(`orchestrate prepare`) and before the first batch (`run build`).

## In a whole-genome run

`orchestrate prepare` freezes the options below into `config.json`. When every task has validated,
`orchestrate run` calls `orchestrate finalize` (merge, then labels). `finalize` can also be run on
its own. It skips finished steps: the merge when `outputs.json` has a `merge` section, and the
labels when every kind has a `labels.manifest.json`.

```
--merge-shard-size 32768      0 = keep the task layout (no merge, no labels)
--keep-sources                keep task_* after the verified merge (default: delete them)
--chr-index TSV               needed for the merge and for --chromosomes other than all
--reference-path DIR          GRCh38 coordinates in the merged summaries; needed for labels
--somatic-vcf --somatic-bed --germline-vcf --germline-bed --reference-fasta --truth-dir   (labels; all or none)
```

`finalize` labels without an SNV AF filter. To apply one to a short-read set, label it again with
`label --snv-min-af 0.07` (`tools/jobs/relabel.sh ... 0.07`; see [Commands](#commands)).

## Merge (`merge_shards.py`)

Input: `<root>/outputs.json` of a completed `orchestrate run`. Output per kind (SNV, INDEL):

```
<tensors>/<kind>/<chrom>_shard_NNNNN_data.npy     chr1..chr22, --shard-size tensors each, last one shorter
<tensors>/<kind>/<chrom>_variant_summary.ndjson   source record, plus chrom, shard_file, new shard_index /
                                                  index_within_shard, source_task / source_shard_index /
                                                  source_index_within_shard, grch38 {chrom, pos0, ref, alt, node_reverse}
<tensors>/<kind>/filtered_candidates.ndjson       all task audit streams, task order
<tensors>/<kind>/unsupported_events.ndjson
<tensors>/<kind>/manifest.json                    layout "chromosome-shards", the shared keys (below), per chromosome
                                                  its shards (file, tensors, SHA-256) and summary, chr_index, sources
                                                  (run root, outputs.json SHA-256, tasks)
<tensors>/<kind>/validation_report.json
<tensors>/non_autosomal/<kind>/...                chrX, chrY, chrM, chrEBV, unplaced (not training data)
<root>/batch_timing.ndjson                        every task's batch timings
```

**One rule set per kind.** Every task manifest of a kind must agree with the first task's on
`tensor_format`, `tensor_storage`, `shape`, `dtype`, `channels`, `encodings`, `row_selection`,
`row_order`, `window_encoding`, `parameters`, `sample_unit`, `site_definition`, `read_cap` and
`debug_rows` (`SHARED_KEYS`). These keys, plus `graph_index` and `gai_version`, are copied into the
merged manifest. The read cap (`max_node_reads` and the rule text) is compared like every other
shared key, so a run whose tasks disagree on it is refused instead of merged under mixed rules.

The copy runs in parallel (`--workers`, default 8): every (kind, chromosome) group is one job that
reads its records' rows with plain sequential file reads and appends them to its shards, and the
audit streams are copied in slices straight to their offsets in the concatenated file. The bytes
are those of a sequential copy.

Records keep source order: task index, then shard, then index within the shard. The tasks of one
run cover contiguous node ranges in order, so its merged shards are node-sorted. A merge over tasks
of two node lists (e.g. a run followed by extra tasks for nodes added later) holds the first list's
tensors followed by the second's. Use `node_id` (or `grch38`) when position order matters.
`grch38` is null when the candidate has no GRCh38 position: the node has no unique GRCh38 visit
(alternative branches), or a multi-node deletion's nodes are not consecutive on GRCh38. `pos0` is
the 0-based REF start (INS: the boundary between `pos0 - 1` and `pos0`). Bases are on the contig's
forward strand, without the anchor base.

Everything is written to a hidden `.merging/` directory. It is published only when all of these
hold:

* every output shard's data SHA-256 (re-read from disk) equals the hash of the source tensors
  copied into it;
* every summary record equals its source record apart from the position fields;
* totals match `outputs.json`;
* shard headers, lengths and file sizes agree;
* random tensors reloaded from the source task shards equal the merged ones.

After publishing, the task directories are deleted unless `--keep-sources` is given. `outputs.json`
gets a `merge` section (the original is kept as `outputs.pre_merge.json`), and `status.json` becomes
`finalized` once the sources are gone. A second merge of the same run is refused, and so is a merge
into a directory whose manifest already has a merged layout.

**Earlier layout name.** Readers accept `chromosome-shards-v1` as `chromosome-shards`
(`common.EARLIER_FORMAT_NAMES`). The six existing merged sets (HG008 PacBio / ONT-UL / Illumina,
COLO829T Illumina / fiberseq / ONT) use that name. `label` reads only `layout` and `chromosomes`
from a merged manifest, so it relabels those sets as they are.

## Labels (`truth_labels.py`)

1. **Truth alleles**: every ALT on chr1–22 is trimmed of shared prefix/suffix bases. Symbolic
   alleles, `*`, and alleles containing N are skipped. MNVs become one SNV per differing base.
   Complex replacements are kept only as "nearby truth". An SNV or DEL whose REF differs from the
   FASTA is dropped (`reference_mismatch` in the truth stats). FILTER `.` counts as PASS; GT is the
   first sample's.
2. **Equivalent positions**: an indel in a repeat is enumerated at every equivalent placement
   on GRCh38, shifted left and right (at most 5,000 placements). VCFs are left-aligned. The builder
   left-normalizes on the node's forward strand, which on a reverse-oriented node is the other end
   of the repeat.
3. **Node keys**: each placement becomes the builder's candidate identity
   `node:start:KIND:REF>ALT` in node-forward coordinates. Reverse-oriented nodes flip the
   position and reverse-complement the bases. An insertion at a node junction gets a key on
   both nodes. Only nodes with one GRCh38 visit get keys. A deletion over several consecutive
   reference nodes of one orientation gets the builder's multi-node key: `start` on the
   forward-first node, REF over all nodes, and the further nodes as an `@n2+n3` suffix
   (`node:start:DEL:REF>@n2+n3`). A span over nodes of mixed orientation has no key, because the
   builder only joins a deletion over mappings of one orientation.
4. **Labels** are decided from the tensor's representative allele (A1, `candidate_id`). The site's
   other alleles (A2, A3, ...) are in the record's `alleles`. A tensor allele "is" a truth allele
   when its key equals one of that truth allele's keys. The confident region is somatic BED ∩
   germline BED.

| value | name | rule |
|---:|---|---|
| 1 | somatic | A1 is a somatic truth allele with FILTER PASS/`.`, inside or outside the BEDs; or A1 overlaps such an allele by more than 45 % (partial, below: another allele of the site anywhere, a truth allele at the same place inside the confident region) |
| 2 | germline | A1 is a germline truth allele with FILTER PASS/`.` or only `GAP1`/`GAP2`, inside or outside the BEDs; or A1 overlaps such an allele by more than 45 % (partial) |
| 0 | non | another allele of the site is such a truth allele and A1 overlaps it by 45 % or less (inside or outside the BEDs); every other tensor with a GRCh38 position inside the confident region: no truth allele, a truth allele nearby that A1 does not overlap enough (errors and artifacts next to real variants), or A1 is a germline allele with another FILTER (dipcall `HET1`/`HET2`) |
| −1 | ignore | an SNV below `--snv-min-af` (whatever the truth); a filtered somatic truth allele; no position; outside the confident region; an off-reference tensor without a partial somatic match |

`GAP_FILTERS` (`GAP1`/`GAP2`, dipcall): the allele is on one assembled haplotype and the other is
uncalled. It is in the normal genome (zygosity unknown), so it counts as germline truth. `HET1`/`HET2`
mean the assembled haplotype itself is ambiguous, so they are not a variant call.

`classify()` applies these checks in order; the first that applies decides. This table lists every
reason the code returns:

| # | check | label | reason |
|---:|---|---:|---|
| 1 | `--snv-min-af` given, SNV tensor with `af` below it | −1 | `below_snv_min_af` |
| 2 | A1 is a somatic truth allele: one with PASS | 1 | `representative_allele_is_somatic_truth` |
| | only filtered ones | −1 | `somatic_truth_filtered` |
| 3 | A1 is a germline truth allele: one with PASS | 2 | `representative_allele_is_germline_truth` |
| | else one with only GAP1/GAP2 | 2 | `representative_allele_is_germline_truth_gap_filtered` |
| | else (e.g. HET1/HET2): noted for step 8 | | |
| 4 | another allele of the site is a PASS somatic or PASS/GAP germline truth allele: A1 overlaps a somatic one by > 45 % | 1 | `allele_partial_somatic_truth` |
| | else a germline one by > 45 % | 2 | `allele_partial_germline_truth` |
| | else | 0 | `truth_matches_non_representative_allele` |
| 5 | another allele of the site is a filtered somatic truth allele | −1 | `somatic_truth_filtered` |
| 6 | no GRCh38 position and no anchor (off-reference tensors, below) | −1 | `not_on_unique_grch38_node` |
| 7 | position outside the confident region | −1 | `outside_confident_region` |
| 8 | A1 is a germline truth allele with another FILTER (step 3) | 0 | `germline_truth_filtered` |
| 9 | a PASS somatic truth allele at the same place that A1 overlaps by > 45 % | 1 | `residual_partial_somatic_truth` |
| 10 | off-reference tensor (anchored), nothing above | −1 | `off_reference_no_truth_match` |
| 11 | a PASS/GAP germline truth allele at the same place that A1 overlaps by > 45 % | 2 | `residual_partial_germline_truth` |
| 12 | a truth allele of either set, any FILTER, within 10 bp (`NEAR_BP`) | 0 | `near_truth_allele_mismatch` |
| 13 | otherwise | 0 | `confident_no_truth_allele` |

Steps 1–5 do not look at the BEDs. A site whose other allele is a PASS somatic or PASS/GAP germline
truth allele is labelled by that allele (1, 2 or 0) wherever it lies. Step 8 never applies to an
off-reference tensor, because only GRCh38 nodes have truth keys. With several truth alleles to
compare, the best overlap decides, and somatic is tried before germline.

**Positions.** A GRCh38 tensor's interval is its REF span `[pos0, pos0 + len(REF))` (one base for an
SNV); an INS uses the two bases around its insertion point, `[pos0 − 1, pos0 + 1)`. The
confident-region test needs the whole interval inside one interval of the region. For an anchored
off-reference tensor, the test uses the middle base of its anchor interval. "The same place" (steps
9, 11) means that the truth allele's span of equivalent placements intersects or touches the
tensor's interval; for an anchored tensor it means within 10 bp of the anchor interval. "Within 10
bp" (step 12) means fewer than 10 bases between the tensor's interval and the truth allele's whole
span of equivalent placements. An insertion in a repeat spans from its leftmost to its rightmost
equivalent boundary, because reads may write it anywhere in the repeat. Both truth sets, every
FILTER, and complex alleles count as nearby truth.

**Partial matches** (`MIN_OVERLAP` = 0.45: "more than 45 %"). The same event written differently by
the graph alignment takes the truth's label. `<chrom>_labels.ndjson` marks it with `partial`,
`overlap` (3 decimals) and `partial_truth`, and `labels.manifest.json` counts it under `partial`
(`allele_somatic`, `residual_germline`, ...).

* `partial: "allele"` (step 4): another allele of the site is the truth allele. For example, A1 +6C
  next to the truth +5C (5 of 6 bases, 83 %), A1 DEL AA next to the truth DEL AAA (67 %), and DEL AA
  next to DEL AAAA (50 %) take the truth's label. DEL A next to DEL AAA (33 %) and +4C next to +9C
  (44 %) are 0 (`truth_matches_non_representative_allele`) unless, for a somatic truth, the
  haplotype overlap (below) is higher.
* `partial: "residual"` (steps 9, 11): a truth allele at the same place that no allele of the site
  equals. chr1:201 DEL AA overlaps the truth chr1:200 DEL AAA (2 of 3 bases, 67 %) and takes its
  label. chr1:210 DEL AA is not at the same place, so it is 0 (`near_truth_allele_mismatch`).
* **Allele overlap** (both truth sets):
  * DEL/DEL: the deleted bases both remove, positions included. The tensor's deleted interval is
    compared against the truth's span of equivalent placements, capped at the shorter deletion.
  * INS/INS: the inserted bases both carry at the same insertion point. The tensor's point must lie
    in the truth's span; the score is the longest common subsequence with the truth insertion
    rotated to that point.
  * Both are divided by the length of the longer allele. Every other pair is 0 (another SNV base;
    different kinds).
* **Haplotype overlap** (somatic truths only, when the allele overlap is 45 % or less): it uses up
  to 20 A1 rows and 10 REF rows (`EVIDENCE_ROWS`, evenly spaced; read bases of channel 0, reads of
  at least 30 bases). Each read is fitted inside the GRCh38 window of the truth allele (90 bases on
  each side) with and without the truth allele. Let d_ref and d_truth be those edit counts minus the
  REF reads' median error, and let e be the event size (1 for an SNV, otherwise the longer allele).
  A read's overlap is max(0, (d_ref + e − d_truth) / 2) / max(d_ref, e), and the tensor's is the
  median over its A1 reads. With fewer than 3 A1 reads there is no haplotype overlap. This catches
  a truth written as a graph branch plus a residual edit (truth +AAA, tensor DEL A after a +AAAA
  branch), or as a skipped node plus a mismatch. It is also the only way an off-reference tensor
  matches, with each read taken in the orientation that fits better.

**Off-reference tensors** (`grch38` null) are placed by node ID (`anchor()`; the interval is
recorded as `anchor` in labels.ndjson). Minigraph-Cactus numbers nodes in topological order, so a
branch node's ID lies between the IDs of its flanks. The anchor interval runs from the nearest
unique reference node below to the nearest one above, each searched within 200 IDs (`ANCHOR_REACH`),
both nodes included; when only one side has such a node, the interval is that node. On HG008 PacBio,
86 % of the branch-node residual edits of missed somatic INDELs lie within 10 bp of that interval.

A tensor has no anchor (−1, `not_on_unique_grch38_node`) when there is no reference node within
200 IDs on either side, when the two neighbours are on different contigs, or when more than 1024
GRCh38 bases lie between them (`ANCHOR_GAP`). The nodes' own lengths (up to 1024 bp each) are not
counted, so an SNV bubble between two long nodes keeps its position. The gap limit matters around
centromeres, where the neighbours by ID lie up to 213 Mb apart and the interval would match any
truth. Without the limit, 120 of the 194 COLO829T fiberseq SNV tensors with an interval over 1 Mb
became partial 1, and the label job spent 45 and 54 min on chr10 and chr2 against minutes for most
chromosomes. Such gaps are rare: 676 of the 106,788 anchored COLO829T fiberseq SNV tensors have one.

An anchored tensor inside the confident region is 1 on a partial somatic match (step 9) and −1
otherwise (step 10). Germline residuals are not tried on off-reference tensors.

**Why −1 is kept narrow.** A tumor-only caller meets every tensor at test time. So, apart from
off-reference tensors (below), −1 is used only for filtered somatic truth and for tensors that the
caller can also leave out without truth: outside the BED it calls in, no position, or an SNV under
the AF filter. Everything else it will meet is 0 or a truth label, errors and artifacts next to
real variants included. On HG008 PacBio, 35,259 of the 86,625 SNV 0s and 205,094 of the 1,518,564
INDEL 0s are `near_truth_allele_mismatch`.

* **`--snv-min-af`** (the short-read sets use 0.07) labels as if the build had used that SNV AF
  filter: a lower-AF SNV tensor is −1 whatever the truth, and this is the first check. `af` is in
  every summary record, so a caller can apply the same filter without truth.
* **Off-reference tensors without a partial somatic match** are −1 for a different reason: they are
  not training negatives, because most branch nodes carry no truth at all (the truth sets are
  GRCh38 VCFs). A partial somatic match still makes one 1.

**Outputs**, next to the merged shards and in summary order:

* `<chrom>_shard_NNNNN_labels.npy`: int8, one per data shard.
* `<chrom>_labels.ndjson`, one line per tensor with:
  * `candidate_id`, `site_id`, `chrom`, `shard_file`, `index_within_shard`, `label`, `label_name`,
    `reason`;
  * `somatic` / `germline`: the truth alleles that an allele of the site equals, each with truth_id,
    VCF fields, GT, filters, in_bed, the matched candidate and whether it is A1;
  * `grch38`;
  * `partial`, plus `overlap` and `partial_truth` when set;
  * `anchor`, for anchored tensors.

The label files are written as temporaries and renamed only after every shard is covered.

`labels.manifest.json` (one per kind) holds:

* `format` `"truth-labels"`;
* `rules_sha256`: the SHA-256 of `truth_labels.py`. The rules are that module and its constants,
  so any edit of the file (comments included) changes the hash;
* `created`, `labels` (name → value), `near_bp`, `tensors`;
* `counts` (per chromosome and label name), `totals`, `reasons`, `partial`;
* `min_overlap`, `anchor_reach`, `anchor_gap`, `snv_min_af` (null when not given);
* the provenance: `truth` (per set: VCF and BED with their SHA-256, parsing stats, truth table
  path), `confident_region` (definition and bp), `fasta`, `reference_path`.

The truth tables `<set>.graph.tsv` go to `--truth-dir`. They list every truth allele with its
trimmed form, kind, GT, filters, in_bed, number of placements and keys, and depend only on the graph
and the truth inputs (VCF, BED, FASTA), not on the run. `<set>.recall.tsv` and `truth_recall.json`
go to `--recall-dir` (default: the tensor directory). Per truth allele, they give one status:

* `tensor_representative`
* `tensor_non_representative_allele`
* `complex_allele`
* `no_unique_grch38_node_key`
* `filtered` (with the builder's reasons)
* `no_candidate`

Recall counts key matches only; partial matches are not in it.

## Commands

From the repository root, with `P="python -m indexed_gam_pipeline_v4"`:

```bash
# once per graph
$P.tools.graph_prep ref-path-scan --gfa G.gfa --output DIR [--reference-sample GRCh38]
$P.tools.graph_prep ref-path-check --path DIR --graph-index DB --fasta FA [--samples 100000]
$P.tools.graph_prep chr-index --components-dir D --reference-path DIR --output PREFIX [--graph-index DB]
# manual re-merge / re-label of a run root (orchestrate finalize does both)
$P.tensor_postprocessing merge --root /path/to/run --chr-index TSV [--shard-size 32768] [--keep-sources] \
    [--workers 8] [--spots 200] [--reference-path DIR]
$P.tensor_postprocessing label --tensors /path/to/tensors --reference-path DIR --fasta FA \
    --somatic-vcf VCF --somatic-bed BED --germline-vcf VCF --germline-bed BED --truth-dir DIR \
    [--kinds SNV INDEL] [--recall-dir DIR] [--snv-min-af 0.07]
```

HPRC v1.1 d9 graph files. The files under `$G` are the current ones; the reference-path and
chr-index JSON record their formats under the earlier key and name (see the top of this file).

```bash
G=/scratch/jshen/data/pansoma_v2_tensors/graph_index
$P.tools.graph_prep ref-path-scan \
    --gfa /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa --output $G/hprc-v1.1-mc-grch38.d9.grch38_path
$P.tools.graph_prep ref-path-check --path $G/hprc-v1.1-mc-grch38.d9.grch38_path \
    --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite \
    --fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
$P.tools.graph_prep chr-index \
    --components-dir /scratch/jshen/data/AF-Filtered_VG_Indexes/chr_component_vs_GRCh38_summary \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path --output $G/hprc-v1.1-mc-grch38.d9.chr_node_ranges \
    --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite
```

`tools/jobs/graph_prep.sh OUTDIR` runs the same three commands as one Slurm job, with these inputs
as its defaults (other graphs: `graph_prep.sh OUTDIR GFA GRAPH_INDEX FASTA COMPONENTS_DIR`); it
names the outputs after the GFA, as above, and does not overwrite them.

Measured on HPRC v1.1 d9:

* ref-path-scan: 7.4 min / 4.0 GiB (49,092,514 GRCh38 nodes, none visited twice);
* ref-path-check: 0 length and 0 sequence mismatches;
* chr-index: 29 s;
* `graph_prep.sh` (all three, 2 CPUs, `--mem=5G`): 10 min, peak 4.0 GiB. Its outputs equal the
  files under `$G`: every `.npy`, `walks.ndjson` and the `.tsv` byte for byte; the JSON only in the
  format key (`format`, plain name, instead of `version`, `-v1` name), `scan_seconds` and paths.

**Relabelling the six existing sets.** `tools/jobs/relabel.sh` labels one merged set with the
package it is in (main README, section 4, "Tools"). On the data side it runs from
`/scratch/jshen/data/pansoma_v2_tensors/pipeline_code/`, a `git archive` of this package plus the
compiled `.so`, with the commit recorded in `pipeline_code/git_head.txt`. Before labelling, the
script does two things:

* checks that `TENSORS/SNV/manifest.json` has a merged layout (either name);
* copies the current `labels.manifest.json` files and recall reports to
  `<sample dir>/labels_backup_<tensors dir>_<time>_<job>/`.

It then runs `tensor_postprocessing label` with the HPRC reference-path directory and the GRCh38
FASTA above (1 CPU, `--mem=19G`, 6 h). The COLO829T `label_job.sh` scripts (step 4 of those runs)
call it with the COLO829T truth files below.

```bash
D=/scratch/jshen/data/pansoma_v2_tensors; H=/scratch/jshen/data/HG008_GIAB; Q=/scratch/qfu/COLO829BL_DSA/dipcall_hg38
J=$D/pipeline_code/indexed_gam_pipeline_v4/tools/jobs
# sbatch -J NAME -o LOG $J/relabel.sh TENSORS SOMATIC_VCF SOMATIC_BED GERMLINE_VCF GERMLINE_BED TRUTH_DIR [SNV_MIN_AF]
HG008="$H/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz
       $H/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed
       $H/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz $H/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed"
COLO="$D/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz $D/COLO829T_truth/SMaHT_v2_easy_difficult_extreme.union.bed
      $Q/dipcall_hg38.dip.vcf.gz $Q/dipcall_hg38.dip.bed"
relabel() { sbatch -J relabel_$1 -o $D/$1/slurm-relabel-%j.out $J/relabel.sh $D/$1/v3_tensors $2 $D/$1/truth $3; }
relabel Liss_lab_PacBio_Revio_20240125 "$HG008"             # HG008 PacBio HiFi
relabel Liss_lab_Northeastern-ONT-UL-20241216 "$HG008"      # HG008 ONT-UL
relabel Liss_lab_BCM_Illumina-WGS_20240313 "$HG008" 0.07    # HG008 Illumina
relabel COLO829T_Illumina "$COLO" 0.07
relabel COLO829T_fiberseq "$COLO"
relabel COLO829T_ONT "$COLO"
```

Truth sets:

* **HG008:** the GIAB HG008-T somatic small-variant draft benchmark (tumor-variants VCF and
  `_all.bed`), and HG008-N dipcall (`dip.vcf.gz`, `dip.bed`).
* **COLO829T:** the validated SNV + INDEL union VCF; as somatic BED, the union of the SMaHT
  easy/difficult/extreme `_v2` BEDs (they tile the genome, so their intersection is empty); and
  COLO829BL dipcall.

Label job, 1 CPU:

| set | wall time | MaxRSS |
|---|---|---|
| HG008 PacBio | 25 min | 13.9 GiB |
| HG008 ONT-UL | 46 min | 14.0 GiB |
| HG008 Illumina | 53 min | 13.9 GiB |
| COLO829T Illumina | 32 min | 15.7 GiB |
| COLO829T fiberseq | 32 min | 15.7 GiB |
| COLO829T ONT (7.7 M tensors) | 2 h 12 min | 15.7 GiB |

The label counts of the six sets are in the [main README](../README.md), section 5.

A change to the label rules (any edit of `truth_labels.py`, through `rules_sha256`) or to the merged
bytes changes the goldens (`tests/golden_hashes.json`).
Record them again from the committed change (main README, "Testing and goldens").
