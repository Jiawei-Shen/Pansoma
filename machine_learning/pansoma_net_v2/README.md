# pansoma_net_v2

PansomaNetV2 trains and predicts directly on the merged Pansoma tensor sets. Each set is a directory with
`SNV/` and `INDEL/` from `tensor_postprocessing` merge + label. The tensors are (8, 200, 101) int8. Labels
are 0 non, 1 somatic, 2 germline and −1 ignore.

It differs from v1 (`../pansoma_net`: one `.npy` file per sample in class folders, 6 channels, z-score with
fixed means) in three input changes, all done on the GPU. The tensors on disk do not change.

1. **One-hot categories.** Read base, site allele, graph base and alignment operation (6 codes each) and
   strand (2) become 0/1 planes, so no base or operation is "larger" than another.
2. **Derived planes and a 1×1 front.**
   - `differs` is read base ≠ graph base.
   - The row's site block is taken from the summary's `row_groups`, as in DeepVariant's
     `read_supports_variant`: A1 (the representative allele), ALT (A2…Ak), REF or OTHER, set on every
     covered cell of that read.
   - Two 1×1 convolutions (64, 64, GELU) combine the planes of each cell before the 4×4 stride-4 stem
     mixes cells.
3. **Normalization that keeps 0.**
   - One-hot planes are not normalized.
   - BQ, MAPQ and path count get a masked z-score. Valid cells are BQ ≥ 0, MAPQ ≥ 0 and path count > 0; the
     z-score is 0 elsewhere.
   - **BQ and MAPQ are z-scored per dataset** (`train --norm dataset`, the default from 2026-10-13): each
     dataset (a set's `SNV/` or `INDEL/`) with the mean and std of its own valid cells, from up to `--stats-samples`
     (20,000) of its tensors, labels not used. Training fits one row per training dataset (a joint run on P100 + P75
     + P50 + P25: four rows; checkpoint `dataset_stats`, `dataset_stats.json`, kept on `--resume`); `predict` fits
     one on each predicted dataset (`.metrics.json` `normalization`). BQ and MAPQ scales depend on the
     instrument, run and aligner (4-level binned BQ on HG008T Illumina, a median of 30 on COLO829T Illumina; MAPQ
     std 2.3 on HG008T PacBio, 8.9 on COLO829T Illumina), so with the training statistics a model saw another
     dataset's values shifted. This aligns datasets' means and spreads, not the shape of their distributions, and
     not the reads of two runs mixed inside one dataset.
   - Path count keeps the training statistics (checkpoint `stats`, fitted once on the training tensors): it is a
     property of the graph, stored relative to its haplotypes, not of the sequencing. `--norm train` (and every
     checkpoint before 2026-10-13, config without `norm`) z-scores BQ and MAPQ with the training statistics too.
   - A `covered` plane (ch4 > 0) and a `bq_missing` plane (a read cell with BQ −1) mark real data, because
     MAPQ 0 and BQ 0 are legal values.
   - Every plane is exactly 0 on padding.

The width is zero-padded from 101 to 104 so the stem also covers column 100. The v1 stem skips it.

## Planes (36, `encode.PLANES`)

| planes | source | values |
|---|---|---|
| `read_base=1..6`, `site_allele=1..6`, `graph_base=1..6`, `operation=1..6` | ch0, ch2, ch5, ch4 | 0/1 (A C G T N gap; M X I D complex aligned-gap) |
| `strand=1..2` | ch7 | 0/1 |
| `base_quality`, `mapping_quality`, `path_count` | ch1, ch3, ch6 | masked z-score, 0 where invalid or no read |
| `covered`, `bq_missing` | ch4 > 0; covered and ch1 = −1 | 0/1 |
| `differs` | covered and ch0 ≠ ch5 | 0/1 |
| `row_A1`, `row_ALT`, `row_REF`, `row_OTHER` | `row_groups` of the summary | 0/1 on the covered cells of the row |

## Use

One model per sample, platform and kind (SNV or INDEL). The test chromosome is chr1; it is never used for
training or selection.

```bash
cd machine_learning
P=/wanglab/jshen/anaconda3/bin/python       # torch 2.8 + timm 1.0 (GPU)
T=/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors

# train on chr2-22 (5 % of their ~1 Mb node blocks validate), chr1 left out; -1 tensors are not used
# (the recipes: see "SNV recipe" and "INDEL recipe" below)
$P -m pansoma_net_v2.train --tensors $T --kinds SNV --output runs/HG008_Illumina_SNV \
    --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75 \
    --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
$P -m pansoma_net_v2.train --tensors $T --kinds INDEL --output runs/HG008_Illumina_INDEL \
    --depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --batch-size 1024 --lr 1e-4 --class-weights sqrt \
    --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
torchrun --nproc_per_node=2 -m pansoma_net_v2.train --ddp ...           # several GPUs, same options

# test: every chr1 tensor (-1 included), the checkpoint's statistics and somatic threshold
$P -m pansoma_net_v2.predict --checkpoint runs/HG008_Illumina_SNV/best.pth --tensors $T --kinds SNV \
    --chroms chr1 --output runs/HG008_Illumina_SNV/test_chr1
```

`machine_learning/pansoma_net_v2/jobs/run.sh NAME TENSORS KIND [train options]` runs both steps as one Slurm job.
The Slurm job scripts in `jobs/` (training + chr1 test, chr1 test again, calls, re-scoring of calls, the INDEL
recipe's submitter, the index and chunk caches; `jobs/README.txt` lists them) are local copies for this cluster and
not in git (`.gitignore`); their logs and `submitted.txt` go to `/scratch/jshen/data/pansoma_net_v2_runs/jobs`. On RIS the jobs are `ris/` (`ris/README.md`).
`python -m pansoma_net_v2.summary [NAME ...]` prints each run's per-epoch table, its current step and, when
done, its chr1 test result (per tensor and against the truth VCF) and the rtg vcfeval reports of its calls
(`<run>/vcf_*/eval_bed`, `eval_nobed`: PASS P/R/F1 and the chr1 curve's best point, raw and after the PoN). The table
columns are:

- training loss, speed, minutes, GPU peak;
- validation per tensor: somatic AP and s.F1 / s.P / s.R at the best threshold, the threshold t, and the precision
  where the somatic recall reaches 0.9 (P@R.9);
- against the truth VCF at the same t: t.AP, t.F1 / t.P / t.R and the ceiling (below);
- germline argmax F1; `*` marks the best epoch so far (by s.F1, as `best.pth`).

Useful `train` options (the recipes, "SNV recipe" and "INDEL recipe" below, set the size, `--epochs`,
`--batch-size`, `--lr`, `--class-weights`, `--non-fraction 1.0`, `--block v2` and `--ignore-reasons`):

- **Split.**
  - `--test-chroms` (default chr1; `--test-chroms` alone: none) leaves chromosomes out of every set;
    `--set-test-chroms SET CHROM ...` (repeatable) out of one set only, e.g. a test sample's chr1 while the other
    sets train on theirs (the checkpoint's `chroms.set_test`).
  - `--val-fraction 0.05` of the `--val-block-nodes 20000`-node blocks of the training chromosomes; whole
    blocks, so validation tensors share no reads with training tensors. `--val-chroms` instead validates on
    whole chromosomes; `--train-chroms` limits training to some chromosomes (a quick test).
- **Imbalance.**
  - `--class-weights balanced|sqrt|pow:E|none|w0,w1,w2`: balanced is n / (3 n_c) over one epoch's tensors;
    sqrt is its square root, pow:E its power E. The loss is normalized by the batch's summed weights, so
    scaling all three weights changes nothing.
  - `--non-fraction 1.0` (the default and the recipes): every tensor every epoch (smaller values: below).
- **Label reasons left out.** `--ignore-reasons R ...` treats the tensors of those label reasons as -1 in training and
  validation; the label files are unchanged. Both recipes leave out
  `residual_partial_somatic_truth allele_partial_somatic_truth`. For SNV: SNV tensors that partially match an INDEL
  truth (mostly read-end SNV shadows of the INDEL) taught the model low-AF artefacts next to germline INDELs, and as
  SNV calls they never match the truth (analysis/v1_vs_v2_20260929).
- **Model.**
  - `--depths`, `--dims`, `--front`, `--drop-path` set the size (defaults `3 3 27 3` / `192 384 768 1536`,
    as in v1's training script).
  - `--block v2` (the default; below, Model).
  - `--drop-planes PLANE ...` leaves encoder planes out of the input (e.g. `path_count`; "Channel 6 storage" below).
- **Selection.** `best.pth` is the epoch with the highest thresholded per-tensor somatic F1 on the validation tensors
  (ties: the lower validation loss); it stores that threshold.
- **Speed.** `--no-compile` and `--no-channels-last` turn off the defaults on CUDA (below). Training reads
  chunks from a zstd cache (`--read-threads 16`, `--chunk-rows 512`, `--window-chunks 32`, `--chunk-cache`;
  "Reading the training tensors" below). `--num-workers` sets the DataLoader workers of the statistics (0 also reads
  the validation without workers), `--val-workers` those of the validation.
- **Memory.** `--micro-batches K` runs the forward and backward of each batch in K parts (gradient accumulation):
  the same update as the whole batch, with the loss normalized over the whole batch, and the GPU memory of one part.
  Batch 1024 of the default model fills a 94 GB H100 NVL (about 94.7 of 95.8 GB); batch 512 takes about 49 GB.
- **Normalization.** `--norm dataset|train` (above, "Normalization that keeps 0"); `--stats-samples 20000`.
- Also: `--amp bf16|off`, `--resume` (keeps the checkpoint's statistics, the per-dataset ones too), `--seed`.

Options kept for older checkpoints and experiments (no recipe uses them):

- `--block v1` and `--decay-all`: the blocks and the weight decay of the checkpoints made before 2026-09-29 (below,
  Model).
- `--scalars` feeds the site scalars to the head (`data.SCALARS`: log coverage, log site coverage, log
  ALT / REF / OTHER counts, AF, second allele AF, allele count, log event length). They are the counts
  before the 200-row cap, z-scored with training statistics that are stored in the model.
- `--non-fraction F` < 1: each epoch takes every somatic and germline tensor and a fresh random fraction F of the non
  tensors. `--keep-non-af A` with it takes the non tensors with AF ≥ A every epoch; the sampled non tensors then get
  loss weight 1/F, so each epoch's loss stays an unbiased estimate of the loss over all non tensors, whatever their AF.
- `--exclude-candidates FILE` (candidate ids, first column) takes those tensors out before the model, as if never
  built: -1 in training and validation; `predict` leaves them out of the test by default (their truth alleles are
  misses) and `graph_vcf` writes no record for them. The INDEL recipe does not use it (below).
- `--epoch-samples N` (a fresh random subset of the training tensors each epoch) and `--val-samples N` (a fixed random
  subset of the validation tensors).
- `--no-chunk-cache`: training reads the chunks from the raw shards instead of the zstd cache.

**Evaluation** (`metrics.py`). The target is the somatic class; germline is reported too.

- **Which tensors count.** Validation and the test chromosome score the tensors labelled 0/1/2, plus the
  off-reference tensors without a truth match (`off_reference_no_truth_match`, −1 for training) as non:
  test-time calling meets them inside the BED, and a somatic call on one is a false positive. The other −1
  are left out, because calling drops them without truth: outside the BED (`outside_confident_region`),
  below the AF floors (`below_snv_min_af`, `below_indel_min_af`), and no GRCh38 position
  (`not_on_unique_grch38_node`). Training uses only 0/1/2.
- **The confident region.** Labels 1 and 2 outside the confident region (somatic BED ∩ germline BED of the
  labels, the labeller's test in `bed.py`) are trained on but not scored: the test drops them with the BED,
  as it drops the outside −1.

- Each validation reports the somatic and germline average precision (PR-AUC), the argmax precision /
  recall / F1 per class, and the somatic threshold t of the best F1.
- With t, a tensor is somatic when p_somatic ≥ t, and otherwise the larger of non and germline. The
  checkpoint stores t (`somatic_threshold`), and `predict` applies it to the test chromosome.

**Against the truth VCF** (`metrics.truth_report`, the `truth` entry of each report). Per tensor, a truth
allele can count twice (two partial INDEL tensors of one truth) or not at all (no tensor). So each report also
scores the truth alleles themselves:

- The truth alleles are the somatic truth rows of the set's `somatic.recall.tsv` that pass and lie in the
  BED, of the model's kind (SNV: SNP; INDEL: DEL, INS), on the scored chromosomes.
- Each labelled-1 tensor stands for the truth alleles its labels record: the representative allele's
  matches, or `partial_truth` for a partial match.
- A truth allele scores the highest p_somatic of its tensors. A truth allele without a tensor is a miss at
  every threshold.
- A false call is a somatic call on a tensor scored 0 or 2 (off-reference no-match tensors included).
- A tensor of the other kind's truth (an SNV tensor that partially matches an INS or DEL truth) is somatic for
  training and per tensor. Against the truth VCF, the truth allele it finds is a true call, counted once
  (`other_tp`), but not in the recall, which is over the model's kind. A labelled-1 tensor whose truth alleles
  count for neither kind (not PASS, outside the BED) is neither (`unmatched_tensors`).
- Precision = distinct truth alleles found (either kind) / (those + false calls); recall = found / all truth
  alleles of the kind; `ceiling` = the recall if every tensor were called; AP over the same ranking;
  `at_threshold` uses t, `best` the best-F1 threshold.
- Validation: a truth allele belongs to the validation block that holds its position (each block starts at
  its first tensor), so validation truth alleles without a tensor count as misses there too.

On HG008T chr1 (the recipe runs of `analysis/hg008t_snv_indel_7fc7396_20261012`; SNV Illumina: its earlier run
`HG008_Illumina_SNV_h90`, the same recipe and set):

| | truth alleles | with a tensor (ceiling): Illumina / PacBio / ONT |
|---|---|---|
| SNV | 697 SNP | 659 (0.945) / 658 (0.944) / 663 (0.951) |
| INDEL | 588 DEL + INS | 317 (0.539) / 413 (0.702) / 408 (0.694) |

Outputs in `--output`:

- `args.json` (the options) and `train.log`.
- `metrics.jsonl` (per epoch: train loss and accuracy, speed, GPU peak memory, and the validation report).
- `stats.json` (the training statistics) and, with `--norm dataset`, `dataset_stats.json` (one row per training
  dataset).
- `last.pth`, and `best.pth` (the epoch of the best thresholded validation somatic F1).

A checkpoint holds `format: pansoma_net_v2`, `config`, `model_state_dict` (the encoder's and scalar
statistics are buffers in it), `stats`, `dataset_stats` (norm dataset: the training datasets' BQ / MAPQ statistics;
config `norm`), `somatic_threshold`, `classes`, `planes`, `scalars`, the training
data (directories, label provenance, counts), the channel-6 storage of those sets (`path_count_storage`,
`path_count_haplotypes`), the chromosome split and `args`. Rebuild it with `PansomaNetV2.from_checkpoint(path)`
(the model only; `predict` and `train --resume` check the storage).

`predict` writes, for every tensor of the chosen chromosomes, `<sample>.<set>.<KIND>.predictions.ndjson.gz`.
Each record has chrom, candidate_id, label (truth), test_label (the scored label, null when left out),
in_test, reason, off_reference (the node is off the GRCh38 path, as in the paper's pangenome-only calls),
truth_ids (the somatic truth alleles of a labelled-1 tensor), p_non, p_somatic, p_germline, and pred (with
the threshold). It also writes `.metrics.json`: the same report over the test tensors, with the truth report
and the counts of test tensors on off-reference nodes (`--threshold` overrides t), and `normalization`: the BQ / MAPQ
z-score used (`norm`; for norm dataset the statistics fitted on this dataset's predicted tensors, `--stats-samples`
of them; `--norm` overrides the checkpoint's). `--ignore-reasons` (default: the
checkpoint's training `--ignore-reasons`) leaves those reasons out of the test like -1; their records keep label,
reason and truth_ids and say `ignored: true`. INDEL sets are the exception to the default: the partial somatic matches
(`residual_partial_somatic_truth`, `allele_partial_somatic_truth`) stay in the INDEL test as 1, so a called INDEL
tensor of another length in the same homopolymer (or a residual edit next to a graph allele) counts its truth as
found in the per-tensor and truth-level numbers. SNV tests leave them out. The calls, and so rtg vcfeval, do not
change: there such a call is an allele of the wrong length.
Off-reference tensors are marked for analysis only; the test scores them like the others.

**Off-reference rescue (SNV).** The model learned that a low path count (ch6: how many HPRC haplotypes run
through the node) is not somatic, and on an off-reference node every read has a low path count, so it hardly ever
calls somatic there. `predict --offref-site90 W` (default 10; 0 turns it off) scores every SNV off-reference tensor
again with its ch6 cells below the every-haplotype code raised to it in the site columns (the middle column ± W;
uncovered cells stay 0), and calls it somatic when somatic is the most probable class of those probabilities. The
every-haplotype code is the set's channel-6 storage's (`data.PATH_COUNT_ALL`): 100 on `int8-count-haplotypes100` (the
count relative to the build's `--haplotypes` H, H → 100), 90 on the older `int8-count-linear100-log2` (kept for the
checkpoints trained on it; "Channel 6 storage" below). Its record adds p_offref_non,
p_offref_somatic, p_offref_germline and offref_call. `.metrics.json` adds `offref_rescue` (tensors, calls before
and after, `floor`) and the rule `pipeline`: offref_call on the off-reference tensors, the checkpoint's validation recall-0.9
threshold on the others.

**SNV recipe** (HG008T Illumina, PacBio, ONT): the 50.4 M model `--depths 3 3 27 3 --dims 96 192 384 768` (the
default widths halved), `--class-weights pow:0.75 --batch-size 1024 --lr 1e-4 --epochs 20 --non-fraction 1.0
--block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth`, every plane kept (ch6
included), one whole GPU; calls at the validation recall-0.9 threshold, off-reference tensors by the rescue above
(graph_vcf's defaults). On HG008T PacBio chr1 (12 epochs) the default 199.6 M model at lr 1e-4 overfits after
epoch ~7 (validation loss 0.51 -> 2.16) and the 50.4 M model does not; with the PoN in the BED both reach P/R/F1
0.473/0.832/0.604 (199.6 M, best epoch 5) and 0.493/0.822/0.616 (50.4 M).

**INDEL recipe** (every platform; chosen on the six P100 sets): the 8.7 M model `--depths 2 2 6 2 --dims
64 128 256 512`, `--class-weights sqrt --batch-size 1024 --lr 1e-4 --epochs 12 --non-fraction 1.0 --block v2
--ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth`, every plane and every candidate kept
(no `--exclude-candidates`); calls at the validation recall-0.9 threshold (graph_vcf's default; no off-reference
rescue for INDELs), then the INDEL PoN rule (gnomAD / CoLoRSdb exact allele: AF ≥ 0.01 in homopolymers of ≥ 7 bases,
else AF ≥ 0.001; below). The chr1 curve's best point is recorded next to the PASS result (vcfeval's report,
`summary`); it is chosen on the test chromosome and never used as the threshold. On a MIG 2g.24gb slice the batch
needs `--micro-batches 4` (about 23 GB as one batch). `machine_learning/pansoma_net_v2/jobs/submit_indel.sh SET ...`
submits the recipe per set from that directory (chunk cache, training + chr1 test, calls + vcfeval); logs and
`submitted.txt` go to `pansoma_net_v2_runs/jobs`. The recipe's chr1 results on HG008T:
`analysis/hg008t_snv_indel_7fc7396_20261012`.

Why (chr1 in the BED after the INDEL PoN, F1 at the validation recall-0.9 threshold; the runs
`<sample>_<platform>_INDEL_small_*` in `pansoma_net_v2_runs`):

- **lr 1e-4.** With the candidate filter (below) HG008T Illumina / PacBio / ONT 0.241 / 0.135 / 0.017 against
  0.226 / 0.130 / 0.016 at 4e-4; without it the two are within 0.01 (`analysis/hg008t_snv_indel_7fc7396_20261012`).
- **sqrt class weights.** pow:0.8 lowers the validation somatic AP in six of the eight long-read pairs and the
  recall-0.9 F1 in seven (more calls at the same recall).
- **No candidate filter.** The filter (INDEL candidates with AF < 0.10, and with AF < 0.15 whose allele is not in the
  other sample's set of the platform) takes out tensors the model already ranks low; on HG008T the unfiltered runs
  are about equal or better (0.239 / 0.132 / 0.015 against 0.226 / 0.130 / 0.016). On COLO829T's long reads the filtered
  runs score higher, but on at most about 33 truth INDELs left after the PoN.
- **The recall-0.9 rule.** The validation best-F1 threshold t scores higher on the HG008T long reads (PacBio 0.214,
  ONT 0.101 against 0.132, 0.015) and lower on Illumina (0.194 against 0.239), and calls only 9–185 tensors on
  COLO829T; argmax called 0–2 INDELs in three Illumina runs.

The 8.7 M model is within noise of the 199.6 M model on chr1 (best F1 0.13–0.15 without the PoN, with or without the
partial labels or ch6) and trains about four times faster (HG008T Illumina: 55 s against 214 s per epoch). What
limits INDELs is outside the model: the PoN rule tags 413 of the 635 HG008T chr1 truth INDELs (81 of 131 on
COLO829T), so the curves' recall after it stops at 0.24–0.26; 46 % of the HG008T Illumina truth INDELs get no
candidate (the allele is a path of the HPRC graph); and
the false calls are the patient's germline homopolymer INDELs at the same VAF as the somatic ones
(`analysis/indel_20260930`, `analysis/tensor_recall_20260930`).

## Calls: graph VCF → GRCh38 VCF → PoN → rtg vcfeval

Three steps turn `predict`'s output into VCFs and score them against the truth VCF. Each module's docstring has
the full rules.

```bash
cd machine_learning
G=/scratch/jshen/data/pansoma_v2_tensors/graph_index; PON=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs
$P -m pansoma_net_v2.graph_vcf --predictions <run>/test_chr1/<name>.SNV.predictions.ndjson.gz --tensors $T \
    --kind SNV --graph-index $G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite --output out/HG008T_Illumina.SNV.graph.vcf.gz
$P -m pansoma_net_v2.linear_vcf --graph-vcf out/HG008T_Illumina.SNV.graph.vcf.gz \
    --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path --fasta GRCh38_no_alt_analysis_set.fasta \
    --output out/HG008T_Illumina.SNV.linear.vcf.gz
$P -m pansoma_net_v2.vcfeval --calls out/HG008T_Illumina.SNV.linear.vcf.gz \
    --pon $PON/af-only-gnomad.hg38.vcf.gz $PON/Homo_sapiens_assembly38.dbsnp138.vcf.gz $PON/1000g_pon.hg38.vcf.gz \
          $PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz \
    --truth <somatic truth>.vcf.gz --bed <benchmark>.bed --sdf GRCh38.sdf --rtg rtg --output out/eval_bed
```

The graph index is the one the set was built with (`manifest.json` `graph_index.path`), the reference path and FASTA
those of its labels (`labels.manifest.json` `reference_path`, `fasta`); the example is HG008T_Illumina's (HPRC v1.1
d9). INDEL calls use the same command; the PoN script applies its INDEL rule (gnomAD / CoLoRSdb only, AF ≥ 0.01 in
homopolymers of ≥ 7 bases, else AF ≥ 0.001; the homopolymer from the calls' `##reference` FASTA)
to them.

`machine_learning/pansoma_net_v2/jobs/vcf.sh RUN KIND [OUT]` runs the three on a run's chr1 predictions into
`<run>/vcf_chr1/` (or OUT), with the PoN for both kinds. The set is `pansoma_v2_tensors/<set>/tensors`, `<set>` the
prefix of the predictions file, for a set of any graph (d9, d46, linear); the truth VCF, BED, FASTA and reference path
come from the set's `labels.manifest.json`, the graph index from its `manifest.json`. A predictions file whose prefix
is no set directory stops ("no tensor set"); `SET_NAME=<set>` picks one set when the run predicted several.
`THRESHOLD=t` (a number, or `predict`) overrides the PASS threshold (needed for checkpoints without validation
thresholds). vcfeval runs twice there:

- `eval_bed`: inside the benchmark BED.
- `eval_nobed`: on the whole chromosome (reusing the PoN-tagged calls through `--pon-vcf`). This is the form
  of the earlier rtg runs of ClairS-TO, DeepSomatic and Pansoma v1 (`--squash-ploidy --sample ALT,ALT`, no BED).

**graph_vcf.** One record per tensor: the site's representative allele A1, the allele the model scores.

- Coordinates: CHROM is the node ID, POS is 1-based on the node's forward strand, and ID is the candidate ID.
- INDEL padding uses the node base before the event. At node offset 0 it uses the base after the event,
  or N (flag NOPAD) when the deletion ends at a node end.
- INFO holds the exact event (START, KIND, EVREF, EVALT, PATH), p_somatic / p_germline / p_non, and the site.
- FORMAT holds A1's read support: GT (always 0/1), DP, AD and AF.
- QUAL = −10 log10(1 − p_somatic). rtg rounds scores to 3 decimals, so p near 1 would collapse on a
  linear scale.
- FILTER:
  - PASS: p_somatic ≥ the checkpoint's validation threshold at recall 0.9 (`--target-recall`; the checkpoint
    is the one in predict's `.metrics.json`). The model favours recall and the PoN raises precision.
    `--threshold T` sets the threshold. `--threshold predict` keeps `predict`'s own call (best-F1 threshold).
    A record with predict's off-reference rescue is called by offref_call instead; its QUAL is from
    p_offref_somatic and INFO OFFREF_P_* keeps the rescored probabilities. INDEL predictions have no rescue.
  - LowQual: not called.
  - LowAF: AF below `--min-af`, by default the AF floor of the checkpoint's training labels (their optional
    `snv_min_af` / `indel_min_af`; none when the training labels had no floor, as on every current set, and then no
    record is LowAF). A relabel of the predicted set, or a new sample without labels, does not change it.
- Tensors with p_somatic < `--min-score` (0.01) that are not called are left out. The LowQual records give
  vcfeval its curve.

**linear_vcf.** Each event is projected with `tensor_postprocessing`'s `ReferencePath.linear`, the projection
merge stored as the summaries' `grch38`.

- The projected REF must equal the FASTA, or the run stops.
- INDELs are left-aligned on the FASTA and padded from it.
- Records of one GRCh38 allele are merged (INFO MERGED, MERGED_ALT; ID lists every candidate). The same
  insertion is often seen at the end of one node and at the start of the next.
- Events with no GRCh38 position go unchanged to `<name>.linear.unplaced.vcf.gz`, in graph coordinates:
  CHROM is the node ID, INFO START the 0-based start offset on the node, and ID `node:start:KIND:REF>ALT`.
  These are off-reference nodes, and deletions over nodes that are not GRCh38 neighbours. Placing them would
  need the graph's edges.

**vcfeval.**

- **PoN** (optional: `--pon` or `--pon-vcf`). `scripts/filter_panel_of_normals.py` with the four PoNs. It matches
  by allele, with the rule of each record's kind: SNV gnomAD / CoLoRSdb at AF ≥ 0.0001, dbSNP only records not
  flagged somatic, 1000G; INDEL gnomAD / CoLoRSdb only, at AF ≥ 0.01 in homopolymers of ≥ 7 bases and AF ≥ 0.001
  otherwise (the repository README, section 5). vcfeval passes the calls' `##reference` FASTA (or `--fasta`).
  It also runs over the truth records, which shows how much recall the PoN leaves to any caller.
- **Truth.** The PASS truth records of the kind on the predicted chromosomes.
- **rtg vcfeval.** Always with `--squash-ploidy --sample ALT,ALT -f QUAL`. For the calls without the PoN (raw)
  and, with a PoN, with it (pon), it runs twice:
  - once on the PASS records: the result at the threshold;
  - once with `--all-records` on PASS + LowQual: the curve (best F1, precision at recall 0.5–0.95, ceiling).
- **Reports.** `report.json` and `report.txt`. Each set's row has the result at the threshold and the curve's best
  F1 with its p_somatic, P, R, TP, FP and FN: a choice on the test chromosomes, recorded for reference, never used as
  the threshold. report.txt also gives the PASS calls without a GRCh38 position: they are not evaluated, so precision
  is also shown with them counted as false calls.

HG008T chr1 with the recipes (`analysis/hg008t_snv_indel_7fc7396_20261012`, whose `tables.md` also has the raw and
no-BED rows and the training; SNV Illumina from its earlier run `HG008_Illumina_SNV_h90`, the same recipe and set):
in the benchmark BED after the PoN of each kind, PASS at the checkpoint's validation recall-0.9 threshold
(`val.somatic_at_recall`). Best F1 is the best point of the chr1 curve, with its p_somatic: a test-set choice,
recorded for reference.

| | truth | PASS threshold | PASS evaluated (unplaced) | PoN TP / FP | PoN P / R / F1 | chr1 best F1 at p |
|---|---|---|---|---|---|---|
| SNV Illumina | 697 | 0.419 | 5,556 (570) | 573 / 666 | 0.463 / 0.822 / 0.592 | 0.612 at 0.563 |
| SNV PacBio | 697 | 0.314 | 5,540 (211) | 577 / 584 | 0.497 / 0.828 / 0.621 | 0.641 at 0.395 |
| SNV ONT | 697 | 0.449 | 5,843 (1,409) | 576 / 661 | 0.466 / 0.826 / 0.596 | 0.620 at 0.612 |
| INDEL Illumina | 587 | 0.119 | 3,596 (563) | 140 / 434 | 0.244 / 0.238 / 0.241 | 0.264 at 0.209 |
| INDEL PacBio | 587 | 0.040 | 11,759 (2,075) | 150 / 2,630 | 0.054 / 0.256 / 0.089 | 0.260 at 0.259 |
| INDEL ONT | 587 | 0.018 | 41,052 (7,288) | 134 / 17,268 | 0.008 / 0.228 / 0.015 | 0.140 at 0.168 |

- **Agreement with the truth-level metrics.** For SNV, the best point of rtg's raw curve is that of
  `metrics.truth_report` in predict's `.metrics.json` (the same p_somatic; TP / FP: PacBio 271 / 780 in both,
  Illumina 288 / 954 against 288 / 953, ONT 259 / 783 against 258 / 779; that analysis's `summary.json`). For
  INDEL the truth report counts the partial matches the INDEL test keeps (above) and finds more (Illumina TP 156
  against rtg's 134); on runs whose test left them out the TP agree (`HG008_<P>_INDEL_small` in the same file: 152 /
  106 / 70 against 152 / 107 / 70).
- **Partial labels.** Both recipes train without the partial-match labels (`train --ignore-reasons
  residual_partial_somatic_truth allele_partial_somatic_truth`). No partial SNV tensor is an rtg true positive:
  its allele is an SNV and its truth an INDEL.
- **Pansoma v1** (e053) with the SNV PoN rule on HG008T Illumina chr1: best F1 0.587, against 0.616 for a v2 SNV run of
  the earlier tensor build, on the scorer of `analysis/v1_vs_v2_20260929` (somatic BED ∩ germline BED;
  `synthesis/snv_pon1e4.out`).
- **The PoN on INDELs.** The SNV rule would tag 544 of the 635 chr1 INDEL truth records
  (`analysis/tensor_recall_20260930`), so after it no caller could recall more than 14 %. The INDEL rule tags 413
  (rtg's recall ceiling after it: 0.24–0.26 in the BED); AF ≥ 0.01 for every INDEL tags 383 and removes fewer false
  calls (HG008T Illumina PASS F1 0.218 against 0.246 on the same calls). For SNV the SNV rule tags 61 of 702.
- **Re-scoring existing calls.** `machine_learning/pansoma_net_v2/jobs/vcfeval_only.sh RUN KIND VCFDIR OUT` runs the
  PoN and rtg steps only, on an existing `<run>/<VCFDIR>/*.linear.vcf.gz`, into OUT (a full path; `<run>/vcf_chr1`
  re-scores in place): calls made with an earlier PoN rule get this checkout's, and runs predicted on an earlier build
  of the set can still be scored (`graph_vcf` checks every candidate ID against the set's
  `chr1_variant_summary.ndjson`, so `vcf.sh` stops on them).

## GPU runs (measured 2026-09-28, node tequila)

`--gres=gpu:24gb:1` is an H100 NVL MIG slice (`2g.24gb`); `gpu:h100` is a whole H100 NVL (94 GB). Measured bf16
training speed of the default size, then with the v1 blocks (199.4 M parameters; with the v2 blocks 199.6 M), a whole
H100, batch 256 unless noted:

| setting | speed | GPU memory (peak) |
|---|---|---|
| MIG 2g.24gb, eager, batch 64 / 96 | 122 / 126 tensors/s | 13.7 / 18.7 GiB (128 does not fit) |
| H100, eager, batch 64 / 256 / 512 | 409 / 478 / 485 tensors/s | 13.6 / 44.1 / 84.6 GiB |
| H100, + channels_last + fused AdamW | 627 tensors/s | 44.8 GiB |
| H100, + torch.compile (the default) | 1,404 tensors/s | 26.9 GiB |
| 2 H100s, DDP, eager, batch 384 per GPU | 926 tensors/s | 65 GiB per GPU, summed RSS 47–50 GiB |

- `cudnn.benchmark` changed nothing. Above batch 256 the speed stays flat, so a larger batch only changes
  the optimization.
- Compile takes ~2 min. In fp32 without TF32 the compiled channels_last model gives the eager model's
  logits and gradients (`tests/test_gpu.py`). In bf16 the logits differ by ≤ 0.04, and the loss over 20
  steps follows the same curve.
- Triton builds a helper with `-lcuda`. tequila has only `libcuda.so.1`, so `env.triton_libcuda()` (called
  by `train`) points `TRITON_LIBCUDA_PATH` at a private `libcuda.so` link. Without it compile fails
  whenever that helper is not in Triton's cache.

## Data and index

`data.KindIndex` reads each `<set>/<KIND>/`:

- `manifest.json` (and from it the channel-6 storage `tensor_storage`, and the build's `parameters.haplotypes`; the
  older spelling, `tensor_storage_version` with a `-v1` suffix, is still read) and the label files `*_labels.npy`;
- from every summary record: `candidate_id`, `node_id`, `shard_file`, `index_within_shard`, `row_groups`,
  and the numbers behind `SCALARS` (the record's own fields, not those of its `alleles[]` entries).

It caches the result in `--cache-dir` (default `<output>/index_cache`) as
`<sample>.<set>.<KIND>.<path hash>.npz` plus `.candidates.txt`, and rebuilds it when `labels.manifest.json`
changes (format, version, rules_sha256, created, tensors, the AF floors). From `<chrom>_labels.ndjson` (the
summary's order, checked by candidate id) it also keeps each tensor's reason, off-reference flag, confident
region flag, GRCh38 position (the anchor's middle off the reference) and truth ids, and derives the
evaluation labels. Jobs can share a cache: each process writes its own temporary files, then
renames them. Reading the summaries takes about a minute per
2.5 M tensors.

**Channel 6 storage.** The path-count channel has two codes: `int8-count-haplotypes100` (what indexed_gam_pipeline_v4
builds, every set on disk: relative to `--haplotypes` H, H → 100, one code per missing haplotype down to 51) and
`int8-count-linear100-log2` (the count up to 100, then 100 + ceil(log2(count − 99)): the HPRC v1.1 d9 sets built before
`--haplotypes` existed, since deleted; kept for the checkpoints trained before 2026-10-01, which read only them).
A model reads one of them: `train` refuses
sets that differ (`data.common_path_count`) and records it in the checkpoint (`path_count_storage`,
`path_count_haplotypes`); `train --resume` and `predict` refuse sets coded otherwise than the checkpoint, before
any prediction (a checkpoint without `path_count_storage` counts as `int8-count-linear100-log2`, the storage it was
trained on). A model trained with `--drop-planes path_count` never reads channel 6 (it is the only plane made from
it): it trains and predicts on sets of any storage, and its `path_count_storage` is null. Sets of the same storage
built with another `--haplotypes` are predicted, with a note, and `.metrics.json` records `path_count` (the set's
storage and H, the checkpoint's H).

Loaders are persistent and their workers' start is retried (`train.retry_workers`, up to 3 times). On
tequila, three of four jobs started together once lost their forkserver workers' semaphores in /dev/shm
(FileNotFoundError in SemLock._rebuild). Running workers hold the semaphores, so only the start is exposed.

Validation, prediction and the statistics read each tensor with one `pread` of 161.6 KB in DataLoader workers.
There is no memory map: mapped pages stay in every worker's RSS, and Slurm's summed-RSS limit counts them once
per worker (a first GPU run with 14 workers and memory maps passed 48 GB within 200 steps).

### Reading the training tensors (`chunks.py`)

Random reads cannot keep an H100 busy (~1,400 tensors/s compiled). Measured on tequila on 2026-09-28:

| read | speed |
|---|---|
| one random 161.6 KB `pread` from BeeGFS, one thread, under load | 38 tensors/s (~26 ms per read) |
| a whole 5.3 GB shard in one sequential read, one thread | ~1,290 tensors/s (210 MB/s) |
| training with one random read per tensor in 48 DataLoader workers (the reader before the chunks; four jobs on the node) | 300–550 tensors/s, GPU busy ~20 % |
| chunks of 512 rows, 16 threads, raw shards (other jobs running) | ~350 tensors/s: all readers of the node together got ~300 MB/s |

The tensors are 45 % zero bytes; zstd level 1 packs a 512-row chunk (83 MB) 28.8× and unpacks it at ~1.1 GB/s
per thread. So the training shards are copied once into a chunk cache: per shard one `.zst` file of compressed
chunks and an `.idx.npy` of their offsets, named by the shard's path, size, mtime and chunk size.
`train` builds missing ones (rank 0) before training; or ahead of time:

```bash
$P -m pansoma_net_v2.chunks --tensors $T --kinds SNV INDEL --index-cache <runs>/index_cache \
    --cache <runs>/index_cache/chunks          # HG008 Illumina: 148 SNV shards, ~440 GB -> ~15 GB, ~15 min
```

Each epoch `chunks.ChunkLoader` takes the epoch's tensors from the `EpochSampler` (so `--non-fraction` and
`--keep-non-af` are unchanged), groups them by chunk, reads the chunks in a random order with a thread pool, and
draws the batches from a shuffled pool of 32 chunks at a time, i.e. from ~32 random places of the genome. With
DDP the chunks are dealt to the ranks by load and every rank yields the same number of tensors.

## Model

`model.PansomaNetV2` = `TensorEncoder` → 1×1 front → `ConvNeXtCBAM` → 3 logits. The backbone comes from v1's
`mynet.ConvNeXtCBAMClassifier` (state-dict names kept), so v2 has no import from v1. Weights use the backbone's
own truncated-normal initialization; v1's training script re-initialized them with Kaiming.

`--block` chooses the blocks (config key `block`):

- `v2` (default): ConvNeXt(-V2) as published.
  - GRN takes each channel's L2 norm over the positions and divides it by the mean over the channels, with
    per-channel gamma / beta.
  - A block returns shortcut + branch, so every stage has an identity path.
  - The stem and the downsampling layers have no GELU.
- `v1`: mynet's blocks, for loading and comparing the checkpoints made before 2026-09-29 (a config without
  `block` loads as v1).
  - Its GRN normalizes each read row by the norm over (positions, channels), with one scalar gamma / beta.
  - A GELU follows every residual sum (no identity path through the 39 blocks), the stem and every downsampling.

AdamW decays every weight matrix but not biases, norm weights and the GRN gamma / beta; `--decay-all` (kept for
comparison with the older checkpoints) decays every parameter, as the runs before 2026-09-29 did.

## Tests

```bash
cd machine_learning && $P -m unittest discover -s pansoma_net_v2/tests -t .
```

78 tests; `test_gpu` (4) needs a CUDA device, `test_real_data` (4) the real tensor set, and the rtg case of `test_vcf`
(1) rtg, otherwise they are skipped.

- `test_encode`: every plane is 0 on padding; plane names and shape; one-hot groups; masked z-score (MAPQ 0 ≠
  padding, BQ −1 → `bq_missing`); `differs`; row blocks fill exactly the covered cells of each read; statistics only
  from valid covered cells, their round trip and the unfitted encoder; dropped planes; bf16 output; per-dataset
  statistics: a dataset with BQ and MAPQ doubled and shifted encodes like the original with its own statistics,
  per sample in a mixed batch, path count on the encoder's statistics, no table = the buffers, not in the state dict.
- `test_data`:
  - the index (node, scalars) against a synthetic merged set;
  - the channel-6 storage of a set (`tensor_storage`, the `-v1` spelling, refusal of an unknown one) and of
    several sets together;
  - the cache is reused and rebuilt when the labels change, and two samples' `tensors` sets share one cache
    without collisions;
  - −1 is never selected for training; evaluation adds only the off-reference no-match tensors and drops
    labels 1 and 2 outside the region;
  - truth ids from the labels, and the validation truth alleles of the blocks;
  - the block split keeps whole blocks and is deterministic;
  - `EpochSampler` gives disjoint evaluation slices and equal padded training slices per rank; with a non fraction
    < 1 it takes the kept non tensors every epoch and samples the rest.
- `test_model_train`:
  - forward and checkpoint round trip, with and without scalars, and with a dropped plane;
  - blocks v1 / v2: a v2 block with a zero branch is the identity, GRN v2 against its formula, checkpoints without
    `block` load as v1 with the same logits; the weight-decay groups;
  - loss without −1; the class weights (balanced, sqrt, pow:E); `--micro-batches 2` gives the whole batch's update;
  - forkserver workers; a worker start is retried, then raised;
  - CPU runs: train 2 epochs (chr1 left out, block validation, threshold stored, one BQ / MAPQ statistics row per
    training dataset), resume (statistics kept, the per-dataset ones too), then predict chr1 with the stored
    threshold and statistics fitted on the predicted dataset (the decision rules, the SNV off-reference rescue);
    `--norm train` and a checkpoint without `norm` predict with the training statistics, the old one as before,
    and `predict --norm dataset` fits the dataset's for it; and the truth
    report recomputed from the predictions; an `--ignore-reasons` run (training and validation counts, predict's
    ignored records with label, reason and truth ids kept); runs with `--scalars` and `--keep-non-af`;
  - `--exclude-candidates`: -1 in place, and out of training, validation and the test (predict's default is the
    checkpoint's file; `--exclude-candidates` without a file excludes none);
  - predict's default test: the reasons trained without stay out of the SNV test, the INDEL test keeps the partial
    somatic matches; an explicit `--ignore-reasons` still applies;
  - channel-6 storage: training on mixed storages and `--resume` across storages refused, the checkpoint
    records the storage, predict refuses another storage before any output (a checkpoint without the record
    counts as `int8-count-linear100-log2`), the rescue floor is 100 on `int8-count-haplotypes100` and 90 on
    `int8-count-linear100-log2`.
- `test_real_data`: runs if `PANSOMA_TEST_TENSORS` exists (default: COLO829T Illumina `tensors`); ~70 s.
  Checks the index counts against the manifests (and the region flag against the labeller's reasons, the truth
  ids against `somatic.recall.tsv`), candidates against a summary, and that the row blocks agree
  with the tensors (A1 rows carry the A1 base at the site column, REF rows the graph base, OTHER rows no
  site allele, no reads after the blocks). It also checks the scalars against a summary and encodes real
  tensors.
- `test_chunks`: every tensor of an epoch once with its own blocks, scalars and label, a new order each epoch;
  equal shares for DDP ranks; the zstd cache gives the same tensors, is reused, and a cache of another chunk
  size is refused; read errors reach the training loop.
- `test_metrics`: average precision (ties together, as sklearn), the best-F1 threshold, and the thresholded
  call; the precision and threshold at a recall against brute force; the truth report against brute force,
  duplicates and other-kind truth counted once.
- `test_vcf`: a small graph (a GRCh38 walk, an off-reference node) and one tensor per case.
  - graph_vcf: node padding (offset 0, node end, NOPAD, several nodes); FILTER from the checkpoint's validation
    threshold at a recall, `--threshold` or predict's calls, and the AF floor of the checkpoint's training
    labels (not the predicted set's), with the errors for a checkpoint without them; `--min-score` and
    `--min-af`; the off-reference rescue called by offref_call; no record for excluded tensors; a summary out of
    line; an event that does not fit the graph.
  - linear_vcf: projection equal to `ReferencePath.linear`, left-alignment, a merged junction insertion
    (PASS kept, A1 reads summed), unplaced events, a FASTA that disagrees.
  - the PoN script: the rules of each kind and the homopolymer run.
  - vcfeval: record kinds, summary parsing (no baseline, NaN), the curve points. With rtg (`$RTG` or on PATH): PoN
    tags on calls and truth (a sites-only truth without ##contig), calls vs curve, `--pon-vcf` reuse and refusal,
    no PoN.
- `test_gpu`: runs on a CUDA node. The GPU encoding equals the CPU one, a bf16 training step learns,
  compile + channels_last gives the eager logits and gradients in fp32, and the compiled model has no graph break
  and does not recompile.
