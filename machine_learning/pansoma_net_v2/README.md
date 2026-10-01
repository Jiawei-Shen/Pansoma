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
   - BQ, MAPQ and path count get a masked z-score. The statistics are fitted once on the covered, valid
     cells of the training tensors and stored in the model (checkpoint `stats`); they are never refitted on
     other data. Valid cells are BQ ≥ 0, MAPQ ≥ 0 and path count > 0; the z-score is 0 elsewhere.
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
# (the SNV recipe: see "SNV recipe" below)
$P -m pansoma_net_v2.train --tensors $T --kinds SNV --output runs/HG008_Illumina_SNV \
    --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75 \
    --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
torchrun --nproc_per_node=2 -m pansoma_net_v2.train --ddp ...           # several GPUs, same options

# test: every chr1 tensor (-1 included), the checkpoint's statistics and somatic threshold
$P -m pansoma_net_v2.predict --checkpoint runs/HG008_Illumina_SNV/best.pth --tensors $T --kinds SNV \
    --chroms chr1 --output runs/HG008_Illumina_SNV/test_chr1
```

`/scratch/jshen/data/pansoma_net_v2_runs/jobs/run.sh NAME SET KIND [train options]` runs both steps as one
Slurm job. `python -m pansoma_net_v2.summary [NAME ...]` prints each run's per-epoch table, its current step and, when
done, its chr1 test result (per tensor and against the truth VCF). The table columns are:

- training loss, speed, minutes, GPU peak;
- validation per tensor: somatic AP and s.F1 / s.P / s.R at the best threshold, and the threshold t;
- against the truth VCF at the same t: t.AP, t.F1 / t.P / t.R and the ceiling (below);
- germline F1.

Useful `train` options:

- **Split.**
  - `--test-chroms` (default chr1).
  - `--val-fraction 0.05` of the `--val-block-nodes 20000`-node blocks of the training chromosomes; whole
    blocks, so validation tensors share no reads with training tensors. `--val-chroms` instead validates on
    whole chromosomes.
- **Imbalance.**
  - `--non-fraction F`: each epoch takes every somatic and germline tensor and a fresh random fraction F of
    the non tensors.
  - `--keep-non-af A` (with `--non-fraction`): the non tensors with AF ≥ A are taken every epoch. The sampled
    non tensors get loss weight 1/F, so each epoch's loss stays an unbiased estimate of the loss over all
    non tensors, whatever their AF.
  - `--class-weights balanced|sqrt|pow:E|none|w0,w1,w2`: balanced is n / (3 n_c) over one epoch's tensors;
    sqrt is its square root, pow:E its power E. The loss is normalized by the batch's summed weights, so
    scaling all three weights changes nothing.
- **Label reasons left out.** `--ignore-reasons R ...` treats the tensors of those label reasons as -1 in training and
  validation; the label files are unchanged. HG008T Illumina SNV runs leave out
  `residual_partial_somatic_truth allele_partial_somatic_truth`: SNV tensors that partially match an INDEL truth
  (mostly read-end SNV shadows of the INDEL) taught the model low-AF artefacts next to germline INDELs, and as SNV
  calls they never match the truth (analysis/v1_vs_v2_20260929).
- **Model.**
  - `--scalars` feeds the site scalars to the head (`data.SCALARS`: log coverage, log site coverage, log
    ALT / REF / OTHER counts, AF, second allele AF, allele count, log event length). They are the counts
    before the 200-row cap, z-scored with training statistics that are stored in the model.
  - `--block v2|v1` (below, Model) and `--decay-all`.
  - `--depths`, `--dims`, `--front`, `--drop-path` set the size (defaults `3 3 27 3` / `192 384 768 1536`,
    as in v1's training script).
- **Selection.** `--select f1|ap|truth_f1` picks the best checkpoint by the thresholded somatic F1 per tensor
  (default), the somatic average precision, or the best F1 against the truth VCF (below). With `truth_f1` the
  stored threshold is the one of that F1.
- **Speed.** `--no-compile` and `--no-channels-last` turn off the defaults on CUDA (below). Training reads
  chunks from a zstd cache (`--read-threads 16`, `--chunk-rows 512`, `--window-chunks 32`, `--chunk-cache`;
  "Reading the training tensors" below); `--no-chunk-cache` reads the chunks from the raw shards, and
  `--random-reads` goes back to one pread per tensor in DataLoader workers (`--num-workers`).
- **Memory.** `--micro-batches K` runs the forward and backward of each batch in K parts (gradient accumulation):
  the same update as the whole batch, with the loss normalized over the whole batch, and the GPU memory of one part.
  Batch 1024 of the default model fills a 94 GB H100 NVL (about 94.7 of 95.8 GB); batch 512 takes about 49 GB.
- Also: `--epoch-samples`, `--val-samples`, `--stats-samples 20000`, `--amp bf16|off`, `--resume` (keeps the
  checkpoint's statistics), `--seed`.

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

On HG008 Illumina chr1:

| | truth alleles | with a tensor (ceiling) | tensors of the other kind's truth |
|---|---|---|---|
| SNV | 697 SNP | 662 (95.0 %) | 188 SNV tensors of 137 INDEL truth alleles |
| INDEL | 588 DEL + INS | 317 (53.9 %) | 1 |

Of those 137 INDEL truth alleles, 74 also have an INDEL tensor. Adding the two models' truth reports would
count them twice, so `combine` scores a sample's SNV and INDEL predictions together against the whole truth
VCF, as a pipeline's output VCF is compared with it:

```bash
$P -m pansoma_net_v2.combine HG008_Illumina_SNV_base HG008_Illumina_INDEL_base --output HG008_Illumina.json
```

- The truth alleles are those of all predicted kinds.
- A truth allele is found when any of its test tensors, in either set, is called by its model (each model's
  own threshold). It counts once, however many calls it has.
- Precision, recall, F1, the ceiling, per truth kind, which sets found each truth allele, and the calls
  repeated on a truth already found. No AP: the two models' scores are not on one scale.
- With both sets, the INDEL truth alleles on chr1 with a tensor are 380 (64.6 %) instead of 317.

Outputs in `--output`:

- `args.json` (the options) and `train.log`.
- `metrics.jsonl` (per epoch: train loss and accuracy, speed, GPU peak memory, and the validation report).
- `stats.json`.
- `last.pth`, and `best.pth` (best validation by `--select`).

A checkpoint holds `format: pansoma_net_v2`, `config`, `model_state_dict` (the encoder's and scalar
statistics are buffers in it), `stats`, `somatic_threshold`, `classes`, `planes`, `scalars`, the training
data (directories, label provenance, counts), the chromosome split and `args`. Rebuild it with
`PansomaNetV2.from_checkpoint(path)`.

`predict` writes, for every tensor of the chosen chromosomes, `<sample>.<set>.<KIND>.predictions.ndjson.gz`.
Each record has chrom, candidate_id, label (truth), test_label (the scored label, null when left out),
in_test, reason, off_reference (the node is off the GRCh38 path, as in the paper's pangenome-only calls),
truth_ids (the somatic truth alleles of a labelled-1 tensor), p_non, p_somatic, p_germline, and pred (with
the threshold). It also writes `.metrics.json`: the same report over the test tensors, with the truth report
and the counts of test tensors on off-reference nodes (`--threshold` overrides t). `--ignore-reasons` (default: the
checkpoint's training `--ignore-reasons`) leaves those reasons out of the test like -1; their records keep label,
reason and truth_ids and say `ignored: true`.
Off-reference tensors are marked for analysis only; the test scores them like the others.

**Off-reference rescue (SNV).** The model learned that a low path count (ch6: how many HPRC haplotypes run
through the node) is not somatic, and on an off-reference node every read has a low path count, so it hardly ever
calls somatic there. `predict --offref-site90 W` (default 10; 0 turns it off) scores every SNV off-reference tensor
again with its ch6 cells below 90 set to 90 in the site columns (the middle column ± W; uncovered cells stay 0), and
calls it somatic when somatic is the most probable class of those probabilities. Its record adds p_offref_non,
p_offref_somatic, p_offref_germline and offref_call. `.metrics.json` adds `offref_rescue` (tensors, calls before
and after) and the rule `pipeline`: offref_call on the off-reference tensors, the checkpoint's validation recall-0.9
threshold on the others.

**SNV recipe** (HG008T Illumina, PacBio, ONT): the 50.4 M model `--depths 3 3 27 3 --dims 96 192 384 768` (the
default widths halved), `--class-weights pow:0.75 --batch-size 1024 --lr 1e-4 --epochs 20 --non-fraction 1.0
--block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth`, every plane kept (ch6
included), one whole GPU; calls at the validation recall-0.9 threshold, off-reference tensors by the rescue above
(graph_vcf's defaults). On HG008T PacBio chr1 (12 epochs) the default 199.6 M model at lr 1e-4 overfits after
epoch ~7 (validation loss 0.51 -> 2.16) and the 50.4 M model does not; with the PoN in the BED both reach P/R/F1
0.473/0.832/0.604 (199.6 M, best epoch 5) and 0.493/0.822/0.616 (50.4 M).

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

INDEL calls go through vcfeval without `--pon` for now.

`/scratch/jshen/data/pansoma_net_v2_runs/jobs/vcf.sh RUN KIND` runs the three on a run's chr1 predictions into
`<run>/vcf_chr1/`, with the truth VCF and BED of the set's labels. It uses the PoN for SNV only. For a checkpoint
without validation thresholds (trained before `val.somatic_at_recall`), set `THRESHOLD=t`. vcfeval runs twice there:

- `eval_bed`: inside the benchmark BED.
- `eval_nobed`: on the whole chromosome (SNV: reusing the PoN-tagged calls through `--pon-vcf`). This is the form
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
  - LowAF: below the AF floor of the checkpoint's training labels (Illumina SNV 0.07; none for INDEL or the
    long-read sets). A relabel of the predicted set, or a new sample without labels, does not change it.
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

- **PoN** (optional: `--pon` or `--pon-vcf`; SNV only for now). `scripts/filter_panel_of_normals.py` with the
  four PoNs. All match by allele; gnomAD and CoLoRSdb only at AF ≥ 0.0001, dbSNP only records not flagged somatic.
  It also runs over the truth records, which shows how much recall the PoN leaves to any caller.
- **Truth.** The PASS truth records of the kind on the predicted chromosomes.
- **rtg vcfeval.** Always with `--squash-ploidy --sample ALT,ALT -f QUAL`. For the calls without the PoN (raw)
  and, with a PoN, with it (pon), it runs twice:
  - once on the PASS records: the result at the threshold;
  - once with `--all-records` on PASS + LowQual: the curve (best F1, precision at recall 0.5–0.95, ceiling).
- **Reports.** `report.json` and `report.txt`. The latter also gives the PASS calls without a GRCh38 position:
  they are not evaluated, so precision is also shown with them counted as false calls.

HG008 Illumina chr1 (2026-09-30; SNV `HG008_Illumina_SNV_allnon_sqrt_b1024`, `<run>/vcf_chr1/`; INDEL
`HG008_Illumina_INDEL_base`, `<run>/vcf_chr1_r09/`). PASS is at the validation recall-0.9 threshold: SNV p ≥ 0.254
(the checkpoint's `val.somatic_at_recall`), INDEL p ≥ 0.034 (`THRESHOLD`, from the run's saved validation
probabilities `val_probs_e2.npz`; its checkpoint predates `val.somatic_at_recall`). INDEL without the PoN. Best F1
is the best point of the chr1 curve, a test-set choice.

| | PASS evaluated (unplaced) | TP | FP | P | R | F1 | best F1 (P, R) |
|---|---|---|---|---|---|---|---|
| SNV raw, in BED (697) | 5,487 (8) | 598 | 4,889 | 0.109 | 0.858 | 0.193 | 0.266 (0.222, 0.333) |
| SNV pon, in BED | 1,212 | 575 | 637 | 0.474 | 0.825 | 0.602 | 0.615 (0.515, 0.763) |
| SNV pon, no BED (702) | 1,259 | 579 | 680 | 0.460 | 0.825 | 0.591 | 0.606 (0.502, 0.763) |
| INDEL raw, in BED (587) | 4,396 (1,281) | 217 | 4,179 | 0.049 | 0.370 | 0.087 | 0.150 (0.106, 0.256) |
| INDEL raw, no BED (635) | 4,929 (1,281) | 234 | 4,695 | 0.048 | 0.368 | 0.084 | 0.143 (0.100, 0.252) |

At predict's best-F1 threshold (`--threshold predict`, `<run>/vcf_chr1/`) INDEL_base gives, raw and no BED:
899 PASS evaluated (324 unplaced), P 0.109, R 0.154, F1 0.128.

- **Agreement with the truth-level metrics** (`HG008_Illumina_SNV_keephard_w100` and INDEL_base at predict's
  threshold). The results agree with `metrics.truth_report`:
  - SNV raw: TP 315 in both.
  - FP: 1,617 against 1,627. rtg does not score the 17 unplaced false calls, counts 4 SNV tensors next to
    INS truths as FP, and 3 calls lie in the somatic BED but outside the germline BED.
  - INDEL TP: 92 against 123. The other 31 are unplaced residual-partial tensors.
  - SNV pon, in BED, of the table: TP 575 and FP 637 here, TP 574 and FP 635 with the scorer of
    `analysis/v1_vs_v2_20260929` (same PoN rule, somatic BED ∩ germline BED).
- **Partial labels.** The SNV run is trained without the SNV partial-match labels (`train --ignore-reasons
  residual_partial_somatic_truth allele_partial_somatic_truth`). No partial SNV tensor is an rtg true positive:
  its allele is an SNV and its truth an INDEL.
- **Pansoma v1** (e053) with the same PoN rule on that scorer: best F1 0.587, against 0.616 for this SNV run.
- **The PoN on INDELs.** It tags 544 of the 635 chr1 INDEL truth records, CoLoRSdb 533 of them, most in
  homopolymers at population AF around 0.1–1 %. So after it, no caller can recall more than 14 %. For SNV it tags
  61 of 702.

## GPU runs (measured 2026-09-28, node tequila)

`--gres=gpu:24gb:1` is an H100 NVL MIG slice (`2g.24gb`); `gpu:h100` is a whole H100 NVL (94 GB). The default
model has 199.4 M parameters. Measured bf16 training speed (a whole H100, batch 256 unless noted):

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

- `manifest.json` and the label files `*_labels.npy`;
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
| training with 48 workers and random reads (four jobs on the node) | 300–550 tensors/s, GPU busy ~20 % |
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

AdamW decays every weight matrix but not biases, norm weights and the GRN gamma / beta; `--decay-all` decays
every parameter, as the runs before 2026-09-29 did.

## Tests

```bash
cd machine_learning && $P -m unittest discover -s pansoma_net_v2/tests -t .
```

- `test_encode`: every plane is 0 on padding; one-hot groups; masked z-score (MAPQ 0 ≠ padding, BQ −1 →
  `bq_missing`); `differs`; row blocks fill exactly the covered cells of each read; statistics only from
  valid covered cells; bf16 output.
- `test_data`:
  - the index (node, scalars) against a synthetic merged set;
  - the cache is reused and rebuilt when the labels change, and two samples' `tensors` sets share one cache
    without collisions;
  - −1 is never selected for training; evaluation adds only the off-reference no-match tensors and drops
    labels 1 and 2 outside the region;
  - truth ids from the labels, and the validation truth alleles of the blocks;
  - the block split keeps whole blocks and is deterministic;
  - `EpochSampler` gives disjoint evaluation slices and equal padded training slices per rank.
- `test_model_train`:
  - forward and checkpoint round trip, with and without scalars;
  - blocks v1 / v2: a v2 block with a zero branch is the identity, GRN v2 against its formula, checkpoints without
    `block` load as v1 with the same logits; the weight-decay groups;
  - loss without −1;
  - forkserver workers;
  - CPU runs: train 2 epochs (chr1 left out, block validation, threshold stored), resume (statistics kept),
    then predict chr1 with the stored threshold, and the truth report recomputed from the predictions; runs
    with `--scalars`, `--select truth_f1` and `--keep-non-af`.
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
  call; the truth report against brute force, duplicates and other-kind truth counted once.
- `test_combine`: one truth found from both sets counts once; repeated calls; other-kind truth with one set;
  two predictions of one kind are refused.
- `test_vcf`: a small graph (a GRCh38 walk, an off-reference node) and one tensor per case.
  - graph_vcf: node padding (offset 0, node end, NOPAD, several nodes); FILTER from the checkpoint's validation
    threshold at a recall, `--threshold` or predict's calls, and the AF floor of the checkpoint's training
    labels (not the predicted set's), with the errors for a checkpoint without them; a summary out of line.
  - linear_vcf: projection equal to `ReferencePath.linear`, left-alignment, a merged junction insertion
    (PASS kept, A1 reads summed), unplaced events, a FASTA that disagrees.
  - vcfeval: summary parsing (no baseline, NaN), the curve points. With rtg (`$RTG` or on PATH): PoN tags on
    calls and truth (a sites-only truth without ##contig), calls vs curve, `--pon-vcf` reuse and refusal, no PoN.
- `test_gpu`: runs on a CUDA node. The GPU encoding equals the CPU one, a bf16 training step learns, and
  compile + channels_last gives the eager logits and gradients in fp32.
