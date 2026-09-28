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

```bash
cd machine_learning
P=/wanglab/jshen/anaconda3/bin/python       # torch 2.8 + timm 1.0 (GPU)
T=/scratch/jshen/data/pansoma_v2_tensors

# train: chromosomes split train / --val-chroms / --test-chroms (left out); -1 tensors are not used
$P -m pansoma_net_v2.train --tensors $T/COLO829T_Illumina/v3_tensors $T/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors \
    --output runs/illumina --val-chroms chr20 --test-chroms chr21 chr22 --epochs 30 --epoch-samples 1000000
torchrun --nproc_per_node=4 -m pansoma_net_v2.train --ddp ...            # several GPUs, same options

# predict: every tensor (-1 included) with the checkpoint's statistics; metrics over the labelled ones
$P -m pansoma_net_v2.predict --checkpoint runs/illumina/best.pth --tensors $T/COLO829T_Illumina/v3_tensors \
    --chroms chr21 chr22 --output runs/illumina/test
```

Useful `train` options:

- `--kinds SNV INDEL` (default both).
- `--epoch-samples` / `--val-samples`: a fresh random subset of the training tensors each epoch, and a
  fixed subset of the validation tensors.
- `--class-weights balanced|none|w0,w1,w2`: balanced is n / (3 n_c) over the training tensors.
- `--depths`, `--dims`, `--front`, `--drop-path`: model size (defaults `3 3 27 3` / `192 384 768 1536`, as
  in v1's training script).
- `--stats-samples 20000`, `--amp bf16|off`, `--resume` (keeps the checkpoint's statistics), `--seed`.

Outputs in `--output`:

- `train.log` and `metrics.jsonl` (per epoch: train loss and accuracy, speed, validation loss, accuracy,
  per-class precision / recall / F1, confusion matrix).
- `stats.json`.
- `last.pth`, and `best.pth` (the best validation somatic F1).

A checkpoint holds `format: pansoma_net_v2`, `config`, `model_state_dict` (the encoder's statistics are
buffers in it), `stats`, `classes`, `planes`, the training data (directories, label versions, counts), the
chromosome split and `args`. Rebuild it with `PansomaNetV2.from_checkpoint(path)`.

`predict` writes `<sample>.<set>.<KIND>.predictions.tsv.gz` (chrom, candidate_id, label, p_non, p_somatic,
p_germline, pred) and `.metrics.json`.

## GPU runs (measured 2026-09-28)

On the `gpu` partition (node tequila), `--gres=gpu:24gb:1` is an H100 NVL MIG slice (`2g.24gb`). The default
model has 199.4 M parameters. A bf16 training step on real tensors measured:

| batch | GPU memory (peak) | speed |
|---|---|---|
| 32 | 8.6 GiB | 111 tensors/s |
| 64 | 13.7 GiB | 122 tensors/s |
| 96 | 18.7 GiB | 126 tensors/s |
| 128 | out of memory | |

A short run on COLO829T Illumina (2 epochs of 20,000 tensors, batch 64, 14 workers) trained at ~120
tensors/s. The job's summed RSS was 24.9 GiB (main process 8.3 GB). At that speed one pass over COLO829T
Illumina's 1.67 M labelled training tensors takes ~3.9 h, so use `--epoch-samples`, the full H100s
(`gpu:h100`), torchrun over several GPUs, or a smaller `--dims`.

```bash
sbatch -p gpu --gres=gpu:24gb:1 -c 16 --mem=30G -t 24:00:00 -o train-%j.out --wrap \
  "cd /scratch/jshen/Github/Pansoma/machine_learning && $P -m pansoma_net_v2.train --tensors ... --output ... --num-workers 14"
```

## Data and index

`data.KindIndex` reads each `<set>/<KIND>/`: `manifest.json`, the label files `*_labels.npy`, and from every
summary record `candidate_id`, `shard_file`, `index_within_shard`, `af` and `row_groups`. It caches the
result in `--cache-dir` (default `<output>/index_cache`) and rebuilds it when `labels.manifest.json`
changes (version, created, tensors, snv_min_af). Reading the summaries takes about a minute per 2.5 M
tensors (COLO829T Illumina: 2.65 M tensors in ~60 s, 1.7 GB).

Each sample is one `pread` of 161.6 KB at its row's offset in the shard. There is no memory map: mapped
pages stay in every DataLoader worker's RSS, and Slurm's summed-RSS limit counts them once per worker (a
first GPU run with 14 workers and memory maps passed 48 GB within 200 steps).

## Model

`model.PansomaNetV2` = `TensorEncoder` → 1×1 front → `ConvNeXtCBAM` → 3 logits. The backbone is v1's
`mynet.ConvNeXtCBAMClassifier`, copied unchanged (state-dict names included) apart from its print, so v2 has
no import from v1. Weights use the backbone's own truncated-normal initialization; v1's training script
re-initialized them with Kaiming.

## Tests

```bash
cd machine_learning && $P -m unittest discover -s pansoma_net_v2/tests -t .
```

- `test_encode`: every plane is 0 on padding; one-hot groups; masked z-score (MAPQ 0 ≠ padding, BQ −1 →
  `bq_missing`); `differs`; row blocks fill exactly the covered cells of each read; statistics only from
  valid covered cells; bf16 output.
- `test_data`: the index against a synthetic merged set; the cache is reused, and rebuilt when the labels
  change; −1 is never selected for training; `EpochSampler` gives disjoint evaluation slices and equal
  padded training slices per rank.
- `test_model_train`: forward, checkpoint round trip, loss without −1, and a CPU run of train (2 epochs),
  resume (statistics kept), then predict (every tensor written).
- `test_real_data`: runs if `PANSOMA_TEST_TENSORS` exists (default: COLO829T Illumina `v3_tensors`); ~70 s.
  Checks the index counts against the manifests, candidates against a summary, and that the row blocks agree
  with the tensors (A1 rows carry the A1 base at the site column, REF rows the graph base, OTHER rows no
  site allele, no reads after the blocks). It also encodes real tensors.
- `test_gpu`: runs on a CUDA node. The GPU encoding equals the CPU one, and a bf16 training step learns.
