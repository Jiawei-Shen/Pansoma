# E12b: multi-cancer models against the single-tumor models

Reviewer points: **R1.7** (training diversity; retrained models, held-out lines, an HCC1395 reference), **R1.m1**,
**R1.6**. Follows E12a (single-tumor models applied across tumors).

**Status (2026-10-19 13:40): the numbers below are the first predictions (predict.py before f344bca); every ONT chr1
test is invalid (see Caveats).** With the fix (f344bca: predict reuses the checkpoint's statistics for a training
dataset and fits other datasets on all their chromosomes) all 12 models are being re-predicted on tequila (array
445940 -> `test_chr{1,20}_normfix`, `vcf_*_normfix`; `tmp/reviewer_cpu_20261019/E12b_repredict/jobs.txt`); the tables
will be regenerated from those. Runs trained and called on tequila (2 x H100) 2026-10-16 to 10-19 by another session
(`machine_learning/pansoma_net_v2/jobs/submit_multicancer.sh`, `pansoma_net_v2_runs/jobs/multicancer_20261016/jobs.txt`:
chunk caches 417715-417717, trains and vcf jobs 417718-417777, all COMPLETED); collected here by `collect.py`
(`tables.md`, `summary.json`).

## Design

All on HPRC v2.1 d46, the SNV and INDEL recipes of the 7fc7396 benchmark (SNV 50.4 M, 20 epochs, pow 0.75; INDEL small
model; `--norm dataset`, the default since 2026-10-13), one seed:

| group | training sets (chromosomes) | tests |
|---|---|---|
| multicancer4 | HCC1395T (without chr1 and chr20), HCC1954T, H1437T, H2009T (all autosomes) | chr1 of HCC1395T, HG008T, COLO829T; HCC1395T chr20 |
| multicancer6 | multicancer4 + HG008T and COLO829T (without chr1) | the same |

- HCC1395T chr1 / chr20 are held-out chromosomes of a training tumor (chr20: ClairS-TO's held-out chromosome).
- For multicancer4, HG008T and COLO829T are tumors never seen in training (held-out lines).
- `vcf.sh` per test set (graph VCF, GRCh38 VCF, PoN, rtg vcfeval in the BED); INDEL with the current PoN rule (truth
  tagged 418 of 635 on HG008T, checked).
- Baselines: the single-tumor d46 model of the same tumor (in-sample; CPU re-predictions
  `reviewer_cpu_20261013/<run>/vcf_chr1/eval_bed`, INDEL re-scored with the current rule) and of the other tumor
  (cross-tumor; E12a `summary.json`).

## Results, SNV (chr1, BED, after the PoN, PASS at the training sets' recall-0.9 validation threshold)

| test | platform | multicancer4 P / R / F1 | multicancer6 P / R / F1 | single-tumor, in-sample | single-tumor, other tumor |
|---|---|---|---|---|---|
| HCC1395T chr1 | Illumina | 0.886 / 0.792 / **0.837** | 0.888 / 0.793 / **0.838** | - | 0.790 (HG008T), 0.269 (COLO829T) |
| HCC1395T chr1 | PacBio | 0.864 / 0.752 / **0.804** | 0.879 / 0.719 / **0.791** | - | 0.698 (HG008T), 0.741 (COLO829T fiber-seq) |
| HCC1395T chr20 | Illumina / PacBio | 0.794 / 0.766 | 0.796 / 0.769 | - | - |
| HG008T chr1 | Illumina | 0.560 / 0.666 / **0.608** | 0.603 / 0.623 / **0.613** | 0.429 / 0.852 / **0.570** | 0.359 (COLO829T) |
| HG008T chr1 | PacBio | 0.546 / 0.703 / **0.615** | 0.582 / 0.608 / **0.595** | 0.472 / 0.835 / **0.603** | 0.557 (COLO829T fiber-seq) |
| COLO829T chr1 | Illumina | 0.825 / 0.527 / **0.643** | 0.886 / 0.639 / **0.742** | 0.867 / 0.714 / **0.783** | 0.708 (HG008T) |
| COLO829T chr1 | fiber-seq | 0.886 / 0.447 / **0.594** | 0.865 / 0.701 / **0.774** | 0.866 / 0.707 / **0.779** | 0.754 (HG008T PacBio) |

- **HCC1395T** (held-out chromosome of a training tumor): 0.837 / 0.804, +0.047 / +0.063 over the best single-tumor
  model of another tumor; chr20 0.794 / 0.766. Adding HG008T and COLO829T changes it by +0.001 / −0.013.
- **HG008T, never seen (multicancer4):** F1 0.608 / 0.615, slightly above HG008T's own model (0.570 / 0.603), with
  higher precision and lower recall (0.666 / 0.703 against 0.852 / 0.835); far above the COLO829T model (0.359 / 0.557).
- **COLO829T, never seen (multicancer4):** 0.643 / 0.594, below both its own model (0.783 / 0.779) and the HG008T model
  (0.708 / 0.754), through recall (0.527 / 0.447). None of the four training lines has COLO829T's UV-dominated spectrum
  (E18: COLO829T truth closest to SBS7a, the CASTLE lines to SBS4 / SBS40a) [interpretation, not tested].
  multicancer6, which includes COLO829T chr2-22, reaches 0.742 / 0.774, close to the single-tumor model.

## Results, INDEL (best F1 of the chr1 curve, an upper bound; PASS F1 in brackets)

| test | platform | multicancer4 | multicancer6 | single-tumor, in-sample | single-tumor, other tumor |
|---|---|---|---|---|---|
| HCC1395T chr1 | Illumina | 0.662 (0.563) | 0.664 (0.425) | - | 0.257 (HG008T), 0.337 (COLO829T) |
| HCC1395T chr1 | PacBio | 0.548 (0.313) | 0.468 (0.191) | - | - |
| HG008T chr1 | Illumina | 0.223 (0.186) | 0.269 (0.262) | 0.256 (0.253) | 0.270 (COLO829T) |
| HG008T chr1 | PacBio | 0.092 (0.080) | 0.241 (0.154) | 0.275 (0.098) | - |
| COLO829T chr1 | Illumina | 0.171 (0.162) | 0.178 (0.162) | 0.168 (0.117) | 0.141 (HG008T) |
| COLO829T chr1 | fiber-seq | 0.065 (0.044) | 0.155 (0.103) | 0.137 (0.021) | - |

- HCC1395T INDELs are scored under rule (a) (unmatched repeat INDELs ignored, `analysis/truth_labels.md` R13) and its
  SEQC2 truth has few homopolymer INDELs, so its numbers are not comparable with HG008T / COLO829T.
- On HG008T PacBio and COLO829T fiber-seq the unseen-tumor model is far lower (0.092 / 0.065); the model with them
  in training is close to the single-tumor one.

## Caveats

- **`--norm dataset` at prediction (all ONT chr1 rows invalid).** `predict.py` fits each test set's BQ / MAPQ statistics
  on the predicted chromosome only; training fitted them on the training chromosomes. ONT chr1's MAPQ mean is about
  0.5 training SD below chr2-22 (RIS summary `analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/README.md`,
  section "`--norm dataset`: predict fits the statistics on the predicted chromosomes only"), and the ONT chr1 rows
  collapse (SNV F1 0.046-0.492 against 0.674 / 0.700 on HCC1395T chr20, whose statistics are typical). The fix (fit on
  the whole set, labels unused) is a pipeline change pending the user; then re-predict the ONT tests (tequila GPU).
  Illumina and PacBio chr1 MAPQ shifts are −0.37 to +0.01 SD; in-sample `--norm dataset` models moved by −0.006 to
  +0.004 there (RIS summary), so their rows are kept, with this caveat.
- One training run each; the seed spread of single-tumor SNV models is up to ±0.08 PASS F1 (E26 §3).
- The PASS threshold is the multi-tumor validation recall-0.9 threshold; recall-precision trade-offs differ from the
  single-tumor models (higher precision, lower recall on HG008T).
- multicancer4 / multicancer6 HCC1395T, HCC1954T, H1437T, H2009T truths: SEQC2 / CASTLE (R2.2); the CASTLE truths were
  made with DeepSomatic and favour it.
