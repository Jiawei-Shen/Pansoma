# E12b: multi-cancer models against the single-tumor models

Reviewer points: **R1.7** (training diversity; retrained models, held-out lines, an HCC1395 reference), **R1.m1**,
**R1.6**. Follows E12a (single-tumor models applied across tumors).

**Status (2026-10-20): all numbers are the re-predictions with the `--norm dataset` fix** (f344bca: predict reuses the
checkpoint's statistics for a training dataset and fits other datasets on all their chromosomes; tequila array 445940
-> `test_chr{1,20}_normfix`, vcf jobs 445941-445990 -> `vcf_*_normfix`, all 48 reports present;
`tmp/reviewer_cpu_20261019/E12b_repredict/jobs.txt`). The first predictions are kept in `summary.json` as
`multicancer{4,6}_prefix`. ClairS-TO and DeepSomatic tumor-only were run on the HCC1395T test data on tequila
(2026-10-19/20, `tmp/reviewer_cpu_20261019/E8_hcc1395/`, scored with E8 `score_comparator.py`) and are now columns of
`tables.md`, next to the re-scored existing runs on HG008T / COLO829T (E8). Runs trained and called on tequila (2 x H100)
2026-10-16 to 10-19 by another session (`machine_learning/pansoma_net_v2/jobs/submit_multicancer.sh`,
`pansoma_net_v2_runs/jobs/multicancer_20261016/jobs.txt`: chunk caches 417715-417717, trains and vcf jobs
417718-417777, all COMPLETED); collected here by `collect.py` (`tables.md`, `summary.json`).

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

## Results, SNV (chr1 and HCC1395T chr20, BED, after the PoN, PASS at the training sets' recall-0.9 validation threshold)

P / R / **F1**; cross-tumor single-tumor for HCC1395T = HG008T model / COLO829T model. Comparators: PASS after the same
PoN; DeepSomatic's models were trained on HCC1395 with chr1 held out, ClairS-TO's never saw HCC1395 (E8b D3, C3).

| test set | platform | multicancer4 (target tumor not in training, except HCC1395T) | multicancer6 | in-sample single-tumor | cross-tumor single-tumor | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|---|---|---|
| HCC1395T chr1 | Illumina | 0.877 / 0.816 / **0.845** | 0.879 / 0.822 / **0.850** | - | 0.801 / 0.778 / **0.790** / 0.934 / 0.157 / **0.269** | 0.781 / 0.927 / **0.848** | 0.806 / 0.859 / **0.832** |
| HCC1395T chr20 | Illumina | 0.834 / 0.764 / **0.797** | 0.839 / 0.781 / **0.809** | - | - | 0.769 / 0.937 / **0.845** | 0.649 / 0.835 / **0.730** |
| HCC1395T chr1 | PacBio | 0.856 / 0.787 / **0.820** | 0.843 / 0.809 / **0.825** | - | 0.792 / 0.624 / **0.698** / 0.868 / 0.646 / **0.741** | 0.766 / 0.927 / **0.839** | 0.805 / 0.822 / **0.814** |
| HCC1395T chr20 | PacBio | 0.803 / 0.750 / **0.775** | 0.804 / 0.773 / **0.788** | - | - | 0.740 / 0.957 / **0.835** | 0.778 / 0.812 / **0.795** |
| HCC1395T chr1 | ONT | 0.650 / 0.668 / **0.658** | 0.767 / 0.697 / **0.730** | - | 0.763 / 0.628 / **0.689** / 0.783 / 0.515 / **0.621** | 0.755 / 0.766 / **0.760** | 0.791 / 0.730 / **0.760** |
| HCC1395T chr20 | ONT | 0.721 / 0.661 / **0.690** | 0.755 / 0.691 / **0.721** | - | - | 0.717 / 0.791 / **0.752** | 0.759 / 0.736 / **0.747** |
| HG008T chr1 | Illumina | 0.543 / 0.684 / **0.605** | 0.592 / 0.611 / **0.601** | 0.429 / 0.852 / **0.570** | 0.603 / 0.255 / **0.359** | 0.376 / 0.889 / **0.529** | 0.410 / 0.845 / **0.552** |
| HG008T chr1 | PacBio | 0.535 / 0.720 / **0.614** | 0.586 / 0.670 / **0.625** | 0.472 / 0.835 / **0.603** | 0.567 / 0.548 / **0.557** | 0.425 / 0.806 / **0.557** | 0.455 / 0.844 / **0.591** |
| HG008T chr1 | ONT | 0.492 / 0.726 / **0.587** | 0.635 / 0.551 / **0.590** | 0.432 / 0.859 / **0.575** | 0.579 / 0.420 / **0.487** | 0.310 / 0.869 / **0.457** | 0.414 / 0.831 / **0.552** |
| COLO829T chr1 | Illumina | 0.827 / 0.564 / **0.671** | 0.880 / 0.661 / **0.755** | 0.867 / 0.714 / **0.783** | 0.792 / 0.640 / **0.708** | 0.427 / 0.835 / **0.565** | 0.431 / 0.764 / **0.551** |
| COLO829T chr1 | PacBio | 0.892 / 0.446 / **0.595** | 0.880 / 0.691 / **0.774** | 0.866 / 0.707 / **0.779** | 0.817 / 0.701 / **0.754** | 0.732 / 0.555 / **0.632** | 0.779 / 0.745 / **0.762** |
| COLO829T chr1 | ONT | 0.497 / 0.308 / **0.380** | 0.817 / 0.670 / **0.737** | 0.787 / 0.692 / **0.736** | 0.678 / 0.589 / **0.630** | 0.582 / 0.695 / **0.633** | 0.749 / 0.684 / **0.715** |

- **HCC1395T** (held-out chromosomes of a training tumor): multicancer4 / multicancer6 chr1 0.845 / 0.850 (Illumina),
  0.820 / 0.825 (PacBio), 0.658 / 0.730 (ONT), above every single-tumor model of another line (at most 0.790 / 0.741 /
  0.689). DeepSomatic-TO, trained on HCC1395 except chr1, reaches 0.848 / 0.839 / 0.760 and ClairS-TO 0.832 / 0.814 /
  0.760: equal on Illumina, DeepSomatic-TO 0.014 higher on PacBio, both 0.030 higher on ONT (multicancer6). On chr20
  DeepSomatic-TO leads (0.845 / 0.835 / 0.752 against 0.809 / 0.788 / 0.721).
- **HG008T, never seen (multicancer4):** 0.605 / 0.614 / 0.587, above HG008T's own models (0.570 / 0.603 / 0.575) and
  both comparators (DeepSomatic-TO 0.529 / 0.557 / 0.457, ClairS-TO 0.552 / 0.591 / 0.552), with higher precision and
  lower recall than the own models.
- **COLO829T, never seen (multicancer4):** 0.671 / 0.595 / 0.380, through recall (0.564 / 0.446 / 0.308); none of the
  four training lines has COLO829T's UV-dominated spectrum (E18) [interpretation, not tested]. With COLO829T chr2-22 in
  training (multicancer6): 0.755 / 0.774 / 0.737, against its own models 0.783 / 0.779 / 0.736 and the comparators
  0.565 / 0.632 / 0.633 (DeepSomatic-TO) and 0.551 / 0.762 / 0.715 (ClairS-TO).
- The fix changed the Illumina and PacBio rows by at most 0.03 (e.g. HCC1395T chr1 0.837 -> 0.845) and restored ONT
  (first predictions 0.046-0.492).

## Results, INDEL (best F1 of the chr1 curve, an upper bound; PASS F1 in brackets)

| test set | platform | multicancer4 (target tumor not in training, except HCC1395T) | multicancer6 | in-sample single-tumor | cross-tumor single-tumor | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|---|---|---|
| HCC1395T chr1 | Illumina | **0.662** (0.565) | **0.690** (0.484) | - | **0.257** (0.176) / **0.337** (0.314) | **0.672** (0.605) | **0.415** (0.641) |
| HCC1395T chr20 | Illumina | **0.699** (0.525) | **0.741** (0.453) | - | - | **0.828** (0.624) | **0.408** (0.565) |
| HCC1395T chr1 | PacBio | **0.537** (0.311) | **0.471** (0.187) | - | - / - | **0.628** (0.612) | **0.498** (0.566) |
| HCC1395T chr20 | PacBio | **0.635** (0.315) | **0.538** (0.182) | - | - | **0.796** (0.733) | **0.457** (0.548) |
| HCC1395T chr1 | ONT | **0.255** (0.013) | **0.173** (0.007) | - | - / - | **0.519** (0.471) | **0.534** (0.569) |
| HCC1395T chr20 | ONT | **0.342** (0.016) | **0.267** (0.009) | - | - | **0.805** (0.707) | **0.522** (0.558) |
| HG008T chr1 | Illumina | **0.222** (0.184) | **0.272** (0.262) | **0.256** (0.253) | **0.270** (0.262) | **0.294** (0.208) | **0.236** (0.132) |
| HG008T chr1 | PacBio | **0.090** (0.079) | **0.247** (0.119) | **0.275** (0.098) | - | **0.248** (0.053) | **0.196** (0.129) |
| HG008T chr1 | ONT | **0.101** (0.047) | **0.126** (0.014) | **0.159** (0.015) | - | **0.154** (0.034) | **0.155** (0.127) |
| COLO829T chr1 | Illumina | **0.170** (0.162) | **0.181** (0.164) | **0.168** (0.117) | **0.141** (0.081) | **0.200** (0.191) | **0.148** (0.122) |
| COLO829T chr1 | PacBio | **0.064** (0.042) | **0.158** (0.098) | **0.137** (0.021) | - | **0.174** (0.133) | **0.169** (0.158) |
| COLO829T chr1 | ONT | **0.009** (0.002) | **0.056** (0.005) | **0.089** (0.002) | - | **0.091** (0.021) | **0.117** (0.103) |

- HCC1395T INDELs are scored under rule (a) (unmatched repeat INDELs ignored, `analysis/truth_labels.md` R13) and its
  SEQC2 truth has few homopolymer INDELs, so its numbers are not comparable with HG008T / COLO829T.
- On HCC1395T DeepSomatic-TO leads on every platform (best F1 0.672 / 0.628 / 0.519 against multicancer4 0.662 / 0.537 /
  0.255), most on ONT, where Pansoma's PASS F1 is near zero (0.013). On HG008T and COLO829T all callers are low
  (best F1 at most 0.294) and close.
- On HG008T PacBio and COLO829T fiber-seq the unseen-tumor model is far lower (0.090 / 0.064); the model with them in
  training is close to the single-tumor one.

## Caveats

- **`--norm dataset` at prediction.** The first predictions fitted each test set's BQ / MAPQ statistics on the predicted
  chromosome only, which broke every ONT chr1 row (ONT chr1's MAPQ mean is about 0.5 training SD below chr2-22). All
  numbers above are from the re-predictions after the fix (f344bca); the first ones are in `summary.json` (`*_prefix`).
- One training run each; the seed spread of single-tumor SNV models is up to ±0.08 PASS F1 (E26 §3).
- The PASS threshold is the multi-tumor validation recall-0.9 threshold; recall-precision trade-offs differ from the
  single-tumor models (higher precision, lower recall on HG008T).
- multicancer4 / multicancer6 HCC1395T, HCC1954T, H1437T, H2009T truths: SEQC2 / CASTLE (R2.2); the CASTLE truths were
  made with DeepSomatic and favour it.
