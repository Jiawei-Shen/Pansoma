# Fact sheet v3: run identity for R1.1 (channel 6), R3.9a / R3.9b, chr2 INDEL confirmation, R3.1 split (2026-10-30)

Tags: **EXISTS** (number already in a repo / run file), **COMPUTED NOW** (computed today from existing files; scripts and
outputs in `SP/runid/`), **MISSING** (no such run or file locally). No new training, prediction or Slurm job was run.

Abbreviations:
- `SP` = `/tmp/claude-10008/-scratch-jshen-Github-Pansoma/2a27ec06-fc25-4d6b-b0b5-2fbc7483d580/scratchpad`
- `RIS` = `/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/<run>/` (local copy of the RIS H100 runs:
  `args.json`, `train.log`, `test_chr1*/…metrics.json`, sometimes `vcf_chr1*/eval_bed/report.json`)
- `CPU` = `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/`: chr1 re-prediction of the RIS
  checkpoint on CPU in fp32 (`tmp/reviewer_cpu_20261013/predict_cpu.sh`: `--amp off`, `best.pth` symlinked to RIS)
- `RISSUM` = `analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/{summary.json,README.md,tables.md}`
  (`summary.json` `seeds` / `no_ch6` = RIS GPU originals)
- `PJ` = `analysis/benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012/summary.json` `earlier` (the
  COLO829T `*_SNV_h90` runs; `runs/COLO829_*_SNV_h90/train.log`)
- Collected run identities: `SP/runid/collect_runs.py` -> `SP/runid/runs.tsv`, `runs.json` (COMPUTED NOW, from args.json
  `seed/norm/ddp/batch_size/drop_planes`, train.log `device cuda x N` and `done; best epoch N (f1 …)`, report.json,
  metrics.json `somatic_ap`).

Common to every SNV run below (args.json, EXISTS): ConvNeXt 50.4 M, 20 epochs, lr 1e-4, wd 0.05, `pow:0.75`,
`--batch-size 1024` **per GPU** (`train.py` line 99 help "per GPU"; `ris/snv_pipeline.sh` header "global batch N x
1024"), val_fraction 0.05, val_block_nodes 20,000, test chr1. Checkpoint rule: epoch with the highest validation
pre-PoN per-tensor somatic F1 (`--best-by f1`, the default; args.json has no `best_by` key = older default f1;
`train.py` line 24). PASS = `p_somatic >=` that checkpoint's validation recall-0.9 threshold. "norm train-stats" = run
trained before `--norm` existed (`args.json "norm": null`): BQ / MAPQ / path count z-scored with training statistics,
equivalent to `--norm train` (RISSUM README lines 17-21; group 0 bit-identical check).

---------------------------------------------------------------------------------------------------------------------

## Q1. Channel-6 table of R1.1: identity of every cell

### 1.1 What the letter's table uses (letter `SP/ver2_parts/R1.1.md` table; COMPUTED NOW `SP/runid/runs.tsv`)

| dataset | graph | arm | run | seed | norm | GPUs x per-GPU batch (global) | epoch (val F1) | PASS thr | P / R / F1 (letter) | AP | single / mean |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HG008T Illumina | v1.1 | with | `HG008_Illumina_SNV_7fc7396_ddp4_b1024` | 0 | train-stats | 4 x 1024 (4,096) | 4 (0.2851) | 0.251 | 0.452 / 0.831 / 0.585 | 0.199 | single |
| HG008T Illumina | v1.1 | without | `HG008_Illumina_SNV_nopc_7fc7396_ddp4_b1024` | 0 | train | 4 x 1024 (4,096) | 3 (0.1467) | 0.178 | 0.407 / 0.809 / 0.542 | 0.075 | single |
| HG008T PacBio | v1.1 | with | `HG008_PacBio_SNV_7fc7396` | 0 | train-stats | **1 x 1024 (1,024)** | 7 (0.3298) | 0.314 | 0.497 / 0.828 / 0.621 | 0.233 | single |
| HG008T PacBio | v1.1 | without | `HG008_PacBio_SNV_nopc_7fc7396` | 0 | train | 1 x 1024 (1,024) | 6 (0.1362) | 0.283 | 0.456 / 0.798 / 0.580 | 0.078 | single |
| HG008T ONT | v1.1 | with | `HG008_ONT_SNV_7fc7396` | 0 | train-stats | **1 x 1024** | 8 (0.3009) | 0.449 | 0.466 / 0.826 / 0.596 | 0.232 | single |
| HG008T ONT | v1.1 | without | `HG008_ONT_SNV_nopc_7fc7396` | 0 | train | 1 x 1024 | 6 (0.1370) | 0.138 | 0.424 / 0.834 / 0.562 | 0.065 | single |
| COLO829T Illumina | v1.1 | with | `COLO829_Illumina_SNV_h90` (2026-10-03, code 9d6bd8d) | 0 | train-stats | **1 x 1024** | 6 (0.6871) | 0.330 | 0.860 / 0.711 / 0.779 (exact 0.7785 = 3476/4465 -> **0.778**) | 0.725 | single |
| COLO829T Illumina | v1.1 | without | `COLO829_Illumina_SNV_nopc_7fc7396` | 0 | train | 1 x 1024 | **2** (0.5209) | 0.182 | 0.850 / 0.687 / 0.760 | 0.544 | single |
| COLO829T Fiber-seq | v1.1 | with | `COLO829_fiberseq_SNV_h90` (2026-10-03) | 0 | train-stats | **1 x 1024** | 12 (0.6917) | 0.394 | 0.857 / 0.717 / 0.781 | blank in letter; **0.724 EXISTS** (PJ `earlier."fiberseq SNV".test.somatic_ap` 0.7243) | single |
| COLO829T Fiber-seq | v1.1 | without | `COLO829_fiberseq_SNV_nopc_7fc7396` | 0 | train | 1 x 1024 | 9 (0.5240) | 0.099 | 0.846 / 0.715 / 0.775 | 0.575 | single |
| COLO829T ONT | v1.1 | with | `COLO829_ONT_SNV_h90` (2026-10-03) | 0 | train-stats | **1 x 1024** | 8 (0.6779) | 0.463 | 0.727 / 0.686 / 0.706 | blank in letter; **0.667 EXISTS** (PJ, 0.6667) | single |
| COLO829T ONT | v1.1 | without | `COLO829_ONT_SNV_nopc_7fc7396` | 0 | train | 1 x 1024 | 5 (0.5047) | 0.157 | 0.713 / 0.675 / 0.693 | 0.480 | single |
| all six | v2.1 | with | `<S>_<P>_SNV_d46_7fc7396{,_s1,_s2}` | 0, 1, 2 | train-stats / train / train | 4 x 1024 (4,096) | see 1.3 | see 1.3 | **three-seed mean** | **range over 3 seeds** | **mean** |
| all six | v2.1 | without | `<S>_<P>_SNV_nopc_d46_7fc7396` | 0 | train | 4 x 1024 (4,096) | 3 / 6 / 5 / 6 / 7 / 7 | see 1.3 | single seed 0 | single | single |

Sources: v1.1 P/R/F1 and v2.1 no-ch6 = RISSUM `summary.json` `no_ch6[*].d9_with / d9_without / d46_without` (EXISTS);
epochs, seeds, GPUs = args.json + train.log (EXISTS; h90 epochs/GPUs from PJ `earlier.*.training` and
`runs/COLO829_*_SNV_h90/train.log` "device cuda x 1", "done; best epoch 6 / 12 / 8"; h90 seed 0 = default `--seed 0` at
9d6bd8d and `COLO829_Illumina_SNV_h90_at9d6bd8d/args.json` seed 0 reproduces it bit for bit, RISSUM README l.76-78);
AP = `test_chr1*/…SNV.metrics.json` `somatic_ap` (EXISTS).

### 1.2 Matching verdict (COMPUTED NOW from 1.1)

- **v2.1, all six pairs: NOT matched.** "With" = three-seed mean (P/R/F1) and three-seed AP range; "without" = one run,
  seed 0. A seed-matched pair exists for all six (with = `<S>_<P>_SNV_d46_7fc7396`, seed 0, the first seed of Table
  R3.1): rebuilt in 1.3. Remaining difference inside the seed-0 pair: code commit (with = 7fc7396, trained 2026-10-12,
  pre-`--norm`; without = RIS code of 2026-10-16, after f2a0e4e; torch.compile last-bit kernel differences, RISSUM README
  l.79-83). A fully code-matched alternative exists for 5 of 6 (normbq pairs, 1.4).
- Small consequence of the mean: the letter's v2.1 with-means are means of 3-decimal per-seed values and use the fp32 CPU
  re-prediction for seed 0 of HG008T PacBio, so they differ from Table R3.1 by 0.001 in places (HG008T PacBio R 0.826 vs
  0.827, F1 0.606 vs 0.607; HG008T ONT F1 0.590 vs 0.589; COLO829T Fiber-seq F1 0.781 vs 0.780). Gone once rebuilt.
- **v1.1 HG008T Illumina: matched** (4 GPU both, seed 0, same norm, both validation-selected). Only code commit differs.
- **v1.1 HG008T PacBio / ONT and COLO829T x 3: matched within the pair on GPUs (1 / 1), global batch (1,024), seed (0),
  norm and checkpoint rule, but not the main-table v1.1 recipe** (Table R3.1 v1.1 = 4 GPUs, global batch 4,096). A
  4-GPU no-ch6 model for these five is **MISSING** (RISSUM README l.22-24: group 2 kept the 1-GPU slots of
  `gpu_nopc_commands.sh`). The letter's footnote discloses the 1 GPU.
- **v1.1 with-arm code differs from the no-ch6 arm**: HG008T 1-GPU controls = 7fc7396 (they reproduced the 9d6bd8d h90
  runs exactly, RISSUM README l.70); COLO829T controls = h90 (9d6bd8d, 2026-10-03); no-ch6 = 2026-10-16 code. For HG008T
  this does not change results (README l.82: HG008T does not show the compile differences). **For COLO829T Illumina it
  does** (Q2): the same seed/recipe with the current code selects epoch 8 and scores F1 0.624. For COLO829T Fiber-seq /
  ONT no current-code rerun of the control exists (**MISSING**), so their sensitivity is unknown.
- AP column: the v1.1 Fiber-seq / ONT blanks can be filled (0.724 / 0.667, EXISTS).

### 1.3 Rebuilt table from seed-matched pairs, validation-selected checkpoints (COMPUTED NOW; `SP/runid/ch6_rebuild.py` -> `SP/runid/ch6_rebuild.md`)

"With" v2.1 = seed 0 RIS GPU values (= Table R3.9a seed 0; RISSUM `seeds`); AP = RIS `test_chr1` metrics.

| dataset | v1.1 with (run, epoch) | v1.1 without (epoch) | ΔF1 / ΔAP | v2.1 with seed 0 (epoch) | v2.1 without (epoch) | ΔF1 / ΔAP |
|---|---|---|---|---|---|---|
| HG008T Illumina | 0.452 / 0.831 / 0.585 (0.199), ep 4, 4 GPU | 0.407 / 0.809 / 0.542 (0.075), ep 3 | −0.043 / −0.124 | 0.429 / 0.852 / 0.570 (0.310), ep 5 | 0.398 / 0.783 / 0.528 (0.089), ep 3 | −0.042 / −0.221 |
| HG008T PacBio | 0.497 / 0.828 / 0.621 (0.233), ep 7, 1 GPU | 0.456 / 0.798 / 0.580 (0.078), ep 6 | −0.041 / −0.155 | 0.473 / 0.839 / 0.605 (0.353), ep 7 | 0.466 / 0.819 / 0.594 (0.075), ep 6 | −0.011 / −0.278 |
| HG008T ONT | 0.466 / 0.826 / 0.596 (0.232), ep 8, 1 GPU | 0.424 / 0.834 / 0.562 (0.065), ep 6 | −0.034 / −0.167 | 0.432 / 0.859 / 0.575 (0.328), ep 9 | 0.447 / 0.841 / 0.584 (0.071), ep 5 | **+0.009** / −0.258 |
| COLO829T Illumina, control h90 (10-03 code) | 0.860 / 0.711 / 0.778 (0.725), ep 6, 1 GPU | 0.850 / 0.687 / 0.760 (0.544), ep 2 | −0.019 / −0.181 | 0.867 / 0.715 / 0.783 (0.818), ep 7 | 0.362 / 0.678 / 0.472 (0.582), ep 6 | −0.311 / −0.236 |
| COLO829T Illumina, control `_h90_repro` (current code, same seed/recipe, `--norm train`) | 0.554 / 0.713 / 0.624 (0.740), ep 8, 1 GPU | same | **+0.136** / −0.196 | | | |
| (COLO829T Illumina 1-GPU seeds 1 / 2, current code: `_h90_s1` / `_s2`) | 0.649 / 0.704 / 0.675 (0.729), ep 6; 0.682 / 0.692 / 0.687 (0.729), ep 5 | same (seed 0) | +0.085 / +0.073 F1 | | | |
| COLO829T Fiber-seq | 0.857 / 0.717 / 0.781 (0.724), ep 12, 1 GPU | 0.846 / 0.715 / 0.775 (0.575), ep 9 | −0.006 / −0.150 | 0.866 / 0.707 / 0.779 (0.829), ep 9 | 0.883 / 0.685 / 0.771 (0.592), ep 7 | −0.007 / −0.237 |
| COLO829T ONT | 0.727 / 0.686 / 0.706 (0.667), ep 8, 1 GPU | 0.713 / 0.675 / 0.693 (0.480), ep 5 | −0.012 / −0.186 | 0.788 / 0.692 / 0.737 (0.755), ep 10 | 0.703 / 0.649 / 0.675 (0.484), ep 7 | −0.062 / −0.271 |

Derived (COMPUTED NOW):
- With the letter's controls (h90 for COLO829T Illumina): F1 lower without ch6 in **11 / 12** seed-0 pairs, median
  **−0.026** (exact ΔF1: −0.0431, −0.0425, −0.0407, −0.0111, −0.0335, +0.0088, −0.0186, −0.3115, −0.0059, −0.0071,
  −0.0125, −0.0617). The letter's "11 of 12 first-seed pairs, median −0.026" is therefore already the seed-0 count; only
  the table shows means.
- With the code-matched control `_h90_repro` for COLO829T Illumina v1.1: **10 / 12, median −0.023**. The no-ch6
  model (0.760) is above all three current-code 1-GPU seeds of the control (0.624 / 0.675 / 0.687) and below only the
  2026-10-03 run (0.778).
- AP lower without ch6 in **12 / 12** (every variant above, incl. `_h90_repro` 0.740 -> 0.544). Ranges: HG008T
  0.199-0.353 -> 0.065-0.089; COLO829T 0.667-0.829 -> 0.480-0.592 (both graphs; v2.1 only: 0.755-0.829 -> 0.484-0.592).
- v2.1 no-ch6 against the with-ch6 three-seed range (`--norm train`, RISSUM seeds): below the range on 5 / 6; HG008T ONT
  inside (0.584 in 0.575-0.598) (EXISTS, RISSUM README l.95-99).
- no-ch6 − GRCh38-only: vs GRCh38-only seed 0 **+0.004 to +0.292** (12 / 12 positive; smallest COLO829T Illumina v2.1
  0.472 vs 0.468); vs GRCh38-only three-seed mean **+0.011 to +0.299** (letter's numbers); two of twelve inside the
  GRCh38-only seed range (HG008T PacBio v1.1 0.580 in 0.513-0.592; COLO829T Illumina v2.1 0.472 in 0.386-0.530) (EXISTS
  `SP/facts/ablation.md` §5, re-checked).

### 1.4 Secondary, fully code-matched v2.1 pairs: normbq (`--norm dataset`, BQ per dataset), seed 0, 4 GPUs (COMPUTED NOW, `SP/runid/runs.tsv`)

| dataset | with `<S>_<P>_SNV_d46_normbq` (epoch) | without `<S>_<P>_SNV_nopc_d46_normbq` (epoch) | ΔF1 | AP with -> without |
|---|---|---|---|---|
| HG008T Illumina | 0.426 / 0.852 / 0.568 (5) | 0.398 / 0.786 / 0.529 (3) | −0.039 | 0.311 -> 0.088 |
| HG008T PacBio | 0.471 / 0.838 / 0.603 (7) | 0.466 / 0.825 / 0.596 (6) | −0.007 | 0.350 -> 0.075 |
| HG008T ONT | 0.433 / 0.859 / 0.576 (9) | 0.447 / 0.842 / 0.584 (5) | +0.008 | 0.326 -> 0.071 |
| COLO829T Illumina | 0.867 / 0.712 / 0.782 (7) | 0.367 / 0.680 / 0.477 (6) | −0.305 | 0.819 -> 0.582 |
| COLO829T Fiber-seq | 0.863 / 0.707 / 0.777 (9) | 0.764 / 0.704 / 0.732 (8) | −0.045 | 0.829 -> 0.590 |
| COLO829T ONT | 0.814 / 0.700 / 0.753 (9) | **MISSING** (`COLO829_ONT_SNV_nopc_d46_normbq` requeued, not synced) | - | - |

Same direction as `--norm train` (4 of 5 lower; HG008T ONT +0.008 again).

### 1.5 What the AP column is (EXISTS, code)

**Tensor-level, before population filtering, threshold-free.** `machine_learning/pansoma_net_v2/predict.py` line 298
`report = metrics.report(labels, probs, threshold)` -> `metrics.py` lines 102-108: keeps tensors with label >= 0 and
computes `somatic_ap = average_precision(labels == SOMATIC, probs[:, SOMATIC])` (sklearn-style AP, `metrics.py`
lines 32-45) over all labelled chr1 candidate tensors (unlabelled −1 tensors, i.e. outside the labels' confident
region, and the ignored partial-label tensors are excluded; PoN-tagged tensors are included). Written to
`test_chr1*/…SNV.metrics.json` `somatic_ap`, which is what the table quotes. P / R / F1 are call-level (rtg vcfeval of
the PASS calls after the PoN, benchmark BED). Note: because AP excludes label −1 tensors, it cannot see the COLO829T
Illumina v1.1 false calls of Q2 (all on label −1 tensors), RISSUM README l.88-91.

---------------------------------------------------------------------------------------------------------------------

## Q2. COLO829T Illumina, HPRC v1.1: 0.779 (channel-6 table) vs 0.572 (main table)

| | channel-6 table "with" (letter 0.779) | main table seed 0 (R3.1 / R3.9a, 0.572) | code-matched ch6 control |
|---|---|---|---|
| run | `COLO829_Illumina_SNV_h90` | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024` | `COLO829_Illumina_SNV_h90_repro` |
| trained | RIS 2026-10-03, code 9d6bd8d | RIS 2026-10-16 | RIS 2026-10-19, current code |
| seed / norm | 0 / train-stats | 0 / train | 0 / train |
| GPUs x batch (global) | 1 x 1024 (1,024) | 4 x 1024 (4,096) | 1 x 1024 (1,024) |
| checkpoint | **validation-selected epoch 6** (val F1 0.6871 vs 0.6870 at ep 8) | validation-selected epoch 8 (0.6806) | validation-selected epoch 8 (0.6874 vs 0.6858 at ep 6) |
| PASS threshold | 0.3303 | 0.3267 | 0.4549 |
| TP / FP (after PoN) | 1,738 / 284 | 1,705 / 1,815 | 1,742 / 1,401 |
| P / R / F1 | 0.8595 / 0.7114 / 0.7785 | 0.4844 / 0.6979 / 0.5719 | 0.5542 / 0.7131 / 0.6237 |
| chr1 somatic AP | 0.725 | 0.736 | 0.740 |

(EXISTS: PJ `earlier."Illumina SNV"`, RISSUM summary `no_ch6` / `seeds`, RIS `args.json` / `train.log` /
`vcf_chr1*/eval_bed/report.json`; RISSUM README l.40-91.)

- **The 0.779 control is a different run, not the exploratory epoch-6 checkpoint.** It is the 1-GPU run of 2026-10-03,
  whose own validation chose epoch 6 by a 0.0001 margin. The exploratory epoch-6 row of R3.9a is
  `COLO829_Illumina_SNV_d9_normbq_sweep_e06` (4 GPUs, normbq, chr1-chosen, F1 0.773). The letter's 0.779 should read
  0.778 (F1 = 2·1738 / (2·1738 + 284 + 705) = 0.77850).
- Why they differ: different recipe (1 vs 4 GPUs; code) and, decisively, a different selected epoch. From epoch 7 these
  models call 1,000-1,600 unlabelled candidates (label −1, outside the labels' confident region) in chr1:121-126 and
  143-144 Mb; validation (and AP) never sees them, so the epoch choice is a near-tie (RISSUM README l.53-91; epoch sweep
  `epoch_sweep.md`: epochs 1-6 0.772-0.784, epoch 6 0.779, epochs 7-9 0.624-0.650). The 4-GPU seeds all choose epoch 8
  (0.572 / 0.626 / 0.566); the current-code 1-GPU seeds choose 8 / 6 / 5 (0.624 / 0.675 / 0.687). Only the 10-03 run
  reaches 0.778; the 10-03 code reproduces it bit for bit (`_h90_at9d6bd8d`), so it is genuine but code-dependent
  (torch.compile last-bit kernels after f2a0e4e).
- **Validation-selected matched value for the ch6 pair**: `_h90_repro` = 0.554 / 0.713 / 0.624 (AP 0.740), same seed,
  GPUs, norm and code generation as the no-ch6 run (0.850 / 0.687 / 0.760, AP 0.544, epoch 2). On F1 the model without
  channel 6 is then **higher by 0.136**; AP is lower by 0.196. For the main-table comparison the validation-selected
  value stays 0.572 (seed 0) / 0.588 ± 0.033 (three seeds).

---------------------------------------------------------------------------------------------------------------------

## Q3. R3.9a vs R3.9b: same runs?

### 3.1 Run IDs of every model in R3.9a / R3.9b / R3.1 (COMPUTED NOW, `SP/runid/r39_runids.py` -> `SP/runid/r39_runids.md`)

All 4 GPUs x 1024 (global 4,096); checkpoint = validation-selected epoch; PASS = validation recall 0.9.

| dataset | reference | run (seeds 1 / 2: `_s1` / `_s2`) | norm s0 / s1 / s2 | epoch s0 / s1 / s2 | PASS thr s0 / s1 / s2 | R3.9b seed-0 call set |
|---|---|---|---|---|---|---|
| HG008T Illumina | GRCh38-only | `HG008_Illumina_SNV_linear_7fc7396` | train-stats / train / train | 6 / 3 / 3 | 0.064 / 0.129 / 0.092 | CPU fp32 re-prediction |
| HG008T Illumina | v1.1 | `HG008_Illumina_SNV_7fc7396_ddp4_b1024` | train-stats / train / train | 4 / 7 / 5 | 0.251 / 0.071 / 0.317 | CPU fp32 re-prediction |
| HG008T Illumina | v2.1 | `HG008_Illumina_SNV_d46_7fc7396` | train-stats / train / train | 5 / 5 / 6 | 0.324 / 0.554 / 0.340 | CPU fp32 re-prediction |
| HG008T PacBio | GRCh38-only | `HG008_PacBio_SNV_linear_7fc7396` | train-stats / train / train | 7 / 3 / 5 | 0.040 / 0.143 / 0.069 | CPU fp32 re-prediction |
| HG008T PacBio | v1.1 | `HG008_PacBio_SNV_7fc7396_ddp4_b1024` | train / train / train | 8 / 11 / 6 | 0.214 / 0.117 / 0.114 | RIS GPU |
| HG008T PacBio | v2.1 | `HG008_PacBio_SNV_d46_7fc7396` | train-stats / train / train | 7 / 9 / 9 | 0.291 / 0.430 / 0.325 | CPU fp32 re-prediction |
| HG008T ONT | GRCh38-only | `HG008_ONT_SNV_linear_7fc7396` | train-stats / train / train | 8 / 8 / 6 | 0.042 / 0.050 / 0.145 | CPU fp32 re-prediction |
| HG008T ONT | v1.1 | `HG008_ONT_SNV_7fc7396_ddp4_b1024` | train / train / train | 9 / 9 / 8 | 0.328 / 0.116 / 0.245 | RIS GPU |
| HG008T ONT | v2.1 | `HG008_ONT_SNV_d46_7fc7396` | train-stats / train / train | 9 / 8 / 6 | 0.231 / 0.482 / 0.592 | CPU fp32 re-prediction |
| COLO829T Illumina | GRCh38-only | `COLO829_Illumina_SNV_linear_7fc7396` | train-stats / train / train | 7 / 6 / 6 | 0.081 / 0.127 / 0.161 | CPU fp32 re-prediction |
| COLO829T Illumina | v1.1 | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024` | train / train / train | 8 / 8 / 8 | 0.327 / 0.435 / 0.373 | RIS GPU (copied) |
| COLO829T Illumina | v2.1 | `COLO829_Illumina_SNV_d46_7fc7396` | train-stats / train / train | 7 / 6 / 7 | 0.640 / 0.703 / 0.303 | CPU fp32 re-prediction |
| COLO829T Fiber-seq | GRCh38-only | `COLO829_fiberseq_SNV_linear_7fc7396` | train-stats / train / train | 10 / 9 / 7 | 0.085 / 0.166 / 0.129 | CPU fp32 re-prediction |
| COLO829T Fiber-seq | v1.1 | `COLO829_fiberseq_SNV_7fc7396_ddp4_b1024` | train / train / train | 9 / 9 / 9 | 0.364 / 0.306 / 0.431 | RIS GPU (copied) |
| COLO829T Fiber-seq | v2.1 | `COLO829_fiberseq_SNV_d46_7fc7396` | train-stats / train / train | 9 / 8 / 10 | 0.481 / 0.437 / 0.344 | CPU fp32 re-prediction |
| COLO829T ONT | GRCh38-only | `COLO829_ONT_SNV_linear_7fc7396` | train-stats / train / train | 6 / 7 / 7 | 0.092 / 0.112 / 0.152 | CPU fp32 re-prediction |
| COLO829T ONT | v1.1 | `COLO829_ONT_SNV_7fc7396_ddp4_b1024` | train / train / train | 9 / 9 / 6 | 0.335 / 0.407 / 0.425 | RIS GPU (copied) |
| COLO829T ONT | v2.1 | `COLO829_ONT_SNV_d46_7fc7396` | train-stats / train / train | 10 / 10 / 8 | 0.568 / 0.604 / 0.535 | CPU fp32 re-prediction |

R3.9a exploratory row (COLO829T Illumina v1.1, chr1-chosen epochs, normbq, 4 GPUs): `COLO829_Illumina_SNV_d9_normbq_sweep_e06`
(seed 0, thr 0.277, 0.832 / 0.723 / 0.773), `…_d9_normbq_s1_sweep_e05` (seed 1, 0.492, 0.705 / 0.705 / 0.705),
`…_d9_normbq_s2_sweep_e06` (seed 2, 0.336, 0.783 / 0.698 / 0.738). (EXISTS)

R3.9a values = RISSUM `summary.json` `seeds` (RIS GPU predictions, bf16 amp) for all 54 runs (EXISTS; matches
`SP/facts/ablation.md` §1.1 exactly). R3.9b = `tmp/assessment_20261030/main_bootstrap/pairs.tsv` (job 462895; every
letter R3.9b HPRC row equals pairs.tsv rounded, checked by `SP/runid/check_r39b.py` -> `check_r39b.out`). pairs.py reads
`reviewer_cpu_20261013/<run>` first, then RIS (`pairs.py` `BASES`).

### 3.2 Finding (COMPUTED NOW, `SP/runid/check_r39b.out`)

- **Same checkpoints, seeds, epochs, norm and thresholds in both tables; the assessment's normbq hypothesis is wrong.**
  The difference is the prediction run: for the 13 seed-0 call sets marked "CPU fp32 re-prediction" above, R3.9b
  scores the fp32 CPU re-prediction of the same `best.pth` (`predict_cpu.sh`, `--amp off`), while R3.9a uses the RIS
  bf16 GPU predictions. Re-predicted vs GPU (TP / FP): HG008T PacBio v2.1 **582 / 650 vs 585 / 652 (recall 0.835 vs
  0.839)**, HG008T PacBio GRCh38-only 576 / 969 vs 576 / 973, HG008T Illumina GRCh38-only 548 / 1,089 vs 548 / 1,091,
  HG008T Illumina v1.1 579 / 702 vs 579 / 703, COLO829T Illumina GRCh38-only 1,423 / 2,221 vs 1,423 / 2,216, COLO829T
  Illumina v2.1 1,745 / 268 vs 1,746 / 268, COLO829T Fiber-seq GRCh38-only 1,525 / 521 vs 1,526 / 522, COLO829T ONT
  GRCh38-only 1,325 / 2,575 vs 1,327 / 2,575, COLO829T ONT v2.1 1,690 / 457 vs 1,690 / 454; the other four identical.
  Max |ΔF1| per call set 0.0016, |ΔR| 0.0043.
- The 0.835 that the assessment links to a normbq model is this same fp32 re-prediction of the `--norm train` seed-0
  checkpoint (it also appears as HG008T PacBio "own model" 0.472 / 0.835 / 0.603 in Supplementary Table R1.7b, letter
  l.849). The normbq seed-0 model `HG008_PacBio_SNV_d46_normbq` has 0.471 / 0.838 / 0.603 (RIS).
- R3.9b point estimates vs R3.9a seed-0 differences: 16 of 54 entries differ by more than rounding, all by about 0.001
  except three: HG008T PacBio v2.1 − GRCh38-only ΔR **+0.009 vs +0.013**, ΔF1 **+0.090 vs +0.092**; HG008T PacBio
  v2.1 − v1.1 ΔR **+0.006 vs +0.010** (its ΔF1 −0.004 vs −0.0025). (v1.1 − GRCh38-only HG008T PacBio ΔF1 +0.094 vs +0.095 is also a
  re-prediction effect, via GRCh38-only.)
- R3.9b's `Pansoma − comparator` rows use the same seed-0 `*_SNV_d46_7fc7396` fp32 re-predictions (`--norm train`)
  (EXISTS: `experiments/E16_bootstrap/paired.tsv`, A = `<S>_<P>_SNV_d46_7fc7396/vcf_chr1`), not the per-dataset-
  normalised (normbq) three-seed models that R3.9a's caption says R1.6 / R2.3 / R3.5 use for caller comparisons.
- Rebuilding R3.9b on the RIS GPU call sets would need the RIS `vcf_chr1/eval_bed/vcfeval/pon_calls` of the 13 seed-0
  runs (or their `predictions.ndjson.gz`); locally the RIS dirs hold only `test_chr1/*.metrics.json`: **MISSING** (RIS
  sync needed).

### 3.3 Point estimate and interval (EXISTS, code)

Point estimate = difference of P / R / F1 computed from the **whole-chr1 rtg totals** of the two call sets (not a
bootstrap mean or median): `tmp/assessment_20261030/main_bootstrap/pairs.py` `pa = E.metrics(point[a])` with `point` =
block counts summed, asserted equal to report.json `pon.at_threshold`; interval = 2.5th / 97.5th percentiles of the
2,000 paired resample differences (`E16_bootstrap/e16_bootstrap.py` `ci`, lines 140-145; blocks = 1-Mb chr1 windows
overlapping the benchmark BED, 231 HG008T / 249 COLO829T; resample seed 20261015, drawn once per block set so both
members of a pair share them).

---------------------------------------------------------------------------------------------------------------------

## Q4. chr2 confirmation of the final INDEL models: how the overall ΔF1 and CI are built (EXISTS)

Script: `/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/final_20261023/analyze.py --rule ppr8` (defaults
`--resamples 2000 --seed 20261023`) -> `results_ppr8.json` `groups."unseen|confirmation"`, `logs/analyze_ppr8_final.out`
l.16; pre-registered in `PREREGISTRATION.md` §4 (+ Amendment 1, PASS rule ppr8).

- Cells: 6 = {HG008T, COLO829T} chr2 x {Illumina, PacBio (COLO829T Fiber-seq), ONT}; models
  `unseen_s{0,1,2}_{Illumina,PacBio,ONT}_INDEL` (four other tumors in training); comparators DeepSomatic-TO 1.9.0 and
  ClairS-TO 0.4.2 run on chr2 (`/scratch/jshen/data/indel_gap_20261020/confirm_chr2/score`).
- Per cell: Pansoma F1 = **mean of the three seeds' F1** (`analyze.py` l.139; not pooled counts); difference = seed-mean
  F1 − comparator F1 (l.144).
- Bootstrap (`Blocks`, l.55-78): 1-Mb windows of chr2 overlapping the tumor's benchmark BED (242 HG008T, 243 COLO829T;
  COMPUTED NOW), 2,000 resamples with replacement, multiplicities drawn once per (chromosome, BED) with RNG
  `[20261023, crc32("chr2|<bed>")]` (l.65). The BED string is identical for all seeds, platforms and both comparators
  of a tumor (checked NOW), so within a tumor everything is paired; HG008T and COLO829T are resampled independently.
  Per resample: each seed's F1 from resampled rtg counts, seed-mean (l.141), minus comparator (l.145-146).
- Overall: **unweighted mean over the six cells** of the per-cell differences, per resample; 95 % CI = 2.5 / 97.5
  percentiles (l.161, `ci` l.88). Seed variation is not resampled (the three models are fixed).
- Result: vs ClairS-TO **−0.0206 [−0.0388, −0.0008]**, vs DeepSomatic-TO **+0.0344 [+0.0164, +0.0526]**; on par 6 / 6;
  verdict success. Per cell vs ClairS-TO (seed-mean − ClairS-TO): HG008T Illumina −0.023 [−0.044, −0.002], PacBio
  −0.004, ONT −0.005; COLO829T Illumina −0.012, Fiber-seq −0.011, ONT −0.069 [−0.107, −0.028] (below in 6 / 6, interval
  excludes zero in 2). Vs DeepSomatic-TO: above in 5 / 6 (HG008T Illumina −0.102).

---------------------------------------------------------------------------------------------------------------------

## Q5. R3.1 validation split overlap and its effect on training data (EXISTS + COMPUTED NOW from `tmp/assessment_20261030/ablation/q3/val_split.json`)

- Rule (EXISTS, `data.py` 307-322): validation = node blocks (`node_id // 20000`) with SHA-1(`"{seed}:{chrom}:{block}"`)
  < 0.05, on chr2-22; node IDs differ between graphs, so the same seed gives different GRCh38 regions per graph.
- Size: validation = **4.5-6.2 % of the chr2-22 GRCh38 span** (115-158 Mb of 2,560 Mb, tile definition) and **4.0-7.3 %
  of the labelled tensors** (train.log counts, `q3/val_counts.txt`; COMPUTED NOW).
- Overlap between graphs, same seed (54 pairs, tile): intersection = **2.5-7.6 % of either arm's validation region**
  (v1.1-v2.1 2.5-5.1 %, GRCh38-only-v1.1 3.3-6.5 %, GRCh38-only-v2.1 4.1-7.6 %); **Jaccard 0.014-0.034** (bins definition
  2.2-8.2 %, Jaccard 0.013-0.033). Between seeds of one graph: 3.3-10.1 %, Jaccard 0.017-0.051. Two independent 5 %
  samples would overlap by about 5 % of each. Same graph + seed across datasets: Jaccard 0.97-0.99 (same node IDs).
  The letter's "2-6 %" understates the upper end slightly (one-direction share 3.7-5.8 % from `SP/facts/ablation.md`
  §3.2; both directions 2.5-7.6 %).
- Effect on training data (COMPUTED NOW, tile): each arm trains on 93.8-95.5 % of the chr2-22 span; **4.4-6.3 % of one
  arm's training territory is validation (held out) in the other arm**; the territories trained on by exactly one of
  the two arms total **8.8-10.7 % of chr2-22**, so about 90 % of chr2-22 is training territory for both. On top of this
  the graphs produce different candidate tensors (training tensors 0.35-4.6 M per set), which is not a split effect.
- Consequence (already stated in R3.1): different validation regions -> different selected epochs and PASS thresholds;
  seed variation = initialisation + split.

---------------------------------------------------------------------------------------------------------------------

## Draft text

### D1. R1.1 item 1 (replacement)

1. *Haplotype-count channel (channel 6).* For every column it gives the number of haplotype paths through the node the
read traverses, relative to the graph's haplotype count. We retrained every HPRC SNV model without it, with alignments,
candidates, recipe, seed and number of GPUs unchanged, and compared each with its control trained the same way (Table
R1.1). Somatic average precision of the chr1 candidates, measured before population filtering, falls in all twelve pairs
(HG008T 0.20–0.35 to 0.07–0.09; COLO829T 0.67–0.83 to 0.48–0.59). SNV F1 after population filtering falls in ten of
twelve pairs (median −0.023). The two exceptions are HG008T ONT on v2.1 (+0.009, within the seed range of the model with
channel 6) and COLO829T Illumina on v1.1. There the control calls about 1,100 unlabelled candidates in two chr1 segments
(R3.1); the model without channel 6, selected at epoch 2, makes 297 false calls in all. The F1 loss is smaller than the AP loss because
population filtering removes many of the extra germline calls, so channel 6 and the filter partly overlap.

(Alternative if the 2026-10-03 COLO829T Illumina control is kept: "…falls in eleven of twelve pairs (median −0.026); the
exception is HG008T ONT on v2.1…", plus footnote † below.)

### D2. Table R1.1 (replacement; seed-0 pairs)

| chr1 SNVs, seed 0 (P / R / F1 after population filtering; somatic AP before filtering) | HPRC v1.1, with channel 6 | HPRC v1.1, without | HPRC v2.1, with channel 6 | HPRC v2.1, without |
|---|---|---|---|---|
| HG008T Illumina | 0.452 / 0.831 / 0.585 (0.199) | 0.407 / 0.809 / 0.542 (0.075) | 0.429 / 0.852 / 0.570 (0.310) | 0.398 / 0.783 / 0.528 (0.089) |
| HG008T PacBio\* | 0.497 / 0.828 / 0.621 (0.233) | 0.456 / 0.798 / 0.580 (0.078) | 0.473 / 0.839 / 0.605 (0.353) | 0.466 / 0.819 / 0.594 (0.075) |
| HG008T ONT\* | 0.466 / 0.826 / 0.596 (0.232) | 0.424 / 0.834 / 0.562 (0.065) | 0.432 / 0.859 / 0.575 (0.328) | 0.447 / 0.841 / 0.584 (0.071) |
| COLO829T Illumina\* | 0.554 / 0.713 / 0.624 (0.740)† | 0.850 / 0.687 / 0.760 (0.544) | 0.867 / 0.715 / 0.783 (0.818) | 0.362 / 0.678 / 0.472 (0.582) |
| COLO829T Fiber-seq\* | 0.857 / 0.717 / 0.781 (0.724) | 0.846 / 0.715 / 0.775 (0.575) | 0.866 / 0.707 / 0.779 (0.829) | 0.883 / 0.685 / 0.771 (0.592) |
| COLO829T ONT\* | 0.727 / 0.686 / 0.706 (0.667) | 0.713 / 0.675 / 0.693 (0.480) | 0.788 / 0.692 / 0.737 (0.755) | 0.703 / 0.649 / 0.675 (0.484) |

Caption: Each cell is one training run with seed 0; the two members of a pair share alignments, candidates, labels,
recipe, number of GPUs and checkpoint rule (the epoch with the highest validation somatic F1 in that run; PASS at its
validation recall-0.9 threshold; base quality, mapping quality and path count standardised with training-set
statistics). HPRC v2.1 pairs and HG008T Illumina v1.1: 4 GPUs (global batch 4,096); the models with channel 6 are the
first seeds of Table R3.1, whose three-seed means are given there. \*HPRC v1.1 columns: both members on 1 GPU (global
batch 1,024); Table R3.1's v1.1 models use 4 GPUs, and no 4-GPU model without channel 6 was trained for these five.
P / R / F1: PASS calls after population filtering in the benchmark regions. Somatic AP: average precision of the somatic
probability over all labelled chr1 candidates, before population filtering and independent of the threshold.
†The same seed and recipe trained earlier (with an earlier code version) selected epoch 6 instead of epoch 8
(validation F1 0.6871 against 0.6870) and scored 0.860 / 0.711 / 0.778 (AP 0.725). Seeds 1 and 2 of this recipe score F1
0.675 and 0.687 (AP 0.729).

### D3. R1.1 item 3 (one clause)

"Retrained without channel 6, both HPRC graphs are still above the GRCh38-only model of the same seed in all twelve
pairs (+0.004 to +0.292 SNV F1; +0.011 to +0.299 against the GRCh38-only three-seed mean, two of twelve within its seed
range)…"

### D4. Supplementary Table R3.9b: caption (numbers of the rows unchanged; they equal `pairs.tsv`)

First-seed models of Table R3.1, with the same checkpoints, validation-selected epochs and PASS thresholds as
Supplementary Table R3.9a. A − B in precision, recall and F1 (benchmark regions, after population filtering). Point
estimates are the differences on the whole of chr1, not bootstrap means. Intervals are the 2.5th and 97.5th percentiles
over 2,000 resamples of the 1-Mb chr1 blocks overlapping the benchmark regions; both members of a pair are scored on the
same resamples. For this analysis the GRCh38-only and HPRC v2.1 models, and HG008T Illumina HPRC v1.1, were re-run
on chr1 in full precision on a second cluster. These call sets differ from the GPU predictions of Table R3.9a by at most
3 true and 5 false positives per model (F1 within 0.002). As a result, a few point estimates differ from the R3.9a
differences in the third decimal. The largest case is HG008T PacBio: v2.1 recall is 0.835 here against 0.839, so
v2.1 − GRCh38-only Δrecall is +0.009 here against +0.013, and ΔF1 +0.090 against +0.092. The Pansoma − comparator rows use the
same first-seed HPRC v2.1 models (training-set standardisation). Comparators: DeepSomatic-TO 1.9.0 and ClairS-TO (0.4.2
on HG008T, 0.4.1 on COLO829T), tumor-only, scored with the same population filter.

(If exact agreement with R3.9a is preferred instead: sync the RIS `vcf_chr1/eval_bed/vcfeval/pon_calls` of the 13
seed-0 runs and rerun `tmp/assessment_20261030/main_bootstrap/pairs.py` with `BASES` = RIS first. The rows that would
move, by 0.001-0.004: HG008T PacBio all three comparisons, COLO829T Illumina v1.1 − GRCh38-only ΔP and v2.1 − v1.1
ΔR/ΔF1, COLO829T Fiber-seq ΔP of both HPRC − GRCh38-only rows, COLO829T ONT ΔP of all three and v2.1 − v1.1 ΔF1.)

### D5. chr2 confirmation: CI method sentence (R3.9 / R3.10 / R1.7 / Supplementary R1.7b)

"For each tumor–platform cell, we take the mean of the three seeds' F1 minus the comparator's F1. The overall difference is
the unweighted mean of this over the six cells. Its 95% interval is the 2.5th–97.5th percentile over 2,000 bootstrap
resamples of the 1-Mb chr2 blocks overlapping each tumor's benchmark regions. All seeds and both comparators of a tumor
are scored on the same resamples, and HG008T and COLO829T are resampled independently. The interval therefore covers
genomic sampling of chr2, not variation between training seeds."

### D6. R3.1 (optional, one parenthesis)

"…so they fall in different GRCh38 regions in each graph and for each seed (overlap 3–8% of each validation region,
as expected for independent 5% samples; about 90% of chr2–22 is training data for both arms)…"
