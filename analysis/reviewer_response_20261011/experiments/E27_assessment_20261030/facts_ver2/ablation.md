# Fact sheet: reference ablation, thresholds, split, ch6, context, seeds (assessment items R1.1, R1.2, R3.1, R3.9, R3.m1)

Compiled 2026-10-30. Status tags: **EXISTS** (number already in a repo/run file), **COMPUTED NOW** (computed today from existing
run outputs; scripts + outputs in `/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/`), **MISSING**.

Abbreviations: `RIS` = `/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/<run>/` (local copies of the RIS runs);
`CPU` = `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/` (tequila re-predictions of the RIS checkpoints);
`RISSUM` = `analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/{tables.md,summary.json,README.md}`.
All SNV numbers: chr1, benchmark BED, rtg vcfeval, PASS = the checkpoint's validation recall-0.9 threshold, "after PoN" = the
common population filter (4 resources, SNV AF >= 1e-4), unless stated.

---------------------------------------------------------------------------------------------------------------------------

## 1. Main per-dataset three-reference table (SNV)

### 1.0 Which setting, and is the recipe the same? (EXISTS)

- Two three-seed settings exist (RISSUM README "Results" and "The reported models with `--norm dataset`"):
  - **`--norm train`** (BQ/MAPQ/path-count z-scored with the training statistics): GRCh38-only (`<S>_<P>_SNV_linear_7fc7396{,_s1,_s2}`),
    HPRC v1.1 d9 matched 4-GPU (`<S>_<P>_SNV_7fc7396_ddp4_b1024{,_s1,_s2}`), HPRC v2.1 d46 (`<S>_<P>_SNV_d46_7fc7396{,_s1,_s2}`).
    **This is the only setting with all three arms x three seeds.** Seed-0 linear/d46 runs were trained before `--norm` existed
    (`args.json` `"norm": null` = training statistics, equivalent per RISSUM README bullet 1); d9 matched and all seed 1/2 runs have
    `"norm": "train"`.
  - **normbq** (`--norm dataset`, BQ per dataset, f344bca): d9 and d46 only (`<S>_<P>_SNV_{d9,d46}_normbq{,_s1,_s2}`). GRCh38-only
    normbq seeds were submitted but are **MISSING** (RISSUM `tables.md` rows "linear" = "-"; no `*_linear_normbq` run in RIS copy).
  - normds (BQ+MAPQ per dataset, 2026-10-16) is superseded (ONT predict bug, RISSUM README "--norm dataset: predict fits ...").
- Letter usage: R3.9a = normbq d9/d46 (COLO829T Illumina d9 at chr1-chosen epochs 6/5/6 + epoch-8 row); R3.9c and the
  "+0.028 to +0.292", "35 of 36" = `--norm train` (validation-selected everywhere); R1.2/R3.1 "−0.015 to +0.021" = normbq with
  COLO829T Illumina d9 at epochs 6/5/6; R1.2/R3.1 "+0.079 to +0.154" (before PoN) = a mix (five sets `--norm train`, COLO829T
  Illumina d9 normbq sweep epochs 6/5/6 vs d46 normds/normbq; R3.1.md evidence 2026-10-20 bullet 1).
- Recipe identity across arms (`args.json` of each run, checked): SNV = ConvNeXt 50.4 M (dims 96-768, depths 3/3/27/3), 20 epochs,
  lr 1e-4, wd 0.05, pow:0.75, batch 1,024/GPU, 4-GPU DDP, `val_fraction` 0.05, `val_block_nodes` 20,000, seeds 0/1/2, PASS at
  validation recall 0.9, checkpoint = best validation tensor-level somatic F1. **Same for all three arms in `--norm train`.**
  Not identical:
  (i) validation blocks are node-ID blocks, so they are different GRCh38 regions in each graph (section 3);
  (ii) numeric PASS thresholds differ (section 2);
  (iii) ch6 encodes path count relative to H = 17 / 90 / 464 (GRCh38-only ch6 is constant 100: E1.1 `results/tables.md` §3-4);
  (iv) Giraffe version: GRCh38-only 1.66.0, HG008T d9 1.66.0, COLO829T Illumina d9 <= 1.62, COLO829T long-read d9 1.64-1.66,
  d46 1.77.0 (points/R3.1.md evidence "vg 版本" and E22 check 2);
  (v) code commit: seed-0 linear/d46 ran earlier code than the d9 matched and seed-1/2 runs (RISSUM README "Code updates";
  torch.compile last-bit differences after f2a0e4e, README "Resolved 2026-10-19").

### 1.1 Main table, `--norm train`, validation-selected checkpoint for every model (COMPUTED NOW from EXISTS files)

P / R / F1 after PoN from RISSUM `summary.json` `seeds` (RIS originals); F1 before PoN = `report.json` `raw.at_threshold.f1`
(seed 0 of linear, d46 and HG008T Illumina d9 from the tequila re-predictions in CPU, the rest from RIS); best epoch from
`train.log` "done; best epoch N". Script `collect_main.py` -> `main_runs.json`, `main_table.py` -> `main_table.md`.
COLO829T Illumina v1.1 here = validation-selected epoch 8 on all three seeds.

| dataset | reference | P s0 / s1 / s2 (mean ± SD) | R | F1 after PoN | F1 before PoN | PASS threshold s0/s1/s2 | epoch |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | GRCh38-only | 0.334 / 0.394 / 0.388 (0.372 ± 0.033) | 0.786 / 0.783 / 0.789 (0.786 ± 0.003) | 0.469 / 0.524 / 0.520 (0.504 ± 0.030) | 0.009 / 0.010 / 0.010 (0.010 ± 0.000) | 0.064 / 0.129 / 0.092 | 6/3/3 |
| HG008T Illumina | HPRC v1.1 (d9) | 0.452 / 0.458 / 0.460 (0.457 ± 0.005) | 0.831 / 0.819 / 0.801 (0.817 ± 0.015) | 0.585 / 0.587 / 0.585 (0.586 ± 0.001) | 0.189 / 0.184 / 0.191 (0.188 ± 0.003) | 0.251 / 0.071 / 0.317 | 4/7/5 |
| HG008T Illumina | HPRC v2.1 (d46) | 0.429 / 0.457 / 0.432 (0.439 ± 0.016) | 0.852 / 0.828 / 0.815 (0.832 ± 0.019) | 0.570 / 0.589 / 0.565 (0.575 ± 0.013) | 0.295 / 0.336 / 0.322 (0.318 ± 0.021) | 0.324 / 0.554 / 0.340 | 5/5/6 |
| HG008T PacBio | GRCh38-only | 0.372 / 0.472 / 0.463 (0.436 ± 0.055) | 0.826 / 0.792 / 0.805 (0.808 ± 0.017) | 0.513 / 0.592 / 0.588 (0.564 ± 0.044) | 0.009 / 0.010 / 0.010 (0.010 ± 0.000) | 0.040 / 0.143 / 0.069 | 7/3/5 |
| HG008T PacBio | HPRC v1.1 (d9) | 0.479 / 0.443 / 0.464 (0.462 ± 0.018) | 0.829 / 0.808 / 0.831 (0.823 ± 0.013) | 0.608 / 0.572 / 0.595 (0.592 ± 0.018) | 0.181 / 0.178 / 0.170 (0.176 ± 0.006) | 0.214 / 0.117 / 0.114 | 8/11/6 |
| HG008T PacBio | HPRC v2.1 (d46) | 0.473 / 0.477 / 0.487 (0.479 ± 0.007) | 0.839 / 0.822 / 0.821 (0.827 ± 0.010) | 0.605 / 0.604 / 0.611 (0.607 ± 0.004) | 0.312 / 0.343 / 0.337 (0.331 ± 0.017) | 0.291 / 0.430 / 0.325 | 7/9/9 |
| HG008T ONT | GRCh38-only | 0.394 / 0.368 / 0.394 (0.385 ± 0.015) | 0.842 / 0.826 / 0.780 (0.816 ± 0.032) | 0.537 / 0.509 / 0.524 (0.523 ± 0.014) | 0.008 / 0.008 / 0.009 (0.009 ± 0.001) | 0.042 / 0.050 / 0.145 | 8/8/6 |
| HG008T ONT | HPRC v1.1 (d9) | 0.473 / 0.451 / 0.447 (0.457 ± 0.014) | 0.821 / 0.816 / 0.822 (0.820 ± 0.003) | 0.600 / 0.581 / 0.579 (0.587 ± 0.012) | 0.175 / 0.176 / 0.165 (0.172 ± 0.006) | 0.328 / 0.116 / 0.245 | 9/9/8 |
| HG008T ONT | HPRC v2.1 (d46) | 0.432 / 0.467 / 0.474 (0.457 ± 0.022) | 0.859 / 0.824 / 0.809 (0.831 ± 0.026) | 0.575 / 0.596 / 0.598 (0.589 ± 0.013) | 0.287 / 0.331 / 0.330 (0.316 ± 0.025) | 0.231 / 0.482 / 0.592 | 9/8/6 |
| COLO829T Illumina | GRCh38-only | 0.391 / 0.289 / 0.499 (0.393 ± 0.105) | 0.583 / 0.583 / 0.564 (0.576 ± 0.011) | 0.468 / 0.386 / 0.530 (0.461 ± 0.072) | 0.073 / 0.076 / 0.094 (0.081 ± 0.011) | 0.081 / 0.127 / 0.161 | 7/6/6 |
| COLO829T Illumina | HPRC v1.1 (d9), epoch 8 (validation) | 0.484 / 0.562 / 0.475 (0.507 ± 0.048) | 0.698 / 0.705 / 0.702 (0.702 ± 0.004) | 0.572 / 0.626 / 0.566 (0.588 ± 0.033) | 0.406 / 0.456 / 0.403 (0.422 ± 0.030) | 0.327 / 0.435 / 0.373 | 8/8/8 |
| COLO829T Illumina | HPRC v2.1 (d46) | 0.867 / 0.704 / 0.820 (0.797 ± 0.084) | 0.715 / 0.705 / 0.717 (0.712 ± 0.007) | 0.783 / 0.704 / 0.765 (0.751 ± 0.042) | 0.684 / 0.601 / 0.647 (0.644 ± 0.042) | 0.640 / 0.703 / 0.303 | 7/6/7 |
| COLO829T Fiber-seq | GRCh38-only | 0.745 / 0.841 / 0.868 (0.818 ± 0.065) | 0.625 / 0.600 / 0.576 (0.600 ± 0.024) | 0.680 / 0.700 / 0.692 (0.691 ± 0.010) | 0.075 / 0.101 / 0.112 (0.096 ± 0.019) | 0.085 / 0.166 / 0.129 | 10/9/7 |
| COLO829T Fiber-seq | HPRC v1.1 (d9) | 0.785 / 0.876 / 0.845 (0.835 ± 0.046) | 0.723 / 0.707 / 0.709 (0.713 ± 0.009) | 0.752 / 0.782 / 0.771 (0.768 ± 0.015) | 0.523 / 0.562 / 0.545 (0.544 ± 0.019) | 0.364 / 0.306 / 0.431 | 9/9/9 |
| COLO829T Fiber-seq | HPRC v2.1 (d46) | 0.866 / 0.846 / 0.849 (0.854 ± 0.011) | 0.707 / 0.713 / 0.737 (0.719 ± 0.016) | 0.779 / 0.774 / 0.789 (0.780 ± 0.008) | 0.683 / 0.657 / 0.664 (0.668 ± 0.013) | 0.481 / 0.437 / 0.344 | 9/8/10 |
| COLO829T ONT | GRCh38-only | 0.340 / 0.353 / 0.475 (0.390 ± 0.075) | 0.543 / 0.540 / 0.569 (0.551 ± 0.016) | 0.418 / 0.427 / 0.518 (0.454 ± 0.055) | 0.066 / 0.066 / 0.080 (0.071 ± 0.008) | 0.092 / 0.112 / 0.152 | 6/7/7 |
| COLO829T ONT | HPRC v1.1 (d9) | 0.793 / 0.816 / 0.833 (0.814 ± 0.020) | 0.701 / 0.693 / 0.673 (0.689 ± 0.014) | 0.744 / 0.749 / 0.745 (0.746 ± 0.003) | 0.503 / 0.533 / 0.537 (0.524 ± 0.019) | 0.335 / 0.407 / 0.425 | 9/9/6 |
| COLO829T ONT | HPRC v2.1 (d46) | 0.788 / 0.734 / 0.744 (0.755 ± 0.029) | 0.692 / 0.703 / 0.718 (0.704 ± 0.013) | 0.737 / 0.718 / 0.730 (0.729 ± 0.009) | 0.626 / 0.588 / 0.595 (0.603 ± 0.020) | 0.568 / 0.604 / 0.535 | 10/10/8 |

(P/R/F1 after PoN of this table = R3.9c of the letter exactly. Seed-0 re-predictions in CPU reproduce RIS within 0.001 F1,
e.g. HG008T Illumina d9 0.5854 vs 0.5851.)

### 1.2 Derived numbers, validation-selected (`--norm train`) (COMPUTED NOW)

HPRC − GRCh38-only, three-seed means:

| dataset | arm | ΔF1 after PoN | ΔP | ΔR | ΔF1 before PoN | seeds HPRC > GRCh38-only | per-seed ΔF1 |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | v1.1 | +0.081 | +0.085 | +0.031 | +0.178 | 3/3 | +0.116 / +0.063 / +0.065 |
| HG008T Illumina | v2.1 | +0.070 | +0.067 | +0.045 | +0.308 | 3/3 | +0.101 / +0.065 / +0.045 |
| HG008T PacBio | v1.1 | +0.028 | +0.027 | +0.015 | +0.167 | 2/3 | +0.095 / −0.019 / +0.008 |
| HG008T PacBio | v2.1 | +0.042 | +0.043 | +0.020 | +0.321 | 3/3 | +0.092 / +0.012 / +0.023 |
| HG008T ONT | v1.1 | +0.064 | +0.072 | +0.003 | +0.163 | 3/3 | +0.064 / +0.072 / +0.056 |
| HG008T ONT | v2.1 | +0.066 | +0.072 | +0.014 | +0.308 | 3/3 | +0.038 / +0.087 / +0.074 |
| COLO829T Illumina | v1.1 (epoch 8) | +0.127 | +0.114 | +0.125 | +0.341 | 3/3 | +0.104 / +0.240 / +0.037 |
| COLO829T Illumina | v2.1 | +0.290 | +0.404 | +0.136 | +0.563 | 3/3 | +0.316 / +0.318 / +0.235 |
| COLO829T Fiber-seq | v1.1 | +0.078 | +0.017 | +0.112 | +0.447 | 3/3 | +0.073 / +0.082 / +0.079 |
| COLO829T Fiber-seq | v2.1 | +0.090 | +0.036 | +0.119 | +0.572 | 3/3 | +0.099 / +0.073 / +0.097 |
| COLO829T ONT | v1.1 | +0.292 | +0.425 | +0.138 | +0.454 | 3/3 | +0.326 / +0.322 / +0.227 |
| COLO829T ONT | v2.1 | +0.274 | +0.366 | +0.154 | +0.532 | 3/3 | +0.319 / +0.291 / +0.213 |

- **HPRC − GRCh38-only ΔF1 after PoN: +0.028 to +0.292** (unchanged by the checkpoint question: this setting already uses epoch 8
  for COLO829T Illumina v1.1; the RIS README "Graph - linear, three-seed means: PASS +0.028 to +0.292" is the same). ΔP +0.017 to
  +0.425, ΔR +0.003 to +0.154 (both positive on all 12). Best F1 (test-tuned): +0.015 to +0.211 (RISSUM README).
- **Seed-matched: 35 / 36** HPRC > GRCh38-only (the exception: HG008T PacBio v1.1 seed 1, 0.572 vs 0.592).
- **v2.1 − v1.1, validation-selected (`--norm train`)**: after PoN −0.018 (COLO829T ONT; RISSUM rounds −0.017) to +0.163 (COLO829T
  Illumina); without COLO829T Illumina −0.018 to +0.015. Before PoN **+0.079 to +0.222** (precision +0.089 to +0.280; recall −0.000 to
  +0.016). Best F1 (test-tuned) d46 − d9: −0.014 to +0.091 (RISSUM README "d46 - d9, three-seed means").
- Per set, d46 vs d9 seeds (`--norm train`): d46 higher on 1/3, 2/3, 2/3, 3/3, 2/3, 0/3 seeds; seed ranges do **not** overlap for
  COLO829T Illumina (d46 above) and COLO829T ONT (d46 0.718-0.737 below d9 0.744-0.749).

### 1.3 Secondary: normbq (BQ per dataset), d9 / d46 only (COMPUTED NOW from EXISTS)

| dataset | reference | P (mean ± SD) | R | F1 after PoN | F1 before PoN | PASS thr. | epoch |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | v1.1 | 0.456/0.459/0.461 (0.459 ± 0.003) | 0.831/0.819/0.803 (0.818 ± 0.014) | 0.588/0.588/0.586 (0.588 ± 0.001) | 0.189/0.185/0.192 (0.189 ± 0.003) | 0.256/0.070/0.317 | 4/7/5 |
| HG008T Illumina | v2.1 | 0.426/0.451/0.433 (0.437 ± 0.013) | 0.852/0.831/0.814 (0.832 ± 0.019) | 0.568/0.585/0.566 (0.573 ± 0.010) | 0.294/0.332/0.323 (0.316 ± 0.020) | 0.317/0.550/0.352 | 5/5/6 |
| HG008T PacBio | v1.1 | 0.484/0.422/0.460 (0.455 ± 0.031) | 0.828/0.811/0.829 (0.823 ± 0.010) | 0.611/0.555/0.592 (0.586 ± 0.028) | 0.182/0.175/0.169 (0.175 ± 0.007) | 0.231/0.101/0.107 | 8/11/6 |
| HG008T PacBio | v2.1 | 0.471/0.487/0.485 (0.481 ± 0.009) | 0.838/0.814/0.814 (0.822 ± 0.014) | 0.603/0.609/0.607 (0.606 ± 0.003) | 0.307/0.351/0.339 (0.332 ± 0.023) | 0.285/0.400/0.345 | 7/7/7 |
| HG008T ONT | v1.1 | 0.472/0.443/0.447 (0.454 ± 0.016) | 0.814/0.819/0.831 (0.821 ± 0.009) | 0.598/0.575/0.581 (0.585 ± 0.012) | 0.175/0.172/0.163 (0.170 ± 0.006) | 0.322/0.091/0.091 | 9/9/7 |
| HG008T ONT | v2.1 | 0.433/0.470/0.481 (0.461 ± 0.025) | 0.859/0.848/0.824 (0.844 ± 0.018) | 0.576/0.605/0.607 (0.596 ± 0.017) | 0.287/0.315/0.335 (0.312 ± 0.024) | 0.221/0.409/0.311 | 9/12/8 |
| COLO829T Illumina | v1.1 epoch 8 (validation) | 0.476/0.565/0.477 (0.506 ± 0.051) | 0.698/0.706/0.703 (0.702 ± 0.004) | 0.566/0.628/0.568 (0.587 ± 0.035) | 0.400/0.458/0.404 (0.421 ± 0.032) | 0.320/0.441/0.373 | 8/8/8 |
| COLO829T Illumina | v1.1 epochs 6/5/6 (chosen on chr1) | 0.832/0.705/0.782 (0.773 ± 0.064) | 0.723/0.705/0.698 (0.708 ± 0.013) | 0.773/0.705/0.738 (0.739 ± 0.034) | 0.531/0.518/0.546 (0.532 ± 0.014) | 0.277/0.492/0.336 | 6/5/6 |
| COLO829T Illumina | v2.1 | 0.867/0.706/0.821 (0.798 ± 0.083) | 0.712/0.705/0.717 (0.712 ± 0.006) | 0.782/0.705/0.765 (0.751 ± 0.040) | 0.684/0.602/0.647 (0.644 ± 0.041) | 0.626/0.701/0.304 | 7/6/7 |
| COLO829T Fiber-seq | v1.1 | 0.798/0.876/0.843 (0.839 ± 0.039) | 0.722/0.708/0.710 (0.713 ± 0.008) | 0.758/0.783/0.771 (0.771 ± 0.012) | 0.528/0.564/0.544 (0.545 ± 0.018) | 0.368/0.309/0.405 | 9/9/9 |
| COLO829T Fiber-seq | v2.1 | 0.863/0.848/0.852 (0.854 ± 0.008) | 0.707/0.710/0.734 (0.717 ± 0.015) | 0.777/0.773/0.788 (0.779 ± 0.008) | 0.680/0.658/0.665 (0.668 ± 0.011) | 0.479/0.470/0.368 | 9/8/10 |
| COLO829T ONT | v1.1 | 0.793/0.822/0.830 (0.815 ± 0.020) | 0.700/0.693/0.675 (0.689 ± 0.013) | 0.744/0.752/0.745 (0.747 ± 0.004) | 0.505/0.534/0.532 (0.524 ± 0.016) | 0.339/0.385/0.445 | 9/9/6 |
| COLO829T ONT | v2.1 | 0.814/0.753/0.740 (0.769 ± 0.039) | 0.700/0.701/0.718 (0.706 ± 0.010) | 0.753/0.726/0.729 (0.736 ± 0.015) | 0.642/0.601/0.593 (0.612 ± 0.027) | 0.521/0.636/0.520 | 9/10/8 |

v2.1 − v1.1 (normbq): after PoN **−0.015 to +0.164** with epoch 8 (−0.015 to +0.021 only with COLO829T Illumina at epochs 6/5/6);
before PoN **+0.088 to +0.224** with epoch 8 (+0.088 to +0.157 with 6/5/6; COLO829T Illumina alone +0.224 vs +0.112), precision
before PoN +0.088 to +0.282. Best F1 d46 − d9: COLO829T Illumina +0.092 at epoch 8 vs +0.012 at 6/5/6 (RISSUM tables.md
"Arm differences"). Seeds: HG008T Illumina d46 below d9 on all 3 seeds (0.568/0.585/0.566 vs 0.588/0.588/0.586; ranges do not
overlap) even in this setting.

Seed-SD ranges over the 12 HPRC rows (for R3.9 point 1): normbq with epoch 8: P SD 0.003-0.083, **R SD 0.004-0.019**, F1 SD
0.001-0.040 (median 0.012); normbq with 6/5/6: R SD 0.006-0.019 (the letter's numbers); `--norm train`: P 0.005-0.084, R 0.003-0.026,
F1 0.001-0.042 (median 0.012); GRCh38-only: P 0.015-0.105, R 0.003-0.032, F1 0.010-0.072.

### 1.4 INDEL, three references (EXISTS: `analysis/reviewer_response_20261011/experiments/E3b_rescue_off_and_indel_threshold/tables.md` §2)

One seed per arm; recipe = the 7fc7396 INDEL recipe for all three arms (8.7 M, sqrt weights, lr 1e-4, 12 epochs, 1 GPU, no site
features, best epoch by tensor-level validation F1; `args.json` of `*_INDEL_{linear_,,d46_}7fc7396`); COLO829T v1.1 = the matched
lr-1e-4 rerun `COLO829_*_INDEL_7fc7396`. PoN = the 71cb51a INDEL rule. **This is not the final INDEL recipe** of R3.10/R1.6
(site features, `--best-by postpon`, `--norm dataset`, non-reasons, PASS ppr8 = post-PoN validation recall 0.8: `args.json` of
`unseen_s0_*_INDEL`, `tmp/indel_gap_20261020/final_20261023/PREREGISTRATION.md` lines 20, 79-81); no three-reference run with the
final recipe exists (MISSING).

| dataset | rule | GRCh38-only | v1.1 | v2.1 |
|---|---|---|---|---|
| HG008T Illumina | r0.9 PASS F1 / valbest F1 / test-tuned best | 0.146 / 0.080 / 0.202 | 0.251 / 0.218 / 0.270 | 0.253 / 0.207 / 0.256 |
| HG008T PacBio | same | 0.133 / 0.036 / 0.176 | 0.092 / 0.227 / 0.262 | 0.098 / 0.200 / 0.275 |
| HG008T ONT | same | 0.023 / 0.048 / 0.085 | 0.015 / 0.117 / 0.144 | 0.015 / 0.154 / 0.159 |
| COLO829T Illumina | same | 0.034 / 0.030 (2 TP/1 FP) / 0.070 | 0.104 / 0.080 (6/14) / 0.161 | 0.117 / 0.044 (3/2) / 0.168 |
| COLO829T Fiber-seq | same | 0.044 / no call / 0.095 | 0.024 / 0.041 (3/11) / 0.142 | 0.021 / 0.093 (8/33) / 0.137 |
| COLO829T ONT | same | 0.004 / 0.000 (0/6) / 0.020 | 0.003 / 0.061 (5/27) / 0.095 | 0.002 / 0.027 (2/14) / 0.089 |

- valbest (validation best truth-level F1 threshold): HPRC higher in 12/12 (COLO829T: 3-41 PASS calls per HPRC run).
- r0.9: GRCh38-only higher on the four long-read datasets (8 of 12 pairs lower for HPRC: HG008T PacBio/ONT, COLO829T Fiber-seq/ONT).
- d46 − d9 INDEL: r0.9 PASS −0.003 to +0.013 (max |Δ| 0.013, COLO829T Illumina 0.104 → 0.117); test-tuned best −0.013 to +0.015
  (RISSUM tables.md "Matched d9"); valbest −0.036 to +0.052 (COLO829T Illumina 0.080 → 0.044, Fiber-seq 0.041 → 0.093).

## 2. Numeric PASS thresholds (EXISTS: each run's `vcf_chr1/eval_bed/report.json` "threshold" = `*.graph.stats.json` "threshold")

Rule (code): `graph_vcf.py` lines 36-38 / 178-182: PASS = `p_somatic >=` the checkpoint's `val.somatic_at_recall["0.9"].threshold`,
i.e. the highest p_somatic at which the per-tensor validation somatic recall reaches 0.9, at the best epoch (verified:
`COLO829_Illumina_SNV_7fc7396_ddp4_b1024/metrics.jsonl` epoch 8 contains 0.3267, = its report threshold). Ranges over seeds 0-2:

| dataset | GRCh38-only (`--norm train`) | v1.1 (`--norm train`) | v2.1 (`--norm train`) | v1.1 (normbq) | v2.1 (normbq) |
|---|---|---|---|---|---|
| HG008T Illumina | 0.064-0.129 | 0.071-0.317 | 0.324-0.554 | 0.070-0.317 | 0.317-0.550 |
| HG008T PacBio | 0.040-0.143 | 0.114-0.214 | 0.291-0.430 | 0.101-0.231 | 0.285-0.400 |
| HG008T ONT | 0.042-0.145 | 0.116-0.328 | 0.231-0.592 | 0.091-0.322 | 0.221-0.409 |
| COLO829T Illumina | 0.081-0.161 | 0.327-0.435 (ep. 8) | 0.303-0.703 | 0.320-0.441 (ep. 8); 0.277-0.492 (ep. 6/5/6) | 0.304-0.701 |
| COLO829T Fiber-seq | 0.085-0.166 | 0.306-0.431 | 0.344-0.481 | 0.309-0.405 | 0.368-0.479 |
| COLO829T ONT | 0.092-0.152 | 0.335-0.425 | 0.535-0.604 | 0.339-0.445 | 0.520-0.636 |

Per seed: section 1.1 / 1.3 tables. GRCh38-only thresholds are 2-10x lower (its validation somatic AP is ~0.01, RISSUM README).
INDEL (one seed, E3b §2 "threshold r0.9 / valbest"): GRCh38-only 0.005-0.048 / 0.116-0.357; v1.1 0.0046-0.120 / 0.215-0.466;
v2.1 0.012-0.107 / 0.277-0.400.

## 3. Validation split

### 3.1 Code (EXISTS; `machine_learning/pansoma_net_v2/data.py` 307-322, `train.py` 549-555)

- `block_split`: tensors of the training chromosomes (chr2-22; train.log "train [...]", test chr1) are grouped into blocks
  `node_id // 20000` per chromosome; a block is validation if `_block_u(seed, chrom, block) < 0.05`, where `_block_u` =
  SHA-1 of the string `"{seed}:{chrom}:{block}"` mapped to [0, 1). **Deterministic, seeded hash; the training seed (`args.seed`)
  is the hash seed, so seeds 0 / 1 / 2 give different validation blocks.** Block membership uses all tensors (labels −1 included,
  `labelled=False`); training keeps label >= 0, validation keeps `eval_label` >= 0.
- Node IDs differ between graphs, so block k is a different GRCh38 interval in each graph (and the GRCh38-only graph has 32-bp
  nodes); same seed does not mean same genomic split.
- Logged labelled sizes (train.log "training tensors ..., validation ..."; validation somatic tensors s0/s1/s2): HG008T Illumina
  linear 484/487/366, d9 469/343/424, d46 357/340/346; COLO829T Illumina linear 1,775/2,039/1,654, d9 2,060/1,579/2,057, d46
  1,737/1,724/1,639 (all 18 arms x 3 seeds in `tmp/assessment_20261030/ablation/q3/val_counts.txt`).

### 3.2 GRCh38 overlap of the validation regions (COMPUTED NOW; Slurm job 462875 on meow,
`tmp/assessment_20261030/ablation/q3/{val_split.py,run.sh,val_split.json,q3_462875.out}`)

Input: the SNV index caches `/scratch/jshen/data/pansoma_net_v2_runs/index_cache/<set>.tensors.SNV.<hash>.npz` (per tensor chrom,
node, pos0, off_reference), the 18 sets (6 datasets x 3 graphs), seeds 0-2, chr2-22, re-implementing `block_split` exactly.
Two definitions of a block's GRCh38 extent: **tile** = the code's own `validation_truth` rule (a block starts at its smallest pos0;
a position belongs to the last block starting at or before it); **bins** = 10-kb GRCh38 bins whose on-reference tensors are > 50 %
in validation blocks (robust to scattered nodes; only bins with tensors count). (A 'hull' [min, max] definition is unusable on the
HPRC graphs: blocks contain far-scattered tensors, hull totals 274-777 Mb.)

| graph | seed | validation blocks / all blocks | validation GRCh38 bp (tile) | (bins) |
|---|---|---|---|---|
| GRCh38-only | 0 / 1 / 2 | 220-221 / 229 / 183 of 3,966-3,975 | 140.8-141.8 / 153.4-153.5 / 116.8-117.2 Mb | 124-136 / 130-143 / 104-114 Mb |
| v1.1 d9 | 0 / 1 / 2 | 150 / 130 / 153-154 of 2,664-2,669 | 142.7-144.5 / 126.0-128.2 / 155.6-157.8 Mb | 73-138 / 65-116 / 82-145 Mb |
| v2.1 d46 | 0 / 1 / 2 | 291-296 / 283-289 / 284-289 of 6,114-6,238 | 115.7-116.7 / 115.2-116.6 / 121.3-123.3 Mb | 59-104 / 57-104 / 61-107 Mb |

(ranges over the six datasets; chr2-22 span covered by tensors ≈ 2,559-2,561 Mb, so validation ≈ 4.5-6.2 % of it.)

Overlap (ranges over the six datasets):

| comparison | intersection (tile) | Jaccard (tile) | share of A's validation bp inside B's (tile) | Jaccard (bins) |
|---|---|---|---|---|
| GRCh38-only vs v1.1, same seed (s0 / s1 / s2) | 5.3-5.6 / 7.5-8.3 / 5.2-5.5 Mb | 0.019-0.020 / 0.027-0.030 / 0.019-0.020 | 3.7-3.9 % / 4.9-5.4 % / 4.5-4.7 % | 0.014-0.028 |
| GRCh38-only vs v2.1, same seed | 5.8-6.0 / 8.6-8.8 / 6.2-6.6 Mb | 0.023-0.024 / 0.033-0.034 / 0.026-0.028 | 4.1-5.8 % | 0.016-0.033 |
| v1.1 vs v2.1, same seed | 4.8-5.0 / 5.8-5.9 / 3.9-4.0 Mb | 0.019-0.020 / 0.025 / 0.014-0.015 | 2.5-4.7 % | 0.013-0.023 |
| same graph, seed pairs (0-1, 0-2, 1-2) | GRCh38-only 6.4-10.2, v1.1 7.7-14.6, v2.1 4.0-6.3 Mb | 0.017-0.051 | 3.5-10.1 % | 0.018-0.048 |
| same graph + seed, different dataset (HG008T Illumina vs COLO829T ONT, s0) | 114.7-142.0 Mb | 0.972-0.991 | 98-99.7 % | 0.71-0.94 |

**Conclusion:** the validation regions of the three graphs are, in GRCh38 terms, essentially independent 5 % samples (overlap
≈ 2-6 % of each, Jaccard 0.014-0.034, what two independent 5 % draws give), and the three seeds of one graph are likewise nearly
disjoint. Same graph + same seed gives the same regions across datasets (as expected: same node IDs). So the reference comparison
pairs models whose validation sets (and so their PASS thresholds and selected epochs) come from different genomic regions; seed
variation is joint initialisation + split variation, not initialisation alone. chr1 test and chr2-22 training are common.

## 4. Paired block bootstrap (E16 / supplementary R3.9)

### 4.1 What the existing tables used (EXISTS)

`experiments/E16_bootstrap/paired.tsv` (rows `arms`, BED, 1 Mb; 2,000 resamples of 1-Mb chr1 blocks, seed 20261015) and
`supplementary/R3.9/make_paired_table.py` (line 5: "first-seed in-sample models (COLO829T Illumina v1.1: epoch 6, E16's
ARM_RUN)") -> `paired_table.tsv` = letter Table R3.9b. **One seed (seed 0) for every pair; no aggregation over seeds.** Runs:
- GRCh38-only and v2.1: `<S>_<P>_SNV_{linear,d46}_7fc7396` seed 0 (`--norm train`; tequila re-predictions in CPU, F1 within
  0.002 of RIS) = main-table seed 0.
- v1.1: HG008T Illumina `HG008_Illumina_SNV_7fc7396_ddp4_b1024` (matched, = main table); **HG008T PacBio / ONT the 1-GPU runs
  `HG008_{PacBio,ONT}_SNV_7fc7396` (F1 0.621 / 0.596), not the matched 4-GPU seed 0 of the main table (0.608 / 0.600)**; COLO829T
  Fiber-seq / ONT the matched `*_7fc7396_ddp4_b1024` (= main table); **COLO829T Illumina `COLO829_Illumina_SNV_d9_normbq_sweep_e06`
  (normbq, chr1-chosen epoch 6, F1 0.773), not the main table's 0.572 (epoch 8)**.
- So the R3.9b point estimates match the main `--norm train` table seed 0 except three v1.1 pairs (HG008T PacBio, HG008T ONT,
  COLO829T Illumina). The earlier E16 version (R3.1.md evidence 2026-10-19) had COLO829T Illumina v1.1 at epoch 8: +0.104
  [+0.022, +0.151].

### 4.2 Recomputed with validation-selected matched checkpoints, all three seeds (COMPUTED NOW)

E16's own functions (`e16_bootstrap.BlockSets/positions/block_counts/metrics/ci`, same resamples) on the `--norm train` main-table
runs, seeds 0/1/2, every pair; Slurm job 462896 (meow, 14 s); `tmp/assessment_20261030/ablation/q4/{e16_pairs_q4.py,run.sh,
e16_pairs_q4.tsv,e16_pairs_q4.md,summary_q4.md}`. rtg TP/FN/FP totals equal each run's report.json for all 54 call sets.
ΔF1 [95 % CI], 1-Mb blocks, BED, after PoN:

| dataset | pair | seed 0 | seed 1 | seed 2 |
|---|---|---|---|---|
| HG008T Illumina | v1.1 − GRCh38-only | +0.116 [+0.062, +0.171] | +0.063 [+0.040, +0.088] | +0.065 [+0.031, +0.102] |
| HG008T Illumina | v2.1 − GRCh38-only | +0.101 [+0.050, +0.152] | +0.065 [+0.036, +0.093] | +0.045 [+0.013, +0.078] |
| HG008T Illumina | v2.1 − v1.1 | −0.015 [−0.030, −0.001] | +0.002 [−0.029, +0.026] | −0.020 [−0.038, −0.001] |
| HG008T PacBio | v1.1 − GRCh38-only | +0.094 [+0.031, +0.174] | −0.019 [−0.088, +0.029] | +0.008 [−0.009, +0.025] |
| HG008T PacBio | v2.1 − GRCh38-only | +0.090 [+0.028, +0.172] | +0.012 [−0.007, +0.029] | +0.023 [+0.006, +0.042] |
| HG008T PacBio | v2.1 − v1.1 | −0.004 [−0.020, +0.014] | +0.031 [−0.014, +0.096] | +0.016 [+0.002, +0.030] |
| HG008T ONT | v1.1 − GRCh38-only | +0.064 [+0.037, +0.095] | +0.072 [+0.032, +0.119] | +0.056 [+0.031, +0.082] |
| HG008T ONT | v2.1 − GRCh38-only | +0.038 [+0.010, +0.069] | +0.087 [+0.046, +0.133] | +0.074 [+0.052, +0.096] |
| HG008T ONT | v2.1 − v1.1 | −0.025 [−0.041, −0.011] | +0.015 [+0.001, +0.027] | +0.018 [+0.001, +0.036] |
| COLO829T Illumina | v1.1 (epoch 8) − GRCh38-only | +0.104 [+0.022, +0.151] | +0.240 [+0.103, +0.274] | +0.037 [+0.008, +0.077] |
| COLO829T Illumina | v2.1 − GRCh38-only | +0.316 [+0.116, +0.453] | +0.318 [+0.118, +0.390] | +0.235 [+0.119, +0.334] |
| COLO829T Illumina | v2.1 − v1.1 (epoch 8) | +0.211 [+0.061, +0.340] | +0.078 [+0.014, +0.131] | +0.198 [+0.048, +0.318] |
| COLO829T Fiber-seq | v1.1 − GRCh38-only | +0.073 [+0.002, +0.132] | +0.082 [+0.056, +0.108] | +0.078 [+0.050, +0.104] |
| COLO829T Fiber-seq | v2.1 − GRCh38-only | +0.099 [+0.062, +0.141] | +0.073 [+0.047, +0.099] | +0.097 [+0.076, +0.117] |
| COLO829T Fiber-seq | v2.1 − v1.1 | +0.026 [−0.014, +0.089] | −0.009 [−0.022, +0.003] | +0.018 [+0.003, +0.040] |
| COLO829T ONT | v1.1 − GRCh38-only | +0.326 [+0.219, +0.422] | +0.322 [+0.154, +0.459] | +0.227 [+0.089, +0.367] |
| COLO829T ONT | v2.1 − GRCh38-only | +0.319 [+0.208, +0.417] | +0.291 [+0.132, +0.418] | +0.213 [+0.072, +0.352] |
| COLO829T ONT | v2.1 − v1.1 | −0.008 [−0.035, +0.010] | −0.031 [−0.064, −0.004] | −0.014 [−0.057, +0.014] |

- **HPRC − GRCh38-only, seed 0 (validation-selected, matched): ΔF1 +0.038 to +0.326; CI excludes zero in 12/12 for F1, 10/12 for
  precision, 8/12 for recall.** Seed 1: −0.019 to +0.322, 10/12 (P 8, R 9); seed 2: +0.008 to +0.235, 11/12 (P 8, R 7). Over the
  36 seed-matched pairs: 33 exclude zero, none significantly below zero.
- **v2.1 − v1.1 with validation-selected v1.1**: seed 0 −0.025 to +0.211 (CI above 0: 1 [COLO829T Illumina], below 0: 2 [HG008T
  Illumina, ONT]); seed 1 −0.031 to +0.078 (2 above, 1 below); seed 2 −0.020 to +0.198 (4 above, 1 below). COLO829T Illumina v2.1 −
  v1.1 excludes zero on all three seeds (+0.211, +0.078, +0.198).
- COLO829T Illumina at the chr1-chosen epochs (normbq; `tmp/assessment_20261030/ablation/q9/e16_pairs_q9.md`): v1.1 − GRCh38-only
  +0.306 [+0.107, +0.441] / +0.319 [+0.109, +0.402] / +0.208 [+0.094, +0.301] (seeds 0/1/2 at epochs 6/5/6); v2.1 − v1.1 +0.010
  [+0.000, +0.019] / −0.001 [−0.021, +0.029] / +0.027 [+0.008, +0.055]. normbq epoch 8: v1.1 − GRCh38-only +0.098 [+0.016, +0.145] /
  +0.242 [+0.104, +0.277] / +0.039 [+0.009, +0.079]; v2.1 − v1.1 +0.217 [+0.064, +0.348] / +0.076 [+0.013, +0.130] / +0.197 [+0.048,
  +0.315].
- With 5-Mb blocks several CIs widen to include zero (e.g. COLO829T Illumina v1.1 epoch 8 − GRCh38-only seed 0 +0.104 [−0.018,
  +0.148]); see `e16_pairs_q4.md`.

## 5. Channel 6 (path count) ablation (COMPUTED NOW from EXISTS; `tmp/assessment_20261030/ablation/q5_ch6.py` -> `q5_ch6.md`)

Same alignments and candidate tensors; no-ch6 = `train --drop-planes path_count` (`args.json` `drop_planes`), one seed. P / R / F1
after PoN (BED, PASS r0.9) from each run's `vcf_chr1/eval_bed/report.json` (COLO829T d9 controls from RISSUM `summary.json`
`no_ch6.d9_with`); somatic AP = chr1 per-tensor, before PoN, `test_chr1/*.metrics.json` `somatic_ap`.

| dataset | graph (control recipe) | with ch6: P / R / F1 (AP) | no ch6: P / R / F1 (AP) | ΔF1 |
|---|---|---|---|---|
| HG008T Illumina | v1.1 (4 GPU) | 0.452 / 0.831 / 0.585 (0.199) | 0.407 / 0.809 / 0.542 (0.075) | −0.043 |
| HG008T PacBio | v1.1 (1 GPU, `HG008_PacBio_SNV_7fc7396`) | 0.497 / 0.828 / 0.621 (0.233) | 0.456 / 0.798 / 0.580 (0.078) | −0.041 |
| HG008T ONT | v1.1 (1 GPU) | 0.466 / 0.826 / 0.596 (0.232) | 0.424 / 0.834 / 0.562 (0.065) | −0.034 |
| COLO829T Illumina | v1.1 (1 GPU `_h90`, its own validation epoch 6) | 0.860 / 0.711 / 0.778 (0.725, RISSUM README l.61) | 0.850 / 0.687 / 0.760 (0.544) | −0.019 |
| COLO829T Fiber-seq | v1.1 (1 GPU `_h90`) | 0.857 / 0.717 / 0.781 (AP MISSING locally) | 0.846 / 0.715 / 0.775 (0.575) | −0.006 |
| COLO829T ONT | v1.1 (1 GPU `_h90`) | 0.727 / 0.686 / 0.706 (AP MISSING locally) | 0.713 / 0.675 / 0.693 (0.480) | −0.012 |
| HG008T Illumina | v2.1, `--norm train`, 3 seeds | 0.439 / 0.832 / 0.575 mean (AP 0.310-0.327) | 0.398 / 0.783 / 0.528 (0.089) | −0.047 vs mean |
| HG008T PacBio | v2.1 | 0.479 / 0.827 / 0.607 (0.322-0.353) | 0.466 / 0.819 / 0.594 (0.075) | −0.012 |
| HG008T ONT | v2.1 | 0.457 / 0.831 / 0.589 (0.328-0.370) | 0.447 / 0.841 / 0.584 (0.071) | −0.006 (inside seed range 0.575-0.598) |
| COLO829T Illumina | v2.1 | 0.797 / 0.712 / 0.751 (0.818-0.824) | 0.362 / 0.678 / 0.472 (0.582) | −0.279 |
| COLO829T Fiber-seq | v2.1 | 0.854 / 0.719 / 0.780 (0.826-0.829) | 0.883 / 0.685 / 0.771 (0.592) | −0.009 |
| COLO829T ONT | v2.1 | 0.755 / 0.704 / 0.729 (0.741-0.755) | 0.703 / 0.649 / 0.675 (0.484) | −0.053 |
| (normbq v2.1) | HG008T Ill / PacBio / ONT, COLO829T Ill / Fiber-seq | with-ch6 means 0.573 / 0.606 / 0.596 / 0.751 / 0.779 | 0.529 / 0.596 / 0.584 / 0.477 / 0.732 | −0.044 / −0.011 / −0.012 / −0.274 / −0.047; COLO829T ONT no-ch6 normbq MISSING |

(ΔF1 of the v2.1 rows = no-ch6 minus the with-ch6 three-seed mean, from 4-decimal report values. Per-seed values: `q5_ch6.md`. With-ch6 v2.1 per-seed AP: HG008T Ill 0.310/0.327/0.327, PacBio 0.353/0.322/0.340, ONT
0.328/0.356/0.370; COLO829T Ill 0.818/0.824/0.818, Fiber-seq 0.829/0.828/0.826, ONT 0.755/0.741/0.748.)
- SNV F1 lower without ch6 in 11/12 seed-0 pairs (`--norm train`; RISSUM tables.md "No ch6": d9 −0.043/−0.041/−0.033/−0.019/−0.006/
  −0.012, d46 vs seed 0 −0.042/−0.011/+0.009/−0.311/−0.007/−0.062; median −0.026). Against the d46 three-seed mean 6/6 lower
  (HG008T ONT −0.006, within its seed range). AP drops on every pair (HG008T 0.199-0.370 → 0.065-0.089; COLO829T d46 0.741-0.829 →
  0.484-0.592).
- The d9 no-ch6 controls are 1-GPU runs (`_7fc7396` / `_h90`) except HG008T Illumina, not the 4-GPU matched recipe of section 1.
- INDEL (best F1, one seed): d9 −0.034 to +0.012, d46 −0.051 to +0.013 (RISSUM tables.md "No ch6"); 9 / 12 lower (RISSUM README: 20 of 24 SNV+INDEL pairs lower).
- no-ch6 − GRCh38-only three-seed mean (`--norm train`): d9 +0.038 / +0.016 / +0.039 / +0.299 / +0.084 / +0.239, d46 +0.024 / +0.030 /
  +0.061 / +0.011 / +0.080 / +0.221 → +0.011 to +0.299; two of 12 inside the GRCh38-only seed range (HG008T PacBio d9 0.580 in
  0.513-0.592; COLO829T Illumina d46 0.472 in 0.386-0.530).

## 6. Cross-node (path-following) context ablation

**Existing controlled ablation before today: MISSING (now run, see 6.1).** No experiment compares node-confined with path-following windows on the same
alignments and candidates (searched `analysis/reviewer_response_20261011/experiments/`, `supplementary/`, `machine_learning/
pansoma_net_v2/` options: only `--drop-planes` exists, no column mask). What exists: E1.7 (`experiments/E1.7_window_node_span/
README.md`, `window_span.tsv`) measures how often the ±50-column window reaches beyond the candidate's node (somatic SNV 56.8 % /
75.5 % / 71.9 %, somatic INDEL 91.7 % / 96.3 % / 97.1 %; **only three sets**: HG008T PacBio d9 and d46, COLO829T Fiber-seq d46; the
letter's "57–76 %" / "92–97 %" come from these three), and Supplementary Fig. R1.m2 (an example). The v1 (node-confined) vs v2
comparison differs in pipeline, candidates, labels and model, so it is not this ablation.

**Can the off-node columns be identified from stored data?** The tensors hold no per-column node ID (8 channels:
`indexed_gam_pipeline_v4/README.md` lines 992-1003; per-column node coordinates are written only by a `--debug-rows` build,
README line 253). But the own-node visit of every cell is reconstructible from the tensor plus three per-tensor numbers already
stored: channel 5 (graph base; 1-5 = a cell that consumes one graph base, 6 = gap/insertion), channel 4 (covered), the summary's
`start` and `site_layout.insertion_slots` (`<chrom>_variant_summary.ndjson`) and the node length (graph-index sqlite named in the
set's `manifest.json`). Rule (per row, rows are in the candidate node's forward orientation): left of the anchor column a cell is
on the own node iff the number of graph-consuming cells between it and the anchor <= `start`; the insertion slots are own; right of
the site a cell is own iff `start` + consumed/empty cells before it < node length. Implemented in
`tmp/assessment_20261030/ablation/q6/node_mask.py` (`own_node_mask`). **Checked today against the per-column node coordinates of
small `--debug-rows` builds of the pipeline's synthetic test inputs** (`q6/g8work/G{2,3,4a,5,7,8}-python/out/out/{SNV,INDEL}`):
0 disagreements in 17,008 covered cells, of which 5,124 lie off the own node visit (G7 INDEL 45/72, G8 SNV 2,668/8,571, G8 INDEL
2,411/7,087). **Not yet checked on a real tensor set** (needs a `--debug-rows` rebuild of a small region; MISSING). The index cache
(`data.py` KindIndex) does not keep `start` / slots / node length, so a side array (3 ints per tensor, in summary = index order)
must be built once per set.

**Masking ablation feasibility.**
- Representation: "a cell without evidence is 0 in every channel" (README line 992; `encode.py` lines 3-5, 17-19: every plane is 0
  where there is no evidence, `covered = ch4 > 0`). Zeroing all 8 channels of the off-node cells makes them identical to
  no-coverage cells, including the derived planes (differs, site-block planes) and ch6; the site-level scalars stay. So "zero" is
  consistent with the encoding.
- Where: (i) `data.py` `TensorDataset.__getitem__` (line ~454, after the `preadv`) for validation and predict (`predict.py` line 228
  builds a `TensorDataset`; the off-reference rescue subclass `OffrefRescue`, line 148, inherits it); (ii) `chunks.py`
  `ChunkLoader._pack` (line 241) for training batches; or once on the GPU in `encode.py` `forward` (line 115) from ch4/ch5 plus the
  3-int side array passed through the dataset like `blocks`.
- (a) **Inference only, existing checkpoints** (FEASIBLE, small): add a `--mask-off-node` flag to predict; re-predict chr1 with the
  12 HPRC SNV seed-0 checkpoints (and INDEL if wanted), then `jobs/vcf.sh`. Command pattern (as the tequila re-predictions,
  `tmp/reviewer_cpu_20261013/task_gpu.sh` lines 20-27): `python -m pansoma_net_v2.predict --checkpoint <run>/best.pth --tensors
  /scratch/jshen/data/pansoma_v2_tensors/<set>/tensors --kinds SNV --chroms chr1 --output <out>/test_chr1 --cache-dir
  <runs>/index_cache --num-workers 6 --batch-size 256 [--mask-off-node]` then `SET_NAME=<set> bash machine_learning/pansoma_net_v2/
  jobs/vcf.sh <run> SNV <out>`. chr1 predict wall time on tequila GPU was 50-2,021 s per run (`tmp/reviewer_cpu_20261013/logs/
  gpu_predict_vcf*-*.out`, "predict N s" lines, e.g. COLO829T ONT d46 SNV 2,021 s, HG008T Illumina d9 SNV 1,254 s). Caveat: an
  inference-only mask puts the model outside its training distribution (it never saw truncated rows), so a drop does not isolate
  the value of the context; it only bounds how much the trained models use it. Also apply the same mask on the GRCh38-only arm
  (32-bp nodes: nearly every column is off-node) only as a sanity check.
- (b) **Retrain with the mask in the loader** (the clean test; FEASIBLE, GPU): same recipe, mask applied in training, validation
  and predict. Training compute (sum of `metrics.jsonl` `train_seconds` x GPUs, one seed; `q6/gpu_hours.tsv` and checked today):
  v2.1 SNV six models 13.2 GPU-h (HG008T Ill 3.26, PacBio 1.44, ONT 1.39; COLO829T Ill 4.64, Fiber-seq 0.97, ONT 1.54; 4 x H100);
  v1.1 matched SNV 14.5 GPU-h; INDEL (1 GPU) v2.1 3.2, v1.1 3.7 GPU-h. Validation and data-loading time are not in
  `train_seconds`, so wall time is longer. Three seeds triple it. No GPU job was launched.

### 6.1 Update 2026-10-30 evening: masking ablation RUN (user asked "为啥还不运行")

**Implementation (isolated copy, no tracked file changed):** `tmp/assessment_20261030/ablation/q6/code/pansoma_net_v2` =
repo `machine_learning/pansoma_net_v2` at 90d7c0f plus `offnode.py` and a `--mask-off-node DIR` flag in `train.py`,
`predict.py`, `data.py` (`TensorDataset`), `chunks.py` (`ChunkLoader._pack`). Every cell outside the candidate's own node
visit is zeroed in all 8 channels; the site columns are always kept. Per-tensor inputs (node, start, slots, span, node length) for
the 12 HPRC SNV sets: `q6/side/<SET>.SNV.side.npz` (`side/build_side.py`; regexes spot-checked against full JSON parses;
node order identical to the index caches except COLO829T Illumina d9, whose local cache predates its 2026-10-02 rebuild,
so `offnode.side_arrays` checks the order against the live index at run time). Unmasked runs with the same code reproduce
the existing chr1 F1 exactly (12/12).

**Rule checked on real data (COMPUTED NOW):** `--debug-rows` builds of ~2,000 chr1 nodes (all somatic-tensor nodes + a
strided sample) and comparison with the per-column node coordinates (`q6/realcheck/`; Slurm jobs 462940, 462951).
- HG008T PacBio v2.1: SNV 149 disagreements in 11.9 M covered cells (48.1 % of cells off the own node); INDEL 27 in 7.8 M
  (65.8 % off).
- HG008T Illumina v1.1: SNV 2,088 in 24.6 M (25.7 % off); INDEL 7,996 in 7.7 M (0.10 %; 62.8 % off), not counting the
  site-span cells, which are kept by design.
- Every disagreement is a pure insertion cell exactly at a node junction, which vg assigns to one or the other neighbouring
  mapping (`q6/realcheck/categorize.py`).

**(a) Inference only (COMPUTED NOW)**
- Runs: seed-0 `--norm train` checkpoints, chr1, same tensors and candidates; tequila 12-GB GPU array 462953, 24 tasks
  (12 masked + 12 unmasked).
- Outputs: `/scratch/jshen/data/pansoma_net_v2_runs/offnode_mask_20261030/<run>__{mask,nomask}/`; table from `q6/infer/collect.py`
  -> `q6/infer/results.md`.
- Metrics: P / R / F1 after PoN at the checkpoint's PASS threshold, BED; (test-tuned best F1; chr1 per-tensor somatic AP).

| set | graph | unmasked P / R / F1 (best; AP) | masked P / R / F1 (best; AP) | ΔF1 PASS | Δbest F1 |
|---|---|---|---|---|---|
| HG008T Illumina | v2.1 | 0.429 / 0.852 / 0.570 (0.613; 0.310) | 0.477 / 0.821 / 0.604 (0.605; 0.313) | +0.033 | −0.008 |
| HG008T PacBio | v2.1 | 0.472 / 0.836 / 0.604 (0.622; 0.352) | 0.482 / 0.679 / 0.563 (0.597; 0.311) | −0.040 | −0.025 |
| HG008T ONT | v2.1 | 0.432 / 0.859 / 0.575 (0.603; 0.328) | 0.469 / 0.815 / 0.595 (0.600; 0.326) | +0.020 | −0.003 |
| COLO829T Illumina | v2.1 | 0.867 / 0.714 / 0.783 (0.786; 0.818) | 0.872 / 0.597 / 0.709 (0.733; 0.748) | −0.074 | −0.053 |
| COLO829T Fiber-seq | v2.1 | 0.866 / 0.707 / 0.779 (0.783; 0.829) | 0.871 / 0.598 / 0.709 (0.752; 0.767) | −0.069 | −0.031 |
| COLO829T ONT | v2.1 | 0.787 / 0.692 / 0.736 (0.739; 0.755) | 0.739 / 0.600 / 0.662 (0.669; 0.688) | −0.074 | −0.069 |
| HG008T Illumina | v1.1 | 0.452 / 0.831 / 0.585 (0.607; 0.199) | 0.492 / 0.798 / 0.609 (0.610; 0.188) | +0.024 | +0.002 |
| HG008T PacBio | v1.1 | 0.478 / 0.829 / 0.607 (0.630; 0.199) | 0.484 / 0.732 / 0.583 (0.588; 0.176) | −0.024 | −0.042 |
| HG008T ONT | v1.1 | 0.473 / 0.821 / 0.600 (0.616; 0.199) | 0.486 / 0.752 / 0.591 (0.597; 0.189) | −0.009 | −0.019 |
| COLO829T Illumina | v1.1 (epoch 8) | 0.484 / 0.698 / 0.572 (0.660; 0.736) | 0.531 / 0.663 / 0.590 (0.631; 0.692) | +0.017 | −0.028 |
| COLO829T Fiber-seq | v1.1 | 0.785 / 0.723 / 0.752 (0.774; 0.721) | 0.776 / 0.676 / 0.723 (0.726; 0.646) | −0.030 | −0.048 |
| COLO829T ONT | v1.1 | 0.793 / 0.701 / 0.744 (0.745; 0.626) | 0.834 / 0.610 / 0.705 (0.718; 0.570) | −0.039 | −0.027 |

Summary of (a):
- Recall falls in 12/12 (−0.031 to −0.157). Test-tuned best F1 falls in 11/12 (−0.003 to −0.069). Somatic AP falls in 11/12
  (up to −0.075).
- PASS F1 falls in 8/12 (−0.009 to −0.074) and rises in 4 (HG008T Illumina both graphs, HG008T ONT v2.1, COLO829T Illumina
  v1.1). The rises come from higher precision at the fixed threshold, as fewer candidates pass.
- Caveat: the models never saw masked tensors, so (a) shows that the trained models use the cross-node columns; it does not
  measure what node-confined training would achieve.

**(b) Retraining with the mask: NOT RUN, waiting for the user's choice of platform.**
- Script ready: `q6/retrain/train_tequila.sh`, 2 x H100 with global batch 4096 = 2 x 2048, 2 micro-batches; `MODE=nomask`
  for a same-setup control.
- RIS (4 x H100, identical to the controls) needs the q6 code copy and `side/` on RIS storage.
- Cost: d46 SNV x 6, one seed, about 13 GPU-h of training on H100.

## 7. Seeds of the cross-tumor (multi-tumor) models (EXISTS)

- **SNV multicancer4 / multicancer6: one seed per platform** (runs `multicancer{4,6}_{Illumina,PacBio,ONT}_SNV` in
  `/scratch/jshen/data/pansoma_net_v2_runs/`, `args.json` `seed` 0, no `_s1`/`_s2` runs; E12b README "Design ... one seed").
  Training chromosomes: HCC1395T without chr1/chr20, **HCC1954T, H1437T, H2009T on all autosomes including chr1**
  (`set_test_chroms` only HCC1395T; `test_chroms` []); multicancer6 adds HG008T/COLO829T without chr1. The single-tumor cross-tumor
  runs (E12a, `xs__*`) are also one seed.
- **INDEL final batch (`tmp/indel_gap_20261020/final_20261023/results_ppr8.md`): three seeds** for both the tumor-unseen (m4:
  HCC1395T, HCC1954T, H1437T, H2009T; runs `unseen_s{0,1,2}_{Illumina,PacBio,ONT}_INDEL`) and release (m6: + HG008T, COLO829T;
  `release_s{0,1,2}_*`) models (PREREGISTRATION.md lines 32-33); here chr1 is held out of **every** tumor (and chr20 of HCC1395T).
  Recipe: 8.7 M, sqrt, 12 epochs, site features, `--best-by postpon`, `--norm dataset`, non-reasons, PASS ppr8.

## 8. Letter statements inconsistent with these numbers (letter = `scratchpad/assemble/response_letter_20261015.md`)

Corrections assume the rule "validation-selected checkpoint for every model" (COLO829T Illumina v1.1 at epoch 8).

**R1.1**
- Point 2, retrained no-ch6 AP "from 0.199–0.353 ... (COLO829T, v2.1: from 0.755–0.829 ...)": these with-ch6 ranges are seed 0
  only; over the three v2.1 seeds HG008T with-ch6 AP is 0.310–0.370 (all HPRC with-ch6 HG008T runs 0.199–0.370) and COLO829T v2.1
  0.741–0.829 (section 5). The no-ch6 side (0.065–0.089; 0.484–0.592) is correct.
- Point 2, "SNV F1 after population filtering falls in 11 of 12 models (median 0.03, at most 0.31)": correct (median −0.026), but
  the six v1.1 pairs use the 1-GPU controls (`HG008_{PacBio,ONT}_SNV_7fc7396`, `COLO829_*_SNV_h90`), not the matched 4-GPU recipe
  of the main table; say so or compare to the matched runs.
- Point 2, inference-time ch6 ablation "12 HPRC SNV models": the three COLO829T v1.1 models there are the purity-joint models
  `COLO829_*_SNV_joint_7fc7396_gpu4` (E1.1 README line 35), not the P100 models of the main table.
- Point 3 "+0.028 to +0.292": correct with validation-selected checkpoints (`--norm train`).

**R1.2**
- "Its gain is mainly before population filtering (SNV F1 +0.079 to +0.154, through precision)": with COLO829T Illumina v1.1 at
  the validation-selected epoch 8: **+0.079 to +0.222** (`--norm train`; precision +0.089 to +0.280) or +0.088 to +0.224 (normbq).
- "afterwards the two graphs are equal on all six datasets (... PASS F1 differ by −0.015 to +0.021 and ... best F1 by at most
  0.012; COLO829T Illumina v1.1 at epochs 5–6)": validation-selected: **−0.015 to +0.164** (normbq) or −0.018 to +0.163
  (`--norm train`); best F1 up to **+0.092** (COLO829T Illumina). Not equal on COLO829T Illumina. Even on the other five sets the
  seed ranges do not overlap for HG008T Illumina in normbq (v2.1 below v1.1 on all three seeds, −0.015) and for COLO829T ONT in
  `--norm train` (v2.1 below on all three seeds, −0.018).
- "indel F1 differs by at most 0.006": with the current matched runs (`COLO829_*_INDEL_7fc7396`, 71cb51a PoN) the r0.9 PASS
  difference is −0.003 to **+0.013** (COLO829T Illumina 0.104 → 0.117), test-tuned best −0.013 to +0.015, valbest −0.036 to
  +0.052 (section 1.4). (0.006 was HG008T only, R1.2.md evidence line 50.)
- "and indel F1 is also higher" (than GRCh38-only): true at the validation best-F1 threshold (12/12, few calls on COLO829T);
  at the default recall-0.9 PASS GRCh38-only is higher on the four long-read datasets (8 of 12 pairs). State the threshold.
- "single-run paired block-bootstrap intervals exclude zero in all twelve pairs": still true with validation-selected matched
  seed-0 runs (ΔF1 +0.038 to +0.326, section 4.2), but the current Table R3.9b behind it uses COLO829T Illumina v1.1 at epoch 6
  and the 1-GPU HG008T PacBio/ONT v1.1 runs. Seeds 1 and 2 give 10/12 and 11/12 (33/36 overall, none significantly negative).

**R3.1**
- "everything downstream is identical": the training recipe and selection rules are identical in `--norm train` (section 1.0),
  but the validation regions are different GRCh38 intervals in each graph (Jaccard 0.014–0.034, section 3), so the numeric PASS
  thresholds (section 2) and the selected epochs differ; ch6 is coded against H = 17/90/464; Giraffe versions differ.
- Point 2 "Before filtering, v2.1 is higher ... (SNV F1 +0.079 to +0.154), mostly through precision (+0.089 to +0.179)":
  validation-selected **+0.079 to +0.222, precision +0.089 to +0.280** (`--norm train`).
- Point 2 "After filtering ... the two graphs are equal within the seed range on all six datasets (means differ by −0.015 to
  +0.021). COLO829T Illumina needs v1.1 at epochs 5–6, chosen on chr1": validation-selected **−0.018 to +0.163** (`--norm train`)
  / −0.015 to +0.164 (normbq); "within the seed range on all six" is false (see R1.2 above). The chr1-chosen epochs contradict
  R3.9 point 3 ("PASS is set on the held-out validation blocks, never on chr1") and R3.m1 "Checkpoint".
- Point 3 "Retrained without the haplotype-count channel, both HPRC graphs still beat GRCh38-only (+0.011 to +0.299)": numbers
  correct (no-ch6 single seed vs GRCh38-only three-seed mean), but 2 of 12 lie inside the GRCh38-only seed range (HG008T PacBio
  v1.1 0.580 vs 0.513–0.592; COLO829T Illumina v2.1 0.472 vs 0.386–0.530) and the v1.1 no-ch6 runs are 1-GPU.
- Point 1 numbers (0.586–0.768 / 0.575–0.780 / 0.454–0.691; P and R higher on all six; INDEL valbest 12/12; r0.9 GRCh38-only higher
  on four) are all correct with validation-selected checkpoints.

**R3.9**
- Point 1 "recall SD is 0.006–0.019, precision SD 0.003–0.083 and F1 SD 0.001–0.040 (median 0.012)": these are the normbq rows with
  COLO829T Illumina v1.1 at epochs 6/5/6. With epoch 8: **recall SD 0.004–0.019** (P and F1 unchanged). With the `--norm train` table
  (R3.9c): P 0.005–0.084, R 0.003–0.026, F1 0.001–0.042 (median 0.012).
- Point 1 "seeds ... also determine the held-out validation blocks": correct (SHA-1 of seed:chrom:block); add that the three seeds'
  validation regions overlap by only 2–10 % and the three graphs' by 2–6 % (section 3), so seed SD = joint initialisation + split.
- Point 1 "35 of 36": correct (validation-selected).
- Point 3 "PASS is set on the held-out validation blocks, never on chr1": the threshold yes, but R3.9a and R3.9b (and R1.2/R3.1)
  use COLO829T Illumina v1.1 checkpoints chosen on chr1 (epochs 6/5/6, or epoch 6 for seed 0 in R3.9b).
- Supplementary Table R3.9a and R3.9c: all per-seed P/R/F1 values match the run files (sections 1.1, 1.3).
- Point 2 "For HPRC versus GRCh38-only, ΔF1 is +0.038 to +0.326 and excludes zero in all 12 SNV comparisons; Δprecision excludes
  zero in 11 and Δrecall in eight": with validation-selected matched seed-0 runs the ΔF1 range and 12/12 hold, **Δprecision excludes
  zero in 10** (COLO829T Illumina v1.1 epoch 8: +0.094 [−0.068, +0.158]), Δrecall 8.
- Table R3.9b caption "COLO829T Illumina on v1.1 at epoch 6": rows to replace — COLO829T Illumina v1.1 − GRCh38-only +0.104 [+0.022,
  +0.151] (ΔP +0.094 [−0.068, +0.158], ΔR +0.115 [+0.095, +0.135]) instead of +0.306 [+0.107, +0.441]; HG008T PacBio v1.1 (matched)
  +0.094 [+0.031, +0.174] instead of +0.108 [+0.046, +0.187]; HG008T ONT v1.1 (matched) +0.064 [+0.037, +0.095] instead of +0.059
  [+0.034, +0.091] (sections 4.1-4.2; full P/R CIs in `q4/e16_pairs_q4.md`). Pansoma-vs-comparator rows use v2.1 seed 0 and are
  unaffected.
- Point 2 should say the intervals are per single seed (seed 0); per-seed intervals for seeds 1-2 now exist (section 4.2).

**R3.m1 (table)**
- "Checkpoint: the epoch with the best validation somatic F1": true for every SNV run, but not for the COLO829T Illumina v1.1
  epochs 6/5/6 where these are used, and the final INDEL models use `--best-by postpon` (best post-PoN validation F1, with site
  features).
- "Calling: PASS at the recall-0.9 validation threshold; INDELs also at the validation best-F1 threshold": the final INDEL models
  (R3.10, R1.6) call at ppr8 (post-PoN validation recall 0.8 of the callable positives, `graph_vcf --threshold postpon-recall
  --target-recall 0.8`); r0.9 / valbest apply to the 7fc7396 INDEL ablation runs only.
- "Labels and split: chr2–22 for training ...; chr1 of every evaluated tumor held out": the multi-tumor SNV models train HCC1954T,
  H1437T and H2009T on all autosomes including chr1 (`args.json` `set_test_chroms` = HCC1395T chr1, chr20 only); the final INDEL
  models hold out chr1 of every tumor. Also say how the 5 % is drawn (seeded hash of node-ID blocks; different GRCh38 regions per
  graph and per seed).
- "Training": the INDEL row omits the final INDEL recipe's site features, `--norm dataset`, non-reasons; the input normalisation
  differs between R3.9a (normbq) and R3.9c (`--norm train`) and is not stated.

## 9. Secondary analyses using COLO829T Illumina v1.1: epoch 6 vs validation-selected epoch 8

Runs: epoch 8 = `COLO829_Illumina_SNV_7fc7396_ddp4_b1024` (`--norm train`, seed 0; `train.log` "done; best epoch 8") or
`COLO829_Illumina_SNV_d9_normbq` (normbq); epoch 6 = `COLO829_Illumina_SNV_d9_normbq_sweep_e06` (chr1-chosen). New outputs under
`tmp/assessment_20261030/ablation/q9/` (Slurm job 462881: E3b/E8 split, E16 COLO pairs, E19 ep8/ep6) and `q9b/strata_ep8.py`.

| analysis | COLO829T Illumina v1.1 run used now | epoch-8 numbers | status |
|---|---|---|---|
| E10 SMaHT strata, call counts (`experiments/E10_region_strata/tables.md` §2) | epoch 8 (`ddp4`) in tables.md; the letter (R1.3, R2.5) substitutes epoch 6 (`colo829t_illumina_d9_epochs.tsv`: extreme 243 calls, P 0.593) | extreme: 1,725 calls, P 0.076, R 0.560, F1 0.134 (tables.md row "COLO829T Illumina SNV d9"). HPRC arms then make **102–1,725 extreme PASS calls at precision 0.076–0.764** (letter: 102–382 at 0.161–0.764); GRCh38-only 351–2,210 at 0.056–0.253 unchanged. Precision still higher than GRCh38-only in 12/12 (0.076 vs 0.056) | EXISTS |
| E16 SMaHT strata paired (`E16_bootstrap/strata_paired.tsv`) and E10b GIAB strata (`E10b_giab_strata/strata_paired.tsv`) | epoch 8 (`ddp4`) in both tsv; `E10b_giab_strata/summary_ep6.py` swaps in epoch 6 (`colo_ill_d9_ep6.json`) to make `summary_ep6.md`, which the letter quotes | recomputed the same summary without the swap (`q9b/strata_ep8.md`, which reproduces `summary_ep6.md` exactly in its ep6 columns). Extreme: ΔF1 +0.033 to +0.531, **CI above 0 in 10/12** (ep6: 11/12; COLO829T Illumina v1.1 +0.033 [−0.074, +0.091] vs ep6 +0.503 [+0.213, +0.588]); difficult −0.005 to +0.259, **> 0 in 10/12** (ep6 11/12; COLO829T Illumina v1.1 −0.001 [−0.080, +0.065]), CI above 0 in 3; easy unchanged (−0.003 to +0.154, 6 above 0). GIAB all difficult regions +0.054 to +0.483, **10/12** above 0 (ep6 11/12); low mappability +0.034 to +0.548, 10 (ep6 11); segmental duplications −0.063 to +0.331, 9 (ep6 10; COLO829T Illumina v1.1 now −0.063 [−0.140, +0.142]); tandem repeats + homopolymers +0.031 to +0.513, 11 (ep6 12); satellites +0.002 to +0.490, 5; never significantly below zero in any stratum (both versions) | EXISTS (tsv) / COMPUTED NOW (summary) |
| E19 germline density (R2.5) (`E19_germline_density/tables.md`) | **no COLO829T Illumina v1.1 row at all** (COLO829T d9 single-set runs were not local on 10-15; only the Fiber-seq joint model); R2.5's density numbers are HG008T plus COLO829T v2.1 | computed now with E19's script unchanged (`q9/e19_colo_ill_d9_q9.py`, `q9/E19_ep8/tables.md`, `q9/E19_ep6/tables.md`): W = 50 bp, F1 by 0 / 1 / 2+ nearby germline variants: epoch 8 0.557 / 0.719 / 0.563 (FP 1,730 / 58 / 27; F1 all 0.572) vs epoch 6 0.786 / 0.737 / 0.565 (FP 263 / 64 / 29); GRCh38-only 0.451 / 0.648 / 0.524, v2.1 0.790 / 0.753 / 0.684 | COMPUTED NOW (not used in letter) |
| E19b dense-germline reads (`E19b_dense_germline_reads/README.md` lines 57-59) | COLO829T columns are "- / v2.1" (no v1.1) | not applicable; "94–97 % of 23–49 such false calls removed" is v2.1 only | unaffected |
| E3b complete-output F1 (`E3b_.../tables.md` §1) | epoch 8 (`ddp4`, seed 0) already | rescue on 186 graph-coordinate PASS / F1 0.5546; off 53 / 0.5668 (EXISTS). Epoch 6 for comparison (`q9/e3b_e8_q9.md`): on 693 / 0.6715, off 38 / 0.7671. Seeds 1/2 epoch 8 (rescue on only): 0.5668 / 0.5029 (`ddp4`), normbq 0.5484 / 0.5693 / 0.5048 | EXISTS / COMPUTED NOW |
| E16 arm pairs (`E16_bootstrap/paired.tsv`, letter Table R3.9b) | epoch 6 (`sweep_e06`, E16 ARM_RUN) | v1.1 − GRCh38-only +0.104 [+0.022, +0.151] (ΔP +0.094 [−0.068, +0.158], ΔR +0.115 [+0.095, +0.135]); v2.1 − v1.1 +0.211 [+0.061, +0.340]; all seeds in section 4.2 | COMPUTED NOW |
| E8 VAF split (`E8_comparators/latest_models.md` line 36) | epoch 6 | recall at truth VAF < 0.08: 0.067 (405) both; VAF >= 0.08: 0.823 (2,038) at epoch 8 vs 0.853 at epoch 6 (`q9/e3b_e8_q9.md`) | COMPUTED NOW |
| E8 `latest_models.md` comparator table (Pansoma COLO829T Illumina v1.1 seed range) | epochs 6/5/6 (line 3, line 12) | epoch-8 seed range = section 1.3 (P 0.476–0.565, R 0.698–0.706, F1 0.566–0.628, normbq) | EXISTS (RISSUM) |
| E9a projection (`E9a_unplaced_projection/results/*`) | both: `ddp4` (epoch 8) rows and a separate epoch-6 projection (README line 465) | in `results/tables.md` / `metrics_by_run.tsv` | EXISTS |
| E9c off-ref filter, E9i normal reads, E9j off-ref COLO829T, E9k DSA somatic off-ref (R3.3 / R3.4) | epoch 6 only (`E9c/tables.md` lines 14, 42; `E9i/README.md` line 11; `E9j/README.md` line 20; `E9k/README.md` line 13) | no epoch-8 row | MISSING (would need E9c's pipeline rerun on the `ddp4` run; not attempted) |
| E21 (`E21_training_composition/README.md` line 669) | "selected epoch 6" refers to an E13 199.6 M model, not COLO829T Illumina | — | unaffected |

Letter consequences (with epoch 8): R1.3 "significantly in 11 of 12 in both [GIAB union of difficult regions, SMaHT extreme]"
→ 10 of 12 in both; "significantly in 10–12 of 12 pairs in low-mappability, segmental-duplication and tandem-repeat regions" →
9–11; "rises in all 12 ... in the GIAB union of difficult regions and in SMaHT's extreme regions" still true (+0.054, +0.033 for
COLO829T Illumina v1.1); "102–382 at 0.161–0.764" → 102–1,725 at 0.076–0.764; "the exception being COLO829T Fiber-seq on v1.1"
→ COLO829T Fiber-seq and COLO829T Illumina on v1.1. R2.5 "raise SNV precision from 0.056–0.253 to 0.161–0.764 in all twelve
comparisons" → to 0.076–0.764 (still all twelve); "in difficult regions F1 changes by −0.005 to +0.259 (higher in eleven of
twelve)" → higher in ten of twelve; "COLO829T Illumina on v1.1 at epochs 5–6" qualifier drops.
