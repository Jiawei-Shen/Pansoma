# E19: does the graph help more where the patient's germline variants cluster?

Reviewer point: **R2.5** (where the precision gain comes from); also R3.1 (reference ablation). Shared code:
`../call_level_common/calls_lib.py`.

**Status (2026-10-15 17:58, job 406570): every call set that exists locally is processed (93).**
- In-sample (45): SNV linear / d9 / d46 for HG008T Illumina / PacBio / ONT and linear / d46 for COLO829T Illumina /
  fiberseq / ONT (15); INDEL the same arms at the standard threshold and at `vcf_chr1_valbest` (30).
- COLO829T fiber-seq joint d9 model (8): SNV and INDEL x P100 / P75 / P50 / P25. Not a matched arm.
- Cross-sample `xs__` (40): 30 SNV (Illumina / PacBio-fiberseq / ONT) and 10 Illumina INDEL, all d46.
- Missing: COLO829T d9 single-set runs (not local), comparators (not re-scored), seed repeats (GPU); see
  "Still missing" after the interpretation.

## Question

A pangenome carries common germline alleles, so reads from a haplotype dense in germline variants should align
better to it than to GRCh38. Does the advantage of the graph arms (d9, d46) over GRCh38-only grow with the local
density of the normal's germline variants? And how do precision and recall depend on that density within each arm?

## Data

- **Call sets**: `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1*`, written by
  `machine_learning/pansoma_net_v2/jobs/vcf.sh`. Used: rtg `eval_bed/vcfeval/pon_calls/` (PASS calls at the run's
  threshold, after the GRCh38 PoN, in the benchmark BED):
  - `tp-baseline` + `fn` = the truth;
  - `tp` + `fp` = the calls.
  - Included: in-sample (`main`, `vcf_chr1_valbest`), joint-model purity sets and cross-sample (`xs__`) call sets.
    `__norescue` re-calls are left out: their placed calls equal the main run's (E4b, E10).
- **Germline truth of the normal** (the set's `labels.manifest.json`): HG008-N dipcall, COLO829BL dipcall, and the
  CASTLE normals' DeepVariant hybrid QUAL10 calls (HCC1395BL, H1437, H2009, HCC1954) for the cross-sample targets.
  Chr1 records with an ALT genotype, any FILTER.
  Pansoma is tumor-only; the normal is used only for this stratification.
- **Truth with a tensor**: the truth allele has an exact somatic match in the call set's own `chr1_labels.ndjson`
  (`calls_lib.truth_tensor_keys`).

## Method

1. Each rtg record is normalized (left-aligned, trimmed). Its reference interval is the substituted or deleted bases,
   or for an insertion the anchor and the next base.
2. Density n_W = the number of germline records that intersect [start − W, end + W) but not the variant's own
   interval, for W = 50 and 500 bp.
3. Strata 0 / 1 / 2+.
4. Recall per stratum = TP-baseline / truth of the stratum; precision = TP calls / calls of the stratum; F1 from the
   two.
5. `tensor share / R with tensor` = the share of the stratum's truth with a tensor, and recall among those. A low
   tensor share means the loss is in candidate generation (no tensor is built), not in the model.
6. **Paired arm differences**: the truth and the BED are the same in every arm of a sample / platform / kind, so each
   truth record has a stratum and a found / missed state per arm. Reported per stratum for d9 − linear,
   d46 − linear and d46 − d9:
   - dR, dP and dF1;
   - FP counts;
   - the discordant truth counts (found by A only / by B only), with an exact two-sided McNemar p.
7. Checks: the stratum totals equal rtg's TP-baseline and FP counts for every call set.

## Scripts and jobs

| file | role |
|---|---|
| `germline_density.py` | per-record densities -> `DATA/E19/<call set>.v2.tsv.gz`; `tables.md`, `summary.json` |
| `collect.sh` | runs `germline_density.py` |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E19_germline_density/job.sh` | Slurm wrapper (`-p general`, 4 CPUs, 12 GB); logs in `logs/` |

`DATA` = `/scratch/jshen/data/reviewer_experiments_20261013/call_level`. Job 406570 (all 93 call sets): 314 s, 1.44 GB max RSS (`/usr/bin/time` in the log).

## Results (all call sets, 2026-10-15 18:00; 93 call sets)

chr1, in the BED, PASS after the PoN; `tables.md` of 17:58. Arms: linear = GRCh38-only, d9 = HPRC v1.1 d9, d46 = HPRC
v2.1 d46. (*) = a pair with the HG008T PacBio / ONT SNV d9 checkpoints, trained with a different recipe (1 GPU; E3
caveat 1), so these pairs are not recipe-matched.

### 1. chr1 in-sample SNV per stratum, W = 50 bp (`tables.md` section 1)

| call set (chr1) | truth 0 / 1 / 2+ | FP 0 / 1 / 2+ | R 0 / 1 / 2+ | P 0 / 1 / 2+ | F1 0 / 1 / 2+ | F1 all | tensor share 0 / 1 / 2+ | R with tensor 0 / 1 / 2+ |
|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV linear | 2091 / 282 / 70 | 2163 / 35 / 23 | 0.592 / 0.539 / 0.471 | 0.364 / 0.813 / 0.589 | 0.451 / 0.648 / 0.524 | 0.468 | 0.835 / 0.819 / 0.843 | 0.709 / 0.658 / 0.559 |
| COLO829T Illumina SNV d46 | 2091 / 282 / 70 | 227 / 36 / 5 | 0.724 / 0.681 / 0.557 | 0.870 / 0.842 / 0.886 | 0.790 / 0.753 / 0.684 | 0.783 | 0.824 / 0.801 / 0.757 | 0.878 / 0.850 / 0.736 |
| COLO829T fiberseq SNV linear | 2091 / 282 / 70 | 440 / 48 / 33 | 0.637 / 0.564 / 0.471 | 0.752 / 0.768 / 0.500 | 0.690 / 0.650 / 0.485 | 0.679 | 0.828 / 0.816 / 0.771 | 0.770 / 0.691 / 0.611 |
| COLO829T fiberseq SNV d46 | 2091 / 282 / 70 | 228 / 37 / 3 | 0.715 / 0.684 / 0.557 | 0.868 / 0.839 / 0.929 | 0.784 / 0.754 / 0.696 | 0.779 | 0.824 / 0.809 / 0.786 | 0.868 / 0.846 / 0.709 |
| COLO829T ONT SNV linear | 2091 / 282 / 70 | 2386 / 140 / 49 | 0.550 / 0.507 / 0.443 | 0.325 / 0.505 / 0.388 | 0.409 / 0.506 / 0.413 | 0.418 | 0.824 / 0.823 / 0.771 | 0.668 / 0.616 / 0.574 |
| COLO829T ONT SNV d46 | 2091 / 282 / 70 | 384 / 62 / 11 | 0.699 / 0.681 / 0.514 | 0.792 / 0.756 / 0.766 | 0.743 / 0.716 / 0.615 | 0.736 | 0.821 / 0.816 / 0.786 | 0.852 / 0.835 / 0.655 |
| HG008T Illumina SNV linear | 561 / 97 / 39 | 693 / 160 / 236 | 0.832 / 0.732 / 0.256 | 0.403 / 0.307 / 0.041 | 0.543 / 0.433 / 0.070 | 0.470 | 0.989 / 0.979 / 0.872 | 0.841 / 0.747 / 0.294 |
| HG008T Illumina SNV d9 | 561 / 97 / 39 | 515 / 114 / 73 | 0.881 / 0.753 / 0.308 | 0.490 / 0.390 / 0.141 | 0.629 / 0.514 / 0.194 | 0.585 | 0.986 / 0.928 / 0.436 | 0.893 / 0.811 / 0.706 |
| HG008T Illumina SNV d46 | 561 / 97 / 39 | 590 / 136 / 66 | 0.897 / 0.814 / 0.308 | 0.460 / 0.367 / 0.154 | 0.608 / 0.506 / 0.205 | 0.570 | 0.982 / 0.948 / 0.590 | 0.913 / 0.859 / 0.522 |
| HG008T PacBio SNV linear | 561 / 97 / 39 | 645 / 144 / 180 | 0.877 / 0.753 / 0.282 | 0.433 / 0.336 / 0.058 | 0.580 / 0.465 / 0.096 | 0.514 | 0.995 / 0.979 / 0.769 | 0.882 / 0.768 / 0.367 |
| HG008T PacBio SNV d9 | 561 / 97 / 39 | 488 / 77 / 18 | 0.879 / 0.763 / 0.256 | 0.503 / 0.490 / 0.357 | 0.639 / 0.597 / 0.299 | 0.621 | 0.984 / 0.928 / 0.436 | 0.893 / 0.822 / 0.588 |
| HG008T PacBio SNV d46 | 561 / 97 / 39 | 525 / 96 / 29 | 0.882 / 0.773 / 0.308 | 0.485 / 0.439 / 0.293 | 0.626 / 0.560 / 0.300 | 0.603 | 0.980 / 0.948 / 0.513 | 0.900 / 0.815 / 0.600 |
| HG008T ONT SNV linear | 561 / 97 / 39 | 650 / 129 / 125 | 0.891 / 0.784 / 0.282 | 0.435 / 0.371 / 0.081 | 0.584 / 0.503 / 0.126 | 0.537 | 0.996 / 0.979 / 0.692 | 0.894 / 0.800 / 0.407 |
| HG008T ONT SNV d9 | 561 / 97 / 39 | 540 / 97 / 23 | 0.881 / 0.722 / 0.308 | 0.478 / 0.419 / 0.343 | 0.619 / 0.530 / 0.324 | 0.596 | 0.989 / 0.928 / 0.487 | 0.890 / 0.778 / 0.632 |
| HG008T ONT SNV d46 | 561 / 97 / 39 | 610 / 126 / 52 | 0.909 / 0.794 / 0.308 | 0.455 / 0.379 / 0.188 | 0.607 / 0.513 / 0.233 | 0.575 | 0.989 / 0.948 / 0.487 | 0.919 / 0.837 / 0.632 |

### 2. chr1 in-sample SNV per stratum, W = 500 bp (`tables.md` section 4)

| call set (chr1) | truth 0 / 1 / 2+ | FP 0 / 1 / 2+ | R 0 / 1 / 2+ | P 0 / 1 / 2+ | F1 0 / 1 / 2+ |
|---|---|---|---|---|---|
| COLO829T Illumina SNV linear | 821 / 635 / 987 | 2054 / 51 / 116 | 0.593 / 0.573 / 0.580 | 0.192 / 0.877 / 0.831 | 0.290 / 0.693 / 0.683 |
| COLO829T Illumina SNV d46 | 821 / 635 / 987 | 107 / 61 / 100 | 0.730 / 0.707 / 0.706 | 0.848 / 0.880 / 0.875 | 0.785 / 0.784 / 0.781 |
| COLO829T fiberseq SNV linear | 821 / 635 / 987 | 260 / 84 / 177 | 0.621 / 0.628 / 0.624 | 0.662 / 0.826 / 0.777 | 0.641 / 0.714 / 0.692 |
| COLO829T fiberseq SNV d46 | 821 / 635 / 987 | 125 / 53 / 90 | 0.706 / 0.713 / 0.704 | 0.823 / 0.895 / 0.885 | 0.760 / 0.794 / 0.784 |
| COLO829T ONT SNV linear | 821 / 635 / 987 | 1754 / 296 / 525 | 0.529 / 0.561 / 0.542 | 0.198 / 0.546 / 0.505 | 0.288 / 0.553 / 0.523 |
| COLO829T ONT SNV d46 | 821 / 635 / 987 | 161 / 111 / 185 | 0.695 / 0.694 / 0.687 | 0.780 / 0.799 / 0.786 | 0.735 / 0.743 / 0.733 |
| HG008T Illumina SNV linear | 192 / 182 / 323 | 291 / 202 / 596 | 0.839 / 0.819 / 0.737 | 0.356 / 0.425 / 0.285 | 0.500 / 0.559 / 0.411 |
| HG008T Illumina SNV d9 | 192 / 182 / 323 | 194 / 178 / 330 | 0.901 / 0.868 / 0.768 | 0.471 / 0.470 / 0.429 | 0.619 / 0.610 / 0.550 |
| HG008T Illumina SNV d46 | 192 / 182 / 323 | 222 / 199 / 371 | 0.901 / 0.874 / 0.811 | 0.438 / 0.444 / 0.414 | 0.589 / 0.589 / 0.548 |
| HG008T PacBio SNV linear | 192 / 182 / 323 | 246 / 180 / 543 | 0.859 / 0.890 / 0.771 | 0.401 / 0.474 / 0.314 | 0.547 / 0.618 / 0.447 |
| HG008T PacBio SNV d9 | 192 / 182 / 323 | 180 / 155 / 248 | 0.875 / 0.868 / 0.777 | 0.483 / 0.505 / 0.503 | 0.622 / 0.638 / 0.611 |
| HG008T PacBio SNV d46 | 192 / 182 / 323 | 195 / 175 / 280 | 0.901 / 0.835 / 0.796 | 0.470 / 0.465 / 0.479 | 0.618 / 0.597 / 0.598 |
| HG008T ONT SNV linear | 192 / 182 / 323 | 262 / 206 / 436 | 0.885 / 0.863 / 0.805 | 0.394 / 0.433 / 0.374 | 0.545 / 0.576 / 0.510 |
| HG008T ONT SNV d9 | 192 / 182 / 323 | 201 / 177 / 282 | 0.896 / 0.874 / 0.759 | 0.461 / 0.473 / 0.465 | 0.609 / 0.614 / 0.576 |
| HG008T ONT SNV d46 | 192 / 182 / 323 | 227 / 206 / 355 | 0.911 / 0.879 / 0.817 | 0.435 / 0.437 / 0.426 | 0.589 / 0.584 / 0.561 |

### 3. chr1 paired arm differences, SNV (`tables.md` sections 7 and 8)

Recall discordance (found by A only / by B only) and its McNemar p are over all strata; per-stratum discordances are
in `tables.md`.

| set (chr1) | A - B | dF1 0 / 1 / 2+ / all (W 50) | dP 0 / 1 / 2+ (W 50) | dR all | found A only / B only (McNemar p) | FP A / B in 2+ (W 50) | dF1 0 / 1 / 2+ (W 500) |
|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | d46 - linear | 0.339 / 0.105 / 0.160 / 0.316 | 0.506 / 0.029 / 0.297 | 0.132 | 396 / 74 (<0.001) | 5 / 23 | 0.495 / 0.091 / 0.098 |
| COLO829T fiberseq SNV | d46 - linear | 0.094 / 0.104 / 0.211 / 0.099 | 0.116 / 0.071 / 0.429 | 0.083 | 297 / 94 (<0.001) | 3 / 33 | 0.119 / 0.080 / 0.092 |
| COLO829T ONT SNV | d46 - linear | 0.334 / 0.210 / 0.202 / 0.319 | 0.467 / 0.251 / 0.378 | 0.149 | 466 / 101 (<0.001) | 11 / 49 | 0.447 / 0.190 / 0.210 |
| HG008T Illumina SNV | d9 - linear | 0.087 / 0.081 / 0.123 / 0.116 | 0.087 / 0.083 / 0.101 | 0.044 | 57 / 26 (<0.001) | 73 / 236 | 0.119 / 0.051 / 0.139 |
| HG008T Illumina SNV | d46 - linear | 0.066 / 0.073 / 0.135 / 0.101 | 0.058 / 0.060 / 0.113 | 0.066 | 64 / 18 (<0.001) | 66 / 236 | 0.089 / 0.030 / 0.137 |
| HG008T Illumina SNV | d46 - d9 | -0.021 / -0.008 / 0.012 / -0.015 | -0.029 / -0.023 / 0.013 | 0.022 | 25 / 10 (0.017) | 66 / 73 | -0.030 / -0.021 / -0.002 |
| HG008T PacBio SNV | d9 - linear (*) | 0.060 / 0.132 / 0.203 / 0.108 | 0.070 / 0.154 / 0.300 | 0.001 | 31 / 30 (1.000) | 18 / 180 | 0.075 / 0.020 / 0.164 |
| HG008T PacBio SNV | d46 - linear | 0.047 / 0.095 / 0.204 / 0.090 | 0.053 / 0.102 / 0.235 | 0.009 | 40 / 34 (0.561) | 29 / 180 | 0.071 / -0.021 / 0.151 |
| HG008T PacBio SNV | d46 - d9 (*) | -0.013 / -0.037 / 0.001 / -0.018 | -0.017 / -0.051 / -0.064 | 0.007 | 26 / 21 (0.560) | 29 / 18 | -0.004 / -0.041 / -0.013 |
| HG008T ONT SNV | d9 - linear (*) | 0.035 / 0.027 / 0.199 / 0.059 | 0.043 / 0.048 / 0.262 | -0.016 | 23 / 34 (0.185) | 23 / 125 | 0.064 / 0.038 / 0.066 |
| HG008T ONT SNV | d46 - linear | 0.022 / 0.010 / 0.107 / 0.038 | 0.021 / 0.009 / 0.107 | 0.017 | 34 / 22 (0.141) | 52 / 125 | 0.044 / 0.008 / 0.050 |
| HG008T ONT SNV | d46 - d9 (*) | -0.013 / -0.017 / -0.091 / -0.021 | -0.022 / -0.040 / -0.155 | 0.033 | 39 / 16 (0.003) | 52 / 23 | -0.020 / -0.030 / -0.016 |

### 4. chr1 INDEL, W = 50 bp (`tables.md` sections 1 and 7)

At the run's standard threshold (validation recall 0.9). The tables do not give the chr1 curve's best F1, so this is a
stratification, not an arm ranking; the `[valbest]` rows are in `tables.md` section 1.

| set (chr1) | A - B | dF1 0 / 1 / 2+ / all | dR 0 / 1 / 2+ | found A only / B only, all (p); 2+ (p) | tensor share 0 / 1 / 2+, arm A | tensor share 0 / 1 / 2+, arm B |
|---|---|---|---|---|---|---|
| COLO829T Illumina INDEL | d46 - linear | 0.116 / 0.016 / 0.014 / 0.077 | 0.065 / -0.091 / -0.059 | 11 / 8 (0.648); 0 / 1 (1.000) | 0.587 / 0.273 / 0.059 | 0.902 / 0.818 / 0.706 |
| COLO829T fiberseq INDEL | d46 - linear | -0.038 / -0.017 / -0.009 / -0.020 | 0.022 / -0.182 / -0.176 | 9 / 14 (0.405); 0 / 3 (0.250) | 0.685 / 0.318 / 0.118 | 0.946 / 1.000 / 0.941 |
| COLO829T ONT INDEL | d46 - linear | -0.002 / -0.001 / -0.004 / -0.002 | -0.033 / -0.045 / -0.176 | 2 / 9 (0.065); 0 / 3 (0.250) | 0.674 / 0.273 / 0.176 | 0.989 / 0.864 / 0.882 |
| HG008T Illumina INDEL | d9 - linear | 0.068 / 0.067 / 0.035 / 0.108 | -0.015 / -0.048 / -0.059 | 14 / 28 (0.044); 1 / 4 (0.375) | 0.458 / 0.274 / 0.176 | 0.905 / 0.774 / 0.627 |
| HG008T Illumina INDEL | d46 - linear | 0.069 / 0.079 / 0.040 / 0.111 | -0.018 / -0.024 / -0.059 | 15 / 28 (0.066); 1 / 4 (0.375) | 0.500 / 0.310 / 0.176 | 0.905 / 0.774 / 0.627 |
| HG008T Illumina INDEL | d46 - d9 | 0.001 / 0.012 / 0.005 / 0.003 | -0.002 / 0.024 / 0.000 | 6 / 5 (1.000); 0 / 0 | 0.500 / 0.310 / 0.176 | 0.458 / 0.274 / 0.176 |
| HG008T PacBio INDEL | d9 - linear | -0.084 / -0.015 / -0.028 / -0.036 | 0.022 / -0.036 / -0.196 | 26 / 29 (0.788); 0 / 10 (0.002) | 0.562 / 0.333 / 0.216 | 0.982 / 0.929 / 0.863 |
| HG008T PacBio INDEL | d46 - linear | -0.065 / -0.018 / -0.027 / -0.030 | 0.024 / -0.048 / -0.176 | 26 / 28 (0.892); 0 / 9 (0.004) | 0.575 / 0.357 / 0.216 | 0.982 / 0.929 / 0.863 |
| HG008T PacBio INDEL | d46 - d9 | 0.019 / -0.003 / 0.001 / 0.006 | 0.002 / -0.012 / 0.020 | 4 / 3 (1.000); 1 / 0 (1.000) | 0.575 / 0.357 / 0.216 | 0.562 / 0.333 / 0.216 |
| HG008T ONT INDEL | d9 - linear | -0.009 / -0.008 / -0.006 / -0.007 | 0.013 / -0.083 / -0.098 | 15 / 21 (0.405); 0 / 5 (0.062) | 0.584 / 0.357 / 0.196 | 0.887 / 0.738 / 0.725 |
| HG008T ONT INDEL | d46 - linear | -0.010 / -0.007 / -0.004 / -0.008 | 0.020 / -0.071 / -0.078 | 20 / 21 (1.000); 1 / 5 (0.219) | 0.582 / 0.381 / 0.216 | 0.887 / 0.738 / 0.725 |
| HG008T ONT INDEL | d46 - d9 | -0.001 / 0.001 / 0.002 / -0.000 | 0.007 / 0.012 / 0.020 | 11 / 6 (0.332); 1 / 0 (1.000) | 0.582 / 0.381 / 0.216 | 0.584 / 0.357 / 0.196 |

INDEL truth in the BED: HG008T 452 / 84 / 51 (W = 50) and COLO829T 92 / 22 / 17. INDEL F1 all at the standard
threshold: HG008T Illumina 0.133 / 0.241 / 0.244 (linear / d9 / d46), PacBio 0.125 / 0.089 / 0.095, ONT 0.022 / 0.015 /
0.015; COLO829T Illumina 0.032 / 0.109, fiberseq 0.040 / 0.020, ONT 0.004 / 0.002 (linear / d46).

### 5. chr1 COLO829T fiber-seq joint d9 model per purity (`tables.md` sections 2 and 5)

Trained on P100 / P75 / P50 / P25 together on the d9 graph: not a matched arm. SNV, W = 50 bp; the fiber-seq in-sample
rows for comparison.

| call set (chr1) | FP 0 / 1 / 2+ | R 0 / 1 / 2+ | F1 0 / 1 / 2+ | F1 all | F1 0 / 1 / 2+ (W 500) |
|---|---|---|---|---|---|
| COLO829T fiberseq SNV d9 joint, P100 | 190 / 34 / 7 | 0.712 / 0.638 / 0.543 | 0.790 / 0.726 / 0.661 | 0.779 | 0.780 / 0.782 / 0.777 |
| COLO829T fiberseq SNV d9 joint, P75 | 451 / 75 / 19 | 0.652 / 0.585 / 0.486 | 0.698 / 0.632 / 0.553 | 0.687 | 0.696 / 0.693 / 0.675 |
| COLO829T fiberseq SNV d9 joint, P50 | 501 / 91 / 14 | 0.631 / 0.585 / 0.486 | 0.675 / 0.613 / 0.576 | 0.665 | 0.679 / 0.660 / 0.656 |
| COLO829T fiberseq SNV d9 joint, P25 | 354 / 74 / 12 | 0.672 / 0.628 / 0.529 | 0.730 / 0.664 / 0.622 | 0.719 | 0.739 / 0.709 / 0.710 |
| COLO829T fiberseq SNV d46 (in-sample) | 228 / 37 / 3 | 0.715 / 0.684 / 0.557 | 0.784 / 0.754 / 0.696 | 0.779 | 0.760 / 0.794 / 0.784 |
| COLO829T fiberseq SNV linear (in-sample) | 440 / 48 / 33 | 0.637 / 0.564 / 0.471 | 0.690 / 0.650 / 0.485 | 0.679 | 0.641 / 0.714 / 0.692 |

Joint INDEL F1 all: 0.007 / 0.005 / 0.004 / 0.003 (P100 / P75 / P50 / P25; standard threshold).

### 6. chr1 cross-sample, SNV d46, W = 50 bp (`tables.md` section 3)

F1 0 / 1 / 2+ (F1 all). PacBio targets of a COLO829T model use the COLO829T fiber-seq model.

| target set (chr1) | model HG008T (same platform) | model COLO829T (same platform) | in-sample d46 of the target |
|---|---|---|---|
| COLO829T Illumina | 0.722 / 0.646 / 0.547 (0.708) | - | 0.790 / 0.753 / 0.684 (0.783) |
| COLO829T fiberseq | 0.768 / 0.683 / 0.641 (0.754) | - | 0.784 / 0.754 / 0.696 (0.779) |
| COLO829T ONT | 0.645 / 0.581 / 0.440 (0.630) | - | 0.743 / 0.716 / 0.615 (0.736) |
| HG008T Illumina | - | 0.388 / 0.286 / 0.109 (0.359) | 0.608 / 0.506 / 0.205 (0.570) |
| HG008T PacBio | - | 0.580 / 0.486 / 0.361 (0.557) | 0.626 / 0.560 / 0.300 (0.603) |
| HG008T ONT | - | 0.512 / 0.408 / 0.208 (0.487) | 0.607 / 0.513 / 0.233 (0.575) |
| H1437T Illumina | 0.856 / 0.847 / 0.635 (0.849) | 0.457 / 0.468 / 0.368 (0.456) | - |
| H1437T PacBio | 0.834 / 0.837 / 0.720 (0.832) | 0.746 / 0.790 / 0.759 (0.753) | - |
| H1437T ONT | 0.789 / 0.791 / 0.675 (0.787) | 0.759 / 0.787 / 0.685 (0.761) | - |
| H2009T Illumina | 0.942 / 0.920 / 0.833 (0.937) | 0.453 / 0.470 / 0.407 (0.454) | - |
| H2009T PacBio | 0.936 / 0.924 / 0.875 (0.933) | 0.808 / 0.818 / 0.784 (0.809) | - |
| H2009T ONT | 0.891 / 0.847 / 0.776 (0.884) | 0.804 / 0.792 / 0.741 (0.801) | - |
| HCC1395T Illumina | 0.799 / 0.765 / 0.529 (0.790) | 0.266 / 0.284 / 0.308 (0.269) | - |
| HCC1395T PacBio | 0.700 / 0.692 / 0.621 (0.698) | 0.741 / 0.752 / 0.659 (0.741) | - |
| HCC1395T ONT | 0.694 / 0.660 / 0.580 (0.689) | 0.623 / 0.616 / 0.571 (0.621) | - |
| HCC1954T Illumina | 0.399 / 0.332 / 0.189 (0.383) | 0.419 / 0.392 / 0.295 (0.412) | - |
| HCC1954T PacBio | 0.319 / 0.282 / 0.183 (0.311) | 0.479 / 0.425 / 0.361 (0.468) | - |
| HCC1954T ONT | 0.367 / 0.329 / 0.288 (0.360) | 0.432 / 0.408 / 0.333 (0.426) | - |

SNV truth 0 / 1 / 2+ (W = 50): H1437T 3324 / 532 / 85, H2009T 9491 / 1394 / 221, HCC1395T 3032 / 356 / 52, HCC1954T
1212 / 201 / 37. Cross-sample INDEL (Illumina d46, 10 call sets) F1 0 / 1 / 2+ (all):

- HG008T model on COLO829T 0.131 / 0.021 / 0.013 (0.076), H1437T 0.057 / 0.046 / 0.034 (0.055), H2009T
  0.118 / 0.127 / 0.103 (0.118), HCC1395T 0.186 / 0.118 / 0.111 (0.174), HCC1954T 0.042 / 0.040 / 0.017 (0.041).
- COLO829T model on HG008T 0.301 / 0.105 / 0.072 (0.254), H1437T 0.186 / 0.202 / 0.000 (0.181), H2009T
  0.275 / 0.288 / 0.148 (0.274), HCC1395T 0.313 / 0.273 / 0.250 (0.307), HCC1954T 0.121 / 0.119 / 0.000 (0.117).
- The CASTLE INDEL 2+ strata hold only 2-6 truth records.

## Interpretation (all call sets)

1. **HG008T SNV: on all three platforms the graph's F1 gain over GRCh38-only is largest where germline variants
   cluster (W = 50 bp, 2+), and it comes from precision.**
   - dF1 2+: Illumina 0.123 (d9) / 0.135 (d46), PacBio 0.203 (*) / 0.204, ONT 0.199 (*) / 0.107. In stratum 0 the
     gain is 0.022-0.087.
   - The linear arm has 236 / 180 / 125 FPs in the 2+ stratum (Illumina / PacBio / ONT; precision 0.041 / 0.058 /
     0.081), against 18-73 for the graph arms.
   - Recall gains significantly only on Illumina (graph finds 57 / 64 truth SNVs that linear misses, linear 26 / 18;
     McNemar p <0.001). On PacBio (31 / 30, 40 / 34) and ONT (23 / 34, 34 / 22) the two arms find the same truth
     (p 0.14-1.0), so the long-read gain is precision only.
   - At W = 500 bp the gain is again largest in 2+ for Illumina (0.139 / 0.137) and PacBio (0.164 / 0.151), but not
     monotonic (stratum 1 is lowest); ONT is flat (0.064 / 0.038 / 0.066 for d9).
2. **COLO829T SNV: Illumina and ONT do not follow the pattern; this is an artefact of the stratification.**
   - Their linear FPs sit in stratum 0 (2,163 and 2,386 at W = 50; d46: 227 and 384), so the gain there is 0.339 /
     0.334 F1.
   - E4b: 1,928 (Illumina) and 1,454 (ONT) of the linear arm's "other" FPs lie outside the COLO829BL germline truth
     BED, where the truth has no germline records, so they count as density 0 by construction. E10 places the
     Illumina ones in SMaHT extreme regions.
   - **fiber-seq**, whose linear arm has only 440 stratum-0 FPs (202 "other" FPs outside that BED), follows the
     HG008T pattern: dF1 0.094 / 0.104 / 0.211, dP 2+ 0.429 (FP 3 vs 33).
   - Recall gains on all three COLO829T platforms (297-466 truth found by d46 only vs 74-101 by linear only;
     p <0.001), with a positive dR in every stratum.
3. **In germline-dense loci the graph loses truth at candidate generation, much more for INDELs.**
   - SNV, HG008T, 2+ at W = 50: the graph arms build an exact tensor for 43.6-59.0 % of the truth (recall among those
     0.522-0.706); linear builds 69.2-87.2 % but recalls only 0.294-0.407 of them. Recall 2+ ends up similar
     (0.256-0.308 for every HG008T arm).
   - COLO829T SNV 2+: tensor share 0.757-0.786 (d46) vs 0.771-0.843 (linear), a much smaller loss.
   - INDEL: graph tensor share is 0.059-0.216 in 2+ and 0.458-0.685 even in stratum 0, against 0.627-0.941 and
     0.887-0.989 for linear (section 4; E2.4 finds the missing graph tensors' AF in the linear tensors of the same
     reads).
   - So the graph INDEL arms lose recall in strata 1 and 2+: HG008T PacBio 2+ dR -0.196 / -0.176, found by the graph
     only / linear only 0 / 10 (p 0.002) and 0 / 9 (p 0.004); COLO829T fiber-seq / ONT dR 2+ -0.176.
   - Presumably a somatic allele inside or next to a germline bubble is not an exact tensor allele; this was not
     checked per record.
4. **INDEL arm differences at the standard threshold come from precision and are not density-specific.**
   - HG008T Illumina: graph - linear dF1 0.108 / 0.111 overall (dP 0.155 / 0.159), largest in strata 0 and 1 (0.067-
     0.079), smallest in 2+ (0.035 / 0.040).
   - COLO829T Illumina: d46 - linear 0.077 overall (0.109 vs 0.032), mostly in stratum 0 (0.116).
   - HG008T PacBio: graph below linear (-0.036 / -0.030); ONT and COLO829T fiber-seq / ONT: every arm is at F1
     <= 0.040, so the strata say little there.
   - The chr1 curve's best F1 (the INDEL comparison rule) is not computed in E19.
5. **d46 vs d9 (HG008T SNV): lower precision in every W = 500 stratum (dP -0.013 to -0.040), extra recall
   concentrated in germline-dense windows.**
   - W = 50 dF1 all -0.015 / -0.018 (*) / -0.021 (*) (Illumina / PacBio / ONT), with no density trend except ONT 2+
     (-0.091; FP 52 vs 23).
   - W = 500, 2+: dR 0.043 / 0.019 / 0.059; found by d46 only / d9 only 18 / 4 (p 0.004), 13 / 7 (p 0.263), 22 / 3
     (p <0.001). Presumably the d46 graph's extra haplotypes recover truth in germline-dense windows; its extra FPs
     are spread over all strata.
   - INDEL d46 - d9: dF1 within ±0.019 in every stratum, no significant discordance.
6. **Joint fiber-seq d9 model (not a matched arm): the same density gradient at every purity.** F1 0 -> 2+ (W = 50):
   0.790 -> 0.661 (P100), 0.698 -> 0.553 (P75), 0.675 -> 0.576 (P50), 0.730 -> 0.622 (P25). At P100 it equals the
   fiber-seq d46 model overall (0.779 vs 0.779) and is 0.028 / 0.035 lower in strata 1 / 2+. At W = 500 the strata
   are flat (P100 0.780 / 0.782 / 0.777).
7. **Cross-sample: the germline-dense stratum is the hardest for a model from another patient too.**
   - F1 2+ is below F1 0 in 38 of 40 cross-sample call sets (exceptions: COLO829T Illumina model on HCC1395T
     Illumina, 0.266 -> 0.308; COLO829T fiber-seq model on H1437T PacBio, 0.746 -> 0.759).
   - Against the target's in-sample d46 model, the HG008T models lose 0.016-0.098 F1 in stratum 0 and 0.055-0.175 in
     2+ on COLO829T; the COLO829T models lose 0.046-0.220 in stratum 0 and -0.061 to +0.096 in 2+ on HG008T, where
     the in-sample 2+ F1 is already 0.205-0.300.
   - The model's source patient matters more than density: on H1437T / H2009T / HCC1395T Illumina the HG008T model is
     0.399-0.533 F1 above the COLO829T model in stratum 0 (the COLO829T Illumina model's recall is 0.155-0.305 there;
     E12a: transferred threshold).

### Still missing

- **COLO829T d9 single-set runs are not local**, so there is no COLO829T d9 - linear or d46 - d9 pair; the only
  COLO829T d9 call sets are the fiber-seq joint model's (section 5, not a matched arm).
- **HG008T PacBio / ONT SNV d9** are not recipe-matched to d46 (*); the matched reruns are pending.
- **Comparators** (ClairS-TO / DeepSomatic-TO) are not yet re-scored on these chr1 truths; E19 has no comparator rows.
- **Seed repeats** are pending GPU: one training seed per model, so the McNemar p values cover chr1 records only.
- **INDEL at the chr1 curve's best F1**: E19 uses the standard threshold and valbest only.
- Cross-sample INDEL exists only for Illumina d46; no cross-sample call sets on the linear or d9 arms.

## Caveats

- chr1 only, so the 2+ strata are small: HG008T has 39 truth SNVs at W = 50 bp.
- Density counts records, not alleles, and includes non-PASS dipcall records.
- A variant's own germline overlap is excluded, so a germline FP of the same allele does not count itself.
- The McNemar p values are per stratum and not corrected for multiple comparisons.

## Rerun

All four call-level analyses: `bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/call_level_common/collect_all.sh`.
This one only:
`sbatch -J E19_germline_density -o /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E19_germline_density/logs/%x-%j.out /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E19_germline_density/job.sh`
(or `bash collect.sh` on the login node). The run rewrites `tables.md` and `summary.json`; this README's numbers are those of 17:58 (job 406570).
