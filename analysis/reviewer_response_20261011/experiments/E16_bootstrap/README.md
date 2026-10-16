# E16: block-bootstrap confidence intervals of the chr1 results (R3.9)

Reviewer point: **R3.9** (stability and statistical uncertainty: every number is one chr1 test without an interval).
Also used by R3.1 / R1.1 / R1.2 (are the reference-arm differences larger than the sampling noise?) and R1.7
(cross-sample vs in-sample).

## Question

1. How wide is the sampling uncertainty of each chr1 P / R / F1 after the PoN, inside the benchmark BED and without it?
2. Are the differences between reference arms (GRCh38-only vs HPRC v1.1 d9 vs HPRC v2.1 d46) of the same sample /
   platform / kind larger than that uncertainty, when both arms are resampled on the same blocks?
3. The same for cross-sample vs in-sample models on the same target set, and for the other paired variants (--norm
   dataset vs not, valbest vs recall-0.9 threshold, no rescue vs rescue, comparators vs Pansoma).

## Method

- **Call sets.** Every finished `vcf_chr1*/` under `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/`
  (`callsets.py`; finished = rtg's `done` marker and report.json):
  - in-sample runs of the three arms; `__norescue`; `vcf_chr1_valbest`;
  - cross-sample `xs__...` (and `__normds`);
  - the COLO829T fiber-seq joint model (`vcf_chr1_<SET>` per purity);
  - comparators listed in `comparators.tsv` (slot; empty).
- **Records.** rtg's `eval_bed/vcfeval/pon_calls/` and `eval_nobed/vcfeval/pon_calls/` from `jobs/vcf.sh`. These are
  the PASS calls after the PoN: tp-baseline and fn (truth side), tp and fp (call side).
  - P = TPc / (TPc + FP), R = TPb / (TPb + FN), F1 = 2PR / (P + R), as rtg defines them.
  - The totals reproduce report.json `pon.at_threshold` for every call set and mode (204 of 204 call-set / mode pairs).
- **Blocks.** chr1 is cut into fixed 1-Mb windows (main analysis) and 5-Mb windows (sensitivity).
  - In the BED, the block set is every window that overlaps the set's benchmark BED: 231 (HG008T), 249 (COLO829T) and
    232 (HCC1395T, HCC1954T, H1437T, H2009T) 1-Mb blocks.
  - Without the BED, it is every window with a non-N GRCh38 base: 232.
  - Each record goes to the block of its POS. No record fell outside its block set.
- **Bootstrap.** 2,000 resamples of n blocks with replacement (n = size of the block set), 95 % percentile intervals.
  - The resample multiplicities are drawn once per block set (seed 20261015), so every call set on the same truth and
    BED uses the same resamples.
  - A paired difference A − B is computed resample by resample. The pair is skipped if the truth or BED differ.
  - `P(dF1 > 0)` = the share of resamples in which A has the higher F1.
- **Pairs.**
  - arms: d9 − linear, d46 − linear, d46 − d9 (in-sample, standard threshold). A pair with a d9 checkpoint of a known
    different recipe (E3 caveat 1: HG008T PacBio / ONT SNV d9, 1 GPU) is marked (*);
  - joint d9: the COLO829T fiber-seq joint model against d46 / linear. This is flagged, because the joint model is
    trained on four purities and is not a matched arm;
  - cross-sample − in-sample on the same target set;
  - `__normds` − old cross-sample;
  - valbest − standard, norescue − standard (same predictions);
  - Pansoma − comparator.

## What the intervals do not cover

- **Unplaced PASS calls** (off-reference, no GRCh38 position; linear_vcf's `*.linear.unplaced.vcf.gz`) are outside
  rtg, so they are outside every interval here. They are listed per call set (`tables.md` §1, last column). Their
  effect is measured by E3 §2 / E9a / E9c and by the `__norescue` re-calls (E3b).
- **One training seed per model.** The intervals are sampling uncertainty over chr1 regions for a fixed model. Seed
  variance needs repeated GPU training (the multi-seed item of the RIS list,
  `overview_20261013.md` §3.1).
- **One test chromosome.** The blocks are chr1 only. Blocks are assumed independent beyond 1 Mb; the 5-Mb blocks
  check this.
- **Fixed thresholds.** The PASS threshold is taken as given (validation recall 0.9, or valbest). Its own uncertainty
  is not resampled.

## Files

| file | content |
|---|---|
| `callsets.py` | discovery of call sets, shared with E17a and E18 (run-name parsing, rtg directories, truth / BED, comparators) |
| `comparators.tsv` | comparator slot (header and a commented example; no active rows) |
| `e16_bootstrap.py` | the analysis (login node; about 30 s now, 0.9 GB) |
| `collect.sh` | reruns it on whatever exists (`--slurm` submits it as one small job instead) |
| `tables.md` | §1a / §1b per call set (in the BED / no BED); paired sections; block-size sensitivity; checks |
| `ci.tsv` | every call set x mode x block size: counts, P / R / F1 with CI, width, SE |
| `paired.tsv` | every pair x mode x block size: dF1, dP, dR with CI, SE, P(d > 0) |
| `summary.json` | all of the above, block-set sizes, report checks |

## Results (all call sets, 2026-10-15 18:03; 102 call sets)

`tables.md` was generated after the last call set finished (17:50). It has 102 call sets:

- **In-sample SNV (24):** linear for HG008T Illumina / PacBio / ONT and COLO829T Illumina / ONT / fiber-seq; d46 for
  the same six; d9 for HG008T Illumina / PacBio / ONT only; 9 `__norescue` copies (the d9 and d46 runs).
- **In-sample INDEL (30):** the same 15 runs (linear 6, d46 6, d9 3), each at the standard threshold and at
  `vcf_chr1_valbest`.
- **COLO829T fiber-seq joint d9 model (8):** SNV and INDEL on P100 / P75 / P50 / P25. It is trained on four purities,
  so it is not a matched arm.
- **Cross-sample (40):**
  - SNV (30): HG008T <-> COLO829T on all three platforms, and both tumors' d46 models -> HCC1395T / HCC1954T /
    H1437T / H2009T on all three platforms;
  - INDEL (10): the same on Illumina only.
- **Checks.** rtg totals = report.json in 204 of 204 call-set / mode pairs. 0 records fall outside the block sets.

**1. Interval widths** (`tables.md` §1a, §7; in the BED, 1-Mb blocks):

| chr1 call set (in the BED) | linear F1 [95 % CI] | d9 F1 [95 % CI] | d46 F1 [95 % CI] |
|---|---|---|---|
| HG008T Illumina SNV | 0.470 [0.409, 0.531] | 0.585 [0.550, 0.621] | 0.570 [0.535, 0.604] |
| HG008T PacBio SNV | 0.514 [0.426, 0.586] | 0.621 [0.587, 0.653] (*) | 0.603 [0.569, 0.635] |
| HG008T ONT SNV | 0.537 [0.492, 0.578] | 0.596 [0.559, 0.630] (*) | 0.575 [0.536, 0.610] |
| COLO829T Illumina SNV | 0.468 [0.327, 0.668] | not local | 0.783 [0.766, 0.799] |
| COLO829T ONT SNV | 0.418 [0.317, 0.532] | not local | 0.736 [0.704, 0.762] |
| COLO829T fiber-seq SNV | 0.679 [0.630, 0.725] | joint model: 0.779 [0.763, 0.796] | 0.779 [0.758, 0.797] |

(*) d9 checkpoint trained on 1 GPU (global batch 1,024), linear / d46 on 4 GPUs: not a matched arm (E3 caveat 1).

- The d9 / d46 in-sample SNV F1 intervals are 0.033-0.074 wide. HG008T's are the widest because it has only 697
  truth SNVs on chr1 in the BED.
- The linear SNV intervals are 0.086-0.341 wide. Their precision intervals are very wide, for example COLO829T
  Illumina P = 0.391 [0.226, 0.777]. A few blocks hold most of their false positives (point 2).
- Other F1 interval widths (§7, in the BED):
  - joint-model SNV: 0.034-0.039;
  - cross-sample SNV: 0.009-0.083;
  - in-sample INDEL: 0.002-0.114;
  - cross-sample INDEL: 0.020-0.180.
- **5-Mb blocks** widen the intervals only slightly. In the BED, the median 5-Mb / 1-Mb width ratio is 1.16 for
  in-sample SNV (range 1.10-1.54), 1.05 for in-sample INDEL, 1.22 for cross-sample SNV and 1.12 for cross-sample
  INDEL.

**2. Reference arms, SNV: both pangenome arms beat GRCh38-only on every set** (§2; standard threshold; in the BED,
1-Mb blocks; P(dF1 > 0) >= 0.998 in every row; every 5-Mb interval also excludes 0):

| chr1 SNV pair (in the BED) | dF1 [95 % CI] | dP [95 % CI] | dR [95 % CI] |
|---|---|---|---|
| COLO829T Illumina d46 − linear | +0.316 [+0.116, +0.453] | +0.476 [+0.086, +0.634] | +0.132 [+0.109, +0.153] |
| COLO829T ONT d46 − linear | +0.319 [+0.208, +0.417] | +0.447 [+0.276, +0.569] | +0.149 [+0.122, +0.175] |
| COLO829T fiber-seq d46 − linear | +0.099 [+0.062, +0.141] | +0.120 [+0.028, +0.211] | +0.083 [+0.064, +0.100] |
| HG008T Illumina d9 − linear | +0.116 [+0.062, +0.171] | +0.117 [+0.061, +0.172] | +0.044 [+0.016, +0.075] |
| HG008T Illumina d46 − linear | +0.101 [+0.050, +0.152] | +0.094 [+0.040, +0.144] | +0.066 [+0.039, +0.095] |
| HG008T PacBio d9 (*) − linear | +0.108 [+0.046, +0.187] | +0.125 [+0.057, +0.205] | +0.001 [−0.024, +0.028] |
| HG008T PacBio d46 − linear | +0.090 [+0.028, +0.172] | +0.100 [+0.032, +0.182] | +0.009 [−0.019, +0.037] |
| HG008T ONT d9 (*) − linear | +0.059 [+0.034, +0.091] | +0.072 [+0.044, +0.107] | −0.016 [−0.039, +0.007] |
| HG008T ONT d46 − linear | +0.038 [+0.010, +0.069] | +0.038 [+0.009, +0.070] | +0.017 [−0.004, +0.040] |

- The gain is mostly precision.
  - Recall also rises on COLO829T (+0.083 to +0.149) and on HG008T Illumina (+0.044 / +0.066).
  - It does not rise on HG008T PacBio / ONT, where the recall intervals cover 0.
- **Without the BED** the HG008T gaps are larger: d9 / d46 − linear = +0.380 / +0.363 on Illumina, +0.221 / +0.198 on
  ONT and +0.255 / +0.238 on PacBio. The reason is false positives outside the HG008T BED: the linear arm has 4,332 /
  1,935 / 1,940 FP on the whole of chr1, against 1,089 / 904 / 969 inside the BED. The COLO829T numbers do not change,
  because its BED already covers these regions.
- **Where the linear false positives are.** This is a login-node check of the rtg fp / tp / fn files; it is not in
  `tables.md`.
  - Many linear FPs fall in chr1:121-126 Mb, the centromere and its flanks. In the BED this holds 1,888 of 2,221 FP
    (COLO829T Illumina), 1,346 of 2,575 (ONT), 254 of 521 (fiber-seq), 227 of 1,089 (HG008T Illumina), 209 of 969
    (PacBio) and 80 of 904 (ONT).
  - Without the BED, the HG008T counts are 3,119 of 4,332 (Illumina), 1,067 of 1,940 (PacBio) and 914 of 1,935 (ONT).
  - The d46 runs have at most 17 FP there; the d9 runs have 0.
  - The pangenome advantage stays when 121-126 Mb is left out (point F1, in the BED, linear vs d46): COLO829T
    Illumina 0.678 vs 0.787, ONT 0.530 vs 0.739, fiber-seq 0.721 vs 0.782; HG008T Illumina 0.520 vs 0.575 (d9
    0.587), PacBio 0.567 vs 0.606 (d9 0.623), ONT 0.557 vs 0.577 (d9 0.598).
  - So the d46 − linear gap shrinks without 121-126 Mb, but stays positive: +0.316 -> +0.109 (COLO829T Illumina),
    +0.319 -> +0.209 (ONT), +0.099 -> +0.061 (fiber-seq), +0.101 -> +0.055 (HG008T Illumina), +0.090 -> +0.039
    (PacBio), +0.038 -> +0.020 (ONT). These are point values without intervals. E10 stratifies the same calls by
    region.
- **d46 − d9 (HG008T only):**
  - Illumina, recipe-matched (`HG008_Illumina_SNV_7fc7396_ddp4_b1024`): −0.015 [−0.030, −0.001], P(dF1 > 0) = 0.021.
  - PacBio (*): −0.018 [−0.030, −0.005].
  - ONT (*): −0.021 [−0.039, −0.004].
  - In all three, d46 has the higher recall (+0.022 / +0.007 / +0.033) and the lower precision (−0.023 / −0.025 /
    −0.034).
  - The intervals barely exclude 0 and cover chr1 sampling with one seed per arm, so they support "d9 and d46 are
    close" (E3: −0.021 to +0.031) rather than a ranking.
- **COLO829T has no single-set d9 run.** Only the joint d9 model exists for it (§3, flagged):
  - joint d9 − linear (fiber-seq) = +0.100 [+0.061, +0.141];
  - d46 − joint d9 = −0.001 [−0.014, +0.010], P(dF1 > 0) = 0.48, a tie.

**3. Reference arms, INDEL** (§2, §5). This is the PASS set at the standard threshold (validation recall 0.9) after the
PoN. The benchmark compares INDEL by the chr1 curve's best F1 instead, which is a threshold sweep and is not resampled
here.

- **Illumina: the pangenome is better.**
  - HG008T d9 − linear = +0.108 [+0.088, +0.129]; d46 − linear = +0.111 [+0.089, +0.133].
  - COLO829T d46 − linear = +0.077 [+0.048, +0.108].
  - HG008T d46 − d9 = +0.003 [−0.009, +0.015], a tie.
- **Long reads: the linear arm is better at this threshold**, and every interval excludes 0. The difference is
  precision. All arms flood at the standard threshold, for example 18,416 FP for HG008T ONT d46 and 28,301 for
  COLO829T ONT d46 in the BED.
  - HG008T PacBio: d9 − linear −0.036 [−0.050, −0.023]; d46 − linear −0.030 [−0.041, −0.020].
  - HG008T ONT: −0.007 [−0.009, −0.005] and −0.008 [−0.010, −0.006].
  - COLO829T fiber-seq d46 − linear: −0.020 [−0.030, −0.011].
  - COLO829T ONT: −0.002 [−0.003, −0.001].
  - The recall intervals of these pairs all cover or touch 0.
- **The valbest threshold** (§5) raises long-read INDEL F1 a lot:
  - HG008T ONT d46: +0.136 [+0.108, +0.164] (to 0.151);
  - HG008T PacBio d9 / d46: +0.137 / +0.104 (to 0.225 / 0.198);
  - COLO829T fiber-seq d46: +0.070 (to 0.090);
  - COLO829T ONT d46: +0.025 (to 0.027).
  - At valbest the pangenome arms are above the linear arm, whose valbest F1 is 0.036 (HG008T PacBio), 0.047
    (HG008T ONT) and undefined (0 TP) for COLO829T fiber-seq / ONT. No arm pair is resampled at valbest.
  - On Illumina, valbest lowers F1: HG008T d46 −0.039 [−0.073, −0.008], linear −0.053, COLO829T d46 −0.065; HG008T
    d9 −0.025 [−0.054, +0.002].

**4. Cross-sample** (§3, §4; in the BED):

| chr1 pair (in the BED): model -> target vs the target's own d46 model | dF1 [95 % CI] | dP | dR |
|---|---|---|---|
| HG008T Illumina -> COLO829T Illumina SNV | −0.075 [−0.088, −0.063] | −0.075 | −0.074 |
| HG008T ONT -> COLO829T ONT SNV | −0.106 [−0.136, −0.072] | −0.109 | −0.102 |
| HG008T PacBio -> COLO829T fiber-seq SNV | −0.024 [−0.038, −0.009] | −0.049 | −0.006 |
| COLO829T Illumina -> HG008T Illumina SNV | −0.211 [−0.252, −0.172] | +0.175 | −0.597 |
| COLO829T ONT -> HG008T ONT SNV | −0.088 [−0.125, −0.054] | +0.147 | −0.439 |
| COLO829T fiber-seq -> HG008T PacBio SNV | −0.046 [−0.076, −0.017] | +0.094 | −0.287 |
| HG008T Illumina -> COLO829T Illumina INDEL | −0.033 [−0.049, −0.019] | −0.026 | +0.015 |
| COLO829T Illumina -> HG008T Illumina INDEL | +0.010 [−0.008, +0.028] | +0.078 | −0.032 |

- All six SNV transfers lose F1 beyond the noise. The two INDEL transfers lose less: HG008T -> COLO829T loses, and
  COLO829T -> HG008T is a tie (P(dF1 > 0) = 0.87).
- The COLO829T SNV models gain precision on HG008T but lose much more recall. That is the transferred-threshold effect
  (E12a); E18 adds a context part (a UV prior).
- On the four CASTLE targets there is no in-sample model to pair with, so these are separate F1 values with their
  intervals (§1a):

| chr1 SNV F1 (in the BED), HG008T model / COLO829T model | Illumina | ONT | PacBio (COLO829T: fiber-seq model) |
|---|---|---|---|
| H1437T | 0.849 / 0.456 | 0.787 / 0.761 | 0.832 / 0.753 |
| H2009T | 0.937 / 0.454 | 0.884 / 0.801 | 0.933 / 0.809 |
| HCC1395T | 0.790 / 0.269 | 0.689 / 0.621 | 0.698 / 0.741 |
| HCC1954T | 0.383 / 0.412 | 0.360 / 0.426 | 0.311 / 0.468 |

  - The HG008T models are ahead on H1437T and H2009T (all platforms) and on HCC1395T Illumina / ONT. The COLO829T
    models are ahead on HCC1954T (all platforms) and on HCC1395T PacBio.
  - The COLO829T Illumina model keeps a high precision (0.908-0.973) but has a recall of 0.157-0.305 there.
  - HCC1954T is low for every model, because precision is 0.222-0.373.
  - Cross-sample INDEL (Illumina; F1 in the BED, HG008T / COLO829T model): H1437T 0.055 / 0.181, H2009T 0.118 /
    0.274, HCC1395T 0.174 / 0.307, HCC1954T 0.041 / 0.117. On H1437T / H2009T / HCC1954T the HG008T INDEL model has a
    recall of 0.56-0.81 at a precision of 0.021-0.064.

**5. COLO829T fiber-seq joint d9 model by purity** (§1a; in the BED):

- SNV F1: P100 0.779 [0.763, 0.796], P75 0.687 [0.666, 0.706], P50 0.665 [0.648, 0.682], P25 0.719 [0.702, 0.737].
- INDEL F1 at the standard threshold: 0.003-0.007, with 9,879-19,537 FP.

**6. No rescue (E3b re-calls).** All 9 `__norescue` call sets are identical to the standard ones inside rtg (0
difference in every resample; §6). Removing the off-reference rescue changes only the unplaced PASS calls:

| chr1 call set | unplaced PASS, standard | unplaced PASS, no rescue |
|---|---:|---:|
| HG008T Illumina d9 | 115 | 12 |
| HG008T Illumina d46 | 754 | 29 |
| HG008T PacBio d9 | 210 | 9 |
| HG008T PacBio d46 | 240 | 6 |
| HG008T ONT d9 | 1,410 | 26 |
| HG008T ONT d46 | 1,258 | 162 |
| COLO829T Illumina d46 | 516 | 18 |
| COLO829T ONT d46 | 4,569 | 98 |
| COLO829T fiber-seq d46 | 257 | 78 |

### Still missing

- **COLO829T d9 single-set runs are not local.** COLO829T has no matched d9 − linear or d46 − d9 pair; only the joint
  d9 model (flagged) stands in for it.
- **HG008T PacBio / ONT SNV d9 are the 1-GPU recipe** (*). Matched d9 reruns are needed for a clean d46 − d9 on those
  two platforms.
- **Comparators are not yet re-scored.** `comparators.tsv` has no active rows, so there are no Pansoma − comparator
  pairs.
- **No `__normds` cross-sample re-calls exist yet**, so there are no `__normds` − old cross-sample pairs.
- **Seed repeats are pending GPU.** Every interval is chr1 sampling noise for one trained model.
- **Unplaced PASS calls stay outside every interval** (counts in §1, last column; E3 §2 / E9a / E9c).

## Rerun when more runs finish

```bash
bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E16_bootstrap/collect.sh
```

All three call-level analyses (E16, E17a, E18) at once:

```bash
cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments && \
  for e in E16_bootstrap E17a_indel_strata E18_sbs96; do bash $e/collect.sh; done
```

`bash collect.sh --slurm` submits the same thing as one small job (`-p general -c 2`) with its log in
`tmp/reviewer_cpu_20261013/<ID>/logs/`. New in-sample arms, valbest, normds and comparator rows (`comparators.tsv`)
are paired automatically.
