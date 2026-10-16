# E10: precision and recall per genomic stratum (difficult regions)

Reviewer point: **R1.3** (performance in difficult regions; ONT graph-coordinate calls). Shared code:
`../call_level_common/calls_lib.py`.

**Status (2026-10-15 18:22, job 406569): every in-sample and joint-model call set that exists locally is
processed.** 62 call sets:

- in-sample, 54: 15 SNV and 15 INDEL at the run's main threshold (validation recall 0.9), 15 INDEL valbest and 9 SNV
  `__norescue` re-calls;
  - linear and d46: HG008T Illumina / PacBio / ONT, COLO829T Illumina / fiberseq / ONT;
  - d9: HG008T Illumina / PacBio / ONT;
- joint, 8: the COLO829T fiberseq joint-purity d9 model on P100 / P75 / P50 / P25, SNV and INDEL (trained on four
  purities, so not a matched d9 arm).

What is still missing is listed at the end of the results.

## Question

Where do the calls fail? This analysis reports precision, recall and F1 per stratum, per reference arm and platform:

- low-mappability and segmental-duplication-like regions (SMaHT easy / difficult / extreme);
- outside the HG008T benchmark BED;
- homopolymers and STRs.

It also counts where the unplaced (graph-coordinate) PASS calls lie.

## Data

- **Call sets**: `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1*`, written by
  `machine_learning/pansoma_net_v2/jobs/vcf.sh`. All numbers are on chr1 (the test chromosome), for PASS calls at the
  run's threshold after the GRCh38 PoN (rtg `pon_calls`). Used here: in-sample runs (`main`, `__norescue`,
  `vcf_chr1_valbest`) and the joint-model purity sets (`vcf_chr1_<SET>`). Cross-sample `xs__` runs are not used.
- **Strata.** GIAB stratifications are not local and were not downloaded. Local BEDs:
  - **SMaHT v2 easy / difficult / extreme**:
    `/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_{easy,difficult,extreme}_v2.bed`. COLO829T's somatic BED
    (`pansoma_v2_tensors/COLO829T_truth/SMaHT_v2_easy_difficult_extreme.union.bed`) is only their union. The README
    there says the three files are disjoint and tile chr1-22,X,Y, so the union BED is every base. The three
    separate files are local and are used here.
    - The COLO829T truth VCF's `RGN` annotation agrees with this assignment for all 2,443 SNV truth records in the
      BED (`tables.md` section 7).
    - On chr1 (`tables.md` section 1): easy 179.2 Mb, difficult 25.3 Mb, extreme 44.4 Mb.
    - The HG008T benchmark BED (224.0 Mb of chr1) covers 99.8 % of easy, 98.8 % of difficult and 45.3 % of extreme.
  - **HG008T benchmark BED**, inside vs outside, from the evaluation without the BED (`eval_nobed`).
  - **Repeat context** from GRCh38:
    - INDELs: the pipeline's `RepeatContext.classify` (`indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py`,
      imported, not modified) on the event bases.
    - SNVs: `RepeatContext` classifies INS / DEL only, so SNVs use the same module's `tract()` with the same
      thresholds: a period-1 run of >= 7 bp touching the base is HP>=7; a period 2-6 run of >= 3 units is STR.
      "Touching" includes a run that ends next to the base.
    - Reported as HP>=7 / STR / other.
- **Truth with a tensor**: a truth allele is "with tensor" when `chr1_labels.ndjson` of the call set's own tensor set
  has an exact somatic match for it (partial matches not counted).
- **Unplaced calls**: E9a projection TSVs `/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection/<run>.projected.tsv`.
  - Classes i / ii carry a GRCh38 position; class iii has none.
  - They exist for the nine in-sample graph SNV runs (HG008T Illumina / PacBio / ONT d9 and d46; COLO829T Illumina /
    fiberseq / ONT d46), not for INDEL, `__norescue` or the joint model.

## Method

1. Read each rtg record of `eval_bed` and `eval_nobed` (`tp-baseline`, `fn` = truth; `tp`, `fp` = calls).
2. Normalize the allele: left-align and trim on the set's FASTA.
3. Assign the stratum that contains the first base of its reference interval.
4. Compute recall and precision per stratum:
   - Recall = TP-baseline / truth of the stratum.
   - Precision = TP calls / (TP + FP calls) of the stratum.
   - F1 from the two.
   - `tensor share / R with tensor` splits recall into candidate generation (does the truth allele have a tensor?) and
     the model (recall among the truth with a tensor).
5. Sanity check: the HG008T "inside" row of the no-BED evaluation reproduces rtg's in-BED counts. For example, HG008T
   PacBio d46 has 697 truth and 1,232 calls, against rtg's TP 582 + FN 115 and TP 582 + FP 650.

## Scripts and jobs

| file | role |
|---|---|
| `region_strata.py` | per-record strata -> `DATA/E10/<call set>.v2.tsv.gz`; tables (`tables.md`, `summary.json`) |
| `collect.sh` | runs `region_strata.py` (login node or inside the Slurm wrapper) |
| `../call_level_common/calls_lib.py` | call-set discovery, VCF reading, normalization, BEDs, repeat classes, truth-with-tensor keys |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E10_region_strata/job.sh` | Slurm wrapper (`-p general`, 4 CPUs, 12 GB); logs in `logs/` |

`DATA` = `/scratch/jshen/data/reviewer_experiments_20261013/call_level`. Jobs 405326 (first run, 39 s, 0.87 GB),
405435, 405459 (snapshot 13:56) and 406569 (all call sets, 18:22, 148 s, 0.87 GB max RSS). Per-record caches are
recomputed only when the call set's rtg files are newer.

## Results (all call sets, 2026-10-15 18:00; 62 call sets; chr1; `tables.md` of 18:22)

### 1. SMaHT strata, SNV, in the BED (`tables.md` sections 2 and 3)

| chr1 SNV call set | easy R / P / F1 | difficult R / P / F1 | extreme calls; R / P / F1 | extreme: truth with tensor / R with tensor |
|---|---|---|---|---|
| HG008T Illumina linear | 0.850 / 0.497 / 0.627 | 0.851 / 0.396 / 0.540 | 562; 0.363 / 0.059 / 0.101 | 0.857 / 0.423 |
| HG008T Illumina d9 | 0.927 / 0.496 / 0.647 | 0.851 / 0.454 / 0.592 | 149; 0.264 / 0.161 / 0.200 | 0.626 / 0.421 |
| HG008T Illumina d46 | 0.934 / 0.483 / 0.637 | 0.862 / 0.426 / 0.570 | 206; 0.374 / 0.165 / 0.229 | 0.670 / 0.557 |
| HG008T PacBio linear | 0.888 / 0.503 / 0.642 | 0.851 / 0.435 / 0.576 | 458; 0.451 / 0.090 / 0.149 | 0.857 / 0.526 |
| HG008T PacBio d9 (*) | 0.909 / 0.516 / 0.659 | 0.839 / 0.507 / 0.632 | 102; 0.352 / 0.314 / 0.332 | 0.604 / 0.582 |
| HG008T PacBio d46 | 0.919 / 0.495 / 0.644 | 0.851 / 0.477 / 0.612 | 114; 0.341 / 0.272 / 0.302 | 0.615 / 0.554 |
| HG008T ONT linear | 0.911 / 0.492 / 0.639 | 0.885 / 0.433 / 0.581 | 351; 0.407 / 0.105 / 0.167 | 0.835 / 0.487 |
| HG008T ONT d9 (*) | 0.913 / 0.490 / 0.638 | 0.816 / 0.467 / 0.594 | 116; 0.341 / 0.267 / 0.300 | 0.670 / 0.508 |
| HG008T ONT d46 | 0.933 / 0.482 / 0.636 | 0.851 / 0.435 / 0.576 | 213; 0.451 / 0.192 / 0.270 | 0.659 / 0.683 |
| COLO829T Illumina linear | 0.584 / 0.917 / 0.714 | 0.615 / 0.848 / 0.713 | 2,210; 0.526 / 0.056 / 0.101 | 0.855 / 0.615 |
| COLO829T Illumina d46 | 0.739 / 0.886 / 0.806 | 0.694 / 0.815 / 0.750 | 165; 0.538 / 0.764 / 0.632 | 0.756 / 0.712 |
| COLO829T fiberseq linear | 0.633 / 0.912 / 0.747 | 0.639 / 0.842 / 0.727 | 494; 0.534 / 0.253 / 0.343 | 0.778 / 0.687 |
| COLO829T fiberseq d46 | 0.727 / 0.890 / 0.801 | 0.670 / 0.848 / 0.749 | 199; 0.590 / 0.693 / 0.637 | 0.756 / 0.780 |
| COLO829T fiberseq d9, joint-purity model, P100 (not a matched arm) | 0.721 / 0.897 / 0.800 | 0.684 / 0.847 / 0.757 | 162; 0.534 / 0.772 / 0.631 | 0.752 / 0.710 |
| COLO829T ONT linear | 0.552 / 0.708 / 0.620 | 0.522 / 0.434 / 0.474 | 2,055; 0.491 / 0.056 / 0.100 | 0.799 / 0.615 |
| COLO829T ONT d46 | 0.714 / 0.841 / 0.772 | 0.667 / 0.815 / 0.733 | 281; 0.543 / 0.452 / 0.493 | 0.786 / 0.690 |

(*) d9 checkpoint of a different recipe (1 GPU; E3 caveat 1). Only HG008T Illumina has a recipe-matched d9 / d46 pair.
The `__norescue` re-calls have exactly the placed metrics of their main runs.

Truth counts in the BED:
- HG008T: easy 519, difficult 87, extreme 91.
- COLO829T: easy 1,918, difficult 291, extreme 234.

- **Extreme.** The six linear arms make 351-2,210 PASS calls there, at precision 0.056-0.253 (F1 0.100-0.343). The
  nine graph arms make 102-281, at 0.161-0.764 (F1 0.200-0.637).
- **Easy.** The linear arm's precision is about equal on HG008T (0.492-0.503 vs 0.482-0.516). On COLO829T it is
  higher for Illumina and fiberseq (0.917 / 0.912 vs 0.886 / 0.890) and lower for ONT (0.708 vs 0.841).
  - Its recall is lower everywhere: HG008T 0.850 / 0.888 / 0.911 (Illumina / PacBio / ONT) vs 0.909-0.934 for the
    graph arms; COLO829T 0.552-0.633 vs 0.714-0.739.
  - COLO829T easy: the share of truth with a tensor is the same in both arms (linear 0.822-0.828, d46 0.821-0.827).
    Recall among the truth with a tensor is lower on linear (0.671-0.764 vs 0.870-0.896). So COLO829T's linear recall
    loss sits in the model and threshold, not in candidate generation.
- **Tensor share.** In HG008T easy and difficult, 98.9-100 % of the truth has a tensor in every arm. COLO829T easy has
  only 82.1-82.8 % in every arm (E2.4: COLO829T has many truth alleles below the tensor build AF).
  - HG008T extreme: the linear arms build tensors for 0.835-0.857 of the truth, the graph arms for 0.604-0.670.
  - Recall among the truth with a tensor: linear 0.423-0.526, graph 0.421-0.683.
- **Joint-purity model** (COLO829T fiberseq d9): P100 easy / difficult / extreme F1 0.800 / 0.757 / 0.631, against
  0.801 / 0.749 / 0.637 for the single-set d46 model. P75 / P50 / P25: easy 0.700 / 0.676 / 0.730, extreme 0.571 /
  0.565 / 0.614.

### 2. SMaHT strata, INDEL, in the BED (`tables.md` sections 2 and 3)

The threshold is each run's own: main = validation recall 0.9. E10 does not give the chr1 curve's best F1 per stratum
(the INDEL comparison rule), so these are threshold-level numbers. Truth in the BED: HG008T easy 116, difficult 249,
extreme 222; COLO829T 28 / 39 / 64.

| chr1 INDEL call set, main threshold | easy R / P / F1 | difficult R / P / F1 | extreme R / P / F1 | truth with tensor: difficult / extreme |
|---|---|---|---|---|
| HG008T Illumina linear | 0.612 / 0.247 / 0.352 | 0.189 / 0.086 / 0.118 | 0.162 / 0.040 / 0.064 | 0.884 / 0.770 |
| HG008T Illumina d9 | 0.655 / 0.388 / 0.487 | 0.161 / 0.261 / 0.199 | 0.108 / 0.107 / 0.107 | 0.341 / 0.270 |
| HG008T Illumina d46 | 0.672 / 0.345 / 0.456 | 0.161 / 0.252 / 0.196 | 0.104 / 0.125 / 0.113 | 0.390 / 0.315 |
| HG008T PacBio linear | 0.517 / 0.193 / 0.281 | 0.197 / 0.102 / 0.135 | 0.198 / 0.041 / 0.068 | 0.976 / 0.932 |
| HG008T PacBio d9 | 0.690 / 0.226 / 0.340 | 0.173 / 0.053 / 0.082 | 0.122 / 0.017 / 0.029 | 0.458 / 0.374 |
| HG008T PacBio d46 | 0.690 / 0.196 / 0.305 | 0.173 / 0.057 / 0.085 | 0.126 / 0.019 / 0.034 | 0.482 / 0.383 |
| HG008T ONT linear | 0.638 / 0.049 / 0.090 | 0.137 / 0.007 / 0.014 | 0.144 / 0.006 / 0.011 | 0.851 / 0.779 |
| HG008T ONT d9 | 0.655 / 0.033 / 0.063 | 0.145 / 0.004 / 0.009 | 0.099 / 0.003 / 0.006 | 0.490 / 0.378 |
| HG008T ONT d46 | 0.672 / 0.026 / 0.050 | 0.141 / 0.004 / 0.008 | 0.117 / 0.004 / 0.007 | 0.482 / 0.396 |
| COLO829T Illumina linear | 0.464 / 0.072 / 0.124 | 0.128 / 0.025 / 0.041 | 0.172 / 0.009 / 0.016 | 0.974 / 0.734 |
| COLO829T Illumina d46 | 0.643 / 0.131 / 0.218 | 0.231 / 0.078 / 0.116 | 0.078 / 0.025 / 0.038 | 0.487 / 0.312 |
| COLO829T fiberseq linear | 0.571 / 0.068 / 0.121 | 0.205 / 0.025 / 0.044 | 0.234 / 0.012 / 0.023 | 0.974 / 0.922 |
| COLO829T fiberseq d46 | 0.607 / 0.049 / 0.091 | 0.231 / 0.008 / 0.016 | 0.125 / 0.005 / 0.009 | 0.590 / 0.375 |
| COLO829T fiberseq d9, joint model P100 (not matched) | 0.679 / 0.007 / 0.014 | 0.231 / 0.003 / 0.005 | 0.125 / 0.002 / 0.004 | 0.513 / 0.359 |
| COLO829T ONT linear | 0.679 / 0.005 / 0.011 | 0.179 / 0.001 / 0.002 | 0.172 / 0.001 / 0.002 | 0.949 / 0.938 |
| COLO829T ONT d46 | 0.607 / 0.004 / 0.007 | 0.179 / 0.001 / 0.001 | 0.094 / 0.000 / 0.001 | 0.564 / 0.391 |

- 80 % of the HG008T INDEL truth (471 of 587) lies in difficult or extreme. There, recall is 0.10-0.20 and
  precision at most 0.26 in every arm.
- Easy: the graph arms have the higher F1 for Illumina (HG008T 0.487 / 0.456 vs 0.352; COLO829T 0.218 vs 0.124) and
  HG008T PacBio (0.340 / 0.305 vs 0.281). They have the lower F1 for HG008T ONT (0.063 / 0.050 vs 0.090), COLO829T
  fiberseq (0.091 vs 0.121) and COLO829T ONT (0.007 vs 0.011), where easy-region precision is at most 0.068.
- The graph arms have an exact somatic tensor for 0.34-0.59 of the difficult and 0.27-0.40 of the extreme INDEL truth.
  The linear arms have one for 0.85-0.98 and 0.73-0.94. Only exact label keys count (section 4 and Caveats).
- valbest (`tables.md` section 2) raises easy-region precision to 0.29-0.48 on the HG008T graph arms (F1 0.326-0.420),
  at recall 0.30-0.47.

### 3. HG008T inside vs outside the benchmark BED (`tables.md` section 4)

- Outside the BED, chr1 has 5 truth SNVs and 47 truth INDELs. Precision there is a lower bound, because the truth is
  not complete.
- SNV: the graph arms make 26-52 PASS calls outside (precision 0.077-0.154). The linear arms make 3,246 (Illumina),
  974 (PacBio) and 1,033 (ONT), at precision 0.001-0.003.
- INDEL, main threshold: outside calls are 183 / 225 / 2,184 for linear (Illumina / PacBio / ONT), 60 / 509 / 1,565
  for d9 and 72 / 2,558 / 3,696 for d46 (precision 0.002-0.133).
- Inside the BED the SNV numbers equal the BED evaluation. The INDEL inside truth is 588 against rtg's 587, because a
  record's stratum is that of its first reference base.

### 4. Repeat context, in the BED (`tables.md` section 5)

SNV, `(truth) R / P`:

| chr1 SNV call set | HP>=7 | STR | other |
|---|---|---|---|
| HG008T Illumina linear | (34) 0.235 / 0.095 | (49) 0.531 / 0.158 | (614) 0.837 / 0.370 |
| HG008T Illumina d9 | 0.206 / 0.137 | 0.531 / 0.236 | 0.889 / 0.487 |
| HG008T Illumina d46 | 0.206 / 0.132 | 0.592 / 0.218 | 0.909 / 0.465 |
| HG008T PacBio linear | 0.265 / 0.153 | 0.633 / 0.267 | 0.873 / 0.391 |
| HG008T PacBio d9 (*) | 0.206 / 0.333 | 0.510 / 0.439 | 0.888 / 0.504 |
| HG008T PacBio d46 | 0.235 / 0.258 | 0.571 / 0.384 | 0.889 / 0.484 |
| HG008T ONT linear | 0.235 / 0.093 | 0.592 / 0.254 | 0.896 / 0.426 |
| HG008T ONT d9 (*) | 0.206 / 0.212 | 0.551 / 0.375 | 0.883 / 0.479 |
| HG008T ONT d46 | 0.235 / 0.148 | 0.592 / 0.299 | 0.915 / 0.455 |
| COLO829T Illumina linear | (31) 0.290 / 0.209 | (127) 0.535 / 0.469 | (2,285) 0.589 / 0.389 |
| COLO829T Illumina d46 | 0.419 / 0.542 | 0.583 / 0.796 | 0.726 / 0.874 |
| COLO829T fiberseq linear | 0.290 / 0.474 | 0.535 / 0.667 | 0.634 / 0.752 |
| COLO829T fiberseq d46 | 0.355 / 0.524 | 0.606 / 0.778 | 0.718 / 0.874 |
| COLO829T fiberseq d9, joint model P100 | 0.387 / 0.600 | 0.606 / 0.802 | 0.708 / 0.888 |
| COLO829T ONT linear | 0.226 / 0.200 | 0.583 / 0.385 | 0.544 / 0.339 |
| COLO829T ONT d46 | 0.387 / 0.375 | 0.614 / 0.729 | 0.700 / 0.797 |

Truth with a tensor, HG008T SNV:
- HP>=7: graph arms 0.529-0.676, linear 0.794-0.941;
- STR: graph 0.735-0.776, linear 0.878-0.939;
- other: 0.982-0.995 in every arm.

INDEL, main threshold (truth in the BED: HG008T HP>=7 435, STR 65, other 87; COLO829T 70 / 27 / 34):
- HP>=7 holds 74 % of the HG008T INDEL truth. Recall there is 0.207-0.274 in every arm, precision at most 0.260.
  The graph arms have a tensor for 0.430-0.572 of it, the linear arms for 0.862-0.984.
- STR: HG008T recall 0.046-0.123 on linear, 0.000-0.031 on the graph arms (truth with tensor 0.123-0.200 vs
  0.662-0.908).
- Other (no repeat): the graph arms have the higher F1. HG008T Illumina 0.332 / 0.332 (d9 / d46) vs 0.123, PacBio
  0.123 / 0.127 vs 0.075, ONT 0.083 / 0.082 vs 0.045; COLO829T Illumina 0.142 vs 0.023, fiberseq 0.059 vs 0.038.
  COLO829T ONT is the exception (0.004 vs 0.008).
  Graph recall among the truth with a tensor is 0.818-0.930 on HG008T.

### 5. Unplaced PASS calls per stratum (`tables.md` section 6)

From the E9a projection of the GRCh38 position, now available for all nine in-sample graph SNV runs. Class iii (no
position) is counted separately.

| chr1 SNV call set | unplaced PASS | easy | difficult | extreme | no position | extreme share of the positioned | in the HG008T BED |
|---|---|---|---|---|---|---|---|
| HG008T Illumina d9 | 115 | 29 | 29 | 53 | 4 | 47.7 % | 105 |
| HG008T Illumina d46 | 754 | 72 | 74 | 352 | 256 | 70.7 % | 422 |
| HG008T PacBio d9 (*) | 210 | 40 | 60 | 105 | 5 | 51.2 % | 193 |
| HG008T PacBio d46 | 240 | 46 | 69 | 125 | 0 | 52.1 % | 227 |
| HG008T ONT d9 (*) | 1,410 | 393 | 461 | 431 | 125 | 33.5 % | 1,219 |
| HG008T ONT d46 | 1,258 | 227 | 350 | 450 | 231 | 43.8 % | 981 |
| COLO829T Illumina d46 | 516 | 48 | 55 | 186 | 227 | 64.4 % | - |
| COLO829T fiberseq d46 | 257 | 39 | 58 | 118 | 42 | 54.9 % | - |
| COLO829T ONT d46 | 4,569 | 278 | 355 | 829 | 3,107 | 56.7 % | - |

- Extreme is 17.8 % of chr1 (44.4 of 248.9 Mb) but holds 33.5-70.7 % of the unplaced calls with a position.
- Not projected:
  - the graph INDEL runs: 368-9,899 unplaced PASS at the main threshold, 0-182 at valbest;
  - the `__norescue` re-calls: 6-162;
  - the joint model: 547-705 (SNV) and 4,195-4,800 (INDEL).
- The linear arms have none.

### Still missing

- COLO829T d9 single-set runs are not local. COLO829T has only linear vs d46 (plus the joint-purity d9 model, which is
  not a matched arm).
- Comparators (ClairS-TO / DeepSomatic-TO) are not yet re-scored, so there are no comparator strata.
- Seed repeats are pending GPU. So are recipe-matched d9 reruns for HG008T PacBio / ONT SNV (marked (*)).
- No E9a projection yet for the INDEL, `__norescue` and joint-model unplaced calls.
- No per-stratum chr1 curve-best F1 for INDEL.

## Interpretation (all call sets)

1. **The GRCh38-only arm fails in SMaHT "extreme" regions, on every platform and both tumors.**
   - The six linear SNV call sets make 351-2,210 PASS calls in extreme at precision 0.056-0.253. The graph arms make
     102-281 at 0.161-0.764.
   - Outside the HG008T BED the linear arms make 974-3,246 PASS SNV calls, the graph arms 26-52.
   - In easy regions the linear precision is as high (HG008T, COLO829T Illumina / fiberseq) but its recall is lower.
   - So the graph's precision advantage sits mostly in extreme regions. E4b finds that the linear arms' extra FPs on
     COLO829T and HG008T Illumina are mostly low-AF calls that are not germline or population alleles (on COLO829T
     mostly outside the assembly-confident germline BED): the signature of reads misplaced on a reference that lacks
     the alternative sequence. On HG008T PacBio / ONT the extra FPs are mostly germline.
   - On COLO829T the graph also has the higher recall in every stratum. In easy regions the tensor share is the same
     in both arms, so that part is the model, not candidate generation.
2. **Within the graph arms, SNV difficulty is still concentrated in "extreme".**
   - HG008T graph F1: easy 0.636-0.659, difficult 0.570-0.632, extreme 0.200-0.332.
   - COLO829T d46: easy 0.772-0.806, difficult 0.733-0.750, extreme 0.493-0.637.
   - In HG008T extreme about a third of the truth has no graph tensor (tensor share 0.604-0.670).
   - d9 vs d46, HG008T Illumina (the matched pair): d46 builds more extreme tensors (0.670 vs 0.626), with recall
     0.374 vs 0.264 at the same precision (0.165 vs 0.161). PacBio and ONT pairs are not recipe-matched.
3. **Homopolymers are the weakest stratum, and they dominate INDEL.**
   - SNV HP>=7 recall is 0.206-0.265 on HG008T in every arm. The graph arms build tensors for 53-68 % of these truths,
     linear for 79-94 %.
   - 74 % of the HG008T INDEL truth is in HP>=7, with recall 0.21-0.27 in every arm. STR INDELs are almost never
     recovered by the graph arms (recall at most 0.031).
   - Outside repeats, the graph arms' INDEL F1 is 1.6-2.7 times the linear one on HG008T.
   - The graph arms' low INDEL tensor share in repeats needs E17a: exact keys only are counted here.
4. **ONT.**
   - Placed SNVs: ONT does not lose more than Illumina in extreme regions. HG008T ONT d46 extreme recall is 0.451
     (Illumina d46 0.374, PacBio d46 0.341). COLO829T ONT extreme precision is lower (0.452 vs 0.764 Illumina d46,
     0.693 fiberseq d46).
   - INDEL: at the main threshold ONT precision is at most 0.049 in every stratum and arm.
5. **Unplaced calls are enriched in extreme regions or have no GRCh38 position.**
   - 33.5-70.7 % of the positioned unplaced calls are in extreme (17.8 % of chr1).
   - COLO829T ONT: 3,107 of 4,569 have no consistent anchor pair (E9a class iii); HG008T Illumina d46: 256 of 754.
   - So the "ambiguity" in R1.3 is mostly calls outside the GRCh38 frame. E9c classifies them as population and
     germline alleles.

## Caveats

- chr1 only. Small strata: HG008T extreme has 91 truth SNVs, HP>=7 34; only 5 truth SNVs and 47 truth INDELs lie
  outside the BED; COLO829T has 131 truth INDELs in all (28 in easy).
- SMaHT strata are a proxy for GIAB difficult-region stratifications (not local; not downloaded).
- All three arms exist for HG008T (three platforms); COLO829T has linear and d46 only (its d9 single-set runs are not
  local). HG008T PacBio / ONT SNV d9 checkpoints are of a different recipe (E3 caveat 1).
- INDEL rows are at the run's own threshold (main or valbest). The INDEL comparison rule (the chr1 curve's best F1) is
  not available per stratum.
- "Truth with tensor" counts exact somatic label keys only. An INDEL whose tensor matches only partially (common for
  repeat representations) counts as without a tensor, so the INDEL tensor shares, above all of the graph arms in
  repeats, can understate the candidates.
- Unplaced calls per stratum need an E9a projection of the same call set. E9a has run only for the nine in-sample
  graph SNV runs.
- A record's stratum is that of its first reference base. An indel that spans a stratum boundary counts once.

## Rerun

All four call-level analyses: `bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/call_level_common/collect_all.sh`.
This one only:
`sbatch -J E10_region_strata -o /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E10_region_strata/logs/%x-%j.out /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E10_region_strata/job.sh`
(or `bash collect.sh` on the login node). The run rewrites `tables.md` and `summary.json`; this README's numbers are
those of 18:22 (job 406569).
