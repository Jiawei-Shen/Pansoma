# RIS INDEL experiments with the INDEL pipeline of 9919d4a ("ppr8"), 2026-10-26

The user (2026-10-26): rerun the training and tests that the revised INDEL pipeline changes (SNV is unchanged): the
COLO829T dilutions, HG008T, the seeds, and any other INDEL task. Decisions (the user, same day): keep each experiment's
training design and swap the recipe; also the no-ch6 INDEL runs, the cross-tumor predictions and the post-PoN best-F1
scoring.

## The recipe

`ris/recipes.sh` `KIND=INDEL` = the "INDEL recipe" of `machine_learning/pansoma_net_v2/README.md` at 9919d4a: the 8.7 M
model, `--epochs 12 --batch-size 1024 --lr 1e-4 --weight-decay 0.05 --warmup-epochs 1 --class-weights sqrt
--non-fraction 1.0 --block v2`, partial somatic labels ignored, and the three changes:

- labels that follow the benchmark: `--non-reasons repeat_indel_without_truth --non-reasons-in-bed
  outside_confident_region` (the HG008T / COLO829T labels, rules `5430533aa100`, have no `repeat_indel_without_truth`;
  `outside_confident_region` tensors inside the somatic BED become negatives);
- the five site features: `--site-features $PANSOMA_RUNS/site_features` (`ris/site_features.sh SET`, one file per set);
- `--best-by postpon`: the epoch with the best post-PoN validation somatic F1; PASS (vcf.sh's default, the checkpoint's
  call rule) = the post-PoN validation recall 0.8 of the callable positives. `vcf.sh RUN INDEL postpon _postpon` scores
  the post-PoN best-F1 threshold beside it (`vcf_chr1*_postpon/`).

`--norm dataset` (the default): every model here is new. One H100 per run as the earlier INDEL runs; the purity joint
models on the GPUs, slots and memory of the 7fc7396 joint INDEL runs (Illumina 1, fiberseq / ONT 2 H100s, DDP with
`--batch-size 1024` per GPU).

## Environment changes on RIS

- `bcftools` (the site-feature builder streams the databases with `bcftools query`) is not in the image: bioconda
  bcftools 1.21 (htslib 1.24) in `$PANSOMA_RUNS/env/bcftools` (`conda create -p ... -c conda-forge -c bioconda
  bcftools=1.21` from the login node's miniconda), `$BCFTOOLS` in `ris/env.sh`.
- `ris/site_features.sh SET` (new), `ris/predict_set.sh NAME SET KIND` (new: a run's checkpoint on chr1 of another
  set into `<run>/test_chr1_<SET>/`), `ris/vcf.sh ... [THRESHOLD [SUFFIX]]` (output suffix), `ris/recipes.sh` (the
  INDEL recipe above), `ris/env.sh` (`BCFTOOLS`, `SITE_FEATURES`).
- CPU unit tests of 9919d4a in the container: 106 tests OK, 8 skipped (GPU, real data, pipeline import)
  (`pnv2_unit_tests_9919d4a`, 425051).

## What runs (`submit.sh`; LSF job ids in `$PANSOMA_RUNS/logs/submit_indel_ppr8_20261026.txt`)

| step | runs | what | replaces (earlier runs) | response |
|---|---:|---|---|---|
| S | 27 | site features: HG008T / COLO829T x platform x d9 / d46 / GRCh38-only, COLO829T P75 / P50 / P25 (d9) | - | - |
| A | 18 | seed 0, one model per set: `<S>_<P>_INDEL_{d46,d9,linear}_ppr8` | `*_INDEL_{d46,d9,linear}_normbq`, `*_INDEL_{,d46_}7fc7396`, `*_INDEL_linear_7fc7396` | R1.2, R1.6, R3.1, R3.9, R3.10 |
| B | 6 | no ch6 on d46: `<S>_<P>_INDEL_nopc_d46_ppr8` | `*_INDEL_nopc_d46_normbq` | R1.1, R2.5 |
| C | 3 | COLO829T purity, one model per platform on P100 + P75 + P50 + P25 (d9), chr1 of each: `COLO829_<P>_INDEL_joint_ppr8` | `COLO829_*_INDEL_joint_7fc7396` | R1.5, R2.4, R3.10 |
| D | 36 | seeds 1 and 2 of A: `..._ppr8_s1`, `_s2` | (INDEL had one seed) | R3.9 |
| E | 18 | cross-tumor: each d46 model (seeds 0-2) on chr1 of the other tumor of its platform (Illumina <-> Illumina, PacBio <-> fiberseq, ONT <-> ONT), `<run>/test_chr1_<SET>/` | E12a / E16 §4 single-tumor INDEL rows | R1.7, R3.9 |

Not here: the multi-cancer INDEL models (HCC1395T, HCC1954T, H1437T, H2009T; tensors on wangcluster only, the user runs
them there) and the readnorm (per-read BQ) dilution variant (the BQ read-source fix is still the user's decision; the
Illumina / ONT dilution scores keep that caveat).

## Results

All numbers: `tables.md` (`collect.py`; per-run train / vcfeval logs under `runs/`, every score in `summary.json`).
chr1, rtg vcfeval after the PoN, in the benchmark BED. **PASS** = the checkpoint's call rule (post-PoN validation recall
0.8 of the callable positives); **post-PoN best** = the post-PoN validation best-F1 threshold (`vcf_chr1*_postpon/`);
**best** = the chr1 curve's best point (an upper bound, chosen on chr1). chr1 truth in the BED: HG008T 587 INDELs (418
PoN-tagged), COLO829T 131 (81 PoN-tagged), so the recall ceiling after the PoN is 0.24-0.28 on the graph references and
0.24-0.38 on GRCh38-only (0.28-0.38 on Illumina / PacBio / fiberseq). COLO829T's 131 make its F1 coarse: one call moves it by
about 0.005.

### Per-set models: three seeds (A + D)

PASS F1, mean [min-max] of seeds 0-2 (post-PoN best, mean); earlier = the seed-0 `*_INDEL_<arm>_normbq` PASS
(validation recall 0.9).

| set | d46 | d9 | GRCh38-only | earlier d46 / d9 / GRCh38-only |
|---|---|---|---|---|
| HG008T Illumina | 0.264 [0.261-0.269] (0.229) | 0.260 [0.249-0.267] (0.226) | 0.272 [0.264-0.279] (0.276) | 0.253 / 0.250 / 0.133 |
| HG008T PacBio | 0.276 [0.269-0.283] (0.274) | 0.275 [0.267-0.285] (0.270) | 0.247 [0.244-0.250] (0.251) | 0.089 / 0.110 / 0.133 |
| HG008T ONT | 0.030 [0.018-0.043] (0.154) | 0.069 [0.049-0.083] (0.165) | 0.051 [0.033-0.068] (0.141) | 0.015 / 0.015 / 0.022 |
| COLO829T Illumina | 0.130 [0.123-0.144] (0.078) | 0.163 [0.138-0.209] (0.145) | 0.172 [0.154-0.209] (0.181) | 0.117 / 0.104 / 0.035 |
| COLO829T fiberseq | 0.127 [0.117-0.138] (0.075) | 0.111 [0.071-0.154] (0.119) | 0.178 [0.177-0.179] (0.197) | 0.022 / 0.022 / 0.043 |
| COLO829T ONT | 0.006 [0.005-0.007] (0.078) | 0.010 [0.007-0.017] (0.112) | 0.011 [0.004-0.015] (0.077) | 0.002 / 0.003 / 0.004 |

- Against the earlier runs (seed 0, PASS), every set and reference gains: Illumina / PacBio / fiberseq +0.010 to
  +0.185, ONT +0.004 to +0.068, mostly precision (the earlier recall-0.9 PASS called at precision 0.01-0.09 on PacBio / fiberseq). The
  chr1 best F1 also rises, most on GRCh38-only (HG008T Illumina 0.183 -> 0.283, PacBio 0.177 -> 0.261, COLO829T
  Illumina 0.068 -> 0.184, fiberseq 0.109 -> 0.224).
- Recall at PASS is 0.17-0.29 in every run, below the PoN ceiling; precision decides F1.
- ONT: PASS of these single-set models stays near zero. Their post-PoN recall-0.8 threshold is very low on ONT
  (p >= 0.007-0.14), and PASS calls 1,800-17,800 false INDELs on chr1. The post-PoN best-F1 threshold reaches 0.14-0.17
  (HG008T) and 0.08-0.11 (COLO829T), near the chr1 best (0.15-0.18 / 0.11-0.14).
  Correction (2026-11-02): this is a property of the single-set models, not of the ppr8 call rule as such. The final
  INDEL models of the pipeline (one model per platform on chr2-22 of every tumor, "release", d46, three seeds;
  `../reviewer_response_20261011/experiments/E8_comparators/latest_models.md`) reach ONT PASS F1 0.168 (HG008T,
  P 0.185 / R 0.155) and 0.090 (COLO829T, P 0.064 / R 0.160) with the same rule, against 0.030 / 0.006 here. Their
  PASS is stricter on every set (recall 0.12-0.17 against 0.20-0.23 here, higher precision). Why the single-set ONT
  threshold falls so low (one tumor's training data and validation) is not established.
- Graph vs GRCh38-only, three-seed means: HG008T PacBio is the one set where the graphs win clearly (d46 / d9 +0.029 /
  +0.028 PASS; every graph seed above every GRCh38-only seed). HG008T Illumina: no graph advantage (d46 - GRCh38-only
  -0.008). COLO829T Illumina and fiberseq: GRCh38-only is better (d46 -0.042 / -0.051, d9 -0.009 / -0.068). d46 - d9:
  -0.039 (HG008T ONT) and -0.033 (COLO829T Illumina), otherwise within +-0.016. Part of the GRCh38-only edge is
  its higher recall ceiling after the PoN (the chr1 details table: COLO829T Illumina 0.359 vs 0.267 on d46).
- Seed spread of PASS (max - min): up to 0.018 on HG008T Illumina / PacBio, up to 0.083 on COLO829T (d9).

### No ch6 on d46 (B)

Seed 0 without ch6 against the three seeds with it: HG008T Illumina 0.260 (with: 0.261-0.269), PacBio 0.274
(0.269-0.283), ONT 0.031 (0.018-0.043), COLO829T Illumina 0.150 (0.123-0.144), fiberseq 0.118 (0.117-0.138), ONT
0.012 (0.005-0.007). Differences to seed 0: within +-0.007, except COLO829T fiberseq (-0.020, still inside its with-ch6
seeds): no measurable effect of the path-count channel on INDEL under this recipe.

### COLO829T purity, one model per platform on P100 + P75 + P50 + P25 (C)

PASS F1 per purity P100 / P75 / P50 / P25 (one threshold per model):

| platform | ppr8 PASS | post-PoN best | best | 7fc7396 joint PASS | 7fc7396 joint best |
|---|---|---|---|---|---|
| Illumina | 0.178 / 0.301 / 0.283 / 0.195 | 0.000 / 0.191 / 0.331 / 0.237 | 0.191 / 0.340 / 0.339 / 0.245 | 0.129 / 0.201 / 0.238 / 0.195 | 0.167 / 0.312 / 0.330 / 0.233 |
| fiberseq | 0.054 / 0.040 / 0.029 / 0.019 | 0.104 / 0.104 / 0.082 / 0.099 | 0.172 / 0.153 / 0.124 / 0.116 | 0.007 / 0.005 / 0.004 / 0.003 | 0.069 / 0.049 / 0.055 / 0.074 |
| ONT | 0.003 / 0.003 / 0.003 / 0.002 | 0.126 / 0.106 / 0.066 / 0.021 | 0.135 / 0.136 / 0.086 / 0.039 | 0.001 / 0.001 / 0.001 / 0.001 | 0.055 / 0.051 / 0.070 / 0.062 |

- Illumina: PASS (p >= 0.41) above the 7fc7396 joint at P100-P50, equal at P25; at P100 0.178, inside the single-set
  d9 seeds (0.138-0.209). Its post-PoN best-F1 threshold is p >= 0.919 (the best epoch's pooled validation), which
  calls 2-55 sites per purity on chr1 (P100: 2 calls, none true).
- fiberseq: PASS at p >= 0.04 calls 1,083-2,432 false INDELs per purity (precision 0.01-0.03), the ONT pattern; the
  ranking is better than 7fc7396's (best 0.12-0.17 vs 0.05-0.07). The best epoch is 4.
- ONT: PASS at p >= 0.018 calls 20,361-21,461 false INDELs per purity (F1 0.002-0.003). The post-PoN best-F1
  threshold (p >= 0.516) gives 0.126 / 0.106 / 0.066 / 0.021, and the ranking beats 7fc7396's at P100 / P75 (best
  0.135 / 0.136 vs 0.055 / 0.051) but not at P50 / P25 (0.086 / 0.039 vs 0.070 / 0.062). The best epoch is 2 (post-PoN
  validation F1 0.106-0.165 over the 12 epochs).
- The Illumina and ONT dilutions still carry the read-source BQ cue (the readnorm fix is the user's decision): their
  dilution scores are likely optimistic.

### Cross-tumor (E): each tumor's d46 models on chr1 of the other tumor

| model trained on | tested on | PASS seeds 0 / 1 / 2 | target's own d46 models PASS (mean) | cross - own |
|---|---|---|---|---|
| HG008T Illumina | COLO829T Illumina | 0.139 / 0.141 / 0.145 | 0.130 | +0.011 |
| COLO829T Illumina | HG008T Illumina | 0.264 / 0.250 / 0.239 | 0.264 | -0.013 |
| HG008T PacBio | COLO829T fiberseq | 0.140 / 0.132 / 0.134 | 0.127 | +0.008 |
| COLO829T fiberseq | HG008T PacBio | 0.241 / 0.160 / 0.242 | 0.276 | -0.062 |
| HG008T ONT | COLO829T ONT | 0.002 / 0.002 / 0.002 | 0.006 | -0.004 |
| COLO829T ONT | HG008T ONT | 0.037 / 0.024 / 0.032 | 0.030 | +0.001 |

A model trained on the other tumor of the platform is as good as the target's own models on Illumina and HG008T PacBio
-> COLO829T fiberseq (within -0.013 to +0.011). COLO829T fiberseq -> HG008T PacBio is -0.062, mostly seed 1 (0.160;
seeds 0 and 2 0.241 / 0.242, -0.035). The cross predictions' post-PoN best-F1 threshold is unstable (COLO829T models on
HG008T: 0.055-0.095 Illumina, 0.000-0.169 PacBio): the threshold comes from the source tumor's validation.

### PASS vs the post-PoN best-F1 threshold

Per-set models on Illumina, PacBio and fiberseq (three-seed means): post-PoN best - PASS is -0.052 to +0.019; it
loses most on the graph references (HG008T Illumina d46 / d9 -0.035 / -0.034, COLO829T Illumina d46 -0.052,
fiberseq d46 -0.052) and gains a little on GRCh38-only (+0.004 to +0.019). The post-PoN best-F1 threshold rests on few
validation positives (somatic validation tensors without a PoN tag: HG008T 94-107, COLO829T 36-39 per d46 set) and is
noisy. On ONT, and on the fiberseq / ONT purity models, the post-PoN best is far better (PASS precision 0.001-0.05);
for ONT this holds for these single-set models, not for the multi-tumor release models (see the ONT correction above).

### Throughput

The 2-GPU DDP joint runs trained at 6,950-7,450 tensors/s while single-GPU runs shared the host, and the ONT joint at
8,300-8,560 once it ran alone (epochs 6-12), against 16,700-17,000 for the 7fc7396 joint runs on the same GPUs / slots /
memory (GPU utilization about 28%; CPU 3-6 of the 14 slots). Single-GPU runs (the 1-GPU Illumina joint too) are
unchanged (8,400-9,000 tensors/s). Not diagnosed (no pipeline change on RIS); the fiberseq joint took 3.3 h (7fc7396:
1.5 h), the ONT joint 5.6 h (2.5 h).

## Follow-ups (2026-10-29)

- `comparators/`: ClairS-TO (the user's runs) and DeepSomatic tumor-only (P100: the user's runs; P75 / P50 / P25: run
  here on the mixture BAMs, exec-419 H100s) on COLO829T chr1 at every purity, scored with E8's protocol (the P100
  scores reproduce E8's exactly), next to the ppr8 joint models and the P100 models. SNV: Pansoma at or above both
  callers at every purity (fiberseq P75 / P50 within 0.02). INDEL: ClairS-TO / DeepSomatic ahead on long reads at
  Pansoma's PASS; at the post-PoN best-F1 threshold Pansoma is close to ClairS-TO on ONT; the Illumina lead comes from
  the joint model (the cue-free P100 model is below DeepSomatic at P75-P25). `comparators/README.md`, `tables.md`.
- `p100_on_dilutions/`: the P100 single-set d9 models (INDEL ppr8, SNV normbq; seeds 0-2) on chr1 of the mixtures,
  `<run>/test_chr1_<SET>/`, `vcf_chr1_<SET>[_postpon]/`; their scores are the "P100 model" columns of
  `comparators/tables.md` (no read-source cue: these models never saw a mixture).
- `fp_context/`: why GRCh38-only beats d46 on COLO829T INDELs. d46's extra PASS false positives are mostly HP>=7
  (66% of d46's FPs vs 28-37% of GRCh38-only's); they are (a) homopolymer candidates both arms have, called by d46's
  lower PASS threshold, and (b) germline homopolymer-length alleles carried on a non-reference path: the graph
  candidate is the read's edit relative to that path, written on GRCh38 as a rare +-1-3 change that the PoN's exact
  match misses (67-81% of these have another, common allele at the tract), while GRCh38-only's candidate there is the
  common allele the PoN removes. Separately, 10-31% of all PASS FPs (both arms) are alleles at AF >= the PoN floor in
  gnomAD / CoLoRSdb after normalisation that the PoN's exact POS / REF / ALT match misses. `fp_context/tables.md`,
  `popaf.md`.
