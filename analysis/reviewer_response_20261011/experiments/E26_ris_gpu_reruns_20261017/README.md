# E26: RIS GPU reruns of 2026-10-17 (no-ch6 retraining, second seeds, matched d9 recipe, per-dataset normalisation, per-read BQ centring, reproduction)

> **Superseded (2026-10-19).** All 87 RIS runs are done: three seeds per arm and E13 (538012c), and the `--norm dataset`
> fix (f344bca; BQ only per dataset, statistics of all chromosomes; normbq retraining running on RIS). Use the RIS
> summary `analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/` for every number. This directory keeps
> the 10-17 two-seed collection (chr1 AP per run, seed-matched arm table) for the record only; its seed statements are
> two-seed and its normds rows predate the fix.

Reviewer points: **R1.1** (retraining without channel 6), **R3.9** (seeds, reproducibility), **R3.1** / **R1.2** (matched
v1.1 recipe, arms by seed), **R3.m1** (normalisation setting), **R2.4** / **R1.7** / **R2.1** (read-source BQ cue of the
dilutions). The runs were trained and called on RIS by the user and synced on 2026-10-17 into
`/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/` (89 run folders; the list is
`tmp/reviewer_cpu_20261017/ris_check/runs_20261017.txt`, selected by ctime because `rsync -a` keeps the RIS mtimes).

**The canonical summary of these runs is the RIS one** (pushed from RIS as 608fc69 on 2026-10-17):
`analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/` (README, `tables.md`, `collect.py`; group
definitions, run names and LSF ids). It plans three seeds per arm (`_s1`, `_s2`; seed 2 and the E13 learning curve
were still running on 2026-10-17) and recommends three-seed means as the final numbers. E26 was written before that
summary was pulled; its numbers agree with it and it adds the chr1 AP of each run and the seed-matched arm table (3b).

`collect.py` (this directory) writes `tables.md` and `summary.json`. All numbers: chr1, in the benchmark BED, PASS at
the validation recall-0.9 threshold, after the PoN (rtg vcfeval of `jobs/vcf.sh` on RIS); "best" = the chr1 curve's
test-tuned best F1 (an upper bound); "AP" = chr1 per-tensor somatic average precision before the PoN
(`test_chr1/*.metrics.json`). INDEL numbers are the `vcf_chr1*_pon71cb51a` reports (current INDEL PoN rule, HP >= 8 bp
AF >= 0.01, else >= 1e-4); the SNV baselines not synced today come from the 7fc7396 benchmark summaries (RIS, same
`vcf.sh`).

## Transfer check

No `transfer/reviewer_gpu_20261017.{txt,sha256}` arrived, so the transfer was checked locally (job 423724,
`tmp/reviewer_cpu_20261017/ris_check/check.sh`): `gzip -t` on every `.gz`, `json.load` on every `.json` and
`torch.load` on every `.pth` that changed locally on 2026-10-17. Result (16 min): 16,824 files; 5,944 `.gz`, 536
`.json` and 58 `.pth` checked, none bad. This checks that every file is complete and readable, not that it equals the
RIS copy (no checksum list).

## Results (`tables.md`)

### 1. Retrained without channel 6 (`--drop-planes path_count`, E1.1 `gpu_nopc_commands.sh`)

Pairs differ only in channel 6 (same recipe and GPU count; d9 HG008T PacBio / ONT on 1 GPU, COLO829T d9 against the
`*_SNV_h90` runs, COLO829T d9 INDEL against the new with-ch6 controls `COLO829_*_INDEL_7fc7396`).

- **SNV:** post-filter PASS F1 is lower without channel 6 in 11 of 12 HPRC models (−0.006 to −0.311; HG008T ONT d46
  +0.009). Before filtering the loss is large everywhere: chr1 somatic AP 0.199–0.353 → 0.065–0.089 on HG008T,
  0.755–0.829 → 0.484–0.592 on COLO829T d46. Population filtering removes most of what channel 6 removes, so after it
  the median loss is 0.03 F1, with COLO829T Illumina d46 (−0.311; AP 0.818 → 0.582) and ONT d46 (−0.062) the
  exceptions. This is the retraining counterpart of E1.1's inference-time ablation (same conclusion: the model uses
  channel 6 as a germline prior that partly overlaps with the filter).
- **INDEL:** PASS F1 −0.030 to +0.016, best F1 −0.051 to +0.013, AP nearly unchanged (0.012–0.123 with channel 6, 0.004–0.113 without):
  channel 6 does little for indels.

### 2. Matched d9 SNV recipe (4 GPUs, batch 1024, the d46 / GRCh38-only recipe)

E3 caveat 1 was that 8 of the 12 d9 models used another recipe. With the matched recipe (`*_SNV_7fc7396_ddp4_b1024`;
HG008T Illumina already was):

| set | E3 d9 | matched d9 | d46 | GRCh38-only |
|---|---:|---:|---:|---:|
| HG008T PacBio | 0.621 | 0.608 | 0.605 | 0.513 |
| HG008T ONT | 0.596 | 0.600 | 0.575 | 0.537 |
| COLO829T Illumina | 0.778 | 0.572 (AP 0.736) | 0.783 | 0.468 |
| COLO829T fiberseq | 0.781 | 0.752 | 0.779 | 0.680 |
| COLO829T ONT | 0.706 | 0.744 | 0.737 | 0.418 |

The matched d9 still beats GRCh38-only on all six (+0.064 to +0.326). COLO829T Illumina's matched d9 has a much lower
PASS F1 (0.572) at a near-normal AP (0.736): 1,815 false calls after the PoN against 284 for the 1-GPU run, most of
them at chr1:121–126 and 143 Mb on tensors outside the labels' confident region (RIS summary).

### 3. Second training seed (`*_s1`), SNV, 17 pairs

- PASS F1 s1 − s0: −0.082 to +0.079 (|median| 0.021). chr1 AP s1 − s0: within ±0.039. So a second seed changes the
  ranking quality little, but the recall-0.9 validation threshold moves enough to change PASS F1 by up to 0.08,
  largest on COLO829T Illumina (d46 −0.080 at AP 0.818 vs 0.824) and the GRCh38-only arms (−0.082 to +0.079).
- HPRC minus GRCh38-only, same seed (`tables.md` 3b): positive in 20 of 21 comparisons; the exception is HG008T PacBio
  d9 seed 1 (−0.019; d46 seed 1 +0.012). Two-seed means: d9 +0.038 to +0.172, d46 +0.052 to +0.317 (COLO829T ONT
  GRCh38-only has no second seed). d46 − d9 two-seed mean: −0.007 to +0.014 on five sets, +0.145 on COLO829T Illumina
  (the low matched-d9 threshold).
- Consequence: single-run differences below about 0.05 PASS F1 are within training variance (R3.9); E16's
  block-bootstrap intervals cover chr1 sampling only.

### 4. Per-dataset normalisation (`*_d46_normds`)

Same as the training statistics on Illumina and PacBio / fiber-seq (SNV −0.006 to +0.000, INDEL −0.001 to −0.000), lower
on ONT: HG008T −0.027 (AP 0.328 → 0.266), COLO829T −0.217 (AP 0.755 → 0.433). This is not a model effect: `predict.py`
fits the dataset statistics on the predicted chromosomes only (chr1), training fitted them on chr2–22, and ONT chr1's
MAPQ mean is 0.50–0.51 training SD lower than chr2–22's (RIS summary, section "`--norm dataset`: predict fits the
statistics on the predicted chromosomes only"). The fix (fit on the whole set in predict and train) is a pipeline change
pending the user; then the 8 models are re-predicted. R3.m1's normalisation line stays [待补] until then.

### 5. Per-read BQ centring, COLO829T Illumina joint INDEL (`joint_readnorm`)

| purity | joint PASS / best | readnorm PASS / best |
|---|---|---|
| P100 | 0.129 / 0.167 | 0.077 / 0.182 |
| P75 | 0.201 / 0.312 | 0.070 / 0.146 |
| P50 | 0.238 / 0.330 | 0.057 / 0.133 |
| P25 | 0.195 / 0.233 | 0.037 / 0.133 |

With the read-source BQ cue removed, the dilutions no longer score above P100 (best F1 0.133–0.146 against 0.182),
which is what the artefact explanation (benchmark README) predicts; their earlier scores were inflated by the cue. PASS
precision collapses at P75–P25 (0.021–0.040), so the recall-0.9 threshold is not usable for these models. SNV and the
ONT joint models were not rerun with readnorm.

### 6. Reproduction (`HG008_Illumina_INDEL_7fc7396_repro`)

Retrained with the same command: identical PASS (0.264 / 0.238 / 0.251) and best F1 (0.269): training is deterministic
for a fixed seed and hardware.

## Caveats

- Two of the three planned seeds so far; seed variance per arm is from two values (the RIS summary will give
  three-seed means once seed 2 is synced).
- The COLO829T Illumina d46 nopc drop (−0.311) and the normds ONT drops are single runs; with the seed variance of 3,
  only the AP changes (−0.236, −0.322) are clearly outside it.
- Section 2's E3 d9 runs and the matched ones are both single seeds.
