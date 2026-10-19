# Benchmarks of the 7fc7396 pipeline on RIS (2026-10-12)

The SNV and INDEL recipes of `machine_learning/pansoma_net_v2/README.md` at 7fc7396 (INDEL PoN homopolymer rule,
INDEL lr 1e-4, INDEL tests counting partial matches), run on RIS compute1's `epigenome` H100s with the `ris/` scripts.
Every model is trained on chr2-22 and tested on chr1 (rtg vcfeval in the benchmark BED and without it, after the PoN
of each kind; PASS at the checkpoint's validation recall-0.9 threshold). Multi-GPU runs use DDP with `--batch-size
1024` per GPU (global batch N x 1024).

| folder | what | status |
|---|---|---|
| `hg008t_snv_indel_7fc7396_20261012/` | HG008T Illumina / PacBio / ONT, SNV and INDEL, HPRC v1.1 d9 graph | done |
| `colo829t_purity_joint_7fc7396_20261012/` | COLO829T purity dilutions: one model per platform and kind on P100 + P75 + P50 + P25, tested per purity | done |
| `hprc_v21_d46_7fc7396_20261012/` | the same 12 sample x platform x kind runs on the HPRC v2.1 d46 graph, against d9 | done |
| `linear_grch38_7fc7396_20261012/` | the same 12 on a GRCh38-only graph, against d9 and d46 | done |
| `indel_pon_71cb51a_20261016/` | every INDEL call set re-scored with the 71cb51a INDEL PoN rule (HP >= 8: AF 0.01, else 1e-4) | done |
| `reviewer_gpu_7fc7396_20261016/` | the RIS GPU jobs of the 2026-10-15 reviewer response: matched d9, no ch6, `--norm dataset`, seeds, E13 | done |

Each folder has a `README.md` (what was run, LSF jobs, results and their reading), a `collect.py` (standard library;
reads `$PANSOMA_RUNS/<run>` and writes `tables.md`, `summary.json` and `runs/<run>/` copies of the training logs and rtg
reports) and `tables.md`.

## Main results (chr1, in the BED, after the PoN)

SNV: PASS F1 at the validation recall-0.9 threshold. INDEL: the best F1 of the chr1 rtg curve, because the recall-0.9
threshold lands far below the best point for INDELs (PASS F1 in brackets). The best point is chosen on the test
chromosome: an upper bound for a threshold rule, not a held-out result. INDEL numbers here use the 71cb51a INDEL PoN
rule (homopolymer >= 8 bp: AF >= 0.01, else 1e-4; `indel_pon_71cb51a_20261016/`, +0.000 to +0.013 against the 7fc7396
rule, no comparison changes direction); the sub-folders' tables keep the 7fc7396 rule they were written with.

**HG008T, v1.1 d9** (`hg008t_snv_indel_7fc7396_20261012`):

| platform | SNV (PASS) | INDEL (best; PASS) |
|---|---:|---:|
| Illumina | 0.585 (4 GPUs) | 0.269 (0.251) |
| PacBio | 0.621 | 0.263 (0.092) |
| ONT | 0.596 | 0.143 (0.015) |

SNV reproduces the 2026-10-03 runs exactly (PacBio, ONT); INDEL lr 1e-4 gives the same as lr 4e-4 under the same PoN,
and the 7fc7396 homopolymer PoN rule raised Illumina INDEL 0.218 -> 0.246 (PASS). On the long reads the recall-0.9
threshold passes 2,780 / 17,402 INDEL calls after the PoN (7fc7396 rule), far from the best points.

**COLO829T purity, joint models** (`colo829t_purity_joint_7fc7396_20261012`), P100 / P75 / P50 / P25:

| platform | SNV (PASS) | INDEL (best) |
|---|---|---|
| Illumina | 0.675 / 0.840 / 0.833 / 0.804 | 0.167 / 0.312 / 0.330 / 0.233 |
| fiberseq | 0.779 / 0.687 / 0.665 / 0.720 | 0.069 / 0.049 / 0.055 / 0.074 |
| ONT | 0.677 / 0.654 / 0.639 / 0.666 | 0.055 / 0.051 / 0.070 / 0.062 |

- Illumina P100 is low because one validation threshold is pooled over the four purities (P100 best point 0.784).
- The Illumina and ONT dilutions carry a read-source artefact: the COLO829BL normal reads have higher base quality
  than the tumor reads (Illumina read groups: mean BQ 30.77 vs 29.43, nothing else differs), so the ALT reads' BQ
  marks somatic variants (AUC of that difference alone up to 0.71 Illumina, 0.69 ONT at P25); their scores are likely
  optimistic. fiberseq has no such cue. A per-read BQ centering removes it in the tensors; the fix is open.
- With BQ normalized per read (`ris/tools/readnorm.py`), Illumina SNV becomes 0.739 / 0.709 / 0.678 / 0.695 (precision
  0.68-0.82 instead of 0.97; INDEL best 0.182 / 0.146 / 0.133 / 0.133): the dilution scores were the read-source cue.
- On every platform P75 and P50 score below P25 (the anomaly of `analysis/reviewer_response_20261011`, reproduced;
  Illumina once the BQ cue is removed): germline SNVs of the tumor's LOH regions appear at AF (1 - p) / 2 = 0.125 /
  0.25 / 0.375 where the tumor lost the ALT copy; at P50 that lies inside the somatic AF (the most germline false
  calls after the PoN, 463-523), at P75 inside the P25 somatic AF the joint model also trains on, at P25 above the
  somatic AF (`loh_fp.md`). A real LOH tumor at intermediate purity has the same germline; a single site cannot
  separate them.
- INDEL: Illumina dilutions above P100 (0.312 / 0.330 / 0.233 vs 0.167, likely helped by the same BQ cue); long reads
  0.05-0.07 at every purity (4-17 true calls at the best points), below the single-set P100 models (0.124 / 0.074).

**HPRC v2.1 d46 against v1.1 d9** (`hprc_v21_d46_7fc7396_20261012`): after the PoN the same (PASS F1 -0.021 to +0.006
on 11 of 12 sets, COLO829T ONT SNV +0.031); before the PoN d46 is much better for SNV (raw PASS F1 +0.10-0.13 on all six
SNV sets), because fewer germline variants become candidates on the larger pangenome, which is what the PoN removes
anyway. INDEL: d46 gives more truth INDELs a tensor, but they are population alleles the PoN tags (the PoN-limited
ceilings are unchanged); by best F1, d46 - d9 is -0.013 / +0.012 / +0.016 on HG008T Illumina / PacBio / ONT and
-0.034 / +0.013 / +0.015 on COLO829T (against the 2026-10-06 lr 4e-4 d9 runs; +0.007 / -0.005 / -0.006 against the
matched lr 1e-4 d9 runs of `reviewer_gpu_7fc7396_20261016`): equal within noise. With three seeds per arm and the d9 SNV
models on 4 GPUs as d46 (`reviewer_gpu_7fc7396_20261016`), d46 - d9 SNV is -0.017 to +0.015 (PASS F1, three-seed means)
on five of six sets, inside the seed spread; COLO829T Illumina +0.163, where the 4-GPU d9 runs are low on every seed.

**Linear GRCh38 graph** (`linear_grch38_7fc7396_20261012`), after the PoN, linear / d9 / d46:

| sample | SNV (PASS) Illumina / PacBio-fiberseq / ONT | INDEL (best) Illumina / PacBio-fiberseq / ONT |
|---|---|---|
| HG008T | 0.469 / 0.513 / 0.537 vs 0.585 / 0.621 / 0.596 (d9) | 0.201 / 0.176 / 0.084 vs 0.269 / 0.263 / 0.143 (d9) |
| COLO829T | 0.468 / 0.680 / 0.418 vs 0.778 / 0.781 / 0.706 (d9) | 0.070 / 0.095 / 0.020 vs 0.202 / 0.124 / 0.074 (d9) |

The pangenome is worth 0.06-0.31 SNV F1 even after the PoN on one seed, +0.03 to +0.29 (PASS F1) and +0.02 to +0.21
(best F1) as three-seed means of every arm (`reviewer_gpu_7fc7396_20261016`; smallest on HG008T PacBio) (linear: every germline variant of the sample is a
candidate, raw F1 0.01-0.07; more false calls survive the PoN, and many more outside the BED). For INDEL the linear
graph keeps the somatic INDELs a pangenome absorbs (truth with a tensor 0.75-0.87 vs 0.54-0.70 on HG008T) but ranks
worse: best F1 below d9 / d46 on every INDEL set (its higher long-read PASS F1 comes from where its recall-0.9
threshold lands). Overall: linear < d9 = d46 after the PoN, for SNV and INDEL; before it, linear << d9 < d46.

## Data and code

- Tensor sets: `/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors/<SET>[_P75|_P50|_P25|_d46|_linear]`
  (labels `5430533aa100`), each checked before use (every manifest shard present with the right size, labels and
  summaries present).
- Graph indexes in `Dataset/pansoma_v2_tensors/graph_index/` (d9; d46 and linear linked from
  `Dataset/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/` and `Dataset/Linear_Reference_VG_Indexes/graph_index/`, both
  checked with `ris/tools/check_graph_index.py`).
- Runs: `/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs/<run>`; scripts `ris/` (`snv_pipeline.sh`,
  `joint_pipeline.sh`, `joint_run.sh`, `recipes.sh`, `vcf.sh [--set]`, `vcfeval_only.sh`, `prebuild.sh`).
