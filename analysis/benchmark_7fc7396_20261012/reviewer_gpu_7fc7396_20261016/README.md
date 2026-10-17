# RIS GPU jobs of the 2026-10-15 reviewer response (submitted 2026-10-16)

The GPU (RIS) list of `analysis/reviewer_response_20261011/overview_20261015.md` §3.1 and of the `[GPU RIS]` items of
its `points/*.md`, as far as the data on RIS allow. Every run is `ris/snv_pipeline.sh` (train with `ris/recipes.sh`,
predict chr1, then `vcf.sh`: graph VCF, GRCh38 VCF, PoN, rtg vcfeval in and outside the BED) on the epigenome H100s:
SNV on 4 GPUs (DDP, `--batch-size 1024` per GPU), INDEL on 1, unless stated. Inference on other sets is for tequila.

| group | runs | what | response | script |
|---|---:|---|---|---|
| 0 repro | 1 | `HG008_Illumina_INDEL_7fc7396` retrained as it was | R3.9 | `submit.sh` |
| 1 matched d9 | 8 | HG008T PacBio / ONT, COLO829T x 3 SNV on 4 GPUs (as d46 / linear); COLO829T x 3 INDEL at lr 1e-4 | R3.1, R1.2, R2.5, R3.2, R1.6 | `submit.sh` |
| 2 no ch6 | 24 | `experiments/E1.1_ch6_counterfactual/gpu_nopc_commands.sh` (`--drop-planes path_count`), d9 and d46, SNV and INDEL | R1.1, R2.5, R3.1, R3.2 | `submit.sh` |
| 3 `--norm dataset` | 8 | d46 SNV x 6 sets, Illumina INDEL x 2, for the cross-sample tests (E12) | R1.7, R1.m1, R3.m2 | `submit.sh` |
| 4 seeds | 36 | the main SNV model of linear / d9 (group 1) / d46 x 6 sets, seeds 1 and 2 | R3.9, R3.1 | `submit.sh` |
| 5 E13 | 10 | learning curve (10 / 25 / 50 / 75 % of the training blocks x 2 replicates) and model size (8.7 M, 199.6 M) on `HG008T_PacBio_d46` | R1.7 | `submit_e13.sh` |

- `--norm train` in every group but 3. Every control (seed 0, with ch6, the 100 % / 50.4 M point, the d46 and linear
  runs) was trained before `--norm` existed, i.e. with the training statistics; with `--norm train` a pair differs only
  in the factor tested.
- Group 0 checks that: the current code with `--norm train` reproduces `HG008_Illumina_INDEL_7fc7396` bit for bit (all
  12 epochs' metrics, best epoch 3, threshold 0.39199385046958923, and all 12,674 chr1 `p_somatic` identical).
- Group 2 keeps the names, GPU counts, slots and memory of `gpu_nopc_commands.sh`: the d9 SNV no-ch6 runs other than
  HG008T Illumina are on 1 GPU, paired with `HG008_*_SNV_7fc7396` / `COLO829_*_SNV_h90`. Its three COLO829T d9 INDEL
  with-ch6 controls (`COLO829_*_INDEL_7fc7396`) are group 1.
- The ch6 test (user, 2026-10-16): on the P100 sets only, main arm HPRC v2.1 d46; COLO829T models are the single-set
  P100 ones, not the purity-joint models (which E1.1's d9 COLO829T rows used). The d9 no-ch6 runs are kept and run to
  the end; the d46 ones were moved to the front of the queue (btop) while d9 ones were still pending.
- Group 5: the exclude files of E21 README §6 are remade on RIS (`e13_subsets.sh` / `e13_subsets_ris.py`, which runs
  E21's `make_e13_subsets.py` unchanged on the RIS copy of the set; SHA-1 draws over the seed-0 split) into
  `$PANSOMA_RUNS/e13_subsets/`. Checked: `e13_subsets_HG008T_PacBio_d46_SNV.tsv` here equals E21's in all 9 rows (kept
  blocks, excluded ids, training tensors per class, truth alleles; every column but `file`). `e13_exclude.sha256`: the
  eight files, for a byte comparison with tequila's
  (`sha256sum -c` in `/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/e13`).
- Run names: group 1 `<S>_<P>_SNV_7fc7396_ddp4_b1024`, `COLO829_<P>_INDEL_7fc7396`; group 2 `<S>_<P>_{SNV,INDEL}_nopc_7fc7396`
  (HG008T Illumina SNV `..._ddp4_b1024`), `<S>_<P>_{SNV,INDEL}_nopc_d46_7fc7396`; group 3 `<S>_<P>_SNV_d46_normds`,
  `<S>_Illumina_INDEL_d46_normds`; group 4 `<d9 / d46 / linear run>_s1`, `_s2`; group 5 `HG008_PacBio_SNV_d46_e13_*`.
  LSF job ids: `$PANSOMA_RUNS/logs/submit_reviewer_gpu_20261016.txt`.

## Results (2026-10-17 12:00; `tables.md` from `collect.py`)

chr1 after the PoN, in the benchmark BED (COLO829T: the full SMaHT union). SNV: PASS F1 at the validation recall-0.9
threshold, best F1 of the chr1 curve beside it; INDEL: best F1, with the 71cb51a INDEL PoN rule. Done: groups 0-3 and
17 of the 18 seed-1 runs; seed 2 and E13 running / queued.

- **Repro**: bit-identical (above).
- **Matched d9** (SNV on 4 GPUs as d46): d46 - d9 on one seed is -0.015 / -0.003 / -0.025 on HG008T Illumina / PacBio /
  ONT SNV and +0.212 / +0.026 / -0.007 on COLO829T (the matched COLO829T Illumina d9 run has 1,815 false calls after
  the PoN against 284 for the earlier 1-GPU run, most of them in chr1:121-126 and 143 Mb, tensors outside the labels'
  confident region); INDEL +0.008 / -0.005 / -0.006 on COLO829T. 1 -> 4 GPUs on the same recipe changes d9 SNV by
  -0.014 / +0.004 (HG008T PacBio / ONT) and -0.206 / -0.029 / +0.038 (COLO829T), as large as the seed spread below.
- **No ch6** (24 pairs, one seed each): SNV d9 -0.006 to -0.043 on all six sets; d46 -0.007 to -0.311 on five, +0.009
  on HG008T ONT; INDEL HG008T -0.001 to -0.034 on all six, COLO829T -0.051 to +0.013 (11-26 true calls at the best
  points). 20 of 24 pairs lower without ch6. The seed spread below (PASS F1 up to 0.035 between two seeds on the
  HG008T graph arms, 0.08 on HG008T linear and on COLO829T Illumina d46) is of the same size as most single pairs,
  so the per-pair numbers need no-ch6 seeds before they are quoted as effects; the direction (lower on 20 of 24) is
  the robust part.
- **`--norm dataset`**: in-sample within -0.006 to +0.004 of `--norm train` on six models; the two ONT SNV models
  lower because predict fitted the statistics on chr1 only (next section). Cross-sample predictions wait for that fix.
- **Seeds** (two of three seeds so far, SNV; `tables.md`):
  - d46 - d9 changes sign between the two seeds on HG008T Illumina (-0.015 / +0.002), PacBio (-0.003 / +0.032), ONT
    (-0.025 / +0.015) and COLO829T fiberseq (+0.027 / -0.008); it keeps its sign on COLO829T Illumina (+0.211 / +0.078;
    both 4-GPU d9 runs low, 0.572 / 0.626 against 0.778 for the 1-GPU run) and ONT (-0.007 / -0.031).
  - PASS F1 moves more than the best F1 between seeds, most on the linear arm (HG008T PacBio 0.513 / 0.592, best
    0.596 / 0.604), whose threshold is chosen on a validation somatic AP of 0.012 / 0.011 (validation F1 0.04).
  - By best F1 the graphs stay above linear on every set and seed so far: HG008T Illumina +0.06-0.07, ONT +0.03-0.05,
    PacBio +0.00-0.04; COLO829T Illumina +0.10-0.23, fiberseq +0.08-0.10, ONT +0.23 (seed 0).
  - The README's "pangenome worth 0.06-0.31 SNV F1 after the PoN" is a one-seed PASS F1 statement; with the seeds the
    HG008T PacBio PASS gap goes from +0.09 (seed 0) to -0.02 / +0.01 (seed 1). Final numbers: three-seed means.
- **E13**: queued.

## Code updates while the jobs ran

The RIS checkout was pulled to 689fc2f / 2c39bad on 2026-10-16 18:00-18:02, so jobs that started later run that code:
- `train --set-test-chroms` (689fc2f): without the option the split, validation and sampling are those of before
  (`chroms_of` returns the chromosomes unchanged), so the runs here are unaffected.
- The INDEL PoN rule of 71cb51a (homopolymer >= 8 bp: AF >= 0.01, else 1e-4; was >= 7 bp and 1e-3; SNV unchanged): the
  `vcf.sh` of an INDEL run here has the old or the new rule depending on when it ran, so every INDEL call set (these and
  the benchmark's) was re-scored with the new rule into `<run>/vcf_chr1_pon71cb51a` (`../indel_pon_71cb51a_20261016/`);
  compare INDEL runs on those.
- To do once no job uses it: the comment in `ris/vcf.sh` (line 45) still states the 7fc7396 INDEL rule.

## `--norm dataset`: predict fits the statistics on the predicted chromosomes only

The 8 `--norm dataset` models score as their `--norm train` counterparts in-sample except the two ONT SNV models
(`tables.md`): COLO829T ONT 0.737 -> 0.520, HG008T ONT 0.575 -> 0.548 (PASS F1), with the same validation F1 (0.737 /
0.739 on COLO829T ONT). `predict.py` fits the dataset statistics on `positions`, the tensors of `--chroms` (here chr1
only), while `train.py` fitted them on the training chromosomes (chr2-22). chr1 is not a typical sample of a set: its
MAPQ mean is 0.50 / 0.51 training SD lower on the two ONT sets (COLO829T ONT 48.6 -> 39.7, std 17.7 -> 19.8), so every
chr1 tensor is z-scored differently from training, and the chr1 somatic AP falls 0.755 -> 0.433 / 0.328 -> 0.266.
The other six have MAPQ shifts of -0.37 to +0.01 SD and change by -0.006 to +0.004.

| model | BQ mean, training -> chr1 (shift in SD) | MAPQ mean, training -> chr1 (shift in SD) | chr1 somatic AP, train -> dataset norm |
|---|---|---|---|
| HG008T Illumina SNV | 33.3 -> 33.3 (-0.01) | 56.1 -> 54.3 (-0.16) | 0.310 -> 0.309 |
| HG008T PacBio SNV | 38.1 -> 38.5 (+0.06) | 42.1 -> 34.5 (-0.37) | 0.353 -> 0.350 |
| HG008T ONT SNV | 41.2 -> 42.9 (+0.16) | 47.1 -> 38.1 (-0.50) | 0.328 -> 0.266 |
| COLO829T Illumina SNV | 28.7 -> 28.5 (-0.04) | 54.8 -> 52.7 (-0.17) | 0.818 -> 0.819 |
| COLO829T fiberseq SNV | 37.8 -> 38.1 (+0.05) | 39.7 -> 35.7 (-0.19) | 0.829 -> 0.823 |
| COLO829T ONT SNV | 29.8 -> 30.1 (+0.03) | 48.6 -> 39.7 (-0.51) | 0.755 -> 0.433 |
| HG008T Illumina INDEL | 34.1 -> 34.0 (-0.02) | 56.3 -> 56.4 (+0.01) | 0.123 -> 0.123 |
| COLO829T Illumina INDEL | 28.9 -> 28.7 (-0.03) | 55.9 -> 55.6 (-0.03) | 0.057 -> 0.058 |

The checkpoints are not affected (training fitted the statistics on chr2-22); their chr1 predictions are, and so would
be any cross-sample prediction of one chromosome of another tumor. A fix is a pipeline change (pending the user):
fit the statistics of a dataset on the whole set (all its chromosomes, labels not used) in predict, and likewise in
train, so in-sample training and test use the same statistics; then re-predict the 8 models.

## Not run here

- HCC1395T / CASTLE (HCC1954T, H1437T, H2009T) in-sample models and leave-one-line-out (R1.7, R1.m1, R1.6): their
  tensor sets are on tequila only. A d46 set is 0.45-0.65 TB (SNV 430 GB, INDEL 22 GB), the twelve about 6 TB; storage2
  fs1 had 7.9 TB free (99 % used) on 2026-10-16.
- The purity joint models after the BQ read-source fix, with per-purity thresholds (E14b; R1.5, R2.4, R3.10): the fix is
  not decided (the per-read BQ experiment is `colo829t_purity_joint_7fc7396_20261012`).
- User decisions: E9f (off-reference tensors labelled by the assembly), E15 low-AF build, E14c / E14d second tumor, the
  three-fold chromosome cross-validation (R3.9, optional), HapMap (E23).
