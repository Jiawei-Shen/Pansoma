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

## Results (all 87 runs, 2026-10-19; `tables.md` from `collect.py`)

chr1 after the PoN, in the benchmark BED (COLO829T: the full SMaHT union). SNV: PASS F1 at the validation recall-0.9
threshold, best F1 of the chr1 curve beside it; INDEL: best F1, with the 71cb51a INDEL PoN rule.

- **Repro**: bit-identical (above).
- **Seeds** (SNV, three seeds per arm and set; mean and range in `tables.md`):
  - Seed spread, PASS F1: 0.002-0.036 on the HG008T graph arms, 0.005-0.079 on the COLO829T graph arms, 0.020-0.144 on
    linear; best F1: at most 0.019 on the HG008T graph arms, 0.065 on COLO829T, 0.061 on linear. The PASS F1 of one
    run is a noisy readout; linear's threshold is chosen on a validation somatic AP of about 0.01.
  - d46 - d9, three-seed means, PASS (best): HG008T Illumina -0.011 (+0.004), PacBio +0.015 (-0.005), ONT +0.002
    (+0.004); COLO829T fiberseq +0.012 (+0.007), ONT -0.017 (-0.014): equal within the seed spread. COLO829T Illumina
    +0.163 (+0.091): the 4-GPU d9 runs are low on all three seeds (PASS 0.566-0.626; 1,340-1,900 false calls after
    the PoN, 1,085-1,659 of them in chr1:121-126 and 143-144 Mb, tensors outside the labels' confident region), while
    the 1-GPU d9 run scored 0.778 with 284 (33 there), as d46 (0.751 mean): a systematic effect of the 4-GPU recipe on
    d9 for this set, not a one-seed accident.
  - Graph - linear, three-seed means: PASS +0.028 to +0.292 (HG008T +0.028 to +0.082, COLO829T +0.077 to +0.292),
    best +0.015 to +0.211 (HG008T +0.015 to +0.065, COLO829T +0.076 to +0.211); positive on every set for both graphs.
    The benchmark README's one-seed "0.06-0.31 after the PoN" becomes +0.03 to +0.29 (PASS) / +0.02 to +0.21 (best);
    the smallest gain is HG008T PacBio (+0.03 / +0.04 PASS, +0.02 / +0.02 best), where linear's seed 0 (0.513) was low.
- **No ch6** (24 pairs, one seed each; on d46 the with-ch6 side has the three seeds, same 4-GPU recipe):
  - On d46 every no-ch6 SNV model is below the with-ch6 three-seed range by best F1 (HG008T Illumina -0.071, PacBio
    -0.015, ONT -0.009; COLO829T Illumina -0.173, fiberseq -0.006, ONT -0.059 against the mean) and on five of six by
    PASS F1 (HG008T ONT inside: 0.584 against 0.575-0.598). HG008T Illumina d9 (4 GPUs) is below its range too (PASS
    0.542 against 0.585-0.587).
  - The other pairs (one seed per side; d9 SNV on 1 GPU): d9 SNV -0.006 to -0.041; INDEL HG008T d9 -0.010 / -0.034 /
    -0.011, d46 -0.004 / -0.030 / -0.001; COLO829T INDEL -0.051 to +0.013 on 11-26 true calls at the best points.
  - 20 of 24 pairs are lower without ch6. The no-ch6 models have one seed each, so the size per set is uncertain by
    about the seed spread; the direction is consistent.
- **Matched d9** (SNV on 4 GPUs as d46): 1 -> 4 GPUs on the same recipe changes d9 SNV by -0.014 / +0.004 (HG008T PacBio
  / ONT) and -0.206 / -0.029 / +0.038 (COLO829T); the three-seed comparison of the arms is the one above. COLO829T d9
  INDEL at lr 1e-4: best 0.161 / 0.142 / 0.095 (Illumina / fiberseq / ONT), d46 - d9 +0.008 / -0.005 / -0.006.
- **`--norm dataset`**: in-sample within -0.006 to +0.004 of `--norm train` on six models; the two ONT SNV models
  lower because predict fitted the statistics on chr1 only (next section). Cross-sample predictions wait for that fix.
- **E13** (HG008T PacBio d46, SNV; PASS F1 (best F1), two replicates per fraction; 100 % = the three seeds):
  - Training blocks 10 % 0.574 (0.579), 25 % 0.582 (0.588), 50 % 0.597 (0.600), 75 % 0.600 (0.616), 100 % 0.607
    (0.616): rising with the data, by +0.03 (PASS) from 10 % to 100 % and +0.010 from 50 % to 100 %; the replicates
    agree within 0.001-0.009.
  - Model size on all blocks: 8.7 M 0.616 (0.632), 50.4 M 0.607 (0.616; seeds 0.604-0.611), 199.6 M 0.603 (0.608): no
    gain from capacity; the 8.7 M model is at or above the 50.4 M seed range (one run).

## Code updates while the jobs ran

The RIS checkout was pulled to 689fc2f / 2c39bad on 2026-10-16 18:00-18:02, so jobs that started later run that code:
- `train --set-test-chroms` (689fc2f): without the option the split, validation and sampling are those of before
  (`chroms_of` returns the chromosomes unchanged), so the runs here are unaffected.
- The INDEL PoN rule of 71cb51a (homopolymer >= 8 bp: AF >= 0.01, else 1e-4; was >= 7 bp and 1e-3; SNV unchanged): the
  `vcf.sh` of an INDEL run here has the old or the new rule depending on when it ran, so every INDEL call set (these and
  the benchmark's) was re-scored with the new rule into `<run>/vcf_chr1_pon71cb51a` (`../indel_pon_71cb51a_20261016/`);
  compare INDEL runs on those.
- The comment in `ris/vcf.sh` stating the INDEL rule was updated once no job used the script (2026-10-19).

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
