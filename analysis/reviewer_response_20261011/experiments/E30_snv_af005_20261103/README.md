# E30: SNV tensor AF floor 0.05 on the CASTLE PacBio HiFi sets (2026-11-03 / 04)

Work directory `tmp/snv_af05_20261103` (job scripts, logs, `code/` = pipeline archive, `code_ml/` = archive of
96e2d6c); this directory keeps the scripts, the counts, the checks and the results.

## Why

E27 `v3/hcc1395_snv_gap`: of the 450 HCC1395T PacBio HiFi chr1 truth SNVs only DeepSomatic-TO finds (against the
six-tumor model's 42), 201 have a Pansoma candidate below the SNV tensor AF floor 0.08, 71 fewer than 3 ALT reads,
11 no candidate; 167 have a tensor. User 2026-11-03: rebuild with the floor at 0.05.

## Candidates the floor removes (`candidates/count.py`, `filtered_candidates.ndjson`; AF = the prefilter's upper bound)

SNV alleles filtered only by the AF floor (>= 3 ALT reads), chr1, AF 0.05-0.08, and the BED somatic truths among them:

| set | chr1 SNV tensors now | alleles at AF 0.05-0.08 | BED somatic truths | truths / alleles |
|---|---|---|---|---|
| HCC1395T PacBio HiFi | 105,386 | 8,123 (+7.7%) | 167 | 2.1% |
| HCC1395T Illumina | 268,264 | 290,746 (+108%) | 144 | 0.05% |
| HCC1395T ONT | 276,601 | 199,177 (+72%) | 114 | 0.06% |
| HG008T PacBio HiFi | 140,610 | 10,032 (+7.1%) | 1 | 0.01% |
| COLO829T Fiber-seq | 241,213 | 13,907 (+5.8%) | 100 | 0.7% |

37 HCC1395T PacBio truths lie below 0.05 (discovery selects nodes at imperfect fraction > 0.05, so 0.05 is the lowest
floor without a new discovery). DeepSomatic-TO's own calls on HCC1395T chr1 at its VAF 0.05-0.08 (after the PoN):
PacBio HiFi 184 TP / 54 FP, Illumina 143 / 203, ONT 100 / 113; below 0.05 (PacBio) 93 / 290.

Decisions (user 2026-11-03): the four CASTLE tumors only (HCC1395T, HCC1954T, H1437T, H2009T PacBio_d46); HG008T and
COLO829T stay as they are (HG008T's assembly-based truth has almost no low-VAF SNVs: 1 truth in 10,032 candidates).
Illumina and ONT not rebuilt (the band is almost all errors there).

## Sets `pansoma_v2_tensors/<SET>_snvaf05`

`make_prepare.py`: the old run's config (GAM, the old discovery's node list and stats, haplotypes 464, task count,
processes, every builder parameter, chromosome index, merge, reference path) and the truth of the set's current SNV
labels, with `--snv-min-af 0.05` (INDEL 0.08). Code: git archive of 0136760 + the native decoder. `check_config.py` at
the end of prepare: no difference but the SNV floor, the (old) node list path and the QUAL >= 10 germline truth of the
current labels. Chains `submit_set.sh` (prepare -> run -> finalize -> SNV site features), `jobs.tsv`; runs 1.5-2.3 h on
48 processes, finalize ~1.5 h. The unused new INDEL tensors were deleted.

Checks (`check_tensors.py SET chr22 chr1`, `results/check_*.txt`): every old site present; same-allele tensors
byte-identical (chr22 all, chr1 a random 30,000: 0 differ in each set); labels of all same-allele sites identical; new
alleles all at AF 0.0500-0.0799. chr1 sites +6.0-6.8% (HCC1395T 105,386 -> 112,600); most added sites are -1 (outside
the labels' confident region, as for the old tensors); HCC1395T chr1 BED somatic truths with a new tensor 169.
Genome SNV tensors: HCC1395T 1,073,061 -> 1,154,325, HCC1954T 1,291,833 -> 1,391,153, H1437T 1,289,743 -> 1,396,847,
H2009T 1,319,176 -> 1,415,166.

## Models (`train.sh`, `submit_train.sh`)

The site-feature recipe of E29 (`train_sitef.sh`: 50.4 M, 20 epochs, pow:0.75, lr 1e-4, batch 1024 on one H100,
`--site-features --best-by postpon --call-recall 0.9`) with `--pon-negative-fraction 0.25` (user: train on part of the
negatives; repository 96e2d6c, reviewed): each epoch every somatic tensor and every non / germline tensor the PoN does
not tag, a fresh 25% of the PoN-tagged non / germline tensors (76% of the training tensors, almost all germline) with
loss weight 4; class weights from all tensors; validation unchanged. An epoch is 43% of the tensors.

| run | data | training time |
|---|---|---|
| multicancer6_PacBio_SNV_sitef_af05_pnf25 | the four `_snvaf05` sets + HG008T_PacBio_d46, COLO829T_fiberseq_d46 | 3 h 07 min |
| multicancer6_PacBio_SNV_sitef_pnf25 (control) | the six sets of sitef_ppr9 (floor 0.08), same sampling | 3 h 24 min |
| multicancer6_PacBio_SNV_sitef_ppr9 (E29) | floor 0.08, every negative | 5 h 56 min |

ppr9 vs pnf25 = the sampling; pnf25 vs af05_pnf25 = the floor. Tests: HCC1395T chr1 / chr20 on each model's own set
(`_snvaf05` for af05), HG008T and COLO829T chr1 on their unchanged sets.

## Results (`compare.py` -> `results/results_PacBio_af05.md`; evaluation fixed before the results)

SNV F1 (P / R) after the PoN, benchmark BED. Rule 1 = the PASS of the checkpoint (post-PoN validation callable recall
0.9), rule 2 = the same at 0.95 (E29, registered before the Illumina results; post hoc on PacBio HiFi):

| model | HCC1395T chr20 | HCC1395T chr1 | HG008T chr1 | COLO829T chr1 |
|---|---|---|---|---|
| AF 0.05 + sampling, rule 1 | 0.825 (0.845 / 0.806) | 0.856 (0.889 / 0.826) | 0.669 (0.617 / 0.730) | 0.755 (0.898 / 0.651) |
| AF 0.05 + sampling, rule 2 | 0.841 (0.830 / 0.854) | 0.865 (0.867 / 0.863) | 0.642 (0.540 / 0.793) | 0.778 (0.860 / 0.710) |
| 0.08 + sampling (control), rule 1 | 0.800 (0.851 / 0.756) | 0.832 (0.888 / 0.783) | 0.653 (0.620 / 0.690) | 0.747 (0.903 / 0.637) |
| 0.08 + sampling (control), rule 2 | 0.812 (0.824 / 0.801) | 0.843 (0.865 / 0.822) | 0.643 (0.549 / 0.776) | 0.774 (0.857 / 0.705) |
| 0.08, every negative (ppr9), rule 1 | 0.794 (0.847 / 0.747) | 0.828 (0.895 / 0.769) | 0.659 (0.614 / 0.710) | 0.777 (0.879 / 0.697) |
| 0.08, every negative (ppr9), rule 2 | 0.813 (0.832 / 0.795) | 0.845 (0.872 / 0.820) | 0.645 (0.549 / 0.780) | 0.792 (0.859 / 0.734) |
| normfix (no site features) | 0.788 (0.804 / 0.773) | 0.825 (0.843 / 0.809) | 0.625 (0.586 / 0.670) | 0.774 (0.880 / 0.691) |
| DeepSomatic-TO | 0.835 (0.740 / 0.957) | 0.839 (0.766 / 0.927) | 0.557 (0.425 / 0.806) | 0.632 (0.732 / 0.555) |
| ClairS-TO | 0.795 (0.778 / 0.812) | 0.814 (0.805 / 0.822) | 0.591 (0.455 / 0.844) | 0.762 (0.779 / 0.745) |

- The floor (af05 vs control): HCC1395T chr20 +0.025 (rule 1) / +0.029 (rule 2), chr1 +0.024 / +0.022, from recall
  (chr1 0.783 -> 0.826, 0.822 -> 0.863; chr20 0.756 -> 0.806, 0.801 -> 0.854) at the same precision (within 0.006).
  The F1 differences are inside the seed range, but the recall comes from truths the floor-0.08 models have no tensor
  for: of the AF 0.05-0.08 truths the af05 model calls 0.713 / 0.808 on chr1 (167) and 0.868 / 0.882 on chr20 (68);
  DeepSomatic-TO 0.994 / 0.985, ClairS-TO 0.162 / 0.221. HG008T and COLO829T (their sets unchanged): +0.016 / -0.001
  and +0.008 / +0.004.
- Against DeepSomatic-TO under rule 2: HCC1395T chr20 0.841 vs 0.835, chr1 0.865 vs 0.839 (under rule 1 0.825 and
  0.856). DeepSomatic's tumor-only models were trained on HCC1395 with this truth, chr1 held out (E8b D2, D10); chr20 is
  in its training chromosomes.
- The sampling (control vs ppr9): HCC1395T and HG008T within 0.006 under both rules; COLO829T -0.030 (rule 1) /
  -0.018 (rule 2), recall 0.697 -> 0.637 / 0.734 -> 0.705. Within the seed range; the only test that moved.
- HCC1395T recall by read AF (rule 2, af05 / control / ppr9): chr1 AF 0.6-0.85 0.846 / 0.832 / 0.871, >= 0.85 0.885 /
  0.885 / 0.901; the truths without any tensor fall from 420 to 252 (chr1) and 144 to 76 (chr20).

One run per model (seed 0). Rule 2 is post hoc for PacBio HiFi (registered in E29 after its PacBio results; the
Illumina site-feature model is its test, E29).

## Illumina and ONT (user 2026-11-04 13:20: test them too); registered 2026-11-04 13:30, before any of their results

The PacBio HiFi design on the other two platforms. Sets `<SET>_snvaf05` for HCC1395T / HCC1954T / H1437T / H2009T
Illumina_d46 and ONT_d46 (`make_prepare.py`, `submit_set.sh`; Illumina on maotai / meow / miao / tequila, each ONT
build after the Illumina one on the same node); HG008T and COLO829T unchanged. ONT site features of the six old sets
(tmp/snv_sitefeat_20261103/build_sitef.sh). Models (`train.sh DATA F RUN SEED PLAT`, `--pon-negative-fraction 0.25`):
multicancer6_<PLAT>_SNV_sitef_pnf25 (control, floor 0.08) and multicancer6_<PLAT>_SNV_sitef_af05_pnf25 (floor 0.05 on the
four CASTLE sets). Evaluation as for PacBio HiFi (`compare.py PLAT`): both rules, HCC1395T chr20 first, chr1 of the
three tumors, recall by read AF; normfix, DeepSomatic-TO, ClairS-TO; Illumina also sitef_ppr9 (E29).
Expectation from the counts (chr1, AF 0.05-0.08: Illumina 290,746 candidates for 144 BED truths, ONT 199,177 for 114,
against PacBio HiFi 8,123 for 167): at most +0.042 / +0.033 recall, and a precision cost; the test is whether the
models call the band without losing precision. One seed; a difference below the seed range is not an effect.
