# E29: SNV models with population site features (2026-11-03)

> **Status (2026-11-03 14:30).** PacBio: features of the six sets built, `multicancer6_PacBio_SNV_sitef_ppr9`
> training on a tequila H100 (job 467902; then predictions and the four vcf jobs 467903-467906). Illumina and ONT not
> started. The scripts run in `tmp/snv_sitefeat_20261103/` (logs, `code/`); copies here.

User 2026-11-03: give the SNV models the population-AF site features of the final INDEL models and select the
checkpoint and the PASS threshold after the PoN (the INDEL "ppr" rule), so that a model can learn that a high-AF
allele absent from the population databases is likely somatic. PacBio first (user), then Illumina and ONT.

Motivation (E27 `v3/hcc1395_snv_gap/recall_by_af.md`, chr1): on HCC1395T the six-tumor model's recall falls with the
read allele fraction (Illumina 0.979 at AF 0.08-0.25, 0.645 at AF >= 0.85; DeepSomatic-TO 0.985 / 0.930); the missed
high-AF somatic SNVs (LOH-like, AF ~1) score like germline. The SNV models have no population input, and most
high-AF germline candidates they learn from are common alleles that the PoN removes anyway.

## Code (`tmp/snv_sitefeat_20261103/code/`)

`git archive` of the repository at ca50ab0 (`code/BASE_COMMIT`), `indexed_gam_pipeline_v4` linked, plus the SNV
builder in `machine_learning/pansoma_net_v2/site_features.py`. The same builder is now in the repository
(`site_features.py build --kind SNV --dbsnp --g1000`, test `SnvBuilderTest`; `ris/site_features.sh` with KIND=SNV);
the copy differs from it only in the module docstring. INDEL files and their checks are unchanged.

- `build --kind SNV --dbsnp --g1000`: the same arrays and model input as the INDEL files (VERSION
  site-features-v1, five inputs; hp_len is 0 for SNVs). VCF allele = linear_vcf's (POS = pos0 + 1).
  `popaf_exact` = largest gnomAD / CoLoRSdb AF of the exact allele; `popaf_pos` = largest AF of any single-base ALT at
  the POS; `pon_tag` = scripts/filter_panel_of_normals.py's SNV rule over the four PoN VCFs (gnomAD and CoLoRSdb
  AF >= 1e-4, dbSNP with SAO != 2, any 1000G record), i.e. the PanelOfNormals tag vcfeval applies.
- `pon_rule(kind)`: the recorded and checked rule is the SNV one for SNV files.
- Check (`check/check_chr1.py`, HCC1395T_PacBio_d46 chr1 against the PoN-tagged linear VCF of
  multicancer6_PacBio_SNV normfix): 13,144 alleles, 8,894 tagged by both, 4,250 by neither, 0 disagreements.

## Premise (`check/premise.py` -> `check/premise.txt`; training chromosomes, PacBio sets)

Share of somatic labels among the labelled candidates by read AF, all / not in the population databases (pon_tag 0)
/ population alleles (pon_tag 1):

| set | AF 0.6-0.85 | AF >= 0.85 |
|---|---|---|
| HCC1395T | 6.2% / 49.0% / 0.27% | 4.1% / 52.5% / 0.11% |
| HG008T | 2.2% / 17.6% / 0.28% | 2.1% / 31.5% / 0.13% |
| COLO829T Fiber-seq | 7.7% / 42.4% / 0.37% | 5.4% / 44.7% / 0.21% |
| H1437T | 10.8% / 65.9% / 0.23% | 8.2% / 62.6% / 0.15% |

Without population information a high-AF candidate is somatic in 2-8% of the training examples; with it, half of
those absent from the databases are. On HG008T the non-population candidates at AF 0.4-0.6 are only 29% somatic
(the patient's private germline variants), which the features cannot separate.

## Runs

- `build_sitef.sh SET` -> `/scratch/jshen/data/snv_sitefeat_20261103/site_features/<SET>.SNV.npz` (general, -c 4).
- `train_sitef.sh PLAT` -> `pansoma_net_v2_runs/multicancer6_<PLAT>_SNV_sitef_ppr9`: the data, chromosomes, seed
  and recipe of `multicancer6_<PLAT>_SNV` (jobs/multicancer_20261016: six tumors, HCC1395T without chr1 / chr20,
  HG008T and COLO829T without chr1, 50.4 M, 20 epochs, pow:0.75, lr 1e-4, batch 1024 on one H100, `--norm dataset`)
  plus `--site-features`, `--best-by postpon`, `--call-recall 0.9`. Then predictions of chr1 (HCC1395T, HG008T,
  COLO829T) and HCC1395T chr20.
- `vcf_code.sh` (vcf.sh on `code/`) per test set: PASS = the checkpoint's call rule (post-PoN validation recall 0.9
  of the callable positives), PoN, rtg in the benchmark BED.
- Jobs: `jobs.tsv`.

## Evaluation (fixed before any result)

Against `multicancer6_<PLAT>_SNV` (normfix predictions) and the comparators (E8), chr1 SNVs after the PoN in the
benchmark regions, PASS at each model's own rule:
1. HCC1395T chr20 (held out of both models' training and validation; not used in the diagnosis above): F1, P, R.
2. chr1 of HCC1395T, HG008T and COLO829T: F1, P, R (HCC1395T chr1 was used in the diagnosis).
3. HCC1395T recall by read AF (the targeted effect; `recall_by_af.py`).
The change is one run against one run (seed 0); differences smaller than the seed range of comparable models
(0.01-0.08 F1) are not taken as an effect.

## Results: PacBio HiFi (2026-11-03; `compare.py`, `results/results_PacBio.md`)

Training: best epoch 6 by post-PoN validation F1 (0.872; epochs 6-20 all 0.866-0.872); PASS = p_somatic >= 0.875
(post-PoN validation recall 0.9 of the callable positives). The normfix model's PASS is p >= 0.622.

SNV F1 (P / R) after the PoN, benchmark BED, PASS at each caller's own rule:

| test | site features (ppr9) | normfix | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|
| HCC1395T chr20 | 0.794 (0.847 / 0.747) | 0.788 (0.804 / 0.773) | 0.835 (0.740 / 0.957) | 0.795 (0.778 / 0.812) |
| HCC1395T chr1 | 0.828 (0.895 / 0.769) | 0.825 (0.843 / 0.809) | 0.839 (0.766 / 0.927) | 0.814 (0.805 / 0.822) |
| HG008T chr1 | 0.659 (0.614 / 0.710) | 0.625 (0.586 / 0.670) | 0.557 (0.425 / 0.806) | 0.591 (0.455 / 0.844) |
| COLO829T chr1 | 0.777 (0.879 / 0.697) | 0.774 (0.880 / 0.691) | 0.632 (0.732 / 0.555) | 0.762 (0.779 / 0.745) |

By the pre-specified rule there is no effect: F1 changes by +0.003 to +0.034, inside the seed range. On HCC1395T the
model gains precision and loses recall (chr1 0.809 -> 0.769, chr20 0.773 -> 0.747), at every read-AF bin at PASS
(AF >= 0.85: chr1 0.801 -> 0.789, 251 truths; chr20 7 truths).

Diagnostic after the results (`matched_recall.py`; not pre-specified, HCC1395T chr1 was used in the diagnosis): the
ranking is better but the PASS rule is stricter. At the cutoff where the new model's recall equals normfix's PASS
recall (chr1 p >= 0.840, chr20 p >= 0.859), precision is 0.879 vs 0.843 (chr1) and 0.842 vs 0.804 (chr20), and chr1
recall at AF >= 0.85 is 0.884 vs 0.801 (AF 0.25-0.6: 0.935-0.937 vs 0.952). Precision at recall 0.8 (vcfeval curve):
HCC1395T chr1 0.884 vs 0.849, chr20 0.826 vs 0.792, HG008T 0.542 vs 0.504. Test-tuned best F1 (reference only):
0.847 vs 0.831 (chr1), 0.817 vs 0.799 (chr20), 0.673 vs 0.638 (HG008T), 0.796 vs 0.791 (COLO829T), the HCC1395T
optima at p ~0.79, near the post-PoN validation best-F1 cutoff (0.795) rather than the recall-0.9 one (0.875).
Truth SNVs without a tensor in the VCF (chr1 433 of 3,440, chr20 151) stay missed, as expected.

## Second PASS rule, registered 2026-11-03 21:57 before any Illumina or ONT test result

The Illumina model (multicancer6_Illumina_SNV_sitef_ppr9, job 467959) was at epoch 7 of 20 and had no test predictions.

- Rule 2: PASS = p_somatic >= the best epoch's post-PoN validation threshold at callable recall 0.95
  (val `somatic_at_recall_postpon_callable["0.95"]`; same checkpoint, no retraining). Rule 1 (ppr9, recall 0.9)
  stays the pre-specified result.
- Why: on PacBio HiFi the recall-0.9 cutoff (0.875) lies above the validation best-F1 cutoff (0.795), where the
  callable validation recall is already 0.967, and the HCC1395T test optimum is at p ~0.79. Chosen after the
  PacBio test results, so its PacBio numbers below are not evidence; Illumina (and ONT, if run) is the test.
- Evaluation: items 1-3 above for both rules (`compare.py` for rule 1; `at_cutoffs.py`, from the vcfeval curve of
  the same VCF, for rule 2), against the same models and comparators, with the same seed-range caveat.

PacBio HiFi under rule 2 (p >= 0.825; post hoc), F1 (P / R): HCC1395T chr20 0.813 (0.832 / 0.795), chr1 0.845
(0.872 / 0.820), HG008T 0.645 (0.549 / 0.780), COLO829T 0.792 (0.859 / 0.734); against normfix +0.025 / +0.020 /
+0.020 / +0.018, against rule 1 +0.019 / +0.017 / -0.014 / +0.015 (HG008T loses precision: private germline).
HCC1395T chr1 recall at AF >= 0.85 0.908 (normfix 0.801, DeepSomatic-TO 0.888, ClairS-TO 0.944), AF 0.6-0.85 0.874
(0.836). The validation best-F1 cutoff (0.795) gives within 0.005 F1 of rule 2 on all four tests.
