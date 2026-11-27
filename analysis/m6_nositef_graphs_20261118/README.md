# Six-tumor models without the population site features as a model input, HPRC v1.1 and v2.1 (2026-11-18)

User 2026-11-18: run the latest pipeline, but without site features, for the six-tumor model on HPRC v2.1 and v1.1,
test and inference on HCC1395T, COLO829T and HG008T, three platforms, SNV and INDEL.

The GRCh38-only six-tumor models (`m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef`, old cluster) never had
the population site features as a model input; the v1.1 / v2.1 six-tumor models do. These 12 models give v1.1 and v2.1
the same treatment, so the three references differ in the graph only.

- Pipeline: repository b7fa8cf (pansoma_net_v2 code last changed in 1798490), `ris/recipes.sh` with `NOSITEF=1`: the
  SNV / INDEL recipe plus `--site-features-pon-only` (the feature files give only the PoN tags for the post-PoN epoch
  selection and the call rule; the five features are not fed to the head; predict reads no feature file).
- Design = the other six-tumor models: HCC1395T chr1 + chr20, HG008T chr1 and COLO829T chr1 left out of training and
  validation, HCC1954T / H1437T / H2009T train on every chromosome, 5 % of the node blocks validate, seed 0, one H100,
  batch 1024 (not DDP: the models they are compared with trained on one GPU).
- Sets = those of the six-tumor models with the features: SNV Illumina / PacBio HiFi on the AF 0.05 builds of the four
  CASTLE lines and COLO829T (`_snvaf05`; HG008T at 0.08), SNV ONT and every INDEL model on the AF 0.08 sets.
- Test: chr1 of the three test tumors (HCC1395T also chr20) predicted with `best.pth` (`test_<chrom>_<SET>/`), then
  `ris/vcf.sh --set` per tumor (`vcf_chr1_<SET>/`: graph VCF, GRCh38 VCF, population filter, rtg vcfeval in the
  benchmark BED and without it) at the checkpoint's call rule (SNV rule 2; INDEL ppr8, ONT pp).
- Runs: `$PANSOMA_RUNS/m6_nositef_graphs_20261118/multicancer6_<P>_<KIND>_{d9,d46}_nositef`.

Files: `submit.sh` (login node; submits every model whose sets are on storage2, skips the others and the ones that
exist), `train_predict.sh` (the GPU job), `to_ris.sh` (old cluster: copies the missing sets), `ready.py` (is a copied
set complete), `collect.py` -> `results.md` (the 18 SNV call sets and the comparison with the figure's six-tumor models). On the lab cluster
`collect.py` reads the copy of the runs pulled from RIS (2026-11-26), where the GRCh38-only and v1.1 ONT cells of the
comparison are no longer "pending".
`strata.py` -> `strata.md` (user 2026-11-24): Response Table 2.7 (HG008T chr1 SNV by the number of HG008-N germline
variants within 50 bp, E19's stratification) with the HPRC rows on these models, so no model of the table takes the
site features; the GRCh38-only rows are those of `E19_germline_density/strata_latest.md` (call sets on the old cluster).
`regions.py` -> `regions.md` (user 2026-11-24): the SNV columns of Response Tables 1.5 = 2.5 = 3.15 (ΔF1 by region,
before and after the filter) with the HPRC side on these models; GRCh38-only counts from
`ver4_20261123/regions/per_callset_strata.tsv`; point estimates, no bootstrap intervals (they need the GRCh38-only
records: `to_ris_grch38only.sh`, run on the old cluster, copies them). Strata BEDs downloaded to
`Dataset/GIAB_stratifications_v3.5` (md5 checked) and `Dataset/SMaHT_Regional_Categorization`.
`r26p2_nositef.py` -> `r26p2_nositef.md`, `.json` (user 2026-11-24): Reviewer 2 comment 6 point 2 (germline SNV
candidates the graph does not absorb: rejected by the model / removed by the filter) on these models, next to the
numbers of the models with the features (`ver4_20261123/r26p2`).

## Status

| model | state | LSF jobs (train; VCF HCC1395T, HG008T, COLO829T) |
|---|---|---|
| SNV v1.1 Illumina | done (best epoch 9) | 929947; 929948, 929949, 929950 |
| SNV v1.1 PacBio HiFi | done (best epoch 14) | 929951; 929952, 929953, 929954 |
| SNV v1.1 ONT | done (best epoch 9) | 929955; 929956, 929957, 929958 |
| SNV v2.1 Illumina | done (best epoch 8) | 932784; 932785, 932786, 932787 |
| SNV v2.1 PacBio HiFi | done (best epoch 14) | 932788; 932789, 932790, 932791 |
| SNV v2.1 ONT | done (best epoch 10) | 932792; 932793, 932794, 932795 |
| INDEL v1.1, v2.1 (6 models) | not run on RIS (no INDEL tensors of the four CASTLE lines on storage2); submitted on the lab cluster 2026-11-24, section below | |

The 18 v1.1 SNV sets were checked before submission: every shard of each manifest is on storage2 with the manifest's
size, label manifests and site-feature files present. The v2.1 SNV sets missing on 2026-11-18 morning
(`{HCC1954T,H1437T,H2009T}_{Illumina,PacBio}_d46_snvaf05`, `{HCC1954T,H1437T,H2009T}_ONT_d46`, eight SNV site-feature
files that were links to /scratch) were copied from the old cluster by the user on 2026-11-18 (12:30-14:30) and passed
the same check (`ready.py`) before their models were submitted; the v2.1 GPU jobs ask for 5 slots (`SLOTS=5`).

SNV v1.1 submitted 2026-11-18 12:05, v2.1 14:36; all six trainings, 18 chr1 VCF / vcfeval jobs and the HCC1395T chr20
predictions done by 2026-11-23 (results.md). The INDEL models were not submitted: the CASTLE lines' INDEL sets were
not copied to storage2.

## INDEL models on the lab cluster (user 2026-11-24)

User 2026-11-24: train and test the six-tumor INDEL models without site features on HPRC v1.1 and v2.1, three platforms.
Run on tequila, where every INDEL set is (`train_indel_lab.sh`, `submit_indel_lab.sh`; job ids in `jobs_indel_lab.tsv`,
logs in `tmp/m6_nositef_indel_20261124/logs`).

- Recipe and code = the GRCh38-only six-tumor INDEL model they pair with
  (`m4m6_graphs_20261110/multicancer6_<P>_INDEL_linear_nositef`): the final INDEL recipe with `code_final`
  (`tmp/indel_gap_20261020`; 8.7 M, 12 epochs, batch 1024, lr 1e-4, sqrt, wd 0.05, rule (a) negatives, `--best-by postpon`,
  `--micro-batches 4`, seed 0) and `--pon-features` in place of `--site-features` (PoN tags for the epoch selection and
  the call rule only; the features are not a model input; predict reads no feature file).
- Sets = those of the six-tumor model with the features on the same graph: v1.1 `multicancer6_<P>_INDEL_d9`
  (`tmp/m4m6_graphs_20261110/sets.sh`), v2.1 `release_s0_<P>_INDEL` (`<tumor>_<P>_d46`, COLO829T PacBio =
  `COLO829T_fiberseq_d46`). Held out of training and validation: chr1 of every tumor and HCC1395T chr20.
- Test: chr1 of HCC1395T, HG008T, COLO829T (and HCC1395T chr20) predicted with `best.pth`, then `vcf_v3.sh` per tumor:
  `vcf_ppr8_<chrom>_<tumor>` (PASS ppr8) and, for ONT, `vcf_pp_<chrom>_<tumor>` (PASS pp, the ONT rule).
- Runs: `/scratch/jshen/data/pansoma_net_v2_runs/m6_nositef_graphs_20261118/multicancer6_<P>_INDEL_{d9,d46}_nositef`.
- `collect_indel_lab.py` -> `results_indel.md`: chr1 before / after the population filter, and PASS after the filter
  next to GRCh38-only without the features and v1.1 / v2.1 with them ("pending" until a call set exists).

| model | training job (GPU, time limit) | training finished (best epoch) |
|---|---|---|
| INDEL v2.1 ONT | 484620 (H100, 2 d) | 2026-11-25 01:14 (5) |
| INDEL v1.1 ONT | 484629 (24 GB slice, 2 d) | 2026-11-25 14:38 (6) |
| INDEL v2.1 PacBio HiFi | 484638 (24 GB slice, 30 h) | 2026-11-25 04:26 (6) |
| INDEL v1.1 PacBio HiFi | 484643 (H100, 30 h) | 2026-11-24 21:29 (7) |
| INDEL v2.1 Illumina | 484648 (H100, 14 h; moved to a 24 GB slice 2026-11-24 13:53, `scontrol update TresPerNode`, so the H100s go to the HCC1395T path-count SNV models) | 2026-11-24 23:03 (3) |
| INDEL v1.1 Illumina | 484653 (24 GB slice, 14 h) | 2026-11-24 14:28 (6) |

All six trainings and their VCF jobs completed; `results_indel.md` has every call set (written 2026-11-26).

Each training asks for 4 CPUs (tequila had 16 idle; the paired trainings used 1.1-1.5 cores). A training cut by its
time limit continues from `last.pth` when `submit_indel_lab.sh` is run again (`ONLY="<P>:<graph>"`).
