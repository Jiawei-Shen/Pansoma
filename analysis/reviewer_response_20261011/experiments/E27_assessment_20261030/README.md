# E27: checks for the response letter ver2 / ver3 (2026-10-30 to 11-01)

Written while the letter was rewritten from the two assessments (`Reviewer_Response_Assessment.md`,
`Reviewer_Response_Assessment_v2.md`) and the ver3 evaluation (`Reviewer_Response_Evaluation_ver3.md`). The letters
(`Reviewer_Letter_Response_ver2.md`, `_ver3.md` and their figure folders) are outside the repository, in
`Pansoma_papers/ver1/Jiawei-Smaht-Paper/Pansoma_Resubmit/`. The scripts ran in `tmp/assessment_20261030/<subdir>/`
(git-ignored), where the job logs, intermediate data and code copies stay; this directory keeps the scripts, the small
result tables and the fact sheets. No model was trained; GPU inference only for `ablation/q6/infer` (tequila).

**Fact sheets.** `facts_ver2/` (ablation, denominators, crosstumor, misc, r2code; for ver2) and `facts_ver3/` (runid,
labels, methods; for ver3) give every number the letters use with its source path, tagged EXISTS / COMPUTED NOW /
MISSING. They contain the draft text of that day; where the letter changed later, the letter is current.

| directory | letter points | contents |
|---|---|---|
| `ablation/` | R3.1, R3.9, R1.1, R1.2 | `main_table.*`, `collect_main.py`: the three-reference chr1 SNV table, validation-selected checkpoints only (`--norm train`, the only setting with all three references x three seeds); `q3/`: overlap of the validation blocks of the three graphs and of the seeds in GRCh38 (Jaccard 0.014-0.034 between graphs); `q4/`: 2,000-resample 1-Mb block bootstrap of the 54 seed-matched pairs; `q5_ch6.*`: channel-6 ablation as seed-0 pairs with the same GPU count (Table R1.1); `q9/`, `q9b/`: COLO829T Illumina v1.1 at the validation-selected epoch 8 (extreme regions, strata) |
| `ablation/q6/` | R1.1, R1.m2 | Cross-node context at inference (Supplementary Table R1.1): `offnode.py` and `offnode_vs_repo.patch` (`--mask-off-node` in train / predict / data / chunks, against `machine_learning/pansoma_net_v2` of 90d7c0f) zero every cell outside the candidate's own node visit; `side/` builds the per-tensor node spans; `realcheck/` checks the own-node rule against `--debug-rows` builds; `infer/results.md`: the 12 seed-0 HPRC SNV models of Table R3.1 with and without the mask on the same chr1 candidates (recall falls in 12/12, somatic AP in 11/12, PASS F1 in 8/12). `retrain/train_tequila.sh` was not used; the retraining is E28. |
| `e9c_ep8/` | R3.3 | E9c / E9d / E9i rerun in an isolated copy with the validation-selected epoch-8 COLO829T Illumina v1.1 run instead of the exploratory epoch-6 one (`tables_epoch8.md`, `summary_epoch8.json`, `ranges.py`); the only changed letter range: calls PASS on all three platforms 37-81 -> 32-81. |
| `denominators/` | R1.4, R1.5, R2.6, R1.8 | Absorption denominators (Table R1.4), label status of the truth alleles per reference, HapMap Mutect2 on the full truth, the R2.6 filter table. |
| `crosstumor/` | R2.4, R1.7, R3.10 | Dilution recall by expected VAF, E15b precision by allele fraction, compact E17a strata, final INDEL model counts and P / R / F1. |
| `main_bootstrap/` | R3.9 | Pairs for the paired bootstrap of the main comparisons. |
| `misc/` | R3.2, R3.7 | False positives and germline alleles split by graph representation (`r32_*`; the ">= 46 of 463 haplotypes" proxy), HG008-T tumor-assembly support (`r37_*`). |
| `v3/candidates/` | R1.2, R3.1 | `candidate_retention.py`: matched-normal germline alleles and somatic truth SNVs that become a candidate on the GRCh38-only, HPRC v1.1 and v2.1 graphs (chr1-22, all truth records; `truth_recall.json` of each tensor set; Supplementary Table R1.2). |
| `v3/region_tables/` | R1.3, R3.3 | `make_tables.py`: SNV calling by region (E16 SMaHT strata, E10b GIAB v3.5 strata, paired intervals; INDEL in homopolymers >= 7 bp), long-read alignment by SMaHT region (E11a G sample) and graph-coordinate calls per run by validation class (epoch-8 E9c / E9d / E9i); the letter's Supplementary Tables R1.3a, R1.3b, R3.3. |
| `v3/labels/` | R1.8, R1.5, R3.2 | Off-path validation negatives (`offref_label_counts`, `offref_val_effect`: share of validation negatives, their false positives at the selected checkpoints), exact-germline / partial-somatic conflicts, the HapMap conditional mask breakdown (1,701 SNVs / 6,482 INDELs), the CHM13 exception of the proxy, the rtg region rule. |
