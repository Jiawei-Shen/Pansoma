# Fact sheet v3: R3.m1 (site features, configs), R3.m2 (PR-curve rules), R2.code (public state), 2026-10-30

Read-only check. Nothing tracked was edited, staged or pushed. Abbreviations:
`REPO` = /scratch/jshen/Github/Pansoma; `ML` = REPO/machine_learning/pansoma_net_v2 (develop working tree);
`CF` = REPO/tmp/indel_gap_20261020/code_final/machine_learning/pansoma_net_v2 (the frozen code the final INDEL runs used);
`RUNS` = /scratch/jshen/data/pansoma_net_v2_runs; `E8W` = /scratch/jshen/data/reviewer_experiments_20261013/E8_comparators;
`E8` = REPO/analysis/reviewer_response_20261011/experiments/E8_comparators; `R2C` = REPO/tmp/assessment_20261030/r2code.

---------------------------------------------------------------------------------------------------------------------

## 1. R3.m1: the five site features of the final INDEL models

### 1a. Where they're switched on (EXISTS)
- `REPO/tmp/indel_gap_20261020/final_20261023/train_final.sh`: `--site-features /scratch/jshen/data/indel_gap_20261020/site_features --best-by postpon`,
  no `--scalars`. Code = `CF` (`code_final/CHANGES.md`: copy of code_v3f, base commit 77cee0e plus INDEL-gap patches).
- All 18 final runs (`RUNS/{release,unseen}_s{0,1,2}_{Illumina,PacBio,ONT}_INDEL/args.json`): `"site_features": ".../indel_gap_20261020/site_features"`,
  `"scalars": false`, `"best_by": "postpon"`, `"norm": "dataset"`, seeds 0/1/2. Checked by script: no other argument differs except
  `tensors` / `set_test_chroms` / `seed` / `micro_batches`, and `resume` in 4 runs (release_s2 and unseen_s2, ONT and PacBio).
- Checkpoint `RUNS/release_s0_Illumina_INDEL/best.pth` (loaded on CPU): `config.site_features = 5`, `config.scalars = 0`;
  `site_features.version = "site-features-v1"`, `names = [log_popaf_exact, log_popaf_pos, pon_tag, hp_len, has_grch38]`,
  transform string as below; weights `scalar_mlp.0.weight (64, 5)`, `backbone.head.weight (3, 576)`.

### 1b. Definitions (EXISTS: `ML/site_features.py` docstring lines 1-55, `transform()`, `build_chrom()`, `_normalise()`, `_stream_db()`;
same `VERSION`, `SITE_FEATURES`, `TRANSFORM`, `transform()` in `CF/site_features.py` lines 40-48, 149-159; the files used were built by
`REPO/tmp/indel_gap_20261020/site_features/site_features.py`, documented in `/scratch/jshen/data/indel_gap_20261020/site_features/README.md`)

The allele every feature refers to: the tensor summary's `grch38` projection {pos0, ref, alt} (the ReferencePath projection
linear_vcf uses), REF checked against the FASTA, then `linear_vcf.Genome.left_align` + `Genome.vcf_allele`, i.e. the allele the
output VCF carries. Off-reference tensors (grch38 null) get raw value -1 (pon_tag -1, hp_len 0).
Inputs: GRCh38 FASTA `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`; databases
`/scratch/jshen/data/Pansoma/panel_of_normal_VCFs/af-only-gnomad.hg38.vcf.gz` and `.../CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz`
(npz `meta.dbs` of HG008T_ONT_d46.INDEL.npz). No read data, no dbSNP, no 1000G PoN.

| # | model input | raw array | definition / computation | source | model scaling | missing / off-reference |
|---|---|---|---|---|---|---|
| 1 | `log_popaf_exact` | `popaf_exact` | max over gnomAD and CoLoRSdb of the AF of the **exact** (POS, REF, ALT). DB records: `bcftools query -r <chrom> -e 'TYPE="snp"' -f '%POS\t%REF\t%ALT\t%AF\n'`, split per ALT, shared bases trimmed (suffix, then prefix, to an empty allele), left-aligned + padded with the same Genome functions; complex alleles (both sides left) dropped; AF '.' ignored; AF read as float32 | population DBs + reference | `1 + log10(clip(AF,0,1) + 1e-6) / 6` → 0 (absent) … 0.333 (AF 1e-4) … 1 (AF 1) | absent from both DBs → AF 0 → 0; off-reference raw -1 → clipped to 0 |
| 2 | `log_popaf_pos` | `popaf_pos` | max normalised AF of **any** pure INDEL allele of either DB whose left-aligned POS equals the tensor's VCF POS (any length/sequence; includes the exact allele) | population DBs + reference | same log map | none at POS → 0; off-reference → 0 |
| 3 | `pon_tag` | `pon_tag` | 1 if `scripts/filter_panel_of_normals.py`'s INDEL rule tags the allele: a gnomAD or CoLoRSdb record with the same POS and REF and the ALT among its ALTs (unnormalised, = `record_matches`), AF ≥ floor·(1-1e-6), floor 0.01 if hp_len ≥ 8 else 1e-4 (`INDEL_MIN_AF_HOMOPOLYMER`, `INDEL_MIN_AF`, `INDEL_HOMOPOLYMER` lines 32-34; script sha256 a9e6e42f…, unchanged since the build) | population DBs + reference (the population filter itself) | `pon_tag == 1` → 1, else 0 | not tagged → 0; off-reference (-1) → 0 |
| 4 | `hp_len` | `hp_len` | `filter_panel_of_normals.homopolymer_run(fasta, chrom, POS, REF, ALT)` (lines 41-62): if the inserted/deleted bases are one repeated base, the run of that base on the reference left + right of the anchor (200-bp window), at least the deleted length; 0 otherwise (incl. multi-base non-homopolymer INDELs) | reference FASTA | `min(max(hp,0), 20) / 20` | off-reference → 0 |
| 5 | `has_grch38` | derived | 1 if `pon_tag >= 0 and popaf_exact >= 0`, i.e. the tensor has a GRCh38 allele | tensor's graph→GRCh38 projection | 0/1 | 0 for off-reference |

- Not z-scored (docstring: "a fixed map to [0, 1], not z-scored"); final `np.clip(out, 0, 1)`; float32. (EXISTS `ML/site_features.py` `transform`)
- Missing feature file or a file that disagrees with the tensor index (chromosomes, counts, candidate-id sha1, off-reference
  consistency, VCF-position agreement ≥ 99 %, PoN-rule constants) → run refused; no imputation. (EXISTS `load()`)
- How they enter the network (EXISTS `ML/model.py` lines 203-264, `ML/data.py` lines 399-434): `TensorDataset` appends the 5 values to
  each sample's scalars; with `scalars=0` the model input is just these 5 columns → `scalar_mlp = Linear(5, 64) + GELU` → concatenated
  with the 512-d pooled ConvNeXt (2-2-6-2, 64-512) features → `head = Linear(576, 3)`. The 9 read-count "SCALARS" (coverage, AF …) are
  NOT used (`--scalars` off).
- Also used outside the network (EXISTS `CF/train.py` lines 82-86, 625-626; train_final.sh): `pon_tag` drives the checkpoint choice
  (`--best-by postpon`: best validation F1 with PoN-tagged tensors never called) and the post-filter validation recall-0.8 PASS threshold.

### 1c. Population prior: YES (explicit)
Features 1, 2 and 3 are population priors: allele frequencies from gnomAD (GATK af-only resource) and CoLoRSdb v1.1.0, and the population
filter's own INDEL tag (gnomAD / CoLoRSdb only). Features 4 and 5 are not (reference sequence, projection status).

### 1d. Consequences for other points (letter ver2, `.../Pansoma_Resubmit/Reviewer_Letter_Response_ver2.md`)
- R3.m2 line 748 says "DeepSomatic-TO receives population allele frequencies as an input channel"; nothing says Pansoma's final INDEL
  models also do. For fairness the Methods / R3.m2 should say so.
- R2.6 point 3 (line ~385) and R1.6 table ("Pansoma PASS, before population filtering", INDEL column): the final INDEL models' "without
  filter" numbers are not free of population information (the model sees gnomAD/CoLoRSdb AF and the filter's tag, and the checkpoint and
  threshold were chosen after the filter). Add a clause, e.g. "(the final INDEL models take population frequencies as input, so this row
  is not population-free)".

---------------------------------------------------------------------------------------------------------------------

## 2. R3.m1: configuration files and split manifests

| item | status | where |
|---|---|---|
| args.json (every CLI argument incl. seed, val_fraction 0.05, val_block_nodes 20000, set_test_chroms) for the 18 final INDEL runs | EXISTS, not in git | `RUNS/{release,unseen}_s{0,1,2}_*_INDEL/args.json` (+ best.pth, train.log, metrics.jsonl, dataset_stats.json) |
| args.json for every Pansoma run read by `E8/latest_models.json`, `purity_comparators.json`, `E17a_indel_strata/final_models.json` (54 RIS runs: SNV normbq 3 seeds × arms, Illumina INDEL normbq, long-read INDEL 7fc7396) | EXISTS (all 72 run dirs have args.json + best.pth) | `RUNS/ris_7fc7396_20261012/<run>/`; the 3 COLO829 Illumina d9 `*_sweep_e05/e06` (chr1-selected epochs) have **no train.log** |
| multicancer4/6 SNV, HapMap_Illumina_{SNV_full,INDEL}, COLO829T_fiberseq_MIX_SNV, COLO829T_ONT_SNV_joint_{control,readnorm}_tq | EXISTS (args.json + best.pth) | `RUNS/<run>/` |
| training job scripts of the final INDEL batch | EXISTS, git-ignored (`.gitignore:13 /tmp/`) | `REPO/tmp/indel_gap_20261020/final_20261023/{train_final.sh,PREREGISTRATION.md,driver.py,state.json}` |
| frozen code of the final INDEL batch | EXISTS, git-ignored | `REPO/tmp/indel_gap_20261020/code_final/` (base 77cee0e + patches; differs from develop's `ML/site_features.py`, `train.py`, `predict.py`, `data.py`, `graph_vcf.py`) |
| site-feature files (18 sets) | EXISTS, not in git | `/scratch/jshen/data/indel_gap_20261020/site_features/*.INDEL.npz` (+ README.md) |
| split recorded in the checkpoint | EXISTS | best.pth keys `chroms` {train, val, test, set_test}, `data` [{directory, labels, train count}], `args` |
| explicit list of validation blocks per model ("split manifest") | **MISSING** as a file. Deterministic: `CF/data.py` `_block_u` = sha1(f"{seed}:{chrom}:{block}") < 0.05, block = node_id // 20000 (lines 276-291), so it can be regenerated from args.json |
| per-set train/val counts | EXISTS (tracked, public) only for the 7fc7396-era sets / 42 runs: `analysis/.../E21_training_composition/set_split_counts.tsv`, `runs_training.tsv`; the normbq SNV models and the 18 final INDEL models are not in them |
| "Supplementary Table R3.m1" | **MISSING** (no file in `Pansoma_Resubmit/` or `analysis/`; `git ls-files` finds no args.json / config / manifest) |
| any args.json / split manifest in git (tracked) | **MISSING** (`git ls-files | grep -iE "args\.json|split|manifest"` → only E21 set_split_counts.tsv and points/R3.m1.md) |
| deposited anywhere public | **MISSING** |

Other R3.m1 table details from args.json that the letter's table omits: final INDEL models use `--norm dataset` (BQ z-scored per dataset),
`--non-fraction 1.0`, `--ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth`, `--non-reasons repeat_indel_without_truth`,
`--non-reasons-in-bed outside_confident_region`; 4 of 18 runs were resumed (`resume` set). These belong in Supplementary Table R3.m1.

---------------------------------------------------------------------------------------------------------------------

## 3. R3.m2: what the comparator scoring actually does

Source: `E8/score_comparator.py` (all comparator reports under `E8W/score/<caller>/<dataset>/<KIND>/report.json`, 2026-10-16),
`ML/vcfeval.py` (functions `vcfeval`, `subset`, `run_pon`, `describe`), `E8/latest_models.py` lines 141-176 (R1.6 matched precision),
`E8/collect.py` lines 120-124, 285-323. All EXISTS.

### 3a. Record preparation (both tools)
- `bcftools norm -m -any -f <labels FASTA> -c w -t chr1 -O v <vcf>`; only records of the kind (SNV: ALT length = REF; INDEL: other).
- Dropped: `if set(filt.split(";")) & NO_CALL_FILTERS: continue` with `NO_CALL_FILTERS = {"RefCall", "NoCall"}`; ALT `*`, `.`, `<…>`.
- One record per (POS, REF, ALT), highest QUAL kept; QUAL '.' → 0. FILTER copied unchanged.
- Counted in the main chr1 reports: DeepSomatic RefCall/NoCall dropped 22,651–983,567 per set; ClairS-TO wrote **no** RefCall/NoCall
  records (counter absent), although its header declares RefCall.

### 3b. Call sets (`score_comparator.evaluate`)
- native operating point `pass`: `subset(calls_all, calls["pass"], lambda f: f == {"PASS"})` → FILTER exactly PASS: every native filter applies.
- harmonized operating point `pass_pon`: our filter tags `calls.all` (`scripts/filter_panel_of_normals.py`, 4 PoNs, adds `PanelOfNormals`),
  then `subset(..., lambda f: f == {"PASS"})` → native PASS **and** not tagged by us.
- curve `all`: every kept record (any FILTER), no PoN; curve `all_pon`: `subset(calls.all.pon, ..., lambda f: "PanelOfNormals" not in f)`.
- rtg: `rtg vcfeval -b truth -c calls -t GRCh38.primary.sdf --squash-ploidy --sample ALT,ALT --bed-regions chr1 -f QUAL [-e BED] [--all-records]`;
  `--all-records` only for `all` / `all_pon`. **Score swept = the VCF QUAL column** (`-f QUAL`).

### 3c. FILTER values kept on the curve (`all_pon`), from the main reports' `counts.filters` and the VCF headers
- ClairS-TO 0.4.1/0.4.2: kept = PASS, **NonSomatic** ("Non-somatic variant tagged by panel of normals"; in every main VCF NonSomatic ⇔ an INFO
  `PoN_1..4` flag: gnomAD r2.1 AF≥0.001 allele, dbSNP 138 non-somatic allele, 1000G PoN position, CoLoRSdb AF≥0.001 position; no
  Verdict flags present), **LowQual, VariantCluster, StrandBias, ReadStartEnd, Realignment, LowSeqEntropy, LowAltBQ, MultiHap, NoAncestry**.
  Dropped: RefCall (none present).
- DeepSomatic 1.9.0: kept = PASS, **GERMLINE** ("Non somatic variants" = the model's germline class, not a database tag), **PON** (default PoN).
  Dropped: RefCall, NoCall. (LowQual declared, never present.)
- So: native population tags AND technical filters are both lifted on the curve; only our PoN applies. The letter's "Technical quality filters
  are kept in both operating points" is true for the two PASS operating points only, not for the curves.
- Size of the effect, SNV, `all_pon` composition (counted from `E8W/score/*/<set>/SNV/calls.all_pon.vcf.gz`):
  DeepSomatic HG008T Illumina / PacBio / ONT: PASS 2,425 / 1,470 / 2,297 + GERMLINE 924 / 274 / 222; COLO829T Illumina / fiberseq / ONT:
  PASS 4,779 / 1,853 / 2,915 + GERMLINE 1,327 / 1,447 / 603 (+10 PON). ClairS-TO HG008T: PASS 1,875 / 1,356 / 1,489 + technical-only
  1,618 / 172 / 597 + NonSomatic 192 / 30 / 41; COLO829T: PASS 4,335 / 2,334 / 2,233 + technical-only 5,453 / 365 / 563 + NonSomatic 493 / 91 / 46.
  E8 README lines 131-132, 208: DeepSomatic GERMLINE records carry high QUAL.
- Pansoma's own curve for comparison (`ML/vcfeval.py` docstring §3, `CURVE_FILTERS = {"PASS","LowQual"}`, line 288): PASS + LowQual records,
  not LowAF, not PanelOfNormals, QUAL = phred(P_SOMATIC); records below p 0.01 are not written (`graph_vcf --min-score 0.01`).

### 3d. Matched-precision recall (R1.6 point 2 numbers)
- `E8/latest_models.py` lines 153-169: comparator side = its **harmonized operating point** `pass_pon` (precision P_c, recall R_c); Pansoma side =
  `max(recall on Pansoma's post-filter curve "pon" with precision >= P_c)`, per seed, then min-max/mean. **The comparators' curves are not
  used.** (`collect.py` also computes the reverse — comparator `all_pon` curve recall at Pansoma's PASS+PoN precision — in `E8/tables.md`;
  the letter does not quote it.)
- PR-curve figure for the comparators: **MISSING** (curve data only in report.json `curve_roc`; not in `Reviewer_Letter_Response_ver2_figures/`).
- INDEL reports use the current PoN rule (truth tagged 418/635 HG008T, 81/131 COLO829T on chr1, all 12 main INDEL reports, 2026-10-16);
  the `score_comparator.py` docstring still describes the old rule (HP ≥ 7, 1e-3): stale comment only.
- HG008T ONT DeepSomatic scored VCF = `E8W/calls/deepsomatic/HG008T_ONT/output.vcf.gz` (`command.txt`: `--model_type=ONT_TUMOR_ONLY`), written
  18:46, scored 18:52/18:54 on 2026-10-16: consistent with the letter.

---------------------------------------------------------------------------------------------------------------------

## 4. R2.code: verifiable public state today (2026-10-30 ~13:30)

| check | result | status |
|---|---|---|
| visibility / default branch (`gh repo view`) | PUBLIC, default = `main` | EXISTS |
| `git ls-remote --symref origin` | HEAD → main = **e1c0447** (2026-08-13, "Explain inference emit output modes", v1 .dat/.idx pipeline); develop = **0bd80aa** (2026-10-30 13:13, analysis commit by a peer session) | EXISTS |
| develop vs main | origin/develop 269 commits ahead of origin/main; local develop = origin/develop (0 ahead/behind) | EXISTS |
| tags (`git tag`) / GitHub releases (`gh release list`) | none / none | **MISSING** |
| native.py self-test fallback (`R2C/patches/native_selftest_fallback.patch`) | in **no** commit (`git log --all -S"building the native decoder failed"` empty); HEAD, origin/develop, indel-recipe-20261023 and the working tree still raise on a failed static self-test (native.py lines 400, 405) | **MISSING** |
| code changes since the verified commit 44821f3 | none in indexed_gam_pipeline_v4/, machine_learning/, scripts/, tests/ (`git diff --stat 44821f3 0bd80aa -- …` empty): the verification still applies to develop's code | EXISTS |
| env files on develop | only retired v1 `machine_learning/pansoma_net/{Dockerfile,environment.yml,requirements.txt}`; no `env/` | **MISSING** |
| env files on main | v1 `environment.yml`, `docker/` (old pipeline) | EXISTS (wrong pipeline) |
| verified pinned envs + lock files | `R2C/env_files/pansoma-{cpu,gpu}.yml`, `*.conda-lock.txt`, `*.pip-lock.txt` (git-ignored tmp/) | EXISTS locally only |
| minimal example | `R2C/example/{run_example.sh,check_example.py,expected/expected.json,data_SHA256SUMS.txt}`, data 345 MB (git-ignored) | EXISTS locally only |
| README install / example section on develop | none (README grep: no install/conda/zenodo/example section) | **MISSING** |
| v2 checkpoints / example data deposited (Zenodo, DOI) | none; `pretrained_model/Pansoma_zenodo` (untracked) = v1 models; no DOI in tracked files | **MISSING** |
| response drafts on the public branch | origin/develop tracks `Reviewer_Response_Assessment.md`, `Reviewer_Response_Assessment_v2.md` and 32 files in `analysis/reviewer_response_20261011/points/` | EXISTS (confidentiality risk) |
| clean-env results (logs in `R2C/logs/`) | final one-command run 3:39.42 wall, max RSS 5,705,548 KiB (= 5.4 GiB / 5.8 GB, so "5.7 GB" is a unit slip); GPU run 2:20.31 on "NVIDIA H100 NVL MIG 1g.12gb" (a 12 GB MIG slice, not a full H100); CPU: v4 207 tests OK, ML 112 OK, viz 14 OK, vg 8 OK; GPU: ML 112 OK; `check_example`: "OK: matches the expected output" | EXISTS |
| example scope | models = `HG008_Illumina_{SNV,INDEL}_d9_normbq` (HPRC v1.1); does **not** exercise the final site-feature INDEL models, giraffe alignment, graph preparation, site-feature building, labelling or training | EXISTS (scope limit) |

---------------------------------------------------------------------------------------------------------------------

# Draft text

## A. R3.m1, Methods definition of the site features (108 words)

Final INDEL models also take five site features of the candidate's left-aligned GRCh38 allele: (1) its highest allele frequency in gnomAD
or CoLoRSdb v1.1.0 (exact allele; database records normalised alike); (2) the highest frequency of any database INDEL allele at the same
position; (3) whether the population filter tags it; (4) the reference homopolymer length of the inserted or deleted base; (5) whether a
GRCh38 allele exists. Frequencies are mapped to 1 + log10(AF + 10⁻⁶)/6 and length to min(L, 20)/20; absent alleles and off-reference
candidates receive 0. A linear layer (5→64, GELU) joins them to the pooled ConvNeXt features before the classifier. Features 1–3 are
population priors.

Table cell replacement (R3.m1 "networks" row): "… (final INDEL models also take five site features, three of them population allele
frequencies from gnomAD and CoLoRSdb; Methods) …". Add to R3.m2 "Population resources": "Pansoma's final INDEL models likewise take gnomAD
and CoLoRSdb allele frequencies as input (Methods)." Add to R2.6 point 3 / R1.6 INDEL rows: "the final INDEL models take population
frequencies as input, so their values before filtering are not free of population information."

## B. R3.m1, configuration and split sentence (honest)

Replace "Supplementary Table R3.m1 records, for every reported model, its configuration and its split manifest (training chromosomes,
validation blocks, test chromosome), which are the sole reference for the results." with:

"Supplementary Table R3.m1 lists, for every reported model, its training arguments, code version and data split: the training tumors and
chromosomes, the chromosomes held out for each tumor, and the validation fraction, block size and seed, from which the validation blocks
follow deterministically. The configuration files, training logs and validation-block lists will be deposited with the model checkpoints at
[Zenodo DOI]."

(To make it true: generate the table from `RUNS/*/args.json` + best.pth `chroms`/`data`; write per-model validation-block lists with
`data.block_split`; deposit args.json, train.log, metrics.jsonl and block lists with the checkpoints; include `--norm`, non/ignore reasons,
resume, GPU count.)

## C. R3.m2, replacement for "The precision–recall curves rank … Technical quality filters are kept in both operating points."

Native performance uses each tool's PASS records, so all of its own filters apply: population and panel-of-normals tags (ClairS-TO
NonSomatic; DeepSomatic-TO PON and its GERMLINE class) and technical filters (for example ClairS-TO LowQual, VariantCluster and StrandBias);
the harmonized operating point keeps those PASS records that our population filter does not tag. The precision–recall curves instead take
every record except RefCall and NoCall, whatever its FILTER, remove the alleles our filter tags, and rank the rest by the tool's QUAL (rtg
vcfeval --all-records -f QUAL), so along the curves no native filter, population or technical, is applied. Recall at matched precision is
Pansoma's highest recall on its own post-filter curve at a precision at least equal to that of each comparator's harmonized operating
point, compared with that operating point's recall; the comparators' curves are not used for it.

(Optional Supplementary Note sentence: "Records are split and left-aligned, one record per allele is kept (highest QUAL), and DeepSomatic-TO
GERMLINE records, which carry high QUAL, enter the curves when our filter does not tag them (222–924 HG008T SNV records next to 1,470–2,425
PASS records).")

## D. R2.code

### D1. Version that is true today
We thank the reviewer for testing Pansoma, and we apologise for the difficulty. We traced two causes: the repository's default branch still
held the original pipeline, and in a fresh environment the compiled extension of the current pipeline could fail its self-test without
falling back to a portable build. We have corrected the build and prepared pinned conda environments (CPU and GPU, with lock files) and a
minimal example: one command runs a sorted GAM of HG008-T Illumina reads over GRCh38 chr1:242.5–243.0 Mb (a chromosome held out from
training) through candidate discovery, tensor construction, the HPRC v1.1 SNV and INDEL models, VCF conversion and population filtering, and
compares the result with the expected output. In freshly created environments without user site-packages, the example finished in 3.7 min on
four CPU threads (5.4 GiB peak memory) and in 2.3 min on a 12 GB partition of an NVIDIA H100; the pipeline and model test suites passed, and
the example's tensors were byte-identical to those of the genome-wide run. These changes will be released on the repository's default
branch as [release tag], with the example data and model checkpoints at [Zenodo DOI]. Read alignment and graph preparation use vg and the
published HPRC graph files, as described in the README.

### D2. Final version (only after the steps below)
… "The fixed code, environments, example and README are in release [release tag] (commit [hash]) of https://github.com/Jiawei-Shen/Pansoma,
now the default branch, and the example data and model checkpoints are at [Zenodo DOI]. From a fresh clone of that release and the deposited
data, the example finished in …" (rest as D1).

### D3. Steps the authors must take for D2 to be true
1. Decide what to do with `Reviewer_Response_Assessment.md`, `Reviewer_Response_Assessment_v2.md` and `analysis/reviewer_response_20261011/points/`
   (32 files) on public `develop` before pointing reviewers at the repo (remove from the release branch; history keeps them unless rewritten).
2. Commit `R2C/patches/native_selftest_fallback.patch` to `indexed_gam_pipeline_v4/native.py` (+ docstring, v4 README §4 "native compile").
3. Commit `env/pansoma-{cpu,gpu}.yml`, `env/*.conda-lock.txt`, `env/*.pip-lock.txt` (from `R2C/env_files/`); omit Dockerfile / `.def` unless built
   and tested; mark `machine_learning/pansoma_net/{Dockerfile,environment.yml,requirements.txt}` as v1-only.
4. Commit `example/run_example.sh`, `check_example.py`, `expected/expected.json`, `data_SHA256SUMS.txt` (optionally `make_example.py`).
5. README: Installation + Minimal example sections (commands of `facts/r2code.md` §3, expected output, runtime), external tools (vg, gbz_graph_index),
   note that `/scratch` / `/wanglab` paths are the authors' cluster.
6. Merge `develop` into `main` or make `develop` the default branch (main is 269 commits behind and holds the v1 pipeline).
7. Create a tag + GitHub release (none exist) and record its commit hash.
8. Deposit on Zenodo: example data tarball (345 MB), the two example checkpoints, and the reported models (at least the 18 final INDEL
   checkpoints, their site-feature files or the build command with gnomAD/CoLoRSdb paths, and the SNV normbq checkpoints), plus args.json /
   train.log / validation-block lists (R3.m1); put the DOI in README and `run_example.sh`.
9. Re-run `run_example.sh` + `check_example.py` from a fresh clone of the tag with the data downloaded from the DOI; update the times if they change.
10. Fill [release tag], [hash] and [Zenodo DOI] in R2.code, R3.m1 and the opening summary (line 15 of ver2 also claims "per-model configuration
    and split manifests"), and remove the [待确认] note.
