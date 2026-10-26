# E8b / E8c: literature and resource facts for the NBT response

Status (2026-10-13, collected 22:54):

- DONE: literature facts (sections 3.1-3.8); local checks 4.1-4.4 (login node + Slurm 387259); section 4.5
  (somatic truth alleles in the four resources) for all six somatic truth sets: HG008T, COLO829T, HCC1395, HCC1954,
  H1437, H2009 (Slurm array 387253 tasks 0-5, all COMPLETED, exit 0; collected with `bash collect.sh`).
- PENDING (held): the HapMap-positives contrast row of 4.5 (chr20; array task 387253_6). It is held (`JobHeldUser`),
  with the other reviewer jobs, after the user's complaint about the lab-cluster footprint; it stays held until the
  user agrees a cap. When it has run: `bash collect.sh` and add the row to 4.5.
  - Update 2026-10-15: 387253_6 COMPLETED at 13:12. Raw output
    `/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts/somatic_truth_in_resources.HapMap_positives_chr20.tsv`:
    SNV n 129,857, Pansoma rule 125,026 (96.3 %), gnomAD AF >= 1e-2 93,049, dbSNP SAO=2 0; INDEL n 41,174, Pansoma
    rule 16,182 (39.3 %), SAO=2 0. `bash collect.sh` has not been re-run, so the 4.5 table does not have the row yet.
- PENDING (user decisions, not computations): HG008 truth v0.2 vs v0.3 (G1, G2); drop the HapMap benchmark or use the
  existing donor-free graph (H10); how to treat the H1437 CASTLE truth (4.5); the ClairS-TO held-out chromosome
  (C4) can only be settled by HKU-BAL, because the comparators are not re-run (user decision 2026-10-13).
- Key result: no somatic truth allele of any set is a dbSNP-somatic (SAO=2) record; Pansoma's PoN rule tags
  2.8-10.2 % of SNV truth alleles on five sets but 30.7 % on H1437, where 27.8 % of the "somatic" SNVs are
  gnomAD alleles with AF >= 1 %; for INDELs, ClairS-TO's default rule tags more truth alleles than Pansoma's in every
  set (e.g. HG008T 87.1 % vs 61.8 %).

Reviewer points served: R1.m1, R3.m2, R2.6, R3.m3, R1.5, R2.1, R2.2 (each row says which one).

## 1. Question

Which facts about the comparators (DeepSomatic, ClairS-TO), the four population / panel-of-normals (PoN) resources,
the two HPRC graphs, and the four truth sources (GIAB HG008, SMaHT COLO829, SEQC2 HCC1395, SMaHT HapMap
mixture) can the response letter state, with a primary source for each? Where the drafts
(`../../response_to_reviewers_detailed.md`, `../../response_to_reviewers.md`) say `[待核实]`, what is the verified value?

## 2. Method

- Primary sources only: journal articles (Nature Biotechnology, Nature Communications, PMC full text of bioRxiv
  preprints), the tools' own GitHub repositories (README, release notes, code at a release tag), consortium READMEs on
  GitHub / NCBI FTP / Zenodo, and NCBI SRA run metadata. Nature HTML was read with a cookie-enabled `curl`; GitHub
  through `gh api`; bioRxiv metadata through `api.biorxiv.org` (bioRxiv full text was rate-limited, HTTP 429/1015).
- Local files were read only (headers, record counts); nothing was downloaded except documentation and three small
  `model.example_info.json` files of DeepSomatic 1.10.0 (a few hundred bytes each) into the session scratchpad.
- Four small local computations (section 4) check what the literature cannot: the gnomAD release behind the GATK
  file, dbSNP's somatic flags, CoLoRSdb's minimum allele count, whether the benchmark donors' private alleles are in
  the resources, and how many somatic truth alleles each resource / rule would tag.
- "[unverified]" marks anything that could not be checked against a primary source; "[inference]" marks a conclusion
  drawn from verified numbers.

Local data (read only):

- Resources: `/scratch/jshen/data/Pansoma/panel_of_normal_VCFs/` `af-only-gnomad.hg38.vcf.gz`,
  `Homo_sapiens_assembly38.dbsnp138.vcf.gz`, `1000g_pon.hg38.vcf.gz`,
  `CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz` (the files of the pipeline's PoN step).
- Somatic truth sets (4.5): HG008T `/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz`;
  COLO829T `/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz`;
  HCC1395 `/scratch/jshen/data/CASTLE/HCC1395/truth/HCC1395_somatic_SEQC2_v1.2.1.vcf.gz`;
  HCC1954 / H1437 / H2009 `/scratch/jshen/data/CASTLE/<LINE>/truth/<LINE>_somatic_CASTLE.vcf.gz` (DeepSomatic's
  training truth, somatic_only PASS); HapMap contrast `/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HapMap_somatic_snv_indel.vcf.gz` (chr20).
- Germline alleles (4.4): SMaHT HAPMAP6 v1.6 SNV benchmark (`/scratch/jshen/data/HapMap/HapMap_trueset/`), HapMap_truth
  HG005 germline, dipcall of the HG008-N and COLO829BL assemblies, CASTLE DeepVariant-hybrid germline VCFs (QUAL >= 10);
  exact paths in `membership_private_alleles.py`.

## 3. Facts table

Columns: claim | verified value | source | exact quote or location | used by.

### 3.1 DeepSomatic (Park, Cook, Chang et al., Nat Biotechnol 44, 1569-1578 (2026); published online 2025-10-16)

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| D1 | Citation | Nat Biotechnol 44, 1569-1578 (2026), doi:10.1038/s41587-025-02839-x; received 2024-08-21, published 2025-10-16; preprint bioRxiv 10.1101/2024.08.16.608331 (v1 2024-08-19) | https://www.nature.com/articles/s41587-025-02839-x ; https://api.biorxiv.org/details/biorxiv/10.1101/2024.08.16.608331 | "Nat Biotechnol 44, 1569–1578 (2026) ... Received: 21 August 2024 Accepted: 01 September 2025 Published: 16 October 2025" | R3.m2 |
| D2 | Train / tune / test chromosomes, all models | chr2-chr20 train, chr21-chr22 tune, chr1 held out, for every model in the paper | Nature article, Results "DeepSomatic overview" | "For training, using reads aligned to GRCh38, we use variants from chromosome 2 to chromosome 20 for training, and tune with chromosomes 21 and 22, while holding out chromosome 1 for validation. In all the models we describe in this paper, we consistently use this training scheme of holding out chromosome 1 for validation." | R1.m1, R3.m2 |
| D3 | Single-cancer models | HCC1395 only (SEQC2 truth, 41,072 SNV+INDEL); Illumina = SEQC2 reads, PacBio = PacBio public data, ONT = new R10.4 | same | "We used this set of 41,072 somatic SNV and insertion and deletion (indel) variants to train three initial sequence technology-specific models of DeepSomatic: single-cancer Illumina, single-cancer PacBio and single-cancer ONT. For sequencing data, we used existing Illumina reads from SEQC2 and PacBio HiFi sequencing data from PacBio ..., and additionally generated ONT v.R10.4 chemistry sequencing data" | R1.m1 |
| D4 | Multi-cancer models: lines in training | HCC1395, HCC1937, HCC1954, H1437, H2009; Hs578 held out completely | same, "Sequencing six tumor–normal cell lines" | "We used five of these lines (HCC1395, HCC1937, HCC1954, H1437 and H2009) to train new DeepSomatic models, which we call multi-cancer models ... We held out the Hs578 sample as a control for testing." | R1.m1, R1.7 |
| D5 | CASTLE = six pairs | HCC1395, HCC1937, HCC1954, Hs578 (breast), H1437, H2009 (lung) | same; https://github.com/CASTLE-Panel/castle README | "we generated sequencing data for a set of four breast cancer cell lines: HCC1395, HCC1937, HCC1954 and Hs578, and two lung cancer cell lines, H1437 and H2009" | R1.m1 |
| D6 | HG008 in DeepSomatic | Not trained on; used as a held-out evaluation line (whole genome) | Nature article, multi-cancer evaluation | "(4) held-out cell lines, evaluating two matching tumor–normal cancer cell lines, Hs578, which is part of CASTLE ..., and HG008, which has been sequenced by the GIAB effort independently. We evaluate these held-out cell lines across the whole genome instead of only chromosome 1" | R1.m1, R3.m2 |
| D7 | COLO829 in DeepSomatic | Not named in the article text (0 occurrences of "COLO" in the HTML full text); CASTLE says its COLO829 data are "not a part of the CASTLE cell lines" | Nature HTML (searched); CASTLE README "Release 1" | "Also includes additional COLO829 data (with ONT R9 and Illumina) that used in the manuscripts, but not a part of the CASTLE cell lines." Training use of COLO829 by DeepSomatic: none stated [unverified beyond the text search] | R3.m2 |
| D8 | Truth of HCC1954 / H1437 / H2009 used here is DeepSomatic's own benchmark, and the paper calls it biased | Our CASTLE truth = `DeepSomatic_HCC1395-model_benchmark` (made by DeepSomatic single-cancer models) | Nature article; `/scratch/jshen/data/CASTLE/README.txt` | "The chromosome 1 test set evaluation directly assesses against the trained single-cancer models; it is strongly biased toward the single-cancer trained models because the set was derived directly from it, and is biased toward DeepSomatic in general." | R1.m1, R3.m2, R2.2 |
| D9 | Purity / contamination training data | Multi-cancer models include tumor-purity and tumor-in-normal titrations; HCC1395 used SEQC2 titration BAMs | Nature article, Results + Methods "Generating titration bam files" | "datasets with various levels of tumor purity and normal purity were generated ... and included to train the DeepSomatic multi-cancer models." / "For the HCC1395 cell line, titration bams published by the SEQC2 consortium were used for training." | R1.m1, R2.4 |
| D10 | Tumor-only models in the paper | Illumina, PacBio, ONT tumor-only trained on HCC1395 (SEQC2 truth), chr1 held out | Nature article, "Extension of DeepSomatic" | "three tumor-only DeepSomatic models, Illumina tumor-only, PacBio tumor-only and ONT tumor-only, were trained using the HCC1395 benchmark from SEQC2 using a panel of normals (Methods), with chromosome 1 held out." | R3.m2 |
| D11 | FFPE and WES models | FFPE_WGS / FFPE_WES trained on SEQC2 HCC1395 FFPE samples, chr1 excluded; WES on SEQC2 samples | Nature article | "These models were trained on FFPE samples of the HCC1395 cell line published by SEQC2, and evaluated against the SEQC2 benchmarking set of variants on chromosome 1, which was excluded from training" / "The WES model was also trained on samples published by SEQC2" | R3.m2 |
| D12 | Tumor-only population input (make_examples) | Illumina: DeepVariant 1000 Genomes calls, AF > 0.05; PacBio / ONT: CoLoRSdb, AF > 0.05, as an extra tensor channel | Nature Methods "DeepSomatic tumor-only variant calling method"; code `run_deepsomatic.py` r1.9; `model.example_info.json` 1.10.0 | "For Illumina data, we use DeepVariant 1KG variant calls with population variants observed above 0.05. For PacBio and ONT, we use the Consortium of Long Read Sequencing Database data similarly filtered for variants observed with allele frequency higher than 0.05." 1.10.0 files: `population_vcfs` = `AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz` (WGS_TUMOR_ONLY), `AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz` (PACBIO_ / ONT_TUMOR_ONLY) | R3.m2, R2.6 |
| D13 | Tumor-only PoN filter (postprocess) | Short reads: dbSNP + gnomAD + 1000 Genomes; long reads: CoLoRSdb + dbSNP + gnomAD. In the released tool it is OFF unless `--use_default_pon_filtering=true` (flag default False); files `PON_dbsnp138_gnomad_ILMN1000g_pon.vcf.gz` (WGS/WES/FFPE TO) and `PON_dbsnp138_gnomad_PB1000g_pon.vcf.gz` (PacBio/ONT TO) | Nature Methods; https://github.com/google/deepvariant/blob/r1.10/scripts/run_deepsomatic.py (lines 282-299, 594-617); Dockerfile.deepsomatic r1.10 | "For short reads, a panel of normals is done using VCFs available from dbSNP, gnomAD and 1000 Genomes project. For long reads, we use the Consortium of Long Read Sequencing Database VCF with dbSNP and gnomAD." Code: `_USE_DEFAULT_PON_FILTERING = flags.DEFINE_boolean('use_default_pon_filtering', False, ...)`; docs: "By using `--use_default_pon_filtering=true` the somatic variants will be filtered using the default PON vcf that contains variant calls from dbSNP, gnomAD and 1000 genomes." | R3.m2, R2.6 |
| D14 | Docs: chr1 is the held-out chromosome for the released tumor-only models | chr1 | https://github.com/google/deepsomatic/blob/r1.10/docs/deepsomatic-case-study-wgs-tumor-only.md (same text in pacbio / ont tumor-only case studies) | "We run the analysis on `chr1` that we hold out during training." | R1.m1, R3.m2 |
| D15 | Release history of the models | 1.7.0 (2024-08-27): first tumor-only (Illumina WGS, PacBio, ONT), FFPE, WES, ONT T/N; 1.8.0 (2024-12-04): no model update (WGS T/N updated 2025-04-02); 1.9.0 (2025-05-13): FFPE tumor-only models, WGS and WGS_TUMOR_ONLY retrained with all manuscript data; 1.10.0 (2026-03-05): flags read from `model.example_info.json`, small models for tumor-normal modes | https://github.com/google/deepsomatic/releases | v1.9.0: "The `WGS` and `WGS_TUMOR_ONLY` models have been retrained with all datasets described in the manuscript, tumor-in-normal and normal contamination datasets." v1.10.0: "`run_deepsomatic` now reads parameters from `model.example_info.json` files" / "Introduced small models for tumor-normal modes in DeepSomatic improving the runtime between 12% to 40%." No 1.10.0 note mentions retraining. | R3.m2 |
| D16 | Which data the released PACBIO_ / ONT_TUMOR_ONLY models (1.10.0) were trained on | Not stated in any release note; the paper's tumor-only models used HCC1395 only | release notes, paper | [unverified]: whether the 1.9/1.10 long-read tumor-only checkpoints include the other CASTLE lines | R1.m1, R3.m2 |
| D18 | CASTLE alignments and depths (relevant for running comparators on the same BAMs) | Illumina BWA-MEM2, PacBio pbmm2, ONT minimap2; PacBio about 60x T / 60x N, ONT about 90x T / 30x N, Illumina downsampled to 100x T / 50x N | Nature article, "Sequencing six tumor–normal cell lines" | "Illumina reads were aligned using BWA-MEM2, PacBio HiFi reads were aligned with pbmm2 and ONT reads were aligned with minimap2. The PacBio alignments have about 60 times tumor and 60 times normal average coverage, and the ONT alignments have about 90 times tumor and 30 times normal coverage ... The Illumina alignments were downsampled to 100 times tumor and 50 times normal average coverage" | R3.m2 |
| D17 | Model names | WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES, WGS_TUMOR_ONLY, WES_TUMOR_ONLY, PACBIO_TUMOR_ONLY, ONT_TUMOR_ONLY, FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY (+ hybrid PacBio-Illumina) | Dockerfile.deepsomatic r1.10; README r1.10 | Dockerfile: "wgs wes pacbio ont ffpe_wgs ffpe_wes wgs_tumor_only wes_tumor_only pacbio_tumor_only ont_tumor_only ffpe_wgs_tumor_only ffpe_wes_tumor_only" | R3.m2 |

### 3.2 ClairS-TO (Chen, Zheng, Su et al., Nat Commun 16, 9630 (2025); GitHub v0.5.1)

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| C1 | Citation, version | Nat Commun 16, 9630 (2025), doi:10.1038/s41467-025-64547-z, published 2025-10-31; latest release v0.5.1 (2026-09-09); v0.4.0 (2025-03-10) "is the version used for the preprint" | https://www.nature.com/articles/s41467-025-64547-z ; https://github.com/HKU-BAL/ClairS-TO (README, releases) | "*v0.5.1 (Sep 9, 2026)* : Added ONT Dorado HAC 5kHz v6 model (`-p ont_r10_dorado_hac_5khz_v6`)." | R3.m2 |
| C2 | SS (synthetic) training data | ONT: GIAB HG002 + HG001 (EPI2ME); PacBio Revio: HG003 + HG004; Illumina: HG003 + HG004 | Nat Commun article | "ClairS-TO SS was trained on synthetic datasets derived from GIAB HG002 and HG001 from the ONT EPI2ME Labs." / "The ClairS-TO SS PacBio model was trained on a synthetic dataset derived from GIAB HG003 and HG004, sequenced on the PacBio Revio system." / "For training, we used GIAB HG003 and HG004 for the SS model" (Illumina) | R3.m2 |
| C3 | SSRS (real tumors) | Fine-tuned on HCC1937, HCC1954, H1437, H2009 (CASTLE, with DeepSomatic multi-cancer-model truth); HCC1395 excluded entirely; COLO829 not trained on | Nat Commun; README; https://github.com/HKU-BAL/ClairS-series-model-training (clairs-to SSRS template lists `*_DeepSomatic_multicancer-model.vcf.gz`) | "ClairS-TO SSRS was augmented from the SS model using four real cancer cell lines—HCC1937, HCC1954, H1437, and H2009—from Park et al." / "left out the entire HCC1395 cell line for benchmarking purposes" | R1.m1, R3.m2 |
| C4 | Held-out chromosome: sources disagree | README: chr20 excluded from all models; paper: chr1 excluded; SS training template: chr20 excluded | README "Pre-trained Models"; Nat Commun Discussion; `clairs-to_ss_model_training_workflow.sh` line 57 | README: "All models were trained with chr20 excluded (including only chr1-19, 21, 22)." Paper: "In addition to excluding chromosome 1 from training, ClairS-TO uses only four real tumor samples". Template: "# Array of chromosomes ... Exclude chromosome 20 for full training." Which applies to the released ssrs checkpoints: [unverified]. If chr20, ClairS-TO ssrs has seen chr1 of HCC1954 / H1437 / H2009. | R1.m1, R3.m2 |
| C5 | Benchmarks in the paper | COLO829 (NYGC truth: 42,993 SNV, 985 INDEL) and HCC1395 (SEQC2 HighConf+MedConf in HC regions: 39,447 SNV, 1,602 INDEL); no HG008 (0 occurrences in full text) | Nat Commun | "We used two cancer cell lines COLO829 and HCC1395 for benchmarking. The truth somatic variants of COLO829 are provided by NYGC, with a total of 42,993 SNVs and 985 Indels." | R3.m2 |
| C6 | Default PoNs | gnomAD (GATK af-only-gnomad.hg38, called "gnomAD r2.1"), AF >= 0.001; dbSNP 138 non-somatic; GATK 1000G PoN, all sites; CoLoRSdb v1.1.0, AF >= 0.001 | README "Tagging non-somatic variant using panels of normals"; Nat Commun Methods | Table: PoN1 source `gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz`, 268,225,276 total, "Sites with AF ≥ 0.001", 35,551,905 used; PoN2 `Homo_sapiens_assembly38.dbsnp138.vcf`, 60,691,395, "Non-Somatic sites", 60,683,019; PoN3 `1000g_pon.hg38.vcf.gz`, 2,609,566, "All sites"; PoN4 CoLoRSdb v1.1.0, 49,550,902, AF ≥ 0.001, 41,175,834 | R3.m2, R2.6 |
| C7 | PoN matching rule | Exact allele for gnomAD and dbSNP; position only for 1000G PoN and CoLoRSdb (default `True,True,False,False`) | Nat Commun Methods; README `--panel_of_normals_require_allele_matching` | "We required exact allele matching for a variant to be tagged as non-somatic for gnomAD and dbSNP, and positional matching for 1000 G PoN and CoLoRSdb." | R3.m2, R2.6 |
| C8 | ClairS-TO "non-somatic" dbSNP = drop SAO 2 and 3 | 60,691,395 - 60,683,019 = 8,376 = SAO=2 (2,753) + SAO=3 (5,623) in our copy (section 4.1). Pansoma drops only SAO=2 | README table; local count | [inference] from the exact count match | R2.6, R3.m2 |
| C9 | Our four PoN files are the same files ClairS-TO uses | 1000G PoN 2,609,566 records (local count = README total); dbSNP 60,691,395 (local = README); CoLoRSdb 49,550,902 (Zenodo stats = README) | local counts; README | gnomAD total not recounted locally [unverified count; same GATK URL] | R2.6 |
| C10 | Platform models (v0.5.1) | ONT R10.4.1: dorado SUP 5 kHz (r1041_e82_400bps_sup_v420; `_ss`, `_ssrs`, plain), dorado SUP 4 kHz (v410), dorado HAC 5 kHz v6.0.0 (new in v0.5.1), dorado HAC 4 kHz, Guppy6 SUP 4 kHz, Guppy6 HAC 5 kHz; PacBio `hifi_revio` (Revio, SMRTbell prep kit 3.0; `_ss`, `_ssrs`); Illumina `ilmn` (NovaSeq/HiSeq X; `_ss`, `_ssrs`). No dorado SUP v5 model. Reference GRCh38_no_alt (ONT, HiFi) / GRCh38 (Illumina) | README "Pre-trained Models" table | e.g. "ONT ... r1041_e82_400bps_sup_v420 ... Dorado SUP ... Nov. 20, 2024 ... `ont_r10_dorado_sup_5khz_ssrs`" | R3.m2 |
| C11 | PoN effect in ClairS-TO itself (context) | COLO829 50x ONT, SSRS F1: no PoN 2.31%, all four 76.49% | Nat Commun Results | "ClairS-TO SSRS achieved F1-scores of 2.31, 39.56, 44.80, 3.87, 71.19, and 76.49% when using no PoNs, gnomAD only, dbSNP only, 1000 G PoN only, CoLoRSdb only, and all four PoNs, respectively." | R2.6 |
| C12 | v0.3.0 introduced CoLoRSdb and ssrs | v0.3.0 (2024-10-11) | README changelog | "Added using CoLoRSdb ... as a PoN for tagging non-somatic variant. The idea was inspired by Park et al., 2024." | R3.m2 |

### 3.3 The four population / PoN resources (`/scratch/jshen/data/Pansoma/panel_of_normal_VCFs`)

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| P1 | gnomAD file = GATK `af-only-gnomad.hg38.vcf.gz`, lifted from GRCh37 in Aug 2017 | Header: `bcftools norm ... af-only-gnomad.raw.sites.b37.vcf.gz; Date=Thu Aug 17 17:16:02 2017` and a `Hg19toGRCh38` step on 2017-08-22; INFO has only AC and AF | local header | `##bcftools_normCommand=norm -m +any --do-not-normalize -Oz -o zeta_af-only-gnomad_Hg19toGRCh38.vcf.gz delta_af-only-gnomad_Hg19toGRCh38.vcf.gz; Date=Tue Aug 22 17:14:23 2017` | R2.6 |
| P2 | Which gnomAD release / samples | gnomAD v2 (GRCh37) era, genome callset only: implied AN = AC/AF is 40,000-41,000 for 86 % of chr1:1-60 Mb records with AC >= 2,000 (273,255 / 316,383), max 42,600; at the coding SNP rs1801133 AC 10,857, AF 0.267 (AN about 40,660). No exome-scale AN (gnomAD v2 exomes would give up to about 246,000). AN 40,628 = 2 x 20,314 at chr1:10,067 | local computation (`implied_an_gnomad.sh`, section 4.2); gnomAD posts https://gnomad.broadinstitute.org/news/2017-02-the-genome-aggregation-database/ ("exome sequence data from 123,136 individuals and whole genome sequencing from 15,496 individuals") and https://gnomad.broadinstitute.org/news/2018-10-gnomad-v2-1/ ("October 17, 2018"; 125,748 exomes, 15,708 genomes) | The release label "r2.1" used by ClairS-TO cannot be exact: the file (Aug 2017) predates gnomAD v2.1 (2018-10-17) [inference]. AFs come from about 20.3 k genomes, more than the 15,496 post-QC v2.0 genomes, consistent with the "raw" (pre-sample-QC) genome callset named in the header [inference]. | R2.6 |
| P3 | gnomAD AF floor of 1e-4 in Pansoma's SNV rule | With AN about 40,600, AF >= 1e-4 means AC >= 5 | arithmetic on P2 | [inference] | R2.6 |
| P4 | dbSNP 138 | 60,691,395 records; SAO 0: 60,665,141; SAO 1: 17,878; SAO 2 (somatic): 2,753; SAO 3 (both): 5,623; dbSNP_BUILD_ID=138, fileDate=20130806 | local count (`dbsnp_sao_counts.sh`); header | `##INFO=<ID=SAO,...,Description="Variant Allele Origin: 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both">` | R2.6 |
| P5 | GATK 1000G PoN | 2,609,566 sites; Mutect2 2.1-beta calls on 1000 Genomes exome BAMs (bwamem GRCh38DH 20150826) [inference from the header path `1kg/exome_pon`]; the header keeps one sample's command (HG02775, PJL); sites-only | local header + count | `##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2  --tumorSampleName HG02775 ... --input shlee-dev/1kg/exome_pon/CreateMutect2Pon1of2/.../HG02775.alt_bwamem_GRCh38DH.20150826.PJL.exome.bam` ; number and identity of the other samples: [unverified, no public list found] | R2.6 |
| P6 | CoLoRSdb v1.1.0 | 2024-07-31, doi:10.5281/zenodo.13145123; 1,381 samples (autosomes), 982 (chrX), 516 (chrY); 49,550,902 small-variant records; PacBio HiFi, DeepVariant + GLnexus | Zenodo README + stats | "Variants were merged from 1381 samples for autosomes and chrM, from 982 samples for chrX ..., and from 516 samples (all male) for chrY." | R2.6 |
| P7 | CoLoRSdb cohorts | Children's Mercy rare-disease cohort (largest, 85 % European), HPRC (127 genomes "from 1,000 genomes"), HGSVC (37), SolveRD (100), HudsonAlpha (50 probands + 30 parents); no per-sample list published | https://colorsdb.org/ "Sample Description" | "Other population control cohorts include genomes from HPRC (127 genomes @30-40x, from 1,000 genomes) and HGSVC (37 genomes @30-40x)." | R2.6, R1.5 |
| P8 | CoLoRSdb has no singletons | chr20: 1,234,398 records, 0 with AC = 1, minimum AC = 2 | local (`colorsdb_chr20_min_ac.txt`) | So a contributor's private heterozygous allele is not in CoLoRSdb; with AN up to 2,762, AF >= 1e-4 (Pansoma's floor) admits every record [inference] | R2.6 |
| P9 | Can HG008 be in these resources? | Not in gnomAD af-only (2017), dbSNP 138 (2013) or the 1000G PoN (2017): HG008 data collection began Fall 2023. CoLoRSdb (2024) cohorts are not GIAB; the private-allele test (4.4) shows HG008-N with the lowest-tier match rates | GIAB HG008 v0.3 README ("Date of data collection: Fall 2023 - Ongoing"); dates in P1, P4, P5; section 4.4 | [inference] for CoLoRSdb | R2.6 |
| P10 | COLO829 / COLO829BL, HCC1395(BL), HCC1954, H1437, H2009 in the resources | Not 1000 Genomes samples; not among CoLoRSdb's listed cohorts; their normals' rare alleles match gnomAD / CoLoRSdb least of all 12 individuals tested (4.4). gnomAD sample list not public | sections 4.4, P7 | [inference]; membership in gnomAD v2 cannot be excluded from a list [unverified] | R2.6 |
| P11 | HapMap donors in the resources | HG00438, HG02257, HG02486, HG02622 are 1000 Genomes samples in HPRC; HG002, HG005 are GIAB/PGP samples in HPRC release 2. CoLoRSdb includes 127 HPRC genomes; their rare alleles match CoLoRSdb 2-5 x more often than the cancer-line normals' (4.4) | HPRC release-2 metadata `hprc_release2_sample_metadata.csv`; P7; 4.4 | consistent with membership in CoLoRSdb [inference]; 1000G PoN sample list [unverified] | R2.6, R1.5, R2.1 |

### 3.4 HPRC v1.1 and v2.1 graphs

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| H1 | v1.1 individuals | 44 individuals (47 assembled minus HG002, HG005, NA19240) + GRCh38 + CHM13 = 90 haplotypes | https://github.com/human-pangenomics/hpp_pangenome_resources/blob/main/README-v1.md | "Of the 47 samples assembled (94 assemblies) in year 1, all but three samples were included in graph constructions (HG002, HG005 and NA19240 were excluded for evaluation purposes). GRCh38 and CHM13 were added to make the total number of haplotypes included 90." | R3.m3, R1.5 |
| H2 | v1.1 populations (44) | AFR 23, AMR 16, EAS 4, SAS 1, EUR 0 | `analysis/graph_absorbed_somatic_20261001/hprc_individuals.tsv` (88 sample haplotypes) and the HPRC year-1 metadata (`tmp/hprc_v21_samples/hprc_year1_sample_metadata.txt`, HPP_Year1_Assemblies/sample_metadata); NA21309 (MKK, HapMap 3) counted as AFR | local count | R3.m3 |
| H3 | d9 meaning | `cactus-pangenome --filter 9`: nodes covered by fewer than 9 of 90 haplotypes removed; reference paths are never filtered, and the GBZ has two of them: GRCh38 and CHM13 (`vg paths -M -S <sample>`, vg 1.77, 2026-10-15: 179 GRCh38 and 109 CHM13 subpaths, all with SENSE = REFERENCE), so CHM13 sequence stays in the graph even where fewer than 9 haplotypes carry it [inference from the quote and the path sense; the ~100 CHM13 subpaths fit clipping, not a 9 / 90 frequency filter]; "about 10 %" of the haplotypes | https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md (HPRC v1.1 section); doc/mc-pangenomes/hprc-v1.1-mc.md (command line with `--filter 9`) | "we use `--filter 9` and `--giraffe` to make giraffe indexes on the subgraph covered by at least 9/90 haplotypes (filter does not apply to the reference)." README-v1: "created by dropping nodes and edges supported by fewer than 10% of haplotypes" | R3.m3 |
| H4 | v2.1 individuals | 232 assembled, 231 in the graph (HG00272 held out, chrX misassembly) + GRCh38 + CHM13 = 464 haplotypes; our d46 GBZ has 462 sample haplotypes from 231 samples + 2 references | hpp_pangenome_resources README; local `tmp/hprc_v21_samples/haplotypes.tsv` (from `vg paths -M` of the d46 GBZ) | "Of the 232 samples assembled (464 assemblies) in year 2, all but one were included in graph construction (HG00272 was held out due to an apparent large-scale misassembly in chrX). GRCh38 and CHM13 were added to make the total number of haplotypes included 464." (The HPRC2 preprint abstract says "HPRC2 contributes 460 haplotypes"; doi:10.64898/2026.07.21.739710.) | R3.m3 |
| H5 | v2.1 populations (231) | AFR 70 (68 in 1000G populations ACB 12, ASW 6, ESN 10, GWD 17, LWK 6, MSL 9, YRI 8; plus NA21309 MKK and HG06807 ASL "African Americans living in St. Louis"), AMR 44, EAS 50, SAS 36, EUR 29 (FIN 8, GBR 10, IBS 3, TSI 8), no population label 2 (HG002, HG005; GIAB/PGP Ashkenazi Jewish and Chinese) | `hprc_release2_sample_metadata.csv` (github.com/human-pangenomics/hprc_intermediate_assembly, data_tables/sample) joined to the GBZ sample list; HPRC year-1 metadata for HG002 / HG005 | The previous note "227 with metadata (AFR 68, EAS 50, AMR 44, SAS 36, EUR 29)" is right; the remaining four are HG002, HG005, HG06807, NA21309. Year-1 metadata: "HG002 HPRC_PLUS male 3140 PGP - Ashkenazim Jewish", "HG005 HPRC_PLUS male 3150 PGP - Chinese" | R3.m3 |
| H6 | d46 meaning | `--filter 46` (46 of 464, about 10 %); as in H3, GRCh38 and CHM13 are reference paths of the GBZ (`vg paths -M`, vg 1.77, 2026-10-15: 195 GRCh38 and 93 CHM13 subpaths, all REFERENCE) and are not filtered [inference, as in H3] | hprc-v2.1 S3 README (cactus-pangenome command lines); hpp_pangenome_resources README | "The "AF-Filtered VG indexes" above were created by dropping nodes and edges supported by fewer than 10% of haplotypes" | R3.m3 |
| H7 | HG002 / HG005 / NA19240 in v2.1 | Included in the default v2.1 graphs; separate "benchmark-graphs" exclude them | hpp_pangenome_resources README | "Versions of the CHM13 and GRCh38-based graphs are available here with `HG002`, `HG005`, and `NA19240` excluded. Use these if, for example, you are running Genome-In-A-Bottle benchmarks." | R1.5, R3.m3 |
| H8 | HapMap donors in the graphs | v1.1: HG00438, HG02257, HG02486, HG02622 (4 of 6); v2.1 (our d46): all 6 (adds HG002 and the 83.5 % background HG005) | local GBZ sample lists (v1.1 `hprc_samples.txt`; d46 `haplotypes.tsv`) | correction: the drafts say "four of its six donors are in both HPRC graphs"; for v2.1 it is six of six | R1.5, R2.1, R3.m3 |
| H9 | Benchmark donors in the graphs | HG008, COLO829, HCC1395, HCC1954, H1437, H2009 (and HCC1937, Hs578T) are in neither sample list | local lists; release-2 metadata | relatives cannot be checked for anonymised cell lines [unverified] | R3.m3 |
| H10 | Donor-removed graph already exists | `mgi_non_HG005.d9`: seqfile = 40 HPRC year-1 individuals + GRCh38 + CHM13 (82 haplotypes); none of the six HapMap donors and not NA19240 | `/scratch/jshen/data/non_HG005_VG_Indexes/mgi_non_HG005.seqfile` (84 lines incl. `_MINIGRAPH_`) | Makes E23 mostly done (graph and indexes exist); the HapMap tensor chain on it was submitted 2026-10-12 per the session memory [unverified here] | R1.5, R2.1 |

### 3.5 GIAB HG008 somatic small-variant benchmark

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| G1 | Version used here | v0.2 draft, local `HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz`: 18,952 records (9,609 SNV, 9,343 INDEL); `_all.bed` 41,835 regions, 2,809,425,965 bp; VCF header dates 2025-09-22 | local counts and header | `##bcftools_mergeCommand=merge --force-samples vcfeval_HG8Tchr19fixvHG8Nz200k_hap1|fn_somatic.vcf.gz; Date=Mon Sep 22 15:20:10 2025` | R2.2 |
| G2 | Public status | v0.2 is not in the public GIAB FTP analysis directory (only smvar V0.3 is); v0.3 (2026-04-25) calls itself the "initial public draft"; 18,946 variants (9,607 SNV, 9,339 INDEL), all.bed 42,480 regions, 2.808 Gb | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-smvar_DraftBenchmark_V0.3-20260425/README.md ; https://www.nist.gov/programs-projects/cancer-genome-bottle | "V0.3, current version, initial public draft, 2026-04-25 ... First public draft of the phased HG008-T somatic small variant (SNV/INDEL) benchmark, derived primarily from a tumor-normal assembly comparison." The analysis FTP lists smvar only as V0.3 (stvar V0.1-V0.3, stvar-CNV V0.4-V0.5) | R2.2 |
| G3 | How it is made | Tumor (HG008-T v3.2) vs normal (HG008-N v6.3, verkko 2.2 with HiFi, herro-corrected ONT, Hi-C) assembly comparison; candidates classified by repun with technology-specific AF/DP patterns from PacBio HiFi, ONT, Illumina / Element alignments to the normal assembly and GRCh38; liftOver to GRCh38; three BED tiers; nine community call sets used to curate FPs/FNs | v0.3 README "METHODOLOGICAL INFORMATION" | "Built from the HG008 normal assembly v6.3 (verkko 2.2 with HiFi, herro-corrected ONT, and Hi-C data) and HG008 tumor assembly v3.2." / "Technologies include PacBio HiFi, Oxford Nanopore long reads, and Illumina/Element short reads aligned to both the normal assembly and GRCh38." | R2.2 |
| G4 | Germline-modifying somatic variants | About 44 % of confirmed somatic variants change representation on GRCh38 | v0.3 README Limitations 3; preprint doi:10.64898/2026.05.01.722316 | "~44% of confirmed somatic variants on GRCh38 have a different representation on GRCh38 than on the assembly" | R2.2, R2.3 |
| G5 | Recommended comparison | aardvark (vcfeval and hap.py also tested); filter VAF < ~5-10 % calls (clonal benchmark) | v0.3 README | "We currently recommend aardvark for benchmarking somatic small variant calls against this draft benchmark, and we have also tested rtg vcfeval and hap.py." | R2.2, R3.m2 |
| G6 | HG008 donor | Female of European genetic ancestry, PDAC cell line + normal pancreatic and duodenal tissue; cell line ATCC CRL-3734 | v0.3 README Background; Sci Data doi:10.1038/s41597-025-05438-2 | "derived from a female of European genetic ancestry with explicit consent for public release of genomic data" | R3.m3, R2.6 |
| G7 | GIAB stratifications | v3.5 on the GIAB FTP (dir dated 2024-09-16; GRCh38 tarball 1.3 GB) and a newer v3.6 (2025-04-25); not local | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/ | directory listing | R1.3 / E10 |

### 3.6 SMaHT COLO829 truth (COLO829BLT50 reference set)

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| S1 | Publication | Ha et al. (SMaHT), bioRxiv 10.1101/2025.10.13.681545 (v1 2025-10-14), PMC12632906 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12632906/ | title "Comprehensive benchmarking of somatic single-nucleotide variant and indel detection at ultra-low allele fractions using short- and long-read data" | R2.2 |
| S2 | Construction | Illumina candidates from four callers (Mutect2, Strelka2, RUFUS, VarNet) on COLO829 200x vs COLO829BL 300x, validated with PacBio (COLO829 300x, COLO829BL 180x): ALT present in COLO829, absent in COLO829BL | PMC12632906 Fig. 1 legend / Results | "Reference set generation involves two steps: candidate variant identification from Illumina short-read data using four complementary somatic callers (Mutect2, Strelka2, RUFUS, and VarNet), followed by validation with PacBio long reads." (The GitHub README says "three variant callers"; the paper says four.) | R2.2 |
| S3 | Size | 44,005 SNVs, 961 insertions, 1,098 deletions; identical to our `COLO829T_somatic_snv_indel.vcf.gz` (44,005 SNV + 2,059 INDEL). The abstract says 44,008 SNVs | PMC12632906 Results; `COLO829T_truth/README.txt` | "Validated SNVs and indels from all methods were then combined into the final reference set, comprising 44,005 SNVs, 961 insertions, and 1,098 deletions" | R2.2 |
| S4 | Region classes (v2 BEDs) | easy = 1000 Genomes strict mask; difficult = in PanMask strict (pm151b) but outside the 1KG strict mask; extreme = outside PanMask; genome shares 74 / 10 / 16 %; files `SMaHT_{easy,difficult,extreme}_v2.bed.gz` | PMC12632906 Methods "Evaluation across different genomic contexts"; github.com/parklab/SMaHT_SNV_COLO829BLT50_HAPMAP `3.negative_control/utils/IntersectRegions.sh` | "Easy regions were defined by the 1000 Genomes strict mask ... Difficult regions corresponded to sites within the PanMask strict (pm151b) set but outside the 1KG strict mask ... Extreme regions encompassed loci outside PanMask" | R2.2, R1.3 |
| S5 | Version | Our files: qfu `Step1_Validated_union_{snv,ind}_VAF_RGN.VAF_added`; region BEDs `_v2`. Public copies: `Resource/COLO829BLT50_reference_set/Truthset_{SNV,Indel}_COLO829.vcf.gz` in the parklab repository; no version number found [unverified whether identical to ours] | local README; GitHub | - | R2.2 |

### 3.7 SEQC2 HCC1395 truth and titrations

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| Q1 | Construction | 21 replicates x 3 aligners (BWA-MEM, Bowtie2, NovoAlign) x 6 callers (MuTect2, SomaticSniper, VarDict, MuSE, Strelka2, TNscope) = 378 call sets, scored by SomaticSeq and NeuSomatic; HighConf / MedConf = reference set; Ion Torrent, AmpliSeq, HiSeq WES and PacBio for validation | Fang et al., Nat Biotechnol 39, 1151-1160 (2021), doi:10.1038/s41587-021-00993-6, PMC8532138 | "six mutation callers (MuTect2, SomaticSniper, VarDict, MuSE, Strelka2 and TNscope). We used SomaticSeq and NeuSomatic to score each mutation call based on evidence collected from the 378 call sets." | R2.2, R1.m1 |
| Q2 | v1.2.1 vs v1.2 | v1.2.1 = v1.2 calls restricted to `High-Confidence_Regions_v1.2.bed` (removes 113 SNVs and 297 INDELs) | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2.1/README.md | "There are 113 SNVs and 297 indels in those subtracted regions ... All the calls there are inside `High-Confidence_Regions_v1.2.bed`." | R2.2 |
| Q3 | Short-read limitation (the paper's own) | Callable regions limited to Illumina-analysable genome; PacBio validation had modest depth | PMC8532138 Discussion | "First, reference somatic mutations were called in regions of the genome limited to those that could be analyzed with the Illumina sequencing platforms." INDEL superset `sINDEL.MDKT.superSet.v1.2` (4 callers) vs SNV `sSNV.MSDUKT` (6) [file names on FTP; the letter mapping to callers is an inference]. Homopolymer share (15.5 % of SEQC2 INDELs in HP >= 7 vs 74.9 % for HG008T) is the project's own number (detailed doc R1.m1), not from SEQC2 | R1.m1, R2.2 |
| Q4 | Titration series | DNA mixtures of HCC1395 into HCC1395BL at 75, 50, 20, 10, 5 % tumor; triplicate HiSeq 4000 WGS by Genentech ("SPP_GT"); library names encode tumor:normal ratio | Xiao et al., Nat Biotechnol 39, 1141-1150 (2021), doi:10.1038/s41587-021-00994-5, PMC8506910; NCBI SRA runinfo (BioProject PRJNA489865) | "DNAs from HCC1395 and HCC1395BL were pooled at various ratios to create mixtures of 75%, 50%, 20%, 10%, and 5%. All libraries from these experiments were sequenced in triplicate on the HiSeq series by Genentech (GT)." | R2.4, R1.5 |
| Q5 | SRA runs per titration | 3:1 (75 %): SPP_GT_3-1_1 SRR7890895; _2 SRR7890894, SRR7955942-SRR7955944; _3 SRR7890864, SRR7961192-SRR7961194. 1:1 (50 %): _1 SRR7890936; _2 SRR7890957, SRR7963260-SRR7963262; _3 SRR7890956, SRR7963263-SRR7963265. 1:4 (20 %): _1 SRR7890903, SRR7947146-48, SRR7947188-90; _2 SRR7890902; _3 SRR7890897, SRR7947149-51. 1:9 (10 %): _1 SRR7890896, SRR7947152-54; _2 SRR7890899; _3 SRR7890898, SRR7947191-93. 1:19 (5 %): _1 SRR7890955; _2 SRR7890901; _3 SRR7890900. Pure tumor 1-0 and normal 0-1 also exist. All HiSeq 4000, 2x150 | NCBI E-utilities `efetch rettype=runinfo` for `PRJNA489865 AND SPP` (51 runs) | ratio-to-purity mapping is from the names (3-1, 1-1, 1-4, 1-9, 1-19 = 75, 50, 20, 10, 5 %) [inference, matches Q4 exactly] | R2.4, R1.5 |
| Q6 | Our local titration run | `SPP_GT_1-1_1.sra` = SRR7890936, 50 % purity, 1,533,849,407 spots x 300 bp (about 460 Gb, about 148x) | SRA runinfo | answers the doc's "purity 档位不明" | R2.4 |

### 3.8 The HapMap mixture benchmark

| # | Claim | Verified value | Source | Exact quote / location | Used by |
|---|---|---|---|---|---|
| M1 | What it is | SMaHT "HAPMAP6" lymphoblastoid mixture: 83.5 % HG005, 10 % HG02622, 2 % each HG002, HG02257, HG02486, 0.5 % HG00438; positives are the minor donors' germline alleles absent from HG005 (VAF 0.25-16.5 %) | SMaHT Network, bioRxiv 2025.10.09.678885, PMC12632347 (Methods); Ha et al. PMC12632906 | "Based on cell counts, volumes of individual cell suspensions were calculated to achieve the target mix ratio: 83.5% HG005, 10% HG02622, 2% each of HG002, HG02257, and HG02486, and 0.5% HG00438." | R1.5, R2.1 |
| M2 | Truth used in the manuscript | Pangenome-graph truth of Kong et al. (WashU), local v1.6 (`SMHTHAPMAP6_GRCh38_v1.6_somatic_benchmark_*`; somatic 7,870,692 records per `HapMap_truth/README.txt`) | Kong, Tang, Ruttenberg et al., bioRxiv 10.1101/2025.09.29.679336 (v2 2025-10-04); public README v1.0.0 (wangftp.wustl.edu/~jmacias/Graph_based_HapMap_Truth_Set) | "The SMaHT Consortium mixed six HapMap cell lines to create artificial somatic variants spanning 0.25% to 16.5% variant allele fractions. We developed a technology-agnostic method that builds pangenome graphs from individual assemblies" | R1.5, R2.1, R2.2 |
| M3 | Donor ancestry | HG00438 CHS, HG02257 ACB, HG02486 ACB, HG02622 GWD (1000 Genomes, HPRC); HG002 Ashkenazi Jewish, HG005 Han Chinese (GIAB/PGP) | HPRC release-2 metadata; year-1 metadata | (the web-search summary that called HG02622 "European/Ashkenazi" is wrong) | R1.5 |
| M4 | Companion paper cited by R1 | Fu, Xin, Miao et al., "Leveraging Human Pangenome for Improved Somatic Variant Detection", bioRxiv doi:10.64898/2026.01.04.697580 (v2 2026-01-05); benchmarks GRCh38, graph pangenomes and pangenome-inferred DSAs on the HapMap mixture and COLO829 | api.biorxiv.org | "we benchmarked somatic variant detection using GRCh38, graph-based pangenomes, and pangenome-inferred DSAs with a HapMap mixture dataset and the COLO829 melanoma cell line." | R1.5 |
| M5 | Donors in graphs | see H8 (4 of 6 in v1.1, 6 of 6 in v2.1) and H10 (donor-free 82-haplotype graph exists) | - | - | R1.5, R2.1 |

## 4. Local computations

All scripts are in this folder; large outputs in `/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts/`;
Slurm scripts and logs in `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E8b_literature_facts/` (`jobs.txt`).

### 4.1 dbSNP 138 by SAO (`dbsnp_sao_counts.sh` -> `dbsnp_sao_counts.tsv`, login node)

| SAO | meaning | records |
|---|---|---|
| 0 | unspecified | 60,665,141 |
| 1 | germline | 17,878 |
| 2 | somatic | 2,753 |
| 3 | both | 5,623 |
| total | | 60,691,395 |

Pansoma drops SAO=2 only (`scripts/filter_panel_of_normals.py`, `non_somatic=True` -> `SAO != 2`); ClairS-TO's
"non-somatic" file keeps 60,683,019 = total minus SAO 2 and 3.

### 4.2 Implied AN of the gnomAD file (`implied_an_gnomad.sh` -> `implied_an_gnomad.tsv`, login node)

chr1:1-60,000,000, records with AC >= 2,000: 316,383; implied AN = AC/AF falls in 40,000-40,999 for 273,255 (86.4 %),
39,000-39,999 for 13,888, 41,000-41,999 for 4,680; maximum 42,600. Spot checks: chr1:10,067 AC 3, AF 7.384e-05
(AN 40,628); rs1801133 (MTHFR coding) AC 10,857, AF 0.267 (AN about 40,660); CoLoRSdb at the same site AN 2,762, AF 0.303.

### 4.3 CoLoRSdb minimum AC (`colorsdb_chr20_min_ac.txt`, login node)

`bcftools query -r chr20 -f '%INFO/AC\n'`: 1,234,398 records, 0 with AC = 1, minimum AC 2.

### 4.4 Private-allele membership test, chr20 SNVs (`membership_private_alleles.py`, Slurm 387259, 5 min, 2.0 GB)

For each individual, its germline SNV alleles on chr20 are matched (POS, REF, ALT) to gnomAD (AC) and CoLoRSdb (AC).
"Rare matched %" = alleles found at AC 1-2 / (alleles absent + alleles at AC 1-2): a contributor's private alleles sit
at AC 1-2 in a resource that contains it, and are absent from one that does not. Full counts:
`membership_private_alleles.chr20.tsv`; ratios: `membership_private_alleles.chr20.summary.tsv`.

| individual | source of alleles | SNV alleles | gnomAD rare matched % | CoLoRSdb rare matched % |
|---|---|---|---|---|
| HG002 | HAPMAP6 v1.6 (alleles not in HG005) | 31,678 | 3.8 | 22.8 |
| HG00438 | same | 23,942 | 15.4 | 16.4 |
| HG02257 | same | 46,032 | 5.6 | 25.9 |
| HG02486 | same | 56,980 | 3.4 | 18.7 |
| HG02622 | same | 54,417 | 6.2 | 22.5 |
| HG005 | HapMap_truth germline | 78,737 | 5.3 | 14.5 |
| HG008-N | dipcall of HG008-N assembly | 82,979 | 2.5 | 7.9 |
| COLO829BL | dipcall of COLO829BL assembly | 79,119 | 1.5 | 6.7 |
| HCC1395BL | CASTLE DeepVariant hybrid, QUAL >= 10 | 93,422 | 0.8 | 4.4 |
| HCC1954 normal | same | 100,187 | 1.8 | 7.2 |
| H1437 normal | same | 96,863 | 0.6 | 6.4 |
| H2009 normal | same | 96,314 | 0.6 | 4.0 |

Interpretation [inference]: in CoLoRSdb the six HPRC / GIAB donors match 14.5-25.9 % against 4.0-7.9 % for HG008-N,
COLO829BL and the four CASTLE normals, even though CoLoRSdb is mostly European and HG008 / COLO829 are European. This
fits CoLoRSdb containing the HPRC donors (its website lists 127 HPRC genomes) and not containing the benchmark
lines. Because CoLoRSdb has no AC=1 records (4.3), the signal comes from AC=2 (a contributor's private homozygous
allele, or two carriers). For gnomAD the split is weaker (HapMap / GIAB donors 3.4-15.4 %, benchmark lines 0.6-2.5 %);
HG00438 stands out (15.4 %); ancestry and call-set differences (assembly dipcall vs DeepVariant hybrid) confound the
comparison, so gnomAD membership stays [unverified]. The dbSNP column of the full table (share of gnomAD-absent
alleles present in dbSNP 138) does not separate the groups (HG008-N 33.2 %, HapMap donors 20.0-46.6 %) and is not used.

### 4.5 Somatic truth alleles found in the resources (`somatic_truth_in_resources.py`, Slurm array 387253)

Six somatic truth sets DONE (tasks 0-5); HapMap contrast DONE (task 6, 2026-10-15; in `somatic_truth_in_resources.tsv`
since then). HapMap chr20 positives (donor alleles absent from HG005, all records, no BED): SNV n 129,857, Pansoma
rule 125,026 (96.3 %), ClairS-TO rule 123,210 (94.9 %); INDEL n 41,174, Pansoma rule 16,182 (39.3 %) under the rule
of this table (HP >= 7: 0.01, else 1e-3), 16,630 (40.4 %) under the 2026-10-16 rule (HP >= 8: 0.01, else 1e-4; job
450562, `/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts/somatic_truth_in_resources.HapMap_positives_chr20.pon_hp8_1e-4.tsv`),
ClairS-TO rule 19,755 (48.0 %). Used by R2.1 point 3 (96.3 % / 40.4 %).

For every ALT of every record of a truth VCF (all chromosomes, no BED; multi-allelic records split, so n can exceed
the record count by a few), the four resources are queried and the script counts: exact allele matches (POS, REF,
ALT) with AF floors for gnomAD / CoLoRSdb; dbSNP matches by SAO; position-only matches for 1000G PoN / CoLoRSdb;
the Pansoma rule, computed by importing the pipeline's own `scripts/filter_panel_of_normals.py` (unchanged since
a01969d; SNV: gnomAD and CoLoRSdb exact allele with AF >= 1e-4, dbSNP exact allele with SAO != 2, 1000G PoN exact
allele; INDEL: gnomAD and CoLoRSdb exact allele with AF >= 0.01 if the ALT changes a GRCh38 homopolymer >= 7 bp, else
AF >= 0.001; dbSNP and 1000G not used for INDELs); and the ClairS-TO default rule (C6, C7: gnomAD AF >= 1e-3 exact,
dbSNP SAO not 2 exact, 1000G PoN position, CoLoRSdb position with AF >= 1e-3). A "tagged" truth allele would be
removed from a call set by that rule, so the share is a recall cost, not a precision gain.

Tables: `somatic_truth_in_resources.tsv` (counts), `somatic_truth_in_resources.pct.tsv` (shares), written by
`bash collect.sh` from `/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts/somatic_truth_in_resources.<SET>.tsv`.
Two columns of the TSV are left out below because they hold no extra information: `CoLoRSdb_AF>=1e-4` equals
`CoLoRSdb_any` in all 12 rows and `CoLoRSdb_AC=1` is 0 in all 12 rows (no singletons, 4.3); `dbSNP_SAO!=2` equals
`dbSNP_any` in all rows.

Counts (share of n). Truth sets: GIAB HG008-T v0.2; SMaHT COLO829 (COLO829BLT50); SEQC2 HCC1395 v1.2.1;
CASTLE = DeepSomatic HCC1395-model benchmark (D8).

SNV

| truth set | n | Pansoma rule | ClairS-TO rule | gnomAD any | gnomAD AF >= 1e-4 | gnomAD AF >= 1e-3 | gnomAD AF >= 1e-2 | gnomAD AC <= 2 | dbSNP any | dbSNP SAO=2 | dbSNP SAO=3 | 1000G allele | 1000G position | CoLoRSdb any | CoLoRSdb AF >= 1e-3 | CoLoRSdb position |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HG008T GIAB v0.2 | 9,609 | 982 (10.2 %) | 906 (9.4 %) | 1,220 (12.7 %) | 666 (6.9 %) | 425 (4.4 %) | 273 (2.8 %) | 424 (4.4 %) | 522 (5.4 %) | 0 | 1 (<0.1 %) | 12 (0.1 %) | 24 (0.2 %) | 648 (6.7 %) | 593 (6.2 %) | 778 (8.1 %) |
| COLO829T SMaHT | 44,005 | 1,606 (3.6 %) | 1,378 (3.1 %) | 3,547 (8.1 %) | 1,149 (2.6 %) | 472 (1.1 %) | 186 (0.4 %) | 1,923 (4.4 %) | 860 (2.0 %) | 0 | 9 (<0.1 %) | 15 (<0.1 %) | 26 (0.1 %) | 551 (1.3 %) | 445 (1.0 %) | 955 (2.2 %) |
| HCC1395 SEQC2 v1.2.1 | 39,447 | 1,102 (2.8 %) | 1,127 (2.9 %) | 2,387 (6.1 %) | 772 (2.0 %) | 397 (1.0 %) | 212 (0.5 %) | 1,344 (3.4 %) | 604 (1.5 %) | 0 | 2 (<0.1 %) | 94 (0.2 %) | 111 (0.3 %) | 604 (1.5 %) | 527 (1.3 %) | 1,009 (2.6 %) |
| HCC1954 CASTLE | 24,033 | 1,956 (8.1 %) | 1,699 (7.1 %) | 2,896 (12.1 %) | 1,504 (6.3 %) | 900 (3.7 %) | 574 (2.4 %) | 1,113 (4.6 %) | 1,166 (4.9 %) | 0 | 1 (<0.1 %) | 170 (0.7 %) | 183 (0.8 %) | 1,199 (5.0 %) | 1,096 (4.6 %) | 1,460 (6.1 %) |
| H1437 CASTLE | 82,529 | 25,323 (30.7 %) | 25,591 (31.0 %) | 26,350 (31.9 %) | 24,438 (29.6 %) | 23,822 (28.9 %) | 22,982 (27.8 %) | 1,578 (1.9 %) | 23,982 (29.1 %) | 0 | 1 (<0.1 %) | 7,664 (9.3 %) | 7,709 (9.3 %) | 24,440 (29.6 %) | 24,329 (29.5 %) | 25,343 (30.7 %) |
| H2009 CASTLE | 161,918 | 11,445 (7.1 %) | 12,188 (7.5 %) | 15,079 (9.3 %) | 10,196 (6.3 %) | 8,880 (5.5 %) | 8,141 (5.0 %) | 4,094 (2.5 %) | 9,560 (5.9 %) | 0 | 1 (<0.1 %) | 2,676 (1.7 %) | 2,790 (1.7 %) | 9,548 (5.9 %) | 9,288 (5.7 %) | 11,760 (7.3 %) |

INDEL

| truth set | n | Pansoma rule | ClairS-TO rule | gnomAD any | gnomAD AF >= 1e-4 | gnomAD AF >= 1e-3 | gnomAD AF >= 1e-2 | gnomAD AC <= 2 | dbSNP any | dbSNP SAO=2 | dbSNP SAO=3 | 1000G allele | 1000G position | CoLoRSdb any | CoLoRSdb AF >= 1e-3 | CoLoRSdb position |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HG008T GIAB v0.2 | 9,362 | 5,783 (61.8 %) | 8,155 (87.1 %) | 1,237 (13.2 %) | 1,146 (12.2 %) | 1,022 (10.9 %) | 667 (7.1 %) | 55 (0.6 %) | 2,253 (24.1 %) | 0 | 0 | 19 (0.2 %) | 196 (2.1 %) | 7,550 (80.6 %) | 7,350 (78.5 %) | 8,154 (87.1 %) |
| COLO829T SMaHT | 2,059 | 1,070 (52.0 %) | 1,524 (74.0 %) | 356 (17.3 %) | 324 (15.7 %) | 280 (13.6 %) | 188 (9.1 %) | 24 (1.2 %) | 497 (24.1 %) | 0 | 0 | 9 (0.4 %) | 48 (2.3 %) | 1,388 (67.4 %) | 1,338 (65.0 %) | 1,521 (73.9 %) |
| HCC1395 SEQC2 v1.2.1 | 1,625 | 119 (7.3 %) | 224 (13.8 %) | 108 (6.6 %) | 97 (6.0 %) | 85 (5.2 %) | 52 (3.2 %) | 10 (0.6 %) | 78 (4.8 %) | 0 | 0 | 8 (0.5 %) | 14 (0.9 %) | 169 (10.4 %) | 156 (9.6 %) | 227 (14.0 %) |
| HCC1954 CASTLE | 2,556 | 528 (20.7 %) | 1,027 (40.2 %) | 654 (25.6 %) | 506 (19.8 %) | 358 (14.0 %) | 165 (6.5 %) | 105 (4.1 %) | 417 (16.3 %) | 0 | 0 | 17 (0.7 %) | 38 (1.5 %) | 810 (31.7 %) | 736 (28.8 %) | 1,039 (40.6 %) |
| H1437 CASTLE | 4,896 | 2,508 (51.2 %) | 2,767 (56.5 %) | 2,183 (44.6 %) | 2,105 (43.0 %) | 1,990 (40.6 %) | 1,843 (37.6 %) | 57 (1.2 %) | 2,297 (46.9 %) | 0 | 0 | 510 (10.4 %) | 532 (10.9 %) | 2,585 (52.8 %) | 2,550 (52.1 %) | 2,749 (56.1 %) |
| H2009 CASTLE | 6,519 | 1,422 (21.8 %) | 1,960 (30.1 %) | 1,217 (18.7 %) | 1,071 (16.4 %) | 946 (14.5 %) | 793 (12.2 %) | 103 (1.6 %) | 1,262 (19.4 %) | 0 | 0 | 208 (3.2 %) | 241 (3.7 %) | 1,717 (26.3 %) | 1,652 (25.3 %) | 1,963 (30.1 %) |

The n agree with the truth record counts elsewhere: HCC1395 39,447 + 1,625 = 41,072 records (CASTLE README); COLO829T
44,005 SNV + 2,059 INDEL (S3); HCC1954 / H1437 / H2009 26,589 / 87,425 / 168,437 alleles from 26,578 / 87,410 /
168,397 records (CASTLE README; the difference is split multi-allelic records).

Reading:

1. No leakage of tumor variants through dbSNP: 0 truth alleles of any set match a dbSNP 138 record flagged somatic
   (SAO=2), and at most 9 (COLO829T SNV) match SAO=3. HG008 data postdate all three 2013-2017 resources (P9); for the
   older lines (COLO829, HCC1395 and the CASTLE lines) the zero SAO=2 matches are the evidence.
2. Pansoma's SNV floor AF >= 1e-4 admits every CoLoRSdb record (the 1e-4 and "any" columns are identical in all rows),
   as predicted in P8.
3. SNV cost of the PoN step, all records [inference: the share of truth SNVs any caller loses after this rule]:
   Pansoma rule 2.8 % (HCC1395), 3.6 % (COLO829T), 7.1 % (H2009), 8.1 % (HCC1954), 10.2 % (HG008T), 30.7 % (H1437);
   ClairS-TO default 2.9 / 3.1 / 7.5 / 7.1 / 9.4 / 31.0 %. Pansoma tags more SNVs than ClairS-TO on HG008T, COLO829T and
   HCC1954 (its gnomAD floor is 1e-4 rather than 1e-3: HG008T gnomAD AF >= 1e-4 6.9 % vs >= 1e-3 4.4 %) and fewer on
   HCC1395, H1437 and H2009 (ClairS-TO's position-only 1000G PoN / CoLoRSdb matches).
4. INDEL cost: ClairS-TO's default tags more truth INDELs than Pansoma's rule in every set: HG008T 87.1 vs 61.8 %,
   COLO829T 74.0 vs 52.0 %, H1437 56.5 vs 51.2 %, HCC1954 40.2 vs 20.7 %, H2009 30.1 vs 21.8 %, HCC1395 13.8 vs 7.3 %.
   The driver is CoLoRSdb: 80.6 % (HG008T) and 67.4 % (COLO829T) of truth INDEL alleles are exact CoLoRSdb alleles and
   87.1 / 73.9 % share a position with a CoLoRSdb record, against 10.4 / 14.0 % for the Illumina-derived SEQC2 truth.
   [inference] the long-read-built truths (HG008T, COLO829T) are rich in homopolymer / STR length changes (Q3: 74.9 %
   of HG008T INDELs in HP >= 7 vs 15.5 % for SEQC2), and such alleles are common in CoLoRSdb.
5. Common population alleles in somatic SNV truth (gnomAD AF >= 1e-2): COLO829T 0.4 %, HCC1395 0.5 %, HCC1954 2.4 %,
   HG008T 2.8 %, H2009 5.0 %, H1437 27.8 %; exact 1000G PoN alleles: H1437 9.3 %, H2009 1.7 %, all others <= 0.7 %.
   [inference] A somatic mutation that happens to create a common germline allele is rare (0.4-0.5 % in the two
   independently validated sets, COLO829T and HCC1395), so most of the 22,982 H1437 truth SNVs with gnomAD AF >= 1 %
   are probably germline alleles labelled somatic in DeepSomatic's CASTLE benchmark (which the DeepSomatic authors
   call biased, D8); H2009 shows a weaker form of the same. Cause not examined here [unverified]; a match of these
   alleles against the H1437 normal's germline VCF would test it. Consequence: on H1437 every tumor-only caller that
   applies a PoN (Pansoma, ClairS-TO, DeepSomatic with `--use_default_pon_filtering=true`) loses about 31 % of the
   SNV truth by construction, so H1437 recall reflects the truth set, not the caller.
6. HG008T (assembly-derived truth) has more population-allele matches than COLO829T / HCC1395 (SNV Pansoma rule
   10.2 % vs 3.6 / 2.8 %); possible contributors are somatic changes at germline-variant sites (G4: about 44 % of
   HG008-T somatic variants change representation on GRCh38) [inference, not tested].
7. HapMap positives (chr20; germline alleles of the minor donors, the expected high-match contrast): PENDING,
   task 387253_6 held.

## 5. Corrections and verified values for the drafts

1. R1.5 / R2.1 / R3.m3: "four of its six donors are in both HPRC graphs" is wrong for v2.1: all six HapMap donors,
   including the 83.5 % background HG005 and HG002, are in v2.1 (H7, H8). Four are in v1.1.
2. R1.5 / E23: a donor-free graph already exists (`mgi_non_HG005.d9`, 40 individuals, 82 haplotypes, no HapMap donor,
   no NA19240) (H10). E23 is not a one-week build; only the HapMap tensors / labels / models on it remain.
3. R1.m1 / R3.m2: DeepSomatic held out chr1 in every model (D2, D14); comparisons on HCC1395 / HCC1954 / H1437 / H2009
   chr1 are on DeepSomatic's held-out chromosome. But the CASTLE truth for HCC1954 / H1437 / H2009 was made by
   DeepSomatic itself and the authors call it "strongly biased" toward it (D8). ClairS-TO ssrs was fine-tuned on
   HCC1937 / HCC1954 / H1437 / H2009 and the held-out chromosome is chr20 by the README, chr1 by the paper (C4).
   The comparators are not re-run (user decision 2026-10-13; the user's existing ClairS-TO / DeepSomatic results are
   used), so the letter should name the ClairS-TO model of those results and, if it is `_ssrs`, state that its
   held-out chromosome on HCC1954 / H1437 / H2009 is ambiguous (C4) unless HKU-BAL clarifies.
4. R3.m2: versions to cite: DeepSomatic 1.10.0 (2026-03-05), ClairS-TO 0.5.1 (2026-09-09). DeepSomatic tumor-only uses
   a population-AF channel by default (1000G DeepVariant for Illumina, CoLoRSdb for PacBio / ONT) and applies its PoN
   filter only with `--use_default_pon_filtering=true` (D12, D13); ClairS-TO tags with four PoNs by default (C6, C7).
5. R2.6: our four PoN files are exactly ClairS-TO's four source files (C9), but the rules differ: Pansoma uses
   exact-allele matching for all four, AF >= 1e-4 for SNVs in gnomAD / CoLoRSdb (= AC >= 5 in gnomAD, any record in
   CoLoRSdb), and keeps dbSNP SAO=3 records; ClairS-TO uses AF >= 1e-3, position-only for 1000G PoN and CoLoRSdb, and
   drops SAO 2 and 3 (C7, C8, P3, P8).
6. R2.6: "gnomAD r2.1" is not the right label for the GATK file; it is a 2017 GRCh37->GRCh38 liftover of gnomAD v2
   genome AFs from about 20.3 k genomes (P1, P2).
7. R2.2: HG008 truth v0.2 was a pre-public draft; v0.3 (2026-04-25) is the first public version with nearly the same
   counts (G1, G2). Either cite v0.2 as "pre-release draft provided by GIAB" or move to v0.3 (needs a 0.3 MB download).
8. R2.4: the local SEQC2 titration run is 50 % purity (Q6); 75 / 20 / 10 / 5 % runs are listed in Q5.
9. R2.2: SMaHT COLO829 truth = four Illumina callers + PacBio validation (S2), not "three callers" (GitHub README).
10. R2.6: no somatic truth allele of HG008T, COLO829T, HCC1395, HCC1954, H1437 or H2009 is a dbSNP-somatic (SAO=2)
    record (4.5), so the resources do not carry the benchmarks' tumor variants; the matches that exist are
    population alleles and cost recall (Pansoma rule: SNV 2.8-10.2 % on five sets; INDEL 7.3-61.8 %). For INDELs
    Pansoma's allele + AF rule removes fewer truth alleles than ClairS-TO's default rule in all six sets.
11. R1.m1 / R2.2 / R3.m2: the H1437 CASTLE truth is unlike the others: 27.8 % of its somatic SNVs are gnomAD alleles
    with AF >= 1 % (0.4-5.0 % elsewhere) and 30.7 % are tagged by Pansoma's PoN rule (31.0 % by ClairS-TO's) (4.5).
    Any H1437 comparison should say that about 31 % of its SNV truth is removed by any PoN, whatever the caller.
    Consistent with this, only about 70 % of the H1437 truth SNVs get a tensor on every platform (0.711 / 0.709 / 0.713,
    Illumina / PacBio / ONT, d46 `truth_recall.json`, 2026-10-15), against 0.915-0.936 for HCC1954 and H2009 SNVs.
12. R3.1 / R3.6 / R3.m3 / R1.5 / R2.1 (2026-10-15): both HPRC GBZs carry CHM13 as a reference-sense path next to GRCh38
    (H3, H6). The frequency filter therefore does not remove CHM13 sequence, and the GRCh38-only vs HPRC contrast
    mixes reference completeness (CHM13) with haplotype diversity. Earlier wording that only the GRCh38 path is
    unfiltered, or that CHM13-only sequence is filtered out, is wrong.

## 6. Status per reviewer point (this experiment's share)

| Point | What E8b had to deliver | Status | Remaining |
|---|---|---|---|
| R1.m1 | DeepSomatic training lines and held-out chromosomes; which lines are fair to compare on | Resolved (D2-D8, D14, C3, C4); H1437 truth caveat added (4.5, item 11 of section 5) | ClairS-TO held-out chromosome (chr20 vs chr1) needs an answer from HKU-BAL (comparators are not re-run, user decision); HCC1937 not downloaded (optional E12b) |
| R3.m2 | Versions, checkpoints, tumor-only inputs, PoN defaults of both comparators | Resolved for the facts (D1, D10-D18, C1-C12) | Exact command lines and models come from the user's existing ClairS-TO / DeepSomatic results (not re-run) |
| R2.6 | Composition of the four resources; whether benchmark donors / tumor variants are in them | Resolved for the six somatic truth sets (P1-P11, 4.1-4.5: 0 dbSNP-somatic matches in any set; matches are population alleles) | HapMap contrast row of 4.5 PENDING: task 387253_6 held until the user agrees a cluster cap; gnomAD v2 and 1000G PoN sample lists are not public (membership only inferred) |
| R3.m3 | HPRC release, individuals, populations, filtering, benchmark donors in graphs | Resolved (H1-H9) | none |
| R1.5 | HapMap composition, donors in graphs, leakage, companion paper | Resolved (M1-M5, H7-H10) | Decision: drop HapMap or use the existing donor-free graph |
| R2.1 | Same facts as R1.5 | Resolved | same decision |
| R2.2 | How each truth set was generated | Resolved (G1-G6, S1-S5, Q1-Q3, D8, M2); H1437 truth has 27.8 % common-allele SNVs (4.5) | Decide v0.2 vs v0.3 for HG008; decide how to present H1437; confirm our COLO829 truth equals the public parklab files [unverified] |

## 7. Caveats

- bioRxiv full text could not be fetched (rate limit); DeepSomatic and ClairS-TO facts come from the journal versions,
  which can differ from the preprints.
- Sample lists of the GATK gnomAD file, the GATK 1000G PoN and CoLoRSdb are not public per individual; section 4.4
  is statistical evidence, not a list. The 4.4 allele lists differ in origin (assembly dipcall, HAPMAP6 graph truth,
  DeepVariant hybrid), which affects absent / present ratios.
- 4.5 counts all truth records on all chromosomes (no BED, not only chr1) and, for HapMap, chr20 only; shares inside
  the benchmark BEDs or on chr1 can differ. The "rule" columns are what a PoN would remove from a perfect call set,
  not measured caller recall.
- 4.5 Pansoma rule = `scripts/filter_panel_of_normals.py` as on disk on 2026-10-13 (last commit a01969d); ClairS-TO
  rule = re-implementation of the documented default (C6, C7), not ClairS-TO's own code.
- The H1437 reading in 4.5 (germline alleles in the truth) is an inference from population AFs; it was not checked
  against the H1437 normal.
- No data were downloaded; GIAB v0.3 HG008 files, GIAB stratifications and the public SMaHT COLO829 truth files would
  each need a download and the user's approval.

## 8. Commands

```
# login node (short)
bash implied_an_gnomad.sh > implied_an_gnomad.tsv
bash dbsnp_sao_counts.sh > dbsnp_sao_counts.tsv
bash colorsdb_min_ac.sh > colorsdb_chr20_min_ac.txt
# Slurm (general, 1 CPU each)
sbatch /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E8b_literature_facts/membership.sh          # 387259
sbatch /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E8b_literature_facts/somatic_truth_array.sh # 387253_[0-6]
# collect (re-run 2026-10-13 22:54; reads finished outputs only, prints the missing HapMap file)
bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8b_literature_facts/collect.sh
```

Jobs (`jobs.txt`: `somres_array=387253 membership=387259`; logs in `.../tmp/reviewer_cpu_20261013/E8b_literature_facts/logs/`):

| job | what | state (sacct, 2026-10-13) | elapsed | MaxRSS |
|---|---|---|---|---|
| 387259 | membership test, chr20 (4.4) | COMPLETED, exit 0 | 0:06:21 | 2.0 GB |
| 387253_0 | HG008T_GIAB_v0.2 (4.5) | COMPLETED, exit 0 | 0:38:44 | 125 MB |
| 387253_1 | COLO829T_SMaHT | COMPLETED, exit 0 | 1:28:40 | 125 MB |
| 387253_2 | HCC1395_SEQC2_v1.2.1 | COMPLETED, exit 0 | 1:26:06 | 126 MB |
| 387253_3 | HCC1954_CASTLE | COMPLETED, exit 0 | 1:14:57 | 125 MB |
| 387253_4 | H1437_CASTLE | COMPLETED, exit 0 | 1:33:16 | 127 MB |
| 387253_5 | H2009_CASTLE | COMPLETED, exit 0 | 0:35:46 | 131 MB |
| 387253_6 | HapMap_positives_chr20 | PENDING (JobHeldUser): held with the other reviewer jobs after the lab-cluster footprint complaint; not to be released until the user agrees a cap | - | - |

The COLO829T log has only header warnings (INFO tags `VAF_Ill`, `RGN`, ... not defined in the truth VCF header).
