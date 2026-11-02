# Fact sheet: full-truth denominators for R1.4, R1.5 / R2.1, R1.6, R2.6, R1.8 (2026-10-30)

Labels used here: **EXISTS** (the number is already in a repo or data file, path given), **COMPUTED NOW** (computed today
by a script in `/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators/`; that directory is git-ignored and
holds the script, its output and the log), **MISSING** (not found and not computable here).
No tracked repo file was edited.

Scripts and outputs computed today (directory `D = /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators`):

| script | what | output | how |
|---|---|---|---|
| `D/absorption_denominators.py` | absorption with full-truth denominators, from E5b `per_truth.tsv.gz` and E5c `per_locus.tsv.gz` | `D/absorption_denominators.md` / `.json` | login node |
| `D/label_status.py` (`ARM=d46` / `d9` / `linear`) | per-truth label status (exact / partial / filtered / no candidate), chr1 and chr2–22, three seed splits; tensors by label × reason × off-reference | `D/label_status.json`, `D/label_status.d9.json`, `D/label_status.linear.json` | Slurm 462879, 462885 (miao); E21 index caches, checked against each set's current `labels.manifest.json` |
| `D/hapmap_m2_fulltruth.sbatch` | Mutect2 HapMap chr1 rerun scored on the **full** chr1 BED (graph-carried truth counted as FN) | `D/hapmap_m2_full/call_vaf_by_value.Mutect2_{pass,pass_cluster}.{SNV,INDEL}.md`, `D/hapmap_m2_full/eval.*` | Slurm 462880 (miao) |
| `D/r26_filter_table.py` | HG008T chr1 with / without population filtering, all arms, from every run's `report.json` | `D/r26_filter_table.md` / `.json` | login node |

Two checks passed. The absorption script reproduces E5b (1,467 / 1,418 / 49; 347 / 341 / 6) and E5c (41.5 / 42.1 / 29.4 %;
net 5.8 / 7.0 / 3.5 %). The seed-0 training-block counts of `label_status.py` equal E21 `set_split_counts.tsv`
(for example HG008T Illumina d46 SNV: somatic 7,274, germline 256,583, non 1,619,483).

---

## 1. R1.4: absorption with full-truth denominators

### 1a. Definitions (each truth allele, per platform; "tagged" = the current population filter: allele-exact, INDEL rule HP ≥ 8 AF ≥ 0.01, else AF ≥ 1e-4)

- **strict**: this is what the letter's 1,467 / 347 means (E5b `s5_tables.py` §2/§3). The allele has a representative
  tensor on the GRCh38-only graph, no candidate on the HPRC graph, and a graph path spells its ALT. The letter's wording
  ("get no candidate on at least one platform although a graph path spells their ALT") leaves out the GRCh38-only
  condition.
- **arm-only**: no candidate on the HPRC graph and the ALT is spelled by the graph, with no GRCh38-only condition. This
  is the definition E8c used for the comparator analysis ("2–108 INDELs", E5b §7).
- **spelled**: the ALT is spelled by a graph path, whatever the reads do (E5b §5).
- **broad (E5c)**: a tensor allele equals the truth on GRCh38-only (representative or not) and none does on the HPRC
  graph. The net figure subtracts the HPRC-only gains.

### 1b. chr1–22, all truth records (the E5b denominators): COMPUTED NOW from EXISTS inputs (`D/absorption_denominators.md` §"chr1-22 all records")

| | HG008T SNV | HG008T INDEL | COLO829T SNV | COLO829T INDEL |
|---|---|---|---|---|
| truth n | 8,690 | 8,496 | 42,035 | 1,912 |
| population-tagged | 910 (10.5 %) | 5,326 (62.7 %) | 1,442 (3.4 %) | 1,005 (52.6 %) |
| **v2.1 strict**, ≥ 1 platform: absorbed / tagged / NOT tagged | 140 (1.6 %) / 140 / 0 | 1,467 (17.3 %) / 1,418 / **49 (0.58 %)** | 54 (0.13 %) / 54 / 0 | 347 (18.1 %) / 341 / **6 (0.31 %)** |
| v2.1 strict per platform: NOT tagged (% truth) | 0 / 0 / 0 | 38 / 31 / 19 (0.45 / 0.36 / 0.22 %) of 1,301 / 839 / 472 absorbed | 0 / 0 / 0 | 5 / 4 / 2 (0.26 / 0.21 / 0.10 %) of 315 / 198 / 128 |
| v2.1 strict, all 3 platforms | 91 / 0 untagged | 345 / 12 untagged (0.14 %) | 35 / 0 | 99 / 2 (0.10 %) |
| v2.1 arm-only, ≥ 1 platform | 141 / 0 untagged | 1,608 (18.9 %) / 61 untagged (0.72 %) | 57 / 0 | 365 (19.1 %) / 6 (0.31 %) |
| v2.1 arm-only per platform: absorbed (NOT tagged) | 111 / 128 / 131 (0) | 1,576 / 998 / 833 (58 / 41 / 36 = 0.68 / 0.48 / 0.42 %) | 48 / 55 / 54 (0) | 363 / 220 / 172 (6 / 4 / 3) |
| v2.1 ALT spelled by graph (any reads) | 153 (0 untagged) | 1,809 (21.3 %; 66 untagged, 0.78 %) | 70 (1 untagged) | 398 (20.8 %; 8 untagged, 0.42 %) |
| **v1.1 strict**, ≥ 1 platform: absorbed / NOT tagged | 184 / 0 | 1,818 (21.4 %) / **55 (0.65 %)** | 64 / 1 | 382 (20.0 %) / **5 (0.26 %)** |
| v1.1 arm-only, ≥ 1 platform | 188 / 1 untagged | 2,039 / 80 untagged (0.94 %) | 66 / 1 | 402 / 5 (0.26 %) |

Per platform = Illumina / PacBio (Fiber-seq for COLO829T) / ONT. Over chr1–22 inside the truth BED, the HG008T
denominators are 8,543 SNVs and 7,982 INDELs (61.96 % tagged). These are in `D/absorption_denominators.md` §"chr1-22
in BED"; the COLO829T truth lies entirely in its BED.

### 1c. chr1 inside the benchmark BED (the vcfeval scope): COMPUTED NOW (`D/absorption_denominators.md` §"chr1 in BED")

| | HG008T SNV | HG008T INDEL | COLO829T SNV | COLO829T INDEL |
|---|---|---|---|---|
| truth n (rtg baseline) | 697 (697) | 588 (rtg 587) | 2,443 | 131 |
| tagged | 60 (8.6 %) | 384 (65.3 %) | 84 (3.4 %) | 81 (61.8 %) |
| v2.1 strict, per platform: absorbed (NOT tagged) | 10 / 10 / 7 (0) | 93 / 66 / 39 (2 / 2 / 2 = 0.34 % each) | 3 / 3 / 2 (0) | 23 / 21 / 18 (1 / 1 / 1 = 0.76 %) |
| v2.1 strict ≥ 1 platform | 11 (0) | 107 (3 untagged, 0.51 %) | 3 (0) | 31 (1, 0.76 %) |
| v2.1 arm-only per platform (= E8c sets) | 10 / 10 / 10 (0) | 108 / 78 / 64 (4 / 3 / 2) | 3 / 3 / 3 (0) | 34 / 26 / 21 (1 / 1 / 1) |
| v1.1 strict per platform | 11 / 11 / 8 (0) | 114 / 64 / 38 (0 / 1 / 1) | 3 / 3 / 2 (0) | 26 / 21 / 16 (1 / 1 / 0) |

The letter's "1–2 INDELs per dataset in the scored region" is the strict count after the GRCh38-only loss (E5b §8:
2 / 2 / 2 and 1 / 1 / 1), and it **EXISTS** there. On the arm-only (E8c) definition the count is **1–4**.

### 1d. Reconciling the broad E5c loss with the strict count on one denominator: COMPUTED NOW (`D/absorption_denominators.md`, broad table)

HG008T INDEL, v2.1, chr1–22 (n = 8,496), per platform Illumina / PacBio / ONT:

| step | Illumina | PacBio | ONT |
|---|---|---|---|
| broad loss vs GRCh38-only (any tensor allele) | 3,528 (41.5 %) | 3,577 (42.1 %) | 2,494 (29.4 %) |
| of which tagged anyway | 3,013 (85.4 %) | 2,974 (83.1 %) | 2,139 (85.8 %) |
| broad loss, NOT tagged (% truth) | 515 (6.06 %) | 603 (7.10 %) | 355 (4.18 %) |
| minus HPRC-only gains (NOT tagged) | 21 | 8 | 59 |
| **net NOT-tagged loss** (E5c's 3.5–7.0 %) | 494 (5.81 %) | 595 (7.00 %) | 296 (3.48 %) |
| of the 515 / 603 / 355: strict absorbed | 38 | 31 | 19 |
| E5c class I1 (ALT reads walk graph alt nodes, no edit) | 44 | 58 | 32 |
| E5c class I2 (alt nodes plus a residual mismatch / indel) | 414 | 483 | 290 |

COLO829T INDEL, v2.1 (n = 1,912):

- Broad loss: 612 / 567 / 457 (32.0 / 29.7 / 23.9 %).
- Not tagged: 118 / 133 / 110 (6.2 / 7.0 / 5.8 %). Net: 114 / 132 / 103 (5.96 / 6.90 / 5.39 %).
- Of these, strict absorbed 5 / 4 / 2; I1 11 / 11 / 11; I2 67 / 66 / 55.

On v1.1 the net NOT-tagged loss is HG008T 6.16 / 7.23 / 3.83 % and COLO829T 6.59 / 6.80 / 5.02 %, so **3.8–7.2 %**.

Reading the two together. On one denominator, the graph-specific candidate loss that survives filtering is 3.5–7.0 %
of truth INDELs (v2.1). About 89–91 % of it is reads that follow a graph allele (I1 + I2 on HG008T). Only 7–9 % of it is
entire absorption (I1), and the strict spelling test catches 0.22–0.45 % of truth per platform (0.58 % as a union).
The letter's "absorption explains only a small part" is therefore true only for entire, spelled absorption. Partial
graph-path alignment (I2) is the bulk of the post-filter loss (E5c README: I2 80–82 % on HG008T).

The letter quotes two "net" figures:

- R1.2: 3.5–7.0 % (E5c: any tensor allele, v2.1).
- R1.4: 1.7–7.7 % (E5b §1: representative tensor; d9 and d46 together; d46 alone 1.7–7.6 %).

Both **EXIST**. They differ only in whether a non-representative tensor allele counts as a tensor.

SNV, v2.1, chr1–22: net NOT-tagged loss 0.21–0.35 % (HG008T) and 0.21–0.65 % (COLO829T). Strict absorbed SNVs are all
tagged (0 untagged).

### 1e. Comparators at the absorbed alleles (chr1, BED; v2.1 arm-only definition): EXISTS

Source: `analysis/reviewer_response_20261011/experiments/E8c_comparators_absorbed/tables_d46.md` and `summary_d46.json`.

- Absorbed per dataset: INDEL 108 / 78 / 64 (HG008T), 34 / 26 / 21 (COLO829T), 4 / 3 / 2 (HCC1395T); SNV 10 / 10 / 10,
  3 / 3 / 3, 4 / 3 / 4.
- Any record emitted, absorbed INDELs: DeepSomatic-TO 36–100 %, ClairS-TO 3–75 %.
- Native PASS: 9 of 780 allele–caller pairs, all DeepSomatic-TO INDELs (HG008T Illumina 5, ONT 1; COLO829T Illumina 1;
  HCC1395T Illumina 2).
- After the common filter: 1 pair (DeepSomatic-TO, HG008T ONT).
- No absorbed SNV is PASS for either caller.
- v1.1 (read-level "perfect bypass" set): 7 of 882 pairs PASS, 0 after filtering (`E8c/tables.md`).

---

## 2. R1.5 / R2.1: the HapMap mixture

Scope: chr1, v1.6 benchmark BED, truth VAF 1–16.5 %, atomized graph alleles of `mgi_non_HG005.d9`, Pansoma 200-row
models trained 10-24 (SNV `HapMap_Illumina_SNV_full`, INDEL `HapMap_Illumina_INDEL`), PASS = validation recall-0.9,
**no population filtering**.

Sources:

- Pansoma, EXISTS: `experiments/E23_hapmap_lowvaf/atomized/call_vaf_by_value.HapMap_Illumina_{SNV_full,INDEL}.md`
  (full BED) and `E23/mutect2/call_vaf_by_value.Pansoma_*.md` (eval BED without graph-carried sites, job 461100).
- Mutect2 on the eval BED, EXISTS: `E23/mutect2/call_vaf_by_value.Mutect2_*.md`.
- Mutect2 on the full BED, COMPUTED NOW: `D/hapmap_m2_full/`.

| | SNV | INDEL |
|---|---|---|
| original truth (chr1, BED, VAF 1–16.5 %) | 392,567 | 124,499 |
| graph-carried (excluded) | 175,773 (44.8 %) | 25,938 (20.8 %) |
| remaining, full-BED split | 216,794 | 98,561 |
| remaining in `eval_chr1.bed` (graph-carried sites ±1 bp masked; the letter's evaluation) | 215,094 | 92,079 |
| **Pansoma recall, conditional** (eval BED) | 0.548 (117,834 / 215,094) | 0.209 (19,219 / 92,079) |
| Pansoma recall, conditional (full-BED split) | 0.546 (118,308 / 216,794) | 0.204 (20,118 / 98,561) |
| **Pansoma recall, full truth** (graph-carried = FN) | **0.302** (118,403 / 392,567) | **0.166** (20,643 / 124,499) |
| Pansoma precision, calls with AF 1–16.5 %: eval BED / full BED | 0.990 / 0.990 | 0.784 / 0.776 |
| Pansoma whole PASS call set, all VAFs (rtg): eval BED / full BED | P 0.990, R 0.491 / P 0.990, R 0.277 | P 0.780, R 0.193 / P 0.772, R 0.155 |
| Mutect2 PASS recall: conditional (eval BED) / full truth | 0.581 / **0.630** | 0.117 / **0.150** |
| Mutect2 PASS + clustered_events / haplotype: conditional / full | 0.691 / 0.764 | 0.157 / 0.197 |
| Mutect2 PASS precision, calls with AF 1–16.5 %: eval / full | 0.975 / 0.986 | 0.860 / 0.900 |
| Mutect2 PASS whole call set (rtg): eval / full | P 0.362, R 0.523 / P 0.527, R 0.580 | P 0.266, R 0.108 / P 0.379, R 0.139 |

The last rows are the summary.txt of each vcfeval; Mutect2's low whole-set precision comes from calls outside the
1–16.5 % AF range.

Per-bin recall, conditional (EXISTS, used in R2.4):

- Pansoma SNV: 0.302 / 0.807 / 0.906 / 0.796 at 1–2 / 2–5 / 5–10 / 10–16.5 %.
- Pansoma SNV on the full truth (COMPUTED from the EXISTS atomized file): 0.225 / 0.232 / 0.466 / 0.105.

Answers:

- **Were the comparators scored on the same denominator?** Yes, for Mutect2 only, on chr1. The eval-BED scoring uses
  the same `eval_chr1.bed`, the same truth and the same rtg options (`tmp/mutect2_hapmap_20261029/evaluate.sh`).
  - On the conditional denominator, Mutect2 PASS has the higher SNV recall (0.581 vs 0.548) and the lower AF-binned
    precision (0.975 vs 0.990).
  - On the full truth, the SNV recall gap widens to 0.630 vs 0.302: a linear caller sees the graph-carried alleles as
    mismatches.
  - INDEL: Pansoma 0.209 vs Mutect2 0.117 conditional; 0.166 vs 0.150 on the full truth.
  - No ClairS-TO or DeepSomatic run on HapMap was found (MISSING).
- **Was the HapMap truth pre-filtered by a population resource?** No such step is documented.
  - Truth = SMHTHAPMAP6 v1.6 without HG005 (qfu), the minor donors' alleles absent from HG005
    (`/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/README.txt`). Its source is Kong et al., a graph-based truth
    built from assemblies (E8b M2).
  - The Pansoma evaluation applies no population filter to the calls (E23 README; R2.1 notes of 10-29).
  - If the filter were applied, it would tag 95.7 % of the chr1 SNV truth records (411,638 / 430,009) and 42.2 % of
    the INDEL records (58,662 / 138,854). Source: `/scratch/jshen/data/pansoma_net_v2_runs/HapMap_Illumina_*/vcf_chr1/eval_bed/report.txt`, EXISTS.
  - On chr20 positives it would tag SNV 96.3 % (125,026 / 129,857) and INDEL 40.4 % on the new rule (E8b task 6).
  - After the filter, recall is SNV 0.007 and INDEL 0.021 (E23 README Result 3).
- **Companion paper (Fu et al., 2026) "little difference across 1–16.5 %, clear drop below 1 %":** MISSING.
  - The only repo source is E8b M4, which gives the paper's existence and its abstract (HapMap mixture + COLO829
    benchmarks). It has no VAF statement.
  - R1.4's evidence records a 10-20 WebFetch of the Methods (surjection to GRCh38; Mutect2 tumor-only on HapMap), again
    nothing on VAF.
  - My one fetch of biorxiv `10.64898/2026.01.04.697580v2.full` today returned HTTP 429. The [待核实] stays.
- **The original submission's HapMap truth:** E23 README Result 6 (EXISTS) says it was a fixed subset of 96,393 truth
  SNVs whose origin and script are unknown, scored for precision against a larger truth. Whether that subset was
  population-filtered: MISSING.

---

## 3. R1.6: one table, chr1, benchmark BED, final models

Models:

- SNV: `--norm dataset` v2.1 d46, three seeds (`<S>_<P>_SNV_d46_normbq{,_s1,_s2}`).
- INDEL: final release models, ppr8 operating threshold, three seeds (`release_s{0,1,2}_<P>_INDEL/vcf_ppr8_chr1_*_d46`).
- Comparators: E8 native PASS → PASS + the common filter.

The recall rows are the rtg baseline (n = 697 / 587 / 2,443 / 131); they EXIST in
`experiments/E8_comparators/latest_models.md` and in `/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score/*/report.json`.
The candidate rows are COMPUTED NOW (`D/label_status.json`, chr1, PASS truth in BED, n = 697 / 588 / 2,443 / 131).

| chr1, BED (Illumina / PacBio or Fiber-seq / ONT) | HG008T SNV | COLO829T SNV | HG008T INDEL | COLO829T INDEL |
|---|---|---|---|---|
| max recall under the filter policy (untagged / truth, in BED) | 0.914 (637 / 697) | 0.966 (2,359 / 2,443) | **0.347** (204 / 588) | 0.382 (50 / 131) |
| same, whole chr1 incl. outside BED (the letter's numbers; EXISTS, report.json `truth_pon_tags`) | 0.913 (641 / 702) | 0.966 | 0.342 (217 / 635) | 0.382 |
| candidate, **exact allele** | 0.956 / 0.950 / 0.956 | 0.820 / 0.822 / 0.819 | **0.432 / 0.480 / 0.417** | **0.466 / 0.519 / 0.427** |
| candidate, **partial only** (of which off-reference) | 0.003 / 0.003 / 0.001 (all off-ref) | ≤ 0.001 | 0.180 / 0.279 / 0.310 (65 / 110 / 101 alleles off-ref) | 0.107 / 0.191 / 0.244 (10 / 12 / 13 off-ref) |
| exact + partial (= letter's row "reaches a candidate"; predict's `with_tensor`) | 0.958 / 0.953 / 0.957 (predict 0.954 / 0.948 / 0.954) | 0.821 / 0.822 / 0.819 (predict 0.815 / 0.817 / 0.812) | 0.612 / 0.759 / 0.726 (predict 0.611 / 0.757 / 0.724) | 0.573 / 0.710 / 0.672 (predict 0.565 / 0.695 / 0.649) |
| exact candidate AND not tagged (attainable after filtering) | 0.898 / 0.894 / 0.904 | 0.795 / 0.797 / 0.795 | 0.259 / 0.259 / 0.247 | 0.275 / 0.282 / 0.260 |
| Pansoma PASS, **before** filtering | 0.847–0.887 / 0.841–0.874 / 0.849–0.900 | 0.722–0.735 / 0.723–0.751 / 0.716–0.734 | 0.111–0.152 / 0.164–0.191 / 0.165–0.184 | 0.130–0.176 / 0.107–0.137 / 0.160–0.168 |
| Pansoma PASS, **after** filtering | 0.814–0.852 / 0.814–0.838 / 0.824–0.859 | 0.705–0.717 / 0.707–0.734 / 0.700–0.718 | 0.111–0.148 / 0.160–0.177 / 0.150–0.162 | 0.130–0.145 / 0.107–0.137 / 0.160 |
| Pansoma post-filter curve maximum | 0.891–0.895 / 0.889–0.891 / 0.897–0.900 | 0.793 / 0.794–0.795 / 0.788–0.790 | 0.254–0.256 / 0.256–0.259 / 0.230–0.237 | 0.275 / 0.260–0.267 / 0.221–0.244 |
| DeepSomatic-TO native PASS → after filter | 0.911→0.889 / 0.816→0.806 / 0.891→0.869 | 0.846→0.835 / 0.562→0.555 / 0.705→0.695 | 0.174→0.141 / 0.036→0.031 / 0.019→0.019 | 0.214→0.183 / 0.092→0.084 / 0.015→0.015 |
| ClairS-TO native PASS → after filter | 0.859→0.845 / 0.857→0.844 / 0.845→0.831 | 0.774→0.764 / 0.753→0.745 / 0.692→0.684 | 0.080→0.080 / 0.075→0.075 / 0.075→0.075 | 0.145→0.145 / 0.107→0.107 / 0.092→0.092 |

Notes:

- The post-filter curve maximum is bounded by "exact candidate AND not tagged" in every dataset. It equals it on
  COLO829T Illumina INDEL (0.275 = 36 / 131).
- The INDEL ceiling is therefore set by exact candidate generation plus the filter, not by the classifier.
- Partial-only INDEL matches cannot be credited by rtg. They are a different allele, and the off-reference ones are
  graph-coordinate calls.
- Predict's `with_tensor` differs from my count by 1–6 alleles because predict leaves out tensors whose position test
  falls outside the BED.
- Same exact-only comparison on GRCh38-only and v1.1 (COMPUTED NOW, `D/label_status.{linear,d9}.json`). HG008T INDEL
  exact candidate:
  - GRCh38-only 0.804 / 0.842 / 0.587 (with partial 0.828 / 0.874 / 0.750: the letter's "0.75–0.87");
  - v1.1 0.401 / 0.471 / 0.440;
  - v2.1 0.417–0.480.
- Exact and untagged, HG008T INDEL:
  - GRCh38-only 168 / 178 / 155 = 0.286 / 0.303 / 0.264;
  - v2.1 152 / 152 / 145. GRCh38-only has 16 / 26 / 10 more callable alleles (= E5b §8).
  - COLO829T: GRCh38-only exact 0.809–0.939, exact untagged 0.336–0.382; v2.1 0.260–0.282.
- SNV exact candidate: GRCh38-only 0.977–0.981 (HG008T) and 0.822–0.833 (COLO829T).
- Rescue: SNV `vcf_chr1` reports were made with predict's off-reference rescue on (graph_vcf "call" string). E3b shows
  that rtg TP / FP / FN are identical without it; only graph-coordinate calls change.

---

## 4. R2.6: HG008T with and without population filtering, chr1, BED

COMPUTED NOW from EXISTS report.json files (`D/r26_filter_table.md`). Cells: range over seeds (mean); P / R / F1 at PASS.

The GRCh38-only arm has only `--norm train` 7fc7396 runs. Its seed-0 values and the seed-0 matched d46 / d9 values are
read from `analysis/benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/runs/*/vcf_chr1.eval_bed.txt`, because
their report.json files are not local.

| | graph / model | without filter P / R / F1 | with filter P / R / F1 |
|---|---|---|---|
| SNV Illumina | v2.1 final (3 seeds) | 0.176–0.206 / 0.847–0.887 / 0.294–0.332 | 0.426–0.451 / 0.814–0.852 / 0.566–0.585 |
| | v1.1 final (3) | 0.104–0.108 / 0.834–0.862 / 0.185–0.192 | 0.456–0.461 / 0.803–0.831 / 0.586–0.588 |
| | GRCh38-only 7fc7396 (3) | 0.005 / 0.836–0.839 / 0.009–0.010 | 0.334–0.394 / 0.783–0.789 / 0.469–0.524 |
| SNV PacBio | v2.1 final | 0.186–0.222 / 0.841–0.874 / 0.307–0.351 | 0.471–0.487 / 0.814–0.838 / 0.603–0.609 |
| | v1.1 final | 0.094–0.102 / 0.845–0.862 / 0.169–0.182 | 0.422–0.484 / 0.811–0.829 / 0.555–0.611 |
| | GRCh38-only | 0.005 / 0.845–0.882 / 0.009–0.010 | 0.372–0.472 / 0.792–0.826 / 0.513–0.592 |
| SNV ONT | v2.1 final | 0.171–0.209 / 0.849–0.900 / 0.287–0.335 | 0.433–0.481 / 0.824–0.859 / 0.576–0.607 |
| | v1.1 final | 0.090–0.097 / 0.845–0.862 / 0.163–0.175 | 0.443–0.472 / 0.814–0.831 / 0.575–0.598 |
| | GRCh38-only | 0.004–0.005 / 0.835–0.892 / 0.008–0.009 | 0.368–0.394 / 0.780–0.842 / 0.509–0.537 |
| INDEL Illumina | v2.1 final release, ppr8 (3) | 0.133–0.460 / 0.111–0.152 / 0.141–0.224 | 0.376–0.460 / 0.111–0.148 / 0.171–0.224 |
| | v2.1 single-set normbq, r0.9 (1) | 0.071 / 0.375 / 0.119 | 0.267 / 0.240 / 0.253 |
| | v1.1 single-set normbq, r0.9 (1) | 0.059 / 0.354 / 0.101 | 0.265 / 0.237 / 0.250 |
| | GRCh38-only 7fc7396, r0.9 (1) | 0.011 / 0.738 / 0.021 | 0.102 / 0.259 / 0.146 |
| INDEL PacBio | v2.1 final release (3) | 0.433–0.504 / 0.164–0.191 / 0.241–0.277 | 0.448–0.510 / 0.160–0.177 / 0.239–0.263 |
| | v2.1 / v1.1 / GRCh38-only 7fc7396, r0.9 | 0.025 / 0.431 / 0.046; 0.021 / 0.431 / 0.041; 0.012 / 0.755 / 0.024 | 0.061 / 0.256 / 0.098; 0.056 / 0.256 / 0.092; 0.090 / 0.256 / 0.133 |
| INDEL ONT | v2.1 final release (3) | 0.100–0.151 / 0.165–0.184 / 0.130–0.158 | 0.153–0.212 / 0.150–0.162 / 0.157–0.178 |
| | v2.1 / v1.1 / GRCh38-only 7fc7396, r0.9 | 0.005 / 0.385 / 0.011; 0.006 / 0.400 / 0.011; 0.005 / 0.538 / 0.009 | 0.008 / 0.237 / 0.015; 0.008 / 0.228 / 0.015; 0.012 / 0.235 / 0.023 |

Test-tuned best F1 without → with filter is in `D/r26_filter_table.md`. For example, SNV v2.1 final 0.398–0.419 →
0.601–0.623; GRCh38-only 0.034–0.050 → 0.544–0.604. **MISSING:** final INDEL models on v1.1 or GRCh38-only; they were
trained and evaluated on d46 only.

The letter's point 3 (HG008T SNV v1.1, F1 0.182–0.191 → 0.585–0.621; INDEL 0.011–0.101 → 0.015–0.251; best
0.090–0.150 → 0.144–0.270) is the seed-0 7fc7396 benchmark. It EXISTS in
`analysis/benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/tables.md` and R2.6 evidence 10-16. It is not the
final three-seed models and not v2.1 / GRCh38-only.

**R2.6 point 4 ("[待补: the authors will revise this paragraph.]") facts:**

1. HapMap donors and the resources: HG00438, HG02257, HG02486 and HG02622 are 1000 Genomes samples, and CoLoRSdb
   includes 127 HPRC genomes (E8b P11). The donors' rare chr20 SNVs are found in CoLoRSdb 14.5–25.9 % of the time and
   in gnomAD 3.4–15.4 %, against 4.0–7.9 % and 0.6–2.5 % for the cancer normals (E8b §4.4). EXISTS.
2. Pre-filtering vs post-filtering: the truth is not pre-filtered by any population resource (§2 above), and the
   Pansoma (and Mutect2) HapMap evaluation is **without** population filtering of calls (E23). The current sentence
   "The HapMap positives were already filtered in that benchmark (R1.5)" has no source. R2.6 evidence line 29 records it
   as the user's phrasing, and line 130 notes it is inconsistent with E23.
3. What filtering would do: it tags 95.7 % of chr1 SNV truth (411,638 / 430,009) and 42.2 % of chr1 INDEL truth
   (58,662 / 138,854); chr20 positives 96.3 % SNV and 40.4 % INDEL. Post-filter recall is SNV 0.007 and INDEL 0.021.
   EXISTS.
4. Whether the original submission's HapMap truth subset (96,393 SNVs) was population-filtered: MISSING (E23 Result 6:
   its origin is unknown).

---

## 5. R1.8: label rules and counts

### 5a. Rules: EXISTS

Source: `indexed_gam_pipeline_v4/tensor_postprocessing/README.md` §Labels and `truth_labels.py` `classify()`. The labels
of every HG008T / COLO829T set have `rules_sha256 5430533aa100` (commit 1e9ff6d). The current file (663a7ff2da37,
a01969d) is a cleanup with byte-identical labels; 17b540137dc3 (4f36af1) only adds `--ignore-repeat-indels`, which is
used for SEQC2 / CASTLE.

Key equality: a truth INDEL is enumerated at up to 5,000 equivalent placements, and node keys are made only on nodes the
GRCh38 path visits once (either orientation). The first check that applies decides:

| # | check | label | reason |
|---|---|---|---|
| 1 | AF floor (unset for these sets) | −1 | below_*_min_af |
| 2 | A1 (representative allele) key = somatic truth, PASS (inside **or outside** the BEDs) / only filtered | 1 / −1 | representative_allele_is_somatic_truth / somatic_truth_filtered |
| 3 | A1 = germline truth, PASS or GAP1/2 (inside or outside BEDs) | 2 | representative_allele_is_germline_truth(_gap_filtered) |
| 4 | another allele of the site = truth: A1 overlaps a somatic one > 45 % → 1; else a germline one > 45 % → 2; else 0 | 1 / 2 / 0 | allele_partial_somatic_truth / allele_partial_germline_truth / truth_matches_non_representative_allele |
| 5 | another allele = filtered somatic truth | −1 | somatic_truth_filtered |
| 6 | no GRCh38 position and no anchor | −1 | not_on_unique_grch38_node |
| 7 | outside confident region (somatic BED ∩ germline BED) | −1 | outside_confident_region |
| 8 | A1 = germline truth with another FILTER | 0 | germline_truth_filtered |
| 9 | PASS somatic truth at the same place, allele overlap > 45 % **or haplotype overlap** (read refit in the truth's GRCh38 window) > 45 % | 1 | residual_partial_somatic_truth |
| 10 | off-reference (anchored) tensor, nothing above | −1 | off_reference_no_truth_match |
| 11 | PASS/GAP germline truth at the same place, allele overlap > 45 % (no haplotype overlap; not tried off-reference) | 2 | residual_partial_germline_truth |
| 12 | any truth within 10 bp | 0 | near_truth_allele_mismatch |
| 13 | otherwise | 0 | confident_no_truth_allele |

Somatic is tried before germline. What training does with these labels:

- Every final model (SNV normbq; INDEL release) trains with `--ignore-reasons residual_partial_somatic_truth
  allele_partial_somatic_truth` (`args.json`). Partial somatic matches therefore become −1 in training **and**
  validation; partial germline matches stay germline.
- Validation and chr1 test score `off_reference_no_truth_match` as non (`eval_as_non` in the index meta; predict.py
  docstring).
- The final INDEL models also use `--non-reasons-in-bed outside_confident_region`, which trains and validates
  outside-confident-region tensors inside the somatic BED as 0, and `--non-reasons repeat_indel_without_truth`
  (SEQC2 / CASTLE). The SNV models use neither.
- INDEL tests (predict) count partial somatic matches as found truth; SNV tests leave them out.

**Off-GRCh38-path tensors, which statement is true.**

- An off-reference tensor can **never** be an exact somatic match: there are 0 exact somatic labels off-reference in
  all 12 v2.1 set × kind (and in the d9 / linear sets). Keys exist only on GRCh38 nodes.
- It **can** be labelled somatic through step 9. It is placed by node-ID anchoring to the flanking reference nodes
  (within 200 IDs, ≤ 1,024 GRCh38 bp), and the haplotype overlap of its reads against the truth allele's GRCh38 window
  can exceed 45 %.
- Seed-0 training blocks hold 461 / 235 / 225 (HG008T SNV), 840 / 1,372 / 1,248 (HG008T INDEL), 97 / 71 / 66
  (COLO829T SNV) and 94 / 209 / 231 (COLO829T INDEL) such tensors. They are all excluded from training by
  `--ignore-reasons`.
- Off-reference INDEL tensors **are** training negatives through step 4 (label 0, truth_matches_non_representative_allele):
  119 / 83 / 332 (HG008T) and 101 / 90 / 717 (COLO829T), seed-0 training blocks. SNV: 0.
- Off-reference tensors without a match are validation negatives: SNV 1,748–4,059 and INDEL 838–12,117 per set (E21
  `set_split_counts.tsv` `val_offref_as_non`).
- So R3.4's "93.6–99.7 % of the somatic examples" counts chr1 test tensors with truth label 1, **including** step-9
  partial matches that are never training or SNV test examples. Of the somatic **training** examples, 100 % are on
  GRCh38-path nodes.

### 5b. Counts per training set (v2.1, chr2–22): COMPUTED NOW (`D/label_status.json`)

Truth = PASS somatic in the truth BED, chr2–22. "Training region" = outside the validation blocks of the seed-0 split;
seeds 1 / 2 differ by ≤ 0.4 percentage points (json `seeds`).

| set | kind | truth chr2–22 | training region (seed 0) | **exact training example** | partial only (excluded) | no tensor: filtered / no candidate / other allele of site (all chr2–22) | PASS truth outside BED with an exact positive |
|---|---|---:|---:|---|---|---|---:|
| HG008T Illumina | SNV | 7,846 | 7,467 | 7,160 (95.9 %) | 28 (0.4 %) | 67 / 227 / 2 | 115 |
| HG008T PacBio | SNV | 7,846 | 7,465 | 7,111 (95.3 %) | 16 | 43 / 318 / 2 | 115 |
| HG008T ONT | SNV | 7,846 | 7,464 | 7,141 (95.7 %) | 17 | 52 / 272 / 2 | 117 |
| COLO829T Illumina | SNV | 39,592 | 37,679 | 34,095 (90.5 %) | 20 | 3,075 / 637 / 0 | 0 |
| COLO829T Fiber-seq | SNV | 39,592 | 37,693 | 33,982 (90.2 %) | 13 | 3,281 / 568 / 0 | 0 |
| COLO829T ONT | SNV | 39,592 | 37,676 | 33,760 (89.6 %) | 14 | 3,597 / 467 / 0 | 0 |
| HG008T Illumina | INDEL | 7,394 | 7,021 | 3,213 (45.8 %) | 1,114 (15.9 %) | 507 / 2,273 / 47 | 161 |
| HG008T PacBio | INDEL | 7,394 | 7,022 | 3,546 (50.5 %) | 1,761 (25.1 %) | 460 / 1,069 / 275 | 188 |
| HG008T ONT | INDEL | 7,394 | 7,024 | 3,066 (43.7 %) | 1,966 (28.0 %) | 622 / 1,009 / 462 | 155 |
| COLO829T Illumina | INDEL | 1,781 | 1,702 | 1,037 (60.9 %) | 92 (5.4 %) | 74 / 501 / 10 | 0 |
| COLO829T Fiber-seq | INDEL | 1,781 | 1,701 | 1,134 (66.7 %) | 185 (10.9 %) | 120 / 238 / 32 | 0 |
| COLO829T ONT | INDEL | 1,781 | 1,701 | 1,007 (59.2 %) | 300 (17.6 %) | 113 / 176 / 118 | 0 |

The partial-only INDEL alleles on chr2–22 split into GRCh38-node and off-reference matches:

- HG008T: 518 / 662, 758 / 1,100, 1,043 / 1,034 (Illumina, PacBio, ONT);
- COLO829T: 53 / 49, 89 / 110, 189 / 125.

"Not labelable" truth in the BED: none. Every exact match is labelled 1. The non-training categories are exactly
partial-only, filtered, no candidate and other allele of the site. Somatic positives also come from PASS truth
**outside** the BED (last column), which the letter does not mention.

Training tensors (seed-0 training blocks; E21 matches):

| set | somatic exact | somatic partial left out (on GRCh38 / off-ref) | germline exact | **germline partial** (% of germline) |
|---|---:|---|---:|---|
| HG008T Illumina SNV / INDEL | 7,274 / 3,383 | 1,207 / 461; 677 / 840 | 256,583 / 46,929 | 0 / 8,331 (15.1 %) |
| HG008T PacBio SNV / INDEL | 7,227 / 3,761 | 229 / 235; 938 / 1,372 | 270,175 / 128,299 | 0 / 34,590 (21.2 %) |
| HG008T ONT SNV / INDEL | 7,259 / 3,249 | 215 / 225; 1,220 / 1,248 | 302,622 / 128,058 | 0 / 47,983 (27.3 %) |
| COLO829T Illumina SNV / INDEL | 34,110 / 1,048 | 259 / 97; 79 / 94 | 280,542 / 51,400 | 0 / 9,710 (15.9 %) |
| COLO829T Fiber-seq SNV / INDEL | 33,998 / 1,141 | 112 / 71; 140 / 209 | 320,431 / 140,716 | 0 / 40,089 (22.2 %) |
| COLO829T ONT SNV / INDEL | 33,779 / 1,018 | 134 / 66; 284 / 231 | 355,824 / 154,132 | 0 / 58,297 (27.4 %) |

Germline partial matches are mostly `residual_partial_germline_truth`, for example 31,107 of 34,590 on HG008T PacBio.

Outside-confident-region INDEL tensors trained as 0 by the final INDEL models (inside the somatic BED, seed-0 training
blocks):

- HG008T: 2,512 / 12,027 / 17,932 of 21,522 / 100,267 / 124,062 (12–14 %).
- COLO829T: 13,299 / 157,899 / 122,958 (100 %; the COLO829T somatic BED, SMaHT union, contains all of them).

---

## 6. Statements in the current letter that conflict with these facts

Letter: `/tmp/claude-10008/-scratch-jshen-Github-Pansoma/2a27ec06-fc25-4d6b-b0b5-2fbc7483d580/scratchpad/assemble/response_letter_20261015.md`,
identical to the repo `points/*.md`.

**R1.4**

1. "1,467 HG008T and 347 COLO829T truth INDELs (and 140 and 54 truth SNVs) get no candidate on at least one platform
   although a graph path spells their ALT"
   - These counts also require a tensor on the GRCh38-only graph.
   - Without that condition: 1,608 / 365 INDELs and 141 / 57 SNVs.
   - Give denominators: 17.3 % of 8,496 and 18.1 % of 1,912.
2. "At the alleles the HPRC v2.1 graph absorbs on chr1 (2–108 INDELs ...)" uses the arm-only definition, not the one
   behind 1,467. The two sentences use different definitions; on chr1 the strict counts are 93 / 66 / 39 and 23 / 21 / 18.
3. "1–2 INDELs per dataset in the scored region" holds for the strict definition only; on the E8c definition it is 1–4.
4. "HPRC v1.1 ... similar (1.2% and 0.5% by a read-level test)" uses another definition and an ALT-split denominator
   (9,362). The same strict test on v1.1 gives **0.65 % (55 / 8,496) and 0.26 % (5 / 1,912)**.
5. "absorption explains only a small part of this (38 / 31 / 19 of the 446 / 542 / 268)" is true only for strict
   entire absorption.
   - On the E5c definition, 89–91 % of the post-filter graph-specific loss (HG008T) is reads following a graph allele,
     mostly with a residual edit (I2: 414 / 483 / 290 of 515 / 603 / 355).
   - This contradicts the framing of R1.2 ("The concern holds for indels ... net loss 3.5–7.0 %").
6. The "1.7–7.7%" net INDEL loss (R1.4, E5b representative tensors) and R1.2's "3.5–7.0%" (E5c any tensor allele)
   measure the same thing two ways. Pick one: on v2.1, E5c net = 3.5–7.0 %; E5b = 1.7–7.6 %.
7. "After filtering, 24–27% of truth INDELs are called on every reference" comes from the E1.6 / 7fc7396 recall-0.9
   runs. With the final release models (ppr8) the post-filter INDEL recall is **0.111–0.177 (HG008T) and 0.107–0.160
   (COLO829T)**.

**R1.5 / R2.1**

8. Recall is reported only on the conditional denominator; the full-truth figures are missing.
   - Full-truth recall (chr1, 1–16.5 %): SNV **0.302**, INDEL **0.166** (conditional 0.548 / 0.209).
   - Graph-carried share: 44.8 % (SNV) and 20.8 % (INDEL) of the truth.
   - Mutect2 on the same denominators: 0.581 / 0.117 conditional, **0.630 / 0.150** on the full truth.
9. The eval-BED masking (±1 bp) also drops non-graph-carried truth. The remaining denominators are 215,094 SNVs and
   92,079 INDELs, not 216,794 and 98,561.

**R1.6**

10. The "maximum attainable recall 0.342 (HG008T)" and the table row "(whole chr1 truth) 0.913 / 0.966 / 0.342 / 0.382"
    use 635 / 702 whole-chr1 records, while every other row uses the BED baseline (587 / 697). In the BED: **0.347**
    (204 / 588), 0.914, 0.966, 0.382. Likewise "418 of the 635" in the BED is **384 of 588**.
11. The row "truth allele reaches a candidate" for INDELs (0.611–0.757 / 0.565–0.695) includes partial matches.
    - Exact only: **0.417–0.480 (HG008T)** and **0.427–0.519 (COLO829T)**.
    - Partial only: 0.180–0.310 and 0.107–0.244, of which 65–110 and 10–13 alleles are off-reference.
12. "A GRCh38-only graph gives more HG008T truth INDELs a candidate (0.75–0.87 against 0.61–0.76)" includes partial
    matches. Exact only: 0.587–0.842 against 0.417–0.480. In the same bullet, "HG008T: 0.43–0.50" (HCC1395 comparison)
    is the exact genome-wide E5b figure, a different definition from the table row.
13. The INDEL operating-point and comparator numbers check out against the files (0.111–0.177 / 0.107–0.160; curve
    0.230–0.259 / 0.221–0.275; comparators 0.015–0.183 after the filter). No conflict.

**R2.6**

14. Point 3 reports only v1.1 seed-0 7fc7396 runs. With the final models (three seeds):
    - SNV F1 without → with filter: v2.1 0.287–0.351 → 0.566–0.609; v1.1 0.163–0.192 → 0.555–0.611; GRCh38-only
      0.008–0.010 → 0.469–0.592.
    - INDEL v2.1 final release: 0.130–0.277 → 0.157–0.263.
    - Two minor value differences: the letter's 0.251 and 0.144 are 0.250 and 0.1435 in the current reports.
15. Point 4: "The HapMap positives were already filtered in that benchmark" has no source and contradicts E23 (no
    filtering is applied to the calls; the truth is not population-filtered).
16. Point 1, "the benchmark does not overstate its effect", is not supported as a certainty. The private-allele test
    shows low membership (0.6–2.5 % / 4.0–7.9 %); it does not exclude membership.

**R1.8**

17. "A candidate on a node off the GRCh38 path has no GRCh38 truth allele to be compared with and is not used as a
    training example" is wrong on two counts.
    - Anchored off-reference tensors are compared through haplotype overlap and labelled somatic (step 9); these are
      excluded from training only as partial matches.
    - Off-reference INDEL tensors are training negatives (83–717 per set) through step 4.
    - Off-reference tensors without a match are validation and test negatives.
    - Correct statement: no somatic training example and no SNV training example is off the GRCh38 path.
18. "Not used for training are tensors without a truth match outside the confident region": false for the final INDEL
    models. They train such tensors as non inside the somatic BED: 12–14 % of them on HG008T, 100 % on COLO829T.
19. "positives require exact allele identity" is correct, but positives also include PASS truth outside the benchmark
    BED: 115–117 SNVs and 155–188 INDELs on HG008T chr2–22.
20. "In the training region, 94.8–95.9% (HG008T) and 89.5–90.5% (COLO829T) of truth SNVs" is consistent (v2.1 seed 0:
    95.3–95.9 % and 89.6–90.5 %). The INDEL counterparts the reviewer asked for (exact 43.7–50.5 % / 59.2–66.7 %;
    partial-only 15.9–28.0 % / 5.4–17.6 %) are not in the letter.
21. The germline partial-match share (15.1–27.4 % of INDEL germline training examples; 0 for SNVs) is not quantified
    in the letter.
22. R3.4's "93.6–99.7% of the somatic examples" conflicts with R1.8's claim only because it counts chr1 test tensors
    with truth label 1, including partial matches that are not training examples. Of the somatic training examples,
    100 % are on GRCh38 nodes. R2.2's "off-reference nodes without a somatic match are excluded" is true for training
    (−1), but such tensors are validation and test negatives.
