# E4b: what the false positives that survive the PoN are

> **2026-10-16, INDEL PoN rule changed** (user decision: gnomAD / CoLoRSdb exact allele, AF >= 0.01 in homopolymers >= 8 bp, else AF >= 1e-4; was >= 7 bp and 1e-3). The 44 local chr1 INDEL call sets were re-scored in place (`tmp/reviewer_cpu_20261013/pon_hp8_20261016/rescore.sh`, array 418042; the old evaluations are kept as `eval_bed.pon_hp7_1e-3` / `eval_nobed.pon_hp7_1e-3`) and this experiment's collector was rerun (job 418071), so `tables.md` / `summary.json` / TSVs now use the new rule. INDEL numbers quoted in the prose below may still be those of the earlier rule (changes are small: e.g. HG008T Illumina v1.1 PASS F1 0.241 -> 0.251, best 0.265 -> 0.270; HG008T chr1 truth INDELs tagged 413 -> 418 of 635). SNV results are unchanged.


Reviewer points:
- **R3.2**: what the graph removes; the composition of the candidates that remain.
- **R2.5**: where the precision gain comes from.

Shared code: `../call_level_common/calls_lib.py`.

**Status (2026-10-15 18:20, job 406568): every call set that exists locally is processed.** 102 call sets, each
in the BED and without it:

- in-sample, 54 call sets:
  - 15 SNV and 15 INDEL at the run's main threshold (validation recall 0.9): linear and d46 for HG008T Illumina /
    PacBio / ONT and COLO829T Illumina / fiberseq / ONT, d9 for the three HG008T platforms;
  - 15 INDEL `vcf_chr1_valbest` re-thresholds of the same runs;
  - 9 SNV `__norescue` re-calls (the six d46 and three d9 runs);
- joint, 8 call sets: the COLO829T fiberseq joint-purity d9 model on P100 / P75 / P50 / P25, SNV and INDEL. It is
  trained on four purities, so it is not a matched d9 arm;
- cross-sample (`xs__`), 40 call sets: 30 SNV (HG008T and COLO829T d46 models of all three platforms, each on the
  other tumor and on HCC1395T / HCC1954T / H1437T / H2009T) and 10 INDEL (Illumina only).

What is still missing is listed at the end of the results.

## Question

After the GRCh38 PoN, which false positives are left? In particular, does the linear (GRCh38-only) arm's extra FP
come from germline or population alleles that a pangenome would have absorbed?

## Data

- **Call sets**: `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1*`, written by
  `machine_learning/pansoma_net_v2/jobs/vcf.sh`. The rtg FP records of the PASS calls after the PoN are
  `eval_bed/vcfeval/pon_calls/fp.vcf.gz` (in the benchmark BED) and the same under `eval_nobed`.
  - Discovery (`calls_lib.discover`) reads the run names:
    - in-sample `<SAMPLE>_<PLATFORM>_<KIND>_{7fc7396 | d46_7fc7396 | linear_7fc7396}`, plus
      `HG008_Illumina_SNV_7fc7396_ddp4_b1024`;
    - `xs__<model>__<set>`;
    - `<run>__norescue`;
    - `vcf_chr1_valbest`;
    - the joint purity sets `vcf_chr1_<SET>`.
  - The arm comes from the tensor set name (`_d46`, `_linear`, else d9). A call set counts once both evaluations and
    every rtg directory are done.
- **Germline truth of the normal** (from the set's `labels.manifest.json`):
  - HG008-N dipcall: `HG008N_GRCh38_dipcall.dip.vcf.gz`, 379,135 chr1 records with an ALT genotype;
  - COLO829BL dipcall: `dipcall_hg38.dip.vcf.gz`, 374,468;
  - for the cross-sample targets, the CASTLE normals' DeepVariant hybrid QUAL10 calls: HCC1395BL 377,806,
    H1437 389,253, H2009 384,898, HCC1954 390,500.
  - Every record whose genotype carries an ALT is used, of any FILTER.
  - Pansoma is tumor-only; the normal is used only for this analysis.
- **HPRC v1.1 haplotypes**: `/scratch/jshen/data/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz` (44 diploid samples +
  CHM13).
  - `../call_level_common/prep_hprc.py` builds `DATA/hprc_v1.1_chr1_alleles.norm.tsv.gz`.
  - It reads 2,147,805 chr1 records and skips 1,402 ALTs with AC 0. Each remaining ALT is left-aligned and trimmed;
    MNPs are also split into single bases. This gives 2,193,706 allele keys (36 s).
  - The HPRC v2.1 VCF is not local, so the d46 runs are also checked against v1.1.
- **Population resources**: the four PoN VCFs in `/scratch/jshen/data/Pansoma/panel_of_normal_VCFs`.
  - The PoN floors are imported from `scripts/filter_panel_of_normals.py`: SNV 1e-4; INDEL 0.01 in homopolymers
    >= 7 bp, else 0.001.

`DATA` = `/scratch/jshen/data/reviewer_experiments_20261013/call_level`.

## Method

Each FP allele is left-aligned and trimmed on the set's FASTA (`calls_lib.normalize`; MNPs also give their single
bases). It then gets one exclusive class, the first that applies:

| # | class | rule |
|---|---|---|
| 1 | germline, same allele | the normalized allele is carried by the normal's germline truth |
| 2 | germline, overlapping other allele | a germline record of another allele intersects the call's reference interval (deleted bases; insertion = anchor + next base) |
| 3 | near germline (<= 10 bp) | a germline record intersects the interval widened by 10 bp (the labels' `near_bp`) |
| 4 | HPRC v1.1 allele | the normalized allele has AC > 0 among the HPRC v1.1 haplotypes |
| 5 | gnomAD / CoLoRSdb below PoN floor | gnomAD or CoLoRSdb carries the allele (normalized match), with AF below the PoN floor for that allele |
| 6 | gnomAD / CoLoRSdb >= floor, other representation | they carry it at AF >= the floor, but in a representation the PoN's exact POS / REF / ALT match misses |
| 7 | HP>=7 / STR | GRCh38 repeat context (pipeline `RepeatContext`; SNVs via its `tract()`; see E10) |
| 8 | other | none of the above |

Non-exclusive flags (`tables.md` section 5):
- any germline record within 10 bp;
- HPRC allele;
- gnomAD / CoLoRSdb at any AF;
- dbSNP (not somatic-only) or 1000G (the PoN uses neither for INDELs);
- HP>=7, STR;
- outside the germline truth BED (where a germline allele can be missing from the truth).

Each table also gives the median FORMAT AF of the call (tumor reads).

Checks:
- The classified rows equal rtg's FP count in every call-set x evaluation row (`summary.json`: `fp` = `fp_rows`).
- The normalizer was unit-tested on shifted homopolymer and STR indels.

## Scripts and jobs

| file | role |
|---|---|
| `fp_decomposition.py` | per-FP classes -> `DATA/E4b/<call set>.<eval>.fp.v2.tsv.gz`; `tables.md`, `summary.json` |
| `collect.sh` | builds the HPRC table if missing, then runs `fp_decomposition.py` |
| `../call_level_common/prep_hprc.py` | HPRC v1.1 chr1 allele table (run once, login node, 36 s, 1.0 GB) |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E4b_fp_decomposition/job.sh` | Slurm wrapper (`-p general`, 4 CPUs, 12 GB); logs in `logs/` |

Jobs 405339 (first full run, 162 s, 1.04 GB max RSS), 405434, 405458 (snapshot 13:55) and 406568 (all call sets,
18:20, 1,639 s, 1.47 GB max RSS). Reruns classify only the call sets whose rtg FP file is newer than the cache. Job
405325 was an earlier version that I cancelled.

## Results (all call sets, 2026-10-15 18:00; 102 call sets; chr1, in the BED, PASS after the PoN; `tables.md` of 18:20)

### 1. In-sample SNV (chr1, `tables.md` sections 1 and 5)

| chr1 call set | TP | FP | germline same | germline overlap | near germline | HPRC | pop < floor | HP>=7 / STR | other | any germline <= 10 bp (flag) | outside germline BED (other) | median AF of "other" |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HG008T Illumina SNV linear | 548 | 1,089 | 642 (59.0 %) | 90 (8.3 %) | 79 (7.3 %) | 6 | 6 | 33 | 224 (20.6 %) | 74.7 % | 196 (128) | 0.230 |
| HG008T Illumina SNV d9 | 579 | 702 | 609 (86.8 %) | 46 (6.6 %) | 38 (5.4 %) | 0 | 1 | 4 | 2 | 98.7 % | 0 (0) | 0.404 |
| HG008T Illumina SNV d46 | 594 | 792 | 655 (82.7 %) | 55 (6.9 %) | 37 (4.7 %) | 2 | 4 | 17 | 18 | 94.7 % | 8 (4) | 0.211 |
| HG008T PacBio SNV linear | 576 | 969 | 816 (84.2 %) | 22 (2.3 %) | 27 (2.8 %) | 16 | 0 | 5 | 83 (8.6 %) | 89.3 % | 219 (69) | 0.680 |
| HG008T PacBio SNV d9 (*) | 577 | 583 | 571 (97.9 %) | 7 (1.2 %) | 4 (0.7 %) | 0 | 0 | 0 | 1 | 99.8 % | 1 (0) | 0.247 |
| HG008T PacBio SNV d46 | 582 | 650 | 619 (95.2 %) | 11 (1.7 %) | 8 (1.2 %) | 0 | 1 | 4 | 4 | 98.3 % | 5 (3) | 0.514 |
| HG008T ONT SNV linear | 587 | 904 | 707 (78.2 %) | 30 (3.3 %) | 59 (6.5 %) | 12 | 1 | 7 | 83 (9.2 %) | 88.4 % | 73 (58) | 0.500 |
| HG008T ONT SNV d9 (*) | 576 | 660 | 618 (93.6 %) | 13 (2.0 %) | 10 (1.5 %) | 1 | 1 | 11 | 6 | 97.1 % | 0 (0) | 0.314 |
| HG008T ONT SNV d46 | 599 | 788 | 674 (85.5 %) | 31 (3.9 %) | 21 (2.7 %) | 10 | 4 | 19 | 24 | 92.3 % | 8 (3) | 0.332 |
| COLO829T Illumina SNV linear | 1,423 | 2,221 | 101 (4.5 %) | 6 (0.3 %) | 7 (0.3 %) | 7 | 5 | 95 | 1,994 (89.8 %) | 5.2 % | 2,031 (1,928) | 0.286 |
| COLO829T Illumina SNV d46 | 1,745 | 268 | 191 (71.3 %) | 4 (1.5 %) | 8 (3.0 %) | 0 | 2 | 10 | 52 (19.4 %) | 75.7 % | 18 (16) | 0.888 |
| COLO829T fiberseq SNV linear | 1,525 | 521 | 133 (25.5 %) | 1 (0.2 %) | 17 (3.3 %) | 13 | 4 | 23 | 326 (62.6 %) | 29.0 % | 243 (202) | 0.195 |
| COLO829T fiberseq SNV d46 | 1,728 | 268 | 172 (64.2 %) | 2 (0.7 %) | 6 (2.2 %) | 2 | 2 | 12 | 71 (26.5 %) | 67.2 % | 44 (37) | 0.200 |
| COLO829T fiberseq SNV d9, joint-purity model, P100 (not a matched arm) | 1,707 | 231 | 158 (68.4 %) | 5 (2.2 %) | 5 (2.2 %) | 0 | 3 | 10 | 50 (21.6 %) | 72.7 % | 18 (16) | 0.484 |
| COLO829T ONT SNV linear | 1,325 | 2,575 | 126 (4.9 %) | 10 (0.4 %) | 44 (1.7 %) | 64 | 132 (5.1 %) | 88 | 2,100 (81.6 %) | 7.2 % | 1,567 (1,454) | 0.153 |
| COLO829T ONT SNV d46 | 1,690 | 457 | 195 (42.7 %) | 5 (1.1 %) | 14 (3.1 %) | 0 | 25 (5.5 %) | 19 | 197 (43.1 %) | 47.3 % | 53 (36) | 0.207 |

(*) d9 checkpoint of a different recipe (1 GPU; E3 caveat 1). Only HG008T Illumina has a recipe-matched d9 / d46 pair.

- Class 6 (">= floor, other representation") holds 0-11 FPs per SNV call set (most: COLO829T ONT linear, 11).
- The 9 `__norescue` re-calls have exactly the placed FPs of their main runs; the rescue changes only unplaced calls.
- The gnomAD / CoLoRSdb flag at any AF covers 2.2-37.2 % of these FPs. Those the patient carries are counted in class 1
  (precedence); class 5 holds only the rest.
- Germline same-allele FPs that are also HPRC v1.1 alleles (common germline that a pangenome carries):
  - linear: 40 / 48 / 41 (HG008T Illumina / PacBio / ONT) and 1 / 0 / 5 (COLO829T Illumina / fiberseq / ONT);
  - graph arms: 0-2 per call set.

### 2. Linear vs graph arms, SNV (chr1, `tables.md` section 4)

`extra` = linear FP minus graph-arm FP; shares are of the extra FPs.

| chr1 SNV, in the BED | FP linear / d9 / d46 | extra vs d9 / d46 | germline classes of the extra (of which HPRC v1.1) | population classes of the extra | "other" of the extra |
|---|---|---|---|---|---|
| HG008T Illumina | 1,089 / 702 / 792 | 387 / 297 | 118 (30.5 %; 58) / 64 (21.5 %; 57) | 18 (4.7 %) / 11 (3.7 %) | +222 / +206 |
| HG008T PacBio | 969 / 583 (*) / 650 | 386 / 319 | 283 (73.3 %; 59) / 227 (71.2 %; 60) | 16 (4.1 %) / 12 (3.8 %) | +82 / +79 |
| HG008T ONT | 904 / 660 (*) / 788 | 244 / 116 | 155 (63.5 %; 47) / 70 (60.3 %; 45) | 16 (6.6 %) / −1 | +77 / +59 |
| COLO829T Illumina | 2,221 / - / 268 | - / 1,953 | −89 (HPRC 2) | 15 (0.8 %) | +1,942 |
| COLO829T fiberseq | 521 / - / 268 | - / 253 | −29 (HPRC 2) | 16 (6.3 %) | +255 |
| COLO829T ONT | 2,575 / - / 457 | - / 2,118 | −34 (HPRC 7) | 180 (8.5 %) | +1,903 |

- HG008T long reads: the linear arm's extra FPs are mostly germline, chiefly the same allele as the normal (PacBio
  +245 / +197, ONT +89 / +33 vs d9 / d46). Only 47 / 47 and 39 / 40 of those are HPRC v1.1 alleles.
- HG008T Illumina and all three COLO829T platforms: the extra FPs are mostly "other".
- d9 vs d46: d46 has more FPs on all three HG008T platforms: Illumina +90 (germline same allele +46), PacBio +67
  (+48), ONT +128 (+56; other +18, HPRC +9).

### 3. In-sample INDEL (chr1, `tables.md` sections 1, 4 and 5)

The threshold is each run's own: main = validation recall 0.9, valbest = the run's validation-best threshold. E4b does
not give the chr1 curve's best F1 (the INDEL comparison rule); it decomposes the PASS FPs at these thresholds.

| chr1 INDEL, main threshold | TP | FP | P | germline classes (1-3) | population classes (4-6) | HP>=7 (flag) |
|---|---|---|---|---|---|---|
| HG008T Illumina linear | 154 | 1,577 | 0.089 | 83.7 % | 9.9 % | 37.3 % |
| HG008T Illumina d9 | 140 | 434 | 0.244 | 94.0 % | 3.9 % | 69.8 % |
| HG008T Illumina d46 | 141 | 428 | 0.248 | 93.9 % | 3.3 % | 67.3 % |
| HG008T PacBio linear | 153 | 1,705 | 0.082 | 89.5 % | 9.0 % | 42.2 % |
| HG008T PacBio d9 | 150 | 2,636 | 0.054 | 70.7 % | 19.9 % | 71.5 % |
| HG008T PacBio d46 | 151 | 2,452 | 0.058 | 79.6 % | 15.3 % | 66.8 % |
| HG008T ONT linear | 140 | 11,892 | 0.012 | 51.1 % | 42.6 % | 77.7 % |
| HG008T ONT d9 | 134 | 17,253 | 0.008 | 42.6 % | 51.1 % | 91.9 % |
| HG008T ONT d46 | 139 | 18,416 | 0.007 | 35.8 % | 57.3 % | 91.9 % |
| COLO829T Illumina linear | 29 | 1,642 | 0.017 | 59.4 % | 6.3 % | 21.0 % |
| COLO829T Illumina d46 | 32 | 422 | 0.070 | 83.4 % | 6.2 % | 55.5 % |
| COLO829T fiberseq linear | 39 | 1,779 | 0.021 | 74.4 % | 10.5 % | 27.2 % |
| COLO829T fiberseq d46 | 34 | 3,172 | 0.011 | 57.1 % | 38.2 % | 80.8 % |
| COLO829T fiberseq d9, joint model P100 (not matched) | 36 | 9,879 | 0.004 | 25.9 % | 55.3 % | 79.6 % |
| COLO829T ONT linear | 37 | 19,117 | 0.002 | 34.6 % | 46.7 % | 69.8 % |
| COLO829T ONT d46 | 30 | 28,301 | 0.001 | 33.2 % | 43.0 % | 74.0 % |

- **Illumina**: at about the same TP the graph arms have 3.6-3.9 times fewer FPs (HG008T 434 / 428 vs 1,577;
  COLO829T 422 vs 1,642).
  - HG008T: of the linear arm's 1,143 / 1,149 extra FPs (vs d9 / d46), germline classes are 912 (79.8 %) / 918
    (79.9 %), with 207 / 236 HPRC v1.1 alleles; population classes 139 (12.2 %) / 142 (12.4 %).
  - COLO829T: 1,220 extra; germline classes 623 (51.1 %; HPRC 190), population 77 (6.3 %), other +449.
- **Long reads**: at about the same TP the graph arms have more FPs than linear: HG008T PacBio 2,636 / 2,452 vs
  1,705, ONT 17,253 / 18,416 vs 11,892; COLO829T fiberseq 3,172 vs 1,779, ONT 28,301 vs 19,117.
  - The graph surplus sits in homopolymers (HP>=7 flag 66.8-91.9 % of graph FPs vs 42.2-77.7 % linear), in
    germline overlapping-other-allele, and in HPRC / population alleles. Example, HG008T ONT d46 − linear: HPRC
    +1,855, below floor +2,180, other representation +1,449.
  - Class 6 ("population allele >= floor, other representation") is large for ONT: 2,773 (15.1 %, HG008T d46), 2,435
    (14.1 %, d9), 1,324 (11.1 %, linear), 3,226 (11.4 %, COLO829T d46). These are population INDELs that the PoN's
    exact POS / REF / ALT match misses.
- **valbest** cuts the FPs to 0-247 per call set. On HG008T Illumina / PacBio graph arms (121-145 FPs) the germline
  classes are 97.5-98.4 %; on HG008T ONT d9 / d46 (114 / 247) they are 58.8 % / 42.1 %, population 36.8 % / 51.4 %.

### 4. Joint-purity model (COLO829T fiberseq d9; chr1, `tables.md` section 2)

| chr1, joint model | SNV FP | SNV germline same allele | median AF of germline-same FPs | INDEL FP | INDEL HP>=7 (flag) |
|---|---|---|---|---|---|
| P100 | 231 | 158 (68.4 %) | 0.985 | 9,879 | 79.6 % |
| P75 | 545 | 469 (86.1 %) | 0.137 | 14,408 | 83.5 % |
| P50 | 606 | 524 (86.5 %) | 0.235 | 16,690 | 83.3 % |
| P25 | 440 | 351 (79.8 %) | 0.317 | 19,537 | 83.7 % |

- In the mixtures the extra SNV FPs are germline same-allele calls at low AF. The median AFs (0.137 / 0.235 / 0.317)
  are close to (1 − purity) / 2 = 0.125 / 0.25 / 0.375. This fits normal alleles lost in the tumor (LOH) that the
  admixed normal brings back. It is an inference from the AFs, not verified per call.

### 5. Cross-sample call sets (chr1, `tables.md` sections 3 and 5)

- **SNV** (30): the target's germline (same allele) is 52.8-99.2 % of the FPs.
  - HCC1954T targets: 1,123-3,229 FPs per call set (P 0.222-0.373), 93.4-99.2 % germline same allele. Their median
    AF is 0.205-0.274 with the COLO829T models and 0.638-0.729 with the HG008T models.
  - Without HCC1954T: 52.8-98.4 %. The lowest are ONT models on other tumors (HG008T ONT -> COLO829T ONT 52.8 %,
    -> H2009T 55.1 %), which add "other" calls at median AF 0.16-0.26.
- **INDEL** (10, Illumina d46 only):
  - HG008T model on H1437T / H2009T / HCC1954T: 4,017 / 4,472 / 4,321 FPs. Of these, 96.2 / 97.0 / 92.9 % are in
    HP>=7 and 94.2 / 95.1 / 92.3 % are gnomAD / CoLoRSdb alleles. Class 6 holds 31.3 / 35.0 / 29.3 %.
  - The COLO829T model on the same targets: 517 / 601 / 687 FPs.
  - COLO829T Illumina with the HG008T model: 730 FPs. Germline same allele 199 (27.3 %), overlapping 271 (37.1 %),
    near 132 (18.1 %); by flags 54.2 % HP>=7, 21.1 % HPRC v1.1, 15.2 % dbSNP / 1000G (not used by the INDEL PoN).

### Still missing

- COLO829T d9 single-set runs are not local. COLO829T has only linear vs d46 (plus the joint-purity d9 model, which is
  not a matched arm).
- Comparators (ClairS-TO / DeepSomatic-TO) are not yet re-scored, so there is no comparator FP decomposition.
- Seed repeats are pending GPU. So are recipe-matched d9 reruns for HG008T PacBio / ONT SNV (current d9 = 1 GPU,
  marked (*)).
- Cross-sample INDEL exists only for Illumina; HPRC v2.1 haplotypes are not local (d46 is checked against v1.1).
- Unplaced PASS calls are not rtg FPs and are not here (E9c).

## Interpretation (all call sets)

1. **On the graph arms, the SNV FPs left after the PoN are mostly the patient's own germline variants.**
   - HG008T, six graph call sets: 82.7-97.9 % carry the same allele as the normal's germline truth; 92.3-99.8 % have
     a germline record within 10 bp.
   - At most 2 of these germline FPs per call set are HPRC v1.1 alleles. Population alleles below the PoN floor are
     at most 4.
   - COLO829T d46: 42.7-71.3 % germline same allele; the rest is mostly "other" (point 3).
   - So the remaining SNV precision limit is germline variation that neither the graph nor the PoN carries. In
     tumor-only mode nothing upstream can see it. This agrees with E9c / E9d.
2. **The linear arm's extra SNV FPs depend on the platform and the tumor (R3.2, R2.5).**
   - HG008T PacBio / ONT: mostly germline (71-73 % / 60-64 % of the extra), chiefly same-allele germline calls, and
     45-60 of them are HPRC v1.1 alleles. Here the pangenome removes germline FPs.
   - HG008T Illumina: germline 21.5-30.5 % of the extra; the largest part is "other" (+206 / +222, median AF 0.230,
     128 of 224 outside the HG008-N germline BED).
   - COLO829T, all three platforms: almost all of the extra is "other": +1,942 (Illumina), +255 (fiberseq), +1,903
     (ONT). Of the linear arms' "other" FPs, 1,928 of 1,994, 202 of 326 and 1,454 of 2,100 lie outside the COLO829BL
     germline BED, at median AF 0.286 / 0.195 / 0.153.
   - E10 places them in SMaHT "extreme" regions. The linear arms make 351-2,210 PASS calls there at precision
     0.056-0.253: HG008T Illumina / PacBio / ONT 562 / 458 / 351 at 0.059 / 0.090 / 0.105; COLO829T 2,210 / 494 /
     2,055 at 0.056 / 0.253 / 0.056. The graph arms make 102-281 at 0.161-0.764.
   - Low AF, extreme regions and outside the assembly-confident BED: this fits reads from paralogous or alternative
     haplotype sequence misplaced on GRCh38-only (reference bias), not germline alleles the PoN missed. It is an
     inference from the strata, not verified per call.
3. **COLO829T graph arms have a larger "other" class.**
   - Illumina: 52 FPs at median AF 0.888, candidates for truth incompleteness (somatic events missing from the
     COLO829T truth); not verified here.
   - fiberseq 71 at median AF 0.200 and ONT 197 at 0.207: low-AF calls that look like sequencing or alignment errors.
   - Only 16 / 52, 37 / 71 and 36 / 197 of them lie outside the COLO829BL germline BED.
4. **INDEL: the graph removes germline FPs for Illumina; long reads are limited by homopolymer errors.**
   - Illumina: the graph arms keep 3.6-3.9 times fewer FPs at about the same TP. About 80 % (HG008T) and 51 %
     (COLO829T) of the linear arm's extra FPs are germline INDELs.
   - PacBio, ONT, fiberseq: at the recall-0.9 threshold most FPs are in HP>=7 (graph 66.8-91.9 %) or are population
     INDELs. The graph arms make more of them than linear at about the same TP. These are the long-read
     homopolymer-length errors (E17a). The FP counts are at different thresholds; the INDEL ranking uses the chr1
     curve's best F1 (benchmark), not these counts.
   - The PoN's exact-representation match misses 11-15 % of the ONT INDEL FPs and 29-35 % of the HG008T Illumina
     model's FPs on H1437T / H2009T / HCC1954T (class 6). Normalizing before the PoN match is a candidate fix; it is not
     applied or tested here.
5. **Cross-sample and joint models.**
   - Cross-sample FPs are again mostly the target's germline. HCC1954T stands out: 93-99 % germline same allele,
     1,123-3,229 FPs per call set.
   - The joint model's mixtures add low-AF germline FPs at about the normal's share of the AF (section 4).

## Caveats

- chr1 only. Several classes hold only a handful of FPs per call set.
- The germline truth includes non-PASS dipcall records. On every in-sample graph SNV call set each same-allele match
  is a PASS (`.`) record (section 5). On the linear SNV arms some are not: HG008T Illumina 638 of 642, PacBio 692 of
  816, ONT 700 of 707; COLO829T Illumina 98 of 101, fiberseq 124 of 133, ONT 117 of 126. INDEL call sets differ by
  0-29.
- The germline truth covers only its confident BED (dipcall; DeepVariant for the CASTLE normals). The COLO829T linear
  arms' "other" FPs lie mostly outside it, so some of them may be germline that the truth does not cover.
- INDEL FPs are at each run's own threshold. The arms reach about the same chr1 TP there, but this is not the INDEL
  comparison rule (the chr1 curve's best F1, `analysis/benchmark_7fc7396_20261012`).
- HPRC v1.1 haplotypes only (v2.1 not local). The HPRC class is checked after the germline classes, so it counts only
  population alleles that the patient does not carry.
- The classes are exclusive by precedence. The repeat class holds only FPs that matched no germline or population
  allele; the HP>=7 and STR flags give the full counts.
- Unplaced (graph-coordinate) PASS calls are not rtg FPs and are not in this table. E9c classifies them.

## Rerun

All four call-level analyses: `bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/call_level_common/collect_all.sh`.
This one only:
`sbatch -J E4b_fp_decomposition -o /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E4b_fp_decomposition/logs/%x-%j.out /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E4b_fp_decomposition/job.sh`.
The run rewrites `tables.md` and `summary.json`; this README's numbers are those of 18:20 (job 406568).
