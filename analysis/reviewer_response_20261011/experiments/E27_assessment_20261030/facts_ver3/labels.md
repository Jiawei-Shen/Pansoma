# Fact sheet v3: labels, denominators, graph representation (R1.5, R1.8, R3.2, Note R2.3) — 2026-10-30

Tags: **EXISTS** (already in a repo/data file, path given), **COMPUTED NOW** (today, script + output in
`V = /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/v3/labels/`, git-ignored; logs in `V/logs/`),
**MISSING** (not stored, not computable here). No tracked repo file was edited; nothing staged.
Paths: `EXP = /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments`,
`RUNS = /scratch/jshen/data/pansoma_net_v2_runs`, `TP = /scratch/jshen/Github/Pansoma/indexed_gam_pipeline_v4/tensor_postprocessing`,
`NET = /scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2`.

Scripts run today:

| script | where | output |
|---|---|---|
| `V/hapmap_mask_breakdown.py` | login node | `V/hapmap_mask_breakdown.json`, `V/logs/hapmap_mask_breakdown.log` |
| `V/rtg_region_rule_check.py`, `V/rtg_region_rule_check2.py` | login node | `V/rtg_region_rule_check.json`, `V/logs/rtg_region_rule_check{,2}.log` |
| `V/offref_val_effect.py` | login node | `V/offref_val_effect.json`, `V/logs/offref_val_effect.log` |
| `V/offref_label_counts.py` (`.sbatch`) | Slurm 462920, miao | `V/offref_label_counts.json`, `V/logs/offref_label_counts-462920.out` |
| `V/germline_somatic_conflict.sbatch` | Slurm 462921, meow | `V/logs/germline_somatic_conflict-462921.out` |
| `V/r32_chm13_exception.py` | login node | `V/r32_chm13_exception.json`, `V/logs/r32_chm13_exception.log` |

---

## Q1. R1.5 HapMap denominators

### Which script made the numbers (EXISTS)

- Mask: `EXP/E23_hapmap_lowvaf/exclude_graph_carried.py` (job 460499, 2026-10-29; submit script
  `/scratch/jshen/Github/Pansoma/tmp/hapmap_nonHG005_20261012/train_nograph/submit.sh`). It writes
  `/scratch/jshen/data/reviewer_experiments_20261013/E23_nograph/eval_chr1.bed` = chr1 benchmark BED minus, for every
  graph-carried truth allele (truth ∩ atomized graph alleles, any VAF, SNV and INDEL), the interval
  `[pos-2, pos-1+len(REF)+1)` (0-based), i.e. the allele's REF span **±1 bp**. Log: 216,860,000 bp → 216,195,801 bp after
  removing 207,432 merged spans (`.../train_nograph/exclude-460499.out`).
- "Graph-carried" = the truth allele (CHROM POS REF ALT) is an allele of `mgi_non_HG005.d9`'s own vg deconstruct VCF,
  all snarl levels, with alleles ≤ 50 bp atomized (`EXP/E23_hapmap_lowvaf/graph_alleles.sh`, ATOMIZE=1;
  `/scratch/jshen/data/reviewer_experiments_20261013/E23_graph_alleles.mgi_non_HG005.d9.atomized.tsv.gz`). This is a real
  graph-content check.
- Scoring: `/scratch/jshen/Github/Pansoma/tmp/mutect2_hapmap_20261029/evaluate.sh` and `pansoma_nograph.sh` (job 461100):
  rtg vcfeval `--squash-ploidy --sample ALT,ALT -e eval_chr1.bed` (Pansoma also `--bed-regions` = all chr1), then
  `EXP/E23_hapmap_lowvaf/call_vaf_by_value.py` (truth records with FORMAT AF 0.0099–0.1651).
  Outputs: Pansoma `RUNS/HapMap_Illumina_{SNV_full,INDEL}/vcf_chr1/eval_bed_nograph_all/`, Mutect2
  `/scratch/jshen/data/reviewer_experiments_20261013/E23_mutect2/eval.pass.{SNV,INDEL}/`; tables
  `EXP/E23_hapmap_lowvaf/mutect2/call_vaf_by_value.*.md`.
- Full-truth scoring: Pansoma = the run's own `vcf_chr1/eval_bed/` (vcf.sh; full v1.6 BED), table
  `EXP/E23_hapmap_lowvaf/atomized/call_vaf_by_value.HapMap_Illumina_{SNV_full,INDEL}.md`; Mutect2 = Slurm 462880,
  `/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators/hapmap_m2_full/` (same truth, same full BED).

### The extra exclusions are masking side effects (COMPUTED NOW, `V/hapmap_mask_breakdown.json`)

| | SNV | INDEL | tag |
|---|---|---|---|
| truth (chr1, BED, VAF 1–16.5 %) | 392,567 | 124,499 | EXISTS (atomized tables) |
| graph-carried | 175,773 (44.8 %) | 25,938 (20.8 %) | EXISTS |
| not graph-carried (full-BED split) | 216,794 | 98,561 | EXISTS |
| not graph-carried but dropped by the mask | **1,701** | **6,482** | COMPUTED NOW |
| – same VCF position as a graph-carried allele (another allele of the same site; the truth VCF is split-allelic) | 383 | 5,009 | COMPUTED NOW |
| – REF span overlaps a graph-carried allele (other position) | 131 | 377 | COMPUTED NOW |
| – within the ±1 bp margin of a graph-carried allele | 1,187 | 823 | COMPUTED NOW |
| – insertion whose insertion point abuts the mask (graph-carried allele 2 bp away) | 0 | 273 | COMPUTED NOW |
| graph-carried truth NOT masked (REF base `N`: chr1:248,755,122 N>G, absent from the label truth keys) | 1 | 0 | COMPUTED NOW |
| **scored in eval BED** = 216,794 − 1,701 + 1 / 98,561 − 6,482 | **215,094** | **92,079** | EXISTS, closes |

- Of the 5,009 same-position INDELs, 4,686 sit at a site whose graph-carried allele is an INDEL (other repeat lengths at
  the same STR/homopolymer, e.g. chr1:744,866 G>GTA×n); the masking allele is at VAF 1–16.5 % for 1,137+363+123 SNVs /
  4,781+353+798 INDELs, else at VAF < 1 % (also masked).
- rtg `-e` rule that explains every record (COMPUTED NOW, `V/rtg_region_rule_check*.log`): a truth record is kept iff its
  anchor-trimmed allele lies in the eval BED — SNV: the base (392,567 / 392,567 agree); DEL: any deleted base
  (56,617 / 56,617); INS: both bases flanking the insertion point (67,882 / 67,882).
- Not normalisation/decomposition: the truth and both evaluations are the same `truth.{SNV,INDEL}.vcf.gz` (byte-identical
  records in `eval_bed` and `eval_bed_nograph_all`); no record is outside F (E − F = 0).

### Same denominator for both tools; full-truth TPs recomputed (COMPUTED NOW / EXISTS)

- Mutect2 and Pansoma eval-BED truth record sets are **identical** (215,094 / 92,079; `mutect2_eval_equals_pansoma_eval: true`).
- Full-truth TPs are NOT conditional TPs carried over; each comes from a separate vcfeval on the unmasked BED with
  graph-carried truth as FN. Closure for Pansoma: SNV 117,834 (eval BED) + 474 (TPs among the 1,701) + 95 (graph-carried
  TPs) = 118,403 → 0.302; INDEL 19,219 + 899 + 525 = 20,643 → 0.166. Mutect2 full: 247,194 / 392,567 = 0.630;
  18,682 / 124,499 = 0.150 (EXISTS, `.../hapmap_m2_full/call_vaf_by_value.Mutect2_pass.*.md`).
- Conditional recall on the unmasked split (216,794 / 98,561), for reference: Pansoma 0.546 (118,308) / 0.204 (20,118);
  Mutect2 0.578 (125,363) / 0.112 (11,058) (EXISTS, atomized + hapmap_m2_full tables). The eval-BED values: Pansoma 0.548
  / 0.209, Mutect2 0.581 / 0.117. The comparison does not change.
- The precision row is from the eval-BED (conditional) evaluation: Pansoma 0.990 (118,388 / 119,549), 0.784
  (17,966 / 22,920); Mutect2 0.975 (125,215 / 128,401), 0.860 (9,303 / 10,818). On the full BED: Pansoma 0.990 / 0.776,
  Mutect2 0.986 / 0.900 (EXISTS).

---

## Q2. R1.8 label precedence (EXISTS: `TP/truth_labels.py` `classify()` lines 589–690, `truth_overlap()` 497–517; `TP/README.md` §Labels)

Order in code (first match decides; label values 1 somatic, 2 germline, 0 non/"other", −1 ignore):

1. AF floors (unset for HG008T/COLO829T) → −1.
2. A1 (representative allele, most ALT reads) **equals a somatic truth key**: PASS → 1; only filtered → −1. Inside or outside BEDs.
3. A1 **equals a germline truth key** (PASS or GAP1/GAP2) → 2. Inside or outside BEDs. (HET1/HET2 noted for step 8.)
4. Another allele of the site equals a truth key: A1 overlap > 45 % with a **somatic** one → 1 (`allele_partial_somatic_truth`);
   else > 45 % with a **germline** one → 2; else → 0. (No BED test.)
5. Another allele is a filtered somatic truth → −1.
6. No GRCh38 position and no anchor → −1.
7. Outside confident region (somatic BED ∩ germline BED) → −1 (`outside_confident_region`).
8. A1 = germline truth with HET1/HET2 → 0.
9. Somatic PASS truth "at the same place", allele or read-haplotype overlap > 45 % → 1 (`residual_partial_somatic_truth`).
10. Off-reference (anchored), nothing above → −1 (`off_reference_no_truth_match`).
11. Germline truth at the same place, allele overlap > 45 % → 2 (`residual_partial_germline_truth`; not tried off-reference).
12–13. → 0.

**Conflict answer:** exact germline (step 3) beats partial somatic (steps 4, 9): such a candidate is **germline**, a
training negative. The letter's sentence "Somatic matches are tried before germline matches, so a candidate that partly
matches a somatic allele is excluded rather than labelled germline" is true only partial-vs-partial (step 4 somatic
before germline; step 9 before 11) and for exact somatic vs anything. Measured conflicts (COMPUTED NOW, Slurm 462921,
chr1–22, d46 sets): A1 exact germline while another allele of its site is a somatic truth allele: SNV 0 / 0 / 1 (HG008T
Illumina / PacBio / ONT), 0 / 0 / 0 (COLO829T); INDEL 49 / 270 / 247 and 9 / 31 / 40 (of 53,483–175,203 exact-germline
INDEL tensors). Residual (step-9) conflicts are not recorded by the labeller (it returns at step 3): MISSING, would need
an instrumented relabel.

Other order facts: exact-somatic tensors with a germline site mate (somatic wins): SNV 0–9, INDEL 8–206 per set
(COMPUTED NOW, same job). Outside-confident-region tensors are −1 before steps 9/11, i.e. their residual partial matches
are never tested.

**Use of the labels** (EXISTS: `RUNS/ris_7fc7396_20261012/*_SNV_d46_normbq*/args.json`, `RUNS/release_s*_*_INDEL/args.json`;
`NET/train.py` 545–572, 419–468; `NET/data.py` 14–20, 56, 225–228; `NET/metrics.py` 70–91, 247–296; `NET/predict.py`
docstring; `NET/vcfeval.py` docstring; `NET/linear_vcf.py` 13–19):

| | training | validation (checkpoint + PASS threshold) | chr1 benchmark |
|---|---|---|---|
| exact somatic (2) | positive | positive | — (benchmark = rtg on calls) |
| germline (3, 4g, 11) | negative | negative | — |
| other (4, 8, 12, 13) | negative | negative | — |
| partial somatic (4s, 9) | left out (`--ignore-reasons`, all final models) | left out | INDEL tensor-level test counts them as found; rtg unaffected |
| off-path, no match (10) | left out (−1) | **negative** (`EVAL_AS_NON`) | calls → `*.linear.unplaced.vcf.gz`, **not scored by rtg**; R3.4 worst case counts them FP |
| outside confident region (7) | left out; **final INDEL models: 0 if pos0/anchor-middle inside the somatic truth BED** (`--non-reasons-in-bed`) | same as training | rtg scores calls inside the benchmark BED (`-e`) |
| filtered somatic / no position (2f, 5, 6) | left out | left out | — |

- Checkpoint: SNV `--best-by f1` = highest validation per-tensor somatic F1 before population filtering (at its F1-optimal
  threshold); final INDEL `--best-by postpon` = highest per-tensor post-filter F1 (PoN-tagged tensors never called).
- PASS threshold: SNV = highest p_somatic at which validation per-tensor recall reaches 0.9 (`somatic_at_recall`);
  INDEL "ppr8" = post-filter recall of the untagged (callable) positives reaches 0.8 (`postpon_report.at_recall_callable`;
  e.g. release_s0_Illumina 0.38437 = `vcf_ppr8_*/...graph.stats.json`). Both depend **only on positives**.
- Benchmark (confirmed): `NET/vcfeval.py` = rtg vcfeval `--squash-ploidy --sample ALT,ALT --bed-regions <chroms> -e <BED>`,
  PASS truth of the kind, population filter applied to calls; "The calls (graph_vcf) and so rtg vcfeval do not depend on"
  tensor labels (`predict.py`). Unplaced PASS calls are counted in `report.json` `unplaced_pass` only.

---

## Q3. Off-path tensors without truth

### Labels (EXISTS, code above)
Training −1; validation 0 (`data.py` `EVAL_AS_NON = ("off_reference_no_truth_match",)`, comment "-1 for training, 0 (non)
when evaluating"); predict chr1 tensor-level test 0; rtg benchmark: not scored (unplaced VCF). Not "other": not used in training.

### Counts (EXISTS: each model's `train.log` line "validation includes N off-reference tensors without a truth match, as non"; parsed in `V/offref_val_effect.json`)

| models | validation off-path negatives | share of validation non |
|---|---|---|
| SNV normbq (18: 6 sets × 3 seeds) | 1,367–4,347 (HG008T Illumina 3,951–4,347; PacBio 1,367–1,784; ONT 1,391–1,748; COLO829T Illumina 1,899–2,239; Fiber-seq 2,321–2,546; ONT 3,209–3,519) | 4.7–42.5 % (PacBio/Fiber-seq highest: 21.9–42.5 %) |
| final INDEL release_s{0,1,2} (9; 6 tumors joint) | Illumina 8,813–9,314; PacBio 31,560–32,686; ONT 53,317–54,260 | Illumina 27.3–29.4 %, PacBio 7.4–7.9 %, ONT 6.6–6.7 % |

Single-tumor sets, seed-0 validation blocks (COMPUTED NOW, `V/offref_label_counts.json`): SNV 4,059 / 1,784 / 1,748
(HG008T) and 1,899 / 2,546 / 3,519 (COLO829T) — equal to train.log; INDEL 838 / 5,205 / 7,504 and 909 / 7,729 / 12,117.

### Effect on threshold and checkpoint

- **PASS threshold: unaffected, exactly** (EXISTS, code): `metrics.somatic_at_recall` and `postpon_report.at_recall_callable`
  pick the highest p_somatic whose cumulative positives reach the recall target; negatives never enter.
- **Checkpoint: cannot be recomputed from stored data** (MISSING): runs keep only `best.pth`, `last.pth` and per-epoch
  aggregate `metrics.jsonl`; no per-epoch checkpoints or per-tensor validation predictions. It would take re-training the
  models (GPU) with train.py logging the selection score with `off_reference_no_truth_match` left out of validation (same
  seeds, otherwise unchanged), or saving per-epoch validation predictions; priority = the 9 INDEL release models.
- **Proxy (COMPUTED NOW, `V/offref_val_effect.json`)**: stored chr1 test predictions of the selected checkpoints
  (`*/test_chr1*/…predictions.ndjson.gz`), same model, same tensor population:

| | SNV (18 tests) | final INDEL (18 tests: HG008T + COLO829T) |
|---|---|---|
| off-path negatives among FPs at the **selection** threshold | 0–3 of 413–907 (≤ 0.4 %) | 0–64 of 10–92 (0–78 %) |
| per-tensor somatic F1 change if they were left out | +0.0000 to +0.0005 | +0.000 to +0.049 (release_s2_Illumina: COLO829T 38/49 FP, +0.049; HG008T 64/92, +0.024) |
| off-path negatives among FPs at the **PASS** threshold | 1–44 of 899–3,013 | 9–177 of 36–704 (7–80 %) |
| validation margin best − runner-up epoch (`metrics.jsonl`) | 0.0001–0.0199 (< 0.001: 3 models) | 0.0002–0.0319 (< 0.001: release_s1_PacBio, release_s2_ONT) |

Reading: for SNV models the off-path negatives cannot matter beyond swapping near-tied epochs (validation F1
difference < 0.001) [inference]. For the final INDEL models they are a material share of post-filter false positives
(population filtering cannot tag an allele without a GRCh38 representation, `pon_tag −1`), so excluding them could
change the selected epoch of several models [inference; not testable without re-training].

### The "83–717 off-path INDEL negatives" (COMPUTED NOW, `V/offref_label_counts.json`)

- Seed-0 training blocks: 119 / 83 / 332 (HG008T Illumina / PacBio / ONT), 101 / 90 / 717 (COLO829T); SNV 0 (reproduces
  the earlier fact sheet). Label 0 by step 4 (`truth_matches_non_representative_allele`).
- What they are (chr2–22, all blocks): **1,518 / 1,518 are deletions spanning several nodes (`@` multi-node key) whose
  walk cannot be projected to GRCh38**, at a site whose other allele (a GRCh38-node allele) is a truth allele:
  germline 1,509, somatic 8, both 1. Example `COLO829T_ONT_d46` chr20: A1 `110194393:85:DEL:TT>@110194397`, site mate
  `110194393:85:DEL:T>` = germline truth chr20:3,988,457 CTT>CT.
- Why 0: with no GRCh38 position the allele overlap is set to 0 (`truth_overlap`: `... if lin else 0.0`), and for germline
  truth no haplotype overlap is tried (for somatic truth chrom is None, so neither) → step 4 "else 0". So they are
  negatives **by site** (a different allele at a truth site), not by a measured mismatch; at the 1,509 germline sites they
  are not somatic either way; the 9 at somatic sites are potential mislabels.
- **Second off-path training-negative route, not in the letter** (final INDEL models only): off-path tensors with reason
  `outside_confident_region` whose anchor middle lies in the somatic truth BED are trained as 0 by
  `--non-reasons-in-bed` (`train.in_somatic_bed` uses pos0 = anchor middle off-reference): seed-0 training blocks
  618 / 3,010 / 4,050 (HG008T) and 1,360 / 9,270 / 12,924 (COLO829T; the COLO829T somatic BED contains all of them);
  validation 22 / 132 / 112 and 36 / 378 / 480. SNV models: none.

---

## Q4. R3.2 "represented in the graph"

- What the analysis did (EXISTS): carrier-count **proxy**, not a graph-path check. `EXP/E7_graph_composition/s4_calibration.py`
  line 12: "carrier count c among the 463 v2.1 haplotypes; 'predicted in d46' = c >= 46". Carriers from
  `/scratch/jshen/data/hprc_v2.1_vcf/hprc-v2.1-mc-grch38.wave.vcf.gz` (HPRC v2.1 MC-GRCh38, vcfbub -l 0 -a 100000,
  vcfwave, vcffixup; 232 columns = CHM13 haploid + 231 diploid = 463 haplotypes; GRCh38 = backbone, not counted);
  alleles left-aligned and split into primitives, a multi-primitive truth allele counts only if ≥ 46 haplotypes carry all
  primitives (AND of masks) (E7 README "Method", "HPRC v2.1 axis"). E7 README: "There is no d46 counterpart of the d9
  deconstruct VCF ... so for v2.1 the graph-content check ... is not possible".
- The letter's strata reproduce exactly (COMPUTED NOW, `V/r32_chm13_exception.json`; also EXISTS in
  `EXP/E7_graph_composition/v2.1/d46_content.tsv`): somatic SNV 143 / 63; somatic INDEL 1,522 (19.2 %) / 340 (17.8 %);
  germline SNV 87.3 % / 87.7 %; germline INDEL 77.8 % / 77.2 % (HG008 / COLO829, PASS, in BED, chr1–22).
- Reference-path exception (COMPUTED NOW): the proxy counts CHM13 as one of the 463 filtered haplotypes, but in the real
  d46 graph the CHM13 (and GRCh38) paths are never filtered (R3.m3; E7 README note 2026-10-15). Alleles carried by CHM13
  with < 46 carriers (in the "no" rows, yet spelled by the CHM13 path): germline SNV 29,896 (0.8 %) / 26,494 (0.7 %);
  germline INDEL 9,388 (1.0 %) / 8,887 (0.9 %); somatic SNV 4 / 4; somatic INDEL 106 (1.3 %) / 23 (1.2 %).
  GRCh38 alleles are the reference and are never variants.
- Calibration that supports the proxy (EXISTS, E7 README §5 / `v2.1/tables.md` G): germline SNV candidate loss on d46 is
  17.0 % at 35–45 carriers and 98.4 % at 46–57 (HG008T PacBio; COLO829T Fiber-seq 14.6 → 98.1 %); 97.8–98.0 % of the
  germline SNV loss falls on predicted alleles. Somatic INDELs with no v2.1 carrier still lose 29.3 % / 26.2 %.
- The v1.1 sentence of R3.2 (87 % / 75 %) uses real d9 content (the d9 graph's deconstruct VCF, `EXP/E7_graph_composition/d9_content.tsv`), a graph-content check; R2.3's "absorbed" uses a real path check
  (an HPRC v2.1 path spells the ALT, E5b). Three different definitions.

---

## Q5. Note R2.3 "61 of the 702" vs 697 (EXISTS + COMPUTED NOW)

- Same truth version: GIAB HG008-T somatic smvar benchmark v0.2 (`/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz`).
- 702 = all PASS chr1 truth SNV records, inside and outside the benchmark BED; 61 tagged by the population filter
  (`RUNS/ris_7fc7396_20261012/HG008_Illumina_SNV_d46_normbq/vcf_chr1/eval_bed/report.json` `truth_pon_tags`:
  records 702, tagged 61; gnomAD 39, dbSNP 31, 1000G 1, CoLoRSdb 42), EXISTS.
- 697 = the same records inside the benchmark BED (rtg baseline tp-baseline + fn); 60 of them tagged, 637 untagged
  (COMPUTED NOW from `truth.SNV.pon.vcf.gz` × `vcfeval/raw_calls`; matches the earlier sheet's 0.914 = 637/697).
  So 5 PASS chr1 truth SNVs lie outside the BED, 1 of them tagged.
- Consistent numbers: **60 of 697 (8.6 %)** in the benchmark regions (the denominator every other chr1 SNV number uses).

---

# Draft text

## R1.5 table rows (replace the row "evaluated after masking graph-carried sites")

| | SNVs | INDELs |
|---|---|---|
| truth alleles | 392,567 | 124,499 |
| carried by the graph (left out of the conditional evaluation) | 175,773 (44.8%) | 25,938 (20.8%) |
| not carried by the graph | 216,794 | 98,561 |
| scored in the conditional evaluationᵃ | 215,094 | 92,079 |
| Pansoma recall: conditional / full truth | 0.548 / 0.302 | 0.209 / 0.166 |
| Mutect2 recall: conditional / full truth | 0.581 / 0.630 | 0.117 / 0.150 |
| precision of calls at VAF 1–16.5% (conditional evaluation): Pansoma / Mutect2 | 0.990 / 0.975 | 0.784 / 0.860 |

ᵃ The conditional evaluation removes each graph-carried allele's reference span ±1 bp from the benchmark regions, so their calls are not scored either; this also removes 1,701 SNVs and 6,482 INDELs that the graph does not carry but that share a position with a graph-carried allele (383 and 5,009, mostly other repeat lengths at the same INDEL site) or overlap or adjoin one (1,318 and 1,473), and it misses one graph-carried SNV whose reference base is N. Both callers are scored on the same 215,094 / 92,079 alleles, and full-truth recall comes from a separate evaluation on the unmasked regions (graph-carried alleles counted as misses).

Optional sentence: "Without these side effects (216,794 / 98,561 alleles), conditional recall is 0.546 / 0.204 for Pansoma and 0.578 / 0.112 for Mutect2."

## R1.8 rules (replace the rule table and the sentence "Somatic matches are tried before germline matches ...")

*Labels* are decided from the candidate's representative allele (the allele with most ALT reads); the first rule that applies decides:

1. It equals a PASS somatic truth allele (inside or outside the benchmark regions): **somatic**. It equals only a filtered somatic truth allele: excluded.
2. It equals a matched-normal germline truth allele: **germline**, also when another allele of its site is a somatic truth allele (9–270 INDEL candidates per dataset, at most one SNV).
3. Another allele of the same candidate site is a truth allele: **partial somatic** if the representative allele overlaps a somatic one by more than 45%; otherwise **germline** if it overlaps a germline one by more than 45%; otherwise **other**.
4. No GRCh38 position and no anchoring reference node: excluded.
5. Outside the confident region (somatic ∩ germline benchmark regions): excluded.
6. A somatic truth allele at the same place overlaps it by more than 45% (allele, or read haplotype in the truth's GRCh38 window): **partial somatic**.
7. Off the GRCh38 path with none of the above: **no label** (unknown).
8. A germline truth allele at the same place overlaps it by more than 45%: **germline**.
9. Everything else: **other** (sequencing and alignment errors, including errors next to true variants).

*Training:* somatic = positives; germline and other = negatives; partial somatic, no label and excluded candidates are not used. The final INDEL models also train candidates of rule 5 that lie inside the somatic truth regions as other.
*Validation* (checkpoint and PASS threshold): the same labels, except that off-path candidates without a label (rule 7) count as negatives and candidates of rule 5 are left out (final INDEL models: as in training). The checkpoint has the highest validation somatic F1 per candidate (SNVs before, INDELs after population filtering). The PASS threshold is the score at which validation recall reaches 0.9 (SNVs) or 0.8 of the positives that population filtering keeps (INDELs); it depends on the positives only.
*Benchmark* (chr1): no candidate labels. rtg vcfeval compares the PASS calls on GRCh38, after population filtering, with the PASS truth alleles inside the benchmark regions (alleles only, --squash-ploidy), so representation differences are resolved by rtg's haplotype matching. Off-path calls without a GRCh38 position are not scored; R3.4 reports them and, as a worst case, counts each as a false positive.

## Off-path negatives (replace "off-path INDEL candidates that match another allele of a truth site are negatives (83–717 per training set)")

Off-path candidates without any truth match (rule 7) carry no evidence either way and are not training examples. In validation we score them as negatives, a conservative assumption: they are 4.7–42.5% of the validation negatives of the SNV models and 6.6–29.4% of those of the final INDEL models. The assumption cannot move the PASS threshold, which depends only on positives. It can affect checkpoint selection: in the chr1 test predictions of the selected checkpoints, off-path candidates make at most 3 of 413–907 false positives at the SNV selection threshold (F1 change ≤ 0.0005), but 0–64 of 10–92 at the INDEL threshold (F1 change ≤ 0.05), because population filtering cannot remove an allele that has no GRCh38 representation. [待补/limitation: per-epoch validation predictions were not kept, so whether leaving them out would select other INDEL epochs needs re-training.] Two kinds of off-path candidates are training negatives. 83–717 per training set are deletions spanning several nodes whose path cannot be projected to GRCh38, at a site whose GRCh38 allele is a truth allele (1,509 of 1,518 a germline allele); they are labelled other because a different allele of the site is the truth, without a measured overlap. The final INDEL models also train off-path candidates outside the confident region but inside the somatic truth regions as other (618–12,924 per training set), as they do on-path ones.

## R3.2 column definition (replace "represented in the graph (carried by at least 46 of the 463 HPRC v2.1 haplotypes, the criterion of the d46 filter)")

Column header: "≥ 46 of 463 HPRC v2.1 haplotypes carry the allele (proxy for d46 graph content)".

Definition: an allele is counted as represented when at least 46 of the 463 non-GRCh38 haplotypes of the HPRC v2.1 release VCF (CHM13 and 231 diploid individuals; alleles left-aligned and split into SNVs and indels, all parts carried by the same haplotypes) carry it, the threshold of the d46 filter. This is a carrier-frequency proxy, not a test that a path of the distributed d46 graph spells the allele (no variant file of the filtered graph is published); it tracks the graph closely for germline SNVs (germline SNVs lose 15–17% of their candidates on HPRC v2.1 at 35–45 carriers and 98% at 46–57). It treats CHM13 as an ordinary haplotype, whereas the graph keeps the CHM13 and GRCh38 reference paths unfiltered; alleles carried by CHM13 but by fewer than 46 haplotypes (about 1% of germline variants, 4 somatic SNVs and 23–106 somatic INDELs per tumor) are therefore in the "no" rows although the graph contains them. The HPRC v1.1 shares (87% / 75%) are exact alleles of the d9 graph itself.

## Note R2.3 fix (replace "The population filter tags 61 of the 702 HG008T truth SNVs (8.7%)")

"The population filter tags 60 of the 697 HG008T chr1 truth SNVs in the benchmark regions (8.6%; 61 of 702 including the five PASS truth SNVs outside them); ..."
