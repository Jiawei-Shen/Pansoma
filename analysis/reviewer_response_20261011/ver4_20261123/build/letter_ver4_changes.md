# Reviewer letter ver4: what was revised (2026-11-23)

`Reviewer_Letter_Response_ver4.docx` = ver3 (`Reviewer_Letter_Response_ver3.docx`, saved 15:37) with the data of the
tables and figures moved to the latest pipeline, and the captions and sentences that restate them updated. Built by
`analysis/reviewer_response_20261011/ver4_20261123/build/build_ver4.py` (edit steps in `steps.py`; every operation is
logged in `Reviewer_Letter_Response_ver4.buildlog.txt`). Use Word's Compare (ver3 against ver4) to see every change in place.

**Models used ("latest pipeline")**: the six-tumor models (HCC1395T, HCC1954T, H1437T, H2009T, HG008T, COLO829T), seed 0
unless stated, tested on chr1 of HG008T / COLO829T / HCC1395T × Illumina / PacBio HiFi (COLO829T: Fiber-seq) / ONT.
- SNV: population site features as a model input; AF 0.05 candidate sets for Illumina and PacBio HiFi / Fiber-seq (0.08 for
  HG008T and ONT); PASS at validation recall 0.95 after the population filter; chr1 of the three benchmark tumors and chr20
  of HCC1395T held out.
- INDEL: the final release models (three seeds); PASS ppr8 on Illumina / PacBio HiFi, pp on ONT; chr1 of all six tumors and
  chr20 of HCC1395T held out.
- Three-reference comparisons: the six-tumor models trained on GRCh38-only (no site features), HPRC v1.1 and HPRC v2.1.

No model was trained for ver4. New computation was run only where the result did not exist on these models (section 4).

## 1. Tables and figures

| ver4 item | ver3 content | ver4 content | source (`ver4_20261123/…`) |
|---|---|---|---|
| Response Figure 1.2 (masking outside the candidate node) | October single-tumor models, HG008T + COLO829T, 6 rows | six-tumor models, 9 rows (HCC1395T added): recall falls in 9 of 9 (0.204–0.344), F1 in 9 of 9 (0.057–0.243), precision in 5 of 9; panel b scale ±0.35 | masking/ |
| Response Figure 1.3 (path-count channel) | October single-tumor models, 6 rows (COLO829T Illumina ΔF1 +0.311) | the only with / without pairs on the current pipeline: COLO829T single-tumor, Illumina and Fiber-seq, HPRC v2.1 and v1.1 (4 rows); seed range of the with-channel models added; panel b scale ±0.15 | pathcount/ |
| Response Tables 1.1 = 2.9 = 3.3 (SNV, no filter) | single-tumor models (HG008T, COLO829T), HCC1395T v2.1 only | six-tumor models on three references, all 9 rows; comparator columns: same call sets, 2 F1 cells move by 0.001 (recomputed from counts) | tables/ |
| Response Tables 1.2 = 2.10 = 3.4 (INDEL, no filter) | single-tumor INDEL models; HCC1395T v2.1 only, ONT at ppr8 | six-tumor INDEL models on three references (v2.1 = mean of 3 seeds, ONT at pp); 2 comparator F1 cells move by 0.001 | tables/ |
| Response Tables 1.5 = 2.5 = 3.15 (regions, before the filter) | October single-tumor models, 12 rows; INDEL column from earlier models | six-tumor models, 18 rows (HCC1395T added), new bootstrap intervals; INDEL column from the six-tumor INDEL models | regions/ |
| Response Tables 1.7 = 2.3 = 3.6 (absorbed INDEL context) | HPRC v1.1 numbers under a v2.1 caption | the HPRC v2.1 table (1,904 / 394 / 46 absorbed), with the first-row percentages the caption promises | tables/ |
| Response Tables 2.6 = 3.8, 3.9 (false positives by germline / graph) | 2026-10-12 single-tumor models, 6 rows; two different headers | six-tumor models (GRCh38-only vs HPRC v2.1), 9 rows; one header wording for the three tables | fpdecomp/ |
| Response Table 2.7 (SNV, nearby germline) | 2026-10-12 HG008T single-tumor models | six-tumor models on three references | E19 `strata_latest.md` |
| Response Table 2.8 (INDEL, nearby germline) | v2.1 rows seed 0 only | v2.1 rows = best seed + mean of three seeds (3 rows added) | E19 `strata_latest.md` |
| Response Table 3.1 (three references, after the filter) | October single-tumor models, per seed, 18 rows | six-tumor models, 45 rows: GRCh38-only, v1.1, v2.1 and v1.1 / v2.1 without site features (the RIS models); ΔF1 with paired 95% intervals | tables/, regions/ |
| Response Table 3.2 (per seed) | INDEL only, HG008T + COLO829T, ONT at ppr8 | INDEL four- and six-tumor, three tumors, ONT at pp; plus four rows of COLO829T single-tumor SNV seeds | tables/ |
| Response Table 3.10 (graph-coordinate validation) | 541 calls of older models (all nine rows) | 935 calls of the six-tumor models; every column recomputed (assemblies, matched-normal reads, cross-platform) | offref/ |
| Response Table 3.11 (CHM13) | 396 of 541 | 690 of 935 (73.8%) | offref/ |
| Response Table 3.17 (settings) | no PASS rule; old checkpoint rule; SNV-only site features; one split for SNV and INDEL | Checkpoint = post-filter validation F1; new Calling row; Tensors row covers SNV and INDEL and the GRCh38-only models; Labels row gives the INDEL split | – |
| Response Figure 3.2 (SBS96, callers) | 2026-10-12 single-tumor models; only HG008T under this caption | three images (HG008T, COLO829T, HCC1395T), six-tumor models | sbs96/ |
| Response Figure 3.3 (SBS96, graph-coordinate calls) | wrong image (COLO829T callers' spectra), October rescue-on caption | the 935 graph-coordinate calls of Table 3.10, 3 × 3 panels | offref/ |

Unchanged: Response Figures 1.1 (caption typo fixed), 1.4 = 2.3 = 3.1, 1.5 = 2.1, 1.6 = 2.2, 2.4 = 3.4; Response Tables
1.3 = 2.4 = 3.7, 1.4, 1.6 = 2.1, 1.8, 2.2, 2.11 = 3.12, 2.12 = 3.13, 3.5, 3.14, 3.16, 3.18.

## 2. Sentences changed

| point | ver3 | ver4 |
|---|---|---|
| R1.1 point 1 | masking "reduced SNV recall in all twelve comparisons" | six-tumor models: recall lower in all nine (0.204–0.344), F1 in all nine (0.057–0.243) |
| R1.1 point 2 | "Removing channel 6 … all twelve … eleven of twelve" | four comparisons on the current pipeline (seed 0): ΔF1 0.061–0.079 on Illumina (1.9× more false calls without the channel), 0.001–0.012 on Fiber-seq; 11 of 12 seeds above the model without it |
| R1.2 point 1 | benefit "especially evident before population filtering" | same, plus "where the HPRC models also take population site features as an input" |
| R1.3 | "same … pipeline"; "all twelve"; INDEL "Illumina gained, whereas long reads lost precision"; "ONT produced the most graph-coordinate PASS calls" | pipeline differs in the site-feature input; "all eighteen"; INDEL F1 before filtering higher on HPRC in 17 of 18; graph-coordinate calls "most numerous on ONT in total" |
| R1.4 | graph-specific loss 0.2–0.7% (SNV), 3.5–7.0% (INDEL); v2.1 "84–91%" | 0.1–0.7% and 1.1–7.0% (with the HCC1395T columns of the table); v2.1 85–91% |
| R1.6 | "sensitivity below VAF 0.08" | below the candidate floor (0.05–0.08) |
| R1.7 | four-tumor models "remain below both on COLO829T Fiber-seq and ONT" | above both on HG008T Illumina / PacBio HiFi (0.623, 0.626) and COLO829T Illumina (0.736); below ClairS-TO on COLO829T Fiber-seq (0.639); on ONT 0.004 below ClairS-TO on HG008T (0.548) and below both on COLO829T (0.365); six-tumor above both on all six HG008T / COLO829T datasets |
| R2.3 point 2 | 110 (6.8%) / 27 (6.7%); "about 92%"; 1,153 of 1,351; 26 carried; 22 same-length sites (7 / 11 / 4); 37–41 (10.7–11.8%) | 120 (7.5%) / 35 (8.7%); "about 91%"; 1,148 of 1,346; 34 carried; 15 sites (7 / 4 / 4); 38–42 (11.0–12.1%) |
| R2.4 point 1 | floor 0.05 | 0.05 (0.08 for ONT and for HG008T) |
| R2.5 point 2 | COLO829T: excess false calls "mainly from SMaHT extreme regions" | before filtering, v2.1 reduces false calls and improves precision in SMaHT extreme regions in all three tumors, all platforms |
| R2.5 point 3 | "using the same pipeline on the GRCh38-only graph" | "on the GRCh38-only graph …, without the population site features" ("over 70%" kept: 71–87%) |
| R2.5 point 4 | largest gains at ≥2 nearby germline variants (+0.107 to +0.204), "mainly through higher precision"; general conclusion | larger gains with nearby germline variants (+0.083 to +0.280 against +0.068 to +0.132), higher precision and higher or equal recall; conclusion limited to HG008T |
| R2.6 point 2 | model rejects 79–88%, filter removes 74–89%, 2.4–3.2% remain, 83–98% of false calls | six-tumor models on v1.1 and v2.1 reject 97.8–98.9%; the filter removes at most 15 more per call set; 1.1–2.2% remain; 85–98% |
| R2.6 point 3 | single-tumor 0.31–0.33 / 0.61–0.67; "cannot be attributed to that filter or to overlap between its resources and the test samples"; six-tumor 0.44–0.45; "INDELs show no consistent difference" | six-tumor before the filter 0.635–0.651 / 0.744–0.802 / 0.712–0.884 against 0.013–0.093 on GRCh38-only; states that the HPRC models take population allele frequencies and the filter tag as input, so the before-filter gap is graph plus this input; gives the matched no-site-feature gain after filtering (0.038–0.295); INDELs "show the same pattern before filtering" |
| R3.1 | HG008T + COLO829T, single-tumor, three seeds; "all six datasets … 0.028–0.292"; intervals cited as Table 3.2 | three tumors, six-tumor models, matched control; HG008T / COLO829T gain 0.075–0.305 (0.038–0.295 without site features); HCC1395T within 0.03 on Illumina / PacBio HiFi, ONT v2.1 0.047 lower (0.079 without site features); intervals in Table 3.1 |
| R3.1 points 1 / 2 | "0.01 to 0.31–0.33 … 0.07–0.10 to 0.61–0.67 … cannot be attributed to downstream population filtering"; "advantage also persists after removing the haplotype-count channel" | 0.013–0.093 against 0.630–0.879 (v1.1) and 0.635–0.884 (v2.1), with the site-feature statement and the pointer to the after-filter matched comparison; the channel sentence removed (not tested on the current models) |
| R3.2 | 71–86%, 0.2–0.6%, 83–88%; INDEL 0.1–0.3% / 15–22%; Table 3.6 and 3.7 references swapped | 71–87%, 0.1–0.7%, 83–90%; INDEL ≤0.2% on Illumina, 15–29% (HG008T, COLO829T) and 0.8–2.6% (HCC1395T) on long reads; references fixed |
| R3.3 | "HPRC v1.1 and v2.1; HCC1395T: v2.1 only"; 34 tumor-only calls | HPRC v2.1, six-tumor models; 40 of 935, 27 of them with no ALT read in the matched normal |
| R3.4 | 396 of 541 (73.2%); 0.5–10.1% / 0.4–2.2% | 690 of 935 (73.8%); 3.6–11.4% / 1.8–3.9% |
| R3.7 | 6.8% / 6.7%; 92.2% / 92.5%; 1,153 of 1,351; 226 of 319 (70.8%) | 7.5% / 8.7%; 90.9% / 90.5%; 1,148 of 1,346; 226 of 311 (72.7%) |
| R3.8 | "allele fraction of at least 0.08" | 0.05 to 0.08 |
| R3.9 point 1 | SNV SD 0.001–0.042 (October), INDEL 0.003–0.027 | SNV SD 0.003–0.048 (COLO829T single-tumor current models), INDEL 0.004–0.064 (three tumors, ONT at pp); three largest spreads named |
| R3.10 points 2 / 4 | "three of the six"; long homopolymers 44–70% | four of the six (mean INDEL F1 of three seeds); 35–61% (seed means, ONT at pp) |

Captions rewritten to name the models, the model input and the PASS rule: Tables 1.1 / 2.9 / 3.3, 1.2 / 2.10 / 3.4,
1.5 / 2.5 / 3.15, 2.6 / 3.8, 2.7, 2.8, 3.1, 3.2, 3.9, 3.10, 3.11; Figures 1.2, 1.3, 3.2, 3.3. Removed: "HCC1395T has
tensors only on HPRC v2.1", "HCC1395T uses six-tumor models", the "Exploratory" epoch paragraph under Table 3.1 (October
models), "non-germline calls account for 43–57% of GRCh38-graph ONT FP" (now 13% / 49% / 6%), the rescue-on sentence of
Figure 3.3. The "Reproducibility" paragraph under Table 3.1 now says its retraining and second-cluster checks were on the
earlier single-tumor models.

## 3. Open points for you to decide (not changed in ver4)

1. **Before-filter comparisons mix graph and site features.** The HPRC models take population allele frequencies and the
   filter tag as input; the GRCh38-only models do not. Their SNV call sets barely change at the filter (0–47 false calls
   removed, against 27,661–114,970 on GRCh38-only), so the before-filter gaps in Tables 1.1 / 1.5 / 2.6 are much larger than
   the graph alone. ver4 states this in the captions and in R1.2, R2.6, R3.1. The matched no-site-feature six-tumor models
   exist only on RIS and only with after-filter numbers (used in Table 3.1).
2. **After the filter the picture is weaker than the before-filter region table shows.** SMaHT extreme ΔF1 is above zero
   in 11 of 18 pairs (17 of 18 positive); HCC1395T is level, and HCC1395T ONT on v2.1 is worse (whole chr1 −0.047
   [−0.104, −0.012]; SMaHT extreme −0.268). INDEL F1 in homopolymers ≥7 bp is higher on HPRC in 9 of 18 after the filter
   (17 of 18 before). The after-filter table is `regions/table_after_filter.md`. The 2026-11-12 choice of a single-tumor
   region table was replaced by the six-tumor one.
3. **Nearby-germline effect (R2.5 point 4) holds on HG008T only**; COLO829T Illumina is reversed (+0.209 / +0.223 with no
   nearby variant against +0.083 / +0.088 at ≥2). Tables for the other tumors: E19 `strata_latest.md` sections 8, 9.
4. **Claims in the text that the six-tumor results do not support everywhere** (left as they are): R1.6 "no significant
   loss of recall" and R2.3 "broadly comparable" recall / "leading F1" (SNV recall HG008T ONT 0.656 against 0.869 / 0.831,
   HG008T Illumina 0.776 against 0.890 / 0.845, HCC1395T ONT 0.676 against 0.766 / 0.730; F1 not the highest on HCC1395T
   ONT SNV and on HCC1395T INDELs).
5. **Path-count channel**: no six-tumor, HG008T or ONT model without the channel exists (about 50 H100-hours to train);
   the Illumina v2.1 model without the channel lies inside the seed range of the models with it.
6. **Six-tumor SNV seeds 1 / 2** are on RIS; the SNV seed SD uses the COLO829T single-tumor models. One INDEL seed is an
   outlier before the filter on Illumina (seed 1: 582 false calls on HG008T against 109 / 102), which pulls the v2.1 means
   of Tables 1.2 / 2.10 / 3.4 down; Table 3.9 uses seed 0.
7. **Four-tumor ONT SNV model**: R1.7 uses the current-recipe model (0.548 / 0.365), the one that fails on chr1; the
   original-recipe model gives 0.587 / 0.380. The letter still has no table with the four-tumor results (old R1.7b).
8. **Left from ver3, outside the data**: "[Figure]" in R3.9 point 4; R3.6 point 5 and R3.9 point 2 cite Table 3.15 for
   read placement, MAPQ and comparator callers, which it does not contain; R1.6 cites Table 1.7 as a recall analysis; empty
   numbered list items after Figure 1.3, Table 1.3 and Table 3.10; old labels (R1.minor2, R1.m2, R3.m1, "Supplementary Note
   R2.3"); R3.m2 "Training data" and R3.10 "for both tumors" still describe the two-tumor design; Table 1.8 caption says
   SNV models have one seed; R2.2 and R2.3 count normal-carried alleles with different rules (0 against 1 COLO829T SNV).
9. **Candidate-level items** (Figures 1.4 = 2.3 = 3.1, 1.6 = 2.2; Tables 1.3 = 2.4 = 3.7, 3.5) are on the AF 0.08
   candidate sets, while the SNV models now use AF 0.05 sets for COLO829T and HCC1395T on Illumina / PacBio HiFi.
10. **Main-text Figure 4a** must be redrawn from `sbs96/sbs96.tsv` to match Response Figure 3.2.

## 4. What was computed for ver4

| item | what ran | check |
|---|---|---|
| masking (Fig. 1.2) | masked predict of the three six-tumor SNV models on chr1 of the nine sets (three 12gb GPU jobs, 28–48 min), VCF and scoring | the same code without the mask on HG008T PacBio HiFi gives TP / FP / FN 567 / 481 / 130 against 567 / 479 / 130 in the existing report: per-tensor probabilities differ by at most 0.005 between GPU types and 6 of 140,610 cross the threshold |
| graph-coordinate calls (Tables 3.10, 3.11, Fig. 3.3) | context extraction, assembly alignment, matched-normal read counts (incl. a GIAB HG008-N fetch), CHM13 projection for the 935 calls; CPU, about 1 h | the copied scripts reproduce all 9 rows of the ver3 Table 3.10 on the ver3 inputs |
| false positives (Tables 2.6 / 3.8, 3.9) | the ver3 script on the six-tumor call sets; CPU minutes | reproduces one October cell |
| regions (Tables 1.5 …, Table 3.1 intervals) | block bootstrap on existing call sets; 40 s | point values equal `report.json` |
| SBS96 (Fig. 3.2) | counting from existing call sets | comparator cosines equal the October file |
| R2.6 point 2, recounts (R2.3, R3.9, R3.10) | statistics on existing predictions and tables | old values reproduced first |

An independent check of the test build against the raw evaluation reports (8 checkers, 89 findings) was run before the
final build; its corrections are in sections 1–2 and its open points in section 3.

## 5. Notes on the Word file

- Edits are plain (untracked) edits; ver3's tracked changes and its six comments are kept. Two exceptions: the six HCC1395T
  rows added to each of Tables 1.5 / 2.5 / 3.15 are copies of rows that ver3 holds as tracked insertions, so they appear
  as insertions dated 2026-11-10; and four small tracked insertions of Li, Daofeng (">=7bp", "regions" in R1.3, "mainly
  from" in R2.5, and one more) were inside sentences that ver4 rewrote, so they are gone.
- Table 1.2 keeps ver3's yellow highlight and tracked insertion; its copies 2.10 and 3.4 have neither.
- Response Figure 3.2 is now three images of 6.5 × 3.9 in (about one and a half pages).
- Not opened in Word here (no Word on the cluster); the package was checked for well-formed XML, unique ids and
  relationships, and read back with pandoc.
