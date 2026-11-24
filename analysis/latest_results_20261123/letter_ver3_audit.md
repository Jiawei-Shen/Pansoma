# Reviewer letter ver3 (Reviewer_Letter_Response_ver3.docx, saved 2026-11-23 15:37) against ver2 and the latest results

Checked: the ver3 Word file converted with pandoc (tracked changes accepted), compared paragraph by paragraph and word by
word with `Reviewer_Letter_Response_ver2_clean.docx`, and its 9 embedded images compared by checksum. Model results quoted
below are chr1 of the test tumor, benchmark regions, rtg vcfeval, PASS calls after the population filter, read from the
runs' `eval_bed/report.json`. The replacement tables on the current models are in `letter_ver2_audit.md` (same folder);
section numbers below (§) refer to that file. A Word lock file was present, so edits made after 15:37 are not covered.

## 1. What changed from ver2 to ver3

No number changed in any table that ver2 already had. The 9 images are byte-identical to ver2's and in the same order.

| change | detail | check |
|---|---|---|
| Renumbering | every figure / table is now "Response Figure x.y" / "Response Table x.y" (mapping in §2 below) | – |
| Wording | "PoN" / "panel of normal" -> "population filter" / "resource VCFs" in text, headers and captions | – |
| Added: Response Table 3.1 | per-seed chr1 SNV results on three references (the old R3.9a) | equals `supplementary/R3.9/three_ref_table.md`; October single-tumor models (see 4.5) |
| Added: Response Table 3.2 | per-seed chr1 INDEL results of the four-tumor and six-tumor models | equals the runs at ppr8 on every platform, ONT included (see 4.1) |
| Copies added | R2.6a as Response Tables 1.1 and 3.3 (was only in R2.6); R2.6b as 3.4; the candidate figure as Response Figure 3.1; the alignment figure and table as Response Figure 3.4 / Table 3.12 | same numbers as the originals |
| Removed | R1.m1: the paragraph with the four-tumor / six-tumor F1 values (0.845 / 0.820 / 0.658 ...) and the CASTLE paragraph (24.5 %, 4.6 %, 99.1–99.7 %); the "Calling" row of the settings table; the post hoc epoch paragraph below it (its content is now in the caption of Response Table 3.1); every citation of a table that was not in the letter (R1.1(c), R1.3b, R1.7b, R1.7c, R3.3b, R3.4a, R3.10, Fig. R1.7); "Supplementary Table [X]" placeholders; "[待补: ...]"; the highlighted "update the data" note | – |
| Changed values | HPRC v2.1 haplotypes: "456" (R1.2) and "460" (R3.1) -> 464; settings table: allele fraction ≥0.05 for SNVs (≥0.08 for HG008T and all ONT datasets), ≥0.08 for INDELs; "allele quality" -> "base quality"; Tensors row: "SNV models also take population site features"; Training row: "other" examples below allele fraction 0.1 sampled at 25 % | correct |
| Rewritten: R1.m2 | new example text: 8-bp deletion over "six graph nodes of 1–2 bp, carried by 19 to 78 of the 464 haplotypes"; SNV example 18 / 72 of 90 | correct (`supplementary/R1.m2/indel/README.md`: 20, 19, 27, 51, 78, 77 of 464) |
| Rewritten: R3.m1 | new sentence: the main Methods listed the thresholds of the original candidate-count comparison, the Supplementary Methods those of tensor generation | not checked against the original submission |
| Rewritten: R3.m2 | software / population-resource paragraphs shortened; "HG008T ONT reads were basecalled with dorado SUP v5.0" added | not rechecked (no model result) |

## 2. Status of every item in ver3

"old models" = not the current pipeline (definition at the top of `letter_ver2_audit.md`). HG008T single-tumor rows count as current (not retrained).

| ver3 label | ver2 label | status in ver3 | replacement |
|---|---|---|---|
| Response Figure 1.1 | Fig. R1.minor2 | OK; caption typo `TITATATA>T` still there (figure: TTATATATA>T) | – |
| Response Figure 1.2 | Fig. R1.1(a) | old models (COLO829T rows); text still says "all twelve comparisons", the figure shows six | none exists (§3.1) |
| Response Figure 1.3 | Fig. R1.1(b) | old models (COLO829T rows); caption label fixed; text still "all twelve ... eleven of twelve" and "Removing channel 6" | §3.1 |
| Response Figure 1.4 = 2.3 = 3.1 | Fig. R1.2 | OK; now cited as a figure | – |
| Response Table 1.1 = 2.9 = 3.3 | Table R2.6a | old models (COLO829T, HCC1395T rows); "HCC1395T has tensors only on HPRC v2.1" | §3.2 |
| Response Table 1.2 = 2.10 = 3.4 | Table R2.6b | HCC1395T rows: only v2.1, ONT at ppr8 | §3.3 |
| Response Table 1.3 = 2.4 = 3.7 | Table R1.4 | OK (two rounding cells, AF 0.08 sets) | – |
| Response Table 1.4 | Table R1.3 SMaHT | OK | – |
| Response Table 1.5 = 2.5 = 3.15 | Table R1.3a | old models in the four COLO829T Illumina / Fiber-seq rows | §3.4 |
| Response Table 1.6, 2.1 | Table R1.5, R2.1 | OK | – |
| Response Figure 1.5 = 2.1 | Fig. R1.6 | OK, current (six-tumor models, INDEL ONT at pp) | – |
| Response Figure 1.6 = 2.2 | Fig. R1.6a | OK | – |
| Response Table 1.7 = 2.3 = 3.6 | Table R1.6b | STILL WRONG: HPRC v1.1 numbers (2,277 / 430 / 52) under an HPRC v2.1 caption | §2.1 |
| Response Table 1.8 | Table R1.7a | OK (caption: "all other models use every example" not true for the current ONT SNV models) | – |
| Response Table 2.2 | Table R2.2 | OK | – |
| Response Table 2.6 = 3.8, 3.9 | Table R3.2b, R3.2c | old models (2026-10-12 single-tumor runs) | none exists (§3.7) |
| Response Table 2.7 | Table R2.5c | old models | §3.5 |
| Response Table 2.8 | Table R2.5d | current, HPRC v2.1 rows seed 0 only; caption does not name the model | §3.5 |
| Response Figure 2.4 = 3.4, Table 2.11 = 3.12 | Fig. / Table R2.7 | OK | – |
| Response Table 2.12 = 3.13, 3.14 | Table R3.6, R3.6a | OK | – |
| Response Table 3.1 | (new; old R3.9a) | equals its source; October models, see 4.5 | §3.8 |
| Response Table 3.2 | (new) | equals the runs; ONT rows at ppr8, see 4.1 | 4.1 below |
| Response Table 3.5 | Table R3.2a | OK | – |
| Response Table 3.10, 3.11 | Table R3.3a, R3.4b | old models (COLO829T, HCC1395T rows) | counts only (§3.7) |
| Response Figure 3.2 | Fig. R3.5a | old models; only the HG008T half sits under this caption | none exists |
| Response Figure 3.3 | Fig. R3.5b | STILL WRONG IMAGE: the COLO829T half of the callers' spectra; caption (6,025 calls) is the October rescue-on output | §2.2 |
| Response Table 3.16 | Table R2.3 | table OK; R2.3 / R3.7 text still has the earlier counts | §2.3 |
| Response Table 3.17 | R3.m1 settings table | partly updated, see 3 and 4.4 | §3.6 |
| Response Table 3.18 | R3.m2 comparator table | OK | – |

## 3. Findings of the ver2 audit: fixed or still open

**Fixed in ver3**

- HPRC v2.1 haplotype count (464 in R1.2 and R3.1).
- Caption label of the path-count figure ("R1.1(c)").
- The candidate figure cited as "Supplementary Table R1.2".
- R1.m1 values of the older four-tumor / six-tumor models (paragraph removed).
- Citations of tables that were not in the letter (removed; R3.9a / R3.9b are now Response Tables 3.1 / 3.2, but see 4.2).
- Placeholders "[待补]", "Supplementary Table [X]", the "update the data" note.
- Channel numbering in the R1.m2 text ("panels 2 and 3 of Response Figure 1.1").
- Settings table: candidate allele-fraction floor, site features, 25 % sampling.

**Still open**

| item | where in ver3 | what is wrong | correct / latest |
|---|---|---|---|
| absorbed-INDEL table | Response Tables 1.7, 2.3, 3.6 | v1.1 numbers under a v2.1 caption | v2.1: 1,904 (22.4%) / 394 (20.6%) / 46 (2.8%); HP + STR 93.3% / 73.1% / 58.7% (§2.1) |
| graph-coordinate SBS96 figure | Response Figure 3.3 | wrong image; caption from the October rescue-on output | `supplementary/R3.5/fig_R3.5b_graph_coordinate.png` exists only for that output (§2.2) |
| normal-haplotype counts in the text | R2.3 point 2, R3.7 | 110 (6.8%), 27 (6.7%), 92.2% / 92.5%, 1,153 of 1,351, 226 of 319 (70.8%) | Table 3.16: 120 (7.5%), 35 (8.7%), 90.9% / 90.5%, 1,148 of 1,346, 226 of 311 (72.7%) (§2.3) |
| "all twelve comparisons" | R1.1 (masking; path-count channel) | the figures show six | – |
| "Removing channel 6" | R1.1 point 2 | the figure numbers the channel 3 | – |
| caption typo | Response Figure 1.1 | `TITATATA>T` | `TTATATATA>T` |
| allele-fraction floor in the text | R3.8 "at least 0.08"; R2.4 "about 0.05" | two values | 0.05 (SNV, Illumina and PacBio HiFi / Fiber-seq of COLO829T and the four lines), 0.08 (HG008T, ONT, INDEL) |
| settings table, Checkpoint row | Response Table 3.17 | "highest validation somatic F1" | highest validation F1 after the population filter (§3.6) |
| placeholder | R3.9 point 4 | "[Figure]" | – |
| "HCC1395T has tensors only on HPRC v2.1" | captions of Response Tables 1.1, 1.2, 2.9, 2.10, 3.3, 3.4 | HCC1395T has GRCh38-only and v1.1 sets and six-tumor models | §3.2, §3.3 |
| six-tumor model before the filter on HG008T | R2.6 point 3 | "0.44–0.45" | 0.635 / 0.651 / 0.640 |
| COLO829T single-tumor before the filter | R2.6 point 3, R3.1 point 1 | "0.07–0.10 to 0.61–0.67" | §3.2 |
| HPRC against GRCh38-only | R3.1 | "+0.028 to +0.292 on all six datasets" | COLO829T level with site features (§3.8) |
| six-tumor INDEL against the comparators | R3.10 point 2 | "three of the six" | 4 of 6 (§3.8) |
| four-tumor model on COLO829T | R1.7 | "below both on COLO829T Fiber-seq and ONT" | Fiber-seq 0.639: above DeepSomatic-TO (0.632), below ClairS-TO (0.762) |
| graph-coordinate fractions | R3.3, R3.4 | 541 calls, 34 tumor-only, 396 CHM13, 0.5–10.1 % / 0.4–2.2 % | older COLO829T / HCC1395T models (§3.7) |
| empty numbered list items | after Response Figure 1.3, Response Table 1.3, Response Table 3.10 | "1. 2. 3." without text | – |

## 4. New problems introduced in ver3

### 4.1 Response Table 3.2 gives ONT at ppr8; Response Figures 1.5 / 2.1 plot the same models at pp

The ONT rule was changed to pp on 2026-11-13 for the four-tumor and six-tumor INDEL models, and the PR figure uses it. The
table's ONT rows (and its caption, "validation recall 0.8 after filtering") are ppr8. Example: six-tumor HG008T ONT, table
P / R 0.192 / 0.150, figure dot 0.522 / 0.080. ONT rows at pp (seeds 0 / 1 / 2; mean in brackets):

| tumor | platform | training set | precision | recall | F1 | F1 in Table 3.2 (ppr8) |
|---|---|---|---|---|---|---|
| HG008T | ONT | 6 tumors | 0.522 / 0.524 / 0.548 | 0.080 / 0.075 / 0.078 | 0.139 / 0.131 / 0.137 (0.136) | 0.168 / 0.178 / 0.157 (0.168) |
| HG008T | ONT | 4 tumors | 0.574 / 0.517 / 0.562 | 0.060 / 0.053 / 0.061 | 0.108 / 0.096 / 0.111 (0.105) | 0.139 / 0.143 / 0.146 (0.143) |
| COLO829T | ONT | 6 tumors | 0.196 / 0.204 / 0.290 | 0.084 / 0.076 / 0.069 | 0.118 / 0.111 / 0.111 (0.113) | 0.109 / 0.088 / 0.074 (0.090) |
| COLO829T | ONT | 4 tumors | 0.071 / 0.056 / 0.105 | 0.084 / 0.069 / 0.061 | 0.077 / 0.061 / 0.077 (0.072) | 0.069 / 0.054 / 0.085 (0.070) |

Illumina and PacBio HiFi / Fiber-seq rows of Table 3.2 are current. If ONT stays at ppr8 in this table, the caption has
to say that the figure uses another threshold for ONT; R3.9's "0.003–0.027 (median 0.014)" SD range was computed with
the ppr8 rows.

### 4.2 R3.1 cites Response Table 3.2 for the paired confidence intervals

"SNV results after population filtering are reported in Response Table 3.1, with paired confidence intervals in Response
Table 3.2": Response Table 3.2 is the INDEL per-seed table. The paired table (old R3.9b,
`analysis/reviewer_response_20261011/supplementary/R3.9/paired_table.md`) is not in the letter.

### 4.3 The settings table now contradicts the captions of the no-filter tables

Response Table 3.17 says "SNV models also take population site features of the full graph as the input". The captions of
Response Tables 1.1 / 2.9 / 3.3 say "Pansoma before our population filter (it has none of its own)", and Response Tables
1.2 / 2.10 / 3.4 say "unlike the SNV models, their values without the filter are not free of population information".
Both captions are true only for the October models that these tables contain.

### 4.4 No PASS rule in the settings table

The "Calling" row was deleted. The current rules (SNV: validation recall 0.95 of the callable positives after the
population filter; INDEL: validation recall 0.8 after the filter on Illumina / PacBio HiFi, the post-filter validation
best-F1 threshold on ONT; off-reference rescue off) are not stated in the table, and the only PASS rule in the letter's
tables is "validation recall 0.9 before filtering" in the caption of Response Table 3.1, the old rule.

### 4.5 Response Table 3.1 is the October models

It agrees with the R3.1 text and with Response Table 3.3, so the R3.1 block is internally consistent, but it describes
single-tumor models at AF 0.08 without site features and with PASS at validation recall 0.9, while Response Table 3.17
describes the current pipeline. On the current pipeline the COLO829T rows are (§3.8 of `letter_ver2_audit.md`):

| COLO829T chr1 SNV, current pipeline | seeds 0 / 1 / 2, F1 | Table 3.1 (October), F1 |
|---|---|---|
| Illumina, HPRC v2.1 | 0.797 / 0.725 / 0.816 | 0.783 / 0.704 / 0.765 |
| Illumina, HPRC v1.1 | 0.774 / 0.714 / 0.722 | 0.572 / 0.626 / 0.566 |
| Illumina, GRCh38-only with site features | 0.773 / 0.752 / – | 0.468 / 0.386 / 0.530 (no site features) |
| Illumina, GRCh38-only without site features | 0.501 (seed 0) | |
| Fiber-seq, HPRC v2.1 | 0.818 / 0.820 / 0.823 | 0.779 / 0.774 / 0.789 |
| Fiber-seq, HPRC v1.1 | 0.819 / 0.819 / 0.825 | 0.752 / 0.782 / 0.771 |
| Fiber-seq, GRCh38-only with site features | 0.816 / 0.820 / 0.812 | 0.680 / 0.700 / 0.692 (no site features) |
| Fiber-seq, GRCh38-only without site features | 0.718 (seed 0) | |
| ONT, HPRC v2.1 | 0.754 (seed 0) | 0.737 / 0.718 / 0.730 |

HG008T rows of Table 3.1 are current (models not retrained).

### 4.6 R1.7 compares models with the comparators without a table

The R1.7 text still says the four-tumor models exceed both comparators on HG008T and stay below both on COLO829T
Fiber-seq and ONT, and that the six-tumor models are above both on HG008T and COLO829T; the table it relied on (R1.7b)
is no longer cited or included. Latest R1.7b and R1.7c: §3.8 of `letter_ver2_audit.md`.

### 4.7 Leftover labels of the old numbering

"R3.1" and "After filtering: R3.1 (point 3, single-tumor models)" in the captions of Response Tables 3.3 / 3.4; "(R3.5)"
and "in R3.1" in the caption of Response Table 3.1; "(versions as in R3.m1)" in the methods under Response Table 3.12;
"(R1.m2)" in R1.8 point 3.

## 5. Word comments in ver3 (not data)

"may remove this sentence?", "what is this?", "a bed file provide by the smaht consortium", "do you have clear definitions
for easy, difficult, extreme regions?", "Yes", "now at Cell Genomics".
