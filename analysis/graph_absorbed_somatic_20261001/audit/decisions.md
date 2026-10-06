# Audit decisions: s9 integration of graph-absorbed HG008T somatic truths (2026-10-01)

> Record of the audit-fix step. The numbers below were taken from that s9 run (Slurm 379498, Illumina read paths only,
> before the chr1-22 fix of the PoN control group). Categories did not change afterwards; read-derived numbers (path
> choice, closest sub-reasons, truth-read flags) and the PoN control percentages did. Current numbers: `../README.md`,
> `../tables.md`. Paths: `$A` = this analysis folder, `$D` = `/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001`,
> `$T` = `tmp/graph_absorbed_somatic_20261001` (its `s9.sbatch` = `$A/s9_integrate.sbatch`). The audit scripts and
> outputs are in `$D/audit/`.

Four audits checked the s9 integration: hg008n (normal-hap membership), hprc (HPRC membership and frequency), logic (rules and totals) and interpretation (biology). This file records the decision on every finding and the numbers after the fixes.

## What changed

- **New step s8b_array_alleles.py** (`$A`; output `$D/array_alleles.tsv`, log `$T/s8b_array_alleles.login.log`; ran on the login node, 10 min, 72 MB).
  - It grows the whole tandem array around each truth: the s1 event core plus every exact periodic stretch of period 1-60, including adjacent and compound arrays.
  - It picks common unique GRCh38 24-mer anchors that are present once in every placed contig (HG008-N hap1/hap2 and the HG008-T contigs).
  - It extracts each haplotype's sequence between those anchors and compares it with R (GRCh38) and with A (GRCh38 with the truth applied).
  - It rebuilds each normal hap independently from the phased dipcall records as a second source.
  - It compares the tumor contigs with the normal haps.
  - It applies the truth INFO event (`asm_event`, normal-assembly frame, 0-based) to the event hap's sequence (`T_info`).
  - Result: all 2,478 loci have GRCh38 anchors. The assembly and dipcall sequences are identical for 4,808 haplotypes and differ for 57; only one of those differences changes an ALT call (17053 hap1).
- **s9_integrate.py rewritten in place.** The docstring gives the method. A copy of the old script is in `audit/s9_integrate.before_fixes.py`, and the old outputs are in `$D/before_fixes/`.
- **Final run:** Slurm 379498, 4 min 42 s, MaxRSS 9.1 GB, 12G requested. Log: `$T/s9_integrate.379498.log`.
- **Read data:** Illumina only. PacBio rest jobs (379468) were at 50/122 loci per task at 21:53. ONT (379469) has 3 tasks running and 17 pending. The script picks up any `read_path_summary_P.tsv` once `s3_read_paths.py merge P` has run, and is re-run unchanged with `sbatch $T/s9.sbatch`.

## Final numbers (union of perfect sets; chr1-22 = all loci, no chrX/Y)

Each cell reads present_broad / present_rare / absent_HPRC_other / other_ambiguous (total).

| set | SNV chr1-22 | INDEL chr1-22 |
|---|---|---|
| PacBio | 2 / 3 / 129 / 36 (170) | 41 / 17 / 1,157 / 148 (1,363) |
| ONT | 2 / 3 / 117 / 26 (148) | 26 / 10 / 805 / 148 (989) |
| Illumina | 2 / 2 / 97 / 21 (122) | 68 / 28 / 1,682 / 138 (1,916) |
| union | 2 / 3 / 150 / 46 (201) | 75 / 30 / 1,934 / 238 (2,277) |

Other scopes for the union:

| scope | SNV | INDEL |
|---|---|---|
| chr2-22 | 2 / 3 / 138 / 41 (184) | 66 / 28 / 1,781 / 217 (2,092) |
| chr1 BED | 0 / 0 / 12 / 5 (17) | 8 / 2 / 140 / 19 (169) |
| GIAB nogermlineoverlap BED | 2 / 2 / 64 / 3 (71) | 43 / 14 / 840 / 6 (903) |

Before the fixes, the union was SNV 25 / 5 / 150 / 21 and INDEL 105 / 37 / 2,043 / 92.

**HG008N_absent_HPRC_other, split by HPRC support (SNV / INDEL):**
- exact_allele: 65 / 1,788. At least one complete HPRC haplotype spells the exact ALT allele over the array.
- element_set_other_allele: 76 / 101. HPRC haplotypes walk all the somatic elements but spell another allele, usually another repeat length.
- recombinant_pieces: 9 / 45. No single haplotype walks a whole ALT path.
- Exact-allele frequency median: SNV 0.0, INDEL 0.167. Set level (any ALT path): 0.322 / 0.193. Element level: 0.320 / 0.197.

**HG008N_present (110 loci = 5 SNV + 105 INDEL):**
- Pattern: other_hap_only 100, other_hap_only:patient_frame 9, no_event_hap 1.
- dipcall says "identical germline allele" at 110/110.
- Every tumor contig equals the carrying hap at 101/110. 9 loci also have a tumor copy that equals neither normal hap.
- The present classes broad/rare use the exact-allele frequency. The class agrees across the element, set, exact and full-VCF frequencies at 92/110.

**other_ambiguous (SNV / INDEL):**
- closest:read_allele_not_truth 27 / 148
- closest:read_allele_nearer_truth 0 / 16
- closest_no_read_element 12 / 23
- with_germline:germline_allele_off_event_hap 3 / 17
- no_HPRC_carrier 1 / 14, plus no_HPRC_carrier:CHM13_carries 0 / 6
- HG008N_hap_unresolved: array_NA 2 / 6, conflict 1 / 0
- evidence_conflict:germline_ALT_vs_truth_INFO_event_novel 0 / 6
- evidence_conflict:ALT_on_both_normal_haps 0 / 1
- evidence_conflict:ALT_on_normal_event_hap 0 / 1

**How loci moved:**
- Of the old 172 present loci, 109 stay present (5255 makes 110).
  - 46 went to absent: 23 SNV and 23 INDEL, all failing the whole-array test.
  - 17 went to other_ambiguous: 6 INFO event novel, 3 closest, 3 array_NA, 2 off-phase germline, 1 both haps, 1 event hap, 1 assembly-vs-dipcall conflict.
- Of the old absent loci, 158 went to other_ambiguous: 143 closest, 12 off-phase germline and 3 with no HPRC carrier.
- 4 old other_ambiguous loci moved out: 5255 to present, and 8980, 12442 and 13628 to absent. The rule `evidence_conflict:dipcall_identical_vs_assembly` was dropped.

## Decisions per finding

### hg008n audit (normal-hap membership)

| # | sev | finding | decision |
|---|---|---|---|
| H1 | major | 15 present INDEL loci wrong at the repeat-array level (865 887 1666 2590 4470 4915 5570 5957 6810 7285 8485 10203 11302 13633 16861) | **Fixed.** s8b compares whole arrays, and s9 takes the array call: a hap is ALT only if its sequence between the common anchors equals GRCh38+truth, or that plus its own dipcall records outside the array (ALT_plus_flank, 7 haps). I re-checked 865, 887, 1666, 4470, 5570, 8485, 13633 and 16861 by hand (sequences in s8b test mode): none of their normal haps equals A over the array. All 15 are now out of present: 13 absent (7 exact_allele, 4 element_set_other_allele, 2 recombinant_pieces) and 2 other_ambiguous (2590 and 10203 are s2 closest). |
| H2 | major | 24 present loci (16 SNVs) are germline only at the coordinate; the tumor allele is novel | **Fixed.** The array rule moves 18 to absent. 3 go to other_ambiguous: 7731 (off-phase germline), 12872 (INFO event novel) and 17053 (assembly and dipcall disagree). 3 stay present: 30, 9593 and 15386. At those three, GIAB's own INFO event applied to the event hap gives exactly the other hap's array allele (`T_info_class eq_ALT=other_hap`). The tumor assembly has an additional novel copy, flagged in `tumor_verdict` (novel_*). For 30, the normal haps are A13 / A14 and the tumor copies are A14 / A15, which an A15 assembly error would explain. Present SNVs: 30 before, 5 now. The auditor's D-class SNVs 8124, 8415, 9872, 9873, 10913, 11941, 12130, 15656 and 15658 also fail the array test. I read 8124, 8415 and 12130: the carrier's array has another repeat count, the INFO event gives a novel tumor allele, and the tumor contig is novel. The auditor checked D only automatically, with a weaker SNV test (site base), so I keep my result. |
| H3 | major | event-hap patterns 'both' / 'event_hap_only' unreliable; LOH reading | **Fixed.** Patterns now come from the array calls. 'both' (1 locus, 7473) and 'event_hap_only' (1 locus, 10740) are other_ambiguous with evidence_conflict sub-reasons, not present. The old 7 'both' and 12 'event_hap_only' are dissolved. The interpretation no longer says LOH; see B1. `tumor_verdict` records whether every tumor copy equals the carrying hap (101/110). |
| H4 | minor | 5255 missed carrier | **Fixed.** No normal hap equals GRCh38+truth: hap1 is +13, with the 12-bp insertion plus a 1-bp insertion that both haps carry inside the array. However, the INFO event applied to the event hap (hap2) gives exactly hap1's array allele, and the tumor contig equals hap1. New pattern `other_hap_only:patient_frame` (9 loci): present_rare. |
| H5 | minor | 5 unverifiable present loci (2353 2356 2363 2364 VNTR, 7473 paralog) | **Fixed.** 2356, 2363 and 2364: no common anchors on the normal haps and outside the dip BED, so HG008N_hap_unresolved:array_NA. 2353 is now closest (read allele not truth). 7473 is ALT on both haps of the copy s6 placed; the INFO event is outside that stretch (`T_info_class event_outside_stretch`), so it goes to evidence_conflict:ALT_on_both_normal_haps. All 5 are other_ambiguous. |
| H6 | minor | HG008T_any_ALT under-calls | **Fixed.** It is renamed `HG008T_any_ALT_narrow` and no longer used in any rule or interpretation. New columns: `tumor_any_ALT_array` (a tumor copy equals A over the array; strict) and `tumor_verdict` (a tumor copy not in the normal = novel_*: 2,079/2,084 absent loci). |
| H7 | note | absent robust; 24 loci with no visible tumor change | **Noted.** Column `tumor_shows_no_change`: 3 absent and 4 ambiguous loci over the whole array. The audit's 24 used compound-run windows. For example, truths 4-7 sit in one somatic cluster, where the array does change. |
| H8 | note | dipcall phase hap1\|hap2 correct; dipcall 'identical' not independent | **Accepted.** dipcall is used as phased haplotype sequence (s8b reconstruction), not as record identity. The rule `evidence_conflict:dipcall_identical_vs_assembly` was dropped. `dipcall_relation` is kept as an annotation. |
| H9 | note | asm_event is 0-based; it can name a neighbouring event | **Accepted and used** (s8b T_info). REF matched at 0-based for every event inside its stretch (0 event_ref_mismatch). 27 events lie outside the anchored stretch (event_outside_stretch) and 30 loci have none. |
| H10 | note | stale chr2-22 numbers in the integration text | **Fixed.** All numbers here come from summary.json / per_variant.tsv. |
| H11 | note | limits of the audit method | **Noted.** |

### hprc audit (HPRC membership and frequency)

| # | sev | finding | decision |
|---|---|---|---|
| P1 | major | element-level frequency counts node walkers, not allele carriers | **Fixed.** RARE_AF (0.20) and the absent split now use `exact_allele_freq`: complete haplotypes whose bracketed sequence over the whole array equals the ALT allele, or truth+germline for with_germline. Element-level and set-level frequencies are kept as columns. I verified truth 378 myself: element 79/81, exact 10/81. Absent SNVs: 85/150 have no exact HPRC carrier and fall in element_set_other_allele or recombinant_pieces. |
| P2 | major | 143 closest absent loci have carriers of a different allele | **Fixed.** Every s2 closest locus without a truth-spelling read is other_ambiguous; no closest locus has one. That is 226 loci (39 SNV / 187 INDEL). Sub-reasons: read_allele_not_truth (175; at 1 locus the kept set could not be applied, so it gets the default), nearer_truth (16) and no reads (35, of which 5 have an s2 primary nearer the truth: `closest_kind s2_primary:*`). |
| P3 | major | representation dependence; set_carriers 0 | **Fixed.** `set_carriers_any_path` is the union over every candidate ALT path. hprc_support `recombinant_pieces` (absent: 9 SNV / 45 INDEL) means every element has carriers but no haplotype walks a whole path; these are shown separately in all tables. 14267: any-path 7, exact 2/23; it is now present_rare (patient frame). 1105: absent, recombinant_pieces. |
| P4 | major | partial-haplotype exclusion biases frequency upward | **Verified and noted.** I recomputed it myself from the full-VCF genotypes against the s4 status: carrier fraction is 0.246 among d9-complete haplotypes vs 0.117 among partial ones (absent category 0.238 vs 0.101). Denominators are unchanged, but `exact_allele_freq_all88` and `set_freq_any_path_all88` (lower bounds), `full_vcf_AF` / `d9_vcf_AF` and `freq_class_agree` were added. The class differs across definitions at 18/110 present and 471/2,084 absent loci. The table note calls these frequencies conditional and upper-biased. |
| P5 | minor | 'carriers follow the panel' understates category-specific deviation | **Fixed.** tables.md has an observed/expected table with z per category. Absent exact allele: AFR 1.108 (z +21.4), AMR 0.891 (z −16.7), EAS 0.885 (z −7.4), SAS 0.917 (z −2.5). Present_broad: AFR 0.916 (z −5.9), AMR 1.098 (z +5.3). The z values ignore pairing and linkage. |
| P6 | minor | d9 floor is ≥7 in the final graph; 'rare' = near the floor; present with 0 carriers | **Noted.** The table note gives the ≥7 figure from the audit's whole-GFA scan. 'Rare' now means an exact-allele frequency below 0.2. The 4 present_rare loci without an exact HPRC carrier show hprc_support none / recombinant / element_set. |
| P7 | note | with_germline compounds | **Noted.** element_alone keeps the per-element class. For with_germline loci the exact allele is truth plus the patient's germline alleles. |
| P8 | note | membership extraction correct | No change. |

### logic audit (rules and totals)

| # | sev | finding | decision |
|---|---|---|---|
| F1 | critical | 129 closest absent loci + 2 present encode germline / near-GRCh38 | **Fixed** (see P2). All of them are other_ambiguous now. Closest elements come from the majority of truth-spelling reads, and only those loci can be categorised (0 occur). The kept set is applied to GRCh38 and Levenshtein-compared with the truth and non-truth haplotypes (`closest_lev_truth` / `closest_lev_nontruth`, `closest_kind`). |
| F6 | major | Illumina 'ALT' read walks that spell germline only | **Fixed as columns and flags; the category is not changed.** Each read walk gets a class: alt / alt+germline / ref / germline_only / other. alt+germline means the 6 nearest PASS alleles, or every PASS record of one normal hap within the span. Verified at 509: of 68 reads, 24 alt and 44 germline_only. Path choice and `{P}_truth_frac_contain_set` use truth-spelling reads, and `{P}_n_truth_reads`, `{P}_n_germline_only_reads` and `{P}_all_reads_majority_class` were added. The locus flag `truth_read_support` is yes / no_truth_spelling_read / no_perfect_platform_read_data. Absent loci: no_truth_spelling_read for 25 SNV / 65 INDEL; no perfect-platform reads (Illumina only for now) for 53 / 252. Illumina-perfect loci where the commonest walk of all good reads spells germline only: 8 SNV + 160 INDEL (110 absent, 57 ambiguous, 1 present); the audit counted 175. A walk that spells both the truth on one hap's background and a germline allele counts as alt+germline here (sequence-identical; the audit's a7 ranked germline_only first). **Why the category stays:** for exact / with_germline loci the graph spells the truth (s2), so a truth-carrying read walks it without edit. The flag says whether the tumor reads in hand show that. |
| F2 | major | 10 present loci rest on ALT calls that contradict hap_diffs | **Fixed.** The narrow call now accepts minimap2 ALT only when the hap's hap_diffs (all, or one or two) give alt_hap, and giraffe-only ALT likewise; spans_agree False gives unresolved. The narrow call is kept as `{h}_narrow_carries_alt` and only confirms a 'no' when the array has no sequence (2 haps). Narrow vs final: ALT→no 10, ALT_plus_other→no 43, ALT→unresolved 5. All 10 listed loci are out of present. |
| F4 | major | 63 absent loci with set_carriers 0 (graph recombination) | **Fixed as a sub-class, kept in absent.** hprc_support `recombinant_pieces` uses the union over all ALT paths. The graph does spell the truth and the patient does not carry it, so these are still pangenome-induced; the mechanism differs and is reported separately (9 SNV / 45 INDEL). |
| F5 | major | 32 absent SNVs whose base is deleted on the event hap | **Fixed in normal_allele_class, category kept.** New class `SNV_site_deleted_on_event_hap` (s6 site_call of the event hap): 30 absent SNVs + 8 ambiguous. The column `asm_event_kind` shows the truth INFO event kind; see B3. The somatic allele as a sequence is still in the graph and, for exact_allele, in HPRC, so the category stands. The interpretation table names the class. |
| F-spans | minor | spans_disagree ignored | **Fixed** in the narrow call. The array call uses common anchors (a wrong copy loses them) and dipcall, which places segdup copies genome-wide. |
| F-numbers | minor | wrong numbers in the report text | **Fixed** (numbers above). |
| F-phase | minor | with_germline germline alleles phased to the non-event hap | **Fixed.** Phase-consistent candidate paths come first. If none exists, the locus is other_ambiguous `with_germline:germline_allele_off_event_hap` (3 SNV / 17 INDEL, 12 of them previously absent, plus 2585/2586). At 1887 the 12 Illumina reads spell alt+germline, but s8b shows the tumor event copy is exactly A, without hap1's +AT. |
| F-wording | note | 'event_hap_only = truth conflict' too strong | **Fixed.** The interpretation text now explains the conflict, naming the INFO event as the change relative to the normal. |
| F-small | note | tie rule, RARE_AF borderline, n_complete < 30, 'walked' semantics, 2585/2586 | Tie rule now in the docstring: the first path in s2 order, the primary when tied; a tie is flagged `:tie`. `HG008N_hap*_walked` in per_node is graph containment only. 2585/2586 are not resolved from partial paths; they are now off-phase ambiguous. |

### interpretation audit (biology)

| # | sev | finding | decision |
|---|---|---|---|
| B1 | major | other_hap_only is not LOH of the event hap | **Fixed (wording).** I checked the audit table b6_local_loh.tsv for the old 141 other_hap_only loci: E_lost 0, O_lost 45, both_retained 91. New text: the tumor event hap now carries the other hap's allele by local conversion or recurrent slippage; a tumor-only caller sees a germline allele. s9 does not test LOH itself. |
| B2 | major | 'a linear caller with PoN loses them too' holds for SNVs, not INDELs | **Fixed.** New table: % tagged in each category vs the 16,493 non-absorbed truths. Absent SNV: repo rule 100% vs 8.3% non-absorbed; pop AF ≥ 0.05: 74.7% vs 1.3%. Absent INDEL: repo rule 99.0% vs 78.7%; AF ≥ 1e-3: 98.7% vs 74.9%; ≥ 0.01: 93.5% vs 47.0%; ≥ 0.05: 69.9% vs 13.5%. The repo INDEL path has no PoN, so for INDELs all 1,934 absent are graph-specific losses (about 30% at a 0.05 floor). |
| B3 | major | absent 'SNVs' are mostly repeat length changes / fragments | **Fixed (columns + table).** `asm_event_kind` gives the normal-frame INFO event kind and `n_truths_same_asm_event` the number of truths sharing it. Absent SNVs in the normal frame: SNV 82, INS 43, DEL 21, none 4. 34 share their event with another truth. normal_allele_class: germline_STR_other_length 91, SNV_site_deleted_on_event_hap 30, both_REF 21, other 8. The CpG count was not recomputed; that stays in the audit. |
| B4 | minor | 19 'truth conflict' loci are assembly-analysis conflicts | **Fixed.** After the array calls, 2 remain (7473 both, 10740 event_hap_only). Both are other_ambiguous evidence_conflict. |
| B5 | minor | superpopulation shifts | Fixed (P5). |
| B6 | minor | GIAB nogermlineoverlap scope | **Fixed.** Column `in_nogermline_bed` and a fourth scope in every category table. Union: SNV 2/2/64/3, INDEL 43/14/840/6. |
| B7 | minor | 29 absent SNVs without tumor-read support | **Noted, list not adopted.** Of the 29, 22 are absent now and 7 ambiguous. Of the 22: 8 no_truth_spelling_read, 10 not Illumina-perfect (no perfect-platform reads yet), 4 with truth reads. tumor_any_ALT_array is False for all 29. The audit's windows are heuristic; my flags (`truth_read_support`, `tumor_verdict`, `T_info_class`) are the per-locus evidence. |
| B8 | note | normal-assembly STR errors bounded at about 6-8% per hap | **Noted.** No HG008-N reads are on disk; the bound is indirect, from the audit. |
| B9 | note | present is small because of the patient genotype, not the panel | **Noted.** Audit figure: 140/196 germline-identical truths are already absorbed. |
| B10 | note | mechanism wording | **Fixed.** Absent exact_allele now reads: somatic slippage or substitution recreating a population allele kept by the d9 graph. |

## Caveats of the fixed version

- **Read data:** only Illumina has read data, so truth_read_support and path choice use Illumina reads even at loci that are not Illumina-perfect. The flag counts only perfect platforms with data.
- **Read classes:**
  - They spell the graph path of the walk over its anchor span; reads may carry edits more than 5 bp from the site.
  - When the read cut ends inside a long interrupted array, the class is representation-dependent. Example: 14267, where 51 reads are 'other'.
  - A truth overlapping a germline record of the same hap cannot be composed edit-wise, so the read is 'other'.
- **Array calls:**
  - Exact periodic runs miss imperfect VNTR copies.
  - There are 16 haps with no array sequence: 10 narrow unresolved, 4 narrow ALT, 2 narrow no. Those 14 that are not 'no' put their locus in other_ambiguous.
- **Tumor assembly:** contig names are not phased to the normal; a 'novel' copy can be an assembly error, for example in A15 homopolymers.
- **Frequencies:** they are conditional on complete traversals and upper-biased (P4). Exact allele is computed over the whole array, so SNVs in arrays rarely have exact HPRC carriers.
- **2585 / 2586:** no complete HPRC haplotype; not resolved from partial paths.
