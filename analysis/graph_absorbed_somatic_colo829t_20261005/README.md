# COLO829T somatic truths absorbed by the HPRC graph: repeat context, and whose alleles they are (COLO829BL or other HPRC individuals)

2026-10-05. COLO829T version of `analysis/graph_absorbed_somatic_20261001` (HG008T): the sequence context of the lost truths
(sections 1-6) and the HPRC haplotypes, individuals and populations that carry them (sections 7-8). Pansoma makes a
candidate (and a tensor) only from the edits a read has relative to the graph nodes it aligns to. When a somatic allele is
already a path of the HPRC v1.1 d9 graph, its reads follow that branch or skip edge with no edit, and the truth gets no
tensor. This analysis takes every COLO829T somatic truth allele lost this way and asks:

1. how many SNVs and INDELs, per platform (fiberseq, ONT, Illumina);
2. in what sequence context they sit: homopolymer, STR, other repeats or none, in the patient's own germline sequence
   (the matched normal COLO829BL, verkko 2.1 hapX / hapY) and in GRCh38;
3. how the absorption rate depends on that context, the SMaHT region, the truth VAF and the normal's own genotype;
4. how this compares with HG008T;
5. which HPRC v1.1 haplotypes, individuals and populations carry each absorbed allele and its graph nodes, and whether the
   allele is the patient's own germline allele (COLO829BL: germline filtering) or other people's (a pangenome-induced false
   negative).

## Short answer

- **Count.** The ALT reads align perfectly (no edit) to a graph branch or skip edge for:
  - fiberseq: 64 SNVs / 218 INDELs; ONT: 58 / 170; Illumina: 53 / 396.
  - That is 8.9-12.1% of each platform's no-candidate SNVs and 49.0-57.6% of its no-candidate INDELs. Every
    such INDEL is class I1 of the miss analysis.
  - Union of the three platforms: **67 SNVs + 430 INDELs = 497 truth alleles** (22.5% of all 1,912 chr1-22 truth
    INDELs, 0.16% of the 42,035 SNVs); chr1: 5 + 38.
  - A d9 path spells the truth (alone or with the normal's nearby germline alleles) at 64 SNVs / 409 INDELs.
- **The absorbed INDELs are repeat changes.** In the patient frame (union, 430):
  - **325 (75.6%) change the length of a tandem repeat by whole units**: homopolymer `HP>=7` 241 (56.0%), `STR` 83 (19.3%),
    `HP4-6` 0, `VNTR>6` 1. 316 of these 325 (97.2%) change exactly one unit (one base of a homopolymer, one
    copy of an STR unit).
  - 19 more sit in or next to a repeat without being a whole-unit change (`in_STR` 17, inside an imperfect VNTR /
    minisatellite 2), and 16 are substitutions inside a repeat (the GRCh38 INDEL record absorbs a germline length
    difference). Only 9 (2.1%) are non-repeat.
  - The homopolymers are mostly 10-29 bp: 7-9 bp 6, 10-14 bp 86, 15-19 bp 94, 20-29 bp 55, >=30 bp 0; 210 of 241 are A/T. STR motifs: (AC)n 34, (AT)n 23, (AGAT)n 13.
  - The rest: the normal already carries the ALT (33, 7.7%), a complex change (27, 6.3%), unresolved 1
    (no sequence, or only one haplotype).
  - **Germline-like or somatic?** At 35 more absorbed INDELs (`SNV in STR2-6>=3copies` 15, `complex` 9, `SNV` 6, `HP>=7` 2, `none` 1, `STR2-6>=3copies` 1, `SNV in HP4-6` 1) a COLO829BL
    haplotype already has exactly the ALT length but other bases than GRCh38 + truth (mostly a germline SNP in the
    array, or a tie where the other haplotype is GRCh38). Either the tumor
    changed the patient's repeat by one unit, or the truth is a germline length allele written against GRCh38; the tumor
    allele was not assembled, so this is not decided. So up to 68 (15.8%) of the absorbed INDELs may be
    germline-like, not only the 33 (7.7%) where the normal carries the ALT exactly (section 2).
- **The absorption rate rises with homopolymer length.** Over all 1,912 truth INDELs, a one-base change of a homopolymer is
  absorbed at 5.2% (7-9 bp), 20.3% (10-14 bp), 25.1% (15-19 bp), 32.5% (20-29 bp); a one-unit STR change at
  46.1%; a non-repeat INDEL at 1.3% (3 / 235; INDELs inside an imperfect VNTR are not
  counted as non-repeat: 2 / 4 absorbed).
- **The patient's germline matters.**
  - At 33 INDELs a COLO829BL haplotype already has the truth allele over the whole tandem array; 32 of them are
    in SMaHT Extreme regions. There the graph acts as germline filtering, as for HG008T's 105 `HG008N_present` INDELs.
    33 (60.0%) of the 55 truth INDELs whose allele the normal carries are absorbed.
  - A germline repeat-length allele in the normal (another array length on >= 1 haplotype) raises the rate for STRs
    (52.2% vs 38.0% when both haplotypes equal GRCh38), hardly for homopolymers (24.1% vs 21.3%).
- **SNVs are graph SNP alleles, not repeat events.** S3 (the ALT is an existing graph SNP allele) is 57 of 64 (fiberseq), 51 of 58 (ONT), 50 of 53 (Illumina).
  38 of the 67 are non-repeat in both frames; 28 are CpG transitions, which are absorbed at
  1.04% against 0.10% for the other truth SNVs.
- **Region, not VAF.** 33 / 562 (5.9%) of the INDELs in SMaHT Easy regions are absorbed, against 186 / 641 (29.0%) (Difficult) and
  211 / 709 (29.8%) (Extreme). Easy regions hold the short homopolymers, but the gap remains at the same length
  (10-14 bp: 13.0% Easy vs 27.1% Difficult). The truth VAF hardly matters: 22.5% (0.1-0.25), 25.3% (0.25-0.4), 21.1% (>=0.4).
- **Whose allele is it? Mostly other people's** (sections 7-8; HPRC v1.1 walks of the d9 graph, complete traversals only).
  - **Pangenome-induced false negative: 427 (85.9%)** (62 SNVs, 365 INDELs). Neither COLO829BL haplotype carries the
    ALT, but HPRC haplotypes walk the same nodes. At 388 of them (58 SNVs / 330 INDELs) at least one HPRC haplotype
    carries exactly the allele: median exact-allele frequency 0.204 / 0.170 (SNV / INDEL), carried by a median of
    15 / 12 HPRC individuals, in >= 2 superpopulations at 333 loci and by a single individual at 11.
    The INDELs are the repeat slippages of section 2 (239 `HP>=7`, 75 `STR`); the SNVs are common graph SNP alleles.
  - **Germline filtering: 33 (6.6%)**, the 'normal carries ALT' loci of section 2 (24 common in HPRC, 9 rare).
  - **Ambiguous: 37 (7.4%)**, mainly the 24 'closest' loci (no d9 path spells the truth).
  - **Populations.** The carriers of the absent alleles lean AFR, the panel's largest and most diverse group: observed / expected
    carrier haplotypes AFR 1.18 (z +17.7), AMR 0.82 (z -13.5), EAS 0.79 (z -6.6), SAS 0.86 (z -2.2). HPRC v1.1 has no
    EUR haplotype; the COLO829 donor is a European (white) male.
  - **Same as HG008T**: 84.9% vs 85.0% of the absorbed INDELs are pangenome-induced false negatives,
    7.4% vs 4.6% germline filtering; each HPRC haplotype carries the exact allele of a median
    19.9% of the absent INDELs it traverses (HG008T 20.0%).
- **Compared with HG008T: the same mechanism, a different truth set** (both sides classed with the HG008 b5 rule here).
  - Per-length rates are close up to 19 bp: one-base homopolymer change 4.7 vs 4.6% (7-9 bp), 20.0 vs 21.5% (10-14 bp), 25.6 vs 28.8% (15-19 bp), 31.0 vs 38.4% (20-29 bp)
    (COLO829T vs HG008T); COLO829T is lower at 20-29 bp. One-unit STR change 45.8% vs 37.8%.
  - Homopolymers are a smaller share of the absorbed INDELs (58.1% vs 81.4%) mainly because the COLO829T truth has fewer
    and shorter homopolymer INDELs: 60.7% of its truth INDELs are `HP>=7` events (HG008T 79.1%), and only 1 is a
    one-base change of a homopolymer >= 30 bp (HG008T 513). Overall INDEL absorption: 22.5% vs 26.8%.
  - Absorbed SNVs: 0.16% of truth SNVs vs 2.3%. HG008T's were mostly homopolymer / STR events written as GRCh38 SNVs
    (180 of 201 in a GRCh38 `HP>=7` or `STR`); COLO829T's are mostly non-repeat graph SNP alleles (18 of 67 in such a repeat).

## Definitions

| Term | Meaning |
|---|---|
| d9 graph | HPRC v1.1 Minigraph-Cactus GRCh38 graph, frequency-filtered (nodes on fewer than ~9 of the 90 haplotypes removed). The COLO829T tensors were built on it |
| perfect bypass | status `no_candidate` on the platform (tensor sets of 2026-09-30) and the majority of the ALT-like reads leave GRCh38 through a non-GRCh38 node run or an edge that skips GRCh38 bases, with no edit within +-5 bp: read-level reason `site_bypassed:branch:no_edit` / `site_bypassed:skip_edge:no_edit` of the COLO829T miss analysis. **Absorbed** = perfect bypass on >= 1 platform (the union) |
| perfect set | the truths that are perfect bypass on one platform; per-set tables count truth alleles |
| d9 match | (c2) a d9 path between the window anchors spells GRCh38 + truth: `exact`; or only with nearby COLO829BL dipcall PASS alleles applied too: `with_germline`; no such path, the nearest graph haplotypes: `closest` |
| element | one departure of the ALT path from GRCh38: `B:x:y:>n1>n2` = a run of non-GRCh38 nodes between GRCh38 nodes x and y; `S:x:y` = an edge x->y that skips GRCh38 bases. Subtypes `ins_branch`, `del_skip`, `snv_branch`, `mnv_branch`, `replacing_branch` |
| GRCh38 frame | the truth VCF record in GRCh38, classed with the HG008 audit b1 rule (c1): unit = minimal period of the inserted / deleted bases, tract = the longest exact period-unit stretch touching the event |
| tandem array | (c3, the HG008 s8b rule) the event window grown by every exact periodic stretch of period 1-60 touching it (cap 3 kb), bracketed by GRCh38 24-mers that are unique and found once in each placed COLO829BL haplotype |
| hap copy | (c3) where a COLO829BL haplotype's sequence is cut from: the best hit of a 10-kb GRCh38 window on the same-chromosome ragtag scaffold (minimap2 asm5); where that copy is suspect (MAPQ < 20, a covering hit on another scaffold, a dipcall phased het SNV contradicted, or no hit) the contig copy dipcall itself aligned there (`dipcall copy`, the paralog fix). Anchors are searched within 400 bp of the array, else within 2 kb (`wide anchor`) so that both haplotypes share them |
| array allele | a COLO829BL haplotype's sequence between the anchors (assembly; the dipcall reconstruction where the assembly cannot be anchored), with germline differences outside the array reverted: `REF` (= GRCh38), `ALT` (= GRCh38 + truth), `germline_len` (another length), `germline_seq` (same length, other bases), `NA`. `flank_len` = the reverted length (>= 4 bp flagged; the hap's total change is `total_len_change`) |
| event hap / patient frame | the haplotype whose array allele is closest to GRCh38 + truth (Levenshtein; tie -> the one equal to GRCh38, else hapX). The **patient-frame event** is the edit from that allele to GRCh38 + truth, left-normalised. COLO829T has no truth INFO event (HG008T used the GIAB `HG008Nv62SOMATICVARIANT`), so the event is derived from the normal |
| `HP>=7` / `HP4-6` | INDEL event with unit 1 in a homopolymer of >= 7 / 4-6 bp (the normal's length, the stretch holding the event's own bases) |
| `STR2-6>=3copies` (`STR`) | unit 2-6 with >= 3 copies |
| `VNTR>6` | unit > 6, the inserted / deleted sequence repeated next to the event |
| `in_STR` | none of these, but the event sits inside or touches a period 1-6 stretch >= 10 bp without being a whole-unit change of it: e.g. 2719, a 2-bp TC insertion inside a (GT)n; most absorbed ones are 1-bp changes of another base next to a homopolymer (the b5 rule calls these `HP>=7`) |
| `imperfect VNTR` | a `none` event inside an imperfect tandem repeat (VNTR / minisatellite with unit 7-150 bp and indels between copies) that reaches >= 20 bp beyond the exact-period array: 12-mers recurring 7-150 bp apart cover >= 60 % of every 31-bp window of a >= 40-bp run holding the truth (GRCh38 +-600 bp; `SNV in imperfect VNTR` for substitutions) |
| `none` | INDEL event not in a repeat |
| `SNV in <class>` / `SNV` | the patient-frame event is a substitution: inside an `HP>=7` / `HP4-6` / `STR` stretch, or not. On a GRCh38 INDEL record this means the record absorbs a germline length difference |
| `complex` | the event changes REF and ALT bases at once (no single INDEL or SNV turns the closest haplotype into the ALT) |
| `normal carries ALT` | a COLO829BL haplotype already has GRCh38 + truth over the whole array |
| `NA` / `single hap` (unresolved) | no haplotype sequence (no common anchor, no assembly hit) / only one haplotype has a sequence and it is not the ALT: the missing one may carry the ALT, so no class is given (the one-hap class is in `pf_label_one_hap`) |
| ALT-length hap | (`alt_len_haps`) an INDEL truth where a COLO829BL haplotype's array allele has exactly the ALT length but is not the ALT (other bases, mostly a germline SNP): the normal may already carry the ALT length |
| group | `HP>=7` / `STR` / other repeat (`HP4-6`, `VNTR>6`, `in_STR`, `imperfect VNTR`) / SNV in repeat / non-repeat (`none`, `SNV`) / normal carries ALT / complex / unresolved (`NA`, `single hap`). Repeat-unit changes = `HP>=7`, `STR`, `HP4-6`, `VNTR>6` events; in or next to a repeat = `in_STR`, `imperfect VNTR` |
| b5 rule | the HG008 audit b5 rule verbatim (`pf_class_b5`, `label_b5`, `pf_tract_b5`): the tract is any period-unit stretch merely touching the event; no imperfect-VNTR check. The default (c3) rule requires the stretch to hold the event's own bases. HG008T numbers use the b5 rule, so every COLO829T vs HG008T rate comparison uses the b5 rule on both sides |
| units changed | abs(length change) / unit, for `HP>=7`, `HP4-6`, `STR`, `VNTR>6` events |
| mechanism class | (background rates) `HP>=7 1-unit, tract <bin>` = one-base change of a homopolymer of that length; `STR 1-unit`; `repeat multi-unit` = `HP>=7` / `STR` with >= 2 units; the other labels as above. COLO829T: c3 rule (`mechanism_class`) and b5 rule (`mechanism_class_b5`, label and tract of the b5 rule); HG008T: the same function on its b5 table |
| germline status | from the two array alleles, first that applies: normal carries ALT; germline other length (`germline_len` on >= 1 hap); germline same length, other bases; both REF; one hap NA, other REF; both NA |
| RGN | SMaHT region of the truth (truth VCF): Easy / Difficult / Extreme |
| VAF bin | the truth's Illumina VAF (`VAF_Ill`): < 0.1, 0.1-0.25, 0.25-0.4, >= 0.4 |
| miss classes | S3 = the ALT is an existing graph SNP allele; S4b = repeat, absorbed by a graph path; S4a = repeat, the SNV becomes an edit on another branch; S2 = ALT not seen in same-length reads; I1 = the INDEL allele is fully in the graph |
| somatic elements | (c7) `exact`: every element of the chosen ALT path (the c2 primary; no read data to choose); `with_germline`: minus the elements that alone spell the path's germline alleles; `closest`: the c2 closest primary's elements minus germline-alone ones (the locus is ambiguous) |
| complete traversal | an HPRC haplotype walk (c5, from the d9 GFA) that visits both window anchors in order. Only these count in HPRC frequencies: 88 HPRC haplotypes (44 samples x 2); CHM13 apart; GRCh38 never |
| exact allele | (c7) a haplotype's sequence over the tandem array, between the nearest GRCh38 nodes it shares with the ALT path (or over the whole window), equals GRCh38 + truth; for `with_germline` also the ALT path's allele (GRCh38 + truth + the patient's germline alleles) |
| element set | a walk contains every element of an ALT path (node level; any candidate path) |
| category | (c7, the HG008 s9 rules, section 8) `normal_present_broad` / `_rare` (germline filtering), `normal_absent_HPRC_other` (pangenome-induced false negative; sub_class `exact_allele` / `element_set_other_allele` / `recombinant_pieces`), `other_ambiguous` |
| RARE_AF | 0.20 on the exact-allele frequency among complete traversals (broad >= 0.20 > rare) |
| O/E | observed / expected carrier haplotypes of a superpopulation: expected = carriers x its share of the locus' complete haplotypes, summed over the loci |

## Data

| Item | Source |
|---|---|
| graph | `/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.*` (d9 sqlite / reference-path index via `indel_graph_paths_20260930/classify.py`) |
| truth | `/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz`; chr1-22 alleles in `analysis/tensor_recall_20260930/per_truth_COLO829T.tsv` (all in the BED, all PASS; VAF_Ill, VAF_PB, RGN) |
| no-candidate statuses | `per_truth_COLO829T.tsv` (tensor sets of 2026-09-30, now in `pansoma_v2_tensors/backup_ch6_linear100_20261001`) |
| read-level reasons, miss classes | `analysis/tensor_recall_20260930/colo829t_miss/{fiberseq,ONT,Illumina}/{snv,indel}_no_candidate.tsv` |
| COLO829BL | SMaHT verkko 2.1 donor-specific assembly, ragtag GRCh38 scaffolds hapX / hapY (`/scratch/qfu/COLO829BL_DSA/ragtag_hg38/`, read only); dipcall vs GRCh38 `dipcall_hg38/dipcall_hg38.dip.vcf.gz` + `.dip.bed`, GT = hapY\|hapX (c3 check on 6,232 distinct isolated phased het SNVs near the loci: 6,068 consistent, 0 swapped, 0 other, 164 undetermined, 0 with different verdicts at two loci); dipcall's own contig alignments `dipcall_hg38.hap{1,2}.paf.gz` (hap1 = hapY) and the raw verkko contigs for the dipcall copy |
| GRCh38 | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` |
| HG008T (comparison) | `analysis/graph_absorbed_somatic_20261001/{per_variant,indel_repeat_context,per_locus_populations,hprc_individuals}.tsv`, `per_truth_HG008T.tsv`, audit tables `b1_truth_context.tsv` / `b5_normal_frame.tsv` |
| HPRC walks | the d9 GFA `hprc-v1.1-mc-grch38.d9.gfa` (P and W lines; c5) |
| HPRC genotypes | the d9 VCF and the full-graph `vg deconstruct` VCF (`hprc-v1.1-mc-grch38.raw.vcf.gz`, read through the HG008 `full_vcf_index`), 44 sample columns + CHM13 (c6) |
| populations | HPRC v1.1 sample -> 1000 Genomes population / superpopulation: a byte-identical copy of the HG008 analysis' `hprc_sample_metadata.tsv` (HPRC Year 1 metadata, checked with IGSR / Coriell). AFR 23 / AMR 16 / EAS 4 / SAS 1 / **EUR 0** samples. The COLO829 donor is a 45-year-old white male (ATCC) |

- The 2026-10-01 tensor rebuild (channel 6, build AF 0.08) cannot change this set: a read with no edit gives no candidate
  in any build.

## 1. How many truths are lost because the reads align perfectly to a graph branch

| platform | SNV no candidate | SNV perfect bypass | of which a d9 path spells the truth | INDEL no candidate | INDEL perfect bypass | of which a d9 path spells the truth | chr1 SNV / INDEL perfect bypass |
|---|---:|---:|---:|---:|---:|---:|---:|
| fiberseq | 528 | 64 (12.1%) | 61 | 436 | 218 (50.0%) | 205 | 5 / 23 |
| ONT | 547 | 58 (10.6%) | 56 | 347 | 170 (49.0%) | 154 | 4 / 19 |
| Illumina | 598 | 53 (8.9%) | 51 | 687 | 396 (57.6%) | 380 | 4 / 34 |
| union |  | 67 | 64 |  | 430 | 409 | 5 / 38 |

All chr1-22 truth alleles: SNV 42,035, INDEL 1,912. Perfect bypass on all three platforms: SNV 46, INDEL 134. Share of all truth: fiberseq SNV 0.15% / INDEL 11.4%; ONT SNV 0.14% / INDEL 8.9%; Illumina SNV 0.13% / INDEL 20.7%; union SNV 0.16% / INDEL 22.5%

- **Platforms.** 134 INDELs and 46 SNVs are perfect bypass on all three platforms.
- **Miss classes and d9 match** (`c2_graph_paths.py`):

| platform | SNV S3 | SNV S4b | SNV S4a | SNV S2 | INDEL I1 | SNV exact | SNV with_germline | SNV closest | INDEL exact | INDEL with_germline | INDEL closest |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fiberseq | 57 | 6 | 0 | 1 | 218 | 61 | 0 | 3 | 195 | 10 | 13 |
| ONT | 51 | 2 | 4 | 1 | 170 | 56 | 0 | 2 | 144 | 10 | 16 |
| Illumina | 50 | 2 | 1 | 0 | 396 | 51 | 0 | 2 | 367 | 13 | 16 |
| union |  |  |  |  |  | 64 | 0 | 3 | 394 | 15 | 21 |

S3 = the ALT is an existing graph SNP allele; S4b = repeat, absorbed by a graph path; S4a = repeat, the SNV becomes an edit on another branch; S2 = ALT not seen in same-length reads; I1 = the INDEL allele is fully in the graph. exact / with_germline / closest: c2 (a d9 path spells GRCh38 + truth; only with COLO829BL dipcall alleles applied too; no such path).

- At the 21 INDEL and 3 SNV 'closest' loci no d9 path spells the truth, so the reads that align perfectly there follow another allele (most likely a germline repeat length; the reads were not re-decoded for COLO829T).

## 2. Repeat context of the absorbed INDELs

Each INDEL truth is classed in the **patient frame**: the change from the closest COLO829BL haplotype to GRCh38 + truth over the tandem array (`c3_normal_frame.py`). This is the frame that matters: the GRCh38 record of a repeat INDEL depends on the GRCh38 repeat length, while the tumor changed the patient's own repeat.

### Patient frame, grouped (per perfect set)

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| HP>=7 | 79 (36.2%) | 51 (30.0%) | 235 (59.3%) | 241 (56.0%) |
| STR | 57 (26.1%) | 44 (25.9%) | 75 (18.9%) | 83 (19.3%) |
| other repeat | 15 (6.9%) | 15 (8.8%) | 17 (4.3%) | 20 (4.7%) |
| SNV in repeat | 11 (5.0%) | 12 (7.1%) | 12 (3.0%) | 16 (3.7%) |
| non-repeat | 8 (3.7%) | 8 (4.7%) | 8 (2.0%) | 9 (2.1%) |
| normal carries ALT | 27 (12.4%) | 22 (12.9%) | 27 (6.8%) | 33 (7.7%) |
| complex | 21 (9.6%) | 18 (10.6%) | 21 (5.3%) | 27 (6.3%) |
| unresolved | 0 (0.0%) | 0 (0.0%) | 1 (0.3%) | 1 (0.2%) |
| total | 218 | 170 | 396 | 430 |

- **Per platform.** Illumina's perfect set is 59.3% homopolymer; fiberseq and ONT are 36.2% / 30.0% homopolymer and have relatively more STR, normal-carries-ALT and complex loci (the same direction as HG008T: PacBio 73.6%, ONT 66.9%, Illumina 84.2% homopolymer).
- **Detailed labels** (`repeat_context.md` 2a): `in_STR` 17, `HP4-6` 0, `VNTR>6` 1; substitution events `SNV in STR2-6>=3copies` 15, `SNV in HP4-6` 1; non-repeat `none` 3, `SNV` 6.
- **Size of the change.** 316 of the 325 repeat INDEL events (97.2%) change one unit, 7 change 2-3 units, 2 more.
- **Homopolymer length** (the normal's): 7-9 bp 6, 10-14 bp 86, 15-19 bp 94, 20-29 bp 55, >=30 bp 0. Bases A/T 210, C/G 31; insertions 112, deletions 129.
- **STR motifs** (canonical unit, both strands): AC 34, AT 23, AGAT 13, AG 5, AAAG 3, AAAT 2, AAT 2, AAATG 1.
- **Patient frame vs GRCh38 record.** Of the 347 patient-frame INDEL events, 278 have the GRCh38 record's length change, 52 another size and 17 the opposite sign (e.g. GRCh38 says a deletion, but the patient's haplotypes are shorter, so the tumor inserted). The event haplotype's array allele is GRCh38 at 270 loci, another germline length at 115, same length other bases at 11 and ALT at 33.
- **Normal carries ALT** (33): the carrying haplotype is hapX 30, hapY 3; 32 are in SMaHT Extreme regions; VAF_Ill >= 0.25 at 32; assembly and dipcall both give ALT on a haplotype at 23; 10 carry a c3 flag (`nf_flags`), 6 a duplicated-region sign (MAPQ < 20 or an other-scaffold hit). The other haplotype's change to the ALT: none 13, HP>=7 8, in_STR 3, HP4-6 3, complex 3, STR2-6>=3copies 2, SNV in HP>=7 1; at the 13 'none' loci the normal is heterozygous for the truth allele, a non-repeat INDEL (one haplotype ALT, the other one non-repeat INDEL away; the other haplotype is GRCh38 at 25 of all 33).
- **Complex** (27; 24 in Extreme regions): the other haplotype gives complex 20, in_STR 6, VNTR>6 1; d9 match exact 14, with_germline 8, closest 5.
- **The odd labels come with a germline difference in the array.** 62 of the 63 other-repeat, SNV-in-repeat and complex INDELs have a germline allele other than GRCh38 in the array. A = GRCh38 + truth keeps GRCh38's bases there, so the change from the patient's haplotype to A also reverts germline differences. Examples: 949 (`complex`; a homozygous germline C>G two bases from a (TG)n insertion, the c2 `with_germline` allele; hapX already has the +14 length), 2719 (`in_STR`; GRCh38 + GTGT keeps GRCh38's C where the patient has G, so the event becomes a 2-bp TC insertion inside (GT)n), 1996 (`SNV in STR`; hapX is (CT)13(AT)8, the length of GRCh38 + truth).
- **Two readings.** At 35 absorbed INDELs a COLO829BL haplotype has exactly the ALT length but other bases than A (`SNV in STR2-6>=3copies` 15, `complex` 9, `SNV` 6, `HP>=7` 2, `none` 1, `STR2-6>=3copies` 1, `SNV in HP4-6` 1; column `alt_len_haps`). (1) The tumor changed the patient's repeat by one unit and the truth record was written against GRCh38; or (2) the truth is a germline length allele the normal already carries, like the 'normal carries ALT' group. No tumor assembly was used, so the two are not told apart: the germline-like share of the absorbed INDELs is between 7.7% (33) and 15.8% (68).

### GRCh38 frame (b1 rule, per perfect set)

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| HP>=7 | 81 (37.2%) | 47 (27.6%) | 238 (60.1%) | 244 (56.7%) |
| STR2-6>=3copies | 86 (39.4%) | 71 (41.8%) | 99 (25.0%) | 120 (27.9%) |
| HP4-6 | 2 (0.9%) | 3 (1.8%) | 3 (0.8%) | 3 (0.7%) |
| VNTR>6 | 2 (0.9%) | 2 (1.2%) | 2 (0.5%) | 2 (0.5%) |
| none | 47 (21.6%) | 47 (27.6%) | 54 (13.6%) | 61 (14.2%) |
| total | 218 | 170 | 396 | 430 |

- 61 of the union (14.2%) look non-repeat in GRCh38. In the patient frame 13 of them are whole-unit repeat changes, 6 sit in or next to a repeat (`in_STR`, imperfect VNTR), 9 are substitutions in a repeat, 15 normal carries ALT, 12 complex, 0 unresolved and 6 non-repeat (normal carries ALT 15, complex 12, SNV in STR2-6>=3copies 8, STR2-6>=3copies 7, HP>=7 6, SNV 4, in_STR 4, none 2, imperfect VNTR 2, SNV in HP4-6 1). The patient's germline repeat length or an interruption hides the repeat in GRCh38. Union crosstab (`repeat_context.md` 2g):

| GRCh38 class \ patient frame | HP>=7 | STR | other repeat | SNV in repeat | non-repeat | normal carries ALT | complex | unresolved | total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 227 | 0 | 5 | 0 | 0 | 9 | 3 | 0 | 244 |
| STR2-6>=3copies | 8 | 76 | 8 | 7 | 3 | 6 | 11 | 1 | 120 |
| HP4-6 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 3 |
| VNTR>6 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 2 |
| none | 6 | 7 | 6 | 9 | 6 | 15 | 12 | 0 | 61 |

GRCh38-frame "none" INDELs (union) in the patient frame: normal carries ALT 15, complex 12, SNV in STR2-6>=3copies 8, STR2-6>=3copies 7, HP>=7 6, SNV 4, in_STR 4, none 2, imperfect VNTR 2, SNV in HP4-6 1. The same INDELs under the b5 rule in GRCh38 (c1 b5rule_class): in_STR 32, none 29.

## 3. Background: absorption rate by context, all chr1-22 truth INDELs

Patient frame (c3 covers all 1,912 truth INDELs):

| class | truth INDELs | absorbed | % |
|---|---:|---:|---:|
| HP>=7 1-unit, tract 7-9 | 116 | 6 | 5.2% |
| HP>=7 1-unit, tract 10-14 | 424 | 86 | 20.3% |
| HP>=7 1-unit, tract 15-19 | 375 | 94 | 25.1% |
| HP>=7 1-unit, tract 20-29 | 163 | 53 | 32.5% |
| HP>=7 1-unit, tract >=30 | 1 | 0 | 0.0% |
| STR 1-unit | 165 | 76 | 46.1% |
| repeat multi-unit | 48 | 9 | 18.8% |
| HP4-6 | 67 | 0 | 0.0% |
| VNTR>6 | 13 | 1 | 7.7% |
| in_STR | 92 | 17 | 18.5% |
| imperfect VNTR | 4 | 2 | 50.0% |
| SNV in repeat | 65 | 16 | 24.6% |
| non-repeat INDEL | 235 | 3 | 1.3% |
| non-repeat SNV | 11 | 6 | 54.5% |
| complex | 72 | 27 | 37.5% |
| normal carries ALT | 55 | 33 | 60.0% |
| unresolved | 6 | 1 | 16.7% |

- Same by label (`repeat_context.md` 3a): `HP>=7` 241 / 1,103 (21.8%), `STR` 83 / 189 (43.9%), `none` 3 / 235 (1.3%), `normal carries ALT` 33 / 55 (60.0%).

GRCh38 frame:

| class | truth INDELs | absorbed | % |
|---|---:|---:|---:|
| HP>=7 | 1,167 | 244 | 20.91% |
| STR2-6>=3copies | 275 | 120 | 43.64% |
| HP4-6 | 97 | 3 | 3.09% |
| VNTR>6 | 4 | 2 | 50.00% |
| none | 369 | 61 | 16.53% |
| total | 1,912 | 430 | 22.49% |

## 4. SNVs

- 67 absorbed SNVs: S3 57 (fiberseq), 51 (ONT), 50 (Illumina); S4b 6, 2, 2; d9 match exact 64, closest 3.
- GRCh38 context: `none` 42, `STR` 14, `HP4-6` 7, `HP>=7` 4. Patient frame: `SNV` (non-repeat) 38, substitution in a repeat 23, `complex` 3, `in_STR` INDEL event 2, normal carries ALT 1.
- CpG transitions: 28 of 67; 22 of the 38 non-repeat ones.

All truth SNVs, GRCh38 frame:

| class | truth SNVs | absorbed | % |
|---|---:|---:|---:|
| HP>=7 | 545 | 4 | 0.73% |
| STR2-6>=3copies | 1,322 | 14 | 1.06% |
| HP4-6 | 5,914 | 7 | 0.12% |
| none | 34,254 | 42 | 0.12% |
| total | 42,035 | 67 | 0.16% |

| class | truth SNVs | absorbed | % |
|---|---:|---:|---:|
| CpG Ti | 2,705 | 28 | 1.04% |
| not CpG Ti | 39,330 | 39 | 0.10% |
| C>A | 5,178 | 6 | 0.12% |
| C>G | 1,338 | 0 | 0.00% |
| C>T | 28,417 | 40 | 0.14% |
| T>A | 2,136 | 9 | 0.42% |
| T>C | 2,762 | 7 | 0.25% |
| T>G | 2,204 | 5 | 0.23% |

- The truth SNVs are UV-dominated (C>T 28,417 of 42,035) and almost never absorbed (0.14% of C>T). The absorbed ones are enriched for CpG transitions, the commonest population SNPs, which the d9 graph keeps as SNP bubbles (S3).

Absorbed SNVs, GRCh38 class x CpG transition x patient frame:

| GRCh38 class | CpG Ti | patient frame | SNVs |
|---|---:|---:|---:|
| none | CpG Ti | SNV | 22 |
| none | other | SNV | 16 |
| STR2-6>=3copies | other | SNV in STR2-6>=3copies | 9 |
| HP4-6 | CpG Ti | SNV in HP4-6 | 4 |
| HP>=7 | other | SNV in HP>=7 | 3 |
| none | other | SNV in imperfect VNTR | 3 |
| HP4-6 | other | SNV in HP4-6 | 3 |
| STR2-6>=3copies | other | complex | 3 |
| STR2-6>=3copies | CpG Ti | in_STR | 1 |
| HP>=7 | other | in_STR | 1 |
| STR2-6>=3copies | CpG Ti | SNV in STR2-6>=3copies | 1 |
| none | other | normal carries ALT | 1 |

## 5. Region, VAF and the normal's genotype

Cells are absorbed (union) / all chr1-22 truth INDELs of the patient-frame group and stratum.

### SMaHT region

| patient-frame group \ RGN | Easy | Difficult | Extreme |
|---|---:|---:|---:|
| HP>=7 | 24 / 247 (9.7%) | 135 / 498 (27.1%) | 82 / 358 (22.9%) |
| STR | 5 / 23 (21.7%) | 35 / 68 (51.5%) | 43 / 98 (43.9%) |
| other repeat | 1 / 78 (1.3%) | 8 / 29 (27.6%) | 11 / 69 (15.9%) |
| SNV in repeat | 0 / 11 (0.0%) | 3 / 11 (27.3%) | 13 / 43 (30.2%) |
| non-repeat | 2 / 198 (1.0%) | 1 / 17 (5.9%) | 6 / 31 (19.4%) |
| normal carries ALT | 0 / 0 | 1 / 6 (16.7%) | 32 / 49 (65.3%) |
| complex | 1 / 5 (20.0%) | 2 / 11 (18.2%) | 24 / 56 (42.9%) |
| unresolved | 0 / 0 | 1 / 1 (100.0%) | 0 / 5 (0.0%) |
| all INDELs | 33 / 562 (5.9%) | 186 / 641 (29.0%) | 211 / 709 (29.8%) |

Cells: absorbed (union) / all chr1-22 truth INDELs of the class and stratum (% absorbed).

Homopolymer one-base changes by length and region:

| tract \ RGN | Easy | Difficult | Extreme |
|---|---:|---:|---:|
| tract 7-9 | 3 / 80 (3.8%) | 0 / 16 (0.0%) | 3 / 20 (15.0%) |
| tract 10-14 | 20 / 154 (13.0%) | 45 / 166 (27.1%) | 21 / 104 (20.2%) |
| tract 15-19 | 1 / 11 (9.1%) | 64 / 252 (25.4%) | 29 / 112 (25.9%) |
| tract 20-29 | 0 / 0 | 24 / 58 (41.4%) | 29 / 105 (27.6%) |
| tract >=30 | 0 / 0 | 0 / 0 | 0 / 1 (0.0%) |

Cells: absorbed (union) / all chr1-22 truth INDELs whose patient-frame event is a 1-base change of a homopolymer of that length.

- SNVs: Easy 33 / 32,945 (0.1%); Difficult 11 / 4,922 (0.2%); Extreme 23 / 4,168 (0.6%).

### Truth VAF (Illumina)

| patient-frame group \ VAF | <0.1 | 0.1-0.25 | 0.25-0.4 | >=0.4 |
|---|---:|---:|---:|---:|
| HP>=7 | 0 / 0 | 19 / 105 (18.1%) | 84 / 355 (23.7%) | 138 / 643 (21.5%) |
| STR | 0 / 0 | 14 / 25 (56.0%) | 33 / 63 (52.4%) | 36 / 101 (35.6%) |
| other repeat | 0 / 0 | 2 / 18 (11.1%) | 1 / 44 (2.3%) | 17 / 114 (14.9%) |
| SNV in repeat | 0 / 0 | 1 / 3 (33.3%) | 4 / 17 (23.5%) | 11 / 45 (24.4%) |
| non-repeat | 0 / 0 | 1 / 15 (6.7%) | 2 / 63 (3.2%) | 6 / 168 (3.6%) |
| normal carries ALT | 0 / 0 | 1 / 5 (20.0%) | 10 / 18 (55.6%) | 22 / 32 (68.8%) |
| complex | 0 / 0 | 3 / 11 (27.3%) | 14 / 27 (51.9%) | 10 / 34 (29.4%) |
| unresolved | 0 / 0 | 0 / 0 | 1 / 3 (33.3%) | 0 / 3 (0.0%) |
| all INDELs | 0 / 0 | 41 / 182 (22.5%) | 149 / 590 (25.3%) | 240 / 1,140 (21.1%) |

Cells: absorbed (union) / all chr1-22 truth INDELs of the class and stratum (% absorbed).

### Normal-frame germline status

| patient-frame group \ germline status | normal carries ALT | germline other length | germline same length, other bases | both REF | one hap NA, other REF | both NA |
|---|---:|---:|---:|---:|---:|---:|
| HP>=7 | 0 / 0 | 56 / 232 (24.1%) | 6 / 30 (20.0%) | 179 / 841 (21.3%) | 0 / 0 | 0 / 0 |
| STR | 0 / 0 | 47 / 90 (52.2%) | 1 / 7 (14.3%) | 35 / 92 (38.0%) | 0 / 0 | 0 / 0 |
| other repeat | 0 / 0 | 12 / 51 (23.5%) | 7 / 22 (31.8%) | 1 / 103 (1.0%) | 0 / 0 | 0 / 0 |
| SNV in repeat | 0 / 0 | 16 / 65 (24.6%) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| non-repeat | 0 / 0 | 8 / 25 (32.0%) | 1 / 6 (16.7%) | 0 / 215 (0.0%) | 0 / 0 | 0 / 0 |
| normal carries ALT | 33 / 55 (60.0%) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| complex | 0 / 0 | 27 / 66 (40.9%) | 0 / 6 (0.0%) | 0 / 0 | 0 / 0 | 0 / 0 |
| unresolved | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 3 (0.0%) | 1 / 3 (33.3%) |
| all INDELs | 33 / 55 (60.0%) | 166 / 529 (31.4%) | 15 / 71 (21.1%) | 215 / 1,251 (17.2%) | 0 / 3 (0.0%) | 1 / 3 (33.3%) |

Cells: absorbed (union) / all chr1-22 truth INDELs of the class and stratum (% absorbed).

- Overall, `germline other length` loci are absorbed more often (31.4%) than `both REF` loci (17.2%), but much of that is class mix: `both REF` holds 215 of the 246 non-repeat INDELs. Within a class the germline length allele matters for STRs (52.2% vs 38.0%) and other repeats (23.5% vs 1.0%), hardly for homopolymers (24.1% vs 21.3%). None of the truth INDELs labelled 'SNV in repeat' or 'complex' has both haplotypes equal to GRCh38 in the array (0 / 0 such loci).

## 6. Side by side with HG008T

HG008T numbers are recomputed by `c4_tables.py` from the HG008 files (checked against `indel_repeat_context.md`). HG008T patient frame = the GIAB INFO event in the HG008-N v6.2 event haplotype, classed with the b5 rule. For a same-rule comparison every rate and share below uses the b5 rule for COLO829T too (the tables show the c3 rule beside it). The events themselves still differ in origin: COLO829T derives them from the closest normal haplotype, HG008T takes the INFO event.

|  | COLO829T | HG008T |
|---|---:|---:|
| truth alleles SNV / INDEL | 42,035 / 1,912 | 8,690 / 8,496 |
| fiberseq / PacBio HiFi: SNV / INDEL perfect (% of no candidate) | 64 (12.1%) / 218 (50.0%) | 170 (42.7%) / 1,363 (45.8%) |
| ONT: SNV / INDEL perfect | 58 (10.6%) / 170 (49.0%) | 148 (43.3%) / 989 (40.3%) |
| Illumina: SNV / INDEL perfect | 53 (8.9%) / 396 (57.6%) | 122 (36.0%) / 1,916 (48.5%) |
| union SNV (% of all truth SNVs) | 67 (0.2%) | 201 (2.3%) |
| union INDEL (% of all truth INDELs) | 430 (22.5%) | 2,277 (26.8%) |

| group | COLO829T (c3 rule) | COLO829T (b5 rule) | HG008T (b5 rule) |
|---|---:|---:|---:|
| HP>=7 | 241 (56.0%) | 250 (58.1%) | 1,853 (81.4%) |
| STR | 83 (19.3%) | 89 (20.7%) | 292 (12.8%) |
| other repeat | 20 (4.7%) | 3 (0.7%) | 6 (0.3%) |
| SNV in repeat | 16 (3.7%) | 16 (3.7%) | 104 (4.6%) |
| non-repeat | 9 (2.1%) | 11 (2.6%) | 1 (0.0%) |
| normal carries ALT | 33 (7.7%) | 33 (7.7%) | 0 (0.0%) |
| complex | 27 (6.3%) | 27 (6.3%) | 0 (0.0%) |
| unresolved | 1 (0.2%) | 1 (0.2%) | 21 (0.9%) |
| total | 430 | 430 | 2,277 |

HG008T patient frame = the GIAB truth INFO event (HG008Nv62SOMATICVARIANT) in the HG008-N v6.2 event haplotype, always a simple event, so HG008T has no 'complex' or 'normal carries ALT' rows; 'unresolved' = no INFO event there. HG008T truths whose tumor allele is the normal's other-haplotype allele (category HG008N_present_*): 105 INDELs (4.6%), classed by their INFO event inside the rows above.

|  | COLO829T (c3 rule): of all absorbed INDELs | COLO829T (c3): of repeat / HP>=7 events | COLO829T (b5 rule): of all | COLO829T (b5): of repeat / HP>=7 events | HG008T (b5): of all absorbed INDELs | HG008T (b5): of repeat / HP>=7 events |
|---|---:|---:|---:|---:|---:|---:|
| 1 unit | 316 (73.5%) | 316 (97.2%) | 331 (77.0%) | 331 (97.4%) | 2,103 (92.4%) | 2,103 (97.9%) |
| 2-3 units | 7 (1.6%) | 7 (2.2%) | 7 (1.6%) | 7 (2.1%) | 38 (1.7%) | 38 (1.8%) |
| >3 units | 2 (0.5%) | 2 (0.6%) | 2 (0.5%) | 2 (0.6%) | 8 (0.4%) | 8 (0.4%) |
| HP 7-9 bp | 6 (1.4%) | 6 (2.5%) | 6 (1.4%) | 6 (2.4%) | 21 (0.9%) | 21 (1.1%) |
| HP 10-14 bp | 86 (20.0%) | 86 (35.7%) | 89 (20.7%) | 89 (35.6%) | 377 (16.6%) | 377 (20.3%) |
| HP 15-19 bp | 94 (21.9%) | 94 (39.0%) | 99 (23.0%) | 99 (39.6%) | 519 (22.8%) | 519 (28.0%) |
| HP 20-29 bp | 55 (12.8%) | 55 (22.8%) | 56 (13.0%) | 56 (22.4%) | 806 (35.4%) | 806 (43.5%) |
| HP >=30 bp | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 130 (5.7%) | 130 (7.0%) |

b5 rule = the HG008 rule on both sides (the same-rule comparison); the c3 rule is the COLO829T default of sections 2-5.

### Background rates

| class | COLO829T | HG008T |
|---|---:|---:|
| HP>=7 | 244 / 1,167 (20.9%) | 1,688 / 6,360 (26.5%) |
| STR2-6>=3copies | 120 / 275 (43.6%) | 369 / 936 (39.4%) |
| HP4-6 | 3 / 97 (3.1%) | 15 / 241 (6.2%) |
| VNTR>6 | 2 / 4 (50.0%) | 2 / 6 (33.3%) |
| complex | 0 / 0 | 1 / 1 (100.0%) |
| none | 61 / 369 (16.5%) | 202 / 952 (21.2%) |
| total | 430 / 1,912 (22.5%) | 2,277 / 8,496 (26.8%) |

| class | COLO829T (c3 rule) | COLO829T (b5 rule) | HG008T (b5 rule) |
|---|---:|---:|---:|
| HP>=7 1-unit, tract 7-9 | 6 / 116 (5.2%) | 6 / 129 (4.7%) | 21 / 453 (4.6%) |
| HP>=7 1-unit, tract 10-14 | 86 / 424 (20.3%) | 89 / 444 (20.0%) | 376 / 1,748 (21.5%) |
| HP>=7 1-unit, tract 15-19 | 94 / 375 (25.1%) | 99 / 387 (25.6%) | 517 / 1,795 (28.8%) |
| HP>=7 1-unit, tract 20-29 | 53 / 163 (32.5%) | 54 / 174 (31.0%) | 795 / 2,072 (38.4%) |
| HP>=7 1-unit, tract >=30 | 0 / 1 (0.0%) | 0 / 1 (0.0%) | 122 / 513 (23.8%) |
| STR 1-unit | 76 / 165 (46.1%) | 82 / 179 (45.8%) | 270 / 715 (37.8%) |
| repeat multi-unit | 9 / 48 (18.8%) | 9 / 49 (18.4%) | 44 / 221 (19.9%) |
| HP4-6 | 0 / 67 (0.0%) | 0 / 87 (0.0%) | 1 / 204 (0.5%) |
| VNTR>6 | 1 / 13 (7.7%) | 1 / 13 (7.7%) | 3 / 15 (20.0%) |
| in_STR | 17 / 92 (18.5%) | 2 / 20 (10.0%) | 2 / 27 (7.4%) |
| imperfect VNTR | 2 / 4 (50.0%) | 0 / 0 | 0 / 0 |
| SNV in repeat | 16 / 65 (24.6%) | 16 / 65 (24.6%) | 104 / 253 (41.1%) |
| non-repeat INDEL | 3 / 235 (1.3%) | 5 / 220 (2.3%) | 0 / 415 (0.0%) |
| non-repeat SNV | 6 / 11 (54.5%) | 6 / 11 (54.5%) | 1 / 24 (4.2%) |
| complex | 27 / 72 (37.5%) | 27 / 72 (37.5%) | 0 / 0 |
| normal carries ALT | 33 / 55 (60.0%) | 33 / 55 (60.0%) | 0 / 0 |
| unresolved | 1 / 6 (16.7%) | 1 / 6 (16.7%) | 21 / 41 (51.2%) |

Compare the two b5 columns (the same rule, tract = the b5 tract). HG008T row counts differ from mechanism.md, which counts SNV and INDEL records together; here INDEL records only.

| class | COLO829T | HG008T |
|---|---:|---:|
| HP>=7 | 4 / 545 (0.7%) | 73 / 399 (18.3%) |
| STR2-6>=3copies | 14 / 1,322 (1.1%) | 107 / 572 (18.7%) |
| HP4-6 | 7 / 5,914 (0.1%) | 4 / 880 (0.5%) |
| none | 42 / 34,254 (0.1%) | 17 / 6,839 (0.2%) |
| CpG Ti | 28 / 2,705 (1.0%) | 11 / 782 (1.4%) |
| total | 67 / 42,035 (0.2%) | 201 / 8,690 (2.3%) |

- **Shares (b5 rule on both sides).** Homopolymers are 58.1% of COLO829T's absorbed INDELs vs 81.4% of HG008T's; STRs 20.7% vs 12.8%. One-unit changes are 97.4% vs 97.9% of the repeat events. COLO829T has no absorbed homopolymer >= 30 bp (7.0% of HG008T's `HP>=7` events) and fewer at 20-29 bp (22.4% vs 43.5%).
- **Rates (b5 rule on both sides).** COLO829T minus HG008T, percentage points: one-base homopolymer changes +0.0 (7-9 bp), -1.5 (10-14 bp), -3.2 (15-19 bp), -7.3 (20-29 bp), one-unit STR changes +8.0. Up to 19 bp the rates are close; at 20-29 bp COLO829T is lower (54 / 174), and STRs are absorbed more often. The lower homopolymer share follows the truth set: `HP>=7` events are 60.7% of COLO829T truth INDELs and 79.1% of HG008T's (c3 rule for COLO829T: 57.7%).
- **SNVs.** The no-candidate SNV counts are similar (528-598 vs 339-398 per platform), but a much smaller share is perfect bypass (8.9-12.1% vs 36.0-43.3%): HG008T's GRCh38 SNVs in homopolymers / STRs are absorbed at 18.3% / 18.7%, COLO829T's at 0.7% / 1.1%.
- **Rules.** With the b5 rule COLO829T has `HP>=7` 250 and `STR` 89 (18 absorbed INDELs change label, 14 of them c3 `in_STR` -> b5 `HP>=7` / `STR`).
- **Germline.** COLO829T 'normal carries ALT' 33 (7.7%) corresponds to HG008T's 105 `HG008N_present_*` INDELs (4.6%), which HG008T's table classes by their INFO event (inside the `HP>=7` / `STR` rows). HG008T has no 'complex' row because its INFO events are simple.

## 7. HPRC haplotypes, individuals, populations

`c5_hprc_walks.py` (the HG008 s4 code) scans the d9 GFA once (P and W lines, 24 processes) and extracts every haplotype's walk between the window anchors. A median of 79 of the 88 HPRC haplotypes per locus have a complete traversal (SNV 85, INDEL 78); 37 loci have fewer than 44, and truth 28014 has none (a 'closest' locus). `c6_hprc_vcf.py` reads the HPRC genotypes of the d9 and the full-graph VCF. `c7_hprc_membership.py` (the HPRC part of HG008 s9) counts the carriers at three levels:

| level | carrier = | normal_absent_HPRC_other SNV / INDEL | normal_present (SNV + INDEL) |
|---|---:|---:|---:|
| **exact allele** (the categories use this) | the walk spells GRCh38 + truth over the tandem array (with_germline: + the germline alleles) | 0.189 / 0.159 (0.204 / 0.170 at loci with >= 1 carrier) | 0.318 |
| element set | the walk contains all elements of any ALT path | 0.205 / 0.191 | 0.368 |
| element (min) | the walk contains each element (minimum over the elements) | 0.205 / 0.198 | 0.364 |

- **Absent loci by HPRC support** (`sub_class`; individuals and superpopulations from `per_locus_populations.tsv`):

| sub_class | kind | loci | median exact-allele freq | median element-set freq | median carrier individuals | >= 2 superpops | one superpop | one population | one individual | none | CHM13 exact |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `exact_allele` | SNV | 58 | 0.204 | 0.206 | 15 | 46 | 10 | 0 | 2 | 0 | 9 |
| `exact_allele` | INDEL | 330 | 0.170 | 0.193 | 12 | 287 | 32 | 2 | 9 | 0 | 65 |
| `element_set_other_allele` | SNV | 4 | 0.000 | 0.158 | 9.5 | 3 | 1 | 0 | 0 | 0 | 0 |
| `element_set_other_allele` | INDEL | 33 | 0.000 | 0.180 | 11 | 29 | 1 | 1 | 2 | 0 | 1 |
| `recombinant_pieces` | INDEL | 2 | 0.000 | 0.000 | 0 | 0 | 0 | 0 | 0 | 2 | 0 |

Carrier individuals / superpopulations (per_locus_populations.tsv): exact-allele carriers; for `element_set_other_allele` / `recombinant_pieces` the haplotypes walking all elements of an ALT path. First match of one individual > one population > one superpopulation.

  - Unlike HG008T (85 of 150 absent SNVs without an exact HPRC carrier), 58 of the 62 absent COLO829T SNVs have one: they are mostly non-repeat graph SNP alleles (section 4), so the population allele is the same allele.
  - `element_set_other_allele`: HPRC haplotypes walk all the somatic nodes but spell another allele over the array, mostly another repeat length.
- **'Rare'** means an exact-allele frequency below 0.20 among complete traversals. 9 present loci are rare; 1 of them (551) has no exact HPRC carrier although HPRC haplotypes walk its nodes. Rarity in large populations is gnomAD / CoLoRSdb (section 8, PoN table).
- **d9 floor.** The d9 graph keeps nodes walked by >= ~9 of the 90 haplotypes. Of the 413 somatic branch elements, 12 have a branch node that fewer than 9 haplotypes visit in the c5 window rows: 10 of them lie on CHM13 (a reference path, which the filter keeps); the other 2 (3346: 8, 37542: 4 haplotypes) are window counts, a lower bound (no whole-GFA scan was run for them).
- **Frequencies are conditional** on complete traversals and on the d9 graph: a haplotype cut inside the window does not count. `hprc_exact_allele_freq_over88` is the lower bound.
- **Cross-checks.**
  - d9 / full-graph VCF (c6 `carriers_any` vs the GFA allele-level carriers, on GFA-complete haplotypes): identical carrier sets at 53 / 65 SNV and 308 / 427 INDEL loci (d9), 53 / 65 and 312 / 427 (full). The VCF counts snarl alleles: most differences are haplotypes with the record's allele but another allele elsewhere in the array, or loci without a matching record (`with_germline`, compound repeats; `hprc_tables.md`).
  - Sequence-level re-derivation (`c7_hprc_membership.py check`: the nearest 16-mers left and right of the allele window that are unique in GRCh38, in GRCh38 + truth and in each walk; carrier = the walk's sequence between them equals GRCh38 + truth): identical carrier sets at 473 of 474 anchored loci (34,053 locus-haplotype pairs compared; 1,467 skipped because a germline variant sits in an anchor; 17 loci have no unique anchor in the window). The exception is 29562 (10 by c7, 14 by sequence, of 43): these haplotypes spell GRCh38 + truth between the anchors but compensate beyond the array bracket, so c7 undercounts there (the category does not change).
  - Window check: none of the 17,311 complete walks that spell GRCh38 over the whole window contains a somatic element or counts as a carrier; all 6,297 walks that spell GRCh38 + truth over the whole window are exact carriers (1 only through the whole-window rule, at 29562).
- **Populations.** The alleles that hide somatic events (absent, exact allele) lean AFR: observed / expected carrier haplotypes AFR 1.18 (z +17.7), AMR 0.82 (z -13.5), EAS 0.79 (z -6.6), SAS 0.86 (z -2.2). AFR is the panel's largest and most diverse group. The patient's own germline alleles (present_broad, 24 loci) lean slightly AFR too (AFR 1.07 (z +2.7), AMR 0.92 (z -2.3)), unlike HG008T's (AFR 0.92, AMR 1.10). HPRC v1.1 has no EUR haplotype, so the donor's ancestry is not in the panel. CHM13 has the exact allele at 75 absent and 6 present loci.

| category | level | loci | carrier haplotypes | AFR O/E (z) | AMR O/E (z) | EAS O/E (z) | SAS O/E (z) | CHM13 carries |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `normal_present_broad` | exact allele | 24 | 747 | 1.07 (z +2.7) | 0.92 (z -2.3) | 0.96 (z -0.5) | 0.94 (z -0.3) | 3 |
| `normal_present_broad` | element set | 24 | 756 | 1.06 (z +2.4) | 0.92 (z -2.2) | 0.97 (z -0.3) | 0.92 (z -0.4) | 3 |
| `normal_present_rare` | exact allele | 9 | 73 | 1.02 (z +0.1) | 0.91 (z -0.6) | 1.24 (z +0.7) | 1.08 (z +0.1) | 3 |
| `normal_present_rare` | element set | 9 | 101 | 1.01 (z +0.1) | 0.93 (z -0.7) | 1.19 (z +0.8) | 1.16 (z +0.3) | 3 |
| `normal_absent_HPRC_other` | exact allele | 427 | 6594 | 1.18 (z +17.7) | 0.82 (z -13.5) | 0.79 (z -6.6) | 0.86 (z -2.2) | 75 |
| `normal_absent_HPRC_other` | element set | 427 | 7876 | 1.17 (z +19.0) | 0.83 (z -14.6) | 0.80 (z -7.0) | 0.88 (z -2.1) | 96 |
| `other_ambiguous` | exact allele | 31 | 44 | 1.07 (z +0.5) | 0.94 (z -0.4) | 1.14 (z +0.3) | 0.66 (z -0.5) | 1 |
| `other_ambiguous` | element set | 31 | 358 | 0.91 (z -2.3) | 1.11 (z +2.2) | 0.97 (z -0.3) | 1.19 (z +0.8) | 7 |
| all | exact allele | 491 | 7458 | 1.16 (z +17.6) | 0.84 (z -13.6) | 0.82 (z -6.3) | 0.86 (z -2.2) | 82 |
| all | element set | 491 | 9091 | 1.15 (z +18.1) | 0.85 (z -14.0) | 0.83 (z -6.7) | 0.90 (z -1.9) | 109 |

Expected per locus = carriers x the superpopulation's share of that locus' complete haplotypes; z hypergeometric, ignoring the pairing of haplotypes within individuals and linkage between loci (it overstates significance). HPRC v1.1: AFR 23, AMR 16, EAS 4, SAS 1 samples, no EUR; the COLO829 donor is a European (white) male.

  Per 1000 Genomes population (`hprc_populations.tsv`; summed over the loci):

| superpop | population | name | HPRC haplotypes | absent INDEL | absent SNV | present INDEL (broad + rare) |
|---|---:|---:|---:|---:|---:|---:|
| AFR | ACB | African Caribbean in Barbados | 14 | 0.237 | 0.277 | 0.350 |
| AFR | ASW | African Ancestry in Southwest US | 2 | 0.230 | 0.291 | 0.412 |
| AFR | ESN | Esan in Nigeria | 2 | 0.237 | 0.298 | 0.500 |
| AFR | GWD | Gambian in Western Division, The Gambia | 16 | 0.230 | 0.274 | 0.392 |
| AFR | MKK | Maasai in Kinyawa, Kenya (HapMap 3) | 2 | 0.227 | 0.270 | 0.396 |
| AFR | MSL | Mende in Sierra Leone | 8 | 0.245 | 0.313 | 0.394 |
| AFR | YRI | Yoruba in Ibadan, Nigeria | 2 | 0.215 | 0.231 | 0.378 |
| AMR | CLM | Colombian in Medellin, Colombia | 8 | 0.162 | 0.207 | 0.274 |
| AMR | PEL | Peruvian in Lima, Peru | 8 | 0.160 | 0.186 | 0.296 |
| AMR | PUR | Puerto Rican in Puerto Rico | 16 | 0.167 | 0.183 | 0.351 |
| EAS | CHS | Southern Han Chinese, China | 6 | 0.163 | 0.197 | 0.327 |
| EAS | KHV | Kinh in Ho Chi Minh City, Vietnam | 2 | 0.138 | 0.171 | 0.382 |
| SAS | PJL | Punjabi in Lahore,Pakistan | 2 | 0.173 | 0.174 | 0.327 |
| CHM13 | CHM13 | reference (T2T CHM13, hydatidiform mole) | 1 | 0.183 | 0.148 | 0.188 |

  The absent-INDEL exact-allele fraction is 0.215-0.245 in the AFR populations, 0.160-0.167 AMR, 0.138-0.163 EAS, 0.173 SAS (`populations.md`: loci per superpopulation, patterns, shared / private, one-superpopulation loci by population).
- **Per individual** (`hprc_individuals.tsv`). Each HPRC haplotype carries the exact allele of a median 61 absent INDELs (19.9% of those it traverses completely; 16.7% of all 365) and 14 absent SNVs. By superpopulation: AFR 69 (23.1%), AMR 53 (16.3%), EAS 51 (15.5%), SAS 56.5 (17.3%). The range is 38 (HG00673#1, CHS) to 89 (HG03453#2, MSL); CHM13 66. `hprc_carriers.tsv.gz` lists the carrier haplotypes of every element and allele.

## 8. Cross-evaluation with the normal: germline filtering or pangenome-induced false negative

Categories (`c7_hprc_membership.py`: the HG008 s9 rules, with the COLO829BL whole-array alleles of c3 in place of HG008-N), RARE_AF = 0.20 on the exact-allele HPRC frequency. Rules apply in this order:

| order | category | rule | SNV | INDEL |
|---|---:|---:|---:|---:|
| 1 | `other_ambiguous` | no d9 path spells the truth (c2 `closest`; no read data to take the elements from) | 3 | 21 |
| 1 | `other_ambiguous` | no somatic element (every path element alone spells a germline allele, or none) | 0 | 0 |
| 1 | `other_ambiguous` | the `with_germline` path needs a germline allele that dipcall phases to the other haplotype than the c3 event hap | 0 | 2 |
| 1 | `other_ambiguous` | no HPRC haplotype traverses the window completely | 0 | 0 |
| 1 | `other_ambiguous` | both COLO829BL haplotypes carry the ALT (homozygous germline; truth conflict) | 0 | 0 |
| 2 | `normal_present_broad` | a COLO829BL haplotype carries the ALT over the whole tandem array ('normal carries ALT'); exact-allele HPRC frequency >= 0.20 | 0 | 24 |
| 2 | `normal_present_rare` | the same, frequency < 0.20 | 1 | 8 |
| 3 | `other_ambiguous` | a COLO829BL haplotype has no array sequence and the other one does not carry the ALT | 0 | 1 |
| 3 | `other_ambiguous` | no complete HPRC haplotype walks every somatic element (CHM13 may) | 1 | 9 |
| 4 | `normal_absent_HPRC_other` | both haplotypes resolved, neither carries the ALT; sub_class `exact_allele`: >= 1 complete HPRC haplotype spells the exact allele | 58 | 330 |
| 4 | `normal_absent_HPRC_other` | `element_set_other_allele`: HPRC haplotypes walk all elements of an ALT path but spell another allele over the array | 4 | 33 |
| 4 | `normal_absent_HPRC_other` | `recombinant_pieces`: every element is walked by some HPRC haplotype, a whole ALT path by none | 0 | 2 |

| category | interpretation |
|---|---|
| `normal_present_*` | the patient's germline allele (a COLO829BL haplotype carries GRCh38 + truth over the whole tandem array), already a graph path. A tumor-only caller sees a germline allele there and would filter it anyway: the graph acts as **germline filtering** |
| `normal_absent_HPRC_other` | **pangenome-induced false negative**: the somatic change recreates an allele of other individuals (`exact_allele`), or is spelled by their nodes (`element_set_other_allele`, `recombinant_pieces`) |
| `other_ambiguous` | not resolved; mostly no d9 path spells the truth ('closest'; the COLO829T reads were not re-decoded) |

Counts per category (present_broad / present_rare / absent_HPRC_other / other_ambiguous):

| set | SNV chr1-22 | INDEL chr1-22 | SNV chr1 | INDEL chr1 |
|---|---:|---:|---:|---:|
| fiberseq | 0 / 1 / 59 / 4 (64) | 18 / 8 / 175 / 17 (218) | 0 / 0 / 3 / 2 (5) | 7 / 3 / 12 / 1 (23) |
| ONT | 0 / 1 / 54 / 3 (58) | 17 / 4 / 129 / 20 (170) | 0 / 0 / 3 / 1 (4) | 7 / 1 / 9 / 2 (19) |
| Illumina | 0 / 1 / 50 / 2 (53) | 20 / 6 / 345 / 25 (396) | 0 / 0 / 3 / 1 (4) | 9 / 2 / 22 / 1 (34) |
| **union** | **0 / 1 / 62 / 4 (67)** | **24 / 8 / 365 / 33 (430)** | 0 / 0 / 3 / 2 (5) | 10 / 3 / 22 / 3 (38) |

Cells: normal_present_broad / normal_present_rare / normal_absent_HPRC_other / other_ambiguous (total). chr1 = the chr1 truths (every COLO829T truth is inside the BED).

- **Sub-classes per platform** (chr1-22 and chr1): `hprc_tables.md`.
- **Present loci** (33): the carrying haplotype is hapX 31, hapY 2; 32 are in SMaHT Extreme regions. They inherit the caveats of section 2 ('Normal carries ALT': dipcall copy, duplicated-region signs). The one other 'normal carries ALT' locus, 38917, is a 'closest' locus (rule 1).
- **Germline-like reading** (section 2, 'Two readings'). At 22 absent loci (`exact_allele` 12, `element_set_other_allele` 9, `recombinant_pieces` 1) a COLO829BL haplotype already has the ALT length with other bases (`germline_like_alt_len`); 13 more such loci are ambiguous. If these are germline length alleles, the pangenome-induced INDEL losses are between 343 and 365.
- **Patient-frame group per category:**

| patient-frame group | SNV normal_present_broad | SNV normal_present_rare | SNV normal_absent_HPRC_other | SNV other_ambiguous | INDEL normal_present_broad | INDEL normal_present_rare | INDEL normal_absent_HPRC_other | INDEL other_ambiguous |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 0 | 0 | 0 | 0 | 0 | 0 | 239 | 2 |
| STR | 0 | 0 | 0 | 0 | 0 | 0 | 75 | 8 |
| other repeat | 0 | 0 | 1 | 1 | 0 | 0 | 17 | 3 |
| SNV in repeat | 0 | 0 | 20 | 3 | 0 | 0 | 10 | 6 |
| non-repeat | 0 | 0 | 38 | 0 | 0 | 0 | 7 | 2 |
| normal carries ALT | 0 | 1 | 0 | 0 | 24 | 8 | 0 | 1 |
| complex | 0 | 0 | 3 | 0 | 0 | 0 | 17 | 10 |
| unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

### Is the graph helping or hurting?

- **Germline filtering: 33 (6.6%).** The tumor allele is a COLO829BL germline allele. COLO829T has no truth INFO event, so whether the tumor's other haplotype changed to it (HG008T's case: local conversion or recurrent slippage) or the truth is a germline allele is not decided here; either way a tumor-only caller sees a germline allele, so the graph costs nothing a tumor-only design could keep.
- **Interference: 427 (85.9%).**
  - INDELs (365), patient frame: `HP>=7` 239, `STR` 75, `other repeat` 17, `complex` 17, `SNV in repeat` 10, `non-repeat` 7. Mostly a one-unit slippage that recreates a length allele found in HPRC, which the d9 graph keeps as a branch or skip edge (section 3: the absorption rate rises with homopolymer length).
  - SNVs (62): `non-repeat` 38, `SNV in repeat` 20, `complex` 3, `other repeat` 1. 58 have an exact HPRC carrier (median frequency 0.204): common population SNP alleles (S3, CpG transitions enriched, section 4).
- **Ambiguous: 37 (7.4%)**: 24 closest (6 with a graph allele nearer the truth than GRCh38), 10 with no complete HPRC haplotype walking the elements (2 CHM13 only), 2 off-phase germline, 1 normal unresolved.
- **What a linear-reference tumor-only caller would lose anyway (population PoN).** SNV: the repo PoN rule tags 100.0% (3.3%) of the absent SNVs (in brackets: the other chr1-22 truth SNVs), pop-AF >= 0.05 83.9% (0.1%). A linear caller would lose most of these SNVs too. INDEL: the repo INDEL path has no PoN; a population-AF filter would tag the absent INDELs 98.6% / 95.1% / 69.3% at AF 1e-3 / 0.01 / 0.05, against 58.4% / 34.2% / 9.0% of the other truth INDELs.

| group | kind | n | repo rule % | pop AF >= 1e-3 % | pop AF >= 0.01 % | pop AF >= 0.05 % |
|---|---:|---:|---:|---:|---:|---:|
| `normal_present_broad` | SNV | 0 |  |  |  |  |
| `normal_present_rare` | SNV | 1 | 100.0 | 100.0 | 100.0 | 100.0 |
| `normal_absent_HPRC_other` | SNV | 62 | 100.0 | 100.0 | 100.0 | 83.9 |
| `other_ambiguous` | SNV | 4 | 75.0 | 75.0 | 75.0 | 25.0 |
| all absorbed | SNV | 67 | 98.5 | 98.5 | 98.5 | 80.6 |
| not absorbed (other COLO829T truths) | SNV | 41,968 | 3.3 | 1.2 | 0.3 | 0.1 |
| `normal_present_broad` | INDEL | 24 | 100.0 | 100.0 | 100.0 | 100.0 |
| `normal_present_rare` | INDEL | 8 | 100.0 | 100.0 | 100.0 | 75.0 |
| `normal_absent_HPRC_other` | INDEL | 365 | 99.5 | 98.6 | 95.1 | 69.3 |
| `other_ambiguous` | INDEL | 33 | 81.8 | 78.8 | 54.5 | 12.1 |
| all absorbed | INDEL | 430 | 98.1 | 97.2 | 92.3 | 66.7 |
| not absorbed (other COLO829T truths) | INDEL | 1,482 | 63.0 | 58.4 | 34.2 | 9.0 |

repo rule = `scripts/filter_panel_of_normals.py` (allele match; gnomAD / CoLoRSdb AF >= 1e-4; dbSNP non-somatic; 1000G), from `analysis/tensor_recall_20260930/pon/truth_pon_COLO829T.tsv`; the repo INDEL path has no PoN. pop AF = max(gnomAD, CoLoRSdb).

### HG008T vs COLO829T

|  | COLO829T | HG008T |
|---|---:|---:|
| absorbed truth alleles, SNV / INDEL | 67 / 430 | 201 / 2,277 |
| germline filtering (the normal carries the ALT: present broad + rare), SNV / INDEL | 1 (1.5%) / 32 (7.4%) | 5 (2.5%) / 105 (4.6%) |
| pangenome-induced false negative (absent from the normal, HPRC carries), SNV / INDEL | 62 (92.5%) / 365 (84.9%) | 150 (74.6%) / 1,935 (85.0%) |
|   of which >= 1 HPRC haplotype has the exact allele, SNV / INDEL | 58 / 330 | 65 / 1,789 |
| ambiguous, SNV / INDEL | 4 (6.0%) / 33 (7.7%) | 46 (22.9%) / 237 (10.4%) |
| median complete HPRC traversals per locus (of 88) | 79 | 76 |
| absent: median exact-allele frequency, SNV / INDEL | 0.189 / 0.159 | 0.000 / 0.167 |
| absent: the same at loci with >= 1 carrier, SNV / INDEL | 0.204 / 0.170 | 0.105 / 0.177 |
| present: median exact-allele frequency | 0.318 | 0.315 |
| absent exact-allele loci carried in >= 2 superpopulations / by one HPRC individual only | 333 (85.8%) / 11 (2.8%) | 1,672 (90.2%) / 65 (3.5%) |
| absent, exact allele: O/E (z) AFR / AMR / EAS / SAS | 1.18 (+17.7) / 0.82 (-13.5) / 0.79 (-6.6) / 0.86 (-2.2) | 1.11 (+21.4) / 0.89 (-16.7) / 0.89 (-7.4) / 0.92 (-2.5) |
| present broad, exact allele: O/E (z) AFR / AMR | 1.07 (+2.7) / 0.92 (-2.3) | 0.92 (-5.9) / 1.10 (+5.3) |
| CHM13 has the exact allele (absent loci) | 75 | 298 |
| per HPRC haplotype: median share of the absent INDELs it traverses whose exact allele it carries (median count) | 0.199 (61 of 365) | 0.200 (315 of 1,935) |

HG008T from analysis/graph_absorbed_somatic_20261001 (s9 / s13; recomputed here with the same functions and checked against its README / tables.md). The same category rules: HG008T `HG008N_*` = COLO829T `normal_*`; closest loci are ambiguous on both sides; HG008T chose among several ALT paths with re-decoded reads, COLO829T takes the c2 primary. % = of the absorbed alleles of that kind.

- **The same picture for INDELs.** 84.9% vs 85.0% of the absorbed INDELs are pangenome-induced false negatives; the exact allele is about as common in HPRC and as widely shared, and each HPRC haplotype carries it at about 20% of the absent INDELs it traverses on both sides.
- **More germline filtering in COLO829T** (7.4% vs 4.6% of the INDELs): the 'normal carries ALT' loci, 32 of 33 in SMaHT Extreme regions.
- **SNVs differ.** COLO829T's absent SNVs are graph SNP alleles that HPRC haplotypes carry exactly (58 of 62); HG008T's were repeat events, 85 of 150 without an exact carrier. COLO829T has fewer ambiguous SNVs (6.0% vs 22.9%; HG008T had 39 closest SNVs).
- **Populations.** The absent alleles lean AFR on both sides (AFR O/E 1.18 vs 1.11). The patient's germline alleles (present_broad) lean AMR for HG008T (1.10, z +5.3) but not for COLO829T (0.92, z -2.3; 24 loci). Neither donor's (European) ancestry is in HPRC v1.1.

## Files

This folder (`columns.tsv` explains every column of `per_variant.tsv` and of the HPRC files, and gives its source file and column):

| File | Content |
|---|---|
| `README.md` | this report (written by `c4_tables.py`) |
| `per_variant.tsv` | one row per absorbed truth allele (497 rows, 180 columns): truth and VAF / RGN; per platform perfect flag, status, miss class, read-level reason, read counts; d9 window, match, primary path and its elements (compact), in-event element subtypes; GRCh38-frame class / unit / tract (and the b5 rule), SNV CpG / Ti-Tv / trinucleotide; COLO829BL array, per-hap placement, assembly / dipcall calls, array allele and its distance to ALT; event hap, patient-frame event, unit, tract, class (c3 and b5 rule), group, units changed, homopolymer bin, germline status, mechanism class, the other hap's event, flags; R / A / hapX / hapY array sequences |
| `columns.tsv` | meaning and source of every column of `per_variant.tsv`, `hprc_per_variant.tsv`, `hprc_per_node.tsv`, `hprc_carriers.tsv.gz`, `hprc_individuals.tsv`, `hprc_populations.tsv`, `per_locus_populations.tsv` |
| `repeat_context.md` | every repeat-context table of `c4_tables.py` (2a-2h, 3a-3i, 4a-4d, 5a-5h), printed by the script too (the tables of sections 7-8 are printed only) |
| `hprc_per_variant.tsv` | one row per absorbed truth allele (497 rows, 80 columns): perfect flags, d9 match, chosen ALT path and its somatic elements, allele window; HPRC complete / partial / absent haplotypes, exact-allele / element-set / element frequencies (and over 88), carrier individuals, d9 / full VCF AF and agreement, per-superpopulation complete / carrier haplotypes, CHM13, the carrier haplotypes; COLO829BL hapX / hapY array alleles, event hap, pattern, germline-like flag; category, sub_class, sub_reason, interpretation |
| `hprc_per_node.tsv` | one row per (truth allele, somatic element), 623 rows: node ids and sequences, GRCh38 interval, subtype, `alone`; HPRC carriers / frequency, per superpopulation, node coverage (d9 floor), CHM13, carrier haplotypes |
| `hprc_carriers.tsv.gz` | one row per (truth_id, level, carrying haplotype with a complete traversal; HPRC and CHM13): `element` rows (= `hprc_per_node` carriers) and `allele` rows (`carries_set_any_path`, `carries_exact_allele` = `hprc_per_variant` `hprc_exact_carrier_haps`); denominators = `hprc_n_complete` / `{SP}_complete` |
| `hprc_individuals.tsv` | per HPRC haplotype (and CHM13): loci per category and kind it traverses completely / carries (exact allele, element set) |
| `hprc_populations.tsv` | per 1000 Genomes population: haplotypes, exact-allele carrier haplotypes and complete traversals per category and kind |
| `per_locus_populations.tsv`, `populations.md` | per truth allele: carrier / complete haplotypes per superpopulation and per 1000G population, superpopulation pattern, shared vs private (one superpopulation / population / individual); `populations.md` = the tables per group (+ O/E, per-haplotype summary) |
| `hprc_tables.md` | every table of `c7_hprc_membership.py` (categories per platform and scope, sub-classes, frequency distributions per level, O/E, crosstabs, VCF cross-check, d9 floor; `$D` column names) |
| `audit_decisions.md` | the two audits of the first version (numbers; independent normal-frame check): every finding, the decision and what changed |
| `c0_*.py` ... `c8_populations.py` | the scripts (below) |

Intermediate data, in `$D = /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005/`:
`variant_set.tsv` (c0), `grch38_context.tsv` (c1, all 43,947 truth alleles), `loci.tsv` (c1b), `graph_paths.tsv` and
`graph_elements.tsv` (c2), `normal_frame.tsv` (c3, 1,979 loci: the absorbed truths + all truth INDELs), `c3_dipcall_phase.tsv`
(GT-order check), `c3_queries.fa`, `c3_hapX.paf` / `c3_hapY.paf` / `c3_*_unloc.paf` and their minimap2 logs;
`prefix_20261005/` (c3 / c4 outputs before the audit fixes; c4 reads its `normal_frame.tsv` for the 'label changed' row of 2h);
`audit/` (the audits' own scripts and tables); `hprc_local_paths.tsv.gz`, `hprc_haplotypes.tsv`, `hprc_window_coverage.tsv` (c5);
`hprc_vcf_alleles.tsv`, `hprc_vcf_haplotypes.tsv.gz`, `hprc_sample_metadata.tsv` (c6); `hprc_per_variant.tsv` (all c7 columns),
`hprc_per_node.tsv`, `hprc_membership.tsv.gz` (every haplotype x locus x element / allele, carriers or not), `hprc_tables.md`,
`hprc_summary.json`, `hprc_spot_checks.txt`, `hprc_check_kmer.tsv` (c7); `hprc_columns.tsv` (c8). Job scripts and logs:
`tmp/graph_absorbed_somatic_colo829t_20261005/`.

Scripts (each has a docstring with inputs, method, outputs and assumptions, and its results at the end):

| Step | Script / job | What it does |
|---|---|---|
| c0 | `c0_variant_set.py` | the 497 perfect-bypass truth alleles (`variant_set.tsv`) |
| c1 | `c1_grch38_context.py` | GRCh38-frame repeat context of every chr1-22 truth allele, b1 rule (`grch38_context.tsv`); `--check-hg008` compares the port with the HG008 audit tables |
| c1b | `c1b_loci.py` | d9 graph window, anchors and event window per locus (`loci.tsv`) |
| c2 | `c2_graph_paths.py`, `c2_graph_paths.sbatch` | d9 paths spelling the truth ALT; their elements (`graph_paths.tsv`, `graph_elements.tsv`) |
| c3 | `c3_normal_frame.py`, `c3_minimap2.sbatch` | COLO829BL array alleles (ragtag copy, or dipcall's copy where the ragtag one is suspect), imperfect-VNTR check, the patient-frame event and class (`normal_frame.tsv`, `c3_dipcall_phase.tsv`) |
| c4 | `c4_tables.py` | `per_variant.tsv`, `columns.tsv`, `repeat_context.md`, this README (run last: it reads the c7 / c8 outputs) |
| c5 | `c5_hprc_walks.py`, `c5_hprc_walks.sbatch` | HPRC haplotype walks through every window from the d9 GFA (`hprc_local_paths.tsv.gz`, `hprc_haplotypes.tsv`, `hprc_window_coverage.tsv`; HG008 s4) |
| c6 | `c6_hprc_vcf.py` | d9 / full-graph VCF carriers by sequence (`hprc_vcf_alleles.tsv`, `hprc_vcf_haplotypes.tsv.gz`); population metadata copy (HG008 s5) |
| c7 | `c7_hprc_membership.py` | somatic elements, HPRC membership (element / set / exact allele), COLO829BL calls, categories (`hprc_per_variant.tsv`, `hprc_per_node.tsv`, `hprc_membership.tsv.gz`, `hprc_tables.md`, `hprc_summary.json`); `check` = the sequence-level re-derivation (`hprc_check_kmer.tsv`) (HG008 s9, HPRC part) |
| c8 | `c8_populations.py` | the HPRC files of this folder, `per_locus_populations.tsv`, `populations.md`, `$D/hprc_columns.tsv` (HG008 s10 HPRC part + s13) |

## Reproduce

All inputs are read only (the COLO829BL files belong to another user: never write or index there). Use
`/wanglab/jshen/anaconda3/bin/python` with `PYTHONDONTWRITEBYTECODE=1`; minimap2 2.28 at `/opt/apps/minimap2/2.28/minimap2`.

1. `python c0_variant_set.py` (login, seconds); `python c1_grch38_context.py` (login, 469 s, 0.19 GB);
   `python c1b_loci.py` (login, 3.4 min, 0.18 GB).
2. `sbatch tmp/graph_absorbed_somatic_colo829t_20261005/c2_graph_paths.sbatch` (32 CPUs, 3G; 6 min, MaxRSS 1.6 GB).
3. `python c3_normal_frame.py queries`; `sbatch tmp/graph_absorbed_somatic_colo829t_20261005/c3_minimap2.sbatch` (16 CPUs,
   24G; 5 min 40 s, peak RSS 12.7 GB by /usr/bin/time, sacct MaxRSS 7.4 GB); the unlocalized check on the login node:
   `minimap2 -x asm5 -c --secondary=yes -N 10 -t 4 <ragtag ..._hap{X,Y}_unlocalized_normalized.fa> $D/c3_queries.fa >
   $D/c3_hap{X,Y}_unloc.paf` (16 s, 1.3 GB); then `python c3_normal_frame.py` (login, 3 min 54 s, 0.26 GB; it also reads dipcall's
   `hap{1,2}.paf.gz` and the raw verkko contigs, read only).
4. `sbatch tmp/graph_absorbed_somatic_colo829t_20261005/c5_hprc_walks.sbatch` (24 CPUs, 16G; 8 min 6 s, MaxRSS 10.1 GB).
5. `python c6_hprc_vcf.py` (login, 2.8 min, 2.9 GB).
6. `python c7_hprc_membership.py 6` (login, 1.5 min, 0.16 GB per process; it imports `c2_graph_paths.py` for the d9 node index);
   `python c7_hprc_membership.py check 6` (login, 43 s).
7. `python c8_populations.py` (login, 2 s, 26 MB).
8. `python c4_tables.py` (login, 11 s, 0.44 GB).

## Assumptions and caveats

- **The patient frame is derived, not given.** HG008T used the GIAB truth INFO event in the HG008-N assembly; COLO829T's
  truth has none, so c3 takes the COLO829BL haplotype closest to GRCh38 + truth and derives the event. Levenshtein
  charges a k-bp indel k, so at some loci the closest haplotype gives `complex` while the other gives a single event
  (flag `other_hap_single_event`, 7 absorbed INDELs; the alternative is in `other_hap_*`). Ties: 15 absorbed truths
  chose the GRCh38 haplotype, 6 fell back to hapX.
- **A = GRCh38 + truth.** The patient's germline differences inside the array are not applied to A (HG008T's INFO event
  was already in the patient's frame), so loci with such differences become `complex`, `in_STR` or `SNV in <class>`
  (section 2, 62 of 63). The `with_germline` d9 match (15 INDELs) is the graph-side view of the same effect.
  At 35 of them a normal haplotype already has the ALT length: these may be germline length alleles rather than
  somatic unit changes (section 2, 'Two readings').
- **The tumor is not checked.** Only the normal's assembly is used. Whether the tumor haplotype carries exactly GRCh38 +
  truth over the array (HG008T: `novel_ALT` vs `novel_other`) is not tested here.
- **Paralog fix.** Where the ragtag copy of a haplotype is suspect, c3 uses the contig copy dipcall aligned there
  (35 of the 1,979 loci, 12 absorbed; the label changed at 4 absorbed truths, `repeat_context.md` 2h).
  That dipcall's long colinear alignment holds the patient's own copy is an assumption. After the fix a copy used
  contradicts an isolated phased het SNV of dipcall at 0 loci (before: 10). No segmental-duplication track was used.
- **Normal carries ALT.** 33 absorbed INDELs, 32 in Extreme regions, mostly on hapX
  (30). 6 of them have a duplicated-region sign (a ragtag hit with MAPQ < 20 or a covering hit on
  another scaffold); of the 6 at chr1:146-149 Mb, 2 do. At 6 of the 6 the haplotypes come from dipcall's copy and agree with
  dipcall by construction; a paralogous copy cannot be excluded. dipcall agreement is not independent evidence (same assembly).
- **Class rules.** `pf_class` requires the repeat to hold the event's own bases; `pf_class_b5` (the HG008 rule) does not. They
  differ at 18 absorbed INDELs; section 6 shows both. SNV events use the b5 rule.
- **Germline differences outside the array** are reverted to GRCh38 by one affine-gap alignment before the event is derived;
  with >= 3 such differences (flag `flank_diffs>=3`) the array / flank split is less reliable. In long or imperfect repeats
  the aligner can put repeat-copy indels at the array edge or the anchor, where they are reverted: the haplotype's array
  allele then understates its repeat-length change (e.g. 40275 hapX: array +96, whole region +475). `total_len_change`
  keeps the whole change and flag `flank_len>=4` marks 14 absorbed INDELs (8 on the event haplotype).
  Event labels at the loci read by hand were not affected.
- **Imperfect repeats.** The arrays are exact-period stretches. A k-mer recurrence check (12-mers 7-150 bp apart; thresholds
  chosen by the audit, not tuned) finds the truth inside a longer imperfect VNTR / minisatellite at 9 absorbed INDELs
  (complex 4, imperfect VNTR 2, normal carries ALT 2, STR2-6>=3copies 1); only `none` events are relabelled (`imperfect VNTR`). The b5 rule and
  the HG008T side have no such check.
- **Wide anchors.** 3 absorbed INDELs needed the 2-kb anchor retry (imperfect VNTR 1, HP>=7 1, normal carries ALT 1); their haplotype
  sequences span longer, often repeat-rich stretches.
- **Unresolved.** 3 loci have no haplotype sequence: 1142 (background; hapX no common anchor, hapY no common anchor); 15472 (background; hapX no common anchor, hapY no common anchor); 41860 (absorbed; hapX no assembly hit, hapY no assembly hit).
  3 have one haplotype only (`single hap`): 110 (background; one-hap label STR2-6>=3copies); 29403 (background; one-hap label imperfect VNTR); 35171 (background; one-hap label STR2-6>=3copies).
- **Truths sharing an array.** 7 absorbed truths have another truth allele inside their tandem array (2115 SNV `SNV in HP4-6` with 2116; 3346 SNV `in_STR` with 3345; 6942 SNV `SNV in imperfect VNTR` with 6943; 6943 SNV `SNV in imperfect VNTR` with 6942; 28711 INDEL `in_STR` with 28712; 40275 INDEL `complex` with 40276; 40276 INDEL `complex` with 40275).
  A applies only the truth itself, so the derived event ignores the other allele (40275 + 40276 applied together still
  match neither haplotype).
- **'closest' loci** (21 INDEL, 3 SNV): no d9 path spells the truth; the perfectly aligned reads were not re-decoded.
- **HPRC frequencies** are over complete traversals, so they are conditional on the d9 graph and biased upward (a haplotype cut
  inside the window does not count); `hprc_exact_allele_freq_over88` is the lower bound. 'Rare' (< 0.20) is relative to HPRC
  v1.1, which has no EUR haplotype; the COLO829 donor is European. The O/E z-scores ignore haplotype pairing and linkage, so
  they overstate significance.
- **No read data for the HPRC part.** HG008T used re-decoded reads to choose among several ALT paths; COLO829T keeps the c2
  primary at the 34 multi-path loci (the element-set level uses every path; the exact-allele level compares
  sequence over the array, so it hardly depends on the path). All 24 'closest' loci stay ambiguous.
- **Event hap in the categories.** The `with_germline` phase check uses c3's derived event hap (the closest haplotype to GRCh38 +
  truth), not a truth INFO event, and 'ALT on the event hap' is no conflict here (the carrying haplotype is the closest one by
  construction). 2 loci fail the phase check.
- **Exact-allele rule.** c7 adds one rule to the HG008 one (a walk that spells GRCh38 + truth over the whole window counts); it
  adds 1 haplotype at 29562. The sequence-level check still finds 4 more carriers there (c7 counts 10 of 43, the sequence 14).
- **Present = c3 'normal carries ALT'**, so the present categories inherit c3's caveats (dipcall copy, duplicated-region signs;
  32 of 33 in SMaHT Extreme regions). The 22 absent loci flagged `germline_like_alt_len` may be germline length
  alleles (section 8).
- **d9 floor.** 2 somatic branch elements not on CHM13 have fewer than 9 visiting haplotypes in the c5 window rows (3346: 8, 37542: 4);
  these are lower bounds, not checked with a whole-GFA scan.
- **GRCh38 frame quirks** (b1 kept unchanged so the classes match HG008T): an INDEL with no stretch found has tract 0; for a
  unit > 6 the tract counts copies only to the right; a deletion such as ATATAT>A gets unit 5.
- **Counting.** One GRCh38 record = one truth allele, as in HG008T. The SNV background is GRCh38 frame only (c3 covers truth
  INDELs and the absorbed SNVs). No COLO829T truth INDEL has VAF_Ill < 0.1.
