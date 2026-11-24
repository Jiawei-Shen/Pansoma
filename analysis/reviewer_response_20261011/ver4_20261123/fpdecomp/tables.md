# False-positive decomposition, six-tumor current models (ver4, 2026-11-23)

chr1, benchmark regions, PASS calls, seed 0. GRCh38-only = six-tumor model without population site features as a model input (m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef); HPRC v2.1 = six-tumor current model (SNV snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, INDEL release_s0_<P>_INDEL). Share = % of the GRCh38-graph FP alleles; candidates retained = % of the alleles of that class with a tensor in the v2.1 test set of the same reads. In the HPRC v2.1 graph = carried by >= 46 of 463 HPRC v2.1 haplotypes. A class with fewer than 20 alleles shows "k of n" instead of the retained percentage; "-" = the class is empty; "<0.1%" = above zero and under 0.05 %.

## Response Table 3.8 (= 2.6), SNV, FP before population filtering

| Dataset | GRCh38-graph FP before population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP before population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 28,178 | 81%; 0.2% | 16%; 87% | 3% | 468 |
| HG008T PacBio HiFi | 84,853 | 85%; 0.3% | 13%; 87% | 2% | 481 |
| HG008T ONT | 50,403 | 83%; 0.5% | 16%; 88% | 1% | 273 |
| COLO829T Illumina | 29,530 | 79%; 0.1% | 11%; 89% | 10% | 550 |
| COLO829T Fiber-seq | 38,677 | 71%; 0.6% | 12%; 83% | 17% | 291 |
| COLO829T ONT | 30,128 | 80%; 0.7% | 12%; 90% | 8% | 477 |
| HCC1395T Illumina | 77,905 | 86%; 0.1% | 14%; 87% | 0.7% | 502 |
| HCC1395T PacBio HiFi | 86,347 | 86%; 0.2% | 13%; 86% | 1% | 416 |
| HCC1395T ONT | 115,580 | 87%; 0.7% | 12%; 86% | 0.9% | 780 |

## Response Table 3.9, INDEL, FP before population filtering

| Dataset | GRCh38-graph FP before population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP before population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 16,016 | 86%; 0.2% | 14%; 65% | 0.8% | 109 |
| HG008T PacBio HiFi | 18,590 | 79%; 15% | 17%; 56% | 5% | 115 |
| HG008T ONT | 3,642 | 71%; 24% | 16%; 66% | 13% | 43 |
| COLO829T Illumina | 5,977 | 85%; <0.1% | 12%; 56% | 4% | 49 |
| COLO829T Fiber-seq | 1,764 | 43%; 19% | 21%; 25% | 36% | 27 |
| COLO829T ONT | 180 | 44%; 29% | 7%; 8 of 12 | 49% | 46 |
| HCC1395T Illumina | 6,594 | 84%; <0.1% | 14%; 73% | 2% | 92 |
| HCC1395T PacBio HiFi | 10,207 | 80%; 2.6% | 15%; 62% | 5% | 126 |
| HCC1395T ONT | 3,713 | 79%; 0.8% | 15%; 75% | 6% | 86 |

## Second pair: SNV, FP after population filtering

| Dataset | GRCh38-graph FP after population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP after population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 517 | 0.2%; 0 of 1 | 40%; 94% | 60% | 460 |
| HG008T PacBio HiFi | 680 | 0.3%; 1 of 2 | 69%; 96% | 31% | 479 |
| HG008T ONT | 429 | 0.9%; 0 of 4 | 79%; 97% | 20% | 273 |
| COLO829T Illumina | 1,480 | 0.1%; 0 of 1 | 7%; 99% | 93% | 536 |
| COLO829T Fiber-seq | 2,637 | 0%; - | 7%; 75% | 93% | 287 |
| COLO829T ONT | 602 | 0%; - | 21%; 98% | 79% | 467 |
| HCC1395T Illumina | 663 | 0%; - | 64%; 100% | 36% | 495 |
| HCC1395T PacBio HiFi | 627 | 0%; - | 73%; 99% | 27% | 413 |
| HCC1395T ONT | 610 | 0%; - | 85%; 100% | 15% | 733 |

## Second pair: INDEL, FP after population filtering

| Dataset | GRCh38-graph FP after population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP after population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 190 | 6%; 0 of 12 | 80%; 78% | 14% | 108 |
| HG008T PacBio HiFi | 541 | 9%; 6.4% | 48%; 59% | 43% | 107 |
| HG008T ONT | 234 | 1%; 0 of 3 | 19%; 91% | 80% | 43 |
| COLO829T Illumina | 79 | 1%; 0 of 1 | 57%; 76% | 42% | 48 |
| COLO829T Fiber-seq | 274 | 4%; 0 of 11 | 28%; 25% | 68% | 26 |
| COLO829T ONT | 55 | 0%; - | 4%; 1 of 2 | 96% | 45 |
| HCC1395T Illumina | 92 | 0%; - | 53%; 98% | 47% | 89 |
| HCC1395T PacBio HiFi | 198 | 0.5%; 0 of 1 | 37%; 92% | 63% | 83 |
| HCC1395T ONT | 168 | 0%; - | 17%; 100% | 83% | 83 |

## Counts behind the four tables

n = GRCh38-only FP alleles of the class; candidate on v2.1 = of them, alleles with a tensor in the v2.1 set; FP on v2.1 = of them, alleles that are also a PASS false positive of the v2.1 model at the same stage.

### SNV, FP before population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 697 | 28,178 (28,178) | 22,744; 56; 0 | 4,569; 3,990; 138 | 865; 352; 21 | 468 (468) | 159 |
| HG008T PacBio HiFi | 697 | 84,853 (84,853) | 72,172; 243; 0 | 11,385; 9,884; 323 | 1,296; 167; 4 | 481 (481) | 327 |
| HG008T ONT | 697 | 50,403 (50,403) | 41,794; 218; 0 | 7,870; 6,889; 150 | 739; 175; 3 | 273 (273) | 153 |
| COLO829T Illumina | 2,443 | 29,530 (29,530) | 23,442; 35; 0 | 3,106; 2,766; 80 | 2,982; 1,759; 146 | 550 (550) | 226 |
| COLO829T Fiber-seq | 2,443 | 38,677 (38,677) | 27,599; 171; 0 | 4,603; 3,810; 77 | 6,475; 1,250; 71 | 291 (291) | 148 |
| COLO829T ONT | 2,443 | 30,128 (30,128) | 24,250; 161; 0 | 3,602; 3,232; 61 | 2,276; 971; 97 | 477 (477) | 158 |
| HCC1395T Illumina | 3,440 | 77,905 (77,905) | 66,643; 75; 0 | 10,718; 9,298; 285 | 544; 288; 82 | 502 (502) | 367 |
| HCC1395T PacBio HiFi | 3,440 | 86,347 (86,347) | 74,001; 168; 0 | 11,520; 9,872; 253 | 826; 158; 79 | 416 (416) | 332 |
| HCC1395T ONT | 3,440 | 115,580 (115,580) | 100,803; 706; 0 | 13,723; 11,858; 285 | 1,054; 361; 15 | 780 (780) | 300 |

### INDEL, FP before population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 587 | 16,016 (16,016) | 13,708; 30; 0 | 2,187; 1,424; 40 | 121; 15; 3 | 109 (109) | 43 |
| HG008T PacBio HiFi | 587 | 18,590 (18,590) | 14,615; 2,215; 0 | 3,118; 1,741; 69 | 857; 303; 5 | 115 (115) | 74 |
| HG008T ONT | 587 | 3,642 (3,642) | 2,586; 631; 0 | 568; 373; 7 | 488; 309; 0 | 43 (43) | 7 |
| COLO829T Illumina | 131 | 5,977 (5,977) | 5,055; 2; 0 | 701; 394; 8 | 221; 40; 1 | 49 (49) | 9 |
| COLO829T Fiber-seq | 131 | 1,764 (1,764) | 751; 145; 0 | 378; 94; 1 | 635; 205; 1 | 27 (27) | 2 |
| COLO829T ONT | 131 | 180 (180) | 79; 23; 0 | 12; 8; 0 | 89; 64; 1 | 46 (46) | 1 |
| HCC1395T Illumina | 133 | 6,594 (6,594) | 5,550; 1; 0 | 940; 683; 34 | 104; 17; 9 | 92 (92) | 43 |
| HCC1395T PacBio HiFi | 133 | 10,207 (10,207) | 8,188; 215; 0 | 1,549; 965; 77 | 470; 119; 26 | 126 (126) | 103 |
| HCC1395T ONT | 133 | 3,713 (3,713) | 2,925; 23; 0 | 568; 427; 18 | 220; 152; 24 | 86 (86) | 42 |

### SNV, FP after population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 697 | 517 (517) | 1; 0; 0 | 208; 196; 137 | 308; 219; 21 | 460 (460) | 158 |
| HG008T PacBio HiFi | 697 | 680 (680) | 2; 1; 0 | 469; 450; 322 | 209; 17; 4 | 479 (479) | 326 |
| HG008T ONT | 697 | 429 (429) | 4; 0; 0 | 339; 329; 150 | 86; 33; 3 | 273 (273) | 153 |
| COLO829T Illumina | 2,443 | 1,480 (1,480) | 1; 0; 0 | 106; 105; 75 | 1,373; 1,083; 145 | 536 (536) | 220 |
| COLO829T Fiber-seq | 2,443 | 2,637 (2,637) | 0; 0; 0 | 174; 131; 75 | 2,463; 505; 70 | 287 (287) | 145 |
| COLO829T ONT | 2,443 | 602 (602) | 0; 0; 0 | 127; 125; 61 | 475; 245; 97 | 467 (467) | 158 |
| HCC1395T Illumina | 3,440 | 663 (663) | 0; 0; 0 | 425; 424; 282 | 238; 233; 82 | 495 (495) | 364 |
| HCC1395T PacBio HiFi | 3,440 | 627 (627) | 0; 0; 0 | 456; 453; 253 | 171; 92; 79 | 413 (413) | 332 |
| HCC1395T ONT | 3,440 | 610 (610) | 0; 0; 0 | 517; 516; 285 | 93; 61; 14 | 733 (733) | 299 |

### INDEL, FP after population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 587 | 190 (190) | 12; 0; 0 | 152; 118; 40 | 26; 8; 3 | 108 (108) | 43 |
| HG008T PacBio HiFi | 587 | 541 (541) | 47; 3; 0 | 260; 153; 62 | 234; 138; 4 | 107 (107) | 66 |
| HG008T ONT | 587 | 234 (234) | 3; 0; 0 | 44; 40; 7 | 187; 152; 0 | 43 (43) | 7 |
| COLO829T Illumina | 131 | 79 (79) | 1; 0; 0 | 45; 34; 8 | 33; 14; 0 | 48 (48) | 8 |
| COLO829T Fiber-seq | 131 | 274 (274) | 11; 0; 0 | 77; 19; 1 | 186; 105; 1 | 26 (26) | 2 |
| COLO829T ONT | 131 | 55 (55) | 0; 0; 0 | 2; 1; 0 | 53; 47; 1 | 45 (45) | 1 |
| HCC1395T Illumina | 133 | 92 (92) | 0; 0; 0 | 49; 48; 31 | 43; 11; 9 | 89 (89) | 40 |
| HCC1395T PacBio HiFi | 133 | 198 (198) | 1; 0; 0 | 73; 67; 37 | 124; 56; 24 | 83 (83) | 61 |
| HCC1395T ONT | 133 | 168 (168) | 0; 0; 0 | 28; 28; 15 | 140; 116; 24 | 83 (83) | 39 |

## PASS recall and precision of the two models (report.json)

The two models are not at the same recall: the FP totals of the first and last column go with these recalls.

| kind | Dataset | truth | GRCh38-only before: TP / FP / recall / precision | HPRC v2.1 before: TP / FP / recall / precision | GRCh38-only after: TP / FP / recall / precision | HPRC v2.1 after: TP / FP / recall / precision |
|---|---|---|---|---|---|---|
| SNV | HG008T Illumina | 697 | 381 / 28,178 / 0.547 / 0.013 | 542 / 468 / 0.778 / 0.537 | 370 / 517 / 0.531 / 0.417 | 541 / 460 / 0.776 / 0.540 |
| SNV | HG008T PacBio HiFi | 697 | 548 / 84,853 / 0.786 / 0.006 | 569 / 481 / 0.816 / 0.542 | 528 / 680 / 0.758 / 0.437 | 567 / 479 / 0.813 / 0.542 |
| SNV | HG008T ONT | 697 | 462 / 50,403 / 0.663 / 0.009 | 457 / 273 / 0.656 / 0.626 | 444 / 429 / 0.637 / 0.509 | 457 / 273 / 0.656 / 0.626 |
| SNV | COLO829T Illumina | 2,443 | 1,554 / 29,530 / 0.636 / 0.050 | 1,856 / 550 / 0.760 / 0.771 | 1,509 / 1,480 / 0.618 / 0.505 | 1,848 / 536 / 0.756 / 0.775 |
| SNV | COLO829T Fiber-seq | 2,443 | 1,784 / 38,677 / 0.730 / 0.044 | 1,832 / 291 / 0.750 / 0.863 | 1,732 / 2,637 / 0.709 / 0.396 | 1,830 / 287 / 0.749 / 0.864 |
| SNV | COLO829T ONT | 2,443 | 1,386 / 30,128 / 0.567 / 0.044 | 1,731 / 477 / 0.709 / 0.784 | 1,349 / 602 / 0.552 / 0.691 | 1,729 / 467 / 0.708 / 0.787 |
| SNV | HCC1395T Illumina | 3,440 | 3,127 / 77,905 / 0.909 / 0.039 | 3,124 / 502 / 0.908 / 0.862 | 3,074 / 663 / 0.894 / 0.823 | 3,113 / 495 / 0.905 / 0.863 |
| SNV | HCC1395T PacBio HiFi | 3,440 | 3,063 / 86,347 / 0.890 / 0.034 | 2,959 / 416 / 0.860 / 0.877 | 3,010 / 627 / 0.875 / 0.828 | 2,950 / 413 / 0.858 / 0.877 |
| SNV | HCC1395T ONT | 3,440 | 2,545 / 115,580 / 0.740 / 0.022 | 2,331 / 780 / 0.678 / 0.749 | 2,498 / 610 / 0.726 / 0.804 | 2,326 / 733 / 0.676 / 0.760 |
| INDEL | HG008T Illumina | 587 | 261 / 16,016 / 0.445 / 0.016 | 65 / 109 / 0.111 / 0.374 | 108 / 190 / 0.184 / 0.362 | 65 / 108 / 0.111 / 0.376 |
| INDEL | HG008T PacBio HiFi | 587 | 279 / 18,590 / 0.475 / 0.015 | 96 / 115 / 0.164 / 0.455 | 120 / 541 / 0.204 / 0.182 | 94 / 107 / 0.160 / 0.468 |
| INDEL | HG008T ONT | 587 | 57 / 3,642 / 0.097 / 0.015 | 47 / 43 / 0.080 / 0.522 | 34 / 234 / 0.058 / 0.127 | 47 / 43 / 0.080 / 0.522 |
| INDEL | COLO829T Illumina | 131 | 39 / 5,977 / 0.298 / 0.006 | 17 / 49 / 0.130 / 0.258 | 18 / 79 / 0.137 / 0.186 | 17 / 48 / 0.130 / 0.262 |
| INDEL | COLO829T Fiber-seq | 131 | 26 / 1,764 / 0.198 / 0.015 | 14 / 27 / 0.107 / 0.341 | 14 / 274 / 0.107 / 0.049 | 14 / 26 / 0.107 / 0.350 |
| INDEL | COLO829T ONT | 131 | 1 / 180 / 0.008 / 0.006 | 11 / 46 / 0.084 / 0.193 | 1 / 55 / 0.008 / 0.018 | 11 / 45 / 0.084 / 0.196 |
| INDEL | HCC1395T Illumina | 133 | 99 / 6,594 / 0.744 / 0.015 | 99 / 92 / 0.744 / 0.518 | 93 / 92 / 0.699 / 0.503 | 99 / 89 / 0.744 / 0.527 |
| INDEL | HCC1395T PacBio HiFi | 133 | 106 / 10,207 / 0.797 / 0.010 | 95 / 126 / 0.714 / 0.430 | 93 / 198 / 0.699 / 0.320 | 91 / 83 / 0.684 / 0.523 |
| INDEL | HCC1395T ONT | 133 | 71 / 3,713 / 0.534 / 0.019 | 80 / 86 / 0.602 / 0.482 | 70 / 168 / 0.526 / 0.294 | 80 / 83 / 0.602 / 0.491 |

## Sets and PASS rules

| Dataset | kind | GRCh38-only set | PASS rule (GRCh38-only) | HPRC v2.1 set | PASS rule (v2.1) |
|---|---|---|---|---|---|
| HG008T Illumina | SNV | HG008T_Illumina_linear | vcf_chr1_HG008T: p_somatic >= 0.09534827619791031 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | HG008T_Illumina_d46 | vcf_chr1_HG008T_Illumina_d46: p_somatic >= 0.7630133032798767 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HG008T PacBio HiFi | SNV | HG008T_PacBio_linear | vcf_chr1_HG008T: p_somatic >= 0.08754626661539078 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | HG008T_PacBio_d46 | vcf_chr1_HG008T_PacBio_d46: p_somatic >= 0.8157289028167725 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HG008T ONT | SNV | HG008T_ONT_linear | vcf_chr1_HG008T: p_somatic >= 0.05439719185233116 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | HG008T_ONT_d46 | vcf_chr1_HG008T_ONT_d46: p_somatic >= 0.8221198916435242 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| COLO829T Illumina | SNV | COLO829T_Illumina_linear_snvaf05 | vcf_chr1_COLO829T: p_somatic >= 0.09534827619791031 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | COLO829T_Illumina_d46_snvaf05 | vcf_chr1_COLO829T_Illumina_d46: p_somatic >= 0.7630133032798767 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| COLO829T Fiber-seq | SNV | COLO829T_fiberseq_linear_snvaf05 | vcf_chr1_COLO829T: p_somatic >= 0.08754626661539078 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | COLO829T_fiberseq_d46_snvaf05 | vcf_chr1_COLO829T_fiberseq_d46: p_somatic >= 0.8157289028167725 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| COLO829T ONT | SNV | COLO829T_ONT_linear | vcf_chr1_COLO829T: p_somatic >= 0.05439719185233116 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | COLO829T_ONT_d46 | vcf_chr1_COLO829T_ONT_d46: p_somatic >= 0.8221198916435242 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HCC1395T Illumina | SNV | HCC1395T_Illumina_linear_snvaf05 | vcf_chr1_HCC1395T: p_somatic >= 0.09534827619791031 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | HCC1395T_Illumina_d46_snvaf05 | vcf_chr1_HCC1395T_Illumina_d46: p_somatic >= 0.7630133032798767 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HCC1395T PacBio HiFi | SNV | HCC1395T_PacBio_linear_snvaf05 | vcf_chr1_HCC1395T: p_somatic >= 0.08754626661539078 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | HCC1395T_PacBio_d46_snvaf05 | vcf_chr1_HCC1395T_PacBio_d46: p_somatic >= 0.8157289028167725 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HCC1395T ONT | SNV | HCC1395T_ONT_linear | vcf_chr1_HCC1395T: p_somatic >= 0.05439719185233116 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) | HCC1395T_ONT_d46 | vcf_chr1_HCC1395T_ONT_d46: p_somatic >= 0.8221198916435242 (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HG008T Illumina | INDEL | HG008T_Illumina_linear | vcf_ppr8_chr1_HG008T: p_somatic >= 0.06674449145793915 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) | HG008T_Illumina_d46 | vcf_ppr8_chr1_HG008T_Illumina_d46: p_somatic >= 0.38437458872795105 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HG008T PacBio HiFi | INDEL | HG008T_PacBio_linear | vcf_ppr8_chr1_HG008T: p_somatic >= 0.06694088131189346 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) | HG008T_PacBio_d46 | vcf_ppr8_chr1_HG008T_PacBio_d46: p_somatic >= 0.6854676604270935 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HG008T ONT | INDEL | HG008T_ONT_linear | vcf_pp_chr1_HG008T: p_somatic >= 0.14897342026233673 (validation post-PoN best F1 (val.somatic_threshold_postpon)) | HG008T_ONT_d46 | vcf_pp_chr1_HG008T_ONT_d46: p_somatic >= 0.7489863634109497 (validation post-PoN best F1 (val.somatic_threshold_postpon)) |
| COLO829T Illumina | INDEL | COLO829T_Illumina_linear | vcf_ppr8_chr1_COLO829T: p_somatic >= 0.06674449145793915 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) | COLO829T_Illumina_d46 | vcf_ppr8_chr1_COLO829T_Illumina_d46: p_somatic >= 0.38437458872795105 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| COLO829T Fiber-seq | INDEL | COLO829T_fiberseq_linear | vcf_ppr8_chr1_COLO829T: p_somatic >= 0.06694088131189346 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) | COLO829T_fiberseq_d46 | vcf_ppr8_chr1_COLO829T_fiberseq_d46: p_somatic >= 0.6854676604270935 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| COLO829T ONT | INDEL | COLO829T_ONT_linear | vcf_pp_chr1_COLO829T: p_somatic >= 0.14897342026233673 (validation post-PoN best F1 (val.somatic_threshold_postpon)) | COLO829T_ONT_d46 | vcf_pp_chr1_COLO829T_ONT_d46: p_somatic >= 0.7489863634109497 (validation post-PoN best F1 (val.somatic_threshold_postpon)) |
| HCC1395T Illumina | INDEL | HCC1395T_Illumina_linear | vcf_ppr8_chr1_HCC1395T: p_somatic >= 0.06674449145793915 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) | HCC1395T_Illumina_d46 | vcf_ppr8_chr1_HCC1395T_Illumina_d46: p_somatic >= 0.38437458872795105 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HCC1395T PacBio HiFi | INDEL | HCC1395T_PacBio_linear | vcf_ppr8_chr1_HCC1395T: p_somatic >= 0.06694088131189346 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) | HCC1395T_PacBio_d46 | vcf_ppr8_chr1_HCC1395T_PacBio_d46: p_somatic >= 0.6854676604270935 (validation post-PoN recall 0.8 of the callable positives (val.somatic_at_recall_postpon_callable)) |
| HCC1395T ONT | INDEL | HCC1395T_ONT_linear | vcf_pp_chr1_HCC1395T: p_somatic >= 0.14897342026233673 (validation post-PoN best F1 (val.somatic_threshold_postpon)) | HCC1395T_ONT_d46 | vcf_pp_chr1_HCC1395T_ONT_d46: p_somatic >= 0.7489863634109497 (validation post-PoN best F1 (val.somatic_threshold_postpon)) |
