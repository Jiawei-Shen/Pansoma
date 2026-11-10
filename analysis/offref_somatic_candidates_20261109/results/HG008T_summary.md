# HG008T off-reference somatic candidates (keys)

| group | kind | T1 | T2 | T3_reads_only | asm_only | rejected_normal | other | chr1 T1 | chr1 T2 |
|---|---|---|---|---|---|---|---|---|---|
| ctrl_germline | SNV | 0 | 0 | 0 | 0 | 2512 | 74 | 0 | 0 |
| ctrl_germline | INDEL | 0 | 9 | 0 | 8 | 1412 | 1177 | 0 | 0 |
| ctrl_non | SNV | 0 | 0 | 0 | 4 | 378 | 935 | 0 | 0 |
| ctrl_non | INDEL | 0 | 4 | 4 | 0 | 369 | 940 | 0 | 0 |
| ctrl_somatic | SNV | 1326 | 267 | 0 | 12 | 13 | 38 | 75 | 15 |
| ctrl_somatic | INDEL | 217 | 1391 | 1 | 52 | 71 | 256 | 11 | 78 |
| offref | SNV | 23 | 167 | 9 | 2378 | 39571 | 100629 | 1 | 21 |
| offref | INDEL | 5 | 882 | 88 | 1692 | 26381 | 229501 | 0 | 78 |

## Off-reference candidates by platform sets with a tensor of the key

| tier | kind | keys | events (anchors within 10 bp merged) | Illumina | ONT | PacBio | Illumina,ONT | Illumina,PacBio | ONT,PacBio | Illumina,ONT,PacBio |
|---|---|---|---|---|---|---|---|---|---|---|
| T1 | SNV | 23 | 22 | 6 | 0 | 1 | 0 | 1 | 2 | 13 |
| T1 | INDEL | 5 | 5 | 1 | 0 | 2 | 0 | 0 | 0 | 2 |
| T2 | SNV | 167 | 146 | 49 | 5 | 12 | 0 | 23 | 19 | 59 |
| T2 | INDEL | 882 | 836 | 145 | 19 | 65 | 7 | 223 | 69 | 354 |
| T3_reads_only | SNV | 9 | 9 | 1 | 1 | 2 | 0 | 0 | 4 | 1 |
| T3_reads_only | INDEL | 88 | 88 | 0 | 16 | 29 | 0 | 0 | 43 | 0 |

## Off-reference keys by label reason (in-BED keys without truth = FDR proxy; residual_partial = next to GRCh38 somatic truth)

| reason | kind | keys | T1 | T2 | T3_reads_only | asm_only | rejected_normal | other |
|---|---|---|---|---|---|---|---|---|
| not_on_unique_grch38_node | SNV | 48686 | 1 | 8 | 0 | 2031 | 8347 | 38299 |
| not_on_unique_grch38_node | INDEL | 13036 | 0 | 3 | 1 | 297 | 1347 | 11388 |
| off_reference_no_truth_match | SNV | 71281 | 5 | 22 | 8 | 76 | 25264 | 45906 |
| off_reference_no_truth_match | INDEL | 222137 | 0 | 52 | 85 | 315 | 20890 | 200795 |
| outside_confident_region | SNV | 22394 | 6 | 36 | 1 | 246 | 5834 | 16271 |
| outside_confident_region | INDEL | 20848 | 0 | 138 | 2 | 958 | 3779 | 15971 |
| residual_partial_somatic_truth | SNV | 416 | 11 | 101 | 0 | 25 | 126 | 153 |
| residual_partial_somatic_truth | INDEL | 2067 | 4 | 687 | 0 | 121 | 244 | 1011 |
| truth_matches_non_representative_allele | INDEL | 461 | 1 | 2 | 0 | 1 | 121 | 336 |

## Known off-reference somatic keys

| known | T1 | T2 | T3_reads_only | asm_only | rejected_normal | other |
|---|---|---|---|---|---|---|
| asm_truth | 10 | 90 | 0 | 6 | 48 | 25 |

## Funnel (rows): set, kind, group, asm_class, normal_side

| set | kind | group | asm_class | normal_side | rows |
|---|---|---|---|---|---|
| HG008T_Illumina_d46 | INDEL | ctrl_germline | neither | no_site_diff | 3 |
| HG008T_Illumina_d46 | INDEL | ctrl_germline | neither | ref_exact | 1 |
| HG008T_Illumina_d46 | INDEL | ctrl_germline | neither | site_plus_other | 4 |
| HG008T_Illumina_d46 | INDEL | ctrl_germline | normal |  | 803 |
| HG008T_Illumina_d46 | INDEL | ctrl_germline | short_flank |  | 59 |
| HG008T_Illumina_d46 | INDEL | ctrl_germline | tumor_only | ref_exact | 9 |
| HG008T_Illumina_d46 | INDEL | ctrl_germline | tumor_only | site_only | 1 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | neither | no_alignment | 17 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | neither | no_site_diff | 13 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | neither | ref_exact | 11 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | neither | site_only | 6 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | neither | site_plus_other | 5 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | normal |  | 335 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | short_flank |  | 46 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | tumor_only | no_site_diff | 2 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | tumor_only | ref_exact | 4 |
| HG008T_Illumina_d46 | INDEL | ctrl_non | tumor_only | site_only | 1 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | neither | no_alignment | 2 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | neither | site_plus_other | 7 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | normal |  | 4 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | short_flank |  | 9 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | tumor_only | no_alignment | 1 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | tumor_only | ref_exact | 787 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | tumor_only | site_only | 44 |
| HG008T_Illumina_d46 | INDEL | ctrl_somatic | tumor_only | site_plus_other | 8 |
| HG008T_Illumina_d46 | INDEL | offref | neither | no_alignment | 547 |
| HG008T_Illumina_d46 | INDEL | offref | neither | no_site_diff | 719 |
| HG008T_Illumina_d46 | INDEL | offref | neither | ref_exact | 184 |
| HG008T_Illumina_d46 | INDEL | offref | neither | site_only | 356 |
| HG008T_Illumina_d46 | INDEL | offref | neither | site_plus_other | 599 |
| HG008T_Illumina_d46 | INDEL | offref | no_query |  | 22 |
| HG008T_Illumina_d46 | INDEL | offref | normal |  | 17626 |
| HG008T_Illumina_d46 | INDEL | offref | short_flank |  | 5277 |
| HG008T_Illumina_d46 | INDEL | offref | tumor_only | no_alignment | 26 |
| HG008T_Illumina_d46 | INDEL | offref | tumor_only | no_site_diff | 23 |
| HG008T_Illumina_d46 | INDEL | offref | tumor_only | ref_exact | 905 |
| HG008T_Illumina_d46 | INDEL | offref | tumor_only | site_only | 519 |
| HG008T_Illumina_d46 | INDEL | offref | tumor_only | site_plus_other | 111 |
| HG008T_Illumina_d46 | SNV | ctrl_germline | neither | no_alignment | 11 |
| HG008T_Illumina_d46 | SNV | ctrl_germline | neither | no_site_diff | 5 |
| HG008T_Illumina_d46 | SNV | ctrl_germline | neither | ref_exact | 2 |
| HG008T_Illumina_d46 | SNV | ctrl_germline | neither | site_only | 1 |
| HG008T_Illumina_d46 | SNV | ctrl_germline | normal |  | 855 |
| HG008T_Illumina_d46 | SNV | ctrl_germline | short_flank |  | 5 |
| HG008T_Illumina_d46 | SNV | ctrl_germline | tumor_only | site_plus_other | 1 |
| HG008T_Illumina_d46 | SNV | ctrl_non | neither | no_alignment | 73 |
| HG008T_Illumina_d46 | SNV | ctrl_non | neither | no_site_diff | 3 |
| HG008T_Illumina_d46 | SNV | ctrl_non | neither | ref_exact | 177 |
| HG008T_Illumina_d46 | SNV | ctrl_non | neither | site_only | 3 |
| HG008T_Illumina_d46 | SNV | ctrl_non | neither | site_plus_other | 160 |
| HG008T_Illumina_d46 | SNV | ctrl_non | normal |  | 21 |
| HG008T_Illumina_d46 | SNV | ctrl_non | short_flank |  | 3 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | neither | site_only | 1 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | neither | site_plus_other | 3 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | normal |  | 14 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | short_flank |  | 2 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | no_alignment | 2 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | no_site_diff | 1 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | ref_exact | 809 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | site_only | 42 |
| HG008T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | site_plus_other | 6 |
| HG008T_Illumina_d46 | SNV | offref | neither | no_alignment | 25642 |
| HG008T_Illumina_d46 | SNV | offref | neither | no_site_diff | 6630 |
| HG008T_Illumina_d46 | SNV | offref | neither | ref_exact | 25308 |
| HG008T_Illumina_d46 | SNV | offref | neither | site_only | 10218 |
| HG008T_Illumina_d46 | SNV | offref | neither | site_plus_other | 37653 |
| HG008T_Illumina_d46 | SNV | offref | normal |  | 53591 |
| HG008T_Illumina_d46 | SNV | offref | short_flank |  | 12065 |
| HG008T_Illumina_d46 | SNV | offref | tumor_only | no_alignment | 149 |
| HG008T_Illumina_d46 | SNV | offref | tumor_only | no_site_diff | 110 |
| HG008T_Illumina_d46 | SNV | offref | tumor_only | ref_exact | 240 |
| HG008T_Illumina_d46 | SNV | offref | tumor_only | site_only | 208 |
| HG008T_Illumina_d46 | SNV | offref | tumor_only | site_plus_other | 101 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | neither | no_alignment | 12 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | neither | no_site_diff | 3 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | neither | ref_exact | 463 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | neither | site_only | 48 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | neither | site_plus_other | 9 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | no_query |  | 1 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | normal |  | 309 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | short_flank |  | 26 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | tumor_only | ref_exact | 8 |
| HG008T_ONT_d46 | INDEL | ctrl_germline | tumor_only | site_only | 1 |
| HG008T_ONT_d46 | INDEL | ctrl_non | neither | no_alignment | 10 |
| HG008T_ONT_d46 | INDEL | ctrl_non | neither | ref_exact | 400 |
| HG008T_ONT_d46 | INDEL | ctrl_non | neither | site_only | 11 |
| HG008T_ONT_d46 | INDEL | ctrl_non | neither | site_plus_other | 11 |
| HG008T_ONT_d46 | INDEL | ctrl_non | normal |  | 7 |
| HG008T_ONT_d46 | INDEL | ctrl_non | short_flank |  | 1 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | neither | no_alignment | 4 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | neither | ref_exact | 54 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | neither | site_only | 88 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | neither | site_plus_other | 9 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | normal |  | 8 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | short_flank |  | 11 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | tumor_only | no_alignment | 1 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | tumor_only | ref_exact | 659 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | tumor_only | site_only | 27 |
| HG008T_ONT_d46 | INDEL | ctrl_somatic | tumor_only | site_plus_other | 2 |
| HG008T_ONT_d46 | INDEL | offref | neither | no_alignment | 4735 |
| HG008T_ONT_d46 | INDEL | offref | neither | no_site_diff | 1216 |
| HG008T_ONT_d46 | INDEL | offref | neither | ref_exact | 97368 |
| HG008T_ONT_d46 | INDEL | offref | neither | site_only | 67643 |
| HG008T_ONT_d46 | INDEL | offref | neither | site_plus_other | 6096 |
| HG008T_ONT_d46 | INDEL | offref | no_query |  | 2288 |
| HG008T_ONT_d46 | INDEL | offref | normal |  | 32183 |
| HG008T_ONT_d46 | INDEL | offref | short_flank |  | 12147 |
| HG008T_ONT_d46 | INDEL | offref | tumor_only | no_alignment | 180 |
| HG008T_ONT_d46 | INDEL | offref | tumor_only | no_site_diff | 25 |
| HG008T_ONT_d46 | INDEL | offref | tumor_only | ref_exact | 1376 |
| HG008T_ONT_d46 | INDEL | offref | tumor_only | site_only | 540 |
| HG008T_ONT_d46 | INDEL | offref | tumor_only | site_plus_other | 30 |
| HG008T_ONT_d46 | SNV | ctrl_germline | neither | no_alignment | 9 |
| HG008T_ONT_d46 | SNV | ctrl_germline | neither | no_site_diff | 14 |
| HG008T_ONT_d46 | SNV | ctrl_germline | neither | ref_exact | 5 |
| HG008T_ONT_d46 | SNV | ctrl_germline | neither | site_only | 10 |
| HG008T_ONT_d46 | SNV | ctrl_germline | neither | site_plus_other | 4 |
| HG008T_ONT_d46 | SNV | ctrl_germline | normal |  | 837 |
| HG008T_ONT_d46 | SNV | ctrl_germline | tumor_only | ref_exact | 1 |
| HG008T_ONT_d46 | SNV | ctrl_non | neither | no_alignment | 60 |
| HG008T_ONT_d46 | SNV | ctrl_non | neither | no_site_diff | 13 |
| HG008T_ONT_d46 | SNV | ctrl_non | neither | ref_exact | 135 |
| HG008T_ONT_d46 | SNV | ctrl_non | neither | site_only | 32 |
| HG008T_ONT_d46 | SNV | ctrl_non | neither | site_plus_other | 60 |
| HG008T_ONT_d46 | SNV | ctrl_non | normal |  | 136 |
| HG008T_ONT_d46 | SNV | ctrl_non | short_flank |  | 4 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | neither | no_alignment | 7 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | neither | ref_exact | 2 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | neither | site_only | 2 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | neither | site_plus_other | 4 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | normal |  | 7 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | tumor_only | no_site_diff | 1 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | tumor_only | ref_exact | 815 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | tumor_only | site_only | 37 |
| HG008T_ONT_d46 | SNV | ctrl_somatic | tumor_only | site_plus_other | 5 |
| HG008T_ONT_d46 | SNV | offref | neither | no_alignment | 12238 |
| HG008T_ONT_d46 | SNV | offref | neither | no_site_diff | 9483 |
| HG008T_ONT_d46 | SNV | offref | neither | ref_exact | 3799 |
| HG008T_ONT_d46 | SNV | offref | neither | site_only | 10301 |
| HG008T_ONT_d46 | SNV | offref | neither | site_plus_other | 8858 |
| HG008T_ONT_d46 | SNV | offref | no_query |  | 64 |
| HG008T_ONT_d46 | SNV | offref | normal |  | 104535 |
| HG008T_ONT_d46 | SNV | offref | short_flank |  | 15563 |
| HG008T_ONT_d46 | SNV | offref | tumor_only | no_alignment | 710 |
| HG008T_ONT_d46 | SNV | offref | tumor_only | no_site_diff | 102 |
| HG008T_ONT_d46 | SNV | offref | tumor_only | ref_exact | 1004 |
| HG008T_ONT_d46 | SNV | offref | tumor_only | site_only | 183 |
| HG008T_ONT_d46 | SNV | offref | tumor_only | site_plus_other | 95 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | neither | no_alignment | 2 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | neither | no_site_diff | 2 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | neither | ref_exact | 509 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | neither | site_only | 18 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | neither | site_plus_other | 1 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | normal |  | 318 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | short_flank |  | 25 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | tumor_only | ref_exact | 4 |
| HG008T_PacBio_d46 | INDEL | ctrl_germline | tumor_only | site_only | 1 |
| HG008T_PacBio_d46 | INDEL | ctrl_non | neither | ref_exact | 432 |
| HG008T_PacBio_d46 | INDEL | ctrl_non | neither | site_plus_other | 2 |
| HG008T_PacBio_d46 | INDEL | ctrl_non | normal |  | 6 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | neither | ref_exact | 19 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | neither | site_only | 75 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | normal |  | 7 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | short_flank |  | 6 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | tumor_only | no_site_diff | 1 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | tumor_only | ref_exact | 733 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | tumor_only | site_only | 27 |
| HG008T_PacBio_d46 | INDEL | ctrl_somatic | tumor_only | site_plus_other | 2 |
| HG008T_PacBio_d46 | INDEL | offref | neither | no_alignment | 1556 |
| HG008T_PacBio_d46 | INDEL | offref | neither | no_site_diff | 953 |
| HG008T_PacBio_d46 | INDEL | offref | neither | ref_exact | 76276 |
| HG008T_PacBio_d46 | INDEL | offref | neither | site_only | 29155 |
| HG008T_PacBio_d46 | INDEL | offref | neither | site_plus_other | 1757 |
| HG008T_PacBio_d46 | INDEL | offref | no_query |  | 614 |
| HG008T_PacBio_d46 | INDEL | offref | normal |  | 37523 |
| HG008T_PacBio_d46 | INDEL | offref | short_flank |  | 9218 |
| HG008T_PacBio_d46 | INDEL | offref | tumor_only | no_alignment | 70 |
| HG008T_PacBio_d46 | INDEL | offref | tumor_only | no_site_diff | 21 |
| HG008T_PacBio_d46 | INDEL | offref | tumor_only | ref_exact | 1684 |
| HG008T_PacBio_d46 | INDEL | offref | tumor_only | site_only | 934 |
| HG008T_PacBio_d46 | INDEL | offref | tumor_only | site_plus_other | 27 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | neither | no_alignment | 4 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | neither | no_site_diff | 7 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | neither | ref_exact | 4 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | neither | site_only | 8 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | neither | site_plus_other | 2 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | normal |  | 850 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | short_flank |  | 1 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | tumor_only | no_site_diff | 1 |
| HG008T_PacBio_d46 | SNV | ctrl_germline | tumor_only | ref_exact | 3 |
| HG008T_PacBio_d46 | SNV | ctrl_non | neither | no_alignment | 10 |
| HG008T_PacBio_d46 | SNV | ctrl_non | neither | no_site_diff | 16 |
| HG008T_PacBio_d46 | SNV | ctrl_non | neither | site_only | 48 |
| HG008T_PacBio_d46 | SNV | ctrl_non | neither | site_plus_other | 4 |
| HG008T_PacBio_d46 | SNV | ctrl_non | normal |  | 350 |
| HG008T_PacBio_d46 | SNV | ctrl_non | short_flank |  | 6 |
| HG008T_PacBio_d46 | SNV | ctrl_non | tumor_only | no_alignment | 2 |
| HG008T_PacBio_d46 | SNV | ctrl_non | tumor_only | ref_exact | 4 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | neither | site_only | 2 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | neither | site_plus_other | 1 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | normal |  | 8 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | tumor_only | no_alignment | 1 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | tumor_only | no_site_diff | 1 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | tumor_only | ref_exact | 821 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | tumor_only | site_only | 42 |
| HG008T_PacBio_d46 | SNV | ctrl_somatic | tumor_only | site_plus_other | 4 |
| HG008T_PacBio_d46 | SNV | offref | neither | no_alignment | 10061 |
| HG008T_PacBio_d46 | SNV | offref | neither | no_site_diff | 8129 |
| HG008T_PacBio_d46 | SNV | offref | neither | ref_exact | 2252 |
| HG008T_PacBio_d46 | SNV | offref | neither | site_only | 8334 |
| HG008T_PacBio_d46 | SNV | offref | neither | site_plus_other | 7758 |
| HG008T_PacBio_d46 | SNV | offref | no_query |  | 48 |
| HG008T_PacBio_d46 | SNV | offref | normal |  | 107765 |
| HG008T_PacBio_d46 | SNV | offref | short_flank |  | 5452 |
| HG008T_PacBio_d46 | SNV | offref | tumor_only | no_alignment | 232 |
| HG008T_PacBio_d46 | SNV | offref | tumor_only | no_site_diff | 81 |
| HG008T_PacBio_d46 | SNV | offref | tumor_only | ref_exact | 274 |
| HG008T_PacBio_d46 | SNV | offref | tumor_only | site_only | 237 |
| HG008T_PacBio_d46 | SNV | offref | tumor_only | site_plus_other | 62 |
