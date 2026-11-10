# COLO829T off-reference somatic candidates (keys)

| group | kind | T1 | T2 | T3_reads_only | asm_only | rejected_normal | other | chr1 T1 | chr1 T2 |
|---|---|---|---|---|---|---|---|---|---|
| ctrl_germline | SNV | 0 | 1 | 0 | 0 | 2524 | 71 | 0 | 0 |
| ctrl_germline | INDEL | 0 | 2 | 0 | 7 | 1418 | 1177 | 0 | 0 |
| ctrl_non | SNV | 9 | 7 | 1 | 6 | 577 | 713 | 0 | 0 |
| ctrl_non | INDEL | 0 | 4 | 2 | 2 | 489 | 821 | 0 | 0 |
| ctrl_somatic | SNV | 1697 | 343 | 169 | 16 | 5 | 99 | 80 | 16 |
| ctrl_somatic | INDEL | 160 | 673 | 35 | 29 | 17 | 273 | 4 | 31 |
| offref | SNV | 110 | 182 | 54 | 3151 | 2827 | 40766 | 7 | 12 |
| offref | INDEL | 1 | 330 | 396 | 1850 | 5377 | 292187 | 0 | 15 |

## Off-reference candidates by platform sets with a tensor of the key

| tier | kind | keys | events (anchors within 10 bp merged) | Illumina | ONT | fiberseq | Illumina,ONT | Illumina,fiberseq | ONT,fiberseq | Illumina,ONT,fiberseq |
|---|---|---|---|---|---|---|---|---|---|---|
| T1 | SNV | 110 | 108 | 5 | 1 | 3 | 3 | 2 | 19 | 77 |
| T1 | INDEL | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| T2 | SNV | 182 | 166 | 45 | 19 | 9 | 2 | 9 | 30 | 68 |
| T2 | INDEL | 330 | 314 | 48 | 7 | 29 | 0 | 51 | 55 | 140 |
| T3_reads_only | SNV | 54 | 52 | 3 | 0 | 8 | 1 | 4 | 27 | 11 |
| T3_reads_only | INDEL | 396 | 396 | 0 | 9 | 204 | 0 | 1 | 182 | 0 |

## Off-reference keys by label reason (in-BED keys without truth = FDR proxy; residual_partial = next to GRCh38 somatic truth)

| reason | kind | keys | T1 | T2 | T3_reads_only | asm_only | rejected_normal | other |
|---|---|---|---|---|---|---|---|---|
| not_on_unique_grch38_node | SNV | 8427 | 20 | 51 | 3 | 1523 | 96 | 6734 |
| not_on_unique_grch38_node | INDEL | 13824 | 0 | 4 | 1 | 391 | 89 | 13339 |
| off_reference_no_truth_match | SNV | 34907 | 66 | 87 | 46 | 557 | 2463 | 31688 |
| off_reference_no_truth_match | INDEL | 278735 | 1 | 238 | 381 | 1273 | 5070 | 271772 |
| outside_confident_region | SNV | 3655 | 8 | 8 | 2 | 1064 | 258 | 2315 |
| outside_confident_region | INDEL | 6553 | 0 | 3 | 13 | 180 | 207 | 6150 |
| residual_partial_somatic_truth | SNV | 101 | 16 | 36 | 3 | 7 | 10 | 29 |
| residual_partial_somatic_truth | INDEL | 291 | 0 | 84 | 1 | 6 | 10 | 190 |
| truth_matches_non_representative_allele | INDEL | 738 | 0 | 1 | 0 | 0 | 1 | 736 |

## Known off-reference somatic keys

| known | T1 | T2 | T3_reads_only | asm_only | rejected_normal | other |
|---|---|---|---|---|---|---|
| dsa_somatic | 32 | 14 | 2 | 32 | 4 | 16 |
| smaht_truth | 7 | 1 | 1 | 2 | 2 | 3 |

## Funnel (rows): set, kind, group, asm_class, normal_side

| set | kind | group | asm_class | normal_side | rows |
|---|---|---|---|---|---|
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | neither | no_alignment | 4 |
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | neither | no_site_diff | 11 |
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | neither | ref_exact | 1 |
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | neither | site_plus_other | 2 |
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | normal |  | 809 |
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | short_flank |  | 51 |
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | tumor_only | ref_exact | 1 |
| COLO829T_Illumina_d46 | INDEL | ctrl_germline | tumor_only | site_only | 1 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | neither | no_alignment | 9 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | neither | no_site_diff | 10 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | neither | ref_exact | 31 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | neither | site_only | 4 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | neither | site_plus_other | 9 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | normal |  | 307 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | short_flank |  | 28 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | tumor_only | no_alignment | 1 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | tumor_only | ref_exact | 36 |
| COLO829T_Illumina_d46 | INDEL | ctrl_non | tumor_only | site_only | 5 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | neither | no_alignment | 1 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | neither | ref_exact | 119 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | neither | site_plus_other | 2 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | normal |  | 5 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | short_flank |  | 7 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | tumor_only | ref_exact | 629 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | tumor_only | site_only | 16 |
| COLO829T_Illumina_d46 | INDEL | ctrl_somatic | tumor_only | site_plus_other | 2 |
| COLO829T_Illumina_d46 | INDEL | offref | neither | no_alignment | 571 |
| COLO829T_Illumina_d46 | INDEL | offref | neither | no_site_diff | 786 |
| COLO829T_Illumina_d46 | INDEL | offref | neither | ref_exact | 245 |
| COLO829T_Illumina_d46 | INDEL | offref | neither | site_only | 280 |
| COLO829T_Illumina_d46 | INDEL | offref | neither | site_plus_other | 692 |
| COLO829T_Illumina_d46 | INDEL | offref | no_query |  | 19 |
| COLO829T_Illumina_d46 | INDEL | offref | normal |  | 17099 |
| COLO829T_Illumina_d46 | INDEL | offref | short_flank |  | 4649 |
| COLO829T_Illumina_d46 | INDEL | offref | tumor_only | no_alignment | 21 |
| COLO829T_Illumina_d46 | INDEL | offref | tumor_only | no_site_diff | 40 |
| COLO829T_Illumina_d46 | INDEL | offref | tumor_only | ref_exact | 524 |
| COLO829T_Illumina_d46 | INDEL | offref | tumor_only | site_only | 404 |
| COLO829T_Illumina_d46 | INDEL | offref | tumor_only | site_plus_other | 130 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | neither | no_alignment | 14 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | neither | no_site_diff | 5 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | neither | ref_exact | 1 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | neither | site_only | 1 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | neither | site_plus_other | 1 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | normal |  | 851 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | short_flank |  | 4 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | tumor_only | ref_exact | 2 |
| COLO829T_Illumina_d46 | SNV | ctrl_germline | tumor_only | site_only | 1 |
| COLO829T_Illumina_d46 | SNV | ctrl_non | neither | no_alignment | 84 |
| COLO829T_Illumina_d46 | SNV | ctrl_non | neither | no_site_diff | 7 |
| COLO829T_Illumina_d46 | SNV | ctrl_non | neither | ref_exact | 147 |
| COLO829T_Illumina_d46 | SNV | ctrl_non | neither | site_only | 18 |
| COLO829T_Illumina_d46 | SNV | ctrl_non | neither | site_plus_other | 135 |
| COLO829T_Illumina_d46 | SNV | ctrl_non | normal |  | 48 |
| COLO829T_Illumina_d46 | SNV | ctrl_non | tumor_only | ref_exact | 1 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | neither | ref_exact | 77 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | neither | site_only | 6 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | normal |  | 8 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | no_alignment | 1 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | no_site_diff | 1 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | ref_exact | 736 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | site_only | 47 |
| COLO829T_Illumina_d46 | SNV | ctrl_somatic | tumor_only | site_plus_other | 4 |
| COLO829T_Illumina_d46 | SNV | offref | neither | no_alignment | 9238 |
| COLO829T_Illumina_d46 | SNV | offref | neither | no_site_diff | 4830 |
| COLO829T_Illumina_d46 | SNV | offref | neither | ref_exact | 7197 |
| COLO829T_Illumina_d46 | SNV | offref | neither | site_only | 3230 |
| COLO829T_Illumina_d46 | SNV | offref | neither | site_plus_other | 12443 |
| COLO829T_Illumina_d46 | SNV | offref | normal |  | 48510 |
| COLO829T_Illumina_d46 | SNV | offref | short_flank |  | 8857 |
| COLO829T_Illumina_d46 | SNV | offref | tumor_only | no_alignment | 68 |
| COLO829T_Illumina_d46 | SNV | offref | tumor_only | no_site_diff | 304 |
| COLO829T_Illumina_d46 | SNV | offref | tumor_only | ref_exact | 560 |
| COLO829T_Illumina_d46 | SNV | offref | tumor_only | site_only | 352 |
| COLO829T_Illumina_d46 | SNV | offref | tumor_only | site_plus_other | 268 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | neither | no_alignment | 7 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | neither | no_site_diff | 11 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | neither | ref_exact | 483 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | neither | site_only | 65 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | neither | site_plus_other | 12 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | normal |  | 283 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | short_flank |  | 16 |
| COLO829T_ONT_d46 | INDEL | ctrl_germline | tumor_only | ref_exact | 3 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | neither | no_alignment | 9 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | neither | no_site_diff | 1 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | neither | ref_exact | 398 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | neither | site_only | 4 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | neither | site_plus_other | 18 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | normal |  | 8 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | short_flank |  | 1 |
| COLO829T_ONT_d46 | INDEL | ctrl_non | tumor_only | ref_exact | 1 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | neither | no_alignment | 3 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | neither | no_site_diff | 1 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | neither | ref_exact | 122 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | neither | site_only | 76 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | neither | site_plus_other | 4 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | normal |  | 14 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | short_flank |  | 10 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | tumor_only | ref_exact | 536 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | tumor_only | site_only | 10 |
| COLO829T_ONT_d46 | INDEL | ctrl_somatic | tumor_only | site_plus_other | 3 |
| COLO829T_ONT_d46 | INDEL | offref | neither | no_alignment | 11053 |
| COLO829T_ONT_d46 | INDEL | offref | neither | no_site_diff | 3542 |
| COLO829T_ONT_d46 | INDEL | offref | neither | ref_exact | 133510 |
| COLO829T_ONT_d46 | INDEL | offref | neither | site_only | 105360 |
| COLO829T_ONT_d46 | INDEL | offref | neither | site_plus_other | 22922 |
| COLO829T_ONT_d46 | INDEL | offref | no_query |  | 2951 |
| COLO829T_ONT_d46 | INDEL | offref | normal |  | 34093 |
| COLO829T_ONT_d46 | INDEL | offref | short_flank |  | 12684 |
| COLO829T_ONT_d46 | INDEL | offref | tumor_only | no_alignment | 35 |
| COLO829T_ONT_d46 | INDEL | offref | tumor_only | no_site_diff | 92 |
| COLO829T_ONT_d46 | INDEL | offref | tumor_only | ref_exact | 816 |
| COLO829T_ONT_d46 | INDEL | offref | tumor_only | site_only | 516 |
| COLO829T_ONT_d46 | INDEL | offref | tumor_only | site_plus_other | 96 |
| COLO829T_ONT_d46 | SNV | ctrl_germline | neither | no_alignment | 5 |
| COLO829T_ONT_d46 | SNV | ctrl_germline | neither | no_site_diff | 16 |
| COLO829T_ONT_d46 | SNV | ctrl_germline | neither | ref_exact | 11 |
| COLO829T_ONT_d46 | SNV | ctrl_germline | neither | site_only | 10 |
| COLO829T_ONT_d46 | SNV | ctrl_germline | neither | site_plus_other | 4 |
| COLO829T_ONT_d46 | SNV | ctrl_germline | normal |  | 830 |
| COLO829T_ONT_d46 | SNV | ctrl_germline | short_flank |  | 4 |
| COLO829T_ONT_d46 | SNV | ctrl_non | neither | no_alignment | 20 |
| COLO829T_ONT_d46 | SNV | ctrl_non | neither | no_site_diff | 21 |
| COLO829T_ONT_d46 | SNV | ctrl_non | neither | ref_exact | 139 |
| COLO829T_ONT_d46 | SNV | ctrl_non | neither | site_only | 22 |
| COLO829T_ONT_d46 | SNV | ctrl_non | neither | site_plus_other | 23 |
| COLO829T_ONT_d46 | SNV | ctrl_non | normal |  | 204 |
| COLO829T_ONT_d46 | SNV | ctrl_non | short_flank |  | 3 |
| COLO829T_ONT_d46 | SNV | ctrl_non | tumor_only | ref_exact | 6 |
| COLO829T_ONT_d46 | SNV | ctrl_non | tumor_only | site_only | 1 |
| COLO829T_ONT_d46 | SNV | ctrl_non | tumor_only | site_plus_other | 1 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | neither | no_alignment | 4 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | neither | ref_exact | 70 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | neither | site_only | 5 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | neither | site_plus_other | 3 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | normal |  | 14 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | tumor_only | no_alignment | 1 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | tumor_only | no_site_diff | 1 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | tumor_only | ref_exact | 749 |
| COLO829T_ONT_d46 | SNV | ctrl_somatic | tumor_only | site_only | 33 |
| COLO829T_ONT_d46 | SNV | offref | neither | no_alignment | 15900 |
| COLO829T_ONT_d46 | SNV | offref | neither | no_site_diff | 21781 |
| COLO829T_ONT_d46 | SNV | offref | neither | ref_exact | 7689 |
| COLO829T_ONT_d46 | SNV | offref | neither | site_only | 17853 |
| COLO829T_ONT_d46 | SNV | offref | neither | site_plus_other | 15036 |
| COLO829T_ONT_d46 | SNV | offref | no_query |  | 33 |
| COLO829T_ONT_d46 | SNV | offref | normal |  | 119818 |
| COLO829T_ONT_d46 | SNV | offref | short_flank |  | 17395 |
| COLO829T_ONT_d46 | SNV | offref | tumor_only | no_alignment | 122 |
| COLO829T_ONT_d46 | SNV | offref | tumor_only | no_site_diff | 287 |
| COLO829T_ONT_d46 | SNV | offref | tumor_only | ref_exact | 996 |
| COLO829T_ONT_d46 | SNV | offref | tumor_only | site_only | 362 |
| COLO829T_ONT_d46 | SNV | offref | tumor_only | site_plus_other | 300 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | neither | no_alignment | 1 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | neither | no_site_diff | 3 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | neither | ref_exact | 425 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | neither | site_only | 12 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | neither | site_plus_other | 2 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | normal |  | 391 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | short_flank |  | 36 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | tumor_only | ref_exact | 8 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_germline | tumor_only | site_only | 2 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_non | neither | no_alignment | 1 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_non | neither | ref_exact | 427 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_non | neither | site_only | 1 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_non | neither | site_plus_other | 3 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_non | normal |  | 8 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | neither | no_alignment | 1 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | neither | ref_exact | 125 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | neither | site_only | 52 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | neither | site_plus_other | 1 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | normal |  | 10 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | short_flank |  | 11 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | tumor_only | ref_exact | 586 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | tumor_only | site_only | 11 |
| COLO829T_fiberseq_d46 | INDEL | ctrl_somatic | tumor_only | site_plus_other | 1 |
| COLO829T_fiberseq_d46 | INDEL | offref | neither | no_alignment | 1979 |
| COLO829T_fiberseq_d46 | INDEL | offref | neither | no_site_diff | 2009 |
| COLO829T_fiberseq_d46 | INDEL | offref | neither | ref_exact | 89199 |
| COLO829T_fiberseq_d46 | INDEL | offref | neither | site_only | 56024 |
| COLO829T_fiberseq_d46 | INDEL | offref | neither | site_plus_other | 4325 |
| COLO829T_fiberseq_d46 | INDEL | offref | no_query |  | 532 |
| COLO829T_fiberseq_d46 | INDEL | offref | normal |  | 47919 |
| COLO829T_fiberseq_d46 | INDEL | offref | short_flank |  | 11148 |
| COLO829T_fiberseq_d46 | INDEL | offref | tumor_only | no_alignment | 56 |
| COLO829T_fiberseq_d46 | INDEL | offref | tumor_only | no_site_diff | 83 |
| COLO829T_fiberseq_d46 | INDEL | offref | tumor_only | ref_exact | 1410 |
| COLO829T_fiberseq_d46 | INDEL | offref | tumor_only | site_only | 1169 |
| COLO829T_fiberseq_d46 | INDEL | offref | tumor_only | site_plus_other | 79 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | neither | no_alignment | 1 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | neither | no_site_diff | 12 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | neither | ref_exact | 6 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | neither | site_only | 14 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | neither | site_plus_other | 2 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | normal |  | 840 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | short_flank |  | 3 |
| COLO829T_fiberseq_d46 | SNV | ctrl_germline | tumor_only | ref_exact | 2 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | neither | no_alignment | 8 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | neither | no_site_diff | 13 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | neither | ref_exact | 2 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | neither | site_only | 26 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | neither | site_plus_other | 4 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | normal |  | 356 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | short_flank |  | 3 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | tumor_only | no_site_diff | 2 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | tumor_only | ref_exact | 22 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | tumor_only | site_only | 1 |
| COLO829T_fiberseq_d46 | SNV | ctrl_non | tumor_only | site_plus_other | 3 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | neither | no_alignment | 1 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | neither | no_site_diff | 1 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | neither | ref_exact | 78 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | neither | site_only | 4 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | neither | site_plus_other | 1 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | normal |  | 14 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | tumor_only | no_alignment | 1 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | tumor_only | ref_exact | 747 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | tumor_only | site_only | 32 |
| COLO829T_fiberseq_d46 | SNV | ctrl_somatic | tumor_only | site_plus_other | 1 |
| COLO829T_fiberseq_d46 | SNV | offref | neither | no_alignment | 8103 |
| COLO829T_fiberseq_d46 | SNV | offref | neither | no_site_diff | 12632 |
| COLO829T_fiberseq_d46 | SNV | offref | neither | ref_exact | 3172 |
| COLO829T_fiberseq_d46 | SNV | offref | neither | site_only | 13601 |
| COLO829T_fiberseq_d46 | SNV | offref | neither | site_plus_other | 8654 |
| COLO829T_fiberseq_d46 | SNV | offref | no_query |  | 45 |
| COLO829T_fiberseq_d46 | SNV | offref | normal |  | 137573 |
| COLO829T_fiberseq_d46 | SNV | offref | short_flank |  | 8132 |
| COLO829T_fiberseq_d46 | SNV | offref | tumor_only | no_alignment | 107 |
| COLO829T_fiberseq_d46 | SNV | offref | tumor_only | no_site_diff | 259 |
| COLO829T_fiberseq_d46 | SNV | offref | tumor_only | ref_exact | 966 |
| COLO829T_fiberseq_d46 | SNV | offref | tumor_only | site_only | 473 |
| COLO829T_fiberseq_d46 | SNV | offref | tumor_only | site_plus_other | 245 |
