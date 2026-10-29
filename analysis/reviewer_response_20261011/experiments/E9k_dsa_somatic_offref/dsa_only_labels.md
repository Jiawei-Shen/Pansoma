# E9k part 4b: labels of the tensors of the DSA-only / germline-allele somatic SNVs (dsa_only_tensors.py labels)

A and B events of dsa_only_tensors.md with a single-base key. GRCh38 position of an A key: ReferencePath.linear, as the
labels. Label region: dip.bed (somatic label 0 region, COLO829T_truth/README.txt). Truth: SMaHT COLO829T (44,005 SNVs).

## d9

| where the edit lands | event class | events | in the SMaHT truth | COLO829T_Illumina: tensor (somatic / non / germline / ignore) | COLO829T_fiberseq: tensor (somatic / non / germline / ignore) | COLO829T_ONT: tensor (somatic / non / germline / ignore) |
|---|---|---:|---:|---:|---:|---:|
| chr1-22, label region | dsa_only | 194 | 17 | 40 (20 / 18 / 2 / 0) | 30 (17 / 11 / 2 / 0) | 25 (15 / 8 / 2 / 0) |
| chr1-22, label region | grch38_diff_ref | 9 | 0 | 9 (0 / 3 / 6 / 0) | 5 (0 / 2 / 3 / 0) | 5 (0 / 2 / 3 / 0) |
| chr1-22, outside the label region | dsa_only | 1,273 | 105 | 226 (95 / 0 / 0 / 131) | 118 (77 / 0 / 0 / 41) | 119 (76 / 0 / 0 / 43) |
| chr1-22, outside the label region | grch38_diff_ref | 7 | 1 | 3 (1 / 0 / 0 / 2) | 4 (1 / 0 / 0 / 3) | 5 (1 / 0 / 0 / 4) |
| chrX / chrY | dsa_only | 30 | 3 | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) |
| other GRCh38 contig (Un / random / M) | dsa_only | 550 | 0 | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) |
| off-reference node | dsa_only | 171 | 0 | 126 (5 / 0 / 0 / 121) | 132 (5 / 0 / 0 / 127) | 131 (5 / 0 / 0 / 126) |
| off-reference node | grch38_diff_ref | 30 | 0 | 27 (1 / 0 / 0 / 26) | 17 (1 / 0 / 0 / 16) | 17 (1 / 0 / 0 / 16) |

Reasons of the tensors' labels (all events of this graph):

- COLO829T_Illumina: ignore / outside_confident_region 151; somatic / representative_allele_is_somatic_truth 113; ignore / off_reference_no_truth_match 103; ignore / not_on_unique_grch38_node 26; non / confident_no_truth_allele 14; somatic / residual_partial_somatic_truth 9; germline / representative_allele_is_germline_truth 8; non / near_truth_allele_mismatch 7
- COLO829T_fiberseq: ignore / off_reference_no_truth_match 100; somatic / representative_allele_is_somatic_truth 92; ignore / outside_confident_region 61; ignore / not_on_unique_grch38_node 26; somatic / residual_partial_somatic_truth 9; non / confident_no_truth_allele 8; non / near_truth_allele_mismatch 5; germline / representative_allele_is_germline_truth 5
- COLO829T_ONT: ignore / off_reference_no_truth_match 97; somatic / representative_allele_is_somatic_truth 89; ignore / outside_confident_region 66; ignore / not_on_unique_grch38_node 26; somatic / residual_partial_somatic_truth 9; non / near_truth_allele_mismatch 5; non / confident_no_truth_allele 5; germline / representative_allele_is_germline_truth 5

Baseline, all COLO829T SMaHT truth SNVs on chr1-22 in the label region (41,230): the representative allele of a tensor labelled somatic

- COLO829T_Illumina: 37,286 (90.4%)
- COLO829T_fiberseq: 37,237 (90.3%)
- COLO829T_ONT: 36,948 (89.6%)

## d46

| where the edit lands | event class | events | in the SMaHT truth | COLO829T_Illumina_d46: tensor (somatic / non / germline / ignore) | COLO829T_fiberseq_d46: tensor (somatic / non / germline / ignore) | COLO829T_ONT_d46: tensor (somatic / non / germline / ignore) |
|---|---|---:|---:|---:|---:|---:|
| chr1-22, label region | dsa_only | 138 | 18 | 40 (19 / 20 / 1 / 0) | 30 (16 / 13 / 1 / 0) | 27 (15 / 11 / 1 / 0) |
| chr1-22, label region | grch38_diff_ref | 5 | 0 | 4 (0 / 0 / 4 / 0) | 1 (0 / 0 / 1 / 0) | 1 (0 / 0 / 1 / 0) |
| chr1-22, outside the label region | dsa_only | 1,078 | 105 | 213 (94 / 0 / 0 / 119) | 137 (77 / 0 / 0 / 60) | 118 (74 / 0 / 0 / 44) |
| chr1-22, outside the label region | grch38_diff_ref | 6 | 0 | 2 (0 / 0 / 0 / 2) | 2 (0 / 0 / 0 / 2) | 4 (0 / 0 / 0 / 4) |
| chrX / chrY | dsa_only | 30 | 3 | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) |
| other GRCh38 contig (Un / random / M) | dsa_only | 516 | 0 | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) | 0 (0 / 0 / 0 / 0) |
| off-reference node | dsa_only | 398 | 0 | 139 (4 / 0 / 0 / 135) | 155 (4 / 0 / 0 / 151) | 142 (4 / 0 / 0 / 138) |
| off-reference node | grch38_diff_ref | 29 | 0 | 27 (1 / 0 / 0 / 26) | 22 (1 / 0 / 0 / 21) | 20 (1 / 0 / 0 / 19) |

Reasons of the tensors' labels (all events of this graph):

- COLO829T_Illumina_d46: ignore / outside_confident_region 134; somatic / representative_allele_is_somatic_truth 111; ignore / off_reference_no_truth_match 98; ignore / not_on_unique_grch38_node 50; non / confident_no_truth_allele 12; non / near_truth_allele_mismatch 8; somatic / residual_partial_somatic_truth 7; germline / representative_allele_is_germline_truth 5
- COLO829T_fiberseq_d46: ignore / off_reference_no_truth_match 98; somatic / representative_allele_is_somatic_truth 91; ignore / outside_confident_region 76; ignore / not_on_unique_grch38_node 60; non / near_truth_allele_mismatch 7; somatic / residual_partial_somatic_truth 7; non / confident_no_truth_allele 6; germline / representative_allele_is_germline_truth 2
- COLO829T_ONT_d46: ignore / off_reference_no_truth_match 94; somatic / representative_allele_is_somatic_truth 87; ignore / outside_confident_region 59; ignore / not_on_unique_grch38_node 52; somatic / residual_partial_somatic_truth 7; non / near_truth_allele_mismatch 6; non / confident_no_truth_allele 5; germline / representative_allele_is_germline_truth 2

Baseline, all COLO829T SMaHT truth SNVs on chr1-22 in the label region (41,230): the representative allele of a tensor labelled somatic

- COLO829T_Illumina_d46: 37,297 (90.5%)
- COLO829T_fiberseq_d46: 37,243 (90.3%)
- COLO829T_ONT_d46: 36,943 (89.6%)

