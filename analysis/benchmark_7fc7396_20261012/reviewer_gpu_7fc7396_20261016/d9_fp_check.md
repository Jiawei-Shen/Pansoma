# False calls of the COLO829T Illumina SNV models (chr1, in the BED, after the PoN, PASS)

(`d9_fp_check.py`; windows chr1:121-126 and 143-144 Mb; reason = the label reason of the call's tensor.)


## 1. Reasons of the false calls

| run | recipe | false calls | in the windows | outside_confident_region (label -1) | of them in the windows | other reasons |
|---|---|---:|---:|---:|---:|---|
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024 | d9, 4 GPUs | 1,815 | 1,569 | 1,499 | 1,482 | representative_allele_is_germline_truth 148, confident_no_truth_allele 145, near_truth_allele_mismatch 19, residual_partial_somatic_truth 4 |
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s1 | d9, 4 GPUs | 1,340 | 1,085 | 1,080 | 1,062 | representative_allele_is_germline_truth 169, confident_no_truth_allele 72, near_truth_allele_mismatch 15, residual_partial_somatic_truth 4 |
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s2 | d9, 4 GPUs | 1,900 | 1,659 | 1,612 | 1,603 | representative_allele_is_germline_truth 160, confident_no_truth_allele 108, near_truth_allele_mismatch 16, residual_partial_somatic_truth 4 |
| COLO829_Illumina_SNV_h90 | d9, 1 GPU | 284 | 33 | 35 | 28 | representative_allele_is_germline_truth 169, confident_no_truth_allele 57, near_truth_allele_mismatch 19, residual_partial_somatic_truth 4 |
| COLO829_Illumina_SNV_h90_s1 | d9, 1 GPU | 928 | 675 | 646 | 634 | representative_allele_is_germline_truth 176, confident_no_truth_allele 94, near_truth_allele_mismatch 8, residual_partial_somatic_truth 4 |
| COLO829_Illumina_SNV_h90_s2 | d9, 1 GPU | 790 | 552 | 561 | 544 | representative_allele_is_germline_truth 148, confident_no_truth_allele 60, near_truth_allele_mismatch 18, residual_partial_somatic_truth 3 |
| COLO829_Illumina_SNV_d46_7fc7396 | d46, 4 GPUs | 268 | 19 | 18 | 15 | representative_allele_is_germline_truth 191, confident_no_truth_allele 45, near_truth_allele_mismatch 11, residual_partial_somatic_truth 3 |
| COLO829_Illumina_SNV_d46_7fc7396_s1 | d46, 4 GPUs | 725 | 494 | 485 | 483 | representative_allele_is_germline_truth 183, confident_no_truth_allele 47, near_truth_allele_mismatch 7, residual_partial_somatic_truth 3 |
| COLO829_Illumina_SNV_d46_7fc7396_s2 | d46, 4 GPUs | 386 | 97 | 99 | 93 | representative_allele_is_germline_truth 229, confident_no_truth_allele 46, near_truth_allele_mismatch 8, residual_partial_somatic_truth 4 |

## 2. The window false calls of the 4-GPU d9 runs on the d46 set

- distinct false-call alleles of the three 4-GPU d9 seeds in the windows: 3,060 (in all three seeds: 235)
- with a d46 tensor of the same GRCh38 allele: 2,444 (80%)
- called by COLO829_Illumina_SNV_d46_7fc7396 (p_somatic >= 0.640): 18 of the 2,444
- called by COLO829_Illumina_SNV_d46_7fc7396_s1 (p_somatic >= 0.703): 312 of the 2,444
- called by COLO829_Illumina_SNV_d46_7fc7396_s2 (p_somatic >= 0.303): 62 of the 2,444
