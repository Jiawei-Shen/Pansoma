# E9i: matched-normal reads at COLO829T graph-coordinate SNV calls (chr1)

ALT / REF 21-mer of each call counted in COLO829BL ONT reads (D001 + D003) and COLO829T ONT reads over the call's GRCh38 locus +- 500 bp. 'absent' = 0 normal ALT reads with >= 5 REF reads; 'somatic candidate' = absent from the normal and ALT in >= 2 tumor reads.

| run | group | calls | ALT in normal | one normal read | absent from normal | low coverage / no k-mer | somatic candidates |
|---|---|---:|---|---|---|---|---|
| COLO829_Illumina_SNV_d46_7fc7396 | unplaced_pass | 516 | 440 (85.3%) | 29 (5.6%) | 35 (6.8%) | 12 (2.3%) | 12 (2.3%) |
| COLO829_Illumina_SNV_d46_7fc7396 | ctrl_tp | 150 | 6 (4.0%) | 5 (3.3%) | 139 (92.7%) | 0 (0.0%) | 139 (92.7%) |
| COLO829_Illumina_SNV_d46_7fc7396 | ctrl_fp | 150 | 126 (84.0%) | 2 (1.3%) | 21 (14.0%) | 1 (0.7%) | 18 (12.0%) |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | unplaced_pass | 693 | 592 (85.4%) | 37 (5.3%) | 43 (6.2%) | 21 (3.0%) | 10 (1.4%) |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | ctrl_tp | 150 | 4 (2.7%) | 6 (4.0%) | 140 (93.3%) | 0 (0.0%) | 139 (92.7%) |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | ctrl_fp | 150 | 121 (80.7%) | 4 (2.7%) | 23 (15.3%) | 2 (1.3%) | 13 (8.7%) |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | unplaced_pass | 3000 | 2334 (77.8%) | 381 (12.7%) | 165 (5.5%) | 120 (4.0%) | 24 (0.8%) |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | ctrl_tp | 150 | 5 (3.3%) | 7 (4.7%) | 138 (92.0%) | 0 (0.0%) | 138 (92.0%) |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | ctrl_fp | 150 | 107 (71.3%) | 2 (1.3%) | 41 (27.3%) | 0 (0.0%) | 38 (25.3%) |
| COLO829_ONT_SNV_d46_7fc7396 | unplaced_pass | 3000 | 1981 (66.0%) | 341 (11.4%) | 480 (16.0%) | 198 (6.6%) | 40 (1.3%) |
| COLO829_ONT_SNV_d46_7fc7396 | ctrl_tp | 150 | 2 (1.3%) | 9 (6.0%) | 139 (92.7%) | 0 (0.0%) | 139 (92.7%) |
| COLO829_ONT_SNV_d46_7fc7396 | ctrl_fp | 150 | 88 (58.7%) | 2 (1.3%) | 59 (39.3%) | 1 (0.7%) | 24 (16.0%) |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | unplaced_pass | 1483 | 1376 (92.8%) | 37 (2.5%) | 33 (2.2%) | 37 (2.5%) | 10 (0.7%) |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | ctrl_tp | 150 | 6 (4.0%) | 10 (6.7%) | 133 (88.7%) | 1 (0.7%) | 133 (88.7%) |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | ctrl_fp | 150 | 117 (78.0%) | 17 (11.3%) | 15 (10.0%) | 1 (0.7%) | 12 (8.0%) |
| COLO829_fiberseq_SNV_d46_7fc7396 | unplaced_pass | 257 | 225 (87.5%) | 21 (8.2%) | 9 (3.5%) | 2 (0.8%) | 8 (3.1%) |
| COLO829_fiberseq_SNV_d46_7fc7396 | ctrl_tp | 150 | 6 (4.0%) | 3 (2.0%) | 141 (94.0%) | 0 (0.0%) | 141 (94.0%) |
| COLO829_fiberseq_SNV_d46_7fc7396 | ctrl_fp | 150 | 131 (87.3%) | 1 (0.7%) | 17 (11.3%) | 1 (0.7%) | 17 (11.3%) |

## Unplaced PASS calls by E9c category

| run | category | calls | ALT in normal | absent from normal | somatic candidates |
|---|---|---:|---:|---:|---:|
| COLO829_Illumina_SNV_d46_7fc7396 | germline | 222 | 177 | 19 | 1 |
| COLO829_Illumina_SNV_d46_7fc7396 | population | 276 | 258 | 4 | 1 |
| COLO829_Illumina_SNV_d46_7fc7396 | tumor_only | 17 | 5 | 11 | 10 |
| COLO829_Illumina_SNV_d46_7fc7396 | unsupported | 1 | 0 | 1 | 0 |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | germline | 208 | 143 | 30 | 1 |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | population | 469 | 446 | 0 | 0 |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | tumor_only | 12 | 0 | 12 | 9 |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | unsupported | 3 | 2 | 1 | 0 |
| COLO829_Illumina_SNV_d9_normbq_sweep_e06 | unverifiable | 1 | 1 | 0 | 0 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | germline | 904 | 622 | 76 | 8 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | population | 1945 | 1588 | 73 | 3 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | tumor_only | 10 | 0 | 10 | 10 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | unsupported | 55 | 45 | 6 | 3 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024 | unverifiable | 86 | 79 | 0 | 0 |
| COLO829_ONT_SNV_d46_7fc7396 | germline | 1361 | 607 | 410 | 10 |
| COLO829_ONT_SNV_d46_7fc7396 | population | 1495 | 1260 | 53 | 17 |
| COLO829_ONT_SNV_d46_7fc7396 | tumor_only | 12 | 1 | 11 | 11 |
| COLO829_ONT_SNV_d46_7fc7396 | unsupported | 60 | 46 | 6 | 2 |
| COLO829_ONT_SNV_d46_7fc7396 | unverifiable | 72 | 67 | 0 | 0 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | germline | 177 | 126 | 18 | 1 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | population | 1193 | 1157 | 5 | 0 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | tumor_only | 6 | 1 | 5 | 5 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | unsupported | 31 | 22 | 5 | 4 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024 | unverifiable | 76 | 70 | 0 | 0 |
| COLO829_fiberseq_SNV_d46_7fc7396 | germline | 27 | 15 | 0 | 0 |
| COLO829_fiberseq_SNV_d46_7fc7396 | population | 200 | 187 | 4 | 3 |
| COLO829_fiberseq_SNV_d46_7fc7396 | tumor_only | 7 | 1 | 5 | 5 |
| COLO829_fiberseq_SNV_d46_7fc7396 | unsupported | 9 | 9 | 0 | 0 |
| COLO829_fiberseq_SNV_d46_7fc7396 | unverifiable | 14 | 13 | 0 | 0 |
