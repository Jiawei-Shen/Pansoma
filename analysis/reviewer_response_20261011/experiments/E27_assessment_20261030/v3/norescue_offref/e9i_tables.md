# E9i: matched-normal reads at COLO829T graph-coordinate SNV calls (chr1)

ALT / REF 21-mer of each call counted in COLO829BL ONT reads (D001 + D003) and COLO829T ONT reads over the call's GRCh38 locus +- 500 bp. 'absent' = 0 normal ALT reads with >= 5 REF reads; 'somatic candidate' = absent from the normal and ALT in >= 2 tumor reads.

| run | group | calls | ALT in normal | one normal read | absent from normal | low coverage / no k-mer | somatic candidates |
|---|---|---:|---|---|---|---|---|
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024__norescue | unplaced_pass | 53 | 45 (84.9%) | 3 (5.7%) | 5 (9.4%) | 0 (0.0%) | 3 (5.7%) |
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024__norescue | ctrl_tp | 150 | 4 (2.7%) | 9 (6.0%) | 137 (91.3%) | 0 (0.0%) | 137 (91.3%) |
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024__norescue | ctrl_fp | 150 | 130 (86.7%) | 6 (4.0%) | 13 (8.7%) | 1 (0.7%) | 4 (2.7%) |
| COLO829_Illumina_SNV_d46_7fc7396__norescue | unplaced_pass | 18 | 10 (55.6%) | 5 (27.8%) | 3 (16.7%) | 0 (0.0%) | 3 (16.7%) |
| COLO829_Illumina_SNV_d46_7fc7396__norescue | ctrl_tp | 150 | 6 (4.0%) | 5 (3.3%) | 139 (92.7%) | 0 (0.0%) | 139 (92.7%) |
| COLO829_Illumina_SNV_d46_7fc7396__norescue | ctrl_fp | 150 | 126 (84.0%) | 2 (1.3%) | 21 (14.0%) | 1 (0.7%) | 18 (12.0%) |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024__norescue | unplaced_pass | 152 | 85 (55.9%) | 34 (22.4%) | 30 (19.7%) | 3 (2.0%) | 6 (3.9%) |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024__norescue | ctrl_tp | 150 | 5 (3.3%) | 7 (4.7%) | 138 (92.0%) | 0 (0.0%) | 138 (92.0%) |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024__norescue | ctrl_fp | 150 | 107 (71.3%) | 2 (1.3%) | 41 (27.3%) | 0 (0.0%) | 38 (25.3%) |
| COLO829_ONT_SNV_d46_7fc7396__norescue | unplaced_pass | 98 | 46 (46.9%) | 15 (15.3%) | 17 (17.3%) | 20 (20.4%) | 6 (6.1%) |
| COLO829_ONT_SNV_d46_7fc7396__norescue | ctrl_tp | 150 | 2 (1.3%) | 9 (6.0%) | 139 (92.7%) | 0 (0.0%) | 139 (92.7%) |
| COLO829_ONT_SNV_d46_7fc7396__norescue | ctrl_fp | 150 | 88 (58.7%) | 2 (1.3%) | 59 (39.3%) | 1 (0.7%) | 24 (16.0%) |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024__norescue | unplaced_pass | 35 | 31 (88.6%) | 2 (5.7%) | 2 (5.7%) | 0 (0.0%) | 2 (5.7%) |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024__norescue | ctrl_tp | 150 | 6 (4.0%) | 10 (6.7%) | 133 (88.7%) | 1 (0.7%) | 133 (88.7%) |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024__norescue | ctrl_fp | 150 | 117 (78.0%) | 17 (11.3%) | 15 (10.0%) | 1 (0.7%) | 12 (8.0%) |
| COLO829_fiberseq_SNV_d46_7fc7396__norescue | unplaced_pass | 78 | 59 (75.6%) | 10 (12.8%) | 9 (11.5%) | 0 (0.0%) | 2 (2.6%) |
| COLO829_fiberseq_SNV_d46_7fc7396__norescue | ctrl_tp | 150 | 6 (4.0%) | 3 (2.0%) | 141 (94.0%) | 0 (0.0%) | 141 (94.0%) |
| COLO829_fiberseq_SNV_d46_7fc7396__norescue | ctrl_fp | 150 | 131 (87.3%) | 1 (0.7%) | 17 (11.3%) | 1 (0.7%) | 17 (11.3%) |

## Unplaced PASS calls by E9c category

| run | category | calls | ALT in normal | absent from normal | somatic candidates |
|---|---|---:|---:|---:|---:|
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024__norescue | germline | 32 | 29 | 1 | 0 |
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024__norescue | population | 17 | 16 | 0 | 0 |
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024__norescue | tumor_only | 3 | 0 | 3 | 3 |
| COLO829_Illumina_SNV_7fc7396_ddp4_b1024__norescue | unsupported | 1 | 0 | 1 | 0 |
| COLO829_Illumina_SNV_d46_7fc7396__norescue | germline | 13 | 9 | 0 | 0 |
| COLO829_Illumina_SNV_d46_7fc7396__norescue | population | 1 | 0 | 0 | 0 |
| COLO829_Illumina_SNV_d46_7fc7396__norescue | tumor_only | 4 | 1 | 3 | 3 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024__norescue | germline | 94 | 38 | 21 | 2 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024__norescue | population | 54 | 47 | 5 | 0 |
| COLO829_ONT_SNV_7fc7396_ddp4_b1024__norescue | tumor_only | 4 | 0 | 4 | 4 |
| COLO829_ONT_SNV_d46_7fc7396__norescue | germline | 61 | 26 | 7 | 0 |
| COLO829_ONT_SNV_d46_7fc7396__norescue | population | 33 | 19 | 7 | 3 |
| COLO829_ONT_SNV_d46_7fc7396__norescue | tumor_only | 4 | 1 | 3 | 3 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024__norescue | germline | 15 | 13 | 0 | 0 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024__norescue | population | 18 | 18 | 0 | 0 |
| COLO829_fiberseq_SNV_7fc7396_ddp4_b1024__norescue | tumor_only | 2 | 0 | 2 | 2 |
| COLO829_fiberseq_SNV_d46_7fc7396__norescue | germline | 48 | 36 | 5 | 0 |
| COLO829_fiberseq_SNV_d46_7fc7396__norescue | population | 28 | 23 | 2 | 0 |
| COLO829_fiberseq_SNV_d46_7fc7396__norescue | tumor_only | 1 | 0 | 1 | 1 |
| COLO829_fiberseq_SNV_d46_7fc7396__norescue | unsupported | 1 | 0 | 1 | 1 |
