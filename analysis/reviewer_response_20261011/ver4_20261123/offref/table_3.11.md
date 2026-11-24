# Response Table 3.11 (ver4): CHM13 representation of graph-coordinate SNV PASS calls

Six-tumor current SNV models (snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, seed 0), HPRC v2.1, chr1, default output. Definitions in table_3_11.py.

| **Tumor** | **Platform** | **Calls** | **CHM13 allele** |
|---|---|---|---|
| HG008T | Illumina | 77 | 50 (64.9%) |
| HG008T | PacBio HiFi | 143 | 113 (79.0%) |
| HG008T | ONT | 90 | 83 (92.2%) |
| COLO829T | Illumina | 110 | 53 (48.2%) |
| COLO829T | Fiber-seq | 79 | 46 (58.2%) |
| COLO829T | ONT | 111 | 76 (68.5%) |
| HCC1395T | Illumina | 76 | 54 (71.1%) |
| HCC1395T | PacBio HiFi | 106 | 95 (89.6%) |
| HCC1395T | ONT | 143 | 120 (83.9%) |
| **Total** | | **935** | **690 (73.8%)** |
| Total excluding HCC1395T | | 610 | 421 (69.0%) |

Per tumor: HG008T 246 (79.4%) of 310; COLO829T 175 (58.3%) of 300; HCC1395T 269 (82.8%) of 325.
Other calls (all nine): CHM13 anchor interval only 171, no CHM13 position 74; GRCh38 allele through the bubble's flanking nodes (E9a class i) 60; GRCh38 (class i) or CHM13 allele 709 (75.8%).
CHM13 VCF check (PASS calls placed on CHM13 but not on GRCh38 = CHM13 allele by projection): HG008T Illumina 50 / 50, HG008T PacBio HiFi 113 / 113, HG008T ONT 83 / 83, COLO829T Illumina 53 / 53, COLO829T Fiber-seq 46 / 46, COLO829T ONT 76 / 76, HCC1395T Illumina 54 / 54, HCC1395T PacBio HiFi 95 / 95, HCC1395T ONT 120 / 120.
