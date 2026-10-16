# Existing whole-GAM `vg stats -a` reports (cross-check of the sampled numbers)

| set | arm | records | unmapped % | perfect % of all | gapless % of all | soft-clip events per 100 mapped reads | mean MAPQ | file |
|---|---|---:|---:|---:|---:|---:|---:|---|
| HG008T_PacBio | d9 | 14,014,144 | 2.64 | 1.39 | 2.40 | 2.13 | 58.19 | `/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008_PacBio_gam_stats.txt` |
| HG008T_ONT | d9 | 3,642,157 | 7.77 | 2.56 | 5.13 | 36.49 | 57.69 | `/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_Northeastern-ONT-UL-20241216/HG008_ONT_gam_stats.txt` |
| HG008T_Illumina | d9 | 2,974,694,440 | 1.42 | 76.86 | 97.67 | 2.92 | 55.14 | `/scratch/jshen/data/HG008_GIAB/vg_giraffe_GAM/Liss_lab_BCM_Illumina-WGS_20240313/HG008-T_p23_BCM_IlluminaWGS_20240313.vg.stats` |
| HG008T_Illumina | linear | 2,974,694,440 | 1.51 | 67.54 | 95.56 | 3.17 | 55.32 | `/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_BCM_Illumina-WGS_20240313/vg_linear_hg38/HG008T_BCM_IlluminaWGS_hg38_linear_vg.gam.stats.txt` |
| COLO829T_fiberseq | linear | 27,846,735 | 5.62 | 0.17 | 0.65 | 4.78 | 57.73 | `/scratch/jshen/data/COLO829T/COLO829T_fiberseq/vg_linear_hg38/SMHTCOLO829T_fiberseq_hg38_linear_vg.gam.stats.txt` |
| COLO829T_ONT | d9 | 14,448,438 | 13.23 | 0.01 | 6.24 | 173.44 | 57.01 | `/scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.vg_stats.txt` |
| COLO829T_ONT | linear | 14,448,438 | 13.53 | 0.01 | 5.83 | 173.45 | 57.26 | `/scratch/jshen/data/COLO829T/COLO829T_ONT/vg_linear_hg38/COLO829T_ONT_hg38_linear_vg.gam.stats.txt` |
