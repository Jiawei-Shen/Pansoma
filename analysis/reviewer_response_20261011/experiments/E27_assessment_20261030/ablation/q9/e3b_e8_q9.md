# q9: COLO829T Illumina SNV d9, E3b complete-output F1 and E8 VAF split (chr1, eval_bed, PASS after PoN)

## E3b complete-output F1 (every graph-coordinate PASS call counted as FP)

| label | run | placed F1 | ON: graph-coord PASS / complete F1 | OFF: graph-coord PASS / complete F1 |
|---|---|---|---|---|
| epoch 8, norm train, s0 | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024` | 0.5719 | 186 / 0.5546 | 53 / 0.5668 |
| epoch 6, normbq, s0 | `COLO829_Illumina_SNV_d9_normbq_sweep_e06` | 0.7734 | 693 / 0.6715 | 38 / 0.7671 |
| d46 s0 (reference) | `COLO829_Illumina_SNV_d46_7fc7396` | 0.7832 | 516 / 0.7019 | 18 / 0.7801 |
| epoch 8, norm train, s1 | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s1` | 0.6259 | 574 / 0.5668 | - (no __norescue run) |
| epoch 8, norm train, s2 | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s2` | 0.5664 | 766 / 0.5029 | - (no __norescue run) |
| epoch 8, normbq, s0 | `COLO829_Illumina_SNV_d9_normbq` | 0.5658 | 191 / 0.5484 | - (no __norescue run) |
| epoch 8, normbq, s1 | `COLO829_Illumina_SNV_d9_normbq_s1` | 0.6277 | 564 / 0.5693 | - (no __norescue run) |
| epoch 8, normbq, s2 | `COLO829_Illumina_SNV_d9_normbq_s2` | 0.5684 | 761 / 0.5048 | - (no __norescue run) |
| epoch 5, normbq, s1 | `COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05` | 0.7052 | 357 / 0.6571 | - (no __norescue run) |
| epoch 6, normbq, s2 | `COLO829_Illumina_SNV_d9_normbq_s2_sweep_e06` | 0.7378 | 474 / 0.6692 | - (no __norescue run) |

## E8 VAF split (truth INFO VAF_Ill; recall below / at or above 0.08)

| label | run | recall VAF < 0.08 (n) | recall VAF >= 0.08 (n) |
|---|---|---|---|
| epoch 8, norm train, s0 | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024` | 0.067 (405) | 0.823 (2038) |
| epoch 6, normbq, s0 | `COLO829_Illumina_SNV_d9_normbq_sweep_e06` | 0.067 (405) | 0.853 (2038) |
| epoch 8, norm train, s1 | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s1` | 0.067 (405) | 0.832 (2038) |
| epoch 8, norm train, s2 | `COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s2` | 0.067 (405) | 0.829 (2038) |
| epoch 8, normbq, s0 | `COLO829_Illumina_SNV_d9_normbq` | 0.067 (405) | 0.823 (2038) |
| epoch 8, normbq, s1 | `COLO829_Illumina_SNV_d9_normbq_s1` | 0.067 (405) | 0.833 (2038) |
| epoch 8, normbq, s2 | `COLO829_Illumina_SNV_d9_normbq_s2` | 0.067 (405) | 0.829 (2038) |
| epoch 5, normbq, s1 | `COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05` | 0.069 (405) | 0.831 (2038) |
| epoch 6, normbq, s2 | `COLO829_Illumina_SNV_d9_normbq_s2_sweep_e06` | 0.069 (405) | 0.823 (2038) |
