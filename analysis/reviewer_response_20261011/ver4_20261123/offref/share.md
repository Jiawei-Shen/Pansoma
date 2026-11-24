# Graph-coordinate share of all PASS SNV calls (six-tumor current models, HPRC v2.1, chr1, default output)

Whole chr1 (no benchmark BED), after the GRCh38 population filter. Definitions in share.py.

| Tumor | Platform | Graph-coordinate PASS (U) | Placed PASS, rtg TP + FP (no BED, after filter) | share U / (TP + FP + U) | Placed PASS records (after filter) | share U / (records + U) | Placed PASS records (before filter) | share before filter |
|---|---|---|---|---|---|---|---|---|
| HG008T | Illumina | 77 | 1,188 | 6.1% | 1,188 | 6.1% | 1,201 | 6.0% |
| HG008T | PacBio HiFi | 143 | 1,112 | 11.4% | 1,112 | 11.4% | 1,117 | 11.3% |
| HG008T | ONT | 90 | 778 | 10.4% | 778 | 10.4% | 778 | 10.4% |
| COLO829T | Illumina | 110 | 2,384 | 4.4% | 2,384 | 4.4% | 2,406 | 4.4% |
| COLO829T | Fiber-seq | 79 | 2,117 | 3.6% | 2,117 | 3.6% | 2,123 | 3.6% |
| COLO829T | ONT | 111 | 2,196 | 4.8% | 2,196 | 4.8% | 2,208 | 4.8% |
| HCC1395T | Illumina | 76 | 4,084 | 1.8% | 4,084 | 1.8% | 4,110 | 1.8% |
| HCC1395T | PacBio HiFi | 106 | 3,733 | 2.8% | 3,733 | 2.8% | 3,750 | 2.7% |
| HCC1395T | ONT | 143 | 3,546 | 3.9% | 3,546 | 3.9% | 3,601 | 3.8% |

Ranges (min-max over the call sets):

- HG008T / COLO829T: U / (TP + FP + U) 3.6%-11.4%; U / (PASS records + U) 3.6%-11.4%
- HCC1395T: U / (TP + FP + U) 1.8%-3.9%; U / (PASS records + U) 1.8%-3.9%
- all nine: U / (TP + FP + U) 1.8%-11.4%; U / (PASS records + U) 1.8%-11.4%
- HG008T: U / (TP + FP + U) 6.1%-11.4%; U / (PASS records + U) 6.1%-11.4%
- COLO829T: U / (TP + FP + U) 3.6%-4.8%; U / (PASS records + U) 3.6%-4.8%
