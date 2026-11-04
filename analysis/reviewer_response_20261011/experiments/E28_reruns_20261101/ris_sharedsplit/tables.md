# E28 (3) on RIS: R3.1 SNV, validation chr20-22 shared by the three references (seed 0)

(`collect.py`; chr1, after the PoN, in the benchmark BED. PASS = the checkpoint's call rule (validation recall 0.9); best = the chr1 curve's best point. Recipe of Table R3.1: 50.4 M, 20 epochs, pow:0.75, lr 1e-4, 4 x H100 x 1024, `--norm train`, `--val-chroms chr20 chr21 chr22`, training chr2-19.)

| set | reference | best epoch | PASS P / R / F1 | best F1 | run |
|---|---|---|---|---|---|
| HG008T Illumina | GRCh38-only | 4 | 0.304 / 0.824 / **0.444** | 0.522 | `HG008_Illumina_SNV_linear_valchr20-22` |
| HG008T Illumina | HPRC v1.1 (d9) | 6 | 0.452 / 0.828 / **0.584** | 0.611 | `HG008_Illumina_SNV_d9_valchr20-22` |
| HG008T Illumina | HPRC v2.1 (d46) | 5 | 0.448 / 0.841 / **0.585** | 0.607 | `HG008_Illumina_SNV_d46_valchr20-22` |
| HG008T PacBio | GRCh38-only | 5 | 0.452 / 0.845 / **0.589** | 0.602 | `HG008_PacBio_SNV_linear_valchr20-22` |
| HG008T PacBio | HPRC v1.1 (d9) | 8 | 0.477 / 0.831 / **0.606** | 0.632 | `HG008_PacBio_SNV_d9_valchr20-22` |
| HG008T PacBio | HPRC v2.1 (d46) | 9 | 0.472 / 0.834 / **0.603** | 0.625 | `HG008_PacBio_SNV_d46_valchr20-22` |
| HG008T ONT | GRCh38-only | 5 | 0.355 / 0.845 / **0.500** | 0.599 | `HG008_ONT_SNV_linear_valchr20-22` |
| HG008T ONT | HPRC v1.1 (d9) | 10 | 0.454 / 0.842 / **0.590** | 0.595 | `HG008_ONT_SNV_d9_valchr20-22` |
| HG008T ONT | HPRC v2.1 (d46) | 7 | 0.479 / 0.834 / **0.608** | 0.619 | `HG008_ONT_SNV_d46_valchr20-22` |
| COLO829T Illumina | GRCh38-only | 7 | 0.209 / 0.557 / **0.304** | 0.488 | `COLO829_Illumina_SNV_linear_valchr20-22` |
| COLO829T Illumina | HPRC v1.1 (d9) | 6 | 0.724 / 0.728 / **0.726** | 0.758 | `COLO829_Illumina_SNV_d9_valchr20-22` |
| COLO829T Illumina | HPRC v2.1 (d46) | 7 | 0.852 / 0.730 / **0.786** | 0.790 | `COLO829_Illumina_SNV_d46_valchr20-22` |
| COLO829T fiberseq | GRCh38-only | 8 | 0.908 / 0.549 / **0.684** | 0.759 | `COLO829_fiberseq_SNV_linear_valchr20-22` |
| COLO829T fiberseq | HPRC v1.1 (d9) | 10 | 0.871 / 0.701 / **0.777** | 0.786 | `COLO829_fiberseq_SNV_d9_valchr20-22` |
| COLO829T fiberseq | HPRC v2.1 (d46) | 11 | 0.865 / 0.722 / **0.787** | 0.792 | `COLO829_fiberseq_SNV_d46_valchr20-22` |
| COLO829T ONT | GRCh38-only | 6 | 0.684 / 0.533 / **0.599** | 0.601 | `COLO829_ONT_SNV_linear_valchr20-22` |
| COLO829T ONT | HPRC v1.1 (d9) | 7 | 0.786 / 0.704 / **0.743** | 0.747 | `COLO829_ONT_SNV_d9_valchr20-22` |
| COLO829T ONT | HPRC v2.1 (d46) | 7 | 0.846 / 0.664 / **0.744** | 0.756 | `COLO829_ONT_SNV_d46_valchr20-22` |
