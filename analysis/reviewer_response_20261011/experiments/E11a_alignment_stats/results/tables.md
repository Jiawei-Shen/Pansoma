# E11a tables: the same reads on four references

Written by compare.py; definitions in README.md. Percentages: perfect of all reads in the row; mismatch, indel, soft clip and MAPQ bins of mapped reads.

## HG008T_PacBio

### G. Genome-wide name-hash sample (crc32(name) % 100 < 10), graph arms

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 14,014,144 | 389,438 | 2.78 | 0.04 | 1.86 | 95.88 | 1,402,384 |
| HPRC v1.1 d9 | 14,014,144 | 369,434 | 2.64 | 0.07 | 2.15 | 95.88 | 1,402,384 |
| HPRC v2.1 d46 | 14,014,144 | 1,592 | 0.01 | 0.07 | 4.80 | 83.32 | 1,402,384 |

Sampled reads (union over graph arms): 1,402,384; absent from an arm (counted unmapped there): GRCh38-only graph 0, HPRC v1.1 d9 0, HPRC v2.1 d46 0.

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph (unmapped 2.79 %) | 1,402,384 | 0.15 | 94.53 | 99.43 | 2.65 | 0.04 | 2.19 | 0.72 | 1.16 | 95.89 | 0.9972 | 1.57 | 1.58 | 0.23 | - | - | - | 0.00 |
| HPRC v1.1 d9 (unmapped 2.64 %) | 1,402,384 | 1.39 | 67.29 | 97.54 | 2.00 | 0.07 | 2.44 | 0.72 | 0.89 | 95.89 | 0.9991 | 0.63 | 1.39 | 0.17 | 92.08 | 8.27 | 98.01 | 0.00 |
| HPRC v2.1 d46 (unmapped 0.01 %) | 1,402,384 | 1.42 | 66.40 | 97.61 | 3.75 | 0.07 | 5.47 | 0.96 | 10.15 | 83.35 | 0.9990 | 1.72 | 1.56 | 0.80 | 89.94 | 7.07 | 93.82 | 0.00 |

Same reads, B relative to A:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,402,384 | 1,362,330 | 1.34 | 1.23 | 0.20 | 0.43 | 0.15 | 1.39 | 1.23 | 0.00 | 2.57 | 1.94 | 98.88 | 99.98 | 307 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,402,384 | 1,363,213 | 1.05 | 12.59 | 0.12 | 1.03 | 0.15 | 1.42 | 1.26 | 0.00 | 2.57 | 3.75 | 98.57 | 99.98 | 291 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,402,384 | 1,365,262 | 1.03 | 12.37 | 0.11 | 0.87 | 1.39 | 1.42 | 0.16 | 0.13 | 1.94 | 3.75 | 98.68 | 99.97 | 454 |

Reads that became perfect / stopped being perfect:

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 17,267 | 95.99 | 11.0 | 15.60 | 84.32 | 0.49 | 0.03 | 9.0 | 99.72 | 1.94 | 0.38 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 17,721 | 96.16 | 12.0 | 15.71 | 84.19 | 0.47 | 0.05 | 9.0 | 99.77 | 0.45 | 32.43 | 3 | 66.67 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,221 | 97.07 | 15.0 | 47.59 | 52.32 | 0.05 | 0.09 | 1.0 | 99.46 | 0.54 | 29.70 | 1,770 | 89.77 |

By region (anchor: linear 1,363,265, d9 2,718, d46 35,097, none 1,304):

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 1,273,579 | GRCh38-only graph | 0.01 | 0.16 | 94.35 | 99.44 | 1.94 | 0.00 | 99.14 | 0.9973 | - |
| SMaHT easy | 1,273,579 | HPRC v1.1 d9 | 0.00 | 1.46 | 66.17 | 97.51 | 1.38 | 0.00 | 99.37 | 0.9991 | 93.41 |
| SMaHT easy | 1,273,579 | HPRC v2.1 d46 | 0.00 | 1.50 | 64.28 | 97.54 | 1.38 | 0.00 | 89.45 | 0.9991 | 93.07 |
| SMaHT difficult | 5,066 | GRCh38-only graph | 9.59 | 0.12 | 96.86 | 98.58 | 7.51 | 0.00 | 94.85 | 0.9961 | - |
| SMaHT difficult | 5,066 | HPRC v1.1 d9 | 8.19 | 1.58 | 74.37 | 96.13 | 2.82 | 0.02 | 97.76 | 0.9990 | 95.38 |
| SMaHT difficult | 5,066 | HPRC v2.1 d46 | 0.00 | 1.76 | 74.40 | 96.09 | 10.74 | 0.00 | 77.48 | 0.9988 | 95.50 |
| SMaHT extreme | 109,232 | GRCh38-only graph | 29.30 | 0.07 | 97.02 | 99.48 | 10.69 | 0.66 | 48.75 | 0.9949 | - |
| SMaHT extreme | 109,232 | HPRC v1.1 d9 | 28.15 | 0.63 | 81.85 | 97.99 | 8.32 | 1.18 | 46.31 | 0.9982 | 76.48 |
| SMaHT extreme | 109,232 | HPRC v2.1 d46 | 0.03 | 0.62 | 86.54 | 98.31 | 23.55 | 0.85 | 22.52 | 0.9948 | 61.38 |
| segdup | 79,263 | GRCh38-only graph | 10.82 | 0.09 | 96.29 | 99.52 | 6.80 | 0.26 | 61.91 | 0.9962 | - |
| segdup | 79,263 | HPRC v1.1 d9 | 10.80 | 0.88 | 76.40 | 97.95 | 4.62 | 0.84 | 59.65 | 0.9988 | 91.24 |
| segdup | 79,263 | HPRC v2.1 d46 | 0.01 | 0.88 | 78.20 | 98.20 | 13.90 | 0.10 | 40.34 | 0.9985 | 85.83 |
| not segdup | 1,308,614 | GRCh38-only graph | 1.83 | 0.16 | 94.41 | 99.43 | 2.21 | 0.03 | 98.15 | 0.9973 | - |
| not segdup | 1,308,614 | HPRC v1.1 d9 | 1.73 | 1.43 | 66.60 | 97.51 | 1.63 | 0.03 | 98.31 | 0.9991 | 92.51 |
| not segdup | 1,308,614 | HPRC v2.1 d46 | 0.00 | 1.47 | 65.34 | 97.56 | 2.51 | 0.06 | 86.79 | 0.9991 | 90.88 |
| in HG008T benchmark BED | 1,325,549 | GRCh38-only graph | 0.24 | 0.16 | 94.40 | 99.43 | 2.09 | 0.01 | 97.98 | 0.9973 | - |
| in HG008T benchmark BED | 1,325,549 | HPRC v1.1 d9 | 0.21 | 1.46 | 66.49 | 97.49 | 1.47 | 0.03 | 98.13 | 0.9991 | 93.31 |
| in HG008T benchmark BED | 1,325,549 | HPRC v2.1 d46 | 0.00 | 1.49 | 64.72 | 97.54 | 1.65 | 0.00 | 87.71 | 0.9991 | 92.91 |
| outside HG008T benchmark BED | 63,403 | GRCh38-only graph | 46.95 | 0.02 | 98.71 | 99.70 | 17.18 | 1.21 | 27.59 | 0.9875 | - |
| outside HG008T benchmark BED | 63,403 | HPRC v1.1 d9 | 45.32 | 0.19 | 91.38 | 98.90 | 14.48 | 1.36 | 25.42 | 0.9932 | 57.94 |
| outside HG008T benchmark BED | 63,403 | HPRC v2.1 d46 | 0.05 | 0.20 | 94.97 | 98.87 | 35.64 | 1.44 | 8.48 | 0.9600 | 40.61 |

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.03 | 0.04 | 1.30 | 0.00 | 99.98 | 4 | 96.05 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.03 | 0.09 | 1.34 | 0.00 | 99.98 | 3 | 96.20 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.01 | 0.07 | 0.17 | 0.13 | 99.98 | 8 | 97.12 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 0.90 | 0.37 | 1.46 | 0.00 | 99.58 | 6 | 97.30 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 0.90 | 1.09 | 1.64 | 0.00 | 99.28 | 11 | 97.59 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.26 | 1.55 | 0.20 | 0.02 | 99.08 | 5 | 100.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 3.05 | 6.39 | 0.55 | 0.00 | 82.74 | 249 | 94.20 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 1.45 | 13.97 | 0.55 | 0.00 | 78.56 | 219 | 94.64 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.47 | 11.54 | 0.08 | 0.09 | 80.75 | 239 | 95.65 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 2.22 | 5.35 | 0.79 | 0.00 | 88.08 | 123 | 94.75 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1.27 | 11.76 | 0.79 | 0.00 | 86.21 | 139 | 94.92 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.37 | 8.73 | 0.11 | 0.11 | 89.39 | 160 | 95.45 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 0.09 | 0.13 | 1.27 | 0.00 | 99.61 | 136 | 96.04 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 0.05 | 0.29 | 1.31 | 0.00 | 99.45 | 94 | 96.20 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.03 | 0.30 | 0.16 | 0.13 | 99.38 | 92 | 97.14 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 0.09 | 0.21 | 1.29 | 0.00 | 99.70 | 85 | 96.01 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 0.07 | 0.52 | 1.33 | 0.00 | 99.65 | 97 | 96.17 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.05 | 0.40 | 0.17 | 0.13 | 99.69 | 111 | 97.13 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 4.47 | 8.20 | 0.17 | 0.00 | 70.96 | 181 | 93.58 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 1.67 | 15.77 | 0.19 | 0.00 | 63.42 | 139 | 93.22 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.84 | 14.12 | 0.04 | 0.03 | 66.84 | 141 | 92.31 |

### W. Window sample (148 GRCh38 windows, flank 25,000 bp; reads placed in a window by at least one arm: 34,332; random 5,984, random_extra 26,125, segdup 2,223)

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 97.88 | 2.12 | 0.00 |
| GRCh38-only graph | genome | 96.23 | 3.35 | 0.43 |
| HPRC v1.1 d9 | genome | 97.04 | 2.66 | 0.30 |
| HPRC v2.1 d46 | genome | 96.76 | 3.23 | 0.00 |

Unmapped records of the whole file (GAM head / BAM unplaced): GRCh38 BAM 1,055

Reads placed in the window by the arm:

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 33,604 | 0.18 | 95.01 | 99.44 | 6.98 | 0.03 | 3.03 | 0.60 | 0.59 | 95.75 | 0.9971 | 1.55 | 1.60 | 0.71 | - | - | - | 0.00 |
| GRCh38-only graph | 33,036 | 0.18 | 95.02 | 99.44 | 2.27 | 0.06 | 2.46 | 0.54 | 1.08 | 95.87 | 0.9972 | 1.36 | 1.58 | 0.16 | - | - | - | 0.00 |
| HPRC v1.1 d9 | 33,316 | 1.52 | 67.20 | 97.51 | 1.61 | 0.23 | 3.28 | 0.43 | 0.72 | 95.33 | 0.9991 | 0.42 | 1.38 | 0.08 | 93.40 | 7.59 | 98.56 | 0.00 |
| HPRC v2.1 d46 | 33,221 | 1.54 | 65.21 | 97.57 | 1.83 | 0.00 | 3.56 | 0.81 | 10.48 | 85.14 | 0.9991 | 0.42 | 1.40 | 0.16 | 92.96 | 4.77 | 98.50 | 0.00 |

Same reads, both arms place them in the window (B relative to A); last columns: reads A places in the window that B places elsewhere / leaves unmapped:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both  A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 32,780 | 32,780 | 1.06 | 1.16 | 0.12 | 0.38 | 0.18 | 1.54 | 1.36 | 0.00 | 2.15 | 1.55 | 99.79 | 100.00 | 1 | 255 | 1 | 523 | 13 |
| HPRC v2.1 d46 vs GRCh38-only graph | 32,692 | 32,692 | 0.88 | 12.17 | 0.06 | 0.80 | 0.18 | 1.56 | 1.37 | 0.00 | 2.17 | 1.54 | 99.75 | 100.00 | 1 | 344 | 0 | 429 | 100 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 32,897 | 32,897 | 1.20 | 11.79 | 0.09 | 0.60 | 1.54 | 1.55 | 0.12 | 0.12 | 1.52 | 1.52 | 99.77 | 100.00 | 1 | 419 | 0 | 247 | 77 |
| GRCh38-only graph vs GRCh38 BAM | 32,880 | 32,880 | 0.86 | 2.83 | 0.04 | 0.96 | 0.18 | 0.18 | 0.00 | 0.00 | 5.40 | 2.16 | 99.62 | 100.00 | 0 | 620 | 104 | 156 | 0 |
| HPRC v1.1 d9 vs GRCh38 BAM | 32,756 | 32,756 | 0.72 | 2.66 | 0.07 | 1.24 | 0.18 | 1.54 | 1.36 | 0.00 | 5.29 | 1.51 | 99.59 | 100.00 | 0 | 779 | 69 | 560 | 0 |
| HPRC v2.1 d46 vs GRCh38 BAM | 32,775 | 32,775 | 0.80 | 12.90 | 0.01 | 1.69 | 0.18 | 1.55 | 1.37 | 0.00 | 5.46 | 1.65 | 99.46 | 100.00 | 0 | 828 | 1 | 446 | 0 |

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 446 | 97.53 | 11.0 | 14.80 | 85.20 | 0.22 | 0.00 | 8.5 | 100.00 | 1.35 | 0.45 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 449 | 97.55 | 11.0 | 15.59 | 84.41 | 0.22 | 0.00 | 9.0 | 100.00 | 0.00 | 31.18 | 0 | n/a |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 41 | 97.56 | 13.5 | 36.59 | 63.41 | 0.00 | 0.00 | 0.0 | 100.00 | 0.00 | 24.39 | 38 | 86.84 |
| GRCh38-only graph vs GRCh38 BAM | 0 | | | | | | | | | | | | |
| HPRC v1.1 d9 vs GRCh38 BAM | 446 | 97.53 | 11.0 | 14.80 | 84.98 | 1.79 | 0.00 | 8.5 | 100.00 | 0.22 | 3.81 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38 BAM | 449 | 97.55 | 11.0 | 15.59 | 84.19 | 1.56 | 0.00 | 9.0 | 100.00 | 0.00 | 31.85 | 0 | n/a |

By window kind and region:

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 5,984 | GRCh38 BAM | 97.84 | 2.16 | 0.00 | 0.10 | 95.05 | 99.57 | 5.33 | 0.02 | 97.83 | 0.9971 |
| windows: random | 5,984 | GRCh38-only graph | 97.66 | 2.34 | 0.00 | 0.10 | 95.09 | 99.57 | 2.19 | 0.00 | 96.59 | 0.9971 |
| windows: random | 5,984 | HPRC v1.1 d9 | 99.08 | 0.92 | 0.00 | 1.33 | 70.16 | 97.81 | 1.62 | 0.88 | 96.04 | 0.9990 |
| windows: random | 5,984 | HPRC v2.1 d46 | 96.86 | 3.14 | 0.00 | 1.38 | 66.65 | 97.83 | 1.55 | 0.00 | 88.30 | 0.9991 |
| windows: random_extra | 26,125 | GRCh38 BAM | 99.18 | 0.82 | 0.00 | 0.20 | 94.96 | 99.41 | 6.95 | 0.02 | 96.40 | 0.9972 |
| windows: random_extra | 26,125 | GRCh38-only graph | 97.25 | 2.27 | 0.48 | 0.20 | 94.95 | 99.40 | 1.94 | 0.05 | 97.34 | 0.9972 |
| windows: random_extra | 26,125 | HPRC v1.1 d9 | 97.16 | 2.49 | 0.35 | 1.62 | 65.59 | 97.36 | 1.44 | 0.08 | 97.75 | 0.9991 |
| windows: random_extra | 26,125 | HPRC v2.1 d46 | 97.48 | 2.51 | 0.00 | 1.63 | 63.88 | 97.44 | 1.72 | 0.00 | 87.05 | 0.9991 |
| windows: segdup | 2,223 | GRCh38 BAM | 82.64 | 17.36 | 0.00 | 0.16 | 95.70 | 99.56 | 12.63 | 0.22 | 79.91 | 0.9960 |
| windows: segdup | 2,223 | GRCh38-only graph | 80.34 | 18.76 | 0.90 | 0.17 | 95.86 | 99.61 | 7.17 | 0.45 | 72.51 | 0.9961 |
| windows: segdup | 2,223 | HPRC v1.1 d9 | 90.19 | 9.31 | 0.49 | 0.80 | 78.85 | 98.40 | 3.74 | 0.15 | 62.59 | 0.9989 |
| windows: segdup | 2,223 | HPRC v2.1 d46 | 88.08 | 11.92 | 0.00 | 0.77 | 78.24 | 98.52 | 4.09 | 0.00 | 51.07 | 0.9988 |
| SMaHT easy | 30,016 | GRCh38 BAM | 99.95 | 0.05 | 0.00 | 0.18 | 94.78 | 99.42 | 4.66 | 0.00 | 100.00 | 0.9973 |
| SMaHT easy | 30,016 | GRCh38-only graph | 99.93 | 0.05 | 0.02 | 0.18 | 94.87 | 99.43 | 1.78 | 0.00 | 99.17 | 0.9973 |
| SMaHT easy | 30,016 | HPRC v1.1 d9 | 100.00 | 0.00 | 0.00 | 1.57 | 66.07 | 97.48 | 1.35 | 0.00 | 99.42 | 0.9991 |
| SMaHT easy | 30,016 | HPRC v2.1 d46 | 99.98 | 0.02 | 0.00 | 1.59 | 63.97 | 97.54 | 1.35 | 0.00 | 89.49 | 0.9991 |
| SMaHT difficult | 160 | GRCh38 BAM | 100.00 | 0.00 | 0.00 | 0.00 | 100.00 | 100.00 | 20.00 | 0.00 | 99.38 | 0.9937 |
| SMaHT difficult | 160 | GRCh38-only graph | 95.00 | 0.00 | 5.00 | 0.00 | 100.00 | 100.00 | 8.55 | 0.00 | 95.39 | 0.9940 |
| SMaHT difficult | 160 | HPRC v1.1 d9 | 95.62 | 0.62 | 3.75 | 1.96 | 77.12 | 96.73 | 1.96 | 0.00 | 99.35 | 0.9992 |
| SMaHT difficult | 160 | HPRC v2.1 d46 | 97.50 | 2.50 | 0.00 | 1.92 | 77.56 | 96.79 | 3.85 | 0.00 | 84.62 | 0.9992 |
| SMaHT extreme | 4,156 | GRCh38 BAM | 82.82 | 17.18 | 0.00 | 0.15 | 96.83 | 99.59 | 26.61 | 0.26 | 58.51 | 0.9944 |
| SMaHT extreme | 4,156 | GRCh38-only graph | 69.49 | 27.33 | 3.18 | 0.17 | 96.36 | 99.52 | 6.99 | 0.69 | 61.57 | 0.9957 |
| SMaHT extreme | 4,156 | HPRC v1.1 d9 | 75.72 | 21.94 | 2.33 | 1.02 | 77.47 | 97.81 | 4.04 | 2.41 | 56.15 | 0.9988 |
| SMaHT extreme | 4,156 | HPRC v2.1 d46 | 73.53 | 26.44 | 0.02 | 1.01 | 76.70 | 97.97 | 6.48 | 0.00 | 42.54 | 0.9987 |
| segdup | 3,507 | GRCh38 BAM | 84.46 | 15.54 | 0.00 | 0.20 | 95.88 | 99.43 | 10.77 | 0.30 | 75.22 | 0.9963 |
| segdup | 3,507 | GRCh38-only graph | 82.52 | 16.62 | 0.86 | 0.21 | 95.96 | 99.41 | 6.01 | 0.66 | 68.76 | 0.9963 |
| segdup | 3,507 | HPRC v1.1 d9 | 90.33 | 9.12 | 0.54 | 1.14 | 75.73 | 97.79 | 2.81 | 2.37 | 62.63 | 0.9989 |
| segdup | 3,507 | HPRC v2.1 d46 | 85.32 | 14.68 | 0.00 | 1.14 | 74.23 | 97.93 | 3.11 | 0.00 | 51.14 | 0.9989 |
| not segdup | 30,825 | GRCh38 BAM | 99.41 | 0.59 | 0.00 | 0.18 | 94.93 | 99.45 | 6.62 | 0.00 | 97.73 | 0.9972 |
| not segdup | 30,825 | GRCh38-only graph | 97.78 | 1.84 | 0.38 | 0.18 | 94.93 | 99.44 | 1.91 | 0.00 | 98.47 | 0.9972 |
| not segdup | 30,825 | HPRC v1.1 d9 | 97.80 | 1.92 | 0.27 | 1.56 | 66.30 | 97.48 | 1.48 | 0.00 | 98.77 | 0.9991 |
| not segdup | 30,825 | HPRC v2.1 d46 | 98.07 | 1.93 | 0.00 | 1.57 | 64.32 | 97.54 | 1.71 | 0.00 | 88.51 | 0.9991 |
| in HG008T benchmark BED | 32,658 | GRCh38 BAM | 99.21 | 0.79 | 0.00 | 0.19 | 94.89 | 99.43 | 5.03 | 0.01 | 98.77 | 0.9972 |
| in HG008T benchmark BED | 32,658 | GRCh38-only graph | 99.03 | 0.83 | 0.13 | 0.19 | 94.96 | 99.44 | 2.06 | 0.05 | 97.63 | 0.9972 |
| in HG008T benchmark BED | 32,658 | HPRC v1.1 d9 | 99.42 | 0.48 | 0.10 | 1.56 | 66.56 | 97.47 | 1.43 | 0.21 | 97.47 | 0.9991 |
| in HG008T benchmark BED | 32,658 | HPRC v2.1 d46 | 99.09 | 0.91 | 0.00 | 1.57 | 64.54 | 97.53 | 1.49 | 0.00 | 87.29 | 0.9991 |
| outside HG008T benchmark BED | 1,674 | GRCh38 BAM | 71.92 | 28.08 | 0.00 | 0.00 | 98.42 | 99.75 | 59.55 | 0.58 | 14.53 | 0.9744 |
| outside HG008T benchmark BED | 1,674 | GRCh38-only graph | 41.46 | 52.45 | 6.09 | 0.00 | 98.13 | 99.57 | 11.82 | 0.58 | 13.98 | 0.9895 |
| outside HG008T benchmark BED | 1,674 | HPRC v1.1 d9 | 50.66 | 45.10 | 4.24 | 0.12 | 91.75 | 98.82 | 8.37 | 0.83 | 13.56 | 0.9964 |
| outside HG008T benchmark BED | 1,674 | HPRC v2.1 d46 | 51.43 | 48.57 | 0.00 | 0.12 | 90.36 | 99.07 | 14.75 | 0.00 | 4.41 | 0.9955 |

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | HPRC v1.1 d9 vs GRCh38-only graph | 5,815 | 0.12 | 0.07 | 1.24 | 0.00 | 99.95 | 0.50 | 0.00 |
| windows: random | HPRC v2.1 d46 vs GRCh38-only graph | 5,765 | 0.03 | 0.14 | 1.27 | 0.00 | 100.00 | 1.35 | 0.00 |
| windows: random | HPRC v2.1 d46 vs HPRC v1.1 d9 | 5,781 | 0.02 | 0.14 | 0.14 | 0.12 | 99.91 | 2.50 | 0.00 |
| windows: random | GRCh38-only graph vs GRCh38 BAM | 5,826 | 0.00 | 0.43 | 0.00 | 0.00 | 99.86 | 0.50 | 0.00 |
| windows: random | HPRC v1.1 d9 vs GRCh38 BAM | 5,810 | 0.00 | 0.33 | 1.24 | 0.00 | 99.79 | 0.77 | 0.00 |
| windows: random | HPRC v2.1 d46 vs GRCh38 BAM | 5,762 | 0.00 | 0.50 | 1.27 | 0.00 | 99.81 | 1.59 | 0.00 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38-only graph | 25,270 | 0.09 | 0.14 | 1.43 | 0.00 | 99.77 | 0.53 | 0.00 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38-only graph | 25,247 | 0.03 | 0.43 | 1.44 | 0.00 | 99.73 | 0.63 | 0.00 |
| windows: random_extra | HPRC v2.1 d46 vs HPRC v1.1 d9 | 25,262 | 0.02 | 0.38 | 0.13 | 0.12 | 99.75 | 0.47 | 0.00 |
| windows: random_extra | GRCh38-only graph vs GRCh38 BAM | 25,323 | 0.04 | 0.35 | 0.00 | 0.00 | 99.57 | 1.93 | 0.35 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38 BAM | 25,262 | 0.08 | 0.41 | 1.43 | 0.00 | 99.54 | 2.27 | 0.24 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38 BAM | 25,329 | 0.01 | 0.71 | 1.44 | 0.00 | 99.38 | 2.25 | 0.00 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 1,695 | 0.47 | 5.13 | 0.77 | 0.00 | 99.47 | 5.10 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1,680 | 0.54 | 8.63 | 0.71 | 0.00 | 99.17 | 5.94 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,854 | 1.24 | 5.12 | 0.00 | 0.05 | 99.57 | 7.53 | 0.00 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 1,731 | 0.17 | 11.73 | 0.00 | 0.00 | 99.60 | 5.01 | 0.76 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 1,684 | 0.24 | 16.75 | 0.77 | 0.00 | 99.64 | 7.95 | 0.38 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 1,684 | 0.06 | 20.55 | 0.71 | 0.00 | 99.47 | 8.33 | 0.00 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 29,996 | 0.02 | 0.01 | 1.39 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 29,995 | 0.02 | 0.03 | 1.40 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 30,009 | 0.01 | 0.03 | 0.13 | 0.12 | 100.00 | 0.02 | 0.00 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 29,996 | 0.00 | 0.05 | 0.00 | 0.00 | 100.00 | 0.00 | 0.02 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 30,002 | 0.00 | 0.04 | 1.39 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 30,001 | 0.00 | 0.06 | 1.40 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 152 | 1.32 | 0.00 | 1.97 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 152 | 1.32 | 0.00 | 1.97 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 153 | 0.00 | 0.00 | 0.00 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | GRCh38-only graph vs GRCh38 BAM | 152 | 0.00 | 1.32 | 0.00 | 0.00 | 100.00 | 0.00 | 5.00 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38 BAM | 153 | 0.00 | 0.00 | 1.96 | 0.00 | 100.00 | 0.62 | 3.75 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38 BAM | 156 | 0.00 | 1.92 | 1.92 | 0.00 | 100.00 | 2.50 | 0.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 2,632 | 1.14 | 4.64 | 0.99 | 0.00 | 97.38 | 8.83 | 0.03 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 2,545 | 0.43 | 9.90 | 0.98 | 0.00 | 96.74 | 11.88 | 0.00 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,735 | 0.99 | 6.91 | 0.04 | 0.07 | 97.22 | 13.09 | 0.00 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 2,732 | 0.44 | 10.94 | 0.00 | 0.00 | 95.42 | 18.01 | 2.61 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 2,601 | 0.88 | 15.07 | 1.00 | 0.00 | 94.85 | 22.60 | 1.83 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 2,618 | 0.15 | 20.40 | 0.95 | 0.00 | 93.24 | 23.91 | 0.03 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 2,725 | 1.06 | 4.44 | 1.06 | 0.00 | 98.42 | 5.84 | 0.00 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 2,660 | 0.45 | 9.02 | 1.02 | 0.00 | 98.53 | 8.09 | 0.00 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,861 | 0.94 | 6.15 | 0.03 | 0.10 | 98.60 | 9.69 | 0.00 |
| segdup | GRCh38-only graph vs GRCh38 BAM | 2,807 | 0.14 | 10.55 | 0.00 | 0.00 | 99.36 | 4.46 | 0.78 |
| segdup | HPRC v1.1 d9 vs GRCh38 BAM | 2,700 | 0.52 | 14.37 | 1.07 | 0.00 | 98.30 | 8.37 | 0.47 |
| segdup | HPRC v2.1 d46 vs GRCh38 BAM | 2,656 | 0.11 | 19.80 | 1.02 | 0.00 | 98.38 | 10.33 | 0.00 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 30,055 | 0.03 | 0.02 | 1.39 | 0.00 | 99.91 | 0.29 | 0.00 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 30,032 | 0.02 | 0.07 | 1.41 | 0.00 | 99.85 | 0.36 | 0.00 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 30,036 | 0.01 | 0.07 | 0.13 | 0.12 | 99.88 | 0.37 | 0.00 |
| not segdup | GRCh38-only graph vs GRCh38 BAM | 30,073 | 0.03 | 0.07 | 0.00 | 0.00 | 99.64 | 1.59 | 0.26 |
| not segdup | HPRC v1.1 d9 vs GRCh38 BAM | 30,056 | 0.03 | 0.06 | 1.39 | 0.00 | 99.71 | 1.73 | 0.18 |
| not segdup | HPRC v2.1 d46 vs GRCh38 BAM | 30,119 | 0.00 | 0.10 | 1.40 | 0.00 | 99.56 | 1.70 | 0.00 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 32,254 | 0.06 | 0.33 | 1.38 | 0.00 | 99.98 | 0.27 | 0.00 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 32,221 | 0.05 | 0.68 | 1.39 | 0.00 | 99.98 | 0.37 | 0.00 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 32,275 | 0.07 | 0.44 | 0.13 | 0.12 | 99.97 | 0.59 | 0.00 |
| in HG008T benchmark BED | GRCh38-only graph vs GRCh38 BAM | 32,295 | 0.01 | 0.59 | 0.00 | 0.00 | 99.97 | 0.22 | 0.10 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38 BAM | 32,258 | 0.02 | 0.88 | 1.38 | 0.00 | 99.95 | 0.37 | 0.07 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38 BAM | 32,231 | 0.01 | 1.26 | 1.39 | 0.00 | 99.93 | 0.52 | 0.00 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 526 | 3.80 | 3.80 | 0.19 | 0.00 | 87.83 | 24.06 | 0.14 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 471 | 0.85 | 9.13 | 0.21 | 0.00 | 83.44 | 32.13 | 0.00 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 622 | 1.13 | 9.16 | 0.00 | 0.00 | 89.39 | 26.65 | 0.00 |
| outside HG008T benchmark BED | GRCh38-only graph vs GRCh38 BAM | 585 | 1.71 | 21.37 | 0.00 | 0.00 | 80.51 | 45.51 | 5.90 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38 BAM | 498 | 3.61 | 24.50 | 0.20 | 0.00 | 76.31 | 54.73 | 3.90 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38 BAM | 544 | 0.37 | 27.39 | 0.18 | 0.00 | 71.32 | 54.82 | 0.00 |

## HG008T_ONT

### G. Genome-wide name-hash sample (crc32(name) % 100 < 10), graph arms

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 3,642,157 | 287,106 | 7.88 | 0.93 | 2.26 | 95.18 | 364,076 |
| HPRC v1.1 d9 | 3,642,157 | 283,090 | 7.77 | 1.00 | 2.55 | 95.02 | 364,076 |
| HPRC v2.1 d46 | 3,642,157 | 292,931 | 8.04 | 0.14 | 1.71 | 94.16 | 364,076 |

Sampled reads (union over graph arms): 364,076; absent from an arm (counted unmapped there): GRCh38-only graph 0, HPRC v1.1 d9 0, HPRC v2.1 d46 0.

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph (unmapped 7.86 %) | 364,076 | 1.32 | 92.47 | 95.19 | 30.37 | 0.93 | 1.79 | 1.10 | 1.00 | 95.19 | 0.9941 | 5.22 | 4.27 | 1.76 | - | - | - | 0.00 |
| HPRC v1.1 d9 (unmapped 7.75 %) | 364,076 | 2.54 | 86.79 | 94.47 | 29.89 | 1.00 | 2.03 | 0.99 | 0.94 | 95.03 | 0.9964 | 4.26 | 4.12 | 1.64 | 81.48 | 8.82 | 96.42 | 0.00 |
| HPRC v2.1 d46 (unmapped 8.01 %) | 364,076 | 2.50 | 87.67 | 95.38 | 29.61 | 0.14 | 2.27 | 1.80 | 1.61 | 94.18 | 0.9963 | 4.35 | 4.33 | 1.71 | 81.27 | 6.88 | 93.46 | 0.00 |

Same reads, B relative to A:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 364,076 | 334,048 | 1.56 | 2.13 | 0.47 | 0.72 | 1.32 | 2.54 | 1.22 | 0.00 | 27.98 | 27.57 | 97.45 | 99.68 | 1,019 |
| HPRC v2.1 d46 vs GRCh38-only graph | 364,076 | 330,369 | 2.06 | 3.58 | 0.72 | 1.68 | 1.32 | 2.50 | 1.20 | 0.02 | 27.98 | 27.24 | 97.23 | 99.38 | 1,953 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 364,076 | 331,036 | 2.36 | 3.17 | 0.79 | 1.48 | 2.54 | 2.50 | 0.06 | 0.10 | 27.57 | 27.24 | 97.40 | 99.39 | 1,915 |

Reads that became perfect / stopped being perfect:

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 4,458 | 95.56 | 2.0 | 67.50 | 32.12 | 0.61 | 0.18 | 1.0 | 98.67 | 1.66 | 1.57 | 18 | 44.44 |
| HPRC v2.1 d46 vs GRCh38-only graph | 4,374 | 95.61 | 2.0 | 69.25 | 30.32 | 0.48 | 0.18 | 1.0 | 98.81 | 1.60 | 1.63 | 88 | 53.41 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 213 | 85.92 | 2.0 | 54.46 | 44.60 | 0.94 | 0.00 | 1.0 | 95.31 | 3.29 | 4.23 | 367 | 66.21 |

By region (anchor: linear 335,447, d9 1,693, d46 2,996, none 23,940):

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 298,748 | GRCh38-only graph | 0.44 | 1.39 | 93.78 | 96.56 | 28.31 | 0.01 | 99.64 | 0.9943 | - |
| SMaHT easy | 298,748 | HPRC v1.1 d9 | 0.26 | 2.62 | 88.03 | 95.89 | 27.97 | 0.02 | 99.54 | 0.9965 | 83.69 |
| SMaHT easy | 298,748 | HPRC v2.1 d46 | 0.18 | 2.59 | 87.97 | 95.94 | 28.14 | 0.00 | 98.86 | 0.9965 | 83.50 |
| SMaHT difficult | 5,168 | GRCh38-only graph | 2.94 | 3.52 | 75.94 | 82.00 | 38.38 | 7.66 | 87.96 | 0.9934 | - |
| SMaHT difficult | 5,168 | HPRC v1.1 d9 | 2.03 | 6.97 | 62.41 | 80.09 | 37.61 | 7.58 | 87.18 | 0.9961 | 64.72 |
| SMaHT difficult | 5,168 | HPRC v2.1 d46 | 10.57 | 6.79 | 65.51 | 86.59 | 32.39 | 0.48 | 87.91 | 0.9957 | 65.36 |
| SMaHT extreme | 31,583 | GRCh38-only graph | 9.52 | 1.46 | 81.22 | 83.31 | 45.40 | 9.27 | 53.85 | 0.9908 | - |
| SMaHT extreme | 31,583 | HPRC v1.1 d9 | 10.50 | 3.23 | 76.93 | 82.03 | 43.52 | 10.09 | 53.10 | 0.9946 | 66.06 |
| SMaHT extreme | 31,583 | HPRC v2.1 d46 | 13.99 | 3.04 | 86.79 | 90.67 | 38.93 | 1.35 | 49.56 | 0.9920 | 64.47 |
| segdup | 19,187 | GRCh38-only graph | 5.26 | 1.29 | 93.54 | 95.31 | 34.44 | 0.99 | 65.07 | 0.9919 | - |
| segdup | 19,187 | HPRC v1.1 d9 | 4.86 | 2.48 | 88.32 | 94.36 | 33.09 | 1.65 | 62.48 | 0.9953 | 80.93 |
| segdup | 19,187 | HPRC v2.1 d46 | 1.34 | 2.23 | 90.30 | 95.41 | 34.72 | 0.87 | 58.23 | 0.9947 | 78.38 |
| not segdup | 316,312 | GRCh38-only graph | 1.09 | 1.43 | 92.36 | 95.19 | 29.67 | 0.92 | 97.28 | 0.9942 | - |
| not segdup | 316,312 | HPRC v1.1 d9 | 1.03 | 2.76 | 86.60 | 94.47 | 29.23 | 0.95 | 97.31 | 0.9965 | 81.95 |
| not segdup | 316,312 | HPRC v2.1 d46 | 1.66 | 2.73 | 87.39 | 95.38 | 28.74 | 0.08 | 96.86 | 0.9964 | 81.88 |
| in HG008T benchmark BED | 322,762 | GRCh38-only graph | 0.71 | 1.46 | 92.43 | 95.41 | 29.48 | 0.78 | 97.14 | 0.9943 | - |
| in HG008T benchmark BED | 322,762 | HPRC v1.1 d9 | 0.62 | 2.81 | 86.56 | 94.67 | 28.99 | 0.81 | 97.08 | 0.9965 | 82.86 |
| in HG008T benchmark BED | 322,762 | HPRC v2.1 d46 | 1.38 | 2.77 | 87.27 | 95.42 | 28.57 | 0.06 | 97.01 | 0.9965 | 82.91 |
| outside HG008T benchmark BED | 13,009 | GRCh38-only graph | 16.78 | 0.56 | 92.35 | 89.01 | 44.01 | 5.19 | 46.43 | 0.9754 | - |
| outside HG008T benchmark BED | 13,009 | HPRC v1.1 d9 | 16.92 | 0.99 | 90.85 | 88.32 | 43.38 | 6.60 | 44.43 | 0.9824 | 51.98 |
| outside HG008T benchmark BED | 13,009 | HPRC v2.1 d46 | 8.13 | 0.85 | 95.52 | 94.13 | 43.31 | 1.87 | 30.85 | 0.9141 | 47.43 |

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.03 | 0.10 | 1.24 | 0.00 | 99.97 | 11 | 95.48 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.11 | 0.48 | 1.22 | 0.02 | 99.78 | 328 | 95.64 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.17 | 0.46 | 0.05 | 0.08 | 99.77 | 323 | 88.89 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 0.82 | 1.41 | 3.44 | 0.00 | 90.56 | 14 | 92.70 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 0.82 | 5.30 | 3.35 | 0.08 | 94.66 | 80 | 92.49 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.34 | 5.22 | 0.23 | 0.41 | 95.07 | 83 | 91.67 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 4.93 | 6.92 | 1.79 | 0.02 | 72.77 | 925 | 96.82 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 7.09 | 14.12 | 1.70 | 0.12 | 70.38 | 1,285 | 96.27 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7.22 | 11.89 | 0.12 | 0.32 | 72.66 | 1,190 | 71.05 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 3.16 | 6.41 | 1.19 | 0.01 | 84.92 | 342 | 96.07 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 4.96 | 10.88 | 1.07 | 0.13 | 76.98 | 846 | 96.12 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 5.40 | 8.06 | 0.10 | 0.35 | 78.56 | 757 | 85.00 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 0.29 | 0.35 | 1.33 | 0.01 | 98.30 | 608 | 95.51 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 0.39 | 1.03 | 1.31 | 0.02 | 98.69 | 847 | 95.57 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.44 | 1.00 | 0.06 | 0.09 | 98.78 | 839 | 85.94 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 0.23 | 0.40 | 1.36 | 0.00 | 98.59 | 215 | 95.59 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 0.34 | 0.92 | 1.33 | 0.02 | 98.81 | 952 | 95.68 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.40 | 0.81 | 0.06 | 0.10 | 98.88 | 924 | 87.56 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 7.30 | 9.82 | 0.46 | 0.03 | 65.42 | 737 | 91.67 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 10.82 | 22.60 | 0.41 | 0.12 | 54.93 | 750 | 88.68 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 11.01 | 19.93 | 0.08 | 0.23 | 58.28 | 684 | 54.55 |

### W. Window sample (148 GRCh38 windows, flank 25,000 bp; reads placed in a window by at least one arm: 9,000; random 1,616, random_extra 6,753, segdup 631)

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 94.20 | 4.07 | 1.73 |
| GRCh38-only graph | genome | 92.29 | 5.21 | 2.50 |
| HPRC v1.1 d9 | genome | 93.62 | 4.07 | 2.31 |
| HPRC v2.1 d46 | genome | 92.42 | 3.70 | 3.88 |

Unmapped records of the whole file (GAM head / BAM unplaced): GRCh38 BAM 141,830

Reads placed in the window by the arm:

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 8,478 | 1.07 | 94.17 | 96.74 | 36.77 | 0.42 | 2.24 | 0.83 | 1.40 | 95.10 | 0.9940 | 5.64 | 5.17 | 2.40 | - | - | - | 0.00 |
| GRCh38-only graph | 8,306 | 1.06 | 93.82 | 96.16 | 28.91 | 0.53 | 1.85 | 1.10 | 0.88 | 95.64 | 0.9941 | 4.95 | 4.41 | 1.35 | - | - | - | 0.00 |
| HPRC v1.1 d9 | 8,426 | 2.44 | 87.10 | 94.48 | 29.44 | 1.51 | 2.33 | 1.03 | 0.95 | 94.18 | 0.9965 | 4.03 | 4.31 | 1.21 | 84.19 | 8.63 | 96.42 | 0.00 |
| HPRC v2.1 d46 | 8,318 | 2.42 | 88.45 | 95.84 | 28.59 | 0.16 | 2.16 | 1.61 | 1.38 | 94.69 | 0.9964 | 3.97 | 4.38 | 1.18 | 83.53 | 5.66 | 96.78 | 0.00 |

Same reads, both arms place them in the window (B relative to A); last columns: reads A places in the window that B places elsewhere / leaves unmapped:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both  A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 8,156 | 8,156 | 1.37 | 1.58 | 0.29 | 0.42 | 1.08 | 2.46 | 1.39 | 0.00 | 28.37 | 28.18 | 99.61 | 99.91 | 7 | 143 | 7 | 249 | 21 |
| HPRC v2.1 d46 vs GRCh38-only graph | 8,055 | 8,055 | 1.51 | 2.15 | 0.42 | 0.72 | 1.09 | 2.46 | 1.37 | 0.00 | 28.08 | 27.73 | 99.45 | 99.94 | 5 | 188 | 63 | 223 | 40 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 8,096 | 8,096 | 2.00 | 2.00 | 0.63 | 0.75 | 2.53 | 2.47 | 0.06 | 0.12 | 27.95 | 27.80 | 99.49 | 99.97 | 2 | 202 | 128 | 197 | 25 |
| GRCh38-only graph vs GRCh38 BAM | 8,170 | 8,170 | 1.65 | 2.47 | 0.27 | 1.10 | 1.09 | 1.08 | 0.00 | 0.01 | 34.93 | 28.37 | 99.62 | 100.00 | 0 | 94 | 214 | 94 | 42 |
| HPRC v1.1 d9 vs GRCh38 BAM | 8,137 | 8,137 | 1.46 | 2.59 | 0.22 | 1.19 | 1.11 | 2.52 | 1.41 | 0.00 | 34.96 | 28.18 | 99.52 | 99.91 | 7 | 142 | 199 | 172 | 117 |
| HPRC v2.1 d46 vs GRCh38 BAM | 8,080 | 8,080 | 1.37 | 2.33 | 0.16 | 1.10 | 1.13 | 2.46 | 1.37 | 0.04 | 34.93 | 27.98 | 99.37 | 99.94 | 5 | 206 | 192 | 237 | 1 |

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 113 | 95.58 | 2.0 | 65.49 | 34.51 | 0.88 | 0.00 | 2.0 | 100.00 | 1.77 | 1.77 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 110 | 94.55 | 2.0 | 66.36 | 33.64 | 0.91 | 0.00 | 1.5 | 100.00 | 0.00 | 2.73 | 0 | n/a |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 5 | 80.00 | 2.0 | 0.00 | 100.00 | 0.00 | 0.00 | 0.0 | 100.00 | 20.00 | 20.00 | 10 | 60.00 |
| GRCh38-only graph vs GRCh38 BAM | 0 | | | | | | | | | | | | |
| HPRC v1.1 d9 vs GRCh38 BAM | 115 | 95.65 | 2.0 | 62.61 | 34.78 | 5.22 | 0.00 | 1.0 | 98.26 | 4.35 | 2.61 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38 BAM | 111 | 94.59 | 2.0 | 63.06 | 34.23 | 5.41 | 0.00 | 1.0 | 100.00 | 2.70 | 2.70 | 3 | 0.00 |

By window kind and region:

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 1,616 | GRCh38 BAM | 92.51 | 3.28 | 4.21 | 1.27 | 94.18 | 96.92 | 37.46 | 0.40 | 95.79 | 0.9939 |
| windows: random | 1,616 | GRCh38-only graph | 90.84 | 6.87 | 2.29 | 1.29 | 93.53 | 96.39 | 29.90 | 0.61 | 95.64 | 0.9940 |
| windows: random | 1,616 | HPRC v1.1 d9 | 95.30 | 2.66 | 2.04 | 2.34 | 84.94 | 92.34 | 31.88 | 3.77 | 92.01 | 0.9966 |
| windows: random | 1,616 | HPRC v2.1 d46 | 89.98 | 4.02 | 6.00 | 2.41 | 88.31 | 96.35 | 28.95 | 0.00 | 96.22 | 0.9964 |
| windows: random_extra | 6,753 | GRCh38 BAM | 95.71 | 3.07 | 1.23 | 1.05 | 94.09 | 96.61 | 35.94 | 0.25 | 96.33 | 0.9942 |
| windows: random_extra | 6,753 | GRCh38-only graph | 93.88 | 3.54 | 2.58 | 1.03 | 93.79 | 96.01 | 28.09 | 0.49 | 96.99 | 0.9943 |
| windows: random_extra | 6,753 | HPRC v1.1 d9 | 94.25 | 3.33 | 2.41 | 2.50 | 87.16 | 94.83 | 28.39 | 1.05 | 96.51 | 0.9965 |
| windows: random_extra | 6,753 | HPRC v2.1 d46 | 94.08 | 2.47 | 3.45 | 2.46 | 88.05 | 95.66 | 28.02 | 0.16 | 96.32 | 0.9965 |
| windows: segdup | 631 | GRCh38 BAM | 82.41 | 16.80 | 0.79 | 0.77 | 95.19 | 97.88 | 45.00 | 2.69 | 77.88 | 0.9918 |
| windows: segdup | 631 | GRCh38-only graph | 78.92 | 18.86 | 2.22 | 0.80 | 95.18 | 97.39 | 36.35 | 0.80 | 78.51 | 0.9922 |
| windows: segdup | 631 | HPRC v1.1 d9 | 82.57 | 15.53 | 1.90 | 2.11 | 92.71 | 96.55 | 35.12 | 0.38 | 72.17 | 0.9960 |
| windows: segdup | 631 | HPRC v2.1 d46 | 80.98 | 16.01 | 3.01 | 1.96 | 93.74 | 96.67 | 34.64 | 0.59 | 70.06 | 0.9957 |
| SMaHT easy | 7,529 | GRCh38 BAM | 98.63 | 0.24 | 1.13 | 0.97 | 94.60 | 97.06 | 35.87 | 0.00 | 99.19 | 0.9942 |
| SMaHT easy | 7,529 | GRCh38-only graph | 96.60 | 1.38 | 2.02 | 0.99 | 94.75 | 96.98 | 28.09 | 0.01 | 99.55 | 0.9943 |
| SMaHT easy | 7,529 | HPRC v1.1 d9 | 98.01 | 0.20 | 1.79 | 2.20 | 87.87 | 95.20 | 29.00 | 1.12 | 98.27 | 0.9966 |
| SMaHT easy | 7,529 | HPRC v2.1 d46 | 96.84 | 0.35 | 2.82 | 2.18 | 88.88 | 96.38 | 28.12 | 0.00 | 99.16 | 0.9965 |
| SMaHT difficult | 160 | GRCh38 BAM | 81.25 | 10.62 | 8.12 | 3.85 | 82.31 | 93.85 | 34.62 | 0.00 | 93.85 | 0.9905 |
| SMaHT difficult | 160 | GRCh38-only graph | 89.38 | 10.00 | 0.62 | 3.50 | 79.02 | 88.11 | 35.66 | 6.99 | 90.21 | 0.9904 |
| SMaHT difficult | 160 | HPRC v1.1 d9 | 85.62 | 13.75 | 0.62 | 9.49 | 64.23 | 84.67 | 33.58 | 5.84 | 92.70 | 0.9951 |
| SMaHT difficult | 160 | HPRC v2.1 d46 | 85.62 | 2.50 | 11.88 | 9.49 | 64.96 | 87.59 | 35.04 | 0.00 | 90.51 | 0.9952 |
| SMaHT extreme | 1,311 | GRCh38 BAM | 70.33 | 25.25 | 4.42 | 1.52 | 92.41 | 94.58 | 44.25 | 3.90 | 62.36 | 0.9913 |
| SMaHT extreme | 1,311 | GRCh38-only graph | 67.89 | 26.62 | 5.49 | 1.24 | 88.65 | 90.79 | 34.49 | 3.71 | 64.61 | 0.9926 |
| SMaHT extreme | 1,311 | HPRC v1.1 d9 | 69.41 | 25.10 | 5.49 | 3.41 | 84.29 | 90.11 | 32.42 | 3.96 | 61.32 | 0.9962 |
| SMaHT extreme | 1,311 | HPRC v2.1 d46 | 67.89 | 23.11 | 9.00 | 3.26 | 88.54 | 92.70 | 31.46 | 1.46 | 58.65 | 0.9954 |
| segdup | 1,038 | GRCh38 BAM | 78.71 | 20.62 | 0.67 | 0.98 | 93.51 | 96.94 | 40.27 | 3.67 | 72.34 | 0.9925 |
| segdup | 1,038 | GRCh38-only graph | 75.63 | 22.06 | 2.31 | 1.02 | 92.87 | 96.43 | 32.48 | 0.89 | 69.68 | 0.9928 |
| segdup | 1,038 | HPRC v1.1 d9 | 78.61 | 19.27 | 2.12 | 2.21 | 89.09 | 96.08 | 29.78 | 1.10 | 64.95 | 0.9961 |
| segdup | 1,038 | HPRC v2.1 d46 | 77.94 | 19.08 | 2.99 | 2.22 | 90.48 | 96.17 | 31.40 | 0.87 | 64.03 | 0.9957 |
| not segdup | 7,962 | GRCh38 BAM | 96.22 | 1.91 | 1.87 | 1.08 | 94.24 | 96.72 | 36.39 | 0.08 | 97.53 | 0.9941 |
| not segdup | 7,962 | GRCh38-only graph | 94.46 | 3.01 | 2.52 | 1.06 | 93.92 | 96.13 | 28.53 | 0.49 | 98.35 | 0.9942 |
| not segdup | 7,962 | HPRC v1.1 d9 | 95.58 | 2.08 | 2.34 | 2.47 | 86.89 | 94.31 | 29.41 | 1.55 | 97.32 | 0.9966 |
| not segdup | 7,962 | HPRC v2.1 d46 | 94.31 | 1.70 | 3.99 | 2.44 | 88.23 | 95.81 | 28.29 | 0.08 | 97.99 | 0.9965 |
| in HG008T benchmark BED | 8,491 | GRCh38 BAM | 96.69 | 1.55 | 1.75 | 1.07 | 94.18 | 96.85 | 36.24 | 0.16 | 97.11 | 0.9941 |
| in HG008T benchmark BED | 8,491 | GRCh38-only graph | 94.94 | 2.81 | 2.25 | 1.09 | 93.82 | 96.34 | 28.74 | 0.46 | 97.22 | 0.9942 |
| in HG008T benchmark BED | 8,491 | HPRC v1.1 d9 | 96.13 | 1.81 | 2.06 | 2.47 | 86.95 | 94.60 | 29.20 | 1.43 | 96.07 | 0.9966 |
| in HG008T benchmark BED | 8,491 | HPRC v2.1 d46 | 94.52 | 1.81 | 3.66 | 2.48 | 88.15 | 95.96 | 28.23 | 0.02 | 97.13 | 0.9965 |
| outside HG008T benchmark BED | 509 | GRCh38 BAM | 52.65 | 45.97 | 1.38 | 1.12 | 94.03 | 93.66 | 52.99 | 8.58 | 33.58 | 0.9868 |
| outside HG008T benchmark BED | 509 | GRCh38-only graph | 48.13 | 45.19 | 6.68 | 0.00 | 93.88 | 90.20 | 34.29 | 2.86 | 43.67 | 0.9868 |
| outside HG008T benchmark BED | 509 | HPRC v1.1 d9 | 51.87 | 41.65 | 6.48 | 1.52 | 91.67 | 90.91 | 37.12 | 3.79 | 35.98 | 0.9926 |
| outside HG008T benchmark BED | 509 | HPRC v2.1 d46 | 57.37 | 35.17 | 7.47 | 0.68 | 96.58 | 92.47 | 38.36 | 3.77 | 27.40 | 0.9875 |

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | HPRC v1.1 d9 vs GRCh38-only graph | 1,445 | 0.76 | 0.28 | 1.18 | 0.00 | 99.93 | 1.57 | 0.00 |
| windows: random | HPRC v2.1 d46 vs GRCh38-only graph | 1,427 | 0.28 | 0.49 | 1.12 | 0.00 | 99.93 | 2.18 | 0.61 |
| windows: random | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,433 | 0.28 | 0.84 | 0.07 | 0.14 | 100.00 | 3.05 | 3.90 |
| windows: random | GRCh38-only graph vs GRCh38 BAM | 1,450 | 0.07 | 0.55 | 0.00 | 0.00 | 99.72 | 0.67 | 2.34 |
| windows: random | HPRC v1.1 d9 vs GRCh38 BAM | 1,441 | 0.21 | 0.42 | 1.18 | 0.00 | 99.72 | 1.40 | 2.21 |
| windows: random | HPRC v2.1 d46 vs GRCh38 BAM | 1,428 | 0.07 | 0.77 | 1.12 | 0.00 | 99.72 | 2.47 | 2.01 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38-only graph | 6,258 | 0.18 | 0.27 | 1.44 | 0.00 | 99.52 | 1.20 | 0.09 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38-only graph | 6,196 | 0.34 | 0.56 | 1.45 | 0.00 | 99.44 | 1.55 | 0.73 |
| windows: random_extra | HPRC v2.1 d46 vs HPRC v1.1 d9 | 6,208 | 0.48 | 0.58 | 0.05 | 0.08 | 99.48 | 1.43 | 1.04 |
| windows: random_extra | GRCh38-only graph vs GRCh38 BAM | 6,243 | 0.27 | 0.69 | 0.00 | 0.02 | 99.57 | 0.82 | 2.58 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38 BAM | 6,243 | 0.19 | 0.70 | 1.46 | 0.00 | 99.46 | 1.01 | 2.40 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38 BAM | 6,217 | 0.14 | 0.66 | 1.46 | 0.05 | 99.39 | 1.50 | 2.31 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 453 | 0.44 | 2.87 | 1.32 | 0.00 | 99.78 | 8.84 | 0.20 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 432 | 2.08 | 3.70 | 0.93 | 0.00 | 98.15 | 11.65 | 1.61 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 455 | 3.74 | 2.86 | 0.22 | 0.66 | 98.02 | 12.28 | 0.38 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 477 | 0.84 | 8.18 | 0.00 | 0.00 | 100.00 | 5.96 | 2.31 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 453 | 0.66 | 10.38 | 1.55 | 0.00 | 99.78 | 10.77 | 2.12 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 435 | 0.69 | 8.51 | 0.92 | 0.00 | 97.93 | 13.85 | 2.50 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 7,272 | 0.07 | 0.07 | 1.24 | 0.00 | 100.00 | 0.00 | 0.01 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 7,248 | 0.18 | 0.28 | 1.20 | 0.00 | 100.00 | 0.19 | 0.15 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7,263 | 0.21 | 0.30 | 0.03 | 0.07 | 100.00 | 0.33 | 1.25 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 7,273 | 0.00 | 0.22 | 0.00 | 0.00 | 100.00 | 0.05 | 2.01 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 7,285 | 0.01 | 0.23 | 1.24 | 0.00 | 100.00 | 0.08 | 1.82 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 7,278 | 0.00 | 0.30 | 1.20 | 0.00 | 100.00 | 0.22 | 1.78 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 132 | 0.00 | 0.76 | 6.06 | 0.00 | 100.00 | 7.69 | 0.00 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 127 | 0.00 | 0.00 | 6.30 | 0.00 | 100.00 | 2.80 | 8.39 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 127 | 0.79 | 0.00 | 0.00 | 0.00 | 100.00 | 2.19 | 5.11 |
| SMaHT difficult | GRCh38-only graph vs GRCh38 BAM | 128 | 0.78 | 0.00 | 0.00 | 0.00 | 100.00 | 0.77 | 0.77 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38 BAM | 128 | 0.78 | 0.78 | 6.25 | 0.00 | 100.00 | 0.77 | 0.77 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38 BAM | 126 | 0.79 | 0.00 | 6.35 | 0.00 | 100.00 | 1.54 | 1.54 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 752 | 2.53 | 3.72 | 1.99 | 0.00 | 95.74 | 14.83 | 0.67 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 680 | 3.09 | 5.59 | 2.21 | 0.00 | 93.53 | 19.10 | 4.49 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 706 | 4.96 | 5.52 | 0.42 | 0.71 | 94.19 | 19.23 | 3.19 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 769 | 2.73 | 9.62 | 0.00 | 0.13 | 95.97 | 9.65 | 6.94 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 724 | 2.21 | 10.91 | 2.35 | 0.00 | 94.61 | 14.64 | 6.83 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 676 | 1.78 | 9.91 | 2.37 | 0.44 | 92.46 | 20.39 | 6.29 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 711 | 2.11 | 3.23 | 0.98 | 0.00 | 97.47 | 9.30 | 0.13 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 652 | 2.61 | 4.60 | 1.07 | 0.00 | 95.55 | 15.54 | 1.40 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 679 | 4.42 | 4.86 | 0.44 | 0.44 | 95.73 | 16.18 | 0.61 |
| segdup | GRCh38-only graph vs GRCh38 BAM | 748 | 1.34 | 9.63 | 0.00 | 0.00 | 98.26 | 5.75 | 2.69 |
| segdup | HPRC v1.1 d9 vs GRCh38 BAM | 705 | 1.13 | 10.92 | 1.28 | 0.00 | 96.60 | 11.14 | 2.57 |
| segdup | HPRC v2.1 d46 vs GRCh38 BAM | 655 | 0.92 | 10.08 | 1.22 | 0.00 | 94.96 | 17.38 | 2.45 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 7,445 | 0.12 | 0.15 | 1.42 | 0.00 | 99.81 | 0.93 | 0.08 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 7,403 | 0.23 | 0.38 | 1.39 | 0.00 | 99.80 | 0.88 | 0.69 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7,417 | 0.28 | 0.38 | 0.03 | 0.09 | 99.84 | 0.92 | 1.62 |
| not segdup | GRCh38-only graph vs GRCh38 BAM | 7,422 | 0.16 | 0.24 | 0.00 | 0.01 | 99.76 | 0.61 | 2.51 |
| not segdup | HPRC v1.1 d9 vs GRCh38 BAM | 7,432 | 0.13 | 0.27 | 1.43 | 0.00 | 99.80 | 0.67 | 2.32 |
| not segdup | HPRC v2.1 d46 vs GRCh38 BAM | 7,425 | 0.09 | 0.31 | 1.39 | 0.04 | 99.76 | 0.84 | 2.25 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 7,986 | 0.19 | 0.24 | 1.40 | 0.00 | 99.86 | 0.88 | 0.05 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 7,909 | 0.33 | 0.53 | 1.39 | 0.00 | 99.81 | 1.23 | 0.66 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7,929 | 0.42 | 0.54 | 0.05 | 0.08 | 99.81 | 1.36 | 1.49 |
| in HG008T benchmark BED | GRCh38-only graph vs GRCh38 BAM | 7,997 | 0.11 | 0.74 | 0.00 | 0.00 | 99.86 | 0.37 | 2.23 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38 BAM | 7,989 | 0.08 | 0.76 | 1.41 | 0.00 | 99.82 | 0.65 | 2.05 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38 BAM | 7,940 | 0.08 | 0.76 | 1.40 | 0.00 | 99.75 | 1.28 | 2.01 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 170 | 5.29 | 8.82 | 0.59 | 0.00 | 87.65 | 29.39 | 1.22 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 146 | 5.48 | 10.96 | 0.00 | 0.00 | 80.14 | 36.33 | 4.08 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 167 | 10.78 | 10.78 | 0.60 | 2.40 | 84.43 | 34.47 | 2.27 |
| outside HG008T benchmark BED | GRCh38-only graph vs GRCh38 BAM | 173 | 7.51 | 17.92 | 0.00 | 0.58 | 88.44 | 23.88 | 11.57 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38 BAM | 148 | 8.11 | 24.32 | 1.35 | 0.00 | 83.11 | 33.21 | 11.57 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38 BAM | 140 | 5.00 | 20.71 | 0.00 | 2.14 | 77.86 | 37.69 | 10.07 |

## HG008T_Illumina

### W. Window sample (48 GRCh38 windows, flank 1,000 bp; reads placed in a window by at least one arm: 1,465,754; random 1,129,415, segdup 336,339)

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 92.17 | 7.79 | 0.04 |
| GRCh38-only graph | window | 91.35 | 7.87 | 0.78 |
| HPRC v1.1 d9 | window | 91.39 | 7.86 | 0.74 |
| HPRC v2.1 d46 | window | 91.36 | 7.89 | 0.75 |

Unmapped records of the whole file (GAM head / BAM unplaced): GRCh38 BAM 2,165,824; GRCh38-only graph 44,910,000 (1.51 % of all reads); HPRC v1.1 d9 42,368,000 (1.42 % of all reads); HPRC v2.1 d46 42,009,000 (1.41 % of all reads)

Reads placed in the window by the arm:

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 1,351,053 | 68.50 | 29.40 | 2.94 | 2.84 | 5.06 | 0.40 | 1.97 | 4.38 | 88.20 | 1.0000 | 4.73 | 0.21 | 0.83 | - | - | - | 0.00 |
| GRCh38-only graph | 1,338,964 | 68.85 | 29.40 | 2.96 | 2.19 | 3.43 | 3.31 | 2.77 | 3.61 | 86.89 | 1.0000 | 4.63 | 0.21 | 0.70 | - | - | - | 0.00 |
| HPRC v1.1 d9 | 1,339,624 | 79.15 | 20.09 | 0.62 | 1.94 | 4.29 | 3.44 | 2.67 | 3.36 | 86.24 | 1.0000 | 3.69 | 0.05 | 0.64 | 14.42 | 5.05 | 98.63 | 0.00 |
| HPRC v2.1 d46 | 1,339,133 | 79.24 | 20.01 | 0.63 | 1.92 | 4.24 | 3.27 | 2.66 | 3.54 | 86.29 | 1.0000 | 3.67 | 0.05 | 0.63 | 14.21 | 3.46 | 98.71 | 0.00 |

Same reads, both arms place them in the window (B relative to A); last columns: reads A places in the window that B places elsewhere / leaves unmapped:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both  A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,277,850 | 1,277,850 | 1.59 | 2.71 | 0.43 | 0.92 | 69.12 | 79.36 | 10.24 | 0.00 | 2.05 | 1.83 | 99.71 | 100.00 | 19 | 60,526 | 588 | 60,967 | 807 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,277,730 | 1,277,730 | 1.54 | 3.26 | 0.41 | 0.92 | 69.07 | 79.41 | 10.34 | 0.00 | 2.04 | 1.84 | 99.64 | 100.00 | 19 | 60,650 | 584 | 60,688 | 715 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,293,604 | 1,293,604 | 0.96 | 1.38 | 0.40 | 0.38 | 79.30 | 79.41 | 0.29 | 0.17 | 1.86 | 1.87 | 99.75 | 100.00 | 0 | 45,757 | 263 | 45,315 | 214 |
| GRCh38-only graph vs GRCh38 BAM | 1,296,578 | 1,296,578 | 3.75 | 7.37 | 0.11 | 2.28 | 69.18 | 69.18 | 0.00 | 0.00 | 2.30 | 1.97 | 99.58 | 100.00 | 25 | 43,333 | 11,142 | 41,880 | 506 |
| HPRC v1.1 d9 vs GRCh38 BAM | 1,270,875 | 1,270,875 | 3.39 | 6.85 | 0.26 | 2.89 | 69.26 | 79.61 | 10.34 | 0.00 | 2.27 | 1.72 | 99.72 | 100.00 | 42 | 69,543 | 10,635 | 68,255 | 494 |
| HPRC v2.1 d46 vs GRCh38 BAM | 1,275,980 | 1,275,980 | 3.39 | 7.07 | 0.24 | 2.89 | 69.24 | 79.65 | 10.41 | 0.00 | 2.26 | 1.72 | 99.60 | 100.00 | 43 | 64,381 | 10,692 | 62,664 | 489 |

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 130,826 | 94.15 | 1.0 | 86.02 | 13.37 | 0.69 | 0.00 | 1.0 | 99.74 | 5.31 | 1.79 | 4 | 25.00 |
| HPRC v2.1 d46 vs GRCh38-only graph | 132,115 | 94.12 | 1.0 | 86.09 | 13.29 | 0.72 | 0.00 | 1.0 | 99.72 | 5.05 | 1.74 | 5 | 20.00 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 3,687 | 93.25 | 1.0 | 86.06 | 12.99 | 1.06 | 0.00 | 1.0 | 99.40 | 8.06 | 1.71 | 2,241 | 17.89 |
| GRCh38-only graph vs GRCh38 BAM | 22 | 0.00 | 0.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 4.55 | 40.91 | 9.09 | 40 | 0.00 |
| HPRC v1.1 d9 vs GRCh38 BAM | 131,452 | 94.16 | 1.0 | 84.74 | 13.93 | 1.43 | 0.00 | 1.0 | 99.72 | 6.35 | 3.89 | 13 | 23.08 |
| HPRC v2.1 d46 vs GRCh38 BAM | 132,810 | 94.14 | 1.0 | 84.84 | 13.83 | 1.44 | 0.00 | 1.0 | 99.72 | 6.17 | 3.82 | 13 | 30.77 |

By window kind and region:

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 1,129,415 | GRCh38 BAM | 96.71 | 3.26 | 0.04 | 69.72 | 28.18 | 2.83 | 2.68 | 2.28 | 95.33 | 1.0000 |
| windows: random | 1,129,415 | GRCh38-only graph | 96.02 | 3.16 | 0.82 | 70.17 | 28.07 | 2.81 | 2.05 | 1.52 | 93.72 | 1.0000 |
| windows: random | 1,129,415 | HPRC v1.1 d9 | 95.97 | 3.25 | 0.78 | 80.25 | 19.03 | 0.45 | 1.79 | 1.76 | 93.78 | 1.0000 |
| windows: random | 1,129,415 | HPRC v2.1 d46 | 95.55 | 3.67 | 0.78 | 80.54 | 18.78 | 0.41 | 1.78 | 1.66 | 93.88 | 1.0000 |
| windows: segdup | 336,339 | GRCh38 BAM | 76.96 | 23.00 | 0.04 | 63.37 | 34.56 | 3.39 | 3.52 | 16.77 | 58.09 | 1.0000 |
| windows: segdup | 336,339 | GRCh38-only graph | 75.66 | 23.69 | 0.65 | 63.20 | 35.03 | 3.60 | 2.82 | 11.56 | 57.75 | 1.0000 |
| windows: segdup | 336,339 | HPRC v1.1 d9 | 76.03 | 23.35 | 0.62 | 74.49 | 24.57 | 1.34 | 2.55 | 15.04 | 54.28 | 1.0000 |
| windows: segdup | 336,339 | HPRC v2.1 d46 | 77.30 | 22.07 | 0.63 | 73.84 | 25.13 | 1.57 | 2.52 | 14.91 | 54.80 | 1.0000 |
| SMaHT easy | 918,984 | GRCh38 BAM | 99.89 | 0.08 | 0.03 | 71.55 | 26.87 | 2.15 | 1.64 | 0.01 | 99.81 | 1.0000 |
| SMaHT easy | 918,984 | GRCh38-only graph | 99.28 | 0.08 | 0.65 | 71.94 | 26.70 | 2.13 | 1.21 | 0.12 | 97.77 | 1.0000 |
| SMaHT easy | 918,984 | HPRC v1.1 d9 | 99.26 | 0.11 | 0.63 | 81.89 | 17.68 | 0.30 | 1.06 | 0.16 | 97.75 | 1.0000 |
| SMaHT easy | 918,984 | HPRC v2.1 d46 | 99.28 | 0.09 | 0.63 | 82.04 | 17.53 | 0.28 | 1.07 | 0.18 | 97.45 | 1.0000 |
| SMaHT difficult | 154,161 | GRCh38 BAM | 99.21 | 0.72 | 0.07 | 61.35 | 34.91 | 5.92 | 5.91 | 0.36 | 97.78 | 1.0000 |
| SMaHT difficult | 154,161 | GRCh38-only graph | 98.46 | 0.60 | 0.95 | 61.77 | 34.97 | 5.94 | 5.02 | 0.52 | 94.55 | 1.0000 |
| SMaHT difficult | 154,161 | HPRC v1.1 d9 | 98.21 | 0.99 | 0.80 | 73.31 | 25.15 | 1.07 | 4.36 | 0.40 | 94.88 | 1.0000 |
| SMaHT difficult | 154,161 | HPRC v2.1 d46 | 98.28 | 0.91 | 0.81 | 73.56 | 24.92 | 1.01 | 4.37 | 0.57 | 94.55 | 1.0000 |
| SMaHT extreme | 392,609 | GRCh38 BAM | 71.36 | 28.59 | 0.05 | 62.44 | 34.69 | 3.90 | 5.11 | 24.17 | 44.92 | 1.0000 |
| SMaHT extreme | 392,609 | GRCh38-only graph | 70.01 | 28.96 | 1.03 | 62.50 | 35.26 | 4.07 | 3.92 | 16.03 | 46.53 | 1.0000 |
| SMaHT extreme | 392,609 | HPRC v1.1 d9 | 70.31 | 28.70 | 0.99 | 73.30 | 25.30 | 1.44 | 3.49 | 20.09 | 43.45 | 1.0000 |
| SMaHT extreme | 392,609 | HPRC v2.1 d46 | 70.11 | 28.90 | 0.99 | 73.07 | 25.52 | 1.59 | 3.41 | 19.70 | 44.76 | 1.0000 |
| segdup | 385,653 | GRCh38 BAM | 72.04 | 27.93 | 0.04 | 64.85 | 33.07 | 3.21 | 3.61 | 23.16 | 50.48 | 1.0000 |
| segdup | 385,653 | GRCh38-only graph | 70.94 | 28.39 | 0.67 | 64.82 | 33.45 | 3.40 | 2.83 | 14.76 | 50.23 | 1.0000 |
| segdup | 385,653 | HPRC v1.1 d9 | 71.42 | 27.94 | 0.64 | 74.77 | 24.27 | 1.29 | 2.59 | 19.04 | 46.78 | 1.0000 |
| segdup | 385,653 | HPRC v2.1 d46 | 71.32 | 28.04 | 0.64 | 74.38 | 24.61 | 1.49 | 2.53 | 18.67 | 48.26 | 1.0000 |
| not segdup | 1,080,101 | GRCh38 BAM | 99.36 | 0.60 | 0.04 | 69.45 | 28.45 | 2.87 | 2.64 | 0.37 | 97.96 | 1.0000 |
| not segdup | 1,080,101 | GRCh38-only graph | 98.64 | 0.54 | 0.82 | 69.88 | 28.35 | 2.84 | 2.03 | 0.52 | 96.30 | 1.0000 |
| not segdup | 1,080,101 | HPRC v1.1 d9 | 98.53 | 0.69 | 0.78 | 80.29 | 19.01 | 0.45 | 1.77 | 0.48 | 96.45 | 1.0000 |
| not segdup | 1,080,101 | HPRC v2.1 d46 | 98.52 | 0.70 | 0.78 | 80.50 | 18.82 | 0.41 | 1.77 | 0.51 | 96.12 | 1.0000 |
| in HG008T benchmark BED | 1,362,171 | GRCh38 BAM | 95.29 | 4.68 | 0.04 | 69.01 | 28.95 | 2.81 | 2.70 | 3.22 | 91.15 | 1.0000 |
| in HG008T benchmark BED | 1,362,171 | GRCh38-only graph | 94.30 | 4.90 | 0.79 | 69.46 | 28.83 | 2.80 | 2.07 | 2.18 | 89.79 | 1.0000 |
| in HG008T benchmark BED | 1,362,171 | HPRC v1.1 d9 | 93.82 | 5.42 | 0.76 | 79.96 | 19.34 | 0.47 | 1.81 | 2.49 | 89.67 | 1.0000 |
| in HG008T benchmark BED | 1,362,171 | HPRC v2.1 d46 | 93.89 | 5.35 | 0.76 | 80.03 | 19.28 | 0.48 | 1.81 | 2.50 | 89.63 | 1.0000 |
| outside HG008T benchmark BED | 103,583 | GRCh38 BAM | 51.26 | 48.71 | 0.03 | 56.06 | 40.49 | 6.05 | 6.34 | 50.13 | 16.05 | 1.0000 |
| outside HG008T benchmark BED | 103,583 | GRCh38-only graph | 52.52 | 46.86 | 0.62 | 54.45 | 42.67 | 6.63 | 5.16 | 32.93 | 18.33 | 1.0000 |
| outside HG008T benchmark BED | 103,583 | HPRC v1.1 d9 | 59.47 | 39.96 | 0.57 | 62.33 | 35.73 | 3.72 | 4.55 | 41.74 | 15.07 | 1.0000 |
| outside HG008T benchmark BED | 103,583 | HPRC v2.1 d46 | 58.16 | 41.25 | 0.59 | 62.55 | 35.56 | 3.81 | 4.37 | 41.10 | 15.47 | 1.0000 |

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | HPRC v1.1 d9 vs GRCh38-only graph | 1,066,836 | 0.21 | 0.26 | 10.10 | 0.00 | 99.89 | 1.59 | 0.04 |
| windows: random | HPRC v2.1 d46 vs GRCh38-only graph | 1,061,925 | 0.18 | 0.41 | 10.36 | 0.00 | 99.91 | 2.04 | 0.04 |
| windows: random | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,063,750 | 0.11 | 0.27 | 0.33 | 0.08 | 99.91 | 1.84 | 0.02 |
| windows: random | GRCh38-only graph vs GRCh38 BAM | 1,068,255 | 0.06 | 1.38 | 0.00 | 0.00 | 99.91 | 1.36 | 0.83 |
| windows: random | HPRC v1.1 d9 vs GRCh38 BAM | 1,063,674 | 0.12 | 1.48 | 10.13 | 0.00 | 99.90 | 1.82 | 0.79 |
| windows: random | HPRC v2.1 d46 vs GRCh38 BAM | 1,062,167 | 0.10 | 1.62 | 10.36 | 0.00 | 99.91 | 1.96 | 0.79 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 211,014 | 1.54 | 4.25 | 10.92 | 0.00 | 98.80 | 17.03 | 0.05 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 215,805 | 1.53 | 3.44 | 10.23 | 0.00 | 98.30 | 15.15 | 0.05 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 229,854 | 1.73 | 0.89 | 0.09 | 0.61 | 99.02 | 10.08 | 0.03 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 228,323 | 0.35 | 6.51 | 0.01 | 0.01 | 98.03 | 10.98 | 0.81 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 207,201 | 1.00 | 10.12 | 11.43 | 0.00 | 98.81 | 19.18 | 0.77 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 213,813 | 0.95 | 9.24 | 10.64 | 0.00 | 98.10 | 16.61 | 0.78 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 911,336 | 0.06 | 0.15 | 9.97 | 0.00 | 100.00 | 0.08 | 0.03 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 911,692 | 0.06 | 0.24 | 10.12 | 0.00 | 100.00 | 0.04 | 0.03 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 911,830 | 0.10 | 0.20 | 0.25 | 0.10 | 100.00 | 0.03 | 0.01 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 911,778 | 0.01 | 0.97 | 0.00 | 0.00 | 100.00 | 0.03 | 0.64 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 911,289 | 0.01 | 1.07 | 9.98 | 0.00 | 100.00 | 0.10 | 0.63 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 911,629 | 0.01 | 1.17 | 10.13 | 0.00 | 100.00 | 0.06 | 0.63 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 151,057 | 0.20 | 0.30 | 11.41 | 0.00 | 99.98 | 0.43 | 0.05 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 151,144 | 0.20 | 0.43 | 11.67 | 0.00 | 99.99 | 0.37 | 0.05 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 151,263 | 0.17 | 0.30 | 0.46 | 0.20 | 99.99 | 0.06 | 0.03 |
| SMaHT difficult | GRCh38-only graph vs GRCh38 BAM | 150,834 | 0.10 | 1.78 | 0.00 | 0.00 | 99.97 | 0.46 | 0.92 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38 BAM | 150,601 | 0.11 | 1.89 | 11.46 | 0.00 | 99.98 | 0.74 | 0.79 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38 BAM | 150,707 | 0.10 | 2.04 | 11.71 | 0.00 | 99.99 | 0.67 | 0.80 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 215,457 | 2.11 | 4.61 | 10.56 | 0.00 | 98.31 | 21.52 | 0.09 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 214,894 | 2.05 | 4.14 | 10.36 | 0.00 | 97.86 | 21.73 | 0.09 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 230,511 | 1.74 | 1.15 | 0.30 | 0.44 | 98.64 | 16.46 | 0.03 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 233,966 | 0.52 | 7.74 | 0.01 | 0.01 | 97.71 | 15.11 | 1.38 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 208,985 | 1.49 | 11.54 | 11.11 | 0.00 | 98.32 | 24.09 | 1.32 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 213,644 | 1.33 | 10.86 | 10.67 | 0.00 | 97.64 | 22.42 | 1.32 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 216,721 | 1.75 | 4.39 | 9.71 | 0.00 | 98.34 | 20.74 | 0.05 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 216,419 | 1.67 | 3.67 | 9.32 | 0.00 | 97.90 | 20.85 | 0.05 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 231,067 | 1.85 | 1.03 | 0.22 | 0.52 | 98.65 | 16.08 | 0.03 |
| segdup | GRCh38-only graph vs GRCh38 BAM | 235,178 | 0.36 | 6.77 | 0.01 | 0.01 | 97.80 | 14.46 | 0.89 |
| segdup | HPRC v1.1 d9 vs GRCh38 BAM | 210,765 | 1.22 | 10.51 | 10.26 | 0.00 | 98.37 | 23.29 | 0.85 |
| segdup | HPRC v2.1 d46 vs GRCh38 BAM | 215,675 | 1.07 | 9.61 | 9.66 | 0.00 | 97.70 | 21.51 | 0.86 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 1,061,129 | 0.15 | 0.21 | 10.35 | 0.00 | 99.99 | 0.36 | 0.04 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1,061,311 | 0.15 | 0.36 | 10.55 | 0.00 | 99.99 | 0.34 | 0.04 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,062,537 | 0.09 | 0.24 | 0.30 | 0.10 | 99.99 | 0.14 | 0.02 |
| not segdup | GRCh38-only graph vs GRCh38 BAM | 1,061,400 | 0.06 | 1.29 | 0.00 | 0.00 | 99.98 | 0.29 | 0.81 |
| not segdup | HPRC v1.1 d9 vs GRCh38 BAM | 1,060,110 | 0.07 | 1.37 | 10.36 | 0.00 | 99.99 | 0.45 | 0.77 |
| not segdup | HPRC v2.1 d46 vs GRCh38 BAM | 1,060,305 | 0.07 | 1.53 | 10.56 | 0.00 | 99.99 | 0.43 | 0.77 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 1,244,554 | 0.40 | 0.84 | 10.38 | 0.00 | 99.76 | 3.07 | 0.04 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 1,245,168 | 0.38 | 0.85 | 10.48 | 0.00 | 99.66 | 3.03 | 0.04 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,246,924 | 0.40 | 0.38 | 0.29 | 0.17 | 99.80 | 2.41 | 0.02 |
| in HG008T benchmark BED | GRCh38-only graph vs GRCh38 BAM | 1,260,986 | 0.09 | 2.10 | 0.00 | 0.00 | 99.60 | 2.03 | 0.81 |
| in HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38 BAM | 1,241,608 | 0.22 | 2.66 | 10.46 | 0.00 | 99.75 | 3.56 | 0.78 |
| in HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38 BAM | 1,246,508 | 0.20 | 2.67 | 10.53 | 0.00 | 99.63 | 3.18 | 0.78 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38-only graph | 33,296 | 1.43 | 3.90 | 4.87 | 0.00 | 98.03 | 38.70 | 0.10 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38-only graph | 32,562 | 1.48 | 3.77 | 4.88 | 0.00 | 98.72 | 40.05 | 0.10 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs HPRC v1.1 d9 | 46,680 | 0.49 | 0.38 | 0.27 | 0.13 | 98.54 | 24.18 | 0.04 |
| outside HG008T benchmark BED | GRCh38-only graph vs GRCh38 BAM | 35,592 | 0.94 | 8.78 | 0.00 | 0.00 | 99.01 | 31.90 | 1.07 |
| outside HG008T benchmark BED | HPRC v1.1 d9 vs GRCh38 BAM | 29,267 | 2.27 | 12.34 | 5.43 | 0.00 | 98.38 | 43.85 | 1.03 |
| outside HG008T benchmark BED | HPRC v2.1 d46 vs GRCh38 BAM | 29,472 | 2.17 | 12.23 | 5.42 | 0.00 | 98.37 | 43.44 | 1.05 |

## COLO829T_fiberseq

### G. Genome-wide name-hash sample (crc32(name) % 100 < 10), graph arms

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 27,846,735 | 1,564,272 | 5.62 | 0.05 | 2.17 | 92.38 | 2,784,787 |
| HPRC v1.1 d9 | 27,846,735 | 1,510,899 | 5.43 | 0.10 | 2.48 | 92.63 | 2,784,787 |
| HPRC v2.1 d46 | 27,846,735 | 8,094 | 0.03 | 0.06 | 7.61 | 76.87 | 2,784,787 |

Sampled reads (union over graph arms): 2,687,506; absent from an arm (counted unmapped there): GRCh38-only graph 0, HPRC v1.1 d9 0, HPRC v2.1 d46 0.

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph (unmapped 3.37 %) | 2,687,506 | 0.18 | 95.81 | 99.31 | 3.26 | 0.04 | 2.48 | 0.86 | 3.50 | 93.12 | 0.9970 | 1.80 | 1.82 | 0.30 | - | - | - | 0.00 |
| HPRC v1.1 d9 (unmapped 3.23 %) | 2,687,506 | 1.50 | 74.69 | 97.12 | 2.58 | 0.08 | 2.71 | 0.84 | 3.03 | 93.34 | 0.9989 | 0.84 | 1.62 | 0.24 | 92.54 | 8.12 | 97.90 | 0.00 |
| HPRC v2.1 d46 (unmapped 0.01 %) | 2,687,506 | 1.54 | 74.05 | 97.23 | 4.96 | 0.04 | 6.34 | 1.17 | 13.14 | 79.31 | 0.9989 | 2.19 | 1.87 | 1.13 | 90.23 | 8.56 | 92.17 | 0.00 |

Same reads, B relative to A:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 2,687,506 | 2,594,452 | 2.42 | 1.36 | 0.25 | 0.44 | 0.18 | 1.50 | 1.32 | 0.00 | 3.15 | 2.50 | 98.72 | 99.97 | 747 |
| HPRC v2.1 d46 vs GRCh38-only graph | 2,687,506 | 2,596,786 | 1.18 | 16.06 | 0.16 | 1.15 | 0.18 | 1.54 | 1.36 | 0.00 | 3.15 | 4.96 | 98.29 | 99.97 | 757 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,687,506 | 2,600,479 | 1.20 | 15.83 | 0.14 | 1.01 | 1.50 | 1.54 | 0.17 | 0.13 | 2.50 | 4.96 | 98.50 | 99.96 | 912 |

Reads that became perfect / stopped being perfect:

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 35,496 | 96.08 | 11.0 | 16.06 | 83.85 | 0.45 | 0.04 | 9.0 | 99.56 | 6.95 | 0.41 | 3 | 100.00 |
| HPRC v2.1 d46 vs GRCh38-only graph | 36,452 | 96.23 | 11.0 | 16.23 | 83.69 | 0.47 | 0.04 | 9.0 | 99.57 | 0.49 | 53.51 | 4 | 100.00 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 4,446 | 97.59 | 15.0 | 48.31 | 51.60 | 0.09 | 0.07 | 1.0 | 99.19 | 0.52 | 49.92 | 3,491 | 90.17 |

By region (anchor: linear 2,596,923, d9 5,831, d46 80,838, none 3,914):

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 2,411,244 | GRCh38-only graph | 0.01 | 0.19 | 95.63 | 99.31 | 2.34 | 0.00 | 96.78 | 0.9971 | - |
| SMaHT easy | 2,411,244 | HPRC v1.1 d9 | 0.00 | 1.60 | 73.71 | 97.08 | 1.76 | 0.00 | 97.22 | 0.9990 | 94.07 |
| SMaHT easy | 2,411,244 | HPRC v2.1 d46 | 0.00 | 1.64 | 72.14 | 97.12 | 1.77 | 0.00 | 86.15 | 0.9990 | 93.73 |
| SMaHT difficult | 9,660 | GRCh38-only graph | 9.49 | 0.04 | 97.47 | 98.68 | 8.81 | 0.01 | 89.87 | 0.9957 | - |
| SMaHT difficult | 9,660 | HPRC v1.1 d9 | 9.03 | 1.33 | 79.27 | 95.88 | 3.24 | 0.00 | 94.24 | 0.9989 | 95.91 |
| SMaHT difficult | 9,660 | HPRC v2.1 d46 | 0.00 | 1.31 | 80.48 | 96.07 | 11.60 | 0.01 | 72.29 | 0.9987 | 95.87 |
| SMaHT extreme | 229,856 | GRCh38-only graph | 32.76 | 0.08 | 97.92 | 99.32 | 12.67 | 0.59 | 45.32 | 0.9942 | - |
| SMaHT extreme | 229,856 | HPRC v1.1 d9 | 31.77 | 0.68 | 85.84 | 97.66 | 10.35 | 1.22 | 43.49 | 0.9977 | 77.01 |
| SMaHT extreme | 229,856 | HPRC v2.1 d46 | 0.04 | 0.70 | 89.94 | 98.14 | 29.72 | 0.51 | 19.21 | 0.9922 | 62.21 |
| segdup | 164,038 | GRCh38-only graph | 14.31 | 0.10 | 97.48 | 99.43 | 8.02 | 0.22 | 57.60 | 0.9958 | - |
| segdup | 164,038 | HPRC v1.1 d9 | 14.45 | 0.92 | 82.14 | 97.67 | 5.68 | 0.87 | 55.99 | 0.9985 | 90.70 |
| segdup | 164,038 | HPRC v2.1 d46 | 0.01 | 0.93 | 84.03 | 98.01 | 18.10 | 0.17 | 35.47 | 0.9980 | 85.38 |
| not segdup | 2,486,722 | GRCh38-only graph | 2.13 | 0.19 | 95.68 | 99.30 | 2.69 | 0.02 | 95.75 | 0.9970 | - |
| not segdup | 2,486,722 | HPRC v1.1 d9 | 2.02 | 1.56 | 74.02 | 97.07 | 2.09 | 0.03 | 96.13 | 0.9989 | 93.18 |
| not segdup | 2,486,722 | HPRC v2.1 d46 | 0.00 | 1.60 | 73.03 | 97.15 | 3.31 | 0.04 | 83.25 | 0.9989 | 91.38 |

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.03 | 0.04 | 1.41 | 0.00 | 99.99 | 13 | 96.03 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.03 | 0.13 | 1.45 | 0.00 | 99.98 | 11 | 96.18 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.01 | 0.11 | 0.17 | 0.14 | 99.98 | 15 | 97.58 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 1.57 | 0.47 | 1.28 | 0.00 | 99.61 | 8 | 96.77 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 1.53 | 1.30 | 1.27 | 0.00 | 99.33 | 9 | 96.75 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.24 | 1.15 | 0.04 | 0.05 | 99.50 | 4 | 100.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 3.40 | 6.26 | 0.60 | 0.00 | 82.11 | 579 | 97.16 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 1.77 | 14.20 | 0.62 | 0.00 | 76.99 | 563 | 97.32 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.78 | 12.16 | 0.10 | 0.08 | 80.27 | 528 | 97.77 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 2.49 | 5.22 | 0.82 | 0.00 | 86.82 | 323 | 96.86 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1.57 | 11.90 | 0.83 | 0.00 | 84.06 | 402 | 96.84 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.66 | 9.18 | 0.13 | 0.12 | 88.22 | 343 | 98.58 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 0.10 | 0.13 | 1.37 | 0.00 | 99.62 | 277 | 96.05 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 0.05 | 0.34 | 1.41 | 0.00 | 99.44 | 181 | 96.20 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.03 | 0.36 | 0.17 | 0.13 | 99.39 | 204 | 97.54 |

### W. Window sample (148 GRCh38 windows, flank 25,000 bp; reads placed in a window by at least one arm: 62,786; random 12,286, random_extra 46,047, segdup 4,453)

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 97.65 | 2.35 | 0.00 |
| GRCh38-only graph | genome | 97.13 | 2.53 | 0.33 |
| HPRC v1.1 d9 | genome | 98.05 | 1.68 | 0.27 |
| HPRC v2.1 d46 | genome | 97.61 | 2.39 | 0.00 |

Unmapped records of the whole file (GAM head / BAM unplaced): GRCh38 BAM 248,616

Reads placed in the window by the arm:

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 61,311 | 0.18 | 96.59 | 99.39 | 8.76 | 0.01 | 1.40 | 0.66 | 0.74 | 97.19 | 0.9969 | 1.44 | 1.79 | 0.30 | - | - | - | 0.00 |
| GRCh38-only graph | 60,986 | 0.18 | 96.42 | 99.40 | 2.79 | 0.03 | 2.63 | 0.63 | 3.36 | 93.36 | 0.9969 | 1.48 | 1.82 | 0.19 | - | - | - | 0.00 |
| HPRC v1.1 d9 | 61,561 | 1.58 | 74.27 | 97.16 | 1.96 | 0.17 | 3.49 | 0.50 | 2.80 | 93.04 | 0.9989 | 0.46 | 1.61 | 0.10 | 94.13 | 7.28 | 98.46 | 0.00 |
| HPRC v2.1 d46 | 61,287 | 1.64 | 71.95 | 97.20 | 2.12 | 0.00 | 3.80 | 0.88 | 13.26 | 82.07 | 0.9990 | 0.45 | 1.63 | 0.18 | 93.71 | 4.85 | 98.46 | 0.00 |

Same reads, both arms place them in the window (B relative to A); last columns: reads A places in the window that B places elsewhere / leaves unmapped:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both  A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 60,488 | 60,488 | 2.17 | 1.18 | 0.24 | 0.38 | 0.18 | 1.59 | 1.41 | 0.00 | 2.68 | 1.91 | 99.79 | 100.00 | 1 | 495 | 3 | 1,052 | 21 |
| HPRC v2.1 d46 vs GRCh38-only graph | 60,312 | 60,312 | 0.87 | 15.14 | 0.09 | 0.91 | 0.18 | 1.65 | 1.47 | 0.00 | 2.63 | 1.91 | 99.75 | 100.00 | 0 | 674 | 0 | 848 | 127 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 60,831 | 60,831 | 1.23 | 14.87 | 0.11 | 0.84 | 1.59 | 1.65 | 0.19 | 0.13 | 1.87 | 1.91 | 99.80 | 100.00 | 0 | 730 | 0 | 345 | 111 |
| GRCh38-only graph vs GRCh38 BAM | 60,720 | 60,720 | 0.56 | 5.67 | 0.09 | 1.69 | 0.18 | 0.18 | 0.00 | 0.00 | 8.28 | 2.68 | 99.71 | 100.00 | 1 | 451 | 140 | 266 | 0 |
| HPRC v1.1 d9 vs GRCh38 BAM | 60,448 | 60,448 | 0.57 | 5.00 | 0.16 | 1.71 | 0.18 | 1.59 | 1.41 | 0.00 | 8.22 | 1.89 | 99.65 | 99.99 | 7 | 758 | 105 | 1,113 | 0 |
| HPRC v2.1 d46 vs GRCh38 BAM | 60,371 | 60,371 | 0.40 | 16.29 | 0.01 | 2.24 | 0.18 | 1.65 | 1.47 | 0.00 | 8.23 | 1.97 | 99.62 | 100.00 | 1 | 940 | 0 | 915 | 1 |

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 851 | 96.12 | 12.0 | 13.75 | 86.25 | 0.35 | 0.00 | 10.0 | 100.00 | 7.76 | 0.82 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 887 | 96.05 | 12.0 | 13.87 | 86.13 | 0.23 | 0.00 | 10.0 | 100.00 | 0.56 | 52.99 | 0 | n/a |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 117 | 96.58 | 16.0 | 43.59 | 56.41 | 0.00 | 0.00 | 1.0 | 100.00 | 1.71 | 49.57 | 81 | 97.53 |
| GRCh38-only graph vs GRCh38 BAM | 0 | | | | | | | | | | | | |
| HPRC v1.1 d9 vs GRCh38 BAM | 851 | 96.12 | 12.0 | 13.63 | 86.02 | 2.35 | 0.00 | 10.0 | 100.00 | 0.00 | 16.57 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38 BAM | 887 | 96.05 | 12.0 | 13.75 | 85.91 | 2.03 | 0.00 | 10.0 | 100.00 | 0.00 | 53.89 | 0 | n/a |

By window kind and region:

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 12,286 | GRCh38 BAM | 97.58 | 2.42 | 0.00 | 0.12 | 96.99 | 99.54 | 8.16 | 0.01 | 97.54 | 0.9969 |
| windows: random | 12,286 | GRCh38-only graph | 97.50 | 2.49 | 0.01 | 0.12 | 96.87 | 99.54 | 2.45 | 0.00 | 93.62 | 0.9969 |
| windows: random | 12,286 | HPRC v1.1 d9 | 98.70 | 1.29 | 0.01 | 1.49 | 74.17 | 97.19 | 1.61 | 0.63 | 93.49 | 0.9990 |
| windows: random | 12,286 | HPRC v2.1 d46 | 96.58 | 3.42 | 0.00 | 1.62 | 71.36 | 97.13 | 1.58 | 0.00 | 84.08 | 0.9990 |
| windows: random_extra | 46,047 | GRCh38 BAM | 99.42 | 0.58 | 0.00 | 0.20 | 96.50 | 99.38 | 8.41 | 0.00 | 98.05 | 0.9969 |
| windows: random_extra | 46,047 | GRCh38-only graph | 98.79 | 0.86 | 0.35 | 0.20 | 96.32 | 99.39 | 2.49 | 0.03 | 95.21 | 0.9970 |
| windows: random_extra | 46,047 | HPRC v1.1 d9 | 98.64 | 1.06 | 0.30 | 1.62 | 73.96 | 97.13 | 1.86 | 0.04 | 95.96 | 0.9989 |
| windows: random_extra | 46,047 | HPRC v2.1 d46 | 98.74 | 1.26 | 0.00 | 1.66 | 71.61 | 97.20 | 2.03 | 0.00 | 84.73 | 0.9990 |
| windows: segdup | 4,453 | GRCh38 BAM | 79.59 | 20.41 | 0.00 | 0.20 | 96.42 | 99.10 | 15.29 | 0.11 | 84.85 | 0.9963 |
| windows: segdup | 4,453 | GRCh38-only graph | 78.98 | 19.99 | 1.03 | 0.20 | 96.19 | 99.06 | 7.85 | 0.14 | 68.44 | 0.9963 |
| windows: segdup | 4,453 | HPRC v1.1 d9 | 90.16 | 9.18 | 0.65 | 1.34 | 78.08 | 97.33 | 4.16 | 0.30 | 58.70 | 0.9987 |
| windows: segdup | 4,453 | HPRC v2.1 d46 | 88.77 | 11.23 | 0.00 | 1.47 | 77.66 | 97.37 | 4.91 | 0.00 | 45.43 | 0.9987 |
| SMaHT easy | 55,625 | GRCh38 BAM | 99.94 | 0.06 | 0.00 | 0.19 | 96.44 | 99.40 | 7.60 | 0.00 | 100.00 | 0.9970 |
| SMaHT easy | 55,625 | GRCh38-only graph | 99.94 | 0.05 | 0.01 | 0.19 | 96.29 | 99.41 | 2.16 | 0.00 | 96.80 | 0.9970 |
| SMaHT easy | 55,625 | HPRC v1.1 d9 | 100.00 | 0.00 | 0.00 | 1.62 | 73.64 | 97.12 | 1.61 | 0.00 | 97.29 | 0.9990 |
| SMaHT easy | 55,625 | HPRC v2.1 d46 | 99.97 | 0.03 | 0.00 | 1.68 | 71.12 | 97.16 | 1.60 | 0.00 | 86.43 | 0.9990 |
| SMaHT difficult | 281 | GRCh38 BAM | 99.64 | 0.36 | 0.00 | 0.36 | 96.07 | 99.64 | 23.57 | 0.00 | 95.71 | 0.9920 |
| SMaHT difficult | 281 | GRCh38-only graph | 92.88 | 0.00 | 7.12 | 0.38 | 95.79 | 99.62 | 7.66 | 0.00 | 93.49 | 0.9930 |
| SMaHT difficult | 281 | HPRC v1.1 d9 | 93.95 | 0.00 | 6.05 | 2.65 | 67.42 | 95.83 | 1.52 | 0.00 | 97.35 | 0.9991 |
| SMaHT difficult | 281 | HPRC v2.1 d46 | 95.73 | 4.27 | 0.00 | 2.60 | 66.54 | 95.91 | 3.72 | 0.00 | 82.16 | 0.9990 |
| SMaHT extreme | 6,880 | GRCh38 BAM | 79.08 | 20.90 | 0.01 | 0.13 | 98.14 | 99.30 | 19.85 | 0.13 | 68.54 | 0.9953 |
| SMaHT extreme | 6,880 | GRCh38-only graph | 74.62 | 22.69 | 2.69 | 0.14 | 97.92 | 99.24 | 9.41 | 0.33 | 56.02 | 0.9958 |
| SMaHT extreme | 6,880 | HPRC v1.1 d9 | 82.46 | 15.35 | 2.19 | 1.16 | 80.73 | 97.55 | 5.41 | 1.85 | 51.21 | 0.9985 |
| SMaHT extreme | 6,880 | HPRC v2.1 d46 | 78.63 | 21.37 | 0.00 | 1.22 | 80.74 | 97.63 | 7.41 | 0.00 | 37.25 | 0.9985 |
| segdup | 6,486 | GRCh38 BAM | 80.97 | 19.03 | 0.00 | 0.19 | 97.60 | 99.26 | 14.64 | 0.13 | 79.09 | 0.9963 |
| segdup | 6,486 | GRCh38-only graph | 79.80 | 19.09 | 1.11 | 0.19 | 97.41 | 99.25 | 7.19 | 0.27 | 64.08 | 0.9964 |
| segdup | 6,486 | HPRC v1.1 d9 | 88.47 | 10.78 | 0.76 | 1.15 | 78.67 | 97.44 | 3.40 | 1.78 | 58.23 | 0.9987 |
| segdup | 6,486 | HPRC v2.1 d46 | 83.49 | 16.51 | 0.00 | 1.14 | 78.15 | 97.53 | 4.04 | 0.00 | 45.48 | 0.9987 |
| not segdup | 56,300 | GRCh38 BAM | 99.57 | 0.43 | 0.00 | 0.18 | 96.50 | 99.41 | 8.21 | 0.00 | 98.88 | 0.9970 |
| not segdup | 56,300 | GRCh38-only graph | 99.13 | 0.63 | 0.25 | 0.18 | 96.33 | 99.41 | 2.38 | 0.01 | 96.07 | 0.9970 |
| not segdup | 56,300 | HPRC v1.1 d9 | 99.15 | 0.64 | 0.21 | 1.62 | 73.82 | 97.13 | 1.81 | 0.01 | 96.62 | 0.9990 |
| not segdup | 56,300 | HPRC v2.1 d46 | 99.24 | 0.76 | 0.00 | 1.69 | 71.35 | 97.17 | 1.94 | 0.00 | 85.62 | 0.9990 |

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | HPRC v1.1 d9 vs GRCh38-only graph | 11,888 | 0.20 | 0.05 | 1.39 | 0.00 | 99.87 | 0.76 | 0.00 |
| windows: random | HPRC v2.1 d46 vs GRCh38-only graph | 11,781 | 0.01 | 0.07 | 1.50 | 0.00 | 99.99 | 1.65 | 0.00 |
| windows: random | HPRC v2.1 d46 vs HPRC v1.1 d9 | 11,824 | 0.00 | 0.12 | 0.21 | 0.11 | 99.92 | 2.49 | 0.00 |
| windows: random | GRCh38-only graph vs GRCh38 BAM | 11,926 | 0.01 | 0.72 | 0.00 | 0.00 | 99.82 | 0.52 | 0.01 |
| windows: random | HPRC v1.1 d9 vs GRCh38 BAM | 11,876 | 0.07 | 0.63 | 1.39 | 0.00 | 99.77 | 0.93 | 0.01 |
| windows: random | HPRC v2.1 d46 vs GRCh38 BAM | 11,760 | 0.00 | 0.39 | 1.51 | 0.00 | 99.88 | 1.91 | 0.00 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38-only graph | 45,296 | 0.20 | 0.11 | 1.43 | 0.00 | 99.78 | 0.42 | 0.01 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38-only graph | 45,261 | 0.06 | 0.52 | 1.47 | 0.00 | 99.73 | 0.50 | 0.00 |
| windows: random_extra | HPRC v2.1 d46 vs HPRC v1.1 d9 | 45,287 | 0.01 | 0.55 | 0.19 | 0.15 | 99.81 | 0.29 | 0.00 |
| windows: random_extra | GRCh38-only graph vs GRCh38 BAM | 45,384 | 0.10 | 0.78 | 0.00 | 0.00 | 99.67 | 0.63 | 0.23 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38 BAM | 45,292 | 0.13 | 0.62 | 1.43 | 0.00 | 99.62 | 0.87 | 0.19 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38 BAM | 45,332 | 0.01 | 1.03 | 1.47 | 0.00 | 99.59 | 0.97 | 0.00 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 3,304 | 1.03 | 5.24 | 1.18 | 0.00 | 99.55 | 6.06 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 3,270 | 0.70 | 9.45 | 1.31 | 0.00 | 99.08 | 7.02 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 3,720 | 1.59 | 6.67 | 0.13 | 0.03 | 99.30 | 7.35 | 0.00 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 3,410 | 0.21 | 17.30 | 0.00 | 0.00 | 99.91 | 2.88 | 0.90 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 3,280 | 0.88 | 20.76 | 1.19 | 0.00 | 99.51 | 6.94 | 0.51 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 3,279 | 0.09 | 25.53 | 1.31 | 0.00 | 99.12 | 7.48 | 0.00 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 55,591 | 0.03 | 0.00 | 1.43 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 55,587 | 0.03 | 0.05 | 1.49 | 0.00 | 100.00 | 0.01 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 55,607 | 0.00 | 0.05 | 0.20 | 0.13 | 100.00 | 0.03 | 0.00 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 55,585 | 0.00 | 0.09 | 0.00 | 0.00 | 100.00 | 0.00 | 0.01 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 55,590 | 0.00 | 0.07 | 1.43 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 55,588 | 0.00 | 0.11 | 1.49 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 261 | 1.53 | 0.00 | 2.30 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 261 | 1.53 | 0.38 | 2.30 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 264 | 0.00 | 0.38 | 0.00 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT difficult | GRCh38-only graph vs GRCh38 BAM | 261 | 0.00 | 1.53 | 0.00 | 0.00 | 100.00 | 0.00 | 6.79 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38 BAM | 264 | 0.00 | 0.00 | 2.27 | 0.00 | 100.00 | 0.00 | 5.71 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38 BAM | 268 | 0.00 | 1.87 | 2.24 | 0.00 | 100.00 | 4.29 | 0.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 4,636 | 2.74 | 4.94 | 1.06 | 0.00 | 97.20 | 9.64 | 0.06 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 4,464 | 0.69 | 11.65 | 1.12 | 0.00 | 96.57 | 13.05 | 0.00 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 4,960 | 1.25 | 9.72 | 0.14 | 0.12 | 97.56 | 12.57 | 0.00 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 4,874 | 1.07 | 19.96 | 0.00 | 0.00 | 96.39 | 8.29 | 2.13 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 4,594 | 2.05 | 21.70 | 1.07 | 0.00 | 95.36 | 13.93 | 1.64 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 4,515 | 0.13 | 28.46 | 1.11 | 0.00 | 94.88 | 17.02 | 0.00 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 4,760 | 1.83 | 4.29 | 0.97 | 0.00 | 98.55 | 8.04 | 0.00 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 4,616 | 0.71 | 9.92 | 0.93 | 0.00 | 98.53 | 10.82 | 0.00 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 5,119 | 1.19 | 7.93 | 0.12 | 0.18 | 99.24 | 10.79 | 0.00 |
| segdup | GRCh38-only graph vs GRCh38 BAM | 5,008 | 0.18 | 18.41 | 0.00 | 0.00 | 99.42 | 3.62 | 1.03 |
| segdup | HPRC v1.1 d9 vs GRCh38 BAM | 4,744 | 1.05 | 20.11 | 0.97 | 0.00 | 98.08 | 8.99 | 0.69 |
| segdup | HPRC v2.1 d46 vs GRCh38 BAM | 4,626 | 0.09 | 26.11 | 0.93 | 0.00 | 98.05 | 11.92 | 0.00 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 55,728 | 0.11 | 0.05 | 1.44 | 0.00 | 99.89 | 0.14 | 0.01 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 55,696 | 0.04 | 0.17 | 1.52 | 0.00 | 99.85 | 0.20 | 0.00 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 55,712 | 0.01 | 0.19 | 0.20 | 0.13 | 99.85 | 0.20 | 0.00 |
| not segdup | GRCh38-only graph vs GRCh38 BAM | 55,712 | 0.08 | 0.19 | 0.00 | 0.00 | 99.74 | 0.47 | 0.15 |
| not segdup | HPRC v1.1 d9 vs GRCh38 BAM | 55,704 | 0.08 | 0.15 | 1.45 | 0.00 | 99.78 | 0.51 | 0.12 |
| not segdup | HPRC v2.1 d46 vs GRCh38 BAM | 55,745 | 0.00 | 0.26 | 1.51 | 0.00 | 99.75 | 0.56 | 0.00 |

## COLO829T_ONT

### G. Genome-wide name-hash sample (crc32(name) % 100 < 10), graph arms

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 14,448,438 | 1,954,638 | 13.53 | 0.88 | 2.72 | 93.85 | 1,444,965 |
| HPRC v1.1 d9 | 14,448,438 | 1,911,847 | 13.23 | 1.10 | 3.19 | 93.46 | 1,444,965 |
| HPRC v2.1 d46 | 14,448,438 | 1,968,240 | 13.62 | 0.38 | 2.80 | 89.57 | 1,444,965 |

Sampled reads (union over graph arms): 1,444,965; absent from an arm (counted unmapped there): GRCh38-only graph 0, HPRC v1.1 d9 0, HPRC v2.1 d46 0.

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph (unmapped 13.53 %) | 1,444,965 | 0.00 | 92.50 | 93.27 | 99.03 | 0.88 | 2.44 | 1.48 | 1.33 | 93.86 | 0.9887 | 9.16 | 8.59 | 2.82 | - | - | - | 0.00 |
| HPRC v1.1 d9 (unmapped 13.23 %) | 1,444,965 | 0.01 | 89.77 | 92.82 | 99.03 | 1.10 | 2.73 | 1.42 | 1.26 | 93.48 | 0.9910 | 8.36 | 8.55 | 2.75 | 73.15 | 8.58 | 94.47 | 0.00 |
| HPRC v2.1 d46 (unmapped 13.61 %) | 1,444,965 | 0.01 | 90.36 | 93.48 | 99.04 | 0.38 | 3.51 | 3.13 | 3.43 | 89.54 | 0.9908 | 8.92 | 9.12 | 2.85 | 73.51 | 9.37 | 90.10 | 0.00 |

Same reads, B relative to A:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,444,965 | 1,246,241 | 2.02 | 2.81 | 0.54 | 0.96 | 0.00 | 0.01 | 0.00 | 0.00 | 85.63 | 85.92 | 96.80 | 99.63 | 4,356 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,444,965 | 1,226,010 | 2.75 | 7.57 | 1.01 | 3.63 | 0.00 | 0.01 | 0.00 | 0.00 | 85.63 | 85.56 | 95.01 | 98.57 | 16,164 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,444,965 | 1,229,957 | 3.16 | 7.05 | 1.14 | 3.33 | 0.01 | 0.01 | 0.00 | 0.00 | 85.92 | 85.56 | 95.20 | 98.60 | 15,892 |

Reads that became perfect / stopped being perfect:

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 39 | 87.18 | 1.0 | 64.10 | 35.90 | 0.00 | 0.00 | 1.0 | 97.44 | 7.69 | 0.00 | 9 | 22.22 |
| HPRC v2.1 d46 vs GRCh38-only graph | 40 | 87.50 | 1.0 | 65.00 | 35.00 | 2.50 | 0.00 | 1.0 | 97.50 | 2.50 | 2.50 | 8 | 12.50 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 4 | 50.00 | 2.0 | 50.00 | 25.00 | 25.00 | 0.00 | 1.0 | 50.00 | 0.00 | 25.00 | 2 | 50.00 |

By region (anchor: linear 1,249,488, d9 7,096, d46 14,331, none 174,050):

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 1,085,973 | GRCh38-only graph | 0.73 | 0.00 | 93.33 | 94.75 | 99.04 | 0.02 | 99.16 | 0.9891 | - |
| SMaHT easy | 1,085,973 | HPRC v1.1 d9 | 0.43 | 0.01 | 90.72 | 94.34 | 99.04 | 0.03 | 98.94 | 0.9913 | 75.88 |
| SMaHT easy | 1,085,973 | HPRC v2.1 d46 | 0.75 | 0.01 | 90.70 | 94.31 | 99.04 | 0.02 | 95.88 | 0.9913 | 76.01 |
| SMaHT difficult | 33,374 | GRCh38-only graph | 2.06 | 0.01 | 75.42 | 74.16 | 98.95 | 6.29 | 88.82 | 0.9901 | - |
| SMaHT difficult | 33,374 | HPRC v1.1 d9 | 1.27 | 0.03 | 67.43 | 71.89 | 98.95 | 7.16 | 88.06 | 0.9928 | 49.03 |
| SMaHT difficult | 33,374 | HPRC v2.1 d46 | 9.86 | 0.03 | 69.54 | 75.89 | 98.93 | 1.24 | 74.60 | 0.9923 | 51.48 |
| SMaHT extreme | 129,271 | GRCh38-only graph | 9.20 | 0.01 | 89.07 | 85.18 | 98.89 | 7.08 | 53.44 | 0.9837 | - |
| SMaHT extreme | 129,271 | HPRC v1.1 d9 | 9.05 | 0.01 | 86.47 | 84.84 | 98.91 | 8.94 | 51.94 | 0.9872 | 60.18 |
| SMaHT extreme | 129,271 | HPRC v2.1 d46 | 9.53 | 0.01 | 91.40 | 90.35 | 98.96 | 2.90 | 43.54 | 0.9823 | 60.42 |
| segdup | 81,169 | GRCh38-only graph | 6.25 | 0.00 | 93.60 | 92.52 | 98.84 | 1.41 | 58.05 | 0.9852 | - |
| segdup | 81,169 | HPRC v1.1 d9 | 5.70 | 0.00 | 91.29 | 92.04 | 98.83 | 2.81 | 55.80 | 0.9887 | 73.38 |
| segdup | 81,169 | HPRC v2.1 d46 | 2.56 | 0.00 | 93.54 | 93.98 | 98.93 | 1.89 | 49.42 | 0.9866 | 73.09 |
| not segdup | 1,167,449 | GRCh38-only graph | 1.32 | 0.01 | 92.37 | 93.34 | 99.04 | 0.82 | 96.92 | 0.9890 | - |
| not segdup | 1,167,449 | HPRC v1.1 d9 | 1.04 | 0.01 | 89.58 | 92.88 | 99.03 | 0.96 | 96.71 | 0.9912 | 73.68 |
| not segdup | 1,167,449 | HPRC v2.1 d46 | 1.86 | 0.01 | 90.02 | 93.44 | 99.04 | 0.21 | 93.18 | 0.9912 | 73.97 |

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.05 | 0.18 | 0.00 | 0.00 | 99.93 | 107 | 92.86 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.23 | 1.82 | 0.00 | 0.00 | 99.08 | 4,736 | 89.66 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.34 | 1.78 | 0.00 | 0.00 | 99.08 | 4,677 | 50.00 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 0.96 | 1.81 | 0.02 | 0.00 | 91.64 | 101 | 66.67 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 1.23 | 13.09 | 0.02 | 0.00 | 87.12 | 1,579 | 66.67 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.50 | 12.82 | 0.00 | 0.00 | 87.61 | 1,563 | n/a |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 4.63 | 7.41 | 0.00 | 0.01 | 72.45 | 3,655 | 100.00 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 7.51 | 17.30 | 0.00 | 0.01 | 62.82 | 7,753 | 100.00 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7.72 | 14.63 | 0.00 | 0.00 | 64.75 | 7,399 | 50.00 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 3.58 | 7.14 | 0.00 | 0.00 | 80.25 | 1,367 | 100.00 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 6.29 | 13.77 | 0.00 | 0.00 | 66.30 | 6,240 | 100.00 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 6.80 | 10.79 | 0.00 | 0.00 | 67.45 | 5,994 | 0.00 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 0.30 | 0.49 | 0.00 | 0.00 | 98.25 | 2,496 | 89.19 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 0.54 | 2.78 | 0.00 | 0.00 | 97.53 | 7,828 | 87.18 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.64 | 2.68 | 0.00 | 0.00 | 97.65 | 7,645 | 66.67 |

### W. Window sample (148 GRCh38 windows, flank 25,000 bp; reads placed in a window by at least one arm: 33,011; random 6,224, random_extra 24,248, segdup 2,539)

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 92.45 | 6.31 | 1.24 |
| GRCh38-only graph | genome | 87.84 | 6.84 | 5.32 |
| HPRC v1.1 d9 | genome | 88.75 | 6.16 | 5.10 |
| HPRC v2.1 d46 | genome | 88.53 | 4.51 | 6.96 |

Unmapped records of the whole file (GAM head / BAM unplaced): GRCh38 BAM 1,070,571

Reads placed in the window by the arm:

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 30,519 | 0.01 | 93.61 | 94.77 | 99.63 | 0.67 | 2.22 | 1.53 | 2.15 | 93.42 | 0.9884 | 12.94 | 12.16 | 3.76 | - | - | - | 0.00 |
| GRCh38-only graph | 28,997 | 0.01 | 92.85 | 93.79 | 99.07 | 0.71 | 2.17 | 1.22 | 1.31 | 94.59 | 0.9889 | 8.84 | 8.99 | 2.59 | - | - | - | 0.00 |
| HPRC v1.1 d9 | 29,296 | 0.01 | 90.04 | 93.19 | 99.07 | 1.13 | 2.72 | 1.32 | 1.22 | 93.60 | 0.9913 | 8.04 | 8.93 | 2.28 | 76.33 | 7.93 | 94.73 | 0.00 |
| HPRC v2.1 d46 | 29,226 | 0.01 | 90.69 | 93.78 | 99.08 | 0.22 | 3.25 | 2.82 | 3.22 | 90.50 | 0.9910 | 8.59 | 9.34 | 2.24 | 75.66 | 5.71 | 94.11 | 0.00 |

Same reads, both arms place them in the window (B relative to A); last columns: reads A places in the window that B places elsewhere / leaves unmapped:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both  A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 28,388 | 28,388 | 1.60 | 1.80 | 0.32 | 0.57 | 0.01 | 0.01 | 0.00 | 0.00 | 99.06 | 99.07 | 99.54 | 99.96 | 12 | 593 | 16 | 824 | 84 |
| HPRC v2.1 d46 vs GRCh38-only graph | 27,550 | 27,550 | 1.71 | 4.92 | 0.51 | 2.15 | 0.01 | 0.01 | 0.00 | 0.00 | 99.05 | 99.06 | 99.38 | 99.88 | 33 | 977 | 470 | 1,473 | 203 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 27,759 | 27,759 | 2.28 | 4.81 | 0.73 | 2.08 | 0.01 | 0.01 | 0.00 | 0.00 | 99.05 | 99.05 | 99.42 | 99.86 | 36 | 1,020 | 517 | 1,335 | 132 |
| GRCh38-only graph vs GRCh38 BAM | 28,417 | 28,417 | 2.40 | 2.85 | 0.62 | 1.06 | 0.01 | 0.01 | 0.00 | 0.00 | 99.62 | 99.06 | 99.55 | 99.99 | 4 | 377 | 1,725 | 372 | 208 |
| HPRC v1.1 d9 vs GRCh38 BAM | 28,355 | 28,355 | 2.30 | 2.83 | 0.64 | 1.31 | 0.01 | 0.01 | 0.00 | 0.00 | 99.62 | 99.06 | 99.48 | 99.97 | 7 | 517 | 1,647 | 704 | 237 |
| HPRC v2.1 d46 vs GRCh38 BAM | 27,662 | 27,662 | 1.91 | 4.99 | 0.47 | 2.36 | 0.01 | 0.01 | 0.00 | 0.00 | 99.62 | 99.05 | 99.33 | 99.86 | 36 | 1,017 | 1,840 | 1,513 | 51 |

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 0 | | | | | | | | | | | | |
| GRCh38-only graph vs GRCh38 BAM | 0 | | | | | | | | | | | | |
| HPRC v1.1 d9 vs GRCh38 BAM | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38 BAM | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |

By window kind and region:

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 6,224 | GRCh38 BAM | 94.33 | 4.37 | 1.30 | 0.00 | 93.32 | 94.38 | 99.64 | 0.53 | 94.77 | 0.9888 |
| windows: random | 6,224 | GRCh38-only graph | 89.88 | 5.03 | 5.09 | 0.00 | 92.60 | 93.76 | 98.95 | 0.54 | 95.46 | 0.9891 |
| windows: random | 6,224 | HPRC v1.1 d9 | 91.07 | 4.26 | 4.68 | 0.00 | 89.61 | 93.08 | 98.98 | 1.39 | 94.72 | 0.9916 |
| windows: random | 6,224 | HPRC v2.1 d46 | 89.03 | 4.50 | 6.47 | 0.00 | 90.15 | 93.61 | 98.95 | 0.11 | 92.53 | 0.9915 |
| windows: random_extra | 24,248 | GRCh38 BAM | 94.04 | 4.69 | 1.27 | 0.01 | 93.67 | 94.82 | 99.61 | 0.45 | 94.79 | 0.9885 |
| windows: random_extra | 24,248 | GRCh38-only graph | 89.31 | 5.18 | 5.51 | 0.01 | 92.87 | 93.75 | 99.06 | 0.76 | 96.23 | 0.9890 |
| windows: random_extra | 24,248 | HPRC v1.1 d9 | 89.54 | 5.12 | 5.34 | 0.02 | 89.97 | 93.13 | 99.06 | 1.02 | 96.01 | 0.9913 |
| windows: random_extra | 24,248 | HPRC v2.1 d46 | 90.02 | 2.70 | 7.28 | 0.02 | 90.59 | 93.75 | 99.10 | 0.17 | 92.58 | 0.9910 |
| windows: segdup | 2,539 | GRCh38 BAM | 72.71 | 26.51 | 0.79 | 0.00 | 93.82 | 95.40 | 99.89 | 3.90 | 72.16 | 0.9861 |
| windows: segdup | 2,539 | GRCh38-only graph | 68.81 | 27.14 | 4.06 | 0.00 | 93.36 | 94.33 | 99.60 | 0.57 | 71.55 | 0.9870 |
| windows: segdup | 2,539 | HPRC v1.1 d9 | 75.50 | 20.72 | 3.78 | 0.00 | 92.07 | 94.16 | 99.48 | 1.62 | 63.07 | 0.9901 |
| windows: segdup | 2,539 | HPRC v2.1 d46 | 73.18 | 21.82 | 5.00 | 0.00 | 93.43 | 94.67 | 99.25 | 1.02 | 59.96 | 0.9896 |
| SMaHT easy | 26,908 | GRCh38 BAM | 98.90 | 0.82 | 0.29 | 0.00 | 94.12 | 95.48 | 99.61 | 0.00 | 98.34 | 0.9886 |
| SMaHT easy | 26,908 | GRCh38-only graph | 93.68 | 1.11 | 5.21 | 0.00 | 93.76 | 95.13 | 99.05 | 0.03 | 99.18 | 0.9891 |
| SMaHT easy | 26,908 | HPRC v1.1 d9 | 94.26 | 0.80 | 4.94 | 0.01 | 90.98 | 94.41 | 99.05 | 0.26 | 98.71 | 0.9913 |
| SMaHT easy | 26,908 | HPRC v2.1 d46 | 93.48 | 0.88 | 5.64 | 0.01 | 91.09 | 94.65 | 99.05 | 0.01 | 96.23 | 0.9913 |
| SMaHT difficult | 1,006 | GRCh38 BAM | 70.97 | 18.89 | 10.14 | 0.00 | 75.49 | 78.85 | 99.44 | 0.00 | 82.77 | 0.9907 |
| SMaHT difficult | 1,006 | GRCh38-only graph | 76.94 | 20.97 | 2.09 | 0.00 | 73.51 | 75.45 | 98.71 | 6.59 | 88.11 | 0.9903 |
| SMaHT difficult | 1,006 | HPRC v1.1 d9 | 76.94 | 21.27 | 1.79 | 0.00 | 66.41 | 73.26 | 98.84 | 7.75 | 88.24 | 0.9929 |
| SMaHT difficult | 1,006 | HPRC v2.1 d46 | 75.94 | 9.05 | 15.01 | 0.00 | 69.63 | 75.00 | 98.95 | 0.39 | 72.64 | 0.9919 |
| SMaHT extreme | 5,097 | GRCh38 BAM | 62.66 | 32.80 | 4.53 | 0.06 | 93.39 | 92.39 | 99.84 | 6.45 | 54.82 | 0.9841 |
| SMaHT extreme | 5,097 | GRCh38-only graph | 59.17 | 34.26 | 6.57 | 0.07 | 90.19 | 87.23 | 99.34 | 4.87 | 57.89 | 0.9869 |
| SMaHT extreme | 5,097 | HPRC v1.1 d9 | 61.98 | 31.47 | 6.55 | 0.06 | 88.26 | 88.26 | 99.27 | 6.49 | 53.94 | 0.9902 |
| SMaHT extreme | 5,097 | HPRC v2.1 d46 | 64.90 | 22.78 | 12.32 | 0.06 | 92.47 | 91.51 | 99.33 | 1.72 | 51.06 | 0.9863 |
| segdup | 4,090 | GRCh38 BAM | 64.74 | 34.35 | 0.90 | 0.00 | 93.58 | 94.52 | 99.92 | 5.70 | 65.86 | 0.9860 |
| segdup | 4,090 | GRCh38-only graph | 61.42 | 34.89 | 3.69 | 0.00 | 92.99 | 93.19 | 99.40 | 1.00 | 61.86 | 0.9871 |
| segdup | 4,090 | HPRC v1.1 d9 | 65.99 | 30.59 | 3.42 | 0.00 | 91.14 | 93.44 | 99.41 | 3.45 | 56.69 | 0.9902 |
| segdup | 4,090 | HPRC v2.1 d46 | 70.73 | 24.21 | 5.06 | 0.00 | 94.23 | 94.47 | 99.38 | 0.93 | 54.68 | 0.9868 |
| not segdup | 28,921 | GRCh38 BAM | 96.37 | 2.34 | 1.29 | 0.01 | 93.61 | 94.79 | 99.61 | 0.20 | 96.04 | 0.9886 |
| not segdup | 28,921 | GRCh38-only graph | 91.58 | 2.87 | 5.55 | 0.01 | 92.83 | 93.84 | 99.04 | 0.68 | 97.70 | 0.9891 |
| not segdup | 28,921 | HPRC v1.1 d9 | 91.96 | 2.70 | 5.33 | 0.02 | 89.92 | 93.16 | 99.04 | 0.90 | 97.35 | 0.9913 |
| not segdup | 28,921 | HPRC v2.1 d46 | 91.05 | 1.73 | 7.22 | 0.02 | 90.30 | 93.70 | 99.05 | 0.14 | 94.43 | 0.9913 |

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | HPRC v1.1 d9 vs GRCh38-only graph | 5,486 | 0.42 | 0.29 | 0.00 | 0.00 | 99.87 | 1.89 | 0.04 |
| windows: random | HPRC v2.1 d46 vs GRCh38-only graph | 5,355 | 0.41 | 1.90 | 0.00 | 0.00 | 99.79 | 2.84 | 1.43 |
| windows: random | HPRC v2.1 d46 vs HPRC v1.1 d9 | 5,393 | 0.32 | 1.97 | 0.00 | 0.00 | 99.81 | 3.07 | 1.78 |
| windows: random | GRCh38-only graph vs GRCh38 BAM | 5,489 | 0.51 | 0.73 | 0.00 | 0.00 | 99.80 | 1.21 | 5.30 |
| windows: random | HPRC v1.1 d9 vs GRCh38 BAM | 5,475 | 0.64 | 0.62 | 0.00 | 0.00 | 99.87 | 1.86 | 4.89 |
| windows: random | HPRC v2.1 d46 vs GRCh38 BAM | 5,373 | 0.32 | 1.99 | 0.00 | 0.00 | 99.68 | 3.05 | 5.43 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38-only graph | 21,335 | 0.21 | 0.29 | 0.00 | 0.00 | 99.55 | 1.42 | 0.06 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38-only graph | 20,816 | 0.39 | 1.88 | 0.00 | 0.00 | 99.49 | 2.20 | 1.68 |
| windows: random_extra | HPRC v2.1 d46 vs HPRC v1.1 d9 | 20,868 | 0.53 | 1.86 | 0.00 | 0.00 | 99.51 | 2.08 | 1.81 |
| windows: random_extra | GRCh38-only graph vs GRCh38 BAM | 21,280 | 0.58 | 0.63 | 0.00 | 0.00 | 99.54 | 0.90 | 5.78 |
| windows: random_extra | HPRC v1.1 d9 vs GRCh38 BAM | 21,312 | 0.61 | 0.77 | 0.00 | 0.00 | 99.43 | 0.97 | 5.56 |
| windows: random_extra | HPRC v2.1 d46 vs GRCh38 BAM | 20,911 | 0.44 | 1.87 | 0.00 | 0.00 | 99.45 | 2.05 | 6.25 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 1,567 | 1.53 | 5.42 | 0.00 | 0.00 | 98.34 | 10.30 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1,379 | 2.76 | 7.18 | 0.00 | 0.00 | 96.23 | 19.58 | 1.49 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,498 | 5.01 | 5.41 | 0.00 | 0.00 | 96.80 | 20.61 | 1.25 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 1,648 | 1.58 | 7.65 | 0.00 | 0.00 | 98.79 | 5.47 | 5.25 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 1,568 | 1.15 | 11.03 | 0.00 | 0.00 | 98.72 | 10.08 | 4.98 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 1,378 | 1.67 | 11.10 | 0.00 | 0.00 | 96.23 | 20.10 | 5.25 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 25,190 | 0.04 | 0.12 | 0.00 | 0.00 | 100.00 | 0.06 | 0.01 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 24,811 | 0.19 | 1.39 | 0.00 | 0.00 | 99.94 | 0.79 | 0.78 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 24,884 | 0.28 | 1.39 | 0.00 | 0.00 | 99.94 | 0.86 | 1.03 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 25,171 | 0.06 | 0.35 | 0.00 | 0.00 | 99.98 | 0.17 | 5.24 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 25,237 | 0.07 | 0.48 | 0.00 | 0.00 | 99.99 | 0.17 | 4.99 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 24,954 | 0.06 | 1.51 | 0.00 | 0.00 | 99.94 | 0.79 | 5.43 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 709 | 0.28 | 1.27 | 0.00 | 0.00 | 99.86 | 8.40 | 0.00 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 622 | 0.64 | 5.31 | 0.00 | 0.00 | 99.84 | 9.69 | 9.95 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 626 | 1.28 | 5.59 | 0.00 | 0.00 | 99.84 | 8.91 | 10.21 |
| SMaHT difficult | GRCh38-only graph vs GRCh38 BAM | 679 | 4.86 | 0.59 | 0.00 | 0.00 | 100.00 | 2.52 | 2.38 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38 BAM | 683 | 4.98 | 1.17 | 0.00 | 0.00 | 100.00 | 2.38 | 1.96 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38 BAM | 611 | 3.93 | 4.91 | 0.00 | 0.00 | 99.84 | 9.24 | 5.18 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 2,489 | 3.13 | 4.90 | 0.00 | 0.00 | 94.82 | 17.01 | 0.46 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 2,117 | 4.20 | 10.06 | 0.00 | 0.00 | 92.73 | 23.28 | 6.53 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,249 | 5.56 | 8.71 | 0.00 | 0.00 | 93.55 | 23.20 | 5.60 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 2,567 | 5.03 | 8.06 | 0.00 | 0.00 | 95.13 | 9.80 | 9.83 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 2,435 | 5.34 | 9.94 | 0.00 | 0.00 | 94.00 | 14.21 | 9.55 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 2,097 | 4.43 | 11.73 | 0.00 | 0.00 | 91.99 | 23.17 | 11.18 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 2,181 | 2.66 | 4.81 | 0.00 | 0.00 | 96.70 | 13.14 | 0.04 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1,833 | 3.38 | 8.02 | 0.00 | 0.00 | 95.25 | 24.96 | 2.07 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,972 | 4.92 | 6.64 | 0.00 | 0.00 | 96.04 | 25.08 | 1.85 |
| segdup | GRCh38-only graph vs GRCh38 BAM | 2,307 | 1.78 | 8.71 | 0.00 | 0.00 | 97.57 | 7.55 | 5.32 |
| segdup | HPRC v1.1 d9 vs GRCh38 BAM | 2,182 | 2.06 | 10.86 | 0.00 | 0.00 | 96.56 | 12.58 | 5.02 |
| segdup | HPRC v2.1 d46 vs GRCh38 BAM | 1,828 | 1.75 | 12.36 | 0.00 | 0.00 | 94.37 | 25.26 | 5.70 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 26,207 | 0.13 | 0.22 | 0.00 | 0.00 | 99.78 | 0.99 | 0.06 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 25,717 | 0.31 | 1.73 | 0.00 | 0.00 | 99.68 | 1.32 | 1.58 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 25,787 | 0.41 | 1.73 | 0.00 | 0.00 | 99.68 | 1.29 | 1.76 |
| not segdup | GRCh38-only graph vs GRCh38 BAM | 26,110 | 0.52 | 0.38 | 0.00 | 0.00 | 99.72 | 0.64 | 5.68 |
| not segdup | HPRC v1.1 d9 vs GRCh38 BAM | 26,173 | 0.52 | 0.52 | 0.00 | 0.00 | 99.72 | 0.66 | 5.43 |
| not segdup | HPRC v2.1 d46 vs GRCh38 BAM | 25,834 | 0.38 | 1.65 | 0.00 | 0.00 | 99.69 | 1.25 | 6.06 |

## COLO829T_Illumina

### W. Window sample (48 GRCh38 windows, flank 1,000 bp; reads placed in a window by at least one arm: 1,267,256; random 972,840, segdup 294,416)

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38-only graph | window | 92.35 | 7.38 | 0.27 |
| HPRC v1.1 d9 | window | 92.23 | 7.76 | 0.01 |
| HPRC v2.1 d46 | window | 92.06 | 7.67 | 0.27 |

Unmapped records of the whole file (GAM head / BAM unplaced): GRCh38-only graph 24,512,000 (0.99 % of all reads); HPRC v1.1 d9 15,568,000 (0.63 % of all reads); HPRC v2.1 d46 23,242,000 (0.94 % of all reads)

Reads placed in the window by the arm:

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 1,170,285 | 70.72 | 26.87 | 2.84 | 2.63 | 3.98 | 3.69 | 2.96 | 3.75 | 85.61 | 1.0000 | 3.95 | 0.21 | 0.68 | - | - | - | 0.00 |
| HPRC v1.1 d9 | 1,168,846 | 81.25 | 17.08 | 0.58 | 2.64 | 4.83 | 3.88 | 2.79 | 3.38 | 85.11 | 1.0000 | 2.97 | 0.04 | 0.88 | 14.54 | 4.89 | 98.98 | 0.00 |
| HPRC v2.1 d46 | 1,166,614 | 81.61 | 17.00 | 0.56 | 2.36 | 4.74 | 3.60 | 2.75 | 3.45 | 85.45 | 1.0000 | 2.94 | 0.04 | 0.61 | 14.42 | 3.54 | 98.96 | 0.00 |

Same reads, both arms place them in the window (B relative to A); last columns: reads A places in the window that B places elsewhere / leaves unmapped:

| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both  A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,102,189 | 1,102,189 | 1.57 | 3.15 | 0.39 | 1.07 | 71.02 | 81.66 | 10.64 | 0.00 | 2.49 | 2.29 | 99.58 | 100.00 | 16 | 68,027 | 69 | 63,261 | 3,396 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,107,463 | 1,107,463 | 1.57 | 3.20 | 0.38 | 0.91 | 71.01 | 81.73 | 10.72 | 0.00 | 2.49 | 2.30 | 99.46 | 100.00 | 22 | 62,585 | 237 | 58,890 | 261 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,114,035 | 1,114,035 | 0.98 | 0.79 | 0.42 | 0.24 | 81.66 | 81.78 | 0.20 | 0.07 | 2.29 | 2.30 | 99.69 | 100.00 | 9 | 51,473 | 3,338 | 52,552 | 27 |

| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 117,286 | 94.11 | 1.0 | 85.76 | 13.63 | 0.70 | 0.00 | 1.0 | 99.60 | 5.17 | 1.93 | 7 | 57.14 |
| HPRC v2.1 d46 vs GRCh38-only graph | 118,738 | 94.02 | 1.0 | 85.76 | 13.63 | 0.72 | 0.00 | 1.0 | 99.54 | 5.25 | 1.55 | 2 | 0.00 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,206 | 92.20 | 1.0 | 83.50 | 15.28 | 1.36 | 0.00 | 1.0 | 99.14 | 10.47 | 0.73 | 783 | 21.97 |

By window kind and region:

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 972,840 | GRCh38-only graph | 96.96 | 2.76 | 0.28 | 71.73 | 25.82 | 2.73 | 2.54 | 1.77 | 93.53 | 1.0000 |
| windows: random | 972,840 | HPRC v1.1 d9 | 96.97 | 3.02 | 0.01 | 82.34 | 16.01 | 0.38 | 2.57 | 2.03 | 93.55 | 1.0000 |
| windows: random | 972,840 | HPRC v2.1 d46 | 96.35 | 3.37 | 0.28 | 82.78 | 15.86 | 0.35 | 2.29 | 1.91 | 94.03 | 1.0000 |
| windows: segdup | 294,416 | GRCh38-only graph | 77.11 | 22.66 | 0.23 | 66.52 | 31.26 | 3.31 | 3.00 | 13.17 | 52.72 | 1.0000 |
| windows: segdup | 294,416 | HPRC v1.1 d9 | 76.58 | 23.41 | 0.01 | 76.70 | 21.53 | 1.40 | 2.94 | 16.56 | 49.84 | 1.0000 |
| windows: segdup | 294,416 | HPRC v2.1 d46 | 77.89 | 21.88 | 0.24 | 76.83 | 21.67 | 1.46 | 2.62 | 16.34 | 50.40 | 1.0000 |
| SMaHT easy | 795,796 | GRCh38-only graph | 99.67 | 0.07 | 0.26 | 73.67 | 24.31 | 2.06 | 1.82 | 0.06 | 98.10 | 1.0000 |
| SMaHT easy | 795,796 | HPRC v1.1 d9 | 99.87 | 0.13 | 0.00 | 83.82 | 14.81 | 0.26 | 1.94 | 0.09 | 97.88 | 1.0000 |
| SMaHT easy | 795,796 | HPRC v2.1 d46 | 99.66 | 0.08 | 0.26 | 84.17 | 14.73 | 0.23 | 1.69 | 0.09 | 98.01 | 1.0000 |
| SMaHT difficult | 126,883 | GRCh38-only graph | 99.24 | 0.33 | 0.43 | 63.11 | 32.50 | 5.62 | 5.60 | 0.65 | 93.83 | 1.0000 |
| SMaHT difficult | 126,883 | HPRC v1.1 d9 | 99.31 | 0.68 | 0.00 | 75.26 | 21.85 | 0.90 | 5.39 | 0.60 | 94.05 | 1.0000 |
| SMaHT difficult | 126,883 | HPRC v2.1 d46 | 98.98 | 0.59 | 0.42 | 75.84 | 21.69 | 0.83 | 4.94 | 0.51 | 94.32 | 1.0000 |
| SMaHT extreme | 344,577 | GRCh38-only graph | 72.90 | 26.87 | 0.23 | 65.22 | 32.16 | 3.92 | 3.68 | 18.06 | 42.08 | 1.0000 |
| SMaHT extreme | 344,577 | HPRC v1.1 d9 | 72.00 | 27.97 | 0.02 | 76.06 | 21.91 | 1.45 | 3.49 | 22.16 | 39.70 | 1.0000 |
| SMaHT extreme | 344,577 | HPRC v2.1 d46 | 71.96 | 27.80 | 0.24 | 76.34 | 21.89 | 1.49 | 3.18 | 21.78 | 40.79 | 1.0000 |
| segdup | 340,128 | GRCh38-only graph | 73.12 | 26.67 | 0.21 | 67.57 | 30.17 | 3.24 | 2.99 | 16.94 | 45.84 | 1.0000 |
| segdup | 340,128 | HPRC v1.1 d9 | 72.48 | 27.52 | 0.00 | 77.07 | 21.13 | 1.32 | 2.98 | 21.28 | 43.15 | 1.0000 |
| segdup | 340,128 | HPRC v2.1 d46 | 72.51 | 27.27 | 0.22 | 77.32 | 21.14 | 1.38 | 2.68 | 20.88 | 44.45 | 1.0000 |
| not segdup | 927,128 | GRCh38-only graph | 99.40 | 0.30 | 0.29 | 71.58 | 25.99 | 2.74 | 2.53 | 0.49 | 96.35 | 1.0000 |
| not segdup | 927,128 | HPRC v1.1 d9 | 99.48 | 0.51 | 0.01 | 82.37 | 16.00 | 0.38 | 2.55 | 0.44 | 96.33 | 1.0000 |
| not segdup | 927,128 | HPRC v2.1 d46 | 99.23 | 0.48 | 0.29 | 82.76 | 15.89 | 0.35 | 2.27 | 0.42 | 96.44 | 1.0000 |

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | HPRC v1.1 d9 vs GRCh38-only graph | 923,112 | 0.19 | 0.29 | 10.87 | 0.00 | 99.86 | 2.13 | 0.01 |
| windows: random | HPRC v2.1 d46 vs GRCh38-only graph | 921,328 | 0.16 | 0.31 | 11.02 | 0.00 | 99.87 | 2.31 | 0.02 |
| windows: random | HPRC v2.1 d46 vs HPRC v1.1 d9 | 923,075 | 0.13 | 0.15 | 0.22 | 0.05 | 99.89 | 1.87 | 0.28 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 179,077 | 1.45 | 5.08 | 9.47 | 0.00 | 98.17 | 21.11 | 0.01 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 186,135 | 1.49 | 3.88 | 9.26 | 0.00 | 97.39 | 17.99 | 0.02 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 190,960 | 1.83 | 0.66 | 0.11 | 0.16 | 98.73 | 15.01 | 0.30 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 792,129 | 0.06 | 0.18 | 10.38 | 0.00 | 100.00 | 0.13 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 792,726 | 0.05 | 0.13 | 10.50 | 0.00 | 100.00 | 0.04 | 0.01 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 792,275 | 0.13 | 0.09 | 0.17 | 0.05 | 100.00 | 0.05 | 0.26 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 125,124 | 0.20 | 0.28 | 12.31 | 0.00 | 99.98 | 0.63 | 0.00 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 125,325 | 0.19 | 0.28 | 12.56 | 0.00 | 100.00 | 0.44 | 0.03 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 125,270 | 0.17 | 0.15 | 0.34 | 0.09 | 99.98 | 0.17 | 0.42 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 184,936 | 1.98 | 5.43 | 10.63 | 0.00 | 97.54 | 26.35 | 0.03 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 189,412 | 1.92 | 4.58 | 10.43 | 0.00 | 96.83 | 24.55 | 0.04 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 196,490 | 1.76 | 0.91 | 0.22 | 0.16 | 98.28 | 20.50 | 0.30 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 184,746 | 1.59 | 5.16 | 9.09 | 0.00 | 97.59 | 25.71 | 0.01 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 189,475 | 1.52 | 4.08 | 8.87 | 0.00 | 96.89 | 23.79 | 0.02 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 195,843 | 1.87 | 0.78 | 0.16 | 0.16 | 98.32 | 20.27 | 0.29 |
| not segdup | HPRC v1.1 d9 vs GRCh38-only graph | 917,443 | 0.15 | 0.25 | 10.95 | 0.00 | 99.98 | 0.44 | 0.01 |
| not segdup | HPRC v2.1 d46 vs GRCh38-only graph | 917,988 | 0.15 | 0.26 | 11.10 | 0.00 | 99.99 | 0.37 | 0.02 |
| not segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 918,192 | 0.11 | 0.12 | 0.21 | 0.05 | 99.99 | 0.16 | 0.29 |

