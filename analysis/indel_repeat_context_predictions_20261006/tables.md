# Repeat context: HG008T truth INDELs vs Pansoma INDEL predictions (chr1 test)

## 1. Truth INDELs

| class | chr1-22 (GRCh38 frame) | chr1 BED (GRCh38 frame) | chr1-22 (patient frame) | chr1 BED (patient frame) |
|---|---:|---:|---:|---:|
| HP>=7 | 6,360 (74.9%) | 436 (74.1%) | 6,722 (79.1%) | 458 (77.9%) |
| STR2-6>=3copies | 936 (11.0%) | 65 (11.1%) | 795 (9.4%) | 66 (11.2%) |
| HP4-6 | 241 (2.8%) | 7 (1.2%) | 204 (2.4%) | 7 (1.2%) |
| VNTR>6 | 6 (0.1%) | 0 (0.0%) | 15 (0.2%) | 0 (0.0%) |
| complex | 1 (0.0%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) |
| none | 952 (11.2%) | 80 (13.6%) | 415 (4.9%) | 34 (5.8%) |
| in_STR |  |  | 27 (0.3%) | 1 (0.2%) |
| substitution in a repeat |  |  | 253 (3.0%) | 18 (3.1%) |
| substitution |  |  | 24 (0.3%) | 2 (0.3%) |
| no INFO event |  |  | 41 (0.5%) | 2 (0.3%) |
| total | 8,496 | 588 | 8,496 | 588 |

## 2. Pansoma PASS INDEL calls on chr1 (inside the GIAB BED, GRCh38-placed; rtg vcfeval TP / FP), by class

raw = PASS at the run threshold before the PoN; pon = after the INDEL PoN (gnomAD / CoLoRSdb AF >= 0.01). Truth chr1 BED for comparison: HP>=7 436, STR2-6>=3copies 65, HP4-6 7, none 80 (total 588).

### HG008_Illumina_INDEL_base

PASS calls on chr1: 899 placed + 324 off-reference (no GRCh38 allele, not classed).

| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 776 (94.3%) | 90 | 686 | 0.116 | 0 | 0 | 0 |  | 436 | 0.206 | 0.000 |
| STR2-6>=3copies | 0 (0.0%) | 0 | 0 |  | 0 | 0 | 0 |  | 65 | 0.000 | 0.000 |
| HP4-6 | 35 (4.3%) | 0 | 35 | 0.000 | 0 | 0 | 0 |  | 7 | 0.000 | 0.000 |
| none | 12 (1.5%) | 2 | 10 | 0.167 | 0 | 0 | 0 |  | 80 | 0.025 | 0.000 |
| total | 823 | 92 | 731 | 0.112 | 0 | 0 | 0 | 0.000 | 588 | 0.156 | 0.000 |

HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): 7-9 14 / 90 (truth 31), 10-14 51 / 328 (truth 137), 15-19 22 / 188 (truth 128), 20-29 3 / 73 (truth 112), >=30 0 / 7 (truth 28)

### HG008_Illumina_INDEL_base_b1024

PASS calls on chr1: 6,181 placed + 1,376 off-reference (no GRCh38 allele, not classed).

| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 2,970 (54.8%) | 179 | 2791 | 0.060 | 0 | 0 | 0 |  | 436 | 0.411 | 0.000 |
| STR2-6>=3copies | 1,004 (18.5%) | 7 | 997 | 0.007 | 0 | 0 | 0 |  | 65 | 0.108 | 0.000 |
| HP4-6 | 295 (5.4%) | 7 | 288 | 0.024 | 0 | 0 | 0 |  | 7 | 1.000 | 0.000 |
| VNTR>6 | 34 (0.6%) | 0 | 34 | 0.000 | 0 | 0 | 0 |  | 0 |  |  |
| none | 1,114 (20.6%) | 35 | 1079 | 0.031 | 0 | 0 | 0 |  | 80 | 0.438 | 0.000 |
| total | 5,417 | 228 | 5189 | 0.042 | 0 | 0 | 0 | 0.000 | 588 | 0.388 | 0.000 |

HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): 7-9 25 / 333 (truth 31), 10-14 75 / 767 (truth 137), 15-19 56 / 685 (truth 128), 20-29 20 / 841 (truth 112), >=30 3 / 165 (truth 28)

### HG008_Illumina_INDEL_nopartial_b1024

PASS calls on chr1: 4,563 placed + 897 off-reference (no GRCh38 allele, not classed).

| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 2,600 (64.2%) | 180 | 2420 | 0.069 | 0 | 0 | 0 |  | 436 | 0.413 | 0.000 |
| STR2-6>=3copies | 463 (11.4%) | 3 | 460 | 0.006 | 0 | 0 | 0 |  | 65 | 0.046 | 0.000 |
| HP4-6 | 263 (6.5%) | 6 | 257 | 0.023 | 0 | 0 | 0 |  | 7 | 0.857 | 0.000 |
| VNTR>6 | 26 (0.6%) | 0 | 26 | 0.000 | 0 | 0 | 0 |  | 0 |  |  |
| none | 700 (17.3%) | 28 | 672 | 0.040 | 0 | 0 | 0 |  | 80 | 0.350 | 0.000 |
| total | 4,052 | 217 | 3835 | 0.054 | 0 | 0 | 0 | 0.000 | 588 | 0.369 | 0.000 |

HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): 7-9 25 / 306 (truth 31), 10-14 78 / 716 (truth 137), 15-19 56 / 603 (truth 128), 20-29 20 / 674 (truth 112), >=30 1 / 121 (truth 28)

### HG008_Illumina_INDEL_nopartial_b1024_nopc

PASS calls on chr1: 4,627 placed + 598 off-reference (no GRCh38 allele, not classed).

| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 2,412 (61.9%) | 178 | 2234 | 0.074 | 0 | 0 | 0 |  | 436 | 0.408 | 0.000 |
| STR2-6>=3copies | 358 (9.2%) | 3 | 355 | 0.008 | 0 | 0 | 0 |  | 65 | 0.046 | 0.000 |
| HP4-6 | 246 (6.3%) | 6 | 240 | 0.024 | 0 | 0 | 0 |  | 7 | 0.857 | 0.000 |
| VNTR>6 | 32 (0.8%) | 0 | 32 | 0.000 | 0 | 0 | 0 |  | 0 |  |  |
| none | 850 (21.8%) | 27 | 823 | 0.032 | 0 | 0 | 0 |  | 80 | 0.338 | 0.000 |
| total | 3,898 | 214 | 3684 | 0.055 | 0 | 0 | 0 | 0.000 | 588 | 0.364 | 0.000 |

HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): 7-9 25 / 279 (truth 31), 10-14 78 / 676 (truth 137), 15-19 55 / 577 (truth 128), 20-29 18 / 596 (truth 112), >=30 2 / 106 (truth 28)

### HG008_Illumina_INDEL_nopartial_b1024_small

PASS calls on chr1: 3,871 placed + 543 off-reference (no GRCh38 allele, not classed).

| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 2,232 (62.7%) | 169 | 2063 | 0.076 | 405 (55.7%) | 102 | 303 | 0.252 | 436 | 0.388 | 0.234 |
| STR2-6>=3copies | 304 (8.5%) | 2 | 302 | 0.007 | 71 (9.8%) | 0 | 71 | 0.000 | 65 | 0.031 | 0.000 |
| HP4-6 | 243 (6.8%) | 7 | 236 | 0.029 | 61 (8.4%) | 7 | 54 | 0.115 | 7 | 1.000 | 1.000 |
| VNTR>6 | 24 (0.7%) | 0 | 24 | 0.000 | 5 (0.7%) | 0 | 5 | 0.000 | 0 |  |  |
| none | 757 (21.3%) | 33 | 724 | 0.044 | 185 (25.4%) | 33 | 152 | 0.178 | 80 | 0.412 | 0.412 |
| total | 3,560 | 211 | 3349 | 0.059 | 727 | 142 | 585 | 0.195 | 588 | 0.359 | 0.241 |

HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): 7-9 25 / 280 (truth 31), 10-14 73 / 664 (truth 137), 15-19 55 / 551 (truth 128), 20-29 16 / 511 (truth 112), >=30 0 / 57 (truth 28)

### HG008_Illumina_INDEL_nopartial_b1024_nopc_small

PASS calls on chr1: 3,881 placed + 374 off-reference (no GRCh38 allele, not classed).

| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 2,242 (63.1%) | 174 | 2068 | 0.078 | 424 (56.7%) | 105 | 319 | 0.248 | 436 | 0.399 | 0.241 |
| STR2-6>=3copies | 286 (8.0%) | 2 | 284 | 0.007 | 70 (9.4%) | 0 | 70 | 0.000 | 65 | 0.031 | 0.000 |
| HP4-6 | 247 (7.0%) | 6 | 241 | 0.024 | 65 (8.7%) | 6 | 59 | 0.092 | 7 | 0.857 | 0.857 |
| VNTR>6 | 28 (0.8%) | 0 | 28 | 0.000 | 6 (0.8%) | 0 | 6 | 0.000 | 0 |  |  |
| none | 750 (21.1%) | 28 | 722 | 0.037 | 183 (24.5%) | 28 | 155 | 0.153 | 80 | 0.350 | 0.350 |
| total | 3,553 | 210 | 3343 | 0.059 | 748 | 139 | 609 | 0.186 | 588 | 0.357 | 0.236 |

HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): 7-9 25 / 272 (truth 31), 10-14 75 / 659 (truth 137), 15-19 55 / 548 (truth 128), 20-29 18 / 518 (truth 112), >=30 1 / 71 (truth 28)

### HG008_Illumina_INDEL_nopartial_b1024_nopc_small_lr1e4

PASS calls on chr1: 4,682 placed + 565 off-reference (no GRCh38 allele, not classed).

| class | raw calls | raw TP | raw FP | raw precision | pon calls | pon TP | pon FP | pon precision | truth chr1 BED | raw recall | pon recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 2,471 (57.9%) | 174 | 2297 | 0.070 | 485 (52.2%) | 106 | 379 | 0.219 | 436 | 0.399 | 0.243 |
| STR2-6>=3copies | 627 (14.7%) | 3 | 624 | 0.005 | 145 (15.6%) | 0 | 145 | 0.000 | 65 | 0.046 | 0.000 |
| HP4-6 | 259 (6.1%) | 7 | 252 | 0.027 | 67 (7.2%) | 7 | 60 | 0.104 | 7 | 1.000 | 1.000 |
| VNTR>6 | 32 (0.7%) | 0 | 32 | 0.000 | 7 (0.8%) | 0 | 7 | 0.000 | 0 |  |  |
| none | 882 (20.7%) | 32 | 850 | 0.036 | 226 (24.3%) | 32 | 194 | 0.142 | 80 | 0.400 | 0.400 |
| total | 4,271 | 216 | 4055 | 0.051 | 930 | 145 | 785 | 0.156 | 588 | 0.367 | 0.247 |

HP>=7 raw calls by homopolymer length (TP / FP; truth chr1 BED): 7-9 25 / 295 (truth 31), 10-14 75 / 692 (truth 137), 15-19 55 / 592 (truth 128), 20-29 18 / 611 (truth 112), >=30 1 / 107 (truth 28)

