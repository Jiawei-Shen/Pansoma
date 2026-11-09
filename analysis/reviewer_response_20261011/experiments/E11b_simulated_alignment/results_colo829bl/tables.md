# E11b tables (simulated reads from the COLO829BL verkko 2.1 assembly, chr1 + chr20 one-to-one segments)

correct = surjected primary alignment on the true chromosome with start or end within 50 bp of the true GRCh38 interval; wrong = placed on GRCh38 but not correct; perfect = graph alignment (GAF) covers the whole read with no mismatch / indel; MAPQ = giraffe's graph MAPQ. GRCh38 linear reference arm: BAM primary record; perfect = no clip, no I/D, NM 0; MAPQ = the aligner's MAPQ.

## illumina: overall

| arm | reads | unmapped | graph-mapped, not surjected | correct (50 bp) | wrong | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | wrong at MAPQ 60 | MAPQ 0 | uses off-ref node |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GRCh38 linear reference (bwa mem / minimap2) | 5,000,000 | 870 (0.017%) | 0 | 98.725% | 62,904 | 1.258% | 63.83% | 34,367 (1.077% of perfect) | 6,355 (0.130%) | 2,364 (0.049%) | 1.37% | 0 (0.00%) |
| GRCh38-only graph (H 17) | 5,000,000 | 3,260 (0.065%) | 0 | 98.719% | 60,769 | 1.216% | 63.79% | 33,212 (1.041% of perfect) | 5,077 (0.104%) | 3,076 (0.064%) | 0.81% | 0 (0.00%) |
| HPRC v1.1 d9 | 5,000,000 | 1,564 (0.031%) | 16,442 | 98.657% | 49,170 | 0.987% | 72.94% | 33,465 (0.918% of perfect) | 1,641 (0.034%) | 1,101 (0.023%) | 1.03% | 613,657 (12.27%) |
| HPRC v2.1 d46 | 5,000,000 | 1,422 (0.028%) | 18,221 | 98.667% | 47,001 | 0.944% | 73.04% | 32,254 (0.883% of perfect) | 1,279 (0.026%) | 1,032 (0.021%) | 1.02% | 615,587 (12.31%) |
| CHM13-only graph (control; CHM13 truth) | 4,982,247 | 2,318 (0.046%) | 0 | 98.895% | 52,722 | 1.059% | 64.83% | 28,969 (0.897% of perfect) | 5,625 (0.116%) | 3,846 (0.080%) | 0.76% | 0 (0.00%) |

Error-free reads only (no simulated sequencing error; any graph mismatch is then a real difference between the read and the graph path it took):

| arm | reads | perfect | perfect but wrong | wrong / placed | wrong at MAPQ>=30 |
|---|---|---|---|---|---|
| GRCh38 linear reference (bwa mem / minimap2) | 3,702,437 | 86.19% | 34,267 (1.074%) | 1.254% | 4,691 |
| GRCh38-only graph (H 17) | 3,702,437 | 86.14% | 33,075 (1.037%) | 1.204% | 3,751 |
| HPRC v1.1 d9 | 3,702,437 | 98.48% | 33,286 (0.913%) | 0.977% | 1,206 |
| HPRC v2.1 d46 | 3,702,437 | 98.62% | 32,072 (0.878%) | 0.929% | 928 |
| CHM13-only graph (control; CHM13 truth) | 3,689,347 | 87.54% | 28,822 (0.893%) | 1.045% | 4,132 |

## illumina: MAPQ calibration (reads placed on GRCh38; error rate = wrong / placed in the bin; expected = mean 10^(-MAPQ/10))

| MAPQ | GRCh38 linear reference (bwa mem / minimap2): placed / wrong / rate / expected | GRCh38-only graph (H 17): placed / wrong / rate / expected | HPRC v1.1 d9: placed / wrong / rate / expected | HPRC v2.1 d46: placed / wrong / rate / expected | CHM13-only graph (control; CHM13 truth): placed / wrong / rate / expected |
|---|---|---|---|---|---|
| 0 | 68,731 / 45,636 / 6.64e-01 / 1.00e+00 | 40,374 / 27,937 / 6.92e-01 / 1.00e+00 | 43,961 / 26,235 / 5.97e-01 / 1.00e+00 | 42,912 / 25,426 / 5.93e-01 / 1.00e+00 | 37,786 / 25,960 / 6.87e-01 / 1.00e+00 |
| 1-9 | 5,970 / 2,883 / 4.83e-01 / 3.49e-01 | 44,533 / 21,274 / 4.78e-01 / 7.72e-01 | 51,359 / 19,256 / 3.75e-01 / 7.77e-01 | 48,543 / 19,078 / 3.93e-01 / 7.79e-01 | 35,534 / 16,497 / 4.64e-01 / 7.68e-01 |
| 10-19 | 6,315 / 2,385 / 3.78e-01 / 3.81e-02 | 29,224 / 5,277 / 1.81e-01 / 3.93e-02 | 28,553 / 1,650 / 5.78e-02 / 3.99e-02 | 28,120 / 908 / 3.23e-02 / 4.02e-02 | 26,710 / 3,874 / 1.45e-01 / 3.92e-02 |
| 20-29 | 30,461 / 5,645 / 1.85e-01 / 3.38e-03 | 10,369 / 1,204 / 1.16e-01 / 3.08e-03 | 10,455 / 388 / 3.71e-02 / 2.60e-03 | 10,490 / 310 / 2.96e-02 / 2.61e-03 | 9,696 / 766 / 7.90e-02 / 2.92e-03 |
| 30-39 | 6,153 / 876 / 1.42e-01 / 4.36e-04 | 17,223 / 1,029 / 5.97e-02 / 8.15e-04 | 16,067 / 264 / 1.64e-02 / 8.34e-04 | 16,162 / 139 / 8.60e-03 / 8.35e-04 | 16,096 / 883 / 5.49e-02 / 8.21e-04 |
| 40-49 | 67,901 / 2,649 / 3.90e-02 / 6.30e-05 | 16,098 / 688 / 4.27e-02 / 3.80e-05 | 15,896 / 195 / 1.23e-02 / 3.90e-05 | 15,723 / 77 / 4.90e-03 / 3.90e-05 | 15,562 / 651 / 4.18e-02 / 3.80e-05 |
| 50-59 | 13,638 / 466 / 3.42e-02 / 4.00e-06 | 7,293 / 284 / 3.89e-02 / 4.00e-06 | 7,092 / 81 / 1.14e-02 / 4.00e-06 | 7,096 / 31 / 4.37e-03 / 4.00e-06 | 7,409 / 245 / 3.31e-02 / 4.00e-06 |
| 60 | 4,799,961 / 2,364 / 4.93e-04 / 1.00e-06 | 4,831,626 / 3,076 / 6.37e-04 / 1.00e-06 | 4,808,611 / 1,101 / 2.29e-04 / 1.00e-06 | 4,811,311 / 1,032 / 2.14e-04 / 1.00e-06 | 4,831,136 / 3,846 / 7.96e-04 / 1.00e-06 |

## illumina: by COLO829BL variants in read

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | GRCh38 linear reference (bwa mem / minimap2) | 4,286,311 | 0.000% | 99.183% | 0.817% | 74.02% | 24,986 | 833 (0.020%) | 1.21% |
| 0 | GRCh38-only graph (H 17) | 4,286,311 | 0.025% | 99.165% | 0.810% | 73.99% | 24,248 | 478 (0.011%) | 0.63% |
| 0 | HPRC v1.1 d9 | 4,286,311 | 0.027% | 99.003% | 0.844% | 74.03% | 26,170 | 110 (0.003%) | 0.89% |
| 0 | HPRC v2.1 d46 | 4,286,311 | 0.026% | 99.005% | 0.833% | 74.04% | 25,679 | 69 (0.002%) | 0.87% |
| 0 | CHM13-only graph (control; CHM13 truth) | 4,273,381 | 0.040% | 99.109% | 0.851% | 69.34% | 22,493 | 3,712 (0.089%) | 0.62% |
| 1 | GRCh38 linear reference (bwa mem / minimap2) | 593,881 | 0.141% | 96.327% | 3.537% | 2.86% | 7,598 | 3,662 (0.643%) | 2.25% |
| 1 | GRCh38-only graph (H 17) | 593,881 | 0.350% | 96.371% | 3.290% | 2.80% | 7,245 | 3,163 (0.557%) | 1.74% |
| 1 | HPRC v1.1 d9 | 593,881 | 0.060% | 96.643% | 1.613% | 66.87% | 5,561 | 1,106 (0.197%) | 1.81% |
| 1 | HPRC v2.1 d46 | 593,881 | 0.043% | 96.700% | 1.439% | 67.51% | 5,086 | 1,002 (0.178%) | 1.82% |
| 1 | CHM13-only graph (control; CHM13 truth) | 590,470 | 0.086% | 97.772% | 2.143% | 38.71% | 5,557 | 1,295 (0.227%) | 1.51% |
| 2-3 | GRCh38 linear reference (bwa mem / minimap2) | 111,468 | 0.028% | 95.531% | 4.443% | 1.20% | 1,326 | 1,239 (1.168%) | 2.45% |
| 2-3 | GRCh38-only graph (H 17) | 111,468 | 0.083% | 95.642% | 4.279% | 1.17% | 1,293 | 932 (0.884%) | 2.16% |
| 2-3 | HPRC v1.1 d9 | 111,468 | 0.033% | 97.165% | 2.230% | 64.49% | 1,244 | 250 (0.238%) | 1.79% |
| 2-3 | HPRC v2.1 d46 | 111,468 | 0.041% | 97.203% | 1.873% | 65.34% | 1,079 | 119 (0.113%) | 1.77% |
| 2-3 | CHM13-only graph (control; CHM13 truth) | 110,475 | 0.055% | 97.454% | 2.492% | 32.34% | 711 | 390 (0.368%) | 1.62% |
| >=4 | GRCh38 linear reference (bwa mem / minimap2) | 8,340 | 0.000% | 76.535% | 23.465% | 5.50% | 457 | 621 (9.319%) | 9.58% |
| >=4 | GRCh38-only graph (H 17) | 8,340 | 0.228% | 77.842% | 21.980% | 5.13% | 426 | 504 (7.620%) | 8.44% |
| >=4 | HPRC v1.1 d9 | 8,340 | 0.084% | 84.077% | 14.300% | 57.43% | 490 | 175 (2.635%) | 8.13% |
| >=4 | HPRC v2.1 d46 | 8,340 | 0.144% | 84.688% | 11.547% | 59.20% | 410 | 89 (1.356%) | 7.81% |
| >=4 | CHM13-only graph (control; CHM13 truth) | 7,921 | 0.215% | 87.426% | 12.386% | 29.42% | 208 | 228 (3.351%) | 5.83% |

## illumina: by sv or insertion

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| read spans an SV (>=50 bp) of its haplotype | GRCh38 linear reference (bwa mem / minimap2) | 17,020 | 5.088% | 27.867% | 70.639% | 20.47% | 3,215 | 2,052 (37.714%) | 37.28% |
| read spans an SV (>=50 bp) of its haplotype | GRCh38-only graph (H 17) | 17,020 | 11.551% | 29.154% | 67.039% | 19.56% | 3,061 | 2,065 (34.642%) | 30.02% |
| read spans an SV (>=50 bp) of its haplotype | HPRC v1.1 d9 | 17,020 | 1.081% | 24.183% | 44.483% | 61.30% | 1,815 | 259 (7.891%) | 29.54% |
| read spans an SV (>=50 bp) of its haplotype | HPRC v2.1 d46 | 17,020 | 0.640% | 24.512% | 40.705% | 62.73% | 1,618 | 233 (7.194%) | 31.18% |
| read spans an SV (>=50 bp) of its haplotype | CHM13-only graph (control; CHM13 truth) | 15,742 | 1.918% | 57.559% | 41.315% | 47.15% | 3,076 | 751 (8.893%) | 28.25% |
| >=10 read bases absent from GRCh38 | GRCh38 linear reference (bwa mem / minimap2) | 26,769 | 3.250% | 48.213% | 50.168% | 16.36% | 3,714 | 3,028 (21.139%) | 25.12% |
| >=10 read bases absent from GRCh38 | GRCh38-only graph (H 17) | 26,769 | 7.471% | 49.793% | 46.187% | 15.80% | 3,568 | 2,983 (19.819%) | 19.93% |
| >=10 read bases absent from GRCh38 | HPRC v1.1 d9 | 26,769 | 0.728% | 47.454% | 25.360% | 58.75% | 2,453 | 1,066 (9.270%) | 20.08% |
| >=10 read bases absent from GRCh38 | HPRC v2.1 d46 | 26,769 | 0.445% | 47.619% | 23.104% | 60.06% | 2,263 | 1,017 (8.703%) | 21.24% |
| >=10 read bases absent from GRCh38 | CHM13-only graph (control; CHM13 truth) | 25,269 | 1.235% | 72.152% | 26.946% | 40.69% | 3,164 | 895 (5.097%) | 18.06% |
| neither | GRCh38 linear reference (bwa mem / minimap2) | 4,971,672 | 0.000% | 99.000% | 1.000% | 64.10% | 30,594 | 3,282 (0.067%) | 1.25% |
| neither | GRCh38-only graph (H 17) | 4,971,672 | 0.025% | 98.987% | 0.989% | 64.07% | 29,592 | 2,038 (0.042%) | 0.70% |
| neither | HPRC v1.1 d9 | 4,971,672 | 0.028% | 98.935% | 0.901% | 73.02% | 30,936 | 568 (0.012%) | 0.93% |
| neither | HPRC v2.1 d46 | 4,971,672 | 0.026% | 98.945% | 0.868% | 73.11% | 29,916 | 252 (0.005%) | 0.91% |
| neither | CHM13-only graph (control; CHM13 truth) | 4,955,439 | 0.040% | 99.034% | 0.926% | 64.96% | 25,764 | 4,700 (0.097%) | 0.67% |

## illumina: by distance to nearest COLO829BL variant

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| inside read | GRCh38 linear reference (bwa mem / minimap2) | 713,689 | 0.121% | 95.971% | 3.912% | 2.63% | 9,381 | 5,522 (0.809%) | 2.37% |
| inside read | GRCh38-only graph (H 17) | 713,689 | 0.307% | 96.040% | 3.664% | 2.57% | 8,964 | 4,599 (0.677%) | 1.88% |
| inside read | HPRC v1.1 d9 | 713,689 | 0.056% | 96.578% | 1.858% | 66.39% | 7,295 | 1,531 (0.227%) | 1.88% |
| inside read | HPRC v2.1 d46 | 713,689 | 0.044% | 96.639% | 1.623% | 67.07% | 6,575 | 1,210 (0.179%) | 1.89% |
| inside read | CHM13-only graph (control; CHM13 truth) | 708,866 | 0.083% | 97.607% | 2.312% | 37.61% | 6,476 | 1,913 (0.280%) | 1.57% |
| 1-50 bp | GRCh38 linear reference (bwa mem / minimap2) | 363,299 | 0.001% | 98.752% | 1.247% | 73.59% | 2,891 | 324 (0.091%) | 1.44% |
| 1-50 bp | GRCh38-only graph (H 17) | 363,299 | 0.052% | 98.713% | 1.236% | 73.52% | 2,790 | 241 (0.068%) | 0.94% |
| 1-50 bp | HPRC v1.1 d9 | 363,299 | 0.043% | 98.710% | 1.036% | 73.87% | 2,631 | 60 (0.017%) | 1.30% |
| 1-50 bp | HPRC v2.1 d46 | 363,299 | 0.039% | 98.697% | 1.006% | 73.88% | 2,558 | 39 (0.011%) | 1.22% |
| 1-50 bp | CHM13-only graph (control; CHM13 truth) | 361,905 | 0.054% | 98.902% | 1.044% | 69.96% | 2,267 | 310 (0.088%) | 0.85% |
| 51-500 bp | GRCh38 linear reference (bwa mem / minimap2) | 1,689,455 | 0.000% | 99.102% | 0.898% | 74.02% | 10,718 | 508 (0.031%) | 1.20% |
| 51-500 bp | GRCh38-only graph (H 17) | 1,689,455 | 0.028% | 99.085% | 0.888% | 73.98% | 10,336 | 231 (0.014%) | 0.68% |
| 51-500 bp | HPRC v1.1 d9 | 1,689,455 | 0.031% | 98.998% | 0.834% | 74.01% | 10,149 | 50 (0.003%) | 0.95% |
| 51-500 bp | HPRC v2.1 d46 | 1,689,455 | 0.029% | 99.010% | 0.809% | 74.01% | 9,788 | 30 (0.002%) | 0.91% |
| 51-500 bp | CHM13-only graph (control; CHM13 truth) | 1,684,127 | 0.039% | 99.145% | 0.816% | 70.20% | 8,807 | 1,062 (0.064%) | 0.64% |
| >500 bp | GRCh38 linear reference (bwa mem / minimap2) | 2,233,557 | 0.000% | 99.314% | 0.686% | 74.09% | 11,377 | 1 (0.000%) | 1.18% |
| >500 bp | GRCh38-only graph (H 17) | 2,233,557 | 0.019% | 99.300% | 0.681% | 74.07% | 11,122 | 6 (0.000%) | 0.54% |
| >500 bp | HPRC v1.1 d9 | 2,233,557 | 0.022% | 99.053% | 0.820% | 74.08% | 13,390 | 0 (0.000%) | 0.78% |
| >500 bp | HPRC v2.1 d46 | 2,233,557 | 0.021% | 99.051% | 0.822% | 74.08% | 13,333 | 0 (0.000%) | 0.79% |
| >500 bp | CHM13-only graph (control; CHM13 truth) | 2,227,349 | 0.040% | 99.115% | 0.846% | 68.59% | 11,419 | 2,340 (0.107%) | 0.57% |

## illumina: by smaht region

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| easy | GRCh38 linear reference (bwa mem / minimap2) | 4,022,421 | 0.002% | 99.995% | 0.003% | 65.02% | 27 | 47 (0.001%) | 0.00% |
| easy | GRCh38-only graph (H 17) | 4,022,421 | 0.005% | 99.992% | 0.003% | 65.02% | 25 | 41 (0.001%) | 0.00% |
| easy | HPRC v1.1 d9 | 4,022,421 | 0.001% | 99.975% | 0.003% | 73.25% | 87 | 22 (0.001%) | 0.03% |
| easy | HPRC v2.1 d46 | 4,022,421 | 0.001% | 99.976% | 0.002% | 73.33% | 50 | 21 (0.001%) | 0.03% |
| easy | CHM13-only graph (control; CHM13 truth) | 4,019,044 | 0.018% | 99.954% | 0.029% | 65.93% | 381 | 861 (0.021%) | 0.00% |
| difficult | GRCh38 linear reference (bwa mem / minimap2) | 576,452 | 0.099% | 99.352% | 0.549% | 58.58% | 824 | 1,090 (0.190%) | 0.23% |
| difficult | GRCh38-only graph (H 17) | 576,452 | 0.209% | 99.351% | 0.441% | 58.56% | 715 | 1,132 (0.197%) | 0.19% |
| difficult | HPRC v1.1 d9 | 576,452 | 0.009% | 99.315% | 0.211% | 72.14% | 745 | 835 (0.146%) | 0.15% |
| difficult | HPRC v2.1 d46 | 576,452 | 0.012% | 99.314% | 0.211% | 72.26% | 732 | 833 (0.146%) | 0.15% |
| difficult | CHM13-only graph (control; CHM13 truth) | 575,734 | 0.059% | 99.676% | 0.265% | 60.17% | 351 | 492 (0.086%) | 0.14% |
| extreme | GRCh38 linear reference (bwa mem / minimap2) | 401,127 | 0.053% | 85.083% | 14.872% | 59.40% | 33,516 | 5,218 (1.788%) | 16.79% |
| extreme | GRCh38-only graph (H 17) | 401,127 | 0.461% | 85.050% | 14.556% | 59.01% | 32,472 | 3,904 (1.408%) | 9.77% |
| extreme | HPRC v1.1 d9 | 401,127 | 0.369% | 84.485% | 12.369% | 70.99% | 32,633 | 784 (0.299%) | 12.33% |
| extreme | HPRC v2.1 d46 | 401,127 | 0.331% | 84.612% | 11.869% | 71.30% | 31,472 | 425 (0.161%) | 12.11% |
| extreme | CHM13-only graph (control; CHM13 truth) | 387,469 | 0.327% | 86.757% | 12.959% | 60.29% | 28,237 | 4,272 (1.532%) | 9.50% |

## illumina: by str fraction of read

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | GRCh38 linear reference (bwa mem / minimap2) | 4,493,663 | 0.018% | 98.755% | 1.227% | 65.08% | 31,693 | 4,697 (0.107%) | 1.40% |
| 0 | GRCh38-only graph (H 17) | 4,493,663 | 0.062% | 98.742% | 1.197% | 65.05% | 30,600 | 3,778 (0.086%) | 0.82% |
| 0 | HPRC v1.1 d9 | 4,493,663 | 0.034% | 98.679% | 1.002% | 73.18% | 30,730 | 1,404 (0.032%) | 1.04% |
| 0 | HPRC v2.1 d46 | 4,493,663 | 0.030% | 98.686% | 0.964% | 73.27% | 29,674 | 1,108 (0.025%) | 1.02% |
| 0 | CHM13-only graph (control; CHM13 truth) | 4,478,245 | 0.049% | 98.884% | 1.068% | 65.98% | 26,778 | 4,884 (0.112%) | 0.78% |
| (0, 0.5) | GRCh38 linear reference (bwa mem / minimap2) | 500,948 | 0.013% | 98.607% | 1.380% | 53.09% | 2,643 | 1,482 (0.301%) | 1.09% |
| (0, 0.5) | GRCh38-only graph (H 17) | 500,948 | 0.090% | 98.630% | 1.281% | 53.08% | 2,589 | 1,151 (0.235%) | 0.62% |
| (0, 0.5) | HPRC v1.1 d9 | 500,948 | 0.007% | 98.584% | 0.828% | 71.02% | 2,698 | 206 (0.043%) | 0.96% |
| (0, 0.5) | HPRC v2.1 d46 | 500,948 | 0.004% | 98.619% | 0.752% | 71.19% | 2,535 | 148 (0.030%) | 0.94% |
| (0, 0.5) | CHM13-only graph (control; CHM13 truth) | 498,924 | 0.023% | 99.028% | 0.950% | 54.93% | 2,172 | 668 (0.136%) | 0.59% |
| >=0.5 | GRCh38 linear reference (bwa mem / minimap2) | 5,389 | 0.000% | 84.227% | 15.773% | 16.44% | 31 | 176 (3.857%) | 5.18% |
| >=0.5 | GRCh38-only graph (H 17) | 5,389 | 0.557% | 88.068% | 11.439% | 16.16% | 23 | 148 (3.109%) | 4.88% |
| >=0.5 | HPRC v1.1 d9 | 5,389 | 0.353% | 86.862% | 3.225% | 50.40% | 37 | 31 (0.779%) | 4.60% |
| >=0.5 | HPRC v2.1 d46 | 5,389 | 0.742% | 87.753% | 2.050% | 52.70% | 45 | 23 (0.568%) | 4.55% |
| >=0.5 | CHM13-only graph (control; CHM13 truth) | 5,078 | 0.453% | 95.884% | 3.679% | 22.19% | 19 | 73 (1.510%) | 0.79% |

## illumina: by homopolymer ge7

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| yes | GRCh38 linear reference (bwa mem / minimap2) | 837,847 | 0.015% | 99.050% | 0.935% | 58.47% | 3,822 | 1,018 (0.124%) | 1.00% |
| yes | GRCh38-only graph (H 17) | 837,847 | 0.064% | 99.059% | 0.878% | 58.44% | 3,660 | 621 (0.076%) | 0.58% |
| yes | HPRC v1.1 d9 | 837,847 | 0.027% | 99.032% | 0.725% | 72.09% | 4,092 | 241 (0.030%) | 0.71% |
| yes | HPRC v2.1 d46 | 837,847 | 0.025% | 99.051% | 0.690% | 72.23% | 3,861 | 187 (0.023%) | 0.68% |
| yes | CHM13-only graph (control; CHM13 truth) | 835,484 | 0.044% | 99.167% | 0.789% | 59.98% | 3,189 | 905 (0.110%) | 0.53% |
| no | GRCh38 linear reference (bwa mem / minimap2) | 4,162,153 | 0.018% | 98.659% | 1.323% | 64.91% | 30,545 | 5,337 (0.131%) | 1.45% |
| no | GRCh38-only graph (H 17) | 4,162,153 | 0.066% | 98.651% | 1.284% | 64.87% | 29,552 | 4,456 (0.110%) | 0.85% |
| no | HPRC v1.1 d9 | 4,162,153 | 0.032% | 98.581% | 1.040% | 73.11% | 29,373 | 1,400 (0.035%) | 1.10% |
| no | HPRC v2.1 d46 | 4,162,153 | 0.029% | 98.590% | 0.995% | 73.20% | 28,393 | 1,092 (0.027%) | 1.08% |
| no | CHM13-only graph (control; CHM13 truth) | 4,146,763 | 0.047% | 98.841% | 1.113% | 65.80% | 25,780 | 4,720 (0.117%) | 0.80% |

## illumina: by graph alignment uses offref node

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| yes | HPRC v1.1 d9 | 613,657 | 0.000% | 95.579% | 1.790% | 71.65% | 6,907 | 411 (0.072%) | 2.55% |
| yes | HPRC v2.1 d46 | 615,587 | 0.000% | 95.579% | 1.506% | 71.86% | 5,821 | 356 (0.062%) | 2.44% |
| no | HPRC v1.1 d9 | 4,384,779 | 0.000% | 99.122% | 0.878% | 73.15% | 26,558 | 1,230 (0.029%) | 0.82% |
| no | HPRC v2.1 d46 | 4,382,991 | 0.000% | 99.133% | 0.867% | 73.23% | 26,433 | 923 (0.022%) | 0.82% |

## illumina: all arms on the reads that have a one-to-one truth on both GRCh38 and CHM13 (4,982,247 of 5,000,000 reads; GRCh38 arms scored on GRCh38, CHM13 arm on CHM13)

| arm | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|
| GRCh38 linear reference (bwa mem / minimap2) | 0.017% | 98.852% | 1.131% | 63.85% | 30,829 | 5,491 | 1.24% |
| GRCh38-only graph (H 17) | 0.064% | 98.845% | 1.091% | 63.81% | 29,708 | 4,386 | 0.76% |
| HPRC v1.1 d9 | 0.031% | 98.798% | 0.882% | 72.96% | 29,821 | 1,368 | 0.97% |
| HPRC v2.1 d46 | 0.028% | 98.800% | 0.853% | 73.06% | 28,990 | 1,216 | 0.94% |
| CHM13-only graph (control; CHM13 truth) | 0.047% | 98.895% | 1.059% | 64.83% | 28,969 | 5,625 | 0.76% |

## illumina: same reads, two arms (correct = within 50 bp; confident wrong = placed wrong with MAPQ>=30)

| A vs B | correct in both | correct only in A | correct only in B | wrong/unplaced in both | confident wrong only in A | confident wrong only in B |
|---|---|---|---|---|---|---|
| linref vs linear | 4,921,807 | 14,419 | 14,164 | 49,610 | 2,351 | 1,073 |
| linref vs d9 | 4,912,497 | 23,729 | 20,327 | 43,447 | 4,977 | 263 |
| linref vs d46 | 4,912,606 | 23,620 | 20,750 | 43,024 | 5,265 | 189 |
| linear vs d9 | 4,915,995 | 19,976 | 16,829 | 47,200 | 3,684 | 248 |
| linear vs d46 | 4,916,332 | 19,639 | 17,024 | 47,005 | 3,968 | 170 |
| d9 vs d46 | 4,918,455 | 14,369 | 14,901 | 52,275 | 568 | 206 |

## hifi: overall

| arm | reads | unmapped | graph-mapped, not surjected | correct (50 bp) | wrong | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | wrong at MAPQ 60 | MAPQ 0 | uses off-ref node |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GRCh38 linear reference (bwa mem / minimap2) | 40,000 | 0 (0.000%) | 0 | 99.555% | 178 | 0.445% | 0.00% | 0 (n/a of perfect) | 36 (0.091%) | 19 (0.048%) | 0.03% | 0 (0.00%) |
| GRCh38-only graph (H 17) | 40,000 | 12 (0.030%) | 0 | 99.550% | 168 | 0.420% | 0.00% | 0 (n/a of perfect) | 16 (0.041%) | 6 (0.016%) | 0.01% | 0 (0.00%) |
| HPRC v1.1 d9 | 40,000 | 3 (0.007%) | 24 | 99.625% | 123 | 0.308% | 0.00% | 0 (n/a of perfect) | 9 (0.023%) | 5 (0.013%) | 0.07% | 36,421 (91.05%) |
| HPRC v2.1 d46 | 40,000 | 0 (0.000%) | 52 | 99.508% | 145 | 0.363% | 0.00% | 0 (n/a of perfect) | 7 (0.018%) | 3 (0.010%) | 0.00% | 36,390 (90.97%) |
| CHM13-only graph (control; CHM13 truth) | 39,679 | 2 (0.005%) | 0 | 99.579% | 165 | 0.416% | 0.00% | 0 (n/a of perfect) | 31 (0.079%) | 19 (0.051%) | 0.03% | 0 (0.00%) |

## hifi: MAPQ calibration (reads placed on GRCh38; error rate = wrong / placed in the bin; expected = mean 10^(-MAPQ/10))

| MAPQ | GRCh38 linear reference (bwa mem / minimap2): placed / wrong / rate / expected | GRCh38-only graph (H 17): placed / wrong / rate / expected | HPRC v1.1 d9: placed / wrong / rate / expected | HPRC v2.1 d46: placed / wrong / rate / expected | CHM13-only graph (control; CHM13 truth): placed / wrong / rate / expected |
|---|---|---|---|---|---|
| 0 | 10 / 5 / 5.00e-01 / 1.00e+00 | 5 / 4 / 8.00e-01 / 1.00e+00 | 16 / 14 / 8.75e-01 / 1.00e+00 | 0 / 0 / n/a / n/a | 11 / 9 / 8.18e-01 / 1.00e+00 |
| 1-9 | 161 / 97 / 6.02e-01 / 5.71e-01 | 586 / 130 / 2.22e-01 / 4.42e-01 | 691 / 94 / 1.36e-01 / 4.49e-01 | 787 / 134 / 1.70e-01 / 4.03e-01 | 537 / 116 / 2.16e-01 / 4.53e-01 |
| 10-19 | 43 / 20 / 4.65e-01 / 3.74e-02 | 48 / 14 / 2.92e-01 / 4.80e-02 | 29 / 3 / 1.03e-01 / 4.29e-02 | 152 / 4 / 2.63e-02 / 4.69e-02 | 35 / 6 / 1.71e-01 / 4.48e-02 |
| 20-29 | 40 / 20 / 5.00e-01 / 5.41e-03 | 56 / 4 / 7.14e-02 / 3.56e-03 | 40 / 3 / 7.50e-02 / 3.27e-03 | 227 / 0 / 0.00e+00 / 3.67e-03 | 43 / 3 / 6.98e-02 / 3.49e-03 |
| 30-39 | 22 / 5 / 2.27e-01 / 4.17e-04 | 151 / 3 / 1.99e-02 / 3.51e-04 | 126 / 3 / 2.38e-02 / 3.52e-04 | 704 / 2 / 2.84e-03 / 3.52e-04 | 138 / 3 / 2.17e-02 / 3.46e-04 |
| 40-49 | 28 / 8 / 2.86e-01 / 4.70e-05 | 485 / 6 / 1.24e-02 / 3.40e-05 | 413 / 0 / 0.00e+00 / 3.40e-05 | 1,984 / 1 / 5.04e-04 / 3.70e-05 | 445 / 6 / 1.35e-02 / 3.40e-05 |
| 50-59 | 37 / 4 / 1.08e-01 / 4.00e-06 | 1,519 / 1 / 6.58e-04 / 4.00e-06 | 1,395 / 1 / 7.17e-04 / 4.00e-06 | 4,797 / 1 / 2.08e-04 / 4.00e-06 | 1,489 / 3 / 2.01e-03 / 4.00e-06 |
| 60 | 39,659 / 19 / 4.79e-04 / 1.00e-06 | 37,138 / 6 / 1.62e-04 / 1.00e-06 | 37,263 / 5 / 1.34e-04 / 1.00e-06 | 31,297 / 3 / 9.60e-05 / 1.00e-06 | 36,979 / 19 / 5.14e-04 / 1.00e-06 |

## hifi: by COLO829BL variants in read

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | GRCh38 linear reference (bwa mem / minimap2) | 1,817 | 0.000% | 99.890% | 0.110% | 0.00% | 0 | 0 (0.000%) | 0.06% |
| 0 | GRCh38-only graph (H 17) | 1,817 | 0.000% | 99.945% | 0.055% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 0 | HPRC v1.1 d9 | 1,817 | 0.000% | 99.890% | 0.110% | 0.00% | 0 | 0 (0.000%) | 0.06% |
| 0 | HPRC v2.1 d46 | 1,817 | 0.000% | 99.890% | 0.055% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 0 | CHM13-only graph (control; CHM13 truth) | 1,812 | 0.000% | 99.614% | 0.386% | 0.00% | 0 | 3 (0.168%) | 0.00% |
| 1 | GRCh38 linear reference (bwa mem / minimap2) | 2,548 | 0.000% | 98.744% | 1.256% | 0.00% | 0 | 8 (0.317%) | 0.08% |
| 1 | GRCh38-only graph (H 17) | 2,548 | 0.079% | 98.783% | 1.139% | 0.00% | 0 | 6 (0.241%) | 0.12% |
| 1 | HPRC v1.1 d9 | 2,548 | 0.000% | 98.744% | 0.632% | 0.00% | 0 | 2 (0.081%) | 0.67% |
| 1 | HPRC v2.1 d46 | 2,548 | 0.000% | 98.744% | 0.475% | 0.00% | 0 | 2 (0.082%) | 0.00% |
| 1 | CHM13-only graph (control; CHM13 truth) | 2,541 | 0.039% | 98.819% | 1.142% | 0.00% | 0 | 0 (0.000%) | 0.08% |
| 2-3 | GRCh38 linear reference (bwa mem / minimap2) | 4,301 | 0.000% | 99.884% | 0.116% | 0.00% | 0 | 1 (0.023%) | 0.00% |
| 2-3 | GRCh38-only graph (H 17) | 4,301 | 0.023% | 99.861% | 0.116% | 0.00% | 0 | 1 (0.024%) | 0.00% |
| 2-3 | HPRC v1.1 d9 | 4,301 | 0.000% | 99.744% | 0.186% | 0.00% | 0 | 0 (0.000%) | 0.02% |
| 2-3 | HPRC v2.1 d46 | 4,301 | 0.000% | 99.582% | 0.303% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 2-3 | CHM13-only graph (control; CHM13 truth) | 4,269 | 0.000% | 99.344% | 0.656% | 0.00% | 0 | 10 (0.240%) | 0.00% |
| >=4 | GRCh38 linear reference (bwa mem / minimap2) | 31,334 | 0.000% | 99.556% | 0.444% | 0.00% | 0 | 27 (0.087%) | 0.02% |
| >=4 | GRCh38-only graph (H 17) | 31,334 | 0.029% | 99.547% | 0.425% | 0.00% | 0 | 9 (0.029%) | 0.01% |
| >=4 | HPRC v1.1 d9 | 31,334 | 0.010% | 99.665% | 0.310% | 0.00% | 0 | 7 (0.023%) | 0.03% |
| >=4 | HPRC v2.1 d46 | 31,334 | 0.000% | 99.537% | 0.380% | 0.00% | 0 | 5 (0.016%) | 0.00% |
| >=4 | CHM13-only graph (control; CHM13 truth) | 31,057 | 0.003% | 99.672% | 0.325% | 0.00% | 0 | 18 (0.059%) | 0.03% |

## hifi: by sv or insertion

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| read spans an SV (>=50 bp) of its haplotype | GRCh38 linear reference (bwa mem / minimap2) | 2,855 | 0.000% | 97.758% | 2.242% | 0.00% | 0 | 19 (0.680%) | 0.07% |
| read spans an SV (>=50 bp) of its haplotype | GRCh38-only graph (H 17) | 2,855 | 0.350% | 97.443% | 2.214% | 0.00% | 0 | 11 (0.403%) | 0.11% |
| read spans an SV (>=50 bp) of its haplotype | HPRC v1.1 d9 | 2,855 | 0.070% | 98.038% | 1.235% | 0.00% | 0 | 5 (0.182%) | 0.70% |
| read spans an SV (>=50 bp) of its haplotype | HPRC v2.1 d46 | 2,855 | 0.000% | 97.793% | 1.238% | 0.00% | 0 | 4 (0.148%) | 0.00% |
| read spans an SV (>=50 bp) of its haplotype | CHM13-only graph (control; CHM13 truth) | 2,695 | 0.037% | 97.848% | 2.116% | 0.00% | 0 | 8 (0.308%) | 0.11% |
| >=10 read bases absent from GRCh38 | GRCh38 linear reference (bwa mem / minimap2) | 10,000 | 0.000% | 98.880% | 1.120% | 0.00% | 0 | 31 (0.313%) | 0.02% |
| >=10 read bases absent from GRCh38 | GRCh38-only graph (H 17) | 10,000 | 0.090% | 98.820% | 1.091% | 0.00% | 0 | 14 (0.143%) | 0.05% |
| >=10 read bases absent from GRCh38 | HPRC v1.1 d9 | 10,000 | 0.010% | 99.130% | 0.651% | 0.00% | 0 | 8 (0.082%) | 0.23% |
| >=10 read bases absent from GRCh38 | HPRC v2.1 d46 | 10,000 | 0.000% | 98.930% | 0.703% | 0.00% | 0 | 5 (0.052%) | 0.00% |
| >=10 read bases absent from GRCh38 | CHM13-only graph (control; CHM13 truth) | 9,805 | 0.000% | 99.184% | 0.816% | 0.00% | 0 | 12 (0.125%) | 0.09% |
| neither | GRCh38 linear reference (bwa mem / minimap2) | 29,445 | 0.000% | 99.793% | 0.207% | 0.00% | 0 | 3 (0.010%) | 0.03% |
| neither | GRCh38-only graph (H 17) | 29,445 | 0.003% | 99.810% | 0.187% | 0.00% | 0 | 1 (0.003%) | 0.00% |
| neither | HPRC v1.1 d9 | 29,445 | 0.003% | 99.800% | 0.187% | 0.00% | 0 | 1 (0.003%) | 0.02% |
| neither | HPRC v2.1 d46 | 29,445 | 0.000% | 99.715% | 0.234% | 0.00% | 0 | 1 (0.003%) | 0.00% |
| neither | CHM13-only graph (control; CHM13 truth) | 29,323 | 0.003% | 99.720% | 0.276% | 0.00% | 0 | 18 (0.062%) | 0.01% |

## hifi: by distance to nearest COLO829BL variant

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| inside read | GRCh38 linear reference (bwa mem / minimap2) | 38,183 | 0.000% | 99.539% | 0.461% | 0.00% | 0 | 36 (0.095%) | 0.02% |
| inside read | GRCh38-only graph (H 17) | 38,183 | 0.031% | 99.531% | 0.438% | 0.00% | 0 | 16 (0.043%) | 0.01% |
| inside read | HPRC v1.1 d9 | 38,183 | 0.008% | 99.612% | 0.317% | 0.00% | 0 | 9 (0.024%) | 0.07% |
| inside read | HPRC v2.1 d46 | 38,183 | 0.000% | 99.489% | 0.378% | 0.00% | 0 | 7 (0.019%) | 0.00% |
| inside read | CHM13-only graph (control; CHM13 truth) | 37,867 | 0.005% | 99.578% | 0.417% | 0.00% | 0 | 28 (0.075%) | 0.03% |
| 1-50 bp | GRCh38 linear reference (bwa mem / minimap2) | 25 | 0.000% | 96.000% | 4.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 1-50 bp | GRCh38-only graph (H 17) | 25 | 0.000% | 96.000% | 4.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 1-50 bp | HPRC v1.1 d9 | 25 | 0.000% | 96.000% | 4.000% | 0.00% | 0 | 0 (0.000%) | 4.00% |
| 1-50 bp | HPRC v2.1 d46 | 25 | 0.000% | 92.000% | 4.167% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 1-50 bp | CHM13-only graph (control; CHM13 truth) | 25 | 0.000% | 92.000% | 8.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | GRCh38 linear reference (bwa mem / minimap2) | 159 | 0.000% | 100.000% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | GRCh38-only graph (H 17) | 159 | 0.000% | 100.000% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | HPRC v1.1 d9 | 159 | 0.000% | 100.000% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | HPRC v2.1 d46 | 159 | 0.000% | 100.000% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | CHM13-only graph (control; CHM13 truth) | 159 | 0.000% | 100.000% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| >500 bp | GRCh38 linear reference (bwa mem / minimap2) | 1,633 | 0.000% | 99.939% | 0.061% | 0.00% | 0 | 0 (0.000%) | 0.06% |
| >500 bp | GRCh38-only graph (H 17) | 1,633 | 0.000% | 100.000% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| >500 bp | HPRC v1.1 d9 | 1,633 | 0.000% | 99.939% | 0.061% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| >500 bp | HPRC v2.1 d46 | 1,633 | 0.000% | 100.000% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| >500 bp | CHM13-only graph (control; CHM13 truth) | 1,628 | 0.000% | 99.693% | 0.307% | 0.00% | 0 | 3 (0.186%) | 0.00% |

## hifi: by smaht region

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| easy | GRCh38 linear reference (bwa mem / minimap2) | 37,639 | 0.000% | 99.997% | 0.003% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| easy | GRCh38-only graph (H 17) | 37,639 | 0.013% | 99.987% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| easy | HPRC v1.1 d9 | 37,639 | 0.000% | 99.997% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| easy | HPRC v2.1 d46 | 37,639 | 0.000% | 99.995% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| easy | CHM13-only graph (control; CHM13 truth) | 37,424 | 0.000% | 99.984% | 0.016% | 0.00% | 0 | 6 (0.016%) | 0.00% |
| difficult | GRCh38 linear reference (bwa mem / minimap2) | 698 | 0.000% | 99.713% | 0.286% | 0.00% | 0 | 2 (0.286%) | 0.00% |
| difficult | GRCh38-only graph (H 17) | 698 | 0.000% | 99.713% | 0.286% | 0.00% | 0 | 2 (0.288%) | 0.00% |
| difficult | HPRC v1.1 d9 | 698 | 0.000% | 99.713% | 0.286% | 0.00% | 0 | 2 (0.287%) | 0.00% |
| difficult | HPRC v2.1 d46 | 698 | 0.000% | 99.570% | 0.287% | 0.00% | 0 | 2 (0.294%) | 0.00% |
| difficult | CHM13-only graph (control; CHM13 truth) | 694 | 0.000% | 99.712% | 0.288% | 0.00% | 0 | 2 (0.289%) | 0.00% |
| extreme | GRCh38 linear reference (bwa mem / minimap2) | 1,663 | 0.000% | 89.477% | 10.523% | 0.00% | 0 | 34 (2.410%) | 0.60% |
| extreme | GRCh38-only graph (H 17) | 1,663 | 0.421% | 89.597% | 10.024% | 0.00% | 0 | 14 (1.399%) | 0.30% |
| extreme | HPRC v1.1 d9 | 1,663 | 0.180% | 91.161% | 7.392% | 0.00% | 0 | 7 (0.784%) | 1.74% |
| extreme | HPRC v2.1 d46 | 1,663 | 0.000% | 88.455% | 8.860% | 0.00% | 0 | 5 (0.791%) | 0.00% |
| extreme | CHM13-only graph (control; CHM13 truth) | 1,561 | 0.128% | 89.814% | 10.071% | 0.00% | 0 | 23 (2.398%) | 0.70% |

## hifi: by str fraction of read

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| <1% | GRCh38 linear reference (bwa mem / minimap2) | 16,626 | 0.000% | 99.687% | 0.313% | 0.00% | 0 | 9 (0.054%) | 0.02% |
| <1% | GRCh38-only graph (H 17) | 16,626 | 0.042% | 99.663% | 0.295% | 0.00% | 0 | 6 (0.037%) | 0.00% |
| <1% | HPRC v1.1 d9 | 16,626 | 0.018% | 99.735% | 0.229% | 0.00% | 0 | 6 (0.037%) | 0.02% |
| <1% | HPRC v2.1 d46 | 16,626 | 0.000% | 99.651% | 0.271% | 0.00% | 0 | 6 (0.037%) | 0.00% |
| <1% | CHM13-only graph (control; CHM13 truth) | 16,523 | 0.012% | 99.667% | 0.321% | 0.00% | 0 | 14 (0.086%) | 0.01% |
| 1-5% | GRCh38 linear reference (bwa mem / minimap2) | 23,248 | 0.000% | 99.501% | 0.499% | 0.00% | 0 | 19 (0.082%) | 0.03% |
| 1-5% | GRCh38-only graph (H 17) | 23,248 | 0.009% | 99.510% | 0.482% | 0.00% | 0 | 8 (0.035%) | 0.02% |
| 1-5% | HPRC v1.1 d9 | 23,248 | 0.000% | 99.553% | 0.366% | 0.00% | 0 | 3 (0.013%) | 0.11% |
| 1-5% | HPRC v2.1 d46 | 23,248 | 0.000% | 99.428% | 0.414% | 0.00% | 0 | 1 (0.004%) | 0.00% |
| 1-5% | CHM13-only graph (control; CHM13 truth) | 23,053 | 0.000% | 99.519% | 0.481% | 0.00% | 0 | 17 (0.075%) | 0.04% |
| >=5% | GRCh38 linear reference (bwa mem / minimap2) | 126 | 0.000% | 92.064% | 7.937% | 0.00% | 0 | 8 (6.504%) | 0.00% |
| >=5% | GRCh38-only graph (H 17) | 126 | 2.381% | 92.064% | 5.691% | 0.00% | 0 | 2 (1.770%) | 0.00% |
| >=5% | HPRC v1.1 d9 | 126 | 0.000% | 98.413% | 0.000% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| >=5% | HPRC v2.1 d46 | 126 | 0.000% | 95.238% | 3.226% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| >=5% | CHM13-only graph (control; CHM13 truth) | 103 | 0.000% | 99.029% | 0.971% | 0.00% | 0 | 0 (0.000%) | 0.00% |

## hifi: by graph alignment uses offref node

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| yes | HPRC v1.1 d9 | 36,421 | 0.000% | 99.649% | 0.286% | 0.00% | 0 | 4 (0.011%) | 0.08% |
| yes | HPRC v2.1 d46 | 36,390 | 0.000% | 99.522% | 0.336% | 0.00% | 0 | 7 (0.020%) | 0.00% |
| no | HPRC v1.1 d9 | 3,576 | 0.000% | 99.469% | 0.531% | 0.00% | 0 | 5 (0.143%) | 0.00% |
| no | HPRC v2.1 d46 | 3,610 | 0.000% | 99.363% | 0.637% | 0.00% | 0 | 0 (0.000%) | 0.00% |

## hifi: all arms on the reads that have a one-to-one truth on both GRCh38 and CHM13 (39,679 of 40,000 reads; GRCh38 arms scored on GRCh38, CHM13 arm on CHM13)

| arm | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|
| GRCh38 linear reference (bwa mem / minimap2) | 0.000% | 99.604% | 0.396% | 0.00% | 0 | 31 | 0.01% |
| GRCh38-only graph (H 17) | 0.023% | 99.602% | 0.376% | 0.00% | 0 | 11 | 0.01% |
| HPRC v1.1 d9 | 0.008% | 99.655% | 0.282% | 0.00% | 0 | 4 | 0.07% |
| HPRC v2.1 d46 | 0.000% | 99.546% | 0.336% | 0.00% | 0 | 4 | 0.00% |
| CHM13-only graph (control; CHM13 truth) | 0.005% | 99.579% | 0.416% | 0.00% | 0 | 31 | 0.03% |

## hifi: same reads, two arms (correct = within 50 bp; confident wrong = placed wrong with MAPQ>=30)

| A vs B | correct in both | correct only in A | correct only in B | wrong/unplaced in both | confident wrong only in A | confident wrong only in B |
|---|---|---|---|---|---|---|
| linref vs linear | 39,803 | 19 | 17 | 161 | 25 | 5 |
| linref vs d9 | 39,783 | 39 | 67 | 111 | 29 | 2 |
| linref vs d46 | 39,740 | 82 | 63 | 115 | 33 | 4 |
| linear vs d9 | 39,781 | 39 | 69 | 111 | 8 | 1 |
| linear vs d46 | 39,742 | 78 | 61 | 119 | 12 | 3 |
| d9 vs d46 | 39,773 | 77 | 30 | 120 | 4 | 2 |

## ont: overall

| arm | reads | unmapped | graph-mapped, not surjected | correct (50 bp) | wrong | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | wrong at MAPQ 60 | MAPQ 0 | uses off-ref node |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GRCh38 linear reference (bwa mem / minimap2) | 40,000 | 0 (0.000%) | 0 | 99.530% | 188 | 0.470% | 1.00% | 1 (0.249% of perfect) | 39 (0.098%) | 23 (0.058%) | 0.02% | 0 (0.00%) |
| GRCh38-only graph (H 17) | 40,000 | 310 (0.775%) | 0 | 98.535% | 276 | 0.695% | 1.00% | 1 (0.249% of perfect) | 163 (0.416%) | 146 (0.375%) | 0.01% | 0 (0.00%) |
| HPRC v1.1 d9 | 40,000 | 243 (0.607%) | 26 | 98.745% | 233 | 0.586% | 8.40% | 4 (0.119% of perfect) | 135 (0.345%) | 132 (0.339%) | 0.07% | 37,208 (93.02%) |
| HPRC v2.1 d46 | 40,000 | 231 (0.578%) | 59 | 97.720% | 622 | 1.566% | 9.01% | 3 (0.083% of perfect) | 323 (0.825%) | 189 (0.486%) | 0.00% | 37,148 (92.87%) |
| CHM13-only graph (control; CHM13 truth) | 39,684 | 296 (0.746%) | 0 | 98.604% | 258 | 0.655% | 1.52% | 0 (0.000% of perfect) | 145 (0.372%) | 140 (0.362%) | 0.01% | 0 (0.00%) |

## ont: MAPQ calibration (reads placed on GRCh38; error rate = wrong / placed in the bin; expected = mean 10^(-MAPQ/10))

| MAPQ | GRCh38 linear reference (bwa mem / minimap2): placed / wrong / rate / expected | GRCh38-only graph (H 17): placed / wrong / rate / expected | HPRC v1.1 d9: placed / wrong / rate / expected | HPRC v2.1 d46: placed / wrong / rate / expected | CHM13-only graph (control; CHM13 truth): placed / wrong / rate / expected |
|---|---|---|---|---|---|
| 0 | 7 / 5 / 7.14e-01 / 1.00e+00 | 2 / 1 / 5.00e-01 / 1.00e+00 | 15 / 12 / 8.00e-01 / 1.00e+00 | 0 / 0 / n/a / n/a | 2 / 2 / 1.00e+00 / 1.00e+00 |
| 1-9 | 170 / 104 / 6.12e-01 / 5.80e-01 | 234 / 99 / 4.23e-01 / 4.09e-01 | 331 / 78 / 2.36e-01 / 4.24e-01 | 335 / 174 / 5.19e-01 / 3.94e-01 | 217 / 99 / 4.56e-01 / 3.85e-01 |
| 10-19 | 46 / 22 / 4.78e-01 / 4.64e-02 | 107 / 5 / 4.67e-02 / 4.50e-02 | 96 / 6 / 6.25e-02 / 4.32e-02 | 133 / 77 / 5.79e-01 / 4.43e-02 | 98 / 9 / 9.18e-02 / 4.69e-02 |
| 20-29 | 44 / 18 / 4.09e-01 / 4.76e-03 | 142 / 8 / 5.63e-02 / 4.28e-03 | 165 / 2 / 1.21e-02 / 4.27e-03 | 98 / 48 / 4.90e-01 / 4.67e-03 | 138 / 3 / 2.17e-02 / 4.15e-03 |
| 30-39 | 34 / 8 / 2.35e-01 / 4.48e-04 | 132 / 11 / 8.33e-02 / 5.17e-04 | 97 / 2 / 2.06e-02 / 5.00e-04 | 100 / 57 / 5.70e-01 / 4.47e-04 | 118 / 2 / 1.69e-02 / 4.92e-04 |
| 40-49 | 33 / 4 / 1.21e-01 / 4.40e-05 | 65 / 4 / 6.15e-02 / 4.60e-05 | 50 / 1 / 2.00e-02 / 4.50e-05 | 97 / 50 / 5.15e-01 / 4.20e-05 | 54 / 2 / 3.70e-02 / 4.10e-05 |
| 50-59 | 24 / 4 / 1.67e-01 / 4.00e-06 | 51 / 2 / 3.92e-02 / 5.00e-06 | 41 / 0 / 0.00e+00 / 5.00e-06 | 77 / 27 / 3.51e-01 / 5.00e-06 | 45 / 1 / 2.22e-02 / 4.00e-06 |
| 60 | 39,642 / 23 / 5.80e-04 / 1.00e-06 | 38,957 / 146 / 3.75e-03 / 1.00e-06 | 38,936 / 132 / 3.39e-03 / 1.00e-06 | 38,870 / 189 / 4.86e-03 / 1.00e-06 | 38,716 / 140 / 3.62e-03 / 1.00e-06 |

## ont: by COLO829BL variants in read

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | GRCh38 linear reference (bwa mem / minimap2) | 1,958 | 0.000% | 99.898% | 0.102% | 20.43% | 0 | 0 (0.000%) | 0.00% |
| 0 | GRCh38-only graph (H 17) | 1,958 | 0.613% | 99.030% | 0.360% | 20.43% | 0 | 7 (0.363%) | 0.00% |
| 0 | HPRC v1.1 d9 | 1,958 | 0.562% | 99.081% | 0.359% | 20.33% | 0 | 7 (0.365%) | 0.00% |
| 0 | HPRC v2.1 d46 | 1,958 | 0.562% | 97.957% | 1.439% | 20.33% | 0 | 18 (0.935%) | 0.00% |
| 0 | CHM13-only graph (control; CHM13 truth) | 1,952 | 0.768% | 98.566% | 0.671% | 4.87% | 0 | 11 (0.573%) | 0.00% |
| 1 | GRCh38 linear reference (bwa mem / minimap2) | 2,620 | 0.000% | 98.740% | 1.260% | 0.00% | 0 | 7 (0.270%) | 0.11% |
| 1 | GRCh38-only graph (H 17) | 2,620 | 0.382% | 98.092% | 1.533% | 0.00% | 0 | 28 (1.089%) | 0.04% |
| 1 | HPRC v1.1 d9 | 2,620 | 0.267% | 98.092% | 0.925% | 12.52% | 1 | 13 (0.511%) | 0.46% |
| 1 | HPRC v2.1 d46 | 2,620 | 0.267% | 96.794% | 2.160% | 12.94% | 1 | 34 (1.331%) | 0.00% |
| 1 | CHM13-only graph (control; CHM13 truth) | 2,614 | 0.344% | 98.164% | 1.497% | 2.45% | 0 | 16 (0.626%) | 0.04% |
| 2-3 | GRCh38 linear reference (bwa mem / minimap2) | 4,312 | 0.000% | 99.884% | 0.116% | 0.00% | 0 | 0 (0.000%) | 0.00% |
| 2-3 | GRCh38-only graph (H 17) | 4,312 | 0.649% | 98.933% | 0.420% | 0.00% | 0 | 13 (0.307%) | 0.00% |
| 2-3 | HPRC v1.1 d9 | 4,312 | 0.580% | 98.910% | 0.467% | 9.04% | 0 | 13 (0.310%) | 0.07% |
| 2-3 | HPRC v2.1 d46 | 4,312 | 0.557% | 97.635% | 1.681% | 9.65% | 0 | 43 (1.021%) | 0.00% |
| 2-3 | CHM13-only graph (control; CHM13 truth) | 4,280 | 0.654% | 98.505% | 0.847% | 1.59% | 0 | 21 (0.501%) | 0.00% |
| >=4 | GRCh38 linear reference (bwa mem / minimap2) | 31,110 | 0.000% | 99.524% | 0.476% | 0.00% | 1 | 32 (0.103%) | 0.01% |
| >=4 | GRCh38-only graph (H 17) | 31,110 | 0.836% | 98.486% | 0.684% | 0.00% | 1 | 115 (0.377%) | 0.00% |
| >=4 | HPRC v1.1 d9 | 31,110 | 0.643% | 98.756% | 0.589% | 7.22% | 3 | 102 (0.335%) | 0.04% |
| >=4 | HPRC v2.1 d46 | 31,110 | 0.607% | 97.795% | 1.509% | 7.88% | 2 | 228 (0.749%) | 0.00% |
| >=4 | CHM13-only graph (control; CHM13 truth) | 30,838 | 0.791% | 98.657% | 0.556% | 1.23% | 0 | 97 (0.321%) | 0.00% |

## ont: by sv or insertion

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| read spans an SV (>=50 bp) of its haplotype | GRCh38 linear reference (bwa mem / minimap2) | 2,855 | 0.000% | 97.583% | 2.417% | 0.00% | 0 | 21 (0.751%) | 0.14% |
| read spans an SV (>=50 bp) of its haplotype | GRCh38-only graph (H 17) | 2,855 | 1.506% | 95.937% | 2.596% | 0.00% | 0 | 52 (1.877%) | 0.03% |
| read spans an SV (>=50 bp) of its haplotype | HPRC v1.1 d9 | 2,855 | 0.911% | 96.778% | 1.532% | 5.32% | 1 | 24 (0.871%) | 0.56% |
| read spans an SV (>=50 bp) of its haplotype | HPRC v2.1 d46 | 2,855 | 0.771% | 95.867% | 2.424% | 5.74% | 0 | 24 (0.879%) | 0.00% |
| read spans an SV (>=50 bp) of its haplotype | CHM13-only graph (control; CHM13 truth) | 2,699 | 1.075% | 96.591% | 2.360% | 0.74% | 0 | 31 (1.185%) | 0.04% |
| >=10 read bases absent from GRCh38 | GRCh38 linear reference (bwa mem / minimap2) | 9,996 | 0.000% | 98.799% | 1.200% | 0.01% | 1 | 34 (0.344%) | 0.04% |
| >=10 read bases absent from GRCh38 | GRCh38-only graph (H 17) | 9,996 | 0.950% | 97.739% | 1.323% | 0.01% | 1 | 75 (0.768%) | 0.02% |
| >=10 read bases absent from GRCh38 | HPRC v1.1 d9 | 9,996 | 0.650% | 98.229% | 0.878% | 5.51% | 2 | 43 (0.441%) | 0.21% |
| >=10 read bases absent from GRCh38 | HPRC v2.1 d46 | 9,996 | 0.620% | 97.339% | 1.667% | 5.91% | 2 | 67 (0.689%) | 0.00% |
| >=10 read bases absent from GRCh38 | CHM13-only graph (control; CHM13 truth) | 9,801 | 0.847% | 98.204% | 0.957% | 0.66% | 0 | 44 (0.459%) | 0.02% |
| neither | GRCh38 linear reference (bwa mem / minimap2) | 29,451 | 0.000% | 99.783% | 0.217% | 1.36% | 0 | 3 (0.010%) | 0.01% |
| neither | GRCh38-only graph (H 17) | 29,451 | 0.703% | 98.825% | 0.475% | 1.36% | 0 | 82 (0.284%) | 0.00% |
| neither | HPRC v1.1 d9 | 29,451 | 0.587% | 98.930% | 0.482% | 9.43% | 2 | 88 (0.305%) | 0.02% |
| neither | HPRC v2.1 d46 | 29,451 | 0.564% | 97.885% | 1.493% | 10.11% | 1 | 248 (0.859%) | 0.00% |
| neither | CHM13-only graph (control; CHM13 truth) | 29,334 | 0.702% | 98.759% | 0.542% | 1.83% | 0 | 96 (0.333%) | 0.00% |

## ont: by distance to nearest COLO829BL variant

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| inside read | GRCh38 linear reference (bwa mem / minimap2) | 38,042 | 0.000% | 99.511% | 0.489% | 0.00% | 1 | 39 (0.103%) | 0.02% |
| inside read | GRCh38-only graph (H 17) | 38,042 | 0.783% | 98.510% | 0.713% | 0.00% | 1 | 156 (0.419%) | 0.01% |
| inside read | HPRC v1.1 d9 | 38,042 | 0.610% | 98.728% | 0.598% | 7.79% | 4 | 128 (0.344%) | 0.07% |
| inside read | HPRC v2.1 d46 | 38,042 | 0.578% | 97.708% | 1.573% | 8.43% | 3 | 305 (0.819%) | 0.00% |
| inside read | CHM13-only graph (control; CHM13 truth) | 37,732 | 0.745% | 98.606% | 0.654% | 1.35% | 0 | 134 (0.362%) | 0.01% |
| 1-50 bp | GRCh38 linear reference (bwa mem / minimap2) | 36 | 0.000% | 97.222% | 2.778% | 19.44% | 0 | 0 (0.000%) | 0.00% |
| 1-50 bp | GRCh38-only graph (H 17) | 36 | 0.000% | 97.222% | 2.778% | 19.44% | 0 | 1 (2.857%) | 0.00% |
| 1-50 bp | HPRC v1.1 d9 | 36 | 0.000% | 97.222% | 2.778% | 19.44% | 0 | 1 (2.857%) | 0.00% |
| 1-50 bp | HPRC v2.1 d46 | 36 | 0.000% | 91.667% | 5.714% | 16.67% | 0 | 0 (0.000%) | 0.00% |
| 1-50 bp | CHM13-only graph (control; CHM13 truth) | 36 | 2.778% | 91.667% | 5.714% | 5.56% | 0 | 2 (5.882%) | 0.00% |
| 51-500 bp | GRCh38 linear reference (bwa mem / minimap2) | 205 | 0.000% | 100.000% | 0.000% | 24.88% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | GRCh38-only graph (H 17) | 205 | 0.976% | 99.024% | 0.000% | 24.88% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | HPRC v1.1 d9 | 205 | 0.976% | 99.024% | 0.000% | 24.88% | 0 | 0 (0.000%) | 0.00% |
| 51-500 bp | HPRC v2.1 d46 | 205 | 0.976% | 98.049% | 0.985% | 24.88% | 0 | 1 (0.498%) | 0.00% |
| 51-500 bp | CHM13-only graph (control; CHM13 truth) | 204 | 0.980% | 99.020% | 0.000% | 7.84% | 0 | 0 (0.000%) | 0.00% |
| >500 bp | GRCh38 linear reference (bwa mem / minimap2) | 1,717 | 0.000% | 99.942% | 0.058% | 19.92% | 0 | 0 (0.000%) | 0.00% |
| >500 bp | GRCh38-only graph (H 17) | 1,717 | 0.582% | 99.068% | 0.351% | 19.92% | 0 | 6 (0.355%) | 0.00% |
| >500 bp | HPRC v1.1 d9 | 1,717 | 0.524% | 99.126% | 0.351% | 19.80% | 0 | 6 (0.356%) | 0.00% |
| >500 bp | HPRC v2.1 d46 | 1,717 | 0.524% | 98.078% | 1.405% | 19.86% | 0 | 17 (1.005%) | 0.00% |
| >500 bp | CHM13-only graph (control; CHM13 truth) | 1,712 | 0.701% | 98.657% | 0.647% | 4.50% | 0 | 9 (0.534%) | 0.00% |

## ont: by smaht region

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| easy | GRCh38 linear reference (bwa mem / minimap2) | 37,591 | 0.000% | 99.992% | 0.008% | 1.01% | 0 | 2 (0.005%) | 0.00% |
| easy | GRCh38-only graph (H 17) | 37,591 | 0.747% | 99.005% | 0.249% | 1.01% | 0 | 93 (0.249%) | 0.00% |
| easy | HPRC v1.1 d9 | 37,591 | 0.583% | 99.146% | 0.270% | 8.53% | 0 | 101 (0.270%) | 0.00% |
| easy | HPRC v2.1 d46 | 37,591 | 0.564% | 99.157% | 0.276% | 9.16% | 0 | 103 (0.276%) | 0.00% |
| easy | CHM13-only graph (control; CHM13 truth) | 37,384 | 0.725% | 99.026% | 0.251% | 1.53% | 0 | 93 (0.251%) | 0.00% |
| difficult | GRCh38 linear reference (bwa mem / minimap2) | 712 | 0.000% | 99.719% | 0.281% | 0.70% | 0 | 2 (0.281%) | 0.00% |
| difficult | GRCh38-only graph (H 17) | 712 | 0.843% | 98.736% | 0.425% | 0.70% | 0 | 3 (0.425%) | 0.00% |
| difficult | HPRC v1.1 d9 | 712 | 0.421% | 99.017% | 0.564% | 7.58% | 1 | 4 (0.565%) | 0.00% |
| difficult | HPRC v2.1 d46 | 712 | 0.421% | 98.315% | 1.130% | 6.88% | 1 | 7 (0.999%) | 0.00% |
| difficult | CHM13-only graph (control; CHM13 truth) | 706 | 0.850% | 98.725% | 0.429% | 1.13% | 0 | 3 (0.429%) | 0.00% |
| extreme | GRCh38 linear reference (bwa mem / minimap2) | 1,697 | 0.000% | 89.216% | 10.784% | 0.94% | 1 | 35 (2.444%) | 0.35% |
| extreme | GRCh38-only graph (H 17) | 1,697 | 1.355% | 88.038% | 10.753% | 0.94% | 1 | 67 (5.616%) | 0.12% |
| extreme | HPRC v1.1 d9 | 1,697 | 1.238% | 89.747% | 7.753% | 6.01% | 3 | 30 (2.841%) | 1.59% |
| extreme | HPRC v2.1 d46 | 1,697 | 0.943% | 65.645% | 31.446% | 6.60% | 2 | 213 (19.925%) | 0.00% |
| extreme | CHM13-only graph (control; CHM13 truth) | 1,594 | 1.192% | 88.645% | 10.286% | 1.51% | 0 | 49 (4.371%) | 0.13% |

## ont: by str fraction of read

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| <1% | GRCh38 linear reference (bwa mem / minimap2) | 16,715 | 0.000% | 99.701% | 0.299% | 1.39% | 0 | 11 (0.066%) | 0.02% |
| <1% | GRCh38-only graph (H 17) | 16,715 | 0.778% | 98.714% | 0.513% | 1.39% | 0 | 49 (0.298%) | 0.00% |
| <1% | HPRC v1.1 d9 | 16,715 | 0.610% | 98.941% | 0.445% | 9.85% | 2 | 50 (0.304%) | 0.04% |
| <1% | HPRC v2.1 d46 | 16,715 | 0.586% | 97.918% | 1.392% | 10.42% | 1 | 133 (0.810%) | 0.00% |
| <1% | CHM13-only graph (control; CHM13 truth) | 16,617 | 0.746% | 98.712% | 0.546% | 2.00% | 0 | 52 (0.318%) | 0.00% |
| 1-5% | GRCh38 linear reference (bwa mem / minimap2) | 23,156 | 0.000% | 99.452% | 0.548% | 0.73% | 1 | 20 (0.087%) | 0.02% |
| 1-5% | GRCh38-only graph (H 17) | 23,156 | 0.756% | 98.454% | 0.796% | 0.73% | 1 | 108 (0.477%) | 0.01% |
| 1-5% | HPRC v1.1 d9 | 23,156 | 0.600% | 98.618% | 0.687% | 7.40% | 2 | 85 (0.376%) | 0.09% |
| 1-5% | HPRC v2.1 d46 | 23,156 | 0.570% | 97.620% | 1.657% | 8.04% | 2 | 186 (0.822%) | 0.00% |
| 1-5% | CHM13-only graph (control; CHM13 truth) | 22,960 | 0.745% | 98.532% | 0.728% | 1.19% | 0 | 91 (0.405%) | 0.01% |
| >=5% | GRCh38 linear reference (bwa mem / minimap2) | 129 | 0.000% | 91.473% | 8.527% | 0.00% | 0 | 8 (6.452%) | 0.00% |
| >=5% | GRCh38-only graph (H 17) | 129 | 3.876% | 89.922% | 6.452% | 0.00% | 0 | 6 (4.959%) | 0.00% |
| >=5% | HPRC v1.1 d9 | 129 | 1.550% | 96.124% | 0.800% | 0.78% | 0 | 0 (0.000%) | 0.00% |
| >=5% | HPRC v2.1 d46 | 129 | 0.775% | 89.922% | 7.937% | 1.55% | 0 | 4 (3.448%) | 0.00% |
| >=5% | CHM13-only graph (control; CHM13 truth) | 107 | 0.935% | 97.196% | 1.887% | 0.00% | 0 | 2 (1.942%) | 0.00% |

## ont: by graph alignment uses offref node

| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|
| yes | HPRC v1.1 d9 | 37,208 | 0.000% | 99.366% | 0.565% | 7.77% | 4 | 120 (0.328%) | 0.07% |
| yes | HPRC v2.1 d46 | 37,148 | 0.000% | 98.288% | 1.556% | 8.42% | 2 | 307 (0.839%) | 0.00% |
| no | HPRC v1.1 d9 | 2,549 | 0.000% | 99.098% | 0.902% | 18.52% | 0 | 15 (0.601%) | 0.00% |
| no | HPRC v2.1 d46 | 2,621 | 0.000% | 98.283% | 1.717% | 18.28% | 1 | 16 (0.622%) | 0.00% |

## ont: all arms on the reads that have a one-to-one truth on both GRCh38 and CHM13 (39,684 of 40,000 reads; GRCh38 arms scored on GRCh38, CHM13 arm on CHM13)

| arm | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | MAPQ 0 |
|---|---|---|---|---|---|---|---|
| GRCh38 linear reference (bwa mem / minimap2) | 0.000% | 99.582% | 0.418% | 1.01% | 1 | 36 | 0.01% |
| GRCh38-only graph (H 17) | 0.761% | 98.601% | 0.642% | 1.01% | 1 | 153 | 0.01% |
| HPRC v1.1 d9 | 0.602% | 98.780% | 0.561% | 8.44% | 4 | 129 | 0.07% |
| HPRC v2.1 d46 | 0.575% | 97.805% | 1.497% | 9.05% | 3 | 315 | 0.00% |
| CHM13-only graph (control; CHM13 truth) | 0.746% | 98.604% | 0.655% | 1.52% | 0 | 145 | 0.01% |

## ont: same reads, two arms (correct = within 50 bp; confident wrong = placed wrong with MAPQ>=30)

| A vs B | correct in both | correct only in A | correct only in B | wrong/unplaced in both | confident wrong only in A | confident wrong only in B |
|---|---|---|---|---|---|---|
| linref vs linear | 39,394 | 418 | 20 | 168 | 17 | 141 |
| linref vs d9 | 39,432 | 380 | 66 | 122 | 30 | 126 |
| linref vs d46 | 39,037 | 775 | 51 | 137 | 35 | 319 |
| linear vs d9 | 39,373 | 41 | 125 | 461 | 51 | 23 |
| linear vs d46 | 38,968 | 446 | 120 | 466 | 68 | 228 |
| d9 vs d46 | 39,044 | 454 | 44 | 458 | 25 | 213 |

