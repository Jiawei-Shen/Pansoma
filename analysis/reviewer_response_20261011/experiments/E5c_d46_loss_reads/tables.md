# E5c tables: truth alleles lost on HPRC v2.1 d46 against GRCh38-only, by how their ALT reads align on d46

Lost = a tensor allele equals the truth on GRCh38-only (status tensor_representative or tensor_non_representative_allele) and none does on d46 (no_candidate or filtered). chr1-22, all truth records. Class = majority representation of the locus's error-free ALT reads (reads.py / report.py).


## INDEL: class of each lost locus (% of lost loci)

| set | d46 status | n lost | I1 | I2 | I3 | I10 | I4 | I5 | I6 | I7 | I8 | I9 | I1 + I2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HG008T_Illumina_d46 | all | 3,528 | 47.5 | 48.0 | 2.2 | 0.4 | 1.2 | 0.3 | 0.0 | 0.4 | 0.0 | 0.0 | 95.5 |
| HG008T_Illumina_d46 | no_candidate | 2,831 | 52.5 | 43.7 | 1.7 | 0.0 | 1.3 | 0.4 | 0.0 | 0.4 | 0.0 | 0.0 | 96.1 |
| HG008T_Illumina_d46 | filtered | 697 | 27.4 | 65.4 | 4.0 | 2.2 | 0.7 | 0.0 | 0.0 | 0.3 | 0.0 | 0.0 | 92.8 |
| HG008T_PacBio_d46 | all | 3,577 | 43.2 | 52.4 | 2.8 | 0.3 | 0.6 | 0.1 | 0.1 | 0.4 | 0.0 | 0.1 | 95.6 |
| HG008T_PacBio_d46 | no_candidate | 2,899 | 41.8 | 53.6 | 3.0 | 0.0 | 0.8 | 0.2 | 0.1 | 0.5 | 0.0 | 0.1 | 95.4 |
| HG008T_PacBio_d46 | filtered | 678 | 49.3 | 47.2 | 2.1 | 1.3 | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 96.5 |
| HG008T_ONT_d46 | all | 2,494 | 44.5 | 51.4 | 2.5 | 0.6 | 0.7 | 0.0 | 0.0 | 0.3 | 0.0 | 0.0 | 95.8 |
| HG008T_ONT_d46 | no_candidate | 1,846 | 43.7 | 52.2 | 2.8 | 0.0 | 0.9 | 0.0 | 0.1 | 0.4 | 0.0 | 0.1 | 95.9 |
| HG008T_ONT_d46 | filtered | 648 | 46.8 | 48.9 | 1.9 | 2.2 | 0.3 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 95.7 |
| COLO829T_Illumina_d46 | all | 612 | 59.2 | 30.2 | 6.5 | 0.0 | 2.8 | 0.5 | 0.0 | 0.7 | 0.0 | 0.2 | 89.4 |
| COLO829T_Illumina_d46 | no_candidate | 523 | 62.7 | 26.2 | 6.5 | 0.0 | 3.1 | 0.6 | 0.0 | 0.8 | 0.0 | 0.2 | 88.9 |
| COLO829T_Illumina_d46 | filtered | 89 | 38.2 | 53.9 | 6.7 | 0.0 | 1.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 92.1 |
| COLO829T_fiberseq_d46 | all | 567 | 53.3 | 32.5 | 9.9 | 0.0 | 1.9 | 1.4 | 0.0 | 0.9 | 0.0 | 0.2 | 85.7 |
| COLO829T_fiberseq_d46 | no_candidate | 424 | 53.3 | 33.0 | 8.0 | 0.0 | 2.4 | 1.9 | 0.0 | 1.2 | 0.0 | 0.2 | 86.3 |
| COLO829T_fiberseq_d46 | filtered | 143 | 53.1 | 30.8 | 15.4 | 0.0 | 0.7 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 83.9 |
| COLO829T_ONT_d46 | all | 457 | 52.1 | 33.3 | 9.6 | 1.8 | 1.3 | 0.4 | 0.0 | 1.5 | 0.0 | 0.0 | 85.3 |
| COLO829T_ONT_d46 | no_candidate | 322 | 52.8 | 33.5 | 9.3 | 0.0 | 1.9 | 0.6 | 0.0 | 1.9 | 0.0 | 0.0 | 86.3 |
| COLO829T_ONT_d46 | filtered | 135 | 50.4 | 32.6 | 10.4 | 5.9 | 0.0 | 0.0 | 0.0 | 0.7 | 0.0 | 0.0 | 83.0 |

## INDEL: lost loci whose ALT reads go through graph alt nodes (entire or partial)

| set | n lost with >= 3 ALT reads | any ALT read on alt nodes | majority of ALT reads on alt nodes | ALT reads pooled: entire / partial / GRCh38 other edit / truth edit / other (%) |
|---|---|---|---|---|
| HG008T_Illumina_d46 | 3,474 | 98.2 | 97.1 | 44.1 / 51.5 / 2.7 / 0.8 / 0.9 |
| HG008T_PacBio_d46 | 3,550 | 97.6 | 96.5 | 45.1 / 49.9 / 2.9 / 1.3 / 0.8 |
| HG008T_ONT_d46 | 2,476 | 97.8 | 96.6 | 36.0 / 57.4 / 3.6 / 2.1 / 0.9 |
| COLO829T_Illumina_d46 | 592 | 93.1 | 92.4 | 53.4 / 37.0 / 7.5 / 0.6 / 1.5 |
| COLO829T_fiberseq_d46 | 548 | 90.1 | 88.9 | 49.4 / 37.2 / 9.5 / 1.8 / 2.1 |
| COLO829T_ONT_d46 | 449 | 89.8 | 86.9 | 40.8 / 43.2 / 11.0 / 2.1 / 2.9 |

## SNV: class of each lost locus (% of lost loci)

| set | d46 status | n lost | I1 | I2 | I3 | I10 | I4 | I5 | I6 | I7 | I8 | I9 | I1 + I2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HG008T_Illumina_d46 | all | 263 | 44.5 | 20.2 | 8.7 | 0.4 | 19.0 | 2.3 | 0.0 | 4.9 | 0.0 | 0.0 | 64.6 |
| HG008T_Illumina_d46 | no_candidate | 207 | 48.8 | 19.3 | 4.3 | 0.0 | 19.8 | 2.9 | 0.0 | 4.8 | 0.0 | 0.0 | 68.1 |
| HG008T_Illumina_d46 | filtered | 56 | 28.6 | 23.2 | 25.0 | 1.8 | 16.1 | 0.0 | 0.0 | 5.4 | 0.0 | 0.0 | 51.8 |
| HG008T_PacBio_d46 | all | 316 | 42.1 | 20.9 | 10.1 | 0.6 | 14.6 | 3.5 | 0.0 | 8.2 | 0.0 | 0.0 | 63.0 |
| HG008T_PacBio_d46 | no_candidate | 277 | 45.1 | 21.3 | 7.9 | 0.0 | 13.4 | 4.0 | 0.0 | 8.3 | 0.0 | 0.0 | 66.4 |
| HG008T_PacBio_d46 | filtered | 39 | 20.5 | 17.9 | 25.6 | 5.1 | 23.1 | 0.0 | 0.0 | 7.7 | 0.0 | 0.0 | 38.5 |
| HG008T_ONT_d46 | all | 290 | 43.8 | 22.8 | 11.0 | 0.3 | 14.1 | 0.3 | 0.0 | 7.6 | 0.0 | 0.0 | 66.6 |
| HG008T_ONT_d46 | no_candidate | 242 | 47.9 | 22.3 | 8.3 | 0.0 | 13.2 | 0.4 | 0.0 | 7.9 | 0.0 | 0.0 | 70.2 |
| HG008T_ONT_d46 | filtered | 48 | 22.9 | 25.0 | 25.0 | 2.1 | 18.8 | 0.0 | 0.0 | 6.2 | 0.0 | 0.0 | 47.9 |
| COLO829T_Illumina_d46 | all | 398 | 11.6 | 3.8 | 2.0 | 4.3 | 32.2 | 45.5 | 0.0 | 0.8 | 0.0 | 0.0 | 15.3 |
| COLO829T_Illumina_d46 | no_candidate | 342 | 12.9 | 2.9 | 1.5 | 0.0 | 30.4 | 51.5 | 0.0 | 0.9 | 0.0 | 0.0 | 15.8 |
| COLO829T_Illumina_d46 | filtered | 56 | 3.6 | 8.9 | 5.4 | 30.4 | 42.9 | 8.9 | 0.0 | 0.0 | 0.0 | 0.0 | 12.5 |
| COLO829T_fiberseq_d46 | all | 281 | 18.1 | 2.1 | 2.8 | 10.3 | 33.1 | 32.0 | 0.0 | 1.1 | 0.0 | 0.4 | 20.3 |
| COLO829T_fiberseq_d46 | no_candidate | 222 | 22.1 | 1.8 | 2.7 | 0.0 | 31.1 | 40.5 | 0.0 | 1.4 | 0.0 | 0.5 | 23.9 |
| COLO829T_fiberseq_d46 | filtered | 59 | 3.4 | 3.4 | 3.4 | 49.2 | 40.7 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 6.8 |
| COLO829T_ONT_d46 | all | 220 | 20.5 | 4.1 | 2.3 | 34.5 | 33.2 | 3.6 | 0.0 | 1.4 | 0.0 | 0.5 | 24.5 |
| COLO829T_ONT_d46 | no_candidate | 102 | 43.1 | 6.9 | 1.0 | 0.0 | 38.2 | 7.8 | 0.0 | 2.0 | 0.0 | 1.0 | 50.0 |
| COLO829T_ONT_d46 | filtered | 118 | 0.8 | 1.7 | 3.4 | 64.4 | 28.8 | 0.0 | 0.0 | 0.8 | 0.0 | 0.0 | 2.5 |

## SNV: lost loci whose ALT reads go through graph alt nodes (entire or partial)

| set | n lost with >= 3 ALT reads | any ALT read on alt nodes | majority of ALT reads on alt nodes | ALT reads pooled: entire / partial / GRCh38 other edit / truth edit / other (%) |
|---|---|---|---|---|
| HG008T_Illumina_d46 | 207 | 85.0 | 81.6 | 55.6 / 29.0 / 8.1 / 0.6 / 6.7 |
| HG008T_PacBio_d46 | 259 | 81.5 | 77.6 | 49.7 / 27.2 / 13.4 / 0.6 / 9.2 |
| HG008T_ONT_d46 | 248 | 81.5 | 78.6 | 47.3 / 33.0 / 11.2 / 1.3 / 7.1 |
| COLO829T_Illumina_d46 | 89 | 69.7 | 68.5 | 61.2 / 20.6 / 7.6 / 5.7 / 5.0 |
| COLO829T_fiberseq_d46 | 98 | 59.2 | 58.2 | 62.9 / 10.0 / 9.8 / 11.9 / 5.4 |
| COLO829T_ONT_d46 | 139 | 42.4 | 39.6 | 56.0 / 11.7 / 8.9 / 18.0 / 5.3 |

## INDEL: class by repeat context (lost loci; % of the class row)

| set | repeat | n lost | I1 | I2 | I3 | I10 | I4 + I5 | rest |
|---|---|---|---|---|---|---|---|---|
| HG008T_Illumina_d46 | HP>=7 | 2,739 | 49.2 | 49.3 | 0.3 | 0.4 | 0.8 | 0.0 |
| HG008T_Illumina_d46 | STR | 379 | 43.3 | 51.2 | 1.1 | 0.5 | 4.0 | 0.0 |
| HG008T_Illumina_d46 | other | 410 | 40.2 | 36.3 | 15.4 | 0.2 | 4.4 | 3.4 |
| HG008T_PacBio_d46 | HP>=7 | 2,448 | 44.6 | 54.9 | 0.0 | 0.2 | 0.2 | 0.0 |
| HG008T_PacBio_d46 | STR | 463 | 43.8 | 54.2 | 1.1 | 0.0 | 0.6 | 0.2 |
| HG008T_PacBio_d46 | other | 666 | 37.7 | 41.9 | 14.1 | 0.6 | 3.0 | 2.7 |
| HG008T_ONT_d46 | HP>=7 | 1,630 | 45.4 | 53.8 | 0.1 | 0.7 | 0.1 | 0.0 |
| HG008T_ONT_d46 | STR | 423 | 46.8 | 52.2 | 0.5 | 0.0 | 0.2 | 0.2 |
| HG008T_ONT_d46 | other | 441 | 38.8 | 41.5 | 13.6 | 0.7 | 3.6 | 1.8 |
| COLO829T_Illumina_d46 | HP>=7 | 345 | 70.4 | 27.2 | 0.9 | 0.0 | 1.2 | 0.3 |
| COLO829T_Illumina_d46 | STR | 117 | 58.1 | 35.9 | 2.6 | 0.0 | 3.4 | 0.0 |
| COLO829T_Illumina_d46 | other | 150 | 34.0 | 32.7 | 22.7 | 0.0 | 8.0 | 2.7 |
| COLO829T_fiberseq_d46 | HP>=7 | 234 | 66.7 | 32.9 | 0.0 | 0.0 | 0.4 | 0.0 |
| COLO829T_fiberseq_d46 | STR | 115 | 65.2 | 31.3 | 0.9 | 0.0 | 2.6 | 0.0 |
| COLO829T_fiberseq_d46 | other | 218 | 32.6 | 32.6 | 25.2 | 0.0 | 6.9 | 2.8 |
| COLO829T_ONT_d46 | HP>=7 | 169 | 60.4 | 35.5 | 0.0 | 4.1 | 0.0 | 0.0 |
| COLO829T_ONT_d46 | STR | 109 | 67.0 | 30.3 | 0.9 | 0.0 | 1.8 | 0.0 |
| COLO829T_ONT_d46 | other | 179 | 35.2 | 33.0 | 24.0 | 0.6 | 3.4 | 3.9 |

## INDEL: lost loci not PoN-tagged (what the graph costs after population filtering), by class

| set | n lost, not tagged | I1 | I2 | I3 | I10 | I4 + I5 | rest |
|---|---|---|---|---|---|---|---|
| HG008T_Illumina_d46 | 515 | 7.8 | 80.4 | 7.8 | 0.8 | 2.7 | 0.6 |
| HG008T_PacBio_d46 | 603 | 9.1 | 80.1 | 8.0 | 0.5 | 1.5 | 0.8 |
| HG008T_ONT_d46 | 355 | 8.2 | 81.7 | 7.9 | 0.8 | 1.4 | 0.0 |
| COLO829T_Illumina_d46 | 118 | 9.3 | 56.8 | 28.0 | 0.0 | 4.2 | 1.7 |
| COLO829T_fiberseq_d46 | 133 | 8.3 | 49.6 | 34.6 | 0.0 | 6.8 | 0.8 |
| COLO829T_ONT_d46 | 110 | 6.4 | 50.0 | 33.6 | 3.6 | 5.5 | 0.9 |

## Controls (tensor on both sets): class of the locus (%), same method

| set | kind | n | I10 (truth edit made) | I1 | I2 | I3 | other |
|---|---|---|---|---|---|---|---|
| HG008T_Illumina_d46 | INDEL | 100 | 91.0 | 2.0 | 7.0 | 0.0 | 0.0 |
| HG008T_Illumina_d46 | SNV | 50 | 98.0 | 0.0 | 0.0 | 2.0 | 0.0 |
| HG008T_PacBio_d46 | INDEL | 100 | 84.0 | 10.0 | 6.0 | 0.0 | 0.0 |
| HG008T_PacBio_d46 | SNV | 50 | 98.0 | 0.0 | 0.0 | 0.0 | 2.0 |
| HG008T_ONT_d46 | INDEL | 100 | 77.0 | 10.0 | 10.0 | 3.0 | 0.0 |
| HG008T_ONT_d46 | SNV | 50 | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| COLO829T_Illumina_d46 | INDEL | 100 | 98.0 | 0.0 | 2.0 | 0.0 | 0.0 |
| COLO829T_Illumina_d46 | SNV | 50 | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| COLO829T_fiberseq_d46 | INDEL | 100 | 90.0 | 6.0 | 4.0 | 0.0 | 0.0 |
| COLO829T_fiberseq_d46 | SNV | 50 | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| COLO829T_ONT_d46 | INDEL | 100 | 77.0 | 13.0 | 7.0 | 3.0 | 0.0 |
| COLO829T_ONT_d46 | SNV | 50 | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
