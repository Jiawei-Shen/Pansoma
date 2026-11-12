# Truth VAF, Illumina and PacBio (chr1-22, all truth alleles of the label truth VCF)

Quartiles and low-VAF shares over the alleles with reads (DP > 0); COLO829T = SMaHT VAF_Ill / VAF_PB (no DP).

| tumor | platform | type | n | DP = 0 | Q25 | median | Q75 | VAF < 0.05 | VAF < 0.1 |
|---|---|---|---|---|---|---|---|---|---|
| HG008T | Illumina | SNV | 8,690 | 30 | 0.467 | 0.536 | 1.000 | 0.8% | 1.0% |
| HG008T | Illumina | INDEL | 8,482 | 13 | 0.227 | 0.378 | 0.555 | 7.8% | 11.6% |
| HG008T | PacBio | SNV | 8,690 | 38 | 0.467 | 0.544 | 1.000 | 0.7% | 0.8% |
| HG008T | PacBio | INDEL | 8,482 | 12 | 0.318 | 0.429 | 0.566 | 2.3% | 2.9% |
| COLO829T | Illumina | SNV | 42,035 | 0 | 0.332 | 0.497 | 0.662 | 6.1% | 10.6% |
| COLO829T | Illumina | INDEL | 1,912 | 0 | 0.321 | 0.439 | 0.605 | 0.0% | 0.0% |
| COLO829T | PacBio | SNV | 42,035 | 0 | 0.324 | 0.491 | 0.657 | 6.2% | 10.6% |
| COLO829T | PacBio | INDEL | 1,912 | 0 | 0.382 | 0.482 | 0.629 | 0.0% | 0.0% |
| HCC1395T | Illumina | SNV | 39,447 | 0 | 0.170 | 0.341 | 0.565 | 5.3% | 13.7% |
| HCC1395T | Illumina | INDEL | 1,625 | 0 | 0.138 | 0.256 | 0.419 | 5.9% | 15.6% |
| HCC1395T | PacBio | SNV | 39,447 | 8 | 0.167 | 0.348 | 0.569 | 6.3% | 15.2% |
| HCC1395T | PacBio | INDEL | 1,625 | 0 | 0.183 | 0.341 | 0.512 | 5.7% | 12.9% |

HCC1395T Illumina library check, chr1 truth alleles with reads in both: SRR7890824 (used above) vs the WGS_FD_T_1 + 2 reads Pansoma uses (E8 HCC1395T_Illumina.chr1_chr20.bam), same counting:
n 3,573, Pearson r 0.984, median VAF 0.352 vs 0.352, median DP 53 vs 105.
