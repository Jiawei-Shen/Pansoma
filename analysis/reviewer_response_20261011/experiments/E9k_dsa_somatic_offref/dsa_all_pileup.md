# E9k part 5b: ALT reads at the DSA somatic SNVs in the COLO829T GRCh38 BAMs (dsa_all_tensors.py pileup)

Reads at each SNV's GRCh38 position in the BAM of the reads the tensor set was built from (same sample):
Fiber-seq `SMHTCOLO829T-X-X-M45-B004-uwsc-MERGED-pbmm2_1.13.0_GRCh38.aligned.sorted.bam`, Illumina `SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIZO17Q9I-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.bam`. Bases as pysam count_coverage, no base
quality threshold; MAPQ >= 10 = the tensor builder's min_mapq. Tensor = a tensor with this SNV as A1 in the v1.1 or
v2.1 set of that platform.

## fiberseq

| SNVs | reads at the SNV | SNVs | with a tensor (v1.1 or v2.1 set) |
|---|---|---:|---:|
| lifts, same base, in the SMaHT GRCh38 truth (random 2,000) (2,000) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 1,993 (99.7%) | 1,988 (100%) |
| lifts, same base, in the SMaHT GRCh38 truth (random 2,000) (2,000) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 6 (0.3%) | 0 (0%) |
| lifts, same base, in the SMaHT GRCh38 truth (random 2,000) (2,000) | ALT only in reads with MAPQ < 10 | 1 (0.1%) | 0 (0%) |
| lifts, same base, not in the truth (4,160) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 3,083 (74.1%) | 1,125 (36%) |
| lifts, same base, not in the truth (4,160) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 214 (5.1%) | 84 (39%) |
| lifts, same base, not in the truth (4,160) | ALT only in reads with MAPQ < 10 | 238 (5.7%) | 22 (9%) |
| lifts, same base, not in the truth (4,160) | no ALT read | 591 (14.2%) | 64 (11%) |
| lifts, same base, not in the truth (4,160) | no reads | 34 (0.8%) | 1 (3%) |
| lifts, on a germline allele (148) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 121 (81.8%) | 31 (26%) |
| lifts, on a germline allele (148) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 3 (2.0%) | 1 (33%) |
| lifts, on a germline allele (148) | ALT only in reads with MAPQ < 10 | 14 (9.5%) | 0 (0%) |
| lifts, on a germline allele (148) | no ALT read | 10 (6.8%) | 1 (10%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 231 (11.3%) | 126 (55%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 1,097 (53.6%) | 36 (3%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | ALT only in reads with MAPQ < 10 | 186 (9.1%) | 17 (9%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | no ALT read | 501 (24.5%) | 17 (3%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | no reads | 32 (1.6%) | 1 (3%) |

## illumina

| SNVs | reads at the SNV | SNVs | with a tensor (v1.1 or v2.1 set) |
|---|---|---:|---:|
| lifts, same base, in the SMaHT GRCh38 truth (random 2,000) (2,000) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 1,995 (99.8%) | 1,992 (100%) |
| lifts, same base, in the SMaHT GRCh38 truth (random 2,000) (2,000) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 5 (0.2%) | 0 (0%) |
| lifts, same base, not in the truth (4,160) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 2,801 (67.3%) | 996 (36%) |
| lifts, same base, not in the truth (4,160) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 201 (4.8%) | 78 (39%) |
| lifts, same base, not in the truth (4,160) | ALT only in reads with MAPQ < 10 | 491 (11.8%) | 72 (15%) |
| lifts, same base, not in the truth (4,160) | no ALT read | 648 (15.6%) | 105 (16%) |
| lifts, same base, not in the truth (4,160) | no reads | 19 (0.5%) | 4 (21%) |
| lifts, on a germline allele (148) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 116 (78.4%) | 38 (33%) |
| lifts, on a germline allele (148) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 4 (2.7%) | 0 (0%) |
| lifts, on a germline allele (148) | ALT only in reads with MAPQ < 10 | 16 (10.8%) | 1 (6%) |
| lifts, on a germline allele (148) | no ALT read | 6 (4.1%) | 0 (0%) |
| lifts, on a germline allele (148) | no reads | 6 (4.1%) | 1 (17%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | ALT >= 3 reads and >= 8% at MAPQ >= 10 | 330 (16.1%) | 226 (68%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | ALT 1-2 reads or < 8% at MAPQ >= 10 | 1,619 (79.1%) | 44 (3%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | ALT only in reads with MAPQ < 10 | 51 (2.5%) | 3 (6%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | no ALT read | 45 (2.2%) | 9 (20%) |
| DSA only, class A on d9 (GRCh38 position of the graph node) (2,047) | no reads | 2 (0.1%) | 0 (0%) |

