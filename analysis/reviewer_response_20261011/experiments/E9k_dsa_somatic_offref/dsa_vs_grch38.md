# E9k part 2: the COLO829BL-DSA somatic SNVs on GRCh38 (dsa_vs_grch38.py)

DSA somatic SNVs: 44,810 (/wanglab/jshen/scratch_jshen/data/COLO829T/somatic_snvs.vcf.gz); lift: /scratch/qfu/COLO829BL_DSA/COLO829BL_hg38.chain; GRCh38 truth: SMaHT COLO829T (42,035 autosomal SNVs).

| class | DSA somatic SNVs | in the SMaHT GRCh38 truth |
|---|---:|---:|
| lifts to GRCh38, same reference base | 39,310 (87.7%) | 35,150 |
| lifts, GRCh38 base differs (on a germline allele) | 148 (0.3%) | 0 |
| no GRCh38 position in the chain (DSA only) | 5,352 (11.9%) | 0 |
| all | 44,810 | 35,150 |

Lifted to an autosome: 37,536; the SMaHT GRCh38 truth has 42,035 autosomal SNVs, 35,150 of which (83.6%) are DSA somatic SNVs by this lift.

## The 35 DSA somatic sites that Pansoma calls in graph coordinates (part 1; * = also PASS without the rescue)

| DSA site | change | class | GRCh38 | in the SMaHT GRCh38 truth | calls |
|---|---|---|---|---|---|
| haplotype1-0000012:642469 | C>T | lifts to GRCh38, same reference base | chr1:1288697:+ | False | v2.1 ONT |
| haplotype1-0000012:952820 | G>T | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina, v2.1 Illumina, v2.1 Fiber-seq, v2.1 ONT |
| haplotype1-0000012:2039002 | G>A | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina |
| haplotype1-0000012:2083800 | C>T | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina, v1.1 Fiber-seq, v1.1 ONT*, v2.1 Illumina |
| haplotype1-0000012:2114688 | C>A | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Fiber-seq, v2.1 Illumina*, v2.1 ONT* |
| haplotype1-0000012:2312940 | G>T | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina, v1.1 ONT*, v2.1 Illumina |
| haplotype1-0000012:3748081 | C>T | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina*, v1.1 Fiber-seq*, v1.1 ONT*, v2.1 Illumina*, v2.1 Fiber-seq, v2.1 ONT* |
| haplotype1-0000012:4253456 | G>A | lifts to GRCh38, same reference base | chr1:4666546:+ | False | v2.1 Illumina, v2.1 Fiber-seq |
| haplotype1-0000012:12579751 | G>A | lifts to GRCh38, same reference base | chr1:12873251:+ | True | v1.1 Illumina |
| haplotype1-0000012:12603243 | G>A | lifts to GRCh38, same reference base | chr1:12896733:+ | True | v1.1 ONT, v2.1 Fiber-seq |
| haplotype1-0000012:12695168 | G>A | lifts to GRCh38, same reference base | chr1:13015854:+ | False | v1.1 ONT |
| haplotype1-0000012:12891506 | G>A | lifts to GRCh38, same reference base | chr1:13239715:+ | False | v1.1 Illumina, v1.1 ONT, v2.1 Illumina, v2.1 ONT |
| haplotype1-0000012:17019427 | G>T | lifts to GRCh38, same reference base | chr1_KI270713v1_random:10988:- | False | v1.1 Fiber-seq, v1.1 ONT, v2.1 Illumina, v2.1 ONT |
| haplotype1-0000012:17207074 | C>T | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina*, v1.1 Fiber-seq, v2.1 Illumina* |
| haplotype1-0000012:17228294 | G>A | lifts to GRCh38, same reference base | chr1_KI270714v1_random:15466:+ | False | v2.1 Fiber-seq |
| haplotype1-0000012:23310365 | G>T | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina*, v1.1 Fiber-seq*, v1.1 ONT, v1.1 ONT*, v2.1 Illumina*, v2.1 Fiber-seq*, v2.1 ONT* |
| haplotype1-0000012:42627451 | G>A | lifts, GRCh38 base differs (on a germline allele) | chr1:42095766:+ | False | v1.1 Illumina, v2.1 Illumina |
| haplotype1-0000012:75785052 | T>G | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina, v1.1 Fiber-seq, v1.1 ONT |
| haplotype1-0000012:104167618 | G>A | lifts to GRCh38, same reference base | chr1:103660049:+ | False | v1.1 ONT, v2.1 ONT |
| haplotype1-0000012:104237577 | G>A | lifts to GRCh38, same reference base | chr1:103730020:+ | False | v1.1 ONT, v2.1 ONT |
| haplotype1-0000012:120885398 | C>T | lifts to GRCh38, same reference base | chr1:120413508:- | True | v2.1 Illumina, v2.1 ONT |
| haplotype1-0000012:122082190 | C>T | lifts to GRCh38, same reference base | chr1:121634506:+ | True | v2.1 Illumina, v2.1 Fiber-seq, v2.1 ONT |
| haplotype1-0000012:128606674 | C>T | no GRCh38 position in the chain (DSA only) | - | False | v2.1 Illumina |
| haplotype1-0000012:134744900 | C>A | no GRCh38 position in the chain (DSA only) | - | False | v2.1 Illumina |
| haplotype1-0000012:144748380 | G>T | no GRCh38 position in the chain (DSA only) | - | False | v1.1 ONT* |
| haplotype1-0000012:144751337 | G>A | no GRCh38 position in the chain (DSA only) | - | False | v2.1 ONT |
| haplotype1-0000012:144756150 | C>A | no GRCh38 position in the chain (DSA only) | - | False | v1.1 ONT*, v2.1 Fiber-seq* |
| haplotype1-0000012:161457702 | G>C | no GRCh38 position in the chain (DSA only) | - | False | v1.1 ONT* |
| haplotype1-0000012:174545094 | G>T | lifts, GRCh38 base differs (on a germline allele) | chr1:174474598:+ | False | v1.1 Illumina, v1.1 ONT |
| haplotype1-0000012:188837834 | G>A | lifts to GRCh38, same reference base | chr1:188764300:+ | False | v1.1 Illumina |
| haplotype1-0000012:190107661 | G>A | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina, v2.1 Illumina, v2.1 ONT |
| haplotype1-0000012:204346846 | G>A | no GRCh38 position in the chain (DSA only) | - | False | v1.1 Illumina, v2.1 Illumina |
| haplotype1-0000012:238371661 | C>A | lifts to GRCh38, same reference base | chr1:238211589:+ | True | v2.1 Illumina |
| haplotype1-0000012:247672629 | G>T | lifts to GRCh38, same reference base | chr1:247487800:+ | False | v2.1 Illumina |
| haplotype1-0000012:248562822 | G>T | lifts to GRCh38, same reference base | chr1:248388669:+ | True | v2.1 ONT* |

Summary of these sites: no GRCh38 position in the chain (DSA only), not in the GRCh38 truth: 17; lifts, GRCh38 base differs (on a germline allele), not in the GRCh38 truth: 2; lifts to GRCh38, same reference base, not in the GRCh38 truth: 10; lifts to GRCh38, same reference base, in the GRCh38 truth: 6
