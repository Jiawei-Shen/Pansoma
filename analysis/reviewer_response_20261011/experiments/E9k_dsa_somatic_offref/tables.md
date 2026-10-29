# E9k: COLO829T graph-coordinate PASS SNV calls against the COLO829BL-DSA somatic SNVs (match_dsa_somatic.py)

DSA somatic VCF: /wanglab/jshen/scratch_jshen/data/COLO829T/somatic_snvs.vcf.gz (44,810 SNV alleles). Calls: E9c contexts aligned to DSA hapX / hapY (E9d PAF).

| graph | platform | group | calls (of unplaced PASS) | DSA somatic (PASS) | normal REF, not in VCF | normal carries ALT | not placed |
|---|---|---|---:|---|---:|---:|---:|
| v1.1 (d9) | Illumina | unplaced_pass | 693 (of 693) | 15 (15) | 342 | 300 | 36 |
| v1.1 (d9) | Illumina | gold | 2 | 0 (0) | 0 | 1 | 1 |
| v1.1 (d9) | Illumina | ctrl_tp | 150 | 123 (123) | 24 | 3 | 0 |
| v1.1 (d9) | Illumina | ctrl_fp | 150 | 13 (13) | 100 | 36 | 1 |
| v1.1 (d9) | Illumina | ctrl_pon | 150 | 1 (1) | 135 | 14 | 0 |
| v1.1 (d9) | Illumina | unplaced_pass, rescue off | 35 (of 38) | 3 (3) | 25 | 7 | 0 |
| v1.1 (d9) | Illumina | norescue_pass | 38 (of 38) | 3 (3) | 28 | 7 | 0 |
| v1.1 (d9) | Fiber-seq | unplaced_pass | 1,483 (of 1,483) | 7 (7) | 541 | 761 | 174 |
| v1.1 (d9) | Fiber-seq | gold | 1 | 0 (0) | 0 | 0 | 1 |
| v1.1 (d9) | Fiber-seq | ctrl_tp | 150 | 123 (123) | 22 | 4 | 1 |
| v1.1 (d9) | Fiber-seq | ctrl_fp | 150 | 12 (12) | 103 | 34 | 1 |
| v1.1 (d9) | Fiber-seq | ctrl_pon | 150 | 2 (2) | 137 | 11 | 0 |
| v1.1 (d9) | Fiber-seq | unplaced_pass, rescue off | 18 (of 35) | 2 (2) | 16 | 0 | 0 |
| v1.1 (d9) | Fiber-seq | norescue_pass | 35 (of 35) | 2 (2) | 32 | 1 | 0 |
| v1.1 (d9) | ONT | unplaced_pass | 3,000 (of 3,421) | 12 (12) | 1733 | 983 | 272 |
| v1.1 (d9) | ONT | gold | 2 | 0 (0) | 0 | 1 | 1 |
| v1.1 (d9) | ONT | ctrl_tp | 150 | 128 (128) | 20 | 1 | 1 |
| v1.1 (d9) | ONT | ctrl_fp | 150 | 16 (16) | 108 | 17 | 9 |
| v1.1 (d9) | ONT | ctrl_pon | 150 | 4 (4) | 130 | 16 | 0 |
| v1.1 (d9) | ONT | unplaced_pass, rescue off | 79 (of 152) | 3 (3) | 75 | 1 | 0 |
| v1.1 (d9) | ONT | norescue_pass | 152 (of 152) | 8 (8) | 137 | 7 | 0 |
| v2.1 (d46) | Illumina | unplaced_pass | 516 (of 516) | 19 (19) | 255 | 236 | 6 |
| v2.1 (d46) | Illumina | gold | 1 | 0 (0) | 0 | 0 | 1 |
| v2.1 (d46) | Illumina | ctrl_tp | 150 | 130 (130) | 18 | 1 | 1 |
| v2.1 (d46) | Illumina | ctrl_fp | 150 | 17 (17) | 115 | 17 | 1 |
| v2.1 (d46) | Illumina | ctrl_pon | 150 | 10 (10) | 133 | 7 | 0 |
| v2.1 (d46) | Illumina | unplaced_pass, rescue off | 17 (of 18) | 4 (4) | 12 | 1 | 0 |
| v2.1 (d46) | Illumina | norescue_pass | 18 (of 18) | 4 (4) | 13 | 1 | 0 |
| v2.1 (d46) | Fiber-seq | unplaced_pass | 257 (of 257) | 7 (7) | 90 | 112 | 48 |
| v2.1 (d46) | Fiber-seq | gold | 2 | 0 (0) | 1 | 0 | 1 |
| v2.1 (d46) | Fiber-seq | ctrl_tp | 150 | 131 (131) | 17 | 2 | 0 |
| v2.1 (d46) | Fiber-seq | ctrl_fp | 150 | 13 (13) | 104 | 32 | 1 |
| v2.1 (d46) | Fiber-seq | ctrl_pon | 150 | 7 (7) | 134 | 8 | 1 |
| v2.1 (d46) | Fiber-seq | unplaced_pass, rescue off | 6 (of 78) | 1 (1) | 4 | 1 | 0 |
| v2.1 (d46) | Fiber-seq | norescue_pass | 78 (of 78) | 2 (2) | 48 | 27 | 1 |
| v2.1 (d46) | ONT | unplaced_pass | 3,000 (of 4,569) | 13 (13) | 1957 | 846 | 184 |
| v2.1 (d46) | ONT | gold | 2 | 0 (0) | 1 | 0 | 1 |
| v2.1 (d46) | ONT | ctrl_tp | 150 | 137 (137) | 13 | 0 | 0 |
| v2.1 (d46) | ONT | ctrl_fp | 150 | 10 (10) | 121 | 15 | 4 |
| v2.1 (d46) | ONT | ctrl_pon | 150 | 2 (2) | 139 | 8 | 1 |
| v2.1 (d46) | ONT | unplaced_pass, rescue off | 61 (of 98) | 3 (3) | 50 | 6 | 2 |
| v2.1 (d46) | ONT | norescue_pass | 98 (of 98) | 4 (4) | 82 | 10 | 2 |

## Graph-coordinate calls that match a DSA somatic SNV (96 call records)

| graph | platform | group | candidate id | rescue off PASS | DSA site (contig:pos:REF>ALT:FILTER) |
|---|---|---|---|---|---|
| v1.1 (d9) | Illumina | unplaced_pass | `51282:154:SNP:C>A` | False | haplotype1-0000012:952820:G>T:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `132995:1:SNP:G>T` | False | haplotype1-0000012:2312940:G>T:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `145516:29:SNP:G>A` | True | haplotype1-0000012:3748081:C>T:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `344973:0:SNP:C>A` | False | haplotype1-0000012:12579751:G>A:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `528467:0:SNP:C>T` | False | haplotype1-0000012:42627451:G>A:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `1091012:4:SNP:T>G` | False | haplotype1-0000012:75785052:T>G:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `2022473:0:SNP:G>T` | False | haplotype1-0000012:174545094:G>T:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `2270420:1:SNP:C>T` | False | haplotype1-0000012:188837834:G>A:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `2298565:74:SNP:G>A` | False | haplotype1-0000012:190107661:G>A:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `2589149:43:SNP:G>A` | False | haplotype1-0000012:2039002:G>A:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `2600934:10:SNP:G>A` | False | haplotype1-0000012:12891506:G>A:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `2846601:11:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `3327532:14:SNP:C>T` | True | haplotype1-0000012:17207074:C>T:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `3383389:0:SNP:G>A` | False | haplotype1-0000012:2083800:C>T:PASS |
| v1.1 (d9) | Illumina | unplaced_pass | `4156612:143:SNP:G>A` | False | haplotype1-0000012:204346846:G>A:PASS |
| v1.1 (d9) | Illumina | norescue_pass | `145516:29:SNP:G>A` | True | haplotype1-0000012:3748081:C>T:PASS |
| v1.1 (d9) | Illumina | norescue_pass | `2846601:11:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v1.1 (d9) | Illumina | norescue_pass | `3327532:14:SNP:C>T` | True | haplotype1-0000012:17207074:C>T:PASS |
| v1.1 (d9) | Fiber-seq | unplaced_pass | `145516:29:SNP:G>A` | True | haplotype1-0000012:3748081:C>T:PASS |
| v1.1 (d9) | Fiber-seq | unplaced_pass | `1091012:4:SNP:T>G` | False | haplotype1-0000012:75785052:T>G:PASS |
| v1.1 (d9) | Fiber-seq | unplaced_pass | `2708908:0:SNP:A>T` | False | haplotype1-0000012:17019427:G>T:PASS |
| v1.1 (d9) | Fiber-seq | unplaced_pass | `2846606:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v1.1 (d9) | Fiber-seq | unplaced_pass | `3328476:77:SNP:C>T` | False | haplotype1-0000012:17207074:C>T:PASS |
| v1.1 (d9) | Fiber-seq | unplaced_pass | `3377740:0:SNP:G>T` | False | haplotype1-0000012:2114688:C>A:PASS |
| v1.1 (d9) | Fiber-seq | unplaced_pass | `3383389:0:SNP:G>A` | False | haplotype1-0000012:2083800:C>T:PASS |
| v1.1 (d9) | Fiber-seq | norescue_pass | `145516:29:SNP:G>A` | True | haplotype1-0000012:3748081:C>T:PASS |
| v1.1 (d9) | Fiber-seq | norescue_pass | `2846606:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `145516:29:SNP:G>A` | True | haplotype1-0000012:3748081:C>T:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `1091012:4:SNP:T>G` | False | haplotype1-0000012:75785052:T>G:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `1592349:90:SNP:C>T` | False | haplotype1-0000012:104237577:G>A:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `2022473:0:SNP:G>T` | False | haplotype1-0000012:174545094:G>T:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `2600934:10:SNP:G>A` | False | haplotype1-0000012:12891506:G>A:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `2708908:0:SNP:A>T` | False | haplotype1-0000012:17019427:G>T:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `2846601:11:SNP:G>T` | False | haplotype1-0000012:23310365:G>T:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `2846606:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `2984582:77:SNP:G>A` | False | haplotype1-0000012:12603243:G>A:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `2986397:44:SNP:G>A` | False | haplotype1-0000012:12695168:G>A:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `3383389:0:SNP:G>A` | True | haplotype1-0000012:2083800:C>T:PASS |
| v1.1 (d9) | ONT | unplaced_pass | `3385904:436:SNP:C>T` | False | haplotype1-0000012:104167618:G>A:PASS |
| v1.1 (d9) | ONT | norescue_pass | `132995:1:SNP:G>T` | True | haplotype1-0000012:2312940:G>T:PASS |
| v1.1 (d9) | ONT | norescue_pass | `145516:29:SNP:G>A` | True | haplotype1-0000012:3748081:C>T:PASS |
| v1.1 (d9) | ONT | norescue_pass | `2846606:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v1.1 (d9) | ONT | norescue_pass | `3383389:0:SNP:G>A` | True | haplotype1-0000012:2083800:C>T:PASS |
| v1.1 (d9) | ONT | norescue_pass | `3878398:18:SNP:G>C` | True | haplotype1-0000012:161457702:G>C:PASS |
| v1.1 (d9) | ONT | norescue_pass | `4007912:315:SNP:G>T` | True | haplotype1-0000012:144748380:G>T:PASS |
| v1.1 (d9) | ONT | norescue_pass | `4008159:7:SNP:C>A` | True | haplotype1-0000012:144756150:C>A:PASS |
| v1.1 (d9) | ONT | norescue_pass | `4008362:21:SNP:C>A` | True | haplotype1-0000012:144756150:C>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `1024740:43:SNP:G>T` | False | haplotype1-0000012:952820:G>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `1279908:247:SNP:C>T` | False | haplotype1-0000012:2083800:C>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `1331451:14:SNP:C>A` | True | haplotype1-0000012:2114688:C>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `1421774:68:SNP:G>T` | False | haplotype1-0000012:2312940:G>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `1651991:2:SNP:C>T` | True | haplotype1-0000012:3748081:C>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `1680255:6:SNP:G>A` | False | haplotype1-0000012:4253456:G>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `1801552:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `2613414:0:SNP:G>A` | False | haplotype1-0000012:42627451:G>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `8189379:74:SNP:G>A` | False | haplotype1-0000012:190107661:G>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `8767546:143:SNP:G>A` | False | haplotype1-0000012:204346846:G>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `10357341:0:SNP:T>A` | False | haplotype1-0000012:238371661:C>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `10793319:10:SNP:G>A` | False | haplotype1-0000012:12891506:G>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `11745251:0:SNP:G>T` | False | haplotype1-0000012:247672629:G>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `11931074:14:SNP:C>T` | True | haplotype1-0000012:17207074:C>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `12066148:76:SNP:C>T` | False | haplotype1-0000012:120885398:C>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `13604241:9:SNP:C>T` | False | haplotype1-0000012:122082190:C>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `15373031:14:SNP:C>T` | False | haplotype1-0000012:128606674:C>T:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `15393890:14:SNP:C>A` | False | haplotype1-0000012:134744900:C>A:PASS |
| v2.1 (d46) | Illumina | unplaced_pass | `19284253:0:SNP:G>T` | False | haplotype1-0000012:17019427:G>T:PASS |
| v2.1 (d46) | Illumina | norescue_pass | `1331451:14:SNP:C>A` | True | haplotype1-0000012:2114688:C>A:PASS |
| v2.1 (d46) | Illumina | norescue_pass | `1651991:2:SNP:C>T` | True | haplotype1-0000012:3748081:C>T:PASS |
| v2.1 (d46) | Illumina | norescue_pass | `1801552:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v2.1 (d46) | Illumina | norescue_pass | `11931074:14:SNP:C>T` | True | haplotype1-0000012:17207074:C>T:PASS |
| v2.1 (d46) | Fiber-seq | unplaced_pass | `1024740:43:SNP:G>T` | False | haplotype1-0000012:952820:G>T:PASS |
| v2.1 (d46) | Fiber-seq | unplaced_pass | `1651991:2:SNP:C>T` | False | haplotype1-0000012:3748081:C>T:PASS |
| v2.1 (d46) | Fiber-seq | unplaced_pass | `1680255:6:SNP:G>A` | False | haplotype1-0000012:4253456:G>A:PASS |
| v2.1 (d46) | Fiber-seq | unplaced_pass | `1801552:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v2.1 (d46) | Fiber-seq | unplaced_pass | `10760510:66:SNP:G>A` | False | haplotype1-0000012:12603243:G>A:PASS |
| v2.1 (d46) | Fiber-seq | unplaced_pass | `11935710:14:SNP:G>A` | False | haplotype1-0000012:17228294:G>A:PASS |
| v2.1 (d46) | Fiber-seq | unplaced_pass | `13604241:9:SNP:C>T` | False | haplotype1-0000012:122082190:C>T:PASS |
| v2.1 (d46) | Fiber-seq | norescue_pass | `1801552:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v2.1 (d46) | Fiber-seq | norescue_pass | `5728710:11:SNP:C>A` | True | haplotype1-0000012:144756150:C>A:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `976969:1:SNP:C>T` | False | haplotype1-0000012:642469:C>T:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `1024740:43:SNP:G>T` | False | haplotype1-0000012:952820:G>T:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `1331451:14:SNP:C>A` | True | haplotype1-0000012:2114688:C>A:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `1651991:2:SNP:C>T` | True | haplotype1-0000012:3748081:C>T:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `1801552:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `5016946:154:SNP:G>A` | False | haplotype1-0000012:104167618:G>A:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `5723807:21:SNP:G>A` | False | haplotype1-0000012:144751337:G>A:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `8189379:74:SNP:G>A` | False | haplotype1-0000012:190107661:G>A:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `10793319:10:SNP:G>A` | False | haplotype1-0000012:12891506:G>A:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `11938772:69:SNP:G>A` | False | haplotype1-0000012:104237577:G>A:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `12066148:76:SNP:C>T` | False | haplotype1-0000012:120885398:C>T:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `13604241:9:SNP:C>T` | False | haplotype1-0000012:122082190:C>T:PASS |
| v2.1 (d46) | ONT | unplaced_pass | `19284253:0:SNP:G>T` | False | haplotype1-0000012:17019427:G>T:PASS |
| v2.1 (d46) | ONT | norescue_pass | `1331451:14:SNP:C>A` | True | haplotype1-0000012:2114688:C>A:PASS |
| v2.1 (d46) | ONT | norescue_pass | `1651991:2:SNP:C>T` | True | haplotype1-0000012:3748081:C>T:PASS |
| v2.1 (d46) | ONT | norescue_pass | `1801552:3:SNP:G>T` | True | haplotype1-0000012:23310365:G>T:PASS |
| v2.1 (d46) | ONT | norescue_pass | `11827337:14:SNP:G>T` | True | haplotype1-0000012:248562822:G>T:PASS |
