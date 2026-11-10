## R1.m2_INDEL_01  HG008T_Illumina  2543944:0:INS:>T  (shard chr1_shard_00000_data.npy index 5869)

- truth chr1:213363821:C>CT; INS ->T; AF 0.283; rows A1 32, REF 78, other 3; H = 90
- candidate node 2543944 (1 bp, 80 paths), site offset 0-0; allele nodes [2543944]; candidate columns [50, 51]
- GRCh38-path nodes under the window: 14: 2543929 (4 bp, 69), 2543930 (4 bp, 67), 2543931 (4 bp, 55), 2543932 (4 bp, 35), 2543933 (4 bp, 22), 2543935 (1 bp, 74), 2543937 (4 bp, 68), 2543938 (4 bp, 71), 2543939 (4 bp, 76), 2543940 (8 bp, 80), 2543941 (3 bp, 44), 2543943 (6 bp, 80), 2543944 (1 bp, 80), 2543945 (254 bp, 89)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 0-3: A1 code 79 = 69 paths (nodes [2543927, 2543928, 2543929]) [91%+, TTCT] | REF code 79 = 69 paths (nodes [2543927, 2543928, 2543929]) [84%+, TTCT]
  - cols 4-7 DIFF: A1 code 79 = 69 paths (nodes [2543927, 2543928, 2543929]) [96%+, TTCT] | REF code 77 = 67 paths (nodes [2543926, 2543930]) [67%+, TTCT]
  - cols 8-11 DIFF: A1 code 77 = 67 paths (nodes [2543926, 2543930]) [83%+, TTCT] | REF code 65 = 55 paths (nodes [2543931]) [67%+, TTCT]
  - cols 12-15 DIFF: A1 code 65 = 55 paths (nodes [2543931]) [83%+, TTCT] | REF code 35 = 35 paths (nodes [2543932]) [64%+, TTCT]
  - cols 16-16 DIFF: A1 code 35 = 35 paths (nodes [2543932]) [83%+, T] | REF code 84 = 74 paths (nodes [2543935]) [75%+, T]
  - cols 17-19 DIFF: A1 code 35 = 35 paths (nodes [2543932]) [80%+, TCT] | REF code 29 = 29 paths (nodes [2543936]) [75%+, TCT]
  - cols 20-20 DIFF: A1 code 84 = 74 paths (nodes [2543935]) [96%+, T] | REF code 29 = 29 paths (nodes [2543936]) [76%+, T]
  - cols 21-24: A1 code 78 = 68 paths (nodes [2543937]) [100%+, TCTT] | REF code 78 = 68 paths (nodes [2543937]) [97%+, TCTT]
  - cols 25-28: A1 code 81 = 71 paths (nodes [2543938]) [100%+, TCTT] | REF code 81 = 71 paths (nodes [2543938]) [97%+, TCTT]
  - cols 29-32: A1 code 86 = 76 paths (nodes [2543939]) [100%+, TCTT] | REF code 86 = 76 paths (nodes [2543939]) [97%+, TCTT]
  - cols 33-40: A1 code 90 = 80 paths (nodes [2543940, 2543943, 2543944]) [100%+, TCTTTCTT] | REF code 90 = 80 paths (nodes [2543940, 2543943, 2543944]) [100%+, TCTTTCTT]
  - cols 41-43 DIFF: A1 code 54 = 44 paths (nodes [2543941]) [100%+, TCT] | REF code 36 = 36 paths (nodes [2543942]) [85%+, CTC]
  - cols 44-49: A1 code 90 = 80 paths (nodes [2543940, 2543943, 2543944]) [100%+, TCTCTC] | REF code 90 = 80 paths (nodes [2543940, 2543943, 2543944]) [100%+, TCTCTC]
  - cols 50-51 SITE: A1 code 90 = 80 paths (nodes [2543940, 2543943, 2543944]) [100%+, -T] | REF code 90 = 80 paths (nodes [2543940, 2543943, 2543944]) [100%+, -T]
  - cols 52-100: A1 code 99 = 89 paths (nodes [2543945]) [100%+, TTTTTTTTTTTTTCTGGTGGCGTGAGGAGATTGAGGGAGGGAAAACCTC] | REF code 99 = 89 paths (nodes [2543945]) [100%+, TTTTTTTTTTTTTCTGGTGGCGTGAGGAGATTGAGGGAGGGAAAACCTC]

## R1.m2_INDEL_02  HG008T_ONT  2264449:0:DEL:T>  (shard chr1_shard_00003_data.npy index 10131)

- truth chr1:188508320:GT>G; DEL T>-; AF 0.565; rows A1 26, REF 13, other 7; H = 90
- candidate node 2264449 (27 bp, 90 paths), site offset 0-1; allele nodes [2264449]; candidate columns [50, 51]
- GRCh38-path nodes under the window: 6: 2264445 (139 bp, 90), 2264449 (27 bp, 90), 2264450 (1 bp, 89), 2264451 (12 bp, 90), 2264452 (1 bp, 31), 2264454 (43 bp, 90)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 45-47: A1 code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [100%+, TTT] | REF code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [54%+, GGT]
  - cols 48-48 DIFF: A1 code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [100%+, G] | REF code 53 = 43 paths (nodes [2264447]) [77%+, T]
  - cols 49-49 DIFF: A1 code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [100%+, G] | REF code 69 = 59 paths (nodes [2264448, 2264453, 2264456]) [92%+, T]
  - cols 77-77: A1 code 99 = 89 paths (nodes [2264444, 2264450]) [100%+, G] | REF code 99 = 89 paths (nodes [2264444, 2264450]) [92%+, G]
  - cols 79-79: A1 code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [100%+, T] | REF code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [92%+, T]
  - cols 90-90 DIFF: A1 code 31 = 31 paths (nodes [2264452, 2264455]) [100%+, T] | REF code 69 = 59 paths (nodes [2264448, 2264453, 2264456]) [92%+, C]
  - cols 92-92: A1 code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [100%+, C] | REF code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [92%+, C]
  - cols 96-96: A1 code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [100%+, C] | REF code 100 = 90 paths (nodes [2264443, 2264445, 2264449, 2264451, 2264454, 2264457]) [100%+, C]

## R1.m2_INDEL_03  HG008T_PacBio  3693578:0:DEL:A>  (shard chr1_shard_00003_data.npy index 31380)

- truth chr1:218103815:CA>C; DEL A>-; AF 0.366; rows A1 45, REF 37, other 16; H = 90
- candidate node 3693578 (13 bp, 88 paths), site offset 0-1; allele nodes [3693578]; candidate columns [50, 53]
- GRCh38-path nodes under the window: 10: 3693574 (82 bp, 90), 3693578 (13 bp, 88), 3693579 (1 bp, 85), 3693580 (4 bp, 88), 3693581 (1 bp, 87), 3693582 (16 bp, 88), 3693583 (1 bp, 85), 3693584 (1 bp, 70), 3693586 (13 bp, 85), 3693587 (4 bp, 82)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 44-45: A1 code 100 = 90 paths (nodes [3693571, 3693574, 3693588]) [100%+, CG] | REF code 100 = 90 paths (nodes [3693571, 3693574, 3693588]) [100%+, TC]
  - cols 46-46 DIFF: A1 code 100 = 90 paths (nodes [3693571, 3693574, 3693588]) [100%+, T] | REF code 24 = 24 paths (nodes [3693575]) [73%+, A]
  - cols 47-49 DIFF: A1 code 100 = 90 paths (nodes [3693571, 3693574, 3693588]) [98%+, CTC] | REF code 69 = 59 paths (nodes [3693576, 3693577]) [84%+, AAA]
  - cols 50-53 SITE: A1 code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, -AAA] | REF code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, -AAA]
  - cols 54-63: A1 code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, AAAAAAAAAA] | REF code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, AAAAAAAAAA]
  - cols 64-64: A1 code 95 = 85 paths (nodes [3693579, 3693583, 3693586]) [100%+, A] | REF code 95 = 85 paths (nodes [3693579, 3693583, 3693586]) [100%+, A]
  - cols 65-68: A1 code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, AAAA] | REF code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, AAAA]
  - cols 69-69: A1 code 97 = 87 paths (nodes [3693581]) [100%+, G] | REF code 97 = 87 paths (nodes [3693581]) [100%+, G]
  - cols 70-85: A1 code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, AAAAGAAAAGAAAAAA] | REF code 98 = 88 paths (nodes [3693578, 3693580, 3693582]) [100%+, AAAAGAAAAGAAAAAA]
  - cols 86-86: A1 code 95 = 85 paths (nodes [3693579, 3693583, 3693586]) [98%+, A] | REF code 95 = 85 paths (nodes [3693579, 3693583, 3693586]) [95%+, A]
  - cols 87-87: A1 code 80 = 70 paths (nodes [3693584]) [98%+, T] | REF code 80 = 70 paths (nodes [3693584]) [92%+, T]
  - cols 88-100: A1 code 95 = 85 paths (nodes [3693579, 3693583, 3693586]) [98%+, ATATATATATATC] | REF code 95 = 85 paths (nodes [3693579, 3693583, 3693586]) [95%+, ATATATATATATC]

## R1.m2_INDEL_04  HG008T_Illumina_d46  19704905:0:INS:>C  (shard chr1_shard_00000_data.npy index 11589)

- truth chr1:159111805:A>AC; INS ->C; AF 0.133; rows A1 19, REF 123, other 1; H = 464
- candidate node 19704905 (1 bp, 446 paths), site offset 0-0; allele nodes [19704905]; candidate columns [50, 51]
- GRCh38-path nodes under the window: 13: 19704893 (58 bp, 457), 19704894 (1 bp, 456), 19704896 (7 bp, 456), 19704901 (1 bp, 419), 19704902 (9 bp, 448), 19704903 (1 bp, 103), 19704905 (1 bp, 446), 19704907 (25 bp, 448), 19704910 (4 bp, 121), 19704911 (4 bp, 317), 19704912 (4 bp, 433), 19704913 (4 bp, 462), 19704914 (34 bp, 464)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 0-30: A1 code 93 = 457 paths (nodes [19704890, 19704893]) [100%+, AGTTCGAGGCCTGGCCAACATGTTGAAACCC] | REF code 93 = 457 paths (nodes [19704890, 19704893]) [65%+, AGTTCGAGGCCTGGCCAACATGTTGAAACCC]
  - cols 31-38: A1 code 92 = 456 paths (nodes [19704891, 19704894, 19704896]) [100%+, GTCTCTAC] | REF code 92 = 456 paths (nodes [19704891, 19704894, 19704896]) [68%+, GTCTCTAC]
  - cols 39-39: A1 code 55 = 419 paths (nodes [19704901]) [100%+, A] | REF code 55 = 419 paths (nodes [19704901]) [99%+, A]
  - cols 40-48: A1 code 84 = 448 paths (nodes [19704902, 19704907]) [100%+, AAAAAAAAA] | REF code 84 = 448 paths (nodes [19704902, 19704907]) [99%+, AAAAAAAAA]
  - cols 49-49 DIFF: A1 code 25 = 96-104 paths (nodes [19704903]) [100%+, A] | REF code 46 = 342-357 paths (nodes [19704904]) [88%+, C]
  - cols 50-51 SITE: A1 code 82 = 446 paths (nodes [19704905]) [100%+, -A] | REF code 82 = 446 paths (nodes [19704905]) [100%+, -A]
  - cols 52-76: A1 code 84 = 448 paths (nodes [19704902, 19704907]) [100%+, AAAAAAAACTATCTATCTATCTATC] | REF code 84 = 448 paths (nodes [19704902, 19704907]) [100%+, AAAAAAAACTATCTATCTATCTATC]
  - cols 77-80 DIFF: A1 code 14 = 28-32 paths (nodes [19704909]) [92%+, ATCT] | REF code 69 = 433 paths (nodes [19704912]) [68%+, ATCT]
  - cols 81-84 DIFF: A1 code 27 = 113-121 paths (nodes [19704910]) [92%+, ATCT] | REF code 98 = 462 paths (nodes [19704913]) [66%+, ATCT]
  - cols 85-88 DIFF: A1 code 44 = 312-326 paths (nodes [19704911]) [86%+, ATCT] | REF code 100 = 464 paths (nodes [19704914, 19704916]) [64%+, ATCT]
  - cols 89-92 DIFF: A1 code 69 = 433 paths (nodes [19704912]) [83%+, ATCT] | REF code 100 = 464 paths (nodes [19704914, 19704916]) [94%+, ATCT]
  - cols 93-96 DIFF: A1 code 98 = 462 paths (nodes [19704913]) [83%+, ATCT] | REF code 100 = 464 paths (nodes [19704914, 19704916]) [98%+, ATCT]
  - cols 97-100: A1 code 100 = 464 paths (nodes [19704914, 19704916]) [100%+, ATCT] | REF code 100 = 464 paths (nodes [19704914, 19704916]) [98%+, ATCT]

## R1.m2_INDEL_05  HG008T_Illumina_d46  10089773:17:INS:>A  (shard chr1_shard_00000_data.npy index 7439)

- truth chr1:233387578:C>CA; INS ->A; AF 0.493; rows A1 99, REF 101, other 0; H = 464
- candidate node 10089773 (19 bp, 464 paths), site offset 17-17; allele nodes [10089773]; candidate columns [50, 51]
- GRCh38-path nodes under the window: 9: 10089768 (112 bp, 464), 10089769 (1 bp, 463), 10089771 (1 bp, 464), 10089772 (3 bp, 463), 10089773 (19 bp, 464), 10089774 (1 bp, 311), 10089776 (2 bp, 464), 10089777 (1 bp, 279), 10089779 (76 bp, 464)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 28-28: A1 code 99 = 463 paths (nodes [10089769, 10089772]) [100%+, C] | REF code 99 = 463 paths (nodes [10089769, 10089772]) [100%+, C]
  - cols 30-32: A1 code 99 = 463 paths (nodes [10089769, 10089772]) [100%+, AAC] | REF code 99 = 463 paths (nodes [10089769, 10089772]) [100%+, AAC]
  - cols 53-53 DIFF: A1 code 31 = 151-161 paths (nodes [10089775]) [100%+, A] | REF code 43 = 298-311 paths (nodes [10089774]) [85%+, C]
  - cols 54-54: A1 code 100 = 464 paths (nodes [10089765, 10089768, 10089771, 10089773, 10089776, 10089779, 10089782]) [100%+, A] | REF code 100 = 464 paths (nodes [10089765, 10089768, 10089771, 10089773, 10089776, 10089779, 10089782]) [99%+, A]
  - cols 56-56 DIFF: A1 code 34 = 183-194 paths (nodes [10089778]) [98%+, A] | REF code 41 = 270-283 paths (nodes [10089777]) [85%+, C]
  - cols 57-88: A1 code 100 = 464 paths (nodes [10089765, 10089768, 10089771, 10089773, 10089776, 10089779, 10089782]) [100%+, AAAACTAGCTATTGTTTTTACAACCAAGCAAA] | REF code 100 = 464 paths (nodes [10089765, 10089768, 10089771, 10089773, 10089776, 10089779, 10089782]) [97%+, AAAACTAGCTATTGTTTTTACAACCAAGCAAA]

## R1.m2_INDEL_06  HG008T_ONT_d46  8181742:24:INS:>A  (shard chr1_shard_00004_data.npy index 296)

- truth chr1:189870594:C>CA; INS ->A; AF 0.304; rows A1 14, REF 29, other 3; H = 464
- candidate node 8181742 (26 bp, 456 paths), site offset 24-24; allele nodes [8181742]; candidate columns [50, 51]
- GRCh38-path nodes under the window: 8: 8181737 (27 bp, 456), 8181738 (1 bp, 453), 8181741 (1 bp, 453), 8181742 (26 bp, 456), 8181743 (1 bp, 344), 8181745 (27 bp, 456), 8181746 (1 bp, 456), 8181747 (20 bp, 464)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 0-23: A1 code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [93%+, ACAATTAATAAAAGATATATTTTT] | REF code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [97%+, ACAATTAATAAAAGATATATTTTT]
  - cols 24-25: A1 code 89 = 453 paths (nodes [8181738, 8181741]) [86%+, TA] | REF code 89 = 453 paths (nodes [8181738, 8181741]) [97%+, TA]
  - cols 26-49: A1 code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [100%+, AAAAAAATAGTGAACAACAACAAC] | REF code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [100%+, AAAAAAATAGTGAACAACAACAAC]
  - cols 50-51 SITE: A1 code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [100%+, -A] | REF code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [100%+, -A]
  - cols 52-52: A1 code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [100%+, A] | REF code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [100%+, A]
  - cols 53-53 DIFF: A1 code 26 = 105-112 paths (nodes [8181744]) [100%+, A] | REF code 46 = 342-357 paths (nodes [8181743]) [97%+, C]
  - cols 54-81: A1 code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [100%+, AAAAAGCACAAAGACTTCCCTCTTAGGA] | REF code 92 = 456 paths (nodes [8181734, 8181737, 8181742, 8181745, 8181746]) [97%+, AAAAAGCACAAAGACTTCCCTCTTAGGA]
  - cols 82-91: A1 code 100 = 464 paths (nodes [8181747, 8181750]) [79%+, AAAAAAAAAA] | REF code 100 = 464 paths (nodes [8181747, 8181750]) [72%+, AAAAAAAAAA]

## R1.m2_INDEL_07  HG008T_ONT_d46  8237854:0:INS:>TA  (shard chr1_shard_00004_data.npy index 1350)

- truth chr1:191234698:T>TTA; INS ->TA; AF 0.292; rows A1 14, REF 33, other 1; H = 464
- candidate node 8237854 (1 bp, 414 paths), site offset 0-0; allele nodes [8237854]; candidate columns [50, 52]
- GRCh38-path nodes under the window: 46: 8237734 (6 bp, 130), 8237738 (1 bp, 113), 8237739 (1 bp, 111), 8237807 (1 bp, 363), 8237808 (1 bp, 340), 8237810 (5 bp, 366), 8237811 (1 bp, 366), 8237812 (3 bp, 366), 8237813 (1 bp, 267), 8237814 (1 bp, 186), 8237817 (1 bp, 199), 8237819 (1 bp, 365), 8237822 (1 bp, 345), 8237823 (7 bp, 353), 8237824 (1 bp, 249), 8237825 (1 bp, 209), 8237826 (1 bp, 264), 8237831 (1 bp, 371), 8237835 (1 bp, 202), 8237837 (1 bp, 203), 8237840 (1 bp, 191), 8237841 (1 bp, 109), 8237843 (1 bp, 204), 8237845 (8 bp, 209), 8237846 (1 bp, 302), 8237847 (1 bp, 300), 8237849 (1 bp, 264), 8237854 (1 bp, 414), 8237855 (1 bp, 415), 8237856 (8 bp, 444), 8237857 (1 bp, 433), 8237861 (1 bp, 432), 8237862 (1 bp, 419), 8237873 (1 bp, 146), 8237874 (1 bp, 414), 8237876 (1 bp, 417), 8237878 (1 bp, 427), 8237879 (6 bp, 428), 8237880 (1 bp, 342), 8237888 (1 bp, 443), 8237928 (4 bp, 463), 8237930 (1 bp, 463), 8237934 (5 bp, 463), 8237935 (10 bp, 464), 8237936 (1 bp, 463), 8237938 (67 bp, 464)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 0-4 DIFF: A1 code 28 = 122-131 paths (nodes [8237734]) [100%+, TATAT] | REF code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [79%+, TATAT]
  - cols 5-5 DIFF: A1 code 27 = 113-121 paths (nodes [8237738, 8237887, 8237891]) [100%+, A] | REF code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [79%+, A]
  - cols 6-6 DIFF: A1 code 26 = 105-112 paths (nodes [8237739, 8237829, 8237830, 8237841, 8237885, 8237890]) [100%+, T] | REF code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [76%+, T]
  - cols 7-7 DIFF: A1 code 47 = 358-372 paths (nodes [8237807, 8237810, 8237811, 8237812, 8237819, 8237831]) [100%+, T] | REF code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [76%+, A]
  - cols 8-8 DIFF: A1 code 45 = 327-341 paths (nodes [8237808]) [100%+, A] | REF code 39 = 243-255 paths (nodes [8237824]) [76%+, T]
  - cols 9-9 DIFF: A1 code 47 = 358-372 paths (nodes [8237807, 8237810, 8237811, 8237812, 8237819, 8237831]) [100%+, T] | REF code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [76%+, A]
  - cols 10-10 DIFF: A1 code 47 = 358-372 paths (nodes [8237807, 8237810, 8237811, 8237812, 8237819, 8237831]) [100%+, A] | REF code 23 = 81-88 paths (nodes [8237815, 8237828, 8237836, 8237865]) [76%+, A]
  - cols 11-17 DIFF: A1 code 47 = 358-372 paths (nodes [8237807, 8237810, 8237811, 8237812, 8237819, 8237831]) [100%+, TATATAT] | REF code 26 = 105-112 paths (nodes [8237739, 8237829, 8237830, 8237841, 8237885, 8237890]) [79%+, TATATAT]
  - cols 18-18 DIFF: A1 code 40 = 256-269 paths (nodes [8237813, 8237826, 8237849, 8237931, 8237933]) [100%+, A] | REF code 26 = 105-112 paths (nodes [8237739, 8237829, 8237830, 8237841, 8237885, 8237890]) [76%+, A]
  - cols 19-19 DIFF: A1 code 34 = 183-194 paths (nodes [8237814, 8237818, 8237840]) [100%+, T] | REF code 26 = 105-112 paths (nodes [8237739, 8237829, 8237830, 8237841, 8237885, 8237890]) [85%+, T]
  - cols 20-20 DIFF: A1 code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [100%+, A] | REF code 47 = 358-372 paths (nodes [8237807, 8237810, 8237811, 8237812, 8237819, 8237831]) [85%+, A]
  - cols 21-21 DIFF: A1 code 47 = 358-372 paths (nodes [8237807, 8237810, 8237811, 8237812, 8237819, 8237831]) [100%+, A] | REF code 18 = 48-53 paths (nodes [8237736, 8237832, 8237833]) [79%+, T]
  - cols 22-24 DIFF: A1 code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [100%+, TAT] | REF code 18 = 48-53 paths (nodes [8237736, 8237832, 8237833]) [76%+, ATA]
  - cols 25-29 DIFF: A1 code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [100%+, ATATA] | REF code 20 = 60-66 paths (nodes [8237834]) [79%+, ATATA]
  - cols 30-30 DIFF: A1 code 39 = 243-255 paths (nodes [8237824]) [100%+, T] | REF code 20 = 60-66 paths (nodes [8237834]) [76%+, T]
  - cols 31-31 DIFF: A1 code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [100%+, A] | REF code 20 = 60-66 paths (nodes [8237834]) [76%+, A]
  - cols 32-32 DIFF: A1 code 40 = 256-269 paths (nodes [8237813, 8237826, 8237849, 8237931, 8237933]) [100%+, T] | REF code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [76%+, T]
  - cols 33-33 DIFF: A1 code 47 = 358-372 paths (nodes [8237807, 8237810, 8237811, 8237812, 8237819, 8237831]) [100%+, A] | REF code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [85%+, A]
  - cols 34-34 DIFF: A1 code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [100%+, T] | REF code 34 = 183-194 paths (nodes [8237814, 8237818, 8237840]) [76%+, T]
  - cols 35-35 DIFF: A1 code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [100%+, A] | REF code 24 = 89-95 paths (nodes [8237842, 8237850, 8237851]) [76%+, A]
  - cols 36-36 DIFF: A1 code 34 = 183-194 paths (nodes [8237814, 8237818, 8237840]) [100%+, T] | REF code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [79%+, A]
  - cols 37-37 DIFF: A1 code 26 = 105-112 paths (nodes [8237739, 8237829, 8237830, 8237841, 8237885, 8237890]) [100%+, T] | REF code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [79%+, T]
  - cols 38-38 DIFF: A1 code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [100%+, A] | REF code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [79%+, A]
  - cols 39-44: A1 code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [100%+, TATATA] | REF code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [91%+, TATATA]
  - cols 45-46 DIFF: A1 code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [100%+, TA] | REF code 43 = 298-311 paths (nodes [8237846, 8237847, 8237939]) [79%+, TA]
  - cols 47-47 DIFF: A1 code 43 = 298-311 paths (nodes [8237846, 8237847, 8237939]) [100%+, T] | REF code 40 = 256-269 paths (nodes [8237813, 8237826, 8237849, 8237931, 8237933]) [79%+, T]
  - cols 48-48 DIFF: A1 code 43 = 298-311 paths (nodes [8237846, 8237847, 8237939]) [100%+, A] | REF code 24 = 89-95 paths (nodes [8237842, 8237850, 8237851]) [79%+, A]
  - cols 49-49 DIFF: A1 code 40 = 256-269 paths (nodes [8237813, 8237826, 8237849, 8237931, 8237933]) [100%+, T] | REF code 24 = 89-95 paths (nodes [8237842, 8237850, 8237851]) [79%+, T]
  - cols 50-52 SITE: A1 code 50 = 406-414 paths (nodes [8237854, 8237874]) [100%+, --T] | REF code 50 = 406-414 paths (nodes [8237854, 8237874]) [100%+, --T]
  - cols 53-53: A1 code 51 = 415 paths (nodes [8237855]) [100%+, A] | REF code 51 = 415 paths (nodes [8237855]) [100%+, A]
  - cols 54-61: A1 code 80 = 444 paths (nodes [8237856]) [100%+, TATATATA] | REF code 80 = 444 paths (nodes [8237856]) [100%+, TATATATA]
  - cols 62-62: A1 code 69 = 433 paths (nodes [8237857]) [100%+, T] | REF code 69 = 433 paths (nodes [8237857]) [100%+, T]
  - cols 63-63: A1 code 68 = 432 paths (nodes [8237861]) [100%+, A] | REF code 68 = 432 paths (nodes [8237861]) [100%+, A]
  - cols 64-64: A1 code 55 = 419 paths (nodes [8237862]) [100%+, T] | REF code 55 = 419 paths (nodes [8237862]) [100%+, T]
  - cols 65-65 DIFF: A1 code 30 = 141-150 paths (nodes [8237873]) [100%+, A] | REF code 42 = 284-297 paths (nodes [8237871]) [85%+, T]
  - cols 66-66: A1 code 50 = 406-414 paths (nodes [8237854, 8237874]) [100%+, A] | REF code 50 = 406-414 paths (nodes [8237854, 8237874]) [94%+, A]
  - cols 67-67: A1 code 53 = 417 paths (nodes [8237876]) [100%+, T] | REF code 53 = 417 paths (nodes [8237876]) [88%+, T]
  - cols 68-68: A1 code 63 = 427 paths (nodes [8237878]) [100%+, A] | REF code 63 = 427 paths (nodes [8237878]) [88%+, A]
  - cols 69-74: A1 code 64 = 428 paths (nodes [8237879]) [100%+, TATATA] | REF code 64 = 428 paths (nodes [8237879]) [88%+, TATATA]
  - cols 75-75: A1 code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [100%+, T] | REF code 46 = 342-357 paths (nodes [8237822, 8237823, 8237880]) [88%+, T]
  - cols 76-76 DIFF: A1 code 79 = 443 paths (nodes [8237888]) [100%+, A] | REF code 25 = 96-104 paths (nodes [8237848, 8237881, 8237882]) [79%+, A]
  - cols 77-79 DIFF: A1 code 99 = 463 paths (nodes [8237928, 8237930, 8237934, 8237936]) [100%+, TAT] | REF code 25 = 96-104 paths (nodes [8237848, 8237881, 8237882]) [79%+, TAT]
  - cols 80-80 DIFF: A1 code 99 = 463 paths (nodes [8237928, 8237930, 8237934, 8237936]) [100%+, A] | REF code 26 = 105-112 paths (nodes [8237739, 8237829, 8237830, 8237841, 8237885, 8237890]) [79%+, T]
  - cols 81-86 DIFF: A1 code 99 = 463 paths (nodes [8237928, 8237930, 8237934, 8237936]) [100%+, AGCTCT] | REF code 27 = 113-121 paths (nodes [8237738, 8237887, 8237891]) [76%+, ATATAT]
  - cols 87-90 DIFF: A1 code 100 = 464 paths (nodes [8237935, 8237938, 8237941]) [100%+, GTCC] | REF code 27 = 113-121 paths (nodes [8237738, 8237887, 8237891]) [82%+, ATAT]
  - cols 91-91 DIFF: A1 code 100 = 464 paths (nodes [8237935, 8237938, 8237941]) [100%+, C] | REF code 79 = 443 paths (nodes [8237888]) [76%+, A]
  - cols 92-96 DIFF: A1 code 100 = 464 paths (nodes [8237935, 8237938, 8237941]) [100%+, ACTAA] | REF code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [73%+, ATATA]
  - cols 97-97 DIFF: A1 code 99 = 463 paths (nodes [8237928, 8237930, 8237934, 8237936]) [100%+, G] | REF code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [79%+, T]
  - cols 98-98 DIFF: A1 code 100 = 464 paths (nodes [8237935, 8237938, 8237941]) [100%+, A] | REF code 35 = 195-205 paths (nodes [8237817, 8237835, 8237837, 8237843]) [79%+, A]
  - cols 99-100 DIFF: A1 code 100 = 464 paths (nodes [8237935, 8237938, 8237941]) [100%+, TA] | REF code 36 = 206-217 paths (nodes [8237825, 8237845, 8237926]) [70%+, TA]

## R1.m2_INDEL_08  HG008T_PacBio_d46  8349226:45:INS:>CA  (shard chr1_shard_00003_data.npy index 6112)

- truth chr1:194115661:G>GCA; INS ->CA; AF 0.290; rows A1 27, REF 65, other 1; H = 464
- candidate node 8349226 (47 bp, 429 paths), site offset 45-45; allele nodes [8349226]; candidate columns [50, 52]
- GRCh38-path nodes under the window: 14: 8349223 (47 bp, 430), 8349224 (1 bp, 429), 8349226 (47 bp, 429), 8349232 (2 bp, 210), 8349233 (2 bp, 243), 8349234 (2 bp, 299), 8349235 (2 bp, 405), 8349236 (2 bp, 463), 8349237 (1 bp, 462), 8349239 (23 bp, 463), 8349240 (1 bp, 447), 8349241 (1 bp, 445), 8349242 (1 bp, 208), 8349244 (198 bp, 464)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 0-3: A1 code 66 = 430 paths (nodes [8349223]) [100%+, TGAC] | REF code 66 = 430 paths (nodes [8349223]) [100%+, TGAC]
  - cols 4-49: A1 code 65 = 429 paths (nodes [8349224, 8349226]) [100%+, ATGTACCCAACATAACCACAAGGTTCATGTCAGGAGAAAAAGCATG] | REF code 65 = 429 paths (nodes [8349224, 8349226]) [100%+, ATGTACCCAACATAACCACAAGGTTCATGTCAGGAGAAAAAGCATG]
  - cols 50-52 SITE: A1 code 65 = 429 paths (nodes [8349224, 8349226]) [100%+, --C] | REF code 65 = 429 paths (nodes [8349224, 8349226]) [100%+, --C]
  - cols 53-53: A1 code 65 = 429 paths (nodes [8349224, 8349226]) [100%+, A] | REF code 65 = 429 paths (nodes [8349224, 8349226]) [100%+, A]
  - cols 54-55 DIFF: A1 code 36 = 206-217 paths (nodes [8349232, 8349242]) [100%+, CA] | REF code 49 = 390-405 paths (nodes [8349235]) [86%+, CA]
  - cols 56-57 DIFF: A1 code 39 = 243-255 paths (nodes [8349233, 8349243]) [100%+, CA] | REF code 99 = 463 paths (nodes [8349236, 8349239]) [86%+, CA]
  - cols 58-58 DIFF: A1 code 43 = 298-311 paths (nodes [8349234]) [100%+, C] | REF code 98 = 462 paths (nodes [8349237, 8349245]) [86%+, C]
  - cols 59-59 DIFF: A1 code 43 = 298-311 paths (nodes [8349234]) [100%+, A] | REF code 99 = 463 paths (nodes [8349236, 8349239]) [92%+, A]
  - cols 60-61 DIFF: A1 code 49 = 390-405 paths (nodes [8349235]) [100%+, CA] | REF code 99 = 463 paths (nodes [8349236, 8349239]) [88%+, CA]
  - cols 62-63: A1 code 99 = 463 paths (nodes [8349236, 8349239]) [100%+, CA] | REF code 99 = 463 paths (nodes [8349236, 8349239]) [100%+, CA]
  - cols 64-64 DIFF: A1 code 98 = 462 paths (nodes [8349237, 8349245]) [100%+, C] | REF code 99 = 463 paths (nodes [8349236, 8349239]) [92%+, C]
  - cols 65-81: A1 code 99 = 463 paths (nodes [8349236, 8349239]) [100%+, ACACACACACACACACA] | REF code 99 = 463 paths (nodes [8349236, 8349239]) [98%+, ACACACACACACACACA]
  - cols 82-82 DIFF: A1 code 99 = 463 paths (nodes [8349236, 8349239]) [100%+, C] | REF code 83 = 447 paths (nodes [8349240]) [86%+, C]
  - cols 83-83 DIFF: A1 code 99 = 463 paths (nodes [8349236, 8349239]) [100%+, A] | REF code 81 = 445 paths (nodes [8349241]) [86%+, A]
  - cols 84-84 DIFF: A1 code 99 = 463 paths (nodes [8349236, 8349239]) [100%+, C] | REF code 39 = 243-255 paths (nodes [8349233, 8349243]) [86%+, T]
  - cols 85-87 DIFF: A1 code 99 = 463 paths (nodes [8349236, 8349239]) [100%+, ACA] | REF code 100 = 464 paths (nodes [8349244, 8349247]) [88%+, ATA]
  - cols 88-88 DIFF: A1 code 83 = 447 paths (nodes [8349240]) [100%+, C] | REF code 100 = 464 paths (nodes [8349244, 8349247]) [92%+, C]
  - cols 89-89 DIFF: A1 code 81 = 445 paths (nodes [8349241]) [100%+, A] | REF code 100 = 464 paths (nodes [8349244, 8349247]) [92%+, A]
  - cols 90-90 DIFF: A1 code 36 = 206-217 paths (nodes [8349232, 8349242]) [100%+, C] | REF code 100 = 464 paths (nodes [8349244, 8349247]) [92%+, C]
  - cols 91-93: A1 code 100 = 464 paths (nodes [8349244, 8349247]) [100%+, ATA] | REF code 100 = 464 paths (nodes [8349244, 8349247]) [97%+, AAT]

## R1.m2_INDEL_09  HG008T_PacBio_d46  11316111:0:DEL:A>  (shard chr1_shard_00004_data.npy index 15274)

- truth chr1:240187269:CA>C; DEL A>-; AF 0.208; rows A1 11, REF 42, other 0; H = 464
- candidate node 11316111 (15 bp, 462 paths), site offset 0-1; allele nodes [11316111]; candidate columns [50, 51]
- GRCh38-path nodes under the window: 9: 11316103 (44 bp, 459), 11316104 (1 bp, 448), 11316106 (42 bp, 458), 11316111 (15 bp, 462), 11316112 (1 bp, 461), 11316114 (1 bp, 462), 11316115 (1 bp, 231), 11316117 (5 bp, 230), 11316118 (38 bp, 464)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 0-3: A1 code 95 = 459 paths (nodes [11316100, 11316103]) [100%+, AGAT] | REF code 95 = 459 paths (nodes [11316100, 11316103]) [100%+, TCGC]
  - cols 4-4 DIFF: A1 code 95 = 459 paths (nodes [11316100, 11316103]) [100%+, C] | REF code 84 = 448 paths (nodes [11316104]) [90%+, G]
  - cols 5-6 DIFF: A1 code 95 = 459 paths (nodes [11316100, 11316103]) [100%+, GC] | REF code 94 = 458 paths (nodes [11316106]) [90%+, CC]
  - cols 7-7 DIFF: A1 code 84 = 448 paths (nodes [11316104]) [100%+, G] | REF code 94 = 458 paths (nodes [11316106]) [95%+, G]
  - cols 8-46: A1 code 94 = 458 paths (nodes [11316106]) [100%+, CCGTTGCACTCCAGCCTGGGTGAAAGAGCGAAACTCCAT] | REF code 94 = 458 paths (nodes [11316106]) [100%+, TTGCACTCCAGCCTGGGTGAAAGAGCGAAACTCCATCTC]
  - cols 47-47 DIFF: A1 code 94 = 458 paths (nodes [11316106]) [100%+, C] | REF code 27 = 113-121 paths (nodes [11316108]) [90%+, A]
  - cols 48-48 DIFF: A1 code 94 = 458 paths (nodes [11316106]) [100%+, T] | REF code 28 = 122-131 paths (nodes [11316109]) [90%+, A]
  - cols 49-49 DIFF: A1 code 94 = 458 paths (nodes [11316106]) [100%+, C] | REF code 39 = 243-255 paths (nodes [11316110]) [95%+, A]
  - cols 50-51 SITE: A1 code 98 = 462 paths (nodes [11316111, 11316114]) [100%+, AA] | REF code 98 = 462 paths (nodes [11316111, 11316114]) [100%+, AA]
  - cols 52-64: A1 code 98 = 462 paths (nodes [11316111, 11316114]) [100%+, AAAAAAAAAAAAA] | REF code 98 = 462 paths (nodes [11316111, 11316114]) [100%+, AAAAAAAAAAAAA]
  - cols 65-65: A1 code 97 = 461 paths (nodes [11316112]) [100%+, A] | REF code 97 = 461 paths (nodes [11316112]) [100%+, A]
  - cols 66-66: A1 code 98 = 462 paths (nodes [11316111, 11316114]) [100%+, A] | REF code 98 = 462 paths (nodes [11316111, 11316114]) [100%+, A]
  - cols 67-67: A1 code 38 = 231-242 paths (nodes [11316115, 11316116]) [100%+, A] | REF code 38 = 231-242 paths (nodes [11316115, 11316116]) [100%+, G]
  - cols 68-72 DIFF: A1 code 37 = 218-230 paths (nodes [11316117]) [100%+, AAAAG] | REF code 100 = 464 paths (nodes [11316118, 11316121]) [95%+, AAAAG]

## R1.m2_INDEL_10  COLO829T_fiberseq  266437:1:INS:>TTGTTGTTGTTGGTTTTTC  (shard chr1_shard_00000_data.npy index 8215)

- truth chr1:9364848:G>GTTGTTGTTGTTGGTTTTTC; INS ->TTGTTGTTGTTGGTTTTTC; AF 0.357; rows A1 30, REF 52, other 2; H = 90
- candidate node 266437 (2 bp, 90 paths), site offset 1-1; allele nodes [266437]; candidate columns [50, 69]
- GRCh38-path nodes under the window: 7: 266432 (63 bp, 90), 266433 (3 bp, 74), 266434 (24 bp, 90), 266435 (1 bp, 56), 266437 (2 bp, 90), 266440 (1 bp, 89), 266441 (99 bp, 90)
- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:
  - cols 13-15: A1 code 100 = 90 paths (nodes [266430, 266432, 266434, 266437, 266441, 266443]) [100%+, GGG] | REF code 100 = 90 paths (nodes [266430, 266432, 266434, 266437, 266441, 266443]) [98%+, GGG]
  - cols 20-20: A1 code 100 = 90 paths (nodes [266430, 266432, 266434, 266437, 266441, 266443]) [100%+, T] | REF code 100 = 90 paths (nodes [266430, 266432, 266434, 266437, 266441, 266443]) [98%+, T]
  - cols 21-23: A1 code 84 = 74 paths (nodes [266433]) [100%+, TTG] | REF code 84 = 74 paths (nodes [266433]) [96%+, TTG]
  - cols 48-48: A1 code 66 = 56 paths (nodes [266435]) [100%+, C] | REF code 66 = 56 paths (nodes [266435]) [98%+, C]
  - cols 70-70 DIFF: A1 code 99 = 89 paths (nodes [266440]) [100%+, T] | REF code 34 = 34 paths (nodes [266436, 266438, 266444]) [94%+, C]
  - cols 71-88 DIFF: A1 code 100 = 90 paths (nodes [266430, 266432, 266434, 266437, 266441, 266443]) [100%+, GTTGTTGTTGGTTTTTCT] | REF code 33 = 33 paths (nodes [266439]) [92%+, G----------------T]
  - cols 89-89 DIFF: A1 code 100 = 90 paths (nodes [266430, 266432, 266434, 266437, 266441, 266443]) [100%+, G] | REF code 99 = 89 paths (nodes [266440]) [92%+, T]

