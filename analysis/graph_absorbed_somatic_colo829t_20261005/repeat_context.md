# COLO829T absorbed truths: repeat context (c4)

Tables of `c4_tables.py` (definitions in its docstring and in README.md). Set = the truths that are perfect bypass on that platform; union = on >= 1 platform. Patient frame = COLO829BL (matched normal) verkko 2.1 hapX / hapY; GRCh38 frame = the truth VCF record. chr1-22 only (per_truth_COLO829T.tsv holds the chr1-22 alleles).

### 1. Perfect bypass per platform (chr1-22; the chr1 column is a subset)

| platform | SNV no candidate | SNV perfect bypass | of which a d9 path spells the truth | INDEL no candidate | INDEL perfect bypass | of which a d9 path spells the truth | chr1 SNV / INDEL perfect bypass |
|---|---:|---:|---:|---:|---:|---:|---:|
| fiberseq | 528 | 64 (12.1%) | 61 | 436 | 218 (50.0%) | 205 | 5 / 23 |
| ONT | 547 | 58 (10.6%) | 56 | 347 | 170 (49.0%) | 154 | 4 / 19 |
| Illumina | 598 | 53 (8.9%) | 51 | 687 | 396 (57.6%) | 380 | 4 / 34 |
| union |  | 67 | 64 |  | 430 | 409 | 5 / 38 |

All chr1-22 truth alleles: SNV 42,035, INDEL 1,912. Perfect bypass on all three platforms: SNV 46, INDEL 134. Share of all truth: fiberseq SNV 0.15% / INDEL 11.4%; ONT SNV 0.14% / INDEL 8.9%; Illumina SNV 0.13% / INDEL 20.7%; union SNV 0.16% / INDEL 22.5%

### 1b. Miss classes and d9 match of the perfect-bypass truths

| platform | SNV S3 | SNV S4b | SNV S4a | SNV S2 | INDEL I1 | SNV exact | SNV with_germline | SNV closest | INDEL exact | INDEL with_germline | INDEL closest |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fiberseq | 57 | 6 | 0 | 1 | 218 | 61 | 0 | 3 | 195 | 10 | 13 |
| ONT | 51 | 2 | 4 | 1 | 170 | 56 | 0 | 2 | 144 | 10 | 16 |
| Illumina | 50 | 2 | 1 | 0 | 396 | 51 | 0 | 2 | 367 | 13 | 16 |
| union |  |  |  |  |  | 64 | 0 | 3 | 394 | 15 | 21 |

S3 = the ALT is an existing graph SNP allele; S4b = repeat, absorbed by a graph path; S4a = repeat, the SNV becomes an edit on another branch; S2 = ALT not seen in same-length reads; I1 = the INDEL allele is fully in the graph. exact / with_germline / closest: c2 (a d9 path spells GRCh38 + truth; only with COLO829BL dipcall alleles applied too; no such path).

### 2a. INDELs, patient frame (COLO829BL): class label per perfect set

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| HP>=7 | 79 (36.2%) | 51 (30.0%) | 235 (59.3%) | 241 (56.0%) |
| STR2-6>=3copies | 57 (26.1%) | 44 (25.9%) | 75 (18.9%) | 83 (19.3%) |
| VNTR>6 | 1 (0.5%) | 1 (0.6%) | 1 (0.3%) | 1 (0.2%) |
| in_STR | 12 (5.5%) | 13 (7.6%) | 15 (3.8%) | 17 (4.0%) |
| imperfect VNTR | 2 (0.9%) | 1 (0.6%) | 1 (0.3%) | 2 (0.5%) |
| SNV in STR2-6>=3copies | 10 (4.6%) | 11 (6.5%) | 11 (2.8%) | 15 (3.5%) |
| SNV in HP4-6 | 1 (0.5%) | 1 (0.6%) | 1 (0.3%) | 1 (0.2%) |
| SNV | 6 (2.8%) | 5 (2.9%) | 5 (1.3%) | 6 (1.4%) |
| none | 2 (0.9%) | 3 (1.8%) | 3 (0.8%) | 3 (0.7%) |
| complex | 21 (9.6%) | 18 (10.6%) | 21 (5.3%) | 27 (6.3%) |
| normal carries ALT | 27 (12.4%) | 22 (12.9%) | 27 (6.8%) | 33 (7.7%) |
| NA | 0 (0.0%) | 0 (0.0%) | 1 (0.3%) | 1 (0.2%) |
| total | 218 | 170 | 396 | 430 |

### 2b. INDELs, patient frame: grouped

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| HP>=7 | 79 (36.2%) | 51 (30.0%) | 235 (59.3%) | 241 (56.0%) |
| STR | 57 (26.1%) | 44 (25.9%) | 75 (18.9%) | 83 (19.3%) |
| other repeat | 15 (6.9%) | 15 (8.8%) | 17 (4.3%) | 20 (4.7%) |
| SNV in repeat | 11 (5.0%) | 12 (7.1%) | 12 (3.0%) | 16 (3.7%) |
| non-repeat | 8 (3.7%) | 8 (4.7%) | 8 (2.0%) | 9 (2.1%) |
| normal carries ALT | 27 (12.4%) | 22 (12.9%) | 27 (6.8%) | 33 (7.7%) |
| complex | 21 (9.6%) | 18 (10.6%) | 21 (5.3%) | 27 (6.3%) |
| unresolved | 0 (0.0%) | 0 (0.0%) | 1 (0.3%) | 1 (0.2%) |
| total | 218 | 170 | 396 | 430 |

### 2c. INDELs, patient frame: units changed (repeat INDEL events only)

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| 1 unit | 130 (59.6%) | 91 (53.5%) | 305 (77.0%) | 316 (73.5%) |
| 2-3 units | 5 (2.3%) | 3 (1.8%) | 5 (1.3%) | 7 (1.6%) |
| >3 units | 2 (0.9%) | 2 (1.2%) | 1 (0.3%) | 2 (0.5%) |
| (not a repeat INDEL event) | 81 (37.2%) | 74 (43.5%) | 85 (21.5%) | 105 (24.4%) |
| total | 218 | 170 | 396 | 430 |

### 2d. INDELs, patient frame: the normal's homopolymer length (HP>=7 events only)

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| 7-9 | 3 (1.4%) | 4 (2.4%) | 6 (1.5%) | 6 (1.4%) |
| 10-14 | 29 (13.3%) | 18 (10.6%) | 83 (21.0%) | 86 (20.0%) |
| 15-19 | 32 (14.7%) | 18 (10.6%) | 93 (23.5%) | 94 (21.9%) |
| 20-29 | 15 (6.9%) | 11 (6.5%) | 53 (13.4%) | 55 (12.8%) |
| (not HP>=7) | 139 (63.8%) | 119 (70.0%) | 161 (40.7%) | 189 (44.0%) |
| total | 218 | 170 | 396 | 430 |

### 2e. INDELs, GRCh38 frame (b1 rule): class per perfect set

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| HP>=7 | 81 (37.2%) | 47 (27.6%) | 238 (60.1%) | 244 (56.7%) |
| STR2-6>=3copies | 86 (39.4%) | 71 (41.8%) | 99 (25.0%) | 120 (27.9%) |
| HP4-6 | 2 (0.9%) | 3 (1.8%) | 3 (0.8%) | 3 (0.7%) |
| VNTR>6 | 2 (0.9%) | 2 (1.2%) | 2 (0.5%) | 2 (0.5%) |
| none | 47 (21.6%) | 47 (27.6%) | 54 (13.6%) | 61 (14.2%) |
| total | 218 | 170 | 396 | 430 |

### 2f. INDELs, patient frame with the b5 rule (HG008 rule): class label per perfect set

| class | fiberseq (INDEL) | ONT (INDEL) | Illumina (INDEL) | union (INDEL) |
|---|---:|---:|---:|---:|
| HP>=7 | 86 (39.4%) | 57 (33.5%) | 244 (61.6%) | 250 (58.1%) |
| STR2-6>=3copies | 62 (28.4%) | 49 (28.8%) | 81 (20.5%) | 89 (20.7%) |
| VNTR>6 | 1 (0.5%) | 1 (0.6%) | 1 (0.3%) | 1 (0.2%) |
| in_STR | 1 (0.5%) | 2 (1.2%) | 1 (0.3%) | 2 (0.5%) |
| SNV in STR2-6>=3copies | 10 (4.6%) | 11 (6.5%) | 11 (2.8%) | 15 (3.5%) |
| SNV in HP4-6 | 1 (0.5%) | 1 (0.6%) | 1 (0.3%) | 1 (0.2%) |
| SNV | 6 (2.8%) | 5 (2.9%) | 5 (1.3%) | 6 (1.4%) |
| none | 3 (1.4%) | 4 (2.4%) | 3 (0.8%) | 5 (1.2%) |
| complex | 21 (9.6%) | 18 (10.6%) | 21 (5.3%) | 27 (6.3%) |
| normal carries ALT | 27 (12.4%) | 22 (12.9%) | 27 (6.8%) | 33 (7.7%) |
| NA | 0 (0.0%) | 0 (0.0%) | 1 (0.3%) | 1 (0.2%) |
| total | 218 | 170 | 396 | 430 |

### 2g. Union INDELs: GRCh38-frame class x patient-frame group

| GRCh38 class \ patient frame | HP>=7 | STR | other repeat | SNV in repeat | non-repeat | normal carries ALT | complex | unresolved | total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HP>=7 | 227 | 0 | 5 | 0 | 0 | 9 | 3 | 0 | 244 |
| STR2-6>=3copies | 8 | 76 | 8 | 7 | 3 | 6 | 11 | 1 | 120 |
| HP4-6 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 3 |
| VNTR>6 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 2 |
| none | 6 | 7 | 6 | 9 | 6 | 15 | 12 | 0 | 61 |

GRCh38-frame "none" INDELs (union) in the patient frame: normal carries ALT 15, complex 12, SNV in STR2-6>=3copies 8, STR2-6>=3copies 7, HP>=7 6, SNV 4, in_STR 4, none 2, imperfect VNTR 2, SNV in HP4-6 1. The same INDELs under the b5 rule in GRCh38 (c1 b5rule_class): in_STR 32, none 29.

### 2h. Union INDELs: miscellany

| item | counts |
|---|---:|
| normal carries ALT: other hap's event | none 13, HP>=7 8, in_STR 3, HP4-6 3, complex 3, STR2-6>=3copies 2, SNV in HP>=7 1 |
| normal carries ALT: hap carrying it | hapX 30, hapY 3 |
| complex: other hap's event | complex 20, in_STR 6, VNTR>6 1 |
| patient-frame INDEL events vs the GRCh38 length change | same 278, other size 52, opposite sign 17 |
| array allele of the event hap | REF 270, germline_len 115, ALT 33, germline_seq 11 |
| HP>=7 events: ins / del | DEL 129, INS 112 |
| HP>=7 events: base | A/T 210, C/G 31 |
| STR events: unit length | 2 62, 3 2, 4 18, 5 1 |
| STR events: motif (canonical, both strands) | AC 34, AT 23, AGAT 13, AG 5, AAAG 3, AAAT 2, AAT 2, AAATG 1 |
| normal carries ALT: SMaHT region | Extreme 32, Difficult 1 |
| normal carries ALT: VAF_Ill | 0.1-0.25 1, 0.25-0.4 10, >=0.4 22 |
| normal carries ALT: assembly and dipcall both ALT on a hap / any c3 flag | 23 / 10 of 33 |
| normal carries ALT: loci per chromosome Mb (>= 2) | chr1:146Mb 3, chr1:148Mb 3 |
| complex: SMaHT region | Extreme 24, Difficult 2, Easy 1 |
| a normal hap has the ALT length but not the ALT (alt_len_haps), by label | SNV in STR2-6>=3copies 15, complex 9, SNV 6, HP>=7 2, none 1, STR2-6>=3copies 1, SNV in HP4-6 1 |
| in an imperfect VNTR / minisatellite (in_ivntr), by label | complex 4, imperfect VNTR 2, normal carries ALT 2, STR2-6>=3copies 1 |
| a hap taken from dipcall's own copy (paralog fix), by label | normal carries ALT 7, SNV 2, none 1, SNV in HP4-6 1 |
| wide anchors (2 kb retry), by label | imperfect VNTR 1, HP>=7 1, normal carries ALT 1 |
| single hap: the one hap's label | (none absorbed) |
| label changed by the audit fixes, all absorbed truths (truth kind: before -> after) | 2561 INDEL: none -> SNV; 2562 INDEL: STR2-6>=3copies -> SNV; 2563 INDEL: HP4-6 -> normal carries ALT; 2570 SNV: SNV -> SNV in imperfect VNTR; 6942 SNV: SNV -> SNV in imperfect VNTR; 6943 SNV: SNV -> SNV in imperfect VNTR; 29401 INDEL: none -> imperfect VNTR; 35166 INDEL: none -> normal carries ALT; 39443 INDEL: none -> imperfect VNTR; 40329 INDEL: none -> normal carries ALT |
| a hap with >= 4 bp of reverted flank length (flank_len) / the event hap | 14 / 8 |

### 3a. All 1,912 chr1-22 truth INDELs, patient frame: share absorbed (union)

| class | truth INDELs | absorbed | % absorbed |
|---|---:|---:|---:|
| HP>=7 | 1,103 | 241 | 21.8% |
| STR2-6>=3copies | 189 | 83 | 43.9% |
| HP4-6 | 67 | 0 | 0.0% |
| VNTR>6 | 13 | 1 | 7.7% |
| in_STR | 92 | 17 | 18.5% |
| imperfect VNTR | 4 | 2 | 50.0% |
| SNV in HP>=7 | 16 | 0 | 0.0% |
| SNV in STR2-6>=3copies | 46 | 15 | 32.6% |
| SNV in HP4-6 | 3 | 1 | 33.3% |
| SNV | 11 | 6 | 54.5% |
| none | 235 | 3 | 1.3% |
| complex | 72 | 27 | 37.5% |
| normal carries ALT | 55 | 33 | 60.0% |
| single hap | 3 | 0 | 0.0% |
| NA | 3 | 1 | 33.3% |
| total | 1,912 | 430 | 22.5% |

### 3b. All chr1-22 truth INDELs, patient-frame mechanism class: share absorbed (union)

| class | truth INDELs | absorbed | % |
|---|---:|---:|---:|
| HP>=7 1-unit, tract 7-9 | 116 | 6 | 5.2% |
| HP>=7 1-unit, tract 10-14 | 424 | 86 | 20.3% |
| HP>=7 1-unit, tract 15-19 | 375 | 94 | 25.1% |
| HP>=7 1-unit, tract 20-29 | 163 | 53 | 32.5% |
| HP>=7 1-unit, tract >=30 | 1 | 0 | 0.0% |
| STR 1-unit | 165 | 76 | 46.1% |
| repeat multi-unit | 48 | 9 | 18.8% |
| HP4-6 | 67 | 0 | 0.0% |
| VNTR>6 | 13 | 1 | 7.7% |
| in_STR | 92 | 17 | 18.5% |
| imperfect VNTR | 4 | 2 | 50.0% |
| SNV in repeat | 65 | 16 | 24.6% |
| non-repeat INDEL | 235 | 3 | 1.3% |
| non-repeat SNV | 11 | 6 | 54.5% |
| complex | 72 | 27 | 37.5% |
| normal carries ALT | 55 | 33 | 60.0% |
| unresolved | 6 | 1 | 16.7% |

### 3c. All chr1-22 truth INDELs, GRCh38 frame (b1 rule): share absorbed (union)

| class | truth INDELs | absorbed | % |
|---|---:|---:|---:|
| HP>=7 | 1,167 | 244 | 20.91% |
| STR2-6>=3copies | 275 | 120 | 43.64% |
| HP4-6 | 97 | 3 | 3.09% |
| VNTR>6 | 4 | 2 | 50.00% |
| none | 369 | 61 | 16.53% |
| total | 1,912 | 430 | 22.49% |

### 3d. All chr1-22 truth SNVs, GRCh38 frame (b1 rule): share absorbed (union)

| class | truth SNVs | absorbed | % |
|---|---:|---:|---:|
| HP>=7 | 545 | 4 | 0.73% |
| STR2-6>=3copies | 1,322 | 14 | 1.06% |
| HP4-6 | 5,914 | 7 | 0.12% |
| none | 34,254 | 42 | 0.12% |
| total | 42,035 | 67 | 0.16% |

### 3e. Absorbed SNVs, GRCh38-frame class per perfect set

| class | fiberseq (SNV) | ONT (SNV) | Illumina (SNV) | union (SNV) |
|---|---:|---:|---:|---:|
| HP>=7 | 4 (6.2%) | 3 (5.2%) | 2 (3.8%) | 4 (6.0%) |
| STR2-6>=3copies | 13 (20.3%) | 9 (15.5%) | 11 (20.8%) | 14 (20.9%) |
| HP4-6 | 7 (10.9%) | 6 (10.3%) | 7 (13.2%) | 7 (10.4%) |
| none | 40 (62.5%) | 40 (69.0%) | 33 (62.3%) | 42 (62.7%) |
| total | 64 | 58 | 53 | 67 |

### 3f. Absorbed SNVs, patient frame (COLO829BL) per perfect set

| class | fiberseq (SNV) | ONT (SNV) | Illumina (SNV) | union (SNV) |
|---|---:|---:|---:|---:|
| in_STR | 2 (3.1%) | 0 (0.0%) | 1 (1.9%) | 2 (3.0%) |
| SNV in HP>=7 | 3 (4.7%) | 3 (5.2%) | 2 (3.8%) | 3 (4.5%) |
| SNV in STR2-6>=3copies | 9 (14.1%) | 7 (12.1%) | 8 (15.1%) | 10 (14.9%) |
| SNV in HP4-6 | 7 (10.9%) | 6 (10.3%) | 7 (13.2%) | 7 (10.4%) |
| SNV in imperfect VNTR | 3 (4.7%) | 3 (5.2%) | 0 (0.0%) | 3 (4.5%) |
| SNV | 36 (56.2%) | 36 (62.1%) | 32 (60.4%) | 38 (56.7%) |
| complex | 3 (4.7%) | 2 (3.4%) | 2 (3.8%) | 3 (4.5%) |
| normal carries ALT | 1 (1.6%) | 1 (1.7%) | 1 (1.9%) | 1 (1.5%) |
| total | 64 | 58 | 53 | 67 |

### 3g. Absorbed SNVs, CpG transition (GRCh38) per perfect set

| class | fiberseq (SNV) | ONT (SNV) | Illumina (SNV) | union (SNV) |
|---|---:|---:|---:|---:|
| CpG Ti | 26 (40.6%) | 24 (41.4%) | 23 (43.4%) | 28 (41.8%) |
| other | 38 (59.4%) | 34 (58.6%) | 30 (56.6%) | 39 (58.2%) |
| total | 64 | 58 | 53 | 67 |

### 3h. All chr1-22 truth SNVs: share absorbed by CpG transition and substitution

| class | truth SNVs | absorbed | % |
|---|---:|---:|---:|
| CpG Ti | 2,705 | 28 | 1.04% |
| not CpG Ti | 39,330 | 39 | 0.10% |
| C>A | 5,178 | 6 | 0.12% |
| C>G | 1,338 | 0 | 0.00% |
| C>T | 28,417 | 40 | 0.14% |
| T>A | 2,136 | 9 | 0.42% |
| T>C | 2,762 | 7 | 0.25% |
| T>G | 2,204 | 5 | 0.23% |

### 3i. Absorbed SNVs (union): GRCh38 class x CpG Ti x patient frame

| GRCh38 class | CpG Ti | patient frame | SNVs |
|---|---:|---:|---:|
| none | CpG Ti | SNV | 22 |
| none | other | SNV | 16 |
| STR2-6>=3copies | other | SNV in STR2-6>=3copies | 9 |
| HP4-6 | CpG Ti | SNV in HP4-6 | 4 |
| HP>=7 | other | SNV in HP>=7 | 3 |
| none | other | SNV in imperfect VNTR | 3 |
| HP4-6 | other | SNV in HP4-6 | 3 |
| STR2-6>=3copies | other | complex | 3 |
| STR2-6>=3copies | CpG Ti | in_STR | 1 |
| HP>=7 | other | in_STR | 1 |
| STR2-6>=3copies | CpG Ti | SNV in STR2-6>=3copies | 1 |
| none | other | normal carries ALT | 1 |

### 4a. SMaHT region x patient-frame group, INDELs

| patient-frame group \ RGN | Easy | Difficult | Extreme |
|---|---:|---:|---:|
| HP>=7 | 24 / 247 (9.7%) | 135 / 498 (27.1%) | 82 / 358 (22.9%) |
| STR | 5 / 23 (21.7%) | 35 / 68 (51.5%) | 43 / 98 (43.9%) |
| other repeat | 1 / 78 (1.3%) | 8 / 29 (27.6%) | 11 / 69 (15.9%) |
| SNV in repeat | 0 / 11 (0.0%) | 3 / 11 (27.3%) | 13 / 43 (30.2%) |
| non-repeat | 2 / 198 (1.0%) | 1 / 17 (5.9%) | 6 / 31 (19.4%) |
| normal carries ALT | 0 / 0 | 1 / 6 (16.7%) | 32 / 49 (65.3%) |
| complex | 1 / 5 (20.0%) | 2 / 11 (18.2%) | 24 / 56 (42.9%) |
| unresolved | 0 / 0 | 1 / 1 (100.0%) | 0 / 5 (0.0%) |
| all INDELs | 33 / 562 (5.9%) | 186 / 641 (29.0%) | 211 / 709 (29.8%) |

Cells: absorbed (union) / all chr1-22 truth INDELs of the class and stratum (% absorbed).

### 4b. Truth VAF (Illumina) x patient-frame group, INDELs

| patient-frame group \ VAF | <0.1 | 0.1-0.25 | 0.25-0.4 | >=0.4 |
|---|---:|---:|---:|---:|
| HP>=7 | 0 / 0 | 19 / 105 (18.1%) | 84 / 355 (23.7%) | 138 / 643 (21.5%) |
| STR | 0 / 0 | 14 / 25 (56.0%) | 33 / 63 (52.4%) | 36 / 101 (35.6%) |
| other repeat | 0 / 0 | 2 / 18 (11.1%) | 1 / 44 (2.3%) | 17 / 114 (14.9%) |
| SNV in repeat | 0 / 0 | 1 / 3 (33.3%) | 4 / 17 (23.5%) | 11 / 45 (24.4%) |
| non-repeat | 0 / 0 | 1 / 15 (6.7%) | 2 / 63 (3.2%) | 6 / 168 (3.6%) |
| normal carries ALT | 0 / 0 | 1 / 5 (20.0%) | 10 / 18 (55.6%) | 22 / 32 (68.8%) |
| complex | 0 / 0 | 3 / 11 (27.3%) | 14 / 27 (51.9%) | 10 / 34 (29.4%) |
| unresolved | 0 / 0 | 0 / 0 | 1 / 3 (33.3%) | 0 / 3 (0.0%) |
| all INDELs | 0 / 0 | 41 / 182 (22.5%) | 149 / 590 (25.3%) | 240 / 1,140 (21.1%) |

Cells: absorbed (union) / all chr1-22 truth INDELs of the class and stratum (% absorbed).

### 4c. Normal-frame germline status x patient-frame group, INDELs

| patient-frame group \ germline status | normal carries ALT | germline other length | germline same length, other bases | both REF | one hap NA, other REF | both NA |
|---|---:|---:|---:|---:|---:|---:|
| HP>=7 | 0 / 0 | 56 / 232 (24.1%) | 6 / 30 (20.0%) | 179 / 841 (21.3%) | 0 / 0 | 0 / 0 |
| STR | 0 / 0 | 47 / 90 (52.2%) | 1 / 7 (14.3%) | 35 / 92 (38.0%) | 0 / 0 | 0 / 0 |
| other repeat | 0 / 0 | 12 / 51 (23.5%) | 7 / 22 (31.8%) | 1 / 103 (1.0%) | 0 / 0 | 0 / 0 |
| SNV in repeat | 0 / 0 | 16 / 65 (24.6%) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| non-repeat | 0 / 0 | 8 / 25 (32.0%) | 1 / 6 (16.7%) | 0 / 215 (0.0%) | 0 / 0 | 0 / 0 |
| normal carries ALT | 33 / 55 (60.0%) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| complex | 0 / 0 | 27 / 66 (40.9%) | 0 / 6 (0.0%) | 0 / 0 | 0 / 0 | 0 / 0 |
| unresolved | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 3 (0.0%) | 1 / 3 (33.3%) |
| all INDELs | 33 / 55 (60.0%) | 166 / 529 (31.4%) | 15 / 71 (21.1%) | 215 / 1,251 (17.2%) | 0 / 3 (0.0%) | 1 / 3 (33.3%) |

Cells: absorbed (union) / all chr1-22 truth INDELs of the class and stratum (% absorbed).

### 4d. HP>=7 1-unit events by tract length x SMaHT region, INDELs

| tract \ RGN | Easy | Difficult | Extreme |
|---|---:|---:|---:|
| tract 7-9 | 3 / 80 (3.8%) | 0 / 16 (0.0%) | 3 / 20 (15.0%) |
| tract 10-14 | 20 / 154 (13.0%) | 45 / 166 (27.1%) | 21 / 104 (20.2%) |
| tract 15-19 | 1 / 11 (9.1%) | 64 / 252 (25.4%) | 29 / 112 (25.9%) |
| tract 20-29 | 0 / 0 | 24 / 58 (41.4%) | 29 / 105 (27.6%) |
| tract >=30 | 0 / 0 | 0 / 0 | 0 / 1 (0.0%) |

Cells: absorbed (union) / all chr1-22 truth INDELs whose patient-frame event is a 1-base change of a homopolymer of that length.

SNVs by RGN: Easy 33 / 32,945 (0.1%); Difficult 11 / 4,922 (0.2%); Extreme 23 / 4,168 (0.6%)
SNVs by VAF: <0.1 6 / 4,439 (0.1%); 0.1-0.25 4 / 2,750 (0.1%); 0.25-0.4 13 / 6,737 (0.2%); >=0.4 44 / 28,109 (0.2%)

### 5a. Perfect bypass, COLO829T vs HG008T (chr1-22)

|  | COLO829T | HG008T |
|---|---:|---:|
| truth alleles SNV / INDEL | 42,035 / 1,912 | 8,690 / 8,496 |
| fiberseq / PacBio HiFi: SNV / INDEL perfect (% of no candidate) | 64 (12.1%) / 218 (50.0%) | 170 (42.7%) / 1,363 (45.8%) |
| ONT: SNV / INDEL perfect | 58 (10.6%) / 170 (49.0%) | 148 (43.3%) / 989 (40.3%) |
| Illumina: SNV / INDEL perfect | 53 (8.9%) / 396 (57.6%) | 122 (36.0%) / 1,916 (48.5%) |
| union SNV (% of all truth SNVs) | 67 (0.2%) | 201 (2.3%) |
| union INDEL (% of all truth INDELs) | 430 (22.5%) | 2,277 (26.8%) |

### 5b. Union absorbed INDELs, patient-frame group: COLO829T vs HG008T

| group | COLO829T (c3 rule) | COLO829T (b5 rule) | HG008T (b5 rule) |
|---|---:|---:|---:|
| HP>=7 | 241 (56.0%) | 250 (58.1%) | 1,853 (81.4%) |
| STR | 83 (19.3%) | 89 (20.7%) | 292 (12.8%) |
| other repeat | 20 (4.7%) | 3 (0.7%) | 6 (0.3%) |
| SNV in repeat | 16 (3.7%) | 16 (3.7%) | 104 (4.6%) |
| non-repeat | 9 (2.1%) | 11 (2.6%) | 1 (0.0%) |
| normal carries ALT | 33 (7.7%) | 33 (7.7%) | 0 (0.0%) |
| complex | 27 (6.3%) | 27 (6.3%) | 0 (0.0%) |
| unresolved | 1 (0.2%) | 1 (0.2%) | 21 (0.9%) |
| total | 430 | 430 | 2,277 |

HG008T patient frame = the GIAB truth INFO event (HG008Nv62SOMATICVARIANT) in the HG008-N v6.2 event haplotype, always a simple event, so HG008T has no 'complex' or 'normal carries ALT' rows; 'unresolved' = no INFO event there. HG008T truths whose tumor allele is the normal's other-haplotype allele (category HG008N_present_*): 105 INDELs (4.6%), classed by their INFO event inside the rows above.

### 5c. Union absorbed INDELs, patient-frame class label: COLO829T vs HG008T

| class | COLO829T (c3 rule) | COLO829T (b5 rule) | HG008T (b5 rule) |
|---|---:|---:|---:|
| HP>=7 | 241 (56.0%) | 250 (58.1%) | 1,853 (81.4%) |
| STR2-6>=3copies | 83 (19.3%) | 89 (20.7%) | 292 (12.8%) |
| HP4-6 | 0 (0.0%) | 0 (0.0%) | 1 (0.0%) |
| VNTR>6 | 1 (0.2%) | 1 (0.2%) | 3 (0.1%) |
| in_STR | 17 (4.0%) | 2 (0.5%) | 2 (0.1%) |
| imperfect VNTR | 2 (0.5%) | 0 (0.0%) | 0 (0.0%) |
| SNV in HP>=7 | 0 (0.0%) | 0 (0.0%) | 14 (0.6%) |
| SNV in STR2-6>=3copies | 15 (3.5%) | 15 (3.5%) | 90 (4.0%) |
| SNV in HP4-6 | 1 (0.2%) | 1 (0.2%) | 0 (0.0%) |
| SNV | 6 (1.4%) | 6 (1.4%) | 1 (0.0%) |
| none | 3 (0.7%) | 5 (1.2%) | 0 (0.0%) |
| complex | 27 (6.3%) | 27 (6.3%) | 0 (0.0%) |
| normal carries ALT | 33 (7.7%) | 33 (7.7%) | 0 (0.0%) |
| NA | 1 (0.2%) | 1 (0.2%) | 0 (0.0%) |
| no event | 0 (0.0%) | 0 (0.0%) | 21 (0.9%) |

### 5d. Union absorbed INDELs: units changed and homopolymer length, COLO829T vs HG008T

|  | COLO829T (c3 rule): of all absorbed INDELs | COLO829T (c3): of repeat / HP>=7 events | COLO829T (b5 rule): of all | COLO829T (b5): of repeat / HP>=7 events | HG008T (b5): of all absorbed INDELs | HG008T (b5): of repeat / HP>=7 events |
|---|---:|---:|---:|---:|---:|---:|
| 1 unit | 316 (73.5%) | 316 (97.2%) | 331 (77.0%) | 331 (97.4%) | 2,103 (92.4%) | 2,103 (97.9%) |
| 2-3 units | 7 (1.6%) | 7 (2.2%) | 7 (1.6%) | 7 (2.1%) | 38 (1.7%) | 38 (1.8%) |
| >3 units | 2 (0.5%) | 2 (0.6%) | 2 (0.5%) | 2 (0.6%) | 8 (0.4%) | 8 (0.4%) |
| HP 7-9 bp | 6 (1.4%) | 6 (2.5%) | 6 (1.4%) | 6 (2.4%) | 21 (0.9%) | 21 (1.1%) |
| HP 10-14 bp | 86 (20.0%) | 86 (35.7%) | 89 (20.7%) | 89 (35.6%) | 377 (16.6%) | 377 (20.3%) |
| HP 15-19 bp | 94 (21.9%) | 94 (39.0%) | 99 (23.0%) | 99 (39.6%) | 519 (22.8%) | 519 (28.0%) |
| HP 20-29 bp | 55 (12.8%) | 55 (22.8%) | 56 (13.0%) | 56 (22.4%) | 806 (35.4%) | 806 (43.5%) |
| HP >=30 bp | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 130 (5.7%) | 130 (7.0%) |

b5 rule = the HG008 rule on both sides (the same-rule comparison); the c3 rule is the COLO829T default of sections 2-5.

### 5e. Background, GRCh38 frame (b1 rule): share of all chr1-22 truth INDELs absorbed

| class | COLO829T | HG008T |
|---|---:|---:|
| HP>=7 | 244 / 1,167 (20.9%) | 1,688 / 6,360 (26.5%) |
| STR2-6>=3copies | 120 / 275 (43.6%) | 369 / 936 (39.4%) |
| HP4-6 | 3 / 97 (3.1%) | 15 / 241 (6.2%) |
| VNTR>6 | 2 / 4 (50.0%) | 2 / 6 (33.3%) |
| complex | 0 / 0 | 1 / 1 (100.0%) |
| none | 61 / 369 (16.5%) | 202 / 952 (21.2%) |
| total | 430 / 1,912 (22.5%) | 2,277 / 8,496 (26.8%) |

### 5f. Background, GRCh38 frame: share of all chr1-22 truth SNVs absorbed

| class | COLO829T | HG008T |
|---|---:|---:|
| HP>=7 | 4 / 545 (0.7%) | 73 / 399 (18.3%) |
| STR2-6>=3copies | 14 / 1,322 (1.1%) | 107 / 572 (18.7%) |
| HP4-6 | 7 / 5,914 (0.1%) | 4 / 880 (0.5%) |
| none | 42 / 34,254 (0.1%) | 17 / 6,839 (0.2%) |
| CpG Ti | 28 / 2,705 (1.0%) | 11 / 782 (1.4%) |
| total | 67 / 42,035 (0.2%) | 201 / 8,690 (2.3%) |

### 5g. Background, patient frame: share of all chr1-22 truth INDELs absorbed

| class | COLO829T (c3 rule) | COLO829T (b5 rule) | HG008T (b5 rule) |
|---|---:|---:|---:|
| HP>=7 | 241 / 1,103 (21.8%) | 250 / 1,160 (21.6%) | 1,853 / 6,722 (27.6%) |
| STR2-6>=3copies | 83 / 189 (43.9%) | 89 / 203 (43.8%) | 292 / 795 (36.7%) |
| HP4-6 | 0 / 67 (0.0%) | 0 / 87 (0.0%) | 1 / 204 (0.5%) |
| VNTR>6 | 1 / 13 (7.7%) | 1 / 13 (7.7%) | 3 / 15 (20.0%) |
| in_STR | 17 / 92 (18.5%) | 2 / 20 (10.0%) | 2 / 27 (7.4%) |
| imperfect VNTR | 2 / 4 (50.0%) | 0 / 0 | 0 / 0 |
| SNV in HP>=7 | 0 / 16 (0.0%) | 0 / 16 (0.0%) | 14 / 60 (23.3%) |
| SNV in STR2-6>=3copies | 15 / 46 (32.6%) | 15 / 46 (32.6%) | 90 / 192 (46.9%) |
| SNV in HP4-6 | 1 / 3 (33.3%) | 1 / 3 (33.3%) | 0 / 1 (0.0%) |
| SNV | 6 / 11 (54.5%) | 6 / 11 (54.5%) | 1 / 24 (4.2%) |
| none | 3 / 235 (1.3%) | 5 / 220 (2.3%) | 0 / 415 (0.0%) |
| complex | 27 / 72 (37.5%) | 27 / 72 (37.5%) | 0 / 0 |
| normal carries ALT | 33 / 55 (60.0%) | 33 / 55 (60.0%) | 0 / 0 |
| single hap | 0 / 3 (0.0%) | 0 / 3 (0.0%) | 0 / 0 |
| NA | 1 / 3 (33.3%) | 1 / 3 (33.3%) | 0 / 0 |
| no event | 0 / 0 | 0 / 0 | 21 / 41 (51.2%) |

Compare the two b5 columns (the same rule). HG008T has no imperfect-VNTR check, no single-hap label and no 'complex' / 'normal carries ALT' events (its INFO events are simple); 'no event' = no INFO event.

### 5h. Background, patient-frame mechanism class (INDEL records): share absorbed

| class | COLO829T (c3 rule) | COLO829T (b5 rule) | HG008T (b5 rule) |
|---|---:|---:|---:|
| HP>=7 1-unit, tract 7-9 | 6 / 116 (5.2%) | 6 / 129 (4.7%) | 21 / 453 (4.6%) |
| HP>=7 1-unit, tract 10-14 | 86 / 424 (20.3%) | 89 / 444 (20.0%) | 376 / 1,748 (21.5%) |
| HP>=7 1-unit, tract 15-19 | 94 / 375 (25.1%) | 99 / 387 (25.6%) | 517 / 1,795 (28.8%) |
| HP>=7 1-unit, tract 20-29 | 53 / 163 (32.5%) | 54 / 174 (31.0%) | 795 / 2,072 (38.4%) |
| HP>=7 1-unit, tract >=30 | 0 / 1 (0.0%) | 0 / 1 (0.0%) | 122 / 513 (23.8%) |
| STR 1-unit | 76 / 165 (46.1%) | 82 / 179 (45.8%) | 270 / 715 (37.8%) |
| repeat multi-unit | 9 / 48 (18.8%) | 9 / 49 (18.4%) | 44 / 221 (19.9%) |
| HP4-6 | 0 / 67 (0.0%) | 0 / 87 (0.0%) | 1 / 204 (0.5%) |
| VNTR>6 | 1 / 13 (7.7%) | 1 / 13 (7.7%) | 3 / 15 (20.0%) |
| in_STR | 17 / 92 (18.5%) | 2 / 20 (10.0%) | 2 / 27 (7.4%) |
| imperfect VNTR | 2 / 4 (50.0%) | 0 / 0 | 0 / 0 |
| SNV in repeat | 16 / 65 (24.6%) | 16 / 65 (24.6%) | 104 / 253 (41.1%) |
| non-repeat INDEL | 3 / 235 (1.3%) | 5 / 220 (2.3%) | 0 / 415 (0.0%) |
| non-repeat SNV | 6 / 11 (54.5%) | 6 / 11 (54.5%) | 1 / 24 (4.2%) |
| complex | 27 / 72 (37.5%) | 27 / 72 (37.5%) | 0 / 0 |
| normal carries ALT | 33 / 55 (60.0%) | 33 / 55 (60.0%) | 0 / 0 |
| unresolved | 1 / 6 (16.7%) | 1 / 6 (16.7%) | 21 / 41 (51.2%) |

Compare the two b5 columns (the same rule, tract = the b5 tract). HG008T row counts differ from mechanism.md, which counts SNV and INDEL records together; here INDEL records only.

