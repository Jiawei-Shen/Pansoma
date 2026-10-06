# Repeat context of the absorbed INDEL truths (s14)

Patient frame = the truth INFO HG008Nv62SOMATICVARIANT event in the HG008-N v6.2 event haplotype; GRCh38 frame = the truth VCF record in GRCh38. Classes in the s14 docstring.

### Patient frame: repeat class of the INDEL truths, per perfect set

| class | PacBio | ONT | Illumina | union |
|---|---:|---:|---:|---:|
| HP>=7 | 1,003 (73.6%) | 662 (66.9%) | 1,614 (84.2%) | 1,853 (81.4%) |
| STR2-6>=3copies | 248 (18.2%) | 221 (22.3%) | 220 (11.5%) | 292 (12.8%) |
| HP4-6 | 1 (0.1%) | 1 (0.1%) | 1 (0.1%) | 1 (0.0%) |
| VNTR>6 | 3 (0.2%) | 3 (0.3%) | 1 (0.1%) | 3 (0.1%) |
| in_STR | 2 (0.1%) | 2 (0.2%) | 0 (0.0%) | 2 (0.1%) |
| SNV in HP>=7 | 10 (0.7%) | 8 (0.8%) | 9 (0.5%) | 14 (0.6%) |
| SNV in STR2-6>=3copies | 75 (5.5%) | 74 (7.5%) | 60 (3.1%) | 90 (4.0%) |
| SNV | 1 (0.1%) | 1 (0.1%) | 0 (0.0%) | 1 (0.0%) |
| no event | 20 (1.5%) | 17 (1.7%) | 11 (0.6%) | 21 (0.9%) |
| total | 1,363 | 989 | 1,916 | 2,277 |

### Patient frame: units changed (repeat INDELs only)

| class | PacBio | ONT | Illumina | union |
|---|---:|---:|---:|---:|
| 1 unit | 1,222 (89.7%) | 855 (86.5%) | 1,808 (94.4%) | 2,103 (92.4%) |
| 2-3 units | 28 (2.1%) | 27 (2.7%) | 23 (1.2%) | 38 (1.7%) |
| >3 units | 5 (0.4%) | 5 (0.5%) | 5 (0.3%) | 8 (0.4%) |
| (not a repeat INDEL) | 108 (7.9%) | 102 (10.3%) | 80 (4.2%) | 128 (5.6%) |
| total | 1,363 | 989 | 1,916 | 2,277 |

### Patient frame: homopolymer length (HP>=7 INDELs only)

| class | PacBio | ONT | Illumina | union |
|---|---:|---:|---:|---:|
| 7-9 | 8 (0.6%) | 9 (0.9%) | 20 (1.0%) | 21 (0.9%) |
| 10-14 | 170 (12.5%) | 128 (12.9%) | 356 (18.6%) | 377 (16.6%) |
| 15-19 | 269 (19.7%) | 183 (18.5%) | 484 (25.3%) | 519 (22.8%) |
| 20-29 | 462 (33.9%) | 271 (27.4%) | 690 (36.0%) | 806 (35.4%) |
| >=30 | 94 (6.9%) | 71 (7.2%) | 64 (3.3%) | 130 (5.7%) |
| (not HP>=7 INDEL) | 360 (26.4%) | 327 (33.1%) | 302 (15.8%) | 424 (18.6%) |
| total | 1,363 | 989 | 1,916 | 2,277 |

### GRCh38 frame: repeat class of the INDEL truths, per perfect set

| class | PacBio | ONT | Illumina | union |
|---|---:|---:|---:|---:|
| HP>=7 | 876 (64.3%) | 546 (55.2%) | 1,486 (77.6%) | 1,688 (74.1%) |
| STR2-6>=3copies | 308 (22.6%) | 284 (28.7%) | 278 (14.5%) | 369 (16.2%) |
| HP4-6 | 10 (0.7%) | 9 (0.9%) | 12 (0.6%) | 15 (0.7%) |
| VNTR>6 | 2 (0.1%) | 2 (0.2%) | 1 (0.1%) | 2 (0.1%) |
| complex | 0 (0.0%) | 0 (0.0%) | 1 (0.1%) | 1 (0.0%) |
| none | 167 (12.3%) | 148 (15.0%) | 138 (7.2%) | 202 (8.9%) |
| total | 1,363 | 989 | 1,916 | 2,277 |

GRCh38-frame "none" INDELs (union) in the patient frame: HP>=7 77, STR2-6>=3copies 66, SNV in STR2-6>=3copies 42, no event 9, SNV in HP>=7 5, VNTR>6 2, in_STR 1

### Patient frame class x category (union INDELs)

| class | HG008N_present_broad | HG008N_present_rare | HG008N_absent_HPRC_other | other_ambiguous |
|---|---:|---:|---:|---:|
| HP>=7 | 67 | 22 | 1630 | 134 |
| STR2-6>=3copies | 8 | 3 | 225 | 56 |
| HP4-6 | 0 | 0 | 0 | 1 |
| VNTR>6 | 0 | 0 | 1 | 2 |
| in_STR | 0 | 1 | 0 | 1 |
| SNV in HP>=7 | 0 | 0 | 10 | 4 |
| SNV in STR2-6>=3copies | 0 | 3 | 57 | 30 |
| SNV | 0 | 0 | 1 | 0 |
| no event | 0 | 1 | 11 | 9 |

### All 8,496 chr1-22 truth INDELs (GRCh38 frame) and the share absorbed (union)

| class | truth INDELs | absorbed | % absorbed |
|---|---:|---:|---:|
| HP>=7 | 6,360 | 1,688 | 26.5% |
| STR2-6>=3copies | 936 | 369 | 39.4% |
| HP4-6 | 241 | 15 | 6.2% |
| VNTR>6 | 6 | 2 | 33.3% |
| complex | 1 | 1 | 100.0% |
| none | 952 | 202 | 21.2% |
| total | 8,496 | 2,277 | 26.8% |

Patient frame, same: HP>=7 1,853/6,722 (27.6%); STR2-6>=3copies 292/795 (36.7%); HP4-6 1/204 (0.5%); VNTR>6 3/15 (20.0%); in_STR 2/27 (7.4%); SNV in HP>=7 14/60 (23.3%); SNV in STR2-6>=3copies 90/192 (46.9%); SNV in HP4-6 0/1 (0.0%); SNV 1/24 (4.2%); none 0/415 (0.0%); no event 21/41 (51.2%)
