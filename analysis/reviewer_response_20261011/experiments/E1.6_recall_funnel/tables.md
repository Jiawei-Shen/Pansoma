# E1.6 recall funnel: tables (written 2026-10-15 19:51 by e16_recall_funnel.py)

chr1, truth records in the benchmark BED (rtg baseline), PASS at each run's validation recall-0.9 threshold; stage = the first stage at which a truth record is lost. Arms: GRCh38-only = linear, HPRC v1.1 d9, HPRC v2.1 d46. `(joint model)` = COLO829T fiberseq P100 called by the joint P100+P75+P50+P25 d9 model (supplementary, not arm-matched). Per-record rows: `/scratch/jshen/data/reviewer_experiments_20261013/E1.6_recall_funnel/per_truth.tsv.gz`.

## 0. Headline: % of the chr1 truth in the BED lost at each stage, per arm

Cell: `no tensor / tensor, not PASS / PASS, PoN-tagged / called` in % of the truth records (n in section 1). `n.l.` = no local single-set run.

| chr1 truth in BED | GRCh38-only | HPRC v1.1 d9 | HPRC v2.1 d46 |
|---|---|---|---|
| HG008T Illumina SNV | 1.9 / 14.5 / 5.0 / 78.6 | 5.3 / 8.5 / 3.2 / 83.1 | 4.4 / 6.9 / 3.4 / 85.2 |
| HG008T PacBio SNV | 2.0 / 9.8 / 5.6 / 82.6 | 5.5 / 8.9 / 2.9 / 82.8 | 5.0 / 7.9 / 3.6 / 83.5 |
| HG008T ONT SNV | 2.3 / 8.5 / 5.0 / 84.2 | 4.7 / 9.8 / 2.9 / 82.6 | 4.4 / 5.6 / 4.0 / 85.9 |
| HG008T Illumina INDEL | 19.4 / 6.8 / 47.5 / 26.2 | 59.8 / 4.3 / 12.1 / 23.9 | 56.7 / 6.0 / 13.3 / 24.0 |
| HG008T PacBio INDEL | 15.7 / 8.9 / 49.4 / 26.1 | 52.8 / 4.1 / 17.5 / 25.6 | 52.0 / 4.9 / 17.4 / 25.7 |
| HG008T ONT INDEL | 41.2 / 4.9 / 30.0 / 23.9 | 55.5 / 4.4 / 17.2 / 22.8 | 58.1 / 3.4 / 14.8 / 23.7 |
| COLO829T Illumina SNV | 16.7 / 23.4 / 1.7 / 58.2 | n.l. | 18.0 / 9.0 / 1.6 / 71.4 |
| COLO829T fiberseq SNV | 17.5 / 18.3 / 1.8 / 62.4 | n.l. | 17.8 / 9.8 / 1.6 / 70.7 |
| COLO829T ONT SNV | 17.8 / 26.6 / 1.4 / 54.2 | n.l. | 18.1 / 11.1 / 1.6 / 69.2 |
| COLO829T Illumina INDEL | 14.5 / 23.7 / 39.7 / 22.1 | n.l. | 53.4 / 7.6 / 14.5 / 24.4 |
| COLO829T fiberseq INDEL | 6.1 / 18.3 / 45.8 / 29.8 | n.l. | 48.1 / 7.6 / 18.3 / 26.0 |
| COLO829T ONT INDEL | 19.1 / 19.1 / 33.6 / 28.2 | n.l. | 57.3 / 4.6 / 15.3 / 22.9 |

## 1. Funnel per call set (chr1, in the BED; n (% of truth))

| call set (chr1) | arm | truth | 1 no tensor | 2 tensor, not PASS | 3 PASS, PoN-tagged | 4 called | threshold | rtg TP / FN after PoN (report.json) | stages reproduce rtg |
|---|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 2443 | 407 (16.7 %) | 572 (23.4 %) | 41 (1.7 %) | 1423 (58.2 %) | 0.0806 | 1423 / 1020 | yes |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 2443 | 440 (18.0 %) | 219 (9.0 %) | 39 (1.6 %) | 1745 (71.4 %) | 0.6395 | 1745 / 698 | yes |
| COLO829T fiberseq SNV | GRCh38-only | 2443 | 427 (17.5 %) | 448 (18.3 %) | 43 (1.8 %) | 1525 (62.4 %) | 0.0852 | 1525 / 918 | yes |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 2443 | 436 (17.8 %) | 239 (9.8 %) | 40 (1.6 %) | 1728 (70.7 %) | 0.4812 | 1728 / 715 | yes |
| COLO829T ONT SNV | GRCh38-only | 2443 | 434 (17.8 %) | 649 (26.6 %) | 35 (1.4 %) | 1325 (54.2 %) | 0.0920 | 1325 / 1118 | yes |
| COLO829T ONT SNV | HPRC v2.1 d46 | 2443 | 442 (18.1 %) | 271 (11.1 %) | 40 (1.6 %) | 1690 (69.2 %) | 0.5677 | 1690 / 753 | yes |
| HG008T Illumina SNV | GRCh38-only | 697 | 13 (1.9 %) | 101 (14.5 %) | 35 (5.0 %) | 548 (78.6 %) | 0.0640 | 548 / 149 | yes |
| HG008T Illumina SNV | HPRC v1.1 d9 | 697 | 37 (5.3 %) | 59 (8.5 %) | 22 (3.2 %) | 579 (83.1 %) | 0.2513 | 579 / 118 | yes |
| HG008T Illumina SNV | HPRC v2.1 d46 | 697 | 31 (4.4 %) | 48 (6.9 %) | 24 (3.4 %) | 594 (85.2 %) | 0.3237 | 594 / 103 | yes |
| HG008T PacBio SNV | GRCh38-only | 697 | 14 (2.0 %) | 68 (9.8 %) | 39 (5.6 %) | 576 (82.6 %) | 0.0397 | 576 / 121 | yes |
| HG008T PacBio SNV | HPRC v1.1 d9 | 697 | 38 (5.5 %) | 62 (8.9 %) | 20 (2.9 %) | 577 (82.8 %) | 0.3142 | 577 / 120 | yes |
| HG008T PacBio SNV | HPRC v2.1 d46 | 697 | 35 (5.0 %) | 55 (7.9 %) | 25 (3.6 %) | 582 (83.5 %) | 0.2912 | 582 / 115 | yes |
| HG008T ONT SNV | GRCh38-only | 697 | 16 (2.3 %) | 59 (8.5 %) | 35 (5.0 %) | 587 (84.2 %) | 0.0418 | 587 / 110 | yes |
| HG008T ONT SNV | HPRC v1.1 d9 | 697 | 33 (4.7 %) | 68 (9.8 %) | 20 (2.9 %) | 576 (82.6 %) | 0.4492 | 576 / 121 | yes |
| HG008T ONT SNV | HPRC v2.1 d46 | 697 | 31 (4.4 %) | 39 (5.6 %) | 28 (4.0 %) | 599 (85.9 %) | 0.2309 | 599 / 98 | yes |
| COLO829T Illumina INDEL | GRCh38-only | 131 | 19 (14.5 %) | 31 (23.7 %) | 52 (39.7 %) | 29 (22.1 %) | 0.0110 | 29 / 102 | yes |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 131 | 70 (53.4 %) | 10 (7.6 %) | 19 (14.5 %) | 32 (24.4 %) | 0.0222 | 32 / 99 | yes |
| COLO829T fiberseq INDEL | GRCh38-only | 131 | 8 (6.1 %) | 24 (18.3 %) | 60 (45.8 %) | 39 (29.8 %) | 0.0081 | 39 / 92 | yes |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 131 | 63 (48.1 %) | 10 (7.6 %) | 24 (18.3 %) | 34 (26.0 %) | 0.0117 | 34 / 97 | yes |
| COLO829T ONT INDEL | GRCh38-only | 131 | 25 (19.1 %) | 25 (19.1 %) | 44 (33.6 %) | 37 (28.2 %) | 0.0050 | 37 / 94 | yes |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 131 | 75 (57.3 %) | 6 (4.6 %) | 20 (15.3 %) | 30 (22.9 %) | 0.0127 | 30 / 101 | yes |
| HG008T Illumina INDEL | GRCh38-only | 587 | 114 (19.4 %) | 40 (6.8 %) | 279 (47.5 %) | 154 (26.2 %) | 0.0480 | 154 / 433 | yes |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 587 | 351 (59.8 %) | 25 (4.3 %) | 71 (12.1 %) | 140 (23.9 %) | 0.1195 | 140 / 447 | yes |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 587 | 333 (56.7 %) | 35 (6.0 %) | 78 (13.3 %) | 141 (24.0 %) | 0.1065 | 141 / 446 | yes |
| HG008T PacBio INDEL | GRCh38-only | 587 | 92 (15.7 %) | 52 (8.9 %) | 290 (49.4 %) | 153 (26.1 %) | 0.0367 | 153 / 434 | yes |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 587 | 310 (52.8 %) | 24 (4.1 %) | 103 (17.5 %) | 150 (25.6 %) | 0.0397 | 150 / 437 | yes |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 587 | 305 (52.0 %) | 29 (4.9 %) | 102 (17.4 %) | 151 (25.7 %) | 0.0270 | 151 / 436 | yes |
| HG008T ONT INDEL | GRCh38-only | 587 | 242 (41.2 %) | 29 (4.9 %) | 176 (30.0 %) | 140 (23.9 %) | 0.0084 | 140 / 447 | yes |
| HG008T ONT INDEL | HPRC v1.1 d9 | 587 | 326 (55.5 %) | 26 (4.4 %) | 101 (17.2 %) | 134 (22.8 %) | 0.0178 | 134 / 453 | yes |
| HG008T ONT INDEL | HPRC v2.1 d46 | 587 | 341 (58.1 %) | 20 (3.4 %) | 87 (14.8 %) | 139 (23.7 %) | 0.0135 | 139 / 448 | yes |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 2443 | 437 (17.9 %) | 251 (10.3 %) | 48 (2.0 %) | 1707 (69.9 %) | 0.5492 | 1707 / 736 | yes |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 131 | 68 (51.9 %) | 6 (4.6 %) | 21 (16.0 %) | 36 (27.5 %) | 0.0073 | 36 / 95 | yes |

Check: stage 4 = rtg TP-baseline after the PoN; 1+2+3 = FN after the PoN; 3 = TP before - TP after the PoN; 1+2 = FN before the PoN (report.json `pon` / `raw` `at_threshold`).

## 2. Stage 1 (no exact tensor) by reason (chr1; n)

| call set (chr1) | arm | stage 1 | AF floor (min_af) | min_variants | absorbed (graph spells ALT) | no candidate | other allele at site | other | inconsistent | with a partial-match tensor (PASS) |
|---|---|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 407 | 276 | 45 | 0 | 86 | 0 | 0 | 0 | 0 (0) |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 440 | 291 | 80 | 3 | 66 | 0 | 0 | 0 | 3 (1) |
| COLO829T fiberseq SNV | GRCh38-only | 427 | 302 | 48 | 0 | 77 | 0 | 0 | 0 | 0 (0) |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 436 | 321 | 62 | 3 | 50 | 0 | 0 | 0 | 1 (1) |
| COLO829T ONT SNV | GRCh38-only | 434 | 229 | 160 | 0 | 45 | 0 | 0 | 0 | 0 (0) |
| COLO829T ONT SNV | HPRC v2.1 d46 | 442 | 226 | 166 | 3 | 47 | 0 | 0 | 0 | 1 (1) |
| HG008T Illumina SNV | GRCh38-only | 13 | 2 | 0 | 0 | 11 | 0 | 0 | 0 | 0 (0) |
| HG008T Illumina SNV | HPRC v1.1 d9 | 37 | 1 | 4 | 11 | 21 | 0 | 0 | 0 | 3 (3) |
| HG008T Illumina SNV | HPRC v2.1 d46 | 31 | 2 | 2 | 10 | 17 | 0 | 0 | 0 | 2 (1) |
| HG008T PacBio SNV | GRCh38-only | 14 | 2 | 4 | 0 | 8 | 0 | 0 | 0 | 0 (0) |
| HG008T PacBio SNV | HPRC v1.1 d9 | 38 | 1 | 4 | 12 | 21 | 0 | 0 | 0 | 3 (0) |
| HG008T PacBio SNV | HPRC v2.1 d46 | 35 | 1 | 2 | 10 | 22 | 0 | 0 | 0 | 2 (1) |
| HG008T ONT SNV | GRCh38-only | 16 | 0 | 12 | 0 | 4 | 0 | 0 | 0 | 0 (0) |
| HG008T ONT SNV | HPRC v1.1 d9 | 33 | 0 | 5 | 12 | 16 | 0 | 0 | 0 | 4 (1) |
| HG008T ONT SNV | HPRC v2.1 d46 | 31 | 1 | 6 | 10 | 14 | 0 | 0 | 0 | 1 (1) |
| COLO829T Illumina INDEL | GRCh38-only | 19 | 2 | 8 | 0 | 8 | 1 | 0 | 0 | 1 (0) |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 70 | 1 | 3 | 34 | 32 | 0 | 0 | 0 | 14 (6) |
| COLO829T fiberseq INDEL | GRCh38-only | 8 | 1 | 1 | 0 | 4 | 2 | 0 | 0 | 2 (2) |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 63 | 6 | 0 | 26 | 27 | 4 | 0 | 0 | 25 (17) |
| COLO829T ONT INDEL | GRCh38-only | 25 | 5 | 1 | 0 | 0 | 19 | 0 | 0 | 12 (8) |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 75 | 10 | 5 | 21 | 24 | 15 | 0 | 0 | 32 (26) |
| HG008T Illumina INDEL | GRCh38-only | 114 | 33 | 14 | 0 | 34 | 33 | 0 | 0 | 14 (3) |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 351 | 20 | 46 | 137 | 145 | 3 | 0 | 0 | 82 (48) |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 333 | 15 | 46 | 108 | 157 | 7 | 0 | 0 | 106 (64) |
| HG008T PacBio INDEL | GRCh38-only | 92 | 7 | 3 | 0 | 11 | 71 | 0 | 0 | 19 (15) |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 310 | 40 | 29 | 81 | 144 | 16 | 0 | 0 | 137 (81) |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 305 | 26 | 19 | 77 | 164 | 19 | 0 | 0 | 164 (124) |
| HG008T ONT INDEL | GRCh38-only | 242 | 22 | 49 | 0 | 16 | 155 | 0 | 0 | 96 (86) |
| HG008T ONT INDEL | HPRC v1.1 d9 | 326 | 13 | 83 | 57 | 128 | 45 | 0 | 0 | 148 (125) |
| HG008T ONT INDEL | HPRC v2.1 d46 | 341 | 14 | 59 | 64 | 144 | 60 | 0 | 0 | 181 (149) |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 437 | 328 | 69 | 4 | 36 | 0 | 0 | 0 | 0 (0) |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 68 | 10 | 1 | 26 | 27 | 4 | 0 | 0 | 20 (12) |

Reasons from `<set>/tensors/somatic.recall.tsv`; absorbed = no_candidate and the arm's graph spells the ALT (E5b spell tables; never on GRCh38-only). Partial-match tensor: labelled somatic for this truth by an overlap match (another allele / representation); its call cannot be an rtg match.

## 3. Stage 2 (exact tensor, no matching PASS call) by reason (chr1; n)

| call set (chr1) | arm | stage 2 | p_somatic below threshold | PASS, no GRCh38 position | PASS, placed, rtg FN |
|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 572 | 572 | 0 | 0 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 219 | 219 | 0 | 0 |
| COLO829T fiberseq SNV | GRCh38-only | 448 | 448 | 0 | 0 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 239 | 239 | 0 | 0 |
| COLO829T ONT SNV | GRCh38-only | 649 | 649 | 0 | 0 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 271 | 271 | 0 | 0 |
| HG008T Illumina SNV | GRCh38-only | 101 | 101 | 0 | 0 |
| HG008T Illumina SNV | HPRC v1.1 d9 | 59 | 59 | 0 | 0 |
| HG008T Illumina SNV | HPRC v2.1 d46 | 48 | 48 | 0 | 0 |
| HG008T PacBio SNV | GRCh38-only | 68 | 68 | 0 | 0 |
| HG008T PacBio SNV | HPRC v1.1 d9 | 62 | 62 | 0 | 0 |
| HG008T PacBio SNV | HPRC v2.1 d46 | 55 | 55 | 0 | 0 |
| HG008T ONT SNV | GRCh38-only | 59 | 59 | 0 | 0 |
| HG008T ONT SNV | HPRC v1.1 d9 | 68 | 68 | 0 | 0 |
| HG008T ONT SNV | HPRC v2.1 d46 | 39 | 39 | 0 | 0 |
| COLO829T Illumina INDEL | GRCh38-only | 31 | 31 | 0 | 0 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 10 | 10 | 0 | 0 |
| COLO829T fiberseq INDEL | GRCh38-only | 24 | 24 | 0 | 0 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 10 | 10 | 0 | 0 |
| COLO829T ONT INDEL | GRCh38-only | 25 | 25 | 0 | 0 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 6 | 6 | 0 | 0 |
| HG008T Illumina INDEL | GRCh38-only | 40 | 40 | 0 | 0 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 25 | 25 | 0 | 0 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 35 | 35 | 0 | 0 |
| HG008T PacBio INDEL | GRCh38-only | 52 | 52 | 0 | 0 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 24 | 24 | 0 | 0 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 29 | 29 | 0 | 0 |
| HG008T ONT INDEL | GRCh38-only | 29 | 29 | 0 | 0 |
| HG008T ONT INDEL | HPRC v1.1 d9 | 26 | 26 | 0 | 0 |
| HG008T ONT INDEL | HPRC v2.1 d46 | 20 | 20 | 0 | 0 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 251 | 251 | 0 | 0 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 6 | 6 | 0 | 0 |

## 4. PoN tags of the truth records per stage (chr1; tagged / stage n)

| call set (chr1) | arm | truth tagged | stage 1 | stage 2 | stage 3 | stage 4 | stage 3: exact tensor's call tagged | untagged truth n | untagged: 1 % | 2 % | 3 % | 4 % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 84 / 2443 (3.4 %) | 11 / 407 | 32 / 572 | 41 / 41 | 0 / 1423 | 41 / 41 | 2359 | 16.8 | 22.9 | 0.0 | 60.3 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 84 / 2443 (3.4 %) | 23 / 440 | 22 / 219 | 39 / 39 | 0 / 1745 | 39 / 39 | 2359 | 17.7 | 8.4 | 0.0 | 74.0 |
| COLO829T fiberseq SNV | GRCh38-only | 84 / 2443 (3.4 %) | 21 / 427 | 20 / 448 | 43 / 43 | 0 / 1525 | 43 / 43 | 2359 | 17.2 | 18.1 | 0.0 | 64.6 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 84 / 2443 (3.4 %) | 24 / 436 | 20 / 239 | 40 / 40 | 0 / 1728 | 40 / 40 | 2359 | 17.5 | 9.3 | 0.0 | 73.3 |
| COLO829T ONT SNV | GRCh38-only | 84 / 2443 (3.4 %) | 22 / 434 | 27 / 649 | 35 / 35 | 0 / 1325 | 35 / 35 | 2359 | 17.5 | 26.4 | 0.0 | 56.2 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 84 / 2443 (3.4 %) | 25 / 442 | 19 / 271 | 40 / 40 | 0 / 1690 | 40 / 40 | 2359 | 17.7 | 10.7 | 0.0 | 71.6 |
| HG008T Illumina SNV | GRCh38-only | 60 / 697 (8.6 %) | 4 / 13 | 21 / 101 | 35 / 35 | 0 / 548 | 35 / 35 | 637 | 1.4 | 12.6 | 0.0 | 86.0 |
| HG008T Illumina SNV | HPRC v1.1 d9 | 60 / 697 (8.6 %) | 23 / 37 | 15 / 59 | 22 / 22 | 0 / 579 | 22 / 22 | 637 | 2.2 | 6.9 | 0.0 | 90.9 |
| HG008T Illumina SNV | HPRC v2.1 d46 | 60 / 697 (8.6 %) | 20 / 31 | 16 / 48 | 24 / 24 | 0 / 594 | 24 / 24 | 637 | 1.7 | 5.0 | 0.0 | 93.2 |
| HG008T PacBio SNV | GRCh38-only | 60 / 697 (8.6 %) | 8 / 14 | 13 / 68 | 39 / 39 | 0 / 576 | 39 / 39 | 637 | 0.9 | 8.6 | 0.0 | 90.4 |
| HG008T PacBio SNV | HPRC v1.1 d9 | 60 / 697 (8.6 %) | 24 / 38 | 16 / 62 | 20 / 20 | 0 / 577 | 20 / 20 | 637 | 2.2 | 7.2 | 0.0 | 90.6 |
| HG008T PacBio SNV | HPRC v2.1 d46 | 60 / 697 (8.6 %) | 21 / 35 | 14 / 55 | 25 / 25 | 0 / 582 | 25 / 25 | 637 | 2.2 | 6.4 | 0.0 | 91.4 |
| HG008T ONT SNV | GRCh38-only | 60 / 697 (8.6 %) | 12 / 16 | 13 / 59 | 35 / 35 | 0 / 587 | 35 / 35 | 637 | 0.6 | 7.2 | 0.0 | 92.2 |
| HG008T ONT SNV | HPRC v1.1 d9 | 60 / 697 (8.6 %) | 23 / 33 | 17 / 68 | 20 / 20 | 0 / 576 | 20 / 20 | 637 | 1.6 | 8.0 | 0.0 | 90.4 |
| HG008T ONT SNV | HPRC v2.1 d46 | 60 / 697 (8.6 %) | 24 / 31 | 8 / 39 | 28 / 28 | 0 / 599 | 28 / 28 | 637 | 1.1 | 4.9 | 0.0 | 94.0 |
| COLO829T Illumina INDEL | GRCh38-only | 81 / 131 (61.8 %) | 16 / 19 | 13 / 31 | 52 / 52 | 0 / 29 | 52 / 52 | 50 | 6.0 | 36.0 | 0.0 | 58.0 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 81 / 131 (61.8 %) | 56 / 70 | 6 / 10 | 19 / 19 | 0 / 32 | 19 / 19 | 50 | 28.0 | 8.0 | 0.0 | 64.0 |
| COLO829T fiberseq INDEL | GRCh38-only | 81 / 131 (61.8 %) | 8 / 8 | 13 / 24 | 60 / 60 | 0 / 39 | 60 / 60 | 50 | 0.0 | 22.0 | 0.0 | 78.0 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 81 / 131 (61.8 %) | 50 / 63 | 7 / 10 | 24 / 24 | 0 / 34 | 24 / 24 | 50 | 26.0 | 6.0 | 0.0 | 68.0 |
| COLO829T ONT INDEL | GRCh38-only | 81 / 131 (61.8 %) | 19 / 25 | 18 / 25 | 44 / 44 | 0 / 37 | 44 / 44 | 50 | 12.0 | 14.0 | 0.0 | 74.0 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 81 / 131 (61.8 %) | 59 / 75 | 2 / 6 | 20 / 20 | 0 / 30 | 20 / 20 | 50 | 32.0 | 8.0 | 0.0 | 60.0 |
| HG008T Illumina INDEL | GRCh38-only | 379 / 587 (64.6 %) | 76 / 114 | 24 / 40 | 279 / 279 | 0 / 154 | 279 / 279 | 208 | 18.3 | 7.7 | 0.0 | 74.0 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 379 / 587 (64.6 %) | 294 / 351 | 14 / 25 | 71 / 71 | 0 / 140 | 71 / 71 | 208 | 27.4 | 5.3 | 0.0 | 67.3 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 379 / 587 (64.6 %) | 278 / 333 | 23 / 35 | 78 / 78 | 0 / 141 | 78 / 78 | 208 | 26.4 | 5.8 | 0.0 | 67.8 |
| HG008T PacBio INDEL | GRCh38-only | 379 / 587 (64.6 %) | 64 / 92 | 25 / 52 | 290 / 290 | 0 / 153 | 290 / 290 | 208 | 13.5 | 13.0 | 0.0 | 73.6 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 379 / 587 (64.6 %) | 255 / 310 | 21 / 24 | 103 / 103 | 0 / 150 | 103 / 103 | 208 | 26.4 | 1.4 | 0.0 | 72.1 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 379 / 587 (64.6 %) | 250 / 305 | 27 / 29 | 102 / 102 | 0 / 151 | 102 / 102 | 208 | 26.4 | 1.0 | 0.0 | 72.6 |
| HG008T ONT INDEL | GRCh38-only | 379 / 587 (64.6 %) | 191 / 242 | 12 / 29 | 176 / 176 | 0 / 140 | 176 / 176 | 208 | 24.5 | 8.2 | 0.0 | 67.3 |
| HG008T ONT INDEL | HPRC v1.1 d9 | 379 / 587 (64.6 %) | 266 / 326 | 12 / 26 | 101 / 101 | 0 / 134 | 99 / 101 | 208 | 28.8 | 6.7 | 0.0 | 64.4 |
| HG008T ONT INDEL | HPRC v2.1 d46 | 379 / 587 (64.6 %) | 278 / 341 | 14 / 20 | 87 / 87 | 0 / 139 | 86 / 87 | 208 | 30.3 | 2.9 | 0.0 | 66.8 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 84 / 2443 (3.4 %) | 24 / 437 | 12 / 251 | 48 / 48 | 0 / 1707 | 48 / 48 | 2359 | 17.5 | 10.1 | 0.0 | 72.4 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 81 / 131 (61.8 %) | 55 / 68 | 5 / 6 | 21 / 21 | 0 / 36 | 21 / 21 | 50 | 26.0 | 2.0 | 0.0 | 72.0 |

Truth tag: the PoN run by vcfeval over the truth records (`eval_bed/truth.<KIND>.pon.vcf.gz`, allele-exact). The untagged columns give the funnel of the truth the PoN does not tag (the share any caller can still report after the PoN).

## 5. Funnel by variant kind (chr1; % of the kind's truth)

| call set (chr1) | arm | kind | n | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | SNV | 2443 | 16.7 | 23.4 | 1.7 | 58.2 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | SNV | 2443 | 18.0 | 9.0 | 1.6 | 71.4 |
| COLO829T fiberseq SNV | GRCh38-only | SNV | 2443 | 17.5 | 18.3 | 1.8 | 62.4 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | SNV | 2443 | 17.8 | 9.8 | 1.6 | 70.7 |
| COLO829T ONT SNV | GRCh38-only | SNV | 2443 | 17.8 | 26.6 | 1.4 | 54.2 |
| COLO829T ONT SNV | HPRC v2.1 d46 | SNV | 2443 | 18.1 | 11.1 | 1.6 | 69.2 |
| HG008T Illumina SNV | GRCh38-only | SNV | 697 | 1.9 | 14.5 | 5.0 | 78.6 |
| HG008T Illumina SNV | HPRC v1.1 d9 | SNV | 697 | 5.3 | 8.5 | 3.2 | 83.1 |
| HG008T Illumina SNV | HPRC v2.1 d46 | SNV | 697 | 4.4 | 6.9 | 3.4 | 85.2 |
| HG008T PacBio SNV | GRCh38-only | SNV | 697 | 2.0 | 9.8 | 5.6 | 82.6 |
| HG008T PacBio SNV | HPRC v1.1 d9 | SNV | 697 | 5.5 | 8.9 | 2.9 | 82.8 |
| HG008T PacBio SNV | HPRC v2.1 d46 | SNV | 697 | 5.0 | 7.9 | 3.6 | 83.5 |
| HG008T ONT SNV | GRCh38-only | SNV | 697 | 2.3 | 8.5 | 5.0 | 84.2 |
| HG008T ONT SNV | HPRC v1.1 d9 | SNV | 697 | 4.7 | 9.8 | 2.9 | 82.6 |
| HG008T ONT SNV | HPRC v2.1 d46 | SNV | 697 | 4.4 | 5.6 | 4.0 | 85.9 |
| COLO829T Illumina INDEL | GRCh38-only | INS | 77 | 14.3 | 26.0 | 42.9 | 16.9 |
| COLO829T Illumina INDEL | GRCh38-only | DEL | 54 | 14.8 | 20.4 | 35.2 | 29.6 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | INS | 77 | 58.4 | 6.5 | 11.7 | 23.4 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | DEL | 54 | 46.3 | 9.3 | 18.5 | 25.9 |
| COLO829T fiberseq INDEL | GRCh38-only | INS | 77 | 3.9 | 19.5 | 49.4 | 27.3 |
| COLO829T fiberseq INDEL | GRCh38-only | DEL | 54 | 9.3 | 16.7 | 40.7 | 33.3 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | INS | 77 | 55.8 | 3.9 | 16.9 | 23.4 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | DEL | 54 | 37.0 | 13.0 | 20.4 | 29.6 |
| COLO829T ONT INDEL | GRCh38-only | INS | 77 | 15.6 | 23.4 | 37.7 | 23.4 |
| COLO829T ONT INDEL | GRCh38-only | DEL | 54 | 24.1 | 13.0 | 27.8 | 35.2 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | INS | 77 | 62.3 | 1.3 | 15.6 | 20.8 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | DEL | 54 | 50.0 | 9.3 | 14.8 | 25.9 |
| HG008T Illumina INDEL | GRCh38-only | INS | 368 | 21.2 | 7.1 | 47.0 | 24.7 |
| HG008T Illumina INDEL | GRCh38-only | DEL | 219 | 16.4 | 6.4 | 48.4 | 28.8 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | INS | 368 | 65.2 | 3.3 | 10.3 | 21.2 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | DEL | 219 | 50.7 | 5.9 | 15.1 | 28.3 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | INS | 368 | 60.3 | 6.2 | 12.0 | 21.5 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | DEL | 219 | 50.7 | 5.5 | 15.5 | 28.3 |
| HG008T PacBio INDEL | GRCh38-only | INS | 368 | 18.2 | 7.3 | 48.4 | 26.1 |
| HG008T PacBio INDEL | GRCh38-only | DEL | 219 | 11.4 | 11.4 | 51.1 | 26.0 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | INS | 368 | 59.8 | 2.4 | 14.4 | 23.4 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | DEL | 219 | 41.1 | 6.8 | 22.8 | 29.2 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | INS | 368 | 58.4 | 2.2 | 16.0 | 23.4 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | DEL | 219 | 41.1 | 9.6 | 19.6 | 29.7 |
| HG008T ONT INDEL | GRCh38-only | INS | 368 | 44.8 | 5.4 | 27.4 | 22.3 |
| HG008T ONT INDEL | GRCh38-only | DEL | 219 | 35.2 | 4.1 | 34.2 | 26.5 |
| HG008T ONT INDEL | HPRC v1.1 d9 | INS | 368 | 62.0 | 2.2 | 15.2 | 20.7 |
| HG008T ONT INDEL | HPRC v1.1 d9 | DEL | 219 | 44.7 | 8.2 | 20.5 | 26.5 |
| HG008T ONT INDEL | HPRC v2.1 d46 | INS | 368 | 63.3 | 3.3 | 13.0 | 20.4 |
| HG008T ONT INDEL | HPRC v2.1 d46 | DEL | 219 | 49.3 | 3.7 | 17.8 | 29.2 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | SNV | 2443 | 17.9 | 10.3 | 2.0 | 69.9 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | INS | 77 | 61.0 | 2.6 | 14.3 | 22.1 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | DEL | 54 | 38.9 | 7.4 | 18.5 | 35.2 |

## 6. Funnel by truth VAF bin (chr1; % of the bin's truth)

VAF source order: own tensor AF, own build filter min_af (< 0.08), other-arm tensor AF (same reads), COLO829T truth VCF (VAF_Ill / VAF_PB), HG008T other-platform tensor AF (median), unknown. Section 7 gives the sources.

| call set (chr1) | arm | VAF bin | n | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 0-0.08 | 402 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T Illumina SNV | GRCh38-only | 0.08-0.2 | 196 | 1.5 | 3.6 | 3.1 | 91.8 |
| COLO829T Illumina SNV | GRCh38-only | 0.2-0.4 | 103 | 0.0 | 17.5 | 1.0 | 81.6 |
| COLO829T Illumina SNV | GRCh38-only | >0.4 | 1742 | 0.1 | 31.4 | 2.0 | 66.5 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 0-0.08 | 401 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 0.08-0.2 | 194 | 4.6 | 3.6 | 2.1 | 89.7 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 0.2-0.4 | 101 | 6.9 | 15.8 | 1.0 | 76.2 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | >0.4 | 1747 | 1.3 | 11.2 | 1.9 | 85.5 |
| COLO829T fiberseq SNV | GRCh38-only | 0-0.08 | 413 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T fiberseq SNV | GRCh38-only | 0.08-0.2 | 188 | 2.1 | 1.6 | 2.1 | 94.1 |
| COLO829T fiberseq SNV | GRCh38-only | 0.2-0.4 | 102 | 3.9 | 10.8 | 2.9 | 82.4 |
| COLO829T fiberseq SNV | GRCh38-only | >0.4 | 1740 | 0.3 | 24.9 | 2.1 | 72.6 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 0-0.08 | 409 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 0.08-0.2 | 188 | 4.8 | 2.1 | 1.6 | 91.5 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 0.2-0.4 | 102 | 3.9 | 13.7 | 2.0 | 80.4 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | >0.4 | 1744 | 0.8 | 12.7 | 2.0 | 84.5 |
| COLO829T ONT SNV | GRCh38-only | 0-0.08 | 422 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T ONT SNV | GRCh38-only | 0.08-0.2 | 217 | 3.2 | 6.0 | 2.3 | 88.5 |
| COLO829T ONT SNV | GRCh38-only | 0.2-0.4 | 149 | 0.7 | 10.7 | 2.7 | 85.9 |
| COLO829T ONT SNV | GRCh38-only | >0.4 | 1655 | 0.2 | 37.5 | 1.6 | 60.7 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 0-0.08 | 421 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 0.08-0.2 | 219 | 5.0 | 12.3 | 1.8 | 80.8 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 0.2-0.4 | 155 | 0.6 | 16.8 | 1.9 | 80.6 |
| COLO829T ONT SNV | HPRC v2.1 d46 | >0.4 | 1648 | 0.5 | 13.2 | 2.0 | 84.2 |
| HG008T Illumina SNV | GRCh38-only | 0-0.08 | 2 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina SNV | GRCh38-only | 0.08-0.2 | 7 | 42.9 | 42.9 | 14.3 | 0.0 |
| HG008T Illumina SNV | GRCh38-only | 0.2-0.4 | 69 | 7.2 | 8.7 | 8.7 | 75.4 |
| HG008T Illumina SNV | GRCh38-only | >0.4 | 618 | 0.3 | 14.9 | 4.5 | 80.3 |
| HG008T Illumina SNV | GRCh38-only | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina SNV | HPRC v1.1 d9 | 0-0.08 | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina SNV | HPRC v1.1 d9 | 0.08-0.2 | 7 | 85.7 | 14.3 | 0.0 | 0.0 |
| HG008T Illumina SNV | HPRC v1.1 d9 | 0.2-0.4 | 64 | 15.6 | 10.9 | 3.1 | 70.3 |
| HG008T Illumina SNV | HPRC v1.1 d9 | >0.4 | 624 | 3.0 | 8.2 | 3.2 | 85.6 |
| HG008T Illumina SNV | HPRC v1.1 d9 | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina SNV | HPRC v2.1 d46 | 0-0.08 | 2 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina SNV | HPRC v2.1 d46 | 0.08-0.2 | 4 | 75.0 | 25.0 | 0.0 | 0.0 |
| HG008T Illumina SNV | HPRC v2.1 d46 | 0.2-0.4 | 66 | 13.6 | 12.1 | 4.5 | 69.7 |
| HG008T Illumina SNV | HPRC v2.1 d46 | >0.4 | 624 | 2.6 | 6.2 | 3.4 | 87.8 |
| HG008T Illumina SNV | HPRC v2.1 d46 | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio SNV | GRCh38-only | 0-0.08 | 2 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio SNV | GRCh38-only | 0.08-0.2 | 12 | 16.7 | 50.0 | 16.7 | 16.7 |
| HG008T PacBio SNV | GRCh38-only | 0.2-0.4 | 76 | 7.9 | 5.3 | 11.8 | 75.0 |
| HG008T PacBio SNV | GRCh38-only | >0.4 | 606 | 0.5 | 9.6 | 4.6 | 85.3 |
| HG008T PacBio SNV | GRCh38-only | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio SNV | HPRC v1.1 d9 | 0-0.08 | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio SNV | HPRC v1.1 d9 | 0.08-0.2 | 14 | 71.4 | 28.6 | 0.0 | 0.0 |
| HG008T PacBio SNV | HPRC v1.1 d9 | 0.2-0.4 | 77 | 15.6 | 11.7 | 3.9 | 68.8 |
| HG008T PacBio SNV | HPRC v1.1 d9 | >0.4 | 604 | 2.3 | 8.1 | 2.8 | 86.8 |
| HG008T PacBio SNV | HPRC v1.1 d9 | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio SNV | HPRC v2.1 d46 | 0-0.08 | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio SNV | HPRC v2.1 d46 | 0.08-0.2 | 10 | 80.0 | 20.0 | 0.0 | 0.0 |
| HG008T PacBio SNV | HPRC v2.1 d46 | 0.2-0.4 | 74 | 16.2 | 12.2 | 2.7 | 68.9 |
| HG008T PacBio SNV | HPRC v2.1 d46 | >0.4 | 611 | 2.1 | 7.2 | 3.8 | 86.9 |
| HG008T PacBio SNV | HPRC v2.1 d46 | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT SNV | GRCh38-only | 0.08-0.2 | 19 | 21.1 | 47.4 | 15.8 | 15.8 |
| HG008T ONT SNV | GRCh38-only | 0.2-0.4 | 100 | 7.0 | 5.0 | 8.0 | 80.0 |
| HG008T ONT SNV | GRCh38-only | >0.4 | 577 | 0.7 | 7.8 | 4.2 | 87.3 |
| HG008T ONT SNV | GRCh38-only | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT SNV | HPRC v1.1 d9 | 0.08-0.2 | 18 | 55.6 | 38.9 | 0.0 | 5.6 |
| HG008T ONT SNV | HPRC v1.1 d9 | 0.2-0.4 | 99 | 10.1 | 10.1 | 3.0 | 76.8 |
| HG008T ONT SNV | HPRC v1.1 d9 | >0.4 | 579 | 2.1 | 8.8 | 2.9 | 86.2 |
| HG008T ONT SNV | HPRC v1.1 d9 | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT SNV | HPRC v2.1 d46 | 0-0.08 | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT SNV | HPRC v2.1 d46 | 0.08-0.2 | 13 | 53.8 | 30.8 | 0.0 | 15.4 |
| HG008T ONT SNV | HPRC v2.1 d46 | 0.2-0.4 | 105 | 9.5 | 4.8 | 4.8 | 81.0 |
| HG008T ONT SNV | HPRC v2.1 d46 | >0.4 | 577 | 2.1 | 5.2 | 4.0 | 88.7 |
| HG008T ONT SNV | HPRC v2.1 d46 | unknown | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T Illumina INDEL | GRCh38-only | 0-0.08 | 2 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T Illumina INDEL | GRCh38-only | 0.08-0.2 | 8 | 0.0 | 50.0 | 25.0 | 25.0 |
| COLO829T Illumina INDEL | GRCh38-only | 0.2-0.4 | 47 | 17.0 | 4.3 | 46.8 | 31.9 |
| COLO829T Illumina INDEL | GRCh38-only | >0.4 | 74 | 12.2 | 33.8 | 37.8 | 16.2 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 0-0.08 | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 0.08-0.2 | 9 | 44.4 | 33.3 | 0.0 | 22.2 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 0.2-0.4 | 44 | 54.5 | 6.8 | 18.2 | 20.5 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | >0.4 | 77 | 53.2 | 5.2 | 14.3 | 27.3 |
| COLO829T fiberseq INDEL | GRCh38-only | 0-0.08 | 1 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T fiberseq INDEL | GRCh38-only | 0.08-0.2 | 4 | 0.0 | 0.0 | 75.0 | 25.0 |
| COLO829T fiberseq INDEL | GRCh38-only | 0.2-0.4 | 22 | 4.5 | 18.2 | 40.9 | 36.4 |
| COLO829T fiberseq INDEL | GRCh38-only | >0.4 | 104 | 5.8 | 19.2 | 46.2 | 28.8 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 0-0.08 | 6 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 0.08-0.2 | 13 | 38.5 | 46.2 | 15.4 | 0.0 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 0.2-0.4 | 23 | 56.5 | 4.3 | 13.0 | 26.1 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | >0.4 | 89 | 43.8 | 3.4 | 21.3 | 31.5 |
| COLO829T ONT INDEL | GRCh38-only | 0-0.08 | 5 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T ONT INDEL | GRCh38-only | 0.08-0.2 | 32 | 50.0 | 25.0 | 9.4 | 15.6 |
| COLO829T ONT INDEL | GRCh38-only | 0.2-0.4 | 50 | 6.0 | 14.0 | 54.0 | 26.0 |
| COLO829T ONT INDEL | GRCh38-only | >0.4 | 44 | 2.3 | 22.7 | 31.8 | 43.2 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 0-0.08 | 10 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 0.08-0.2 | 37 | 62.2 | 13.5 | 16.2 | 8.1 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 0.2-0.4 | 42 | 50.0 | 0.0 | 26.2 | 23.8 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | >0.4 | 42 | 50.0 | 2.4 | 7.1 | 40.5 |
| HG008T Illumina INDEL | GRCh38-only | 0-0.08 | 33 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina INDEL | GRCh38-only | 0.08-0.2 | 82 | 26.8 | 9.8 | 51.2 | 12.2 |
| HG008T Illumina INDEL | GRCh38-only | 0.2-0.4 | 224 | 17.4 | 3.1 | 56.2 | 23.2 |
| HG008T Illumina INDEL | GRCh38-only | >0.4 | 234 | 2.6 | 10.7 | 47.4 | 39.3 |
| HG008T Illumina INDEL | GRCh38-only | unknown | 14 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 0-0.08 | 20 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 0.08-0.2 | 98 | 85.7 | 11.2 | 3.1 | 0.0 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 0.2-0.4 | 218 | 59.6 | 2.3 | 16.1 | 22.0 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | >0.4 | 236 | 43.2 | 3.8 | 14.0 | 39.0 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | unknown | 15 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 0-0.08 | 15 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 0.08-0.2 | 103 | 81.6 | 11.7 | 4.9 | 1.9 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 0.2-0.4 | 223 | 57.8 | 5.4 | 16.6 | 20.2 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | >0.4 | 231 | 39.0 | 4.8 | 15.6 | 40.7 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | unknown | 15 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio INDEL | GRCh38-only | 0-0.08 | 7 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio INDEL | GRCh38-only | 0.08-0.2 | 38 | 57.9 | 0.0 | 34.2 | 7.9 |
| HG008T PacBio INDEL | GRCh38-only | 0.2-0.4 | 218 | 21.6 | 2.3 | 51.8 | 24.3 |
| HG008T PacBio INDEL | GRCh38-only | >0.4 | 313 | 1.6 | 15.0 | 52.4 | 31.0 |
| HG008T PacBio INDEL | GRCh38-only | unknown | 11 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 0-0.08 | 40 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 0.08-0.2 | 79 | 57.0 | 19.0 | 22.8 | 1.3 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 0.2-0.4 | 188 | 55.3 | 3.2 | 21.3 | 20.2 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | >0.4 | 266 | 40.2 | 1.1 | 16.9 | 41.7 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | unknown | 14 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 0-0.08 | 26 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 0.08-0.2 | 69 | 66.7 | 23.2 | 10.1 | 0.0 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 0.2-0.4 | 202 | 55.0 | 5.4 | 20.8 | 18.8 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | >0.4 | 276 | 39.1 | 0.7 | 19.2 | 40.9 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | unknown | 14 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT INDEL | GRCh38-only | 0-0.08 | 22 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT INDEL | GRCh38-only | 0.08-0.2 | 220 | 66.4 | 3.2 | 23.2 | 7.3 |
| HG008T ONT INDEL | GRCh38-only | 0.2-0.4 | 199 | 27.6 | 6.5 | 40.7 | 25.1 |
| HG008T ONT INDEL | GRCh38-only | >0.4 | 132 | 3.8 | 6.8 | 33.3 | 56.1 |
| HG008T ONT INDEL | GRCh38-only | unknown | 14 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT INDEL | HPRC v1.1 d9 | 0-0.08 | 14 | 92.9 | 0.0 | 7.1 | 0.0 |
| HG008T ONT INDEL | HPRC v1.1 d9 | 0.08-0.2 | 226 | 69.5 | 7.1 | 17.3 | 6.2 |
| HG008T ONT INDEL | HPRC v1.1 d9 | 0.2-0.4 | 201 | 51.7 | 4.5 | 22.4 | 21.4 |
| HG008T ONT INDEL | HPRC v1.1 d9 | >0.4 | 131 | 28.2 | 0.8 | 12.2 | 58.8 |
| HG008T ONT INDEL | HPRC v1.1 d9 | unknown | 15 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT INDEL | HPRC v2.1 d46 | 0-0.08 | 14 | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T ONT INDEL | HPRC v2.1 d46 | 0.08-0.2 | 225 | 77.8 | 5.3 | 12.0 | 4.9 |
| HG008T ONT INDEL | HPRC v2.1 d46 | 0.2-0.4 | 211 | 48.3 | 3.8 | 21.8 | 26.1 |
| HG008T ONT INDEL | HPRC v2.1 d46 | >0.4 | 122 | 28.7 | 0.0 | 11.5 | 59.8 |
| HG008T ONT INDEL | HPRC v2.1 d46 | unknown | 15 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 0-0.08 | 410 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 0.08-0.2 | 186 | 4.3 | 3.8 | 1.6 | 90.3 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 0.2-0.4 | 102 | 4.9 | 15.7 | 3.9 | 75.5 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | >0.4 | 1745 | 0.8 | 13.1 | 2.3 | 83.8 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 0-0.08 | 10 | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 0.08-0.2 | 8 | 50.0 | 25.0 | 12.5 | 12.5 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 0.2-0.4 | 26 | 50.0 | 7.7 | 15.4 | 26.9 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | >0.4 | 87 | 47.1 | 2.3 | 18.4 | 32.2 |

## 7. VAF sources (chr1; truth records)

| call set (chr1) | arm | truth | tensor | build_min_af | other_arm | truth_vcf | other_platform | unknown |
|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 2443 | 2036 | 276 | 0 | 131 | 0 | 0 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 2443 | 2003 | 291 | 31 | 118 | 0 | 0 |
| COLO829T fiberseq SNV | GRCh38-only | 2443 | 2016 | 302 | 2 | 123 | 0 | 0 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 2443 | 2007 | 321 | 12 | 103 | 0 | 0 |
| COLO829T ONT SNV | GRCh38-only | 2443 | 2009 | 229 | 1 | 204 | 0 | 0 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 2443 | 2001 | 226 | 6 | 210 | 0 | 0 |
| HG008T Illumina SNV | GRCh38-only | 697 | 684 | 2 | 1 | 0 | 9 | 1 |
| HG008T Illumina SNV | HPRC v1.1 d9 | 697 | 660 | 1 | 25 | 0 | 10 | 1 |
| HG008T Illumina SNV | HPRC v2.1 d46 | 697 | 666 | 2 | 19 | 0 | 9 | 1 |
| HG008T PacBio SNV | GRCh38-only | 697 | 683 | 2 | 3 | 0 | 8 | 1 |
| HG008T PacBio SNV | HPRC v1.1 d9 | 697 | 659 | 1 | 27 | 0 | 9 | 1 |
| HG008T PacBio SNV | HPRC v2.1 d46 | 697 | 662 | 1 | 25 | 0 | 8 | 1 |
| HG008T ONT SNV | GRCh38-only | 697 | 681 | 0 | 5 | 0 | 10 | 1 |
| HG008T ONT SNV | HPRC v1.1 d9 | 697 | 664 | 0 | 22 | 0 | 10 | 1 |
| HG008T ONT SNV | HPRC v2.1 d46 | 697 | 666 | 1 | 19 | 0 | 10 | 1 |
| COLO829T Illumina INDEL | GRCh38-only | 131 | 113 | 2 | 0 | 16 | 0 | 0 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 131 | 61 | 1 | 51 | 18 | 0 | 0 |
| COLO829T fiberseq INDEL | GRCh38-only | 131 | 125 | 1 | 0 | 5 | 0 | 0 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 131 | 72 | 6 | 47 | 6 | 0 | 0 |
| COLO829T ONT INDEL | GRCh38-only | 131 | 125 | 5 | 0 | 1 | 0 | 0 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 131 | 71 | 10 | 48 | 2 | 0 | 0 |
| HG008T Illumina INDEL | GRCh38-only | 587 | 506 | 33 | 1 | 0 | 33 | 14 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 587 | 239 | 20 | 256 | 0 | 57 | 15 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 587 | 261 | 15 | 238 | 0 | 58 | 15 |
| HG008T PacBio INDEL | GRCh38-only | 587 | 566 | 7 | 0 | 0 | 3 | 11 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 587 | 293 | 40 | 235 | 0 | 5 | 14 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 587 | 301 | 26 | 241 | 0 | 5 | 14 |
| HG008T ONT INDEL | GRCh38-only | 587 | 500 | 22 | 5 | 0 | 46 | 14 |
| HG008T ONT INDEL | HPRC v1.1 d9 | 587 | 304 | 14 | 199 | 0 | 55 | 15 |
| HG008T ONT INDEL | HPRC v2.1 d46 | 587 | 306 | 14 | 199 | 0 | 53 | 15 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 2443 | 2006 | 328 | 12 | 97 | 0 | 0 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 131 | 67 | 10 | 48 | 6 | 0 | 0 |

## 8. COLO829T funnel by the truth VCF's VAF for every record (chr1; VAF_Ill for Illumina, VAF_PB for fiberseq / ONT; % of the bin)

| call set (chr1) | arm | truth-VCF VAF bin | n | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 0-0.08 | 405 | 91.1 | 2.0 | 0.2 | 6.7 |
| COLO829T Illumina SNV | GRCh38-only | 0.08-0.2 | 186 | 18.3 | 1.6 | 1.6 | 78.5 |
| COLO829T Illumina SNV | GRCh38-only | 0.2-0.4 | 107 | 1.9 | 10.3 | 3.7 | 84.1 |
| COLO829T Illumina SNV | GRCh38-only | >0.4 | 1745 | 0.1 | 31.5 | 1.9 | 66.5 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 0-0.08 | 405 | 93.1 | 0.2 | 0.2 | 6.4 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 0.08-0.2 | 186 | 20.4 | 2.7 | 1.1 | 75.8 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 0.2-0.4 | 107 | 2.8 | 14.0 | 2.8 | 80.4 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | >0.4 | 1745 | 1.3 | 11.3 | 1.9 | 85.5 |
| COLO829T fiberseq SNV | GRCh38-only | 0-0.08 | 411 | 96.4 | 0.0 | 0.2 | 3.4 |
| COLO829T fiberseq SNV | GRCh38-only | 0.08-0.2 | 188 | 10.1 | 1.1 | 2.1 | 86.7 |
| COLO829T fiberseq SNV | GRCh38-only | 0.2-0.4 | 102 | 4.9 | 10.8 | 2.0 | 82.4 |
| COLO829T fiberseq SNV | GRCh38-only | >0.4 | 1742 | 0.4 | 25.0 | 2.1 | 72.6 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 0-0.08 | 411 | 95.9 | 0.5 | 0.2 | 3.4 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 0.08-0.2 | 188 | 11.7 | 3.2 | 1.1 | 84.0 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 0.2-0.4 | 102 | 3.9 | 10.8 | 2.9 | 82.4 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | >0.4 | 1742 | 0.9 | 12.6 | 2.0 | 84.5 |
| COLO829T ONT SNV | GRCh38-only | 0-0.08 | 411 | 90.8 | 1.0 | 0.0 | 8.3 |
| COLO829T ONT SNV | GRCh38-only | 0.08-0.2 | 188 | 28.2 | 1.6 | 1.1 | 69.1 |
| COLO829T ONT SNV | GRCh38-only | 0.2-0.4 | 102 | 2.9 | 9.8 | 2.0 | 85.3 |
| COLO829T ONT SNV | GRCh38-only | >0.4 | 1742 | 0.3 | 36.3 | 1.8 | 61.7 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 0-0.08 | 411 | 91.0 | 1.9 | 0.0 | 7.1 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 0.08-0.2 | 188 | 29.3 | 8.0 | 1.1 | 61.7 |
| COLO829T ONT SNV | HPRC v2.1 d46 | 0.2-0.4 | 102 | 2.9 | 14.7 | 2.9 | 79.4 |
| COLO829T ONT SNV | HPRC v2.1 d46 | >0.4 | 1742 | 0.6 | 13.4 | 2.0 | 84.0 |
| COLO829T Illumina INDEL | GRCh38-only | 0.2-0.4 | 41 | 22.0 | 9.8 | 43.9 | 24.4 |
| COLO829T Illumina INDEL | GRCh38-only | >0.4 | 90 | 11.1 | 30.0 | 37.8 | 21.1 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 0.2-0.4 | 41 | 58.5 | 12.2 | 14.6 | 14.6 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | >0.4 | 90 | 51.1 | 5.6 | 14.4 | 28.9 |
| COLO829T fiberseq INDEL | GRCh38-only | 0.2-0.4 | 18 | 5.6 | 16.7 | 33.3 | 44.4 |
| COLO829T fiberseq INDEL | GRCh38-only | >0.4 | 113 | 6.2 | 18.6 | 47.8 | 27.4 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 0.2-0.4 | 18 | 55.6 | 5.6 | 11.1 | 27.8 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | >0.4 | 113 | 46.9 | 8.0 | 19.5 | 25.7 |
| COLO829T ONT INDEL | GRCh38-only | 0.2-0.4 | 18 | 27.8 | 11.1 | 22.2 | 38.9 |
| COLO829T ONT INDEL | GRCh38-only | >0.4 | 113 | 17.7 | 20.4 | 35.4 | 26.5 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 0.2-0.4 | 18 | 61.1 | 11.1 | 5.6 | 22.2 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | >0.4 | 113 | 56.6 | 3.5 | 16.8 | 23.0 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 0-0.08 | 411 | 96.4 | 0.2 | 0.2 | 3.2 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 0.08-0.2 | 188 | 11.7 | 4.8 | 1.6 | 81.9 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 0.2-0.4 | 102 | 3.9 | 15.7 | 2.9 | 77.5 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | >0.4 | 1742 | 0.9 | 12.9 | 2.4 | 83.9 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 0.2-0.4 | 18 | 55.6 | 0.0 | 11.1 | 33.3 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | >0.4 | 113 | 51.3 | 5.3 | 16.8 | 26.5 |

## 9. Funnel by repeat context (chr1; % of the context's truth)

HP>=7: homopolymer tract >= 7 bp; STR: unit 2-6, >= 3 units; other: the rest (calls_lib.Repeats: pipeline RepeatContext for INS / DEL, its tract() for SNVs).

| call set (chr1) | arm | repeat | n | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | HP>=7 | 31 | 16.1 | 45.2 | 9.7 | 29.0 |
| COLO829T Illumina SNV | GRCh38-only | STR | 127 | 18.9 | 26.0 | 1.6 | 53.5 |
| COLO829T Illumina SNV | GRCh38-only | other | 2285 | 16.5 | 23.0 | 1.6 | 58.9 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | HP>=7 | 31 | 22.6 | 29.0 | 6.5 | 41.9 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | STR | 127 | 22.8 | 18.1 | 0.8 | 58.3 |
| COLO829T Illumina SNV | HPRC v2.1 d46 | other | 2285 | 17.7 | 8.2 | 1.6 | 72.6 |
| COLO829T fiberseq SNV | GRCh38-only | HP>=7 | 31 | 22.6 | 48.4 | 0.0 | 29.0 |
| COLO829T fiberseq SNV | GRCh38-only | STR | 127 | 22.8 | 22.0 | 1.6 | 53.5 |
| COLO829T fiberseq SNV | GRCh38-only | other | 2285 | 17.1 | 17.7 | 1.8 | 63.4 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | HP>=7 | 31 | 25.8 | 32.3 | 6.5 | 35.5 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | STR | 127 | 22.0 | 16.5 | 0.8 | 60.6 |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | other | 2285 | 17.5 | 9.1 | 1.6 | 71.8 |
| COLO829T ONT SNV | GRCh38-only | HP>=7 | 31 | 19.4 | 48.4 | 9.7 | 22.6 |
| COLO829T ONT SNV | GRCh38-only | STR | 127 | 19.7 | 20.5 | 1.6 | 58.3 |
| COLO829T ONT SNV | GRCh38-only | other | 2285 | 17.6 | 26.6 | 1.3 | 54.4 |
| COLO829T ONT SNV | HPRC v2.1 d46 | HP>=7 | 31 | 22.6 | 29.0 | 9.7 | 38.7 |
| COLO829T ONT SNV | HPRC v2.1 d46 | STR | 127 | 20.5 | 17.3 | 0.8 | 61.4 |
| COLO829T ONT SNV | HPRC v2.1 d46 | other | 2285 | 17.9 | 10.5 | 1.6 | 70.0 |
| HG008T Illumina SNV | GRCh38-only | HP>=7 | 34 | 5.9 | 20.6 | 50.0 | 23.5 |
| HG008T Illumina SNV | GRCh38-only | STR | 49 | 6.1 | 30.6 | 10.2 | 53.1 |
| HG008T Illumina SNV | GRCh38-only | other | 614 | 1.3 | 12.9 | 2.1 | 83.7 |
| HG008T Illumina SNV | HPRC v1.1 d9 | HP>=7 | 34 | 44.1 | 29.4 | 5.9 | 20.6 |
| HG008T Illumina SNV | HPRC v1.1 d9 | STR | 49 | 26.5 | 10.2 | 10.2 | 53.1 |
| HG008T Illumina SNV | HPRC v1.1 d9 | other | 614 | 1.5 | 7.2 | 2.4 | 88.9 |
| HG008T Illumina SNV | HPRC v2.1 d46 | HP>=7 | 34 | 32.4 | 38.2 | 8.8 | 20.6 |
| HG008T Illumina SNV | HPRC v2.1 d46 | STR | 49 | 22.4 | 8.2 | 10.2 | 59.2 |
| HG008T Illumina SNV | HPRC v2.1 d46 | other | 614 | 1.5 | 5.0 | 2.6 | 90.9 |
| HG008T PacBio SNV | GRCh38-only | HP>=7 | 34 | 14.7 | 17.6 | 41.2 | 26.5 |
| HG008T PacBio SNV | GRCh38-only | STR | 49 | 8.2 | 10.2 | 18.4 | 63.3 |
| HG008T PacBio SNV | GRCh38-only | other | 614 | 0.8 | 9.3 | 2.6 | 87.3 |
| HG008T PacBio SNV | HPRC v1.1 d9 | HP>=7 | 34 | 47.1 | 29.4 | 2.9 | 20.6 |
| HG008T PacBio SNV | HPRC v1.1 d9 | STR | 49 | 26.5 | 16.3 | 6.1 | 51.0 |
| HG008T PacBio SNV | HPRC v1.1 d9 | other | 614 | 1.5 | 7.2 | 2.6 | 88.8 |
| HG008T PacBio SNV | HPRC v2.1 d46 | HP>=7 | 34 | 38.2 | 26.5 | 11.8 | 23.5 |
| HG008T PacBio SNV | HPRC v2.1 d46 | STR | 49 | 22.4 | 14.3 | 6.1 | 57.1 |
| HG008T PacBio SNV | HPRC v2.1 d46 | other | 614 | 1.8 | 6.4 | 2.9 | 88.9 |
| HG008T ONT SNV | GRCh38-only | HP>=7 | 34 | 20.6 | 17.6 | 38.2 | 23.5 |
| HG008T ONT SNV | GRCh38-only | STR | 49 | 12.2 | 14.3 | 14.3 | 59.2 |
| HG008T ONT SNV | GRCh38-only | other | 614 | 0.5 | 7.5 | 2.4 | 89.6 |
| HG008T ONT SNV | HPRC v1.1 d9 | HP>=7 | 34 | 47.1 | 26.5 | 5.9 | 20.6 |
| HG008T ONT SNV | HPRC v1.1 d9 | STR | 49 | 24.5 | 14.3 | 6.1 | 55.1 |
| HG008T ONT SNV | HPRC v1.1 d9 | other | 614 | 0.8 | 8.5 | 2.4 | 88.3 |
| HG008T ONT SNV | HPRC v2.1 d46 | HP>=7 | 34 | 44.1 | 20.6 | 11.8 | 23.5 |
| HG008T ONT SNV | HPRC v2.1 d46 | STR | 49 | 24.5 | 8.2 | 8.2 | 59.2 |
| HG008T ONT SNV | HPRC v2.1 d46 | other | 614 | 0.7 | 4.6 | 3.3 | 91.5 |
| COLO829T Illumina INDEL | GRCh38-only | HP>=7 | 70 | 5.7 | 27.1 | 45.7 | 21.4 |
| COLO829T Illumina INDEL | GRCh38-only | STR | 27 | 18.5 | 18.5 | 55.6 | 7.4 |
| COLO829T Illumina INDEL | GRCh38-only | other | 34 | 29.4 | 20.6 | 14.7 | 35.3 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | HP>=7 | 70 | 44.3 | 7.1 | 22.9 | 25.7 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | STR | 27 | 77.8 | 11.1 | 7.4 | 3.7 |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | other | 34 | 52.9 | 5.9 | 2.9 | 38.2 |
| COLO829T fiberseq INDEL | GRCh38-only | HP>=7 | 70 | 5.7 | 22.9 | 44.3 | 27.1 |
| COLO829T fiberseq INDEL | GRCh38-only | STR | 27 | 7.4 | 14.8 | 70.4 | 7.4 |
| COLO829T fiberseq INDEL | GRCh38-only | other | 34 | 5.9 | 11.8 | 29.4 | 52.9 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | HP>=7 | 70 | 35.7 | 10.0 | 25.7 | 28.6 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | STR | 27 | 77.8 | 0.0 | 18.5 | 3.7 |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | other | 34 | 50.0 | 8.8 | 2.9 | 38.2 |
| COLO829T ONT INDEL | GRCh38-only | HP>=7 | 70 | 31.4 | 15.7 | 27.1 | 25.7 |
| COLO829T ONT INDEL | GRCh38-only | STR | 27 | 3.7 | 22.2 | 66.7 | 7.4 |
| COLO829T ONT INDEL | GRCh38-only | other | 34 | 5.9 | 23.5 | 20.6 | 50.0 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | HP>=7 | 70 | 52.9 | 4.3 | 18.6 | 24.3 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | STR | 27 | 77.8 | 0.0 | 18.5 | 3.7 |
| COLO829T ONT INDEL | HPRC v2.1 d46 | other | 34 | 50.0 | 8.8 | 5.9 | 35.3 |
| HG008T Illumina INDEL | GRCh38-only | HP>=7 | 435 | 14.3 | 4.8 | 54.3 | 26.7 |
| HG008T Illumina INDEL | GRCh38-only | STR | 65 | 40.0 | 7.7 | 47.7 | 4.6 |
| HG008T Illumina INDEL | GRCh38-only | other | 87 | 29.9 | 16.1 | 13.8 | 40.2 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | HP>=7 | 435 | 57.5 | 3.7 | 15.2 | 23.7 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | STR | 65 | 89.2 | 4.6 | 4.6 | 1.5 |
| HG008T Illumina INDEL | HPRC v1.1 d9 | other | 87 | 49.4 | 6.9 | 2.3 | 41.4 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | HP>=7 | 435 | 53.8 | 5.7 | 17.2 | 23.2 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | STR | 65 | 87.7 | 10.8 | 1.5 | 0.0 |
| HG008T Illumina INDEL | HPRC v2.1 d46 | other | 87 | 48.3 | 3.4 | 2.3 | 46.0 |
| HG008T PacBio INDEL | GRCh38-only | HP>=7 | 435 | 13.8 | 4.1 | 54.7 | 27.4 |
| HG008T PacBio INDEL | GRCh38-only | STR | 65 | 27.7 | 12.3 | 50.8 | 9.2 |
| HG008T PacBio INDEL | GRCh38-only | other | 87 | 16.1 | 29.9 | 21.8 | 32.2 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | HP>=7 | 435 | 48.7 | 4.4 | 22.1 | 24.8 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | STR | 65 | 84.6 | 4.6 | 9.2 | 1.5 |
| HG008T PacBio INDEL | HPRC v1.1 d9 | other | 87 | 49.4 | 2.3 | 1.1 | 47.1 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | HP>=7 | 435 | 47.8 | 6.0 | 21.4 | 24.8 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | STR | 65 | 87.7 | 0.0 | 10.8 | 1.5 |
| HG008T PacBio INDEL | HPRC v2.1 d46 | other | 87 | 46.0 | 3.4 | 2.3 | 48.3 |
| HG008T ONT INDEL | GRCh38-only | HP>=7 | 435 | 46.9 | 5.1 | 27.4 | 20.7 |
| HG008T ONT INDEL | GRCh38-only | STR | 65 | 24.6 | 3.1 | 60.0 | 12.3 |
| HG008T ONT INDEL | GRCh38-only | other | 87 | 25.3 | 5.7 | 20.7 | 48.3 |
| HG008T ONT INDEL | HPRC v1.1 d9 | HP>=7 | 435 | 52.4 | 5.3 | 21.1 | 21.1 |
| HG008T ONT INDEL | HPRC v1.1 d9 | STR | 65 | 83.1 | 1.5 | 12.3 | 3.1 |
| HG008T ONT INDEL | HPRC v1.1 d9 | other | 87 | 50.6 | 2.3 | 1.1 | 46.0 |
| HG008T ONT INDEL | HPRC v2.1 d46 | HP>=7 | 435 | 56.6 | 4.4 | 17.0 | 22.1 |
| HG008T ONT INDEL | HPRC v2.1 d46 | STR | 65 | 83.1 | 0.0 | 15.4 | 1.5 |
| HG008T ONT INDEL | HPRC v2.1 d46 | other | 87 | 47.1 | 1.1 | 3.4 | 48.3 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | HP>=7 | 31 | 25.8 | 22.6 | 12.9 | 38.7 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | STR | 127 | 23.6 | 15.0 | 0.8 | 60.6 |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | other | 2285 | 17.5 | 9.8 | 1.9 | 70.8 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | HP>=7 | 70 | 41.4 | 5.7 | 22.9 | 30.0 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | STR | 27 | 85.2 | 3.7 | 7.4 | 3.7 |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | other | 34 | 47.1 | 2.9 | 8.8 | 41.2 |

## 10. Pooled over the three platforms per sample, kind and arm (in-sample runs; sum of per-run records; % of n)

| sample kind (chr1) | arm |  | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | all | 7329 (100.0 %) | 17.3 | 22.8 | 1.6 | 58.3 |
| COLO829T SNV | HPRC v2.1 d46 | all | 7329 (100.0 %) | 18.0 | 9.9 | 1.6 | 70.4 |
| COLO829T INDEL | GRCh38-only | all | 393 (100.0 %) | 13.2 | 20.4 | 39.7 | 26.7 |
| COLO829T INDEL | HPRC v2.1 d46 | all | 393 (100.0 %) | 52.9 | 6.6 | 16.0 | 24.4 |
| HG008T SNV | GRCh38-only | all | 2091 (100.0 %) | 2.1 | 10.9 | 5.2 | 81.8 |
| HG008T SNV | HPRC v1.1 d9 | all | 2091 (100.0 %) | 5.2 | 9.0 | 3.0 | 82.8 |
| HG008T SNV | HPRC v2.1 d46 | all | 2091 (100.0 %) | 4.6 | 6.8 | 3.7 | 84.9 |
| HG008T INDEL | GRCh38-only | all | 1761 (100.0 %) | 25.4 | 6.9 | 42.3 | 25.4 |
| HG008T INDEL | HPRC v1.1 d9 | all | 1761 (100.0 %) | 56.0 | 4.3 | 15.6 | 24.1 |
| HG008T INDEL | HPRC v2.1 d46 | all | 1761 (100.0 %) | 55.6 | 4.8 | 15.2 | 24.5 |

### 10a. by variant kind

| sample kind (chr1) | arm | kind | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | SNV | 7329 (100.0 %) | 17.3 | 22.8 | 1.6 | 58.3 |
| COLO829T SNV | HPRC v2.1 d46 | SNV | 7329 (100.0 %) | 18.0 | 9.9 | 1.6 | 70.4 |
| COLO829T INDEL | GRCh38-only | INS | 231 (58.8 %) | 11.3 | 22.9 | 43.3 | 22.5 |
| COLO829T INDEL | GRCh38-only | DEL | 162 (41.2 %) | 16.0 | 16.7 | 34.6 | 32.7 |
| COLO829T INDEL | HPRC v2.1 d46 | INS | 231 (58.8 %) | 58.9 | 3.9 | 14.7 | 22.5 |
| COLO829T INDEL | HPRC v2.1 d46 | DEL | 162 (41.2 %) | 44.4 | 10.5 | 17.9 | 27.2 |
| HG008T SNV | GRCh38-only | SNV | 2091 (100.0 %) | 2.1 | 10.9 | 5.2 | 81.8 |
| HG008T SNV | HPRC v1.1 d9 | SNV | 2091 (100.0 %) | 5.2 | 9.0 | 3.0 | 82.8 |
| HG008T SNV | HPRC v2.1 d46 | SNV | 2091 (100.0 %) | 4.6 | 6.8 | 3.7 | 84.9 |
| HG008T INDEL | GRCh38-only | INS | 1104 (62.7 %) | 28.1 | 6.6 | 40.9 | 24.4 |
| HG008T INDEL | GRCh38-only | DEL | 657 (37.3 %) | 21.0 | 7.3 | 44.6 | 27.1 |
| HG008T INDEL | HPRC v1.1 d9 | INS | 1104 (62.7 %) | 62.3 | 2.6 | 13.3 | 21.7 |
| HG008T INDEL | HPRC v1.1 d9 | DEL | 657 (37.3 %) | 45.5 | 7.0 | 19.5 | 28.0 |
| HG008T INDEL | HPRC v2.1 d46 | INS | 1104 (62.7 %) | 60.7 | 3.9 | 13.7 | 21.7 |
| HG008T INDEL | HPRC v2.1 d46 | DEL | 657 (37.3 %) | 47.0 | 6.2 | 17.7 | 29.1 |

### 10b. by truth VAF bin (sources as in section 6)

| sample kind (chr1) | arm | VAF bin | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | 0-0.08 | 1237 (16.9 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T SNV | GRCh38-only | 0.08-0.2 | 601 (8.2 %) | 2.3 | 3.8 | 2.5 | 91.3 |
| COLO829T SNV | GRCh38-only | 0.2-0.4 | 354 (4.8 %) | 1.4 | 12.7 | 2.3 | 83.6 |
| COLO829T SNV | GRCh38-only | >0.4 | 5137 (70.1 %) | 0.2 | 31.2 | 1.9 | 66.7 |
| COLO829T SNV | HPRC v2.1 d46 | 0-0.08 | 1231 (16.8 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T SNV | HPRC v2.1 d46 | 0.08-0.2 | 601 (8.2 %) | 4.8 | 6.3 | 1.8 | 87.0 |
| COLO829T SNV | HPRC v2.1 d46 | 0.2-0.4 | 358 (4.9 %) | 3.4 | 15.6 | 1.7 | 79.3 |
| COLO829T SNV | HPRC v2.1 d46 | >0.4 | 5139 (70.1 %) | 0.9 | 12.4 | 2.0 | 84.8 |
| COLO829T INDEL | GRCh38-only | 0-0.08 | 8 (2.0 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T INDEL | GRCh38-only | 0.08-0.2 | 44 (11.2 %) | 36.4 | 27.3 | 18.2 | 18.2 |
| COLO829T INDEL | GRCh38-only | 0.2-0.4 | 119 (30.3 %) | 10.1 | 10.9 | 48.7 | 30.3 |
| COLO829T INDEL | GRCh38-only | >0.4 | 222 (56.5 %) | 7.2 | 24.8 | 40.5 | 27.5 |
| COLO829T INDEL | HPRC v2.1 d46 | 0-0.08 | 17 (4.3 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T INDEL | HPRC v2.1 d46 | 0.08-0.2 | 59 (15.0 %) | 54.2 | 23.7 | 13.6 | 8.5 |
| COLO829T INDEL | HPRC v2.1 d46 | 0.2-0.4 | 109 (27.7 %) | 53.2 | 3.7 | 20.2 | 22.9 |
| COLO829T INDEL | HPRC v2.1 d46 | >0.4 | 208 (52.9 %) | 48.6 | 3.8 | 15.9 | 31.7 |
| HG008T SNV | GRCh38-only | 0-0.08 | 4 (0.2 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | GRCh38-only | 0.08-0.2 | 38 (1.8 %) | 23.7 | 47.4 | 15.8 | 13.2 |
| HG008T SNV | GRCh38-only | 0.2-0.4 | 245 (11.7 %) | 7.3 | 6.1 | 9.4 | 77.1 |
| HG008T SNV | GRCh38-only | >0.4 | 1801 (86.1 %) | 0.5 | 10.8 | 4.4 | 84.2 |
| HG008T SNV | GRCh38-only | unknown | 3 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v1.1 d9 | 0-0.08 | 2 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v1.1 d9 | 0.08-0.2 | 39 (1.9 %) | 66.7 | 30.8 | 0.0 | 2.6 |
| HG008T SNV | HPRC v1.1 d9 | 0.2-0.4 | 240 (11.5 %) | 13.3 | 10.8 | 3.3 | 72.5 |
| HG008T SNV | HPRC v1.1 d9 | >0.4 | 1807 (86.4 %) | 2.5 | 8.4 | 3.0 | 86.2 |
| HG008T SNV | HPRC v1.1 d9 | unknown | 3 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v2.1 d46 | 0-0.08 | 4 (0.2 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v2.1 d46 | 0.08-0.2 | 27 (1.3 %) | 66.7 | 25.9 | 0.0 | 7.4 |
| HG008T SNV | HPRC v2.1 d46 | 0.2-0.4 | 245 (11.7 %) | 12.7 | 9.0 | 4.1 | 74.3 |
| HG008T SNV | HPRC v2.1 d46 | >0.4 | 1812 (86.7 %) | 2.3 | 6.2 | 3.7 | 87.8 |
| HG008T SNV | HPRC v2.1 d46 | unknown | 3 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | GRCh38-only | 0-0.08 | 62 (3.5 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | GRCh38-only | 0.08-0.2 | 340 (19.3 %) | 55.9 | 4.4 | 31.2 | 8.5 |
| HG008T INDEL | GRCh38-only | 0.2-0.4 | 641 (36.4 %) | 22.0 | 3.9 | 49.9 | 24.2 |
| HG008T INDEL | GRCh38-only | >0.4 | 679 (38.6 %) | 2.4 | 11.9 | 47.0 | 38.7 |
| HG008T INDEL | GRCh38-only | unknown | 39 (2.2 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | HPRC v1.1 d9 | 0-0.08 | 74 (4.2 %) | 98.6 | 0.0 | 1.4 | 0.0 |
| HG008T INDEL | HPRC v1.1 d9 | 0.08-0.2 | 403 (22.9 %) | 71.0 | 10.4 | 14.9 | 3.7 |
| HG008T INDEL | HPRC v1.1 d9 | 0.2-0.4 | 607 (34.5 %) | 55.7 | 3.3 | 19.8 | 21.3 |
| HG008T INDEL | HPRC v1.1 d9 | >0.4 | 633 (35.9 %) | 38.9 | 2.1 | 14.8 | 44.2 |
| HG008T INDEL | HPRC v1.1 d9 | unknown | 44 (2.5 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | HPRC v2.1 d46 | 0-0.08 | 55 (3.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | HPRC v2.1 d46 | 0.08-0.2 | 397 (22.5 %) | 76.8 | 10.1 | 9.8 | 3.3 |
| HG008T INDEL | HPRC v2.1 d46 | 0.2-0.4 | 636 (36.1 %) | 53.8 | 4.9 | 19.7 | 21.7 |
| HG008T INDEL | HPRC v2.1 d46 | >0.4 | 629 (35.7 %) | 37.0 | 2.1 | 16.4 | 44.5 |
| HG008T INDEL | HPRC v2.1 d46 | unknown | 44 (2.5 %) | 100.0 | 0.0 | 0.0 | 0.0 |

### 10c. COLO829T by the truth VCF's VAF (every record)

| sample kind (chr1) | arm | truth-VCF VAF bin | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | 0-0.08 | 1227 (16.7 %) | 92.7 | 1.0 | 0.2 | 6.1 |
| COLO829T SNV | GRCh38-only | 0.08-0.2 | 562 (7.7 %) | 18.9 | 1.4 | 1.6 | 78.1 |
| COLO829T SNV | GRCh38-only | 0.2-0.4 | 311 (4.2 %) | 3.2 | 10.3 | 2.6 | 83.9 |
| COLO829T SNV | GRCh38-only | >0.4 | 5229 (71.3 %) | 0.3 | 30.9 | 1.9 | 66.9 |
| COLO829T SNV | HPRC v2.1 d46 | 0-0.08 | 1227 (16.7 %) | 93.3 | 0.9 | 0.2 | 5.6 |
| COLO829T SNV | HPRC v2.1 d46 | 0.08-0.2 | 562 (7.7 %) | 20.5 | 4.6 | 1.1 | 73.8 |
| COLO829T SNV | HPRC v2.1 d46 | 0.2-0.4 | 311 (4.2 %) | 3.2 | 13.2 | 2.9 | 80.7 |
| COLO829T SNV | HPRC v2.1 d46 | >0.4 | 5229 (71.3 %) | 0.9 | 12.4 | 2.0 | 84.7 |
| COLO829T INDEL | GRCh38-only | 0.2-0.4 | 77 (19.6 %) | 19.5 | 11.7 | 36.4 | 32.5 |
| COLO829T INDEL | GRCh38-only | >0.4 | 316 (80.4 %) | 11.7 | 22.5 | 40.5 | 25.3 |
| COLO829T INDEL | HPRC v2.1 d46 | 0.2-0.4 | 77 (19.6 %) | 58.4 | 10.4 | 11.7 | 19.5 |
| COLO829T INDEL | HPRC v2.1 d46 | >0.4 | 316 (80.4 %) | 51.6 | 5.7 | 17.1 | 25.6 |

### 10d. by repeat context

| sample kind (chr1) | arm | repeat | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | HP>=7 | 93 (1.3 %) | 19.4 | 47.3 | 6.5 | 26.9 |
| COLO829T SNV | GRCh38-only | STR | 381 (5.2 %) | 20.5 | 22.8 | 1.6 | 55.1 |
| COLO829T SNV | GRCh38-only | other | 6855 (93.5 %) | 17.1 | 22.4 | 1.6 | 58.9 |
| COLO829T SNV | HPRC v2.1 d46 | HP>=7 | 93 (1.3 %) | 23.7 | 30.1 | 7.5 | 38.7 |
| COLO829T SNV | HPRC v2.1 d46 | STR | 381 (5.2 %) | 21.8 | 17.3 | 0.8 | 60.1 |
| COLO829T SNV | HPRC v2.1 d46 | other | 6855 (93.5 %) | 17.7 | 9.3 | 1.6 | 71.5 |
| COLO829T INDEL | GRCh38-only | HP>=7 | 210 (53.4 %) | 14.3 | 21.9 | 39.0 | 24.8 |
| COLO829T INDEL | GRCh38-only | STR | 81 (20.6 %) | 9.9 | 18.5 | 64.2 | 7.4 |
| COLO829T INDEL | GRCh38-only | other | 102 (26.0 %) | 13.7 | 18.6 | 21.6 | 46.1 |
| COLO829T INDEL | HPRC v2.1 d46 | HP>=7 | 210 (53.4 %) | 44.3 | 7.1 | 22.4 | 26.2 |
| COLO829T INDEL | HPRC v2.1 d46 | STR | 81 (20.6 %) | 77.8 | 3.7 | 14.8 | 3.7 |
| COLO829T INDEL | HPRC v2.1 d46 | other | 102 (26.0 %) | 51.0 | 7.8 | 3.9 | 37.3 |
| HG008T SNV | GRCh38-only | HP>=7 | 102 (4.9 %) | 13.7 | 18.6 | 43.1 | 24.5 |
| HG008T SNV | GRCh38-only | STR | 147 (7.0 %) | 8.8 | 18.4 | 14.3 | 58.5 |
| HG008T SNV | GRCh38-only | other | 1842 (88.1 %) | 0.9 | 9.9 | 2.4 | 86.9 |
| HG008T SNV | HPRC v1.1 d9 | HP>=7 | 102 (4.9 %) | 46.1 | 28.4 | 4.9 | 20.6 |
| HG008T SNV | HPRC v1.1 d9 | STR | 147 (7.0 %) | 25.9 | 13.6 | 7.5 | 53.1 |
| HG008T SNV | HPRC v1.1 d9 | other | 1842 (88.1 %) | 1.2 | 7.6 | 2.5 | 88.7 |
| HG008T SNV | HPRC v2.1 d46 | HP>=7 | 102 (4.9 %) | 38.2 | 28.4 | 10.8 | 22.5 |
| HG008T SNV | HPRC v2.1 d46 | STR | 147 (7.0 %) | 23.1 | 10.2 | 8.2 | 58.5 |
| HG008T SNV | HPRC v2.1 d46 | other | 1842 (88.1 %) | 1.3 | 5.3 | 2.9 | 90.4 |
| HG008T INDEL | GRCh38-only | HP>=7 | 1305 (74.1 %) | 25.0 | 4.7 | 45.4 | 24.9 |
| HG008T INDEL | GRCh38-only | STR | 195 (11.1 %) | 30.8 | 7.7 | 52.8 | 8.7 |
| HG008T INDEL | GRCh38-only | other | 261 (14.8 %) | 23.8 | 17.2 | 18.8 | 40.2 |
| HG008T INDEL | HPRC v1.1 d9 | HP>=7 | 1305 (74.1 %) | 52.9 | 4.4 | 19.5 | 23.2 |
| HG008T INDEL | HPRC v1.1 d9 | STR | 195 (11.1 %) | 85.6 | 3.6 | 8.7 | 2.1 |
| HG008T INDEL | HPRC v1.1 d9 | other | 261 (14.8 %) | 49.8 | 3.8 | 1.5 | 44.8 |
| HG008T INDEL | HPRC v2.1 d46 | HP>=7 | 1305 (74.1 %) | 52.7 | 5.4 | 18.5 | 23.4 |
| HG008T INDEL | HPRC v2.1 d46 | STR | 195 (11.1 %) | 86.2 | 3.6 | 9.2 | 1.0 |
| HG008T INDEL | HPRC v2.1 d46 | other | 261 (14.8 %) | 47.1 | 2.7 | 2.7 | 47.5 |

### 10e. stage 1 reasons (n (% of all truth n)), PoN tags and the funnel of the untagged truth

| sample kind (chr1) | arm | stage 1 | af_floor | min_variants | absorbed | no_candidate | other_allele | other | stage 1 with a PASS partial-match tensor | stage 1 PoN-tagged | truth PoN-tagged | untagged n | untagged: 1 % | 2 % | 3 % | 4 % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | 1268 (17.3 %) | 807 (11.0 %) | 253 (3.5 %) | 0 (0.0 %) | 208 (2.8 %) | 0 (0.0 %) | 0 (0.0 %) | 0 | 54 / 1268 | 252 / 7329 (3.4 %) | 7077 | 17.2 | 22.5 | 0.0 | 60.4 |
| COLO829T SNV | HPRC v2.1 d46 | 1318 (18.0 %) | 838 (11.4 %) | 308 (4.2 %) | 9 (0.1 %) | 163 (2.2 %) | 0 (0.0 %) | 0 (0.0 %) | 3 | 72 / 1318 | 252 / 7329 (3.4 %) | 7077 | 17.6 | 9.4 | 0.0 | 73.0 |
| COLO829T INDEL | GRCh38-only | 52 (13.2 %) | 8 (2.0 %) | 10 (2.5 %) | 0 (0.0 %) | 12 (3.1 %) | 22 (5.6 %) | 0 (0.0 %) | 10 | 43 / 52 | 243 / 393 (61.8 %) | 150 | 6.0 | 24.0 | 0.0 | 70.0 |
| COLO829T INDEL | HPRC v2.1 d46 | 208 (52.9 %) | 17 (4.3 %) | 8 (2.0 %) | 81 (20.6 %) | 83 (21.1 %) | 19 (4.8 %) | 0 (0.0 %) | 49 | 165 / 208 | 243 / 393 (61.8 %) | 150 | 28.7 | 7.3 | 0.0 | 64.0 |
| HG008T SNV | GRCh38-only | 43 (2.1 %) | 4 (0.2 %) | 16 (0.8 %) | 0 (0.0 %) | 23 (1.1 %) | 0 (0.0 %) | 0 (0.0 %) | 0 | 24 / 43 | 180 / 2091 (8.6 %) | 1911 | 1.0 | 9.5 | 0.0 | 89.5 |
| HG008T SNV | HPRC v1.1 d9 | 108 (5.2 %) | 2 (0.1 %) | 13 (0.6 %) | 35 (1.7 %) | 58 (2.8 %) | 0 (0.0 %) | 0 (0.0 %) | 4 | 70 / 108 | 180 / 2091 (8.6 %) | 1911 | 2.0 | 7.4 | 0.0 | 90.6 |
| HG008T SNV | HPRC v2.1 d46 | 97 (4.6 %) | 4 (0.2 %) | 10 (0.5 %) | 30 (1.4 %) | 53 (2.5 %) | 0 (0.0 %) | 0 (0.0 %) | 3 | 65 / 97 | 180 / 2091 (8.6 %) | 1911 | 1.7 | 5.4 | 0.0 | 92.9 |
| HG008T INDEL | GRCh38-only | 448 (25.4 %) | 62 (3.5 %) | 66 (3.7 %) | 0 (0.0 %) | 61 (3.5 %) | 259 (14.7 %) | 0 (0.0 %) | 104 | 331 / 448 | 1137 / 1761 (64.6 %) | 624 | 18.8 | 9.6 | 0.0 | 71.6 |
| HG008T INDEL | HPRC v1.1 d9 | 987 (56.0 %) | 73 (4.1 %) | 158 (9.0 %) | 275 (15.6 %) | 417 (23.7 %) | 64 (3.6 %) | 0 (0.0 %) | 254 | 815 / 987 | 1137 / 1761 (64.6 %) | 624 | 27.6 | 4.5 | 0.0 | 67.9 |
| HG008T INDEL | HPRC v2.1 d46 | 979 (55.6 %) | 55 (3.1 %) | 124 (7.0 %) | 249 (14.1 %) | 465 (26.4 %) | 86 (4.9 %) | 0 (0.0 %) | 337 | 806 / 979 | 1137 / 1761 (64.6 %) | 624 | 27.7 | 3.2 | 0.0 | 69.1 |

## 11. Invariants (all call sets)

- Stage-3 truth records all PoN-tagged themselves: 32 of 32 call sets.
- No stage-4 (called) truth record is PoN-tagged: 32 of 32 call sets.
- Stage 2 is only `below_threshold` (no PASS exact tensor unplaced or unmatched by rtg): 32 of 32 call sets.
- Stage-1 records whose recall status says `tensor_representative` but without an exact tensor in the predictions: 0.
- Stage 3 / 4 records without an exact tensor (rtg matched another tensor's call), per call set: {'HG008_ONT_INDEL_7fc7396__vcf_chr1': (2, 0), 'HG008_ONT_INDEL_d46_7fc7396__vcf_chr1': (1, 0)}.
- Stage-3 records whose exact tensor's call is not PoN-tagged (the match came from another call): {'HG008_ONT_INDEL_7fc7396__vcf_chr1': 2, 'HG008_ONT_INDEL_d46_7fc7396__vcf_chr1': 1}.

All 32 call sets reproduce rtg exactly: True. Runtime 30 s.
