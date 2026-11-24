# ver4 offref: checks (checks.py)

## 1. Truth counts (TP + FN) per call set

| call set | in BED: baseline, raw TP+FN, after-filter TP+FN | no BED: baseline, raw TP+FN, after-filter TP+FN | E18 truth spectrum n |
|---|---|---|---|
| HG008T Illumina | 697, 697, 697 | 702, 702, 702 | 697 |
| HG008T PacBio HiFi | 697, 697, 697 | 702, 702, 702 | 697 |
| HG008T ONT | 697, 697, 697 | 702, 702, 702 | 697 |
| COLO829T Illumina | 2443, 2443, 2443 | 2443, 2443, 2443 | 2443 |
| COLO829T Fiber-seq | 2443, 2443, 2443 | 2443, 2443, 2443 | 2443 |
| COLO829T ONT | 2443, 2443, 2443 | 2443, 2443, 2443 | 2443 |
| HCC1395T Illumina | 3440, 3440, 3440 | 3440, 3440, 3440 | 3440 |
| HCC1395T PacBio HiFi | 3440, 3440, 3440 | 3440, 3440, 3440 | 3440 |
| HCC1395T ONT | 3440, 3440, 3440 | 3440, 3440, 3440 | 3440 |

All equal within a tumor (three platforms) and equal to the E18 truth spectrum counts: OK.

## 2. Graph-coordinate PASS call counts per source

| call set | unplaced VCF PASS | report.json (BED / no BED) | E9c run.json | per_call rows | E9i rows | E9a projection PASS rows | expected |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 77 | 77 / 77 | 77 | 77 | 77 | 77 | 77 |
| HG008T PacBio HiFi | 143 | 143 / 143 | 143 | 143 | 143 | 143 | 143 |
| HG008T ONT | 90 | 90 / 90 | 90 | 90 | 90 | 90 | 90 |
| COLO829T Illumina | 110 | 110 / 110 | 110 | 110 | 110 | 110 | 110 |
| COLO829T Fiber-seq | 79 | 79 / 79 | 79 | 79 | 79 | 79 | 79 |
| COLO829T ONT | 111 | 111 / 111 | 111 | 111 | 111 | 111 | 111 |
| HCC1395T Illumina | 76 | 76 / 76 | 76 | 76 | 76 | 76 | 76 |
| HCC1395T PacBio HiFi | 106 | 106 / 106 | 106 | 106 | 106 | 106 | 106 |
| HCC1395T ONT | 143 | 143 / 143 | 143 | 143 | 143 | 143 | 143 |

All equal: OK. Total 935 calls.

## 3. E9i read counts against the ver3 runs (calls counted in both, same candidate key and 21-mers)

| tumor | calls in both (any group) | identical counts and status | different | of them graph-coordinate PASS calls: in both / identical |
|---|---|---|---|---|
| HG008T | 422 | 422 | 0 | 23 / 23 |
| COLO829T | 281 | 281 | 0 | 24 / 24 |
| HCC1395T | 362 | 362 | 0 | 181 / 181 |

## 4. HG008-N-D remote fetch

regions.bed: 822 regions; fetch log: `done 822 regions 85756 reads`; retries logged: 0; index present: True.

## 5. Call sets (from the VCF headers and report.json)

| call set | tensor set | checkpoint | PASS threshold | call rule | records rescored / called by the off-reference rescue | OFFREF_P_SOMATIC records in the graph-coordinate VCF |
|---|---|---|---|---|---|---|
| HG008T Illumina | HG008T_Illumina_d46 | snv_newpipe_20261109/multicancer6_Illumina_SNV_newpipe/best.pth | 0.7630 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| HG008T PacBio HiFi | HG008T_PacBio_d46 | snv_newpipe_20261109/multicancer6_PacBio_SNV_newpipe/best.pth | 0.8157 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| HG008T ONT | HG008T_ONT_d46 | snv_newpipe_20261109/multicancer6_ONT_SNV_newpipe/best.pth | 0.8221 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| COLO829T Illumina | COLO829T_Illumina_d46_snvaf05 | snv_newpipe_20261109/multicancer6_Illumina_SNV_newpipe/best.pth | 0.7630 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| COLO829T Fiber-seq | COLO829T_fiberseq_d46_snvaf05 | snv_newpipe_20261109/multicancer6_PacBio_SNV_newpipe/best.pth | 0.8157 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| COLO829T ONT | COLO829T_ONT_d46 | snv_newpipe_20261109/multicancer6_ONT_SNV_newpipe/best.pth | 0.8221 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| HCC1395T Illumina | HCC1395T_Illumina_d46_snvaf05 | snv_newpipe_20261109/multicancer6_Illumina_SNV_newpipe/best.pth | 0.7630 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| HCC1395T PacBio HiFi | HCC1395T_PacBio_d46_snvaf05 | snv_newpipe_20261109/multicancer6_PacBio_SNV_newpipe/best.pth | 0.8157 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |
| HCC1395T ONT | HCC1395T_ONT_d46 | snv_newpipe_20261109/multicancer6_ONT_SNV_newpipe/best.pth | 0.8221 | validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable) | 0 / 0 | 0 |

All nine: HPRC v2.1 (d46) graph index, kind SNV, chr1, no off-reference rescue: OK.

## 6. E9c + E9d reproduced on a ver3 run with this batch's code (HCC1395T Illumina six-tumor no-rescue run of 2026-11-04)

E9c calls.tsv (13 graph-coordinate PASS calls + controls + gold; 465 rows): byte-identical to the 2026-11-04 file.
E9d (queries, minimap2 on the five HCC1395 assemblies, classify, collect): 465 rows in the old per_call.tsv.gz, 465 in the reproduction, 465 with the same assembly class, flag, category and population tag, 0 different; graph-coordinate PASS classes in the reproduction {'germline': 9, 'tumor_only': 3, 'unsupported': 1} (ver3 letter row: 9 normal / 3 tumor-only / 1 neither).
