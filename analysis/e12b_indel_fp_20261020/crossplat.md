# HCC1395T chr1: multicancer4's HP>=7 / STR false INDEL calls on the other platforms

(`crossplat.py`. Exact = same POS / REF / ALT; near = same kind (DEL / INS) within 5 bp. Pansoma = the other platform's multicancer4 PASS calls after the PoN (TP + FP); DS / CT = DeepSomatic-TO / ClairS-TO records of any FILTER, RefCall excluded.)

| platform | HP>=7 / STR false calls | in Pansoma of the other platforms: exact / near | in DS of the other platforms: exact / near | in CT of the other platforms: exact / near | in any of them, near |
|---|---:|---|---|---|---|
| Illumina | 58 | 34 (58.6%) / 34 (58.6%) | 35 (60.3%) / 38 (65.5%) | 29 (50.0%) / 30 (51.7%) | 43 (74.1%) |
| PacBio | 288 | 54 (18.8%) / 54 (18.8%) | 49 (17.0%) / 51 (17.7%) | 48 (16.7%) / 49 (17.0%) | 71 (24.7%) |
| ONT | 13,061 | 55 (0.4%) / 57 (0.4%) | 51 (0.4%) / 78 (0.6%) | 51 (0.4%) / 75 (0.6%) | 97 (0.7%) |
