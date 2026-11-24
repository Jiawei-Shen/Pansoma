# Response Table 3.10 (ver4): validation of graph-coordinate SNV PASS calls

Six-tumor current SNV models (snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, seed 0), HPRC v2.1, chr1, default output. Every call was evaluated. Definitions in table_3_10.py.

| **Dataset (chr1 SNV)** | **Graph-coordinate PASS calls** | **In a normal assembly** | **Tumor-only** | **Neither** | **ALT in ≥2 normal reads** | **No normal ALT read, ALT in tumor reads** | **Tumor-only and no normal ALT read** | **Also PASS on both other platforms** |
|---|---|---|---|---|---|---|---|---|
| HG008T Illumina | 77 | 71 (92.2%) | 3 (3.9%) | 3 (3.9%) | 64 (83.1%) | 2 (2.6%) | 1 | 3 |
| HG008T PacBio HiFi | 143 | 127 (88.8%) | 0 (0.0%) | 16 (11.2%) | 136 (95.1%) | 1 (0.7%) | 0 | 3 |
| HG008T ONT | 90 | 90 (100.0%) | 0 (0.0%) | 0 (0.0%) | 82 (91.1%) | 0 (0.0%) | 0 | 3 |
| COLO829T Illumina | 110 | 98 (89.1%) | 8 (7.3%) | 4 (3.6%) | 81 (73.6%) | 6 (5.5%) | 4 | 5 |
| COLO829T Fiber-seq | 79 | 73 (92.4%) | 5 (6.3%) | 1 (1.3%) | 63 (79.7%) | 7 (8.9%) | 5 | 5 |
| COLO829T ONT | 111 | 94 (84.7%) | 15 (13.5%) | 2 (1.8%) | 98 (88.3%) | 12 (10.8%) | 11 | 5 |
| HCC1395T Illumina | 76 | 71 (93.4%) | 5 (6.6%) | 0 (0.0%) | 55 (72.4%) | 3 (3.9%) | 3 | 7 |
| HCC1395T PacBio HiFi | 106 | 102 (96.2%) | 3 (2.8%) | 1 (0.9%) | 97 (91.5%) | 2 (1.9%) | 2 | 7 |
| HCC1395T ONT | 143 | 130 (90.9%) | 1 (0.7%) | 12 (8.4%) | 139 (97.2%) | 1 (0.7%) | 1 | 7 |

Totals, nine call sets: 935 calls; in a normal assembly 856 (91.6%), tumor-only 40 (4.3%), neither 39 (4.2%); ALT in >= 2 normal reads 815 (87.2%); no normal ALT read and ALT in tumor reads 34 (3.6%); neither read class 86 (9.2%); tumor-only assembly and no normal ALT read 27 (2.9%).

Per tumor (three platforms summed):

| Tumor | Calls | In a normal assembly | Tumor-only | Neither | ALT in ≥2 normal reads | No normal ALT read, ALT in tumor reads | Tumor-only and no normal ALT read | Keys PASS on all three platforms |
|---|---|---|---|---|---|---|---|---|
| HG008T | 310 | 288 (92.9%) | 3 (1.0%) | 19 (6.1%) | 282 (91.0%) | 3 (1.0%) | 1 | 3 |
| COLO829T | 300 | 265 (88.3%) | 28 (9.3%) | 7 (2.3%) | 242 (80.7%) | 25 (8.3%) | 20 | 5 |
| HCC1395T | 325 | 303 (93.2%) | 9 (2.8%) | 13 (4.0%) | 291 (89.5%) | 6 (1.8%) | 6 | 7 |

Keys PASS on all three platforms of a tumor (assembly class; read class per platform):

| Tumor | Key (node:offset:SNP:REF>ALT) | Assembly | Reads Illumina / PacBio HiFi or Fiber-seq / ONT |
|---|---|---|---|
| HG008T | 1280181:0:SNP:C>T | normal | alt_normal / alt_normal / alt_normal |
| HG008T | 1519960:39:SNP:C>A | normal | alt_normal / alt_normal / alt_normal |
| HG008T | 5736980:22:SNP:C>T | normal | other / other / other |
| HG008T | (3 keys on all three; 24 on two or more; 283 distinct keys) | | |
| COLO829T | 10588020:0:SNP:G>A | normal | alt_normal / alt_normal / alt_normal |
| COLO829T | 10793319:10:SNP:G>A | tumor | somatic / somatic / somatic |
| COLO829T | 1249536:1:SNP:C>A | normal | alt_normal / alt_normal / alt_normal |
| COLO829T | 12743663:20:SNP:C>A | normal | alt_normal / alt_normal / alt_normal |
| COLO829T | 1801552:3:SNP:G>T | tumor | somatic / somatic / somatic |
| COLO829T | (5 keys on all three; 30 on two or more; 265 distinct keys) | | |
| HCC1395T | 11567240:0:SNP:G>A | normal | alt_normal / alt_normal / alt_normal |
| HCC1395T | 12744951:9:SNP:C>A | normal | alt_normal / alt_normal / alt_normal |
| HCC1395T | 1519960:39:SNP:C>A | normal | alt_normal / alt_normal / alt_normal |
| HCC1395T | 19532684:0:SNP:C>T | normal | alt_normal / alt_normal / alt_normal |
| HCC1395T | 4379492:231:SNP:G>A | tumor | somatic / somatic / somatic |
| HCC1395T | 6556831:0:SNP:C>G | normal | other / other / other |
| HCC1395T | 6557052:124:SNP:T>C | normal | other / other / other |
| HCC1395T | (7 keys on all three; 34 on two or more; 284 distinct keys) | | |

'Neither' by the E9d class (unsupported = only the REF context is found in an assembly; unverifiable = neither context):
HG008T Illumina {'unverifiable': 2, 'unsupported': 1}; HG008T PacBio HiFi {'unsupported': 15, 'unverifiable': 1}; HG008T ONT {}; COLO829T Illumina {'unsupported': 3, 'unverifiable': 1}; COLO829T Fiber-seq {'unverifiable': 1}; COLO829T ONT {'unsupported': 2}; HCC1395T Illumina {}; HCC1395T PacBio HiFi {'unsupported': 1}; HCC1395T ONT {'unverifiable': 3, 'unsupported': 9}.

Read class 'other' (neither read column) by E9i status:
HG008T Illumina {'low normal coverage': 7, 'absent from normal; tumor ALT < 2': 3, 'no k-mer': 1}; HG008T PacBio HiFi {'low normal coverage': 4, 'absent from normal; tumor ALT < 2': 2}; HG008T ONT {'low normal coverage': 5, 'one normal read': 1, 'absent from normal; tumor ALT < 2': 2}; COLO829T Illumina {'absent from normal; tumor ALT < 2': 5, 'no k-mer': 2, 'one normal read': 13, 'low normal coverage': 3}; COLO829T Fiber-seq {'one normal read': 8, 'absent from normal; tumor ALT < 2': 1}; COLO829T ONT {'one normal read': 1}; HCC1395T Illumina {'no k-mer': 2, 'low normal coverage': 12, 'absent from normal; tumor ALT < 2': 4}; HCC1395T PacBio HiFi {'low normal coverage': 6, 'one normal read': 1}; HCC1395T ONT {'low normal coverage': 3}.

Controls (placed PASS calls of the same call sets, up to 150 each; ctrl TP = rtg true positives, ctrl FP = rtg false positives after the population filter, whole chr1): expected tumor-only / no normal ALT read for TP, normal for FP.

| Dataset | ctrl TP: tumor-only assembly | ctrl TP: in a normal assembly | ctrl TP: no normal ALT read, ALT in tumor reads | ctrl TP: ALT in ≥2 normal reads | ctrl FP: in a normal assembly | ctrl FP: tumor-only assembly | ctrl FP: ALT in ≥2 normal reads | ctrl FP: no normal ALT read, ALT in tumor reads |
|---|---|---|---|---|---|---|---|---|
| HG008T Illumina | 147/150 | 2/150 | 150/150 | 0/150 | 137/150 | 1/150 | 122/150 | 9/150 |
| HG008T PacBio HiFi | 148/150 | 2/150 | 150/150 | 0/150 | 138/150 | 0/150 | 137/150 | 11/150 |
| HG008T ONT | 149/150 | 1/150 | 150/150 | 0/150 | 131/150 | 4/150 | 127/150 | 16/150 |
| COLO829T Illumina | 128/150 | 2/150 | 133/150 | 3/150 | 120/150 | 15/150 | 114/150 | 21/150 |
| COLO829T Fiber-seq | 120/150 | 0/150 | 130/150 | 7/150 | 111/150 | 24/150 | 109/150 | 31/150 |
| COLO829T ONT | 129/150 | 3/150 | 138/150 | 3/150 | 121/150 | 8/150 | 97/150 | 26/150 |
| HCC1395T Illumina | 111/150 | 0/150 | 150/150 | 0/150 | 101/150 | 18/150 | 85/150 | 33/150 |
| HCC1395T PacBio HiFi | 111/150 | 0/150 | 149/150 | 1/150 | 83/150 | 32/150 | 80/150 | 70/150 |
| HCC1395T ONT | 130/150 | 0/150 | 150/150 | 0/150 | 116/150 | 22/150 | 50/150 | 35/150 |
