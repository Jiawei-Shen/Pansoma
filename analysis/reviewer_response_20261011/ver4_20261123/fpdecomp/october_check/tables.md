# False-positive decomposition, six-tumor current models (ver4, 2026-11-23)

chr1, benchmark regions, PASS calls, seed 0. GRCh38-only = six-tumor model without population site features as a model input (m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef); HPRC v2.1 = six-tumor current model (SNV snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, INDEL release_s0_<P>_INDEL). Share = % of the GRCh38-graph FP alleles; candidates retained = % of the alleles of that class with a tensor in the v2.1 test set of the same reads. In the HPRC v2.1 graph = carried by >= 46 of 463 HPRC v2.1 haplotypes. A class with fewer than 20 alleles shows "k of n" instead of the retained percentage; "-" = the class is empty; "<0.1%" = above zero and under 0.05 %.

## Response Table 3.8 (= 2.6), SNV, FP before population filtering

| Dataset | GRCh38-graph FP before population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP before population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 121,139 | 84%; 0.3% | 13%; 85% | 3% | 2,870 |

## Response Table 3.9, INDEL, FP before population filtering

| Dataset | GRCh38-graph FP before population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP before population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 39,930 | 79%; 0.3% | 16%; 44% | 5% | 2,867 |

## Second pair: SNV, FP after population filtering

| Dataset | GRCh38-graph FP after population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP after population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 1,089 | 2%; 1 of 17 | 57%; 93% | 41% | 792 |

## Second pair: INDEL, FP after population filtering

| Dataset | GRCh38-graph FP after population filtering | Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1) | Non-germline FP (% of GRCh38-graph FP) | HPRC v2.1 FP after population filtering |
|---|---|---|---|---|---|
| HG008T Illumina | 1,348 | 5%; 0% | 26%; 56% | 69% | 385 |

## Counts behind the four tables

n = GRCh38-only FP alleles of the class; candidate on v2.1 = of them, alleles with a tensor in the v2.1 set; FP on v2.1 = of them, alleles that are also a PASS false positive of the v2.1 model at the same stage.

### SNV, FP before population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 697 | 121,139 (121,139) | 101,866; 322; 5 | 16,200; 13,797; 2,057 | 3,073; 713; 71 | 2,870 (2,870) | 2,133 |

### INDEL, FP before population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 587 | 39,930 (39,930) | 31,450; 87; 40 | 6,316; 2,804; 1,887 | 2,164; 98; 26 | 2,867 (2,867) | 1,953 |

### SNV, FP after population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 697 | 1,089 (1,089) | 17; 1; 0 | 625; 584; 516 | 447; 229; 11 | 792 (792) | 527 |

### INDEL, FP after population filtering

| Dataset | truth (TP + FN) | GRCh38-only FP records (allele keys) | germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1 | germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1 | not germline: n; candidate on v2.1; FP on v2.1 | HPRC v2.1 FP records (allele keys) | v2.1 FP that are also GRCh38-only FP |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 587 | 1,348 (1,348) | 73; 0; 0 | 346; 193; 146 | 929; 54; 13 | 385 (385) | 159 |

## PASS recall and precision of the two models (report.json)

The two models are not at the same recall: the FP totals of the first and last column go with these recalls.

| kind | Dataset | truth | GRCh38-only before: TP / FP / recall / precision | HPRC v2.1 before: TP / FP / recall / precision | GRCh38-only after: TP / FP / recall / precision | HPRC v2.1 after: TP / FP / recall / precision |
|---|---|---|---|---|---|---|
| SNV | HG008T Illumina | 697 | 583 / 121,139 / 0.836 / 0.005 | 618 / 2,870 / 0.887 / 0.177 | 548 / 1,089 / 0.786 / 0.335 | 594 / 792 / 0.852 / 0.429 |
| INDEL | HG008T Illumina | 587 | 433 / 39,930 / 0.738 / 0.011 | 219 / 2,867 / 0.373 / 0.071 | 152 / 1,348 / 0.259 / 0.101 | 141 / 385 / 0.240 / 0.268 |

## Sets and PASS rules

| Dataset | kind | GRCh38-only set | PASS rule (GRCh38-only) | HPRC v2.1 set | PASS rule (v2.1) |
|---|---|---|---|---|---|
| HG008T Illumina | SNV | HG008T_Illumina_linear | vcf_chr1: p_somatic >= 0.06399419158697128 (validation recall 0.9) | HG008T_Illumina_d46 | vcf_chr1: p_somatic >= 0.32374927401542664 (validation recall 0.9); off-reference tensors rescored by predict --offref-site90: offref_call (argmax) |
| HG008T Illumina | INDEL | HG008T_Illumina_linear | vcf_chr1: p_somatic >= 0.04801846295595169 (validation recall 0.9) | HG008T_Illumina_d46 | vcf_chr1: p_somatic >= 0.10646037012338638 (validation recall 0.9) |
