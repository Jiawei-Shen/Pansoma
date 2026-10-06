# c7: HPRC membership and cross-evaluation of the graph-absorbed COLO829T somatic truth alleles (d9 HPRC v1.1)
RARE_AF = 0.2 on the exact-allele HPRC frequency (complete haplotypes spelling the ALT allele over the tandem array). COLO829BL calls: c3 whole-array alleles. No read data.

### Categories per perfect set, chr1-22

| category | fiberseq SNV chr1-22 | fiberseq INDEL chr1-22 | ONT SNV chr1-22 | ONT INDEL chr1-22 | Illumina SNV chr1-22 | Illumina INDEL chr1-22 | union SNV chr1-22 | union INDEL chr1-22 |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | 0 | 19 | 0 | 18 | 0 | 20 | 0 | 25 |
| normal_present_rare | 1 | 8 | 1 | 4 | 1 | 6 | 1 | 8 |
| normal_absent_HPRC_other | 60 | 176 | 55 | 130 | 50 | 347 | 63 | 367 |
| other_ambiguous | 3 | 15 | 2 | 18 | 2 | 23 | 3 | 30 |
| total | 64 | 218 | 58 | 170 | 53 | 396 | 67 | 430 |

### Categories per perfect set, chr1

| category | fiberseq SNV chr1 | fiberseq INDEL chr1 | ONT SNV chr1 | ONT INDEL chr1 | Illumina SNV chr1 | Illumina INDEL chr1 | union SNV chr1 | union INDEL chr1 |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | 0 | 8 | 0 | 8 | 0 | 9 | 0 | 11 |
| normal_present_rare | 0 | 3 | 0 | 1 | 0 | 2 | 0 | 3 |
| normal_absent_HPRC_other | 4 | 12 | 4 | 9 | 3 | 22 | 4 | 22 |
| other_ambiguous | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 2 |
| total | 5 | 23 | 4 | 19 | 4 | 34 | 5 | 38 |

### Categories per perfect set, chr2-22

| category | fiberseq SNV chr2-22 | fiberseq INDEL chr2-22 | ONT SNV chr2-22 | ONT INDEL chr2-22 | Illumina SNV chr2-22 | Illumina INDEL chr2-22 | union SNV chr2-22 | union INDEL chr2-22 |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | 0 | 11 | 0 | 10 | 0 | 11 | 0 | 14 |
| normal_present_rare | 1 | 5 | 1 | 3 | 1 | 4 | 1 | 5 |
| normal_absent_HPRC_other | 56 | 164 | 51 | 121 | 47 | 325 | 59 | 345 |
| other_ambiguous | 2 | 15 | 2 | 17 | 1 | 22 | 2 | 28 |
| total | 59 | 195 | 54 | 151 | 49 | 362 | 62 | 392 |

### Categories with sub-class / sub-reason, chr1-22

| category | sub_class / sub_reason | fiberseq SNV chr1-22 | fiberseq INDEL chr1-22 | ONT SNV chr1-22 | ONT INDEL chr1-22 | Illumina SNV chr1-22 | Illumina INDEL chr1-22 | union SNV chr1-22 | union INDEL chr1-22 |
|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | hapX_only | 0 | 18 | 0 | 17 | 0 | 19 | 0 | 23 |
| normal_present_broad | hapX_only:path_allele | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| normal_present_broad | hapY_only | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| normal_present_rare | hapX_only | 1 | 7 | 1 | 4 | 1 | 5 | 1 | 7 |
| normal_present_rare | hapY_only | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| normal_absent_HPRC_other | element_set_other_allele | 5 | 24 | 3 | 20 | 0 | 27 | 5 | 33 |
| normal_absent_HPRC_other | exact_allele | 55 | 150 | 52 | 108 | 50 | 318 | 58 | 332 |
| normal_absent_HPRC_other | recombinant_pieces | 0 | 2 | 0 | 2 | 0 | 2 | 0 | 2 |
| other_ambiguous | closest:graph_allele_nearer_truth | 1 | 2 | 1 | 3 | 0 | 5 | 1 | 5 |
| other_ambiguous | closest:graph_allele_not_truth | 2 | 11 | 1 | 13 | 2 | 11 | 2 | 16 |
| other_ambiguous | no_HPRC_carrier | 0 | 2 | 0 | 1 | 0 | 3 | 0 | 4 |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 2 |
| other_ambiguous | normal_hap_unresolved:both_NA | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 |

### Categories with sub-class / sub-reason, chr1

| category | sub_class / sub_reason | fiberseq SNV chr1 | fiberseq INDEL chr1 | ONT SNV chr1 | ONT INDEL chr1 | Illumina SNV chr1 | Illumina INDEL chr1 | union SNV chr1 | union INDEL chr1 |
|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | hapX_only | 0 | 7 | 0 | 7 | 0 | 9 | 0 | 10 |
| normal_present_broad | hapX_only:path_allele | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| normal_present_broad | hapY_only | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| normal_present_rare | hapX_only | 0 | 3 | 0 | 1 | 0 | 2 | 0 | 3 |
| normal_present_rare | hapY_only | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| normal_absent_HPRC_other | element_set_other_allele | 1 | 3 | 1 | 2 | 0 | 3 | 1 | 3 |
| normal_absent_HPRC_other | exact_allele | 3 | 9 | 3 | 7 | 3 | 19 | 3 | 19 |
| normal_absent_HPRC_other | recombinant_pieces | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | closest:graph_allele_nearer_truth | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | closest:graph_allele_not_truth | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 |
| other_ambiguous | no_HPRC_carrier | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| other_ambiguous | normal_hap_unresolved:both_NA | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Interpretation per category and sub-class (union, chr1-22)

| category | sub_class / sub_reason | SNV | INDEL | interpretation |
|---|---|---|---|---|
| normal_present_broad | hapX_only | 0 | 23 | COLO829BL germline allele already a graph path (common exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapX); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_present_broad | hapX_only:path_allele | 0 | 1 | COLO829BL germline allele already a graph path (common exact allele in HPRC): the d9 ALT path (GRCh38 + truth + the germline alleles it needs) spells over the whole tandem array the COLO829BL hapX allele; the perfectly aligned reads carry the patient's own germline allele: germline filtering by the graph (patient frame) |
| normal_present_broad | hapY_only | 0 | 1 | COLO829BL germline allele already a graph path (common exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapY); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_present_rare | hapX_only | 1 | 7 | COLO829BL germline allele already a graph path (rare exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapX); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_present_rare | hapY_only | 0 | 1 | COLO829BL germline allele already a graph path (rare exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapY); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_absent_HPRC_other | element_set_other_allele | 5 | 33 | the patient does not carry it; HPRC haplotypes walk all somatic elements but spell another allele over the array (e.g. another repeat length): the graph spells the somatic allele from population nodes, not a population allele: pangenome-induced false negative |
| normal_absent_HPRC_other | exact_allele | 58 | 332 | the patient does not carry it; >= 1 complete HPRC haplotype spells the exact ALT allele over the array (somatic slippage / substitution recreating a population allele kept by the d9 graph): pangenome-induced false negative |
| normal_absent_HPRC_other | recombinant_pieces | 0 | 2 | the patient does not carry it; every somatic element is walked by some HPRC haplotype but none walks a whole ALT path: the ALT path is a recombination of population nodes: pangenome-induced false negative |
| other_ambiguous | closest:graph_allele_nearer_truth | 1 | 5 | unresolved: no d9 path spells the truth; the closest graph allele is nearer the truth, but not the truth (no read data) |
| other_ambiguous | closest:graph_allele_not_truth | 2 | 16 | unresolved: no d9 path spells the truth; the closest graph allele is nearer GRCh38 / the patient's germline than the truth (no read data to say what the reads walk) |
| other_ambiguous | no_HPRC_carrier | 0 | 4 | unresolved: no complete HPRC haplotype walks the somatic elements |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 2 | unresolved: no complete HPRC haplotype walks the somatic elements; CHM13 does |
| other_ambiguous | normal_hap_unresolved:both_NA | 0 | 1 | unresolved: no COLO829BL hap has an array sequence |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 2 | unresolved: the ALT path needs a germline allele that dipcall phases to the other hap than the c3 event hap |

### HPRC frequency per category, exact allele (spells the ALT allele over the array; RARE_AF level) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 25 | 0 | 0 | 0 | 0 | 19 | 6 | 0 | 0.398 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.125 |
| normal_present_rare | INDEL | 8 | 1 | 1 | 1 | 5 | 0 | 0 | 0 | 0.152 |
| normal_absent_HPRC_other | SNV | 63 | 5 | 6 | 1 | 21 | 24 | 5 | 1 | 0.188 |
| normal_absent_HPRC_other | INDEL | 367 | 35 | 31 | 28 | 138 | 112 | 23 | 0 | 0.159 |
| other_ambiguous | SNV | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 |
| other_ambiguous | INDEL | 27 | 24 | 0 | 0 | 1 | 2 | 0 | 0 | 0.0 |

### HPRC frequency per category, exact allele over all 88 HPRC haplotypes (lower bound) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 25 | 0 | 0 | 0 | 6 | 15 | 4 | 0 | 0.33 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.114 |
| normal_present_rare | INDEL | 8 | 1 | 2 | 0 | 5 | 0 | 0 | 0 | 0.108 |
| normal_absent_HPRC_other | SNV | 63 | 5 | 6 | 1 | 23 | 22 | 5 | 1 | 0.182 |
| normal_absent_HPRC_other | INDEL | 367 | 35 | 41 | 29 | 152 | 95 | 15 | 0 | 0.136 |
| other_ambiguous | SNV | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 |
| other_ambiguous | INDEL | 27 | 24 | 0 | 0 | 2 | 1 | 0 | 0 | 0.0 |

### HPRC frequency per category, set level, any ALT path (walks all somatic elements) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 25 | 0 | 0 | 0 | 0 | 19 | 6 | 0 | 0.398 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.125 |
| normal_present_rare | INDEL | 8 | 0 | 1 | 1 | 5 | 0 | 1 | 0 | 0.172 |
| normal_absent_HPRC_other | SNV | 63 | 0 | 4 | 1 | 25 | 25 | 6 | 2 | 0.205 |
| normal_absent_HPRC_other | INDEL | 367 | 2 | 17 | 25 | 143 | 148 | 29 | 3 | 0.193 |
| other_ambiguous | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.159 |
| other_ambiguous | INDEL | 27 | 9 | 3 | 0 | 5 | 6 | 4 | 0 | 0.114 |

### HPRC frequency per category, element level (min over elements) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 25 | 0 | 0 | 0 | 1 | 18 | 6 | 0 | 0.398 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.125 |
| normal_present_rare | INDEL | 8 | 0 | 0 | 1 | 6 | 0 | 1 | 0 | 0.172 |
| normal_absent_HPRC_other | SNV | 63 | 0 | 4 | 1 | 25 | 25 | 6 | 2 | 0.205 |
| normal_absent_HPRC_other | INDEL | 367 | 1 | 16 | 18 | 149 | 153 | 28 | 2 | 0.198 |
| other_ambiguous | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.159 |
| other_ambiguous | INDEL | 27 | 7 | 3 | 0 | 6 | 6 | 5 | 0 | 0.175 |

### HPRC frequency per category, d9 VCF AF (called haplotypes, matching record) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 21 | 0 | 0 | 0 | 1 | 15 | 5 | 0 | 0.393 |
| normal_present_rare | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_rare | INDEL | 7 | 0 | 0 | 1 | 3 | 3 | 0 | 0 | 0.19 |
| normal_absent_HPRC_other | SNV | 62 | 0 | 3 | 1 | 25 | 25 | 6 | 2 | 0.202 |
| normal_absent_HPRC_other | INDEL | 302 | 0 | 5 | 9 | 138 | 120 | 27 | 3 | 0.198 |
| other_ambiguous | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| other_ambiguous | INDEL | 3 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0.173 |

### HPRC frequency per category, full-graph VCF AF (all called haplotypes, matching record) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 21 | 0 | 0 | 0 | 3 | 14 | 4 | 0 | 0.393 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.123 |
| normal_present_rare | INDEL | 7 | 0 | 2 | 0 | 4 | 1 | 0 | 0 | 0.135 |
| normal_absent_HPRC_other | SNV | 63 | 0 | 2 | 1 | 27 | 25 | 7 | 1 | 0.202 |
| normal_absent_HPRC_other | INDEL | 312 | 0 | 13 | 9 | 156 | 115 | 18 | 1 | 0.18 |
| other_ambiguous | SNV | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0.011 |
| other_ambiguous | INDEL | 5 | 0 | 3 | 0 | 1 | 1 | 0 | 0 | 0.045 |

### Superpopulation of HPRC carriers per category: observed / expected (union, chr1-22)

expected per locus = carriers x the superpopulation's share of the complete haplotypes at that locus (hypergeometric variance for z; ignores the pairing of haplotypes within individuals and linkage between loci, so z overstates significance); after z: 95 % interval of O/E from 2,000 locus bootstrap resamples (seed 20261005), which keeps the loci's carriers together. COLO829 donor: European (white male); HPRC v1.1 has no EUR sample (AFR 23, AMR 16, EAS 4, SAS 1).

| category | level | loci | carrier haps | AFR O/E (z; 95% CI) | AMR O/E (z; 95% CI) | EAS O/E (z; 95% CI) | SAS O/E (z; 95% CI) | CHM13 carries |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | exact allele | 25 | 761 | 1.076 (z +3.0; 0.96-1.19) | 0.902 (z -2.8; 0.77-1.04) | 0.964 (z -0.4; 0.72-1.18) | 0.977 (z -0.1; 0.59-1.41) | 4 |
| normal_present_broad | set, any ALT path | 25 | 770 | 1.070 (z +2.8; 0.95-1.18) | 0.908 (z -2.7; 0.78-1.05) | 0.979 (z -0.2; 0.74-1.19) | 0.959 (z -0.2; 0.58-1.38) | 4 |
| normal_present_rare | exact allele | 9 | 73 | 1.015 (z +0.1; 0.71-1.32) | 0.914 (z -0.6; 0.63-1.28) | 1.240 (z +0.7; 0.16-2.98) | 1.082 (z +0.1; 0.00-3.07) | 3 |
| normal_present_rare | set, any ALT path | 9 | 101 | 1.011 (z +0.1; 0.77-1.25) | 0.928 (z -0.7; 0.70-1.19) | 1.188 (z +0.8; 0.34-2.48) | 1.161 (z +0.3; 0.34-2.52) | 3 |
| normal_absent_HPRC_other | exact allele | 430 | 6656 | 1.175 (z +17.4; 1.14-1.21) | 0.828 (z -13.3; 0.79-0.87) | 0.798 (z -6.5; 0.72-0.88) | 0.861 (z -2.1; 0.73-0.99) | 77 |
| normal_absent_HPRC_other | set, any ALT path | 430 | 7975 | 1.168 (z +18.7; 1.13-1.20) | 0.834 (z -14.4; 0.80-0.87) | 0.807 (z -6.9; 0.74-0.88) | 0.885 (z -2.0; 0.77-1.01) | 97 |
| other_ambiguous | exact allele | 28 | 47 | 1.107 (z +0.8; 0.71-1.37) | 0.923 (z -0.5; 0.73-1.20) | 0.891 (z -0.2; 0.00-1.43) | 0.651 (z -0.5; 0.00-2.00) | 1 |
| other_ambiguous | set, any ALT path | 28 | 373 | 0.915 (z -2.1; 0.75-1.05) | 1.090 (z +1.9; 0.95-1.24) | 1.009 (z +0.1; 0.70-1.36) | 1.118 (z +0.5; 0.57-1.70) | 7 |
| all | exact allele | 492 | 7537 | 1.163 (z +17.5; 1.13-1.20) | 0.837 (z -13.6; 0.79-0.88) | 0.819 (z -6.2; 0.74-0.90) | 0.873 (z -2.1; 0.75-0.99) | 85 |
| all | set, any ALT path | 492 | 9219 | 1.149 (z +17.9; 1.12-1.18) | 0.852 (z -14.0; 0.82-0.89) | 0.833 (z -6.5; 0.77-0.90) | 0.905 (z -1.8; 0.80-1.01) | 111 |

### HPRC support x category (union, chr1-22)

| hprc_support | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| element_set_other_allele | 0 | 0 | 0 | 1 | 5 | 33 | 1 | 15 |
| exact_allele | 0 | 25 | 1 | 7 | 58 | 332 | 0 | 3 |
| none | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 10 |
| recombinant_pieces | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 |

### exact-allele frequency class x category (union, chr1-22)

| exact_freq_class | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 3 |
| 0.2-0.5 | 0 | 19 | 0 | 0 | 24 | 112 | 0 | 2 |
| <0.2 | 0 | 0 | 1 | 8 | 33 | 232 | 1 | 25 |
| >=0.5 | 0 | 6 | 0 | 0 | 6 | 23 | 0 | 0 |

### element-level frequency class x category (union, chr1-22)

| elem_freq_class | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| 0.2-0.5 | 0 | 18 | 0 | 0 | 25 | 153 | 0 | 6 |
| <0.2 | 0 | 1 | 1 | 7 | 30 | 184 | 1 | 16 |
| >=0.5 | 0 | 6 | 0 | 1 | 8 | 30 | 0 | 5 |
| NA | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 3 |

### frequency class agrees across element / set / exact / full-VCF x category (union, chr1-22)

| freq_class_agree | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 3 | 0 | 2 | 4 | 81 | 0 | 10 |
| True | 0 | 22 | 1 | 6 | 59 | 286 | 3 | 20 |

### COLO829BL hap carrying the ALT (c3 array allele, or the with_germline path allele) x category (union, chr1-22)

| normal_pattern | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| hapX_only | 0 | 24 | 1 | 7 | 0 | 0 | 0 | 0 |
| hapY_only | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| neither | 0 | 0 | 0 | 0 | 63 | 367 | 3 | 28 |
| unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

### COLO829BL carries: the ALT (GRCh38 + truth) / the with_germline ALT path allele x category (union, chr1-22)

| normal_carries_via | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 63 | 367 | 3 | 29 |
| ALT | 0 | 24 | 1 | 8 | 0 | 0 | 0 | 1 |
| path_allele | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |

### every exact carrier has a GRCh38 copy too (paralog-like) x category (union, chr1-22)

| multicopy_psv_like | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 24 | 0 | 8 | 63 | 367 | 3 | 30 |
| True | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |

### COLO829BL allele class (locus) x category (union, chr1-22)

| normal_allele_class | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| both_NA | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| both_REF | 0 | 0 | 0 | 0 | 51 | 214 | 1 | 1 |
| carries_ALT | 0 | 24 | 1 | 8 | 0 | 0 | 0 | 1 |
| germline_other_length | 0 | 1 | 0 | 0 | 7 | 139 | 2 | 26 |
| germline_same_length_other_seq | 0 | 0 | 0 | 0 | 5 | 14 | 0 | 1 |

### c3 event hap x category (union, chr1-22)

| event_hap | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| hapX | 0 | 24 | 1 | 7 | 58 | 327 | 3 | 19 |
| hapY | 0 | 1 | 0 | 1 | 5 | 40 | 0 | 10 |

### a COLO829BL hap has the ALT length, other bases (no hap ALT) x category (union, chr1-22)

| germline_like_alt_len | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 25 | 1 | 8 | 63 | 343 | 3 | 20 |
| True | 0 | 0 | 0 | 0 | 0 | 24 | 0 | 10 |

### with_germline: germline alleles on the event hap x category (union, chr1-22)

| germline_phase | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 24 | 1 | 8 | 63 | 356 | 3 | 27 |
| consistent | 0 | 1 | 0 | 0 | 0 | 11 | 0 | 1 |
| off_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |

### path choice x category (union, chr1-22)

| path_choice | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| closest_primary | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 21 |
| no_read_data:alternative:phase | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| no_read_data:primary | 0 | 1 | 0 | 0 | 1 | 30 | 0 | 2 |
| single_path | 0 | 24 | 1 | 8 | 62 | 336 | 0 | 7 |

### element source x category (union, chr1-22)

| element_source | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| graph | 0 | 25 | 1 | 8 | 63 | 367 | 0 | 9 |
| graph_closest_primary | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 21 |

### closest: the graph allele x category (union, chr1-22)

| closest_kind | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 25 | 1 | 8 | 63 | 367 | 0 | 9 |
| graph_allele_nearer_truth | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 5 |
| graph_allele_not_truth | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 16 |

### c2 match x category (union, chr1-22)

| match | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| closest | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 21 |
| exact | 0 | 24 | 1 | 8 | 63 | 356 | 0 | 6 |
| with_germline | 0 | 1 | 0 | 0 | 0 | 11 | 0 | 3 |

### tandem array reaches beyond the graph window x category (union, chr1-22)

| allele_window_clipped | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 3 |
| False | 0 | 25 | 1 | 8 | 61 | 358 | 0 | 26 |
| True | 0 | 0 | 0 | 0 | 2 | 9 | 1 | 1 |

### HPRC carriers: GFA allele level vs VCF carriers_any (on GFA-complete haplotypes; loci with somatic elements)

| VCF | kind | loci | with matching record | identical carrier sets | identical, with record | hap-level disagreements |
|---|---|---|---|---|---|---|
| d9 | SNV | 65 | 62 | 52 | 50 | vcf_only:gfa_set_other_allele 143; gfa_only:no_matching_record 10 |
| d9 | INDEL | 427 | 333 | 308 | 277 | vcf_only:gfa_set_other_allele 809; gfa_only:no_matching_record 776; vcf_only:gfa_other_walk 43; gfa_only:vcf_other_allele 28 |
| full | SNV | 65 | 65 | 52 | 52 | vcf_only:gfa_set_other_allele 156 |
| full | INDEL | 427 | 345 | 312 | 285 | vcf_only:gfa_set_other_allele 802; gfa_only:no_matching_record 716; gfa_only:vcf_other_allele 48; vcf_only:gfa_other_walk 43 |

### d9 frequency-filter floor (somatic elements)

node coverage = HPRC + CHM13 + GRCh38 haplotypes whose c5 rows (complete or partial path) visit the branch node (min over the nodes of the element); a lower bound (c5 partial paths stop 500 nodes from an anchor). full_window = loci where all 88 HPRC haplotypes have a complete traversal, where the counts are exact.

| measure | value |
|---|---|
| full_windows | 46 |
| full_window_elements | 46 |
| full_window_branch_node_coverage_min | 2 |
| full_window_node_coverage_hist | {2: 1, 4: 1, 5: 1, 9: 2, 10: 1, 12: 2, 15: 23} |
| full_window_hprc_carriers_min | 0 |
| full_window_hprc_carriers_hist | {0: 1, 1: 1, 3: 1, 4: 1, 6: 1, 9: 2, 10: 1, 11: 2, 12: 3, 13: 1, 14: 4, 15: 28} |
| full_window_not_on_CHM13_hprc_carriers_min | 9 |
| full_window_on_CHM13_below9 | 5 |
| all_not_on_CHM13_node_coverage_lt9 | 2 |
| n_elements | 623 |
| branch_elements | 413 |
| branch_elements_node_coverage_min | 2 |
| node_coverage_lt7 | 9 |
| node_coverage_lt9 | 12 |
| node_coverage_hist | {2: 1, 3: 2, 4: 2, 5: 1, 6: 3, 7: 1, 8: 2, 9: 11, 10: 20, 11: 14, 12: 10, 13: 15, 14: 10, 15: 321} |
| complete_carriers_min | 0 |
| complete_carriers_lt9 | 76 |
| any_row_carriers_min | 0 |
| any_row_carriers_lt9 | 32 |
