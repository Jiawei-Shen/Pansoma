# c7: HPRC membership and cross-evaluation of the graph-absorbed COLO829T somatic truth alleles (d9 HPRC v1.1)
RARE_AF = 0.2 on the exact-allele HPRC frequency (complete haplotypes spelling the ALT allele over the tandem array). COLO829BL calls: c3 whole-array alleles. No read data.

### Categories per perfect set, chr1-22

| category | fiberseq SNV chr1-22 | fiberseq INDEL chr1-22 | ONT SNV chr1-22 | ONT INDEL chr1-22 | Illumina SNV chr1-22 | Illumina INDEL chr1-22 | union SNV chr1-22 | union INDEL chr1-22 |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | 0 | 18 | 0 | 17 | 0 | 20 | 0 | 24 |
| normal_present_rare | 1 | 8 | 1 | 4 | 1 | 6 | 1 | 8 |
| normal_absent_HPRC_other | 59 | 175 | 54 | 129 | 50 | 345 | 62 | 365 |
| other_ambiguous | 4 | 17 | 3 | 20 | 2 | 25 | 4 | 33 |
| total | 64 | 218 | 58 | 170 | 53 | 396 | 67 | 430 |

### Categories per perfect set, chr1

| category | fiberseq SNV chr1 | fiberseq INDEL chr1 | ONT SNV chr1 | ONT INDEL chr1 | Illumina SNV chr1 | Illumina INDEL chr1 | union SNV chr1 | union INDEL chr1 |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | 0 | 7 | 0 | 7 | 0 | 9 | 0 | 10 |
| normal_present_rare | 0 | 3 | 0 | 1 | 0 | 2 | 0 | 3 |
| normal_absent_HPRC_other | 3 | 12 | 3 | 9 | 3 | 22 | 3 | 22 |
| other_ambiguous | 2 | 1 | 1 | 2 | 1 | 1 | 2 | 3 |
| total | 5 | 23 | 4 | 19 | 4 | 34 | 5 | 38 |

### Categories per perfect set, chr2-22

| category | fiberseq SNV chr2-22 | fiberseq INDEL chr2-22 | ONT SNV chr2-22 | ONT INDEL chr2-22 | Illumina SNV chr2-22 | Illumina INDEL chr2-22 | union SNV chr2-22 | union INDEL chr2-22 |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | 0 | 11 | 0 | 10 | 0 | 11 | 0 | 14 |
| normal_present_rare | 1 | 5 | 1 | 3 | 1 | 4 | 1 | 5 |
| normal_absent_HPRC_other | 56 | 163 | 51 | 120 | 47 | 323 | 59 | 343 |
| other_ambiguous | 2 | 16 | 2 | 18 | 1 | 24 | 2 | 30 |
| total | 59 | 195 | 54 | 151 | 49 | 362 | 62 | 392 |

### Categories with sub-class / sub-reason, chr1-22

| category | sub_class / sub_reason | fiberseq SNV chr1-22 | fiberseq INDEL chr1-22 | ONT SNV chr1-22 | ONT INDEL chr1-22 | Illumina SNV chr1-22 | Illumina INDEL chr1-22 | union SNV chr1-22 | union INDEL chr1-22 |
|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | hapX_only | 0 | 18 | 0 | 17 | 0 | 19 | 0 | 23 |
| normal_present_broad | hapY_only | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| normal_present_rare | hapX_only | 1 | 7 | 1 | 4 | 1 | 5 | 1 | 7 |
| normal_present_rare | hapY_only | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| normal_absent_HPRC_other | element_set_other_allele | 4 | 24 | 2 | 21 | 0 | 26 | 4 | 33 |
| normal_absent_HPRC_other | exact_allele | 55 | 149 | 52 | 106 | 50 | 317 | 58 | 330 |
| normal_absent_HPRC_other | recombinant_pieces | 0 | 2 | 0 | 2 | 0 | 2 | 0 | 2 |
| other_ambiguous | closest:graph_allele_nearer_truth | 1 | 2 | 1 | 3 | 0 | 5 | 1 | 5 |
| other_ambiguous | closest:graph_allele_not_truth | 2 | 11 | 1 | 13 | 2 | 11 | 2 | 16 |
| other_ambiguous | no_HPRC_carrier | 1 | 4 | 1 | 3 | 0 | 5 | 1 | 7 |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 2 |
| other_ambiguous | normal_hap_unresolved:both_NA | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 |

### Categories with sub-class / sub-reason, chr1

| category | sub_class / sub_reason | fiberseq SNV chr1 | fiberseq INDEL chr1 | ONT SNV chr1 | ONT INDEL chr1 | Illumina SNV chr1 | Illumina INDEL chr1 | union SNV chr1 | union INDEL chr1 |
|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | hapX_only | 0 | 7 | 0 | 7 | 0 | 9 | 0 | 10 |
| normal_present_broad | hapY_only | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| normal_present_rare | hapX_only | 0 | 3 | 0 | 1 | 0 | 2 | 0 | 3 |
| normal_present_rare | hapY_only | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| normal_absent_HPRC_other | element_set_other_allele | 0 | 3 | 0 | 2 | 0 | 3 | 0 | 3 |
| normal_absent_HPRC_other | exact_allele | 3 | 9 | 3 | 7 | 3 | 19 | 3 | 19 |
| normal_absent_HPRC_other | recombinant_pieces | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | closest:graph_allele_nearer_truth | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | closest:graph_allele_not_truth | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 |
| other_ambiguous | no_HPRC_carrier | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 2 |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| other_ambiguous | normal_hap_unresolved:both_NA | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Interpretation per category and sub-class (union, chr1-22)

| category | sub_class / sub_reason | SNV | INDEL | interpretation |
|---|---|---|---|---|
| normal_present_broad | hapX_only | 0 | 23 | COLO829BL germline allele already a graph path (common exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapX); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_present_broad | hapY_only | 0 | 1 | COLO829BL germline allele already a graph path (common exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapY); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_present_rare | hapX_only | 1 | 7 | COLO829BL germline allele already a graph path (rare exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapX); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_present_rare | hapY_only | 0 | 1 | COLO829BL germline allele already a graph path (rare exact allele in HPRC): the ALT allele over the whole tandem array is a COLO829BL germline allele (hapY); a tumor-only caller sees a germline allele: germline filtering by the graph |
| normal_absent_HPRC_other | element_set_other_allele | 4 | 33 | the patient does not carry it; HPRC haplotypes walk all somatic elements but spell another allele over the array (e.g. another repeat length): the graph spells the somatic allele from population nodes, not a population allele: pangenome-induced false negative |
| normal_absent_HPRC_other | exact_allele | 58 | 330 | the patient does not carry it; >= 1 complete HPRC haplotype spells the exact ALT allele over the array (somatic slippage / substitution recreating a population allele kept by the d9 graph): pangenome-induced false negative |
| normal_absent_HPRC_other | recombinant_pieces | 0 | 2 | the patient does not carry it; every somatic element is walked by some HPRC haplotype but none walks a whole ALT path: the ALT path is a recombination of population nodes: pangenome-induced false negative |
| other_ambiguous | closest:graph_allele_nearer_truth | 1 | 5 | unresolved: no d9 path spells the truth; the closest graph allele is nearer the truth, but not the truth (no read data) |
| other_ambiguous | closest:graph_allele_not_truth | 2 | 16 | unresolved: no d9 path spells the truth; the closest graph allele is nearer GRCh38 / the patient's germline than the truth (no read data to say what the reads walk) |
| other_ambiguous | no_HPRC_carrier | 1 | 7 | unresolved: no complete HPRC haplotype walks the somatic elements |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 2 | unresolved: no complete HPRC haplotype walks the somatic elements; CHM13 does |
| other_ambiguous | normal_hap_unresolved:both_NA | 0 | 1 | unresolved: no COLO829BL hap has an array sequence |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 2 | unresolved: the ALT path needs a germline allele that dipcall phases to the other hap than the c3 event hap |

### HPRC frequency per category, exact allele (spells the ALT allele over the array; RARE_AF level) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 24 | 0 | 0 | 0 | 0 | 18 | 6 | 0 | 0.399 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.125 |
| normal_present_rare | INDEL | 8 | 1 | 1 | 1 | 5 | 0 | 0 | 0 | 0.152 |
| normal_absent_HPRC_other | SNV | 62 | 4 | 6 | 1 | 21 | 24 | 5 | 1 | 0.189 |
| normal_absent_HPRC_other | INDEL | 365 | 35 | 32 | 28 | 136 | 110 | 24 | 0 | 0.159 |
| other_ambiguous | SNV | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 |
| other_ambiguous | INDEL | 29 | 26 | 0 | 0 | 1 | 2 | 0 | 0 | 0.0 |

### HPRC frequency per category, exact allele over all 88 HPRC haplotypes (lower bound) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 24 | 0 | 0 | 0 | 5 | 15 | 4 | 0 | 0.335 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.114 |
| normal_present_rare | INDEL | 8 | 1 | 2 | 0 | 5 | 0 | 0 | 0 | 0.108 |
| normal_absent_HPRC_other | SNV | 62 | 4 | 6 | 1 | 23 | 22 | 5 | 1 | 0.182 |
| normal_absent_HPRC_other | INDEL | 365 | 35 | 42 | 30 | 150 | 93 | 15 | 0 | 0.136 |
| other_ambiguous | SNV | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 |
| other_ambiguous | INDEL | 30 | 27 | 0 | 0 | 2 | 1 | 0 | 0 | 0.0 |

### HPRC frequency per category, set level, any ALT path (walks all somatic elements) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 24 | 0 | 0 | 0 | 0 | 18 | 6 | 0 | 0.399 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.125 |
| normal_present_rare | INDEL | 8 | 0 | 1 | 1 | 5 | 0 | 1 | 0 | 0.172 |
| normal_absent_HPRC_other | SNV | 62 | 0 | 4 | 1 | 24 | 25 | 6 | 2 | 0.205 |
| normal_absent_HPRC_other | INDEL | 365 | 2 | 18 | 26 | 142 | 144 | 30 | 3 | 0.191 |
| other_ambiguous | SNV | 2 | 1 | 0 | 0 | 1 | 0 | 0 | 0 | 0.079 |
| other_ambiguous | INDEL | 29 | 11 | 3 | 0 | 5 | 6 | 4 | 0 | 0.1 |

### HPRC frequency per category, element level (min over elements) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 24 | 0 | 0 | 0 | 1 | 17 | 6 | 0 | 0.399 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.125 |
| normal_present_rare | INDEL | 8 | 0 | 0 | 1 | 6 | 0 | 1 | 0 | 0.172 |
| normal_absent_HPRC_other | SNV | 62 | 0 | 4 | 1 | 24 | 25 | 6 | 2 | 0.205 |
| normal_absent_HPRC_other | INDEL | 365 | 0 | 17 | 19 | 147 | 150 | 30 | 2 | 0.198 |
| other_ambiguous | SNV | 2 | 1 | 0 | 0 | 1 | 0 | 0 | 0 | 0.079 |
| other_ambiguous | INDEL | 29 | 10 | 3 | 0 | 6 | 5 | 5 | 0 | 0.113 |

### HPRC frequency per category, d9 VCF AF (called haplotypes, matching record) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 21 | 0 | 0 | 0 | 1 | 15 | 5 | 0 | 0.393 |
| normal_present_rare | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_rare | INDEL | 7 | 0 | 0 | 1 | 3 | 3 | 0 | 0 | 0.19 |
| normal_absent_HPRC_other | SNV | 61 | 0 | 3 | 1 | 24 | 25 | 6 | 2 | 0.202 |
| normal_absent_HPRC_other | INDEL | 301 | 0 | 5 | 9 | 138 | 119 | 27 | 3 | 0.197 |
| other_ambiguous | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.112 |
| other_ambiguous | INDEL | 4 | 0 | 1 | 0 | 1 | 2 | 0 | 0 | 0.216 |

### HPRC frequency per category, full-graph VCF AF (all called haplotypes, matching record) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| normal_present_broad | SNV | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |
| normal_present_broad | INDEL | 21 | 0 | 0 | 0 | 3 | 14 | 4 | 0 | 0.393 |
| normal_present_rare | SNV | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.123 |
| normal_present_rare | INDEL | 7 | 0 | 2 | 0 | 4 | 1 | 0 | 0 | 0.135 |
| normal_absent_HPRC_other | SNV | 62 | 0 | 2 | 1 | 26 | 25 | 7 | 1 | 0.202 |
| normal_absent_HPRC_other | INDEL | 311 | 0 | 13 | 9 | 156 | 114 | 18 | 1 | 0.18 |
| other_ambiguous | SNV | 2 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0.062 |
| other_ambiguous | INDEL | 6 | 0 | 3 | 0 | 1 | 2 | 0 | 0 | 0.102 |

### Superpopulation of HPRC carriers per category: observed / expected (union, chr1-22)

expected per locus = carriers x the superpopulation's share of the complete haplotypes at that locus (hypergeometric variance for z; ignores the pairing of haplotypes within individuals and linkage between loci, so z overstates significance). COLO829 donor: European (white male); HPRC v1.1 has no EUR sample (AFR 23, AMR 16, EAS 4, SAS 1).

| category | level | loci | carrier haps | AFR O/E (z) | AMR O/E (z) | EAS O/E (z) | SAS O/E (z) | CHM13 carries |
|---|---|---|---|---|---|---|---|---|
| normal_present_broad | exact allele | 24 | 747 | 1.067 (z +2.7) | 0.918 (z -2.3) | 0.957 (z -0.5) | 0.941 (z -0.3) | 3 |
| normal_present_broad | set, any ALT path | 24 | 756 | 1.062 (z +2.4) | 0.924 (z -2.2) | 0.973 (z -0.3) | 0.923 (z -0.4) | 3 |
| normal_present_rare | exact allele | 9 | 73 | 1.015 (z +0.1) | 0.914 (z -0.6) | 1.240 (z +0.7) | 1.082 (z +0.1) | 3 |
| normal_present_rare | set, any ALT path | 9 | 101 | 1.011 (z +0.1) | 0.928 (z -0.7) | 1.188 (z +0.8) | 1.161 (z +0.3) | 3 |
| normal_absent_HPRC_other | exact allele | 427 | 6594 | 1.178 (z +17.7) | 0.825 (z -13.5) | 0.795 (z -6.6) | 0.855 (z -2.2) | 75 |
| normal_absent_HPRC_other | set, any ALT path | 427 | 7876 | 1.172 (z +19.0) | 0.831 (z -14.6) | 0.804 (z -7.0) | 0.877 (z -2.1) | 96 |
| other_ambiguous | exact allele | 31 | 44 | 1.065 (z +0.5) | 0.937 (z -0.4) | 1.144 (z +0.3) | 0.660 (z -0.5) | 1 |
| other_ambiguous | set, any ALT path | 31 | 358 | 0.906 (z -2.3) | 1.106 (z +2.2) | 0.969 (z -0.3) | 1.192 (z +0.8) | 7 |
| all | exact allele | 491 | 7458 | 1.164 (z +17.6) | 0.835 (z -13.6) | 0.816 (z -6.3) | 0.864 (z -2.2) | 82 |
| all | set, any ALT path | 491 | 9091 | 1.150 (z +18.1) | 0.851 (z -14.0) | 0.828 (z -6.7) | 0.898 (z -1.9) | 109 |

### HPRC support x category (union, chr1-22)

| hprc_support | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| element_set_other_allele | 0 | 0 | 0 | 1 | 4 | 33 | 1 | 14 |
| exact_allele | 0 | 24 | 1 | 7 | 58 | 330 | 0 | 3 |
| none | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 14 |
| recombinant_pieces | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 |

### exact-allele frequency class x category (union, chr1-22)

| exact_freq_class | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 3 |
| 0.2-0.5 | 0 | 18 | 0 | 0 | 24 | 110 | 0 | 2 |
| <0.2 | 0 | 0 | 1 | 8 | 32 | 231 | 2 | 27 |
| >=0.5 | 0 | 6 | 0 | 0 | 6 | 24 | 0 | 0 |
| NA | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

### element-level frequency class x category (union, chr1-22)

| elem_freq_class | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| 0.2-0.5 | 0 | 17 | 0 | 0 | 25 | 150 | 0 | 5 |
| <0.2 | 0 | 1 | 1 | 7 | 29 | 183 | 2 | 19 |
| >=0.5 | 0 | 6 | 0 | 1 | 8 | 32 | 0 | 5 |
| NA | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 4 |

### frequency class agrees across element / set / exact / full-VCF x category (union, chr1-22)

| freq_class_agree | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 3 | 0 | 2 | 4 | 81 | 0 | 11 |
| True | 0 | 21 | 1 | 6 | 58 | 284 | 4 | 22 |

### COLO829BL hap carrying the ALT (c3 array allele) x category (union, chr1-22)

| normal_pattern | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| hapX_only | 0 | 23 | 1 | 7 | 0 | 0 | 0 | 0 |
| hapY_only | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| neither | 0 | 0 | 0 | 0 | 62 | 365 | 4 | 31 |
| unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

### COLO829BL allele class (locus) x category (union, chr1-22)

| normal_allele_class | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| both_NA | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| both_REF | 0 | 0 | 0 | 0 | 51 | 214 | 1 | 1 |
| carries_ALT | 0 | 24 | 1 | 8 | 0 | 0 | 0 | 1 |
| germline_other_length | 0 | 0 | 0 | 0 | 6 | 137 | 3 | 29 |
| germline_same_length_other_seq | 0 | 0 | 0 | 0 | 5 | 14 | 0 | 1 |

### c3 event hap x category (union, chr1-22)

| event_hap | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| hapX | 0 | 23 | 1 | 7 | 58 | 325 | 3 | 22 |
| hapY | 0 | 1 | 0 | 1 | 4 | 40 | 1 | 10 |

### a COLO829BL hap has the ALT length, other bases (no hap ALT) x category (union, chr1-22)

| germline_like_alt_len | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 24 | 1 | 8 | 62 | 343 | 4 | 20 |
| True | 0 | 0 | 0 | 0 | 0 | 22 | 0 | 13 |

### with_germline: germline alleles on the event hap x category (union, chr1-22)

| germline_phase | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 24 | 1 | 8 | 62 | 354 | 4 | 29 |
| consistent | 0 | 0 | 0 | 0 | 0 | 11 | 0 | 2 |
| off_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |

### path choice x category (union, chr1-22)

| path_choice | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| closest_primary | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 21 |
| no_read_data:alternative:phase | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| no_read_data:primary | 0 | 1 | 0 | 0 | 1 | 29 | 0 | 3 |
| single_path | 0 | 23 | 1 | 8 | 61 | 335 | 1 | 9 |

### element source x category (union, chr1-22)

| element_source | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| graph | 0 | 24 | 1 | 8 | 62 | 365 | 1 | 12 |
| graph_closest_primary | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 21 |

### closest: the graph allele x category (union, chr1-22)

| closest_kind | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 24 | 1 | 8 | 62 | 365 | 1 | 12 |
| graph_allele_nearer_truth | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 5 |
| graph_allele_not_truth | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 16 |

### c2 match x category (union, chr1-22)

| match | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| closest | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 21 |
| exact | 0 | 24 | 1 | 8 | 62 | 354 | 1 | 8 |
| with_germline | 0 | 0 | 0 | 0 | 0 | 11 | 0 | 4 |

### tandem array reaches beyond the graph window x category (union, chr1-22)

| allele_window_clipped | normal_present_broad SNV | normal_present_broad INDEL | normal_present_rare SNV | normal_present_rare INDEL | normal_absent_HPRC_other SNV | normal_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 3 |
| False | 0 | 24 | 1 | 8 | 61 | 356 | 0 | 29 |
| True | 0 | 0 | 0 | 0 | 1 | 9 | 2 | 1 |

### HPRC carriers: GFA allele level vs VCF carriers_any (on GFA-complete haplotypes; loci with somatic elements)

| VCF | kind | loci | with matching record | identical carrier sets | identical, with record | hap-level disagreements |
|---|---|---|---|---|---|---|
| d9 | SNV | 65 | 62 | 53 | 51 | vcf_only:gfa_set_other_allele 133; gfa_only:no_matching_record 10 |
| d9 | INDEL | 427 | 333 | 308 | 275 | vcf_only:gfa_set_other_allele 791; gfa_only:no_matching_record 737; vcf_only:gfa_other_walk 54; gfa_only:vcf_other_allele 28 |
| full | SNV | 65 | 65 | 53 | 53 | vcf_only:gfa_set_other_allele 146 |
| full | INDEL | 427 | 345 | 312 | 283 | vcf_only:gfa_set_other_allele 784; gfa_only:no_matching_record 677; vcf_only:gfa_other_walk 54; gfa_only:vcf_other_allele 48 |

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
| complete_carriers_lt9 | 83 |
| any_row_carriers_min | 0 |
| any_row_carriers_lt9 | 32 |
