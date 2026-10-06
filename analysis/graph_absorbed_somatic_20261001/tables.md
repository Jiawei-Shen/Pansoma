# s9 integration: graph-absorbed HG008T somatic truth alleles (d9 HPRC v1.1)
RARE_AF = 0.2 on the exact-allele HPRC frequency (complete haplotypes spelling the ALT allele over the tandem array). HG008-N calls: whole-array alleles (s8b). Read data: PacBio 2478 loci, ONT 1365 loci, Illumina 2478 loci. Fixes: audit/decisions.md.

### Categories per perfect set, chr1-22

| category | PacBio SNV chr1-22 | PacBio INDEL chr1-22 | ONT SNV chr1-22 | ONT INDEL chr1-22 | Illumina SNV chr1-22 | Illumina INDEL chr1-22 | union SNV chr1-22 | union INDEL chr1-22 |
|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | 2 | 41 | 2 | 26 | 2 | 68 | 2 | 75 |
| HG008N_present_rare | 3 | 17 | 3 | 10 | 2 | 28 | 3 | 30 |
| HG008N_absent_HPRC_other | 129 | 1157 | 117 | 806 | 97 | 1682 | 150 | 1935 |
| other_ambiguous | 36 | 148 | 26 | 147 | 21 | 138 | 46 | 237 |
| total | 170 | 1363 | 148 | 989 | 122 | 1916 | 201 | 2277 |

### Categories per perfect set, chr2-22

| category | PacBio SNV chr2-22 | PacBio INDEL chr2-22 | ONT SNV chr2-22 | ONT INDEL chr2-22 | Illumina SNV chr2-22 | Illumina INDEL chr2-22 | union SNV chr2-22 | union INDEL chr2-22 |
|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | 2 | 36 | 2 | 24 | 2 | 59 | 2 | 66 |
| HG008N_present_rare | 3 | 16 | 3 | 9 | 2 | 26 | 3 | 28 |
| HG008N_absent_HPRC_other | 118 | 1064 | 105 | 740 | 87 | 1555 | 138 | 1781 |
| other_ambiguous | 33 | 138 | 24 | 136 | 19 | 128 | 41 | 217 |
| total | 156 | 1254 | 134 | 909 | 110 | 1768 | 184 | 2092 |

### Categories per perfect set, chr1 BED

| category | PacBio SNV chr1 BED | PacBio INDEL chr1 BED | ONT SNV chr1 BED | ONT INDEL chr1 BED | Illumina SNV chr1 BED | Illumina INDEL chr1 BED | union SNV chr1 BED | union INDEL chr1 BED |
|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | 0 | 5 | 0 | 2 | 0 | 8 | 0 | 8 |
| HG008N_present_rare | 0 | 1 | 0 | 1 | 0 | 2 | 0 | 2 |
| HG008N_absent_HPRC_other | 11 | 86 | 12 | 60 | 10 | 117 | 12 | 140 |
| other_ambiguous | 3 | 10 | 2 | 11 | 2 | 9 | 5 | 19 |
| total | 14 | 102 | 14 | 74 | 12 | 136 | 17 | 169 |

### Categories per perfect set, chr1-22 nogermline BED

| category | PacBio SNV chr1-22 nogermline BED | PacBio INDEL chr1-22 nogermline BED | ONT SNV chr1-22 nogermline BED | ONT INDEL chr1-22 nogermline BED | Illumina SNV chr1-22 nogermline BED | Illumina INDEL chr1-22 nogermline BED | union SNV chr1-22 nogermline BED | union INDEL chr1-22 nogermline BED |
|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | 2 | 16 | 2 | 9 | 2 | 42 | 2 | 43 |
| HG008N_present_rare | 2 | 4 | 2 | 1 | 1 | 13 | 2 | 14 |
| HG008N_absent_HPRC_other | 55 | 294 | 49 | 198 | 49 | 803 | 64 | 840 |
| other_ambiguous | 2 | 3 | 3 | 3 | 1 | 4 | 3 | 6 |
| total | 61 | 317 | 56 | 211 | 53 | 862 | 71 | 903 |

### Categories with sub-class / sub-reason, chr1-22

| category | sub_class / sub_reason | PacBio SNV chr1-22 | PacBio INDEL chr1-22 | ONT SNV chr1-22 | ONT INDEL chr1-22 | Illumina SNV chr1-22 | Illumina INDEL chr1-22 | union SNV chr1-22 | union INDEL chr1-22 |
|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | other_hap_only | 2 | 41 | 2 | 26 | 2 | 68 | 2 | 75 |
| HG008N_present_rare | no_event_hap | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 |
| HG008N_present_rare | other_hap_only | 1 | 10 | 1 | 4 | 0 | 21 | 1 | 22 |
| HG008N_present_rare | other_hap_only:patient_frame | 2 | 6 | 2 | 5 | 2 | 6 | 2 | 7 |
| HG008N_absent_HPRC_other | element_set_other_allele | 66 | 74 | 54 | 70 | 36 | 67 | 76 | 101 |
| HG008N_absent_HPRC_other | exact_allele | 55 | 1045 | 57 | 706 | 58 | 1601 | 65 | 1789 |
| HG008N_absent_HPRC_other | recombinant_pieces | 8 | 38 | 6 | 30 | 3 | 14 | 9 | 45 |
| other_ambiguous | HG008N_hap_unresolved:array_NA | 2 | 5 | 2 | 5 | 1 | 3 | 2 | 6 |
| other_ambiguous | HG008N_hap_unresolved:conflict | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 |
| other_ambiguous | closest:read_allele_nearer_truth | 2 | 12 | 2 | 13 | 0 | 12 | 3 | 21 |
| other_ambiguous | closest:read_allele_not_truth | 24 | 99 | 18 | 100 | 17 | 97 | 32 | 166 |
| other_ambiguous | closest_no_read_element | 4 | 0 | 0 | 0 | 0 | 0 | 4 | 0 |
| other_ambiguous | evidence_conflict:ALT_on_both_normal_haps | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| other_ambiguous | evidence_conflict:ALT_on_normal_event_hap | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| other_ambiguous | evidence_conflict:germline_ALT_vs_truth_INFO_event_novel | 0 | 3 | 0 | 2 | 0 | 6 | 0 | 6 |
| other_ambiguous | no_HPRC_carrier | 1 | 13 | 1 | 10 | 0 | 6 | 1 | 13 |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 3 | 0 | 3 | 0 | 5 | 0 | 6 |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 3 | 12 | 2 | 12 | 2 | 9 | 3 | 17 |

### Categories with sub-class / sub-reason, chr1 BED

| category | sub_class / sub_reason | PacBio SNV chr1 BED | PacBio INDEL chr1 BED | ONT SNV chr1 BED | ONT INDEL chr1 BED | Illumina SNV chr1 BED | Illumina INDEL chr1 BED | union SNV chr1 BED | union INDEL chr1 BED |
|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | other_hap_only | 0 | 5 | 0 | 2 | 0 | 8 | 0 | 8 |
| HG008N_present_rare | no_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HG008N_present_rare | other_hap_only | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 |
| HG008N_present_rare | other_hap_only:patient_frame | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| HG008N_absent_HPRC_other | element_set_other_allele | 6 | 6 | 6 | 5 | 4 | 5 | 6 | 7 |
| HG008N_absent_HPRC_other | exact_allele | 5 | 76 | 6 | 52 | 6 | 111 | 6 | 128 |
| HG008N_absent_HPRC_other | recombinant_pieces | 0 | 4 | 0 | 3 | 0 | 1 | 0 | 5 |
| other_ambiguous | HG008N_hap_unresolved:array_NA | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | HG008N_hap_unresolved:conflict | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | closest:read_allele_nearer_truth | 0 | 1 | 0 | 2 | 0 | 3 | 0 | 4 |
| other_ambiguous | closest:read_allele_not_truth | 3 | 7 | 2 | 7 | 2 | 6 | 5 | 13 |
| other_ambiguous | closest_no_read_element | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | evidence_conflict:ALT_on_both_normal_haps | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | evidence_conflict:ALT_on_normal_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | evidence_conflict:germline_ALT_vs_truth_INFO_event_novel | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | no_HPRC_carrier | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |

### Categories with sub-class / sub-reason, chr1-22 nogermline BED

| category | sub_class / sub_reason | PacBio SNV chr1-22 nogermline BED | PacBio INDEL chr1-22 nogermline BED | ONT SNV chr1-22 nogermline BED | ONT INDEL chr1-22 nogermline BED | Illumina SNV chr1-22 nogermline BED | Illumina INDEL chr1-22 nogermline BED | union SNV chr1-22 nogermline BED | union INDEL chr1-22 nogermline BED |
|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | other_hap_only | 2 | 16 | 2 | 9 | 2 | 42 | 2 | 43 |
| HG008N_present_rare | no_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HG008N_present_rare | other_hap_only | 1 | 4 | 1 | 1 | 0 | 12 | 1 | 13 |
| HG008N_present_rare | other_hap_only:patient_frame | 1 | 0 | 1 | 0 | 1 | 1 | 1 | 1 |
| HG008N_absent_HPRC_other | element_set_other_allele | 24 | 6 | 17 | 5 | 16 | 8 | 26 | 10 |
| HG008N_absent_HPRC_other | exact_allele | 31 | 287 | 32 | 192 | 33 | 795 | 38 | 829 |
| HG008N_absent_HPRC_other | recombinant_pieces | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| other_ambiguous | HG008N_hap_unresolved:array_NA | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 |
| other_ambiguous | HG008N_hap_unresolved:conflict | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 |
| other_ambiguous | closest:read_allele_nearer_truth | 0 | 2 | 0 | 2 | 0 | 1 | 0 | 3 |
| other_ambiguous | closest:read_allele_not_truth | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | closest_no_read_element | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | evidence_conflict:ALT_on_both_normal_haps | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | evidence_conflict:ALT_on_normal_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| other_ambiguous | evidence_conflict:germline_ALT_vs_truth_INFO_event_novel | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| other_ambiguous | no_HPRC_carrier | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 1 | 0 | 1 | 0 | 2 | 0 | 2 |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Interpretation per category and sub-class (union, chr1-22)

| category | sub_class / sub_reason | SNV | INDEL | interpretation |
|---|---|---|---|---|
| HG008N_present_broad | other_hap_only | 2 | 75 | HG008-N germline allele already a graph path (common exact allele in HPRC): the ALT allele over the whole tandem array is the patient's germline allele on the non-event hap; the tumor event hap now carries it (local conversion or recurrent slippage to the other hap's allele; a tumor-only caller sees a germline allele): germline filtering by the graph |
| HG008N_present_rare | no_event_hap | 0 | 1 | HG008-N germline allele already a graph path (rare exact allele in HPRC): the ALT allele over the array is a germline allele of one normal hap; no event hap in the truth INFO |
| HG008N_present_rare | other_hap_only | 1 | 22 | HG008-N germline allele already a graph path (rare exact allele in HPRC): the ALT allele over the whole tandem array is the patient's germline allele on the non-event hap; the tumor event hap now carries it (local conversion or recurrent slippage to the other hap's allele; a tumor-only caller sees a germline allele): germline filtering by the graph |
| HG008N_present_rare | other_hap_only:patient_frame | 2 | 7 | HG008-N germline allele already a graph path (rare exact allele in HPRC): the tumor event-hap allele (truth INFO event applied to the normal event hap) equals the patient's other-hap germline allele over the array; the GRCh38 ALT differs from it only by germline differences shared by both haps: germline filtering by the graph |
| HG008N_absent_HPRC_other | element_set_other_allele | 76 | 101 | the patient does not carry it; HPRC haplotypes walk all somatic elements but spell another allele over the array (e.g. another repeat length): the graph spells the somatic allele from population nodes, not a population allele: pangenome-induced false negative |
| HG008N_absent_HPRC_other | exact_allele | 65 | 1789 | the patient does not carry it; >= 1 complete HPRC haplotype spells the exact ALT allele over the array (somatic slippage / substitution recreating a population allele kept by the d9 graph): pangenome-induced false negative |
| HG008N_absent_HPRC_other | recombinant_pieces | 9 | 45 | the patient does not carry it; every somatic element is walked by some HPRC haplotype but none walks a whole ALT path: the ALT path is a recombination of population nodes: pangenome-induced false negative |
| other_ambiguous | HG008N_hap_unresolved:array_NA | 2 | 6 | unresolved: a normal hap has no array sequence and the narrow call is not a clear no |
| other_ambiguous | HG008N_hap_unresolved:conflict | 1 | 0 | unresolved: a normal hap: assembly and dipcall disagree over the array |
| other_ambiguous | closest:read_allele_nearer_truth | 3 | 21 | unresolved: no graph path spells the truth; the reads walk a graph allele nearer the truth, but not the truth |
| other_ambiguous | closest:read_allele_not_truth | 32 | 166 | unresolved: no graph path spells the truth; the reads walk a graph allele nearer GRCh38 / the patient's germline than the truth |
| other_ambiguous | closest_no_read_element | 4 | 0 | unresolved: no graph path spells the truth and no read walk to take the elements from |
| other_ambiguous | evidence_conflict:ALT_on_both_normal_haps | 0 | 1 | unresolved: both normal haps carry the GRCh38 ALT over the array: homozygous germline, truth conflict |
| other_ambiguous | evidence_conflict:ALT_on_normal_event_hap | 0 | 1 | unresolved: the GRCh38 ALT is the normal event hap's own allele; the somatic change relative to the normal is the INFO event: truth representation conflict |
| other_ambiguous | evidence_conflict:germline_ALT_vs_truth_INFO_event_novel | 0 | 6 | unresolved: the GRCh38 ALT equals a normal hap's allele, but the truth INFO event gives a tumor allele that neither normal hap has |
| other_ambiguous | no_HPRC_carrier | 1 | 13 | unresolved: no complete HPRC haplotype walks the somatic elements |
| other_ambiguous | no_HPRC_carrier:CHM13_carries | 0 | 6 | unresolved: no complete HPRC haplotype walks the somatic elements; CHM13 does |
| other_ambiguous | with_germline:germline_allele_off_event_hap | 3 | 17 | unresolved: the ALT path needs a germline allele that dipcall phases to the non-event hap |

### HPRC frequency per category, exact allele (spells the ALT allele over the array; RARE_AF level) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | SNV | 2 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0.374 |
| HG008N_present_broad | INDEL | 75 | 0 | 0 | 0 | 0 | 54 | 20 | 1 | 0.392 |
| HG008N_present_rare | SNV | 3 | 0 | 0 | 1 | 2 | 0 | 0 | 0 | 0.128 |
| HG008N_present_rare | INDEL | 30 | 4 | 2 | 5 | 19 | 0 | 0 | 0 | 0.126 |
| HG008N_absent_HPRC_other | SNV | 150 | 85 | 13 | 16 | 25 | 9 | 2 | 0 | 0.0 |
| HG008N_absent_HPRC_other | INDEL | 1935 | 146 | 141 | 140 | 739 | 714 | 55 | 0 | 0.167 |
| other_ambiguous | SNV | 35 | 33 | 0 | 1 | 0 | 1 | 0 | 0 | 0.0 |
| other_ambiguous | INDEL | 194 | 177 | 3 | 5 | 2 | 5 | 2 | 0 | 0.0 |

### HPRC frequency per category, set level, any ALT path (walks all somatic elements) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | SNV | 2 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0.374 |
| HG008N_present_broad | INDEL | 75 | 0 | 0 | 0 | 0 | 52 | 21 | 2 | 0.392 |
| HG008N_present_rare | SNV | 3 | 0 | 0 | 0 | 0 | 2 | 1 | 0 | 0.282 |
| HG008N_present_rare | INDEL | 30 | 3 | 0 | 3 | 18 | 4 | 0 | 2 | 0.171 |
| HG008N_absent_HPRC_other | SNV | 150 | 9 | 2 | 13 | 29 | 42 | 30 | 25 | 0.322 |
| HG008N_absent_HPRC_other | INDEL | 1935 | 45 | 61 | 101 | 796 | 824 | 101 | 7 | 0.193 |
| other_ambiguous | SNV | 35 | 2 | 2 | 0 | 16 | 11 | 2 | 2 | 0.175 |
| other_ambiguous | INDEL | 194 | 40 | 19 | 14 | 35 | 65 | 16 | 5 | 0.165 |

### HPRC frequency per category, element level (min over elements) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | SNV | 2 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0.374 |
| HG008N_present_broad | INDEL | 75 | 0 | 0 | 0 | 1 | 51 | 21 | 2 | 0.392 |
| HG008N_present_rare | SNV | 3 | 0 | 0 | 0 | 0 | 2 | 1 | 0 | 0.282 |
| HG008N_present_rare | INDEL | 30 | 1 | 0 | 4 | 20 | 4 | 0 | 1 | 0.165 |
| HG008N_absent_HPRC_other | SNV | 150 | 0 | 3 | 14 | 34 | 45 | 29 | 25 | 0.32 |
| HG008N_absent_HPRC_other | INDEL | 1935 | 0 | 52 | 104 | 817 | 862 | 94 | 6 | 0.198 |
| other_ambiguous | SNV | 35 | 2 | 1 | 0 | 14 | 14 | 2 | 2 | 0.218 |
| other_ambiguous | INDEL | 194 | 34 | 11 | 12 | 40 | 76 | 16 | 5 | 0.199 |

### HPRC frequency per category, full-graph VCF AF (all called haplotypes, matching record) (union, chr1-22)

| category | kind | n | 0 | (0,0.05) | [0.05,0.1) | [0.1,0.2) | [0.2,0.5) | [0.5,0.8) | [0.8,1] | median |
|---|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | SNV | 2 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0.376 |
| HG008N_present_broad | INDEL | 74 | 0 | 0 | 0 | 7 | 50 | 15 | 2 | 0.343 |
| HG008N_present_rare | SNV | 3 | 0 | 0 | 0 | 0 | 2 | 1 | 0 | 0.315 |
| HG008N_present_rare | INDEL | 24 | 0 | 0 | 2 | 18 | 3 | 0 | 1 | 0.164 |
| HG008N_absent_HPRC_other | SNV | 114 | 0 | 6 | 5 | 27 | 33 | 25 | 18 | 0.335 |
| HG008N_absent_HPRC_other | INDEL | 1695 | 1 | 68 | 52 | 904 | 613 | 55 | 2 | 0.169 |
| other_ambiguous | SNV | 3 | 0 | 1 | 0 | 2 | 0 | 0 | 0 | 0.135 |
| other_ambiguous | INDEL | 34 | 0 | 14 | 8 | 3 | 8 | 1 | 0 | 0.067 |

### Superpopulation of HPRC carriers per category: observed / expected (union, chr1-22)

expected per locus = carriers x the superpopulation's share of the complete haplotypes at that locus (hypergeometric variance for z; ignores the pairing of haplotypes within individuals and linkage between loci, so z overstates significance). HG008 donor: European; HPRC v1.1 has no EUR sample (AFR 23, AMR 16, EAS 4, SAS 1).

| category | level | loci | carrier haps | AFR O/E (z) | AMR O/E (z) | EAS O/E (z) | SAS O/E (z) | CHM13 carries |
|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | exact allele | 77 | 2544 | 0.916 (z -5.9) | 1.098 (z +5.3) | 1.057 (z +1.3) | 1.007 (z +0.1) | 27 |
| HG008N_present_broad | set, any ALT path | 77 | 2610 | 0.917 (z -6.0) | 1.095 (z +5.3) | 1.059 (z +1.4) | 1.030 (z +0.3) | 27 |
| HG008N_present_rare | exact allele | 33 | 300 | 0.912 (z -1.7) | 1.073 (z +1.1) | 1.123 (z +0.7) | 1.296 (z +0.9) | 7 |
| HG008N_present_rare | set, any ALT path | 33 | 552 | 0.998 (z -0.1) | 1.025 (z +0.6) | 0.975 (z -0.2) | 0.744 (z -1.3) | 13 |
| HG008N_absent_HPRC_other | exact allele | 2085 | 28208 | 1.107 (z +21.4) | 0.891 (z -16.7) | 0.885 (z -7.4) | 0.917 (z -2.5) | 298 |
| HG008N_absent_HPRC_other | set, any ALT path | 2085 | 36096 | 1.090 (z +21.1) | 0.910 (z -16.5) | 0.908 (z -7.1) | 0.916 (z -3.1) | 408 |
| other_ambiguous | exact allele | 229 | 238 | 0.792 (z -3.9) | 1.223 (z +3.4) | 1.042 (z +0.3) | 1.515 (z +1.5) | 6 |
| other_ambiguous | set, any ALT path | 229 | 3030 | 0.908 (z -6.7) | 1.110 (z +6.5) | 1.045 (z +1.1) | 0.896 (z -1.2) | 55 |
| all | exact allele | 2424 | 31290 | 1.088 (z +18.5) | 0.912 (z -14.4) | 0.902 (z -6.7) | 0.932 (z -2.2) | 338 |
| all | set, any ALT path | 2424 | 42288 | 1.065 (z +16.8) | 0.938 (z -12.6) | 0.928 (z -6.1) | 0.919 (z -3.2) | 503 |

### HPRC support x category (union, chr1-22)

| hprc_support | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| element_set_other_allele | 0 | 0 | 0 | 1 | 76 | 101 | 31 | 134 |
| exact_allele | 2 | 75 | 3 | 26 | 65 | 1789 | 2 | 16 |
| none | 0 | 0 | 0 | 1 | 0 | 0 | 13 | 77 |
| recombinant_pieces | 0 | 0 | 0 | 2 | 9 | 45 | 0 | 10 |

### exact-allele frequency class x category (union, chr1-22)

| exact_freq_class | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 0 | 0 | 9 | 43 |
| 0.2-0.5 | 2 | 54 | 0 | 0 | 9 | 714 | 1 | 5 |
| <0.2 | 0 | 0 | 3 | 30 | 139 | 1166 | 34 | 187 |
| >=0.5 | 0 | 21 | 0 | 0 | 2 | 55 | 0 | 2 |
| NA | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 |

### element-level frequency class x category (union, chr1-22)

| elem_freq_class | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| 0.2-0.5 | 2 | 51 | 2 | 4 | 45 | 862 | 14 | 76 |
| <0.2 | 0 | 1 | 0 | 25 | 51 | 973 | 17 | 97 |
| >=0.5 | 0 | 23 | 1 | 1 | 54 | 100 | 4 | 21 |
| NA | 0 | 0 | 0 | 0 | 0 | 0 | 11 | 43 |

### frequency class agrees across element / set / exact / full-VCF x category (union, chr1-22)

| freq_class_agree | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 8 | 3 | 7 | 91 | 379 | 17 | 92 |
| True | 2 | 67 | 0 | 23 | 59 | 1556 | 29 | 145 |

### Event-hap pattern x category (union, chr1-22)

| event_hap_pattern | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| both | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| event_hap_only | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| neither | 0 | 0 | 0 | 0 | 150 | 1935 | 43 | 221 |
| no_event_hap | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| other_hap_only | 2 | 75 | 1 | 22 | 0 | 0 | 0 | 6 |
| other_hap_only:patient_frame | 0 | 0 | 2 | 7 | 0 | 0 | 0 | 1 |
| unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 7 |

### HG008-N carries_alt hap1/hap2 (array) x category (union, chr1-22)

| normal_evidence | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| ALT/ALT | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| ALT/no | 1 | 36 | 0 | 12 | 0 | 0 | 0 | 5 |
| no/ALT | 1 | 39 | 1 | 11 | 0 | 0 | 0 | 2 |
| no/no | 0 | 0 | 2 | 7 | 150 | 1935 | 43 | 222 |
| no/unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 2 |
| unresolved/no | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
| unresolved/unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 4 |

### HG008-N narrow (s6/s7 event window) call hap1/hap2 x category (union, chr1-22)

| normal_narrow_evidence | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| ALT/ALT | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| ALT/no | 1 | 36 | 0 | 12 | 1 | 1 | 1 | 6 |
| ALT/unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| ALT_plus_other/ALT_plus_other | 0 | 0 | 0 | 0 | 3 | 1 | 0 | 0 |
| ALT_plus_other/no | 0 | 0 | 1 | 4 | 11 | 5 | 1 | 0 |
| no/ALT | 1 | 39 | 1 | 12 | 0 | 5 | 0 | 2 |
| no/ALT_plus_other | 0 | 0 | 1 | 1 | 8 | 2 | 0 | 1 |
| no/no | 0 | 0 | 0 | 1 | 127 | 1904 | 42 | 215 |
| no/unresolved | 0 | 0 | 0 | 0 | 0 | 7 | 1 | 4 |
| unresolved/ALT | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
| unresolved/no | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 1 |
| unresolved/unresolved | 0 | 0 | 0 | 0 | 0 | 5 | 0 | 5 |

### HG008-N allele class (locus) x category (union, chr1-22)

| normal_allele_class | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| SNV_site_deleted_on_event_hap | 0 | 0 | 0 | 0 | 30 | 0 | 8 | 0 |
| both_REF | 0 | 0 | 0 | 0 | 21 | 656 | 1 | 3 |
| carries_ALT | 2 | 75 | 3 | 30 | 0 | 0 | 0 | 9 |
| germline_STR_other_length | 0 | 0 | 0 | 0 | 91 | 1179 | 30 | 196 |
| other | 0 | 0 | 0 | 0 | 8 | 100 | 4 | 22 |
| unresolved | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 7 |

### truth INFO event applied to the event hap (s8b) x category (union, chr1-22)

| T_info_class | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| NA | 0 | 0 | 0 | 1 | 5 | 12 | 1 | 11 |
| eq_ALT | 0 | 4 | 0 | 1 | 45 | 1642 | 1 | 76 |
| eq_ALT=other_hap | 2 | 71 | 1 | 21 | 0 | 0 | 0 | 0 |
| eq_other_hap | 0 | 0 | 2 | 7 | 0 | 0 | 0 | 1 |
| event_outside_stretch | 0 | 0 | 0 | 0 | 10 | 10 | 3 | 4 |
| novel | 0 | 0 | 0 | 0 | 90 | 271 | 41 | 145 |

### HG008-T contigs vs normal haps over the array (s8b) x category (union, chr1-22)

| tumor_verdict | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| NA | 0 | 0 | 0 | 0 | 1 | 1 | 2 | 10 |
| all_eq_hap1 | 1 | 33 | 1 | 16 | 0 | 1 | 0 | 0 |
| all_eq_hap2 | 1 | 34 | 2 | 13 | 0 | 0 | 0 | 2 |
| no_change_both | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| no_change_hom | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| novel_ALT | 0 | 4 | 0 | 1 | 45 | 1652 | 1 | 76 |
| novel_other | 0 | 4 | 0 | 0 | 104 | 279 | 43 | 147 |

### HG008-T any contig == ALT over the array x category (union, chr1-22)

| tumor_any_ALT_array | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 9 |
| False | 0 | 0 | 2 | 7 | 104 | 282 | 44 | 146 |
| True | 2 | 75 | 1 | 23 | 45 | 1652 | 1 | 82 |

### HG008-T shows no change at the array x category (union, chr1-22)

| tumor_shows_no_change | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 2 | 75 | 3 | 30 | 0 | 0 | 0 | 0 |
| False | 0 | 0 | 0 | 0 | 150 | 1932 | 46 | 233 |
| True | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 4 |

### kind of the truth INFO event (normal-assembly frame) x category (union, chr1-22)

| asm_event_kind | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| DEL | 0 | 41 | 0 | 15 | 21 | 746 | 21 | 81 |
| INS | 0 | 34 | 0 | 11 | 43 | 1110 | 18 | 113 |
| SNV | 2 | 0 | 3 | 3 | 82 | 68 | 6 | 34 |
| none | 0 | 0 | 0 | 1 | 4 | 11 | 1 | 9 |

### truth-spelling reads on a perfect platform x category (union, chr1-22)

| truth_read_support | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| no_truth_spelling_read | 0 | 0 | 0 | 4 | 48 | 79 | 40 | 209 |
| yes | 2 | 75 | 3 | 26 | 102 | 1856 | 6 | 28 |

### Dipcall relation x category (union, chr1-22)

| dipcall_relation | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| germline INDEL within 10 bp | 0 | 0 | 0 | 0 | 45 | 3 | 17 | 0 |
| germline SNV within 10 bp | 0 | 0 | 0 | 0 | 10 | 54 | 13 | 18 |
| germline variant at the same site (other allele) | 0 | 0 | 0 | 0 | 40 | 1195 | 10 | 201 |
| identical germline allele | 2 | 75 | 3 | 30 | 21 | 5 | 4 | 9 |
| none within 10 bp | 0 | 0 | 0 | 0 | 34 | 678 | 2 | 9 |

### in_dip_bed x category (union, chr1-22)

| in_dip_bed | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 5 |
| True | 2 | 75 | 3 | 30 | 150 | 1932 | 46 | 232 |

### in GIAB v0.2 nogermlineoverlap BED x category (union, chr1-22)

| in_nogermline_bed | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 32 | 1 | 16 | 86 | 1095 | 43 | 231 |
| True | 2 | 43 | 2 | 14 | 64 | 840 | 3 | 6 |

### PoN tagged (repo rule) x category (union, chr1-22)

| PoN_tagged | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| False | 0 | 0 | 0 | 1 | 0 | 20 | 2 | 41 |
| True | 2 | 75 | 3 | 29 | 150 | 1915 | 44 | 196 |

### PoN rule tag x category (union, chr1-22)

| PoN_rule | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| 1000G+CoLoRSdb | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| CoLoRSdb | 0 | 9 | 0 | 6 | 15 | 867 | 5 | 108 |
| dbSNP | 0 | 0 | 0 | 0 | 0 | 1 | 3 | 6 |
| dbSNP+1000G+CoLoRSdb | 0 | 2 | 0 | 0 | 1 | 6 | 0 | 0 |
| dbSNP+CoLoRSdb | 0 | 61 | 3 | 19 | 39 | 826 | 6 | 13 |
| gnomAD | 0 | 0 | 0 | 0 | 3 | 7 | 3 | 7 |
| gnomAD+1000G+CoLoRSdb | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| gnomAD+CoLoRSdb | 0 | 1 | 0 | 1 | 14 | 112 | 12 | 46 |
| gnomAD+dbSNP | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 |
| gnomAD+dbSNP+1000G+CoLoRSdb | 1 | 0 | 0 | 0 | 3 | 0 | 0 | 0 |
| gnomAD+dbSNP+CoLoRSdb | 1 | 2 | 0 | 3 | 74 | 93 | 15 | 16 |
| none | 0 | 0 | 0 | 1 | 0 | 20 | 2 | 41 |

### element source x category (union, chr1-22)

| element_source | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| graph | 2 | 75 | 3 | 30 | 150 | 1935 | 7 | 50 |
| graph_closest_no_reads | 0 | 0 | 0 | 0 | 0 | 0 | 4 | 0 |
| reads | 0 | 0 | 0 | 0 | 0 | 0 | 35 | 187 |

### closest: what the reads walk x category (union, chr1-22)

| closest_kind | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 2 | 75 | 3 | 30 | 150 | 1935 | 7 | 51 |
| read_allele_nearer_truth | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 21 |
| read_allele_not_truth | 0 | 0 | 0 | 0 | 0 | 0 | 32 | 165 |
| s2_primary:read_allele_not_truth | 0 | 0 | 0 | 0 | 0 | 0 | 4 | 0 |

### with_germline: germline alleles on the event hap x category (union, chr1-22)

| germline_phase | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
|  | 2 | 75 | 3 | 29 | 138 | 1883 | 43 | 217 |
| consistent | 0 | 0 | 0 | 1 | 12 | 52 | 0 | 3 |
| off_event_hap | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 17 |

### path choice (read counts dropped) x category (union, chr1-22)

| path_choice | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| closest_primary | 0 | 0 | 0 | 0 | 0 | 0 | 4 | 0 |
| no_truth_read_support:alternative:phase | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 |
| no_truth_read_support:primary | 0 | 0 | 0 | 2 | 13 | 19 | 0 | 11 |
| pooled_majority | 0 | 0 | 0 | 0 | 0 | 0 | 35 | 187 |
| single_path | 2 | 72 | 3 | 28 | 121 | 1804 | 7 | 34 |
| truth_read_supported:alternative | 0 | 1 | 0 | 0 | 2 | 54 | 0 | 2 |
| truth_read_supported:alternative:tie | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| truth_read_supported:primary | 0 | 2 | 0 | 0 | 13 | 51 | 0 | 2 |
| truth_read_supported:primary:tie | 0 | 0 | 0 | 0 | 1 | 4 | 0 | 1 |

### s2 match x category (union, chr1-22)

| match | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| closest | 0 | 0 | 0 | 0 | 0 | 0 | 39 | 187 |
| exact | 2 | 75 | 3 | 29 | 138 | 1883 | 4 | 30 |
| with_germline | 0 | 0 | 0 | 1 | 12 | 52 | 3 | 20 |

### PoN: % tagged per category vs the non-absorbed truths (chr1-22)

repo rule = allele match, gnomAD / CoLoRSdb AF >= 1e-4, dbSNP non-somatic, 1000G (SNV; the repo INDEL path has no PoN). pop AF = max(gnomAD, CoLoRSdb). The graph-specific loss = absorbed truths a linear tumor-only caller with that filter would still report.

| group | kind | n | repo rule % | pop AF >= 1e-3 % | pop AF >= 0.01 % | pop AF >= 0.05 % |
|---|---|---|---|---|---|---|
| HG008N_present_broad | SNV | 2 | 100.0 | 100.0 | 100.0 | 100.0 |
| HG008N_present_rare | SNV | 3 | 100.0 | 100.0 | 100.0 | 100.0 |
| HG008N_absent_HPRC_other | SNV | 150 | 100.0 | 100.0 | 96.0 | 74.7 |
| other_ambiguous | SNV | 46 | 95.7 | 87.0 | 69.6 | 45.7 |
| all absorbed | SNV | 201 | 99.0 | 97.0 | 90.0 | 68.7 |
| not absorbed (other HG008T truths) | SNV | 8489 | 8.4 | 5.0 | 2.7 | 1.3 |
| HG008N_present_broad | INDEL | 75 | 100.0 | 100.0 | 100.0 | 100.0 |
| HG008N_present_rare | INDEL | 30 | 96.7 | 96.7 | 96.7 | 90.0 |
| HG008N_absent_HPRC_other | INDEL | 1935 | 99.0 | 98.7 | 93.5 | 69.9 |
| other_ambiguous | INDEL | 237 | 82.7 | 78.9 | 51.1 | 17.7 |
| all absorbed | INDEL | 2277 | 97.3 | 96.7 | 89.3 | 65.7 |
| not absorbed (other HG008T truths) | INDEL | 6219 | 78.2 | 74.1 | 45.3 | 12.0 |

### PoN rule and population AF per category (union, chr1-22)

| category | kind | n | PoN tagged | gnomAD AF n | gnomAD median | gnomAD >=1e-4 | CoLoRSdb AF n | CoLoRSdb median | CoLoRSdb >=1e-4 |
|---|---|---|---|---|---|---|---|---|---|
| HG008N_present_broad | SNV | 2 | 2 | 2 | 0.322 | 2 | 2 | 0.2862 | 2 |
| HG008N_present_broad | INDEL | 75 | 75 | 3 | 0.371 | 3 | 75 | 0.3186 | 75 |
| HG008N_present_rare | SNV | 3 | 3 | 0 |  | 0 | 3 | 0.1227 | 3 |
| HG008N_present_rare | INDEL | 30 | 29 | 4 | 0.1245 | 4 | 29 | 0.1408 | 29 |
| HG008N_absent_HPRC_other | SNV | 150 | 150 | 96 | 0.094 | 95 | 146 | 0.1146 | 146 |
| HG008N_absent_HPRC_other | INDEL | 1935 | 1915 | 215 | 0.065 | 214 | 1906 | 0.0874 | 1906 |
| other_ambiguous | SNV | 46 | 44 | 31 | 0.024 | 30 | 38 | 0.0353 | 38 |
| other_ambiguous | INDEL | 237 | 196 | 70 | 0.0089 | 69 | 183 | 0.0152 | 183 |

### Read-vs-graph agreement per platform (perfect loci)

| platform | kind | perfect loci | with read data | with no-edit spanning ALT reads | with truth-spelling reads | majority of all reads spells germline only | no truth read can test an element | truth majority walk = graph elements | different | median frac truth reads containing the set | frac >= 0.9 | frac < 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PacBio | SNV | 170 | 170 | 166 | 93 | 16 | 0 | 55 | 38 | 1.0 | 89 | 2 |
| PacBio | INDEL | 1363 | 1363 | 1363 | 1174 | 184 | 0 | 1084 | 90 | 1.0 | 1123 | 20 |
| ONT | SNV | 148 | 148 | 148 | 85 | 14 | 0 | 59 | 26 | 1.0 | 84 | 0 |
| ONT | INDEL | 989 | 989 | 989 | 797 | 190 | 0 | 724 | 73 | 1.0 | 754 | 16 |
| Illumina | SNV | 122 | 122 | 121 | 78 | 8 | 0 | 51 | 27 | 1.0 | 76 | 1 |
| Illumina | INDEL | 1916 | 1916 | 1912 | 1724 | 160 | 0 | 1642 | 82 | 1.0 | 1701 | 19 |

### HPRC carriers: GFA allele level vs VCF carriers_any (on GFA-complete haplotypes; loci with somatic elements)

| VCF | kind | loci | with matching record | identical carrier sets | identical, with record | hap-level disagreements |
|---|---|---|---|---|---|---|
| d9 | SNV | 192 | 124 | 72 | 27 | vcf_only:gfa_set_other_allele 2518; gfa_only:no_matching_record 154; gfa_only:vcf_other_allele 15; vcf_only:gfa_other_walk 4 |
| d9 | INDEL | 2234 | 1756 | 1667 | 1486 | vcf_only:gfa_set_other_allele 3501; gfa_only:no_matching_record 2992; vcf_only:gfa_other_walk 309; gfa_only:vcf_other_allele 173; gfa_only:vcf_gt_missing 1 |
| full | SNV | 192 | 122 | 73 | 25 | vcf_only:gfa_set_other_allele 2419; gfa_only:no_matching_record 150; gfa_only:vcf_other_allele 19; vcf_only:gfa_other_walk 4 |
| full | INDEL | 2234 | 1818 | 1688 | 1522 | vcf_only:gfa_set_other_allele 3217; gfa_only:no_matching_record 2521; gfa_only:vcf_other_allele 487; vcf_only:gfa_other_walk 342; gfa_only:vcf_gt_missing 3 |

### d9 frequency-filter floor (somatic elements)

node coverage = HPRC + CHM13 + GRCh38 haplotypes whose s4 rows (complete or partial path) visit the branch node (min over the nodes of the element); a lower bound (s4 partial paths stop 500 nodes from an anchor). full_window = loci where all 88 HPRC haplotypes have a complete traversal, where the counts are exact. A whole-GFA scan (audit hprc_membership a5) found >= 7 HPRC haplotypes for every branch node not on a reference path.

| measure | value |
|---|---|
| full_windows | 114 |
| full_window_elements | 120 |
| full_window_branch_node_coverage_min | 1 |
| full_window_node_coverage_hist | {1: 1, 2: 1, 6: 1, 8: 1, 9: 5, 10: 8, 11: 9, 12: 4, 13: 4, 14: 6, 15: 40} |
| full_window_hprc_carriers_min | 0 |
| full_window_hprc_carriers_hist | {0: 1, 1: 1, 5: 2, 7: 1, 9: 8, 10: 10, 11: 12, 12: 5, 13: 8, 14: 8, 15: 64} |
| full_window_not_on_CHM13_hprc_carriers_min | 9 |
| full_window_on_CHM13_below9 | 5 |
| all_not_on_CHM13_node_coverage_lt9 | 17 |
| n_elements | 3368 |
| branch_elements_node_coverage_min | 1 |
| node_coverage_lt9 | 60 |
| node_coverage_hist | {1: 5, 2: 4, 3: 5, 4: 7, 5: 5, 6: 12, 7: 6, 8: 16, 9: 70, 10: 64, 11: 73, 12: 61, 13: 61, 14: 54, 15: 1889} |
| complete_carriers_min | 0 |
| complete_carriers_lt9 | 542 |
| any_row_carriers_min | 0 |
| any_row_carriers_lt9 | 200 |

### Read data per platform

| platform | loci with read data | loci with no-edit spanning ALT reads | perfect loci | perfect loci with read data |
|---|---|---|---|---|
| PacBio | 2478 | 2459 | 1533 | 1533 |
| ONT | 1365 | 1361 | 1137 | 1137 |
| Illumina | 2478 | 2351 | 2038 | 2038 |
