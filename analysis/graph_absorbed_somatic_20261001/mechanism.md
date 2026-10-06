# Mechanism tables (s11; categories from s9)

### Patient-frame event class x category (GRCh38 kind in the column; union, chr1-22)

|  | HG008N_present_broad SNV | HG008N_present_broad INDEL | HG008N_present_rare SNV | HG008N_present_rare INDEL | HG008N_absent_HPRC_other SNV | HG008N_absent_HPRC_other INDEL | other_ambiguous SNV | other_ambiguous INDEL |
|---|---|---|---|---|---|---|---|---|
| HP>=7 1-unit | 0 | 67 | 0 | 21 | 53 | 1612 | 32 | 131 |
| STR 1-unit | 0 | 7 | 0 | 3 | 9 | 209 | 7 | 51 |
| repeat substitution | 1 | 0 | 3 | 3 | 75 | 67 | 5 | 34 |
| repeat multi-unit / other | 0 | 1 | 0 | 1 | 2 | 35 | 0 | 9 |
| no patient-frame event | 0 | 0 | 0 | 1 | 4 | 11 | 1 | 9 |
| non-repeat SNV (CpG Ti) | 1 | 0 | 0 | 0 | 6 | 0 | 1 | 0 |
| non-repeat SNV | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 |
| in_STR(p4) 1-unit | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| VNTR>6 1-unit | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| in_STR(p2) 1-unit | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| HP4-6 1-unit | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

### Absorption rate by patient-frame context, all HG008T truth alleles on chr1-22

|  | truth alleles | absorbed | % |
|---|---|---|---|
| HP 1-unit, tract 7-9 | 459 | 21 | 4.6 |
| HP 1-unit, tract 10-14 | 1792 | 393 | 21.9 |
| HP 1-unit, tract 15-19 | 1858 | 543 | 29.2 |
| HP 1-unit, tract 20-29 | 2136 | 830 | 38.9 |
| HP 1-unit, tract >=30 | 524 | 129 | 24.6 |
| HP4-6 1-unit | 196 | 1 | 0.5 |
| STR 1-unit | 765 | 286 | 37.4 |
| repeat multi-unit / other | 250 | 48 | 19.2 |
| repeat substitution | 1844 | 188 | 10.2 |
| non-repeat INDEL | 450 | 0 | 0.0 |
| non-repeat SNV (CpG Ti) | 730 | 8 | 1.1 |
| non-repeat SNV | 6074 | 2 | 0.0 |
| VNTR>6 1-unit | 21 | 1 | 4.8 |
| in_STR(p1) 1-unit | 1 | 0 | 0.0 |
| in_STR(p2) 1-unit | 7 | 1 | 14.3 |
| in_STR(p3) 1-unit | 3 | 0 | 0.0 |
| in_STR(p4) 1-unit | 6 | 1 | 16.7 |
| in_STR(p5) 1-unit | 8 | 0 | 0.0 |
| in_STR(p6) 1-unit | 3 | 0 | 0.0 |
| no patient-frame event | 59 | 26 | 44.1 |

### Local LOH around the locus (Illumina allele balance at phased het SNVs +-20 kb) x category (union, chr1-22)

|  | HG008N_present_broad | HG008N_present_rare | HG008N_absent_HPRC_other | other_ambiguous |
|---|---|---|---|---|
| O_lost | 17 | 12 | 722 | 88 |
| both_retained | 59 | 17 | 1110 | 159 |
| imbalanced | 1 | 0 | 4 | 2 |
| no_event_hap | 0 | 1 | 15 | 10 |
| no_hets | 0 | 3 | 234 | 24 |

Controls (800 random non-absorbed truths): both_retained 463, O_lost 243, no_hets 86, imbalanced 6, no_event_hap 1, E_lost 1
