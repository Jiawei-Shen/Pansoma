# Truth alleles already carried by a matched-normal haplotype (chr1-22, PASS truth; R2.2)

Strict rule (`strict_normal_alt.py`): a normal haplotype, reconstructed from the assembly-based dipcall calls, equals
GRCh38 + the truth allele over the site and any tandem repeat it lies in (array, anchors and flank reversion of
`analysis/graph_absorbed_somatic_colo829t_20261005/c3_normal_frame.py`). Looser rule (`overlap.py`): the split,
left-aligned truth allele equals a dipcall germline allele.

| tumor | truth alleles (SNV / INDEL) | strict: SNV / INDEL (in BED) | strict, normal genotype | looser exact allele: SNV / INDEL | no dipcall sequence on either hap (SNV / INDEL) |
|---|---|---|---|---|---|
| HG008T (HG008-N) | 8,690 / 8,496 | 8 / 136 = 144, 0.8% (8 / 122) | 143 het (other hap REF 89, other length 52, other bases 2), 1 hom | 48 / 161 = 209 | 60 / 63 |
| COLO829T (COLO829BL) | 42,035 / 1,912 | 0 / 45 = 45; 2.4% of INDELs (all in BED) | 45 het (REF 28, other length 15, other bases 2) | 1 / 69 = 70 | 848 / 52 |

Checks:
- COLO829T, against the assembly-first c3 call (`normal_frame.tsv`, label `normal carries ALT`, 55 INDELs): 44 agree; the
  other 11 are outside dipcall's dip.bed (chr1:148 Mb, chr14:18-19 Mb, chr15:22 Mb, chr18:14 Mb; the assembly is placed
  there, dipcall gives no sequence); 1 more here (truth 29389: dipcall ALT, assembly 3 bp longer).
- HG008T, against the 110 graph-absorbed `HG008N_present` truths of `analysis/graph_absorbed_somatic_20261001`: 101 are
  strict carriers; the other 9 differ from GRCh38 + truth elsewhere in the array on both haplotypes (s8b array call
  `other`). 8 further absorbed truths are strict carriers.
