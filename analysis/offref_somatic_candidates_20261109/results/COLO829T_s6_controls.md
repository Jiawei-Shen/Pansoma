# COLO829T: s6 rules applied to the on-reference controls (keys)

| group (truth) | kind | keys | offref_somatic_t1 | offref_somatic_t2 | offref_germline | offref_non | offref_non_normal_alt | uncertain |
|---|---|---|---|---|---|---|---|---|
| ctrl_somatic | SNV | 2329 | 1697 (72.9 %) | 343 (14.7 %) | 0 (0.0 %) | 2 (0.1 %) | 0 (0.0 %) | 287 (12.3 %) |
| ctrl_somatic | INDEL | 1187 | 160 (13.5 %) | 673 (56.7 %) | 4 (0.3 %) | 65 (5.5 %) | 0 (0.0 %) | 285 (24.0 %) |
| ctrl_germline | SNV | 2596 | 0 (0.0 %) | 1 (0.0 %) | 2126 (81.9 %) | 15 (0.6 %) | 9 (0.3 %) | 445 (17.1 %) |
| ctrl_germline | INDEL | 2604 | 0 (0.0 %) | 2 (0.1 %) | 1292 (49.6 %) | 370 (14.2 %) | 17 (0.7 %) | 923 (35.4 %) |
| ctrl_non | SNV | 1313 | 9 (0.7 %) | 7 (0.5 %) | 36 (2.7 %) | 432 (32.9 %) | 181 (13.8 %) | 648 (49.4 %) |
| ctrl_non | INDEL | 1318 | 0 (0.0 %) | 4 (0.3 %) | 136 (10.3 %) | 483 (36.6 %) | 103 (7.8 %) | 592 (44.9 %) |

Control keys get the platforms of every target set with a tensor of them, as off-reference keys do. The snvaf05-only SNV keys (AF 0.05-0.08, COLO829T) have no controls of their own.

Uncertain, by note:

- ctrl_somatic: neither_multiplatform 182, reads_only 35, tumor_only_unconfirmed 33, no_query 12, tumor_only_normal_alt 7, asm_normal_multicopy 5, asm_normal_no_read_check 5, asm_normal_reads_absent 4, asm_normal_other_chrom 1, neither_supported 1
- ctrl_germline: neither_multiplatform 647, no_query 100, asm_normal_reads_absent 65, asm_normal_multicopy 46, asm_normal_other_chrom 39, asm_normal_no_read_check 11, tumor_only_unconfirmed 7, tumor_only_normal_alt 6, neither_supported 2
- ctrl_non: neither_multiplatform 306, asm_normal_multicopy 98, asm_normal_other_chrom 74, neither_supported 30, tumor_only_normal_alt 14, asm_normal_no_read_check 11, tumor_only_unconfirmed 6, no_query 5, asm_normal_reads_absent 3, reads_only 2
