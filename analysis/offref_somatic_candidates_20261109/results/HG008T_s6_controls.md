# HG008T: s6 rules applied to the on-reference controls (keys)

| group (truth) | kind | keys | offref_somatic_t1 | offref_somatic_t2 | offref_germline | offref_non | offref_non_normal_alt | uncertain |
|---|---|---|---|---|---|---|---|---|
| ctrl_somatic | SNV | 1656 | 1326 (80.1 %) | 267 (16.1 %) | 4 (0.2 %) | 2 (0.1 %) | 1 (0.1 %) | 56 (3.4 %) |
| ctrl_somatic | INDEL | 1988 | 217 (10.9 %) | 1391 (70.0 %) | 11 (0.6 %) | 104 (5.2 %) | 6 (0.3 %) | 259 (13.0 %) |
| ctrl_germline | SNV | 2586 | 0 (0.0 %) | 0 (0.0 %) | 2236 (86.5 %) | 14 (0.5 %) | 10 (0.4 %) | 326 (12.6 %) |
| ctrl_germline | INDEL | 2606 | 0 (0.0 %) | 9 (0.3 %) | 1335 (51.2 %) | 428 (16.4 %) | 13 (0.5 %) | 821 (31.5 %) |
| ctrl_non | SNV | 1317 | 0 (0.0 %) | 0 (0.0 %) | 51 (3.9 %) | 661 (50.2 %) | 94 (7.1 %) | 511 (38.8 %) |
| ctrl_non | INDEL | 1317 | 0 (0.0 %) | 4 (0.3 %) | 172 (13.1 %) | 646 (49.1 %) | 52 (3.9 %) | 443 (33.6 %) |

Control keys get the platforms of every target set with a tensor of them, as off-reference keys do. The snvaf05-only SNV keys (AF 0.05-0.08, COLO829T) have no controls of their own.

Uncertain, by note:

- ctrl_somatic: tumor_only_unconfirmed 23, neither_multiplatform 8, tumor_only_normal_alt 5, no_query 2, asm_normal_other_chrom 2, asm_normal_reads_absent 2, asm_normal_multicopy 1, asm_normal_no_read_check 1, reads_only 1
- ctrl_germline: asm_normal_multicopy 230, neither_multiplatform 49, asm_normal_other_chrom 14, tumor_only_unconfirmed 8, asm_normal_reads_absent 7, no_query 6, tumor_only_normal_alt 5, asm_normal_no_read_check 5, neither_supported 1
- ctrl_non: asm_normal_multicopy 298, asm_normal_other_chrom 93, asm_normal_no_read_check 42, neither_multiplatform 35, asm_normal_reads_absent 20, no_query 13, neither_supported 4, reads_only 4, tumor_only_unconfirmed 4, tumor_only_normal_alt 2
