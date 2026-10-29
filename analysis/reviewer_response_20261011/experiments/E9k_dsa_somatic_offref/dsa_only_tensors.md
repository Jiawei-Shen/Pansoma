# E9k part 4: tensors for the DSA somatic SNVs without a GRCh38 position (dsa_only_tensors.py)

Events: DSA somatic SNVs of class dsa_only (no GRCh38 position in the chain) and grch38_diff_ref (on a germline allele),
placed on each graph with E9j's method (giraffe of the DSA windows, s3_classify.py). Tensor = a tensor of the set whose
representative allele is the event's key (exact), or a tensor at the same site with another representative allele.

## d9

| event class | events | A (GRCh38 node) | B (off-reference node) | C (graph branch) | D (unmapped / MAPQ < 20 / same) |
|---|---:|---:|---:|---:|---:|
| dsa_only | 5,352 | 2,134 | 183 | 43 | 2,992 |
| grch38_diff_ref | 148 | 18 | 30 | 48 | 52 |

| set | A or B SNP keys | tensor, exact A1 | tensor at the site, other A1 | of them B (off-reference) exact |
|---|---:|---:|---:|---:|
| COLO829T_Illumina | 2,264 | 462 | 32 | 155 |
| COLO829T_fiberseq | 2,264 | 306 | 22 | 149 |
| COLO829T_ONT | 2,264 | 302 | 31 | 148 |

## d46

| event class | events | A (GRCh38 node) | B (off-reference node) | C (graph branch) | D (unmapped / MAPQ < 20 / same) |
|---|---:|---:|---:|---:|---:|
| dsa_only | 5,352 | 1,834 | 422 | 41 | 3,055 |
| grch38_diff_ref | 148 | 13 | 29 | 52 | 54 |

| set | A or B SNP keys | tensor, exact A1 | tensor at the site, other A1 | of them B (off-reference) exact |
|---|---:|---:|---:|---:|
| COLO829T_Illumina_d46 | 2,200 | 425 | 21 | 166 |
| COLO829T_fiberseq_d46 | 2,200 | 347 | 45 | 177 |
| COLO829T_ONT_d46 | 2,200 | 312 | 36 | 162 |

