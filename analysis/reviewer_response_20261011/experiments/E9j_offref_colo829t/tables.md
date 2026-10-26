# E9j: COLO829T assembly off-reference somatic SNV keys (h50 strict read check)

| set | keys (chr1) | tumor-only, read-supported (chr1) | base / rescue / all_H on them | rescue calls on all keys: tumor-only / all | classes (all keys) |
|---|---|---|---|---|---|
| `COLO829T_Illumina` (d9) | 12 (2) | 7 (0) | 1 / 2 / 4 | 2 / 2 | {'tumor_only': 7, 'neither': 2, 'normal': 3} |
| `COLO829T_fiberseq` (d9) | 8 (0) | 7 (0) | 1 / 4 / 6 | 4 / 4 | {'tumor_only': 7, 'neither': 1} |
| `COLO829T_ONT` (d9) | 9 (1) | 7 (0) | 1 / 3 / 3 | 3 / 3 | {'tumor_only': 7, 'neither': 1, 'normal': 1} |
| `COLO829T_Illumina_d46` (d46) | 16 (6) | 12 (3) | 0 / 7 / 5 | 7 / 8 | {'tumor_only': 12, 'normal': 2, 'neither': 2} |
| `COLO829T_fiberseq_d46` (d46) | 12 (4) | 10 (3) | 0 / 4 / 4 | 4 / 4 | {'tumor_only': 10, 'normal': 1, 'neither': 1} |
| `COLO829T_ONT_d46` (d46) | 11 (3) | 10 (3) | 0 / 4 / 5 | 4 / 4 | {'tumor_only': 10, 'normal': 1} |

base = the set's in-sample SNV model at its PASS threshold (validation recall 0.9); rescue = predict's off-reference rescue (site +- 10 raised to the every-haplotype code, argmax); all_H = ch6 set to H everywhere, at the threshold. All autosomes; chr1 in brackets.
