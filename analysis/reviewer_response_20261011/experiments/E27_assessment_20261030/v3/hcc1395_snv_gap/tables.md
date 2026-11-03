# HCC1395T chr1 SNVs: where Pansoma's recall is lost relative to the comparators

| chr1 SNVs in the BED | Illumina | PacBio | ONT |
|---|---|---|---|
| truth SNVs | 3440 | 3440 | 3440 |
| tensor (scored by the model) | 3148 (91.5%) | 3024 (87.9%) | 2652 (77.1%) |
| candidate below allele fraction 0.08 | 244 (7.1%) | 204 (5.9%) | 227 (6.6%) |
| candidate with fewer than 3 ALT reads | 25 (0.7%) | 143 (4.2%) | 356 (10.3%) |
| no candidate | 23 (0.7%) | 69 (2.0%) | 205 (6.0%) |
| recall: six-tumor Pansoma / DeepSomatic-TO | 0.822 / 0.927 | 0.809 / 0.927 | 0.697 / 0.766 |
| found by DeepSomatic-TO only / by Pansoma only | 391 / 28 | 450 / 42 | 290 / 53 |
| — DeepSomatic-TO only: tensor (scored by the model) | 241 | 167 | 179 |
| — DeepSomatic-TO only: candidate below allele fraction 0.08 | 141 | 201 | 101 |
| — DeepSomatic-TO only: candidate with fewer than 3 ALT reads | 0 | 71 | 9 |
| — DeepSomatic-TO only: no candidate | 9 | 11 | 1 |
| recall: six-tumor Pansoma / ClairS-TO | 0.822 / 0.859 | 0.809 / 0.822 | 0.697 / 0.730 |
| found by ClairS-TO only / by Pansoma only | 301 / 172 | 193 / 146 | 181 / 66 |
| — ClairS-TO only: tensor (scored by the model) | 241 | 159 | 167 |
| — ClairS-TO only: candidate below allele fraction 0.08 | 55 | 27 | 13 |
| — ClairS-TO only: candidate with fewer than 3 ALT reads | 0 | 3 | 0 |
| — ClairS-TO only: no candidate | 5 | 4 | 1 |
