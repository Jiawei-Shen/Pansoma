# COLO829T dilutions on chr1, before calling (`purity.py`)

Truth: chr1 somatic truth alleles that pass and lie in the BED; with a tensor = recall ceiling of any model. AF: median (quartiles) of the read AF of the chr1 tensors labelled somatic / germline; measured purity = the somatic SNV median relative to the platform's P100 (biased upward at low purity: tensors need AF >= 0.08).

| platform | purity | SNV truth | SNV with a tensor (ceiling) | SNV filtered / no candidate | INDEL truth | INDEL with a tensor (ceiling) | somatic SNV tensors: AF median (q1-q3) | measured purity | germline SNV AF median | somatic INDEL tensors: AF median |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Illumina | P100 | 2443 | 2008 (0.822) | 384 / 51 | 131 | 57 (0.435) | 0.973 (0.493-1.000), n 2,043 | 1.00 | 0.990 | 0.412 (n 81) |
| Illumina | P75 | 2443 | 1936 (0.792) | 438 / 69 | 131 | 55 (0.420) | 0.640 (0.404-0.747), n 1,968 | 0.66 | 0.789 | 0.317 (n 80) |
| Illumina | P50 | 2443 | 1864 (0.763) | 468 / 111 | 131 | 53 (0.405) | 0.387 (0.298-0.508), n 1,891 | 0.40 | 0.636 | 0.240 (n 71) |
| Illumina | P25 | 2443 | 1748 (0.716) | 482 / 213 | 131 | 36 (0.275) | 0.191 (0.150-0.264), n 1,770 | 0.20 | 0.532 | 0.155 (n 51) |
| fiberseq | P100 | 2443 | 2006 (0.821) | 397 / 40 | 131 | 67 (0.511) | 0.961 (0.493-0.988), n 2,022 | 1.00 | 0.986 | 0.450 (n 97) |
| fiberseq | P75 | 2443 | 1940 (0.794) | 412 / 91 | 131 | 65 (0.496) | 0.632 (0.402-0.742), n 1,956 | 0.66 | 0.785 | 0.373 (n 92) |
| fiberseq | P50 | 2443 | 1870 (0.765) | 444 / 129 | 131 | 61 (0.466) | 0.385 (0.291-0.512), n 1,886 | 0.40 | 0.645 | 0.293 (n 86) |
| fiberseq | P25 | 2443 | 1737 (0.711) | 489 / 217 | 131 | 58 (0.443) | 0.188 (0.150-0.270), n 1,751 | 0.20 | 0.542 | 0.189 (n 79) |
| ONT | P100 | 2443 | 1998 (0.818) | 393 / 52 | 131 | 69 (0.527) | 0.823 (0.449-0.902), n 2,017 | 1.00 | 0.865 | 0.236 (n 110) |
| ONT | P75 | 2443 | 1922 (0.787) | 462 / 59 | 131 | 69 (0.527) | 0.576 (0.375-0.681), n 1,941 | 0.70 | 0.667 | 0.228 (n 104) |
| ONT | P50 | 2443 | 1838 (0.752) | 510 / 95 | 131 | 64 (0.488) | 0.360 (0.277-0.466), n 1,853 | 0.44 | 0.537 | 0.220 (n 90) |
| ONT | P25 | 2443 | 1736 (0.711) | 530 / 177 | 131 | 63 (0.481) | 0.180 (0.143-0.247), n 1,749 | 0.22 | 0.470 | 0.175 (n 77) |
