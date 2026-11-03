# HCC1395T chr1 SNVs: Pansoma six-tumor false negatives and recall by read allele fraction

```
Illumina multicancer6 truth 3440 TP 2827 FN 613 threshold 0.735
  no tensor: filtered:min_af: 244
  tensor, below threshold (somatic is the top class): 166
  tensor, below threshold (germline is the top class): 90
  tensor, population allele (PoN): 53
  no tensor: filtered:min_variants: 25
  no tensor: no_candidate: 23
  tensor, P_SOMATIC < 0.01 (not in VCF): 12
  TP AF median 0.36 [IQR 0.23-0.54] DP median 108.00 [IQR 90.00-126.00]
  FN below thr AF median 0.66 [IQR 0.51-0.98] DP median 113.00 [IQR 92.00-128.00] P_SOMATIC median 0.58 [IQR 0.44-0.68]
  AF [0,0.15): TP 370, FN below thr 1 (germline > somatic in 0); recall of tensor truth 1.00
  AF [0.15,0.25): TP 409, FN below thr 5 (germline > somatic in 5); recall of tensor truth 0.99
  AF [0.25,0.4): TP 850, FN below thr 28 (germline > somatic in 8); recall of tensor truth 0.97
  AF [0.4,0.6): TP 671, FN below thr 64 (germline > somatic in 16); recall of tensor truth 0.91
  AF [0.6,0.85): TP 371, FN below thr 77 (germline > somatic in 21); recall of tensor truth 0.83
  AF [0.85,1.01): TP 156, FN below thr 81 (germline > somatic in 40); recall of tensor truth 0.66
```

Illumina: truth SNVs with a Pansoma tensor in the VCF, by Pansoma read AF; recall per caller
| AF | truth | Pansoma | DeepSomatic-TO | ClairS-TO |
| [0.08, 0.25) | 796 | 0.979 | 0.985 | 0.835 |
| [0.25, 0.4) | 896 | 0.949 | 0.978 | 0.953 |
| [0.4, 0.6) | 749 | 0.896 | 0.960 | 0.940 |
| [0.6, 0.85) | 453 | 0.819 | 0.951 | 0.967 |
| [0.85, 1.01) | 242 | 0.645 | 0.930 | 0.971 |
| no tensor in the VCF | 304 | 0.000 | 0.510 | 0.197 |

PacBio: truth SNVs with a Pansoma tensor in the VCF, by Pansoma read AF; recall per caller
| AF | truth | Pansoma | DeepSomatic-TO | ClairS-TO |
| [0.08, 0.25) | 673 | 0.961 | 0.979 | 0.854 |
| [0.25, 0.4) | 775 | 0.952 | 0.970 | 0.957 |
| [0.4, 0.6) | 881 | 0.952 | 0.966 | 0.949 |
| [0.6, 0.85) | 427 | 0.836 | 0.963 | 0.939 |
| [0.85, 1.01) | 251 | 0.801 | 0.888 | 0.944 |
| no tensor in the VCF | 433 | 0.000 | 0.679 | 0.088 |

ONT: truth SNVs with a Pansoma tensor in the VCF, by Pansoma read AF; recall per caller
| AF | truth | Pansoma | DeepSomatic-TO | ClairS-TO |
| [0.08, 0.25) | 318 | 0.925 | 0.959 | 0.909 |
| [0.25, 0.4) | 886 | 0.953 | 0.970 | 0.960 |
| [0.4, 0.6) | 832 | 0.917 | 0.972 | 0.954 |
| [0.6, 0.85) | 350 | 0.851 | 0.960 | 0.937 |
| [0.85, 1.01) | 246 | 0.805 | 0.841 | 0.943 |
| no tensor in the VCF | 808 | 0.000 | 0.146 | 0.022 |

