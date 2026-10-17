# INDEL call sets under the 7fc7396 and the 71cb51a INDEL PoN rule

(`collect.py`; chr1 after the PoN. 7fc7396 rule: gnomAD / CoLoRSdb exact allele at AF >= 0.01 when the ALT changes a homopolymer of >= 7 bp, else AF >= 1e-3; 71cb51a rule: >= 8 bp, else AF >= 1e-4. Same calls and PASS threshold; best = the chr1 curve's best point. * = the run's own vcf.sh already had the 71cb51a rule.)


## in the benchmark BED

| sample | group | run (set) | PASS F1 7fc7396 -> 71cb51a | best F1 7fc7396 -> 71cb51a | best P / R (71cb51a) | truth tagged 7fc7396 -> 71cb51a |
|---|---|---|---|---|---|---|
| HG008T | --norm dataset, d46 | HG008_Illumina_INDEL_d46_normds | 0.244 -> 0.253 | 0.249 -> 0.256 | 0.287 / 0.232 | 413 -> 418 |
| HG008T | d46 | HG008_Illumina_INDEL_d46_7fc7396 | 0.244 -> 0.253 | 0.249 -> 0.256 | 0.287 / 0.232 | 413 -> 418 |
| HG008T | d46 | HG008_ONT_INDEL_d46_7fc7396 | 0.015 -> 0.015 | 0.155 -> 0.159 | 0.210 / 0.128 | 413 -> 418 |
| HG008T | d46 | HG008_PacBio_INDEL_d46_7fc7396 | 0.095 -> 0.098 | 0.269 -> 0.275 | 0.377 / 0.216 | 413 -> 418 |
| HG008T | d9 | HG008_Illumina_INDEL_7fc7396 | 0.241 -> 0.251 | 0.264 -> 0.269 | 0.367 / 0.213 | 413 -> 418 |
| HG008T | d9 | HG008_ONT_INDEL_7fc7396 | 0.015 -> 0.015 | 0.140 -> 0.143 | 0.224 / 0.106 | 413 -> 418 |
| HG008T | d9 | HG008_PacBio_INDEL_7fc7396 | 0.089 -> 0.092 | 0.260 -> 0.263 | 0.368 / 0.204 | 413 -> 418 |
| HG008T | d9, 2026-10-03 h90 | HG008_Illumina_INDEL_h90 | no PoN -> 0.252 | no PoN -> 0.265 | 0.414 / 0.194 | no PoN -> 418 |
| HG008T | d9, 2026-10-03 h90 | HG008_ONT_INDEL_h90 | no PoN -> 0.012 | no PoN -> 0.133 | 0.176 / 0.107 | no PoN -> 418 |
| HG008T | d9, 2026-10-03 h90 | HG008_PacBio_INDEL_h90 | no PoN -> 0.107 | no PoN -> 0.244 | 0.337 / 0.191 | no PoN -> 418 |
| HG008T | d9, 45d56ad lr 4e-4 | HG008_Illumina_INDEL_small | 0.246 -> 0.255 | 0.263 -> 0.272 | 0.339 / 0.227 | 413 -> 418 |
| HG008T | d9, 45d56ad lr 4e-4 | HG008_ONT_INDEL_small | 0.016 -> 0.017 | 0.147 -> 0.149 | 0.243 / 0.107 | 413 -> 418 |
| HG008T | d9, 45d56ad lr 4e-4 | HG008_PacBio_INDEL_small | 0.097 -> 0.100 | 0.231 -> 0.234 | 0.337 / 0.179 | 413 -> 418 |
| HG008T | d9, repro | HG008_Illumina_INDEL_7fc7396_repro | 0.241 -> 0.251 | 0.264 -> 0.269 | 0.367 / 0.213 | 413 -> 418 |
| HG008T | linear | HG008_Illumina_INDEL_linear_7fc7396 | 0.133 -> 0.146 | 0.199 -> 0.201 | 0.228 / 0.181 | 413 -> 418 |
| HG008T | linear | HG008_ONT_INDEL_linear_7fc7396 | 0.022 -> 0.023 | 0.082 -> 0.084 | 0.074 / 0.097 | 413 -> 418 |
| HG008T | linear | HG008_PacBio_INDEL_linear_7fc7396 | 0.125 -> 0.133 | 0.173 -> 0.176 | 0.160 / 0.196 | 413 -> 418 |
| HG008T | no ch6, d46 | HG008_Illumina_INDEL_nopc_d46_7fc7396 | 0.215 -> 0.223 | 0.246 -> 0.252 | 0.316 / 0.209 | 413 -> 418 |
| HG008T | no ch6, d46 | HG008_ONT_INDEL_nopc_d46_7fc7396 | 0.014 -> 0.014 | 0.154 -> 0.158 | 0.179 / 0.141 | 413 -> 418 |
| HG008T | no ch6, d46 | HG008_PacBio_INDEL_nopc_d46_7fc7396 | 0.102 -> 0.105 | 0.240 -> 0.245 | 0.292 / 0.211 | 413 -> 418 |
| HG008T | no ch6, d9 | HG008_Illumina_INDEL_nopc_7fc7396 | 0.216 -> 0.225 | 0.255 -> 0.259 | 0.377 / 0.198 | 413 -> 418 |
| HG008T | no ch6, d9 | HG008_ONT_INDEL_nopc_7fc7396 | 0.013 -> 0.013 | 0.131 -> 0.133 | 0.209 / 0.097 | 413 -> 418 |
| HG008T | no ch6, d9 | HG008_PacBio_INDEL_nopc_7fc7396 | 0.104 -> 0.108 | 0.226 -> 0.229 | 0.295 / 0.187 | 413 -> 418 |
| COLO829T | --norm dataset, d46 | COLO829_Illumina_INDEL_d46_normds | 0.108 -> 0.116 | 0.167 -> 0.172 | 0.196 / 0.153 | 81 -> 81 |
| COLO829T | d46 | COLO829_Illumina_INDEL_d46_7fc7396 | 0.109 -> 0.117 | 0.163 -> 0.168 | 0.187 / 0.153 | 81 -> 81 |
| COLO829T | d46 | COLO829_ONT_INDEL_d46_7fc7396 | 0.002 -> 0.002 | 0.085 -> 0.089 | 0.094 / 0.084 | 81 -> 81 |
| COLO829T | d46 | COLO829_fiberseq_INDEL_d46_7fc7396 | 0.020 -> 0.021 | 0.129 -> 0.137 | 0.105 / 0.199 | 81 -> 81 |
| COLO829T | d9 | COLO829_Illumina_INDEL_7fc7396 | 0.097 -> 0.104 | 0.153 -> 0.161 | 0.135 / 0.199 | 81 -> 81 |
| COLO829T | d9 | COLO829_ONT_INDEL_7fc7396 | 0.003 -> 0.003 | 0.088 -> 0.095 | 0.085 / 0.107 | 81 -> 81 |
| COLO829T | d9 | COLO829_fiberseq_INDEL_7fc7396 | 0.023 -> 0.024 | 0.135 -> 0.142 | 0.116 / 0.183 | 81 -> 81 |
| COLO829T | d9, 2026-10-03 h90 | COLO829_Illumina_INDEL_h90 | no PoN -> 0.089 | no PoN -> 0.163 | 0.200 / 0.137 | no PoN -> 81 |
| COLO829T | d9, 2026-10-03 h90 | COLO829_ONT_INDEL_h90 | no PoN -> 0.002 | no PoN -> 0.057 | 0.054 / 0.061 | no PoN -> 81 |
| COLO829T | d9, 2026-10-03 h90 | COLO829_fiberseq_INDEL_h90 | no PoN -> 0.008 | no PoN -> 0.144 | 0.265 / 0.099 | no PoN -> 81 |
| COLO829T | d9, 45d56ad lr 4e-4 | COLO829_Illumina_INDEL_small | 0.115 -> 0.125 | 0.195 -> 0.202 | 0.215 / 0.191 | 81 -> 81 |
| COLO829T | d9, 45d56ad lr 4e-4 | COLO829_ONT_INDEL_small | 0.003 -> 0.003 | 0.071 -> 0.074 | 0.188 / 0.046 | 81 -> 81 |
| COLO829T | d9, 45d56ad lr 4e-4 | COLO829_fiberseq_INDEL_small | 0.016 -> 0.017 | 0.120 -> 0.124 | 0.167 / 0.099 | 81 -> 81 |
| COLO829T | linear | COLO829_Illumina_INDEL_linear_7fc7396 | 0.032 -> 0.034 | 0.069 -> 0.070 | 0.056 / 0.092 | 81 -> 81 |
| COLO829T | linear | COLO829_ONT_INDEL_linear_7fc7396 | 0.004 -> 0.004 | 0.019 -> 0.020 | 0.015 / 0.030 | 81 -> 81 |
| COLO829T | linear | COLO829_fiberseq_INDEL_linear_7fc7396 | 0.040 -> 0.044 | 0.092 -> 0.095 | 0.091 / 0.099 | 81 -> 81 |
| COLO829T | no ch6, d46 | COLO829_Illumina_INDEL_nopc_d46_7fc7396 | 0.105 -> 0.112 | 0.176 -> 0.178 | 0.232 / 0.145 | 81 -> 81 |
| COLO829T | no ch6, d46 | COLO829_ONT_INDEL_nopc_d46_7fc7396 * | 0.002 -> 0.002 | 0.102 -> 0.102 | 0.097 / 0.107 | 81 -> 81 |
| COLO829T | no ch6, d46 | COLO829_fiberseq_INDEL_nopc_d46_7fc7396 | 0.017 -> 0.017 | 0.082 -> 0.086 | 0.069 / 0.115 | 81 -> 81 |
| COLO829T | no ch6, d9 | COLO829_Illumina_INDEL_nopc_7fc7396 | 0.091 -> 0.096 | 0.166 -> 0.172 | 0.231 / 0.137 | 81 -> 81 |
| COLO829T | no ch6, d9 | COLO829_ONT_INDEL_nopc_7fc7396 * | 0.002 -> 0.002 | 0.082 -> 0.082 | 0.079 / 0.084 | 81 -> 81 |
| COLO829T | no ch6, d9 | COLO829_fiberseq_INDEL_nopc_7fc7396 | 0.015 -> 0.016 | 0.108 -> 0.113 | 0.131 / 0.099 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina) | 0.121 -> 0.129 | 0.163 -> 0.167 | 0.167 / 0.168 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina_P25) | 0.192 -> 0.195 | 0.233 -> 0.233 | 0.594 / 0.145 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina_P50) | 0.232 -> 0.238 | 0.330 -> 0.330 | 0.544 / 0.237 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina_P75) | 0.191 -> 0.201 | 0.308 -> 0.312 | 0.527 / 0.221 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT) | 0.001 -> 0.001 | 0.054 -> 0.055 | 0.096 / 0.038 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT_P25) | 0.001 -> 0.001 | 0.060 -> 0.062 | 0.167 / 0.038 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT_P50) | 0.001 -> 0.001 | 0.068 -> 0.070 | 0.146 / 0.046 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT_P75) | 0.001 -> 0.001 | 0.051 -> 0.051 | 0.148 / 0.030 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq) | 0.007 -> 0.007 | 0.068 -> 0.069 | 0.080 / 0.061 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq_P25) | 0.003 -> 0.003 | 0.073 -> 0.074 | 0.119 / 0.053 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq_P50) | 0.004 -> 0.004 | 0.053 -> 0.055 | 0.050 / 0.061 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq_P75) | 0.004 -> 0.005 | 0.047 -> 0.049 | 0.030 / 0.130 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina) | 0.074 -> 0.077 | 0.174 -> 0.182 | 0.210 / 0.160 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina_P25) | 0.035 -> 0.037 | 0.132 -> 0.133 | 0.500 / 0.076 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina_P50) | 0.054 -> 0.057 | 0.125 -> 0.133 | 0.113 / 0.160 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina_P75) | 0.066 -> 0.070 | 0.136 -> 0.146 | 0.147 / 0.145 | 81 -> 81 |

## whole chr1, no BED

| sample | group | run (set) | PASS F1 7fc7396 -> 71cb51a | best F1 7fc7396 -> 71cb51a | best P / R (71cb51a) | truth tagged 7fc7396 -> 71cb51a |
|---|---|---|---|---|---|---|
| HG008T | --norm dataset, d46 | HG008_Illumina_INDEL_d46_normds | 0.235 -> 0.242 | 0.240 -> 0.245 | 0.312 / 0.202 | 413 -> 418 |
| HG008T | d46 | HG008_Illumina_INDEL_d46_7fc7396 | 0.235 -> 0.243 | 0.239 -> 0.244 | 0.314 / 0.200 | 413 -> 418 |
| HG008T | d46 | HG008_ONT_INDEL_d46_7fc7396 | 0.013 -> 0.013 | 0.153 -> 0.156 | 0.209 / 0.124 | 413 -> 418 |
| HG008T | d46 | HG008_PacBio_INDEL_d46_7fc7396 | 0.056 -> 0.057 | 0.254 -> 0.260 | 0.341 / 0.209 | 413 -> 418 |
| HG008T | d9 | HG008_Illumina_INDEL_7fc7396 | 0.233 -> 0.242 | 0.252 -> 0.257 | 0.342 / 0.206 | 413 -> 418 |
| HG008T | d9 | HG008_ONT_INDEL_7fc7396 | 0.015 -> 0.015 | 0.134 -> 0.137 | 0.215 / 0.101 | 413 -> 418 |
| HG008T | d9 | HG008_PacBio_INDEL_7fc7396 | 0.082 -> 0.084 | 0.247 -> 0.250 | 0.338 / 0.198 | 413 -> 418 |
| HG008T | d9, 2026-10-03 h90 | HG008_Illumina_INDEL_h90 | no PoN -> 0.241 | no PoN -> 0.255 | 0.393 / 0.189 | no PoN -> 418 |
| HG008T | d9, 2026-10-03 h90 | HG008_ONT_INDEL_h90 | no PoN -> 0.011 | no PoN -> 0.130 | 0.173 / 0.104 | no PoN -> 418 |
| HG008T | d9, 2026-10-03 h90 | HG008_PacBio_INDEL_h90 | no PoN -> 0.102 | no PoN -> 0.235 | 0.314 / 0.187 | no PoN -> 418 |
| HG008T | d9, 45d56ad lr 4e-4 | HG008_Illumina_INDEL_small | 0.237 -> 0.245 | 0.253 -> 0.261 | 0.318 / 0.221 | 413 -> 418 |
| HG008T | d9, 45d56ad lr 4e-4 | HG008_ONT_INDEL_small | 0.015 -> 0.015 | 0.144 -> 0.145 | 0.243 / 0.104 | 413 -> 418 |
| HG008T | d9, 45d56ad lr 4e-4 | HG008_PacBio_INDEL_small | 0.089 -> 0.092 | 0.224 -> 0.227 | 0.311 / 0.178 | 413 -> 418 |
| HG008T | d9, repro | HG008_Illumina_INDEL_7fc7396_repro | 0.233 -> 0.242 | 0.252 -> 0.257 | 0.342 / 0.206 | 413 -> 418 |
| HG008T | linear | HG008_Illumina_INDEL_linear_7fc7396 | 0.127 -> 0.140 | 0.192 -> 0.194 | 0.211 / 0.179 | 413 -> 418 |
| HG008T | linear | HG008_ONT_INDEL_linear_7fc7396 | 0.020 -> 0.020 | 0.077 -> 0.079 | 0.070 / 0.091 | 413 -> 418 |
| HG008T | linear | HG008_PacBio_INDEL_linear_7fc7396 | 0.121 -> 0.128 | 0.163 -> 0.166 | 0.166 / 0.165 | 413 -> 418 |
| HG008T | no ch6, d46 | HG008_Illumina_INDEL_nopc_d46_7fc7396 | 0.204 -> 0.211 | 0.234 -> 0.239 | 0.288 / 0.205 | 413 -> 418 |
| HG008T | no ch6, d46 | HG008_ONT_INDEL_nopc_d46_7fc7396 | 0.012 -> 0.012 | 0.147 -> 0.150 | 0.172 / 0.134 | 413 -> 418 |
| HG008T | no ch6, d46 | HG008_PacBio_INDEL_nopc_d46_7fc7396 | 0.094 -> 0.097 | 0.227 -> 0.231 | 0.265 / 0.205 | 413 -> 418 |
| HG008T | no ch6, d9 | HG008_Illumina_INDEL_nopc_7fc7396 | 0.181 -> 0.188 | 0.247 -> 0.251 | 0.356 / 0.194 | 413 -> 418 |
| HG008T | no ch6, d9 | HG008_ONT_INDEL_nopc_7fc7396 | 0.013 -> 0.013 | 0.127 -> 0.128 | 0.207 / 0.093 | 413 -> 418 |
| HG008T | no ch6, d9 | HG008_PacBio_INDEL_nopc_7fc7396 | 0.099 -> 0.103 | 0.214 -> 0.217 | 0.263 / 0.184 | 413 -> 418 |
| COLO829T | --norm dataset, d46 | COLO829_Illumina_INDEL_d46_normds | 0.108 -> 0.116 | 0.167 -> 0.172 | 0.196 / 0.153 | 81 -> 81 |
| COLO829T | d46 | COLO829_Illumina_INDEL_d46_7fc7396 | 0.109 -> 0.117 | 0.163 -> 0.168 | 0.187 / 0.153 | 81 -> 81 |
| COLO829T | d46 | COLO829_ONT_INDEL_d46_7fc7396 | 0.002 -> 0.002 | 0.085 -> 0.089 | 0.094 / 0.084 | 81 -> 81 |
| COLO829T | d46 | COLO829_fiberseq_INDEL_d46_7fc7396 | 0.020 -> 0.021 | 0.129 -> 0.137 | 0.105 / 0.199 | 81 -> 81 |
| COLO829T | d9 | COLO829_Illumina_INDEL_7fc7396 | 0.097 -> 0.104 | 0.153 -> 0.161 | 0.135 / 0.199 | 81 -> 81 |
| COLO829T | d9 | COLO829_ONT_INDEL_7fc7396 | 0.003 -> 0.003 | 0.088 -> 0.095 | 0.085 / 0.107 | 81 -> 81 |
| COLO829T | d9 | COLO829_fiberseq_INDEL_7fc7396 | 0.023 -> 0.024 | 0.135 -> 0.142 | 0.116 / 0.183 | 81 -> 81 |
| COLO829T | d9, 2026-10-03 h90 | COLO829_Illumina_INDEL_h90 | no PoN -> 0.089 | no PoN -> 0.163 | 0.200 / 0.137 | no PoN -> 81 |
| COLO829T | d9, 2026-10-03 h90 | COLO829_ONT_INDEL_h90 | no PoN -> 0.002 | no PoN -> 0.057 | 0.054 / 0.061 | no PoN -> 81 |
| COLO829T | d9, 2026-10-03 h90 | COLO829_fiberseq_INDEL_h90 | no PoN -> 0.008 | no PoN -> 0.144 | 0.265 / 0.099 | no PoN -> 81 |
| COLO829T | d9, 45d56ad lr 4e-4 | COLO829_Illumina_INDEL_small | 0.115 -> 0.125 | 0.195 -> 0.202 | 0.215 / 0.191 | 81 -> 81 |
| COLO829T | d9, 45d56ad lr 4e-4 | COLO829_ONT_INDEL_small | 0.003 -> 0.003 | 0.071 -> 0.074 | 0.188 / 0.046 | 81 -> 81 |
| COLO829T | d9, 45d56ad lr 4e-4 | COLO829_fiberseq_INDEL_small | 0.016 -> 0.017 | 0.120 -> 0.124 | 0.167 / 0.099 | 81 -> 81 |
| COLO829T | linear | COLO829_Illumina_INDEL_linear_7fc7396 | 0.032 -> 0.034 | 0.069 -> 0.070 | 0.056 / 0.092 | 81 -> 81 |
| COLO829T | linear | COLO829_ONT_INDEL_linear_7fc7396 | 0.004 -> 0.004 | 0.019 -> 0.020 | 0.015 / 0.030 | 81 -> 81 |
| COLO829T | linear | COLO829_fiberseq_INDEL_linear_7fc7396 | 0.040 -> 0.044 | 0.092 -> 0.095 | 0.091 / 0.099 | 81 -> 81 |
| COLO829T | no ch6, d46 | COLO829_Illumina_INDEL_nopc_d46_7fc7396 | 0.105 -> 0.112 | 0.176 -> 0.178 | 0.232 / 0.145 | 81 -> 81 |
| COLO829T | no ch6, d46 | COLO829_ONT_INDEL_nopc_d46_7fc7396 * | 0.002 -> 0.002 | 0.102 -> 0.102 | 0.097 / 0.107 | 81 -> 81 |
| COLO829T | no ch6, d46 | COLO829_fiberseq_INDEL_nopc_d46_7fc7396 | 0.017 -> 0.017 | 0.082 -> 0.086 | 0.069 / 0.115 | 81 -> 81 |
| COLO829T | no ch6, d9 | COLO829_Illumina_INDEL_nopc_7fc7396 | 0.091 -> 0.096 | 0.166 -> 0.172 | 0.231 / 0.137 | 81 -> 81 |
| COLO829T | no ch6, d9 | COLO829_ONT_INDEL_nopc_7fc7396 * | 0.002 -> 0.002 | 0.082 -> 0.082 | 0.079 / 0.084 | 81 -> 81 |
| COLO829T | no ch6, d9 | COLO829_fiberseq_INDEL_nopc_7fc7396 | 0.015 -> 0.016 | 0.108 -> 0.113 | 0.131 / 0.099 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina) | 0.121 -> 0.129 | 0.163 -> 0.167 | 0.167 / 0.168 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina_P25) | 0.192 -> 0.195 | 0.233 -> 0.233 | 0.594 / 0.145 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina_P50) | 0.232 -> 0.238 | 0.330 -> 0.330 | 0.544 / 0.237 | 81 -> 81 |
| COLO829T | purity joint | COLO829_Illumina_INDEL_joint_7fc7396 (COLO829T_Illumina_P75) | 0.191 -> 0.201 | 0.308 -> 0.312 | 0.527 / 0.221 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT) | 0.001 -> 0.001 | 0.054 -> 0.055 | 0.096 / 0.038 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT_P25) | 0.001 -> 0.001 | 0.060 -> 0.062 | 0.167 / 0.038 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT_P50) | 0.001 -> 0.001 | 0.068 -> 0.070 | 0.146 / 0.046 | 81 -> 81 |
| COLO829T | purity joint | COLO829_ONT_INDEL_joint_7fc7396 (COLO829T_ONT_P75) | 0.001 -> 0.001 | 0.051 -> 0.051 | 0.148 / 0.030 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq) | 0.007 -> 0.007 | 0.068 -> 0.069 | 0.080 / 0.061 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq_P25) | 0.003 -> 0.003 | 0.073 -> 0.074 | 0.119 / 0.053 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq_P50) | 0.004 -> 0.004 | 0.053 -> 0.055 | 0.050 / 0.061 | 81 -> 81 |
| COLO829T | purity joint | COLO829_fiberseq_INDEL_joint_7fc7396 (COLO829T_fiberseq_P75) | 0.004 -> 0.005 | 0.047 -> 0.049 | 0.030 / 0.130 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina) | 0.074 -> 0.077 | 0.174 -> 0.182 | 0.210 / 0.160 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina_P25) | 0.035 -> 0.037 | 0.132 -> 0.133 | 0.500 / 0.076 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina_P50) | 0.054 -> 0.057 | 0.125 -> 0.133 | 0.113 / 0.160 | 81 -> 81 |
| COLO829T | purity joint, BQ per read | COLO829_Illumina_INDEL_joint_readnorm (COLO829T_Illumina_P75) | 0.066 -> 0.070 | 0.136 -> 0.146 | 0.147 / 0.145 | 81 -> 81 |
