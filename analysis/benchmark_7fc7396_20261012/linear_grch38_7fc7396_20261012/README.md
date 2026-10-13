# GRCh38-only (linear) graph: HG008T and COLO829T, three platforms, SNV and INDEL, chr1 test

2026-10-12, on RIS compute1, pipeline of 7fc7396 (the v2 README's SNV recipe and INDEL recipe). The same samples,
platforms, recipes and test as the HPRC v1.1 d9 runs (`../hg008t_snv_indel_7fc7396_20261012` and the earlier
COLO829T runs) and the HPRC v2.1 d46 runs (`../hprc_v21_d46_7fc7396_20261012`), on the tensor sets built against a
graph of GRCh38 alone (`<SET>_linear`, `analysis/graph_ablation_tensors_20261005`): the no-pangenome control of the
graph comparison. One model per set, trained on chr2-22, tested on chr1, PASS at the validation recall-0.9 threshold,
the PoN of each kind, rtg vcfeval in the BED and without it.

Files: `collect.py` (standard library; writes `tables.md` with linear / d9 / d46 side by side after the PoN, before it,
and per tensor, `summary.json`, and `runs/<run>/` copies of `train.log` and the rtg reports), `tables.md`.

## Data

`Dataset/pansoma_v2_tensors/<SET>_linear/tensors` (build of 9d6bd8d on the GRCh38-only graph `hg38_linear.giraffe`:
96,872,672 nodes, `--haplotypes 17`, ch6 `int8-count-haplotypes100`, build AF 0.08; labels `5430533aa100`, the same
truth files as the d9 and d46 sets). Checked on 2026-10-12: every shard of every manifest present with the size its
tensor count gives, labels and summaries present. Without population alleles in the graph, every germline variant of
the sample becomes a candidate, so the SNV sets are 2-3x the d9 sets.

| set | SNV tensors (somatic) | INDEL tensors (somatic) |
|---|---:|---:|
| HG008T_Illumina_linear | 5,854,337 (10,711) | 913,786 (7,259) |
| HG008T_PacBio_linear | 3,977,546 (8,962) | 2,444,557 (7,861) |
| HG008T_ONT_linear | 5,182,542 (9,003) | 2,926,299 (6,573) |
| COLO829T_Illumina_linear | 5,007,009 (38,697) | 918,053 (1,859) |
| COLO829T_fiberseq_linear | 4,525,734 (38,066) | 2,229,342 (2,001) |
| COLO829T_ONT_linear | 5,690,176 (37,805) | 4,709,597 (1,829) |

**Graph index.** `graph_vcf` and `linear_vcf` need `hg38_linear.giraffe.graph_index.sqlite` and
`hg38_linear.giraffe.grch38_path/` (the sets' manifests: `/scratch/jshen/data/Linear_Reference_VG_Indexes/graph_index/`
on the old cluster). Copied by the user to `Dataset/Linear_Reference_VG_Indexes/graph_index/` (2026-10-12, complete by
22:09), checked once unchanged for 10 min with `ris/tools/check_graph_index.py` (job 47013): 96,872,672 nodes in the
sqlite, as the linear manifests record; the reference path's arrays consistent with its `meta.json` (195 contigs,
3,315 walks). Linked into `Dataset/pansoma_v2_tensors/graph_index/`, where `ris/vcf.sh` takes it by name.

## What was run

Per set and kind: `ris/prebuild.sh` (CPU), then `ris/run.sh` (train with the recipe, predict chr1) on the `epigenome`
queue's H100s (SNV on 4 under DDP with `--batch-size 1024` per GPU, INDEL on 1), after every COLO829T purity and d46
job had started; then `ris/vcf.sh` (submitted once the graph index was on RIS, each after its GPU job).

```bash
cd /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma
SETS="HG008T_Illumina HG008T_PacBio HG008T_ONT COLO829T_Illumina COLO829T_fiberseq COLO829T_ONT"
for k in SNV INDEL; do for b in $SETS; do s=${b}_linear
  ris/submit.sh -g 0 -n 8 -m 64 pnv2_prebuild_${s}_$k prebuild.sh $s $k; done; done
for k in SNV INDEL; do for b in $SETS; do s=${b}_linear; samp=${b%%T_*}; plat=${b#*T_}
  if [ $k = SNV ]; then g=4 sl=16 m=320 vw=3; else g=1 sl=5 m=256 vw=6; fi
  NOVCF=1 KIND=$k GPUS=$g SLOTS=$sl MEM=$m VAL_WORKERS=$vw WAIT="ended(pnv2_prebuild_${s}_$k)" \
      ris/snv_pipeline.sh ${samp}_${plat}_${k}_linear_7fc7396 $s
done; done
```

| run | GPU job | prebuild job | vcf job |
|---|---|---|---|
| `HG008_Illumina_SNV_linear_7fc7396` | 46916 | 46904 | 47035 |
| `HG008_PacBio_SNV_linear_7fc7396` | 46917 | 46905 | 47036 |
| `HG008_ONT_SNV_linear_7fc7396` | 46918 | 46906 | 47037 |
| `COLO829_Illumina_SNV_linear_7fc7396` | 46919 | 46907 | 47038 |
| `COLO829_fiberseq_SNV_linear_7fc7396` | 46920 | 46908 | 47039 |
| `COLO829_ONT_SNV_linear_7fc7396` | 46921 | 46909 | 47040 |
| `HG008_Illumina_INDEL_linear_7fc7396` | 46922 | 46910 | 47041 |
| `HG008_PacBio_INDEL_linear_7fc7396` | 46924 | 46911 | 47042 |
| `HG008_ONT_INDEL_linear_7fc7396` | 46925 | 46912 | 47043 |
| `COLO829_Illumina_INDEL_linear_7fc7396` | 46926 | 46913 | 47044 |
| `COLO829_fiberseq_INDEL_linear_7fc7396` | 46927 | 46914 | 47045 |
| `COLO829_ONT_INDEL_linear_7fc7396` | 46928 | 46915 | 47046 |

**Comparability.** The linear and d46 runs use the same settings (SNV on 4 GPUs, INDEL on 1, lr 1e-4); the d9 runs
differ for some sets (`../hprc_v21_d46_7fc7396_20261012/README.md`: SNV on 1 GPU except HG008T Illumina, COLO829T INDEL
at lr 4e-4), so linear vs d46 is the cleanest of the three comparisons.

## Results (chr1, all 12 runs)

All numbers: `tables.md`. In the BED. SNV rows: PASS F1 at each model's validation recall-0.9 threshold (the chr1
curve's best F1 in brackets). INDEL rows: the chr1 curve's best F1 after the PoN (PASS F1 in brackets), because the
recall-0.9 threshold lands far below the best point for INDELs; the best point is chosen on the test chromosome, an
upper bound for a threshold rule, not a held-out result. The raw and no-BED columns are PASS F1 for SNV and best F1 for
INDEL:

| sample | platform | kind | after the PoN: linear / d9 / d46 | before the PoN (raw): linear / d9 / d46 | without the BED, after the PoN: linear / d9 / d46 |
|---|---|---|---|---|---|
| HG008T | Illumina | SNV | **0.469** (0.544) / 0.585 / 0.570 | 0.009 / 0.188 / 0.295 | 0.197 / 0.577 / 0.561 |
| HG008T | PacBio | SNV | **0.513** (0.596) / 0.621 / 0.605 | 0.009 / 0.191 / 0.310 | 0.359 / 0.614 / 0.600 |
| HG008T | ONT | SNV | **0.537** (0.570) / 0.596 / 0.575 | 0.008 / 0.182 / 0.287 | 0.365 / 0.586 / 0.563 |
| HG008T | Illumina | INDEL | **0.199** (0.133) / 0.264 (0.241) / 0.249 (0.244) | 0.048 / 0.147 / 0.173 | 0.192 / 0.252 / 0.239 |
| HG008T | PacBio | INDEL | **0.173** (0.125) / 0.260 (0.089) / 0.269 (0.095) | 0.056 / 0.150 / 0.159 | 0.163 / 0.247 / 0.254 |
| HG008T | ONT | INDEL | **0.082** (0.022) / 0.140 (0.015) / 0.155 (0.015) | 0.030 / 0.090 / 0.104 | 0.077 / 0.134 / 0.153 |
| COLO829T | Illumina | SNV | **0.468** (0.554) / 0.778 / 0.783 | 0.073 / 0.571 / 0.685 | 0.468 / 0.778 / 0.783 |
| COLO829T | fiberseq | SNV | **0.680** (0.682) / 0.781 / 0.779 | 0.075 / 0.555 / 0.683 | 0.680 / 0.781 / 0.779 |
| COLO829T | ONT | SNV | **0.418** (0.511) / 0.706 / 0.737 | 0.066 / 0.494 / 0.626 | 0.418 / 0.706 / 0.737 |
| COLO829T | Illumina | INDEL | **0.069** (0.032) / 0.195 (0.115) / 0.163 (0.109) | 0.071 / 0.125 / 0.105 | 0.069 / 0.195 / 0.163 |
| COLO829T | fiberseq | INDEL | **0.092** (0.040) / 0.120 (0.016) / 0.129 (0.020) | 0.027 / 0.069 / 0.081 | 0.092 / 0.120 / 0.129 |
| COLO829T | ONT | INDEL | **0.019** (0.004) / 0.071 (0.003) / 0.085 (0.002) | 0.018 / 0.051 / 0.056 | 0.019 / 0.071 / 0.085 |

**SNV: the pangenome matters, also after the PoN.** The linear graph is below both HPRC graphs on every SNV set after
the PoN: HG008T -0.06 to -0.12, COLO829T -0.10 to -0.31 against d9. Without population alleles in the graph every germline
variant of the sample is a candidate: chr1 test tensors 2.2-7.6x d9's (HG008T PacBio 251,303 vs 32,880), per-tensor somatic
AP 0.012-0.016 against 0.20-0.23 (HG008T), raw PASS F1 0.01 (HG008T) / 0.07 (COLO829T). The PoN removes most of the
germline calls but not all: false calls after it 1.4-1.7x d9's on HG008T, 1.8-7.8x on COLO829T (Illumina 2,216 vs
284). Outside the HG008T benchmark BED it is worse still (0.20-0.37 without the BED against 0.58-0.61 for d9).

**INDEL: the linear graph keeps more truth INDELs, and loses them to germline.** Without HPRC alleles in the graph, the
somatic INDELs that the pangenome absorbs (`analysis/graph_absorbed_somatic_20261001`) get a candidate: truth INDELs
with a tensor 0.83 / 0.87 / 0.75 on HG008T Illumina / PacBio / ONT against 0.54 / 0.70 / 0.69 (d9) and 0.61 / 0.76 / 0.72
(d46); the rtg curve ceilings before the PoN 0.80 / 0.84 / 0.54 against 0.40 / 0.47 / 0.42 (d9). Part survives the PoN
(ceilings after it 0.284 / 0.301 / 0.238 against 0.257 / 0.259 / 0.238; COLO829T fiberseq / ONT 0.298 / 0.282 against
0.244 / 0.237). Compared by the best point of the chr1 curve after the PoN, linear is below both HPRC graphs on every
INDEL set: HG008T 0.199 / 0.173 / 0.082 against 0.264 / 0.260 / 0.140 (d9) and 0.249 / 0.269 / 0.155 (d46) on Illumina /
PacBio / ONT; COLO829T 0.069 / 0.092 / 0.019 against 0.195 / 0.120 / 0.071 (d9) and 0.163 / 0.129 / 0.085 (d46). At its
best point linear has 1.6-2.9x d9's false calls on HG008T (1.9-8.9x on COLO829T) for no more true calls (HG008T PacBio
97 TP / 437 FP against 120 / 217 on d9; ONT 57 / 747 against 64 / 262). Only at the PASS threshold does linear score higher on the long reads (HG008T PacBio
0.125 vs 0.089), because its recall-0.9 threshold happens to land at fewer false calls than d9's (1,703 vs 2,630 at
about the same 150 true calls), not because it separates somatic from germline better.

**So**: linear < d9 = d46 for SNV (PASS F1) and for INDEL (best F1), after the PoN; before it linear << d9 < d46. For
INDEL the linear graph gives more candidates and a slightly higher ceiling after the PoN, but ranks them worse. Training: the linear
SNV sets are 2.9-4.6 M training tensors (HG008T 7.4 k somatic), 76-266 min on 4 H100s; validation somatic AP 0.011-0.012
on HG008T (0.17-0.26 COLO829T).

