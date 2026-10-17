# HPRC v2.1 d46 graph: HG008T and COLO829T, three platforms, SNV and INDEL, chr1 test

> 2026-10-16: the INDEL numbers here use the 7fc7396 INDEL PoN rule (homopolymer >= 7 bp: AF >= 0.01, else 1e-3). With
> the 71cb51a rule (>= 8 bp, else 1e-4) every INDEL call set was re-scored: `../indel_pon_71cb51a_20261016/tables.md`
> (+0.000 to +0.013 F1, no comparison changes direction).

2026-10-12, on RIS compute1, pipeline of 7fc7396 (the v2 README's SNV recipe and INDEL recipe). The same samples,
platforms, recipes and test as the HPRC v1.1 d9 runs, on the tensor sets built against the HPRC v2.1 Minigraph-Cactus
GRCh38 graph filtered at d46 (`<SET>_d46`, `analysis/graph_ablation_tensors_20261005`): one model per set, trained on
chr2-22, tested on chr1, PASS at the validation recall-0.9 threshold, the PoN of each kind, rtg vcfeval in the BED and
without it.

Files: `collect.py` (standard library; writes `tables.md`, `summary.json` and `runs/<run>/` copies of `train.log` and the
rtg reports, for the d46 runs and the d9 runs compared), `tables.md`.

## Data

`Dataset/pansoma_v2_tensors/<SET>_d46/tensors` (build of 9d6bd8d on the d46 graph: `--haplotypes 90`, ch6
`int8-count-haplotypes100`, build AF 0.08; labels `5430533aa100`, the same truth files as the d9 sets). Checked on
2026-10-12: every shard of every manifest present with the size its tensor count gives, labels and summaries present.

| set | SNV tensors (somatic) | INDEL tensors (somatic) |
|---|---:|---:|
| HG008T_Illumina_d46 | 2,794,337 (10,221) | 138,850 (5,560) |
| HG008T_PacBio_d46 | 1,294,697 (8,780) | 1,970,759 (6,902) |
| HG008T_ONT_d46 | 1,710,957 (8,787) | 2,524,440 (6,507) |
| COLO829T_Illumina_d46 | 1,654,891 (38,262) | 134,186 (1,369) |
| COLO829T_fiberseq_d46 | 2,407,332 (37,946) | 1,908,317 (1,677) |
| COLO829T_ONT_d46 | 2,348,867 (37,743) | 4,477,120 (1,718) |

**Graph index.** `graph_vcf` (node sequences) and `linear_vcf` (the GRCh38 reference path) need the d46 graph's
`hprc-v2.1-mc-grch38.d46.graph_index.sqlite` and `hprc-v2.1-mc-grch38.d46.grch38_path/` (the sets' manifests:
`/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/` on the old cluster). Copied by the user to
`Dataset/AF-Filtered_VG_Indexes_HPRC_V2/` (scp, 2026-10-12, its `graph_index/` complete by 16:49) and checked once
unchanged for 10 min with `ris/tools/check_graph_index.py` (job 43996): 94,468,093 nodes in the sqlite, as the d46
manifests record; the reference path's arrays consistent with its `meta.json` (195 contigs, 83,303,684 path nodes,
212,465,662 walks). The same check passes on the d9 index (60,118,570 nodes). Both are linked into
`Dataset/pansoma_v2_tensors/graph_index/`, where `ris/vcf.sh` now takes each set's graph index by the names in its
manifests. Training and predict do not need it.

## What was run

Per set and kind: `ris/prebuild.sh` (index and chunk caches, CPU), then `ris/run.sh` (train with the recipe, predict
chr1) on the `epigenome` queue's H100s, then `ris/vcf.sh` (submitted once the graph index was on RIS, each after its GPU job). GPUs: SNV on 4 under DDP with
`--batch-size 1024` per GPU (global batch 4,096), INDEL on 1. The GPU jobs wait until every COLO829T purity job has
started (`started()` dependencies), so they never go before one.

```bash
cd /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma
PAIRS="HG008T_PacBio_d46:SNV HG008T_PacBio_d46:INDEL HG008T_ONT_d46:SNV HG008T_ONT_d46:INDEL HG008T_Illumina_d46:INDEL
       COLO829T_Illumina_d46:SNV COLO829T_Illumina_d46:INDEL COLO829T_fiberseq_d46:SNV COLO829T_fiberseq_d46:INDEL"
for pk in $PAIRS; do s=${pk%%:*} k=${pk#*:}; ris/submit.sh -g 0 -n 8 -m 48 pnv2_prebuild_${s}_$k prebuild.sh $s $k; done
MIX="started(41421) && started(41440) && started(41445)"   # the COLO829T purity GPU jobs still pending then
for pk in $PAIRS; do s=${pk%%:*} k=${pk#*:}; samp=${s%%T_*}; plat=${s#*T_}; plat=${plat%_d46}
  if [ $k = SNV ]; then g=4 sl=16 vw=3; else g=1 sl=5 vw=6; fi
  NOVCF=1 KIND=$k GPUS=$g SLOTS=$sl MEM=192 VAL_WORKERS=$vw WAIT="$MIX && ended(pnv2_prebuild_${s}_$k)" \
      ris/snv_pipeline.sh ${samp}_${plat}_${k}_d46_7fc7396 $s
done
```

| run | GPU job | prebuild job | vcf job |
|---|---|---|---|
| `HG008_PacBio_SNV_d46_7fc7396` | 41543 | 41534 | 44040 |
| `HG008_PacBio_INDEL_d46_7fc7396` | 41544 | 41535 | 44035 |
| `HG008_ONT_SNV_d46_7fc7396` | 41545 | 41536 | 44041 |
| `HG008_ONT_INDEL_d46_7fc7396` | 41546 | 41537 | 44038 |
| `HG008_Illumina_INDEL_d46_7fc7396` | 41547 | 41538 | 44036 |
| `COLO829_Illumina_SNV_d46_7fc7396` | 41548 | 41539 | 44042 |
| `COLO829_Illumina_INDEL_d46_7fc7396` | 41549 | 41540 | 44037 |
| `COLO829_fiberseq_SNV_d46_7fc7396` | 41550 | 41541 | 44043 |
| `COLO829_fiberseq_INDEL_d46_7fc7396` | 41551 | 41542 | 44039 |
| `HG008_Illumina_SNV_d46_7fc7396` | 42232 (after its transfer, 2026-10-12 15:42) | 42231 | 44044 |
| `COLO829_ONT_SNV_d46_7fc7396` | 45059 (after its transfer, 2026-10-12 17:50) | 45052 | 45060 |
| `COLO829_ONT_INDEL_d46_7fc7396` | 45061 (after its transfer) | 45053 | 45062 |

**The d9 runs compared** (`collect.py` `d9_run`): HG008T the 7fc7396 runs of
`../hg008t_snv_indel_7fc7396_20261012` (the same recipes; SNV PacBio / ONT on one GPU, Illumina on 4); COLO829T SNV
`COLO829_<platform>_SNV_h90` (2026-10-03, the same SNV recipe on one GPU); COLO829T INDEL `COLO829_<platform>_INDEL_small`
(2026-10-06, lr 4e-4, re-scored with the 7fc7396 PoN). So the d46 - d9 difference is the graph alone only for HG008T
Illumina SNV and the HG008T INDELs; elsewhere it also holds the global batch (SNV: 4,096 against 1,024) or the
learning rate (COLO829T INDEL: 1e-4 against 4e-4).

## Results (chr1, all 12 runs)

All numbers: `tables.md`. In the GIAB / SMaHT BED; PASS = each model's validation recall-0.9 threshold; "best" = the
chr1 curve's best F1 (a test-set choice, for reference).

| sample | platform | kind | raw PASS F1, d9 -> d46 | raw best F1, d9 -> d46 | **PoN PASS F1, d9 -> d46** | PoN best F1, d9 -> d46 | PoN recall ceiling, d9 -> d46 |
|---|---|---|---:|---:|---:|---:|---:|
| HG008T | Illumina | SNV | 0.188 -> 0.295 | 0.285 -> 0.406 | **0.585 -> 0.570** | 0.607 -> 0.613 | 0.888 -> 0.892 |
| HG008T | PacBio | SNV | 0.191 -> 0.310 | 0.310 -> 0.414 | **0.621 -> 0.605** | 0.641 -> 0.622 | 0.888 -> 0.889 |
| HG008T | ONT | SNV | 0.182 -> 0.287 | 0.298 -> 0.404 | **0.596 -> 0.575** | 0.620 -> 0.604 | 0.897 -> 0.900 |
| HG008T | Illumina | INDEL | 0.101 -> 0.119 | 0.147 -> 0.173 | **0.241 -> 0.244** | 0.264 -> 0.249 | 0.257 -> 0.257 |
| HG008T | PacBio | INDEL | 0.041 -> 0.046 | 0.150 -> 0.159 | **0.089 -> 0.095** | 0.260 -> 0.269 | 0.259 -> 0.259 |
| HG008T | ONT | INDEL | 0.011 -> 0.011 | 0.090 -> 0.104 | **0.015 -> 0.015** | 0.140 -> 0.155 | 0.238 -> 0.237 |
| COLO829T | Illumina | SNV | 0.571 -> 0.685 | 0.628 -> 0.696 | **0.778 -> 0.783** | 0.780 -> 0.786 | 0.795 -> 0.793 |
| COLO829T | fiberseq | SNV | 0.555 -> 0.683 | 0.628 -> 0.700 | **0.781 -> 0.779** | 0.784 -> 0.783 | 0.792 -> 0.795 |
| COLO829T | ONT | SNV | 0.494 -> 0.626 | 0.586 -> 0.652 | **0.706 -> 0.737** | 0.717 -> 0.739 | 0.790 -> 0.792 |
| COLO829T | Illumina | INDEL | 0.031 -> 0.037 | 0.125 -> 0.105 | **0.115 -> 0.109** | 0.195 -> 0.163 | 0.275 -> 0.267 |
| COLO829T | fiberseq | INDEL | 0.010 -> 0.010 | 0.069 -> 0.081 | **0.016 -> 0.020** | 0.120 -> 0.129 | 0.244 -> 0.260 |
| COLO829T | ONT | INDEL | 0.002 -> 0.002 | 0.051 -> 0.056 | **0.003 -> 0.002** | 0.071 -> 0.085 | 0.237 -> 0.244 |

**After the PoN the two graphs are equivalent.** PASS F1 moves by -0.021 to +0.006 on 11 of 12 sets (COLO829T ONT SNV
+0.031). The only matched SNV pair (HG008T Illumina, both on 4 GPUs) gives -0.015 at PASS and +0.006 at the best
point; the HG008T PacBio / ONT drops (-0.016 / -0.021) compare 4-GPU d46 with 1-GPU d9 models, and the 1 -> 4 GPU change
alone cost HG008T Illumina d9 0.007. The PoN-limited recall ceilings are the same for both graphs (INDEL 0.257 / 0.259
/ 0.237 HG008T; SNV 0.89-0.90).

**Before the PoN, d46 is clearly better for SNV.** Raw PASS F1 rises by 0.10-0.13 on every SNV set (HG008T 0.18-0.19
-> 0.29-0.31, COLO829T 0.49-0.57 -> 0.63-0.69), the curves' best points by 0.07-0.12, the per-tensor somatic AP by
0.09-0.12 (HG008T 0.20-0.23 -> 0.31-0.35), with the same recall ceiling (HG008T 0.94-0.95). The d46 graph has more of
the population's variation, so fewer germline variants of the sample become candidates (HG008T PacBio chr1 SNV tensors
32,880 -> 27,501, COLO829T fiberseq 33,589 -> 26,608) and fewer reach the somatic call. Those are the variants the PoN
removes, so after it the advantage is gone: the graph and the PoN do the same job here.

**INDEL, by the chr1 curve's best F1 after the PoN** (`tables.md`, "INDEL: the chr1 curve's best F1"; the PASS threshold falls far below the best point for INDELs; the best point is chosen on chr1, an upper bound, not a held-out result): d9 -> d46 HG008T 0.264 -> 0.249 (-0.015), 0.260 -> 0.269 (+0.009), 0.140 -> 0.155 (+0.015) on Illumina / PacBio / ONT; COLO829T 0.195 -> 0.163 (-0.033), 0.120 -> 0.129 (+0.009), 0.071 -> 0.085 (+0.013) (COLO829T d9 at lr 4e-4; 6-25 true calls at its best points). The two graphs are equal within noise for INDEL too.

**INDEL: more truth INDELs get a tensor, which the PoN takes away again.** On HG008T the truth INDELs with a tensor
rise from 0.539 / 0.702 / 0.694 to 0.611 / 0.757 / 0.724 (predict's truth report, partial matches counted), raw best F1
by 0.009-0.026, but the PoN tags 413 of the 635 chr1 truth INDELs whatever the graph and the ceilings after it are
unchanged (0.257 / 0.259 / 0.237 vs 0.257 / 0.259 / 0.238): the INDELs d46 recovers are population alleles of gnomAD /
CoLoRSdb. (The COLO829T d9 INDEL per-tensor ceilings, 0.43 / 0.46 / 0.40, are from the 2026-10-06 predict that left
the partial matches out, so they are not comparable to d46's 0.57 / 0.70 / 0.65.)

**So**: with this pipeline's PoN, HPRC v2.1 d46 gives the same end result as v1.1 d9; its value shows only without a
PoN (raw SNV F1 +0.10-0.13), i.e. as a graph that does part of the PoN's work. Without the BED the PoN F1 differences
range from -0.026 (HG008T PacBio INDEL, 0.082 -> 0.056) to +0.031 (COLO829T ONT SNV) (`tables.md`).

