# HG008T SNV and INDEL, three platforms: training and chr1 test with the pipeline of 7fc7396

> 2026-10-16: the INDEL numbers here use the 7fc7396 INDEL PoN rule (homopolymer >= 7 bp: AF >= 0.01, else 1e-3). With
> the 71cb51a rule (>= 8 bp, else 1e-4) every INDEL call set was re-scored: `../indel_pon_71cb51a_20261016/tables.md`
> (+0.000 to +0.013 F1, no comparison changes direction).

2026-10-12, on RIS compute1. The SNV recipe and the INDEL recipe of `machine_learning/pansoma_net_v2/README.md` at
7fc7396 (INDEL PoN homopolymer rule, INDEL recipe lr 1e-4, INDEL tests keep the partial somatic matches), one model
per platform and kind on HG008T Illumina, PacBio and ONT, trained on chr2-22 and tested on chr1.

Files in this folder:

- `collect.py`: reads the run directories (`$PANSOMA_RUNS/<run>`) and writes `tables.md`, `summary.json` and
  `runs/<run>/` (copies of `train.log` and the rtg `report.txt` files). Standard library only:
  `PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py`.
- `tables.md`: training, chr1 per tensor / per truth allele, rtg vcfeval in the BED and without it, comparison with
  the earlier runs.
- `summary.json`: the numbers behind the tables.
- `runs/`: per-run copies of the logs and reports, for the six runs of this analysis and the earlier runs they are
  compared with.

## What was run

Data: the RIS copies of the HG008T sets (`Dataset/pansoma_v2_tensors/HG008T_{Illumina,PacBio,ONT}/tensors`, build
9d6bd8d: `--haplotypes 90`, ch6 `int8-count-haplotypes100`, build AF 0.08; labels `5430533aa100` of 2026-10-02, no AF
floors, `ignore_repeat_indels` off). The sets are unchanged since the 2026-10-03 / 10-06 runs.

Code: the repository at 7fc7396 plus `ris/`. CPU unit tests of 7fc7396 on RIS (`pnv2_setup_7fc7396`): 79 tests OK,
8 skipped (GPU, real data, rtg).

Per run, two LSF jobs (`ris/snv_pipeline.sh`; `KIND=INDEL` for INDEL):

1. `run.sh` on one H100 80GB (`epigenome` queue, compute1-exec-419): `pansoma_net_v2.train` with the recipe, chr1
   left out, 5 % of the chr2-22 node blocks validate; then `pansoma_net_v2.predict` on every chr1 tensor with
   `best.pth`.
2. `vcf.sh` on CPU: `graph_vcf` (PASS = the checkpoint's validation recall-0.9 threshold; SNV off-reference rescue
   on, the default) → `linear_vcf` → `vcfeval` with the four PoNs (the PoN script's rule of each record's kind) and
   rtg vcfeval, in the GIAB BED (`eval_bed`) and on the whole chromosome (`eval_nobed`).

Recipes (`ris/snv_pipeline.sh`, as in the v2 README at 7fc7396; `--read-threads 8 --val-workers 6` added):

| kind | model | options |
|---|---|---|
| SNV | 50.4 M (`--depths 3 3 27 3 --dims 96 192 384 768`) | `--epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75 --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth` |
| INDEL | 8.7 M (`--depths 2 2 6 2 --dims 64 128 256 512`) | `--epochs 12 --batch-size 1024 --lr 1e-4 --class-weights sqrt --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth` (no `--exclude-candidates`) |

The INDEL recipe in `ris/snv_pipeline.sh` was changed from lr 4e-4 (45d56ad) to 1e-4 for these runs.

PoN rules (`scripts/filter_panel_of_normals.py` at 7fc7396): SNV gnomAD / CoLoRSdb exact allele at AF ≥ 0.0001,
dbSNP records not flagged somatic, 1000G; INDEL gnomAD / CoLoRSdb exact allele at AF ≥ 0.01 when the ALT changes the
length of a homopolymer of ≥ 7 bases on GRCh38, AF ≥ 0.001 otherwise.

| run | set | GPU job | vcf job | resources (GPU job) |
|---|---|---|---|---|
| `HG008_Illumina_SNV_7fc7396_ddp4_b1024` | HG008T_Illumina | 40970 | 40971 | 4 H100s (DDP), 16 slots, 256 GB |
| `HG008_PacBio_SNV_7fc7396` | HG008T_PacBio | 40873 | 40874 | 7 slots, 192 GB |
| `HG008_ONT_SNV_7fc7396` | HG008T_ONT | 40877 | 40878 | 7 slots, 192 GB |
| `HG008_Illumina_INDEL_7fc7396` | HG008T_Illumina | 40871 | 40872 | 5 slots, 192 GB |
| `HG008_PacBio_INDEL_7fc7396` | HG008T_PacBio | 40875 | 40876 | 5 slots, 192 GB |
| `HG008_ONT_INDEL_7fc7396` | HG008T_ONT | 40879 | 40880 | 5 slots, 192 GB |

Submitted 2026-10-12 12:22 CDT from the login node:

```bash
cd /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma
ris/submit.sh -g 0 -n 4 -m 16 pnv2_setup_7fc7396 setup_env.sh
for p in Illumina PacBio ONT; do
  MEM=192 ris/snv_pipeline.sh HG008_${p}_SNV_7fc7396 HG008T_$p
  KIND=INDEL SLOTS=5 MEM=192 ris/snv_pipeline.sh HG008_${p}_INDEL_7fc7396 HG008T_$p
done
```

**HG008T Illumina SNV on 4 GPUs.** It was first submitted on one H100 like the others (`HG008_Illumina_SNV_7fc7396`,
40869): 1.9 M training tensors at ~2,300 tensors/s in epoch 1, about 3-4 h for 20 epochs with five other runs on the
host. It was stopped in epoch 1 and resubmitted on 4 H100s under DDP:
`GPUS=4 SLOTS=16 MEM=256 VAL_WORKERS=3 ris/snv_pipeline.sh HG008_Illumina_SNV_7fc7396_ddp4_b1024 HG008T_Illumina`.
`--batch-size 1024` is per GPU under DDP (the chunk loader splits each epoch over the ranks), so this run has a global
batch of 4,096, a quarter of the steps per epoch (465 instead of 1,860) and the same lr 1e-4: it is the SNV recipe
except for the global batch. A first 4-GPU attempt with `--batch-size 256` per GPU (`HG008_Illumina_SNV_7fc7396_ddp4`,
40955) was stopped during torch.compile when the per-GPU batch was set back to 1,024. Neither partial directory is used.

**Earlier runs, for comparison.** SNV: `HG008_<P>_SNV_h90` (2026-10-03, repo 9d6bd8d, the same SNV recipe, same
sets; the SNV PoN rule has not changed since). INDEL: `HG008_<P>_INDEL_small` (2026-10-06, repo 45d56ad, the INDEL
recipe at lr 4e-4). Their chr1 calls were re-scored with the PoN script of 7fc7396
(`ris/vcfeval_only.sh RUN INDEL vcf_chr1 vcf_chr1_pon7fc7396`, jobs 40898-40900), so the INDEL comparison isolates the
learning rate; `vcf_chr1` keeps their scores under the earlier INDEL PoN rule (AF ≥ 0.01 for every INDEL).

## Results (chr1)

All tables, both BED settings and the training columns: `tables.md`. PASS = the checkpoint's validation recall-0.9
threshold; "best" = the best point of the chr1 rtg curve (chosen on the test chromosome: for reference only, never a
threshold). In the GIAB benchmark BED after the PoN of each kind:

| run | PASS threshold | PASS evaluated (unplaced) | PoN TP / FP | PoN P / R / **F1** | chr1 best F1 (at p) | no BED, PoN F1 | earlier run, same PoN: PASS F1 (best) |
|---|---:|---:|---:|---:|---:|---:|---:|
| SNV Illumina (4 GPUs) | 0.251 | 5,682 (116) | 579 / 703 | 0.452 / 0.831 / **0.585** | 0.607 (0.388) | 0.577 | 0.592 (0.612) |
| SNV PacBio | 0.314 | 5,540 (211) | 577 / 584 | 0.497 / 0.828 / **0.621** | 0.641 (0.395) | 0.614 | 0.621 (0.641) |
| SNV ONT | 0.449 | 5,843 (1,409) | 576 / 661 | 0.466 / 0.826 / **0.596** | 0.620 (0.612) | 0.586 | 0.596 (0.620) |
| INDEL Illumina | 0.119 | 3,596 (563) | 140 / 434 | 0.244 / 0.238 / **0.241** | 0.264 (0.209) | 0.233 | 0.246 (0.263) |
| INDEL PacBio | 0.040 | 11,759 (2,075) | 150 / 2,630 | 0.054 / 0.256 / **0.089** | 0.260 (0.259) | 0.082 | 0.097 (0.231) |
| INDEL ONT | 0.018 | 41,052 (7,288) | 134 / 17,268 | 0.008 / 0.228 / **0.015** | 0.140 (0.168) | 0.015 | 0.016 (0.147) |

**INDELs are compared by the chr1 best F1 column** (0.264 / 0.260 / 0.140 for Illumina / PacBio / ONT): the PASS
threshold, the validation recall-0.9 point, lands far below the best point on the long reads (PacBio 0.089 vs 0.260,
ONT 0.015 vs 0.140). The best point is chosen on chr1: an upper bound for a threshold rule, not a held-out result.

Truth on chr1: 697 SNP / 587 INDEL records in the BED (702 / 635 without it). "PASS evaluated" is before the PoN;
unplaced = PASS calls without a GRCh38 position (not scored by rtg; for SNV mostly the off-reference rescue).

**SNV (PacBio, ONT).** The runs reproduce the 2026-10-03 runs (`HG008_<P>_SNV_h90`) exactly: the same best epoch (7,
8), the same PASS thresholds (0.3141991198 / 0.4492301941), the same rtg TP / FP raw and after the PoN and the same 211 /
1,409 unplaced PASS calls. Between 9d6bd8d and 7fc7396 the train / predict / graph_vcf / data code changed only in
7fc7396 (`--exclude-candidates`, the INDEL partial default), which an SNV run without an exclusion file does not
touch, the SNV PoN rule is unchanged, and training is deterministic for a given seed and set. PoN F1 0.621 / 0.596, best 0.641 / 0.620; the PoN tags 61 of the 702 chr1 truth SNVs (recall
bound 0.913).

**SNV Illumina on 4 GPUs.** Global batch 4,096 (`--batch-size 1024` per GPU), so a quarter of the 1-GPU steps per
epoch at the same lr: best epoch 4, PoN F1 0.585 (best 0.607) against 0.592 (0.612) for the 1-GPU run of 2026-10-03,
the same recipe otherwise. Training took 125 min (20 epochs, ~6 min each, 5,280 tensors/s on 4 H100s, other runs sharing
the host throughout) against ~3 h on one GPU.

**INDEL, lr 1e-4 (the 7fc7396 recipe) against lr 4e-4 (45d56ad), both scored with the 7fc7396 PoN.** On RIS the lower
learning rate does not change the PASS F1: Illumina 0.241 vs 0.246, PacBio 0.089 vs 0.097, ONT 0.015 vs 0.016; the
best points are 0.264 vs 0.263, 0.260 vs 0.231, 0.140 vs 0.147. Upstream's figures for the unfiltered recipe
(v2 README at 7fc7396) are 0.239 / 0.132 / 0.015: Illumina and ONT agree, PacBio's PASS F1 is lower here while its
best point agrees (0.260 vs 0.254), so the difference is the threshold, not the ranking.

**INDEL PoN rule (7fc7396).** On the same lr 4e-4 calls, the homopolymer rule (AF ≥ 0.01 in homopolymers ≥ 7 bp,
≥ 0.001 elsewhere) against the earlier rule (AF ≥ 0.01 for every INDEL): Illumina 0.218 → 0.246, PacBio 0.089 →
0.097, ONT 0.016 → 0.016 (best 0.243 → 0.263, 0.220 → 0.231, 0.138 → 0.147). The new rule tags 413 of the 635 chr1
truth INDELs (the old one 383), so after it no caller recalls more than 0.350 of them; the rtg curve ceilings after
the PoN are 0.24–0.26 (0.40–0.47 before).

**What limits INDEL on the long reads is the decision rule.** The validation recall-0.9 threshold falls to 0.040
(PacBio) and 0.018 (ONT), so 11,759 / 41,052 calls PASS (2,780 / 17,402 after the PoN) at precision 0.054 / 0.008,
while the chr1 curves peak at F1 0.260 / 0.140 near p ≈ 0.26 / 0.17. On Illumina the PASS point and the best point are
close (0.241 / 0.264). The chr2-chosen ensemble + rescorer of 2026-10-07 (`ris/README.md`, "INDEL: a general method",
under the earlier PoN rule) reached 0.275 / 0.274 / 0.180 on the same sets.

**Per-tensor INDEL numbers are not comparable with the earlier runs.** From 7fc7396, predict keeps the partial
somatic matches in the INDEL test as 1 (somatic test tensors: Illumina 364, earlier 236; truth ceiling 0.539, earlier
0.400). The rtg results do not depend on it.

**Run times** (H100 80GB, training only, six runs on one host): SNV Illumina 125 min on 4 GPUs, PacBio 36 min, ONT 48 min (20 epochs, ~3,610
tensors/s); INDEL Illumina 3 min, PacBio 33 min, ONT 44 min (12 epochs, ~8,800 tensors/s). Then predict chr1 and
`vcf.sh` (a few minutes each).
