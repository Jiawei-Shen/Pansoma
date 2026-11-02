# E28: reruns of 2026-11-01 (masked retraining, INDEL checkpoint selection, shared validation split)

> **Status (2026-11-02 12:06).** (1) done; (2) five of nine reruns done and scored, two running, two waiting for an H100;
> (3) waiting for an H100 (`jobs.tsv`). Results go to
> `/scratch/jshen/data/pansoma_net_v2_runs/{offnode_retrain_20261101,selcheck_20261101,sharedsplit_20261101}/`; the
> tables are made by `results/collect_*.py` (they write next to themselves).

The user asked to run the three experiments that the letter ver3 lists as not done. `submit.sh` submits all of them:

1. **R1.1, cross-node context, retrained** (`snv/snv_train.sh` MODE `mask` / `nomask`). The six HPRC v2.1 SNV sets,
   seed 0, the recipe of Table R3.1 (50.4 M, 20 epochs, pow:0.75, lr 1e-4, `--norm train`, global batch 4096) on one
   H100 (4 micro-batches of 1024; the RIS controls ran on 4 GPUs x 1024), trained with `--mask-off-node` (every cell
   outside the candidate's own node visit zeroed in training, validation and prediction; E27 `ablation/q6/offnode.py`)
   and without it on the same setup as the control. Code: the E27 q6 copy (repository `pansoma_net_v2` of 90d7c0f +
   `offnode_vs_repo.patch`).
2. **R1.8, checkpoint selection without the off-path validation negatives** (`indel_select/`). The nine final INDEL
   release models (`release_s{0,1,2}_{Illumina,PacBio,ONT}_INDEL`, 6 tumors, chr1 held out) retrained with the same
   arguments, seed, micro-batches and GPU type as the originals (H100 for micro-batches 1, a 24gb MIG slice for 4),
   by the frozen `tmp/indel_gap_20261020/code_final` (77cee0e) plus `train_select.patch`: each epoch also logs the
   selection score (post-filter validation F1) with the off-path candidates without a truth match left out
   (`val.excl_offpath` in `metrics.jsonl`) and keeps that selection's checkpoint (`best_excl_offpath.pth`). Training
   itself is unchanged, so a rerun that reproduces the original per-epoch metrics shows which epoch the other
   selection would have kept for the reported model.
3. **R3.1, a validation split shared by the three references** (`snv/snv_train.sh` MODE `shared`). HG008T and
   COLO829T, three platforms, GRCh38-only, HPRC v1.1 and HPRC v2.1, seed 0, the same recipe and setup as (1) without
   the mask, with `--val-chroms chr20 chr21 chr22` instead of the hash-selected node blocks (which differ between
   graphs), so every reference selects its checkpoint and threshold on the same GRCh38 region; training chr2-19, test
   chr1.

Every SNV training has a CPU job for `jobs/vcf.sh` (chr1, PASS, population filter, rtg vcfeval).

## Results

### (1) Cross-node context, retrained (`results/offnode.md`, `collect_offnode.py`)

chr1 SNVs, benchmark BED, after the PoN, PASS at the checkpoint's validation recall-0.9 threshold; best = best F1 of the
chr1 curve (test-tuned); AP = chr1 per-tensor somatic AP before the PoN; epoch = validation-selected epoch.

| chr1 SNV, v2.1, seed 0 | no mask P / R / F1 (best; AP; epoch) | mask P / R / F1 (best; AP; epoch) | ΔF1 PASS | Δbest | ΔAP | seed range (Table R3.1, 3 seeds, no mask) |
|---|---|---|---|---|---|---|
| HG008_Illumina | 0.462 / 0.852 / 0.599 (0.619; 0.313; 9) | 0.446 / 0.858 / 0.587 (0.610; 0.323; 6) | -0.012 | -0.009 | +0.010 | 0.024 |
| HG008_PacBio | 0.475 / 0.835 / 0.606 (0.621; 0.383; 14) | 0.467 / 0.831 / 0.598 (0.606; 0.303; 9) | -0.008 | -0.015 | -0.080 | 0.008 |
| HG008_ONT | 0.451 / 0.862 / 0.592 (0.612; 0.339; 14) | 0.457 / 0.864 / 0.598 (0.608; 0.329; 17) | +0.005 | -0.004 | -0.010 | 0.023 |
| COLO829_Illumina | 0.866 / 0.714 / 0.783 (0.785; 0.810; 8) | 0.849 / 0.710 / 0.773 (0.778; 0.807; 8) | -0.009 | -0.006 | -0.003 | 0.079 |
| COLO829_fiberseq | 0.872 / 0.701 / 0.777 (0.793; 0.823; 19) | 0.840 / 0.718 / 0.774 (0.782; 0.790; 17) | -0.003 | -0.012 | -0.032 | 0.015 |
| COLO829_ONT | 0.816 / 0.700 / 0.753 (0.754; 0.760; 8) | 0.781 / 0.698 / 0.737 (0.740; 0.740; 8) | -0.016 | -0.014 | -0.020 | 0.018 |

Trained with the mask, PASS F1 is lower in 5 of 6 (-0.003 to -0.016), the best F1 in 6 of 6 and AP in 5 of 6; each
PASS F1 difference is smaller than the three-seed range in 5 of 6 and the same size in HG008 PacBio. At inference only
(E27 `ablation/q6/infer`), the same mask cost up to 0.074 PASS F1: models trained with the context use it, but a model
trained without it recovers nearly all of the accuracy. The one-GPU controls are within 0.03 F1 of the four-GPU seed-0
models of Table R3.1 (HG008 Illumina +0.028, HG008 ONT +0.018, COLO829 ONT +0.017; the others within 0.003).

### (2) INDEL checkpoint selection without the off-path validation negatives (`collect_select.py`, `collect_select_eval.py`)

`results/select.md` (all reruns: reproduction of the original per-epoch metrics, epoch kept by each selection) and
`results/select_eval.md` (test F1 of both checkpoints where they differ; `indel_select/predict_ckpt.sh`,
`submit_eval.sh`, jobs in `results/eval_jobs.tsv`). Five of nine done (Illumina seeds 0-2, PacBio seeds 0-1):

- Reproduction: only Illumina seed 0 is bit-identical to the original; the others differ from epoch 1 on (train loss
  up to 4e-4) on the same GPU type (H100 NVL 24gb MIG slice), so comparisons are within one rerun. With the original
  rule, 4 of 5 keep the reported epoch (PacBio seed 0: 6 instead of 7), and their test results are within 0.011 F1 of
  the reported ones.
- Without the off-path negatives, selection keeps a later epoch in 4 of 5 (2->3, 1->6, 4->6, 6->8; PacBio seed 0
  same). Test F1 after the PoN (ppr8): higher in 10 of 12 chr1 sets (up to +0.142) and in 4 of 4 HCC1395T chr20
  (+0.046 to +0.171), lower in two COLO829T sets (-0.009, -0.031); HG008T recall higher in all four. The
  off-path validation negatives make the selection favour earlier, more conservative checkpoints.
- Letter: ver3 R1.8 and Supplementary Table R1.8 (with a [待补] for the four ONT / PacBio seed-2 reruns).
