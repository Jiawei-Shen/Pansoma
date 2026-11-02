# E28: reruns of 2026-11-01 (masked retraining, INDEL checkpoint selection, shared validation split)

> **Running (submitted 2026-11-01 19:08 on tequila; `jobs.tsv`).** Results go to
> `/scratch/jshen/data/pansoma_net_v2_runs/{offnode_retrain_20261101,selcheck_20261101,sharedsplit_20261101}/`; this
> README gets the tables when the runs finish.

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
