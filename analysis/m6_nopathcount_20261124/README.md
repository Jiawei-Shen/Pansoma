# Six-tumor SNV models on HPRC v2.1 without the haplotype-count channel (2026-11-24)

User 2026-11-24: the with / without haplotype-count comparison (Response Figure 1.3; in ver3 on older models, in the
cluster's ver4 on the COLO829T single-tumor models) on the six-tumor HPRC v2.1 models of the latest pipeline, tested
on the three platforms.

- With the channel, existing: seed 0 `snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe` (old cluster; the six-tumor
  model of Supplementary Fig. R1.6 / Table R1.7b), seeds 1 and 2 `m6_seeds_20261123/multicancer6_<P>_SNV_newpipe_s{1,2}`
  (RIS).
- Without the channel, this batch: `$PANSOMA_RUNS/m6_nopathcount_20261124/multicancer6_<P>_SNV_newpipe_nopc`, P =
  Illumina / PacBio (COLO829T: Fiber-seq) / ONT. Seed 0 with only `--drop-planes path_count` added (the encoder plane
  made from channel 6 is not an input, and the channel is never read); everything else as `../m6_seeds_20261123`: the
  `ris/recipes.sh` SNV recipe with the site features as a model input, the six d46 sets (Illumina / PacBio HiFi on the
  AF 0.05 builds of the four CASTLE lines and COLO829T, HG008T and ONT at 0.08), HCC1395T chr1 + chr20, HG008T chr1
  and COLO829T chr1 left out of training and validation, 20 epochs, batch 1024, one H100 (not DDP). pansoma_net_v2
  code as for the seeds (last changed in 1798490).
- Test: chr1 of the three test tumors (HCC1395T also chr20) predicted with `best.pth`, then `ris/vcf.sh --set` per
  tumor (`vcf_chr1_<SET>/`; PASS = rule 2; population filter; rtg vcfeval in the benchmark BED).
- One model without the channel per platform (seed 0) against three seeds with it: a difference inside the seed
  spread of the models with the channel (`results.md`, last table) is not evidence of an effect of the channel.

Files: `submit.sh` (login node; uses `../m6_nositef_graphs_20261118/train_predict.sh` and `ready.py`), `collect.py` ->
`results.md` (per model, before and after the filter; with − without against seed 0 and against the three seeds).

## Status

| model | submitted | LSF jobs (train; VCF HCC1395T, HG008T, COLO829T) |
|---|---|---|
| Illumina | 2026-11-24 14:29 | 4359; 4360, 4361, 4362 |
| PacBio HiFi | 2026-11-24 14:29 | 4363; 4364, 4365, 4366 |
| ONT | 2026-11-24 14:29 | 4367; 4368, 4369, 4370 |

At submission all eight H100s of exec-419 were in use (one by the Illumina seed 1 training of `m6_seeds_20261123`),
so the three trainings wait on the epigenome queue.

All three trained and scored by 2026-11-26 06:19 (best epochs 10 / 11 / 8 for Illumina / PacBio HiFi / ONT); `results.md`
regenerated 2026-11-26. chr1 F1 after the filter, mean of the three seeds with the channel minus the model without:
HG008T +0.047 / +0.007 / +0.045, COLO829T +0.074 / +0.017 / +0.005, HCC1395T +0.007 / +0.002 / +0.035 (Illumina /
PacBio HiFi or Fiber-seq / ONT); seed 0 with the channel is above the model without in all nine cells. The largest
gap, COLO829T Illumina, is precision: 871 false calls without the channel against 335-358 for seeds 1 and 2.
