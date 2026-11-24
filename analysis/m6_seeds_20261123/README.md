# Six-tumor SNV models of the current pipeline on HPRC v2.1: seeds 1 and 2 (2026-11-23)

User 2026-11-23: a set of multi-seed six-tumor v2.1 models of the latest pipeline, tested on the three datasets
(HCC1395T, COLO829T, HG008T).

- Seed 0 exists: `snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe` (old cluster, one H100), the six-tumor model of
  Supplementary Fig. R1.6 and Table R1.7b. This batch trains seeds 1 and 2 with only `--seed` changed, so the three
  seeds give the spread of that model.
- Recipe: `ris/recipes.sh` SNV (repository 8e6fb71; pansoma_net_v2 code last changed in 1798490, which only adds
  `--site-features-pon-only`, not used here): the SNV site features as a model input, non tensors with AF <= 0.1
  sampled 25 % per epoch at loss weight 1, best epoch by the post-PoN validation F1, PASS at the post-PoN validation
  recall 0.95 of the callable positives (rule 2), 20 epochs, batch 1024, one H100 (not DDP, as seed 0).
- Sets: the six d46 sets of seed 0: Illumina / PacBio HiFi on the AF 0.05 builds of the four CASTLE lines and COLO829T
  (`_d46_snvaf05`), HG008T and ONT at AF 0.08. HCC1395T chr1 + chr20, HG008T chr1 and COLO829T chr1 left out of
  training and validation; HCC1954T, H1437T and H2009T train on every chromosome.
- Test: chr1 of the three test tumors (HCC1395T also chr20) predicted with `best.pth`, then `ris/vcf.sh --set` per
  tumor (`vcf_chr1_<SET>/`: graph VCF, GRCh38 VCF, population filter, rtg vcfeval in the benchmark BED and without it).
- Runs: `$PANSOMA_RUNS/m6_seeds_20261123/multicancer6_<P>_SNV_newpipe_s{1,2}`.
- INDEL: not part of this batch. The six-tumor INDEL model of the final recipe already has three seeds
  (`release_s{0,1,2}_<P>_INDEL`, old cluster), and the CASTLE lines' INDEL tensors are not on storage2.

Files: `submit.sh` (login node; uses `../m6_nositef_graphs_20261118/train_predict.sh` and `ready.py`), `collect.py` ->
`results.md` (seeds 1 and 2 from their reports, seed 0 from `r17b_latest.md`; mean and SD of the PASS F1).

## Status

| model | submitted | LSF jobs (train; VCF HCC1395T, HG008T, COLO829T) |
|---|---|---|
| Illumina seed 1 | 2026-11-23 11:24 | 989814; 989815, 989816, 989817 |
| Illumina seed 2 | 2026-11-23 11:24 | 989818; 989819, 989820, 989821 |
| PacBio HiFi seed 1 | 2026-11-23 11:24 | 989822; 989823, 989824, 989825 |
| PacBio HiFi seed 2 | 2026-11-23 11:24 | 989826; 989827, 989828, 989829 |
| ONT seed 1 | 2026-11-23 11:24 | 989830; 989831, 989832, 989833 |
| ONT seed 2 | 2026-11-23 11:24 | 989834; 989835, 989836, 989837 |

Memory (2026-11-23): LSF's job memory includes the page cache of the chunk files the job reads, which grows by
1-3 GB/min until the job's limit and is then reclaimed; the trainings of 2026-11-18 reached their limit that way,
not by need. Measured with `bjobs -l` during epoch 1: PacBio HiFi 33-41 GB in training, 63 GB right after the first
validation; Illumina 55-75 GB and growing with the cache. The four Illumina / PacBio jobs keep the requests they
started with (384 / 256 GB). The two ONT jobs, pending on memory, were lowered to 128 GB with bmod (job group now
/jiawei.s/default); ONT seed 1 started at once, seed 2 when another user's job freed an H100 (~12:30; the user chose to wait on
epigenome rather than use a general-queue A100). All six run on exec-419, one H100 each. LSF refuses to
change a running job ("Cannot change job options for a running job"), and killing the four to resubmit would hand
their GPUs to the other user's queued jobs. submit.sh now asks for 192 GB (Illumina) / 128 GB (PacBio HiFi, ONT).
