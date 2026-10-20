# E21: what the models are trained on (training-set composition)

Reviewer points: **R1.7** (training scale and tumor diversity: exact parameter count, positive and negative examples
per platform and variant type, whether the data volume supports the capacity, comparison with DeepSomatic's training)
and **R1.8** (how the training set is built: positives, negatives, what "determinable" means, germline, artefacts, VAF
and purity). E21 also prepares **E13** (learning curve and model size), which needs GPUs.

## Summary

Status **PENDING** (mostly done). Done: the per-run table for all 42 models (section 1), the label rules (section 4),
the comparator scale (section 5), the E13 input files and job commands (section 6), and sections 2-3b for the 21 sets
whose index was built: the 18 reference-arm sets (HG008T / COLO829T x 3 platforms x d9 / d46 / linear) and the three
COLO829T Illumina dilutions (P75 / P50 / P25). The replay matches `train.log` exactly for 38 of the 42 runs (0 differ).
Pending: the index builds of the other 12 sets (COLO829T fiberseq / ONT dilutions, HCC1395T x 3, the CASTLE PacBio
sets; Slurm array 387137 tasks 22-33) and the final collect job 387193 are **held** (`JobHeldUser`) after the user's
complaint about the lab-cluster footprint and stay held until the user agrees a cap; E13 itself needs RIS GPUs and the
exclude files copied to RIS by the user. `tables.md` and the TSVs were regenerated on 2026-10-13 23:01 from the 21
finished caches (`collect_runs.py`, `summarize.py`, `make_tables.py` on the login node; no new jobs).

- **Parameters** (from `train.log`): SNV model 50.4 M (ConvNeXt depths 3-3-27-3, widths 96-768), INDEL model 8.7 M
  (2-2-6-2, 64-512).
- **Positives per model** (chr2-22 minus the validation blocks; partial matches left out, so every positive is an exact
  allele match): HG008T SNV 7,098-7,391 somatic tensors, HG008T INDEL 3,121-6,224; COLO829T SNV 33,472-34,341,
  COLO829T INDEL 971-1,641. In every single-set SNV model each positive is a different truth allele (e.g. HG008T PacBio
  d9: 7,098 tensors, 7,098 alleles). Over the three platforms the distinct SNV alleles behind training tensors are
  7,203 / 7,337 / 7,448 for HG008T (d9 / d46 / linear; 7,935 over all three arms) and 34,508 / 34,853 / 34,971 for
  COLO829T (36,880 over all arms). The COLO829T joint models see 124,070-125,916 SNV positives, but they are four
  purities of one truth set: the Illumina joint SNV model's 125,916 positives are 33,797 distinct alleles (P100 alone:
  33,782), its 3,439 INDEL positives 973 alleles (P100: 967). The fiberseq / ONT joint unions are pending (held jobs).
- **Negatives per model**: germline 256,583-381,159 SNV tensors on the two pangenomes and 2.71-3.01 M on GRCh38-only
  (every germline variant of the sample becomes a candidate there); non 67,735-1,832,960 SNV and 20,551-3,289,717 INDEL
  tensors in the single-set models, depending mostly on the platform (Illumina SNV and long-read INDEL errors).
- **Capacity**: about 7,000 parameters per HG008T SNV positive. The SNV networks can memorise the training set: the
  training loss at epoch 20 is 0.002-0.139 while the validation loss reaches its minimum at epoch 2-11 and ends 1.4-9.5x
  higher; the checkpoints (best validation somatic F1) are from epochs 4-12 of 20 (INDEL: 2-7 of 12; validation loss
  ends 1.0-2.5x its minimum). Early stopping on the held-out validation blocks is what keeps the model general; the
  learning curve (E13) is needed to show whether more data would still help. Inputs for it are ready (section 6).
- **Which truth reaches training** (share of the PASS somatic truth in the training region that stands behind a training
  tensor): SNV 94.8-95.9 % (HG008T) and 89.5-90.5 % (COLO829T) on the pangenome arms, 98.5-98.7 % and 90.1-91.4 % on
  GRCh38-only. The COLO829T loss is mostly the build filters (AF 0.08, 3 ALT reads: 2,542-3,430 alleles, its low-VAF
  subclonal SNVs; HG008T 45-71). INDEL: only 42.6-50.5 % (HG008T) and 57.3-66.7 % (COLO829T) on the pangenome arms, where
  1,991-3,209 HG008T INDELs have no candidate at all; on GRCh38-only 57.4-83.7 % and 75.4-96.7 %, with 116-272 (HG008T)
  and 4-49 (COLO829T) without a candidate. Most of the missing INDEL positives are lost on the pangenome arms, before
  labelling.
- **VAF**: the somatic SNV training tensors are mostly clonal: median AF 0.53-0.54 (HG008T) and 0.46-0.51 (COLO829T),
  all arms; AF < 0.20 for 0.4-1.8 % (HG008T) and 6.5-7.4 % (COLO829T); 31-33 % of the HG008T ones have AF >= 0.80. On the
  pangenome arms the germline SNV tensors have nearly the same AF distribution (median 0.48-0.60), so AF alone does not
  separate the two classes; on GRCh38-only the germline median is 0.76-0.99. Non tensors sit at low AF (SNV median
  0.10-0.24). The COLO829T Illumina dilutions move the somatic SNV median to 0.40 / 0.28 / 0.15 (P75 / P50 / P25;
  AF < 0.20: 8.8 / 25.0 / 82.9 %): only the joint models see low-VAF positives. HCC1395T: pending (held jobs).
- **Tumor diversity of what is learned**: the COLO829T positives are 71 % C>T, 65 % C>T at a dipyrimidine (UV), 0.8 %
  C>G at TCW; the HG008T positives 40-41 % C>T, 28-29 % at a dipyrimidine, 19 % C>A, 7.9-8.5 % C>G at TCW
  (APOBEC-like), 9.2-9.6 % C>T at CpG (all arms). The two training lines have clearly different spectra; HCC1395T and the
  CASTLE lines: pending (held jobs).
- **Replay check**: the training / validation counts re-derived from the label files equal the `train.log` counts
  (per class, partial left-out, validation off-reference) of 38 of the 42 runs, 0 differ; the 4 unchecked runs are the
  fiberseq and ONT joint models, whose P75 / P50 / P25 caches are in the held array tasks.

## Question

For every model of the 2026-10-12 benchmark (36 single-set runs: HG008T and COLO829T x Illumina / PacBio (COLO829T:
fiberseq) / ONT x SNV / INDEL x three reference arms; plus the 6 COLO829T purity joint runs):

1. How many training and validation tensors per class (non / somatic / germline), how many ignored, epochs, best epoch,
   validation somatic AP and PASS threshold?
2. How many distinct somatic truth alleles stand behind the training tensors?
3. What is the VAF (tensor AF) distribution of the somatic and germline training tensors, per set, and how do HG008T,
   COLO829T and HCC1395T compare?
4. What are the label rules, in plain English, with the code lines?
5. How does this compare with the comparators' training scale?
6. Which GPU jobs give the learning curve (E13)?

## Data

- Models: `/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012/<run>/` (`train.log`, `metrics.jsonl`,
  `args.json`; read-only) for 36 runs. The six COLO829T HPRC v1.1 d9 single-set models
  (`COLO829_<platform>_SNV_h90`, 2026-10-03; `COLO829_<platform>_INDEL_small`, 2026-10-06, lr 4e-4) are not on this cluster:
  their `train.log` copies in `analysis/benchmark_7fc7396_20261012/*/runs/<run>/` and the training block of that
  folder's `summary.json` (validation AP and recall-0.9 threshold, collected on RIS from `metrics.jsonl`).
- Tensor sets (read-only): `/scratch/jshen/data/pansoma_v2_tensors/<SET>/tensors/{SNV,INDEL}/` (`<chrom>_labels.ndjson`,
  `<chrom>_variant_summary.ndjson`, `*_labels.npy`, `labels.manifest.json`), `tensors/somatic.recall.tsv`,
  `truth/somatic.graph.tsv`. 33 sets: the 18 reference-arm sets (HG008T / COLO829T x 3 platforms x d9 / d46 / linear),
  the 9 COLO829T dilutions (P75 / P50 / P25, d9), HCC1395T x 3 platforms (d46) and, as a snapshot, the three CASTLE
  PacBio sets whose labels were final on 2026-10-13 (HCC1954T, H1437T, H2009T; labels of 2026-10-12; their ONT /
  Illumina sets were being relabelled or built and were not read). **Read so far: 21 of the 33** (the 18 reference-arm
  sets and `COLO829T_Illumina_P75` / `_P50` / `_P25`); the index caches of the other 12 (array tasks 22-33, in the
  order of `sets.txt`) are not built because those tasks are held. Before the summarize run every one of the 21 caches
  was checked to be current (index version, labels meta and directory equal to the set's labels today), so the run
  read them and built nothing.
- Every set read here has the labels the models were trained on: `labels.manifest.json` `rules_sha256` `5430533aa100`
  and the same `created` time stamp as in each `train.log` (checked for all 60 set x kind pairs of the 42 runs:
  `runs_training.labels_created.json` against the local manifests, 0 differ), and the replayed counts equal the logged
  ones (below). HCC1395T and CASTLE: `17b540137dc3` (the same file plus the
  `ignore_repeat_indels` option; `git diff 1e9ff6d 4f36af1 -- indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py`
  changes only `RepeatContext`, the `non()` wrapper and the option's plumbing).

## Method

1. `collect_runs.py` (standard library): parses every `train.log` (class counts, `--ignore-reasons` counts, parameter
   count, GPUs, per-epoch train / validation loss and somatic AP, best epoch) and the best epoch's record of
   `metrics.jsonl` (validation somatic AP; `somatic_at_recall["0.9"].threshold` = the PASS threshold of graph_vcf).
   Writes `runs_training.tsv` / `runs_training.md`.
2. `build_index.py` (Slurm array `build_index.sbatch`, job 387137, 33 tasks, `-c 2 --mem 24G`, `-p general`): builds
   `pansoma_net_v2.data.KindIndex` (imported, unchanged) of each set into this experiment's own cache
   (`/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/index_cache`). The index holds, per
   tensor, chromosome, node, label, evaluation label, label reason, A1 AF (`scalars[:, 5]` = the summary's `af`,
   A1 reads / site reads) and the somatic truth ids a label-1 tensor stands for.
3. `summarize.py`: replays `train.py`'s split with `data.block_split` (chr1 left out; chr2-22 in 20,000-node blocks,
   5 % to validation, seed 0: the `args.json` of all 36 local runs and the headers of the 6 RIS-only logs say so);
   training = label >= 0 in the training blocks, validation = evaluation label >= 0 in the validation blocks, both
   without the recipes' `--ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth`. The replayed
   per-class counts, the partial left-out counts and the validation off-reference counts are compared with every
   `train.log` (`runs_replay_check.tsv`). It then counts distinct somatic truth alleles (`data.somatic_truth`,
   `data.validation_truth`, `somatic.recall.tsv`, `somatic.graph.tsv`), the AF of the training tensors and the SBS6
   spectrum of the SNV truth alleles behind them (the trinucleotide around each VCF position, from the labels' GRCh38
   FASTA, each chromosome read once).
4. `make_e13_subsets.py`: nested subsets of the training blocks as `--exclude-candidates` files for E13.
5. `make_tables.py`: `tables.md`.

Commands (repository root `/scratch/jshen/Github/Pansoma`; jobs in `tmp/reviewer_cpu_20261013/E21_training_composition/`):

```bash
E=analysis/reviewer_response_20261011/experiments/E21_training_composition
J=tmp/reviewer_cpu_20261013/E21_training_composition
/wanglab/jshen/anaconda3/bin/python $E/collect_runs.py
sbatch $J/build_index.sbatch                                   # 387137, sets in $J/sets.txt
sbatch --dependency=afterany:387137 $J/collect.sbatch          # 387193: runs $E/collect.sh
# collect.sh = collect_runs.py; summarize.py; make_e13_subsets.py (HG008T_PacBio_d46, HG008T_PacBio); make_tables.py
# 2026-10-13 23:01, login node, from the 21 finished caches (summarize.py skips a set without a cache):
/wanglab/jshen/anaconda3/bin/python $E/collect_runs.py > /dev/null
/wanglab/jshen/anaconda3/bin/python $E/summarize.py      # "replay check: 38 of 42 runs match their train.log, 0 differ, 4 not checked"
/wanglab/jshen/anaconda3/bin/python $E/make_tables.py    # tables.md: 435 lines
```

Job states (2026-10-13, `squeue` / `sacct`; logs in `$J/logs/`):

| job | what | state |
|---|---|---|
| 387137 tasks 1-21 | index caches of the 18 reference-arm sets and the three COLO829T Illumina dilutions | COMPLETED (`logs/index_<task>.log`; per-kind builds up to 3,713 s, printed max RSS up to 5.5 GB) |
| 387137 tasks 22-33 | COLO829T fiberseq P75 / P50 / P25, ONT P75 / P50 / P25, HCC1395T Illumina / PacBio / ONT d46, HCC1954T / H1437T / H2009T PacBio d46 | PENDING, held (`JobHeldUser`); not released (waiting for the user's cluster cap) |
| 387193 | `collect.sbatch` (`collect.sh`, `--dependency=afterany:387137`) | PENDING, held (`JobHeldUser`) |
| 387270 | `E21_fastatest`: BeeGFS read-speed test behind `contexts()` in `summarize.py` (`logs/fastatest_387270.log`: 46,709,983 bytes in 43.8 s, 242,193,529 bytes in 213.5 s) | COMPLETED |

`make_e13_subsets.py` was not rerun at 23:01: the exclude files and `e13_subsets_*_SNV.tsv` are those written at
10:16 (`HG008T_PacBio`) and 10:32 (`HG008T_PacBio_d46`). `runs_training.*` came out byte-identical to the earlier run.

## 1. Every model's training data (R1.7; `tables.md` T1, `runs_training.tsv`)

36 single-set runs (12 per arm) and the 6 COLO829T purity joint runs. train / val = tensors per class (non / somatic /
germline) from `train.log`; left out = partial somatic matches removed by `--ignore-reasons` (training / validation);
replay = the split re-derived from the label files gives exactly these counts and left-out counts (yes; "-" = not
checked: the fiberseq / ONT joint runs, whose P75 / P50 / P25 caches are in the held array tasks); losses from the
per-epoch log lines; val somatic AP at the best epoch (the checkpoint); PASS thr = the checkpoint's validation
recall-0.9 threshold. The six `COLO829_*_h90` /
`COLO829_*_INDEL_small` rows are the COLO829T d9 models (run directories on RIS only; AP and threshold from the
benchmark `summary.json`).

| run | arm | params (M) | train non / somatic / germline | val non / somatic / germline | left out train / val | replay | epochs (best) | train loss best / last | val loss best / min (ep) / last | val somatic AP | PASS thr |
|---|---|---:|---|---|---|---|---|---|---|---:|---:|
| `HG008_Illumina_SNV_7fc7396_ddp4_b1024` | d9 | 50.4 | 1,620,829 / 7,134 / 275,839 | 103,012 / 469 / 14,737 | 1,720 / 97 | yes | 20/20 (4) | 0.188 / 0.022 | 0.240 / 0.220 (5) / 2.065 | 0.190 | 0.251 |
| `HG008_Illumina_SNV_d46_7fc7396` | d46 | 50.4 | 1,619,483 / 7,274 / 256,583 | 80,578 / 357 / 12,553 | 1,668 / 92 | yes | 20/20 (5) | 0.145 / 0.020 | 0.179 / 0.173 (6) / 1.385 | 0.275 | 0.324 |
| `HG008_Illumina_SNV_linear_7fc7396` | linear | 50.4 | 1,832,960 / 7,371 / 2,710,668 | 107,131 / 484 / 156,406 | 1,893 / 114 | yes | 20/20 (6) | 0.242 / 0.034 | 0.302 / 0.296 (4) / 2.801 | 0.011 | 0.064 |
| `HG008_Illumina_INDEL_7fc7396` | d9 | 8.7 | 27,729 / 3,121 / 61,258 | 2,363 / 175 / 3,471 | 1,458 / 91 | yes | 12/12 (3) | 0.689 / 0.548 | 0.825 / 0.739 (6) / 0.766 | 0.110 | 0.119 |
| `HG008_Illumina_INDEL_d46_7fc7396` | d46 | 8.7 | 20,551 / 3,383 / 55,260 | 1,767 / 162 / 2,622 | 1,517 / 79 | yes | 12/12 (3) | 0.715 / 0.566 | 0.864 / 0.832 (4) / 0.891 | 0.130 | 0.106 |
| `HG008_Illumina_INDEL_linear_7fc7396` | linear | 8.7 | 122,745 / 5,893 / 616,538 | 7,225 / 322 / 35,205 | 460 / 29 | yes | 12/12 (2) | 0.359 / 0.145 | 0.345 / 0.321 (5) / 0.581 | 0.026 | 0.048 |
| `HG008_PacBio_SNV_7fc7396` | d9 | 50.4 | 67,735 / 7,098 / 306,828 | 3,339 / 465 / 15,047 | 464 / 28 | yes | 20/20 (7) | 0.358 / 0.139 | 0.504 / 0.503 (5) / 2.053 | 0.227 | 0.314 |
| `HG008_PacBio_SNV_d46_7fc7396` | d46 | 50.4 | 71,573 / 7,227 / 270,175 | 4,509 / 350 / 12,631 | 464 / 28 | yes | 20/20 (7) | 0.304 / 0.130 | 0.445 / 0.350 (9) / 0.895 | 0.297 | 0.291 |
| `HG008_PacBio_SNV_linear_7fc7396` | linear | 50.4 | 165,119 / 7,370 / 2,745,963 | 13,161 / 477 / 157,349 | 369 / 14 | yes | 20/20 (7) | 0.333 / 0.047 | 0.420 / 0.356 (5) / 2.847 | 0.012 | 0.040 |
| `HG008_PacBio_INDEL_7fc7396` | d9 | 8.7 | 1,307,921 / 3,588 / 165,779 | 88,590 / 194 / 9,313 | 2,401 / 138 | yes | 12/12 (4) | 0.311 / 0.222 | 0.303 / 0.303 (4) / 0.366 | 0.093 | 0.040 |
| `HG008_PacBio_INDEL_d46_7fc7396` | d46 | 8.7 | 1,329,154 / 3,761 / 162,889 | 70,488 / 181 / 7,511 | 2,310 / 124 | yes | 12/12 (4) | 0.307 / 0.225 | 0.325 / 0.323 (6) / 0.423 | 0.129 | 0.027 |
| `HG008_PacBio_INDEL_linear_7fc7396` | linear | 8.7 | 1,326,439 / 6,224 / 682,475 | 77,410 / 349 / 38,369 | 643 / 36 | yes | 12/12 (4) | 0.203 / 0.135 | 0.199 / 0.194 (5) / 0.267 | 0.029 | 0.037 |
| `HG008_ONT_SNV_7fc7396` | d9 | 50.4 | 175,573 / 7,124 / 326,341 | 9,067 / 470 / 15,230 | 452 / 31 | yes | 20/20 (8) | 0.286 / 0.128 | 0.425 / 0.425 (8) / 1.464 | 0.210 | 0.449 |
| `HG008_ONT_SNV_d46_7fc7396` | d46 | 50.4 | 183,510 / 7,259 / 302,622 | 9,488 / 355 / 12,816 | 440 / 24 | yes | 20/20 (9) | 0.234 / 0.100 | 0.364 / 0.318 (6) / 1.115 | 0.284 | 0.231 |
| `HG008_ONT_SNV_linear_7fc7396` | linear | 50.4 | 295,571 / 7,391 / 2,786,375 | 19,218 / 479 / 158,792 | 386 / 17 | yes | 20/20 (8) | 0.291 / 0.060 | 0.434 / 0.412 (5) / 2.604 | 0.012 | 0.042 |
| `HG008_ONT_INDEL_7fc7396` | d9 | 8.7 | 1,696,693 / 3,203 / 190,869 | 113,648 / 181 / 10,756 | 2,343 / 150 | yes | 12/12 (5) | 0.345 / 0.240 | 0.346 / 0.337 (7) / 0.437 | 0.055 | 0.018 |
| `HG008_ONT_INDEL_d46_7fc7396` | d46 | 8.7 | 1,710,450 / 3,249 / 176,041 | 91,378 / 156 / 8,081 | 2,468 / 135 | yes | 12/12 (5) | 0.345 / 0.254 | 0.351 / 0.351 (5) / 0.504 | 0.080 | 0.014 |
| `HG008_ONT_INDEL_linear_7fc7396` | linear | 8.7 | 1,683,055 / 4,242 / 653,473 | 97,780 / 251 / 36,938 | 1,489 / 63 | yes | 12/12 (5) | 0.268 / 0.177 | 0.282 / 0.274 (6) / 0.425 | 0.014 | 0.008 |
| `COLO829_Illumina_SNV_h90` | d9 | 50.4 | 807,850 / 33,782 / 301,623 | 49,537 / 2,060 / 17,069 | 374 / 26 | yes | 20/20 (6) | 0.178 / 0.014 | 0.214 / 0.212 (5) / 1.313 | 0.720 | 0.330 |
| `COLO829_Illumina_SNV_d46_7fc7396` | d46 | 50.4 | 803,764 / 34,110 / 280,542 | 40,055 / 1,737 / 13,382 | 356 / 15 | yes | 20/20 (7) | 0.141 / 0.019 | 0.178 / 0.169 (5) / 0.807 | 0.812 | 0.640 |
| `COLO829_Illumina_SNV_linear_7fc7396` | linear | 50.4 | 1,115,204 / 34,341 / 2,892,361 | 65,265 / 1,775 / 163,915 | 426 / 18 | yes | 20/20 (7) | 0.227 / 0.037 | 0.250 / 0.247 (5) / 2.207 | 0.227 | 0.081 |
| `COLO829_Illumina_SNV_joint_7fc7396_gpu4` | d9 joint P100+P75+P50+P25 | 50.4 | 3,047,060 / 125,916 / 1,356,321 | 187,856 / 7,711 / 76,352 | 1,141 / 59 | yes | 20/20 (7) | 0.075 / 0.002 | 0.132 / 0.118 (3) / 0.621 | 0.961 | 0.804 |
| `COLO829_Illumina_INDEL_small` | d9 | 8.7 | 30,775 / 971 / 67,027 | 2,648 / 65 / 4,023 | 218 / 11 | yes | 12/12 (5) | 0.570 / 0.292 | 0.715 / 0.689 (4) / 1.188 | 0.065 | 0.020 |
| `COLO829_Illumina_INDEL_d46_7fc7396` | d46 | 8.7 | 23,386 / 1,048 / 61,110 | 1,928 / 52 / 2,957 | 173 / 14 | yes | 12/12 (7) | 0.581 / 0.492 | 0.775 / 0.771 (1) / 0.838 | 0.057 | 0.022 |
| `COLO829_Illumina_INDEL_linear_7fc7396` | linear | 8.7 | 120,648 / 1,555 / 655,606 | 6,866 / 91 / 37,569 | 92 / 2 | yes | 12/12 (4) | 0.291 / 0.132 | 0.309 / 0.294 (5) / 0.725 | 0.007 | 0.011 |
| `COLO829_Illumina_INDEL_joint_7fc7396` | d9 joint P100+P75+P50+P25 | 8.7 | 127,515 / 3,439 / 297,626 | 11,193 / 233 / 17,725 | 666 / 29 | yes | 12/12 (6) | 0.403 / 0.282 | 0.591 / 0.579 (4) / 0.802 | 0.576 | 0.025 |
| `COLO829_fiberseq_SNV_h90` | d9 | 50.4 | 111,360 / 33,669 / 350,948 | 6,285 / 2,066 / 17,676 | 203 / 9 | yes | 20/20 (12) | 0.242 / 0.105 | 0.487 / 0.386 (2) / 1.313 | 0.711 | 0.394 |
| `COLO829_fiberseq_SNV_d46_7fc7396` | d46 | 50.4 | 152,721 / 33,998 / 320,431 | 8,016 / 1,725 / 13,912 | 183 / 9 | yes | 20/20 (9) | 0.231 / 0.134 | 0.341 / 0.288 (11) / 0.402 | 0.822 | 0.481 |
| `COLO829_fiberseq_SNV_linear_7fc7396` | linear | 50.4 | 221,989 / 34,100 / 2,950,677 | 13,452 / 1,774 / 165,480 | 107 / 4 | yes | 20/20 (10) | 0.258 / 0.068 | 0.347 / 0.286 (7) / 1.680 | 0.256 | 0.085 |
| `COLO829_fiberseq_SNV_joint_7fc7396_gpu4` | d9 joint P100+P75+P50+P25 | 50.4 | 421,331 / 125,734 / 1,523,106 | 22,975 / 7,713 / 78,121 | 699 / 33 | - | 20/20 (5) | 0.228 / 0.013 | 0.463 / 0.348 (2) / 3.275 | 0.789 | 0.549 |
| `COLO829_fiberseq_INDEL_small` | d9 | 8.7 | 1,124,518 / 1,078 / 181,407 | 78,510 / 73 / 10,651 | 374 / 14 | yes | 12/12 (7) | 0.286 / 0.224 | 0.327 / 0.324 (6) / 0.432 | 0.035 | 0.006 |
| `COLO829_fiberseq_INDEL_d46_7fc7396` | d46 | 8.7 | 1,149,976 / 1,141 / 180,805 | 63,255 / 58 / 8,477 | 349 / 24 | yes | 12/12 (2) | 0.345 / 0.245 | 0.356 / 0.351 (7) / 0.428 | 0.039 | 0.012 |
| `COLO829_fiberseq_INDEL_linear_7fc7396` | linear | 8.7 | 1,096,037 / 1,641 / 730,009 | 62,991 / 92 / 41,613 | 124 / 3 | yes | 12/12 (6) | 0.165 / 0.102 | 0.178 / 0.178 (6) / 0.271 | 0.011 | 0.008 |
| `COLO829_fiberseq_INDEL_joint_7fc7396` | d9 joint P100+P75+P50+P25 | 8.7 | 5,835,645 / 3,965 / 759,274 | 394,134 / 263 / 44,526 | 1,191 / 32 | - | 12/12 (3) | 0.266 / 0.172 | 0.298 / 0.289 (4) / 0.466 | 0.052 | 0.007 |
| `COLO829_ONT_SNV_h90` | d9 | 50.4 | 363,230 / 33,472 / 381,159 | 21,065 / 2,037 / 18,504 | 222 / 6 | yes | 20/20 (8) | 0.213 / 0.050 | 0.354 / 0.332 (6) / 1.610 | 0.690 | 0.463 |
| `COLO829_ONT_SNV_d46_7fc7396` | d46 | 50.4 | 476,300 / 33,779 / 355,824 | 29,549 / 1,717 / 14,494 | 200 / 11 | yes | 20/20 (10) | 0.167 / 0.068 | 0.242 / 0.229 (9) / 0.690 | 0.794 | 0.568 |
| `COLO829_ONT_SNV_linear_7fc7396` | linear | 50.4 | 505,321 / 33,852 / 3,008,184 | 28,773 / 1,750 / 167,861 | 122 / 5 | yes | 20/20 (6) | 0.272 / 0.052 | 0.305 / 0.302 (7) / 2.111 | 0.170 | 0.092 |
| `COLO829_ONT_SNV_joint_7fc7396_gpu4` | d9 joint P100+P75+P50+P25 | 50.4 | 1,280,987 / 124,070 / 1,665,260 | 75,779 / 7,614 / 81,641 | 746 / 25 | - | 20/20 (7) | 0.165 / 0.006 | 0.456 / 0.349 (6) / 2.161 | 0.801 | 0.500 |
| `COLO829_ONT_INDEL_small` | d9 | 8.7 | 3,285,172 / 999 / 220,610 | 218,747 / 65 / 12,957 | 555 / 24 | yes | 12/12 (5) | 0.233 / 0.162 | 0.245 / 0.238 (6) / 0.331 | 0.024 | 0.012 |
| `COLO829_ONT_INDEL_d46_7fc7396` | d46 | 8.7 | 3,289,717 / 1,018 / 212,429 | 173,128 / 53 / 10,104 | 515 / 24 | yes | 12/12 (4) | 0.250 / 0.185 | 0.260 / 0.260 (4) / 0.339 | 0.046 | 0.013 |
| `COLO829_ONT_INDEL_linear_7fc7396` | linear | 8.7 | 3,204,105 / 1,276 / 688,095 | 185,957 / 66 / 39,055 | 329 / 16 | yes | 12/12 (5) | 0.178 / 0.113 | 0.184 / 0.184 (5) / 0.281 | 0.007 | 0.005 |
| `COLO829_ONT_INDEL_joint_7fc7396` | d9 joint P100+P75+P50+P25 | 8.7 | 10,987,344 / 3,747 / 897,514 | 742,249 / 248 / 52,698 | 1,694 / 67 | - | 12/12 (6) | 0.196 / 0.116 | 0.300 / 0.262 (3) / 0.491 | 0.044 | 0.000 |

Reading:

- Parameters: 50.4 M (SNV), 8.7 M (INDEL), printed by `train.py` for every run.
- Graph arms (d9, d46): HG008T SNV 7,098-7,274 somatic, 256,583-326,341 germline, 67,735-1,620,829 non training
  tensors; COLO829T SNV 33,472-34,110 somatic. INDEL: HG008T 3,121-3,761 somatic, COLO829T 971-1,141. GRCh38-only arm:
  more somatic INDELs reach a tensor (HG008T 4,242-6,224) and 8-11x the germline SNV tensors (2.71-3.01 M).
- Validation (5 % of the node blocks): 350-484 somatic SNV tensors for HG008T and 1,717-2,066 for COLO829T; 52-194 somatic
  INDEL tensors on the graph arms. The PASS thresholds of the INDEL models and of the HG008T models rest on these few
  hundred validation positives.
- Joint models (COLO829T P100 + P75 + P50 + P25, d9): 124,070-125,916 somatic SNV and 3,439-3,965 somatic INDEL training
  tensors, 3.7x (SNV) and 3.5-3.8x (INDEL) the single-set P100 model of the platform, but not more alleles: the Illumina joint SNV positives are 33,797
  distinct truth alleles against 33,782 for P100 alone (`runs_replay_check.tsv`, section 2).
- SNV: the validation loss is lowest at epoch 2-11 and ends 1.4-9.5x higher at epoch 20, while the training loss falls
  to 0.002-0.139. INDEL: lowest at epoch 1-7, ends 1.0-2.5x higher. Best epochs (validation somatic F1): SNV 4-12 of 20,
  INDEL 2-7 of 12.
- The COLO829T joint Illumina models have validation somatic AP 0.961 (SNV) and 0.576 (INDEL) against 0.720-0.812 and
  0.057-0.065 for the single-set models: the in-silico dilutions mix tumor reads with COLO829BL reads of higher base
  quality, and that read-source cue inflates their validation (benchmark README, `colo829t_purity_joint_7fc7396_20261012`).

## 2. Somatic truth alleles behind the training tensors (R1.7, R1.8; `tables.md` T3-T4, `somatic_truth_alleles.tsv`)

Computed for the 21 sets with an index cache (the 18 reference-arm sets and the three COLO829T Illumina dilutions); the
fiberseq / ONT dilutions, HCC1395T and CASTLE are pending (held array tasks 22-33). Truth in BED = PASS somatic truth
alleles of the kind on chr2-22 in the labels' confident region; training region = outside the validation blocks
(`data.validation_truth`); behind training tensors = distinct truth alleles a training tensor stands for; status of the
training-region alleles from `somatic.recall.tsv` (rep = A1 of a tensor; other = another allele of a site; filtered =
removed by the build's AF 0.08 / 3-ALT-read filters; none = no candidate at all; why, for HG008T, is in
`analysis/hg008_somatic_miss_20260926`, not re-checked here).

| set | kind | truth in BED | in validation region | training region | behind training tensors | % of training region | rep | other | filtered | none |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HG008T_Illumina | SNV | 7,846 | 487 | 7,359 | 7,134 | 95.3 % | 7,028 | 3 | 62 | 266 |
| HG008T_Illumina | INDEL | 7,394 | 404 | 6,990 | 3,115 | 42.6 % | 2,987 | 71 | 723 | 3,209 |
| HG008T_Illumina_d46 | SNV | 7,846 | 379 | 7,467 | 7,274 | 95.9 % | 7,164 | 3 | 71 | 229 |
| HG008T_Illumina_d46 | INDEL | 7,394 | 373 | 7,021 | 3,377 | 45.8 % | 3,219 | 92 | 655 | 3,055 |
| HG008T_Illumina_linear | SNV | 7,846 | 488 | 7,358 | 7,371 | 98.5 % | 7,250 | 6 | 50 | 52 |
| HG008T_Illumina_linear | INDEL | 7,394 | 412 | 6,982 | 5,880 | 79.6 % | 5,561 | 568 | 581 | 272 |
| HG008T_PacBio | SNV | 7,846 | 483 | 7,363 | 7,098 | 94.8 % | 6,992 | 2 | 48 | 321 |
| HG008T_PacBio | INDEL | 7,394 | 399 | 6,995 | 3,557 | 48.5 % | 3,399 | 318 | 840 | 2,438 |
| HG008T_PacBio_d46 | SNV | 7,846 | 381 | 7,465 | 7,227 | 95.3 % | 7,115 | 2 | 45 | 303 |
| HG008T_PacBio_d46 | INDEL | 7,394 | 372 | 7,022 | 3,728 | 50.5 % | 3,551 | 330 | 565 | 2,576 |
| HG008T_PacBio_linear | SNV | 7,846 | 488 | 7,358 | 7,370 | 98.5 % | 7,245 | 7 | 46 | 60 |
| HG008T_PacBio_linear | INDEL | 7,394 | 414 | 6,980 | 6,185 | 83.7 % | 5,845 | 895 | 124 | 116 |
| HG008T_ONT | SNV | 7,846 | 484 | 7,362 | 7,124 | 95.1 % | 7,020 | 4 | 61 | 277 |
| HG008T_ONT | INDEL | 7,394 | 399 | 6,995 | 3,163 | 43.1 % | 3,023 | 709 | 1,272 | 1,991 |
| HG008T_ONT_d46 | SNV | 7,846 | 382 | 7,464 | 7,259 | 95.7 % | 7,145 | 2 | 53 | 264 |
| HG008T_ONT_d46 | INDEL | 7,394 | 370 | 7,024 | 3,214 | 43.6 % | 3,071 | 793 | 876 | 2,284 |
| HG008T_ONT_linear | SNV | 7,846 | 488 | 7,358 | 7,391 | 98.7 % | 7,265 | 6 | 47 | 40 |
| HG008T_ONT_linear | INDEL | 7,394 | 414 | 6,980 | 4,220 | 57.4 % | 4,005 | 1,652 | 1,051 | 272 |
| COLO829T_Illumina | SNV | 39,592 | 2,274 | 37,318 | 33,782 | 90.4 % | 33,799 | 1 | 2,992 | 526 |
| COLO829T_Illumina | INDEL | 1,781 | 103 | 1,678 | 967 | 57.3 % | 964 | 12 | 112 | 590 |
| COLO829T_Illumina_d46 | SNV | 39,592 | 1,913 | 37,679 | 34,110 | 90.5 % | 34,129 | 1 | 2,935 | 614 |
| COLO829T_Illumina_d46 | INDEL | 1,781 | 79 | 1,702 | 1,042 | 60.9 % | 1,041 | 14 | 89 | 558 |
| COLO829T_Illumina_linear | SNV | 39,592 | 2,030 | 37,562 | 34,341 | 91.4 % | 34,341 | 1 | 2,542 | 678 |
| COLO829T_Illumina_linear | INDEL | 1,781 | 100 | 1,681 | 1,544 | 91.7 % | 1,542 | 24 | 66 | 49 |
| COLO829T_fiberseq | SNV | 39,592 | 2,265 | 37,327 | 33,669 | 90.0 % | 33,688 | 0 | 3,166 | 473 |
| COLO829T_fiberseq | INDEL | 1,781 | 101 | 1,680 | 1,072 | 63.6 % | 1,071 | 33 | 206 | 370 |
| COLO829T_fiberseq_d46 | SNV | 39,592 | 1,899 | 37,693 | 33,998 | 90.1 % | 34,014 | 0 | 3,130 | 549 |
| COLO829T_fiberseq_d46 | INDEL | 1,781 | 80 | 1,701 | 1,134 | 66.7 % | 1,134 | 35 | 138 | 394 |
| COLO829T_fiberseq_linear | SNV | 39,592 | 2,031 | 37,561 | 34,100 | 90.8 % | 34,100 | 1 | 2,902 | 558 |
| COLO829T_fiberseq_linear | INDEL | 1,781 | 98 | 1,683 | 1,627 | 96.7 % | 1,627 | 20 | 24 | 12 |
| COLO829T_ONT | SNV | 39,592 | 2,244 | 37,348 | 33,472 | 89.5 % | 33,513 | 1 | 3,364 | 470 |
| COLO829T_ONT | INDEL | 1,781 | 101 | 1,680 | 986 | 58.6 % | 985 | 169 | 229 | 297 |
| COLO829T_ONT_d46 | SNV | 39,592 | 1,916 | 37,676 | 33,779 | 89.6 % | 33,793 | 1 | 3,430 | 452 |
| COLO829T_ONT_d46 | INDEL | 1,781 | 80 | 1,701 | 1,007 | 59.2 % | 1,007 | 207 | 158 | 329 |
| COLO829T_ONT_linear | SNV | 39,592 | 2,031 | 37,561 | 33,852 | 90.1 % | 33,852 | 0 | 3,404 | 305 |
| COLO829T_ONT_linear | INDEL | 1,781 | 98 | 1,683 | 1,269 | 75.4 % | 1,269 | 295 | 115 | 4 |
| COLO829T_Illumina_P75 | SNV | 39,592 | 2,273 | 37,319 | 33,098 | 88.5 % | 33,115 | 1 | 3,533 | 670 |
| COLO829T_Illumina_P75 | INDEL | 1,781 | 104 | 1,677 | 953 | 56.4 % | 949 | 9 | 122 | 597 |
| COLO829T_Illumina_P50 | SNV | 39,592 | 2,271 | 37,321 | 32,073 | 85.8 % | 32,089 | 0 | 4,249 | 983 |
| COLO829T_Illumina_P50 | INDEL | 1,781 | 105 | 1,676 | 894 | 53.0 % | 889 | 13 | 147 | 627 |
| COLO829T_Illumina_P25 | SNV | 39,592 | 2,259 | 37,333 | 26,963 | 72.1 % | 26,991 | 1 | 8,526 | 1,815 |
| COLO829T_Illumina_P25 | INDEL | 1,781 | 104 | 1,677 | 614 | 36.4 % | 611 | 12 | 384 | 670 |

Unions over platforms and arms, and the joint models' sets (distinct alleles by VCF chrom / pos / REF / ALT;
`somatic_truth_unions.tsv`, which also has the per-platform rows). The fiberseq / ONT "joint" rows of
`somatic_truth_unions.tsv` and `tables.md` T4 contain only the P100 set until the held caches exist:

| sample | kind | group | distinct truth alleles |
|---|---|---|---:|
| HG008T | SNV | all platforms, d9 | 7,203 |
| HG008T | SNV | all platforms, d46 | 7,337 |
| HG008T | SNV | all platforms, linear | 7,448 |
| HG008T | SNV | all platforms, all three arms | 7,935 |
| HG008T | INDEL | all platforms, d9 | 4,037 |
| HG008T | INDEL | all platforms, d46 | 4,154 |
| HG008T | INDEL | all platforms, linear | 6,706 |
| HG008T | INDEL | all platforms, all three arms | 7,107 |
| COLO829T | SNV | all platforms, d9 | 34,508 |
| COLO829T | SNV | all platforms, d46 | 34,853 |
| COLO829T | SNV | all platforms, linear | 34,971 |
| COLO829T | SNV | all platforms, all three arms | 36,880 |
| COLO829T | INDEL | all platforms, d9 | 1,179 |
| COLO829T | INDEL | all platforms, d46 | 1,236 |
| COLO829T | INDEL | all platforms, linear | 1,653 |
| COLO829T | INDEL | all platforms, all three arms | 1,730 |
| COLO829T | SNV | Illumina d9 | 33,782 |
| COLO829T | SNV | Illumina d9 P75 | 33,098 |
| COLO829T | SNV | Illumina d9 P50 | 32,073 |
| COLO829T | SNV | Illumina d9 P25 | 26,963 |
| COLO829T | SNV | Illumina joint P100+P75+P50+P25 | 33,797 |
| COLO829T | SNV | fiberseq d9 | 33,669 |
| COLO829T | SNV | fiberseq joint P100+P75+P50+P25 (only P100 read: P75 / P50 / P25 caches pending) | 33,669 |
| COLO829T | SNV | ONT d9 | 33,472 |
| COLO829T | SNV | ONT joint P100+P75+P50+P25 (only P100 read: P75 / P50 / P25 caches pending) | 33,472 |
| COLO829T | INDEL | Illumina d9 | 967 |
| COLO829T | INDEL | Illumina d9 P75 | 953 |
| COLO829T | INDEL | Illumina d9 P50 | 894 |
| COLO829T | INDEL | Illumina d9 P25 | 614 |
| COLO829T | INDEL | Illumina joint P100+P75+P50+P25 | 973 |
| COLO829T | INDEL | fiberseq d9 | 1,072 |
| COLO829T | INDEL | fiberseq joint P100+P75+P50+P25 (only P100 read: P75 / P50 / P25 caches pending) | 1,072 |
| COLO829T | INDEL | ONT d9 | 986 |
| COLO829T | INDEL | ONT joint P100+P75+P50+P25 (only P100 read: P75 / P50 / P25 caches pending) | 986 |

- SNV: on the pangenome arms (d9, d46) 94.8-95.9 % (HG008T) and 89.5-90.5 % (COLO829T) of the training-region truth
  stands behind a training tensor; on GRCh38-only 98.5-98.7 % and 90.1-91.4 %. The COLO829T shortfall is mostly the
  build's AF / read-count filter (2,542-3,430 alleles over the nine P100 sets: subclonal SNVs below AF 0.08 or with fewer
  than 3 ALT reads), not the graph. HG008T loses 229-321 SNVs without a candidate on the pangenome arms and 40-60 on
  GRCh38-only.
- INDEL: on the pangenome arms 42.6-50.5 % (HG008T) and 57.3-66.7 % (COLO829T); 1,991-3,209 of the ~7,000 HG008T
  training-region INDELs have no candidate and 565-1,272 are filtered. On GRCh38-only the same truth reaches 57.4-83.7 %
  (HG008T) and 75.4-96.7 % (COLO829T), and only 116-272 (HG008T) and 4-49 (COLO829T) INDELs have no candidate. So most
  of the missing INDEL positives are lost on the pangenome arms, before labelling. The earlier miss analysis attributes
  most HG008T INDEL misses to alleles that are graph paths or that the alignments write differently
  (`analysis/hg008_somatic_miss_20260926`; see also `analysis/graph_absorbed_somatic_20261001`; R1.4 / E5; not re-checked
  here). HG008T ONT on GRCh38-only (57.4 %) is the exception among the linear sets: 1,652 of its INDELs are another allele
  of a site and 1,051 are filtered (long-read INDEL errors split the site's reads) [interpretation, unverified].
- Dilutions (COLO829T Illumina d9): the share falls from 90.4 % (P100) to 88.5 / 85.8 / 72.1 % (P75 / P50 / P25) for SNVs
  and from 57.3 % to 56.4 / 53.0 / 36.4 % for INDELs, because the build filters remove more of the diluted alleles
  (filtered SNVs 2,992 at P100, 3,533 / 4,249 / 8,526 at P75 / P50 / P25).
- The positives of the three platforms are nearly the same alleles (HG008T SNV: 7,098-7,134 per platform on d9, 7,203 in
  the union; 7,935 over all platforms and arms). The dilution sets add almost no alleles to the joint model: Illumina
  joint 33,797 SNV alleles against 33,782 at P100 (INDEL 973 against 967). The fiberseq / ONT joint unions: pending.

## 3. VAF of the training tensors (R1.8; `tables.md` T5-T6, `vaf_quantiles.tsv`, `vaf_hist.tsv`)

AF = the A1 allele fraction of the tensor (A1 reads / site reads, `af` in `variant_summary.ndjson`; `scalars[:, 5]` of the
index). Training tensors only (chr2-22, training blocks, the recipes' split). Somatic and germline rows (non rows in
`vaf_quantiles.tsv`):

| set | kind | class | n | 5 % | 25 % | median | 75 % | 95 % | AF < 0.20 | AF >= 0.40 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HG008T_Illumina | SNV | somatic | 7,134 | 0.317 | 0.471 | 0.537 | 0.990 | 1.000 | 0.6 % | 89.1 % |
| HG008T_Illumina | SNV | germline | 275,839 | 0.301 | 0.475 | 0.564 | 1.000 | 1.000 | 2.1 % | 88.2 % |
| HG008T_Illumina | INDEL | somatic | 3,121 | 0.169 | 0.332 | 0.431 | 0.690 | 0.885 | 7.5 % | 58.3 % |
| HG008T_Illumina | INDEL | germline | 61,258 | 0.109 | 0.238 | 0.391 | 0.587 | 0.912 | 19.2 % | 48.4 % |
| HG008T_Illumina_d46 | SNV | somatic | 7,274 | 0.318 | 0.471 | 0.534 | 0.990 | 1.000 | 0.5 % | 89.1 % |
| HG008T_Illumina_d46 | SNV | germline | 256,583 | 0.299 | 0.473 | 0.558 | 1.000 | 1.000 | 2.2 % | 88.0 % |
| HG008T_Illumina_d46 | INDEL | somatic | 3,383 | 0.160 | 0.326 | 0.428 | 0.667 | 0.883 | 8.1 % | 57.9 % |
| HG008T_Illumina_d46 | INDEL | germline | 55,260 | 0.121 | 0.276 | 0.418 | 0.629 | 0.916 | 14.3 % | 53.7 % |
| HG008T_Illumina_linear | SNV | somatic | 7,371 | 0.308 | 0.469 | 0.534 | 0.989 | 1.000 | 1.0 % | 88.3 % |
| HG008T_Illumina_linear | SNV | germline | 2,710,668 | 0.377 | 0.527 | 0.991 | 1.000 | 1.000 | 0.7 % | 94.1 % |
| HG008T_Illumina_linear | INDEL | somatic | 5,893 | 0.140 | 0.290 | 0.404 | 0.602 | 0.848 | 11.1 % | 51.2 % |
| HG008T_Illumina_linear | INDEL | germline | 616,538 | 0.155 | 0.350 | 0.535 | 0.786 | 0.922 | 8.2 % | 67.4 % |
| HG008T_PacBio | SNV | somatic | 7,098 | 0.318 | 0.466 | 0.539 | 0.982 | 1.000 | 0.5 % | 88.6 % |
| HG008T_PacBio | SNV | germline | 306,828 | 0.300 | 0.477 | 0.604 | 1.000 | 1.000 | 2.2 % | 88.6 % |
| HG008T_PacBio | INDEL | somatic | 3,588 | 0.122 | 0.333 | 0.453 | 0.615 | 0.927 | 10.7 % | 62.6 % |
| HG008T_PacBio | INDEL | germline | 165,779 | 0.088 | 0.132 | 0.227 | 0.442 | 0.948 | 44.2 % | 28.6 % |
| HG008T_PacBio_d46 | SNV | somatic | 7,227 | 0.319 | 0.466 | 0.538 | 0.982 | 1.000 | 0.4 % | 88.7 % |
| HG008T_PacBio_d46 | SNV | germline | 270,175 | 0.309 | 0.476 | 0.586 | 1.000 | 1.000 | 1.8 % | 89.0 % |
| HG008T_PacBio_d46 | INDEL | somatic | 3,761 | 0.149 | 0.348 | 0.454 | 0.608 | 0.915 | 8.4 % | 63.6 % |
| HG008T_PacBio_d46 | INDEL | germline | 162,889 | 0.089 | 0.136 | 0.230 | 0.443 | 0.948 | 43.4 % | 28.7 % |
| HG008T_PacBio_linear | SNV | somatic | 7,370 | 0.304 | 0.462 | 0.535 | 0.980 | 1.000 | 1.1 % | 87.5 % |
| HG008T_PacBio_linear | SNV | germline | 2,745,963 | 0.378 | 0.528 | 0.982 | 1.000 | 1.000 | 0.7 % | 94.0 % |
| HG008T_PacBio_linear | INDEL | somatic | 6,224 | 0.198 | 0.350 | 0.448 | 0.586 | 0.865 | 5.0 % | 63.2 % |
| HG008T_PacBio_linear | INDEL | germline | 682,475 | 0.170 | 0.400 | 0.566 | 0.869 | 1.000 | 6.6 % | 75.2 % |
| HG008T_ONT | SNV | somatic | 7,124 | 0.288 | 0.437 | 0.535 | 0.931 | 1.000 | 1.2 % | 83.1 % |
| HG008T_ONT | SNV | germline | 326,341 | 0.247 | 0.444 | 0.575 | 0.938 | 1.000 | 3.3 % | 83.2 % |
| HG008T_ONT | INDEL | somatic | 3,203 | 0.115 | 0.203 | 0.311 | 0.463 | 0.906 | 23.4 % | 34.2 % |
| HG008T_ONT | INDEL | germline | 190,869 | 0.093 | 0.133 | 0.189 | 0.303 | 0.800 | 53.0 % | 17.6 % |
| HG008T_ONT_d46 | SNV | somatic | 7,259 | 0.289 | 0.436 | 0.533 | 0.926 | 1.000 | 1.2 % | 82.8 % |
| HG008T_ONT_d46 | SNV | germline | 302,622 | 0.217 | 0.431 | 0.551 | 0.927 | 1.000 | 4.2 % | 80.7 % |
| HG008T_ONT_d46 | INDEL | somatic | 3,249 | 0.125 | 0.210 | 0.312 | 0.462 | 0.900 | 21.1 % | 34.0 % |
| HG008T_ONT_d46 | INDEL | germline | 176,041 | 0.097 | 0.143 | 0.200 | 0.320 | 0.800 | 48.3 % | 18.5 % |
| HG008T_ONT_linear | SNV | somatic | 7,391 | 0.271 | 0.431 | 0.531 | 0.923 | 1.000 | 1.8 % | 81.6 % |
| HG008T_ONT_linear | SNV | germline | 2,786,375 | 0.345 | 0.522 | 0.930 | 0.978 | 1.000 | 1.0 % | 91.4 % |
| HG008T_ONT_linear | INDEL | somatic | 4,242 | 0.140 | 0.219 | 0.316 | 0.460 | 0.867 | 18.4 % | 33.9 % |
| HG008T_ONT_linear | INDEL | germline | 653,473 | 0.133 | 0.240 | 0.413 | 0.793 | 0.976 | 16.0 % | 51.8 % |
| COLO829T_Illumina | SNV | somatic | 33,782 | 0.159 | 0.377 | 0.507 | 0.664 | 1.000 | 6.7 % | 72.4 % |
| COLO829T_Illumina | SNV | germline | 301,623 | 0.274 | 0.433 | 0.540 | 0.886 | 1.000 | 2.2 % | 79.5 % |
| COLO829T_Illumina | INDEL | somatic | 971 | 0.159 | 0.269 | 0.378 | 0.506 | 0.851 | 11.2 % | 44.8 % |
| COLO829T_Illumina | INDEL | germline | 67,027 | 0.106 | 0.215 | 0.354 | 0.513 | 0.876 | 22.0 % | 42.0 % |
| COLO829T_Illumina_d46 | SNV | somatic | 34,110 | 0.160 | 0.378 | 0.508 | 0.664 | 1.000 | 6.7 % | 72.6 % |
| COLO829T_Illumina_d46 | SNV | germline | 280,542 | 0.274 | 0.430 | 0.536 | 0.761 | 1.000 | 2.2 % | 79.2 % |
| COLO829T_Illumina_d46 | INDEL | somatic | 1,048 | 0.162 | 0.264 | 0.375 | 0.504 | 0.844 | 11.1 % | 44.9 % |
| COLO829T_Illumina_d46 | INDEL | germline | 61,110 | 0.118 | 0.248 | 0.378 | 0.535 | 0.877 | 16.6 % | 46.0 % |
| COLO829T_Illumina_linear | SNV | somatic | 34,341 | 0.160 | 0.377 | 0.507 | 0.661 | 1.000 | 6.7 % | 72.4 % |
| COLO829T_Illumina_linear | SNV | germline | 2,892,361 | 0.314 | 0.500 | 0.964 | 1.000 | 1.000 | 0.7 % | 87.9 % |
| COLO829T_Illumina_linear | INDEL | somatic | 1,555 | 0.165 | 0.265 | 0.370 | 0.500 | 0.827 | 10.7 % | 43.1 % |
| COLO829T_Illumina_linear | INDEL | germline | 655,606 | 0.151 | 0.313 | 0.477 | 0.738 | 0.912 | 9.4 % | 61.6 % |
| COLO829T_fiberseq | SNV | somatic | 33,669 | 0.161 | 0.376 | 0.502 | 0.658 | 1.000 | 6.5 % | 72.7 % |
| COLO829T_fiberseq | SNV | germline | 350,948 | 0.279 | 0.451 | 0.564 | 0.944 | 1.000 | 2.6 % | 80.9 % |
| COLO829T_fiberseq | INDEL | somatic | 1,078 | 0.133 | 0.325 | 0.424 | 0.544 | 0.877 | 9.0 % | 56.4 % |
| COLO829T_fiberseq | INDEL | germline | 181,407 | 0.088 | 0.127 | 0.215 | 0.416 | 0.889 | 46.8 % | 26.3 % |
| COLO829T_fiberseq_d46 | SNV | somatic | 33,998 | 0.161 | 0.376 | 0.502 | 0.658 | 1.000 | 6.5 % | 72.8 % |
| COLO829T_fiberseq_d46 | SNV | germline | 320,431 | 0.287 | 0.452 | 0.560 | 0.970 | 1.000 | 2.1 % | 81.3 % |
| COLO829T_fiberseq_d46 | INDEL | somatic | 1,141 | 0.132 | 0.321 | 0.423 | 0.545 | 0.877 | 8.2 % | 55.4 % |
| COLO829T_fiberseq_d46 | INDEL | germline | 180,805 | 0.089 | 0.130 | 0.220 | 0.418 | 0.893 | 45.9 % | 26.5 % |
| COLO829T_fiberseq_linear | SNV | somatic | 34,100 | 0.161 | 0.378 | 0.502 | 0.656 | 1.000 | 6.5 % | 72.9 % |
| COLO829T_fiberseq_linear | SNV | germline | 2,950,677 | 0.318 | 0.500 | 0.952 | 0.993 | 1.000 | 0.7 % | 88.1 % |
| COLO829T_fiberseq_linear | INDEL | somatic | 1,641 | 0.222 | 0.336 | 0.432 | 0.551 | 0.881 | 3.8 % | 58.2 % |
| COLO829T_fiberseq_linear | INDEL | germline | 730,009 | 0.167 | 0.359 | 0.513 | 0.788 | 0.987 | 7.0 % | 69.1 % |
| COLO829T_ONT | SNV | somatic | 33,472 | 0.151 | 0.346 | 0.460 | 0.598 | 0.926 | 7.3 % | 65.2 % |
| COLO829T_ONT | SNV | germline | 381,159 | 0.184 | 0.377 | 0.500 | 0.707 | 0.951 | 5.6 % | 71.5 % |
| COLO829T_ONT | INDEL | somatic | 999 | 0.114 | 0.195 | 0.274 | 0.411 | 0.717 | 25.4 % | 26.5 % |
| COLO829T_ONT | INDEL | germline | 220,610 | 0.088 | 0.126 | 0.174 | 0.260 | 0.600 | 60.1 % | 12.9 % |
| COLO829T_ONT_d46 | SNV | somatic | 33,779 | 0.150 | 0.344 | 0.457 | 0.595 | 0.920 | 7.4 % | 64.5 % |
| COLO829T_ONT_d46 | SNV | germline | 355,824 | 0.173 | 0.364 | 0.484 | 0.667 | 0.940 | 6.2 % | 68.9 % |
| COLO829T_ONT_d46 | INDEL | somatic | 1,018 | 0.121 | 0.197 | 0.274 | 0.404 | 0.698 | 25.5 % | 25.4 % |
| COLO829T_ONT_d46 | INDEL | germline | 212,429 | 0.091 | 0.132 | 0.179 | 0.267 | 0.589 | 58.0 % | 13.0 % |
| COLO829T_ONT_linear | SNV | somatic | 33,852 | 0.153 | 0.348 | 0.461 | 0.597 | 0.927 | 7.2 % | 65.5 % |
| COLO829T_ONT_linear | SNV | germline | 3,008,184 | 0.275 | 0.454 | 0.763 | 0.917 | 0.976 | 1.6 % | 83.1 % |
| COLO829T_ONT_linear | INDEL | somatic | 1,276 | 0.129 | 0.208 | 0.288 | 0.420 | 0.765 | 21.5 % | 27.5 % |
| COLO829T_ONT_linear | INDEL | germline | 688,095 | 0.116 | 0.202 | 0.329 | 0.583 | 0.903 | 24.0 % | 41.2 % |
| COLO829T_Illumina_P75 | SNV | somatic | 33,098 | 0.152 | 0.290 | 0.397 | 0.505 | 0.712 | 8.8 % | 49.3 % |
| COLO829T_Illumina_P75 | SNV | germline | 349,262 | 0.136 | 0.375 | 0.504 | 0.661 | 0.992 | 12.8 % | 71.6 % |
| COLO829T_Illumina_P75 | INDEL | somatic | 957 | 0.125 | 0.202 | 0.290 | 0.387 | 0.582 | 24.1 % | 23.2 % |
| COLO829T_Illumina_P75 | INDEL | germline | 74,579 | 0.098 | 0.181 | 0.327 | 0.473 | 0.769 | 28.4 % | 37.1 % |
| COLO829T_Illumina_P50 | SNV | somatic | 32,073 | 0.121 | 0.199 | 0.277 | 0.343 | 0.455 | 25.0 % | 12.2 % |
| COLO829T_Illumina_P50 | SNV | germline | 352,514 | 0.232 | 0.402 | 0.500 | 0.614 | 0.992 | 3.2 % | 75.5 % |
| COLO829T_Illumina_P50 | INDEL | somatic | 896 | 0.097 | 0.141 | 0.200 | 0.265 | 0.373 | 49.7 % | 2.9 % |
| COLO829T_Illumina_P50 | INDEL | germline | 77,513 | 0.103 | 0.198 | 0.324 | 0.452 | 0.697 | 25.1 % | 35.1 % |
| COLO829T_Illumina_P25 | SNV | somatic | 26,963 | 0.088 | 0.119 | 0.151 | 0.185 | 0.240 | 82.9 % | 0.1 % |
| COLO829T_Illumina_P25 | SNV | germline | 352,922 | 0.306 | 0.436 | 0.500 | 0.565 | 0.992 | 2.4 % | 84.8 % |
| COLO829T_Illumina_P25 | INDEL | somatic | 615 | 0.083 | 0.102 | 0.122 | 0.151 | 0.204 | 93.5 % | 0.3 % |
| COLO829T_Illumina_P25 | INDEL | germline | 78,507 | 0.105 | 0.209 | 0.338 | 0.439 | 0.667 | 23.2 % | 34.9 % |

Somatic AF histogram (coarse bins; 0.02 bins in `vaf_hist.tsv`):

| set | kind | n | 0.00-0.10 | 0.10-0.20 | 0.20-0.30 | 0.30-0.40 | 0.40-0.60 | 0.60-0.80 | 0.80-1.00 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HG008T_Illumina | SNV | 7,134 | 0.1 % | 0.5 % | 3.1 % | 7.2 % | 49.3 % | 6.9 % | 32.9 % |
| HG008T_Illumina | INDEL | 3,121 | 0.9 % | 6.6 % | 12.4 % | 21.8 % | 27.8 % | 15.3 % | 15.2 % |
| HG008T_Illumina_d46 | SNV | 7,274 | 0.1 % | 0.5 % | 3.1 % | 7.3 % | 49.6 % | 6.5 % | 33.0 % |
| HG008T_Illumina_d46 | INDEL | 3,383 | 1.2 % | 6.9 % | 13.0 % | 21.0 % | 28.6 % | 14.5 % | 14.7 % |
| HG008T_Illumina_linear | SNV | 7,371 | 0.1 % | 0.8 % | 3.5 % | 7.3 % | 49.0 % | 7.0 % | 32.3 % |
| HG008T_Illumina_linear | INDEL | 5,893 | 1.7 % | 9.4 % | 15.6 % | 22.1 % | 26.1 % | 15.8 % | 9.4 % |
| HG008T_PacBio | SNV | 7,098 | 0.1 % | 0.5 % | 3.2 % | 7.7 % | 48.6 % | 7.5 % | 32.5 % |
| HG008T_PacBio | INDEL | 3,588 | 2.4 % | 8.3 % | 10.0 % | 16.7 % | 36.2 % | 14.5 % | 11.9 % |
| HG008T_PacBio_d46 | SNV | 7,227 | 0.0 % | 0.4 % | 3.2 % | 7.6 % | 48.9 % | 7.1 % | 32.7 % |
| HG008T_PacBio_d46 | INDEL | 3,761 | 1.5 % | 7.0 % | 10.1 % | 17.9 % | 37.7 % | 14.5 % | 11.5 % |
| HG008T_PacBio_linear | SNV | 7,370 | 0.1 % | 0.9 % | 3.7 % | 7.8 % | 48.3 % | 7.5 % | 31.6 % |
| HG008T_PacBio_linear | INDEL | 6,224 | 0.6 % | 4.5 % | 11.1 % | 20.7 % | 39.8 % | 15.2 % | 8.1 % |
| HG008T_ONT | SNV | 7,124 | 0.1 % | 1.1 % | 4.5 % | 11.2 % | 43.2 % | 7.6 % | 32.3 % |
| HG008T_ONT | INDEL | 3,203 | 2.6 % | 20.8 % | 24.1 % | 18.3 % | 20.8 % | 5.4 % | 8.0 % |
| HG008T_ONT_d46 | SNV | 7,259 | 0.1 % | 1.1 % | 4.5 % | 11.5 % | 43.3 % | 7.3 % | 32.2 % |
| HG008T_ONT_d46 | INDEL | 3,249 | 1.7 % | 19.4 % | 26.8 % | 18.1 % | 20.7 % | 5.4 % | 7.9 % |
| HG008T_ONT_linear | SNV | 7,391 | 0.1 % | 1.6 % | 5.0 % | 11.6 % | 42.8 % | 7.7 % | 31.2 % |
| HG008T_ONT_linear | INDEL | 4,242 | 1.0 % | 17.5 % | 28.3 % | 19.4 % | 21.3 % | 6.0 % | 6.6 % |
| COLO829T_Illumina | SNV | 33,782 | 1.5 % | 5.1 % | 6.9 % | 14.0 % | 39.3 % | 16.1 % | 17.0 % |
| COLO829T_Illumina | INDEL | 971 | 0.8 % | 10.4 % | 21.7 % | 22.2 % | 28.1 % | 8.3 % | 8.3 % |
| COLO829T_Illumina_d46 | SNV | 34,110 | 1.6 % | 5.1 % | 6.8 % | 14.0 % | 39.4 % | 16.0 % | 17.2 % |
| COLO829T_Illumina_d46 | INDEL | 1,048 | 0.5 % | 10.6 % | 21.9 % | 22.1 % | 28.5 % | 8.9 % | 7.4 % |
| COLO829T_Illumina_linear | SNV | 34,341 | 1.6 % | 5.2 % | 6.9 % | 13.9 % | 39.7 % | 16.0 % | 16.8 % |
| COLO829T_Illumina_linear | INDEL | 1,555 | 0.5 % | 10.3 % | 22.3 % | 23.9 % | 27.7 % | 8.9 % | 6.4 % |
| COLO829T_fiberseq | SNV | 33,669 | 1.5 % | 5.0 % | 6.4 % | 14.5 % | 40.3 % | 15.4 % | 17.0 % |
| COLO829T_fiberseq | INDEL | 1,078 | 1.9 % | 7.1 % | 11.2 % | 23.4 % | 37.0 % | 12.2 % | 7.1 % |
| COLO829T_fiberseq_d46 | SNV | 33,998 | 1.5 % | 5.0 % | 6.3 % | 14.4 % | 40.3 % | 15.3 % | 17.1 % |
| COLO829T_fiberseq_d46 | INDEL | 1,141 | 1.9 % | 6.3 % | 11.8 % | 24.5 % | 36.2 % | 12.4 % | 6.7 % |
| COLO829T_fiberseq_linear | SNV | 34,100 | 1.5 % | 5.0 % | 6.4 % | 14.2 % | 40.8 % | 15.4 % | 16.7 % |
| COLO829T_fiberseq_linear | INDEL | 1,641 | 0.7 % | 3.0 % | 12.7 % | 25.3 % | 38.3 % | 12.5 % | 7.4 % |
| COLO829T_ONT | SNV | 33,472 | 1.8 % | 5.6 % | 10.2 % | 17.3 % | 40.3 % | 9.7 % | 15.3 % |
| COLO829T_ONT | INDEL | 999 | 2.6 % | 22.8 % | 29.8 % | 18.2 % | 18.5 % | 4.0 % | 4.0 % |
| COLO829T_ONT_d46 | SNV | 33,779 | 1.8 % | 5.6 % | 10.4 % | 17.7 % | 40.1 % | 9.3 % | 15.1 % |
| COLO829T_ONT_d46 | INDEL | 1,018 | 1.8 % | 23.8 % | 30.6 % | 18.5 % | 17.9 % | 3.9 % | 3.6 % |
| COLO829T_ONT_linear | SNV | 33,852 | 1.7 % | 5.5 % | 10.0 % | 17.2 % | 40.7 % | 9.6 % | 15.1 % |
| COLO829T_ONT_linear | INDEL | 1,276 | 1.6 % | 19.9 % | 31.5 % | 19.5 % | 17.7 % | 5.6 % | 4.2 % |
| COLO829T_Illumina_P75 | SNV | 33,098 | 1.6 % | 7.2 % | 18.0 % | 23.9 % | 33.4 % | 15.4 % | 0.4 % |
| COLO829T_Illumina_P75 | INDEL | 957 | 1.9 % | 22.3 % | 29.3 % | 23.4 % | 19.5 % | 3.7 % | 0.0 % |
| COLO829T_Illumina_P50 | SNV | 32,073 | 2.3 % | 22.8 % | 34.3 % | 28.6 % | 12.0 % | 0.2 % | 0.0 % |
| COLO829T_Illumina_P50 | INDEL | 896 | 6.2 % | 43.4 % | 34.2 % | 13.3 % | 2.8 % | 0.1 % | 0.0 % |
| COLO829T_Illumina_P25 | SNV | 26,963 | 12.0 % | 70.8 % | 16.4 % | 0.6 % | 0.1 % | 0.0 % | 0.0 % |
| COLO829T_Illumina_P25 | INDEL | 615 | 22.4 % | 71.1 % | 5.9 % | 0.3 % | 0.2 % | 0.2 % | 0.0 % |

- HG008T somatic SNV positives (all three arms): median AF 0.53-0.54, 82-89 % at AF >= 0.40 and 31-33 % at >= 0.80
  (high purity and large LOH regions); 0.4-1.8 % below AF 0.20. COLO829T P100 (all arms): median 0.46-0.51, 6.5-7.4 %
  below 0.20. The reference arm hardly changes the somatic AF distribution of a set.
- On the pangenome arms the germline SNV tensors have almost the same AF distribution as the somatic ones (median
  0.48-0.60, 69-89 % >= 0.40): the model cannot use AF alone to tell them apart. On GRCh38-only the germline SNV median
  is 0.76-0.99 (83-94 % >= 0.40): there every germline variant, homozygous ones included, becomes a candidate, while on
  a pangenome the reads of common germline alleles follow graph paths [interpretation of the 8-11x more germline
  tensors, unverified]. Non tensors are mostly low AF in every set (SNV median 0.10-0.24; 43-93 % below 0.20).
- So the P100 training positives contain very few low-VAF somatic variants. The COLO829T dilutions are the only source of
  them: Illumina P75 / P50 / P25 somatic SNV median AF 0.40 / 0.28 / 0.15, below 0.20 for 8.8 / 25.0 / 82.9 % (INDEL
  median 0.29 / 0.20 / 0.12), while their germline SNV median stays at 0.50 (the diluting reads are COLO829BL's, which
  carry the same germline variants). Only the three joint models train on these sets. The fiberseq / ONT dilutions and
  HCC1395T: pending (held jobs).
- Long-read INDEL germline tensors on the pangenome arms are lower in AF (PacBio / fiberseq median 0.22-0.23, ONT
  0.17-0.20) than the somatic ones (0.27-0.45): many germline INDEL tensors are probably the A1 of a homopolymer site
  whose reads split between lengths [interpretation, unverified]. On GRCh38-only the long-read germline INDEL median is
  0.33-0.57.

### 3b. Mutation spectrum of the SNV positives (R1.7 tumor diversity; `tables.md` T7, `somatic_spectrum.tsv`)

| set | SNV alleles | C>A | C>G | C>T | T>A | T>C | T>G | UV-type C>T | APOBEC-type TCW | TCW C>G | CpG C>T |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HG008T_Illumina | 7,134 | 18.6 % | 13.0 % | 40.6 % | 12.4 % | 10.0 % | 5.5 % | 29.2 % | 20.7 % | 8.4 % | 9.3 % |
| HG008T_Illumina_d46 | 7,274 | 18.7 % | 12.8 % | 40.6 % | 12.3 % | 10.1 % | 5.5 % | 29.2 % | 20.3 % | 8.3 % | 9.5 % |
| HG008T_Illumina_linear | 7,371 | 18.9 % | 12.2 % | 39.8 % | 13.1 % | 10.1 % | 6.0 % | 28.4 % | 19.7 % | 7.9 % | 9.2 % |
| HG008T_PacBio | 7,098 | 18.6 % | 13.1 % | 40.7 % | 12.3 % | 9.9 % | 5.3 % | 29.3 % | 20.8 % | 8.5 % | 9.4 % |
| HG008T_PacBio_d46 | 7,227 | 18.8 % | 13.0 % | 40.8 % | 12.3 % | 9.9 % | 5.3 % | 29.3 % | 20.4 % | 8.4 % | 9.6 % |
| HG008T_PacBio_linear | 7,370 | 18.8 % | 12.2 % | 39.7 % | 13.2 % | 10.0 % | 5.9 % | 28.3 % | 19.6 % | 7.9 % | 9.2 % |
| HG008T_ONT | 7,124 | 18.6 % | 13.0 % | 40.8 % | 12.3 % | 10.0 % | 5.3 % | 29.4 % | 20.9 % | 8.5 % | 9.4 % |
| HG008T_ONT_d46 | 7,259 | 18.8 % | 13.0 % | 40.8 % | 12.2 % | 10.0 % | 5.4 % | 29.3 % | 20.4 % | 8.3 % | 9.6 % |
| HG008T_ONT_linear | 7,391 | 18.9 % | 12.3 % | 39.8 % | 13.1 % | 10.0 % | 5.9 % | 28.4 % | 19.7 % | 7.9 % | 9.2 % |
| COLO829T_Illumina | 33,782 | 11.4 % | 2.4 % | 71.0 % | 4.6 % | 6.1 % | 4.6 % | 64.6 % | 20.9 % | 0.8 % | 6.7 % |
| COLO829T_Illumina_d46 | 34,110 | 11.4 % | 2.3 % | 71.0 % | 4.6 % | 6.2 % | 4.6 % | 64.6 % | 20.8 % | 0.8 % | 6.8 % |
| COLO829T_Illumina_linear | 34,341 | 11.4 % | 2.3 % | 70.9 % | 4.6 % | 6.1 % | 4.6 % | 64.5 % | 20.8 % | 0.8 % | 6.8 % |
| COLO829T_fiberseq | 33,669 | 11.5 % | 2.3 % | 71.0 % | 4.6 % | 6.1 % | 4.6 % | 64.6 % | 21.0 % | 0.8 % | 6.7 % |
| COLO829T_fiberseq_d46 | 33,998 | 11.4 % | 2.3 % | 71.1 % | 4.6 % | 6.1 % | 4.5 % | 64.6 % | 20.9 % | 0.8 % | 6.8 % |
| COLO829T_fiberseq_linear | 34,100 | 11.5 % | 2.3 % | 71.0 % | 4.6 % | 6.1 % | 4.6 % | 64.6 % | 20.8 % | 0.8 % | 6.8 % |
| COLO829T_ONT | 33,472 | 11.3 % | 2.2 % | 71.5 % | 4.4 % | 6.0 % | 4.5 % | 65.1 % | 21.1 % | 0.8 % | 6.7 % |
| COLO829T_ONT_d46 | 33,779 | 11.3 % | 2.2 % | 71.6 % | 4.4 % | 6.1 % | 4.4 % | 65.2 % | 21.0 % | 0.8 % | 6.8 % |
| COLO829T_ONT_linear | 33,852 | 11.3 % | 2.2 % | 71.5 % | 4.4 % | 6.0 % | 4.5 % | 65.1 % | 20.9 % | 0.8 % | 6.9 % |
| COLO829T_Illumina_P75 | 33,098 | 11.2 % | 2.2 % | 71.8 % | 4.4 % | 6.0 % | 4.4 % | 65.4 % | 21.1 % | 0.8 % | 6.8 % |
| COLO829T_Illumina_P50 | 32,073 | 11.0 % | 1.8 % | 73.0 % | 4.2 % | 5.8 % | 4.1 % | 66.6 % | 21.3 % | 0.6 % | 6.9 % |
| COLO829T_Illumina_P25 | 26,963 | 10.8 % | 1.4 % | 74.8 % | 4.0 % | 5.5 % | 3.6 % | 68.3 % | 21.4 % | 0.4 % | 7.1 % |

- COLO829T positives are UV-dominated (71 % C>T, 65 % C>T at a dipyrimidine, 0.8 % C>G at TCW). HG008T positives are
  mixed (40-41 % C>T, 28-29 % at a dipyrimidine, 19 % C>A, 7.9-8.5 % C>G at TCW, 9.2-9.6 % C>T at CpG). The reference
  arm changes these fractions by at most about one point. These are context counts, not a signature fit (E18 does the
  SBS96 decomposition).
- The alleles that survive dilution are slightly more UV-like (C>T 71.0 % at P100, 74.8 % at P25; UV-type 64.6 % ->
  68.3 %), i.e. the low-VAF COLO829T SNVs lost at low purity carry less of the UV pattern [interpretation, unverified].
- HCC1395T and the CASTLE lines (spectra of the added lines): pending (held jobs).


## 4. Label rules in plain English (R1.8)

The rules are `indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py` (module docstring lines 1-58,
`classify` lines 589-688), explained with real examples in `analysis/truth_labels.md` (rule table lines 129-149, notes
lines 151-160). The first rule that applies decides.

1. **Truth alleles become graph keys, all placements of an indel included.** Every ALT is split off and trimmed, MNVs
   become SNVs; an indel in a repeat is expanded to every equivalent placement (up to `MAX_SHIFTS` = 5,000, line 86);
   each placement is converted along the graph's GRCh38 path into the builder's candidate key `node:start:KIND:REF>ALT`
   on the node's forward strand (reverse nodes flipped and reverse-complemented; an insertion at a node junction gets a
   key on both nodes; a deletion over consecutive nodes gets the builder's multi-node key). Docstring lines 3-15;
   `truth_labels.md` section 2 (lines 60-102). Because every placement is enumerated, it does not matter at which end
   of a repeat the builder (left-aligned on the node's forward strand) or the VCF (left-aligned on GRCh38) writes the
   indel: no separate normalisation step is needed on either side.
2. **Somatic (1)**: the tensor's representative allele A1 (the allele with most ALT reads of its site) equals a PASS
   somatic truth key: `classify` lines 622-624 (`representative_allele_is_somatic_truth`; `truth_labels.md` R2, line
   134). This holds inside or outside the BED; when A1 is both a somatic and a germline truth allele, somatic wins.
   Partial matches (allele overlap or read-haplotype overlap > 0.45, `MIN_OVERLAP` line 79; R4 lines 647-651, R8 lines
   673-676) are labelled 1 in the files but **left out of training and validation** by the recipes
   (`--ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth`: `ris/recipes.sh` lines 6 and 9;
   `machine_learning/pansoma_net_v2/train.py` lines 337-346). So every training positive is an exact allele match.
3. **Germline (2)**: A1 equals a PASS germline truth key, or one whose only FILTER is dipcall's GAP1 / GAP2
   (`classify` lines 626-629, `GAP_FILTERS` line 76, `counts_as_truth` lines 335-338); INDELs that overlap a germline
   allele by > 0.45 in the same site or at the same place also get 2 (R4, lines 647-651; R10, lines 681-684). The germline truth is
   the matched normal's: the dipcall calls of the HG008-N assembly (HG008T) and of the COLO829BL donor-specific
   assembly (COLO829T); DeepVariant calls on the matched normal with QUAL >= 10 for HCC1395T and the CASTLE lines (the
   `truth.germline.vcf` of each `labels.manifest.json`). The normal is used only to label the training data and to
   score; calling is tumor-only.
4. **Non (0)**: a tensor on a GRCh38 node inside the confident region (somatic BED ∩ germline BED, line 802) that is
   no truth allele: nothing within 10 bp (R12 `confident_no_truth_allele`, line 688), a truth allele within 10 bp
   that A1 does not match (R11 `near_truth_allele_mismatch`, line 687: errors and artefacts next to real variants),
   or a site whose truth allele is another of its alleles that A1 does not overlap enough (R4, line 652). These are
   the sequencing and alignment artefacts; the model also sees the real germline variants as class 2.
5. **Ignore (-1)**, never trained on:
   - outside the confident region (R6, lines 667-668);
   - no GRCh38 position: not on a unique GRCh38 node and no reference neighbours within 200 node IDs / 1,024 bp
     (R5, lines 659-662);
   - an off-reference (branch-node) tensor without a partial somatic match (R9, lines 679-680), because the germline
     truth has keys only on GRCh38 nodes; validation and the chr1 test count these as non
     (`machine_learning/pansoma_net_v2/data.py` `EVAL_AS_NON`, line 51, applied at lines 197-199, which also leave out every tensor outside the confident region, labels 1 and 2 included);
   - the partial somatic matches of point 2 (training only; they are 1 in the label files);
   - rule (a), only for HCC1395T and the CASTLE lines (`ignore_repeat_indels: true` in their manifests): an INDEL that
     would be 0 in a homopolymer >= 7 bp or an STR (R13, `non()` lines 605-611, `IGNORED_REPEAT_CLASSES` line 87;
     `truth_labels.md` line 160), because SEQC2 / CASTLE truth sets omit most such INDELs;
   - a filtered somatic truth allele (R2', R4'): 0 tensors in every set (`truth_labels.md` section 7).
   "Determinable" (R1.8) therefore means: inside the confident region, on a GRCh38 node (a unique GRCh38 position), and
   either an exact truth allele or none.
6. **No AF floor at labelling; AF 0.08 at the build.** The label runs set no `snv_min_af` / `indel_min_af` (null in
   every `labels.manifest.json` and in `run/config.json` `postprocess.labels`; R1 lines 616-619 never fire;
   `truth_labels.md` lines 42-43; `pansoma_v2_tensors/README.txt` line 11). The build writes a tensor only for an
   allele with A1 AF >= 0.08 (`--snv-min-af 0.08 --indel-min-af 0.08`), at least 3 ALT-supporting reads
   (`--min-variants 3`), alignments with MAPQ > 10 (`--min-mapq 10` is exclusive) and ALT observations with base
   quality >= 10 (`--min-allele-bq 10`; `indexed_gam_pipeline_v4/run.py` lines 77, 80-81) (the build command in `<SET>/run/logs/task_*.resources.txt`; the same values in each
   record's `parameters` in `variant_summary.ndjson`). Every labelled tensor is used each epoch (`--non-fraction 1.0`)
   with class weights `pow:0.75` (SNV) or `sqrt` (INDEL) (`ris/recipes.sh`). VAF and purity are not resampled: each
   line's natural VAF distribution, plus the COLO829T in-silico dilutions in the joint models.

## 5. Comparator training scale (R1.7)

| caller | training data | source |
|---|---|---|
| Pansoma (these models) | one tumor line per model, chr2-22 (HG008T: 8,690 SNV + 8,496 INDEL somatic truth alleles on chr1-22; COLO829T: 42,035 + 1,912), Illumina / PacBio-fiberseq / ONT trained separately; the COLO829T joint models add three in-silico purities of the same line | this folder; truth counts from the `truth.somatic.stats` of the labels manifests |
| DeepSomatic | the five CASTLE tumor / normal cell-line pairs HCC1395, HCC1937, HCC1954, H1437, H2009, on Illumina, PacBio HiFi and ONT, with augmentation; the truth sets of HCC1954 / H1437 / H2009 are its training truth | the detailed response (lines 632, 809-812; the reviewer's statement and our reading); chromosomes held out, number of training examples and the augmentation [unverified: E8b, literature] |
| ClairS-TO | [unverified: not in the detailed response; another agent checks the literature] | - |

For scale on our side: the four added lines carry 41,072 (HCC1395T, SEQC2 v1.2.1), 25,014 (HCC1954T), 84,744 (H1437T)
and 161,091 (H2009T) somatic truth alleles on chr1-22 (311,921 together; 323,516 with the other contigs, which is
the "about 323 k" of the detailed response), against 43,947 (46,064) for COLO829T and 17,186 for HG008T (labels
manifests, `truth.somatic.stats`).

## 6. E13: GPU jobs for the learning curve and the model size (not run; GPUs)

`train.py` has no option to subsample node blocks, but `--exclude-candidates FILE` takes listed tensors out before
the model ("as if never built": -1 in training and validation; predict defaults to the checkpoint's file and graph_vcf
writes no record for them). `make_e13_subsets.py` writes, per kept fraction f and replicate r, the candidate ids of all
labelled tensors in the *dropped training blocks*; validation blocks and chr1 are untouched, so the validation set, its
PASS threshold rule and the chr1 test are those of the full run, and nothing on chr1 is excluded (the script checks that
no excluded id occurs outside the training blocks). Subsets are nested within a replicate (10 % ⊂ 25 % ⊂ 50 % ⊂ 75 %).
`--epoch-samples` is not a substitute: it takes a fresh random subset every epoch, so the model still sees every
training tensor; `--train-chroms` would confound data size with chromosome. Changing `--seed` would also move the
validation blocks (`data.block_split` uses it), so replicates use a different subset draw instead.

Exclude files (16, 70 MB): `/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/e13/<SET>.SNV.keep<10|25|50|75>.rep<1|2>.exclude.txt`
for `HG008T_PacBio_d46` (the d46 arm, as asked) and `HG008T_PacBio` (d9, in case the main figure stays on d9). What
each subset trains on (`e13_subsets_<SET>_SNV.tsv`; the 100 % row equals the `train.log` of the full run):

| HG008T_PacBio_d46 SNV | kept blocks (of 5,662) | train non | train somatic (= distinct truth alleles) | train germline |
|---|---:|---:|---:|---:|
| 100 % (`HG008_PacBio_SNV_d46_7fc7396`) | 5,662 | 71,573 | 7,227 | 270,175 |
| 75 % rep 1 / rep 2 | 4,257 / 4,245 | 57,025 / 53,933 | 5,289 / 5,364 | 202,590 / 195,385 |
| 50 % rep 1 / rep 2 | 2,810 / 2,859 | 42,021 / 26,972 | 3,512 / 3,655 | 136,122 / 133,033 |
| 25 % rep 1 / rep 2 | 1,431 / 1,455 | 32,289 / 13,584 | 1,773 / 1,867 | 66,634 / 64,920 |
| 10 % rep 1 / rep 2 | 556 / 567 | 8,410 / 5,352 | 703 / 779 | 26,642 / 23,645 |

The somatic and germline counts scale with the kept fraction; the non count does not (25 %: 32,289 against 13,584):
the PacBio SNV non tensors are concentrated in few blocks. Two replicates show how much that matters.

**Jobs** (RIS, `ris/` scripts and the SNV recipe of `ris/recipes.sh`, set up exactly like the d46 run they extend:
4 H100 under DDP, `--batch-size 1024` per GPU). First copy the exclude files to RIS (needs the user: scp), then:

```bash
cd /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma
E13=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs/e13_subsets   # copy of .../E21_training_composition/e13
# learning curve: 4 fractions x 2 replicates = 8 runs (the 100 % point is HG008_PacBio_SNV_d46_7fc7396)
for rep in 1 2; do for f in 10 25 50 75; do
  KIND=SNV GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 ris/snv_pipeline.sh HG008_PacBio_SNV_d46_e13_keep${f}_rep${rep} \
      HG008T_PacBio_d46 --exclude-candidates $E13/HG008T_PacBio_d46.SNV.keep${f}.rep${rep}.exclude.txt
done; done
# model size on 100 % of the blocks: 8.7 M and 199.6 M (the 50.4 M point is HG008_PacBio_SNV_d46_7fc7396); the extra
# options come after the recipe's in run.sh, and argparse keeps the last --depths / --dims
KIND=SNV GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 ris/snv_pipeline.sh HG008_PacBio_SNV_d46_e13_8p7M HG008T_PacBio_d46 \
    --depths 2 2 6 2 --dims 64 128 256 512
KIND=SNV GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 ris/snv_pipeline.sh HG008_PacBio_SNV_d46_e13_199p6M HG008T_PacBio_d46 \
    --dims 192 384 768 1536 --micro-batches 2
```

- `--micro-batches 2` for 199.6 M: batch 1024 of that model "fills a 94 GB H100 NVL" (`pansoma_net_v2/README.md`, line
  112); the RIS H100s have 80 GB. Same update as the whole batch.
- Cost: the 100 % d46 run trained 21.5 min on 4 GPUs (`metrics.jsonl`, 20 epochs); the subsets are proportionally
  shorter per epoch, the 199.6 M model about 4x slower per epoch [unverified for SNV; the README measured 214 s vs 55 s
  per epoch for 199.6 M vs 8.7 M on HG008T Illumina INDEL]. All 10 runs: roughly 3-5 GPU-hours on 4 H100
  [estimate]. Each `snv_pipeline.sh` call also submits `vcf.sh` (graph VCF, GRCh38 VCF, PoN, rtg vcfeval on chr1).
- The d9 alternative: the same loops with `HG008T_PacBio`, run name `HG008_PacBio_SNV_e13_*`, files
  `HG008T_PacBio.SNV.keep*.rep*.exclude.txt`, and `GPUS=1 SLOTS=7 MEM=128` like `HG008_PacBio_SNV_7fc7396` (1 GPU,
  35.8 min training).
- Fixed 20 epochs means fewer optimizer steps on smaller subsets; the full model's best epoch is 7, so a 10 % model
  gets about 2 full-data epochs' worth of steps. If that confounds the curve, a matched-steps variant adds
  `--epochs 200 / 80 / 40 / 27` for 10 / 25 / 50 / 75 % (validation every epoch is cheap) [optional].
- Readout per run: chr1 PoN-after PASS F1 and best F1 in the BED (`vcf_chr1/eval_bed/report.json`), per-tensor chr1
  somatic AP (`test_chr1/*.metrics.json`), validation somatic AP, best epoch; the same `collect.py` pattern as
  `analysis/benchmark_7fc7396_20261012/hprc_v21_d46_7fc7396_20261012`.

## 7. Model size (loss curves) and the multi-tumor training sets (R1.7; 2026-10-20, `capacity_multitumor.md`)

`capacity_multitumor.py` (login node, reads only `train.log` / `metrics.jsonl` / `report.json`) writes
`capacity_multitumor.md`, `capacity_loss.tsv` (train and validation loss per epoch) and `multitumor_training.tsv`.

- **A. Model size** (E13, HG008T PacBio d46 SNV, lr 1e-4, 20 epochs; validation loss = class-weighted cross-entropy):
  199.6 M: minimum 0.399 at epoch 3, 5.76x at epoch 20, training loss 0.038, selected epoch 6, PASS F1 0.603; 50.4 M
  (three seeds): minimum at epochs 9-10, 2.56-3.43x at epoch 20, training loss 0.112-0.130, selected 7 / 9 / 9, PASS F1
  0.604-0.611; 8.7 M: minimum at epoch 13, 1.04x at epoch 20, training loss 0.252, selected 16, PASS F1 0.616. The
  INDEL pair of 2026-10-06 (HG008T Illumina d9, lr 4e-4, 12 epochs, sqrt weights; PoN rule of that day for both):
  199.6 M training loss 0.0004, validation loss 5.33x its minimum at epoch 12, PASS F1 0.1865; 8.7 M 0.225, 1.91x,
  0.2426. Larger networks overfit earlier and further; with early stopping the SNV F1 is the same within 0.013.
- **B. Multi-tumor models** (12 runs; the per-set training counts add up to each run's total, checked by the script):
  SNV somatic training tensors 237,577-240,469 (multicancer4) and 278,615-281,853 (multicancer6), INDEL 10,390-10,584
  and 14,809-15,486; selected epochs SNV 10-19, INDEL 6-9 (single-set models: 4-12, 2-7). Validation counts are per run
  (`train.log` does not split them by set).

## Interpretation for the response

- **R1.7, example counts** (fills "Supplementary Table [待补]"): use `tables.md` T1 (per model) and T3 (per set). The
  detailed response's figures for HG008T SNV sets on chr1-22 (8.8-10.3 k somatic, 0.31-0.37 M germline, 0.07-1.9 M non)
  are the label-manifest totals over all tensors; what the models actually train on is smaller: 7.1-7.4 k somatic,
  0.26-0.33 M germline (pangenomes; 2.7-2.8 M on GRCh38-only) and 0.07-1.8 M non per model (chr2-22, minus 5 % for
  validation, minus chr1, minus the -1 tensors). The response should quote the training counts.
- **R1.7, capacity**: honest wording: the 50.4 M network has ~7,000 parameters per HG008T SNV positive and can fit its
  training set (training loss -> 0.002-0.14), so it is regularised by early stopping on held-out node blocks
  (checkpoints at epochs 4-12); the 199.6 M model gave no better test result and the 8.7 M INDEL model is within noise of
  it (`pansoma_net_v2/README.md`). Whether accuracy still rises with data is the E13 learning curve (section 6).
- **R1.7, diversity**: the two training lines differ in spectrum (UV-dominated COLO829T vs mixed HG008T, T7) but each
  model is trained on one line. The cross-line evidence the reviewer asks for is E12 (leave-one-line-out); T7 gives the
  spectra to stratify it, and the HCC1395T / CASTLE rows (held jobs) the spectra of the added lines.
- **R1.7, joint models**: the purity-joint models have 3.5-3.8x the positives of a single-set model but the same
  alleles (Illumina: 33,797 against 33,782 SNV alleles); they add VAF diversity, not mutation diversity. The response
  should not count their tensors as extra training examples of different variants.
- **R1.8, positives**: every training positive is an exact allele match of the representative allele (partial matches
  are left out), every placement of an indel in a repeat is tried, and off-reference tensors are never given a coordinate
  through neighbouring nodes (they are -1 unless their reads match a somatic truth haplotype, and then they are partial
  and left out of training). The "is the positive set distorted" table: T3 (on the pangenome arms SNV 89.5-95.9 % of the
  training-region truth reaches training, INDEL 43-67 %, the INDEL loss being alleles without any candidate, R1.4; on
  GRCh38-only SNV 90.1-98.7 %, INDEL 57-97 %, so the INDEL loss comes with the pangenome, not with the labelling).
- **R1.8, negatives / VAF / purity**: germline is an explicit class (0.26-0.38 M SNV tensors per pangenome model);
  artefacts = everything in the confident region on a GRCh38 node without a truth allele (non); VAF follows each line's
  natural distribution with no AF floor at labelling (AF 0.08 at the build), and positives below AF 0.20 are rare in the
  single-set models (0.4-1.8 % HG008T, 6.5-7.4 % COLO829T): purity is only varied through the COLO829T dilutions, which
  bring the share below AF 0.20 to 8.8 / 25.0 / 82.9 % (Illumina P75 / P50 / P25) in the joint models' inputs.
  This fills the "[待补: VAF histogram of training positives per dataset]" of R1.8 with T5 / T6.

## Caveats

- Counts are what the recipe's split gives (seed 0, 5 % of the 20,000-node blocks): another seed moves the validation
  blocks. All 36 local `args.json` use seed 0, `--val-fraction 0.05`, `--val-block-nodes 20000`, the same
  `--ignore-reasons`, no `--exclude-candidates`, `--non-fraction 1.0`; the six RIS-only logs state the same split.
- The six COLO829T d9 single-set models are older (SNV 2026-10-03, 1 GPU; INDEL 2026-10-06, lr 4e-4); their split and
  labels are the same (replay matches), their validation AP and threshold come from the benchmark `summary.json`.
- "Distinct truth alleles behind training tensors" counts the representative matches only (partial matches are left out
  of training), by truth id; the unions across sets match alleles by VCF chrom / pos / REF / ALT.
- The spectrum uses the GRCh38 bases around the VCF position (no signature fitting); C>T at TCA / TCT counts both as
  UV-type and as APOBEC-type.
- CASTLE rows are a snapshot of the PacBio sets' labels of 2026-10-12 (`17b540137dc3`, germline QUAL >= 10); their ONT /
  Illumina sets were being relabelled or built and are not included. HCC1395T / CASTLE use rule (a), so their INDEL
  negatives are not comparable with HG008T / COLO829T.
- Comparator training data (section 5) is [unverified] beyond what the detailed response states; the comparators are not
  re-run here (user decision: their ClairS-TO / DeepSomatic results exist).
- Until the held caches exist, the `fiberseq joint ...` and `ONT joint ...` rows of `somatic_truth_unions.tsv` /
  `tables.md` T4 hold only the P100 set (the label says joint), and `runs_replay_check.tsv` leaves those four runs
  unchecked (`match` empty; its `replay_*` columns there are the P100 set alone, not a mismatch).
- The tables of today were made by running `collect_runs.py`, `summarize.py` and `make_tables.py` by hand on the login
  node, not by job 387193. When 387193 runs it rewrites every TSV and `tables.md` (and reruns `make_e13_subsets.py`);
  the 21 sets read today should come out unchanged as long as their labels do not change (the caches are keyed on the
  labels meta).
- HCC1395T / CASTLE numbers will be those of the labels present when the held index tasks run, not necessarily the
  2026-10-12 snapshot named above.
- AF is the tensor's A1 fraction over the reads at the site, after the build's read filters (MAPQ > 10, ALT base quality >= 10),
  not the truth VAF.

## Status

**PENDING** (mostly done; 2026-10-13 23:01).

DONE:

- Section 1: all 42 runs (36 single-set, 6 purity-joint) from `train.log` / `metrics.jsonl`; replay exact for 38 of
  42, 0 differ.
- Sections 2, 3 and 3b for the 21 sets with an index cache: all 18 reference-arm sets (HG008T / COLO829T x Illumina /
  PacBio-fiberseq / ONT x d9 / d46 / linear) and the COLO829T Illumina P75 / P50 / P25 dilutions, so also the exact
  allele union of the Illumina joint models.
- Sections 4 (label rules with code lines) and 5 (comparator scale; no comparator re-run, per the user).
- Section 6: the 16 E13 exclude files and their subset tables (`e13_subsets_*_SNV.tsv`), and the RIS commands.

PENDING, with the reason:

- Index caches of 12 sets (Slurm array 387137 tasks 22-33: COLO829T fiberseq / ONT P75 / P50 / P25, HCC1395T Illumina /
  PacBio / ONT d46, HCC1954T / H1437T / H2009T PacBio d46) and the collect job 387193: **held** (`JobHeldUser`) by a
  peer session after the user's complaint about the lab-cluster footprint; they stay held until the user agrees a
  cluster cap (not released here). They give: the replay check of the 4 fiberseq / ONT joint runs, the fiberseq / ONT
  joint allele unions and dilution VAF, and the HCC1395T / CASTLE truth coverage, VAF and spectrum (R1.7 diversity
  rows). Once released, 387193 runs `collect.sh` by itself; by hand afterwards:
  `bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E21_training_composition/collect.sh`.
- E13 (learning curve and model size): not submitted. Needs the exclude files copied to RIS by the user (scp; RIS data
  sync) and RIS GPU time (10 runs, section 6); which arm (d46 or d9) is the user's decision.
- ClairS-TO training data (section 5): [unverified], left to E8b (literature).

## Files

| file | what |
|---|---|
| `build_index.py` | index cache of one set (SNV + INDEL in two processes) |
| `collect_runs.py` | per-run table from `train.log` / `metrics.jsonl` -> `runs_training.tsv`, `runs_training.md`, `runs_training.labels_created.json` |
| `summarize.py` | replayed split, counts, truth alleles, VAF, spectrum -> `set_split_counts.tsv`, `runs_replay_check.tsv`, `somatic_truth_alleles.tsv`, `somatic_truth_unions.tsv`, `vaf_quantiles.tsv`, `vaf_hist.tsv`, `somatic_spectrum.tsv` |
| `make_e13_subsets.py` | E13 exclude files -> `e13_subsets_<SET>_SNV.tsv` here, the files in `/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/e13/` |
| `make_tables.py` | `tables.md` (T1-T7) |
| `capacity_multitumor.py` | section 7: `capacity_multitumor.md`, `capacity_loss.tsv`, `multitumor_training.tsv` |
| `collect.sh` | all of the above after the index builds |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E21_training_composition/` | `build_index.sbatch`, `sets.txt`, `collect.sbatch`, logs |
| `/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/index_cache/` | index caches (large, not in git) |
