# Pansoma on RIS compute1 (WashU)

The other READMEs were written on a Slurm cluster (`/scratch/jshen`, GPU node tequila,
`/wanglab/jshen/anaconda3/bin/python`). RIS compute1 runs LSF, and every job runs in a Docker container. This
directory holds the job scripts for RIS and records how the model steps (`machine_learning/pansoma_net_v2`) run
there. The tensor pipeline itself (`indexed_gam_pipeline_v4`) has not been run on RIS, apart from relabelling sets
(`relabel.sh`).

RIS documentation: <https://washu.atlassian.net/wiki/spaces/RUD/overview> (compute1 is retiring in July 2027).

## Storage

| what | path |
|---|---|
| repository | `/storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma` |
| tensor sets | `/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors/<SET>/tensors` (`HG008T_Illumina`, `HG008T_ONT`, `HG008T_PacBio`, `COLO829T_Illumina`, `COLO829T_ONT`, `COLO829T_fiberseq`) |
| graph index | `Dataset/pansoma_v2_tensors/graph_index` (`hprc-v1.1-mc-grch38.d9.graph_index.sqlite`, `.grch38_path`; links to the HPRC v2.1 d46 and GRCh38-only indexes copied to `Dataset/AF-Filtered_VG_Indexes_HPRC_V2/graph_index` and `Dataset/Linear_Reference_VG_Indexes/graph_index`, checked with `tools/check_graph_index.py`) |
| truth of the labels | `Dataset/HG008_GIAB/draft_v02_benchmark`, `Dataset/HG008_GIAB/dipcall_HG008N_GRCh38`, `Dataset/pansoma_v2_tensors/COLO829T_truth`, `Dataset/COLO829BL_DSA/dipcall_hg38` |
| PoNs | `Dataset/panel_of_normal_VCFs` (gnomAD, dbSNP, 1000G, CoLoRSdb, with `.tbi`) |
| GRCh38 | `Dataset/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa` (`.fai`) and its rtg SDF `….sdf` |
| runs, environment, logs | `/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs` (`$PANSOMA_RUNS`) |

`$PANSOMA_RUNS` holds `env/` (Python user base, rtg, Triton and inductor caches), `index_cache/` (the index
and zstd chunk caches, shared by all runs), `logs/<JOBNAME>.log` and one directory per run. `~/jiawei.s` is a link
to `/storage2/fs1/epigenome/Active/jiawei.s`. There is no `/scratch`, and `/scratch1/fs1` is not readable for us.

**Truth paths.** `submit.sh` mounts (`MOUNTS`): `Dataset/HG008_GIAB/{draft_v02_benchmark,dipcall_HG008N_GRCh38}` →
`/scratch/jshen/data/HG008_GIAB/...`, `Dataset/pansoma_v2_tensors/COLO829T_truth` →
`/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth`, `Dataset/COLO829BL_DSA/dipcall_hg38` →
`/scratch/qfu/COLO829BL_DSA/dipcall_hg38` (all sha256-checked against the manifests on 2026-10-01).

**How the paths are read.** Each `labels.manifest.json` names its truth VCFs and BEDs by their old-cluster paths
(`/scratch/jshen/data/HG008_GIAB/...`), and `data.KindIndex` reads the BEDs from there (the confident region).
The copies in `Dataset/HG008_GIAB` have the manifests' sha256. `submit.sh` bind-mounts the two directories at the
old paths in the container, so the manifests stay unchanged. Without them the index stops with
`FileNotFoundError: /scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/...`. GIAB's FTP has only the v0.3 draft of
the HG008-T somatic benchmark, not the v0.2 the labels used.

## GPUs and queues

Measured 2026-09-29 (`bhosts -gpu -w`, `lshosts -gpu -w`, `bqueues -l`):

| queue | hosts | GPUs | notes |
|---|---|---|---|
| `epigenome` | compute1-exec-419 (48 slots, 1.9 TB) | 8× H100 80GB HBM3 | our condo; `-G compute-epigenome-condo`; on 2026-09-29 all 8 were held by long jobs of another user, free on 2026-10-01 |
| `epigenome-interactive` | the same host | | `bsub -Is`, 24 h limit |
| `general` | A40 48GB (exec-357…360), A100 40GB (281, 282), A100 80GB (393…401), V100 32GB (203…216) | | `-G compute-epigenome`; the A100 80GB hosts were `closed_Adm`; 28 d limit |
| `general-interactive` | V100 only (197, 199, 205) | | 24 h limit |

- The GPU driver is 550.54.15 (CUDA 12.4) on every GPU host.
- `gpuhost` is an exclusive resource (`!gpuhost` in `lshosts`): a GPU job needs `-R gpuhost`, or it pends with
  "There are no suitable hosts for the job".
- `-gpu "num=1:gmem=40G:j_exclusive=yes"` takes any A40, A100 or H100 and no V100 (V100s have no bf16).
  `gmodel=NVIDIAH10080GBHBM3` / `NVIDIAA100_SXM4_80GB` / `NVIDIAA40` pins one model.
- LSF sets `CUDA_VISIBLE_DEVICES` to the job's GPU, so torch sees one device; `nvidia-smi` still lists all of them.
- Memory: every DataLoader worker holds its own copy of the dataset index (several GB for the 5 M tensors of
  the earlier HG008T_Illumina SNV build). `predict` with 12 workers passed 64 GB and swapped; `predict.sh` uses 4
  (`PREDICT_WORKERS`). LSF's `-M` does not kill the job at the limit: it slows down. Watch `MEM` in `bjobs -l`.
- Free GPUs usually come with few free slots (the hosts' CPUs are taken by CPU jobs): a 16-slot GPU job pended
  where a 1-slot job started in minutes. Ask for few slots unless the host is free.
- See what is free: `bhosts -gpu -w`, and `bhosts -w <host>` for the host's state (`ok`, `closed_Full`,
  `closed_Adm`).

## Container and environment

The jobs run in `pytorch/pytorch:2.8.0-cuda12.6-cudnn9-devel` (torch 2.8.0+cu126, Python 3.11, numpy,
zstandard). Its `NVIDIA_REQUIRE_CUDA` accepts the 550 datacenter driver (forward compatibility), and torch 2.8 is
the version the code was developed on. timm and pysam are installed once with `pip install --user` into
`$PANSOMA_RUNS/env/pyuser` (`setup_env.sh`); `env.sh` sets `PYTHONUSERBASE` for every job. RTG Tools 3.13 (with
its own JRE) is in `$PANSOMA_RUNS/env/rtg-tools-3.13` (talkback off in its `rtg.cfg`; give memory as
`rtg RTG_MEM=4g ...`, its default is 90 % of the host's RAM).

- `jiaweiwustk/pansoma` (Docker Hub) is the image of the earlier pipeline: torch 2.1.2+cu121, vg 1.74.1, an
  entrypoint of its own. `pansoma_net_v2` has not been tested with it.
- A host pulls the image (~7.8 GB) the first time it runs it: 1–3 minutes.
- The image has no `git`; `run.sh` reads the commit from `.git`.
- `LSF_DOCKER_VOLUMES` (set in `~/.bashrc`) mounts `/storage2/fs1/epigenome/Active`; `submit.sh` sets it itself
  and adds the truth mounts. `LSF_DOCKER_SHM_SIZE` defaults to 16g in `submit.sh` (DataLoader workers).

## Scripts

`ris/submit.sh [-q QUEUE] [-G GROUP] [-g GPUS] [-n SLOTS] [-m MEM_GB] [-s GPU_SELECT] [-w DEPENDENCY] [-j JOB_GROUP]
JOBNAME SCRIPT [args]` submits `ris/SCRIPT` with the image, the mounts and `-oo $PANSOMA_RUNS/logs/JOBNAME.log`.
Defaults: `-q general -g 1 -n 1 -m 64 -s gmem=40G`; `-g 0` is a CPU job (no `-gpu`, no `gpuhost`). `-w` is bsub's
dependency (e.g. `"done(pnv2_NAME)"`); `-j` the LSF job group, default `$PANSOMA_JOB_GROUP`, else `/<user>/pansoma` when
it exists ("Job groups" below).

| script | what | example |
|---|---|---|
| `setup_env.sh` | once: timm and pysam into the user base, then the CPU unit tests | `ris/submit.sh -g 0 -n 4 -m 16 pnv2_setup setup_env.sh` |
| `gpu_check.sh` | `tests/test_gpu` on a GPU | `ris/submit.sh pnv2_gpu_check gpu_check.sh` |
| `prebuild.sh SET KIND [chunks options]` | index cache and zstd chunk cache on CPU, before a GPU run | `ris/submit.sh -g 0 -n 8 -m 32 pnv2_prebuild prebuild.sh HG008T_Illumina SNV` |
| `run.sh NAME SET KIND [train options]` | train (chr1 left out), then predict chr1 with `best.pth` | below |
| `predict.sh NAME SET KIND [predict options]` | predict chr1 with `<run>/best.pth` (run.sh's last step) | `ris/submit.sh -n 1 -m 64 pnv2_predict_SNV_test predict.sh HG008_Illumina_SNV_test HG008T_Illumina SNV` |
| `relabel.sh SET SNV_MIN_AF` | relabel a set (both kinds) with this checkout's rules and the truth of its current labels; `''` for no SNV floor (the current labels have none); backs up the current labels to `<set>/labels_backup_<date>/` first (CPU, 9 slots: label's recall scan reads in 8 processes besides the main one; ~22 GB, hours) | `ris/submit.sh -g 0 -n 9 -m 32 pnv2_relabel_HG008T_PacBio relabel.sh HG008T_PacBio ''` |
| `snv_pipeline.sh NAME SET [train options]` (login node) | submits `run.sh` with the recipe of `KIND` (`SNV`, the default, or `INDEL`; `recipes.sh`) on an H100 of `epigenome` and `vcf.sh` after it (`-w done(...)`); `SLOTS` / `MEM` (default 7 / 128 GB; the runs of 2026-10-03 to 10-12 used `MEM=192`, `SLOTS=5` for INDEL, and 16 / 256 on 4 GPUs), `GPUS=N` trains under DDP (`--batch-size` per GPU), `VAL_WORKERS` (default 6, per rank), `QUEUE` / `SELECT`, `WAIT=` a dependency of the GPU job, `NOVCF=1` the GPU job only (vcf.sh later) | `KIND=INDEL SLOTS=5 MEM=192 ris/snv_pipeline.sh HG008_PacBio_INDEL_7fc7396 HG008T_PacBio` |
| `joint_pipeline.sh NAME SET[,SET...] [train options]` (login node) | one model on several sets: submits `joint_run.sh` (train on all sets, chr1 left out; predict chr1 of each set into `<run>/test_chr1_<SET>/`) and one `vcf.sh --set SET` per set after it; the same `KIND`, `GPUS`, `SLOTS`, `MEM`, `WAIT`, `VAL_WORKERS` as `snv_pipeline.sh` | `KIND=INDEL ris/joint_pipeline.sh COLO829_Illumina_INDEL_joint COLO829T_Illumina,COLO829T_Illumina_P75,COLO829T_Illumina_P50,COLO829T_Illumina_P25` |
| `recipes.sh` | sourced by the two pipelines: `RECIPE` = the v2 README's train options for `KIND` | |
| `vcfeval_only.sh RUN KIND SRC OUT` | the PoN + rtg steps of `vcf.sh` only, on `<run>/SRC/*.linear.vcf.gz` into `<run>/OUT/`: re-scores earlier calls with this checkout's PoN script | `ris/submit.sh -g 0 -n 4 -m 32 pnv2_rescore vcfeval_only.sh HG008_Illumina_INDEL_small INDEL vcf_chr1 vcf_chr1_pon7fc7396` |
| `vcf.sh [--set SET] RUN KIND [THRESHOLD]` | graph VCF → GRCh38 VCF (the set's graph index and GRCh38 path, by the names in its manifests: HPRC v1.1 d9, v2.1 d46 or the GRCh38-only graph, in `Dataset/pansoma_v2_tensors/graph_index/`) → PoN (the rule of each kind) + rtg vcfeval (in the BED and without it) into `<run>/vcf_chr1/` (`--set`: a joint run's `test_chr1_<SET>/` into `vcf_chr1_<SET>/`); THRESHOLD (a number, or `predict` for predict's own call) overrides the PASS threshold, the checkpoint's validation recall-0.9 threshold | `ris/submit.sh -g 0 -n 4 -m 32 pnv2_vcf vcf.sh HG008_Illumina_SNV_test SNV` |

`run.sh` passes the train options as given; the recipes are in the v2 README ("SNV recipe", "INDEL recipe";
`snv_pipeline.sh` passes them). With `-g 2` or more it trains under `torchrun --standalone --ddp`.

## Workflow

```bash
cd /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma

# 0. once per environment: packages and CPU tests; then the GPU tests
ris/submit.sh -g 0 -n 4 -m 16 pnv2_setup setup_env.sh
ris/submit.sh pnv2_gpu_check gpu_check.sh

# 1. a short test: train on chr21 + chr22, 2 epochs, predict chr1 (any GPU with >= 40 GB, 1 slot)
ris/submit.sh -n 1 -m 64 pnv2_HG008_Illumina_SNV_test run.sh HG008_Illumina_SNV_test HG008T_Illumina SNV \
    --train-chroms chr21 chr22 --epochs 2 --read-threads 4 --num-workers 4 --val-workers 4

# 2. a full run on chr2-22: build the caches on CPU first, then train on an H100 of the epigenome queue
ris/submit.sh -g 0 -n 8 -m 32 pnv2_prebuild prebuild.sh HG008T_Illumina SNV
ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 7 -m 192 pnv2_HG008_Illumina_SNV_recipe \
    run.sh HG008_Illumina_SNV_recipe HG008T_Illumina SNV \
    --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75 \
    --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
    --read-threads 8 --val-workers 6

# 3. calls and rtg vcfeval of the chr1 predictions (a third argument THRESHOLD overrides the PASS threshold: a number
#    or predict); step 2's GPU job and step 3 in one command: MEM=192 ris/snv_pipeline.sh HG008_Illumina_SNV_recipe HG008T_Illumina
ris/submit.sh -g 0 -n 4 -m 32 pnv2_vcf_SNV_recipe vcf.sh HG008_Illumina_SNV_recipe SNV

# watch
bjobs -w; bpeek <jobid>; tail -f $PANSOMA_RUNS/logs/<JOBNAME>.log $PANSOMA_RUNS/<NAME>/train.log
bkill <jobid>
```

`python3 -m pansoma_net_v2.summary [NAME ...]` (run from `machine_learning/` with `PANSOMA_RUNS` set) needs only
the standard library, so it runs on the login node too.

`bmod` of a pending job fails ("You must supply an 'application'"): the esub wants `-a` again. Kill and
resubmit instead.

An interactive shell in the same container, on a GPU (the form of the earlier ConvNeXt jobs; not yet run with this
image):

```bash
LSF_DOCKER_SHM_SIZE=16g bsub -Is -G compute-epigenome-condo -q epigenome-interactive -R 'gpuhost' \
    -gpu "num=1:j_exclusive=yes" -a 'docker(pytorch/pytorch:2.8.0-cuda12.6-cudnn9-devel)' /bin/bash
source /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma/ris/env.sh
cd $REPO/machine_learning && python -m pansoma_net_v2.train --help
```

## Measured on RIS

The test run of step 1 (`HG008_Illumina_SNV_test`, 2026-09-29, A40 on compute1-exec-359, 1 slot, commit 0cde04a,
with `--non-fraction 0.25 --class-weights sqrt --batch-size 256 --lr 2e-4` and the default 199.6 M model; the PoN
rule and graph_vcf's PASS threshold of that commit, both changed since):

| step | time |
|---|---|
| CPU unit tests (`setup_env.sh`, 4 slots) | 65 tests in 329 s, 9 skipped (GPU, real data, rtg) |
| GPU tests (`gpu_check.sh`) | 4 tests in 27–81 s |
| index cache of HG008T_Illumina SNV (5.0 M tensors, 16 GB of summaries), first run only | 4.5 min |
| zstd chunk cache of chr21 + chr22 (7 shards, 33 GB) | 17 s |
| statistics from 20,000 training tensors | 20 s |
| epoch 1 (30,662 tensors, with ~2 min of torch.compile) / epoch 2 | 338 s / 65 s (472 tensors/s) |
| GPU memory, batch 256, compiled | 24.1 GiB allocated, 25.4 GiB reserved |
| `best.pth`, `last.pth` (with the optimizer state) | 2.4 GB each |
| rtg SDF of GRCh38 (`rtg format`, 4 slots) | 75 s |
| `predict.sh`: chr1, 468,721 tensors, 4 workers, 1 slot | 17 min (~450 tensors/s), 30 GB peak memory |
| `vcf.sh`: graph VCF, GRCh38 VCF, PoN, rtg vcfeval with and without the BED (4 slots) | 4 min, 1.9 GB |

The A40 trains at about a third of an H100's speed (1,404 tensors/s on tequila). storage2 read the 7 shards
for the chunk cache in 17 s, much faster than the old cluster's BeeGFS.

The numbers that do not depend on the model agree with the old cluster's HG008 Illumina chr1 values of that build (v2
README at 0cde04a): 697 SNP truth
alleles in the BED and 702 without it, a tensor ceiling of 0.950, the PoN tagging 100 of the 702 truth records. So
the truth files, PoNs and rtg on RIS give what they gave on the old cluster. The model scores themselves (chr1
somatic F1 0.044 per tensor; rtg with the PoN in the BED: TP 38, FP 309) are those of a 2-epoch model trained on
187 somatic tensors, a smoke test only.

`vcfeval` prints `[W::hts_idx_load3] The index file is older than the data file` for the four PoNs: their copies
got mtimes a few seconds apart. It is harmless (`touch` the `.tbi` files to silence it).

The same test on an H100 of the `epigenome` queue (`HG008_Illumina_SNV_test_H100`, 2026-10-01, compute1-exec-419,
`-n 16 -m 128`, caches already built; `gpu_check.sh` passed there too):

| step | A40, 1 slot | H100, 16 slots |
|---|---|---|
| epoch 1 (with torch.compile) / epoch 2 | 338 s / 65 s (472 tensors/s) | 202 s / 20 s (1,561 tensors/s) |
| predict chr1 (4 workers) | 17 min | 5 min |
| whole `run.sh` job | | 14 min |

The two runs give the same losses (epoch 1: 0.4072 / 0.4071) and the same chr1 result (somatic F1 0.044). The
H100 of exec-419 trains a little faster than tequila's H100 NVL (1,404 tensors/s).

After pulling 320e1dd (2026-10-01): the CPU unit tests pass (72, 8 skipped), and `vcf.sh HG008_Illumina_SNV_test SNV
predict` runs with the new PoN rule (allele match for all four PoNs, gnomAD and CoLoRSdb at AF ≥ 1e-4): it tags 61 of
the 702 chr1 truth records instead of 100, so the recall bound after the PoN is 0.913 instead of 0.858.

## Runs of 2026-10-01: the SNV recipe on five sets

Submitted 2026-10-01 16:27 with `snv_pipeline.sh` (repo at 320e1dd plus `ris/`): train (SNV recipe, chr1 left out) and
predict chr1 on one H100 each, then `vcf.sh`. The HG008T sets are relabelled first with `relabel.sh SET 0.08` (the
upstream 2026-09-30 relabel; the RIS copies had the 09-28 labels: Illumina 0.07, PacBio / ONT none). COLO829T_ONT waits
for its copy.

| run | set | SNV floor | jobs (LSF) |
|---|---|---|---|
| `COLO829_Illumina_SNV_recipe` | COLO829T_Illumina (2,505,600 tensors) | 0.07 | 789173 GPU → 789174 vcf |
| `COLO829_fiberseq_SNV_recipe` | COLO829T_fiberseq (1,715,829) | none | 789175 → 789176 |
| `HG008_Illumina_SNV_recipe` | HG008T_Illumina (5,000,888) | 0.08 | 789264 relabel → 790068 (256 GB) → 790069 |
| `HG008_PacBio_SNV_recipe` | HG008T_PacBio (1,180,343) | 0.08 | 789268 → 789269 → 789270 |
| `HG008_ONT_SNV_recipe` | HG008T_ONT (2,372,506) | 0.08 | 789271 → 789272 → 789273 |

The first relabels (789177, 789180, 789183) stopped with exit 127 before labelling (`/usr/bin/time` is not in the image;
`relabel.sh` now uses bash's `time`); their dependants were killed and the chains resubmitted. LSF binds `done(NAME)`
to the newest job of that name (`bjdepinfo -p JOBID` shows the parent).

Results: `$PANSOMA_RUNS/<run>/` (`train.log`, `metrics.jsonl`, `best.pth`, `test_chr1/`, `vcf_chr1/eval_bed/report.txt`,
`vcf_chr1/eval_nobed/report.txt`); logs `$PANSOMA_RUNS/logs/pnv2_<run>.log`, `pnv2_vcf_<run>.log`.

`relabel.sh HG008T_PacBio 0.08` (789268, 2026-10-01, rules 197b5d25bbf4) reproduces the upstream 2026-09-30 labels: all
8,234 SNV tensors that represent a truth allele in `analysis/tensor_recall_20260930/labels/HG008T_PacBio.tsv` have the
same label and reason on RIS (8,231 label 1, 3 `below_snv_min_af`). SNV totals: somatic 8,763, germline 348,027,
non 73,461, ignore 750,092 (09-28: 8,836 / 350,827 / 86,677 / 734,003).
`relabel.sh HG008T_ONT 0.08` (789271) also matches upstream: 8,275 of 8,275 truth-representing SNV tensors (7
`below_snv_min_af`). SNV totals: somatic 8,789, germline 369,101, non 197,581, ignore 1,797,035.
`relabel.sh HG008T_Illumina 0.08` (789264) matches upstream too: 8,279 of 8,279 (6 `below_snv_min_af`). SNV totals:
somatic 10,251, germline 314,363, non 1,876,103, ignore 2,800,171. Its GPU job waited for an H100: at 17:00 all 8 on
exec-419 were taken (4 by these runs, 4 by another user).

**Memory (2026-10-01 17:15, `bjobs -l` MEM; LSF's SWAP column here always equals MEM, so it says nothing).** With 7 slots,
`--val-workers 6` and the 4-worker predict, the training jobs used: COLO829T Illumina 121 GB (2.5 M tensors), HG008T ONT
92 GB (max 128), HG008T PacBio 50 GB (max 128, and the slowest run at ~1,500 tensors/s against ~3,600 for the others,
maybe throttled at the limit), COLO829T fiberseq 46 GB. GPU: 48.8 GiB allocated / 54.1–54.3 GiB reserved on every run
(nvidia: ~55.5 GB of 80). HG008T Illumina (5.0 M tensors, not started) was resubmitted with `MEM=256` (790068 → 790069).
For the 128 GB default of `snv_pipeline.sh`, use `MEM=256` on the large sets.

Results of the 2026-10-01 runs (old tensor build: ch6 int8-count-linear100-log2, SNV build AF 0.06, labels 197b5d25bbf4
with HG008T floor 0.08 / COLO829T Illumina 0.07). chr1, PASS at the validation recall-0.9 threshold, rtg vcfeval; for
COLO829T the somatic BED covers every chromosome whole, so in-BED = no BED. Deleted 2026-10-03 with the run directories
when the sets were rebuilt.

| run | best epoch | truth | raw P / R / F1 | PoN P / R / F1 | PoN, no BED P / R / F1 |
|---|---|---|---|---|---|
| COLO829_fiberseq | 9 | 2,443 | 0.510 / 0.745 / 0.605 | 0.881 / 0.726 / 0.796 | same |
| COLO829_Illumina | 6 | 2,443 | 0.464 / 0.737 / 0.570 | 0.827 / 0.719 / 0.769 | same |
| HG008_PacBio | 9 | 697 | 0.111 / 0.859 / 0.196 | 0.495 / 0.829 / 0.620 | 0.484 / 0.829 / 0.611 |
| HG008_ONT | 8 | 697 | 0.101 / 0.857 / 0.181 | 0.474 / 0.829 / 0.603 | 0.462 / 0.829 / 0.594 |
| HG008_Illumina | 5 | 697 | 0.110 / 0.842 / 0.195 | 0.469 / 0.816 / 0.596 | 0.458 / 0.816 / 0.587 |

## Runs of 2026-10-03: the SNV recipe on the six rebuilt sets

The sets were rebuilt upstream (9d6bd8d: channel 6 relative to `--haplotypes 90`, storage `int8-count-haplotypes100`;
build AF 0.08 for SNV and INDEL; labels `5430533aa100` without AF floors) and copied to RIS on 2026-10-02/03; every
set's truth files match the manifests' sha256. The earlier run directories, index and chunk caches and logs were
deleted (only `env/` kept). Submitted 2026-10-03 18:35 with `MEM=192 ris/snv_pipeline.sh <run> <set>` (repo at 9d6bd8d
plus `ris/`), one H100 each, all eight free; CPU unit tests of 9d6bd8d in `pnv2_setup_9d6bd8d` (914624).

| run | set (SNV tensors) | jobs (GPU → vcf) |
|---|---|---|
| `COLO829_fiberseq_SNV_h90` | COLO829T_fiberseq (1,576,803) | 914625 → 914644 |
| `COLO829_ONT_SNV_h90` | COLO829T_ONT (2,464,433) | 914627 → 914645 |
| `COLO829_Illumina_SNV_h90` | COLO829T_Illumina (1,701,348) | 914629 → 914646 |
| `HG008_PacBio_SNV_h90` | HG008T_PacBio (1,100,923) | 914631 → 914647 |
| `HG008_ONT_SNV_h90` | HG008T_ONT (2,093,336) | 914648 → 914649 |
| `HG008_Illumina_SNV_h90` | HG008T_Illumina (2,851,868) | 914650 → 914651 |

`COLO829T_fiberseq/tensors/SNV` on RIS also holds five shards of the previous build (chr1_00003, chr8_00003,
chr9_00003, chr13_00003, chr16_00002 and their `_labels.npy`, dated 2026-10-01; 32,481 tensors) that the new
`manifest.json` does not list. `data.KindIndex` and the chunk cache take the shards from the manifest, so they are not
read.

**Job groups.** RIS puts jobs in `/<user>/default`, which runs at most 5 at once (`bjgroup`): with the tests and four
GPU jobs, the other two pended ("The specified job group has reached its job limit"). `bgadd -L 10 /jiawei.s/pansoma`
made a group of our own, and `submit.sh` now submits there (`-j`, `PANSOMA_JOB_GROUP`, default `/<user>/pansoma` when it
exists). `bmod -g` is refused by the esub (it wants `-a` and `-G`, and would take the Docker mounts from the current
shell), so the eight pending jobs were killed and resubmitted (the vcf jobs of the four running ones with
`-w "done(pnv2_<run>)"`).

INDEL on the same six sets, submitted 2026-10-03 ~20:00 while the last two SNV runs trained (`KIND=INDEL SLOTS=5
MEM=192 ris/snv_pipeline.sh <run> <set>`): the same SNV recipe (no INDEL recipe is documented; the newest upstream INDEL
runs, `HG008_Illumina_INDEL_nopartial_b1024*`, also leave the partial labels out). Partial labels are a large part of the
INDEL positives (HG008T PacBio 2,742 of 6,840 somatic; `analysis/indel_partial_tensors_20260930`: on HG008T Illumina
chr1, 97 of 128 partial tensors are in homopolymers ≥ 6, median AF 0.17). No PoN for INDEL (`vcf.sh` then).

| run | set (INDEL tensors; somatic) | jobs (GPU → vcf) |
|---|---|---|
| `COLO829_fiberseq_INDEL_h90` | COLO829T_fiberseq (1,782,626; 1,639) | 915178 → 915179 |
| `COLO829_ONT_INDEL_h90` | COLO829T_ONT (4,493,849; 1,755) | 915180 → 915181 |
| `COLO829_Illumina_INDEL_h90` | COLO829T_Illumina (146,555; 1,349) | 915182 → 915183 |
| `HG008_PacBio_INDEL_h90` | HG008T_PacBio (1,923,786; 6,840) | 915184 → 915185 |
| `HG008_ONT_INDEL_h90` | HG008T_ONT (2,511,131; 6,383) | 915186 → 915187 |
| `HG008_Illumina_INDEL_h90` | HG008T_Illumina (151,741; 5,239) | 915188 → 915189 |

**SNV results on the rebuilt sets** (chr1, PASS at the validation recall-0.9 threshold, rtg vcfeval; COLO829T's BED
covers every chromosome whole, so in-BED = no BED):

| run | best epoch | truth (BED) | raw P / R / F1 | PoN P / R / F1 | PoN best F1 | PoN, no BED P / R / F1 | PASS without GRCh38 (PoN P if false) |
|---|---|---|---|---|---|---|---|
| COLO829_fiberseq | 12 | 2,443 | 0.445 / 0.736 / 0.555 | 0.857 / 0.717 / 0.781 | 0.784 | 0.857 / 0.717 / 0.781 | 1,018 (0.572) |
| COLO829_ONT | 8 | 2,443 | 0.381 / 0.705 / 0.494 | 0.727 / 0.686 / 0.706 | 0.717 | 0.727 / 0.686 / 0.706 | 1,675 (0.421) |
| COLO829_Illumina | 6 | 2,443 | 0.469 / 0.729 / 0.571 | 0.860 / 0.711 / 0.778 | 0.780 | 0.860 / 0.711 / 0.778 | 425 (0.710) |
| HG008_PacBio | 7 | 697 | 0.108 / 0.857 / 0.191 | 0.497 / 0.828 / 0.621 | 0.641 | 0.488 / 0.828 / 0.614 | 211 (0.421) |
| HG008_ONT | 8 | 697 | 0.102 / 0.855 / 0.182 | 0.466 / 0.826 / 0.596 | 0.620 | 0.454 / 0.826 / 0.586 | 1,409 (0.218) |
| HG008_Illumina | 6 | 697 | 0.107 / 0.851 / 0.190 | 0.463 / 0.822 / 0.592 | 0.612 | 0.453 / 0.822 / 0.584 | 570 (0.317) |

- Against the old build (2026-10-01): with the PoN, HG008T PacBio 0.620 → 0.621, ONT 0.603 → 0.596, Illumina 0.596 →
  0.592; COLO829T fiberseq 0.796 → 0.781 (recall bound by the lower tensor ceiling, 0.839 → 0.816: build AF 0.08 drops
  low-VAF subclonal truth), Illumina 0.769 → 0.778. COLO829T ONT is new.
- The off-reference rescue (`predict --offref-site90 10`, ch6 raised to 100 on `int8-count-haplotypes100`) makes many
  more PASS calls without a GRCh38 position than on the old build (COLO829T fiberseq 263 → 1,018); rtg does not score
  them, and if they are false the PoN precision falls to the last column.
- `HG008_PacBio` ran at ~3,600 tensors/s with `MEM=192` (1,500/s at 128 GB on 2026-10-01: throttled at the limit).

**INDEL results on the rebuilt sets** (chr1, no PoN, PASS at the checkpoint's validation recall-0.9 threshold, rtg
vcfeval in the BED; HG008T no-BED differs by ≤ 0.002):

| run | best epoch (val F1) | PASS threshold | truth | PASS evaluated | TP | P / R / F1 | best F1 (its P) | ceiling | PASS without GRCh38 |
|---|---|---|---|---|---|---|---|---|---|
| COLO829_fiberseq | 2 (0.082) | 0.0841 | 131 | 20,957 | 58 | 0.003 / 0.443 / 0.005 | 0.094 (0.077) | 0.481 | 5,237 |
| COLO829_ONT | 9 (0.069) | 0.0018 | 131 | 61,998 | 50 | 0.001 / 0.382 / 0.002 | 0.035 (0.025) | 0.382 | 11,839 |
| COLO829_Illumina | 2 (0.114) | 0.0757 | 131 | 3,860 | 52 | 0.013 / 0.397 / 0.026 | 0.075 (0.056) | 0.435 | 403 |
| HG008_PacBio | 4 (0.193) | 0.1444 | 587 | 9,322 | 241 | 0.026 / 0.411 / 0.049 | 0.149 (0.114) | 0.469 | 1,424 |
| HG008_ONT | 6 (0.150) | 0.0231 | 587 | 50,694 | 237 | 0.005 / 0.404 / 0.009 | 0.090 (0.079) | 0.434 | 9,760 |
| HG008_Illumina | 3 (0.217) | 0.1831 | 587 | 3,507 | 213 | 0.061 / 0.363 / 0.104 | 0.150 (0.107) | 0.402 | 544 |

- The SNV recipe does not work for INDEL. The models peak early (epochs 2–9) at validation somatic F1 0.07–0.22, with
  65–194 somatic validation tensors among 6,000–232,000.
- The recall-0.9 validation threshold is then tiny (0.0018 on COLO829T ONT), so tens of thousands of tensors PASS.
  Even the best point of each curve is F1 ≤ 0.15.
- Recall is bound by the tensors (ceiling 0.38–0.48: a truth INDEL without a tensor is a miss at every threshold).
- HG008T Illumina agrees with upstream (`nopartial_b1024_nopc`: R 0.365, ceiling 0.400).

## Runs of 2026-10-06: the INDEL recipe (45d56ad) on the six sets

After pulling 45d56ad: `snv_pipeline.sh` takes the INDEL recipe for `KIND=INDEL` (8.7 M model `--depths 2 2 6 2 --dims
64 128 256 512 --epochs 12 --batch-size 1024 --lr 4e-4 --class-weights sqrt --non-fraction 1.0 --block v2`, partials
ignored), and `vcf.sh` passes the PoN for both kinds (the PoN script's INDEL rule at 45d56ad: a flat gnomAD /
CoLoRSdb exact-allele AF ≥ 0.01). Submitted with `KIND=INDEL SLOTS=5 MEM=192 ris/snv_pipeline.sh <S>_<P>_INDEL_small <set>`:
COLO829_fiberseq 938392 → 938393, COLO829_ONT 938394 → 938395, COLO829_Illumina 938396 → 938397,
HG008_PacBio 938398 → 938399, HG008_ONT 938400 → 938401, HG008_Illumina 938402 → 938403.

## INDEL: a general method to F1 >= 0.15 (2026-10-06/07)

Goal set by the user: INDEL F1 >= 0.15 on chr1 for all six rebuilt sets with a general method (one recipe and one
decision rule for all sets, nothing tuned on chr1, nothing tailored to these datasets). Everything below was chosen
on chr2 (held out of training with `--val-chroms chr2`) and written to `$PANSOMA_RUNS/final_indel/PREREGISTRATION.md`
before the single chr1 evaluation.

**Method.**
1. CNN: the upstream INDEL recipe (45d56ad, 8.7 M model), trained per set with chr2 held out; three seeds (0/1/2), fused
   by the normalized geometric mean of their class probabilities (`tools/fuse_preds.py`, `ens_fuse.sh`; runs `ixe_<set>`).
2. Call-level rescorer (`tools/rescore_indel.py`, Strelka2 EVS-style second stage): one logistic model over
   label-free per-call evidence: logit p_somatic / p_germline / p_non; population AF of the exact allele and of any
   length-changing allele within 6 bp (gnomAD, CoLoRSdb); repeat tract length and homopolymer >= 6; calls within 5 kb /
   20 kb (Mutect2's clustered events); LEFTSHIFT. Fitted once on chr2 of all six sets (equal weight per set),
   PoN floor (0.02) and ONE shared threshold (0.145) chosen on chr2; applied unchanged to chr1.
3. PoN: one flat gnomAD / CoLoRSdb exact-allele AF floor for every INDEL, at that value (`tools/pon_tag.py`; SNV rule
   unchanged), then rtg vcfeval.

**How it was chosen (chr2 only).** Upstream recall-0.9 threshold -> chr2-tuned threshold (`exp_calls.sh` general) ->
rescorer (block CV mean 0.169 -> 0.256 on the seed-0 models, wins in 92-100% of CV repeats in every set and family).
Rejected on chr2 by rules fixed beforehand: focal loss (ixf), pooling per platform (ixp) or over all six sets (ixp_all;
chr2 CV -0.005), per-set PoN/read-AF floors (worse on chr1 than the general rule for COLO829T Illumina), rescorer
feature set v2 (ClairS-TO sequence entropy, Mutect2 repeat-unit/slippage terms, population-absent flag: CV +0.011 but
leave-one-sample-out mean 0.233 -> 0.222). Adopted: the 3-seed ensemble (chr2 CV 0.256 vs 0.243 for the mean single
seed, leave-one-sample-out 0.233 vs 0.204). Code reviewed adversarially (56 agents): no critical/major defect; fixed
the minor ones (PoN float32 tolerance in tune_indel, entropy normalization, errexit in $(...)).

**chr1 results** (`python3 ris/tools/final_summary.py rscF_ixe_v1 ...`; PoN, in the BED; 95% CI by 5 Mb block bootstrap):

| set | upstream recipe | ensemble, threshold only | **ensemble + rescorer** | 95% CI | rescorer from the other sample | single seeds |
|---|---|---|---|---|---|---|
| HG008T Illumina | 0.218 | 0.239 | **0.275** | 0.236-0.311 | 0.210 | 0.223-0.285 |
| HG008T PacBio | 0.089 | 0.237 | **0.274** | 0.236-0.310 | 0.198 | 0.221-0.286 |
| HG008T ONT | 0.016 | 0.149 | **0.180** | 0.150-0.210 | 0.162 | 0.152-0.185 |
| COLO829T Illumina | 0.092 | 0.175 | **0.239** | 0.159-0.320 | 0.244 | 0.179-0.255 |
| COLO829T fiberseq | 0.015 | 0.114 | **0.204** | 0.135-0.277 | 0.195 | 0.179-0.210 |
| COLO829T ONT | 0.003 | 0.057 | **0.111** | 0.047-0.178 | 0.108 | 0.096-0.113 |

Five of six sets reach F1 >= 0.15; COLO829T ONT does not (about 35 PoN-pass true calls against tens of thousands of
homopolymer-error candidates). Without the BED the F1 is 0.000-0.012 lower. The rescorer fitted only on the other
sample's chr2 (its CNN is still the set's own) keeps 5 of 6 at >= 0.15. Caveats: chr2 also selected the CNN checkpoints;
early exploration looked at chr1 calls of the earlier `*_INDEL_small` models; COLO829T has 131 chr1 truth INDELs.
Feature set v2 (rejected on chr2) gives 0.290/0.295/0.179/0.250/0.199/0.111 on chr1, within the intervals above.

Runs: `ixb_<set>` / `ixb_s1_<set>` / `ixb_s2_<set>` (CNNs), `ixe_<set>` (ensemble), `rsc3_*` (chr2 fits),
`rscF_*` (chr1 evaluation); scripts `exp_submit.sh`, `exp_train.sh`, `exp_calls.sh`, `ens_fuse.sh`, `exp_rescore.sh`;
tools `rescore_indel.py`, `fuse_preds.py`, `tune_indel.py`, `pon_tag.py`, `train_focal.py`, `fit_compare.py`,
`exp_summary.py`, `final_summary.py`.
