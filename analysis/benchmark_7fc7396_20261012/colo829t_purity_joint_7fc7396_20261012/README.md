# COLO829T purity dilutions: one joint model per platform and kind, tested on chr1 of each purity

> 2026-10-16: the INDEL numbers here use the 7fc7396 INDEL PoN rule (homopolymer >= 7 bp: AF >= 0.01, else 1e-3). With
> the 71cb51a rule (>= 8 bp, else 1e-4) every INDEL call set was re-scored: `../indel_pon_71cb51a_20261016/tables.md`
> (+0.000 to +0.013 F1, no comparison changes direction).

2026-10-12, on RIS compute1, pipeline of 7fc7396 (the v2 README's SNV recipe and INDEL recipe). For each platform and
kind (SNV, INDEL), one model is trained on the undiluted COLO829T set and its three in-silico tumor / normal (COLO829BL)
dilutions together: P100 + P75 + P50 + P25, chr2-22 of every set. chr1 of each purity is then predicted, called and
scored on its own. This is the design of `analysis/reviewer_response_20261011` (R1.5, R1.6: "trained on COLO829T
dilutions together with the undiluted tumor, results per purity").

Files in this folder:

- `collect.py`: reads the run directories and writes `tables.md`, `summary.json` and `runs/<run>/` (copies of
  `train.log` and of the rtg `report.txt` per purity). Standard library only:
  `PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py`.
- `purity.py` → `purity.md`, `purity.json`: per purity on chr1 before any calling, the truth alleles with a tensor
  (the recall ceiling) and the read AF of the somatic / germline tensors (measured purity).
- `tables.md`: per-purity F1 matrix, training, per-purity rtg details, P100 joint vs the earlier single-set model.
- Read-source artefact checks: `bq_hist.py`, `bq_alt_ref.py` (`bq_alt_ref.json`), `read_source.py` (`read_source.json`:
  the Illumina mixture BAMs' read groups), `bq_relative_check.py` (`bq_relative_check.json`: raw / per-read
  centered / per-read rank BQ), `read_quality_dist.py` (`read_quality_dist.json`: wholly low-quality reads).
- `purity_loh_fp.py` → `loh_fp.md`, `loh_fp.json`: the germline false calls per purity against the AFs of tumor-LOH
  germline variants in a dilution (section 6).

## Data

The twelve COLO829T sets on RIS (`Dataset/pansoma_v2_tensors/COLO829T_<platform>[_P75|_P50|_P25]/tensors`): build
9d6bd8d (`--haplotypes 90`, ch6 `int8-count-haplotypes100`, build AF 0.08 for SNV and INDEL), labels `5430533aa100`
without AF floors, the same truth for every purity (somatic: `COLO829T_somatic_snv_indel.vcf.gz` and the SMaHT v2
easy / difficult / extreme union BED; germline: the COLO829BL dipcall). Checked on 2026-10-12 before training: every
shard of every manifest present with the size its tensor count gives, a labels file for each shard and a summary for
each chromosome.

| platform | purity | SNV tensors (somatic) | INDEL tensors (somatic) |
|---|---|---:|---:|
| Illumina | P100 | 1,701,348 (38,299) | 146,555 (1,349) |
| Illumina | P75 | 1,721,375 (37,434) | 158,624 (1,309) |
| Illumina | P50 | 1,713,275 (36,209) | 165,520 (1,195) |
| Illumina | P25 | 1,720,352 (30,610) | 167,812 (808) |
| fiberseq | P100 | 1,576,803 (37,984) | 1,782,626 (1,639) |
| fiberseq | P75 | 1,370,058 (37,178) | 2,088,886 (1,544) |
| fiberseq | P50 | 1,372,428 (36,016) | 2,275,937 (1,439) |
| fiberseq | P25 | 1,355,692 (30,666) | 2,473,555 (1,195) |
| ONT | P100 | 2,464,433 (37,770) | 4,493,849 (1,755) |
| ONT | P75 | 2,430,644 (36,872) | 3,947,195 (1,621) |
| ONT | P50 | 2,305,979 (35,656) | 3,619,896 (1,469) |
| ONT | P25 | 2,131,060 (29,773) | 3,373,605 (1,300) |

The dilutions share reads with each other and with P100. Training and test never share a chromosome (chr1 is held out
of every set), and the validation blocks are the same node blocks in every set (`data.block_split` depends only on
the seed, chromosome and node block), so no validation region is a training region of another purity.

## What was run

Per platform and kind (`ris/joint_pipeline.sh`, repo 7fc7396 + `ris/`):

1. `ris/prebuild.sh` on CPU for each dilution set and kind: the index cache and the zstd chunk cache of chr2-22 (the
   P100 caches existed from the earlier runs).
2. `ris/joint_run.sh` on H100 80GB GPUs (`epigenome` queue, compute1-exec-419): `pansoma_net_v2.train --tensors <P100>
   <P75> <P50> <P25>` with the recipe (chr1 left out, 5 % of the chr2-22 node blocks validate); then
   `pansoma_net_v2.predict` on chr1 of each set with `best.pth` into `<run>/test_chr1_<SET>/`.
3. `ris/vcf.sh --set <SET>` on CPU, one job per purity: `graph_vcf` (PASS = the checkpoint's validation recall-0.9
   threshold, one threshold for all purities; SNV off-reference rescue on) → `linear_vcf` → `vcfeval` with the four
   PoNs (each kind's rule of 7fc7396) and rtg vcfeval in the BED and without it, into `<run>/vcf_chr1_<SET>/`.

Recipes (`ris/recipes.sh`; `--read-threads 8 --val-workers 3` per rank added):

| kind | model | options |
|---|---|---|
| SNV | 50.4 M (`--depths 3 3 27 3 --dims 96 192 384 768`) | `--epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75 --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth` |
| INDEL | 8.7 M (`--depths 2 2 6 2 --dims 64 128 256 512`) | `--epochs 12 --batch-size 1024 --lr 1e-4 --class-weights sqrt --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth` |

GPUs: only the `epigenome` queue's 8 H100s. The SNV runs train on 4 GPUs each under DDP with `--batch-size 1024` per
GPU (global batch 4,096: a quarter of the steps per epoch of one GPU, the same lr), fiberseq INDEL on 2 (global 2,048),
Illumina INDEL (about 0.5 M training tensors) on one; LSF queues what does not fit.

| run | sets | GPUs, slots, memory | GPU job | vcf jobs (P100, P75, P50, P25) |
|---|---|---|---|---|
| `COLO829_Illumina_SNV_joint_7fc7396_gpu4` | Illumina P100 + P75 + P50 + P25 | 4 H100, 16, 320 GB | 41416 | 41417-41420 |
| `COLO829_fiberseq_SNV_joint_7fc7396_gpu4` | fiberseq P100 + P75 + P50 + P25 | 4 H100, 16, 256 GB | 41421 | 41422-41425 |
| `COLO829_fiberseq_INDEL_joint_7fc7396` | fiberseq P100 + P75 + P50 + P25 | 2 H100, 14, 384 GB | 41440 | 41441-41444 |
| `COLO829_Illumina_INDEL_joint_7fc7396` | Illumina P100 + P75 + P50 + P25 | 1 H100, 5, 128 GB | 41445 | 41446-41449 |
| `COLO829_ONT_INDEL_joint_7fc7396` | ONT P100 + P75 + P50 + P25 | 2 H100, 14, 512 GB | 45054 (after the ONT P75 INDEL transfer; `btop` ahead of the d46 jobs) | 45055-45058 |
| `COLO829_ONT_SNV_joint_7fc7396_gpu4` | ONT P100 + P75 + P50 + P25 | 4 H100, 16, 384 GB | 45407 (after the ONT P75 SNV transfer; `btop`) | 45408, 45409, 45411, 45412 |

The SNV runs first started on 2 GPUs each (`COLO829_<platform>_SNV_joint_7fc7396`, 41314 / 41391) and were stopped
in epoch 1 to run on 4; those partial directories are not used.

Prebuild jobs 41302-41313 (`pnv2_prebuild_COLO829T_<platform>_<P>_<KIND>`); each GPU job waits for its three.
Two prebuilds (Illumina P75 INDEL, fiberseq P75 SNV) built their whole cache ("21 / 50 shards cached") and then
exited 1 in the script's last line, a `du` of the shared chunk directory that hit the other builds' temporary files
(`prebuild.sh` now ignores that). The GPU jobs waiting on them with `done()` (first submission: 41319, 41324, 41329
and their vcf jobs) could never start; they were killed while pending and resubmitted with `ended()` dependencies
(Illumina INDEL without one: its caches were complete). The final submission (the table) needs no dependency.
The three fiberseq INDEL prebuilds (70 / 75 / 79 shards cached) then exited 2 with "error reading input file: Stale
file handle": `prebuild.sh` had been replaced on storage2 (NFS) while they ran. Their caches were complete, and the
`ended()` dependencies let the GPU jobs start.
Submitted 2026-10-12 from the login node:

```bash
cd /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma
for p in Illumina fiberseq; do for k in SNV INDEL; do for f in P75 P50 P25; do
  ris/submit.sh -g 0 -n 8 -m 48 pnv2_prebuild_COLO829T_${p}_${f}_$k prebuild.sh COLO829T_${p}_$f $k
done; done; done
# final submission (the caches built; the SNV run names with _gpu4):
sub() {  # PLATFORM KIND SUFFIX, env GPUS SLOTS MEM
  local p=$1 k=$2 sets="COLO829T_$1,COLO829T_$1_P75,COLO829T_$1_P50,COLO829T_$1_P25"
  KIND=$k VAL_WORKERS=3 ris/joint_pipeline.sh COLO829_${p}_${k}_joint_7fc7396$3 "$sets"
}
GPUS=4 SLOTS=16 MEM=320 sub Illumina SNV _gpu4; GPUS=4 SLOTS=16 MEM=256 sub fiberseq SNV _gpu4
GPUS=2 SLOTS=14 MEM=384 sub fiberseq INDEL;    GPUS=1 SLOTS=5 MEM=128 sub Illumina INDEL
```

**Comparison at P100.** The earlier single-set models of the same P100 sets: SNV `COLO829_<platform>_SNV_h90`
(2026-10-03, the same SNV recipe on one GPU; the SNV PoN rule is unchanged); INDEL `COLO829_<platform>_INDEL_small`
(2026-10-06, the INDEL recipe at lr 4e-4), its chr1 calls re-scored with the 7fc7396 PoN
(`ris/vcfeval_only.sh RUN INDEL vcf_chr1 vcf_chr1_pon7fc7396`, jobs 41370-41372).

## Results (chr1)

All numbers: `tables.md` (per purity: tensors, truth ceiling, rtg raw and after the PoN, in the BED and without it;
COLO829T's BED covers each chromosome whole, so the two agree). PASS = the joint checkpoint's validation recall-0.9
threshold, one threshold for all four purities. In the BED after the PoN. SNV: PASS F1 (the chr1 curve's best F1 in
brackets). INDEL: the chr1 curve's best F1 (PASS F1 in brackets), because the recall-0.9 threshold lands far below the
best point for INDELs; the best point is chosen on the test chromosome, an upper bound for a threshold rule, not a
held-out result:

| platform | kind | P100 | P75 | P50 | P25 | P100 single-set model (earlier run) |
|---|---|---:|---:|---:|---:|---:|
| Illumina | SNV | **0.675** (0.784) | **0.840** (0.841) | **0.833** (0.834) | **0.804** (0.805) | 0.778 (0.780) |
| Illumina | INDEL | **0.163** (0.121) | **0.308** (0.191) | **0.330** (0.232) | **0.233** (0.192) | 0.195 (0.115), lr 4e-4 |
| fiberseq | SNV | **0.779** (0.790) | **0.687** (0.703) | **0.665** (0.671) | **0.720** (0.729) | 0.781 (0.784) |
| fiberseq | INDEL | **0.068** (0.007) | **0.047** (0.004) | **0.053** (0.004) | **0.073** (0.003) | 0.120 (0.016), lr 4e-4 |
| ONT | SNV | **0.677** (0.681) | **0.654** (0.665) | **0.639** (0.654) | **0.666** (0.711) | 0.706 (0.717) |
| ONT | INDEL | **0.054** (0.001) | **0.051** (0.001) | **0.068** (0.001) | **0.060** (0.001) | 0.071 (0.003), lr 4e-4 |

Truth on chr1: 2,443 SNVs and 131 INDELs for every purity. **What each purity gives any model** (`purity.md`, before
calling): the SNV recall ceiling falls from 0.82 (P100) to 0.72 / 0.71 (P25, Illumina / fiberseq): the build's AF 0.08
floor removes more low-VAF truth (filtered 384 -> 482 Illumina). The median somatic SNV tensor AF relative to P100
(the measured purity) is 0.66 / 0.40 / 0.20 for P75 / P50 / P25 on Illumina and fiberseq, 0.70 / 0.44 / 0.22 on ONT
(the reviewer notes: 0.65 / 0.38 / 0.17). chr1 of COLO829T is largely LOH: at P100 the somatic and the germline SNV tensors both have median AF 0.97-0.99.

**1. Illumina SNV P100: the shared threshold, not the ranking.** The pooled validation recall-0.9 threshold is 0.804.
In the dilutions 97-99 % of the somatic truth tensors score above it, at P100 only 66 % (high-AF and low-AF somatic
SNVs alike: 0.659 / 0.668), so P100 recall falls to 0.530 at precision 0.931. At its best chr1 point the P100 curve
reaches 0.784, the single-set P100 model's 0.778. Three quarters of the validation tensors are dilution tensors, so
the pooled threshold follows them; a threshold per purity (from that purity's validation blocks) would be needed.

**2. The Illumina and ONT dilutions carry a read-source artefact (not a biological mechanism; to be fixed).** The base quality of the Illumina reads shifts with purity (BQ
mode 30 at P100, 32-33 at P25: the COLO829BL normal reads have higher BQ than the COLO829T tumor reads;
`bq_hist.py`). In a dilution every ALT read of a somatic SNV comes from the tumor read set, while germline ALT reads
come from both, so the ALT reads' BQ marks somatic variants. On chr2 (training data, `bq_alt_ref.py` ->
`bq_alt_ref.json`), the mean BQ of the A1 (ALT) rows minus the REF rows at the site:

| set | somatic median | germline median | AUC of that difference alone (somatic vs germline) |
|---|---:|---:|---:|
| Illumina P100 / P75 / P50 / P25 | +0.65 / -0.08 / -0.50 / -0.84 | +0.28 / +0.23 / +0.21 / +0.19 | 0.436 / 0.567 / 0.649 / 0.706 |
| fiberseq P100 / P75 / P50 / P25 | +0.77 / +0.72 / +0.77 / +0.90 | +0.44 / +0.38 / +0.37 / +0.33 | 0.456 / 0.459 / 0.453 / 0.439 |
| ONT P100 / P75 / P50 / P25 | +1.48 / -0.30 / -1.26 / -2.13 | +1.32 / +1.29 / +1.22 / +1.21 | 0.493 / 0.610 / 0.657 / 0.689 |

(ONT from `bq_relative_check.py` -> `bq_relative_check.json`, the same test, "raw".) fiberseq reads use the same 7
binned BQ values in tumor and normal, so its dilutions have no such cue; ONT has it as strongly as Illumina. The
Illumina mixture BAMs (`read_source.py` -> `read_source.json`) hold two read groups, `tumor.sub` (NovaSeq X flowcell
LH00266_46_22GM3HLT3) and `normal.sub`, that differ only in BQ (mean 29.43 vs 30.77; the same 30-value alphabet,
MAPQ 59.81 vs 59.84, read length 151, soft-clipped 3.0 % vs 2.6 %). The Illumina (and ONT) dilution scores (Illumina
precision 0.97, chr1 somatic AP 0.98-0.99 against 0.70 at P100; validation somatic AP 0.961) are therefore likely
optimistic for a real low-purity tumor: part of what the model separates in an in-silico mix is the sequencing run of
the reads. A per-read BQ normalization (each read's BQ minus its own median in the window) removes the cue in these
tensors: the AUC is 0.42-0.43 at every Illumina purity and 0.49 at every ONT purity (`bq_relative_check.json`,
"read_centered"), but it also hides wholly low-quality reads (13-33 % of ONT reads sit > 5 Q below the set's typical
read, `read_quality_dist.json`); how to fix it is open.

**3. Long reads (fiberseq, ONT) SNV: P75 and P50 below P25, reproduced and explained.** The pattern of the earlier fiberseq joint model
(0.785 / 0.699 / 0.663 / 0.722 in `analysis/reviewer_response_20261011`) is reproduced (0.779 / 0.687 / 0.665 /
0.720). The extra false calls are germline SNVs in the tumor's LOH regions. Per tensor (chr1, before the PoN), germline
tensors called somatic at the PASS threshold (0.549):

| purity | germline called | where their AF peaks | AF of a germline het whose ALT the tumor lost, (1 - purity) / 2 |
|---|---:|---|---:|
| P100 | 1,802 | 0.8-1.0 (1,794; the tumor kept the ALT: AF ~1) | - |
| P75 | 4,538 | 0.1-0.2 (2,107), 0.0-0.1 (817); 0.8-1.0 (1,186) | 0.125 |
| P50 | 4,470 | 0.2-0.3 (1,776), 0.1-0.2 (997); 0.6-0.8 (878) | 0.25 |
| P25 | 2,480 | 0.3-0.4 (1,107), 0.4-0.6 (519) | 0.375 |

Where the tumor lost the ALT copy of a germline het SNV, the ALT reads come only from the admixed normal, at AF
(1 - purity) / 2. At P75 and P50 that is 0.125 / 0.25, inside the somatic AF range of those mixtures (medians 0.63 /
0.39); at P25 it is 0.375, near ordinary het germline (0.5) and far from the P25 somatic median (0.19). The PoN removes
most of these calls (false calls after it: 230 / 545 / 605 / 440), but the order P100 > P25 > P75 > P50 stays.
ONT shows the same pattern (0.677 / 0.654 / 0.639 / 0.666; germline tensors called at its threshold 0.500: 1,396 at
P100, mostly AF 0.8-1.0; 3,239 at P75, peak 0.1-0.2 (1,401) and 0.0-0.1 (593); 4,011 at P50, 0.1-0.3 (2,582); 2,594
at P25, 0.2-0.6 (2,148)), although ONT also carries the BQ cue (2.). Illumina does not show the dip, plausibly because
there the BQ cue marks such ALT reads as normal-derived.

**4. INDEL, by the best point of each purity's chr1 curve after the PoN** (`tables.md`, "INDEL: the chr1 curve's best
F1": F1, P / R, TP / FP and the p_somatic of the point; 131 truth INDELs, recall ceiling after the PoN 0.24-0.28,
0.145-0.23 at P25):

| platform | P100 | P75 | P50 | P25 | P100 single-set model |
|---|---|---|---|---|---|
| Illumina | 0.163 (P 0.200, R 0.137; 18 TP / 72 FP) | 0.308 (0.509, 0.221; 29 / 28) | 0.330 (0.544, 0.237; 31 / 26) | 0.233 (0.594, 0.145; 19 / 13) | 0.195 (0.200, 0.191; 25 / 100) |
| fiberseq | 0.068 (0.077, 0.061; 8 / 96) | 0.047 (0.029, 0.130; 17 / 573) | 0.053 (0.047, 0.061; 8 / 163) | 0.073 (0.115, 0.053; 7 / 54) | 0.120 (0.151, 0.099; 13 / 73) |
| ONT | 0.054 (0.093, 0.038; 5 / 49) | 0.051 (0.148, 0.030; 4 / 23) | 0.068 (0.133, 0.046; 6 / 39) | 0.060 (0.143, 0.038; 5 / 30) | 0.071 (0.162, 0.046; 6 / 31) |

- Illumina: the dilutions score above P100 (0.308 / 0.330 / 0.233 against 0.163), with precision 0.51-0.59 against
  0.20, like the Illumina SNVs; the read-source BQ cue (2.) applies to INDEL reads as well, so the dilution numbers are
  likely optimistic. At P100 the joint model is below the single-set model (0.163 vs 0.195).
- Long reads: 0.047-0.073 (fiberseq) and 0.051-0.068 (ONT) at every purity, below the single-set P100 models (0.120,
  0.071); the best points hold 4-17 true calls, so the purity trend is within noise. At PASS (the validation recall-0.9
  threshold, 0.007 fiberseq, 8.4e-6 ONT) 20,000-37,500 / 96,000-107,000 calls pass per purity's chr1 and F1 is
  0.001-0.007.

**5. Illumina with BQ normalized per read** (2026-10-13; `ris/tools/readnorm.py`, `ris/joint_run.sh --readnorm`: each
read's BQ minus its own median over its valid cells in the window, over the dataset's BQ std; MAPQ and path count as
before; with `--norm dataset`). Runs `COLO829_Illumina_SNV_joint_readnorm_gpu8` (8 H100, global batch 8,192, jobs
56580-56584) and `COLO829_Illumina_INDEL_joint_readnorm` (1 H100, 56585-56589). In the BED after the PoN:

| | P100 | P75 | P50 | P25 |
|---|---|---|---|---|
| SNV PASS F1 (P / R), joint as before | 0.675 (0.931 / 0.530) | 0.840 (0.968 / 0.742) | 0.833 (0.966 / 0.733) | 0.804 (0.967 / 0.689) |
| SNV PASS F1 (P / R), BQ per read | **0.739** (0.823 / 0.670) | **0.709** (0.737 / 0.683) | **0.678** (0.682 / 0.673) | **0.695** (0.723 / 0.668) |
| SNV best F1: before -> per read | 0.784 -> 0.744 | 0.841 -> 0.709 | 0.834 -> 0.684 | 0.805 -> 0.733 |
| SNV chr1 somatic AP: before -> per read | 0.697 -> 0.706 | 0.979 -> 0.472 | 0.986 -> 0.529 | 0.985 -> 0.872 |
| INDEL best F1: before -> per read | 0.163 -> 0.174 | 0.308 -> 0.136 | 0.330 -> 0.125 | 0.233 -> 0.132 |

Validation somatic AP 0.961 -> 0.813 (SNV), 0.576 -> 0.178 (INDEL). Without the read-source cue the Illumina dilutions
lose their precision of 0.97 (now 0.68-0.74) and their INDEL advantage, and show the long reads' P75 / P50 < P25 pattern
(section 6); P100's PASS F1 rises (the pooled recall-0.9 threshold falls from 0.804 to 0.550 once the dilutions no
longer score near 1), its best point falls (0.784 -> 0.744). Not a clean comparison: the per-dataset z-score and, for
SNV, 8 GPUs (global batch 8,192 against 4,096) changed too; a control (`--norm dataset` without the per-read BQ, 8 GPUs)
would separate them. Training ran at 8,262 tensors/s on 8 GPUs against 13,491 on 4 for the earlier run (186 against
113 min): the per-read median (a sort over 200 x 101 cells per tensor) is slow. Per-read centering also hides wholly
low-quality reads (section 2), so it is a diagnostic here, not yet a method.

**6. Why P25 scores above P75 and P50 on every platform** (`purity_loh_fp.py` -> `loh_fp.md`; Illumina with BQ per read,
fiberseq and ONT the joint models as they are). A heterozygous germline SNV of the patient in a tumor LOH region shows
in a dilution at AF (1 - p) / 2 where the tumor lost its ALT copy ("lost") and (1 + p) / 2 where it kept it ("kept"),
p the tumor read fraction: lost 0.125 / 0.25 / 0.375 and kept 0.875 / 0.75 / 0.625 at P75 / P50 / P25. Germline chr1
tensors called somatic at the PASS threshold (SNV), before the PoN (after it):

| platform | purity | PASS F1 | somatic AF q1 / median / q3 | germline called | of them lost | kept | other |
|---|---|---:|---|---:|---:|---:|---:|
| Illumina | P100 | 0.739 | 0.50 / 0.97 / 1.00 | 1,556 (160) | - | 1,510 (155) | 46 (5) |
| Illumina | P75 | 0.709 | 0.41 / 0.64 / 0.75 | 3,297 (372) | 2,137 (251) | 494 (50) | 666 (71) |
| Illumina | P50 | 0.678 | 0.30 / 0.39 / 0.51 | 4,061 (510) | 1,693 (205) | 912 (87) | 1,456 (218) |
| Illumina | P25 | 0.695 | 0.15 / 0.19 / 0.26 | 2,576 (355) | 883 (112) | 724 (64) | 969 (179) |
| fiberseq | P100 | 0.779 | 0.49 / 0.96 / 0.99 | 1,802 (158) | - | 1,732 (151) | 70 (7) |
| fiberseq | P75 | 0.687 | 0.40 / 0.63 / 0.74 | 4,538 (469) | 2,773 (316) | 1,075 (83) | 690 (70) |
| fiberseq | P50 | 0.665 | 0.29 / 0.39 / 0.51 | 4,470 (523) | 2,114 (214) | 811 (64) | 1,545 (245) |
| fiberseq | P25 | 0.720 | 0.15 / 0.19 / 0.27 | 2,480 (350) | 1,110 (108) | 399 (34) | 971 (208) |
| ONT | P100 | 0.677 | 0.45 / 0.82 / 0.90 | 1,396 (140) | - | 217 (19) | 1,179 (121) |
| ONT | P75 | 0.654 | 0.38 / 0.58 / 0.68 | 3,239 (351) | 1,904 (221) | 303 (22) | 1,032 (108) |
| ONT | P50 | 0.639 | 0.28 / 0.36 / 0.47 | 4,011 (463) | 1,526 (161) | 439 (33) | 2,046 (269) |
| ONT | P25 | 0.666 | 0.14 / 0.18 / 0.25 | 2,594 (337) | 719 (76) | 445 (37) | 1,430 (224) |

("lost" / "kept": within 0.06 of the expected AF; ONT's noisier AF puts more into "other", its P100 kept group at AF ~1
too.) The same on all three platforms, so a property of the sample and the dilution, not of the sequencing:
- P100: the germline false calls are the kept group at AF ~1 (homozygous-looking germline of the LOH regions).
- P50 has the most germline false calls after the PoN (463-523): the lost group's 0.25 lies inside P50's own somatic AF
  (q1 0.28-0.30, median 0.36-0.39).
- P75: the lost group's 0.125 lies below P75's somatic AF but inside the somatic AF of P25 (0.14-0.27) in the same
  training data: the joint model does not know a tensor's purity, and low AF is somatic-like across the four sets.
- P25: the lost group's 0.375 and the kept group's 0.625 straddle ordinary heterozygous germline (0.5) and lie above
  P25's somatic AF (median 0.18-0.19), so fewer pass (76-112 lost after the PoN, against 161-316 at P75 / P50). P25's
  recall ceiling is the lowest (0.69 after the PoN) but its precision the highest of the dilutions on fiberseq and ONT
  (0.786 / 0.670 against 0.741 / 0.659 at P75), hence its F1. Illumina (BQ per read) is the exception at PASS: P75 0.709
  (precision 0.737) against P25 0.695 (0.723); at the curves' best points P25 is ahead there too (0.733 vs 0.709).
- This is not an artefact of in-silico mixing: a real tumor with LOH at intermediate purity has the same lost / kept
  germline. The PoN cannot remove them (rare germline variants, not in the population panels). What separates them
  is outside a single site: a whole LOH segment's heterozygous germline sits at one AF ((1 - p) / 2), somatic variants
  do not; a sample-level purity estimate; or a matched normal.
- INDEL (at each curve's best point; 4-34 true calls): the lost group is again the main germline false call at P75 /
  P50 (Illumina 569 / 293 before the PoN), but the counts are small and the best points' thresholds differ by purity,
  so INDEL does not resolve the purity ordering.

**Training** (`tables.md`): Illumina SNV 4.53 M training tensors (125,916 somatic), 113 min on 4 H100s, best epoch 7,
validation somatic AP 0.961; fiberseq SNV 2.07 M (125,734), 133 min on 4, best epoch 5, AP 0.789; Illumina INDEL 0.43
M (3,439), 10 min on 1, best epoch 6, AP 0.576; fiberseq INDEL 6.60 M (3,965), 80 min on 2, best epoch 3, AP 0.052;
ONT SNV 3.07 M (124,070), 190 min on 4, best epoch 7, AP 0.801; ONT INDEL 11.89 M (3,747), 142 min on 2, best epoch 6,
AP 0.044.

