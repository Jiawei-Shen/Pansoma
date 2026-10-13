# COLO829T purity dilutions: one joint model per platform and kind, tested on chr1 of each purity

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
threshold, one threshold for all four purities. In the BED after the PoN, PASS F1 (the chr1 curve's best F1, a
test-set choice, in brackets):

| platform | kind | P100 | P75 | P50 | P25 | P100 single-set model (earlier run) |
|---|---|---:|---:|---:|---:|---:|
| Illumina | SNV | **0.675** (0.784) | **0.840** (0.841) | **0.833** (0.834) | **0.804** (0.805) | 0.778 (0.780) |
| Illumina | INDEL | **0.121** (0.163) | **0.191** (0.308) | **0.232** (0.330) | **0.192** (0.233) | 0.115 (0.195), lr 4e-4 |
| fiberseq | SNV | **0.779** (0.790) | **0.687** (0.703) | **0.665** (0.671) | **0.720** (0.729) | 0.781 (0.784) |
| fiberseq | INDEL | **0.007** (0.068) | **0.004** (0.047) | **0.004** (0.053) | **0.003** (0.073) | 0.016 (0.120), lr 4e-4 |
| ONT | SNV | **0.677** (0.681) | **0.654** (0.665) | **0.639** (0.654) | **0.666** (0.711) | 0.706 (0.717) |
| ONT | INDEL | **0.001** (0.054) | **0.001** (0.051) | **0.001** (0.068) | **0.001** (0.060) | 0.003 (0.071), lr 4e-4 |

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

**2. Illumina dilutions carry a read-source signature.** The base quality of the Illumina reads shifts with purity (BQ
mode 30 at P100, 32-33 at P25: the COLO829BL normal reads have higher BQ than the COLO829T tumor reads;
`bq_hist.py`). In a dilution every ALT read of a somatic SNV comes from the tumor read set, while germline ALT reads
come from both, so the ALT reads' BQ marks somatic variants. On chr2 (training data, `bq_alt_ref.py` ->
`bq_alt_ref.json`), the mean BQ of the A1 (ALT) rows minus the REF rows at the site:

| set | somatic median | germline median | AUC of that difference alone (somatic vs germline) |
|---|---:|---:|---:|
| Illumina P100 / P75 / P50 / P25 | +0.65 / -0.08 / -0.50 / -0.84 | +0.28 / +0.23 / +0.21 / +0.19 | 0.436 / 0.567 / 0.649 / 0.706 |
| fiberseq P100 / P75 / P50 / P25 | +0.77 / +0.72 / +0.77 / +0.90 | +0.44 / +0.38 / +0.37 / +0.33 | 0.456 / 0.459 / 0.453 / 0.439 |

fiberseq reads use the same 7 binned BQ values in tumor and normal, so its dilutions have no such cue. The Illumina
dilution scores (precision 0.97, chr1 somatic AP 0.98-0.99 against 0.70 at P100; validation somatic AP 0.961) are
therefore likely optimistic for a real low-purity tumor: part of what the model separates in an in-silico mix is the
sequencing run of the reads.

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
at P25, 0.2-0.6 (2,148)). Illumina does not show the dip, plausibly because there the BQ cue (2.) marks such ALT reads
as normal-derived.

**4. INDEL.** Illumina: 0.121 / 0.191 / 0.232 / 0.192 on 131 truth INDELs (19-34 true calls after the PoN), the
recall ceiling after the PoN 0.24-0.26 (0.145 at P25). The long-read INDEL models fail at every purity, like the
single-set long-read COLO829T INDEL models: the validation recall-0.9 threshold is 0.007 (fiberseq) and 8.4e-6 (ONT),
so 20,000-37,500 (fiberseq) and 96,000-107,000 (ONT) calls PASS per purity's chr1 (9,900-19,500 and 54,000-65,500
false calls after the PoN), and even the best chr1 point stays at F1 0.047-0.073 (fiberseq), 0.051-0.068 (ONT).

**Training** (`tables.md`): Illumina SNV 4.53 M training tensors (125,916 somatic), 113 min on 4 H100s, best epoch 7,
validation somatic AP 0.961; fiberseq SNV 2.07 M (125,734), 133 min on 4, best epoch 5, AP 0.789; Illumina INDEL 0.43
M (3,439), 10 min on 1, best epoch 6, AP 0.576; fiberseq INDEL 6.60 M (3,965), 80 min on 2, best epoch 3, AP 0.052;
ONT SNV 3.07 M (124,070), 190 min on 4, best epoch 7, AP 0.801; ONT INDEL 11.89 M (3,747), 142 min on 2, best epoch 6,
AP 0.044.

