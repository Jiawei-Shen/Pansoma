# E9a / E9b / E9e: placing the graph-coordinate ("unplaced") calls and scoring the complete output

> **2026-10-19, COLO829T v1.1 added** (12 graph runs now): the COLO829T HPRC v1.1 (d9) runs of the matched recipe, seed 0 (RIS, synced 2026-10-19: SNV `COLO829_*_SNV_7fc7396_ddp4_b1024`, 4 GPUs; INDEL `COLO829_*_INDEL_7fc7396`, lr 1e-4), placed under `reviewer_cpu_20261013/<run>/` (best.pth -> the RIS copy, test_chr1 copied, the predictions' metrics.json pointed at the local checkpoint, the RIS path kept as `checkpoint_ris`). The local `jobs/vcf.sh` reproduces RIS's PASS / best F1, TP, FP and unplaced counts of all six exactly. Jobs: `tmp/reviewer_cpu_20261019/E3b_colo_d9` (446053) and `tmp/reviewer_cpu_20261019/E9a_colo_d9` (446057). Which v1.1 run COLO829T Illumina should use is still open (RIS is training the 1-GPU recipe with seeds 1 and 2). Class i (GRCh38 allele) 4.3 / 3.2 / 2.2 % (Illumina / Fiber-seq / ONT), no class i or class ii true positive; projection moves complete-output F1 (no BED, after the PoN) by 0.0006-0.0037. `collect.py` also re-read the GRCh38-only INDEL reports, which were re-scored with the 2026-10-16 INDEL PoN rule after the 10-15 tables were written.

Reviewer points: **R3.4** (metrics exclude graph-coordinate calls; the "about 2 %"; CHM13), **R3.3** (validation of
graph-coordinate variants: here the placement step only), **R1.3** (are the many ONT graph-coordinate calls improved
resolution or extra ambiguity?). The design is in `../../response_to_reviewers_detailed.md`, sections R3.3, R3.4 and
R1.3, and rows E9a, E9b and E9e of the E-table (section 6). E9c (graph-native population filter) and E9d (assembly
validation) are the separate experiment `../E9c_offref_filter_validation/`. It reads the per-call TSVs written here.

## Status (2026-10-15)

- **DONE: E9a and E9e on all nine in-sample SNV runs on a pangenome graph that exist on this cluster.**
  - HPRC v2.1 d46: HG008T Illumina / PacBio / ONT and COLO829T Illumina / fiber-seq / ONT.
  - HPRC v1.1 d9: HG008T Illumina / PacBio / ONT.
  - Arrays 387435 (3 d46 runs, 2026-10-13) and 405787 (the other 6, 2026-10-15). `collect.py` rewrote `results/` at
    2026-10-15 17:52, after the last call-set output (17:50).
- **GRCh38-only (linear) arm.** That graph has no off-reference nodes, so it has no unplaced calls and nothing to
  project. `collect.py` reads the six SNV and six INDEL linear runs' own reports (Results §7).
- **E9b (CHM13): dropped by the user on 2026-10-13, restarted by the user on 2026-10-19.** The CHM13 step of
  `prep_graph.py` had failed in both graph-prep jobs (every call of the 2026-10-15 outputs has `chm13_class =
  not_built`). The CHM13 coordinates are now the pipeline's own: `graph_prep ref-path-scan --reference-sample CHM13`
  + `ref-path-check` against chm13v2.0.fa (the scan reads a small GFA from `chm13_mini_gfa.py`: the GRCh38 path's segment lengths + the CHM13 W lines of the leftover `graph/<g>/chm13.gaf`; d9 446492, d46 446491. BeeGFS read the full GFAs at ~10 MB/s on 2026-10-19, and the full d9 GFA scan 446268 stalled and was cancelled. The 109 d9 walks equal the d9 GFA's CHM13 W lines (contig, start, end, node count, min / max node; the GRCh38 scan's walks.ndjson). The first check against chm13v2.0.fa failed on both graphs; comparing every CHM13 base with the FASTAs (`tmp/reviewer_cpu_20261019/E9b_chm13/compare_chm13_paths.py`, jobs 446746, 446758, 446760) showed: d9's two chrY walks are GRCh38's chrY (0 differences over 23.9 Mb) and its other 107 walks equal chm13v2.0.fa; d46 equals chm13v2.0_maskedY_rCRS.fa (rCRS chrM) over all 93 walks. So d9 is rebuilt with `--exclude-contigs chrY` and checked against chm13v2.0.fa (447886), and d46 is checked against chm13v2.0_maskedY_rCRS.fa (447959); CHM13 VCFs 447960 and projection array 447961 follow; `tmp/reviewer_cpu_20261019/E9b_chm13/`), written
  to `<graph_index>/<graph>.chm13_path`. `<scratch>/graph/{d9,d46}/chm13_path` (empty directories) were replaced by
  symlinks to them. The twelve SNV runs are projected again into `<scratch>/runs_e9b` (array 447961), so `runs/` stays
  as it is until the GRCh38 columns are checked equal. pansoma_net_v2's CHM13 VCFs of the 24 SNV / INDEL runs (447960)
  go to `/scratch/jshen/data/reviewer_experiments_20261013/E9b_chm13_vcf/`.
- **E9b result (2026-10-19, `<scratch>/runs_e9b`, twelve SNV runs; `runs/` unchanged).** The GRCh38 columns of every
  `projected.tsv` and all GRCh38 metrics equal those of `runs/` row for row; only the `chm13_*` fields changed (they
  were `not_built`). Unplaced PASS SNV calls (no GRCh38 position), summed over the twelve runs: 14,419. Of these,
  6,860 (47.6 %; 36.7-63.1 % per run) have a direct CHM13 position (the event's node is on the CHM13 path), 5,070
  (35.2 %) lie in a CHM13 anchor interval, and 2,489 (17.3 %) have neither; 2,054 of the last are COLO829T ONT d46.
  By GRCh38 class: i 246 / 252 / 15, ii 3,519 / 4,446 / 71, iii 3,095 / 372 / 2,403 (direct / anchor interval /
  none).
- **CHM13 VCFs (pipeline step, 24 SNV / INDEL runs, `collect_chm13_vcf.py` -> `results/chm13_vcf_by_run.tsv`).** For
  each of the twelve SNV runs, the PASS calls placed on CHM13 but not on GRCh38 equal the E9b direct count exactly
  (111, 312, 1,981, 1,677, 663, 115, 57, 476, 612, 660, 85, 111). PASS calls placed on GRCh38 but not on CHM13 are
  0-101 per SNV run except COLO829T Illumina d9 (2,581). Of those 2,581, 2,554 are at GRCh38 chr1 120-145 Mb (1,820 at
  120-125 Mb, around the modelled centromere), where the graphs' CHM13 path is clipped. Most CHM13-unplaced records of
  a run are non-PASS candidates there (HG008T ONT d9: 47,096 of 51,625 at 120-145 Mb; none of its 809 CHM13-unplaced
  PASS calls is).
- **NOT PROCESSED: INDEL runs, the `__norescue` re-calls, the COLO829T fiber-seq joint model and the cross-sample runs.**
  They are finished and `run_all.sh` would queue them, but that needs a new Slurm array (see "Still missing").
- **PENDING: the final complete-output number, which needs E9c / E9d.** The numbers here count every call that is not
  placed as a GRCh38 allele as a false positive.
- No E9a job is queued or held.

## Question

`linear_vcf` writes every graph-VCF record whose event has no unique GRCh38 position to
`<set>.<KIND>.linear.unplaced.vcf.gz`. These records are not scored by rtg. For SNVs most of them are calls of the
off-reference rescue (ch6 of the site set to the every-haplotype code). In the nine SNV runs here, all of them are. In
each run the number of unplaced PASS calls (115 to 4,569; Results §1) equals `offref_called` in the run's
`*.graph.stats.json`. The questions:

1. How many of the unplaced PASS calls can be written as a GRCh38 variant after all? The bubble's flanking GRCh38
   nodes give a frame (E9a).
2. How many have a T2T-CHM13 position, and how many are on neither reference, i.e. specific to the graph (E9b)?
   **Dropped** (user decision, 2026-10-13). Not answered here.
3. What are P, R and F1 of the complete output when the placed calls are scored together with the unplaced ones, with
   and without the PoN, next to the published placed-only numbers (E9e)?

## Data

**Runs.** These are chr1 re-predictions of the RIS 7fc7396 checkpoints, made on this cluster. The RIS run directories
copied here (`pansoma_net_v2_runs/ris_7fc7396_20261012/<RUN>/`) hold the checkpoint, the training logs and
`test_chr1/*.metrics.json`, but no prediction or VCF files. Each run directory is
`/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<RUN>/vcf_chr1`. Its VCFs come from
`machine_learning/pansoma_net_v2/jobs/vcf.sh` (graph_vcf -> linear_vcf -> PoN -> rtg vcfeval).

| run | set | graph | threshold (validation recall 0.9) | prediction | VCF + eval |
|---|---|---|---:|---|---|
| `HG008_Illumina_SNV_7fc7396_ddp4_b1024` | HG008T_Illumina | HPRC v1.1 d9 | 0.25127 | tequila GPU array 405185 (`tasks2.tsv`) | same task |
| `HG008_PacBio_SNV_7fc7396` | HG008T_PacBio | HPRC v1.1 d9 | 0.3142 | 405185 | same task |
| `HG008_ONT_SNV_7fc7396` | HG008T_ONT | HPRC v1.1 d9 | 0.44923 | 405185 | same task |
| `HG008_Illumina_SNV_d46_7fc7396` | HG008T_Illumina_d46 | HPRC v2.1 d46 | 0.32375 | 405185 | same task |
| `HG008_PacBio_SNV_d46_7fc7396` | HG008T_PacBio_d46 | HPRC v2.1 d46 | 0.29122 | job 387135 (CPU, fp32) | 387323_0 (job 387359) |
| `HG008_ONT_SNV_d46_7fc7396` | HG008T_ONT_d46 | HPRC v2.1 d46 | 0.23086 | 387379_1 (job 387380; tequila GPU, bf16 as on RIS) | same task |
| `COLO829_Illumina_SNV_d46_7fc7396` | COLO829T_Illumina_d46 | HPRC v2.1 d46 | 0.6395 | 405185 | same task |
| `COLO829_fiberseq_SNV_d46_7fc7396` | COLO829T_fiberseq_d46 | HPRC v2.1 d46 | 0.48119 | 405185 | same task |
| `COLO829_ONT_SNV_d46_7fc7396` | COLO829T_ONT_d46 | HPRC v2.1 d46 | 0.56769 | 387379_2 (job 387381; tequila GPU) | same task |

The d9 checkpoints of HG008T PacBio and ONT SNV were trained on 1 GPU (global batch 1,024); the d46 and linear ones and
HG008T Illumina d9 (`ddp4_b1024`) on 4 GPUs (E3 caveat 1). COLO829T has no single-set d9 run on this cluster.

**Development runs (earlier d9 / HPRC v1.1 models, not the paper's 7fc7396 checkpoints).** These used the runs'
existing `vcf_chr1`:

- `pansoma_net_v2_runs/HG008_ONT_SNV_pow075_S_lr1e4_e20` (SNV);
- `pansoma_net_v2_runs/HG008_Illumina_INDEL_small_nofilt` (INDEL).

**Truth and BED.**

- HG008T: GIAB draft v0.2 somatic small variants (`HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz`). chr1
  has 702 SNV truth records, 697 of them in `HG008-T_somatic_smvar_benchmark_v0.2_all.bed` (TP + FN in
  `metrics.json`).
- COLO829T: `pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz`, 2,443 chr1 SNV records. The BED is
  the SMaHT easy + difficult + extreme union, which covers chr1, so for COLO829T the BED and no-BED evaluations are
  identical.

**PoN.** The four PoN VCFs (`af-only-gnomad.hg38`, `dbsnp138`, `1000g_pon.hg38`, `CoLoRSdb v1.1.0`) with the run's own
rules (`scripts/filter_panel_of_normals.py`, called by `pansoma_net_v2.vcfeval`).

**Graphs.**

- d9: `AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.{gfa,gbz}`.
- d46: `AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.{gfa,gbz}`.
- Both: the GRCh38 reference path of each graph index (`*.grch38_path`) and the graph-index sqlite (node sequences).

**Published reference numbers.** `analysis/benchmark_7fc7396_20261012/*/summary.json`: the RIS GPU predictions of the
same checkpoints. `collect.py` copies them into `results/summary.json`.

## Method

1. **Graph prep** (`prep_graph.py`, job script `tmp/reviewer_cpu_20261013/E9a_unplaced_projection/prep_graph.sh`).
   - Builds a CSR of oriented node-side successors from the GFA's L lines. Only the L-line byte range is read: d9 bytes
     47,712,856,890-49,579,776,396; d46 bytes 4,316,535,066-7,216,737,213.
   - Check: every step of the GRCh38 path must be an edge. d9: 49,092,335 steps, 0 without an edge. d46: 83,303,489
     steps, 0 without an edge.
   - Sizes: d9 has 60,118,570 nodes, 70,105,101 L lines and 140,210,202 directed side edges. d46 has 94,468,093 nodes,
     103,555,446 L lines and 207,110,892 directed side edges (`graph/<g>/edges/meta.json`).
   - The CHM13 coordinates (`vg paths -A -S CHM13` then `chm13_path/`) were not built. The step stopped at
     `unexpected GAF name b'@HD'`.
2. **Projection** (`project_unplaced.py`; the full definition is in its docstring). This runs per unplaced record.
   - **Anchors.** A Dijkstra search in bp from both ends of the event's node path. It stops at the first node with
     exactly one GRCh38 visit. Limits: 10,000 bp, 5,000 nodes, a GRCh38 span of 100,000 and 4 x 10^7 alignment cells.
   - **Bubble.** The consistent left + right anchor pair with the smallest span is the bubble. The branch through the
     event node is aligned globally to the GRCh38 interval between the anchors. Scoring: match 2, mismatch -4, gap 4
     plus 2 per base; 32 flanking bases on each side.
   - **Class i (written as a GRCh38 allele).**
     - SNV: the base is aligned to a GRCh38 base, the ALT differs from it, and there is no gap within 10 bases. The
       detail is `exact`, or `branch_base_differs` when the branch base itself differs from GRCh38.
     - INDEL: the event is moved to GRCh38 and left-aligned. The detail is `exact` or `composite`.
   - **Class ii (an interval only).** The detail is one of `in_insertion`, `near_branch_indel` (within 10 bases of a
     branch INDEL), `alt_is_grch38_base` (the read shows GRCh38's own base on a branch node, i.e. no GRCh38 variant) or
     `too_long`.
   - **Class iii (no consistent anchor pair within the limits).** The detail is `no_left_anchor`, `no_right_anchor` or
     `inconsistent_anchors`.
   - **Merging.** Class-i records are added to `linear.vcf.gz` to make `complete.linear.vcf.gz`. When a class-i record
     has the same (CHROM, POS, REF, ALT) as a placed record, one is kept (PASS first, then higher p_somatic). This
     happened 11 / 11 / 7 times.
3. **Complete-output evaluation** (`evaluate_complete.py`).
   - rtg vcfeval of `complete.linear.vcf.gz` through the repository's `pansoma_net_v2.vcfeval`, unmodified. It uses the
     run's truth, BED, SDF, the four PoNs and the run's threshold.
   - **Class ii.** A PASS class-ii call is a TP when a truth record of the kind that rtg left as FN lies within its
     interval (+-1 bp). Each truth record is claimed once, and calls are taken by decreasing p_somatic. Every other
     class-ii call is an FP, and `alt_is_grch38_base` calls are always FPs. The PoN cannot tag class ii or iii (no
     GRCh38 allele).
   - **Class iii.** Counted as FPs (worst case) or dropped.
   - **Scoring variants.**
     - `placed`: the run's own eval.
     - `placed_unplaced_fp`: every unplaced PASS call counts as an FP.
     - `projected_rtg`: placed + class i, with class ii and iii dropped.
     - `projected_ii`: adds class ii scored by interval.
     - `projected_worst`: adds class iii as FPs.
4. **`collect.py`** gathers every `runs/<RUN>/<VCFDIR>/{projection,metrics}.json` into `results/` and writes
   `<scratch>/<RUN>.projected.tsv`, the input of E9c.

## Scripts and jobs

Scripts are in this folder. Logs are in `tmp/reviewer_cpu_20261013/E9a_unplaced_projection/logs/`. Data is in
`/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection/` (below: `<scratch>`).

| step | script | job(s) | state | output |
|---|---|---|---|---|
| graph prep d9 / d46 (1st try) | `prep_graph.sh` -> `prep_graph.py` | 387144 / 387145 | cancelled after 7 min | - |
| graph prep d9 / d46 | same | 387183 (meow, 1:18 h, 17.6 GB) / 387184 (miao, 0:54 h, 32.4 GB) | edges done and checked; FAILED in the CHM13 step | `<scratch>/graph/{d9,d46}/edges/` |
| dev VCF of HG008T PacBio d46 | (inline job) | 387355 | completed | `<scratch>/selfvcf/` (superseded, see Caveats) |
| dev: two d9 runs (`tasks.dev1.txt`, `E9A_OUT=<scratch>/dev`) | `task.sh` | 387367_0 / 387367_1 | _0 completed; _1 metrics written, then the final `mv` of the coordination TSV failed | `<scratch>/dev/<RUN>/vcf_chr1/` |
| E9a + E9e, three d46 runs (`tasks.d46_insample_20261013.txt`, %1) | `task.sh` -> `project_unplaced.py`, `evaluate_complete.py` | 387435_0 (PacBio), _1 (HG008T ONT), _2 (COLO829T ONT) | completed (4:35, 5:56, 5:10) | `<scratch>/runs/<RUN>/vcf_chr1/`, `<scratch>/<RUN>.projected.tsv` |
| E9a + E9e, the other six SNV runs (`tasks.snv_rest_20261015.txt`: HG008T Illumina / PacBio / ONT d9, HG008T Illumina d46, COLO829T Illumina / fiber-seq d46) | same | 405787_0-_5 | completed 2026-10-15 15:08-15:33 | same |
| tables | `collect.py` (login node) | - | last run 2026-10-15 17:52 | `results/{tables.md,projection_by_run.tsv,metrics_by_run.tsv,summary.json}` |

Per run of the first three, projection took 79 / 223 / 195 s and evaluation 195 / 132 / 115 s, on 4 CPUs with at most
1.1 GB RSS (PacBio / HG008T ONT / COLO829T ONT).

`run_all.sh [--dry-run]` builds the task list of every finished run under `reviewer_cpu_20261013` that has no
`metrics.json` yet, and submits it as an array of at most 6 tasks at a time.

## Results (2026-10-15; 9 runs: chr1, SNV, HPRC v1.1 d9 and v2.1 d46)

Rows are ordered by tumor, platform and graph. The GRCh38-only arm has no unplaced calls; its own numbers are in §7.

### 1. chr1 unplaced PASS calls and their GRCh38 classes (from `results/tables.md`)

All filters: unplaced and placed records per run are in `results/projection_by_run.tsv` (for example 7,752 unplaced
next to 20,581 placed for COLO829T ONT d46). The PASS calls are below. The CHM13 columns of `tables.md` are left out
because they were not built.

| run | set | graph | unplaced PASS | i (GRCh38 allele) | ii (interval) | iii (no anchor) | chr1 unplaced share of PASS, raw: before -> after placement | after the PoN: before -> after placement |
|---|---|---|---:|---:|---:|---:|---|---|
| `COLO829_Illumina_SNV_d46_7fc7396` | COLO829T_Illumina_d46 | d46 | 516 | 39 (7.6%) | 250 (48.4%) | 227 (44.0%) | 0.157 -> 0.145 | 0.204 -> 0.190 |
| `COLO829_fiberseq_SNV_d46_7fc7396` | COLO829T_fiberseq_d46 | d46 | 257 | 12 (4.7%) | 203 (79.0%) | 42 (16.3%) | 0.086 -> 0.082 | 0.114 -> 0.109 |
| `COLO829_ONT_SNV_d46_7fc7396` | COLO829T_ONT_d46 | d46 | 4,569 | 128 (2.8%) | 1,334 (29.2%) | 3,107 (68.0%) | 0.597 -> 0.580 | 0.680 -> 0.670 |
| `HG008_Illumina_SNV_7fc7396_ddp4_b1024` | HG008T_Illumina | d9 | 115 | 22 (19.1%) | 89 (77.4%) | 4 (3.5%) | 0.019 -> 0.016 | 0.080 -> 0.066 |
| `HG008_Illumina_SNV_d46_7fc7396` | HG008T_Illumina_d46 | d46 | 754 | 58 (7.7%) | 440 (58.4%) | 256 (34.0%) | 0.174 -> 0.161 | 0.345 -> 0.326 |
| `HG008_PacBio_SNV_7fc7396` | HG008T_PacBio | d9 | 210 | 14 (6.7%) | 191 (91.0%) | 5 (2.4%) | 0.036 -> 0.034 | 0.150 -> 0.141 |
| `HG008_PacBio_SNV_d46_7fc7396` | HG008T_PacBio_d46 | d46 | 240 | 19 (7.9%) | 221 (92.1%) | 0 (0.0%) | 0.069 -> 0.063 | 0.160 -> 0.149 |
| `HG008_ONT_SNV_7fc7396` | HG008T_ONT | d9 | 1,410 | 47 (3.3%) | 1,238 (87.8%) | 125 (8.9%) | 0.191 -> 0.185 | 0.525 -> 0.515 |
| `HG008_ONT_SNV_d46_7fc7396` | HG008T_ONT_d46 | d46 | 1,258 | 43 (3.4%) | 984 (78.2%) | 231 (18.4%) | 0.250 -> 0.241 | 0.466 -> 0.457 |

The unplaced share is computed over the whole of chr1 (no BED): unplaced PASS / (placed PASS + unplaced PASS).

### 2. Details of the chr1 PASS classes (counts; from `results/tables.md`)

| run | i exact | i branch base differs | ii in insertion | ii ALT = GRCh38 base | ii near branch INDEL | ii too long | iii no left anchor | iii no right anchor | class-i PASS tagged by the PoN | class i merged with a placed call |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| COLO829T Illumina d46 | 21 | 18 | 98 | 125 | 12 | 15 | 213 | 14 | 16 of 39 | 3 |
| COLO829T fiber-seq d46 | 8 | 4 | 88 | 86 | 28 | 1 | 30 | 12 | 9 of 12 | 9 |
| COLO829T ONT d46 | 94 | 34 | 589 | 608 | 116 | 21 | 2,836 | 271 | 88 of 128 | 11 |
| HG008T Illumina d9 | 4 | 18 | 39 | 43 | 7 | 0 | 3 | 1 | 20 of 22 | 56 |
| HG008T Illumina d46 | 22 | 36 | 235 | 155 | 42 | 8 | 249 | 7 | 48 of 58 | 32 |
| HG008T PacBio d9 | 4 | 10 | 76 | 100 | 15 | 0 | 5 | 0 | 12 of 14 | 33 |
| HG008T PacBio d46 | 2 | 17 | 108 | 89 | 24 | 0 | 0 | 0 | 18 of 19 | 7 |
| HG008T ONT d9 | 26 | 21 | 336 | 849 | 53 | 0 | 86 | 39 | 40 of 47 | 37 |
| HG008T ONT d46 | 19 | 24 | 383 | 517 | 84 | 0 | 222 | 9 | 38 of 43 | 11 |

`i composite` and `iii inconsistent anchors` are 0 everywhere (SNV runs).

### 3. chr1 PASS P / R / F1 in the BED, after the PoN: placed-only vs complete output (from `results/tables.md`)

| run | published (RIS) | placed (local re-prediction) | placed + all unplaced as FP | + class i (rtg) | + class ii by interval | + class iii as FP (worst) | class ii TP / calls | class iii |
|---|---|---|---|---|---|---|---|---:|
| COLO829T Illumina d46 | 0.867 / 0.715 / **0.783** | 0.867 / 0.714 / **0.783** | 0.690 / 0.714 / **0.702** | 0.858 / 0.715 / **0.780** | 0.765 / 0.715 / **0.739** | 0.695 / 0.715 / **0.705** | 1 / 250 | 227 |
| COLO829T fiber-seq d46 | 0.866 / 0.707 / **0.779** | 0.866 / 0.707 / **0.779** | 0.767 / 0.707 / **0.736** | 0.864 / 0.707 / **0.778** | 0.785 / 0.707 / **0.744** | 0.770 / 0.707 / **0.737** | 0 / 203 | 42 |
| COLO829T ONT d46 | 0.788 / 0.692 / **0.737** | 0.787 / 0.692 / **0.736** | 0.252 / 0.692 / **0.369** | 0.773 / 0.692 / **0.730** | 0.480 / 0.692 / **0.567** | 0.255 / 0.692 / **0.373** | 0 / 1334 | 3107 |
| HG008T Illumina d9 | 0.452 / 0.831 / **0.585** | 0.452 / 0.831 / **0.585** | 0.415 / 0.831 / **0.553** | 0.452 / 0.832 / **0.586** | 0.425 / 0.832 / **0.562** | 0.423 / 0.832 / **0.561** | 0 / 83 | 4 |
| HG008T Illumina d46 | 0.429 / 0.852 / **0.570** | 0.429 / 0.852 / **0.570** | 0.278 / 0.852 / **0.419** | 0.428 / 0.855 / **0.571** | 0.337 / 0.855 / **0.484** | 0.295 / 0.855 / **0.438** | 0 / 375 | 256 |
| HG008T PacBio d9 | 0.497 / 0.828 / **0.621** | 0.497 / 0.828 / **0.621** | 0.421 / 0.828 / **0.558** | 0.497 / 0.829 / **0.622** | 0.431 / 0.829 / **0.567** | 0.429 / 0.829 / **0.566** | 0 / 179 | 5 |
| HG008T PacBio d46 | 0.473 / 0.839 / **0.605** | 0.472 / 0.835 / **0.603** | 0.395 / 0.835 / **0.537** | 0.472 / 0.835 / **0.603** | 0.404 / 0.835 / **0.544** | 0.404 / 0.835 / **0.544** | 0 / 208 | 0 |
| HG008T ONT d9 | 0.466 / 0.826 / **0.596** | 0.466 / 0.826 / **0.596** | 0.218 / 0.826 / **0.345** | 0.465 / 0.828 / **0.595** | 0.239 / 0.828 / **0.370** | 0.227 / 0.828 / **0.356** | 0 / 1177 | 125 |
| HG008T ONT d46 | 0.432 / 0.859 / **0.575** | 0.432 / 0.859 / **0.575** | 0.226 / 0.859 / **0.358** | 0.430 / 0.859 / **0.573** | 0.257 / 0.859 / **0.396** | 0.234 / 0.859 / **0.367** | 0 / 940 | 231 |

The same comparison as chr1 TP / FP / FN (`results/metrics_by_run.tsv`, `bed_pon`):

| run | placed | + all unplaced as FP | + class i (rtg) | + class ii by interval | + class iii as FP |
|---|---|---|---|---|---|
| COLO829T Illumina d46 | 1745 / 268 / 698 | 1745 / 784 / 698 | 1746 / 289 / 697 | 1747 / 538 / 696 | 1747 / 765 / 696 |
| COLO829T fiber-seq d46 | 1728 / 268 / 715 | 1728 / 525 / 715 | 1728 / 271 / 715 | 1728 / 474 / 715 | 1728 / 516 / 715 |
| COLO829T ONT d46 | 1690 / 457 / 753 | 1690 / 5026 / 753 | 1690 / 496 / 753 | 1690 / 1830 / 753 | 1690 / 4937 / 753 |
| HG008T Illumina d9 | 579 / 702 / 118 | 579 / 817 / 118 | 580 / 703 / 117 | 580 / 786 / 117 | 580 / 790 / 117 |
| HG008T Illumina d46 | 594 / 792 / 103 | 594 / 1546 / 103 | 596 / 795 / 101 | 596 / 1170 / 101 | 596 / 1426 / 101 |
| HG008T PacBio d9 | 577 / 583 / 120 | 577 / 793 / 120 | 578 / 584 / 119 | 578 / 763 / 119 | 578 / 768 / 119 |
| HG008T PacBio d46 | 582 / 650 / 115 | 582 / 890 / 115 | 582 / 651 / 115 | 582 / 859 / 115 | 582 / 859 / 115 |
| HG008T ONT d9 | 576 / 660 / 121 | 576 / 2070 / 121 | 577 / 664 / 120 | 577 / 1841 / 120 | 577 / 1966 / 120 |
| HG008T ONT d46 | 599 / 788 / 98 | 599 / 2046 / 98 | 599 / 793 / 98 | 599 / 1733 / 98 | 599 / 1964 / 98 |

The class-i PASS calls under rtg give these chr1 TP / FP counts in the BED (`<scratch>/runs/<RUN>/vcf_chr1/metrics.json`,
`class_i_rtg`); the TP alleles are from the runs' `eval_bed/vcfeval/*/tp.vcf.gz`:

| run | before the PoN | after the PoN | class-i TP after the PoN |
|---|---|---|---|
| COLO829T Illumina d46 | 1 / 37 | 1 / 21 | chr1:238,211,589 C>A |
| COLO829T fiber-seq d46 | 0 / 12 | 0 / 3 | - |
| COLO829T ONT d46 | 0 / 127 | 0 / 39 | - |
| HG008T Illumina d9 | 1 / 21 | 1 / 1 | chr1:78,932,081 A>T |
| HG008T Illumina d46 | 3 / 44 | 2 / 3 | chr1:78,932,081 A>T, chr1:234,995,813 G>T |
| HG008T PacBio d9 | 1 / 13 | 1 / 1 | chr1:78,932,081 A>T |
| HG008T PacBio d46 | 1 / 18 | 0 / 1 | - |
| HG008T ONT d9 | 1 / 40 | 1 / 4 | chr1:78,932,081 A>T |
| HG008T ONT d46 | 1 / 40 | 0 / 5 | - |

All class-i TPs are `branch_base_differs`. The raw class-i TP of the three HG008T d46 runs includes chr1:183,643,262 T>A,
which the PoN tags (`PanelOfNormals` in `eval_bed/complete.linear.pon.vcf.gz`).

### 4. chr1 PASS P / R / F1 in the BED, before the PoN (from `results/tables.md`)

| run | published (RIS) | placed (local re-prediction) | placed + all unplaced as FP | + class i (rtg) | + class ii by interval | + class iii as FP (worst) |
|---|---|---|---|---|---|---|
| COLO829T Illumina d46 | 0.644 / 0.731 / **0.685** | 0.644 / 0.730 / **0.684** | 0.543 / 0.730 / **0.623** | 0.636 / 0.731 / **0.680** | 0.584 / 0.731 / **0.649** | 0.544 / 0.731 / **0.624** |
| COLO829T fiber-seq d46 | 0.647 / 0.724 / **0.683** | 0.646 / 0.724 / **0.683** | 0.591 / 0.724 / **0.650** | 0.643 / 0.724 / **0.681** | 0.599 / 0.724 / **0.655** | 0.591 / 0.724 / **0.650** |
| COLO829T ONT d46 | 0.561 / 0.708 / **0.626** | 0.561 / 0.708 / **0.626** | 0.226 / 0.708 / **0.343** | 0.539 / 0.708 / **0.612** | 0.381 / 0.708 / **0.495** | 0.226 / 0.708 / **0.343** |
| HG008T Illumina d9 | 0.106 / 0.862 / **0.188** | 0.106 / 0.862 / **0.188** | 0.104 / 0.862 / **0.185** | 0.106 / 0.864 / **0.188** | 0.104 / 0.864 / **0.186** | 0.104 / 0.864 / **0.186** |
| HG008T Illumina d46 | 0.177 / 0.887 / **0.295** | 0.177 / 0.887 / **0.295** | 0.146 / 0.887 / **0.250** | 0.176 / 0.891 / **0.293** | 0.159 / 0.891 / **0.270** | 0.149 / 0.891 / **0.255** |
| HG008T PacBio d9 | 0.108 / 0.857 / **0.191** | 0.108 / 0.857 / **0.192** | 0.104 / 0.857 / **0.185** | 0.108 / 0.858 / **0.191** | 0.104 / 0.858 / **0.186** | 0.104 / 0.858 / **0.186** |
| HG008T PacBio d46 | 0.189 / 0.875 / **0.310** | 0.190 / 0.871 / **0.312** | 0.177 / 0.871 / **0.294** | 0.189 / 0.872 / **0.311** | 0.178 / 0.872 / **0.295** | 0.178 / 0.872 / **0.295** |
| HG008T ONT d9 | 0.102 / 0.855 / **0.182** | 0.102 / 0.855 / **0.182** | 0.082 / 0.855 / **0.150** | 0.101 / 0.857 / **0.181** | 0.085 / 0.857 / **0.154** | 0.083 / 0.857 / **0.151** |
| HG008T ONT d46 | 0.171 / 0.900 / **0.287** | 0.171 / 0.900 / **0.287** | 0.127 / 0.900 / **0.223** | 0.169 / 0.901 / **0.285** | 0.135 / 0.901 / **0.235** | 0.129 / 0.901 / **0.225** |

### 5. chr1 PASS P / R / F1 on the whole of chr1 (no BED), after the PoN (from `results/tables.md`)

For COLO829T the BED covers chr1, so these rows equal §3.

| run | published (RIS) | placed (local re-prediction) | placed + all unplaced as FP | + class i (rtg) | + class ii by interval | + class iii as FP (worst) | class ii TP / calls |
|---|---|---|---|---|---|---|---|
| HG008T Illumina d9 | 0.442 / 0.831 / **0.577** | 0.443 / 0.830 / **0.578** | 0.407 / 0.830 / **0.546** | 0.443 / 0.832 / **0.578** | 0.415 / 0.832 / **0.554** | 0.414 / 0.832 / **0.553** | 0 / 89 |
| HG008T Illumina d46 | 0.418 / 0.852 / **0.561** | 0.418 / 0.852 / **0.561** | 0.274 / 0.852 / **0.414** | 0.416 / 0.855 / **0.560** | 0.319 / 0.855 / **0.465** | 0.281 / 0.855 / **0.423** | 0 / 440 |
| HG008T PacBio d9 | 0.488 / 0.828 / **0.614** | 0.489 / 0.828 / **0.614** | 0.415 / 0.828 / **0.553** | 0.489 / 0.829 / **0.615** | 0.421 / 0.829 / **0.559** | 0.420 / 0.829 / **0.557** | 0 / 191 |
| HG008T PacBio d46 | 0.466 / 0.839 / **0.600** | 0.466 / 0.835 / **0.598** | 0.391 / 0.835 / **0.533** | 0.465 / 0.835 / **0.598** | 0.396 / 0.835 / **0.537** | 0.396 / 0.835 / **0.537** | 0 / 221 |
| HG008T ONT d9 | 0.454 / 0.826 / **0.586** | 0.454 / 0.826 / **0.586** | 0.216 / 0.826 / **0.342** | 0.453 / 0.828 / **0.585** | 0.230 / 0.828 / **0.361** | 0.220 / 0.828 / **0.347** | 0 / 1238 |
| HG008T ONT d46 | 0.419 / 0.859 / **0.563** | 0.419 / 0.859 / **0.563** | 0.224 / 0.859 / **0.355** | 0.418 / 0.859 / **0.562** | 0.248 / 0.859 / **0.385** | 0.227 / 0.859 / **0.359** | 0 / 984 |

### 6. Development runs on d9 (HPRC v1.1), earlier models (from `<scratch>/dev/<RUN>/vcf_chr1/{projection,metrics}.json`)

These models are not the paper's 7fc7396 checkpoints. The runs show the method on d9 and on INDELs.

| run | kind | unplaced PASS | i | ii | iii | details (PASS) | class-i PASS PoN-tagged | unplaced share of PASS raw / after the PoN, before -> after placement |
|---|---|---:|---:|---:|---:|---|---:|---|
| `HG008_ONT_SNV_pow075_S_lr1e4_e20` | SNV | 855 | 24 | 761 | 70 | i: 8 exact, 16 branch base differs; ii: 168 in insertion, 562 ALT = GRCh38 base, 31 near branch INDEL; iii: 53 no left, 17 no right | 22 | 0.125 -> 0.121 / 0.399 -> 0.392 |
| `HG008_Illumina_INDEL_small_nofilt` | INDEL | 534 | 49 | 485 | 0 | i: 3 exact, 46 composite; ii: 414 in insertion, 71 near branch INDEL | 14 | 0.122 -> 0.111 / 0.462 -> 0.427 |

| run | eval | placed | + all unplaced as FP | + class i (rtg) | + class ii by interval | + class iii as FP | class ii TP / calls |
|---|---|---|---|---|---|---|---|
| ONT SNV | BED, after the PoN | 0.467 / 0.835 / **0.599** | 0.277 / 0.835 / **0.416** | 0.467 / 0.836 / **0.599** | 0.293 / 0.836 / **0.434** | 0.283 / 0.836 / **0.423** | 0 / 742 |
| ONT SNV | BED, before the PoN | 0.103 / 0.865 / **0.183** | 0.090 / 0.865 / **0.162** | 0.102 / 0.867 / **0.183** | 0.091 / 0.867 / **0.165** | 0.090 / 0.867 / **0.163** | 0 / 742 |
| Illumina INDEL | BED, after the PoN | 0.245 / 0.233 / **0.239** | 0.125 / 0.233 / **0.163** | 0.233 / 0.233 / **0.233** | 0.163 / 0.283 / **0.207** | 0.163 / 0.283 / **0.207** | 29 / 428 |
| Illumina INDEL | BED, before the PoN | 0.059 / 0.353 / **0.100** | 0.051 / 0.353 / **0.089** | 0.058 / 0.353 / **0.099** | 0.058 / 0.395 / **0.101** | 0.058 / 0.395 / **0.101** | 25 / 428 |

In the d9 ONT SNV run, one class-i call is a TP even after the PoN (class-i rtg TP / FP 1 / 1). The INDEL run's class-i
calls give 0 TP before and after the PoN.

### 7. GRCh38-only (linear) arm: the runs' own chr1 reports, in the BED (from `results/tables.md`)

The linear graph has no off-reference nodes, so every PASS call is placed (unplaced PASS = 0) and the placed numbers
are the complete output. Standard threshold (validation recall 0.9); INDEL also at valbest.

| tumor | platform | SNV raw P / R / F1 | SNV after the PoN | INDEL after the PoN | INDEL valbest, after the PoN |
|---|---|---|---|---|---|
| COLO829T | Illumina | 0.039 / 0.599 / 0.073 | 0.391 / 0.583 / **0.468** | 0.017 / 0.221 / 0.032 | 0.667 / 0.015 / 0.030 |
| COLO829T | fiber-seq | 0.040 / 0.642 / 0.075 | 0.745 / 0.624 / **0.679** | 0.021 / 0.298 / 0.040 | - / 0.000 / - |
| COLO829T | ONT | 0.035 / 0.557 / 0.066 | 0.340 / 0.542 / **0.418** | 0.002 / 0.282 / 0.004 | 0.000 / 0.000 / 0.000 |
| HG008T | Illumina | 0.005 / 0.836 / 0.009 | 0.335 / 0.786 / **0.470** | 0.089 / 0.262 / 0.133 | 0.426 / 0.044 / 0.080 |
| HG008T | PacBio | 0.005 / 0.882 / 0.009 | 0.373 / 0.826 / **0.514** | 0.082 / 0.261 / 0.125 | 0.423 / 0.019 / 0.036 |
| HG008T | ONT | 0.004 / 0.892 / 0.008 | 0.394 / 0.842 / **0.537** | 0.012 / 0.238 / 0.022 | 0.182 / 0.027 / 0.047 |

chr1 SNV F1 after the PoN in the BED, linear against the pangenome arms' complete output (§3; d9 = HG008T only):

| tumor / platform | linear | d9 placed / all unplaced as FP / class ii by interval | d46 placed / all unplaced as FP / class ii by interval |
|---|---|---|---|
| COLO829T Illumina | 0.468 | - | 0.783 / 0.702 / 0.739 |
| COLO829T fiber-seq | 0.679 | - | 0.779 / 0.736 / 0.744 |
| COLO829T ONT | 0.418 | - | 0.736 / 0.369 / 0.567 |
| HG008T Illumina | 0.470 | 0.585 / 0.553 / 0.562 | 0.570 / 0.419 / 0.484 |
| HG008T PacBio | 0.514 | 0.621 / 0.558 / 0.567 | 0.603 / 0.537 / 0.544 |
| HG008T ONT | 0.537 | 0.596 / 0.345 / 0.370 | 0.575 / 0.358 / 0.396 |

## Interpretation

1. **Placement onto GRCh38 recovers few calls (R3.4, item 1).** Of the unplaced PASS SNVs on chr1, 2.8 % (COLO829T
   ONT d46) to 19.1 % (HG008T Illumina d9, 22 of 115) can be written as a GRCh38 allele through the bubble anchors.
   - d46: 2.8 / 4.7 / 7.6 % (COLO829T ONT / fiber-seq / Illumina), 3.4 / 7.7 / 7.9 % (HG008T ONT / Illumina / PacBio).
   - d9 (HG008T): 3.3 / 6.7 / 19.1 % (ONT / PacBio / Illumina).
   - This is the number for the "fraction placed" gap in the R3.3 / R3.4 response text.
2. **The placed calls are mostly population alleles; they add almost no recall.**
   - Once they have a GRCh38 allele, the PoN tags 289 of the 382 class-i PASS calls of the nine runs, from 16 of 39
     (COLO829T Illumina) to 18 of 19 (HG008T PacBio d46).
   - After the PoN, class i adds 6 TPs in 5 runs, i.e. 3 distinct truth SNVs: chr1:78,932,081 A>T in four HG008T runs
     (all three d9 runs and Illumina d46), chr1:234,995,813 G>T (HG008T Illumina d46) and chr1:238,211,589 C>A
     (COLO829T Illumina). Recall rises by at most 0.003 (HG008T Illumina d46, 0.852 -> 0.855).
   - Adding class i changes chr1 F1 after the PoN, in the BED, by -0.006 (COLO829T ONT) to +0.001.
3. **Most unplaced calls have no GRCh38 allele at all.**
   - **Class ii, interval only.** 29.2 % (COLO829T ONT) to 92.1 % (HG008T PacBio d46) of the unplaced PASS calls.
     - Calls inside sequence inserted relative to GRCh38: 39 (HG008T Illumina d9) to 589 (COLO829T ONT).
     - Calls whose ALT is GRCh38's own base on a branch node: 43 to 849 per run, i.e. 13.3 % (COLO829T ONT) to 60.2 %
       (HG008T ONT d9, 849 of 1,410) of all unplaced PASS calls. These are no variant with respect to GRCh38.
     - The interval rule matched one missed truth SNV in all nine runs (COLO829T Illumina, 1 of 250 class-ii calls).
   - **Class iii, no anchor.** 0 (HG008T PacBio d46) to 68.0 % (COLO829T ONT, 3,107); 44.0 % for COLO829T Illumina.
     - On HG008T, d46 has more of them than d9: Illumina 256 against 4, ONT 231 against 125 (PacBio 0 against 5).
     - Almost all are `no_left_anchor` (2,836 of 3,107 COLO829T ONT; 249 of 256 HG008T Illumina d46).
     - In the two ONT d46 runs they cluster in node IDs 11.85 M to 12.10 M (1,959 of 3,107 and 215 of 231 class-iii
       PASS calls in `<scratch>/<RUN>.projected.tsv`; class ii: 9 and 8). In the Illumina d46 runs only 50 of 227
       (COLO829T) and 23 of 256 (HG008T) do. The node IDs are not ordered along GRCh38, so this does not locate them.
     - The median AF of the class-iii PASS calls is 0.95 to 1.0 in the Illumina, fiber-seq and HG008T ONT d46 runs
       (germline-like) and 0.37 in COLO829T ONT.
4. **Complete output (R3.4, item 3; E9e before E9c / E9d).** Recall changes by at most +0.003 (class i). The rest of the
   change is in precision. chr1 F1 after the PoN, in the BED (§3):
   - All unplaced PASS calls counted as FPs: COLO829T d46 0.783 / 0.779 / 0.736 -> 0.702 / 0.736 / 0.369 (Illumina /
     fiber-seq / ONT); HG008T d9 0.585 / 0.621 / 0.596 -> 0.553 / 0.558 / 0.345 (Illumina / PacBio / ONT); HG008T d46
     0.570 / 0.603 / 0.575 -> 0.419 / 0.537 / 0.358.
   - Class iii dropped instead (class ii scored by interval): COLO829T 0.739 / 0.744 / 0.567; HG008T d9 0.562 / 0.567 /
     0.370; HG008T d46 0.484 / 0.544 / 0.396.
   - The loss is largest for ONT (-0.217 to -0.367 with every unplaced call false) and smallest for HG008T Illumina d9
     (-0.032) and COLO829T fiber-seq (-0.043).
   - These are upper bounds on the loss. E9c (population haplotypes) and E9d (assemblies) decide how many of the
     remaining calls are population alleles, germline or confirmed somatic.
5. **Against the GRCh38-only arm (R3.1 / R3.4; §7).** The linear arm has no unplaced calls, so its chr1 F1 after the PoN
   is its complete output: 0.468 / 0.679 / 0.418 (COLO829T Illumina / fiber-seq / ONT), 0.470 / 0.514 / 0.537 (HG008T
   Illumina / PacBio / ONT).
   - Placed only, every pangenome run is above linear (0.570-0.783).
   - With every unplaced PASS call false, the pangenome run stays above linear for COLO829T Illumina (0.702) and
     fiber-seq (0.736), HG008T PacBio (d9 0.558, d46 0.537) and HG008T Illumina d9 (0.553). It falls below linear for
     HG008T Illumina d46 (0.419 against 0.470), COLO829T ONT (0.369 against 0.418) and HG008T ONT (d9 0.345, d46 0.358
     against 0.537).
   - With class iii dropped and class ii scored by interval, COLO829T ONT (0.567) and HG008T Illumina d46 (0.484) are
     above linear again; HG008T ONT is not (0.370 / 0.396).
   - The `__norescue` re-calls (E3b) keep the placed calls and cut the unplaced PASS calls to 6-162 per run (E16
     `tables.md` §1a: e.g. HG008T ONT d46 1,258 -> 162, COLO829T ONT 4,569 -> 98). They are not projected here yet.
6. **"About 2 %" (R3.4, item 4).**
   - Before the PoN, unplaced calls are 1.9 % of the PASS SNVs on chr1 for HG008T Illumina d9, the only run near
     "about 2 %". The others: HG008T d9 PacBio 3.6 %, ONT 19.1 %; HG008T d46 PacBio 6.9 %, Illumina 17.4 %, ONT 25.0 %;
     COLO829T d46 fiber-seq 8.6 %, Illumina 15.7 %, ONT 59.7 %.
   - After the PoN they are 8.0 % (HG008T Illumina d9) to 68.0 % (COLO829T ONT), because the PoN removes placed calls
     only.
   - Placement lowers the after-PoN shares by at most 1.9 points (HG008T Illumina d46, 34.5 -> 32.6 %).
   - The d9 development SNV run (an earlier model) is 12.5 % before and 39.9 % after the PoN.
7. **ONT: resolution or ambiguity (R1.3).** From E9a alone, the many ONT graph-coordinate calls are not newly resolved
   somatic variants.
   - ONT has the most unplaced PASS calls: HG008T 1,410 (d9) / 1,258 (d46) against PacBio 210 / 240 and Illumina
     115 / 754; COLO829T 4,569 against Illumina 516 and fiber-seq 257.
   - Only 2.8-3.4 % can be placed. After the PoN, the only placed TP is the shared chr1:78,932,081 A>T (HG008T ONT d9).
   - Most of the rest are calls on GRCh38 bases or inserted sequence (class ii; HG008T ONT d9: 849 of 1,410 on GRCh38
     bases), or sit in a few anchor-less graph regions (class iii; COLO829T ONT).
   - Whether they are population alleles, germline variants or errors is the question of E9c / E9d.
8. **INDEL (d9 development model only).** The 7fc7396 INDEL runs are not projected yet (see "Still missing").
   - 49 of 534 unplaced PASS INDELs can be placed, mostly as `composite` alleles. None is a TP.
   - Most class-ii INDELs lie in inserted sequence. 25 (raw) and 29 (after the PoN) of the 428 in the BED match a
     missed truth INDEL within their interval, a lenient match (Caveat 4).
   - After the PoN, the complete-output F1 is 0.207 with class ii scored by interval and 0.163 with every unplaced call
     counted as an FP, against 0.239 placed-only.

## Caveats

1. **Coverage.** Only chr1, one training seed per model, and only the nine SNV runs of the paper's 7fc7396 models on a
   pangenome graph: HG008T on d9 and d46 (three platforms each), COLO829T on d46 (three platforms).
   - No INDEL run of those models is projected yet. COLO829T has no single-set d9 run on this cluster.
   - HG008T PacBio / ONT d9 are 1-GPU checkpoints (E3 caveat 1), so their d9-d46 differences are not a matched
     comparison; HG008T Illumina d9 (`ddp4_b1024`) is matched.
   - The two d9 development runs use earlier models.
2. **Re-predictions, not the RIS predictions.**
   - "placed" comes from local re-predictions of the RIS checkpoints. HG008T PacBio d46 was predicted on CPU in fp32,
     the other eight on tequila's GPU in bf16, so the tables' column label "placed (CPU)" really means "local".
   - They differ slightly from the published RIS evaluation. In the BED after the PoN, TP / FP / FN is 582 / 650 / 115
     against 585 / 652 / 112 (HG008T PacBio d46), 1745 / 268 / 698 against 1746 / 268 / 697 (COLO829T Illumina),
     1690 / 457 / 753 against 1690 / 454 / 753 (COLO829T ONT), and one FP fewer locally in the three HG008T d9 runs
     (702 / 583 / 660 against 703 / 584 / 661). COLO829T fiber-seq and HG008T Illumina / ONT d46 are identical.
3. **`collect.py` key bug.** The `ii INDEL context` column of `results/tables.md` (and `ii_indel_context_gap` in
   `results/projection_by_run.tsv`) reads the key `ii:indel_context_gap`. `project_unplaced.py` writes
   `ii:near_branch_indel`, so the column showed 0 (correct: 116 / 84 / 24 for COLO829T ONT / HG008T ONT / HG008T PacBio
   d46). The class-ii totals and every metric were unaffected. **Fixed 2026-10-13**: `collect.py` now reads
   `ii:near_branch_indel`; every table here has the correct counts.
4. **Class-ii interval matching is lenient.**
   - Any FN truth record of the kind within the interval +-1 bp counts, whatever its allele.
   - After the PoN, truth records whose placed call was PoN-tagged become FN and can then be claimed. This is why the
     d9 INDEL run has 29 TPs after the PoN against 25 before.
   - For the SNV runs here the rule gave 1 TP in nine runs (COLO829T Illumina, 1 of 250 class-ii calls), so the
     leniency hardly matters for them.
5. **Class iii is counted inside the BED in the worst case** although it has no position. In the `+ class ii by
   interval` variant it is dropped.
6. **No assembly evidence here.** Every unplaced call not placed as a GRCh38 allele is an FP unless the class-ii interval
   rule matches a truth record. COLO829T has no tumor assembly, so even E9d can only rule out germline there.
7. **CHM13 leftovers.** `<scratch>/graph/{d9,d46}/chm13.gaf` (437 MB / 793 MB) are leftovers of the failed CHM13 step,
   and the `chm13_path/` directories are empty. `check_chm13_sqlite.py` (a check of the CHM13 coordinates against the d9
   path-index sqlite) is not needed any more.
8. **Development leftovers in `<scratch>`.**
   - `HG008_Illumina_INDEL_small_nofilt.projected.tsv.tmp.3046608` is the coordination TSV of 387367_1, whose final
     `mv` failed with "Device or resource busy". The task's `metrics.json` is complete.
   - `dev/HG008_ONT_SNV_pow075/` is an earlier interactive test (run name `DEV_HG008_ONT_SNV_pow075`) on the same VCF as
     387367_0. It has the same projection counts and no `metrics.json`.
   - `selfvcf/` is a VCF built here from the PacBio CPU predictions (job 387355) before the run's own `vcf_chr1`
     existed. It has the same counts (3,486 PASS graph records, 240 unplaced PASS; P / R / F1 after the PoN in the BED
     0.472 / 0.835 / 0.603) and is superseded.

## Still missing (2026-10-15)

1. **E9a / E9e on the other finished d9 / d46 vcf directories (75).** `run_all.sh --dry-run` would list them:
   - 18 INDEL in-sample directories: `vcf_chr1` and `vcf_chr1_valbest` of HG008T Illumina / PacBio / ONT on d9 and d46
     and COLO829T Illumina / fiber-seq / ONT on d46;
   - 9 `__norescue` SNV re-calls (E3b; same predictions, few unplaced calls);
   - 8 `vcf_chr1_<SET>` directories of the COLO829T fiber-seq joint d9 model (SNV and INDEL, P100 / P75 / P50 / P25);
   - 40 cross-sample `xs__` directories (E12a). Whether to project them is the user's decision; otherwise submit
     `task.sh` with an explicit task list, as was done with `tasks.snv_rest_20261015.txt`.
   - Each needs a Slurm array task (CPU, `general`, `-c 4 --mem 16G`, about 5 min per directory, at most 6 at a time),
     then `python collect.py`. Not submitted here.
2. **COLO829T d9 single-set runs**: added 2026-10-19 (matched recipe, seed 0; top of this README). COLO829T Illumina v1.1
   at epoch 6 (the user's choice for R3.1 / R3.9, 2026-10-20) is projected separately on 2026-10-21 (see the 2026-10-21
   block below).
3. **Comparators** (E8: ClairS-TO, DeepSomatic on GRCh38) are not yet re-scored. They are GRCh38 callers with no
   unplaced calls, so they enter only the linear-style comparison of §7.
4. **Seed repeats** of the 7fc7396 models are pending GPU (R3.9); every number here is one seed.
5. **Final E9e.** E9c / E9d now cover the same nine SNV runs (`../E9c_offref_filter_validation/`, its §2 and the E9a
   class cross-tab of its §6). The combined number (class i placed by E9a, E9c / E9d on the rest: precision =
   TP_rtg / (TP_rtg + FP_rtg + remaining unplaced PASS - confirmed)) is not assembled yet.
6. **E9b (CHM13).** Done 2026-10-19 (restarted by the user; Status above, `runs_e9b`); E9b x E9c by category:
   `chm13_by_category.py` -> `results/chm13_by_category.tsv` (2026-10-21).
