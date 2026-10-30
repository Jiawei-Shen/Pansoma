# E8: ClairS-TO and DeepSomatic tumor-only under Pansoma's protocol (chr1)

Status (2026-10-13): **CANCELLED by the user.** 2026-10-16: the user's existing VCFs are being re-scored (section 9). The user already has ClairS-TO / DeepSomatic results, so the comparators
are not re-run. No comparator was run in E8.

What is left:

- the prepared scorer, dataset table and job scripts (section 3);
- a check of the scorer on two existing DeepSomatic VCFs (section 4);
- 5.7G of partial preparation data, which the user decides whether to delete (section 5).

The user's existing comparator VCFs can be re-scored under the Pansoma protocol with `score_comparator.py` once the
user gives their paths (section 6).

## 1. Question and reviewer points

The plan was to re-run both tumor-only callers on the same reads as Pansoma and score them with the same truth, BED,
PoN rule and rtg options on chr1. That would show whether Pansoma's chr1 results hold against them, and where each
caller loses recall.

| point | what E8 was to provide |
|---|---|
| R3.m2 comparator setup | versions, models, command lines, BAMs, PoN, PASS criteria; each caller with its own filters and with our PoN |
| R1.6 low sensitivity / R2.3 recall gap to DeepSomatic | comparators' recall before / after our PoN, curve ceilings, recall at matched precision, one protocol |
| R3.5 signatures | SBS96 of truth, comparators and Pansoma PASS + PoN SNVs |
| R3.9 statistics | comparator call sets for the paired bootstrap (E16) |
| R3.10 INDEL claims | INDEL tables and strata (HP>=7, HP4-6, STR, none; INS / DEL) for every call set |
| R1.m1 additional lines | HCC1395 and CASTLE HCC1954 / H1437 / H2009 x 3 platforms |

These experiments were to take their comparator calls from E8: E1.6 (recall funnel), E5c (comparators at absorbed
loci), E10 (strata), E16 (bootstrap) and E17a (INDEL strata). They can use the re-scored user VCFs instead (section 6).

## 2. Decision and what happened to the jobs

- 2026-10-13 from 10:23 CDT: E8 submitted (387336 / 387338 about 11:15, replacing the slow index jobs). There were
  6 preparation jobs (`jobs/jobs.tsv`: bwa, merge, index, chr1 extraction) and 55 jobs 387194-387248
  (`jobs/lanes.tsv`: the callers plus the CASTLE index job 387210) in 6 dependency lanes.
- 2026-10-13 11:46:40-11:47:25: all E8 jobs cancelled (user decision).
  - The comparator image pull 387136 was cancelled too, after 1:49:05. No ClairS-TO or DeepSomatic SIF exists in
    `/scratch/jshen/bin/`.
- Final states (`sacct`): all 61 E8 jobs are CANCELLED. 57 never started. Four had started:

| job | what | ran | how far it got |
|---|---|---|---|
| 387161 | bwa mem COLO829T Illumina FASTQs -> chr1 BAM (`bwa_chr1.sh`) | 1:23:59 | 178,145,788 reads aligned (sum of `Processed` lines in `bams/COLO829T_Illumina.bwa.log`); nothing written (sort still in memory) |
| 387162 | bwa mem HCC1395 WGS_FD_T_1 -> chr1 BAM | 1:23:59 | 128,666,924 reads aligned (`bams/HCC1395T_Illumina_FD_T_1.bwa.log`); nothing written |
| 387164 | index HCC1395 PacBio / ONT raw BAMs via symlinks (`index_bams.sh`) | 0:52:03 | cancelled earlier (11:15) as too slow; replaced by 387336 |
| 387336 | chr1 extraction of HCC1395 PacBio / ONT (`extract_chr1.sh`) | 0:27:07 | partial `bams/HCC1395T_PacBio.chr1.tmp.bam`, truncated (no BGZF EOF block) |

Job logs: `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E8_comparators/logs/`. The four logs above are
there; each ends with the slurmstepd CANCELLED line.

## 3. What was prepared (kept here)

| file | what it is |
|---|---|
| `datasets.tsv` | 27 datasets: HG008T x 3, COLO829T x 3, HCC1395T x 3, CASTLE HCC1954T / H1437T / H2009T x 3, COLO829T purity mixtures P75 / P50 / P25 x 3. Columns: dataset, priority, tensor `set` (truth / BED / FASTA source), BAM, ClairS-TO platform / model (`ss` instead of `ssrs` on the CASTLE lines), ClairS-TO extra options (`--disable_verdict` for the mixtures), DeepSomatic model type, reads and aligner |
| `score_comparator.py` | scores one comparator VCF (or ClairS-TO's snv + indel pair) for one kind on chr1 with Pansoma's protocol (below). Tested (section 4) |
| `collect.py` | builds `tables.md`, `summary.json`, `sbs96.tsv`, `indel_strata.tsv` from `score/<caller>/<dataset>/<KIND>/report.json`, next to Pansoma's 7fc7396 benchmark numbers (arms d9, d46, linear; purity joint models) |
| `tables.md`, `summary.json` | generated 11:27 before the cancellation: Pansoma rows only. The comparator rows, matched precision, SBS96, INDEL strata and Commands tables are empty; `summary.json` lists 108 pending caller / dataset / kind reports |
| `jobs/` | copies of the job scripts (`bwa_chr1.sh`, `hcc1395_illumina_merge.sh`, `index_bams.sh`, `extract_chr1.sh`, `run_clairs_to.sh`, `run_deepsomatic.sh`, `submit_lanes.sh`) and `jobs.tsv` / `lanes.tsv`. These are a record only; none completed |

Scoring protocol of `score_comparator.py`:

- It uses `machine_learning/pansoma_net_v2/vcfeval.py` as a library, unmodified.
- Truth and BED: the tensor set's `tensors/<KIND>/labels.manifest.json`, the same files as Pansoma's chr1 test.
  - CASTLE sets without a manifest fall back to `CASTLE/<line>/truth/`.
  - `--truth` / `--bed` override both.
- Call preparation:
  - `bcftools norm -m -any -f <labels FASTA>`, chr1;
  - only records of the kind, by the same SNV / INDEL rule that selects the truth;
  - RefCall / NoCall records dropped;
  - one record per (POS, REF, ALT);
  - sites-only VCF that keeps the tool's QUAL and FILTER.
- PoN: `scripts/filter_panel_of_normals.py` with the four PoN VCFs of `jobs/vcf.sh`, applied to the calls and to the
  truth. The truth run gives the PoN bound on recall.
- rtg 3.12.1: `--squash-ploidy --sample ALT,ALT -f QUAL`, chr1 as `--bed-regions`. Each call set is scored inside the
  BED (`eval_bed`) and on the whole chromosome (`eval_nobed`), with `--all-records` for the curves.
- Four call sets per run:

| set | content | reviewer setting |
|---|---|---|
| `pass` | the tool's PASS records | (1) tool PASS |
| `pass_pon` | `pass` minus the records our PoN tags | (2) tool PASS + our PoN |
| `all` | every non-RefCall record, curve by QUAL | curve without our PoN |
| `all_pon` | `all` minus PoN-tagged records, curve by QUAL | (3) all calls + our PoN, by the tool's quality |

## 4. Scorer check on two existing DeepSomatic VCFs (COLO829T fiber-seq, SNV, chr1)

These are not E8 runs. They only test the scorer on real DeepSomatic output. Outputs are in
`/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/test/`.

- `smaht_ds160_fiberseq_SNV`: SMaHT's
  `/scratch/jshen/data/COLO829T/COLO829T_fiberseq/SMHTCOLO829T-X-X-M45-B004-uwsc-SMAVCQ2SLU5S-deepsomatic_1.6.0_GRCh38_snv.vcf.gz`.
  It is genome-wide and SNV only. The header does not record the mode, which may be tumor-normal [unverified].
  - Full report: `report.txt` / `report.json`.
- `local_ds190_TO_fiberseq_SNV`: `/scratch/jshen/data/COLO829T/COLO829T_fiberseq/deeepsomatic_TO_results/COLO829T_Pacbio_deepsomatic.vcf.gz`.
  This is DeepSomatic 1.9.0 tumor-only on chr1, with a PON filter.
  - The four `eval_bed` rtg runs finished; the numbers below come from their `summary.txt`.
  - `eval_nobed` finished only `pass`. The other three nobed runs were stopped to free the login node, so there is no
    `report.json`.

The truth has 2,443 chr1 SNV records in the BED. The PoN tags 84 of them, so recall after the PoN is at most 0.966. Cells
are P / R / F1 (TP / FP / FN).

| call set | tool PASS | tool PASS + our PoN | curve, all records, no PoN: best F1 (P, R) | curve, all records + our PoN: best F1 (P, R) at QUAL | curve ceiling + PoN |
|---|---|---|---|---|---|
| SMaHT DeepSomatic 1.6.0 (mode not recorded) | 0.972 / 0.745 / 0.844 (1821 / 53 / 622) | 0.974 / 0.723 / 0.830 (1765 / 47 / 678) | 0.066 (0.034, 0.746) | 0.779 (0.836, 0.728) at 3.1 | 0.728 |
| local DeepSomatic 1.9.0 tumor-only, with PON | 0.592 / 0.564 / 0.578 (1378 / 951 / 1065) | 0.733 / 0.557 / 0.633 (1361 / 496 / 1082) | 0.236 (0.233, 0.240) | 0.792 (0.724, 0.873) at 15.2 | 0.894 |
| Pansoma d9 `COLO829_fiberseq_SNV_h90` (`tables.md`, from the 7fc7396 benchmark) | 0.445 / 0.736 / 0.555 (raw PASS) | 0.857 / 0.717 / 0.781 | - | 0.784 (0.849, 0.728) | 0.792 |

Prepared chr1 SNV records (non-RefCall, `calls.all.vcf.gz`), by FILTER:

- SMaHT 1.6.0: 1,874 PASS, 51,995 GERMLINE;
- local 1.9.0: 2,329 PASS, 3,124 PON, 209,047 GERMLINE.

Interpretation:

- The scorer runs end to end on DeepSomatic VCFs: multi-allelic splitting, RefCall removal, PoN on about 2e5 records,
  and rtg with `--all-records`.
  - The SMaHT VCF took about 7 min.
  - The 1.9.0 VCF took 20 min for the PoN (10:43 -> 11:03) plus 4-5 min per rtg run on the busy login node.
- The 1.9.0 tumor-only run loses most of its SNV recall to its own germline classification.
  - Its PASS recall is 0.564.
  - After our PoN, all its records ranked by QUAL reach recall 0.873 at precision 0.724 (best F1 0.792). The ceiling
    is 0.894.
  - A likely reason [unverified, not analysed]: COLO829T is near-pure, so clonal somatic SNVs (in LOH regions
    especially) sit at germline-like VAF.
- The SMaHT 1.6.0 VCF has far fewer GERMLINE records than the 1.9.0 tumor-only VCF, and PASS precision 0.972. That fits
  a tumor-normal run [unverified]. It should not be used as a tumor-only comparator unless its mode is confirmed.
- Before our PoN, ranking all DeepSomatic records by QUAL is useless (best F1 0.066 / 0.236), because GERMLINE records
  carry high QUAL. Setting (3) is meaningful only after the PoN.

## 5. Leftover data (user decides whether to delete)

Directory: `/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/`. It is 5.7G by `du -sh`, almost all of it
one unusable partial BAM.

| path (under the directory) | size | content | useful? |
|---|---|---|---|
| `bams/HCC1395T_PacBio.chr1.tmp.bam` | 6,046,506,688 bytes | partial chr1 extraction of the HCC1395 PacBio BAM (job 387336), truncated | no |
| `bams/HCC1395T_PacBio.bam` | symlink | -> `/scratch/jshen/data/CASTLE/HCC1395/PacBio/HCC1395.GRCh38.bam` (raw data) | no. Removing the link does not touch the raw BAM; never delete through it |
| `bams/COLO829T_Illumina.bwa.log`, `bams/HCC1395T_Illumina_FD_T_1.bwa.log` | 279,331 / 207,045 bytes | bwa stderr of 387161 / 387162 | record only |
| `bams/COLO829T_Illumina.sorttmp/`, `bams/HCC1395T_Illumina_FD_T_1.sorttmp/` | empty | samtools sort temp dirs | no |
| `calls/`, `score/` | empty | planned outputs | `score/` is where re-scored reports go (section 6) |
| `test/` | 9.7M | the scorer check of section 4 | yes: the source of the section 4 numbers |

Nothing was written next to the raw data.

## 6. How to re-score the user's existing comparator VCFs

Once the user gives the VCF paths, run one scorer call per caller, dataset and kind. ClairS-TO's `snv.vcf.gz` and
`indel.vcf.gz` can be passed together, because the scorer keeps only the records of `--kind`.

- `--set`: the `set` column of `datasets.tsv`, for example `HCC1395T_PacBio_d46`. It decides truth, BED and FASTA.
- `--output`: `score/<caller>/<dataset>/<KIND>`, where `<dataset>` is the first column of `datasets.tsv`. This is the
  path `collect.py` reads.
- Use `--caller clairs_to` or `--caller deepsomatic`; `collect.py` looks only at those two directory names.

Example (ClairS-TO, HG008T PacBio, SNV):

    cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8_comparators
    /wanglab/jshen/anaconda3/bin/python score_comparator.py \
        --vcf /path/given/by/user/snv.vcf.gz --vcf /path/given/by/user/indel.vcf.gz \
        --set HG008T_PacBio --kind SNV --caller clairs_to --threads 4 \
        --output /scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score/clairs_to/HG008T_PacBio/SNV

Repeat with `--kind INDEL`, then build the tables:

    /wanglab/jshen/anaconda3/bin/python collect.py      # tables.md, summary.json, sbs96.tsv, indel_strata.tsv

Practical notes:

- Only chr1 records are scored (`--chrom chr1`, the default), so genome-wide VCFs can be passed as they are.
- The VCFs must be on GRCh38 with `chr` names. The scorer normalizes against the GCA no-alt FASTA of the labels.
  chr1 has the same name and length in GIAB's GRCh38-GIABv3; for other references this is [unverified]. `bcftools norm`
  runs with `-c w`, so a REF mismatch is only warned about.
- The scorer refuses to overwrite an existing `eval_bed/<set>` or `eval_nobed/<set>`. Give a new `--output`, or move the
  old one aside.
- Run time (section 4, busy login node):
  - about 7 min for the SMaHT VCF;
  - 42 min (10:40 -> 11:22) for the 1.9.0 tumor-only VCF, with about 2e5 GERMLINE records on chr1, up to the end of
    the four `eval_bed` runs. The `eval_nobed` runs add roughly as much again [estimate].
  - Memory was not measured.
  - Large VCFs belong in a Slurm job once the user agrees a cluster cap.
- `collect.py` labels the two callers "ClairS-TO 0.5.1" and "DeepSomatic 1.10.0" (`CALLERS`). If the user's runs used
  other versions, edit that constant before collecting.
- `collect.py` fills the Commands table (R3.m2) from `calls/<caller>/<dataset>/command.txt`. Put each user run's
  command line (version, model, BAM) there, or the table stays empty.
- The section 4 DeepSomatic 1.9.0 tumor-only VCF has no `report.json` (`test/` only). To include it, re-score it into
  `score/deepsomatic/COLO829T_fiberseq/{SNV,INDEL}`.

## 7. Caveats for the comparison with the user's results

- Reads: E8 was meant to give the comparators the same reads as the Pansoma tensors. The user's runs may differ in reads,
  aligner and depth (for example COLO829T Illumina, where E8 planned a new bwa alignment of the tensor FASTQs, and the
  purity mixtures). Record the BAM of every scored VCF.
- Mode: confirm that each VCF is tumor-only. The SMaHT DeepSomatic 1.6.0 VCF looks tumor-normal (section 4).
- Training data:
  - Pansoma models are trained on chr2-22 of the same tumor and tested on its chr1. The comparators are general
    models. This asymmetry favours Pansoma and must be stated next to every table.
  - E8b (`../E8b_literature_facts/README.md`) found that DeepSomatic holds out chr1 in every model, including the
    tumor-only models trained on HCC1395 (rows D2, D10, D14). Its chr1 numbers on the CASTLE lines are therefore
    held-out.
  - For ClairS-TO the sources disagree: chr20 or chr1 excluded (E8b row C4). If chr20, the `ssrs` models have seen
    chr1 of HCC1954 / H1437 / H2009, and the `ss` models should be used there.
  - CASTLE truth is DeepSomatic's own benchmark, which the DeepSomatic paper calls biased toward DeepSomatic (E8b row D8).
- QUAL ranking (setting 3) uses the tool's QUAL. For DeepSomatic, GERMLINE records carry high QUAL, so only the curve
  after our PoN (`all_pon`) is meaningful (section 4).
- HG008T ONT is dorado sup v5.0.0. ClairS-TO 0.5.1 has no sup v5 model; E8 would have used
  `ont_r10_dorado_sup_5khz_ssrs` (sup v4.2.0). Check which model the user's run used.

## 8. Status

| item | status |
|---|---|
| scorer `score_comparator.py`, collector `collect.py`, `datasets.tsv` | DONE; scorer tested on two existing DeepSomatic VCFs (section 4) |
| comparator calling (ClairS-TO 0.5.1, DeepSomatic 1.10.0) on 27 datasets | CANCELLED by the user (existing results are used instead); all 61 E8 jobs and image pull 387136 cancelled, nothing called |
| re-scoring the user's ClairS-TO / DeepSomatic VCFs under the Pansoma protocol | PENDING: the user has to give the VCF paths (and the BAM / version / command of each run) |
| leftover data, 5.7G | PENDING user decision: delete or keep (section 5); only `test/` is worth keeping |
| `tables.md` comparator rows | PENDING on the re-scoring above; currently Pansoma rows only |

## 9. Re-scoring the user's existing VCFs (2026-10-16)

The user pointed to `/scratch/jshen/Pansoma_testing_results_V2/` (the comparison data of the original manuscript). The
P100 COLO829T VCFs are in the older `/scratch/jshen/Pansoma_testing_results/VCF_results/`; the old rtg results in V2
name them in their `vcfeval.log` (`-c`). All runs are tumor-only and chr1 only; provenance is from the VCF headers
(`##cmdline`, `##clairs_to_version`, `##DeepVariant_version`) and the DeepSomatic logs (model `example_info.json`,
`--reads_tumor`, PON files). One line per run is in `calls/<caller>/<dataset>/command.txt` (Commands table).

| dataset | ClairS-TO (default PoNs on) | DeepSomatic tumor-only (PON filtering on) | BAM (same sequencing as the tensors) |
|---|---|---|---|
| HG008T_Illumina | 0.4.2 `ilmn_ssrs`, `V2/HG008_GIAB/AF_HPRC/clairs-to_HG008T_WGS_Chr1/` | 1.9.0 WGS_TUMOR_ONLY, `.../deepsomatic_HG008T_WGS_Chr1/HG008T_WGS.vcf.gz` | GIAB `HG008-T_Illumina_195x_GRCh38-GIABv3.bam` (BCM 2024-03-13, the FASTQs of the GAM) |
| HG008T_PacBio | 0.4.2 `hifi_revio_ssrs`, `.../clairs-to_HG008_PacBio_Chr1/` | 1.9.0 PACBIO_TUMOR_ONLY, `.../deepsomatic_HG008T_PacBio_Chr1/HG008T_PacBio.vcf.gz` | GIAB Revio 116x |
| HG008T_ONT | 0.4.2 `ont_r10_dorado_sup_5khz_ssrs`, `.../clairs-to_HG008_ONT_Chr1/` | 1.9.0 **PACBIO_TUMOR_ONLY on the ONT BAM** (logs load `pacbio_tumor_only/example_info.json`), `.../deepsomatic_HG008T_ONT_Chr1/HG008T_ONTvcf.gz` | GIAB ONT-UL R10.4.1 dorado 0.8.1 sup, 54x |
| COLO829T_Illumina | 0.4.1 `ilmn`, `VCF_results/clairs-to_COLO829T_WGS/` | 1.9.0 WGS_TUMOR_ONLY, `VCF_results/deepsomatic_COLO829T_WGS/COLO829T_WGS.chr1.vcf.gz` | SMaHT `SMHTCOLO829T-X-X-M45-A001-uwsc-...-sentieon_bwamem` |
| COLO829T_fiberseq | 0.4.1 `hifi_revio_ssrs`, `VCF_results/clairs-to_COLO829T_Pacbio/` | 1.9.0 PACBIO_TUMOR_ONLY, `VCF_results/deepsomatic_COLO829T_Pacbio/COLO829T_Pacbio.chr1.vcf.gz` | SMaHT `...-M45-B004-uwsc-MERGED-pbmm2` (the fiber-seq reads of the tensors) |
| COLO829T_ONT | 0.4.1 `ont_r10_dorado_sup_5khz_ssrs`, `VCF_results/clairs-to_COLO829T_ONT/` | 1.9.0 ONT_TUMOR_ONLY, `VCF_results/deepsomatic_COLO829T_ONT/COLO829T_ONT.chr1.vcf.gz` | `COLO829T_ONT_merged.bam` (the reads of the tensors) |

- The DeepSomatic files scored are the full VCFs (PASS, PON, GERMLINE and RefCall records), not the PASS / per-kind
  subsets the old rtg runs used, so all four call sets of section 6 (`pass`, `pass_pon`, `all`, `all_pon`) exist.
- **HG008T ONT DeepSomatic used the PacBio model.** Its numbers are not a valid ONT comparator. User (2026-10-16): re-run
  it with `ONT_TUMOR_ONLY` and restrict the evaluation to the BED. Re-run on a tequila GPU (user: not on CPU;
  the first CPU chain, pull 417710 + run 417711, was cancelled): DeepSomatic `1.9.0-gpu` image pulled inside the job
  (`/scratch/jshen/bin/deepsomatic_1.9.0-gpu.sif`), the same BAM, `GRCh38.primary.fasta`, `--regions=chr1`,
  `--use_default_pon_filtering=true` as the user's other 1.9.0 runs (job 417841, `-p gpu -w tequila --gres=gpu:12gb:1`,
  16 CPUs, 80G; `tmp/reviewer_cpu_20261016/E8_rescore/run_deepsomatic_hg008t_ont_gpu.sh`; make_examples 50 min, 15 of 16
  shards in minutes; call_variants on the GPU finished, postprocess_variants was killed at 98 GB, one worker per node core
  each loading the PON VCF; postprocess re-run alone with `--cpus=4`, 64G, job 418308, `postprocess_hg008t_ont_gpu.sh`), calls in `calls/deepsomatic/HG008T_ONT/`, then
  scored like the others; the reported numbers are `eval_bed` (inside the HG008T benchmark BED), `eval_nobed` is whole
  chr1. Array tasks 10 / 12 (the wrong-model VCF) were cancelled so they do not take the same `score/` path.
- The old numbers of the original manuscript differ from this protocol: no BED (`--region chr1` only), the tools' own
  filters only (no common PoN), and for DeepSomatic HG008T PacBio SNV a baseline of SNVs + INDELs (FN 765 > 702 chr1 SNVs).
- Jobs: `tmp/reviewer_cpu_20261016/E8_rescore/` (`tasks.tsv`, `score.sh`; array 417688, 24 tasks, 4 at a time, 4 CPUs,
  32G each; the first submission 417661 failed in 1 s on anaconda's `bcftools` (missing `libcrypto.so.1.0.0`), the
  script now loads `bcftools/1.18`). Outputs: `score/<caller>/<dataset>/<KIND>/`. Then `collect.py` (labels set to
  ClairS-TO 0.4.1 / 0.4.2 and DeepSomatic 1.9.0).
- The COLO829T P25 / P50 / P75 comparator VCFs: checked and scored on 2026-10-29 (section 10). Not scored yet: the HapMap
  runs (`V2/HapMap/`, `V2/clairs-to_HapMap_WGS/`).
- **INDEL PoN rule.** The INDEL PoN rule changed during the re-scoring (user decision 2026-10-16: HP >= 8 bp AF >= 0.01,
  else >= 1e-4; `scripts/filter_panel_of_normals.py` rewritten 17:33, commit 71cb51a). The two HG008T Illumina INDEL
  reports (tasks 3 / 4) had run with the old rule (truth tagged 413 of 635; new rule 418): moved to
  `score/{clairs_to,deepsomatic}/HG008T_Illumina/INDEL.oldPoN_HP7_1e-3` and re-scored (array 418385). Every other INDEL
  report uses the new rule. Pansoma's INDEL rows: `collect.py` (`indel_new_pon_rule`) takes them from the CPU
  re-predictions re-scored with the new rule (`pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1*/eval_bed`,
  old reports kept as `eval_bed.pon_hp7_1e-3`; tmp/reviewer_cpu_20261013/pon_hp8_20261016) instead of the 7fc7396
  summaries; rows without a local re-scored call set (COLO829T d9 single-set runs) are marked "old INDEL PoN rule".
  SNV is unaffected.

## 10. COLO829T purity mixtures P75 / P50 / P25 (2026-10-29)

The user asked whether the ClairS-TO / DeepSomatic evaluations of the COLO829T purity mixtures are correct, and to redo
them for every purity except P100 (no re-calling if the existing calls are sound). Result: `purity_comparators.py` ->
`purity_comparators.md` / `.json`.

**Calls (the user's runs of the original manuscript, `V2/COLO829T_{WGS,Pacbio,ONT}/{clairs-to,deepsomatic-to}_INDEL_COLO829_*_chr1_P*/`;
each directory holds the run's full output despite the "INDEL" in its name).** All 18 are sound; none was re-called.
- ClairS-TO: every VCF header has `##clairs_to_version=0.4.2` and the command line: tumor-only, `--ctg_name chr1`,
  default PoNs, `--platform` `ilmn_ssrs` / `hifi_revio_ssrs` / `ont_r10_dorado_sup_5khz_ssrs`, and the mixture BAM of
  that purity (`Tumor_Normal_inSilico_Mix/COLO829_{WGS,Pacbio,ONT}/MIX_p{75,50,25}_*.mixed.sorted.bam`).
- DeepSomatic 1.9.0: the VCF header does not name the model (identical header lines for the WGS / PacBio / ONT models)
  and no logs were kept. Check by re-calling (`tmp/E8_purity_20261029/model_check.sh`, array 461256, tequila gpu:12gb,
  1.9.0-gpu image, default PON filtering): chr1:20,000,001-25,000,000 of each mixture BAM with the expected model
  (WGS_TUMOR_ONLY / PACBIO_TUMOR_ONLY / ONT_TUMOR_ONLY), records matched against the existing VCF in that window
  (`compare.py`). All nine reproduce the existing calls: PASS sets identical (Jaccard 1.000, 78-200 PASS records per
  window), all records Jaccard >= 0.9998, at most two records with another FILTER / GT, |QUAL diff| <= 0.9 (CPU vs GPU
  inference). Controls with the other long-read model on P75 differ clearly (PASS Jaccard 0.21 ONT with the PacBio model,
  0.48 fiber-seq with the ONT model). So each run used the right model on the right purity's BAM.
- Mixture BAMs (`/scratch/jshen/Scripts/mix_tumor_normal_{WGS,COLO829-Pacbio,COLO829-ONT}.sh`: `samtools view -s` of the
  SMaHT tumor and COLO829BL BAMs, seed 88, merged; read groups `tumor.sub` / `normal.sub`). Tumor sources are the BAMs of
  the P100 comparator runs (Illumina SMaHT sentieon bwa; fiber-seq MERGED pbmm2; ONT `COLO829T_ONT_merged.bam`). Bases
  per read group on chr1:20-22 Mb: tumor : normal = ONT 54 : 25 / 35 : 51 / 16 : 78 x, fiber-seq P50 30 : 48 x, Illumina
  P50 31 : 45 x; reads: tumor P75 / P25 = 3.0 and normal P25 / P75 = 3.0 on all three platforms. The tumor share is
  p x k / (p x k + 1 - p) with one local factor k = 0.63-0.70 (tumor copy number there) -> purities as labelled.

**Scoring.** The original manuscript's evaluation (`*_SNV_rtg_chr1_*_p*.somatic`) used `rtg vcfeval --region chr1`
without the BED and without our PoN (tools' own filters only), SNVs only, so it is not the protocol of the P100 rows.
Re-scored here exactly like P100 (section 9): `score_comparator.py`, set `COLO829T_<platform>_P<q>` (the same COLO829T
truth and BED as P100), SNV and INDEL, all four call sets (`tmp/E8_purity_20261029/score.sh`, array 461262, 36 tasks, all
COMPLETED; outputs `score/<caller>/COLO829T_<platform>_P<q>/<KIND>/`). PoN tags on the truth: SNV 84, INDEL 81 of chr1
(as P100).

**Findings** (PASS after our PoN, chr1, BED; `purity_comparators.md`):
- SNV F1 falls with purity for ClairS-TO on fiber-seq and ONT (0.762 / 0.688 / 0.659 / 0.614; 0.715 / 0.664 / 0.630 /
  0.565 for P100 / P75 / P50 / P25) and for DeepSomatic on Illumina and ONT.
- DeepSomatic fiber-seq and ONT recall is LOWER at P100 than at P75 (fiber-seq 0.555 vs 0.835; ONT 0.695 vs 0.756): at
  P100 it labels 858 (fiber-seq) / 292 (ONT) truth SNVs GERMLINE, against 103 / 72 at P75 and 6 / 4 at P25 (truth
  alleles by the FILTER of the caller's record). Clonal somatic variants at VAF ~0.5 in a pure tumor look germline to a
  tumor-only model; this is the caller's behaviour, not a scoring error.
- ClairS-TO Illumina P100 (0.551) is below its P75 (0.661): the P100 run is ClairS-TO 0.4.1 with `--platform ilmn`
  (`Pansoma_testing_results/VCF_results/clairs-to_COLO829T_WGS/`), the mixtures 0.4.2 with `ilmn_ssrs`; recall is about
  the same (0.764 vs 0.752), precision differs (0.431 vs 0.589). The fiber-seq / ONT P100 runs are 0.4.1 with the same
  model names as the mixtures. So the P100 column is not the same caller configuration for ClairS-TO Illumina.
- INDEL: low everywhere (F1 0.010-0.191); ClairS-TO collapses at P25 on long reads (fiber-seq 0.028, ONT 0.010).
- Against the original manuscript's SNV numbers, our PoN mostly raises precision: e.g. DeepSomatic Illumina P75 F1 0.324
  -> 0.500, ClairS-TO ONT P75 0.611 -> 0.664.
- Reads versus Pansoma's mixture tensors: ONT and fiber-seq mixtures come from the same mixing scripts and source BAMs
  as Pansoma's (fiber-seq: `giraffe_depthx100`, 100x BAM mode); Pansoma's Illumina mixtures are the FASTQ-mode
  `depth120x` files, a different random subset of the same tumor / normal libraries.

