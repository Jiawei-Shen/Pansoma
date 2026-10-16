# E11a: per-read alignment statistics of the same reads on four references

Status: complete (2026-10-15); see "Status" at the end. Tables: `results/tables.md` (written by `compare.py`, last run
2026-10-15, 18:36).

In short: both tumors on all three platforms. HG008T PacBio / ONT / Illumina and COLO829T fiber-seq / ONT are on all
four arms (GRCh38 BAM, GRCh38-only graph, HPRC v1.1 d9, HPRC v2.1 d46). COLO829T Illumina has the three graph arms only,
because no GRCh38 BAM is listed for it.

## Question

For identical reads, how do unmapped / perfect (full-length, edit-free) / mismatch / indel / soft-clip proportions,
MAPQ, and GRCh38 placement change between

- the GRCh38 BAMs used by the comparators (GIAB / SMaHT alignments),
- the GRCh38-only graph ("linear", vg 1.66 autoindex, H 17),
- HPRC v1.1 d9 (H 90),
- HPRC v2.1 d46 (H 464)?

Do reads that become perfect on a pangenome stay where they were (same GRCh38 locus)? How much of the gain comes from reads
routed over off-reference (branch) nodes? Does anything change differently in difficult regions?

## Reviewer points

- R2.7 (alignment statistics overstated; give counts and proportions).
- R3.6 (more perfect alignment is not more correct alignment; MAPQ distributions; difficult regions).
- R1.3 (difficult regions; ONT ambiguity; region-specific statistics).

This experiment answers the part about the same reads on several references: counts, MAPQ, placement concordance and
per-region statistics. It does not measure correctness against a known origin, which needs simulated reads (E11b). It
does not use the GIAB v3.5 stratifications either, because they are not on the server (E10).

## Data

GAMs (one sorted GAM + `.gai` per arm, all giraffe output of the same FASTQ; paths from
`pansoma_v2_tensors/<SET>[_d46|_linear]/run/config.json`, listed in `tasks.tsv`):

The three graph GAMs of each long-read set hold the same number of records (whole-GAM scans): HG008T PacBio 14,014,144,
HG008T ONT 3,642,157, COLO829T fiber-seq 27,846,735 and COLO829T ONT 14,448,438. Where a `vg stats -a` report exists
it gives the same count (section 0).

| set | arm | file | used by | state |
|---|---|---|---|---|
| HG008T_PacBio | linear | `.../Liss_lab_PacBio_Revio_20240125/vg_linear_hg38/HG008T_PacBio_HiFi_hg38_linear_vg.sorted.gam` | G scan (task 0) | done |
| HG008T_PacBio | d9 | `.../Liss_lab_PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x.AF-HPRC.sorted.gam` | G scan (task 1) | done |
| HG008T_PacBio | d46 | `.../Liss_lab_PacBio_Revio_20240125/vg_hprc_v2.1_d46/HG008T_PacBio_HiFi_hprc-v2.1-d46_vg.sorted.gam` | G scan (task 2) | done |
| HG008T_PacBio | grch38_bam | `.../Liss_lab_PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam` | W scan (wtask 4) | done |
| HG008T_ONT | linear | `.../Liss_lab_Northeastern-ONT-UL-20241216/vg_linear_hg38/HG008T_ONT_hg38_linear_vg.sorted.gam` | G scan (task 4) | done |
| HG008T_ONT | d9 | `.../Liss_lab_Northeastern-ONT-UL-20241216/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.AF-HPRC.sorted.gam` | G scan (task 5) | done |
| HG008T_ONT | d46 | `.../Liss_lab_Northeastern-ONT-UL-20241216/vg_hprc_v2.1_d46/HG008T_ONT_hprc-v2.1-d46_vg.sorted.gam` | G scan (task 6) | done |
| HG008T_ONT | grch38_bam | `.../Liss_lab_Northeastern-ONT-UL-20241216/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.bam` | W scan (wtask 5) | done |
| HG008T_Illumina | linear | `.../Liss_lab_BCM_Illumina-WGS_20240313/vg_linear_hg38/HG008T_BCM_IlluminaWGS_hg38_linear_vg.sorted.gam` | W scan (wtask 0) | done |
| HG008T_Illumina | d9 | `/scratch/jshen/data/HG008_GIAB/vg_giraffe_GAM/Liss_lab_BCM_Illumina-WGS_20240313/HG008-T_p23_BCM_IlluminaWGS_20240313.sorted.gam` | W scan (wtask 1) | done |
| HG008T_Illumina | d46 | `.../Liss_lab_BCM_Illumina-WGS_20240313/vg_hprc_v2.1_d46/HG008T_BCM_IlluminaWGS_hprc-v2.1-d46_vg.sorted.gam` | W scan (wtask 2) | done |
| HG008T_Illumina | grch38_bam | `.../Liss_lab_BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam` | W scan (wtask 3) | done |
| COLO829T_fiberseq | linear / d9 / d46 | `tasks.tsv` rows 12-14 | G scans | done |
| COLO829T_fiberseq | grch38_bam | `SMHTCOLO829T-X-X-M45-B004-uwsc-MERGED-pbmm2_1.13.0_GRCh38.aligned.sorted.bam` | W scan (wtask 9) | done |
| COLO829T_ONT | linear / d9 / d46 | `tasks.tsv` rows 16-18 | G scans | done |
| COLO829T_ONT | grch38_bam | `COLO829T_ONT_merged.bam` | W scan (wtask 10) | done |
| COLO829T_Illumina | linear / d9 / d46 | `wtasks.tsv` rows 6-8 (linear = the paired `vg_linear_hg38_pe` GAM) | W scans (wtasks 6-8) | done |

`...` = `/scratch/jshen/data/HG008_GIAB/raw_sequencing_data`. Full paths are in `tasks.tsv` / `wtasks.tsv`. No GRCh38 BAM
is listed for COLO829T Illumina, so that set has the three graph arms only.

Other inputs:

- GRCh38 placement of graph nodes: the pipeline's reference-path arrays (`chrom.npy`, `start0.npy`, `lengths.npy`,
  `path_nodes.npy`, `path_starts.npy`, `meta.json`) in
  `pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path`,
  `pansoma_v2_tensors/graph_index_hprc-v2.1-d46/hprc-v2.1-mc-grch38.d46.grch38_path` and
  `Linear_Reference_VG_Indexes/graph_index/hg38_linear.giraffe.grch38_path`.
- Region BEDs (`compare.py`, `BEDS`): SMaHT v2 easy / difficult / extreme
  (`/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_{easy,difficult,extreme}_v2.bed`), segmental duplications
  (`/scratch/qfu/real_mixture/repeatmasker/GRCh38_segdups.bed`), and the HG008T somatic benchmark v0.2 BED
  (`/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed`, HG008T only).
- Fixed GRCh38 windows (`choose_windows.py`, seed 20261013, chr1-22 + chrX):
  - `windows.tsv`, 48 core windows: 24 random windows of 500 kb and 24 segdup-centred windows of 100 kb (14.4 Mb in
    total).
  - `windows_extra.tsv`, 100 extra random windows of 500 kb (50 Mb), seed 20261014, used for the long-read BAM arm only.
  - `window_nodes.json` holds each window's node-ID range per graph (design check).
  - A copy of the core windows is in the data directory as `windows_core_v1.tsv` (identical to `windows.tsv`).
- Existing whole-GAM `vg stats -a` reports for some arms (`vg_stats_crosscheck.py`), used as a cross-check.

## Method

Every arm aligned the same reads, so a read is identified by (name, mate). Mate is 0 for unpaired reads and 1 / 2 from
`fragment_next` / `fragment_prev` (GAM) or flags 0x40 / 0x80 (BAM). Reads are sampled by a name hash,
`zlib.crc32(name) % MOD < KEEP`. The same reads are therefore selected in every arm, and both mates of a pair are kept.

Two samples (definitions in the `compare.py` docstring):

- **G, genome-wide (long reads, graph arms).** `scan_gam.py` walks a whole sorted GAM in parallel segments.
  - It counts every record: records, records with a path, and the MAPQ histogram.
  - It decodes the sampled reads (`MOD 100, KEEP 10`, a 10 % sample) with the v4 `gam_reader`.
  - Per read it records: mapped, MAPQ, mismatch / inserted / deleted bases and events, soft clips (an insertion at the
    alignment ends), off-reference node visits and their edits, and the GRCh38 placement. The placement is the contig
    holding most GRCh38-path bases of the read's nodes; its span runs from the first to the last of those nodes (node
    granularity).
  - A sampled read absent from an arm's GAM counts as unmapped there. None were absent in any arm of the four
    long-read sets.
- **W, fixed windows (all arms; the only design for short reads and for the BAM arm).** `window_scan.py` collects the
  reads placed in the windows, plus each arm's unmapped reads:
  - GAM: the GRCh38-path nodes in window ± flank, plus the node IDs between consecutive reference IDs at most 64 apart
    (bubble interiors). The GAI runs holding those nodes are walked, and every sampled record visiting one of them is
    kept. Unmapped reads come from the sorted GAM's head (GAI bin [0, 1]).
  - BAM: `fetch(chrom, start - flank, end + flank)`. Primary records only. Supplementary records seen in the same
    fetch are counted. Unmapped reads are the flag 0x4 records in the windows plus `fetch('*')`.
  - Short reads: 48 core windows, flank 1 kb, `MOD 10, KEEP 1` (10 % of reads within the windows).
  - Long reads: 148 windows (core + extra), flank 25 kb, `MOD 100, KEEP 10`. The long-read graph arms reuse the G
    full-scan rows, so only the BAM arm needs a window scan.
  - The read universe is every read that at least one arm places in a window. Each arm then has a status per read:
    `in` (placed overlapping the read's window ± flank), `elsewhere` (mapped but not near the window) or `unmapped`.

Per-read metrics (`compare.derived`, `arm_metrics`, `pair_metrics`):

- perfect: mapped, with no mismatch, no indel and no soft clip (BAM: no S/H either).
- mismatch / indel / soft clip: at least one such base, given as a percentage of mapped reads.
- identity: match / (match + mismatch + inserted + deleted bases).
- Placement concordant: same contig, and the spans overlap or lie within 1 kb of each other.
- Pairs (B relative to A, the same reads): MAPQ up / down, MAPQ ≥ 30 gained / lost, perfect gained / lost, and placement
  concordance. Concordance is reported overall and for reads with MAPQ ≥ 30 in both arms.
- "Became perfect" breakdown: the share on off-reference nodes in B, and what A had (mismatches only, an indel, a soft
  clip, unmapped).

Regions: each read gets one anchor placement, the same for every arm.

- G: the linear-graph placement; if absent, d9; then d46.
- W: the BAM placement when it is in the window; if not, linear, d9, then d46.
- BED strata need at least 50 % of the anchor span. SMaHT easy / difficult / extreme are disjoint (the one covering most
  of the span is used); segdup and the HG008T benchmark BED are separate, overlapping cuts.

Scripts (this folder):

| file | role |
|---|---|
| `scan_gam.py` | G full scan of one sorted GAM (`ROW` per-read schema) |
| `scan_bam.py` | full scan of one BAM (`BROW`); written but the full BAM scan was cancelled (see Status); `primary_row` is reused by `window_scan.py` |
| `choose_windows.py` | writes `windows.tsv`, `windows_extra.tsv`, `window_nodes.json` |
| `window_scan.py` | W scan of one GAM or BAM |
| `compare.py` | reads the finished scans and writes `results/<SET>.json` and `results/tables.md` |
| `vg_stats_crosscheck.py` | writes `results/vg_stats_crosscheck.md` from existing `vg stats -a` reports |
| `tasks.tsv`, `wtasks.tsv` | array task tables of the two job scripts |

Job scripts: `tmp/reviewer_cpu_20261013/E11a_alignment_stats/scan_array.sh` (G; general, 16 CPUs, 24 GB, 10 h) and
`window_array.sh` (W; general, 8 CPUs, 16 GB, 8 h). Logs are in `.../logs/scan_<task>_<array>.out` and
`win_<task>_<array>.out`. Outputs are in `/scratch/jshen/data/reviewer_experiments_20261013/E11a_alignment_stats/{scans,windows}/`.
`test/` there holds the development runs (segment or window subsets), not results.

| job | tasks | content | state |
|---|---|---|---|
| 387166 | 0, 1, 2, 4 | G scans HG008T PacBio linear / d9 / d46, HG008T ONT linear | COMPLETED 2026-10-13 (57 min - 1 h 23 min each; peak RSS 0.54-0.64 GB; CPU use 85-117 % of one core) |
| 387166 | 3 | full scan of the HG008T PacBio GRCh38 BAM | CANCELLED 2026-10-13 10:37 after 13 min; the BAM arm moved to the window design (387287_4) |
| 387166 | 5, 6 | G scans HG008T ONT d9 / d46 | COMPLETED 2026-10-15 13:05-13:37 (15 min, 17 min) |
| 387287 | 0, 1, 3 | W scans HG008T Illumina linear / d9 / GRCh38 BAM | COMPLETED 2026-10-13 (21 min, 17 min, 5 min; peak RSS <= 2.5 GB) |
| 387287 | 4, 5 | W scans HG008T PacBio / ONT GRCh38 BAM (148 windows) | COMPLETED 2026-10-13 (11 min each) |
| 387287 | 2 | W scan HG008T Illumina d46 | requeued 2026-10-13 11:22 and held; COMPLETED 2026-10-15 13:37-13:44 (7 min) |
| 387287 | 6-10 | W scans COLO829T Illumina linear / d9 / d46, COLO829T fiber-seq BAM, COLO829T ONT BAM | COMPLETED 2026-10-15 18:01-18:20 (3-5 min each) |
| 387288 | 12-14, 16-18 | G scans COLO829T fiber-seq and ONT, linear / d9 / d46 | COMPLETED 2026-10-15 13:44-18:01 (35-53 min each, one at a time) |
| 406873 | - | `compare.py` on the COLO829T sets (`sbatch --wrap`, general, 2 CPUs, 32 GB) | COMPLETED 2026-10-15 18:31-18:35 |

The jobs held on 2026-10-13 (lab-cluster cap) were released on 2026-10-15 at about 13:05. The `tasks.tsv` rows 7, 8-11,
15, 19 and 20-22 (full BAM scans and full Illumina GAM scans) were never submitted. Short reads and the BAM arms use the
W design only.

`compare.py` was last run on 2026-10-15: the HG008T JSONs were rewritten at 18:25-18:29, the COLO829T JSONs by job
406873, and `results/tables.md` at 18:36 from all six JSONs. All numbers below come from `results/<SET>.json`,
`results/tables.md`, `results/vg_stats_crosscheck.md` (2026-10-13) and the scan summaries in the data directory.

## Results

Percentages: perfect is a share of all reads in the row. Mismatch, indel, soft clip and the MAPQ bins are shares of
mapped reads. The full tables, including every stratum and pair, are in `results/tables.md`; the sections below copy
the main tables and selected strata. In the W pair tables the headers follow the JSON fields, because the header row of
these tables in `tables.md` is off by one cell (see Caveats).

### 0. Cross-check against whole-file `vg stats -a`

From `results/vg_stats_crosscheck.md`:

| set | arm | records | unmapped % | perfect % of all | gapless % of all | soft-clip events per 100 mapped reads | mean MAPQ |
|---|---|---:|---:|---:|---:|---:|---:|
| HG008T_PacBio | d9 | 14,014,144 | 2.64 | 1.39 | 2.40 | 2.13 | 58.19 |
| HG008T_ONT | d9 | 3,642,157 | 7.77 | 2.56 | 5.13 | 36.49 | 57.69 |
| HG008T_Illumina | d9 | 2,974,694,440 | 1.42 | 76.86 | 97.67 | 2.92 | 55.14 |
| HG008T_Illumina | linear | 2,974,694,440 | 1.51 | 67.54 | 95.56 | 3.17 | 55.32 |
| COLO829T_fiberseq | linear | 27,846,735 | 5.62 | 0.17 | 0.65 | 4.78 | 57.73 |
| COLO829T_ONT | d9 | 14,448,438 | 13.23 | 0.01 | 6.24 | 173.44 | 57.01 |
| COLO829T_ONT | linear | 14,448,438 | 13.53 | 0.01 | 5.83 | 173.45 | 57.26 |

The scans agree with these reports:

- The whole-GAM scans give the same record counts and unmapped rates: HG008T PacBio d9 2.64 %, HG008T ONT d9 7.77 %,
  COLO829T fiber-seq linear 5.62 %, COLO829T ONT d9 13.23 % and linear 13.53 % (sections 1, 4-6).
- HG008T PacBio d9, 10 % G sample: unmapped 2.64 %, perfect 1.39 %, mean MAPQ of mapped reads 58.196 (JSON
  `per_arm.d9.mapq_mean`).
- HG008T Illumina GAM heads: 44,910,000 unmapped records on the GRCh38-only graph (1.51 % of 2,974,694,440) and
  42,368,000 on d9 (1.42 %).

### 1. HG008T PacBio HiFi, G sample (1,402,384 reads, 10 % of 14,014,144), graph arms

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 14,014,144 | 389,438 | 2.78 | 0.04 | 1.86 | 95.88 | 1,402,384 |
| HPRC v1.1 d9 | 14,014,144 | 369,434 | 2.64 | 0.07 | 2.15 | 95.88 | 1,402,384 |
| HPRC v2.1 d46 | 14,014,144 | 1,592 | 0.01 | 0.07 | 4.80 | 83.32 | 1,402,384 |

Per arm, the same 1,402,384 reads:

| arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 2.79 | 0.15 | 94.53 | 99.43 | 2.65 | 0.04 | 2.19 | 0.72 | 1.16 | 95.89 | 0.9972 | 1.57 | 1.58 | 0.23 | - | - | - |
| HPRC v1.1 d9 | 2.64 | 1.39 | 67.29 | 97.54 | 2.00 | 0.07 | 2.44 | 0.72 | 0.89 | 95.89 | 0.9991 | 0.63 | 1.39 | 0.17 | 92.08 | 8.27 | 98.01 |
| HPRC v2.1 d46 | 0.01 | 1.42 | 66.40 | 97.61 | 3.75 | 0.07 | 5.47 | 0.96 | 10.15 | 83.35 | 0.9990 | 1.72 | 1.56 | 0.80 | 89.94 | 7.07 | 93.82 |

Same reads, B relative to A:

| B vs A | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,362,330 | 1.34 | 1.23 | 0.20 | 0.43 | 0.15 | 1.39 | 1.23 | 0.00 | 2.57 | 1.94 | 98.88 | 99.98 | 307 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,363,213 | 1.05 | 12.59 | 0.12 | 1.03 | 0.15 | 1.42 | 1.26 | 0.00 | 2.57 | 3.75 | 98.57 | 99.98 | 291 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,365,262 | 1.03 | 12.37 | 0.11 | 0.87 | 1.39 | 1.42 | 0.16 | 0.13 | 1.94 | 3.75 | 98.68 | 99.97 | 454 |

The reads with MAPQ ≥ 30 in both arms number 1,316,930 (d9 vs linear), 1,308,877 (d46 vs linear) and 1,309,406 (d46 vs
d9). Reads unmapped only on the GRCh38-only graph number 2,986 against d9 and 38,985 against d46; reads unmapped only on
d46 against the GRCh38-only graph number 52 (JSON `pairs`).

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 17,267 | 95.99 | 11.0 | 15.60 | 84.32 | 0.49 | 0.03 | 9.0 | 99.72 | 1.94 | 0.38 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 17,721 | 96.16 | 12.0 | 15.71 | 84.19 | 0.47 | 0.05 | 9.0 | 99.77 | 0.45 | 32.43 | 3 | 66.67 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,221 | 97.07 | 15.0 | 47.59 | 52.32 | 0.05 | 0.09 | 1.0 | 99.46 | 0.54 | 29.70 | 1,770 | 89.77 |

By region. The anchor placement comes from the GRCh38-only graph for 1,363,265 reads, d9 for 2,718, d46 for 35,097, and
none for 1,304. Selected rows:

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 1,273,579 | GRCh38-only graph | 0.01 | 0.16 | 94.35 | 99.44 | 1.94 | 0.00 | 99.14 | 0.9973 | - |
| SMaHT easy | 1,273,579 | HPRC v1.1 d9 | 0.00 | 1.46 | 66.17 | 97.51 | 1.38 | 0.00 | 99.37 | 0.9991 | 93.41 |
| SMaHT easy | 1,273,579 | HPRC v2.1 d46 | 0.00 | 1.50 | 64.28 | 97.54 | 1.38 | 0.00 | 89.45 | 0.9991 | 93.07 |
| SMaHT difficult | 5,066 | GRCh38-only graph | 9.59 | 0.12 | 96.86 | 98.58 | 7.51 | 0.00 | 94.85 | 0.9961 | - |
| SMaHT difficult | 5,066 | HPRC v1.1 d9 | 8.19 | 1.58 | 74.37 | 96.13 | 2.82 | 0.02 | 97.76 | 0.9990 | 95.38 |
| SMaHT difficult | 5,066 | HPRC v2.1 d46 | 0.00 | 1.76 | 74.40 | 96.09 | 10.74 | 0.00 | 77.48 | 0.9988 | 95.50 |
| SMaHT extreme | 109,232 | GRCh38-only graph | 29.30 | 0.07 | 97.02 | 99.48 | 10.69 | 0.66 | 48.75 | 0.9949 | - |
| SMaHT extreme | 109,232 | HPRC v1.1 d9 | 28.15 | 0.63 | 81.85 | 97.99 | 8.32 | 1.18 | 46.31 | 0.9982 | 76.48 |
| SMaHT extreme | 109,232 | HPRC v2.1 d46 | 0.03 | 0.62 | 86.54 | 98.31 | 23.55 | 0.85 | 22.52 | 0.9948 | 61.38 |
| segdup | 79,263 | GRCh38-only graph | 10.82 | 0.09 | 96.29 | 99.52 | 6.80 | 0.26 | 61.91 | 0.9962 | - |
| segdup | 79,263 | HPRC v1.1 d9 | 10.80 | 0.88 | 76.40 | 97.95 | 4.62 | 0.84 | 59.65 | 0.9988 | 91.24 |
| segdup | 79,263 | HPRC v2.1 d46 | 0.01 | 0.88 | 78.20 | 98.20 | 13.90 | 0.10 | 40.34 | 0.9985 | 85.83 |
| outside HG008T benchmark BED | 63,403 | GRCh38-only graph | 46.95 | 0.02 | 98.71 | 99.70 | 17.18 | 1.21 | 27.59 | 0.9875 | - |
| outside HG008T benchmark BED | 63,403 | HPRC v1.1 d9 | 45.32 | 0.19 | 91.38 | 98.90 | 14.48 | 1.36 | 25.42 | 0.9932 | 57.94 |
| outside HG008T benchmark BED | 63,403 | HPRC v2.1 d46 | 0.05 | 0.20 | 94.97 | 98.87 | 35.64 | 1.44 | 8.48 | 0.9600 | 40.61 |

Per-base rates by region, from `results/HG008T_PacBio.json` (`G.strata.<stratum>.per_arm`):

| stratum | arm | mismatch bp / kb | indel events / kb | soft-clipped bp % | mean MAPQ | MAPQ 30-59 % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | GRCh38-only graph | 1.020 | 1.562 | 0.089 | 59.923 | 0.782 | - | - |
| SMaHT easy | HPRC v1.1 d9 | 0.155 | 1.365 | 0.053 | 59.930 | 0.536 | 7.681 | 98.611 |
| SMaHT easy | HPRC v2.1 d46 | 0.150 | 1.368 | 0.054 | 59.263 | 10.401 | 4.647 | 98.637 |
| SMaHT extreme | GRCh38-only graph | 9.159 | 1.896 | 1.570 | 34.937 | n/r | - | - |
| SMaHT extreme | HPRC v1.1 d9 | 6.857 | 1.666 | 1.061 | 33.256 | n/r | 11.615 | 95.499 |
| SMaHT extreme | HPRC v2.1 d46 | 17.635 | 3.300 | 6.528 | 20.498 | n/r | 14.493 | 88.911 |
| segdup | GRCh38-only graph | 4.367 | 1.867 | 0.870 | 40.852 | n/r | - | - |
| segdup | HPRC v1.1 d9 | 2.413 | 1.590 | 0.573 | 39.081 | n/r | 10.852 | 95.740 |
| segdup | HPRC v2.1 d46 | 8.193 | 2.985 | 2.311 | 31.579 | n/r | 11.577 | 91.318 |

(n/r = in the JSON, not extracted here.)

Pairs by region (B relative to A):

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.03 | 0.04 | 1.30 | 0.00 | 99.98 | 4 | 96.05 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.03 | 0.09 | 1.34 | 0.00 | 99.98 | 3 | 96.20 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.01 | 0.07 | 0.17 | 0.13 | 99.98 | 8 | 97.12 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 0.90 | 0.37 | 1.46 | 0.00 | 99.58 | 6 | 97.30 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 0.90 | 1.09 | 1.64 | 0.00 | 99.28 | 11 | 97.59 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.26 | 1.55 | 0.20 | 0.02 | 99.08 | 5 | 100.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 3.05 | 6.39 | 0.55 | 0.00 | 82.74 | 249 | 94.20 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 1.45 | 13.97 | 0.55 | 0.00 | 78.56 | 219 | 94.64 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.47 | 11.54 | 0.08 | 0.09 | 80.75 | 239 | 95.65 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 2.22 | 5.35 | 0.79 | 0.00 | 88.08 | 123 | 94.75 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1.27 | 11.76 | 0.79 | 0.00 | 86.21 | 139 | 94.92 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.37 | 8.73 | 0.11 | 0.11 | 89.39 | 160 | 95.45 |

### 2. HG008T PacBio HiFi, W sample (148 windows, flank 25 kb; 34,332 reads: random 5,984, random_extra 26,125, segdup 2,223), all four arms

The graph arms use their genome-wide G rows, so their status is exact. The BAM arm comes from the window scan, and the
BAM has 1,055 unplaced unmapped records in the whole file.

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 97.88 | 2.12 | 0.00 |
| GRCh38-only graph | genome | 96.23 | 3.35 | 0.43 |
| HPRC v1.1 d9 | genome | 97.04 | 2.66 | 0.30 |
| HPRC v2.1 d46 | genome | 96.76 | 3.23 | 0.00 |

Reads the arm places in the window:

| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 33,604 | 0.18 | 95.01 | 99.44 | 6.98 | 0.03 | 3.03 | 0.60 | 0.59 | 95.75 | 0.9971 | 1.55 | 1.60 | 0.71 | - | - | - |
| GRCh38-only graph | 33,036 | 0.18 | 95.02 | 99.44 | 2.27 | 0.06 | 2.46 | 0.54 | 1.08 | 95.87 | 0.9972 | 1.36 | 1.58 | 0.16 | - | - | - |
| HPRC v1.1 d9 | 33,316 | 1.52 | 67.20 | 97.51 | 1.61 | 0.23 | 3.28 | 0.43 | 0.72 | 95.33 | 0.9991 | 0.42 | 1.38 | 0.08 | 93.40 | 7.59 | 98.56 |
| HPRC v2.1 d46 | 33,221 | 1.54 | 65.21 | 97.57 | 1.83 | 0.00 | 3.56 | 0.81 | 10.48 | 85.14 | 0.9991 | 0.42 | 1.40 | 0.16 | 92.96 | 4.77 | 98.50 |

The BAM arm has at least one supplementary record (inside the same fetch) for 0.94 % of its in-window reads (JSON
`W.per_arm.grch38_bam.supplementary_pct`).

Same reads, both arms place them in the window (B relative to A). The last four columns count the reads one arm places in
the window that the other places elsewhere or leaves unmapped. The headers here follow the JSON fields; the header row
of this table in `tables.md` is off by one cell (see Caveats).

| B vs A | both in window | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | discordant with MAPQ>=30 in both | A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 32,780 | 1.06 | 1.16 | 0.12 | 0.38 | 0.18 | 1.54 | 1.36 | 0.00 | 2.15 | 1.55 | 99.79 | 1 | 255 | 1 | 523 | 13 |
| HPRC v2.1 d46 vs GRCh38-only graph | 32,692 | 0.88 | 12.17 | 0.06 | 0.80 | 0.18 | 1.56 | 1.37 | 0.00 | 2.17 | 1.54 | 99.75 | 1 | 344 | 0 | 429 | 100 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 32,897 | 1.20 | 11.79 | 0.09 | 0.60 | 1.54 | 1.55 | 0.12 | 0.12 | 1.52 | 1.52 | 99.77 | 1 | 419 | 0 | 247 | 77 |
| GRCh38-only graph vs GRCh38 BAM | 32,880 | 0.86 | 2.83 | 0.04 | 0.96 | 0.18 | 0.18 | 0.00 | 0.00 | 5.40 | 2.16 | 99.62 | 0 | 620 | 104 | 156 | 0 |
| HPRC v1.1 d9 vs GRCh38 BAM | 32,756 | 0.72 | 2.66 | 0.07 | 1.24 | 0.18 | 1.54 | 1.36 | 0.00 | 5.29 | 1.51 | 99.59 | 0 | 779 | 69 | 560 | 0 |
| HPRC v2.1 d46 vs GRCh38 BAM | 32,775 | 0.80 | 12.90 | 0.01 | 1.69 | 0.18 | 1.55 | 1.37 | 0.00 | 5.46 | 1.65 | 99.46 | 0 | 828 | 1 | 446 | 0 |

Became perfect:

- d9 vs the BAM: 446 reads, 97.53 % of them on off-ref nodes. On the BAM, 84.98 % had an indel and 14.80 % had mismatches
  only. Placement concordance is 100.00 %.
- d46 vs the BAM: 449 reads, 97.55 % on off-ref nodes, concordance 100.00 %.
- No read lost its perfect status in either comparison. No read became perfect on the GRCh38-only graph relative to the
  BAM (0).

By window kind and region (selected):

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 30,016 | GRCh38 BAM | 99.95 | 0.05 | 0.00 | 0.18 | 94.78 | 99.42 | 4.66 | 0.00 | 100.00 | 0.9973 |
| SMaHT easy | 30,016 | GRCh38-only graph | 99.93 | 0.05 | 0.02 | 0.18 | 94.87 | 99.43 | 1.78 | 0.00 | 99.17 | 0.9973 |
| SMaHT easy | 30,016 | HPRC v1.1 d9 | 100.00 | 0.00 | 0.00 | 1.57 | 66.07 | 97.48 | 1.35 | 0.00 | 99.42 | 0.9991 |
| SMaHT easy | 30,016 | HPRC v2.1 d46 | 99.98 | 0.02 | 0.00 | 1.59 | 63.97 | 97.54 | 1.35 | 0.00 | 89.49 | 0.9991 |
| windows: segdup | 2,223 | GRCh38 BAM | 82.64 | 17.36 | 0.00 | 0.16 | 95.70 | 99.56 | 12.63 | 0.22 | 79.91 | 0.9960 |
| windows: segdup | 2,223 | GRCh38-only graph | 80.34 | 18.76 | 0.90 | 0.17 | 95.86 | 99.61 | 7.17 | 0.45 | 72.51 | 0.9961 |
| windows: segdup | 2,223 | HPRC v1.1 d9 | 90.19 | 9.31 | 0.49 | 0.80 | 78.85 | 98.40 | 3.74 | 0.15 | 62.59 | 0.9989 |
| windows: segdup | 2,223 | HPRC v2.1 d46 | 88.08 | 11.92 | 0.00 | 0.77 | 78.24 | 98.52 | 4.09 | 0.00 | 51.07 | 0.9988 |
| SMaHT extreme | 4,156 | GRCh38 BAM | 82.82 | 17.18 | 0.00 | 0.15 | 96.83 | 99.59 | 26.61 | 0.26 | 58.51 | 0.9944 |
| SMaHT extreme | 4,156 | GRCh38-only graph | 69.49 | 27.33 | 3.18 | 0.17 | 96.36 | 99.52 | 6.99 | 0.69 | 61.57 | 0.9957 |
| SMaHT extreme | 4,156 | HPRC v1.1 d9 | 75.72 | 21.94 | 2.33 | 1.02 | 77.47 | 97.81 | 4.04 | 2.41 | 56.15 | 0.9988 |
| SMaHT extreme | 4,156 | HPRC v2.1 d46 | 73.53 | 26.44 | 0.02 | 1.01 | 76.70 | 97.97 | 6.48 | 0.00 | 42.54 | 0.9987 |

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 1,695 | 0.47 | 5.13 | 0.77 | 0.00 | 99.47 | 5.10 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1,680 | 0.54 | 8.63 | 0.71 | 0.00 | 99.17 | 5.94 | 0.00 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 1,731 | 0.17 | 11.73 | 0.00 | 0.00 | 99.60 | 5.01 | 0.76 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 1,684 | 0.24 | 16.75 | 0.77 | 0.00 | 99.64 | 7.95 | 0.38 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 1,684 | 0.06 | 20.55 | 0.71 | 0.00 | 99.47 | 8.33 | 0.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 2,632 | 1.14 | 4.64 | 0.99 | 0.00 | 97.38 | 8.83 | 0.03 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 2,545 | 0.43 | 9.90 | 0.98 | 0.00 | 96.74 | 11.88 | 0.00 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 2,732 | 0.44 | 10.94 | 0.00 | 0.00 | 95.42 | 18.01 | 2.61 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 2,601 | 0.88 | 15.07 | 1.00 | 0.00 | 94.85 | 22.60 | 1.83 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 2,618 | 0.15 | 20.40 | 0.95 | 0.00 | 93.24 | 23.91 | 0.03 |

### 3. HG008T Illumina, W sample (48 windows, flank 1 kb; 1,465,754 reads: random 1,129,415, segdup 336,339), all four arms

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 92.17 | 7.79 | 0.04 |
| GRCh38-only graph | window | 91.35 | 7.87 | 0.78 |
| HPRC v1.1 d9 | window | 91.39 | 7.86 | 0.74 |
| HPRC v2.1 d46 | window | 91.36 | 7.89 | 0.75 |

Unmapped records of the whole file: GRCh38 BAM 2,165,824 (unplaced, `fetch('*')`); GRCh38-only graph 44,910,000
(1.51 % of all reads); HPRC v1.1 d9 42,368,000 (1.42 %); HPRC v2.1 d46 42,009,000 (1.41 %).

The read universe grew from 1,443,548 reads (2026-10-13, three arms) to 1,465,754 with the d46 arm, because the
universe is every read that at least one arm places in a window. The in-window rows of the first three arms are
unchanged; the status percentages above are over the larger universe.

Reads the arm places in the window:

| arm | reads in window | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 1,351,053 | 68.50 | 29.40 | 2.94 | 2.84 | 5.06 | 0.40 | 1.97 | 4.38 | 88.20 | 1.0000 | 4.73 | 0.21 | 0.83 | - | - | - |
| GRCh38-only graph | 1,338,964 | 68.85 | 29.40 | 2.96 | 2.19 | 3.43 | 3.31 | 2.77 | 3.61 | 86.89 | 1.0000 | 4.63 | 0.21 | 0.70 | - | - | - |
| HPRC v1.1 d9 | 1,339,624 | 79.15 | 20.09 | 0.62 | 1.94 | 4.29 | 3.44 | 2.67 | 3.36 | 86.24 | 1.0000 | 3.69 | 0.05 | 0.64 | 14.42 | 5.05 | 98.63 |
| HPRC v2.1 d46 | 1,339,133 | 79.24 | 20.01 | 0.63 | 1.92 | 4.24 | 3.27 | 2.66 | 3.54 | 86.29 | 1.0000 | 3.67 | 0.05 | 0.63 | 14.21 | 3.46 | 98.71 |

Same reads, both arms place them in the window (B relative to A). The last four columns count the reads one arm places in the window that the other places elsewhere or leaves unmapped:

| B vs A | both in window | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both | A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,277,850 | 1.59 | 2.71 | 0.43 | 0.92 | 69.12 | 79.36 | 10.24 | 0.00 | 2.05 | 1.83 | 99.71 | 100.00 | 19 | 60,526 | 588 | 60,967 | 807 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,277,730 | 1.54 | 3.26 | 0.41 | 0.92 | 69.07 | 79.41 | 10.34 | 0.00 | 2.04 | 1.84 | 99.64 | 100.00 | 19 | 60,650 | 584 | 60,688 | 715 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,293,604 | 0.96 | 1.38 | 0.40 | 0.38 | 79.30 | 79.41 | 0.29 | 0.17 | 1.86 | 1.87 | 99.75 | 100.00 | 0 | 45,757 | 263 | 45,315 | 214 |
| GRCh38-only graph vs GRCh38 BAM | 1,296,578 | 3.75 | 7.37 | 0.11 | 2.28 | 69.18 | 69.18 | 0.00 | 0.00 | 2.30 | 1.97 | 99.58 | 100.00 | 25 | 43,333 | 11,142 | 41,880 | 506 |
| HPRC v1.1 d9 vs GRCh38 BAM | 1,270,875 | 3.39 | 6.85 | 0.26 | 2.89 | 69.26 | 79.61 | 10.34 | 0.00 | 2.27 | 1.72 | 99.72 | 100.00 | 42 | 69,543 | 10,635 | 68,255 | 494 |
| HPRC v2.1 d46 vs GRCh38 BAM | 1,275,980 | 3.39 | 7.07 | 0.24 | 2.89 | 69.24 | 79.65 | 10.41 | 0.00 | 2.26 | 1.72 | 99.60 | 100.00 | 43 | 64,381 | 10,692 | 62,664 | 489 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 130,826 | 94.15 | 1.0 | 86.02 | 13.37 | 0.69 | 0.00 | 1.0 | 99.74 | 5.31 | 1.79 | 4 | 25.00 |
| HPRC v2.1 d46 vs GRCh38-only graph | 132,115 | 94.12 | 1.0 | 86.09 | 13.29 | 0.72 | 0.00 | 1.0 | 99.72 | 5.05 | 1.74 | 5 | 20.00 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 3,687 | 93.25 | 1.0 | 86.06 | 12.99 | 1.06 | 0.00 | 1.0 | 99.40 | 8.06 | 1.71 | 2,241 | 17.89 |
| GRCh38-only graph vs GRCh38 BAM | 22 | 0.00 | 0.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 4.55 | 40.91 | 9.09 | 40 | 0.00 |
| HPRC v1.1 d9 vs GRCh38 BAM | 131,452 | 94.16 | 1.0 | 84.74 | 13.93 | 1.43 | 0.00 | 1.0 | 99.72 | 6.35 | 3.89 | 13 | 23.08 |
| HPRC v2.1 d46 vs GRCh38 BAM | 132,810 | 94.14 | 1.0 | 84.84 | 13.83 | 1.44 | 0.00 | 1.0 | 99.72 | 6.17 | 3.82 | 13 | 30.77 |

Reads with MAPQ >= 30 in both arms (JSON `W.pairs`): 1,192,066 (d9 vs linear), 1,195,562 (d46 vs linear), 1,194,100
(d46 vs d9), 1,207,683 (linear vs BAM), 1,192,470 (d9 vs BAM), 1,196,020 (d46 vs BAM).

By window kind and region (selected; all strata in `tables.md`):

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 1,129,415 | GRCh38 BAM | 96.71 | 3.26 | 0.04 | 69.72 | 28.18 | 2.83 | 2.68 | 2.28 | 95.33 | 1.0000 |
| windows: random | 1,129,415 | GRCh38-only graph | 96.02 | 3.16 | 0.82 | 70.17 | 28.07 | 2.81 | 2.05 | 1.52 | 93.72 | 1.0000 |
| windows: random | 1,129,415 | HPRC v1.1 d9 | 95.97 | 3.25 | 0.78 | 80.25 | 19.03 | 0.45 | 1.79 | 1.76 | 93.78 | 1.0000 |
| windows: random | 1,129,415 | HPRC v2.1 d46 | 95.55 | 3.67 | 0.78 | 80.54 | 18.78 | 0.41 | 1.78 | 1.66 | 93.88 | 1.0000 |
| windows: segdup | 336,339 | GRCh38 BAM | 76.96 | 23.00 | 0.04 | 63.37 | 34.56 | 3.39 | 3.52 | 16.77 | 58.09 | 1.0000 |
| windows: segdup | 336,339 | GRCh38-only graph | 75.66 | 23.69 | 0.65 | 63.20 | 35.03 | 3.60 | 2.82 | 11.56 | 57.75 | 1.0000 |
| windows: segdup | 336,339 | HPRC v1.1 d9 | 76.03 | 23.35 | 0.62 | 74.49 | 24.57 | 1.34 | 2.55 | 15.04 | 54.28 | 1.0000 |
| windows: segdup | 336,339 | HPRC v2.1 d46 | 77.30 | 22.07 | 0.63 | 73.84 | 25.13 | 1.57 | 2.52 | 14.91 | 54.80 | 1.0000 |
| SMaHT easy | 918,984 | GRCh38 BAM | 99.89 | 0.08 | 0.03 | 71.55 | 26.87 | 2.15 | 1.64 | 0.01 | 99.81 | 1.0000 |
| SMaHT easy | 918,984 | GRCh38-only graph | 99.28 | 0.08 | 0.65 | 71.94 | 26.70 | 2.13 | 1.21 | 0.12 | 97.77 | 1.0000 |
| SMaHT easy | 918,984 | HPRC v1.1 d9 | 99.26 | 0.11 | 0.63 | 81.89 | 17.68 | 0.30 | 1.06 | 0.16 | 97.75 | 1.0000 |
| SMaHT easy | 918,984 | HPRC v2.1 d46 | 99.28 | 0.09 | 0.63 | 82.04 | 17.53 | 0.28 | 1.07 | 0.18 | 97.45 | 1.0000 |
| SMaHT difficult | 154,161 | GRCh38 BAM | 99.21 | 0.72 | 0.07 | 61.35 | 34.91 | 5.92 | 5.91 | 0.36 | 97.78 | 1.0000 |
| SMaHT difficult | 154,161 | GRCh38-only graph | 98.46 | 0.60 | 0.95 | 61.77 | 34.97 | 5.94 | 5.02 | 0.52 | 94.55 | 1.0000 |
| SMaHT difficult | 154,161 | HPRC v1.1 d9 | 98.21 | 0.99 | 0.80 | 73.31 | 25.15 | 1.07 | 4.36 | 0.40 | 94.88 | 1.0000 |
| SMaHT difficult | 154,161 | HPRC v2.1 d46 | 98.28 | 0.91 | 0.81 | 73.56 | 24.92 | 1.01 | 4.37 | 0.57 | 94.55 | 1.0000 |
| SMaHT extreme | 392,609 | GRCh38 BAM | 71.36 | 28.59 | 0.05 | 62.44 | 34.69 | 3.90 | 5.11 | 24.17 | 44.92 | 1.0000 |
| SMaHT extreme | 392,609 | GRCh38-only graph | 70.01 | 28.96 | 1.03 | 62.50 | 35.26 | 4.07 | 3.92 | 16.03 | 46.53 | 1.0000 |
| SMaHT extreme | 392,609 | HPRC v1.1 d9 | 70.31 | 28.70 | 0.99 | 73.30 | 25.30 | 1.44 | 3.49 | 20.09 | 43.45 | 1.0000 |
| SMaHT extreme | 392,609 | HPRC v2.1 d46 | 70.11 | 28.90 | 0.99 | 73.07 | 25.52 | 1.59 | 3.41 | 19.70 | 44.76 | 1.0000 |

Pairs by region (B relative to A):

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 211,014 | 1.54 | 4.25 | 10.92 | 0.00 | 98.80 | 17.03 | 0.05 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 215,805 | 1.53 | 3.44 | 10.23 | 0.00 | 98.30 | 15.15 | 0.05 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 229,854 | 1.73 | 0.89 | 0.09 | 0.61 | 99.02 | 10.08 | 0.03 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 228,323 | 0.35 | 6.51 | 0.01 | 0.01 | 98.03 | 10.98 | 0.81 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 207,201 | 1.00 | 10.12 | 11.43 | 0.00 | 98.81 | 19.18 | 0.77 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 213,813 | 0.95 | 9.24 | 10.64 | 0.00 | 98.10 | 16.61 | 0.78 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 911,336 | 0.06 | 0.15 | 9.97 | 0.00 | 100.00 | 0.08 | 0.03 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 911,692 | 0.06 | 0.24 | 10.12 | 0.00 | 100.00 | 0.04 | 0.03 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 911,830 | 0.10 | 0.20 | 0.25 | 0.10 | 100.00 | 0.03 | 0.01 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 911,778 | 0.01 | 0.97 | 0.00 | 0.00 | 100.00 | 0.03 | 0.64 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 911,289 | 0.01 | 1.07 | 9.98 | 0.00 | 100.00 | 0.10 | 0.63 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 911,629 | 0.01 | 1.17 | 10.13 | 0.00 | 100.00 | 0.06 | 0.63 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 215,457 | 2.11 | 4.61 | 10.56 | 0.00 | 98.31 | 21.52 | 0.09 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 214,894 | 2.05 | 4.14 | 10.36 | 0.00 | 97.86 | 21.73 | 0.09 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 230,511 | 1.74 | 1.15 | 0.30 | 0.44 | 98.64 | 16.46 | 0.03 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 233,966 | 0.52 | 7.74 | 0.01 | 0.01 | 97.71 | 15.11 | 1.38 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 208,985 | 1.49 | 11.54 | 11.11 | 0.00 | 98.32 | 24.09 | 1.32 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 213,644 | 1.33 | 10.86 | 10.67 | 0.00 | 97.64 | 22.42 | 1.32 |

### 4. HG008T ONT

**G sample (364,076 reads, 10 % of 3,642,157), graph arms.**

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 3,642,157 | 287,106 | 7.88 | 0.93 | 2.26 | 95.18 | 364,076 |
| HPRC v1.1 d9 | 3,642,157 | 283,090 | 7.77 | 1.00 | 2.55 | 95.02 | 364,076 |
| HPRC v2.1 d46 | 3,642,157 | 292,931 | 8.04 | 0.14 | 1.71 | 94.16 | 364,076 |

Per arm, the same reads:

| arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 7.86 | 1.32 | 92.47 | 95.19 | 30.37 | 0.93 | 1.79 | 1.10 | 1.00 | 95.19 | 0.9941 | 5.22 | 4.27 | 1.76 | - | - | - |
| HPRC v1.1 d9 | 7.75 | 2.54 | 86.79 | 94.47 | 29.89 | 1.00 | 2.03 | 0.99 | 0.94 | 95.03 | 0.9964 | 4.26 | 4.12 | 1.64 | 81.48 | 8.82 | 96.42 |
| HPRC v2.1 d46 | 8.01 | 2.50 | 87.67 | 95.38 | 29.61 | 0.14 | 2.27 | 1.80 | 1.61 | 94.18 | 0.9963 | 4.35 | 4.33 | 1.71 | 81.27 | 6.88 | 93.46 |

Same reads, B relative to A:

| B vs A | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 334,048 | 1.56 | 2.13 | 0.47 | 0.72 | 1.32 | 2.54 | 1.22 | 0.00 | 27.98 | 27.57 | 97.45 | 99.68 | 1,019 |
| HPRC v2.1 d46 vs GRCh38-only graph | 330,369 | 2.06 | 3.58 | 0.72 | 1.68 | 1.32 | 2.50 | 1.20 | 0.02 | 27.98 | 27.24 | 97.23 | 99.38 | 1,953 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 331,036 | 2.36 | 3.17 | 0.79 | 1.48 | 2.54 | 2.50 | 0.06 | 0.10 | 27.57 | 27.24 | 97.40 | 99.39 | 1,915 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 4,458 | 95.56 | 2.0 | 67.50 | 32.12 | 0.61 | 0.18 | 1.0 | 98.67 | 1.66 | 1.57 | 18 | 44.44 |
| HPRC v2.1 d46 vs GRCh38-only graph | 4,374 | 95.61 | 2.0 | 69.25 | 30.32 | 0.48 | 0.18 | 1.0 | 98.81 | 1.60 | 1.63 | 88 | 53.41 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 213 | 85.92 | 2.0 | 54.46 | 44.60 | 0.94 | 0.00 | 1.0 | 95.31 | 3.29 | 4.23 | 367 | 66.21 |

Region anchor: linear 335,447, d9 1,693, d46 2,996, none 23,940. Reads unmapped only on the GRCh38-only graph number
1,814 against d9 and 4,537 against d46; reads unmapped only on d9 / d46 against the GRCh38-only graph number 1,399 /
5,078 (JSON `G.pairs`). By region (selected):

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 298,748 | GRCh38-only graph | 0.44 | 1.39 | 93.78 | 96.56 | 28.31 | 0.01 | 99.64 | 0.9943 | - |
| SMaHT easy | 298,748 | HPRC v1.1 d9 | 0.26 | 2.62 | 88.03 | 95.89 | 27.97 | 0.02 | 99.54 | 0.9965 | 83.69 |
| SMaHT easy | 298,748 | HPRC v2.1 d46 | 0.18 | 2.59 | 87.97 | 95.94 | 28.14 | 0.00 | 98.86 | 0.9965 | 83.50 |
| SMaHT difficult | 5,168 | GRCh38-only graph | 2.94 | 3.52 | 75.94 | 82.00 | 38.38 | 7.66 | 87.96 | 0.9934 | - |
| SMaHT difficult | 5,168 | HPRC v1.1 d9 | 2.03 | 6.97 | 62.41 | 80.09 | 37.61 | 7.58 | 87.18 | 0.9961 | 64.72 |
| SMaHT difficult | 5,168 | HPRC v2.1 d46 | 10.57 | 6.79 | 65.51 | 86.59 | 32.39 | 0.48 | 87.91 | 0.9957 | 65.36 |
| SMaHT extreme | 31,583 | GRCh38-only graph | 9.52 | 1.46 | 81.22 | 83.31 | 45.40 | 9.27 | 53.85 | 0.9908 | - |
| SMaHT extreme | 31,583 | HPRC v1.1 d9 | 10.50 | 3.23 | 76.93 | 82.03 | 43.52 | 10.09 | 53.10 | 0.9946 | 66.06 |
| SMaHT extreme | 31,583 | HPRC v2.1 d46 | 13.99 | 3.04 | 86.79 | 90.67 | 38.93 | 1.35 | 49.56 | 0.9920 | 64.47 |
| segdup | 19,187 | GRCh38-only graph | 5.26 | 1.29 | 93.54 | 95.31 | 34.44 | 0.99 | 65.07 | 0.9919 | - |
| segdup | 19,187 | HPRC v1.1 d9 | 4.86 | 2.48 | 88.32 | 94.36 | 33.09 | 1.65 | 62.48 | 0.9953 | 80.93 |
| segdup | 19,187 | HPRC v2.1 d46 | 1.34 | 2.23 | 90.30 | 95.41 | 34.72 | 0.87 | 58.23 | 0.9947 | 78.38 |
| outside HG008T benchmark BED | 13,009 | GRCh38-only graph | 16.78 | 0.56 | 92.35 | 89.01 | 44.01 | 5.19 | 46.43 | 0.9754 | - |
| outside HG008T benchmark BED | 13,009 | HPRC v1.1 d9 | 16.92 | 0.99 | 90.85 | 88.32 | 43.38 | 6.60 | 44.43 | 0.9824 | 51.98 |
| outside HG008T benchmark BED | 13,009 | HPRC v2.1 d46 | 8.13 | 0.85 | 95.52 | 94.13 | 43.31 | 1.87 | 30.85 | 0.9141 | 47.43 |

Pairs by region (B relative to A):

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.03 | 0.10 | 1.24 | 0.00 | 99.97 | 11 | 95.48 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.11 | 0.48 | 1.22 | 0.02 | 99.78 | 328 | 95.64 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.17 | 0.46 | 0.05 | 0.08 | 99.77 | 323 | 88.89 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 0.82 | 1.41 | 3.44 | 0.00 | 90.56 | 14 | 92.70 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 0.82 | 5.30 | 3.35 | 0.08 | 94.66 | 80 | 92.49 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.34 | 5.22 | 0.23 | 0.41 | 95.07 | 83 | 91.67 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 4.93 | 6.92 | 1.79 | 0.02 | 72.77 | 925 | 96.82 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 7.09 | 14.12 | 1.70 | 0.12 | 70.38 | 1,285 | 96.27 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7.22 | 11.89 | 0.12 | 0.32 | 72.66 | 1,190 | 71.05 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 3.16 | 6.41 | 1.19 | 0.01 | 84.92 | 342 | 96.07 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 4.96 | 10.88 | 1.07 | 0.13 | 76.98 | 846 | 96.12 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 5.40 | 8.06 | 0.10 | 0.35 | 78.56 | 757 | 85.00 |

**W sample (148 windows, flank 25 kb; 9,000 reads: random 1,616, random_extra 6,753, segdup 631), all four arms.** The
BAM has 141,830 unplaced unmapped records in the whole file, and at least one supplementary record (inside the same
fetch) for 1.42 % of its in-window reads (JSON `W.per_arm.grch38_bam.supplementary_pct`). The universe was 8,614 reads
on 2026-10-13, with two arms.

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 94.20 | 4.07 | 1.73 |
| GRCh38-only graph | genome | 92.29 | 5.21 | 2.50 |
| HPRC v1.1 d9 | genome | 93.62 | 4.07 | 2.31 |
| HPRC v2.1 d46 | genome | 92.42 | 3.70 | 3.88 |

Reads the arm places in the window:

| arm | reads in window | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 8,478 | 1.07 | 94.17 | 96.74 | 36.77 | 0.42 | 2.24 | 0.83 | 1.40 | 95.10 | 0.9940 | 5.64 | 5.17 | 2.40 | - | - | - |
| GRCh38-only graph | 8,306 | 1.06 | 93.82 | 96.16 | 28.91 | 0.53 | 1.85 | 1.10 | 0.88 | 95.64 | 0.9941 | 4.95 | 4.41 | 1.35 | - | - | - |
| HPRC v1.1 d9 | 8,426 | 2.44 | 87.10 | 94.48 | 29.44 | 1.51 | 2.33 | 1.03 | 0.95 | 94.18 | 0.9965 | 4.03 | 4.31 | 1.21 | 84.19 | 8.63 | 96.42 |
| HPRC v2.1 d46 | 8,318 | 2.42 | 88.45 | 95.84 | 28.59 | 0.16 | 2.16 | 1.61 | 1.38 | 94.69 | 0.9964 | 3.97 | 4.38 | 1.18 | 83.53 | 5.66 | 96.78 |

Same reads, both arms place them in the window (B relative to A). The last four columns count the reads one arm places in the window that the other places elsewhere or leaves unmapped:

| B vs A | both in window | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both | A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 8,156 | 1.37 | 1.58 | 0.29 | 0.42 | 1.08 | 2.46 | 1.39 | 0.00 | 28.37 | 28.18 | 99.61 | 99.91 | 7 | 143 | 7 | 249 | 21 |
| HPRC v2.1 d46 vs GRCh38-only graph | 8,055 | 1.51 | 2.15 | 0.42 | 0.72 | 1.09 | 2.46 | 1.37 | 0.00 | 28.08 | 27.73 | 99.45 | 99.94 | 5 | 188 | 63 | 223 | 40 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 8,096 | 2.00 | 2.00 | 0.63 | 0.75 | 2.53 | 2.47 | 0.06 | 0.12 | 27.95 | 27.80 | 99.49 | 99.97 | 2 | 202 | 128 | 197 | 25 |
| GRCh38-only graph vs GRCh38 BAM | 8,170 | 1.65 | 2.47 | 0.27 | 1.10 | 1.09 | 1.08 | 0.00 | 0.01 | 34.93 | 28.37 | 99.62 | 100.00 | 0 | 94 | 214 | 94 | 42 |
| HPRC v1.1 d9 vs GRCh38 BAM | 8,137 | 1.46 | 2.59 | 0.22 | 1.19 | 1.11 | 2.52 | 1.41 | 0.00 | 34.96 | 28.18 | 99.52 | 99.91 | 7 | 142 | 199 | 172 | 117 |
| HPRC v2.1 d46 vs GRCh38 BAM | 8,080 | 1.37 | 2.33 | 0.16 | 1.10 | 1.13 | 2.46 | 1.37 | 0.04 | 34.93 | 27.98 | 99.37 | 99.94 | 5 | 206 | 192 | 237 | 1 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 113 | 95.58 | 2.0 | 65.49 | 34.51 | 0.88 | 0.00 | 2.0 | 100.00 | 1.77 | 1.77 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 110 | 94.55 | 2.0 | 66.36 | 33.64 | 0.91 | 0.00 | 1.5 | 100.00 | 0.00 | 2.73 | 0 | n/a |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 5 | 80.00 | 2.0 | 0.00 | 100.00 | 0.00 | 0.00 | 0.0 | 100.00 | 20.00 | 20.00 | 10 | 60.00 |
| GRCh38-only graph vs GRCh38 BAM | 0 |  |  |  |  |  |  |  |  |  |  |  |  |
| HPRC v1.1 d9 vs GRCh38 BAM | 115 | 95.65 | 2.0 | 62.61 | 34.78 | 5.22 | 0.00 | 1.0 | 98.26 | 4.35 | 2.61 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38 BAM | 111 | 94.59 | 2.0 | 63.06 | 34.23 | 5.41 | 0.00 | 1.0 | 100.00 | 2.70 | 2.70 | 3 | 0.00 |

By window kind and region (selected):

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | 631 | GRCh38 BAM | 82.41 | 16.80 | 0.79 | 0.77 | 95.19 | 97.88 | 45.00 | 2.69 | 77.88 | 0.9918 |
| windows: segdup | 631 | GRCh38-only graph | 78.92 | 18.86 | 2.22 | 0.80 | 95.18 | 97.39 | 36.35 | 0.80 | 78.51 | 0.9922 |
| windows: segdup | 631 | HPRC v1.1 d9 | 82.57 | 15.53 | 1.90 | 2.11 | 92.71 | 96.55 | 35.12 | 0.38 | 72.17 | 0.9960 |
| windows: segdup | 631 | HPRC v2.1 d46 | 80.98 | 16.01 | 3.01 | 1.96 | 93.74 | 96.67 | 34.64 | 0.59 | 70.06 | 0.9957 |
| SMaHT easy | 7,529 | GRCh38 BAM | 98.63 | 0.24 | 1.13 | 0.97 | 94.60 | 97.06 | 35.87 | 0.00 | 99.19 | 0.9942 |
| SMaHT easy | 7,529 | GRCh38-only graph | 96.60 | 1.38 | 2.02 | 0.99 | 94.75 | 96.98 | 28.09 | 0.01 | 99.55 | 0.9943 |
| SMaHT easy | 7,529 | HPRC v1.1 d9 | 98.01 | 0.20 | 1.79 | 2.20 | 87.87 | 95.20 | 29.00 | 1.12 | 98.27 | 0.9966 |
| SMaHT easy | 7,529 | HPRC v2.1 d46 | 96.84 | 0.35 | 2.82 | 2.18 | 88.88 | 96.38 | 28.12 | 0.00 | 99.16 | 0.9965 |
| SMaHT extreme | 1,311 | GRCh38 BAM | 70.33 | 25.25 | 4.42 | 1.52 | 92.41 | 94.58 | 44.25 | 3.90 | 62.36 | 0.9913 |
| SMaHT extreme | 1,311 | GRCh38-only graph | 67.89 | 26.62 | 5.49 | 1.24 | 88.65 | 90.79 | 34.49 | 3.71 | 64.61 | 0.9926 |
| SMaHT extreme | 1,311 | HPRC v1.1 d9 | 69.41 | 25.10 | 5.49 | 3.41 | 84.29 | 90.11 | 32.42 | 3.96 | 61.32 | 0.9962 |
| SMaHT extreme | 1,311 | HPRC v2.1 d46 | 67.89 | 23.11 | 9.00 | 3.26 | 88.54 | 92.70 | 31.46 | 1.46 | 58.65 | 0.9954 |

Pairs by region (B relative to A):

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 453 | 0.44 | 2.87 | 1.32 | 0.00 | 99.78 | 8.84 | 0.20 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 432 | 2.08 | 3.70 | 0.93 | 0.00 | 98.15 | 11.65 | 1.61 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 455 | 3.74 | 2.86 | 0.22 | 0.66 | 98.02 | 12.28 | 0.38 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 477 | 0.84 | 8.18 | 0.00 | 0.00 | 100.00 | 5.96 | 2.31 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 453 | 0.66 | 10.38 | 1.55 | 0.00 | 99.78 | 10.77 | 2.12 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 435 | 0.69 | 8.51 | 0.92 | 0.00 | 97.93 | 13.85 | 2.50 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 7,272 | 0.07 | 0.07 | 1.24 | 0.00 | 100.00 | 0.00 | 0.01 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 7,248 | 0.18 | 0.28 | 1.20 | 0.00 | 100.00 | 0.19 | 0.15 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7,263 | 0.21 | 0.30 | 0.03 | 0.07 | 100.00 | 0.33 | 1.25 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 7,273 | 0.00 | 0.22 | 0.00 | 0.00 | 100.00 | 0.05 | 2.01 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 7,285 | 0.01 | 0.23 | 1.24 | 0.00 | 100.00 | 0.08 | 1.82 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 7,278 | 0.00 | 0.30 | 1.20 | 0.00 | 100.00 | 0.22 | 1.78 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 752 | 2.53 | 3.72 | 1.99 | 0.00 | 95.74 | 14.83 | 0.67 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 680 | 3.09 | 5.59 | 2.21 | 0.00 | 93.53 | 19.10 | 4.49 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 706 | 4.96 | 5.52 | 0.42 | 0.71 | 94.19 | 19.23 | 3.19 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 769 | 2.73 | 9.62 | 0.00 | 0.13 | 95.97 | 9.65 | 6.94 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 724 | 2.21 | 10.91 | 2.35 | 0.00 | 94.61 | 14.64 | 6.83 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 676 | 1.78 | 9.91 | 2.37 | 0.44 | 92.46 | 20.39 | 6.29 |

### 5. COLO829T fiber-seq

**G sample (2,687,506 distinct reads from 2,784,787 sampled records; whole GAM 27,846,735 records), graph arms.**

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 27,846,735 | 1,564,272 | 5.62 | 0.05 | 2.17 | 92.38 | 2,784,787 |
| HPRC v1.1 d9 | 27,846,735 | 1,510,899 | 5.43 | 0.10 | 2.48 | 92.63 | 2,784,787 |
| HPRC v2.1 d46 | 27,846,735 | 8,094 | 0.03 | 0.06 | 7.61 | 76.87 | 2,784,787 |

Per arm, the same reads:

| arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 3.37 | 0.18 | 95.81 | 99.31 | 3.26 | 0.04 | 2.48 | 0.86 | 3.50 | 93.12 | 0.9970 | 1.80 | 1.82 | 0.30 | - | - | - |
| HPRC v1.1 d9 | 3.23 | 1.50 | 74.69 | 97.12 | 2.58 | 0.08 | 2.71 | 0.84 | 3.03 | 93.34 | 0.9989 | 0.84 | 1.62 | 0.24 | 92.54 | 8.12 | 97.90 |
| HPRC v2.1 d46 | 0.01 | 1.54 | 74.05 | 97.23 | 4.96 | 0.04 | 6.34 | 1.17 | 13.14 | 79.31 | 0.9989 | 2.19 | 1.87 | 1.13 | 90.23 | 8.56 | 92.17 |

Same reads, B relative to A:

| B vs A | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 2,594,452 | 2.42 | 1.36 | 0.25 | 0.44 | 0.18 | 1.50 | 1.32 | 0.00 | 3.15 | 2.50 | 98.72 | 99.97 | 747 |
| HPRC v2.1 d46 vs GRCh38-only graph | 2,596,786 | 1.18 | 16.06 | 0.16 | 1.15 | 0.18 | 1.54 | 1.36 | 0.00 | 3.15 | 4.96 | 98.29 | 99.97 | 757 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,600,479 | 1.20 | 15.83 | 0.14 | 1.01 | 1.50 | 1.54 | 0.17 | 0.13 | 2.50 | 4.96 | 98.50 | 99.96 | 912 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 35,496 | 96.08 | 11.0 | 16.06 | 83.85 | 0.45 | 0.04 | 9.0 | 99.56 | 6.95 | 0.41 | 3 | 100.00 |
| HPRC v2.1 d46 vs GRCh38-only graph | 36,452 | 96.23 | 11.0 | 16.23 | 83.69 | 0.47 | 0.04 | 9.0 | 99.57 | 0.49 | 53.51 | 4 | 100.00 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 4,446 | 97.59 | 15.0 | 48.31 | 51.60 | 0.09 | 0.07 | 1.0 | 99.19 | 0.52 | 49.92 | 3,491 | 90.17 |

Region anchor: linear 2,596,923, d9 5,831, d46 80,838, none 3,914. Reads unmapped only on the GRCh38-only graph number
6,174 against d9 and 90,413 against d46; reads unmapped only on d9 / d46 number 2,471 / 137 (JSON `G.pairs`). By region
(selected; no HG008T benchmark BED for COLO829T):

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 2,411,244 | GRCh38-only graph | 0.01 | 0.19 | 95.63 | 99.31 | 2.34 | 0.00 | 96.78 | 0.9971 | - |
| SMaHT easy | 2,411,244 | HPRC v1.1 d9 | 0.00 | 1.60 | 73.71 | 97.08 | 1.76 | 0.00 | 97.22 | 0.9990 | 94.07 |
| SMaHT easy | 2,411,244 | HPRC v2.1 d46 | 0.00 | 1.64 | 72.14 | 97.12 | 1.77 | 0.00 | 86.15 | 0.9990 | 93.73 |
| SMaHT difficult | 9,660 | GRCh38-only graph | 9.49 | 0.04 | 97.47 | 98.68 | 8.81 | 0.01 | 89.87 | 0.9957 | - |
| SMaHT difficult | 9,660 | HPRC v1.1 d9 | 9.03 | 1.33 | 79.27 | 95.88 | 3.24 | 0.00 | 94.24 | 0.9989 | 95.91 |
| SMaHT difficult | 9,660 | HPRC v2.1 d46 | 0.00 | 1.31 | 80.48 | 96.07 | 11.60 | 0.01 | 72.29 | 0.9987 | 95.87 |
| SMaHT extreme | 229,856 | GRCh38-only graph | 32.76 | 0.08 | 97.92 | 99.32 | 12.67 | 0.59 | 45.32 | 0.9942 | - |
| SMaHT extreme | 229,856 | HPRC v1.1 d9 | 31.77 | 0.68 | 85.84 | 97.66 | 10.35 | 1.22 | 43.49 | 0.9977 | 77.01 |
| SMaHT extreme | 229,856 | HPRC v2.1 d46 | 0.04 | 0.70 | 89.94 | 98.14 | 29.72 | 0.51 | 19.21 | 0.9922 | 62.21 |
| segdup | 164,038 | GRCh38-only graph | 14.31 | 0.10 | 97.48 | 99.43 | 8.02 | 0.22 | 57.60 | 0.9958 | - |
| segdup | 164,038 | HPRC v1.1 d9 | 14.45 | 0.92 | 82.14 | 97.67 | 5.68 | 0.87 | 55.99 | 0.9985 | 90.70 |
| segdup | 164,038 | HPRC v2.1 d46 | 0.01 | 0.93 | 84.03 | 98.01 | 18.10 | 0.17 | 35.47 | 0.9980 | 85.38 |

Pairs by region (B relative to A):

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.03 | 0.04 | 1.41 | 0.00 | 99.99 | 13 | 96.03 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.03 | 0.13 | 1.45 | 0.00 | 99.98 | 11 | 96.18 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.01 | 0.11 | 0.17 | 0.14 | 99.98 | 15 | 97.58 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 1.57 | 0.47 | 1.28 | 0.00 | 99.61 | 8 | 96.77 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 1.53 | 1.30 | 1.27 | 0.00 | 99.33 | 9 | 96.75 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.24 | 1.15 | 0.04 | 0.05 | 99.50 | 4 | 100.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 3.40 | 6.26 | 0.60 | 0.00 | 82.11 | 579 | 97.16 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 1.77 | 14.20 | 0.62 | 0.00 | 76.99 | 563 | 97.32 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.78 | 12.16 | 0.10 | 0.08 | 80.27 | 528 | 97.77 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 2.49 | 5.22 | 0.82 | 0.00 | 86.82 | 323 | 96.86 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1.57 | 11.90 | 0.83 | 0.00 | 84.06 | 402 | 96.84 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.66 | 9.18 | 0.13 | 0.12 | 88.22 | 343 | 98.58 |

**W sample (148 windows, flank 25 kb; 62,786 reads: random 12,286, random_extra 46,047, segdup 4,453), all four
arms.** The BAM has 248,616 unplaced unmapped records in the whole file; supplementary records for 0.16 % of its
in-window reads.

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 97.65 | 2.35 | 0.00 |
| GRCh38-only graph | genome | 97.13 | 2.53 | 0.33 |
| HPRC v1.1 d9 | genome | 98.05 | 1.68 | 0.27 |
| HPRC v2.1 d46 | genome | 97.61 | 2.39 | 0.00 |

Reads the arm places in the window:

| arm | reads in window | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 61,311 | 0.18 | 96.59 | 99.39 | 8.76 | 0.01 | 1.40 | 0.66 | 0.74 | 97.19 | 0.9969 | 1.44 | 1.79 | 0.30 | - | - | - |
| GRCh38-only graph | 60,986 | 0.18 | 96.42 | 99.40 | 2.79 | 0.03 | 2.63 | 0.63 | 3.36 | 93.36 | 0.9969 | 1.48 | 1.82 | 0.19 | - | - | - |
| HPRC v1.1 d9 | 61,561 | 1.58 | 74.27 | 97.16 | 1.96 | 0.17 | 3.49 | 0.50 | 2.80 | 93.04 | 0.9989 | 0.46 | 1.61 | 0.10 | 94.13 | 7.28 | 98.46 |
| HPRC v2.1 d46 | 61,287 | 1.64 | 71.95 | 97.20 | 2.12 | 0.00 | 3.80 | 0.88 | 13.26 | 82.07 | 0.9990 | 0.45 | 1.63 | 0.18 | 93.71 | 4.85 | 98.46 |

Same reads, both arms place them in the window (B relative to A). The last four columns count the reads one arm places in the window that the other places elsewhere or leaves unmapped:

| B vs A | both in window | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both | A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 60,488 | 2.17 | 1.18 | 0.24 | 0.38 | 0.18 | 1.59 | 1.41 | 0.00 | 2.68 | 1.91 | 99.79 | 100.00 | 1 | 495 | 3 | 1,052 | 21 |
| HPRC v2.1 d46 vs GRCh38-only graph | 60,312 | 0.87 | 15.14 | 0.09 | 0.91 | 0.18 | 1.65 | 1.47 | 0.00 | 2.63 | 1.91 | 99.75 | 100.00 | 0 | 674 | 0 | 848 | 127 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 60,831 | 1.23 | 14.87 | 0.11 | 0.84 | 1.59 | 1.65 | 0.19 | 0.13 | 1.87 | 1.91 | 99.80 | 100.00 | 0 | 730 | 0 | 345 | 111 |
| GRCh38-only graph vs GRCh38 BAM | 60,720 | 0.56 | 5.67 | 0.09 | 1.69 | 0.18 | 0.18 | 0.00 | 0.00 | 8.28 | 2.68 | 99.71 | 100.00 | 1 | 451 | 140 | 266 | 0 |
| HPRC v1.1 d9 vs GRCh38 BAM | 60,448 | 0.57 | 5.00 | 0.16 | 1.71 | 0.18 | 1.59 | 1.41 | 0.00 | 8.22 | 1.89 | 99.65 | 99.99 | 7 | 758 | 105 | 1,113 | 0 |
| HPRC v2.1 d46 vs GRCh38 BAM | 60,371 | 0.40 | 16.29 | 0.01 | 2.24 | 0.18 | 1.65 | 1.47 | 0.00 | 8.23 | 1.97 | 99.62 | 100.00 | 1 | 940 | 0 | 915 | 1 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 851 | 96.12 | 12.0 | 13.75 | 86.25 | 0.35 | 0.00 | 10.0 | 100.00 | 7.76 | 0.82 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 887 | 96.05 | 12.0 | 13.87 | 86.13 | 0.23 | 0.00 | 10.0 | 100.00 | 0.56 | 52.99 | 0 | n/a |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 117 | 96.58 | 16.0 | 43.59 | 56.41 | 0.00 | 0.00 | 1.0 | 100.00 | 1.71 | 49.57 | 81 | 97.53 |
| GRCh38-only graph vs GRCh38 BAM | 0 |  |  |  |  |  |  |  |  |  |  |  |  |
| HPRC v1.1 d9 vs GRCh38 BAM | 851 | 96.12 | 12.0 | 13.63 | 86.02 | 2.35 | 0.00 | 10.0 | 100.00 | 0.00 | 16.57 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38 BAM | 887 | 96.05 | 12.0 | 13.75 | 85.91 | 2.03 | 0.00 | 10.0 | 100.00 | 0.00 | 53.89 | 0 | n/a |

By window kind and region (selected):

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | 4,453 | GRCh38 BAM | 79.59 | 20.41 | 0.00 | 0.20 | 96.42 | 99.10 | 15.29 | 0.11 | 84.85 | 0.9963 |
| windows: segdup | 4,453 | GRCh38-only graph | 78.98 | 19.99 | 1.03 | 0.20 | 96.19 | 99.06 | 7.85 | 0.14 | 68.44 | 0.9963 |
| windows: segdup | 4,453 | HPRC v1.1 d9 | 90.16 | 9.18 | 0.65 | 1.34 | 78.08 | 97.33 | 4.16 | 0.30 | 58.70 | 0.9987 |
| windows: segdup | 4,453 | HPRC v2.1 d46 | 88.77 | 11.23 | 0.00 | 1.47 | 77.66 | 97.37 | 4.91 | 0.00 | 45.43 | 0.9987 |
| SMaHT easy | 55,625 | GRCh38 BAM | 99.94 | 0.06 | 0.00 | 0.19 | 96.44 | 99.40 | 7.60 | 0.00 | 100.00 | 0.9970 |
| SMaHT easy | 55,625 | GRCh38-only graph | 99.94 | 0.05 | 0.01 | 0.19 | 96.29 | 99.41 | 2.16 | 0.00 | 96.80 | 0.9970 |
| SMaHT easy | 55,625 | HPRC v1.1 d9 | 100.00 | 0.00 | 0.00 | 1.62 | 73.64 | 97.12 | 1.61 | 0.00 | 97.29 | 0.9990 |
| SMaHT easy | 55,625 | HPRC v2.1 d46 | 99.97 | 0.03 | 0.00 | 1.68 | 71.12 | 97.16 | 1.60 | 0.00 | 86.43 | 0.9990 |
| SMaHT extreme | 6,880 | GRCh38 BAM | 79.08 | 20.90 | 0.01 | 0.13 | 98.14 | 99.30 | 19.85 | 0.13 | 68.54 | 0.9953 |
| SMaHT extreme | 6,880 | GRCh38-only graph | 74.62 | 22.69 | 2.69 | 0.14 | 97.92 | 99.24 | 9.41 | 0.33 | 56.02 | 0.9958 |
| SMaHT extreme | 6,880 | HPRC v1.1 d9 | 82.46 | 15.35 | 2.19 | 1.16 | 80.73 | 97.55 | 5.41 | 1.85 | 51.21 | 0.9985 |
| SMaHT extreme | 6,880 | HPRC v2.1 d46 | 78.63 | 21.37 | 0.00 | 1.22 | 80.74 | 97.63 | 7.41 | 0.00 | 37.25 | 0.9985 |

Pairs by region (B relative to A):

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 3,304 | 1.03 | 5.24 | 1.18 | 0.00 | 99.55 | 6.06 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 3,270 | 0.70 | 9.45 | 1.31 | 0.00 | 99.08 | 7.02 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 3,720 | 1.59 | 6.67 | 0.13 | 0.03 | 99.30 | 7.35 | 0.00 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 3,410 | 0.21 | 17.30 | 0.00 | 0.00 | 99.91 | 2.88 | 0.90 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 3,280 | 0.88 | 20.76 | 1.19 | 0.00 | 99.51 | 6.94 | 0.51 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 3,279 | 0.09 | 25.53 | 1.31 | 0.00 | 99.12 | 7.48 | 0.00 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 55,591 | 0.03 | 0.00 | 1.43 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 55,587 | 0.03 | 0.05 | 1.49 | 0.00 | 100.00 | 0.01 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 55,607 | 0.00 | 0.05 | 0.20 | 0.13 | 100.00 | 0.03 | 0.00 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 55,585 | 0.00 | 0.09 | 0.00 | 0.00 | 100.00 | 0.00 | 0.01 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 55,590 | 0.00 | 0.07 | 1.43 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 55,588 | 0.00 | 0.11 | 1.49 | 0.00 | 100.00 | 0.00 | 0.00 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 4,636 | 2.74 | 4.94 | 1.06 | 0.00 | 97.20 | 9.64 | 0.06 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 4,464 | 0.69 | 11.65 | 1.12 | 0.00 | 96.57 | 13.05 | 0.00 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 4,960 | 1.25 | 9.72 | 0.14 | 0.12 | 97.56 | 12.57 | 0.00 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 4,874 | 1.07 | 19.96 | 0.00 | 0.00 | 96.39 | 8.29 | 2.13 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 4,594 | 2.05 | 21.70 | 1.07 | 0.00 | 95.36 | 13.93 | 1.64 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 4,515 | 0.13 | 28.46 | 1.11 | 0.00 | 94.88 | 17.02 | 0.00 |

### 6. COLO829T ONT

**G sample (1,444,965 reads, 10 % of 14,448,438), graph arms.**

Whole GAM (every record, exact):

| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 14,448,438 | 1,954,638 | 13.53 | 0.88 | 2.72 | 93.85 | 1,444,965 |
| HPRC v1.1 d9 | 14,448,438 | 1,911,847 | 13.23 | 1.10 | 3.19 | 93.46 | 1,444,965 |
| HPRC v2.1 d46 | 14,448,438 | 1,968,240 | 13.62 | 0.38 | 2.80 | 89.57 | 1,444,965 |

Per arm, the same reads:

| arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 13.53 | 0.00 | 92.50 | 93.27 | 99.03 | 0.88 | 2.44 | 1.48 | 1.33 | 93.86 | 0.9887 | 9.16 | 8.59 | 2.82 | - | - | - |
| HPRC v1.1 d9 | 13.23 | 0.01 | 89.77 | 92.82 | 99.03 | 1.10 | 2.73 | 1.42 | 1.26 | 93.48 | 0.9910 | 8.36 | 8.55 | 2.75 | 73.15 | 8.58 | 94.47 |
| HPRC v2.1 d46 | 13.61 | 0.01 | 90.36 | 93.48 | 99.04 | 0.38 | 3.51 | 3.13 | 3.43 | 89.54 | 0.9908 | 8.92 | 9.12 | 2.85 | 73.51 | 9.37 | 90.10 |

Same reads, B relative to A:

| B vs A | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,246,241 | 2.02 | 2.81 | 0.54 | 0.96 | 0.00 | 0.01 | 0.00 | 0.00 | 85.63 | 85.92 | 96.80 | 99.63 | 4,356 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,226,010 | 2.75 | 7.57 | 1.01 | 3.63 | 0.00 | 0.01 | 0.00 | 0.00 | 85.63 | 85.56 | 95.01 | 98.57 | 16,164 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,229,957 | 3.16 | 7.05 | 1.14 | 3.33 | 0.01 | 0.01 | 0.00 | 0.00 | 85.92 | 85.56 | 95.20 | 98.60 | 15,892 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 39 | 87.18 | 1.0 | 64.10 | 35.90 | 0.00 | 0.00 | 1.0 | 97.44 | 7.69 | 0.00 | 9 | 22.22 |
| HPRC v2.1 d46 vs GRCh38-only graph | 40 | 87.50 | 1.0 | 65.00 | 35.00 | 2.50 | 0.00 | 1.0 | 97.50 | 2.50 | 2.50 | 8 | 12.50 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 4 | 50.00 | 2.0 | 50.00 | 25.00 | 25.00 | 0.00 | 1.0 | 50.00 | 0.00 | 25.00 | 2 | 50.00 |

Region anchor: linear 1,249,488, d9 7,096, d46 14,331, none 174,050. Reads unmapped only on the GRCh38-only graph number
7,529 against d9 and 22,295 against d46; reads unmapped only on d9 / d46 number 3,247 / 23,478 (JSON `G.pairs`). By
region (selected):

| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | 1,085,973 | GRCh38-only graph | 0.73 | 0.00 | 93.33 | 94.75 | 99.04 | 0.02 | 99.16 | 0.9891 | - |
| SMaHT easy | 1,085,973 | HPRC v1.1 d9 | 0.43 | 0.01 | 90.72 | 94.34 | 99.04 | 0.03 | 98.94 | 0.9913 | 75.88 |
| SMaHT easy | 1,085,973 | HPRC v2.1 d46 | 0.75 | 0.01 | 90.70 | 94.31 | 99.04 | 0.02 | 95.88 | 0.9913 | 76.01 |
| SMaHT difficult | 33,374 | GRCh38-only graph | 2.06 | 0.01 | 75.42 | 74.16 | 98.95 | 6.29 | 88.82 | 0.9901 | - |
| SMaHT difficult | 33,374 | HPRC v1.1 d9 | 1.27 | 0.03 | 67.43 | 71.89 | 98.95 | 7.16 | 88.06 | 0.9928 | 49.03 |
| SMaHT difficult | 33,374 | HPRC v2.1 d46 | 9.86 | 0.03 | 69.54 | 75.89 | 98.93 | 1.24 | 74.60 | 0.9923 | 51.48 |
| SMaHT extreme | 129,271 | GRCh38-only graph | 9.20 | 0.01 | 89.07 | 85.18 | 98.89 | 7.08 | 53.44 | 0.9837 | - |
| SMaHT extreme | 129,271 | HPRC v1.1 d9 | 9.05 | 0.01 | 86.47 | 84.84 | 98.91 | 8.94 | 51.94 | 0.9872 | 60.18 |
| SMaHT extreme | 129,271 | HPRC v2.1 d46 | 9.53 | 0.01 | 91.40 | 90.35 | 98.96 | 2.90 | 43.54 | 0.9823 | 60.42 |
| segdup | 81,169 | GRCh38-only graph | 6.25 | 0.00 | 93.60 | 92.52 | 98.84 | 1.41 | 58.05 | 0.9852 | - |
| segdup | 81,169 | HPRC v1.1 d9 | 5.70 | 0.00 | 91.29 | 92.04 | 98.83 | 2.81 | 55.80 | 0.9887 | 73.38 |
| segdup | 81,169 | HPRC v2.1 d46 | 2.56 | 0.00 | 93.54 | 93.98 | 98.93 | 1.89 | 49.42 | 0.9866 | 73.09 |

Pairs by region (B relative to A):

| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 0.05 | 0.18 | 0.00 | 0.00 | 99.93 | 107 | 92.86 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 0.23 | 1.82 | 0.00 | 0.00 | 99.08 | 4,736 | 89.66 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 0.34 | 1.78 | 0.00 | 0.00 | 99.08 | 4,677 | 50.00 |
| SMaHT difficult | HPRC v1.1 d9 vs GRCh38-only graph | 0.96 | 1.81 | 0.02 | 0.00 | 91.64 | 101 | 66.67 |
| SMaHT difficult | HPRC v2.1 d46 vs GRCh38-only graph | 1.23 | 13.09 | 0.02 | 0.00 | 87.12 | 1,579 | 66.67 |
| SMaHT difficult | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1.50 | 12.82 | 0.00 | 0.00 | 87.61 | 1,563 | n/a |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 4.63 | 7.41 | 0.00 | 0.01 | 72.45 | 3,655 | 100.00 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 7.51 | 17.30 | 0.00 | 0.01 | 62.82 | 7,753 | 100.00 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 7.72 | 14.63 | 0.00 | 0.00 | 64.75 | 7,399 | 50.00 |
| segdup | HPRC v1.1 d9 vs GRCh38-only graph | 3.58 | 7.14 | 0.00 | 0.00 | 80.25 | 1,367 | 100.00 |
| segdup | HPRC v2.1 d46 vs GRCh38-only graph | 6.29 | 13.77 | 0.00 | 0.00 | 66.30 | 6,240 | 100.00 |
| segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 6.80 | 10.79 | 0.00 | 0.00 | 67.45 | 5,994 | 0.00 |

**W sample (148 windows, flank 25 kb; 33,011 reads: random 6,224, random_extra 24,248, segdup 2,539), all four arms.**
The BAM has 1,070,571 unplaced unmapped records in the whole file; supplementary records for 2.01 % of its in-window
reads.

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38 BAM | window | 92.45 | 6.31 | 1.24 |
| GRCh38-only graph | genome | 87.84 | 6.84 | 5.32 |
| HPRC v1.1 d9 | genome | 88.75 | 6.16 | 5.10 |
| HPRC v2.1 d46 | genome | 88.53 | 4.51 | 6.96 |

Reads the arm places in the window:

| arm | reads in window | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38 BAM | 30,519 | 0.01 | 93.61 | 94.77 | 99.63 | 0.67 | 2.22 | 1.53 | 2.15 | 93.42 | 0.9884 | 12.94 | 12.16 | 3.76 | - | - | - |
| GRCh38-only graph | 28,997 | 0.01 | 92.85 | 93.79 | 99.07 | 0.71 | 2.17 | 1.22 | 1.31 | 94.59 | 0.9889 | 8.84 | 8.99 | 2.59 | - | - | - |
| HPRC v1.1 d9 | 29,296 | 0.01 | 90.04 | 93.19 | 99.07 | 1.13 | 2.72 | 1.32 | 1.22 | 93.60 | 0.9913 | 8.04 | 8.93 | 2.28 | 76.33 | 7.93 | 94.73 |
| HPRC v2.1 d46 | 29,226 | 0.01 | 90.69 | 93.78 | 99.08 | 0.22 | 3.25 | 2.82 | 3.22 | 90.50 | 0.9910 | 8.59 | 9.34 | 2.24 | 75.66 | 5.71 | 94.11 |

Same reads, both arms place them in the window (B relative to A). The last four columns count the reads one arm places in the window that the other places elsewhere or leaves unmapped:

| B vs A | both in window | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both | A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 28,388 | 1.60 | 1.80 | 0.32 | 0.57 | 0.01 | 0.01 | 0.00 | 0.00 | 99.06 | 99.07 | 99.54 | 99.96 | 12 | 593 | 16 | 824 | 84 |
| HPRC v2.1 d46 vs GRCh38-only graph | 27,550 | 1.71 | 4.92 | 0.51 | 2.15 | 0.01 | 0.01 | 0.00 | 0.00 | 99.05 | 99.06 | 99.38 | 99.88 | 33 | 977 | 470 | 1,473 | 203 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 27,759 | 2.28 | 4.81 | 0.73 | 2.08 | 0.01 | 0.01 | 0.00 | 0.00 | 99.05 | 99.05 | 99.42 | 99.86 | 36 | 1,020 | 517 | 1,335 | 132 |
| GRCh38-only graph vs GRCh38 BAM | 28,417 | 2.40 | 2.85 | 0.62 | 1.06 | 0.01 | 0.01 | 0.00 | 0.00 | 99.62 | 99.06 | 99.55 | 99.99 | 4 | 377 | 1,725 | 372 | 208 |
| HPRC v1.1 d9 vs GRCh38 BAM | 28,355 | 2.30 | 2.83 | 0.64 | 1.31 | 0.01 | 0.01 | 0.00 | 0.00 | 99.62 | 99.06 | 99.48 | 99.97 | 7 | 517 | 1,647 | 704 | 237 |
| HPRC v2.1 d46 vs GRCh38 BAM | 27,662 | 1.91 | 4.99 | 0.47 | 2.36 | 0.01 | 0.01 | 0.00 | 0.00 | 99.62 | 99.05 | 99.33 | 99.86 | 36 | 1,017 | 1,840 | 1,513 | 51 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38-only graph | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 0 |  |  |  |  |  |  |  |  |  |  |  |  |
| GRCh38-only graph vs GRCh38 BAM | 0 |  |  |  |  |  |  |  |  |  |  |  |  |
| HPRC v1.1 d9 vs GRCh38 BAM | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |
| HPRC v2.1 d46 vs GRCh38 BAM | 1 | 100.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0.00 | 1.0 | 100.00 | 0.00 | 0.00 | 0 | n/a |

By window kind and region (selected):

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | 2,539 | GRCh38 BAM | 72.71 | 26.51 | 0.79 | 0.00 | 93.82 | 95.40 | 99.89 | 3.90 | 72.16 | 0.9861 |
| windows: segdup | 2,539 | GRCh38-only graph | 68.81 | 27.14 | 4.06 | 0.00 | 93.36 | 94.33 | 99.60 | 0.57 | 71.55 | 0.9870 |
| windows: segdup | 2,539 | HPRC v1.1 d9 | 75.50 | 20.72 | 3.78 | 0.00 | 92.07 | 94.16 | 99.48 | 1.62 | 63.07 | 0.9901 |
| windows: segdup | 2,539 | HPRC v2.1 d46 | 73.18 | 21.82 | 5.00 | 0.00 | 93.43 | 94.67 | 99.25 | 1.02 | 59.96 | 0.9896 |
| SMaHT easy | 26,908 | GRCh38 BAM | 98.90 | 0.82 | 0.29 | 0.00 | 94.12 | 95.48 | 99.61 | 0.00 | 98.34 | 0.9886 |
| SMaHT easy | 26,908 | GRCh38-only graph | 93.68 | 1.11 | 5.21 | 0.00 | 93.76 | 95.13 | 99.05 | 0.03 | 99.18 | 0.9891 |
| SMaHT easy | 26,908 | HPRC v1.1 d9 | 94.26 | 0.80 | 4.94 | 0.01 | 90.98 | 94.41 | 99.05 | 0.26 | 98.71 | 0.9913 |
| SMaHT easy | 26,908 | HPRC v2.1 d46 | 93.48 | 0.88 | 5.64 | 0.01 | 91.09 | 94.65 | 99.05 | 0.01 | 96.23 | 0.9913 |
| SMaHT extreme | 5,097 | GRCh38 BAM | 62.66 | 32.80 | 4.53 | 0.06 | 93.39 | 92.39 | 99.84 | 6.45 | 54.82 | 0.9841 |
| SMaHT extreme | 5,097 | GRCh38-only graph | 59.17 | 34.26 | 6.57 | 0.07 | 90.19 | 87.23 | 99.34 | 4.87 | 57.89 | 0.9869 |
| SMaHT extreme | 5,097 | HPRC v1.1 d9 | 61.98 | 31.47 | 6.55 | 0.06 | 88.26 | 88.26 | 99.27 | 6.49 | 53.94 | 0.9902 |
| SMaHT extreme | 5,097 | HPRC v2.1 d46 | 64.90 | 22.78 | 12.32 | 0.06 | 92.47 | 91.51 | 99.33 | 1.72 | 51.06 | 0.9863 |

Pairs by region (B relative to A):

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 1,567 | 1.53 | 5.42 | 0.00 | 0.00 | 98.34 | 10.30 | 0.00 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 1,379 | 2.76 | 7.18 | 0.00 | 0.00 | 96.23 | 19.58 | 1.49 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,498 | 5.01 | 5.41 | 0.00 | 0.00 | 96.80 | 20.61 | 1.25 |
| windows: segdup | GRCh38-only graph vs GRCh38 BAM | 1,648 | 1.58 | 7.65 | 0.00 | 0.00 | 98.79 | 5.47 | 5.25 |
| windows: segdup | HPRC v1.1 d9 vs GRCh38 BAM | 1,568 | 1.15 | 11.03 | 0.00 | 0.00 | 98.72 | 10.08 | 4.98 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38 BAM | 1,378 | 1.67 | 11.10 | 0.00 | 0.00 | 96.23 | 20.10 | 5.25 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 25,190 | 0.04 | 0.12 | 0.00 | 0.00 | 100.00 | 0.06 | 0.01 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 24,811 | 0.19 | 1.39 | 0.00 | 0.00 | 99.94 | 0.79 | 0.78 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 24,884 | 0.28 | 1.39 | 0.00 | 0.00 | 99.94 | 0.86 | 1.03 |
| SMaHT easy | GRCh38-only graph vs GRCh38 BAM | 25,171 | 0.06 | 0.35 | 0.00 | 0.00 | 99.98 | 0.17 | 5.24 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38 BAM | 25,237 | 0.07 | 0.48 | 0.00 | 0.00 | 99.99 | 0.17 | 4.99 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38 BAM | 24,954 | 0.06 | 1.51 | 0.00 | 0.00 | 99.94 | 0.79 | 5.43 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 2,489 | 3.13 | 4.90 | 0.00 | 0.00 | 94.82 | 17.01 | 0.46 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 2,117 | 4.20 | 10.06 | 0.00 | 0.00 | 92.73 | 23.28 | 6.53 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,249 | 5.56 | 8.71 | 0.00 | 0.00 | 93.55 | 23.20 | 5.60 |
| SMaHT extreme | GRCh38-only graph vs GRCh38 BAM | 2,567 | 5.03 | 8.06 | 0.00 | 0.00 | 95.13 | 9.80 | 9.83 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38 BAM | 2,435 | 5.34 | 9.94 | 0.00 | 0.00 | 94.00 | 14.21 | 9.55 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38 BAM | 2,097 | 4.43 | 11.73 | 0.00 | 0.00 | 91.99 | 23.17 | 11.18 |

### 7. COLO829T Illumina, W sample (48 windows, flank 1 kb; 1,267,256 reads: random 972,840, segdup 294,416), graph arms only

| arm | scope | in the window % | placed elsewhere % | unmapped % |
|---|---|---:|---:|---:|
| GRCh38-only graph | window | 92.35 | 7.38 | 0.27 |
| HPRC v1.1 d9 | window | 92.23 | 7.76 | 0.01 |
| HPRC v2.1 d46 | window | 92.06 | 7.67 | 0.27 |

Unmapped records of the whole file (GAM heads; 2,479,603,422 reads, `compare.py` `TOTAL_READS`): GRCh38-only graph
24,512,000 (0.99 %), HPRC v1.1 d9 15,568,000 (0.63 %), HPRC v2.1 d46 23,242,000 (0.94 %).

Reads the arm places in the window:

| arm | reads in window | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GRCh38-only graph | 1,170,285 | 70.72 | 26.87 | 2.84 | 2.63 | 3.98 | 3.69 | 2.96 | 3.75 | 85.61 | 1.0000 | 3.95 | 0.21 | 0.68 | - | - | - |
| HPRC v1.1 d9 | 1,168,846 | 81.25 | 17.08 | 0.58 | 2.64 | 4.83 | 3.88 | 2.79 | 3.38 | 85.11 | 1.0000 | 2.97 | 0.04 | 0.88 | 14.54 | 4.89 | 98.98 |
| HPRC v2.1 d46 | 1,166,614 | 81.61 | 17.00 | 0.56 | 2.36 | 4.74 | 3.60 | 2.75 | 3.45 | 85.45 | 1.0000 | 2.94 | 0.04 | 0.61 | 14.42 | 3.54 | 98.96 |

Same reads, both arms place them in the window (B relative to A). The last four columns count the reads one arm places in the window that the other places elsewhere or leaves unmapped:

| B vs A | both in window | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both | A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 1,102,189 | 1.57 | 3.15 | 0.39 | 1.07 | 71.02 | 81.66 | 10.64 | 0.00 | 2.49 | 2.29 | 99.58 | 100.00 | 16 | 68,027 | 69 | 63,261 | 3,396 |
| HPRC v2.1 d46 vs GRCh38-only graph | 1,107,463 | 1.57 | 3.20 | 0.38 | 0.91 | 71.01 | 81.73 | 10.72 | 0.00 | 2.49 | 2.30 | 99.46 | 100.00 | 22 | 62,585 | 237 | 58,890 | 261 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 1,114,035 | 0.98 | 0.79 | 0.42 | 0.24 | 81.66 | 81.78 | 0.20 | 0.07 | 2.29 | 2.30 | 99.69 | 100.00 | 9 | 51,473 | 3,338 | 52,552 | 27 |

Reads that became perfect / stopped being perfect:

| B vs A | perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | placement concordant % | MAPQ up % | MAPQ down % | perfect in A, not in B | of these on off-ref in B % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HPRC v1.1 d9 vs GRCh38-only graph | 117,286 | 94.11 | 1.0 | 85.76 | 13.63 | 0.70 | 0.00 | 1.0 | 99.60 | 5.17 | 1.93 | 7 | 57.14 |
| HPRC v2.1 d46 vs GRCh38-only graph | 118,738 | 94.02 | 1.0 | 85.76 | 13.63 | 0.72 | 0.00 | 1.0 | 99.54 | 5.25 | 1.55 | 2 | 0.00 |
| HPRC v2.1 d46 vs HPRC v1.1 d9 | 2,206 | 92.20 | 1.0 | 83.50 | 15.28 | 1.36 | 0.00 | 1.0 | 99.14 | 10.47 | 0.73 | 783 | 21.97 |

Reads with MAPQ >= 30 in both arms (JSON `W.pairs`): 1,026,464 (d9 vs linear), 1,030,949 (d46 vs linear), 1,029,601
(d46 vs d9).

By window kind and region (selected):

| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: random | 972,840 | GRCh38-only graph | 96.96 | 2.76 | 0.28 | 71.73 | 25.82 | 2.73 | 2.54 | 1.77 | 93.53 | 1.0000 |
| windows: random | 972,840 | HPRC v1.1 d9 | 96.97 | 3.02 | 0.01 | 82.34 | 16.01 | 0.38 | 2.57 | 2.03 | 93.55 | 1.0000 |
| windows: random | 972,840 | HPRC v2.1 d46 | 96.35 | 3.37 | 0.28 | 82.78 | 15.86 | 0.35 | 2.29 | 1.91 | 94.03 | 1.0000 |
| windows: segdup | 294,416 | GRCh38-only graph | 77.11 | 22.66 | 0.23 | 66.52 | 31.26 | 3.31 | 3.00 | 13.17 | 52.72 | 1.0000 |
| windows: segdup | 294,416 | HPRC v1.1 d9 | 76.58 | 23.41 | 0.01 | 76.70 | 21.53 | 1.40 | 2.94 | 16.56 | 49.84 | 1.0000 |
| windows: segdup | 294,416 | HPRC v2.1 d46 | 77.89 | 21.88 | 0.24 | 76.83 | 21.67 | 1.46 | 2.62 | 16.34 | 50.40 | 1.0000 |
| SMaHT easy | 795,796 | GRCh38-only graph | 99.67 | 0.07 | 0.26 | 73.67 | 24.31 | 2.06 | 1.82 | 0.06 | 98.10 | 1.0000 |
| SMaHT easy | 795,796 | HPRC v1.1 d9 | 99.87 | 0.13 | 0.00 | 83.82 | 14.81 | 0.26 | 1.94 | 0.09 | 97.88 | 1.0000 |
| SMaHT easy | 795,796 | HPRC v2.1 d46 | 99.66 | 0.08 | 0.26 | 84.17 | 14.73 | 0.23 | 1.69 | 0.09 | 98.01 | 1.0000 |
| SMaHT difficult | 126,883 | GRCh38-only graph | 99.24 | 0.33 | 0.43 | 63.11 | 32.50 | 5.62 | 5.60 | 0.65 | 93.83 | 1.0000 |
| SMaHT difficult | 126,883 | HPRC v1.1 d9 | 99.31 | 0.68 | 0.00 | 75.26 | 21.85 | 0.90 | 5.39 | 0.60 | 94.05 | 1.0000 |
| SMaHT difficult | 126,883 | HPRC v2.1 d46 | 98.98 | 0.59 | 0.42 | 75.84 | 21.69 | 0.83 | 4.94 | 0.51 | 94.32 | 1.0000 |
| SMaHT extreme | 344,577 | GRCh38-only graph | 72.90 | 26.87 | 0.23 | 65.22 | 32.16 | 3.92 | 3.68 | 18.06 | 42.08 | 1.0000 |
| SMaHT extreme | 344,577 | HPRC v1.1 d9 | 72.00 | 27.97 | 0.02 | 76.06 | 21.91 | 1.45 | 3.49 | 22.16 | 39.70 | 1.0000 |
| SMaHT extreme | 344,577 | HPRC v2.1 d46 | 71.96 | 27.80 | 0.24 | 76.34 | 21.89 | 1.49 | 3.18 | 21.78 | 40.79 | 1.0000 |

Pairs by region (B relative to A):

| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| windows: segdup | HPRC v1.1 d9 vs GRCh38-only graph | 179,077 | 1.45 | 5.08 | 9.47 | 0.00 | 98.17 | 21.11 | 0.01 |
| windows: segdup | HPRC v2.1 d46 vs GRCh38-only graph | 186,135 | 1.49 | 3.88 | 9.26 | 0.00 | 97.39 | 17.99 | 0.02 |
| windows: segdup | HPRC v2.1 d46 vs HPRC v1.1 d9 | 190,960 | 1.83 | 0.66 | 0.11 | 0.16 | 98.73 | 15.01 | 0.30 |
| SMaHT easy | HPRC v1.1 d9 vs GRCh38-only graph | 792,129 | 0.06 | 0.18 | 10.38 | 0.00 | 100.00 | 0.13 | 0.00 |
| SMaHT easy | HPRC v2.1 d46 vs GRCh38-only graph | 792,726 | 0.05 | 0.13 | 10.50 | 0.00 | 100.00 | 0.04 | 0.01 |
| SMaHT easy | HPRC v2.1 d46 vs HPRC v1.1 d9 | 792,275 | 0.13 | 0.09 | 0.17 | 0.05 | 100.00 | 0.05 | 0.26 |
| SMaHT extreme | HPRC v1.1 d9 vs GRCh38-only graph | 184,936 | 1.98 | 5.43 | 10.63 | 0.00 | 97.54 | 26.35 | 0.03 |
| SMaHT extreme | HPRC v2.1 d46 vs GRCh38-only graph | 189,412 | 1.92 | 4.58 | 10.43 | 0.00 | 96.83 | 24.55 | 0.04 |
| SMaHT extreme | HPRC v2.1 d46 vs HPRC v1.1 d9 | 196,490 | 1.76 | 0.91 | 0.22 | 0.16 | 98.28 | 20.50 | 0.30 |

## Interpretation

1. **The short-read perfect-read gain is real, and it comes from graph content, not from the aligner.**
   - HG008T Illumina, the same reads in the same windows: 69.18 % perfect on the GRCh38 BAM and 69.18 % on the
     GRCh38-only graph (0.00 % gained or lost), against 79.61 % on d9 and 79.65 % on d46 (+10.34 / +10.41 points vs the
     BAM; +10.24 / +10.34 vs the GRCh38-only graph). COLO829T Illumina (no BAM): 71.02 % on the GRCh38-only graph,
     81.66 % on d9 and 81.73 % on d46 (+10.64 / +10.72 points).
   - d46 adds little over d9 for short reads: perfect gained / lost 0.29 / 0.17 % (HG008T), 0.20 / 0.07 % (COLO829T).
   - The reads that became perfect were mostly one-mismatch reads: in the six HPRC-vs-GRCh38 comparisons 84.74-86.09 %
     had mismatches only (median 1 mismatched base), 94.02-94.16 % now pass over off-reference nodes, and 99.54-99.74 %
     keep their GRCh38 placement.
   - For long reads, "perfect" is the wrong scale: almost every HiFi / fiber-seq read carries an indel (97.12-99.43 % of
     mapped reads, G), and ONT reads more so. Mismatched bases per kb show the effect instead (GRCh38-only graph -> d9,
     G genome-wide; SMaHT easy from the JSON in brackets): HG008T PacBio 1.57 -> 0.63 (1.020 -> 0.155), COLO829T
     fiber-seq 1.80 -> 0.84 (1.073 -> 0.198), HG008T ONT 5.22 -> 4.26 (4.425 -> 3.582), COLO829T ONT 9.16 -> 8.36
     (8.261 -> 7.576). Indel events per kb change less (1.58 -> 1.39, 1.82 -> 1.62, 4.27 -> 4.12, 8.59 -> 8.55).
   - On d9, 92.08 % (PacBio), 92.54 % (fiber-seq), 81.48 % (HG008T ONT) and 73.15 % (COLO829T ONT) of the mapped reads
     visit at least one off-reference node, and 98.01 / 97.90 / 96.42 / 94.47 % of those visits are edit-free.
   - This answers R2.7 with counts and proportions. The "more perfect reads" statement holds for short reads (about +10
     points on both tumors). For long reads it should be put as fewer mismatches per kb rather than more perfect reads.
2. **Reads that improve stay where they were, and confident placements agree.**
   - Reads with MAPQ >= 30 on both the GRCh38-only graph and d9 (G): concordant 99.98 % for HG008T PacBio (307 discordant
     of 1,316,930), 99.97 % for COLO829T fiber-seq (747 of 2,497,017), 99.68 % for HG008T ONT (1,019 of 319,630) and
     99.63 % for COLO829T ONT (4,356 of 1,175,820). Illumina windows: 19 of 1,192,066 (HG008T) and 16 of 1,026,464
     (COLO829T).
   - The same for d46 vs the GRCh38-only graph: 99.98 % (291), 99.97 % (757), 99.38 % (1,953 of 315,825) and 98.57 %
     (16,164 of 1,134,188); Illumina 19 of 1,195,562 and 22 of 1,030,949.
   - Reads that became perfect keep their locus: 99.72 % (PacBio d9), 99.56 % (fiber-seq d9), 98.67 % (HG008T ONT d9).
   - So in the real data the higher alignment identity is not bought by moving reads to other loci (R3.6); ONT on d46 is
     the exception (point 6). Whether the placements are *correct* is measured in E11b.
3. **MAPQ changes little between the GRCh38-only graph and d9; d46 MAPQ is scaled differently for HiFi-type reads.**
   - d9 vs the GRCh38-only graph, MAPQ >= 30 gained / lost: PacBio 0.20 / 0.43 %, fiber-seq 0.25 / 0.44 %, HG008T ONT
     0.47 / 0.72 %, COLO829T ONT 0.54 / 0.96 % (G); Illumina windows 0.43 / 0.92 % (HG008T), 0.39 / 1.07 % (COLO829T).
   - MAPQ 60 share of mapped records (whole GAM; GRCh38-only / d9 / d46): PacBio 95.88 / 95.88 / 83.32 %, fiber-seq
     92.38 / 92.63 / 76.87 %, HG008T ONT 95.18 / 95.02 / 94.16 %, COLO829T ONT 93.85 / 93.46 / 89.57 %. Illumina in-window
     reads: 86.89 / 86.24 / 86.29 % (HG008T), 85.61 / 85.11 / 85.45 % (COLO829T).
   - The d46 drop is largest for PacBio and fiber-seq and appears in SMaHT easy regions too (MAPQ 60 99.14 -> 89.45 % and
     96.78 -> 86.15 %), where placement concordance with the GRCh38-only graph stays at 99.98 %. It is a change of the
     MAPQ scale (different vg version / graph), not misplacement; the HiFi simulation of E11b shows the same shift without
     extra wrong placements. MAPQ values are not directly comparable between d46 and the other arms.
4. **Difficult regions (R1.3).**
   - SMaHT extreme, placement concordance with the GRCh38-only graph (G): d9 82.74 % (PacBio), 82.11 % (fiber-seq),
     72.77 % (HG008T ONT), 72.45 % (COLO829T ONT); d46 78.56 / 76.99 / 70.38 / 62.82 %. More reads lose MAPQ >= 30 than
     gain it on d9: 6.39 vs 3.05 %, 6.26 vs 3.40 %, 6.92 vs 4.93 %, 7.41 vs 4.63 %; on d46 13.97 vs 1.45 %, 14.20 vs
     1.77 %, 14.12 vs 7.09 %, 17.30 vs 7.51 %.
   - Illumina segdup windows: d9 has more MAPQ 0 reads than the GRCh38-only graph (HG008T 15.04 vs 11.56 %, BAM 16.77 %;
     COLO829T 16.56 vs 13.17 %) and fewer at MAPQ 60 (54.28 vs 57.75 %; 49.84 vs 52.72 %). MAPQ >= 30 lost exceeds
     gained (d9 4.25 vs 1.54 % and 5.08 vs 1.45 %; d46 3.44 vs 1.53 % and 3.88 vs 1.49 %), while the perfect-read gain
     stays (+10.92 / +10.23 points HG008T, +9.47 / +9.26 COLO829T).
   - The pangenome therefore does not reduce mapping ambiguity in segdups for short reads; it slightly increases it,
     consistent with paralogous haplotypes in the graph.
   - Long reads in segdup windows: d9 places more reads in the window than the BAM or the GRCh38-only graph: PacBio
     90.19 % vs 82.64 / 80.34 %, fiber-seq 90.16 % vs 79.59 / 78.98 %, HG008T ONT 82.57 % vs 82.41 / 78.92 %, COLO829T
     ONT 75.50 % vs 72.71 / 68.81 %.
5. **d46 aligns almost every HiFi / fiber-seq read, but not more ONT or Illumina reads, and its additional placements
   are low-confidence.**
   - Whole-GAM unmapped (GRCh38-only / d9 / d46): PacBio 2.78 / 2.64 / 0.01 %, fiber-seq 5.62 / 5.43 / 0.03 %; HG008T ONT
     7.88 / 7.77 / 8.04 %, COLO829T ONT 13.53 / 13.23 / 13.62 %. Illumina whole files: 1.51 / 1.42 / 1.41 % (HG008T),
     0.99 / 0.63 / 0.94 % (COLO829T).
   - The rescued HiFi-type reads sit mostly in SMaHT extreme regions (unmapped 29.30 -> 0.03 % of 109,232 PacBio reads,
     32.76 -> 0.04 % of 229,856 fiber-seq reads). On d46 they align with low MAPQ (MAPQ 60 22.52 vs 48.75 % and 19.21
     vs 45.32 %), more soft clipping (23.55 vs 10.69 % and 29.72 vs 12.67 % of reads) and more mismatches (17.635 vs
     9.159 and 20.957 vs 11.369 bp per kb, JSON).
   - Genome-wide this lifts d46's mismatch rate above d9's (PacBio 1.72 vs 0.63, fiber-seq 2.19 vs 0.84 bp per kb), while
     the two are equal in SMaHT easy regions (0.150 vs 0.155, 0.192 vs 0.198).
   - A lower unmapped rate on d46 is therefore not a quality gain across the board. E11b's HiFi simulation (one-to-one
     regions only) finds no extra confident misplacement on d46, but it cannot test the regions where most rescued reads
     sit.
6. **ONT on d46: where HPRC v2.1 does worse in these statistics.**
   - d46 leaves more ONT reads unmapped than the GRCh38-only graph (whole GAM 8.04 vs 7.88 % for HG008T, 13.62 vs
     13.53 % for COLO829T), most visibly in SMaHT difficult regions (G: 10.57 vs 2.94 % and 9.86 vs 2.06 %).
   - Confident placements disagree more often: discordant with MAPQ >= 30 in both, d46 vs the GRCh38-only graph, 1,953
     (HG008T) and 16,164 (COLO829T) reads against 1,019 / 4,356 for d9; in SMaHT easy regions alone 328 and 4,736
     against 11 / 107 for d9.
   - In SMaHT easy regions d46 removes fewer mismatched ONT bases than d9 (HG008T 4.425 -> 3.769 vs 3.582; COLO829T
     8.261 -> 8.305 vs 7.576 bp per kb, JSON), and indel events per kb rise (4.303 -> 4.340; 8.662 -> 9.197).
   - Which placement is correct is not known: the ONT simulation of E11b was not aligned. The d46 GAMs were also made
     with a newer giraffe (vg 1.77), so the aligner version may explain part of this.
7. **The comparators' GRCh38 BAMs vs giraffe on the GRCh38-only graph** (same in-window reads).
   - Perfect rates are equal: PacBio 0.18 / 0.18 %, fiber-seq 0.18 / 0.18 %, HG008T Illumina 69.18 / 69.18 %, HG008T ONT
     1.09 / 1.08 %, COLO829T ONT 0.01 / 0.01 %. Placement concordance is 99.55-99.71 %.
   - The long-read BAMs soft-clip more reads (PacBio 5.40 vs 2.16 %, fiber-seq 8.28 vs 2.68 %, HG008T ONT 34.93 vs
     28.37 %, COLO829T ONT 99.62 vs 99.06 %), and the BAM rows count hard clips too.
   - Giraffe leaves more reads unmapped. Of the reads the BAM places in a window, the GRCh38-only graph leaves 104
     (PacBio), 140 (fiber-seq), 214 (HG008T ONT), 1,725 (COLO829T ONT) and 11,142 (HG008T Illumina) unmapped; for ONT
     that is 2.01 % (HG008T) and 5.24 % (COLO829T) even in SMaHT easy regions. Whole-GAM unmapped rates on the
     GRCh38-only graph are 2.78 % (PacBio), 5.62 % (fiber-seq), 7.88 % (HG008T ONT), 13.53 % (COLO829T ONT) and 1.51 %
     (HG008T Illumina).
8. **COLO829T ONT reads are clipped almost everywhere.** 99.03-99.04 % of mapped reads carry a soft clip on every graph
   arm and 99.63 % in the BAM, and perfect reads are 0.00-0.01 %. Perfect-read comparisons say nothing for this set;
   the per-kb rates apply. The cause (for example untrimmed adapter sequence) was not investigated.

## Caveats

- **Aligners and versions differ by arm**, so MAPQ and soft-clip comparisons across arms are descriptive.
  - The GRCh38 BAMs are the GIAB / SMaHT alignments. Their `@PG` header lines (`samtools view -H`, samtools 1.18,
    checked 2026-10-15) give the aligners of the HG008T BAMs: Illumina BWA-MEM 0.7.17-r1188 (`bwa mem`, then
    samtools sort / merge), PacBio pbmm2 1.10.0 (`--preset HIFI`; also ccs 8.0.0, jasmine 2.2.0 and hiphase upstream),
    ONT minimap2 2.26-r1175 (`-x map-ont`), all against the GIAB v3 GRCh38 no-alt analysis set. The COLO829T fiber-seq
    BAM is pbmm2 1.13.0 (name and header); the COLO829T ONT BAM header lists minimap2 2.26-r1175 for its alignment steps
    and also a minimap2 2.24-r1122 line of unclear origin (see `../../points/R2.7.md`).
  - The GRCh38-only graph GAMs were made with vg 1.66; the d46 GAMs with vg 1.77 on indexes rebuilt with vg 1.77 (project
    notes on the graph ablation). The HG008T d9 GAMs were made with vg giraffe 1.66.0 (production logs; see the evidence
    in `../../points/R3.1.md`), the same version as the GRCh38-only GAMs; the COLO829T d9 version is unverified.
- **"Perfect" is strict** (no mismatch, indel or soft clip over the full read). It is near zero for long reads by
  construction, so the per-kb rates are the better long-read measure.
- **The GAM placement has node granularity.** It is the contig holding most GRCh38-anchored bases, from the first to the
  last anchored node. A read that touches only off-reference nodes has no placement.
  - Concordance allows a 1 kb gap. Two placements in the same 500 kb window can still be "discordant" (e.g. the 22
    HG008T Illumina reads that became perfect on the GRCh38-only graph relative to the BAM, 4.55 % concordant).
- **The W design conditions on placement.** The universe is the set of reads that at least one arm places in a window,
  so it grows when an arm is added (HG008T Illumina 1,443,548 -> 1,465,754 reads, HG008T ONT 8,614 -> 9,000 between
  2026-10-13 and 2026-10-15).
  - `elsewhere` for a window-scanned GAM means the read visits none of the window's nodes. For the BAM it means the read
    was not seen within window ± flank, which would also count a read missing from the BAM.
  - Supplementary records are seen only inside the same fetch, so the BAM's multi-contig percentage is 0 by construction.
  - The Illumina BAM's whole-file count (2,165,824) covers only unplaced records (`fetch('*')`). Unmapped reads with a
    placed mate are seen only inside windows.
- **The windows are small.** The segdup windows are 24 × 100 kb, and some strata are thin (W "SMaHT difficult": 160
  reads for HG008T PacBio and for HG008T ONT, 281 for COLO829T fiber-seq; HG008T ONT W "SMaHT extreme" 1,311 reads).
- **One anchor per read for the strata.** G uses the linear placement first, so reads unmapped on the GRCh38-only graph
  but mapped on d46 are stratified by their d46 placement.
- **The regions are SMaHT v2 easy / difficult / extreme, segdups and the HG008T benchmark BED** (HG008T only), not GIAB
  v3.5 (not on the server; E10).
- **COLO829T fiber-seq G sample.** The whole-GAM scans count 2,784,787 sampled records but 2,687,506 distinct reads, and
  report 44,480 duplicate (name, mate) keys (JSON `G.whole_file.*.duplicate_keys`). The sampled distinct reads are
  unmapped at 3.37 / 3.23 / 0.01 % against 5.62 / 5.43 / 0.03 % of all records. Not investigated.
- **COLO829T Illumina d9** leaves only 0.01 % of the W universe unmapped, against 0.27 % on the other two graph arms
  (whole-file unmapped 0.63 % against 0.99 / 0.94 %). Not investigated.
- **No correctness against a known origin here** (E11b, Illumina and HiFi only), and no MAPQ calibration.
- MAPQ distributions are reported as bins in the tables; no figure has been made.
- **`results/tables.md` formatting bug.** In the W "Same reads, both arms place them in the window" tables, the header
  lacks a `|` after "discordant with MAPQ>=30 in both". The header row has one cell fewer than the data rows, so the last
  five header labels are shifted. The data rows are correct, and the tables above were relabelled from the JSON fields
  (`concordant_mapq30_both_pct`, `discordant_mapq30_both`, `A_in_B_elsewhere`, `A_in_B_unmapped`, `B_in_A_elsewhere`,
  `B_in_A_unmapped`). `compare.py` was not changed.
- `compare.py` rewrites `results/*.json` and `tables.md` on every run. The numbers here are from the 2026-10-15 run.

## Status

**DONE (2026-10-15).** Every scan job finished and `compare.py` has run on all six sets.

- Scripts and design: `scan_gam.py`, `scan_bam.py`, `window_scan.py`, `choose_windows.py` (windows fixed), `compare.py`,
  `vg_stats_crosscheck.py` (`results/vg_stats_crosscheck.md`).
- HG008T PacBio, four arms: G full scans of the GRCh38-only graph / d9 / d46 (387166_0-2); BAM window scan (387287_4).
- HG008T ONT, four arms: G full scans (387166_4, _5, _6); BAM window scan (387287_5).
- HG008T Illumina, four arms: window scans of the GRCh38-only graph, d9, d46 and the BAM (387287_0, _1, _2, _3).
- COLO829T fiber-seq and ONT, four arms each: G full scans (387288_12-14, _16-18); BAM window scans (387287_9, _10).
- COLO829T Illumina, three graph arms: window scans (387287_6-8).
- `compare.py`: `results/{HG008T,COLO829T}_*.json` (six) and `results/tables.md`, 2026-10-15.

Next steps outside this folder:

1. Fill the open `[待补]` placeholders of R2.7 / R1.3 / R3.6 in `../../response_to_reviewers_detailed.md` from the
   numbers above (a user decision on the wording).
2. The ONT alignment of E11b would show whether the more frequent confident disagreement of ONT reads on d46
   (Interpretation 6) is misplacement.

Not planned:

- Full BAM scans (387166_3 cancelled; the BAM arm uses the window design).
- Full Illumina GAM scans (tasks 8-10 and 20-22 never submitted).
- A GRCh38 BAM arm for COLO829T Illumina (none listed).
- CHM13: dropped by the user on 2026-10-13. It was never part of E11a; it belonged to E11b.
- GIAB v3.5 strata (E10).
- The comparators are not re-run; the user has the ClairS-TO / DeepSomatic results.
