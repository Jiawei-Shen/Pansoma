# ver4 item `offref`: Response Tables 3.10 / 3.11, graph-coordinate share, Response Figure 3.3 (2026-11-23)

Graph-coordinate (GRCh38-unprojectable, "unplaced") SNV PASS calls of the nine six-tumor current call sets:
`/scratch/jshen/data/pansoma_net_v2_runs/snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe/vcf_chr1_<SET>_d46`
(seed 0; HPRC v2.1; chr1; default output, no off-reference rescue), listed in `runs.txt`
(`NAME VCFDIR TUMOR PLATFORM`). 935 calls: HG008T 77 / 143 / 90, COLO829T 110 / 79 / 111, HCC1395T 76 / 106 / 143
(Illumina / PacBio HiFi or Fiber-seq / ONT). No model was trained or re-predicted; no existing file was changed.

## Results (this directory)

| file | content | script |
|---|---|---|
| `table_3.10.md`, `table_3.10.json` | Response Table 3.10 (letter layout), totals, per-tumor sums, the keys PASS on all three platforms, controls | `table_3_10.py` |
| `table_3.11.md`, `table_3.11.json` | Response Table 3.11 (letter layout), totals, CHM13 VCF check | `table_3_11.py` |
| `share.md`, `share.json` | graph-coordinate share of all PASS calls per call set, two denominators | `share.py` |
| `fig_3.3_sbs96_graph_coordinate.png` (3000 x 1800 px), `fig_3.3_sbs96.tsv`, `fig_3.3_cosine.tsv` | Response Figure 3.3 and its numbers | `fig_3_3_sbs96.py` (sambcfenv python) |
| `checks.md` | truth counts, call counts per source, E9i against the ver3 runs, fetch, call-set headers, E9c / E9d reproduction | `checks.py` |
| `EXP/` | `collect.py` output of this batch (`per_call.tsv.gz` = assembly class per call; `tables.md`, `summary.json`) | `code/collect.py` |
| `repro/EXP/` | the same for the reproduction run (one ver3 call set) | |

Every builder has `--repro` (except `checks.py`): the same code on the ver3 inputs, compared with the ver3 numbers.

## Steps

| step | code | job script (`tmp/ver4_20261123/offref/`) | output |
|---|---|---|---|
| E9c contexts / GRCh38 locus / controls per call set | `code/e9c_run.py`, `code/e9c_lib.py` (copies of `tmp/r33_hcc1395_20261104/code`; only `DATA` / `EXP`, with an `E9C_DATA` / `E9C_EXP` override for the reproduction) | `e9c.sh` | `<X4>/DATA/runs/m6np__<SET>__vcf_chr1/{calls.tsv,run.json}` |
| E9d assembly classes | `code/e9d_asm.py`, `code/collect.py` (unchanged copies) | `asm.sh HG008T / COLO829T / HCC1395T`, `collect.sh` | `<X4>/DATA/asm/<tumor>/`, `EXP/per_call.tsv.gz` |
| E9i matched-normal / tumor reads | `e9i/normal_reads.py` (one copy for the three tumors; counting code unchanged), `e9i/regions.py`, `e9i/fetch_normal_hg008n.py` | `e9i_fetch_hg008n.sh` (one remote reader), `e9i_count.sh RUN NSHARDS` | `<X4>/e9i/HG008T/normal_subset.bam`, `<X4>/e9i/out/<run>[.<shard>].tsv` |
| E9a / E9b projection + CHM13 VCF | `experiments/E9a_unplaced_projection/project_unplaced.py` and `pansoma_net_v2.linear_vcf` (both unmodified, run in place) | `r34_task.sh` (array over `runs.txt`; the complete-output rtg evaluation of the old task is left out) | `<X4>/r34/{e9a,chm13_vcf}/<run>/` |

`<X4>` = `/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/offref` (0.64 GB on disk). `DATA/edges_d46`,
`DATA/gold_d46_mapping`, `DATA/ref` and the two chr1 cache files are symlinks into
`reviewer_experiments_20261013/E9c_offref_filter_validation` (read only). Job ids and notes: `tmp/ver4_20261123/offref/jobs.txt`;
logs in `tmp/ver4_20261123/offref/logs/`.

Inputs per tumor (unchanged from ver3):

| tumor | normal assemblies | tumor assembly | normal reads | tumor reads |
|---|---|---|---|---|
| HG008T | HG008-N v6.2 (both haplotypes) | HG008-T v3.2 | HG008-N-D PacBio HiFi Revio 68x (GIAB, remote; 822 regions fetched) | HG008-T PacBio HiFi Revio 116x |
| COLO829T | COLO829BL donor-specific assembly hapX / hapY | hifiasm of COLO829T Fiber-seq HiFi + ONT reads, hap1 / hap2 | COLO829BL ONT (SMaHT D001 + D003) | COLO829T ONT |
| HCC1395T | HCC1395-BL hifiasm hap1 / hap2, HCC1395BL v1.0 | hifiasm of HCC1395 HiFi reads, hap1 / hap2 | HCC1395-BL PacBio HiFi (CASTLE; local full BAM) | HCC1395 PacBio HiFi |

## Definitions

See the docstrings of `table_3_10.py`, `table_3_11.py`, `share.py`, `fig_3_3_sbs96.py`. They are the definitions of the
ver3 letter (`tmp/e9i_hg008t_20261104/r33_results.py`, `tmp/r33_hcc1395_20261104/{r33_results.py,r34/r34_results.py}`,
`supplementary/R3.5/make_sbs96.py` part b), which the `--repro` runs confirm on the ver3 inputs.
