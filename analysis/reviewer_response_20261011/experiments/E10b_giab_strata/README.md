# E10b: calls per GIAB v3.5 stratification (chr1), HPRC graphs against the GRCh38-only (linear) graph

Reviewer point: **R1.3** (performance in difficult regions; a region-specific benchmark). Follows E10 (SMaHT easy /
difficult / extreme) and E16's per-stratum bootstrap (`../E16_bootstrap/strata.md`), which used the SMaHT classes only
because the GIAB stratifications were not on the server.

**Status (2026-10-20): done.** `giab_strata.py` (job 449747, 1 min 53 s; `tmp/reviewer_cpu_20261020/E10b_giab/`) wrote
`strata.md`, `strata_ci.tsv`, `strata_paired.tsv`, `summary.json`.

## Data

- **Strata**: GIAB genome stratifications v3.5, GRCh38 (`genome-stratifications-GRCh38@all.tar.gz`, downloaded
  2026-10-20 with the user's approval from ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/,
  md5 b70ac852472cac7747676b6b6b0c3465 checked against `genome-stratifications-md5s.txt`; the 12 files used are extracted to
  `/scratch/jshen/data/GIAB_stratifications_v3.5/GRCh38@all/`, list in `members.txt` there). Chosen for R1.3 and present
  on chr1: the union of all difficult regions and its complement, low mappability or segdup, low mappability, segmental
  duplications, tandem repeats and homopolymers (all; homopolymers >= 7 bp; tandem repeats <= 50 bp), satellites,
  L1H > 500 bp, GC < 25 % or > 65 %, CMRG falsely duplicated genes (no chr1 records). MHC / KIR / VDJ are not on chr1.
- **Records**: E10's per-record tables (`/scratch/jshen/data/reviewer_experiments_20261013/call_level/E10/`), eval_bed rows
  (PASS after the PoN, in the benchmark BED, rtg parts tp-baseline / fn / tp / fp), in-sample call sets at the main
  threshold, SNV and INDEL, linear / d9 / d46 (as `e16_strata.py`). A record is in a stratum when its VCF position lies in
  one of the stratum's chr1 intervals; strata overlap, so each is scored on its own.
- **Bootstrap**: E16's 1-Mb blocks of the set's BED, 2,000 resamples, seed 20261015; HPRC - linear differences on the
  same resampled blocks.

## Results (`strata.md` section 1; 12 pairs = HG008T / COLO829T x three platforms x v1.1 d9 / v2.1 d46)

| stratum | chr1 Mb | SNV dF1 range | SNV: > 0 / CI above 0 / CI below 0 | INDEL dF1 range | INDEL: > 0 / above / below |
|---|---:|---|---|---|---|
| not in any difficult region | 183.7 | -0.001 to +0.176 | 11 / 7 / 0 | -0.001 to +0.280 | 11 / 9 / 0 |
| all difficult regions | 46.8 | +0.054 to +0.483 | 12 / 10 / 0 | -0.053 to +0.090 | 4 / 4 / 7 |
| low mappability or segdup | 23.2 | +0.033 to +0.528 | 12 / 10 / 0 | -0.078 to +0.125 | 4 / 2 / 3 |
| low mappability | 19.8 | +0.034 to +0.548 | 12 / 10 / 0 | -0.080 to +0.123 | 4 / 2 / 3 |
| segmental duplications | 14.8 | -0.063 to +0.331 | 11 / 9 / 0 | -0.077 to +0.056 | 2 / 1 / 3 |
| tandem repeats and homopolymers | 17.4 | +0.031 to +0.513 | 12 / 11 / 0 | -0.053 to +0.089 | 4 / 4 / 7 |
| homopolymers >= 7 bp | 6.7 | +0.025 to +0.160 | 12 / 5 / 0 | -0.065 to +0.083 | 4 / 4 / 7 |
| tandem repeats <= 50 bp | 2.3 | -0.080 to +0.342 | 9 / 4 / 0 | -0.033 to +0.130 | 5 / 2 / 0 |
| satellites (8 SNV pairs) | 3.4 | +0.002 to +0.490 | 8 / 5 / 0 | - | - |
| GC < 25 % or > 65 % | 14.2 | +0.008 to +0.232 | 12 / 8 / 0 | -0.040 to +0.088 | 5 / 2 / 7 |

- **SNVs**: the HPRC graphs are better than the linear graph in every difficult stratum and never significantly worse.
  In the union of difficult regions all 12 pairs gain (interval above zero in 10); low mappability / segdups / tandem
  repeats 9-11 of 12. Truth SNVs per set in the union: 180 (HG008T) and 443 (COLO829T).
- **INDELs**: outside the difficult regions the graphs gain (9 of 12 above zero). In homopolymers and tandem repeats,
  where most truth INDELs lie (HG008T 465 of 552 difficult-region truth INDELs are in homopolymers >= 7 bp), Illumina
  gains on both tumors (HG008T +0.080 / +0.083, COLO829T +0.041 / +0.056, all above zero) and the long reads lose: HG008T
  PacBio -0.065 / -0.056 and COLO829T fiber-seq -0.029 / -0.032 through precision (PacBio PASS INDELs in homopolymers
  1,180 on linear against 2,101 / 1,864 on d9 / d46 at the same 109-123 true positives), ONT about -0.001 to -0.007 at
  F1 near zero on every arm. Every reference scores at most 0.238 INDEL F1 in these regions.

## COLO829T Illumina v1.1 at epoch 6 (2026-10-20)

The user decided on 2026-10-20 (RIS b662f00) that the COLO829T Illumina HPRC v1.1 (d9) arm uses epochs 6 / 5 / 6 instead
of the validation-selected epoch 8, which calls unlabelled tensors in two chr1 segments (R3.1, R3.9). E10 / E16 /
`giab_strata.py` above still hold the epoch-8 run (`COLO829_Illumina_SNV_7fc7396_ddp4_b1024`).
`colo_ill_d9_ep6.py` (job 449749) rescores that pair with seed 0 at epoch 6
(`ris_7fc7396_20261012/COLO829_Illumina_SNV_d9_normbq_sweep_e06`, `--norm dataset`; PASS F1 0.773) against the same
linear and d46 runs, per SMaHT and GIAB stratum (`colo_ill_d9_ep6.md`, `.json`); the epoch-8 rows reproduce E16's values.
`summary_ep6.py` (`summary_ep6.md`) is the SNV summary with that row replaced:

| stratum | dF1 range | > 0 / CI above 0 / CI below 0 (of 12) |
|---|---|---|
| SMaHT extreme | +0.098 to +0.531 | 12 / 11 / 0 |
| SMaHT difficult | -0.005 to +0.259 | 11 / 3 / 0 |
| GIAB all difficult regions | +0.091 to +0.483 | 12 / 11 / 0 |
| low mappability | +0.106 to +0.548 | 12 / 11 / 0 |
| segmental duplications | +0.018 to +0.331 | 12 / 10 / 0 |
| tandem repeats and homopolymers | +0.091 to +0.513 | 12 / 12 / 0 |

COLO829T Illumina v1.1 itself: SMaHT extreme 243 PASS calls at precision 0.593 (epoch 8: 1,725 at 0.076), dF1 +0.503
[+0.213, +0.588]; GIAB all difficult +0.473 [+0.154, +0.559]. The INDEL rows are unchanged (the decision concerns the SNV
model).

## Caveats

- chr1 only, one training run per arm (the three-seed spread of SNV PASS F1 is up to about 0.08, RIS summary); the
  intervals cover block sampling of chr1, not training variance.
- The COLO829T d9 runs are the call sets present locally in `reviewer_cpu_20261013` (as E10 / E16); E16 and the RIS
  summary discuss which d9 recipe each is.
- Small strata (L1H, satellites, tandem repeats <= 50 bp) have few truth records; read their rows as descriptive.
