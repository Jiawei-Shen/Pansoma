# E19b: the reads and the false calls at germline-dense sites, GRCh38-only against the HPRC graphs

Reviewer point: **R2.5** (does graph-aware calling help when several germline variants are nearby?). Follows E19, which
showed that the HPRC graphs' SNV F1 gain is largest where the normal's germline variants cluster, through precision.
E19b shows why, on the same strata.

**Status (2026-10-16, job 418014 on meow, 51 min): done** for HG008T Illumina / PacBio / ONT and COLO829T Illumina /
fiberseq / ONT, arms linear (GRCh38-only graph), d9 (HPRC v1.1) and d46 (HPRC v2.1).

## Data and method

`dense_reads.py` (this directory; tables in `tables.md`, numbers in `summary.json`). chr1, SNV, inside the benchmark
BED. Germline density n50 = the normal's germline records within ±50 bp, read from E19's per-record files
(`/scratch/jshen/data/reviewer_experiments_20261013/call_level/E19/<call set>.v2.tsv.gz`); strata 0 / 1 / 2+.

- **A. Reads at the truth somatic SNVs** (no model involved). For every truth SNV (E19 `tp-baseline` + `fn` of the
  linear call set) with an exact somatic tensor (`representative_allele_is_somatic_truth` in `chr1_labels.ndjson`) in
  **every** arm of the dataset, each covered read row of the tensor gives:
  - edits = mismatch columns (channel 4 op X) + runs of I / D / complex ops, in the 101-column window, site column 50
    excluded;
  - MAPQ = channel 3.

  Per tensor: mean edits per read, share of reads with ≥ 2 edits, mean MAPQ, share with MAPQ < 20; then the mean over
  the paired sites of the stratum. Tensors are read with one read at the record's byte offset (`np.load` with
  `mmap_mode` faults page by page over BeeGFS; checked identical to `np.load`).
- **B. The GRCh38-only false calls** (PASS at the run's threshold, after the PoN, in the BED; E19 `fp` rows, E4b class
  of each FP from `DATA/E4b/<call set>.eval_bed.fp.v2.tsv.gz`) and the graph arm at the same allele (GRCh38 POS / ALT):
  also a false call there / a tensor that is not called / another allele's tensor at that position / no tensor.

## Results (`tables.md`)

### A. Read context (paired sites)

| dataset | n50 = 0: edits per read lin / d9 / d46 | n50 = 2+: edits per read lin / d9 / d46 | 2+: reads with ≥ 2 edits lin / d9 / d46 | 2+ sites paired / truth |
|---|---|---|---|---|
| HG008T Illumina | 0.29 / 0.29 / 0.29 | 1.94 / 0.80 / 0.72 | 0.54 / 0.16 / 0.13 | 15 / 39 |
| HG008T PacBio | 0.20 / 0.19 / 0.19 | 1.78 / 0.44 / 0.37 | 0.64 / 0.08 / 0.05 | 15 / 39 |
| HG008T ONT | 0.76 / 0.76 / 0.79 | 2.16 / 0.94 / 0.98 | 0.67 / 0.17 / 0.18 | 16 / 39 |
| COLO829T Illumina | 0.21 / 0.21 / 0.21 | 1.81 / 0.55 / 0.57 | 0.53 / 0.10 / 0.10 | 51 / 70 |
| COLO829T fiberseq | 0.20 / 0.20 / 0.20 | 1.95 / 0.37 / 0.38 | 0.67 / 0.05 / 0.05 | 52 / 70 |
| COLO829T ONT | 1.58 / 1.60 / 1.72 | 2.98 / 1.61 / 1.75 | 0.74 / 0.29 / 0.30 | 51 / 70 |

- With no germline variant within 50 bp the three references give the same reads (difference ≤ 0.14 edits per read).
- With two or more, the GRCh38-only reads carry 1.78–2.98 edits per read (53–74 % of reads ≥ 2), the HPRC graph reads
  0.37–1.75 (5–30 %): the germline alleles are graph paths, so the reads follow them without a mismatch. What remains
  on the graphs is about the platform's background (the n50 = 0 row), most clearly on ONT.
- Mean MAPQ is 58.3–60.0 in every arm and stratum, and reads with MAPQ < 20 are ≤ 2 %: the difference is the read
  context of the candidate, not mapping quality.

### B. GRCh38-only false calls, n50 = 2+

| dataset | linear FP | germline-related (same allele / overlapping / ≤ 10 bp) | gone on d9 / d46 | no tensor on d9 / d46 |
|---|---|---|---|---|
| HG008T Illumina | 236 | 193 (81.8 %) | 92.8 % / 91.5 % | 48.7 % / 47.5 % |
| HG008T PacBio | 180 | 175 (97.2 %) | 94.4 % / 93.3 % | 41.1 % / 86.7 % |
| HG008T ONT | 125 | 121 (96.8 %) | 91.2 % / 85.6 % | 68.8 % / 73.6 % |
| COLO829T Illumina | 23 | 5 (21.7 %) | - / 95.7 % | - / 21.7 % |
| COLO829T fiberseq | 33 | 15 (45.5 %) | - / 97.0 % | - / 78.8 % |
| COLO829T ONT | 49 | 28 (57.1 %) | - / 93.9 % | - / 46.9 % |

- The ones still called on the graph are mostly the normal's own germline alleles (e.g. HG008T Illumina d46: 17 of 20).
- For contrast, at n50 = 0 the graph arms still make 58.7–75.5 % of the HG008T linear false calls (isolated rare
  germline alleles, which neither reference contains; E4b), so the graph's removal of linear false calls is specific
  to germline-dense sites.

## Caveats

- A is on the sites with an exact tensor in every arm: 15–16 of the 39 HG008T 2+ truth SNVs (the HPRC graphs have
  exact tensors for fewer of them, E19 `tensor share`), 51–52 of 70 for COLO829T. Small numbers for HG008T.
- The d9 checkpoints of HG008T PacBio / ONT are of another recipe (E3 caveat 1); part A does not depend on models.
- "Gone" = not a false call of the graph arm at the same GRCh38 allele; "no tensor" = no tensor with a GRCh38
  position at that base (an off-reference tensor nearby is not looked for).
