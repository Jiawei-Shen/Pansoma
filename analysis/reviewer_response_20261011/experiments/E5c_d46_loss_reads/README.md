# E5c: why truth alleles lose their tensor on HPRC v2.1 d46 (read level)

Reviewer points: **R1.2 / R1.4** (does the graph absorb somatic alleles?). User request 2026-10-16: count a truth allele
as having a tensor when any tensor allele equals it (representative or not), and show whether the reads at the lost
loci align entirely or partly to graph alt nodes, so that no tensor is generated.

**Status (2026-10-19): done, all six d46 sets.** The ONT sets ran out of memory at 24 GB on 2026-10-17 and were rerun
with 64 GB (423550, 423551, 423552); every set is complete (`tables.md`, `summary.json`, `DATA/per_locus.tsv.gz`).

## Results (`tables.md`; chr1-22, all truth records)

Lost INDELs = a tensor allele equals the truth on GRCh38-only and none does on d46. Class = how the majority of the
locus's error-free ALT reads align on d46.

| set | lost INDELs | I1: ALT reads walk graph alt nodes, no edit | I2: alt nodes plus a residual edit | I1 + I2 | I3: GRCh38 nodes, another edit | rest |
|---|---|---|---|---|---|---|
| HG008T Illumina | 3,528 | 47.5 % | 48.0 % | 95.5 % | 2.2 % | 2.3 % |
| HG008T PacBio | 3,577 | 43.2 % | 52.4 % | 95.6 % | 2.8 % | 1.6 % |
| HG008T ONT | 2,494 | 44.5 % | 51.4 % | 95.8 % | 2.5 % | 1.7 % |
| COLO829T Illumina | 612 | 59.2 % | 30.2 % | 89.4 % | 6.5 % | 4.1 % |
| COLO829T Fiber-seq | 567 | 53.3 % | 32.5 % | 85.7 % | 9.9 % | 4.4 % |
| COLO829T ONT | 457 | 52.1 % | 33.3 % | 85.3 % | 9.6 % | 5.1 % |

- **The loss is the graph.** 85-96 % of the lost loci are I1 or I2 (86.9-97.1 % of those with >= 3 ALT reads have
  most ALT reads on graph alt nodes; pooled ALT reads: 84-96 % entire or partial). The same method on 150 control loci per set (a tensor on both arms) finds the
  truth edit made at 77-98 % of the control INDELs (ONT lowest) and 98-100 % of the control SNVs.
- **Entire and partial:** about equal on HG008T (43-48 % / 48-52 %); entire dominates on COLO829T (52-59 % / 30-33 %).
  Partial = the read follows a graph allele of the same repeat (e.g. another homopolymer length) and makes a residual
  edit; that edit becomes a candidate with another key or falls below the AF / min-variant cut, so the truth key gets
  no tensor.
- **By repeat context:** HP >= 7 and STR make up 81-88 % of the lost INDELs on HG008T and 61-75 % on COLO829T; these
  losses are 94.0-99.6 % I1 + I2; non-repeat INDEL losses 65.2-80.3 %, with
  13.6-25.2 % I3 (reads stay on GRCh38 and spell the event differently).
- **After population filtering** (lost and not PoN-tagged under the current rule: HG008T 515 / 603 / 355, COLO829T
  118 / 133 / 110 for Illumina / PacBio or Fiber-seq / ONT): mostly partial alignments, I2 80-82 % on HG008T and
  50-57 % on COLO829T; entire absorption only 6-9 %; I3 8 % (HG008T) and 28-35 % (COLO829T). HP >= 7 or STR: 78-85 %
  of these on HG008T, 40-55 % on COLO829T.
- **Size of the loss** (same criterion, `tmp` recount from the recall tables and the E5b PoN tags): 24-42 % of truth
  INDELs lose their candidate on d46 (HG008T 41.5 / 42.1 / 29.4 %, COLO829T 32.0 / 29.7 / 23.9 %); 76-86 % of these are
  PoN-tagged anyway; net of d46-only tensors, the loss after the PoN is 3.5-7.0 % of truth INDELs.
- **SNVs** lose far less (HG008T 263-316 lost, 63-67 % I1 + I2). On COLO829T only 15-25 % of the 220-398 lost SNVs are
  graph-related; most have too few ALT reads (I4, 32-33 %), no read with MAPQ > 10 (I5, Illumina and Fiber-seq
  32-46 %) or, on ONT, the truth edit is made but the candidate is filtered (I10, 34.5 %).

- Loci per d46 set: every truth SNV / INDEL (chr1-22, all records) with a tensor on the GRCh38-only set and none on d46
  (status no_candidate or filtered), plus 150 controls with a tensor on both.
- `reads.py`: the read-level classifier of `analysis/hg008_somatic_miss_20260926` (I1-I10), unchanged in method, run
  with each d46 run's frozen v4 source on its GAM. Per locus it counts the ALT reads that walk graph alt nodes / a skip
  edge with no edit (entire), with a residual edit (partial), stay on GRCh38 with another edit, or make the truth edit.
- `report.py`: tables by d46 status, repeat class (E7: HP>=7 / STR / other) and PoN tag (E5b).
- Jobs (`tmp/reviewer_cpu_20261013/E5c/`, `jobs.txt`): test 417646 (passed: 6 lost loci all walk alt nodes, the
  3 controls make the truth edit); HG008T PacBio 417634 and COLO829T Fiber-seq 417635 first (first conclusion), then
  COLO829T ONT / Illumina, HG008T Illumina, HG008T ONT (417636-417639). 8 CPUs, 24 GB per part; maotai excluded
  (BeeGFS stall).
- Data: `/scratch/jshen/data/reviewer_experiments_20261013/E5c_d46_loss_reads/<set>/reads_part*.jsonl`.
