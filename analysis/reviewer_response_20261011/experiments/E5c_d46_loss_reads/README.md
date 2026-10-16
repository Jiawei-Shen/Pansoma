# E5c: why truth alleles lose their tensor on HPRC v2.1 d46 (read level)

Reviewer points: **R1.2 / R1.4** (does the graph absorb somatic alleles?). User request 2026-10-16: count a truth allele
as having a tensor when any tensor allele equals it (representative or not), and show whether the reads at the lost
loci align entirely or partly to graph alt nodes, so that no tensor is generated.

**Status (2026-10-16 18:10): running.** Results, tables and the conclusion will be added here.

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
