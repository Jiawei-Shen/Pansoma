# E8c: linear tumor-only callers at the alleles the HPRC v1.1 graph absorbs (R1.4)

**Question.** The HPRC graph absorbs some somatic alleles: their ALT reads align without edit to an existing path, so
Pansoma gets no candidate (R1.4, `analysis/graph_absorbed_somatic_20261001`, `..._colo829t_20261005`). Do linear
tumor-only callers keep these alleles, or do they lose them too?

**Inputs (read only).** Absorbed alleles: the "perfect bypass" truths of those two analyses (`per_variant.tsv`,
`<platform>_perfect`), chr1, per platform. Comparators: the user's existing tumor-only runs, re-scored by E8
(`/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score/<caller>/<set>/<KIND>/eval_bed/`; DeepSomatic
1.9.0, ClairS-TO 0.4.2 on HG008T and 0.4.1 on COLO829T, chr1 only, rtg vcfeval in the benchmark BED) at four stages; this
experiment uses `all` (any FILTER), `pass` and `pass_pon` (PASS after Pansoma's PoN rule).

**Method.** `absorbed_comparators.py`: an allele is kept at a stage when it is in that stage's tp-baseline; it counts
as in the BED when it is in tp-baseline or fn of the `all` stage. Writes `tables.md`, `summary.json`.

**Result (2026-10-20, `tables.md`).** Per dataset 4-14 absorbed SNVs and 19-136 absorbed INDELs on chr1.
- Both callers emit a record for many absorbed INDELs (DeepSomatic-TO 28-85 %, ClairS-TO 4-74 %), but their own
  filters pass 7 of the 882 absorbed allele-caller pairs (0 of 106 SNV pairs, 7 of 776 INDEL pairs: DeepSomatic-TO
  HG008T Illumina 6 / 136, COLO829T Illumina 1 / 34), against 56-91 % of all chr1 truth SNVs and 2-21 % of all truth
  INDELs at PASS.
- After the same population filtering, none of the absorbed alleles is called by either caller.
- So the linear callers lose these alleles too, at their germline filtering instead of at alignment.
