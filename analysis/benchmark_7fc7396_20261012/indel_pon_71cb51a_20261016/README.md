# The RIS INDEL call sets re-scored with the 71cb51a INDEL PoN rule (2026-10-16)

71cb51a changed the INDEL rule of `scripts/filter_panel_of_normals.py`: a gnomAD / CoLoRSdb exact allele tags an INDEL
at AF >= 0.01 when its ALT changes a homopolymer of >= 8 bp, else at AF >= 1e-4 (7fc7396: >= 7 bp, else 1e-3; chosen
on chr1 calls with one tumor left out, `analysis/reviewer_response_20261011/experiments/E17c_pon_rule_offchr1`). The
SNV rule is unchanged. Every INDEL call set on RIS (59: the 7fc7396 benchmark runs, the purity joint runs per purity,
the 2026-10-03 / 10-06 d9 runs, and the reviewer GPU runs of `../reviewer_gpu_7fc7396_20261016`) was re-scored with it.

- `rescore.sh`: `ris/vcfeval_only.sh` (PoN + rtg vcfeval only) on each run's existing GRCh38 calls
  `<run>/<SRC>/*.INDEL.linear.vcf.gz` into `<run>/<SRC>_pon71cb51a/eval_{bed,nobed}`; the calls and their PASS
  threshold are unchanged, the earlier evaluations are kept. LSF jobs 141864 (a test: HG008_Illumina_INDEL_7fc7396,
  PASS F1 0.241 -> 0.251, as E17c's 0.251) and 141869-141932.
- `collect.py` (standard library): `tables.md` and `summary.json`, the 7fc7396-rule and the 71cb51a-rule scores of each
  call set side by side, in the BED and without it. The 7fc7396-rule scores of the 2026-10-06 runs (`*_INDEL_small`)
  are their `vcf_chr1_pon7fc7396`; the 2026-10-03 `*_INDEL_h90` runs had no INDEL PoN before.

## Results (chr1, in the BED, after the PoN; `tables.md`)

- Small and in one direction: PASS F1 +0.000 to +0.013, best F1 +0.000 to +0.010 on every call set. The new rule tags
  a few more calls (HG008T Illumina d9: 6,604 -> 6,736) and 5 more HG008T truth INDELs (413 -> 418 of the chr1 truth;
  COLO829T unchanged at 81), so the PoN-limited ceilings are unchanged to three decimals.
- No comparison changes direction (best F1, 7fc7396 rule -> 71cb51a rule):
  - d46 - d9: HG008T -0.015 / +0.009 / +0.015 -> -0.013 / +0.012 / +0.016 (Illumina / PacBio / ONT); COLO829T against
    the 2026-10-06 d9 runs -0.033 / +0.009 / +0.013 -> -0.034 / +0.013 / +0.015; against the matched lr 1e-4 d9 runs
    (`COLO829_*_INDEL_7fc7396`) +0.007 / -0.005 / -0.006.
  - linear below d9 and d46 on all six sets (HG008T 0.201 / 0.176 / 0.084, COLO829T 0.070 / 0.095 / 0.020).
  - no ch6 against with ch6 (`../reviewer_gpu_7fc7396_20261016`): HG008T d9 -0.010 / -0.034 / -0.010, d46 -0.004 /
    -0.030 / -0.001; COLO829T Illumina +0.011 (d9) / +0.010 (d46), fiberseq -0.029 / -0.051 and ONT -0.013 / +0.013,
    on fewer than 30 true calls at the best points.
  - purity joint: Illumina 0.167 / 0.312 / 0.330 / 0.233 (P100 / P75 / P50 / P25), fiberseq 0.069 / 0.049 / 0.055 /
    0.074, ONT 0.055 / 0.051 / 0.070 / 0.062; Illumina with BQ per read 0.182 / 0.146 / 0.133 / 0.133.
- The sub-folders' `tables.md` and READMEs keep the 7fc7396-rule INDEL numbers they were written with; the benchmark
  README's main INDEL results now quote the 71cb51a rule.
- The two COLO829T ONT INDEL no-ch6 runs (d9, d46) were re-scored after their own `vcf.sh` (jobs 141931-141932);
  their `vcf.sh` had already run with the 71cb51a rule (the * rows of `tables.md`, same scores).
