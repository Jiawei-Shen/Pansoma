# E3b: SNV without the off-reference rescue; INDEL at a validation-chosen threshold

> **2026-10-16, INDEL PoN rule changed** (user decision: gnomAD / CoLoRSdb exact allele, AF >= 0.01 in homopolymers >= 8 bp, else AF >= 1e-4; was >= 7 bp and 1e-3). The 44 local chr1 INDEL call sets were re-scored in place (`tmp/reviewer_cpu_20261013/pon_hp8_20261016/rescore.sh`, array 418042; the old evaluations are kept as `eval_bed.pon_hp7_1e-3` / `eval_nobed.pon_hp7_1e-3`) and this experiment's collector was rerun (job 418071), so `tables.md` / `summary.json` / TSVs now use the new rule. INDEL numbers quoted in the prose below may still be those of the earlier rule (changes are small: e.g. HG008T Illumina v1.1 PASS F1 0.241 -> 0.251, best 0.265 -> 0.270; HG008T chr1 truth INDELs tagged 413 -> 418 of 635). SNV results are unchanged.


Reviewer points:
- **R3.4 / R3.1 / R1.3 / R2.5**: complete-output metrics, i.e. what the graph-coordinate calls do to the comparison.
- **R1.6 / R3.10 / R3.m1**: an INDEL operating point chosen without the test chromosome.

No new inference. Each existing chr1 prediction set is re-called through the pipeline's own
`machine_learning/pansoma_net_v2/jobs/vcf.sh` (graph_vcf → linear_vcf → PoN → rtg vcfeval, in the BED and without it).
- The prediction sets are the tequila GPU re-predictions of the RIS checkpoints under
  `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/`.
- Script: `variant_calls.sh` (modes `norescue`, `valbest`).
- Driver: `tmp/reviewer_cpu_20261013/E3b/submit.sh`, with `loop.sh` running it as runs finished; small CPU jobs on tequila, 4 CPUs, 2 at a time.
- Tables: `collect.py` → `tables.md`, `summary.json`.

## Method

- **norescue (SNV, graph runs).**
  - predict's off-reference rescue (`--offref-site90 10`) adds `offref_call` / `p_offref_*` to off-reference tensors, and graph_vcf calls those tensors by `offref_call` (argmax after raising the site's path count).
  - The fields are removed from a copy of the predictions (`<run>__norescue/test_chr1/`), so graph_vcf calls every tensor by the checkpoint's validation recall-0.9 threshold.
  - Everything else is unchanged.
- **complete-output F1** = 2TP / (2TP + FP + FN + U).
  - TP / FP / FN are rtg's after the PoN.
  - U = graph-coordinate (unplaced) PASS calls, all counted as false positives. They have no GRCh38 position, so neither the BED nor the PoN applies to them.
- **valbest (INDEL, all arms).** The PASS threshold is the checkpoint's validation best truth-level F1 threshold (`metrics.jsonl` `val.truth.best.threshold` at the best epoch). It is chosen on the 5 % validation node blocks of chr2-22, never on chr1. The default PASS threshold is the validation recall-0.9 threshold.

## Results (`tables.md`; chr1, after the PoN)

**SNV.** Without the rescue the placed calls are identical (rtg TP / FP / FN unchanged); only the graph-coordinate
PASS calls drop. No BED, complete-output F1, rescue on → off, against the GRCh38-only graph (E3 §2, no graph-coordinate
calls):

| set | HPRC v1.1 (d9) | HPRC v2.1 (d46) | GRCh38-only |
|---|---|---|---|
| HG008T Illumina | 0.546 → 0.574 (115 → 12 calls) | 0.414 → 0.553 (754 → 29) | 0.197 |
| HG008T PacBio | 0.553 → 0.612 (210 → 9) | 0.533 → 0.596 (240 → 6) | 0.359 |
| HG008T ONT | 0.342 → 0.579 (1,410 → 26) | 0.355 → 0.524 (1,258 → 162) | 0.365 |
| COLO829T Illumina | (not local) | 0.702 → 0.780 (516 → 18) | 0.468 |
| COLO829T Fiber-seq | (not local) | 0.736 → 0.765 (257 → 78) | 0.680 |
| COLO829T ONT | (not local) | 0.369 → 0.721 (4,569 → 98) | 0.418 |

- With the rescue off, every HPRC arm is above the GRCh38-only graph on the complete output (9 of 9 local
  comparisons).
- With the rescue on, 3 of the 9 are below: HG008T ONT on both graphs and COLO829T ONT on d46. A fourth reversal in
  E3 §2, COLO829T Fiber-seq on d9, uses a run whose predictions are not local.
- The rescue's benefit, the assembly-derived off-reference somatic SNVs it recovers (E1.1 §5), is invisible to rtg.
  Whether to keep it on by default is open; every SNV result can be reported both ways.
- The COLO829T d9 single-set runs (`*_SNV_h90`) were not transferred from RIS, so their rows are missing.

**INDEL.** PASS at the validation-chosen threshold (valbest) against the validation recall-0.9 threshold (r0.9), in
the BED:

| run | r0.9 F1 | valbest F1 | chr1 best F1 (test-tuned) |
|---|---:|---:|---:|
| HG008T Illumina d9 / d46 / linear | 0.241 / 0.244 / 0.133 | 0.216 / 0.204 / 0.080 | 0.265 / 0.249 / 0.200 |
| HG008T PacBio d9 / d46 / linear | 0.089 / 0.095 / 0.125 | 0.226 / 0.199 / 0.036 | 0.259 / 0.269 / 0.173 |
| HG008T ONT d9 / d46 / linear | 0.015 / 0.015 / 0.022 | 0.116 / 0.151 / 0.047 | 0.140 / 0.155 / 0.083 |
| COLO829T Illumina d46 / linear | 0.109 / 0.032 | 0.044 / 0.030 | 0.163 / 0.069 |
| COLO829T Fiber-seq d46 / linear | 0.020 / 0.040 | 0.090 / - (no call) | 0.129 / 0.092 |
| COLO829T ONT d46 / linear | 0.002 / 0.004 | 0.027 / 0.000 | 0.085 / 0.019 |

- The validation-chosen threshold rescues the long-read HPRC INDEL models on HG008T, reaching 74-97 % of the
  test-tuned best F1: PacBio d9 0.089 → 0.226; ONT d46 0.015 → 0.151.
- It costs a little on Illumina (0.241 → 0.216).
- It is too strict for the GRCh38-only models and for COLO829T, where only 0-46 calls pass. COLO829T has 131 truth
  INDELs on chr1 and the validation blocks hold very few somatic INDELs (52-92), so its validation optimum is noisy.
- No single held-out rule is uniformly good, so the INDEL results should be reported at a stated held-out operating
  point, with the test-tuned best F1 labelled as an upper bound.

## Caveats

- One seed and chr1 only; confidence intervals are in E16 (the norescue call sets are identical inside rtg, and only
  the graph-coordinate count changes).
- COLO829T d9 single-set models are missing locally; the COLO829T Fiber-seq joint model (P100 + dilutions) is not
  re-called here.
- The rescue-off counts are exact re-calls; they replace the earlier estimate (6 / 162 / 98 unplaced at p ≥ t for the
  three d46 runs of 2026-10-13, now confirmed).
