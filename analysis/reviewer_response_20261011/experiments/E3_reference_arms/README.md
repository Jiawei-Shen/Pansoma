# E3: reference-level ablation (GRCh38-only vs HPRC v1.1 d9 vs HPRC v2.1 d46)

Reviewer points: **R3.1** (reference ablation), **R1.1** (is the method pangenome-native), **R1.2** (graph size and
composition), **R2.5** (what drives the precision), **R3.2** (what the graph removes).

## Data

The 2026-10-12 RIS runs (`analysis/benchmark_7fc7396_20261012`; pipeline 7fc7396). Each of the 12 sets (HG008T and
COLO829T × Illumina / PacBio (COLO829T: fiberseq) / ONT × SNV / INDEL) is trained and tested once per reference:

| arm | graph | H (ch6) | aligner |
|---|---|---:|---|
| linear | GRCh38 only (`hg38_linear.giraffe.gbz`, vg autoindex) | 17 | vg giraffe 1.66 |
| d9 | HPRC v1.1 Minigraph-Cactus GRCh38, clipped at d9 | 90 | vg giraffe 1.66.0 for HG008T (logs, `points/R3.1.md`); COLO829T [unverified] |
| d46 | HPRC v2.1 Minigraph-Cactus GRCh38, clipped at d46 | 464 | vg giraffe 1.77 |

The same reads (read counts identical across arms), the same tensor build (9d6bd8d, build AF 0.08), the same labels
(`5430533aa100`) and truth, the same recipes, and chr1 held out. `arms_table.py` writes `tables.md` from
`linear_grch38_7fc7396_20261012/summary.json`.

## Results (chr1; `tables.md`)

**SNV: both pangenomes beat GRCh38-only after the PoN on all six sets.** PASS F1 in the BED:

| sample | Illumina | PacBio / fiberseq | ONT |
|---|---|---|---|
| HG008T: linear / d9 / d46 | 0.469 / 0.585 / 0.570 | 0.513 / 0.621 / 0.605 | 0.537 / 0.596 / 0.575 |
| COLO829T: linear / d9 / d46 | 0.468 / 0.778 / 0.783 | 0.680 / 0.781 / 0.779 | 0.418 / 0.706 / 0.737 |

- The gain over linear is +0.04 to +0.32 F1. It holds at the test-tuned best point too: linear is below d9 and d46 on
  all six sets.
- d9 and d46 are equivalent after the PoN: −0.021 to +0.031.
- Mechanism: on the linear graph every germline variant of the sample is a candidate.
  - HG008T Illumina chr1 has 414,453 SNV tensors against 190,697 (d9).
  - Per-tensor somatic AP is 0.012–0.016 against 0.20–0.35 (HG008T).
  - Raw PASS F1 is 0.008–0.075 against 0.18–0.69.
  - The PoN removes most of these germline calls, but 1.4–7.8× more false calls survive it than on d9.

**INDEL: no pangenome advantage after the PoN, except on Illumina.**

- Illumina: d9 / d46 0.24 against linear 0.13 (HG008T), and 0.11 against 0.03 (COLO829T).
- Long reads: linear PASS F1 is higher (HG008T PacBio 0.125 vs 0.089 / 0.095), but its test-tuned best F1 is lower on
  every long-read set (0.173 vs 0.260 / 0.269).
- So the long-read INDEL PASS differences come from where each model's validation recall-0.9 threshold falls, not from
  the ranking.
- The linear arm keeps the INDELs that the pangenome absorbs: truth with a tensor is 0.75–0.87 against 0.54–0.76.
  Its PoN-limited ceiling is higher (HG008T 0.284 / 0.301 / 0.238 against 0.257 / 0.259 / 0.238 for d9), but its
  ranking is worse.

**Complete-output bound (`tables.md` §2).** Unplaced PASS calls (off-reference, no GRCh38 position) are not scored by
rtg. Counting every one of them as a false positive (worst case, no BED):

- SNV d9 / d46 stay above linear on Illumina (both samples) and on HG008T PacBio.
- They fall below linear in four cases:
  - COLO829T fiberseq d9: 0.636 against 0.680; d46 stays above at 0.736;
  - HG008T ONT: 0.342 (d9) / 0.354 (d46) against linear 0.365;
  - COLO829T ONT d46: 0.369 against linear 0.418, with 4,566 unplaced PASS calls against 1,690 true positives.
    d9 stays above at 0.522.
- Every unplaced SNV PASS call is made by predict's off-reference rescue (`offref_rescue.calls` equals the unplaced
  count in the nine runs checked; `points/R3.4.md`).
- Without the rescue, the complete-output F1 of the three d46 runs re-predicted on tequila would be at least 0.596 /
  0.524 / 0.721 (HG008T PacBio / HG008T ONT / COLO829T ONT), above linear.
- E9a (projection) and E9c / E9d (graph-native population filter, assembly validation) replace this bound with measured
  numbers.

## Caveats (stated in the response)

1. **d9 comparability.** The d9 runs match the linear and d46 runs exactly only for HG008T Illumina SNV and the HG008T
   INDELs.
   - HG008T PacBio / ONT SNV d9 trained on 1 GPU (global batch 1,024); linear and d46 trained on 4 GPUs (global batch
     4,096). The batch change alone cost HG008T Illumina d9 0.007.
   - COLO829T d9 SNV are the 2026-10-03 one-GPU runs (`*_SNV_h90`).
   - COLO829T d9 INDEL were trained at lr 4e-4.
   - The matched d9 reruns are GPU jobs, in the RIS list.
2. **One seed per arm.** The variance needs the seed repeats (R3.9; GPU jobs) and a bootstrap over chr1 blocks (CPU,
   once the call sets are local).
3. **The linear arm is Pansoma on a linear graph, trained with the graph recipe.**
   - Its validation somatic AP is about 0.01 on HG008T, and its PASS thresholds are 0.04–0.09.
   - Whether a dedicated linear caller does better is answered by E8 (ClairS-TO, DeepSomatic on GRCh38 BAMs, same
     protocol).
   - The claim is therefore "graph candidates + ch6 vs the same network on GRCh38", not "vs any linear method".
4. **vg versions differ between arms**: 1.66 for linear and the HG008T d9 GAMs, 1.77 for d46, and the COLO829T d9
   GAMs unverified. So HG008T linear vs d9 has no aligner-version effect. A version control (the linear graph with vg 1.77)
   needs a re-index and re-alignment. It is listed, not run.
5. **No-ch6 arm.** The retraining without the path-count plane (`--drop-planes path_count`) on d9 / d46 is a GPU job
   (R1.1). E1.1 (CPU counterfactual on the trained models) measures how much the trained models use ch6.

6. **CHM13 is in both HPRC graphs as an unfiltered reference path** (`vg paths -M`: REFERENCE sense, 109 subpaths in
   d9, 93 in d46; checked 2026-10-15). So "GRCh38-only vs HPRC" mixes two things: a more complete reference (CHM13
   sequence) and population haplotype diversity. A CHM13+GRCh38-only graph would separate them; CHM13 work is dropped
   for now (user decision 2026-10-13), so this stays a stated limitation.

## Status

- Done: the three-arm table, mechanism, and the complete-output worst-case bound.
- Pending: matched d9 reruns and seeds (GPU, RIS); E9 measured complete-output metrics; bootstrap CIs; E8 comparators.
