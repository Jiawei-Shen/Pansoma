# E17a: INDEL results by type, length and repeat context (R3.10)

> **2026-10-29, final INDEL models** (`e17a_final_models.py` -> `final_models.md` / `.tsv` / `.json`; R3.10, R1.6, R2.3).
> The same strata for the final models of the INDEL pipeline (`tmp/indel_gap_20261020/final_20261023/`: PASS = ppr8,
> release = chr2-22 of six tumors, tumor-unseen = the four cell lines, three seeds each, seeds pooled) and for both
> comparators (PASS + our PoN, `../E16_bootstrap/comparators.tsv`). All 36 Pansoma call sets match report.json. chr1, in
> the BED; release models:
> - Overall F1 (seed means): HG008T 0.199 / 0.249 / 0.168 (Illumina / PacBio / ONT), COLO829T 0.174 / 0.170 / 0.090; ClairS-TO 0.132 /
>   0.129 / 0.127 and 0.122 / 0.158 / 0.103; DeepSomatic-TO 0.208 / 0.053 / 0.034 and 0.191 / 0.133 / 0.021. Recall
>   0.117-0.169 (seed means); in the BED the PoN tags 384 of 587 HG008T truth INDELs (bound 0.346) and 81 of 131 on COLO829T
>   (0.382); recall among the truth INDELs the filter keeps is 0.31-0.49.
> - Type: DEL recall above INS in every Pansoma set on HG008T (0.19-0.21 against 0.10-0.14), as for both comparators.
> - Length: 1 bp is recalled best (HG008T 0.19-0.25); 2-4 bp worst (0.02-0.06 on both tumors).
> - Context: STR truth is almost never called (0.00-0.04; 57 of 65 HG008T STR truths are PoN-tagged); outside repeats recall
>   is 0.31-0.44.
> - Long homopolymers: on HG008T, 44 / 61 / 70 % of Pansoma's FPs are HP>=7 calls (Illumina / PacBio / ONT), precision of
>   HP>=7 calls 0.44 / 0.46 / 0.14. COLO829T ONT FPs spread over HP>=7 (105 per seed), HP4-6 (about 99) and no context (96).
> The sections below describe the 2026-10-15 call sets (earlier models, recall-0.9 / valbest thresholds).

> **2026-10-16, INDEL PoN rule changed** (user decision: gnomAD / CoLoRSdb exact allele, AF >= 0.01 in homopolymers >= 8 bp, else AF >= 1e-4; was >= 7 bp and 1e-3). The 44 local chr1 INDEL call sets were re-scored in place (`tmp/reviewer_cpu_20261013/pon_hp8_20261016/rescore.sh`, array 418042; the old evaluations are kept as `eval_bed.pon_hp7_1e-3` / `eval_nobed.pon_hp7_1e-3`) and this experiment's collector was rerun (job 418071), so `tables.md` / `summary.json` / TSVs now use the new rule. INDEL numbers quoted in the prose below may still be those of the earlier rule (changes are small: e.g. HG008T Illumina v1.1 PASS F1 0.241 -> 0.251, best 0.265 -> 0.270; HG008T chr1 truth INDELs tagged 413 -> 418 of 635). SNV results are unchanged.


Reviewer point: **R3.10** (the INDEL claims: where do Pansoma's INDEL calls work and where do they fail?). Also feeds
R1.6 (sensitivity) and R2.6 (the PoN).

## Question

For every chr1 INDEL call set (PASS after the PoN, inside the benchmark BED): TP / FP / FN, recall and precision by

- **type**: INS / DEL (complex = REF and ALT both keep bases after trimming);
- **length** (|len ALT − len REF|): 1, 2-4, 5-10, >10;
- **repeat context**: HP>=7 (homopolymer tract >= 7 bp), HP4-6, STR (unit 2-6, >= 3 units), VNTR (unit > 6, >= 2 units),
  none.

The PoN-tagged truth per stratum gives the recall bound that no caller can exceed after the PoN.

## Method

- Call sets: every finished `vcf_chr1*/` of `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/` with kind
  INDEL (discovery in `../E16_bootstrap/callsets.py`): in-sample runs of the three reference arms, cross-sample `xs__`
  runs, the COLO829T fiber-seq joint model, and `vcf_chr1_valbest` (the PASS threshold chosen on the chr2-22 validation
  blocks; E3b) next to the standard recall-0.9 threshold. Comparators enter through `../E16_bootstrap/comparators.tsv`
  (slot; empty).
- Records: the rtg output of `jobs/vcf.sh`, `eval_bed/vcfeval/pon_calls/` (tp-baseline and fn = truth side; tp and fp =
  call side). Recall uses the truth record's class, precision the call's class. rtg can match a call and a truth record
  of different representation, so the two sides are classified separately.
- Classes come from the label pipeline's own code, imported and not modified
  (`indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py`):
  - `split_alleles` trims the shared prefix and suffix;
  - `RepeatContext.classify` gives the context. It uses the minimal period of the inserted / deleted bases and the
    longest exact run of that unit on GRCh38 that touches the event (labels' FASTA, GRCh38 no-alt analysis set).
- PoN-tagged truth: `eval_bed/truth.INDEL.pon.vcf.gz` (FILTER PanelOfNormals), joined to the truth records in the BED.
- Check: for every call set, the overall TPb / FN / TPc / FP equal report.json `pon.at_threshold` (44 of 44 in the 2026-10-15 18:05 run).

## Files

| file | content |
|---|---|
| `e17a_indel_strata.py` | the analysis (login node; about 15 s now, 40 MB) |
| `collect.sh` | reruns it on whatever exists (`--slurm` submits it as one small job instead) |
| `tables.md` | overall, by type, by length, by context, INS / DEL x context |
| `strata.tsv` | every call set x axis x class: TPb, FN, recall, TPc, FP, precision, PoN-tagged truth, recall bound; also the type x context and type x length crosses |
| `summary.json` | the same counts per call set, with the threshold and checks |
| `e17a_final_models.py` | the same strata for the final INDEL models (release / tumor-unseen, three seeds pooled) and both comparators |
| `final_models.md` / `.tsv` / `.json` | its tables (also FP in HP>=7, retained recall, test-tuned best F1, ceiling), per call set x axis x class, and per group |

## Results (all call sets, 2026-10-15 18:00; 44 INDEL call sets)

`tables.md` was generated at 18:05, after the last call set had finished (17:50). It holds 44 INDEL call sets:

- **in-sample, 30.** HG008T Illumina / PacBio / ONT × GRCh38-only (linear) / HPRC v1.1 d9 / HPRC v2.1 d46, and
  COLO829T Illumina / ONT / fiber-seq × linear / d46. Each appears twice: at the standard threshold and at `[valbest]`.
  There is no COLO829T d9 single-set run (not local).
- **COLO829T fiber-seq joint d9 model, 4.** Trained on P100 + P75 + P50 + P25 and scored on each purity. It is not a
  matched arm.
- **cross-sample, 10, all Illumina d46.** The HG008T model -> COLO829T, H1437T, H2009T, HCC1395T, HCC1954T, and the
  COLO829T model -> HG008T and the four CASTLE lines.

**Thresholds.** Every number here is at a fixed PASS threshold:

- standard = validation recall 0.9 (0.005-0.120);
- `[valbest]` = the threshold chosen on the chr2-22 validation blocks (0.116-0.400).

The chr1 curve-best F1, which is the metric used to compare INDEL arms (report.json `pon.best`), is not in these tables.
E17a shows where the calls work. It does not rank the arms.

**1. Truth sets and the PoN recall bound** (chr1 truth in the BED; `tables.md` §1, §4):

| sample (chr1, in the BED) | truth INDELs | PoN-tagged | recall bound | HP>=7 (tagged) | HP4-6 | STR (tagged) | no context (tagged) |
|---|---:|---:|---:|---:|---:|---:|---:|
| HG008T | 587 | 379 | 0.354 | 435 (294) | 7 | 65 (54) | 80 (31) |
| COLO829T | 131 | 81 | 0.382 | 70 (42) | 6 | 27 (24) | 28 (12) |
| H1437T | 158 | 16 | 0.899 | 38 (9) | 24 | 12 (7) | 84 (0) |
| H2009T | 376 | 35 | 0.907 | 78 (25) | 78 | 23 (9) | 197 (1) |
| HCC1395T | 133 | 14 | 0.895 | 24 (12) | 18 | 4 (2) | 87 (0) |
| HCC1954T | 167 | 29 | 0.826 | 77 (23) | 32 | 16 (6) | 42 (0) |

- On HG008T and COLO829T, the PoN alone caps INDEL recall at 0.354 / 0.382, mostly through HP>=7 tracts and STRs.
- The CASTLE lines have bounds of 0.83-0.91 and mostly INDELs outside repeats. HCC1395T has 59 deletions > 10 bp out of
  133.

**2. In-sample, overall** (`tables.md` §1; F1 with P / R in brackets):

| call set (chr1, in the BED) | linear | d9 | d46 | valbest: linear / d9 / d46 F1 |
|---|---|---|---|---|
| HG008T Illumina | 0.133 (0.089 / 0.262) | 0.241 (0.244 / 0.239) | 0.244 (0.248 / 0.240) | 0.080 / 0.216 / 0.205 |
| HG008T PacBio | 0.125 (0.082 / 0.261) | 0.089 (0.054 / 0.256) | 0.095 (0.058 / 0.257) | 0.036 / 0.225 / 0.198 |
| HG008T ONT | 0.022 (0.012 / 0.239) | 0.015 (0.008 / 0.228) | 0.015 (0.007 / 0.237) | 0.047 / 0.116 / 0.151 |
| COLO829T Illumina | 0.032 (0.017 / 0.221) | not local | 0.109 (0.070 / 0.244) | 0.030 / - / 0.044 |
| COLO829T ONT | 0.004 (0.002 / 0.282) | not local | 0.002 (0.001 / 0.229) | 0 TP / - / 0.027 |
| COLO829T fiber-seq | 0.040 (0.021 / 0.298) | not local | 0.020 (0.011 / 0.260) | 0 calls / - / 0.090 |
| COLO829T fiber-seq, joint d9 (P100 / P75 / P50 / P25) | | 0.007 / 0.005 / 0.004 / 0.003 | | |

- **Recall barely depends on the arm.** At the standard threshold, recall is 0.22-0.30 in every in-sample call set,
  which is 58-78 % of the PoN bound. The arms differ in precision.
- **Illumina: the pangenome arms cut the FPs by 3.6-3.9×.** HG008T: 1,577 FPs (linear) against 434 (d9) and 428 (d46).
  COLO829T: 1,642 against 422 (d46). The cut is mostly outside homopolymers:
  - HG008T STR calls: 546 (linear), 38 (d9), 26 (d46);
  - HG008T calls with no repeat context: 403 / 95 / 119;
  - HG008T HP>=7 calls: 704 / 406 / 389;
  - COLO829T: STR 310 -> 53, no context 865 -> 119, HP>=7 359 -> 252.
- **Long reads: HP>=7 FPs dominate at the standard threshold**, which is 0.005-0.040 for every long-read set. The
  pangenome arms have more FPs than the linear arm here:

  | call set | chr1 HP>=7 PASS calls / all PASS calls |
  |---|---|
  | HG008T ONT | d46 17,017 / 18,555; d9 15,954 / 17,387; linear 9,336 / 12,032 |
  | HG008T PacBio | d9 1,993 / 2,786; d46 1,747 / 2,603; linear 838 / 1,858 |
  | COLO829T ONT | d46 20,966 / 28,331 |
  | COLO829T fiber-seq | d46 2,582 / 3,206; linear 502 / 1,818 |

- **valbest is a 1-bp, high-precision call set.**
  - HG008T d9 / d46: precision 0.39-0.41 (ONT 0.22-0.27).
  - The Illumina and PacBio valbest sets are (almost) all 1 bp: 210 / 210, 205 / 205, 235 / 238 and 199 / 199 calls.
  - It beats the standard threshold on long reads: HG008T PacBio d9 0.225 vs 0.089, d46 0.198 vs 0.095; ONT d9 0.116 vs
    0.015, d46 0.151 vs 0.015; COLO829T fiber-seq d46 0.090 vs 0.020; ONT d46 0.027 vs 0.002.
  - It is worse on Illumina: HG008T d9 0.216 vs 0.241, d46 0.205 vs 0.244.
  - With the linear arm, valbest collapses: recall 0.000-0.044, and 0 TP for COLO829T ONT and fiber-seq. COLO829T
    valbest is nearly empty in every arm (0-46 calls).
- **Joint d9 fiber-seq model.** Its recall (0.229-0.275) is like that of the single-set models. At its threshold of
  0.0073 it has 9,879 (P100) to 19,537 (P25) FPs. Of its PASS calls, 7,887 / 9,915 (P100) and 16,372 / 19,567 (P25)
  are in HP>=7. The single-set d46 fiber-seq model has 3,172 FPs on P100.

**3. Strata, in-sample at the standard threshold** (`tables.md` §2-5):

- **1-bp INDELs carry most of the recovered truth.** Recall by length:

  | length | chr1 recall |
  |---|---|
  | 1 bp | 0.30-0.39 (joint model 0.35-0.42) |
  | 2-4 bp | 0.05-0.19 |
  | 5-10 bp | 0.10-0.20 |
  | > 10 bp | 0.09-0.47 |

- **The linear arm recalls more 2-4 bp INDELs, at far lower precision.**
  - HG008T 2-4 bp recall: 0.13-0.14 (linear) against 0.07-0.08 (d9 / d46). COLO829T ONT / fiber-seq: 0.16 / 0.19
    against 0.05 / 0.08.
  - HG008T Illumina 2-4 bp precision: 0.02 on 723 calls (linear) against 0.08 / 0.10 on 108 / 99 calls.
- **STR INDELs are not recovered, and the PoN removes most of them first.**
  - STR recall is 0.00-0.12 and STR precision 0.00-0.03 in every in-sample call set.
  - The PoN tags 54 / 65 HG008T and 24 / 27 COLO829T STR truth INDELs.
- **HP>=7 recall is 0.21-0.29**, against a PoN bound inside HP>=7 of 0.324 (HG008T) and 0.400 (COLO829T).
  - HP>=7 precision of HG008T Illumina: d9 0.25, d46 0.26, linear 0.16.
  - On long reads it is 0.00-0.06 in the pangenome arms and 0.00-0.14 in the linear arm.
- **INDELs outside repeats are recalled best**: 0.29-0.57.
  - HG008T Illumina precision on them: d9 0.31, d46 0.28, linear 0.07.

**4. Cross-sample, all Illumina d46** (`tables.md` §1, §3, §4):

| model -> target (chr1, in the BED) | P | R | F1 | TP / FP / FN | recall bound |
|---|---:|---:|---:|---|---:|
| HG008T -> COLO829T | 0.045 | 0.260 | 0.076 | 34 / 730 / 97 | 0.382 |
| HG008T -> H1437T | 0.029 | 0.753 | 0.055 | 119 / 4,017 / 39 | 0.899 |
| HG008T -> H2009T | 0.064 | 0.809 | 0.118 | 304 / 4,472 / 72 | 0.907 |
| HG008T -> HCC1395T | 0.134 | 0.248 | 0.174 | 33 / 214 / 100 | 0.895 |
| HG008T -> HCC1954T | 0.021 | 0.563 | 0.041 | 94 / 4,321 / 73 | 0.826 |
| COLO829T -> HG008T | 0.326 | 0.208 | 0.254 | 122 / 252 / 465 | 0.354 |
| COLO829T -> H1437T | 0.115 | 0.424 | 0.181 | 67 / 517 / 91 | 0.899 |
| COLO829T -> H2009T | 0.205 | 0.412 | 0.274 | 155 / 601 / 221 | 0.907 |
| COLO829T -> HCC1395T | 0.368 | 0.263 | 0.307 | 35 / 60 / 98 | 0.895 |
| COLO829T -> HCC1954T | 0.072 | 0.317 | 0.117 | 53 / 687 / 114 | 0.826 |

- **HG008T model (threshold 0.1065) on H1437T / H2009T / HCC1954T: recall 0.56-0.81, but most of its calls are HP>=7
  FPs.**
  - HP>=7 calls: 3,891 / 4,136, 4,387 / 4,776 and 4,060 / 4,415 of its PASS calls, at HP>=7 precision 0.01.
  - On COLO829T, 416 of its 764 PASS calls are HP>=7 (precision 0.05) and 120 are STR (0.01).
  - HCC1395T is the exception: 247 calls, precision 0.134.
- **COLO829T model (threshold 0.0222): fewer calls.** Recall is 0.26-0.42 and precision 0.07-0.37.
  - On the CASTLE INDELs outside repeats: recall 0.29-0.49, precision 0.12-0.74.
  - On HCC1395T it recalls 34 of the 59 deletions > 10 bp at precision 0.79. The HG008T model recalls 6 of 59.
  - On HG008T, at these fixed thresholds, its F1 is 0.254 (P 0.326, R 0.208), against 0.244 for HG008T's own d46 model
    (0.248 / 0.240). Its higher threshold trades recall for precision. The curve-best comparison is not in this table.

## Still missing

- **COLO829T d9 single-set INDEL runs** are not local. COLO829T d9 is covered only by the fiber-seq joint model, which
  is not a matched arm.
- **Comparators** (ClairS-TO / DeepSomatic) are not yet re-scored. The `../E16_bootstrap/comparators.tsv` slot is
  empty.
- **Seed repeats** are pending GPU, so there is one seed per model.
- **Cross-sample INDEL** is Illumina only: there are no long-read `xs__` INDEL runs.
- **Chr1 curve-best F1** is not given per stratum; the strata use fixed thresholds.

## Caveats

1. **Fixed thresholds.** The strata need a fixed call set, so every row uses either the recall-0.9 threshold or
   valbest. The INDEL arm comparison itself uses the chr1 curve-best F1 (report.json `pon.best`).
2. **Small strata.** Many cells have fewer than 20 truth records (e.g. HG008T HP4-6: 7; COLO829T HP4-6: 6; HCC1395T
   STR: 4). E16 gives block-bootstrap intervals for the overall numbers; per-stratum intervals are not computed.
3. **The PoN is part of the result.** Recall is after the PoN. The PoN-tagged column shows how much of each stratum any
   caller loses to it (E17c / R2.6 treat the PoN choice itself).
4. **Unplaced PASS calls** (off-reference, no GRCh38 position) are outside rtg and therefore outside these strata.
5. **One seed** per model; seed repeats are GPU jobs (R3.9).

## Rerun when more runs finish

```bash
bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E17a_indel_strata/collect.sh
```

All three call-level analyses (E16, E17a, E18) at once:

```bash
cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments && \
  for e in E16_bootstrap E17a_indel_strata E18_sbs96; do bash $e/collect.sh; done
```

It processes whatever exists. When the user's ClairS-TO / DeepSomatic VCFs have been re-scored with
`../E8_comparators/score_comparator.py`, add their `eval_bed/pass_pon` directories to
`../E16_bootstrap/comparators.tsv` and rerun.
