# E18: SBS96 mutational spectra of the truth and of Pansoma's calls (R3.5)

Reviewer point: **R3.5** (do the calls reproduce the tumours' mutational signatures?). Also R2.m (minor: SBS96) and,
through the cross-sample rows, R1.7 (does a model carry its training tumour's signature?).

## Question

On chr1 inside each benchmark BED, compare the SBS96 spectrum of the truth SNVs with the spectra of Pansoma's PASS
calls after the PoN (true positives, false positives) and of the missed truth (false negatives). Then ask:

- Do the call sets show the expected hallmarks: UV (C>T at TpC / CpC, CC>TT doublets) for COLO829T, and clock-like
  SBS1 / SBS5 / SBS40 features for HG008T?
- Do the false positives look like the tumour's own mutations or like germline variants?

## Method

- Call sets: every finished SNV `vcf_chr1*/` under `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/`
  (discovery: `../E16_bootstrap/callsets.py`): in-sample runs of the three reference arms, `__norescue`, cross-sample
  `xs__`, the COLO829T fiber-seq joint model; comparators from `../E16_bootstrap/comparators.tsv` (slot; empty).
- Records: rtg's `eval_bed/vcfeval/pon_calls/`.
  - truth = tp-baseline + fn (every truth SNV in the BED; checked to be the same count for every call set of a sample);
  - TP = `tp.vcf.gz` (call side), FP = `fp.vcf.gz`, FN = `fn.vcf.gz`.
- Context: the trinucleotide around POS on the labels' GRCh38 FASTA (no-alt analysis set), pyrimidine-centred. A G or A
  reference base is reverse-complemented together with its context. The chr1 truth VCFs contain no MNP records, so all
  truth SNVs are single-base.
- Germline reference spectrum: the set's germline truth (labels manifest; dipcall of HG008-N and COLO829BL, the CASTLE /
  SEQC2 germline VCF of the HCC1395, HCC1954, H1437 and H2009 normals) on chr1 in the same BED. FILTER PASS / . / GAP1 / GAP2 (the labels'
  germline rule), single-base ALTs carried by the genotype.
- Signatures: COSMIC v3.4 SBS GRCh38 from the installed SigProfilerAssignment package
  (`/wanglab/jshen/anaconda3/lib/python3.11/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/`).
  Nothing was downloaded.
  - Reported: cosine to SBS1, SBS5, SBS7a, SBS7b, SBS18 and SBS40a, and the best three single signatures.
  - No signature refitting was done.
- Hallmarks:
  - UV share = (T[C>T]N + C[C>T]N) / all SNVs;
  - CpG share = N[C>T]G / all SNVs;
  - CC>TT doublets = adjacent SNV pairs (POS, POS+1) with REF CC -> TT or GG -> AA. These are DBS, outside SBS96.
- FP truth-like share: non-negative least squares of the FP spectrum (as proportions) on the sample's truth and germline
  spectra. It is the truth weight's share, and the fit cosine is reported with it. It can only separate the two where
  they differ. cos(truth, germline) is 0.428 for COLO829T, 0.585 for HG008T, 0.520 for H1437T, 0.469 for H2009T, 0.524
  for HCC1395T and 0.802 for HCC1954T, so the share is barely informative on HCC1954T.
- Recall by substitution class: TPb / (TPb + FN), with the PoN-tagged truth of each class (`eval_bed/truth.SNV.pon.vcf.gz`).

## Files

| file | content |
|---|---|
| `e18_sbs96.py` | the analysis (login node; about 25 s, 0.9 GB) |
| `collect.sh` | reruns it on whatever exists (`--slurm` submits it as one small job instead) |
| `tables.md` | §1 truth and germline per sample; §2 calls vs truth (cosines, FP truth-like share, hallmarks); §3 recall by class; §4 cosines to named COSMIC signatures; §5 comparators |
| `sbs96.tsv` | every spectrum (truth, germline, TP / FP / FN per call set), 96 channels in COSMIC order |
| `summary.json` | all of the above per call set |

The germline spectra are cached in `tmp/reviewer_cpu_20261013/E18_sbs96/germline_cache.json`.

## Results (all call sets, 2026-10-15 18:00; 58 SNV call sets)

`tables.md` was generated at 18:06, after the last call set had finished (17:50). It holds 58 SNV call sets:

- **in-sample, 24.**
  - HG008T Illumina / PacBio / ONT × GRCh38-only (linear) / HPRC v1.1 d9 / HPRC v2.1 d46;
  - COLO829T Illumina / ONT / fiber-seq × linear / d46;
  - a `__norescue` copy of every d9 / d46 run. Its placed calls, and so all its spectra, are identical to the standard
    run's.
  - There is no COLO829T d9 single-set run (not local).
- **COLO829T fiber-seq joint d9 model, 4.** P100 / P75 / P50 / P25. It is not a matched arm.
- **cross-sample d46, 30.** The HG008T models -> COLO829T and the four CASTLE lines, and the COLO829T models -> HG008T
  and the four CASTLE lines, on each platform.
- The truth count of each call set matched its truth set in 58 of 58 cases.

**1. The truth spectra show the expected signatures** (`tables.md` §1; chr1, in the BED):

| sample | chr1 truth SNVs | UV share | CpG C>T | CC>TT doublets (adjacent pairs) | best COSMIC (cosine) | cos(truth, germline) |
|---|---:|---:|---:|---|---|---:|
| COLO829T | 2,443 | 58.5 % | 6.2 % | 26 (46) | SBS7a 0.917, SBS7b 0.857 | 0.428 |
| HG008T | 697 | 28.8 % | 7.3 % | 2 (19) | SBS40a 0.780, SBS40c 0.685; SBS5 0.625, SBS1 0.241 | 0.585 |
| H1437T | 3,941 | 12.8 % | 2.1 % | 0 (74) | SBS4 0.870, SBS95 0.842 | 0.520 |
| H2009T | 11,106 | 18.3 % | 1.6 % | 6 (183) | SBS40a 0.754, SBS4 0.716 | 0.469 |
| HCC1395T | 3,440 | 24.6 % | 3.3 % | 15 (39) | SBS40a 0.737, SBS13 0.685 | 0.524 |
| HCC1954T | 1,450 | 19.9 % | 8.6 % | 5 (9) | SBS40a 0.843, SBS40c 0.833, SBS5 0.821 | 0.802 |

- COLO829T is UV (SBS7a), and the two lung lines H1437T and H2009T carry tobacco SBS4.
- HG008T matches no single signature well (best 0.780). It is a flat SBS5 / SBS40-like mixture, and its CpG C>T (7.3 %)
  is below its own germline's (12.6 %).
- The germline spectra are alike in all six samples: best match SBS5 (0.912-0.914), UV share 15.8-15.9 %, CpG share
  12.6-12.9 %.

**2. In-sample calls by reference arm** (§2, §3; PASS after the PoN, chr1 in the BED):

| call set (chr1, in the BED) | arm | TP / FP | cos TP+FP to truth | UV share TP / FP | recall UV C>T | recall T>C | FP truth-like share (fit cos) | cos FP-germline |
|---|---|---|---:|---|---:|---:|---|---:|
| COLO829T Illumina | linear | 1,423 / 2,221 | 0.947 | 65.7 / 29.1 % | 0.654 | 0.228 | 58 % (0.788) | 0.528 |
| COLO829T Illumina | d46 | 1,745 / 268 | 0.996 | 67.6 / 42.5 % | 0.826 | 0.332 | 57 % (0.850) | 0.577 |
| COLO829T ONT | linear | 1,325 / 2,575 | 0.913 | 64.7 / 36.6 % | 0.600 | 0.250 | 53 % (0.779) | 0.552 |
| COLO829T ONT | d46 | 1,690 / 457 | 0.996 | 69.1 / 41.1 % | 0.817 | 0.272 | 61 % (0.956) | 0.622 |
| COLO829T fiber-seq | linear | 1,525 / 521 | 0.994 | 65.9 / 37.6 % | 0.703 | 0.337 | 62 % (0.921) | 0.596 |
| COLO829T fiber-seq | d46 | 1,728 / 268 | 0.997 | 68.5 / 44.0 % | 0.828 | 0.397 | 63 % (0.875) | 0.558 |
| COLO829T fiber-seq P100 | joint d9 | 1,707 / 231 | 0.996 | 68.8 / 47.6 % | 0.822 | 0.326 | 69 % (0.896) | 0.536 |
| HG008T Illumina | linear | 548 / 1,089 | 0.923 | 31.2 / 19.6 % | 0.851 | 0.403 | 34 % (0.899) | 0.844 |
| HG008T Illumina | d9 | 579 / 702 | 0.949 | 30.9 / 20.2 % | 0.891 | 0.613 | 28 % (0.923) | 0.885 |
| HG008T Illumina | d46 | 594 / 792 | 0.940 | 30.8 / 18.8 % | 0.910 | 0.629 | 25 % (0.941) | 0.910 |
| HG008T ONT | linear | 587 / 904 | 0.933 | 30.7 / 20.0 % | 0.896 | 0.661 | 26 % (0.932) | 0.897 |
| HG008T ONT | d9 (*) | 576 / 660 | 0.951 | 31.2 / 21.1 % | 0.896 | 0.613 | 26 % (0.928) | 0.894 |
| HG008T ONT | d46 | 599 / 788 | 0.941 | 31.1 / 21.1 % | 0.925 | 0.645 | 26 % (0.946) | 0.913 |
| HG008T PacBio | linear | 576 / 969 | 0.932 | 31.4 / 19.6 % | 0.900 | 0.484 | 28 % (0.927) | 0.888 |
| HG008T PacBio | d9 (*) | 577 / 583 | 0.958 | 31.4 / 24.0 % | 0.900 | 0.468 | 34 % (0.900) | 0.845 |
| HG008T PacBio | d46 | 582 / 650 | 0.941 | 31.3 / 21.2 % | 0.905 | 0.613 | 24 % (0.934) | 0.905 |

(*) d9 checkpoint of a different recipe (1 GPU; E3 caveat 1).

- **Every arm reproduces the truth spectrum.** The cosine of the TPs alone is 0.982-0.995. With the FPs included
  (TP + FP), it is lowest in the linear arm, which has the most FPs:

  | arm | COLO829T, chr1 cos TP+FP | HG008T, chr1 cos TP+FP |
  |---|---|---|
  | d46 | 0.996-0.997 | 0.940-0.941 |
  | d9 | - | 0.949-0.958 |
  | linear | 0.913-0.994 | 0.923-0.933 |

- **UV hallmark in COLO829T.** The UV share of the TPs is 67.6-69.1 % at d46, against 58.5 % in the truth.
  - Recall of UV-context C>T at d46 is 0.817-0.828, against 0.692-0.714 overall (E16).
  - The linear arm loses this class most: 0.600-0.703, against 0.542-0.624 overall.
  - By class, the linear arm loses C>T (0.582-0.678 against 0.786-0.800 at d46) and T>A (0.378-0.417 against
    0.472-0.638). C>A and C>G change little.
- **CC>TT doublets.** The TPs contain 20 / 21 / 21 of the 26 truth doublets at d46, 18 / 16 / 20 in the linear arm and
  21 for the joint model. Doublets among the FPs:
  - none in any d9 / d46 in-sample call set;
  - 3 in COLO829T ONT linear and 1 in fiber-seq linear;
  - one each in 9 cross-sample call sets (3 on HCC1395T, all 6 on HCC1954T).
- **The weakest class is T>C.** In-sample T>C recall:
  - COLO829T: linear 0.228 / 0.250 / 0.337, d46 0.332 / 0.272 / 0.397 (Illumina / ONT / fiber-seq);
  - HG008T: 0.403-0.661 in all arms.
  - The PoN tags 13 of HG008T's 62 T>C truth SNVs and 10 of COLO829T's 184.

**3. False positives** (§2, §4):

- **HG008T: the FPs look like germline variants in every arm.**
  - Their cosine to the HG008-N germline spectrum is 0.844-0.913.
  - Their best COSMIC match is SBS5 in 7 of 9 call sets. The exceptions are Illumina linear (SBS40a 0.856, SBS5 0.851)
    and PacBio d9 (SBS40a 0.810).
  - The FP truth-like share is 24-34 %.
  - So the linear arm's extra FPs (1,089 / 904 / 969, against 792 / 788 / 650 at d46) are germline-like as well. This
    fits E9d (placed FPs after the PoN are mostly germline).
- **COLO829T, pangenome arms: about 60 % of the FP spectrum looks like the tumour's UV spectrum.**
  - At d46: FP UV share 41.1-44.0 %, between the truth (58.5 %) and the germline (15.8 %). Truth-like share 57-63 %
    (fit cosine 0.850-0.956).
  - The joint d9 model on P100: 47.6 % and 69 %.
  - Some of these FPs may be real UV mutations missing from the truth set, for example subclonal ones. This analysis
    does not prove it.
- **COLO829T, linear arm: the extra FPs are neither UV nor germline.**
  - The linear Illumina / ONT FPs (2,221 / 2,575) have a UV share of 29.1 / 36.6 %.
  - The truth + germline fit describes them poorly: fit cosine 0.788 / 0.779 (fiber-seq: 0.921).
  - Their best COSMIC match is a flat SBS40a (0.821 / 0.802).
  - The linear Illumina FPs have a CpG share of 3.8 %, against 12.7 % in the germline.
- **Joint d9 model at lower purity.** The FP truth-like share is 54 % (P25), 35 % (P50) and 34 % (P75). The cosine of
  the FPs to the germline is 0.641 / 0.768 / 0.788. The TP UV share stays at 69.1-72.6 %, and the UV C>T recall is
  0.770-0.794.

**4. Cross-sample: the COLO829T models carry a UV prior; the HG008T models do not drift** (§2, §3):

| target (chr1 truth UV share) | TP cosine to the target truth: COLO829T models (Illumina / ONT / fiber-seq) | HG008T models (Illumina / ONT / PacBio) | TP UV share: COLO829T models | HG008T models |
|---|---|---|---|---|
| HG008T (28.8 %) | 0.793 / 0.879 / 0.917 | in-sample | 55.6 / 48.8 / 44.2 % | - |
| COLO829T (58.5 %) | in-sample | 0.994 / 0.996 / 0.996 | - | 69.7 / 68.1 / 67.3 % |
| H1437T (12.8 %) | 0.874 / 0.975 / 0.964 | 0.997 / 0.994 / 0.995 | 24.8 / 16.3 / 17.5 % | 13.0 / 13.6 / 13.2 % |
| H2009T (18.3 %) | 0.852 / 0.973 / 0.974 | 0.999 / 0.998 / 0.999 | 39.1 / 24.1 / 24.2 % | 18.3 / 18.7 / 18.2 % |
| HCC1395T (24.6 %) | 0.849 / 0.972 / 0.978 | 0.994 / 0.988 / 0.986 | 40.2 / 34.0 / 30.9 % | 26.7 / 26.5 / 26.5 % |
| HCC1954T (19.9 %) | 0.904 / 0.953 / 0.966 | 0.984 / 0.969 / 0.967 | 30.9 / 25.1 / 23.0 % | 21.0 / 23.1 / 21.8 % |

- **The HG008T models' TPs match each target's truth**: cosine 0.967-0.999.
  - On the CASTLE lines, the UV share of their TPs is within 3.2 points of the truth.
  - On COLO829T it is 67.3-69.7 %, the same as COLO829T's own d46 models (67.6-69.1 %).
- **The COLO829T models over-select UV context everywhere.** The Illumina model drifts most (cosine 0.793-0.904).
  Recall of UV-context C>T against T>C for that model:

  | target (chr1, in the BED) | UV C>T | T>C | overall recall |
  |---|---:|---:|---:|
  | HG008T | 0.493 | 0.032 | 0.255 |
  | H1437T | 0.589 | 0.168 | 0.305 |
  | H2009T | 0.635 | 0.106 | 0.296 |
  | HCC1395T | 0.256 | 0.112 | 0.157 |
  | HCC1954T | 0.713 | 0.368 | 0.460 |

- **The COLO829T models' FPs on other tumours are UV-like.**
  - FP UV share: 41.8-53.0 % on HG008T and 26.6-54.0 % on the CASTLE lines.
  - On HG008T the truth-like share is 72-96 %, but the fit is poor (cosine 0.760-0.784), so that share only means
    "more UV than germline", not "like HG008T".
  - The HG008T models' FPs on other tumours are germline-like: cos FP-germline 0.665-0.936. Their best match is SBS5
    in 12 of 15 call sets, and SBS40a for the three on COLO829T.
- So a model trained on one UV-dominated tumour ranks UV-context C>T higher on every other tumour. This supports
  training on several tumours (R1.7; the `--norm dataset` / LOLO models).

## Still missing

- **COLO829T d9 single-set SNV runs** are not local. COLO829T d9 is covered only by the fiber-seq joint model, which is
  not a matched arm.
- **Comparators** (ClairS-TO / DeepSomatic) are not yet re-scored. The `../E16_bootstrap/comparators.tsv` slot is
  empty, so §5 is empty.
- **Seed repeats** are pending GPU, so there is one seed per model. The HG008T PacBio / ONT d9 checkpoints are of a
  different recipe (marked (*)).
- **Multi-tumour models.** Cross-sample `__normds` / LOLO call sets are not in yet. They would test whether training on
  several tumours removes the UV prior.

## Caveats

1. **Small counts.** HG008T has 697 truth SNVs on chr1 and about 600 TPs per call set over 96 channels. No intervals
   on the cosines are computed, so small differences between HG008T call sets (for example the cosines of the d9 and
   d46 arms) are not interpreted here.
2. **The cosines are spectrum-level.** The FP truth-like share is a two-component fit, not a classification of
   individual calls. It is uninformative where the truth and germline spectra are alike (HCC1954T, 0.802) or where the
   fit cosine is low. COSMIC GRCh38 signatures are normalised to the whole-genome trinucleotide content; this analysis
   counts chr1 inside the BED.
3. **One seed** per model; unplaced PASS calls (no GRCh38 position) are not in rtg and not in these spectra.

## Rerun when more runs finish

```bash
bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E18_sbs96/collect.sh
```

All three call-level analyses (E16, E17a, E18) at once:

```bash
cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments && \
  for e in E16_bootstrap E17a_indel_strata E18_sbs96; do bash $e/collect.sh; done
```

**Comparator slot.** When the user's ClairS-TO / DeepSomatic VCFs have been re-scored with
`../E8_comparators/score_comparator.py`, add one row per call set to `../E16_bootstrap/comparators.tsv`. The row gives
the label, the set, the kind, `<out>/eval_bed/pass_pon`, `<out>/eval_nobed/pass_pon` and `<out>/truth.SNV.pon.vcf.gz`.
Then rerun `collect.sh`. The comparators then appear in §2-4 with role `comparator`.
