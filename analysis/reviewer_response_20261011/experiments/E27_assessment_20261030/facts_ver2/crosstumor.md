# Fact sheet: R1.7, R1.m1, R2.4, R3.10 (+ cross-tumor part of R3.9)

Compiled 2026-10-30. Every number has a source. Status tags: **EXISTS** (already in a repo or run file),
**COMPUTED NOW** (computed today from existing run outputs; scripts and outputs are under
`/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/crosstumor/`, which git ignores), **MISSING** (no result exists).
No tracked file was edited.

Abbreviations: `R` = `/scratch/jshen/data/pansoma_net_v2_runs`; `RR` =
`/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011`; `CT` =
`/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/crosstumor`; `FIN` =
`/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/final_20261023`.

Model names used below:
- **m4-SNV / m6-SNV**: `R/multicancer4_*_SNV` and `R/multicancer6_*_SNV` (E12b). One seed (seed 0). Trained on HCC1395T
  (without chr1 and chr20) and HCC1954T, H1437T and H2009T on **all autosomes, chr1 and chr20 included** (each train.log
  line 2 lists chr1–22; only HCC1395T chr1/chr20 are left out, plus HG008T/COLO829T chr1 in m6). m6 adds the chr2–22 of
  HG008T and COLO829T. PASS = validation recall-0.9 threshold. `--norm dataset` with the f344bca fix ("normfix").
- **m4-INDEL-final ("unseen") / m6-INDEL-final ("release")**: `R/unseen_s{0,1,2}_*_INDEL` and
  `R/release_s{0,1,2}_*_INDEL`, the pre-registered final INDEL batch (`FIN/PREREGISTRATION.md` §2). Three seeds each.
  **chr1 is held out of every training tumor** (and chr20 of HCC1395T). Release = the six tumors, including the chr2–22
  of HG008T and COLO829T. PASS = ppr8 (Amendment 1). Best epoch = best post-PoN validation F1. The model has five site
  features as inputs. Repeat INDELs without a truth allele, and candidates outside the confident region inside the
  somatic BED, are trained as non.
- **m4/m6-INDEL-dev**: `R/multicancer{4,6}_*_INDEL` (E12b INDEL rows). An older recipe with one seed and the CASTLE
  chr1 in training. Superseded by the final batch; do not use for the letter.

Exposure labels used in the tables: **[W]** = within-tumor (other chromosomes of the target tumor are in training);
**[H]** = held-out tumor (the target tumor was never in training); **[T+]** = target added (a multi-tumor model that
includes the target tumor's chr2–22, i.e. also within-tumor).

---

## Q1. R1.7 training and validation example counts per class

**Most complete table: E21 plus one new file.**
- Single-tumor and joint-purity runs (42): `RR/experiments/E21_training_composition/runs_training.tsv` (and `.md`,
  `tables.md` T1). **EXISTS**
- Multi-tumor E12b runs (12, per set and per run): `RR/experiments/E21_training_composition/multitumor_training.tsv` and
  `capacity_multitumor.md` §B. **EXISTS**
- **Final INDEL runs (18; these are the INDEL models used in R1.6 and R3.10): `CT/final_indel_counts.tsv`.** It is parsed
  from each `R/{release,unseen}_s*_*_INDEL/train.log` by `CT/final_indel_counts.py`. **COMPUTED NOW.** E21 does not cover
  these runs.

Validation "non" always includes the off-reference tensors that have no truth match (train.log: "validation includes N
off-reference tensors ... as non"). Training does not use them.

### Single-tumor HPRC v2.1 (d46) models (E21 runs_training.tsv; train non / somatic / germline | val non / somatic / germline) — EXISTS

| model | train | val |
|---|---|---|
| HG008T Illumina SNV | 1,619,483 / 7,274 / 256,583 | 80,578 / 357 / 12,553 |
| HG008T PacBio SNV | 71,573 / 7,227 / 270,175 | 4,509 / 350 / 12,631 |
| HG008T ONT SNV | 183,510 / 7,259 / 302,622 | 9,488 / 355 / 12,816 |
| COLO829T Illumina SNV | 803,764 / 34,110 / 280,542 | 40,055 / 1,737 / 13,382 |
| COLO829T Fiber-seq SNV | 152,721 / 33,998 / 320,431 | 8,016 / 1,725 / 13,912 |
| COLO829T ONT SNV | 476,300 / 33,779 / 355,824 | 29,549 / 1,717 / 14,494 |
| HG008T Illumina INDEL | 20,551 / 3,383 / 55,260 | 1,767 / 162 / 2,622 |
| HG008T PacBio INDEL | 1,329,154 / 3,761 / 162,889 | 70,488 / 181 / 7,511 |
| HG008T ONT INDEL | 1,710,450 / 3,249 / 176,041 | 91,378 / 156 / 8,081 |
| COLO829T Illumina INDEL | 23,386 / 1,048 / 61,110 | 1,928 / 52 / 2,957 |
| COLO829T Fiber-seq INDEL | 1,149,976 / 1,141 / 180,805 | 63,255 / 58 / 8,477 |
| COLO829T ONT INDEL | 3,289,717 / 1,018 / 212,429 | 173,128 / 53 / 10,104 |

### Compact summary (ranges)

| group | train somatic | val somatic | train germline | train non | source |
|---|---|---|---|---|---|
| single-tumor SNV, d9 + d46 | 7,098–34,110 | 350–2,066 | 256,583–381,159 | 67,735–1,620,829 | E21 (EXISTS) |
| single-tumor INDEL, d9 + d46 | 971–3,761 | 52–194 | 55,260–220,610 | 20,551–3,289,717 | E21 (EXISTS) |
| joint COLO829T 4-purity SNV (d9) | 124,070–125,916 | 7,614–7,713 | 1,356,321–1,665,260 | 421,331–3,047,060 | E21 (EXISTS) |
| joint COLO829T 4-purity INDEL (d9) | 3,439–3,965 | 233–263 | 297,626–897,514 | 127,515–10,987,344 | E21 (EXISTS) |
| m4-SNV | 237,577–240,469 | 11,592–11,752 | 1,376,494–1,590,228 | 222,944–6,491,750 | E21 multitumor (EXISTS) |
| m6-SNV | 278,615–281,853 | 13,664–13,829 | 1,914,065–2,248,674 | 447,238–8,914,997 | E21 multitumor (EXISTS) |
| m4/m6-INDEL-dev (superseded) | 10,390–10,584 / 14,809–15,486 | 467–469 / 676–708 | 348,777–885,703 / 465,147–1,274,173 | 19,628–5,670,343 / 63,565–10,670,510 | E21 multitumor (EXISTS) |
| **m4-INDEL-final (unseen), 9 runs** | **9,717–9,930** | **431–462** | 324,181–827,146 | 408,782–10,465,888 | CT/final_indel_counts.tsv (COMPUTED NOW) |
| **m6-INDEL-final (release), 9 runs** | **14,132–14,854** | **620–682** | 440,263–1,215,616 | 467,482–15,606,962 | same (COMPUTED NOW) |

Per-set training counts of the final INDEL runs (seed 0; identical sets in the other seeds up to the validation
blocks), as somatic / germline / non:
- Illumina: HCC1395T 1,142 / 46,593 / 23,024; HCC1954T 1,823 / 99,985 / 128,189; H1437T 2,376 / 92,498 / 128,188;
  H2009T 4,378 / 86,073 / 130,949; release also HG008T 3,383 / 55,260 / 23,063 and COLO829T 1,048 / 61,110 / 36,685.
- PacBio and ONT: see the tsv.
- The non counts are much larger than in the E12b dev runs: for example, HCC1395T Illumina non is 23,024 here against
  3,341 in E12b, because of `--non-reasons`.

Other facts for this table:
- Partial somatic matches left out: the final INDEL runs drop 768–4,711 training tensors (tsv `left_out_train`).
- Final INDEL best epochs: unseen 3–12, release 1–8, out of 12 (tsv `best_epoch`).
- Parameters: 50.4 M (SNV) and 8.7 M (INDEL); the final INDEL logs also print 8.7 M (`train.log` line "model ...
  'site_features': 5, 8.7 M parameters").
- Joint-purity training sets are counted per dataset: one truth allele appears up to four times. For example, the
  Illumina joint SNV model's 125,916 positives are 33,797 distinct alleles (E21 README §Summary).

---

## Q2. Exposure-labelled results table (chr1 BED, after the shared population filter)

### SNV (P / R / F1). Source: `RR/experiments/E12b_multicancer_models/README.md` and `tables.md` (EXISTS)

- m4 and m6 are one seed each.
- The "single" columns use the seed-0 d46 models (E12b "in-sample" column). Cross-tumor single-tumor numbers come from
  E12a (old `--norm train` checkpoints).
- The 3-seed means of the in-sample single-tumor models are in R3.9a.

| target | platform | single-tumor own [W] | single-tumor other [H] | m4-SNV | m6-SNV | DeepSomatic-TO | ClairS-TO |
|---|---|---|---|---|---|---|---|
| HG008T | Illumina | 0.429/0.852/0.570 | COLO: 0.603/0.255/0.359 | [H] 0.543/0.684/0.605 | [T+] 0.592/0.611/0.601 | [H] 0.376/0.889/0.529 | [H] 0.410/0.845/0.552 |
| HG008T | PacBio | 0.472/0.835/0.603 | COLO fs: 0.567/0.548/0.557 | [H] 0.535/0.720/0.614 | [T+] 0.586/0.670/0.625 | 0.425/0.806/0.557 | 0.455/0.844/0.591 |
| HG008T | ONT | 0.432/0.859/0.575 | COLO: 0.579/0.420/0.487 | [H] 0.492/0.726/0.587 | [T+] 0.635/0.551/0.590 | 0.310/0.869/0.457 | 0.414/0.831/0.552 |
| COLO829T | Illumina | 0.867/0.714/0.783 | HG008: 0.792/0.640/0.708 | [H] 0.827/0.564/0.671 | [T+] 0.880/0.661/0.755 | [H] 0.427/0.835/0.565 | [H] 0.431/0.764/0.551 |
| COLO829T | Fiber-seq | 0.866/0.707/0.779 | HG008: 0.817/0.701/0.754 | [H] 0.892/0.446/0.595 | [T+] 0.880/0.691/0.774 | 0.732/0.555/0.632 | 0.779/0.745/0.762 |
| COLO829T | ONT | 0.787/0.692/0.736 | HG008: 0.678/0.589/0.630 | [H] 0.497/0.308/0.380 | [T+] 0.817/0.670/0.737 | 0.582/0.695/0.633 | 0.749/0.684/0.715 |
| HCC1395T chr1 | Illumina | MISSING (no HCC1395 in-sample model) | HG008: 0.801/0.778/0.790; COLO: 0.934/0.157/0.269 | [W] 0.877/0.816/0.845 | [W] 0.879/0.822/0.850 | [W] 0.781/0.927/0.848 | [H] 0.806/0.859/0.832 |
| HCC1395T chr1 | PacBio | MISSING | 0.792/0.624/0.698; 0.868/0.646/0.741 | [W] 0.856/0.787/0.820 | [W] 0.843/0.809/0.825 | [W] 0.766/0.927/0.839 | [H] 0.805/0.822/0.814 |
| HCC1395T chr1 | ONT | MISSING | 0.763/0.628/0.689; 0.783/0.515/0.621 | [W] 0.650/0.668/0.658 | [W] 0.767/0.697/0.730 | [W] 0.755/0.766/0.760 | [H] 0.791/0.730/0.760 |
| HCC1395T chr20 | Ill / PB / ONT | – | – | [W] F1 0.797 / 0.775 / 0.690 | [W] 0.809 / 0.788 / 0.721 | **training data** 0.845 / 0.835 / 0.752 | [H] 0.730 / 0.795 / 0.747 |

Reading: on COLO829T, never seen by m4, m4-SNV is **below both comparators on Fiber-seq and ONT** (0.595 vs
0.632 / 0.762; 0.380 vs 0.633 / 0.715). It is above them only on Illumina. On HG008T, never seen, m4 is above both
comparators on all three platforms (one seed, no bootstrap: **MISSING**).

### INDEL, final batch (ppr8, 3-seed mean P / R / F1; seed F1 range). COMPUTED NOW: `CT/final_indel_prf.tsv` (script `CT/final_indel_prf.py`)

The script reads the `report.json` files listed in `FIN/results_ppr8.json`. The F1 values and bootstraps equal
`FIN/results_ppr8.md`, which EXISTS.

| target | platform | m4-INDEL-final (unseen) | m6-INDEL-final (release) | DeepSomatic-TO | ClairS-TO | truth n |
|---|---|---|---|---|---|---|
| HG008T chr1 | Illumina | [H] 0.524/0.070/0.123 (0.113–0.137) | [T+] 0.405/0.132/0.199 (0.171–0.224) | [H] 0.393/0.141/0.208 | [H] 0.379/0.080/0.132 | 587 |
| HG008T chr1 | PacBio | [H] 0.479/0.084/0.143 (0.135–0.154) | [T+] 0.475/0.169/0.250 (0.239–0.263) | 0.198/0.031/0.053 | 0.449/0.075/0.129 | 587 |
| HG008T chr1 | ONT | [H] 0.424/0.086/0.143 (0.139–0.146) | [T+] 0.185/0.155/0.168 (0.157–0.178) | 0.172/0.019/0.034 | 0.411/0.075/0.127 | 587 |
| COLO829T chr1 | Illumina | [H] 0.439/0.104/0.168 (0.163–0.179) | [T+] 0.264/0.135/0.174 (0.159–0.190) | 0.198/0.183/0.191 | 0.106/0.145/0.122 | 131 |
| COLO829T chr1 | Fiber-seq | [H] 0.187/0.104/0.134 (0.115–0.153) | [T+] 0.317/0.117/0.170 (0.154–0.194) | 0.314/0.084/0.133 | 0.304/0.107/0.158 | 131 |
| COLO829T chr1 | ONT | [H] 0.050/0.120/0.070 (0.054–0.085) | [T+] 0.064/0.160/0.090 (0.074–0.109) | 0.032/0.015/0.021 | 0.119/0.092/0.103 | 131 |
| HCC1395T chr1 | Illumina | [W] 0.698/0.682/0.690 (0.657–0.710) | [W] 0.518/0.729/0.603 (0.534–0.660) | [W] 0.500/0.767/0.605 | [H] 0.597/0.692/0.641 | 133 |
| HCC1395T chr1 | PacBio | [W] 0.612/0.657/0.634 | [W] 0.553/0.684/0.611 | 0.607/0.617/0.612 | 0.622/0.519/0.566 | 133 |
| HCC1395T chr1 | ONT | [W] 0.392/0.617/0.478 (0.429–0.503) | [W] 0.201/0.634/0.304 (0.270–0.340) | 0.460/0.481/0.471 | 0.583/0.556/0.569 | 133 |
| HCC1395T chr20 | Ill / PB / ONT | [W] F1 0.747 / 0.613 / 0.413 | [W] 0.553 / 0.564 / 0.278 | **training data** 0.624 / 0.733 / 0.707 | [H] 0.565 / 0.548 / 0.558 | 42 |
| **HG008T chr2 (pre-registered confirmation)** | Ill / PB / ONT | [H] F1 0.173 / 0.181 / 0.160 | not tested (chr2 is in its training) | 0.275 / 0.090 / 0.058 | 0.196 / 0.185 / 0.164 | 771 |
| **COLO829T chr2 (confirmation)** | Ill / PB / ONT | [H] F1 0.260 / 0.243 / 0.141 | – | 0.239 / 0.197 / 0.092 | 0.272 / 0.254 / 0.210 | 188 |

The confirmation verdict (`FIN/results_ppr8.md`) is **success**:
- 6 of 6 cells are on par.
- The family mean difference is +0.034 [+0.016, +0.053] against DeepSomatic-TO and −0.021 [−0.039, −0.001] against
  ClairS-TO.
- Pansoma is below ClairS-TO in all six chr2 cells (−0.004 to −0.069) and above DeepSomatic-TO in 5 of 6. HG008T
  Illumina is the exception: −0.102 [−0.129, −0.075].
- Recall denominators are the full truth in the BED; population-tagged truth stays in the denominator.

**Are the release INDEL models in R1.6 / R3.10 the m6 type? Yes. EXISTS.** `FIN/PREREGISTRATION.md` §2 defines release
as the six tumors including HG008T and COLO829T, with chr1 of every tumor (and HCC1395T chr20) held out. R1.6's ranges
reproduce from the release rows: HG008T 0.111–0.177 and COLO829T 0.107–0.160 (per-seed recall in
`CT/final_indel_prf.tsv`). R3.10's "0.09–0.25" also comes from the release rows.

The R3.10 INDEL statements are therefore **[T+], target-added results**. For held-out tumors, use the unseen rows:
- chr1 (development cells): F1 0.070–0.168.
- chr2 (confirmation): F1 0.141–0.260.

---

## Q3. Seeds (EXISTS)

- **m4-SNV: one seed** (seed 0; `R/multicancer4_{Illumina,PacBio,ONT}_SNV/args.json`; E12b README "one seed", Caveats
  "One training run each"). The same holds for m6-SNV. No paired bootstrap exists for any m4/m6 SNV result (E16 has none;
  **MISSING**).
- **m4-INDEL: three seeds** (s0, s1, s2) for the final "unseen" models (`FIN/PREREGISTRATION.md` §2). The seed also
  draws the validation blocks and the threshold, so seed spread means joint training and split variation.
- The superseded E12b multicancer4_*_INDEL models have one seed. Development runs `mc4feat_*` and `mc4featS1_Illumina_INDEL`
  (seed 1) exist but are recipe development runs.

---

## Q4. R1.m1 HCC1395 comparator table and training exposure

SNV chr1 and chr20 come from E12b (EXISTS; table above). INDEL chr1 and chr20 come from the final batch (COMPUTED NOW
P/R; F1 in `FIN/results_ppr8.md`; table above). There is no in-sample single-tumor HCC1395 model (**MISSING**; E12a
Caveat 5).

Exposure to HCC1395, with sources in `RR/experiments/E8b_literature_facts/README.md`:
- **Pansoma m4 / m6 (SNV and INDEL)**: HCC1395T chr2–19, 21 and 22 are in training; chr1 and chr20 are held out. The
  HCC1395 results are therefore **within-tumor (held-out chromosome)**.
- **DeepSomatic-TO**: trained on HCC1395 (SEQC2 truth) with chr1 held out (D10, D14).
  - D2: "chr2–chr20 train, chr21–chr22 tune, chr1 held out, for every model in the paper".
  - **HCC1395 chr20 is therefore DeepSomatic training data**, and its chr20 numbers are in-sample. FIN/PREREGISTRATION
    §3.1 also calls chr20 "DeepSomatic-TO training data".
  - Whether the released 1.9 TO checkpoints also contain the other CASTLE lines is **[unverified]** (D16). The paper's TO
    models used HCC1395 only (D10). Only the tumor–normal multi-cancer models are documented with HCC1954, H1437 and
    H2009 (D4).
- **ClairS-TO (ssrs)**: HCC1395 is excluded entirely (C3: "left out the entire HCC1395 cell line"). The model is
  fine-tuned on HCC1937, HCC1954, H1437 and H2009 (C3), so it is held-out on HCC1395 chr1 and chr20.
  - The held-out chromosome for its other lines is disputed: README chr20 versus paper chr1 (C4, [unverified]).
  - COLO829 is not in its training (C3), and HG008 is absent from the ClairS-TO paper (C5).
- **DeepSomatic on HG008**: not trained on it; HG008 is used as a held-out evaluation line (D6).

---

## Q5. R2.4: COLO829T dilutions with real somatic truth

### Setup (EXISTS)

- **Model**: one joint model per platform and kind, trained on chr2–22 of P100 + P75 + P50 + P25 together. The tumor is
  **not held out**; the dilutions share reads with P100.
- **Graph**: the **HPRC v1.1 d9 graph, not v2.1** (`analysis/benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012/README.md`
  "Data": build 9d6bd8d, `--haplotypes 90`).
- **Run**: RIS 7fc7396, one seed, with one PASS threshold (validation recall 0.9) shared by all purities.
- **Tensor depth**: 200 rows in every dilution set, in the low-AF rebuilds and in the HapMap set (each
  `tensors/SNV/manifest.json`: shape [8, 200, 101]; checked now).
- **Measured purity** (median somatic SNV tensor AF relative to P100; `.../purity.md`): Illumina and Fiber-seq
  0.66 / 0.40 / 0.20, ONT 0.70 / 0.44 / 0.22 (nominal 75 / 50 / 25%).
- **Truth on chr1**: 2,443 SNVs and 131 INDELs at every purity; population-tagged truth is SNV 84 and INDEL 81
  (`E8_comparators/purity_comparators.md`).

### Candidate stage (truth with a tensor; ceiling). EXISTS: `purity.md`

| platform | P100 | P75 | P50 | P25 |
|---|---|---|---|---|
| Illumina SNV | 2008 (0.822) | 1936 (0.792) | 1864 (0.763) | 1748 (0.716) |
| Fiber-seq SNV | 2006 (0.821) | 1940 (0.794) | 1870 (0.765) | 1737 (0.711) |
| ONT SNV | 1998 (0.818) | 1922 (0.787) | 1838 (0.752) | 1736 (0.711) |
| Illumina / Fiber-seq / ONT INDEL | 57 / 67 / 69 | 55 / 65 / 69 | 53 / 61 / 64 | 36 / 58 / 63 (of 131) |

### Truth and candidate stage by expected truth VAF. COMPUTED NOW

Outputs: `CT/dilution_vaf.md` (the requested 4 bins) and `CT/dilution_vaf_fine.md` (5 bins); script
`CT/dilution_vaf.py`. The expected VAF is the truth VCF VAF (VAF_Ill for Illumina, VAF_PB for long reads) times the
measured purity. Cells show truth n, then tensor / filtered / no candidate. SNV:

| set | <0.05 | 0.05–0.08 | 0.08–0.10 | 0.10–0.20 | ≥0.20 |
|---|---|---|---|---|---|
| Illumina P100 (default 0.08 floor) | 303; 10/267/26 | 102; 20/81/1 | 44; 28/14/2 | 142; 120/19/3 | 1852; 1830/3/19 |
| Illumina P75 | 391; 6/340/45 | 106; 36/67/3 | 57; 39/16/2 | 121; 107/13/1 | 1768; 1748/2/18 |
| Illumina P50 | 506; 19/398/89 | 85; 30/53/2 | 45; 38/6/1 | 266; 253/10/3 | 1541; 1524/1/16 |
| Illumina P25 | 636; 18/427/191 | 62; 22/36/4 | 204; 194/8/2 | 899; 879/9/11 | 642; 635/2/5 |
| Fiber-seq P25 | 651; 24/430/197 | 50; 18/27/5 | 263; 245/16/2 | 783; 762/11/10 | 696; 688/5/3 |
| ONT P25 | 621; 10/446/165 | 77; 28/46/3 | 32; 27/4/1 | 441; 413/23/5 | 1272; 1258/11/3 |
| **Illumina P100, low-AF rebuild (floor 0.015)** | 303; **211**/77/15 | 102; 92/9/1 | 44; 41/1/2 | 142; 138/1/3 | 1852; 1830/3/19 |
| **Illumina P75 low-AF** | 391; 224/139/28 | 106; 97/6/3 | 57; 55/0/2 | 121; 119/1/1 | 1768; 1748/2/18 |
| **Illumina P50 low-AF** | 506; 247/182/77 | 85; 83/1/1 | 45; 44/0/1 | 266; 262/1/3 | 1541; 1524/1/16 |
| **Illumina P25 low-AF** | 636; 233/234/169 | 62; 56/4/2 | 204; 202/0/2 | 899; 888/0/11 | 642; 637/0/5 |

The other purities and platforms are in `CT/dilution_vaf_fine.md`. INDEL truth by expected VAF: at P100 all 131 are
≥ 0.20; at Illumina P25 the bins hold 12 / 60 / 58 / 1 (<0.05 / 0.05–0.10 / 0.10–0.20 / ≥0.20).

### Final PASS after population filtering

Overall SNV P / R / F1, all EXISTS: `.../colo829t_purity_joint_7fc7396_20261012/tables.md`, README §5 and
`RR/experiments/E14b_ont_readnorm/tables.md`.

| model (all [W]: COLO829T chr2–22 of all purities) | P100 | P75 | P50 | P25 |
|---|---|---|---|---|
| Illumina joint (has the read-source BQ cue) | 0.931/0.530/0.675 | 0.968/0.742/0.840 | 0.966/0.733/0.833 | 0.967/0.689/0.804 |
| Illumina joint, per-read BQ (cue removed; also `--norm dataset` and 8 GPUs) | 0.823/0.670/0.739 | 0.737/0.683/0.709 | 0.682/0.673/0.678 | 0.723/0.668/0.695 |
| Fiber-seq joint (no cue) | 0.881/0.698/0.779 | 0.741/0.639/0.687 | 0.715/0.621/0.665 | 0.786/0.663/0.720 |
| ONT joint (RIS 7fc7396) | 0.723/0.637/0.677 | 0.659/0.648/0.654 | 0.635/0.643/0.639 | 0.670/0.661/0.666 |
| ONT control (E14b, tequila) | 0.858/0.639/0.733 | 0.798/0.658/0.721 | 0.740/0.648/0.691 | 0.779/0.659/0.714 |
| ONT per-read BQ (E14b) | 0.804/0.659/0.724 | 0.731/0.659/0.693 | 0.690/0.643/0.666 | 0.742/0.662/0.700 |

INDEL PASS F1, then test-tuned best F1 in brackets, under the current rule (source:
`analysis/benchmark_7fc7396_20261012/indel_pon_71cb51a_20261016/tables.md`, EXISTS):

| model | P100 | P75 | P50 | P25 |
|---|---|---|---|---|
| Illumina joint | 0.129 (0.167) | 0.201 (0.312) | 0.238 (0.330) | 0.195 (0.233) |
| Illumina per-read BQ | 0.077 (0.182) | 0.070 (0.146) | 0.057 (0.133) | 0.037 (0.133) |
| Fiber-seq | 0.007 (0.069) | 0.005 (0.049) | 0.004 (0.055) | 0.003 (0.074) |
| ONT | 0.001 (0.055) | 0.001 (0.051) | 0.001 (0.070) | 0.001 (0.062) |

### PASS recall by expected VAF (COMPUTED NOW; rtg pon_calls; R (TP/truth))

Bin order: <0.05 / 0.05–0.08 / 0.08–0.10 / 0.10–0.20 / ≥0.20.
- **Fiber-seq joint**:
  - P100: 0.013 / 0.087 / 0.516 / 0.879 / 0.835
  - P75: 0.012 / 0.393 / 0.634 / 0.828 / 0.787
  - P50: 0.026 / 0.400 / 0.788 / 0.898 / 0.776
  - P25: 0.029 (19/651) / 0.300 (15/50) / 0.882 / 0.917 / 0.912
- **ONT control**:
  - P100: 0.033 / 0.125 / 0.258 / 0.567 / 0.782
  - P75: 0.018 / 0.205 / 0.412 / 0.654 / 0.834
  - P50: 0.013 / 0.125 / 0.552 / 0.633 / 0.862
  - P25: 0.011 / 0.299 / 0.656 / 0.873 / 0.924
  - The ONT per-read BQ model is within ±0.05 of the control in every bin.
- **Illumina joint, by bin: MISSING.** The call sets are only on RIS.
- **COLO829T P100 with v2.1 models** (<0.05 / 0.05–0.08):

  | platform | single-tumor own | m6-SNV | m4-SNV |
  |---|---|---|---|
  | Illumina | 0.026 / 0.176 | 0.026 / 0.167 | 0.026 / 0.176 |
  | Fiber-seq | 0.013 / 0.096 | 0.013 / 0.096 | 0.013 / 0.096 |
  | ONT | 0.046 / 0.144 | 0.046 / 0.135 | 0.046 / 0.135 |

  The multi-tumor model therefore gains nothing below VAF 0.08. m6 and m4 were **never run on the dilutions
  (MISSING)**.

### Comparators on the dilutions (EXISTS: `RR/experiments/E8_comparators/purity_comparators.md`, README §10)

Same chr1 BED and the same population filter. SNV F1 for P100 / P75 / P50 / P25:

| platform | ClairS-TO | DeepSomatic-TO |
|---|---|---|
| Illumina | 0.551 / 0.661 / 0.606 / 0.550 | 0.565 / 0.500 / 0.479 / 0.450 |
| Fiber-seq | 0.762 / 0.688 / 0.659 / 0.614 | 0.632 / 0.701 / 0.676 / 0.623 |
| ONT | 0.715 / 0.664 / 0.630 / 0.565 | 0.633 / 0.621 / 0.609 / 0.572 |

INDEL F1 is 0.010–0.191 everywhere.

Caveats from E8 §10:
- ClairS-TO Illumina P100 is version 0.4.1 with `--platform ilmn`; the mixtures used 0.4.2 with `ilmn_ssrs`.
- The comparators' Illumina mixtures are a different random subset of reads from Pansoma's.
- The ONT and Fiber-seq mixtures were made with the same scripts and source BAMs as Pansoma's.

Against the no-cue Pansoma rows:
- Fiber-seq: Pansoma is below DeepSomatic-TO at P75 / P50 (0.687 / 0.665 against 0.701 / 0.676) and above both
  comparators at P25 (0.720).
- Illumina: Pansoma's per-read-BQ model is above both comparators at P75–P25.

### E15b: the HapMap model on COLO829T. Model identity first (EXISTS)

- **Every HapMap number in the revision** (E23, R2.4 point 3, E15b) comes from `R/HapMap_Illumina_SNV_full`. It is the
  50.4 M SNV network, trained on `HapMap_Illumina_nonHG005` tensors with **200 rows** (manifest), floor 0.015 and
  node-alt 0.01, on the donor-free `mgi_non_HG005.d9` graph (H 82), chr2–22, best epoch 8.
- The graph-carried positives **were in its training**, at 2.0% of SNV positives with a tensor. The no-graph retrain was
  cancelled (E23 README Status 2026-10-29).
- **A 500-row HapMap model does not exist.** E23 README: "Rebuild at `--rows 500` not submitted yet (user)".
- The "500 reads" in R2.4 point 2 describes the **submitted manuscript's legacy model** (legacy builder
  `TENSOR_MAX_READ_ROWS = 499` plus a reference row; R2.4 evidence 2026-10-26).

### E15b status (EXISTS)

- Only **P100 SNV** was predicted: `R/HapMap_Illumina_SNV_full/test_chr1_colo829t_illumina/COLO829T_Illumina_lowaf_chr1.tensors.SNV.metrics.json`.
  It uses the low-AF rebuild (d9, floor 0.015) and is scored at tensor level, without a VCF or population filtering.
- P75 / P50 / P25 and INDEL: **MISSING**. vcf jobs 452399–452406 are held (`JobHeldAdmin`), and prediction stopped
  because the user asked for a 500-row model (`tmp/e15b_hapmap_model_colo829t_20261026/jobs.txt`, 2026-10-29 14:15).

### E15b P100 results

- **EXISTS** (metrics.json): truth-level P / R / F1 is 0.325 / 0.132 / 0.188 at the pipeline threshold (validation
  recall 0.9, p ≥ 0.9876) and 0.202 / 0.189 / 0.195 at threshold_t. The model labels 1,808 of 2,296 somatic tensors
  germline.
- **COMPUTED NOW** (`CT/dilution_vaf_fine.md` §3; `CT/e15b_precision_by_af.md`, Slurm 462877 on miao):

| | <0.05 | 0.05–0.08 | 0.08–0.10 | 0.10–0.20 | ≥0.20 |
|---|---|---|---|---|---|
| truth recall, pipeline threshold | **0.432** (131/303) | **0.755** (77/102) | 0.727 | 0.521 | **0.005** (9/1852) |
| truth recall, threshold_t | 0.647 | 0.853 | 0.886 | 0.852 | 0.010 |

- Precision by the call's own AF, at the pipeline threshold, by AF band:

  | call AF | precision (TP/calls) |
  |---|---|
  | 0–0.02 | 0.167 (2/12) |
  | 0.02–0.05 | 0.179 (96/537) |
  | 0.05–0.08 | 0.353 (109/309) |
  | 0.08–0.10 | 0.781 |
  | 0.10–0.20 | 0.949 |

  At threshold_t the same bands give 0.065 / 0.106 / 0.264 / 0.638 / 0.865.
- For comparison, the standard d46 Illumina model's PASS recall at <0.05 / 0.05–0.08 is 0.026 / 0.176.

### Conclusion from data only

- **No configuration has shown useful recall below VAF 0.08 on the real COLO829T dilutions.**
- Standard 0.08 floor: PASS recall at expected VAF < 0.05 is 0.011–0.046 and at 0.05–0.08 is 0.087–0.400 (Fiber-seq and
  ONT dilutions; P100 for all platforms).
- Low-AF rebuild: it gives candidates to 37–70% of the <0.05 truth (Illumina P100–P25), but no classifier was run on the
  diluted sets.
- The only low-AF classifier result (P100, tensor level, before population filtering) recalls 0.43–0.75 below 0.08 at
  precision 0.18–0.35, and it misses almost all clonal SNVs (recall 0.005 at VAF ≥ 0.20).
- The multi-tumor model (R2.4 point 4) has the same floor and the same <0.08 recall as the single-tumor models.

---

## Q6. R3.10 INDEL strata (final release models, ppr8, chr1 BED, after population filtering)

- **EXISTS**: `RR/experiments/E17a_indel_strata/final_models.md` (§1 overall, §2 type, §3 length, §4 context, §5 type ×
  context) and `final_models.tsv`. It covers release, unseen, DeepSomatic-TO and ClairS-TO.
- **COMPUTED NOW**: `CT/e17a_compact.md` (script `CT/e17a_compact.py`). It gives, per cell, truth N (population-tagged),
  TP, FP and FN as seed means, plus P, R and F1, for all four call sets.

Letter-sized version, **release** models. Cells show truth N (population-tagged) and R / P.

| stratum | HG008T Ill | HG008T PacBio | HG008T ONT | COLO Ill | COLO Fiber-seq | COLO ONT |
|---|---|---|---|---|---|---|
| all | 587 (384): 0.13 / 0.40 | 0.17 / 0.47 | 0.16 / 0.18 | 131 (81): 0.13 / 0.24 | 0.12 / 0.32 | 0.16 / 0.06 |
| INS | 368 (248): 0.10 / 0.50 | 0.14 / 0.53 | 0.13 / 0.17 | 77 (52): 0.11 / 0.41 | 0.08 / 0.38 | 0.16 / 0.09 |
| DEL | 219 (136): 0.19 / 0.35 | 0.21 / 0.42 | 0.19 / 0.20 | 54 (29): 0.17 / 0.18 | 0.17 / 0.29 | 0.16 / 0.04 |
| length 1 | 320 (188): 0.19 / 0.42 | 0.25 / 0.47 | 0.22 / 0.18 | 69 (37): 0.19 / 0.25 | 0.15 / 0.28 | 0.21 / 0.06 |
| length 2–4 | 135 (102): 0.04 / 0.22 | 0.06 / 0.37 | 0.06 / 0.13 | 37 (28): 0.04 / 0.14 | 0.02 / 0.25 | 0.05 / 0.04 |
| length 5–10 | 62 (47): 0.09 / 0.52 | 0.10 / 0.64 | 0.10 / 0.23 | 10 (8): 0.03 / 0.07 | 0.10 / 0.60 | 0.10 / 0.08 |
| length >10 | 70 (47): 0.10 / 0.47 | 0.09 / 0.70 | 0.10 / 0.30 | 15 (8): 0.18 / 0.40 | 0.20 / 0.45 | 0.22 / 0.19 |
| HP ≥ 7 | 435 (294): 0.09 / 0.44 | 0.13 / 0.46 | 0.11 / 0.14 | 70 (42): 0.09 / 0.23 | 0.04 / 0.22 | 0.11 / 0.07 |
| HP 4–6 | 7 (0): 0.90 / 0.33 | 1.00 / 0.46 | 0.95 / 0.26 | 6 (3): 0.39 / 0.24 | 0.44 / 0.24 | 0.33 / 0.02 |
| STR | 65 (57): 0.01 / 0.09 | 0.02 / 0.25 | 0.02 / 0.04 | 27 (24): 0.01 / 0.14 | 0.00 / – | 0.04 / 0.04 |
| no repeat | 80 (33): 0.38 / 0.41 | 0.42 / 0.52 | 0.44 / 0.30 | 28 (12): 0.31 / 0.25 | 0.35 / 0.41 | 0.36 / 0.09 |

Counts behind the table:
- Mean TP / FP / FN: HG008T 77.7 / 114.3 / 509.3 (Illumina), 99.3 / 110.0 / 487.7 (PacBio), 91.0 / 410.7 / 496.0 (ONT);
  COLO829T 17.7 / 56.0 / 113.3, 15.3 / 33.3 / 115.7 and 21.0 / 324.7 / 110.0.
- FP in HP ≥ 7, as a share of all FP: HG008T 50.7 / 114.3 (44%) on Illumina, 67.0 / 110.0 (61%) on PacBio and
  286.3 / 410.7 (70%) on ONT; COLO829T 38%, 32% and 32%. COLO829T ONT also has 99.0 of its FP in HP 4–6.
- Population-tagged truth is concentrated in HP ≥ 7 (294 of 435) and STR (57 of 65) on HG008T.

Comparators, all-INDEL R / P (E17a §1):

| | HG008T Ill | HG008T PacBio | HG008T ONT | COLO Ill | COLO Fiber-seq | COLO ONT |
|---|---|---|---|---|---|---|
| DeepSomatic-TO | 0.141 / 0.393 | 0.031 / 0.198 | 0.019 / 0.172 | 0.183 / 0.198 | 0.084 / 0.314 | 0.015 / 0.032 |
| ClairS-TO | 0.080 / 0.379 | 0.075 / 0.449 | 0.075 / 0.411 | 0.145 / 0.106 | 0.107 / 0.304 | 0.092 / 0.119 |

The strata are in `CT/e17a_compact.md`. These are [T+] models; the unseen ([H]) strata are in the same file.

---

## Q7. Conflicts in the current letter (`.../assemble/response_letter_20261015.md`), with corrected values

### R1.7

1. "multi-tumor models see 237,577–281,853 and 10,390–15,486" (INDEL): the INDEL range comes from the superseded E12b
   dev models. The final INDEL models (R1.6 / R3.10) see **9,717–9,930 (unseen) and 14,132–14,854 (release)** somatic
   training examples, with 431–462 / 620–682 validation examples (`CT/final_indel_counts.tsv`). The INDEL models reported
   in R1.6 / R3.10 are multi-tumor, not single-tumor (971–3,761 refers to single-tumor models).
2. "Training uses chr2–22 (HCC1395: chr2–19, 21, 22); ... chr1 is used in neither": this is **false for m4/m6-SNV** (and
   the E12b INDEL models). HCC1954T, H1437T and H2009T were trained on all autosomes, chr1 and chr20 included
   (`R/multicancer*_SNV/train.log` line 2; E12b Design). It is true for the final INDEL models.
3. "Not used: ... for the four added lines, INDEL candidates without a truth allele in homopolymers ≥ 7 bp or STR" and
   "Other examples ... inside the benchmark region": the final INDEL models **train those repeat INDELs as non**
   (`--non-reasons repeat_indel_without_truth`). They also train out-of-confident-region candidates inside the somatic
   BED as non, and they take five site features as input (`FIN/PREREGISTRATION.md` §1). The letter does not mention
   these.
4. Early stopping "best validation somatic F1", "multi-tumor ... 10–19 and 6–9" epochs, and PASS "recall-0.9": the
   final INDEL models use the best **post-population-filtering** validation F1, best epochs **1–12**, and PASS **ppr8**.
5. "Multi-tumor models (one run each)": this holds only for SNV. The final INDEL multi-tumor models have three seeds.
   Multi-tumor INDEL results are absent from R1.7 / R1.m1 (Q2 table).
6. "reaches SNV F1 0.605 / 0.614 / 0.587 on HG008T, above HG008T's own models": that comparison is one seed against one
   seed. Against the 3-seed means of HG008T's own v2.1 models (0.573 / 0.606 / 0.596, R3.9a), **ONT 0.587 is below the
   mean**; the seed range is 0.576–0.607. No CI exists.
7. The own-model PacBio value is written 0.605 (E12a / RIS) in one sentence and 0.603 (E12b CPU re-prediction) in another
   source. Use one value.
8. "On COLO829T ... it reaches 0.671 / 0.595 / 0.380" (m4): add that this is **below both comparators on Fiber-seq and
   ONT** (DeepSomatic-TO 0.632 / 0.633, ClairS-TO 0.762 / 0.715).

### R1.m1

9. "[待补: INDEL results of the multi-tumor models]": fill from the Q2 / Q4 INDEL table.
   - HCC1395 chr1, m4-final: 0.690 / 0.634 / 0.478, against DeepSomatic-TO 0.605 / 0.612 / 0.471 and ClairS-TO
     0.641 / 0.566 / 0.569.
   - HCC1395 chr1, m6-final: 0.603 / 0.611 / 0.304.
10. "these lines are in both comparators' training data": this is verified for ClairS-TO ssrs (C3). It is **not verified
    for DeepSomatic-TO**: the paper's TO models used HCC1395 only (D10), and the released checkpoints are unverified
    (D16).
11. Any HCC1395 chr20 comparison must say that **chr20 is DeepSomatic training data** (D2). ClairS-TO never saw HCC1395
    (C3).

### R2.4

12. Point 2 ("deeper tensors, up to 500 reads instead of 200"): this is the legacy submitted model. **Every HapMap / E15b
    number in the revision is from the 200-row `HapMap_Illumina_SNV_full`, and no 500-row model exists** (E23 README).
    The model identity is inconsistent with point 3.
13. Point 1 says recall "is zero below 0.05, where no truth SNV had a candidate". This holds only under E2.4's mixed VAF
    definition: a truth with a tensor gets its tensor AF, which is ≥ 0.08 by construction. With truth VCF VAF, 5–20 of
    303–307 truth SNVs below 0.05 have a tensor at P100, and the d46 models recall **0.013–0.046** (4–14 of 303–307)
    (`CT/dilution_vaf.md`).
14. Point 4 ("we recommend the multi-tumor model"): m6 has the same floor and the same <0.08 recall as the single-tumor
    models (Q5). It was never tested on dilutions, and it includes COLO829T chr2–22. It is not a low-VAF answer.
15. There is no dilution table (purity, truth n, candidate recall, P/R/F1, same-tumor training); it is "[待补]" and
    points to R1.5. R1.5 points back to R2.4 ("Results by purity and allele fraction are in R2.4"; "its behaviour on
    somatic variants is assessed on the COLO829T dilutions (R2.4)"). This is a **circular reference**. E15b exists only
    for P100 SNV at tensor level.
16. The graph of the dilution models (**HPRC v1.1 d9**, not v2.1) is not stated.

### R3.10

17. Point 2, INDEL claim:
    - The release models are **[T+], with the target tumor's chr2–22 in training**, on development cells.
    - Of the "higher than both" cells, only HG008T PacBio and ONT have CIs that exclude zero against both comparators.
      COLO829T Fiber-seq has +0.038 [−0.035, +0.111] and +0.012 [−0.041, +0.068].
    - The "lower" Illumina cells are not significant against DeepSomatic-TO: −0.009 [−0.027, +0.010] and
      −0.016 [−0.081, +0.044].
    - The held-out-tumor pre-registered result is missing from the letter: unseen models on chr2, Pansoma below
      ClairS-TO in 6 / 6 cells, family −0.021 [−0.039, −0.001]; above DeepSomatic-TO, +0.034 [+0.016, +0.053]
      (`FIN/results_ppr8.md`).
18. Point 2, SNV claim ("exceeds ... in all six tested settings"): this rests on in-sample single-tumor models. Against
    ClairS-TO, 4 of 6 CIs include zero (R3.9b). With the held-out m4, Pansoma is below both comparators on COLO829T
    Fiber-seq and ONT.
19. Point 3 ("recall reported by truth VAF (R2.4)"): R2.4 has per-VAF recall only for undiluted COLO829T. The
    per-purity × VAF recall is in `CT/dilution_vaf*.md` (new) and E2.4 §3 (Fiber-seq only).
20. Point 4 ("too few true INDELs ... (R1.6)"): R1.6 has no dilution INDEL numbers. The source is the purity benchmark:
    131 truth INDELs per purity and 4–34 TP at the best points (README §4).
21. The assessment's "R1.5 base-quality controls" citation is **not in the current R3.10 text**. Base-quality centring is
    now in R1.5 (and done: Illumina SNV and INDEL, ONT SNV).

### R3.9, cross-tumor part

22. "INDELs at the recall-0.9 operating point, Pansoma trades higher recall for lower precision" (R3.9b INDEL rows): these
    rows are from superseded single-tumor d46 INDEL models.
    - The final release models at ppr8 often have **comparable or higher precision**. For example, HG008T PacBio:
      P 0.475 against 0.198 / 0.449.
    - The final batch has its own paired bootstrap (`FIN/results_ppr8.md`).
    - "Cross-tumor HG008T ... one training seed" is true for SNV (E12a, m4-SNV); the final INDEL unseen models have
      three seeds.
