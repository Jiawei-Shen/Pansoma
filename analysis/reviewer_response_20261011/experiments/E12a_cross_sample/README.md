# E12a: cross-sample application of the d46 models (preliminary)

> **2026-10-16, INDEL PoN rule changed** (user decision: gnomAD / CoLoRSdb exact allele, AF >= 0.01 in homopolymers >= 8 bp, else AF >= 1e-4; was >= 7 bp and 1e-3). The 44 local chr1 INDEL call sets were re-scored in place (`tmp/reviewer_cpu_20261013/pon_hp8_20261016/rescore.sh`, array 418042; the old evaluations are kept as `eval_bed.pon_hp7_1e-3` / `eval_nobed.pon_hp7_1e-3`) and this experiment's collector was rerun (job 418071), so `tables.md` / `summary.json` / TSVs now use the new rule. INDEL numbers quoted in the prose below may still be those of the earlier rule (changes are small: e.g. HG008T Illumina v1.1 PASS F1 0.241 -> 0.251, best 0.265 -> 0.270; HG008T chr1 truth INDELs tagged 413 -> 418 of 635). SNV results are unchanged.


Reviewer points: **R1.7** (training scale and tumor diversity: does a model trained on one tumor generalize?),
**R1.m1** (additional benchmarks: HCC1395), **R3.9** (stability), and the in-sample caveat behind every Pansoma number:
the RIS models are trained on chr2-22 of the same tumor they are tested on (chr1).

## What was run

- Models: the 2026-10-12 RIS checkpoints on the HPRC v2.1 d46 graph:
  - `HG008_{Illumina,PacBio,ONT}_SNV_d46_7fc7396`, `COLO829_{Illumina,fiberseq,ONT}_SNV_d46_7fc7396`;
  - the two Illumina INDEL d46 models.
- Each model predicts chr1 of another tumor's d46 tensor set, on the same platform (PacBio ↔ fiberseq), and of HCC1395T
  (SEQC2 v1.2.1 truth, labels with rule (a)).
- Prediction ran on tequila's GPUs (bf16 as on RIS; arrays 387379 and 387323, `tmp/reviewer_cpu_20261013/task_gpu.sh`).
- Then `machine_learning/pansoma_net_v2/jobs/vcf.sh`:
  - graph_vcf at the model's OWN validation recall-0.9 threshold, chosen on its training tumor (no retuning on the
    target);
  - then linear_vcf, the 7fc7396 PoN, and rtg vcfeval in the target's benchmark BED and without it.
- Outputs: `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/xs__<model>__<target set>/`.
- `collect.py` writes `tables.md` and `summary.json`.

## Results (chr1, in the BED, after the PoN; `tables.md`)

SNV, PASS F1 (best F1). The in-sample model of the target is shown for comparison:

| platform | HG008T model → COLO829T | COLO829T in-sample | COLO829T model → HG008T | HG008T in-sample | HG008T model → HCC1395T | COLO829T model → HCC1395T |
|---|---|---|---|---|---|---|
| Illumina | 0.708 (0.717) | 0.783 (0.786) | 0.359 (0.554) | 0.570 (0.613) | 0.790 (0.809) | 0.269 (0.810) |
| PacBio / fiberseq | 0.754 (0.764) | 0.779 (0.783) | 0.557 (0.597) | 0.605 (0.622) | 0.698 (0.775) | 0.741 (0.819) |
| ONT | 0.630 (0.639) | 0.737 (0.739) | 0.487 (0.579) | 0.575 (0.604) | 0.689 (0.699) | 0.621 (0.691) |

Illumina INDEL: HG008T → COLO829T 0.076 (0.135); COLO829T → HG008T 0.254 (0.262), against in-sample 0.244 (0.249);
HG008T → HCC1395T 0.174 (0.247); COLO829T → HCC1395T 0.307 (0.331).

## Reading

- **The HG008T models transfer.** They lose 0.02–0.11 F1 on COLO829T against COLO829T's own models, and reach
  0.69–0.79 on HCC1395T.
  - HG008T (pancreatic) has only about 7,200 somatic SNV training tensors on d46, against about 34,000 for COLO829T, and a
    different mutational spectrum (COLO829 is UV-dominated).
  - So a model trained on one tumor is not tied to that tumor's signature.
- **The COLO829T models lose more at PASS, mostly through the threshold.**
  - Illumina on HCC1395T: recall 0.157 at precision 0.934 at the transferred threshold of 0.640, while the curve
    reaches 0.810 at its best point. On HG008T: 0.359, best 0.554.
  - The ranking transfers; the calibration of p_somatic does not.
  - These checkpoints z-score BQ / MAPQ with the training tumor's statistics. COLO829T Illumina has a median BQ of 30
    against HG008T's 4-level binned BQ.
  - That is the shift that commit f2a0e4e (was cbc94a7; `train --norm dataset`, now the default) removes for new
    models.
- **The long-read COLO829T models transfer better than the Illumina one** (fiberseq → HCC1395T PacBio 0.741, best
  0.819).

## CASTLE lines (added 2026-10-15; same procedure, array 405185)

The same d46 models applied to chr1 of the three CASTLE lines (CASTLE truth and BED of each set's labels; germline
labels QUAL >= 10). SNV PASS F1 after the PoN, in the BED (best F1):

| model (trained on) | HCC1954T | H1437T | H2009T |
|---|---|---|---|
| HG008T Illumina | 0.383 (0.384) | 0.849 (0.852) | 0.937 (0.940) |
| HG008T PacBio | 0.311 (0.369) | 0.832 (0.852) | 0.933 (0.947) |
| HG008T ONT | 0.360 (0.365) | 0.787 (0.792) | 0.884 (0.895) |
| COLO829T Illumina | 0.412 (0.431) | 0.456 (0.828) | 0.454 (0.932) |
| COLO829T Fiber-seq | 0.468 (0.474) | 0.753 (0.874) | 0.809 (0.950) |
| COLO829T ONT | 0.426 (0.426) | 0.761 (0.840) | 0.801 (0.930) |

Illumina INDEL: HG008T model 0.041 / 0.055 / 0.118, COLO829T model 0.117 / 0.181 / 0.274 (HCC1954T / H1437T / H2009T).

Reading:
- On H1437T and H2009T the HG008T models reach 0.79-0.94 without ever seeing these tumors. The COLO829T models
  again lose mostly at their transferred threshold (Illumina recall 0.30, best F1 0.83-0.93).
- HCC1954T is the exception for every model: precision 0.22-0.37 at recall 0.46-0.78 (e.g. HG008T Illumina 1,108 TP /
  3,229 FP). E22 finds no germline problem in the HCC1954 truth (no truth SNV is a normal allele; 2.39 % common gnomAD
  alleles; no normal-LOH segments), while E4b classifies 93-99 % of these false calls as the normal's own germline
  alleles. So on HCC1954T the transferred models leak the patient's private germline variants: a model limitation
  (germline / somatic separation without the training tumor's statistics), not a truth artefact. Why HCC1954 in
  particular (e.g. copy-number or allele-fraction shifts of its germline hets) is not resolved here [unverified].
- The CASTLE truth sets are DeepSomatic's training truth; Pansoma never trained on these lines, so for Pansoma they
  are fully held out. H1437T has only about 70 % tensor coverage on all three platforms: E22 finds that about 30 % of its
  truth SNVs (chr9q, chr1:9-12 Mb) are probably germline alleles the B-lymphoblastoid normal lost through LOH; on chr1
  (the test) only 19 of 4,266 H1437 truth SNVs are common alleles, so the chr1 numbers here are little affected.
- Same caveats as above: old checkpoints (training-tumor BQ / MAPQ statistics), one seed, graph-coordinate calls not
  scored (column "unplaced PASS" in `tables.md`: up to 3,972 for the COLO829T ONT model).

## Caveats

1. **Preliminary.**
   - Old checkpoints, trained with the training tumor's BQ / MAPQ statistics.
   - The final cross-sample evaluation (E12) needs models trained with `--norm dataset` on RIS, ideally on several
     tumors (leave-one-line-out once the CASTLE sets are finished).
2. **One seed, one test chromosome.** No confidence intervals yet (R3.9).
3. **Unplaced PASS calls** (off-reference, not scored by rtg) are listed per row. They are not in the F1 (E9a / E9c).
4. **HCC1395T truth (SEQC2) is short-read derived.**
   - Its INDEL truth misses many homopolymer INDELs, which is why the labels use rule (a).
   - The PoN-limited INDEL ceiling there is 0.43–0.79, against 0.25–0.27 on HG008T / COLO829T. The INDEL numbers are
     not comparable across targets.
5. **No in-sample HCC1395T model yet.** The HCC1395T sets are on this cluster; training them is a GPU job for RIS.

## Status

- Done: the 16 cross-sample predictions of 2026-10-13 and the 24 CASTLE ones of 2026-10-15, with their evaluation.
- Pending, all GPU work on RIS:
  - models with per-dataset normalization;
  - leave-one-line-out;
  - an HCC1395T in-sample model as the reference.
