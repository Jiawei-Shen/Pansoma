# E1.6: recall funnel of the chr1 somatic truth (where each missed truth variant is lost)

> **2026-10-16, INDEL PoN rule changed** (user decision: gnomAD / CoLoRSdb exact allele, AF >= 0.01 in homopolymers >= 8 bp, else AF >= 1e-4; was >= 7 bp and 1e-3). The 44 local chr1 INDEL call sets were re-scored in place (`tmp/reviewer_cpu_20261013/pon_hp8_20261016/rescore.sh`, array 418042; the old evaluations are kept as `eval_bed.pon_hp7_1e-3` / `eval_nobed.pon_hp7_1e-3`) and this experiment's collector was rerun (job 418071), so `tables.md` / `summary.json` / TSVs now use the new rule. INDEL numbers quoted in the prose below may still be those of the earlier rule (changes are small: e.g. HG008T Illumina v1.1 PASS F1 0.241 -> 0.251, best 0.265 -> 0.270; HG008T chr1 truth INDELs tagged 413 -> 418 of 635). SNV results are unchanged.


Reviewer points: **R1.6** (low sensitivity) and **R2.3** (explain the lower recall). Shared code:
`../call_level_common/calls_lib.py` (call sets, rtg files, repeat classes) and `../E2.4_vaf_bins/vaf_bins.py` (truth AF
caches, other-arm sets). Both are imported; nothing in them or in the pipeline is changed.

**Status (2026-10-15 19:51): complete for every call set that exists locally (32).**
- In-sample (30): SNV and INDEL x GRCh38-only / HPRC v1.1 d9 / HPRC v2.1 d46 for HG008T Illumina / PacBio / ONT (18);
  SNV and INDEL x GRCh38-only / d46 for COLO829T Illumina / fiberseq / ONT (12).
- Supplementary (2): COLO829T fiberseq P100 called by the joint d9 model (`COLO829_fiberseq_{SNV_joint_7fc7396_gpu4,
  INDEL_joint_7fc7396}/vcf_chr1_COLO829T_fiberseq`). It is the only COLO829T d9 call set here; the model was trained on
  P100 + P75 + P50 + P25 jointly, so it is not arm-matched.
- Not included: COLO829T d9 single-set runs (not local), cross-sample (`xs__`) and `__norescue` call sets, INDEL
  `vcf_chr1_valbest`, the P75 / P50 / P25 mixtures.

## Question

For every chr1 truth somatic variant in the benchmark BED, at which step does Pansoma lose it?

1. **No tensor**: no candidate tensor carries the truth allele as its representative allele.
2. **Tensor, not PASS**: a tensor exists, but p_somatic is below the PASS threshold (the checkpoint's validation
   threshold at recall 0.9).
3. **PASS, PoN-tagged**: called PASS, then removed by the population filter (PoN).
4. **Called**: rtg true positive after the PoN.

The funnel is summarized per call set and per arm, and split by variant kind (SNV / INS / DEL), truth VAF bin and
repeat context.

## Data

- **Call sets**: `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/`.
  - `test_chr1/<set>.tensors.<KIND>.predictions.ndjson.gz`: one record per chr1 tensor with `p_somatic`,
    `offref_call` (SNV off-reference rescue) and `truth_ids` (the somatic truth alleles a tensor labelled 1 stands
    for: representative or partial match).
  - `vcf_chr1*/eval_bed/` (written by `machine_learning/pansoma_net_v2/jobs/vcf.sh`):
    - `report.json`: threshold and rtg counts;
    - `vcfeval/pon_calls/` (after the PoN) and `vcfeval/raw_calls/` (before it): `tp-baseline` / `fn`;
    - `<set>.<KIND>.linear.pon.vcf.gz`: the calls with their PoN tags;
    - `truth.<KIND>.pon.vcf.gz`: the truth records with their PoN tags (same PoN script and rule).
  - `vcf_chr1*/<set>.<KIND>.linear.unplaced.vcf.gz`: calls without a GRCh38 position.
- **Tensor sets**: `/scratch/jshen/data/pansoma_v2_tensors/<SET>/`.
  - `truth/somatic.graph.tsv`: each truth allele and its graph node keys (the candidate IDs that would be the exact
    allele).
  - `tensors/somatic.recall.tsv`: the per-allele status behind `tensors/truth_recall.json`:
    `tensor_representative`, `tensor_non_representative_allele`, `filtered` (with the builder's reasons `min_af` /
    `min_variants`) or `no_candidate`.
- **Absorption**: E5b's `/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46/per_truth.tsv.gz`, columns
  `d9_spelled` / `d46_spelled` (a graph path from the window's first to last GRCh38 anchor spells the truth ALT).
- **Truth AF**: E2.4's caches `/scratch/jshen/data/reviewer_experiments_20261013/call_level/E2.4/truth_af/<set>.<KIND>.chr1.tsv.gz`
  (the `af` of the matched allele in `chr1_variant_summary.ndjson`).
- **Truth in the BED, chr1** (rtg baseline): HG008T (GIAB draft v0.2) 697 SNV and 587 INDEL records; COLO829T (SMaHT
  union BED) 2,443 SNV and 131 INDEL records (`tables.md` section 1).

## Method

1. **Universe**: the rtg baseline of `pon_calls`, i.e. `tp-baseline` + `fn`: the PASS truth records of the kind on chr1
   in the BED. There is one row per (call set, truth record) in
   `/scratch/jshen/data/reviewer_experiments_20261013/E1.6_recall_funnel/per_truth.tsv.gz` (29,574 rows).
2. **Stage**: the first one that applies.
   - **4 called**: in `pon_calls/tp-baseline`.
   - **3 PoN**: in `raw_calls/tp-baseline` but not in `pon_calls/tp-baseline`. A PASS call matched it, and the PoN
     removed that call.
   - **2 tensor, not PASS**: rtg FN before the PoN, and some tensor's representative allele is exactly the truth
     allele (its `candidate_id` is one of the allele's node keys). Sub-reasons:
     - `below_threshold`: no exact tensor is called, i.e. `p_somatic` < the run's threshold (`report.json`
       `threshold`); for an SNV off-reference tensor, `offref_call` is false;
     - `pass_unplaced`: called, but the call has no GRCh38 position;
     - `pass_not_matched`: called and placed, but rtg FN.
   - **1 no tensor**: no exact representative tensor. Reasons, from `somatic.recall.tsv`:
     - `af_floor`: `filtered` with `min_af`, i.e. the tumour AF upper bound is below the build AF floor 0.08;
     - `min_variants`: `filtered` for read support only;
     - `absorbed`: `no_candidate`, and the arm's graph spells the ALT (E5b). It cannot occur on GRCh38-only;
     - `no_candidate`: `no_candidate`, and the ALT is not spelled;
     - `other_allele`: `tensor_non_representative_allele`. The site has a tensor whose A1 is another allele, and the
       call VCF holds A1 only;
     - `other`: `complex_allele` or `no_unique_grch38_node_key`.

     A flag counts the stage-1 records that have a **partial-match tensor** (labelled somatic for this truth by
     overlap: another length in the same homopolymer, a residual edit) and whether that tensor is PASS. Such a call is
     another allele, so it is an rtg FP, not a match.
3. **PoN tags**: the truth record's own tag (`truth.<KIND>.pon.vcf.gz`) and the tag of its exact tensor's call
   (`linear.pon.vcf.gz`). The "untagged" funnel restricts the truth to records that the PoN does not tag, i.e. the truth
   that any caller using this PoN can still report.
4. **Splits**:
   - **Variant kind**: from `truth_labels.split_alleles` (imported).
   - **Truth VAF bin**: [0, 0.08), [0.08, 0.2), [0.2, 0.4), [0.4, 1]. The VAF source is recorded per record and counted
     in `tables.md` section 7. Sources, in order:
     1. **tensor**: the tumour read AF of the truth allele in the call set's own tensor set (exact match,
        representative first);
     2. **build_min_af**: the own set filtered the allele for `min_af`, so its AF is below 0.08 (bin 0-0.08, no value);
     3. **other_arm**: the AF in the same reads' tensor set on another reference arm;
     4. **truth_vcf**: COLO829T truth VCF `VAF_Ill` for Illumina sets, `VAF_PB` for fiberseq and ONT sets;
     5. **other_platform**: HG008T only, the median tensor AF over the other two platforms' sets;
     6. **unknown**.

     For COLO829T, every record also gets a bin from the truth VCF VAF alone (`tables.md` section 8 and 10c). That
     bin does not depend on whether a tensor exists.
   - **Repeat context**: `calls_lib.Repeats`, which is the pipeline's `RepeatContext` (indexed_gam_pipeline_v4
     `truth_labels`, imported) for INS / DEL and its `tract()` for SNVs. Groups: HP>=7 (homopolymer tract >= 7 bp), STR
     (unit 2-6, >= 3 units), other.
5. **Check against rtg** (each call set, `report.json` `at_threshold`):
   - stage 4 = `pon.tp_baseline`;
   - stages 1 + 2 + 3 = `pon.fn`;
   - stage 3 = `raw.tp_baseline` - `pon.tp_baseline`;
   - stages 1 + 2 = `raw.fn`.

## Scripts and jobs

| file | role |
|---|---|
| `e16_recall_funnel.py` | per-truth stages, reasons, VAF, repeat class; `summary.json`, `tables.md`, `per_truth.tsv.gz` |
| `collect.sh` | runs it (login node or Slurm) |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E1.6/job.sh` | Slurm wrapper (`-p general`, 2 CPUs, 4 GB); not needed so far |

Run: on the login node, `bash collect.sh`, 2026-10-15 19:51: 30 s, 765,588 KB max RSS
(`tmp/reviewer_cpu_20261013/E1.6/logs/login_collect_20261015.log`). No Slurm job was submitted.

## Results

All numbers are from `tables.md` (section numbers as there). The data are chr1, the truth is in the BED, and calls are
PASS at each run's validation recall-0.9 threshold.

### 1. Headline per arm (tables.md section 0)

Cell: `no tensor / tensor, not PASS / PASS, PoN-tagged / called`, in % of the chr1 truth records in the BED (n in
section 2). `n.l.` = no local single-set run (COLO829T d9).

| chr1 truth in BED | GRCh38-only | HPRC v1.1 d9 | HPRC v2.1 d46 |
|---|---|---|---|
| HG008T Illumina SNV | 1.9 / 14.5 / 5.0 / 78.6 | 5.3 / 8.5 / 3.2 / 83.1 | 4.4 / 6.9 / 3.4 / 85.2 |
| HG008T PacBio SNV | 2.0 / 9.8 / 5.6 / 82.6 | 5.5 / 8.9 / 2.9 / 82.8 | 5.0 / 7.9 / 3.6 / 83.5 |
| HG008T ONT SNV | 2.3 / 8.5 / 5.0 / 84.2 | 4.7 / 9.8 / 2.9 / 82.6 | 4.4 / 5.6 / 4.0 / 85.9 |
| HG008T Illumina INDEL | 19.4 / 6.8 / 47.5 / 26.2 | 59.8 / 4.3 / 12.1 / 23.9 | 56.7 / 6.0 / 13.3 / 24.0 |
| HG008T PacBio INDEL | 15.7 / 8.9 / 49.4 / 26.1 | 52.8 / 4.1 / 17.5 / 25.6 | 52.0 / 4.9 / 17.4 / 25.7 |
| HG008T ONT INDEL | 41.2 / 4.9 / 30.0 / 23.9 | 55.5 / 4.4 / 17.2 / 22.8 | 58.1 / 3.4 / 14.8 / 23.7 |
| COLO829T Illumina SNV | 16.7 / 23.4 / 1.7 / 58.2 | n.l. | 18.0 / 9.0 / 1.6 / 71.4 |
| COLO829T fiberseq SNV | 17.5 / 18.3 / 1.8 / 62.4 | n.l. | 17.8 / 9.8 / 1.6 / 70.7 |
| COLO829T ONT SNV | 17.8 / 26.6 / 1.4 / 54.2 | n.l. | 18.1 / 11.1 / 1.6 / 69.2 |
| COLO829T Illumina INDEL | 14.5 / 23.7 / 39.7 / 22.1 | n.l. | 53.4 / 7.6 / 14.5 / 24.4 |
| COLO829T fiberseq INDEL | 6.1 / 18.3 / 45.8 / 29.8 | n.l. | 48.1 / 7.6 / 18.3 / 26.0 |
| COLO829T ONT INDEL | 19.1 / 19.1 / 33.6 / 28.2 | n.l. | 57.3 / 4.6 / 15.3 / 22.9 |

For example, HG008T PacBio SNV on d9: 5.5 % no tensor, 8.9 % tensor below the threshold, 2.9 % PASS but PoN-tagged,
82.8 % called.

### 2. Reconciliation with rtg (tables.md sections 1 and 11)

The stage totals reproduce `report.json` exactly in **all 32 call sets**:
- stage 4 = TP after the PoN;
- stages 1-3 = FN after the PoN;
- stage 3 = TP before the PoN minus TP after it;
- stages 1-2 = FN before the PoN.

There is no mismatch to explain. Invariants (section 11):
- Every stage-3 truth record carries a PoN tag itself (32 of 32 call sets).
- No called (stage-4) truth record is PoN-tagged (32 of 32).
- Stage 2 is always `below_threshold` (32 of 32). No exact tensor that was called PASS ended up without a GRCh38
  position or unmatched by rtg.
- Three stage-3 records have no exact tensor: 2 in HG008T ONT INDEL d9 and 1 in HG008T ONT INDEL d46. rtg matched each
  to a partial-match tensor's PASS call, and the truth record itself is PoN-tagged.

| call set (chr1) | arm | truth | 1 no tensor | 2 tensor, not PASS | 3 PASS, PoN-tagged | 4 called | threshold | rtg TP / FN after PoN (report.json) | stages reproduce rtg |
|---|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV | GRCh38-only | 2443 | 407 (16.7 %) | 572 (23.4 %) | 41 (1.7 %) | 1423 (58.2 %) | 0.0806 | 1423 / 1020 | yes |
| COLO829T Illumina SNV | HPRC v2.1 d46 | 2443 | 440 (18.0 %) | 219 (9.0 %) | 39 (1.6 %) | 1745 (71.4 %) | 0.6395 | 1745 / 698 | yes |
| COLO829T fiberseq SNV | GRCh38-only | 2443 | 427 (17.5 %) | 448 (18.3 %) | 43 (1.8 %) | 1525 (62.4 %) | 0.0852 | 1525 / 918 | yes |
| COLO829T fiberseq SNV | HPRC v2.1 d46 | 2443 | 436 (17.8 %) | 239 (9.8 %) | 40 (1.6 %) | 1728 (70.7 %) | 0.4812 | 1728 / 715 | yes |
| COLO829T ONT SNV | GRCh38-only | 2443 | 434 (17.8 %) | 649 (26.6 %) | 35 (1.4 %) | 1325 (54.2 %) | 0.0920 | 1325 / 1118 | yes |
| COLO829T ONT SNV | HPRC v2.1 d46 | 2443 | 442 (18.1 %) | 271 (11.1 %) | 40 (1.6 %) | 1690 (69.2 %) | 0.5677 | 1690 / 753 | yes |
| HG008T Illumina SNV | GRCh38-only | 697 | 13 (1.9 %) | 101 (14.5 %) | 35 (5.0 %) | 548 (78.6 %) | 0.0640 | 548 / 149 | yes |
| HG008T Illumina SNV | HPRC v1.1 d9 | 697 | 37 (5.3 %) | 59 (8.5 %) | 22 (3.2 %) | 579 (83.1 %) | 0.2513 | 579 / 118 | yes |
| HG008T Illumina SNV | HPRC v2.1 d46 | 697 | 31 (4.4 %) | 48 (6.9 %) | 24 (3.4 %) | 594 (85.2 %) | 0.3237 | 594 / 103 | yes |
| HG008T PacBio SNV | GRCh38-only | 697 | 14 (2.0 %) | 68 (9.8 %) | 39 (5.6 %) | 576 (82.6 %) | 0.0397 | 576 / 121 | yes |
| HG008T PacBio SNV | HPRC v1.1 d9 | 697 | 38 (5.5 %) | 62 (8.9 %) | 20 (2.9 %) | 577 (82.8 %) | 0.3142 | 577 / 120 | yes |
| HG008T PacBio SNV | HPRC v2.1 d46 | 697 | 35 (5.0 %) | 55 (7.9 %) | 25 (3.6 %) | 582 (83.5 %) | 0.2912 | 582 / 115 | yes |
| HG008T ONT SNV | GRCh38-only | 697 | 16 (2.3 %) | 59 (8.5 %) | 35 (5.0 %) | 587 (84.2 %) | 0.0418 | 587 / 110 | yes |
| HG008T ONT SNV | HPRC v1.1 d9 | 697 | 33 (4.7 %) | 68 (9.8 %) | 20 (2.9 %) | 576 (82.6 %) | 0.4492 | 576 / 121 | yes |
| HG008T ONT SNV | HPRC v2.1 d46 | 697 | 31 (4.4 %) | 39 (5.6 %) | 28 (4.0 %) | 599 (85.9 %) | 0.2309 | 599 / 98 | yes |
| COLO829T Illumina INDEL | GRCh38-only | 131 | 19 (14.5 %) | 31 (23.7 %) | 52 (39.7 %) | 29 (22.1 %) | 0.0110 | 29 / 102 | yes |
| COLO829T Illumina INDEL | HPRC v2.1 d46 | 131 | 70 (53.4 %) | 10 (7.6 %) | 19 (14.5 %) | 32 (24.4 %) | 0.0222 | 32 / 99 | yes |
| COLO829T fiberseq INDEL | GRCh38-only | 131 | 8 (6.1 %) | 24 (18.3 %) | 60 (45.8 %) | 39 (29.8 %) | 0.0081 | 39 / 92 | yes |
| COLO829T fiberseq INDEL | HPRC v2.1 d46 | 131 | 63 (48.1 %) | 10 (7.6 %) | 24 (18.3 %) | 34 (26.0 %) | 0.0117 | 34 / 97 | yes |
| COLO829T ONT INDEL | GRCh38-only | 131 | 25 (19.1 %) | 25 (19.1 %) | 44 (33.6 %) | 37 (28.2 %) | 0.0050 | 37 / 94 | yes |
| COLO829T ONT INDEL | HPRC v2.1 d46 | 131 | 75 (57.3 %) | 6 (4.6 %) | 20 (15.3 %) | 30 (22.9 %) | 0.0127 | 30 / 101 | yes |
| HG008T Illumina INDEL | GRCh38-only | 587 | 114 (19.4 %) | 40 (6.8 %) | 279 (47.5 %) | 154 (26.2 %) | 0.0480 | 154 / 433 | yes |
| HG008T Illumina INDEL | HPRC v1.1 d9 | 587 | 351 (59.8 %) | 25 (4.3 %) | 71 (12.1 %) | 140 (23.9 %) | 0.1195 | 140 / 447 | yes |
| HG008T Illumina INDEL | HPRC v2.1 d46 | 587 | 333 (56.7 %) | 35 (6.0 %) | 78 (13.3 %) | 141 (24.0 %) | 0.1065 | 141 / 446 | yes |
| HG008T PacBio INDEL | GRCh38-only | 587 | 92 (15.7 %) | 52 (8.9 %) | 290 (49.4 %) | 153 (26.1 %) | 0.0367 | 153 / 434 | yes |
| HG008T PacBio INDEL | HPRC v1.1 d9 | 587 | 310 (52.8 %) | 24 (4.1 %) | 103 (17.5 %) | 150 (25.6 %) | 0.0397 | 150 / 437 | yes |
| HG008T PacBio INDEL | HPRC v2.1 d46 | 587 | 305 (52.0 %) | 29 (4.9 %) | 102 (17.4 %) | 151 (25.7 %) | 0.0270 | 151 / 436 | yes |
| HG008T ONT INDEL | GRCh38-only | 587 | 242 (41.2 %) | 29 (4.9 %) | 176 (30.0 %) | 140 (23.9 %) | 0.0084 | 140 / 447 | yes |
| HG008T ONT INDEL | HPRC v1.1 d9 | 587 | 326 (55.5 %) | 26 (4.4 %) | 101 (17.2 %) | 134 (22.8 %) | 0.0178 | 134 / 453 | yes |
| HG008T ONT INDEL | HPRC v2.1 d46 | 587 | 341 (58.1 %) | 20 (3.4 %) | 87 (14.8 %) | 139 (23.7 %) | 0.0135 | 139 / 448 | yes |
| COLO829T fiberseq SNV | HPRC v1.1 d9 (joint model) | 2443 | 437 (17.9 %) | 251 (10.3 %) | 48 (2.0 %) | 1707 (69.9 %) | 0.5492 | 1707 / 736 | yes |
| COLO829T fiberseq INDEL | HPRC v1.1 d9 (joint model) | 131 | 68 (51.9 %) | 6 (4.6 %) | 21 (16.0 %) | 36 (27.5 %) | 0.0073 | 36 / 95 | yes |

### 3. Pooled over the three platforms (tables.md sections 10 and 10a)

The pooled tables sum the per-run records, so each truth record counts once per platform.

| sample kind (chr1) | arm |  | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | all | 7329 (100.0 %) | 17.3 | 22.8 | 1.6 | 58.3 |
| COLO829T SNV | HPRC v2.1 d46 | all | 7329 (100.0 %) | 18.0 | 9.9 | 1.6 | 70.4 |
| COLO829T INDEL | GRCh38-only | all | 393 (100.0 %) | 13.2 | 20.4 | 39.7 | 26.7 |
| COLO829T INDEL | HPRC v2.1 d46 | all | 393 (100.0 %) | 52.9 | 6.6 | 16.0 | 24.4 |
| HG008T SNV | GRCh38-only | all | 2091 (100.0 %) | 2.1 | 10.9 | 5.2 | 81.8 |
| HG008T SNV | HPRC v1.1 d9 | all | 2091 (100.0 %) | 5.2 | 9.0 | 3.0 | 82.8 |
| HG008T SNV | HPRC v2.1 d46 | all | 2091 (100.0 %) | 4.6 | 6.8 | 3.7 | 84.9 |
| HG008T INDEL | GRCh38-only | all | 1761 (100.0 %) | 25.4 | 6.9 | 42.3 | 25.4 |
| HG008T INDEL | HPRC v1.1 d9 | all | 1761 (100.0 %) | 56.0 | 4.3 | 15.6 | 24.1 |
| HG008T INDEL | HPRC v2.1 d46 | all | 1761 (100.0 %) | 55.6 | 4.8 | 15.2 | 24.5 |

| sample kind (chr1) | arm | kind | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | SNV | 7329 (100.0 %) | 17.3 | 22.8 | 1.6 | 58.3 |
| COLO829T SNV | HPRC v2.1 d46 | SNV | 7329 (100.0 %) | 18.0 | 9.9 | 1.6 | 70.4 |
| COLO829T INDEL | GRCh38-only | INS | 231 (58.8 %) | 11.3 | 22.9 | 43.3 | 22.5 |
| COLO829T INDEL | GRCh38-only | DEL | 162 (41.2 %) | 16.0 | 16.7 | 34.6 | 32.7 |
| COLO829T INDEL | HPRC v2.1 d46 | INS | 231 (58.8 %) | 58.9 | 3.9 | 14.7 | 22.5 |
| COLO829T INDEL | HPRC v2.1 d46 | DEL | 162 (41.2 %) | 44.4 | 10.5 | 17.9 | 27.2 |
| HG008T SNV | GRCh38-only | SNV | 2091 (100.0 %) | 2.1 | 10.9 | 5.2 | 81.8 |
| HG008T SNV | HPRC v1.1 d9 | SNV | 2091 (100.0 %) | 5.2 | 9.0 | 3.0 | 82.8 |
| HG008T SNV | HPRC v2.1 d46 | SNV | 2091 (100.0 %) | 4.6 | 6.8 | 3.7 | 84.9 |
| HG008T INDEL | GRCh38-only | INS | 1104 (62.7 %) | 28.1 | 6.6 | 40.9 | 24.4 |
| HG008T INDEL | GRCh38-only | DEL | 657 (37.3 %) | 21.0 | 7.3 | 44.6 | 27.1 |
| HG008T INDEL | HPRC v1.1 d9 | INS | 1104 (62.7 %) | 62.3 | 2.6 | 13.3 | 21.7 |
| HG008T INDEL | HPRC v1.1 d9 | DEL | 657 (37.3 %) | 45.5 | 7.0 | 19.5 | 28.0 |
| HG008T INDEL | HPRC v2.1 d46 | INS | 1104 (62.7 %) | 60.7 | 3.9 | 13.7 | 21.7 |
| HG008T INDEL | HPRC v2.1 d46 | DEL | 657 (37.3 %) | 47.0 | 6.2 | 17.7 | 29.1 |

### 4. Why there is no tensor, the PoN, and the truth the PoN leaves (tables.md section 10e)

| sample kind (chr1) | arm | stage 1 | af_floor | min_variants | absorbed | no_candidate | other_allele | other | stage 1 with a PASS partial-match tensor | stage 1 PoN-tagged | truth PoN-tagged | untagged n | untagged: 1 % | 2 % | 3 % | 4 % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | 1268 (17.3 %) | 807 (11.0 %) | 253 (3.5 %) | 0 (0.0 %) | 208 (2.8 %) | 0 (0.0 %) | 0 (0.0 %) | 0 | 54 / 1268 | 252 / 7329 (3.4 %) | 7077 | 17.2 | 22.5 | 0.0 | 60.4 |
| COLO829T SNV | HPRC v2.1 d46 | 1318 (18.0 %) | 838 (11.4 %) | 308 (4.2 %) | 9 (0.1 %) | 163 (2.2 %) | 0 (0.0 %) | 0 (0.0 %) | 3 | 72 / 1318 | 252 / 7329 (3.4 %) | 7077 | 17.6 | 9.4 | 0.0 | 73.0 |
| COLO829T INDEL | GRCh38-only | 52 (13.2 %) | 8 (2.0 %) | 10 (2.5 %) | 0 (0.0 %) | 12 (3.1 %) | 22 (5.6 %) | 0 (0.0 %) | 10 | 43 / 52 | 243 / 393 (61.8 %) | 150 | 6.0 | 24.0 | 0.0 | 70.0 |
| COLO829T INDEL | HPRC v2.1 d46 | 208 (52.9 %) | 17 (4.3 %) | 8 (2.0 %) | 81 (20.6 %) | 83 (21.1 %) | 19 (4.8 %) | 0 (0.0 %) | 49 | 165 / 208 | 243 / 393 (61.8 %) | 150 | 28.7 | 7.3 | 0.0 | 64.0 |
| HG008T SNV | GRCh38-only | 43 (2.1 %) | 4 (0.2 %) | 16 (0.8 %) | 0 (0.0 %) | 23 (1.1 %) | 0 (0.0 %) | 0 (0.0 %) | 0 | 24 / 43 | 180 / 2091 (8.6 %) | 1911 | 1.0 | 9.5 | 0.0 | 89.5 |
| HG008T SNV | HPRC v1.1 d9 | 108 (5.2 %) | 2 (0.1 %) | 13 (0.6 %) | 35 (1.7 %) | 58 (2.8 %) | 0 (0.0 %) | 0 (0.0 %) | 4 | 70 / 108 | 180 / 2091 (8.6 %) | 1911 | 2.0 | 7.4 | 0.0 | 90.6 |
| HG008T SNV | HPRC v2.1 d46 | 97 (4.6 %) | 4 (0.2 %) | 10 (0.5 %) | 30 (1.4 %) | 53 (2.5 %) | 0 (0.0 %) | 0 (0.0 %) | 3 | 65 / 97 | 180 / 2091 (8.6 %) | 1911 | 1.7 | 5.4 | 0.0 | 92.9 |
| HG008T INDEL | GRCh38-only | 448 (25.4 %) | 62 (3.5 %) | 66 (3.7 %) | 0 (0.0 %) | 61 (3.5 %) | 259 (14.7 %) | 0 (0.0 %) | 104 | 331 / 448 | 1137 / 1761 (64.6 %) | 624 | 18.8 | 9.6 | 0.0 | 71.6 |
| HG008T INDEL | HPRC v1.1 d9 | 987 (56.0 %) | 73 (4.1 %) | 158 (9.0 %) | 275 (15.6 %) | 417 (23.7 %) | 64 (3.6 %) | 0 (0.0 %) | 254 | 815 / 987 | 1137 / 1761 (64.6 %) | 624 | 27.6 | 4.5 | 0.0 | 67.9 |
| HG008T INDEL | HPRC v2.1 d46 | 979 (55.6 %) | 55 (3.1 %) | 124 (7.0 %) | 249 (14.1 %) | 465 (26.4 %) | 86 (4.9 %) | 0 (0.0 %) | 337 | 806 / 979 | 1137 / 1761 (64.6 %) | 624 | 27.7 | 3.2 | 0.0 | 69.1 |

Per call set: `tables.md` sections 2 (stage-1 reasons), 3 (stage-2 reasons) and 4 (PoN tags per stage and the untagged
funnel).

### 5. Truth VAF (tables.md sections 10b and 10c; per call set sections 6-8)

| sample kind (chr1) | arm | VAF bin | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | 0-0.08 | 1237 (16.9 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T SNV | GRCh38-only | 0.08-0.2 | 601 (8.2 %) | 2.3 | 3.8 | 2.5 | 91.3 |
| COLO829T SNV | GRCh38-only | 0.2-0.4 | 354 (4.8 %) | 1.4 | 12.7 | 2.3 | 83.6 |
| COLO829T SNV | GRCh38-only | >0.4 | 5137 (70.1 %) | 0.2 | 31.2 | 1.9 | 66.7 |
| COLO829T SNV | HPRC v2.1 d46 | 0-0.08 | 1231 (16.8 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T SNV | HPRC v2.1 d46 | 0.08-0.2 | 601 (8.2 %) | 4.8 | 6.3 | 1.8 | 87.0 |
| COLO829T SNV | HPRC v2.1 d46 | 0.2-0.4 | 358 (4.9 %) | 3.4 | 15.6 | 1.7 | 79.3 |
| COLO829T SNV | HPRC v2.1 d46 | >0.4 | 5139 (70.1 %) | 0.9 | 12.4 | 2.0 | 84.8 |
| COLO829T INDEL | GRCh38-only | 0-0.08 | 8 (2.0 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T INDEL | GRCh38-only | 0.08-0.2 | 44 (11.2 %) | 36.4 | 27.3 | 18.2 | 18.2 |
| COLO829T INDEL | GRCh38-only | 0.2-0.4 | 119 (30.3 %) | 10.1 | 10.9 | 48.7 | 30.3 |
| COLO829T INDEL | GRCh38-only | >0.4 | 222 (56.5 %) | 7.2 | 24.8 | 40.5 | 27.5 |
| COLO829T INDEL | HPRC v2.1 d46 | 0-0.08 | 17 (4.3 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| COLO829T INDEL | HPRC v2.1 d46 | 0.08-0.2 | 59 (15.0 %) | 54.2 | 23.7 | 13.6 | 8.5 |
| COLO829T INDEL | HPRC v2.1 d46 | 0.2-0.4 | 109 (27.7 %) | 53.2 | 3.7 | 20.2 | 22.9 |
| COLO829T INDEL | HPRC v2.1 d46 | >0.4 | 208 (52.9 %) | 48.6 | 3.8 | 15.9 | 31.7 |
| HG008T SNV | GRCh38-only | 0-0.08 | 4 (0.2 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | GRCh38-only | 0.08-0.2 | 38 (1.8 %) | 23.7 | 47.4 | 15.8 | 13.2 |
| HG008T SNV | GRCh38-only | 0.2-0.4 | 245 (11.7 %) | 7.3 | 6.1 | 9.4 | 77.1 |
| HG008T SNV | GRCh38-only | >0.4 | 1801 (86.1 %) | 0.5 | 10.8 | 4.4 | 84.2 |
| HG008T SNV | GRCh38-only | unknown | 3 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v1.1 d9 | 0-0.08 | 2 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v1.1 d9 | 0.08-0.2 | 39 (1.9 %) | 66.7 | 30.8 | 0.0 | 2.6 |
| HG008T SNV | HPRC v1.1 d9 | 0.2-0.4 | 240 (11.5 %) | 13.3 | 10.8 | 3.3 | 72.5 |
| HG008T SNV | HPRC v1.1 d9 | >0.4 | 1807 (86.4 %) | 2.5 | 8.4 | 3.0 | 86.2 |
| HG008T SNV | HPRC v1.1 d9 | unknown | 3 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v2.1 d46 | 0-0.08 | 4 (0.2 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T SNV | HPRC v2.1 d46 | 0.08-0.2 | 27 (1.3 %) | 66.7 | 25.9 | 0.0 | 7.4 |
| HG008T SNV | HPRC v2.1 d46 | 0.2-0.4 | 245 (11.7 %) | 12.7 | 9.0 | 4.1 | 74.3 |
| HG008T SNV | HPRC v2.1 d46 | >0.4 | 1812 (86.7 %) | 2.3 | 6.2 | 3.7 | 87.8 |
| HG008T SNV | HPRC v2.1 d46 | unknown | 3 (0.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | GRCh38-only | 0-0.08 | 62 (3.5 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | GRCh38-only | 0.08-0.2 | 340 (19.3 %) | 55.9 | 4.4 | 31.2 | 8.5 |
| HG008T INDEL | GRCh38-only | 0.2-0.4 | 641 (36.4 %) | 22.0 | 3.9 | 49.9 | 24.2 |
| HG008T INDEL | GRCh38-only | >0.4 | 679 (38.6 %) | 2.4 | 11.9 | 47.0 | 38.7 |
| HG008T INDEL | GRCh38-only | unknown | 39 (2.2 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | HPRC v1.1 d9 | 0-0.08 | 74 (4.2 %) | 98.6 | 0.0 | 1.4 | 0.0 |
| HG008T INDEL | HPRC v1.1 d9 | 0.08-0.2 | 403 (22.9 %) | 71.0 | 10.4 | 14.9 | 3.7 |
| HG008T INDEL | HPRC v1.1 d9 | 0.2-0.4 | 607 (34.5 %) | 55.7 | 3.3 | 19.8 | 21.3 |
| HG008T INDEL | HPRC v1.1 d9 | >0.4 | 633 (35.9 %) | 38.9 | 2.1 | 14.8 | 44.2 |
| HG008T INDEL | HPRC v1.1 d9 | unknown | 44 (2.5 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | HPRC v2.1 d46 | 0-0.08 | 55 (3.1 %) | 100.0 | 0.0 | 0.0 | 0.0 |
| HG008T INDEL | HPRC v2.1 d46 | 0.08-0.2 | 397 (22.5 %) | 76.8 | 10.1 | 9.8 | 3.3 |
| HG008T INDEL | HPRC v2.1 d46 | 0.2-0.4 | 636 (36.1 %) | 53.8 | 4.9 | 19.7 | 21.7 |
| HG008T INDEL | HPRC v2.1 d46 | >0.4 | 629 (35.7 %) | 37.0 | 2.1 | 16.4 | 44.5 |
| HG008T INDEL | HPRC v2.1 d46 | unknown | 44 (2.5 %) | 100.0 | 0.0 | 0.0 | 0.0 |

| sample kind (chr1) | arm | truth-VCF VAF bin | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | 0-0.08 | 1227 (16.7 %) | 92.7 | 1.0 | 0.2 | 6.1 |
| COLO829T SNV | GRCh38-only | 0.08-0.2 | 562 (7.7 %) | 18.9 | 1.4 | 1.6 | 78.1 |
| COLO829T SNV | GRCh38-only | 0.2-0.4 | 311 (4.2 %) | 3.2 | 10.3 | 2.6 | 83.9 |
| COLO829T SNV | GRCh38-only | >0.4 | 5229 (71.3 %) | 0.3 | 30.9 | 1.9 | 66.9 |
| COLO829T SNV | HPRC v2.1 d46 | 0-0.08 | 1227 (16.7 %) | 93.3 | 0.9 | 0.2 | 5.6 |
| COLO829T SNV | HPRC v2.1 d46 | 0.08-0.2 | 562 (7.7 %) | 20.5 | 4.6 | 1.1 | 73.8 |
| COLO829T SNV | HPRC v2.1 d46 | 0.2-0.4 | 311 (4.2 %) | 3.2 | 13.2 | 2.9 | 80.7 |
| COLO829T SNV | HPRC v2.1 d46 | >0.4 | 5229 (71.3 %) | 0.9 | 12.4 | 2.0 | 84.7 |
| COLO829T INDEL | GRCh38-only | 0.2-0.4 | 77 (19.6 %) | 19.5 | 11.7 | 36.4 | 32.5 |
| COLO829T INDEL | GRCh38-only | >0.4 | 316 (80.4 %) | 11.7 | 22.5 | 40.5 | 25.3 |
| COLO829T INDEL | HPRC v2.1 d46 | 0.2-0.4 | 77 (19.6 %) | 58.4 | 10.4 | 11.7 | 19.5 |
| COLO829T INDEL | HPRC v2.1 d46 | >0.4 | 316 (80.4 %) | 51.6 | 5.7 | 17.1 | 25.6 |

### 6. Repeat context (tables.md section 10d; per call set section 9)

| sample kind (chr1) | arm | repeat | n (3 platforms; % of all) | 1 no tensor % | 2 not PASS % | 3 PoN % | 4 called % |
|---|---|---|---|---|---|---|---|
| COLO829T SNV | GRCh38-only | HP>=7 | 93 (1.3 %) | 19.4 | 47.3 | 6.5 | 26.9 |
| COLO829T SNV | GRCh38-only | STR | 381 (5.2 %) | 20.5 | 22.8 | 1.6 | 55.1 |
| COLO829T SNV | GRCh38-only | other | 6855 (93.5 %) | 17.1 | 22.4 | 1.6 | 58.9 |
| COLO829T SNV | HPRC v2.1 d46 | HP>=7 | 93 (1.3 %) | 23.7 | 30.1 | 7.5 | 38.7 |
| COLO829T SNV | HPRC v2.1 d46 | STR | 381 (5.2 %) | 21.8 | 17.3 | 0.8 | 60.1 |
| COLO829T SNV | HPRC v2.1 d46 | other | 6855 (93.5 %) | 17.7 | 9.3 | 1.6 | 71.5 |
| COLO829T INDEL | GRCh38-only | HP>=7 | 210 (53.4 %) | 14.3 | 21.9 | 39.0 | 24.8 |
| COLO829T INDEL | GRCh38-only | STR | 81 (20.6 %) | 9.9 | 18.5 | 64.2 | 7.4 |
| COLO829T INDEL | GRCh38-only | other | 102 (26.0 %) | 13.7 | 18.6 | 21.6 | 46.1 |
| COLO829T INDEL | HPRC v2.1 d46 | HP>=7 | 210 (53.4 %) | 44.3 | 7.1 | 22.4 | 26.2 |
| COLO829T INDEL | HPRC v2.1 d46 | STR | 81 (20.6 %) | 77.8 | 3.7 | 14.8 | 3.7 |
| COLO829T INDEL | HPRC v2.1 d46 | other | 102 (26.0 %) | 51.0 | 7.8 | 3.9 | 37.3 |
| HG008T SNV | GRCh38-only | HP>=7 | 102 (4.9 %) | 13.7 | 18.6 | 43.1 | 24.5 |
| HG008T SNV | GRCh38-only | STR | 147 (7.0 %) | 8.8 | 18.4 | 14.3 | 58.5 |
| HG008T SNV | GRCh38-only | other | 1842 (88.1 %) | 0.9 | 9.9 | 2.4 | 86.9 |
| HG008T SNV | HPRC v1.1 d9 | HP>=7 | 102 (4.9 %) | 46.1 | 28.4 | 4.9 | 20.6 |
| HG008T SNV | HPRC v1.1 d9 | STR | 147 (7.0 %) | 25.9 | 13.6 | 7.5 | 53.1 |
| HG008T SNV | HPRC v1.1 d9 | other | 1842 (88.1 %) | 1.2 | 7.6 | 2.5 | 88.7 |
| HG008T SNV | HPRC v2.1 d46 | HP>=7 | 102 (4.9 %) | 38.2 | 28.4 | 10.8 | 22.5 |
| HG008T SNV | HPRC v2.1 d46 | STR | 147 (7.0 %) | 23.1 | 10.2 | 8.2 | 58.5 |
| HG008T SNV | HPRC v2.1 d46 | other | 1842 (88.1 %) | 1.3 | 5.3 | 2.9 | 90.4 |
| HG008T INDEL | GRCh38-only | HP>=7 | 1305 (74.1 %) | 25.0 | 4.7 | 45.4 | 24.9 |
| HG008T INDEL | GRCh38-only | STR | 195 (11.1 %) | 30.8 | 7.7 | 52.8 | 8.7 |
| HG008T INDEL | GRCh38-only | other | 261 (14.8 %) | 23.8 | 17.2 | 18.8 | 40.2 |
| HG008T INDEL | HPRC v1.1 d9 | HP>=7 | 1305 (74.1 %) | 52.9 | 4.4 | 19.5 | 23.2 |
| HG008T INDEL | HPRC v1.1 d9 | STR | 195 (11.1 %) | 85.6 | 3.6 | 8.7 | 2.1 |
| HG008T INDEL | HPRC v1.1 d9 | other | 261 (14.8 %) | 49.8 | 3.8 | 1.5 | 44.8 |
| HG008T INDEL | HPRC v2.1 d46 | HP>=7 | 1305 (74.1 %) | 52.7 | 5.4 | 18.5 | 23.4 |
| HG008T INDEL | HPRC v2.1 d46 | STR | 195 (11.1 %) | 86.2 | 3.6 | 9.2 | 1.0 |
| HG008T INDEL | HPRC v2.1 d46 | other | 261 (14.8 %) | 47.1 | 2.7 | 2.7 | 47.5 |

## Interpretation (R1.6 / R2.3)

The pooled numbers below are sums over the three platforms (sections 10 and 10a-e).

1. **SNVs, HG008T: about five in six truth SNVs are called on every arm.**
   - Called: 81.8 % (GRCh38-only), 82.8 % (d9), 84.9 % (d46).
   - Losses by stage:
     - no tensor: 2.1 / 5.2 / 4.6 %;
     - below the threshold: 10.9 / 9.0 / 6.8 %;
     - PoN: 5.2 / 3.0 / 3.7 %.
   - Of the truth the PoN does not tag, 89.5 / 90.6 / 92.9 % are called.
   - The graph arms' extra stage-1 loss is mainly absorbed and no-candidate alleles: 35 (d9) and 30 (d46) absorbed of
     2,091. Most of these are PoN-tagged anyway: 70 of 108 stage-1 records on d9, 65 of 97 on d46.
2. **SNVs, COLO829T: the cap is the build AF floor.**
   - No tensor: 17.3 % (GRCh38-only) and 18.0 % (d46).
     - Below the build AF floor 0.08: 807 (11.0 %) and 838 (11.4 %).
     - Too few supporting reads: 253 (3.5 %) and 308 (4.2 %).
   - By the truth VCF's own VAF, 16.7 % of the truth (1,227 records) lies below 0.08. Of these, 92.7 % / 93.3 % have
     no tensor.
   - The second loss is the threshold: 22.8 % on GRCh38-only against 9.9 % on d46. It sits in VAF > 0.4 (30.9 %
     against 12.4 %, by truth-VCF VAF).
   - Called: 58.3 % / 70.4 %. Of the untagged truth: 60.4 % / 73.0 %.
3. **INDELs: about two thirds of the truth is a population allele.**
   - The PoN tags the truth records themselves: HG008T 1,137 of 1,761 (64.6 %), COLO829T 243 of 393 (61.8 %). No caller
     that uses this PoN can report them.
   - Where they are lost depends on the arm:
     - **GRCh38-only**: they are called and then removed, so stage 3 holds 42.3 % (HG008T) and 39.7 % (COLO829T).
     - **Graph arms**: they never become tensors. Stage 1 holds 56.0 % (d9) and 55.6 % (d46) of HG008T and 52.9 % of
       COLO829T d46. On HG008T d9, 815 of the 987 stage-1 records are PoN-tagged. Absorbed alleles (ALT spelled by the
       graph) are 15.6 % (d9) and 14.1 % (d46) of the HG008T truth.
   - The net called share is the same on every arm: HG008T 25.4 / 24.1 / 24.5 %, COLO829T 26.7 / 24.4 %.
   - Of the untagged INDEL truth:
     - called: HG008T 71.6 / 67.9 / 69.1 %, COLO829T 70.0 / 64.0 %;
     - no tensor: HG008T 18.8 % on GRCh38-only against 27.6 % (d9) and 27.7 % (d46); COLO829T 6.0 % against 28.7 %
       (d46). This difference is the graph-specific candidate loss after the PoN.
4. **INDEL context.**
   - STR INDELs are almost never called: HG008T 8.7 / 2.1 / 1.0 %, COLO829T 7.4 / 3.7 %.
   - HP>=7 INDELs (74.1 % of the HG008T INDEL truth): 24.9 / 23.2 / 23.4 % called.
   - Other INDELs: 40.2 / 44.8 / 47.5 % called.
   - Insertions are called less often than deletions on the graph arms: HG008T INS 21.7 % against DEL 28.0 % (d9) and
     29.1 % (d46).
5. **Low VAF.**
   - HG008T SNVs at VAF 0.2-0.4: 77.1 / 72.5 / 74.3 % called, against 84.2 / 86.2 / 87.8 % at VAF > 0.4.
   - HG008T INDELs at VAF 0.08-0.2: 8.5 / 3.7 / 3.3 % called.
   - Below 0.08, every in-sample run calls nothing (section 10b), with one exception: on COLO829T, 6.1 % / 5.6 % of the
     truth-VCF-VAF < 0.08 SNVs are called (section 10c). Their own reads give them a tensor AF of at least 0.08.

## Caveats

1. **VAF source mix.** The VAF bin of a record without its own tensor comes from another arm, another platform or the
   truth VCF (section 7), so for HG008T the bins are not independent of stage 1. Example: on the graph arms, the
   HG008T INDEL VAF 0.08-0.2 bin is 71.0 % (d9) and 76.8 % (d46) no tensor. These are mostly alleles without a tensor of
   their own (absorbed, no candidate, min_variants), binned by the AF of another arm's or another platform's tensor
   (column `vaf_source` of `per_truth.tsv.gz`).

   For COLO829T, the truth-VCF bin (section 10c) avoids this dependence. However, its VAF_Ill / VAF_PB come from other
   sequencing data, so 8 / 17 COLO829T INDELs fall in the 0-0.08 bin by the build filter (10b) while the truth VCF puts
   none below 0.2 (10c).
2. **Threshold.** INDELs are evaluated at the standard threshold (validation recall 0.9), not at the chr1 curve's best
   F1, which is the project's INDEL comparison rule. Stage 2 of the INDEL rows therefore depends on that threshold.
   The `vcf_chr1_valbest` call sets exist and could be added to `call_sets()`.
3. **d9 recipe.** The HG008T PacBio / ONT SNV d9 checkpoints were trained with another recipe (1 GPU; E3 caveat 1), so
   their d9 rows are not recipe-matched to d46. COLO829T has no single-set d9 call set; the joint-model rows are
   supplementary.
4. **Absorbed.** `absorbed` uses E5b's spelling search: an exact path through the window subgraph. A no-candidate truth
   that the graph absorbs in another form (a nearby germline background, a different representation) is counted as
   `no_candidate`.
5. **Partial-match tensors.** 254 (d9) and 337 (d46) of the pooled HG008T INDEL stage-1 records have a PASS
   partial-match tensor, i.e. the model called the event as another allele (section 10e). rtg counts the truth as FN
   and that call as FP. The funnel keeps them in stage 1, because no tensor carries the exact allele.
6. **Pooling.** Pooled rows sum per-run records, so each truth record counts once per platform. They are not unique
   truths.

## Files

- `e16_recall_funnel.py`, `collect.sh`: this folder.
- `tables.md`:
  - sections 0-11;
  - per call set: funnel, stage-1 / stage-2 reasons, PoN tags, kind, VAF, VAF sources, COLO829T truth-VCF VAF,
    repeat;
  - pooled: section 10;
  - invariants: section 11.
- `summary.json`: the same counts per call set (funnel, reasons, splits, checks against `report.json`), pooled per
  sample x kind x arm, and the invariants.
- `/scratch/jshen/data/reviewer_experiments_20261013/E1.6_recall_funnel/per_truth.tsv.gz`: one row per call set x truth
  record. Columns: stage, reason, recall status, exact / partial tensors, max p_somatic, threshold, PoN tags, VAF and
  its source, truth-VCF VAF, repeat class.
- Log: `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E1.6/logs/login_collect_20261015.log`.
