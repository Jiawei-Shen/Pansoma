# E2.4: call-level recall per truth VAF bin (low-VAF detection)

> **2026-10-16, INDEL PoN rule changed** (user decision: gnomAD / CoLoRSdb exact allele, AF >= 0.01 in homopolymers >= 8 bp, else AF >= 1e-4; was >= 7 bp and 1e-3). The 44 local chr1 INDEL call sets were re-scored in place (`tmp/reviewer_cpu_20261013/pon_hp8_20261016/rescore.sh`, array 418042; the old evaluations are kept as `eval_bed.pon_hp7_1e-3` / `eval_nobed.pon_hp7_1e-3`) and this experiment's collector was rerun (job 418071), so `tables.md` / `summary.json` / TSVs now use the new rule. INDEL numbers quoted in the prose below may still be those of the earlier rule (changes are small: e.g. HG008T Illumina v1.1 PASS F1 0.241 -> 0.251, best 0.265 -> 0.270; HG008T chr1 truth INDELs tagged 413 -> 418 of 635). SNV results are unchanged.


Reviewer points: **R2.4** (low-VAF detection); also R3.10 (low-VAF conclusions). Shared code:
`../call_level_common/calls_lib.py`.

**Status (2026-10-15 18:04, job 406571): every call set that exists locally is processed (53).**
- In-sample (45): SNV linear / d9 / d46 for HG008T Illumina / PacBio / ONT and linear / d46 for COLO829T Illumina /
  fiberseq / ONT (15); INDEL the same arms at the standard threshold and at `vcf_chr1_valbest` (30).
- COLO829T fiber-seq joint d9 model (8): SNV and INDEL x P100 / P75 / P50 / P25. Not a matched arm.
- Not included by design: cross-sample (`xs__`) and `__norescue` call sets.
- Missing: COLO829T d9 single-set runs (not local), comparators (not re-scored), seed repeats (GPU); see
  "Still missing" after the interpretation.

## Question

How does call-level recall (and precision) depend on the variant allele fraction? What share of the truth lies below
the tensor build AF floor of 0.08, where no tensor is built at all?

## Data

- **Call sets**: `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1*`, written by
  `machine_learning/pansoma_net_v2/jobs/vcf.sh`.
  - Truth in the BED: rtg `tp-baseline` + `fn`.
    - `pon_calls`: PASS at the run's threshold after the GRCh38 PoN (main result).
    - `raw_calls`: the same PASS calls before the PoN.
  - Included: in-sample runs (`main`, `vcf_chr1_valbest`) and the COLO829T fiberseq joint model per purity
    (`COLO829_fiberseq_{SNV_joint_7fc7396_gpu4,INDEL_joint_7fc7396}/vcf_chr1_<SET>`, when present).
- **Truth VAF**, in order of preference; the source of every record is counted in the tables:
  1. **Own tensor**: the tumor read AF of the truth allele in the call set's own tensor set.
     - It is the `af` of the matched allele in `<set>/tensors/<KIND>/chr1_variant_summary.ndjson`, for the tensor that
       `chr1_labels.ndjson` links to the truth record (the representative tensor first).
     - Exact allele matches only; partial matches are not used.
  2. **Other-arm tensor**: the same AF in the tensor set of another reference arm built from the same reads (linear,
     d9 or d46 of the same sample, platform and purity). A truth that the graph absorbed, for example, has no graph
     tensor but a linear one.
  3. **Truth VCF**: COLO829T `VAF_Ill` for Illumina sets, `VAF_PB` for fiberseq and ONT sets, multiplied by the purity
     for the P75 / P50 / P25 mixtures (an expected VAF).
  4. **Unknown**: HG008T's truth VCF has no read VAF (its `AD` comes from assembly haplotypes).

  The AF lookups are cached in `DATA/E2.4/truth_af/<set>.<KIND>.chr1.tsv.gz`.
- **Bins**: [0, 0.05), [0.05, 0.1), [0.1, 0.2), [0.2, 0.3), [0.3, 0.5), [0.5, 1], plus `unknown`.
- **Build AF floor**: 0.08 for SNV and INDEL (`pansoma_v2_tensors/README.txt`).

`DATA` = `/scratch/jshen/data/reviewer_experiments_20261013/call_level`.

## Method

1. Each truth record in the BED gets a VAF (above) and a bin.
2. Recall per bin = TP-baseline / truth of the bin, after the PoN and before it.
3. Precision per bin of the call's own read AF (FORMAT `AF` of the PASS call): TP / (TP + FP) after the PoN, in the
   BED.
4. Floor share per tensor set: the truth in the BED split by source, with the counts below 0.08.
5. Calibration of the truth-VCF fallback (COLO829T): Pearson r and the median ratio between the own-tensor AF and the
   truth-VCF VAF, over the truth records that have both.

## Scripts and jobs

| file | role |
|---|---|
| `vaf_bins.py` | truth AF caches, per-bin recall and precision, floor shares, calibration; `tables.md`, `summary.json` |
| `collect.sh` | runs `vaf_bins.py` |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E2.4_vaf_bins/job.sh` | Slurm wrapper (`-p general`, 4 CPUs, 12 GB); logs in `logs/` |

Jobs: 405328 (first run, 320 s; it builds the AF caches by scanning the chr1 `variant_summary` files), 405437,
405461 (snapshot 13:55) and 406571 (all 53 call sets, 18:04; 397 s with the new sets' AF caches, 79 MB max RSS).
`DATA` = `/scratch/jshen/data/reviewer_experiments_20261013/call_level`.

## Results (all call sets, 2026-10-15 18:00; 53 call sets)

chr1, in the BED; `tables.md` of 18:04. Arms: linear = GRCh38-only, d9 = HPRC v1.1 d9, d46 = HPRC v2.1 d46. The HG008T
PacBio / ONT SNV d9 checkpoints were trained with a different recipe (1 GPU; E3 caveat 1), so their d9 rows are not
recipe-matched to d46. INDEL rows are at the run's standard threshold (validation recall 0.9); the `[valbest]` rows
are in `tables.md`. The tables do not give the chr1 curve's best F1 (the INDEL comparison rule), so the INDEL rows
describe VAF dependence, not an arm ranking.

### 1. chr1 recall per truth VAF bin, in-sample SNV, PASS after the PoN (`tables.md` section 1)

Cells: `R (TP/truth)`. The last column gives the VAF source of the call set's truth records.

| call set (chr1) | 0-0.05 | 0.05-0.1 | 0.1-0.2 | 0.2-0.3 | 0.3-0.5 | >0.5 | unknown | all | VAF source: own tensor / other-arm tensor / truth VCF / unknown |
|---|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV linear | 0.000 (0/288) | 0.363 (57/157) | 0.815 (123/151) | 0.883 (68/77) | 0.842 (223/265) | 0.633 (952/1505) | - | 0.582 (1423/2443) | 2036 / 1 / 406 / 0 |
| COLO829T Illumina SNV d46 | 0.000 (0/288) | 0.348 (55/158) | 0.799 (119/149) | 0.853 (64/75) | 0.808 (215/266) | 0.857 (1292/1507) | - | 0.714 (1745/2443) | 2003 / 34 / 406 / 0 |
| COLO829T fiberseq SNV linear | 0.000 (0/301) | 0.232 (32/138) | 0.918 (145/158) | 0.852 (75/88) | 0.872 (238/273) | 0.697 (1035/1485) | - | 0.624 (1525/2443) | 2016 / 4 / 423 / 0 |
| COLO829T fiberseq SNV d46 | 0.000 (0/301) | 0.230 (32/139) | 0.897 (140/156) | 0.818 (72/88) | 0.855 (230/269) | 0.842 (1254/1490) | - | 0.707 (1728/2443) | 2007 / 13 / 423 / 0 |
| COLO829T ONT SNV linear | 0.000 (0/284) | 0.382 (66/173) | 0.708 (126/178) | 0.864 (51/59) | 0.874 (341/390) | 0.545 (741/1359) | - | 0.542 (1325/2443) | 2009 / 5 / 429 / 0 |
| COLO829T ONT SNV d46 | 0.000 (0/284) | 0.330 (58/176) | 0.672 (119/177) | 0.763 (45/59) | 0.857 (335/391) | 0.836 (1133/1356) | - | 0.692 (1690/2443) | 2001 / 13 / 429 / 0 |
| HG008T Illumina SNV linear | - | - | 0.000 (0/4) | 0.200 (2/10) | 0.885 (207/234) | 0.774 (339/438) | 0.000 (0/11) | 0.786 (548/697) | 684 / 2 / 0 / 11 |
| HG008T Illumina SNV d9 | - | 0.000 (0/1) | 0.000 (0/2) | 0.222 (2/9) | 0.868 (204/235) | 0.850 (373/439) | 0.000 (0/11) | 0.831 (579/697) | 660 / 26 / 0 / 11 |
| HG008T Illumina SNV d46 | - | - | 0.000 (0/1) | 0.250 (2/8) | 0.888 (207/233) | 0.867 (385/444) | 0.000 (0/11) | 0.852 (594/697) | 666 / 20 / 0 / 11 |
| HG008T PacBio SNV linear | - | 0.000 (0/2) | 0.250 (2/8) | 0.571 (4/7) | 0.900 (216/240) | 0.823 (354/430) | 0.000 (0/10) | 0.826 (576/697) | 683 / 4 / 0 / 10 |
| HG008T PacBio SNV d9 | - | 0.000 (0/2) | 0.000 (0/10) | 0.375 (3/8) | 0.877 (206/235) | 0.852 (368/432) | 0.000 (0/10) | 0.828 (577/697) | 659 / 28 / 0 / 10 |
| HG008T PacBio SNV d46 | - | 0.000 (0/1) | 0.000 (0/6) | 0.286 (2/7) | 0.850 (204/240) | 0.868 (376/433) | 0.000 (0/10) | 0.835 (582/697) | 662 / 25 / 0 / 10 |
| HG008T ONT SNV linear | - | 0.000 (0/5) | 0.231 (3/13) | 0.615 (8/13) | 0.937 (237/253) | 0.843 (339/402) | 0.000 (0/11) | 0.842 (587/697) | 681 / 5 / 0 / 11 |
| HG008T ONT SNV d9 | - | 0.000 (0/4) | 0.077 (1/13) | 0.667 (8/12) | 0.892 (224/251) | 0.845 (343/406) | 0.000 (0/11) | 0.826 (576/697) | 664 / 22 / 0 / 11 |
| HG008T ONT SNV d46 | - | 0.000 (0/4) | 0.222 (2/9) | 0.667 (8/12) | 0.922 (237/257) | 0.871 (352/404) | 0.000 (0/11) | 0.859 (599/697) | 666 / 20 / 0 / 11 |

### 2. chr1 recall per truth VAF bin, in-sample INDEL, standard threshold, after the PoN (`tables.md` section 1)

| call set (chr1) | 0-0.05 | 0.05-0.1 | 0.1-0.2 | 0.2-0.3 | 0.3-0.5 | >0.5 | unknown | all | VAF source: own tensor / other-arm tensor / truth VCF / unknown |
|---|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina INDEL linear | - | - | 0.250 (2/8) | 0.391 (9/23) | 0.292 (14/48) | 0.077 (4/52) | - | 0.221 (29/131) | 113 / 0 / 18 / 0 |
| COLO829T Illumina INDEL d46 | - | - | 0.222 (2/9) | 0.208 (5/24) | 0.261 (12/46) | 0.250 (13/52) | - | 0.244 (32/131) | 61 / 52 / 18 / 0 |
| COLO829T fiberseq INDEL linear | - | - | 0.250 (1/4) | 0.400 (4/10) | 0.341 (15/44) | 0.260 (19/73) | - | 0.298 (39/131) | 125 / 0 / 6 / 0 |
| COLO829T fiberseq INDEL d46 | - | 0.000 (0/5) | 0.000 (0/8) | 0.455 (5/11) | 0.275 (11/40) | 0.269 (18/67) | - | 0.260 (34/131) | 72 / 53 / 6 / 0 |
| COLO829T ONT INDEL linear | - | 0.000 (0/4) | 0.172 (5/29) | 0.250 (7/28) | 0.390 (16/41) | 0.310 (9/29) | - | 0.282 (37/131) | 125 / 1 / 5 / 0 |
| COLO829T ONT INDEL d46 | - | 0.000 (0/5) | 0.088 (3/34) | 0.214 (6/28) | 0.361 (13/36) | 0.286 (8/28) | - | 0.229 (30/131) | 71 / 55 / 5 / 0 |
| HG008T Illumina INDEL linear | - | 0.273 (3/11) | 0.115 (7/61) | 0.169 (15/89) | 0.359 (70/195) | 0.381 (59/155) | 0.000 (0/76) | 0.262 (154/587) | 506 / 5 / 0 / 76 |
| HG008T Illumina INDEL d9 | - | 0.000 (0/12) | 0.000 (0/65) | 0.138 (12/87) | 0.335 (63/188) | 0.409 (65/159) | 0.000 (0/76) | 0.239 (140/587) | 239 / 272 / 0 / 76 |
| HG008T Illumina INDEL d46 | - | 0.000 (0/14) | 0.030 (2/67) | 0.101 (9/89) | 0.363 (66/182) | 0.403 (64/159) | 0.000 (0/76) | 0.240 (141/587) | 261 / 250 / 0 / 76 |
| HG008T PacBio INDEL linear | - | 0.000 (0/4) | 0.086 (3/35) | 0.187 (17/91) | 0.305 (78/256) | 0.304 (55/181) | 0.000 (0/20) | 0.261 (153/587) | 566 / 1 / 0 / 20 |
| HG008T PacBio INDEL d9 | - | 0.000 (0/7) | 0.014 (1/71) | 0.093 (8/86) | 0.276 (66/239) | 0.457 (75/164) | 0.000 (0/20) | 0.256 (150/587) | 293 / 274 / 0 / 20 |
| HG008T PacBio INDEL d46 | - | 0.000 (0/8) | 0.000 (0/59) | 0.085 (8/94) | 0.281 (68/242) | 0.457 (75/164) | 0.000 (0/20) | 0.257 (151/587) | 301 / 266 / 0 / 20 |
| HG008T ONT INDEL linear | - | 0.050 (1/20) | 0.081 (15/185) | 0.250 (31/124) | 0.422 (46/109) | 0.610 (47/77) | 0.000 (0/72) | 0.239 (140/587) | 500 / 15 / 0 / 72 |
| HG008T ONT INDEL d9 | - | 0.040 (1/25) | 0.072 (13/180) | 0.214 (27/126) | 0.396 (42/106) | 0.654 (51/78) | 0.000 (0/72) | 0.228 (134/587) | 304 / 211 / 0 / 72 |
| HG008T ONT INDEL d46 | - | 0.000 (0/24) | 0.061 (11/179) | 0.260 (33/127) | 0.420 (47/112) | 0.658 (48/73) | 0.000 (0/72) | 0.237 (139/587) | 306 / 209 / 0 / 72 |

### 3. chr1 recall per truth VAF bin, COLO829T fiber-seq joint d9 model per purity (`tables.md` section 3)

Trained on P100 / P75 / P50 / P25 together on the d9 graph: not a matched arm. For P75 / P50 / P25 the truth-VCF
fallback VAF is the pure-sample VAF times the purity, and no other-arm tensor exists.

| call set (chr1) | 0-0.05 | 0.05-0.1 | 0.1-0.2 | 0.2-0.3 | 0.3-0.5 | >0.5 | unknown | all | VAF source: own tensor / other-arm tensor / truth VCF / unknown |
|---|---|---|---|---|---|---|---|---|---|
| COLO829T fiberseq SNV d9 (joint model) | 0.000 (0/301) | 0.212 (29/137) | 0.885 (139/157) | 0.759 (66/87) | 0.842 (230/273) | 0.835 (1243/1488) | - | 0.699 (1707/2443) | 2006 / 14 / 423 / 0 |
| COLO829T fiberseq SNV d9 P75 (joint model) | 0.000 (0/383) | 0.289 (35/121) | 0.779 (109/140) | 0.797 (55/69) | 0.826 (370/448) | 0.775 (994/1282) | - | 0.640 (1563/2443) | 1940 / 0 / 503 / 0 |
| COLO829T fiberseq SNV d9 P50 (joint model) | 0.000 (0/440) | 0.270 (40/148) | 0.798 (99/124) | 0.918 (335/365) | 0.820 (710/866) | 0.668 (334/500) | - | 0.621 (1518/2443) | 1870 / 0 / 573 / 0 |
| COLO829T fiberseq SNV d9 P25 (joint model) | 0.000 (0/592) | 0.420 (50/119) | 0.899 (829/922) | 0.921 (493/535) | 0.898 (247/275) | - | - | 0.663 (1619/2443) | 1737 / 0 / 706 / 0 |
| COLO829T fiberseq INDEL d9 (joint model) | - | 0.500 (1/2) | 0.000 (0/6) | 0.462 (6/13) | 0.268 (11/41) | 0.261 (18/69) | - | 0.275 (36/131) | 67 / 58 / 6 / 0 |
| COLO829T fiberseq INDEL d9 P75 (joint model) | - | 0.333 (1/3) | 0.125 (1/8) | 0.250 (4/16) | 0.275 (14/51) | 0.245 (13/53) | - | 0.252 (33/131) | 65 / 0 / 66 / 0 |
| COLO829T fiberseq INDEL d9 P50 (joint model) | - | 0.000 (0/3) | 0.158 (3/19) | 0.250 (9/36) | 0.235 (16/68) | 0.600 (3/5) | - | 0.237 (31/131) | 61 / 0 / 70 / 0 |
| COLO829T fiberseq INDEL d9 P25 (joint model) | - | 0.222 (4/18) | 0.203 (13/64) | 0.184 (7/38) | 0.545 (6/11) | - | - | 0.229 (30/131) | 58 / 0 / 73 / 0 |

### 4. Before the PoN (`tables.md` sections 2 and 4)

- **SNV**: recall before the PoN is 0.000-0.056 higher per bin (bins with >= 20 truth records). Over all truth it is
  0.014-0.018 higher for COLO829T in-sample (joint model 0.017-0.020) and 0.029-0.056 for HG008T (linear
  0.050-0.056, graph 0.029-0.040). The VAF pattern is the same.
- **INDEL**: the PoN removes much more. Before the PoN, linear recall is 0.538-0.756 and graph recall 0.359-0.443;
  after it, 0.221-0.298 and 0.228-0.260 (sections 1 and 2). The linear arm loses 0.300-0.494, the graph arms
  0.121-0.183, the joint model 0.137-0.176. So the linear arm's extra INDEL recall before the PoN is truth that the
  PoN removes (population alleles, gnomAD / CoLoRSdb rule).

### 5. chr1 recall per bin among the truth with an own-set tensor, SNV (`tables.md` section 5)

The model's sensitivity once a tensor exists: cells `R (TP/truth with a tensor)`. The 0-0.05 bin has no tensor in
any set and is left out.

| call set (chr1) | 0.05-0.1 | 0.1-0.2 | 0.2-0.3 | 0.3-0.5 | >0.5 |
|---|---|---|---|---|---|
| COLO829T Illumina SNV linear | 0.934 (57/61) | 0.932 (123/132) | 0.907 (68/75) | 0.842 (223/265) | 0.633 (952/1503) |
| COLO829T Illumina SNV d46 | 0.932 (55/59) | 0.944 (119/126) | 0.889 (64/72) | 0.840 (215/256) | 0.867 (1292/1490) |
| COLO829T fiberseq SNV linear | 0.970 (32/33) | 0.960 (145/151) | 0.882 (75/85) | 0.881 (238/270) | 0.701 (1035/1477) |
| COLO829T fiberseq SNV d46 | 0.941 (32/34) | 0.966 (140/145) | 0.847 (72/85) | 0.861 (230/267) | 0.850 (1254/1476) |
| COLO829T ONT SNV linear | 0.971 (66/68) | 0.887 (126/142) | 0.911 (51/56) | 0.877 (341/389) | 0.547 (741/1354) |
| COLO829T ONT SNV d46 | 0.853 (58/68) | 0.850 (119/140) | 0.804 (45/56) | 0.861 (335/389) | 0.841 (1133/1348) |
| HG008T Illumina SNV linear | - | 0.000 (0/4) | 0.222 (2/9) | 0.885 (207/234) | 0.776 (339/437) |
| HG008T Illumina SNV d9 | 0.000 (0/1) | - | 0.500 (2/4) | 0.879 (204/232) | 0.882 (373/423) |
| HG008T Illumina SNV d46 | - | 0.000 (0/1) | 0.500 (2/4) | 0.892 (207/232) | 0.897 (385/429) |
| HG008T PacBio SNV linear | 0.000 (0/2) | 0.250 (2/8) | 0.571 (4/7) | 0.904 (216/239) | 0.829 (354/427) |
| HG008T PacBio SNV d9 | - | 0.000 (0/4) | 0.600 (3/5) | 0.896 (206/230) | 0.876 (368/420) |
| HG008T PacBio SNV d46 | - | 0.000 (0/2) | 0.500 (2/4) | 0.879 (204/232) | 0.887 (376/424) |
| HG008T ONT SNV linear | 0.000 (0/4) | 0.273 (3/11) | 0.615 (8/13) | 0.937 (237/253) | 0.848 (339/400) |
| HG008T ONT SNV d9 | 0.000 (0/1) | 0.143 (1/7) | 0.800 (8/10) | 0.900 (224/249) | 0.864 (343/397) |
| HG008T ONT SNV d46 | 0.000 (0/1) | 0.400 (2/5) | 0.727 (8/11) | 0.937 (237/253) | 0.889 (352/396) |
| COLO829T fiberseq SNV d9 (joint model) | 0.935 (29/31) | 0.946 (139/147) | 0.786 (66/84) | 0.852 (230/270) | 0.843 (1243/1474) |
| COLO829T fiberseq SNV d9 P75 (joint model) | 0.897 (35/39) | 0.901 (109/121) | 0.833 (55/66) | 0.835 (370/443) | 0.782 (994/1271) |
| COLO829T fiberseq SNV d9 P50 (joint model) | 0.909 (40/44) | 0.884 (99/112) | 0.928 (335/361) | 0.828 (710/857) | 0.673 (334/496) |
| COLO829T fiberseq SNV d9 P25 (joint model) | 0.962 (50/52) | 0.928 (829/893) | 0.954 (493/517) | 0.898 (247/275) | - |

INDEL, same measure (standard threshold): at VAF > 0.5 the HG008T graph arms recall 0.703-0.889 of the truth with a
tensor (Illumina 0.722 / 0.703, PacBio 0.750 / 0.714, ONT 0.879 / 0.889 for d9 / d46), the linear arm 0.304-0.610
(0.381, 0.304, 0.610). At VAF 0.1-0.2 every in-sample INDEL set is at 0.000-0.400 (joint P25 0.591, 13 / 22).

### 6. chr1 truth without an own tensor and below the tensor build AF floor 0.08 (`tables.md` section 6)

No set has an own-tensor AF below 0.08 (that column of `tables.md` is 0 everywhere and is left out).

| set (chr1) | truth | with tensor | no tensor | no tensor: other-arm tensor AF (< 0.08) | no tensor: truth-VCF VAF (< 0.08) | no tensor: VAF unknown | known below floor |
|---|---|---|---|---|---|---|---|
| COLO829T_Illumina_linear SNV | 2443 | 2036 | 407 (16.7 %) | 1 (0) | 406 (368) | 0 | 368 (15.1 %) |
| COLO829T_Illumina_d46 SNV | 2443 | 2003 | 440 (18.0 %) | 34 (0) | 406 (368) | 0 | 368 (15.1 %) |
| COLO829T_fiberseq_linear SNV | 2443 | 2016 | 427 (17.5 %) | 4 (0) | 423 (394) | 0 | 394 (16.1 %) |
| COLO829T_fiberseq_d46 SNV | 2443 | 2007 | 436 (17.8 %) | 13 (0) | 423 (394) | 0 | 394 (16.1 %) |
| COLO829T_ONT_linear SNV | 2443 | 2009 | 434 (17.8 %) | 5 (0) | 429 (370) | 0 | 370 (15.1 %) |
| COLO829T_ONT_d46 SNV | 2443 | 2001 | 442 (18.1 %) | 13 (0) | 429 (370) | 0 | 370 (15.1 %) |
| HG008T_Illumina_linear SNV | 697 | 684 | 13 (1.9 %) | 2 (0) | 0 (0) | 11 | 0 (0.0 %) |
| HG008T_Illumina SNV | 697 | 660 | 37 (5.3 %) | 26 (0) | 0 (0) | 11 | 0 (0.0 %) |
| HG008T_Illumina_d46 SNV | 697 | 666 | 31 (4.4 %) | 20 (0) | 0 (0) | 11 | 0 (0.0 %) |
| HG008T_PacBio_linear SNV | 697 | 683 | 14 (2.0 %) | 4 (0) | 0 (0) | 10 | 0 (0.0 %) |
| HG008T_PacBio SNV | 697 | 659 | 38 (5.5 %) | 28 (0) | 0 (0) | 10 | 0 (0.0 %) |
| HG008T_PacBio_d46 SNV | 697 | 662 | 35 (5.0 %) | 25 (0) | 0 (0) | 10 | 0 (0.0 %) |
| HG008T_ONT_linear SNV | 697 | 681 | 16 (2.3 %) | 5 (0) | 0 (0) | 11 | 0 (0.0 %) |
| HG008T_ONT SNV | 697 | 664 | 33 (4.7 %) | 22 (0) | 0 (0) | 11 | 0 (0.0 %) |
| HG008T_ONT_d46 SNV | 697 | 666 | 31 (4.4 %) | 20 (0) | 0 (0) | 11 | 0 (0.0 %) |
| COLO829T_Illumina_linear INDEL | 131 | 113 | 18 (13.7 %) | 0 (0) | 18 (0) | 0 | 0 (0.0 %) |
| COLO829T_Illumina_d46 INDEL | 131 | 61 | 70 (53.4 %) | 52 (0) | 18 (0) | 0 | 0 (0.0 %) |
| COLO829T_fiberseq_linear INDEL | 131 | 125 | 6 (4.6 %) | 0 (0) | 6 (0) | 0 | 0 (0.0 %) |
| COLO829T_fiberseq_d46 INDEL | 131 | 72 | 59 (45.0 %) | 53 (0) | 6 (0) | 0 | 0 (0.0 %) |
| COLO829T_ONT_linear INDEL | 131 | 125 | 6 (4.6 %) | 1 (0) | 5 (0) | 0 | 0 (0.0 %) |
| COLO829T_ONT_d46 INDEL | 131 | 71 | 60 (45.8 %) | 55 (0) | 5 (0) | 0 | 0 (0.0 %) |
| HG008T_Illumina_linear INDEL | 587 | 506 | 81 (13.8 %) | 5 (0) | 0 (0) | 76 | 0 (0.0 %) |
| HG008T_Illumina INDEL | 587 | 239 | 348 (59.3 %) | 272 (0) | 0 (0) | 76 | 0 (0.0 %) |
| HG008T_Illumina_d46 INDEL | 587 | 261 | 326 (55.5 %) | 250 (0) | 0 (0) | 76 | 0 (0.0 %) |
| HG008T_PacBio_linear INDEL | 587 | 566 | 21 (3.6 %) | 1 (0) | 0 (0) | 20 | 0 (0.0 %) |
| HG008T_PacBio INDEL | 587 | 293 | 294 (50.1 %) | 274 (0) | 0 (0) | 20 | 0 (0.0 %) |
| HG008T_PacBio_d46 INDEL | 587 | 301 | 286 (48.7 %) | 266 (0) | 0 (0) | 20 | 0 (0.0 %) |
| HG008T_ONT_linear INDEL | 587 | 500 | 87 (14.8 %) | 15 (0) | 0 (0) | 72 | 0 (0.0 %) |
| HG008T_ONT INDEL | 587 | 304 | 283 (48.2 %) | 211 (0) | 0 (0) | 72 | 0 (0.0 %) |
| HG008T_ONT_d46 INDEL | 587 | 306 | 281 (47.9 %) | 209 (0) | 0 (0) | 72 | 0 (0.0 %) |
| COLO829T_fiberseq SNV | 2443 | 2006 | 437 (17.9 %) | 14 (0) | 423 (394) | 0 | 394 (16.1 %) |
| COLO829T_fiberseq_P75 SNV | 2443 | 1940 | 503 (20.6 %) | 0 (0) | 503 (438) | 0 | 438 (17.9 %) |
| COLO829T_fiberseq_P50 SNV | 2443 | 1870 | 573 (23.5 %) | 0 (0) | 573 (538) | 0 | 538 (22.0 %) |
| COLO829T_fiberseq_P25 SNV | 2443 | 1737 | 706 (28.9 %) | 0 (0) | 706 (657) | 0 | 657 (26.9 %) |
| COLO829T_fiberseq INDEL | 131 | 67 | 64 (48.9 %) | 58 (0) | 6 (0) | 0 | 0 (0.0 %) |
| COLO829T_fiberseq_P75 INDEL | 131 | 65 | 66 (50.4 %) | 0 (0) | 66 (0) | 0 | 0 (0.0 %) |
| COLO829T_fiberseq_P50 INDEL | 131 | 61 | 70 (53.4 %) | 0 (0) | 70 (0) | 0 | 0 (0.0 %) |
| COLO829T_fiberseq_P25 INDEL | 131 | 58 | 73 (55.7 %) | 0 (0) | 73 (4) | 0 | 4 (3.1 %) |

### 7. Calibration of the truth-VCF fallback, COLO829T (`tables.md` section 7)

Over the truth records that have both an own-set tensor AF and a truth-VCF VAF; the `[valbest]` rows duplicate the
standard ones and are left out.

| call set (chr1) | VCF field | n | Pearson r | median tensor AF | median VCF VAF | median ratio tensor / VCF | VCF VAF < 0.08 with a tensor | their median tensor AF |
|---|---|---|---|---|---|---|---|---|
| COLO829T Illumina SNV linear | VAF_Ill | 2036 | 0.985 | 0.973 | 0.976 | 1.000 | 36 | 0.097 |
| COLO829T Illumina SNV d46 | VAF_Ill | 2003 | 0.993 | 0.974 | 0.976 | 1.000 | 28 | 0.090 |
| COLO829T fiberseq SNV linear | VAF_PB | 2016 | 0.987 | 0.961 | 0.977 | 0.992 | 15 | 0.087 |
| COLO829T fiberseq SNV d46 | VAF_PB | 2007 | 0.989 | 0.961 | 0.977 | 0.992 | 17 | 0.089 |
| COLO829T ONT SNV linear | VAF_PB | 2009 | 0.964 | 0.823 | 0.977 | 0.900 | 38 | 0.097 |
| COLO829T ONT SNV d46 | VAF_PB | 2001 | 0.962 | 0.815 | 0.978 | 0.892 | 37 | 0.097 |
| COLO829T Illumina INDEL linear | VAF_Ill | 113 | 0.833 | 0.432 | 0.481 | 0.904 | 0 | - |
| COLO829T Illumina INDEL d46 | VAF_Ill | 61 | 0.934 | 0.429 | 0.471 | 0.897 | 0 | - |
| COLO829T fiberseq INDEL linear | VAF_PB | 125 | 0.838 | 0.515 | 0.616 | 0.979 | 0 | - |
| COLO829T fiberseq INDEL d46 | VAF_PB | 72 | 0.522 | 0.467 | 0.576 | 0.960 | 0 | - |
| COLO829T ONT INDEL linear | VAF_PB | 125 | 0.604 | 0.315 | 0.629 | 0.561 | 0 | - |
| COLO829T ONT INDEL d46 | VAF_PB | 71 | 0.462 | 0.257 | 0.576 | 0.492 | 0 | - |
| COLO829T fiberseq SNV d9 (joint model) | VAF_PB | 2006 | 0.991 | 0.961 | 0.977 | 0.992 | 15 | 0.091 |
| COLO829T fiberseq SNV d9 P75 (joint model) | VAF_PB | 1940 | 0.935 | 0.632 | 0.734 | 0.988 | 17 | 0.099 |
| COLO829T fiberseq SNV d9 P50 (joint model) | VAF_PB | 1870 | 0.782 | 0.385 | 0.490 | 1.015 | 41 | 0.093 |
| COLO829T fiberseq SNV d9 P25 (joint model) | VAF_PB | 1737 | 0.511 | 0.188 | 0.246 | 1.032 | 38 | 0.102 |
| COLO829T fiberseq INDEL d9 (joint model) | VAF_PB | 67 | 0.628 | 0.474 | 0.576 | 0.960 | 0 | - |
| COLO829T fiberseq INDEL d9 P75 (joint model) | VAF_PB | 65 | 0.471 | 0.394 | 0.432 | 0.962 | 0 | - |
| COLO829T fiberseq INDEL d9 P50 (joint model) | VAF_PB | 61 | 0.470 | 0.303 | 0.283 | 1.074 | 0 | - |
| COLO829T fiberseq INDEL d9 P25 (joint model) | VAF_PB | 58 | 0.252 | 0.197 | 0.141 | 1.359 | 5 | 0.108 |

### 8. chr1 precision per bin of the call's own read AF, SNV, PASS after the PoN (`tables.md` section 8)

No call has AF < 0.05 (the 0-0.05 column is left out).

| call set (chr1) | 0.05-0.1 | 0.1-0.2 | 0.2-0.3 | 0.3-0.5 | >0.5 |
|---|---|---|---|---|---|
| COLO829T Illumina SNV linear | 0.326 (57/175) | 0.184 (123/670) | 0.136 (68/500) | 0.314 (223/710) | 0.599 (952/1589) |
| COLO829T Illumina SNV d46 | 0.932 (55/59) | 0.944 (119/126) | 0.901 (64/71) | 0.960 (215/224) | 0.843 (1292/1533) |
| COLO829T fiberseq SNV linear | 0.395 (32/81) | 0.487 (145/298) | 0.540 (75/139) | 0.824 (238/289) | 0.835 (1035/1239) |
| COLO829T fiberseq SNV d46 | 0.667 (32/48) | 0.838 (140/167) | 0.973 (72/74) | 0.970 (230/237) | 0.853 (1254/1470) |
| COLO829T ONT SNV linear | 0.116 (66/569) | 0.101 (126/1248) | 0.093 (51/548) | 0.550 (341/620) | 0.810 (741/915) |
| COLO829T ONT SNV d46 | 0.523 (58/111) | 0.598 (119/199) | 0.417 (45/108) | 0.903 (335/371) | 0.834 (1133/1358) |
| HG008T Illumina SNV linear | 0.000 (0/28) | 0.000 (0/160) | 0.017 (2/117) | 0.396 (207/523) | 0.419 (339/809) |
| HG008T Illumina SNV d9 | - | 0.000 (0/7) | 0.091 (2/22) | 0.488 (204/418) | 0.447 (373/834) |
| HG008T Illumina SNV d46 | 0.000 (0/2) | 0.000 (0/24) | 0.062 (2/32) | 0.473 (207/438) | 0.433 (385/890) |
| HG008T PacBio SNV linear | 0.000 (0/3) | 0.050 (2/40) | 0.103 (4/39) | 0.468 (216/462) | 0.354 (354/1001) |
| HG008T PacBio SNV d9 | - | - | 0.500 (3/6) | 0.552 (206/373) | 0.471 (368/781) |
| HG008T PacBio SNV d46 | - | 0.000 (0/1) | 0.286 (2/7) | 0.533 (204/383) | 0.447 (376/841) |
| HG008T ONT SNV linear | 0.000 (0/4) | 0.042 (3/72) | 0.113 (8/71) | 0.458 (237/517) | 0.410 (339/827) |
| HG008T ONT SNV d9 | - | 0.200 (1/5) | 0.258 (8/31) | 0.516 (224/434) | 0.448 (343/766) |
| HG008T ONT SNV d46 | 0.000 (0/1) | 0.250 (2/8) | 0.151 (8/53) | 0.480 (237/494) | 0.424 (352/831) |
| COLO829T fiberseq SNV d9 (joint model) | 0.829 (29/35) | 0.897 (139/155) | 0.943 (66/70) | 0.966 (230/238) | 0.863 (1243/1440) |
| COLO829T fiberseq SNV d9 P75 (joint model) | 0.265 (35/132) | 0.288 (109/378) | 0.573 (55/96) | 0.971 (370/381) | 0.887 (994/1121) |
| COLO829T fiberseq SNV d9 P50 (joint model) | 0.741 (40/54) | 0.332 (99/298) | 0.626 (335/535) | 0.893 (710/795) | 0.756 (334/442) |
| COLO829T fiberseq SNV d9 P25 (joint model) | 0.746 (50/67) | 0.918 (829/903) | 0.781 (493/631) | 0.608 (247/406) | 0.000 (0/52) |

INDEL at the standard threshold: precision is <= 0.047 in every call-AF bin below 0.2 for every in-sample and joint
set, and at most 0.354 in any bin (HG008T Illumina d46, AF > 0.5). ONT INDEL calls pile up at AF 0.1-0.3: 9,000-22,900
calls per COLO829T / HG008T ONT set, with precision 0.000-0.007.

## Interpretation (all call sets)

1. **COLO829T low-VAF recall is limited by the tensor build floor, not by the model.**
   - 15.1-16.1 % of COLO829T's chr1 SNV truth in the BED has a known VAF below 0.08 (truth-VCF VAF, no tensor), on
     every in-sample arm and platform; no call can be made for these.
   - The 0-0.05 bin (284-301 truth records in-sample) has recall 0 in every call set.
   - In the 0.05-0.1 bin only 33-68 of 138-176 truth records have an own-set tensor (in-sample). Recall among those
     is 0.853-0.971, recall over the bin 0.230-0.382.
   - Once a tensor exists, recall at VAF 0.05-0.2 is 0.850-0.971 on every COLO829T in-sample arm (linear and d46,
     Illumina / fiber-seq / ONT) and 0.884-0.962 for the joint model at every purity.
2. **Purity series (joint d9 model, not a matched arm): the main loss with purity is the truth below the floor.**
   - Known below 0.08: 16.1 / 17.9 / 22.0 / 26.9 % (P100 / P75 / P50 / P25); the 0-0.05 bin grows from 301 to 383 /
     440 / 592 truth records, all missed.
   - Recall over all SNV truth: 0.699 / 0.640 / 0.621 / 0.663. Among the truth with a tensor, recall at VAF
     0.05-0.3 stays 0.786-0.962 at every purity (P25: 0.962 / 0.928 / 0.954).
   - P25 calls with AF > 0.5 are all false (0 / 52), as expected when every truth VAF is scaled by 0.25.
3. **HG008T cannot test low VAF.**
   - Its chr1 SNV truth has 9-31 records below VAF 0.3 per call set, and none known below 0.08.
   - 1.9-2.3 % (linear) and 4.4-5.5 % (d9, d46) of its SNV truth has no own-set tensor. Another arm's tensor of the
     same reads gives most of them a VAF (2-28 per set; graph-absorbed truth has a linear tensor, for example); 10-11
     have no VAF.
4. **The graph arms' INDEL tensor gap is not a low-VAF effect.**
   - 45.0-59.3 % of the INDEL truth has no own-set tensor in the graph arms (d9, d46, joint P100), against 3.6-14.8 %
     in the linear arm.
   - The linear tensors of the same reads give these a VAF, and none of them is below 0.08 (only P25 has 4 INDEL
     truth records known below the floor, from the scaled truth-VCF VAF).
   - The linear arm's extra recall from these tensors does not survive the PoN (section 4): after the PoN the graph
     and linear arms have the same INDEL recall (0.228-0.260 vs 0.221-0.298), so presumably most of this truth is
     population alleles.
5. **The linear arm loses high-VAF SNV truth.**
   - VAF > 0.5, linear vs d46: COLO829T Illumina 0.633 vs 0.857, fiber-seq 0.697 vs 0.842, ONT 0.545 vs 0.836; HG008T
     Illumina 0.774 vs 0.867, PacBio 0.823 vs 0.868, ONT 0.843 vs 0.871 (d9: 0.850 / 0.852 / 0.845).
   - The COLO829T loss is not in candidate generation (both arms have a tensor for >= 98.9 % of this bin) nor in the PoN
     (before it: 0.651 / 0.716 / 0.560 vs 0.877 / 0.862 / 0.855). Among the truth with a tensor: linear 0.633 /
     0.701 / 0.547 vs d46 0.867 / 0.850 / 0.841.
   - At VAF 0.3-0.5 the arms are similar or linear is slightly higher (COLO829T 0.842 / 0.872 / 0.874 vs 0.808 /
     0.855 / 0.857; HG008T 0.885 / 0.900 / 0.937 vs 0.888 / 0.850 / 0.922).
6. **The linear arm's false calls are low-AF, on every platform.**
   - HG008T linear, call AF 0.05-0.3: precision 0.000-0.017 (Illumina, 305 calls), 0.000-0.103 (PacBio, 82 calls),
     0.000-0.113 (ONT, 147 calls).
   - COLO829T linear, call AF 0.05-0.3: 0.136-0.326 (Illumina, 1,345 calls), 0.395-0.540 (fiber-seq, 518), 0.093-0.116
     (ONT, 2,365); d46: 0.901-0.944 (256), 0.667-0.973 (289), 0.417-0.598 (418).
   - This agrees with E4b (median AF 0.286 / 0.195 / 0.153 of the COLO829T linear "other" FPs) and E10 (extreme
     regions).
7. **ONT low-AF calls are less precise, and ONT tensor AFs run low.**
   - COLO829T ONT d46, call AF 0.05-0.3: precision 0.417-0.598, against 0.901-0.944 (Illumina d46) and 0.667-0.973
     (fiber-seq d46).
   - Tensor AF / VAF_PB median ratio: 0.892-0.900 for ONT SNVs and 0.492-0.561 for ONT INDELs, against 0.960-0.992
     for fiber-seq.
8. **The truth-VCF fallback is calibrated for pure samples; for mixtures it is an unbiased but noisy expected value.**
   - SNV: r 0.985-0.993 (VAF_Ill, Illumina), 0.987-0.991 (VAF_PB, fiber-seq) and 0.962-0.964 (ONT).
   - Joint P75 / P50 / P25: r 0.935 / 0.782 / 0.511, median ratio 0.988 / 1.015 / 1.032.
   - INDEL: r 0.252-0.934 (n 58-125), so the INDEL fallback VAFs are rough.

### Still missing

- **COLO829T d9 single-set runs are not local**, so COLO829T has no in-sample d9 rows; the only COLO829T d9 call sets
  are the fiber-seq joint model's (section 3, not a matched arm).
- **Comparators** (ClairS-TO / DeepSomatic-TO) are not yet re-scored on these chr1 truths; E2.4 has no comparator rows.
- **Seed repeats** are pending GPU: one training seed per model.
- **INDEL at the chr1 curve's best F1**: E2.4 reports the standard threshold and valbest only.
- No purity series for HG008T or for COLO829T Illumina / ONT: those dilutions are not among the
  `reviewer_cpu_20261013` call sets. Cross-sample (`xs__`) call sets are not included in E2.4.

## Caveats

- chr1 only. HG008T has too few low-VAF truth records for a VAF curve, so the COLO829T purity mixtures (joint
  model, section 3) are the low-VAF test.
- Bins mix VAF sources. The fallback VAFs are not the reads the model saw, but they are calibrated against them
  (section 7). For mixtures, the truth-VCF fallback is the pure-sample VAF times the purity: an expected value.
- Partial tensor matches (labels `residual_partial_somatic_truth`) are not used for the VAF.
- Recall per bin counts a truth record found if rtg matched it after the PoN. PoN-tagged truth (population alleles)
  is lost in every bin; the before-PoN tables (section 4) separate this.



## Rerun

All four call-level analyses: `bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/call_level_common/collect_all.sh`.
This one only:
`sbatch -J E2.4_vaf_bins -o /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E2.4_vaf_bins/logs/%x-%j.out /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E2.4_vaf_bins/job.sh`
(or `bash collect.sh` on the login node; once the AF caches exist a run takes seconds). The run rewrites `tables.md`
and `summary.json`; this README's numbers are those of 18:04 (job 406571).
