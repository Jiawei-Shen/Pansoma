# E17c: the PoN rule re-selected on chr2-22 (R2.6, panels of normals)

> **Decision 2026-10-16 (user):** the leave-one-tumor-out rule, HP >= 8: AF >= 0.01, else AF >= 1e-4, is adopted as the INDEL PoN rule (`scripts/filter_panel_of_normals.py`, commit 71cb51a) and used from now on. In the tables below "current" still means the earlier rule (HP >= 7: 0.01, else 1e-3); the purity re-scoring of the new rule was not pursued (user: no more purity tests).


## Question

The population filters were chosen by looking at chr1, the test chromosome:

- SNV: gnomAD / CoLoRSdb allele AF >= 1e-4, dbSNP non-somatic, 1000G (2026-09-30).
- INDEL: gnomAD / CoLoRSdb exact allele, AF >= 0.01 in homopolymers >= 7 bp, else >= 0.001 (2026-10-10).

A reviewer can say the rule was tuned on the test set. This experiment asks three things:

1. Which rule does a model-free selection on the training chromosomes (chr2-22) pick?
2. Would chr1 have picked the same rule?
3. For INDEL, which rule does a selection on chr1 calls pick when the scored tumor is left out?

## Method

**Inputs.** The labelled tensors of the six HPRC v1.1 d9 sets (HG008T PacBio / ONT / Illumina; COLO829T Illumina /
fiberseq / ONT), and the same six on HPRC v2.1 d46 as a replication. Labels 1 (somatic), 2 (germline) and 0 (non) are
used. Also the somatic truth records in the benchmark BED (from the sets' labels.manifest.json), and the four PoN VCFs.

**What the PoN sees.** For each tensor, the allele a call of that tensor would carry in the GRCh38 VCF:

- the labels' `grch38` projection, left-aligned and padded with `pansoma_net_v2.linear_vcf.Genome`;
- tensors without a projection ("unplaced") cannot be tagged.

**PoN features.** Matching uses `scripts/filter_panel_of_normals.py` itself. Its `record_matches`, `homopolymer_run`
and `record_kind` are imported; the script is not modified. Each PoN is streamed once per chromosome, and its records
at the positions of interest are handed to `record_matches`. For each allele the scan stores:

- the highest AF floor at which gnomAD / CoLoRSdb match;
- dbSNP non-somatic, and 1000G;
- the homopolymer run of the script;
- the pipeline's `RepeatContext` class (indexed_gam_pipeline_v4, imported).

**Checks** (`tables.md` §0):

- **Tabix sample:** 13,200 random alleles re-queried on the tabix files: 0 differences.
- **Script on chr22:** the unmodified script run on 3,973 chr22 tensor alleles: 0 per-PoN differences from the
  current rule computed from the features.
- **Pipeline truth tags:** the chr1 truth records the pipeline's vcfeval tags are reproduced exactly:
  - HG008T SNV 702 / 61; COLO829T SNV 2,443 / 84;
  - HG008T INDEL 635 / 413; COLO829T INDEL 131 / 81.

**Rule grid** (`e17c_tables.py`):

- **SNV:** AF floor 1e-5 / 1e-4 / 1e-3 / 0.01, each with and without dbSNP and with and without 1000G (16 rules).
- **INDEL:** one floor for every INDEL (1e-3, 0.01), or a homopolymer split (11 rules). The splits use HP >= 6 / 7 / 8
  with (other, HP) floors (1e-3, 0.01), (1e-3, 0.05) and (1e-4, 0.01).

**Objective.** Defined in the `e17c_tables.py` docstring before any grid number was computed; the same code is then
run on chr1.

- s = the fraction of somatic tensors a rule keeps; n = the fraction of negative tensors (germline + non) it keeps.
- F1 proxy = 2s / (1 + s + n). This is the post-PoN F1 of a caller that, before the PoN, finds every somatic
  candidate and makes as many false calls as true ones (precision 0.5). Its calls are tagged like the tensors.
- Selected rule = the highest mean F1 proxy over the six d9 sets. Ties (within 1e-4) go to the rule losing fewer
  somatic tensors.
- Reported, not used for the selection:
  - the F1 proxy with germline tensors only as negatives;
  - the candidate-level F1 with real counts (a caller that calls every candidate);
  - the d46 sets.

**chr1 calls** (`e17c_calls_chr1.py`). Every finished chr1 INDEL call set is used (vcf.sh done), 44 in all:

- in-sample, at the run's threshold and at the validation-chosen one (valbest);
- the joint COLO829T fiberseq model at each purity;
- cross-sample.

Each call set is tagged under a rule and scored with the pipeline's own functions (`pansoma_net_v2.vcfeval.vcfeval`
and `describe`, imported). Scoring uses the run's `truth.INDEL.vcf.gz`, `regions.bed` and benchmark BED with the
pipeline's rtg options:

- PASS: the PASS records the rule does not tag;
- curve: the PASS + LowQual records it does not tag, for best F1.

Rules scored:

- all 11 INDEL rules, PASS only, on the in-sample and joint P100 call sets;
- the current rule and the chr2-22 tensor rule, PASS and curve, on all 44.

The re-scored current rule must reproduce each run's own report. SNV, whose selection is the current rule, gets only
that check, on its in-sample d9 call sets.

**Leave-one-tumor-out selection on calls.** Defined in the `e17c_calls_chr1.py` docstring before any grid result
existed.

- A tumor's call sets: in-sample, at the run's threshold, P100.
  - HG008T: Illumina / PacBio / ONT × linear / d9 / d46 (9 call sets).
  - COLO829T: Illumina / fiberseq / ONT × linear / d46, plus the joint fiberseq model on d9 (7 call sets).
- Selected on tumor A: the INDEL rule with the highest mean chr1 PASS F1 over A's call sets. Ties within 0.0005 go to
  the higher mean recall.
- That rule is then scored on tumor B's call sets and compared with the current rule, in both directions.
- Repeated per platform (Illumina; PacBio-type = HG008T PacBio with COLO829T fiberseq; ONT) and on the valbest call
  sets.

## Result

**Short answer.**

- **SNV: yes.** The chr2-22 selection is the current rule, AF >= 1e-4 + dbSNP + 1000G. chr1 picks the same floor and
  dbSNP (with and without 1000G tie there, within 1e-4 of the F1 proxy), so the chr1 PoN-after numbers do not change.
- **INDEL: no.** chr2-22 selects HP >= 8: AF >= 0.05, else AF >= 1e-3; the current rule
  (HP >= 7: 0.01, else 1e-3) ranks 5 of 11.
  - chr1 under the same objective picks the same rule as chr2-22. So it is the objective, not the test chromosome, that
    separates it from the current rule.
  - On the 44 chr1 INDEL call sets, the tensor rule gives a lower PASS F1 in 29 and a higher one in 12 (median
    -0.002). Best F1 is lower in 34 (median -0.010). Every cross-sample set is worse (median -0.051).
- **INDEL, leave one tumor out on chr1 calls: not the current rule either, but a close neighbour.**
  - At the run's threshold, both directions pick HP >= 8: AF >= 0.01, else AF >= 1e-4. This keeps the current
    homopolymer floor; the HP threshold is 8 instead of 7 and the other floor 1e-4 instead of 1e-3.
  - The held-out tumor gains mean chr1 PASS F1: COLO829T 0.0307 → 0.0328 (+0.0021) and HG008T 0.1087 → 0.1139
    (+0.0052).
  - No held-out call set gets worse: 6 better and 1 equal of 7 for COLO829T; 9 better of 9 for HG008T.
  - The current rule ranks 5 of 11 on HG008T and 6 of 11 on COLO829T.

### SNV (`tables.md` §1-2, §5, §7)

**chr2-22, six d9 sets.**

- AF >= 1e-4 + dbSNP + 1000G has the highest mean F1 proxy, 0.8293 (rank 1 of 16). The runners-up are the same rule
  without 1000G (0.8288) and AF >= 1e-3 + dbSNP + 1000G (0.8275).
- Means over the six sets for this rule:
  - truth lost: 10.4 % (HG008T), 3.4 % (COLO829T);
  - somatic tensors lost: 6.3 %;
  - removed: 96.2 % of the germline tensors, 39.3 % of the non tensors.
- Per set, it is the best rule on the three COLO829T sets. The three HG008T sets each prefer another rule:
  AF >= 1e-3 + 1000G (PacBio), AF >= 1e-3 + dbSNP + 1000G (ONT), AF >= 0.01 + 1000G (Illumina).

**chr1.**

- The current rule is again rank 1 (0.8240).
- It ties within 1e-4 with the same rule without 1000G, and the tie-break takes that one.

**Other checks.**

- d46 sets: the current rule is selected on chr2-22 and on chr1.
- The choice depends on the objective:
  - the germline-only proxy picks AF >= 1e-3 (current rank 7);
  - the candidate-level F1 picks AF >= 1e-5 + dbSNP + 1000G (current rank 3).
- chr1 calls: the rule is unchanged, so only the reproduction check was run. The 3 d9 in-sample SNV call sets
  re-score identically to their pipeline reports (`calls_chr1.md`).

### INDEL (`tables.md` §3-7)

**chr2-22, six d9 sets.**

- Selected: HP >= 8: AF >= 0.05, else AF >= 1e-3 (mean F1 proxy 0.6760). The current rule scores 0.5972 (rank 5 of 11).
- What decides it is the homopolymer floor, 0.05 vs 0.01. The HP threshold hardly matters:
  - at 0.05, HP >= 6 / 7 / 8 give 0.6749-0.6760;
  - at 0.01, they give 0.5970-0.5972.

Selected vs current, means over the six sets:

| | selected | current |
|---|---:|---:|
| truth lost, HG008T | 35.7 % | 60.9 % |
| truth lost, COLO829T | 31.2 % | 50.6 % |
| somatic tensors lost | 14.1 % | 32.9 % |
| germline removed | 77.4 % | 90.0 % |
| non removed | 15.7 % | 23.9 % |

- **Where the rules differ** (§6). Only the RepeatContext HP>=7 class changes:
  - truth lost: 29.9 vs 63.6 % (HG008T), 25.0 vs 56.6 % (COLO829T);
  - germline removed: 73.2 vs 90.4 %; non removed: 5.3 vs 20.0 %.
  - The script's `homopolymer_run` >= 7 and RepeatContext HP>=7 agree on the truth INDEL alleles: 6,938 both, 0 run
    only, 185 RepeatContext only.
- **Robustness.**
  - The selected rule wins 5 of the 6 sets; HP >= 7 with 0.05 wins the sixth (HG008T ONT).
  - The d46 sets and the germline-only proxy pick the same rule.
  - The candidate-level F1 picks a stricter one: HP >= 8: AF >= 0.01, else AF >= 1e-4 (current rank 5).

**chr1.** The same rule is selected; the current rule ranks 9 of 11.

### chr1 INDEL calls under the chr2-22 tensor rule (`calls_chr1.md`)

There are 44 finished INDEL call sets. The re-scored current rule equals the pipeline report in 44 of 44.

| call sets | n | chr1 PASS F1 better / worse / equal | median change | chr1 best F1 better / worse | median change |
|---|---:|---|---:|---|---:|
| all | 44 | 12 / 29 / 3 | -0.002 | 10 / 34 | -0.010 |
| in-sample, run threshold | 15 | 3 / 12 / 0 | -0.002 | 3 / 12 | -0.010 |
| in-sample, valbest threshold | 15 | 6 / 6 / 3 | +0.000 | 3 / 12 | -0.010 |
| joint fiberseq, 4 purities | 4 | 3 / 1 / 0 | +0.000 | 4 / 0 | +0.006 |
| cross-sample | 10 | 0 / 10 / 0 | -0.051 | 0 / 10 | -0.068 |

- **Range.** PASS F1 changes run from -0.122 (HCC1395T with the COLO829T model) to +0.044.
- **Example, HG008T Illumina INDEL d9.** PASS F1 0.241 → 0.217 (recall 0.238 → 0.334, precision 0.244 → 0.161).
- **Against no PoN.** The current rule beats no PoN on chr1 PASS F1 in 40 of 44 call sets.

### Leave-one-tumor-out selection of the INDEL rule on chr1 calls (`heldout.md`)

**Run threshold, all platforms** (selection on 9 HG008T / 7 COLO829T call sets):

| selected on | rule | held-out tumor | held-out mean chr1 PASS F1 | current rule there | change | held-out call sets better / worse / equal |
|---|---|---|---:|---:|---:|---|
| HG008T | HP >= 8: 0.01, else 1e-4 | COLO829T | 0.0328 | 0.0307 | +0.0021 | 6 / 0 / 1 |
| COLO829T | HP >= 8: 0.01, else 1e-4 | HG008T | 0.1139 | 0.1087 | +0.0052 | 9 / 0 / 0 |

- **A tie on COLO829T.** There, AF >= 1e-3 for every INDEL is tied with the selected rule within 0.0005 (0.0330 vs
  0.0328); the recall tie-break takes HP >= 8: 0.01, else 1e-4.
- **The ranking orders the grid the same way on both tumors.**
  - Best: the three "HP floor 0.01, other floor 1e-4" rules.
  - Then the three "0.01, 1e-3" rules, which include the current one.
  - Then the three HP-floor-0.05 rules, which include the chr2-22 tensor rule.
  - HP >= 8 member of each group: 0.1139 / 0.1094 / 0.1039 on HG008T, 0.0328 / 0.0310 / 0.0272 on COLO829T.
  - AF >= 1e-3 for every INDEL is the exception: last on HG008T, first on COLO829T.

**Per platform, run threshold.**

| platform | selected on HG008T | change on COLO829T | selected on COLO829T | change on HG008T |
|---|---|---:|---|---:|
| Illumina (3 / 2 call sets) | HP >= 8: 0.01, else 1e-4 | +0.0049 | HP >= 8: 0.01, else 1e-4 | +0.0106 |
| PacBio-type (3 / 3) | HP >= 8: 0.01, else 1e-4 | +0.0016 | AF >= 1e-3 every INDEL | -0.0280 |
| ONT (3 / 2) | AF >= 1e-3 every INDEL | +0.0013 | AF >= 1e-3 every INDEL | +0.0039 |

The choices differ for long reads. In one direction this costs F1: COLO829T PacBio-type picks AF >= 1e-3, which drops
HG008T PacBio from 0.1029 to 0.0750 (worse in 3 of 3). COLO829T's long-read PASS F1 values are very small under the
current rule (fiberseq 0.007-0.040, ONT 0.002-0.004), so its long-read selections rest on few true calls.

**Validation-chosen threshold, all platforms** (9 / 6 call sets).

- Both directions pick HP >= 8: AF >= 0.05, else AF >= 1e-3, the chr2-22 tensor rule.
- COLO829T: 0.0318 → 0.0369 (+0.0051; better / worse / equal 2 / 1 / 3).
- HG008T: 0.1416 → 0.1471 (+0.0056; 4 / 5 / 0).
- The current rule ranks 7 and 8 of 11.
- Per platform, the HP-floor-0.05 rules are selected for Illumina and PacBio-type. ONT selects HP >= 8: 0.01, else 1e-4
  on HG008T (+0.0002 on COLO829T) and AF >= 1e-3 on COLO829T (-0.0250 on HG008T).

**In plain terms.**

- The current INDEL rule is not what a leave-one-tumor-out selection picks.
- At the run's threshold, the held-out selection picks a neighbour that keeps its 0.01 homopolymer floor and lowers the
  other floor to 1e-4. That neighbour raises the held-out tumor's mean chr1 PASS F1 by 0.0021 (COLO829T) and 0.0052
  (HG008T), with no held-out call set worse.
- So the current rule is not overfit to chr1 in a way that helps it. A rule chosen without the test tumor does at
  least as well on every held-out call set.
- At the validation-chosen threshold, the held-out selection moves to the 0.05 homopolymer floor. Its gain is similar
  in size (+0.005 to +0.006), but it is mixed per call set.
- Per platform, the long-read selections are unstable, driven by COLO829T's very low long-read F1.

### What this means for the reviewer's point

- **SNV:** the rule chosen on chr1 is the rule chr2-22 chooses.
- **INDEL, model-free:** the tensor objective on chr2-22 picks a more lenient rule than the current one.
  - On chr1 calls that rule is slightly worse at the run's threshold (median -0.002 PASS F1 in-sample, -0.051
    cross-sample).
  - chr1 tensors pick the same rule, so the gap comes from the call-level operating point, not from chr1.
- **INDEL, on calls:** with the scored tumor left out, the selection picks a neighbour of the current rule at the run's
  threshold. Both held-out tumors would gain a little (+0.002, +0.005) rather than lose. Reporting the current rule
  therefore does not flatter chr1.

## Caveats

- **The objective assumes precision 0.5 before the PoN.** The chr1 call sets have much lower precision without the PoN
  (`calls_chr1.md`): 0.102-0.108 for SNV, 0.000-0.333 for INDEL.
  - For F1, a PoN that removes a fraction r of the false calls and a fraction l of the true calls helps when
    r > l (1 + T / FP), with T the truth and FP the false calls. At precision 0.5 and recall 1 that is r > 2 l; with
    many false calls it approaches r > l.
  - So the proxy undervalues the false calls a PoN removes and leans to lenient rules. For SNV this does not change
    the choice.
  - For INDEL it explains both results: the proxy prefers the 0.05 homopolymer floor, and it scores the current rule
    (0.597) below no PoN (0.667). The chr1 calls disagree: the current rule beats no PoN in 40 of 44 sets.
- **The SNV choice depends on the objective** (germline-only proxy: AF >= 1e-3; candidate-level F1: AF >= 1e-5), but
  the primary objective picks the current rule on chr2-22, on chr1 and on d46.
- **Tensors, not calls.** The selection uses labelled candidates, not model calls. It cannot see what the model
  already rejects, for example common germline through the path-count channel.
- **Two truth sets.** Truth lost is per sample. HG008T and COLO829T are the only two samples, so the six sets are two
  truth sets seen through three platforms.
- **Leave-one-tumor-out uses chr1 calls.** It selects on one tumor's chr1 and scores on the other tumor's chr1. That
  removes the scored tumor, not the chromosome.
- **Two tumors, unequal F1.** Equal weights per call set let the Illumina call sets dominate both means. COLO829T
  long-read call sets have PASS F1 near 0, so its selections carry little information there.

## Files

| File | Content |
|---|---|
| `e17c_scan.py` | per chromosome: PoN features of every labelled tensor allele and every truth record → `/scratch/jshen/data/reviewer_experiments_20261013/E17c_pon_rule_offchr1/scan/chr*.json` (Slurm array, `tmp/reviewer_cpu_20261013/E17c/scan_array*.sh`; 2-5 min and ≤ 2.9 GB MaxRSS per chromosome) |
| `e17c_check_script.py` | the unmodified PoN script vs the features on chr22 → `.../check_script_chr22.json` |
| `e17c_tables.py` | objective, rule grid, selection → `tables.md`, `summary.json` |
| `e17c_calls_chr1.py` | `run`: chr1 INDEL call sets re-scored under the grid (Slurm `tmp/reviewer_cpu_20261013/E17c/calls_grid.sh`, 2 shards × 4 CPUs, 24-31 min, 367 rtg runs; rtg outputs in `.../calls_chr1_v2/`, first-version outputs reused from `.../calls_chr1/`). `report`: → `calls_chr1.md`, `heldout.md`, `calls_chr1.json` |
| `tables.md` | §0 checks and counts; §1-4 the grid on chr2-22 and chr1 (SNV, INDEL); §5 selections under the three objectives, d9 and d46, exchange rates; §6 INDEL by repeat context; §7 d46 replication |
| `calls_chr1.md` | chr1 PASS and best-curve F1 of every INDEL call set under the current and the chr2-22 tensor rule (and the SNV reproduction check) |
| `heldout.md` | leave-one-tumor-out INDEL rule selection on chr1 calls: every rule's mean PASS F1 per tumor, held-out scores, per platform, at the run and the validation-chosen threshold |

Rerun `calls_grid.sh` and then `python e17c_calls_chr1.py report` when more chr1 call sets finish: finished rtg
outputs are reused.
