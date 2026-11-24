# Why the comparators' INDEL F1 is far above the multi-cancer models' on HCC1395T (E12b), 2026-10-20

The user: "why are the other two models' INDELs so strong? It is not reasonable." E12b (`../reviewer_response_20261011/
experiments/E12b_multicancer_models/`) has, on HCC1395T chr1 INDEL, best F1 multicancer4 0.662 / 0.537 / 0.255
(Illumina / PacBio / ONT) and PASS F1 0.565 / 0.311 / 0.013, against DeepSomatic-TO 0.672 / 0.628 / 0.519 (PASS 0.605 /
0.612 / 0.471) and ClairS-TO PASS 0.641 / 0.566 / 0.569.

Data, copied from tequila by the user on 2026-10-20 into `$PANSOMA_RUNS/from_tequila_20261020/`: `runs/` the six
E12b INDEL runs (multicancer4 / multicancer6 x Illumina / PacBio / ONT; small files, `test_chr{1,20}_normfix`,
`vcf_*_normfix`; no checkpoints), `E8_score/` the comparators' scoring (`score_comparator.py` output: their calls,
rtg vcfeval per call set), `tensors/` the INDEL label manifests of the training lines and the HCC1395T chr1 / chr20
labels. `indel_fp.py` -> `tables.md` (all test sets, the comparators through E12b `summary.json`); `hcc1395_indel.py`
-> `hcc1395.md` (HCC1395T, every caller call by call).

## HCC1395T (`hcc1395.md`)

Same truth (SEQC2 v1.2.1, 133 INDELs on chr1, 42 on chr20), BED, PoN and rtg vcfeval for every caller; PASS after the
PoN.

| chr1 | multicancer4 TP / FP | DeepSomatic-TO | ClairS-TO | false calls in HP>=7 / STR (mc4 / DS / CT) | the others (mc4 / DS / CT) |
|---|---|---|---|---|---|
| Illumina | 98 / 116 | 102 / 102 | 92 / 62 | 58 / 37 / 11 | 58 / 65 / 51 |
| PacBio | 93 / 372 | 82 / 53 | 69 / 42 | 288 / 11 / 6 | 84 / 42 / 36 |
| ONT | 88 / 13,728 | 64 / 75 | 74 / 53 | 13,061 / 10 / 5 | 667 / 65 / 48 |

- Truth found: as many as the comparators or more, in every repeat class (ONT: HP4-6 12 of 18 against 4 / 12, HP>=7
  6 of 24 against 0 / 3, none 69 of 87 against 60 / 58).
- The gap is the false calls in homopolymers >= 7 bp and STRs on the long reads: 13,061 of 13,728 (ONT) and 288 of 372
  (PacBio). Of those, 13,024 / 13,037 (ONT) and 241 / 257 (PacBio) are not among DeepSomatic's non-RefCall records /
  ClairS-TO's records of any FILTER (DeepSomatic may have them as RefCall; the scorer drops those); their false calls
  in that class are 5-11. Median AF of
  the ONT false calls 0.18 against 0.44 for the true calls (`tables.md`).
- They are the class R13 takes out of the training labels of all four lines: R13 tensors are 46-50 % of the would-be
  negatives (label 0 + R13) on ONT, 63-71 % on PacBio, 75-96 % on Illumina (section 1); the test chromosome's R13
  tensors pass the threshold 12.0 % of the time against 0.1 % for label 0 (ONT, multicancer4).
- Outside HP>=7 / STR: Illumina level with the comparators (58 against 65 / 51), PacBio twice (84 against 42 / 36),
  ONT ten times (667 against 65 / 48; PASS threshold 0.042, validation precision at recall 0.9 only 0.094).
- chr20: the same pattern (ONT 3,742 of 3,955 false calls in HP>=7 / STR); DeepSomatic's chr20 is its training data.
- `indel_fp.repeat_class` agrees with the labels' own `repeat_class` on all 41,233 chr1 tensors compared.
- Real variants the truth lacks, or platform errors? (`crossplat.py` -> `crossplat.md`): a real variant is in every
  platform's reads. Of multicancer4's HP>=7 / STR false calls on chr1, the share present (same kind within 5 bp) in
  any other platform's call set (Pansoma multicancer4 PASS, DeepSomatic-TO or ClairS-TO records of any FILTER): ONT 97
  of 13,061 (0.7 %), PacBio 71 of 288 (24.7 %), Illumina 43 of 58 (74.1 %). So the ONT and PacBio ones are platform
  errors; the sites that may be real variants missing from the truth (or germline) are a few dozen on chr1, against
  about 1.4 million R13 tensors per line on ONT that are, by this test, almost all errors.
- Note on "no record": the scorer drops RefCall / NoCall records, so "not a record" for DeepSomatic-TO means no
  non-RefCall record (its raw output has 388,042 RefCall / NoCall records on ONT chr1). ClairS-TO writes no record when
  its INDEL network calls the site reference (`--print_ref_calls` off). So for neither caller does "no record" show that
  the site was never a candidate.
- Support of these false calls: on ONT 91.5 % have AF >= 0.1 and all have >= 3 alt reads (Pansoma's own pileup); PacBio
  95.1 % AF >= 0.1. A comparator-style candidate floor (DeepSomatic long-read tumor-only: each exact INDEL allele >= 0.10
  of the reads and >= 2 reads; ClairS-TO ONT: >= 0.10, >= 3 reads) would remove 8.5 % / 4.9 % of them and none of the
  true calls. (Their pileups differ: ClairS-TO ONT reads need MAPQ >= 20 and BQ >= 20; DeepSomatic drops insertions with
  mean BQ < 10.)

## How DeepSomatic-TO and ClairS-TO handle INDELs (source code and papers, checked 2026-10-20)

- DeepSomatic 1.9.0 / 1.10.0 (code in google/deepvariant; the 1.10 model JSONs carry the flags of the 1.9 script):
  candidates per exact allele, `vsc_min_fraction_indels` 0.07 (WGS_TUMOR_ONLY) / 0.10 (PACBIO / ONT_TUMOR_ONLY),
  `vsc_min_count_indels` 2, no read realignment on long reads, no homopolymer channel; every candidate is written
  (RefCall, GERMLINE = the network's 0/1, PON, LowQual, NoCall, PASS). Training: within the confident BED every candidate
  not in the truth is a negative (DeepVariant labeler); no homopolymer / STR exclusion in the preprint; tumor-only
  models trained on HCC1395 (SEQC2) chr2-20, chr21-22 tuning, chr1 held out (preprint p.4-5, p.11-12), so HCC1395 chr20
  is most likely its training data (not documented for the released PacBio / ONT tumor-only weights). Its own metrics
  page: tumor-only INDEL F1 on HCC1395 chr1 0.511 / 0.640 / 0.502 (WGS / PacBio / ONT).
- ClairS-TO 0.4.2 / 0.5.1: separate INDEL networks; candidates per exact allele with AF >= 0.10 (ONT) / 0.05 (ilmn,
  hifi_revio), >= 3 alt reads, depth >= 5, MAPQ >= 20, ONT BQ >= 20; no homopolymer handling at the candidate level
  (the early restriction to non-difficult regions was removed in 0.3.0); reference calls not written; LowSeqEntropy,
  VariantCluster, NonSomatic etc. change FILTER only. Training: synthetic tumors from GIAB pairs (plain platform names
  load these SS models), negatives = candidates that are not germline, no repeat exclusion; all models chr20 excluded;
  HCC1395 never used.
- Both train with every non-truth candidate as a negative. Pansoma's multi-cancer INDEL labels (rule (a), R13) are the
  only ones of the three that leave the homopolymer / STR non-truth candidates out.

## Decision (2026-10-20; the user: "you decide")

1. Drop rule (a) (R13, `--ignore-repeat-indels`) for the four multi-cancer lines on all three platforms and retrain
   the six multi-cancer INDEL models (multicancer4 / multicancer6 x Illumina / PacBio / ONT; same recipe, same seed).
   Reasons: (a) it labels as "not a negative" what the scoring counts as a false call, the only such mismatch in the
   pipeline; (b) by the cross-platform test the long-read candidates it removes are almost all platform errors (the
   possibly-real sites the truth lacks are a few dozen on chr1, against about 1.4 million R13 tensors per line on ONT);
   (c) HG008T and COLO829T, and both comparators, train with these candidates as negatives.
2. Scoring unchanged and the same for every caller. Rule (a) is not applied in the scoring: it would remove about
   13,000 of Pansoma's ONT false calls and about 10 of each comparator's, a rule shaped by our own training choice.
3. Candidate generation unchanged for now: a comparator-style AF floor would remove 5-9 % of these false calls. To be
   revisited if, after retraining, the ONT false calls outside HP>=7 / STR stay about ten times the comparators'.
4. Reporting: the retrained runs become the multi-cancer INDEL results; the rule-(a) runs stay documented (this
   README) as the first version. HCC1395T chr20 is held out for Pansoma and ClairS-TO but most likely training data
   for DeepSomatic-TO; say so wherever chr20 is compared. The response's rule-(a) paragraph (R2.2) and E12b README's
   "scored under rule (a)" are to be rewritten accordingly.

What to check after retraining: HP>=7 / STR false calls and the recall of the HP>=7 truth (6-9 of 24 on chr1 now), on
HCC1395T chr1 / chr20; HG008T / COLO829T as now.

## Steps (tequila; labels and tensors are there)

1. Back up the INDEL and SNV label files of the 12 sets `{HCC1395T,HCC1954T,H1437T,H2009T}_{Illumina,PacBio,ONT}_d46`
   (`tensors/{SNV,INDEL}/*labels*`, `labels.manifest.json`): `relabel.sh` keeps no backup and the current runs were
   trained on them.
2. `relabel.sh` without `IGNORE_REPEAT_INDELS`, with the same truth VCF / BED, germline truth and AF settings as the
   current `labels.manifest.json` (`truth`, `snv_min_af`, `indel_min_af`, `reference_path`).
3. Check: SNV labels byte-identical to the backup; INDEL labels differ only where `repeat_indel_without_truth` was -1.
4. Rebuild whatever index / chunk caches key on the labels; retrain the six INDEL models as `submit_multicancer.sh` did,
   under new names; predict and score the same test sets with the current code.
5. Optional, for the mechanism: DeepSomatic's raw `calls/deepsomatic/HCC1395T_ONT/output.vcf.gz` would show whether
   Pansoma's ONT repeat false-call sites are its RefCall candidates (seen and rejected) or not candidates at all.

## Findings

1. **Not recall: the same true calls, many more false ones.** HCC1395T chr1 (133 truth INDELs), PASS: multicancer4
   TP / FP 98 / 116 (Illumina), 93 / 372 (PacBio), 88 / 13,728 (ONT); DeepSomatic-TO about 102 / 102, 82 / 53,
   64 / 75; ClairS-TO about 92 / 62, 69 / 42, 74 / 53. On the curve at the comparators' TP the multi-cancer model has
   (multicancer4) 190 / 53 false calls at DeepSomatic's / ClairS-TO's TP on Illumina (102 / 62 theirs), 118 / 56 on
   PacBio (53 / 42), 719 / 1,844 on ONT (75 / 53): Illumina is level, PacBio somewhat and ONT far behind, so on long
   reads it is the ranking, not only the low PASS threshold (0.042 multicancer4, 0.0062 multicancer6 on ONT).
2. **The false calls are the tensors the training labels leave out.** The four training lines' labels use R13
   (`--ignore-repeat-indels`, `analysis/truth_labels.md`): an INDEL tensor in a GRCh38 homopolymer >= 7 bp or an STR
   without a truth allele is -1 (`repeat_indel_without_truth`), neither trained on nor validated, because SEQC2 / CASTLE
   truths miss many such INDELs. Scoring (rtg vcfeval, `pansoma_net_v2/vcfeval.py`, the same for the comparators) has
   no such rule: every call counts. On HCC1395T chr1:
   - ONT: 13,061 of the 13,728 PASS false calls are HP>=7 (11,905) or STR (1,156); 12,400 come from R13 tensors.
     Of the chr1 tensors, the R13 ones (126,114) are at or above the PASS threshold 12.0 % of the time, the label-0 ones
     (146,278) 0.1 % (multicancer6: 23.6 % against 0.2 %).
   - PacBio: 288 of 372 false calls HP>=7 / STR, 223 from R13 tensors; R13 tensors 79,700, 0.3 % above the threshold
     (label 0: 0.1 %; multicancer6 0.9 % / 0.1 %).
   - Illumina: only 985 R13 tensors; the false calls are mostly germline INDELs (69 of 116 tensors with germline truth).
   - The same on chr20 (ONT 3,558 of 3,955 false calls from R13 tensors).
   So the model never learnt to reject the long-read homopolymer / STR error INDELs of these lines, and validation
   (which excludes them as well) could not show it.
3. **The comparators' training.** DeepSomatic's tumor-only models were trained on HCC1395 with the same SEQC2 truth,
   chr2-20 training, chr21-22 tuning, chr1 held out (E8b D2, D3, D10): chr20 is one of its training chromosomes, so
   its HCC1395T chr20 numbers are not held-out. ClairS-TO never saw HCC1395 (E8b C3) and still has 42-62 false calls,
   so their precision on these lines is not only training exposure. That DeepSomatic learnt SEQC2's homopolymer
   convention (non-truth calls as negatives) is an inference, not checked.
4. **E12b README wording.** "HCC1395T INDELs are scored under rule (a)" holds for the labels only; the scoring counts
   the repeat INDEL calls without truth as false calls for every caller (point 2).
5. **HG008T / COLO829T.** The multi-cancer models' false calls there are label-0 tensors (`confident_no_truth_allele`,
   mostly HP>=7; no R13 on these sets), and the recall ceiling is 0.15-0.26 (HG008T: 397-404 of its 474 HP>=7 truth
   INDELs never called at any threshold); the comparators are low there too.

