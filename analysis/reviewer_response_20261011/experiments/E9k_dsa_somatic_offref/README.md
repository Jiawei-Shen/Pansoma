# E9k: COLO829T graph-coordinate SNV calls against the COLO829BL-DSA somatic SNVs

Question (user 2026-10-29): which COLO829T off-reference (graph-coordinate, "unplaced") PASS SNV predictions are in the
somatic SNV set called on the COLO829BL donor-specific assembly (DSA),
`/wanglab/jshen/scratch_jshen/data/COLO829T/somatic_snvs.vcf.gz` (DeepVariant 1.6.0 on COLO829T PassageB Fiber-seq reads
aligned to the DSA, per its header, so not independent of the Fiber-seq data; 44,810 PASS SNVs; contigs
`haplotype1-*` = SMaHT verkko 2.1 hapX, `haplotype2-*` = hapY, `/scratch/qfu/COLO829BL_DSA/`).

**Status (2026-10-29): parts 1-5b done; used in R3.3 point 5 (rescue-off calls, part 3).**

## Method

E9c already holds, for six COLO829T in-sample SNV runs on chr1 (v1.1 d9: Illumina epoch 6, Fiber-seq, ONT; v2.1 d46: the
three), every unplaced PASS call (ONT: the first 3,000 of 3,421 / 4,569) with its ALT read context, and E9d aligned every
context to DSA hapX and hapY (minimap2 -cx sr --cs, -N 50; 2026-10-23). `match_dsa_somatic.py` follows each alignment
(>= 90% of the context, NM <= 3, no indel within 5 bp of the site) to the site and reads the DSA base there:
- **DSA somatic**: the DSA base differs from the call's ALT and the VCF has that exact SNV (contig, position, REF = DSA
  base, ALT = the call's ALT in DSA orientation);
- **normal REF, not in VCF**: the DSA base differs from the ALT, no such SNV in the VCF;
- **normal carries ALT**: the DSA has the call's ALT (germline);
- **not placed**: no usable alignment.
Controls of the same runs: ctrl_tp (placed calls matching the GRCh38 SMaHT truth after the PoN), ctrl_fp (placed false
positives), ctrl_pon (placed calls removed by the PoN). Rescue-off PASS = PASS in the run's `__norescue` unplaced VCF
(E3b; the pipeline default since 91a8261).

## Results (`tables.md`, `per_call.tsv.gz`, `summary.json`)

- The method recovers the DSA truth: ctrl_tp 123-137 of 150 per run are DSA somatic SNVs (82-91%).
- Unplaced PASS calls in the DSA somatic set: v1.1 Illumina 15 / 693, Fiber-seq 7 / 1,483, ONT 12 / 3,000; v2.1 Illumina
  19 / 516, Fiber-seq 7 / 257, ONT 13 / 3,000. Rescue-off PASS calls (those E9c processed): 1-4 per run.
- 73 call records = **31 distinct DSA sites**, all on `haplotype1-0000012` (hapX chr1); many found by several runs
  (e.g. haplotype1-0000012:23310365 G>T by all six, five of them without the rescue). Only 4 of the 31 sit at the DSA
  position of a SMaHT GRCh38 truth SNV (E9j `events.tsv`, event haplotype only), so most are somatic SNVs of the DSA set
  that the GRCh38 benchmark does not score.
- Most unplaced PASS calls are not DSA somatic: the normal carries the ALT (germline, 112-983 per run) or the normal is
  REF and the VCF has no SNV there (90-1,957 per run; mostly ONT).
- ctrl_fp: 10-17 of 150 placed false positives (against the GRCh38 truth) are DSA somatic SNVs.

## Part 2: the DSA somatic SNVs on GRCh38 (`dsa_vs_grch38.py`, `dsa_vs_grch38.md` / `.json`, `dsa_snvs.tsv.gz`)

Each DSA somatic SNV is lifted with qfu's DSA -> hg38 chain/net (`/scratch/qfu/COLO829BL_DSA/COLO829BL_hg38.chain`,
parsed in Python) and compared with the GRCh38 base and the SMaHT GRCh38 truth (`COLO829T_somatic_snv_indel.vcf.gz`).

| class | DSA somatic SNVs | in the SMaHT GRCh38 truth |
|---|---:|---:|
| lifts to GRCh38, same reference base | 39,310 (87.7%) | 35,150 |
| lifts, but GRCh38 has another base (the SNV sits on a germline allele of the normal) | 148 (0.3%) | 0 |
| no GRCh38 position in the chain (DSA only) | 5,352 (11.9%) | 0 |

35,150 of the 42,035 autosomal SMaHT truth SNVs (83.6%) are DSA somatic SNVs by this lift.

The 31 DSA somatic sites that Pansoma calls in graph coordinates (part 1): 14 DSA only, 2 on a germline allele, 15 with
the same GRCh38 base (2 of them lift to chr1 random contigs, chr1_KI270713v1_random / chr1_KI270714v1_random), and only 5
of the 31 are in the SMaHT GRCh38 truth (part 1's E9j-based count of 4 missed one placed on the other haplotype).

Caveats: "DSA only" means outside every chain block; a block can also break at indels, structural variants or poorly
aligned sequence, so this is sequence the chain does not align to GRCh38, not necessarily sequence absent from GRCh38.

## Part 3: every unplaced PASS call without the rescue (the default since 91a8261)

`norescue_contexts.py` ran E9c's `e9c_run.py` on the six `<model>__norescue/vcf_chr1` re-calls (E3b) with E9c's DATA
redirected to `/scratch/jshen/data/reviewer_experiments_20261013/E9k_dsa_somatic_offref` (no controls; E9c's folder
untouched), and the 328 alt_read contexts were aligned to DSA hapX / hapY as E9d (minimap2 2.28 -cx sr --cs --secondary=yes
-N 50); jobs 460540 / 460541 (`tmp/reviewer_cpu_20261029/E9k/`). Group `norescue_pass` in `tables.md`:

| graph | Illumina | Fiber-seq | ONT |
|---|---|---|---|
| v1.1 (d9) | 3 of 38 DSA somatic | 2 of 35 | 8 of 152 |
| v2.1 (d46) | 4 of 18 | 2 of 78 | 4 of 98 |

The rest: normal REF and no DSA somatic SNV there 13-137 per run (mostly ONT), normal carries the ALT 1-27. With these,
the graph-coordinate calls match **35 distinct DSA somatic sites** (part 2 table updated): 17 DSA only, 2 on a germline
allele, 16 with the same GRCh38 base, 6 of the 35 in the SMaHT GRCh38 truth. Ten of the 35 sites are called without the
rescue, nine of them DSA only.

## Part 4: can the DSA-only SNVs get a tensor? (`dsa_only_tensors.py`, `dsa_only_tensors.md`, `dsa_only_labels.md`)

The 5,352 DSA-only SNVs and the 148 on a germline allele were placed on each graph as in E9j: the DSA window of the SNV's
haplotype (100-bp flanks) without (N) and with (T) the somatic base, giraffe (vg 1.65.0 d9, vg 1.77.0 d46), E9j's
`s3_classify.py` unchanged (A = the somatic edit on a GRCh38 node, B = on an off-reference node, C = the allele is a graph
branch, D = unmapped / MAPQ < 20 / same alignment); jobs 460547 / 460548 / 460555 (`tmp/reviewer_cpu_20261029/E9k/`). The key
of each A / B single-base edit (node:start:SNP:REF>ALT, E9j's rule) was looked up in the six COLO829T SNV tensor sets
(representative allele A1 exact), and the GRCh38 position of each A key was computed as the labels do (ReferencePath.linear).

DSA-only SNVs (5,352):

| | d9 | d46 |
|---|---:|---:|
| D: unmapped / MAPQ < 20 / same alignment (MAPQ < 20: 2,754 / 2,823) | 2,992 | 3,055 |
| C: the allele is a graph branch | 43 | 41 |
| A, B with a complex (not single-base) edit | 99 | 96 |
| A on a non-autosome (chrUn / random 550 / 516, chrX / chrY 30): no tensor, the sets cover chr1-22 | 580 | 546 |
| A on chr1-22, a SMaHT GRCh38 truth SNV (the chain misses the position) | 122 | 123 |
| A on chr1-22, not in the truth | 1,345 | 1,093 |
| B: off-reference node | 171 | 398 |

Tensors (exact A1) per set, Illumina / Fiber-seq / ONT:

| | d9 | d46 |
|---|---|---|
| A, SMaHT truth SNV (122 / 123) | 112 / 91 / 88, labelled somatic | 111 / 91 / 87, labelled somatic |
| A, not in the truth, label region (dip.bed; 177 / 120) | 23 / 16 / 13, labelled non (germline 2) | 23 / 16 / 14, labelled non (germline 1) |
| A, not in the truth, outside dip.bed (1,168 / 973) | 131 / 41 / 43, ignore | 119 / 60 / 44, ignore |
| B, off-reference node (171 / 398) | 126 / 132 / 131, ignore (5 somatic) | 139 / 155 / 142, ignore (4 somatic) |
| DSA-only SNVs with a tensor in at least one set | 428 | 467 |

Baseline: 89.6-90.5% of the 41,230 SMaHT truth SNVs on chr1-22 in dip.bed are the A1 of a somatic-labelled tensor. The
DSA-only SNVs that the chain misses but the GRCh38 truth has reach about the same rate (71-92%); the other A events on
GRCh38 nodes rarely have a tensor (4-13% per set; why not checked); B events have one in 74-77% (d9) and 35-39% (d46) of the
cases. The off-reference tensors are labelled ignore (no GRCh38 truth to match), so they are neither trained nor scored;
the DSA VCF could label them somatic (126-155 per set). Of the 148 SNVs on a germline allele, 48 (d9) / 42 (d46) are
A or B; their tensors are mostly ignore or germline (`dsa_only_labels.md`).

## Part 5: tensors on the HPRC graphs for all 44,810 DSA somatic SNVs (`dsa_all_tensors.py`, `dsa_all_tensors.md` / `.json`)

Part 4's steps on every SNV (windows of all 44,810, giraffe to d9 / d46, `s3_classify.py`; jobs 461133 / 461134 / 461135).
A set has the SNV when the representative allele A1 of one of its tensors is the SNV: by the candidate key of the window's
A / B edit, or by the SNV's GRCh38 lift (the labels' `grch38` position and ALT); the two agree for about 35,500 per set.

| DSA class | SNVs | per set (6 sets) | any v1.1 set | any v2.1 set | any of the six |
|---|---:|---|---:|---:|---:|
| lifts, same base, in the SMaHT GRCh38 truth | 35,150 | 99.3-99.6% | 99.7% | 99.8% | 35,091 (99.8%) |
| lifts, same base, not in the truth | 4,160 | 29.1-33.8% | 36.1% | 36.5% | 1,574 (37.8%) |
| lifts, on a germline allele | 148 | 18.9-26.4% | 27.7% | 23.6% | 42 (28.4%) |
| no GRCh38 position in the chain | 5,352 | 5.1-7.3% | 8.0% | 8.7% | 516 (9.6%) |
| all | 44,810 | 81.3-81.9% | 82.6% | 82.8% | 37,223 (83.1%) |

- The 59 truth SNVs without a tensor in any set: 35 are class C on both graphs (the somatic allele is a graph path, so no
  candidate; graph_absorbed_somatic_colo829t_20261005), 13 A, 10 D.
- The 4,160 lifted SNVs not in the SMaHT truth: their tensors are labelled non in 894-1,038 cases per set (somatic
  105-123, ignore 154-252, germline 8-11), i.e. trained and scored as negatives.
- DSA-only SNVs of class D (2,992 / 3,055) have no key and no lift, so this method cannot find their tensors (0 by
  construction, not shown to be absent); for the lifted ones the GRCh38 lookup does (e.g. 301-411 of the ~430 truth SNVs of
  class D have a tensor).
- "Not in the truth" is against the truth's chr1-22 SNVs (as part 2; the tensor sets are chr1-22): 1,664 of the 4,160
  lift to chrX / chrY (1,557 of them are in the SMaHT truth's chrX / chrY part) and 241 to other GRCh38 contigs, so 1,905
  cannot have a tensor in these sets (part 5b).

## Part 5b: why many have no tensor: ALT reads in the COLO829T GRCh38 BAMs (`dsa_all_tensors.py pileup`, `dsa_all_pileup.md` / `.json`)

Reads at each SNV's GRCh38 position (the lift; for the DSA-only SNVs of class A on d9, the GRCh38 position of the graph
node, part 4b) in the BAMs of the reads the tensor sets were built from (Fiber-seq B004, Illumina A001); job 461357
(461142 timed out: a fetch per site re-read the long reads, 757 GB in 2 h; Fiber-seq is now read once per chromosome).
"ALT supported" = ALT in >= 3 reads and >= 8% of reads at MAPQ >= 10 (the builder's min_mapq and SNV AF). Fiber-seq /
Illumina; tensor = A1 is the SNV in the v1.1 or v2.1 set of that platform.

| SNVs | n | ALT supported | tensor, ALT supported | tensor, all | why no tensor (main) |
|---|---:|---:|---:|---:|---|
| lifts, in the truth (random 2,000, control) | 2,000 | 99.7% / 99.8% | 100% / 100% | 99% / 100% | - |
| lifts, not in the truth, chrX / chrY | 1,664 | 96.7% / 95.3% | 0 | 0 | outside the chr1-22 tensor sets (1,557 are SMaHT truth on chrX / chrY) |
| lifts, not in the truth, other GRCh38 contigs | 241 | 68.9% / 58.9% | 4% / 5% | 4% / 5% | outside the tensor sets |
| lifts, not in the truth, chr1-22 | 2,255 | 58.0% / 47.6% | 86% / 92% | 57% / 55% | no ALT read 23.4% / 26.0% (376 in both); ALT only at MAPQ < 10 9.6% / 18.2%; graph class D 1,023 |
| germline allele, ALT = the GRCh38 base | 115 | (reference reads) | - | 10% / 16% | the mutation turns the germline allele back to the GRCh38 base: no variant on the GRCh38 path (d9 class C 48, D 42) |
| germline allele, another ALT | 33 | 66.7% / 60.6% | 91% / 100% | 67% / 67% | no ALT read 9 / 6 |
| DSA only, class A on d9, node on chr1-22 | 1,467 | 11.8% / 17.9% | 72% / 86% | 13% / 19% | ALT >= 3 reads but < 8%: 38.9% / 70.6% (below); Fiber-seq no ALT read at the node 30.0% |
| DSA only, class A on d9, node elsewhere | 580 | 10.0% / 11.6% | 2% / 1% | 1% / 1% | node on chrUn / random / chrX / chrY, outside the tensor sets |

- The DSA-only nodes on chr1-22 sit in collapsed repeat copies (part 4: 90% TRF): median depth 305 / 4,968 against 188 /
  121 at the truth controls (82% of the Illumina sites >= 2x), median ALT fraction 0.3% / 0.8%. The ALT is in the reads
  but diluted by the reads of the other copies, below the builder's 8%.
- So where the SNV is on chr1-22 and our reads support it, 86-92% of the lifted non-truth SNVs and 72-86% of the DSA-only
  ones have a tensor. The low rates of part 5 come from SNVs outside chr1-22 (1,905 lifts, 580 DSA-only nodes), DSA calls
  without ALT reads in these BAMs, ALT diluted in collapsed repeats, reversions to the GRCh38 base, and windows Giraffe
  does not place uniquely (class D: 2,992 DSA-only, 1,023 chr1-22 lifts).

## Caveats

- Rescue-on calls: ONT runs 3,000 of 3,421 / 4,569 unplaced PASS calls (E9c's sample); the rescue-off calls (the
  default) are complete (part 3).
- The GRCh38-truth comparison uses E9j's event haplotype only; a truth SNV placed on hapY by E9j but called on hapX in the
  DSA VCF is missed, so "4 of 31" is a lower bound.
- The DSA VCF is used as given (no FILTER other than PASS reported; all matched records are PASS).
