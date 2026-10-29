# E9k: COLO829T graph-coordinate SNV calls against the COLO829BL-DSA somatic SNVs

Question (user 2026-10-29): which COLO829T off-reference (graph-coordinate, "unplaced") PASS SNV predictions are in the
somatic SNV set called on the COLO829BL donor-specific assembly (DSA),
`/wanglab/jshen/scratch_jshen/data/COLO829T/somatic_snvs.vcf.gz` (DeepVariant 1.6.0 on COLO829T PassageB Fiber-seq reads
aligned to the DSA, per its header, so not independent of the Fiber-seq data; 44,810 PASS SNVs; contigs
`haplotype1-*` = SMaHT verkko 2.1 hapX, `haplotype2-*` = hapY, `/scratch/qfu/COLO829BL_DSA/`).

**Status (2026-10-29): parts 1-4 done; used in R3.3 point 5 (rescue-off calls, part 3).**

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

## Caveats

- Rescue-on calls: ONT runs 3,000 of 3,421 / 4,569 unplaced PASS calls (E9c's sample); the rescue-off calls (the
  default) are complete (part 3).
- The GRCh38-truth comparison uses E9j's event haplotype only; a truth SNV placed on hapY by E9j but called on hapX in the
  DSA VCF is missed, so "4 of 31" is a lower bound.
- The DSA VCF is used as given (no FILTER other than PASS reported; all matched records are PASS).
