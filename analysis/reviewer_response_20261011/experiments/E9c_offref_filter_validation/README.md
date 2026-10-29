# E9c / E9d: graph-native population filter and assembly validation of unplaced (graph-coordinate) calls

Reviewer points: **R3.3** (validation of graph-coordinate variants), **R3.4** (metrics exclude graph-coordinate calls;
"about 2 %"), **R1.3** (ONT graph-coordinate calls: resolution or ambiguity). Design: `../../response_to_reviewers_detailed.md`,
sections R3.3 / R3.4 / R1.3 and rows E9c / E9d / E9e of the E-table (section 6).

> **2026-10-23: COLO829T assemblies changed.** Normal = the COLO829BL donor-specific assembly (SMaHT verkko 2.1 hapX / hapY,
> all contigs; before: their ragtag-hg38 primary + unlocalized parts); tumor = qfu's hifiasm 0.25.0 of the COLO829T Fiber-seq
> HiFi reads with the COLO829T ONT D003 reads (`asm_targets/`; not independent of the COLO829T Fiber-seq / ONT call sets).
> COLO829T now has the HG008T categories (tumor_only / unsupported instead of not_in_normal); asm 451718, collect 451719
> (`tmp/reviewer_cpu_20261023/E9i/`). `tables.md` / `summary.json` / `per_call.tsv.gz` are the new ones; COLO829T numbers in
> the text below that say "not in normal" are the old classification. Controls: COLO829T ctrl_tp tumor-only 130-137 / 150,
> ctrl_fp germline 94-134 / 150; unplaced PASS tumor-only 0.3-3.3 % (population first), 0.6-4.0 % by the assemblies alone.

**Status (2026-10-15, tables of collect job 405796, 15:27): every in-sample SNV run on d9 / d46 done.** 13 runs: the 9
RIS 7fc7396 in-sample SNV re-predictions on a pangenome graph (HG008T Illumina / PacBio / ONT on HPRC v1.1 d9 and on
HPRC v2.1 d46; COLO829T Illumina / fiber-seq / ONT on d46) and 4 of the 12 legacy d9 runs; the 6 gold-set calibrations
(d9 and d46), both assembly checks and the tables. The tables' inputs (`DATA/runs/*/calls.tsv`, newest 15:18, and the
gold sets) have not changed since. Not processed: the INDEL, `__norescue`, joint-model and cross-sample
re-predictions, and the 8 other legacy runs (see "Status" at the end). The GRCh38-only arm has no off-reference nodes,
so it has no unplaced calls and is outside this experiment.

## Question

A PASS call on a node without a GRCh38 position ("unplaced", `*.linear.unplaced.vcf.gz`) escapes both the GRCh38
panel of normals (PoN) and rtg vcfeval. What are these calls?

1. **Population alleles** that the GRCh38 PoN cannot see (E9c: a graph-native PoN).
2. **Germline** variants of the patient (in the HG008-N normal assembly / the COLO829BL assembly).
3. **Assembly-supported somatic** variants (ALT context in the HG008-T tumor assembly, absent from HG008-N).
4. **Unsupported** (in no assembly: likely sequencing or alignment errors), or not verifiable.

From these classes follows a precision over the complete output (R3.4), per platform (R1.3), together with
cross-platform reproducibility and the fate of the assembly-verified off-reference somatic gold set (R3.3).

## Data

- **Runs (local, "legacy")**: the local d9 runs of the current recipes, chr1 test, PASS at the checkpoint's
  validation recall-0.9 threshold, GRCh38 PoN of the run's time:
  - SNV `{HG008,COLO829T}_{Illumina,PacBio|fiberseq,ONT}_SNV_pow075_S_lr1e4_e20` (2026-09-30 / 10-01; earlier than
    the 2026-10-03 `*_SNV_h90` runs, which the RIS 7fc7396 SNV runs reproduce exactly for HG008T PacBio / ONT according
    to `analysis/benchmark_7fc7396_20261012/README.md`);
  - INDEL `{HG008,COLO829T}_*_INDEL_small_nofilt` (2026-10-10: small model, sqrt, no filter, lr 4e-4 (`args.json`);
    the 7fc7396 RIS recipe uses lr 1e-4, which the benchmark README reports as equal under the same PoN).
  - Processed: `HG008_Illumina_SNV`, `HG008_Illumina_INDEL`, `HG008_PacBio_SNV`, `COLO829T_ONT_SNV`.
- **Runs (RIS 7fc7396 re-predictions, `rcpu`)**: the RIS 7fc7396 checkpoints re-predicted by the main session into
  `pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1*` (CPU array 387323, tequila GPU arrays 387379 and 405185).
  `discover.py` lists every finished one on d9 / d46; `run_all.sh` processes them. Processed: the nine in-sample SNV
  runs, `HG008_{Illumina,PacBio,ONT}_SNV_d46_7fc7396`, `COLO829_{Illumina,fiberseq,ONT}_SNV_d46_7fc7396`,
  `HG008_Illumina_SNV_7fc7396_ddp4_b1024` and `HG008_{PacBio,ONT}_SNV_7fc7396` (d9; the PacBio / ONT d9 checkpoints
  were trained on 1 GPU, E3 caveat 1). COLO829T has no single-set d9 run on this cluster.
- **Graphs**: HPRC v1.1 d9 and v2.1 d46: node sequences and path counts (`graph_index.sqlite`), GRCh38 reference path
  (`grch38_path`), and an edge table built here from the GFA L lines (`build_edges.py`: d9 70,105,101 L lines,
  d46 103,555,446; CSR over oriented node handles in `DATA/edges_{d9,d46}`).
- **Population haplotypes**: `hprc-v1.1-mc-grch38.raw.vcf.gz` (vg deconstruct of the full, not frequency-filtered
  HPRC v1.1 graph: top-level sites, 44 diploid samples + CHM13 = 89 haplotypes, every allele kept). Indexed here through
  a symlink (`DATA/ref/`, tabix). The HPRC v2.1 full-graph VCF is not local (see "Needs"); the d46 runs are therefore
  also checked against the v1.1 haplotypes.
- **Assemblies**: HG008-T v3.2 (`HG008T_v3.2.fasta.gz`, both haplotypes), HG008-N v6.2
  (`HG008N_curatedv6_250714_bothhaps_polished6.2.fasta.gz`); COLO829BL verkko 2.1 hapX / hapY (ragtag GRCh38 scaffolds,
  primary + unlocalized; `/scratch/qfu/COLO829BL_DSA/ragtag_hg38`, read only). No COLO829T tumor assembly is local.
- **Gold set**: the 450 HG008-T somatic events that the HG008-N v6.2 assembly frame puts on a non-reference d9 node,
  233 with an exact HG008T_Illumina tensor key (`pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv`);
  for d46 the same events re-aligned to d46 (`gold_graph.py`, job `gold_d46.sh`: 18,211 events, class A 15,447,
  B (somatic edit off the GRCh38 path = the d46 gold set) 345, C 2,144, D 275; 179 class-B events with an exact
  tensor key). Keys with a tensor: d9 Illumina 231, PacBio 134, ONT 103 of 233; d46 Illumina 175, PacBio 119, ONT 97
  of 179 (job logs `gold_387266.out`, `gold_d46_387271.out`).

`DATA` = `/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation`.

## Method

Per run (`e9c_run.py VCFDIR NAME`), rows = the unplaced PASS calls (a random sample of 3,000, seed 0, when there are
more), 150 random placed controls of each kind (rtg TP after the PoN = somatic truth; rtg FP after the PoN; PASS calls
the GRCh38 PoN tags = population alleles), and the gold keys on chr1.

1. **Contexts from the tensor** (node-forward, up to 50 bases each side): `alt_read` = majority read bases of the A1
   rows (the sample's ALT haplotype); `alt_graph` = graph bases of the A1 rows' path (channel 5) with the site allele
   (channel 2); `ref_graph` = the graph path with REF; `ref_read` = majority read bases of the REF rows. The A1 / REF
   row blocks and the site columns are read from the tensor itself (channel 2 is non-zero exactly over the site
   columns); this equals the `variant_summary` layout on 783 of 783 checked tensors (383 SNV, 400 INDEL).
2. **Bubble locus**: from the call's node, a bp-weighted search along the graph edges in each direction, not through
   GRCh38 nodes, to the nearest node with one GRCh38 visit; the GRCh38 interval between the left and right anchors is
   the locus (`bubble`; one side only `one_side`; span > 100 kb `wide`, not used).
3. **Graph PoN (E9c)**. The ALT query = the event plus `f` bases each side, from `alt_graph` and from `alt_read`;
   `f` = 5, 10, 15 (11 / 21 / 31-mers for an SNV) or adaptive (10, widened by 5 up to 40 while the query occurs in the
   REF contexts, i.e. inside a repeat). The call is a **population allele** when
   - *in-graph*: a walk of the run's own graph within the call's neighbourhood (all nodes within the anchor distance
     plus the query length) spells the query (every d9 / d46 node and edge is carried by >= 10 % of the HPRC
     haplotypes), or
   - *HPRC haplotype*: GRCh38 or at least `min AC` of the 89 HPRC v1.1 haplotypes spell it inside the bubble interval
     +- 100 bp (haplotype sequences = GRCh38 with each haplotype's top-level alleles applied; phased GTs).

   **Main rule, fixed before looking at chr1: adaptive flank, min AC 1, graph or read context.** The grid of flanks and
   AC thresholds is reported for calibration (gold set on chr2-22) and on every group.
4. **Assembly check (E9d, `e9d_asm.py`)**: every context (>= 25 bp) mapped with `minimap2 -cx sr --cs --secondary=yes -N 50`
   to each assembly. A hit covers >= 90 % of the context, has no indel, at most one mismatch, and no mismatch within one
   base of the event. Primary query = `alt_read` (>= 40 bp, else `alt_graph`). Classes: **germline** (ALT in the normal
   assembly), **tumor-only** (ALT in HG008-T, not in HG008-N), **unsupported** (ALT in neither, REF context found),
   **unverifiable** (neither ALT nor REF context found); multi-copy flag when a context hits > 2 loci in one assembly.
   COLO829T: germline / not in normal / unverifiable (no tumor assembly).
5. **Categories and corrected precision (E9e part)**: exclusive order population > germline > tumor-only >
   unverifiable > unsupported. Precision on chr1 after the GRCh38 PoN at the PASS threshold, with TP / FP from the
   run's own rtg report (`eval_nobed`, `eval_bed`) and U unplaced PASS calls:
   `placed` TP/(TP+FP); `all FP` TP/(TP+FP+U); `graph PoN` TP/(TP+FP+U_untagged); `+tumor-only TP`
   (TP+S)/(TP+FP+U_untagged) with S the tumor-only unplaced calls (an upper bound: they are not truth-matched).
   Recall stays rtg's. In the BED, an unplaced call counts when its bubble interval lies inside the BED (no interval:
   counted, conservative). When the unplaced calls were sampled, the sample's fractions are scaled to all U calls.
6. **Cross-platform**: the same `node:offset:allele` key PASS among the unplaced calls of another platform's run of the
   same graph, sample and kind.
7. **Repeat context (R1.3)**: homopolymer >= 7 at the site, else a 2-6 bp unit tandem repeat of >= 3 copies, else other.
8. **E9a cross-tab** (`rcpu` runs that E9a has projected: the nine SNV runs): `collect.py` joins the E9a projection
   class of each unplaced PASS call (`../E9a_unplaced_projection`, `projected.tsv`: i = written as a GRCh38 variant,
   ii = a GRCh38 position or interval only, iii = no consistent anchor pair) to its E9c / E9d category.
9. **GRCh38 PoN reference point for the gold set** (`gold_grch38_pon.py`): the repository's PoN script
   (`scripts/filter_panel_of_normals.py`, unmodified, four PoN VCFs) on the GRCh38 truth records of the same gold events.

## Scripts and jobs

Code in this folder; job scripts and logs in `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E9c_offref_filter_validation/`
(`logs/`); outputs in `DATA`.

| step | script | job script | Slurm job | state | output |
|---|---|---|---|---|---|
| edge tables | `build_edges.py` | `edges.sh` | 387146 (d9), 387147 (d46) | COMPLETED | `DATA/edges_{d9,d46}` |
| pilot run, HG008 PacBio SNV d9 | `e9c_run.py` | none (run directly, 10:36) | none | done | `DATA/runs/legacy__HG008_PacBio_SNV_*`, log `DATA/test_pacbio.log` |
| legacy runs (11, `targets_20261013_104354.tsv`) | `e9c_run.py`, `e9c_lib.py` | `e9c_task.sh` | array 387265 | `_0`, `_1`, `_9` COMPLETED; `_2`-`_7` held, then CANCELLED 2026-10-15 13:05 without running; `_8`, `_10` FAILED | `DATA/runs/legacy__*` |
| d46 in-sample runs (3, `targets_d46_insample_20261013.tsv`) | `e9c_run.py` | `e9c_task.sh` (name `e9c_runs_d46`) | array 387431 | `_0`-`_2` COMPLETED | `DATA/runs/rcpu__*` |
| the other in-sample SNV runs (6, `targets_snv_rest_20261015.tsv`: HG008T Illumina / PacBio / ONT d9, HG008T Illumina d46, COLO829T Illumina / fiber-seq d46) | `e9c_run.py` | `e9c_task.sh` (name `e9c_runs_rest`) | array 405788 | `_0`-`_5` COMPLETED 2026-10-15 15:03-15:18 | `DATA/runs/rcpu__*` |
| gold set d9 | `gold.py` | `gold.sh` | 387266 | COMPLETED | `DATA/gold/HG008T_{Illumina,PacBio,ONT}` |
| gold set d46 | `gold_graph.py` (vg 1.77 giraffe of the assembly windows), `gold.py --graph d46` | `gold_d46.sh` | 387271 | COMPLETED | `DATA/gold_d46_mapping`, `DATA/gold/HG008T_*_d46` |
| GRCh38 PoN on the gold truth records | `gold_grch38_pon.py` | none (run directly, 11:34-11:36) | none | done | `DATA/gold_grch38_pon` |
| assembly check pilot | `e9d_asm.py` | `asm.sh` | 387268 (`e9d_asm_pilot`) | CANCELLED (superseded) | |
| assembly check HG008T (12,059 contexts; 7 runs) | `e9d_asm.py` (minimap2 2.28) | `asm.sh HG008T` | 387432 | COMPLETED (superseded) | `DATA/asm/HG008T` |
| assembly check COLO829T (15,509 contexts; 7 runs) | `e9d_asm.py` | `asm.sh COLO829T` | 387433 | COMPLETED (superseded) | `DATA/asm/COLO829T` |
| tables (7 runs) | `collect.py` | `collect.sh` | 387434 | COMPLETED 2026-10-13 16:07 (superseded) | |
| assembly check HG008T (18,217 contexts; 13 runs + gold) | `e9d_asm.py` | `asm.sh HG008T` | 405793 | COMPLETED 2026-10-15 15:22 | `DATA/asm/HG008T` |
| assembly check COLO829T (18,327 contexts) | `e9d_asm.py` | `asm.sh COLO829T` | 405794 | COMPLETED 15:27 | `DATA/asm/COLO829T` |
| tables (13 runs, 6 gold sets) | `collect.py` | `collect.sh` | 405796 | COMPLETED 15:27 | `tables.md`, `summary.json`, `per_call.tsv.gz` |

`run_all.sh` is the driver (discover, submit the runs, gold, assembly checks and collect with dependencies); `discover.py`
lists the vcf directories that are ready.

## Results (2026-10-15; 13 runs: 9 RIS 7fc7396 in-sample SNV on d9 / d46, 4 legacy d9)

All numbers are chr1 (the test chromosome) unless a table says otherwise, after the run's GRCh38 PoN, at the run's PASS
threshold. Tables are copied from `tables.md` (generated by `collect.py`, job 405796); the full flank x AC grid is in
`tables.md` section 6, per-call rows in `per_call.tsv.gz`, everything in `summary.json`. Run names: `legacy__` = local
pre-RIS d9 runs, `rcpu__` = RIS 7fc7396 re-predictions (d9 = HPRC v1.1, d46 = HPRC v2.1). The GRCh38-only arm has no
unplaced calls; its numbers are in `../E9a_unplaced_projection/` (Results §7).

### 1. What the chr1 unplaced PASS calls are

Exclusive categories (population > germline > tumor-only > unverifiable > unsupported), fractions of the calls with a
tensor context (`tables.md` section 1):

| run | graph | kind | sample | platform | unplaced PASS | processed | with ctx | population (graph PoN) | in-graph | HPRC hap | germline (normal asm) | tumor-only | not in normal | unverifiable | unsupported | not checked |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy__COLO829T_ONT_SNV_pow075_S_lr1e4_e20 | d9 | SNV | COLO829T | ONT | 2515 | 2515 | 2505 | 1892 (75.5 %) | 564 (22.5 %) | 1886 (75.3 %) | 445 (17.8 %) | 0 (0.0 %) | 59 (2.4 %) | 109 (4.4 %) | 0 (0.0 %) | 0 (0.0 %) |
| legacy__HG008_Illumina_INDEL_small_nofilt | d9 | INDEL | HG008T | Illumina | 534 | 534 | 534 | 478 (89.5 %) | 39 (7.3 %) | 477 (89.3 %) | 37 (6.9 %) | 9 (1.7 %) | 0 (0.0 %) | 4 (0.7 %) | 6 (1.1 %) | 0 (0.0 %) |
| legacy__HG008_Illumina_SNV_pow075_S_lr1e4_e20 | d9 | SNV | HG008T | Illumina | 595 | 595 | 595 | 469 (78.8 %) | 143 (24.0 %) | 461 (77.5 %) | 122 (20.5 %) | 2 (0.3 %) | 0 (0.0 %) | 1 (0.2 %) | 1 (0.2 %) | 0 (0.0 %) |
| legacy__HG008_PacBio_SNV_pow075_S_lr1e4_e20 | d9 | SNV | HG008T | PacBio | 105 | 105 | 105 | 88 (83.8 %) | 15 (14.3 %) | 88 (83.8 %) | 15 (14.3 %) | 1 (1.0 %) | 0 (0.0 %) | 1 (1.0 %) | 0 (0.0 %) | 0 (0.0 %) |
| rcpu__COLO829_Illumina_SNV_d46_7fc7396 | d46 | SNV | COLO829T | Illumina | 516 | 516 | 516 | 276 (53.5 %) | 90 (17.4 %) | 233 (45.2 %) | 222 (43.0 %) | 0 (0.0 %) | 16 (3.1 %) | 2 (0.4 %) | 0 (0.0 %) | 0 (0.0 %) |
| rcpu__COLO829_ONT_SNV_d46_7fc7396 | d46 | SNV | COLO829T | ONT | 4569 | 3000 | 3000 | 1495 (49.8 %) | 272 (9.1 %) | 1466 (48.9 %) | 1362 (45.4 %) | 0 (0.0 %) | 72 (2.4 %) | 71 (2.4 %) | 0 (0.0 %) | 0 (0.0 %) |
| rcpu__COLO829_fiberseq_SNV_d46_7fc7396 | d46 | SNV | COLO829T | fiberseq | 257 | 257 | 257 | 200 (77.8 %) | 35 (13.6 %) | 197 (76.7 %) | 27 (10.5 %) | 0 (0.0 %) | 15 (5.8 %) | 15 (5.8 %) | 0 (0.0 %) | 0 (0.0 %) |
| rcpu__HG008_Illumina_SNV_7fc7396_ddp4_b1024 | d9 | SNV | HG008T | Illumina | 115 | 115 | 115 | 95 (82.6 %) | 10 (8.7 %) | 95 (82.6 %) | 19 (16.5 %) | 1 (0.9 %) | 0 (0.0 %) | 0 (0.0 %) | 0 (0.0 %) | 0 (0.0 %) |
| rcpu__HG008_Illumina_SNV_d46_7fc7396 | d46 | SNV | HG008T | Illumina | 754 | 754 | 754 | 496 (65.8 %) | 209 (27.7 %) | 425 (56.4 %) | 250 (33.2 %) | 4 (0.5 %) | 0 (0.0 %) | 2 (0.3 %) | 2 (0.3 %) | 0 (0.0 %) |
| rcpu__HG008_ONT_SNV_7fc7396 | d9 | SNV | HG008T | ONT | 1410 | 1410 | 1410 | 1289 (91.4 %) | 189 (13.4 %) | 1285 (91.1 %) | 59 (4.2 %) | 3 (0.2 %) | 0 (0.0 %) | 52 (3.7 %) | 7 (0.5 %) | 0 (0.0 %) |
| rcpu__HG008_ONT_SNV_d46_7fc7396 | d46 | SNV | HG008T | ONT | 1258 | 1258 | 1258 | 932 (74.1 %) | 174 (13.8 %) | 925 (73.5 %) | 242 (19.2 %) | 3 (0.2 %) | 0 (0.0 %) | 70 (5.6 %) | 11 (0.9 %) | 0 (0.0 %) |
| rcpu__HG008_PacBio_SNV_7fc7396 | d9 | SNV | HG008T | PacBio | 210 | 210 | 210 | 186 (88.6 %) | 28 (13.3 %) | 186 (88.6 %) | 15 (7.1 %) | 1 (0.5 %) | 0 (0.0 %) | 7 (3.3 %) | 1 (0.5 %) | 0 (0.0 %) |
| rcpu__HG008_PacBio_SNV_d46_7fc7396 | d46 | SNV | HG008T | PacBio | 240 | 240 | 240 | 214 (89.2 %) | 47 (19.6 %) | 213 (88.8 %) | 19 (7.9 %) | 1 (0.4 %) | 0 (0.0 %) | 5 (2.1 %) | 1 (0.4 %) | 0 (0.0 %) |

`in-graph` and `HPRC hap` overlap; `population` is their union.

Assembly classes alone, before the population precedence (`tables.md` section 3, `unplaced_pass` rows; calls tagged by
the graph PoN counted from `per_call.tsv.gz`; multi-copy from `summary.json`):

| run | with ctx | germline (asm) | tumor-only (asm) | of which graph-PoN tagged | not in normal (COLO829T) | of which tagged | unsupported | unverifiable | not checked | multi-copy context |
|---|---|---|---|---|---|---|---|---|---|---|
| legacy COLO829T ONT SNV d9 | 2505 | 1814 (72.4 %) | - | - | 235 | 176 | - | 456 | 0 | 1298 |
| legacy HG008T Illumina INDEL d9 | 534 | 434 (81.3 %) | 58 (10.9 %) | 49 | - | - | 22 | 20 | 0 | 13 |
| legacy HG008T Illumina SNV d9 | 595 | 555 (93.3 %) | 17 (2.9 %) | 15 | - | - | 8 | 15 | 0 | 157 |
| legacy HG008T PacBio SNV d9 | 105 | 78 (74.3 %) | 9 (8.6 %) | 8 | - | - | 2 | 16 | 0 | 8 |
| rcpu COLO829T Illumina SNV d46 | 516 | 487 (94.4 %) | - | - | 23 | 7 | - | 6 | 0 | 196 |
| rcpu COLO829T fiber-seq SNV d46 | 257 | 137 (53.3 %) | - | - | 38 | 23 | - | 82 | 0 | 48 |
| rcpu COLO829T ONT SNV d46 | 3000 | 2493 (83.1 %) | - | - | 197 | 125 | - | 309 | 1 | 1594 |
| rcpu HG008T Illumina SNV d9 | 115 | 111 (96.5 %) | 4 (3.5 %) | 3 | - | - | 0 | 0 | 0 | 4 |
| rcpu HG008T Illumina SNV d46 | 754 | 700 (92.8 %) | 24 (3.2 %) | 20 | - | - | 20 | 10 | 0 | 245 |
| rcpu HG008T PacBio SNV d9 | 210 | 128 (61.0 %) | 11 (5.2 %) | 10 | - | - | 13 | 58 | 0 | 6 |
| rcpu HG008T PacBio SNV d46 | 240 | 158 (65.8 %) | 9 (3.8 %) | 8 | - | - | 7 | 66 | 0 | 8 |
| rcpu HG008T ONT SNV d9 | 1410 | 877 (62.2 %) | 16 (1.1 %) | 13 | - | - | 76 | 440 | 1 | 173 |
| rcpu HG008T ONT SNV d46 | 1258 | 724 (57.6 %) | 22 (1.7 %) | 19 | - | - | 101 | 410 | 1 | 279 |

### 2. Corrected chr1 precision (after the GRCh38 PoN, at the run's PASS threshold)

`tables.md` section 2. `placed`: rtg on placed calls; `all FP`: every unplaced PASS call false; `graph PoN`: tagged
unplaced calls removed, the rest false; `+tumor-only TP`: the tumor-only remainder counted as true (upper bound). Recall
is rtg's. COLO829T rows have no HG008-T-style tumor assembly, so `+tumor-only` equals `graph PoN`, and their BED covers
chr1, so `bed` equals `nobed`.

| run | eval | TP | FP | unplaced PASS | unplaced share of PASS | recall | P placed | F1 placed | P all FP | F1 all FP | P graph PoN | F1 graph PoN | P +tumor-only | F1 +tumor-only | unplaced left | unplaced share left |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy__COLO829T_ONT_SNV_pow075_S_lr1e4_e20 | bed | 1765 | 691 | 2515.0 | 0.5059 | 0.7225 | 0.7186 | 0.7206 | 0.3551 | 0.4761 | 0.5746 | 0.6401 | 0.5746 | 0.6401 | 615.4 | 0.2004 |
| legacy__COLO829T_ONT_SNV_pow075_S_lr1e4_e20 | nobed | 1765 | 691 | 2515 | 0.5059 | 0.7225 | 0.7186 | 0.7206 | 0.3551 | 0.4761 | 0.5746 | 0.6401 | 0.5746 | 0.6401 | 615.4 | 0.2004 |
| legacy__HG008_Illumina_INDEL_small_nofilt | bed | 137 | 423 | 475.0 | 0.4589 | 0.2334 | 0.2446 | 0.2389 | 0.1324 | 0.1689 | 0.225 | 0.2291 | 0.2381 | 0.2357 | 49.0 | 0.0805 |
| legacy__HG008_Illumina_INDEL_small_nofilt | nobed | 145 | 477 | 534 | 0.4619 | 0.2283 | 0.2331 | 0.2307 | 0.1254 | 0.1619 | 0.2139 | 0.2209 | 0.2271 | 0.2277 | 56.0 | 0.0826 |
| legacy__HG008_Illumina_SNV_pow075_S_lr1e4_e20 | bed | 584 | 702 | 425.0 | 0.2484 | 0.8379 | 0.4541 | 0.589 | 0.3413 | 0.485 | 0.4339 | 0.5717 | 0.4354 | 0.573 | 60.0 | 0.0446 |
| legacy__HG008_Illumina_SNV_pow075_S_lr1e4_e20 | nobed | 588 | 745 | 595 | 0.3086 | 0.8376 | 0.4411 | 0.5779 | 0.305 | 0.4471 | 0.403 | 0.5442 | 0.4044 | 0.5454 | 126.0 | 0.0864 |
| legacy__HG008_PacBio_SNV_pow075_S_lr1e4_e20 | bed | 579 | 588 | 88.0 | 0.0701 | 0.8307 | 0.4961 | 0.6212 | 0.4614 | 0.5932 | 0.4907 | 0.6169 | 0.4915 | 0.6176 | 13.0 | 0.011 |
| legacy__HG008_PacBio_SNV_pow075_S_lr1e4_e20 | nobed | 583 | 617 | 105 | 0.0805 | 0.8305 | 0.4858 | 0.613 | 0.4467 | 0.581 | 0.479 | 0.6076 | 0.4799 | 0.6083 | 17.0 | 0.014 |
| rcpu__COLO829_Illumina_SNV_d46_7fc7396 | bed | 1745 | 268 | 516.0 | 0.204 | 0.7143 | 0.8669 | 0.7832 | 0.69 | 0.7019 | 0.7745 | 0.7432 | 0.7745 | 0.7432 | 240.0 | 0.1065 |
| rcpu__COLO829_Illumina_SNV_d46_7fc7396 | nobed | 1745 | 268 | 516 | 0.204 | 0.7143 | 0.8669 | 0.7832 | 0.69 | 0.7019 | 0.7745 | 0.7432 | 0.7745 | 0.7432 | 240.0 | 0.1065 |
| rcpu__COLO829_ONT_SNV_d46_7fc7396 | bed | 1690 | 457 | 4569.0 | 0.6803 | 0.6918 | 0.7871 | 0.7364 | 0.2516 | 0.369 | 0.3807 | 0.4911 | 0.3807 | 0.4911 | 2292.1 | 0.5163 |
| rcpu__COLO829_ONT_SNV_d46_7fc7396 | nobed | 1690 | 457 | 4569 | 0.6803 | 0.6918 | 0.7871 | 0.7364 | 0.2516 | 0.369 | 0.3807 | 0.4911 | 0.3807 | 0.4911 | 2292.1 | 0.5163 |
| rcpu__COLO829_fiberseq_SNV_d46_7fc7396 | bed | 1728 | 268 | 257.0 | 0.1141 | 0.7073 | 0.8657 | 0.7786 | 0.767 | 0.7359 | 0.8417 | 0.7687 | 0.8417 | 0.7687 | 57.0 | 0.0278 |
| rcpu__COLO829_fiberseq_SNV_d46_7fc7396 | nobed | 1728 | 268 | 257 | 0.1141 | 0.7073 | 0.8657 | 0.7786 | 0.767 | 0.7359 | 0.8417 | 0.7687 | 0.8417 | 0.7687 | 57.0 | 0.0278 |
| rcpu__HG008_Illumina_SNV_7fc7396_ddp4_b1024 | bed | 579 | 702 | 104.0 | 0.0751 | 0.8307 | 0.452 | 0.5854 | 0.4181 | 0.5562 | 0.4464 | 0.5807 | 0.4472 | 0.5814 | 16.0 | 0.0123 |
| rcpu__HG008_Illumina_SNV_7fc7396_ddp4_b1024 | nobed | 583 | 734 | 115 | 0.0803 | 0.8305 | 0.4427 | 0.5775 | 0.4071 | 0.5464 | 0.4361 | 0.5718 | 0.4368 | 0.5725 | 20.0 | 0.015 |
| rcpu__HG008_Illumina_SNV_d46_7fc7396 | bed | 594 | 792 | 634.0 | 0.3139 | 0.8522 | 0.4286 | 0.5703 | 0.2941 | 0.4372 | 0.3738 | 0.5197 | 0.3763 | 0.5221 | 203.0 | 0.1278 |
| rcpu__HG008_Illumina_SNV_d46_7fc7396 | nobed | 598 | 833 | 754 | 0.3451 | 0.8519 | 0.4179 | 0.5607 | 0.2737 | 0.4143 | 0.3541 | 0.5002 | 0.3564 | 0.5026 | 258.0 | 0.1528 |
| rcpu__HG008_ONT_SNV_7fc7396 | bed | 576 | 660 | 1219.0 | 0.4965 | 0.8264 | 0.466 | 0.596 | 0.2346 | 0.3655 | 0.4331 | 0.5683 | 0.4353 | 0.5703 | 94.0 | 0.0707 |
| rcpu__HG008_ONT_SNV_7fc7396 | nobed | 580 | 697 | 1410 | 0.5247 | 0.8262 | 0.4542 | 0.5862 | 0.2159 | 0.3423 | 0.4149 | 0.5524 | 0.417 | 0.5543 | 121.0 | 0.0866 |
| rcpu__HG008_ONT_SNV_d46_7fc7396 | bed | 599 | 788 | 1182.0 | 0.4601 | 0.8594 | 0.4319 | 0.5749 | 0.2332 | 0.3668 | 0.3513 | 0.4988 | 0.3531 | 0.5005 | 318.0 | 0.1865 |
| rcpu__HG008_ONT_SNV_d46_7fc7396 | nobed | 603 | 836 | 1258 | 0.4664 | 0.859 | 0.419 | 0.5633 | 0.2236 | 0.3548 | 0.3416 | 0.4889 | 0.3433 | 0.4906 | 326.0 | 0.1847 |
| rcpu__HG008_PacBio_SNV_7fc7396 | bed | 577 | 583 | 192.0 | 0.142 | 0.8278 | 0.4974 | 0.6214 | 0.4268 | 0.5632 | 0.489 | 0.6148 | 0.4898 | 0.6155 | 20.0 | 0.0169 |
| rcpu__HG008_PacBio_SNV_7fc7396 | nobed | 581 | 608 | 210 | 0.1501 | 0.8276 | 0.4886 | 0.6145 | 0.4153 | 0.5531 | 0.479 | 0.6068 | 0.4798 | 0.6074 | 24.0 | 0.0198 |
| rcpu__HG008_PacBio_SNV_d46_7fc7396 | bed | 582 | 650 | 226.0 | 0.155 | 0.835 | 0.4724 | 0.6034 | 0.3992 | 0.5401 | 0.4626 | 0.5954 | 0.4634 | 0.5961 | 26.0 | 0.0207 |
| rcpu__HG008_PacBio_SNV_d46_7fc7396 | nobed | 586 | 672 | 240 | 0.1602 | 0.8348 | 0.4658 | 0.598 | 0.3912 | 0.5327 | 0.4564 | 0.5901 | 0.4572 | 0.5908 | 26.0 | 0.0202 |

### 3. Controls: graph PoN tag rate (main rule) and assembly classes

`tables.md` section 3, condensed. `ctrl_tp`: placed PASS calls rtg matched to the somatic truth; `ctrl_fp`: placed false
calls after the PoN; `ctrl_pon`: placed PASS calls the GRCh38 PoN tags; random samples of up to 150 each. The assembly
column gives the main class of each group (tumor-only / not in normal for ctrl_tp, germline for ctrl_fp and ctrl_pon).

| run | ctrl_tp tagged | ctrl_tp: tumor-only / not in normal | ctrl_fp tagged | ctrl_fp: germline | ctrl_pon tagged | ctrl_pon: germline |
|---|---|---|---|---|---|---|
| legacy COLO829T ONT SNV d9 | 0 / 146 | 144 | 10 / 125 (8.0 %) | 87 | 2 / 137 (1.5 %) | 132 |
| legacy HG008T Illumina INDEL d9 | 35 / 145 (24.1 %) | 144 | 71 / 150 (47.3 %) | 137 | 127 / 150 (84.7 %) | 141 |
| legacy HG008T Illumina SNV d9 | 0 / 150 | 148 | 7 / 150 (4.7 %) | 141 | 9 / 150 (6.0 %) | 145 |
| legacy HG008T PacBio SNV d9 | 0 / 150 | 146 | 4 / 150 (2.7 %) | 144 | 8 / 150 (5.3 %) | 149 |
| rcpu COLO829T Illumina SNV d46 | 0 / 150 | 148 | 2 / 150 (1.3 %) | 125 | 7 / 150 (4.7 %) | 136 |
| rcpu COLO829T fiber-seq SNV d46 | 0 / 150 | 148 | 1 / 150 (0.7 %) | 131 | 7 / 150 (4.7 %) | 140 |
| rcpu COLO829T ONT SNV d46 | 0 / 150 | 150 | 3 / 150 (2.0 %) | 110 | 11 / 150 (7.3 %) | 144 |
| rcpu HG008T Illumina SNV d9 | 0 / 150 | 147 | 18 / 150 (12.0 %) | 145 | 10 / 150 (6.7 %) | 147 |
| rcpu HG008T Illumina SNV d46 | 0 / 150 | 149 | 18 / 150 (12.0 %) | 140 | 24 / 150 (16.0 %) | 142 |
| rcpu HG008T PacBio SNV d9 | 0 / 150 | 147 | 4 / 150 (2.7 %) | 145 | 15 / 150 (10.0 %) | 148 |
| rcpu HG008T PacBio SNV d46 | 0 / 150 | 149 | 4 / 150 (2.7 %) | 146 | 17 / 150 (11.3 %) | 147 |
| rcpu HG008T ONT SNV d9 | 0 / 150 | 146 | 1 / 150 (0.7 %) | 139 | 12 / 150 (8.0 %) | 145 |
| rcpu HG008T ONT SNV d46 | 1 / 150 (0.7 %) | 145 | 15 / 150 (10.0 %) | 135 | 13 / 150 (8.7 %) | 140 |

Main rule against the most conservative grid point (adaptive flank, AC >= 9), tagged / n (`tables.md` section 6):

| run | unplaced PASS, AC>=1 | unplaced PASS, AC>=9 | ctrl_tp, AC>=1 | ctrl_tp, AC>=9 | ctrl_pon, AC>=1 | ctrl_pon, AC>=9 | gold, AC>=1 | gold, AC>=9 |
|---|---|---|---|---|---|---|---|---|
| legacy COLO829T ONT SNV d9 | 1892/2505 | 1557/2505 | 0/146 | 0/146 | 2/137 | 1/137 | 2/2 | 0/2 |
| legacy HG008T Illumina INDEL d9 | 478/534 | 85/534 | 35/145 | 1/145 | 127/150 | 10/150 | 3/3 | 0/3 |
| legacy HG008T Illumina SNV d9 | 469/595 | 336/595 | 0/150 | 0/150 | 9/150 | 5/150 | 6/10 | 3/10 |
| legacy HG008T PacBio SNV d9 | 88/105 | 60/105 | 0/150 | 0/150 | 8/150 | 1/150 | 3/7 | 1/7 |
| rcpu COLO829T Illumina SNV d46 | 276/516 | 219/516 | 0/150 | 0/150 | 7/150 | 0/150 | 1/1 | 1/1 |
| rcpu COLO829T fiber-seq SNV d46 | 200/257 | 159/257 | 0/150 | 0/150 | 7/150 | 1/150 | 2/2 | 1/2 |
| rcpu COLO829T ONT SNV d46 | 1495/3000 | 1146/3000 | 0/150 | 0/150 | 11/150 | 2/150 | 2/2 | 1/2 |
| rcpu HG008T Illumina SNV d9 | 95/115 | 46/115 | 0/150 | 0/150 | 10/150 | 1/150 | 6/10 | 3/10 |
| rcpu HG008T Illumina SNV d46 | 496/754 | 376/754 | 0/150 | 0/150 | 24/150 | 5/150 | 9/12 | 7/12 |
| rcpu HG008T PacBio SNV d9 | 186/210 | 131/210 | 0/150 | 0/150 | 15/150 | 0/150 | 3/7 | 1/7 |
| rcpu HG008T PacBio SNV d46 | 214/240 | 160/240 | 0/150 | 0/150 | 17/150 | 1/150 | 9/11 | 6/11 |
| rcpu HG008T ONT SNV d9 | 1289/1410 | 1155/1410 | 0/150 | 0/150 | 12/150 | 1/150 | 3/7 | 1/7 |
| rcpu HG008T ONT SNV d46 | 932/1258 | 803/1258 | 1/150 | 0/150 | 13/150 | 4/150 | 7/9 | 5/9 |

### 4. Gold set: false-tag rate of the graph PoN on assembly-verified off-reference somatic keys (all chromosomes)

Main rule (adaptive flank, min AC 1), from `tables.md` section 4; chr2-22 is the calibration part, chr1 is held out:

| set | part | n | tagged | of which in-graph |
|---|---|---|---|---|
| HG008T_Illumina (d9) | chr2-22 | 218 | 195 (89.4 %) | 44 |
| HG008T_Illumina (d9) | chr1 | 13 | 9 (69.2 %) | 1 |
| HG008T_Illumina (d9) | all | 231 | 204 (88.3 %) | 45 |
| HG008T_PacBio (d9) | chr2-22 | 126 | 107 (84.9 %) | 20 |
| HG008T_PacBio (d9) | chr1 | 8 | 4 (50.0 %) | 1 |
| HG008T_PacBio (d9) | all | 134 | 111 (82.8 %) | 21 |
| HG008T_ONT (d9) | chr2-22 | 96 | 81 (84.4 %) | 16 |
| HG008T_ONT (d9) | chr1 | 7 | 3 (42.9 %) | 1 |
| HG008T_ONT (d9) | all | 103 | 84 (81.6 %) | 17 |
| HG008T_Illumina_d46 | chr2-22 | 160 | 127 (79.4 %) | 37 |
| HG008T_Illumina_d46 | chr1 | 15 | 12 (80.0 %) | 3 |
| HG008T_Illumina_d46 | all | 175 | 139 (79.4 %) | 40 |
| HG008T_PacBio_d46 | chr2-22 | 105 | 78 (74.3 %) | 21 |
| HG008T_PacBio_d46 | chr1 | 14 | 12 (85.7 %) | 4 |
| HG008T_PacBio_d46 | all | 119 | 90 (75.6 %) | 25 |
| HG008T_ONT_d46 | chr2-22 | 87 | 64 (73.6 %) | 12 |
| HG008T_ONT_d46 | chr1 | 10 | 8 (80.0 %) | 4 |
| HG008T_ONT_d46 | all | 97 | 72 (74.2 %) | 16 |

The lowest false-tag rate anywhere in the grid on chr2-22 is at adaptive flank, AC >= 9, for every set: d9 Illumina
146/218 (67.0 %), PacBio 75/126 (59.5 %), ONT 57/96 (59.4 %); d46 Illumina 89/160 (55.6 %), PacBio 48/105 (45.7 %),
ONT 37/87 (42.5 %). The largest is at flank 5, AC >= 1 (d9 Illumina 209/218, 95.9 %).

Graph PoN (main rule) against the GRCh38 PoN (repository rule) on the same events' GRCh38 truth records, and on the
gold keys whose ALT context is tumor-only in the assemblies (`tables.md` section 4, second table):

| set | part | n | graph PoN | GRCh38 PoN | both | neither | tumor-only gold n | graph PoN | GRCh38 PoN |
|---|---|---|---|---|---|---|---|---|---|
| HG008T_Illumina | chr2-22 | 218 | 195 (89.4 %) | 159 (72.9 %) | 152 | 16 | 143 | 123 (86.0 %) | 102 (71.3 %) |
| HG008T_Illumina | chr1 | 13 | 9 (69.2 %) | 9 (69.2 %) | 9 | 4 | 12 | 8 (66.7 %) | 8 (66.7 %) |
| HG008T_Illumina | all | 231 | 204 (88.3 %) | 168 (72.7 %) | 161 | 20 | 155 | 131 (84.5 %) | 110 (71.0 %) |
| HG008T_Illumina_d46 | chr2-22 | 160 | 127 (79.4 %) | 47 (29.4 %) | 45 | 31 | 123 | 91 (74.0 %) | 32 (26.0 %) |
| HG008T_Illumina_d46 | chr1 | 15 | 12 (80.0 %) | 4 (26.7 %) | 4 | 3 | 14 | 11 (78.6 %) | 3 (21.4 %) |
| HG008T_Illumina_d46 | all | 175 | 139 (79.4 %) | 51 (29.1 %) | 49 | 34 | 137 | 102 (74.5 %) | 35 (25.5 %) |
| HG008T_ONT | chr2-22 | 96 | 81 (84.4 %) | 61 (63.5 %) | 59 | 13 | 61 | 49 (80.3 %) | 36 (59.0 %) |
| HG008T_ONT | chr1 | 7 | 3 (42.9 %) | 3 (42.9 %) | 3 | 4 | 6 | 3 (50.0 %) | 3 (50.0 %) |
| HG008T_ONT | all | 103 | 84 (81.6 %) | 64 (62.1 %) | 62 | 17 | 67 | 52 (77.6 %) | 39 (58.2 %) |
| HG008T_ONT_d46 | chr2-22 | 87 | 64 (73.6 %) | 23 (26.4 %) | 21 | 21 | 72 | 49 (68.1 %) | 17 (23.6 %) |
| HG008T_ONT_d46 | chr1 | 10 | 8 (80.0 %) | 3 (30.0 %) | 3 | 2 | 9 | 7 (77.8 %) | 2 (22.2 %) |
| HG008T_ONT_d46 | all | 97 | 72 (74.2 %) | 26 (26.8 %) | 24 | 23 | 81 | 56 (69.1 %) | 19 (23.5 %) |
| HG008T_PacBio | chr2-22 | 126 | 107 (84.9 %) | 86 (68.3 %) | 80 | 13 | 88 | 71 (80.7 %) | 55 (62.5 %) |
| HG008T_PacBio | chr1 | 8 | 4 (50.0 %) | 4 (50.0 %) | 4 | 4 | 8 | 4 (50.0 %) | 4 (50.0 %) |
| HG008T_PacBio | all | 134 | 111 (82.8 %) | 90 (67.2 %) | 84 | 17 | 96 | 75 (78.1 %) | 59 (61.5 %) |
| HG008T_PacBio_d46 | chr2-22 | 105 | 78 (74.3 %) | 26 (24.8 %) | 24 | 25 | 91 | 64 (70.3 %) | 21 (23.1 %) |
| HG008T_PacBio_d46 | chr1 | 14 | 12 (85.7 %) | 5 (35.7 %) | 5 | 2 | 12 | 10 (83.3 %) | 3 (25.0 %) |
| HG008T_PacBio_d46 | all | 119 | 90 (75.6 %) | 31 (26.1 %) | 29 | 27 | 103 | 74 (71.8 %) | 24 (23.3 %) |

Assembly classes and repeat context of the gold keys (all chromosomes; the gold events are somatic by construction):

| set | n | assembly classes | repeat context |
|---|---|---|---|
| HG008T_Illumina (d9) | 231 | germline 60, tumor_only 155, unsupported 8, unverifiable 8 | STR 58, homopolymer>=7 107, other 66 |
| HG008T_PacBio (d9) | 134 | germline 26, tumor_only 96, unsupported 3, unverifiable 9 | STR 26, homopolymer>=7 73, other 35 |
| HG008T_ONT (d9) | 103 | germline 19, tumor_only 67, unsupported 7, unverifiable 10 | STR 21, homopolymer>=7 54, other 28 |
| HG008T_Illumina_d46 | 175 | germline 33, tumor_only 137, unsupported 4, unverifiable 1 | STR 40, homopolymer>=7 101, other 34 |
| HG008T_PacBio_d46 | 119 | germline 11, tumor_only 103, unsupported 3, unverifiable 2 | STR 25, homopolymer>=7 82, other 12 |
| HG008T_ONT_d46 | 97 | germline 9, tumor_only 81, unsupported 5, unverifiable 2 | STR 18, homopolymer>=7 67, other 12 |

### 5. Gold keys on chr1 in each run (graph-only somatic sensitivity)

`tables.md` section 5:

| run | gold keys with a tensor | PASS (unplaced) | tagged | PASS and tagged | assembly classes |
|---|---|---|---|---|---|
| legacy__COLO829T_ONT_SNV_pow075_S_lr1e4_e20 | 2 | 1 | 2 | 1 | unverifiable 2 |
| legacy__HG008_Illumina_INDEL_small_nofilt | 3 | 3 | 3 | 3 | tumor_only 3 |
| legacy__HG008_Illumina_SNV_pow075_S_lr1e4_e20 | 10 | 4 | 6 | 3 | germline 1, tumor_only 9 |
| legacy__HG008_PacBio_SNV_pow075_S_lr1e4_e20 | 7 | 3 | 3 | 2 | tumor_only 7 |
| rcpu__COLO829_Illumina_SNV_d46_7fc7396 | 1 | 0 | 1 | 0 | unverifiable 1 |
| rcpu__COLO829_ONT_SNV_d46_7fc7396 | 2 | 1 | 2 | 1 | germline 1, unverifiable 1 |
| rcpu__COLO829_fiberseq_SNV_d46_7fc7396 | 2 | 0 | 2 | 0 | germline 1, unverifiable 1 |
| rcpu__HG008_Illumina_SNV_7fc7396_ddp4_b1024 | 10 | 2 | 6 | 1 | germline 1, tumor_only 9 |
| rcpu__HG008_Illumina_SNV_d46_7fc7396 | 12 | 9 | 9 | 7 | germline 1, tumor_only 11 |
| rcpu__HG008_ONT_SNV_7fc7396 | 7 | 4 | 3 | 2 | tumor_only 6, unsupported 1 |
| rcpu__HG008_ONT_SNV_d46_7fc7396 | 9 | 5 | 7 | 5 | germline 1, tumor_only 8 |
| rcpu__HG008_PacBio_SNV_7fc7396 | 7 | 3 | 3 | 2 | tumor_only 7 |
| rcpu__HG008_PacBio_SNV_d46_7fc7396 | 11 | 3 | 9 | 3 | germline 2, tumor_only 9 |

(The gold keys are HG008-T events; in the COLO829T runs they are the same graph keys found among the COLO829T
candidates.) The d9 gold key in the bubble chr1:78,932,048-78,932,236 is an untagged tumor-only unplaced PASS call in
all three HG008T d9 runs; E9a writes it as the truth SNV chr1:78,932,081 A>T (class i, an rtg TP after the PoN).

### 6. E9a projection class of the unplaced PASS calls against the E9c / E9d category (`rcpu` runs)

From `summary.json` (`e9a_proj_class_of_unplaced_pass`, `e9a_category_by_class`); i = written as a GRCh38 variant by
E9a, ii = GRCh38 position or interval only, iii = no consistent anchor pair:

| run | class i | class ii | class iii |
|---|---|---|---|
| rcpu COLO829T Illumina SNV d46 | 39: germline 20, population 17, not_in_normal 1, unverifiable 1 | 250: population 195, germline 44, not_in_normal 10, unverifiable 1 | 227: germline 158, population 64, not_in_normal 5 |
| rcpu COLO829T fiber-seq SNV d46 | 12: population 8, germline 2, not_in_normal 2 | 203: population 164, unverifiable 15, germline 13, not_in_normal 11 | 42: population 28, germline 12, not_in_normal 2 |
| rcpu COLO829T ONT SNV d46 (3,000 sampled) | 85: population 51, germline 33, unverifiable 1 | 885: population 717, unverifiable 70, germline 57, not_in_normal 41 | 2030: germline 1272, population 727, not_in_normal 31 |
| rcpu HG008T Illumina SNV d9 | 22: population 16, germline 5, tumor_only 1 | 89: population 78, germline 11 | 4: germline 3, population 1 |
| rcpu HG008T Illumina SNV d46 | 58: population 36, germline 20, tumor_only 2 | 440: population 386, germline 49, tumor_only 2, unsupported 2, unverifiable 1 | 256: germline 181, population 74, unverifiable 1 |
| rcpu HG008T PacBio SNV d9 | 14: population 10, germline 3, tumor_only 1 | 191: population 174, germline 9, unverifiable 7, unsupported 1 | 5: germline 3, population 2 |
| rcpu HG008T PacBio SNV d46 | 19: population 13, germline 6 | 221: population 201, germline 13, unverifiable 5, tumor_only 1, unsupported 1 | 0 |
| rcpu HG008T ONT SNV d9 | 47: population 36, germline 6, tumor_only 2, unverifiable 2, unsupported 1 | 1238: population 1154, unverifiable 49, germline 28, unsupported 6, tumor_only 1 | 125: population 99, germline 25, unverifiable 1 |
| rcpu HG008T ONT SNV d46 | 43: population 32, germline 8, unverifiable 2, unsupported 1 | 984: population 880, unverifiable 68, germline 25, unsupported 10, tumor_only 1 | 231: germline 209, population 20, tumor_only 2 |

### 7. Repeat context of the chr1 unplaced PASS calls (R1.3)

`tables.md` section 7, condensed: per stratum, n: population / germline / other categories.

| run | homopolymer >= 7 | STR | other |
|---|---|---|---|
| legacy COLO829T ONT SNV d9 | 443: 334 / 10 / 99 | 291: 235 / 28 / 28 | 1771: 1323 / 407 / 41 |
| legacy HG008T Illumina INDEL d9 | 411: 369 / 25 / 17 | 114: 103 / 10 / 1 | 9: 6 / 2 / 1 |
| legacy HG008T Illumina SNV d9 | 87: 74 / 12 / 1 | 173: 151 / 20 / 2 | 335: 244 / 90 / 1 |
| legacy HG008T PacBio SNV d9 | 41: 37 / 3 / 1 | 23: 18 / 5 / 0 | 41: 33 / 7 / 1 |
| rcpu COLO829T Illumina SNV d46 | 41: 34 / 5 / 2 | 77: 61 / 13 / 3 | 398: 181 / 204 / 13 |
| rcpu COLO829T fiber-seq SNV d46 | 113: 92 / 2 / 19 | 39: 32 / 3 / 4 | 105: 76 / 22 / 7 |
| rcpu COLO829T ONT SNV d46 (3,000 sampled) | 303: 211 / 19 / 73 | 289: 209 / 67 / 13 | 2408: 1075 / 1276 / 57 |
| rcpu HG008T Illumina SNV d9 | 36: 31 / 5 / 0 | 26: 24 / 2 / 0 | 53: 40 / 12 / 1 |
| rcpu HG008T Illumina SNV d46 | 106: 92 / 12 / 2 | 242: 209 / 31 / 2 | 406: 195 / 207 / 4 |
| rcpu HG008T PacBio SNV d9 | 111: 98 / 6 / 7 | 33: 29 / 4 / 0 | 66: 59 / 5 / 2 |
| rcpu HG008T PacBio SNV d46 | 97: 85 / 6 / 6 | 60: 55 / 5 / 0 | 83: 74 / 8 / 1 |
| rcpu HG008T ONT SNV d9 | 336: 282 / 11 / 43 | 227: 206 / 11 / 10 | 847: 801 / 37 / 9 |
| rcpu HG008T ONT SNV d46 | 387: 322 / 7 / 58 | 219: 185 / 20 / 14 | 652: 425 / 215 / 12 |

### 8. Cross-platform reproducibility of the chr1 unplaced PASS keys (same `node:offset:allele`, same graph)

`tables.md` section 8, condensed: the shared keys by category (calls also PASS on at least one other platform). The
`rcpu` groups have all three platforms (HG008T d9, HG008T d46, COLO829T d46); the legacy d9 group has two.

| group | run | platforms | unplaced PASS keys | also PASS on another platform | shared keys by category |
|---|---|---|---|---|---|
| legacy/d9/HG008T/SNV | legacy HG008T Illumina | 2 | 595 | 50 (8.4 %) | population 36, germline 13, tumor-only 1 |
| legacy/d9/HG008T/SNV | legacy HG008T PacBio | 2 | 105 | 50 (47.6 %) | population 36, germline 13, tumor-only 1 |
| rcpu/d9/HG008T/SNV | rcpu HG008T Illumina d9 | 3 | 115 | 61 (53.0 %) | population 49, germline 11, tumor-only 1 |
| rcpu/d9/HG008T/SNV | rcpu HG008T PacBio d9 | 3 | 210 | 158 (75.2 %) | population 141, germline 13, unverifiable 3, tumor-only 1 |
| rcpu/d9/HG008T/SNV | rcpu HG008T ONT d9 | 3 | 1410 | 169 (12.0 %) | population 147, germline 16, unverifiable 3, unsupported 2, tumor-only 1 |
| rcpu/d46/HG008T/SNV | rcpu HG008T Illumina d46 | 3 | 754 | 163 (21.6 %) | population 134, germline 27, tumor-only 1, unverifiable 1 |
| rcpu/d46/HG008T/SNV | rcpu HG008T PacBio d46 | 3 | 240 | 208 (86.7 %) | population 184, germline 19, unverifiable 4, tumor-only 1 |
| rcpu/d46/HG008T/SNV | rcpu HG008T ONT d46 | 3 | 1258 | 254 (20.2 %) | population 217, germline 29, unverifiable 6, unsupported 1, tumor-only 1 |
| rcpu/d46/COLO829T/SNV | rcpu COLO829T Illumina d46 | 3 | 516 | 187 (36.2 %) | germline 92, population 81, not in normal 13, unverifiable 1 |
| rcpu/d46/COLO829T/SNV | rcpu COLO829T fiber-seq d46 | 3 | 257 | 175 (68.1 %) | population 136, germline 18, not in normal 11, unverifiable 10 |
| rcpu/d46/COLO829T/SNV | rcpu COLO829T ONT d46 | 3 | 4569 | 303 (6.6 %) | of the 3,000 processed: population 109, germline 73, not in normal 14, unverifiable 7 |

The one tumor-only key shared by all three HG008T platforms is, on d9, the gold key of §5 (chr1:78,932,081 A>T, E9a
class i) and, on d46, a non-gold class-ii call at chr1:192,857,712-192,857,718 (homopolymer 7, CC repeat; AF 0.38-0.52;
from `per_call.tsv.gz`).

## Interpretation

1. **The unplaced PASS calls are almost all population or germline alleles (R3.3, R3.4).**
   - The graph PoN tags 49.8 % (COLO829T ONT d46) to 91.4 % (HG008T ONT d9) of them in the 13 runs (§1). By platform
     in the `rcpu` runs: Illumina 53.5-82.6 %, PacBio / fiber-seq 77.8-89.2 %, ONT 49.8-91.4 %.
   - The assembly check alone places 53.3 % (COLO829T fiber-seq) to 96.5 % (HG008T Illumina d9) of them in the
     patient's normal assembly (§1, second table).
   - After the population precedence, the nine `rcpu` runs split as: germline 4.2 % (HG008T ONT d9) to 33.2 % (HG008T
     Illumina d46) on HG008T and 10.5 % (fiber-seq) to 45.4 % (ONT) on COLO829T; unverifiable 0.0-5.8 %; unsupported
     0.0-0.9 % (§1).
   - Assembly-supported tumor-only calls are 0.2-0.9 % of the unplaced PASS calls of the `rcpu` HG008T runs after the
     population precedence (legacy 0.3-1.7 %), and 1.1-5.2 % before it (legacy 2.9-10.9 %).
   - COLO829T `not in normal` (not confirmed somatic): 2.4-5.8 % after the precedence, 4.5 % (Illumina), 6.6 % (ONT,
     sampled) and 14.8 % (fiber-seq) before it.
   - This agrees with the internal assessment in the design (R3.3: graph-coordinate calls are mostly false) and supports
     counting every unconfirmed graph-coordinate call as a false positive in the complete-output metrics (R3.4).
2. **Complete-output precision (R3.4), chr1, no BED, `rcpu` runs** (§2). Placed -> every unplaced call false -> graph
   PoN applied to the unplaced calls:
   - COLO829T d46: Illumina 0.8669 -> 0.69 -> 0.7745, fiber-seq 0.8657 -> 0.767 -> 0.8417, ONT 0.7871 -> 0.2516 ->
     0.3807.
   - HG008T d9: Illumina 0.4427 -> 0.4071 -> 0.4361, PacBio 0.4886 -> 0.4153 -> 0.479, ONT 0.4542 -> 0.2159 -> 0.4149.
   - HG008T d46: Illumina 0.4179 -> 0.2737 -> 0.3541, PacBio 0.4658 -> 0.3912 -> 0.4564, ONT 0.419 -> 0.2236 -> 0.3416.
   - The graph PoN recovers most of the precision loss for the HG008T d9 runs, HG008T PacBio d46 and COLO829T fiber-seq,
     about half of it for HG008T Illumina / ONT d46 and COLO829T Illumina, and little of it for COLO829T ONT.
   - The unplaced share of the PASS calls is 8.0-68.0 % before the graph PoN (the manuscript says "about 2 %") and
     1.5 % (HG008T Illumina d9) to 51.6 % (COLO829T ONT d46) after it.
   - On HG008T, d9 leaves less than d46 after the graph PoN: 1.5 / 2.0 / 8.7 % against 15.3 / 2.0 / 18.5 % (Illumina /
     PacBio / ONT). The d46 remainder is mostly germline class iii (§6). Only Illumina is a matched d9-d46 pair.
   - Counting the tumor-only remainder as true changes precision by at most 0.0025 in the `rcpu` runs (0.0132 for the
     legacy Illumina INDELs).
   - Corrected chr1 F1 of the nine `rcpu` runs (§2; recall is rtg's and does not change):

     | run | F1 placed (nobed) | F1 all unplaced FP (nobed) | F1 graph PoN (nobed) | F1 graph PoN (bed) |
     |---|---|---|---|---|
     | COLO829T Illumina d46 | 0.7832 | 0.7019 | 0.7432 | 0.7432 |
     | COLO829T fiber-seq d46 | 0.7786 | 0.7359 | 0.7687 | 0.7687 |
     | COLO829T ONT d46 | 0.7364 | 0.369 | 0.4911 | 0.4911 |
     | HG008T Illumina d9 | 0.5775 | 0.5464 | 0.5718 | 0.5807 |
     | HG008T PacBio d9 | 0.6145 | 0.5531 | 0.6068 | 0.6148 |
     | HG008T ONT d9 | 0.5862 | 0.3423 | 0.5524 | 0.5683 |
     | HG008T Illumina d46 | 0.5607 | 0.4143 | 0.5002 | 0.5197 |
     | HG008T PacBio d46 | 0.598 | 0.5327 | 0.5901 | 0.5954 |
     | HG008T ONT d46 | 0.5633 | 0.3548 | 0.4889 | 0.4988 |

   - In the BED (E9c's rule: an unplaced call counts when its bubble interval is inside), chr1 F1 after the graph PoN is
     0.7432 / 0.7687 / 0.4911 (COLO829T Illumina / fiber-seq / ONT d46), 0.5807 / 0.6148 / 0.5683 (HG008T d9) and
     0.5197 / 0.5954 / 0.4988 (HG008T d46). The GRCh38-only arm has 0.468 / 0.679 / 0.418 and 0.470 / 0.514 / 0.537
     (E9a §7): only HG008T ONT d46 stays below linear.
3. **The graph PoN is not specific for off-reference somatic variants, so it cannot be used as a filter that keeps
   graph-only somatic calls (R3.3).**
   - On the assembly-verified gold set the main rule tags 74.2-88.3 % of the keys (all chromosomes; §4): 81.6-88.3 % on
     d9 and 74.2-79.4 % on d46. On the held-out chr1 keys it tags 42.9-69.2 % (d9; 3 of 7 to 9 of 13) and 80.0-85.7 %
     (d46; 8 of 10 to 12 of 14). On the tumor-only gold keys (clean positives) it tags 77.6-84.5 % (d9) and 69.1-74.5 %
     (d46), against 58.2-71.0 % and 23.3-25.5 % for the GRCh38 PoN on the same events' GRCh38 records. The calibration on
     chr2-22 finds no grid point with a small loss (lowest 42.5-67.0 %, adaptive flank, AC >= 9).
   - It tags most tumor-only unplaced calls of the HG008T runs: `rcpu` 3 of 4, 20 of 24, 10 of 11, 8 of 9, 13 of 16 and
     19 of 22 (Illumina d9 / d46, PacBio d9 / d46, ONT d9 / d46); legacy 15 of 17, 49 of 58 and 8 of 9.
   - On placed somatic SNV truths it is specific (ctrl_tp 0-0.7 % in all 12 SNV runs), but it tags 35 of 145 (24.1 %)
     placed INDEL truths of HG008T Illumina (homopolymer and STR contexts; 1 of 145 at AC >= 9).
   - The gold events themselves are mostly in homopolymers >= 7 or STRs (165 of 231 on d9 Illumina), and the GRCh38 PoN
     tags 62.1-72.7 % of the d9 gold events' GRCh38 records too (26.1-29.1 % on d46). Many off-reference somatic events
     coincide with alleles seen in other genomes; a population filter on the allele context cannot separate them.
4. **The graph PoN finds common HPRC alleles, not rare population alleles.**
   - It tags only 1.5-16.0 % of the placed SNV calls that the GRCh38 PoN tags (ctrl_pon; 84.7 % for the legacy Illumina
     INDELs), since 89 haplotypes cannot represent the rare alleles of the four GRCh38 PoN VCFs (gnomAD, dbSNP, 1000G,
     CoLoRSdb).
   - That the unplaced calls are still 49.8-91.4 % tagged says they are mostly common alleles. The remainder that the
     assemblies call germline are probably rarer germline alleles, the same class as the placed FPs after the PoN
     (87-146 of 125-150 ctrl_fp are germline in the assemblies).
5. **Graph-only somatic sensitivity is small and not measurable on chr1 alone** (§5).
   - The HG008T SNV runs carry 7-12 gold keys on chr1 (the legacy Illumina INDEL run 3). In the six `rcpu` HG008T runs,
     26 of them are unplaced PASS calls and the graph PoN tags 20 of these 26 (legacy: 10 and 8).
   - One gold key is PASS, untagged and tumor-only on all three d9 platforms; E9a places it as a truth SNV (§5).
6. **Cross-platform reproducibility tracks population / germline status, not somatic status** (§8).
   - Reproduced on another platform: PacBio 75.2 % (d9) / 86.7 % (d46), Illumina 53.0 % / 21.6 %, ONT 12.0 % / 20.2 %
     (HG008T); COLO829T d46 fiber-seq 68.1 %, Illumina 36.2 %, ONT 6.6 %.
   - The shared keys are population (e.g. 184 of 208 HG008T PacBio d46) or germline. One tumor-only key per HG008T graph
     is shared by all three platforms (§8); on COLO829T 11-14 shared keys are `not in normal`.
7. **R1.3 (ONT).** ONT has many more unplaced PASS calls: HG008T 1,410 (d9) / 1,258 (d46) against PacBio 210 / 240 and
   Illumina 115 / 754; COLO829T d46 4,569 (68.0 % of PASS) against Illumina 516 and fiber-seq 257.
   - They are not mostly homopolymer errors: 336 of 1,410 (d9) and 387 of 1,258 (d46) HG008T ONT calls are in
     homopolymers >= 7, a smaller share than PacBio's (111 of 210, 97 of 240). They are population or germline alleles.
   - 80-88 % of HG008T ONT's unplaced keys are not reproduced on another platform (12.0 % / 20.2 % are), against 13-25 %
     for PacBio.
   - These data give no support to "improved resolution"; the extra calls are alleles the PoN cannot reach.
   - COLO829T ONT d46 is the extreme case: after the graph PoN 51.6 % of the output is still unplaced, 1,362 of 3,000
     sampled calls are untagged and germline in COLO829BL (2,493 germline before the population precedence), and 2,030
     of the 3,000 have no consistent anchor pair (E9a class iii; 1,272 of those germline).
   - Pansoma is tumor-only, so the assembly classes are a diagnosis, not a filter.

## Caveats

- **chr1 only, one training seed per model**, and the gold keys on chr1 are few (1-12 per run, 7-15 per gold set); the
  gold calibration uses chr2-22.
- **Run sources and graph pairs.** The legacy runs are the local pre-RIS d9 models (INDEL lr 4e-4); the `rcpu` runs are
  the RIS 7fc7396 models on d9 (HG008T only) and d46. The `rcpu` HG008T d9 / d46 rows are a graph comparison, but only
  Illumina is recipe-matched: the PacBio / ONT d9 checkpoints were trained on 1 GPU (E3 caveat 1). COLO829T has `rcpu`
  runs on d46 only. Legacy and `rcpu` d9 rows of the same set are different models.
- **HG008T PacBio SNV d9 (legacy)** is the 10:36 pilot (run directly, log `DATA/test_pacbio.log`, contexts through the
  `variant_summary` reader); it was not re-run with the 11:05 `e9c_run.py`. Its `calls.tsv` has the same 68 columns as
  the other runs, and the tensor reader matched the `variant_summary` layout on 783 of 783 checked tensors.
- **Sampling.** COLO829T ONT d46 has 4,569 unplaced PASS calls; 3,000 (seed 0) were processed and the precision rows
  scale the sample's fractions to 4,569 (`unplaced left` 2292.1 is therefore an estimate). Every other run was processed
  in full.
- **HPRC v1.1 haplotypes only** (89), also for the d46 runs, because the v2.1 full-graph VCF is not local. The in-graph
  test uses the run's own graph (d46 for the d46 runs). This may leave more d46 calls untagged than a v2.1 haplotype set
  would.
- **The assembly check is context-level.** No indel is allowed in a hit (repeat-length differences never match), NM <= 1,
  50-base contexts. Multi-copy contexts are frequent for COLO829T ONT (1,298 of 2,505 d9, 1,594 of 3,000 d46), COLO829T
  Illumina (196 of 516) and HG008T Illumina / ONT d46 (245 of 754, 279 of 1,258), so their germline class can include
  paralogous copies. On the gold set (somatic by construction) the check calls 9-60 keys per set germline (e.g. 60 of
  231, d9 Illumina): the germline class is not error-free at off-reference loci, while placed somatic truths come out
  tumor-only / not in normal in 144-150 of 145-150 (ctrl_tp).
- **COLO829T has no tumor assembly**: `not_in_normal` is not confirmed somatic, and `+tumor-only TP` equals `graph PoN`.
  The COLO829T BED rows equal the nobed rows in the run's own rtg report.
- **`+tumor-only TP` is an upper bound** (tumor-only calls are not truth-matched).
- **The BED rule for unplaced calls differs from E9a's.** Here an unplaced call counts in the BED when its bubble interval
  lies inside it (no interval: counted); E9a counts every unplaced call. The `all FP` rows therefore differ slightly from
  E9a's `placed + all unplaced as FP` for HG008T.
- **Gold-set purity.** One or two HG008-T gold keys are candidates in the tensors of each COLO829T run (d9 and d46); in
  the ONT and fiber-seq d46 runs one of them is germline in COLO829BL. Together with the GRCh38 PoN tag rate of the d9
  gold events (62.1-72.7 %) this shows that the gold set contains alleles that also occur in other genomes.
- The main rule was fixed before chr1 was looked at and was not re-tuned after the calibration failed; the full grid is
  reported instead.

## Status

### Done

- Edge tables d9 / d46 (387146, 387147).
- E9c + E9d on 13 runs:
  - legacy d9 HG008T Illumina SNV (387265_0), Illumina INDEL (387265_1), PacBio SNV (pilot), COLO829T ONT SNV
    (387265_9);
  - RIS 7fc7396 d46 HG008T PacBio / ONT SNV and COLO829T ONT SNV (387431_0-2);
  - RIS 7fc7396 HG008T Illumina / PacBio / ONT SNV d9, HG008T Illumina SNV d46, COLO829T Illumina / fiber-seq SNV d46
    (405788_0-5).
- Gold-set calibration d9 (387266) and d46 (387271, with the d46 re-alignment of the gold events); GRCh38 PoN on the gold
  truth records (`gold_grch38_pon.py`).
- Assembly checks HG008T (405793, 18,217 contexts) and COLO829T (405794, 18,327 contexts) on every context of the 13
  runs and the gold sets.
- Tables (`collect.py`, 405796, 2026-10-15 15:27). Its inputs have not changed since, so a rerun of `collect.py` would
  reproduce them. No E9c job is queued.

### Still missing (2026-10-15)

| item | state | what it needs |
|---|---|---|
| 18 RIS 7fc7396 INDEL in-sample directories: `vcf_chr1` and `vcf_chr1_valbest` of HG008T Illumina / PacBio / ONT on d9 and d46, COLO829T Illumina / fiber-seq / ONT on d46 | finished, listed by `discover.py`, not processed | an `e9c_runs` array, then `asm.sh HG008T`, `asm.sh COLO829T` and `collect.sh` (`run_all.sh` chains them; Slurm, not run here). E9a projects INDELs too (`../E9a_unplaced_projection/`, also pending). |
| 9 `__norescue` SNV re-calls (E3b) | finished, not processed | the same chain; they carry only 6-162 unplaced PASS calls each (E16 `tables.md` §1a), so they are cheap |
| 8 `vcf_chr1_<SET>` directories of the COLO829T fiber-seq joint d9 model (SNV and INDEL, P100 / P75 / P50 / P25) | finished, not processed | the same chain; the only d9 coverage of COLO829T by a 7fc7396 model |
| 40 cross-sample `xs__` directories (E12a) | finished, not processed | the user's decision whether to project them; `run_all.sh` as written would include them (submit an explicit target list otherwise) |
| COLO829T d9 single-set runs | not on this cluster | the RIS 7fc7396 COLO829T d9 models exist only as the joint model here |
| 8 legacy d9 runs (HG008T PacBio INDEL, ONT SNV, ONT INDEL; COLO829T Illumina SNV / INDEL, fiber-seq SNV / INDEL, ONT INDEL) | 387265_2-_7 cancelled 2026-10-15 13:05 while held; 387265_8 / _10 FAILED on the empty-query bug (fixed in `e9c_lib.py`) | a resubmission; the HG008T SNV ones are superseded by the `rcpu` d9 runs |
| comparators (E8: ClairS-TO, DeepSomatic on GRCh38) | not yet re-scored | they are GRCh38 callers without unplaced calls; they enter only through the complete-output comparison |
| seed repeats | pending GPU (R3.9) | every run here is one seed |
| final E9e chain (E9a placement first, then E9c / E9d on the rest) | the E9a class cross-tab exists for all nine `rcpu` runs (§6) | the combined precision (class i placed by E9a, E9c / E9d on classes ii / iii) is not assembled; §2 is over all unplaced calls |

### Needs (user decisions or data)

- **Use of the graph PoN.** The pre-registered criterion (small loss on the gold set) fails for every rule in the grid.
  The user decides whether the graph PoN is reported only as a classification of the unplaced calls (as here) or used
  as a filter despite the loss of 74-88 % of the gold keys.
- **HPRC v2.1 full-graph VCF** (download; also listed under E7) for the HPRC-haplotype test of the d46 runs.
- **A COLO829T tumor assembly** (optional SMaHT download in the design) to confirm COLO829T `not_in_normal` calls; the
  COLO829BL ONT read evidence of design step E9d-4 is not done.
- CHM13 projection is dropped (user decision 2026-10-13).
