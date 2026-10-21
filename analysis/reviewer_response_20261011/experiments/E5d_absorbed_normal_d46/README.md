# E5d: the matched normal and the tumor at absorbed truth alleles, on HPRC v2.1 (d46) and HPRC v1.1 (d9)

Reviewer point: **R2.3** ("for false negatives on existing graph nodes, check whether the normal's reads align to the
alternative path"), third paragraph ("The normal at missed sites on graph nodes"), so that it can lead with the HPRC v2.1
graph as R1.4 does. The user asked for it on 2026-10-21 ("v2.1的补上").

**Status (2026-10-21): done.** Jobs 450304 (sets), 450312 (minimap2 + s6 parse + s8b), 450314 (tumor reads); `tables.md`,
`summary.json` written by `e4_summary.py` on the login node.

## 1. Question

The R2.3 paragraph quotes, for the HPRC v1.1 (d9) graph, from `analysis/graph_absorbed_somatic_20261001` (HG008T) and
`analysis/graph_absorbed_somatic_colo829t_20261005` (COLO829T):

- truth alleles whose tumor reads align perfectly to an existing graph path: HG008T 2,478, COLO829T 497;
- carried by the matched normal: 110 (4.4%) / 34 (6.8%); absent from both normal haplotypes, on nodes walked by other
  HPRC individuals: 2,085 (84.1%) / 430 (86.5%);
- HG008T, at the 1,935 absent INDELs: the tumor assembly carries a sequence in neither normal haplotype at 1,931 and
  exactly the truth allele at 1,652; tumor Illumina reads show the expected allele at 96.2% (95.7% at random non-absorbed
  truth INDELs);
- the GRCh38-only arm calls at most one of the 38-114 HG008T truth INDELs that v1.1 absorbs (chr1, BED; E5b);
- COLO829T: up to 15.8% of the absorbed INDELs may be germline-like.

E5d gives the same quantities on d46, and on d9 with the same set rule.

## 2. Data and set definition

- **Absorbed set (both graphs, one rule): the E5b set that R1.4 quotes.** A truth allele is absorbed on graph X when, on
  at least one platform, it has no candidate on X (`no_candidate`), a tensor on the GRCh38-only arm, and a path of X
  spells GRCh38 + the truth allele between the window anchors (E5b s4, exact spelling). Statuses and spelling are read
  from `E5b_absorbed_d46` (`per_truth.tsv.gz`); `e1_sets.py` re-derives the spelling path with E5b's own `spells()` and
  asserts that every set allele is spelled again (all were). Platforms: HG008T Illumina / PacBio / ONT, COLO829T
  Illumina / fiberseq / ONT; chr1-22, all truth records.
- **This is not the d9 set of the 2026-10-01 / 10-05 analyses.** That set is read-level: the majority of the ALT-like
  reads leave GRCh38 through a graph branch or skip edge with no edit (re-decoded tensor-set reads, 2026-09-30 sets). It
  is not available on d46 without re-decoding every d46 GAM at every locus (in 2026-10-01 PacBio alone took 7-16
  worker-minutes per locus), so d46 uses the read-independent rule, and d9 is recomputed with the same rule for the
  side-by-side table. Overlap with the read-level d9 sets: HG008T d46 1,258 of 1,607 (SNV 98 of 140, INDEL 1,160 of
  1,467); COLO829T d46 325 of 401; d9 (same rule) 1,826 of 2,002 and 427 of 446. E5b §5.4 showed that the spelling rule
  recovers 80.8% of the read-level HG008T INDELs; the rest are the 'closest' / 'with_germline' classes, which the exact
  spelling cannot see.
- **Set sizes** (`$DATA/sets.tsv`; equal to E5b `tables.md` §3 "absorbed on >= 1 platform"):

| sample | d9 SNV / INDEL | d46 SNV / INDEL |
|---|---|---|
| HG008T | 184 / 1,818 | 140 / 1,467 |
| COLO829T | 64 / 382 | 54 / 347 |

- Normal and tumor data are those of the two earlier analyses: HG008-N v6.2 assembly (both haplotypes) and dipcall;
  HG008-T v3.2 assembly; HG008-T GRCh38 BAMs (Illumina 195x, PacBio HiFi 116x); COLO829BL verkko 2.1 hapX / hapY and its
  dipcall (`/scratch/qfu/COLO829BL_DSA`, read only).
- Large data: `/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46/` (`$DATA`).

## 3. Method

1. **e1 `e1_sets.py`** (job 450304, 8 min, 7.5 GB): the sets above, with the node path of each ALT; writes
   `$DATA/sets.tsv`, `colo829t_d46_absorbed_loci.tsv` (401 loci, for E20) and `hg008t_d46_absorbed_loci.tsv` (1,607).
2. **e2 HG008T allele calls for the 438 set alleles outside the 2,478 of 2026-10-01.** `e2_hg008t_loci.py` writes their
   loci in that analysis' s1 format: the window is E5b's d9 window, which equals s1's window at all 2,478 old loci (first /
   last node, window and event window identical); asm_event = the truth INFO HG008Nv62SOMATICVARIANT. Then the
   2026-10-01 steps unchanged except for the data directory: `e2_s6_assembly_seq.py` (minimap2 2.28 sr of REF / ALT
   contexts against both assemblies) and `e2_s8b_array_alleles.py` (whole tandem array between common unique anchors,
   assembly and dipcall per normal haplotype, tumor contig verdict) (job 450312, 12 min).
3. **e3 `e3_hg008t_tumor_reads.py`** (job 450314, 29 min): the biology audit's `b2_tumor_reads.py` `work()`, imported, on
   the 438 new loci (linear GRCh38 tumor BAMs, <= 300 reads with MAPQ >= 20, repeat window); then the audit's b4 rule on
   old + new rows (expected tumor allele = normal event-hap length + the normal-frame event length; 'seen' at >= max(3,
   15%) of the reads; loci with >= 10 reads and a length-changing event). Controls = the audit's 400 + 400 random
   non-absorbed truths (reused, not redrawn).
4. **e4 `e4_summary.py`**: the classification and tables.
   - HG008T, per allele: the s9 hap call (s8b array call; assembly and dipcall disagreeing on ALT -> unresolved; both NA
     -> the s6 minimap2 call, without s9's giraffe walk, which only confirms a 'no'); pattern by the truth INFO event hap;
     patient frame as s9. **present** = a normal hap carries the ALT (other hap only, no event hap, or patient frame; not
     when the INFO event gives a novel tumor allele); **ambiguous** = ALT on both haps or on the event hap only, the
     INFO-novel conflict, or a hap unresolved; **absent** = both haps resolved, neither carries the ALT. s9's rules that
     need reads or HPRC walks ('closest', 'with_germline', no somatic element, no complete HPRC traversal) are dropped:
     every set allele is exactly spelled.
   - COLO829T, per allele: c3 `normal_frame.tsv` (every truth INDEL and the 67 SNVs of 2026-10-05) hap alleles; the 4
     d46 / 2 d9 SNVs outside it: COLO829BL dipcall at the site (in dip.bed, PASS, GT hapY|hapX). present / absent /
     ambiguous as for HG008T; germline-like = a hap has exactly the ALT length with other bases and no hap carries the ALT
     (c7 `germline_like_alt_len`).
   - "On nodes walked by other HPRC individuals": both graphs keep only nodes and edges supported by >= 10% of their
     haplotypes (90 / 464), and neither contains HG008 or COLO829, so every spelled ALT path is made of population
     nodes and edges. `alt_min_paths` (E5b) gives the fewest haplotypes on a non-GRCh38 node of the path.
   - **Check** (tables.md §4): on the read-level d9 sets these rules give present 110 of the published 110 HG008T and 33
     of 34 COLO829T (the missing one, 2562, is present only through its `with_germline` path allele), and absent for all
     2,085 / 430 published absent loci; the published ambiguous loci become mostly absent here (264 / 31), because their
     ambiguity came from the dropped graph / read rules.

## 4. Results (`tables.md`, `summary.json`)

### 4.1 The matched normal

| sample | graph | set | absorbed | normal carries the ALT | absent from both normal haplotypes | ambiguous |
|---|---|---|---:|---|---|---|
| HG008T | v1.1 (d9) | read-level, 2026-10-01 | 2,478 | 110 (4.4%) | 2,085 (84.1%) | 283 (11.4%) |
| HG008T | v1.1 (d9) | spelled (this rule) | 2,002 | 104 (5.2%) | 1,881 (94.0%) | 17 (0.8%) |
| HG008T | **v2.1 (d46)** | spelled | **1,607** | **110 (6.8%)** | **1,481 (92.2%)** | 16 (1.0%) |
| COLO829T | v1.1 (d9) | read-level, 2026-10-05 | 497 | 34 (6.8%) | 430 (86.5%) | 33 (6.6%) |
| COLO829T | v1.1 (d9) | spelled | 446 | 26 (5.8%) | 419 (93.9%) | 1 (0.2%) |
| COLO829T | **v2.1 (d46)** | spelled | **401** | **27 (6.7%)** | **371 (92.5%)** | 3 (0.7%) |

- Per kind, d46: HG008T SNV 6 present / 130 absent / 4 ambiguous, INDEL 104 / 1,351 / 12; COLO829T SNV 1 / 52 / 1,
  INDEL 26 / 319 / 2.
- PoN-tagged (truth tags of E5b): 97.0% of the d46 sets (HG008T), 98.5% (COLO829T); d9 97.3% / 98.7%.
- The ALT paths' non-GRCh38 nodes are walked by a median of 138 (HG008T) and 110 (COLO829T) of the 464 d46 haplotypes
  (d9: 27 and 23.5 of 90); 634 / 157 d46 ALT paths use skip edges only.

### 4.2 HG008T: tumor assembly and tumor reads at the absorbed INDELs absent from the normal

| graph | set | absent INDELs | tumor contig in neither normal hap | tumor contig = GRCh38 + truth | Illumina reads show the expected allele | PacBio reads |
|---|---|---:|---|---|---|---|
| v1.1 (d9) | read-level, 2026-10-01 | 1,935 | 1,931 | 1,652 | quoted 96.2%; recomputed from the audit's b4 table 1,804 / 1,890 (95.4%) | 1,885 / 1,945 (96.9%) |
| v1.1 (d9) | spelled | 1,706 | 1,702 (99.8%) | 1,478 (86.6%) | 1,586 (97.8%) | 1,605 (98.0%) |
| **v2.1 (d46)** | spelled | **1,351** | **1,347 (99.7%)** | **1,153 (85.3%)** | **1,240 (97.4%)** | 1,260 (97.5%) |
| random non-absorbed truth INDELs (audit controls) | | | | | 353 / 369 (95.7%) | 372 / 386 (96.4%) |

Read-percentage denominators are the loci with >= 10 reads and a length-changing normal-frame event.

### 4.3 GRCh38-only arm at the absorbed HG008T INDELs (chr1, benchmark BED, after the PoN; E5b `linear_recovery.md`, new INDEL rule)

| platform | absorbed on d9: n / GRCh38-only TP (curve / PASS) | absorbed on d46: n / GRCh38-only TP (curve / PASS) |
|---|---|---|
| Illumina | 114 / 0 / 0 | 93 / 2 / 2 |
| PacBio | 64 / 1 / 1 | 66 / 2 / 2 |
| ONT | 38 / 1 / 1 | 39 / 2 / 2 |

On v2.1 the GRCh38-only arm calls at most 2 of the 39-93 absorbed truth INDELs per platform (v1.1: at most 1 of 38-114).
Not recomputed here; the d46 rows already existed in E5b.

### 4.4 COLO829T: germline-like absorbed INDELs

| graph | set | absorbed INDELs | normal carries the ALT | a normal hap has the ALT length, other bases | up to |
|---|---|---:|---:|---:|---|
| v1.1 (d9) | read-level, 2026-10-05 | 430 | 33 | 35 | 68 (15.8%) |
| v1.1 (d9) | spelled | 382 | 25 | 24 | 49 (12.8%) |
| **v2.1 (d46)** | spelled | **347** | **26** | **22** | **48 (13.8%)** |

## 5. Reading

- On v2.1, as on v1.1, the absorbed alleles are mostly not the patient's: 92.2% (HG008T) and 92.5% (COLO829T) are absent
  from both normal haplotypes and lie on population nodes; the normal carries the ALT at 6.8% / 6.7%.
- At the HG008T INDELs absent from the normal, the tumor assembly differs from both normal haplotypes at 99.7% and equals
  the truth allele at 85.3%; tumor reads show the expected allele as often as at random non-absorbed truth INDELs (97.4%
  vs 95.7%, Illumina). These are real somatic events that the graph hides.
- The spelled rule gives a cleaner split than the read-level rule (fewer ambiguous loci) because it keeps only exactly
  spelled alleles; the shares absent / present are of the same order on both graphs.

## 6. Caveats

1. **Set rule.** d46 and the side-by-side d9 use the read-independent spelled rule (E5b, R1.4), not the read-level rule of
   the paragraph's current d9 numbers; the two d9 rows of each table show the difference. A strict test: alleles needing
   nearby germline alleles or spelled only approximately are not in the sets (E5b Caveat 1), and the d46 subgraph has no
   completeness check (E5b Caveat 2).
2. **Reduced classification.** The s9 / c7 rules that need re-decoded reads or HPRC walks are not applied (section 3,
   item 4); no HPRC exact-allele frequencies or populations were computed on d46 (the d46 GFA is 375 GB and there is no
   d46 deconstruct VCF). "On population nodes" rests on the d46 frequency filter; by the pipeline graph index, 49 of 973
   HG008T and 13 of 244 COLO829T d46 ALT paths with a branch node have a node below 46 haplotypes (fewest 2 / 6): the
   index counts haplotype path fragments (d46 cuts every haplotype into ~0.3 M pieces), so these counts are approximate.
3. **The published 96.2%.** Recomputed from the audit's own `b4_tumor_alleles.tsv` with its rule, the d9 read-level absent
   INDELs give 1,804 / 1,890 = 95.4% (Illumina); the control 95.7% reproduces exactly. The 96.2% in the 2026-10-01 README
   could not be traced to a file; the R2.3 sentence should use 95.4% or the E5d numbers.
4. **HG008T normal calls of the 438 new alleles** use the s6 minimap2 call (no s7 giraffe walk) when the array call is NA;
   at the 2,478 old alleles this changes no present / absent call (tables.md §4).
5. COLO829T has no tumor assembly; its read support was not recomputed (the 2026-10-05 analysis did not have it either).

## 7. Files

| file | content |
|---|---|
| `e1_sets.py` | absorbed sets on d9 / d46 with the ALT node path -> `$DATA/sets.tsv`, `*_d46_absorbed_loci.tsv` |
| `colo829t_d46_absorbed_loci.tsv` | 401 COLO829T d46 absorbed loci (chrom, pos, ref, alt, kind, truth_id, platforms, PoN tag, elements, node path `<id><+/->`) |
| `hg008t_d46_absorbed_loci.tsv` | the same for HG008T (1,607) |
| `e2_hg008t_loci.py` | s1-format loci of the 438 new HG008T alleles -> `$DATA/hg008t/{variant_set,loci}.tsv` |
| `e2_s6_assembly_seq.py`, `e2_s8b_array_alleles.py` | copies of the 2026-10-01 s6 / s8b with only the data directory changed |
| `e3_hg008t_tumor_reads.py` | audit b2 on the new loci + b4 rule on all -> `$DATA/hg008t/{b2_new,b4_all}.tsv` |
| `e4_summary.py` | classification, `tables.md`, `summary.json`, `$DATA/per_allele.tsv` |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261021/E5d/` | Slurm scripts and logs (`e1_sets`, `e2_s6_minimap2`, `e3_tumor_reads`) |
