# HG008T somatic truths absorbed by the HPRC graph: HG008-N germline, or other people's alleles?

2026-10-01. Pansoma makes a candidate (and a tensor) only from the edits a read has relative to the graph nodes it aligns
to. When a somatic allele is already a path of the HPRC v1.1 d9 graph, its reads follow that branch or skip edge with no
edit, and the truth gets no tensor. This analysis takes every HG008T somatic truth allele lost this way and asks:

1. how many SNVs and INDELs, per platform;
2. which graph nodes / edges the reads align to;
3. whether those nodes are the patient's own germline (HG008-N hap1 / hap2), or alleles of other HPRC individuals and
   populations that coincide with the somatic site.

## Short answer

- **Count.** The ALT reads align perfectly (no edit) to a graph branch or skip edge for:
  - PacBio: 170 SNVs / 1,363 INDELs; ONT: 148 / 989; Illumina: 122 / 1,916.
  - That is 36-43% of each platform's no-candidate SNVs and 40-48% of its no-candidate INDELs. For INDELs this is exactly
    class I1 of the earlier miss analysis.
  - Union of the three platforms: **201 SNVs + 2,277 INDELs = 2,478 truth alleles**; chr1 inside the BED: 17 + 169.
  - For 162 SNVs / 2,090 INDELs, a d9 path spells the truth allele itself. For the other 39 / 187 ('closest'), the reads
    that align perfectly follow another allele.
- **Nodes.** `per_node.tsv` has one row per (truth, graph element): 3,368 rows on 2,426 loci, using 4,409 distinct
  non-GRCh38 nodes and 1,025 skip edges.
  - INDEL: mostly an insertion branch (the inserted bases as non-GRCh38 nodes between two adjacent GRCh38 nodes) or a
    deletion skip edge.
  - SNV: mostly a same-length SNP branch node.
- **Mostly not germline.** For 2,085 of the 2,478 (84.1%), neither HG008-N haplotype carries the ALT allele, but HPRC
  haplotypes walk the same nodes. This is **pangenome-induced interference**: the graph hides a somatic change.
  - At 1,854 of them, at least one HPRC haplotype carries exactly the same allele. Mostly a one-unit slippage in a long
    homopolymer / STR reproduces a population length allele (median HPRC frequency 0.17) that the d9 graph keeps.
  - At 231 only the nodes are population nodes; no HPRC haplotype has the exact allele.
- **Germline filtering is rare.** Only 110 (4.4%) are an HG008-N germline allele: 5 SNVs, 105 INDELs.
  - In 109 of the 110, the tumor allele is the germline allele of the *other* normal haplotype. The tumor's event
    haplotype changed to it at the site (local conversion or recurrent slippage).
  - A tumor-only caller sees a germline allele there, so absorbing it is germline filtering.
- **Ambiguous.** 283 (11.4%) are unresolved, mainly the 226 'closest' loci.
- **The graph decides for INDELs.** The repo's linear PoN rule would remove the absent-from-HG008-N SNVs anyway (100%
  tagged, vs 8.4% of the other chr1-22 truth SNVs). The repo INDEL path has no PoN, so the INDEL losses come from the
  graph alone; even at a population-AF floor of 0.05, 30% would be kept by a linear caller.

## Definitions

| Term | Meaning |
|---|---|
| d9 graph | HPRC v1.1 Minigraph-Cactus GRCh38 graph, frequency-filtered: nodes on fewer than ~9 of the 90 haplotypes removed (>= 7 HPRC haplotypes per non-reference node in the final GFA). 44 samples x 2 haplotypes + CHM13 + GRCh38. The tensors were built on it |
| GRCh38 node | a node the GRCh38 path visits once, forward (`classify.py ref()`) |
| perfect bypass | the set definition: status `no_candidate` on the platform, and the majority of the ALT-like reads (reads.py: closer to ALT than to REF) leave GRCh38 through a non-GRCh38 node run or an edge that skips GRCh38 bases, with no edit within ±5 bp. Read-level reason `site_bypassed:branch:no_edit` / `site_bypassed:skip_edge:no_edit` |
| window | the GRCh38 nodes overlapping truth ±60 bp, widened to the event window ±20 bp; its first and last node are the **anchors** |
| event window | (s1) union of the periodic stretches (period 1-6, length >= max(2p, p+3)) touching the event, ±5 bp |
| tandem array | (s8b, used for allele calls) the event window grown by every periodic stretch of period 1-60 touching it, including adjacent arrays (up to 3 kb), bracketed by GRCh38 24-mers that are unique and present once in every placed contig |
| ALT path | a d9 path from anchor to anchor that spells GRCh38 + the truth allele: `exact`; or that needs nearby HG008-N germline alleles too: `with_germline`; or none exists: `closest` |
| element | one departure of the ALT path from GRCh38. `B:x:y:>n1>n2` = a run of non-GRCh38 nodes between GRCh38 nodes x and y; `S:x:y` = an edge x->y that skips GRCh38 bases |
| complete traversal | a haplotype walk that visits both anchors in order; only these count in HPRC frequencies |
| exact allele | a haplotype's sequence over the tandem array (between the nearest GRCh38 nodes it shares with the ALT path) equals the ALT allele; for with_germline, ALT + the patient's germline alleles |
| event hap | the HG008-N v6.2 haplotype named by the truth INFO `HG008Nv62SOMATICVARIANT` (`chrN_hapK:pos-REF-ALT`, 0-based): the haplotype the somatic event is on |
| patient frame | the coordinates and sequence of the event hap; the patient-frame event is that INFO event |
| truth-spelling read | a no-edit read whose walk, spelled over its GRCh38 anchor span, equals GRCh38 + truth (or + the patient's germline alleles). A germline-only read spells only HG008-N germline alleles |
| perfect platform | a platform on which the truth is perfect bypass |
| PoN | population panel of normals: the repo rule (`scripts/filter_panel_of_normals.py`: allele match; gnomAD / CoLoRSdb AF >= 1e-4; dbSNP non-somatic; 1000G) |

## Data

| Item | Source |
|---|---|
| graph | `/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.*` (GBZ, 49.6 GB GFA, giraffe indexes, `vg deconstruct` VCF) |
| truth | GIAB `HG008-T_somatic_smvar_benchmark_v0.2` chr1-22 alleles, BED `_all.bed` (168 of the 2,478 lie outside it; the chr1 column is inside); `nogermlineoverlap.bed` as a stricter scope |
| no-candidate statuses | `analysis/tensor_recall_20260930/per_truth_HG008T.tsv` (tensor sets of 2026-09-30, now in `pansoma_v2_tensors/backup_ch6_linear100_20261001`) |
| read-level reasons and miss classes | `analysis/hg008_somatic_miss_20260926/{pacbio_v6,ont,illumina}/{snv,indel}_no_candidate.tsv` |
| HG008-N | v6.2 curated assembly `HG008N_curatedv6_250714_bothhaps_polished6.2.fasta.gz` (`chrN_hap1` / `chrN_hap2`); dipcall `HG008N_GRCh38_dipcall.dip.vcf.gz` (phased hap1\|hap2, checked) + `.dip.bed` |
| HG008-T | v3.2 assembly `HG008T_v3.2.fasta.gz` |
| HPRC genotypes | d9 VCF and the full-graph `vg deconstruct` VCF (`hprc-v1.1-mc-grch38.raw.vcf.gz`), 44 sample columns + CHM13 |
| populations | HPRC Year 1 sample metadata, cross-checked with IGSR (42/44) and Coriell; NA21309 (Maasai, MKK, AFR) only from the HPRC note. HPRC v1.1 = AFR 23 / AMR 16 / EAS 4 / SAS 1 / **EUR 0** samples |
| external AF | gnomAD / CoLoRSdb AF per truth from `analysis/tensor_recall_20260930/pon/truth_pon_HG008T.tsv` |

- The HG008 donor is of European genetic ancestry (McDaniel et al. 2025, *Sci Data* 12:1195,
  doi 10.1038/s41597-025-05438-2).
- The 2026-10-01 tensor rebuild (channel 6, build AF 0.08) cannot change this set: a read with no edit gives no candidate
  in any build.

## 1. How many truths are lost because the reads align perfectly to a graph branch

| Platform | SNV no candidate | **SNV perfect bypass** | of which a d9 path spells the truth | INDEL no candidate | **INDEL perfect bypass** | of which a d9 path spells the truth | chr1 BED SNV / INDEL perfect bypass |
|---|---:|---:|---:|---:|---:|---:|---:|
| PacBio HiFi | 398 | **170** (42.7%) | 140 | 2,976 | **1,363** (45.8%) | 1,252 | 14 / 102 |
| ONT-UL | 342 | **148** (43.3%) | 128 | 2,454 | **989** (40.3%) | 876 | 14 / 74 |
| Illumina | 339 | **122** (36.0%) | 105 | 3,954 | **1,916** (48.5%) | 1,807 | 12 / 136 |
| union | | **201** | 162 | | **2,277** | 2,090 | 17 / 169 |

- **Share of all truth.** The perfect-bypass INDELs are 11.6-22.6% of all 8,496 chr1-22 truth INDEL alleles; the SNVs
  are 1.4-2.0% of all 8,690.
- **Platforms.** 704 INDELs and 93 SNVs are perfect bypass on all three platforms.
- **SNV miss classes** (PacBio / ONT / Illumina):
  - S3 (ALT is an existing graph SNP allele): 55 / 44 / 61;
  - S4b (repeat, absorbed by a graph path): 89 / 84 / 45;
  - S4a (repeat, the SNV becomes an edit on another branch): 22 / 18 / 11;
  - S2 (ALT not seen in same-length reads): 4 / 2 / 4;
  - S1 (low MAPQ): 0 / 0 / 1.

  Set membership follows the read-level reason, not the class. The 21 SNVs whose classes are only S2 / S4a are in the set
  because the majority of their ALT-like reads bypass without edit (16 end up absent, 5 ambiguous).
- **Does a d9 path spell the truth ALT?** Checked by sequence (`s2_graph_paths.py`):
  - INDEL: exactly 2,017; with HG008-N germline alleles applied too 73; no path, only a closest allele 187.
  - SNV: 147 / 15 / 39.
  - At the 226 'closest' loci the reads that align perfectly are not spelling the truth.
- **Tumor reads that spell the truth.** The reads were re-decoded (`s3_read_paths.py`): Illumina and PacBio all 2,478 loci,
  ONT all 1,137 ONT-perfect loci (1,365 loci in total).
  - At least one perfect-platform read spells the truth (or truth + germline) at 1,985 INDEL and 113 SNV loci.
  - With the platform's own reads: PacBio 1,174 / 93, ONT 797 / 85, Illumina 1,724 / 78 (INDEL / SNV).
- **Caveat on the read-level class.** In long repeats reads.py's ALT-like reads can be the patient's own germline
  STR-length reads. The most common read walk spells germline only at 184 PacBio-, 190 ONT- and 160 Illumina-perfect
  INDEL loci. So the read-level class slightly over-counts truths whose *tumor* reads were absorbed. Per-locus counts are
  `{P}_n_truth_reads` and `{P}_n_germline_only_reads` in `per_variant.tsv`.

### Repeat context of the absorbed INDELs (`s14_indel_repeat_context.py`, `indel_repeat_context.md` / `.tsv`)

Each INDEL truth is classed in the **patient frame**: its truth INFO `HG008Nv62SOMATICVARIANT` event, in the sequence of the
HG008-N v6.2 event haplotype. This is the frame that matters, because the GRCh38 record of a repeat INDEL depends on the
GRCh38 repeat length. Definitions:
- unit u = the minimal period of the inserted / deleted bases;
- tract L = the longest exact period-u stretch touching the event;
- `HP>=7` = homopolymer, u = 1 and L >= 7;
- `STR` = u 2-6 with >= 3 copies;
- other repeat = `HP4-6`, `VNTR>6` (unit > 6, repeated next to the event) or `in_STR` (the unit does not repeat, but the
  event sits inside a period 1-6 stretch >= 10 bp);
- `SNV in repeat` = a GRCh38 INDEL record whose patient-frame event is a substitution inside a repeat (the GRCh38 record
  absorbs a germline length difference);
- `none` = not in a repeat; `no event` = the truth has no INFO event.

| patient frame | PacBio | ONT | Illumina | union |
|---|---:|---:|---:|---:|
| homopolymer `HP>=7` | 1,003 (73.6%) | 662 (66.9%) | 1,614 (84.2%) | 1,853 (81.4%) |
| `STR` (period 2-6, >= 3 copies) | 248 (18.2%) | 221 (22.3%) | 220 (11.5%) | 292 (12.8%) |
| other repeat (`HP4-6`, `VNTR>6`, `in_STR`) | 6 | 6 | 2 | 6 |
| `SNV in repeat` (HP or STR) | 85 (6.2%) | 82 (8.3%) | 69 (3.6%) | 104 (4.6%) |
| non-repeat (INDEL or SNV) | 1 | 1 | 0 | 1 |
| no INFO event | 20 | 17 | 11 | 21 |
| total | 1,363 | 989 | 1,916 | 2,277 |

- **Size of the change.** 92.4% of the union (2,103) change a repeat by exactly one unit (one base of a homopolymer, one
  copy of an STR unit); 46 change it by 2 or more units.
- **Homopolymer length.** The homopolymers are long: 7-9 bp 21, 10-14 bp 377, 15-19 bp 519, 20-29 bp 806, >= 30 bp 130.
- **GRCh38 frame.** 202 of the union (8.9%) look non-repeat in GRCh38, but in the patient frame they are repeat events
  (77 HP, 66 STR, 47 SNV in repeat, 3 other, 9 no event). The patient's germline repeat length or an interruption hides
  the repeat in GRCh38.
- **Background.** Of all chr1-22 truth INDELs, the share absorbed is:
  - patient frame: `HP>=7` 1,853 / 6,722 (27.6%), `STR` 292 / 795 (36.7%), `none` 0 / 415;
  - GRCh38 frame: 1,688 / 6,360 HP, 369 / 936 STR, 202 / 952 'none'.

## 2. Which nodes the reads align to

`s2_graph_paths.py` enumerates the d9 paths between the window anchors that spell GRCh38 + truth, and lists their
elements. Locus counts (loci with at least one element of the subtype):

| element subtype | meaning | absent-from-HG008-N loci, INDEL / SNV | present loci, INDEL / SNV |
|---|---|---:|---:|
| `ins_branch` | a run of non-GRCh38 nodes between two adjacent GRCh38 nodes (the inserted sequence) | 1,202 / 16 | 48 / 0 |
| `del_skip` | an edge between two GRCh38 nodes that skips the deleted bases | 759 / 17 | 57 / 0 |
| `snv_branch` | a same-length non-GRCh38 node (SNP bubble) | 146 / 133 | 3 / 5 |
| `replacing_branch` | a branch replacing GRCh38 bases with another length | 74 / 4 | 3 / 0 |
| `mnv_branch` | a same-length branch of >= 2 bases | 12 / 6 | 1 / 0 |

- **Elements per locus.** 1,960 loci have one element and 466 have 2-14. Several elements, or an `snv_branch` at an INDEL
  locus, appear in compound repeats or where the graph spells the INDEL as a substitution plus a shorter
  insertion/deletion. 52 'closest' loci have none: their reads' walk departs from GRCh38 only through germline-alone
  elements.
- **What each element is alone** (`alone` in `per_node.tsv`: GRCh38 with only that element applied). At absent loci:
  1,760 elements spell the truth alone, 756 are one piece of a multi-element path, 134 spell a germline allele of the
  patient alone, and 14 spell truth + germline.
- **Path choice.** 181 loci have several ALT paths. At 134 the path walked by the most truth-spelling reads (all platforms
  pooled) is used. At the other 47 no truth read decides: s2's primary path (fewest elements) is kept at 45, a
  phase-consistent alternative at 2.
- **Read check.** Truth-spelling reads contain the whole element set (median fraction 1.0 on every platform). The reads'
  majority walk equals the graph elements at 1,084 / 1,174 PacBio, 724 / 797 ONT and 1,642 / 1,724 Illumina INDEL loci
  with truth reads (tables.md, "Read-vs-graph").
- **Examples.**
  - Present: chr8:76047168 CA>C (truth 7917, `HG008N_present_broad`). The A16 homopolymer becomes A15; the graph puts the
    skip edge `S:54543053:54543055` at the 3' end (`grch38_start0` 76047183). 51 of 86 HPRC haplotypes carry it, and
    HG008-N hap2 has A15.
  - Absent: chr7:17609297 CAAAAA>C (truth 6624, `HG008N_absent_HPRC_other`). A33 becomes A28 through skip edge
    `S:51443495:51443488`.
    - HG008-N has A27 / A34; the somatic +1 on hap1 makes A28.
    - 36 of 76 HPRC haplotypes carry A28.

## 3. HG008-N haplotype membership

**Allele level (the categories use this; `s8b_array_alleles.py`).**
- A haplotype **carries the ALT** if its sequence over the tandem array is identical to GRCh38 + truth, or identical
  apart from that haplotype's own dipcall variants outside the array.
- The assembly sequence is used first; the phased dipcall reconstruction is used where the assembly cannot be anchored.

| HG008-N hap1 / hap2 vs GRCh38 + truth over the array (union, chr1-22) | SNV | INDEL |
|---|---:|---:|
| one hap carries the ALT | 3 (all present) | 105 (98 present; 7 ambiguous: 6 where the truth INFO event gives a tumor allele neither hap has, 1 on the event hap) |
| both haps carry it | 0 | 1 (ambiguous: homozygous germline, truth conflict) |
| neither carries it | 195 (150 absent, 43 ambiguous, 2 present via the patient frame) | 2,164 (1,935 absent, 222 ambiguous, 7 present via the patient frame) |
| a hap unresolved | 3 | 7 |

- **Present via the patient frame** (9 loci, sub_class `other_hap_only:patient_frame`): no haplotype equals GRCh38 +
  truth exactly, because both haplotypes share other germline differences inside the array. The truth INFO event
  applied to the event hap gives exactly the other hap's array sequence, so the tumor allele is still the patient's
  other-hap germline allele.
- **Present loci.**
  - In all 109 with an event hap, the carrying hap is the *non-event* hap. A carrier on the event hap itself appears
    only once: truth 10740, an evidence conflict in `other_ambiguous`.
  - Dipcall reports the identical germline allele at 110/110, and every tumor contig equals the carrying normal hap at
    101/110.
  - Flanking phased SNPs show the event hap is retained in the tumor: 76 both retained, 29 other hap lost, 0 event hap
    lost (`mechanism.md`). So this is a local change of the event hap to the other hap's allele, not regional LOH.
- **Absent loci, the normal's own allele** (`HG008N_allele_class`, INDEL / SNV):
  - germline STR allele of another length on >= 1 hap: 1,179 / 91;
  - both haps = GRCh38: 656 / 21;
  - SNV whose site base is deleted on the event hap: 0 / 30;
  - other: 100 / 8.
- **Absent loci, the tumor assembly.** Its copy over the array equals GRCh38 + truth at 1,697 (`novel_ALT`). At 383 it is
  another sequence found in neither normal hap (`novel_other`; 356 of them also give a novel allele from the INFO event).
  At those 383 the tumor's own array sequence differs from GRCh38 + truth by germline differences, so the HPRC
  exact-allele carriers hold the GRCh38-frame ALT, not exactly the tumor haplotype.

**Node level (giraffe walks of the assembly windows on d9, `s7`).**
- `per_node.tsv` `HG008N_hap1/2_giraffe_walk_contains` and `per_variant.tsv` `HG008N_hap1/2_walks_elements` say whether
  each normal haplotype's own walk through the graph uses the node / edge. This answers "is the node present in
  hap1 / hap2". Walking a node is not carrying the allele.
- At present loci 113 of 119 elements are walked by a normal hap.
- At absent loci 753 of 2,664 elements are walked by a normal hap. At 253 absent loci (178 INDEL / 75 SNV) a normal hap
  walks *all* the somatic elements but has another allele: 212 germline STR alleles of another length, 21 SNVs whose base
  is deleted, 20 other. A germline repeat of another length often uses the same skip edge or insertion branch plus
  another element or edit. So these nodes are partly the patient's, but the somatic allele is not.

## 4. HPRC paths, individuals, populations and frequency

`s4_hprc_walks.py` scans the d9 GFA once (P and W lines, 24 processes, 5 min) and extracts every haplotype's walk between
the window anchors. A median of 76 of the 88 HPRC haplotypes per locus have a complete traversal. Membership at three
levels:

| level | carrier = | absent-from-HG008-N median frequency, SNV / INDEL | present median |
|---|---|---:|---:|
| **exact allele** (categories use this) | walk spells the exact allele over the tandem array | 0.000 / 0.167 (0.105 / 0.177 at loci with >= 1 carrier) | 0.315 |
| element set | walk contains all elements of any ALT path | 0.322 / 0.193 | |
| element (min) | walk contains each element (minimum over the elements) | 0.320 / 0.198 | |

- **Absent loci by HPRC support** (`sub_class`):
  - `exact_allele`, at least one HPRC haplotype has exactly this allele: 65 SNV / 1,789 INDEL.
    - Carried by a median of 11 individuals.
    - Seen in 2-4 superpopulations at 1,672 of these 1,854 loci; in a single individual at 65.
  - `element_set_other_allele`, HPRC haplotypes walk all the somatic nodes but spell another allele (mostly another
    repeat length): 76 / 101.
  - `recombinant_pieces`, each element is walked by some haplotype but no haplotype walks a whole ALT path: 9 / 45.
  - Absent SNVs sit in repeats (section 5), so the exact array allele rarely recurs: 85 of 150 have no exact HPRC carrier.
- **"Rare".** Rare means an exact-allele frequency below 0.20 among complete traversals: 0-17 carrier haplotypes, and 4
  present_rare loci have no exact HPRC carrier. Separately, every somatic branch node not on a reference path is walked by
  >= 7 HPRC haplotypes in the whole GFA. That was checked on the elements at audit time; 6 elements added later were not
  scanned. For rarity in large populations use `gnomAD_AF` / `CoLoRSdb_AF`.
- **Frequencies are conditional.** They are conditional on complete traversals and biased upward: haplotypes cut by the
  d9 filter inside the window carry the allele at about half the rate. `hprc_exact_allele_freq_over88` is the lower bound.
- **Cross-checks.**
  - The d9 / full-graph VCF genotypes give the same carrier set as the graph walks at 1,667 / 1,688 INDEL loci. The VCF
    counts snarl alleles, so it agrees with the element level, not the array level.
  - An audit re-derivation from the d9 GBZ (`vg paths`, vg 1.65) reproduced every walk and every element carrier for
    8 haplotypes x 2,478 loci.
- **Populations** (observed / expected carrier haplotypes, against each locus's complete haplotypes; z ignores
  haplotype pairing and linkage):
  - absent, exact allele: AFR 1.11 (z +21), AMR 0.89 (z −17), EAS 0.89 (z −7), SAS 0.92 (z −2.5);
  - present_broad: AFR 0.92 (z −6), AMR 1.10 (z +5).
  - The patient's (European) germline alleles lean toward the admixed AMR haplotypes. The alleles hiding somatic events
    lean AFR, the panel's largest and most diverse group. HPRC v1.1 has no EUR haplotype. CHM13 carries the exact allele
    at 298 absent and 34 present loci.

  Per 1000 Genomes population (`hprc_populations.tsv`), exact-allele carrier haplotypes / complete traversals, summed
  over the loci:

  | population | haplotypes | absent INDEL | absent SNV | present_broad INDEL |
  |---|---:|---:|---:|---:|
  | AFR: ACB / ASW / ESN / GWD / MKK / MSL / YRI | 14 / 2 / 2 / 16 / 2 / 8 / 2 | 0.200-0.232 | 0.052-0.112 | 0.316-0.455 |
  | AMR: CLM / PEL / PUR | 8 / 8 / 16 | 0.168-0.176 | 0.050-0.062 | 0.444-0.505 |
  | EAS: CHS / KHV | 6 / 2 | 0.171-0.174 | 0.053-0.073 | 0.440-0.458 |
  | SAS: PJL | 2 | 0.179 | 0.058 | 0.426 |
  | CHM13 | 1 | 0.151 | 0.054 | 0.360 |

- **Per individual** (`hprc_individuals.tsv`). Each HPRC haplotype carries the exact allele of a median 315 absent INDELs
  (20% of those it traverses completely; 16% of all 1,935) and 7 absent SNVs. The range is 244 (HG01123#2) to 366
  (HG01891#2). `hprc_carriers.tsv.gz` lists the carrier haplotypes of every element and allele.

## 5. Cross-evaluation

Categories (`s9_integrate.py`), RARE_AF = 0.20 on the exact-allele HPRC frequency. Rules apply in this order:

| order | category | rule |
|---|---|---|
| 1 | `other_ambiguous` | no d9 path spells the truth and no truth-spelling read ('closest', 226); the path needs a germline allele that dipcall phases to the non-event hap (20); the ALT is on both normal haps or on the event hap (2) |
| 2 | `HG008N_present_broad` / `_rare` | the non-event hap carries the ALT over the array, or (patient frame) the INFO event turns the event hap into the other hap's array allele, or a hap carries it and the truth has no event hap; exact-allele frequency >= / < 0.20. If the INFO event instead gives a novel tumor allele: `other_ambiguous` (6) |
| 3 | `other_ambiguous` | an HG008-N hap unresolved (assembly vs dipcall conflict, or no array sequence; 9); no complete HPRC haplotype walks the elements (20, 6 of them CHM13 only) |
| 4 | `HG008N_absent_HPRC_other` | everything else: both haps resolved, neither carries the ALT, HPRC haplotypes walk the elements; sub_class `exact_allele` / `element_set_other_allele` / `recombinant_pieces` |

| category | interpretation |
|---|---|
| `HG008N_present_*` | the patient's germline allele, already a graph path. The tumor's event hap now carries the other hap's allele. A tumor-only caller would call it germline anyway: the graph acts as **germline filtering** |
| `HG008N_absent_HPRC_other` | **pangenome-induced false negative**: the somatic change recreates an allele of other individuals (`exact_allele`), or is spelled by their nodes (`element_set_other_allele`, `recombinant_pieces`) |
| `other_ambiguous` | not resolved; mostly the reads that align perfectly do not spell the truth ('closest') |

Counts per category (present_broad / present_rare / absent_HPRC_other / other_ambiguous):

| set | SNV chr1-22 | INDEL chr1-22 | SNV chr1 BED | INDEL chr1 BED | SNV chr2-22 | INDEL chr2-22 |
|---|---|---|---|---|---|---|
| PacBio | 2 / 3 / 129 / 36 (170) | 41 / 17 / 1,157 / 148 (1,363) | 0 / 0 / 11 / 3 (14) | 5 / 1 / 86 / 10 (102) | 2 / 3 / 118 / 33 | 36 / 16 / 1,064 / 138 |
| ONT | 2 / 3 / 117 / 26 (148) | 26 / 10 / 806 / 147 (989) | 0 / 0 / 12 / 2 (14) | 2 / 1 / 60 / 11 (74) | 2 / 3 / 105 / 24 | 24 / 9 / 740 / 136 |
| Illumina | 2 / 2 / 97 / 21 (122) | 68 / 28 / 1,682 / 138 (1,916) | 0 / 0 / 10 / 2 (12) | 8 / 2 / 117 / 9 (136) | 2 / 2 / 87 / 19 | 59 / 26 / 1,555 / 128 |
| **union** | **2 / 3 / 150 / 46 (201)** | **75 / 30 / 1,935 / 237 (2,277)** | 0 / 0 / 12 / 5 (17) | 8 / 2 / 140 / 19 (169) | 2 / 3 / 138 / 41 | 66 / 28 / 1,781 / 217 |

- **Stricter GIAB BED.** Inside GIAB's `nogermlineoverlap` BED (union): SNV 2 / 2 / 64 / 3 and INDEL 43 / 14 / 840 / 6.
  The direction is the same. GIAB's own stricter BED excludes most of these loci.
- **Sub-classes.** The per-sub-class table (chr1-22, chr1 BED, nogermline BED) is in `tables.md`.

### Is the graph helping or hurting?

- **Germline filtering: 110 loci (4.4%).** These are real somatic events under GIAB's definition, but the tumor allele
  is the patient's other-haplotype germline allele. Any tumor-only caller would call or filter them as germline, so
  losing them costs nothing a tumor-only design could have kept.
- **Interference: 2,085 loci (84.1%).**
  - **Mechanism: repeat slippage** (`mechanism.md`, patient-frame events):
    - 1,612 of the 1,935 absent INDELs are one-base changes in homopolymers of >= 7 bp; 209 are one-unit STR changes.
    - Over all chr1-22 truths, the absorption rate of one-base homopolymer changes rises with tract length: 4.6% at
      7-9 bp, 21.9% at 10-14, 29.2% at 15-19, 38.9% at 20-29, then 24.6% at >= 30.
    - None of 450 non-repeat INDEL truths is absorbed.
    - The tumor's new length equals a length found in HPRC, which the d9 graph keeps as a branch or skip edge.
  - **Absent SNVs are mostly repeat events.** 75 of the 150 are substitutions inside repeats, and 62 are patient-frame
    homopolymer / STR length changes that GIAB wrote as GRCh38 SNVs. Only 7 are non-repeat SNVs: 6 CpG transitions at
    common SNP sites (gnomAD AF 0.02-0.20) and 1 T>G (gnomAD 0.45).
  - **The tumor reads support these somatic calls** as well as other truths: Illumina reads show the expected tumor
    allele at 96.2% of absent INDEL loci vs 95.7% of random non-absorbed truth INDELs (audit table `b4`).
- **What a linear-reference tumor-only caller would lose anyway (population PoN):**
  - SNV: the repo PoN rule tags 100% of absent SNVs, vs 8.4% of the other 8,489 chr1-22 truth SNVs; pop-AF >= 0.05 tags
    74.7% vs 1.3%. A linear caller would lose these SNVs too.
  - INDEL: the repo INDEL path has no PoN. With a population-AF filter, the absent INDELs are tagged 98.7% / 93.5% /
    69.9% at AF 1e-3 / 0.01 / 0.05, against 74.1% / 45.3% / 12.0% of the other 6,219 chr1-22 truth INDELs. Somatic
    slippage alleles are population alleles by nature; only a high AF floor keeps them in a linear pipeline.
- **The panel's composition does not explain the small present group.** Of the truth alleles that phased dipcall says
  HG008-N carries, 140 / 196 are already absorbed (biology audit). The present category is small because of the
  patient's genotype, not because HPRC v1.1 lacks EUR haplotypes.
- **Scope.** This analysis covers only the somatic truths the graph loses. The benefit side, the patient's germline
  variants that the graph absorbs so they never become candidates, is not measured here. The audit's approximate figure
  from an earlier normalised intersection is 81% of HG008-N SNVs and 73% of INDELs present in the d9 VCF.

## Verification

Four independent audits checked the first integration; each finding was then fixed or documented
(`audit/decisions.md`, whose numbers are from that run; the audit scripts and outputs are in `$D/audit/`).
- **HG008-N membership.** Independent anchoring with minimap2 asm5 on 10-kb windows, own aligner, tumor-contig
  comparison, about 60 loci read by hand.
  - 0 errors in 4,351 absent haplotypes.
  - It found 15 present INDELs wrong at the whole-array level and 24 present only at the GRCh38 coordinate. The array rule
    in s8b came from this, and present fell from 172 to 110 (present SNVs from 30 to 5).
- **HPRC membership.**
  - The d9 GBZ re-derivation agreed at 19,824 / 19,824 locus-haplotype pairs, and node coverage agreed with a whole-GFA
    scan.
  - Element-level frequencies count node walkers, not allele carriers, so the exact-allele level was added.
  - The 'closest' loci are ambiguous.
- **Logic and totals.** Every table of the first integration was recomputed from the raw files (0 differences). It found
  closest loci whose read-derived elements spell germline (now ambiguous) and germline-only Illumina reads (now separate
  read classes). The final tables were checked number by number against the files in a second pass. It found and fixed
  the chrX rows in the PoN control group and the element-vs-locus counts.
- **Interpretation.** Tumor reads support the absent INDELs; the flanking-SNP check rules out regional LOH of the event
  hap. The PoN and patient-frame mechanism statements above come from this audit.

## Assumptions and caveats

- **Definitions.** As in the Definitions section. The truth INFO coordinate is 0-based (REF matches at 2,452 / 2,452).
- **HG008-N membership.**
  - It is decided on the tandem array between unique anchors. Imperfect VNTR copies are not grown into the array.
  - 16 haplotypes have no array sequence; their loci are ambiguous unless the narrow call is a clear 'no'.
  - No HG008-N reads are on disk. Normal-assembly errors in long homopolymers are bounded only indirectly: the tumor's
    retained other haplotype disagrees with the assembly at ~6-8% of loci, the same at absorbed and other truths.
- **HPRC frequencies.** They are over complete traversals (upper-biased) and conditional on the d9 graph. "Rare"
  (< 0.20) is relative to HPRC, not to gnomAD.
- **Read paths.**
  - Complete: Illumina and PacBio all 2,478 loci; ONT all 1,137 ONT-perfect loci (1,365 in total; the ONT-only rest jobs
    covered the perfect loci only).
  - Reads drive the path choice at multi-path loci, the elements of 'closest' loci and the read-support flags, so they
    can move a few loci. Over five s9 runs (Illumina only; + PacBio 2,126 / 2,464 / 2,478 loci; + ONT 518 / 1,365) one
    locus changed category: truth 511 (chr1:148067150, outside the BED). ONT truth reads chose another ALT path there,
    which one HPRC haplotype carries, so it moved from `no_HPRC_carrier` to absent. The complete ONT reads changed no
    category (only read-support numbers).
- **Read-level reasons.** Recomputed from the re-decoded reads, they agree with the earlier miss analysis at
  1,618 / 1,642 PacBio and 1,162 / 1,162 ONT loci with a reason. The PacBio SNV reasons there came from the v5 decoder.
- **Representation duplicates.** 34 absent SNVs share one patient-frame event with another truth (one +A written as
  several GRCh38 records). Counts are per truth allele, as everywhere else.

## Files

This folder (results; `columns.tsv` explains every column of the two main tables):

| File | Content |
|---|---|
| `README.md` | this report |
| `per_variant.tsv` | one row per truth allele (2,478, 88 columns):<br>- platforms and miss classes;<br>- graph match and somatic elements;<br>- read support per platform;<br>- HPRC complete / partial / absent haplotypes, exact-allele, element-set and element frequencies, d9 / full VCF AF, per-superpopulation carriers and the carrier haplotypes;<br>- HG008-N event hap and INFO event, hap1 / hap2 allele calls (array call, length change, node-level walk), allele class, event-hap pattern;<br>- tumor verdict, dipcall, PoN / gnomAD / CoLoRSdb;<br>- category, sub_class, sub_reason, interpretation |
| `per_node.tsv` | one row per (truth allele, somatic element), 3,368 rows / 3,309 distinct elements:<br>- node ids and sequences, GRCh38 interval, type and subtype, `alone`;<br>- reads per platform;<br>- HPRC carriers / frequency, superpopulation carriers and denominators, carrier haplotypes;<br>- CHM13;<br>- node-level HG008-N hap1 / hap2 membership (`HG008N_hap*_giraffe_walk_contains`) |
| `hprc_carriers.tsv.gz` | one row per (truth_id, level, carrying haplotype with a complete traversal; HPRC and CHM13). Columns: `element_id`, `sample`, `hap`, `population_code`, `superpopulation`, `carries_element`, `carries_set_any_path`, `carries_exact_allele`.<br>- Element carriers = `level == element` (= `per_node` `hprc_carriers`).<br>- Exact-allele carriers = `level == allele & carries_exact_allele == True & sample != CHM13` (= `per_variant` `hprc_exact_allele_carriers`).<br>- Denominators = `per_variant` `hprc_n_complete` / `{SP}_complete` |
| `hprc_individuals.tsv` | per HPRC haplotype (and CHM13): loci per category it traverses completely / carries (exact allele, element set) |
| `hprc_populations.tsv` | per 1000 Genomes population: haplotypes, exact-allele carrier haplotypes and complete traversals per category |
| `indel_repeat_context.tsv`, `indel_repeat_context.md` | per absorbed INDEL truth: GRCh38- and patient-frame repeat class, unit, tract, units changed; tables per platform and category |
| `per_locus_populations.tsv`, `populations.md` | per truth allele: carrier haplotypes / complete haplotypes per superpopulation and per 1000G population, the superpopulation pattern, shared vs private (one superpopulation / population / individual); `populations.md` = the tables per category |
| `columns.tsv` | meaning of every column of `per_variant.tsv` / `per_node.tsv`, with its name in `$D/per_variant.tsv` (`tables.md` uses those `$D` names) |
| `tables.md`, `summary.json` | every table of s9: categories in four scopes, sub-classes, frequencies, superpopulations, evidence x category, PoN, read support, VCF cross-check, d9 floor |
| `mechanism.md` | patient-frame event class per category, absorption rate by repeat context, local LOH per category (class labels: `HP>=7 1-unit` = one-base change of a homopolymer >= 7 bp; `STR` = period 2-6; `VNTR>6` = period > 6; `O_lost` / `E_lost` = other / event hap lost by Illumina allele balance at phased heterozygous SNPs ±20 kb) |
| `audit/decisions.md` | every audit finding and what was done |
| `audit/b1_truth_context.py`, `b5_normal_frame.py`, `b6_local_loh.py` | the verification scripts whose tables `s11_mechanism.py` reads |

Scripts (each has a docstring with inputs, method, outputs and assumptions):

| Step | Script / job | What it does |
|---|---|---|
| s0 | `s0_variant_set.py` | the 2,478 perfect-bypass truth alleles (`variant_set.tsv`) |
| s1 | `s1_loci.py` | graph window, anchors, event window and INFO event per locus (`loci.tsv`) |
| s2 | `s2_graph_paths.py`, `s2_graph_paths.sbatch` | d9 paths spelling the truth ALT; their elements (`graph_paths.tsv`, `graph_elements.tsv`) |
| s3 | `s3_read_paths.py`, `s3_job.sh`, `s3_sub_job.sh`, `s3_rest_job.sh`, `s3_merge_s9.sbatch` | re-decodes the reads with the earlier miss analysis' reads.py (frozen run decoders) and records each read's node walk (`read_paths_P.tsv.gz`, `read_path_summary_P.tsv`) |
| s4 | `s4_hprc_walks.py`, `s4_hprc_walks.sh` | HPRC haplotype walks through every window from the d9 GFA (`hprc_local_paths.tsv.gz`, `hprc_haplotypes.tsv`, `hprc_window_coverage.tsv`) |
| s5 | `s5_hprc_vcf.py`, `s5_hprc_vcf.sh`, `s5_index_full_vcf.sh`, `s5b_populations.py` | d9 / full-graph VCF carriers by sequence (`hprc_vcf_alleles.tsv`, `hprc_vcf_haplotypes.tsv.gz`); population metadata (`hprc_sample_metadata.tsv`, raw downloads in `metadata_sources/`); HG008 ancestry (`hg008_donor_ancestry.tsv`) |
| s6 | `s6_assembly_seq.py`, `s6_minimap2.sh` | minimap2 sr REF/ALT contexts vs HG008-N / HG008-T (`asm_seq_status.tsv`, `asm_seq_summary.tsv`) |
| s7 | `s7_assembly_paths.py`, `s7_giraffe.sh`, `s7_parse.sh` | assembly windows mapped to d9 with giraffe (vg 1.65); their node walks (`asm_local_paths.tsv`) |
| s8 | `s8_dipcall.py`, `s8_dipcall.sh` | HG008-N dipcall relation (`dipcall_relation.tsv`) |
| s8b | `s8b_array_alleles.py` | whole-array HG008-N / HG008-T alleles between common unique anchors (`array_alleles.tsv`) |
| s9 | `s9_integrate.py`, `s9_integrate.sbatch` | integration and categories (`per_variant.tsv` 266 columns, `per_node.tsv`, `hprc_membership.tsv.gz`, `tables.md`, `summary.json`) |
| s10 | `s10_deliverables.py` | the tables of this folder and `columns.tsv` |
| s11 | `s11_mechanism.py` | `mechanism.md` (reads three audit tables; their scripts are in `audit/`) |
| s12 | `s12_report_numbers.py` | every number the README quotes, labelled (`report_numbers.txt`) |
| s13 | `s13_populations.py` | `per_locus_populations.tsv`, `populations.md` (which populations carry each locus' graph allele) |
| s14 | `s14_indel_repeat_context.py` | `indel_repeat_context.tsv`, `indel_repeat_context.md` (repeat context of the absorbed INDELs) |

Intermediate data, in `$D = /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/`:
- the files named above, plus `collector_reports.json` (each step's own report, with numbers and caveats) and
  `read_paths_parts/` (raw read-path parts);
- `before_fixes/` and `s9_illumina_only_20261001/`: s9 outputs before the audit fixes, and before the PacBio / ONT reads;
- `audit/`: the four audits (`hg008n_hap_membership/`, `hprc_membership/`, `logic_totals/`, `biology/`) and the report
  check (`report_check/`).

Job scripts and logs are in `tmp/graph_absorbed_somatic_20261001/`; the final run (ONT complete, merge + s9-s13) is
`final3_ont_complete.login.log`. Earlier s9 outputs are kept in `$D/s9_illumina_only_20261001/`,
`$D/s9_pacbio2464_20261001/` and `$D/s9_ont518_20261001/`.

## Reproduce

All inputs are read-only. Use `/wanglab/jshen/anaconda3/bin/python`. vg is `/scratch/wzhang/bin/vg_v1.65.0`, because the
d9 minimizer index is format 10.

1. `python s0_variant_set.py && python s1_loci.py` (login node, minutes).
2. `sbatch s2_graph_paths.sbatch` (32 CPUs, 3G; 1-15 min).
3. Read paths:
   1. `sbatch --array=0-(N-1) -c 40 --mem=75G s3_job.sh Illumina N`. Illumina took 13-15 min per part and 60-63 GB.
      PacBio: `-c 40 --mem=56G`, 1.5 h and 36 GB for 413 loci.
   2. Then, until nothing is left: `python s3_read_paths.py rest-list P TAG` and
      `sbatch --array=0-(N-1) -c 12 --mem=16G s3_rest_job.sh P TAG N` (ONT 22G).
      - Chunks are written as they finish, so a job that hits its time limit keeps what it did.
      - Cost: PacBio 7-16 worker-minutes per locus (I/O-bound); ONT 9-40.
   3. `python s3_read_paths.py merge P` (PacBio 8-10 min, 0.8 GB).
4. `sbatch s4_hprc_walks.sh` (24 processes, 5 min, 13.4 GB).
5. `sbatch s5_index_full_vcf.sh; sbatch s5_hprc_vcf.sh` (5 min, 2.1 GB); `python s5b_populations.py` (fetches the HPRC /
   IGSR / Coriell metadata into `metadata_sources/`).
6. Assemblies:
   1. `python s6_assembly_seq.py queries`, `sbatch s6_minimap2.sh` (21 GB peak, 4 min), `python s6_assembly_seq.py parse`.
   2. `python s7_assembly_paths.py windows`, `sbatch s7_giraffe.sh` (50.8 GB, 8 min), `sbatch s7_parse.sh`
      (parse + crosscheck).
   3. `sbatch s8_dipcall.sh`; `python s8b_array_alleles.py` (login, 10 min, 72 MB).
7. `python s9_integrate.py 6` (login, 3-5 min, 0.6 GB per process) or `sbatch s9_integrate.sbatch` (16 workers, 9 GB).
8. `python s10_deliverables.py && python s11_mechanism.py && python s12_report_numbers.py > report_numbers.txt`.

After any re-run, refresh the README numbers from `report_numbers.txt` (`s12_report_numbers.py` prints every number the
README quotes, labelled by section).
