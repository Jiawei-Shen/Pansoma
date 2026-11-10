# Off-reference somatic candidates (HPRC v2.1 d46, HG008T and COLO829T, 3 platforms)

Goal: a tiered list of somatic candidates among the off-reference tensors (labels with `"grch38": null`) of the six d46
sets (HG008T Illumina / PacBio / ONT, COLO829T Illumina / Fiber-seq / ONT), which the GRCh38 truth cannot label. The
evidence is the tumor and normal assemblies plus matched-normal and tumor reads on three platforms.

Data: `/scratch/jshen/data/pansoma_net_v2_runs/offref_somatic_candidates_20261109` (s1-s5). Job scripts, `jobs.txt` (job
ids, nodes, changes made while running), `targets.tsv` (assemblies, GRCh38), `bams.tsv` (reads) and the s1 task lists:
`jobs/` (copies; the scripts run from, and the code reads `targets.tsv` / `bams.tsv` at,
`tmp/offref_somatic_candidates_20261109`, where the logs are).

## Steps

| step | script | what |
|---|---|---|
| s1 | `s1_extract.py` | Every off-reference tensor, plus 40 somatic / 40 germline / 20 no-variant on-reference controls per table: A1 (ALT-carrying) read consensus over the whole tensor (`alt_q`, site `[a, b)`) and the same with the graph REF at the site (`ref_q`); QC `n_a1`, `min_share`, `n_split`. |
| s2 | `s2_queries.py`, `s2_map.sh` | Unique contexts (`common.usable`: >= 30 bp, >= 10 bases on each side of the site and >= 5 beyond the event window) mapped with minimap2 2.28 `-cx sr --cs -N 50` to the two normal haplotypes, the two tumor haplotypes and GRCh38. |
| s3 | `s3_classify.py` | `asm_class` (normal / tumor_only / neither / short_flank / no_query) by the E9d / E9h h50 strict rule; `normal_side` (ref_exact / site_only / site_plus_other / no_site_diff / no_alignment); multi-copy flags; the patient's own normal allele at the site (`norm_q`); GRCh38 anchor from the alignment of `ref_q`. |
| s4 | `s4_loci.py`, `s4_count.py`, `s4_fetch_hg008n.py` | Rows: controls, known keys, every tumor_only row, and neither + ref_exact / site_only rows that can still reach the reads-only tier (n_a1 >= 4, n_split 0, not a homopolymer / STR INDEL, key on >= 2 platforms or AF >= 0.2; this gate cuts HG008T from about 330k to 15k such rows, most of them long-read homopolymer INDEL errors). ALT and patient-REF k-mers per locus and their read counts in each tumor and normal BAM (k-mer searched only in the part of a read aligned near the site; fragments counted once). HG008-N reads are fetched from GIAB (HiFi, Illumina). |
| s5 | `s5_tiers.py` | Key-level tiers, controls and known off-reference somatic keys through the same funnel. |

Assemblies: HG008-N v6.2 (normal), HG008-T v3.2 (tumor, the assembly GIAB v0.2 was built from); COLO829BL verkko DSA
hapX / hapY (normal), COLO829T hifiasm (tumor; built from the Fiber-seq and ONT D003 reads, so not independent of those
two platforms).
Normal reads: COLO829BL Illumina, Fiber-seq, ONT (D001 + D003, local); HG008-N PacBio HiFi and Illumina (remote).

## Rules

- Event window: site +- 1, widened to the tandem runs (period 1-6) touching it and to every equivalent placement of an
  INDEL, so classes do not depend on contig orientation.
- Read rules (s5): normal present = ALT >= 2 and ALT fraction >= 0.05; absent = patient-REF k-mer >= 8 (HiFi / Fiber-seq)
  or >= 10 (Illumina / ONT), ALT <= 1 and fraction < 0.03. For homopolymer (>= 6) and STR INDELs only Illumina reads count
  (long-read normals and tumors show the same ALT fraction from errors there). Tumor support = ALT >= 3 and fraction >= 0.05.
- Tiers per key: rejected_normal; T1 = tumor_only + ref_exact / site_only, normal absent, evidence on >= 2 platforms
  (read support or a tensor of the key), clean (n_a1 >= 4, n_split 0, single copy, specific k-mer); T2 = the same with
  one platform or not clean; T3_reads_only = neither + ref_exact / site_only, normal absent, >= 2 platforms, clean (not
  in the tumor assembly: possible subclonal or assembly miss); asm_only = tumor_only without a read check.

## Review before the full run

A 26-agent adversarial review of s1-s4 (2026-11-09) found s1 correct (encoding re-derived on 42k rows, 0 mismatches) and
led to these changes before s2-s4 ran: the short-flank filter; the event window (normal_side of STR deletions had
depended on contig orientation); no NM <= 10 cap on the PAF (it removed the normal alignment of every INDEL > 10 bp);
normal alignments ranked by differences outside the window; GRCh38 anchors from the alignments (27-60 % of off-reference
SNV rows have no label anchor); patient-REF k-mers from the normal assembly for site_only rows; window-local k-mer search
(long reads carried far copies); fragment-level counting (the GIAB HG008 Illumina BAMs have no duplicate flags);
supplementary records counted; GRCh38 copy check of the k-mers; resumable remote fetch.

## Results (2026-11-10)

Lists: `results/<TUMOR>_offref_candidates.tsv` (off-reference keys in T1 / T2 / T3_reads_only, with platform sets,
read counts per sample `alt/ref/span`, label reason, anchor); full tables `results/<TUMOR>_summary.md`; every key with
its tier in `DATA/s5/<TUMOR>/keys.tsv.gz`.

Off-reference keys (unit = d46 key; events = keys of a tier and kind with anchors within 10 bp merged). Numbers of the
second s5 run (2026-11-10, with the COLO829T snvaf05 SNV keys and the germline-check rows: 2 HG008T SNV and 1 INDEL T2 keys
became rejected_normal, 3 COLO829T SNV T2 keys were added):

| tumor | kind | T1 | T2 | T3_reads_only | chr1 T1 | chr1 T2 | T1 events | T2 events |
|---|---|---|---|---|---|---|---|---|
| HG008T | SNV | 23 | 167 | 9 | 1 | 21 | 22 | 146 |
| HG008T | INDEL | 5 | 882 | 88 | 0 | 78 | 5 | 836 |
| COLO829T | SNV | 110 | 185 | 54 | 7 | 13 | 108 | 168 |
| COLO829T | INDEL | 1 | 330 | 396 | 0 | 15 | 1 | 314 |

T1 SNVs by platform sets with a tensor of the key: HG008T 13 in all three, 3 in two, 6 Illumina only, 1 PacBio only;
COLO829T 77 in all three, 19 ONT + Fiber-seq, 5 in another pair, 9 in one set.

Controls through the same funnel (keys, both tumors' own control draws):

| tumor | group | kind | keys | T1 + T2 | T3 | rejected_normal |
|---|---|---|---|---|---|---|
| HG008T | ctrl_somatic | SNV | 1656 | 1593 (96.2 %) | 0 | 13 |
| HG008T | ctrl_somatic | INDEL | 1988 | 1608 (80.9 %) | 1 | 71 |
| HG008T | ctrl_germline | both | 5192 | 9 (0.17 %) | 0 | 3924 |
| HG008T | ctrl_non | both | 2634 | 4 (0.15 %) | 4 | 747 |
| COLO829T | ctrl_somatic | SNV | 2329 | 2040 (87.6 %) | 169 | 5 |
| COLO829T | ctrl_somatic | INDEL | 1187 | 833 (70.2 %) | 35 | 17 |
| COLO829T | ctrl_germline | both | 5200 | 3 (0.06 %) | 0 | 3942 |
| COLO829T | ctrl_non | both | 2631 | 20 (0.76 %; 9 T1) | 3 | 1066 |

COLO829T somatic controls in T3 (169 SNVs) are the hifiasm tumor-assembly misses (7-13 % of the clonal truth is absent
from it, E9). Homopolymer / STR rows reach T2 at most (Illumina is their only informative platform), which is why most
somatic INDEL controls are T2.

Off-reference keys by label reason (T1 + T2):

| tumor | next to GRCh38 somatic truth (residual_partial) | in BED, no truth within 10 bp | outside BED | no unique GRCh38 node |
|---|---|---|---|---|
| HG008T | 806 (T1 15) | 79 (T1 5) | 180 (T1 6) | 12 (T1 1) |
| COLO829T | 136 (T1 16) | 392 (T1 67) | 19 (T1 8) | 75 (T1 20) |

- HG008T: most candidates are the GIAB somatic events in an off-reference representation; in-BED keys without truth are
  79 of 885 in-BED T1 + T2 keys (an upper bound on the false share; GIAB can miss off-reference events). T3 is 93 of 97
  in-BED no-truth keys: weak, report separately.
- COLO829T: in-BED no-truth keys dominate T1 (67 of 111). Their reads look like clonal somatic SNVs (median max tumor
  VAF 0.59 vs 0.61 for the somatic controls; normal absent on Illumina, Fiber-seq and ONT), and 9 on-reference
  no-variant controls pass T1 with the same profile (median VAF 0.98), so the SMaHT GRCh38 truth probably lacks some
  clonal SNVs; reads alone cannot settle this.
- Known off-reference somatic keys: HG008T E1.1 asm_d46 179 keys with a tensor: T1 10, T2 90, rejected_normal 48 (19
  have the A1 context in the HG008-N assembly, others have ~100 % or ~30-50 % ALT in both normals: germline in this
  tensor's representation). COLO829T SMaHT (E9j) 16: T1 7, T2 1, T3 1; E9k DSA keys 100: T1 32, T2 14, T3 2.

Calibration that set the read rules (Illumina normal ALT fraction in homopolymer / STR context): >= 0.15 in 98-99 % of
germline controls and 0.6-3.6 % of somatic controls; a first version that let long-read normals and a 0.05 threshold
decide graph-frame SNVs in runs rejected 96 of the 179 HG008T known keys (HiFi homopolymer errors).

Caveats: the COLO829T tumor assembly is built from the Fiber-seq and ONT reads; HG008-T v3.2 has known homopolymer
errors (GIAB v0.4 polished 14,370 regions); T1 keys whose GRCh38 locus has no Illumina reads rest on the HiFi normal alone
(e.g. `39191406:47:SNP:G>A`, HiFi normal 0 / 9); the d46 node key is the unit, so one event in two node frames counts
twice (events column).

## Training datasets (s6, 2026-11-10)

`s6_dataset.py` turns the evidence into training labels and writes `/scratch/jshen/data/pansoma_v2_tensors/<SET>_offref`
for HG008T_{Illumina,PacBio,ONT}_d46, COLO829T_{Illumina,fiberseq,ONT}_d46 and COLO829T_{Illumina,fiberseq}_d46_snvaf05
(SNV only: their INDEL directories were deleted 2026-11-10 by the snv_newpipe session; INDEL training uses d46). Before
s6, the evidence was extended: s1 on the snvaf05 SNV keys absent from d46 (86,845 rows; `s1_extract.py ... BASE_SET`),
COLO829T s2-s5 rerun with them, and a germline read check: off-reference rows whose A1 context is in the normal
assembly at <= 2 loci got loci (purpose `germ`) counted in the Illumina normal (HG008-N: 87,868 new loci fetched from
GIAB; COLO829BL local). Counts of earlier runs are reused (locus ids are content hashes; `s4_count.py ... SELECT TAG`).

Labels (per key, from every platform set of the tumor):

| label | reason | rule |
|---|---|---|
| 1 | `offref_somatic_t1` / `offref_somatic_t2` | s5 tier T1 / T2 |
| 2 | `offref_germline` | A1 context in a normal haplotype at <= 2 loci on the tensor's chromosome (HG008-N contig names; COLO829BL DSA by the SMaHT GRCh38ToDSA chain) and Illumina normal ALT >= 3 reads, fraction >= 0.15 |
| 0 | `offref_non` / `offref_non_normal_alt` | ALT context in neither assembly, a tensor in one platform only, tumor read support on <= 1 platform (normal_alt: the normal reads carry it too) |
| -1 | unchanged | everything else (reads-only, unconfirmed tumor_only, normal-assembly hits without read support, multi-copy or other-chromosome hits, ALT in neither assembly on several platforms, ...) |

Scope: only tensors labelled -1 with reason `off_reference_no_truth_match` (anchored in the confident region, no truth
allele), and not on chr1 (the test chromosome; training never reads it, and the per-tensor chr1 test stays comparable
with the source sets; `--relabel-test` labels it too). Uncertain keys keep -1 and `off_reference_no_truth_match`, so
evaluation counts them as non as before. Every off-reference record carries `offref` (proposed label, class, s5 tier,
platform sets); relabelled ones also `reason_original`. Shards, summaries and manifests are links to the source set;
site features are linked as `<SET>_offref.<KIND>.npz` in `snv_af05_20261103/site_features` and
`indel_gap_20261020/site_features`. Train with `--tensors .../<SET>_offref/tensors`; drop a class with e.g.
`--ignore-reasons offref_somatic_t2 offref_non_normal_alt`.

Relabelled tensors (chr1 excluded) next to the source set's labels (all checked with pansoma_net_v2.data.KindIndex and
site_features.load):

| set | kind | somatic (t1 / t2) | germline | non (normal_alt) | source non / somatic / germline |
|---|---|---|---|---|---|
| HG008T_Illumina_d46_offref | SNV | 24 (4 / 20) | 11,809 | 65,197 (202) | 1,852,495 / 10,221 / 290,341 |
| HG008T_Illumina_d46_offref | INDEL | 38 (0 / 38) | 11,766 | 824 (0) | 23,240 / 5,560 / 62,975 |
| HG008T_PacBio_d46_offref | SNV | 21 (5 / 16) | 10,662 | 5,296 (155) | 77,478 / 8,780 / 305,305 |
| HG008T_PacBio_d46_offref | INDEL | 35 (0 / 35) | 13,402 | 42,314 (233) | 1,522,059 / 6,902 / 186,128 |
| HG008T_ONT_d46_offref | SNV | 19 (5 / 14) | 10,709 | 8,479 (81) | 206,573 / 8,787 / 340,596 |
| HG008T_ONT_d46_offref | INDEL | 21 (0 / 21) | 10,078 | 96,070 (22) | 1,967,176 / 6,507 / 201,363 |
| COLO829T_Illumina_d46_offref | SNV | 112 (53 / 59) | 10,941 | 20,705 (100) | 920,461 / 38,262 / 309,226 |
| COLO829T_Illumina_d46_offref | INDEL | 170 (0 / 170) | 12,222 | 798 (0) | 26,088 / 1,369 / 67,967 |
| COLO829T_fiberseq_d46_offref | SNV | 108 (57 / 51) | 14,364 | 6,538 (131) | 162,510 / 37,946 / 351,891 |
| COLO829T_fiberseq_d46_offref | INDEL | 207 (1 / 206) | 15,418 | 45,000 (832) | 1,321,934 / 1,677 / 204,927 |
| COLO829T_ONT_d46_offref | SNV | 109 (55 / 54) | 17,610 | 18,378 (216) | 546,170 / 37,743 / 391,974 |
| COLO829T_ONT_d46_offref | INDEL | 153 (1 / 152) | 12,188 | 142,789 (451) | 3,779,537 / 1,718 / 241,783 |
| COLO829T_Illumina_d46_snvaf05_offref | SNV | 113 (53 / 60) | 11,488 | 35,892 (100) | 2,386,324 / 39,771 / 312,651 |
| COLO829T_fiberseq_d46_snvaf05_offref | SNV | 109 (58 / 51) | 15,407 | 8,435 (131) | 217,158 / 39,225 / 356,190 |

Most off-reference somatic keys are not in this scope: next to GRCh38 somatic truth they are already labelled 1
(residual_partial_somatic_truth, unchanged), and the novel ones lie outside the confident region or have no anchor.

The same rules on the on-reference controls (keys; each control key gets the platforms of every set with a tensor of it,
as off-reference keys do; `DATA/s6/<TUMOR>/controls.md`):

| tumor | truth | kind | -> somatic | -> germline | -> non | uncertain |
|---|---|---|---|---|---|---|
| HG008T | somatic | SNV | 96.2 % | 0.2 % | 0.2 % | 3.4 % |
| HG008T | somatic | INDEL | 80.9 % | 0.6 % | 5.5 % | 13.0 % |
| HG008T | germline | SNV | 0 % | 86.5 % | 0.9 % | 12.6 % |
| HG008T | germline | INDEL | 0.3 % | 51.2 % | 16.9 % | 31.5 % |
| HG008T | no truth | SNV | 0 % | 3.9 % | 57.3 % | 38.8 % |
| HG008T | no truth | INDEL | 0.3 % | 13.1 % | 53.0 % | 33.6 % |
| COLO829T | somatic | SNV | 87.6 % | 0 % | 0.1 % | 12.3 % |
| COLO829T | somatic | INDEL | 70.2 % | 0.3 % | 5.5 % | 24.0 % |
| COLO829T | germline | SNV | 0 % | 81.9 % | 0.9 % | 17.1 % |
| COLO829T | germline | INDEL | 0.1 % | 49.6 % | 14.9 % | 35.4 % |
| COLO829T | no truth | SNV | 1.2 % | 2.7 % | 46.7 % | 49.4 % |
| COLO829T | no truth | INDEL | 0.3 % | 10.3 % | 44.4 % | 44.9 % |

- Germline INDEL controls that become non: in the review, the Illumina normal had no ALT at ~93 % of them (fraction
  < 0.03 at most): the long-read germline INDEL labels are largely homopolymer read errors at germline positions.
- No-truth INDEL controls that become germline: their A1 context is in the normal assembly and >= 15 % of the Illumina
  normal reads carry it, so they look like germline variants the on-reference labels call `confident_no_truth_allele`.
- No-truth COLO829T SNVs that become somatic look clonal (median VAF 0.98, normal absent on three platforms).
- The snvaf05-only SNV keys (AF 0.05-0.08) have no controls of their own; somatic is 0.1 % (Illumina) and 2.3 %
  (Fiber-seq) of that band's on-reference non + somatic tensors.
- A second adversarial review (18 agents) checked the dataset format against data.py / train.py (lossless JSON rewrite,
  only intended tensors change, caches, links) and led to: site-feature links, the same-chromosome germline rule, the
  manifest totals / reasons, control platforms, and keeping chr1 unchanged.
- For RIS training, the `_offref` sets and their site-feature links still have to be copied there.
