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

Off-reference keys (unit = d46 key; events = keys of a tier and kind with anchors within 10 bp merged):

| tumor | kind | T1 | T2 | T3_reads_only | chr1 T1 | chr1 T2 | T1 events | T2 events |
|---|---|---|---|---|---|---|---|---|
| HG008T | SNV | 23 | 169 | 9 | 1 | 21 | 22 | 148 |
| HG008T | INDEL | 5 | 883 | 88 | 0 | 78 | 5 | 837 |
| COLO829T | SNV | 110 | 182 | 54 | 7 | 12 | 108 | 166 |
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
