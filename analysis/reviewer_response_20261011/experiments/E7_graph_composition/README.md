# E7: in-silico graph composition (haplotype count x frequency threshold)

Reviewer points:

- **R1.2** (composition and size of the pangenome). The reviewer says there is no systematic analysis of how the
  number and diversity of haplotypes affect node pileups and germline identification. The graph decides which alleles
  are represented, and so which variants are suppressed. A larger HPRC release could absorb more somatic mutations that
  overlap population alleles.
- **R3.m3** (composition of the reference). The reviewer asks whether benchmark individuals or their close relatives
  contributed to the graph. Step s6 answers this.

The results also feed R1.4 / R3.2 (which somatic alleles the graph absorbs). The design is in
`../../response_to_reviewers_detailed.md`, section R1.2, "E7 in-silico 图组成分析".

## Summary

Status: **DONE for HPRC v1.1** (steps s1-s7 finished on 2026-10-13). **PENDING**: the HPRC v2.1 haplotype axis
(k up to 462), which needs the v2.1 full-graph VCF; downloading it is a user decision. Also pending: the optional E7b
measured point and a figure. No E7 job is queued or held.

- **The frequency threshold decides what the graph contains. The haplotype count does not.** At a relative threshold of
  10 %, going from k = 16 to all 89 HPRC v1.1 haplotypes changes:
  - the share of the donor's germline SNVs in the graph: HG008 84.5 -> 86.9 %, COLO829 84.7 -> 87.2 %,
    HCC1395 85.6 -> 88.1 %;
  - the share of somatic truth INDELs in the graph: HG008 25.9 -> 25.6 %, COLO829 22.7 -> 22.3 %,
    HCC1395 3.3 -> 3.3 %.
- **Without a frequency filter (d = 1), more haplotypes absorb more somatic INDELs.** Germline coverage saturates
  (HG008 germline SNV: 90.7 / 94.0 / 96.0 % at k = 16 / 32 / 89). Somatic INDEL absorption keeps rising
  (40.4 / 49.0 / 59.0 %).
- On the full v1.1 set, lowering the filter from 9 haplotypes (the d9 graph) to 1 haplotype changes HG008 as follows:
  - germline SNV coverage gains 9.1 points (86.9 -> 96.0 %);
  - somatic INDEL absorption gains 33.4 points (25.6 -> 59.0 %);
  - somatic SNV absorption goes from 2.4 to 4.6 %.
- **Absorption is concentrated in repeat INDELs.** At k = 89, d = 1 the graph contains, for HG008:
  - 70.9 % of homopolymer >= 7 bp somatic INDELs;
  - 66.9 % of STR somatic INDELs;
  - 8.5 % of the other somatic INDELs.

  HCC1395 absorbs much less (7.4 %). This follows its truth composition: only 234 of its 1,559 INDELs are HP >= 7.
- **Ancestry matters more than count.** At equal k the non-AFR haplotypes contain more of each donor's germline alleles
  than the AFR haplotypes. Example at k = 42, f = 10 %, germline SNV: 85.7 vs 79.4 % for HG008. The AFR haplotypes absorb
  more somatic INDELs (25.2 vs 20.7 %).
- **Calibration against the real d9 graph:**
  - The d9 filter is sharp on the v1.1 carrier count. Of the HPRC SNV alleles with >= 9 carriers, 96.2-98.5 % are in
    the d9 VCF. Of those with 1-8 carriers, 2.0-7.2 % are.
  - The model predicts the d9 content (HG008 germline SNV: 86.9 % predicted, 87.0 % found).
  - On identical reads the model explains germline SNV loss almost entirely: 98.1-98.2 % of the germline SNV candidates
    lost on d9 are predicted.
  - It explains only part of the somatic loss: 40-54 % of the somatic INDEL loss and 14-54 % of the somatic SNV loss.
    Somatic INDELs with no HPRC carrier at all still lose their tensor on d9 in 34.0 % (HG008 PacBio) and 26.5 %
    (COLO829 fiber-seq) of cases.
- **Relatedness (R3.m3): no HPRC v1.1 individual is a close relative of the HG008, COLO829 or HCC1395 donor.** The most
  similar individual carries 7.7-8.2 % of the donor's rare germline SNVs. That is at most 1.23x the median of its
  superpopulation.

## Question

1. Build a graph from k haplotypes and keep the alleles carried by at least d of them. What fraction of a donor's
   germline alleles is then in the graph? These alleles stop being candidates (quantity (a)). What fraction of the
   somatic truth alleles is in the graph? These are absorbed (quantity (b)). Which axis drives (a) and (b): k or d?
2. How do (a) and (b) split by SNV / INDEL, by repeat context (HP >= 7, STR, other) and by the ancestry of the drawn
   haplotypes?
3. Does the model predict the real d9 graph: its allele content, and the candidate loss measured on identical reads
   (GRCh38-only vs d9 / d46 tensor sets)?
4. (R3.m3) Is any HPRC v1.1 individual a close relative of a benchmark donor?

## Data

- **HPRC v1.1 haplotype alleles**
  - `/scratch/jshen/data/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz`: vg deconstruct -> vcfbub -> vcfwave. It has 45
    sample columns: CHM13 (haploid) plus 44 diploid samples, giving 89 non-GRCh38 haplotypes. GRCh38 is the backbone
    and is always in the graph.
  - `/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/full_vcf_index/hprc-v1.1-mc-grch38.raw.vcf.gz`:
    the deconstruct VCF before vcfbub / vcfwave. It is used to repair the genotypes of vcfwave-decomposed records (see
    Method).
  - Superpopulation of each haplotype: `analysis/graph_absorbed_somatic_20261001/hprc_individuals.tsv`. It lists AFR 23,
    AMR 16, EAS 4 and SAS 1 individuals, plus CHM13. This confirms the v1.1 counts marked [待核实] in R3.m3.
- **Real d9 graph content**: `/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz`, the
  deconstruct VCF of the d9 graph Pansoma aligns to.
- **Sample alleles**: the truth tables the labels were built from, `truth/{somatic,germline}.graph.tsv` of the tensor
  sets HG008T_PacBio, COLO829T_fiberseq and HCC1395T_PacBio_d46 (`/scratch/jshen/data/pansoma_v2_tensors`):

| Sample | Somatic truth (BED) | Germline truth (BED) |
|---|---|---|
| HG008 | GIAB HG008-T somatic smvar benchmark v0.2 (`..._all.bed`) | HG008-N dipcall (GRCh38) |
| COLO829 | `COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz` (SMaHT v2 easy/difficult/extreme union BED) | COLO829BL dipcall (hg38) |
| HCC1395 | SEQC2 v1.2.1 (SEQC2 high-confidence regions v1.2) | HCC1395BL DeepVariant hybrid, QUAL >= 10 (QUAL10 BED) |

- **Measured candidate status**: `tensors/{somatic,germline}.recall.tsv`. Status `tensor_representative` counts as
  "tensor".
  - HG008T and COLO829T: three platforms x three arms (`<set>_linear` = GRCh38-only, `<set>` = HPRC v1.1 d9,
    `<set>_d46` = HPRC v2.1 d46).
  - HCC1395T: d46 only.
  - The truth tables are identical across the sets: one md5 per (sample, role), checked in the s3_s4_s6 log.
- **Reference**: GRCh38 no-alt analysis set
  (`/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`). Jobs read a node-local copy
  (`local_fasta.sh`) because random access over BeeGFS stalled on busy nodes.
- **Scope**: chr1-22. The main scope "bed" is FILTER PASS and inside the truth's BED. Also kept: "bed_chr1", and "all"
  (every chr1-22 allele).

## Method

1. **s1 (`s1_hprc_alleles.py`, per chromosome)**: computes per-haplotype carriers of every HPRC v1.1 primitive allele.
   - vcfwave GT repair. vcfwave keeps incomplete genotypes for the records it decomposed (INFO ORIGIN set). For every
     decomposed bubble, each raw ALT allele is re-expressed as a minimum exact cover by that bubble's wave primitives. A
     primitive's carriers are the union of the raw carriers of its parent alleles. Untouched records keep their wave GT.
   - All alleles then go through `bcftools norm -f GRCh38 -c w` (left-align) and `common.primitives`: MNPs are split
     into SNVs; INS / DEL keep one anchor base; anything else is COMPLEX.
   - Identical keys are merged as the union of carriers. Each allele stores an 89-bit carrier mask plus its repeat
     context in the GRCh38 frame:
     - HP >= 7: unit 1, tract >= 7 bp;
     - HP2-6: unit 1, tract 2-6 bp;
     - STR: unit 2-6, >= 3 copies;
     - other: everything else;
     - SV: > 50 bp.
2. **s2 (`s2_sample_alleles.py`, per sample)**: normalizes the somatic and germline truth alleles in the same way.
   A truth allele with several primitives is in a graph only if one haplotype set carries all of them (AND of the
   masks). Such alleles go to stratum "other" and are not in the SNV / INDEL strata.
3. **s3 (`s3_curves.py`)**: draws in-silico graphs.
   - Pools: all 89 haplotypes; AFR (46); non-AFR (AMR + EAS + SAS, 42). CHM13 is left out of both population pools.
   - Note (2026-10-15): in the real d9 / d46 graphs CHM13 is an unfiltered REFERENCE path (`vg paths -M`), whereas this
     model counts CHM13 as one of the 89 haplotypes that the carrier threshold filters. The effect on the curves is small
     (one haploid genome), but CHM13-private alleles are always in the real graphs.
   - k in {2, 4, 8, 16, 32, 64} plus the whole pool; for AFR also k = 42 to match non-AFR. There are 20 seeded draws
     per k (seed 20261013) and one draw for the whole pool: 343 draws in total.
   - An allele is in the graph if at least d of the drawn haplotypes carry it. Thresholds:
     - absolute d in {1, 2, 5, 9, 18};
     - relative f in {>= 1 haplotype, 2 %, 5 %, 10 %, 20 %}, with d = max(1, ceil(f k)).
   - Outputs: `curves.tsv` (mean / sd / min / max over draws) and `graph_size.tsv` (HPRC alleles kept per graph).
4. **s5 (`s5_d9_alleles.py`, per chromosome)**: lists the alleles of the real d9 graph from its own deconstruct VCF,
   using the same normalization (`bcftools norm -m -any`, then `primitives`). This is a lower bound: d9 COMPLEX alleles
   are not decomposed.
5. **s4 (`s4_calibration.py`)**: compares the real graphs with the model, allele by allele.
   - "Predicted in d9" means >= 9 of the 89 haplotypes carry the allele (the model for d9 = 9 of 90 haplotypes).
   - Measured loss on d9 (or d46) = alleles with a tensor on the GRCh38-only arm and none on d9 (d46), divided by the
     alleles with a GRCh38-only tensor.
   - Outputs: `calibration.tsv`, `calibration_bins.tsv` (loss by carrier count), `d9_retention.tsv` (d9 VCF presence by
     carrier count), `d9_content.tsv`.
6. **s6 (`s6_relatedness.py`)**: takes each normal's PASS, in-BED germline SNVs that are rare in HPRC (1-2 of 89
   carriers). It counts how many of them each HPRC individual carries, as a share of the total and as a ratio to the
   superpopulation median. The script's premise: a first-degree relative would share about half of them. Output:
   `relatedness.tsv`.
7. **s7 (`s7_tables.py`, run by `collect.sh`)**: writes `tables.md` (sections A-J) and `summary.json`.

## Scripts and jobs

Job scripts and logs are in `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E7_graph_composition/`. Large
outputs are in `/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/` (WORK).

| Step | Job script | Slurm job(s) | State | Output |
|---|---|---|---|---|
| s1 | `s1_hprc_alleles.sh` (general, -c 3, 6 G) | 387149_1-9 (chr1-9); 387261_10-20, _22 (chr10-20, chr22). chr21 ran outside Slurm at 10:17-10:19, before the array (no job log; `chr21.json` complete). | COMPLETED. The 387149 chr21 task was cancelled before it started. | WORK/hprc_alleles/chr*.tsv.gz, .json |
| s2 | `s2_sample_alleles.sh` (general, -c 2, 10 G) | 387260_0-2 | COMPLETED. Earlier attempts 387172 and 387186 were cancelled and superseded. | WORK/sample_alleles/{HG008,COLO829,HCC1395}.{somatic,germline}.tsv.gz |
| s5 | `s5_d9_alleles.sh` (general, -c 3, 4 G) | 387262_1-22 | COMPLETED | WORK/d9_alleles/chr*.tsv.gz, .json |
| s3 + s4 + s6 | `s3_s4_s6.sh` (general, -c 2, 24 G) | 387305 (miao, 28 min) | COMPLETED | `curves.tsv`, `graph_size.tsv`, WORK/joined.npz, WORK/hist_draws.npy, `calibration*.tsv`, `d9_*.tsv`, `relatedness.tsv` |
| s7 | `collect.sh` -> `s7_tables.py` | login node, 2026-10-13 22:55 | done | `tables.md`, `summary.json` |

`collect.sh` submits s3/s4/s6 (`sbatch --wait`) only if `calibration.tsv` is missing. It exists now, so a rerun only
rebuilds `tables.md` and `summary.json`.

Input bookkeeping, from the logs and the WORK json files:

- **HPRC v1.1, chr1-22**: 27,409,686 wave records, of which 6,947,944 were decomposed by vcfwave. All 1,731,073
  decomposed bubbles were found in the raw VCF. The result is 28,057,183 allele keys, 27,915,130 of them with at least
  one carrier.
- **Sample alleles**: every table allele was normalized:

| Sample | Somatic | Germline |
|---|---|---|
| HG008 | 17,186 | 5,061,892 |
| COLO829 | 43,947 | 5,062,184 |
| HCC1395 | 41,072 | 4,828,024 |

- **Joined**: 15,054,305 sample alleles, of which 13,600,577 have all their keys in the HPRC table.
- **d9 deconstruct VCF**: 9,222,127 records, giving 9,703,007 keys (106,088 COMPLEX primitives).

## Results

All percentages are copied from `tables.md`, `curves.tsv` and the other TSVs in this folder. Scope: PASS, in BED,
chr1-22, unless stated.

### 1. Allele counts (tables.md section A; chr1 in brackets)

| Sample | germline SNV | germline INDEL | somatic SNV | somatic INDEL | somatic INDEL_HP>=7 | somatic INDEL_not_HP>=7 |
|---|---|---|---|---|---|---|
| HG008 | 3,737,022 (294,552) | 964,490 (80,280) | 8,543 (697) | 7,933 (582) | 5,840 (427) | 2,093 (155) |
| COLO829 | 3,729,123 (290,501) | 977,293 (80,420) | 42,035 (2,443) | 1,912 (131) | 1,110 (62) | 802 (69) |
| HCC1395 | 3,135,770 (257,092) | 753,706 (66,143) | 39,447 (3,440) | 1,559 (124) | 234 (23) | 1,325 (101) |

### 2. The two axes: % of alleles in the in-silico graph (all-haplotype pool; mean ± sd over 20 draws)

Columns: (a) germline SNV, (a) germline INDEL, (b) somatic SNV, (b) somatic INDEL, (b) somatic INDEL HP >= 7,
(b) somatic INDEL not HP >= 7.

**2a. Haplotype count at a fixed relative threshold f = 10 %** (`curves.tsv`, threshold_type rel, 0.1; tables.md C for
HG008)

| Sample | k | d | germ SNV | germ INDEL | som SNV | som INDEL | som HP>=7 | som not HP>=7 |
|---|---|---|---|---|---|---|---|---|
| HG008 | 8 | 1 | 85.1 ± 1.2 | 74.3 ± 1.0 | 2.7 ± 0.1 | 30.6 ± 0.6 | 37.0 ± 0.9 | 12.7 ± 0.5 |
| HG008 | 16 | 2 | 84.5 ± 0.9 | 73.3 ± 0.8 | 2.4 ± 0.1 | 25.9 ± 0.6 | 31.4 ± 0.8 | 10.7 ± 0.5 |
| HG008 | 32 | 4 | 84.7 ± 0.5 | 73.3 ± 0.4 | 2.3 ± 0.1 | 23.1 ± 0.4 | 28.0 ± 0.6 | 9.4 ± 0.4 |
| HG008 | 64 | 7 | 86.0 ± 0.2 | 75.0 ± 0.2 | 2.4 ± 0.0 | 24.5 ± 0.2 | 29.8 ± 0.3 | 9.8 ± 0.2 |
| HG008 | 89 | 9 | 86.9 | 76.1 | 2.4 | 25.6 | 31.1 | 10.2 |
| COLO829 | 8 | 1 | 85.3 ± 1.2 | 73.8 ± 1.0 | 0.3 ± 0.0 | 26.4 ± 1.1 | 34.0 ± 1.5 | 15.9 ± 0.9 |
| COLO829 | 16 | 2 | 84.7 ± 0.9 | 72.8 ± 0.8 | 0.2 ± 0.0 | 22.7 ± 0.7 | 28.4 ± 1.1 | 14.8 ± 0.4 |
| COLO829 | 32 | 4 | 85.0 ± 0.5 | 72.8 ± 0.4 | 0.2 ± 0.0 | 20.1 ± 0.7 | 24.9 ± 1.0 | 13.6 ± 0.6 |
| COLO829 | 64 | 7 | 86.3 ± 0.2 | 74.5 ± 0.2 | 0.2 ± 0.0 | 21.6 ± 0.4 | 26.6 ± 0.5 | 14.7 ± 0.4 |
| COLO829 | 89 | 9 | 87.2 | 75.6 | 0.2 | 22.3 | 27.2 | 15.5 |
| HCC1395 | 8 | 1 | 86.2 ± 1.3 | 77.7 ± 1.0 | 0.6 ± 0.1 | 3.9 ± 0.3 | 18.6 ± 1.8 | 1.3 ± 0.2 |
| HCC1395 | 16 | 2 | 85.6 ± 1.0 | 77.0 ± 0.8 | 0.5 ± 0.1 | 3.3 ± 0.2 | 14.9 ± 0.9 | 1.3 ± 0.1 |
| HCC1395 | 32 | 4 | 85.9 ± 0.5 | 77.2 ± 0.4 | 0.5 ± 0.0 | 3.0 ± 0.2 | 13.5 ± 0.9 | 1.2 ± 0.1 |
| HCC1395 | 64 | 7 | 87.2 ± 0.2 | 78.8 ± 0.2 | 0.5 ± 0.0 | 3.2 ± 0.1 | 14.5 ± 0.7 | 1.2 ± 0.1 |
| HCC1395 | 89 | 9 | 88.1 | 79.9 | 0.5 | 3.3 | 14.5 | 1.3 |

**2b. Haplotype count without a frequency filter (d = 1)** (tables.md B)

| Sample | k | germ SNV | germ INDEL | som SNV | som INDEL | som HP>=7 | som not HP>=7 |
|---|---|---|---|---|---|---|---|
| HG008 | 8 | 85.1 ± 1.2 | 74.3 ± 1.0 | 2.7 ± 0.1 | 30.6 ± 0.6 | 37.0 ± 0.9 | 12.7 ± 0.5 |
| HG008 | 16 | 90.7 ± 0.6 | 81.3 ± 0.5 | 3.3 ± 0.1 | 40.4 ± 0.6 | 48.8 ± 0.8 | 17.2 ± 0.3 |
| HG008 | 32 | 94.0 ± 0.2 | 85.4 ± 0.2 | 3.8 ± 0.1 | 49.0 ± 0.5 | 58.9 ± 0.6 | 21.4 ± 0.3 |
| HG008 | 64 | 95.5 ± 0.1 | 87.3 ± 0.0 | 4.4 ± 0.0 | 56.3 ± 0.2 | 67.6 ± 0.2 | 24.6 ± 0.2 |
| HG008 | 89 | 96.0 | 87.9 | 4.6 | 59.0 | 70.9 | 25.8 |
| COLO829 | 8 | 85.3 ± 1.2 | 73.8 ± 1.0 | 0.3 ± 0.0 | 26.4 ± 1.1 | 34.0 ± 1.5 | 15.9 ± 0.9 |
| COLO829 | 16 | 90.9 ± 0.6 | 80.9 ± 0.5 | 0.4 ± 0.0 | 34.3 ± 0.8 | 44.9 ± 1.4 | 19.7 ± 0.5 |
| COLO829 | 32 | 94.1 ± 0.2 | 85.0 ± 0.1 | 0.5 ± 0.0 | 41.1 ± 0.6 | 54.6 ± 1.1 | 22.5 ± 0.4 |
| COLO829 | 64 | 95.5 ± 0.1 | 87.0 ± 0.1 | 0.6 ± 0.0 | 46.8 ± 0.3 | 63.1 ± 0.6 | 24.3 ± 0.2 |
| COLO829 | 89 | 96.0 | 87.6 | 0.7 | 49.2 | 66.7 | 24.9 |
| HCC1395 | 8 | 86.2 ± 1.3 | 77.7 ± 1.0 | 0.6 ± 0.1 | 3.9 ± 0.3 | 18.6 ± 1.8 | 1.3 ± 0.2 |
| HCC1395 | 16 | 91.7 ± 0.6 | 84.3 ± 0.5 | 0.7 ± 0.0 | 4.9 ± 0.2 | 23.6 ± 1.3 | 1.6 ± 0.1 |
| HCC1395 | 32 | 94.9 ± 0.2 | 88.0 ± 0.1 | 0.8 ± 0.0 | 6.0 ± 0.2 | 29.2 ± 1.4 | 1.8 ± 0.1 |
| HCC1395 | 64 | 96.3 ± 0.1 | 89.6 ± 0.0 | 0.9 ± 0.0 | 6.9 ± 0.1 | 34.7 ± 0.7 | 2.0 ± 0.1 |
| HCC1395 | 89 | 96.8 | 90.1 | 1.0 | 7.4 | 38.0 | 2.0 |

**2c. Threshold at k = 89 (all v1.1 haplotypes; d = 9 is the model of the d9 graph)** (tables.md B)

| Sample | d | germ SNV | germ INDEL | som SNV | som INDEL | som HP>=7 | som not HP>=7 |
|---|---|---|---|---|---|---|---|
| HG008 | 1 | 96.0 | 87.9 | 4.6 | 59.0 | 70.9 | 25.8 |
| HG008 | 2 | 94.6 | 86.3 | 3.9 | 51.3 | 61.6 | 22.4 |
| HG008 | 5 | 91.3 | 82.0 | 3.0 | 36.8 | 44.2 | 16.3 |
| HG008 | 9 | 86.9 | 76.1 | 2.4 | 25.6 | 31.1 | 10.2 |
| HG008 | 18 | 75.7 | 62.2 | 1.6 | 12.6 | 15.4 | 4.6 |
| COLO829 | 1 | 96.0 | 87.6 | 0.7 | 49.2 | 66.7 | 24.9 |
| COLO829 | 2 | 94.7 | 85.9 | 0.5 | 42.5 | 56.6 | 23.1 |
| COLO829 | 5 | 91.6 | 81.5 | 0.3 | 31.5 | 40.5 | 19.2 |
| COLO829 | 9 | 87.2 | 75.6 | 0.2 | 22.3 | 27.2 | 15.5 |
| COLO829 | 18 | 76.0 | 61.6 | 0.1 | 11.0 | 13.2 | 8.0 |
| HCC1395 | 1 | 96.8 | 90.1 | 1.0 | 7.4 | 38.0 | 2.0 |
| HCC1395 | 2 | 95.5 | 88.8 | 0.8 | 5.7 | 27.8 | 1.8 |
| HCC1395 | 5 | 92.4 | 85.2 | 0.6 | 4.6 | 20.9 | 1.7 |
| HCC1395 | 9 | 88.1 | 79.9 | 0.5 | 3.3 | 14.5 | 1.3 |
| HCC1395 | 18 | 77.0 | 66.5 | 0.4 | 1.9 | 9.0 | 0.7 |

**2d. chr1 only (scope bed_chr1, k = 89; d = 1 / d = 9)** (`curves.tsv`)

| Sample | germ SNV chr1 | germ INDEL chr1 | som SNV chr1 | som INDEL chr1 | som HP>=7 chr1 |
|---|---|---|---|---|---|
| HG008 | 96.2 / 87.1 | 87.9 / 76.1 | 3.7 / 2.4 | 60.8 / 26.8 | 73.8 / 32.1 |
| COLO829 | 95.7 / 87.0 | 87.4 / 75.6 | 0.5 / 0.2 | 55.0 / 29.0 | 69.4 / 30.6 |
| HCC1395 | 96.9 / 88.9 | 89.8 / 80.1 | 0.6 / 0.1 | 12.1 / 5.6 | 56.5 / 26.1 |

chr1 follows chr1-22. The chr1 somatic INDEL counts are small (582 / 131 / 124).

### 3. Repeat context of the absorbed somatic INDELs (k = 89; % in graph at d = 1 / d = 9; `curves.tsv`)

| Sample | HP>=7 | STR | HP2-6 + other | INS | DEL |
|---|---|---|---|---|---|
| HG008 | 70.9 / 31.1 (n 5,840) | 66.9 / 27.3 (622) | 8.5 / 3.0 (1,471) | 61.9 / 27.0 (4,585) | 55.2 / 23.7 (3,348) |
| COLO829 | 66.7 / 27.2 (1,110) | 71.4 / 41.7 (206) | 8.9 / 6.4 (596) | 53.8 / 26.7 (888) | 45.1 / 18.5 (1,024) |
| HCC1395 | 38.0 / 14.5 (234) | 37.0 / 14.8 (27) | 1.3 / 1.0 (1,298) | 22.3 / 10.8 (157) | 5.8 / 2.4 (1,402) |

### 4. Ancestry of the drawn haplotypes: AFR / non-AFR at equal k (tables.md D; mean of 20 draws, k = 42 = whole non-AFR pool)

| Sample | threshold | k | germ SNV | germ INDEL | som SNV | som INDEL | som HP>=7 |
|---|---|---|---|---|---|---|---|
| HG008 | d = 1 | 8 | 80.7 / 86.3 | 70.5 / 75.1 | 2.7 / 2.5 | 31.1 / 26.6 | 37.8 / 31.5 |
| HG008 | d = 1 | 42 | 91.7 / 95.0 | 84.1 / 85.9 | 4.0 / 3.7 | 52.1 / 46.2 | 63.0 / 54.9 |
| HG008 | f = 10 % | 16 | 79.6 / 85.5 | 69.0 / 73.9 | 2.4 / 2.2 | 27.0 / 22.3 | 33.1 / 26.5 |
| HG008 | f = 10 % | 42 | 79.4 / 85.7 | 68.7 / 74.0 | 2.4 / 2.0 | 25.2 / 20.7 | 30.9 / 24.5 |
| COLO829 | d = 1 | 8 | 80.8 / 86.6 | 70.0 / 74.6 | 0.3 / 0.2 | 28.2 / 21.0 | 37.0 / 25.5 |
| COLO829 | d = 1 | 42 | 91.7 / 95.1 | 83.8 / 85.5 | 0.6 / 0.4 | 44.7 / 36.6 | 60.5 / 47.0 |
| COLO829 | f = 10 % | 16 | 79.6 / 85.8 | 68.5 / 73.3 | 0.2 / 0.2 | 24.8 / 17.6 | 32.2 / 20.8 |
| COLO829 | f = 10 % | 42 | 79.6 / 86.1 | 68.2 / 73.5 | 0.2 / 0.1 | 23.2 / 16.0 | 29.7 / 18.7 |
| HCC1395 | d = 1 | 8 | 81.6 / 87.6 | 73.9 / 78.5 | 0.5 / 0.6 | 3.9 / 3.3 | 19.0 / 14.2 |
| HCC1395 | d = 1 | 42 | 92.4 / 95.9 | 86.7 / 88.4 | 0.8 / 0.7 | 6.1 / 5.2 | 31.4 / 24.4 |
| HCC1395 | f = 10 % | 16 | 80.4 / 86.8 | 72.6 / 77.5 | 0.5 / 0.5 | 3.5 / 3.0 | 16.9 / 12.1 |
| HCC1395 | f = 10 % | 42 | 80.3 / 87.1 | 72.5 / 77.7 | 0.4 / 0.5 | 3.3 / 2.8 | 16.1 / 11.1 |

### 5. Graph size: HPRC v1.1 alleles kept, chr1-22, millions (tables.md E)

| k | d | SNV | INDEL <= 50 bp | of which HP>=7 | SV > 50 bp |
|---|---|---|---|---|---|
| 2 | 1 | 4.536 | 1.273 | 0.588 | 0.036 |
| 2 | 2 | 1.474 | 0.313 | 0.135 | 0.006 |
| 8 | 1 | 8.898 | 2.749 | 1.270 | 0.096 |
| 8 | 2 | 5.061 | 1.432 | 0.678 | 0.031 |
| 32 | 1 | 14.697 | 4.638 | 2.065 | 0.226 |
| 32 | 2 | 9.427 | 3.058 | 1.441 | 0.086 |
| 32 | 9 | 3.841 | 0.969 | 0.455 | 0.018 |
| 32 | 18 | 1.785 | 0.336 | 0.141 | 0.006 |
| 89 | 1 | 21.077 | 6.402 | 2.704 | 0.422 |
| 89 | 2 | 13.694 | 4.490 | 2.073 | 0.168 |
| 89 | 9 | 6.870 | 2.124 | 1.015 | 0.046 |
| 89 | 18 | 4.698 | 1.271 | 0.605 | 0.025 |

### 6. Calibration against the real graphs

**6a. How sharp is the d9 filter?** The table gives the % of HPRC v1.1 alleles found as a simple allele in the d9
deconstruct VCF, by carrier count of 89 (chr1-22; tables.md H).

| Stratum | 1 | 2-4 | 5-6 | 7 | 8 | 9 | 10-12 | 13-20 | 21-44 | 45-89 |
|---|---|---|---|---|---|---|---|---|---|---|
| SNV | 2.0 | 2.1 | 4.9 | 5.6 | 7.2 | 96.2 | 96.6 | 96.5 | 97.6 | 98.5 |
| INDEL | 1.2 | 2.8 | 5.5 | 7.4 | 9.6 | 87.1 | 88.7 | 91.2 | 93.4 | 95.4 |
| INDEL_HP>=7 | 0.9 | 2.1 | 4.3 | 6.0 | 7.9 | 91.0 | 92.3 | 94.1 | 95.4 | 96.5 |
| INDEL_not_HP>=7 | 1.3 | 3.4 | 6.5 | 8.7 | 11.1 | 83.6 | 85.4 | 88.4 | 91.4 | 94.6 |
| SV>50 | 2.9 | 6.9 | 14.2 | 19.3 | 26.1 | 62.9 | 64.0 | 69.1 | 76.0 | 88.6 |

**6b. Predicted vs actual d9 content of the sample alleles** (tables.md I; %)

| Sample | Role | Stratum | n | any HPRC carrier | predicted d9 (>= 9) | in d9 VCF | in d9 VCF if predicted / if not |
|---|---|---|---|---|---|---|---|
| HG008 | germline | SNV | 3,737,022 | 96.0 | 86.9 | 87.0 | 98.9 / 7.7 |
| HG008 | germline | INDEL | 964,490 | 87.9 | 76.1 | 75.1 | 97.2 / 5.0 |
| HG008 | somatic | SNV | 8,543 | 4.6 | 2.4 | 2.5 | 91.2 / 0.3 |
| HG008 | somatic | INDEL | 7,933 | 59.0 | 25.6 | 25.2 | 94.0 / 1.5 |
| COLO829 | germline | SNV | 3,729,123 | 96.0 | 87.2 | 87.2 | 98.9 / 7.3 |
| COLO829 | germline | INDEL | 977,293 | 87.6 | 75.6 | 74.6 | 97.2 / 4.6 |
| COLO829 | somatic | SNV | 42,035 | 0.7 | 0.2 | 0.2 | 95.0 / 0.0 |
| COLO829 | somatic | INDEL | 1,912 | 49.2 | 22.3 | 21.1 | 89.4 / 1.5 |
| HCC1395 | germline | SNV | 3,135,770 | 96.8 | 88.1 | 88.6 | 99.5 / 7.7 |
| HCC1395 | germline | INDEL | 753,706 | 90.1 | 79.9 | 79.4 | 98.2 / 4.7 |
| HCC1395 | somatic | SNV | 39,447 | 1.0 | 0.5 | 0.5 | 100.0 / 0.0 |
| HCC1395 | somatic | INDEL | 1,559 | 7.4 | 3.3 | 3.3 | 98.0 / 0.1 |

**6c. Predicted vs measured candidate loss on identical reads** (tables.md F; %).

Columns:

- "tensor" = status `tensor_representative`.
- Loss = has a GRCh38-only tensor but none on d9 (d46).
- "if in / if out" = d9 loss among alleles the model puts in / out of d9.
- "share predicted" = the fraction of the d9 loss that falls on predicted-in alleles.

The HP >= 7 sub-strata and the HCC1395 rows are in `tables.md`. HCC1395 has only d46 sets, so it has no measured loss.

| Sample | Platform | Role | Stratum | n | predicted in d9 | tensor linear / d9 / d46 | loss d9 | loss d46 | d9 loss if in / if out | share predicted |
|---|---|---|---|---|---|---|---|---|---|---|
| HG008 | PacBio | germline | SNV | 3,737,022 | 86.9 | 82.8 / 8.4 / 7.6 | 90.0 | 90.9 | 99.4 / 15.3 | 98.1 |
| HG008 | ONT | germline | SNV | 3,737,022 | 86.9 | 83.4 / 8.6 / 8.0 | 89.8 | 90.6 | 99.3 / 15.0 | 98.1 |
| HG008 | Illumina | germline | SNV | 3,737,022 | 86.9 | 82.4 / 8.1 / 7.5 | 90.2 | 91.0 | 99.5 / 15.3 | 98.1 |
| HG008 | PacBio | germline | INDEL | 964,490 | 76.1 | 75.7 / 15.0 / 14.9 | 82.0 | 82.0 | 86.4 / 63.5 | 84.9 |
| HG008 | ONT | germline | INDEL | 964,490 | 76.1 | 64.4 / 15.9 / 14.7 | 81.9 | 83.1 | 86.9 / 61.5 | 85.5 |
| HG008 | Illumina | germline | INDEL | 964,490 | 76.1 | 70.6 / 5.6 / 5.4 | 92.8 | 93.0 | 99.5 / 61.0 | 88.6 |
| HG008 | PacBio | somatic | SNV | 8,543 | 2.4 | 98.4 / 95.0 / 95.2 | 3.7 | 3.5 | 84.4 / 1.8 | 51.8 |
| HG008 | ONT | somatic | SNV | 8,543 | 2.4 | 98.7 / 95.4 / 95.6 | 3.5 | 3.2 | 84.6 / 1.7 | 53.5 |
| HG008 | Illumina | somatic | SNV | 8,543 | 2.4 | 98.6 / 95.4 / 95.8 | 3.3 | 2.9 | 73.4 / 1.6 | 52.1 |
| HG008 | PacBio | somatic | INDEL | 7,933 | 25.6 | 84.3 / 48.8 / 50.6 | 45.2 | 43.2 | 71.8 / 36.0 | 40.7 |
| HG008 | ONT | somatic | INDEL | 7,933 | 25.6 | 58.0 / 43.6 / 43.7 | 37.5 | 35.6 | 69.7 / 28.6 | 40.4 |
| HG008 | Illumina | somatic | INDEL | 7,933 | 25.6 | 80.1 / 42.9 / 45.9 | 49.2 | 45.8 | 95.6 / 31.9 | 52.7 |
| COLO829 | PacBio_fiberseq | germline | SNV | 3,729,123 | 87.2 | 87.2 / 9.2 / 8.4 | 89.5 | 90.5 | 99.1 / 14.7 | 98.1 |
| COLO829 | ONT | germline | SNV | 3,729,123 | 87.2 | 88.3 / 9.7 / 8.9 | 89.1 | 90.0 | 98.8 / 14.3 | 98.2 |
| COLO829 | Illumina | germline | SNV | 3,729,123 | 87.2 | 86.5 / 8.7 / 8.0 | 90.0 | 90.8 | 99.5 / 14.4 | 98.2 |
| COLO829 | PacBio_fiberseq | germline | INDEL | 977,293 | 75.6 | 78.7 / 15.8 / 15.8 | 81.7 | 81.6 | 86.9 / 60.9 | 85.1 |
| COLO829 | ONT | germline | INDEL | 977,293 | 75.6 | 65.2 / 18.1 / 17.2 | 79.5 | 80.5 | 84.6 / 58.8 | 85.5 |
| COLO829 | Illumina | germline | INDEL | 977,293 | 75.6 | 72.6 / 5.9 / 5.7 | 92.6 | 92.9 | 99.5 / 59.8 | 88.9 |
| COLO829 | PacBio_fiberseq | somatic | SNV | 42,035 | 0.2 | 90.3 / 89.8 / 89.8 | 0.6 | 0.7 | 77.5 / 0.5 | 22.5 |
| COLO829 | ONT | somatic | SNV | 42,035 | 0.2 | 89.6 / 89.3 / 89.2 | 0.5 | 0.6 | 73.9 / 0.4 | 25.0 |
| COLO829 | Illumina | somatic | SNV | 42,035 | 0.2 | 90.9 / 90.1 / 90.1 | 1.0 | 1.0 | 72.2 / 0.9 | 14.4 |
| COLO829 | PacBio_fiberseq | somatic | INDEL | 1,912 | 22.3 | 96.3 / 63.3 / 65.9 | 34.8 | 32.1 | 72.8 / 24.1 | 46.3 |
| COLO829 | ONT | somatic | INDEL | 1,912 | 22.3 | 75.4 / 58.0 / 58.4 | 30.8 | 29.8 | 68.1 / 21.5 | 44.1 |
| COLO829 | Illumina | somatic | INDEL | 1,912 | 22.3 | 91.4 / 57.1 / 60.4 | 38.1 | 34.9 | 93.5 / 22.4 | 54.4 |

**6d. Measured loss by HPRC v1.1 carrier count.** Each cell gives n alleles with a GRCh38-only tensor, then the loss in
% on d9 / d46 (tables.md G; long-read sets).

HG008T PacBio:

| carriers (of 89) | germline SNV | germline INDEL | somatic SNV | somatic INDEL |
|---|---|---|---|---|
| 0 | 95,344; 17.4 / 19.2 | 71,469; 81.9 / 82.0 | 8,047; 1.4 / 1.4 | 2,826; 34.0 / 33.8 |
| 1 | 37,464; 10.0 / 11.4 | 9,196; 45.1 / 44.6 | 53; 18.9 / 18.9 | 491; 33.0 / 33.4 |
| 2-4 | 88,682; 12.5 / 18.6 | 25,463; 44.1 / 46.6 | 67; 23.9 / 23.9 | 945; 39.6 / 40.3 |
| 5-8 | 123,841; 17.2 / 50.9 | 36,484; 45.6 / 61.9 | 48; 31.2 / 39.6 | 714; 41.5 / 43.7 |
| 9-17 | 319,645; 98.5 / 91.0 | 92,746; 88.9 / 83.6 | 61; 80.3 / 54.1 | 861; 71.8 / 56.9 |
| 18-44 | 969,507; 99.2 / 99.3 | 252,475; 87.0 / 86.7 | 72; 80.6 / 75.0 | 752; 72.6 / 69.5 |
| 45-89 | 1,458,896; 99.7 / 99.7 | 242,530; 84.9 / 84.6 | 59; 93.2 / 84.8 | 100; 66.0 / 64.0 |

COLO829T fiber-seq:

| carriers (of 89) | germline SNV | germline INDEL | somatic SNV | somatic INDEL |
|---|---|---|---|---|
| 0 | 103,510; 17.3 / 18.8 | 76,229; 79.6 / 79.9 | 37,679; 0.5 / 0.6 | 932; 26.5 / 26.7 |
| 1 | 40,261; 10.1 / 11.3 | 10,314; 39.8 / 40.2 | 94; 1.1 / 0.0 | 126; 19.1 / 15.9 |
| 2-4 | 92,988; 12.4 / 19.2 | 27,875; 42.2 / 44.7 | 68; 4.4 / 5.9 | 205; 19.0 / 17.6 |
| 5-8 | 133,424; 15.7 / 49.4 | 39,538; 43.5 / 60.0 | 33; 6.1 / 9.1 | 171; 20.5 / 22.8 |
| 9-17 | 347,092; 98.2 / 90.6 | 99,772; 89.5 / 83.7 | 29; 79.3 / 48.3 | 208; 73.1 / 54.3 |
| 18-44 | 1,033,674; 98.9 / 99.1 | 266,450; 87.6 / 87.2 | 35; 71.4 / 71.4 | 162; 72.2 / 66.7 |
| 45-89 | 1,502,253; 99.5 / 99.6 | 248,995; 85.0 / 84.5 | 7; 100.0 / 100.0 | 38; 73.7 / 68.4 |

The other platforms are in `calibration_bins.tsv`.

### 7. Relatedness check (R3.m3) (tables.md J and `relatedness.tsv`)

The test set is each normal's rare germline SNVs (PASS, in BED, 1-2 of 89 HPRC carriers). The table counts how many of
them each HPRC v1.1 individual carries (either haplotype) and shows the top three per sample:

| Sample | rare SNVs | HPRC individual | superpopulation | shared | share % | ratio to superpopulation median |
|---|---|---|---|---|---|---|
| HG008 | 95,232 | HG00733 | AMR | 7,674 | 8.06 | 1.23 |
| HG008 | 95,232 | HG01106 | AMR | 7,444 | 7.82 | 1.19 |
| HG008 | 95,232 | HG01175 | AMR | 7,189 | 7.55 | 1.15 |
| COLO829 | 92,076 | HG00733 | AMR | 7,126 | 7.74 | 1.19 |
| COLO829 | 92,076 | HG01175 | AMR | 6,971 | 7.57 | 1.16 |
| COLO829 | 92,076 | HG01123 | AMR | 6,969 | 7.57 | 1.16 |
| HCC1395 | 74,604 | HG01175 | AMR | 6,112 | 8.19 | 1.23 |
| HCC1395 | 74,604 | HG00733 | AMR | 5,944 | 7.97 | 1.2 |
| HCC1395 | 74,604 | HG00741 | AMR | 5,599 | 7.50 | 1.13 |

Highest ratio in each other superpopulation (`relatedness.tsv`):

- AFR: 2.2 (HG008, NA21309, share 3.10 %), 2.01 (COLO829, HG02145, 2.54 %) and 2.3 (HCC1395, NA20129, 2.84 %).
- EAS: at most 1.26.
- SAS and CHM13 are single individuals, so their ratio is 1.0 by construction.

## Interpretation

1. **Size vs threshold (R1.2).**
   - At a relative frequency filter, the haplotype count barely changes either what the graph removes or what it
     absorbs. From k = 16 to 89 at f = 10 %, germline SNV coverage moves by +2.4 to +2.5 points and somatic INDEL
     absorption by -0.4 to 0 points (section 2a). This is the in-silico counterpart of the measured result: d46
     (464 haplotypes, about 10 %) and d9 (90 haplotypes, about 10 %) lose about the same (section 6c, loss d9 vs d46).
   - The haplotype count matters only when the filter is absolute, as in an unfiltered graph (d = 1). There germline
     coverage saturates early (94.0-94.9 % of germline SNVs at k = 32) while somatic INDEL absorption keeps growing
     (section 2b).
   - So "a larger graph absorbs more somatic variants" holds for an unfiltered graph. It does not hold at the
     frequency filter Pansoma uses.
2. **Where the threshold sits.** On the full v1.1 set, moving from d = 9 to d = 1 has a lopsided cost:
   - germline SNV coverage gains 8.7-9.1 points;
   - somatic INDEL absorption gains 26.9 (COLO829) and 33.4 (HG008) points;
   - moving to d = 18 instead costs 11.1-11.2 points of germline SNV coverage, for 11.3-13.0 points less HG008 /
     COLO829 somatic INDEL absorption (section 2c).

   The d9 graph keeps 6.87 M of the 21.08 M SNV alleles (section 5) and still contains 86.9-88.1 % of the donors'
   germline SNVs.
3. **What gets absorbed.** Mostly repeat-length INDELs: HP >= 7 and STR, at 27-42 % at d = 9 and 67-71 % at d = 1 for
   HG008 / COLO829 (section 3). Other INDELs are at 3.0-6.4 % (d = 9) and 8.5-8.9 % (d = 1), and somatic SNVs at 0.2-2.4 % at d = 9 (section 2c).
   The HCC1395 truth has few HP >= 7 INDELs (234 of 1,559), so its absorption is low.
4. **Ancestry (composition) matters more than count.** For all three donors the 42 non-AFR haplotypes (32 of them AMR)
   contain more of the donor's germline alleles than 42 AFR haplotypes:
   - +3.3 to +3.5 points of germline SNV at d = 1;
   - +6.3 to +6.8 points at f = 10 %.

   The AFR haplotypes absorb more somatic INDELs (section 4). A graph whose haplotypes match the patient's ancestry
   suppresses more germline variants at the same size. The relatedness check fits this: the AMR individuals share the
   most rare alleles with all three donors.
5. **Calibration.**
   - The d9 filter behaves as ">= 9 of the 89 v1.1 haplotypes" (section 6a).
   - The model predicts the d9 content of the sample alleles to within 1.2 points for the SNV / INDEL rows of 6b (up
     to 1.8 points in the INDEL sub-strata of tables.md I).
   - For germline SNVs it explains the measured loss on identical reads almost completely: alleles predicted in d9 lose
     98.8-99.5 % of their candidates, and 98.1-98.2 % of all germline SNV losses are predicted alleles.
   - It explains only about half of the somatic SNV loss on HG008 (51.8-53.5 %), 14.4-25.0 % on COLO829, and 40-54 % of
     the somatic INDEL loss. Somatic INDELs absent from HPRC v1.1 (0 carriers) lose 34.0 % (HG008 PacBio) and 26.5 %
     (COLO829 fiber-seq) of their candidates on d9, and about the same on d46 (33.8 % / 26.7 %).
   - So about half of the graph-arm INDEL loss is not "the allele is in the graph". It is something else (see Caveats
     and Status). This matters for R1.4: graph composition is only part of the somatic INDEL story.
6. **d46 on the v1.1 axis.** d46 removes 49.4-50.9 % of the germline SNVs that have 5-8 v1.1 carriers, against
   15.7-17.2 % for d9. It removes 90.6-91.0 % of those with 9-17 carriers (d9: 98.2-98.5 %) (section 6d). The d46
   threshold is set on v2.1 frequencies, so its step is blurred on the v1.1 carrier axis. Net germline SNV loss is
   slightly higher on d46 (90.0-91.0 % vs 89.1-90.2 %).
7. **Relatedness (R3.m3).** No HPRC v1.1 individual is a close relative of the HG008, COLO829 or HCC1395 donor. The
   maximum share of rare alleles is 7.74-8.19 %, against about 50 % expected for a first-degree relative. Within each
   superpopulation, every individual is at most 1.23x (AMR), 1.26x (EAS) or 2.3x (AFR, at shares of only 2.5-3.1 %) the median.
   The modest AFR outliers are consistent with admixture, but this was not tested.

Suggested text for the open "curves" placeholder of the R1.2 response (draft, numbers from this README):

> We quantified both axes in silico with random subsets of 2-89 HPRC v1.1 haplotypes and carrier thresholds from one
> haplotype to 20 %. At a 10 % threshold, increasing the subset from 16 to 89 haplotypes changed the fraction of
> matched-normal germline SNVs contained in the graph from 84.5 % to 86.9 % (HG008) and of somatic INDELs from 25.9 % to
> 25.6 %. Without a frequency filter, germline coverage saturated (94.0 % at 32 haplotypes, 96.0 % at 89) while somatic
> INDEL absorption kept rising (49.0 % to 59.0 %). This growth was concentrated in homopolymer (>= 7 bp) and STR indels
> (70.9 % and 66.9 % absorbed, versus 8.5 % for other indels). The model reproduced the real graph: 96-98.5 % of HPRC SNV
> alleles carried by at least 9 haplotypes, and at most 7.2 % of those carried by fewer, are in the v1.1 graph. Predicted
> alleles account for 98 % of the germline SNV candidates the graph removes on identical reads.

## Caveats

- **"In the graph" is not "absorbed".**
  - The calibration (6c, 6d) shows the model is close for germline SNVs and partial for somatic alleles.
  - About half of the measured somatic INDEL loss falls on alleles the model puts outside d9, including alleles with
    no HPRC carrier. Plausible causes are other alleles of the same repeat in the graph, representation differences,
    or node-space candidate generation. **None was tested here.**
  - Germline INDELs with 0 v1.1 carriers lose 79.6-81.9 % of their candidates on d9. This suggests many are in the graph
    in a form the exact key match does not see, but E7 cannot separate the two explanations.
- **Allele matching** is exact after left-alignment (`bcftools norm -c w`). The design proposed enumerating equivalent
  placements; this was not done.
  - COMPLEX alleles match only exactly.
  - MNP-like truth alleles need all their primitives on one haplotype and are excluded from the SNV / INDEL strata.
  - The "0 carriers" class therefore mixes truly novel alleles with alleles in another representation.
- **vcfwave GT repair.** For 6,947,944 decomposed primitives, the repaired carrier count:
  - equals the wave INFO AC for 6,435,168 (92.6 %);
  - is above it for 415,452;
  - is below it for 97,324.

  The repaired total is 68.84 M carriers, against 64.33 M AC and 56.05 M in the wave GT. 10,219 raw alleles
  (13,240 carriers) could not be decomposed and are missing.
- **d9 VCF presence is a lower bound.** The d9 deconstruct VCF is not vcfwave-decomposed, which is why only 96-98.5 % of
  SNV alleles with >= 9 carriers are found (6a, 6b).
- **The model covers HPRC v1.1 only** (89 haplotypes; GRCh38 is always the backbone). It cannot extrapolate beyond
  k = 89. d46 is compared only through measured loss, not modelled.
- **Population pools are small.** AFR has 46 haplotypes. Non-AFR has 42 (AMR 32, EAS 8, SAS 2), so "non-AFR" mostly
  means AMR. There is no European pool in v1.1.
- **HCC1395 has no linear or d9 tensor sets**, so it has no measured loss (empty columns in tables.md F).
- **Scope differs from the R1.2 response table** (here chr1-22, PASS within the BED). For example, HG008T PacBio somatic
  SNV tensor rates are 98.4 / 95.0 / 95.2 % here vs 0.983 / 0.947 / 0.949 there.
- **Relatedness** was checked against HPRC v1.1 only, for three donors. A first-degree relative sharing "about half" is
  the script's premise, not a measured control.

## Status

**DONE (2026-10-13)**

- s1-s6 (Slurm jobs above, all COMPLETED).
- s7 `tables.md` / `summary.json` (collect.sh, 22:55).

They cover the haplotype-count x threshold curves for HG008, COLO829 and HCC1395, the ancestry pools, graph size,
calibration against d9 / d46 on HG008T and COLO829T x 3 platforms, and the R3.m3 relatedness check.

**PENDING**

1. **HPRC v2.1 axis** (k = 128, 256, 462; calibration of k = 462, f = 10 % against d46). This needs the HPRC v2.1
   full-graph vcfbub / wave VCF. The download is on the user-decision list (`response_to_reviewers_detailed.md` §7, item
   4). Nothing has been downloaded or submitted.
2. **E7b** (optional, P2): one measured point on a v2.1 d23 graph. It needs the v2.1 full GBZ, new indexes and
   alignment, and is not started (user decision).
3. **Figure.** `curves.tsv` is figure-ready (`WORK/hist_draws.npy` keeps the full per-draw histograms), but no plot has
   been made yet.
4. **The non-composition half of the graph-arm INDEL loss** (6c / 6d). This is outside E7's design. It is a candidate
   follow-up for R1.4 / E5b; no job is planned.
5. **Relatedness** against v2.1 and for the other CASTLE lines (HCC1954, H1437, H2009): not done. The s2 inputs cover
   only HG008, COLO829 and HCC1395.

No E7 job is running, queued or held, so nothing here waits on the lab-cluster cap. The comparator re-runs and the CHM13
work dropped today do not affect E7.

## Files

This folder:

- Scripts: `common.py`, `s1`-`s7` `*.py`, `collect.sh`.
- Results: `curves.tsv`, `graph_size.tsv`, `calibration.tsv`, `calibration_bins.tsv`, `d9_retention.tsv`,
  `d9_content.tsv`, `relatedness.tsv`, `tables.md`, `summary.json`.

Job scripts and logs: `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E7_graph_composition/`.

Data directory `/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/`:

- `hprc_alleles/` (s1; `chr21test.*` is the 10:17 test run);
- `sample_alleles/` (s2);
- `d9_alleles/` (s5);
- `joined.npz` and `hist_draws.npy` (s3).
