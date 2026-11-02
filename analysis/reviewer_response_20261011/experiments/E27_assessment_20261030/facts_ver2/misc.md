# Fact sheet: R2.2, R3.2, R3.3–R3.5 ("22–66"), R3.7, R3.m2 (2026-10-30)

Legend: **EXISTS** = number is in a committed/working-tree file (path given); **COMPUTED NOW** = computed today
(read-only, scripts + outputs in `/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/`); **MISSING** = not
available anywhere I could find. No tracked file was edited. "Letter" = the current `points/*.md` Responses (identical to
the 10-30 assembled copy `.../scratchpad/assemble/points/`).

Short paths: `RR` = `/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011`; `EXP` = `RR/experiments`;
`MISC` = `/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc`; `WORK` =
`/scratch/jshen/data/reviewer_experiments_20261013`.

---

## 1. R2.2: truth sources, and the H1437 numbers

### 1a. Truth-source table (all EXISTS unless marked)

| tumor | somatic truth (name / version) | how generated | validation | region | counts (records; in-BED chr1–22 / chr1 in BED) | germline truth (labels) | known limitations | source |
|---|---|---|---|---|---|---|---|---|
| HG008T (PDAC) | GIAB HG008-T somatic smvar **draft v0.2** (pre-release; not on public FTP; public = v0.3, 2026-04-25, 18,946 variants) | tumor (HG008-T) vs normal (HG008-N **v6.2**) assembly comparison (v0.2 VCF header: `HG008Nv62SOMATICVARIANT`, vcfeval dipHG8N vs dipHG8T) | v0.3 README: repun with PacBio HiFi / ONT / Illumina-Element read AF/DP + 9 community call sets; **whether v0.2 had the same read curation: MISSING (unverified)** | `_all.bed` 41,835 regions, 2,809,425,965 bp (2.81 Gb) | 18,952 (9,609 SNV, 9,343 INDEL); SNV 8,543 / 697, INDEL 7,982 / 588 | HG008-N dipcall (`HG008N_GRCh38_dipcall.dip.vcf.gz`) | draft; INDELs 74.9 % in HP ≥ 7 (chr1–22, 8,496); ~44 % of somatic variants change representation on GRCh38 (v0.3 README); 209 truth alleles = an HG008-N germline allele | E8b G1–G6 (`EXP/E8b_literature_facts/README.md`); `EXP/E22.../graph_overlap.tsv` (HG008 bed / bed_chr1); `analysis/indel_repeat_context_predictions_20261006/README.md` |
| COLO829T (melanoma) | SMaHT COLO829 reference set (Ha et al. 2025, bioRxiv 10.1101/2025.10.13.681545) = SMaHT `COLO829BLT50_reference_set` exactly (memory check 10-29) | Illumina candidates from Mutect2, Strelka2, RUFUS, VarNet (COLO829 200× vs COLO829BL 300×) | PacBio (300× / 180×): ALT in tumor, absent in normal | SMaHT easy ∪ difficult ∪ extreme v2 = whole genome ("scored genome-wide") | 44,005 SNV + 2,059 INDEL (961 ins, 1,098 del); SNV 42,035 / 2,443, INDEL 1,912 / 131 | COLO829BL DSA dipcall (`dipcall_hg38.dip.vcf.gz`; dip.bed 2.656 Gb chr1–22) | Illumina-derived candidates (short-read callable); includes subclonal SNVs: 4,616 SNVs with VAF_Ill < ~0.1 that the SMaHT "challenge" version drops (450 of 2,443 chr1); 55 truth INDELs where a COLO829BL haplotype already carries the ALT over the repeat | E8b S1–S5; `pansoma_v2_tensors/COLO829T_truth/README.txt`; memory `colo829t-truth-vs-smaht-release.md`; `analysis/graph_absorbed_somatic_colo829t_20261005/README.md` |
| HCC1395 (TNBC) | SEQC2 **v1.2.1** HighConf + MedConf sSNV + sINDEL | 21 replicates × 3 aligners × 6 callers = 378 Illumina call sets scored by SomaticSeq / NeuSomatic | Ion Torrent, AmpliSeq, WES, PacBio | `SEQC2_High-Confidence_Regions_v1.2.bed`, 2,485,467,324 bp | 41,072; SNV 39,447 / 3,440, INDEL 1,564 / 124 | HCC1395BL DeepVariant hybrid (DeepSomatic training VCF, FILTER GERMLINE), QUAL ≥ 10 | Illumina-analysable regions only (paper's own limitation); few HP ≥ 7 INDELs (15.5 %) → rule (a) | E8b Q1–Q3; `/scratch/jshen/data/CASTLE/README.txt`; graph_overlap.tsv |
| HCC1954 (breast) | CASTLE = DeepSomatic `DeepSomatic_HCC1395-model_benchmark`, `vcfs/somatic_only` PASS | DeepSomatic model calls on tumor–normal (multi-technology) | DeepSomatic paper: "strongly biased toward the single-cancer trained models … biased toward DeepSomatic in general" (D8) | `..._high-confidence-regions_DeepSomatic-HCC1395-model.bed`, ~2.51 Gb | 26,578; SNV 20,846 / 1,450, INDEL 2,272 / 169 | DeepVariant hybrid of the normal, QUAL ≥ 10 | DeepSomatic-derived; in both comparators' training data; 2.39 % common alleles (IG/TCR, chr7:90–91 Mb) but no normal LOH | CASTLE README; E8b D4/D8; E22 T1/T2 |
| H1437 (lung adeno) | same CASTLE benchmark | same | same | same, ~2.51 Gb | 87,410; SNV 73,623 / 3,941, INDEL 4,439 / 158 | same | **~30 % germline-origin SNVs (normal-LOH chr9q)**, see 1b | E22 |
| H2009 (lung adeno) | same CASTLE benchmark | same | same | same, ~2.51 Gb | 168,397; SNV 144,177 / 11,106, INDEL 5,843 / 378 | same | ~5 % germline-origin (normal-LOH chr7q etc.) | E22 |
| HapMap mixture (Supplement only) | SMaHT HAPMAP6 (83.5 % HG005 + 5 minor donors); pangenome-graph truth **v1.6** of Kong et al. (v1.6_wo_HG005) | minor donors' germline alleles absent from HG005, built from individual assemblies in pangenome graphs | technology-agnostic (Kong et al.) | `SMHTHAPMAP6_GRCh38_v1.6_somatic_benchmark_autosomes.bed`, 2,592,293,110 bp (COMPUTED NOW, awk sum) | 7,870,692 (6,037,703 SNV + 1,832,989 INDEL); evaluated only at expected VAF 1–16.5 % | HG005 alleles (4,920,296 records) | germline alleles, not somatic; donors share alleles with the graph (R2.1) | E8b M1–M3; `pansoma_v2_tensors/HapMap_truth/README.txt`; E23 README |

Counts "in-BED chr1–22 / chr1" are E7/E22 normalized primitive alleles (`EXP/E22_truth_and_version_checks/graph_overlap.tsv`,
scopes `bed` / `bed_chr1`); vcfeval-level chr1 counts differ slightly (e.g. HG008T chr1 truth INDELs 635 in R1.6 vs 588
primitives here).

### 1b. H1437 numbers: exact definitions (all EXISTS; chr1:9–12 Mb split COMPUTED NOW)

Sources: `EXP/E22_truth_and_version_checks/{README.md, tables.md (T1–T4), graph_overlap.tsv, graph_overlap_training.tsv,
label_impact.tsv}`; per-allele table `WORK/E22_truth_and_version_checks/H1437_CASTLE.alleles.tsv.gz`.

| letter figure | where | numerator / denominator | region | definition |
|---|---|---|---|---|
| **24.5 %** | R1.m1 | 18,028 / 73,623 truth SNVs | FILTER PASS, **inside the CASTLE BED, chr1–22** (E22 check 3 `bed` scope) | allele carried by ≥ 46 of the 463 HPRC v2.1 haplotypes (= predicted d46 graph allele, E7 v2.1) |
| **27.8 %** | R2.2 (also R2.6, R3.8) | 22,982 / 82,529 truth SNVs (27.85 %) | **all** truth SNV alleles, genome-wide incl. outside the BED and chrX (in BED: 21,180 / 73,623 = 28.77 %) | gnomAD AF ≥ 1 % (exact allele) |
| **"about 30 %"** | R2.2 | 24,914 / 82,529 (30.19 %) truth SNVs lying in the two normal-LOH segments; independently, the common-allele mixture estimate (27.85 − 0.54) / (90.90 − 0.54) = 30.2 % (~24,942) | all truth SNVs | estimate, not a per-variant call ("a minority of the segment SNVs may be real somatic mutations", E22) |
| **97.5 %** | R2.2 | 22,397 / 22,975 **autosomal** common (gnomAD ≥ 1 %) truth SNVs (97.48 %) | genome-wide autosomes | inside the normal-LOH segments ≥ 3 Mb with normal het fraction < 0.2: chr1:9–12 Mb and chr9:68–139 Mb |
| 99 % / 91 % | R1.m1 | H1437 17,859 / 18,028; H2009 6,055 / 6,637 | in BED, chr1–22 | share of the ≥ 46/463 truth SNVs in **1-Mb bins** with normal het < 0.2 (not the ≥ 3 Mb segments) |
| **19 of 4,266** | R2.2 | 19 common (gnomAD ≥ 1 %) / 4,266 chr1 truth SNVs (0.45 %) | **chr1, all records** (in BED: 17 / 3,941 per graph_overlap `bed_chr1`) | none of the 19 is in a normal-LOH bin (T3) |
| "population filtering tags all but two" / "cap … about 75 %" | R1.m1 | H1437: **18,028 / 18,028 tagged**; H2009: 6,635 / 6,637; cap = 1 − 0.2449 = 75.5 % | in BED, chr1–22 | gnomAD or CoLoRSdb AF ≥ 1e-4, exact allele |

COMPUTED NOW (`H1437_CASTLE.alleles.tsv.gz`, python one-off on login node): the chr1:9–12 Mb segment holds only **30**
truth SNVs (29 in BED), **0** of them common. All 22,397 common SNVs in segments are in **chr9:68–139 Mb** (24,884 SNVs
there in total; 24,884 + 30 = 24,914). The 19 chr1 common SNVs lie at chr1:51.7–240.0 Mb, all in bins with normal het
≥ 0.2 (18 in het 0.4–1.0, 1 in 0.2–0.4); 17 are in BED.

**Do the 97.5 % and the chr1 19 / 4,266 apply to the same set?** No. 97.5 % is over the 22,975 autosomal common truth
SNVs genome-wide (all in chr9q); 19 / 4,266 is over chr1 truth SNVs (all records), whose common alleles all lie outside
the LOH bins. The chr1:9–12 Mb LOH segment contributes 30 truth SNVs to the "about 30 %" (24,914) but no common allele
to the 97.5 %.

Other lines (EXISTS, E22): H2009 8,141 / 161,917 common (5.03 %), 87.23 % of autosomal common in segments
(chr5:47–51, chr7:55–62, chr7:64–67, chr7:122–143 Mb), segments hold 9,824 (6.07 %), mixture 5.0 %; HCC1954 574 / 24,033
(2.39 %), 0 in normal-LOH segments, mixture 2.1 % (~500), common SNVs clustered at chr7:90–91 Mb and IGH/TRA/IGK;
SEQC2 HCC1395 212 / 39,447 (0.54 %).

---

## 2. R3.2: candidate retention, graph representation, final FPs

### 2a. Somatic and germline candidate retention by graph representation (COMPUTED NOW)

Script `MISC/r32_germline_split.py` (Slurm 462872, miao, ~10 min); output `MISC/r32_germline_split.tsv` (336 rows: also
`bed_chr1` scope and the `v11_vcf` definition) and `MISC/r32_germline_split_compact.md`. Same inputs/definitions as E7
`s4_calibration.py`: truth alleles = PASS, in BED, chr1–22 (E7 `joined.npz`, v1.1 and v2.1, row order checked identical);
"candidate" = status `tensor_representative` in `<set>/tensors/{somatic,germline}.recall.tsv`. Regression: HG008T PacBio
germline SNV linear / d9 / d46 = 0.828 / 0.084 / 0.076 and predicted-in-d9 0.8694 / -d46 0.8733 equal E7
`calibration.tsv` / `v2.1/calibration.tsv`.

"Represented in the graph": **v1.1** = ≥ 9 of 89 HPRC v1.1 haplotypes carry the allele (E7 "predicted in d9"; the real
d9 deconstruct-VCF membership `v11_vcf` gives the same picture, in the TSV); **v2.1** = ≥ 46 of 463 HPRC v2.1 haplotypes
(E7 "predicted in d46"). **No d46 deconstruct VCF exists (HPRC publishes only GBZ/indexes), so for v2.1 the
representation is carrier-count based only** (E7 README "There is no d46 counterpart of the d9 deconstruct VCF").

Retention = % of truth alleles with a representative tensor. v1.1 columns compare linear → d9; v2.1 columns linear → d46.

| tumor / platform | class | all alleles: n; linear / d9 / d46 | v1.1-represented: n; linear → d9 | not v1.1-repr.: n; linear → d9 | v2.1-represented: n; linear → d46 | not v2.1-repr.: n; linear → d46 |
|---|---|---|---|---|---|---|
| HG008T PacBio | somatic SNV | 8,543; 98.4 / 95.0 / 95.2 | 205; 93.7 → 15.1 | 8,338; 98.5 → 96.9 | 143; 96.5 → 18.2 | 8,400; 98.4 → 96.5 |
| | somatic INDEL | 7,933; 84.3 / 48.8 / 50.6 | 2,030; 84.4 → 26.3 | 5,903; 84.3 → 56.5 | 1,522; 85.9 → 25.0 | 6,411; 83.9 → 56.7 |
| | germline SNV | 3,737,022; 82.8 / 8.4 / 7.6 | 3,249,012; 84.6 → 0.6 | 488,010; 70.8 → 60.4 | 3,263,667; 84.6 → 0.4 | 473,355; 70.1 → 57.7 |
| | germline INDEL | 964,490; 75.7 / 15.0 / 14.9 | 734,132; 80.1 → 12.0 | 230,358; 61.9 → 24.6 | 750,735; 80.0 → 11.6 | 213,755; 60.7 → 26.4 |
| HG008T ONT | somatic SNV | 8,543; 98.7 / 95.4 / 95.6 | 205; 91.7 → 14.1 | 8,338; 98.8 → 97.4 | 143; 96.5 → 13.3 | 8,400; 98.7 → 97.0 |
| | somatic INDEL | 7,933; 58.0 / 43.6 / 43.7 | 2,030; 49.3 → 27.4 | 5,903; 61.0 → 49.2 | 1,522; 52.1 → 26.2 | 6,411; 59.4 → 47.9 |
| | germline SNV | 3,737,022; 83.4 / 8.6 / 8.0 | 3,249,012; 85.1 → 0.6 | 488,010; 72.0 → 61.6 | 3,263,667; 85.2 → 0.4 | 473,355; 71.4 → 59.9 |
| | germline INDEL | 964,490; 64.4 / 15.9 / 14.7 | 734,132; 68.2 → 13.7 | 230,358; 52.3 → 23.0 | 750,735; 68.1 → 11.9 | 213,755; 51.5 → 24.3 |
| HG008T Illumina | somatic SNV | 8,543; 98.5 / 95.4 / 95.8 | 205; 97.1 → 25.9 | 8,338; 98.6 → 97.1 | 143; 97.9 → 23.1 | 8,400; 98.6 → 97.1 |
| | somatic INDEL | 7,933; 80.1 / 42.9 / 45.9 | 2,030; 84.9 → 4.0 | 5,903; 78.5 → 56.2 | 1,522; 84.4 → 1.6 | 6,411; 79.1 → 56.4 |
| | germline SNV | 3,737,022; 82.4 / 8.1 / 7.5 | 3,249,012; 84.3 → 0.4 | 488,010; 69.7 → 59.3 | 3,263,667; 84.3 → 0.3 | 473,355; 69.3 → 57.1 |
| | germline INDEL | 964,490; 70.6 / 5.6 / 5.4 | 734,132; 76.7 → 0.7 | 230,358; 51.4 → 21.4 | 750,735; 76.4 → 0.4 | 213,755; 50.3 → 23.2 |
| COLO829T Fiber-seq | somatic SNV | 42,035; 90.3 / 89.8 / 89.8 | 80; 88.8 → 20.0 | 41,955; 90.3 → 90.0 | 63; 88.9 → 17.5 | 41,972; 90.3 → 89.9 |
| | somatic INDEL | 1,912; 96.3 / 63.3 / 65.9 | 426; 95.8 → 26.5 | 1,486; 96.5 → 73.8 | 340; 94.1 → 23.8 | 1,572; 96.8 → 75.0 |
| | germline SNV | 3,729,123; 87.2 / 9.2 / 8.4 | 3,253,239; 88.6 → 0.8 | 475,884; 77.8 → 66.7 | 3,270,210; 88.6 → 0.5 | 458,913; 77.3 → 64.2 |
| | germline INDEL | 977,293; 78.7 / 15.8 / 15.8 | 738,833; 83.3 → 12.0 | 238,460; 64.6 → 27.6 | 754,006; 83.3 → 11.8 | 223,287; 63.1 → 29.6 |
| COLO829T ONT | somatic SNV | 42,035; 89.6 / 89.3 / 89.2 | 80; 86.2 → 22.5 | 41,955; 89.6 → 89.4 | 63; 85.7 → 19.0 | 41,972; 89.6 → 89.4 |
| | somatic INDEL | 1,912; 75.4 / 57.9 / 58.4 | 426; 67.6 → 34.3 | 1,486; 77.6 → 64.7 | 340; 70.9 → 31.5 | 1,572; 76.3 → 64.2 |
| | germline SNV | 3,729,123; 88.3 / 9.7 / 8.9 | 3,253,239; 89.6 → 1.1 | 475,884; 79.6 → 68.6 | 3,270,210; 89.6 → 0.8 | 458,913; 79.2 → 66.8 |
| | germline INDEL | 977,293; 65.2 / 18.1 / 17.3 | 738,833; 69.3 → 15.9 | 238,460; 52.6 → 25.2 | 754,006; 69.4 → 14.5 | 223,287; 51.3 → 26.4 |
| COLO829T Illumina | somatic SNV | 42,035; 90.9 / 90.1 / 90.1 | 80; 98.8 → 27.5 | 41,955; 90.9 → 90.2 | 63; 95.2 → 25.4 | 41,972; 90.9 → 90.2 |
| | somatic INDEL | 1,912; 91.4 / 57.1 / 60.4 | 426; 90.8 → 6.3 | 1,486; 91.5 → 71.7 | 340; 89.7 → 0.6 | 1,572; 91.7 → 73.3 |
| | germline SNV | 3,729,123; 86.5 / 8.7 / 8.0 | 3,253,239; 88.0 → 0.4 | 475,884; 75.8 → 65.2 | 3,270,210; 88.0 → 0.3 | 458,913; 75.4 → 62.8 |
| | germline INDEL | 977,293; 72.6 / 5.9 / 5.7 | 738,833; 79.5 → 0.7 | 238,460; 51.2 → 22.0 | 754,006; 79.3 → 0.4 | 223,287; 50.0 → 23.8 |

Readings: (i) graph-represented germline SNVs are almost fully removed (0.3–1.1 % retained on the graph that represents
them); unrepresented germline SNVs keep 57–69 % (vs 69–80 % on linear). (ii) Graph-represented germline INDELs are
removed on Illumina (0.4–0.7 % retained) but **not on long reads (11.6–15.9 % retained)**. (iii) Somatic INDELs lose
candidates even when **not** represented (e.g. HG008T PacBio 83.9 → 56.7 %, COLO829T Illumina 91.7 → 73.3 % on v2.1);
represented somatic INDELs fall to 0.6–34 %. (iv) 19.2–25.6 % of HG008T and 17.8–22.3 % of COLO829T in-BED somatic
INDELs are graph-represented (1,522–2,030 / 7,933; 340–426 / 1,912). The R3.2 letter's 87 % / 75 % (germline SNV /
INDEL) and 0.2–2.5 % / 21–25 % (somatic SNV / INDEL) use the real d9 VCF (`v11_vcf`; EXISTS, `EXP/E7_graph_composition/d9_content.tsv`) and are reproduced here (e.g. HG008 germline SNV 3,250,159 / 3,737,022 = 0.8697).

### 2b. Background (label 0 = in somatic ∩ germline BED, no truth allele) candidates (EXISTS; re-collected)

`analysis/graph_ablation_tensors_20261005/collect.py` rerun to `MISC/ablation_summary_20261030.json` (reads the 18
sets' `labels.manifest.json`); identical to the R3.2 evidence where it lists them.

| set | SNV non: linear → d9 → d46 | INDEL non: linear → d9 → d46 |
|---|---|---|
| HG008T Illumina | 2,114,925 → 1,876,969 → 1,852,495 | 141,083 → 31,433 → 23,240 |
| HG008T PacBio | 188,910 → 73,454 → 77,478 | 1,531,401 → 1,518,657 → 1,522,059 |
| HG008T ONT | 339,711 → 197,575 → 206,573 | 1,951,739 → 1,976,637 → 1,967,176 |
| COLO829T Illumina | 1,285,779 → 934,450 → 920,461 | 137,944 → 34,660 → 26,088 |
| COLO829T Fiber-seq | 250,913 → 120,489 → 162,510 | 1,272,815 → 1,310,694 → 1,321,934 |
| COLO829T ONT | 577,417 → 410,794 → 546,170 | 3,714,293 → 3,822,962 → 3,779,537 |

### 2c. Final FPs after population filtering, chr1, in BED, split by graph representation (SNV EXISTS in E4b; v2.1 split and INDEL COMPUTED NOW)

Script `MISC/r32_fp_split.py` (login node, 6 s) on E4b's per-FP tables
(`WORK/call_level/E4b/<call set>.eval_bed.fp.v2.tsv.gz`); output `MISC/r32_fp_split.tsv`, `MISC/r32_fp_split_compact.md`.
Call sets = E4b in-sample, seed 0, 7fc7396 recipe, **main threshold (validation recall 0.9)**, current INDEL PoN rule
(HP ≥ 8 AF ≥ 0.01, else 1e-4). Caveats: HG008T PacBio / ONT d9 are the 1-GPU recipe (E4b "(*)"); COLO829T d9
single-set runs are not local, so COLO829T has linear vs d46 only; **the INDEL rows are the 7fc7396 INDEL models at the
val-recall-0.9 threshold, not the final INDEL models (ppr8) of R1.6/R3.10.** E4b FP totals equal `EXP/E4b_fp_decomposition/summary.json`
(note: the INDEL FP numbers in the E4b README prose, e.g. 1,577 for HG008T Illumina linear, are pre-10-16 PoN; tables/summary give 1,348).

"Germline" = class 1 of E4b (FP allele equals an allele of the normal's dipcall genotype). Representation: v2.1 = FP
allele carried by ≥ 46 / 463 v2.1 haplotypes (joined to `WORK/E7_graph_composition/v2.1/hprc_alleles/chr1.tsv.gz`);
E4b's own flag = AC > 0 among the 89 v1.1 haplotypes (any carrier).

| call set | arm | FP | germline same, v2.1 ≥46 / not | germline same, v1.1 AC>0 / AC=0 | germline overlap or ≤ 10 bp | HPRC v1.1 allele (not germline) | population classes | HP≥7/STR | other |
|---|---|---|---|---|---|---|---|---|---|
| HG008T Illumina SNV | linear / d9 / d46 | 1,089 / 702 / 792 | 17/625 · 0/609 · 0/655 | 40/602 · 1/608 · 2/653 | 169 · 84 · 92 | 6 · 0 · 2 | 15 · 3 · 8 | 33 · 4 · 17 | 224 · 2 · 18 |
| HG008T PacBio SNV | linear / d9 / d46 | 969 / 583 / 650 | 28/788 · 0/571 · 0/619 | 48/768 · 1/570 · 1/618 | 49 · 11 · 19 | 16 · 0 · 0 | 0 · 0 · 4 | 5 · 0 · 4 | 83 · 1 · 4 |
| HG008T ONT SNV | linear / d9 / d46 | 904 / 660 / 788 | 23/684 · 0/618 · 0/674 | 41/666 · 2/616 · 1/673 | 89 · 23 · 52 | 12 · 1 · 10 | 6 · 1 · 9 | 7 · 11 · 19 | 83 · 6 · 24 |
| HG008T Illumina INDEL | linear / d9 / d46 | 1,348 / 390 / 385 | 73/346 · 4/170 · 0/178 | 162/257 · 50/124 · 46/132 | 738 · 191 · 183 | 47 · 6 · 4 | 43 · 10 · 8 | 78 · 8 · 9 | 23 · 1 · 3 |
| HG008T PacBio INDEL | linear / d9 / d46 | 1,523 / 2,534 / 2,324 | 117/519 · 12/260 · 3/290 | 246/390 · 89/183 · 99/194 | 720 · 1,499 · 1,544 | 69 · 227 · 177 | 72 · 289 · 185 | 13 · 85 · 48 | 13 · 162 · 77 |
| HG008T ONT INDEL | linear / d9 / d46 | 11,373 / 17,023 / 18,178 | 164/502 · 9/232 · 5/244 | 271/395 · 76/165 · 76/173 | 4,937 · 6,914 · 6,140 | 1,890 · 3,147 · 3,750 | 3,131 · 5,638 · 6,759 | 528 · 948 · 1,058 | 221 · 135 · 222 |
| COLO829T Illumina SNV | linear / d46 | 2,221 / 268 | 2/99 · 0/191 | 1/100 · 0/191 | 13 · 12 | 7 · 0 | 11 · 3 | 95 · 10 | 1,994 · 52 |
| COLO829T Fiber-seq SNV | linear / d46 | 521 / 268 | 0/133 · 0/172 | 0/133 · 0/172 | 18 · 8 | 13 · 2 | 8 · 3 | 23 · 12 | 326 · 71 |
| COLO829T ONT SNV | linear / d46 | 2,575 / 457 | 0/126 · 0/195 | 5/121 · 0/195 | 54 · 19 | 64 · 0 | 143 · 27 | 88 · 19 | 2,100 · 197 |
| COLO829T Illumina INDEL | linear / d46 | 1,527 / 384 | 89/191 · 0/145 | 130/150 · 45/100 | 610 · 171 | 31 · 4 | 42 · 20 | 88 · 17 | 476 · 27 |
| COLO829T Fiber-seq INDEL | linear / d46 | 1,601 / 3,069 | 158/431 · 4/249 | 241/348 · 96/157 | 577 · 1,465 | 88 · 421 | 77 · 781 | 40 · 62 | 230 · 87 |
| COLO829T ONT INDEL | linear / d46 | 18,380 / 27,857 | 159/489 · 11/304 | 268/380 · 127/188 | 5,415 · 8,796 | 2,788 · 3,692 | 5,963 · 8,322 | 1,150 · 1,637 | 2,416 · 5,095 |

HG008T SNV sums (three platforms): germline-same **v1.1 AC>0: 129 → 4 (d9) / 4 (d46)**; AC=0: **2,036 → 1,794 / 1,944**
(= the letter; EXISTS, E4b §5). Under stricter "in the graph" definitions: v1.1 ≥ 9/89 carriers 87 → 0 / 0
(unrepresented 2,078 → 1,798 / 1,948); v2.1 ≥ 46/463 (the d46 graph) 68 → 0 on d46 (unrepresented 2,097 → 1,948).
COLO829T SNV: represented germline FPs are 0–5 per call set on linear (v1.1 AC>0 6 in total, v2.1 2), and unrepresented
germline FPs **rise** on d46 (354 → 558); the linear excess is "other" (2,100 / 1,994 / 326). HG008T INDEL (Illumina):
represented germline 73 → 0 (v2.1), unrepresented 346 → 178, overlap/near 738 → 183; long reads: total INDEL FPs rise on
the graph (PacBio 1,523 → 2,324; ONT 11,373 → 18,178), driven by overlap/near-germline, HPRC-allele and population
classes. **MISSING:** the same split on the matched three-seed `normbq` / final INDEL models (E4b was not rerun on them).

---

## 3. "22–66", "2–8", "0–4", "57–79 %" (all EXISTS, `EXP/E9h_offref_key_reads/{README.md, tables.md §0/§3, results/key_classes.tsv}`; event dedup COMPUTED NOW from key_classes.tsv)

- **Unit.** Tensor keys (one SNP-type tensor key on an off-reference node per key) taken from
  `pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv` (d9) and `WORK/E1.1_ch6_counterfactual/asm_d46/offref_somatic_d46.tsv` (d46):
  GIAB HG008T v0.2 **truth events** (truth INFO `HG008Nv62SOMATICVARIANT`) whose somatic edit lies on an off-reference
  node of the graph (class B). COMPUTED NOW: within every run keys map 1 : 1 to events, so 22–66 are **deduplicated
  events per run**, not repeated tensors. In the normal-assembly frame most are **INDEL events** (repeat-length
  changes) seen as SNV tensors in graph coordinates: rescued tumor-only keys = INDEL / SNV events 15/7 (PacBio d9),
  15/14 (ONT d9), 46/20 (Illumina d46), 20/13 (PacBio d46), 19/13 (ONT d46). Union over the five runs: 82 distinct events
  (35 on d9 runs, 69 on d46 runs). Labels: 21–55 per run carry `residual_partial_somatic_truth` (R8: overlaps a PASS
  GRCh38 truth allele > 0.45; ignored in training), the rest −1 (`off_reference_no_truth_match` / outside region).
- **Runs.** Five runs with E1.1 scores (seed-0 7fc7396 SNV models): HG008T PacBio d9, ONT d9 (HPRC v1.1), Illumina,
  PacBio, ONT d46 (HPRC v2.1). HG008T Illumina d9 has keys (184) but no model scores. So "per HG008T dataset" = per run,
  across two graphs (v1.1: 22, 29; v2.1: 66, 33, 32).
- **Region.** All autosomes chr1–22 ("genome-wide"), including the training chromosomes chr2–22 (these tensors are not
  training examples: partial-somatic/−1).
- **Validation.** `tumor_only` under `h50 strict`: the consensus of the key's ALT-carrying reads (whole 101-column
  tensor, no indel, no mismatch within 1 bp of the site) has a hit in HG008-T v3.2 and none in HG008-N v6.2. Controls:
  T contexts tumor-only 218 / 256, N contexts normal 253 / 256.
- **Numbers.** Keys 112 / 89 / 135 / 82 / 67; tumor-only (read-supported) 72 / 51 / 93 / 65 / 50 (+ Illumina d9 115 /
  184) → **57–79 %** (57.3 % ONT d9 … 79.3 % PacBio d46; includes Illumina d9 62.5 %). Base model at its PASS threshold t
  on these keys: **0 / 0 / 3 / 1 / 4** ("0–4"). Rescue (`site10_high`, **argmax**): **22 / 29 / 66 / 33 / 32**
  ("22–66"; at threshold t: 33 / 28 / 70 / 40 / 35). Tensor level, no AF floor, so an upper bound on PASS calls.
  Precision of the rescue on all keys: 22/32, 29/46, 66/80, 33/38, 32/40 tumor-only (63–87 %).
- **chr1.** Keys 8 / 8 / 12 / 11 / 9; tumor-only 7 / 6 / 10 / 9 / 7; rescued tumor-only keys **3 / 4 / 8 / 2 / 4
  ("2–8")**, same five runs, argmax, tensor level. The comparison "against 210–1,410 graph-coordinate calls" (R3.4) is
  E9a chr1 graph-coordinate PASS of the same runs (210 / 1,410 / 754 / 240 / 1,258; per R3.3 evidence).

---

## 4. R3.7: COLO829T tumor assembly and the 35 ambiguous INDELs

- **Exists (EXISTS):** hifiasm 0.25.0 assembly by qfu of the COLO829T **Fiber-seq HiFi** reads with the COLO829T **ONT
  D003** reads as ultra-long: `/scratch/qfu/COLO829T_fiberseq/COLO829T_hifiasm/COLO829T.asm.bp.hap{1,2}.p_ctg.gfa`;
  FASTA copies `WORK/E9c_offref_filter_validation/asm_targets/COLO829T_hifiasm_hap{1,2}.fa.gz` (hap1 368 contigs,
  2,677,016,425 bp; hap2 188 contigs, 2,809,641,207 bp). In use since 2026-10-23 (E9c / E9i / E9j for R3.3). Built from
  the same long reads Pansoma analyses → not independent (`EXP/E9j_offref_colo829t/README.md` Data and Caveats).
- **Not used for R3.7's counts.** `analysis/graph_absorbed_somatic_colo829t_20261005/README.md` (2026-10-05): "No tumor
  assembly was used, so the two are not told apart: … between 7.7% (33) and 15.8% (68)". E5d was not rerun with it
  (R2.3 evidence 10-26).
- **R3.7 counts are the v1.1 read-level definition.** "Absorbed = perfect bypass on ≥ 1 platform (the union)":
  `no_candidate` on the 2026-09-30 d9 tensor sets and the majority of ALT-like reads leave GRCh38 through a non-GRCh38
  node run / skip edge with no edit within ±5 bp (read-level reasons `site_bypassed:{branch,skip_edge}:no_edit`), both
  READMEs' Definitions. HG008T 2,478 = 201 SNV + 2,277 INDEL; COLO829T 497 = 67 SNV + 430 INDEL. (E5b's tensor-level
  "no candidate and ALT spelled by the graph" is a different definition, d9 and d46.)
- **The 35 loci:** `per_variant.tsv` rows with `kind2 = INDEL`, `alt_len_haps` non-empty and `germline_status =
  germline other length` (exactly 35; a 36th ALT-length row is "normal carries ALT").
- **Tumor-assembly check: see 4b (COMPUTED NOW).**

### 4b. Tumor assembly at the 35 loci (COMPUTED NOW)

Script `MISC/r37_tumor_asm_35.py` (Slurm 462876, meow, 11 min 44 s, 2.9 GB); raw hits `MISC/r37_tumor_asm_35.tsv`, per
locus `MISC/r37_tumor_asm_35_summary.tsv`. Method: the analysis' own GRCh38 24-mer anchors (left = GRCh38[a−24, a),
right = GRCh38[b, b+24), from `anchor_a0` / `anchor_b0`; `R_seq` re-checked against GRCh38 for all 35); every exact
anchor-pair occurrence (either strand, within len(R) + 2 kb) in each hifiasm haplotype is extracted and compared with
`A_seq` (GRCh38 + truth), `R_seq`, and the COLO829BL DSA `hapX_seq` / `hapY_seq` of `per_variant.tsv`.

| outcome | loci | truth VAF median (Illumina / PacBio) | reading |
|---|---|---|---|
| a tumor haplotype carries exactly the truth allele (absent from both COLO829BL haplotypes) | **13** (1996, 4077, 15132, 18759, 24082, 28540, 35691, 38059, 39047, 40058, 41015, 41448, 43547) | 0.590 / 0.511 | somatic change on a REF-like background (reading 1) |
| tumor haplotypes equal only COLO829BL haplotype sequences (incl. the ALT-length one); no truth allele | **14** (1166, 2190, 2561, 2562, 11008, 14592, 27786, 28830, 28866, 29389, 29579, 36311, 39235, 39468; 12 in SMaHT Extreme) | 0.428 / 0.781 (several PacBio VAF ≈ 1) | germline-like (reading 2: truth = the normal's ALT-length allele written against GRCh38), or the event is missing from the assembly |
| tumor has a sequence absent from the normal but ≠ truth allele (4 with the other tumor hap = a normal hap; 4 novel only / one hap no hit) | **8** (949, 8883, 11351, 13579, 30763, 35971, 36369, 39393) | — | unresolved; 5 of them (949, 11351, 13579, 36369, 39393) have exactly the truth's length change and differ from the truth allele by 1 base (event on a SNP-carrying haplotype, or assembly error) |

Consequence for R3.7's "between 7.7% (33) and 15.8% (68)" (of 430 absorbed INDELs): with the tumor assembly the
germline-like share is **47 (10.9 %) to 55 (12.8 %)** (33 + 14, plus the 8 unresolved); 13 of the 35 resolve as somatic.
Caveats: the assembly was made from the same COLO829T long reads (consistency, not independent); a two-haplotype
assembly of an aneuploid tumor can miss subclonal or third haplotypes; 12 / 14 germline-like and 10 / 13 somatic loci
lie in SMaHT Extreme regions (repeats where assemblies err).

---

## 5. R3.m2 comparators

### 5a. ClairS-TO INDEL calling region and options (EXISTS: run logs / VCF headers; code read today from the 0.4.2 image)

| run | version | platform (model dirs) | ctg | `--bed_fn` / `--call_indels_only_in_these_regions` | INDEL min AF | SNV min AF | source |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | 0.4.2 | `ilmn_ssrs` (`ilmn_ssrs/indel/…`) | chr1 | none / none | 0.05 | 0.05 | `/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC/clairs-to_HG008T_WGS_Chr1/run_clairs_to.log` |
| HG008T PacBio | 0.4.2 | `hifi_revio_ssrs` | chr1 | none / none | 0.05 | 0.05 | `.../clairs-to_HG008_PacBio_Chr1/run_clairs_to.log` |
| HG008T ONT | 0.4.2 | `ont_r10_dorado_sup_5khz_ssrs` | chr1 | none / none | **0.1** | 0.05 | `.../clairs-to_HG008_ONT_Chr1/run_clairs_to.log` |
| COLO829T Illumina | 0.4.1 | **`ilmn`** (= `ilmn_ss`, synthetic) | chr1 | none / none | 0.05 | 0.05 | `/scratch/jshen/Pansoma_testing_results/VCF_results/clairs-to_COLO829T_WGS/run_clairs_to.log` |
| COLO829T Fiber-seq | 0.4.1 | `hifi_revio_ssrs` | chr1 | none / none | 0.05 | 0.05 | `.../clairs-to_COLO829T_Pacbio/run_clairs_to.log` |
| COLO829T ONT | 0.4.1 | `ont_r10_dorado_sup_5khz_ssrs` | chr1 | none / none | **0.1** | 0.05 | `.../clairs-to_COLO829T_ONT/run_clairs_to.log` |
| HCC1395 ×3 | 0.4.2 | `ilmn_ssrs` / `hifi_revio_ssrs` / `ont_r10_dorado_sup_5khz_ssrs` | chr1,chr20 | none / none | 0.05 / 0.05 / 0.1 | 0.05 | `WORK/E8_comparators/calls/clairs_to/HCC1395T_*/run_clairs_to.log` |

- All runs: `DISABLE INDEL CALLING: False`, `BED FILE PATH: None`; the extract step is given
  `--call_indels_only_in_these_regions <out>/tmp/split_indel_beds/chr1`, but that directory is **empty** in every run
  (no BED supplied). In `extract_candidates_calling.py` (0.4.2 image, line ~399) `pass_indel_bed_region = len(indel_tree)
  == 0 or …`, so with no BED every INDEL candidate on the contig passes. ClairS-TO changelog (image README): v0.1.0
  restricted INDELs to "whole genome − GIAB v3.3 all difficult + CMRG"; **v0.3.0 removed this restriction**. So for
  0.4.1 and 0.4.2, INDELs were called on all of chr1 (chr1 + chr20 for HCC1395) on every platform; all three model
  families are in the INDEL-supported platform list (`run_clairs_to` line ~539). Saved copies: `MISC/run_clairs_to_v0.4.2.txt`,
  `MISC/ec_v0.4.2.txt`, `MISC/clairs_to_readme_v0.4.2.txt`. (0.4.1 code not inspected; its run logs show the same empty
  INDEL BED.)
- ClairS-TO PoN tags (VCF header, all nine runs): PoN_1 gnomAD r2.1 AF ≥ 0.001 (allele matching **True**), PoN_2 dbSNP
  138 non-somatic (**True**), PoN_3 1000G PoN (**False** = position only), PoN_4 CoLoRSdb v1.1.0 AF ≥ 0.001 (**False**).
  Tagged records get FILTER `NonSomatic`.
- DeepSomatic: `--regions=chr1` (HCC1395: "chr1 chr20"), no BED; `--use_default_pon_filtering=true` (new runs explicit;
  original RIS runs inferred from FILTER=PON records and PON VCFs in the postprocess logs); `calls/deepsomatic/*/command.txt`.

### 5b. Matching and filter handling in scoring (EXISTS: `EXP/E8_comparators/score_comparator.py`, `machine_learning/pansoma_net_v2/vcfeval.py`)

- Calls: `bcftools norm -m -any -f <labels FASTA>` (split multiallelics, left-align), chr1, records of the kind only,
  RefCall / NoCall dropped, one record per (POS, REF, ALT) keeping the highest QUAL; sites-only VCF with tool QUAL/FILTER.
- rtg vcfeval 3.12.1: `--squash-ploidy --sample ALT,ALT -f QUAL --bed-regions chr1` (+ `-e <benchmark BED>` for
  eval_bed), `--all-records` only for the curves; **no `--decompose`, no `--ref-overlap`**. So matching is rtg's
  haplotype-aware allele matching with genotype ignored, the same function Pansoma's `vcf.sh` uses.
- Truth = all PASS truth records of the kind on chr1 (not PoN-filtered); the PoN-tagged truth is only used for the
  "maximum attainable recall" bound.
- Sets: `pass` = FILTER exactly PASS ("native"); `pass_pon` = `pass` minus records our PoN tags ("harmonized");
  `all` / `all_pon` = every non-RefCall/NoCall record (any FILTER) without / minus our PoN, ranked by QUAL (curves).
  Consequences: **native population tags (ClairS-TO `NonSomatic`, DeepSomatic `PON` and `GERMLINE`) and quality filters
  are treated as filtered in both the native and the harmonized operating point** (harmonized = native PASS ∩ not our PoN,
  i.e. the union of both filters, not an identical total filter); **in the PR curves (`all_pon`) they are not applied**,
  only our PoN. E8 README §7: for DeepSomatic, GERMLINE records carry high QUAL, so only `all_pon` curves are meaningful.
- Our PoN (`scripts/filter_panel_of_normals.py`, the four resources of `jobs/vcf.sh`): SNV gnomAD / CoLoRSdb AF ≥ 1e-4
  exact allele + dbSNP non-somatic + 1000G; INDEL gnomAD / CoLoRSdb exact allele, AF ≥ 0.01 in HP ≥ 8, else ≥ 1e-4
  (71cb51a; HG008T Illumina INDEL reports re-scored with it, array 418385).

### 5c. Corrected HG008T ONT DeepSomatic run is what every table reads (EXISTS, checked today)

- `WORK/E8_comparators/score/deepsomatic/HG008T_ONT/{SNV,INDEL}/report.json` → `input_vcfs` =
  `WORK/E8_comparators/calls/deepsomatic/HG008T_ONT/output.vcf.gz`, whose `command.txt` is `--model_type=ONT_TUMOR_ONLY`
  (jobs 417841 + 418308). There is only this one HG008T_ONT directory under `score/deepsomatic/`.
- Consumers (all read `WORK/E8_comparators/score/<caller>/<dataset>/<KIND>/report.json` or its rtg dirs):
  `EXP/E8_comparators/latest_models.py` → `latest_models.md/json` (rows HG008T_ONT DeepSomatic-TO SNV 0.310 / 0.869 /
  0.457 = pass_pon of that report; INDEL 0.172 / 0.019 / 0.034); `EXP/E8_comparators/collect.py` → `tables.md`,
  `summary.json`, `sbs96.tsv` (→ `RR/supplementary/R3.5/make_sbs96.py`); `EXP/E16_bootstrap/comparators.tsv` (E16 / R3.9
  paired tables, `RR/supplementary/R3.9/make_paired_table.py`); `EXP/E12b_multicancer_models/collect.py`;
  `EXP/E8c_comparators_absorbed/absorbed_comparators.py`; `EXP/E8_comparators/purity_comparators.py` (COLO829T only).
- The wrong-model VCF (`/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC/deepsomatic_HG008T_ONT_Chr1/`,
  PACBIO_TUMOR_ONLY) is referenced by **no** script in `analysis/` (grep). Old numbers (submission Supplementary: 618 TP /
  7,549 FP, P 0.0757) appear nowhere in the current tables.
- **MISSING:** a separate PR-curve figure file; the curves exist only as rtg `--all-records` outputs (`all_pon`) inside
  the same score dirs, so they inherit the corrected run.

### 5d. Training exposure of the comparators (EXISTS, `EXP/E8b_literature_facts/README.md` rows D/C)

| tumor | DeepSomatic 1.9.0 tumor-only | ClairS-TO `ssrs` (0.4.x) |
|---|---|---|
| HG008T | not trained on in the paper; held-out whole-genome evaluation line (D6). Whether 1.9.0 WGS_TUMOR_ONLY "retrained with all datasets described in the manuscript" (D15) includes HG008: unverified | not in the paper (0 occurrences, C5); not trained on |
| COLO829T | not named in the article (D7); CASTLE ships COLO829 data "used in the manuscripts" (ONT R9, Illumina): training use unverified | not trained on (C3); a paper benchmark (NYGC truth, C5). COLO829T Illumina run used `ilmn` = SS model (synthetic HG003/HG004, C2) |
| HCC1395 | paper's tumor-only models trained on HCC1395 (SEQC2 truth), chr1 held out (D10, D2, D14) | excluded entirely (C3) |
| HCC1954, H1437, H2009 | in the multi-cancer training lines (D4) and CASTLE truth was made by DeepSomatic (D8); paper's tumor-only models = HCC1395 only (D10); 1.9.0 WGS_TUMOR_ONLY retrained on "all datasets" (D15), long-read TO checkpoints' data unstated (D16); chr1 held out in every model (D2) | SSRS fine-tuned on HCC1937, HCC1954, H1437, H2009 with DeepSomatic multi-cancer truth (C3); held-out chromosome chr20 (README) or chr1 (paper) — unresolved (C4) |

---

## 6. Statements in the letter that conflict with these facts

| section | letter text | correction | basis |
|---|---|---|---|
| R2.2 | "27.8% are common population alleles, 97.5% of them where the B-lymphoblastoid normal lost heterozygosity (chr9q, chr1:9–12 Mb)" | 97.5 % = 22,397 / 22,975 **autosomal** common SNVs; **all of them in chr9:68–139 Mb**; chr1:9–12 Mb has 30 truth SNVs and none common → drop "chr1:9–12 Mb" (or say "chr9q") | 1b |
| R2.2 | "On chr1 … negligible (19 of 4,266 truth SNVs)" | correct as a count (all chr1 records; 17 / 3,941 in BED), but it is a different set from the 97.5 %: none of the 19 lies in a normal-LOH bin; say "19 of 4,266 chr1 truth SNVs are common alleles" | 1b |
| R2.2 | "about 30% of H1437 truth SNVs are probably germline" | keep, but it is the segment count 24,914 / 82,529 (30.2 %), equal to the mixture estimate, not a per-variant call | 1b |
| R2.2 | "(H2009 about 5%, HCC1954 none)" | HCC1954 has no normal-LOH excess but a 2.1 % mixture estimate (~500; IG/TCR, chr7:90–91 Mb); "none" only for the LOH mechanism | E22 T4 |
| R1.m1 (feeds R2.2 reconciliation) | "population filtering tags all but two of them" | H1437: all 18,028 of 18,028; "all but two" is H2009 (6,635 / 6,637) | E22 README check 3 |
| R1.m1 | "99% and 91% of these alleles lie in segments where … lost heterozygosity" | they lie in 1-Mb **bins** with normal het < 0.2 (not the ≥ 3 Mb segments of R2.2) | graph_overlap.tsv |
| R3.2 | "carrying a germline allele of the normal that is in HPRC v1.1 fall from 129 … to 4 … whereas those absent from the graph fall only from 2,036 to 1,794–1,944" | "in HPRC v1.1" = carried by ≥ 1 of 89 v1.1 haplotypes, and "absent" = 0 carriers, not graph membership; with the graph's own criterion: ≥ 9/89 87 → 0 / 0, ≥ 46/463 (d46) 68 → 0; unrepresented 2,078 → 1,798 / 1,948. Seed-0 7fc7396 runs at val-recall-0.9; PacBio/ONT d9 = 1-GPU recipe | 2c |
| R3.2 | SNV only ("Stratified by graph content …") | INDEL and COLO829T now available (2a, 2c): represented somatic INDELs fall to 0.6–34 %; unrepresented somatic INDELs also lose 12–28 points; long-read germline INDELs represented in the graph keep 11.6–15.9 % candidates; long-read INDEL FPs increase on the graph | 2a, 2c |
| R3.2 | "Non-truth SNV candidates fall by 5–61% (non-truth INDEL … 75–84% on Illumina …)" | consistent (2b) | 2b |
| R3.4 | "read-supported off-reference somatic SNVs that rtg cannot credit (22–66 per HG008T dataset genome-wide, 0–4 without it; R3.3)" | per **run** (5 runs: PacBio and ONT on v1.1; Illumina, PacBio, ONT on v2.1; v1.1 22–29, v2.1 32–66), all autosomes; distinct GIAB truth events, mostly INDEL events in the normal frame (15–46 per run) seen as SNV tensors; argmax at tensor level without AF floor (upper bound on PASS); base 0–4 at t | 3 |
| R3.4 | "on chr1 these are 2–8 per HG008T run" | correct (3 / 4 / 8 / 2 / 4, same five runs, tensor level) | 3 |
| R3.5 | "recovers 22–66 assembly-validated somatic SNVs per HG008T dataset that GRCh38 cannot represent" | these are events of the GRCh38 truth (GIAB v0.2 records) whose graph representation is off-reference, so "GRCh38 cannot represent" is inaccurate → "whose graph-coordinate calls cannot be scored on GRCh38"; validation is read-consensus vs assemblies (h50 strict), per run as above | 3 |
| R1.4 | "57–79% have read support, 0–4 are called" | correct; 57–79 % includes HG008T Illumina d9 (no model scores), 0–4 over five runs | 3 |
| R3.7 | "Without a tumor assembly we cannot resolve 35 further INDELs … between 7.7% and 15.8%" | a COLO829T tumor assembly exists (qfu hifiasm 0.25.0, Fiber-seq HiFi + ONT D003; used in R3.3) but was not used in this analysis and is not independent. Checked today: at 13 of the 35 it carries exactly the truth allele, at 14 only normal haplotype sequences, 8 unresolved → germline-like share 10.9–12.8 % (47–55 of 430), with the "same reads" caveat | 4, 4b |
| R3.7 | "haplotype-resolved assemblies of the matched normals … and, for HG008T, of the tumor" | consistent with the analysis as run (COLO829T tumor assembly not used); add the v1.1 read-level definition ("perfect bypass on ≥ 1 platform") | 4 |
| R3.m2 | "ClairS-TO used its default tagging (gnomAD and CoLoRSdb AF ≥ 0.001, dbSNP 138 non-somatic sites, the GATK 1000 Genomes panel of normals)" | add matching: gnomAD and dbSNP exact allele, 1000G PoN and CoLoRSdb position only | 5a |
| R3.m2 | "We report native PASS calls and PASS calls after applying the same additional population filter" | state that both operating points keep the tools' own FILTERs (NonSomatic, PON, GERMLINE, LowQual …), so the harmonized set = native ∪ our PoN, not an identical total filter; PR curves use every non-RefCall record minus our PoN, i.e. native population tags are not applied there | 5b |
| R3.m2 | "same truth, benchmark regions and rtg vcfeval options" | specify: bcftools norm -m -any (left-align), rtg vcfeval `--squash-ploidy --sample ALT,ALT -f QUAL`, no decomposition, full (unfiltered) truth | 5b |
| R3.m2 | (INDEL calling region not stated) | ClairS-TO called INDELs on all of chr1 (no BED; whole-genome default since v0.3.0); INDEL min AF 0.1 on ONT, 0.05 otherwise | 5a |
| R3.m2 | "ClairS-TO's ssrs models exclude HCC1395" | correct; add that ssrs were fine-tuned on HCC1954 / H1437 / H2009 (CASTLE, DeepSomatic truth) and DeepSomatic's training lines include them, so CASTLE results are not held-out for either comparator (chr1 is held out by DeepSomatic; ClairS-TO hold-out chr20 or chr1 unresolved) | 5d |
| R3.3 | "COLO829T: a hifiasm assembly of the tumor's own HiFi and ONT reads" | consistent (Fiber-seq HiFi + ONT D003, qfu) | 4 |
| R2.3 | "COLO829T has no tumor assembly independent of the sequencing data analysed here" | consistent | 4 |
