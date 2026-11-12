# Figure data: truth VAF (Illumina, PacBio) and germline overlap with GRCh38-only, HPRC v1.1, HPRC v2.1

Data for the three R scripts the user supplied on 2026-11-11 (`COLO829T_VAF_violin.R`, `HG008T_VAF_violin.R`,
`overlapping_HPRC_germline_HG008.R`), now for all three tumors (HG008T, COLO829T, HCC1395T) and all three references.

## 1. Truth VAF (`truth_vaf/`)

Every somatic truth allele of the truth VCF that the labels use (chr1-22; SNV / INS / DEL, multi-ALT split, other
alleles skipped), with its VAF in Illumina and in PacBio reads:

| tumor | truth VCF | Illumina | PacBio |
|---|---|---|---|
| HG008T | GIAB HG008-T somatic smvar benchmark v0.2 (tumorvariants), BED v0.2_all | read counts, `HG008-T_Illumina_195x_GRCh38-GIABv3.bam` | read counts, `HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam` |
| COLO829T | SMaHT COLO829T reference set (`COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz`), BED SMaHT v2 easy+difficult+extreme | SMaHT `VAF_Ill` (truth VCF INFO) | SMaHT `VAF_PB` (truth VCF INFO) |
| HCC1395T | SEQC2 v1.2.1 (HighConf + MedConf), BED SEQC2 High-Confidence Regions v1.2 | read counts, SEQC2 SRR7890824 tumor WGS (bwa-mem GRCh38, `SEQC2/SRR7890824_Tumor/SRR7890824_Tumor.sorted.bam`) | read counts, CASTLE HCC1395 HiFi (`CASTLE/HCC1395/PacBio/HCC1395.GRCh38.bam`) |

- `count_vaf.py`: the counting rules of the earlier HG008T truth-VAF table (`/scratch/jshen/tmp/HG008_trueset_germline_overlapping/count_vaf_from_vcf_with_pysam.py`)
  unchanged: pysam pileup at POS, primary alignments, MAPQ >= 20, base quality >= 20; SNV = REF / ALT base;
  INS / DEL = indel of the truth length right after the anchor base; VAF = ALT reads / reads with a base at POS.
  Checked: the first chr1 records reproduce that table exactly (DP, REF and ALT counts).
- COLO829T keeps the SMaHT values, the source of the earlier COLO829T violin (no whole-genome COLO829T Illumina tumor BAM
  on this cluster; E8 has chr1 only).
- HCC1395T Illumina: Pansoma's HCC1395T Illumina sets use the SEQC2 WGS_FD_T_1 + 2 reads, of which only chr1 / chr20 are
  aligned to GRCh38 (E8). SRR7890824 is a whole-genome SEQC2 HCC1395 tumor library; `summary.md` compares the two on chr1.
- `make_tables.py` writes:
  - `<T>_truth_VAF_Illumina_PacBio.csv` (CHROM, POS, REF, ALT, TYPE, PASSED, IN_BED, [RGN,] DP_Ill, ALT_COUNT_Ill,
    VAF_Ill, DP_PB, ALT_COUNT_PB, VAF_PB): input of `COLO829T_VAF_violin.R` (reads `VAF_Ill`, `VAF_PB`; the earlier
    COLO829T input was SNV only, so filter `TYPE == "SNV"` for the same plot);
  - `counts/<T>_<Illumina|PacBio>.truth_vaf.tsv` (TYPE, VAF, ...): input of `HG008T_VAF_violin.R` (reads `TYPE`, `VAF`);
  - `summary.md`: n, quartiles and low-VAF shares.
- IN_BED: the REF span (anchor base included) inside one merged interval of the truth BED.

## 2. Germline overlap (`germline_overlap/`)

`overlap.py`: germline alleles of the three normals (HG008N, COLO829BL, HCC1395BL) against the alleles of each graph,
per chromosome, SNV and INDEL (<= 50 bp), counted as distinct normalized allele keys. Inputs are the E7 allele tables
(`analysis/reviewer_response_20261011/experiments/E7_graph_composition`; WORK under
`/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition`), where the normal's germline truth and every
HPRC allele are left-normalized against GRCh38 and decomposed the same way:

| graph | alleles |
|---|---|
| `GRCh38` (GRCh38-only) | none: the linear graph is vg autoindex of the GRCh38 FASTA alone |
| `HPRC_v1.1_d9` | HPRC v1.1 alleles carried by >= 9 of the 89 haplotypes (the d9 filter) |
| `HPRC_v2.1_d46` | HPRC v2.1 alleles carried by >= 46 of the 463 haplotypes (the d46 filter) |
| `HPRC_v1.1_d9_VCF` (check) | the alleles of the real d9 graph's own deconstruct VCF (complex ALTs not decomposed) |

There is no deconstruct VCF of the d46 graph (HPRC publishes the d46 GBZ only), so both HPRC graphs use the same
carrier-count rule; the d9 VCF row shows how close the rule is to the real d9 graph.

Normal germline truth (the labels' germline VCF / BED): HG008-N dipcall (HG008N_GRCh38_dipcall.dip.vcf.gz / .dip.bed),
COLO829BL dipcall (`/scratch/qfu/COLO829BL_DSA/dipcall_hg38`), HCC1395BL DeepVariant hybrid QUAL >= 10 with its
high-confidence BED.

Scopes: `all` = every chr1-22 allele of both sides (as the earlier HG008N figure, which intersected whole VCFs);
`bed` = normal alleles with FILTER PASS inside the normal's germline BED, and graph alleles whose REF span lies in that BED.

Outputs: `overlap_long.tsv` (all numbers); `R_input/isec_<graph>_vs_<normal>.<scope>.per_chr.snv_indel.tsv` in the
layout of the earlier HG008N table (A = graph, B = normal), read by `overlapping_HPRC_germline_HG008.R` as
AF-HPRC = A and HG008N = B; `excluded.tsv` (alleles > 50 bp, COMPLEX, multi-primitive normal alleles).
`summarize.py` writes `summary.md` (chr1-22 and chr1 totals per scope).

Results, chr1-22, share of the normal's germline alleles that the graph carries (`summary.md`; GRCh38-only is 0 for
every normal by construction):

| normal | graph | SNV, all | INDEL, all | SNV, bed | INDEL, bed |
|---|---|---|---|---|---|
| HG008N | HPRC v1.1 d9 | 82.4% | 74.8% | 86.9% | 76.1% |
| HG008N | HPRC v2.1 d46 | 82.7% | 76.5% | 87.3% | 77.8% |
| COLO829BL | HPRC v1.1 d9 | 82.7% | 74.3% | 87.2% | 75.6% |
| COLO829BL | HPRC v2.1 d46 | 83.1% | 75.8% | 87.7% | 77.2% |
| HCC1395BL | HPRC v1.1 d9 | 84.2% | 75.8% | 88.1% | 79.8% |
| HCC1395BL | HPRC v2.1 d46 | 84.4% | 76.9% | 88.6% | 80.9% |

- Checks: the bed-scope numbers are those of E7 Table I (germline truth carried by each graph); the v1.1 graph sizes
  (all scope: 6,869,974 SNV, 2,123,519 INDEL) are those of E7 Table E; the real d9 VCF gives the same shares within
  1 point (HG008N bed 87.0% / 75.1%).
- v2.1 d46 carries fewer alleles than v1.1 d9 (6.30 M vs 6.87 M SNV, 1.80 M vs 2.12 M INDEL, all scope) but more of each
  normal's germline alleles (+0.2 to +0.5 points SNV, +1.1 to +1.7 points INDEL). Both filters keep alleles on about 10 %
  of the haplotypes (9 / 89, 46 / 463).

## Jobs

`tmp/vaf_germline_overlap_20261111/` (`count_vaf.sh`, `overlap.sh`, `jobs.txt`, `logs/`).
