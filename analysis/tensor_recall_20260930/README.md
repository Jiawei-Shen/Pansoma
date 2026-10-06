# Somatic truth VCF → tensors, and the four PoNs: HG008T and COLO829T, three platforms each

2026-09-30. For every tensor set: how many somatic truth VCF records are built into tensors, how many are missing and
why, SNV and INDEL separately. A second part compares the truth and the tensors with the four PoN VCFs.
All numbers come from the files listed at the end; `tables.md` and `pon/pon_tables.md` have every table in full.

## Data

| Sample | Platform | Tensor set (`/scratch/jshen/data/pansoma_v2_tensors/`) | Truth |
|---|---|---|---|
| HG008T | PacBio HiFi Revio 116× | `HG008T_PacBio/tensors` | GIAB `HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz` (18,952 records, all PASS), BED `_all.bed` |
| HG008T | ONT-UL R10.4.1 54× | `HG008T_ONT/tensors` | same |
| HG008T | Illumina WGS (BCM) | `HG008T_Illumina/tensors` | same |
| COLO829T | PacBio fiberseq | `COLO829T_fiberseq/tensors` | `COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz` (46,064 records, all PASS; INFO VAF_Ill, VAF_PB, RGN); BED = whole genome |
| COLO829T | ONT | `COLO829T_ONT/tensors` | same |
| COLO829T | Illumina WGS | `COLO829T_Illumina/tensors` | same |

- The tensors cover chr1–22 only. Build thresholds on every set: SNV AF 0.06, INDEL AF 0.08, at least 3 supporting reads,
  MAPQ > 10.
- Labels:
  - Rules `rules_sha256 197b5d25bbf4`.
  - SNV label AF floor: HG008T 0.08 (relabelled 2026-09-30); COLO829T Illumina 0.07; none on the others.
  - The HG008T `README.txt` files in the data directory still say `snv_min_af None`. The `labels.manifest.json` files
    say 0.08.
- A truth allele's status comes from `<set>/tensors/somatic.recall.tsv` (key match; partial matches do not count).
  - **built**: a tensor contains the allele, as its representative allele A1 or as another allele of the site.
  - **filtered**: the builder made the candidate but dropped it. `AF` = below the AF floor; `<3 reads` = fewer than 3
    supporting reads.
  - **no candidate**: no read produced this allele as an edit on a target node.

## Recent updates and test results this builds on

- **Tensors:** unchanged since 2026-09-27.
- **Labels:**
  - Haplotype-overlap fix, relabelled 2026-09-28 (69bdb2a).
  - HG008T SNV relabelled with the 0.08 AF floor (label manifests dated 2026-09-30 UTC).
  - The recall statuses do not depend on the labels.
- **PoN rule** (`scripts/filter_panel_of_normals.py`):
  - Allele match for all four PoNs (a3365c7).
  - gnomAD and CoLoRSdb only at ALT AF ≥ 1e-4 (79afb23).
  - dbSNP only if not somatic (SAO ≠ 2).
  - 1000G: any matching allele.
- **Latest chr1 tests:** pansoma_net_v2 `HG008_*_SNV_allnon_sqrt_b1024` with the PoN at 1e-4, `eval_bed`, 697 chr1 truth
  records. The rtg "ceiling" is the recall of all predictions at the lowest score.

  | chr1 in BED, HG008T SNV | P | R | F1 | ceiling raw / after PoN | tensor recall here: built / built and not PoN-tagged |
  |---|---:|---:|---:|---:|---:|
  | Illumina | 0.474 | 0.825 | 0.602 | 0.931 / 0.884 | 94.7% / 89.4% |
  | PacBio | 0.485 | 0.826 | 0.612 | 0.928 / 0.885 | 94.5% / 89.4% |
  | ONT | 0.462 | 0.839 | 0.596 | 0.937 / 0.891 | 95.3% / 90.0% |

  HG008T Illumina INDEL (`nopartial_b1024_nopc`, no PoN, 587 chr1 records in BED): R 0.365, ceiling 0.400. On chr1 in
  the BED, 40.6% of Illumina INDEL truth alleles are built.

  The test excludes tensors labelled −1 below the AF floor, so the test ceiling sits a little under the tensor recall.
- **Still running:** `HG008_{Illumina,PacBio,ONT}_SNV_allnon_sqrt_nopc_lr1e4`, with their VCF jobs queued.

## 1. Truth VCF records built into tensors

chr1–22 truth alleles. Multi-allelic records are split: HG008T has 15 INDEL records = 30 alleles. chrX/Y records get no
tensors, so they are all missing. Percentages are of the chr1–22 alleles.

| Sample | Platform | Kind | VCF records (on chrX/Y) | chr1–22 alleles | **built** | filtered: AF / <3 reads / both | no candidate |
|---|---|---|---:|---:|---:|---:|---:|
| HG008T | PacBio HiFi | SNV | 9,609 (919) | 8,690 | **8,236 (94.8%)** | 18 / 38 / 0 | 398 (4.6%) |
| HG008T | ONT | SNV | 9,609 (919) | 8,690 | **8,280 (95.3%)** | 2 / 66 / 0 | 342 (3.9%) |
| HG008T | Illumina | SNV | 9,609 (919) | 8,690 | **8,283 (95.3%)** | 13 / 55 / 0 | 339 (3.9%) |
| HG008T | PacBio HiFi | INDEL | 9,343 (862) | 8,496 | **4,445 (52.3%)** | 463 / 500 / 112 | 2,976 (35.0%) |
| HG008T | ONT | INDEL | 9,343 (862) | 8,496 | **4,490 (52.8%)** | 245 / 1,228 / 79 | 2,454 (28.9%) |
| HG008T | Illumina | INDEL | 9,343 (862) | 8,496 | **3,633 (42.8%)** | 239 / 637 / 33 | 3,954 (46.5%) |
| COLO829T | PacBio fiberseq | SNV | 44,005 (1,970) | 42,035 | **38,550 (91.7%)** | 2,540 / 417 / 0 | 528 (1.3%) |
| COLO829T | ONT | SNV | 44,005 (1,970) | 42,035 | **38,357 (91.3%)** | 1,671 / 1,460 / 0 | 547 (1.3%) |
| COLO829T | Illumina | SNV | 44,005 (1,970) | 42,035 | **38,629 (91.9%)** | 2,106 / 702 / 0 | 598 (1.4%) |
| COLO829T | PacBio fiberseq | INDEL | 2,059 (147) | 1,912 | **1,249 (65.3%)** | 131 / 67 / 29 | 436 (22.8%) |
| COLO829T | ONT | INDEL | 2,059 (147) | 1,912 | **1,303 (68.1%)** | 124 / 109 / 29 | 347 (18.1%) |
| COLO829T | Illumina | INDEL | 2,059 (147) | 1,912 | **1,104 (57.7%)** | 19 / 99 / 3 | 687 (35.9%) |

`tables.md` has more tables:

| Table | Content |
|---|---|
| §1 | representative vs other allele, and % of all VCF records |
| §2 | INS vs DEL |
| §6 | labels |
| §7 | chr1 |
| §8 | the same truth allele across platforms |

Headlines:
- **INS vs DEL:** INS is worse than DEL everywhere. HG008T built INS 36–44% vs DEL 52–65%; COLO829T INS 50–56% vs DEL
  64–79%.
- **Labels:** every truth allele built as a representative allele is labelled 1. The exceptions are the SNV AF floor
  (label −1 `below_snv_min_af`): HG008T 3 / 7 / 6 (PacBio / ONT / Illumina) and COLO829T Illumina 389.
- **chr1** (the test chromosome):
  - HG008T in the BED: SNV 94.5 / 95.3 / 94.7% and INDEL 49.8 / 51.7 / 40.6% built (PacBio / ONT / Illumina).
  - COLO829T: SNV 84.4 / 84.4 / 85.0% and INDEL 51.1 / 52.7 / 43.5% built (fiberseq / ONT / Illumina). COLO829T chr1 has
    relatively more low-VAF SNV truth: 13–14% filtered, vs 7% genome-wide.
- **Across platforms:**
  - HG008T: SNV built on all three 94.0%, on ≥ 1 96.3%; INS on ≥ 1 50.0%; DEL on ≥ 1 73.8%.
  - COLO829T: SNV on all three 89.0%, on ≥ 1 94.2%; INS on ≥ 1 60.1%; DEL on ≥ 1 82.6%.

## 2. Why truth alleles are missing

### No candidate: read-level classes (tables.md §3)

Every no-candidate truth allele was examined at read level with the method of `analysis/hg008_somatic_miss_20260926`:
- the reads over the site are fetched from the GAM and decoded with the run's frozen decoder;
- the ALT reads are compared with the REF and ALT haplotypes;
- the allele is classed by how the ALT reads represent it.

HG008T uses the 2026-09-26 per-truth classes, which still cover every current miss. COLO829T was analysed here, in
`colo829t_miss/`. Controls (200 built INDELs per set) check the method: on COLO829T 187/200 (Illumina), 188/200
(fiberseq) and 179/200 (ONT) show the truth edit on a built node.

| Class | Meaning |
|---|---|
| S1 | no MAPQ>10 read spans the site |
| S2 | < 3 ALT reads |
| S3 | ALT is an existing graph SNP allele (reads take that branch, so there is no edit) |
| S4a/b/c | repeat: the SNV becomes an edit on another branch / is absorbed by a graph path / no clear ALT |
| I1 | INDEL allele fully in the graph (reads take the branch or skip edge, so there is no edit) |
| I2 | INDEL allele partly in the graph (only a residual edit with another spelling) |
| I3 | reads stay on GRCh38 but spell another allele |
| I4 | < 3 ALT reads |
| I5 | no MAPQ>10 reads |
| I6 | read indel > 50 bp |

| No candidate | HG008T PacBio | HG008T ONT | HG008T Illumina | COLO829T fiberseq | COLO829T ONT | COLO829T Illumina |
|---|---:|---:|---:|---:|---:|---:|
| SNV total | 398 | 342 | 339 | 528 | 547 | 598 |
| S1 low MAPQ / no reads | 48 | 1 | 55 | 227 | 45 | 237 |
| S2 < 3 ALT reads | 46 | 36 | 40 | 109 | 356 | 200 |
| S3 ALT = graph SNP allele | 65 | 46 | 76 | 58 | 52 | 65 |
| S4a+b+c repeats | 239 | 259 | 168 | 134 | 94 | 96 |
| INDEL total (INS + DEL) | 2,976 | 2,454 | 3,954 | 436 | 347 | 687 |
| I1 + I2 allele fully / partly in the graph | 2,739 (92.0%) | 2,237 (91.2%) | 3,420 (86.5%) | 354 (81.2%) | 280 (80.7%) | 565 (82.2%) |
| I3 GRCh38 path, other spelling | 100 | 113 | 106 | 43 | 43 | 59 |
| I4 < 3 ALT reads | 67 | 60 | 306 | 14 | 13 | 33 |
| I5 no MAPQ>10 reads | 29 | 1 | 101 | 19 | 4 | 26 |
| I6 > 50 bp, I7 unclear, I9 | 41 | 43 | 21 | 6 | 7 | 4 |

**SNV:**
- HG008T misses are mostly the graph's own alleles and repeats: S3 + S4 are 72–89% of the misses.
- COLO829T has two reasons HG008T does not:
  - S1 sits in SMaHT "Extreme" regions (fiberseq 202/227, Illumina 208/237) at normal truth VAF (median 0.34–0.37).
  - S2 is low-VAF truth: Illumina median truth VAF 0.026, 168/200 below 0.10. ONT has the most S2 (356 of 547 misses,
    median truth VAF 0.031, 291 below 0.10, 236 in Easy regions).
- ONT reads map in the segdups: S1 is 1 on HG008T ONT and 45 on COLO829T ONT, vs 48–55 (HG008T) and 227–237
  (COLO829T) on the other platforms.

**INDEL:** in all six sets the allele is already fully or partly in the HPRC graph (I1 + I2 = 81–92%). The ALT reads
then follow a graph branch or skip edge, and there is no edit to make a candidate from. COLO829T INDEL misses sit mostly
in Difficult/Extreme regions.

### Filtered (tables.md §4, §5)

These candidates were built but dropped:

| Set | Kind | Reason | Median ALT reads | Median AF |
|---|---|---|---:|---:|
| all | any | `<3 reads` | 1–2 | – |
| HG008T PacBio | INDEL | AF | 4 | 0.054 |
| COLO829T Illumina | SNV | AF | 5 | 0.040 |
| COLO829T fiberseq | SNV | AF | 6 | 0.033 |

- **COLO829T SNV:** the AF floor removes low-VAF, subclonal truth. 3,062 truth SNVs have VAF_PB < 0.06 (3,037 by VAF_Ill),
  and 74–82% of those are filtered. From VAF ≥ 0.10 upward, 89–99.5% are built.
- **HG008T ONT INDEL:** `<3 reads` dominates (1,228). ONT spreads one allele over several edit spellings.

## 3. The four PoNs (pon/pon_tables.md)

- **Rule:** the repo rule (`scripts/filter_panel_of_normals.py`), applied by streaming each PoN over each chromosome
  once.
- **Cross-check:** exact agreement with the script on:
  - the chr22 truth (713 alleles);
  - the PoN outputs of the three chr1 SNV vcfeval runs (9,999 / 11,368 / 15,515 records);
  - one chr1 INDEL linear VCF (7,544 records), run through the script here.
- **Tensor alleles:** built as `pansoma_net_v2.linear_vcf` builds them (GRCh38 projection, left-aligned, padded).
  Off-reference tensors have no GRCh38 allele, so no PoN reaches them.

**Truth alleles in the PoNs** (chr1–22, tagged by the rule):

| Sample | Kind | truth alleles | allele in any PoN (any AF) | **tagged** | main PoN |
|---|---|---:|---:|---:|---|
| HG008T | SNV | 8,690 | 1,383 (15.9%) | **910 (10.5%)** | gnomAD 618, CoLoRSdb 605, dbSNP 501 |
| HG008T | INDEL | 8,496 | 7,160 (84.3%) | **7,078 (83.3%)** | CoLoRSdb 6,869 (4,225 only there), dbSNP 2,088 |
| COLO829T | SNV | 42,035 | 3,629 (8.6%) | **1,442 (3.4%)** | gnomAD 1,098, dbSNP 737, CoLoRSdb 537 |
| COLO829T | INDEL | 1,912 | 1,383 (72.3%) | **1,356 (70.9%)** | CoLoRSdb 1,284, dbSNP 464 |

- **Truth SNVs in gnomAD:** most match at very low population AF (HG008T median 1.5e-4, COLO829T 4.9e-5). The 1e-4
  floor spares 516 HG008T and 2,332 COLO829T gnomAD matches.
- **Truth INDELs:** they are common population alleles, CoLoRSdb median AF 0.022 (homopolymer/STR alleles; see also
  `analysis/indel_20260930`). An INDEL PoN would remove 71–83% of the true somatic INDELs, which is why INDELs run
  without a PoN.
- **chr1 cross-check:** the chr1 HG008T numbers match vcfeval's truth tags (SNV 61/702; INDEL 544 of 636 alleles =
  635 records).

**Recall ceiling before / after the PoN** (built, and built and not tagged; chr1–22):

| Set | SNV built → after PoN | INDEL built → after PoN |
|---|---:|---:|
| HG008T PacBio | 94.8% → 88.4% | 52.3% → 12.9% |
| HG008T ONT | 95.3% → 88.9% | 52.8% → 12.8% |
| HG008T Illumina | 95.3% → 88.5% | 42.8% → 12.6% |
| COLO829T fiberseq | 91.7% → 88.7% | 65.3% → 23.6% |
| COLO829T ONT | 91.3% → 88.3% | 68.1% → 24.0% |
| COLO829T Illumina | 91.9% → 88.9% | 57.7% → 23.0% |

On HG008T, the truth alleles the tensors miss are mostly PoN alleles anyway (pon_tables §5b):
- HG008T SNV no candidate: 77–88% are PoN-tagged, vs 6.7–7.2% of the built ones. S3 and S4b are 91–100% tagged.
- INDEL I1 is 95.8–98.1% PoN-tagged on HG008T and 96.2–99.4% on COLO829T.
- After the PoN the missing-but-untagged SNV truth is small. For HG008T PacBio, 454 SNVs are missing; 357 of them are
  tagged.

COLO829T SNV misses are different: only 20–25% of the no-candidate SNVs are PoN-tagged (S3 still 77–98%). Most are
low-VAF (S2) or Extreme-region (S1) truth, so they stay missing after the PoN.

**PoN effect on the tensors** (chr1–22, % of placed tensors tagged by the rule; pon_tables §6, chr1 there too):

| Set | SNV label 2 | SNV label 0 | SNV label 1 | INDEL label 2 | INDEL label 0 | INDEL label 1 |
|---|---:|---:|---:|---:|---:|---:|
| HG008T PacBio | 96.4% | 60.4% | 7.7% | 97.7% | 44.7% | 73.7% |
| HG008T ONT | 95.8% | 37.1% | 7.9% | 96.7% | 37.9% | 70.0% |
| HG008T Illumina | 96.4% | 18.1% | 12.8% | 94.7% | 55.1% | 70.2% |
| COLO829T fiberseq | 96.1% | 55.0% | 3.4% | 97.4% | 54.2% | 62.8% |
| COLO829T ONT | 95.5% | 28.7% | 3.4% | 96.2% | 25.8% | 57.0% |
| COLO829T Illumina | 96.6% | 21.7% | 3.5% | 94.9% | 59.0% | 59.8% |

- **Germline (label 2):** the PoN removes 94.7–97.7% of germline tensors on every platform.
- **Label 0:** the tagged share depends on the platform:
  - SNV: 37–60% on the long reads vs 18–22% on Illumina. Illumina's SNV label 0 is mostly low-AF sequencing errors
    (`analysis/hg008_somatic_miss_20260926`), which no PoN contains.
  - INDEL: 26–59%.
- **Label 1 SNV:** the partial labels (`residual_partial_somatic_truth`) are 23–51% tagged. This fits the finding that
  partial SNV labels are germline-like.
- **Unplaced:** label-1 tensors without a GRCh38 allele (off-reference) cannot be reached by the PoN: SNV 0.3–6.1%,
  INDEL 8.7–19.8%.

**INDEL PoN rules compared** (`pon/indel_rules/README.md`, 2026-10-01):
- Other callers' rules on our data:
  - ClairS-TO (position matching for 1000G/CoLoRSdb) tags even more truth INDELs (87% / 73%).
  - DeepSomatic-like exact variant keys tag 19% / 23%.
  - Mutect2-like (1000G PoN + gnomAD AF ≥ 0.01) tags 7% / 9%.
- The overlap is real: 84% / 72% of the truth INDEL alleles are exactly in gnomAD or CoLoRSdb (homopolymer/STR alleles).
- A population-AF floor of 0.01 on the exact allele (gnomAD or CoLoRSdb) gives the best chr1 INDEL F1:
  - PASS 0.186 vs 0.094 with no PoN and 0.114 with the repo rule;
  - best 0.235 vs 0.150 / 0.126;
  - it costs recall: the ceiling drops to 28–43%.
- Chosen on chr1, the test chromosome. The genome-wide AF distributions support the 0.01–0.05 region.

## Files

| Path | Content |
|---|---|
| `tables.md` | all recall tables (§1–§8), from `recall_tables.py` |
| `recall_tables.py` | builds `tables.md`, `per_truth_HG008T.tsv`, `per_truth_COLO829T.tsv` |
| `finalize_job.sh` | Slurm job that reruns `recall_tables.py` and `pon/pon_tables.py` (used after the last COLO829T ONT report, 379062) |
| `per_truth_<sample>.tsv` | one row per truth allele (COLO829T with VAF_Ill / VAF_PB / RGN). Per platform: status, miss class, label of its tensor, ALT reads and AF of a filtered candidate |
| `collect_labels.py`, `labels/<set>.tsv` | label of the tensor that represents each built truth allele (Slurm 378983) |
| `filtered_truth/` | `patterns.py` + `filtered_scan.sh` (Slurm 378981): the builder's `filtered_candidates.ndjson` records of the filtered truth keys (`*.hits.ndjson`). The `*.patterns` inputs are not in git; `patterns.py` rewrites them |
| `colo829t_miss/scripts/` | read-level miss analysis for COLO829T (`reads.py`, `snv_detail.py`, `snv_nodes.py`, `report.py` = the HG008T scripts with `SAMPLE` as a parameter; job scripts) |
| `colo829t_miss/<platform>/` | per platform: `indel_no_candidate.tsv`, `snv_no_candidate.tsv`, `controls.tsv`, `residual_edits_labelled_germline.tsv`, `summary.json`, the raw `reads_part*.jsonl` and `snv_detail.jsonl`, job logs |
| `colo829t_miss/jobs.txt` | Slurm job ids |
| `pon/pon_scan.py`, `pon/pon_scan_job.sh` | PoN overlap per chromosome (Slurm array 379003) → `pon/scan/<chrom>.json`. `chr1.keys.tsv` (34 MB, input of the check below) is not in git; `DUMP_KEYS=1` rewrites it |
| `pon/check_vs_vcfeval.py` | the cross-check against vcfeval PoN outputs |
| `pon/pon_tables.py` → `pon/pon_tables.md`, `pon/truth_pon_<sample>.tsv` | PoN tables; one row per truth allele: overlaps per PoN, gnomAD/CoLoRSdb AF, status per platform |

## Reproduce

All inputs are read-only.

1. `sbatch collect_labels_job.sh`
2. `sbatch filtered_truth/filtered_scan.sh`
3. COLO829T read-level classes, per platform:
   1. `sbatch --array=0-(N-1) … reads_job.sh`, then `snv_job.sh`, then `report_job.sh` (commands in the job-script
      headers).
   2. The class report for each platform needs `SAMPLE=COLO829T_<platform>`.
4. `sbatch --array=0-23 pon/pon_scan_job.sh`
5. `python recall_tables.py`
6. `python pon/pon_tables.py` (needs `per_truth_*.tsv` from step 5).

Resources seen:
- read-level jobs 37–54 GB per job (24–40 processes);
- class reports ≈ 3 GB;
- PoN scan ≤ 2.2 GB per chromosome;
- label scan 0.24 GB.
