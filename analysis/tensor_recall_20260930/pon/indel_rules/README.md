# INDEL PoN: how other callers do it, and which rule works on our data (2026-10-01)

**Question.** The repo PoN rule tags 83% of the HG008T and 71% of the COLO829T somatic truth INDELs
(`../pon_tables.md` §1). Is that because of how we match, and would the rules of other tumor-only callers do better?

**Short answer.** The matching is not the cause: 84% (HG008T) and 72% (COLO829T) of the truth INDEL alleles are in
gnomAD or CoLoRSdb as the same allele. ClairS-TO's own rule tags even more. What helps is a **population-AF floor on
the exact allele**: AF ≥ 0.01 in gnomAD or CoLoRSdb, with no dbSNP and no position matching. On the chr1 calls it
roughly doubles the INDEL F1 (PASS 0.094 → 0.186, best 0.150 → 0.235); the repo rule gives 0.114 / 0.126.

## 1. What other callers do

Checked in their code and papers on 2026-10-01.

| Caller | PoN sources | Matching | AF cutoff | INDELs | Source |
|---|---|---|---|---|---|
| **ClairS-TO** (v0.4.4) | gnomAD r2.1 (GATK af-only), dbSNP 138 non-somatic, GATK 1000G PoN, CoLoRSdb v1.1.0 | exact allele for gnomAD and dbSNP; **position only** for 1000G and CoLoRSdb (`--panel_of_normals_require_allele_matching` default `True,True,False,False`) | gnomAD and CoLoRSdb sites with AF ≥ 0.001 | yes, same tagging. Also `LowSeqEntropy` for INDELs in low-entropy sequence; `--indel_min_af` 0.1 | README options; `src/nonsomatic_tagging.py` (`apply_one`: exact `POS REF ALT` key or `POS` only); paper (Nat Commun 2025): "exact allele matching … for gnomAD and dbSNP, and positional matching for 1000 G PoN and CoLoRSdb"; COLO829 INDEL recall "capped at ~50%", AUPRC 0.20 |
| **DeepSomatic** tumor-only | one merged PoN per platform: `PON_dbsnp138_gnomad_ILMN1000g_pon.vcf.gz` (Illumina) or `PON_dbsnp138_gnomad_PB1000g_pon.vcf.gz` (PacBio, ONT), shipped in the Docker image | exact variant key: contig, start, REF and the **whole ALT list** (`postprocess_variants.should_filter`) | not stated in the code (the PoN files are built in advance) | yes, every PASS call | `scripts/run_deepsomatic.py`, `deepvariant/postprocess_variants.py`; docs: "PON vcf that contains variant calls from dbSNP, gnomAD and 1000 genomes"; issue #25 (developer): "For PON of DeepSomatic, we use dbSNP, clinvar and 1KGP variants". **CoLoRSdb is not a hard filter there**: the long-read models read `AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz` (CoLoRSdb alleles at AF ≥ 0.05) through `--population_vcfs`, i.e. as the population-AF input channel of the network (issue #44 command line). Published tumor-only INDEL results with the default PoN (HCC1395, `docs/metrics.md`): recall 0.848 / 0.758 / 0.574, F1 0.511 / 0.640 / 0.502 (Illumina WGS / PacBio / ONT) |
| **Mutect2** (GATK) | PoN = sites seen in normals (the 1000G PoN we have is the GATK one); gnomAD = germline resource | PoN: hard filter; gnomAD: prior on the germline hypothesis, combined with the tumor AF in `FilterMutectCalls` | none (probabilistic) | yes | GATK best practices. Not re-read here: the GATK docs page returned 403, so this row is from general knowledge |

Our runs of these rules use our four files (`/scratch/jshen/data/Pansoma/panel_of_normal_VCFs`). DeepSomatic's own PoN
files are not public outside its image, so the "DeepSomatic-like" rule (exact variant key in dbSNP, gnomAD or 1000G) is
an approximation. Mutect2's model is approximated by "1000G PoN + gnomAD AF ≥ 0.01".

## 2. Results

- Rule definitions: `pon_features.py`.
- Truth and tensor columns: chr2–22 (`rules_tables.md` §1–3; chr1 there too).
- chr1 call columns: `HG008_Illumina_INDEL_nopartial_b1024`, chr1 in the BED, 587 truth records (`calls_chr1.md`; the
  `_nopc` run ranks the same).

| Rule | truth tagged HG008T / COLO829T | INDEL recall ceiling after: HG008T Illumina / PacBio / COLO fiberseq | germline (label 2) removed: HG008T PacBio / Illumina | chr1 calls: PASS F1 | chr1 calls: best F1 |
|---|---:|---:|---:|---:|---:|
| no PoN (current INDEL setting) | 0 / 0 | 43.0 / 52.5 / 66.4% | 0 / 0 | 0.094 | 0.150 |
| repo rule (= SNV rule) | 83.1 / 70.5% | 12.8 / 13.1 / 24.0% | 97.6 / 94.7% | 0.114 | 0.126 |
| ClairS-TO | 86.9 / 73.1% | 11.6 / 11.7 / 23.0% | 98.5 / 96.2% | 0.123 | 0.131 |
| DeepSomatic-like | 18.5 / 22.6% | 36.0 / 42.6 / 55.1% | 56.9 / 53.1% | 0.134 | 0.178 |
| Mutect2-like | 7.3 / 8.8% | 40.3 / 49.3 / 61.5% | 18.4 / 32.0% | 0.118 | 0.164 |
| repo rule without CoLoRSdb | 33.5 / 36.0% | 31.1 / 36.1 / 47.0% | 74.5 / 63.8% | 0.142 | 0.180 |
| population AF ≥ 1e-3 (gnomAD or CoLoRSdb allele) | 79.9 / 66.7% | 14.9 / 15.3 / 27.0% | 96.8 / 92.8% | 0.127 | 0.142 |
| **population AF ≥ 0.01** | 56.9 / 46.8% | 27.9 / 29.7 / 42.4% | 90.1 / 82.0% | **0.186** | **0.235** |
| population AF ≥ 0.05 | 26.1 / 21.3% | 40.3 / 45.1 / 59.5% | 71.7 / 52.6% | 0.154 | 0.215 |
| population AF ≥ 0.1 | 14.2 / 12.7% | 42.5 / 49.0 / 63.1% | 57.8 / 30.0% | 0.120 | 0.177 |

Adding the 1000G PoN and dbSNP `COMMON` (≥ 1% in a 1000G population) to a floor changes almost nothing (0.186 / 0.232
at 0.01).

### Why: population AF of the exact allele (max of gnomAD, CoLoRSdb), chr1–22 (`rules_tables.md` §4)

| | not in either | < 1e-3 | 1e-3 – 0.01 | 0.01 – 0.05 | 0.05 – 0.1 | ≥ 0.1 |
|---|---:|---:|---:|---:|---:|---:|
| HG008T somatic truth INDELs | 16.1% | 3.8% | 23.0% | 30.7% | 12.0% | 14.4% |
| COLO829T somatic truth INDELs | 27.9% | 5.0% | 19.8% | 25.3% | 8.7% | 13.2% |
| HG008T PacBio germline INDEL tensors | 2.3% | 0.9% | 6.7% | 18.3% | 13.9% | 58.0% |
| HG008T Illumina germline INDEL tensors | 5.1% | 2.1% | 10.8% | 29.3% | 22.7% | 30.1% |

- **Somatic INDELs are population alleles.** Most truth INDELs are homopolymer/STR length changes, which recur at the
  same alleles that are polymorphic in the population. On chr1, `analysis/indel_20260930` finds:
  - 74% of the truth INDELs in ≥ 5 bp homopolymers;
  - 43% "composite": on GRCh38 the somatic step is written together with the patient's own STR length difference.
- **Matching tighter does not help.** The overlap is with the exact allele, so trimming REF/ALT or position matching
  barely changes it (`../pon_tables.md` §1–2). Position matching (ClairS-TO) only makes it worse.
- **Where the truth and germline differ.** The truth sits mostly at AF 1e-3–0.05; the individual's germline INDELs
  mostly at ≥ 0.05.
  - A low floor (1e-4, 1e-3: the repo rule, ClairS-TO) removes both.
  - A floor of 0.01–0.05 keeps most of the somatic INDELs at AF 1e-3–0.01, and still removes 69–90% of the long-read
    germline.
  - dbSNP has no AF in this file and contributes STR alleles at any frequency.
- **COLO829T** has the same shape, so this is not specific to HG008T.

## 2b. DeepSomatic on COLO829T: the same problem, called `GERMLINE` by the model (`deepsomatic_colo829t.py`)

A DeepSomatic tumor-only run already on disk: v1.9.0, `PACBIO_TUMOR_ONLY` model on the COLO829T fiberseq BAM, default
PoN, chr1 only (`/scratch/jshen/data/COLO829T/COLO829T_fiberseq/deeepsomatic_TO_results/`, 2025-08). This is where the
COLO829T chr1 truth alleles end up:

| chr1 truth | PASS | GERMLINE (model) | RefCall | PON | other / no record |
|---|---:|---:|---:|---:|---:|
| INDEL (131) | 12 (9.2%) | 72 (55.0%) | 36 (27.5%) | 0 | 11 |
| SNV (2,443) | 1,378 (56.4%) | 854 (35.0%) | 102 (4.2%) | 16 (0.7%) | 93 |

- **INDELs:** DeepSomatic's hard PoN removes no truth INDEL. Its model labels 55% of them `GERMLINE`, and 55 of those 72
  are alleles our PoN rule tags.
  - The GERMLINE truth INDELs: median population AF 0.018 (42 of 72 at ≥ 0.01), VAF 0.71.
  - The 12 PASS truth INDELs: 11 of 12 absent from gnomAD/CoLoRSdb, VAF 0.47.
  - So DeepSomatic loses the population-allele somatic INDELs too. Its population-AF input (CoLoRSdb AF ≥ 0.05) and the
    high VAF make the model call them germline. Its chr1 INDEL recall here is 9%.
- **SNVs:** the GERMLINE truth SNVs have VAF 1.00 and no population AF. These are COLO829's clonal LOH SNVs, a VAF
  effect, not a PoN one.
- **Caveat:** fiberseq is not the HiFi data the PacBio model was trained on, so treat DeepSomatic's absolute numbers on
  this run with care. The GERMLINE-vs-population pattern does not depend on that.

## 3. Caveats

- **Test-chromosome tuning.** The call-level columns come from chr1, the test chromosome. The INDEL runs keep no
  validation predictions (`val_fraction` 0.05, nothing saved), so the floor was compared on chr1.
  - The genome-wide distributions above (chr1–22, tensors and truth of both samples) show the same separation. That
    supports 0.01–0.05 as the region, but not one exact value.
  - Before using a floor, fix it once (e.g. from validation predictions) rather than from further chr1 results.
- **Recall cost.** Every floor still costs recall. At 0.01 the INDEL recall ceiling drops from 43–69% to 28–43%; at
  0.05 to 40–60%. Which is preferable depends on how much the downstream steps favour recall.
- **Approximated rules.** DeepSomatic's own PoN and Mutect2's germline model are approximated as described in §1.

## Files

| File | Content |
|---|---|
| `pon_features.py` | per-allele PoN features and the rules |
| `indel_scan.py`, `indel_scan_job.sh` | features of every truth INDEL and INDEL tensor, chr1–22 (Slurm array 379130, ≤ 0.6 GB, ≤ 3 min per chromosome) → `scan/chr*.json` |
| `rules_tables.py` → `rules_tables.md` | truth tagged, recall ceilings, tensors by label, population-AF distributions |
| `calls_chr1.py` → `calls_chr1.md`, `calls_chr1.tsv` | chr1 calls of the two HG008T Illumina INDEL runs under every rule (from their rtg TP/FP) |
| `deepsomatic_colo829t.py` → `deepsomatic_colo829t.txt` | COLO829T chr1 truth in the local DeepSomatic tumor-only run |

The repo rule rebuilt from the features equals `../pon_scan.py` (= `scripts/filter_panel_of_normals.py`) on chr22:
94 truth INDELs and 72,037 tagged tensors, 0 differences.
