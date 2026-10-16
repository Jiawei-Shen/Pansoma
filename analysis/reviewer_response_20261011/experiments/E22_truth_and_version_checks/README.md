# E22: CASTLE truth germline check and COLO829T d9 giraffe version

Reviewer points: check 1 serves **R2.2, R2.6 and R1.m1**; check 2 serves **R3.1 and R3.m1**.

Status (2026-10-15): both checks are done. Every number below comes from the outputs listed in "Files". All data were
opened read-only. Compute was light: one Slurm job (408124, `-p general -c 4`) made the bcftools subsets, and the
rest ran on the login node. No GAM was modified.

## Short answers

- **Check 1, H1437.** About 30 % of the H1437 CASTLE truth SNVs are probably germline in origin. They are not germline
  calls of the normal: 0 of 82,529 truth SNVs are an allele in the normal's genotype. Instead, they are alleles that
  the B-lymphoblastoid normal has lost through LOH. 22,982 truth SNVs (27.85 %) are gnomAD AF >= 1 % alleles, and
  97.5 % of the autosomal ones lie in two segments where the normal is almost entirely homozygous: chr9:68-139 Mb and
  chr1:9-12 Mb. In chr9:68-139 Mb the normal's het fraction is 0.061; it is 0.610 on the other autosomes. These two
  segments hold 24,914 truth SNVs (30.2 %). The common-allele mixture estimate gives about 24,900 (30.2 %).
- **Check 1, other lines.** H2009 has the same pattern on a smaller scale: 8,141 common SNVs (5.03 %), 87 % of the
  autosomal ones in chr7:122-143 Mb and three smaller normal-LOH segments. HCC1954 shows no normal LOH: 574 common
  SNVs (2.39 %), clustered at chr7:90-91 Mb and at the IGH, TRA and IGK loci. The SEQC2 HCC1395 control has 212
  (0.54 %).
- **Check 1, chr1 and training.** On chr1, the test chromosome, the effect is negligible for every line, for example
  19 common SNVs among 4,266 H1437 chr1 truth SNVs. The main effect is on training labels: about 10 % of H1437's
  somatic SNV training tensors lie in the normal-LOH segments.
- **Check 2.** The vg version cannot be read from any surviving log or script. The GAMs themselves fix it to a range:
  - COLO829T_fiberseq: vg **1.64.1–1.66.0**.
  - COLO829T_ONT: vg **1.64.0–1.66.0**.
  - COLO829T_Illumina: vg **≤ 1.62**. This is not 1.66.0, so it differs from the HG008T d9 GAMs. The exact version
    cannot be recovered.

---

## Check 1: do the CASTLE somatic truth sets contain germline variants?

### (a) Which truth files were used, and is the germline FILTER respected?

All three platform sets of a line use the same files, according to `tensors/SNV/labels.manifest.json` of
`/scratch/jshen/data/pansoma_v2_tensors/<LINE>T_{PacBio,ONT,Illumina}_d46`:

| line | somatic truth (labels) | derived from | germline truth (labels) | BEDs |
|---|---|---|---|---|
| H1437 | `CASTLE/H1437/truth/H1437_somatic_CASTLE.vcf.gz` (sha256 98fc46873dca…; 84,729 records loaded) | DeepSomatic `vcfs/somatic_only/H1437_…_somatic-only.vcf.gz`, PASS records (**somatic-only file**) | `H1437_germline_DeepVariant_hybrid_QUAL10.vcf.gz` (552d5ae01046…) | somatic: `H1437_high-confidence-regions_DeepSomatic-HCC1395-model.bed`; germline: `…_germline_high-confidence-regions_QUAL10.bed` |
| HCC1954 | `HCC1954_somatic_CASTLE.vcf.gz` (c74e5e629a79…; 25,004) | same, somatic-only PASS | `HCC1954_germline_DeepVariant_hybrid_QUAL10.vcf.gz` (f5bbdf6f0433…) | same pattern |
| H2009 | `H2009_somatic_CASTLE.vcf.gz` (c5b67d578f1c…; 161,057) | same, somatic-only PASS | `H2009_germline_DeepVariant_hybrid_QUAL10.vcf.gz` (b6a512adedac…) | same pattern |
| HCC1395 | `HCC1395_somatic_SEQC2_v1.2.1.vcf.gz` (fc8cfb9b6382…; 41,072) | SEQC2 v1.2.1 high-confidence sSNV + sINDEL (HighConf + MedConf) | `HCC1395BL_germline_DeepVariant_hybrid_QUAL10.vcf.gz` (43ed18e08ed0…) | somatic: `SEQC2_High-Confidence_Regions_v1.2.bed`; germline: `HCC1395_germline_high-confidence-regions_QUAL10.bed` |

(The record counts are the manifests' `truth.somatic.stats.records`, contigs of the graph only.)

- **The somatic truth is the somatic-only file.** Pansoma did not use the mixed file for somatic truth.
  `make_truth_castle.py` (`tmp/castle_hcc1395_20261006/`) keeps the PASS records of `somatic_only/`. Those files
  contain only PASS records (`filter_counts.tsv`): H1437 87,410, HCC1954 26,578, H2009 168,397.
- **The germline truth comes from the mixed file.** It is the `FILTER GERMLINE` records of the mixed benchmark VCF
  (`vcfs/<LINE>_DeepSomatic_HCC1395-model.vcf.gz`). Only the normal's genotype column is used, and only the ALTs in
  that genotype are kept. Since 2026-10-12, records with QUAL < 10 are moved out of the label VCF and their spans are
  removed from the germline BED. The mixed files contain only PASS and GERMLINE values: for H1437, 87,292 PASS and
  5,241,303 GERMLINE.
- **The GERMLINE FILTER is respected.**
  - No GERMLINE record enters the somatic truth.
  - The somatic-only file has 118 (H1437), 46 (HCC1954) and 217 (H2009) more PASS records than the mixed file. In
    each case the somatic allele sits at a site that bcftools merged into a GERMLINE record, but the normal's genotype
    carries a different allele at that site. The per-allele join counts 66 + 49 (H1437), 11 + 33 (HCC1954) and
    161 + 52 (H2009) such SNV + INDEL alleles. Keeping them as somatic is correct.
- **Few alleles are in both truth sets.** Only 7 somatic truth alleles are also an allele of the normal's genotype,
  all with QUAL >= 10. The label code puts somatic first, so these tensors are labelled somatic:
  - H1437: 2 INDELs.
  - HCC1954: 1 SNV, chr14:105734718 T>C. This is at the IGH locus, tumor VAF 0.785.
  - H2009: 3 SNVs, chr12:76402110, chr18:3608769 and chrX:232067, with tumor VAF 0.92-0.98 (classic tumor LOH of
    a normal het), and 1 INDEL.

### (b) Germline and common-allele content of the truth SNVs

Method (`truth_germline_check.py`):

- Every ALT of every record of the truth VCF is joined with three sources by exact allele (POS, REF, ALT):
  - the normal's GERMLINE calls in the mixed DeepSomatic VCF, at any QUAL;
  - the label germline VCF;
  - gnomAD and CoLoRSdb from `/scratch/jshen/data/Pansoma/panel_of_normal_VCFs`.
- In each normal's genome-wide germline SNVs (QUAL >= 10), the het fraction het / (het + hom-alt) is computed per
  1-Mb bin. A bin with a het fraction below 0.1 is homozygous in the normal, i.e. LOH in the normal cell line.
- Controls:
  - H2009 and HCC1395 (SEQC2), as asked;
  - two CASTLE truth sets that Pansoma does not use, HCC1395 and HCC1937 (DeepSomatic somatic-only). These show the
    behaviour of the same truth pipeline on other lines.

T1 (SNV alleles; `tables.md`):

| truth set | used | SNVs (all / in BED) | same allele in the normal's genotype (any QUAL / QUAL>=10) | normal call at same POS, other allele | gnomAD AF >= 1 %: all / in BED | CoLoRSdb AF >= 1 % |
|---|---|---|---|---|---|---|
| H1437 CASTLE | yes | 82,529 / 73,623 | 0 / 0 | 77 | **22,982 (27.85 %)** / 21,180 (28.77 %) | 23,491 (28.46 %) |
| HCC1954 CASTLE | yes | 24,033 / 20,846 | 1 / 1 | 17 | 574 (2.39 %) / 478 (2.29 %) | 693 (2.88 %) |
| H2009 CASTLE | yes | 161,917 / 144,176 | 3 / 3 | 185 | 8,141 (5.03 %) / 7,662 (5.31 %) | 8,244 (5.09 %) |
| HCC1395 SEQC2 | yes | 39,447 / 39,447 | 0 / 0 | 22 | 212 (0.54 %) / 212 (0.54 %) | 330 (0.84 %) |
| HCC1395 CASTLE | no | 120,941 / 35,040 | 2 / 2 | 66 | 72,264 (59.75 %) / 193 (0.55 %) | 73,784 (61.01 %) |
| HCC1937 CASTLE | no | 52,685 / 16,631 | 0 / 0 | 23 | 31,991 (60.72 %) / 271 (1.63 %) | 32,933 (62.51 %) |

Baseline: on chr20 inside the BED, 89.6-91.0 % of the H1437, HCC1954, H2009 and HCC1937 normals' germline SNVs are
gnomAD AF >= 1 % alleles. For example, H1437 has 67,919 of 74,719 (T4). HCC1395BL has 83.3 % inside the SEQC2 BED and
92.6 % inside the DeepSomatic BED. A true somatic set has about 0.5 % (SEQC2).

T2: where the common truth SNVs lie (autosomes):

| truth set | common SNVs | in bins where the normal's het fraction < 0.1 | normal-LOH segments >= 3 Mb (het < 0.2) | truth SNVs in those segments | tumor VAF median: common / rare |
|---|---|---|---|---|---|
| H1437 | 22,975 | 21,821 (94.98 %) | chr1:9-12, chr9:68-139 Mb | 24,914 (30.19 % of all truth SNVs), 22,397 of them common | 0.336 / 0.483 |
| H2009 | 8,122 | 6,852 (84.36 %) | chr5:47-51, chr7:55-62, chr7:64-67, chr7:122-143 Mb | 9,824 (6.07 %), 7,085 common | 0.674 / 0.501 |
| HCC1954 | 552 | 1 (0.18 %) | none | 0 | 0.478 / 0.243 |
| HCC1395 SEQC2 | 212 | 48 (22.64 %) | chr6:1-62, chr16:47-90 Mb and others | 207 (0.52 %) | 0.596 / 0.364 |

Notes on T2:

- The normal's het fraction is 0.061 in H1437 chr9:68-139 Mb (4,382 het / 67,904 hom-alt) and 0.057 in H2009
  chr7:122-143 Mb. It is 0.610-0.613 on the rest of the autosomes.
- Top bins of common SNVs:
  - H1437: all on chr9, for example chr9:92-93 Mb with 739.
  - H2009: all on chr7:122-142 Mb.
  - HCC1954: chr7:90-91 Mb (197); chr14:105-107 Mb, IGH (64); chr14:22-23 Mb, TRA/TRD (32); chr2:88-89 Mb, IGK (5).
- **[inference]** The normals are EBV B-lymphoblastoid lines. Lost germline haplotypes, or V(D)J deletions at the
  IG/TCR loci, make germline alleles look absent from the normal.

Other measures:

- **Mixture estimate (T4).** The germline-like share is (observed − 0.54 %) / (normal baseline − 0.54 %):
  - H1437: 30.2 %, about 24,942;
  - H2009: 5.0 %, about 8,047;
  - HCC1954: 2.1 %, about 500.
- **chr1, the test chromosome (T3).** Common truth SNVs are 19 of 4,266 (H1437), 33 of 11,792 (H2009), 27 of 1,596
  (HCC1954) and 9 of 3,440 (HCC1395). None of them lies in a normal-LOH bin.
- **Effect on Pansoma's labels** (`label_impact.tsv`, SNV, d46 sets):
  - **H1437.** 5,609 / 5,688 / 5,649 somatic-labelled tensors (PacBio / ONT / Illumina; 9.8-9.9 % of 56,937-57,371)
    carry a truth allele in the normal-LOH segments. 3,759-3,857 carry a common allele. On chr1 only 11-14 of about
    4,000 do.
  - **H1437 tensor coverage.** Only 22.4-22.7 % of the 24,914 segment truth SNVs get a tensor, compared with
    89.3-89.6 % of the truth SNVs outside the segments and common alleles. This explains H1437's ~69 % overall
    coverage (68.8 / 69.1 / 69.0 %). **[inference]** Most common alleles are paths in the HPRC v2.1 d46 graph, so
    reads carrying them align without a mismatch and make no candidate.
  - **H2009.** 3,099-3,167 somatic-labelled tensors (2.1-2.2 %) are in segments.
  - **HCC1954.** None are in segments, and 159-170 carry a common allele.
- **INDELs (T5).** The common-allele share is H1437 37.64 %, H2009 12.16 %, HCC1954 6.46 % and SEQC2 3.20 %. Repeat
  INDELs make this rate less specific than for SNVs.
- **DeepSomatic's own BEDs.** For HCC1395 and HCC1937, DeepSomatic's BEDs leave out the normal-LOH regions: within
  the BED, only 0.55 % and 1.63 % of truth SNVs are common, against 59.75 % and 60.72 % overall. The H1437 and H2009
  BEDs do not leave them out: 28.77 % and 5.31 % are common within the BED.

### (c) Conclusion

**H1437.** The truth plausibly contains a large block of germline-origin variants. These are not germline calls
that slipped past the FILTER, since no truth SNV is an allele in the normal's genotype. They are constitutional
heterozygous alleles that the matched normal cell line lost through LOH, mainly on chr9q (chr9:68-139 Mb), and that
the tumor retains. The tumor-normal pipeline therefore called them "somatic":

- 27.85 % of the truth SNVs are common gnomAD alleles, against 0.54 % for SEQC2.
- 97.5 % of the autosomal common ones lie in the two segments where the normal's het fraction is 0.06.
- Their tumor VAF (median 0.34) is that of retained heterozygous alleles.
- Both the segment count and the common-allele mixture estimate put the share at about 30 % (about 24,900 SNVs).
  This is an estimate, not a per-variant call: a minority of the segment SNVs may be real somatic mutations.

**H2009.** The same mechanism affects H2009 at about 5-6 % (chr7q).

**HCC1954.** This check gives no sign of a germline problem: 1 normal-genotype allele, 2.39 % common alleles that
cluster at chr7:90-91 Mb and the IG/TCR loci, and no normal-LOH segment. The E18 germline-like spectrum and the E4b
germline FPs therefore are not explained by germline alleles inside the HCC1954 truth.

**Tumor LOH.** Variants mislabelled through LOH in the tumor (normal het, tumor homozygous) are negligible: 3 SNVs in
H2009 and 1 in HCC1954.

**What this cannot detect.** Private germline variants that are absent from gnomAD and not called in the normal for
coverage reasons. That would need the normal BAMs, which are not local.

**Effect on Pansoma.** The chr1 test truth is essentially unaffected, at 0.45 % or less common alleles for H1437. The
effect is on chr2-22 training labels: about 10 % of H1437's somatic SNV tensors lie in these segments.

---

## Check 2: vg giraffe version of the COLO829T HPRC v1.1 (d9) GAMs

### The GAMs behind the tensor sets

Each tensor set reads a `.sorted.gam` (`run/config.json` `inputs.gam.path`). The sorted GAMs were made by
`vg_v1.77.0 gamsort` from the giraffe outputs (`/scratch/jshen/Pansoma_scripts/sort_gam_COLO829*.sh`, 2026-09).
Sorting does not realign reads.

| tensor set | sorted GAM | giraffe GAM (file time) |
|---|---|---|
| COLO829T_Illumina | `/scratch/jshen/data/COLO829T/illumina/GAM/COLO829T.sorted.gam` | `COLO829T.gam` (2025-05-09 10:43) |
| COLO829T_fiberseq | `/scratch/jshen/data/COLO829T/COLO829T_fiberseq/GAM/COLO829T_fiberseq.sorted.gam` | `COLO829T_fiberseq.gam` (2025-06-11 09:10) |
| COLO829T_ONT | `/scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.sorted.gam` | `COLO829T_ONT_merged.gam` (2025-09-04 18:50) |

### Logs and scripts: nothing records the binary

- `/scratch/jshen/Log_info/giraffe_COLO829T.log` (2025-05-08, empty) belongs to the Illumina job. Its `.err` was
  overwritten on 2025-08-26 by a later long-read job with the same name: it now contains reads SRR8955953.*.
- `COLO829T_fiberseq_AF-Filtered_pipeline_result.err` (2025-06-11) is the fiberseq giraffe stderr. Its read names
  (m84046_…) match the GAM, but it has no version line.
- No ONT giraffe log exists. `giraffe_COLO829T_ONT.err` (2025-09-04 21:15, empty) was written after the GAM, by the
  parse step.
- No script that writes these GAMs was found in `/scratch/jshen/Scripts`, `/scratch/jshen/Pansoma_scripts`,
  `~/.bash_history` or the Pansoma working tree.

### Evidence in the GAMs (read-only scan, `gam_tag_scan.py`, `gam_params.sh` → `vg_fingerprint.tsv`)

**How the fingerprint works.** Giraffe of vg >= 1.63 writes a `PARAMS_JSON` group into the GAM.
`vg_version_fingerprint.sh` mapped simulated reads on a 300-kb test graph with every local vg binary from 1.50.0 to
1.66.0, plus 1.77.0. Results:

- **Short-read default:** vg 1.50-1.62 write no PARAMS_JSON. vg 1.63.x, 1.64.0-1.65.0, 1.66.0 and 1.77.0 each write
  a distinct payload: md5 5ab31e48…, 8c9d0b57…, d261df87… and 9f454471….
- **`-b hifi`:** 1.64.0, 1.64.1, 1.65.0 and 1.66.0 all give fd3a0842…. 1.63.x and 1.77.0 differ.
- **`-b r10`:** 1.64.0-1.66.0 all give cf46e333…. 1.63.x and 1.77.0 differ.
- **Old-format minimizer index:** with a v9 index (no zipcodes), vg 1.63.x maps 0 of 4,000 reads ("MinimizerHeader:
  Expected v10, got v9"), and vg 1.64-1.66 abort.

**Reference GAMs.** These match their scripts' vg_v1.66.0 and so validate the method:

- HG008T Illumina d9: d261df87…, unique to 1.66.0.
- HG008T PacBio d9 and ONT d9: fd3a0842… and cf46e333….
- COLO829 ONT P25 mixture: cf46e333….

**COLO829T_fiberseq: vg 1.64.1, 1.65.0 or 1.66.0.**

- Its PARAMS_JSON is fd3a08426553, byte-identical to the 1.64.0-1.66.0 `-b hifi` payload and to HG008T PacBio d9.
- Its giraffe stderr contains "Cannot find an acceptable path from … with max path length of …" 5 times. That
  message string exists in vg 1.64.1, 1.65.0 and 1.66.0, but not in 1.63.1 or 1.64.0.

**COLO829T_ONT: vg 1.64.0-1.66.0.** Its PARAMS_JSON is cf46e333ee44, identical to the `-b r10` payload of
1.64.0-1.66.0 and to HG008T ONT d9.

**COLO829T_Illumina: vg ≤ 1.62. It is not 1.66.0.**

- **No PARAMS_JSON.** There is none in the first 3,000 groups (2,374,254 alignments). The 1.66.0 HG008T Illumina d9
  GAM has it at group 11, and every ≥1.63 test run writes it.
- **The index excludes ≥1.63.** On 2025-05-09 the only d9 minimizer index found was the 2023 HPRC
  `hprc-v1.1-mc-grch38.d9.min`, which has header v9. The search covered `/scratch/jshen/data` to depth 3, plus
  `AF-Filtered_VG_Indexes`. The v10 `shortread.withzip.min` appeared on 2025-06-07. ≥1.63
  cannot map with a v9 index, yet the GAM has mapped reads: 1,982 of the first 2,000 records have a path.
- **Same pattern in other runs.** The June 2025 Illumina d9 GAMs (COLO829BL.gam and SEQC2 SRR7890824.gam) also have
  no PARAMS_JSON.
- **Exact version unrecoverable.** The annotation keys are identical across 1.50-1.66. The candidates installed at
  the time are:
  - `ml load vg`, whose default module is 1.59.0;
  - `/scratch/wzhang/bin/vg`, which points to 1.61.0;
  - `vg_v1.43`-`vg_v1.62`.

### Context, not proof

- `/scratch/wzhang/bin/vg_v1.66.0` was installed on 2025-06-02, the d9 long-read indexes are dated 2025-06-05, and
  vg_v1.65.0 was also present.
- All of the user's later giraffe scripts use vg_v1.66.0: the COLO829 mixtures in 2025-10/11, HG008T in 2026-02/04
  and the linear runs in 2026-05.
- 1.66.0 is therefore the most likely long-read version, but the files only prove the ranges above.

To separate 1.64.1 / 1.65.0 / 1.66.0 one would re-map a few reads on the full d9 long-read indexes (~16 GB per
version) and compare alignments; this was not run under the light-compute rule.

**For the response.** Only the COLO829T Illumina d9 GAM differs in version (≤ 1.62) from the HG008T d9 GAMs (1.66.0).
The long-read COLO829T d9 GAMs share 1.66.0's parameter set exactly.

---

## Files

In this folder:

| file | what |
|---|---|
| `sets.tsv` | truth sets, BEDs, DeepSomatic VCFs, normal sample names |
| `subset_resources.sh` | bcftools subsets at the truth positions plus normal germline tables (ran as Slurm 408124) |
| `truth_germline_check.py`, `run_analysis.sh` | per-allele join, normal het-fraction bins, `summary.json` |
| `make_tables.py` | `tables.md` (T1-T5) and `summary_tables.json` |
| `label_impact.py` | `label_impact.tsv`: somatic-labelled tensors and tensor coverage by class |
| `filter_counts.sh` | `filter_counts.tsv`: FILTER values of the DeepSomatic mixed and somatic-only VCFs |
| `gam_tag_scan.py` | read-only GAM group/tag scanner |
| `vg_version_fingerprint.sh` | test-graph giraffe runs per vg version |
| `gam_params.sh` | `vg_fingerprint.tsv`: PARAMS_JSON of the real GAMs vs the test runs |

Intermediate files are in `/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks/`:

- `<set>.alleles.tsv.gz`: one row per truth allele;
- the resource subsets;
- `gam_params/`: payloads and scans;
- `vg_fingerprint/`: test graphs, reads, GAMs and logs per version.

Logs are in `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E22_truth_and_version_checks/logs/`.

To reproduce, run `bash subset_resources.sh --slurm`, which also runs `run_analysis.sh`. Then run
`python label_impact.py` and `bash filter_counts.sh`. For check 2, run `bash vg_version_fingerprint.sh` and then
`bash gam_params.sh`.
