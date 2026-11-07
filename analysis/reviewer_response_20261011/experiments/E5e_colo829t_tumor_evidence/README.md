# E5e: COLO829T tumor assembly and tumor reads at the normal-absent absorbed alleles (HPRC v2.1, R2.3)

## Question

E5d gives, for HG008T on HPRC v2.1 (d46), two kinds of tumor evidence at the absorbed truth INDELs that neither
matched-normal haplotype carries (1,351 loci): the HG008-T tumor assembly (a tumor contig differs from both normal
haplotypes at 99.7%, equals GRCh38 + truth at 85.3%) and the tumor reads (the expected tumor allele length is seen at
97.4% of the evaluable loci with Illumina, 97.5% with PacBio; random non-absorbed truth INDELs 95.7% / 96.4%). The
user's absorbed-allele table (R2.3, 2026-11-05) left both cells empty for COLO829T. E5e computes them for the 311
normal-absent COLO829T INDELs (and the tumor assembly also for the normal-absent SNVs).

**The COLO829T tumor assembly is built from the tumor's own reads** (hifiasm of COLO829T HiFi + ONT, the assembly of
R3.7 and E9c). The truth set was called from tumor reads too, so a tumor contig that equals the truth allele shows
consistency between the assembly and the reads, not independent support. The HG008-T v3.2 assembly of E5d is also an
assembly of HG008-T reads.

## Data

| what | where |
|---|---|
| alleles | E5d COLO829T d46 set (401 = 54 SNV + 347 INDEL) with E5d's normal category (`$E5D/per_allele.tsv`): absent 311 INDEL + 52 SNV (E5d after its step 5, 2026-11-05; 319 before) |
| patient frame | `analysis/graph_absorbed_somatic_colo829t_20261005` c3 `normal_frame.tsv` (data dir): anchors, R / A, COLO829BL hapX / hapY sequences (whole tandem array between common unique anchors; s8b port), event hap, normal-frame event length; every chr1-22 truth INDEL and the 67 SNVs of 2026-10-05 (49 of the 52 absent SNVs) |
| tumor assembly | `E9c_offref_filter_validation/asm_targets/COLO829T_hifiasm_hap{1,2}.fa.gz` |
| normal dipcall (reads step) | `/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.vcf.gz` (GT = hapY\|hapX), cleaned like the HG008-N copy of the audit: records with only ACGTN alleles (drops the 52,500 records with `*`), bcftools 1.18 `norm -f GRCh38 -m -any`, sort, index -> `$DATA/COLO829BL_dipcall.clean.norm.sorted.vcf.gz` (5,198,467 records) |
| tumor reads | linear GRCh38 BAMs: Illumina `/scratch/qfu/COLO829T_WGS/SMHTCOLO829T-X-X-M45-A001-...-sentieon_bwamem_202308.01_GRCh38` (read only), Fiber-seq `COLO829T_fiberseq/SMHTCOLO829T-X-X-M45-B004-uwsc-MERGED-pbmm2_1.13.0_GRCh38`, ONT `COLO829T_ONT/COLO829T_ONT_merged.bam` (minimap2 2.26 map-ont) |
| controls | 400 random COLO829T chr1-22 truth INDELs of `normal_frame.tsv` in neither the d9 absorbed set nor the d46 set (`random.Random(9).sample` over the truth_id-sorted pool of 1,410) |

Large data: `/scratch/jshen/data/reviewer_experiments_20261013/E5e_colo829t_tumor_evidence/` (`$DATA`). Jobs and logs:
`/scratch/jshen/Github/Pansoma/tmp/e5e_colo829t_20261105/` (`jobs.txt`).

## Method

1. `e1_loci.py`: `$DATA/loci.tsv` (311 absent INDELs, 52 absent SNVs, 49 of them in normal_frame, 400 controls). Event
   hap for the audit's dipcall code: hapY -> GT index 0 ('hap1'), hapX -> index 1 ('hap2').
2. `e2_tumor_asm.py` (Slurm array, 8 chunks): the R3.7 / E27 `altlen_tumor_asm.py` search on every locus. The c3
   anchors (24-mers ending at `anchor_a0`, starting at `anchor_b0`) are searched exactly in every tumor contig, both
   strands; T = the tumor bases between them, oriented to GRCh38. Per locus every hit is a tumor copy.
3. `e3_tumor_reads.py`: the audit's `b2_tumor_reads.py` `work()` imported unchanged, with its DIP / BAMS paths set to the
   COLO829 inputs; the b4 rule from E5d `e3_hg008t_tumor_reads.b4_row` (imported): EXP = NE_len + nf_len (nf_len = c3
   `pf_len_change`, the normal-frame event length, as the audit's b5 `nf_len`); EXP seen at >= max(3, 15 %) of the reads;
   loci with >= 10 reads (MAPQ >= 20, <= 300) and a length-changing event.
4. `e4_summary.py`: tumor verdict = E5d s8b step 6 line for line (each copy vs hapX / hapY: hom / hapX / hapY / novel;
   novel_ALT / novel_other / no_change_* / all_eq_* / NA); E5d's HG008T columns are novel_ALT + novel_other ("in neither
   normal hap") and novel_ALT ("= GRCh38 + truth") over all normal-absent INDELs. Also the R3.7 per-tumor-hap outcome,
   the 22-locus check, and a strict reads variant (EXP differs from both normal lengths) for both tumors.

Difference from E5d's HG008T assembly step: E5d places each tumor contig by minimap2 (s6) and requires the anchors once
in the placed segment; here every exact anchor pair in the tumor assembly is a copy (as in R3.7), so a paralogous copy
can add a 'novel' copy. The table reports the loci with at most one copy per tumor hap separately.

## Results

`tables.md` and `summary.json` hold every number; per-locus verdicts in `$DATA/per_locus.tsv`.

**Tumor assembly, 311 normal-absent INDELs** (E5d HG008T: 1,342 / 1,346 = 99.7% and 1,148 / 1,346 = 85.3%):

| COLO829T, d46 | n | tumor copy in neither normal hap | tumor copy = GRCh38 + truth |
|---|---:|---|---|
| INDEL | 311 | 258 (83.0%) | 226 (72.7%) |
| SNV (49 of the 52 normal-absent SNVs; 3 lack c3 anchors) | 49 | 41 (83.7%) | 34 (69.4%) |

INDEL verdicts: novel_ALT 226, novel_other 32, no_change_hom 38, all_eq_hapX 10, all_eq_hapY 2, no_change_both 2, NA 1.
Loci with at most one copy per tumor hap (296): 82.1% / 72.0%, so paralogous copies do not drive the numbers.
Check: at the 22 ALT-length loci of E27 v3 `colo829t_altlen_tumor_asm` the 51 hit rows are identical and all 22 outcomes
reproduce (computed independently of the E5d class). After E5d step 5, 7 of the 22 are present (the normal hap carries
the truth plus another germline difference; the tumor assembly carries only normal sequences there, as expected); of the
15 still absent the tumor assembly carries the truth allele at 7, only normal sequences at 4, and is unresolved at 4.

At the 52 loci where the tumor assembly carries only normal-haplotype sequence (38 of them HP>=7), the tumor reads still
show the expected tumor allele (Illumina 45 / 45, Fiber-seq 46 / 47 evaluable loci): the assembly consensus misses
these somatic alleles (homopolymers, subclonal fractions); it does not show them to be germline.

**Tumor reads, 311 normal-absent INDELs** (audit b2 / b4 rule; HG008T E5d: Illumina 1,236 / 1,269 = 97.4%, PacBio
1,256 / 1,288 = 97.5%, controls 95.7% / 96.4%):

| platform | absorbed: expected allele seen / evaluable loci | 400 random non-absorbed truth INDELs |
|---|---|---|
| Illumina | 282 / 289 (97.6%) | 365 / 371 (98.4%) |
| Fiber-seq | 289 / 297 (97.3%) | 366 / 375 (97.6%) |
| ONT | 237 / 297 (79.8%) | 308 / 375 (82.1%) |

Not evaluable: Illumina 12 loci < 10 reads (MAPQ >= 20) and 10 without a length-changing normal-frame event (Fiber-seq /
ONT 3 and 11). Strict variant (EXP differs from both normal dipcall lengths): COLO829T 97.5% / 97.3% / 79.4%, controls
98.4% / 97.6% / 81.9%; HG008T 97.4% / 97.5%, controls 95.6% / 96.4% (E5d rows) - the rule's germline-length overlap
changes nothing material. ONT is lower for absorbed and control loci alike (homopolymer length spread; 218 of the 311
are HP>=7), so it measures the platform, not the alleles.

## Caveats

- Consistency, not independent support: the COLO829T assembly and the truth set both come from COLO829T reads.
- EXP can equal a germline length of the other normal haplotype (e.g. a +1 somatic change on a hap whose partner is
  already +1); then a tumor read of that length does not show the somatic allele. The E5d rule does not exclude these
  loci; the strict columns do, for both tumors.
- ONT reads carry homopolymer-length errors; at HP>=7 loci (218 of the 311) the 15 % threshold can be passed by error
  tails or missed by spread-out lengths. Illumina reads are often below MAPQ 20 in segmental duplications.

## Status

Done 2026-11-05. e2 470987 (array 0-7, meow; 22-24 min each, MaxRSS 2.8-2.9 GB), e3 470988 (miao, 12 workers, 10 min,
2.8 GB), e1 / e4 on the login node. e1 / e4 rerun 2026-11-05 after E5d step 5 (e2 / e3 outputs
cover the 319 old loci, a superset; e4 now restricts the read tables to the current loci). Not committed.
