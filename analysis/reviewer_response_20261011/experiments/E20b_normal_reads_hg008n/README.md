# E20b: do the matched normal's reads take the alternative graph path at the HG008T absorbed truth alleles?

Reviewer point: **R2.3** (Reviewer 2, major 3): "for false negatives at existing graph nodes, check whether the normal
sample's reads preferentially align to the alternative path". E20 answered it for COLO829T (COLO829BL ONT reads); the
user asked on 2026-11-05 to fill the HG008T column of the R2.3 absorbed-allele table the same way. HG008-N reads are not
local (the reason E20 was COLO829BL only, user 2026-10-21); they are read remotely from GIAB at the loci, as the R3.3
analysis does (`tmp/e9i_hg008t_20261104/fetch_normal.py`).

**Status (2026-11-05 15:40): done.** fetch 470970 (15.7 min, 0.9 GB), extract 470971 / 470972 (1-8 min), giraffe 470973 / 470974 (12 / 16 min, MaxRSS 41.1 / 41.7 GB), classify 470975 / 470976 (7 / 9 min, 6.7 / 8.3 GB); `report.py` and `control_snv.py` on the login node (`tables.md`, `summary.json`, `$D/control/snv_pileup.tsv`).

## Question

At a truth allele the graph absorbs (a d46 path spells GRCh38 + truth, so tumor reads that carry it walk that path with
no edit and give no candidate), where do the normal's reads go? If the graph pulls reads onto the alternative path,
normal reads whose own bases are GRCh38 (REF) would also walk the ALT path; if not, the normal's REF reads stay on GRCh38
and only reads that carry the allele take the branch.

## Data

| what | where |
|---|---|
| loci, HPRC v2.1 (d46) | the 1,607 HG008T absorbed alleles (140 SNVs, 1,467 INDELs) of E5d `hg008t_d46_absorbed_loci.tsv` (E5b's 'absorbed' rule, the set R1.4 / R2.3 quote: no candidate on d46 and a tensor on the GRCh38-only graph on >= 1 platform, and a d46 path spells GRCh38 + truth); truth records of `pansoma_v2_tensors/HG008T_PacBio_d46/truth/somatic.graph.tsv` |
| normal reads | HG008-N-D PacBio HiFi, GIAB BCM Revio 68x, GRCh38-GIABv3 pbmm2 BAM, read remotely over the regions (`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Revio_20240313/HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam`; local copy of its .bai) |
| tumor reads | HG008-T PacBio HiFi 116x `HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam` (the reads of Pansoma's HG008T PacBio GAMs) |
| graph | `AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.*` (GBZ, vg 1.77 dist, longread minimizers / zipcodes); graph index / GRCh38 path as the tensor sets |
| check | Pansoma's own HG008T PacBio HiFi d46 GAM, through the tumor counts E5c took from it (`E5c_d46_loss_reads/per_locus.tsv.gz`, set HG008T_PacBio_d46) |

The normal is PacBio HiFi, not ONT as for COLO829BL: GIAB's HG008-N HiFi BAM is the one the R3.3 analysis already reads
remotely, and HiFi has fewer homopolymer length errors, so single-read labels at INDELs are less noisy than in E20.
Normal and tumor are on the same platform and go through the same selection and alignment.

## Method

E20's steps; code reused unchanged where it is sample-independent.

1. `loci.py`: loci table and `regions.bed` (every locus +- 500 bp, merged: 1,580 regions, 1.58 Mb).
2. `fetch_normal.py` (`fetch.sh`): HG008-N records overlapping the regions, one remote connection with retries, each
   record once -> `normal_subset.bam`. `extract.sh`: primary records (`-F 0x900`) overlapping the regions, each read
   once (`samtools view -M -L`), to FASTQ, from the normal subset and from the tumor BAM.
3. `giraffe.sh`: both FASTQs to d46 with the settings of Pansoma's HG008T PacBio d46 GAM
   (`/scratch/jshen/Scripts/hprc_v21_d46_vg177_giraffe.sh` task 1, job 379456: vg 1.77, `-b hifi`, longread minimizers +
   zipcodes), then `vg gamsort -i`.
4. `classify.sh`: E20's `classify.py` unchanged (`TEST_DIR` = the E20b data dir), stream mode: per locus, every record
   with MAPQ > 10 that spans the allele (matched GRCh38 bases within 2 kb on both sides) gets a sequence label (its own
   bases: REF-like / ALT-like / unclear) and a path label (GRCh38 nodes / a path spelling the ALT / closer to the ALT /
   other).
5. `report.py`: E20's pooled definitions (imported), normal against tumor, per kind, per E5d assembly class of the
   allele (present / absent / ambiguous in the HG008-N haplotypes), and the realigned-tumor vs GAM check.

## Results

See `tables.md` for every number. Normal = HG008-N-D PacBio HiFi (119,141 reads extracted), tumor = HG008-T PacBio
HiFi (141,730), both realigned to d46; 1,601 / 1,602 of the 1,607 loci have a spanning read, median 72 (normal) and
86 (tumor) spanning reads per locus.

| | normal (HG008-N) | tumor (HG008-T) |
|---|---|---|
| REF-like reads on the ALT path, all loci | 427 / 58,736 (0.7%) | 345 / 38,545 (0.9%) |
| same, SNV / INDEL | 27 / 3,559 (0.8%) / 400 / 55,177 (0.7%) | 32 / 1,555 (2.1%) / 313 / 36,990 (0.8%) |
| reads carrying the ALT, SNV loci | 2,758 / 9,805 (28.1%) | 7,317 / 11,707 (62.5%) |
| ALT-like reads, INDEL loci | 41,449 / 104,295 (39.7%) | 82,105 / 124,752 (65.8%) |
| ALT reads on the ALT path | 42,753 / 44,207 (96.7%) | 87,320 / 89,422 (97.6%) |

- **The graph does not pull REF reads onto the ALT path** (0.7% of REF-like normal reads; COLO829BL ONT 1.1%, E20),
  while 96.7% of ALT-carrying reads take it.
- **The normal reads carry the ALT at many absorbed SNVs, unlike COLO829BL (2.2%); this is mostly how the HG008T truth
  is represented, plus a few SNVs E5d had misclassified.** Per locus (>= 5 reads, graph labels): normal ALT fraction
  < 0.1 at 61 of the 140 SNVs. The linear pileup agrees (`control_snv.py`, MAPQ >= 20: 23.9% normal, 65.5% tumor), and
  random non-absorbed truth SNVs show 0.2% in the normal (293 of 298 loci < 0.1), so it is not tumor-in-normal
  contamination. Two causes (`snv_event_frame.py`, `snv_event_frame.md`, 2026-11-05):
  1. E5d called a normal haplotype ALT only when its whole array equalled GRCh38 + truth; a haplotype with the ALT plus
     another germline difference in the same array counted as not carrying it. E5d step 5 (`e5_alt_with_germline.py`)
     fixed this: 5 SNVs moved to present and 10 to ambiguous (d46). By the corrected class: present 11 SNVs, normal
     47.1% ALT-like (heterozygous germline, as expected); ambiguous 14, 41.3%; absent 115, 1,997 / 8,073 (24.7%).
  2. The HG008T truth is called on the HG008-N assembly and projected to GRCh38. At the absent SNVs whose normal-frame
     event is an INDEL (41) the normal is 39.4% ALT-like: the event is a homopolymer / STR length change that appears on
     GRCh38 as an SNV, and the normal's own repeat length makes its reads match the projected base. At the 72 absent SNVs
     that are SNVs in the normal frame too, 16.0% of normal reads are ALT-like (pileup 8.9%, < 0.1 at 57 of 70 loci;
     tumor 61.9%); the 12 with pileup >= 0.1 all lie at repeat boundaries where HG008-N differs from GRCh38 (e.g. five
     GRCh38 SNVs at chr13:84,388,899-907 from one normal-frame T>G). The COLO829T truth was called from GRCh38 reads
     against the normal, so its SNVs do not have this projection.
  The whole-array comparison (E5d) is the measure of whether the normal carries the tumor allele; the per-read SNV
  fraction is not comparable between the two truth sets.
- INDELs: 39.7% of normal reads are ALT-like (tumor 65.8%); E5d present INDELs median 0.556 per locus, absent 0.257.
  HiFi length errors in homopolymers / STRs make single-read INDEL labels less reliable than SNV labels; no INDEL
  control was run.
- Check: realigned tumor reads against Pansoma's own HG008T PacBio d46 GAM (E5c counts), 1,289 loci: median |difference|
  of the ALT-read fraction 0.000, all 1,289 within 0.1; median spanning reads 87 and 87.

## Caveats

- Reads are selected by their GRCh38 (pbmm2) placement; a read that pbmm2 places elsewhere but giraffe would place at
  the locus is missed, for normal and tumor alike.
- Normal platform differs from E20 (HiFi here, ONT for COLO829BL); compare normal against tumor within a sample.
- The normal is HG008-N-D (duodenal) HiFi from GIAB; the HG008-N v6.2 assembly that E5d classifies against may come
  from other HG008-N data.
- `report.py` reads E5d's per-allele class at run time; rerun 2026-11-05 after E5d step 5 (class counts above).
- `control_snv.py` uses linear pileups (MAPQ >= 20, all base qualities), not the graph labels.

## Files

| file | what |
|---|---|
| `loci.py`, `loci_summary.json` | loci and regions |
| `fetch_normal.py` | remote HG008-N read fetch |
| `control_snv.py` | linear pileup at the 140 absorbed and 300 random truth SNVs (normal remote, tumor local) |
| `report.py`, `tables.md`, `summary.json` | tables |
| `snv_event_frame.py`, `snv_event_frame.md` | the 140 SNVs by E5d class and normal-frame event kind, graph labels and linear pileup |
| `../E20_normal_reads_colo829bl/classify.py` | per-locus read labels (unchanged) |
| `/scratch/jshen/Github/Pansoma/tmp/e20b_hg008n_20261105/` | `fetch.sh`, `extract.sh`, `giraffe.sh`, `classify.sh`, `jobs.txt`, logs |
| `/scratch/jshen/data/reviewer_experiments_20261013/E20b_normal_reads_hg008n/` | `loci_d46.tsv`, `regions.bed`, `idx/`, `normal_subset.bam`, `fastq/`, `gam/`, `classify/*.jsonl`, `per_locus.tsv` (not in git) |
