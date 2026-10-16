# E11b: alignment correctness with simulated reads of known origin (R3.6, R1.3)

## Question

Are more perfect alignments on a pangenome graph also more correct alignments? Reads were simulated from the HG008-N
v6.2 diploid assembly, which is in none of the references, so every read's true GRCh38 position is known. The same reads
are aligned to three graph references (GRCh38-only, HPRC v1.1 d9, HPRC v2.1 d46). For each one we measure:

- the misplacement rate;
- "perfect but wrong" alignments;
- MAPQ calibration;
- unmapped reads.

Results are split by HG008-N germline variants in the read, repeat content, SMaHT region class and off-reference path
use.

Scope as run: Illumina, HiFi and ONT on all three arms. The ONT-like reads were first left out by a submission error
(see Status) and were aligned and scored on 2026-10-15 (jobs 407745-407748). The CHM13-only control arm was planned and
then dropped by the user on 2026-10-13.

## Reviewer points answered

- **R3.6.** More perfect alignment is not more correct alignment. Compare GRCh38, T2T-CHM13 and HPRC, use
  assembly-based truth and MAPQ. This experiment covers point 3 of the R3.6 response ("Correctness with known
  origin") for GRCh38 and HPRC. The T2T-CHM13 comparison is no longer part of it, because the CHM13 work was dropped.
- **R1.3.** Pangenome ambiguity in difficult regions and with noisy ONT reads. This experiment covers the simulated-read
  part of point 1 of the R1.3 response (misplacement and MAPQ calibration per stratum on each reference) for Illumina,
  HiFi and ONT. For ONT the answer is unfavourable to HPRC v2.1 (Interpretation 7).

## Data

| item | path |
|---|---|
| HG008-N v6.2 haplotypes | `/scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_polished6.2_hap1.fasta`, `..._hap2.fasta.gz` |
| hap -> GRCh38 truth | dipcall: `/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.hap{1,2}.paf.gz` (minimap2 2.20 `-c --paf-no-hit -xasm5 --cs`), `.hap{1,2}.bed` (paftools call R lines = GRCh38 covered by exactly one contig), `.hap{1,2}.var.gz` (per-haplotype variants) |
| GRCh38 | `/scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta` |
| region classes | SMaHT easy / difficult / extreme v2 (`/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_{difficult,extreme}_v2.bed`; easy = the rest; the same BEDs as the COLO829T truth). GIAB v3.5 stratifications are not on the server. |
| GRCh38-only graph ("linear", H 17) | `/scratch/jshen/data/Linear_Reference_VG_Indexes/hg38_linear.{giraffe.gbz,dist,shortread.withzip.min,shortread.zipcodes,longread.withzip.min,longread.zipcodes}` (vg 1.66 autoindex; minimizer format v10) |
| HPRC v1.1 d9 (H 90) | `/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.{gbz,dist,shortread.withzip.min,shortread.zipcodes,longread.withzip.min,longread.zipcodes}` (minimizer format v10) |
| HPRC v2.1 d46 (H 464) | `/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.{gbz,dist,shortread.withzip.min,...}` (vg 1.77-rebuilt dist / min, format v11) |
| off-reference nodes | the pipeline's `grch38_path/chrom.npy` per graph (`chrom == -1`: node not on the GRCh38 path) |
| simulated reads + truth | `/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/sim/{illumina,hifi,ont}/` (`sim_R1.fq.gz` + `sim_R2.fq.gz` or `sim.fq.gz`, `truth.npz`, `sim_meta.json`); segments in `sim/segments.tsv`, `sim/alignments.{tsv,npz}`, `sim/segments_meta.json` |
| alignments | `.../E11b_simulated_alignment/aln/{linear,d9,d46}/{illumina,hifi,ont}.{gam,gaf.gz,bam}` with `*.giraffe.log`, `*.surject.log` (`ont.*` written 2026-10-15 19:25-19:38) |
| per-read scores | `.../E11b_simulated_alignment/eval/<arm>.<platform>.{json,reads.npz}` (9 sets) |
| CHM13 truth (dropped control; built, unused) | `/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/chm13_truth/HG008N_CHM13_dipcall.hap{1,2}.*` (dipcall commands against `/scratch/qfu/reference/chm13/chm13v2.0_maskedY_rCRS.fa`); `segments.tsv` (350 segments, 583,229,106 bases); `illumina.chm13_truth.npz`, `hifi.chm13_truth.npz` (2026-10-13) and `ont.chm13_truth.npz` (written by collect 387349 on 2026-10-15) |
| CHM13-only graph (dropped control) | not built: `chm13_linear/` holds only the partial `tmp-autoindex/` of the cancelled job 387281 |

HG008 is in neither HPRC graph. HPRC v1.1 has no HG008, and the v2.1 sample list was checked on 2026-10-01
(memory note `hprc-v2-graph`).

## Method

1. **One-to-one segments** (`build_segments.py`). A segment is a haplotype interval with these properties:
   - it lies in one primary dipcall alignment of `chrN_hapK -> chrN`, with N = chr1 or chr20, MAPQ >= 5 and block
     >= 50 kb;
   - its GRCh38 image is inside one interval of `hapK.bed`, so the reference side is covered once;
   - no haplotype base is covered by a second primary or inversion alignment, so the query side is covered once;
   - it is at least 1 kb long.

   The CIGAR (`cg:Z`) gives a base-level hap -> GRCh38 map. For inserted bases, the map gives the insertion point.
   Result (`sim/segments_meta.json`): 358 segments, 570,972,967 haplotype bases (hap1 175 alignments,
   285,835,588 bases; hap2 172 alignments, 285,137,379 bases).
2. **Simulation** (`simulate.py`; seed 20261013). Templates are drawn uniformly over the segment bases of both
   haplotypes, with 50/50 strand and pure-ACGT templates only. A read's truth is its template's GRCh38 interval.
   - **Illumina:** 2,500,000 pairs (5,000,000 reads), 2 x 150, fragment N(400, 50) clipped to [250, 600], substitutions
     0.2 % per base, constant Q37 qualities. This simulator is our own: wgsim's read names do not say which mate comes
     from the left end.
   - **HiFi-like:** 40,000 reads, length N(15 kb, 3 kb) clipped to [5, 30] kb, 0.1 % errors (0.04 % substitutions,
     0.03 % 1-bp insertions, 0.03 % 1-bp deletions).
   - **ONT-like:** 40,000 templates, length N(15 kb, 4 kb) clipped to [3, 30] kb. Errors are added by pbsim3
     template mode with the R10.4 `ERRHMM-ONT-HQ` model (`--accuracy-mean 0.98`). The local pbsim3 is
     `/scratch/wzhang/bin/pbsim`. Qualities were set to a constant Q20.
     - pbsim3's own report (`sim/ont/pbsim/pbsim.log`): 40,000 reads, 600,367,967 template bases, read length mean
       14,932 (SD 3,963, min 2,871, max 29,950), read accuracy mean 0.9592 (SD 0.0470), substitution 1.13 %,
       insertion 1.20 %, deletion 1.72 %.
     - minimap2 NM / block length of the first 2,000 reads against their own templates (`results/ont_error_profile.txt`;
       1,801 reads with a primary hit to their own template): median 2.17 %, IQR 0.82-4.96 %, 90th percentile 8.60 %.
3. **Truth check** (`check_truth.py`, first 20,000 Illumina rows). For reads with no HG008-N indel or SV and no inserted
   bases, two things must hold:
   - the read length equals `ref_end - ref_start`;
   - the mismatches against GRCh38 equal the haplotype's SNVs in the read plus the simulated errors.
4. **Alignment** (`align_arm.sh`, Slurm `general`, 16 threads). Each arm was meant to use the vg version of the
   production GAMs for that graph:
   - linear: `/scratch/wzhang/bin/vg_v1.66.0`;
   - d9: `/scratch/wzhang/bin/vg_v1.65.0` (the job logs confirm `vg version v1.65.0 "Carfon"`, e.g.
     `align_E11b_align_d9_lr_387346.out`). Correction 2026-10-15: the HG008T production d9 GAMs were made with vg
     1.66.0, not 1.65.0 (their giraffe logs in `/scratch/jshen/Log_info/`, `vg version v1.66.0 "Navetta"`; see
     `../../points/R3.1.md`). The d9 arm here is therefore one minor version older than production. The COLO829T d9 GAMs'
     version is in no log; their embedded `PARAMS_JSON` fingerprint gives 1.64.1-1.66.0 (Fiber-seq), 1.64.0-1.66.0 (ONT)
     and <= 1.62 (Illumina, exact version unrecoverable) (`../E22_truth_and_version_checks/README.md`, check 2);
   - d46: `/scratch/jshen/bin/vg_v1.77.0`.

   ```
   vg giraffe -Z <gbz> -d <dist> -m <shortread|longread>.withzip.min -z <...>.zipcodes -f R1 -f R2 -t 16        # Illumina
   vg giraffe ... -f sim.fq.gz -b hifi | -b r10 -t 16                                                            # HiFi / ONT
   vg convert -t 16 -G <arm>/<platform>.gam <gbz> | gzip -1 > <platform>.gaf.gz                                  # graph view
   vg surject -x <gbz> [-n GRCh38 for d9/d46] [-D short|long for vg 1.77] [-i for pairs] -t 16 -b <gam> > <platform>.bam
   ```
5. **Scoring** (`evaluate.py`, primary records only).
   - *correct*: the surjected alignment is on the true chromosome, and its start or end is within 50 bp of the true
     interval.
   - *wrong*: placed on GRCh38 but not correct.
   - *graph-mapped, not surjected*: aligned in the graph but without a GRCh38 projection; counted neither correct nor
     wrong (so it lowers the correct share without raising the wrong share).
   - *perfect*: the graph alignment (GAF `cs`) covers the whole read with no mismatch or indel.
   - *MAPQ*: giraffe's graph MAPQ, the one the Pansoma pipeline reads.
   - *Calibration*: wrong / placed in each MAPQ bin, against the mean 10^(-MAPQ/10).
   - *Strata*, all from the truth:
     - HG008-N variants of the read's own haplotype inside its true interval (from `hapK.var.gz`; SV = >= 50 bp);
     - distance to the nearest such variant;
     - inserted (non-GRCh38) bases;
     - STR content: period 1-6, tract >= max(12, 3 x period);
     - homopolymer >= 7;
     - SMaHT class by majority of the read's GRCh38 bases;
     - for d9 / d46, whether the graph alignment uses a node that is off the GRCh38 path.
   - *Same reads, two arms*: reads correct in one arm and not the other, and confident wrong (MAPQ >= 30) only in one
     arm.
   - GAF header lines (`@HD VN:Z:1.0`, written by vg 1.77 `convert`) are skipped. This fix was made on 2026-10-15 after
     the first collect failed on the d46 GAF; the GRCh38-only and d9 GAFs (vg 1.66 / 1.65) have no header.
   - `summarize.py` turns the per-arm JSONs into `results/tables.md`, `results/summary.{tsv,json}` and the generated
     block of this README.
   - The scoring chain (`align_arm.sh` steps -> `evaluate.py` -> `summarize.py`) was run end to end on a 400-read test
     against a tiny autoindexed graph (`.../tinytest/base/`). That test only checks the code; its numbers mean nothing.
6. **CHM13-only control: dropped (user decision 2026-10-13).** The plan was to align the same reads to a CHM13-only
   linear graph and score them against the template lifted to CHM13. The hap -> CHM13 alignments (387282, 387283) and
   the CHM13 segments / truth had already been made. The graph index (387281) and the CHM13 alignment job (387348) were
   cancelled. `align_arm.sh`, `chm13_truth.py`, `collect.sh` and `summarize.py` still contain the `chm13` arm; with no
   `aln/chm13/` output it is skipped.

## Exact commands

Job scripts and logs: `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E11b_simulated_alignment/`; job ids in
`jobs.txt` there, Slurm logs in `logs/`.

```
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E11b_simulated_alignment
D=/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment
J=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E11b_simulated_alignment
python $E/build_segments.py --out $D/sim                                     # login node, 12 s
sbatch --export=ALL,PLATFORM=illumina,N=2500000 $J/sim.sbatch                # 387278
sbatch --export=ALL,PLATFORM=hifi,N=40000 $J/sim.sbatch                      # 387340
sbatch --export=ALL,PLATFORM=ont,N=40000 $J/sim.sbatch                       # 387341
sbatch $J/chm13_index.sbatch                                                 # 387281 (cancelled: CHM13 dropped)
sbatch --export=ALL,HAP=1 $J/chm13_paf.sbatch; sbatch --export=ALL,HAP=2 $J/chm13_paf.sbatch   # 387282, 387283
sbatch --mem=72G --export=ALL,ARM=linear,PLATFORMS=illumina $J/align.sbatch  # 387342 (d9 387343 at 72G, d46 387344 at 84G)
sbatch --mem=40G --export=ALL,ARM=linear,PLATFORMS=hifi,ont $J/align.sbatch  # 387345 (d9 387346, d46 387347); memory raised
                                                                             # at 11:40 to 64G / 64G / 72G (jobs.txt).
                                                                             # sbatch cut PLATFORMS at the comma: HiFi only
sbatch --mem=72G --export=ALL,ARM=chm13,PLATFORMS=illumina,hifi,ont $J/align.sbatch   # 387348 (cancelled: CHM13 dropped)
sbatch $J/collect.sbatch                                                     # 387349 = bash $E/collect.sh (afterany on 387342-387347); FAILED
sbatch -o logs/collect_%j.out collect.sbatch                                 # 406576, 2026-10-15 17:53, after the evaluate.py fix
sbatch --mem=24G --export=ALL,ARM=linear,PLATFORMS=ont -J E11b_align_linear_ont align.sbatch                 # 407745 (ONT, one platform per job)
sbatch --mem=36G --dependency=afterany:407745 --export=ALL,ARM=d9,PLATFORMS=ont -J E11b_align_d9_ont align.sbatch     # 407746
sbatch --mem=56G --dependency=afterany:407746 --export=ALL,ARM=d46,PLATFORMS=ont -J E11b_align_d46_ont align.sbatch   # 407747
sbatch --dependency=afterany:407747 -o logs/collect_%j.out collect.sbatch    # 407748 (scores only the sets without an eval/*.json)
```

Pre-alignment checks, run on the login node on 2026-10-13 (22:55-22:58) with the commands of `collect.sh`:

```
python $E/check_truth.py $D/sim/illumina /scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta > $E/results/truth_check.txt
gzip -dc $D/sim/ont/sim.fq.gz | head -n 8000 > $D/sim/ont/first2000.fq; minimap2 2.28 -c ... (collect.sh) > $E/results/ont_error_profile.txt
python $E/segment_stats.py $D/sim /scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta $E/results
python $E/summarize.py --evals $D/eval --out $E/results                                                            # 0 alignment sets then
```

`collect.sh` builds the CHM13 truth of every platform whenever both hap -> CHM13 PAFs exist, so collect 387349 wrote
`chm13_truth/ont.chm13_truth.npz` on 2026-10-15. It is unused.

## Job status (final, 2026-10-15)

| job | what | state |
|---|---|---|
| 387278 | simulate Illumina | COMPLETED |
| 387340 | simulate HiFi | COMPLETED |
| 387341 | simulate ONT (pbsim3) | COMPLETED |
| 387282, 387283 | hap1 / hap2 -> CHM13 dipcall alignments | COMPLETED (minimap2 peak RSS 27.8 / 28.3 GB) |
| 387281 | CHM13-only graph index | CANCELLED 2026-10-13 11:51 during graph construction (CHM13 dropped) |
| 387348 | align all platforms to CHM13 | CANCELLED (CHM13 dropped) |
| 387342 / 387343 / 387344 | align Illumina: linear / d9 / d46 (16 CPUs; 72G / 72G / 84G) | COMPLETED 2026-10-15 13:05-14:17 (15 min 51 s / 17 min 38 s / 38 min 36 s) |
| 387345 / 387346 / 387347 | align long reads: linear / d9 / d46 (16 CPUs; 64G / 64G / 72G) | COMPLETED 2026-10-15 13:05-13:38 (6 min 50 s / 8 min 38 s / 16 min 58 s), **HiFi only** |
| 387349 | collect (`collect.sh`) | FAILED 2026-10-15 14:22 in `evaluate.py` on the d46 Illumina GAF header; the linear and d9 Illumina evaluations it wrote are kept |
| 406576 | collect rerun after the `evaluate.py` fix | COMPLETED 2026-10-15 17:53-17:57; scored d46 Illumina and the three HiFi sets, rebuilt the tables (6 arm x platform sets) |
| 407745 / 407746 / 407747 | align ONT: linear / d9 / d46 (16 CPUs; 24G / 36G / 56G; chained) | COMPLETED 2026-10-15 19:23-19:38 (4 min 12 s / 4 min 21 s / 7 min 9 s; MaxRSS 16.3 / 24.2 / 39.5 GB); logs `align_E11b_align_{linear,d9,d46}_ont_40774{5,6,7}.out` confirm vg v1.66.0 / v1.65.0 / v1.77.0 and `-b r10` |
| 407748 | collect for the three ONT sets | COMPLETED 2026-10-15 19:38-19:39; 40,000 truth rows, 0 GAF duplicates per arm; rebuilt the tables (9 arm x platform sets) |

The six alignment jobs were held from 2026-10-13 until about 13:05 on 2026-10-15 (lab-cluster cap) and then released.

## Results

The first collect (387349) failed on the d46 Illumina GAF: vg 1.77 `convert` writes a header line (`@HD VN:Z:1.0`) that
`evaluate.py` read as a read name. With the fix (Method step 5) collect was rerun (406576). The ONT sets were aligned
and scored afterwards (407745-407748). Nine (arm, platform) sets are scored: Illumina, HiFi and ONT, each on the
GRCh38-only graph, HPRC v1.1 d9 and HPRC v2.1 d46. In every set each simulated read has exactly one primary record in
the GAF and in the BAM (0 duplicates), and the GAF and BAM MAPQs agree for every read (`results/summary.json` and
`eval/<arm>.<platform>.json`: `gaf_records_seen` = `bam_records_seen` = reads, `gaf_duplicates` 0, `mapq_gaf_eq_bam`
1.0).

The block below is regenerated by `summarize.py` (`collect.sh`): per platform the overall table, error-free reads
(Illumina), MAPQ calibration and same reads on two arms, then the pre-alignment checks. The per-stratum tables are in
`results/tables.md`; the main rows are copied after the block.

<!-- RESULTS:BEGIN -->
Generated by summarize.py from the evaluate.py outputs (per-stratum tables: `results/tables.md`).


correct = surjected primary alignment on the true chromosome with start or end within 50 bp of the true GRCh38 interval; wrong = placed on GRCh38 but not correct; perfect = graph alignment (GAF) covers the whole read with no mismatch / indel; MAPQ = giraffe's graph MAPQ.

### illumina: overall

| arm | reads | unmapped | graph-mapped, not surjected | correct (50 bp) | wrong | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | wrong at MAPQ 60 | MAPQ 0 | uses off-ref node |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GRCh38-only graph (H 17) | 5,000,000 | 3,404 (0.068%) | 0 | 98.698% | 61,705 | 1.235% | 63.69% | 33,847 (1.063% of perfect) | 4,538 (0.093%) | 2,408 (0.050%) | 0.84% | 0 (0.00%) |
| HPRC v1.1 d9 | 5,000,000 | 1,500 (0.030%) | 17,436 | 98.642% | 48,959 | 0.983% | 72.94% | 33,204 (0.910% of perfect) | 904 (0.019%) | 344 (0.007%) | 1.07% | 619,803 (12.40%) |
| HPRC v2.1 d46 | 5,000,000 | 1,553 (0.031%) | 19,340 | 98.659% | 46,144 | 0.927% | 73.06% | 32,025 (0.877% of perfect) | 631 (0.013%) | 294 (0.006%) | 1.07% | 622,386 (12.45%) |

Error-free reads only (no simulated sequencing error; any graph mismatch is then a real difference between the read and the graph path it took):

| arm | reads | perfect | perfect but wrong | wrong / placed | wrong at MAPQ>=30 |
|---|---|---|---|---|---|
| GRCh38-only graph (H 17) | 3,702,448 | 86.00% | 33,682 (1.058%) | 1.224% | 3,348 |
| HPRC v1.1 d9 | 3,702,448 | 98.48% | 32,983 (0.905%) | 0.967% | 659 |
| HPRC v2.1 d46 | 3,702,448 | 98.64% | 31,800 (0.871%) | 0.917% | 464 |

### illumina: MAPQ calibration (reads placed on GRCh38; error rate = wrong / placed in the bin; expected = mean 10^(-MAPQ/10))

| MAPQ | GRCh38-only graph (H 17): placed / wrong / rate / expected | HPRC v1.1 d9: placed / wrong / rate / expected | HPRC v2.1 d46: placed / wrong / rate / expected |
|---|---|---|---|
| 0 | 41,978 / 29,500 / 7.03e-01 / 1.00e+00 | 45,494 / 27,828 / 6.12e-01 / 1.00e+00 | 44,253 / 26,433 / 5.97e-01 / 1.00e+00 |
| 1-9 | 43,093 / 21,024 / 4.88e-01 / 7.71e-01 | 49,532 / 18,353 / 3.71e-01 / 7.78e-01 | 46,648 / 17,938 / 3.85e-01 / 7.78e-01 |
| 10-19 | 28,286 / 5,335 / 1.89e-01 / 3.96e-02 | 27,169 / 1,511 / 5.56e-02 / 3.99e-02 | 27,301 / 908 / 3.33e-02 / 4.01e-02 |
| 20-29 | 10,165 / 1,308 / 1.29e-01 / 3.09e-03 | 10,156 / 363 / 3.57e-02 / 2.65e-03 | 10,220 / 234 / 2.29e-02 / 2.63e-03 |
| 30-39 | 17,540 / 1,030 / 5.87e-02 / 8.21e-04 | 16,416 / 254 / 1.55e-02 / 8.38e-04 | 16,238 / 149 / 9.18e-03 / 8.38e-04 |
| 40-49 | 16,423 / 794 / 4.83e-02 / 3.80e-05 | 16,206 / 222 / 1.37e-02 / 3.90e-05 | 16,175 / 144 / 8.90e-03 / 3.90e-05 |
| 50-59 | 7,390 / 306 / 4.14e-02 / 4.00e-06 | 7,109 / 84 / 1.18e-02 / 4.00e-06 | 7,312 / 44 / 6.02e-03 / 4.00e-06 |
| 60 | 4,831,721 / 2,408 / 4.98e-04 / 1.00e-06 | 4,808,982 / 344 / 7.20e-05 / 1.00e-06 | 4,810,960 / 294 / 6.10e-05 / 1.00e-06 |

### illumina: same reads, two arms (correct = within 50 bp; confident wrong = placed wrong with MAPQ>=30)

| A vs B | correct in both | correct only in A | correct only in B | wrong/unplaced in both | confident wrong only in A | confident wrong only in B |
|---|---|---|---|---|---|---|
| linear vs d9 | 4,915,989 | 18,902 | 16,116 | 48,993 | 3,906 | 272 |
| linear vs d46 | 4,915,522 | 19,369 | 17,441 | 47,668 | 4,093 | 186 |
| d9 vs d46 | 4,918,494 | 13,611 | 14,469 | 53,426 | 490 | 217 |

### hifi: overall

| arm | reads | unmapped | graph-mapped, not surjected | correct (50 bp) | wrong | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | wrong at MAPQ 60 | MAPQ 0 | uses off-ref node |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GRCh38-only graph (H 17) | 40,000 | 20 (0.050%) | 0 | 99.455% | 198 | 0.495% | 0.00% | 0 (n/a of perfect) | 13 (0.033%) | 1 (0.003%) | 0.01% | 0 (0.00%) |
| HPRC v1.1 d9 | 40,000 | 8 (0.020%) | 40 | 99.552% | 131 | 0.328% | 0.00% | 0 (n/a of perfect) | 6 (0.015%) | 2 (0.005%) | 0.04% | 36,664 (91.66%) |
| HPRC v2.1 d46 | 40,000 | 0 (0.000%) | 51 | 99.490% | 153 | 0.383% | 0.00% | 0 (n/a of perfect) | 3 (0.008%) | 0 (0.000%) | 0.00% | 36,672 (91.68%) |

### hifi: MAPQ calibration (reads placed on GRCh38; error rate = wrong / placed in the bin; expected = mean 10^(-MAPQ/10))

| MAPQ | GRCh38-only graph (H 17): placed / wrong / rate / expected | HPRC v1.1 d9: placed / wrong / rate / expected | HPRC v2.1 d46: placed / wrong / rate / expected |
|---|---|---|---|
| 0 | 3 / 3 / 1.00e+00 / 1.00e+00 | 7 / 3 / 4.29e-01 / 1.00e+00 | 0 / 0 / n/a / n/a |
| 1-9 | 612 / 167 / 2.73e-01 / 4.47e-01 | 707 / 110 / 1.56e-01 / 4.46e-01 | 775 / 146 / 1.88e-01 / 4.10e-01 |
| 10-19 | 35 / 9 / 2.57e-01 / 5.32e-02 | 46 / 10 / 2.17e-01 / 5.22e-02 | 172 / 4 / 2.33e-02 / 5.24e-02 |
| 20-29 | 50 / 6 / 1.20e-01 / 3.74e-03 | 45 / 2 / 4.44e-02 / 4.07e-03 | 243 / 0 / 0.00e+00 / 3.54e-03 |
| 30-39 | 145 / 5 / 3.45e-02 / 3.55e-04 | 125 / 2 / 1.60e-02 / 3.40e-04 | 719 / 0 / 0.00e+00 / 3.49e-04 |
| 40-49 | 474 / 7 / 1.48e-02 / 3.40e-05 | 409 / 1 / 2.45e-03 / 3.30e-05 | 1,971 / 1 / 5.07e-04 / 3.70e-05 |
| 50-59 | 1,515 / 0 / 0.00e+00 / 4.00e-06 | 1,387 / 1 / 7.21e-04 / 4.00e-06 | 4,766 / 2 / 4.20e-04 / 4.00e-06 |
| 60 | 37,146 / 1 / 2.70e-05 / 1.00e-06 | 37,226 / 2 / 5.40e-05 / 1.00e-06 | 31,303 / 0 / 0.00e+00 / 1.00e-06 |

### hifi: same reads, two arms (correct = within 50 bp; confident wrong = placed wrong with MAPQ>=30)

| A vs B | correct in both | correct only in A | correct only in B | wrong/unplaced in both | confident wrong only in A | confident wrong only in B |
|---|---|---|---|---|---|---|
| linear vs d9 | 39,721 | 61 | 100 | 118 | 11 | 4 |
| linear vs d46 | 39,710 | 72 | 86 | 132 | 12 | 2 |
| d9 vs d46 | 39,736 | 85 | 60 | 119 | 4 | 1 |

### ont: overall

| arm | reads | unmapped | graph-mapped, not surjected | correct (50 bp) | wrong | wrong / placed | perfect | perfect but wrong | wrong at MAPQ>=30 | wrong at MAPQ 60 | MAPQ 0 | uses off-ref node |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GRCh38-only graph (H 17) | 40,000 | 282 (0.705%) | 0 | 98.528% | 307 | 0.773% | 0.85% | 1 (0.295% of perfect) | 170 (0.433%) | 158 (0.405%) | 0.00% | 0 (0.00%) |
| HPRC v1.1 d9 | 40,000 | 234 (0.585%) | 46 | 98.690% | 244 | 0.614% | 7.94% | 1 (0.032% of perfect) | 131 (0.335%) | 125 (0.321%) | 0.07% | 37,387 (93.47%) |
| HPRC v2.1 d46 | 40,000 | 218 (0.545%) | 66 | 97.817% | 589 | 1.483% | 8.62% | 2 (0.058% of perfect) | 299 (0.764%) | 181 (0.466%) | 0.00% | 37,349 (93.37%) |

### ont: MAPQ calibration (reads placed on GRCh38; error rate = wrong / placed in the bin; expected = mean 10^(-MAPQ/10))

| MAPQ | GRCh38-only graph (H 17): placed / wrong / rate / expected | HPRC v1.1 d9: placed / wrong / rate / expected | HPRC v2.1 d46: placed / wrong / rate / expected |
|---|---|---|---|
| 0 | 1 / 1 / 1.00e+00 / 1.00e+00 | 9 / 7 / 7.78e-01 / 1.00e+00 | 0 / 0 / n/a / n/a |
| 1-9 | 259 / 112 / 4.32e-01 / 4.12e-01 | 336 / 96 / 2.86e-01 / 4.22e-01 | 329 / 173 / 5.26e-01 / 4.07e-01 |
| 10-19 | 110 / 17 / 1.55e-01 / 4.08e-02 | 86 / 7 / 8.14e-02 / 4.08e-02 | 133 / 63 / 4.74e-01 / 4.21e-02 |
| 20-29 | 127 / 7 / 5.51e-02 / 4.15e-03 | 173 / 3 / 1.73e-02 / 4.00e-03 | 128 / 54 / 4.22e-01 / 4.54e-03 |
| 30-39 | 112 / 2 / 1.79e-02 / 4.85e-04 | 90 / 1 / 1.11e-02 / 4.70e-04 | 88 / 46 / 5.23e-01 / 4.65e-04 |
| 40-49 | 64 / 4 / 6.25e-02 / 5.00e-05 | 53 / 2 / 3.77e-02 / 4.30e-05 | 88 / 40 / 4.55e-01 / 4.20e-05 |
| 50-59 | 46 / 6 / 1.30e-01 / 4.00e-06 | 43 / 3 / 6.98e-02 / 4.00e-06 | 67 / 32 / 4.78e-01 / 5.00e-06 |
| 60 | 38,999 / 158 / 4.05e-03 / 1.00e-06 | 38,930 / 125 / 3.21e-03 / 1.00e-06 | 38,883 / 181 / 4.66e-03 / 1.00e-06 |

### ont: same reads, two arms (correct = within 50 bp; confident wrong = placed wrong with MAPQ>=30)

| A vs B | correct in both | correct only in A | correct only in B | wrong/unplaced in both | confident wrong only in A | confident wrong only in B |
|---|---|---|---|---|---|---|
| linear vs d9 | 39,341 | 70 | 135 | 454 | 63 | 24 |
| linear vs d46 | 38,997 | 414 | 130 | 459 | 78 | 207 |
| d9 vs d46 | 39,058 | 418 | 69 | 455 | 29 | 197 |

### Where the simulation samples (one-to-one segments, GRCh38 non-N bases covered)

| chrom | non-N bases | hap1 covered | hap2 covered | hap1 easy / difficult / extreme covered |
|---|---|---|---|---|
| chr1 | 230,481,014 | 97.7% | 97.5% | 99.9% / 99.6% / 81.0% |
| chr20 | 63,944,257 | 94.7% | 94.4% | 99.9% / 98.8% / 56.9% |

### truth check

```
{'checked': 19271, 'exact': 19267, 'mismatch': 4, 'skipped_with_indel': 729}
('i4835/1', 'chr1', 101930920, 2, 1, 2)
('i8358/1', 'chr1', 108947738, 1, 1, 1)
('i11178/1', 'chr1', 72909987, 2, 1, 2)
('i15106/1', 'chr20', 22608854, 2, 2, 1)
```

### ont error profile

```
ONT-like reads vs their own templates (minimap2 -c, primary blocks): n=1801, NM/block median 0.0217, q25 0.0082, q75 0.0496, q90 0.0860
```

<!-- RESULTS:END -->

hap2 SMaHT coverage (`results/segment_stats.json`; the generated table shows hap1 only):

| chrom | hap2 easy / difficult / extreme covered | SMaHT non-N bases easy / difficult / extreme |
|---|---|---|
| chr1 | 99.85% / 99.63% / 78.86% | 179,223,254 / 25,309,436 / 25,948,324 |
| chr20 | 99.69% / 98.54% / 55.01% | 49,278,286 / 7,258,874 / 7,407,097 |

### Selected strata (copied from `results/tables.md`)

Each cell lists GRCh38-only graph / HPRC v1.1 d9 / HPRC v2.1 d46. Percentages as in `tables.md`: correct and unmapped
of all reads in the stratum, wrong of the reads placed on GRCh38.

**Illumina** (5,000,000 reads):

| stratum | reads | correct (50 bp) % | wrong / placed % | wrong at MAPQ>=30 (reads) | perfect % | perfect but wrong (reads) | unmapped % |
|---|---:|---|---|---|---|---|---|
| HG008-N variants in read: 0 | 4,279,026 | 99.177 / 99.022 / 99.028 | 0.799 / 0.827 / 0.804 | 466 / 98 / 68 | 73.97 / 74.02 / 74.02 | 23,739 / 25,367 / 24,839 | 0.025 / 0.027 / 0.026 |
| variants in read: 1 | 600,460 | 96.212 / 96.475 / 96.561 | 3.430 / 1.622 / 1.388 | 2,669 / 389 / 236 | 2.84 / 67.04 / 67.78 | 7,974 / 5,698 / 5,200 | 0.370 / 0.052 / 0.065 |
| variants in read: 2-3 | 112,202 | 95.400 / 96.930 / 97.002 | 4.516 / 2.514 / 2.226 | 849 / 203 / 129 | 1.42 / 64.44 / 65.53 | 1,574 / 1,525 / 1,403 | 0.088 / 0.038 / 0.039 |
| variants in read: >=4 | 8,312 | 76.239 / 82.736 / 82.868 | 23.595 / 15.713 / 14.371 | 554 / 214 / 198 | 6.80 / 59.44 / 60.64 | 560 / 614 / 583 | 0.217 / 0.072 / 0.120 |
| read spans an SV (>=50 bp) of its haplotype | 18,851 | 25.139 / 21.240 / 21.686 | 71.670 / 51.015 / 46.047 | 2,350 / 316 / 239 | 23.50 / 62.74 / 64.29 | 4,173 / 2,431 / 2,203 | 11.262 / 0.833 / 1.358 |
| >=10 read bases absent from GRCh38 | 28,075 | 47.316 / 45.817 / 46.169 | 48.691 / 25.354 / 21.916 | 2,512 / 338 / 247 | 17.37 / 59.73 / 60.99 | 4,256 / 2,520 / 2,271 | 7.783 / 0.609 / 0.994 |
| neither SV nor insertion | 4,970,368 | 98.992 / 98.942 / 98.958 | 0.984 / 0.897 / 0.855 | 1,968 / 548 / 364 | 63.97 / 73.02 / 73.13 | 29,551 / 30,631 / 29,696 | 0.024 / 0.027 / 0.026 |
| nearest variant: inside read | 720,974 | 95.856 / 96.388 / 96.472 | 3.832 / 1.925 / 1.668 | 4,072 / 806 / 563 | 2.67 / 66.54 / 67.34 | 10,108 / 7,837 / 7,186 | 0.325 / 0.050 / 0.062 |
| nearest variant: 1-50 bp | 366,880 | 98.727 / 98.688 / 98.710 | 1.215 / 1.075 / 0.998 | 220 / 44 / 30 | 73.62 / 73.98 / 73.99 | 2,803 / 2,745 / 2,585 | 0.059 / 0.044 / 0.045 |
| nearest variant: 51-500 bp | 1,709,231 | 99.136 / 99.046 / 99.060 | 0.840 / 0.799 / 0.759 | 241 / 50 / 30 | 74.01 / 74.04 / 74.05 | 9,857 / 9,779 / 9,386 | 0.025 / 0.028 / 0.027 |
| nearest variant: >500 bp | 2,202,915 | 99.284 / 99.059 / 99.056 | 0.697 / 0.807 / 0.807 | 5 / 4 / 8 | 74.00 / 74.01 / 74.01 | 11,079 / 12,843 / 12,868 | 0.019 / 0.023 / 0.022 |
| SMaHT easy | 4,021,943 | 99.994 / 99.982 / 99.979 | 0.002 / 0.003 / 0.002 | 25 / 19 / 12 | 64.90 / 73.22 / 73.30 | 17 / 71 / 43 | 0.004 / 0.001 / 0.001 |
| SMaHT difficult | 577,525 | 99.423 / 99.405 / 99.395 | 0.339 / 0.081 / 0.083 | 341 / 51 / 57 | 58.53 / 72.24 / 72.36 | 326 / 215 / 228 | 0.239 / 0.015 / 0.037 |
| SMaHT extreme | 400,532 | 84.639 / 84.086 / 84.349 | 14.967 / 12.561 / 11.893 | 4,172 / 834 / 562 | 58.92 / 71.17 / 71.59 | 33,504 / 32,918 / 31,754 | 0.463 / 0.346 / 0.329 |
| STR fraction of read >= 0.5 | 5,219 | 89.251 / 88.139 / 88.159 | 10.199 / 2.315 / 2.480 | 169 / 20 / 33 | 17.07 / 53.38 / 54.90 | 13 / 28 / 37 | 0.613 / 0.326 / 0.709 |
| homopolymer >= 7 in read | 838,028 | 99.084 / 99.079 / 99.105 | 0.852 / 0.691 / 0.647 | 475 / 99 / 52 | 58.48 / 72.28 / 72.45 | 3,507 / 3,871 / 3,678 | 0.065 / 0.022 / 0.026 |

Illumina reads whose graph alignment uses an off-reference node (d9 619,803, d46 622,386 reads; `tables.md`): correct
95.499 / 95.474 %, wrong / placed 1.736 / 1.464 %, wrong at MAPQ >= 30 290 / 229 reads (0.050 / 0.040 %). The other
reads: correct 99.121 / 99.147 %, wrong / placed 0.879 / 0.853 %. Every "graph-mapped, not surjected" read is in the
off-reference group (17,436 of the d9 group, 19,340 of the d46 group; `summary.json`).

Placement counts behind the SV row (`summary.json`, `strata.sv_or_insertion`): of the 18,851 SV-spanning reads,

| arm | placed on GRCh38 | correct | wrong | graph-mapped, not surjected | unmapped |
|---|---:|---:|---:|---:|---:|
| GRCh38-only graph | 16,728 | 4,739 | 11,989 | 0 | 2,123 |
| HPRC v1.1 d9 | 8,174 | 4,004 | 4,170 | 10,520 | 157 |
| HPRC v2.1 d46 | 7,577 | 4,088 | 3,489 | 11,018 | 256 |

For the >500 bp stratum the wrong reads are 15,360 / 17,755 / 17,752 (`summary.json`).

**HiFi** (40,000 reads; no read is perfect on any arm, see Interpretation 6):

| stratum | reads | correct (50 bp) % | wrong / placed % | wrong at MAPQ>=30 (reads) | unmapped % |
|---|---:|---|---|---|---|
| HG008-N variants in read: 0 | 1,633 | 99.939 / 99.755 / 99.755 | 0.061 / 0.184 / 0.245 | 0 / 0 / 0 | 0.000 / 0.000 / 0.000 |
| variants in read: 1 | 2,478 | 98.588 / 98.426 / 98.547 | 1.253 / 0.813 / 0.772 | 4 / 0 / 0 | 0.161 / 0.000 / 0.000 |
| variants in read: >=4 | 31,790 | 99.487 / 99.616 / 99.544 | 0.472 / 0.296 / 0.359 | 9 / 5 / 3 | 0.041 / 0.025 / 0.000 |
| read spans an SV (>=50 bp) of its haplotype | 2,892 | 96.853 / 97.545 / 97.545 | 2.608 / 1.467 / 1.605 | 10 / 6 / 3 | 0.553 / 0.138 / 0.000 |
| >=10 read bases absent from GRCh38 | 10,466 | 98.643 / 99.025 / 99.025 | 1.215 / 0.614 / 0.614 | 11 / 6 / 3 | 0.143 / 0.029 / 0.000 |
| neither SV nor insertion | 28,963 | 99.748 / 99.745 / 99.658 | 0.238 / 0.225 / 0.297 | 0 / 0 / 0 | 0.014 / 0.014 / 0.000 |
| nearest variant: >500 bp | 1,424 | 100.000 / 99.719 / 99.789 | 0.000 / 0.211 / 0.211 | 0 / 0 / 0 | 0.000 / 0.000 / 0.000 |
| SMaHT easy | 37,644 | 99.997 / 99.997 / 100.000 | 0.000 / 0.003 / 0.000 | 0 / 1 / 0 | 0.003 / 0.000 / 0.000 |
| SMaHT extreme | 1,658 | 86.972 / 89.264 / 87.696 | 12.073 / 8.075 / 9.521 | 13 / 5 / 3 | 1.086 / 0.483 / 0.000 |
| STR fraction of read >= 5 % | 106 | 84.906 / 93.396 / 95.283 | 9.091 / 3.884 / 0.000 | 2 / 2 / 0 | 6.604 / 0.000 / 0.000 |

HiFi reads using an off-reference node: 36,664 (91.66 %, d9) / 36,672 (91.68 %, d46).

**ONT** (40,000 reads; `tables.md` "## ont" sections, 2026-10-15 19:39):

| stratum | reads | correct (50 bp) % | wrong / placed % | wrong at MAPQ>=30 (reads) | unmapped % |
|---|---:|---|---|---|---|
| HG008-N variants in read: 0 | 1,754 | 99.088 / 98.860 / 98.233 | 0.230 / 0.459 / 1.204 | 4 / 5 / 15 | 0.684 / 0.570 / 0.570 |
| variants in read: 1 | 2,514 | 97.733 / 97.772 / 96.619 | 1.562 / 0.807 / 2.135 | 25 / 8 / 23 | 0.716 / 0.517 / 0.477 |
| variants in read: >=4 | 31,635 | 98.536 / 98.736 / 97.895 | 0.761 / 0.605 / 1.429 | 130 / 109 / 230 | 0.708 / 0.597 / 0.553 |
| read spans an SV (>=50 bp) of its haplotype | 2,874 | 95.442 / 96.277 / 96.033 | 3.517 / 2.053 / 2.439 | 66 / 35 / 28 | 1.079 / 0.731 / 0.522 |
| >=10 read bases absent from GRCh38 | 10,467 | 97.707 / 98.175 / 97.487 | 1.607 / 0.935 / 1.658 | 94 / 51 / 78 | 0.697 / 0.525 / 0.459 |
| neither SV nor insertion | 28,972 | 98.830 / 98.882 / 97.936 | 0.459 / 0.486 / 1.414 | 71 / 73 / 213 | 0.715 / 0.611 / 0.580 |
| SMaHT easy | 37,633 | 99.086 / 99.184 / 99.208 | 0.235 / 0.248 / 0.251 | 88 / 93 / 94 | 0.680 / 0.569 / 0.542 |
| SMaHT difficult | 678 | 98.230 / 98.525 / 97.935 | 0.448 / 0.447 / 0.747 | 3 / 3 / 5 | 1.327 / 1.032 / 1.032 |
| SMaHT extreme | 1,689 | 86.205 / 87.744 / 66.785 | 12.919 / 9.080 / 30.284 | 79 / 35 / 200 | 1.006 / 0.770 / 0.414 |
| STR fraction of read >= 5 % | 121 | 84.297 / 91.736 / 92.562 | 11.304 / 4.310 / 2.609 | 8 / 3 / 0 | 4.959 / 1.653 / 0.826 |

Wrong reads by SMaHT class (`eval/<arm>.ont.json`, `strata.smaht_region`; GRCh38-only / d9 / d46): easy 88 / 93 / 94,
difficult 3 / 3 / 5, extreme 216 / 148 / 490 (of 1,672 / 1,630 / 1,618 placed). Of these, the reads whose placement does
not overlap the true interval at all (`wrong_not_overlapping`) are 1 / 0 / 0 in easy and 181 / 124 / 464 in extreme
regions (overall 182 / 124 / 466). The easy-region "wrong" reads therefore overlap their true locus on every arm: start
and end are both more than 50 bp off (presumably clipped noisy read ends; not inspected), not misplacements. ONT reads
using an off-reference node: 37,387 (93.47 %, d9) / 37,349 (93.37 %, d46). Graph-mapped, not surjected: 0 / 46 / 66 (64
of the d46 ones in SMaHT extreme).

## Interpretation

1. **Illumina: the HPRC graphs misplace fewer reads, above all fewer confident ones (R3.6).**
   - Wrong / placed: 1.235 % (GRCh38-only graph) -> 0.983 % (d9) / 0.927 % (d46).
   - Wrong at MAPQ >= 30: 4,538 -> 904 / 631 reads; at MAPQ 60: 2,408 -> 344 / 294.
   - On the same reads, 3,906 are confidently wrong only on the GRCh38-only graph against 272 only on d9 (4,093
     against 186 for d46).
   - Correct within 50 bp is flat: 98.698 / 98.642 / 98.659 %. The HPRC arms are 0.056 / 0.039 points lower because
     17,436 / 19,340 reads align in the graph without a GRCh38 projection (counted neither correct nor wrong). On the same
     reads, 18,902 are correct only on the GRCh38-only graph and 16,116 only on d9 (19,369 / 17,441 for d46).
2. **More perfect alignments come without more perfect-but-wrong ones.**
   - Perfect: 63.69 -> 72.94 / 73.06 % of all reads; 86.00 -> 98.48 / 98.64 % of the 3,702,448 reads without a simulated
     error.
   - Perfect but wrong: 33,847 (1.063 % of perfect) -> 33,204 (0.910 %) / 32,025 (0.877 %). The count is lower and so is
     the share.
   - Almost all perfect-but-wrong reads come from SMaHT extreme regions (33,504 / 32,918 / 31,754). Each is an exact match
     to another copy of a repeat, and the GRCh38-only graph has the same problem. In SMaHT easy regions there are
     17 / 71 / 43 of them among 4,021,943 reads.
3. **Where the HPRC graphs help.** The gain grows with the number of HG008-N variants in the read:
   - wrong / placed: 3.430 -> 1.622 / 1.388 % with one variant, 4.516 -> 2.514 / 2.226 % with 2-3, 23.595 -> 15.713 /
     14.371 % with four or more (correct 76.239 -> 82.736 / 82.868 %);
   - SMaHT difficult: wrong / placed 0.339 -> 0.081 / 0.083 %, confident wrong 341 -> 51 / 57;
   - SMaHT extreme: 14.967 -> 12.561 / 11.893 %, confident wrong 4,172 -> 834 / 562;
   - STR-rich reads (STR fraction >= 0.5): 10.199 -> 2.315 / 2.480 %.
4. **Where the HPRC graphs do worse.**
   - **SV-spanning Illumina reads** (18,851): correct 25.139 % on the GRCh38-only graph, 21.240 % on d9 and 21.686 % on
     d46. The HPRC graphs align almost all of them (unmapped 11.262 % -> 0.833 / 1.358 %), but 10,520 / 11,018 align over
     SV-allele nodes and have no GRCh38 projection. Of the reads that are placed, 4,004 of 8,174 (d9) and 4,088 of 7,577
     (d46) are correct, against 4,739 of 16,728 on the GRCh38-only graph, and confident wrong falls from 2,350 to
     316 / 239. The loss is therefore in projecting non-GRCh38 alignments back to GRCh38, not in misplacement. Reads with
     >= 10 inserted bases behave the same way (correct 47.316 -> 45.817 / 46.169 %).
   - **Reads far from any HG008-N variant** (> 500 bp; 2,202,915 reads): slightly more misplacement on the HPRC graphs.
     Wrong / placed rises 0.697 -> 0.807 / 0.807 % (15,360 -> 17,755 / 17,752 reads), correct falls 99.284 -> 99.059 /
     99.056 %, and perfect but wrong rises 11,079 -> 12,843 / 12,868. These extra errors have low MAPQ: wrong at
     MAPQ >= 30 is 5 / 4 / 8 reads, and MAPQ 0 rises 0.55 -> 0.80 / 0.79 % (`tables.md`). Reads with no HG008-N variant at
     all show the same pattern: wrong / placed 0.799 / 0.827 / 0.804 %, correct 99.177 / 99.022 / 99.028 %. A plausible
     reading (not tested here) is that the population haplotypes add near-identical paths in repeats, so that a read
     that matches GRCh38 exactly now has more equally good placements.
   - SMaHT extreme: correct 84.639 -> 84.086 / 84.349 % even though wrong / placed falls; 13,978 (d9) and 15,765 (d46)
     extreme-region reads are not surjected (`summary.json`), and MAPQ 0 rises 10.15 -> 12.81 / 12.70 %.
   - HiFi reads with neither an SV nor an insertion (28,963): d46 is slightly worse (wrong / placed 0.238 / 0.225 /
     0.297 %, correct 99.748 / 99.745 / 99.658 %), but none of these wrong reads has MAPQ >= 30.
   - **ONT reads on d46** (Interpretation 7): about twice the misplacement of the GRCh38-only graph (1.483 % against
     0.773 % of placed reads), almost all in SMaHT extreme regions, with many of the errors at MAPQ >= 30.
5. **MAPQ is overconfident on every arm, least on the HPRC graphs** (Illumina, generated calibration table).
   - Observed error at MAPQ 60: 4.98e-04 / 7.20e-05 / 6.10e-05 (nominal 1.00e-06).
   - MAPQ 30-39: 5.87e-02 / 1.55e-02 / 9.18e-03 (nominal 8.21e-04 / 8.38e-04 / 8.38e-04).
   - MAPQ 10-19: 1.89e-01 / 5.56e-02 / 3.33e-02 (nominal 3.96e-02 / 3.99e-02 / 4.01e-02), so the HPRC arms are close to
     nominal in this bin.
   - Below MAPQ 10 the error is lower than nominal (MAPQ 0: 7.03e-01 / 6.12e-01 / 5.97e-01).
   - So no reference gives a calibrated MAPQ, but a MAPQ threshold removes a larger share of the wrong reads on the HPRC
     graphs. MAPQ 0 reads are 0.84 / 1.07 / 1.07 % of all reads.
   - HiFi: the wrong reads sit at MAPQ 1-9 (167 / 110 / 146 of 198 / 131 / 153). At MAPQ 60, 1 of 37,146, 2 of 37,226 and
     0 of 31,303 reads are wrong. d46 gives MAPQ 60 to fewer reads and MAPQ 30-59 to more (7,456 against 2,134 / 1,921
     reads, 3 of them wrong): the same d46 MAPQ shift as on the real reads (E11a).
6. **HiFi** (40,000 reads): correct 99.455 / 99.552 / 99.490 %, wrong / placed 0.495 / 0.328 / 0.383 %, wrong at
   MAPQ >= 30 13 / 6 / 3 reads, unmapped 20 / 8 / 0.
   - No HiFi read is perfect on any arm: 0.1 % simulated error over reads of about 15 kb leaves essentially no error-free
     read, so "perfect but wrong" does not apply to long reads (as in E11a, where per-kb rates replace it).
   - SV-spanning HiFi reads are placed correctly 96.853 / 97.545 / 97.545 % of the time, unlike Illumina: the reads are
     long enough to anchor in the flanks.
   - SMaHT extreme: 86.972 / 89.264 / 87.696 % correct.
7. **ONT** (40,000 noisy reads, mean accuracy 95.9 %): d9 is slightly better than the GRCh38-only graph, **d46 is
   worse** (R1.3's noisy-ONT question).
   - Wrong / placed: 0.773 % (GRCh38-only) / 0.614 % (d9) / 1.483 % (d46); wrong at MAPQ >= 30: 170 (0.433 %) / 131
     (0.335 %) / 299 (0.764 %) of the MAPQ >= 30 placements; at MAPQ 60: 158 / 125 / 181. Correct within 50 bp:
     98.528 / 98.690 / 97.817 %. Unmapped falls slightly: 0.705 / 0.585 / 0.545 %.
   - The d46 excess sits in SMaHT extreme regions: wrong 216 / 148 / 490 reads (12.919 / 9.080 / 30.284 % of placed),
     correct 86.205 / 87.744 / 66.785 %, confident wrong 79 / 35 / 200. That is 274 of the 282 extra wrong reads and 121
     of the 129 extra confident ones on d46 against the GRCh38-only graph (my subtraction of the per-class counts above).
     SMaHT easy regions are equal on the three arms (0.235 / 0.248 / 0.251 %; these 88-94 reads overlap their true
     interval, see the ONT strata note).
   - d46 MAPQ is badly calibrated for ONT: in every bin from MAPQ 10 to 59 about half of the d46 placements are wrong
     (63 / 133, 54 / 128, 46 / 88, 40 / 88, 32 / 67; rates 0.42-0.52), against 0.018-0.155 (GRCh38-only) and
     0.011-0.081 (d9). HiFi on d46 shows the opposite (MAPQ 10-59 nearly error-free, Interpretation 5), so the effect is
     specific to the noisy reads (or to the `-b r10` preset of vg 1.77), not to the d46 MAPQ scale in general.
   - Same reads: 418 are correct only on d9 and 69 only on d46; 197 are confidently wrong only on d46 and 29 only on d9
     (GRCh38-only vs d46: 78 against 207).
   - By HG008-N variants in the read (wrong / placed, GRCh38-only -> d9 / d46): one variant 1.562 -> 0.807 / 2.135 %,
     >= 4 variants 0.761 -> 0.605 / 1.429 %, none 0.230 -> 0.459 / 1.204 % (1,754 reads). So d9 helps where the read
     carries a variant, as for Illumina, and d46 is worse in every variant-count class (2-3 variants: 0.614 -> 0.638 / 1.621 %). STR-rich reads gain on both graphs (11.304 ->
     4.310 / 2.609 %, 121 reads), and SV-spanning reads are placed better on both (wrong / placed 3.517 -> 2.053 /
     2.439 %, correct 95.442 -> 96.277 / 96.033 %).
   - This matches E11a's real-data result (E11a README, Interpretation 6): real ONT reads have more confidently
     discordant placements on d46 than on d9 (MAPQ >= 30 in both arms, vs the GRCh38-only graph: 1,953 against 1,019 for
     HG008T, 16,164 against 4,356 for COLO829T). E11b now shows that, at least in simulation, these d46 disagreements
     are mostly misplacements.
   - Two qualifications: d46 is aligned with vg 1.77 and the other arms with 1.66 / 1.65, so the aligner version is
     confounded with the graph; and the pbsim3 reads are noisier than R10.4.1 sup reads (Caveats), so absolute rates are
     pessimistic.
8. **For the response (R3.6 point 3, R1.3 point 1).** With known origin, the HPRC graphs give more perfect alignments
   and fewer misplaced short and HiFi reads; confident Illumina misplacement drops 5x (d9) to 7x (d46). The costs are:
   - 0.35 % (d9) / 0.39 % (d46) of Illumina reads, most of them on SV or inserted sequence, align without a GRCh38
     projection;
   - slightly more low-MAPQ misplacement of reads that carry no variant;
   - for ONT, d46 (aligned with vg 1.77) misplaces about twice as many reads as the GRCh38-only graph (1.483 % against
     0.773 %; at MAPQ >= 30 0.764 % against 0.433 %), almost all in SMaHT extreme regions, while d9 is slightly better
     (0.614 %; 0.335 %). This should be stated plainly in R3.6 and R1.3.

## Caveats

- **ONT reads are a pessimistic case.** pbsim3 ERRHMM-ONT-HQ at `--accuracy-mean 0.98` gave a mean accuracy of 95.9 %
  and a median NM / block of 2.2 % with a long tail (90th percentile 8.6 %). This is noisier than R10.4.1 sup reads, so
  absolute ONT misplacement rates are likely higher than for real data; the comparison between arms uses the same reads.
  Constant Q20 qualities. Per-stratum ONT counts are small (1,689 SMaHT extreme reads, 121 STR-rich reads). In SMaHT
  easy regions about 90 ONT reads per arm count as "wrong" although they overlap their true interval (both ends more than
  50 bp off), so the ONT "wrong" rate includes a constant end-point component.
- **The d9 arm used vg 1.65.0, not the production 1.66.0.** The HG008T production d9 GAMs were made with vg 1.66.0
  (2026-10-15 correction, `../../points/R3.1.md`). The COLO829T d9 GAMs are 1.64.1-1.66.0 (Fiber-seq), 1.64.0-1.66.0 (ONT)
  and <= 1.62 (Illumina) by their `PARAMS_JSON` fingerprint (E22 check 2). The other arms match production
  (GRCh38-only 1.66.0, d46 1.77.0). Because the three arms also use three giraffe versions, part of any difference
  between them can come from the aligner version rather than the graph: the d46 MAPQ scale, and the d46 ONT
  misplacement (Interpretation 7), in particular.
- **CHM13 is inside both HPRC arms.** Both HPRC GBZs carry CHM13 as a reference-sense path next to GRCh38
  (`vg paths -M -S CHM13`, vg 1.77: 109 subpaths in d9, 93 in d46, all with SENSE = REFERENCE; E8b H3 / H6). The
  frequency filter does not apply to reference paths, so CHM13 sequence, including sequence that GRCh38 lacks, is
  unfiltered in both HPRC graphs and reads can align to it. The HPRC arms therefore differ from the GRCh38-only arm in
  reference completeness as well as in population haplotypes. The CHM13-only arm that would separate the two was dropped
  by the user (2026-10-13), so this experiment cannot attribute the gains to one or the other, and the response to R3.6
  point 3 / R1.3 point 1 must be narrowed to GRCh38-only graph vs HPRC v1.1 / v2.1.
- **HG008-N is in no graph.** HPRC v1.1 has no HG008, and HG008 is not in the v2.1 sample list (checked 2026-10-01).
  The simulated reads are therefore out of sample for every arm, as real patient reads are; HG008-N alleles that are
  common in HPRC are still in the graphs.
- **Simulated reads.** The error models are simple for Illumina and HiFi (uniform; constant qualities). Simulation
  covers chr1 + chr20 only, and only where the assembly-to-GRCh38 truth is one-to-one: unresolved segdups, centromeres
  and regions covered twice are under-sampled by construction (55-81 % of SMaHT extreme bases are covered). The results
  are a lower bound on misplacement in difficult regions, not a genome-wide rate. Reads come from the normal (HG008-N)
  genome; tumor-specific rearrangements are not modelled.
- The truth is the dipcall minimap2 asm5 alignment. Within one-to-one blocks of >= 50 kb this is reliable, but small
  local misalignments (e.g., indel placement in STRs) are inside the 50 bp tolerance only if they are short.
- The metrics are per read, so depth does not matter. The read counts set the precision instead: about 1.3x of the
  segment bases for Illumina (5,000,000 x 150 bp over 570,972,967 bases), and 40,000 reads each for HiFi and ONT. HiFi
  and ONT strata are small (e.g. 106 / 121 STR-rich reads, 1,658 / 1,689 SMaHT extreme reads).
- Surjection to GRCh38 is needed to score graph alignments. A read that aligns in the graph but has no GRCh38 projection
  is counted as "graph-mapped, not surjected", not as wrong. This lowers the HPRC arms' correct share for SV-spanning
  reads (Interpretation 4) without counting a misplacement.
- **No linear-aligner arm.** There is no bwa-mem2 / minimap2 / pbmm2 GRCh38 BAM of the simulated reads. Comparators are
  not re-run, and the GRCh38 arm here is the GRCh38-only vg graph, which isolates the graph effect within giraffe.
- Stratification uses SMaHT classes and sequence-derived STR / homopolymer masks, not GIAB v3.5 (not on the server).

## Status

**DONE for Illumina, HiFi and ONT (2026-10-15 19:39).**

DONE:
- one-to-one segments, the Illumina / HiFi / ONT simulation (387278, 387340, 387341), the truth check, the segment
  coverage per SMaHT class, the ONT error profile and the code smoke test on a tiny graph (2026-10-13);
- Illumina on the three arms (387342-387344) and HiFi on the three arms (387345-387347), 2026-10-15;
- ONT on the three arms (407745-407747, 2026-10-15 19:23-19:38). The first long-read jobs (387345-387347) had been
  submitted with `--export=ALL,ARM=<arm>,PLATFORMS=hifi,ont`; `sbatch` splits `--export` at commas, so `PLATFORMS`
  arrived as `hifi` and ONT was skipped. The ONT jobs were submitted one platform per job;
- scoring and tables: collect 406576 (2026-10-15 17:57) after the `evaluate.py` GAF-header fix, then collect 407748
  (19:39) for the ONT sets; `results/tables.md`, `results/summary.{tsv,json}` and the generated block above hold 9
  (arm, platform) sets.

DROPPED (user decision 2026-10-13): the CHM13-only control. 387281 and 387348 were cancelled. The CHM13 truth files in
`chm13_truth/` (including `ont.chm13_truth.npz`, written by collect 387349) are unused.

Test leftovers in the data directory, not used by any result:
- `test_illumina/` (20,000-pair login-node test);
- `test_hifi/`, `test_ont/` (empty; their login-node tests stopped after loading the haplotypes);
- `tinytest/` (pbsim3 and giraffe smoke tests).
