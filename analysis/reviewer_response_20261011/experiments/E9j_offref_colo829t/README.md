# E9j: COLO829T off-reference somatic SNVs from the patient's assemblies (R3.3 point 1, R1.1)

The COLO829T counterpart of the HG008T analysis behind R3.3 point 1: `pansoma_net_v2_runs/assembly_offref_20260929`
(events in the normal assembly, giraffe, classes A / B / C / D, keys), E1.1 section 5 (model scores on the keys: base,
off-reference rescue, all_H) and E9h (read support of the keys: the A1-read consensus in the tumor assembly and not in
the normal). Possible since 2026-10-23 with a COLO829T tumor assembly (qfu's hifiasm) and the COLO829BL donor-specific
assembly paths the user gave.

**Status (2026-10-23): done** (jobs and logs `tmp/reviewer_cpu_20261023/E9j/`, `jobs.txt`; data
`/scratch/jshen/data/reviewer_experiments_20261013/E9j_offref_colo829t/`; code snapshot `code_8580a9f` there).

## Data

| what | where |
|---|---|
| somatic truth | SMaHT COLO829T, `pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz` (the sets' labels.manifest.json); 42,035 SNVs on chr1-22 (INDELs not done) |
| normal assembly | COLO829BL DSA, SMaHT verkko 2.1 hapX / hapY (`/scratch/qfu/COLO829BL_DSA/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFTH1B2J6-verkko_2.1_hapX.fa`, `...SMASFK4V33V8-verkko_2.1_hapY.fa`) |
| tumor assembly | hifiasm 0.25.0 of the COLO829T Fiber-seq HiFi reads with the COLO829T ONT D003 reads as ultra-long (qfu, `/scratch/qfu/COLO829T_fiberseq/COLO829T_hifiasm/COLO829T.asm.bp.hap{1,2}.p_ctg.gfa`), as FASTA in `E9c_offref_filter_validation/asm_targets/COLO829T_hifiasm_hap{1,2}.fa.gz` (hap1 368 contigs, 2.68 Gb) |
| graphs | d9 HPRC v1.1 (vg 1.65.0, shortread minimizers v10), d46 HPRC v2.1 (vg 1.77.0); GRCh38 paths and graph indexes of the tensor sets |
| models | the in-sample COLO829T SNV models of R3.1 / R3.4: d9 `COLO829_Illumina_SNV_d9_normbq_sweep_e06` (epoch 6, seed 0), `COLO829_{fiberseq,ONT}_SNV_7fc7396_ddp4_b1024`; d46 `COLO829_{Illumina,fiberseq,ONT}_SNV_d46_7fc7396` (`pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/best.pth`) |

## Method

1. `s1_windows.py`. The SMaHT truth is on GRCh38 only, so each truth SNV is placed in the patient frame: its GRCh38
   window (100-bp flanks) is aligned to DSA hapX and hapY (minimap2 2.28 -cx sr, `s1_map.sh MODE=grch38`); the primary
   alignment (>= 90 % aligned, MAPQ >= 20, the site in a match or mismatch block) gives the DSA position and base. Per
   haplotype: N = the DSA window (100-bp flanks), T = the same with the ALT. The event haplotype is the one whose T window
   has a strict hit in the tumor assembly (`s1_map.sh MODE=tumor`; >= 90 % aligned, NM <= 1, no indel, the SNV base
   matched); if both do, the one with an exact (NM 0) hit; identical N windows on both haplotypes -> hapX.
2. `s2_giraffe.sh`: giraffe the N / T windows of the event haplotype to d9 and d46 (settings of the HG008T runs).
3. `s3_classify.py` (the code of `assembly_offref_20260929/s3_classify.py` + `E1.1/asm_d46_classify.py`): A = the somatic
   edit (T-only difference) on a GRCh38 node, B = on an off-reference node (an off-reference somatic SNV), C = the allele
   is a graph branch, D = unmapped / MAPQ < 20 / same alignment. B keys in node-forward coordinates.
4. `s4_scores.py` (GPU, tequila, `s4_scores_gpu.sh`): every B key with a tensor (exact candidate id, all autosomes) scored by
   the set's in-sample model: base (at the PASS threshold), rescue (predict's off-reference rescue: ch6 below the
   every-haplotype code raised to it at site +- 10, argmax) and all_H (ch6 = H on every valid cell, at the threshold),
   with the checkpoint's own normalization (current predict code, frozen copy).
5. `s5_reads.py` + `s5_map.sh`: the E9h read check, h50 strict rule: the A1-read consensus of the key's tensor (up to 101
   bp) aligned to the four haplotypes (minimap2 -ax sr); tumor_only = a strict hit in the tumor assembly and none in the
   normal (read-supported somatic); normal / ambiguous / neither as E9h (locus from the event windows' alignments).

## Results

**Event haplotype** (42,035 truth SNVs): placed on both DSA haplotypes with the REF base 40,902 / 40,951 (hapX / hapY);
not placed 882. Identical haplotype windows 27,488 (5,060 of them without a tumor hit), tumor hit on one haplotype 2,306,
exact tumor hit on one 3,441, ambiguous 5,393, unresolved 2,525 (hapX kept). Events: 41,153.

**Classes** (`classes_<graph>.tsv`; all 41,153 events / the 28,175 with a resolved event haplotype; chr1):

| graph | A | B (off-reference) | C (graph branch) | D | chr1 B |
|---|---|---|---|---|---|
| d9 (v1.1) | 40,843 / 27,960 | 26 (25 SNP + 1 INS edit) / 15 | 54 / 35 | 230 / 165 | 3 |
| d46 (v2.1) | 40,837 / 27,957 | 35 / 19 | 46 / 30 | 235 / 169 | 6 |

20 events are B on both graphs; 6 d9 B events are A on d46, 12 d46 B events are A on d9. The normal haplotype itself is
off-reference at the site for 3 (d9) and 12 (d46) B events. For comparison, HG008T (truth from the assemblies): 77 of
9,169 SNV events are B on d9 (0.84 %, 53 with the normal haplotype off-reference) against 26 of 41,153 here (0.06 %).

**Keys, scores and read support** (`tables.md`, `summary.json`, `results/key_classes.tsv`; all autosomes, chr1 in brackets):

| set | keys with a tensor (chr1) | read-supported somatic (tumor-only, h50 strict) (chr1) | base / rescue / all_H on them | rescue calls on all keys: read-supported / all |
|---|---|---|---|---|
| `COLO829T_Illumina` (d9) | 12 (2) | 7 (0) | 1 / 2 / 4 | 2 / 2 |
| `COLO829T_fiberseq` (d9) | 8 (0) | 7 (0) | 1 / 4 / 6 | 4 / 4 |
| `COLO829T_ONT` (d9) | 9 (1) | 7 (0) | 1 / 3 / 3 | 3 / 3 |
| `COLO829T_Illumina_d46` | 16 (6) | 12 (3) | 0 / 7 / 5 | 7 / 8 |
| `COLO829T_fiberseq_d46` | 12 (4) | 10 (3) | 0 / 4 / 4 | 4 / 4 |
| `COLO829T_ONT_d46` | 11 (3) | 10 (3) | 0 / 4 / 5 | 4 / 4 |

Of the 25 / 35 SNP keys, 8-16 per set have a tensor. Their labels: residual_partial_somatic_truth (an ignored SNV label)
4-11 per set, off_reference_no_truth_match (-1) 2-4, not_on_unique_grch38_node (-1) 1; none is a training positive. The
only base call is one key (`18425560:27:SNP:C>T`, read-supported) on all three d9 sets.

Same pattern as HG008T (E1.1 section 5, E9h): the base models call 0-1 of the read-supported off-reference somatic SNVs,
the rescue and all_H call several (2-7 and 3-6), and nearly every rescue call on these keys is read-supported
(24 of 25 across the six sets).

## Caveats

- Few events: the SMaHT COLO829T truth is a GRCh38 call set, so somatic SNVs in sequence where the patient's haplotype
  differs from GRCh38 are rare in it (0.06 % B against 0.84 % for the assembly-based HG008T truth). The numbers are
  small (8-16 keys per set) and describe the same 25-35 events on every platform.
- The tumor assembly was built from the COLO829T Fiber-seq HiFi and ONT D003 reads, so its read support is not
  independent of the Fiber-seq and ONT tensors (it is for Illumina). The event-haplotype choice uses the same assembly.
- 7,918 events have an ambiguous or unresolved event haplotype (hapX kept); restricted to resolved ones the counts are in
  the classes table (15 / 19 B).
- SNVs only; INDEL truths were not placed (repeat INDELs need the patient-frame array logic of
  `analysis/graph_absorbed_somatic_colo829t_20261005`).
- One seed per model; d46 models and the matched d9 fiber-seq / ONT models use norm train, the Illumina d9 model norm
  dataset (checkpoint statistics, matched by set name since the checkpoint records the RIS path).
