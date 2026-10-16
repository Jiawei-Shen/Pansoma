# E5b: somatic truth alleles absorbed by the graph, on HPRC v2.1 d46 and HPRC v1.1 d9, against the GRCh38-only arm; HCC1395T

**Status (2026-10-13, 23:00).**

- DONE: absorption classification on d9 and d46 for HG008T and COLO829T (3 platforms each), and on d46 for HCC1395T
  (s1-s5); PoN cross-tables; check of the graph-spelling search against the complete d9 edge table (t1/t2).
  `tables.md`, `summary.json` and `lost_alleles.tsv.gz` are final. They were written at 22:53, after both spelling
  outputs, and carry no PRELIMINARY banner.
- PENDING:
  - `linear_recovery.md` (s6): which lost / absorbed truth INDELs each arm calls after the PoN. It needs the 18 INDEL
    in-sample re-predictions (linear, d9 / joint d9, d46). They have not been submitted: they wait for the user's
    lab-cluster cap.
  - One formatting bug in `s5_tables.py` must be fixed before that collect (see Caveats, item 6).

## 1. Question and reviewer points

- **R1.4** (candidate generation removes somatic variants already encoded in the graph). A somatic allele that is
  already a graph path gives reads that align to it without edit. Such an allele produces no candidate, no tensor and
  no call. The earlier analyses quantified this on HPRC v1.1 d9 only:
  - `analysis/graph_absorbed_somatic_20261001` (HG008T);
  - `analysis/graph_absorbed_somatic_colo829t_20261005`;
  - `tmp/graph_absorbed_pon_20261012`.
- **R1.2** (composition and size of the pangenome). With HPRC2 (many more haplotypes), are more somatic alleles
  absorbed?
- E5b, as designed in `../../response_to_reviewers_detailed.md` (R1.4, E-table row E5b), redoes the classification
  as follows:
  - on the HPRC v2.1 d46 graph, next to d9;
  - against the GRCh38-only ("linear") arm of the same pipeline (7fc7396), which shows what the graph loses relative
    to a linear reference;
  - on HCC1395T, which has no assembly and only d46 tensor sets. Here the question is whether the truth ALT is a
    graph allele, and the PoN cross-table.
- Questions:
  1. How many truth alleles have a tensor on linear but none on d9 / d46? Of those, how many are absorbed, i.e. no
     candidate while a graph path spells the ALT?
  2. How many of them does the PoN tag anyway? The PoN tags truth alleles allele-exactly; the after-PoN loss is the
     loss specific to the graph.
  3. Does d46 absorb more or less than d9?
  4. Is absorption on HCC1395T (SEQC2 truth) like that on HG008T / COLO829T?
  5. (pending, s6) Which of the lost alleles does the linear arm actually call after the PoN?
- Not here: E5c, comparator recall at the absorbed loci. Comparators are not re-run (user decision 2026-10-13);
  E5c would use the user's existing ClairS-TO / DeepSomatic results.

## 2. Data

- Tensor sets (read-only): `/scratch/jshen/data/pansoma_v2_tensors/<SET>/`.
  - `tensors/somatic.recall.tsv` gives the per truth allele status: `tensor_representative` (T),
    `tensor_non_representative_allele` (A), `filtered` (F) or `no_candidate` (N). No other status occurs.
  - `truth/somatic.graph.tsv` gives the truth allele.
  - Arms:
    - d9 = `<SET>`;
    - d46 = `<SET>_d46`;
    - linear = `<SET>_linear`.
  - HG008T: Illumina / PacBio / ONT × (d9, d46, linear).
  - COLO829T: Illumina / fiberseq / ONT × (d9, d46, linear).
  - HCC1395T: Illumina / PacBio / ONT × d46 only.
- Truth alleles, chr1-22, all records (not only PASS / in the BED):

| sample | SNV | INS | DEL | records |
|---|---:|---:|---:|---:|
| HG008T | 8,690 | 4,951 | 3,545 | 17,186 |
| COLO829T | 42,035 | 888 | 1,024 | 43,947 |
| HCC1395T (SEQC2 v1.2.1) | 39,447 | 157 | 1,468 | 41,072 |

  Source: `status_matrix.tsv.gz`.
- Graphs (`e5b_common.py`):
  - d9 = `AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz`, H = 90;
  - d46 = `AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.gbz`, H = 464.
  - Both use the pipeline's `grch38_path` arrays for GRCh38 node placement and the `graph_index.sqlite` for haplotype
    counts.
- PoN tags of the truth: `scripts/filter_panel_of_normals.py`, rule of 7fc7396. A truth allele is tagged when
  `PANSOMA_PON` is set on its record.
  - HG008T: `tmp/graph_absorbed_pon_20261012/truth.pon.vcf.gz`.
  - COLO829T: `tmp/graph_absorbed_pon_20261012/colo.truth.pon.vcf.gz` (sites only).
  - HCC1395T: `hcc1395.truth.pon.vcf.gz`, made for this experiment (data dir, `hcc1395_pon.log`). It tagged 1,221 of
    41,072 records (SNV 1,102, INDEL 119): PoN1_gnomAD 828, PoN2_dbSNP 604, PoN3_1000G 94, PoN4_CoLoRSdb 701.
- Read-level d9 absorbed set, used for validation: `per_variant.tsv` of the two earlier analyses (HG008T 201 SNV +
  2,277 INDEL; COLO829T 67 SNV + 430 INDEL).
- Large data: `/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46/`.

## 3. Method

1. **s1 `s1_status_matrix.py`.** One row per (sample, truth allele) with its status on every platform × arm. The
   script asserts that truth_id, position and alleles are identical on every arm of a sample. Output:
   `status_matrix.tsv.gz`.
2. **s2 `s2_windows.py`.** One GRCh38 window per truth allele, as in `graph_absorbed_somatic_20261001` s1:
   - start from the GRCh38 nodes overlapping pos ± 60 bp;
   - widen to the periodic stretch around the event (period 1-6) ± 20 bp;
   - the first and last GRCh38 node are the anchors, separately on d9 and d46.

   Outputs: `loci.tsv.gz` (102,205 loci) and the seed node lists `window_nodes_{d9,d46}.txt` (370,316 and 531,796
   nodes).
3. **s3 `s3_extract.sh`.** `vg find -N <seeds> -c 250 -L` cuts out every node within 250 bp of a window node and the
   edges among them. vg 1.65.0 is used on both graphs: vg 1.77 crashes in `find` on d9, and also crashed on d46
   (job 387310). Subgraphs:
   - d9: 693,448 nodes, 672,127 edges;
   - d46: 817,172 nodes, 781,542 edges.
4. **s4 `s4_spell.py d9|d46`.** Does a graph path spell the truth ALT?
   - `ref_hap` = GRCh38 from the first to the last anchor; `alt_hap` = `ref_hap` with the truth allele applied.
   - The ALT is **spelled** when a path from (first,+) to (last,+), whose GRCh38 nodes all lie inside the window,
     spells `alt_hap` exactly. This is a prefix-pruned DFS: the "exact" class of `graph_absorbed_somatic_20261001`
     s2, without its germline-background and closest-allele searches.
   - `ref_ok` checks that the GRCh38 path spells `ref_hap`, i.e. that the subgraph is complete.
   - For spelled ALTs it records the departures from GRCh38 (`elements`) and the smallest haplotype count on the
     branch nodes (`alt_min_paths`).
   - Output: `spell_{d9,d46}.tsv.gz`.
5. **s5 `s5_tables.py`.** Classes per platform, for arm X ∈ {d9, d46} against linear:
   - `lost_X` = tensor (T) on linear, none on X. Lost alleles split by their status on X:
     - **absorbed** = no candidate on X and the ALT spelled by X;
     - **no candidate, ALT not spelled**;
     - **filtered**;
     - **other allele**.
   - `gained_X` = tensor on X, none on linear.
   - d46 vs d9 the same way (`d46_only` / `d9_only`).
   - "After PoN" = the truth allele is not PoN-tagged. It is the tensor-level ceiling of what a caller on that arm can
     still report after the PoN.
   - Outputs: `tables.md`, `summary.json`, `lost_alleles.tsv.gz` (alleles without a tensor on d9 or d46 on at least
     one platform) and `per_truth.tsv.gz` (data dir).
6. **Checks.**
   - **t1 `t1_check_spell_d9_sqlite.py`** runs the same search over the complete d9 node / edge table
     (`hprc-v1.1-mc-grch38.d9.GRCh38_CHM13.path_index.sqlite`). It covers the first 300 HG008T loci of the
     read-level absorbed set (`t1_spell_d9_sqlite.tsv`).
   - **t2 `t2_compare_d9.py`** compares t1 with s4's subgraph answer (`t2_compare_d9.txt`).
7. **s6 `s6_linear_recovery.py` (pending).** It joins `per_truth.tsv.gz` with the rtg `tp-baseline.vcf.gz` of the
   after-PoN evaluation in the benchmark BED:
   - `vcf_chr1*/eval_bed/vcfeval/pon_curve` (any score);
   - `vcf_chr1*/eval_bed/vcfeval/pon_calls` (PASS).

   The inputs are the 18 INDEL runs (linear, d9 or COLO829T joint d9, d46) under
   `/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/`. Scope: chr1 in-BED truth INDELs.
8. **`collect.sh`** runs t2, s5 and s6 on the login node. It needs both `spell_*.tsv.gz` files.

## 4. Scripts and jobs

Scripts: this folder. Slurm scripts and logs: `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E5b_absorbed_d46/`.

| step | job | resources / outcome |
|---|---|---|
| s1 status matrix | no Slurm job (no log; output 10:09) | - |
| s2 windows | 387156 (`s2_windows.sbatch`) | general, -c 2, 20 GB; 29 min, MaxRSS 6.7 GB; COMPLETED |
| vg test | 387157 | cancelled on maotai |
| vg test | 387252 (`vg_test.sbatch`) | vg 1.65 on d9 works (5 nodes / 4 edges). The d46 part was uninformative: test nodes 100-102 are not in d46. State FAILED from that part. |
| s3 d9 subgraph | 387290 (`s3_extract_d9.sbatch`) | vg 1.65.0; 50 min, 13.3 GB; COMPLETED |
| s3 d46 subgraph | 387310 | vg 1.77.0 crashed: `Path NA18505#2#JBHIJZ010000004.1#33675 is not indexed for get_position_of_step`. FAILED. |
| s3 d46 subgraph | 387428 (`s3_extract_d46.sbatch`) | vg 1.65.0; 7 h 12 min, 25.7 GB; COMPLETED |
| s4 spell d46 | 387429 (`s4_spell_d46.sbatch`, after 387428) | 3 min, 5.8 GB; COMPLETED |
| s4 spell d9 | 387430 (`s4_spell_d9.sbatch`, after 387429) | 2.7 min, 2.6 GB; COMPLETED |
| HCC1395T truth PoN | no Slurm file in the job dir | `hcc1395_pon.log` (data dir), 34 min |
| t1 / t2, s5, s6 | login node (`collect.sh`) | `t2_compare_d9.txt`, `tables.md`, `summary.json` and `linear_recovery.md` written 22:53 |

`s5_prelim.stderr.txt` is the stderr of a preliminary s5 run at 11:49, before the spelling outputs existed. It is
superseded.

## 5. Results

All tables are copied from `tables.md` (section numbers as there); `summary.json` holds the same counts. "Tagged" =
PoN-tagged truth allele; "after PoN" = not tagged.

### 5.1 Truth alleles with a tensor, per arm (tables.md §1, chr1-22, all truth records)

| sample | platform | kind | n (not tagged) | tensor: linear / d9 / d46 | tensor and not tagged: linear / d9 / d46 | d9 - linear, d46 - linear (not tagged) |
|---|---|---|---|---|---|---|
| HG008T | Illumina | SNV | 8,690 (7,780) | 8,552 (0.984) / 8,273 (0.952) / 8,307 (0.956) | 7,714 (0.992) / 7,684 (0.988) / 7,686 (0.988) | -30 (-0.3% of truth), -28 (-0.3%) |
| HG008T | Illumina | INDEL | 8,496 (3,251) | 6,720 (0.791) / 3,555 (0.418) / 3,821 (0.450) | 2,683 (0.825) / 2,276 (0.700) / 2,306 (0.709) | -407 (-4.8% of truth), -377 (-4.4%) |
| HG008T | PacBio | SNV | 8,690 (7,780) | 8,545 (0.983) / 8,231 (0.947) / 8,249 (0.949) | 7,704 (0.990) / 7,681 (0.987) / 7,680 (0.987) | -23 (-0.3% of truth), -24 (-0.3%) |
| HG008T | PacBio | INDEL | 8,496 (3,251) | 7,077 (0.833) / 4,062 (0.478) / 4,226 (0.497) | 2,807 (0.863) / 2,294 (0.706) / 2,313 (0.711) | -513 (-6.0% of truth), -494 (-5.8%) |
| HG008T | ONT | SNV | 8,690 (7,780) | 8,567 (0.986) / 8,268 (0.951) / 8,290 (0.954) | 7,741 (0.995) / 7,718 (0.992) / 7,725 (0.993) | -23 (-0.3% of truth), -16 (-0.2%) |
| HG008T | ONT | INDEL | 8,496 (3,251) | 4,848 (0.571) / 3,628 (0.427) / 3,642 (0.429) | 2,241 (0.689) / 2,074 (0.638) / 2,075 (0.638) | -167 (-2.0% of truth), -166 (-2.0%) |
| COLO829T | Illumina | SNV | 42,035 (40,593) | 38,222 (0.909) / 37,864 (0.901) / 37,861 (0.901) | 36,871 (0.908) / 36,608 (0.902) / 36,598 (0.902) | -263 (-0.6% of truth), -273 (-0.6%) |
| COLO829T | Illumina | INDEL | 1,912 (929) | 1,747 (0.914) / 1,092 (0.571) / 1,155 (0.604) | 869 (0.935) / 739 (0.795) / 751 (0.808) | -130 (-6.8% of truth), -118 (-6.2%) |
| COLO829T | fiberseq | SNV | 42,035 (40,593) | 37,945 (0.903) / 37,756 (0.898) / 37,736 (0.898) | 36,635 (0.902) / 36,514 (0.900) / 36,491 (0.899) | -121 (-0.3% of truth), -144 (-0.3%) |
| COLO829T | fiberseq | INDEL | 1,912 (929) | 1,842 (0.963) / 1,210 (0.633) / 1,260 (0.659) | 910 (0.980) / 753 (0.811) / 755 (0.813) | -157 (-8.2% of truth), -155 (-8.1%) |
| COLO829T | ONT | SNV | 42,035 (40,593) | 37,665 (0.896) / 37,523 (0.893) / 37,514 (0.892) | 36,354 (0.896) / 36,276 (0.894) / 36,263 (0.893) | -78 (-0.2% of truth), -91 (-0.2%) |
| COLO829T | ONT | INDEL | 1,912 (929) | 1,441 (0.754) / 1,108 (0.579) / 1,116 (0.584) | 789 (0.849) / 701 (0.755) / 695 (0.748) | -88 (-4.6% of truth), -94 (-4.9%) |
| HCC1395T | Illumina | SNV | 39,447 (38,345) | - / - / 35,266 (0.894) | - / - / 34,455 (0.899) | - |
| HCC1395T | Illumina | INDEL | 1,625 (1,506) | - / - / 1,346 (0.828) | - / - / 1,280 (0.850) | - |
| HCC1395T | PacBio | SNV | 39,447 (38,345) | - / - / 34,374 (0.871) | - / - / 33,604 (0.876) | - |
| HCC1395T | PacBio | INDEL | 1,625 (1,506) | - / - / 1,375 (0.846) | - / - / 1,307 (0.868) | - |
| HCC1395T | ONT | SNV | 39,447 (38,345) | - / - / 30,869 (0.783) | - / - / 30,125 (0.786) | - |
| HCC1395T | ONT | INDEL | 1,625 (1,506) | - / - / 1,239 (0.762) | - / - / 1,171 (0.778) | - |

### 5.2 INDELs lost on a pangenome arm against linear, by class (tables.md §2, INDEL rows; SNV rows in tables.md)

Each cell: alleles (PoN-tagged, share tagged; not tagged = graph-specific loss after the PoN).

| sample | platform | kind | arm | lost (all) | absorbed (no candidate, ALT spelled by the graph) | no candidate, ALT not spelled | filtered | other allele | gained (tensor on arm, not on linear) |
|---|---|---|---|---|---|---|---|---|---|
| HG008T | Illumina | INDEL | d9 | 3,341 (2,860 tagged, 85.6 %; 481 after PoN) | 1,613 (1,578 tagged, 97.8 %; 35 after PoN) | 987 (713 tagged, 72.2 %; 274 after PoN) | 681 (528 tagged, 77.5 %; 153 after PoN) | 60 (41 tagged, 68.3 %; 19 after PoN) | 176 (102 tagged, 58.0 %; 74 after PoN) |
| HG008T | Illumina | INDEL | d46 | 3,104 (2,640 tagged, 85.1 %; 464 after PoN) | 1,301 (1,262 tagged, 97.0 %; 39 after PoN) | 1,139 (859 tagged, 75.4 %; 280 after PoN) | 577 (454 tagged, 78.7 %; 123 after PoN) | 87 (65 tagged, 74.7 %; 22 after PoN) | 205 (118 tagged, 57.6 %; 87 after PoN) |
| HG008T | PacBio | INDEL | d9 | 3,229 (2,648 tagged, 82.0 %; 581 after PoN) | 935 (901 tagged, 96.4 %; 34 after PoN) | 1,203 (864 tagged, 71.8 %; 339 after PoN) | 835 (710 tagged, 85.0 %; 125 after PoN) | 256 (173 tagged, 67.6 %; 83 after PoN) | 214 (146 tagged, 68.2 %; 68 after PoN) |
| HG008T | PacBio | INDEL | d46 | 3,078 (2,507 tagged, 81.4 %; 571 after PoN) | 839 (807 tagged, 96.2 %; 32 after PoN) | 1,433 (1,072 tagged, 74.8 %; 361 after PoN) | 546 (451 tagged, 82.6 %; 95 after PoN) | 260 (177 tagged, 68.1 %; 83 after PoN) | 227 (150 tagged, 66.1 %; 77 after PoN) |
| HG008T | ONT | INDEL | d9 | 1,847 (1,548 tagged, 83.8 %; 299 after PoN) | 441 (421 tagged, 95.5 %; 20 after PoN) | 583 (438 tagged, 75.1 %; 145 after PoN) | 591 (507 tagged, 85.8 %; 84 after PoN) | 232 (182 tagged, 78.4 %; 50 after PoN) | 627 (495 tagged, 78.9 %; 132 after PoN) |
| HG008T | ONT | INDEL | d46 | 1,752 (1,462 tagged, 83.4 %; 290 after PoN) | 472 (452 tagged, 95.8 %; 20 after PoN) | 687 (541 tagged, 78.7 %; 146 after PoN) | 364 (304 tagged, 83.5 %; 60 after PoN) | 229 (165 tagged, 72.1 %; 64 after PoN) | 546 (422 tagged, 77.3 %; 124 after PoN) |
| COLO829T | Illumina | INDEL | d9 | 666 (533 tagged, 80.0 %; 133 after PoN) | 354 (350 tagged, 98.9 %; 4 after PoN) | 199 (107 tagged, 53.8 %; 92 after PoN) | 103 (70 tagged, 68.0 %; 33 after PoN) | 10 (6 tagged, 60.0 %; 4 after PoN) | 11 (8 tagged, 72.7 %; 3 after PoN) |
| COLO829T | Illumina | INDEL | d46 | 609 (485 tagged, 79.6 %; 124 after PoN) | 315 (309 tagged, 98.1 %; 6 after PoN) | 194 (105 tagged, 54.1 %; 89 after PoN) | 88 (65 tagged, 73.9 %; 23 after PoN) | 12 (6 tagged, 50.0 %; 6 after PoN) | 17 (11 tagged, 64.7 %; 6 after PoN) |
| COLO829T | fiberseq | INDEL | d9 | 642 (483 tagged, 75.2 %; 159 after PoN) | 183 (180 tagged, 98.4 %; 3 after PoN) | 203 (112 tagged, 55.2 %; 91 after PoN) | 219 (173 tagged, 79.0 %; 46 after PoN) | 37 (18 tagged, 48.6 %; 19 after PoN) | 10 (8 tagged, 80.0 %; 2 after PoN) |
| COLO829T | fiberseq | INDEL | d46 | 591 (434 tagged, 73.4 %; 157 after PoN) | 198 (193 tagged, 97.5 %; 5 after PoN) | 214 (123 tagged, 57.5 %; 91 after PoN) | 142 (96 tagged, 67.6 %; 46 after PoN) | 37 (22 tagged, 59.5 %; 15 after PoN) | 9 (7 tagged, 77.8 %; 2 after PoN) |
| COLO829T | ONT | INDEL | d9 | 444 (330 tagged, 74.3 %; 114 after PoN) | 115 (114 tagged, 99.1 %; 1 after PoN) | 142 (75 tagged, 52.8 %; 67 after PoN) | 144 (116 tagged, 80.6 %; 28 after PoN) | 43 (25 tagged, 58.1 %; 18 after PoN) | 111 (85 tagged, 76.6 %; 26 after PoN) |
| COLO829T | ONT | INDEL | d46 | 429 (313 tagged, 73.0 %; 116 after PoN) | 128 (126 tagged, 98.4 %; 2 after PoN) | 149 (84 tagged, 56.4 %; 65 after PoN) | 104 (71 tagged, 68.3 %; 33 after PoN) | 48 (32 tagged, 66.7 %; 16 after PoN) | 104 (82 tagged, 78.8 %; 22 after PoN) |

SNV absorbed alleles (tables.md §2) are PoN-tagged in every row: 0 left after the PoN, except 1 on COLO829T Illumina
d9. The SNVs lost after the PoN fall in two classes, per platform and arm:

- "no candidate, ALT not spelled": HG008T 15-27, COLO829T 47-257;
- "filtered": HG008T 6-13, COLO829T 30-105.

### 5.3 Union over the three platforms (tables.md §3)

| sample | kind | arm | truth n (tagged) | lost on >= 1 platform | absorbed on >= 1 platform | lost on all 3 | absorbed on all 3 | lost on all 3, not absorbed anywhere | gained on all 3 (tensor on the arm, not on linear) |
|---|---|---|---|---|---|---|---|---|---|
| HG008T | SNV | d9 | 8,690 (910 tagged, 10.5 %; 7,780 after PoN) | 426 (359 tagged, 84.3 %; 67 after PoN) | 184 (184 tagged, 100.0 %; 0 after PoN) | 200 (193 tagged, 96.5 %; 7 after PoN) | 113 (113 tagged, 100.0 %; 0 after PoN) | 63 (56 tagged, 88.9 %; 7 after PoN) | 0 (0 tagged, -; 0 after PoN) |
| HG008T | SNV | d46 | 8,690 (910 tagged, 10.5 %; 7,780 after PoN) | 396 (333 tagged, 84.1 %; 63 after PoN) | 140 (140 tagged, 100.0 %; 0 after PoN) | 180 (172 tagged, 95.6 %; 8 after PoN) | 91 (91 tagged, 100.0 %; 0 after PoN) | 74 (66 tagged, 89.2 %; 8 after PoN) | 0 (0 tagged, -; 0 after PoN) |
| HG008T | INDEL | d9 | 8,496 (5,245 tagged, 61.7 %; 3,251 after PoN) | 4,275 (3,555 tagged, 83.2 %; 720 after PoN) | 1,818 (1,759 tagged, 96.8 %; 59 after PoN) | 1,188 (1,013 tagged, 85.3 %; 175 after PoN) | 315 (309 tagged, 98.1 %; 6 after PoN) | 601 (443 tagged, 73.7 %; 158 after PoN) | 30 (18 tagged, 60.0 %; 12 after PoN) |
| HG008T | INDEL | d46 | 8,496 (5,245 tagged, 61.7 %; 3,251 after PoN) | 4,039 (3,335 tagged, 82.6 %; 704 after PoN) | 1,467 (1,416 tagged, 96.5 %; 51 after PoN) | 1,126 (956 tagged, 84.9 %; 170 after PoN) | 345 (333 tagged, 96.5 %; 12 after PoN) | 628 (473 tagged, 75.3 %; 155 after PoN) | 39 (25 tagged, 64.1 %; 14 after PoN) |
| COLO829T | SNV | d9 | 42,035 (1,442 tagged, 3.4 %; 40,593 after PoN) | 489 (113 tagged, 23.1 %; 376 after PoN) | 64 (63 tagged, 98.4 %; 1 after PoN) | 119 (51 tagged, 42.9 %; 68 after PoN) | 41 (41 tagged, 100.0 %; 0 after PoN) | 72 (4 tagged, 5.6 %; 68 after PoN) | 2 (0 tagged, 0.0 %; 2 after PoN) |
| COLO829T | SNV | d46 | 42,035 (1,442 tagged, 3.4 %; 40,593 after PoN) | 552 (114 tagged, 20.7 %; 438 after PoN) | 54 (54 tagged, 100.0 %; 0 after PoN) | 101 (45 tagged, 44.6 %; 56 after PoN) | 35 (35 tagged, 100.0 %; 0 after PoN) | 65 (9 tagged, 13.8 %; 56 after PoN) | 3 (1 tagged, 33.3 %; 2 after PoN) |
| COLO829T | INDEL | d9 | 1,912 (983 tagged, 51.4 %; 929 after PoN) | 799 (615 tagged, 77.0 %; 184 after PoN) | 382 (376 tagged, 98.4 %; 6 after PoN) | 323 (251 tagged, 77.7 %; 72 after PoN) | 89 (89 tagged, 100.0 %; 0 after PoN) | 150 (80 tagged, 53.3 %; 70 after PoN) | 5 (5 tagged, 100.0 %; 0 after PoN) |
| COLO829T | INDEL | d46 | 1,912 (983 tagged, 51.4 %; 929 after PoN) | 757 (577 tagged, 76.2 %; 180 after PoN) | 347 (340 tagged, 98.0 %; 7 after PoN) | 298 (227 tagged, 76.2 %; 71 after PoN) | 99 (97 tagged, 98.0 %; 2 after PoN) | 148 (80 tagged, 54.1 %; 68 after PoN) | 6 (6 tagged, 100.0 %; 0 after PoN) |

### 5.4 Truth ALT spelled by a graph path at its locus, independent of the reads (tables.md §5)

| sample | kind | n | spelled by d9 | spelled by d46 | d9 only | d46 only | both | spelled by d46, PoN-tagged |
|---|---|---:|---|---|---:|---:|---:|---|
| HG008T | SNV | 8,690 | 215 (2.5 %) | 153 (1.8 %) | 87 | 25 | 128 | 153 (153 tagged, 100.0 %; 0 after PoN) |
| HG008T | INDEL | 8,496 | 2,319 (27.3 %) | 1,809 (21.3 %) | 716 | 206 | 1,603 | 1,809 (1,740 tagged, 96.2 %; 69 after PoN) |
| COLO829T | SNV | 42,035 | 83 (0.2 %) | 70 (0.2 %) | 18 | 5 | 65 | 70 (69 tagged, 98.6 %; 1 after PoN) |
| COLO829T | INDEL | 1,912 | 457 (23.9 %) | 398 (20.8 %) | 112 | 53 | 345 | 398 (389 tagged, 97.7 %; 9 after PoN) |
| HCC1395T | SNV | 39,447 | 202 (0.5 %) | 228 (0.6 %) | 9 | 35 | 193 | 228 (228 tagged, 100.0 %; 0 after PoN) |
| HCC1395T | INDEL | 1,625 | 52 (3.2 %) | 47 (2.9 %) | 8 | 3 | 44 | 47 (47 tagged, 100.0 %; 0 after PoN) |

Validation against the read-level d9 absorbed set of the earlier analyses. That set = ALT reads align with no edit to a
d9 branch / skip edge, union of 3 platforms, on the 2026-09-30 tensor sets.

| sample | kind | read-level absorbed (d9) | of which ALT spelled by d9 (this check) | of which no candidate on d9 on >= 1 platform now | this check: no candidate on d9 on >= 1 platform and spelled | of which in the read-level set |
|---|---|---:|---:|---:|---:|---:|
| HG008T | SNV | 201 | 145 | 201 | 188 | 145 |
| HG008T | INDEL | 2,277 | 1,839 | 2,277 | 2,039 | 1,839 |
| COLO829T | SNV | 67 | 64 | 67 | 66 | 64 |
| COLO829T | INDEL | 430 | 382 | 430 | 402 | 382 |

Spelling check on the complete d9 edge table (`t2_compare_d9.txt`; 300 HG008T read-level loci):

```
(ref_ok sqlite, ref_ok vg-find, alt sqlite, alt vg-find): count over 300 loci
('1', '1', '0', '0') 44
('1', '1', '1', '0') 22
('1', '1', '1', '1') 234
agree on alt_spelled: 278 of 300
```

By the earlier analysis' graph_match class (`t1_spell_d9_sqlite.tsv`), the complete table spells:

- all 256 "exact" loci;
- none of the 32 "closest" loci;
- none of the 12 "with_germline" loci.

### 5.5 HCC1395T on d46, and the same absorbed share on every sample (tables.md §7)

| sample | platform | SNV absorbed (share) | INDEL absorbed (share) | INDEL absorbed after PoN (share of truth INDELs) |
|---|---|---|---|---|
| HG008T | Illumina | 111 (1.3 %) | 1,576 (18.5 %) | 61 (0.7 %) |
| HG008T | PacBio | 128 (1.5 %) | 998 (11.7 %) | 44 (0.5 %) |
| HG008T | ONT | 131 (1.5 %) | 833 (9.8 %) | 38 (0.4 %) |
| COLO829T | Illumina | 48 (0.1 %) | 363 (19.0 %) | 7 (0.4 %) |
| COLO829T | fiberseq | 55 (0.1 %) | 220 (11.5 %) | 5 (0.3 %) |
| COLO829T | ONT | 54 (0.1 %) | 172 (9.0 %) | 4 (0.2 %) |
| HCC1395T | Illumina | 225 (0.6 %) | 46 (2.8 %) | 0 (0.0 %) |
| HCC1395T | PacBio | 225 (0.6 %) | 33 (2.0 %) | 0 (0.0 %) |
| HCC1395T | ONT | 226 (0.6 %) | 26 (1.6 %) | 0 (0.0 %) |

Absorbed here = no candidate on d46 and the ALT spelled by d46, as a share of all truth alleles of the kind. Unlike
§5.2, it does not require a tensor on linear. HCC1395T per class (tables.md §7):

- no candidate: SNV 436 / 824 / 2,481 and INDEL 77 / 75 / 130 (Illumina / PacBio / ONT);
- filtered: SNV 3,745 / 4,248 / 6,095 and INDEL 198 / 166 / 235.

Most HCC1395T alleles without a tensor are therefore filtered, not absorbed.

### 5.6 chr1 in the benchmark BED, INDEL, after the PoN (tables.md §8; the vcfeval test scope)

| sample | platform | truth INDEL (not tagged) | tensor & not tagged: linear / d9 / d46 | lost on d9 vs linear, not tagged (absorbed) | lost on d46 vs linear, not tagged (absorbed) |
|---|---|---|---|---|---|
| HG008T | Illumina | 588 (209) | 170 / 151 / 153 | 27 (0) | 25 (2) |
| HG008T | PacBio | 588 (209) | 181 / 153 / 153 | 30 (1) | 31 (2) |
| HG008T | ONT | 588 (209) | 157 / 148 / 145 | 18 (1) | 20 (2) |
| COLO829T | Illumina | 131 (50) | 47 / 36 / 36 | 11 (1) | 11 (1) |
| COLO829T | fiberseq | 131 (50) | 50 / 37 / 37 | 13 (1) | 13 (1) |
| COLO829T | ONT | 131 (50) | 44 / 34 / 34 | 10 (0) | 10 (1) |

### 5.7 d46 against d9 (tables.md §4)

Net d46 − d9 tensors after the PoN (column "net d46 - d9 after PoN"):

| sample | kind | Illumina | PacBio / fiberseq | ONT |
|---|---|---:|---:|---:|
| HG008T | SNV | +2 | -1 | +7 |
| HG008T | INDEL | +30 | +19 | +1 |
| COLO829T | SNV | -10 | -23 | -13 |
| COLO829T | INDEL | +12 | +2 | -6 |

HG008T INDELs with a tensor on d46 but not on d9 that are absorbed on d9: 266 / 132 / 92. After the PoN, 20 / 12 / 8
of them remain. The reverse direction (tensor on d9 only, absorbed on d46) has 64 / 76 / 106, and after the PoN
11 / 9 / 4.

### 5.8 Lost truth INDELs called after the PoN (s6, `linear_recovery.md`): PENDING

The populations are counted: chr1, in the BED, lost / absorbed on d9 / d46 against linear. All 36 rtg outputs are
missing (18 runs × `pon_curve` / `pon_calls`), so every TP column reads "pending". For COLO829T, the "absorbed" n of
the current file is wrong (see Caveats, item 6).

## 6. Interpretation

1. **Absorption exists on d46 as on d9, and is almost entirely PoN-tagged.**
   - Absorbed truth alleles = no candidate on the arm, a tensor on linear, and the ALT spelled by the graph (§5.2,
     tables.md §2). Of these, 95.5-97.8 % of HG008T INDELs, 97.5-99.1 % of COLO829T INDELs and 98.2-100 % of SNVs
     are PoN-tagged.
   - After the PoN, per platform:

     | sample | arm | absorbed INDELs left (Illumina / PacBio or fiberseq / ONT) |
     |---|---|---|
     | HG008T | d9 | 35 / 34 / 20 |
     | HG008T | d46 | 39 / 32 / 20 |
     | COLO829T | d9 | 4 / 3 / 1 |
     | COLO829T | d46 | 6 / 5 / 2 |

     No absorbed SNV is left, except 1 on COLO829T Illumina d9.
   - Union over platforms: 59 (d9) and 51 (d46) HG008T INDELs, 6 and 7 COLO829T INDELs (§5.3).
   - On d46, absorbed-after-PoN INDELs are 0.4-0.7 % of HG008T truth INDELs, 0.2-0.4 % of COLO829T's and 0.0 % of
     HCC1395T's (§5.5).
   - The reviewer's mechanism is real. The loss specific to the graph, beyond what the population filter removes
     anyway, is a fraction of a percent of truth INDELs.
2. **The larger graph does not absorb more (R1.2).**
   - d46 spells fewer truth ALTs than d9 (§5.4):
     - HG008T INDEL: 21.3 % against 27.3 %;
     - COLO829T INDEL: 20.8 % against 23.9 %;
     - HG008T SNV: 1.8 % against 2.5 %.
   - d46 has more INDEL tensors after the PoN than d9 on 5 of 6 sets (net +1 to +30; §5.7).
   - SNVs are within ±23 tensors, and slightly lower on d46 for COLO829T.
   - Both graphs are AF-filtered (H = 90 and 464). The filter, not the haplotype count, appears to set what is
     absorbed. That reading comes from the design document; this folder does not test it.
3. **Absorption explains only a small part of what the graph arms lose against linear.**
   - Tensor-level ceiling after the PoN, d9 − linear (§5.1):
     - HG008T INDEL: −407 / −513 / −167 (−4.8 / −6.0 / −2.0 % of truth);
     - COLO829T INDEL: −130 / −157 / −88;
     - d46 is similar.
   - Most of the after-PoN loss is "no candidate, ALT not spelled" (HG008T 145-361 per platform and arm), "filtered"
     (60-153) and "other allele" (19-83); absorbed is 20-39 (§5.2).
   - In the vcfeval scope (chr1, BED, §5.6), the graph arms lose 18-31 HG008T and 10-13 COLO829T INDEL tensors after
     the PoN against linear, of which 0-2 are absorbed.
   - The remainder is a different representation / alignment of repeat INDELs on the graph (or partial absorption
     that the exact-spelling test cannot see; Caveats 1-2), not absorption in the strict sense. s6 will show whether
     the linear arm actually calls these alleles.
4. **HCC1395T absorbs little.**
   - d46 spells 2.9 % of SEQC2 truth INDELs and 0.6 % of SNVs, against 21.3 % and 20.8 % of HG008T / COLO829T truth
     INDELs.
   - Every d46-spelled HCC1395T allele is PoN-tagged (228 / 228 SNV, 47 / 47 INDEL), so none of them is a
     graph-specific loss after the PoN.
   - The HCC1395T truth is also far less PoN-tagged: INDELs 7.3 %, against 61.7 % (HG008T) and 51.4 % (COLO829T).
   - The absorbed fraction therefore depends on the truth set's composition. Why the SEQC2 truth differs (for
     example fewer long-homopolymer INDELs) is not tested here.
5. **Validation.**
   - The read-independent spelling test recovers 80.8 % (1,839 / 2,277) of the HG008T read-level absorbed INDELs and
     88.8 % (382 / 430) of COLO829T's.
   - On the t1 sample, the complete d9 edge table spells exactly the "exact" class (256 / 256) and none of the
     "closest" / "with_germline" alleles.
   - The rest of the gap comes from the non-exact classes (Caveat 1) and the incomplete subgraph (Caveat 2).

## 7. Caveats

1. **"Spelled" is a strict test.** It asks for the exact GRCh38 + truth haplotype between the window anchors. It
   misses:
   - alleles whose graph path needs nearby germline alleles of the patient ("with_germline");
   - near matches ("closest").

   These alleles are counted as "no candidate, ALT not spelled". On the 300 t1 loci that is 44 of 300 read-level
   absorbed alleles.
2. **The vg-find subgraph (250 bp context) is incomplete.** Of 256 loci spelled with the complete d9 edge table, 22
   (8.6 %) are not spelled on the subgraph (t2). The d9 "absorbed" counts are therefore lower bounds. d46 has no
   equivalent check: no complete d46 edge table was used.
   - The GRCh38 sanity check (`ref_ok`) fails at only 2 loci per graph, and the same 2 on both: HG008T truth_id 9463
     and COLO829T truth_id 28403, both SNVs on chr10.
   - Most "not spelled" loci are `first_mismatch` (d9 57,053, d46 40,058 of 102,205): the allele lies inside the
     first anchor node, so no graph path can carry it. This is correct, not a failure.
3. **The classes come from tensor-set statuses.** "Lost after PoN" is a tensor-level ceiling, not called TPs. Called
   TPs are s6 (pending).
4. **Arms differ in more than the graph.** Linear, d9 and d46 use separate alignments (GAMs) of the same reads. The
   "no candidate, ALT not spelled" class mixes alignment and representation differences in repeats; it is not broken
   down here.
5. **Scope.**
   - §5.1-5.5 cover all truth records on chr1-22, not only PASS / in the BED.
   - §5.6 is the vcfeval scope (chr1, in the BED).
   - The read-level validation set was defined on the 2026-09-30 tensor sets. Some of its alleles have a tensor on
     d9 now: HG008T PacBio 352 and ONT 423 INDELs (tables.md §2b).
6. **Bug in the TSV outputs (counts in `tables.md` / `summary.json` are unaffected).**
   - The cause: `s5_tables.py` `pon_tags()` keeps the trailing newline of the INFO column, and the COLO829T PoN VCF
     is sites-only.
   - The effect: in `per_truth.tsv.gz` and `lost_alleles.tsv.gz`, every PoN-tagged COLO829T record (2,425 =
     1,442 SNV + 983 INDEL) is split over two lines (22 + 11 fields).
   - `tables.md` / `summary.json` are computed in memory, and the tag stays non-empty, so they are correct.
   - `s6_linear_recovery.py` reads `per_truth.tsv.gz`, so the COLO829T "absorbed on d9 / d46" n in
     `linear_recovery.md` misses every tagged allele. Recount by rejoining the split lines (ad-hoc check, same
     definition):

     | platform | absorbed on d9: file / rejoined | absorbed on d46: file / rejoined |
     |---|---|---|
     | Illumina | 1 / 26 | 1 / 23 |
     | fiberseq | 1 / 21 | 1 / 21 |
     | ONT | 0 / 16 | 1 / 18 |

     The "not PoN-tagged" column and the "lost" rows are right.
   - **Fixed 2026-10-13**: `pon_tags()` strips the newline (`line.rstrip('\n')`) and `collect.sh` was rerun.
     `linear_recovery.md` now has the rejoined counts (26 / 23, 21 / 21, 16 / 18), and every line of
     `lost_alleles.tsv.gz` has 32 fields.
7. Only the main PoN rule of 7fc7396 is used. Other floors are not tested.

## 8. Status

DONE (2026-10-13):

- s1 `status_matrix.tsv.gz`, s2 `loci.tsv.gz` (job 387156), s3 subgraphs d9 (387290) and d46 (387428, after the vg
  1.77 crash in 387310).
- s4 spelling on d46 (387429) and d9 (387430).
- HCC1395T truth PoN.
- t1 / t2 checks.
- s5 final `tables.md`, `summary.json`, `lost_alleles.tsv.gz` and `per_truth.tsv.gz`.
- No E5b job is running or held.

PENDING:

1. **s6 `linear_recovery.md`: called TPs among the lost / absorbed INDELs.**
   - It needs the 18 in-sample INDEL re-predictions on tequila GPUs (`tmp/reviewer_cpu_20261013/task_gpu.sh`; rows
     `insample_linear`, `insample_d9`, `joint_d9`, `insample_d46` of `tasks_full_plan.tsv`):
     - `{HG008,COLO829}_*_INDEL_linear_7fc7396`;
     - `HG008_*_INDEL_7fc7396`;
     - `COLO829_*_INDEL_joint_7fc7396`;
     - `{HG008,COLO829}_*_INDEL_d46_7fc7396`.
   - They are not in the submitted first small batch (`tasks.tsv`, GPU array 387379, which holds SNV d46 and
     cross-sample tasks only).
   - Reason: the lab-cluster footprint. New work goes in small batches under a cap the user has not yet agreed;
     the held reviewer jobs of other experiments stay held until then.
   - Then run `./collect.sh` again (login node, a few minutes; the `pon_tags()` bug of Caveat 6 is fixed).
2. **E5c** (comparator calls at the absorbed loci) is a separate item. It uses the user's existing ClairS-TO /
   DeepSomatic-TO results; comparators are not re-run.
3. **Optional, not planned:** a d46 completeness check like t1, which needs a complete d46 edge table.
