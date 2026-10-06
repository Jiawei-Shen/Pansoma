# Repeat context of the HG008T truth INDELs and of Pansoma's INDEL calls (chr1 test)

2026-10-06. Question: how many of the HG008T truth INDELs are homopolymer / STR events, and what share of Pansoma's INDEL
predictions fall in the same classes. Script `repeat_context.py` (one run on the login node, ~6 min), tables in `tables.md`,
one row per evaluated call in `calls_context.tsv`.

- **Class rule** (the b1 rule of `analysis/graph_absorbed_somatic_20261001`, applied in GRCh38 to truth and calls alike):
  unit = minimal period of the inserted / deleted bases; tract = longest exact period-unit stretch touching the event;
  `HP>=7` homopolymer of >= 7 bp (unit 1), `HP4-6`, `STR` = unit 2-6 with >= 3 copies, `VNTR>6`, `none` = not in a repeat.
  The truth is also given in the **patient frame** (the GIAB INFO event on the HG008-N v6.2 event haplotype, b5 rule),
  where GRCh38-frame 'none' records mostly turn out to be repeat events or substitutions inside repeats.
- **Calls:** the seven current HG008T Illumina INDEL runs (`/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_INDEL_*`;
  there is no PacBio / ONT INDEL model). PASS calls on chr1 inside the GIAB BED, GRCh38-placed, split into TP / FP by the
  runs' rtg vcfeval outputs; `raw` = before the PoN, `pon` = after the 2026-10-06 INDEL PoN (gnomAD / CoLoRSdb exact allele
  AF >= 0.01). Off-reference PASS calls (no GRCh38 allele) are counted but not classed.
- **Status of the PoN numbers:** complete for the three `*_small*` runs (rtg outputs of 2026-10-06 13:09). The
  `vcf_chr1_pon` / `vcf_chr1_r09_pon` directories of `base`, `base_b1024`, `nopartial_b1024` and `nopartial_b1024_nopc` were
  removed by another session at 13:07 (the INDEL-PoN evaluations are being redone there); re-run `repeat_context.py` when
  they are back. The raw numbers of all seven runs are complete.

## Truth INDELs (chr1-22: 8,496 alleles; chr1 BED: 588)

| class | chr1-22, GRCh38 frame | chr1 BED, GRCh38 frame | chr1-22, patient frame | chr1 BED, patient frame |
|---|---:|---:|---:|---:|
| homopolymer `HP>=7` | 6,360 (74.9%) | 436 (74.1%) | 6,722 (79.1%) | 458 (77.9%) |
| `STR` (unit 2-6, >= 3 copies) | 936 (11.0%) | 65 (11.1%) | 795 (9.4%) | 66 (11.2%) |
| `HP4-6` | 241 (2.8%) | 7 (1.2%) | 204 (2.4%) | 7 (1.2%) |
| `VNTR>6` / `in_STR` | 6 | 0 | 15 + 27 | 0 + 1 |
| substitution inside a repeat (GRCh38 wrote it as an INDEL) | | | 253 (3.0%) | 18 (3.1%) |
| non-repeat (`none`, substitution) | 952 (11.2%) | 80 (13.6%) | 439 (5.2%) | 36 (6.1%) |
| no INFO event | | | 41 | 2 |

Homopolymer + STR = 86% of the truth INDELs in GRCh38, 89% in the patient frame.

## Pansoma INDEL calls on chr1 (in the BED)

Share of the PASS calls per class, and precision / recall per class (`nopartial_b1024_nopc_small` shown; all runs in `tables.md`):

| class | raw calls (share) | raw TP / FP | raw precision | after PoN (share) | pon TP / FP | pon precision | truth chr1 BED | recall raw / pon |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `HP>=7` | 2,242 (63.1%) | 174 / 2,068 | 0.078 | 424 (56.7%) | 105 / 319 | 0.248 | 436 | 0.40 / 0.24 |
| `STR` | 286 (8.0%) | 2 / 284 | 0.007 | 70 (9.4%) | 0 / 70 | 0.000 | 65 | 0.03 / 0.00 |
| `HP4-6` | 247 (7.0%) | 6 / 241 | 0.024 | 65 (8.7%) | 6 / 59 | 0.092 | 7 | 0.86 / 0.86 |
| `VNTR>6` | 28 (0.8%) | 0 / 28 | 0.000 | 6 (0.8%) | 0 / 6 | 0.000 | 0 | |
| `none` | 750 (21.1%) | 28 / 722 | 0.037 | 183 (24.5%) | 28 / 155 | 0.153 | 80 | 0.35 / 0.35 |
| total | 3,553 | 210 / 3,343 | 0.059 | 748 | 139 / 609 | 0.186 | 588 | 0.36 / 0.24 |

Plus 374 off-reference PASS calls (no GRCh38 allele; 324-1,376 across the runs).

- **The calls have the same shape as the truth, but more non-repeat and `HP4-6`.** Homopolymer `HP>=7` is 55-64% of the
  calls (truth 74%), `STR` 8-19% (truth 11%), `HP4-6` 5-7% (truth 1%), non-repeat 15-22% (truth 14%); `base` (R0.9 threshold
  on a weaker model) is 94% `HP>=7`.
- **Precision is poor in every class, worst in STRs.** Raw precision 0.06-0.08 for `HP>=7` and 0.03-0.04 for non-repeat; STRs
  are nearly all false (2-7 TP among 286-1,004 calls, precision <= 0.008), as are `HP4-6` and `VNTR>6`.
- **The PoN removes mostly homopolymer calls.** It cuts `HP>=7` calls by 81% (2,242 -> 424) and their TP by 40% (174 -> 105), so
  `HP>=7` precision rises to 0.22-0.25; STR calls fall by 75% and lose their last TPs; non-repeat calls lose no TP (28 -> 28,
  precision 0.04 -> 0.15). After the PoN the call set is 52-57% `HP>=7`, 9-16% STR, 24-25% non-repeat.
- **Recall by class** (raw / after PoN): `HP>=7` 0.39-0.41 / 0.23-0.24, `STR` 0.03-0.11 / 0.00, non-repeat 0.34-0.44 / unchanged.
  Long homopolymers are the worst: raw TP / truth 25 / 31 (7-9 bp), 75 / 137 (10-14), 55 / 128 (15-19), 18 / 112 (20-29),
  1 / 28 (>= 30), while FP calls grow with length (272, 659, 548, 518, 71). The STR and long-homopolymer truth is largely the
  part absorbed by the graph (`analysis/graph_absorbed_somatic_20261001`: 1,916 of the Illumina truth INDELs have no
  candidate because the reads align perfectly to a graph branch), so it is missing before the model sees anything.

Files: `repeat_context.py`, `tables.md` (truth tables; per run: class table, HP length bins), `calls_context.tsv`.

## Other somatic INDEL truth sets, same class rule (chr1-22, PASS alleles)

| truth set | INDELs | `HP>=7` | `HP4-6` | `STR` | non-repeat | 1-bp share | HP tract median |
|---|---:|---:|---:|---:|---:|---:|---:|
| GIAB HG008T v0.2 (PDAC cell line; HiFi + ONT + assembly) | 8,496 | 74.9% | 2.8% | 11.0% | 11.2% | 55% | 17 bp |
| SMaHT COLO829T union (ours; melanoma; multi-tech) | 1,912 | 61.0% | 5.1% | 14.4% | 19.3% | 65% | 15 bp |
| NYGC COLO829 v6 (Illumina tumor-normal, Strelka2 + Mutect2 ensemble) | 953 | 16.3% | 13.0% | 10.0% | 60.4% | 44% | 9 bp |
| SEQC2 HCC1395 v1.2.1 sINDEL, HC regions (breast cancer; Illumina multi-site; DeepSomatic's training truth) | 1,625 | 15.5% | 14.4% | 1.8% | 68.2% | 43% | 10 bp |

The two long-read / assembly-based sets (HG008T, SMaHT COLO829T) are 61-75% long-homopolymer INDELs; the two
short-read ensemble sets (NYGC COLO829, SEQC2 HCC1395) are 60-68% non-repeat with short homopolymers (median tract 9-10 bp),
because short-read callers and the consensus / high-confidence filtering remove long-homopolymer calls. The same tumor
(COLO829) gives 16% vs 61% `HP>=7` depending on the truth-set method. Computed with the classifier of `repeat_context.py`
on the local VCFs (`/scratch/jshen/data/SEQC2/`, `/scratch/jshen/data/COLO829T/`).
