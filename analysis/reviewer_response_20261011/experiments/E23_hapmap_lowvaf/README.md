# E23: HapMap mixture on the donor-free graph, low-VAF detection limit

Reviewer points: **R2.4** (VAF-decreasing sensitivity, the low-VAF model), **R1.5** / **R2.1** (HapMap as a benchmark,
leakage), **R3.10** (low-VAF claims). User request 2026-10-22: rerun HapMap with the current pipeline and use it to
answer these points, with a reasoned account of low-VAF detection.

## Question

How far down in allele fraction do candidates, tensors and calls reach when the reads are deep enough, and which step
sets the limit? HapMap is the only data set we have with a dense truth below the 0.08 candidate floor of the default
build (90 % of its truth alleles have VAF < 0.08; COLO829T dilutions reach VAF ~0.05 at the lowest purity).

What HapMap cannot measure (R1.5, R2.1): its positives are germline alleles of five minor donors, not somatic
mutations. Many are common population alleles: the graph absorbs those it carries, and population filtering tags
96.3 % of its chr20 SNV positives (E8b task 6). E23 therefore reports the detection curve **without** population
filtering as a read-evidence limit, and with it only to show why the mixture is not a tumor-only somatic benchmark.

## Data

- Reads: the 8 HapMap libraries (LIB027514-17 + DIL re-runs; 10.5e9 Illumina reads, ~510x), Giraffe GAMs on
  `mgi_non_HG005.d9` (HPRC year-1 minus the six HapMap donors and NA19240: 40 individuals + GRCh38 + CHM13, 82 haplotypes,
  d9 filter; the graph of the submitted HapMap results). Sorted per library and read as one (multi-GAM, f3c159d).
- Tensors: `pansoma_v2_tensors/HapMap_Illumina_nonHG005` (code f3c159d; discovery `--node-alt 0.01`; builder SNV / INDEL AF
  floor **0.015**, `--min-variants 3`, `--haplotypes 82`, autosomes): 43,201,211 tensors (SNV 42,415,603, INDEL 785,608).
- Truth: `pansoma_v2_tensors/HapMap_truth/` (README there): somatic = SMHTHAPMAP6 v1.6 without HG005 (7.87 M; FORMAT AF =
  expected VAF, 0.0025-0.165); germline = alleles of the background HG005; BED = v1.6 autosome benchmark BED.
- Labels: `tmp/hapmap_nonHG005_20261012/label_parallel.py` (job 450540), the run's frozen rules, one process per
  chromosome (the serial finalize would have exceeded its time limit); checked byte for byte against the serial
  run's finished chromosomes (SNV chr1, chr10-13; `tensors/SNV/.serial_label_tmp_20261022`).
- Models (tequila H100s, `tmp/hapmap_nonHG005_20261012/train/`): SNV recipe with 4 M sampled tensors per epoch
  (`HapMap_Illumina_SNV_sample4M`, job 450542), INDEL recipe (`HapMap_Illumina_INDEL`, 450543), then the SNV recipe on
  all training tensors on two GPUs (`HapMap_Illumina_SNV_full`, 450544); chr2-22 train, chr1 test; `vcf.sh` afterwards.

## Analyses

1. `tensor_vaf.py` (after the labels): per truth VAF bin, the share of truth alleles with a tensor, filtered at the AF
   floor, filtered for < 3 ALT reads, or without any candidate (absorbed by the graph, or no imperfect read).
2. Call level (after the models): recall per truth VAF bin and precision per call AF bin on chr1, without and with
   population filtering (rtg vcfeval; `vcf.sh` gives the PoN arm, a second `vcfeval` without `--pon` the other).

## Results (2026-10-26)

Files: `tables.md` / `tensor_vaf.tsv` (tensor level) and `call_vaf.<run>.md` (call level) use the exact graph alleles
(`graph_alleles.sh`); `atomized/` repeats them with multi-base graph alleles of <= 50 bp split into SNVs and indels
(`ATOMIZE=1`), the definition the response uses. `donors_in_graph.chr1.tsv` (`donors_in_graph.py`): in-graph share by
the number of minor donors carrying an allele. `neighbours.txt` (`neighbours.sh`): truth positions with another one
within 50 bp.

1. **The graph absorbs shared alleles.** On chr1 (BED), SNVs carried by 1 / 2 / 3 / 4 / 5 minor donors are paths of
   the graph in 23.8 / 65.4 / 90.5 / 97.4 / 99.0 % of cases (atomized); about 64 % of the in-graph SNVs are shared
   by two or more donors. In-graph SNVs have a tensor in about 1 % of cases and are never called (recall <= 0.001),
   in-graph INDELs in about 2 %. Atomizing matters in the bins of shared alleles only: the SNVs it adds have no tensor
   either (the not-in-graph tensor share at 0.08-0.095 rises from 0.43 to 0.62, at 0.03-0.045 from 0.68 to 0.80).
2. **Tensor level, SNVs not carried by the graph** (autosomes, atomized): with a tensor 0.032 / 0.084 / 0.336 / 0.845
   / 0.933 at VAF 0.0025 / 0.005 / 0.01 / 0.02-0.025 / 0.05-0.075. At 0.01, 45.6 % are filtered (AF floor or < 3 ALT
   reads) and 20.8 % have no candidate (the `--node-alt 0.01` discovery threshold equals the VAF). The curve is flatter
   than the Poisson prediction below: 0.34 against 0.13 at VAF 0.01, 0.80 against 0.98 at 0.03-0.045, although the
   median depth at the true chr1 calls is 531 (FORMAT DP, p10 328, p90 604) [inference: ALT counts more dispersed than
   Poisson, e.g. PCR duplicates (not removed) and uneven mixing].
3. **Call level, chr1, BED, before population filtering, truth not carried by the graph** (atomized):

   | truth VAF | SNV full | SNV sampled | INDEL |
   |---|---|---|---|
   | 0.0025 | 0.006 | 0.007 | 0.047 |
   | 0.005 | 0.033 | 0.032 | 0.059 |
   | 0.01 | 0.295 | 0.299 | 0.116 |
   | 0.0125-0.015 | 0.525 | 0.520 | 0.108 |
   | 0.02-0.025 | 0.809 | 0.801 | 0.189 |
   | 0.03-0.045 | 0.751 | 0.743 | 0.127 |
   | 0.05-0.075 | 0.905 | 0.901 | 0.365 |
   | 0.08-0.095 | 0.497 | 0.497 | 0.080 |
   | 0.10-0.165 | 0.791 | 0.769 | 0.175 |

   Whole call sets (chr1, BED, PASS = validation recall-0.9 threshold): SNV full P 0.990, R 0.277 (ceiling 0.308),
   sampled P 0.990, R 0.277: the full model adds nothing. INDEL P 0.772, R 0.155 (ceiling 0.176). SNV precision by call
   AF is 0.981-0.994 in every bin; INDEL 0.62-0.88. Where SNVs reach the classifier (VAF 0.01, 0.02-0.025, 0.05-0.075),
   recall is 85-95 % of the tensor share. After population filtering, SNV recall is 0.007 and INDEL 0.021: the PoN
   tags 95.7 % of the chr1 SNV truth records (411,638 of 430,009).
4. **Precision is likely optimistic.** 35.2 % of the HapMap chr1 truth positions have another truth position within
   50 bp, against 4.4 % for COLO829T: reads that carry one donor allele often carry others [inference: this makes a
   true allele easier to tell from an error than for an isolated somatic mutation; not tested].
5. **INDELs** stay low even without graph absorption: 45 % of the not-in-graph INDELs at 0.05-0.075 have no
   candidate (chr1, atomized) [not investigated; probably repeat indels represented or aligned differently, as on
   HG008T / COLO829T (R1.2)].

## Expected limits (reasoning before the results)

Read sampling alone sets a floor. With a candidate floor `a` and `--min-variants 3`, a site of depth D needs
max(3, ceil(a D)) ALT reads. At D = 510 and a = 0.015 that is 8 reads; for Poisson read sampling the chance of a
candidate is 0.005 at VAF 0.005, 0.14 at 0.01, 0.50 at 0.015, 0.80 at 0.02 and 0.985 at 0.03. With the default floor
0.08 it is 0.51 at VAF 0.08 and below 0.003 at 0.05 for any depth >= 100 (the 0.08 floor, not the depth, is the limit
there). The table of `tensor_vaf.py` should follow this curve where the graph does not absorb the allele.

Further limits, by step:
- **Discovery** (`--node-alt`, 0.05 by default, 0.01 here): a node is examined only if more than this share of its
  mappings carry an edit. Sequencing errors add to the share, so it rarely binds on short reads.
- **Tensor rows**: 200 rows, sampled uniformly over the allele blocks, so an allele at VAF f fills about min(D, 200) x f
  rows (3 rows at VAF 0.015, 16 at 0.08). Beyond 200x the extra ALT reads count toward the AF filters but are not shown
  to the network. A low-VAF configuration could keep every ALT row (not done).
- **Classifier**: the default models never saw tensors below AF 0.08; the HapMap models are trained at 0.015.
- **Sequencing errors**: at VAF 0.0025-0.005 (1-3 ALT reads at 510x) a true allele is not separable from errors without
  UMIs or duplex reads; no tumor-only caller resolves these from one sample.
- **Graph absorption** (HapMap-specific): alleles carried by >= 9 of the 82 haplotypes are graph paths; reads that carry
  them align perfectly and give no candidate. For real somatic truth this matters little: on HG008T the absorbed INDEL
  truth that population filtering would not tag anyway is 1.3 % / 0.6 % of the INDEL truth
  (`tmp/graph_absorbed_pon_20261012`).
- **Population filtering** (HapMap-specific): the positives are population alleles (96.3 % of chr20 SNVs tagged).

## Status

- 2026-10-22 16:26: labels 450540 running (maotai), chunk cache 450541 (yuengling, ~5 s per shard); models
  450542 (SNV, 4 M per epoch) and 450543 (INDEL) wait for both, then 450544 (SNV, all tensors, two H100s); `vcf.sh`
  450545-450547.
- 2026-10-22 16:30: E8b truth-in-resources for the HapMap chr20 positives rerun with the 10-16 INDEL rule (job 450562;
  `somatic_truth_in_resources.HapMap_positives_chr20.pon_hp8_1e-4.tsv`), for R2.1 point 3.
- 2026-10-23: labels 450540 done 10-22 23:54 (7.5 h; MaxRSS 292 GiB). Check against the serial run: the 354 finished
  serial files (SNV chr1, chr10-13: labels.ndjson and every shard's labels.npy) are byte-identical; its interrupted
  chr14_labels.ndjson.tmp (957,238 complete lines) equals the first 957,238 lines of the new file.
- Response drafts with [待补: E23 ...] placeholders: R2.4 ("Below the candidate floor", the low-VAF model, usage),
  R1.5 (Supplement role), R2.1 (Supplement role, filtering share), R3.10 (Low VAF).
- 2026-10-26: models done (tequila H100: SNV sampled 1 GPU 16.9 h, INDEL 1 GPU 31 min, SNV full 2 GPUs 36.1 h); `vcf.sh`
  of SNV full rerun at 32G (452207; 450547 OOM at 8G). Atomized graph alleles: `norm -a` on everything exceeded 64 GB
  (452208, on the 1.6 Mb REF / 1.28 Mb ALT allele at chr5:69,765,910, SMN region), so only alleles <= 50 bp are
  atomized (452298, 5 min, 1.1 GB); tables 452299, SNV full call table 452300. Numbers filled into R2.4 ("Below the
  candidate floor") and R3.10 (Low VAF).
