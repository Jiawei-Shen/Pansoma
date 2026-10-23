# E20: do the matched normal's reads take the alternative graph path at the COLO829T absorbed truth alleles?

Reviewer point: **R2.3** (Reviewer 2, major 3): "for false negatives at existing graph nodes, check whether the normal
sample's reads preferentially align to the alternative path". User decision 2026-10-21: COLO829BL only (its ONT reads
are local; HG008-N reads would have to be downloaded).

**Status (2026-10-22): done.** Classify jobs 450335 / 450437 (d9 normal / tumor) and 450435 / 450436 (d46) finished; `report.py` run on the login node on 2026-10-22 (`tables.md`, `summary.json`). The d46 loci are E5d's 401 absorbed alleles.

## Question

At a truth allele the graph absorbs (the tumor's ALT reads walk a graph branch or skip edge that already spells the ALT,
so no edit and no candidate), where do the normal's reads go? If the graph pulls reads onto the alternative path, normal
reads whose own bases are GRCh38 (REF) would also walk the ALT path; if not, the normal's REF reads stay on GRCh38 and only
reads that carry the allele take the branch.

## Data

| what | where |
|---|---|
| loci, HPRC v1.1 (d9) | the 497 absorbed alleles (67 SNVs, 430 INDELs) of `analysis/graph_absorbed_somatic_colo829t_20261005` (`variant_set.tsv`), with its categories (`hprc_per_variant.tsv`: germline filtering 34 / pangenome-induced false negative 430 / ambiguous 33) and patient-frame classes (`normal_frame.tsv`) |
| loci, HPRC v2.1 (d46) | the 401 absorbed alleles (54 SNVs, 347 INDELs) of E5d `colo829t_d46_absorbed_loci.tsv`: E5b's 'absorbed' rule, the set R1.4 / R2.3 quote (no candidate on d46 and a tensor on the GRCh38-only graph on >= 1 platform, and a d46 path spells GRCh38 + truth); 325 of them are also in the d9 set. (The regions were extracted for a wider set, the 599 E5c I1 loci and every E5c COLO829T lost locus; all 401 lie inside them.) |
| normal reads | COLO829BL ONT, SMaHT, GRCh38 minimap2 BAMs `SMHTCOLO829BL-X-X-M45-D001-...` and `...-D003-...` (`/scratch/jshen/data/COLO829T/COLO829BL_ONT/`, read only) |
| tumor reads | COLO829T ONT `COLO829T_ONT_merged.bam` (the reads of Pansoma's COLO829T ONT GAMs), and Pansoma's GAMs themselves: d9 `COLO829T_ONT_merged.sorted.gam`, d46 `vg_hprc_v2.1_d46/COLO829T_ONT_hprc-v2.1-d46_vg.sorted.gam` |
| graphs | d9 `AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.*`, d46 `AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.*`; graph index / GRCh38 path as the tensor sets |

## Method

1. `loci.py`: loci tables and `regions.bed` (every d9 / d46 locus and every E5c COLO829T lost locus +- 500 bp, merged:
   1,371 loci, 1,330 regions, 1.34 Mb).
2. `extract.sh`: primary records (`-F 0x900`) overlapping the regions, each once (`samtools view -M -L`), to FASTQ, from
   the two normal BAMs and the tumor BAM. Normal = D001 + D003.
3. `giraffe.sh`: the FASTQs to each graph with the giraffe settings of the COLO829T ONT GAMs (d9: vg 1.66, `-b r10`,
   d9 longread minimizers / zipcodes / dist; d46: vg 1.77, the same with the d46 indexes), `vg gamsort -i` (vg 1.77).
   The tumor reads are realigned too, so normal and tumor go through the same selection and alignment; Pansoma's own
   tumor GAMs are classified as well as a check that the realigned subset aligns as in the whole-genome run.
4. `classify.py`: per locus, every record with MAPQ > 10 that spans the allele (matched GRCh38 bases within 2 kb on both
   sides; E5c's read selection and decoding with the frozen pipeline decoder) gets two labels: the **sequence** label
   (its own bases between the anchors: REF-like / ALT-like / unclear, by edit distance to GRCh38 and GRCh38 + truth) and
   the **path** label (the graph bases it walks there: GRCh38 nodes with consecutive positions; a path spelling the ALT;
   one closer to the ALT than to GRCh38; another path).
5. `report.py`: pooled and per-locus fractions, normal against tumor; per d9 category.

## Results

See `tables.md` for every number. Main points (normal = COLO829BL ONT, tumor = COLO829T ONT, both realigned the same way):

- The graph does not pull reads onto the alternative path: reads whose own bases are REF-like walk the ALT path in 1.1% (normal, both graphs) and 1.6-1.8% (tumor) of cases, while 94-97% of ALT-like reads do.
- SNVs: 4.0% (d9) / 2.2% (d46) of normal reads carry the ALT, against 45% of tumor reads.
- INDELs (mostly homopolymer / STR length changes): 41-42% of normal and 58-59% of tumor reads are ALT-like; ONT length errors in such repeats make single-read REF / ALT labels unreliable, so germline status there rests on the assemblies.
- Realigned tumor reads agree with Pansoma's own GAMs (per-locus ALT fraction within 0.1 at 249 / 250 (d9) and 256 / 258 (d46) loci).

## Caveats

- Reads are selected by their GRCh38 (minimap2) placement; a read that minimap2 places elsewhere but giraffe would place
  at the locus is missed, for normal and tumor alike.
- ONT reads carry homopolymer-length errors, so at homopolymer loci the sequence label of single reads is noisy; the
  comparison that matters is normal against tumor through the same pipeline.

## Files

| file | what |
|---|---|
| `loci.py`, `loci_summary.json` | loci and regions |
| `classify.py` | per-locus read labels |
| `report.py`, `tables.md`, `summary.json` | tables |
| `/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261021/E20/` | `extract.sh`, `giraffe.sh`, `classify.sh`, `jobs.txt`, logs |
| `/scratch/jshen/data/reviewer_experiments_20261013/E20_normal_reads_colo829bl/` | `loci_*.tsv`, `regions.bed`, `fastq/`, `gam/`, `classify/*.jsonl`, `per_locus.tsv` (not in git) |
