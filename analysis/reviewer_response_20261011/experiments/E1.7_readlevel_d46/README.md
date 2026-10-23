# E1.7: read-level classes of the truth alleles without a candidate, HPRC v2.1 (d46) tensor sets

Reviewer point: **R2.3** (explain the false negatives; replaces the categories of the original Supplementary Fig. 2).
The method of `analysis/tensor_recall_20260930/colo829t_miss` (itself that of `analysis/hg008_somatic_miss_20260926`, HPRC
v1.1), run on the six current HPRC v2.1 sets (HG008T / COLO829T x Illumina / PacBio HiFi or Fiber-seq / ONT, chr1-22).

**Status (2026-10-22): done.** `classes.md` / `classes.tsv`, `supp_fig2.png` / `.pdf` (the revised Supplementary Fig. 2).

## Method

For every truth allele with status `no_candidate` in the set's `tensors/somatic.recall.tsv`, plus 250 controls that got a
tensor: the GAM records over the site are fetched and decoded, the reads that span it are compared with the REF and ALT
haplotypes, and the ALT reads' alignment decides the class (`scripts/report.py`; SNV classes S1-S4c, INDEL classes I1-I9;
definitions in `analysis/tensor_recall_20260930/README.md` section 2). `make_supp_fig2.py` merges them into five groups:

| group | SNV classes | INDEL classes |
|---|---|---|
| ALT already in the graph | S3 | I1 (fully), I2 (partly, residual edit with another spelling) |
| repeat: absorbed or written differently | S4a, S4b, S4c | I3 (reads on GRCh38, another spelling) |
| fewer than 3 ALT reads or low allele fraction | S2 | I4 |
| no read with MAPQ > 10 | S1 | I5 |
| other | - | I6 (> 50 bp), I7 (unclear), I9 |

Changes from the 2026-09-30 scripts (the classification is unchanged):
- the set's package, graph index and GRCh38 path come from its `run/config.json` (v4 / d46);
- target nodes from the single `discovery/target_nodes.txt` (the run's thousands of part files are slow on BeeGFS);
- the reader is the current repository's `indexed_gam_pipeline_v4` (`E17_SOURCE`): it reads the long-read GAMs through their
  per-record index (`<gam>.gri`); the runs' frozen 2026-10-03 copy predates it and did not finish one ONT locus in 70 min;
- `report.py`: the truth VAF fields (`VAF_Ill`, `VAF_PB`) exist only in the COLO829T VCF.

## Results (`classes.md`)

| set | INDELs without a candidate, ALT already in the graph | SNVs without a candidate: main groups |
|---|---|---|
| HG008T Illumina | 3,236 of 3,761 (86.0%) | repeat 52.2%, graph 20.3%, MAPQ 17.3% (n 295) |
| HG008T PacBio HiFi | 2,939 of 3,169 (92.7%) | repeat 64.0%, MAPQ 16.3%, graph 13.2% (n 386) |
| HG008T ONT | 2,590 of 2,800 (92.5%) | repeat 81.2%, graph 10.9% (n 329) |
| COLO829T Illumina | 523 of 643 (81.3%) | low ALT 38.6%, MAPQ 33.8%, repeat 20.0% (n 721) |
| COLO829T Fiber-seq | 385 of 463 (83.2%) | MAPQ 42.3%, repeat 30.8%, low ALT 18.8% (n 626) |
| COLO829T ONT | 332 of 390 (85.1%) | low ALT 70.2%, repeat 18.5% (n 523) |

The same picture as on HPRC v1.1 on 2026-09-30 (INDEL I1 + I2 81-92%). Controls: see each set's `controls.tsv` and
`summary.json` (`control:*`).

## Jobs

`tmp/reviewer_cpu_20261021/E1.7/` (`tasks.tsv`, `reads_array.sh`, `post_array.sh`, `jobs.txt`, logs): reads 450393
(tasks 1-3), 450438 (4-18), 450559 (11 and 18 again at 128 GB after OOM at 64 GB; the script resumes from the loci
already written); post 450563_5, 450573 (COLO829T ONT / Fiber-seq at 128 GB); reports 450587-450589.

## Files

| file | what |
|---|---|
| `scripts/reads.py`, `snv_nodes.py`, `snv_detail.py`, `report.py` (+ job scripts) | per-set analysis (SAMPLE = set) |
| `<set>/reads_part*.jsonl`, `snv_detail.jsonl`, `snv_nodes.json` | raw per-locus records |
| `<set>/indel_no_candidate.tsv`, `snv_no_candidate.tsv`, `controls.tsv`, `summary.json` | per-set classes |
| `make_supp_fig2.py`, `classes.md`, `classes.tsv`, `supp_fig2.png` / `.pdf` | the revised Supplementary Fig. 2 and its table |
