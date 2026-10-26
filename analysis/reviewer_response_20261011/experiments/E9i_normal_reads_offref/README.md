# E9i: matched-normal read evidence for COLO829T graph-coordinate SNV calls (R3.3)

Reviewer 3, major 3, asks for orthogonal validation of the graph-coordinate calls, "matched-normal alignment evidence"
first. E9c / E9d check the calls against assemblies; this experiment checks them against the matched normal's reads.

**Status (2026-10-23): done.** Array 451636 (60 shards, `tmp/reviewer_cpu_20261023/E9i/`; 451619 cancelled at once:
`-w` with four nodes allocated all four to every task), `report.py` -> `tables.md`, `summary.json`.

## Method

- Calls: E9c's `calls.tsv` of the six COLO829T in-sample SNV runs (HPRC v1.1 d9: Illumina at epoch 6, Fiber-seq, ONT;
  HPRC v2.1 d46: Illumina, Fiber-seq, ONT): the unplaced (graph-coordinate) PASS calls (COLO829T ONT: E9c's random
  3,000) and the placed controls `ctrl_tp` (true somatic after the PoN, 150) and `ctrl_fp` (false positives after the
  PoN, 150).
- Reads: COLO829BL ONT (SMaHT D001 + D003, minimap2 to GRCh38) = matched normal; COLO829T ONT
  (`COLO829T_ONT_merged.bam`) = tumor check. Primary records overlapping the call's GRCh38 locus (E9c `locus_*`: the
  GRCh38 interval of the call's bubble) +- 500 bp. No COLO829BL Illumina or Fiber-seq BAM on GRCh38 is local (Illumina:
  FASTQ only; Fiber-seq: aligned to the COLO829BL assembly); HG008-N reads are not local (HG008T uses the HG008-N
  assembly, E9d).
- Test (`normal_reads.py`): the 21-mer centred on the call in its ALT context (E9c `alt_read`, the consensus of the
  call's ALT reads; `alt_graph` when missing) and in its REF context (`ref_graph`), counted in each read on both strands.
  'ALT in normal' = >= 2 normal reads with the ALT 21-mer; 'absent from normal' = 0 such reads and >= 5 with the REF
  21-mer; 'somatic candidate' = absent from the normal and >= 2 tumor reads with the ALT 21-mer.

## Results (`tables.md`)

| group | ALT in normal | absent from normal | somatic candidate |
|---|---|---|---|
| placed true somatic (ctrl_tp), 6 x 150 | 1.3-4.0 % | 88.7-94.0 % | 88.7-94.0 % |
| placed false positives after the PoN (ctrl_fp), 6 x 150 | 58.7-87.3 % | 10.0-39.3 % | 8.0-25.3 % |
| graph-coordinate PASS calls (257-3,000 per run) | 66.0-92.8 % | 2.2-16.0 % | 0.7-3.1 % |

- The assay separates the controls: true somatic calls lack the ALT in the normal's reads (89-94 %), and false positives
  carry it in most runs (59-87 %; the rest are mostly ONT errors without tumor support either).
- Most graph-coordinate calls are in the normal's reads (66-93 % per run); 0.7-3.1 % per run are somatic candidates.
- Agreement with the assemblies (E9c, COLO829T with the hifiasm tumor assembly and the COLO829BL DSA): of the 64
  assembly tumor-only calls, 54 lack the ALT in the normal's reads and 50 are somatic candidates; of the 5,578
  population alleles, 88 % carry the ALT in >= 2 normal reads.

## Caveats

- ONT reads and exact 21-mers: a read with an error inside the 21-mer is missed, so 'one normal read' and the low
  'ALT in normal' share of some false-positive controls (ONT d46 58.7 %) include errors.
- Repeats: an ALT 21-mer that also occurs elsewhere in the locus window is counted; the controls bound this
  (1.3-4.0 % of true somatic calls show the ALT in the normal).
- COLO829T only; ONT normal reads only.
