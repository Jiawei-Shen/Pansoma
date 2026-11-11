# Tensors of the two R1.m2 figures

| Tensor | Figure | Source (read-only) |
|---|---|---|
| `R1.m2_tensor_across_nodes.npy` | `../fig_R1.m2_tensor_across_nodes.png` | `HG008T_PacBio/tensors/SNV/chr1_shard_00000_data.npy`, index 32257: candidate 3130240:147:SNP:G>A, somatic truth SNV chr1:242,529,453 G>A (HPRC v1.1, d9) |
| `R1.m2_INDEL_span_01.npy` | `../indel/fig_R1.m2_INDEL_span_01.png` | `HG008T_Illumina_d46/tensors/INDEL/chr1_shard_00000_data.npy`, index 271: candidate 1811743:0:DEL, somatic truth deletion chr1:23,017,177 TTATATATA>T (HPRC v2.1, d46) |

Sources are under `/scratch/jshen/data/pansoma_v2_tensors/`. Each `.npy` is a byte copy of the shard row, shape
(1, 8, 200, 101) int8 (channels x read rows x window columns; channel order as in the pipeline's tensor format).
`<name>.summary.ndjson` is the tensor's `chr1_variant_summary.ndjson` record, with `shard_index` / `index_within_shard`
set to 0 and the source values in `source_*`. `<name>.label.json` holds its `chr1_labels.ndjson` record (label 1,
`representative_allele_is_somatic_truth`) and its value in `chr1_shard_00000_labels.npy`.

`extract.py` writes these files. `render.sh [outdir]` redraws both figures from them with `scripts/visualize_tensor.py`;
the output is pixel-identical to the committed PNGs (checked 2026-11-11).
