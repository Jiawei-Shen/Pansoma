# E1.7: how often a tensor window reaches beyond the candidate's own node

Reviewer points: **R1.1** (path-following tensors), **R1.m2** (tensor windows truncated at node boundaries).

For every chr1 tensor of a set, `window_span.py` reads the site's node, start and end from
`tensors/<KIND>/chr1_variant_summary.ndjson`, and takes the node's length from the set's graph-index sqlite.
- The 101-column window covers 50 columns on each side of the site.
- It stays inside the candidate's node only if the node has at least 50 bases on both sides of the site.
- Somatic means label 1 in `chr1_labels.ndjson`.
- The INDEL rule uses the same ±50 bases around the event on the node, an approximation, since INDEL windows also hold
  insertion slots.

The script is read-only and was run on the login node on 2026-10-16:
`python window_span.py HG008T_PacBio HG008T_PacBio_d46 HG008T_PacBio_linear COLO829T_fiberseq_d46 > window_span.tsv`.

## Results (`window_span.tsv`, chr1, all tensors and the somatic ones)

| set | graph | kind | tensors: window beyond the node | somatic tensors: window beyond the node | median node length (bp) |
|---|---|---|---:|---:|---:|
| HG008T PacBio | HPRC v1.1 (d9) | SNV | 51,560 / 85,133 (60.6 %) | 399 / 703 (56.8 %) | 95 |
| HG008T PacBio | HPRC v1.1 (d9) | INDEL | 123,144 / 161,836 (76.1 %) | 465 / 507 (91.7 %) | 110 |
| HG008T PacBio | HPRC v2.1 (d46) | SNV | 56,304 / 140,610 (40.0 %) | 533 / 706 (75.5 %) | 1,024 |
| HG008T PacBio | HPRC v2.1 (d46) | INDEL | 143,857 / 167,382 (85.9 %) | 496 / 515 (96.3 %) | 45 |
| COLO829T Fiber-seq | HPRC v2.1 (d46) | SNV | 71,828 / 241,213 (29.8 %) | 1,456 / 2,025 (71.9 %) | 1,024 |
| COLO829T Fiber-seq | HPRC v2.1 (d46) | INDEL | 139,541 / 164,198 (85.0 %) | 102 / 105 (97.1 %) | 34 |
| HG008T PacBio | GRCh38-only | SNV | 302,200 / 302,200 (100 %) | 722 / 722 (100 %) | 32 |
| HG008T PacBio | GRCh38-only | INDEL | 203,810 / 203,810 (100 %) | 591 / 591 (100 %) | 32 |

## Reading

- On the HPRC graphs the window of 57–76 % of the somatic SNV tensors and 92–97 % of the somatic INDEL tensors reaches
  beyond the candidate's node.
- A tensor confined to its node, as in the first version, would have truncated the context of most somatic candidates.
- INDEL candidates sit on short nodes (median 34–110 bp), since they are mostly in homopolymer / STR bubbles.
- Many SNV candidates sit on long reference nodes: d46 nodes are cut at 1,024 bp, the median here.
- Somatic SNVs lie near node boundaries more often than all SNV tensors do (d46: 72–76 % against 30–40 %).
- The GRCh38-only graph (vg autoindex) has 32-bp nodes, so every window crosses nodes there.
