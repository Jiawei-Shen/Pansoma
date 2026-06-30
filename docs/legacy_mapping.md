# Legacy Mapping

This file records where old scripts were placed during the first migration pass. The copied files are preserved for reference; they have not all been refactored into clean package modules yet.

## Reference Pipeline Scripts

| Old file | New file |
| --- | --- |
| `utility/Pansoma_build_nodes_from_idx_5channels.py` | `scripts/build_node_json.py` |
| `parseGAM_find_unperfect_nodes_V6.py` | `scripts/find_unperfect_nodes.py` |
| `parseGAM_filtered_unperfect_nodes_Ver6CPP.py` | `scripts/build_dat_idx.py` |
| `find_GRCh38_position.py` | `scripts/build_grch38_path_json.py` |
| `utility/node_filter_LongRead.py` | `scripts/filter_node_json.py` |
| `node_pileup_tensors/node_5channel_batchP_pileup_testing_tensor.py` | `scripts/generate_testing_tensors.py` |
| `node_pileup_tensors/node_6channel_batchP_pileup_tensor_Ver7_label.py` | `scripts/label_tensors.py` |
| `node_pileup_tensors/tensor_visualization_V8.py` | `scripts/visualize_tensor.py` |
| `node_pileup_tensors/node_tensor_classify_true_false.py` | `scripts/classify_tensors.py` |
| cluster-only chromosome component script | `scripts/build_chr_node_filters.sh` |

## Legacy Archives

| Old area | New area |
| --- | --- |
| root-level `*.py` processing scripts | `experiments/legacy/processing/` |
| `utility/` | `experiments/legacy/utilities/utility/` |
| `node_pileup_tensors/` | `experiments/legacy/tensors/node_pileup_tensors/` |
| `Figure_drawing/` | `experiments/legacy/figures/Figure_drawing/` |
| `test_new_dataformat/` | `experiments/legacy/test_new_dataformat/test_new_dataformat/` |
| `CPP_NPU_format/` | `cpp/legacy/CPP_NPU_format/` |

## Package Seeds

Several files were copied into `src/pangenome_ml_data_generation/` as starting points for future refactors. These files may still contain command-line parsing, old names, or direct path assumptions. Treat them as seeded code, not finished library modules.
