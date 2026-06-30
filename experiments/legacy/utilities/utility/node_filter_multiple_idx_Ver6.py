#!/usr/bin/env python3
import json
import struct
import argparse
import sys
import time
import copy
import re
import os
import pysam  # for VCF queries via pysam.VariantFile
from typing import Dict, Set, List, Optional

# ─────────────────────────────────────────────────────────────────────────────
# Latest format constants ONLY (kept for reference)
GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")   # major, minor, block_count, reserved[16]

# Latest per-block header in .dat (18B): nid, nrec, flags, R, C
BLOCK_HDR_PACK = struct.Struct("<I I H I I")   # 18 bytes

# Latest .idx entry (30B): nid, offset, block_size, n_records, flags, R, C
IDX_ENTRY_PACK = struct.Struct("<I Q I I H I I")
IDX_ENTRY_SIZE = IDX_ENTRY_PACK.size  # 30

# ─────────────────────────────────────────────────────────────────────────────
# Strict latest-only IDX loader
# ─────────────────────────────────────────────────────────────────────────────
def load_index_latest(idx_path: str) -> Dict[int, dict]:
    """
    Strictly parse latest .idx format (30 bytes/entry):
      <I Q I I H I I> → nid, offset, block_size, n_records, flags, R, C

    Returns dict[nid] = {
        "start": offset,
        "block_size": block_size,
        "n_records": n_records,
        "flags": flags,
        "max_read_len": R,
        "max_cigar_len": C,
    }
    """
    node_index: Dict[int, dict] = {}

    try:
        with open(idx_path, "rb") as f:
            hdr = f.read(4)
            if len(hdr) != 4:
                raise RuntimeError("idx header too short")
            (count,) = struct.unpack("<I", hdr)

            # file size must match exactly: 4 + count*30
            f.seek(0, os.SEEK_END)
            size = f.tell()
            expected = 4 + count * IDX_ENTRY_SIZE
            if size != expected:
                raise RuntimeError(
                    f"idx size mismatch for latest format: file={size}, expected={expected} "
                    f"(count={count}, entry={IDX_ENTRY_SIZE})"
                )

            f.seek(4)
            for i in range(count):
                rec = f.read(IDX_ENTRY_SIZE)
                if len(rec) != IDX_ENTRY_SIZE:
                    raise RuntimeError(f"truncated idx at entry {i+1}")
                nid, off, blk_sz, nrec, flg, R, C = IDX_ENTRY_PACK.unpack(rec)
                node_index[nid] = {
                    "start": off,
                    "block_size": blk_sz,
                    "n_records": nrec,
                    "flags": flg,
                    "max_read_len": int(R),
                    "max_cigar_len": int(C),
                }
    except FileNotFoundError:
        print(f"Error: idx not found: {idx_path}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"Error parsing latest idx {idx_path}: {e}", file=sys.stderr)
        return {}

    return node_index

# ─────────────────────────────────────────────────────────────────────────────
# JSON helpers
# ─────────────────────────────────────────────────────────────────────────────
def load_json_data_ids_and_map(json_filepath: str):
    main_json_data = {}
    node_ids_set = set()
    json_nodes_map = {}
    try:
        with open(json_filepath, 'r') as f_json:
            main_json_data = json.load(f_json)
            if not isinstance(main_json_data, dict):
                print(f"Error: Root JSON content in {json_filepath} is not an object/dictionary.", file=sys.stderr)
                return None, set(), {}
            node_list_from_json = main_json_data.get("nodes")
            if not isinstance(node_list_from_json, list):
                print(f"Warning: 'nodes' key not found or not a list in {json_filepath}.", file=sys.stderr)
                return main_json_data, set(), {}
            for item_index, item in enumerate(node_list_from_json):
                if not isinstance(item, dict):
                    continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None:
                        continue
                    node_id_int = int(node_id_str)
                    node_ids_set.add(node_id_int)
                    json_nodes_map[node_id_int] = item
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Node ID format error for item {item_index} in 'nodes' (JSON): '{node_id_str}'. Error: {e}. Skipping.",
                        file=sys.stderr)
            return main_json_data, node_ids_set, json_nodes_map
    except FileNotFoundError:
        print(f"Error: JSON file not found: {json_filepath}", file=sys.stderr)
        return None, set(), {}
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_filepath}.", file=sys.stderr)
        return None, set(), {}
    except Exception as e:
        print(f"Unexpected error reading JSON {json_filepath}: {e}", file=sys.stderr)
        return None, set(), {}

def extract_chromosome_from_path_pattern(path_pattern: Optional[str]) -> Optional[str]:
    if not path_pattern or not isinstance(path_pattern, str):
        return None
    match = re.search(r'(chr([0-9A-Za-z_]+|X|Y|M|MT))', path_pattern)
    if match:
        return match.group(1)
    return None

# ─────────────────────────────────────────────────────────────────────────────
# Prefix file expansion
# ─────────────────────────────────────────────────────────────────────────────
def read_idx_paths_from_prefix_file(prefix_file: str) -> List[str]:
    """
    Read lines from prefix file; for each non-empty, non-comment line:
      - strip whitespace
      - if it ends with '.idx' (case-sensitive), use as-is
      - else append '.idx'
    """
    paths: List[str] = []
    try:
        with open(prefix_file, 'r') as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith('#'):
                    continue
                if line.endswith('.idx'):
                    paths.append(line)
                else:
                    paths.append(line + '.idx')
    except FileNotFoundError:
        print(f"Error: --npu_paths_file not found: {prefix_file}", file=sys.stderr)
    except Exception as e:
        print(f"Error reading --npu_paths_file {prefix_file}: {e}", file=sys.stderr)
    return paths

# ─────────────────────────────────────────────────────────────────────────────
# Main filter (union across multiple IDX files)
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath: str,
                                idx_filepaths: List[str],
                                output_json_filepath: str,
                                vcf_file: Optional[str] = None,
                                txt_output_path: Optional[str] = None):
    print("Step 1: Loading JSON data and node map...")
    main_json_structure, target_node_ids_from_json, json_nodes_map = load_json_data_ids_and_map(json_filepath)
    if main_json_structure is None:
        print("Critical error: Could not load main JSON. Cannot proceed.")
        return
    print(f"Loaded {len(json_nodes_map)} nodes from JSON for potential processing.")

    chromosome_for_vcf = None
    if vcf_file:
        print(f"Using pysam for VCF queries on: {vcf_file}")
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_vcf = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_vcf:
            print(
                f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'. Using 'chr1' as fallback.",
                file=sys.stderr)
            chromosome_for_vcf = "chr1"
        else:
            print(f"Extracted chromosome '{chromosome_for_vcf}' for VCF queries.")

    print("\nStep 2: Loading node indices (LATEST IDX only) from multiple files...")
    union_idx_node_ids: Set[int] = set()
    total_found_files = 0
    for idx_path in idx_filepaths:
        if not os.path.exists(idx_path):
            print(f"Warning: idx path does not exist, skipping: {idx_path}", file=sys.stderr)
            continue
        idx_data = load_index_latest(idx_path)
        if idx_data:
            total_found_files += 1
        before_union = len(union_idx_node_ids)
        union_idx_node_ids |= set(idx_data.keys())
        after_union = len(union_idx_node_ids)
        print(f"  Loaded {len(idx_data)} nodes from: {idx_path}")
        print(f"    Union grew: {before_union} → {after_union}")

    print(f"Total unique node IDs across {total_found_files} readable idx file(s): {len(union_idx_node_ids)}")
    if not union_idx_node_ids:
        print("Warning: No nodes were loaded from any idx. Output will be empty unless VCF-only filtering keeps none anyway.", file=sys.stderr)

    print("\nStep 3: Identifying node IDs present in JSON and in ANY idx (union)...")
    common_node_ids_int = target_node_ids_from_json.intersection(union_idx_node_ids)
    num_initial_common_nodes = len(common_node_ids_int)
    if num_initial_common_nodes == 0:
        print("No common node IDs found between JSON and the union of IDX files. Output 'nodes' list will be empty.")
    else:
        print(f"Found {num_initial_common_nodes} common node IDs between JSON and union of IDXs for potential processing.")

    ultimate_filtered_nodes_list = []
    total_nodes_queried_with_pysam = 0
    total_nodes_filtered_out_by_vcf_prereq = 0
    total_nodes_filtered_out_by_empty_vcf = 0
    total_nodes_kept_with_vcf_results = 0

    if num_initial_common_nodes > 0:
        print(
            f"\nStep 4: Processing {num_initial_common_nodes} common nodes for filtering and VCF query (if applicable)...")

        try:
            with pysam.VariantFile(vcf_file) if vcf_file else open(os.devnull, 'w') as vcf:
                sorted_common_ids = sorted(list(common_node_ids_int))

                batch_node_counter_step4 = 0
                batch_queried_for_vcf_step4 = 0
                batch_kept_with_vcf_results_step4 = 0
                batch_start_time_step4 = time.time()

                for node_idx_step4, node_id_int in enumerate(sorted_common_ids):
                    original_node_item = json_nodes_map.get(node_id_int)
                    if not original_node_item:
                        sys.stderr.write(f"Internal Warning: Common node ID {node_id_int} not in JSON map. Skipping.\n")
                        continue

                    current_node_copy = copy.deepcopy(original_node_item)
                    include_node_in_final_output = True

                    if vcf_file and chromosome_for_vcf:
                        try:
                            start_0based_str = current_node_copy.get("grch38_position_start")
                            length_str = current_node_copy.get("length")

                            if start_0based_str is None or length_str is None:
                                include_node_in_final_output = False
                                total_nodes_filtered_out_by_vcf_prereq += 1
                            else:
                                start_0based = int(start_0based_str)
                                length = int(length_str)
                                if length <= 0:
                                    include_node_in_final_output = False
                                    total_nodes_filtered_out_by_vcf_prereq += 1
                                else:
                                    query_end_0based = start_0based + length
                                    vcf_records_iterator = vcf.fetch(chromosome_for_vcf, start_0based, query_end_0based)
                                    vcf_results = [str(record).strip() for record in vcf_records_iterator]

                                    current_node_copy["vcf_query_results"] = vcf_results
                                    total_nodes_queried_with_pysam += 1
                                    batch_queried_for_vcf_step4 += 1

                                    if not vcf_results:
                                        include_node_in_final_output = False
                                        total_nodes_filtered_out_by_empty_vcf += 1
                                    else:
                                        total_nodes_kept_with_vcf_results += 1
                                        batch_kept_with_vcf_results_step4 += 1
                        except (ValueError, KeyError) as e:
                            include_node_in_final_output = False
                            total_nodes_filtered_out_by_vcf_prereq += 1

                    if include_node_in_final_output:
                        ultimate_filtered_nodes_list.append(current_node_copy)

                    batch_node_counter_step4 += 1
                    if (batch_node_counter_step4 % 100000 == 0 and batch_node_counter_step4 > 0) or \
                            (node_idx_step4 + 1 == num_initial_common_nodes and batch_node_counter_step4 > 0):
                        batch_time_step4 = time.time() - batch_start_time_step4
                        print(
                            f"  Step 4 Batch: Processed {batch_node_counter_step4} nodes "
                            f"({node_idx_step4 + 1}/{num_initial_common_nodes} total common). "
                            f"VCF queried in batch: {batch_queried_for_vcf_step4}. "
                            f"Kept with VCF results in batch: {batch_kept_with_vcf_results_step4}. "
                            f"Time: {batch_time_step4:.2f}s."
                        )
                        batch_node_counter_step4 = 0
                        batch_queried_for_vcf_step4 = 0
                        batch_kept_with_vcf_results_step4 = 0
                        batch_start_time_step4 = time.time()
        except FileNotFoundError:
            print(f"Error: VCF file not found at '{vcf_file}'. Cannot perform VCF filtering.", file=sys.stderr)
            ultimate_filtered_nodes_list = [json_nodes_map[node_id] for node_id in common_node_ids_int]
        except ValueError as e:
            print(f"Error processing VCF file '{vcf_file}'. Is it a valid bgzipped VCF with a .tbi index? Details: {e}",
                  file=sys.stderr)
            ultimate_filtered_nodes_list = [json_nodes_map[node_id] for node_id in common_node_ids_int]

    output_json_structure = copy.deepcopy(main_json_structure)
    output_json_structure["nodes"] = ultimate_filtered_nodes_list
    num_ultimate_nodes = len(ultimate_filtered_nodes_list)
    print(f"\nStep 5: Final processing complete. {num_ultimate_nodes} nodes will be in output JSON.")

    if vcf_file:
        print(f"VCF Query Summary:")
        print(f"  Total nodes eligible and queried with pysam: {total_nodes_queried_with_pysam}")
        print(f"  Nodes filtered out due to VCF prerequisites: {total_nodes_filtered_out_by_vcf_prereq}")
        print(f"  Nodes filtered out due to empty VCF query results: {total_nodes_filtered_out_by_empty_vcf}")
        print(f"  Nodes kept that had non-empty VCF results: {total_nodes_kept_with_vcf_results}")

    if txt_output_path:
        print(f"\nStep 6: Writing {num_ultimate_nodes} ultimate filtered node ID(s) to {txt_output_path}...")
        try:
            with open(txt_output_path, 'w') as f_txt:
                ids_to_write = [str(node.get("node_id")) for node in ultimate_filtered_nodes_list if
                                node.get("node_id") is not None]
                try:
                    sorted_ids = sorted(ids_to_write, key=int)
                except ValueError:
                    sorted_ids = sorted(ids_to_write)
                for node_id_str in sorted_ids:
                    f_txt.write(f"{node_id_str}\n")
            print(f"Successfully wrote ultimate filtered node IDs to {txt_output_path}")
        except Exception as e:
            print(f"Error writing ultimate filtered node IDs to TXT {txt_output_path}: {e}", file=sys.stderr)

    print(f"\nStep 7: Writing final JSON (with {num_ultimate_nodes} nodes) to {output_json_filepath}...")
    try:
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(output_json_structure, f_out_json, indent=4)
        print(f"Successfully wrote resulting JSON to {output_json_filepath}")
    except Exception as e:
        print(f"Error writing output JSON {output_json_filepath}: {e}", file=sys.stderr)

# ─────────────────────────────────────────────────────────────────────────────
# CLI (mutually exclusive idx source)
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter 'nodes' in JSON based on union of LATEST .idx files (30B entries) "
                    "OR a single npu_paths_file (prefix list; '.idx' appended if missing)."
    )
    parser.add_argument("json_path", help="Path to the input JSON file.")
    parser.add_argument("output_json_path", help="Path for the output filtered JSON file.")
    parser.add_argument("idx_paths", nargs='*',
                        help="Zero or more explicit .idx files. Cannot be used with --npu_paths_file.")
    parser.add_argument("--npu_paths_file",
                        help="File containing idx prefixes (one per line). "
                             "Lines without '.idx' will have '.idx' appended.",
                        default=None)
    parser.add_argument("--vcf_file",
                        help="Optional: bgzipped & tabix-indexed VCF (.vcf.gz). Filters nodes whose region has no variants.",
                        default=None)
    parser.add_argument("--txt", dest="txt_output_path",
                        help="Optional: write final node IDs (one per line).", default=None)
    args = parser.parse_args()

    # Require exactly one source of idx input
    if (not args.idx_paths and not args.npu_paths_file) or (args.idx_paths and args.npu_paths_file):
        parser.error("Provide EITHER one or more idx_paths OR --npu_paths_file, but not both.")

    if args.npu_paths_file:
        all_idx_paths = read_idx_paths_from_prefix_file(args.npu_paths_file)
    else:
        all_idx_paths = args.idx_paths

    if not all_idx_paths:
        sys.exit("No valid idx files found.")

    filter_json_nodes_and_write(
        args.json_path,
        all_idx_paths,
        args.output_json_path,
        args.vcf_file,
        args.txt_output_path
    )

if __name__ == "__main__":
    main()
