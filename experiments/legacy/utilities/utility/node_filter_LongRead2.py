#!/usr/bin/env python3
import json
import struct
import argparse
import sys
import time
import copy
import re
import os
import pysam  # Replaces subprocess and shutil for VCF operations

# ─────────────────────────────────────────────────────────────────────────────
# New-format constants (only used if you pass dat_path to load_index)
GLOBAL_MAGIC = b"MYFMT\x01"                       # from the writer
GLOBAL_VER_PACK = struct.Struct("<BBI16s")        # major, minor, block_count, reserved[16]
BLOCK_HDR_PACK = struct.Struct("<I I H I")        # node_id, n_records, flags, node_length
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size              # 14 bytes

# ─────────────────────────────────────────────────────────────────────────────
# Functions for reading IDX file (UPDATED to support old/new formats)
# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path, dat_path=None):
    """
    Read .idx produced by your new writer. Supports:
      - NEW: per-entry 26 bytes: <I Q I I H I>
      - OLD: per-entry 22 bytes: <I Q I I H>
    Returns dict[node_id] = {
        "start": offset,
        "block_size": block_size,
        "n_records": n_records,
        "flags": flags,
        "node_length": node_length_or_0
    }
    If dat_path is provided and node_length == 0 (old idx), fill node_length by reading .dat.
    """
    node_index = {}

    try:
        with open(idx_path, "rb") as f:
            raw = f.read(4)
            if len(raw) != 4:
                print(f"Error: Could not read block count from {idx_path}", file=sys.stderr)
                return {}
            (count,) = struct.unpack("<I", raw)

            # determine per-entry size: total_remaining // count
            f.seek(0, os.SEEK_END)
            remaining = f.tell() - 4
            if count == 0 or remaining <= 0:
                print(f"Error: {idx_path} is empty or corrupt.", file=sys.stderr)
                return {}
            entry_size = remaining // count

            # old (22) vs new (26)
            if entry_size not in (22, 26):
                # Be generous: default to new format if not cleanly divisible
                entry_size = 26

            f.seek(4)
            for i in range(count):
                data = f.read(entry_size)
                if len(data) != entry_size:
                    print(f"Error: Truncated .idx at entry {i} in {idx_path}", file=sys.stderr)
                    return {}

                if entry_size == 26:
                    node_id, offset, block_size, n_records, flags, node_len = struct.unpack("<I Q I I H I", data)
                else:  # 22 bytes (old)
                    node_id, offset, block_size, n_records, flags = struct.unpack("<I Q I I H", data)
                    node_len = 0

                node_index[node_id] = {
                    "start": offset,
                    "block_size": block_size,
                    "n_records": n_records,
                    "flags": flags,
                    "node_length": node_len,
                }

    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except struct.error as e:
        print(f"Error: Could not unpack data from IDX file {idx_path}. Details: {e}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}

    # If we have old idx (node_length == 0) but the .dat is available, fill lengths from .dat
    if dat_path and any(info["node_length"] <= 0 for info in node_index.values()):
        try:
            with open(dat_path, "rb") as df:
                magic = df.read(len(GLOBAL_MAGIC))
                if magic != GLOBAL_MAGIC:
                    print(f"Warning: .dat magic mismatch for {dat_path}; cannot fill node_length.", file=sys.stderr)
                else:
                    maj, minr, dat_count, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))
                    # Optional informational note about writer version
                    if minr < 3:
                        print(f"Warning: .dat minor version {maj}.{minr} predates node-length CIGAR; "
                              f"index-filling still works but records may be old-format.", file=sys.stderr)
                    # read each block header by seeking to offset
                    for node_id, info in node_index.items():
                        if info["node_length"] > 0:
                            continue
                        df.seek(info["start"], os.SEEK_SET)
                        hdr = df.read(BLOCK_HDR_SIZE)
                        if len(hdr) != BLOCK_HDR_SIZE:
                            print(f"Warning: Cannot read block header for node {node_id} in {dat_path}", file=sys.stderr)
                            continue
                        nid2, nrec2, flg2, node_len = BLOCK_HDR_PACK.unpack(hdr)
                        if nid2 != node_id or nrec2 != info["n_records"]:
                            print(f"Warning: .dat/.idx mismatch for node {node_id} (idx n={info['n_records']}, dat n={nrec2})", file=sys.stderr)
                        info["node_length"] = node_len
        except FileNotFoundError:
            print(f"Note: .dat not found at {dat_path}; node_length will remain 0 for old .idx.", file=sys.stderr)
        except Exception as e:
            print(f"Note: Could not read .dat to fill node_length: {e}", file=sys.stderr)

    return node_index

# ─────────────────────────────────────────────────────────────────────────────
# Function for reading JSON data (Unchanged)
# ─────────────────────────────────────────────────────────────────────────────
def load_json_data_ids_and_map(json_filepath):
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
                if not isinstance(item, dict): continue
                try:
                    node_id_str = item.get("node_id")
                    if node_id_str is None: continue
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

# ─────────────────────────────────────────────────────────────────────────────
# Function to extract chromosome from path pattern (Unchanged)
# ─────────────────────────────────────────────────────────────────────────────
def extract_chromosome_from_path_pattern(path_pattern):
    if not path_pattern or not isinstance(path_pattern, str):
        return None
    match = re.search(r'(chr([0-9A-Za-z_]+|X|Y|M|MT))', path_pattern)
    if match:
        return match.group(1)
    return None

# ─────────────────────────────────────────────────────────────────────────────
# Main processing function (unchanged, except calling load_index with optional dat)
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath, vcf_file=None, txt_output_path=None, dat_filepath=None):
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

    print("\nStep 2: Loading node index (IDX)...")
    # NEW: pass dat_filepath (optional) so we can fill node_length if idx is old
    idx_data = load_index(idx_filepath, dat_path=dat_filepath)
    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs...")
    common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)
    num_initial_common_nodes = len(common_node_ids_int)
    if num_initial_common_nodes == 0:
        print("No common node IDs found between JSON and IDX. Output 'nodes' list will be empty.")
    else:
        print(f"Found {num_initial_common_nodes} common node IDs between JSON and IDX for potential processing.")

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
                    if (batch_node_counter_step4 % 5000 == 0 and batch_node_counter_step4 > 0) or \
                            (node_idx_step4 + 1 == num_initial_common_nodes and batch_node_counter_step4 > 0):
                        batch_time_step4 = time.time() - batch_start_time_step4
                        print(
                            f"  Step 4 Batch: Processed {batch_node_counter_step4} nodes ({node_idx_step4 + 1}/{num_initial_common_nodes} total common). "
                            f"VCF queried in batch: {batch_queried_for_vcf_step4}. Kept with VCF results in batch: {batch_kept_with_vcf_results_step4}. "
                            f"Time: {batch_time_step4:.2f}s.")
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
        print(
            f"  Nodes filtered out due to VCF prerequisites (missing coords/length, etc.): {total_nodes_filtered_out_by_vcf_prereq}")
        print(f"  Nodes filtered out due to empty VCF query results: {total_nodes_filtered_out_by_empty_vcf}")
        print(f"  Nodes kept that had non-empty VCF results: {total_nodes_kept_with_vcf_results}")

    if txt_output_path:
        print(f"\nStep 6: Writing {num_ultimate_nodes} ultimate filtered node ID(s) to {txt_output_path}...")
        if ultimate_filtered_nodes_list:
            try:
                with open(txt_output_path, 'w') as f_txt:
                    ids_to_write = [str(node.get("node_id")) for node in ultimate_filtered_nodes_list if
                                    node.get("node_id") is not None]
                    try:
                        sorted_ids = sorted(ids_to_write, key=int)
                    except ValueError:
                        sorted_ids = sorted(ids_to_write)
                    for node_id_str in sorted_ids: f_txt.write(f"{node_id_str}\n")
                print(f"Successfully wrote ultimate filtered node IDs to {txt_output_path}")
            except Exception as e:
                print(f"Error writing ultimate filtered node IDs to TXT {txt_output_path}: {e}", file=sys.stderr)
        else:
            print(f"No ultimate filtered node IDs to write to {txt_output_path}.")
            try:
                with open(txt_output_path, 'w') as f_txt:
                    pass
                print(f"Created empty TXT file at {txt_output_path}.")
            except Exception as e:
                print(f"Error creating empty TXT file {txt_output_path}: {e}", file=sys.stderr)

    print(f"\nStep 7: Writing final JSON (with {num_ultimate_nodes} nodes) to {output_json_filepath}...")
    try:
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(output_json_structure, f_out_json, indent=4)
        print(f"Successfully wrote resulting JSON to {output_json_filepath}")
    except Exception as e:
        print(f"Error writing output JSON {output_json_filepath}: {e}", file=sys.stderr)

# ─────────────────────────────────────────────────────────────────────────────
# CLI (added --dat for optional node_length fill)
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter 'nodes' list in a JSON based on IDX file and VCF content using pysam. Supports new/old .idx formats."
    )
    parser.add_argument("json_path", help="Path to the input JSON file.")
    parser.add_argument("idx_path", help="Path to the .idx file.")
    parser.add_argument("output_json_path", help="Path for the output filtered JSON file.")
    parser.add_argument("--vcf_file",
                        help="Optional: Path to a bgzipped and tabix-indexed VCF file (.vcf.gz). Nodes will be filtered if they lack VCF prerequisites OR if their region yields no variants from this VCF.",
                        default=None)
    parser.add_argument("--txt", dest="txt_output_path",
                        help="Optional: Path to output filtered node IDs (one per line).", default=None)
    parser.add_argument("--dat", dest="dat_filepath",
                        help="Optional: Path to the matching .dat file (used to fill node_length for old .idx).",
                        default=None)
    args = parser.parse_args()

    filter_json_nodes_and_write(
        args.json_path,
        args.idx_path,
        args.output_json_path,
        args.vcf_file,
        args.txt_output_path,
        args.dat_filepath
    )

if __name__ == "__main__":
    main()
