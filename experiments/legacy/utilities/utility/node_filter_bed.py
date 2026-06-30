#!/usr/bin/env python3
import json
import struct
import argparse
import sys
import time
import copy
import subprocess
import shutil  # For checking tool availability
import re  # For extracting chromosome name


# ─────────────────────────────────────────────────────────────────────────────
# Functions for reading IDX file (UNCHANGED)
# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path):
    """
    Loads the node index from a .idx file into a dictionary.
    """
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}. File might be empty or corrupted.",
                      file=sys.stderr)
                return {}
            blocks_num, = struct.unpack("<I", blocks_num_bytes)
            for i in range(blocks_num):
                header_data = f.read(4 + 8 + 4 + 4 + 2)
                if len(header_data) < (4 + 8 + 4 + 4 + 2):
                    # End of file reached prematurely
                    break
                block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    skipped_bytes = f.read(metadata_len)
                    if len(skipped_bytes) < metadata_len:
                        # End of file reached prematurely
                        break
                node_index[block_id] = {"start": block_start, "size": block_size, "n_records": n_records}
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except struct.error as e:
        print(f"Error: Could not unpack data from IDX file {idx_path}. Details: {e}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}
    return node_index


# ─────────────────────────────────────────────────────────────────────────────
# Function for reading JSON data, extracting Node IDs, and creating a node map (UNCHANGED)
# ─────────────────────────────────────────────────────────────────────────────
def load_json_data_ids_and_map(json_filepath):
    """
    Loads the main JSON file, extracts all node IDs, and creates a map for quick lookups.
    """
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


# ─────────────────────────────────────────────────────────────────────────────
# Function to extract chromosome from path pattern (UNCHANGED)
# ─────────────────────────────────────────────────────────────────────────────
def extract_chromosome_from_path_pattern(path_pattern):
    """
    Extracts a chromosome name (e.g., 'chr1', 'chrX') from a string pattern.
    """
    if not path_pattern or not isinstance(path_pattern, str):
        return None
    match = re.search(r'(chr([0-9A-Za-z_]+|X|Y|M|MT))', path_pattern)
    if match:
        return match.group(1)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# OPTIMIZED Main processing function
# ─────────────────────────────────────────────────────────────────────────────
def filter_json_nodes_and_write(json_filepath, idx_filepath, output_json_filepath, bed_file=None, txt_output_path=None):
    """
    Main function to coordinate the loading, filtering, and writing of data.
    This version is OPTIMIZED to call bedtools only once for performance.
    """
    bedtools_path = None
    if bed_file:
        bedtools_path = shutil.which("bedtools")
        if not bedtools_path:
            print("Warning: --bed_file provided, but bedtools not found in PATH. BED file filtering will be disabled.",
                  file=sys.stderr)
        else:
            print(f"Using bedtools found at: {bedtools_path}")

    print("Step 1: Loading JSON data and creating node map...")
    main_json_structure, target_node_ids_from_json, json_nodes_map = load_json_data_ids_and_map(json_filepath)
    if main_json_structure is None:
        print("Critical error: Could not load main JSON structure. Cannot proceed.", file=sys.stderr)
        return
    print(f"Loaded {len(json_nodes_map)} nodes from JSON for potential processing.")

    chromosome_for_filter = None
    if bed_file and bedtools_path:
        path_pattern = main_json_structure.get("path_name_input_pattern")
        chromosome_for_filter = extract_chromosome_from_path_pattern(path_pattern)
        if not chromosome_for_filter:
            print(f"Warning: Could not extract chromosome from 'path_name_input_pattern': '{path_pattern}'. "
                  "Filtering will fail for nodes unless they have required coordinate information.", file=sys.stderr)
        else:
            print(f"Extracted chromosome '{chromosome_for_filter}' for use in filtering.")

    print("\nStep 2: Loading node index from IDX file...")
    idx_data = load_index(idx_filepath)
    available_idx_node_ids = set(idx_data.keys())
    print(f"Found {len(available_idx_node_ids)} unique node IDs in {idx_filepath}")

    print("\nStep 3: Identifying common node IDs between JSON and IDX...")
    common_node_ids_int = target_node_ids_from_json.intersection(available_idx_node_ids)
    num_initial_common_nodes = len(common_node_ids_int)
    if num_initial_common_nodes == 0:
        print("No common node IDs found. Output 'nodes' list will be empty.")
    else:
        print(f"Found {num_initial_common_nodes} common node IDs to process.")

    ultimate_filtered_nodes_list = []

    # If no BED file is provided, just keep all common nodes.
    if not (bed_file and bedtools_path and chromosome_for_filter):
        print("\nStep 4: No BED filtering requested or possible. Keeping all common nodes.")
        for node_id_int in sorted(list(common_node_ids_int)):
            if node_id_int in json_nodes_map:
                ultimate_filtered_nodes_list.append(json_nodes_map[node_id_int])
    else:
        # --- OPTIMIZED BED FILTERING LOGIC ---
        print(
            f"\nStep 4: Preparing all {num_initial_common_nodes} common nodes for a single, efficient bedtools query...")

        all_nodes_bed_string = ""
        nodes_for_query_count = 0
        nodes_failed_prereq_count = 0

        # Create a single string containing all node regions in BED format.
        # This is much faster than writing to a temporary file.
        for node_id_int in sorted(list(common_node_ids_int)):
            node_item = json_nodes_map.get(node_id_int)
            if not node_item: continue

            try:
                start_0based_str = node_item.get("grch38_position_start")
                length_str = node_item.get("length")

                if start_0based_str is None or length_str is None:
                    nodes_failed_prereq_count += 1
                    continue

                start_0based = int(start_0based_str)
                length = int(length_str)

                if length <= 0:
                    nodes_failed_prereq_count += 1
                    continue

                # BED format: chrom, start, end, name (we use node_id as the name)
                end_0based = start_0based + length
                all_nodes_bed_string += f"{chromosome_for_filter}\t{start_0based}\t{end_0based}\t{node_id_int}\n"
                nodes_for_query_count += 1

            except (ValueError, TypeError):
                nodes_failed_prereq_count += 1
                continue

        print(f"Prepared {nodes_for_query_count} valid node regions for bedtools.")
        print(f"{nodes_failed_prereq_count} nodes were skipped due to missing or invalid coordinate/length data.")

        kept_node_ids = set()
        if nodes_for_query_count > 0:
            print("\nStep 5: Running a single 'bedtools intersect' command on all nodes at once...")
            # Use '-wa' to write the original entry from 'stdin' for each overlap.
            command = [bedtools_path, 'intersect', '-a', 'stdin', '-b', bed_file, '-wa']

            start_time = time.time()
            process = subprocess.run(
                command,
                input=all_nodes_bed_string,
                capture_output=True,
                text=True,
                check=False
            )
            end_time = time.time()
            print(f"Bedtools command finished in {end_time - start_time:.2f} seconds.")

            if process.returncode != 0 and process.stderr:
                print(f"Error/Warning from bedtools: {process.stderr.strip()}", file=sys.stderr)

            # The output of bedtools will be the lines from our input string that had an overlap.
            # We just need to parse the node ID (the 4th column) from each line.
            for line in process.stdout.strip().split('\n'):
                if not line: continue
                fields = line.split('\t')
                if len(fields) >= 4:
                    try:
                        kept_node_ids.add(int(fields[3]))
                    except (ValueError, IndexError):
                        print(f"Warning: Could not parse node ID from bedtools output line: '{line}'", file=sys.stderr)

        print(f"Identified {len(kept_node_ids)} nodes that overlap with regions in the BED file.")

        # Build the final list of node objects from the set of kept IDs.
        for node_id_int in kept_node_ids:
            if node_id_int in json_nodes_map:
                ultimate_filtered_nodes_list.append(json_nodes_map[node_id_int])

        # Sort the final list by node ID for consistent output.
        ultimate_filtered_nodes_list.sort(key=lambda x: int(x['node_id']))

        # --- Final Summary ---
        print("\nBED Filter Summary:")
        print(f"  - Nodes eligible and queried with bedtools: {nodes_for_query_count}")
        print(f"  - Nodes filtered out due to missing prerequisites: {nodes_failed_prereq_count}")
        print(f"  - Nodes kept that overlapped with a BED region: {len(kept_node_ids)}")
        print(f"  - Nodes filtered out due to no overlap in BED file: {nodes_for_query_count - len(kept_node_ids)}")

    num_ultimate_nodes = len(ultimate_filtered_nodes_list)
    final_step_number = 6 if bed_file and bedtools_path else 5
    print(
        f"\nStep {final_step_number}: Final processing complete. A total of {num_ultimate_nodes} nodes will be written to the output.")

    if txt_output_path:
        print(
            f"\nStep {final_step_number + 1}: Writing {num_ultimate_nodes} filtered node ID(s) to {txt_output_path}...")
        try:
            with open(txt_output_path, 'w') as f_txt:
                if ultimate_filtered_nodes_list:
                    ids_to_write = [str(node.get("node_id")) for node in ultimate_filtered_nodes_list if
                                    node.get("node_id") is not None]
                    f_txt.write('\n'.join(ids_to_write) + '\n')
            print(f"Successfully wrote node IDs to {txt_output_path}")
        except Exception as e:
            print(f"Error writing filtered node IDs to TXT file {txt_output_path}: {e}", file=sys.stderr)

    print(
        f"\nStep {final_step_number + 2}: Writing final JSON (with {num_ultimate_nodes} nodes) to {output_json_filepath}...")
    try:
        # Rebuild the final JSON structure
        output_json_structure = copy.deepcopy(main_json_structure)
        output_json_structure["nodes"] = ultimate_filtered_nodes_list
        with open(output_json_filepath, 'w') as f_out_json:
            json.dump(output_json_structure, f_out_json, indent=4)
        print(f"Successfully wrote resulting JSON to {output_json_filepath}")
    except Exception as e:
        print(f"Error writing output JSON {output_json_filepath}: {e}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# REVISED Command-line interface
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Filter 'nodes' in a JSON file based on an IDX file and, optionally, a BED file of high-confidence regions. "
                    "Outputs a modified JSON and optionally a TXT list of the kept node IDs.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("json_path", help="Path to the input JSON file.")
    parser.add_argument("idx_path", help="Path to the corresponding .idx file.")
    parser.add_argument("output_json_path", help="Path for the output filtered JSON file.")
    parser.add_argument("--bed_file",
                        help="Optional: Path to a BED file containing regions for filtering.\n"
                             "Nodes will be kept only if their genomic coordinates overlap with a region in this file.\n"
                             "Requires 'bedtools' to be installed and in the system's PATH.",
                        default=None)
    parser.add_argument("--txt", dest="txt_output_path",
                        help="Optional: Path to write the final list of filtered node IDs (one per line).",
                        default=None)
    args = parser.parse_args()

    filter_json_nodes_and_write(
        args.json_path,
        args.idx_path,
        args.output_json_path,
        bed_file=args.bed_file,
        txt_output_path=args.txt_output_path
    )


if __name__ == "__main__":
    main()
