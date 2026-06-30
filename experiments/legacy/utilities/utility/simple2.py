import argparse
import struct
import sys


def load_index(idx_path):
    """
    Loads node IDs from a binary .idx file.

    The structure of the .idx file is assumed to be:
    - 4 bytes: Number of blocks (unsigned int)
    - For each block:
        - 4 bytes: block_id (unsigned int)
        - 8 bytes: block_start (unsigned long long)
        - 4 bytes: block_size (unsigned int)
        - 4 bytes: n_records (unsigned int)
        - 2 bytes: metadata_len (unsigned short)
        - metadata_len bytes: metadata (to be skipped)

    Args:
        idx_path (str): The path to the .idx file.

    Returns:
        dict: A dictionary mapping node IDs (block_id) to their info.
              Returns an empty dictionary if an error occurs.
    """
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            # Read the total number of blocks
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}. File might be empty or corrupted.",
                      file=sys.stderr)
                return {}

            blocks_num, = struct.unpack("<I", blocks_num_bytes)

            # Read each block's header
            for i in range(blocks_num):
                header_data = f.read(4 + 8 + 4 + 4 + 2)
                if len(header_data) < (4 + 8 + 4 + 4 + 2):
                    print(f"Warning: Reached end of file unexpectedly while reading block {i + 1}/{blocks_num}.",
                          file=sys.stderr)
                    break

                block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", header_data)

                # Skip metadata if it exists
                if metadata_len > 0:
                    skipped_bytes = f.read(metadata_len)
                    if len(skipped_bytes) < metadata_len:
                        print(f"Warning: Reached end of file unexpectedly while skipping metadata for block {i + 1}.",
                              file=sys.stderr)
                        break

                node_index[block_id] = {"start": block_start, "size": block_size, "n_records": n_records}
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except struct.error as e:
        print(f"Error: Could not unpack data from IDX file {idx_path}. The file may be corrupted. Details: {e}",
              file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}

    return node_index


def find_matches_in_txt(txt_path, idx_nodes):
    """
    Finds all node IDs from a text file that are present in a given set of nodes.
    Each line in the text file is assumed to contain one node ID.

    Args:
        txt_path (str): The path to the input text file.
        idx_nodes (set): A set of node IDs loaded from the index file.

    Returns:
        set: A set of all unique matching node IDs found. Returns None on file error.
    """
    matching_nodes = set()
    line_num = 0
    try:
        with open(txt_path, 'r') as f:
            for line in f:
                line_num += 1
                node_id_str = line.strip()

                # Skip empty lines
                if not node_id_str:
                    continue

                try:
                    # Convert the line to an integer and check for its existence
                    node_id = int(node_id_str)
                    if node_id in idx_nodes:
                        matching_nodes.add(node_id)
                except ValueError:
                    print(f"Warning: Skipping non-integer value on line {line_num}: '{node_id_str}'", file=sys.stderr)
                    continue
    except FileNotFoundError:
        print(f"Error: TXT file not found at {txt_path}", file=sys.stderr)
        return None  # Return None on error
    except Exception as e:
        print(f"An unexpected error occurred while reading TXT file {txt_path}: {e}", file=sys.stderr)
        return None

    return matching_nodes


def main():
    """
    Main function to orchestrate the node finding process.
    """
    parser = argparse.ArgumentParser(
        description="Find all unique node IDs from a text file that also exist in a binary .idx file."
    )
    parser.add_argument("txt_file", help="Path to the input text file containing one node ID per line.")
    parser.add_argument("idx_file", help="Path to the .idx index file.")
    args = parser.parse_args()

    # 1. Load node IDs from the index file
    print(f"Loading nodes from index file: {args.idx_file}...", file=sys.stderr)
    index_data = load_index(args.idx_file)
    if not index_data:
        print("Could not load index file. Exiting.", file=sys.stderr)
        sys.exit(1)

    idx_nodes = set(index_data.keys())
    print(f"Loaded {len(idx_nodes)} unique nodes from the index.", file=sys.stderr)

    # 2. Find all matching nodes in the text file
    print(f"Finding all matching nodes from TXT file: {args.txt_file}...", file=sys.stderr)
    found_nodes = find_matches_in_txt(args.txt_file, idx_nodes)

    # 3. Output the final result
    if found_nodes is not None:
        print("\n--- Matching Nodes Found ---")
        if found_nodes:
            # Sort the nodes for consistent output
            # for node in sorted(list(found_nodes)):
            #     print(node)
            print(f"\nTotal unique matching nodes: {len(found_nodes)}", file=sys.stderr)
        else:
            print("No matching nodes were found.")


if __name__ == "__main__":
    main()
