import os
import json
import shutil
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


def process_single_node(node_dir, tensors_path):
    """
    Processes a single node directory: moves .npy files and reads the summary JSON.
    This function is designed to be run in a separate thread.

    Args:
        node_dir (Path): The path to the individual node_id directory.
        tensors_path (Path): The destination directory for tensor files.

    Returns:
        tuple: A tuple containing (node_id, loaded_json_data, num_npy_files_moved).
               Returns (node_id, None, num_npy_files_moved) if JSON is missing or invalid.
    """
    node_id = node_dir.name
    npy_files_moved = 0
    json_data = None

    # --- Move and rename .npy files ---
    try:
        for npy_file in node_dir.glob("*.npy"):
            new_filename = f"{node_id}_{npy_file.name}"
            destination_path = tensors_path / new_filename
            shutil.move(str(npy_file), str(destination_path))
            npy_files_moved += 1
    except Exception as e:
        # Return the error to be handled by the main thread
        raise IOError(f"Failed to move .npy files for node {node_id}: {e}")

    # --- Read the variant_summary.json ---
    summary_file = node_dir / "variant_summary.json"
    if summary_file.exists():
        try:
            with open(summary_file, 'r') as f:
                json_data = json.load(f)
        except (json.JSONDecodeError, Exception) as e:
            # Log the error but don't stop the whole process
            print(f"\nWarning: Could not parse or read JSON for node {node_id}: {e}")
            # json_data remains None

    return node_id, json_data, npy_files_moved


def organize_node_data(source_dir, output_dir, num_workers):
    """
    Organizes node data in parallel by moving and renaming .npy files
    and merging variant_summary.json files.

    Args:
        source_dir (str): The path to the directory containing node_id subdirectories.
        output_dir (str): The path to the directory where organized data will be saved.
        num_workers (int): The number of parallel threads to use.
    """
    # --- 1. Setup Output Directories ---
    source_path = Path(source_dir)
    output_path = Path(output_dir)
    tensors_path = output_path / "tensors"

    if not source_path.is_dir():
        print(f"Error: Source directory not found at '{source_dir}'")
        return

    output_path.mkdir(exist_ok=True)
    tensors_path.mkdir(exist_ok=True)
    print(f"Output will be saved in: '{output_path}'")
    print(f"Tensors will be moved to: '{tensors_path}'")

    # --- 2. Discover Tasks ---
    node_dirs = [d for d in source_path.iterdir() if d.is_dir() and d.name.isdigit()]
    if not node_dirs:
        print(f"Warning: No node_id subdirectories found in '{source_dir}'")
        return

    print(f"\nFound {len(node_dirs)} node directories to process with {num_workers} workers...")

    # --- 3. Process Data in Parallel ---
    merged_summary_data = {}
    total_npy_files_moved = 0
    nodes_processed_count = 0  # <-- ADDED: Counter for milestone printing

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all tasks to the thread pool
        future_to_node = {
            executor.submit(process_single_node, node_dir, tensors_path): node_dir.name
            for node_dir in node_dirs
        }

        # Use tqdm to create a progress bar as tasks complete
        for future in tqdm(as_completed(future_to_node), total=len(node_dirs), desc="Processing nodes"):
            nodes_processed_count += 1  # <-- ADDED: Increment counter
            node_id_str = future_to_node[future]
            try:
                node_id, data, npy_count = future.result()
                total_npy_files_moved += npy_count
                if data:
                    merged_summary_data[node_id] = data
                else:
                    # This case is handled by the warning inside process_single_node
                    pass
            except Exception as exc:
                print(f"\nError processing node {node_id_str}: {exc}")

            # --- ADDED: Milestone printing logic ---
            if nodes_processed_count % 10000 == 0 and nodes_processed_count > 0:
                tqdm.write(f"--- Milestone: {nodes_processed_count}/{len(node_dirs)} nodes processed. ---")

    # --- 4. Write the Final Merged JSON ---
    output_json_path = output_path / "merged_variant_summary.json"
    print("\nWriting final merged JSON file...")
    try:
        # Sort keys for consistent output, which is good practice
        sorted_merged_data = dict(sorted(merged_summary_data.items(), key=lambda item: int(item[0])))
        with open(output_json_path, 'w') as f:
            json.dump(sorted_merged_data, f, indent=2)
    except Exception as e:
        print(f"\nError: Failed to write merged JSON file: {e}")
        return

    print("\n--- Organization Complete ---")
    print(f"Total .npy files moved: {total_npy_files_moved}")
    print(f"Total nodes merged in JSON: {len(merged_summary_data)}")
    print(f"Merged summary saved to: '{output_json_path}'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Organize node data in parallel by moving .npy files and merging summary JSONs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("source_directory", help="The source directory containing the node_id subfolders.")
    parser.add_argument("output_directory", help="The directory where the organized data will be saved.")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Number of parallel worker threads to use.")

    args = parser.parse_args()

    organize_node_data(args.source_directory, args.output_directory, args.workers)
