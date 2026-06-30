import os
import shutil
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


def restore_single_npy_file(npy_file_path, restore_base_path):
    """
    Restores a single .npy file to its original node directory.
    It parses the node_id from the filename prefix.

    Args:
        npy_file_path (Path): The full path to the prefixed .npy file.
        restore_base_path (Path): The base path where the restored node directory will be created.
    """
    filename = npy_file_path.name

    # --- 1. Parse node_id and original filename ---
    # Assumes the format is "nodeid_originalfilename.npy"
    try:
        node_id, original_filename = filename.split('_', 1)
    except ValueError:
        tqdm.write(f"\nWarning: Skipping file with unexpected format: {filename}")
        return

    # --- 2. Create destination directory (thread-safe) ---
    node_restore_path = restore_base_path / node_id
    node_restore_path.mkdir(exist_ok=True)

    # --- 3. Move and rename the file ---
    destination_path = node_restore_path / original_filename
    try:
        shutil.move(str(npy_file_path), str(destination_path))
    except Exception as e:
        tqdm.write(f"\nWarning: Could not move file '{filename}': {e}")


def restore_npy_files(source_dir, output_dir, num_workers):
    """
    Restores .npy files from a central 'tensors' directory back to
    individual node_id subdirectories in parallel.

    Args:
        source_dir (str): The directory containing the 'tensors' folder.
        output_dir (str): The directory where the original structure will be restored.
        num_workers (int): The number of parallel threads to use.
    """
    # --- 1. Validate Paths ---
    source_path = Path(source_dir)
    output_path = Path(output_dir)
    tensors_source_path = source_path / "tensors"

    if not all([source_path.is_dir(), tensors_source_path.is_dir()]):
        print(f"Error: Source directory '{source_dir}' is not valid.")
        print("It must contain a 'tensors' subdirectory.")
        return

    output_path.mkdir(exist_ok=True)
    print(f"Restored data will be saved in: '{output_path}'")

    # --- 2. Discover .npy files to move ---
    npy_files = list(tensors_source_path.glob("*.npy"))
    if not npy_files:
        print("Warning: No .npy files found in the 'tensors' directory.")
        return

    print(f"\nFound {len(npy_files)} .npy files to restore with {num_workers} workers...")

    # --- 3. Process Files in Parallel ---
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_file = {
            executor.submit(restore_single_npy_file, npy_file, output_path): npy_file.name
            for npy_file in npy_files
        }

        for future in tqdm(as_completed(future_to_file), total=len(npy_files), desc="Restoring .npy files"):
            try:
                future.result()  # Check for exceptions from the worker thread
            except Exception as exc:
                filename = future_to_file[future]
                tqdm.write(f"\nError processing file {filename}: {exc}")

    print("\n--- Restoration Complete ---")
    print(f"Total .npy files moved back: {len(npy_files)}")
    print(f"Restored data is located in: '{output_path}'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Restore .npy files from an organized 'tensors' directory back to node_id subfolders.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("source_directory", help="The source directory containing the 'tensors' folder.")
    parser.add_argument("output_directory", help="The directory where the original folder structure will be restored.")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Number of parallel worker threads to use.")

    args = parser.parse_args()

    restore_npy_files(args.source_directory, args.output_directory, args.workers)
