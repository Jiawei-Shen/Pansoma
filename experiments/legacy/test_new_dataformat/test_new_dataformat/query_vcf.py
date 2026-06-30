import json
import subprocess
import shlex  # For safely formatting command arguments for display
import argparse
import os
import concurrent.futures  # For parallel processing


def execute_bcftools_query(vcf_file_path, region_string):
    """
    Executes a bcftools view query for a given region.

    Args:
        vcf_file_path (str): Path to the VCF file.
        region_string (str): The region to query (e.g., "chr1:100-200").

    Returns:
        tuple: (list_of_variant_strings, bcftools_stderr_str)
               Returns (None, error_message) if bcftools command itself fails.
    """
    # Command to execute: bcftools view <file.vcf.gz> <chr:from-to>
    command = ["bcftools", "view", vcf_file_path, region_string]
    print(f"  [Worker PID:{os.getpid()}] Executing: {' '.join(shlex.quote(arg) for arg in command)}")
    try:
        # Run the command
        process = subprocess.run(
            command,
            capture_output=True,  # Capture stdout and stderr
            text=True,  # Decode output as text
            check=False  # Do not raise an exception for non-zero exit codes automatically
            # We will check stderr manually
        )

        variants_list = []
        if process.stdout:
            # Split stdout into lines, removing any empty lines that might result from split
            variants_list = [line for line in process.stdout.strip().split('\n') if line]

        return variants_list, process.stderr.strip()  # Return list of variants and any stderr output

    except FileNotFoundError:
        # This occurs if 'bcftools' executable is not found
        error_msg = f"Error: 'bcftools' command not found. Please ensure it is installed and in your PATH."
        print(error_msg)  # Print to console for immediate user feedback
        return None, error_msg  # Indicate critical failure
    except Exception as e:
        # Catch any other unexpected errors during subprocess execution
        error_msg = f"An unexpected error occurred while running bcftools for region {region_string}: {e}"
        print(error_msg)  # Print to console
        return None, error_msg


def process_item_for_bcftools(item, vcf_file_path):
    """
    Processes a single item from the input JSON, constructs the region,
    and prepares it for bcftools query.

    Args:
        item (dict): A dictionary from the input JSON (expected to have fields like
                     'top_level_node_id', 'haplotype_id', 'node_sequence_length',
                     'grch38_chromosome_region', 'grch38_start_position').
        vcf_file_path (str): Path to the VCF file.

    Returns:
        dict: A dictionary containing original item info, the query region string,
              a list of found variants, and any error messages from bcftools.
    """
    # Initialize the result entry with all expected fields from the input item
    base_result_entry = {
        "top_level_node_id": item.get("top_level_node_id"),
        "haplotype_id": item.get("haplotype_id"),
        "node_sequence_length": item.get("node_sequence_length"),
        "grch38_chromosome_region": item.get("grch38_chromosome_region"),
        "grch38_start_position": item.get("grch38_start_position"),
        "query_region": None,  # Will be populated if query proceeds
        "variants": [],  # Initialize as empty list for variants
        "bcftools_error": None  # Will store any errors
    }

    # Extract necessary fields from the input item
    start_pos_str = item.get("grch38_start_position")
    seq_len_val = item.get("node_sequence_length")
    chromosome = item.get("grch38_chromosome_region")

    # Validate required fields for constructing the query region
    if not chromosome:
        base_result_entry["bcftools_error"] = "Skipped: Missing 'grch38_chromosome_region'."
        print(
            f"  Skipping item (HapID: {base_result_entry.get('haplotype_id', 'N/A')}): Missing 'grch38_chromosome_region'.")
        return base_result_entry

    if start_pos_str is None:  # Check if start position is null (from previous script error)
        base_result_entry["bcftools_error"] = "Skipped: 'grch38_start_position' is null."
        print(
            f"  Skipping item (HapID: {base_result_entry.get('haplotype_id', 'N/A')}): 'grch38_start_position' is null.")
        return base_result_entry

    if seq_len_val is None:  # Check if sequence length is null
        base_result_entry["bcftools_error"] = "Skipped: 'node_sequence_length' is null."
        print(
            f"  Skipping item (HapID: {base_result_entry.get('haplotype_id', 'N/A')}): 'node_sequence_length' is null.")
        return base_result_entry

    try:
        # Convert start position and sequence length to integers
        start_pos = int(start_pos_str)
        seq_len = int(seq_len_val)
        if seq_len < 0:  # Sequence length cannot be negative
            base_result_entry[
                "bcftools_error"] = f"Skipped: Invalid 'node_sequence_length' ({seq_len}). Must be non-negative."
            print(
                f"  Skipping item (HapID: {base_result_entry.get('haplotype_id', 'N/A')}): Invalid 'node_sequence_length' ({seq_len}).")
            return base_result_entry
    except ValueError:
        # Handle cases where conversion to int fails
        base_result_entry[
            "bcftools_error"] = "Skipped: 'grch38_start_position' or 'node_sequence_length' is not a valid integer."
        print(
            f"  Skipping item (HapID: {base_result_entry.get('haplotype_id', 'N/A')}): 'grch38_start_position' or 'node_sequence_length' not a valid int.")
        return base_result_entry

    # Calculate end position based on the formula: grch38_start_position + node_sequence_length
    # This defines the end of the window.
    # For VCF queries (1-based, inclusive), if node_sequence_length is the number of bases,
    # a common query would be start_pos to (start_pos + seq_len - 1).
    # However, the user specified the formula for the *end coordinate* as start + length.
    # If seq_len is 0, region is start-start.
    # If seq_len is 1, region is start-(start+1).
    end_pos = start_pos + seq_len

    # Sanity check for calculated end position
    if end_pos < start_pos:  # This should only occur if seq_len was negative, which is checked above.
        base_result_entry[
            "bcftools_error"] = f"Skipped: Calculated end position ({end_pos}) is less than start position ({start_pos}) with length ({seq_len})."
        print(
            f"  Skipping item (HapID: {base_result_entry.get('haplotype_id', 'N/A')}): Calculated end_pos {end_pos} < start_pos {start_pos} with length {seq_len}.")
        return base_result_entry

    # Construct the region string for bcftools
    query_region_str = f"{chromosome}:{start_pos}-{end_pos}"
    base_result_entry["query_region"] = query_region_str

    # Execute the bcftools query
    variants_list, bcftools_stderr_output = execute_bcftools_query(vcf_file_path, query_region_str)

    if variants_list is not None:  # Indicates bcftools command was attempted (not a FileNotFoundError for bcftools itself)
        base_result_entry["variants"] = variants_list  # Will be an empty list if no variants found

    if bcftools_stderr_output:  # Store any stderr output from bcftools (warnings or errors)
        # If bcftools itself wasn't found, the error is critical and should have been raised by execute_bcftools_query
        if "bcftools command not found" in bcftools_stderr_output:
            # This case should ideally be caught and raised in execute_bcftools_query
            raise FileNotFoundError(bcftools_stderr_output)
        base_result_entry["bcftools_error"] = bcftools_stderr_output

    return base_result_entry


def main():
    """
    Main function to parse arguments and run the VCF querying process.
    """
    parser = argparse.ArgumentParser(
        description="Query a VCF file using bcftools view based on regions derived from an input JSON file. "
                    "The input JSON should be an array of objects, each with fields like "
                    "'grch38_chromosome_region', 'grch38_start_position', and 'node_sequence_length'.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # Shows default values in help
    )
    parser.add_argument(
        "input_json",
        help="Path to the input JSON file (output from the previous script containing regions)."
    )
    parser.add_argument(
        "vcf_file",
        help="Path to the gzipped and indexed VCF file (e.g., SMaHT_COLO829_SNV_truth_set_v1.0.vcf.gz)."
    )
    parser.add_argument(
        "-o", "--output_json",
        default="bcftools_query_results.json",  # Default output file name
        help="Path to save the output JSON file with bcftools query results."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=min(4, os.cpu_count() if os.cpu_count() else 1),  # Default to a modest number of threads
        help="Number of parallel worker threads to use for bcftools queries."
    )

    args = parser.parse_args()

    # --- Validate input files ---
    if not os.path.exists(args.input_json):
        print(f"Error: Input JSON file not found at '{args.input_json}'")
        return

    if not os.path.exists(args.vcf_file):
        print(f"Error: VCF file not found at '{args.vcf_file}'")
        return
    # Check for VCF index (.tbi or .csi)
    if not (os.path.exists(args.vcf_file + ".tbi") or os.path.exists(args.vcf_file + ".csi")):
        print(f"Error: Index file (.tbi or .csi) not found for '{args.vcf_file}'. "
              "Please index the VCF file using tabix or bcftools index.")
        return

    # --- Load input JSON data ---
    try:
        with open(args.input_json, 'r') as f:
            input_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{args.input_json}'. Make sure it's a valid JSON file.")
        return

    if not isinstance(input_data, list):
        print(f"Error: Input JSON is not a list as expected. Found type: {type(input_data)}")
        return

    all_results = []  # List to store results from all processed items

    print(f"Processing {len(input_data)} items from '{args.input_json}'...")
    print(f"Querying VCF: '{args.vcf_file}' with bcftools view.")
    print(f"Using up to {args.threads} parallel threads.")
    print(f"Output will be saved to: '{args.output_json}'")

    # Using ThreadPoolExecutor because bcftools is an external I/O-bound command.
    # For CPU-bound Python functions, ProcessPoolExecutor is usually better.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        # Create a dictionary to map futures to some identifiable information for logging
        future_to_item_info = {
            executor.submit(process_item_for_bcftools, item, args.vcf_file):
                item.get('haplotype_id', f'index_{i}')  # Use haplotype_id or index for logging
            for i, item in enumerate(input_data)
        }

        processed_count = 0
        for future in concurrent.futures.as_completed(future_to_item_info):
            item_id_for_log = future_to_item_info[future]  # Get the identifier for logging
            processed_count += 1
            try:
                result_entry = future.result()  # Get the result from the completed future
                all_results.append(result_entry)
                if processed_count % 50 == 0 or processed_count == len(input_data):  # Log progress periodically
                    print(f"  Processed {processed_count}/{len(input_data)} items...")
            except FileNotFoundError as e:  # Critical error: bcftools not found
                print(f"CRITICAL ERROR: {e}. Aborting further processing.")
                # Attempt to cancel remaining futures (best-effort)
                for f_cancel in future_to_item_info:
                    if not f_cancel.done() and not f_cancel.cancelled():
                        f_cancel.cancel()
                break  # Stop processing if bcftools itself is not found
            except Exception as exc:
                # Catch any other unexpected exceptions from future.result() or process_item_for_bcftools
                print(f"An error occurred processing item associated with ID '{item_id_for_log}': {exc}")
                # Create a basic error entry for this item in the results
                all_results.append({
                    "haplotype_id": item_id_for_log,  # Or more context from original item if available
                    "query_region": "N/A - processing error",
                    "variants": [],
                    "bcftools_error": str(exc)
                })

    # --- Save all results to the output JSON file ---
    try:
        with open(args.output_json, 'w') as outfile:
            json.dump(all_results, outfile, indent=4)  # Use indent for pretty printing
        print(f"\nSuccessfully saved {len(all_results)} results to '{args.output_json}'")
    except IOError:
        print(f"Error: Could not write results to '{args.output_json}'.")
    except Exception as e:
        print(f"An unexpected error occurred while writing the output JSON: {e}")


if __name__ == "__main__":
    main()
