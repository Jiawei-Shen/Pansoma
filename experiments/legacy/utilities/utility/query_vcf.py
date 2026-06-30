import json
import subprocess
import shlex  # For safely formatting command arguments for display
import argparse
import os
import concurrent.futures  # For parallel processing


def execute_bcftools_query(vcf_file_path, region_string):
    """
    Executes a bcftools view query for a given region, excluding header lines.

    Args:
        vcf_file_path (str): Path to the VCF file.
        region_string (str): The region to query (e.g., "chr1:100-200").

    Returns:
        tuple: (list_of_variant_strings, bcftools_stderr_str)
               Returns (None, error_message) if bcftools command itself fails.
    """
    # Command to execute: bcftools view -H <file.vcf.gz> <chr:from-to>
    # The -H or --no-header option suppresses VCF header output
    command = ["bcftools", "view", "-H", vcf_file_path, region_string]
    # Worker PID is useful for debugging parallel processes
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
            # Since -H is used, no header lines should be present in stdout
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
        item (dict): A dictionary from the input JSON.
        vcf_file_path (str): Path to the VCF file.

    Returns:
        dict: A dictionary containing original item info, the query region string,
              a list of found variants, and an internal 'bcftools_error' field.
    """
    base_result_entry = {
        "top_level_node_id": item.get("top_level_node_id"),
        "haplotype_id": item.get("haplotype_id"),
        "node_sequence_length": item.get("node_sequence_length"),
        "grch38_chromosome_region": item.get("grch38_chromosome_region"),
        "grch38_start_position": item.get("grch38_start_position"),
        "query_region": None,
        "variants": [],
        "bcftools_error": None  # Internal field for error tracking
    }

    start_pos_str = item.get("grch38_start_position")
    seq_len_val = item.get("node_sequence_length")
    chromosome = item.get("grch38_chromosome_region")
    item_hap_id_for_log = base_result_entry.get('haplotype_id', 'N/A')

    if not chromosome:
        base_result_entry["bcftools_error"] = "Skipped: Missing 'grch38_chromosome_region'."
        print(f"  Skipping item (HapID: {item_hap_id_for_log}): Missing 'grch38_chromosome_region'.")
        return base_result_entry

    if start_pos_str is None:
        base_result_entry["bcftools_error"] = "Skipped: 'grch38_start_position' is null."
        print(f"  Skipping item (HapID: {item_hap_id_for_log}): 'grch38_start_position' is null.")
        return base_result_entry

    if seq_len_val is None:
        base_result_entry["bcftools_error"] = "Skipped: 'node_sequence_length' is null."
        print(f"  Skipping item (HapID: {item_hap_id_for_log}): 'node_sequence_length' is null.")
        return base_result_entry

    try:
        start_pos = int(start_pos_str)
        seq_len = int(seq_len_val)
        if seq_len < 0:
            base_result_entry[
                "bcftools_error"] = f"Skipped: Invalid 'node_sequence_length' ({seq_len}). Must be non-negative."
            print(f"  Skipping item (HapID: {item_hap_id_for_log}): Invalid 'node_sequence_length' ({seq_len}).")
            return base_result_entry
    except ValueError:
        base_result_entry[
            "bcftools_error"] = "Skipped: 'grch38_start_position' or 'node_sequence_length' is not a valid integer."
        print(
            f"  Skipping item (HapID: {item_hap_id_for_log}): 'grch38_start_position' or 'node_sequence_length' not a valid int.")
        return base_result_entry

    end_pos = start_pos + seq_len

    if end_pos < start_pos:
        base_result_entry[
            "bcftools_error"] = f"Skipped: Calculated end position ({end_pos}) is less than start position ({start_pos}) with length ({seq_len})."
        print(
            f"  Skipping item (HapID: {item_hap_id_for_log}): Calculated end_pos {end_pos} < start_pos {start_pos} with length {seq_len}.")
        return base_result_entry

    query_region_str = f"{chromosome}:{start_pos}-{end_pos}"
    base_result_entry["query_region"] = query_region_str

    variants_list, bcftools_stderr_output = execute_bcftools_query(vcf_file_path, query_region_str)

    if variants_list is not None:
        base_result_entry["variants"] = variants_list

    if bcftools_stderr_output:
        if "bcftools command not found" in bcftools_stderr_output:
            raise FileNotFoundError(bcftools_stderr_output)
        # Store bcftools stderr output. If there was a pre-skip error, append.
        if base_result_entry["bcftools_error"]:  # e.g., from pre-skip
            base_result_entry["bcftools_error"] += f"; bcftools stderr: {bcftools_stderr_output}"
        else:
            base_result_entry["bcftools_error"] = f"bcftools stderr: {bcftools_stderr_output}"

    return base_result_entry


def main():
    """
    Main function to parse arguments and run the VCF querying process.
    """
    parser = argparse.ArgumentParser(
        description="Query a VCF file using bcftools view based on regions from an input JSON file. "
                    "The 'bcftools_error' field will not be included in the final output JSON. "
                    "Use --save-empty-records to include records with no variants (but no errors).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
        default="bcftools_query_results.json",
        help="Path to save the output JSON file."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=None,
        help=(f"Number of parallel worker threads to use for bcftools queries. "
              f"(default: number of CPU cores, or 1 if CPU count cannot be determined. "
              f"Currently detected cores: {os.cpu_count() if os.cpu_count() else 'N/A'})")
    )
    parser.add_argument(
        "--save-empty-records",
        action="store_true",
        default=False,
        help="If specified, save records to output even if no variants were found (as long as there was no bcftools error during query)."
    )

    args = parser.parse_args()

    max_workers = args.threads if args.threads is not None else (os.cpu_count() if os.cpu_count() else 1)
    if max_workers <= 0:
        print("Warning: Number of threads must be positive. Defaulting to 1.")
        max_workers = 1

    print("--- Effective Configuration ---")
    print(f"Input JSON: {args.input_json}")
    print(f"VCF File: {args.vcf_file}")
    print(f"Output JSON: {args.output_json}")
    print(f"Max Workers (Threads): {max_workers}")
    print(f"Save Empty Records: {args.save_empty_records}")
    print("-----------------------------")

    if not os.path.exists(args.input_json):
        print(f"Error: Input JSON file not found at '{args.input_json}'")
        return
    if not os.path.exists(args.vcf_file):
        print(f"Error: VCF file not found at '{args.vcf_file}'")
        return
    if not (os.path.exists(args.vcf_file + ".tbi") or os.path.exists(args.vcf_file + ".csi")):
        print(f"Error: Index file (.tbi or .csi) not found for '{args.vcf_file}'. Please index the VCF file.")
        return

    try:
        with open(args.input_json, 'r') as f:
            input_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{args.input_json}'. Make sure it's a valid JSON file.")
        return

    if not isinstance(input_data, list):
        print(f"Error: Input JSON is not a list as expected. Found type: {type(input_data)}")
        return

    final_results_to_save = []

    print(f"\nProcessing {len(input_data)} items from '{args.input_json}'...")
    print(f"Querying VCF: '{args.vcf_file}' with bcftools view (no header).")

    items_to_process = []
    skipped_due_to_input_validation = 0

    for i, item in enumerate(input_data):
        item_hap_id_for_log = item.get('haplotype_id', f'index_{i}')
        start_pos_str = item.get("grch38_start_position")
        seq_len_val = item.get("node_sequence_length")
        chromosome = item.get("grch38_chromosome_region")

        if not chromosome or start_pos_str is None or seq_len_val is None:
            print(
                f"  [PRE-SKIP] Item (HapID: {item_hap_id_for_log}): Missing one or more required fields (chromosome, start_position, node_sequence_length).")
            skipped_due_to_input_validation += 1
            continue
        try:
            int(start_pos_str)
            if int(seq_len_val) < 0:
                print(
                    f"  [PRE-SKIP] Item (HapID: {item_hap_id_for_log}): Invalid 'node_sequence_length' ({seq_len_val}).")
                skipped_due_to_input_validation += 1
                continue
        except ValueError:
            print(
                f"  [PRE-SKIP] Item (HapID: {item_hap_id_for_log}): 'grch38_start_position' or 'node_sequence_length' not a valid int.")
            skipped_due_to_input_validation += 1
            continue
        items_to_process.append(item)

    if not items_to_process:
        print("No valid items to process after initial input validation.")
        if skipped_due_to_input_validation > 0:
            print(f"Skipped {skipped_due_to_input_validation} items due to missing/invalid fields in input JSON.")
        return

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_item_info = {
            executor.submit(process_item_for_bcftools, item, args.vcf_file):
                item.get('haplotype_id', f'original_index_{input_data.index(item) if item in input_data else "N/A"}')
            for item in items_to_process
        }

        processed_count = 0
        skipped_from_output_count = 0
        for future in concurrent.futures.as_completed(future_to_item_info):
            item_id_for_log = future_to_item_info[future]
            processed_count += 1
            try:
                result_entry = future.result()

                # Condition for saving to output file
                save_this_record = False
                # Check if there was a bcftools error reported by the worker or a pre-skip error
                # The process_item_for_bcftools function populates result_entry["bcftools_error"]
                # if any issue occurs, including pre-skips or actual bcftools stderr.
                actual_bcftools_error = result_entry.get("bcftools_error")

                if actual_bcftools_error is None:  # Only proceed if the query and pre-check were clean
                    if result_entry.get("variants"):  # Variants were found
                        save_this_record = True
                    elif args.save_empty_records:  # No variants, but user wants to save empty records
                        save_this_record = True

                if save_this_record:
                    # Create a new dict for output, excluding the internal "bcftools_error" field
                    output_record = {k: v for k, v in result_entry.items() if k != "bcftools_error"}
                    final_results_to_save.append(output_record)
                else:
                    skipped_from_output_count += 1
                    # Log why it was skipped from output
                    query_region_log = result_entry.get('query_region', 'N/A')
                    if actual_bcftools_error:
                        # Don't re-log simple "Skipped:" messages from pre-validation if they somehow got here
                        if not actual_bcftools_error.startswith("Skipped:"):
                            print(
                                f"  Item (HapID: {item_id_for_log}): Query for region '{query_region_log}' had issues. Error: '{actual_bcftools_error}'. Not saved.")
                        # else: pre-skip error was already logged by process_item_for_bcftools
                    elif not result_entry.get("variants") and not args.save_empty_records:
                        print(
                            f"  Item (HapID: {item_id_for_log}): No variants found in region '{query_region_log}'. Not saved (as per --save-empty-records=False).")
                    # else: Other reasons for not saving (e.g. variants empty and save_empty_records is False)

                if processed_count % 50 == 0 or processed_count == len(items_to_process):
                    print(f"  Processed {processed_count}/{len(items_to_process)} submitted tasks...")

            except FileNotFoundError as e:
                print(f"CRITICAL ERROR: {e}. Aborting further processing.")
                for f_cancel in future_to_item_info:
                    if not f_cancel.done() and not f_cancel.cancelled():
                        f_cancel.cancel()
                break
            except Exception as exc:
                print(f"An error occurred processing future for item associated with ID '{item_id_for_log}': {exc}")
                skipped_from_output_count += 1

    print(f"\n--- Summary ---")
    print(f"Total items in input JSON: {len(input_data)}")
    print(f"Items skipped due to initial validation: {skipped_due_to_input_validation}")
    print(f"Items submitted for bcftools query: {len(items_to_process)}")
    print(f"Items processed by workers: {processed_count}")
    print(f"Records saved to output: {len(final_results_to_save)}")
    print(f"Records not saved to output (no variants/error/skipped by flag): {skipped_from_output_count}")

    try:
        with open(args.output_json, 'w') as outfile:
            json.dump(final_results_to_save, outfile, indent=4)
        print(f"\nSuccessfully saved {len(final_results_to_save)} records to '{args.output_json}'")
    except IOError:
        print(f"Error: Could not write results to '{args.output_json}'.")
    except Exception as e:
        print(f"An unexpected error occurred while writing the output JSON: {e}")


if __name__ == "__main__":
    main()
