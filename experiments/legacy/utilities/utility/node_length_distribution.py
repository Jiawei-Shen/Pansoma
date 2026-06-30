import json
import argparse
import collections
import statistics
import math
import sys


def analyze_node_lengths(json_filepath, num_bins_requested=10, top_n_freq=10):
    """
    Analyzes and prints the distribution of node lengths from a JSON file.
    The JSON file is expected to have a top-level "nodes" key containing a list
    of objects, each with a "length" field.
    """
    try:
        with open(json_filepath, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: JSON file not found at {json_filepath}", file=sys.stderr)
        return
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_filepath}. Check for syntax errors.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An unexpected error occurred while reading the JSON file: {e}", file=sys.stderr)
        return

    if not isinstance(data, dict):
        print("Error: Root JSON content is not an object/dictionary.", file=sys.stderr)
        return

    node_list = data.get("nodes")
    if not isinstance(node_list, list):
        print(f"Error: The key 'nodes' was not found in '{json_filepath}', or its value is not a list.",
              file=sys.stderr)
        return

    print(f"Found {len(node_list)} node entries in the 'nodes' list of '{json_filepath}'.")

    lengths = []
    nodes_with_invalid_format_or_length = 0
    for i, node_obj in enumerate(node_list):
        if not isinstance(node_obj, dict):
            nodes_with_invalid_format_or_length += 1
            continue

        length_val = node_obj.get("length")
        if length_val is None:
            nodes_with_invalid_format_or_length += 1
            continue

        try:
            lengths.append(int(length_val))
        except (ValueError, TypeError):
            nodes_with_invalid_format_or_length += 1
            continue

    if nodes_with_invalid_format_or_length > 0:
        print(
            f"Note: {nodes_with_invalid_format_or_length} node entries were skipped due to being non-dictionary items, or missing/non-numeric 'length' fields.")

    if not lengths:
        print("No valid node lengths found to analyze.")
        return

    print(f"Successfully extracted {len(lengths)} valid node lengths.")

    # --- Summary Statistics ---
    print("\nSummary Statistics for Node Lengths:")
    print("-----------------------------------")
    min_len = min(lengths)
    max_len = max(lengths)
    mean_len = statistics.mean(lengths)
    median_len = statistics.median(lengths)
    print(f"{'Count:':<15} {len(lengths)}")
    print(f"{'Minimum:':<15} {min_len}")
    print(f"{'Maximum:':<15} {max_len}")
    print(f"{'Mean:':<15} {mean_len:.2f}")
    print(f"{'Median:':<15} {median_len}")
    if len(lengths) > 1:
        stdev_len = statistics.stdev(lengths)
        print(f"{'Std Dev:':<15} {stdev_len:.2f}")
    else:
        print(f"{'Std Dev:':<15} N/A (single data point)")
    print("-----------------------------------")

    # --- Top N Most Frequent Lengths ---
    if top_n_freq > 0:
        print(f"\nTop {top_n_freq} Most Frequent Node Lengths:")
        print("--------------------------------------")
        print(f"{'Length':<10} | {'Count':<10} | {'Percentage':<10}")
        print("--------------------------------------")
        length_counts_counter = collections.Counter(lengths)
        total_valid_lengths = len(lengths)
        for length, count in length_counts_counter.most_common(top_n_freq):
            percentage = (count / total_valid_lengths) * 100 if total_valid_lengths > 0 else 0
            print(f"{length:<10} | {count:<10} | {percentage:>9.2f}%")
        print("--------------------------------------")

    # --- Binned Distribution ---
    print("\nLength Distribution (Binned):")
    if max_len == min_len:  # All lengths are the same
        print("------------------------------------------------------")
        print(f"{'Bin Range (inclusive)':<25} | {'Count':<10} | {'Percentage':<10}")
        print("------------------------------------------------------")
        print(f"{min_len:<25} | {len(lengths):<10} | {100.0:>9.2f}%")
        print("------------------------------------------------------")
    else:
        actual_num_bins = max(1, num_bins_requested)  # Ensure at least one bin

        # Calculate bin_width, ensuring it's at least 1.
        # This width aims to distribute the range (max_len - min_len + 1) across actual_num_bins.
        data_span = max_len - min_len + 1
        bin_width = math.ceil(data_span / actual_num_bins)
        bin_width = max(1, bin_width)  # Ensure bin_width is at least 1

        # Determine the number of bins we will actually use based on this width
        num_effective_bins = math.ceil(data_span / bin_width)

        bin_counts_list = [0] * num_effective_bins

        for length_val in lengths:
            # Determine which bin length_val falls into
            # The index is based on how many bin_widths away the length_val is from min_len
            bin_index = math.floor((length_val - min_len) / bin_width)
            # Clamp bin_index to be within the valid range [0, num_effective_bins - 1]
            bin_index = max(0, min(bin_index, num_effective_bins - 1))
            bin_counts_list[bin_index] += 1

        print(f"(Using {num_effective_bins} bins, approximate width {bin_width})")
        print("------------------------------------------------------")
        print(f"{'Bin Range (inclusive)':<25} | {'Count':<10} | {'Percentage':<10}")
        print("------------------------------------------------------")
        total_valid_lengths = len(lengths)
        for i in range(num_effective_bins):
            bin_start_val = min_len + i * bin_width
            # The end of the bin is one less than the start of the next bin
            bin_end_val = bin_start_val + bin_width - 1

            # Ensure the last bin's displayed end does not exceed max_len
            # and the bin_start_val itself does not exceed max_len (already handled by num_effective_bins effectively)
            bin_end_val = min(bin_end_val, max_len)

            range_str = f"{bin_start_val} - {bin_end_val}"
            count = bin_counts_list[i]
            percentage = (count / total_valid_lengths) * 100 if total_valid_lengths > 0 else 0
            print(f"{range_str:<25} | {count:<10} | {percentage:>9.2f}%")
        print("------------------------------------------------------")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate and display the distribution of node lengths from a JSON file. "
                    "The JSON file must contain a top-level 'nodes' key, which is a list of objects, "
                    "each object having a 'length' field."
    )
    parser.add_argument("json_filepath", help="Path to the input JSON file.")
    parser.add_argument(
        "--num_bins",
        type=int,
        default=10,
        help="Approximate number of bins for the length distribution histogram (default: 10)."
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=10,
        help="Number of most frequent lengths to display (default: 10; use 0 to disable)."
    )

    args = parser.parse_args()

    analyze_node_lengths(args.json_filepath, args.num_bins, args.top_n)


if __name__ == "__main__":
    main()