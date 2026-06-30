#!/usr/bin/env python3
import argparse
import pickle
import json

def main():
    """
    Converts a stats pickle file to either TSV or JSON format based on command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Convert a stats pickle file to TSV or JSON.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("pickle_file", help="Path to the input stats.pkl file.")
    parser.add_argument("output_file", help="Path for the output file.")
    parser.add_argument(
        "--format",
        choices=['tsv', 'json'],
        default='tsv',
        help="Output format."
    )
    args = parser.parse_args()

    # --- 1. Load data and apply filter ---
    print(f"Loading stats from {args.pickle_file}...")
    with open(args.pickle_file, 'rb') as fh:
        stats_data = pickle.load(fh)

    filtered_records = []
    for node_id, stats in stats_data.items():
        perfect = int(stats.get("perfect", 0))
        not_perfect = int(stats.get("not_perfect", 0))
        total_records = perfect + not_perfect

        # Filtering logic from the original script
        if total_records > 0 and not_perfect > 1 and (not_perfect / total_records) > 0.05:
            max_read = int(stats.get("max_read_length", 1) or 1)
            max_cigar = int(stats.get("max_cigar_length", 1) or 1)
            filtered_records.append({
                "node_id": int(node_id),
                "n_records": total_records,
                "max_read_len": max_read,
                "max_cigar_len": max_cigar
            })

    print(f"Found {len(filtered_records)} records passing the filter.")

    # --- 2. Write output based on selected format ---
    if args.format == 'tsv':
        with open(args.output_file, 'w') as fh:
            # Write header
            fh.write("node_id\tn_records\tmax_read_len\tmax_cigar_len\n")
            # Write data rows
            for record in filtered_records:
                fh.write(
                    f"{record['node_id']}\t"
                    f"{record['n_records']}\t"
                    f"{record['max_read_len']}\t"
                    f"{record['max_cigar_len']}\n"
                )

    elif args.format == 'json':
        with open(args.output_file, 'w') as fh:
            # indent makes the JSON file human-readable
            json.dump(filtered_records, fh, indent=2)

    print(f"Successfully generated {args.output_file}")

if __name__ == "__main__":
    main()