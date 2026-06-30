import subprocess
import json
import argparse


def extract_reads(gam_file, segment_ids):
    # Run vg view to convert the GAM file to JSON
    result = subprocess.run(['vg', 'view', '-a', gam_file], stdout=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"Error running vg view on {gam_file}")
        return

    matching_reads = []

    # Process each read in the GAM file
    for line in result.stdout.strip().split('\n'):
        read = json.loads(line)

        # Check if any of the path segments match the provided segment IDs
        if 'path' in read:
            for mapping in read['path']['mapping']:
                if 'position' in mapping and 'node_id' in mapping['position']:
                    if str(mapping['position']['node_id']) in segment_ids:
                        matching_reads.append(read['name'])
                        break

    return matching_reads


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find reads in a GAM file that contain specific segment IDs.")
    parser.add_argument("gam_file", help="Path to the input GAM file.")
    parser.add_argument("segment_ids", nargs='+', help="List of segment IDs to search for.")
    args = parser.parse_args()

    reads = extract_reads(args.gam_file, args.segment_ids)

    if reads:
        print(f"Reads containing specified segments ({', '.join(args.segment_ids)}):")
        for read in reads:
            print(read)
    else:
        print("No reads found containing the specified segment IDs.")
