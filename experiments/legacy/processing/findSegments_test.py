import json
import subprocess
import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


def extract_segments_from_json(gbz_file, path_range):
    # Run vg find and pipe to vg view to get JSON output
    cmd = f"vg find -x {gbz_file} -p {path_range} | vg view -v -j -"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    segment_info = {}

    for line in process.stdout:
        data = json.loads(line)

        # Extract node (segment) IDs
        if 'node' in data:
            for node in data['node']:
                node_id = node['id']

                # Check if this node is part of a path
                if 'path' in data:
                    for path in data['path']:
                        for mapping in path['mapping']:
                            if mapping['position']['node_id'] == node_id:
                                if node_id not in segment_info:
                                    segment_info[node_id] = {'segment_id': node_id, 'paths': []}
                                if path['name'] not in segment_info[node_id]['paths']:
                                    segment_info[node_id]['paths'].append(path['name'])

    process.stdout.close()
    process.wait()

    return list(segment_info.values())


def call_find_segments(gbz_file, chromosome, position, record, windowLength=100, ref="GRCh38#0"):
    """Call extract_segments_from_json with the given chromosome and position."""
    formatted_position = f"{ref}#{chromosome}:{position - windowLength}-{position + windowLength}"
    segments = extract_segments_from_json(gbz_file, formatted_position)
    print(f"Output for {formatted_position}: {segments}")
    return {formatted_position: {"segments": segments, "variants": record}}


def process_record(gbz_file, record):
    chromosome = record.get('chromosome')
    position = record.get('position')

    if not chromosome or not position:
        return None, f"Skipping record due to missing chromosome or position: {record}"

    print(f"Processing {chromosome}:{position}")
    segment_data = call_find_segments(gbz_file, chromosome, position, record)
    return segment_data, None


def main(json_file, gbz_file, threads, output_file):
    # Load JSON data
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)

    results = {}

    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_record = {executor.submit(process_record, gbz_file, record): record for record in data}

        for i, future in enumerate(as_completed(future_to_record)):
            record = future_to_record[future]
            try:
                segment_data, error = future.result()
                if error:
                    print(error)
                    continue
                results.update(segment_data)
            except Exception as e:
                print(f"Error processing record {record}: {e}")

            # Save results every 80 iterations
            if (i + 1) % 80 == 0:
                try:
                    with open(output_file, 'w') as outfile:
                        json.dump(results, outfile, indent=4)
                    print(f"Intermediate results saved to {output_file} after {i + 1} records.")
                except Exception as e:
                    print(f"Error saving intermediate results: {e}")

    # Save final results
    try:
        with open(output_file, 'w') as outfile:
            json.dump(results, outfile, indent=4)
        print(f"Final results saved to {output_file}.")
    except Exception as e:
        print(f"Error saving final results: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process JSON file to extract segment IDs and paths.')
    parser.add_argument('-j', '--json_file', type=str,
                        help='Path to the JSON file containing chromosome and position data.')
    parser.add_argument('-x', '--gbz_file', type=str, required=True, help='Path to the GBZ index or graph file.')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads to use for parallel processing.')
    parser.add_argument('-o', '--output_path', default="output_seg.json", type=str,
                        help='Path to the output file containing segment and position data.')

    args = parser.parse_args()
    main(args.json_file, args.gbz_file, args.threads, output_file=args.output_path)
