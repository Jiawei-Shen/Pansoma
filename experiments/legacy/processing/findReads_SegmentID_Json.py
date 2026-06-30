import json
import argparse
import os
import subprocess
import threading
import queue

def extract_nodes(json_file_path):
    """Extracts node (segment) IDs for each genomic region from the given JSON file."""
    with open(json_file_path, "r") as f:
        data = json.load(f)

    node_dict = {}

    for region, details in data.items():
        if "segments" in details:
            node_ids = [segment["segment_id"] for segment in details["segments"]]
            node_dict[region] = node_ids

    return node_dict

def extract_reads_from_gam(gam_file, segment_ids):
    """Extracts reads from a GAM file that map to a given region's segment IDs."""
    result = subprocess.run(['vg', 'view', '-a', gam_file], stdout=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"Error running vg view on {gam_file}")
        return []

    matching_reads = set()

    for line in result.stdout.strip().split('\n'):
        read = json.loads(line)

        if 'path' in read:
            for mapping in read['path']['mapping']:
                if 'position' in mapping and 'node_id' in mapping['position']:
                    if str(mapping['position']['node_id']) in segment_ids:
                        matching_reads.add(read['name'])
                        break

    return matching_reads

def extract_reads_from_gaf(gaf_file, segment_ids):
    """Extracts reads from a GAF file that map to a given region's segment IDs."""
    matching_reads = set()

    with open(gaf_file, 'r') as f:
        for line in f:
            columns = line.strip().split("\t")
            if len(columns) > 5:
                path = columns[5]  # Segment/path information
                for node in segment_ids:
                    if node in path:
                        matching_reads.add(columns[0])  # Read name
                        break

    return matching_reads

def process_region(region, segment_ids, alignment_file, file_extension, result_queue):
    """Processes each region individually and stores results in a queue."""
    if file_extension == ".gam":
        reads = extract_reads_from_gam(alignment_file, segment_ids)
    elif file_extension == ".gaf":
        reads = extract_reads_from_gaf(alignment_file, segment_ids)
    else:
        raise ValueError("Unsupported file format. Please provide a GAM or GAF file.")

    result_queue.put((region, reads))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find reads mapped to nodes from a JSON file in a GAM or GAF file.")
    parser.add_argument("json_file", type=str, help="Input JSON file containing segment IDs.")
    parser.add_argument("alignment_file", type=str, help="Input GAM or GAF file for read extraction.")

    args = parser.parse_args()
    file_extension = os.path.splitext(args.alignment_file)[-1].lower()

    # Extract nodes from JSON
    nodes_dict = extract_nodes(args.json_file)
    print(f"Extracted {len(nodes_dict)} regions from {args.json_file}")

    # Initialize thread queue
    result_queue = queue.Queue()
    threads = []

    for region, segment_ids in nodes_dict.items():
        thread = threading.Thread(target=process_region, args=(region, segment_ids, args.alignment_file, file_extension, result_queue))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

    # Print results
    while not result_queue.empty():
        region, reads = result_queue.get()
        print(f"\nRegion: {region}")
        if reads:
            print(f"  Reads mapped ({len(reads)} total):")
            for read in reads:
                print(f"    {read}")
        else:
            print("  No reads found for this region.")
