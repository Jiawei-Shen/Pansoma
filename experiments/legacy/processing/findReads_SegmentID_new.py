import json
import argparse
import subprocess
import threading
import queue

def load_json(json_file):
    """Loads the JSON file and maps segment IDs to regions."""
    with open(json_file, "r") as f:
        data = json.load(f)

    segment_to_region = {}  # Map segment_id to region
    for region, details in data.items():
        if "segments" in details:
            for segment in details["segments"]:
                segment_id = str(segment["segment_id"])  # Convert to string for lookup
                segment_to_region[segment_id] = region

        # Initialize read storage for each region
        if "aligned_reads" not in details:
            details["aligned_reads"] = []

    return data, segment_to_region

def process_gam_worker(read_queue, segment_to_region, json_data, lock, processed_count):
    """Worker thread to process reads from the queue as soon as they arrive."""
    while True:
        try:
            read = read_queue.get(timeout=3)  # Get a read from the queue
        except queue.Empty:
            break  # Exit when no more reads are being added

        processed = False  # Track if the read was assigned to any region
        if 'path' in read:
            for mapping in read['path']['mapping']:
                if 'position' in mapping and 'node_id' in mapping['position']:
                    segment_id = str(mapping['position']['node_id'])

                    if segment_id in segment_to_region:
                        region = segment_to_region[segment_id]

                        # Check if 'offset' exists in position, otherwise default to 0
                        mapping_position = mapping["position"].get("offset", 0)

                        # Construct read details
                        read_info = {
                            "read_name": read['name'],
                            "sequence_length": len(read["sequence"]),
                            "mapping_position": mapping_position
                        }

                        # Lock required for safe multi-threaded dictionary updates
                        with lock:
                            json_data[region]["aligned_reads"].append(read_info)
                        processed = True
                        break  # Avoid duplicate reads in the same region

        # Update and print progress every 1000 reads
        if processed:
            with lock:
                processed_count[0] += 1
                if processed_count[0] % 1000 == 0:
                    print(f"Processed {processed_count[0]} reads...")
                # print(f"Processed {processed_count[0]} reads...")

        read_queue.task_done()


def process_gam_file(gam_file, segment_to_region, json_data, num_threads):
    """Reads GAM file, adds reads to queue, and processes them in parallel."""
    read_queue = queue.Queue()
    lock = threading.Lock()
    processed_count = [0]  # Shared counter for processed reads

    # Start worker threads **before** reading the GAM file
    threads = []
    for _ in range(num_threads):
        thread = threading.Thread(target=process_gam_worker, args=(read_queue, segment_to_region, json_data, lock, processed_count))
        thread.start()
        threads.append(thread)

    # Stream the GAM file as JSON **while workers are running**
    process = subprocess.Popen(['vg', 'view', '-a', gam_file], stdout=subprocess.PIPE, text=True)
    print("Started processing GAM file...")

    for line in process.stdout:
        try:
            read = json.loads(line.strip())  # Read and parse JSON line by line
            read_queue.put(read)  # Add to queue while workers are processing
        except json.JSONDecodeError as e:
            print(f"⚠️ Skipping invalid JSON read: {e}")

    # Wait for all threads to finish
    read_queue.join()
    for thread in threads:
        thread.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find reads mapped to segments from a JSON file and update the JSON with aligned reads.")
    parser.add_argument("json_file", type=str, help="Input JSON file containing segment IDs.")
    parser.add_argument("gam_file", type=str, help="Input GAM file.")
    parser.add_argument("output_json", type=str, help="Output JSON file with updated aligned reads.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4)")

    args = parser.parse_args()

    # Load JSON and segment-region mapping
    json_data, segment_to_region = load_json(args.json_file)
    print(f"Loaded {len(segment_to_region)} segment IDs from {args.json_file}")

    # Process the GAM file with user-defined thread count
    process_gam_file(args.gam_file, segment_to_region, json_data, args.threads)

    # Save to new JSON file
    with open(args.output_json, "w") as f:
        json.dump(json_data, f, indent=4)

    print(f"\nUpdated JSON saved to {args.output_json}")
