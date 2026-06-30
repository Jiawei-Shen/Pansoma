import json
import argparse
import threading
import queue
import vg_pb2  # Import VG's Protobuf schema

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
    """Worker thread to process reads from the queue."""
    while True:
        try:
            read = read_queue.get(timeout=3)  # Get a read from the queue
        except queue.Empty:
            break

        processed = False
        for mapping in read.path.mapping:
            segment_id = str(mapping.position.node_id)

            if segment_id in segment_to_region:
                region = segment_to_region[segment_id]

                # Construct read details
                read_info = {
                    "read_name": read.name,
                    "sequence_length": len(read.sequence),
                    "mapping_position": mapping.position.offset
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

        read_queue.task_done()

def iterate_gam_binary(gam_file, read_queue):
    """Reads a GAM file in binary format and adds reads to a queue."""
    with open(gam_file, "rb") as f:
        while True:
            try:
                read = vg_pb2.Alignment()

                # Read the size of the next message (4 bytes)
                size_bytes = f.read(4)
                if not size_bytes:
                    break  # End of file

                # Convert the 4-byte size header to an integer
                size = int.from_bytes(size_bytes, "little")

                # Read the exact number of bytes expected for the message
                message_data = f.read(size)
                if len(message_data) != size:
                    print(f"⚠️ Skipping corrupted read (expected {size} bytes, got {len(message_data)})")
                    continue  # Skip to the next read

                # Parse the protobuf message
                read.ParseFromString(message_data)
                read_queue.put(read)  # Add valid read to the queue

            except Exception as e:
                print(f"⚠️ Skipping a corrupted read: {e}")
                continue  # Skip the problematic read and move to the next one


def process_gam_file(gam_file, segment_to_region, json_data, num_threads):
    """Distributes binary GAM parsing across multiple threads."""
    read_queue = queue.Queue()
    lock = threading.Lock()
    processed_count = [0]  # Shared counter for processed reads

    # Start a thread to read the GAM file
    reader_thread = threading.Thread(target=iterate_gam_binary, args=(gam_file, read_queue))
    reader_thread.start()

    # Start worker threads
    threads = []
    for _ in range(num_threads):
        thread = threading.Thread(target=process_gam_worker, args=(read_queue, segment_to_region, json_data, lock, processed_count))
        thread.start()
        threads.append(thread)

    # Wait for all threads to finish
    reader_thread.join()
    read_queue.join()
    for thread in threads:
        thread.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find reads mapped to segments from a JSON file and update the JSON with aligned reads.")
    parser.add_argument("json_file", type=str, help="Input JSON file containing segment IDs.")
    parser.add_argument("gam_file", type=str, help="Input GAM file (binary).")
    parser.add_argument("output_json", type=str, help="Output JSON file with updated aligned reads.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")

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
