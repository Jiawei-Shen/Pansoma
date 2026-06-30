import subprocess
import json
import argparse
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def extract_nodes(graph_xg, path_name, output_json="nodes.json"):
    """Extract nodes with ID, sequence, and length from a given path."""
    path_vg = graph_xg.replace(".xg", ".vg")

    # Run vg find to extract nodes from the path
    subprocess.run(f"vg find -x {graph_xg} -p {path_name} > {path_vg}", shell=True, check=True)

    # Convert nodes to JSON format including sequence
    subprocess.run(f"vg view -v {path_vg} | jq -c '{{nodes: [.node[] | {{id: .id, sequence: .sequence}}]}}' > {output_json}", shell=True, check=True)

    print(f"[✔] Extracted node IDs, sequences, and saved to {output_json}")

def load_nodes(nodes_json):
    """Load node IDs, sequences, and compute their lengths."""
    with open(nodes_json, "r") as f:
        data = json.load(f)

    node_info = {}
    for node in data["nodes"]:
        node_id = str(node["id"])
        sequence = node["sequence"]
        node_info[node_id] = {"sequence": sequence, "length": len(sequence)}

    return node_info

def process_read(line, node_info, node_read_map, lock):
    """Check if the read aligns to any node and store it grouped by node."""
    try:
        read = json.loads(line)
        read_info = {
            "read_name": read["name"],
            "sequence": read["sequence"],
            "mapping_quality": read.get("mapping_quality", 0),
            "score": read.get("score", None),
            "quality": read.get("quality", ""),
            "path": read.get("path", {})
        }

        mapped_nodes = set()  # Store unique node IDs the read maps to
        for mapping in read.get("path", {}).get("mapping", []):
            node_id = str(mapping["position"].get("node_id", ""))
            if node_id in node_info:
                mapped_nodes.add(node_id)

        # Immediately store the read into the dictionary
        with lock:
            for node_id in mapped_nodes:
                if node_id not in node_read_map:
                    node_read_map[node_id] = {
                        "sequence": node_info[node_id]["sequence"],
                        "length": node_info[node_id]["length"],
                        "reads": []
                    }
                node_read_map[node_id]["reads"].append(read_info)

    except json.JSONDecodeError:
        return  # Skip invalid JSON reads

def filter_reads(input_gam, nodes_file, output_json, threads=4):
    """Filter reads from GAM that align to extracted nodes and group them by node in real-time."""
    # Load node information (id, sequence, length)
    node_info = load_nodes(nodes_file)

    # Shared data structure for storing grouped reads
    node_read_map = {}
    lock = threading.Lock()  # Ensures thread-safe writing to the dictionary

    # Start GAM file processing
    process = subprocess.Popen(["vg", "view", "-a", input_gam], stdout=subprocess.PIPE, text=True)

    processed_count = 0
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for _ in executor.map(lambda line: process_read(line, node_info, node_read_map, lock), process.stdout):
            processed_count += 1
            if processed_count % 5000 == 0:
                print(f"[INFO] Processed {processed_count} reads...")
                with lock:
                    save_json(node_read_map, output_json)

    # Final save
    save_json(node_read_map, output_json)
    print(f"[✔] Filtered reads grouped by node saved to {output_json}")
    return output_json

def save_json(data, output_file):
    """Saves the current JSON state to a file atomically."""
    temp_file = output_file + ".tmp"
    with open(temp_file, "w") as f:
        json.dump(data, f, indent=2)
    subprocess.run(["mv", temp_file, output_file])  # Atomic move
    print(f"[✔] Saved progress to {output_file}")

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description="Extract reads from a GAM file that align to a specific path and group them by node.")

    parser.add_argument("-x", "--xg", required=True, help="Graph XG file")
    parser.add_argument("-g", "--gam", required=True, help="Input GAM file")
    parser.add_argument("-p", "--path", required=True, help="Path name to filter reads")
    parser.add_argument("-n", "--nodes", default="hg38_chr5_nodes.json", help="Output JSON file for extracted node data")
    parser.add_argument("-j", "--json", default="grouped_reads_by_hg38_chr5_nodes.json", help="Output JSON file grouping reads by node")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for processing")

    args = parser.parse_args()

    # Extract node data (ID, sequence, length)
    extract_nodes(args.xg, args.path, args.nodes)

    # Process reads from GAM and group them by node
    filter_reads(args.gam, args.nodes, args.json, args.threads)

if __name__ == "__main__":
    main()
