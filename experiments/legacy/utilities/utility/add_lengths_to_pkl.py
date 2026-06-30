#!/usr/bin/env python3
import argparse
import pickle
import gzip

def parse_gfa_lengths(gfa_path):
    """
    Parse GFA v1 segment lines and return {node_id:int -> length:int}.
    Prefer LN:i:<len> tag; fallback = len(sequence).
    """
    node_lengths = {}
    open_func = gzip.open if gfa_path.endswith(".gz") else open
    with open_func(gfa_path, "rt") as f:
        for line in f:
            if not line.startswith("S"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[1]
            seq = parts[2]

            # Extract length
            length = None
            for field in parts[3:]:
                if field.startswith("LN:i:"):
                    try:
                        length = int(field.split(":", 2)[2])
                    except ValueError:
                        pass
                    break
            if length is None and seq != "*":
                length = len(seq)
            if length is None:
                length = 0

            try:
                node_id = int(name)
            except ValueError:
                continue  # skip non-numeric node names

            node_lengths[node_id] = length
    return node_lengths

def main():
    parser = argparse.ArgumentParser(description="Add node lengths from GFA to an existing PKL file.")
    parser.add_argument("input_pkl", help="Input PKL file (with perfect/not_perfect counts)")
    parser.add_argument("gfa", help="GFA file (can be .gz)")
    parser.add_argument("output_pkl", help="Output PKL file with lengths added")
    args = parser.parse_args()

    # Load existing PKL
    with open(args.input_pkl, "rb") as f:
        data = pickle.load(f)

    # Parse GFA lengths
    lengths = parse_gfa_lengths(args.gfa)

    # Merge
    missing = 0
    for nid, rec in data.items():
        ln = lengths.get(nid, 0)
        rec["length"] = ln
        if ln == 0 and nid not in lengths:
            missing += 1

    print(f"Updated {len(data)} nodes. Missing lengths for {missing} nodes.")

    # Save updated PKL
    with open(args.output_pkl, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    main()
