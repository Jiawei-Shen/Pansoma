#!/usr/bin/env python3
import struct
import argparse
import subprocess

def extract_node_ids(idx_path, out_txt):
    node_ids = []
    with open(idx_path, "rb") as f:
        # Read block count (first 4 bytes)
        block_count_bytes = f.read(4)
        if not block_count_bytes or len(block_count_bytes) < 4:
            raise RuntimeError("Invalid .idx file: cannot read block count")
        block_count = struct.unpack("<I", block_count_bytes)[0]

        # Each entry = <I Q I I H> (22 bytes)
        entry_struct = struct.Struct("<I Q I I H")

        for _ in range(block_count):
            entry_data = f.read(entry_struct.size)
            if len(entry_data) < entry_struct.size:
                raise RuntimeError("Unexpected EOF in .idx file")
            node_id, _, _, _, _ = entry_struct.unpack(entry_data)
            node_ids.append(node_id)

    with open(out_txt, "w") as out_f:
        for nid in node_ids:
            out_f.write(f"{nid}\n")

    print(f"Extracted {len(node_ids)} node IDs to {out_txt}")
    return out_txt

def main():
    parser = argparse.ArgumentParser(description="Extract node IDs from .idx and run gbztool find-batch.")
    parser.add_argument("idx_path", help="Path to .idx file")
    parser.add_argument("gbz_path", help="Path to .gbz file")
    parser.add_argument("json_out", help="Output JSON file")
    parser.add_argument("--tmp_txt", default="tmp_node_id.txt", help="Temporary node ID list file")
    parser.add_argument("--gbztool", default="/scratch/xsong/project/gbz-tool/gbztool",
                        help="Path to gbztool executable")
    args = parser.parse_args()

    txt_path = extract_node_ids(args.idx_path, args.tmp_txt)

    # Run gbztool
    cmd = [args.gbztool, "find-batch", args.gbz_path, txt_path, args.json_out]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
