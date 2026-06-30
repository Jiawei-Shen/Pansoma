#!/usr/bin/env python3
"""
Aggregate perfect / un‑perfect read‑segments per VG node
with strand, and save as JSON / Pickle / both.
"""

import argparse, gzip, json, pickle, time, base64
import concurrent.futures, vg_pb2

# ─────────────────────────── helpers ─────────────────────────────────────────
def read_varint(stream):
    value = shift = 0
    while True:
        b = stream.read(1)
        if not b: raise EOFError
        byte = b[0]
        value |= (byte & 0x7F) << shift
        if not (byte & 0x80): return value
        shift += 7

def file_is_gzip(path: str) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

def parse_gam_groups(path: str, tag="GAM"):
    opener = gzip.open if file_is_gzip(path) else open
    with opener(path, "rb") as fh:
        while True:
            try:
                group_count = read_varint(fh)
            except EOFError:
                break
            if group_count == 0:
                continue
            try:
                tag_len = read_varint(fh)
                group_tag = fh.read(tag_len).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if group_tag != tag:
                for _ in range(group_count - 1):      # skip whole group
                    try:
                        skip = read_varint(fh)
                        fh.seek(skip, 1)
                    except EOFError:
                        break
                continue
            messages = []
            for _ in range(group_count - 1):
                try:
                    size = read_varint(fh)
                    msg = fh.read(size)
                except EOFError:
                    break
                if len(msg) == size:
                    messages.append(msg)
            if messages:
                yield messages

# ─────────────────────────── worker ──────────────────────────────────────────
def process_group(args):
    message_list, wanted_nodes, chrom_filter = args
    node_segments, read_count = {}, 0
    for raw in message_list:
        aln = vg_pb2.Alignment()
        aln.ParseFromString(raw)

        if chrom_filter and not any(rp.name == chrom_filter for rp in aln.refpos):
            continue

        read_count += 1
        seq = aln.sequence
        qual_bytes = aln.quality
        mapq = aln.mapping_quality
        read_offset = 0

        for mapping in aln.path.mapping:
            if not mapping.position.node_id:
                continue
            node_id = mapping.position.node_id
            if node_id not in wanted_nodes:
                for e in mapping.edit:
                    read_offset += max(e.from_length, len(e.sequence))
                continue
            node_offset = mapping.position.offset
            strand_char = "-" if mapping.position.is_reverse else "+"

            part_seq, part_q = [], bytearray()
            for edit in mapping.edit:
                if edit.from_length:
                    part_seq.append(
                        seq[read_offset : read_offset + edit.from_length].lower()
                        if edit.sequence else
                        seq[read_offset : read_offset + edit.from_length]
                    )
                    part_q.extend(
                        qual_bytes[read_offset : read_offset + edit.from_length]
                    )
                    read_offset += edit.from_length
                elif edit.sequence:   # insertion
                    ins = len(edit.sequence)
                    part_seq.append(seq[read_offset : read_offset + ins].lower())
                    part_q.extend(qual_bytes[read_offset : read_offset + ins])
                    read_offset += ins

            node_segments.setdefault(node_id, []).append({
                "offset": node_offset,
                "sequence": "".join(part_seq),
                "base_quality": base64.b64encode(part_q).decode(),
                "read_quality": mapq,
                "strand": strand_char
            })

    return node_segments, read_count

def merge_partial(parts):
    merged = {}; total = 0
    for part, cnt in parts:
        total += cnt
        for nid, segs in part.items():
            merged.setdefault(nid, []).extend(segs)
    return merged, total

# ─────────────────────────── pipeline ────────────────────────────────────────
def run_pipeline(
    gam_path, stats_pkl, prefix, fmt, threads, max_pending, milestone, chrom
):
    print(f"Loading stats: {stats_pkl}")
    with open(stats_pkl, "rb") as fh:
        stats = pickle.load(fh)

    total_nodes = len(stats)
    nodes_with_unperfect = sum(1 for s in stats.values() if s["not_perfect"] > 0)
    filtered_nodes = {
        int(nid)
        for nid, s in stats.items()
        if s["not_perfect"] > 1 and s["not_perfect"] / (s["perfect"] + s["not_perfect"]) > 0.1
    }

    print("\nNode‑level overview")
    print(f"  Total nodes               : {total_nodes}")
    print(f"  Nodes with ≥1 un‑perfect  : {nodes_with_unperfect} "
          f"({nodes_with_unperfect/total_nodes*100:.2f} %)")
    print(f"  Nodes passing filter      : {len(filtered_nodes)} "
          f"({len(filtered_nodes)/total_nodes*100:.2f} %)\n")

    partials = []; reads_total = 0; next_milestone = milestone
    start = time.perf_counter()

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        pending = []
        for batch in parse_gam_groups(gam_path):
            pending.append(ex.submit(process_group, (batch, filtered_nodes, chrom)))
            while len(pending) >= max_pending:
                done, not_done = concurrent.futures.wait(
                    pending, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    partials.append(fut.result())
                    reads_total += fut.result()[1]
                pending = list(not_done)
                if reads_total >= next_milestone:
                    print(f"Milestone {reads_total} reads | "
                          f"{time.perf_counter()-start:.1f}s")
                    next_milestone += milestone

        for fut in concurrent.futures.as_completed(pending):
            partials.append(fut.result()); reads_total += fut.result()[1]
            if reads_total >= next_milestone:
                print(f"Milestone {reads_total} reads | "
                      f"{time.perf_counter()-start:.1f}s")
                next_milestone += milestone

    merged, _ = merge_partial(partials)

    # ----- pre‑write summary -----
    print("\n--- pre‑write summary ---")
    print(f"  Reads processed : {reads_total}")
    print(f"  Nodes saved     : {len(merged)}")
    print(f"  Elapsed so far  : {time.perf_counter()-start:.2f}s")

    if fmt in ("json", "both"):
        jfile = prefix + ".json"
        with open(jfile, "w") as jf:
            json.dump({str(k): v for k, v in merged.items()}, jf, indent=2)
        print(f"  JSON written    : {jfile}")
    if fmt in ("pkl", "both"):
        pfile = prefix + ".pkl"
        with open(pfile, "wb") as pf:
            pickle.dump(merged, pf, pickle.HIGHEST_PROTOCOL)
        print(f"  Pickle written  : {pfile}")

    # ----- final summary -----
    print("\n--- final summary ---")
    print(f"  Reads processed : {reads_total}")
    print(f"  Nodes saved     : {len(merged)}")
    print(f"  Total time      : {time.perf_counter()-start:.2f}s\n")

# ─────────────────────────── CLI ─────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Extract node‑mapped read segments with strand and high un‑perfect ratio"
    )
    ap.add_argument("gam_file")
    ap.add_argument("stats_pickle")
    ap.add_argument("output_prefix")
    ap.add_argument("--save-format", choices=["json", "pkl", "both"],
                    default="json", help="File format(s) to write")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max_pending", type=int, default=16)
    ap.add_argument("--milestone", type=int, default=100_000_000)
    ap.add_argument("--chr", default="", help="Filter by chromosome name")
    args = ap.parse_args()

    run_pipeline(
        args.gam_file, args.stats_pickle, args.output_prefix, args.save_format,
        args.threads, args.max_pending, args.milestone, args.chr
    )

if __name__ == "__main__":
    main()
