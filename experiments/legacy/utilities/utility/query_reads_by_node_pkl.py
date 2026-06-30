#!/usr/bin/env python3
import argparse
import pickle
import sys
from typing import Dict, Tuple, Optional, List


def load_pickle(pkl_path: str) -> Dict[int, Dict[str, int]]:
    with open(pkl_path, "rb") as f:
        obj = pickle.load(f)

    if not isinstance(obj, dict):
        raise TypeError(f"Pickle root object is {type(obj)}, expected dict")

    # Normalize keys to int (in case they are str)
    out: Dict[int, Dict[str, int]] = {}
    for k, v in obj.items():
        try:
            ik = int(k)
        except Exception:
            continue
        if isinstance(v, dict):
            out[ik] = v
    return out


def compute_stats(rec: Dict[str, int]) -> Tuple[int, int, int, float, float]:
    perfect = int(rec.get("perfect", 0))
    not_perfect = int(rec.get("not_perfect", 0))
    total = perfect + not_perfect
    perfect_ratio = (perfect / total) if total > 0 else 0.0
    notperf_ratio = (not_perfect / total) if total > 0 else 0.0
    return perfect, not_perfect, total, perfect_ratio, notperf_ratio


def print_node(node_id: int, rec: Optional[Dict[str, int]]) -> None:
    if rec is None:
        print(f"[NOT FOUND] node_id={node_id}")
        return

    perfect, not_perfect, total, perfect_ratio, notperf_ratio = compute_stats(rec)
    max_read_length = int(rec.get("max_read_length", 0))
    max_cigar_length = int(rec.get("max_cigar_length", 0))

    print(f"node_id: {node_id}")
    print(f"  perfect:        {perfect}")
    print(f"  not_perfect:    {not_perfect}")
    print(f"  total:          {total}")
    print(f"  perfect_ratio:  {perfect_ratio:.6f}")
    print(f"  notperf_ratio:  {notperf_ratio:.6f}")
    print(f"  max_read_length (sum to_length on node): {max_read_length}")
    print(f"  max_cigar_length (chars):                {max_cigar_length}")


def read_ids_file(path: str) -> List[int]:
    ids: List[int] = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # allow either "123" or "123 ..." per line
            tok = s.split()[0]
            ids.append(int(tok))
    return ids


def top_k(data: Dict[int, Dict[str, int]], by: str, k: int, min_total: int) -> None:
    rows = []
    for nid, rec in data.items():
        p, np, total, pr, npr = compute_stats(rec)
        if total < min_total:
            continue

        if by == "not_perfect":
            key = np
        elif by == "perfect":
            key = p
        elif by == "total":
            key = total
        elif by == "notperf_ratio":
            key = npr
        elif by == "perfect_ratio":
            key = pr
        else:
            raise ValueError(by)

        rows.append((key, nid, p, np, total, pr, npr, rec))

    rows.sort(key=lambda x: x[0], reverse=True)
    rows = rows[:k]

    print(
        "rank\tkey\tnode_id\tperfect\tnot_perfect\ttotal\tperfect_ratio\tnotperf_ratio\t"
        "max_read_length\tmax_cigar_length"
    )
    for i, (key, nid, p, np, total, pr, npr, rec) in enumerate(rows, 1):
        mrl = int(rec.get("max_read_length", 0))
        mcl = int(rec.get("max_cigar_length", 0))
        key_str = f"{key:.6f}" if isinstance(key, float) else str(key)
        print(f"{i}\t{key_str}\t{nid}\t{p}\t{np}\t{total}\t{pr:.6f}\t{npr:.6f}\t{mrl}\t{mcl}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Query reads_by_node.pkl (node_id -> {perfect, not_perfect, max_read_length, max_cigar_length})."
    )
    ap.add_argument("pkl", help="Path to .pkl produced by your GAM parser")

    # Multiple IDs:
    ap.add_argument(
        "--nodes",
        type=int,
        nargs="+",
        default=[],
        help="One or more node IDs to query (space-separated). Example: --nodes 18 21 315301",
    )
    ap.add_argument(
        "--nodes_file",
        default="",
        help="Text file containing node IDs (one per line; lines starting with # ignored).",
    )

    # Keep your top-K feature
    ap.add_argument("--top", type=int, default=0, help="List top K nodes (0 disables).")
    ap.add_argument(
        "--by",
        choices=["not_perfect", "perfect", "total", "notperf_ratio", "perfect_ratio"],
        default="not_perfect",
        help="Metric to sort by for --top (default: not_perfect).",
    )
    ap.add_argument("--min_total", type=int, default=1, help="Min total segments for --top (default: 1).")

    args = ap.parse_args()

    data = load_pickle(args.pkl)
    print(f"[INFO] Loaded {len(data)} nodes from {args.pkl}")

    ids: List[int] = []
    if args.nodes_file:
        ids.extend(read_ids_file(args.nodes_file))
    ids.extend(args.nodes)

    if not ids and args.top <= 0:
        print("[ERROR] Provide --nodes/--nodes_file and/or --top.")
        sys.exit(2)

    # Dedup while preserving order
    seen = set()
    ids = [x for x in ids if not (x in seen or seen.add(x))]

    if ids:
        for nid in ids:
            print_node(nid, data.get(nid))
            print()  # blank line between records

    if args.top and args.top > 0:
        top_k(data, by=args.by, k=args.top, min_total=args.min_total)


if __name__ == "__main__":
    main()