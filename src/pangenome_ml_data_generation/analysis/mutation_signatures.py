#!/usr/bin/env python3

import argparse
import gzip
import os
from collections import Counter


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1].upper()


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def parse_gfa_nodes(gfa_path):
    nodes = {}

    with open_maybe_gz(gfa_path) as f:
        for line in f:
            if not line.startswith("S\t"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            node_id = fields[1]
            seq = fields[2].upper()

            if seq == "*":
                continue

            nodes[node_id] = seq

    return nodes


def load_neighbor_map(path):
    """
    TSV format:

    node_id    left_node    right_node

    Missing value can be "." or empty.

    Example:
    1001       1000         1002
    1005       .            1006
    1009       1008         .
    """
    neighbor_map = {}

    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        header_lc = [h.lower() for h in header]

        required = ["node_id", "left_node", "right_node"]
        for col in required:
            if col not in header_lc:
                raise ValueError(f"--neighbor-map missing required column: {col}")

        i_node = header_lc.index("node_id")
        i_left = header_lc.index("left_node")
        i_right = header_lc.index("right_node")

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            node_id = fields[i_node].strip()
            left_raw = fields[i_left].strip() if i_left < len(fields) else "."
            right_raw = fields[i_right].strip() if i_right < len(fields) else "."

            left_node = left_raw if left_raw not in ("", ".") else None
            right_node = right_raw if right_raw not in ("", ".") else None

            neighbor_map[str(node_id)] = {
                "left": str(left_node) if left_node is not None else None,
                "right": str(right_node) if right_node is not None else None,
            }

    return neighbor_map


def all_sbs96_channels():
    bases = ["A", "C", "G", "T"]

    mutation_types = [
        ("C", "A"),
        ("C", "G"),
        ("C", "T"),
        ("T", "A"),
        ("T", "C"),
        ("T", "G"),
    ]

    channels = []

    for left in bases:
        for ref, alt in mutation_types:
            for right in bases:
                channels.append(f"{left}[{ref}>{alt}]{right}")

    return channels


def get_graph_position(fields):
    """
    Return:
        node_id, pos_1based

    For graph VCF:
        POS is node-local 1-based.
        INFO may contain NID and NSO.

    If NID/NSO exists:
        POS = NSO + 1
    """
    chrom = fields[0]
    pos = fields[1]
    info = parse_info(fields[7])

    if "NID" in info and "NSO" in info:
        node_id = str(info["NID"])
        pos_1based = int(info["NSO"]) + 1
        return node_id, pos_1based

    node_id = str(chrom)
    pos_1based = int(pos)
    return node_id, pos_1based


def get_context_same_node_only(nodes, node_id, pos_1based):
    """
    Strict default:
    only use trinucleotide context from the same node.

    If variant is at node boundary, return unresolved.
    """
    node_id = str(node_id)

    if node_id not in nodes:
        return None, "missing_node"

    seq = nodes[node_id]
    idx = pos_1based - 1

    if idx < 0 or idx >= len(seq):
        return None, "position_out_of_range"

    if idx == 0:
        return None, "left_boundary_context_missing"

    if idx == len(seq) - 1:
        return None, "right_boundary_context_missing"

    context = seq[idx - 1: idx + 2].upper()

    if len(context) != 3:
        return None, "bad_context_length"

    if any(base not in "ACGT" for base in context):
        return None, "non_acgt_context"

    return context, "ok"


def get_neighbor_node_from_map(neighbor_map, node_id, side):
    if neighbor_map is None:
        return None

    node_id = str(node_id)

    if node_id not in neighbor_map:
        return None

    return neighbor_map[node_id].get(side)


def get_neighbor_node_by_id_plus_minus_one(node_id, side):
    """
    Heuristic fallback:
    left  = node_id - 1
    right = node_id + 1

    This is not graph-topology-safe unless node IDs are path ordered.
    """
    try:
        nid = int(node_id)
    except ValueError:
        return None

    if side == "left":
        return str(nid - 1)

    if side == "right":
        return str(nid + 1)

    return None


def get_context_with_neighbor_rescue(
    nodes,
    node_id,
    pos_1based,
    neighbor_map=None,
    allow_node_id_plus_minus_one=False,
):
    """
    Optional rescue:
    If trinucleotide context crosses node boundary, use left/right neighbor node.

    Priority:
      1. --neighbor-map if provided
      2. node_id +/- 1 only if --allow-node-id-plus-minus-one is set

    Default behavior of the whole script remains no rescue.
    """
    node_id = str(node_id)

    if node_id not in nodes:
        return None, "missing_node"

    seq = nodes[node_id]
    idx = pos_1based - 1

    if idx < 0 or idx >= len(seq):
        return None, "position_out_of_range"

    center = seq[idx]

    # Left base
    if idx > 0:
        left = seq[idx - 1]
    else:
        left_node = get_neighbor_node_from_map(neighbor_map, node_id, "left")

        if left_node is None and allow_node_id_plus_minus_one:
            left_node = get_neighbor_node_by_id_plus_minus_one(node_id, "left")

        if left_node is None:
            return None, "left_boundary_neighbor_not_available"

        left_node = str(left_node)

        if left_node not in nodes:
            return None, "left_boundary_neighbor_missing_in_gfa"

        if len(nodes[left_node]) == 0:
            return None, "left_boundary_neighbor_empty"

        left = nodes[left_node][-1]

    # Right base
    if idx < len(seq) - 1:
        right = seq[idx + 1]
    else:
        right_node = get_neighbor_node_from_map(neighbor_map, node_id, "right")

        if right_node is None and allow_node_id_plus_minus_one:
            right_node = get_neighbor_node_by_id_plus_minus_one(node_id, "right")

        if right_node is None:
            return None, "right_boundary_neighbor_not_available"

        right_node = str(right_node)

        if right_node not in nodes:
            return None, "right_boundary_neighbor_missing_in_gfa"

        if len(nodes[right_node]) == 0:
            return None, "right_boundary_neighbor_empty"

        right = nodes[right_node][0]

    context = (left + center + right).upper()

    if len(context) != 3:
        return None, "bad_context_length"

    if any(base not in "ACGT" for base in context):
        return None, "non_acgt_context"

    return context, "ok_neighbor_rescued"


def canonical_sbs96(context, ref, alt):
    context = context.upper()
    ref = ref.upper()
    alt = alt.upper()

    if len(context) != 3:
        return None

    if len(ref) != 1 or len(alt) != 1:
        return None

    if ref not in "ACGT" or alt not in "ACGT":
        return None

    if ref == alt:
        return None

    if ref in {"C", "T"}:
        left, center, right = context
        return f"{left}[{ref}>{alt}]{right}"

    rc_context = revcomp(context)
    rc_ref = revcomp(ref)
    rc_alt = revcomp(alt)

    left, center, right = rc_context
    return f"{left}[{rc_ref}>{rc_alt}]{right}"


def sample_name_from_path(path):
    base = os.path.basename(path)

    for suffix in [".vcf.gz", ".vcf", ".gz"]:
        if base.endswith(suffix):
            base = base[: -len(suffix)]

    return base


def count_sbs96_for_vcf(
    vcf_path,
    nodes,
    min_prob=None,
    min_af=None,
    require_pass=True,
    allow_neighbor_rescue=False,
    neighbor_map=None,
    allow_node_id_plus_minus_one=False,
):
    counts = Counter()
    skipped = Counter()
    mismatch_examples = []

    with open_maybe_gz(vcf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                skipped["malformed_vcf_line"] += 1
                continue

            ref = fields[3].upper()
            alt = fields[4].upper()
            filt = fields[6]
            info = parse_info(fields[7])

            if require_pass and filt != "PASS":
                skipped["not_pass"] += 1
                continue

            vtype = info.get("TYPE")
            if vtype is not None and vtype != "X":
                skipped["not_type_x"] += 1
                continue

            if len(ref) != 1 or len(alt) != 1:
                skipped["not_snv"] += 1
                continue

            if "," in alt:
                skipped["multi_allelic"] += 1
                continue

            if min_prob is not None:
                if "PROB" not in info:
                    skipped["missing_prob"] += 1
                    continue

                try:
                    prob = float(info["PROB"])
                except ValueError:
                    skipped["bad_prob"] += 1
                    continue

                if prob <= min_prob:
                    skipped["low_or_equal_prob"] += 1
                    continue

            if min_af is not None:
                if "AF" not in info:
                    skipped["missing_af"] += 1
                    continue

                try:
                    af = float(info["AF"])
                except ValueError:
                    skipped["bad_af"] += 1
                    continue

                if af < min_af:
                    skipped["low_af"] += 1
                    continue

            try:
                node_id, pos_1based = get_graph_position(fields)
            except Exception:
                skipped["bad_graph_position"] += 1
                continue

            if allow_neighbor_rescue:
                context, status = get_context_with_neighbor_rescue(
                    nodes=nodes,
                    node_id=node_id,
                    pos_1based=pos_1based,
                    neighbor_map=neighbor_map,
                    allow_node_id_plus_minus_one=allow_node_id_plus_minus_one,
                )
            else:
                context, status = get_context_same_node_only(
                    nodes=nodes,
                    node_id=node_id,
                    pos_1based=pos_1based,
                )

            if context is None:
                skipped[status] += 1
                continue

            if status == "ok_neighbor_rescued":
                skipped["neighbor_rescued_context_used"] += 1

            if context[1] != ref:
                skipped["ref_context_mismatch"] += 1

                if len(mismatch_examples) < 10:
                    mismatch_examples.append(
                        {
                            "node_id": node_id,
                            "pos_1based": pos_1based,
                            "ref": ref,
                            "alt": alt,
                            "context": context,
                            "center": context[1],
                        }
                    )

                continue

            channel = canonical_sbs96(context, ref, alt)

            if channel is None:
                skipped["cannot_canonicalize"] += 1
                continue

            counts[channel] += 1

    return counts, skipped, mismatch_examples


def write_cosmic_matrix(sample_to_counts, out_path):
    channels = all_sbs96_channels()
    sample_names = list(sample_to_counts.keys())

    with open(out_path, "w") as out:
        out.write("Mutation Types\t" + "\t".join(sample_names) + "\n")

        for channel in channels:
            row = [channel]
            for sample in sample_names:
                row.append(str(sample_to_counts[sample].get(channel, 0)))
            out.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert graph VCF(s) + GFA into COSMIC-style SBS96 mutation matrix."
    )

    parser.add_argument(
        "--gfa",
        required=True,
        help="GFA file containing node sequences",
    )

    parser.add_argument(
        "--vcf",
        nargs="+",
        required=True,
        help="One or more graph VCF/VCF.GZ files",
    )

    parser.add_argument(
        "--sample-names",
        nargs="+",
        default=None,
        help="Optional sample names, same number/order as --vcf",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output COSMIC-style SBS96 matrix TSV",
    )

    parser.add_argument(
        "--min-prob",
        type=float,
        default=None,
        help="Only keep variants with PROB > this value",
    )

    parser.add_argument(
        "--min-af",
        type=float,
        default=None,
        help="Only keep variants with AF >= this value",
    )

    parser.add_argument(
        "--allow-non-pass",
        action="store_true",
        help="Do not require FILTER=PASS",
    )

    parser.add_argument(
        "--allow-neighbor-rescue",
        action="store_true",
        help=(
            "Enable left/right neighbor rescue for boundary trinucleotide context. "
            "Default: OFF."
        ),
    )

    parser.add_argument(
        "--neighbor-map",
        default=None,
        help=(
            "Optional TSV with columns: node_id, left_node, right_node. "
            "Used only when --allow-neighbor-rescue is enabled."
        ),
    )

    parser.add_argument(
        "--allow-node-id-plus-minus-one",
        action="store_true",
        help=(
            "If --allow-neighbor-rescue is enabled and --neighbor-map does not provide "
            "a neighbor, use node_id-1/node_id+1 as a fallback. Not recommended unless "
            "node IDs are path ordered."
        ),
    )

    args = parser.parse_args()

    if args.sample_names is not None and len(args.sample_names) != len(args.vcf):
        raise ValueError("--sample-names must have the same number of entries as --vcf")

    neighbor_map = None
    if args.neighbor_map is not None:
        if not args.allow_neighbor_rescue:
            print("[WARNING] --neighbor-map provided but --allow-neighbor-rescue is OFF; neighbor map will not be used.")
        else:
            print(f"[INFO] Loading neighbor map: {args.neighbor_map}")
            neighbor_map = load_neighbor_map(args.neighbor_map)
            print(f"[INFO] Loaded neighbor map records: {len(neighbor_map):,}")

    print(f"[INFO] Loading GFA nodes: {args.gfa}")
    nodes = parse_gfa_nodes(args.gfa)
    print(f"[INFO] Loaded nodes: {len(nodes):,}")

    print(f"[INFO] Neighbor rescue enabled: {args.allow_neighbor_rescue}")
    print(f"[INFO] node_id +/- 1 fallback enabled: {args.allow_node_id_plus_minus_one}")

    sample_to_counts = {}

    for i, vcf_path in enumerate(args.vcf):
        if args.sample_names is not None:
            sample = args.sample_names[i]
        else:
            sample = sample_name_from_path(vcf_path)

        print("============================================================")
        print(f"[INFO] Sample: {sample}")
        print(f"[INFO] VCF: {vcf_path}")

        if args.min_prob is not None:
            print(f"[INFO] Applying PROB filter: PROB > {args.min_prob}")

        if args.min_af is not None:
            print(f"[INFO] Applying AF filter: AF >= {args.min_af}")

        counts, skipped, mismatch_examples = count_sbs96_for_vcf(
            vcf_path=vcf_path,
            nodes=nodes,
            min_prob=args.min_prob,
            min_af=args.min_af,
            require_pass=not args.allow_non_pass,
            allow_neighbor_rescue=args.allow_neighbor_rescue,
            neighbor_map=neighbor_map,
            allow_node_id_plus_minus_one=args.allow_node_id_plus_minus_one,
        )

        sample_to_counts[sample] = counts

        print(f"[INFO] Counted SBS96 SNVs: {sum(counts.values()):,}")

        if skipped:
            print("[INFO] Skipped/rescue records:")
            for k, v in skipped.most_common():
                print(f"  {k}: {v:,}")

        if mismatch_examples:
            print("[WARNING] First REF/context mismatch examples:")
            for ex in mismatch_examples:
                print(
                    f"  node={ex['node_id']} pos={ex['pos_1based']} "
                    f"REF={ex['ref']} ALT={ex['alt']} "
                    f"context={ex['context']} center={ex['center']}"
                )

    write_cosmic_matrix(sample_to_counts, args.out)

    print("============================================================")
    print(f"[INFO] Wrote COSMIC-style SBS96 matrix: {args.out}")
    print(f"[INFO] Samples: {len(sample_to_counts)}")
    print("[INFO] Rows: 96")


if __name__ == "__main__":
    main()