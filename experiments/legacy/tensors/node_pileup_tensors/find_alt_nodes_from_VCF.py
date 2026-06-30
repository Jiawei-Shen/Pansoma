#!/usr/bin/env python3
"""
extract_snps_and_insertions.py
------------------------------
Emit SNPs (one row per single-base difference) and/or simple insertions
(ALT longer than REF and containing one new node) from a vg-deconstructed VCF.

Output columns (TAB):
  CHROM  POS  ID  TYPE  REF_BASE  REF_NODE  ALT_STR  ALT_NODE(S)
"""

import argparse, re, sys, pysam
from typing import List

_SPLIT = re.compile(r"[><]+")


def split_nodes(trav: str) -> List[int]:
    """'>15051>15050>15048' → [15051, 15050, 15048]"""
    return [int(tok) for tok in _SPLIT.split(trav) if tok]


# ------------------------------------------------------------------  SNP  ----
def snps_for_alt(ref_seq, ref_nodes, alt_seq, alt_nodes):
    """Yield (offset, ref_base, ref_node, alt_base, alt_node) for each SNP."""
    if len(ref_seq) != len(alt_seq):
        return
    for i, (rb, ab) in enumerate(zip(ref_seq, alt_seq)):
        if rb == ab:
            continue
        ref_node = next((n for n in ref_nodes[i:] if n not in alt_nodes),
                        ref_nodes[-1])
        alt_node = next((n for n in alt_nodes[i:] if n not in ref_nodes),
                        alt_nodes[-1])
        yield i, rb, ref_node, ab, alt_node


# -------------------------------------------------------------  INSERTION ----
def insertion_info(ref_seq, ref_nodes, alt_seq, alt_nodes):
    """
    Identifies an insertion event corresponding to one or more new nodes.

    This function defines an insertion as a variant where the alternate
    sequence is longer than the reference and its path contains at least
    one node not present in the reference path.

    Args:
        ref_seq (str): The reference sequence.
        ref_nodes (list): The list of node IDs for the reference path.
        alt_seq (str): The alternate sequence, which includes the insertion.
        alt_nodes (list): The list of node IDs for the alternate path.

    Returns:
        A tuple containing:
        (anchor_base, anchor_node, inserted_sequence, insertion_nodes_str, pos_offset)
        or None if the criteria for an insertion are not met.
    """
    # For an insertion, the alternate sequence must be longer than the reference.
    if len(alt_seq) <= len(ref_seq):
        return None

    # Find the set of nodes that are unique to the alternate path.
    insertion_nodes = set(alt_nodes) - set(ref_nodes)

    # **Key Change 1**: If there are no unique nodes, it's not a graph-based
    # insertion. This now allows for one or more nodes.
    if not insertion_nodes:
        return None

    # **Key Change 2**: Format the insertion nodes into a deterministic,
    # comma-separated string. Sorting ensures the output is consistent.
    insertion_nodes_str = ",".join(map(str, sorted(list(insertion_nodes))))

    # The inserted_sequence is the portion of the alt_seq that extends
    # beyond the ref_seq.
    inserted_sequence = alt_seq[len(ref_seq):]

    # The anchor position is defined by the end of the reference path.
    anchor_base = ref_seq[-1]
    anchor_node = ref_nodes[-1]
    position_offset = len(ref_seq) - 1

    # Return the details of the insertion.
    return anchor_base, anchor_node, inserted_sequence, insertion_nodes_str, position_offset


# ------------------------------------------------------------------- main ----
def main(args):
    """Main logic: open VCF, iterate records, and print selected variants."""
    vcf = pysam.VariantFile(args.vcf)
    if "AT" not in vcf.header.info:
        sys.exit("ERROR: VCF header lacks AT tag")

    print("CHROM\tPOS\tID\tTYPE\tREF_BASE\tREF_NODE\tALT_STR\tALT_NODE(S)")

    for rec in vcf:
        ref_seq   = rec.ref
        ref_nodes = split_nodes(rec.info["AT"][0])

        for alt_idx, alt_seq in enumerate(rec.alts):
            alt_nodes = split_nodes(rec.info["AT"][alt_idx + 1])

            # Conditionally process SNPs
            if args.type in ("snp", "both"):
                for off, rb, rn, ab, an in snps_for_alt(
                    ref_seq, ref_nodes, alt_seq, alt_nodes
                ):
                    print(
                        f"{rec.contig}\t{rec.pos + off}\t{rec.id}\tSNP\t"
                        f"{rb}\t{rn}\t{ab}\t{an}"
                    )

            # Conditionally process insertions
            if args.type in ("insertion", "both"):
                ins = insertion_info(ref_seq, ref_nodes, alt_seq, alt_nodes)
                if ins:
                    anchor_b, anchor_n, ins_seq, ins_nodes, pos_off = ins
                    pos = rec.pos + pos_off
                    print(
                        f"{rec.contig}\t{pos}\t{rec.id}\tINS\t"
                        f"{anchor_b}\t{anchor_n}\t{ins_seq}\t{ins_nodes}"
                    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Extract SNPs and/or insertions from a vg-deconstructed VCF."
    )
    ap.add_argument("vcf", help="Input .vcf or .vcf.gz produced by vg deconstruct")
    ap.add_argument(
        "-t", "--type",
        help="Type of variant to extract.",
        choices=["snp", "insertion", "both"],
        default="both"
    )
    args = ap.parse_args()
    main(args)