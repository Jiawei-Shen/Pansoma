#!/usr/bin/env python3

import argparse
import struct
import sys
import pysam  # We need pysam to read VCF files


def load_index(idx_path):
    """Loads node IDs from a binary .idx file. (Unchanged)"""
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}.", file=sys.stderr)
                return {}
            blocks_num, = struct.unpack("<I", blocks_num_bytes)
            for i in range(blocks_num):
                # Read only the relevant parts of the header to get the block_id
                header_data = f.read(22)  # Read the full header
                if len(header_data) < 22:
                    break
                block_id, _, _, _, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    f.read(metadata_len)  # Skip metadata
                node_index[block_id] = True  # We only care about the existence of the node ID
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}
    return node_index


def verify_variants_in_vcf(tsv_path, idx_nodes, vcf_path):
    """
    For variants in the TSV file whose nodes are in the index, verifies
    they exist at the correct locus in a standard VCF file.
    """
    print("Verifying TSV records against index and VCF...", file=sys.stderr)
    results = {
        "variants_to_check": 0,
        "matches_in_vcf": 0,
        "mismatches_in_vcf": 0
    }
    line_num = 0

    try:
        # Open the VCF file. It must be indexed (e.g., .vcf.gz.tbi) for fetch to work.
        vcf_file = pysam.VariantFile(vcf_path)
    except Exception as e:
        print(f"Error opening VCF file {vcf_path}. Is it a valid, indexed VCF? Details: {e}", file=sys.stderr)
        return None

    try:
        with open(tsv_path, 'r') as f:
            next(f, None)  # Skip header line
            for line in f:
                line_num += 1
                if not line.strip():
                    continue

                parts = line.strip().split()
                if len(parts) < 8:
                    continue

                # Extract variant data from the TSV line
                tsv_chrom, tsv_pos, _, _, tsv_ref, _, tsv_alt, alt_nodes_str = parts[:8]
                tsv_pos = int(tsv_pos)

                # First, check if any node on this line is a "matching node" from the index
                try:
                    node_ids = [int(n) for n in alt_nodes_str.split(',')]
                    is_node_in_idx = any(nid in idx_nodes for nid in node_ids)
                except ValueError:
                    is_node_in_idx = False

                # If no nodes from this line are in the index, we skip the VCF check for it.
                if not is_node_in_idx:
                    continue

                results["variants_to_check"] += 1
                is_variant_in_vcf = False

                # Now, query the VCF by coordinate to find the variant
                try:
                    # pysam is 0-based, VCF is 1-based. Fetching [pos-1, pos] gets records at that position.
                    for rec in vcf_file.fetch(tsv_chrom, tsv_pos - 1, tsv_pos + len(tsv_alt) + 1):
                        # We need an exact match on position, reference allele, and one of the alternate alleles
                        # if rec.pos == tsv_pos and rec.ref == tsv_ref and tsv_alt in rec.alts:
                        if rec:
                            is_variant_in_vcf = True
                            break  # Found it
                except ValueError:
                    # This error occurs if the chromosome from the TSV (e.g., "chr1") is not in the VCF header
                    print(f"Warning: Chromosome '{tsv_chrom}' from TSV line {line_num} not found in VCF.",
                          file=sys.stderr)

                if is_variant_in_vcf:
                    results["matches_in_vcf"] += 1
                else:
                    results["mismatches_in_vcf"] += 1
                    # print(
                    #     f"Mismatch: Variant {tsv_chrom}:{tsv_pos} {tsv_ref}>{tsv_alt} (from TSV line with valid node) not found in VCF.",
                    #     file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: TSV file not found at {tsv_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading TSV file {tsv_path} on line {line_num}: {e}",
              file=sys.stderr)
        return None

    return results


def main():
    """
    Main function to orchestrate the node verification process.
    """
    parser = argparse.ArgumentParser(
        description="For TSV variants whose nodes are in an index, verify they exist in a standard VCF file."
    )
    parser.add_argument("tsv_file", help="Path to the input TSV file.")
    parser.add_argument("idx_file", help="Path to the .idx index file.")
    parser.add_argument("vcf_file", help="Path to the standard, indexed VCF file (.vcf.gz).")
    args = parser.parse_args()

    # 1. Load node IDs from the index file
    idx_data = load_index(args.idx_file)
    if not idx_data:
        sys.exit(1)
    idx_nodes = set(idx_data.keys())
    print(f"Loaded {len(idx_nodes)} unique nodes from the index.", file=sys.stderr)

    # 2. Verify the relevant TSV records against the VCF
    results = verify_variants_in_vcf(args.tsv_file, idx_nodes, args.vcf_file)

    # 3. Output the final summary
    if results:
        print("\n--- Verification Summary ---")
        print(f"Total variants from TSV with nodes in index: {results['variants_to_check']}")
        print(f"  - Variants found in VCF at correct locus:    {results['matches_in_vcf']}")
        print(f"  - Variants NOT found in VCF:                 {results['mismatches_in_vcf']}")


if __name__ == "__main__":
    main()