#!/usr/bin/env python3
import argparse
import json
import os
import sys
import time
from collections import defaultdict
import re  # <-- added

import numpy as np
import pysam

# ─────────────────────────────────────────────────────────────────────────────
# Helpers

def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def format_time(seconds: float) -> str:
    seconds = int(seconds)
    h, r = divmod(seconds, 3600)
    m, s = divmod(r, 60)
    return f"{h:02}:{m:02}:{s:02}"

def _parse_key(variant_key: str):
    """
    Expected canonical forms:
      SNP/other:   "{offset}_{type}_{ref}_{alt}"
      Deletion D:  "{anchor}_{D}_{anchor+deleted}_{anchor}"
      Insertion I: "{anchor}_{I}_{anchorBase}_{anchorBase+inserted}"

    Returns: (offset, vtype, REF, ALT) all uppercased.
    """
    parts = variant_key.split("_")
    if len(parts) < 4:
        raise ValueError(f"Unexpected variant_key: {variant_key}")
    offset = int(parts[0])
    vtype = parts[1]
    ref = parts[2].upper()
    alt = parts[3].upper()
    return offset, vtype, ref, alt

def _vcf_has_partial_match(vcf, chrom: str, pos: int, v_ref: str, v_alt: str, v_type: str) -> bool:
    """
    Fetch records at POS and test:
      - Deletion (D):     VCF.REF contains v_ref  AND  any VCF.ALT contains v_alt
      - Insertion (I):    any VCF.ALT contains v_alt  AND  VCF.REF contains v_ref
      - Other (e.g. X):   exact REF==v_ref and ALT==v_alt
    """
    for rec in vcf.fetch(chrom, max(0, pos - 1), pos + 1):
        ref_truth = (rec.ref or "").upper()
        alts_truth = [(a or "").upper() for a in (rec.alts or [])]
        if v_type == "D":
            if ref_truth and (v_ref[len(v_alt):] in ref_truth):
                return True
        elif v_type == "I":
            if ref_truth and (v_ref in ref_truth):
                for a in alts_truth:
                    if v_alt in a:
                        return True
        else:
            for a in alts_truth:
                if ref_truth == v_ref and a == v_alt:
                    return True
    return False

def load_node_positions(ref_json_path: str):
    """Return dict[int node_id] -> int grch38_position_start."""
    pos = {}
    with open(ref_json_path, "r") as f:
        data = json.load(f)
    nodes = data if isinstance(data, list) else data.get("nodes", [])
    for node in nodes:
        nid_raw = node.get("node_id")
        if nid_raw is None:
            continue
        try:
            nid = int(nid_raw)
        except (TypeError, ValueError):
            continue
        p = node.get("grch38_position_start")
        if isinstance(p, int):
            pos[nid] = p
    return pos

def detect_shard_bases(data_dir: str):
    """
    Scan data_dir for files like:
      PREFIX00000_data.npy
      shard_00000_data.npy
      COLO829T_ONT_chr3_00000_data.npy
    and build a mapping:
      shard_idx (int) -> base path WITHOUT '_data.npy'/'_labels.npy'

    Example:
      'COLO829T_ONT_chr3_00000_data.npy'
      -> idx = 0
      -> base = '/.../COLO829T_ONT_chr3_00000'
    """
    pattern = re.compile(r"^(.*?)(\d{5})_data\.npy$")
    mapping = {}
    for fname in os.listdir(data_dir):
        m = pattern.match(fname)
        if not m:
            continue
        prefix, idx_str = m.group(1), m.group(2)
        try:
            idx = int(idx_str)
        except ValueError:
            continue
        base = os.path.join(data_dir, prefix + idx_str)
        if idx in mapping and mapping[idx] != base:
            log(f"[WARN] Multiple data files found for shard {idx:05d}: "
                f"{os.path.basename(mapping[idx])} vs {fname}; using {fname}")
        mapping[idx] = base
    log(f"Detected data prefixes for {len(mapping)} shard(s) in {data_dir}")
    return mapping

# ─────────────────────────────────────────────────────────────────────────────
# Main classification stage

def main():
    ap = argparse.ArgumentParser(
        description="Classify variant tensors from variant_summary.ndjson using VCF + ref JSON."
    )
    ap.add_argument("variant_summary_ndjson",
                    help="variant_summary.ndjson produced by tensor-generation stage")
    ap.add_argument("ref_json",
                    help="Reference node JSON with grch38_position_start per node")
    ap.add_argument("vcf_file",
                    help="Truth VCF (bgzipped) with candidate variants")
    ap.add_argument("--chr", required=True,
                    help="Chromosome name for VCF fetch (e.g., chr1)")
    ap.add_argument("--data-dir", required=True,
                    help="Directory containing *data.npy; labels will be written here")
    ap.add_argument("--unknown-label", type=int, default=-1,
                    help="Integer label to assign when no decision can be made (default: -1)")
    args = ap.parse_args()

    vs_path = args.variant_summary_ndjson
    if not os.path.isfile(vs_path):
        sys.exit(f"variant_summary_ndjson not found: {vs_path}")
    if not os.path.isfile(args.ref_json):
        sys.exit(f"ref_json not found: {args.ref_json}")
    if not os.path.isfile(args.vcf_file):
        sys.exit(f"VCF not found: {args.vcf_file}")
    if not os.path.isdir(args.data_dir):
        sys.exit(f"data-dir not found or not a directory: {args.data_dir}")

    log("Loading node positions from ref JSON...")
    node_pos = load_node_positions(args.ref_json)
    log(f"Loaded positions for {len(node_pos):,} nodes from ref JSON.")

    log("Opening VCF with pysam...")
    vcf = pysam.VariantFile(args.vcf_file)

    # labels_by_shard[shard_idx] -> list of labels
    labels_by_shard = defaultdict(list)

    out_classified_path = os.path.join(
        os.path.dirname(vs_path),
        "variant_summary_classified.ndjson",
    )
    out_f = open(out_classified_path, "w")

    total_lines = 0
    total_true = total_false = total_unknown = 0
    start = time.time()

    log(f"Starting classification from {vs_path}...")
    with open(vs_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            total_lines += 1
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue

            node_id = rec.get("node_id")
            key = rec.get("variant_key")
            shard_idx = rec.get("shard_index")
            idx_within = rec.get("index_within_shard")

            if node_id is None or key is None or shard_idx is None or idx_within is None:
                continue

            try:
                node_id_int = int(node_id)
            except (TypeError, ValueError):
                continue

            # parse key -> offset, type, ref, alt
            try:
                offset, vtype, v_ref, v_alt = _parse_key(key)
            except ValueError:
                # cannot interpret key; mark unknown
                genomic_pos = None
                label_int = args.unknown_label
                classification = "unknown"
            else:
                start_pos = node_pos.get(node_id_int)
                if start_pos is None:
                    genomic_pos = None
                    label_int = args.unknown_label
                    classification = "unknown"
                else:
                    genomic_pos = start_pos + offset  # 1-based coordinate

                    is_match = _vcf_has_partial_match(
                        vcf=vcf,
                        chrom=args.chr,
                        pos=genomic_pos,
                        v_ref=v_ref,
                        v_alt=v_alt,
                        v_type=vtype,
                    )
                    if is_match:
                        label_int = 1
                        classification = "true"
                        total_true += 1
                    else:
                        # If we want a strict binary classification:
                        label_int = 0
                        classification = "false"
                        total_false += 1

            # Place label into labels_by_shard
            try:
                si = int(shard_idx)
                ii = int(idx_within)
            except (TypeError, ValueError):
                continue

            arr = labels_by_shard[si]
            # Grow list if needed, then assign at index
            if len(arr) <= ii:
                arr.extend([args.unknown_label] * (ii + 1 - len(arr)))
            arr[ii] = label_int
            if label_int == args.unknown_label:
                total_unknown += 1

            # Write updated record
            out_rec = dict(rec)
            out_rec["genomic_position"] = genomic_pos
            out_rec["classification"] = classification
            out_rec["label_int"] = label_int
            out_f.write(json.dumps(out_rec) + "\n")

            # Periodic progress
            if total_lines % 100000 == 0:
                dt = time.time() - start
                rate = total_lines / dt if dt > 0 else 0.0
                log(
                    f"Classified {total_lines:,} variants "
                    f"(true={total_true:,}, false={total_false:,}, unknown={total_unknown:,}) "
                    f"→ {rate:.1f} variants/s"
                )

    out_f.close()
    log(f"Classification NDJSON written to: {out_classified_path}")

    # Detect data prefixes from existing *_data.npy files
    shard_bases = detect_shard_bases(args.data_dir)

    # Now write *_labels.npy for each shard, matching the detected prefix
    log("Writing shard label arrays...")
    for shard_idx, labels_list in sorted(labels_by_shard.items()):
        # Ensure no None; replace with unknown_label if any
        labels_array = np.array(
            [lbl if lbl is not None else args.unknown_label for lbl in labels_list],
            dtype=np.int8,
        )

        base = shard_bases.get(shard_idx)
        if base is None:
            # Fallback to old naming scheme
            base = os.path.join(args.data_dir, f"shard_{shard_idx:05d}")
            log(
                f"[WARN] No data file detected for shard {shard_idx:05d}; "
                f"using default base name {os.path.basename(base)}"
            )

        labels_path = base + "_labels.npy"
        log(f"  shard {shard_idx:05d}: {labels_array.shape[0]} labels -> {labels_path}")
        np.save(labels_path, labels_array)

    elapsed = time.time() - start
    log("All done.")
    log(f"Total variants classified: {total_lines:,}")
    log(f"True:    {total_true:,}")
    log(f"False:   {total_false:,}")
    log(f"Unknown: {total_unknown:,}")
    log(f"Elapsed time: {format_time(elapsed)}")

if __name__ == "__main__":
    main()
