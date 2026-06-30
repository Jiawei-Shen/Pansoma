#!/usr/bin/env python3
"""
Merge node JSONs:
- Input files look like:
  {
    "path_name_input_pattern": "^W\\tGRCh38\\t0\\tchr1",
    "path_identifier_gfa": "W:GRCh38/0/chr1",
    "path_source_type": "W",
    "nodes": [ ... node objects ... ]
  }

- Output: ONE JSON array (list) of node objects. For each node:
  * Add a chromosome field (default key: "chrom"), e.g. "chr1"
  * Convert AF list -> digit string with 8 bins:
      0: [0, 1e-6)
      1: [1e-6, 1e-5)
      2: [1e-5, 1e-4)
      3: [1e-4, 1e-3)
      4: [1e-3, 1e-2)
      5: [1e-2, 0.1)
      6: [0.1, 0.5)
      7: [0.5, 1.0]   (1.0 inclusive)

Usage:
  python merge_nodes.py -o merged.json file1.json [file2.json ...]
Options:
  --af-field GENOMEAD field name to transform (default: genomead_af)
  --chrom-key KEY       field name to store chromosome (default: chrom)
"""

import argparse
import gzip
import json
import os
import re
import sys
from glob import glob
from typing import Any, Dict, Iterable, List, Optional

CHR_RE = re.compile(r"(chr[0-9XYM]+)", re.IGNORECASE)

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def extract_chr(meta: Dict[str, Any], fallback_name: str) -> Optional[str]:
    """Try to find chrN from header fields or file name."""
    for key in ("path_identifier_gfa", "path_name_input_pattern"):
        val = meta.get(key)
        if isinstance(val, str):
            m = CHR_RE.search(val)
            if m:
                return m.group(1)
    # fallback to file name
    m = CHR_RE.search(os.path.basename(fallback_name))
    return m.group(1) if m else None

def af_to_bin_digit(af: float) -> str:
    """Map AF to bin digit '0'..'7' per specified ranges."""
    # Handle non-finite
    try:
        x = float(af)
    except Exception:
        return "0"
    if x < 0:
        x = 0.0
    # edges define [edge[i], edge[i+1]) except last is inclusive
    edges = [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5]
    for i in range(7):
        if edges[i] <= x < edges[i+1]:
            return str(i)
    # last bin: [0.5, 1.0] inclusive (and anything >=1 clamps to 7)
    return "7"

def af_list_to_digit_string(vals: Iterable[Any]) -> str:
    return "".join(af_to_bin_digit(v) for v in vals)

def load_nodes_from_file(path: str) -> (List[Dict[str, Any]], Optional[str]):
    """Read one JSON; return (nodes, chr). Accepts dict-with-nodes or bare list."""
    with open_maybe_gzip(path, "rt") as f:
        data = json.load(f)

    if isinstance(data, dict) and "nodes" in data:
        chrom = extract_chr(data, path)
        nodes = data["nodes"] if isinstance(data["nodes"], list) else []
        return nodes, chrom
    elif isinstance(data, list):
        # Already a list of nodes (no header)
        chrom = extract_chr({}, path)
        return data, chrom
    else:
        return [], extract_chr(data if isinstance(data, dict) else {}, path)

def main():
    ap = argparse.ArgumentParser(description="Merge node JSONs; add chromosome; convert AF list to digit string.")
    ap.add_argument("inputs", nargs="+", help="Input JSON/JSON.GZ files (shell globs OK)")
    ap.add_argument("-o", "--out", required=True, help="Output JSON path")
    ap.add_argument("--af-field", default="genomead_af", help="Node field containing AF list to convert")
    ap.add_argument("--chrom-key", default="chrom", help="Key name to store chromosome inside each node")
    args = ap.parse_args()

    # Expand any globs passed as literals (in case shell didn't expand)
    file_list: List[str] = []
    for pat in args.inputs:
        matched = glob(pat)
        file_list.extend(matched if matched else [pat])

    out_nodes: List[Dict[str, Any]] = []

    for fp in file_list:
        nodes, chrom = load_nodes_from_file(fp)

        # If we couldn't parse chr from header/file name, leave None (or set explicitly)
        for n in nodes:
            if chrom is not None:
                n[args.chrom_key] = chrom

            # Convert AF list -> number string
            if args.af_field in n:
                af_val = n[args.af_field]
                if isinstance(af_val, list):
                    n[args.af_field] = af_list_to_digit_string(af_val)
                else:
                    # handle scalar or other types
                    n[args.af_field] = af_list_to_digit_string([af_val])

        out_nodes.extend(nodes)

    # Write a JSON *list* of nodes
    with open(args.out, "w") as fo:
        json.dump(out_nodes, fo, ensure_ascii=False)

if __name__ == "__main__":
    main()
