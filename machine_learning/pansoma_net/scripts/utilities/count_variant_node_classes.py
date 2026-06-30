#!/usr/bin/env python3
"""
Count variants by node class (reference / ref-alt / alt-alt) from a node/offset VCF.

Classification rule (per user spec):
  - reference node: JSON record HAS `genomead_af` (case-insensitive; nested allowed)
  - ref-alt node : JSON record exists but DOES NOT have `genomead_af`
  - alt-alt node : node_id not present in JSON at all

Inputs:
  --in_vcf    : VCF with CHROM = numeric node_id (non-numeric CHROM lines are ignored but tallied)
  --map_json  : Node JSON (array-of-objects OR dict-of-objects), or JSONL with --map_jsonl
  --map_jsonl : Interpret --map_json as JSON Lines (one dict per line)
  --out_tsv   : (optional) write counts as a TSV

Usage example:
  python count_variant_node_classes.py \
      --in_vcf nodes.vcf.gz \
      --map_json nodes.json \
      --out_tsv node_counts.tsv
"""

import argparse
import gzip
import json
import sys
from typing import Dict, Iterable, Iterator, Optional, Tuple


# ---------------- JSON helpers ----------------

ID_KEYS = {"node_id", "node", "id", "nid"}

def _get_first(d: dict, keys: Iterable[str]):
    for k in keys:
        if k in d:
            return d[k]
    return None

def _node_id_from_obj(obj: dict) -> Optional[int]:
    nid = _get_first(obj, ID_KEYS)
    if nid is None:
        return None
    try:
        return int(nid)
    except Exception:
        return None

def _has_key_recursive(obj: dict, wanted_keys_lower: set) -> bool:
    """Return True if any of the wanted keys (lowercased) appear anywhere in (possibly nested) dict."""
    for k, v in obj.items():
        if str(k).lower() in wanted_keys_lower:
            return True
        if isinstance(v, dict) and _has_key_recursive(v, wanted_keys_lower):
            return True
        # If lists may contain dicts, scan them too
        if isinstance(v, (list, tuple)):
            for elt in v:
                if isinstance(elt, dict) and _has_key_recursive(elt, wanted_keys_lower):
                    return True
    return False

def iter_json_records(path: str, jsonl: bool) -> Iterator[dict]:
    """Yield dict records from JSON (list/dict) or JSONL (one dict per line)."""
    if jsonl:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or not s.startswith("{"):
                    continue
                try:
                    obj = json.loads(s)
                except Exception:
                    continue
                if isinstance(obj, dict):
                    yield obj
    else:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            for v in data.values():
                if isinstance(v, dict):
                    yield v
        elif isinstance(data, list):
            for v in data:
                if isinstance(v, dict):
                    yield v
        else:
            raise ValueError("Unsupported JSON structure: expected list/dict of dicts, or use --map_jsonl")

def build_node_class_map(json_path: str, jsonl: bool) -> Dict[int, str]:
    """
    Returns node_id -> {'reference','ref-alt'}
    Nodes not present in the map are 'alt-alt'.
    """
    cls_map: Dict[int, str] = {}
    wanted = {"genomead_af"}  # case-insensitive match

    for obj in iter_json_records(json_path, jsonl):
        nid = _node_id_from_obj(obj)
        if nid is None:
            continue
        if _has_key_recursive(obj, wanted):
            cls_map[nid] = "reference"
        else:
            cls_map[nid] = "ref-alt"
    return cls_map


# ---------------- VCF helpers ----------------

def open_text_auto(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")

def count_variants_by_class(vcf_path: str, node_class: Dict[int, str]) -> Tuple[int, int, int, int]:
    """
    Returns: (n_reference, n_ref_alt, n_alt_alt, n_non_numeric_chrom)
    Only lines with numeric CHROM (node ids) are classified; others are tallied separately.
    """
    n_ref = n_refalt = n_altalt = n_nonnum = 0
    with open_text_auto(vcf_path) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            chrom = cols[0]
            try:
                nid = int(chrom)
            except Exception:
                n_nonnum += 1
                continue

            cls = node_class.get(nid)
            if cls == "reference":
                n_ref += 1
            elif cls == "ref-alt":
                n_refalt += 1
            else:
                n_altalt += 1
    return n_ref, n_refalt, n_altalt, n_nonnum

def pct(n: int, d: int) -> str:
    return f"{(100.0 * n / d):.2f}%" if d else "0.00%"


# ---------------- CLI ----------------

def main():
    ap = argparse.ArgumentParser(description="Count variants on reference / ref-alt / alt-alt nodes (uses 'genomead_af' to identify reference nodes).")
    ap.add_argument("--in_vcf", required=True, help="Input VCF (.vcf or .vcf.gz) with CHROM=node_id.")
    ap.add_argument("--map_json", required=True, help="Node metadata JSON (array/dict) or JSONL (use --map_jsonl).")
    ap.add_argument("--map_jsonl", action="store_true", help="Interpret --map_json as JSONL (one JSON object per line).")
    ap.add_argument("--out_tsv", default=None, help="Optional: write counts to this TSV.")
    args = ap.parse_args()

    print("=== Node class counting ===")
    print(f"VCF       : {args.in_vcf}")
    print(f"Node JSON : {args.map_json} (jsonl={bool(args.map_jsonl)})")

    node_class = build_node_class_map(args.map_json, args.map_jsonl)
    n_ref_nodes  = sum(1 for c in node_class.values() if c == "reference")
    n_ralt_nodes = sum(1 for c in node_class.values() if c == "ref-alt")
    print(f"Node classes from JSON -> reference:{n_ref_nodes}  ref-alt:{n_ralt_nodes}  (others => alt-alt)")

    n_ref, n_refalt, n_altalt, n_nonnum = count_variants_by_class(args.in_vcf, node_class)
    total = n_ref + n_refalt + n_altalt
    grand = total + n_nonnum

    print("\n=== Variant counts (node-CHROM only) ===")
    print(f"reference : {n_ref:10d}  ({pct(n_ref, total)})")
    print(f"ref-alt   : {n_refalt:10d}  ({pct(n_refalt, total)})")
    print(f"alt-alt   : {n_altalt:10d}  ({pct(n_altalt, total)})")
    print(f"----------------------------------------")
    print(f"TOTAL     : {total:10d}")
    if n_nonnum:
        print(f"Non-numeric CHROM (ignored in classes): {n_nonnum}  (grand total lines = {grand})")

    if args.out_tsv:
        try:
            with open(args.out_tsv, "w", encoding="utf-8") as w:
                w.write("class\tcount\tpercent_of_node_records\n")
                w.write(f"reference\t{n_ref}\t{pct(n_ref, total)}\n")
                w.write(f"ref-alt\t{n_refalt}\t{pct(n_refalt, total)}\n")
                w.write(f"alt-alt\t{n_altalt}\t{pct(n_altalt, total)}\n")
                w.write(f"TOTAL_NODE_CHROM\t{total}\t100.00%\n")
                if n_nonnum:
                    w.write(f"NON_NUMERIC_CHROM_IGNORED\t{n_nonnum}\t-\n")
            print(f"\nWrote TSV: {args.out_tsv}")
        except Exception as e:
            print(f"WARNING: failed to write --out_tsv: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
