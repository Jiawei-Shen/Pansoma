#!/usr/bin/env python3
"""
Convert a node/offset VCF to linear-reference coordinates and write TWO outputs:

  1) --out_linear_vcf.gz : records successfully converted to linear (CHROM=chr*, POS=linear)
  2) --out_nodes_vcf.gz  : records that stayed in node form (CHROM=node_id, POS=offset+1 for VCF)

Additionally prints:
  - Variant counts by class (reference / ref-alt / alt-alt)
  - DISTINCT NODE COUNTS by class (reference / ref-alt / alt-alt) observed in the VCF
  - TSV coverage of ref-alt nodes **against the TSV ALT column**:
      how many ref-alt nodes seen in the VCF are/aren't present in the TSV ALT column

Classification rule for nodes (JSON-driven):
  - reference: JSON record contains 'genomead_af' (case-insensitive, nested allowed)
  - ref-alt  : JSON record exists but lacks 'genomead_af'
  - alt-alt  : node id not present in the JSON

VCF COMPLIANCE:
- Only standard 8 columns are written.
- Original node/offset are preserved in INFO as: NID=<node_id>;NSO=<offset>
- For node-form output, POS is written as (offset + 1) to satisfy VCF/tabix (1-based).
"""

import argparse
import csv
import gzip
import json
import os
import re
import sys
from typing import Dict, Iterable, List, Optional, Tuple, Union, Set

from tqdm import tqdm

# pysam for bgzip + tabix
try:
    import pysam
except Exception:
    pysam = None


# ------------------------ TSV (ALT->REF) mapping ------------------------ #

def load_alt_to_ref_tsv(tsv_path: Optional[str]) -> Tuple[Dict[int, int], Set[int], Set[int]]:
    """
    Returns:
      alt_to_ref     : dict mapping ALT node -> REF node
      tsv_ref_nodes  : set of REF node ids present in TSV (even if no ALT listed)
      tsv_alt_nodes  : set of ALT node ids present in TSV
    """
    if not tsv_path:
        return {}, set(), set()

    alt_to_ref: Dict[int, int] = {}
    tsv_ref_nodes: Set[int] = set()
    tsv_alt_nodes: Set[int] = set()

    with open(tsv_path, "r", encoding="utf-8") as f:
        try:
            csv.field_size_limit(sys.maxsize)
        except OverflowError:
            csv.field_size_limit(2**31 - 1)
        reader = csv.DictReader(f, delimiter="\t")
        hdr = [c.strip() for c in (reader.fieldnames or [])]
        ref_node_key = None
        alt_nodes_key = None
        for k in hdr:
            lk = k.lower()
            if ref_node_key is None and lk in ("ref_node", "refnode", "ref"):
                ref_node_key = k
            if alt_nodes_key is None and lk in ("alt_node(s)", "alt_nodes", "altnodes", "alt"):
                alt_nodes_key = k
        if ref_node_key is None or alt_nodes_key is None:
            raise ValueError("TSV must have REF_NODE and ALT_NODE(S) columns (case-insensitive).")

        for row in reader:
            # collect REF node
            try:
                ref_node = int(str(row[ref_node_key]).strip())
                tsv_ref_nodes.add(ref_node)
            except Exception:
                continue

            # collect ALT nodes and ALT->REF mapping
            alts_raw = str(row[alt_nodes_key]).strip()
            for token in re.findall(r"\d+", alts_raw):
                try:
                    alt_node = int(token)
                except Exception:
                    continue
                alt_to_ref[alt_node] = ref_node
                tsv_alt_nodes.add(alt_node)

    return alt_to_ref, tsv_ref_nodes, tsv_alt_nodes


# ------------------------ Node map (JSON/JSONL) ------------------------ #

def _first_present(d: dict, keys: Iterable[str], default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default

def _to_int(x) -> Optional[int]:
    try:
        return int(x)
    except Exception:
        return None

def iter_node_map_records(json_path: str, jsonl: bool) -> Iterable[dict]:
    if jsonl:
        with open(json_path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or not s.startswith("{"):
                    continue
                try:
                    obj = json.loads(s)
                    if isinstance(obj, dict):
                        yield obj
                except Exception:
                    continue
    else:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, list):
            for obj in data:
                if isinstance(obj, dict):
                    yield obj
        elif isinstance(data, dict):
            for _, v in data.items():
                if isinstance(v, dict):
                    yield v
        else:
            raise ValueError("Unsupported JSON structure for node map.")

def _has_key_case_insensitive(obj: dict, key_name: str) -> bool:
    """True if key_name (case-insensitive) appears anywhere in obj (recursively)."""
    kn = key_name.lower()
    for k, v in obj.items():
        if str(k).lower() == kn:
            return True
        if isinstance(v, dict) and _has_key_case_insensitive(v, key_name):
            return True
        if isinstance(v, (list, tuple)):
            for elt in v:
                if isinstance(elt, dict) and _has_key_case_insensitive(elt, key_name):
                    return True
    return False

def load_node_map(json_path: str,
                  jsonl: bool,
                  node_start_is_1_based: bool = True) -> Tuple[
                      Dict[int, Tuple[str, int, str, int]],  # anchors: node -> (chrom, start0, strand, length)
                      Set[int],                              # nodes_seen in JSON
                      Set[int],                              # nodes_with_genomead_af (reference nodes by rule)
                  ]:
    anchors: Dict[int, Tuple[str, int, str, int]] = {}
    nodes_seen: Set[int] = set()
    nodes_with_genomead_af: Set[int] = set()

    for rec in iter_node_map_records(json_path, jsonl):
        nid = _first_present(rec, ["node_id", "node", "id", "nid"])
        node_id = _to_int(nid)
        if node_id is None:
            continue
        nodes_seen.add(node_id)

        if _has_key_case_insensitive(rec, "genomead_af"):
            nodes_with_genomead_af.add(node_id)

        chrom = _first_present(rec, ["chrom", "chr", "chromosome", "contig"])
        start1 = _first_present(rec, ["grch38_position_start", "start_1based", "pos1"])
        start0 = _first_present(rec, ["start_0based", "pos0"])
        strand = str(_first_present(rec, ["strand_in_path", "strand"], "+") or "+")
        length = _to_int(_first_present(rec, ["length", "len", "node_length"]))

        # decide start0
        if start0 is not None:
            s0 = _to_int(start0)
        elif start1 is not None:
            s1 = _to_int(start1)
            s0 = (s1 - 1) if s1 is not None else None
        else:
            generic = _first_present(rec, ["start", "pos", "start_pos"])
            if generic is not None:
                sg = _to_int(generic)
                if sg is not None:
                    s0 = (sg - 1) if node_start_is_1_based else sg
                else:
                    s0 = None
            else:
                s0 = None

        if chrom is not None and s0 is not None and length is not None:
            anchors[node_id] = (str(chrom), int(s0), strand if strand in ("+", "-") else "+", int(length))

    return anchors, nodes_seen, nodes_with_genomead_af


# ------------------------ VCF helpers ------------------------ #

def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")

def parse_vcf_header_and_count_records(vcf_path: str) -> Tuple[List[str], int]:
    header: List[str] = []
    n = 0
    with open_maybe_gzip(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            else:
                n += 1
    return header, n

def append_info(info: str, kv_pairs: List[Tuple[str, str]]) -> str:
    parts = [] if info.strip() == "." or info.strip() == "" else [info.strip()]
    for k, v in kv_pairs:
        if v is None or v == "":
            continue
        parts.append(f"{k}={v}")
    return ";".join(parts) if parts else "."


# ------------------------ Conversion core ------------------------ #

def convert_vcf(
    in_vcf: str,
    out_linear_vcf_gz: str,
    out_nodes_vcf_gz: str,
    anchors: Dict[int, Tuple[str, int, str, int]],
    nodes_seen_in_json: Set[int],
    ref_nodes_with_af: Set[int],
    alt_to_ref: Dict[int, int],
    tsv_ref_nodes: Set[int],
    tsv_alt_nodes: Set[int],
    offset_is_0_based: bool,
    sort_records: bool,
    compress_level: int,
    do_index: bool,
) -> None:
    header_lines, total = parse_vcf_header_and_count_records(in_vcf)

    kept_header = [
        h for h in header_lines
        if not h.startswith("##contig=")
           and not h.startswith("#CHROM")
           and not h.startswith("##fileformat")
    ]

    converted: List[Tuple[str, int, str, str, str, str, str]] = []
    unconverted: List[Tuple[Union[str,int], Union[str,int], str, str, str, str, str]] = []

    # Variant counters
    cnt_ref = cnt_refalt = cnt_altalt = cnt_nonnum = 0
    # DISTINCT node sets observed in VCF
    nodes_ref_seen: Set[int] = set()
    nodes_refalt_seen: Set[int] = set()
    nodes_altalt_seen: Set[int] = set()

    with open_maybe_gzip(in_vcf) as fin, \
         tqdm(total=total if total > 0 else None, desc="Convert", unit="rec",
              dynamic_ncols=True, leave=True) as bar:

        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                if total:
                    bar.update(1)
                continue

            chrom_node = parts[0]
            pos_offset = parts[1]
            vid  = parts[2] if len(parts) > 2 else "."
            ref  = parts[3] if len(parts) > 3 else "N"
            alt  = parts[4] if len(parts) > 4 else "N"
            qual = parts[5] if len(parts) > 5 else "."
            filt = parts[6] if len(parts) > 6 else "PASS"
            info = parts[7] if len(parts) > 7 else "."

            # Parse numeric node/offset
            try:
                node_id = int(chrom_node)
                offset  = int(pos_offset)
                is_numeric = True
            except Exception:
                is_numeric = False

            # Classify for counts + distinct node sets
            if not is_numeric:
                cnt_nonnum += 1
                info2 = append_info(info, [("NID", str(chrom_node)), ("NSO", str(pos_offset))])
                unconverted.append((chrom_node, pos_offset, vid, ref, alt, qual, filt, info2))
                if total:
                    bar.update(1)
                continue
            else:
                if node_id in nodes_seen_in_json:
                    if node_id in ref_nodes_with_af:
                        cnt_ref += 1
                        nodes_ref_seen.add(node_id)
                    else:
                        cnt_refalt += 1
                        nodes_refalt_seen.add(node_id)
                else:
                    cnt_altalt += 1
                    nodes_altalt_seen.add(node_id)

            base_node = node_id

            # Map ALT node to its REF node if present
            ref_node_candidate = alt_to_ref.get(base_node)
            anchor = None
            if ref_node_candidate is not None:
                anchor = anchors.get(ref_node_candidate)
            if anchor is None:
                anchor = anchors.get(base_node)
            if anchor is None and base_node in nodes_seen_in_json and ref_node_candidate is not None:
                anchor = anchors.get(ref_node_candidate)

            # Preserve original node/offset in INFO
            info2 = append_info(info, [("NID", str(node_id)), ("NSO", str(offset))])

            if anchor is None:
                # Keep in node coords; write POS as offset+1 (VCF is 1-based)
                pos_for_nodes_vcf = (offset if not offset_is_0_based else offset + 1)
                unconverted.append((str(node_id), int(pos_for_nodes_vcf), vid, ref, alt, qual, filt, info2))
                if total:
                    bar.update(1)
                continue

            chrom_lin, start0, strand, length = anchor
            off0 = offset if offset_is_0_based else (offset - 1)
            if strand == "+":
                pos_lin = start0 + off0 + 1
            else:
                pos_lin = start0 + (length - 1 - off0) + 1

            converted.append((chrom_lin, int(pos_lin), vid, ref, alt, qual, filt, info2))
            if total:
                bar.update(1)

    # Sorting (needed for tabix)
    need_sort = sort_records or do_index
    if need_sort:
        def chrom_key_linear(c: str):
            m = re.fullmatch(r"(?:chr)?(\d+)", c, flags=re.IGNORECASE)
            if m:
                return (0, int(m.group(1)))
            cl = c.lower()
            if cl in ("chrx", "x"): return (1, 23)
            if cl in ("chry", "y"): return (1, 24)
            if cl in ("chrm", "chrmt", "m", "mt"): return (1, 25)
            return (2, cl)
        converted.sort(key=lambda r: (chrom_key_linear(r[0]), int(r[1])))

        def chrom_key_nodes(c: Union[str,int]):
            try:
                return (0, int(c))
            except Exception:
                return (1, str(c))
        unconverted.sort(key=lambda r: (chrom_key_nodes(r[0]), int(r[1]) if str(r[1]).isdigit() else 0))

    # Collect contigs
    lin_contigs, seen_lin = [], set()
    for c, *_ in converted:
        if c not in seen_lin:
            seen_lin.add(c)
            lin_contigs.append(c)

    node_contigs, seen_node = [], set()
    for c, *_ in unconverted:
        sc = str(c)
        if sc not in seen_node:
            seen_node.add(sc)
            node_contigs.append(sc)

    def ensure_pysam_or_die():
        if pysam is None:
            print("ERROR: pysam is required for bgzip + tabix. Install: pip install pysam (or conda -c bioconda pysam)", file=sys.stderr)
            sys.exit(2)

    def write_plain_and_bgzip(final_gz_path: str,
                              records: List[Tuple],
                              contigs: List[str],
                              meta_desc: str):
        if not final_gz_path.endswith(".vcf.gz"):
            final_gz_path = final_gz_path + ".vcf.gz"
        os.makedirs(os.path.dirname(final_gz_path), exist_ok=True)
        tmp_vcf = final_gz_path[:-3]

        with open(tmp_vcf, "w", encoding="utf-8") as out:
            # MUST be first per spec
            out.write('##fileformat=VCFv4.2\n')

            # Then any kept header lines (skip duplicates)
            for h in kept_header:
                if h.startswith("#CHROM") or h.lower().startswith("##fileformat"):
                    continue
                out.write(h + "\n")

            # Then your custom META/INFO/contigs
            out.write(
                f'##META=<ID={meta_desc},Description="Converted using node map and optional ALT->REF TSV; NID/NSO store original node and offset.">\n')
            out.write('##INFO=<ID=NID,Number=1,Type=String,Description="Original node_id (CHROM in input)">\n')
            out.write('##INFO=<ID=NSO,Number=1,Type=String,Description="Original starting offset (POS in input)">\n')
            for c in contigs:
                out.write(f"##contig=<ID={c}>\n")

            # Finally the column header
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            with tqdm(total=len(records), desc=f"Write {meta_desc}", unit="rec",
                      dynamic_ncols=True, leave=True) as wbar:
                for chrom, pos, vid, ref, alt, qual, filt, info in records:
                    out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n")
                    wbar.update(1)

        ensure_pysam_or_die()
        if do_index:
            produced_gz = pysam.tabix_index(tmp_vcf, preset="vcf", force=True, keep_original=False)
            produced_tbi = produced_gz + ".tbi"
            if os.path.abspath(produced_gz) != os.path.abspath(final_gz_path):
                try:
                    if os.path.exists(final_gz_path):
                        os.remove(final_gz_path)
                except Exception:
                    pass
                os.replace(produced_gz, final_gz_path)
                if os.path.exists(produced_tbi):
                    os.replace(produced_tbi, final_gz_path + ".tbi")
            return final_gz_path, final_gz_path + ".tbi"
        else:
            if os.path.exists(final_gz_path):
                try:
                    os.remove(final_gz_path)
                except Exception:
                    pass
            try:
                pysam.tabix_compress(tmp_vcf, final_gz_path, force=True)
            except TypeError:
                pysam.tabix_compress(tmp_vcf, final_gz_path)
            try:
                os.remove(tmp_vcf)
            except Exception:
                pass
            return final_gz_path, None

    # Write outputs
    linear_gz, linear_tbi = write_plain_and_bgzip(out_linear_vcf_gz, converted, lin_contigs, meta_desc="CONVERTED")
    nodes_gz, nodes_tbi   = write_plain_and_bgzip(out_nodes_vcf_gz,  unconverted, node_contigs, meta_desc="UNCONVERTED")

    # PRINT distinct node counts (from VCF)
    print("=== Distinct node counts (observed in VCF) ===")
    print(f"reference nodes : {len(nodes_ref_seen)}")
    print(f"ref-alt nodes   : {len(nodes_refalt_seen)}")
    print(f"alt-alt nodes   : {len(nodes_altalt_seen)}")

    # TSV coverage of ref-alt nodes (observed in VCF) — check TSV **ALT** column
    refalt_missing_in_tsv = {n for n in nodes_refalt_seen if n not in tsv_alt_nodes}
    refalt_present_in_tsv = nodes_refalt_seen - refalt_missing_in_tsv
    print("\n=== TSV coverage of ref-alt nodes (VS TSV ALT column) ===")
    print(f"ref-alt nodes seen in VCF             : {len(nodes_refalt_seen)}")
    print(f"found in TSV ALT column               : {len(refalt_present_in_tsv)}")
    print(f"NOT found in TSV ALT column           : {len(refalt_missing_in_tsv)}")

    # (Optional) TSV presence of VCF nodes in general
    tsv_ref_seen = len({n for n in (nodes_ref_seen | nodes_refalt_seen | nodes_altalt_seen) if n in tsv_ref_nodes})
    tsv_alt_seen = len({n for n in (nodes_ref_seen | nodes_refalt_seen | nodes_altalt_seen) if n in tsv_alt_nodes})
    print("\n=== TSV node presence in VCF (sanity) ===")
    print(f"TSV REF nodes (unique)                : {len(tsv_ref_nodes)}")
    print(f"TSV REF nodes seen in VCF (any class) : {tsv_ref_seen}")
    print(f"TSV ALT nodes (unique)                : {len(tsv_alt_nodes)}")
    print(f"TSV ALT nodes seen in VCF (any class) : {tsv_alt_seen}")

    # PRINT variant counts
    total_classified = cnt_ref + cnt_refalt + cnt_altalt
    grand_total = total_classified + cnt_nonnum
    print("\n=== Variant class counts (by input node id) ===")
    print(f"reference : {cnt_ref}")
    print(f"ref-alt   : {cnt_refalt}")
    print(f"alt-alt   : {cnt_altalt}")
    print(f"----------------------------------")
    print(f"TOTAL (node CHROM only): {total_classified}")
    if cnt_nonnum:
        print(f"Non-numeric CHROM (ignored in classes): {cnt_nonnum}  (grand total lines = {grand_total})")

    print(f"\nDone.\n  Linear : {os.path.abspath(linear_gz)} (records={len(converted)})"
          f"{'' if not linear_tbi else f'  [index: {os.path.abspath(linear_tbi)}]'}\n"
          f"  Nodes  : {os.path.abspath(nodes_gz)} (records={len(unconverted)})"
          f"{'' if not nodes_tbi else f'  [index: {os.path.abspath(nodes_tbi)}]'}")


# ------------------------ CLI ------------------------ #

def main():
    ap = argparse.ArgumentParser(description="Convert node/offset VCF to linear-reference coordinates; output linear + node VCF.GZ (tabix-indexed) and print class & node counts.")
    ap.add_argument("--in_vcf",             required=True, help="Input VCF (.vcf or .vcf.gz) with CHROM=node_id, POS=offset.")
    ap.add_argument("--map_json",           required=True, help="Node map JSON (array) or JSONL.")
    ap.add_argument("--map_jsonl",          action="store_true", help="Interpret --map_json as JSONL.")
    ap.add_argument("--tsv",                default=None, help="Optional TSV mapping ALT_NODE(S) -> REF_NODE; coverage checks use TSV ALT column.")
    ap.add_argument("--out_linear_vcf",     required=True, help="Output path for linear VCF.GZ (will be bgzip compressed and indexed).")
    ap.add_argument("--out_nodes_vcf",      required=True, help="Output path for node-form VCF.GZ (will be bgzip compressed and indexed).")
    ap.add_argument("--offset_is_0_based",  type=lambda s: str(s).lower() in ("1","true","t","yes","y"), default=True,
                        help="Offsets in input VCF are 0-based (default: true). Set false if 1-based.")
    ap.add_argument("--node_start_is_1_based", type=lambda s: str(s).lower() in ("1","true","t","yes","y"), default=True,
                        help="Node starts in JSON are 1-based (default: true).")
    ap.add_argument("--sort",               action="store_true", help="Sort outputs by CHROM and POS (linear) / by node and offset (nodes).")
    ap.add_argument("--no-index",           action="store_true", help="Do not create Tabix indexes (.tbi).")
    ap.add_argument("--compress-level",     type=int, default=6, help="(Kept for compatibility; pysam.tabix_index does not accept level.)")
    args = ap.parse_args()

    def normalize_vcfgz(p: str) -> str:
        if p.endswith(".vcf.gz"):
            return p
        if p.endswith(".vcf"):
            return p + ".gz"
        return p + ".vcf.gz"

    out_linear_vcf_gz = normalize_vcfgz(args.out_linear_vcf)
    out_nodes_vcf_gz  = normalize_vcfgz(args.out_nodes_vcf)

    print("=== Convert VCF Configuration ===")
    print(f"in_vcf                : {os.path.abspath(args.in_vcf)}")
    print(f"map_json              : {os.path.abspath(args.map_json)} (jsonl={bool(args.map_jsonl)})")
    print(f"tsv (ALT->REF mapping): {os.path.abspath(args.tsv) if args.tsv else '(none)'}")
    print(f"out_linear_vcf.gz     : {os.path.abspath(out_linear_vcf_gz)}")
    print(f"out_nodes_vcf.gz      : {os.path.abspath(out_nodes_vcf_gz)}")
    print(f"offset_is_0_based     : {args.offset_is_0_based}")
    print(f"node_start_is_1_based : {args.node_start_is_1_based}")
    print(f"sort                  : {args.sort}")
    print(f"index                 : {not args.no_index}")
    print("=" * 70)

    alt_to_ref, tsv_ref_nodes, tsv_alt_nodes = load_alt_to_ref_tsv(args.tsv) if args.tsv else ({}, set(), set())
    if args.tsv:
        print(f"ALT->REF mappings loaded: {len(alt_to_ref)} "
              f"(unique TSV REF nodes={len(tsv_ref_nodes)}, unique TSV ALT nodes={len(tsv_alt_nodes)})")

    anchors, nodes_seen, ref_nodes_with_af = load_node_map(
        args.map_json,
        jsonl=args.map_jsonl,
        node_start_is_1_based=args.node_start_is_1_based
    )
    # JSON-level node summary
    n_ref_nodes_json  = len(ref_nodes_with_af)
    n_refalt_nodes_json = len(nodes_seen - ref_nodes_with_af)
    print(f"Node map: anchors={len(anchors)} (with chrom+start+strand+length), nodes_seen={len(nodes_seen)} "
          f"(reference nodes={n_ref_nodes_json}, ref-alt nodes={n_refalt_nodes_json})")

    convert_vcf(
        in_vcf=args.in_vcf,
        out_linear_vcf_gz=out_linear_vcf_gz,
        out_nodes_vcf_gz=out_nodes_vcf_gz,
        anchors=anchors,
        nodes_seen_in_json=nodes_seen,
        ref_nodes_with_af=ref_nodes_with_af,
        alt_to_ref=alt_to_ref,
        tsv_ref_nodes=tsv_ref_nodes,
        tsv_alt_nodes=tsv_alt_nodes,
        offset_is_0_based=args.offset_is_0_based,
        sort_records=args.sort,
        compress_level=args.compress_level,
        do_index=(not args.no_index),
    )


if __name__ == "__main__":
    main()
