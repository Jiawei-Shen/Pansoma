#!/usr/bin/env python3
"""
Create a node/offset VCF from predictions + per-node variant_summary.json,
then output BGZF-compressed VCF (.vcf.gz) and a Tabix index (.tbi).

Filename patterns (auto-detected):
  • 5-part (new):    <nodeID>_<offset>_<X>_<REF>_<ALT>.npy  -> CHROM from filename
  • 4-part (legacy): <offset>_<X>_<REF>_<ALT>.npy           -> CHROM from parent dir

Custom semantics:
  - CHROM = node_id
  - POS   = starting offset
  - REF/ALT parsed from filename
  - Classification:
      * pred_class == "true" -> keep as FILTER=PASS (optionally gated by --min_true_prob)
      * otherwise, if probability >= --refcall_prob -> keep as FILTER=RefCall
      * else drop

INFO fields written (when available from inputs):
  PROB   : probability for the 'true' class (or provided prob field)
  AF     : alt_allele_frequency (variant_summary.json)
  DP     : coverage_at_locus (variant_summary.json)
  AD     : "ref,alt" counts (ref_allele_count_at_locus, alt_allele_count)
  OTHER  : other_allele_count_at_locus
  BQ     : mean_alt_allele_base_quality
"""

import argparse
import json
import os
import re
import sys
from typing import Dict, Iterator, List, Optional, Tuple

from tqdm import tqdm

# pysam for bgzip + tabix
try:
    import pysam
except Exception:
    pysam = None

# Accept DNA letters (any case) OR the special '*' allele (spanning deletion)
_BASE_OR_STAR = r"(?:[ACGTNacgtn]+|\*)"

# 5-part: <nodeID>_<offset>_<X>_<REF>_<ALT>.npy
NAME_RE5 = re.compile(
    rf"""^(?P<node_id>-?\d+)
        _(?P<offset>-?\d+)
        _(?P<chrom_hint>[^_]+)
        _(?P<ref>{_BASE_OR_STAR})
        _(?P<alt>{_BASE_OR_STAR})
        \.[Nn][Pp][Yy]$
    """,
    re.VERBOSE,
)

# 4-part: <offset>_<X>_<REF>_<ALT>.npy
NAME_RE4 = re.compile(
    rf"""^(?P<offset>-?\d+)
        _(?P<chrom_hint>[^_]+)
        _(?P<ref>{_BASE_OR_STAR})
        _(?P<alt>{_BASE_OR_STAR})
        \.[Nn][Pp][Yy]$
    """,
    re.VERBOSE,
)

def _strip_node_prefix(bn: str) -> str:
    """
    If basename looks like "<nodeID>_<rest>", drop the leading "<nodeID>_".
    """
    i = bn.find("_")
    if i > 0:
        maybe = bn[:i]
        if maybe.lstrip("-").isdigit():
            return bn[i+1:]
    return bn

def read_json_array(path: str) -> List[dict]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list):
        raise ValueError("Predictions JSON must be a list when not using --jsonl.")
    return [x for x in data if isinstance(x, dict)]

def iter_jsonl(path: str) -> Iterator[dict]:
    with open(path, "r", encoding="utf-8") as f:
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

def true_prob(rec: dict) -> Optional[float]:
    """
    Extract a probability for the positive/'true' class.
    Priority:
      1) rec['probs']['true'] (case-insensitive key match)
      2) rec['pred_prob']
    """
    p = rec.get("probs")
    if isinstance(p, dict):
        for k, v in p.items():
            if str(k).lower() == "true":
                try:
                    return float(v)
                except Exception:
                    break
    try:
        return float(rec.get("pred_prob"))
    except Exception:
        return None

def load_variant_summary_for_node(node_dir: str) -> Dict[str, dict]:
    """
    Build an index: basename(.npy) -> row dict from node_dir/variant_summary.json
    Stores both 4-part and 5-part keys to match predictions reliably.
    """
    idx: Dict[str, dict] = {}
    path = os.path.join(node_dir, "variant_summary.json")
    if not os.path.exists(path):
        return idx
    try:
        with open(path, "r", encoding="utf-8") as f:
            summary = json.load(f)
    except Exception:
        return idx

    rows: List[dict] = []
    if isinstance(summary, dict) and isinstance(summary.get("variants_passing_af_filter"), list):
        rows = [r for r in summary["variants_passing_af_filter"] if isinstance(r, dict)]
    elif isinstance(summary, list):
        rows = [r for r in summary if isinstance(r, dict)]

    node_id_guess = os.path.basename(node_dir.rstrip("/"))

    for r in rows:
        tf = r.get("tensor_file") or r.get("variant_key")
        if isinstance(tf, str):
            bn = os.path.basename(tf)
            if not bn.lower().endswith(".npy"):
                bn = bn + ".npy"
            # 4-part key
            idx[bn] = r
            # 5-part key "<nodeID>_<4part>"
            if node_id_guess and node_id_guess.lstrip("-").isdigit():
                idx[f"{node_id_guess}_{bn}"] = r
    return idx

def escape_info(v) -> str:
    return (
        str(v)
        .replace(" ", "_")
        .replace(",", "%2C")
        .replace(";", "%3B")
        .replace("=", "%3D")
    )

def ensure_pysam_or_die():
    if pysam is None:
        msg = (
            "ERROR: pysam is required to bgzip and index the VCF.\n"
            "Install with:  pip install pysam\n"
            "Or in conda:   conda install -c bioconda pysam"
        )
        print(msg, file=sys.stderr)
        sys.exit(2)

def parse_filename(path: str) -> Optional[Tuple[str, int, str, str]]:
    """
    Parse tensor filename to (node_id, offset, ref, alt).
    Supports:
      - <nodeID>_<offset>_<X>_<REF>_<ALT>.npy  (node_id from filename)
      - <offset>_<X>_<REF>_<ALT>.npy           (node_id from parent directory)
    """
    base = os.path.basename(path)

    m5 = NAME_RE5.match(base)
    if m5:
        node_id = m5.group("node_id")
        try:
            offset = int(m5.group("offset"))
        except Exception:
            return None
        ref = m5.group("ref")
        alt = m5.group("alt")
        return node_id, offset, ref, alt

    m4 = NAME_RE4.match(base)
    if m4:
        node_id = os.path.basename(os.path.dirname(path.rstrip("/")))
        try:
            offset = int(m4.group("offset"))
        except Exception:
            return None
        ref = m4.group("ref")
        alt = m4.group("alt")
        return node_id, offset, ref, alt

    return None

def scan_predictions(
    iterable: Iterator[dict],
    min_true_prob: Optional[float],
    refcall_prob: Optional[float],
) -> Tuple[List[Tuple[str, int, str, str, dict, str]], int, int, int]:
    """
    Returns:
      records: list of (node_id, offset, ref, alt, info_dict, filter_label)
      total_in, kept_out, missing_pattern
    """
    records: List[Tuple[str, int, str, str, dict, str]] = []
    total_in = 0
    kept_out = 0
    missing_pattern = 0
    by_node_cache: Dict[str, Dict[str, dict]] = {}

    bar = tqdm(total=None, desc="Scan preds", unit="rec", dynamic_ncols=True, leave=True)
    for obj in iterable:
        total_in += 1
        bar.update(1)

        # Extract essential fields
        path = obj.get("path") or obj.get("file") or obj.get("src") or obj.get("tensor") or obj.get("npy")
        if not path:
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        parsed = parse_filename(path)
        if not parsed:
            missing_pattern += 1
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        node_id, offset, ref, alt = parsed
        prob = true_prob(obj)
        pred_class = str(obj.get("pred_class", "")).strip().lower()

        # Decide keep and FILTER label
        keep = False
        filter_label = "PASS"
        if pred_class == "true":
            if (min_true_prob is None) or (prob is None) or (prob >= min_true_prob):
                keep = True
                filter_label = "PASS"
        else:
            if (refcall_prob is not None) and (prob is not None) and (prob >= refcall_prob):
                keep = True
                filter_label = "RefCall"

        if not keep:
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        # Load node's summary (cached)
        node_dir = os.path.dirname(path)
        if node_id not in by_node_cache:
            by_node_cache[node_id] = load_variant_summary_for_node(node_dir)

        base = os.path.basename(path)
        info = {
            "PROB": f"{prob:.6g}" if prob is not None else None,
        }

        # Try exact base (5-part), then the stripped 4-part key.
        row = by_node_cache[node_id].get(base)
        if row is None:
            base4 = _strip_node_prefix(base)
            if base4 != base:
                row = by_node_cache[node_id].get(base4)

        if row:
            af = row.get("alt_allele_frequency")
            cov = row.get("coverage_at_locus")
            rc  = row.get("ref_allele_count_at_locus")
            ac  = row.get("alt_allele_count")
            oc  = row.get("other_allele_count_at_locus")
            bq  = row.get("mean_alt_allele_base_quality")

            if af is not None:
                try:
                    info["AF"] = f"{float(af):.6g}"
                except Exception:
                    info["AF"] = escape_info(af)
            if cov is not None:
                if isinstance(cov, (int, float)):
                    try:
                        info["DP"] = str(int(cov))
                    except Exception:
                        info["DP"] = escape_info(cov)
                else:
                    info["DP"] = escape_info(cov)
            if rc is not None and ac is not None:
                try:
                    info["AD"] = f"{int(rc)},{int(ac)}"
                except Exception:
                    info["AD"] = f"{rc},{ac}"
            if oc is not None:
                try:
                    info["OTHER"] = str(int(oc))
                except Exception:
                    info["OTHER"] = escape_info(oc)
            if bq is not None:
                try:
                    info["BQ"] = f"{float(bq):.3f}"
                except Exception:
                    info["BQ"] = escape_info(bq)

        records.append((node_id, offset, ref, alt, info, filter_label))
        kept_out += 1
        bar.set_postfix_str(f"kept={kept_out}")
    bar.close()

    return records, total_in, kept_out, missing_pattern

def main():
    ap = argparse.ArgumentParser(description="Generate node/offset VCF from predictions + variant_summary.json, then bgzip+index. Supports 4- and 5-part filenames.")
    ap.add_argument("--predictions", required=True, help="predictions.json (array) or .jsonl")
    ap.add_argument("--jsonl", action="store_true", help="Interpret --predictions as JSONL")
    ap.add_argument("--output_vcf", required=True, help="Output path (.vcf or .vcf.gz). Final output will be .vcf.gz")
    ap.add_argument("--sort", action="store_true", help="Sort by node_id (CHROM) then POS")
    ap.add_argument("--min_true_prob", type=float, default=None, help="Keep only true records with prob>=this (if prob present)")
    ap.add_argument("--refcall_prob", type=float, default=None,
                    help="Also include non-true predictions with prob>=this, labeled FILTER=RefCall")
    ap.add_argument("--no-index", action="store_true", help="Do not create a Tabix index (.tbi)")
    args = ap.parse_args()

    print("=== VCF Build Configuration ===")
    print(f"Predictions: {os.path.abspath(args.predictions)}  (jsonl={bool(args.jsonl)})")
    print(f"Output VCF : {os.path.abspath(args.output_vcf)}")
    print(f"min_true_prob: {args.min_true_prob} | refcall_prob: {args.refcall_prob} | sort: {args.sort} | index: {not args.no_index}")
    print("=" * 40)

    # Build iterator and scan
    if args.jsonl:
        iterable = iter_jsonl(args.predictions)
    else:
        data = read_json_array(args.predictions)
        iterable = iter(data)

    records, total_in, kept_out, missing_pattern = scan_predictions(
        iterable, args.min_true_prob, args.refcall_prob
    )

    if missing_pattern:
        tqdm.write(
            "WARNING: {n} filenames did not match either "
            "'<nodeID>_<offset>_<X>_<REF>_<ALT>.npy' or '<offset>_<X>_<REF>_<ALT>.npy'; skipped."
            .format(n=missing_pattern)
        )

    if kept_out == 0:
        print("No qualifying predictions to write. Nothing to do.", file=sys.stderr)
        sys.exit(0)

    # Sorting: required for tabix indexing to be valid
    want_index = (not args.no_index)
    if args.sort or want_index:
        def chrom_key(c: str):
            try:
                return (0, int(c))
            except Exception:
                return (1, c)
        records.sort(key=lambda r: (chrom_key(r[0]), r[1]))
        if not args.sort:
            print("Note: Sorting enabled automatically to allow Tabix indexing.", file=sys.stderr)

    # Resolve output paths
    out_user = os.path.abspath(args.output_vcf)
    out_dir = os.path.dirname(out_user)
    os.makedirs(out_dir, exist_ok=True)

    # We'll write a temporary plain VCF, then use pysam.tabix_index to auto-compress+index.
    base_no_ext = out_user
    if out_user.endswith(".vcf.gz"):
        base_no_ext = out_user[:-3]   # strip ".gz" -> ".vcf"
    elif out_user.endswith(".vcf"):
        base_no_ext = out_user
    else:
        base_no_ext = out_user + ".vcf"

    tmp_vcf_path = base_no_ext

    # Write plain VCF
    with open(tmp_vcf_path, "w", encoding="utf-8") as out:
        # Header
        out.write("##fileformat=VCFv4.2\n")
        out.write("##META=<ID=SCHEMA,Description=\"Custom: CHROM=node_id, POS=starting_offset; supports 4- and 5-part filenames\">\n")
        out.write("##INFO=<ID=PROB,Number=1,Type=Float,Description=\"Model probability for positive/'true' class\">\n")
        out.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Alt allele frequency from variant_summary.json\">\n")
        out.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Coverage at locus from variant_summary.json\">\n")
        out.write("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths: ref,alt (from variant_summary.json)\">\n")
        out.write("##INFO=<ID=OTHER,Number=1,Type=Integer,Description=\"Other allele count at locus\">\n")
        out.write("##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Mean alt allele base quality\">\n")
        out.write('##FILTER=<ID=RefCall,Description="Non-true prediction retained because prob >= --refcall_prob">\n')

        # Declare contigs (node IDs)
        seen = set()
        for chrom, _, _, _, _, _ in records:
            if chrom not in seen:
                out.write(f"##contig=<ID={chrom}>\n")
                seen.add(chrom)

        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        with tqdm(total=len(records), desc="Write VCF", unit="rec", dynamic_ncols=True, leave=True) as wbar:
            for chrom, pos, ref, alt, info, filt in records:
                info_str = ";".join(f"{k}={escape_info(v)}" for k, v in info.items() if v is not None) or "."
                out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info_str}\n")
                wbar.update(1)

    print(f"Plain VCF written: {tmp_vcf_path} (records={len(records)}, scanned={total_in}, kept_out={kept_out})")

    # Compress + (optionally) index using pysam.tabix_index (auto-compress)
    ensure_pysam_or_die()
    if want_index:
        print(f"Compressing + indexing with pysam.tabix_index: {tmp_vcf_path} ...")
        final_vcf_gz = pysam.tabix_index(
            tmp_vcf_path,
            preset="vcf",
            force=True,
            keep_original=False
        )
        final_tbi = final_vcf_gz + ".tbi"
    else:
        # Only compress, no index
        final_vcf_gz = base_no_ext if base_no_ext.endswith(".vcf.gz") else base_no_ext + ".gz"
        print(f"Compressing with pysam.tabix_compress (no index): {final_vcf_gz} ...")
        try:
            pysam.tabix_compress(tmp_vcf_path, final_vcf_gz, force=True)
        except TypeError:
            pysam.tabix_compress(tmp_vcf_path, final_vcf_gz)
        final_tbi = None
        try:
            os.remove(tmp_vcf_path)
        except Exception:
            pass

    # If caller explicitly asked for .vcf.gz, ensure final path matches exactly
    out_user_is_gz = out_user.endswith(".vcf.gz")
    if out_user_is_gz and os.path.abspath(final_vcf_gz) != os.path.abspath(out_user):
        try:
            os.replace(final_vcf_gz, out_user)
            if final_tbi and os.path.exists(final_tbi):
                os.replace(final_tbi, out_user + ".tbi")
            final_vcf_gz = out_user
            final_tbi = out_user + ".tbi" if want_index else None
        except Exception:
            pass

    print(f"VCF (BGZF) written: {final_vcf_gz}")
    if want_index:
        print(f"VCF index (.tbi) written: {final_tbi}")

if __name__ == "__main__":
    main()
