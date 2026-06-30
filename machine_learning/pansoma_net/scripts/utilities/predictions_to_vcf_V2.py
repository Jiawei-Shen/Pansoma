#!/usr/bin/env python3
"""
Create a node/offset VCF from predictions produced on the *latest shard format*
(shard_XXXXX_data.npy + variant_summary.ndjson).

LATEST SHARD FORMAT (your current pipeline):
  - shard files:        shard_00000_data.npy, shard_00001_data.npy, ...
  - predictions JSON:   records contain:
        {'shard_path', 'index_in_shard', 'pred_class', 'pred_prob', 'probs', ...}
  - variant_summary:    NDJSON file (one line per tensor) containing at least:
        shard_index, index_within_shard,
        node_id, v_pos, v_ref, v_alt, v_type,
        alt_allele_frequency, coverage_at_locus, ref_allele_count_at_locus,
        alt_allele_count, other_allele_count_at_locus, mean_alt_allele_base_quality, ...

VCF semantics (same as before):
  - CHROM = node_id
  - POS   = v_pos (0-based in your meta; VCF is 1-based so we write POS=v_pos+1 by default)
  - REF/ALT = v_ref/v_alt from variant_summary
  - FILTER:
      * pred_class == "true" -> PASS (optionally gated by --min_true_prob)
      * otherwise, if --refcall_prob provided and prob >= it -> RefCall
      * else drop

INFO fields written (when available):
  PROB, AF, DP, AD, OTHER, BQ, TYPE, SHARD (shard_index), IDX (index_within_shard)

Output:
  - writes plain .vcf then bgzip+tabix index (requires pysam) unless --no-index.

Notes:
  - This script supports JSON array predictions (default) or JSONL via --jsonl.
  - This script does NOT require per-node variant_summary.json nor filename parsing.
"""

import argparse
import json
import os
import re
import sys
from typing import Dict, Iterator, List, Optional, Tuple

from tqdm import tqdm

try:
    import pysam
except Exception:
    pysam = None


# ----------------------------- helpers ------------------------------------- #

def ensure_pysam_or_die():
    if pysam is None:
        msg = (
            "ERROR: pysam is required to bgzip and index the VCF.\n"
            "Install with:  pip install pysam\n"
            "Or in conda:   conda install -c bioconda pysam"
        )
        print(msg, file=sys.stderr)
        sys.exit(2)


def escape_info(v) -> str:
    return (
        str(v)
        .replace(" ", "_")
        .replace(",", "%2C")
        .replace(";", "%3B")
        .replace("=", "%3D")
    )


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
    Extract probability for the positive/'true' class.
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


_SHARD_BASENAME_RE = re.compile(r"(?:^|/)(?:shard_)?(?P<idx>\d+)(?:_data)?\.npy$", re.IGNORECASE)

def shard_index_from_path(shard_path: str) -> Optional[int]:
    """
    Extract shard_index from paths like:
      /.../shard_00007_data.npy -> 7
      shard_12_data.npy -> 12
      shard_00007.npy -> 7
    """
    if not shard_path:
        return None
    bn = os.path.basename(str(shard_path))
    m = _SHARD_BASENAME_RE.search(bn)
    if not m:
        # fallback: try find first integer group in filename
        m2 = re.search(r"(\d+)", bn)
        if not m2:
            return None
        try:
            return int(m2.group(1))
        except Exception:
            return None
    try:
        return int(m.group("idx"))
    except Exception:
        return None


def iter_variant_summary_ndjson(path: str) -> Iterator[dict]:
    """
    Stream NDJSON lines (one JSON dict per line).
    """
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


def load_variant_summary_index(
    ndjson_path: str,
    shard_prefix: str = "shard_",
    shard_digits: int = 5,
) -> Dict[Tuple[int, int], dict]:
    """
    Build a mapping: (shard_index, index_within_shard) -> meta dict

    This is fast enough for ~0.5M records, and avoids any filename parsing.
    """
    idx: Dict[Tuple[int, int], dict] = {}
    missing = 0
    total = 0
    with tqdm(desc="Load variant_summary", unit="rec", dynamic_ncols=True, leave=True) as bar:
        for row in iter_variant_summary_ndjson(ndjson_path):
            total += 1
            bar.update(1)
            try:
                sidx = int(row.get("shard_index"))
                iws  = int(row.get("index_within_shard"))
            except Exception:
                missing += 1
                continue
            idx[(sidx, iws)] = row
    if missing:
        tqdm.write(f"WARNING: {missing} variant_summary rows missing shard_index/index_within_shard; skipped.")
    tqdm.write(f"Loaded variant_summary index: {len(idx)} / {total} rows")
    return idx


# ----------------------------- core ---------------------------------------- #

def scan_predictions_shard_mode(
    iterable: Iterator[dict],
    vs_index: Dict[Tuple[int, int], dict],
    min_true_prob: Optional[float],
    refcall_prob: Optional[float],
    pos_is_1based: bool,
) -> Tuple[List[Tuple[str, int, str, str, dict, str]], int, int, int, int]:
    """
    Returns:
      records: list of (node_id, pos, ref, alt, info_dict, filter_label)
      total_in, kept_out, missing_shard_fields, missing_meta, missing_prob
    """
    records: List[Tuple[str, int, str, str, dict, str]] = []
    total_in = 0
    kept_out = 0
    missing_shard_fields = 0
    missing_meta = 0
    missing_prob = 0

    bar = tqdm(total=None, desc="Scan preds", unit="rec", dynamic_ncols=True, leave=True)
    for obj in iterable:
        total_in += 1
        bar.update(1)

        shard_path = obj.get("shard_path")
        idx_in_shard = obj.get("index_in_shard")

        if shard_path is None or idx_in_shard is None:
            missing_shard_fields += 1
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        sidx = shard_index_from_path(str(shard_path))
        try:
            iws = int(idx_in_shard)
        except Exception:
            sidx = None

        if sidx is None:
            missing_shard_fields += 1
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        meta = vs_index.get((sidx, iws))
        if meta is None:
            missing_meta += 1
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        prob = true_prob(obj)
        if prob is None:
            missing_prob += 1
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

        # Extract CHROM/POS/REF/ALT from meta
        try:
            node_id = str(meta["node_id"])
        except Exception:
            missing_meta += 1
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        # Your meta uses v_pos (0-based). VCF is 1-based.
        # Default: write v_pos+1 unless user says it's already 1-based.
        try:
            vpos0 = int(meta.get("v_pos"))
        except Exception:
            missing_meta += 1
            bar.set_postfix_str(f"kept={kept_out}")
            continue

        pos = vpos0 if pos_is_1based else (vpos0 + 1)

        ref = str(meta.get("v_ref", "N"))
        alt = str(meta.get("v_alt", "N"))

        info: Dict[str, Optional[str]] = {
            "PROB": f"{prob:.6g}" if prob is not None else None,
            "TYPE": str(meta.get("v_type")) if meta.get("v_type") is not None else None,
            "SHARD": str(sidx),
            "IDX": str(iws),
        }

        af = meta.get("alt_allele_frequency")
        cov = meta.get("coverage_at_locus")
        rc  = meta.get("ref_allele_count_at_locus")
        ac  = meta.get("alt_allele_count")
        oc  = meta.get("other_allele_count_at_locus")
        bq  = meta.get("mean_alt_allele_base_quality")

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

        records.append((node_id, pos, ref, alt, info, filter_label))
        kept_out += 1
        bar.set_postfix_str(f"kept={kept_out}")

    bar.close()
    return records, total_in, kept_out, missing_shard_fields, missing_meta, missing_prob


# ----------------------------- main ---------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description="Generate VCF from shard-mode predictions + variant_summary.ndjson (latest format), then bgzip+tabix."
    )
    ap.add_argument("--predictions", required=True, help="predictions.json (array) or .jsonl")
    ap.add_argument("--jsonl", action="store_true", help="Interpret --predictions as JSONL")
    ap.add_argument("--variant_summary", required=True,
                    help="Path to variant_summary.ndjson produced by tensor generation (contains shard_index/index_within_shard).")
    ap.add_argument("--output_vcf", required=True,
                    help="Output path (.vcf or .vcf.gz). Final output will be .vcf.gz if indexing/compressing.")
    ap.add_argument("--sort", action="store_true", help="Sort by node_id (CHROM) then POS")
    ap.add_argument("--min_true_prob", type=float, default=None,
                    help="Keep only true records with prob>=this (if prob present)")
    ap.add_argument("--refcall_prob", type=float, default=None,
                    help="Also include non-true predictions with prob>=this, labeled FILTER=RefCall")
    ap.add_argument("--pos_is_1based", action="store_true",
                    help="If set, treat meta v_pos as already 1-based; otherwise write POS=v_pos+1 (default).")
    ap.add_argument("--no-index", action="store_true", help="Do not create a Tabix index (.tbi)")
    args = ap.parse_args()

    print("=== VCF Build Configuration ===")
    print(f"Predictions     : {os.path.abspath(args.predictions)}  (jsonl={bool(args.jsonl)})")
    print(f"Variant summary : {os.path.abspath(args.variant_summary)}")
    print(f"Output VCF      : {os.path.abspath(args.output_vcf)}")
    print(f"min_true_prob: {args.min_true_prob} | refcall_prob: {args.refcall_prob} | sort: {args.sort} | index: {not args.no_index}")
    print(f"POS convention  : {'1-based (no +1)' if args.pos_is_1based else 'VCF POS = v_pos + 1 (default)'}")
    print("=" * 60)

    if not os.path.isfile(args.variant_summary):
        raise SystemExit(f"variant_summary not found: {args.variant_summary}")

    # Load summary index
    vs_index = load_variant_summary_index(args.variant_summary)

    # Build iterator and scan predictions
    if args.jsonl:
        iterable = iter_jsonl(args.predictions)
    else:
        data = read_json_array(args.predictions)
        iterable = iter(data)

    records, total_in, kept_out, miss_fields, miss_meta, miss_prob = scan_predictions_shard_mode(
        iterable,
        vs_index=vs_index,
        min_true_prob=args.min_true_prob,
        refcall_prob=args.refcall_prob,
        pos_is_1based=args.pos_is_1based,
    )

    if miss_fields:
        tqdm.write(f"WARNING: {miss_fields} prediction records missing shard_path/index_in_shard; skipped.")
    if miss_meta:
        tqdm.write(f"WARNING: {miss_meta} prediction records could not be matched to variant_summary; skipped.")
    if miss_prob:
        tqdm.write(f"NOTE: {miss_prob} kept/checked records had no usable probability (PROB omitted).")

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
    out_dir = os.path.dirname(out_user) or "."
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
        out.write("##META=<ID=SCHEMA,Description=\"Custom: CHROM=node_id; POS from v_pos (+1 by default); source=variant_summary.ndjson\">\n")
        out.write("##INFO=<ID=PROB,Number=1,Type=Float,Description=\"Model probability for positive/'true' class\">\n")
        out.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Alt allele frequency from variant_summary.ndjson\">\n")
        out.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Coverage at locus from variant_summary.ndjson\">\n")
        out.write("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths: ref,alt\">\n")
        out.write("##INFO=<ID=OTHER,Number=1,Type=Integer,Description=\"Other allele count at locus\">\n")
        out.write("##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Mean alt allele base quality\">\n")
        out.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type from variant_summary (X/I/D)\">\n")
        out.write("##INFO=<ID=SHARD,Number=1,Type=Integer,Description=\"Shard index\">\n")
        out.write("##INFO=<ID=IDX,Number=1,Type=Integer,Description=\"Index within shard\">\n")
        out.write('##FILTER=<ID=RefCall,Description="Non-true prediction retained because prob >= --refcall_prob">\n')

        # Declare contigs (node IDs)
        seen = set()
        for chrom, *_ in records:
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
